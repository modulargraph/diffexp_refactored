// Opt-in experiment only. No production recurrence or conditioning changes.
#include "diffexp/adjoint_transport.hpp"
#include <chrono>
#include <ctime>
#include <iostream>
#include <memory>
using namespace diffexp;
using B=Jet::Ball;
struct AcbArray {
  slong size;acb_ptr values;
  explicit AcbArray(slong n):size(n),values(_acb_vec_init(n)){}
  ~AcbArray(){_acb_vec_clear(values,size);}
  AcbArray(const AcbArray&)=delete;
};
Boundary experimental_chart(const std::vector<adjoint_detail::Entry>& entries,const Boundary& initial,
    const B& center,const B& step,unsigned order,unsigned grouped) {
  const unsigned d=initial.size(),width=initial[0].size();const auto bits=B::precision();
  Jet x(0,order+1,bits);x.set(0,center);x.set(1,B(1));
  auto imaginary=x.constant(0);imaginary.set(0,B::from_strings("0","1"));
  AcbArray coefficients(static_cast<slong>(entries.size())*order);
  for(unsigned e=0;e<entries.size();++e){auto jet=evaluate(entries[e].coefficient,x,{{"x",x},{"I",imaginary}});
    for(unsigned n=0;n<order;++n){auto value=jet.at(n);acb_swap(coefficients.values+e*order+n,value.raw());}}
  const slong stride=static_cast<slong>(d)*width;
  AcbArray values(static_cast<slong>(order+1)*stride);
  const auto at=[&](unsigned n,unsigned i,unsigned k){return values.values+static_cast<slong>(n)*stride+i*width+k;};
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<width;++k)acb_set(at(0,i,k),initial[i][k].raw());
  B sum;
  for(unsigned n=0;n<order;++n){
    for(unsigned e=0;e<entries.size();++e)for(const auto& [row,column]:entries[e].positions)
      for(unsigned k=entries[e].epsilon;k<width;++k){
        auto target=at(n+1,row,k);auto coefficient=coefficients.values+e*order;
        if(grouped){
          // Real Acb arrays permit the documented signed stride. This does
          // not reinterpret C++ wrapper objects as C arrays or copy the lag
          // vectors into temporary gather buffers.
          if(grouped==2)acb_dot_precise(sum.raw(),target,0,coefficient,1,at(n,column,k-entries[e].epsilon),-stride,n+1,bits);
          else acb_dot(sum.raw(),target,0,coefficient,1,at(n,column,k-entries[e].epsilon),-stride,n+1,bits);
          acb_swap(target,sum.raw());
        }else for(unsigned m=0;m<=n;++m)
          if(!acb_is_zero(coefficient+m) && !acb_is_zero(at(n-m,column,k-entries[e].epsilon)))
            acb_addmul(target,coefficient+m,at(n-m,column,k-entries[e].epsilon),bits);
      }
    for(unsigned i=0;i<d;++i)for(unsigned k=0;k<width;++k)acb_div_ui(at(n+1,i,k),at(n+1,i,k),n+1,bits);
  }
  Boundary out(d,std::vector<B>(width,B(0)));
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<width;++k){
    acb_zero(out[i][k].raw());
    for(unsigned n=order+1;n-->0;){
      acb_mul(out[i][k].raw(),out[i][k].raw(),step.raw(),bits);
      acb_add(out[i][k].raw(),out[i][k].raw(),at(n,i,k),bits);
    }
  }
  return out;
}
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
struct Timing {double wall,cpu;};
template<class F>Timing time_charts(unsigned repeats,F run){
  auto start=std::chrono::steady_clock::now();auto cpu=std::clock();
  for(unsigned i=0;i<repeats;++i){auto out=run();require(out.back().back().is_finite(),"nonfinite timed result");}
  return {std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/repeats,double(std::clock()-cpu)/CLOCKS_PER_SEC/repeats};
}
int main(int argc,char** argv){try{
  if(argc<2 || std::string(argv[1])!="--run"){std::cout<<"Opt-in: benchmark_adjoint_dot --run [repetitions,1..20] [precise]\n";return 0;}
  unsigned repeats=argc>2?std::stoul(argv[2]):3;require(repeats>0 && repeats<=20,"repetition budget");
  const unsigned grouped_mode=argc>3 && std::string(argv[3])=="precise"?2:1;
  constexpr unsigned order=80,high=5,masters=10,observables=4,d=masters*observables+1;
  B::set_precision(384);ExactField field({"x","eps","I"});Exact x(field,"x"),one(field,"1");
  std::vector<RationalLineEntry> sparse;
  for(unsigned block=0;block<observables;++block)for(unsigned row=0;row<masters;++row)
    for(unsigned edge=0;edge<3;++edge){
      auto coefficient=(edge%2?-one:one)/(one.constant(2+edge)+x);
      sparse.push_back({block*masters+row,block*masters+(row+edge)%masters,edge,coefficient});
    }
  for(unsigned block=0;block<observables;++block)sparse.push_back({block*masters+3,d-1,0,one/(one+x*x)});
  auto compiled=adjoint_detail::compile(sparse);
  const auto center=B::from_strings("1/8","1/16"),step=B::from_strings("1/32","1/64");
  for(bool dense:{false,true}){
    Boundary initial(d,std::vector<B>(high+1,B(0)));initial.back()[0]=B(1);
    for(unsigned i=0;i<d-1;++i)for(unsigned k=0;k<=high;++k)if(dense || (k==0 && i%masters==0))
      initial[i][k]=B::from_strings(std::to_string(1+(i+3*k)%7)+"/16",std::to_string(static_cast<int>(i%3)-1)+"/32");
    auto reference=adjoint_detail::chart(compiled,initial,center,step,order);
    auto scalar=experimental_chart(compiled,initial,center,step,order,false);
    auto grouped=experimental_chart(compiled,initial,center,step,order,grouped_mode);
    std::size_t contained=0,nonzero=0;double worst_radius_ratio=0;
    slong minimum_reference_accuracy=1000000,minimum_grouped_accuracy=1000000;
    for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=high;++k){
      require(acb_equal(reference[i][k].raw(),scalar[i][k].raw()),"raw scalar experiment changed the reference recurrence");
      require(grouped[i][k].is_finite() && acb_overlaps(reference[i][k].raw(),grouped[i][k].raw()),"grouped retained enclosure does not overlap reference");
      contained+=acb_contains(reference[i][k].raw(),grouped[i][k].raw());
      const auto oldrad=std::max(mag_get_d(arb_radref(acb_realref(reference[i][k].raw()))),mag_get_d(arb_radref(acb_imagref(reference[i][k].raw()))));
      const auto newrad=std::max(mag_get_d(arb_radref(acb_realref(grouped[i][k].raw()))),mag_get_d(arb_radref(acb_imagref(grouped[i][k].raw()))));
      if(oldrad)worst_radius_ratio=std::max(worst_radius_ratio,newrad/oldrad);
      if(!reference[i][k].is_zero()){
        ++nonzero;minimum_reference_accuracy=std::min(minimum_reference_accuracy,acb_rel_accuracy_bits(reference[i][k].raw()));
        minimum_grouped_accuracy=std::min(minimum_grouped_accuracy,acb_rel_accuracy_bits(grouped[i][k].raw()));
      }
    }
    auto current=time_charts(repeats,[&]{return adjoint_detail::chart(compiled,initial,center,step,order);});
    auto raw=time_charts(repeats,[&]{return experimental_chart(compiled,initial,center,step,order,false);});
    auto dot=time_charts(repeats,[&]{return experimental_chart(compiled,initial,center,step,order,grouped_mode);});
    std::cout<<"kernel="<<(grouped_mode==2?"precise-dot":"dot")<<" case="<<(dense?"dense":"sparse")<<" order="<<order<<" eps_high="<<high<<" bits=384 dimension="<<d
      <<" current_cpu="<<current.cpu<<" raw_scalar_cpu="<<raw.cpu<<" grouped_cpu="<<dot.cpu
      <<" speedup="<<current.cpu/dot.cpu<<" current_wall="<<current.wall<<" grouped_wall="<<dot.wall
      <<" contained="<<contained<<"/"<<d*(high+1)<<" nonzero="<<nonzero<<" worst_radius_ratio="<<worst_radius_ratio
      <<" min_accuracy="<<minimum_reference_accuracy<<"/"<<minimum_grouped_accuracy<<std::endl;
  }
  std::cout<<"All retained enclosures overlap; raw scalar outputs match bitwise. No omitted-tail certification or full-example speed claim.\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
