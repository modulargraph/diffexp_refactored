#include "diffexp/level_preparation.hpp"
#include "diffexp/families.hpp"
#include <iostream>
using namespace diffexp;
void check(bool c,const char* why){if(!c)throw std::runtime_error(why);}
void hermetic() {
  ExactField field({"x","d"});Exact x(field,"x"),d(field,"d");
  ibp::ScalarProducts space(1,{},x);auto denominator=space.zero();denominator.constant=x;denominator.linear[0]=-x.constant(1);
  ibp::PropagatorBasis basis({space,{denominator},{}});level::Options options;std::size_t calls=0,last=0;
  auto provider=[&](const auto& batch,const auto&) {
    check(batch.size()>last,"closure demand must grow between unresolved passes");last=batch.size();++calls;
    fire::Result out;out.success=true;
    for(const auto& a:batch){if(a==ibp::Integral{1})out.reductions.emplace(a,ibp::Relation{{{1},x.constant(1)}});
      else if(a==ibp::Integral{2})out.reductions.emplace(a,ibp::Relation{{{1},Exact(field,"(1-d/2)/x")}});
      else throw std::runtime_error("unexpected hermetic demand");}
    return out;
  };
  auto closed=level::prepare(basis,d,field,0,{{1},{1}},options,provider);
  check(closed.success&&calls==2&&closed.ordered_basis==std::vector<ibp::Integral>{{1}},"batched tadpole closure");
  check(closed.matrix[0][0]==Exact(field,"(d/2-1)/x")&&closed.target_rows.size()==2&&closed.target_rows[0][0]==x.constant(1),"exact differential and duplicate observable coordinates");
  auto identity=[&](const auto& batch,const auto&){fire::Result out;out.success=true;for(const auto& a:batch)out.reductions.emplace(a,ibp::Relation{{a,x.constant(1)}});return out;};
  options.max_passes=2;auto limited=level::prepare(basis,d,field,0,{{1}},options,identity);
  check(!limited.success&&limited.passes==2&&limited.reason.find("pass budget")!=std::string::npos,"unresolved master expansion is finitely bounded");
  options.max_passes=8;options.max_demands=1;limited=level::prepare(basis,d,field,0,{{1}},options,identity);
  check(!limited.success&&limited.passes==1&&limited.reason.find("demand budget")!=std::string::npos,"demand cap checked before provider");
  options.max_demands=2000;
  auto zero=[&](const auto& batch,const auto&){fire::Result out;out.success=true;for(const auto& a:batch)out.reductions.emplace(a,ibp::Relation{});return out;};
  auto vanishing=level::prepare(basis,d,field,0,{{1}},options,zero);
  check(vanishing.success&&vanishing.ordered_basis.empty()&&vanishing.target_rows.size()==1,"zero observable span");
}
int main(int argc,char** argv) {try {
  hermetic();
  ExactField field({"x","d"});Exact x(field,"x"),d(field,"d");
  auto family=ibp::merge(ibp::quadratic_family(feynman::banana(2,{Rational(1),Rational(1),Rational(1)}),x),0,1,x);
  ibp::PropagatorBasis basis(family);level::Options options;
  auto invalid=level::prepare(basis,d,field,0,{},options);check(!invalid.success&&!invalid.reason.empty(),"empty closure demand rejected");
  if(argc<2){std::cout<<"level preparation validation passed; supply FIRE executable for integration\n";return 0;}
  options.provider.executable=argv[1];
  auto result=level::prepare(basis,d,field,0,{{3,1,0,0,0}},options);
  std::cout<<"sunrise passes="<<result.passes<<" masters="<<result.ordered_basis.size()<<" demands="<<result.demands<<" seconds="<<result.elapsed_seconds<<" reason="<<result.reason<<"\n";
  check(result.success,"merged sunrise closure failed");
  ibp::Generator generator(basis,d);ibp::ExactReducer reducer(x,300);
  ibp::for_each_seed(2,5,{1,1,50},[&](const ibp::Integral& seed){for(auto& row:generator.relations(seed))reducer.insert(std::move(row));});
  for(std::size_t i=0;i<result.ordered_basis.size();++i) {
    auto residual=generator.derivative(result.ordered_basis[i],0);
    for(std::size_t j=0;j<result.ordered_basis.size();++j)ibp::add(residual,result.ordered_basis[j],-result.matrix[i][j]);
    check(reducer.reduce(residual).remainder.empty(),"closed differential row disagrees with native IBPs");
  }
  ibp::Relation residual{{{3,1,0,0,0},x.constant(1)}};
  for(std::size_t j=0;j<result.ordered_basis.size();++j)ibp::add(residual,result.ordered_basis[j],-result.target_rows[0][j]);
  check(reducer.reduce(residual).remainder.empty(),"closed observable row disagrees with native IBPs");
  check(result.differential_witnesses.size()==result.matrix.size()&&result.target_witnesses.size()==result.target_rows.size(),"durable closure certificates");
  if(argc>2&&std::string(argv[2])=="--banana4") {
    auto banana=ibp::quadratic_family(feynman::banana(4,{Rational(1),Rational(1),Rational(1),Rational(1),Rational(1)}),x);
    banana=ibp::merge(banana,0,1,x.constant(Rational("1/3")));
    banana=ibp::merge(banana,0,1,x.constant(Rational("2/5")));
    banana=ibp::merge(banana,0,1,x);
    ibp::PropagatorBasis bb(banana);ibp::Integral target(bb.denominators.size(),0);target[0]=4;target[1]=1;
    options.provider.timeout_seconds=45;options.total_timeout_seconds=90;
    auto br=level::prepare(bb,d,field,0,{target},options);
    std::cout<<"banana4 success="<<br.success<<" slots="<<bb.denominators.size()<<" passes="<<br.passes<<" masters="<<br.ordered_basis.size()<<" demands="<<br.demands<<" seconds="<<br.elapsed_seconds<<" reason="<<br.reason<<"\n";
    check(br.success,"two-propagator Banana4 closure failed");
    for(const auto& dir:br.fire_directories)std::cout<<dir<<"\n";
  }
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
