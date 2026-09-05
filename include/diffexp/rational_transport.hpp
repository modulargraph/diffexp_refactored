#pragma once
#include "diffexp/canonical.hpp"
#include "diffexp/geometry.hpp"

namespace diffexp {
struct RationalLineEntry {unsigned row,column,epsilon;Exact coefficient;};
struct PreparedRationalEntry {unsigned row,column,epsilon;data::Expr coefficient;};
inline std::vector<PreparedRationalEntry> compile_rational_entries(const std::vector<RationalLineEntry>& entries) {
  std::vector<PreparedRationalEntry> result;
  for(const auto& entry:entries) {
    const auto& names=entry.coefficient.variables();
    for(const auto& terms:{entry.coefficient.numerator_terms(),entry.coefficient.denominator_terms()})
      for(const auto& term:terms) for(std::size_t i=0;i<term.powers.size();++i)
        if(term.powers[i] && names[i]!="x" && names[i]!="I")
          throw std::invalid_argument("rational line has an unsubstituted parameter: "+names[i]);
    result.push_back({entry.row,entry.column,entry.epsilon,data::Reader(entry.coefficient.str()).read()});
  }
  return result;
}

// Ordinary physical-basis evolution for an epsilon-polynomial connection. No
// homogeneous basis or spectral transform is constructed for a regular step.
inline Boundary rational_chart(const std::vector<PreparedRationalEntry>& entries,
    const Boundary& boundary,const Jet::Ball& center,const Jet::Ball& step,unsigned order) {
  using B=Jet::Ball;
  if(boundary.empty() || boundary[0].empty() || !order)throw std::invalid_argument("empty rational-chart boundary or Taylor order");
  const unsigned d=boundary.size(),kmax=boundary.at(0).size()-1;
  for(const auto& row:boundary)if(row.size()!=kmax+1)throw std::invalid_argument("rational-chart boundary window mismatch");
  const auto bits=B::precision();Jet x(0,order+1,bits);x.set(0,center);x.set(1,B(1));
  struct Entry {unsigned row,col,epsilon;std::vector<B> coefficients;};std::vector<Entry> matrix;
  auto imaginary=x.constant(0);imaginary.set(0,B::from_strings("0","1"));
  for(const auto& entry:entries) {
    if(entry.row>=d || entry.column>=d)throw std::invalid_argument("rational matrix index exceeds boundary size");
    auto jet=evaluate(entry.coefficient,x,{{"x",x},{"I",imaginary}});
    Entry e{entry.row,entry.column,entry.epsilon,{}};
    for(unsigned n=0;n<order;++n)e.coefficients.push_back(jet.at(n));matrix.push_back(std::move(e));
  }
  std::vector<B> values(static_cast<std::size_t>(order+1)*d*(kmax+1),B(0));
  const auto at=[&](unsigned n,unsigned i,unsigned k)->B&{return values[(static_cast<std::size_t>(n)*d+i)*(kmax+1)+k];};
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=kmax;++k)at(0,i,k)=boundary[i][k];
  for(unsigned n=0;n<order;++n) {
    for(const auto& e:matrix)for(unsigned k=e.epsilon;k<=kmax;++k)
      for(unsigned m=0;m<=n;++m)
        if(!e.coefficients[m].is_zero() && !at(n-m,e.col,k-e.epsilon).is_zero())
          acb_addmul(at(n+1,e.row,k).raw(),e.coefficients[m].raw(),at(n-m,e.col,k-e.epsilon).raw(),bits);
    for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=kmax;++k)
      acb_div_ui(at(n+1,i,k).raw(),at(n+1,i,k).raw(),n+1,bits);
  }
  Boundary out(d,std::vector<B>(kmax+1,B(0)));
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=kmax;++k)
    for(unsigned n=order+1;n-->0;)out[i][k]=out[i][k]*step+at(n,i,k);
  return out;
}
inline Boundary rational_chart(const std::vector<RationalLineEntry>& entries,
    const Boundary& boundary,const Jet::Ball& center,const Jet::Ball& step,unsigned order) {
  return rational_chart(compile_rational_entries(entries),boundary,center,step,order);
}

inline Boundary rational_line(const std::vector<RationalLineEntry>& entries,Boundary boundary,unsigned order=64) {
  using B=Jet::Ball;
  const auto compiled=compile_rational_entries(entries);
  std::set<std::string> denominators;std::vector<B> roots;
  for(const auto& e:entries) {
    auto p=e.coefficient.denominator();
    const auto& names=e.coefficient.variables();
    auto imaginary=std::find(names.begin(),names.end(),"I");
    if(imaginary!=names.end())p=polynomial_norm(p,imaginary-names.begin(),e.coefficient.constant(-1));
    if(p.is_rational())continue;
    p=p/p.constant(p.numerator_terms()[0].coefficient);
    if(!denominators.insert(p.str()).second)continue;
    auto parameter=std::find(names.begin(),names.end(),"x");
    if(parameter==names.end())throw std::invalid_argument("nonconstant rational-line denominator has no path variable");
    auto found=polynomial_roots(p,parameter-names.begin(),B::precision());roots.insert(roots.end(),found.begin(),found.end());
  }
  double center=0;unsigned charts=0;
  while(center<1) {
    if(++charts>20000)throw std::runtime_error("physical-basis transport chart budget exhausted");
    double next=clearance_endpoint(center,roots);
    B c,end;acb_set_d(c.raw(),center);acb_set_d(end.raw(),next);
    boundary=rational_chart(compiled,boundary,c,end-c,order);center=next;
  }
  return boundary;
}

inline data::Expr read_named_array(const std::string& path,const std::string& name) {
  std::ifstream input(path);if(!input)throw std::runtime_error("cannot read reference data: "+path);
  std::ostringstream buffer;buffer<<input.rdbuf();auto text=buffer.str();
  auto location=text.find(name+" = {");
  if(location==std::string::npos)throw std::invalid_argument("named numeric array not found: "+name);
  auto start=text.find('{',location);std::size_t depth=0,end=start;
  for(;end<text.size();++end) {if(text[end]=='{')++depth;if(text[end]=='}' && --depth==0)break;}
  if(end==text.size())throw std::invalid_argument("unterminated numeric array");
  return data::Reader(text.substr(start,end-start+1)).read();
}

inline int run_original_banana_equal(const std::string& directory) {
  using B=Jet::Ball;B::set_precision(384);
  const auto start=std::chrono::steady_clock::now();
  auto a0=data::read_file(directory+"/Data/Banana/EqualMass/dt_0.m");
  auto a1=data::read_file(directory+"/Data/Banana/EqualMass/dt_1.m");
  auto current=read_boundary(read_named_array(directory+"/BananaEqualMass.wl","boundaryAtMinusOne"),4,4,384);
  auto reference=read_boundary(read_named_array(directory+"/BananaEqualMass.wl","referenceAt20"),4,4,384);
  const std::vector<std::pair<std::string,std::string>> legs{{"-1","-1+5*I"},{"-1+5*I","20+5*I"},{"20+5*I","20"}};
  for(const auto& [from,to]:legs) {
    ExactField field({"x","I"});Exact x(field,"x"),a(field,from),b(field,to);auto line=a+(b-a)*x;
    std::vector<RationalLineEntry> entries;
    for(unsigned epsilon=0;epsilon<2;++epsilon) {
      const auto& matrix=epsilon?a1:a0;
      if(matrix.args.size()!=4)throw std::invalid_argument("banana matrix dimension");
      for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j) {
        auto c=evaluate_exact(matrix.args[i].args.at(j),x,{{"t",line}})*(b-a);
        // The exact field models I algebraically; reduce its powers in each
        // polynomial before numerical evaluation and singularity discovery.
        c=reduce_square(c.numerator(),1,x.constant(-1))/reduce_square(c.denominator(),1,x.constant(-1));
        if(!c.is_zero())entries.push_back({i,j,epsilon,std::move(c)});
      }
    }
    current=rational_line(entries,std::move(current));
  }
  double maximum=0;
  for(unsigned i=0;i<4;++i)for(unsigned k=0;k<5;++k) {
    maximum=std::max(maximum,finite_reference_error(current[i][k],reference[i][k],384));
  }
  std::cout<<"Original equal-mass banana at t=20: maximum discrepancy = "<<maximum<<", seconds = "
    <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
  std::cout<<"Homotopic upper-half-plane contour; numerical comparison, without a Taylor-tail certificate.\n";
  return std::isfinite(maximum) && maximum<1e-10?0:1;
}

inline int run_original_banana_unequal(const std::string& directory) {
  using B=Jet::Ball;B::set_precision(384);
  auto start=std::chrono::steady_clock::now();
  auto a0=data::read_file(directory+"/Data/Banana/EqualMass/dt_0.m");
  auto a1=data::read_file(directory+"/Data/Banana/EqualMass/dt_1.m");
  auto current=read_boundary(data::read_file(directory+"/Reference/BananaBoundaryAtMinusOneEps7.m"),4,7,384);
  auto reference=read_boundary(data::read_file(directory+"/Reference/BananaUnequalMassAt50.m"),15,4,384);
  for(const auto& [from,to]:std::vector<std::pair<std::string,std::string>>{{"-1","-1+10*I"},{"-1+10*I","50+10*I"},{"50+10*I","50"}}) {
    ExactField field({"x","I"});Exact x(field,"x"),a(field,from),b(field,to);auto line=a+(b-a)*x;
    std::vector<RationalLineEntry> entries;
    for(unsigned epsilon=0;epsilon<2;++epsilon)for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j) {
      auto c=evaluate_exact((epsilon?a1:a0).args.at(i).args.at(j),x,{{"t",line}})*(b-a);
      c=reduce_square(c.numerator(),1,x.constant(-1))/reduce_square(c.denominator(),1,x.constant(-1));
      if(!c.is_zero())entries.push_back({i,j,epsilon,std::move(c)});
    }
    current=rational_line(entries,std::move(current));
  }
  Boundary lifted;
  for(unsigned index:{0,0,0,0,0,0,1,1,1,1,2,3,3,3,3})lifted.push_back(current[index]);
  ExactField field({"x","I"});Exact x(field,"x"),one(field,1);
  std::map<std::string,Exact> point{{"psq",x.constant(50)},{"mm1",one+x},
    {"mm2",one+x/x.constant(2)},{"mm3",one+x/x.constant(3)},{"mm4",one}};
  std::map<std::tuple<unsigned,unsigned,unsigned>,Exact> combined;
  for(unsigned mass=1;mass<=3;++mass)for(unsigned epsilon=0;epsilon<2;++epsilon) {
    std::cout<<"Unequal banana: compiling mass "<<mass<<", epsilon degree "<<epsilon<<std::endl;
    auto matrix=data::read_file(directory+"/Data/Banana/UnequalMass/dmm"+std::to_string(mass)+"_"+std::to_string(epsilon)+".m");
    if(matrix.args.size()!=15)throw std::invalid_argument("unequal banana matrix dimension");
    for(unsigned i=0;i<15;++i)for(unsigned j=0;j<15;++j) {
      auto c=evaluate_exact(matrix.args[i].args.at(j),x,point)/x.constant(mass);
      if(c.is_zero())continue;
      auto [p,inserted]=combined.try_emplace({i,j,epsilon},x.constant(0));p->second=p->second+c;
    }
  }
  std::vector<RationalLineEntry> entries;
  for(auto& [key,c]:combined)if(!c.is_zero()) {
    auto [i,j,e]=key;entries.push_back({i,j,e,std::move(c)});
  }
  auto output=rational_line(entries,std::move(lifted));
  double maximum=0;
  for(unsigned i=0;i<15;++i)for(unsigned k=0;k<5;++k) {
    maximum=std::max(maximum,finite_reference_error(output[i][k],reference[i][k],384));
  }
  std::cout<<"Original unequal-mass banana: 15 masters through epsilon 7, comparison through epsilon 4; maximum discrepancy = "<<maximum
    <<", seconds = "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
  std::cout<<"Numerical reference comparison; Taylor remainders not yet certified.\n";
  return std::isfinite(maximum) && maximum<1e-10?0:1;
}
} // namespace diffexp
