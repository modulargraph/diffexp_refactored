#pragma once
#include "diffexp/canonical.hpp"
#include "diffexp/geometry.hpp"
#include <cmath>
#include <filesystem>
#include <limits>

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

// Heuristic size of the omitted tail of one retained Taylor polynomial: the
// last four |y_n h^n| relative to 1+|value|, as a base-2 logarithm maximised
// over all components. Ball midpoints are used so inherited input radii do not
// masquerade as slow convergence; nonfinite enclosures are left to the
// callers' own finiteness checks. This is a step-control signal, never a
// remainder certificate.
class TaylorTailEstimate {
 public:
  TaylorTailEstimate(const kernel::ComplexBall& step,unsigned order)
      :step_(log2_midpoint(step.raw())),first_(order>3?order-3:1),order_(order) {}
  // coefficient(n) returns the acb_struct pointer of this component's y_n.
  template<class Coefficient> void component(const acb_struct* value,Coefficient&& coefficient) {
    const double scale=log2_midpoint(value);
    if(scale==std::numeric_limits<double>::infinity())return;
    double terms[4],top=-std::numeric_limits<double>::infinity();unsigned count=0;
    for(unsigned n=first_;n<=order_;++n) {
      const double term=log2_midpoint(coefficient(n))+static_cast<double>(n)*step_;
      if(!std::isfinite(term))continue; // exact zero terms carry no tail; nonfinite ones no information
      terms[count++]=term;top=std::max(top,term);
    }
    if(!count)return;
    double sum=0;for(unsigned i=0;i<count;++i)sum+=std::exp2(terms[i]-top);
    const double denominator=scale>0?scale+std::log2(1+std::exp2(-scale)):std::log2(1+std::exp2(scale));
    worst_=std::max(worst_,top+std::log2(sum)-denominator);
  }
  double log2_relative()const {return worst_;}
 private:
  static double log2_midpoint(const acb_struct* value) {
    mag_t real,imaginary;mag_init(real);mag_init(imaginary);
    arf_get_mag(real,arb_midref(acb_realref(value)));arf_get_mag(imaginary,arb_midref(acb_imagref(value)));
    mag_hypot(real,real,imaginary);
    const double result=mag_is_zero(real)?-std::numeric_limits<double>::infinity():
      mag_is_inf(real)?std::numeric_limits<double>::infinity():mag_get_d_log2_approx(real);
    mag_clear(imaginary);mag_clear(real);return result;
  }
  double step_,worst_=-std::numeric_limits<double>::infinity();unsigned first_,order_;
};
// Ordinary charts step a quarter of the distance to the nearest pole, so a
// Fuchsian solution has late terms of relative size about 4^-order. A tail
// estimate far above that (the constant 2^32 absorbs polynomial growth and
// component scale separation) means the pole distance does not control this
// chart, e.g. an entire but stiff direction, and its step must be shortened.
// Terms below the working precision cannot be resolved and are not required.
inline double taylor_tail_excess(double log2_relative,unsigned order,slong bits) {
  return log2_relative-(32.0-std::min(2.0*static_cast<double>(order),static_cast<double>(bits)));
}
inline bool taylor_tail_consistent(double log2_relative,unsigned order,slong bits) {
  return !(taylor_tail_excess(log2_relative,order,bits)>0);
}

namespace rational_transport_detail {
// One Taylor-expanded coefficient and every (row,column) position sharing it.
struct ConvolutionTerm {
  unsigned epsilon;const std::pair<unsigned,unsigned>* positions;std::size_t count;std::vector<Jet::Ball> coefficients;
};
// (n+1) y_{n+1} = sum_m A_m y_{n-m} with the complete Cauchy product, followed
// by Horner evaluation at step. The physical-basis chart and the adjoint
// scalar reference are this one recurrence.
inline Boundary convolution_chart(const std::vector<ConvolutionTerm>& terms,const Boundary& initial,
    const Jet::Ball& step,unsigned order,double* truncation_log2) {
  using B=Jet::Ball;const unsigned d=initial.size(),kmax=initial[0].size()-1;const auto bits=B::precision();
  std::vector<B> values(static_cast<std::size_t>(order+1)*d*(kmax+1),B(0));
  const auto at=[&](unsigned n,unsigned i,unsigned k)->B&{return values[(static_cast<std::size_t>(n)*d+i)*(kmax+1)+k];};
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=kmax;++k)at(0,i,k)=initial[i][k];
  for(unsigned n=0;n<order;++n) {
    for(const auto& e:terms)for(std::size_t p=0;p<e.count;++p) {
      const auto [row,column]=e.positions[p];
      for(unsigned k=e.epsilon;k<=kmax;++k)for(unsigned m=0;m<=n;++m)
        if(!e.coefficients[m].is_zero() && !at(n-m,column,k-e.epsilon).is_zero())
          acb_addmul(at(n+1,row,k).raw(),e.coefficients[m].raw(),at(n-m,column,k-e.epsilon).raw(),bits);
    }
    for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=kmax;++k)
      acb_div_ui(at(n+1,i,k).raw(),at(n+1,i,k).raw(),n+1,bits);
  }
  Boundary out(d,std::vector<B>(kmax+1,B(0)));
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=kmax;++k)
    for(unsigned n=order+1;n-->0;)out[i][k]=out[i][k]*step+at(n,i,k);
  if(truncation_log2) {
    TaylorTailEstimate tail(step,order);
    for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=kmax;++k)
      tail.component(out[i][k].raw(),[&](unsigned n){return at(n,i,k).raw();});
    *truncation_log2=tail.log2_relative();
  }
  return out;
}
} // namespace rational_transport_detail

// Ordinary physical-basis evolution for an epsilon-polynomial connection. No
// homogeneous basis or spectral transform is constructed for a regular step.
// truncation_log2, when requested, receives the TaylorTailEstimate of this chart.
inline Boundary rational_chart(const std::vector<PreparedRationalEntry>& entries,
    const Boundary& boundary,const Jet::Ball& center,const Jet::Ball& step,unsigned order,
    double* truncation_log2=nullptr) {
  using B=Jet::Ball;
  if(boundary.empty() || boundary[0].empty() || !order)throw std::invalid_argument("empty rational-chart boundary or Taylor order");
  const unsigned d=boundary.size(),kmax=boundary.at(0).size()-1;
  for(const auto& row:boundary)if(row.size()!=kmax+1)throw std::invalid_argument("rational-chart boundary window mismatch");
  const auto bits=B::precision();Jet x(0,order+1,bits);x.set(0,center);x.set(1,B(1));
  std::vector<std::pair<unsigned,unsigned>> positions;positions.reserve(entries.size());
  std::vector<rational_transport_detail::ConvolutionTerm> matrix;matrix.reserve(entries.size());
  auto imaginary=x.constant(0);imaginary.set(0,B::from_strings("0","1"));
  for(const auto& entry:entries) {
    if(entry.row>=d || entry.column>=d)throw std::invalid_argument("rational matrix index exceeds boundary size");
    auto jet=evaluate(entry.coefficient,x,{{"x",x},{"I",imaginary}});
    positions.emplace_back(entry.row,entry.column);
    rational_transport_detail::ConvolutionTerm e{entry.epsilon,&positions.back(),1,{}};
    for(unsigned n=0;n<order;++n)e.coefficients.push_back(jet.at(n));matrix.push_back(std::move(e));
  }
  return rational_transport_detail::convolution_chart(matrix,boundary,step,order,truncation_log2);
}
inline Boundary rational_chart(const std::vector<RationalLineEntry>& entries,
    const Boundary& boundary,const Jet::Ball& center,const Jet::Ball& step,unsigned order,
    double* truncation_log2=nullptr) {
  return rational_chart(compile_rational_entries(entries),boundary,center,step,order,truncation_log2);
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
    B c;acb_set_d(c.raw(),center);
    // The geometric step only bounds late terms for pole-controlled growth.
    // Reject a chart whose retained tail contradicts it and halve its step.
    for(unsigned halvings=0;;++halvings) {
      B end;acb_set_d(end.raw(),next);double tail=-std::numeric_limits<double>::infinity();
      auto candidate=rational_chart(compiled,boundary,c,end-c,order,&tail);
      if(taylor_tail_consistent(tail,order,B::precision())){boundary=std::move(candidate);break;}
      const double smaller=center+(next-center)/2;
      if(halvings>=40 || !(smaller>center) || !(smaller<next))
        throw std::runtime_error("physical-basis transport Taylor truncation step budget exhausted");
      next=smaller;
    }
    center=next;
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

inline data::Expr original_banana_reference(const std::string& directory,const std::string& name,unsigned epsilon_order=4) {
  const auto path=directory+"/Reference/"+(name=="boundaryAtMinusOne"?"BananaBoundaryAtMinusOneEps7.m":"BananaEqualMassAt20.m");
  if(std::filesystem::exists(path)) {
    auto data=data::read_file(path);for(auto& row:data.args){if(row.args.size()<epsilon_order+1)throw std::invalid_argument("banana reference epsilon window");row.args.resize(epsilon_order+1);}return data;
  }
  return read_named_array(directory+"/BananaEqualMass.wl",name);
}

inline int run_original_banana_equal(const std::string& directory) {
  using B=Jet::Ball;B::set_precision(384);
  const auto start=std::chrono::steady_clock::now();
  auto a0=data::read_file(directory+"/Data/Banana/EqualMass/dt_0.m");
  auto a1=data::read_file(directory+"/Data/Banana/EqualMass/dt_1.m");
  auto current=read_boundary(original_banana_reference(directory,"boundaryAtMinusOne"),4,4,384);
  auto reference=read_boundary(original_banana_reference(directory,"referenceAt20"),4,4,384);
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
