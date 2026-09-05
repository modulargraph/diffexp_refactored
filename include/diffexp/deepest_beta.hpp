#pragma once
#include "diffexp/feynman.hpp"
#include "diffexp/tadpole.hpp"
#include "diffexp/rational_transport.hpp"

namespace diffexp::feynman {
struct LaurentValue {
  int epsilon_low;
  std::vector<Jet::Ball> coefficients;
  bool taylor_tail_certified=false;
  Jet::Ball at(int order) const {
    if(order<epsilon_low || static_cast<std::size_t>(order-epsilon_low)>=coefficients.size())
      throw std::out_of_range("Laurent value outside computed window");
    return coefficients[order-epsilon_low];
  }
};
struct OriginFactor {long power;Exact regular;Rational value;};
inline OriginFactor factor_origin(const Exact& expression,std::size_t variable) {
  if(expression.is_zero())throw std::invalid_argument("identically zero deepest Symanzik polynomial");
  auto valuation=[&](const std::vector<Exact::Term>& terms) {
    auto lowest=std::numeric_limits<unsigned long>::max();
    for(const auto& term:terms)lowest=std::min(lowest,term.powers.at(variable));
    if(lowest>1000000)throw std::invalid_argument("origin valuation exceeds resource limit");
    return static_cast<long>(lowest);
  };
  const long power=valuation(expression.numerator_terms())-valuation(expression.denominator_terms());
  auto x=expression.variable(variable);
  auto regular=power<0?expression*x.pow(-power):expression/x.pow(power);
  std::vector<Exact> origin;
  for(std::size_t i=0;i<expression.variable_count();++i)origin.push_back(i==variable?expression.constant(0):expression.variable(i));
  return {power,regular,regular.substitute(origin).rational()};
}

// Integrate one beta functional of the deepest scalar boundary. Its local
// x^(a+b eps) factor is kept symbolic through the primitive, including the
// epsilon pole when a+n=-1. Only the ordinary analytic factor is Taylor
// expanded. This implements dimensional analytic continuation at the lower
// endpoint; the upper endpoint must be regular. It returns a numerical series
// approximation until its omitted Taylor tail has a certificate.
inline LaurentValue deepest_beta(const Symanzik& geometry,unsigned loops,int d0,
    unsigned left,unsigned right,unsigned epsilon_top,unsigned taylor_order=80) {
  using B=Jet::Ball;
  if(!loops || loops>64 || d0< -64 || d0>64 || d0%2 || !left || !right || left>10000 || right>10000 || taylor_order<8 || taylor_order>1000 || epsilon_top>100)
    throw std::invalid_argument("invalid deepest beta request");
  auto request=merge_request({static_cast<int>(left),static_cast<int>(right)},0,1);
  const auto& names=geometry.U.variables();auto found=std::find(names.begin(),names.end(),"x");
  if(found==names.end())throw std::invalid_argument("deepest beta needs the path variable x");
  const auto variable=static_cast<std::size_t>(found-names.begin());
  auto uf=factor_origin(geometry.U,variable),ff=factor_origin(geometry.F,variable);
  const long v=left+right,a=v-static_cast<long>(loops)*d0/2,c=v-(static_cast<long>(loops)+1)*d0/2;
  const long q=c*uf.power-a*ff.power+left-1;
  const long slope=(static_cast<long>(loops)+1)*uf.power-static_cast<long>(loops)*ff.power;
  if(q<=-1 && -q-1>static_cast<long>(taylor_order))
    throw std::invalid_argument("Taylor order does not reach the resonant primitive coefficient");
  const unsigned primitive_poles=q<=-1 && slope!=0?1:0;
  auto gamma=tadpole(uf.value,ff.value,v,loops,d0,epsilon_top+primitive_poles);
  const int low=gamma.epsilon_low-static_cast<int>(primitive_poles),top=static_cast<int>(epsilon_top);
  const auto bits=B::precision();Jet x(0,taylor_order+1,bits);x.set(1,B(1));
  // Rational constants are evaluated as expressions, not decimal literals.
  auto constant=[&](const Rational& value){return evaluate(data::Reader(value.str()).read(),x,{});};
  auto ur=evaluate(data::Reader(uf.regular.str()).read(),x,{{"x",x}})/constant(uf.value);
  auto fr=evaluate(data::Reader(ff.regular.str()).read(),x,{{"x",x}})/constant(ff.value);
  auto prefactor=ur.pow(c)*fr.pow(-a)*(x.constant(1)-x).pow(right-1)*constant(request.normalization);
  auto logarithm=x.constant(loops+1)*ur.log()-x.constant(loops)*fr.log();
  const auto maximum=static_cast<unsigned>(top+primitive_poles-gamma.epsilon_low);
  std::vector<Jet> analytic(maximum+1,x.constant(0));auto log_power=x.constant(1);
  for(unsigned m=0;m<=maximum;++m) {
    for(unsigned j=0;j+m<=maximum && j<gamma.coefficients.size();++j) {
      auto coefficient=x.constant(0);coefficient.set(0,gamma.coefficients[j]);
      analytic[j+m]=analytic[j+m]+prefactor*log_power*coefficient;
    }
    log_power=log_power*logarithm/x.constant(m+1);
  }
  // Choose a dyadic overlap inside the proved root clearance of the regular
  // analytic factors. This establishes convergence, not a bound on the tail.
  std::vector<B> singularities;
  for(const auto& p:{uf.regular.numerator(),uf.regular.denominator(),ff.regular.numerator(),ff.regular.denominator()}) {
    if(p.is_rational())continue;
    auto roots=polynomial_roots(p,variable,bits);singularities.insert(singularities.end(),roots.begin(),roots.end());
  }
  const double allowed=clearance_endpoint(0,singularities);long divisor=8;
  while(1.0/divisor>allowed) {
    if(divisor>(1L<<50))throw std::runtime_error("deepest overlap requires a dedicated endpoint rescaling");
    divisor*=2;
  }
  const Rational overlap=Rational(1)/Rational(divisor);const auto h=B::from_strings(overlap.str());
  std::vector<B> primitive(top-low+1,B(0));B h_power(1);
  for(unsigned n=0;n<=taylor_order;++n) {
    const long denominator=static_cast<long>(n)+q+1;
    for(unsigned k=0;k<analytic.size();++k) {
      auto coefficient=analytic[k].at(n)*h_power;const int exponent=gamma.epsilon_low+static_cast<int>(k);
      if(coefficient.is_zero())continue;
      if(!denominator) {
        if(!slope)throw std::domain_error("deepest endpoint has a logarithmic divergence not regulated by epsilon");
        if(exponent-1>=low && exponent-1<=top)primitive[exponent-1-low]+=coefficient/B(slope);
      } else {
        auto term=coefficient/B(denominator);
        for(int target=exponent;target<=top;++target) {
          if(target>=low)primitive[target-low]+=term;
          term=term*B(-slope)/B(denominator);
        }
      }
    }
    h_power=h_power*h;
  }
  B log_h,scale;acb_log(log_h.raw(),h.raw(),bits);acb_pow_si(scale.raw(),h.raw(),q+1,bits);
  Jet exponential(0,top-low+1,bits);if(exponential.length()>1)exponential.set(1,B(slope)*log_h);
  exponential=exponential.exp();std::vector<B> matched(top-low+1,B(0));
  for(unsigned k=0;k<matched.size();++k)for(unsigned j=0;j<=k;++j)matched[k]+=scale*primitive[k-j]*exponential.at(j);
  std::vector<Exact> point;const auto sample=geometry.U;
  for(std::size_t i=0;i<sample.variable_count();++i)point.push_back(sample.variable(i));
  point[variable]=sample.constant(overlap);
  auto seed=tadpole(geometry.U.substitute(point).rational(),geometry.F.substitute(point).rational(),v,loops,d0,epsilon_top);
  Boundary boundary{std::vector<B>(top-low+1,B(0)),std::move(matched)};
  for(int k=std::max(seed.epsilon_low,low);k<=top;++k)boundary[0][k-low]=seed.coefficients[k-seed.epsilon_low];
  auto parameter=sample.variable(variable);
  auto a0=sample.constant(c)*geometry.U.derivative(variable)/geometry.U-sample.constant(a)*geometry.F.derivative(variable)/geometry.F;
  auto a1=sample.constant(loops+1)*geometry.U.derivative(variable)/geometry.U-sample.constant(loops)*geometry.F.derivative(variable)/geometry.F;
  auto weight=sample.constant(request.normalization)*parameter.pow(left-1)*(sample.constant(1)-parameter).pow(right-1);
  const auto width=sample.constant(Rational(1)-overlap);
  point[variable]=sample.constant(overlap)+width*parameter;
  // A regular terminal point is a deliberate contract of this scalar path.
  auto endpoint=point;endpoint[variable]=sample.constant(1);
  if(geometry.U.substitute(endpoint).is_zero() || geometry.F.substitute(endpoint).is_zero())
    throw std::domain_error("deepest beta needs a separate singular upper-endpoint operation");
  std::vector<RationalLineEntry> entries{{0,0,0,width*a0.substitute(point)},
    {0,0,1,width*a1.substitute(point)},{1,0,0,width*weight.substitute(point)}};
  auto output=rational_line(entries,std::move(boundary),taylor_order);
  return {low,std::move(output[1]),false};
}
} // namespace diffexp::feynman
