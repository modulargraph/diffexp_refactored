#pragma once
#include "diffexp/certified_deepest_beta.hpp"
#include "diffexp/epsilon_gauge.hpp"
#include "diffexp/certified_rational_transport.hpp"

namespace diffexp::feynman {
namespace scalar_functional_detail {
struct Weight {
  int x_power,epsilon_power;
  std::optional<std::size_t> epsilon_variable;
  Exact regular;
};
inline Weight prepare_weight(const Exact& weight,std::size_t x_variable) {
  if(weight.is_zero())return {0,0,std::nullopt,weight};
  std::optional<std::size_t> epsilon;
  for(std::size_t i=0;i<weight.variables().size();++i) {
    const auto& name=weight.variables()[i];
    if(name=="eps")epsilon=i;
    else if(name!="x" && name!="I") {
      for(const auto& terms:{weight.numerator_terms(),weight.denominator_terms()})for(const auto& term:terms)
        if(term.powers[i])throw std::invalid_argument("unbound scalar functional weight parameter: "+name);
    }
  }
  const auto valuation=[&](std::size_t variable) {
    unsigned long n=100001,d=100001;
    for(const auto& term:weight.numerator_terms())n=std::min(n,term.powers[variable]);
    for(const auto& term:weight.denominator_terms())d=std::min(d,term.powers[variable]);
    if(n>100000 || d>100000)throw std::invalid_argument("scalar weight valuation budget");
    return static_cast<int>(n)-static_cast<int>(d);
  };
  int xp=valuation(x_variable),ep=epsilon?valuation(*epsilon):0;
  if(std::abs(ep)>100)throw std::invalid_argument("scalar weight epsilon valuation budget");
  auto regular=weight;
  const auto deflate=[&](std::size_t variable,int power){if(power<0)regular=regular*weight.variable(variable).pow(-power);else regular=regular/weight.variable(variable).pow(power);};
  deflate(x_variable,xp);if(epsilon)deflate(*epsilon,ep);
  std::vector<Exact> origin;for(std::size_t i=0;i<weight.variable_count();++i)origin.push_back(weight.variable(i));
  origin[x_variable]=weight.constant(0);if(epsilon)origin[*epsilon]=weight.constant(0);
  auto denominator=regular.denominator().substitute(origin);
  auto origin_ball=native_tail_detail::eval(data::Reader(denominator.str()).read(),kernel::ComplexBall(0),"x").value;
  if(origin_ball.contains_zero() || !origin_ball.is_finite())
    throw std::domain_error("scalar weight is not jointly analytic after x/epsilon monomial deflation");
  origin[x_variable]=weight.constant(1);
  auto upper=regular.denominator().substitute(origin);
  auto upper_ball=native_tail_detail::eval(data::Reader(upper.str()).read(),kernel::ComplexBall(0),"x").value;
  if(upper_ball.contains_zero() || !upper_ball.is_finite())
    throw std::domain_error("scalar functional requires a regular upper-endpoint weight");
  return {xp,ep,epsilon,std::move(regular)};
}
inline std::vector<Exact> epsilon_series(const Exact& rational,std::optional<std::size_t> epsilon,unsigned top) {
  std::vector<Exact> result(top+1,rational.constant(0));
  if(!epsilon){result[0]=rational;return result;}
  const auto polynomial=[&](const std::vector<Exact::Term>& terms) {
    std::vector<Exact> out(top+1,rational.constant(0));
    for(const auto& term:terms) {
      auto n=term.powers[*epsilon];if(n>top)continue;
      auto value=rational.constant(term.coefficient);
      for(std::size_t i=0;i<term.powers.size();++i)if(i!=*epsilon && term.powers[i])value=value*rational.variable(i).pow(term.powers[i]);
      out[n]=out[n]+value;
    }
    return out;
  };
  auto numerator=polynomial(rational.numerator_terms()),denominator=polynomial(rational.denominator_terms());
  if(denominator[0].is_zero())throw std::domain_error("scalar weight epsilon denominator has no regular constant");
  for(unsigned n=0;n<=top;++n) {
    auto value=numerator[n];for(unsigned j=1;j<=n;++j)value=value-denominator[j]*result[n-j];
    result[n]=value/denominator[0];
  }
  return result;
}
} // namespace scalar_functional_detail

// Integrate weight(x,eps)*J_power(x,eps) with exact rational weight combined
// before deflation/expansion. The retained affine-epsilon endpoint power is
// integrated meromorphically; the omitted analytic tail is integrated with
// absolute log-moment bounds. Negative weight epsilon valuation shifts the
// source and requests gamma lookahead. No quadrature enters the solver.
// Geometry and weight must belong to the same ExactField, with geometry
// independent of eps and other parameters specialized. The lower weight must
// be jointly analytic after monomial deflation; the upper endpoint is regular.
namespace scalar_functional_detail {
inline LaurentValue integrate_one_sided(const Symanzik& geometry,unsigned loops,int d0,
    unsigned power,const Exact& weight,unsigned epsilon_top,const CertifiedDeepestBetaOptions& options={}) {
  const unsigned taylor_order=options.taylor_order;
  using B=Jet::Ball;
  if(!loops || loops>64 || d0< -64 || d0>64 || d0%2 || !power || power>20000 || taylor_order<8 || taylor_order>1000 || epsilon_top>100)
    throw std::invalid_argument("invalid scalar functional request");
  if(options.working_bits<64 || options.working_bits>1000000 || !options.requested_digits || options.requested_digits>10000 || !options.max_charts || !options.max_disk_attempts || (options.f_rim!=1 && options.f_rim!=-1))
    throw std::invalid_argument("invalid certified scalar functional precision or proof budget");
  const Rational requested_radius(options.initial_radius);
  if(!(Rational(0)<requested_radius))throw std::invalid_argument("certified scalar functional radius must be positive");
  struct PrecisionScope {slong previous;explicit PrecisionScope(slong bits):previous(B::precision()){B::set_precision(bits);}~PrecisionScope(){B::set_precision(previous);}} precision(options.working_bits);

  const auto& names=geometry.U.variables();auto found=std::find(names.begin(),names.end(),"x");
  if(found==names.end())throw std::invalid_argument("scalar functional needs the path variable x");
  const auto variable=static_cast<std::size_t>(found-names.begin());
  (void)(geometry.U+geometry.F+weight); // enforce a single owned exact field
  for(std::size_t i=0;i<names.size();++i)if(names[i]=="eps" && (!geometry.U.derivative(i).is_zero() || !geometry.F.derivative(i).is_zero()))
    throw std::invalid_argument("scalar functional geometry must be epsilon independent");
  if(weight.is_zero())return {0,std::vector<B>(epsilon_top+1,B(0)),true};
  auto uf=factor_origin(geometry.U,variable),ff=factor_origin(geometry.F,variable);
  const long v=power,a=v-static_cast<long>(loops)*d0/2,c=v-(static_cast<long>(loops)+1)*d0/2;
  auto weights=scalar_functional_detail::prepare_weight(weight,variable);
  const long q=c*uf.power-a*ff.power+weights.x_power;
  const long slope=(static_cast<long>(loops)+1)*uf.power-static_cast<long>(loops)*ff.power;
  if(q<=-1 && -q-1>static_cast<long>(taylor_order))
    throw std::invalid_argument("Taylor order does not reach the resonant primitive coefficient");
  const unsigned primitive_poles=q<=-1 && slope!=0?1:0;
  const unsigned gamma_top=static_cast<unsigned>(std::max<long>(0,static_cast<long>(epsilon_top)-weights.epsilon_power+primitive_poles));
  auto gamma=tadpole(uf.value,ff.value,v,loops,d0,gamma_top,options.f_rim);
  const int source_low=gamma.epsilon_low+weights.epsilon_power;
  const int low=std::min(0,source_low-static_cast<int>(primitive_poles)),top=static_cast<int>(epsilon_top);
  const unsigned maximum=static_cast<unsigned>(std::max(0,top+static_cast<int>(primitive_poles)-source_low));
  const auto weight_coefficients=scalar_functional_detail::epsilon_series(weights.regular,weights.epsilon_variable,maximum);
  const auto bits=B::precision();Jet x(0,taylor_order+1,bits);x.set(1,B(1));
  // Rational constants are evaluated as expressions, not decimal literals.
  auto constant=[&](const Rational& value){return evaluate(data::Reader(value.str()).read(),x,{});};
  auto ur=evaluate(data::Reader(uf.regular.str()).read(),x,{{"x",x}})/constant(uf.value);
  auto fr=evaluate(data::Reader(ff.regular.str()).read(),x,{{"x",x}})/constant(ff.value);
  auto prefactor=ur.pow(c)*fr.pow(-a);
  auto logarithm=x.constant(loops+1)*ur.log()-x.constant(loops)*fr.log();
  std::vector<Jet> weight_jets;
  for(const auto& coefficient:weight_coefficients)
    weight_jets.push_back(evaluate(data::Reader(coefficient.str()).read(),x,{{"x",x}}));
  std::vector<Jet> analytic(maximum+1,x.constant(0));auto log_power=x.constant(1);
  for(unsigned m=0;m<=maximum;++m) {
    for(unsigned j=0;j+m<=maximum && j<gamma.coefficients.size();++j) {
      auto coefficient=x.constant(0);coefficient.set(0,gamma.coefficients[j]);
      for(unsigned w=0;w+j+m<=maximum;++w)
        analytic[w+j+m]=analytic[w+j+m]+prefactor*log_power*coefficient*weight_jets[w];
    }
    log_power=log_power*logarithm/x.constant(m+1);
  }
  // Whole-disk proof of the ordinary analytic factor, independent of its
  // stored Taylor coefficients. Principal logs apply to ratios equal to 1
  // at the origin, so the explicit endpoint x^(q+slope*epsilon) is untouched.
  using M=NativeTailMagnitude;
  auto ur_exact=uf.regular/uf.regular.constant(uf.value);
  auto fr_exact=ff.regular/ff.regular.constant(ff.value);
  const auto signed_power=[](const Exact& value,long power){return power<0?value.constant(1)/value.pow(-power):value.pow(power);};
  auto p_exact=signed_power(ur_exact,c)*signed_power(fr_exact,-a);
  auto p_expression=data::Reader(p_exact.str()).read();
  auto l_expression=data::Reader(std::to_string(loops+1)+"*Log["+ur_exact.str()+"]-"+std::to_string(loops)+"*Log["+fr_exact.str()+"]").read();
  Rational radius=requested_radius;M p_bound,l_bound;std::vector<M> weight_bounds;bool disk_proved=false;
  for(unsigned attempt=0;attempt<options.max_disk_attempts;++attempt) {
    B disk(0);M::upper_abs(B::from_strings(radius.str())).add_error_to(disk);
    try {
      auto p=native_tail_detail::eval(p_expression,disk,"x").value;
      auto l=native_tail_detail::eval(l_expression,disk,"x").value;
      native_tail_detail::require_finite(p);native_tail_detail::require_finite(l);
      std::vector<M> candidate_weights;
      for(const auto& coefficient:weight_coefficients) {
        auto value=native_tail_detail::eval(data::Reader(coefficient.str()).read(),disk,"x").value;
        native_tail_detail::require_finite(value);candidate_weights.push_back(M::upper_abs(value));
      }
      p_bound=M::upper_abs(p);l_bound=M::upper_abs(l);weight_bounds=std::move(candidate_weights);disk_proved=true;break;
    }catch(const std::domain_error&) {radius=radius/Rational(2);}
  }
  if(!disk_proved)throw std::runtime_error("certified scalar functional endpoint disk-proof budget exhausted");
  const Rational quarter=radius/Rational(4);
  const Rational overlap=quarter<Rational("1/8")?quarter:Rational("1/8");
  const auto h=B::from_strings(overlap.str());
  const long first_denominator=static_cast<long>(taylor_order)+q+2;
  if(first_denominator<=0)throw std::domain_error("certified scalar functional omitted tail is not integrable");
  std::vector<M> analytic_bounds(maximum+1);M logarithm_factor=M::one();
  for(unsigned m=0;m<=maximum;++m) {
    for(unsigned j=0;j+m<=maximum && j<gamma.coefficients.size();++j)
      for(unsigned w=0;w+j+m<=maximum;++w)analytic_bounds[w+j+m]+=p_bound*logarithm_factor*M::upper_abs(gamma.coefficients[j])*weight_bounds[w];
    logarithm_factor=logarithm_factor*l_bound/M::lower_abs(B::from_strings(std::to_string(m+1)));
  }
  std::vector<B> primitive(top-low+1,B(0));B h_power(1);
  for(unsigned n=0;n<=taylor_order;++n) {
    const long denominator=static_cast<long>(n)+q+1;
    for(unsigned k=0;k<analytic.size();++k) {
      auto coefficient=analytic[k].at(n)*h_power;const int exponent=source_low+static_cast<int>(k);
      if(coefficient.is_zero())continue;
      if(!denominator) {
        if(!slope)throw std::domain_error("scalar functional endpoint has a logarithmic divergence not regulated by epsilon");
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
  // The unseen powers all have n+q>-1; integrate their log moments
  // absolutely. Unlike the retained resonant primitive they have no DR pole.
  // For D=N+q+2 and H=-log(h), normalized moments obey
  // T_0=1/D, T_j=H^j/(j! D)+T_(j-1)/D.
  const auto radius_lower=M::lower_abs(B::from_strings(radius.str()));
  const auto ratio=M::upper_abs(h)/radius_lower;
  const auto gap=M::positive_difference_lower(M::one(),ratio);
  if(gap.is_zero())throw std::domain_error("certified scalar functional overlap outside endpoint disk");
  const auto common=M::upper_abs(scale)*ratio.power_upper(taylor_order+1)/gap;
  const auto d_lower=M::lower_abs(B::from_strings(std::to_string(first_denominator)));
  const auto h_log=M::upper_abs(log_h);
  const auto slope_abs=M::upper_abs(B(slope));
  const unsigned tail_depth=static_cast<unsigned>(std::max(0,top-source_low));
  std::vector<M> moments(tail_depth+1);M hlog_factor=M::one(),slope_power=M::one();
  for(unsigned j=0;j<=tail_depth;++j) {
    moments[j]=hlog_factor/d_lower;
    if(j)moments[j]+=moments[j-1]/d_lower;
    hlog_factor=hlog_factor*h_log/M::lower_abs(B::from_strings(std::to_string(j+1)));
  }
  for(unsigned j=0;j<=tail_depth;++j) {
    for(unsigned k=0;k+j<=tail_depth;++k) {
      auto error=common*analytic_bounds[k]*slope_power*moments[j];
      if(!error.is_finite())throw std::domain_error("nonfinite certified endpoint primitive tail");
      if(source_low+static_cast<int>(k+j)<=top)error.add_error_to(matched[source_low+static_cast<int>(k+j)-low]);
    }
    slope_power=slope_power*slope_abs;
  }
  std::vector<Exact> point;const auto sample=geometry.U;
  for(std::size_t i=0;i<sample.variable_count();++i)point.push_back(sample.variable(i));
  point[variable]=sample.constant(overlap);
  auto seed=tadpole(geometry.U.substitute(point).rational(),geometry.F.substitute(point).rational(),v,loops,d0,static_cast<unsigned>(std::max<long>(0,static_cast<long>(epsilon_top)-weights.epsilon_power)),options.f_rim);
  Boundary boundary{std::vector<B>(top-low+1,B(0)),std::move(matched)};
  for(int k=std::max(seed.epsilon_low+weights.epsilon_power,low);k<=top;++k)boundary[0][k-low]=seed.coefficients[k-weights.epsilon_power-seed.epsilon_low];
  auto parameter=sample.variable(variable);
  auto a0=sample.constant(c)*geometry.U.derivative(variable)/geometry.U-sample.constant(a)*geometry.F.derivative(variable)/geometry.F;
  auto a1=sample.constant(loops+1)*geometry.U.derivative(variable)/geometry.U-sample.constant(loops)*geometry.F.derivative(variable)/geometry.F;
  const auto width=sample.constant(Rational(1)-overlap);
  point[variable]=sample.constant(overlap)+width*parameter;
  // A regular terminal point is a deliberate contract of this scalar path.
  auto endpoint=point;endpoint[variable]=sample.constant(1);
  if(geometry.U.substitute(endpoint).is_zero() || geometry.F.substitute(endpoint).is_zero())
    throw std::domain_error("scalar functional needs a separate singular upper-endpoint operation");
  std::vector<RationalLineEntry> entries{{0,0,0,width*a0.substitute(point)},
    {0,0,1,width*a1.substitute(point)}};
  for(unsigned e=0;e<weight_coefficients.size();++e) {
    auto regular_weight=signed_power(parameter,weights.x_power)*weight_coefficients[e];
    if(!regular_weight.is_zero())entries.push_back({1,0,e,width*regular_weight.substitute(point)});
  }
  CertifiedRationalOptions ordinary;
  ordinary.working_bits=options.working_bits;ordinary.requested_digits=options.requested_digits;
  ordinary.taylor_order=taylor_order;ordinary.max_charts=options.max_charts;
  ordinary.max_disk_attempts=options.max_disk_attempts;ordinary.initial_radius=options.initial_radius;
  auto result=certified_rational_line(entries,boundary,ordinary);
  return {low,std::move(result.boundary[1]),true};
}
} // namespace scalar_functional_detail

// A singular upper endpoint is handled by reversing the second half of the
// same scalar functional. Both lower-endpoint primitives must independently
// certify their full omitted tails; their final sum is checked again. This is
// the same meromorphic DR endpoint prescription on each arm, with an explicit
// common F rim. No arbitrary real-interval integration across an interior zero
// of U/F is admitted by the ordinary disk proofs.
inline LaurentValue certified_scalar_functional(const Symanzik& geometry,unsigned loops,int d0,
    unsigned power,const Exact& weight,unsigned epsilon_top,const CertifiedDeepestBetaOptions& options={}) {
  using B=kernel::ComplexBall;
  if(options.working_bits<64 || options.working_bits>1000000 || !options.requested_digits || options.requested_digits>10000 || (options.f_rim!=1 && options.f_rim!=-1))
    throw std::invalid_argument("invalid scalar functional precision or F rim");
  struct PrecisionScope {slong previous;explicit PrecisionScope(slong bits):previous(B::precision()){B::set_precision(bits);}~PrecisionScope(){B::set_precision(previous);}} precision(options.working_bits);
  (void)(geometry.U.constant(0)+geometry.F+weight);
  const auto& names=weight.variables();auto xi=std::find(names.begin(),names.end(),"x"),ei=std::find(names.begin(),names.end(),"eps");
  if(xi==names.end())throw std::invalid_argument("scalar functional requires x");
  const auto variable=static_cast<std::size_t>(xi-names.begin());auto x=weight.variable(variable);
  std::vector<Exact> upper;for(std::size_t i=0;i<weight.variable_count();++i)upper.push_back(weight.variable(i));
  upper[variable]=weight.constant(1);
  bool split=false;
  try {split=geometry.U.substitute(upper).is_zero() || geometry.F.substitute(upper).is_zero();}
  catch(const std::domain_error&){split=true;}
  auto denominator=weight.denominator();
  if(ei!=names.end()) {
    const auto epsilon=static_cast<std::size_t>(ei-names.begin());
    auto valuation=*exact_epsilon_valuation(denominator,epsilon);
    if(valuation)denominator=denominator/weight.variable(epsilon).pow(valuation);
    upper[epsilon]=weight.constant(0);
  }
  auto upper_denominator=denominator.substitute(upper);
  auto separation=native_tail_detail::eval(data::Reader(upper_denominator.str()).read(),B(0),"x").value;
  if(separation.contains_zero() || !separation.is_finite())split=true;
  if(!split)return scalar_functional_detail::integrate_one_sided(geometry,loops,d0,power,weight,epsilon_top,options);
  const Rational anchor(options.scalar_split_anchor);
  if(anchor<=Rational(0) || anchor>=Rational(1) || options.requested_digits>=10000)
    throw std::invalid_argument("invalid scalar two-endpoint split anchor or precision budget");
  auto inner=options;++inner.requested_digits;
  std::vector<LaurentValue> halves;
  for(bool reverse:{false,true}) {
    const Rational length=reverse?Rational(1)-anchor:anchor;
    std::vector<Exact> point;for(std::size_t i=0;i<weight.variable_count();++i)point.push_back(weight.variable(i));
    point[variable]=reverse?weight.constant(1)-weight.constant(length)*x:weight.constant(length)*x;
    Symanzik mapped{geometry.U.substitute(point),geometry.F.substitute(point)};
    auto mapped_weight=weight.constant(length)*weight.substitute(point);
    halves.push_back(scalar_functional_detail::integrate_one_sided(mapped,loops,d0,power,mapped_weight,epsilon_top,inner));
  }
  const int low=std::min(halves[0].epsilon_low,halves[1].epsilon_low);
  LaurentValue result{low,std::vector<B>(static_cast<int>(epsilon_top)-low+1,B(0)),false};
  for(const auto& half:halves) {
    if(!half.taylor_tail_certified)throw std::logic_error("uncertified arm in scalar two-endpoint functional");
    for(unsigned k=0;k<half.coefficients.size();++k)result.coefficients[half.epsilon_low+static_cast<int>(k)-low]+=half.coefficients[k];
  }
  const std::string tolerance="1/1"+std::string(options.requested_digits,'0');
  if(!native_enclosure_meets_tolerance({result.coefficients},tolerance))throw std::runtime_error("scalar two-endpoint sum does not meet requested absolute digits");
  result.taylor_tail_certified=true;return result;
}
} // namespace diffexp::feynman
