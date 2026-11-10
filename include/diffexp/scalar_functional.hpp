#pragma once
#include "diffexp/deepest_beta.hpp"
#include "diffexp/certified_tail.hpp"
#include "diffexp/epsilon_gauge.hpp"
#include <optional>

namespace diffexp::feynman {
struct CertifiedDeepestBetaOptions {
  slong working_bits=256;
  unsigned requested_digits=20;
  unsigned taylor_order=112;
  unsigned max_charts=512;
  unsigned max_disk_attempts=32;
  std::string initial_radius="1/4";
  int f_rim=1;
  std::string scalar_split_anchor="1/2";
};
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
using B=Jet::Ball;using M=NativeTailMagnitude;
// A rational function of x alone, kept as its two polynomials.
struct Fraction {data::Expr numerator,denominator;unsigned numerator_degree,denominator_degree;};
inline Fraction compile_fraction(const Exact& value,std::size_t variable) {
  const auto degree=[&](const std::vector<Exact::Term>& terms) {
    unsigned long d=0;for(const auto& term:terms)d=std::max(d,term.powers.at(variable));
    if(d>100000)throw std::invalid_argument("scalar functional polynomial degree budget");
    return static_cast<unsigned>(d);
  };
  return {data::Reader(value.numerator().str()).read(),data::Reader(value.denominator().str()).read(),
    degree(value.numerator_terms()),degree(value.denominator_terms())};
}
// Both polynomials Taylor-shifted to a chart centre. The shift is exact
// polynomial arithmetic on balls, so evaluating the shifted coefficients on
// |z|<=R encloses the function on the whole disk without the dependency of an
// unshifted evaluation, and the same coefficients give the quotient's series.
struct Local {
  Jet numerator,denominator;B centre_value,enclosure;
  Jet series(unsigned length)const {
    const auto resized=[&](const Jet& p){Jet out(0,length,p.bits());for(unsigned n=0;n<std::min(length,p.length());++n)out.set(n,p.at(n));return out;};
    return resized(numerator)/resized(denominator);
  }
};
inline Local localize(const Fraction& f,const B& centre,const B& disk,slong bits) {
  const auto shifted=[&](const data::Expr& polynomial,unsigned degree) {
    Jet x(0,std::max(degree+1,2u),bits);x.set(0,centre);x.set(1,B(1));return evaluate(polynomial,x,{{"x",x}});
  };
  auto numerator=shifted(f.numerator,f.numerator_degree),denominator=shifted(f.denominator,f.denominator_degree);
  const auto bound=denominator.evaluate_polynomial(disk);
  if(bound.contains_zero())throw std::domain_error("scalar functional disk meets a pole");
  auto enclosure=numerator.evaluate_polynomial(disk)/bound;native_tail_detail::require_finite(enclosure);
  auto centre_value=numerator.at(0)/denominator.at(0);
  return {std::move(numerator),std::move(denominator),std::move(centre_value),std::move(enclosure)};
}
inline B ball_power(const B& value,long power) {
  if(power<0 && value.contains_zero())throw std::domain_error("scalar functional disk meets a zero");
  B out;acb_pow_si(out.raw(),value.raw(),power,B::precision());native_tail_detail::require_finite(out);return out;
}
inline B ball_log(const B& value) {
  if(!native_tail_detail::avoids_principal_cut(value))throw std::domain_error("scalar functional disk leaves the principal logarithmic sheet");
  B out;acb_log(out.raw(),value.raw(),B::precision());native_tail_detail::require_finite(out);return out;
}
// With U=x^p*u, F=x^r*f and weight=eps^s*x^t*W the integrand about a centre
// c is, exactly,
//   T_c(eps) * c^t * X^(q+slope*eps) * ur^C * fr^(-A) * exp(eps*L) * W,
// where T_c is the gamma boundary at c (at the origin: of the regular parts),
// X=(c+z)/c, ur=u(c+z)/u(c), fr=f(c+z)/f(c), L=(loops+1)*log(ur)-loops*log(fr)
// and q=C*p-A*r+t. Integer C,A and principal logarithms of ratios equal to one
// at the centre make every factor single valued on a disk where the ratios
// avoid the cut; on the real segment through c this is the selected F rim.
// At the origin X^(q+slope*eps) is kept symbolic by the caller.
struct Integrand {
  Fraction u,f;std::vector<Fraction> weights;
  long c,a,q,slope;unsigned loops;
};
// Whole-disk proof of one chart. bounds[k] majorizes on |z|<=R the
// coefficient of eps^k (counted from the lowest source order) of the analytic
// factor above, so Cauchy's inequality bounds its Taylor coefficients.
struct Chart {
  Local u,f;std::vector<Local> weights;B centre;bool monomial;std::vector<M> bounds;
};
inline Chart prove(const Integrand& g,const std::vector<B>& gamma,const B& centre,bool monomial,
    const M& radius,unsigned depth,slong bits) {
  B disk(0);radius.add_error_to(disk);
  Chart chart{localize(g.u,centre,disk,bits),localize(g.f,centre,disk,bits),{},centre,monomial,{}};
  const auto ur=chart.u.enclosure/chart.u.centre_value,fr=chart.f.enclosure/chart.f.centre_value;
  auto prefactor=ball_power(ur,g.c)*ball_power(fr,-g.a);
  auto logarithm=B(static_cast<long>(g.loops)+1)*ball_log(ur)-B(static_cast<long>(g.loops))*ball_log(fr);
  if(monomial) {
    const auto ratio=B(1)+disk/centre;
    prefactor=prefactor*ball_power(ratio,g.q);
    if(g.slope)logarithm=logarithm+B(g.slope)*ball_log(ratio);
  }
  const auto p=M::upper_abs(prefactor),l=M::upper_abs(logarithm);
  // source[k]=sum_{j+m=k} |gamma_j| p l^m/m!
  std::vector<M> source(depth+1);M power=M::one();
  for(unsigned m=0;m<=depth;++m) {
    for(unsigned j=0;j+m<=depth && j<gamma.size();++j)source[j+m]+=p*power*M::upper_abs(gamma[j]);
    power=power*l/M::lower_abs(B(static_cast<long>(m)+1));
  }
  chart.bounds.assign(depth+1,M());
  for(unsigned e=0;e<=depth && e<g.weights.size();++e) {
    chart.weights.push_back(localize(g.weights[e],centre,disk,bits));
    const auto w=M::upper_abs(chart.weights.back().enclosure);
    for(unsigned k=e;k<=depth;++k)chart.bounds[k]+=w*source[k-e];
  }
  for(const auto& bound:chart.bounds)if(!bound.is_finite())throw std::domain_error("nonfinite scalar functional disk majorant");
  return chart;
}
// Taylor series about the chart centre of the same analytic factor.
inline std::vector<Jet> expand(const Integrand& g,const Chart& chart,const std::vector<B>& gamma,
    unsigned depth,unsigned length,slong bits) {
  const Jet zero(0,length,bits);
  const auto constant=[&](const B& value){auto j=zero;j.set(0,value);return j;};
  auto ur=chart.u.series(length)*constant(B(1)/chart.u.centre_value),fr=chart.f.series(length)*constant(B(1)/chart.f.centre_value);
  auto prefactor=ur.pow(g.c)*fr.pow(-g.a);
  auto logarithm=constant(B(static_cast<long>(g.loops)+1))*ur.log()-constant(B(static_cast<long>(g.loops)))*fr.log();
  if(chart.monomial) {
    auto ratio=constant(B(1));ratio.set(1,B(1)/chart.centre);
    prefactor=prefactor*ratio.pow(g.q);
    if(g.slope)logarithm=logarithm+constant(B(g.slope))*ratio.log();
  }
  // source[k]=sum_{j+m=k} gamma_j prefactor logarithm^m/m!
  std::vector<Jet> source(depth+1,zero);auto power=prefactor;
  for(unsigned m=0;m<=depth;++m) {
    for(unsigned j=0;j+m<=depth && j<gamma.size();++j)source[j+m]=source[j+m]+power*constant(gamma[j]);
    if(m<depth)power=power*logarithm*constant(B(1)/B(static_cast<long>(m)+1));
  }
  std::vector<Jet> out(depth+1,zero);
  for(unsigned e=0;e<chart.weights.size();++e) {
    if(chart.weights[e].enclosure.is_zero())continue;
    const auto weight=chart.weights[e].series(length);
    for(unsigned k=e;k<=depth;++k)out[k]=out[k]+weight*source[k-e];
  }
  return out;
}
} // namespace scalar_functional_detail

// Integrate weight(x,eps)*J_power(x,eps) over [0,1] with the exact rational
// weight combined before deflation/expansion. J is known in closed form at
// every point, so the interval is covered by charts on which the integrand is
// expanded exactly and integrated term by term; nothing is transported from
// chart to chart. The endpoint chart keeps x^(q+slope*eps) symbolic and is
// integrated meromorphically, with absolute log-moment bounds on its omitted
// tail; every other chart is a symmetric interval about its centre inside a
// disk of proved analyticity, with a Cauchy bound on its omitted tail.
// Negative weight epsilon valuation shifts the source and requests gamma
// lookahead. No quadrature enters the solver. Geometry and weight must belong
// to the same ExactField, with geometry independent of eps and other
// parameters specialized. The lower weight must be jointly analytic after
// monomial deflation; the upper endpoint is regular.
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
  kernel::PrecisionScope precision(options.working_bits);

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
  const auto bits=B::precision();
  Integrand g{compile_fraction(uf.regular,variable),compile_fraction(ff.regular,variable),{},c,a,q,slope,loops};
  for(const auto& coefficient:weight_coefficients)g.weights.push_back(compile_fraction(coefficient,variable));
  Rational radius=requested_radius;std::optional<Chart> origin;
  for(unsigned attempt=0;attempt<options.max_disk_attempts && !origin;++attempt) {
    try {origin.emplace(prove(g,gamma.coefficients,B(0),false,M::upper_abs(B::from_strings(radius.str())),maximum,bits));}
    catch(const std::domain_error&) {radius=radius/Rational(2);}
  }
  if(!origin)throw std::runtime_error("certified scalar functional endpoint disk-proof budget exhausted");
  const auto analytic=expand(g,*origin,gamma.coefficients,maximum,taylor_order+1,bits);
  const auto& analytic_bounds=origin->bounds;
  const Rational quarter=radius/Rational(4);
  const Rational overlap=quarter<Rational("1/8")?quarter:Rational("1/8");
  const auto h=B::from_strings(overlap.str());
  const long first_denominator=static_cast<long>(taylor_order)+q+2;
  if(first_denominator<=0)throw std::domain_error("certified scalar functional omitted tail is not integrable");
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
        // x^(D-1) with D<0 and no epsilon in the exponent has no analytic
        // continuation: h^D/D would be a Hadamard finite part, not the integral.
        if(denominator<0 && !slope)throw std::domain_error("scalar functional endpoint has a power divergence not regulated by epsilon");
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
  // A regular terminal point is a deliberate contract of this scalar path.
  point[variable]=sample.constant(1);
  if(geometry.U.substitute(point).is_zero() || geometry.F.substitute(point).is_zero())
    throw std::domain_error("scalar functional needs a separate singular upper-endpoint operation");
  // Charts from the overlap to the upper endpoint. A chart is the interval
  // [c-s,c+s] inside the proved disk of radius R=s/rho; with |A_n|<=M/R^n
  // its odd powers cancel and the omitted even ones are bounded by
  //   2 M s rho^(N+1) / ((N+2)(1-rho)).
  // A failed disk proof halves R. rho=sixtyfourths/64 follows the largest
  // value this bound allows for the proved M, so neither the tail limit nor
  // the chart length is wasted; the estimate only proposes, the bound decides.
  if(top>=source_low) {
    const unsigned depth=static_cast<unsigned>(top-source_low);
    const auto gamma_order=static_cast<unsigned>(std::max<long>(0,static_cast<long>(epsilon_top)-weights.epsilon_power));
    const bool monomial=q!=0 || slope!=0;
    const auto limit=M::lower_abs(B::from_strings("1/1"+std::string(options.requested_digits+20,'0')));
    const double log_limit=-std::log(10.)*(options.requested_digits+20.);
    Rational position=overlap,half=(Rational(1)-overlap)/Rational(2);unsigned charts=0;long sixtyfourths=16;
    while(position<Rational(1)) {
      if(charts++>=options.max_charts)throw std::runtime_error("certified scalar functional chart budget exhausted before the upper endpoint");
      const Rational remaining=(Rational(1)-position)/Rational(2);if(remaining<half)half=remaining;
      std::optional<Chart> chart;TadpoleBoundary seed;B step,scale;std::vector<M> tails;long next_sixtyfourths=sixtyfourths;
      for(unsigned attempt=0;attempt<options.max_disk_attempts && !chart;++attempt) {
        const Rational middle=position+half;point[variable]=sample.constant(middle);
        const auto u_middle=geometry.U.substitute(point).rational(),f_middle=geometry.F.substitute(point).rational();
        if(f_middle.is_zero() || u_middle<=Rational(0))
          throw std::domain_error("certified scalar functional real path meets a zero of U or F");
        try {
          seed=tadpole(u_middle,f_middle,v,loops,d0,gamma_order,options.f_rim);
          const auto centre=B::from_strings(middle.str()),disk_radius=B::from_strings((half*Rational(64)/Rational(sixtyfourths)).str());
          step=B::from_strings(half.str());
          auto candidate=prove(g,seed.coefficients,centre,monomial,M::upper_abs(disk_radius),depth,bits);
          acb_pow_si(scale.raw(),centre.raw(),weights.x_power,bits);
          const auto ratio=M::upper_abs(step)/M::lower_abs(disk_radius);
          const auto gap=M::positive_difference_lower(M::one(),ratio);
          if(gap.is_zero())throw std::domain_error("certified scalar functional chart outside its disk");
          const auto common=M::from_ui(2)*M::upper_abs(scale)*M::upper_abs(step)*ratio.power_upper(taylor_order+1)/
            (gap*M::lower_abs(B(static_cast<long>(taylor_order)+2)));
          bool small=true;tails.clear();double largest=0;
          for(const auto& bound:candidate.bounds){tails.push_back(common*bound);if(!(tails.back()<=limit))small=false;largest=std::max(largest,tails.back().approximate_upper());}
          // tail ~ rho^(N+2) at fixed R: the admissible rho for this M.
          const double current=sixtyfourths/64.;
          const long admissible=largest>0?std::clamp(static_cast<long>(64*0.97*current*std::exp((log_limit-std::log(largest))/(taylor_order+2.))),1L,24L):24L;
          if(small){chart.emplace(std::move(candidate));next_sixtyfourths=admissible;break;}
          const long shorter=std::min(admissible,sixtyfourths-1);
          if(shorter<1)throw std::domain_error("certified scalar functional tail cannot be reduced on this disk");
          half=half*Rational(shorter)/Rational(sixtyfourths);sixtyfourths=shorter;
        }catch(const std::domain_error&) {half=half/Rational(2);}
      }
      if(!chart)throw std::runtime_error("certified scalar functional chart disk/tail-proof budget exhausted");
      const auto series=expand(g,*chart,seed.coefficients,depth,taylor_order+1,bits);
      for(unsigned k=0;k<=depth;++k) {
        B sum(0),power=step;const auto square=step*step;
        for(unsigned n=0;n<=taylor_order;n+=2){sum+=series[k].at(n)*power/B(static_cast<long>(n)+1);power=power*square;}
        sum=B(2)*scale*sum;tails[k].add_error_to(sum);
        matched[source_low+static_cast<int>(k)-low]+=sum;
      }
      // Start the next chart on a disk twice as large, at its admissible rho.
      position+=half*Rational(2);half=half*Rational(2)*Rational(next_sixtyfourths)/Rational(sixtyfourths);sixtyfourths=next_sixtyfourths;
    }
  }
  const std::string tolerance="1/1"+std::string(options.requested_digits,'0');
  if(!native_enclosure_meets_tolerance({matched},tolerance))
    throw std::runtime_error("certified scalar functional final enclosure does not meet requested absolute digits; arithmetic precision or Taylor order is insufficient");
  return {low,std::move(matched),true};
}
} // namespace scalar_functional_detail

// A singular upper endpoint is handled by reversing the second half of the
// same scalar functional. Both arms certify their full omitted tails; their
// final sum is checked again. This is the same meromorphic DR endpoint
// prescription on each arm, with an explicit common F rim. No chart crosses an
// interior zero of U or F: its disk proof cannot succeed there.
inline LaurentValue certified_scalar_functional(const Symanzik& geometry,unsigned loops,int d0,
    unsigned power,const Exact& weight,unsigned epsilon_top,const CertifiedDeepestBetaOptions& options={}) {
  using B=kernel::ComplexBall;
  if(options.working_bits<64 || options.working_bits>1000000 || !options.requested_digits || options.requested_digits>10000 || (options.f_rim!=1 && options.f_rim!=-1))
    throw std::invalid_argument("invalid scalar functional precision or F rim");
  kernel::PrecisionScope precision(options.working_bits);
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
// Beta functional of the deepest scalar boundary: the certified scalar
// functional of power left+right with weight
// normalization*x^(left-1)*(1-x)^(right-1). A regular upper endpoint is
// required; the lower F constant uses options.f_rim.
inline LaurentValue certified_deepest_beta(const Symanzik& geometry,unsigned loops,int d0,
    unsigned left,unsigned right,unsigned epsilon_top,const CertifiedDeepestBetaOptions& options={}) {
  if(!left || !right || left>10000 || right>10000)throw std::invalid_argument("invalid deepest beta request");
  const auto request=merge_request({static_cast<int>(left),static_cast<int>(right)},0,1);
  const auto& names=geometry.U.variables();auto found=std::find(names.begin(),names.end(),"x");
  if(found==names.end())throw std::invalid_argument("deepest beta needs the path variable x");
  const auto x=geometry.U.variable(static_cast<std::size_t>(found-names.begin()));
  return scalar_functional_detail::integrate_one_sided(geometry,loops,d0,left+right,
    geometry.U.constant(request.normalization)*x.pow(left-1)*(geometry.U.constant(1)-x).pow(right-1),epsilon_top,options);
}
} // namespace diffexp::feynman
