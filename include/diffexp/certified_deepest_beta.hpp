#pragma once
#include "diffexp/deepest_beta.hpp"
#include "diffexp/certified_tail.hpp"
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
// Certificate-only variant: the lower primitive's omitted analytic Taylor
// terms are bounded after keeping its affine-epsilon endpoint power exact.
// A certified gamma overlap boundary then starts ordinary transport. The flag
// is set only after the final balls meet 10^(-requested_digits) absolute radius.
// A regular upper endpoint is required. The lower F constant uses f_rim;
// normalized analytic U/F ratios start at one on that sheet. Unsupported or
// inconclusive whole-disk proofs fail.
inline LaurentValue certified_deepest_beta(const Symanzik& geometry,unsigned loops,int d0,
    unsigned left,unsigned right,unsigned epsilon_top,const CertifiedDeepestBetaOptions& options={}) {
  const unsigned taylor_order=options.taylor_order;
  using B=Jet::Ball;
  if(!loops || loops>64 || d0< -64 || d0>64 || d0%2 || !left || !right || left>10000 || right>10000 || taylor_order<8 || taylor_order>1000 || epsilon_top>100)
    throw std::invalid_argument("invalid deepest beta request");
  if(options.working_bits<64 || options.working_bits>1000000 || !options.requested_digits || options.requested_digits>10000 || !options.max_charts || !options.max_disk_attempts || (options.f_rim!=1 && options.f_rim!=-1))
    throw std::invalid_argument("invalid certified deepest precision or proof budget");
  const Rational requested_radius(options.initial_radius);
  if(!(Rational(0)<requested_radius))throw std::invalid_argument("certified deepest radius must be positive");
  struct PrecisionScope {slong previous;explicit PrecisionScope(slong bits):previous(B::precision()){B::set_precision(bits);}~PrecisionScope(){B::set_precision(previous);}} precision(options.working_bits);
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
  auto gamma=tadpole(uf.value,ff.value,v,loops,d0,epsilon_top+primitive_poles,options.f_rim);
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
  // Whole-disk proof of the ordinary analytic factor, independent of its
  // stored Taylor coefficients. Principal logs apply to ratios equal to 1
  // at the origin, so the explicit endpoint x^(q+slope*epsilon) is untouched.
  using M=NativeTailMagnitude;
  auto ur_exact=uf.regular/uf.regular.constant(uf.value);
  auto fr_exact=ff.regular/ff.regular.constant(ff.value);
  const auto signed_power=[](const Exact& value,long power){return power<0?value.constant(1)/value.pow(-power):value.pow(power);};
  auto p_exact=signed_power(ur_exact,c)*signed_power(fr_exact,-a)*(geometry.U.constant(1)-geometry.U.variable(variable)).pow(right-1)*geometry.U.constant(request.normalization);
  auto p_expression=data::Reader(p_exact.str()).read();
  auto l_expression=data::Reader(std::to_string(loops+1)+"*Log["+ur_exact.str()+"]-"+std::to_string(loops)+"*Log["+fr_exact.str()+"]").read();
  Rational radius=requested_radius;M p_bound,l_bound;bool disk_proved=false;
  for(unsigned attempt=0;attempt<options.max_disk_attempts;++attempt) {
    B disk(0);M::upper_abs(B::from_strings(radius.str())).add_error_to(disk);
    try {
      auto p=native_tail_detail::eval(p_expression,disk,"x").value;
      auto l=native_tail_detail::eval(l_expression,disk,"x").value;
      native_tail_detail::require_finite(p);native_tail_detail::require_finite(l);
      p_bound=M::upper_abs(p);l_bound=M::upper_abs(l);disk_proved=true;break;
    }catch(const std::domain_error&) {radius=radius/Rational(2);}
  }
  if(!disk_proved)throw std::runtime_error("certified deepest endpoint disk-proof budget exhausted");
  const Rational quarter=radius/Rational(4);
  const Rational overlap=quarter<Rational("1/8")?quarter:Rational("1/8");
  const auto h=B::from_strings(overlap.str());
  const long first_denominator=static_cast<long>(taylor_order)+q+2;
  if(first_denominator<=0)throw std::domain_error("certified deepest omitted tail is not integrable");
  std::vector<M> analytic_bounds(maximum+1);M logarithm_factor=M::one();
  for(unsigned m=0;m<=maximum;++m) {
    for(unsigned j=0;j+m<=maximum && j<gamma.coefficients.size();++j)
      analytic_bounds[j+m]+=p_bound*logarithm_factor*M::upper_abs(gamma.coefficients[j]);
    logarithm_factor=logarithm_factor*l_bound/M::lower_abs(B::from_strings(std::to_string(m+1)));
  }
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
  // The unseen powers all have n+q>-1; integrate their log moments
  // absolutely. Unlike the retained resonant primitive they have no DR pole.
  // For D=N+q+2 and H=-log(h), normalized moments obey
  // T_0=1/D, T_j=H^j/(j! D)+T_(j-1)/D.
  const auto radius_lower=M::lower_abs(B::from_strings(radius.str()));
  const auto ratio=M::upper_abs(h)/radius_lower;
  const auto gap=M::positive_difference_lower(M::one(),ratio);
  if(gap.is_zero())throw std::domain_error("certified deepest overlap outside endpoint disk");
  const auto common=M::upper_abs(scale)*ratio.power_upper(taylor_order+1)/gap;
  const auto d_lower=M::lower_abs(B::from_strings(std::to_string(first_denominator)));
  const auto h_log=M::upper_abs(log_h);
  const auto slope_abs=M::upper_abs(B(slope));
  const unsigned tail_depth=static_cast<unsigned>(top-gamma.epsilon_low);
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
      error.add_error_to(matched[gamma.epsilon_low+static_cast<int>(k+j)-low]);
    }
    slope_power=slope_power*slope_abs;
  }
  std::vector<Exact> point;const auto sample=geometry.U;
  for(std::size_t i=0;i<sample.variable_count();++i)point.push_back(sample.variable(i));
  point[variable]=sample.constant(overlap);
  auto seed=tadpole(geometry.U.substitute(point).rational(),geometry.F.substitute(point).rational(),v,loops,d0,epsilon_top,options.f_rim);
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
  const auto compiled=compile_rational_entries(entries);
  const auto local_limit=M::lower_abs(B::from_strings("1/1"+std::string(options.requested_digits+20,'0')));
  Rational center(0);unsigned charts=0;
  while(center<Rational(1)) {
    if(charts++>=options.max_charts)throw std::runtime_error("certified deepest ordinary chart budget exhausted");
    Rational witness=requested_radius;std::optional<NativeAnalyticDisk> disk;
    const B center_ball=B::from_strings(center.str());
    // Translate the exact rational disk model before interval evaluation.
    // This removes catastrophic dependency in expanded denominators near a
    // small overlap; the Taylor recurrence still reuses its compiled matrix.
    std::vector<Exact> shifted_point;for(std::size_t i=0;i<sample.variable_count();++i)shifted_point.push_back(sample.variable(i));
    shifted_point[variable]=parameter+sample.constant(center);
    std::vector<NativeAnalyticEntry> centered_matrix;
    for(const auto& e:entries)centered_matrix.push_back({e.row,e.column,e.epsilon,data::Reader(e.coefficient.substitute(shifted_point).str()).read()});
    Rational advance(0),step_divisor(4);B step;std::optional<NativeTaylorTail> tail;
    for(unsigned attempt=0;attempt<options.max_disk_attempts;++attempt) {
      bool shrink_radius=true;
      try {
        auto candidate=NativeAnalyticDisk::certify(centered_matrix,2,top-low,B(0),witness.str());
        const Rational remaining=Rational(1)-center,quarter_step=witness/step_divisor;
        advance=remaining<quarter_step?remaining:quarter_step;
        if(!(Rational(0)<advance))throw std::runtime_error("certified deepest transport made no progress");
        step=B::from_strings(advance.str());
        auto candidate_tail=certify_native_taylor_tail(candidate,boundary,step,taylor_order);
        bool small=true;for(const auto& error:candidate_tail.absolute)if(!(error<=local_limit))small=false;
        if(small){disk.emplace(std::move(candidate));tail.emplace(std::move(candidate_tail));break;}
        // Reduce h/R independently of the witness radius; otherwise a fixed
        // geometric factor creates an artificial accuracy floor at this N.
        if(candidate.radius_upper()*candidate.ordinary_growth_norm()<=M::from_ui(8)) {
          step_divisor=step_divisor*Rational(2);shrink_radius=false;
        }
      }catch(const std::domain_error&) {}
      if(shrink_radius)witness=witness/Rational(2);
    }
    if(!disk || !tail)throw std::runtime_error("certified deepest ordinary disk/tail-proof budget exhausted");
    auto polynomial=rational_chart(compiled,boundary,center_ball,step,taylor_order);
    for(auto& row:polynomial)for(unsigned k=0;k<row.size();++k)tail->absolute[k].add_error_to(row[k]);
    boundary=std::move(polynomial);center+=advance;
  }
  const std::string tolerance="1/1"+std::string(options.requested_digits,'0');
  if(!native_enclosure_meets_tolerance({boundary[1]},tolerance)) {
    double largest=0;for(const auto& value:boundary[1])largest=std::max(largest,std::max(mag_get_d(arb_radref(acb_realref(value.raw()))),mag_get_d(arb_radref(acb_imagref(value.raw())))));
    std::ostringstream message;message<<"certified deepest final enclosure does not meet requested absolute digits; largest Cartesian radius="<<largest<<", charts="<<charts;
    throw std::runtime_error(message.str());
  }
  return {low,std::move(boundary[1]),true};
}
} // namespace diffexp::feynman
