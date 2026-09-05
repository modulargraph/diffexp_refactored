#pragma once
#include "diffexp/jet.hpp"
#include "diffexp/exact.hpp"

namespace diffexp::feynman {
struct TadpoleBoundary {int epsilon_low;std::vector<Jet::Ball> coefficients;};

// Gamma(v-L*d/2)/Gamma(v) U^(v-(L+1)d/2) F^(-v+Ld/2).
// The affine dimension is d0-2 eps. Integer dimension and explicit negative-F
// rims cover the current examples; unsupported U phases fail explicitly.
inline TadpoleBoundary tadpole(const Rational& U,const Rational& F,unsigned power,
    unsigned loops,int d0,unsigned epsilon_order,int rim=1) {
  using B=Jet::Ball;
  if(!power || !loops || d0%2 || U<=Rational(0) || F.is_zero() || (rim!=1 && rim!=-1) || epsilon_order>1000)
    throw std::invalid_argument("unsupported or invalid tadpole boundary parameters");
  const long a=static_cast<long>(power)-static_cast<long>(loops)*d0/2;
  const long c=static_cast<long>(power)-(static_cast<long>(loops)+1)*d0/2;
  const int low=a<=0?-1:0;
  Jet eps(0,epsilon_order-low+1,B::precision());if(eps.length()>1)eps.set(1,B(1));
  auto b=eps.constant(loops),gamma=(eps.constant(1)+b*eps).gamma();
  auto shift=eps.constant(1);
  if(a>0)for(long k=1;k<a;++k)shift=shift*(eps.constant(k)+b*eps);
  else {
    // Factor the sole epsilon zero before any series division.
    shift=eps.constant(1)/b;
    for(long k=a;k<0;++k)shift=shift/(eps.constant(k)+b*eps);
  }
  B log_u,log_f;auto u=B::from_strings(U.str()),f=B::from_strings(F.str());
  acb_log(log_u.raw(),u.raw(),B::precision());
  if(F<Rational(0)) {
    auto abs_f=B::from_strings((-F).str());acb_log(log_f.raw(),abs_f.raw(),B::precision());
    B pi;acb_const_pi(pi.raw(),B::precision());log_f+=pi*B::from_strings("0",rim==1?"1":"-1");
  } else acb_log(log_f.raw(),f.raw(),B::precision());
  auto lambda=eps.constant(0);lambda.set(0,B(loops+1)*log_u-B(loops)*log_f);
  auto exponential=(lambda*eps).exp();
  B numerator,denominator,gamma_v;
  acb_pow_si(numerator.raw(),u.raw(),c,B::precision());acb_pow_si(denominator.raw(),f.raw(),a,B::precision());
  acb_gamma(gamma_v.raw(),B(power).raw(),B::precision());
  auto constant=eps.constant(0);constant.set(0,numerator/(denominator*gamma_v));
  auto result=constant*gamma*shift*exponential;
  TadpoleBoundary out{low,{}};for(unsigned n=0;n<eps.length();++n)out.coefficients.push_back(result.at(n));
  return out;
}
} // namespace diffexp::feynman
