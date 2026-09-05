#pragma once
#include "diffexp/kernel/scalar.hpp"
#include <array>

namespace diffexp::oracle {
struct BoxReference {
  static constexpr int epsilon_low=-2,epsilon_high=0;
  std::array<kernel::ComplexBall,3> coefficients;
  kernel::ComplexBall at(int epsilon_order) const {
    if(epsilon_order< -2 || epsilon_order>0)throw std::out_of_range("box reference covers epsilon -2..0 only");
    return coefficients[epsilon_order+2];
  }
};

// Independent analytic reference: Ellis and Zanderighi, arXiv:0712.1851v4,
// eqs. (2.1), (2.2), (4.18), https://arxiv.org/pdf/0712.1851 .
// Their I4 divides out rGamma; the native scalar restores
// rGamma=Gamma(1+eps)Gamma(1-eps)^2/Gamma(1-2eps), with mu^2=1.
// Four D_i=-q_i^2 signs multiply to +1. No exp(gamma*eps) is included.
// For S=-s>0,T=-t>0:
// I_native = rGamma/(s*t) * {2(S^-eps+T^-eps)/eps^2
//                          -log(S/T)^2-pi^2} + O(eps).
// The returned Arb balls rigorously enclose these three Laurent coefficients;
// this is NOT a bound on the finite-epsilon O(eps) remainder or an FT solver.
// Only the Euclidean massless on-shell domain s,t<0 is supported; both rims
// agree there. Kinematics and scalar normalization match the native box fixture.
inline BoxReference massless_box_reference(const kernel::Rational& s=kernel::Rational(-1),
    const kernel::Rational& t=kernel::Rational("-1/3"),slong bits=192) {
  if(s.sign()>=0 || t.sign()>=0 || bits<64 || bits>4096)
    throw std::invalid_argument("box oracle requires Euclidean s,t<0 and 64..4096 bits");
  using B=kernel::ComplexBall;B sb,tb,Ls,Lt,gamma,pi,denominator,sum,work;
  auto rational=[&](B& b,const kernel::Rational& q) {
    fmpq_t exact;fmpq_init(exact);fmpq_set_str(exact,q.str().c_str(),10);
    arb_set_fmpq(acb_realref(b.raw()),exact,bits);fmpq_clear(exact);
  };
  rational(sb,-s);rational(tb,-t);
  arb_log(acb_realref(Ls.raw()),acb_realref(sb.raw()),bits);
  arb_log(acb_realref(Lt.raw()),acb_realref(tb.raw()),bits);
  arb_mul(acb_realref(denominator.raw()),acb_realref(sb.raw()),acb_realref(tb.raw()),bits);
  arb_const_euler(acb_realref(gamma.raw()),bits);arb_const_pi(acb_realref(pi.raw()),bits);
  arb_add(acb_realref(sum.raw()),acb_realref(Ls.raw()),acb_realref(Lt.raw()),bits);
  BoxReference result;
  auto out=[&](unsigned i){return acb_realref(result.coefficients[i].raw());};
  rational(result.coefficients[0],kernel::Rational(4)/(s*t));
  arb_mul_ui(out(1),acb_realref(gamma.raw()),2,bits);
  arb_add(out(1),out(1),acb_realref(sum.raw()),bits);arb_mul_si(out(1),out(1),-2,bits);
  arb_div(out(1),out(1),acb_realref(denominator.raw()),bits);
  arb_mul(out(2),acb_realref(Ls.raw()),acb_realref(Lt.raw()),bits);
  arb_mul(acb_realref(work.raw()),acb_realref(gamma.raw()),acb_realref(sum.raw()),bits);arb_add(out(2),out(2),acb_realref(work.raw()),bits);
  arb_mul(acb_realref(work.raw()),acb_realref(gamma.raw()),acb_realref(gamma.raw()),bits);arb_add(out(2),out(2),acb_realref(work.raw()),bits);
  arb_mul_ui(out(2),out(2),2,bits);
  arb_mul(acb_realref(work.raw()),acb_realref(pi.raw()),acb_realref(pi.raw()),bits);
  arb_mul_ui(acb_realref(work.raw()),acb_realref(work.raw()),4,bits);arb_div_ui(acb_realref(work.raw()),acb_realref(work.raw()),3,bits);
  arb_sub(out(2),out(2),acb_realref(work.raw()),bits);arb_div(out(2),out(2),acb_realref(denominator.raw()),bits);
  for(const auto& c:result.coefficients)if(!c.is_finite())throw std::runtime_error("box coefficient enclosure failed");
  return result;
}
} // namespace diffexp::oracle
