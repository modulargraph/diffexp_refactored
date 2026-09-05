#pragma once
#include "diffexp/jet.hpp"
namespace diffexp::oracle {
struct BoxBubbleReference {
  static constexpr int epsilon_low=-3;
  int epsilon_high=0;
  std::vector<kernel::ComplexBall> coefficients;
  kernel::ComplexBall at(int k)const{if(k<epsilon_low||k>epsilon_high)throw std::out_of_range("box-bubble coefficient outside oracle window");return coefficients.at(k-epsilon_low);}
};
// Independent loop-by-loop gamma reference. Native D=-q^2, d=4-2eps,
// measure product d^d l/(i*pi^(d/2)), no exp(2*gamma*eps), s<0.
// Native routing leaves bubble denominators (l1+p12)^2,(l1-l2)^2.
// Its integration yields Gamma(eps)Gamma(1-eps)^2/Gamma(2-2eps)
// times [-(l2+p12)^2]^-eps. The remaining triangle has shifts
// 0,p12,p123 and powers (1,1+eps,1), only p12^2=s nonzero.
// The Feynman-parametric convention is BDK hep-ph/9306240 eq2.2:
// https://arxiv.org/pdf/hep-ph/9306240 . Dirichlet integration of
// F=(-s)*a*b gives Gamma(1+2eps)Gamma(-2eps)Gamma(-eps)/
// [Gamma(1+eps)Gamma(1-3eps)] * (-s)^(-1-2eps).
// Cancelling shifted gamma factors gives the analytic germ
// eps^3*J = (-s)^(-1-2eps) Gamma(1-eps)^3 Gamma(1+2eps)/
// [2(1-2eps)Gamma(1-3eps)]. No extra (-1)^5 factor belongs here.
// Rigorous Laurent coefficient balls only; no finite-eps remainder claim.
inline BoxBubbleReference box_bubble_reference(const kernel::Rational& s=kernel::Rational(-1),unsigned top=0,slong bits=192) {
  using B=kernel::ComplexBall;
  if(s.sign()>=0||top>20||bits<64||bits>4096)throw std::invalid_argument("box-bubble oracle requires s<0, top<=20, bits64..4096");
  const auto previous=B::precision();struct Restore{slong bits;~Restore(){B::set_precision(bits);}}restore{previous};B::set_precision(bits);
  Jet eps(0,top+4,bits);eps.set(1,B(1));auto one=eps.constant(1),two=eps.constant(2);
  auto scale=eps.constant(0);scale.set(0,B::from_strings((-s).str()));
  auto regular=((one-eps).gamma().pow(3)*(one+two*eps).gamma())/(two*(one-two*eps)*(one-eps.constant(3)*eps).gamma());
  regular=regular*((-one-two*eps)*scale.log()).exp();
  BoxBubbleReference out;out.epsilon_high=top;
  for(unsigned k=0;k<top+4;++k){auto c=regular.at(k);if(!c.is_finite())throw std::runtime_error("box-bubble gamma coefficient failed");out.coefficients.push_back(c);}return out;
}
}
