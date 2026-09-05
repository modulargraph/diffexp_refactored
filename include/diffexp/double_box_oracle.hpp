#pragma once
#include "diffexp/jet.hpp"
#include <array>
namespace diffexp::oracle {
struct DoubleBoxReference {
  static constexpr int epsilon_low=-4,epsilon_high=0;
  std::array<kernel::ComplexBall,5> coefficients;
  kernel::ComplexBall at(int k)const{if(k< -4||k>0)throw std::out_of_range("double-box oracle covers eps[-4,0]");return coefficients[k+4];}
};
// Smirnov hep-ph/9905323v2 eqs1,22,23,25,26:
// https://arxiv.org/pdf/hep-ph/9905323 . Fixed native massless planar fixture
// s=-1,t=-1/3,d=4-2eps, seven unit D=-q^2 propagators. His momentum
// measure and e^-2gammaeps prefactor give J_native=-3*exp(-2gammaeps)*K.
// This includes the odd seven-propagator sign; no further exp normalization.
// All logs are real Euclidean. Coefficient balls are rigorous, not an O(eps)
// remainder bound. Nielsen series follows by expanding eq26:
// S_{a,2}(z)=sum_{k>=1} H_{k-1} z^k/k^(a+1).
inline DoubleBoxReference double_box_planar_reference(slong bits=192) {
  using B=kernel::ComplexBall;
  if(bits<64||bits>4096)throw std::invalid_argument("double-box oracle bits64..4096");
  const auto previous=B::precision();struct Restore{slong p;~Restore(){B::set_precision(p);}}restore{previous};B::set_precision(bits);
  B x=B::from_strings("1/3"),z=-x,power(1),h(0),li2(0),li3(0),li4(0),s12(0),s22(0);
  for(slong k=1;k<=bits;++k) {
    power=power*z;const auto n=B(k),n2=n*n;li2+=power/n2;li3+=power/(n2*n);li4+=power/(n2*n2);
    s12+=h*power/n2;s22+=h*power/(n2*n);h+=B(1)/n;
  }
  // H_{k-1}<=k, hence every omitted absolute summand <=(1/3)^k.
  // A common geometric remainder bounds all five independent series.
  B tail;acb_pow_ui(tail.raw(),x.raw(),bits+1,bits);tail=tail/(B(1)-x);
  for(auto* value:{&li2,&li3,&li4,&s12,&s22})arb_add_error(acb_realref(value->raw()),acb_realref(tail.raw()));
  B l=x,m=B(1)+x,pi,zeta,gamma;acb_log(l.raw(),l.raw(),bits);acb_log(m.raw(),m.raw(),bits);
  acb_const_pi(pi.raw(),bits);arb_zeta_ui(acb_realref(zeta.raw()),3,bits);arb_const_euler(acb_realref(gamma.raw()),bits);
  const auto l2=l*l,l3=l2*l,l4=l2*l2,p2=pi*pi,p4=p2*p2;
  std::array<B,5> k;
  k[0]=B(-4);k[1]=B(5)*l;k[2]=-B(2)*l2+B(5)*p2/B(2);
  k[3]=-B(2)*l3/B(3)-B(11)*p2*l/B(2)+B(65)*zeta/B(3)
    -B(2)*(B(2)*li3-B(2)*l*li2-(l2+p2)*m);
  k[4]=B(4)*l4/B(3)+B(6)*p2*l2-B(88)*zeta*l/B(3)+B(29)*p4/B(30)
    -B(4)*(s22-l*s12)+B(44)*li4-B(4)*(m+B(6)*l)*li3
    +B(2)*(l2+B(2)*l*m+B(10)*p2/B(3))*li2+(l2+p2)*m*m
    -B(2)*(B(4)*l3+B(5)*p2*l-B(6)*zeta)*m/B(3);
  DoubleBoxReference out;B exponent(1);
  for(unsigned n=0;n<5;++n){for(unsigned j=n;j<5;++j)out.coefficients[j]+=B(-3)*exponent*k[j-n];exponent=exponent*(-B(2)*gamma)/B(n+1);}
  for(const auto& value:out.coefficients)if(!value.is_finite())throw std::runtime_error("double-box reference not finite");return out;
}
}
