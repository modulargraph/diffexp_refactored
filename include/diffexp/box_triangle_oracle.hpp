#pragma once
#include "diffexp/kernel/scalar.hpp"
#include <array>
namespace diffexp::oracle {
namespace box_triangle_detail {
using B=kernel::ComplexBall;
inline B gamma(const B& z){B r;acb_gamma(r.raw(),z.raw(),B::precision());return r;}
inline B psi(const B& z){B r;acb_digamma(r.raw(),z.raw(),B::precision());return r;}
inline B exponential(const B& z){B r;acb_exp(r.raw(),z.raw(),B::precision());return r;}
// Right-hand residues of Gamma(A+z)Gamma(c+z)Gamma(d+z)
// Gamma(-z)Gamma(b-z)^2*x^z, A=1+2eps,b=-2eps,x=1/3.
// Single poles z=n and double poles z=b+n are combined before summation.
inline B mb(const B& eps,const B& c,unsigned terms) {
 const B one(1),x=B::from_strings("1/3"),A=one+B(2)*eps,b=-B(2)*eps,d(1);
 B logx=x;acb_log(logx.raw(),logx.raw(),B::precision());
 B single=gamma(A)*gamma(b)*gamma(b)*gamma(c)*gamma(d);
 B paired=gamma(-b)*gamma(c+b)*gamma(d+b)*exponential(b*logx);
 B q=logx-psi(-b)+psi(c+b)+psi(d+b)-psi(one),sum(0);
 for(unsigned n=0;n<terms;++n) {
  sum+=single-paired*q;B k(n),next(n+1),back=-b-k-one;
  single=single*(-x)*(A+k)*(c+k)*(d+k)/(next*(b-k-one)*(b-k-one));
  paired=paired*x*(c+b+k)*(d+b+k)/(next*back);
  q+=one/back+one/(c+b+k)+one/(d+b+k)-one/next;
 }
 return sum;
}
inline B scalar(const B& e,unsigned terms) {
 const B one(1),two(2),three(3);auto e2=e*e,e3=e2*e;
 auto gm=gamma(one-e),g2m=gamma(one-two*e),g3m=gamma(one-three*e);
 auto sunset=gamma(two*e-one)*gm*gm*gm/gamma(three-three*e);
 B logthird=B::from_strings("1/3");acb_log(logthird.raw(),logthird.raw(),B::precision());
 auto sunset_t=exponential((one-two*e)*logthird)*sunset;
 auto bubble_triangle=gamma(e)*gm*gm*gamma(two*e)*g2m*g2m/(gamma(two-two*e)*gamma(two-three*e));
 auto insertion=gm*gm/(gamma(two-two*e)*g3m)*mb(e,e,terms);
 auto diagonal=gamma(-e)*gamma(-e)*gm/(g3m*g2m*g2m)*mb(e,one,terms);
 // Exact original, unmerged native-family IBP reduction, scalar target
 // {0,1,1,1,1,1,1,0,0}; see docs/box-triangle-original-ibp.txt.
 return ((-B(162)*e3+B(243)*e2-B(117)*e+B(18))/(two*e3))*sunset
  +((-B(486)*e3+B(729)*e2-B(351)*e+B(54))/(two*e3))*sunset_t
  +((-B(54)*e2+B(45)*e-B(9))/(two*e2))*bubble_triangle
  +((-B(18)*e+B(9))/e)*insertion+B(12)*diagonal;
}
}
struct BoxTriangleReference {
 static constexpr int epsilon_low=-4,epsilon_high=0;
 std::array<kernel::ComplexBall,5> coefficients;
 std::array<kernel::ComplexBall,5> resolution_differences;
 bool coefficient_enclosures_certified=false;
 kernel::ComplexBall at(int k)const{if(k< -4||k>0)throw std::out_of_range("boxtriangle reference eps[-4,0]");return coefficients[k+4];}
};
// Independent numerical reference, NOT a rigorous error certificate.
// Native D=-q^2, s=-1,t=-1/3,d=4-2eps, no exp(2gammaeps).
// Smirnov/Veretin hep-ph/9907385v3 eq21 gives diagonal MB integral;
// https://arxiv.org/pdf/hep-ph/9907385 . Insertion uses the one-loop
// box MB formula obtained from the same Feynman-parameter identity (eq20).
// Remove eps^-4 before Cauchy coefficient extraction. Two different radii,
// quadrature sizes and residue truncations supply a convergence diagnostic;
// agreement is an estimated numerical check, never proof of omitted tails.
inline BoxTriangleReference box_triangle_reference(slong bits=256) {
 using namespace box_triangle_detail;
 if(bits<192||bits>1024)throw std::invalid_argument("boxtriangle reference bits192..1024");
 const auto previous=B::precision();struct Restore{slong p;~Restore(){B::set_precision(p);}}restore{previous};B::set_precision(bits);
 auto extract=[&](unsigned samples,unsigned terms,const char* radius) {
  std::array<B,5> out;B r=B::from_strings(radius),pi;acb_const_pi(pi.raw(),bits);
  for(unsigned j=0;j<samples;++j) {
   // Half-step avoids real epsilon points; conjugate samples are retained.
   auto angle=pi*B(2*j+1)/B(samples);B phase;acb_mul_onei(phase.raw(),angle.raw());phase=exponential(phase);
   auto e=r*phase,regular=e*e*e*e*scalar(e,terms),inverse=B(1);
   for(unsigned k=0;k<5;++k){out[k]+=regular*inverse/B(samples);inverse=inverse/e;}
  }
  return out;
 };
 auto low=extract(48,128,"1/16"),high=extract(64,160,"1/20");BoxTriangleReference out;out.coefficients=high;
 for(unsigned k=0;k<5;++k){out.resolution_differences[k]=high[k]-low[k];if(!high[k].is_finite())throw std::runtime_error("nonfinite MB/Cauchy reference");}
 return out;
}
}
