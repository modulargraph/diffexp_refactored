#pragma once
#include "diffexp/box_oracle.hpp"
#include "diffexp/families.hpp"
#include <flint/acb_hypgeom.h>

namespace diffexp::oracle {
namespace pentagon_detail {
using B=kernel::ComplexBall;
struct Precision {slong previous;explicit Precision(slong bits):previous(B::precision()){if(bits<64||bits>4096)throw std::invalid_argument("oracle requires 64..4096 bits");B::set_precision(bits);}~Precision(){B::set_precision(previous);}};
inline B rational(const Rational& r){return B::from_strings(r.str());}
inline B logarithm(const Rational& r){B out=rational(r);acb_log(out.raw(),out.raw(),B::precision());return out;}
}
// Ellis/Zanderighi arXiv:0712.1851v4 eqs2.1,2.2,4.19:
// https://arxiv.org/pdf/0712.1851 . The paper removes rGamma; restore it here.
// Four massless propagators, external squares (0,0,0,m2), diagonals s,t.
// Native D=-q^2, measure d^d l/(i*pi^(d/2)), d=4-2eps, mu^2=1,
// no exp(gamma*eps). Only s,t,m2<0; real branches, either F rim agrees.
// These are rigorous coefficient balls, not a finite-epsilon remainder bound.
inline BoxReference one_mass_box_reference(const Rational& s,const Rational& t,const Rational& m2,slong bits=192) {
  using namespace pentagon_detail;Precision precision(bits);
  if(s.sign()>=0 || t.sign()>=0 || m2.sign()>=0)throw std::invalid_argument("one-mass box requires s,t,m2<0");
  const auto ls=logarithm(-s),lt=logarithm(-t),lm=logarithm(-m2);
  B gamma,pi,ds,dt;arb_const_euler(acb_realref(gamma.raw()),bits);acb_const_pi(pi.raw(),bits);
  auto zs=rational(Rational(1)-m2/s),zt=rational(Rational(1)-m2/t);
  acb_hypgeom_dilog(ds.raw(),zs.raw(),bits);acb_hypgeom_dilog(dt.raw(),zt.raw(),bits);
  const auto denominator=rational(s*t),sum=ls+lt-lm;
  BoxReference out;
  out.coefficients[0]=rational(Rational(2)/(s*t));
  out.coefficients[1]=-B(2)*(sum+gamma)/denominator;
  out.coefficients[2]=(B(2)*ls*lt-lm*lm-B(2)*(ds+dt)-pi*pi/B(2)+B(2)*gamma*sum+gamma*gamma)/denominator;
  for(const auto& c:out.coefficients)if(!c.is_finite())throw std::runtime_error("one-mass box coefficient evaluation failed");
  return out;
}
struct PentagonReference {
  static constexpr int epsilon_low=-2,epsilon_high=0;
  std::array<kernel::ComplexBall,3> coefficients;
  std::array<Rational,5> box_weights;
  std::array<std::array<Rational,3>,5> box_kinematics; // s,t,m2 for each omitted denominator
  kernel::ComplexBall at(int k)const {if(k< -2||k>0)throw std::out_of_range("pentagon oracle only covers eps -2..0");return coefficients[k+2];}
};
// Bern/Dixon/Kosower hep-ph/9306240 eqs2.2–2.4 and3.23–3.24:
// https://arxiv.org/pdf/hep-ph/9306240 . Set their auxiliary alpha_i=1.
// F=a^T S a, Sij=-(q_i-q_j)^2/2, b=(1/2)S^-1*1.
// Native J5=sum b_i*J4(pinchi)+O(eps). The D=6-2eps pentagon is finite,
// so the dimension-shift term proportional to eps contributes no eps[-2,0].
// Their In[1] is already the positive Gamma(n-d/2)*F^(d/2-n)
// parameter integral: DO NOT add another (-1)^5 conversion.
// Contract: five ordered cyclic massless external legs, massless internal
// lines, strictly negative nonadjacent squared distances. Exceptional/singular
// kinematics or non-unit loop coefficient are explicitly unsupported.
inline PentagonReference massless_pentagon_reference(const feynman::Family& family,slong bits=192) {
  using namespace pentagon_detail;Precision precision(bits);family.validate();
  if(family.loops!=1 || family.lines.size()!=5)throw std::invalid_argument("pentagon oracle needs one loop and five lines");
  for(const auto& line:family.lines)if(!line.mass_squared.is_zero() || line.loop_coefficients.size()!=1 || line.loop_coefficients[0]!=Rational(1))
    throw std::invalid_argument("pentagon oracle needs massless unit loop-coefficient lines");
  std::array<std::array<Rational,5>,5> distances;
  for(unsigned i=0;i<5;++i)for(unsigned j=0;j<5;++j) {
    Rational square(0);for(unsigned a=0;a<family.external_gram.size();++a)for(unsigned b=0;b<family.external_gram.size();++b)
      square=square+(family.lines[i].external_coefficients[a]-family.lines[j].external_coefficients[a])*family.external_gram[a][b]*(family.lines[i].external_coefficients[b]-family.lines[j].external_coefficients[b]);
    distances[i][j]=square;
    const bool adjacent=i==j || (i+1)%5==j || (j+1)%5==i;
    if(adjacent?!square.is_zero():square.sign()>=0)throw std::invalid_argument("pentagon oracle needs cyclic on-shell legs and Euclidean nonadjacent invariants");
  }
  // Exact elimination for S*b=1/2. No fitted numerical reduction coefficients.
  std::array<std::array<Rational,6>,5> matrix;
  for(unsigned i=0;i<5;++i){for(unsigned j=0;j<5;++j)matrix[i][j]=-distances[i][j]/Rational(2);matrix[i][5]=Rational("1/2");}
  for(unsigned column=0;column<5;++column) {
    unsigned pivot=column;while(pivot<5&&matrix[pivot][column].is_zero())++pivot;
    if(pivot==5)throw std::domain_error("singular pentagon kinematic S matrix");
    std::swap(matrix[column],matrix[pivot]);const auto divisor=matrix[column][column];for(unsigned j=column;j<6;++j)matrix[column][j]=matrix[column][j]/divisor;
    for(unsigned row=0;row<5;++row)if(row!=column){const auto multiplier=matrix[row][column];for(unsigned j=column;j<6;++j)matrix[row][j]=matrix[row][j]-multiplier*matrix[column][j];}
  }
  PentagonReference out;
  for(unsigned omitted=0;omitted<5;++omitted) {
    out.box_weights[omitted]=matrix[omitted][5];
    std::array<unsigned,4> remaining;unsigned k=0;for(unsigned i=0;i<5;++i)if(i!=omitted)remaining[k++]=i;
    auto s=distances[remaining[0]][remaining[2]],t=distances[remaining[1]][remaining[3]];Rational m2(0);
    unsigned massive=0;for(unsigned i=0;i<4;++i){auto edge=distances[remaining[i]][remaining[(i+1)%4]];if(!edge.is_zero()){m2=edge;++massive;}}
    if(massive!=1)throw std::logic_error("pentagon pinch did not produce a one-mass box");
    out.box_kinematics[omitted]={s,t,m2};auto box=one_mass_box_reference(s,t,m2,bits);
    for(unsigned e=0;e<3;++e)out.coefficients[e]+=rational(out.box_weights[omitted])*box.coefficients[e];
  }
  return out;
}
inline PentagonReference massless_pentagon_reference(slong bits=192){return massless_pentagon_reference(feynman::example_family("pentagon").momenta,bits);}
} // namespace diffexp::oracle
