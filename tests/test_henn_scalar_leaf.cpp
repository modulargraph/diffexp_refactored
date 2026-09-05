#include "diffexp/recursion_leaf.hpp"
#include "diffexp/causal.hpp"
#include <iostream>
using namespace diffexp;using B=Jet::Ball;
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
int main(){try {
  B::set_precision(384);ExactField field({"x","eps","I"});Exact x(field,"x"),eps(field,"eps"),dimension(field,"4-2*eps");
  const auto family=feynman::example_family("henn_double_pentagon_x0");
  auto raw=ibp::quadratic_family(family.momenta,x,family.physical_count);
  const auto anchors=causal::henn_anchors();auto prescription=causal::current_example(family.name,7,anchors,true);
  for(unsigned level=0;level<6;++level)raw=recursion::at_anchor(ibp::PropagatorBasis(ibp::merge(raw,0,1,x)),0,anchors[level]);
  require(ibp::PropagatorBasis(raw).physical_count==2,"Henn frozen ladder must reach exactly two physical denominators");
  const ibp::PropagatorBasis final(ibp::merge(raw,0,1,x));
  ibp::Integral target(final.denominators.size(),0);target[0]=2;
  auto unit=gaussian::reduce(final,target,dimension);
  require(unit.U==Exact(field,"x*(43750-27391*x)/62500") && unit.F==Exact(field,"x^2*(87680943*x-243432000)/7812500000"),"prescribed-anchor native Henn scalar geometry changed");
  // These exact linear factors prove all nonzero roots lie outside [0,1],
  // and fix U>0,F<0 on (0,1); root plots alone are not the sign proof.
  const Rational uroot("43750/27391"),froot("243432000/87680943");
  require(uroot>Rational(1) && froot>Rational(1),"Henn leaf has an open-unit U/F zero");
  feynman::CertifiedDeepestBetaOptions options;options.working_bits=384;options.requested_digits=20;options.f_rim=prescription.f_rim;
  auto integrate=[&](const gaussian::Boundary& source,const Exact& weight){
    auto value=feynman::certified_scalar_functional({source.U,source.F},source.loops,4,source.power,source.multiplier*weight,4,options);
    require(value.taylor_tail_certified && native_enclosure_meets_tolerance({value.coefficients},"1/100000000000000000000"),"Henn scalar probe must enclose both arithmetic and Taylor tails at20digits");return value;
  };
  auto scalar=integrate(unit,x.constant(1));require(scalar.epsilon_low==-1,"unit Henn scalar gamma pole");
  // The fully accumulated unit-power beta request would have left power7,
  // right power1. This tests that scalar functional, not the complete ladder.
  target[0]=8;auto accumulated=gaussian::reduce(final,target,dimension);
  auto beta=integrate(accumulated,x.constant(7)*x.pow(6));
  require(beta.epsilon_low==0,"Henn accumulated unit-power scalar should start at epsilon zero");
  // A genuine auxiliary numerator produces an exact rational Gaussian weight.
  target[0]=2;target[1]=-1;auto numerator=gaussian::reduce(final,target,dimension);
  require(!numerator.multiplier.is_zero() && numerator.U==unit.U && numerator.F==unit.F,"Henn auxiliary numerator geometry/weight");
  auto numerator_value=integrate(numerator,x.constant(1));
  auto rational=integrate(numerator,x.constant(1)/(eps.pow(2)*(x.constant(2)+x+eps)));
  require(rational.epsilon_low==numerator_value.epsilon_low-2,"Henn rational-weight epsilon pole lookahead");
  options.f_rim=1;auto opposite=feynman::certified_scalar_functional({unit.U,unit.F},unit.loops,4,unit.power,x.constant(1),4,options);
  for(int k=scalar.epsilon_low;k<=4;++k){B conjugate;acb_conj(conjugate.raw(),scalar.at(k).raw());require(acb_overlaps(conjugate.raw(),opposite.at(k).raw()),"Henn leaf F rims must be conjugate");}
  std::cout<<"Henn scalar leaf: exact prescribed-anchor geometry, no open-unit roots, lower F rim, inherited beta power, auxiliary numerator and rational epsilon lookahead certified; no FIRE/full-ladder claim\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
