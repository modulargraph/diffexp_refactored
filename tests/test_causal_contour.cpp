#include "diffexp/scalar_functional.hpp"
#include "diffexp/families.hpp"
#include <iostream>
using namespace diffexp;
using B=kernel::ComplexBall;
static void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
int main(){try {
  B::set_precision(256);ExactField field({"x","eps","I"});Exact x(field,"x");
  // Gamma(1+eps)*(-1 +/- i0)^(-1-eps): finite term -1,
  // epsilon coefficient EulerGamma +/- i*pi. Opposite sheets first differ here.
  auto upper=feynman::tadpole(Rational(1),Rational(-1),2,1,2,1,1);
  auto lower=feynman::tadpole(Rational(1),Rational(-1),2,1,2,1,-1);
  require(arb_is_positive(acb_imagref(upper.coefficients[1].raw())),"upper F rim must yield +i*pi");
  require(arb_is_negative(acb_imagref(lower.coefficients[1].raw())),"lower F rim must yield -i*pi");
  auto scalar=feynman::certified_scalar_functional({x.constant(1),x.constant(-1)},1,2,2,x.constant(1),1);
  require(scalar.taylor_tail_certified,"constant negative-F scalar should enclose its selected sheet");
  require(acb_overlaps(scalar.at(1).raw(),upper.coefficients[1].raw()),"current scalar default uses upper F rim");
  require(!acb_overlaps(scalar.at(1).raw(),lower.coefficients[1].raw()),"tail certificate must not be mistaken for lower-F prescription");
  auto rejected=[&](const Exact& f,const char* label){
    feynman::CertifiedDeepestBetaOptions options;options.max_charts=12;options.max_disk_attempts=12;options.taylor_order=24;options.requested_digits=8;
    bool failed=false;try{feynman::certified_scalar_functional({x.constant(1),f},1,2,2,x.constant(1),0,options);}catch(const std::exception& e){failed=true;std::cout<<label<<": "<<e.what()<<'\n';}
    require(failed,"unsupported real scalar path unexpectedly returned a certificate");
  };
  auto endpoint=feynman::certified_scalar_functional({x.constant(1),x.constant(1)-x},1,2,2,x.constant(1),0);
  B euler;arb_const_euler(acb_realref(euler.raw()),256);
  require(endpoint.taylor_tail_certified && acb_contains(endpoint.at(-1).raw(),B(-1).raw()) && acb_contains(endpoint.at(0).raw(),euler.raw()),"upper endpoint DR primitive -Gamma(1+eps)/eps");
  rejected(x.constant(1)-x.constant(2)*x,"interior F zero");
  // Same positive imaginary x gives opposite F rims for opposite simple-root slopes.
  auto z=B::from_strings("1/2","1/10");auto increasing=B(2)*z-B(1),decreasing=B(1)-B(2)*z;
  require(arb_is_positive(acb_imagref(increasing.raw())) && arb_is_negative(acb_imagref(decreasing.raw())),"x orientation alone is not a universal F rim");
  // Exact physical Henn simplex slice: alpha_2=x, all other alpha=1.
  // Homogeneous rescaling puts every x>0 strictly inside the projective simplex.
  auto henn=feynman::example_family("henn_double_pentagon_x0");henn.momenta.lines.resize(8);
  std::vector<Exact> parameters(8,x.constant(1));parameters[1]=x;
  auto geometry=feynman::symanzik(henn.momenta,parameters);
  require(geometry.U==Exact(field,"5*x+16") && geometry.F==Exact(field,"3*x-29"),"Henn timelike interior threshold slice");
  std::cout<<"causal contour audit probes PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
