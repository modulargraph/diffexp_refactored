#include "diffexp/banana_oracle.hpp"
#include <iostream>
using namespace diffexp::oracle;
void require(bool yes,const char* why) { if(!yes) throw std::runtime_error(why); }
int main() { try {
  Banana4Options options;options.target_bits=64;options.working_bits=160;
  auto equal=banana4_equal_bessel(options);
  Ball reference;
  arb_set_str(acb_realref(reference.raw()),"39.655526834297652529992823046933581156446060218710",160);
  require(acb_contains(equal.value.raw(),reference.raw()),"equal-mass reference not enclosed");
  auto unequal=banana4_unequal_bessel(options);
  require(arb_lt(acb_realref(unequal.value.raw()),acb_realref(equal.value.raw())),"unequal result should be smaller");
  require(equal.analytic_queries>0 && unequal.analytic_queries>0,"analytic error-bound queries not exercised");
  std::cout<<"equal ";acb_printn(equal.value.raw(),40,0);std::cout<<"\nunequal ";acb_printn(unequal.value.raw(),40,0);
  std::cout<<"\nevaluations "<<equal.evaluations<<" "<<unequal.evaluations<<"\n";
  Banana4Options refined=options;refined.target_bits=112;refined.working_bits=224;
  auto unequal_refined=banana4_unequal_bessel(refined);
  require(acb_contains(unequal.value.raw(),unequal_refined.value.raw()),"unequal result inconsistent on refinement");
  Ball unequal_reference;
  arb_set_str(acb_realref(unequal_reference.raw()),"28.647167500233570163016706",224);
  require(acb_contains(unequal.value.raw(),unequal_reference.raw()),"unequal numerical reference not enclosed");
  // Deliberately cross exp(u)'s right-half-plane boundary in an analytic
  // query. A midpoint-only callback would incorrectly return a finite value.
  detail::BananaContext context;context.max_evaluations=100;
  for(auto& m:context.masses)acb_one(m.raw());
  Ball input,output;arb_add_error_2exp_si(acb_imagref(input.raw()),1);
  detail::banana_density(output.raw(),input.raw(),&context,1,160);
  require(!output.is_finite(),"analytic query across K0 branch domain must fail");
  bool rejected=false;
  try {banana4_bessel({Rational(1),Rational(1),Rational(1),Rational(1),Rational("1/2")},options);}
  catch(const std::invalid_argument&) {rejected=true;}
  require(rejected,"unsupported lower-tail mass domain must fail");
  // N=2 bubble has an elementary Feynman-parameter integral, independently
  // checking normalization and the rescaling used for sub-unit masses.
  auto elementary_bubble=[](long four_m2_plus_one) {
    Ball root(four_m2_plus_one),inverse,value;
    acb_sqrt(root.raw(),root.raw(),256);acb_inv(inverse.raw(),root.raw(),256);
    acb_atanh(value.raw(),inverse.raw(),256);acb_mul(value.raw(),value.raw(),inverse.raw(),256);
    acb_mul_2exp_si(value.raw(),value.raw(),2);
    // The elementary expression is real on 0<1/sqrt(4m²+1)<1.
    require(arb_contains_zero(acb_imagref(value.raw())),"elementary bubble imaginary enclosure");
    arb_zero(acb_imagref(value.raw()));return value;
  };
  auto bubble=equal_banana_bessel(1,refined);
  require(acb_contains(bubble.value.raw(),elementary_bubble(5).raw()),"bubble elementary reference not enclosed");
  auto small_mass=banana_bessel({Rational("1/4"),Rational("1/4")},options);
  require(acb_contains(small_mass.value.raw(),elementary_bubble(2).raw()),"sub-unit mass rescaling reference not enclosed");
  auto sunrise=equal_banana_bessel(2,refined);
  Ball sunrise_reference;
  arb_set_str(acb_realref(sunrise_reference.raw()),"2.2367927002126465105292636268841307993871175",224);
  require(acb_contains(sunrise.value.raw(),sunrise_reference.raw()),"full sunrise independently refined reference not enclosed");
  auto three_loop=equal_banana_bessel(3,options);
  Ball three_loop_reference;
  arb_set_str(acb_realref(three_loop_reference.raw()),"8.2681045358689687315430153454799888687286185",224);
  require(acb_contains(three_loop.value.raw(),three_loop_reference.raw()),"three-loop normalization reference not enclosed");
  for(const auto& masses:std::vector<std::vector<Rational>>{{Rational(1)},
      {Rational(1),Rational(0)},{Rational(1),Rational(-1)},std::vector<Rational>(6,Rational(1))}) {
    rejected=false;try{banana_bessel(masses,options);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"unsupported propagator count or nonpositive masses accepted");
  }
  std::cout<<"sunrise ";acb_printn(sunrise.value.raw(),48,0);std::cout<<"\n";
  options.max_evaluations=1;rejected=false;
  try {banana4_equal_bessel(options);}catch(const std::runtime_error&) {rejected=true;}
  require(rejected,"exhausted budget must fail explicitly");
  std::cout<<"unequal refined ";acb_printn(unequal_refined.value.raw(),48,0);std::cout<<"\n";
  std::cout<<"massive banana independent Bessel oracle PASS\n";
} catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 1;} }
