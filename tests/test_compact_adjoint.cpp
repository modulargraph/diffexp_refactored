#include "diffexp/adjoint_checkpoint.hpp"
#include <iostream>
using namespace diffexp;
using B=Jet::Ball;
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
int main(){try {
  B::set_precision(256);ExactField field({"x","eps","I"});
  auto e=[&](const char* text){return Exact(field,text);};
  const auto zero=e("0"),one=e("1");
  // Signed Laurent forcing, independent observable rows and complex legs.
  ExactEpsilonMatrix a{{e("eps/(2+x+eps)"),zero},{e("1/(3+x)"),zero}};
  ExactEpsilonMatrix forcing{{e("x/eps"),zero},{zero,e("2*x/eps")}};
  LaurentRows input{-2,2,{{{B(1),B(2),B(0),B(0),B(0)},{B(2),B(0),B(1),B(0),B(0)}},
                         {{B(3),B(0),B(0),B(0),B(0)},{B(4),B(1),B(0),B(0),B(0)}}}};
  std::vector<Exact> path{zero,e("I/32"),e("1/16+I/32"),e("1/16")};
  AdjointOptions legacy;legacy.taylor_order=32;legacy.rational_circuit_recurrence=false;
  auto reference=transport_adjoint_rows(a,input,forcing,path,legacy);
  AdjointOptions options=legacy;options.rational_circuit_recurrence=true;
  require(adjoint_checkpoint::identity(a,input,forcing,path,options)!=
          adjoint_checkpoint::identity(a,input,forcing,path,legacy),"different kernels share checkpoint identity");
  for(unsigned batch:{0,1}) {
    AdjointConditioningStats stats;options.conditioning_stats=&stats;options.max_rows_per_batch=batch;
    auto output=transport_adjoint_rows(a,input,forcing,path,options);
    require(output.low==-2 && output.high==2,"changed legitimate negative Laurent orders");
    require(stats.circuit_charts>0 && stats.circuit_preparations>0,"compact kernel was not used");
    for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j)for(unsigned k=0;k<5;++k) {
      const auto delta=output.coefficients[i][j][k]-reference.coefficients[i][j][k];
      mag_t bound;mag_init(bound);acb_get_mag(bound,delta.raw());
      require(mag_cmp_2exp_si(bound,-170)<0,"same-order legacy/compact mismatch");mag_clear(bound);
    }
  }
  // Nilpotent flow has a closed polynomial solution but independent input
  // radii need the full map action; this tests prepared sharing under fallback.
  a={{e("-1000"),e("-1000")},{e("1000"),e("1000")}};
  B first(1),second(2);arb_add_error_2exp_si(acb_realref(first.raw()),-220);
  arb_add_error_2exp_si(acb_realref(second.raw()),-220);
  input={0,0,{{{first},{second}}}};forcing={{zero,zero}};path={zero,one};
  options={};options.centered_before_rational=true;options.centered_map_workers=2;
  AdjointConditioningStats stats;options.conditioning_stats=&stats;
  auto result=transport_adjoint_rows(a,input,forcing,path,options);
  require(stats.circuit_homogeneous_columns==2 && stats.circuit_preparations==2,
          "homogeneous columns did not share chart preparation");
  require(stats.rational_compilations==0,"compact centered route expanded fallback expression trees");
  Boundary expected{{B(1001)*first-B(1000)*second},{B(1000)*first-B(999)*second}};
  for(unsigned j=0;j<2;++j)require(acb_contains(result.coefficients[0][j][0].raw(),expected[j][0].raw()),
                                  "nilpotent exact uncertain flow excluded");
  require(B::precision()==256,"map workers leaked precision");
  // Inherited uncertainty exceeds the half-precision reserve, but a nilpotent
  // chart grows it by <2x. The old predicate skipped centering and needlessly
  // invoked a full legacy recurrence at every step. The exact polynomial flow
  // is an independent oracle; no omitted Taylor tail enters this case.
  a={{zero,e("-1/4")},{zero,zero}};
  first=B(1);second=B(1);
  arb_add_error_2exp_si(acb_realref(first.raw()),-80);
  arb_add_error_2exp_si(acb_realref(second.raw()),-80);
  input={0,0,{{{first},{second}}}};
  options={};options.taylor_order=16;options.centered_before_rational=true;
  options.compact_centered_recovery=false;
  AdjointConditioningStats old_stats;options.conditioning_stats=&old_stats;
  auto old_result=transport_adjoint_rows(a,input,forcing,path,options);
  const auto old_identity=adjoint_checkpoint::identity(a,input,forcing,path,options);
  options.compact_centered_recovery=true;
  AdjointConditioningStats new_stats;options.conditioning_stats=&new_stats;
  auto recovered=transport_adjoint_rows(a,input,forcing,path,options);
  require(old_stats.rational_cross_checks>0,"regression failed to exercise inherited-reserve reference");
  require(new_stats.rational_cross_checks==0 && new_stats.rational_compilations==0 &&
          new_stats.centered_charts>0,"compact recovery repeated legacy uncertainty comparison");
  require(new_stats.conditioning_subdivisions==0,"inherited noise triggered pointless subdivisions");
  require(old_identity!=adjoint_checkpoint::identity(a,input,forcing,path,options),
          "changed recovery policy reused old checkpoint identity");
  expected={{first},{second+B::from_strings("0.25")*first}};
  for(unsigned j=0;j<2;++j) {
    require(acb_contains(recovered.coefficients[0][j][0].raw(),expected[j][0].raw()),
            "compact inherited-noise recovery excluded exact uncertain solution");
    require(acb_overlaps(recovered.coefficients[0][j][0].raw(),old_result.coefficients[0][j][0].raw()),
            "changed recovery policy disagrees with legacy retained enclosure");
  }
  // Direct FLINT conversions preserve arbitrary-size signed rational values.
  auto large=e("-123456789012345678901234567890123456789/987654321");
  require(large.constant(large.rational())==large,"direct exact rational conversion changed value");
  require(*exact_epsilon_valuation(e("eps^7/(eps^11*(1+x))"),1)==-4,
          "direct exponent valuation changed signed order");
  std::cout<<"compact adjoint integration passed\n";return 0;
}catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}}
