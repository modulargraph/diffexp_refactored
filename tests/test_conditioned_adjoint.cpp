#include "diffexp/adjoint_transport.hpp"
#include <iostream>
using namespace diffexp;using B=Jet::Ball;
void check(bool b,const char* why){if(!b)throw std::runtime_error(why);}
int main(){try {
  B::set_precision(384);ExactField f({"x","eps","I"});Exact zero(f,"0"),one(f,"1");
  // y'=y/(1+x)^20 has analytic solution exp((1-(1+x)^-19)/19).
  // Compare its SAME retained N80 polynomial, including uncertain initial data.
  std::vector<RationalLineEntry> entries{{0,0,0,Exact(f,"1/(1+x)^20")}};
  auto polynomial=polynomial_transport::compile(entries,1,80,0);auto rational=adjoint_detail::compile(entries);
  check(polynomial.polynomial_rows==1,"adversarial test did not exercise polynomial recurrence");
  B value(1);arb_add_error_2exp_si(acb_realref(value.raw()),-300);Boundary input{{value}};B h=B::from_strings("1/4");
  const auto fast=[&](const Boundary& b){return polynomial_transport::chart(polynomial,b,B(0),h,80);};
  const auto old=[&](const Boundary& b){return adjoint_detail::chart(rational,b,B(0),h,80);};
  const auto unsafe=fast(input),reference=old(input);AdjointConditioningStats stats;
  const auto guarded=adjoint_detail::conditioned_chart(input,fast,old,&stats);
  check(stats.polynomial_charts==1 && stats.rational_cross_checks==1,"large denominator cancellation did not trigger bounded fallback");
  check(adjoint_detail::enclosure_quality(unsafe).approximate_upper()>1e-30 && adjoint_detail::enclosure_quality(guarded).approximate_upper()<1e-80,"conditioning fallback did not repair arithmetic wrapping");
  // Independent exact rational recurrence supplies the retained polynomial.
  std::vector<Rational> a(80,Rational(0)),c(81,Rational(0));a[0]=c[0]=Rational(1);
  for(unsigned n=1;n<80;++n)a[n]=-a[n-1]*Rational(n+19)/Rational(n);
  for(unsigned n=0;n<80;++n){for(unsigned m=0;m<=n;++m)c[n+1]+=a[m]*c[n-m];c[n+1]=c[n+1]/Rational(n+1);}
  Rational exact(0);for(unsigned n=81;n-->0;)exact=exact/Rational(4)+c[n];
  const B expected=B::from_strings(exact.str())*value;
  check(acb_contains(guarded[0][0].raw(),expected.raw()),"guarded retained polynomial excludes independent exact recurrence enclosure");
  check(acb_overlaps(guarded[0][0].raw(),reference[0][0].raw()),"intersection lost rational enclosure");
  // The same high-degree denominator also makes a polynomial homogeneous
  // column lose its accuracy reserve. Its rational fallback must retain the
  // independent exact N80 value, including the inherited input uncertainty.
  AdjointOptions map_options;AdjointConditioningStats map_stats;
  map_options.conditioning_stats=&map_stats;
  auto map_guarded=transport_adjoint_rows({{Exact(f,"-1/(1+x)^20")}},
      LaurentRows{0,0,{{{value}}}},{{zero}},{zero,Exact(f,"1/4")},map_options);
  check(map_stats.polynomial_homogeneous_columns>0 && map_stats.rational_homogeneous_columns>0,
      "inaccurate polynomial homogeneous column did not fall back");
  check(acb_contains(map_guarded.coefficients[0][0][0].raw(),expected.raw()),
      "homogeneous map fallback excluded independent exact retained polynomial");
  // Fixed reserve catches cumulative small losses, even below the 256x trigger.
  B broad(1);arb_add_error_2exp_si(acb_realref(broad.raw()),-190);
  check(adjoint_detail::needs_rational_cross_check({{broad}},{{broad}}),"fixed accuracy reserve depends on per-chart loss");
  bool rejected=false;try{adjoint_detail::intersect_retained_enclosures({{B(0)}},{{B(1)}});}catch(const std::logic_error&){rejected=true;}
  check(rejected,"disjoint retained-polynomial enclosures accepted");
  B invalid;acb_indeterminate(invalid.raw());auto recovered=adjoint_detail::intersect_retained_enclosures({{invalid}},{{B(1)}});
  check(acb_equal_si(recovered[0][0].raw(),1),"nonfinite polynomial candidate did not use finite reference");
  rejected=false;try{adjoint_detail::intersect_retained_enclosures({{invalid}},{{invalid}});}catch(const std::domain_error&){rejected=true;}
  check(rejected,"two nonfinite candidates accepted");
  AdjointOptions options;options.conditioning_stats=&stats;stats={};
  auto initial=exact_laurent_rows({{one}},zero,0);
  auto benign=transport_adjoint_rows({{Exact(f,"1/(1+x)")}},initial,{{zero}},{zero,one},options);
  check(stats.polynomial_charts>0 && !stats.rational_cross_checks && !stats.rational_compilations,"well-conditioned path lost lazy polynomial fast path");
  check(NativeTailMagnitude::upper_abs(benign.coefficients[0][0][0]-B::from_strings("1/2")).approximate_upper()<1e-40,"benign analytic solution mismatch");
  // Exercise lazy reference compilation and the independent-row batch path.
  options.max_rows_per_batch=1;stats={};
  auto many=exact_laurent_rows({{one},{one},{one}},zero,0);
  auto batched=transport_adjoint_rows({{Exact(f,"-1/(1+x)^20")}},many,{{zero},{zero},{zero}},{zero,Exact(f,"1/4")},options);
  check(stats.rational_cross_checks>0 && stats.rational_compilations==1,"batched fallback did not reuse one lazy compilation per leg");
  check(acb_overlaps(batched.coefficients[0][0][0].raw(),batched.coefficients[2][0][0].raw()),"batched fallback changed independent row result");
  std::cout<<"Conditioned adjoint: exact adversarial retained polynomial, fixed reserve, intersections, nonfinite rejection, lazy fast path and batching passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
