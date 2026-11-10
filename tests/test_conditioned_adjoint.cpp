#include "diffexp/adjoint_transport.hpp"
#include <iostream>
using namespace diffexp;using B=Jet::Ball;
void check(bool b,const char* why){if(!b)throw std::runtime_error(why);}
double quality(const B& value){return adjoint_detail::enclosure_quality(value).approximate_upper();}
int main(){try {
  B::set_precision(384);ExactField f({"x","eps","I"});Exact zero(f,"0"),one(f,"1");
  // y'=y/(1+x)^20 has analytic solution exp((1-(1+x)^-19)/19). Its N80
  // Taylor polynomial at a quarter of the pole distance is an arithmetic
  // fixture, not an approximation (its tail is O(0.1)): an independent exact
  // rational recurrence supplies the SAME retained polynomial.
  std::vector<RationalLineEntry> entries{{0,0,0,Exact(f,"1/(1+x)^20")}};
  B value(1);arb_add_error_2exp_si(acb_realref(value.raw()),-300);Boundary input{{value}};B h=B::from_strings("1/4");
  std::vector<Rational> a(80,Rational(0)),c(81,Rational(0));a[0]=c[0]=Rational(1);
  for(unsigned n=1;n<80;++n)a[n]=-a[n-1]*Rational(n+19)/Rational(n);
  for(unsigned n=0;n<80;++n){for(unsigned m=0;m<=n;++m)c[n+1]+=a[m]*c[n-m];c[n+1]=c[n+1]/Rational(n+1);}
  Rational exact(0);for(unsigned n=81;n-->0;)exact=exact/Rational(4)+c[n];
  const B expected=B::from_strings(exact.str())*value;
  // The quotient by a pole of multiplicity 20 is where an interval recurrence
  // wraps; the majorized finite-lag quotient must keep the input's accuracy.
  double tail=0;
  const auto compiled=rational_circuit::compile(entries,1,0,1,0,80);
  check(compiled.data().denominators.at(0).majorized && compiled.data().denominators[0].poles.at(0).multiplicity==20,"pole multiplicity of the adversarial denominator");
  rational_circuit::Statistics work;
  const auto circuit=rational_circuit::chart(compiled,input,B(0),h,80,&work);tail=work.truncation_log2;
  check(acb_contains(circuit[0][0].raw(),expected.raw()),"finite-lag retained polynomial excludes the independent exact recurrence");
  check(quality(circuit[0][0])<1e-80,"multiple-pole quotient lost its arithmetic accuracy");
  check(acb_overlaps(circuit[0][0].raw(),rational_chart(entries,input,B(0),h,80)[0][0].raw()),"finite-lag and convolution enclosures are disjoint");
  check(!taylor_tail_consistent(tail,80,384),"unconverged fixture reported a consistent Taylor tail");
  // The transport must not accept that chart: it shortens steps until the
  // retained tail meets the target and reproduces the analytic solution,
  // carrying the input uncertainty through.
  {
    AdjointConditioningStats stats;AdjointOptions options;options.conditioning_stats=&stats;
    auto accurate=transport_adjoint_rows({{Exact(f,"-1/(1+x)^20")}},LaurentRows{0,0,{{{value}}}},{{zero}},{zero,Exact(f,"1/4")},options);
    B power,analytic;acb_pow_si(power.raw(),B::from_strings("5/4").raw(),-19,384);
    const B exponent=(B(1)-power)/B(19);acb_exp(analytic.raw(),exponent.raw(),384);
    check(stats.tail_limited_charts>0 && stats.charts>1,"slowly convergent chart near a high-order pole was accepted");
    check(NativeTailMagnitude::upper_abs(accurate.coefficients[0][0][0]-analytic*value).approximate_upper()<1e-40,
        "transport near a pole of order 20 missed the analytic solution");
    check(quality(accurate.coefficients[0][0][0])<1e-80,"transport near a pole of order 20 lost the input accuracy");
    // Independent rows in separate batches give the same solution.
    options.max_rows_per_batch=1;stats={};
    auto many=exact_laurent_rows({{one},{one},{one}},zero,0);
    auto batched=transport_adjoint_rows({{Exact(f,"-1/(1+x)^20")}},many,{{zero},{zero},{zero}},{zero,Exact(f,"1/4")},options);
    for(unsigned row=0;row<3;++row)
      check(NativeTailMagnitude::upper_abs(batched.coefficients[row][0][0]-analytic).approximate_upper()<1e-40,"batched rows near a pole of order 20");
  }
  // A well-conditioned connection: one direct recurrence per chart, no map.
  {
    AdjointConditioningStats stats;AdjointOptions options;options.conditioning_stats=&stats;
    auto benign=transport_adjoint_rows({{Exact(f,"1/(1+x)")}},exact_laurent_rows({{one}},zero,0),{{zero}},{zero,one},options);
    check(stats.charts>0 && stats.recurrences==stats.charts && !stats.centered_charts && !stats.conditioning_subdivisions,"well-conditioned path did extra work");
    check(NativeTailMagnitude::upper_abs(benign.coefficients[0][0][0]-B::from_strings("1/2")).approximate_upper()<1e-45,"benign analytic solution mismatch");
  }
  // The fixed reserve is half the working precision, whatever the step.
  B broad(1);arb_add_error_2exp_si(acb_realref(broad.raw()),-190);
  check(adjoint_detail::beyond_reserve({{broad}}) && !adjoint_detail::beyond_reserve({{value}}),"fixed accuracy reserve");
  B invalid;acb_indeterminate(invalid.raw());
  check(adjoint_detail::beyond_reserve({{invalid}}),"nonfinite coefficient within reserve");
  // A known input radius can flow into an exactly zero accumulator. Its
  // exact solution is y1=y0(0), y2=x*y0(0), including the inherited radius.
  B inherited(1);arb_add_error_2exp_si(acb_realref(inherited.raw()),-130);
  AdjointConditioningStats transfer_stats;AdjointOptions transfer_options;
  transfer_options.conditioning_stats=&transfer_stats;
  auto transferred=transport_adjoint_rows({{zero,-one},{zero,zero}},
      LaurentRows{0,0,{{{inherited},{B(0)}}}},{{zero,zero}},{zero,one},transfer_options);
  for(unsigned j=0;j<2;++j)check(acb_contains(transferred.coefficients[0][j][0].raw(),inherited.raw()),
      "coupled uncertainty transfer excluded the exact interval solution");
  check(transfer_stats.conditioning_subdivisions==0,"existing uncertainty caused futile conditioning subdivision");
  B widened(0);arb_add_error_2exp_si(acb_realref(widened.raw()),-110);
  check(adjoint_detail::needs_conditioning_subdivision({{inherited},{B(0)}},{{inherited},{widened}}),
      "uncertainty transfer disabled the bounded arithmetic growth check");
  check(!adjoint_detail::needs_conditioning_subdivision({{inherited},{B(0)}},{{inherited},{inherited}}),
      "plain uncertainty transfer was treated as lost accuracy");
  std::cout<<"Conditioned adjoint: exact adversarial retained polynomial, multiple-pole quotient, tail-selected steps, fixed reserve and uncertainty transfer passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
