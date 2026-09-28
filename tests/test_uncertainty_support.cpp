#include "diffexp/adjoint_transport.hpp"
#include <iostream>
using namespace diffexp;using B=Jet::Ball;
static void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
int main(){try {
  B::set_precision(256);ExactField field({"x","eps","I"});
  auto e=[&](const char* s){return Exact(field,s);};const auto z=e("0");
  // A^2=0, hence the exact flow is polynomial and independent input boxes
  // can be enclosed without any appeal to an omitted Taylor-tail estimate.
  ExactEpsilonMatrix a{{z,e("-eps/4"),z},{z,z,z},{z,e("-eps^2/8"),z}};
  ExactEpsilonMatrix forcing(2,std::vector<Exact>(3,z));
  LaurentRows initial{-3,2,std::vector(2,std::vector(3,std::vector<B>(6,B(0))))};
  initial.coefficients[0][0][4]=B(1);initial.coefficients[1][2][2]=B(2);
  arb_add_error_2exp_si(acb_realref(initial.coefficients[0][0][4].raw()),-80);
  arb_add_error_2exp_si(acb_imagref(initial.coefficients[1][2][2].raw()),-80);
  const std::vector<Exact> path{z,e("1")};
  AdjointOptions options;options.taylor_order=16;options.centered_map_working_bits=128;
  options.centered_map_workers=2;options.max_rows_per_batch=1;
  AdjointConditioningStats baseline_stats;options.conditioning_stats=&baseline_stats;
  options.compact_centered_only=false;
  auto baseline=transport_adjoint_rows(a,initial,forcing,path,options);
  options.compact_centered_only=true;
  for(unsigned subdivisions:{0,12}) {
    AdjointConditioningStats stats;options.conditioning_stats=&stats;options.max_conditioning_halvings=subdivisions;
    auto output=transport_adjoint_rows(a,initial,forcing,path,options);
    require(output.low==-3 && output.high==2,"physical Laurent window changed");
    require(stats.centered_only_charts==2 && stats.centered_only_fallbacks==0,"centered-only route not exercised for both batches");
    require(stats.circuit_charts==2 && stats.polynomial_charts==0,"redundant whole-input source recurrence remains");
    require(stats.circuit_homogeneous_columns==2,"uncertainty support union lost later batch");
    require(stats.rational_cross_checks==0 && stats.conditioning_subdivisions==0,"closed flow needed legacy/subdivision");
    auto expected=initial;
    for(unsigned r=0;r<2;++r)for(unsigned k=0;k<6;++k) {
      if(k>=1)expected.coefficients[r][1][k]+=B::from_strings("0.25")*initial.coefficients[r][0][k-1];
      if(k>=2)expected.coefficients[r][1][k]+=B::from_strings("0.125")*initial.coefficients[r][2][k-2];
    }
    for(unsigned r=0;r<2;++r)for(unsigned j=0;j<3;++j)for(unsigned k=0;k<6;++k) {
      require(acb_contains(output.coefficients[r][j][k].raw(),expected.coefficients[r][j][k].raw()),"closed uncertain solution excluded");
      require(acb_overlaps(output.coefficients[r][j][k].raw(),baseline.coefficients[r][j][k].raw()),"centered result disagrees with reference");
    }
  }
  // Exact input retains the ordinary fast path, requiring no uncertainty map.
  for(auto& row:initial.coefficients)for(auto& series:row)for(auto& value:series)
    acb_get_mid(value.raw(),value.raw());
  AdjointConditioningStats exact_stats;options.conditioning_stats=&exact_stats;
  transport_adjoint_rows(a,initial,forcing,path,options);
  require(exact_stats.centered_only_charts==0 && exact_stats.homogeneous_chart_maps==0,
          "exact input needlessly constructed an uncertainty map");
  // Noise confined to a constant accumulator needs no active map columns.
  initial.coefficients[0][1][3]=B(3);
  arb_add_error_2exp_si(acb_realref(initial.coefficients[0][1][3].raw()),-80);
  AdjointConditioningStats accumulator_stats;options.conditioning_stats=&accumulator_stats;
  transport_adjoint_rows(a,initial,forcing,path,options);
  require(accumulator_stats.centered_only_charts==1 && accumulator_stats.circuit_homogeneous_columns==0 &&
          accumulator_stats.exact_input_map_columns_skipped==2,
          "exact active input columns were needlessly solved");
  // A rational nilpotent connection has a known finite Taylor map, but
  // reduced-precision interval cancellation consumes its arithmetic reserve.
  // Retry the compact map at full precision before any legacy column solves.
  B::set_precision(512);
  ExactEpsilonMatrix cancellation{{e("-1000/(3+x)"),e("-1000/(3+x)")},
                                  {e("1000/(3+x)"),e("1000/(3+x)")}};
  B first(1),second(2);
  arb_add_error_2exp_si(acb_realref(first.raw()),-80);
  arb_add_error_2exp_si(acb_realref(second.raw()),-80);
  LaurentRows noisy{0,0,{{{first},{second}}}};
  ExactEpsilonMatrix no_forcing{{z,z}};
  AdjointOptions retry_options;retry_options.taylor_order=16;
  retry_options.centered_map_working_bits=128;retry_options.centered_map_workers=2;
  AdjointConditioningStats retry_stats;retry_options.conditioning_stats=&retry_stats;
  const std::vector<Exact> retry_path{z,e("1/2")};
  auto retried=transport_adjoint_rows(cancellation,noisy,no_forcing,retry_path,retry_options);
  require(retry_stats.full_precision_map_retries==1 && retry_stats.rational_homogeneous_columns==0 &&
          retry_stats.rational_compilations==0 && retry_stats.reference_seconds==0,
          "deficient low-precision map entered legacy before compact full-precision retry");
  require(retry_stats.circuit_homogeneous_columns==4 && retry_stats.polynomial_homogeneous_columns==4,
          "discarded low-precision map work missing from counters");
  require(retry_stats.conditioning_subdivisions==0 && retry_stats.centered_only_fallbacks==0,
          "map retry weakened acceptance or changed chart subdivision");
  retry_options.centered_map_working_bits=512;
  AdjointConditioningStats full_stats;retry_options.conditioning_stats=&full_stats;
  auto full=transport_adjoint_rows(cancellation,noisy,no_forcing,retry_path,retry_options);
  require(full_stats.full_precision_map_retries==0 && full_stats.rational_homogeneous_columns==0,
          "full-precision control unexpectedly needed fallback");
  B logarithm(0),power(1),ratio=B(1)/B(6);
  for(unsigned n=1;n<=retry_options.taylor_order;++n) {
    power*=ratio;
    logarithm+=(n%2?B(1):B(-1))*power/B(n);
  }
  const auto correction=B(1000)*logarithm*(first-second);
  const std::vector<B> exact_retained{first+correction,second+correction};
  for(unsigned j=0;j<2;++j) {
    require(acb_overlaps(retried.coefficients[0][j][0].raw(),exact_retained[j].raw()),
            "retried map disagrees with closed retained polynomial");
    require(acb_equal(retried.coefficients[0][j][0].raw(),full.coefficients[0][j][0].raw()),
            "retry changed the full-precision retained enclosure");
  }
  // Deliberately misuse a shortened map: future callers must not silently
  // interpret uncomputed high coefficients as structural zero coefficients.
  adjoint_detail::SparseChartMap map{1,3,std::vector<adjoint_detail::SparseChartMap::Column>(1)};
  map.columns[0].minimum_input_order=2;map.columns[0].entries.push_back({0,{B(1)}});
  Boundary input{{B(1),B(0),B(0)},{B(1),B(0),B(0)}};
  arb_add_error_2exp_si(acb_realref(input[0][0].raw()),-80);
  bool rejected=false;try{adjoint_detail::centered_action(input,[](const Boundary& b){return b;},map,1);}catch(const std::logic_error&){rejected=true;}
  require(rejected,"shortened map accepted unknown coefficients");
  std::cout<<"uncertainty support tests passed\n";return 0;
}catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}}
