#include "diffexp/adjoint_transport.hpp"
#include <iostream>
using namespace diffexp;
using B=Jet::Ball;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
int main(){try {
  require(AdjointOptions{}.centered_map_workers==1,"centered map threading changed the default");
  B::set_precision(384);const auto caller=std::this_thread::get_id();
  const int original_flint_threads=flint_get_num_threads();
  if(FLINT_USES_TLS && mpfr_buildopt_tls_p())flint_set_num_threads(3);
  std::vector<unsigned> visited(8,0);std::atomic<bool> worker_ok{true};
  const auto used=adjoint_detail::homogeneous_columns(8,4,128,[&](unsigned column) {
    if(B::precision()!=128 || flint_get_num_threads()>1)worker_ok=false;
    const auto thread=std::this_thread::get_id();
    const auto nested=adjoint_detail::homogeneous_columns(2,4,192,[&](unsigned) {
      if(std::this_thread::get_id()!=thread || B::precision()!=192)worker_ok=false;
    });
    if(nested!=1 || B::precision()!=128)worker_ok=false;
    ++visited[column];
  });
  require(worker_ok && B::precision()==384,"homogeneous workers leaked precision or nested workers");
  require(used==(FLINT_USES_TLS && mpfr_buildopt_tls_p()?4:1),"homogeneous worker limit ignored");
  for(auto count:visited)require(count==1,"homogeneous column was skipped or repeated");
  require(flint_get_num_threads()==(FLINT_USES_TLS && mpfr_buildopt_tls_p()?3:original_flint_threads),
      "homogeneous workers mutated the caller FLINT pool");
  flint_set_num_threads(original_flint_threads);
  bool propagated=false;std::atomic<unsigned> completed{0};
  try {adjoint_detail::homogeneous_columns(2,2,128,[&](unsigned column) {
    ++completed;if(column==0)throw std::logic_error("first column");throw std::runtime_error("second column");
  });}catch(const std::logic_error& error){propagated=std::string(error.what())=="first column";}
  require(propagated && completed==(FLINT_USES_TLS && mpfr_buildopt_tls_p()?2u:1u) && B::precision()==384,
      "homogeneous worker failures did not join and propagate deterministically");

  for(bool centered_first:{false,true})for(const slong map_bits:{0,128}) {
  B::set_precision(384);ExactField field({"x","eps","I"});
  // A^2=0 exactly. Taylor propagation is exactly (I+h*A)y for N>=1,
  // but propagating independent input intervals through every Taylor
  // coefficient loses this cancellation and creates fictitious higher powers.
  Exact a(field,"1000"),minus(field,"-1000");
  std::vector<RationalLineEntry> entries{{0,0,0,a},{0,1,0,minus},{1,0,0,a},{1,1,0,minus}};
  auto compiled=adjoint_detail::compile(entries);
  B first(1),second(2);arb_add_error_2exp_si(acb_realref(first.raw()),-300);
  arb_add_error_2exp_si(acb_realref(second.raw()),-300);
  Boundary input{{first},{second}};
  auto naive=adjoint_detail::chart(compiled,input,B(0),B(1),80);
  Boundary expected{{B(1001)*first-B(1000)*second},{B(1000)*first-B(999)*second}};
  Boundary mapped(2,std::vector<B>(1,B(0)));
  for(unsigned column=0;column<2;++column) {
    Boundary basis(2,std::vector<B>(1,B(0)));basis[column][0]=B(1);
    auto map=adjoint_detail::chart(compiled,basis,B(0),B(1),80);
    for(unsigned row=0;row<2;++row)mapped[row][0]+=map[row][0]*input[column][0];
  }
  require(adjoint_detail::enclosure_quality(naive).approximate_upper()>1e20,"nilpotent wrapping reproduction did not fail");
  require(adjoint_detail::enclosure_quality(mapped).approximate_upper()<1e-80,"linear chart action did not preserve cancellation");
  for(unsigned row=0;row<2;++row)require(acb_contains(mapped[row][0].raw(),expected[row][0].raw()),"independent exact nilpotent flow excluded");
  std::cout<<"Same N80/384 bits: naive quality="<<adjoint_detail::enclosure_quality(naive).approximate_upper()
    <<", linear-action quality="<<adjoint_detail::enclosure_quality(mapped).approximate_upper()<<'\n';
  Exact zero(field,"0"),one(field,"1"),eps(field,"eps"),x(field,"x");
  ExactEpsilonMatrix connection{{minus,minus},{a,a}},forcing{{zero,zero}};
  LaurentRows start{0,0,{{{first},{second}}}};
  AdjointConditioningStats stats;AdjointOptions options;options.centered_map_working_bits=map_bits;options.centered_before_rational=centered_first;options.conditioning_stats=&stats;
  auto centered=transport_adjoint_rows(connection,start,forcing,{zero,one},options);
  const auto serial_stats=stats;AdjointConditioningStats parallel_stats;
  auto parallel_options=options;parallel_options.centered_map_workers=4;parallel_options.conditioning_stats=&parallel_stats;
  parallel_options.chart_observer=[&](unsigned,double,double,const LaurentRows&) {
    require(std::this_thread::get_id()==caller,"parallel map invoked chart observer on a worker");
  };
  auto parallel=transport_adjoint_rows(connection,start,forcing,{zero,one},parallel_options);
  for(unsigned j=0;j<2;++j)require(acb_equal(centered.coefficients[0][j][0].raw(),parallel.coefficients[0][j][0].raw()),
      "parallel centered map changed a numeric ball");
  require(parallel_stats.polynomial_homogeneous_columns==serial_stats.polynomial_homogeneous_columns &&
      parallel_stats.rational_homogeneous_columns==serial_stats.rational_homogeneous_columns &&
      parallel_stats.homogeneous_chart_maps==serial_stats.homogeneous_chart_maps,
      "parallel columns changed deterministic conditioning counters");
  require(parallel_stats.max_homogeneous_map_workers==(FLINT_USES_TLS && mpfr_buildopt_tls_p()?2u:1u),
      "centered map did not use the available independent columns");
  // The compact kernel counts coefficient/pivot/result storage as well as Y.
  // This fits one complete workspace, but still cannot fit two map workers.
  parallel_options.max_taylor_cells=3*(parallel_options.taylor_order+1)+32;parallel_stats={};
  auto limited=transport_adjoint_rows(connection,start,forcing,{zero,one},parallel_options);
  require(parallel_stats.max_homogeneous_map_workers==1,"parallel map exceeded total Taylor workspace budget");
  for(unsigned j=0;j<2;++j)require(acb_equal(centered.coefficients[0][j][0].raw(),limited.coefficients[0][j][0].raw()),
      "workspace-limited parallel map changed a numeric ball");
  require(stats.centered_charts==1 && stats.homogeneous_chart_maps==1,"integrated nilpotent fallback was not used");
  if(centered_first)require(stats.centered_before_rational_charts==1 && stats.rational_cross_checks==0,
    "centered-first route repeated the unnecessary uncertain-input rational recurrence");
  require(stats.polynomial_homogeneous_columns==2 && stats.rational_homogeneous_columns==0,"exact nilpotent map did not use finite-lag columns");
  require(stats.polynomial_midpoint_charts==1 && stats.rational_midpoint_charts==0,"exact nilpotent midpoint did not use guarded finite-lag recurrence");
  for(unsigned row=0;row<2;++row)
    require(acb_contains(centered.coefficients[0][row][0].raw(),expected[row][0].raw()),"integrated chart excluded exact interval flow");
  require(adjoint_detail::enclosure_quality(Boundary{centered.coefficients[0][0],centered.coefficients[0][1]}).approximate_upper()<1e-80,"integrated chart retained fictitious interval growth");
  options.max_centered_map_cells=1;options.max_conditioning_halvings=0;stats={};
  auto bounded=transport_adjoint_rows(connection,start,forcing,{zero,one},options);
  require(stats.centered_budget_skips==1 && !stats.homogeneous_chart_maps,"centered map ignored explicit memory cap");
  require(adjoint_detail::enclosure_quality(Boundary{bounded.coefficients[0][0]}).approximate_upper()>1e20,"budget test did not retain original broad enclosure");
  // Laurent forcing and epsilon coupling: M=(1+eps)N, N^2=0,
  // b(x)=x*(1,2)/eps. The complete solution is the cubic polynomial
  // (I+x*M)y0 + x^2*b/2 + x^3*M*b/6, including complex intermediate legs.
  for(auto& row:connection)for(auto& entry:row)entry=entry*(one+eps);
  auto laurent_forcing=ExactEpsilonMatrix{{x/eps,one.constant(2)*x/eps},{x/eps,one.constant(2)*x/eps}};
  LaurentRows data{0,1,{{{first,B(3)},{second,B(4)}},{{first,B(3)},{second,B(4)}}}};
  const std::vector<Exact> path{zero,Exact(field,"I/8"),Exact(field,"1/4+I/8"),Exact(field,"1/4")};
  const B h=B::from_strings("1/4"),six=B::from_strings("1/6"),half=B::from_strings("1/2");
  Boundary analytic(2,std::vector<B>(3,B(0)));
  const B n0=B(1000)*(first-second),n1=B(-1000),nb=B(-1000);
  for(unsigned j=0;j<2;++j) {
    analytic[j][0]=h*h*half*B(j+1)+h*h*h*six*nb;
    analytic[j][1]=input[j][0]+h*n0+h*h*h*six*nb;
    analytic[j][2]=B(j+3)+h*n1+h*n0;
  }
  for(bool polynomial:{false,true})for(unsigned batch:{0,1}) {
    options=AdjointOptions{};options.centered_map_workers=4;options.centered_map_working_bits=map_bits;options.centered_before_rational=centered_first;options.polynomial_recurrence=polynomial;options.max_rows_per_batch=batch;
    stats={};options.conditioning_stats=&stats;
    auto result=transport_adjoint_rows(connection,data,laurent_forcing,path,options);
    require(result.low==-1 && result.high==1 && stats.centered_charts>0,"Laurent/complex test did not exercise centered action");
    require(polynomial?stats.polynomial_homogeneous_columns>0:stats.rational_homogeneous_columns>0,"homogeneous map ignored recurrence selection");
    require(polynomial?stats.polynomial_midpoint_charts>0:
      stats.rational_midpoint_charts>0 && stats.polynomial_midpoint_charts==0,"centered midpoint ignored recurrence selection");
    for(unsigned observable=0;observable<2;++observable)for(unsigned j=0;j<2;++j)for(unsigned k=0;k<3;++k) {
      const auto& value=result.coefficients[observable][j][k];
      require(acb_overlaps(value.raw(),analytic[j][k].raw()),"centered epsilon convolution/forcing disagrees with exact cubic solution");
      const auto error=NativeTailMagnitude::upper_abs(value-analytic[j][k]).approximate_upper();
      if(error>=1e-40)std::cerr<<"mode="<<polynomial<<" batch="<<batch<<" row="<<j<<" coefficient="<<k<<" error="<<error<<'\n';
      require(error<1e-40,"centered epsilon/complex path failed forty-digit analytic comparison");
    }
    if(map_bits)require(stats.full_precision_map_retries>0,"inadequate reduced map did not retry full precision");
    require(stats.conditioning_subdivisions>0,"complex forcing did not exercise chart rejection and rollback");
  }
  require(B::precision()==384,"centered fallback changed working precision");
  std::cout<<"Centered adjoint passed exact nilpotent and cubic Laurent/complex solutions, both recurrences, row batching and memory cap\n";
  }
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
