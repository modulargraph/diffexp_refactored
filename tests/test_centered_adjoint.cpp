#include "diffexp/adjoint_transport.hpp"
#include <iostream>
using namespace diffexp;
using B=Jet::Ball;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
int main(){try {
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
  AdjointConditioningStats stats;AdjointOptions options;options.conditioning_stats=&stats;
  auto centered=transport_adjoint_rows(connection,start,forcing,{zero,one},options);
  require(stats.centered_charts==1 && stats.homogeneous_chart_maps==1,"integrated nilpotent fallback was not used");
  require(stats.polynomial_homogeneous_columns==2 && stats.rational_homogeneous_columns==0,"exact nilpotent map did not use finite-lag columns");
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
    options=AdjointOptions{};options.polynomial_recurrence=polynomial;options.max_rows_per_batch=batch;
    stats={};options.conditioning_stats=&stats;
    auto result=transport_adjoint_rows(connection,data,laurent_forcing,path,options);
    require(result.low==-1 && result.high==1 && stats.centered_charts>0,"Laurent/complex test did not exercise centered action");
    require(polynomial?stats.polynomial_homogeneous_columns>0:stats.rational_homogeneous_columns>0,"homogeneous map ignored recurrence selection");
    for(unsigned observable=0;observable<2;++observable)for(unsigned j=0;j<2;++j)for(unsigned k=0;k<3;++k) {
      const auto& value=result.coefficients[observable][j][k];
      require(acb_overlaps(value.raw(),analytic[j][k].raw()),"centered epsilon convolution/forcing disagrees with exact cubic solution");
      const auto error=NativeTailMagnitude::upper_abs(value-analytic[j][k]).approximate_upper();
      if(error>=1e-40)std::cerr<<"mode="<<polynomial<<" batch="<<batch<<" row="<<j<<" coefficient="<<k<<" error="<<error<<'\n';
      require(error<1e-40,"centered epsilon/complex path failed forty-digit analytic comparison");
    }
    require(stats.conditioning_subdivisions>0,"complex forcing did not exercise chart rejection and rollback");
  }
  require(B::precision()==384,"centered fallback changed working precision");
  std::cout<<"Centered adjoint passed exact nilpotent and cubic Laurent/complex solutions, both recurrences, row batching and memory cap\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
