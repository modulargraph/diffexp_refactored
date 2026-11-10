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
  AdjointOptions options;options.taylor_order=16;options.max_rows_per_batch=1;
  for(unsigned subdivisions:{0,12}) {
    AdjointConditioningStats stats;options.conditioning_stats=&stats;options.max_conditioning_halvings=subdivisions;
    auto output=transport_adjoint_rows(a,initial,forcing,path,options);
    require(output.low==-3 && output.high==2,"physical Laurent window changed");
    // Both batches carry uncertainty beyond the half-precision reserve.
    require(stats.centered_charts==2 && stats.conditioning_subdivisions==0,"uncertain batches were not carried by the homogeneous map");
    auto expected=initial;
    for(unsigned r=0;r<2;++r)for(unsigned k=0;k<6;++k) {
      if(k>=1)expected.coefficients[r][1][k]+=B::from_strings("0.25")*initial.coefficients[r][0][k-1];
      if(k>=2)expected.coefficients[r][1][k]+=B::from_strings("0.125")*initial.coefficients[r][2][k-2];
    }
    for(unsigned r=0;r<2;++r)for(unsigned j=0;j<3;++j)for(unsigned k=0;k<6;++k) {
      require(acb_contains(output.coefficients[r][j][k].raw(),expected.coefficients[r][j][k].raw()),"closed uncertain solution excluded");
      mag_t radius;mag_init(radius);mag_hypot(radius,arb_radref(acb_realref(output.coefficients[r][j][k].raw())),arb_radref(acb_imagref(output.coefficients[r][j][k].raw())));
      require(mag_cmp_2exp_si(radius,-77)<0,"closed uncertain solution inflated");mag_clear(radius);
    }
  }
  // Exact input retains the direct recurrence, requiring no uncertainty map.
  auto exact=initial;
  for(auto& row:exact.coefficients)for(auto& series:row)for(auto& value:series)
    acb_get_mid(value.raw(),value.raw());
  AdjointConditioningStats exact_stats;options.conditioning_stats=&exact_stats;
  transport_adjoint_rows(a,exact,forcing,path,options);
  require(exact_stats.centered_charts==0 && exact_stats.recurrences==2,"exact input needlessly constructed an uncertainty map");
  // Noise confined to a column no row depends on needs no map recurrence
  // beyond its trivial unit column.
  exact.coefficients[0][1][3]=B(3);
  arb_add_error_2exp_si(acb_realref(exact.coefficients[0][1][3].raw()),-80);
  auto accumulated=transport_adjoint_rows(a,exact,forcing,path,options);
  require(acb_contains(accumulated.coefficients[0][1][3].raw(),exact.coefficients[0][1][3].raw()),"accumulator noise was not retained exactly");
  // A rational nilpotent connection: the exact flow is
  // y + 1000*log((3+x)/3)*(y1-y2)*(1,1). Interval propagation of the
  // uncertain input through the recurrence would lose the cancellation.
  B::set_precision(512);
  ExactEpsilonMatrix cancellation{{e("-1000/(3+x)"),e("-1000/(3+x)")},
                                  {e("1000/(3+x)"),e("1000/(3+x)")}};
  B first(1),second(2);
  arb_add_error_2exp_si(acb_realref(first.raw()),-80);
  arb_add_error_2exp_si(acb_realref(second.raw()),-80);
  LaurentRows noisy{0,0,{{{first},{second}}}};
  ExactEpsilonMatrix no_forcing{{z,z}};
  AdjointConditioningStats stats;AdjointOptions defaults;defaults.conditioning_stats=&stats;
  auto carried=transport_adjoint_rows(cancellation,noisy,no_forcing,{z,e("1/2")},defaults);
  require(stats.centered_charts>0 && stats.conditioning_subdivisions==0,"uncertain rational nilpotent flow was not carried by the map");
  B logarithm;acb_log(logarithm.raw(),B::from_strings("7/6").raw(),512);
  const auto correction=B(1000)*logarithm*(first-second);
  const std::vector<B> flow{first+correction,second+correction};
  for(unsigned j=0;j<2;++j) {
    const auto& value=carried.coefficients[0][j][0];
    require(acb_overlaps(value.raw(),flow[j].raw()),"carried uncertainty disagrees with the closed flow");
    B middle;acb_get_mid(middle.raw(),flow[j].raw());
    require(NativeTailMagnitude::upper_abs(value-middle).approximate_upper()<1e-20,"carried uncertainty was inflated");
  }
  std::cout<<"uncertainty support tests passed\n";return 0;
}catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}}
