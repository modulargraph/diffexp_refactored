#include "diffexp/adjoint_transport.hpp"
#include <iostream>
using namespace diffexp;
using B=Jet::Ball;
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
// Independent oracle. u=g^T solves u'=-A^T u+f^T*source, here integrated row
// by row with the epsilon-expanded full-convolution executor at fixed
// quarter-radius steps: no finite-lag recurrence, no step selection, no
// shared code with the adjoint transport below the exact rational input.
LaurentRows forward_reference(const ExactEpsilonMatrix& a,const LaurentRows& input,const ExactEpsilonMatrix& forcing,
    const std::vector<Exact>& path,unsigned order) {
  const auto d=input.columns();const auto [xi,ei]=path_epsilon_variables(path[0]);(void)xi;
  int low=std::min(0,input.low);
  for(const auto& row:forcing)for(const auto& c:row)if(!c.is_zero())low=std::min(low,static_cast<int>(*exact_epsilon_valuation(c,ei)));
  const auto zero=path[0].constant(0),eps=zero.variable(ei);
  LaurentRows out{low,input.high,{}};
  for(std::size_t i=0;i<input.coefficients.size();++i) {
    ExactEpsilonMatrix system(d+1,std::vector<Exact>(d+1,zero));
    for(std::size_t j=0;j<d;++j) {
      for(std::size_t k=0;k<d;++k)system[j][k]=-a[k][j];
      system[j][d]=forcing[i][j]*eps.pow(-low);
    }
    LaurentBoundary start{low,Boundary(d+1,std::vector<B>(input.high-low+1,B(0))),false};
    for(std::size_t j=0;j<d;++j)for(int k=input.low;k<=input.high;++k)start.values[j][k-low]=input.coefficients[i][j][k-input.low];
    start.values[d][0]=B(1);
    auto end=transport_laurent(system,std::move(start),path,order);
    end.values.pop_back();out.coefficients.push_back(std::move(end.values));
  }
  return out;
}
int main(){try {
  B::set_precision(256);ExactField field({"x","eps","I"});
  auto e=[&](const char* text){return Exact(field,text);};
  const auto zero=e("0"),one=e("1");
  // Signed Laurent forcing, a moving epsilon-dependent pole, independent
  // observable rows and complex legs.
  ExactEpsilonMatrix a{{e("eps/(2+x+eps)"),zero},{e("1/(3+x)"),zero}};
  ExactEpsilonMatrix forcing{{e("x/eps"),zero},{zero,e("2*x/eps")}};
  LaurentRows input{-2,2,{{{B(1),B(2),B(0),B(0),B(0)},{B(2),B(0),B(1),B(0),B(0)}},
                         {{B(3),B(0),B(0),B(0),B(0)},{B(4),B(1),B(0),B(0),B(0)}}}};
  std::vector<Exact> path{zero,e("I/32"),e("1/16+I/32"),e("1/16")};
  const auto reference=forward_reference(a,input,forcing,path,80);
  for(unsigned batch:{0,1}) {
    AdjointConditioningStats stats;AdjointOptions options;options.conditioning_stats=&stats;options.max_rows_per_batch=batch;
    auto output=transport_adjoint_rows(a,input,forcing,path,options);
    require(output.low==-2 && output.high==2,"changed legitimate negative Laurent orders");
    require(stats.charts>0 && stats.recurrences>0 && stats.circuit_preparations>0,"transport work was not counted");
    for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j)for(unsigned k=0;k<5;++k) {
      const auto delta=output.coefficients[i][j][k]-reference.coefficients[i][j][k];
      mag_t bound;mag_init(bound);acb_get_mag(bound,delta.raw());
      require(mag_cmp_2exp_si(bound,-140)<0,"adjoint transport disagrees with the independent forward convolution");mag_clear(bound);
    }
  }
  // Inherited uncertainty beyond the half-precision reserve in a nilpotent
  // flow: the exact polynomial solution is an independent oracle, no omitted
  // Taylor tail enters, and no step may be wasted on subdivision.
  a={{zero,e("-1/4")},{zero,zero}};
  B first(1),second(1);
  arb_add_error_2exp_si(acb_realref(first.raw()),-80);
  arb_add_error_2exp_si(acb_realref(second.raw()),-80);
  input={0,0,{{{first},{second}}}};forcing={{zero,zero}};path={zero,one};
  AdjointConditioningStats stats;AdjointOptions options;options.taylor_order=16;options.conditioning_stats=&stats;
  auto recovered=transport_adjoint_rows(a,input,forcing,path,options);
  require(stats.centered_charts>0 && stats.conditioning_subdivisions==0,"inherited noise was not carried by the homogeneous map");
  const Boundary expected{{first},{second+B::from_strings("0.25")*first}};
  for(unsigned j=0;j<2;++j) {
    require(acb_contains(recovered.coefficients[0][j][0].raw(),expected[j][0].raw()),
            "inherited-noise chart excluded the exact uncertain solution");
    mag_t radius;mag_init(radius);mag_hypot(radius,arb_radref(acb_realref(recovered.coefficients[0][j][0].raw())),arb_radref(acb_imagref(recovered.coefficients[0][j][0].raw())));
    require(mag_cmp_2exp_si(radius,-77)<0,"inherited-noise chart inflated the exact uncertain solution");mag_clear(radius);
  }
  // Direct FLINT conversions preserve arbitrary-size signed rational values.
  auto large=e("-123456789012345678901234567890123456789/987654321");
  require(large.constant(large.rational())==large,"direct exact rational conversion changed value");
  require(*exact_epsilon_valuation(e("eps^7/(eps^11*(1+x))"),1)==-4,
          "direct exponent valuation changed signed order");
  std::cout<<"compact adjoint integration passed\n";return 0;
}catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}}
