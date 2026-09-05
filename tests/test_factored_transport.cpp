#include "diffexp/factored_transport.hpp"
#include <iostream>
using namespace diffexp;
namespace ft=factored_transport;
namespace lb=linear_boundary;
using B=Jet::Ball;
void check(bool yes,const char* why){if(!yes)throw std::runtime_error(why);}
void near(const B& a,const B& b){auto delta=a-b;mag_t m;mag_init(m);acb_get_mag(m,delta.raw());bool good=mag_cmp_2exp_si(m,-90)<0;mag_clear(m);check(good,"factored analytic map mismatch");}
double radius(const B& b){return mag_get_d(arb_radref(acb_realref(b.raw())));}
LaurentRows rows(unsigned r,unsigned c,int low,int high){return {low,high,std::vector(r,std::vector(c,std::vector<B>(high-low+1,B(0))))};}
int main(){try{
  B::set_precision(384);ExactField field({"x","eps"});auto e=[&](const char* text){return Exact(field,text);};
  auto source=std::make_shared<const LaurentBoundary>(LaurentBoundary{0,{{B(2),B(0),B(0)}},false});
  auto initial=rows(1,1,0,2);for(auto& v:initial.coefficients[0][0])v=B(1);
  ft::Options options;options.transport.taylor_order=64;
  auto result=ft::evolve({{e("0")}},lb::Expression{initial,source},{{e("1/eps")}},
      {e("0"),e("1/4")},1,options);
  check(result.required_initial_high==2 && result.integrated.transform.low==-1,"factored negative forcing demand/lower bound");
  check(result.physical.leaf_source==source && result.integrated.leaf_source==source,"factored transport replaced shared source");
  for(const auto& value:result.integrated.transform.coefficients[0][0])near(value,B::from_strings("1/4"));
  initial.high=1;initial.coefficients[0][0].resize(2);bool demanded=false;
  try{(void)ft::evolve({{e("0")}},lb::Expression{initial,source},{{e("1/eps")}}, {e("0"),e("1/4")},1,options);}
  catch(const ft::MapDemand& demand){demanded=demand.required_high==2;}
  check(demanded,"factored missing upper map coefficients invented as zero");
  // Dense rank-one nilpotent connection: A^2=0, so the forward map and
  // integrated observables have exact linear/quadratic analytic formulas.
  std::vector<int> u{1,2,3},v{1,1,-1};
  ExactEpsilonMatrix a(3,std::vector<Exact>(3,e("0"))),b{{e("1"),e("2"),e("1")},{e("-1"),e("1"),e("3")}};
  auto m=rows(3,2,0,0);for(unsigned i=0;i<3;++i)for(unsigned j=0;j<2;++j){a[i][j]=e("0");m.coefficients[i][j][0]=B(2*i+j+1);}
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)a[i][j]=e("0").constant(u[i]*v[j]);
  auto leaf2=std::make_shared<const LaurentBoundary>(LaurentBoundary{0,{{B(1)},{B(1)}},false});
  auto dense=ft::evolve(a,{m,leaf2},b,{e("0"),e("1/4")},0,options);
  int weights[2][3]={{1,2,1},{-1,1,3}};
  for(unsigned j=0;j<2;++j){B contraction(0);for(unsigned k=0;k<3;++k)contraction+=B(v[k])*m.coefficients[k][j][0];
    for(unsigned i=0;i<3;++i)near(dense.physical.transform.coefficients[i][j][0],m.coefficients[i][j][0]+B(u[i])*contraction/B(4));
    for(unsigned i=0;i<2;++i){B expected(0);for(unsigned k=0;k<3;++k)expected+=B(weights[i][k])*(m.coefficients[k][j][0]/B(4)+B(u[k])*contraction/B(32));near(dense.integrated.transform.coefficients[i][j][0],expected);}}
  auto terminal=rows(2,3,0,0);auto minus_b=b;for(auto& row:minus_b)for(auto& value:row)value=-value;
  auto adjoint=transport_adjoint_rows(a,terminal,minus_b,{e("1/4"),e("0")},options.transport);
  auto direct_comparison=lb::compose(adjoint,m,0);
  for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j)near(dense.integrated.transform.coefficients[i][j][0],direct_comparison.coefficients[i][j][0]);
  // Logarithmic observables share a pole-bearing uncertain source. Combine
  // their maps first; materializing them separately discards the correlation.
  B uncertain(2);arb_add_error_2exp_si(acb_realref(uncertain.raw()),-10);
  auto shared=std::make_shared<const LaurentBoundary>(LaurentBoundary{-1,{{uncertain,B(0)}},false});
  auto repeated=rows(2,1,0,1);repeated.coefficients[0][0][0]=B(1);repeated.coefficients[1][0][0]=B(1);
  auto logs=ft::evolve({{e("1/(1-x)"),e("0")},{e("0"),e("1/(1-x)")}},
      {repeated,shared},{{e("1"),e("0")},{e("0"),e("-1")}}, {e("0"),e("1/4")},1,options);
  B logarithm;acb_log(logarithm.raw(),B::from_strings("3/4").raw(),384);
  near(logs.integrated.transform.coefficients[0][0][0],-logarithm);
  near(logs.integrated.transform.coefficients[1][0][0],logarithm);
  auto sum=rows(1,2,0,1);sum.coefficients[0][0][0]=B(1);sum.coefficients[0][1][0]=B(1);
  auto combined=lb::compose(sum,logs.integrated,1);
  auto final=lb::materialize(combined,0);
  auto separately=lb::materialize(logs.integrated,0);
  auto independent=separately.values[0][0]+separately.values[1][0];
  check(radius(final.values[0][0])<radius(independent)*1e-20,"factored map failed to preserve source correlation");
  check(!logs.omitted_tail_certified && combined.leaf_source==shared,"factored source/certificate contract");
  std::cout<<"Factored forward/integrated maps, epsilon lookahead, dense analytic/adjoint equivalence and logarithmic shared-source cancellation passed\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
