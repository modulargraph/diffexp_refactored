#include "diffexp/recursion_pipeline.hpp"
#include <iostream>
using namespace diffexp;
namespace ed=epsilon_demand;
namespace ft=factored_transport;
using B=Jet::Ball;
void check(bool b,const char* why){if(!b)throw std::runtime_error(why);}
template<class E,class F>void rejects(F f,const char* why){try{f();}catch(const E&){return;}throw std::runtime_error(why);}
LaurentRows rows(unsigned r,unsigned c,int low,int high){return {low,high,std::vector(r,std::vector(c,std::vector<B>(high-low+1,B(0))))};}
B at(const LaurentRows& a,unsigned i,unsigned j,int k){return k<a.low?B(0):a.coefficients.at(i).at(j).at(k-a.low);}
void agrees(const B& a,const B& b){check(a.is_finite()&&b.is_finite()&&acb_overlaps(a.raw(),b.raw()),"component/analytic retained coefficient mismatch");}
level::Result native_sunrise(const ibp::PropagatorBasis& basis,const Exact& dimension,
    const ExactField&,std::size_t parameter,const std::vector<ibp::Integral>& requested,const level::Options&) {
  check(basis.physical_count==2&&basis.denominators.size()==5,"test closure supports merged sunrise only");
  ibp::Generator generator(basis,dimension);ibp::ExactReducer reducer(dimension,300);
  ibp::for_each_seed(2,5,{1,1,50},[&](const ibp::Integral& seed){for(auto& row:generator.relations(seed))reducer.insert(std::move(row));});
  ibp::BasisReduction coordinates(reducer,{{1,0,0,0,0},{1,1,0,0,0},{1,1,-1,0,0}},dimension);
  level::Result result;result.success=true;result.ordered_basis=coordinates.ordered_basis();result.requested_integrals=requested;
  for(const auto& master:result.ordered_basis)result.matrix.push_back(coordinates.resolve(generator.derivative(master,parameter)));
  for(const auto& target:requested)result.target_rows.push_back(coordinates.resolve({{target,dimension.constant(1)}}));
  return result;
}
int main(){try{
 B::set_precision(384);ExactField field({"x","eps"});auto e=[&](const char*s){return Exact(field,s);};
 ExactEpsilonMatrix a{{e("0"),e("eps^2")},{e("0"),e("0")}};
 auto closed=ed::close(a,{2,-1},1);check(closed==ed::Highs({2,0}),"epsilon closure orientation");
 check(ed::close(a,{std::nullopt,-1},1)==ed::Highs({std::nullopt,-1}),"backward reachability direction");
 rejects<std::domain_error>([&]{ed::close({{e("0"),e("eps^-1")},{e("1"),e("0")}}, {0,std::nullopt},1);},"reachable negative epsilon cycle accepted");
 check(ed::close({{e("eps^-1"),e("0")},{e("0"),e("0")}}, {std::nullopt,-2},1)==ed::Highs({std::nullopt,-2}),"unreachable cycle demanded");
 check(ed::close({{e("0"),e("eps^-2")},{e("eps^2"),e("0")}}, {-3,std::nullopt},1)==ed::Highs({-3,-1}),"zero-weight cycle or negative Laurent closure");
 auto source=std::make_shared<const LaurentBoundary>(LaurentBoundary{-1,{{B(2),B(0),B(0),B(0),B(0)},{B(3),B(0),B(0),B(0),B(0)}},false});
 auto full=rows(2,2,-2,3);
 for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j)for(int k=-2;k<=3;++k)
   full.coefficients[i][j][k+2]=B(3*(i+1)+2*j+k);
 // A zero-containing, nonzero-radius first coefficient may not be stripped.
 full.coefficients[0][0][0]=B(0);arb_add_error_2exp_si(acb_realref(full.coefficients[0][0][0].raw()),-220);
 ed::Expression input{{full,source},{2,0}};
 // These stored cells are explicitly UNKNOWN, and must never enter arithmetic.
 for(unsigned j=0;j<2;++j){input.value.transform.coefficients[0][j][5]=B(9999);
   for(int k=1;k<=3;++k)input.value.transform.coefficients[1][j][k+2]=B(9999);}
 ExactEpsilonMatrix b{{e("eps^-1"),e("0")},{e("0"),e("eps^-2")}};
 ft::Options options;options.transport.taylor_order=20;
 auto component=ft::evolve_components(a,input,b,{e("0"),e("1/4")},{2,-1},{0,-2},options);
 check(component.initial_high==ed::Highs({2,0}),"component source demand raised unnecessarily");
 check(component.physical.high==ed::Highs({2,-1})&&component.integrated.high==ed::Highs({0,-2}),"component contracts lost");
 check(component.physical.value.leaf_source==source&&component.integrated.value.leaf_source==source,"immutable source replaced");
 auto legacy=ft::evolve(a,{full,source},b,{e("0"),e("1/4")},2,0,options);
 const B t=B(1)/B(4);
 for(unsigned j=0;j<2;++j) {
   for(unsigned i=0;i<2;++i)for(int k=component.physical.value.transform.low;k<=*component.physical.high[i];++k) {
     B expected=at(full,i,j,k);if(i==0)expected+=t*at(full,1,j,k-2);
     agrees(at(component.physical.value.transform,i,j,k),expected);
     agrees(at(component.physical.value.transform,i,j,k),at(legacy.physical.transform,i,j,k));
   }
   for(unsigned i=0;i<2;++i)for(int k=component.integrated.value.transform.low;k<=*component.integrated.high[i];++k) {
     B expected=i==0?t*at(full,0,j,k+1)+t*t*at(full,1,j,k-1)/B(2):t*at(full,1,j,k+2);
     agrees(at(component.integrated.value.transform,i,j,k),expected);
     agrees(at(component.integrated.value.transform,i,j,k),at(legacy.integrated.transform,i,j,k));
   }
 }
 const auto uncertain=at(component.physical.value.transform,0,0,-2);
 check(!uncertain.is_zero()&&acb_contains(uncertain.raw(),full.coefficients[0][0][0].raw()),"initial zero-containing ball was stripped or narrowed");
 auto missing=input;missing.high[1]=-1;
 rejects<ft::MapDemand>([&]{ft::evolve_components(a,missing,b,{e("0"),e("1/4")},{2,-1},{0,-2},options);},"unknown source upper tail silently zeroed");
 // A retained all-zero operator has UNKNOWN tail at high+1, not infinity.
 auto zero_map=rows(1,2,0,1);
 check(ed::pullback(zero_map,{3})==ed::Highs({1,1}),"zero retained operator treated as exact zero");
 // An all-zero retained input through -1 starts its unknown tail at zero.
 // A map known zero through 1 therefore has zero product through 1 without
 // needing its unknown coefficient at 2.
 auto padded_input=rows(2,2,-1,1);
 ed::Expression padded{{padded_input,source},{-1,-1}};
 auto padded_product=ed::compose(zero_map,padded,1);
 for(const auto& column:padded_product.transform.coefficients[0])for(const auto& value:column)
   check(value.is_zero(),"zero-prefix composition produced a nonzero coefficient");
 rejects<linear_boundary::CompositionDemand>([&]{ed::compose(zero_map,padded,2);},"stored zeros beyond component high treated as known");
 padded.value.transform.coefficients[0][0][0]=B(0);
 arb_add_error_2exp_si(acb_realref(padded.value.transform.coefficients[0][0][0].raw()),-220);
 rejects<linear_boundary::CompositionDemand>([&]{ed::compose(zero_map,padded,1);},"zero-containing input coefficient treated as exact zero");
 auto short_input=input;short_input.high={0,0};
 rejects<linear_boundary::CompositionDemand>([&]{ed::compose(zero_map,short_input,3);},"unknown operator tail accepted");
 auto consumer=rows(1,2,-1,5);consumer.coefficients[0][1][0]=B(1);
 rejects<linear_boundary::CompositionDemand>([&]{ed::compose(consumer,input,0);},"component compose consumed unknown upper coefficient");
 auto composed=ed::compose(consumer,input,-1);
 for(unsigned j=0;j<2;++j)for(int k=composed.transform.low;k<=-1;++k)
   agrees(at(composed.transform,0,j,k),at(full,1,j,k+1));
 // Fully exact zero prefixes may be dropped; preserving a nonzero epsilon
 // grade then requires undoing the diagonal shift on extraction.
 auto prefix=rows(2,2,-2,3);prefix.coefficients[0][0][3]=B(2);prefix.coefficients[1][1][2]=B(3);
 auto stripped=ft::evolve_components(a,{{prefix,source},{2,0}},b,{e("0"),e("1/4")},{2,-1},{0,-2},options);
 agrees(at(stripped.physical.value.transform,0,0,1),B(2));
 agrees(at(stripped.integrated.value.transform,1,1,-2),B(3)/B(4));
 // Negative physical connection valuations are legal when the demanded
 // subgraph admits the diagonal regularization; no epsilon grades are sampled.
 ExactEpsilonMatrix singular{{e("0"),e("eps^-2")},{e("0"),e("0")}};
 auto negative=ft::evolve_components(singular,{{full,source},{std::nullopt,-1}},{},
   {e("0"),e("1/4")},{-3,std::nullopt},{},options);
 for(unsigned j=0;j<2;++j)for(int k=-4;k<=-3;++k)
   agrees(at(negative.physical.value.transform,0,j,k),t*at(full,1,j,k+2));
 auto count_input=rows(1,2,0,2);count_input.coefficients[0][0][0]=B(1);
 auto count=ft::evolve_components({{e("0")}},{{count_input,source},{2}},{{e("eps^-2")}},
   {e("0"),e("1/4")},{0},{0},options);
 check(count.retained_slots==12,"packed storage counter omitted negative accumulator grades");
 // Complete nonleaf recursion comparison, including exact endpoint operators,
 // direct requests, dotted requests and separate observable epsilon poles.
 auto graph=recursion::prepare(feynman::example_family("sunrise"),
    {{1,1,1,0,0},{2,1,1,0,0},{0,1,1,0,0}},{},native_sunrise);
 const auto ei=path_epsilon_variables(graph.dimension).second;
 for(std::size_t i=0;i<graph.nodes[0].observable_rows.size();++i)
   for(auto& value:graph.nodes[0].observable_rows[i])
     value=epsilon_gauge_detail::multiply_power(value,ei,-static_cast<int>(i%3));
 recursion::NumericalOptions settings;settings.working_bits=256;settings.endpoint_order=24;
 settings.ordinary_order=64;settings.leaf_digits=18;settings.automatic_basis_scan=false;
 settings.linear_method=recursion::LinearMethod::factored;settings.component_epsilon_demands=false;
 recursion::Evaluator rectangular(graph,settings);auto old=rectangular.evaluate(0);
 settings.component_epsilon_demands=true;
 recursion::Evaluator scheduled(graph,settings);auto now=scheduled.evaluate(0);
 check(old.values.size()==now.values.size()&&scheduled.statistics().component_transports==2,
       "pipeline comparison bypassed component transport");
 for(std::size_t i=0;i<now.values.size();++i)for(int k=std::max(now.low,old.low);k<=0;++k)
   check(NativeTailMagnitude::upper_abs(now.values[i][k-now.low]-old.values[i][k-old.low]).approximate_upper()<1e-13,
         "component/rectangular recursive epsilon-pole observable disagreement");
 check(!now.taylor_tail_certified,"component pipeline asserted an omitted spatial tail certificate");
 std::cout<<"Component demands: asymmetric analytic/legacy agreement, negative Laurent poles, uncertainty, missing tails and cycle guards passed\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
