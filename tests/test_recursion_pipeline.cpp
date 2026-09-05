#include "diffexp/recursion_pipeline.hpp"
#include "diffexp/banana_oracle.hpp"
#include <iostream>
using namespace diffexp;using B=kernel::ComplexBall;
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
double upper(const B& b){return NativeTailMagnitude::upper_abs(b).approximate_upper();}
// Hermetic closure provider: native generated IBPs at the one nonleaf sunrise
// level. It constructs exactly the same ordered-coordinate contract as FIRE,
// without invoking an external binary or installing prepared coefficients.
level::Result native_sunrise(const ibp::PropagatorBasis& basis,const Exact& dimension,
    const ExactField&,std::size_t parameter,const std::vector<ibp::Integral>& requested,const level::Options&) {
  if(basis.physical_count!=2 || basis.denominators.size()!=5)throw std::invalid_argument("test closure only supports merged sunrise");
  ibp::Generator generator(basis,dimension);ibp::ExactReducer reducer(dimension,300);
  ibp::for_each_seed(2,5,{1,1,50},[&](const ibp::Integral& seed){for(auto& row:generator.relations(seed))reducer.insert(std::move(row));});
  ibp::BasisReduction coordinates(reducer,{{1,0,0,0,0},{1,1,0,0,0},{1,1,-1,0,0}},dimension);
  level::Result result;result.success=true;result.ordered_basis=coordinates.ordered_basis();
  for(const auto& master:result.ordered_basis)result.matrix.push_back(coordinates.resolve(generator.derivative(master,parameter)));
  for(const auto& target:requested)result.target_rows.push_back(coordinates.resolve({{target,dimension.constant(1)}}));
  return result;
}
template<class F>void fails(F f,const std::string& reason){try{f();}catch(const std::exception& e){require(std::string(e.what()).find(reason)!=std::string::npos,"unexpected recursive pipeline failure");return;}throw std::runtime_error("invalid recursive pipeline request accepted");}
int main(){try {
  B::set_precision(384);
  const auto endpoint_cache=std::filesystem::temp_directory_path()/("de3-pipeline-endpoints-"+std::to_string(::getpid()));
  struct Cleanup {std::filesystem::path path;~Cleanup(){std::error_code error;std::filesystem::remove_all(path,error);}} cleanup{endpoint_cache};
  require(recursion::NumericalOptions{}.observable_adjoint,"shared-source adjoint recursion must be the default");
  auto bubble=recursion::prepare(feynman::example_family("bubble"),{{1,1},{2,1},{1,0},{0,1},{-1,1},{0,0}});
  recursion::Evaluator bubble_solver(bubble);auto b=bubble_solver.evaluate(4);
  require(b.taylor_tail_certified && b.values.size()==6,"bubble graph must execute all scalar-leaf typed requests with certificates");
  oracle::BananaOptions oracle_options;oracle_options.target_bits=160;oracle_options.working_bits=320;
  auto bubble_reference=oracle::equal_banana_bessel(1,oracle_options).value;
  require(acb_contains(b.values[0][-b.low].raw(),bubble_reference.raw()),"bubble recursive graph must contain independent Bessel oracle");
  auto dotted=(bubble_reference+B(1))/B(5);
  require(acb_contains(b.values[1][-b.low].raw(),dotted.raw()),"dotted bubble request must contain mass-derivative identity");
  auto tadpole=feynman::tadpole(Rational(1),Rational(1),1,1,2,4);
  for(int k=-1;k<=4;++k)for(unsigned row:{2,3})
    require(acb_overlaps(b.values[row][k-b.low].raw(),tadpole.coefficients[k+1].raw()),"early one-sided bubble endpoint must equal tadpole boundary");
  for(const auto& value:b.values[5])require(value.is_zero(),"empty bubble sector must be structurally zero");
  auto before=bubble_solver.statistics();(void)bubble_solver.evaluate(2);
  require(bubble_solver.statistics().leaf_evaluations==before.leaf_evaluations && bubble_solver.statistics().cache_hits>before.cache_hits,"numeric leaf high-window cache not reused");

  auto graph=recursion::prepare(feynman::example_family("sunrise"),
    {{1,1,1,0,0},{2,1,1,0,0},{0,1,1,0,0}}, {},native_sunrise);
  std::vector<std::string> phases;recursion::NumericalOptions observed;
  observed.observable_adjoint=false;
  observed.endpoint_cache_directory=endpoint_cache/"sunrise";
  observed.progress=[&](std::size_t depth,const std::string& phase,int high){require(depth<graph.nodes.size() && high<=100,"bounded progress callback context");phases.push_back(phase);};
  // Geometry must retain a pole hidden by a numerator cancellation at eps=0.
  auto near=graph;near.nodes[0].operations.resize(1);near.nodes[0].operations[0]=pullback::Plan{};
  near.nodes[0].operations[0].operation=feynman::Operation::BetaIntegral;
  near.nodes[0].requested.resize(1);near.nodes[0].observable_rows={std::vector<Exact>(3,near.dimension.constant(0))};
  near.nodes[0].observable_rows[0][0]=near.dimension.constant(1);
  for(auto& row:near.nodes[0].closure.matrix)for(auto& entry:row)entry=near.dimension.constant(0);
  auto path=near.dimension.variable(0),epsilon=near.dimension.variable(1);
  near.nodes[0].closure.matrix[0][0]=(path.constant(1)-path.constant(100)*path+epsilon)/(path.constant(1)-path.constant(100)*path);
  recursion::NumericalOptions geometry_only;geometry_only.endpoint_order=8;
  geometry_only.endpoint_cache_directory=endpoint_cache/"geometry";
  geometry_only.progress=[](std::size_t,const std::string& phase,int){if(phase.starts_with("endpoint overlap="))throw std::runtime_error("geometry ready");};
  recursion::Evaluator near_solver(near,geometry_only);
  fails([&]{near_solver.evaluate(0);},"geometry ready");
  require(near_solver.statistics().endpoint_series_built>0,"interrupted endpoint preparation must save exact work");
  recursion::Evaluator resumed_near(near,geometry_only);
  fails([&]{resumed_near.evaluate(0);},"geometry ready");
  require(resumed_near.statistics().endpoint_series_built==0 && resumed_near.statistics().endpoint_series_reused==2,
      "new evaluator must verify and reuse both exact endpoints after later interruption");
  const auto& clearance=near_solver.endpoint_geometry(0);
  require(clearance.overlap<=Rational("1/1600") && clearance.nearest_pole_lower.approximate_upper()>0.009999 && clearance.nearest_pole_lower.approximate_upper()<0.010001,"endpoint overlap must resolve poles hidden at epsilon zero");
  recursion::Evaluator solver(graph,observed);auto result=solver.evaluate(0);
  const auto& geometry=solver.endpoint_geometry(0);
  require(geometry.overlap<=observed.overlap && NativeTailMagnitude::upper_abs(B::from_strings(geometry.overlap.str()))*NativeTailMagnitude::from_ui(16)<=geometry.nearest_pole_lower,"geometric endpoint margin must be proven from all nonzero pole enclosures");
  require(!result.taylor_tail_certified && result.values.size()==3,"nonleaf finite endpoint/ordinary series must remain explicitly uncertified");
  auto reference=oracle::equal_banana_bessel(2).value;
  require(upper(result.values[0][-result.low]-reference)<1e-20,"generic sunrise FT recursion disagrees with independent Bessel oracle");
  for(int k=result.low;k<0;++k)require(upper(result.values[0][k-result.low])<1e-20,"physical sunrise spurious Laurent poles must cancel");
  recursion::NumericalOptions adjoint;adjoint.observable_adjoint=true;
  adjoint.endpoint_cache_directory=observed.endpoint_cache_directory;
  recursion::Evaluator adjoint_solver(graph,adjoint);auto adjoint_result=adjoint_solver.evaluate(0);
  require(solver.statistics().endpoint_series_built>0 && adjoint_solver.statistics().endpoint_series_built==0 &&
      adjoint_solver.statistics().endpoint_series_reused==2,"exact series cache must survive evaluator and numerical-method changes");
  require(!adjoint_result.taylor_tail_certified,"adjoint finite-series result must remain uncertified");
  for(unsigned i=0;i<result.values.size();++i)for(int k=std::max(result.low,adjoint_result.low);k<=0;++k)
    require(upper(result.values[i][k-result.low]-adjoint_result.values[i][k-adjoint_result.low])<1e-20,"adjoint and value recursion must agree for scalar/dotted/endpoint sunrise requests");
  require(upper(adjoint_result.values[0][-adjoint_result.low]-reference)<1e-20,"adjoint sunrise must agree with independent Bessel oracle");
  auto interrupted_options=adjoint;
  interrupted_options.ordinary_cache_directory=endpoint_cache/"ordinary";
  bool upper_arm=false;
  interrupted_options.progress=[&](std::size_t,const std::string& phase,int) {upper_arm=phase=="upper adjoint transport";};
  interrupted_options.adjoint.chart_observer=[&](unsigned,double,double,const LaurentRows&) {
    if(upper_arm)throw std::runtime_error("intentional upper arm interruption");
  };
  recursion::Evaluator interrupted(graph,interrupted_options);
  fails([&]{interrupted.evaluate(0);},"intentional upper arm interruption");
  require(interrupted.statistics().ordinary_checkpoints.saved>1,"evaluator failed to persist accepted ordinary charts");
  interrupted_options.progress={};interrupted_options.adjoint.chart_observer={};
  recursion::Evaluator resumed(graph,interrupted_options);const auto resumed_result=resumed.evaluate(0);
  require(resumed.statistics().ordinary_checkpoints.loaded==2 && resumed.statistics().ordinary_checkpoints.completed_reused==1,
    "new evaluator did not reuse completed lower and partial upper arms");
  for(unsigned i=0;i<adjoint_result.values.size();++i)for(unsigned k=0;k<adjoint_result.values[i].size();++k)
    require(acb_equal(adjoint_result.values[i][k].raw(),resumed_result.values[i][k].raw()),"evaluator continuation changed retained result bits");
  require(resumed.linear_expression().leaf_source==resumed.linear_expression(1).leaf_source,
    "evaluator continuation lost common immutable leaf identity");
  // Deleting the first propagator factorizes into two gamma tadpoles.
  Jet eps(0,3,384);eps.set(1,B(1));auto gamma=(eps.constant(1)+eps).gamma();auto product=gamma*gamma;
  for(int k=-2;k<=0;++k)require(upper(result.values[2][k-result.low]-product.at(k+2))<1e-20,"nonleaf lower endpoint must reproduce factorized gamma product");
  require(result.values[1][-result.low].is_finite(),"dotted sunrise shared accumulator must be finite");
  require(std::find(phases.begin(),phases.end(),"endpoint preparation")!=phases.end() && std::find(phases.begin(),phases.end(),"child refinement")!=phases.end() && std::find(phases.begin(),phases.end(),"upper matching")!=phases.end(),"progress callback must expose preparation, refinement and matching");
  require(solver.statistics().exact_plans==1 && solver.statistics().refinements>0,"endpoint plans must survive numeric-only epsilon refinement");
  before=solver.statistics();(void)solver.evaluate(0);
  require(solver.statistics().exact_plans==before.exact_plans && solver.statistics().numeric_evaluations==before.numeric_evaluations,"completed numerical result cache not reused");
  // Direct identities do not need an endpoint expansion at all.
  auto direct=graph;direct.nodes[0].operations.resize(1);direct.nodes[0].operations[0].operation=feynman::Operation::Direct;
  direct.nodes[0].requested.resize(1);direct.nodes[0].observable_rows={std::vector<Exact>(3,direct.dimension.constant(0))};
  recursion::Evaluator direct_solver(direct);auto direct_result=direct_solver.evaluate(0);
  require(direct_solver.statistics().exact_plans==0 && direct_result.values[0][0].is_zero(),"direct identity should bypass endpoint planning");
  // Three-level direct graph: the middle node duplicates one leaf source.
  // Materializing those two rows independently destroys a shared cancellation.
  auto correlated=graph;correlated.nodes.insert(correlated.nodes.begin()+1,graph.nodes[0]);
  auto& first=correlated.nodes[0];auto& middle_node=correlated.nodes[1];
  first.requested.resize(1);first.operations.resize(1);first.operations[0].operation=feynman::Operation::Direct;
  const auto big=correlated.dimension.constant(Rational("1000000000000000000000000000000000000000000000000000000000000"));
  auto exact_zero=correlated.dimension.constant(0),exact_one=correlated.dimension.constant(1);
  first.observable_rows={{big+exact_one,-big,exact_zero}};
  middle_node.requested=first.closure.ordered_basis;
  middle_node.operations=std::vector<pullback::Plan>(3);for(auto& operation:middle_node.operations)operation.operation=feynman::Operation::Direct;
  middle_node.observable_rows={{exact_one,exact_zero,exact_zero},{exact_one,exact_zero,exact_zero},{exact_zero,exact_zero,exact_one}};
  recursion::NumericalOptions values_options;values_options.observable_adjoint=false;
  recursion::Evaluator global(correlated,adjoint),independent(correlated,values_options);
  auto correlated_value=global.evaluate(0),independent_value=independent.evaluate(0),leaf_value=global.evaluate(2,0);
  const auto expected_leaf=leaf_value.values[0][-leaf_value.low];
  const double global_error=upper(correlated_value.values[0][-correlated_value.low]-expected_leaf);
  const double independent_error=upper(independent_value.values[0][-independent_value.low]-expected_leaf);
  require(global_error<1e-25 && independent_error>global_error*1e15,"global expression must preserve correlation across distinct intermediate output rows");
  auto global_before=global.statistics();(void)global.evaluate(0);
  require(global.statistics().numeric_evaluations==global_before.numeric_evaluations,"global expression materialized cache remains synchronized");
  (void)global.evaluate(1);auto refreshed=global.evaluate(0);
  require(upper(refreshed.values[0][-refreshed.low]-expected_leaf)<1e-25,"global expression remains synchronized after source-window refinement");
  auto zero_chain=correlated;zero_chain.nodes[1].closure.ordered_basis.clear();
  recursion::Evaluator zero_expression(zero_chain,adjoint);auto zero_result=zero_expression.evaluate(0);
  for(const auto& coefficient:zero_result.values[0])require(coefficient.is_zero(),"exact-zero child must retain a composable synchronized expression");
  recursion::NumericalOptions limits=values_options;limits.max_epsilon=1;
  fails([&]{recursion::Evaluator limited(graph,limits);limited.evaluate(0);},"epsilon demand exceeds finite budget");
  fails([&]{solver.evaluate(101);},"epsilon budget");
  require(B::precision()==384,"recursive evaluator must restore ambient precision");
  std::cout<<"Generic FT pipeline: certified bubble leaf/typed requests; sunrise scalar/dotted/endpoint recursion versus independent oracle; exact-plan reuse and finite epsilon demands passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
