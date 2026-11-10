#include "diffexp/recursion_pipeline.hpp"
#include "diffexp/banana_oracle.hpp"
#include "diffexp/family_config.hpp"
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
  const auto scratch_directory=std::filesystem::temp_directory_path()/("de3-pipeline-endpoints-"+std::to_string(::getpid()));
  struct Cleanup {std::filesystem::path path;~Cleanup(){std::error_code error;std::filesystem::remove_all(path,error);}} cleanup{scratch_directory};
  require(recursion::NumericalOptions{}.linear_method==recursion::LinearMethod::adjoint,"shared-source adjoint recursion must be the default");
  {
    // The option contract is checkable without a graph, i.e. before the exact
    // preparation a front end would otherwise run first.
    recursion::validate_numerical_options(recursion::NumericalOptions{});
    recursion::NumericalOptions invalid;invalid.working_bits=32;
    fails([&]{recursion::validate_numerical_options(invalid);},"working precision");
    invalid=recursion::NumericalOptions{};invalid.ordinary_order=4;
    fails([&]{recursion::validate_numerical_options(invalid);},"ordinary order");
    invalid=recursion::NumericalOptions{};invalid.overlap=Rational("3/4");
    fails([&]{recursion::validate_numerical_options(invalid);},"overlap");
    // Only a document that names epsilon_order has chosen it; a bare example
    // name carries the template default.
    auto document=family_config::describe(feynman::example_family("bubble")).as_object();
    require(family_config::parse(document).epsilon_order_given && family_config::parse(document).epsilon_order==4,"configured epsilon order not recognised");
    document.erase("epsilon_order");
    require(!family_config::parse(document).epsilon_order_given && !family_config::load("bubble").epsilon_order_given &&
            family_config::load("bubble").epsilon_order==4,"absent epsilon order reported as chosen");
  }
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
  geometry_only.automatic_basis_scan=false;
  // Abort after the immutable plan commits, not from a construction callback
  // before endpoint_geometry() is available.
  geometry_only.progress=[](std::size_t,const std::string& phase,int){if(phase=="child boundary")throw std::runtime_error("geometry ready");};
  recursion::Evaluator near_solver(near,geometry_only);
  fails([&]{near_solver.evaluate(0);},"geometry ready");
  require(near_solver.statistics().endpoint_series_built==2,"both endpoint series are built exactly once per plan");
  const auto& clearance=near_solver.endpoint_geometry(0);
  require(clearance.overlap<=Rational("1/1600") && clearance.nearest_pole_lower.approximate_upper()>0.009999 && clearance.nearest_pole_lower.approximate_upper()<0.010001,"endpoint overlap must resolve poles hidden at epsilon zero");
  // Compare the new endpoint route in original transport coordinates, including
  // the epsilon gauge, coupled modes, and automatic fallback on the other end.
  auto functional_graph=near;
  for(auto& row:functional_graph.nodes[0].closure.matrix)for(auto& entry:row)entry=path.constant(0);
  for(unsigned i=0;i<3;++i)functional_graph.nodes[0].closure.matrix[i][i]=path.constant(i+1)*epsilon/path;
  functional_graph.nodes[0].closure.matrix[0][1]=path.constant(1)/path+path/(path.constant(1)-path);
  functional_graph.nodes[0].closure.matrix[1][2]=epsilon/(path.constant(1)-path);
  recursion::NumericalOptions full_endpoint;full_endpoint.endpoint_order=24;
  recursion::Evaluator full_endpoint_solver(functional_graph,full_endpoint);
  auto full_stage=full_endpoint_solver.prepare_adjoint_stage(0,4);
  auto direct_endpoint=full_endpoint;direct_endpoint.direct_endpoint_functionals=true;
  recursion::Evaluator direct_endpoint_solver(functional_graph,direct_endpoint);
  auto direct_stage=direct_endpoint_solver.prepare_adjoint_stage(0,4);
  require(direct_endpoint_solver.statistics().direct_endpoints==1 &&
      direct_endpoint_solver.statistics().endpoint_series_built==1,
      "direct endpoint should skip the lower basis and retain the fixed-sector upper basis");
  require(full_stage.lower_path==direct_stage.lower_path && full_stage.upper_path==direct_stage.upper_path,
      "direct endpoint route changed contour geometry");
  for(const auto pair:{std::pair{&full_stage.lower_endpoint,&direct_stage.lower_endpoint},
      std::pair{&full_stage.upper_endpoint,&direct_stage.upper_endpoint}})
    for(unsigned i=0;i<pair.first->coefficients.size();++i)for(unsigned j=0;j<pair.first->columns();++j)
      for(int k=std::min(pair.first->low,pair.second->low);k<=4;++k){
        const auto a=k<pair.first->low?B(0):pair.first->coefficients[i][j][k-pair.first->low];
        const auto b=k<pair.second->low?B(0):pair.second->coefficients[i][j][k-pair.second->low];
        require(upper(a-b)<1e-24,"direct and full endpoint maps disagree in original coordinates");
      }
  {
    // Mixed endpoints: the lower one is direct (one row per operation), the
    // upper series contributes a zero-constraint row for the fixed 1/(1-x)
    // weight. Both operators must expose the same row count to the demand
    // propagation and to transport.
    auto constrained=functional_graph;
    constrained.nodes[0].observable_rows[0][2]=path.constant(1)/(path.constant(1)-path);
    auto mixed=direct_endpoint;mixed.automatic_basis_scan=false;unsigned constraints=0;
    mixed.progress=[&](std::size_t,const std::string& phase,int){if(phase.rfind("physical endpoint zero constraints",0)==0)++constraints;};
    recursion::Evaluator mixed_solver(constrained,mixed);
    const auto profile=mixed_solver.epsilon_demand_profile(0);
    require(constraints==1 && mixed_solver.statistics().direct_endpoints==1 && mixed_solver.statistics().endpoint_series_built==1,
        "mixed direct/series endpoint fixture lost its constraint row");
    require(profile.child_component_loss.size()==3 && mixed_solver.statistics().demand_preflights==1,
        "direct endpoint operator lacks the other endpoint's constraint rows");
  }
  recursion::Evaluator solver(graph,observed);auto result=solver.evaluate(0);
  const auto& geometry=solver.endpoint_geometry(0);
  require(geometry.overlap<=observed.overlap && NativeTailMagnitude::upper_abs(B::from_strings(geometry.overlap.str()))*NativeTailMagnitude::from_ui(16)<=geometry.nearest_pole_lower,"geometric endpoint margin must be proven from all nonzero pole enclosures");
  require(!result.taylor_tail_certified && result.values.size()==3,"nonleaf finite endpoint/ordinary series must remain explicitly uncertified");
  auto reference=oracle::equal_banana_bessel(2).value;
  require(upper(result.values[0][-result.low]-reference)<1e-20,"generic sunrise FT recursion disagrees with independent Bessel oracle");
  for(int k=result.low;k<0;++k)require(upper(result.values[0][k-result.low])<1e-20,"physical sunrise spurious Laurent poles must cancel");
  {
    // The same recursion expanded around epsilon=1/64: its order-zero value is
    // the Laurent series at epsilon=0 summed at 1/64, and no level needs more
    // than order zero of its child.
    recursion::NumericalOptions plain;plain.working_bits=384;
    recursion::Evaluator series_solver(graph,plain);const auto series=series_solver.evaluate(10);
    const Rational shift("1/64");
    const auto shifted_graph=recursion::shift_epsilon(graph,shift);
    auto shifted_options=plain;shifted_options.shifted_leaf_orders=24;
    recursion::Evaluator shifted_solver(shifted_graph,shifted_options);const auto sample=shifted_solver.evaluate(0);
    require(sample.low==0 && sample.high()==0,"a shifted recursion has no epsilon poles and needs no wider window");
    B sum(0),power(1),at=B::from_strings(shift.str());
    for(int k=series.low;k<0;++k)power=power/at;
    for(int k=series.low;k<=series.high();++k){sum=sum+series.values[0][k-series.low]*power;power=power*at;}
    require(upper(sample.values[0][0]-sum)<1e-15,"shifted-epsilon sunrise disagrees with its summed Laurent series");
  }
  // A second, independent linear route: the forward (factored) maps must
  // reproduce the adjoint result for scalar, dotted and endpoint requests.
  recursion::NumericalOptions adjoint;
  auto factored=adjoint;factored.linear_method=recursion::LinearMethod::factored;
  recursion::Evaluator adjoint_solver(graph,adjoint),factored_solver(graph,factored);
  auto adjoint_result=adjoint_solver.evaluate(0);const auto factored_result=factored_solver.evaluate(0);
  require(adjoint_solver.statistics().endpoint_series_built==2 && adjoint_solver.statistics().adjoint_selections==1 &&
      factored_solver.statistics().factored_selections==1,"each route builds its two endpoint series once");
  require(!adjoint_result.taylor_tail_certified,"adjoint finite-series result must remain uncertified");
  for(unsigned i=0;i<adjoint_result.values.size();++i)for(int k=std::max(factored_result.low,adjoint_result.low);k<=0;++k)
    require(upper(factored_result.values[i][k-factored_result.low]-adjoint_result.values[i][k-adjoint_result.low])<1e-20,"adjoint and factored recursion must agree for scalar/dotted/endpoint sunrise requests");
  require(upper(adjoint_result.values[0][-adjoint_result.low]-reference)<1e-20,"adjoint sunrise must agree with independent Bessel oracle");
  require(upper(factored_result.values[0][-factored_result.low]-reference)<1e-20,"factored sunrise must agree with independent Bessel oracle");
  auto coefficient_options=adjoint;
  unsigned progress_events=0;
  coefficient_options.ordinary_fraction_progress=[&](double fraction) {
    require(fraction>=0 && fraction<=1,"ordinary progress fraction outside original path");++progress_events;
  };
  recursion::Evaluator coefficient_solver(graph,coefficient_options);
  const auto coefficient_result=coefficient_solver.evaluate(0);
  require(progress_events>2,"ordinary transport reported no path progress");
  for(unsigned i=0;i<adjoint_result.values.size();++i)for(unsigned k=0;k<adjoint_result.values[i].size();++k)
    require(acb_equal(adjoint_result.values[i][k].raw(),coefficient_result.values[i][k].raw()),"progress reporting changed retained result bits");
  auto interrupted_options=adjoint;
  interrupted_options.ordinary_cache_directory=scratch_directory/"ordinary";
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
  require(std::find(phases.begin(),phases.end(),"endpoint preparation")!=phases.end() && std::find(phases.begin(),phases.end(),"child boundary")!=phases.end() && std::find(phases.begin(),phases.end(),"upper adjoint transport")!=phases.end(),"progress callback must expose preparation, child demand and transport");
  // The structural preflight makes the first child request sufficient, and a
  // later higher epsilon order extends the same exact plan.
  require(solver.statistics().exact_plans==1 && solver.statistics().refinements==0,"child window preflight did not avoid a refinement");
  const auto higher=solver.evaluate(1);
  require(solver.statistics().exact_plans==1 && solver.statistics().endpoint_series_built==2 &&
      upper(higher.values[0][-higher.low]-reference)<1e-20,"a higher epsilon request must reuse the exact endpoint plan");
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
  recursion::Evaluator global(correlated,adjoint);
  auto correlated_value=global.evaluate(0),leaf_value=global.evaluate(2,0);
  const auto expected_leaf=leaf_value.values[0][-leaf_value.low];
  // Rows evaluated independently would carry the leaf's own uncertainty times
  // 1e60; the shared expression cancels the two huge coefficients exactly.
  const double global_error=upper(correlated_value.values[0][-correlated_value.low]-expected_leaf);
  require(global_error<1e-25,"global expression must preserve correlation across distinct intermediate output rows");
  auto global_before=global.statistics();(void)global.evaluate(0);
  require(global.statistics().numeric_evaluations==global_before.numeric_evaluations,"global expression materialized cache remains synchronized");
  (void)global.evaluate(1);auto refreshed=global.evaluate(0);
  require(upper(refreshed.values[0][-refreshed.low]-expected_leaf)<1e-25,"global expression remains synchronized after source-window refinement");
  auto zero_chain=correlated;zero_chain.nodes[1].closure.ordered_basis.clear();
  recursion::Evaluator zero_expression(zero_chain,adjoint);auto zero_result=zero_expression.evaluate(0);
  for(const auto& coefficient:zero_result.values[0])require(coefficient.is_zero(),"exact-zero child must retain a composable synchronized expression");
  // A Direct pole of depth two needs two child orders beyond the budget.
  recursion::NumericalOptions limits;limits.max_epsilon=1;auto deep=direct;
  deep.nodes[0].observable_rows[0][0]=exact_one/correlated.dimension.variable(1).pow(2);
  fails([&]{recursion::Evaluator limited(deep,limits);limited.evaluate(0);},"epsilon demand exceeds finite budget");
  fails([&]{solver.evaluate(101);},"epsilon budget");
  require(B::precision()==384,"recursive evaluator must restore ambient precision");
  std::cout<<"Generic FT pipeline: certified bubble leaf/typed requests; sunrise scalar/dotted/endpoint recursion versus independent oracle; exact-plan reuse and finite epsilon demands passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
