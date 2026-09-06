#include "diffexp/recursion_pipeline.hpp"
#include <iostream>
using namespace diffexp;using B=Jet::Ball;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
double upper(const B& value){return NativeTailMagnitude::upper_abs(value).approximate_upper();}
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
int main(){try{
  B::set_precision(256);
  auto graph=recursion::prepare(feynman::example_family("sunrise"),
      {{1,1,1,0,0},{2,1,1,0,0},{0,1,1,0,0}},{},native_sunrise);
  recursion::NumericalOptions settings;settings.working_bits=256;settings.endpoint_order=24;
  settings.ordinary_order=64;settings.leaf_digits=18;
  require(settings.linear_method==recursion::LinearMethod::adjoint,"linear method default changed");
  recursion::Evaluator baseline(graph,settings);auto expected=baseline.evaluate(0);
  settings.linear_method=recursion::LinearMethod::factored;
  recursion::Evaluator factored(graph,settings);auto actual=factored.evaluate(0);
  for(unsigned i=0;i<actual.values.size();++i)for(int k=std::max(actual.low,expected.low);k<=0;++k)
    require(upper(actual.values[i][k-actual.low]-expected.values[i][k-expected.low])<1e-14,
        "factored/adjoint scalar, dotted or endpoint sunrise mismatch");
  require(factored.statistics().factored_selections>0 && factored.statistics().adjoint_selections==0,
      "forced factored route was not selected");
  require(factored.statistics().exact_plans==1 && !actual.taylor_tail_certified,"factored plan reuse/tail contract");
  auto before=factored.statistics();(void)factored.evaluate(0);
  require(factored.statistics().numeric_evaluations==before.numeric_evaluations,"factored expression/value cache desynchronized");
  settings.linear_method=recursion::LinearMethod::automatic;
  recursion::Evaluator automatic(graph,settings);auto chosen=automatic.evaluate(0);
  require(automatic.statistics().factored_selections+automatic.statistics().adjoint_selections==1,
      "automatic method selection not recorded");
  require(upper(chosen.values[0][-chosen.low]-actual.values[0][-actual.low])<1e-14,"automatic route mismatch");
  // Duplicate a nonleaf scalar observable, then cancel large coefficients at
  // its parent. The factored middle retains one immutable leaf expression.
  auto correlated=graph;
  correlated.nodes[0].requested[1]=correlated.nodes[0].requested[0];
  correlated.nodes[0].operations[1]=correlated.nodes[0].operations[0];
  correlated.nodes[0].observable_rows[1]=correlated.nodes[0].observable_rows[0];
  auto parent=correlated.nodes[0];parent.closure.ordered_basis=parent.requested;
  parent.requested.resize(1);parent.operations.resize(1);parent.operations[0].operation=feynman::Operation::Direct;
  auto zero=graph.dimension.constant(0),one=graph.dimension.constant(1),big=graph.dimension.constant(Rational("100000000000000000000"));
  parent.observable_rows={{big+one,-big,zero}};
  correlated.nodes.insert(correlated.nodes.begin(),std::move(parent));
  settings.linear_method=recursion::LinearMethod::factored;
  const auto cache=std::filesystem::temp_directory_path()/("de3-factored-continuation-"+std::to_string(::getpid()));
  struct Cleanup {std::filesystem::path path;~Cleanup(){std::error_code ec;std::filesystem::remove_all(path,ec);}} cleanup{cache};
  // Explicitly exercise local-chart interruption and exact resumption.
  settings.ordinary_method=recursion::OrdinaryMethod::taylor;
  settings.endpoint_cache_directory=cache/"endpoints";settings.ordinary_cache_directory=cache/"ordinary";
  settings.progress=[](std::size_t,const std::string& phase,int) {
    if(phase=="middle factored transport")throw std::runtime_error("intentional factored interruption");
  };
  recursion::Evaluator interrupted(correlated,settings);bool stopped=false;
  try{interrupted.evaluate(0);}catch(const std::runtime_error& e){stopped=std::string(e.what())=="intentional factored interruption";}
  require(stopped && interrupted.statistics().ordinary_checkpoints.saved>0,"factored ordinary adapter failed to save completed lower map");
  settings.progress={};
  recursion::Evaluator shared(correlated,settings);auto cancellation=shared.evaluate(0);
  require(shared.statistics().ordinary_checkpoints.completed_reused==1,"factored evaluator repeated completed lower map");
  require(shared.linear_expression().leaf_source==shared.linear_expression(2).leaf_source,"factored continuation lost common leaf identity");
  require(upper(cancellation.values[0][-cancellation.low]-actual.values[0][-actual.low])<1e-14,
      "factored recursion lost correlated child-source cancellation");
  require(shared.statistics().factored_selections>0,"correlated regression bypassed factored integration");
  settings.ordinary_cache_directory.clear();settings.endpoint_cache_directory.clear();
  // A cost estimate can select a poorly conditioned direction. Inject that
  // bounded failure after its lower transport, with a completed cached leaf.
  // Automatic mode must recover once; forced mode and unrelated errors must
  // remain visible. Repeated order requests must not restart the failed route.
  auto wide=graph;
  for(unsigned i=0;i<16;++i) {
    wide.nodes[0].requested.push_back(graph.nodes[0].requested[0]);
    wide.nodes[0].operations.push_back(graph.nodes[0].operations[0]);
    wide.nodes[0].observable_rows.push_back(graph.nodes[0].observable_rows[0]);
  }
  unsigned failures=0;
  settings.linear_method=recursion::LinearMethod::automatic;
  settings.progress=[&](std::size_t,const std::string& phase,int) {
    if(phase=="middle factored transport"){++failures;throw ArithmeticConditioningFailure("injected bounded conditioning failure");}
  };
  recursion::Evaluator recovered(wide,settings);auto recovery=recovered.evaluate(0);
  require(failures==1 && recovered.statistics().conditioning_method_fallbacks==1 &&
    recovered.statistics().factored_selections==1 && recovered.statistics().adjoint_selections==1,
    "automatic method did not perform exactly one conditioning fallback");
  require(recovered.statistics().exact_plans==1 && recovered.statistics().leaf_evaluations==1,
    "conditioning fallback repeated exact preparation or the leaf evaluation");
  require(upper(recovery.values[0][-recovery.low]-expected.values[0][-expected.low])<1e-14,
    "conditioning fallback changed the independent adjoint result");
  (void)recovered.evaluate(1);
  require(failures==1 && recovered.statistics().factored_selections==1,
    "higher-order request retried a failed method on the same plan");
  settings.linear_method=recursion::LinearMethod::factored;bool forced_failed=false;
  try {recursion::Evaluator forced(wide,settings);forced.evaluate(0);}catch(const ArithmeticConditioningFailure&){forced_failed=true;}
  require(forced_failed,"forced factored mode silently changed method");
  settings.linear_method=recursion::LinearMethod::automatic;
  settings.progress=[](std::size_t,const std::string& phase,int) {
    if(phase=="middle factored transport")throw std::runtime_error("unrelated failure");
  };
  bool unrelated=false;
  try {recursion::Evaluator invalid(wide,settings);invalid.evaluate(0);}catch(const ArithmeticConditioningFailure&){throw;}
  catch(const std::runtime_error& error){unrelated=std::string(error.what())=="unrelated failure";}
  require(unrelated,"automatic method hid an unrelated failure");
  settings.progress=[](std::size_t,const std::string& phase,int) {
    if(phase=="middle factored transport") {
      ArithmeticConditioningFailure failure("failed descendant");failure.recursion_depth=1;throw failure;
    }
  };
  recursion::Evaluator dependent(wide,settings);bool descendant=false;
  try {dependent.evaluate(0);}catch(const ArithmeticConditioningFailure& failure){descendant=failure.recursion_depth==1;}
  require(descendant && dependent.statistics().conditioning_method_fallbacks==0 &&
    dependent.statistics().adjoint_selections==0,"parent retried transport after a descendant failed");
  std::cout<<"Factored recursion: scalar/dotted/endpoint sunrise equivalence, automatic choice, plan/cache reuse and shared-source cancellation passed\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
