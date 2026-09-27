#include "diffexp/recursion_pipeline.hpp"
#include <iostream>
using namespace diffexp;
namespace bp=recursion::basis_persistence;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
template<class F>void rejects(F f,const char* message){try{f();}catch(const std::exception&){return;}throw std::runtime_error(message);}
int main(){try{
  const auto root=std::filesystem::temp_directory_path()/("diffexp-basis-persistence-"+std::to_string(::getpid()));
  struct Cleanup{std::filesystem::path root;~Cleanup(){std::error_code e;std::filesystem::remove_all(root,e);}}cleanup{root};
  auto graph=recursion::prepare(feynman::example_family("bubble"),{{1,1}});
  auto& node=graph.nodes.front();
  // Synthetic single-master closure: no external reduction or ordinary solve.
  node.scalar_leaf=false;
  ibp::Integral master(node.merged.denominators.size(),0);master[0]=1;
  node.closure.ordered_basis={master};node.closure.success=true;
  const auto one=graph.dimension.constant(1),zero=graph.dimension.constant(0),x=graph.dimension.variable(0);
  node.closure.matrix={{one/(one+x)}};
  node.observable_rows={{one}};node.operations={pullback::Plan{}};
  causal::Prescription prescription;prescription.f_rim=-1;prescription.levels={{1}};prescription.provenance="synthetic exact cache test";
  recursion::NumericalOptions options;options.endpoint_order=4;options.working_bits=192;
  options.basis_scan_seconds=0.001;options.causal_prescription=prescription;
  options.endpoint_cache_directory=root/"endpoints";
  require(bp::directory(options)==options.endpoint_cache_directory/"physical-bases-v1","derived basis directory");
  options.basis_cache_directory=root/"selected";
  require(bp::directory(options)==options.basis_cache_directory,"explicit basis directory");
  auto id=bp::identity(graph,0,0,1,prescription);
  artifacts::Store store(options.basis_cache_directory);
  require(!bp::load(store,id,node.closure.matrix,0),"unexpected initial cache hit");
  recursion::PhysicalBasisTransform transform{node.closure.matrix,{{one/(one+x)}},{{one+x}},{{zero}}};
  bp::Choice choice{transform,true};
  bp::save(store,id,choice,node.closure.matrix,0);
  auto loaded=bp::load(store,id,node.closure.matrix,0);
  require(loaded && loaded->selected && loaded->basis.to_basis==transform.to_basis &&
          loaded->basis.from_basis==transform.from_basis && loaded->basis.connection==transform.connection,"exact basis round trip");
  const auto unit=fuchsify::detail::identity(1,one);
  bp::save(store,id,{{node.closure.matrix,unit,unit,node.closure.matrix},false},node.closure.matrix,0);
  require(bp::load(store,id,node.closure.matrix,0)->selected,"completed selection was replaced");
  auto corrupted=bp::payload(choice);
  corrupted["from_basis"]=cached_level::detail::matrix_json(ExactEpsilonMatrix{{one}});
  rejects([&]{bp::decode(corrupted,node.closure.matrix,0);},"false inverse accepted");
  corrupted=bp::payload(choice);corrupted["connection"]=cached_level::detail::matrix_json(ExactEpsilonMatrix{{one}});
  rejects([&]{bp::decode(corrupted,node.closure.matrix,0);},"false gauge identity accepted");
  corrupted=bp::payload(choice);corrupted["selected"]=false;
  rejects([&]{bp::decode(corrupted,node.closure.matrix,0);},"nonidentity original decision accepted");
  corrupted=bp::payload(choice);corrupted["candidate_child_loss"]=0;
  rejects([&]{bp::decode(corrupted,node.closure.matrix,0);},"cached order-specific proof accepted");
  // A valid artifact checksum/certificate is not an algebraic verifier.
  auto hostile_id=id;hostile_id.algorithm_version="tampered-payload-test";
  auto bad_payload=bp::payload(choice);
  bad_payload["connection"]=cached_level::detail::matrix_json(ExactEpsilonMatrix{{one}});
  store.put(hostile_id,Demand{0,0,0,64,0},bad_payload,
    artifacts::Certificate{"exact",bp::verifier,bp::scope,{{"claim","deliberately false gauge for test"}}});
  rejects([&]{bp::load(store,hostile_id,node.closure.matrix,0);},"cache load trusted a false exact certificate");
  auto changed=graph;changed.nodes[0].closure.matrix={{zero}};
  require(bp::identity(changed,0,0,1,prescription).key()!=id.key(),"original connection missing from identity");
  changed=graph;changed.nodes[0].closure.ordered_basis[0][0]=2;
  require(bp::identity(changed,0,0,1,prescription).key()!=id.key(),"ordered master identity missing");
  changed=graph;changed.nodes[0].observable_rows={{one+one}};
  require(bp::identity(changed,0,0,1,prescription).key()!=id.key(),"observable identity missing");
  auto changed_prescription=prescription;changed_prescription.levels[0].x_detour_sign=-1;
  require(bp::identity(graph,0,0,1,changed_prescription).key()!=id.key(),"domain identity missing");
  auto revision=id;revision.algorithm_version="different-scan-revision";
  require(!bp::load(store,revision,node.closure.matrix,0),"scan revision ignored");
  recursion::Evaluator resumed(graph,options);
  (void)resumed.epsilon_demand_profile(0);
  require(resumed.statistics().basis_scans==0 && resumed.statistics().basis_cache_hits==1 &&
          resumed.statistics().basis_selections==1 && resumed.statistics().demand_preflights>=2,
          "resumed selected basis was rescanned or acceptance gates skipped");
  require(resumed.basis_selections().at(0).cache_hit && resumed.basis_selections().at(0).accepted,
          "persisted selection reporting");
  auto changed_resources=options;changed_resources.endpoint_order=6;changed_resources.basis_scan_seconds=0.000001;
  recursion::Evaluator refined(graph,changed_resources);(void)refined.epsilon_demand_profile(0);
  require(refined.statistics().basis_scans==0 && refined.statistics().basis_cache_hits==1 &&
          refined.statistics().demand_preflights>=2,"changed receiving resources reused demand scores or rescanned");
  auto explicit_options=options;
  explicit_options.physical_bases.emplace(0,recursion::PhysicalBasisTransform{node.closure.matrix,unit,unit,node.closure.matrix});
  recursion::Evaluator explicit_solver(graph,explicit_options);(void)explicit_solver.epsilon_demand_profile(0);
  require(explicit_solver.statistics().basis_scans==0 && explicit_solver.statistics().basis_cache_hits==0 &&
          explicit_solver.basis_selections().empty(),"explicit basis did not override cache");
  auto disabled=options;disabled.automatic_basis_scan=false;
  recursion::Evaluator disabled_solver(graph,disabled);(void)disabled_solver.epsilon_demand_profile(0);
  require(disabled_solver.statistics().basis_cache_hits==0 && disabled_solver.basis_selections().empty(),"disabled scan loaded automatic basis");
  // A completed no-improvement decision is also durable, avoiding a future
  // wall-clock scan selecting different coordinates during a continuation.
  changed=graph;changed.nodes[0].closure.matrix={{zero}};
  auto original_options=options;original_options.basis_cache_directory=root/"original";
  recursion::Evaluator first(changed,original_options);(void)first.epsilon_demand_profile(0);
  require(first.statistics().basis_scans==1 && first.statistics().basis_cache_writes==1,"original decision was not saved");
  recursion::Evaluator second(changed,original_options);(void)second.epsilon_demand_profile(0);
  require(second.statistics().basis_scans==0 && second.statistics().basis_cache_hits==1 &&
          !second.basis_selections().at(0).accepted,"original decision was rescanned");
  // A persisted coordinate remains pinned even when its cost ties the original.
  auto pinned=options;pinned.basis_cache_directory=root/"equal-cost-pinned";
  artifacts::Store pinned_store(pinned.basis_cache_directory);
  const auto eps=graph.dimension.variable(1);
  bp::save(pinned_store,id,{{node.closure.matrix,{{eps}},{{one/eps}},node.closure.matrix},true},node.closure.matrix,0);
  recursion::Evaluator tied(graph,pinned);(void)tied.epsilon_demand_profile(0);
  require(tied.basis_selections().front().cache_hit && tied.basis_selections().front().accepted,
    "persisted equal-cost coordinates silently reverted to the original basis");
  std::cout<<"Exact physical basis persistence, identity binding, fresh acceptance and explicit overrides passed\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
