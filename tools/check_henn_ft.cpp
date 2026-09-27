// Opt-in native boundary reconstruction, not transport from a supplied boundary.
// Usage: diffexp_henn_ft component1|all CACHE FIRE_OR_DASH ANCILLARY_DIRECTORY
//                    [FIRE_SECONDS [GRAPH_SECONDS [adjoint|factored|auto [regular-anchor-trial]]]]
// Ancillary data are parsed as data; no Wolfram evaluation is used.
#include "stage_probe_io.hpp"
#include "diffexp/level_cache.hpp"
#include "diffexp/fire_modular.hpp"
#include "diffexp/ibp_solver_provider.hpp"
#include "diffexp/henn_boundary.hpp"
#include "diffexp/canonical.hpp"
#include <charconv>
#include <iostream>

using namespace diffexp;
namespace json=boost::json;
using B=kernel::ComplexBall;
using Clock=std::chrono::steady_clock;
static double elapsed(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now()-start).count();
}
static unsigned budget(const char* text) {
  std::string s(text);unsigned result=0;
  auto parsed=std::from_chars(s.data(),s.data()+s.size(),result);
  if(parsed.ec!=std::errc{} || parsed.ptr!=s.data()+s.size() || !result || result>86400)
    throw std::invalid_argument("time budgets must be integers in 1..86400");
  return result;
}
static json::object ball(const B& value) {
  auto part=[](const arb_t x){char* text=arb_get_str(x,45,0);std::string result(text);flint_free(text);return result;};
  return {{"real",part(acb_realref(value.raw()))},{"imaginary",part(acb_imagref(value.raw()))}};
}
int main(int argc,char** argv) {
  std::ostream output(std::cout.rdbuf());
  struct Redirect {
    std::streambuf* old=std::cout.rdbuf(std::cerr.rdbuf());
    ~Redirect(){std::cout.rdbuf(old);}
  } redirect;
  auto start=Clock::now();std::string phase="arguments";
  json::object report{{"schema","DiffExp.HennFTAcceptance/v1"},{"status","error"},
    {"absolute_threshold","1e-20"},{"full_evaluator_end_to_end",false},
    {"reference_is_used_as_boundary",false},{"reference_coefficient_enclosures_certified",false}};
  json::array events;unsigned built=0,reused=0;AdjointConditioningStats conditioning;
  try {
    if(argc<5 || argc>10 || (std::string(argv[1])!="component1" && std::string(argv[1])!="all"))
      throw std::invalid_argument("usage: diffexp_henn_ft component1|all CACHE FIRE_OR_IBP_SOLVER_OR_DASH ANCILLARY_DIRECTORY [PROBE_SECONDS [GRAPH_SECONDS [adjoint|factored|auto [regular-anchor-trial [SETTINGS_JSON]]]]]");
    const bool all=std::string(argv[1])=="all",cache_only=std::string(argv[3])=="-";
    const bool native_ibp=std::string(argv[3])=="ibp-solver";
    json::object settings;if(argc>9)settings=json::parse(fire::read_text(argv[9])).as_object();
    const std::set<std::string> allowed{"compact_centered_only","circuit_grouped_dot","compact_centered_recovery","rational_circuit_recurrence","endpoint_order","ordinary_order","centered_map_working_bits","centered_before_rational","working_bits","leaf_digits","canonical_high","specialize_exact_point","ordinary_method","prepare_only","ibp_numerators","ibp_dots","ibp_max_degree","dense_endpoint_verification","univariate_endpoint_projection","cleared_endpoint_projection","finite_lag_endpoint_projection","preflight_depth","evaluate_depth","evaluate_high","direct_endpoint_functionals","spectral_endpoint_clustering","spectral_max_block_size","spectral_max_nodes","spectral_max_block_nodes","spectral_seconds_budget","spectral_max_subdivisions","transport_capture_directory","transport_capture_stop_after"};
    for(const auto& item:settings)if(!allowed.contains(std::string(item.key())) &&
        item.key()!="basis_scan_artifacts" && item.key()!="epsilon_demand_depths" && item.key()!="centered_map_workers" &&
        item.key()!="endpoint_power_scan_directory" && item.key()!="automatic_basis_scan" &&
        item.key()!="basis_scan_seconds" && item.key()!="ordinary_cache_directory" && item.key()!="component_epsilon_demands" &&
        item.key()!="component_stage_directory" &&
        item.key()!="certify_operator_principal_part" && item.key()!="principal_part_endpoint_milliseconds")
      throw std::invalid_argument("unknown acceptance setting: "+std::string(item.key()));
    auto number=[&](const char* key,unsigned fallback,unsigned low,unsigned high){auto value=settings.contains(key)?settings.at(key).to_number<unsigned>():fallback;if(value<low||value>high)throw std::invalid_argument(std::string("invalid acceptance setting: ")+key);return value;};
    const bool prepare_only=settings.contains("prepare_only")&&settings.at("prepare_only").as_bool();
    if(settings.contains("evaluate_high")&&!settings.contains("evaluate_depth"))
      throw std::invalid_argument("evaluate_high requires evaluate_depth");
    if(settings.contains("evaluate_depth")&&(prepare_only||settings.contains("preflight_depth")))
      throw std::invalid_argument("subtree evaluation cannot be combined with preparation-only modes");
    report["settings"]=settings;
    const std::string method=argc>7?argv[7]:"adjoint";
    const bool regular_anchor_trial=argc>8&&std::string(argv[8])=="regular-anchor-trial";
    if(argc>8&&!regular_anchor_trial)throw std::invalid_argument("unknown Henn acceptance contour trial");
    if(method!="adjoint" && method!="factored" && method!="auto")
      throw std::invalid_argument("Henn linear method must be adjoint, factored or auto");
    report["linear_transport"]=method;
    const std::filesystem::path cache=argv[2],data_directory=argv[4];
    const auto definitions=data_directory/"dlogBasisXB.txt",reference_file=data_directory/"XB_Boundary_values_X0.txt";
    report["scope"]=all?"all108":"component1";report["cache_only"]=cache_only;
    report["basis_path"]=std::filesystem::absolute(definitions).string();
    report["basis_sha256"]=artifacts::detail::sha256(fire::read_text(definitions,2000000));
    report["reference_path"]=std::filesystem::absolute(reference_file).string();
    report["reference_sha256"]=artifacts::detail::sha256(fire::read_text(reference_file,2000000));
    report["reference_kind"]="published independent canonical boundary decimal table; rounding balls are not reference error certificates";
    phase="canonical definitions";auto canonical=henn::read_x0(definitions.string());
    if(!all) {
      // This explicitly selected observable has its own demand identity. Keep
      // the full 257-target union unchanged in the default all-component mode.
      canonical.components.resize(1);canonical.scalar_targets.clear();
      for(const auto& [index,c]:canonical.components.front())if(!c.is_zero())canonical.scalar_targets.push_back(index);
      canonical.nonzero_targets=canonical.scalar_targets;
      if(canonical.scalar_targets.empty())throw std::invalid_argument("component 1 unexpectedly has no scalar source");
    }
    report["canonical_rows"]=canonical.components.size();report["scalar_target_count"]=canonical.scalar_targets.size();
    json::array targets;
    for(const auto& index:canonical.scalar_targets){json::array row;for(auto n:index)row.emplace_back(n);targets.push_back(std::move(row));}
    report["ordered_scalar_targets"]=std::move(targets);
    const int top=number("canonical_high",4,0,4);
    const int raw_high=std::max(0,henn::needed_scalar_high(canonical,top));
    report["requested_canonical_high"]=top;report["requested_scalar_high"]=raw_high;
    recursion::Options exact;recursion::NumericalOptions numerical;
    numerical.automatic_basis_scan=!settings.contains("automatic_basis_scan") || settings.at("automatic_basis_scan").as_bool();
    numerical.basis_scan_seconds=number("basis_scan_seconds",30,1,60);
    numerical.component_epsilon_demands=!settings.contains("component_epsilon_demands") || settings.at("component_epsilon_demands").as_bool();
    report["basis_selections"]=json::array{};
    numerical.basis_selection_observer=[&](const recursion::BasisSelection& selected) {
      json::object item{{"depth",selected.depth},{"masters",selected.masters},{"accepted",selected.accepted},
        {"original_expression_characters",selected.original_characters},{"candidate_expression_characters",selected.candidate_characters},
        {"candidates_checked",selected.candidates_checked},{"original_child_loss",selected.original_child_loss},{"candidate_child_loss",selected.candidate_child_loss},
        {"candidates_screened",selected.candidates_screened},{"full_candidates_checked",selected.full_candidates_checked},
        {"scan_seconds",selected.scan_seconds},{"cache_hit",selected.cache_hit},{"total_selection_seconds",selected.total_selection_seconds},{"reason",selected.reason},{"transport_trials",0}};
      report["basis_selections"].as_array().push_back(item);
      std::cerr<<"HENN_BASIS_SELECTION "<<json::serialize(item)<<'\n';
    };
    std::size_t active_depth=0;
    numerical.ordinary_fraction_progress=[&](double fraction) {
      std::cerr<<"HENN_PROGRESS "<<json::serialize(json::object{{"kind","arm"},{"depth",active_depth},
        {"fraction",fraction},{"seconds",elapsed(start)}})<<'\n';
    };
    numerical.endpoint_order=number("endpoint_order",32,4,128);
    numerical.direct_endpoint_functionals=settings.contains("direct_endpoint_functionals")&&settings.at("direct_endpoint_functionals").as_bool();
    numerical.spectral.endpoint_clustering=settings.contains("spectral_endpoint_clustering")&&settings.at("spectral_endpoint_clustering").as_bool();
    numerical.spectral.max_block_size=number("spectral_max_block_size",4,1,32);
    numerical.spectral.max_nodes=number("spectral_max_nodes",128,16,128);
    numerical.spectral.max_block_nodes=number("spectral_max_block_nodes",256,16,2048);
    numerical.spectral.seconds_budget=number("spectral_seconds_budget",15,1,1200);
    numerical.spectral.max_subdivisions=number("spectral_max_subdivisions",0,0,256);
    numerical.ordinary_progress=[&](const auto& text){std::cerr<<elapsed(start)<<"s "<<text<<'\n';};
    numerical.ordinary_order=number("ordinary_order",80,8,1000);
    numerical.adjoint.rational_circuit_recurrence=!settings.contains("rational_circuit_recurrence") || settings.at("rational_circuit_recurrence").as_bool();
    numerical.adjoint.compact_centered_only=!settings.contains("compact_centered_only") || settings.at("compact_centered_only").as_bool();
    numerical.adjoint.circuit_grouped_dot=!settings.contains("circuit_grouped_dot") || settings.at("circuit_grouped_dot").as_bool();
    numerical.adjoint.compact_centered_recovery=!settings.contains("compact_centered_recovery") || settings.at("compact_centered_recovery").as_bool();
    // Captures must record effective execution policy even when defaults were used.
    settings["rational_circuit_recurrence"]=numerical.adjoint.rational_circuit_recurrence;
    settings["compact_centered_recovery"]=numerical.adjoint.compact_centered_recovery;
    settings["compact_centered_only"]=numerical.adjoint.compact_centered_only;
    settings["circuit_grouped_dot"]=numerical.adjoint.circuit_grouped_dot;
    report["settings"]=settings;
    numerical.adjoint.centered_before_rational=settings.contains("centered_before_rational")&&settings.at("centered_before_rational").as_bool();
    numerical.adjoint.centered_map_working_bits=number("centered_map_working_bits",0,0,4096);
    numerical.adjoint.centered_map_workers=number("centered_map_workers",1,1,4);
    if(numerical.adjoint.centered_map_working_bits && numerical.adjoint.centered_map_working_bits<64)
      throw std::invalid_argument("centered_map_working_bits must be zero or at least64");
    numerical.working_bits=number("working_bits",384,128,4096);
    std::size_t captured_problems=0;
    if(settings.contains("transport_capture_directory")) {
      auto directory=std::filesystem::path(std::string(settings.at("transport_capture_directory").as_string()));
      std::filesystem::create_directories(directory);
      numerical.ordinary_problem_observer=[&,directory](const auto& matrix,const auto& initial,const auto& forcing,const auto& vertices) {
        json::array variables,path;for(const auto& name:vertices.front().variables())variables.emplace_back(name);
        for(const auto& point:vertices)path.emplace_back(point.str());
        json::object problem{{"schema","DiffExp.FTTransportProblem/v1"},{"variables",variables},{"matrix",numerical_rows_io::exact_matrix(matrix)},
          {"initial",numerical_rows_io::exact_rows(initial)},{"forcing",numerical_rows_io::exact_matrix(forcing)},
          {"path",path},{"working_bits",B::precision()},{"settings",settings}};
        auto file=directory/("problem-"+std::to_string(++captured_problems)+".json");
        if(std::filesystem::exists(file))throw std::runtime_error("transport capture refuses to overwrite a problem");
        std::ofstream stream(file);stream<<json::serialize(problem)<<'\n';stream.close();
        if(!stream)throw std::runtime_error("transport capture write failed");
        if(captured_problems>=number("transport_capture_stop_after",1000,1,1000))throw std::runtime_error("requested transport capture completed; numerical acceptance not run");
      };
    }
    numerical.leaf_digits=number("leaf_digits",28,20,500);
    numerical.matching.specialize_exact_point=settings.contains("specialize_exact_point")&&settings.at("specialize_exact_point").as_bool();
    numerical.matching.certify_operator_principal_part=settings.contains("certify_operator_principal_part")&&settings.at("certify_operator_principal_part").as_bool();
    numerical.principal_part_endpoint_milliseconds=number("principal_part_endpoint_milliseconds",5000,0,60000);
    if(settings.contains("ordinary_method")){
      auto ordinary=std::string(settings.at("ordinary_method").as_string());
      if(ordinary!="auto"&&ordinary!="spectral"&&ordinary!="taylor"&&ordinary!="ultraspherical")throw std::invalid_argument("ordinary method");
      numerical.ordinary_method=ordinary=="ultraspherical"?recursion::OrdinaryMethod::ultraspherical:
        ordinary=="spectral"?recursion::OrdinaryMethod::spectral:ordinary=="taylor"?recursion::OrdinaryMethod::taylor:recursion::OrdinaryMethod::automatic;
    }
    exact.anchors=causal::henn_anchors();
    numerical.causal_prescription=causal::current_example("henn_double_pentagon_x0",7,exact.anchors,true);
    if(regular_anchor_trial) {
      exact.anchors=causal::henn_anchors();exact.anchors.front()=Rational("1/3");
      auto prescription=causal::current_example("henn_double_pentagon_x0",7,causal::henn_anchors(),true);
      prescription.provenance+="; explicit acceptance trial: first matching anchor1/3, unchanged lower F rim and level signs. Requires independent canonical reference comparison; no homotopy certificate claimed.";
      numerical.causal_prescription=std::move(prescription);
    }
    report["regular_anchor_contour_trial"]=regular_anchor_trial;
    numerical.linear_method=method=="auto"?recursion::LinearMethod::automatic:
      method=="factored"?recursion::LinearMethod::factored:recursion::LinearMethod::adjoint;
    exact.reduction.provider.timeout_seconds=argc>5?budget(argv[5]):600;
    exact.total_timeout_seconds=argc>6?budget(argv[6]):1800;
    exact.reduction.total_timeout_seconds=exact.total_timeout_seconds;
    exact.reduction.provider.threads=4;exact.reduction.provider.simplifier_threads=1;
    exact.reduction.provider.memory_bytes=6ull*1024*1024*1024;
    if(!cache_only&&!native_ibp)exact.reduction.provider.executable=argv[3];
    exact.reduction.batch_cache_directory=cache/"fire-batches";
    exact.reduction.pending_cache_directory=cache/"fire-pending";
    numerical.endpoint_cache_directory=cache/"affine-series";numerical.ordinary_cache_directory=cache/"ordinary-transport";numerical.adjoint.conditioning_stats=&conditioning;
    if(settings.contains("ordinary_cache_directory"))
      numerical.ordinary_cache_directory=artifacts::detail::string(settings.at("ordinary_cache_directory"));
    numerical.endpoint.univariate_epsilon_recurrence=true;
    numerical.endpoint.finite_lag_projection=settings.contains("finite_lag_endpoint_projection")&&settings.at("finite_lag_endpoint_projection").as_bool();
    numerical.endpoint.cleared_epsilon_projection=settings.contains("cleared_endpoint_projection")&&settings.at("cleared_endpoint_projection").as_bool();
    numerical.endpoint.univariate_epsilon_projection=settings.contains("univariate_endpoint_projection")&&settings.at("univariate_endpoint_projection").as_bool();
    // The 69-column N32 endpoint exceeds the generic two-billion-product
    // receiving cap. This finite resource budget keeps every exact check.
    numerical.endpoint_cache_verification.dense_epsilon_polynomials=settings.contains("dense_endpoint_verification")&&settings.at("dense_endpoint_verification").as_bool();
    numerical.endpoint_cache_verification.max_term_products=20000000000ULL;
    report["resources"]=json::object{{"endpoint_order",numerical.endpoint_order},{"ordinary_order",numerical.ordinary_order},
      {"working_bits",numerical.working_bits},{"leaf_digits",numerical.leaf_digits},
      {"endpoint_coefficient_domain","exact univariate rational epsilon"},{"endpoint_cache_layout","incremental columns v1; legacy read compatibility"},
      {"endpoint_verification_product_budget",numerical.endpoint_cache_verification.max_term_products},
      {"physical_endpoint_constraint_tolerance",numerical.endpoint_constraint_tolerance.str()},
      {"endpoint_projection_product_budget",numerical.endpoint.max_projection_products},
      {"endpoint_limit_projection","fixed nonpositive-power sectors; full configured-order series retained for integration"},
      {"fire_seconds",exact.reduction.provider.timeout_seconds},{"graph_seconds",exact.total_timeout_seconds},
      {"fire_workers",exact.reduction.provider.threads},{"simplifier_workers",exact.reduction.provider.simplifier_threads},
      {"fire_memory_bytes",exact.reduction.provider.memory_bytes}};
    numerical.progress=[&](std::size_t depth,const std::string& text,int high) {
      active_depth=depth;
      std::cerr<<"HENN_PROGRESS "<<json::serialize(json::object{{"kind","level"},{"depth",depth},
        {"phase",text},{"epsilon_high",high},{"seconds",elapsed(start)}})<<'\n';
      std::cerr<<elapsed(start)<<"s numerical level "<<depth<<": "<<text<<", epsilon through "<<high<<'\n';
    };
    artifacts::Store store(cache);
    double modular_seconds=0;std::size_t modular_batches=0;
    report["uncached_reduction_provider"]=cache_only?"disabled":native_ibp?"IBP Solver":"FIRE finite-field reduction and rational reconstruction";
    auto provider=[&](const auto& b,const auto& d,const auto& f,auto p,const auto& sources,const auto& options) {
      level::Provider fallback;
      if(cache_only)fallback=[](const auto&,const auto&){fire::Result result;result.reason="Henn acceptance cache-only mode forbids FIRE";return result;};
      else if(native_ibp){
        ibp_solver::Options native;native.cache_directory=cache/"ibp-solver";native.degree_slices=true;native.probe_timeout_seconds=options.provider.timeout_seconds;native.sample_checkpoint_interval=32;
        native.closure_parameter=p;native.selected_source_rows=true;native.target_sector_seeds=true;
        native.numerators=number("ibp_numerators",1,0,8);native.dots=number("ibp_dots",1,0,8);native.max_degree=number("ibp_max_degree",64,1,64);
        native.progress=[&](const auto& text){std::cerr<<elapsed(start)<<"s "<<text<<'\n';};
        auto session=std::make_shared<ibp_solver::Session>(b,d,f,native);
        fallback=[session,&modular_seconds,&modular_batches](const auto& batch,const auto& limits){const auto began=Clock::now();auto result=(*session)(batch,limits);modular_seconds+=elapsed(began);++modular_batches;return result;};
      }else {
        fire_modular::Options modular;modular.executable=std::filesystem::path(argv[3]).parent_path()/"FIRE7p";
        modular.cache_directory=cache/"fire-modular";
        modular.progress=[&](const std::string& text){std::cerr<<elapsed(start)<<"s "<<text<<'\n';};
        auto session=std::make_shared<fire_modular::Session>(b,d,f,modular,options.batch_cache_directory);
        fallback=[session,&modular_seconds,&modular_batches](const auto& batch,const auto& limits){
          const auto began=Clock::now();auto result=(*session)(batch,limits);
          modular_seconds+=elapsed(began);++modular_batches;return result;
        };
      }
      auto prepared=cached_level::prepare(store,b,d,f,p,sources,options,fallback);
      if(prepared.cache_hit)++reused;else if(prepared.result.success)++built;
      return std::move(prepared.result);
    };
    phase="exact graph";auto began=Clock::now();
    auto graph=recursion::prepare(feynman::example_family("henn_double_pentagon_x0"),canonical.scalar_targets,exact,provider,
      [&](const recursion::Event& event) {
        std::cerr<<"exact level "<<event.depth<<": "<<event.masters<<" masters, "<<event.elapsed_seconds<<" s\n";
        events.push_back(json::object{{"depth",event.depth},{"masters",event.masters},{"requests",event.requests},
          {"sources",event.sources},{"provider_passes",event.provider_passes},{"seconds",event.elapsed_seconds},{"scalar_leaf",event.scalar_leaf}});
      });
    report["exact_seconds"]=elapsed(began);report["modular_provider_seconds"]=modular_seconds;
    report["modular_provider_batches"]=modular_batches;phase="native FT evaluation";began=Clock::now();
    if(prepare_only){report["status"]="prepared";report["exact_events"]=events;report["exact_systems_built"]=built;report["exact_systems_reused"]=reused;report["total_seconds"]=elapsed(start);output<<json::serialize(report)<<'\n';return 0;}
    json::array anchors;for(const auto& node:graph.nodes){anchors.emplace_back(node.anchor.str());std::cerr<<"matching anchor "<<anchors.size()<<": "<<node.anchor.str()<<'\n';}
    report["matching_anchors"]=std::move(anchors);
    if(settings.contains("basis_scan_artifacts")) {
      json::array loaded;
      for(const auto& selected:settings.at("basis_scan_artifacts").as_array()) {
        const auto& selection=selected.as_object();artifacts::detail::keys(selection,{"depth","file"});
        const auto depth=selection.at("depth").to_number<unsigned>();
        if(depth>=graph.nodes.size() || graph.nodes[depth].scalar_leaf)
          throw std::invalid_argument("physical basis selection requires a nonleaf depth");
        const auto file=artifacts::detail::string(selection.at("file"));
        const auto contents=fire::read_text(file,64*1024*1024);
        const auto artifact=json::parse(contents).as_object();
        if(artifact.at("schema")!="DiffExp.FTEpsilonBasisScan/v1")
          throw std::invalid_argument("physical basis scan schema");
        std::vector<std::string> names;
        for(const auto& name:artifact.at("variables").as_array())names.push_back(artifacts::detail::string(name));
        const auto& graph_names=graph.dimension.variables();
        if(names.size()>graph_names.size() || !std::equal(names.begin(),names.end(),graph_names.begin()))
          throw std::invalid_argument("physical basis symbol ordering mismatch");
        const auto& original=graph.nodes[depth].closure.ordered_basis;
        std::vector<diffexp::ibp::Integral> ids;
        for(const auto& row:artifact.at("original_ordered_basis").as_array()) {
          ids.emplace_back();for(const auto& power:row.as_array())ids.back().push_back(power.to_number<int>());
        }
        if(ids!=original)throw std::invalid_argument("physical basis original integral ordering mismatch");
        const auto matrix=[&](const char* key) {
          const auto& rows=artifact.at(key).as_array();
          if(rows.size()!=original.size())throw std::invalid_argument("physical basis row count");
          ExactEpsilonMatrix result;
          for(const auto& row:rows) {
            if(row.as_array().size()!=original.size())throw std::invalid_argument("physical basis column count");
            result.emplace_back();for(const auto& entry:row.as_array())
              result.back().push_back(graph.dimension.parse(artifacts::detail::string(entry)));
          }
          return result;
        };
        recursion::PhysicalBasisTransform basis{matrix("original_connection"),matrix("S"),matrix("S_inverse"),matrix("new_connection")};
        // Recheck the algebra against the actual immutable graph; producer
        // flags and a sample-level complexity score are not acceptance proof.
        recursion::verify_physical_basis(basis,graph.nodes[depth].closure.matrix,path_epsilon_variables(graph.dimension).first);
        if(!numerical.physical_bases.emplace(depth,std::move(basis)).second)
          throw std::invalid_argument("duplicate physical basis depth");
        loaded.push_back(json::object{{"depth",depth},{"file",file},
          {"sha256",artifacts::detail::sha256(contents)},{"exact_identities_verified",true}});
      }
      report["physical_bases"]=std::move(loaded);
    }
    recursion::Evaluator evaluator(graph,numerical);
    if(settings.contains("epsilon_demand_depths")) {
      if(settings.contains("preflight_depth") || settings.contains("evaluate_depth"))
        throw std::invalid_argument("epsilon demand inspection cannot be combined with another diagnostic mode");
      phase="epsilon demand inspection";json::array profiles;std::set<unsigned> seen;
      int new_total=raw_high,old_total=raw_high;
      for(const auto& item:settings.at("epsilon_demand_depths").as_array()) {
        const auto depth=item.to_number<unsigned>();
        if(!seen.insert(depth).second)throw std::invalid_argument("duplicate epsilon inspection depth");
        const auto profile_start=Clock::now();auto p=evaluator.epsilon_demand_profile(depth);
        json::object profile{{"depth",p.depth},{"level",p.depth+1},{"masters",p.masters},
          {"lower_endpoint_low",p.lower_endpoint_low},{"upper_endpoint_low",p.upper_endpoint_low},
          {"beta_low",p.beta_low},{"anchor_inverse_pole_loss",p.anchor_inverse_pole_loss},
          {"adjoint_child_loss",p.adjoint_child_loss},{"factored_child_loss",p.factored_child_loss},
          {"old_common_high_child_loss",p.common_high_child_loss},{"seconds",elapsed(profile_start)}};
        json::array component_losses;for(const auto& loss:p.child_component_loss)
          if(loss)component_losses.emplace_back(*loss);else component_losses.emplace_back(nullptr);
        profile["child_component_loss"]=std::move(component_losses);
        new_total+=p.factored_child_loss;old_total+=p.common_high_child_loss;
        std::cerr<<"HENN_EPSILON_DEMAND "<<json::serialize(profile)<<'\n';profiles.push_back(std::move(profile));
        if(settings.contains("component_stage_directory")) {
          const auto directory=std::filesystem::path(artifacts::detail::string(settings.at("component_stage_directory")));
          std::filesystem::create_directories(directory);
          const auto destination=directory/("level-"+std::to_string(depth+1)+".json");
          if(std::filesystem::exists(destination))throw std::runtime_error("component fixture refuses to overwrite an artifact");
          auto stage=evaluator.prepare_factored_consumers(depth,20);
          json::array names;for(const auto& name:graph.dimension.variables())names.emplace_back(name);
          json::object fixture{{"schema","DiffExp.ComponentChartConsumers/v1"},{"variables",names},
            {"A",stage_probe::exact_matrix(stage.connection)},{"beta",stage_probe::exact_matrix(stage.beta)},
            {"right",numerical_rows_io::exact_rows(stage.right)},{"output_high",10},{"depth",depth},
            {"settings",settings},{"endpoint_map_kind","retained numerical operator; only exact zero coefficients used for demand bounds"},
            {"basis_sha256",report.at("basis_sha256")},{"omitted_tail_certified",false}};
          std::ofstream stream(destination);stream<<json::serialize(fixture)<<'\n';stream.close();
          if(!stream)throw std::runtime_error("component fixture write failed");
        }
        if(settings.contains("endpoint_power_scan_directory")) {
          auto proposal=evaluator.optimize_endpoint_epsilon_basis(depth);
          json::array ids,names,shifts;
          for(const auto& integral:graph.nodes[depth].closure.ordered_basis) {
            json::array row;for(auto power:integral)row.emplace_back(power);ids.push_back(std::move(row));
          }
          for(const auto& name:graph.dimension.variables())names.emplace_back(name);
          for(auto shift:proposal.additional_shifts)shifts.emplace_back(shift);
          json::object artifact{{"schema","DiffExp.FTEpsilonBasisScan/v1"},{"variables",names},
            {"original_ordered_basis",ids},{"original_connection",stage_probe::exact_matrix(proposal.basis.original_connection)},
            {"S",stage_probe::exact_matrix(proposal.basis.to_basis)},{"S_inverse",stage_probe::exact_matrix(proposal.basis.from_basis)},
            {"new_connection",stage_probe::exact_matrix(proposal.basis.connection)},
            {"additional_epsilon_shifts",shifts},{"original_child_loss",proposal.original_child_loss},
            {"estimated_child_loss",proposal.estimated_child_loss},{"exact_identities_verified",true},
            {"method","integer difference constraints using endpoint, integrand and anchor maps"},
            {"production_adoption_authorized_by_this_artifact",false}};
          const auto directory=std::filesystem::path(artifacts::detail::string(settings.at("endpoint_power_scan_directory")));
          std::filesystem::create_directories(directory);const auto destination=directory/("level-"+std::to_string(depth+1)+".json");
          if(std::filesystem::exists(destination))throw std::runtime_error("endpoint epsilon scan refuses to overwrite an artifact");
          std::ofstream stream(destination);stream<<json::serialize(artifact)<<'\n';stream.close();
          if(!stream)throw std::runtime_error("endpoint epsilon scan artifact write failed");
          std::cerr<<"HENN_ENDPOINT_POWER_SCAN "<<json::serialize(json::object{{"depth",depth},
            {"original_child_loss",proposal.original_child_loss},{"estimated_child_loss",proposal.estimated_child_loss},
            {"additional_epsilon_shifts",shifts},{"file",destination.string()}})<<'\n';
        }
      }
      report["epsilon_demand_profiles"]=std::move(profiles);
      report["complete_chain_profile"]=seen.size()+1==graph.nodes.size();
      if(seen.size()+1==graph.nodes.size()) {
        report["predicted_factored_leaf_high"]=new_total;
        report["previous_common_high_leaf_high"]=old_total;
      }
      report["status"]="epsilon_demand_inspected";report["exact_events"]=events;
      const auto& stats=evaluator.statistics();
      report["principal_part_certificates"]=json::object{{"attempts",stats.principal_part_attempts},
        {"certified",stats.principal_part_certificates},{"reused",stats.principal_part_reused},
        {"zero_coefficients",stats.principal_part_zero_coefficients},{"full_formal_tail_certificate",false}};
      report["total_seconds"]=elapsed(start);output<<json::serialize(report)<<'\n';return 0;
    }
    if(settings.contains("preflight_depth")){auto depth=number("preflight_depth",0,0,5);phase="endpoint operator preflight";auto stage=evaluator.prepare_adjoint_stage(depth,0);report["status"]="preflight_completed";report["preflight_depth"]=depth;report["full_evaluator_end_to_end"]=false;report["total_seconds"]=elapsed(start);output<<json::serialize(report)<<'\n';return 0;}
    if(settings.contains("evaluate_depth")) {
      const auto depth=number("evaluate_depth",0,0,static_cast<unsigned>(graph.nodes.size()-1));
      const auto high=number("evaluate_high",0,0,numerical.max_epsilon);
      phase="diagnostic subtree evaluation";
      report["evaluated_depth"]=depth;report["evaluated_high"]=high;
      auto partial=evaluator.evaluate(depth,high);
      report["status"]="subtree_completed";
      report["numerical_seconds"]=elapsed(began);
      report["epsilon_low"]=partial.low;report["epsilon_high"]=partial.high();
      report["reported_ft_radii_include_omitted_tails"]=partial.taylor_tail_certified;
      json::array values;
      for(const auto& row:partial.values){json::array coefficients;for(const auto& value:row)coefficients.push_back(ball(value));values.push_back(std::move(coefficients));}
      report["subtree_coefficients"]=std::move(values);
      const auto& s=evaluator.statistics();
      report["numerical_statistics"]=json::object{{"leaf_evaluations",s.leaf_evaluations},{"direct_endpoints",s.direct_endpoints},
        {"ordinary_seconds",s.ordinary_seconds},{"spectral_seconds",s.spectral_seconds},
        {"spectral_attempts",s.spectral_attempts},{"spectral_accepted",s.spectral_accepted},
        {"conditioning_method_fallbacks",s.conditioning_method_fallbacks},
        {"completed_checkpoints_reused",s.ordinary_checkpoints.completed_reused}};
      report["conditioning"]=json::object{{"centered_charts",conditioning.centered_charts},
      {"circuit_charts",conditioning.circuit_charts},{"circuit_homogeneous_columns",conditioning.circuit_homogeneous_columns},
      {"circuit_preparations",conditioning.circuit_preparations},{"circuit_addmul_operations",conditioning.circuit_addmul_operations},
      {"circuit_scalar_operations",conditioning.circuit_scalar_operations},{"circuit_dot_products",conditioning.circuit_dot_products},
      {"centered_only_charts",conditioning.centered_only_charts},{"centered_only_fallbacks",conditioning.centered_only_fallbacks},{"exact_input_map_columns_skipped",conditioning.exact_input_map_columns_skipped},{"circuit_peak_live_cells",conditioning.circuit_peak_live_cells},
      {"compilation_seconds",conditioning.compilation_seconds},{"preparation_seconds",conditioning.preparation_seconds},
      {"source_seconds",conditioning.source_seconds},{"homogeneous_seconds",conditioning.homogeneous_seconds},
      {"reference_seconds",conditioning.reference_seconds},
      {"reduced_precision_homogeneous_maps",conditioning.reduced_precision_homogeneous_maps},
      {"full_precision_map_retries",conditioning.full_precision_map_retries},
      {"centered_before_rational_charts",conditioning.centered_before_rational_charts},
      {"centered_first_fallbacks",conditioning.centered_first_fallbacks},
        {"conditioning_subdivisions",conditioning.conditioning_subdivisions},
        {"rational_cross_checks",conditioning.rational_cross_checks}};
      report["exact_events"]=events;report["exact_systems_built"]=built;report["exact_systems_reused"]=reused;
      report["total_seconds"]=elapsed(start);output<<json::serialize(report)<<'\n';return 0;
    }
    auto raw=evaluator.evaluate(raw_high);
    report["full_evaluator_end_to_end"]=true;
    report["numerical_seconds"]=elapsed(began);
    // Preserve the common leaf and coefficient map for further diagnosis; the
    // comparison below still uses the complete public evaluator result.
    report["linear_expression"]=stage_probe::exact_expression(evaluator.linear_expression());
    report["linear_expression_sha256"]=artifacts::detail::sha256(artifacts::detail::canonical(report.at("linear_expression")));
    phase="canonical projection";auto result=henn::project_boundary(canonical,graph.nodes.front().requested,raw,top,numerical.working_bits);
    report["epsilon_low"]=result.low;report["epsilon_high"]=result.high();
    report["reported_ft_radii_include_omitted_tails"]=result.taylor_tail_certified;
    phase="independent comparison";B::set_precision(numerical.working_bits);
    auto reference=read_boundary(data::read_file(reference_file.string()),108,top,numerical.working_bits);
    const auto limit=NativeTailMagnitude::lower_abs(B::from_strings("1/100000000000000000000"));
    bool pass=true;json::array comparisons,poles;
    for(std::size_t row=0;row<result.values.size();++row)for(int k=0;k<=top;++k) {
      const B actual=k<result.low?B(0):result.values.at(row).at(k-result.low);
      auto difference=actual-reference.at(row).at(k);
      const bool good=actual.is_finite() && reference[row][k].is_finite() && NativeTailMagnitude::upper_abs(difference)<=limit;
      pass&=good;comparisons.push_back(json::object{{"component",row+1},{"epsilon_order",k},
        {"ft",ball(actual)},{"reference",ball(reference[row][k])},{"difference",ball(difference)},{"pass",good}});
    }
    for(std::size_t row=0;row<result.values.size();++row)for(int k=result.low;k<0;++k) {
      const auto& value=result.values[row][k-result.low];const bool good=value.is_finite() && NativeTailMagnitude::upper_abs(value)<=limit;
      pass&=good;poles.push_back(json::object{{"component",row+1},{"epsilon_order",k},{"coefficient",ball(value)},{"pass",good}});
    }
    report["comparisons"]=std::move(comparisons);report["forbidden_pole_audit"]=std::move(poles);
    const auto& s=evaluator.statistics();report["numerical_statistics"]=json::object{{"exact_plans",s.exact_plans},
      {"component_transports",s.component_transports},{"component_packed_slots",s.component_slots},
      {"component_rectangular_slots",s.component_rectangular_slots},
      {"functional_reduction_attempts",s.functional_reduction_attempts},
      {"functional_reduction_certificates",s.functional_reduction_certificates},
      {"basis_scans",s.basis_scans},{"basis_selections",s.basis_selections},{"basis_cache_hits",s.basis_cache_hits},{"basis_cache_writes",s.basis_cache_writes},
      {"principal_part_attempts",s.principal_part_attempts},{"principal_part_certificates",s.principal_part_certificates},
      {"principal_part_reused",s.principal_part_reused},{"principal_part_zero_coefficients",s.principal_part_zero_coefficients},
      {"leaf_evaluations",s.leaf_evaluations},{"direct_endpoints",s.direct_endpoints},{"refinements",s.refinements},{"endpoint_series_built",s.endpoint_series_built},
      {"endpoint_series_reused",s.endpoint_series_reused},{"demand_preflights",s.demand_preflights},{"local_operator_reuses",s.local_operator_reuses},
      {"adjoint_selections",s.adjoint_selections},{"factored_selections",s.factored_selections},
      {"conditioning_method_fallbacks",s.conditioning_method_fallbacks},
      {"endpoint_constraint_rows",s.endpoint_constraint_rows},{"endpoint_constraint_coefficients",s.endpoint_constraint_coefficients},
      {"maximum_endpoint_constraint_residual",s.maximum_endpoint_constraint_residual}};
    report["numerical_statistics"].as_object()["spectral_attempts"]=s.spectral_attempts;
    report["numerical_statistics"].as_object()["spectral_accepted"]=s.spectral_accepted;
    report["numerical_statistics"].as_object()["spectral_reused"]=s.spectral_reused;
    report["numerical_statistics"].as_object()["spectral_subdivisions"]=s.spectral_subdivisions;
    report["numerical_statistics"].as_object()["spectral_accepted_subsegments"]=s.spectral_accepted_subsegments;
    report["numerical_statistics"].as_object()["spectral_seconds"]=s.spectral_seconds;
    report["numerical_statistics"].as_object()["ordinary_seconds"]=s.ordinary_seconds;
    json::array rejections;for(const auto& reason:s.spectral_rejections)rejections.emplace_back(reason);report["spectral_rejections"]=std::move(rejections);
    report["ordinary_checkpoints"]=json::object{{"loaded",s.ordinary_checkpoints.loaded},{"saved",s.ordinary_checkpoints.saved},{"completed_reused",s.ordinary_checkpoints.completed_reused},{"omitted_tail_certified",false}};
    report["conditioning"]=json::object{{"centered_charts",conditioning.centered_charts},
      {"circuit_charts",conditioning.circuit_charts},{"circuit_homogeneous_columns",conditioning.circuit_homogeneous_columns},
      {"circuit_preparations",conditioning.circuit_preparations},{"circuit_addmul_operations",conditioning.circuit_addmul_operations},
      {"circuit_scalar_operations",conditioning.circuit_scalar_operations},{"circuit_dot_products",conditioning.circuit_dot_products},
      {"centered_only_charts",conditioning.centered_only_charts},{"centered_only_fallbacks",conditioning.centered_only_fallbacks},{"exact_input_map_columns_skipped",conditioning.exact_input_map_columns_skipped},{"circuit_peak_live_cells",conditioning.circuit_peak_live_cells},
      {"compilation_seconds",conditioning.compilation_seconds},{"preparation_seconds",conditioning.preparation_seconds},
      {"source_seconds",conditioning.source_seconds},{"homogeneous_seconds",conditioning.homogeneous_seconds},
      {"reference_seconds",conditioning.reference_seconds},
      {"reduced_precision_homogeneous_maps",conditioning.reduced_precision_homogeneous_maps},
      {"full_precision_map_retries",conditioning.full_precision_map_retries},
      {"centered_before_rational_charts",conditioning.centered_before_rational_charts},
      {"centered_first_fallbacks",conditioning.centered_first_fallbacks},
      {"conditioning_subdivisions",conditioning.conditioning_subdivisions},{"rational_cross_checks",conditioning.rational_cross_checks}};
    report["status"]=pass?"pass":"fail";
  } catch(const std::exception& error) {
    report["phase"]=phase;report["error"]=error.what();std::cerr<<phase<<": "<<error.what()<<'\n';
  }
  // Preserve phase costs on failures as well as completed calculations.
  report["transport_work"]=json::object{
    {"circuit_charts",conditioning.circuit_charts},{"circuit_homogeneous_columns",conditioning.circuit_homogeneous_columns},
    {"circuit_preparations",conditioning.circuit_preparations},{"circuit_addmul_operations",conditioning.circuit_addmul_operations},
    {"circuit_scalar_operations",conditioning.circuit_scalar_operations},{"circuit_dot_products",conditioning.circuit_dot_products},
      {"centered_only_charts",conditioning.centered_only_charts},{"centered_only_fallbacks",conditioning.centered_only_fallbacks},{"exact_input_map_columns_skipped",conditioning.exact_input_map_columns_skipped},{"circuit_peak_live_cells",conditioning.circuit_peak_live_cells},
    {"compilation_seconds",conditioning.compilation_seconds},{"preparation_seconds",conditioning.preparation_seconds},
    {"source_seconds",conditioning.source_seconds},{"homogeneous_seconds",conditioning.homogeneous_seconds},
    {"reference_seconds",conditioning.reference_seconds},{"rational_cross_checks",conditioning.rational_cross_checks},
    {"conditioning_subdivisions",conditioning.conditioning_subdivisions}};
  report["exact_events"]=std::move(events);report["exact_systems_built"]=built;report["exact_systems_reused"]=reused;
  report["total_seconds"]=elapsed(start);output<<json::serialize(report)<<'\n';
  return report.at("status")=="pass"?0:1;
}
