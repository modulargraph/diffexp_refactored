// Opt-in native boundary reconstruction, not transport from a supplied boundary.
// Usage: diffexp_henn_ft component1|all CACHE FIRE_OR_DASH ANCILLARY_DIRECTORY
//                    [FIRE_SECONDS [GRAPH_SECONDS [adjoint|factored|auto [regular-anchor-trial]]]]
// Ancillary data are parsed as data; no Wolfram evaluation is used.
#include "stage_probe_io.hpp"
#include "diffexp/level_cache.hpp"
#include "diffexp/fire_modular.hpp"
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
  json::object report{{"schema","DiffExp3.HennFTAcceptance/v1"},{"status","error"},
    {"absolute_threshold","1e-20"},{"full_evaluator_end_to_end",false},
    {"reference_is_used_as_boundary",false},{"reference_coefficient_enclosures_certified",false}};
  json::array events;unsigned built=0,reused=0;
  try {
    if(argc<5 || argc>9 || (std::string(argv[1])!="component1" && std::string(argv[1])!="all"))
      throw std::invalid_argument("usage: diffexp_henn_ft component1|all CACHE FIRE_OR_DASH ANCILLARY_DIRECTORY [FIRE_SECONDS [GRAPH_SECONDS [adjoint|factored|auto]]]");
    const bool all=std::string(argv[1])=="all",cache_only=std::string(argv[3])=="-";
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
    constexpr int top=4;
    const int raw_high=std::max(0,henn::needed_scalar_high(canonical,top));
    report["requested_canonical_high"]=top;report["requested_scalar_high"]=raw_high;
    recursion::Options exact;recursion::NumericalOptions numerical;AdjointConditioningStats conditioning;
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
    if(!cache_only)exact.reduction.provider.executable=argv[3];
    exact.reduction.batch_cache_directory=cache/"fire-batches";
    exact.reduction.pending_cache_directory=cache/"fire-pending";
    numerical.endpoint_cache_directory=cache/"affine-series";numerical.ordinary_cache_directory=cache/"ordinary-transport";numerical.adjoint.conditioning_stats=&conditioning;
    numerical.endpoint.univariate_epsilon_recurrence=true;
    // The 69-column N32 endpoint exceeds the generic two-billion-product
    // receiving cap. This finite resource budget keeps every exact check.
    numerical.endpoint_cache_verification.max_term_products=20000000000ULL;
    report["resources"]=json::object{{"endpoint_order",numerical.endpoint_order},{"ordinary_order",numerical.ordinary_order},
      {"working_bits",numerical.working_bits},{"leaf_digits",numerical.leaf_digits},
      {"endpoint_coefficient_domain","exact univariate rational epsilon"},{"endpoint_cache_layout","incremental columns v1; legacy read compatibility"},
      {"endpoint_verification_product_budget",numerical.endpoint_cache_verification.max_term_products},
      {"physical_endpoint_constraint_tolerance",numerical.endpoint_constraint_tolerance.str()},
      {"endpoint_projection_product_budget",numerical.endpoint.max_projection_products},
      {"endpoint_limit_projection","fixed nonpositive-power sectors; full N32 series retained for integration"},
      {"fire_seconds",exact.reduction.provider.timeout_seconds},{"graph_seconds",exact.total_timeout_seconds},
      {"fire_workers",exact.reduction.provider.threads},{"simplifier_workers",exact.reduction.provider.simplifier_threads},
      {"fire_memory_bytes",exact.reduction.provider.memory_bytes}};
    numerical.progress=[&](std::size_t depth,const std::string& text,int high) {
      std::cerr<<elapsed(start)<<"s numerical level "<<depth<<": "<<text<<", epsilon through "<<high<<'\n';
    };
    artifacts::Store store(cache);
    double modular_seconds=0;std::size_t modular_batches=0;
    report["uncached_reduction_provider"]=cache_only?"disabled":"native finite-field reduction and rational reconstruction";
    auto provider=[&](const auto& b,const auto& d,const auto& f,auto p,const auto& sources,const auto& options) {
      level::Provider fallback;
      if(cache_only)fallback=[](const auto&,const auto&){fire::Result result;result.reason="Henn acceptance cache-only mode forbids FIRE";return result;};
      else {
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
    json::array anchors;for(const auto& node:graph.nodes){anchors.emplace_back(node.anchor.str());std::cerr<<"matching anchor "<<anchors.size()<<": "<<node.anchor.str()<<'\n';}
    report["matching_anchors"]=std::move(anchors);
    recursion::Evaluator evaluator(graph,numerical);auto raw=evaluator.evaluate(raw_high);
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
      {"leaf_evaluations",s.leaf_evaluations},{"refinements",s.refinements},{"endpoint_series_built",s.endpoint_series_built},
      {"endpoint_series_reused",s.endpoint_series_reused},{"demand_preflights",s.demand_preflights},{"local_operator_reuses",s.local_operator_reuses},
      {"adjoint_selections",s.adjoint_selections},{"factored_selections",s.factored_selections},
      {"conditioning_method_fallbacks",s.conditioning_method_fallbacks},
      {"endpoint_constraint_rows",s.endpoint_constraint_rows},{"endpoint_constraint_coefficients",s.endpoint_constraint_coefficients},
      {"maximum_endpoint_constraint_residual",s.maximum_endpoint_constraint_residual}};
    report["ordinary_checkpoints"]=json::object{{"loaded",s.ordinary_checkpoints.loaded},{"saved",s.ordinary_checkpoints.saved},{"completed_reused",s.ordinary_checkpoints.completed_reused},{"omitted_tail_certified",false}};
    report["conditioning"]=json::object{{"centered_charts",conditioning.centered_charts},
      {"conditioning_subdivisions",conditioning.conditioning_subdivisions},{"rational_cross_checks",conditioning.rational_cross_checks}};
    report["status"]=pass?"pass":"fail";
  } catch(const std::exception& error) {
    report["phase"]=phase;report["error"]=error.what();std::cerr<<phase<<": "<<error.what()<<'\n';
  }
  report["exact_events"]=std::move(events);report["exact_systems_built"]=built;report["exact_systems_reused"]=reused;
  report["total_seconds"]=elapsed(start);output<<json::serialize(report)<<'\n';
  return report.at("status")=="pass"?0:1;
}
