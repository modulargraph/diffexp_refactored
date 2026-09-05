// Opt-in, potentially expensive acceptance executable. No automatic family sweep.
// Usage: test_ft_examples FAMILY [CACHE_DIR [FIRE_PATH_OR_DASH [ENDPOINT_N [ORDINARY_N [BITS [GRAPH_TIMEOUT_SECONDS [METHOD]]]]]]]
// Supported: bubble, sunrise, banana, banana_unequal, banana4,
// banana4_unequal, box, pentagon, pentagon_massive, kite, box_bubble, box_triangle, double_box_planar. Defaults N32/N80/384 bits and graph timeout 600s.
// METHOD is adjoint (default), factored, auto or values.
// '-' leaves FIRE unavailable; verified durable cache hits still work.
// Output: one Boost.JSON report on stdout, phase progress on stderr.
#include "diffexp/recursion_pipeline.hpp"
#include "diffexp/level_cache.hpp"
#include "diffexp/banana_oracle.hpp"
#include "diffexp/box_oracle.hpp"
#include "diffexp/pentagon_oracle.hpp"
#include "diffexp/box_bubble_oracle.hpp"
#include "diffexp/double_box_oracle.hpp"
#include "diffexp/box_triangle_oracle.hpp"
#include <boost/json.hpp>
#include <charconv>
#include <iostream>
#include <chrono>
using namespace diffexp;
namespace json=boost::json;
using B=kernel::ComplexBall;
using Clock=std::chrono::steady_clock;
static double seconds(Clock::time_point start){return std::chrono::duration<double>(Clock::now()-start).count();}
static unsigned number(const char* value){unsigned out=0;std::string s(value);auto p=std::from_chars(s.data(),s.data()+s.size(),out);if(p.ec!=std::errc{} || p.ptr!=s.data()+s.size() || !out)throw std::invalid_argument("resource argument must be a positive integer");return out;}
static std::string real_text(const arb_t value){char* s=arb_get_str(value,45,0);std::string out(s);flint_free(s);return out;}
static json::object ball(const B& value){return {{"real",real_text(acb_realref(value.raw()))},{"imaginary",real_text(acb_imagref(value.raw()))}};}
int main(int argc,char** argv) {
  // Core transport diagnostics currently use cout; reserve stdout for JSON.
  std::ostream json_output(std::cout.rdbuf());
  struct DiagnosticRedirect {std::streambuf* original;DiagnosticRedirect():original(std::cout.rdbuf(std::cerr.rdbuf())){}~DiagnosticRedirect(){std::cout.rdbuf(original);}} redirect;
  const auto start=Clock::now();json::object report{{"schema","DiffExp.FTAcceptance/v1"},{"status","error"},{"absolute_threshold","1e-20"},{"comparison_is_numerical_acceptance_not_certification",true}};
  std::string phase="arguments";unsigned reused=0,built=0;json::array exact_events;
  try {
    if(argc<2 || argc>9)throw std::invalid_argument("usage: test_ft_examples FAMILY [CACHE [FIRE_OR_DASH [ENDPOINT_N [ORDINARY_N [BITS [GRAPH_TIMEOUT_SECONDS [METHOD]]]]]]]");
    const std::string name=argv[1];report["family"]=name;
    const std::vector<std::string> supported{"bubble","sunrise","banana","banana_unequal","banana4","banana4_unequal","box","pentagon","pentagon_massive","box_bubble","double_box_planar","kite","box_triangle"};
    if(std::find(supported.begin(),supported.end(),name)==supported.end())throw std::invalid_argument("family has no independent reference in this runner");
    const bool numerical_pins=name=="pentagon_massive",estimated_reference=name=="kite";
    const unsigned requested_top=numerical_pins?2:0;
    report["requested_epsilon_high"]=requested_top;
    if(numerical_pins)report["absolute_threshold"]="1e-18";
    if(estimated_reference)report["absolute_threshold"]="5e-7";
    const auto cache=argc>2?std::filesystem::path(argv[2]):std::filesystem::temp_directory_path()/"diffexp-exact-cache";
    recursion::Options exact;recursion::NumericalOptions numerical;
    AdjointConditioningStats conditioning;
    numerical.adjoint.conditioning_stats=&conditioning;
    exact.reduction.batch_cache_directory=cache/"fire-batches";
    exact.reduction.pending_cache_directory=cache/"fire-pending";
    numerical.endpoint_cache_directory=cache/"affine-series";numerical.ordinary_cache_directory=cache/"ordinary-transport";
    if(argc>3 && std::string(argv[3])!="-")exact.reduction.provider.executable=argv[3];
    if(argc>4)numerical.endpoint_order=number(argv[4]);
    if(argc>5)numerical.ordinary_order=number(argv[5]);
    if(argc>6)numerical.working_bits=number(argv[6]);
    if(argc>7)exact.total_timeout_seconds=number(argv[7]);
    const std::string method=argc>8?argv[8]:"adjoint";
    if(method!="values" && method!="adjoint" && method!="factored" && method!="auto")
      throw std::invalid_argument("METHOD must be adjoint, factored, auto or values");
    numerical.observable_adjoint=method!="values";report["method"]=method;
    numerical.linear_method=method=="factored"?recursion::LinearMethod::factored:
      method=="auto"?recursion::LinearMethod::automatic:recursion::LinearMethod::adjoint;
    if(numerical.endpoint_order>1000 || numerical.ordinary_order>1000 || numerical.working_bits<64 || numerical.working_bits>1000000 || exact.total_timeout_seconds>86400)
      throw std::invalid_argument("resource argument outside bounded acceptance limits");
    report["resources"]={{"endpoint_order",numerical.endpoint_order},{"ordinary_order",numerical.ordinary_order},{"working_bits",numerical.working_bits},{"leaf_digits",numerical.leaf_digits},{"graph_timeout_seconds",exact.total_timeout_seconds},{"max_refinements",numerical.max_refinements}};
    report["algorithms"]={{"shared_leaf_sources",numerical.observable_adjoint},{"frobenius_finite_lag",numerical.endpoint.finite_lag_recurrence},{"numeric_operator_convolution",numerical.matching.numeric_operator_convolution},{"polynomial_adjoint",numerical.adjoint.polynomial_recurrence},{"adjoint_conditioning_guard",true},{"verified_endpoint_cache",!numerical.endpoint_cache_directory.empty()}};
    report["algorithms"].as_object()["linear_transport"]=method;
    report["algorithms"].as_object()["centered_chart_action"]=numerical.adjoint.centered_input;
    report["resources"].as_object()["max_conditioning_halvings"]=numerical.adjoint.max_conditioning_halvings;
    report["cache_directory"]=std::filesystem::absolute(cache).string();
    artifacts::Store store(cache);
    auto provider=[&](const auto& basis,const auto& dimension,const auto& field,auto parameter,const auto& sources,const auto& budget) {
      level::Provider fallback;
      if(argc>3 && std::string(argv[3])=="-")fallback=[](const auto&,const auto&){fire::Result result;result.reason="acceptance cache-only mode forbids FIRE";return result;};
      auto result=cached_level::prepare(store,basis,dimension,field,parameter,sources,budget,fallback);
      if(result.cache_hit)++reused;else if(result.result.success)++built;
      return std::move(result.result);
    };
    phase="exact graph";auto t=Clock::now();
    std::cerr<<"exact graph: "<<name<<'\n';
    auto family=feynman::example_family(name);
    auto graph=recursion::prepare(family,{},exact,provider,[&](const auto& event){
      std::cerr<<"exact level "<<event.depth<<": "<<event.masters<<" masters, "<<event.elapsed_seconds<<" s\n";
      exact_events.push_back(json::object{{"depth",event.depth},{"physical_propagators",event.physical_propagators},{"requests",event.requests},{"sources",event.sources},{"masters",event.masters},{"provider_passes",event.provider_passes},{"seconds",event.elapsed_seconds},{"scalar_leaf",event.scalar_leaf}});
    });
    report["exact_seconds"]=seconds(t);report["levels"]=graph.nodes.size();
    phase="FT numerical evaluation";t=Clock::now();
    numerical.progress=[&](std::size_t depth,const std::string& text,int high){std::cerr<<seconds(start)<<"s numerical level "<<depth<<": "<<text<<", epsilon through "<<high<<'\n';};
    recursion::Evaluator evaluator(graph,numerical);auto result=evaluator.evaluate(requested_top);
    report["numerical_seconds"]=seconds(t);
    const auto& stats=evaluator.statistics();report["numerical_statistics"]={{"exact_plans",stats.exact_plans},{"leaf_evaluations",stats.leaf_evaluations},{"numeric_evaluations",stats.numeric_evaluations},{"refinements",stats.refinements},{"cache_hits",stats.cache_hits},{"endpoint_series_built",stats.endpoint_series_built},{"endpoint_series_reused",stats.endpoint_series_reused}};
    report["numerical_statistics"].as_object()["adjoint_selections"]=stats.adjoint_selections;
    report["numerical_statistics"].as_object()["factored_selections"]=stats.factored_selections;
    report["numerical_statistics"].as_object()["demand_preflights"]=stats.demand_preflights;
    report["numerical_statistics"].as_object()["local_operator_reuses"]=stats.local_operator_reuses;
    report["numerical_statistics"].as_object()["conditioning_method_fallbacks"]=stats.conditioning_method_fallbacks;
    report["ordinary_checkpoints"]={{"loaded",stats.ordinary_checkpoints.loaded},{"saved",stats.ordinary_checkpoints.saved},{"completed_reused",stats.ordinary_checkpoints.completed_reused},{"omitted_tail_certified",false}};
    report["adjoint_conditioning"]={{"polynomial_charts",conditioning.polynomial_charts},{"rational_cross_checks",conditioning.rational_cross_checks},{"rational_compilations",conditioning.rational_compilations},
      {"centered_charts",conditioning.centered_charts},{"homogeneous_chart_maps",conditioning.homogeneous_chart_maps},
      {"centered_budget_skips",conditioning.centered_budget_skips},{"conditioning_subdivisions",conditioning.conditioning_subdivisions},
      {"polynomial_homogeneous_columns",conditioning.polynomial_homogeneous_columns},{"rational_homogeneous_columns",conditioning.rational_homogeneous_columns}};
    report["ft_taylor_tail_certified"]=result.taylor_tail_certified;
    report["reported_ft_radii_include_omitted_tails"]=result.taylor_tail_certified;
    if(result.values.size()!=1 || result.high()<static_cast<int>(requested_top))throw std::runtime_error("unexpected scalar FT output shape/window");
    phase="independent reference";t=Clock::now();std::cerr<<phase<<": "<<name<<'\n';B::set_precision(numerical.working_bits);
    std::map<int,B> references;
    if(estimated_reference) {
      references.emplace(0,B::from_strings("223983917019259/1000000000000000"));
      report["reference_kind"]="independent numerical quadrature with estimated error";
      report["reference_source"]="Docs/Results.md: Massive kite";
      report["reference_reliability"]="heuristic/conservative absolute error estimate, not a rigorous enclosure; no high-digit validation";
      report["reference_estimated_absolute_error"]="1.8e-7";
      report["reference_ball_meaning"]="rounding enclosure of supplied midpoint only; independent quadrature error is recorded separately";
    } else if(numerical_pins) {
      references.emplace(0,B::from_strings("1813378668630195764/100000000000000000000"));
      references.emplace(1,B::from_strings("7613115416144053565/1000000000000000000000"));
      references.emplace(2,B::from_strings("5214475578477681142/1000000000000000000000"));
      report["reference_kind"]="independent convergent numerical Feynman-parameter integration pins";
      report["reference_source"]="Scripts/pentagon_massive_oracle.m: pmOracleReference[]";
      report["reference_reliability"]="numerical pins only, no rigorous error balls; comparison capped at18 absolute decimal digits";
      report["reference_comparison_digits"]=18;
      report["reference_ball_meaning"]="rounding enclosure of each supplied decimal pin, not a certified enclosure of the true coefficient";
    } else if(name=="box_triangle") {
      auto ref=oracle::box_triangle_reference(256);
      for(int k=-4;k<=0;++k)references.emplace(k,ref.at(k));
      report["reference_kind"]="independent original-IBP and Mellin-Barnes/Cauchy numerical reference";
      report["reference_source"]="Smirnov-Veretin hep-ph/9907385 eq21 and docs/box-triangle-original-ibp";
      report["reference_reliability"]="estimated convergence, not a rigorous coefficient enclosure";
      report["reference_estimated_absolute_error"]="1e-30";
      report["reference_ball_meaning"]="arithmetic error only; MB truncation and Cauchy alias error estimated separately";
      json::array diagnostics;for(const auto& difference:ref.resolution_differences)diagnostics.push_back(ball(difference));
      report["reference_resolution_differences"]=std::move(diagnostics);
    } else if(name=="double_box_planar") {
      auto ref=oracle::double_box_planar_reference(256);
      for(int k=-4;k<=0;++k)references.emplace(k,ref.at(k));
      report["reference_kind"]="independent analytic Smirnov planar double-box Laurent coefficients";
    } else if(name=="box_bubble") {
      auto ref=oracle::box_bubble_reference(Rational(-1),0,256);
      for(int k=-3;k<=0;++k)references.emplace(k,ref.at(k));
      report["reference_kind"]="independent loop-by-loop box-bubble gamma Laurent coefficients";
    } else if(name=="pentagon") {
      auto ref=oracle::massless_pentagon_reference(family.momenta,256);
      for(int k=-2;k<=0;++k)references.emplace(k,ref.at(k));
      report["reference_kind"]="independent analytic dimension-shift pentagon Laurent coefficients";
    } else if(name=="box") {
      auto ref=oracle::massless_box_reference(Rational(-1),Rational("-1/3"),256);
      for(int k=-2;k<=0;++k)references.emplace(k,ref.at(k));
      report["reference_kind"]="independent analytic massless box Laurent coefficients";
    } else {
      std::vector<Rational> masses;for(const auto& line:family.momenta.lines)masses.push_back(line.mass_squared);
      oracle::BananaOptions options;options.target_bits=100;options.working_bits=256;
      auto ref=oracle::banana_bessel(masses,options);references.emplace(0,ref.value);
      report["reference_kind"]="independent certified coordinate-space Bessel epsilon-zero coefficient";
      report["oracle_evaluations"]=ref.evaluations;report["oracle_analytic_queries"]=ref.analytic_queries;
    }
    report["reference_seconds"]=seconds(t);report["reference_coefficient_enclosures_certified"]=!numerical_pins && !estimated_reference && name!="box_triangle";
    phase="comparison";const B tolerance=B::from_strings(estimated_reference?"1/2000000":numerical_pins?"1/1000000000000000000":"1/100000000000000000000");
    const auto limit=NativeTailMagnitude::lower_abs(tolerance);bool pass=true;json::array comparisons,poles;
    auto coefficient=[&](int k){if(k<result.low)return B(0);if(k>result.high())throw std::runtime_error("FT upper coefficient missing");return result.values[0].at(k-result.low);};
    for(const auto& [order,ref]:references) {
      auto actual=coefficient(order),difference=actual-ref;
      bool good=actual.is_finite() && ref.is_finite() && NativeTailMagnitude::upper_abs(difference)<=limit;pass=pass&&good;
      comparisons.push_back(json::object{{"epsilon_order",order},{"ft",ball(actual)},
        {"solver_real_arithmetic_radius_exponent",actual.real_radius_exponent()},
        {"solver_imaginary_arithmetic_radius_exponent",actual.imag_radius_exponent()},
        {"radius_exponent_meaning","radius is strictly less than2^exponent, or exactly zero; omitted tails excluded"},
        {"reference",ball(ref)},{"difference",ball(difference)},{"pass",good}});
    }
    // Finite bananas must have no negative Laurent coefficient. The box permits
    // eps^-2 and eps^-1 (as does the pentagon), but no more singular terms. Declared low is structural.
    const int first_allowed=(name=="double_box_planar" || name=="box_triangle")?-4:name=="box_bubble"?-3:(name=="box" || name=="pentagon")?-2:0;
    report["forbidden_pole_absolute_threshold"]="1e-20";
    const auto pole_limit=NativeTailMagnitude::lower_abs(B::from_strings("1/100000000000000000000"));
    for(int k=result.low;k<first_allowed;++k) {
      auto value=coefficient(k);bool good=value.is_finite() && NativeTailMagnitude::upper_abs(value)<=pole_limit;pass=pass&&good;
      poles.push_back(json::object{{"epsilon_order",k},{"coefficient",ball(value)},{"pass",good}});
    }
    report["comparisons"]=std::move(comparisons);report["forbidden_pole_audit"]=std::move(poles);
    report["epsilon_low"]=result.low;report["epsilon_high"]=result.high();report["status"]=pass?"pass":"fail";
    report["exact_systems_built"]=built;report["exact_systems_reused"]=reused;report["exact_events"]=std::move(exact_events);report["total_seconds"]=seconds(start);
    json_output<<json::serialize(report)<<'\n';return pass?0:1;
  } catch(const std::exception& e) {
    report["phase"]=phase;report["error"]=e.what();report["exact_systems_built"]=built;report["exact_systems_reused"]=reused;report["exact_events"]=std::move(exact_events);report["total_seconds"]=seconds(start);
    std::cerr<<phase<<": "<<e.what()<<'\n';json_output<<json::serialize(report)<<'\n';return 1;
  }
}
