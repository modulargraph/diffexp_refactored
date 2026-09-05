// Opt-in diagnostic, not an independent integral acceptance test.
#include "diffexp/recursion_pipeline.hpp"
#include "diffexp/level_cache.hpp"
#include "stage_probe_io.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
using namespace diffexp;
namespace json=boost::json;
using B=Jet::Ball;
using namespace diffexp::stage_probe;
int main(int argc,char** argv) {
  json::object report{{"schema","DiffExp3.RecursionStageProbe/v1"},{"status","error"},{"independent_reference_comparison",false}};
  const auto begin=std::chrono::steady_clock::now();const auto cpu=std::clock();
  const auto elapsed=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();};
  json::array observations;AdjointConditioningStats conditioning;int code=0;
  try {
    if(argc<5 || argc>6)throw std::invalid_argument("usage: probe_recursion_stage CACHE FAMILY ZERO_BASED_DEPTH EPSILON_HIGH [adjoint|factored|auto|endpoints]");
    const std::filesystem::path cache=argv[1];const std::string family=argv[2],method=argc>5?argv[5]:"adjoint";
    const auto depth=std::stoul(argv[3]),high=std::stoul(argv[4]);
    if(depth>15 || high>100 || (method!="adjoint" && method!="factored" && method!="auto" && method!="endpoints"))throw std::invalid_argument("stage/method budget");
    report["family"]=family;report["depth"]=depth;report["requested_high"]=high;report["method"]=method;
    artifacts::Store store(cache);recursion::Options exact;
    auto provider=[&](const auto& basis,const auto& dimension,const auto& field,auto xi,const auto& sources,const auto& budget) {
      auto hit=cached_level::prepare(store,basis,dimension,field,xi,sources,budget,[](const auto&,const auto&)->fire::Result{throw std::runtime_error("stage probe forbids FIRE and exact cache misses");});
      if(!hit.cache_hit || !hit.result.success)throw std::runtime_error("verified exact graph unavailable");return hit.result;
    };
    auto graph=recursion::prepare(feynman::example_family(family),{},exact,provider);
    if(depth>=graph.nodes.size() || (method=="endpoints" && graph.nodes[depth].scalar_leaf))
      throw std::invalid_argument("stage depth is unavailable or a scalar leaf has no endpoint problem");
    json::array requested;
    for(const auto& integral:graph.nodes[depth].requested){json::array indices;for(auto power:integral)indices.emplace_back(power);requested.push_back(std::move(indices));}
    report["requested_basis"]=std::move(requested);report["dimension"]=graph.dimension.str();
    recursion::NumericalOptions options;options.endpoint_cache_directory=cache/"affine-series";options.ordinary_cache_directory=cache/"ordinary-transport";
    options.linear_method=method=="factored"?recursion::LinearMethod::factored:method=="auto"?recursion::LinearMethod::automatic:recursion::LinearMethod::adjoint;
    options.adjoint.conditioning_stats=&conditioning;
    options.progress=[&](auto d,const auto& phase,int order){std::cerr<<elapsed()<<"s depth="<<d<<" "<<phase<<" high="<<order<<'\n'<<std::flush;};
    options.operator_observer=[&](auto d,const auto& phase,const auto& rows) {
      auto item=quality(rows);item["depth"]=d;item["phase"]=phase;item["native_elapsed_seconds"]=elapsed();
      std::cerr<<json::serialize(item)<<'\n'<<std::flush;observations.push_back(std::move(item));
    };
    report["resources"]={{"endpoint_order",options.endpoint_order},{"ordinary_order",options.ordinary_order},{"working_bits",options.working_bits},{"leaf_digits",options.leaf_digits}};
    recursion::Evaluator solver(graph,options);
    if(method=="endpoints") {
      auto stage=solver.prepare_adjoint_stage(depth,high);
      B::set_precision(options.working_bits);
      report["lower_endpoint_arithmetic"]=quality(stage.lower_endpoint);
      report["upper_endpoint_arithmetic"]=quality(stage.upper_endpoint);
      report["epsilon_high_meaning"]="operator coefficient top before final inverse epsilon gauge; no child boundary evaluated";
      auto payload=stage_payload(stage,graph,depth);
      report["prepared_stage_sha256"]=artifacts::detail::sha256(artifacts::detail::canonical(payload));
      report["prepared_stage"]=std::move(payload);
    } else {
      auto value=solver.evaluate(depth,high);
      LaurentRows result{value.low,value.high(),{}};for(const auto& row:value.values)result.coefficients.push_back({row});
      B::set_precision(options.working_bits);report["output_arithmetic"]=quality(result);
      report["output_rows"]=exact_rows(result);
      auto expression=exact_expression(solver.linear_expression(depth));
      report["linear_expression_sha256"]=artifacts::detail::sha256(artifacts::detail::canonical(expression));
      report["linear_expression"]=std::move(expression);
    }
    const auto& statistics=solver.statistics();
    report["statistics"]={{"plans",statistics.exact_plans},{"leaf_evaluations",statistics.leaf_evaluations},{"refinements",statistics.refinements},{"preflights",statistics.demand_preflights},{"local_map_reuses",statistics.local_operator_reuses},{"conditioning_method_fallbacks",statistics.conditioning_method_fallbacks},{"factored_selections",statistics.factored_selections},{"adjoint_selections",statistics.adjoint_selections}};
    report["status"]="completed_diagnostic";
  }catch(const std::exception& error){report["error"]=error.what();code=1;}
  report["observations"]=std::move(observations);
  report["conditioning"]={{"polynomial_charts",conditioning.polynomial_charts},{"rational_cross_checks",conditioning.rational_cross_checks},{"centered_charts",conditioning.centered_charts},{"homogeneous_maps",conditioning.homogeneous_chart_maps},{"subdivisions",conditioning.conditioning_subdivisions},{"polynomial_homogeneous_columns",conditioning.polynomial_homogeneous_columns},{"rational_homogeneous_columns",conditioning.rational_homogeneous_columns}};
  report["native_elapsed_seconds"]=elapsed();report["process_cpu_seconds"]=static_cast<double>(std::clock()-cpu)/CLOCKS_PER_SEC;
  std::cout<<json::serialize(report)<<'\n';return code;
}
