// Explicit replay of an uncertified checkpoint, not an automatic trusted cache.
#include "stage_probe_io.hpp"
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
using namespace diffexp;using namespace diffexp::stage_probe;
namespace json=boost::json;
int main(int argc,char** argv) {
  const auto begin=std::chrono::steady_clock::now();const auto cpu=std::clock();
  json::object output{{"schema","DiffExp3.AdjointStageReplay/v1"},{"status","error"},{"independent_reference_comparison",false}};
  AdjointConditioningStats stats;int code=0;
  try {
    if(argc<3 || argc>5)throw std::invalid_argument("usage: replay_adjoint_stage PREPARED_REPORT lower|upper [current|baseline [CHECKPOINT]]");
    const std::filesystem::path file=argv[1];const std::string side=argv[2],mode=argc>3?argv[3]:"current";
    if((side!="lower" && side!="upper") || (mode!="current" && mode!="baseline"))throw std::invalid_argument("replay side/method");
    if(std::filesystem::file_size(file)>64*1024*1024)throw std::length_error("prepared report file budget");
    std::ifstream input(file);std::string text((std::istreambuf_iterator<char>(input)),{});
    auto report=json::parse(text).as_object();auto stage=read_stage(report);
    const auto& resources=report.at("resources").as_object();const auto bits=artifacts::detail::integer(resources.at("working_bits"));
    const auto order=artifacts::detail::integer(resources.at("ordinary_order"));
    if(bits<64 || order<8 || order>1000)throw std::invalid_argument("replay numerical resource budget");
    B::set_precision(bits);AdjointOptions options;options.taylor_order=order;options.conditioning_stats=&stats;
    if(mode=="baseline"){options.centered_input=false;options.max_conditioning_halvings=0;}
    auto initial=side=="lower"?stage.lower_endpoint:stage.upper_endpoint;
    auto vertices=side=="lower"?stage.lower_path:stage.upper_path;
    const auto& forcing=side=="lower"?stage.lower_forcing:stage.upper_forcing;
    const auto [xi,ei]=path_epsilon_variables(vertices.front());(void)xi;
    int expected_low=std::min(0,initial.low);
    for(const auto& row:forcing)for(const auto& coefficient:row)if(!coefficient.is_zero())
      expected_low=std::min(expected_low,static_cast<int>(*exact_epsilon_valuation(coefficient,ei)));
    const std::filesystem::path checkpoint=argc>4?argv[4]:"";
    const auto source_hash=artifacts::detail::sha256(text);
    constexpr auto algorithm="conditioned-centered-adjoint-v1";
    output["resumed_checkpoint"]=false;
    if(!checkpoint.empty() && std::filesystem::exists(checkpoint)) {
      if(std::filesystem::file_size(checkpoint)>64*1024*1024)throw std::length_error("arm checkpoint file budget");
      std::ifstream saved(checkpoint);std::string bytes((std::istreambuf_iterator<char>(saved)),{});
      auto state=json::parse(bytes).as_object();auto hash=artifacts::detail::string(state.at("checkpoint_sha256"));state.erase("checkpoint_sha256");
      if(artifacts::detail::sha256(artifacts::detail::canonical(state))!=hash || state.at("schema")!="DiffExp3.AdjointArmCheckpoint/v1" ||
          artifacts::detail::string(state.at("source_report_sha256"))!=source_hash ||
          artifacts::detail::string(state.at("side"))!=side || artifacts::detail::string(state.at("mode"))!=mode ||
          state.at("algorithm")!=algorithm || state.at("resources")!=resources || state.at("omitted_tails_included")!=false)
        throw std::invalid_argument("arm checkpoint identity/configuration mismatch");
      auto restored=read_rows(state.at("rows"));
      if(restored.low!=expected_low || restored.high!=initial.high || restored.columns()!=initial.columns() || restored.coefficients.size()!=initial.coefficients.size())
        throw std::invalid_argument("arm checkpoint shape/window mismatch");
      const auto& path=state.at("remaining_path").as_array();
      if(path.empty() || path.size()>32)throw std::invalid_argument("arm checkpoint remaining path budget");
      // Parse through the existing field so equation and path contexts agree.
      const auto parse=[&](const std::string& expression) {
        return stage.connection[0][0].parse(expression);
      };
      std::vector<Exact> remaining;
      for(const auto& vertex:path)remaining.push_back(parse(artifacts::detail::string(vertex)));
      if(!(remaining.back()==vertices.back()))throw std::invalid_argument("arm checkpoint changed the destination");
      vertices=std::move(remaining);initial=std::move(restored);output["resumed_checkpoint"]=true;
    }
    if(!checkpoint.empty()) {
      options.chart_observer=[&](unsigned leg,double,double to,const LaurentRows& rows) {
        auto rest=remaining_path(vertices,leg,to);json::array path;for(const auto& vertex:rest)path.emplace_back(vertex.str());
        json::object state{{"schema","DiffExp3.AdjointArmCheckpoint/v1"},{"algorithm",algorithm},
          {"source_report_sha256",source_hash},{"side",side},{"mode",mode},{"resources",resources},
          {"remaining_path",path},{"rows",exact_rows(rows)},{"omitted_tails_included",false}};
        state["checkpoint_sha256"]=artifacts::detail::sha256(artifacts::detail::canonical(state));
        auto temporary=checkpoint;temporary+=".new";
        {std::ofstream out(temporary);out<<artifacts::detail::canonical(state)<<'\n';if(!out)throw std::runtime_error("cannot save completed adjoint chart");}
        std::filesystem::rename(temporary,checkpoint);
        auto diagnostic=quality(rows);diagnostic["leg"]=leg;diagnostic["parameter"]=to;
        diagnostic["centered_charts"]=stats.centered_charts;diagnostic["subdivisions"]=stats.conditioning_subdivisions;
        std::cerr<<json::serialize(diagnostic)<<'\n'<<std::flush;
      };
      output["checkpoint_path"]=std::filesystem::absolute(checkpoint).string();
    }
    output["initial_arithmetic"]=quality(initial);output["side"]=side;output["mode"]=mode;output["resources"]=resources;
    output["prepared_stage_sha256"]=report.at("prepared_stage_sha256");output["input_report_sha256"]=artifacts::detail::sha256(text);
    output["input_report"]=std::filesystem::absolute(file).string();
    std::cerr<<"Replaying "<<side<<" arm at N"<<order<<" and "<<bits<<" bits; no endpoint or child recomputation\n"<<std::flush;
    auto result=transport_adjoint_rows(stage.connection,std::move(initial),forcing,vertices,options);
    output["result_arithmetic"]=quality(result);output["result_rows"]=exact_rows(result);
    output["status"]="completed_diagnostic";
  }catch(const std::exception& error){output["error"]=error.what();code=1;}
  output["conditioning"]={{"polynomial_charts",stats.polynomial_charts},{"rational_cross_checks",stats.rational_cross_checks},{"centered_charts",stats.centered_charts},{"homogeneous_maps",stats.homogeneous_chart_maps},{"subdivisions",stats.conditioning_subdivisions},{"polynomial_homogeneous_columns",stats.polynomial_homogeneous_columns},{"rational_homogeneous_columns",stats.rational_homogeneous_columns}};
  output["native_elapsed_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
  output["process_cpu_seconds"]=static_cast<double>(std::clock()-cpu)/CLOCKS_PER_SEC;
  std::cout<<json::serialize(output)<<'\n';return code;
}
