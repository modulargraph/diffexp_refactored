// Explicit staged composition; saved numeric reports are uncertified inputs.
#include "stage_probe_io.hpp"
#include "diffexp/level_cache.hpp"
#include "diffexp/banana_oracle.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
using namespace diffexp;using namespace diffexp::stage_probe;
namespace json=boost::json;
struct Input {
  std::string bytes;json::object report;
  explicit Input(const std::filesystem::path& path) {
    if(std::filesystem::file_size(path)>64*1024*1024)throw std::length_error("stage report file budget");
    std::ifstream in(path);bytes.assign(std::istreambuf_iterator<char>(in),{});report=json::parse(bytes).as_object();
  }
  std::string hash()const{return artifacts::detail::sha256(bytes);}
};
std::string string(const json::value& value){return artifacts::detail::string(value);}
json::array basis(const std::vector<ibp::Integral>& values) {
  json::array result;for(const auto& integral:values){json::array row;for(auto index:integral)row.emplace_back(index);result.push_back(std::move(row));}return result;
}
bool equal(const ExactEpsilonMatrix& a,const ExactEpsilonMatrix& b) {
  if(a.size()!=b.size())return false;
  for(unsigned i=0;i<a.size();++i){if(a[i].size()!=b[i].size())return false;
    for(unsigned j=0;j<a[i].size();++j)if(!(a[i][j]==a[i][j].parse(b[i][j].str())))return false;}
  return true;
}
std::string part(const arb_t value){char* raw=arb_get_str(value,50,0);std::string result(raw);flint_free(raw);return result;}
json::object display(const B& value){return {{"real",part(acb_realref(value.raw()))},{"imaginary",part(acb_imagref(value.raw()))}};}
int main(int argc,char** argv) {
  const auto started=std::chrono::steady_clock::now();int code=0;
  json::object result{{"schema","DiffExp3.StagedRecursionComposition/v1"},{"status","error"},
    {"full_evaluator_end_to_end",false},{"omitted_tails_included",false}};
  try {
    if(argc!=6)throw std::invalid_argument("usage: compose_recursion_stage CACHE PREPARED_REPORT LOWER_REPORT UPPER_REPORT CHILD_REPORT");
    Input prepared(argv[2]),lower(argv[3]),upper(argv[4]),child(argv[5]);auto stage=read_stage(prepared.report);
    const auto family=string(prepared.report.at("family"));const auto depth=artifacts::detail::integer(prepared.report.at("depth"));
    const auto& resources=prepared.report.at("resources");
    if(depth<0 || depth>15 || prepared.report.at("status")!="completed_diagnostic" ||
        child.report.at("status")!="completed_diagnostic" || string(child.report.at("family"))!=family ||
        artifacts::detail::integer(child.report.at("depth"))!=depth+1 || child.report.at("resources")!=resources)
      throw std::invalid_argument("prepared/child report identity or numerical resources mismatch");
    B::set_precision(artifacts::detail::integer(resources.as_object().at("working_bits")));
    const auto arm=[&](const Input& input,const char* side) {
      const auto& r=input.report;
      if(r.at("schema")!="DiffExp3.AdjointStageReplay/v1" || r.at("status")!="completed_diagnostic" ||
          r.at("side")!=side || r.at("resources")!=resources || string(r.at("input_report_sha256"))!=prepared.hash() ||
          r.at("prepared_stage_sha256")!=prepared.report.at("prepared_stage_sha256"))
        throw std::invalid_argument("completed arm does not belong to this prepared stage");
      auto rows=read_rows(r.at("result_rows"));
      const auto& endpoint=std::string(side)=="lower"?stage.lower_endpoint:stage.upper_endpoint;
      const auto& forcing=std::string(side)=="lower"?stage.lower_forcing:stage.upper_forcing;
      const auto [unused,epsilon]=path_epsilon_variables(stage.lower_path.front());(void)unused;
      int low=std::min(0,endpoint.low);
      for(const auto& row:forcing)for(const auto& entry:row)if(!entry.is_zero())
        low=std::min(low,static_cast<int>(*exact_epsilon_valuation(entry,epsilon)));
      if(rows.low!=low || rows.high!=endpoint.high || rows.columns()!=endpoint.columns() || rows.coefficients.size()!=endpoint.coefficients.size())
        throw std::invalid_argument("completed arm dimensions or coefficient top mismatch");
      return rows;
    };
    auto left=arm(lower,"lower"),right=arm(upper,"upper");
    artifacts::Store store(argv[1]);
    auto provider=[&](const auto& basis,const auto& dimension,const auto& field,auto xi,const auto& sources,const auto& budget) {
      auto hit=cached_level::prepare(store,basis,dimension,field,xi,sources,budget,[](const auto&,const auto&)->fire::Result{throw std::runtime_error("composition forbids FIRE and exact cache misses");});
      if(!hit.cache_hit || !hit.result.success)throw std::runtime_error("verified exact graph unavailable");return hit.result;
    };
    auto original=feynman::example_family(family);auto graph=recursion::prepare(original,{},recursion::Options{},provider);
    if(depth+1>=graph.nodes.size())throw std::invalid_argument("prepared stage has no child in the exact graph");
    const auto& node=graph.nodes[depth];auto [xi,ei]=path_epsilon_variables(graph.dimension);
    if(prepared.report.at("prepared_stage").as_object().at("ordered_master_basis")!=basis(node.closure.ordered_basis) ||
        child.report.at("requested_basis")!=basis(graph.nodes[depth+1].requested) ||
        child.report.at("requested_basis")!=basis(node.closure.ordered_basis) ||
        string(child.report.at("dimension"))!=graph.dimension.str())
      throw std::invalid_argument("ordered basis or dimension does not match the verified graph");
    auto gauge=epsilon_diagonal_gauge(node.closure.matrix,ei);
    if(gauge.shifts!=stage.epsilon_gauge_shifts || !equal(gauge.matrix,stage.connection))
      throw std::invalid_argument("prepared connection or epsilon gauge differs from the verified graph");
    auto zero=graph.dimension.constant(0),x=graph.dimension.variable(xi),one=x.constant(1),anchor=x.constant(node.anchor);
    ExactEpsilonMatrix forcing(node.operations.size(),std::vector<Exact>(gauge.matrix.size(),zero)),direct=forcing;
    for(unsigned i=0;i<node.operations.size();++i) {
      const auto& op=node.operations[i];
      if(op.operation==feynman::Operation::BetaIntegral) {
        forcing[i]=node.observable_rows[i];
        for(auto& entry:forcing[i])entry=entry*x.constant(op.normalization)*x.pow(op.left_power)*(one-x).pow(op.right_power);
      } else if(op.operation==feynman::Operation::Direct)direct[i]=node.observable_rows[i];
    }
    forcing=gauge_observable_columns(gauge,forcing);auto negative=forcing;for(auto& row:negative)for(auto& entry:row)entry=-entry;
    if(!equal(forcing,stage.lower_forcing) || !equal(negative,stage.upper_forcing) ||
        stage.lower_path.back().str()!=anchor.str() || stage.upper_path.back().str()!=anchor.str())
      throw std::invalid_argument("prepared forcing or destination differs from the verified graph");
    auto local=shift_laurent_columns(add_laurent_rows(left,right),stage.epsilon_gauge_shifts,true);
    local=add_laurent_rows(local,exact_laurent_rows(direct,anchor,local.high));
    auto inner=read_expression(child.report.at("linear_expression"),string(child.report.at("linear_expression_sha256")));
    if(artifacts::detail::integer(child.report.at("requested_high"))< -local.low)
      throw std::invalid_argument("child report does not cover the parent's epsilon demand");
    auto expression=linear_boundary::compose(local,inner,-inner.leaf_source->low);
    auto values=linear_boundary::materialize(expression,0);
    LaurentRows output{values.low,values.high(),{}};for(const auto& row:values.values)output.coefficients.push_back({row});
    result["family"]=family;result["depth"]=depth;result["resources"]=resources;
    result["local_arithmetic"]=quality(local);result["composed_arithmetic"]=quality(expression.transform);
    result["output_arithmetic"]=quality(output);result["output_rows"]=exact_rows(output);
    result["inputs"]={{"prepared_sha256",prepared.hash()},{"lower_sha256",lower.hash()},{"upper_sha256",upper.hash()},{"child_sha256",child.hash()}};
    result["verified_graph_binding"]=true;result["shared_leaf_preserved"]=true;result["status"]="completed_diagnostic";
    const std::vector<std::string> finite_bananas{"sunrise","banana","banana_unequal","banana4","banana4_unequal"};
    if(depth==0 && std::find(finite_bananas.begin(),finite_bananas.end(),family)!=finite_bananas.end()) {
      if(values.values.size()!=1 || values.low>0)throw std::invalid_argument("finite banana output dimensions");
      std::vector<Rational> masses;for(const auto& line:original.momenta.lines)masses.push_back(line.mass_squared);
      oracle::BananaOptions options;options.target_bits=100;options.working_bits=256;
      auto reference=oracle::banana_bessel(masses,options);
      const auto difference=values.values[0][-values.low]-reference.value;
      const auto threshold=NativeTailMagnitude::lower_abs(B::from_strings("1/100000000000000000000"));
      bool pass=NativeTailMagnitude::upper_abs(difference)<=threshold;
      result["finite_part"]=display(values.values[0][-values.low]);result["reference"]=display(reference.value);result["difference"]=display(difference);
      result["absolute_threshold"]="1e-20";result["reference_kind"]="independent certified coordinate-space Bessel epsilon-zero coefficient";
      result["finite_part_pass"]=pass;json::array poles;
      for(int k=values.low;k<0;++k) {
        auto value=values.values[0][k-values.low];const bool good=NativeTailMagnitude::upper_abs(value)<=threshold;
        pass&=good;poles.push_back(json::object{{"epsilon_order",k},{"value",display(value)},{"pass",good}});
      }
      result["forbidden_poles"]=std::move(poles);result["status"]=pass?"pass":"fail";
      if(!pass)code=1;
    }
  }catch(const std::exception& e){result["error"]=e.what();code=1;}
  result["native_elapsed_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
  std::cout<<json::serialize(result)<<'\n';return code;
}
