#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace json = boost::json;

namespace {

json::object request(const std::string& text) {
  return json::parse(diffexp2::run_recurrence_json(text)).as_object();
}

std::string prepare_scalar_chart(const std::string& session,
                                 const std::string& key,
                                 const std::string& identity) {
  const auto response = request(std::string(R"json({
    "schema":2,"op":"chart.prepare","session":")json") + session +
    R"json(","key":")json" + key + R"json(","identity":")json" +
    identity + R"json(",
    "analytic":{
      "geometry":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,"prescriptions":[]},
      "principal_matrix":[[{"exact":"0","proven_zero":true}]],
      "native_scc_capabilities":{"regular":true,"identity_gauge":true,
        "identity_v":true,"no_pseudo":true}},
    "scc":{"components":[[0]],"structural_edges":[],
      "condensation_edges":[],"topological_order":[0],"coupling_depth":0},
    "problem":{"domain":"rational","d":1,"fb":0,"w":3,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],"val":[null]}],
      "d0_inverse":"1","blocks":[[0]],
      "assembly":{"identity":true,"poly":[],"rat":[],"val":[0]},
      "chop_digits":0}
  })json");
  if (response.at("status") != "ok")
    throw std::runtime_error(json::serialize(response));
  return std::string(response.at("chart").as_string());
}

std::string regular_run(bool seed) {
  return std::string(R"json({
    "nmax":7,"p":0,"has_initial":)json") + (seed ? "true" : "false") +
    R"json(,"adaptive_probe":false,
    "a_target":"0","b_target":"0","a_shift_min":0,
    "a_shifts":["0","1","2","3","4","5","6","7"],
    "schedule":[
      [{"case":"R","da":"0","db":"0"}],
      [{"case":"T","da":"1","db":"0"}],
      [{"case":"T","da":"2","db":"0"}],
      [{"case":"T","da":"3","db":"0"}],
      [{"case":"T","da":"4","db":"0"}],
      [{"case":"T","da":"5","db":"0"}],
      [{"case":"T","da":"6","db":"0"}],
      [{"case":"T","da":"7","db":"0"}]],
    "initial":)json" + (seed ? R"json(["1","0","0"])json"
                                  : R"json(["0","0","0"])json") +
    R"json(,"initial_validity":)json" + (seed ? "[2]" : "[null]") +
    R"json(,"source":null,"return_u":false})json";
}

std::string metadata(const std::string& checkpoint) {
  return std::string(R"json({
    "chart":{"center_exact":"0","scale_exact":"1",
      "radius":"2","infinite_radius":false},
    "tag":{"a":{"domain":"rational","canonical":"0"},
      "b":{"domain":"rational","canonical":"0"}},
    "prescriptions":[],"checkpoint_identity":")json") + checkpoint +
    R"json("})json";
}

double real_midpoint(const json::value& coefficient) {
  return std::stod(std::string(
      coefficient.as_array().front().as_string()));
}

}  // namespace

int main() {
  const auto created = request(R"json({
    "schema":2,"op":"session.create","domain":"rational",
    "precision_bits":256,"output_digits":30,"scc_capacity":2,
    "local_capacity":2
  })json");
  const auto session = std::string(created.at("session").as_string());
  const auto first = prepare_scalar_chart(
      session, "column-principal-0-key", "column-principal-0");
  const auto second = prepare_scalar_chart(
      session, "column-principal-1-key", "column-principal-1");

  const auto prepared = request(std::string(R"json({
    "schema":2,"op":"scc.prepare","session":")json") + session +
    R"json(","key":"two-block-column-key",
    "identity":"two-block-column-parent-v1",
    "parent":{
      "dimension":2,
      "exact_system_record":[
        [{"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true}],
        [{"exact":"eps","proven_zero":false},
         {"exact":"0","proven_zero":true}]],
      "exact_theta_record":[
        [{"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true}],
        [{"exact":"eps*t","proven_zero":false},
         {"exact":"0","proven_zero":true}]],
      "chart":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,"prescriptions":[]},
      "scc":{"components":[[0],[1]],"structural_edges":[[0,1]],
        "condensation_edges":[[0,1]],"topological_order":[0,1],
        "coupling_depth":1},
      "execution":{"mode":"BlockSequentialStrict","work_t_order":7},
      "work_contract":{"work_min":0,"requested_min":0,
        "requested_max":1,"work_complete_max":2,"public_t_order":1,
        "wolfram_coupling_depth":2}},
    "blocks":[
      {"block":0,"vertices":[0],"chart":")json" + first +
    R"json(","principal_identity":"column-principal-0","regular":true,
       "identity_gauge":true,"identity_v":true,"no_pseudo":true},
      {"block":1,"vertices":[1],"chart":")json" + second +
    R"json(","principal_identity":"column-principal-1","regular":true,
       "identity_gauge":true,"identity_v":true,"no_pseudo":true}],
    "couplings":[{
      "source_block":0,"target_block":1,
      "source_vertices":[0],"target_vertices":[1],
      "rows":1,"columns":1,"exact_identity":"block-0-to-1-column",
      "domain":"rational","symbols":[],
      "entries":[{"row":0,"column":0,
        "source_vertex":0,"target_vertex":1,
        "exact_original_entry":"eps","exact_theta_entry":"eps*t",
        "multiplier":{"epsilon_shift":1,"center_pole_order":0,
          "kernels":[
            ["0","1","0","0","0","0","0","0"],
            ["0","0","0","0","0","0","0","0"],
            ["0","0","0","0","0","0","0","0"]],
          "exact_identity":"eps*t","proven_zero":false}}]}]
  })json");
  if (prepared.at("status") != "ok") {
    std::cerr << "scc.prepare: " << json::serialize(prepared) << '\n';
    return EXIT_FAILURE;
  }
  const auto scc = std::string(prepared.at("scc").as_string());

  const auto solved = request(std::string(R"json({
    "schema":2,"op":"scc.solve_column","session":")json") + session +
    R"json(","scc":")json" + scc +
    R"json(","checkpoint_identity":"two-block-column-checkpoint",
    "seed":{"block":0,"run":)json" + regular_run(true) +
    R"json(,"metadata":)json" + metadata("seed-work") +
    R"json(},"targets":[{"block":1,"run":)json" + regular_run(false) +
    R"json(,"metadata":)json" + metadata("target-work") +
    R"json(}]})json");

  json::object evaluated;
  json::object scc_stats;
  json::object session_stats;
  if (solved.at("status") == "ok") {
    const auto local = std::string(solved.at("local").as_string());
    evaluated = request(std::string(R"json({
      "schema":2,"op":"local.evaluate","session":")json") + session +
      R"json(","local":")json" + local +
      R"json(","point":{"exact":"1/2"},
      "options":{"tail_estimate":false}})json");
    scc_stats = request(std::string(R"json({
      "schema":2,"op":"scc.stats","session":")json") + session +
      R"json(","scc":")json" + scc + R"json("})json");
    session_stats = request(std::string(R"json({
      "schema":2,"op":"session.stats","session":")json") + session +
      R"json("})json");
  }

  bool propagated_value = false;
  if (evaluated.if_contains("status") != nullptr &&
      evaluated.at("status") == "ok") {
    const auto& coefficients = evaluated.at("value").as_object()
                                   .at("coefficients").as_array();
    propagated_value = coefficients.size() == 4 &&
        std::abs(real_midpoint(coefficients[0]) - 1.0) < 1e-20 &&
        std::abs(real_midpoint(coefficients[1])) < 1e-20 &&
        std::abs(real_midpoint(coefficients[2])) < 1e-20 &&
        std::abs(real_midpoint(coefficients[3]) - 0.5) < 1e-20;
  }

  bool provenance_ok = false;
  bool diagnostics_ok = false;
  if (solved.at("status") == "ok") {
    const auto& provenance = solved.at("column_provenance").as_object();
    const auto& diagnostics = solved.at("block_diagnostics").as_array();
    provenance_ok =
        std::string(provenance.at("scc").as_string()) == scc &&
        provenance.at("scc_exact_identity") ==
            "two-block-column-parent-v1" &&
        provenance.at("seed_block") == 0 &&
        provenance.at("basis_index") == 0 &&
        std::string(provenance.at("exact_column_identity").as_string())
            .find("\"targets\"") != std::string::npos;
    diagnostics_ok = diagnostics.size() == 2 &&
        diagnostics[1].as_object().at("source_epsilon_min") == 1 &&
        diagnostics[1].as_object().at("source_epsilon_max") == 2 &&
        diagnostics[1].as_object().at("predecessors").as_array().size() == 1 &&
        diagnostics[1].as_object().at("predecessors").as_array()[0] == 0;
  }

  const bool ok = solved.at("status") == "ok" &&
      std::string(solved.at("chart").as_string()) == scc &&
      solved.at("dimension") == 2 &&
      solved.at("epsilon_min") == 0 && solved.at("epsilon_max") == 1 &&
      solved.at("taylor_complete_max") == 1 &&
      solved.at("native_retained") == true &&
      solved.at("json_coefficients") == 0 && provenance_ok &&
      diagnostics_ok && propagated_value &&
      scc_stats.at("execution_implemented") == true &&
      scc_stats.at("general_scc_execution") == false &&
      scc_stats.at("scc_column_solves") == 1 &&
      session_stats.at("locals") == 1 &&
      session_stats.at("local_solves") == 1 &&
      session_stats.at("scc_column_solves") == 1;

  if (!ok) {
    std::cerr << "scc.solve_column: " << json::serialize(solved) << '\n'
              << "local.evaluate: " << json::serialize(evaluated) << '\n'
              << "scc.stats: " << json::serialize(scc_stats) << '\n'
              << "session.stats: " << json::serialize(session_stats) << '\n';
  }

  (void)request(std::string(R"json({
    "schema":2,"op":"session.close","session":")json") + session +
    R"json("})json");
  std::cout << (ok ? "PASS" : "FAIL")
            << ": persistent two-block scalar SCC column smoke\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
