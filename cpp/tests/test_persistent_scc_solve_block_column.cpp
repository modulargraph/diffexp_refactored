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

json::object request(json::object value) {
  return request(json::serialize(value));
}

std::string prepare_scalar_chart(const std::string& session) {
  const auto response = request(std::string(R"json({
    "schema":2,"op":"chart.prepare","session":")json") + session +
    R"json(","key":"regular-block-scalar-key",
    "identity":"regular-block-scalar-v1",
    "analytic":{
      "geometry":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,"prescriptions":[]},
      "principal_matrix":[[{"exact":"0","proven_zero":true}]],
      "native_scc_capabilities":{"regular":true,"identity_gauge":true,
        "identity_v":true,"no_pseudo":true}},
    "scc":{"components":[[0]],"structural_edges":[],
      "condensation_edges":[],"topological_order":[0],"coupling_depth":0},
    "problem":{"domain":"rational","d":1,"fb":-1,"w":4,
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

std::string prepare_vector_chart(const std::string& session) {
  const auto response = request(std::string(R"json({
    "schema":2,"op":"chart.prepare","session":")json") + session +
    R"json(","key":"regular-block-vector-key",
    "identity":"regular-block-vector-v1",
    "analytic":{
      "geometry":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,"prescriptions":[]},
      "principal_matrix":[
        [{"exact":"0","proven_zero":true},
         {"exact":"1","proven_zero":false}],
        [{"exact":"1","proven_zero":false},
         {"exact":"0","proven_zero":true}]],
      "native_scc_capabilities":{"regular":true,"identity_gauge":true,
        "identity_v":true,"no_pseudo":true}},
    "scc":{"components":[[0,1]],"structural_edges":[[0,1],[1,0]],
      "condensation_edges":[],"topological_order":[0],"coupling_depth":0},
    "problem":{"domain":"rational","d":2,"fb":-1,"w":4,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[
        {"poly":[],"rat":[],"val":[null,null,null,null]},
        {"poly":[{"s":0,"e":[[0,1,"1"],[1,0,"1"]]}],
         "rat":[],"val":[null,0,0,null]}],
      "d0_inverse":"1","blocks":[[0],[1]],
      "assembly":{"identity":true,"poly":[],"rat":[],
        "val":[0,null,null,0]},"chop_digits":0}
  })json");
  if (response.at("status") != "ok")
    throw std::runtime_error(json::serialize(response));
  return std::string(response.at("chart").as_string());
}

json::object regular_run(std::uint32_t dimension, bool seed,
                         std::uint32_t seed_component = 0) {
  constexpr std::uint32_t nmax = 8;
  constexpr std::uint32_t frame_width = 4;
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    shifts.push_back(json::string(std::to_string(n)));
    json::array row;
    for (std::uint32_t component = 0; component < dimension; ++component) {
      row.push_back(json::object{{"case", n == 0 ? "R" : "T"},
                                 {"da", std::to_string(n)}, {"db", "0"}});
    }
    schedule.push_back(std::move(row));
  }
  json::array initial;
  json::array validity;
  for (std::uint32_t component = 0; component < dimension; ++component) {
    for (std::uint32_t epsilon = 0; epsilon < frame_width; ++epsilon)
      initial.push_back(seed && component == seed_component && epsilon == 1
                            ? "1" : "0");
    validity.push_back(seed ? json::value(2) : json::value(nullptr));
  }
  return json::object{{"nmax", nmax}, {"p", 0},
                      {"has_initial", seed}, {"adaptive_probe", false},
                      {"a_target", "0"}, {"b_target", "0"},
                      {"a_shift_min", 0}, {"a_shifts", std::move(shifts)},
                      {"schedule", std::move(schedule)},
                      {"initial", std::move(initial)},
                      {"initial_validity", std::move(validity)},
                      {"source", nullptr}, {"return_u", false}};
}

json::object metadata(const std::string& checkpoint) {
  return json::object{
      {"chart", json::object{{"center_exact", "0"}, {"scale_exact", "1"},
                              {"radius", "2"},
                              {"infinite_radius", false}}},
      {"tag", json::object{
          {"a", json::object{{"domain", "rational"}, {"canonical", "0"}}},
          {"b", json::object{{"domain", "rational"}, {"canonical", "0"}}}}},
      {"prescriptions", json::array{}},
      {"checkpoint_identity", checkpoint}};
}

double real_midpoint(const json::value& coefficient) {
  return std::stod(std::string(coefficient.as_array().front().as_string()));
}

bool close(double left, double right) {
  return std::abs(left - right) < 1e-20;
}

}  // namespace

int main() {
  const auto created = request(R"json({
    "schema":2,"op":"session.create","domain":"rational",
    "precision_bits":256,"output_digits":30,"scc_capacity":2,
    "local_capacity":2
  })json");
  const auto session = std::string(created.at("session").as_string());
  const auto scalar = prepare_scalar_chart(session);
  const auto vector = prepare_vector_chart(session);

  const auto prepared = request(std::string(R"json({
    "schema":2,"op":"scc.prepare","session":")json") + session +
    R"json(","key":"regular-block-dag-key",
    "identity":"regular-block-dag-parent-v2",
    "parent":{
      "dimension":3,
      "exact_system_record":[
        [{"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true}],
        [{"exact":"1/eps","proven_zero":false},
         {"exact":"0","proven_zero":true},
         {"exact":"1","proven_zero":false}],
        [{"exact":"0","proven_zero":true},
         {"exact":"1","proven_zero":false},
         {"exact":"0","proven_zero":true}]],
      "exact_theta_record":[
        [{"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true}],
        [{"exact":"t/eps","proven_zero":false},
         {"exact":"0","proven_zero":true},
         {"exact":"t","proven_zero":false}],
        [{"exact":"0","proven_zero":true},
         {"exact":"t","proven_zero":false},
         {"exact":"0","proven_zero":true}]],
      "chart":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,"prescriptions":[]},
      "scc":{"components":[[0],[1,2]],
        "structural_edges":[[0,1],[1,2],[2,1]],
        "condensation_edges":[[0,1]],"topological_order":[0,1],
        "coupling_depth":1},
      "execution":{"mode":"BlockSequentialStrict","work_t_order":8},
      "work_contract":{"work_min":-1,"requested_min":-1,
        "requested_max":1,"work_complete_max":2,"public_t_order":2,
        "wolfram_coupling_depth":2}},
    "blocks":[
      {"block":0,"vertices":[0],"chart":")json" + scalar +
    R"json(","principal_identity":"regular-block-scalar-v1",
       "regular":true,"identity_gauge":true,"identity_v":true,
       "no_pseudo":true},
      {"block":1,"vertices":[1,2],"chart":")json" + vector +
    R"json(","principal_identity":"regular-block-vector-v1",
       "regular":true,"identity_gauge":true,"identity_v":true,
       "no_pseudo":true}],
    "couplings":[{
      "source_block":0,"target_block":1,
      "source_vertices":[0],"target_vertices":[1,2],
      "rows":2,"columns":1,"exact_identity":"block-0-to-1-vector",
      "domain":"rational","symbols":[],
      "entries":[
        {"row":0,"column":0,"source_vertex":0,"target_vertex":1,
         "exact_original_entry":"1/eps","exact_theta_entry":"t/eps",
         "multiplier":{"epsilon_shift":-1,"center_pole_order":0,
           "kernels":[
             ["0","1","0","0","0","0","0","0","0"],
             ["0","0","0","0","0","0","0","0","0"],
             ["0","0","0","0","0","0","0","0","0"],
             ["0","0","0","0","0","0","0","0","0"]],
           "exact_identity":"t/eps","proven_zero":false}},
        {"row":1,"column":0,"source_vertex":0,"target_vertex":2,
         "exact_original_entry":"0","exact_theta_entry":"0",
         "multiplier":{"epsilon_shift":1,"center_pole_order":0,
           "kernels":[],"exact_identity":"0","proven_zero":true}}
      ]}]
  })json");
  if (prepared.at("status") != "ok") {
    std::cerr << "scc.prepare: " << json::serialize(prepared) << '\n';
    return EXIT_FAILURE;
  }
  const auto scc = std::string(prepared.at("scc").as_string());

  json::array targets;
  targets.push_back(json::object{{"block", 1},
                                 {"run", regular_run(2, false)},
                                 {"metadata", metadata("vector-target")}});
  const auto propagated = request(json::object{
      {"schema", 2}, {"op", "scc.solve_column"}, {"session", session},
      {"scc", scc}, {"checkpoint_identity", "regular-block-propagated"},
      {"seed", json::object{{"block", 0}, {"run", regular_run(1, true)},
                             {"metadata", metadata("scalar-seed")}}},
      {"targets", std::move(targets)}});

  json::object propagated_value;
  if (propagated.at("status") == "ok") {
    propagated_value = request(json::object{
        {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
        {"local", propagated.at("local")},
        {"point", json::object{{"exact", "1/2"}}},
        {"options", json::object{{"tail_estimate", false}}}});
  }

  const auto vector_seed = request(json::object{
      {"schema", 2}, {"op", "scc.solve_column"}, {"session", session},
      {"scc", scc}, {"checkpoint_identity", "regular-block-vector-seed"},
      {"seed", json::object{{"block", 1}, {"run", regular_run(2, true, 1)},
                             {"metadata", metadata("vector-seed")}}},
      {"targets", json::array{}}});

  json::object vector_seed_value;
  if (vector_seed.at("status") == "ok") {
    vector_seed_value = request(json::object{
        {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
        {"local", vector_seed.at("local")},
        {"point", json::object{{"exact", "1/2"}}},
        {"options", json::object{{"tail_estimate", false}}}});
  }
  const auto stats = request(json::object{
      {"schema", 2}, {"op", "scc.stats"}, {"session", session},
      {"scc", scc}});

  bool propagated_ok = false;
  if (propagated_value.if_contains("status") != nullptr &&
      propagated_value.at("status") == "ok") {
    const auto& coefficients = propagated_value.at("value").as_object()
                                   .at("coefficients").as_array();
    propagated_ok = coefficients.size() == 9 &&
        close(real_midpoint(coefficients[0]), 0.0) &&
        close(real_midpoint(coefficients[1]), 0.5) &&
        close(real_midpoint(coefficients[2]), 0.125) &&
        close(real_midpoint(coefficients[3]), 1.0) &&
        close(real_midpoint(coefficients[4]), 0.0) &&
        close(real_midpoint(coefficients[5]), 0.0);
  }

  bool vector_seed_ok = false;
  if (vector_seed_value.if_contains("status") != nullptr &&
      vector_seed_value.at("status") == "ok") {
    const auto& coefficients = vector_seed_value.at("value").as_object()
                                   .at("coefficients").as_array();
    const auto& provenance = vector_seed.at("column_provenance").as_object();
    vector_seed_ok = coefficients.size() == 6 &&
        close(real_midpoint(coefficients[0]), 0.0) &&
        close(real_midpoint(coefficients[1]), 0.5) &&
        close(real_midpoint(coefficients[2]), 1.125) &&
        provenance.at("seed_block") == 1 &&
        provenance.at("basis_index") == 2;
  }

  const bool ok = propagated.at("status") == "ok" && propagated_ok &&
      vector_seed.at("status") == "ok" && vector_seed_ok &&
      propagated.at("execution_capability") ==
          "exact-rational-regular-block-dag-column-v2" &&
      stats.at("execution_implemented") == true &&
      stats.at("execution_scope") ==
          "exact-rational-regular-block-dag-column-v2" &&
      stats.at("scalar_block_dag_column_execution") == false &&
      stats.at("regular_block_dag_column_execution") == true &&
      stats.at("active_coupling_entries") == 1 &&
      stats.at("proven_zero_coupling_entries") == 1 &&
      stats.at("min_coupling_shift") == -1 &&
      stats.at("max_coupling_shift") == -1 &&
      stats.at("scc_column_solves") == 2;

  if (!ok) {
    std::cerr << "propagated: " << json::serialize(propagated) << '\n'
              << "propagated value: " << json::serialize(propagated_value)
              << '\n' << "vector seed: " << json::serialize(vector_seed)
              << '\n' << "vector seed value: "
              << json::serialize(vector_seed_value) << '\n'
              << "scc.stats: " << json::serialize(stats) << '\n';
  }
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  std::cout << (ok ? "PASS" : "FAIL")
            << ": persistent 2+1-dimensional regular SCC block columns\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
