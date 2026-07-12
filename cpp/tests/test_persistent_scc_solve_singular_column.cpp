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

std::string prepare_singular_chart(const std::string& session,
                                   const std::string& key,
                                   const std::string& identity) {
  const auto response = request(std::string(R"json({
    "schema":2,"op":"chart.prepare","session":")json") + session +
    R"json(","key":")json" + key + R"json(","identity":")json" +
    identity + R"json(",
    "analytic":{
      "geometry":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,"prescriptions":[]},
      "principal_matrix":[[{"exact":"((1/2)+(eps/3))/t",
        "proven_zero":false}]],
      "native_scc_capabilities":{"regular":false,"identity_gauge":true,
        "identity_v":true,"no_pseudo":true}},
    "scc":{"components":[[0]],"structural_edges":[[0,0]],
      "condensation_edges":[],"topological_order":[0],"coupling_depth":0},
    "problem":{"domain":"rational","d":1,"fb":-3,"w":6,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[
          {"s":0,"e":[[0,0,"1/2"]]},
          {"s":1,"e":[[0,0,"1/3"]]}],
        "rat":[],"val":[0]}],
      "d0_inverse":"1","blocks":[[0]],
      "assembly":{"identity":true,"poly":[],"rat":[],"val":[0]},
      "chop_digits":0}
  })json");
  if (response.at("status") != "ok")
    throw std::runtime_error(json::serialize(response));
  return std::string(response.at("chart").as_string());
}

json::object sector_run(bool seed, std::uint32_t log_max) {
  constexpr std::uint32_t nmax = 6;
  constexpr std::uint32_t frame_width = 6;
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    const auto numerator = 1 + 2 * n;
    shifts.push_back(json::string(
        numerator == 1 ? "1/2" : std::to_string(numerator) + "/2"));
    schedule.push_back(json::array{json::object{
        {"case", n == 0 ? "R" : "T"},
        {"da", std::to_string(n)}, {"db", "0"}}});
  }
  json::array initial;
  json::array validity;
  for (std::uint32_t log = 0; log <= log_max; ++log) {
    for (std::uint32_t epsilon = 0; epsilon < frame_width; ++epsilon)
      initial.push_back(seed && log == 0 && epsilon == 3 ? "1" : "0");
    validity.push_back(seed && log == 0 ? json::value(2)
                                        : json::value(nullptr));
  }
  return json::object{{"nmax", nmax}, {"p", log_max},
      {"has_initial", seed}, {"adaptive_probe", false},
      {"a_target", "1/2"}, {"b_target", "1/3"},
      {"a_shift_min", 0}, {"a_shifts", std::move(shifts)},
      {"schedule", std::move(schedule)}, {"initial", std::move(initial)},
      {"initial_validity", std::move(validity)}, {"source", nullptr},
      {"return_u", false}};
}

json::object metadata(const std::string& checkpoint) {
  return json::object{
      {"chart", json::object{{"center_exact", "0"}, {"scale_exact", "1"},
                              {"radius", "2"},
                              {"infinite_radius", false}}},
      {"tag", json::object{
          {"a", json::object{{"domain", "rational"},
                              {"canonical", "1/2"}}},
          {"b", json::object{{"domain", "rational"},
                              {"canonical", "1/3"}}}}},
      {"prescriptions", json::array{}},
      {"checkpoint_identity", checkpoint}};
}

double real_midpoint(const json::value& coefficient) {
  return std::stod(
      std::string(coefficient.as_array().front().as_string()));
}

bool close(double left, double right) {
  return std::abs(left - right) < 1e-18;
}

json::object column_request(const std::string& session,
                            const std::string& scc,
                            std::uint32_t target_log_max,
                            const std::string& checkpoint) {
  json::array targets;
  targets.push_back(json::object{
      {"block", 1}, {"run", sector_run(false, target_log_max)},
      {"metadata", metadata(checkpoint + ":target")}});
  return json::object{
      {"schema", 2}, {"op", "scc.solve_column"}, {"session", session},
      {"scc", scc}, {"checkpoint_identity", checkpoint},
      {"seed", json::object{{"block", 0}, {"run", sector_run(true, 0)},
                              {"metadata", metadata(checkpoint + ":seed")}}},
      {"targets", std::move(targets)}};
}

}  // namespace

int main() {
  const auto created = request(R"json({
    "schema":2,"op":"session.create","domain":"rational",
    "precision_bits":256,"output_digits":30,"scc_capacity":2,
    "local_capacity":2
  })json");
  const auto session = std::string(created.at("session").as_string());
  const auto first = prepare_singular_chart(
      session, "singular-column-principal-0-key",
      "singular-column-principal-0");
  const auto second = prepare_singular_chart(
      session, "singular-column-principal-1-key",
      "singular-column-principal-1");

  const auto prepared = request(std::string(R"json({
    "schema":2,"op":"scc.prepare","session":")json") + session +
    R"json(","key":"singular-two-block-column-key",
    "identity":"singular-two-block-column-parent-v1",
    "parent":{
      "dimension":2,
      "exact_system_record":[
        [{"exact":"((1/2)+(eps/3))/t","proven_zero":false},
         {"exact":"0","proven_zero":true}],
        [{"exact":"1/(eps*t)","proven_zero":false},
         {"exact":"((1/2)+(eps/3))/t","proven_zero":false}]],
      "exact_theta_record":[
        [{"exact":"(1/2)+(eps/3)","proven_zero":false},
         {"exact":"0","proven_zero":true}],
        [{"exact":"1/eps","proven_zero":false},
         {"exact":"(1/2)+(eps/3)","proven_zero":false}]],
      "chart":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,"prescriptions":[]},
      "scc":{"components":[[0],[1]],
        "structural_edges":[[0,0],[0,1],[1,1]],
        "condensation_edges":[[0,1]],"topological_order":[0,1],
        "coupling_depth":1},
      "execution":{"mode":"BlockSequentialStrict","work_t_order":6},
      "work_contract":{"work_min":-3,"requested_min":-2,
        "requested_max":0,"work_complete_max":2,"public_t_order":0,
        "wolfram_coupling_depth":2}},
    "blocks":[
      {"block":0,"vertices":[0],"chart":")json" + first +
    R"json(","principal_identity":"singular-column-principal-0",
       "regular":false,"identity_gauge":true,"identity_v":true,
       "no_pseudo":true},
      {"block":1,"vertices":[1],"chart":")json" + second +
    R"json(","principal_identity":"singular-column-principal-1",
       "regular":false,"identity_gauge":true,"identity_v":true,
       "no_pseudo":true}],
    "couplings":[{
      "source_block":0,"target_block":1,
      "source_vertices":[0],"target_vertices":[1],
      "rows":1,"columns":1,"exact_identity":"singular-0-to-1",
      "domain":"rational","symbols":[],
      "entries":[{"row":0,"column":0,
        "source_vertex":0,"target_vertex":1,
        "exact_original_entry":"1/(eps*t)",
        "exact_theta_entry":"1/eps",
        "multiplier":{"epsilon_shift":-1,"center_pole_order":0,
          "kernels":[
            ["1","0","0","0","0","0","0"],
            ["0","0","0","0","0","0","0"],
            ["0","0","0","0","0","0","0"],
            ["0","0","0","0","0","0","0"],
            ["0","0","0","0","0","0","0"],
            ["0","0","0","0","0","0","0"]],
          "exact_identity":"1/eps","proven_zero":false}}]}]
  })json");
  if (prepared.at("status") != "ok") {
    std::cerr << "scc.prepare: " << json::serialize(prepared) << '\n';
    return EXIT_FAILURE;
  }
  const auto scc = std::string(prepared.at("scc").as_string());

  // p=0 cannot absorb the exact n=0 resonance.  This must be an E5, not a
  // truncated success or a fallback to another solver.
  const auto insufficient = request(column_request(
      session, scc, 0, "singular-column-insufficient-log"));
  const auto solved = request(column_request(
      session, scc, 1, "singular-column-checkpoint"));

  json::object evaluated;
  if (solved.at("status") == "ok") {
    evaluated = request(json::object{
        {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
        {"local", solved.at("local")},
        {"point", json::object{{"exact", "1/4"}}},
        {"options", json::object{{"tail_estimate", false}}}});
  }
  const auto stats = request(json::object{
      {"schema", 2}, {"op", "scc.stats"}, {"session", session},
      {"scc", scc}});

  bool value_ok = false;
  if (evaluated.if_contains("status") != nullptr &&
      evaluated.at("status") == "ok") {
    const auto& value = evaluated.at("value").as_object();
    const auto& coefficients = value.at("coefficients").as_array();
    const double logarithm = std::log(0.25);
    value_ok = value.at("min") == -1 &&
        value.at("max") == 0 && coefficients.size() == 4 &&
        close(real_midpoint(coefficients[0]), 0.0) &&
        close(real_midpoint(coefficients[1]), 0.5 * logarithm) &&
        close(real_midpoint(coefficients[2]), 0.5) &&
        close(real_midpoint(coefficients[3]),
              logarithm * logarithm / 6.0);
  }

  const bool ok = insufficient.at("status") == "error" &&
      insufficient.at("id") == "E5" &&
      solved.at("status") == "ok" &&
      solved.at("execution_capability") ==
          "exact-rational-regular-singular-scalar-block-dag-column-v1" &&
      solved.at("epsilon_min") == -2 && solved.at("epsilon_max") == 0 &&
      solved.at("dimension") == 2 && value_ok &&
      stats.at("execution_implemented") == true &&
      stats.at("execution_scope") ==
          "exact-rational-regular-singular-scalar-block-dag-column-v1" &&
      stats.at("regular_singular_scalar_block_dag_column_execution") == true &&
      stats.at("min_coupling_shift") == -1 &&
      stats.at("max_coupling_shift") == -1 &&
      stats.at("scc_column_solves") == 1;

  if (!ok) {
    std::cerr << "insufficient: " << json::serialize(insufficient) << '\n'
              << "solved: " << json::serialize(solved) << '\n'
              << "evaluated: " << json::serialize(evaluated) << '\n'
              << "stats: " << json::serialize(stats) << '\n';
  }
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  std::cout << (ok ? "PASS" : "FAIL")
            << ": persistent scalar regular-singular SCC resonance\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
