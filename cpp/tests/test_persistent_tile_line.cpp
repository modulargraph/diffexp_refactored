#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cmath>
#include <cstdlib>
#include <future>
#include <iostream>
#include <stdexcept>
#include <string>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
      .as_object();
}

std::string prepare_chart(const std::string& session,
                          const std::string& key,
                          const std::string& identity,
                          const std::string& center) {
  auto value = json::parse(R"json({
    "schema":2,"op":"chart.prepare","session":"placeholder",
    "key":"placeholder","identity":"placeholder",
    "analytic":{
      "geometry":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,
        "prescriptions":[{"factor_exact":"t","sign":-1,
          "multiplicity":1,"leading_coefficient_sign":1}]},
      "principal_matrix":[[{"exact":"0","proven_zero":true}]],
      "native_scc_capabilities":{"regular":true,"identity_gauge":true,
        "identity_v":true,"no_pseudo":true}},
    "scc":{"components":[[0]],"structural_edges":[],
      "condensation_edges":[],"topological_order":[0],
      "coupling_depth":0},
    "problem":{"domain":"rational","d":1,"fb":0,"w":3,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],"val":[null]}],
      "d0_inverse":"1","blocks":[[0]],
      "assembly":{"identity":true,"poly":[],"rat":[],"val":[0]},
      "chop_digits":0}
  })json").as_object();
  value["session"] = session;
  value["key"] = key;
  value["identity"] = identity;
  value.at("analytic").as_object().at("geometry").as_object()
      ["center_exact"] = center;
  const auto response = request(std::move(value));
  if (response.at("status") != "ok")
    throw std::runtime_error(json::serialize(response));
  return std::string(response.at("chart").as_string());
}

json::object solve_anchor_local(const std::string& session,
                                const std::string& chart) {
  auto value = json::parse(R"json({
    "schema":2,"op":"local.solve","session":"placeholder",
    "chart":"placeholder",
    "run":{"nmax":0,"p":0,"has_initial":true,
      "adaptive_probe":false,"a_target":"0","b_target":"0",
      "a_shift_min":0,"a_shifts":["0"],
      "schedule":[[{"case":"R","da":"0","db":"0"}]],
      "initial":["1","0","0"],"initial_validity":[2],
      "source":null,"return_u":false},
    "metadata":{"chart":{"center_exact":"0","scale_exact":"1",
        "radius":"2","infinite_radius":false},
      "tag":{"a":{"domain":"rational","canonical":"0"},
        "b":{"domain":"rational","canonical":"0"}},
      "prescriptions":[{"factor_exact":"t","sign":-1,
        "multiplicity":1,"leading_coefficient_sign":1}],
      "checkpoint_identity":"anchor-local-v1"}
  })json").as_object();
  value["session"] = session;
  value["chart"] = chart;
  const auto response = request(std::move(value));
  if (response.at("status") != "ok")
    throw std::runtime_error(json::serialize(response));
  return response;
}

json::object topology() {
  return json::object{
      {"singular_points", json::array{}},
      {"boundary_points", json::array{}},
      {"complex_projections", json::array{}},
      {"branch_sheets", json::array{
           json::object{{"factor_exact", "t"}, {"sign", -1}}}}};
}

json::object arm(const std::string& endpoint,
                 const std::string& anchor,
                 const std::string& endpoint_chart) {
  return json::object{
      {"from_exact", "0"}, {"to_exact", endpoint},
      {"charts", json::array{anchor, endpoint_chart}},
      {"topology", topology()}};
}

double exported_midpoint(const json::object& response) {
  const auto& coefficients = response.at("value").as_object()
      .at("coefficients").as_array();
  return std::stod(std::string(
      coefficients.front().as_array().front().as_string()));
}

}  // namespace

int main() {
  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"},
      {"domain", "rational"}, {"output_digits", 40}});
  const auto session = std::string(created.at("session").as_string());
  const auto anchor = prepare_chart(
      session, "tile-anchor", "tile-anchor-v1", "0");
  const auto lower_chart = prepare_chart(
      session, "tile-lower", "tile-lower-v1", "-2/3");
  const auto upper_chart = prepare_chart(
      session, "tile-upper", "tile-upper-v1", "2/3");
  const auto local = solve_anchor_local(session, anchor);
  const auto local_handle = std::string(local.at("local").as_string());

  const auto plan = request(json::object{
      {"schema", 2}, {"op", "tile.plan"}, {"session", session},
      {"checkpoint_identity", "independent-arms-v1"},
      {"division_order", 3},
      {"lower", arm("-2/3", anchor, lower_chart)},
      {"upper", arm("2/3", anchor, upper_chart)}});
  const auto plan_handle = std::string(plan.at("tile_plan").as_string());
  const auto lower_match = request(json::object{
      {"schema", 2}, {"op", "tile.match_interval"},
      {"session", session}, {"tile_plan", plan_handle},
      {"arm", "lower"}, {"match", 0}});
  const auto upper_tile = request(json::object{
      {"schema", 2}, {"op", "tile.integration_interval"},
      {"session", session}, {"tile_plan", plan_handle},
      {"arm", "upper"}, {"tile", 0}});

  auto integrate = [&](const char* selected_arm, const char* checkpoint) {
    return request(json::object{
        {"schema", 2}, {"op", "integration.line"},
        {"session", session}, {"tile_plan", plan_handle},
        {"local", local_handle}, {"arm", selected_arm}, {"tile", 0},
        {"epsilon", json::object{{"min", 0}, {"max", 0}}},
        {"source_checkpoint_identity", "anchor-local-v1"},
        {"tile_plan_checkpoint_identity", "independent-arms-v1"},
        {"checkpoint_identity", checkpoint}});
  };
  auto lower_future = std::async(std::launch::async, integrate,
                                 "lower", "lower-line-v1");
  auto upper_future = std::async(std::launch::async, integrate,
                                 "upper", "upper-line-v1");
  const auto lower_line = lower_future.get();
  const auto upper_line = upper_future.get();
  const auto lower_line_handle =
      std::string(lower_line.at("line").as_string());
  const auto upper_line_handle =
      std::string(upper_line.at("line").as_string());

  // Published results own immutable plan and local snapshots.  Public tokens
  // may therefore be released without invalidating inspection or export.
  (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                             {"session", session}, {"local", local_handle}});
  (void)request(json::object{{"schema", 2}, {"op", "tile.release"},
                             {"session", session},
                             {"tile_plan", plan_handle}});
  const auto lower_stats = request(json::object{
      {"schema", 2}, {"op", "integration.stats"},
      {"session", session}, {"line", lower_line_handle}});
  const auto lower_export = request(json::object{
      {"schema", 2}, {"op", "integration.export"},
      {"session", session}, {"line", lower_line_handle},
      {"checkpoint_identity", "lower-line-v1"}, {"output_digits", 40}});
  const auto upper_export = request(json::object{
      {"schema", 2}, {"op", "integration.export"},
      {"session", session}, {"line", upper_line_handle},
      {"checkpoint_identity", "upper-line-v1"}, {"output_digits", 40}});
  const auto checkpoint_rejection = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
      {"path", "/tmp/diffexp2-live-line-result.de2cp"},
      {"checkpoint_identity", "must-not-drop-lines"}});
  const auto session_stats = request(json::object{
      {"schema", 2}, {"op", "session.stats"}, {"session", session}});

  const auto lower_value = exported_midpoint(lower_export);
  const auto upper_value = exported_midpoint(upper_export);
  const bool ok = created.at("tile_plan_capability") ==
          "retained-exact-independent-arm-tile-plan-v1" &&
      created.at("line_integration_capability") ==
          "retained-native-stored-truncation-physical-tile-integral-v1" &&
      plan.at("status") == "ok" && plan.at("independent_arms") == true &&
      plan.at("division_order") == 3 &&
      lower_match.at("physical_exact") == "-1/3" &&
      lower_match.at("producing_local_exact") == "-1/3" &&
      lower_match.at("receiving_local_exact") == "1/3" &&
      upper_tile.at("physical_begin_exact") == "0" &&
      upper_tile.at("physical_end_exact") == "1/3" &&
      upper_tile.at("local_begin_exact") == "0" &&
      upper_tile.at("local_end_exact") == "1/3" &&
      lower_line.at("status") == "ok" &&
      upper_line.at("status") == "ok" &&
      lower_line.at("scope") == "stored_truncation" &&
      lower_line.at("error").as_object().at("guarantee") == "none" &&
      lower_line.at("effective_rim") == -1 &&
      lower_stats.at("status") == "ok" &&
      lower_stats.at("exports") == 0 &&
      lower_export.at("json_coefficients") == 1 &&
      upper_export.at("json_coefficients") == 1 &&
      std::abs(lower_value + 1.0 / 3.0) < 1e-30 &&
      std::abs(upper_value - 1.0 / 3.0) < 1e-30 &&
      checkpoint_rejection.at("status") == "error" &&
      std::string(checkpoint_rejection.at("detail").as_string()).find(
          "line-result") != std::string::npos &&
      session_stats.at("locals") == 0 &&
      session_stats.at("tile_plans") == 0 &&
      session_stats.at("line_results") == 2 &&
      session_stats.at("line_integrations") == 2 &&
      session_stats.at("pending_line_integrations") == 0;

  if (!ok) {
    std::cerr << "plan: " << json::serialize(plan) << '\n'
              << "lower match: " << json::serialize(lower_match) << '\n'
              << "upper tile: " << json::serialize(upper_tile) << '\n'
              << "lower line: " << json::serialize(lower_line) << '\n'
              << "upper line: " << json::serialize(upper_line) << '\n'
              << "lower export: " << json::serialize(lower_export) << '\n'
              << "upper export: " << json::serialize(upper_export) << '\n'
              << "checkpoint: " << json::serialize(checkpoint_rejection)
              << '\n' << "session: " << json::serialize(session_stats)
              << '\n';
  }

  (void)request(json::object{{"schema", 2}, {"op", "integration.release"},
                             {"session", session},
                             {"line", lower_line_handle}});
  (void)request(json::object{{"schema", 2}, {"op", "integration.release"},
                             {"session", session},
                             {"line", upper_line_handle}});
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  std::cout << (ok ? "PASS" : "FAIL")
            << ": retained independent tile-plan and line-result protocol\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
