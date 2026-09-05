#include "diffexp/kernel/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include <unistd.h>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp::kernel::run_recurrence_json(json::serialize(value)))
      .as_object();
}

std::string prepare_chart(const std::string& session) {
  auto response = request(json::parse(std::string(R"json({
    "schema":2,"op":"chart.prepare","session":")json") + session +
      R"json(","key":"retained-endpoint-chart",
    "identity":"retained-endpoint-chart-v1",
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
  })json").as_object());
  if (response.at("status") != "ok")
    throw std::runtime_error(json::serialize(response));
  return std::string(response.at("chart").as_string());
}

std::string solve_regulated_local(const std::string& session,
                                  const std::string& chart) {
  auto response = request(json::parse(std::string(R"json({
    "schema":2,"op":"local.solve","session":")json") + session +
      R"json(","chart":")json" + chart + R"json(",
    "run":{"nmax":0,"p":0,"has_initial":true,
      "adaptive_probe":false,"a_target":"-3","b_target":"2",
      "a_shift_min":0,"a_shifts":["-3"],
      "schedule":[[{"case":"R","da":"0","db":"0"}]],
      "initial":["13","17","19"],"initial_validity":[2],
      "source":null,"return_u":false},
    "metadata":{"chart":{"center_exact":"0","scale_exact":"1",
        "radius":"2","infinite_radius":false},
      "tag":{"a":{"domain":"rational","canonical":"-3"},
        "b":{"domain":"rational","canonical":"2"}},
      "prescriptions":[{"factor_exact":"t","sign":-1,
        "multiplicity":1,"leading_coefficient_sign":1}],
      "checkpoint_identity":"regulated-local-v1"}
  })json").as_object());
  if (response.at("status") != "ok")
    throw std::runtime_error(json::serialize(response));
  return std::string(response.at("local").as_string());
}

}  // namespace

int main() {
  const std::string checkpoint_path =
      "/tmp/diffexp2_live_endpoint_checkpoint_" +
      std::to_string(static_cast<long long>(::getpid())) + ".de2cp";
  std::filesystem::remove(checkpoint_path);
  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", "rational"},
      {"output_digits", 30}, {"local_capacity", 1},
      {"endpoint_capacity", 1}});
  const auto session = std::string(created.at("session").as_string());
  const auto chart = prepare_chart(session);
  const auto local = solve_regulated_local(session, chart);

  const auto endpoint = request(json::object{
      {"schema", 2}, {"op", "local.endpoint_limit"},
      {"session", session}, {"local", local},
      {"source_checkpoint_identity", "regulated-local-v1"},
      {"checkpoint_identity", "regulated-endpoint-v1"},
      {"approach_direction", -1}, {"rim", -1},
      {"cancellation", json::object{
           {"mode", "exact-coefficient-field"}}}});
  const auto endpoint_handle =
      std::string(endpoint.at("endpoint").as_string());

  // The endpoint owns only its result and provenance.  Releasing the source
  // after successful publication must not affect inspection or export.
  (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                             {"session", session}, {"local", local}});
  const auto before_export = request(json::object{
      {"schema", 2}, {"op", "endpoint.stats"}, {"session", session},
      {"endpoint", endpoint_handle}});
  const auto exported = request(json::object{
      {"schema", 2}, {"op", "endpoint.export"}, {"session", session},
      {"endpoint", endpoint_handle},
      {"checkpoint_identity", "regulated-endpoint-v1"},
      {"output_digits", 30}});
  const auto after_export = request(json::object{
      {"schema", 2}, {"op", "endpoint.stats"}, {"session", session},
      {"endpoint", endpoint_handle}});
  const auto stats = request(json::object{
      {"schema", 2}, {"op", "session.stats"}, {"session", session}});
  const auto saved_checkpoint = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
      {"path", checkpoint_path},
      {"checkpoint_identity", "must-not-drop-live-endpoint"}});

  const auto& analytic = endpoint.at("analytic_regularization").as_object();
  const auto& metadata = analytic.at("metadata").as_object();
  const auto& sector = metadata.at("sectors").as_array().front().as_object();
  const auto& prescription =
      metadata.at("prescriptions").as_array().front().as_object();
  const auto& exported_value = exported.at("value").as_object();
  const auto& exported_zero =
      exported_value.at("coefficients").as_array().front().as_array();
  const bool ok = created.at("status") == "ok" &&
      created.at("endpoint_limit_capability") ==
          "retained-native-endpoint-sector-limit-v1" &&
      endpoint.at("status") == "ok" &&
      endpoint.at("native_retained") == true &&
      endpoint.at("json_coefficients") == 0 &&
      endpoint.at("effective_rim") == -1 &&
      endpoint.at("approach_direction") == -1 &&
      analytic.at("dropped_regulated_sectors") == 1 &&
      sector.at("b").as_object().at("canonical") == "2" &&
      prescription.at("sign") == -1 &&
      before_export.at("exports") == 0 &&
      before_export.if_contains("value") == nullptr &&
      exported.at("status") == "ok" &&
      exported.at("compatibility_export") == true &&
      exported.at("json_coefficients") == 1 &&
      exported_zero.size() >= 2 && exported_zero[0] == "0" &&
      exported_zero[1] == "0" && after_export.at("exports") == 1 &&
      stats.at("locals") == 0 && stats.at("endpoints") == 1 &&
      stats.at("endpoint_limits") == 1 &&
      stats.at("endpoint_exports") == 1 &&
      stats.at("pending_endpoint_limits") == 0 &&
      saved_checkpoint.at("status") == "ok" &&
      saved_checkpoint.at("endpoints") == 1 &&
      std::filesystem::exists(checkpoint_path) &&
      stats.at("endpoint_limit_capability") ==
          "retained-native-endpoint-sector-limit-v1";

  if (!ok) {
    std::cerr << "endpoint: " << json::serialize(endpoint) << '\n'
              << "before export: " << json::serialize(before_export) << '\n'
              << "export: " << json::serialize(exported) << '\n'
              << "after export: " << json::serialize(after_export) << '\n'
              << "checkpoint save: "
              << json::serialize(saved_checkpoint) << '\n'
              << "session stats: " << json::serialize(stats) << '\n';
  }

  (void)request(json::object{{"schema", 2}, {"op", "endpoint.release"},
                             {"session", session},
                             {"endpoint", endpoint_handle}});
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  std::filesystem::remove(checkpoint_path);
  std::cout << (ok ? "PASS" : "FAIL")
            << ": retained endpoint result protocol smoke\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
