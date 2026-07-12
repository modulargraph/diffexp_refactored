#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
      .as_object();
}

json::array prescriptions() {
  return json::array{json::object{
      {"factor_exact", "x"}, {"sign", -1}, {"multiplicity", 1},
      {"leading_coefficient_sign", 1}}};
}

std::string prepare_chart(const std::string& session,
                          const std::string& domain,
                          const std::string& suffix,
                          const std::string& center) {
  auto value = json::parse(R"json({
    "schema":2,"op":"chart.prepare","session":"placeholder",
    "key":"placeholder","identity":"placeholder",
    "analytic":{
      "geometry":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,
        "prescriptions":[]},
      "principal_matrix":[
        [{"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true}],
        [{"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true}]],
      "native_scc_capabilities":{"regular":true,"identity_gauge":true,
        "identity_v":true,"no_pseudo":true}},
    "scc":{"components":[[0],[1]],"structural_edges":[],
      "condensation_edges":[],"topological_order":[0,1],
      "coupling_depth":0},
    "problem":{"domain":"rational","precision_bits":256,
      "d":2,"fb":0,"w":3,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],
        "val":[null,null,null,null]}],
      "d0_inverse":"1","blocks":[[0],[1]],
      "assembly":{"identity":true,"poly":[],"rat":[],
        "val":[0,null,null,0]},"chop_digits":0}
  })json").as_object();
  value["session"] = session;
  value["key"] = "rational-row-" + domain + "-" + suffix;
  value["identity"] = "rational-row-chart-" + domain + "-" + suffix;
  value.at("problem").as_object()["domain"] = domain;
  auto& geometry = value.at("analytic").as_object()
      .at("geometry").as_object();
  geometry["center_exact"] = center;
  geometry["prescriptions"] = prescriptions();
  const auto response = request(std::move(value));
  if (response.at("status") != "ok")
    throw std::runtime_error("prepare: " + json::serialize(response));
  return std::string(response.at("chart").as_string());
}

json::object run() {
  json::array schedule_row{
      json::object{{"case", "R"}, {"da", "0"}, {"db", "0"}},
      json::object{{"case", "R"}, {"da", "0"}, {"db", "0"}}};
  json::array schedule;
  schedule.push_back(std::move(schedule_row));
  return json::object{
      {"nmax", 0}, {"p", 0}, {"has_initial", true},
      {"adaptive_probe", false}, {"a_target", "1/2"},
      {"b_target", "0"}, {"a_shift_min", 0},
      {"a_shifts", json::array{"1/2"}},
      {"schedule", std::move(schedule)},
      {"initial", json::array{"1", "0", "0", "2", "0", "0"}},
      {"initial_validity", json::array{2, 2}},
      {"source", nullptr}, {"return_u", false}};
}

json::object metadata(const std::string& checkpoint) {
  return json::object{
      {"chart", json::object{{"center_exact", "0"},
                              {"scale_exact", "1"}, {"radius", "2"},
                              {"infinite_radius", false}}},
      {"tag", json::object{
          {"a", json::object{{"domain", "rational"},
                              {"canonical", "1/2"}}},
          {"b", json::object{{"domain", "rational"},
                              {"canonical", "0"}}}}},
      {"prescriptions", prescriptions()},
      {"checkpoint_identity", checkpoint}};
}

json::object solve_local(const std::string& session,
                         const std::string& chart,
                         const std::string& checkpoint) {
  const auto response = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart}, {"run", run()},
      {"metadata", metadata(checkpoint)}});
  if (response.at("status") != "ok")
    throw std::runtime_error("solve: " + json::serialize(response));
  return response;
}

json::array zero_kernel() {
  return json::array{"0"};
}

json::object multiplier(std::int32_t shift, std::uint32_t pole,
                        json::array leading,
                        std::string identity) {
  return json::object{
      {"epsilon_shift", shift}, {"center_pole_order", pole},
      {"kernels", json::array{
           std::move(leading), zero_kernel(), zero_kernel()}},
      {"exact_identity", std::move(identity)},
      {"proven_zero", false}};
}

json::object rational_row() {
  return json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", 2}, {"exact_identity", "[1/t,3eps]"},
      {"entries", json::array{
           json::object{{"column", 0},
                        {"multiplier", multiplier(
                             0, 1, json::array{"1"},
                             "1/t")}},
           json::object{{"column", 1},
                        {"multiplier", multiplier(
                             1, 0, json::array{"3"},
                             "3eps")}}}}};
}

json::object zero_row() {
  return json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", 2}, {"exact_identity", "[0,0]"},
      {"entries", json::array{}}};
}

json::object apply_row(const std::string& session,
                       const std::string& local,
                       json::object row,
                       const std::string& source_checkpoint,
                       const std::string& checkpoint) {
  return request(json::object{
      {"schema", 2}, {"op", "local.apply_rational_row"},
      {"session", session}, {"local", local}, {"row", std::move(row)},
      {"source_checkpoint_identity", source_checkpoint},
      {"checkpoint_identity", checkpoint}});
}

json::object topology() {
  return json::object{
      {"singular_points", json::array{}},
      {"boundary_points", json::array{}},
      {"complex_projections", json::array{}},
      {"branch_sheets", json::array{
           json::object{{"factor_exact", "x"}, {"sign", -1}}}}};
}

json::object arm(const std::string& endpoint,
                 const std::string& anchor,
                 const std::string& endpoint_chart) {
  return json::object{
      {"from_exact", "0"}, {"to_exact", endpoint},
      {"charts", json::array{anchor, endpoint_chart}},
      {"topology", topology()}};
}

double coefficient_midpoint(const json::object& response,
                            std::size_t index) {
  const auto& encoded = response.at("value").as_object()
      .at("coefficients").as_array().at(index).as_array();
  return std::stod(std::string(encoded.at(0).as_string()));
}

bool rational_protocol() {
  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"},
      {"domain", "rational"}, {"output_digits", 50},
      {"local_capacity", 8}});
  const auto session = std::string(created.at("session").as_string());
  const auto anchor = prepare_chart(session, "rational", "anchor", "0");
  const auto lower = prepare_chart(session, "rational", "lower", "-2/3");
  const auto upper = prepare_chart(session, "rational", "upper", "2/3");
  const auto source = solve_local(session, anchor, "row-source-rational");
  const auto source_handle = std::string(source.at("local").as_string());

  auto malformed_row = rational_row();
  malformed_row.at("entries").as_array().front().as_object()
      .at("multiplier").as_object().at("kernels").as_array()
      .front() = json::array{"1", "0"};
  const auto malformed = apply_row(
      session, source_handle, std::move(malformed_row),
      "row-source-rational", "malformed-row");
  const auto projected = apply_row(
      session, source_handle, rational_row(), "row-source-rational",
      "row-result-rational");
  const auto projected_handle =
      std::string(projected.at("local").as_string());
  const auto zero = apply_row(
      session, source_handle, zero_row(), "row-source-rational",
      "zero-row-result");
  const auto zero_handle = std::string(zero.at("local").as_string());

  const auto plan = request(json::object{
      {"schema", 2}, {"op", "tile.plan"}, {"session", session},
      {"checkpoint_identity", "row-plan"}, {"division_order", 3},
      {"lower", arm("-2/3", anchor, lower)},
      {"upper", arm("2/3", anchor, upper)}});
  const auto plan_handle = std::string(plan.at("tile_plan").as_string());

  // The projected local strongly owns its source and remains fully usable
  // after the public vector-local token is released.
  (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                             {"session", session},
                             {"local", source_handle}});
  const auto positive = request(json::object{
      {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
      {"local", projected_handle},
      {"point", json::object{{"exact", "1/4"}}},
      {"options", json::object{{"tail_estimate", false}}}});
  const auto negative = request(json::object{
      {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
      {"local", projected_handle},
      {"point", json::object{{"exact", "-1/4"}}},
      {"options", json::object{{"tail_estimate", false}}}});
  const auto zero_value = request(json::object{
      {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
      {"local", zero_handle},
      {"point", json::object{{"exact", "1/4"}}},
      {"options", json::object{{"tail_estimate", false}}}});

  const auto line = request(json::object{
      {"schema", 2}, {"op", "integration.line"},
      {"session", session}, {"tile_plan", plan_handle},
      {"local", projected_handle}, {"arm", "upper"}, {"tile", 0},
      {"epsilon", json::object{{"min", 0}, {"max", 2}}},
      {"source_checkpoint_identity", "row-result-rational"},
      {"tile_plan_checkpoint_identity", "row-plan"},
      {"checkpoint_identity", "row-line"}});
  const auto line_handle = std::string(line.at("line").as_string());
  const auto exported = request(json::object{
      {"schema", 2}, {"op", "integration.export"},
      {"session", session}, {"line", line_handle},
      {"checkpoint_identity", "row-line"}, {"output_digits", 50}});
  const auto stats = request(json::object{
      {"schema", 2}, {"op", "session.stats"}, {"session", session}});

  const auto projected_json = json::serialize(projected);
  const auto expected_eps0 = 2.0 / std::sqrt(3.0);
  const auto expected_eps1 = 4.0 / (3.0 * std::sqrt(3.0));
  const bool ok =
      created.at("rational_row_application_capability") ==
          "retained-native-rational-row-local-application-v1" &&
      malformed.at("status") == "error" &&
      std::string(malformed.at("detail").as_string()).find(
          "exact source Taylor width") != std::string::npos &&
      projected.at("status") == "ok" &&
      projected.at("application_capability") ==
          "retained-native-rational-row-local-application-v1" &&
      projected.at("native_retained") == true &&
      projected.at("json_coefficients") == 0 &&
      projected.at("dimension") == 1 &&
      projected.at("epsilon_min") == 0 &&
      projected.at("epsilon_max") == 2 &&
      projected.at("checkpoint_identity") == "row-result-rational" &&
      projected.at("strong_derivation_ownership") == true &&
      projected_json.find("\"kernels\"") == std::string::npos &&
      projected_json.find("\"coefficients\"") == std::string::npos &&
      zero.at("dimension") == 1 &&
      zero.at("epsilon_min") == 0 && zero.at("epsilon_max") == 2 &&
      positive.at("status") == "ok" && positive.at("imaginary_sign").is_null() &&
      std::abs(coefficient_midpoint(positive, 0) - 2.0) < 1e-30 &&
      std::abs(coefficient_midpoint(positive, 1) - 3.0) < 1e-30 &&
      negative.at("status") == "ok" && negative.at("imaginary_sign") == -1 &&
      std::abs(coefficient_midpoint(zero_value, 0)) < 1e-30 &&
      line.at("status") == "ok" && line.at("dimension") == 1 &&
      line.at("json_coefficients") == 0 &&
      exported.at("json_coefficients") == 3 &&
      std::abs(coefficient_midpoint(exported, 0) - expected_eps0) < 1e-14 &&
      std::abs(coefficient_midpoint(exported, 1) - expected_eps1) < 1e-14 &&
      stats.at("locals") == 2 && stats.at("pending_local_solves") == 0;

  if (!ok) {
    std::cerr << "malformed: " << json::serialize(malformed) << '\n'
              << "projected: " << projected_json << '\n'
              << "positive: " << json::serialize(positive) << '\n'
              << "negative: " << json::serialize(negative) << '\n'
              << "zero: " << json::serialize(zero_value) << '\n'
              << "line: " << json::serialize(line) << '\n'
              << "exported: " << json::serialize(exported) << '\n'
              << "stats: " << json::serialize(stats) << '\n';
  }

  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  return ok;
}

bool acb_protocol() {
  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", "acb"},
      {"precision_bits", 256}, {"output_digits", 50},
      {"local_capacity", 4}});
  const auto session = std::string(created.at("session").as_string());
  const auto chart = prepare_chart(session, "acb", "anchor", "0");
  const auto source = solve_local(session, chart, "row-source-acb");
  const auto source_handle = std::string(source.at("local").as_string());
  json::array acb_leading;
  acb_leading.push_back(json::array{"[1 +/- 1e-40]", "0"});
  auto uncertain = multiplier(
      0, 0, std::move(acb_leading),
      "uncertain-unit");
  const auto row = json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", 2}, {"exact_identity", "[uncertain-one,0]"},
      {"entries", json::array{json::object{
           {"column", 0}, {"multiplier", std::move(uncertain)}}}}};
  const auto projected = apply_row(
      session, source_handle, row, "row-source-acb", "row-result-acb");
  if (projected.at("status") != "ok") {
    std::cerr << "Acb apply failed: " << json::serialize(projected) << '\n';
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    return false;
  }
  const auto projected_handle =
      std::string(projected.at("local").as_string());
  const auto evaluated = request(json::object{
      {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
      {"local", projected_handle},
      {"point", json::object{{"exact", "1/4"}}},
      {"options", json::object{{"tail_estimate", false}}}});
  const auto& first = evaluated.at("value").as_object()
      .at("coefficients").as_array().front().as_array();
  const bool ok =
      created.at("rational_row_application_capability") ==
          "retained-native-rational-row-local-application-v1" &&
      projected.at("status") == "ok" && projected.at("dimension") == 1 &&
      projected.at("json_coefficients") == 0 &&
      projected.at("metadata").as_object()
          .at("prescriptions").as_array().size() == 1 &&
      evaluated.at("status") == "ok" && evaluated.at("arithmetic_enclosed") == true &&
      first.size() == 4 && first.at(2) != "zero";
  if (!ok) {
    std::cerr << "Acb projected: " << json::serialize(projected) << '\n'
              << "Acb evaluated: " << json::serialize(evaluated) << '\n';
  }
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  return ok;
}

}  // namespace

int main() {
  const bool ok = rational_protocol() && acb_protocol();
  std::cout << (ok ? "PASS" : "FAIL")
            << ": retained Rational/Acb rational-row projection and tile integration\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
