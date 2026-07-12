#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
      .as_object();
}

json::array nested(json::array inner) {
  json::array outer;
  outer.push_back(std::move(inner));
  return outer;
}

json::object session_create(const std::string& domain) {
  json::object value{{"schema", 2}, {"op", "session.create"},
                     {"domain", domain}, {"output_digits", 40}};
  if (domain == "acb") value["precision_bits"] = 256;
  return request(std::move(value));
}

json::object prepared_problem(const std::string& domain,
                              bool exponential) {
  json::array nhat{
      json::object{{"poly", json::array{}}, {"rat", json::array{}},
                   {"val", json::array{nullptr}}}};
  if (exponential)
    nhat.push_back(json::object{
        {"poly", json::array{json::object{
             {"s", 0},
             {"e", nested(json::array{0, 0, "1"})}}}},
        {"rat", json::array{}}, {"val", json::array{0}}});
  json::object result{
      {"domain", domain}, {"d", 1}, {"fb", 0}, {"w", 1},
      {"d_lags", nested(json::array{
           json::object{{"s", 0}, {"v", "1"}}})},
      {"denominators", json::array{}}, {"nhat_lags", std::move(nhat)},
      {"d0_inverse", "1"},
      {"blocks", nested(json::array{0})},
      {"assembly", json::object{
           {"identity", true}, {"poly", json::array{}},
           {"rat", json::array{}}, {"val", json::array{0}}}},
      {"chop_digits", 50}};
  if (domain == "acb") result["precision_bits"] = 256;
  return result;
}

std::string prepare_chart(const std::string& session,
                          const std::string& domain,
                          const std::string& key,
                          const std::string& identity,
                          const std::string& center,
                          bool exponential) {
  const auto principal = exponential
      ? json::object{{"exact", "t"}, {"proven_zero", false}}
      : json::object{{"exact", "0"}, {"proven_zero", true}};
  json::object chart_request{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", key}, {"identity", identity},
      {"analytic", json::object{
           {"geometry", json::object{
                {"center_exact", center}, {"scale_exact", "1"},
                {"radius_exact", "2"}, {"infinite_radius", false},
                {"prescriptions", json::array{}}}},
           {"principal_matrix", nested(json::array{principal})},
           {"native_scc_capabilities", json::object{
                {"regular", true}, {"identity_gauge", true},
                {"identity_v", true}, {"no_pseudo", true}}}}},
      {"scc", json::object{
           {"components", nested(json::array{0})},
           {"structural_edges", json::array{}},
           {"condensation_edges", json::array{}},
           {"topological_order", json::array{0}},
           {"coupling_depth", 0}}},
      {"problem", prepared_problem(domain, exponential)}};
  auto response = request(std::move(chart_request));
  if (response.at("status") != "ok")
    throw std::runtime_error(json::serialize(response));
  return std::string(response.at("chart").as_string());
}

json::object regular_run(std::uint32_t nmax, bool sourced) {
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    shifts.push_back(json::string(std::to_string(n)));
    schedule.push_back(json::array{json::object{
        {"case", n == 0 ? "R" : "T"},
        {"da", json::string(std::to_string(n))}, {"db", "0"}}});
  }
  json::value source = nullptr;
  if (sourced)
    source = json::object{{"frames", json::array{"0"}},
                          {"validity", json::array{0}},
                          {"present", json::array{true}}};
  return json::object{
      {"nmax", nmax}, {"p", 0}, {"has_initial", true},
      {"adaptive_probe", false}, {"a_target", "0"},
      {"b_target", "0"}, {"a_shift_min", 0},
      {"a_shifts", std::move(shifts)},
      {"schedule", std::move(schedule)},
      {"initial", json::array{"1"}},
      {"initial_validity", json::array{0}},
      {"source", std::move(source)}, {"return_u", false}};
}

json::object local_metadata(const std::string& checkpoint) {
  return json::object{
      {"chart", json::object{{"center_exact", "0"},
                              {"scale_exact", "1"}, {"radius", "2"},
                              {"infinite_radius", false}}},
      {"tag", json::object{
           {"a", json::object{{"domain", "rational"},
                               {"canonical", "0"}}},
           {"b", json::object{{"domain", "rational"},
                               {"canonical", "0"}}}}},
      {"prescriptions", json::array{}},
      {"checkpoint_identity", checkpoint}};
}

json::object solve_local(const std::string& session,
                         const std::string& chart,
                         std::uint32_t nmax, bool sourced,
                         const std::string& checkpoint) {
  return request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart}, {"run", regular_run(nmax, sourced)},
      {"metadata", local_metadata(checkpoint)}});
}

json::object evaluate_certified(const std::string& session,
                                const std::string& local,
                                const std::string& point) {
  return request(json::object{
      {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
      {"local", local}, {"point", json::object{{"exact", point}}},
      {"options", json::object{{"tail_estimate", false},
                                {"certified_tail_radius_exact", "1"}}}});
}

json::object empty_topology() {
  return json::object{{"singular_points", json::array{}},
                      {"boundary_points", json::array{}},
                      {"complex_projections", json::array{}},
                      {"branch_sheets", json::array{}}};
}

json::object arm(const std::string& endpoint,
                 const std::string& anchor,
                 const std::string& endpoint_chart) {
  return json::object{{"from_exact", "0"}, {"to_exact", endpoint},
                      {"charts", json::array{anchor, endpoint_chart}},
                      {"topology", empty_topology()}};
}

json::object integrate(const std::string& session,
                       const std::string& plan,
                       const std::string& local,
                       const std::string& source_checkpoint,
                       const std::string& result_checkpoint) {
  return request(json::object{
      {"schema", 2}, {"op", "integration.line"},
      {"session", session}, {"tile_plan", plan}, {"local", local},
      {"arm", "upper"}, {"tile", 0},
      {"epsilon", json::object{{"min", 0}, {"max", 0}}},
      {"source_checkpoint_identity", source_checkpoint},
      {"tile_plan_checkpoint_identity", "tail-plan-v1"},
      {"checkpoint_identity", result_checkpoint},
      {"certify_tail", true}});
}

bool run_rational_protocol() {
  const auto created = session_create("rational");
  const auto session = std::string(created.at("session").as_string());
  const auto anchor = prepare_chart(
      session, "rational", "tail-anchor", "tail-anchor-v1", "0", true);
  const auto lower = prepare_chart(
      session, "rational", "tail-lower", "tail-lower-v1", "-2/3", false);
  const auto upper = prepare_chart(
      session, "rational", "tail-upper", "tail-upper-v1", "2/3", false);
  const auto eligible = solve_local(
      session, anchor, 4, false, "tail-local-v1");
  const auto unsupported = solve_local(
      session, anchor, 0, true, "sourced-local-v1");
  const auto eligible_handle =
      std::string(eligible.at("local").as_string());
  const auto unsupported_handle =
      std::string(unsupported.at("local").as_string());
  const auto evaluated = evaluate_certified(
      session, eligible_handle, "1/2");
  const auto unsupported_evaluation = evaluate_certified(
      session, unsupported_handle, "1/2");
  const auto local_stats = request(json::object{
      {"schema", 2}, {"op", "local.stats"}, {"session", session},
      {"local", eligible_handle}});
  const auto plan = request(json::object{
      {"schema", 2}, {"op", "tile.plan"}, {"session", session},
      {"checkpoint_identity", "tail-plan-v1"}, {"division_order", 3},
      {"lower", arm("-2/3", anchor, lower)},
      {"upper", arm("2/3", anchor, upper)}});
  const auto plan_handle = std::string(plan.at("tile_plan").as_string());
  const auto certified_line = integrate(
      session, plan_handle, eligible_handle, "tail-local-v1",
      "certified-line-v1");
  const auto unsupported_line = integrate(
      session, plan_handle, unsupported_handle, "sourced-local-v1",
      "unsupported-line-v1");
  const auto certified_handle =
      std::string(certified_line.at("line").as_string());
  const auto exported = request(json::object{
      {"schema", 2}, {"op", "integration.export"},
      {"session", session}, {"line", certified_handle},
      {"checkpoint_identity", "certified-line-v1"},
      {"output_digits", 40}});

  const auto& value = evaluated.at("value").as_object();
  const auto& theta = evaluated.at("theta").as_object();
  const auto& certified_diagnostics =
      certified_line.at("diagnostics").as_object();
  const auto& unsupported_diagnostics =
      unsupported_line.at("diagnostics").as_object();
  const bool ok =
      created.at("certified_tail_capability") ==
          "retained-regular-homogeneous-gronwall-cauchy-tail-v1" &&
      eligible.at("tail_majorant").as_object().at("status") ==
          "certified" &&
      eligible.at("tail_majorant").as_object().at("attached") == true &&
      eligible.at("tail_majorant").as_object().at("operator_identity") ==
          "tail-anchor-v1" &&
      eligible.at("tail_majorant").as_object().at(
          "local_checkpoint_identity") == "tail-local-v1" &&
      eligible.at("tail_majorant").as_object().at(
          "checkpoint_serialized") == true &&
      unsupported.at("tail_majorant").as_object().at("status") ==
          "unsupported" &&
      evaluated.at("tail_certificate").as_object().at("status") ==
          "certified" &&
      value.at("error").as_object().at("guarantee") == "certified" &&
      theta.at("error").as_object().at("guarantee") == "certified" &&
      unsupported_evaluation.at("tail_certificate").as_object().at(
          "status") == "unsupported" &&
      unsupported_evaluation.at("value").as_object().if_contains("error") ==
          nullptr &&
      local_stats.at("tail_certificate_requests") == 1 &&
      local_stats.at("tail_certificate_certified") == 1 &&
      certified_line.at("scope") == "full_local_with_certified_tail" &&
      certified_line.at("capability") ==
          "retained-native-certified-full-local-physical-tile-integral-v1" &&
      certified_line.at("error").as_object().at("guarantee") ==
          "certified" &&
      std::string(certified_line.at("error").as_object().at(
          "provenance").as_string()).find("physical_jacobian_exact=1") !=
          std::string::npos &&
      certified_diagnostics.at("tail_certificate_requested") == true &&
      certified_diagnostics.at("tail_certificate_status") == "certified" &&
      !certified_diagnostics.at("tail_witness_radius_exact").is_null() &&
      unsupported_line.at("scope") == "stored_truncation" &&
      unsupported_line.at("error").as_object().at("guarantee") == "none" &&
      unsupported_diagnostics.at("tail_certificate_requested") == true &&
      unsupported_diagnostics.at("tail_certificate_status") ==
          "unsupported" &&
      exported.at("scope") == "full_local_with_certified_tail" &&
      exported.at("error_guarantee") == "certified" &&
      exported.at("value").as_object().at("error").as_object().at(
          "guarantee") == "certified";

  if (!ok) {
    std::cerr << "eligible: " << json::serialize(eligible) << '\n'
              << "unsupported: " << json::serialize(unsupported) << '\n'
              << "evaluate: " << json::serialize(evaluated) << '\n'
              << "unsupported evaluate: "
              << json::serialize(unsupported_evaluation) << '\n'
              << "local stats: " << json::serialize(local_stats) << '\n'
              << "certified line: " << json::serialize(certified_line)
              << '\n' << "unsupported line: "
              << json::serialize(unsupported_line) << '\n'
              << "export: " << json::serialize(exported) << '\n';
  }
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  return ok;
}

bool run_acb_attachment() {
  const auto created = session_create("acb");
  const auto session = std::string(created.at("session").as_string());
  const auto chart = prepare_chart(
      session, "acb", "acb-tail", "acb-tail-v1", "0", true);
  const auto local = solve_local(
      session, chart, 4, false, "acb-tail-local-v1");
  const auto evaluated = evaluate_certified(
      session, std::string(local.at("local").as_string()), "1/3");
  const bool ok = local.at("tail_majorant").as_object().at("status") ==
          "certified" &&
      local.at("tail_majorant").as_object().at("attached") == true &&
      evaluated.at("tail_certificate").as_object().at("status") ==
          "certified" &&
      evaluated.at("value").as_object().at("error").as_object().at(
          "guarantee") == "certified";
  if (!ok)
    std::cerr << "Acb local: " << json::serialize(local) << '\n'
              << "Acb evaluate: " << json::serialize(evaluated) << '\n';
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  return ok;
}

}  // namespace

int main() {
  const bool rational = run_rational_protocol();
  const bool acb = run_acb_attachment();
  const bool ok = rational && acb;
  std::cout << (ok ? "PASS" : "FAIL")
            << ": persistent certified regular tail protocol\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
