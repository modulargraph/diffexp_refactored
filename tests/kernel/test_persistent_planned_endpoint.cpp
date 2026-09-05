#include "diffexp/kernel/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

#include <unistd.h>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp::kernel::run_recurrence_json(json::serialize(value)))
      .as_object();
}

json::array prescriptions(const std::string& factor,
                          std::optional<std::int32_t> sign) {
  if (!sign.has_value()) return {};
  return json::array{json::object{
      {"factor_exact", factor}, {"sign", *sign}, {"multiplicity", 1},
      {"leading_coefficient_sign", 1}}};
}

std::string prepare_chart(const std::string& session,
                          const std::string& domain,
                          const std::string& identity,
                          const std::string& center,
                          const std::string& scale,
                          const std::string& factor,
                          std::optional<std::int32_t> sign) {
  auto value = json::parse(R"json({
    "schema":2,"op":"chart.prepare","session":"placeholder",
    "key":"placeholder","identity":"placeholder",
    "analytic":{
      "geometry":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,
        "prescriptions":[]},
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
  value["key"] = identity;
  value["identity"] = identity;
  value.at("problem").as_object()["domain"] = domain;
  if (domain == "acb")
    value.at("problem").as_object()["precision_bits"] = 256;
  auto& geometry = value.at("analytic").as_object()
      .at("geometry").as_object();
  geometry["center_exact"] = center;
  geometry["scale_exact"] = scale;
  geometry["prescriptions"] = prescriptions(factor, sign);
  const auto response = request(std::move(value));
  if (response.at("status") != "ok")
    throw std::runtime_error("chart.prepare: " + json::serialize(response));
  return std::string(response.at("chart").as_string());
}

json::object solve_local(const std::string& session,
                         const std::string& chart,
                         const std::string& checkpoint,
                         const std::string& center,
                         const std::string& scale,
                         const std::string& factor,
                         std::optional<std::int32_t> sign) {
  json::array schedule_row{json::object{
      {"case", "R"}, {"da", "0"}, {"db", "0"}}};
  json::array schedule;
  schedule.push_back(std::move(schedule_row));
  const auto response = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart},
      {"run", json::object{
          {"nmax", 0}, {"p", 0}, {"has_initial", true},
          {"adaptive_probe", false}, {"a_target", "0"},
          {"b_target", "0"}, {"a_shift_min", 0},
          {"a_shifts", json::array{"0"}},
          {"schedule", std::move(schedule)},
          {"initial", json::array{"1", "0", "0"}},
          {"initial_validity", json::array{2}},
          {"source", nullptr}, {"return_u", false}}},
      {"metadata", json::object{
          {"chart", json::object{{"center_exact", center},
                                  {"scale_exact", scale},
                                  {"radius", "2"},
                                  {"infinite_radius", false}}},
          {"tag", json::object{
              {"a", json::object{{"domain", "rational"},
                                  {"canonical", "0"}}},
              {"b", json::object{{"domain", "rational"},
                                  {"canonical", "0"}}}}},
          {"prescriptions", prescriptions(factor, sign)},
          {"checkpoint_identity", checkpoint}}}});
  if (response.at("status") != "ok")
    throw std::runtime_error("local.solve: " + json::serialize(response));
  return response;
}

json::object topology(const std::string& factor,
                      std::optional<std::int32_t> sign) {
  json::array sheets;
  if (sign.has_value())
    sheets.push_back(json::object{{"factor_exact", factor},
                                  {"sign", *sign}});
  return json::object{{"singular_points", json::array{}},
                      {"boundary_points", json::array{}},
                      {"complex_projections", json::array{}},
                      {"branch_sheets", std::move(sheets)}};
}

json::object arm(const std::string& endpoint,
                 const std::string& anchor,
                 const std::string& final_chart,
                 const std::string& factor,
                 std::optional<std::int32_t> sign) {
  return json::object{{"from_exact", "0"}, {"to_exact", endpoint},
                      {"charts", json::array{anchor, final_chart}},
                      {"topology", topology(factor, sign)}};
}

json::object planned_endpoint(const std::string& session,
                              const json::object& plan,
                              const std::string& arm_name,
                              const json::object& local,
                              const std::string& checkpoint,
                              const std::string& cancellation_mode) {
  return request(json::object{
      {"schema", 2}, {"op", "tile.endpoint_limit"},
      {"session", session},
      {"tile_plan", plan.at("tile_plan")},
      {"tile_plan_checkpoint_identity", plan.at("checkpoint_identity")},
      {"arm", arm_name}, {"local", local.at("local")},
      {"source_checkpoint_identity", local.at("checkpoint_identity")},
      {"checkpoint_identity", checkpoint},
      {"cancellation", json::object{{"mode", cancellation_mode}}}});
}

bool run_domain(const std::string& domain, bool reversed_scale) {
  const auto checkpoint_path =
      (std::filesystem::temp_directory_path() /
       ("diffexp3-plan-endpoint-" + domain + "-" +
        std::to_string(::getpid()) + ".de2cp")).string();
  std::filesystem::remove(checkpoint_path);
  std::string session;
  std::string restored_session;
  bool ok = false;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"}, {"domain", domain},
        {"precision_bits", 256}, {"output_digits", 40},
        {"chart_capacity", 8}, {"local_capacity", 8},
        {"endpoint_capacity", 8}, {"tile_plan_capacity", 4}});
    session = std::string(created.at("session").as_string());

    const std::string scale = reversed_scale ? "-1" : "1";
    const std::string lower_factor = domain + "-lower-factor";
    const std::string upper_factor = domain + "-upper-factor";
    const std::optional<std::int32_t> lower_sign = 1;
    const std::optional<std::int32_t> upper_sign =
        domain == "rational" ? std::optional<std::int32_t>(-1)
                             : std::nullopt;
    const auto anchor = prepare_chart(
        session, domain, domain + "-endpoint-anchor", "0", scale,
        domain + "-anchor-factor", std::nullopt);
    const auto lower_chart = prepare_chart(
        session, domain, domain + "-endpoint-lower", "-2/3", scale,
        lower_factor, lower_sign);
    const auto upper_chart = prepare_chart(
        session, domain, domain + "-endpoint-upper", "2/3", scale,
        upper_factor, upper_sign);
    const auto lower_local = solve_local(
        session, lower_chart, domain + "-lower-local", "-2/3", scale,
        lower_factor, lower_sign);
    const auto upper_local = solve_local(
        session, upper_chart, domain + "-upper-local", "2/3", scale,
        upper_factor, upper_sign);
    const auto plan = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", domain + "-endpoint-plan"},
        {"division_order", 3},
        {"lower", arm("-2/3", anchor, lower_chart,
                       lower_factor, lower_sign)},
        {"upper", arm("2/3", anchor, upper_chart,
                       upper_factor, upper_sign)}});
    if (plan.at("status") != "ok")
      throw std::runtime_error("tile.plan: " + json::serialize(plan));

    const auto cancellation = domain == "rational"
        ? "exact-coefficient-field" : "exact-or-acb-singleton";
    const auto lower = planned_endpoint(
        session, plan, "lower", lower_local, domain + "-lower-endpoint",
        cancellation);
    const auto upper = planned_endpoint(
        session, plan, "upper", upper_local, domain + "-upper-endpoint",
        cancellation);

    auto mismatched_request = json::object{
        {"schema", 2}, {"op", "tile.endpoint_limit"},
        {"session", session}, {"tile_plan", plan.at("tile_plan")},
        {"tile_plan_checkpoint_identity", plan.at("checkpoint_identity")},
        {"arm", "lower"}, {"local", upper_local.at("local")},
        {"source_checkpoint_identity", upper_local.at("checkpoint_identity")},
        {"checkpoint_identity", domain + "-mismatch"},
        {"cancellation", json::object{{"mode", cancellation}}}};
    const auto mismatch = request(std::move(mismatched_request));
    auto stale_request = json::object{
        {"schema", 2}, {"op", "tile.endpoint_limit"},
        {"session", session}, {"tile_plan", plan.at("tile_plan")},
        {"tile_plan_checkpoint_identity", "stale-plan"},
        {"arm", "upper"}, {"local", upper_local.at("local")},
        {"source_checkpoint_identity", upper_local.at("checkpoint_identity")},
        {"checkpoint_identity", domain + "-stale"},
        {"cancellation", json::object{{"mode", cancellation}}}};
    const auto stale = request(std::move(stale_request));
    auto stale_local_request = json::object{
        {"schema", 2}, {"op", "tile.endpoint_limit"},
        {"session", session}, {"tile_plan", plan.at("tile_plan")},
        {"tile_plan_checkpoint_identity", plan.at("checkpoint_identity")},
        {"arm", "lower"}, {"local", lower_local.at("local")},
        {"source_checkpoint_identity", "stale-local"},
        {"checkpoint_identity", domain + "-stale-local"},
        {"cancellation", json::object{{"mode", cancellation}}}};
    const auto stale_local = request(std::move(stale_local_request));
    auto invalid_arm_request = json::object{
        {"schema", 2}, {"op", "tile.endpoint_limit"},
        {"session", session}, {"tile_plan", plan.at("tile_plan")},
        {"tile_plan_checkpoint_identity", plan.at("checkpoint_identity")},
        {"arm", "middle"}, {"local", lower_local.at("local")},
        {"source_checkpoint_identity", lower_local.at("checkpoint_identity")},
        {"checkpoint_identity", domain + "-invalid-arm"},
        {"cancellation", json::object{{"mode", cancellation}}}};
    const auto invalid_arm = request(std::move(invalid_arm_request));
    auto forbidden_direction = json::object{
        {"schema", 2}, {"op", "tile.endpoint_limit"},
        {"session", session}, {"tile_plan", plan.at("tile_plan")},
        {"tile_plan_checkpoint_identity", plan.at("checkpoint_identity")},
        {"arm", "upper"}, {"local", upper_local.at("local")},
        {"source_checkpoint_identity", upper_local.at("checkpoint_identity")},
        {"checkpoint_identity", domain + "-forbidden-direction"},
        {"cancellation", json::object{{"mode", cancellation}}},
        {"approach_direction", 1}};
    const auto caller_direction = request(std::move(forbidden_direction));

    const auto plan_handle = std::string(plan.at("tile_plan").as_string());
    const auto lower_handle =
        std::string(lower.at("endpoint").as_string());
    const auto upper_handle =
        std::string(upper.at("endpoint").as_string());
    (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                               {"session", session},
                               {"local", lower_local.at("local")}});
    (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                               {"session", session},
                               {"local", upper_local.at("local")}});
    (void)request(json::object{{"schema", 2}, {"op", "tile.release"},
                               {"session", session},
                               {"tile_plan", plan_handle}});
    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
        {"path", checkpoint_path},
        {"checkpoint_identity", domain + "-plan-endpoint-session"}});
    const auto before_close = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    session.clear();

    const auto restored = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint_path},
        {"expected_identity", domain + "-plan-endpoint-session"}});
    if (restored.at("status") == "ok")
      restored_session = std::string(restored.at("session").as_string());
    const auto lower_restored = restored.at("status") == "ok"
        ? request(json::object{{"schema", 2}, {"op", "endpoint.stats"},
                               {"session", restored_session},
                               {"endpoint", lower_handle}})
        : json::object{{"status", "error"}};
    const auto upper_restored = restored.at("status") == "ok"
        ? request(json::object{{"schema", 2}, {"op", "endpoint.stats"},
                               {"session", restored_session},
                               {"endpoint", upper_handle}})
        : json::object{{"status", "error"}};

    const auto expected_lower_approach = reversed_scale ? -1 : 1;
    const auto expected_upper_approach = reversed_scale ? 1 : -1;
    const auto expected_lower_rim = reversed_scale ? -1 : 1;
    const bool expected_upper_null = domain == "acb";
    const bool upper_rim_ok = expected_upper_null
        ? upper.at("effective_rim").is_null()
        : upper.at("effective_rim") == -1;
    const bool restored_upper_rim_ok = expected_upper_null
        ? upper_restored.at("effective_rim").is_null()
        : upper_restored.at("effective_rim") == -1;
    ok = created.at("planned_endpoint_limit_capability") ==
             "retained-native-plan-bound-endpoint-sector-limit-v1" &&
        lower.at("status") == "ok" && upper.at("status") == "ok" &&
        lower.at("execution_scope") == "plan-bound-final-arm-endpoint" &&
        lower.at("approach_direction") == expected_lower_approach &&
        upper.at("approach_direction") == expected_upper_approach &&
        lower.at("effective_rim") == expected_lower_rim && upper_rim_ok &&
        lower.at("requested_rim").is_null() &&
        upper.at("requested_rim").is_null() &&
        lower.at("source").as_object().at("endpoint_exact") == "-2/3" &&
        upper.at("source").as_object().at("endpoint_exact") == "2/3" &&
        mismatch.at("status") == "error" &&
        std::string(mismatch.at("detail").as_string()).find(
            "does not name the retained final chart") != std::string::npos &&
        stale.at("status") == "error" &&
        std::string(stale.at("detail").as_string()).find(
            "stale or mismatched") != std::string::npos &&
        stale_local.at("status") == "error" &&
        std::string(stale_local.at("detail").as_string()).find(
            "stale or mismatched") != std::string::npos &&
        invalid_arm.at("status") == "error" &&
        std::string(invalid_arm.at("detail").as_string()).find(
            "must be lower or upper") != std::string::npos &&
        caller_direction.at("status") == "error" &&
        std::string(caller_direction.at("detail").as_string()).find(
            "unknown or missing fields") != std::string::npos &&
        saved.at("status") == "ok" && saved.at("endpoints") == 2 &&
        before_close.at("locals") == 0 && before_close.at("tile_plans") == 0 &&
        restored.at("status") == "ok" &&
        restored.at("locals").as_array().empty() &&
        restored.at("tile_plans").as_array().empty() &&
        lower_restored.at("status") == "ok" &&
        upper_restored.at("status") == "ok" &&
        lower_restored.at("approach_direction") == expected_lower_approach &&
        upper_restored.at("approach_direction") == expected_upper_approach &&
        lower_restored.at("effective_rim") == expected_lower_rim &&
        restored_upper_rim_ok;

    if (!ok) {
      std::cerr << domain << " lower: " << json::serialize(lower) << '\n'
                << domain << " upper: " << json::serialize(upper) << '\n'
                << domain << " mismatch: " << json::serialize(mismatch) << '\n'
                << domain << " stale: " << json::serialize(stale) << '\n'
                << domain << " stale local: "
                << json::serialize(stale_local) << '\n'
                << domain << " invalid arm: "
                << json::serialize(invalid_arm) << '\n'
                << domain << " caller direction: "
                << json::serialize(caller_direction) << '\n'
                << domain << " save: " << json::serialize(saved) << '\n'
                << domain << " restore: " << json::serialize(restored) << '\n'
                << domain << " restored lower: "
                << json::serialize(lower_restored) << '\n'
                << domain << " restored upper: "
                << json::serialize(upper_restored) << '\n';
    }
  } catch (const std::exception& error) {
    std::cerr << domain << " FAIL: " << error.what() << '\n';
  }
  if (!session.empty())
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
  if (!restored_session.empty())
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", restored_session}});
  std::remove(checkpoint_path.c_str());
  return ok;
}

}  // namespace

int main() {
  const bool rational = run_domain("rational", false);
  const bool acb = run_domain("acb", true);
  const bool ok = rational && acb;
  std::cout << (ok ? "PASS" : "FAIL")
            << ": plan-bound Rational/Acb endpoint direction, rim, ownership, and checkpoint\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
