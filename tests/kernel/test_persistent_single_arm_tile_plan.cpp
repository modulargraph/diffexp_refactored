#include "diffexp/kernel/checkpoint.hpp"
#include "diffexp/kernel/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace json = boost::json;

namespace {

constexpr const char* kSingleCapability =
    "retained-exact-single-arm-tile-plan-v1";
constexpr const char* kTwoCapability =
    "retained-exact-independent-arm-tile-plan-v1";
constexpr const char* kSingleCheckpointSchema =
    "diffexp3-retained-single-arm-tile-plan-v1";

json::object request(json::object value) {
  return json::parse(diffexp::kernel::run_recurrence_json(json::serialize(value)))
      .as_object();
}

std::string prepare_chart(const std::string& session,
                          const std::string& key,
                          const std::string& center) {
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
    "problem":{"domain":"rational","d":1,"fb":0,"w":2,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],"val":[null]}],
      "d0_inverse":"1","blocks":[[0]],
      "assembly":{"identity":true,"poly":[],"rat":[],"val":[0]},
      "chop_digits":0}
  })json").as_object();
  value["session"] = session;
  value["key"] = key;
  value["identity"] = key + "-identity";
  value.at("analytic").as_object().at("geometry").as_object()
      ["center_exact"] = center;
  const auto response = request(std::move(value));
  if (response.at("status") != "ok")
    throw std::runtime_error("chart.prepare: " + json::serialize(response));
  return std::string(response.at("chart").as_string());
}

json::object topology() {
  return json::object{{"singular_points", json::array{}},
                      {"boundary_points", json::array{}},
                      {"complex_projections", json::array{}},
                      {"branch_sheets", json::array{}}};
}

json::object arm(const std::string& endpoint,
                 const std::string& anchor,
                 const std::string& receiver) {
  return json::object{{"from_exact", "0"}, {"to_exact", endpoint},
                      {"charts", json::array{anchor, receiver}},
                      {"topology", topology()}};
}

json::object degenerate_arm(const std::string& anchor) {
  return json::object{{"from_exact", "0"}, {"to_exact", "0"},
                      {"charts", json::array{anchor}},
                      {"topology", topology()}};
}

json::object plan_arm(const std::string& session,
                      const std::string& checkpoint,
                      json::object arm_request) {
  return request(json::object{
      {"schema", 2}, {"op", "tile.plan_arm"}, {"session", session},
      {"checkpoint_identity", checkpoint}, {"division_order", 3},
      {"arm", std::move(arm_request)}});
}

json::object interval(const std::string& session,
                      const std::string& plan,
                      const std::string& arm_name) {
  return request(json::object{
      {"schema", 2}, {"op", "tile.match_interval"},
      {"session", session}, {"tile_plan", plan},
      {"arm", arm_name}, {"match", 0}});
}

}  // namespace

int main() {
  const std::string path =
      "/tmp/diffexp3-persistent-single-arm-plan.de2cp";
  std::remove(path.c_str());
  std::string session;
  std::string restored_session;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "rational"}, {"output_digits", 30},
        {"chart_capacity", 3}, {"tile_plan_capacity", 4}});
    session = std::string(created.at("session").as_string());
    const auto anchor = prepare_chart(session, "single-plan-anchor", "0");
    const auto lower_chart = prepare_chart(
        session, "single-plan-lower", "-2/3");
    const auto upper_chart = prepare_chart(
        session, "single-plan-upper", "2/3");

    const auto lower = plan_arm(
        session, "single-plan-lower-checkpoint",
        arm("-2/3", anchor, lower_chart));
    const auto upper = plan_arm(
        session, "single-plan-upper-checkpoint",
        arm("2/3", anchor, upper_chart));
    const auto two = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "two-plan-checkpoint"},
        {"division_order", 3},
        {"lower", arm("-2/3", anchor, lower_chart)},
        {"upper", arm("2/3", anchor, upper_chart)}});
    if (lower.at("status") != "ok" || upper.at("status") != "ok" ||
        two.at("status") != "ok")
      throw std::runtime_error(
          "plan creation failed: " + json::serialize(lower) + " / " +
          json::serialize(upper) + " / " + json::serialize(two));
    const auto lower_handle = std::string(lower.at("tile_plan").as_string());
    const auto upper_handle = std::string(upper.at("tile_plan").as_string());
    const auto two_handle = std::string(two.at("tile_plan").as_string());

    const auto lower_match = interval(session, lower_handle, "lower");
    const auto upper_match = interval(session, upper_handle, "upper");
    const auto two_lower_match = interval(session, two_handle, "lower");
    const auto two_upper_match = interval(session, two_handle, "upper");
    const auto absent_upper = interval(session, lower_handle, "upper");
    const auto absent_lower = interval(session, upper_handle, "lower");
    if (lower_match.at("status") != "ok" ||
        upper_match.at("status") != "ok" ||
        absent_upper.at("status") != "error" ||
        absent_lower.at("status") != "error")
      throw std::runtime_error(
          "single-arm lookup contract failed: " +
          json::serialize(lower_match) + " / " +
          json::serialize(upper_match) + " / " +
          json::serialize(absent_upper) + " / " +
          json::serialize(absent_lower));

    const auto before_malformed = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto malformed = plan_arm(
        session, "single-plan-malformed", degenerate_arm(anchor));
    const auto after_malformed = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const bool malformed_atomic = malformed.at("status") == "error" &&
        before_malformed.at("tile_plans") ==
            after_malformed.at("tile_plans") &&
        before_malformed.at("tile_plans_created") ==
            after_malformed.at("tile_plans_created") &&
        after_malformed.at("pending_tile_plans") == 0;
    if (!malformed_atomic)
      throw std::runtime_error(
          "malformed single arm published state: " +
          json::serialize(malformed) + " / " +
          json::serialize(after_malformed));

    const auto released_two = request(json::object{
        {"schema", 2}, {"op", "tile.release"}, {"session", session},
        {"tile_plan", two_handle}});
    if (released_two.at("status") != "ok")
      throw std::runtime_error(
          "two-arm tile.release: " + json::serialize(released_two));

    for (const auto& chart : {anchor, lower_chart, upper_chart}) {
      const auto released = request(json::object{
          {"schema", 2}, {"op", "chart.release"}, {"session", session},
          {"chart", chart}});
      if (released.at("status") != "ok")
        throw std::runtime_error(
            "chart.release: " + json::serialize(released));
    }
    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
        {"path", path}, {"checkpoint_identity", "single-plan-roundtrip"}});
    if (saved.at("status") != "ok")
      throw std::runtime_error("checkpoint.save: " +
                               json::serialize(saved));
    const auto payload = json::parse(
        diffexp::kernel::checkpoint::read(path).payload_json).as_object();
    std::size_t single_records = 0;
    for (const auto& raw_plan :
         payload.at("retained_tile_plans").as_array()) {
      const auto& record = raw_plan.as_object();
      const auto schema = std::string(record.at("schema").as_string());
      if (schema == kSingleCheckpointSchema) {
        ++single_records;
        if (record.at("arm_name") != "lower" &&
            record.at("arm_name") != "upper")
          throw std::runtime_error(
              "single checkpoint lost its arm name");
      } else {
        throw std::runtime_error(
            "unexpected tile checkpoint schema: " + schema);
      }
    }
    const auto& visibility = payload.at("session").as_object()
                                 .at("registry_visibility").as_object();
    const bool closure_ok =
        payload.at("prepared_charts").as_array().size() == 3 &&
        visibility.at("charts").as_array().empty() &&
        visibility.at("tile_plans").as_array().size() == 2 &&
        single_records == 2;
    if (!closure_ok)
      throw std::runtime_error(
          "checkpoint did not retain the hidden chart closure: " +
          json::serialize(payload));

    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    session.clear();
    const auto restored = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"}, {"path", path},
        {"expected_identity", "single-plan-roundtrip"}});
    if (restored.at("status") == "ok")
      restored_session = std::string(restored.at("session").as_string());
    if (restored.at("status") != "ok")
      throw std::runtime_error("checkpoint.restore: " +
                               json::serialize(restored));

    const auto restored_lower = request(json::object{
        {"schema", 2}, {"op", "tile.stats"},
        {"session", restored_session}, {"tile_plan", lower_handle}});
    const auto restored_upper = request(json::object{
        {"schema", 2}, {"op", "tile.stats"},
        {"session", restored_session}, {"tile_plan", upper_handle}});
    const auto restored_lower_match = interval(
        restored_session, lower_handle, "lower");
    const auto restored_upper_match = interval(
        restored_session, upper_handle, "upper");

    const bool ok =
        created.at("single_arm_tile_plan_capability") == kSingleCapability &&
        lower.at("capability") == kSingleCapability &&
        upper.at("capability") == kSingleCapability &&
        lower.at("arm_name") == "lower" &&
        upper.at("arm_name") == "upper" &&
        lower.at("matches") == 1 && lower.at("tiles") == 2 &&
        upper.at("matches") == 1 && upper.at("tiles") == 2 &&
        lower_match.at("physical_exact") == "-1/3" &&
        upper_match.at("physical_exact") == "1/3" &&
        two.at("capability") == kTwoCapability &&
        two.at("independent_arms") == true &&
        two.at("lower_matches") == 1 && two.at("upper_matches") == 1 &&
        two_lower_match.at("status") == "ok" &&
        two_lower_match.at("physical_exact") == "-1/3" &&
        two_upper_match.at("status") == "ok" &&
        two_upper_match.at("physical_exact") == "1/3" &&
        released_two.at("checkpoint_identity") == "two-plan-checkpoint" &&
        restored.at("charts").as_array().empty() &&
        restored.at("tile_plans").as_array().size() == 2 &&
        restored_lower.at("status") == "ok" &&
        restored_lower.at("capability") == kSingleCapability &&
        restored_lower.at("arm_name") == "lower" &&
        restored_lower.at("match_interval_queries") == 2 &&
        restored_upper.at("status") == "ok" &&
        restored_upper.at("capability") == kSingleCapability &&
        restored_upper.at("arm_name") == "upper" &&
        restored_upper.at("match_interval_queries") == 2 &&
        restored_lower_match.at("physical_exact") == "-1/3" &&
        restored_upper_match.at("physical_exact") == "1/3";
    if (!ok)
      std::cerr << "lower: " << json::serialize(lower) << '\n'
                << "upper: " << json::serialize(upper) << '\n'
                << "two: " << json::serialize(two) << '\n'
                << "restored: " << json::serialize(restored) << '\n'
                << "restored lower: " << json::serialize(restored_lower)
                << '\n' << "restored upper: "
                << json::serialize(restored_upper) << '\n';

    for (const auto& plan : {lower_handle, upper_handle}) {
      const auto released = request(json::object{
          {"schema", 2}, {"op", "tile.release"},
          {"session", restored_session}, {"tile_plan", plan}});
      if (released.at("status") != "ok")
        throw std::runtime_error(
            "tile.release: " + json::serialize(released));
    }
    const auto final_stats = request(json::object{
        {"schema", 2}, {"op", "session.stats"},
        {"session", restored_session}});
    if (final_stats.at("tile_plans") != 0)
      throw std::runtime_error(
          "tile.release retained a public plan: " +
          json::serialize(final_stats));
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", restored_session}});
    restored_session.clear();
    std::remove(path.c_str());
    std::cout << (ok ? "PASS" : "FAIL")
              << ": retained one-or-two-arm tile-plan foundation\n";
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    if (!restored_session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored_session}});
    std::remove(path.c_str());
    std::cerr << "FAIL: retained single-arm tile plan: "
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
