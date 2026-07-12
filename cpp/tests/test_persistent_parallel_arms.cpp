#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
      .as_object();
}

json::array prescriptions() {
  return json::array{json::object{
      {"factor_exact", "parallel-arm-f"}, {"sign", -1},
      {"multiplicity", 1}, {"leading_coefficient_sign", 1}}};
}

std::string prepare_chart(const std::string& session,
                          const std::string& name,
                          const std::string& center) {
  json::array principal_row;
  principal_row.push_back(
      json::object{{"exact", "0"}, {"proven_zero", true}});
  json::array principal_matrix;
  principal_matrix.push_back(std::move(principal_row));
  json::array component;
  component.emplace_back(0);
  json::array components;
  components.push_back(component);
  json::array d_lag;
  d_lag.push_back(json::object{{"s", 0}, {"v", "1"}});
  json::array d_lags;
  d_lags.push_back(std::move(d_lag));
  json::array blocks;
  blocks.push_back(std::move(component));
  const auto response = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", name}, {"identity", name + "-identity"},
      {"analytic", json::object{
           {"geometry", json::object{
                {"center_exact", center}, {"scale_exact", "1"},
                {"radius_exact", "2"}, {"infinite_radius", false},
                {"prescriptions", prescriptions()}}},
           {"principal_matrix", std::move(principal_matrix)},
           {"native_scc_capabilities", json::object{
                {"regular", true}, {"identity_gauge", true},
                {"identity_v", true}, {"no_pseudo", true}}}}},
      {"scc", json::object{
           {"components", std::move(components)},
           {"structural_edges", json::array{}},
           {"condensation_edges", json::array{}},
           {"topological_order", json::array{0}},
           {"coupling_depth", 0}}},
      {"problem", json::object{
           {"domain", "rational"}, {"d", 1}, {"fb", 0}, {"w", 3},
           {"d_lags", std::move(d_lags)},
           {"denominators", json::array{}},
           {"nhat_lags", json::array{json::object{
                {"poly", json::array{}}, {"rat", json::array{}},
                {"val", json::array{nullptr}}}}},
           {"d0_inverse", "1"},
           {"blocks", std::move(blocks)},
           {"assembly", json::object{
                {"identity", true}, {"poly", json::array{}},
                {"rat", json::array{}}, {"val", json::array{0}}}},
           {"chop_digits", 0}}}});
  if (response.at("status") != "ok")
    throw std::runtime_error("prepare: " + json::serialize(response));
  return std::string(response.at("chart").as_string());
}

std::string solve_local(const std::string& session,
                        const std::string& chart,
                        const std::string& center,
                        const std::string& checkpoint,
                        const std::string& value) {
  json::array schedule_row;
  schedule_row.push_back(json::object{
      {"case", "R"}, {"da", "0"}, {"db", "0"}});
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
           {"initial", json::array{value, "0", "0"}},
           {"initial_validity", json::array{2}}, {"source", nullptr},
           {"return_u", false}}},
      {"metadata", json::object{
           {"chart", json::object{
                {"center_exact", center}, {"scale_exact", "1"},
                {"radius", "2"}, {"infinite_radius", false}}},
           {"tag", json::object{
                {"a", json::object{{"domain", "rational"},
                                     {"canonical", "0"}}},
                {"b", json::object{{"domain", "rational"},
                                     {"canonical", "0"}}}}},
           {"prescriptions", prescriptions()},
           {"checkpoint_identity", checkpoint}}}});
  if (response.at("status") != "ok")
    throw std::runtime_error("solve: " + json::serialize(response));
  return std::string(response.at("local").as_string());
}

json::object topology() {
  return json::object{
      {"singular_points", json::array{}},
      {"boundary_points", json::array{}},
      {"complex_projections", json::array{}},
      {"branch_sheets", json::array{json::object{
           {"factor_exact", "parallel-arm-f"}, {"sign", -1}}}}};
}

json::object arm(const std::string& endpoint,
                 const std::vector<std::string>& charts) {
  json::array handles;
  for (const auto& chart : charts) handles.emplace_back(chart);
  return json::object{{"from_exact", "0"}, {"to_exact", endpoint},
                      {"charts", std::move(handles)},
                      {"topology", topology()}};
}

json::object arm_execution(const std::string& first,
                           const std::string& second) {
  return json::object{
      {"receiving_basis",
       json::array{json::array{first}, json::array{second}}},
      {"match_policies", json::array{json::object{}, json::object{}}}};
}

json::object run_arms_request(const std::string& session,
                              const std::string& plan,
                              const std::string& anchor,
                              json::object lower, json::object upper,
                              const std::string& root) {
  return json::object{
      {"schema", 2}, {"op", "integration.run_arms"},
      {"session", session}, {"tile_plan", plan}, {"anchor", anchor},
      {"tile_plan_checkpoint_identity", "parallel-plan"},
      {"anchor_checkpoint_identity", "parallel-anchor"},
      {"epsilon", json::object{{"min", 0}, {"max", 2},
                                {"required_complete_max", 2}}},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
           {"root", root}}},
      {"lower", std::move(lower)}, {"upper", std::move(upper)}};
}

double exported_coefficient(const json::object& response,
                            std::size_t epsilon_index = 0) {
  const auto& raw = response.at("value").as_object()
                        .at("coefficients").as_array().at(epsilon_index);
  return raw.is_string()
      ? std::stod(std::string(raw.as_string()))
      : std::stod(std::string(raw.as_array().front().as_string()));
}

}  // namespace

int main() {
  const std::string checkpoint_path =
      "/tmp/diffexp2-persistent-parallel-arms.de2cp";
  std::remove(checkpoint_path.c_str());
  std::string session;
  std::string restored_session;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "rational"}, {"output_digits", 40},
        {"chart_capacity", 8}, {"local_capacity", 16},
        {"match_capacity", 8}, {"tile_plan_capacity", 2}});
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_chart(session, "pair-anchor", "0");
    const auto lower_1_chart = prepare_chart(session, "pair-lower-1", "-2/3");
    const auto lower_2_chart = prepare_chart(session, "pair-lower-2", "-4/3");
    const auto upper_1_chart = prepare_chart(session, "pair-upper-1", "2/3");
    const auto upper_2_chart = prepare_chart(session, "pair-upper-2", "4/3");
    const auto anchor = solve_local(session, anchor_chart, "0",
                                    "parallel-anchor", "2");
    const auto lower_1 = solve_local(session, lower_1_chart, "-2/3",
                                     "parallel-lower-basis-1", "1");
    const auto lower_2 = solve_local(session, lower_2_chart, "-4/3",
                                     "parallel-lower-basis-2", "1");
    const auto upper_1 = solve_local(session, upper_1_chart, "2/3",
                                     "parallel-upper-basis-1", "1");
    const auto upper_2 = solve_local(session, upper_2_chart, "4/3",
                                     "parallel-upper-basis-2", "1");
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "parallel-plan"}, {"division_order", 3},
        {"lower", arm("-4/3", {anchor_chart, lower_1_chart, lower_2_chart})},
        {"upper", arm("4/3", {anchor_chart, upper_1_chart, upper_2_chart})}});
    if (planned.at("status") != "ok")
      throw std::runtime_error("plan: " + json::serialize(planned));
    const auto plan = std::string(planned.at("tile_plan").as_string());

    const auto before_failure = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto failed = request(run_arms_request(
        session, plan, anchor, arm_execution(lower_1, lower_2),
        // The first upper basis is deliberately bound to the lower chart.
        // Admission succeeds, one worker may complete, but publication must
        // still be all-or-nothing after the upper plan-binding failure.
        arm_execution(lower_1, upper_2), "parallel-failure"));
    const auto after_failure = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    if (failed.at("status") != "error" ||
        before_failure.at("locals") != after_failure.at("locals") ||
        before_failure.at("matches") != after_failure.at("matches") ||
        before_failure.at("line_results") !=
            after_failure.at("line_results") ||
        after_failure.at("pending_matches") != 0 ||
        after_failure.at("pending_local_solves") != 0 ||
        after_failure.at("pending_line_integrations") != 0)
      throw std::runtime_error(
          "failed arm pair published partial state: " +
          json::serialize(failed) + " / " + json::serialize(after_failure));

    const auto marched = request(run_arms_request(
        session, plan, anchor, arm_execution(lower_1, lower_2),
        arm_execution(upper_1, upper_2), "parallel-success"));
    if (marched.at("status") != "ok")
      throw std::runtime_error("march: " + json::serialize(marched));
    const auto& arms = marched.at("arms").as_object();
    const auto& lower = arms.at("lower").as_object();
    const auto& upper = arms.at("upper").as_object();
    const auto& lower_local = lower.at("final_local").as_object();
    const auto& upper_local = upper.at("final_local").as_object();
    const auto& lower_line = lower.at("line_result").as_object();
    const auto& upper_line = upper.at("line_result").as_object();
    const auto& combined_line =
        marched.at("combined_line_result").as_object();

    // Every public source token can disappear after the single call.  The
    // returned locals and aggregates retain the complete owner chains.
    for (const auto& local : {anchor, lower_1, lower_2, upper_1, upper_2})
      (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                                 {"session", session}, {"local", local}});
    (void)request(json::object{{"schema", 2}, {"op", "tile.release"},
                               {"session", session}, {"tile_plan", plan}});
    for (const auto& chart : {anchor_chart, lower_1_chart, lower_2_chart,
                              upper_1_chart, upper_2_chart})
      (void)request(json::object{{"schema", 2}, {"op", "chart.release"},
                                 {"session", session}, {"chart", chart}});

    const auto lower_evaluation = request(json::object{
        {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
        {"local", lower_local.at("local")},
        {"point", json::object{{"exact", "0"}}},
        {"options", json::object{{"tail_estimate", false}}}});
    const auto upper_evaluation = request(json::object{
        {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
        {"local", upper_local.at("local")},
        {"point", json::object{{"exact", "0"}}},
        {"options", json::object{{"tail_estimate", false}}}});
    const auto export_line = [&](const json::object& line) {
      return request(json::object{
          {"schema", 2}, {"op", "integration.export"},
          {"session", session}, {"line", line.at("line")},
          {"checkpoint_identity", line.at("checkpoint_identity")},
          {"output_digits", 40}});
    };
    const auto lower_export = export_line(lower_line);
    const auto upper_export = export_line(upper_line);
    const auto combined_export = export_line(combined_line);

    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
        {"path", checkpoint_path},
        {"checkpoint_identity", "parallel-arm-roundtrip"}});
    if (saved.at("status") != "ok")
      throw std::runtime_error("checkpoint save: " + json::serialize(saved));
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    session.clear();
    const auto restored = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint_path},
        {"expected_identity", "parallel-arm-roundtrip"}});
    if (restored.at("status") != "ok")
      throw std::runtime_error(
          "checkpoint restore: " + json::serialize(restored));
    restored_session = std::string(restored.at("session").as_string());
    const auto restored_combined = request(json::object{
        {"schema", 2}, {"op", "integration.export"},
        {"session", restored_session}, {"line", combined_line.at("line")},
        {"checkpoint_identity", combined_line.at("checkpoint_identity")},
        {"output_digits", 40}});

    const auto local_value = [&](const json::object& evaluation) {
      const auto& raw = evaluation.at("value").as_object()
                            .at("coefficients").as_array().front();
      return raw.is_string()
          ? std::stod(std::string(raw.as_string()))
          : std::stod(std::string(raw.as_array().front().as_string()));
    };
    const bool ok =
        created.at("parallel_arm_march_capability") ==
            "retained-native-concurrent-two-arm-march-v1" &&
        marched.at("capability") ==
            "retained-native-concurrent-two-arm-march-v1" &&
        marched.at("atomic_publication") == true &&
        marched.at("workers") == 2 &&
        marched.at("max_parallel_arms") == 2 &&
        marched.at("worker_overlap") == true &&
        marched.at("json_coefficients") == 0 &&
        lower.at("matches") == 2 && lower.at("tiles") == 3 &&
        upper.at("matches") == 2 && upper.at("tiles") == 3 &&
        lower_line.at("capability") == "retained-native-line-aggregate-v1" &&
        upper_line.at("capability") == "retained-native-line-aggregate-v1" &&
        combined_line.at("capability") ==
            "retained-native-line-aggregate-v1" &&
        lower_evaluation.at("status") == "ok" &&
        upper_evaluation.at("status") == "ok" &&
        std::abs(local_value(lower_evaluation) - 2.0) < 1e-30 &&
        std::abs(local_value(upper_evaluation) - 2.0) < 1e-30 &&
        std::abs(exported_coefficient(lower_export) + 8.0 / 3.0) < 1e-30 &&
        std::abs(exported_coefficient(upper_export) - 8.0 / 3.0) < 1e-30 &&
        std::abs(exported_coefficient(combined_export) - 16.0 / 3.0) <
            1e-30 &&
        restored_combined.at("status") == "ok" &&
        std::abs(exported_coefficient(restored_combined) - 16.0 / 3.0) <
            1e-30;
    if (!ok) {
      std::cerr << "failed response: " << json::serialize(failed) << '\n'
                << "marched: " << json::serialize(marched) << '\n'
                << "lower export: " << json::serialize(lower_export) << '\n'
                << "upper export: " << json::serialize(upper_export) << '\n'
                << "combined export: " << json::serialize(combined_export)
                << '\n' << "restored: " << json::serialize(restored)
                << '\n';
    }
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", restored_session}});
    restored_session.clear();
    std::remove(checkpoint_path.c_str());
    std::cout << (ok ? "PASS" : "FAIL")
              << ": persistent concurrent two-arm native march\n";
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    if (!restored_session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored_session}});
    std::remove(checkpoint_path.c_str());
    std::cerr << "FAIL: persistent parallel arms: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
