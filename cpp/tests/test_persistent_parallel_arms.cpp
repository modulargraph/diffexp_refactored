#include "diffexp2/json_codec.hpp"
#include "diffexp2/matching.hpp"

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
                          const std::string& domain,
                          const std::string& name,
                          const std::string& center,
                          bool nontrivial = false,
                          std::int32_t frame_base = 0) {
  json::array principal_row;
  principal_row.push_back(
      json::object{{"exact", nontrivial ? "t" : "0"},
                   {"proven_zero", !nontrivial}});
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
  json::array nhat_lags;
  nhat_lags.push_back(json::object{
      {"poly", json::array{}}, {"rat", json::array{}},
      {"val", json::array{nullptr}}});
  if (nontrivial) {
    json::array exponent;
    exponent.push_back(json::array{0, 0, "1"});
    nhat_lags.push_back(json::object{
        {"poly", json::array{json::object{
             {"s", 0}, {"e", std::move(exponent)}}}},
        {"rat", json::array{}}, {"val", json::array{0}}});
  }
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
           {"domain", domain}, {"precision_bits", 256},
           {"d", 1}, {"fb", frame_base}, {"w", 3},
           {"d_lags", std::move(d_lags)},
           {"denominators", json::array{}},
           {"nhat_lags", std::move(nhat_lags)},
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
                        const std::string& value,
                        std::uint32_t nmax = 0,
                        const std::string& next_epsilon_value = "0") {
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    shifts.emplace_back(std::to_string(n));
    schedule.push_back(json::array{json::object{
        {"case", n == 0 ? "R" : "T"},
        {"da", std::to_string(n)}, {"db", "0"}}});
  }
  const auto response = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart},
      {"run", json::object{
           {"nmax", nmax}, {"p", 0}, {"has_initial", true},
           {"adaptive_probe", false}, {"a_target", "0"},
           {"b_target", "0"}, {"a_shift_min", 0},
           {"a_shifts", std::move(shifts)},
           {"schedule", std::move(schedule)},
           {"initial", json::array{value, next_epsilon_value, "0"}},
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

json::object integrand_row(std::int32_t shift,
                           const std::string& identity) {
  const json::array one_kernel{"1", "0", "0", "0", "0"};
  const json::array zero_kernel{"0", "0", "0", "0", "0"};
  return json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", 1}, {"exact_identity", identity},
      {"entries", json::array{json::object{
           {"column", 0},
           {"multiplier", json::object{
                {"epsilon_shift", shift}, {"center_pole_order", 0},
                {"kernels", json::array{one_kernel, zero_kernel,
                                        zero_kernel}},
                {"exact_identity", identity},
                {"proven_zero", false}}}}}}};
}

json::object arm_execution(const std::string& first,
                           const std::string& second,
                           const std::vector<std::int32_t>& shifts) {
  if (shifts.size() != 3)
    throw std::invalid_argument("test arm requires three row shifts");
  return json::object{
      {"receiving_basis",
       json::array{json::array{first}, json::array{second}}},
      {"integrand_rows", json::array{
           integrand_row(shifts[0], "test-row-1:" +
                                      std::to_string(shifts[0])),
           integrand_row(shifts[1], "test-row-2:" +
                                      std::to_string(shifts[1])),
           integrand_row(shifts[2], "test-row-3:" +
                                      std::to_string(shifts[2]))}}};
}

json::object single_hop_arm_execution(const std::string& basis) {
  json::array basis_set;
  basis_set.emplace_back(basis);
  json::array receiving_basis;
  receiving_basis.push_back(std::move(basis_set));
  return json::object{
      {"receiving_basis", std::move(receiving_basis)},
      {"integrand_rows", json::array{
           integrand_row(0, "positive-leading-row-1"),
           integrand_row(0, "positive-leading-row-2")}}};
}

json::object run_arms_request(const std::string& session,
                              const std::string& plan,
                              const std::string& anchor,
                              json::object lower, json::object upper,
                              const std::string& root,
                              std::int32_t match_required_complete_max = 2,
                              std::int32_t work_min = -1) {
  return json::object{
      {"schema", 2}, {"op", "integration.run_arms"},
      {"session", session}, {"tile_plan", plan}, {"anchor", anchor},
      {"tile_plan_checkpoint_identity", "parallel-plan"},
      {"anchor_checkpoint_identity", "parallel-anchor"},
      {"epsilon", json::object{{"min", work_min}, {"max", 2},
                                {"required_complete_max", 1},
                                {"match_required_complete_max",
                                 match_required_complete_max}}},
      {"refinement", json::object{
           {"relative_tolerance", "1e-30"}, {"max_steps", 2}}},
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

bool exact_zero_acb_coefficient(const json::value& raw) {
  if (!raw.is_array()) return false;
  const auto& ball = raw.as_array();
  return ball.size() == 4 && ball[0].is_string() && ball[1].is_string() &&
      ball[2].is_string() && ball[3].is_string() &&
      std::stod(std::string(ball[0].as_string())) == 0.0 &&
      std::stod(std::string(ball[1].as_string())) == 0.0 &&
      ball[2].as_string() == "zero" && ball[3].as_string() == "zero";
}

double acb_real_midpoint(const json::value& raw) {
  return std::stod(std::string(raw.as_array().front().as_string()));
}

bool acb_full_rank_pivot_smoke() {
  diffexp2::ComplexBall::set_precision(256);
  using Matrix =
      diffexp2::matching_detail::DenseScalarMatrix<diffexp2::ComplexBall>;
  const auto uncertain = [] {
    return diffexp2::ComplexBall::from_strings("[0 +/- 1]");
  };
  const Matrix triangular{{diffexp2::ComplexBall(1), uncertain()},
                          {diffexp2::ComplexBall(0),
                           diffexp2::ComplexBall(1)}};
  if (diffexp2::matching_detail::certify_full_rank_by_nonzero_pivots(
          triangular, "parallel-arm triangular Acb rank proof") != 2)
    return false;

  const Matrix ambiguous{{uncertain(), diffexp2::ComplexBall(0)},
                         {diffexp2::ComplexBall(0),
                          diffexp2::ComplexBall(1)}};
  try {
    (void)diffexp2::matching_detail::certify_full_rank_by_nonzero_pivots(
        ambiguous, "parallel-arm ambiguous Acb rank proof");
  } catch (const diffexp2::MatchingArithmeticError& error) {
    return error.code == diffexp2::MatchingArithmeticErrorCode::AmbiguousZero;
  }
  return false;
}

bool positive_leading_match_smoke(const std::string& domain) {
  std::string session;
  try {
    json::object create{
        {"schema", 2}, {"op", "session.create"}, {"domain", domain},
        {"output_digits", 30}, {"chart_capacity", 4},
        {"local_capacity", 12}, {"match_capacity", 4},
        {"tile_plan_capacity", 2}};
    if (domain == "acb") create["precision_bits"] = 256;
    const auto created = request(std::move(create));
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_chart(
        session, domain, domain + "-positive-anchor", "0", false, 1);
    const auto lower_chart = prepare_chart(
        session, domain, domain + "-positive-lower", "-2/3");
    const auto upper_chart = prepare_chart(
        session, domain, domain + "-positive-upper", "2/3");
    const auto anchor = solve_local(
        session, anchor_chart, "0", "parallel-anchor", "1");
    const auto lower_basis = solve_local(
        session, lower_chart, "-2/3", domain + "-positive-lower-basis",
        "1");
    const auto upper_basis = solve_local(
        session, upper_chart, "2/3", domain + "-positive-upper-basis",
        "1");
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "parallel-plan"}, {"division_order", 3},
        {"lower", arm("-2/3", {anchor_chart, lower_chart})},
        {"upper", arm("2/3", {anchor_chart, upper_chart})}});
    if (planned.at("status") != "ok")
      throw std::runtime_error(
          "positive-leading plan: " + json::serialize(planned));
    const auto marched = request(run_arms_request(
        session, std::string(planned.at("tile_plan").as_string()), anchor,
        single_hop_arm_execution(lower_basis),
        single_hop_arm_execution(upper_basis),
        domain + "-positive-leading", 2, 0));
    if (marched.at("status") != "ok")
      throw std::runtime_error(
          "positive-leading march: " + json::serialize(marched));
    const auto& arms = marched.at("arms").as_object();
    const auto& lower_final =
        arms.at("lower").as_object().at("final_local").as_object();
    const auto& upper_final =
        arms.at("upper").as_object().at("final_local").as_object();
    const auto& combined = marched.at("combined_line_result").as_object();
    const auto evaluate_final = [&](const json::object& local) {
      return request(json::object{
          {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
          {"local", local.at("local")},
          {"point", json::object{{"exact", "0"}}},
          {"options", json::object{{"tail_estimate", false}}},
          {"output_digits", 40}});
    };
    const auto combined_export = request(json::object{
        {"schema", 2}, {"op", "integration.export"},
        {"session", session}, {"line", combined.at("line")},
        {"checkpoint_identity", combined.at("checkpoint_identity")},
        {"output_digits", 40}});
    bool coefficient_ok = false;
    if (domain == "rational") {
      coefficient_ok = lower_final.at("epsilon_min") == 1 &&
          upper_final.at("epsilon_min") == 1 &&
          combined.at("epsilon_min") == 1 &&
          combined_export.at("status") == "ok" &&
          combined_export.at("value").as_object().at("min") == 1 &&
          std::abs(exported_coefficient(combined_export) - 4.0 / 3.0) <
              1e-30;
    } else {
      const auto lower_value = evaluate_final(lower_final);
      const auto upper_value = evaluate_final(upper_final);
      const auto& lower_coefficients = lower_value.at("value").as_object()
                                           .at("coefficients").as_array();
      const auto& upper_coefficients = upper_value.at("value").as_object()
                                           .at("coefficients").as_array();
      const auto& combined_coefficients =
          combined_export.at("value").as_object()
              .at("coefficients").as_array();
      coefficient_ok = lower_final.at("epsilon_min") == 0 &&
          upper_final.at("epsilon_min") == 0 &&
          combined.at("epsilon_min") == 0 &&
          lower_value.at("status") == "ok" &&
          upper_value.at("status") == "ok" &&
          combined_export.at("status") == "ok" &&
          lower_value.at("value").as_object().at("min") == 1 &&
          upper_value.at("value").as_object().at("min") == 1 &&
          lower_coefficients.size() >= 1 &&
          upper_coefficients.size() >= 1 &&
          combined_coefficients.size() >= 2 &&
          exact_zero_acb_coefficient(combined_coefficients[0]) &&
          std::abs(acb_real_midpoint(lower_coefficients[0]) - 1.0) <
              1e-30 &&
          std::abs(acb_real_midpoint(upper_coefficients[0]) - 1.0) <
              1e-30 &&
          std::abs(acb_real_midpoint(combined_coefficients[1]) -
                   4.0 / 3.0) < 1e-30;
      if (!coefficient_ok)
        std::cerr << "Acb positive-leading coefficients: lower="
                  << json::serialize(lower_value.at("value"))
                  << " upper=" << json::serialize(upper_value.at("value"))
                  << " combined="
                  << json::serialize(combined_export.at("value")) << '\n';
    }
    const bool ok = coefficient_ok &&
        arms.at("lower").as_object().at("matches") == 1 &&
        arms.at("upper").as_object().at("matches") == 1 &&
        combined_export.at("status") == "ok";
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    session.clear();
    if (!ok)
      std::cerr << domain << " positive-leading match failed: mins="
                << lower_final.at("epsilon_min") << ','
                << upper_final.at("epsilon_min") << ','
                << combined.at("epsilon_min") << '\n';
    return ok;
  } catch (const std::exception& error) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    std::cerr << domain << " positive-leading match failed: "
              << error.what() << '\n';
    return false;
  }
}

bool acb_unit_leading_proof_smoke() {
  const std::string checkpoint_path =
      "/tmp/diffexp2-persistent-parallel-arms-acb.de2cp";
  std::remove(checkpoint_path.c_str());
  std::string session;
  std::string restored_session;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"}, {"domain", "acb"},
        {"precision_bits", 256}, {"output_digits", 30},
        {"chart_capacity", 8}, {"local_capacity", 16},
        {"match_capacity", 8}, {"tile_plan_capacity", 2}});
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_chart(
        session, "acb", "acb-pair-anchor", "0", true);
    const auto lower_1_chart = prepare_chart(
        session, "acb", "acb-pair-lower-1", "-2/3", true);
    const auto lower_2_chart = prepare_chart(
        session, "acb", "acb-pair-lower-2", "-4/3", true);
    const auto upper_1_chart = prepare_chart(
        session, "acb", "acb-pair-upper-1", "2/3", true);
    const auto upper_2_chart = prepare_chart(
        session, "acb", "acb-pair-upper-2", "4/3", true);
    const auto anchor = solve_local(session, anchor_chart, "0",
                                    "parallel-anchor", "2", 4);
    const auto lower_1 = solve_local(session, lower_1_chart, "-2/3",
                                     "acb-lower-basis-1", "1", 4);
    const auto lower_2 = solve_local(session, lower_2_chart, "-4/3",
                                     "acb-lower-basis-2", "1", 4);
    const auto defective_lower_1 = solve_local(
        session, lower_1_chart, "-2/3", "acb-lower-defective-1",
        "[0 +/- 1]", 4);
    const auto upper_1 = solve_local(session, upper_1_chart, "2/3",
                                     "acb-upper-basis-1", "1", 4);
    const auto upper_2 = solve_local(session, upper_2_chart, "4/3",
                                     "acb-upper-basis-2", "1", 4);
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "parallel-plan"}, {"division_order", 3},
        {"lower", arm("-4/3", {anchor_chart, lower_1_chart, lower_2_chart})},
        {"upper", arm("4/3", {anchor_chart, upper_1_chart, upper_2_chart})}});
    if (planned.at("status") != "ok")
      throw std::runtime_error("Acb plan: " + json::serialize(planned));
    const auto plan = std::string(planned.at("tile_plan").as_string());
    const auto lower_match = request(json::object{
        {"schema", 2}, {"op", "tile.match_interval"},
        {"session", session}, {"tile_plan", plan},
        {"arm", "lower"}, {"match", 0}});
    const auto evaluated_basis = request(json::object{
        {"schema", 2}, {"op", "local.evaluate"},
        {"session", session}, {"local", lower_1},
        {"point", json::object{{"exact",
             lower_match.at("receiving_local_exact")}}},
        {"options", json::object{{"tail_estimate", false}}}});
    const auto& basis_coefficient = evaluated_basis.at("value").as_object()
        .at("coefficients").as_array().front();
    const auto basis_midpoint = std::stod(std::string(
        basis_coefficient.as_array().front().as_string()));
    const auto before_defective = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto defective = request(run_arms_request(
        session, plan, anchor,
        arm_execution(defective_lower_1, lower_2, {0, 0, 0}),
        arm_execution(upper_1, upper_2, {0, 0, 0}),
        "acb-parallel-defective"));
    const auto after_defective = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto* defective_detail = defective.if_contains("detail");
    const bool defective_rejected = defective.at("status") == "error" &&
        defective_detail != nullptr && defective_detail->is_string() &&
        std::string(defective_detail->as_string())
            .find("overlap zero") != std::string::npos &&
        before_defective.at("locals") == after_defective.at("locals") &&
        before_defective.at("matches") == after_defective.at("matches") &&
        before_defective.at("line_results") ==
            after_defective.at("line_results") &&
        after_defective.at("pending_matches") == 0 &&
        after_defective.at("pending_local_solves") == 0 &&
        after_defective.at("pending_line_integrations") == 0;
    const auto marched = request(run_arms_request(
        session, plan, anchor,
        arm_execution(lower_1, lower_2, {0, -1, -1}),
        arm_execution(upper_1, upper_2, {0, -1, -1}),
        "acb-parallel-success"));
    if (marched.at("status") != "ok")
      throw std::runtime_error(
          "Acb success march: " + json::serialize(marched));
    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", session}, {"path", checkpoint_path},
        {"checkpoint_identity", "acb-parallel-roundtrip"}});
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    session.clear();
    const auto restored = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint_path},
        {"expected_identity", "acb-parallel-roundtrip"}});
    if (restored.at("status") == "ok")
      restored_session = std::string(restored.at("session").as_string());
    const auto restored_export = restored_session.empty()
        ? json::object{{"status", "error"}}
        : request(json::object{
              {"schema", 2}, {"op", "integration.export"},
              {"session", restored_session},
              {"line", marched.at("combined_line_result").as_object()
                           .at("line")},
              {"checkpoint_identity",
               marched.at("combined_line_result").as_object()
                   .at("checkpoint_identity")},
              {"output_digits", 30}});
    const bool ok = defective_rejected &&
        evaluated_basis.at("status") == "ok" &&
        std::abs(basis_midpoint - 1.0) > 1e-8 &&
        marched.at("status") == "ok" && saved.at("status") == "ok" &&
        restored.at("status") == "ok" &&
        restored_export.at("status") == "ok" &&
        marched.at("worker_overlap") == true &&
        marched.at("max_parallel_arms") == 2 &&
        marched.at("json_coefficients") == 0 &&
        marched.at("epsilon").as_object().at("required_complete_max") == 1 &&
        marched.at("epsilon").as_object().at(
            "match_required_complete_max") == 2 &&
        marched.at("combined_line_result").as_object().at("epsilon_min") == -1 &&
        marched.at("combined_line_result").as_object().at("epsilon_max") == 1;
    if (!ok)
      std::cerr << "Acb proof smoke: basis_midpoint=" << basis_midpoint
                << " defective=" << defective.at("status")
                << " detail="
                << (defective.if_contains("detail") != nullptr
                        ? std::string(defective.at("detail").as_string())
                        : std::string("<none>"))
                << " march=" << marched.at("status")
                << " save=" << saved.at("status")
                << " restore=" << restored.at("status")
                << " export=" << restored_export.at("status")
                << " combined_epsilon="
                << marched.at("combined_line_result").as_object()
                       .at("epsilon_min")
                << ":"
                << marched.at("combined_line_result").as_object()
                       .at("epsilon_max")
                << " defective_rejected=" << defective_rejected
                << " maps=" << before_defective.at("locals") << "/"
                << after_defective.at("locals") << ","
                << before_defective.at("matches") << "/"
                << after_defective.at("matches") << ","
                << before_defective.at("line_results") << "/"
                << after_defective.at("line_results") << " pending="
                << after_defective.at("pending_matches") << ","
                << after_defective.at("pending_local_solves") << ","
                << after_defective.at("pending_line_integrations") << '\n';
    if (!restored_session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored_session}});
    restored_session.clear();
    std::remove(checkpoint_path.c_str());
    return ok;
  } catch (const std::exception& error) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    if (!restored_session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored_session}});
    std::remove(checkpoint_path.c_str());
    std::cerr << "Acb unit-leading proof smoke failed: " << error.what()
              << '\n';
    return false;
  }
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
        {"chart_capacity", 8}, {"local_capacity", 32},
        {"match_capacity", 8}, {"tile_plan_capacity", 2}});
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_chart(
        session, "rational", "pair-anchor", "0");
    const auto lower_1_chart = prepare_chart(
        session, "rational", "pair-lower-1", "-2/3");
    const auto lower_guard_chart = prepare_chart(
        session, "rational", "pair-lower-guard", "-2/3", false, -1);
    const auto lower_2_chart = prepare_chart(
        session, "rational", "pair-lower-2", "-4/3");
    const auto upper_1_chart = prepare_chart(
        session, "rational", "pair-upper-1", "2/3");
    const auto upper_2_chart = prepare_chart(
        session, "rational", "pair-upper-2", "4/3");
    const auto anchor = solve_local(session, anchor_chart, "0",
                                    "parallel-anchor", "2");
    const auto lower_1 = solve_local(session, lower_1_chart, "-2/3",
                                     "parallel-lower-basis-1", "1");
    const auto lower_guard = solve_local(
        session, lower_guard_chart, "-2/3",
        "parallel-lower-nonzero-discarded-frame", "1");
    const auto lower_guard_zero = solve_local(
        session, lower_guard_chart, "-2/3",
        "parallel-lower-exact-zero-discarded-frame", "0", 0, "1");
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

    const auto guard_planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "parallel-plan"}, {"division_order", 3},
        {"lower", arm("-4/3", {anchor_chart, lower_guard_chart,
                                  lower_2_chart})},
        {"upper", arm("4/3", {anchor_chart, upper_1_chart,
                                 upper_2_chart})}});
    if (guard_planned.at("status") != "ok")
      throw std::runtime_error(
          "lower-frame guard plan: " + json::serialize(guard_planned));
    const auto guard_plan = std::string(
        guard_planned.at("tile_plan").as_string());
    const auto before_lower_guard = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto discarded_lower = request(run_arms_request(
        session, guard_plan, anchor,
        arm_execution(lower_guard, lower_2, {0, 0, 0}),
        arm_execution(upper_1, upper_2, {0, 0, 0}),
        "parallel-nonzero-discarded-lower", 1, 0));
    const auto after_lower_guard = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    if (discarded_lower.at("status") != "error" ||
        std::string(discarded_lower.at("detail").as_string())
                .find("discard a nonzero lower epsilon coefficient") ==
            std::string::npos ||
        before_lower_guard.at("locals") != after_lower_guard.at("locals") ||
        before_lower_guard.at("matches") !=
            after_lower_guard.at("matches") ||
        before_lower_guard.at("line_results") !=
            after_lower_guard.at("line_results") ||
        after_lower_guard.at("pending_matches") != 0 ||
        after_lower_guard.at("pending_local_solves") != 0 ||
        after_lower_guard.at("pending_line_integrations") != 0)
      throw std::runtime_error(
          "nonzero discarded Rational lower frame was not rejected atomically: " +
          json::serialize(discarded_lower) + " / " +
          json::serialize(after_lower_guard));

    const auto zero_discarded_lower = request(run_arms_request(
        session, guard_plan, anchor,
        arm_execution(lower_guard_zero, lower_2, {0, 0, 0}),
        arm_execution(upper_1, upper_2, {0, 0, 0}),
        "parallel-zero-discarded-lower", 1, 0));
    if (zero_discarded_lower.at("status") != "ok")
      throw std::runtime_error(
          "exact-zero discarded Rational lower frame was rejected: " +
          json::serialize(zero_discarded_lower));
    const auto& zero_arms = zero_discarded_lower.at("arms").as_object();
    for (const auto& arm_name : {"lower", "upper"}) {
      const auto& arm_result = zero_arms.at(arm_name).as_object();
      (void)request(json::object{
          {"schema", 2}, {"op", "local.release"}, {"session", session},
          {"local", arm_result.at("final_local").as_object().at("local")}});
      (void)request(json::object{
          {"schema", 2}, {"op", "integration.release"},
          {"session", session},
          {"line", arm_result.at("line_result").as_object().at("line")}});
    }
    (void)request(json::object{
        {"schema", 2}, {"op", "integration.release"},
        {"session", session},
        {"line", zero_discarded_lower.at("combined_line_result")
                     .as_object().at("line")}});

    const auto before_failure = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto failed = request(run_arms_request(
        session, plan, anchor,
        arm_execution(lower_1, lower_2, {0, 0, 0}),
        // The first upper basis is deliberately bound to the lower chart.
        // Admission succeeds, one worker may complete, but publication must
        // still be all-or-nothing after the upper plan-binding failure.
        arm_execution(lower_1, upper_2, {0, 0, 0}),
        "parallel-failure"));
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
        session, plan, anchor,
        arm_execution(lower_1, lower_2, {0, -1, -1}),
        arm_execution(upper_1, upper_2, {0, -1, -1}),
        "parallel-success"));
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
    for (const auto& local : {anchor, lower_1, lower_guard,
                              lower_guard_zero, lower_2, upper_1, upper_2})
      (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                                 {"session", session}, {"local", local}});
    (void)request(json::object{{"schema", 2}, {"op", "tile.release"},
                               {"session", session}, {"tile_plan", plan}});
    (void)request(json::object{{"schema", 2}, {"op", "tile.release"},
                               {"session", session},
                               {"tile_plan", guard_plan}});
    for (const auto& chart : {anchor_chart, lower_1_chart, lower_2_chart,
                              lower_guard_chart, upper_1_chart,
                              upper_2_chart})
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
    const bool ok = positive_leading_match_smoke("rational") &&
        positive_leading_match_smoke("acb") &&
        acb_full_rank_pivot_smoke() &&
        acb_unit_leading_proof_smoke() &&
        created.at("parallel_arm_march_capability") ==
            "retained-native-concurrent-two-arm-march-v1" &&
        marched.at("capability") ==
            "retained-native-concurrent-two-arm-march-v1" &&
        marched.at("atomic_publication") == true &&
        marched.at("workers") == 2 &&
        marched.at("max_parallel_arms") == 2 &&
        marched.at("worker_overlap") == true &&
        marched.at("json_coefficients") == 0 &&
        marched.at("epsilon").as_object().at("required_complete_max") == 1 &&
        marched.at("epsilon").as_object().at(
            "match_required_complete_max") == 2 &&
        lower.at("matches") == 2 && lower.at("tiles") == 3 &&
        upper.at("matches") == 2 && upper.at("tiles") == 3 &&
        lower_line.at("capability") == "retained-native-line-aggregate-v1" &&
        lower_line.at("epsilon_min") == -1 &&
        lower_line.at("epsilon_max") == 1 &&
        upper_line.at("capability") == "retained-native-line-aggregate-v1" &&
        upper_line.at("epsilon_min") == -1 &&
        upper_line.at("epsilon_max") == 1 &&
        combined_line.at("capability") ==
            "retained-native-line-aggregate-v1" &&
        combined_line.at("epsilon_min") == -1 &&
        combined_line.at("epsilon_max") == 1 &&
        lower_evaluation.at("status") == "ok" &&
        upper_evaluation.at("status") == "ok" &&
        std::abs(local_value(lower_evaluation) - 2.0) < 1e-30 &&
        std::abs(local_value(upper_evaluation) - 2.0) < 1e-30 &&
        std::abs(exported_coefficient(lower_export, 0) + 2.0) <
            1e-30 &&
        std::abs(exported_coefficient(lower_export, 1) + 2.0 / 3.0) <
            1e-30 &&
        std::abs(exported_coefficient(upper_export, 0) - 2.0) < 1e-30 &&
        std::abs(exported_coefficient(upper_export, 1) - 2.0 / 3.0) <
            1e-30 &&
        std::abs(exported_coefficient(combined_export, 0) - 4.0) <
            1e-30 &&
        std::abs(exported_coefficient(combined_export, 1) - 4.0 / 3.0) <
            1e-30 &&
        restored_combined.at("status") == "ok" &&
        std::abs(exported_coefficient(restored_combined, 0) - 4.0) <
            1e-30 &&
        std::abs(exported_coefficient(restored_combined, 1) - 4.0 / 3.0) <
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
