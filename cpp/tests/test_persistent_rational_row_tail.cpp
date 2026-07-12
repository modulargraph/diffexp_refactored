#include "diffexp2/checkpoint.hpp"
#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include <unistd.h>

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

void require_ok(const json::object& response, const char* label) {
  if (response.at("status") != "ok")
    throw std::runtime_error(std::string(label) + ": " +
                             json::serialize(response));
}

json::array prescriptions() {
  return json::array{json::object{
      {"factor_exact", "t"}, {"sign", -1}, {"multiplicity", 1},
      {"leading_coefficient_sign", 1}}};
}

std::string prepare_chart(const std::string& session,
                          const std::string& domain,
                          const std::string& suffix) {
  json::array principal_row;
  principal_row.push_back(json::object{
      {"exact", "0"}, {"proven_zero", true}});
  json::array principal_matrix;
  principal_matrix.push_back(std::move(principal_row));
  json::object problem{
      {"domain", domain}, {"d", 1}, {"fb", 0}, {"w", 1},
      {"d_lags", nested(json::array{json::object{
           {"s", 0}, {"v", "1"}}})},
      {"denominators", json::array{}},
      {"nhat_lags", json::array{json::object{
           {"poly", json::array{}}, {"rat", json::array{}},
           {"val", json::array{nullptr}}}}},
      {"d0_inverse", "1"},
      {"blocks", nested(json::array{0})},
      {"assembly", json::object{
           {"identity", true}, {"poly", json::array{}},
           {"rat", json::array{}}, {"val", json::array{0}}}},
      {"chop_digits", 50}};
  if (domain == "acb") problem["precision_bits"] = 256;
  auto prepared = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", "row-tail-chart-" + suffix},
      {"identity", "row-tail-operator-" + suffix},
      {"analytic", json::object{
           {"geometry", json::object{
                {"center_exact", "0"}, {"scale_exact", "1"},
                {"radius_exact", "2"}, {"infinite_radius", false},
                {"prescriptions", prescriptions()}}},
           {"principal_matrix", std::move(principal_matrix)},
           {"native_scc_capabilities", json::object{
                {"regular", true}, {"identity_gauge", true},
                {"identity_v", true}, {"no_pseudo", true}}}}},
      {"scc", json::object{
           {"components", nested(json::array{0})},
           {"structural_edges", json::array{}},
           {"condensation_edges", json::array{}},
           {"topological_order", json::array{0}},
           {"coupling_depth", 0}}},
      {"problem", std::move(problem)}});
  require_ok(prepared, "chart.prepare");
  return std::string(prepared.at("chart").as_string());
}

json::object solve_anchor(const std::string& session,
                          const std::string& chart,
                          const std::string& checkpoint) {
  auto solved = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart},
      {"run", json::object{
           {"nmax", 0}, {"p", 0}, {"has_initial", true},
           {"adaptive_probe", false}, {"a_target", "0"},
           {"b_target", "0"}, {"a_shift_min", 0},
           {"a_shifts", json::array{"0"}},
           {"schedule", nested(json::array{json::object{
                {"case", "R"}, {"da", "0"}, {"db", "0"}}})},
           {"initial", json::array{"1"}},
           {"initial_validity", json::array{0}},
           {"source", nullptr}, {"return_u", false}}},
      {"metadata", json::object{
           {"chart", json::object{
                {"center_exact", "0"}, {"scale_exact", "1"},
                {"radius", "2"}, {"infinite_radius", false}}},
           {"tag", json::object{
                {"a", json::object{{"domain", "rational"},
                                     {"canonical", "0"}}},
                {"b", json::object{{"domain", "rational"},
                                     {"canonical", "0"}}}}},
           {"prescriptions", prescriptions()},
           {"checkpoint_identity", checkpoint}}}});
  require_ok(solved, "local.solve");
  if (solved.at("tail_majorant").as_object().at("status") != "certified")
    throw std::runtime_error(
        "fixture source did not acquire a certified regular tail: " +
        json::serialize(solved));
  return solved;
}

json::object topology() {
  return json::object{
      {"singular_points", json::array{}},
      {"boundary_points", json::array{}},
      {"complex_projections", json::array{}},
      {"branch_sheets", json::array{json::object{
           {"factor_exact", "t"}, {"sign", -1}}}}};
}

json::object prepare_plan(const std::string& session,
                          const std::string& chart,
                          const std::string& checkpoint) {
  auto planned = request(json::object{
      {"schema", 2}, {"op", "tile.plan_arm"}, {"session", session},
      {"checkpoint_identity", checkpoint}, {"division_order", 3},
      {"arm", json::object{
           {"from_exact", "0"}, {"to_exact", "1/2"},
           {"charts", json::array{chart}}, {"topology", topology()}}}});
  require_ok(planned, "tile.plan_arm");
  return planned;
}

json::object run_state(const std::string& session,
                       const json::object& plan,
                       const std::string& plan_checkpoint,
                       const json::object& anchor,
                       const std::string& root) {
  auto state = request(json::object{
      {"schema", 2}, {"op", "transport.run_arm"},
      {"session", session}, {"tile_plan", plan.at("tile_plan")},
      {"anchor", anchor.at("local")},
      {"tile_plan_checkpoint_identity", plan_checkpoint},
      {"anchor_checkpoint_identity", anchor.at("checkpoint_identity")},
      {"arm", "upper"}, {"receiving_basis", json::array{}},
      {"epsilon", json::object{
           {"min", 0}, {"max", 0}, {"required_complete_max", 0},
           {"match_required_complete_max", 0}}},
      {"refinement", json::object{
           {"relative_tolerance", "1e-30"}, {"max_steps", 0}}},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
           {"root", root}}}});
  require_ok(state, "transport.run_arm");
  if (state.at("tiles") != 1 || state.at("matches") != 0)
    throw std::runtime_error(
        "rational-row tail fixture did not retain one match-free tile");
  return state;
}

json::object rational_row(const std::string& identity) {
  return json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", 1}, {"exact_identity", identity},
      {"entries", json::array{json::object{
           {"column", 0},
           {"multiplier", json::object{
                {"epsilon_shift", 0}, {"center_pole_order", 0},
                {"kernels", nested(json::array{"1"})},
                {"analytic_coefficients", json::array{json::object{
                     {"numerator", json::array{"1"}},
                     {"denominator", json::array{"1", "-1/2"}}}}},
                {"exact_identity", identity},
                {"proven_zero", false}}}}}}};
}

json::object contract_require(const std::string& session,
                              const json::object& state,
                              const std::string& root,
                              const std::string& checkpoint) {
  auto contracted = request(json::object{
      {"schema", 2}, {"op", "transport.contract"},
      {"session", session},
      {"transport_state", state.at("transport_state")},
      {"transport_state_checkpoint_identity", state.at("checkpoint_identity")},
      {"transport_state_provenance_identity", state.at("provenance_identity")},
      {"checkpoint_policy", json::object{
           {"schema",
            "diffexp2-deterministic-transport-contraction-checkpoints-v1"},
           {"root", root}}},
      {"observables", json::array{json::object{
           {"identity", root + ":observable"},
           {"checkpoint_identity", checkpoint},
           {"integrand_rows", json::array{
                rational_row(root + ":row")}},
           {"epsilon", json::object{
                {"min", 0}, {"max", 0},
                {"required_complete_max", 0}}},
           {"tail_policy", "require"}}}}});
  require_ok(contracted, "transport.contract require");
  if (contracted.at("lines").as_array().size() != 1)
    throw std::runtime_error(
        "required rational-row contraction returned the wrong line count");
  return contracted.at("lines").as_array().front().as_object();
}

json::object export_line(const std::string& session,
                         const json::object& line) {
  auto exported = request(json::object{
      {"schema", 2}, {"op", "integration.export"},
      {"session", session}, {"line", line.at("line")},
      {"checkpoint_identity", line.at("checkpoint_identity")},
      {"output_digits", 50}});
  require_ok(exported, "integration.export");
  return exported;
}

double coefficient_midpoint(const json::object& exported) {
  const auto& encoded = exported.at("value").as_object()
      .at("coefficients").as_array().front().as_array();
  return std::stod(std::string(encoded.front().as_string()));
}

const json::object& retained_line(const json::object& payload,
                                  const std::string& handle) {
  for (const auto& raw : payload.at("retained_line_results").as_array()) {
    const auto& line = raw.as_object();
    if (line.at("handle").as_string() == handle) return line;
  }
  throw std::runtime_error("checkpoint omitted certified rational-row line");
}

void run_domain(const std::string& domain) {
  const auto suffix = domain + "-" + std::to_string(::getpid());
  const auto checkpoint_path =
      (std::filesystem::temp_directory_path() /
       ("diffexp2-rational-row-tail-" + suffix + ".de2cp")).string();
  std::filesystem::remove(checkpoint_path);
  std::string session;
  std::string restored_session;
  try {
    json::object create{
        {"schema", 2}, {"op", "session.create"}, {"domain", domain},
        {"output_digits", 50}, {"chart_capacity", 1},
        {"local_capacity", 4}, {"match_capacity", 1},
        {"tile_plan_capacity", 1}, {"transport_state_capacity", 1},
        {"line_result_capacity", 2}};
    if (domain == "acb") create["precision_bits"] = 256;
    const auto created = request(std::move(create));
    require_ok(created, "session.create");
    session = std::string(created.at("session").as_string());
    const auto chart = prepare_chart(session, domain, suffix);
    const auto anchor_checkpoint = "row-tail-anchor-" + suffix;
    const auto anchor = solve_anchor(session, chart, anchor_checkpoint);
    const auto plan_checkpoint = "row-tail-plan-" + suffix;
    const auto plan = prepare_plan(session, chart, plan_checkpoint);
    const auto state = run_state(
        session, plan, plan_checkpoint, anchor, "row-tail-state-" + suffix);
    const auto line_checkpoint = "row-tail-line-" + suffix;
    const auto line = contract_require(
        session, state, "row-tail-contract-" + suffix, line_checkpoint);
    const auto exported = export_line(session, line);
    const auto& error = exported.at("value").as_object().at("error").as_object();
    const auto error_upper = error.at("absolute_upper_approx")
        .as_array().front().as_double();
    const auto actual_tail = 2.0 * std::log(4.0 / 3.0) - 0.5;
    if (exported.at("scope") != "full_local_with_certified_tail" ||
        exported.at("error_guarantee") != "certified" ||
        error.at("guarantee") != "certified" ||
        std::abs(coefficient_midpoint(exported) - 0.5) > 1e-30 ||
        !std::isfinite(error_upper) || error_upper < actual_tail)
      throw std::runtime_error(
          "required rational-row certificate is not a valid full-line bound: " +
          json::serialize(exported));

    const auto save_identity = "row-tail-save-" + suffix;
    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", session}, {"path", checkpoint_path},
        {"checkpoint_identity", save_identity}});
    require_ok(saved, "checkpoint.save");
    const auto payload = json::parse(
        diffexp2::checkpoint::read(checkpoint_path).payload_json).as_object();
    const auto& saved_line = retained_line(
        payload, std::string(line.at("line").as_string()));
    const auto& saved_result = saved_line.at("result").as_object();
    const auto& saved_error = saved_result.at("value").as_object()
        .at("error").as_object();
    if (saved_result.at("scope") != "full_local_with_certified_tail" ||
        saved_error.at("guarantee") != "certified" ||
        saved_error.at("absolute_exact").as_array().empty())
      throw std::runtime_error(
          "checkpoint did not retain the exact rational-row tail bound");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"}, {"session", session}}),
        "session.close");
    session.clear();

    const auto restored = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint_path}, {"expected_identity", save_identity}});
    require_ok(restored, "checkpoint.restore");
    restored_session = std::string(restored.at("session").as_string());
    const auto restored_export = export_line(restored_session, line);
    if (restored_export.at("scope") != "full_local_with_certified_tail" ||
        restored_export.at("error_guarantee") != "certified" ||
        restored_export.at("value") != exported.at("value"))
      throw std::runtime_error(
          "restored rational-row line changed its certified value/bound");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"},
        {"session", restored_session}}), "restored session.close");
    restored_session.clear();
    std::filesystem::remove(checkpoint_path);
  } catch (...) {
    if (!session.empty())
      (void)request(json::object{
          {"schema", 2}, {"op", "session.close"}, {"session", session}});
    if (!restored_session.empty())
      (void)request(json::object{
          {"schema", 2}, {"op", "session.close"},
          {"session", restored_session}});
    std::filesystem::remove(checkpoint_path);
    throw;
  }
}

}  // namespace

int main() {
  try {
    run_domain("rational");
    run_domain("acb");
    std::cout << "persistent Rational/Acb rational-row certified tails passed\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
