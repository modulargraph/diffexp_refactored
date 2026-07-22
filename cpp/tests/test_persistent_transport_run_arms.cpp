#include "diffexp2/checkpoint.hpp"
#include "diffexp2/json_codec.hpp"
#include "diffexp2/local_solution.hpp"
#include "diffexp2/scalar.hpp"

#include <boost/json.hpp>

#include <array>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
      .as_object();
}
void require_ok(const json::object& response, const char* label) {
  if (response.at("status") != "ok")
    throw std::runtime_error(std::string(label) + ": " +
                             json::serialize(response));
}

std::uint64_t counter(const json::object& object, const char* key) {
  const auto& value = object.at(key);
  if (value.is_uint64()) return value.as_uint64();
  if (value.is_int64() && value.as_int64() >= 0)
    return static_cast<std::uint64_t>(value.as_int64());
  throw std::runtime_error(std::string("expected nonnegative counter: ") + key);
}

diffexp2::ComplexBall decode_encoded_ball(const json::value& raw) {
  const auto& encoded = raw.as_array();
  if (encoded.size() != 4 || !encoded[0].is_string() ||
      !encoded[1].is_string() || !encoded[2].is_string() ||
      !encoded[3].is_string())
    throw std::runtime_error("unexpected encoded Acb coefficient");
  auto ball = diffexp2::ComplexBall::from_strings(
      std::string(encoded[0].as_string()),
      std::string(encoded[1].as_string()));
  const auto add_radius = [](arb_t value, const json::value& exponent) {
    const auto text = std::string(exponent.as_string());
    if (text != "zero")
      arb_add_error_2exp_si(
          value, static_cast<slong>(std::stol(text)));
  };
  add_radius(acb_realref(ball.raw()), encoded[2]);
  add_radius(acb_imagref(ball.raw()), encoded[3]);
  return ball;
}

bool epsilon_vectors_overlap(const json::value& left_raw,
                             const json::value& right_raw) {
  const auto& left = left_raw.as_object();
  const auto& right = right_raw.as_object();
  if (left.at("min") != right.at("min") ||
      left.at("max") != right.at("max") ||
      left.at("dimension") != right.at("dimension"))
    return false;
  const auto& left_coefficients = left.at("coefficients").as_array();
  const auto& right_coefficients = right.at("coefficients").as_array();
  if (left_coefficients.size() != right_coefficients.size()) return false;
  for (std::size_t index = 0; index < left_coefficients.size(); ++index) {
    const auto left_ball =
        decode_encoded_ball(left_coefficients[index]);
    const auto right_ball =
        decode_encoded_ball(right_coefficients[index]);
    if (!acb_overlaps(left_ball.raw(), right_ball.raw())) return false;
  }
  return true;
}

json::array prescriptions(bool branch_sensitive = true,
                          std::int32_t sign = -1) {
  if (!branch_sensitive) return json::array{};
  if (sign != -1 && sign != 1)
    throw std::invalid_argument("prescription sign must be +/-1");
  return json::array{json::object{
      {"factor_exact", "paired-state-f"}, {"sign", sign},
      {"multiplicity", 1}, {"leading_coefficient_sign", 1}}};
}

json::object epsilon_rational_one() {
  return json::object{{"zero", false}, {"valuation", 0},
                      {"numerator", json::array{"1"}},
                      {"denominator", json::array{"1"}}};
}

json::object epsilon_rational_minus_one() {
  return json::object{{"zero", false}, {"valuation", 0},
                      {"numerator", json::array{"-1"}},
                      {"denominator", json::array{"1"}}};
}

std::string prepare_chart(const std::string& session,
                          const std::string& name,
                          const std::string& center,
                          bool branch_sensitive = true,
                          const std::string& domain = "rational",
                          int chop_digits = 0,
                          std::int32_t prescription_sign = -1,
                          const std::string& relative_accuracy_max =
                              "1/100",
                          bool near_boundary_q_zero = false) {
  json::array principal_row{
      json::object{{"exact", "0"}, {"proven_zero", true}}};
  json::array component{0};
  json::array d_lag{json::object{{"s", 0}, {"v", "1"}}};
  json::array principal_matrix;
  principal_matrix.push_back(std::move(principal_row));
  json::array components;
  components.push_back(component);
  json::array d_lags;
  d_lags.push_back(std::move(d_lag));
  if (near_boundary_q_zero)
    d_lags.push_back(json::array{
        json::object{{"s", 0}, {"v", "-1"}}});
  json::array blocks;
  blocks.push_back(std::move(component));
  const auto owner = "de2-operator-" + name;
  const auto one = epsilon_rational_one();
  json::array physical_q{one};
  if (near_boundary_q_zero)
    physical_q.push_back(epsilon_rational_minus_one());
  json::array physical_c;
  physical_c.push_back(json::array{});
  json::object problem{
           {"domain", domain},
           {"d", 1}, {"fb", 0}, {"w", 3},
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
           {"physical_ode", json::object{
                {"schema", "diffexp2-physical-cleared-ode-v1"},
                {"basis", "physical-original-master"},
                {"theta_coordinate", "local-t"},
                {"owner_signature_identity", owner},
                {"payload_identity", "de2-physical-ode-" + name},
                {"q", std::move(physical_q)},
                {"c", std::move(physical_c)}}},
           {"chop_digits", chop_digits}};
  if (domain == "acb") problem["precision_bits"] = 256;
  json::object analytic{
      {"geometry", json::object{
           {"center_exact", center}, {"scale_exact", "1"},
           {"radius_exact", "2"}, {"infinite_radius", false},
           {"prescriptions", prescriptions(
                branch_sensitive, prescription_sign)}}},
      {"principal_matrix", std::move(principal_matrix)},
      {"native_scc_capabilities", json::object{
           {"regular", true}, {"identity_gauge", true},
           {"identity_v", true}, {"no_pseudo", true}}}};
  if (!relative_accuracy_max.empty())
    analytic["regular_value_relative_accuracy_max_exact"] =
        relative_accuracy_max;
  const auto prepared = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", name}, {"identity", owner},
      {"analytic", std::move(analytic)},
      {"scc", json::object{
           {"components", std::move(components)},
           {"structural_edges", json::array{}},
           {"condensation_edges", json::array{}},
           {"topological_order", json::array{0}},
           {"coupling_depth", 0}}},
      {"problem", std::move(problem)}});
  require_ok(prepared, "chart.prepare");
  return std::string(prepared.at("chart").as_string());
}

std::string solve_local(const std::string& session,
                        const std::string& chart,
                        const std::string& center,
                        const std::string& checkpoint,
                        const std::string& value,
                        std::uint32_t nmax = 0,
                        bool branch_sensitive = true,
                        std::int32_t prescription_sign = -1,
                        std::int32_t initial_validity = 2,
                        const std::string& equation_owner = {}) {
  json::array schedule;
  json::array shifts;
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    shifts.emplace_back(std::to_string(n));
    schedule.push_back(json::array{json::object{
        {"case", n == 0 ? "R" : "T"},
        {"da", std::to_string(n)}, {"db", "0"}}});
  }
  json::object solve_request{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart},
      {"run", json::object{
           {"nmax", nmax}, {"p", 0}, {"has_initial", true},
           {"adaptive_probe", false}, {"a_target", "0"},
           {"b_target", "0"}, {"a_shift_min", 0},
           {"a_shifts", std::move(shifts)},
           {"schedule", std::move(schedule)},
           {"initial", json::array{value, "0", "0"}},
           {"initial_validity", json::array{initial_validity}},
           {"source", nullptr},
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
           {"prescriptions", prescriptions(
                branch_sensitive, prescription_sign)},
           {"checkpoint_identity", checkpoint}}}};
  if (!equation_owner.empty())
    solve_request["equation_owner"] = equation_owner;
  const auto solved = request(std::move(solve_request));
  require_ok(solved, "local.solve");
  return std::string(solved.at("local").as_string());
}

std::string prepare_regular_equation_owner(
    const std::string& session, const std::string& name,
    const std::string& center, bool branch_sensitive = true,
    std::int32_t prescription_sign = -1,
    const std::string& relative_accuracy_max = "1/100") {
  const auto identity = "de2-equation-" + name;
  json::array physical_c;
  physical_c.emplace_back(json::array{});
  json::array physical_q;
  physical_q.emplace_back(epsilon_rational_one());
  const auto prepared = request(json::object{
      {"schema", 2}, {"op", "regular_equation.prepare"},
      {"session", session},
      {"capability",
       "frame-independent-regular-physical-equation-owner-v1"},
      {"key", "regular-equation:" + name}, {"identity", identity},
      {"dimension", 1},
      {"relative_accuracy_max_exact", relative_accuracy_max},
      {"geometry", json::object{
           {"center_exact", center}, {"scale_exact", "1"},
           {"radius_exact", "2"}, {"infinite_radius", false},
           {"prescriptions", prescriptions(
                branch_sensitive, prescription_sign)}}},
      {"physical_ode", json::object{
           {"schema", "diffexp2-physical-cleared-ode-v1"},
           {"basis", "physical-original-master"},
           {"theta_coordinate", "local-t"},
           {"owner_signature_identity", identity},
           {"payload_identity", "de2-physical-ode-" + name},
           {"q", std::move(physical_q)},
           {"c", std::move(physical_c)}}}});
  require_ok(prepared, "regular_equation.prepare");
  return std::string(prepared.at("equation_owner").as_string());
}

std::string solve_log_local(const std::string& session,
                            const std::string& chart,
                            const std::string& center,
                            const std::string& checkpoint,
                            std::int32_t prescription_sign = -1) {
  constexpr std::uint32_t nmax = 4;
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    shifts.emplace_back(std::to_string(n));
    schedule.push_back(json::array{json::object{
        {"case", n == 0 ? "R" : "T"},
        {"da", std::to_string(n)}, {"db", "0"}}});
  }
  const auto solved = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart},
      {"run", json::object{
           {"nmax", nmax}, {"p", 1}, {"has_initial", true},
           {"adaptive_probe", false}, {"a_target", "0"},
           {"b_target", "0"}, {"a_shift_min", 0},
           {"a_shifts", std::move(shifts)},
           {"schedule", std::move(schedule)},
           {"initial", json::array{"2", "0", "0", "1", "0", "0"}},
           {"initial_validity", json::array{2, 2}}, {"source", nullptr},
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
           {"prescriptions", prescriptions(true, prescription_sign)},
           {"checkpoint_identity", checkpoint}}}});
  require_ok(solved, "log local.solve");
  return std::string(solved.at("local").as_string());
}

std::string prepare_multiblock_chart(const std::string& session,
                                     const std::string& name,
                                     const std::string& center,
                                     const std::string& domain =
                                         "rational") {
  constexpr std::uint32_t dimension = 2;
  json::array principal_matrix;
  json::array components;
  json::array blocks;
  json::array order;
  json::array null_matrix;
  json::array assembly_matrix;
  json::array d_lag{json::object{{"s", 0}, {"v", "1"}}};
  json::array d_lags;
  d_lags.push_back(std::move(d_lag));
  json::array physical_c;
  physical_c.push_back(json::array{});
  for (std::uint32_t row = 0; row < dimension; ++row) {
    json::array principal_row;
    for (std::uint32_t column = 0; column < dimension; ++column) {
      principal_row.push_back(
          json::object{{"exact", "0"}, {"proven_zero", true}});
      null_matrix.push_back(nullptr);
      assembly_matrix.push_back(
          row == column ? json::value(0) : json::value(nullptr));
    }
    principal_matrix.push_back(std::move(principal_row));
    components.push_back(json::array{row});
    blocks.push_back(json::array{row});
    order.push_back(row);
  }
  const auto owner = "de2-operator-" + name;
  json::object problem{
      {"domain", domain}, {"d", dimension}, {"fb", 0},
      {"w", 3},
      {"d_lags", std::move(d_lags)},
      {"denominators", json::array{}},
      {"nhat_lags", json::array{json::object{
           {"poly", json::array{}}, {"rat", json::array{}},
           {"val", std::move(null_matrix)}}}},
      {"d0_inverse", "1"}, {"blocks", std::move(blocks)},
      {"assembly", json::object{
           {"identity", true}, {"poly", json::array{}},
           {"rat", json::array{}},
           {"val", std::move(assembly_matrix)}}},
      {"physical_ode", json::object{
           {"schema", "diffexp2-physical-cleared-ode-v1"},
           {"basis", "physical-original-master"},
           {"theta_coordinate", "local-t"},
           {"owner_signature_identity", owner},
           {"payload_identity", "de2-physical-ode-" + name},
           {"q", json::array{epsilon_rational_one()}},
           {"c", std::move(physical_c)}}},
      {"chop_digits", 0}};
  if (domain == "acb") problem["precision_bits"] = 256;
  const auto prepared = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", name}, {"identity", owner},
      {"analytic", json::object{
           {"regular_value_relative_accuracy_max_exact", "1/100"},
           {"geometry", json::object{
                {"center_exact", center}, {"scale_exact", "1"},
                {"radius_exact", "2"}, {"infinite_radius", false},
                {"prescriptions", prescriptions(true)}}},
           {"principal_matrix", std::move(principal_matrix)},
           {"native_scc_capabilities", json::object{
                {"regular", true}, {"identity_gauge", true},
                {"identity_v", true}, {"no_pseudo", true}}}}},
      {"scc", json::object{
           {"components", std::move(components)},
           {"structural_edges", json::array{}},
           {"condensation_edges", json::array{}},
           {"topological_order", std::move(order)},
           {"coupling_depth", 0}}},
      {"problem", std::move(problem)}});
  require_ok(prepared, "multiblock chart.prepare");
  return std::string(prepared.at("chart").as_string());
}

std::string solve_multiblock_local(
    const std::string& session, const std::string& chart,
    const std::string& center, const std::string& checkpoint,
    const std::array<std::string, 2>& values,
    std::uint32_t nmax = 4) {
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    shifts.emplace_back(std::to_string(n));
    json::array row;
    for (std::uint32_t block = 0; block < 2; ++block)
      row.push_back(json::object{
          {"case", n == 0 ? "R" : "T"},
          {"da", std::to_string(n)}, {"db", "0"}});
    schedule.push_back(std::move(row));
  }
  json::array initial;
  for (const auto& value : values) {
    initial.emplace_back(value);
    initial.push_back("0");
    initial.push_back("0");
  }
  const auto solved = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart},
      {"run", json::object{
           {"nmax", nmax}, {"p", 0}, {"has_initial", true},
           {"adaptive_probe", false}, {"a_target", "0"},
           {"b_target", "0"}, {"a_shift_min", 0},
           {"a_shifts", std::move(shifts)},
           {"schedule", std::move(schedule)},
           {"initial", std::move(initial)},
           {"initial_validity", json::array{2, 2}}, {"source", nullptr},
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
           {"prescriptions", prescriptions(true)},
           {"checkpoint_identity", checkpoint}}}});
  require_ok(solved, "multiblock local.solve");
  return std::string(solved.at("local").as_string());
}

json::object topology(bool branch_sensitive = true,
                      std::int32_t sign = -1) {
  return json::object{
      {"singular_points", json::array{}},
      {"boundary_points", json::array{}},
      {"complex_projections", json::array{}},
      {"branch_sheets", branch_sensitive
           ? json::array{json::object{
                 {"factor_exact", "paired-state-f"}, {"sign", sign}}}
           : json::array{}}};
}

json::object arm(const std::string& endpoint,
                 const std::string& anchor_chart,
                 const std::string& receiving_chart,
                 bool branch_sensitive = true,
                 std::int32_t prescription_sign = -1) {
  return json::object{
      {"from_exact", "0"}, {"to_exact", endpoint},
      {"charts", json::array{anchor_chart, receiving_chart}},
      {"topology", topology(branch_sensitive, prescription_sign)}};
}

json::object receiving_basis(const std::string& local) {
  json::array basis_set;
  basis_set.emplace_back(local);
  json::array basis;
  basis.push_back(std::move(basis_set));
  return json::object{
      {"receiving_basis", std::move(basis)}};
}

json::object epsilon_contract() {
  return json::object{{"min", 0}, {"max", 2},
                      {"required_complete_max", 1},
                      {"match_required_complete_max", 2}};
}

json::object refinement() {
  return json::object{{"relative_tolerance", "1e-30"},
                      {"max_steps", 2}};
}

json::object run_pair(const std::string& session,
                      const std::string& plan,
                      const std::string& anchor,
                      const std::string& lower_basis,
                      const std::string& upper_basis,
                      const std::string& root) {
  return request(json::object{
      {"schema", 2}, {"op", "transport.run_arms"},
      {"session", session}, {"tile_plan", plan}, {"anchor", anchor},
      {"tile_plan_checkpoint_identity", "paired-state-plan"},
      {"anchor_checkpoint_identity", "paired-state-anchor"},
      {"epsilon", epsilon_contract()}, {"refinement", refinement()},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
           {"root", root}}},
      {"lower", receiving_basis(lower_basis)},
      {"upper", receiving_basis(upper_basis)}});
}

json::object run_consuming_pair(const std::string& session,
                                const std::string& plan,
                                const std::string& anchor,
                                const std::string& lower_basis,
                                const std::string& upper_basis,
                                const std::string& root) {
  return request(json::object{
      {"schema", 2}, {"op", "transport.run_arms_consuming"},
      {"session", session}, {"tile_plan", plan}, {"anchor", anchor},
      {"tile_plan_checkpoint_identity", "consuming-state-plan"},
      {"anchor_checkpoint_identity", "consuming-state-anchor"},
      {"epsilon", epsilon_contract()}, {"refinement", refinement()},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
           {"root", root}}},
      {"lower", receiving_basis(lower_basis)},
      {"upper", receiving_basis(upper_basis)}});
}

json::object consume_hop(const std::string& session,
                         const std::string& plan,
                         const std::string& arm_name,
                         const std::string& incoming,
                         const std::string& incoming_checkpoint,
                         const std::string& basis,
                         const std::string& root,
                         const std::string& plan_checkpoint =
                             "streaming-state-plan") {
  return request(json::object{
      {"schema", 2}, {"op", "transport.consume_hop"},
      {"session", session}, {"tile_plan", plan},
      {"tile_plan_checkpoint_identity", plan_checkpoint},
      {"arm", arm_name}, {"match", 0},
      {"receiving_basis", json::array{basis}},
      {"incoming", incoming},
      {"incoming_checkpoint_identity", incoming_checkpoint},
      {"epsilon", json::object{{"min", 0}, {"max", 2},
                                  {"required_complete_max", 2}}},
      {"refinement", refinement()},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
           {"root", root}}}});
}

json::object value_solver(const std::string& center,
                          bool branch_sensitive = false,
                          std::uint32_t nmax = 4,
                          const std::string& relative_error_max = "1/100",
                          const std::string& tail_proxy_max = "1/100",
                          std::int32_t prescription_sign = -1) {
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    shifts.emplace_back(std::to_string(n));
    schedule.push_back(json::array{json::object{
        {"case", n == 0 ? "R" : "T"},
        {"da", std::to_string(n)}, {"db", "0"}}});
  }
  return json::object{
      {"schema", "diffexp2-native-regular-value-solver-prototype-v1"},
      {"run", json::object{
           {"nmax", nmax}, {"p", 0}, {"has_initial", true},
           {"adaptive_probe", false}, {"a_target", "0"},
           {"b_target", "0"}, {"a_shift_min", 0},
           {"a_shifts", std::move(shifts)},
           {"schedule", std::move(schedule)},
           {"initial", json::array{"1", "0", "0"}},
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
           {"prescriptions", prescriptions(
                branch_sensitive, prescription_sign)},
           {"checkpoint_identity", "value-solver-prototype"}}},
      {"tail_proxy_max_exact", tail_proxy_max},
      {"relative_accuracy_max_exact", relative_error_max}};
}

json::object consume_value_hop(
    const std::string& session, const std::string& plan,
    const std::string& arm_name, const std::string& incoming,
    const std::string& incoming_checkpoint, json::object solver,
    const std::string& root,
    const std::string& plan_checkpoint = "streaming-state-plan") {
  return request(json::object{
      {"schema", 2}, {"op", "transport.consume_value_hop"},
      {"session", session}, {"tile_plan", plan},
      {"tile_plan_checkpoint_identity", plan_checkpoint},
      {"arm", arm_name}, {"match", 0},
      {"value_solver", std::move(solver)},
      {"incoming", incoming},
      {"incoming_checkpoint_identity", incoming_checkpoint},
      {"epsilon", json::object{{"min", 0}, {"max", 2},
                                  {"required_complete_max", 2}}},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
           {"root", root}}}});
}

json::object physical_value_solver(
    const std::string& center, std::uint32_t nmax = 8,
    const std::string& relative_error_max = "1/100") {
  return json::object{
      {"schema", "diffexp2-native-ordinary-physical-value-solver-v1"},
      {"taylor_complete_max", nmax},
      {"metadata", json::object{
           {"chart", json::object{
                {"center_exact", center}, {"scale_exact", "1"},
                {"radius", "2"}, {"infinite_radius", false}}},
           {"tag", json::object{
                {"a", json::object{{"domain", "rational"},
                                     {"canonical", "0"}}},
                {"b", json::object{{"domain", "rational"},
                                     {"canonical", "0"}}}}},
           {"prescriptions", prescriptions(true)},
           {"checkpoint_identity", "physical-value-solver-prototype"}}},
      {"relative_accuracy_max_exact", relative_error_max}};
}

json::object consume_physical_value_hop(
    const std::string& session, const std::string& plan,
    const std::string& arm_name, std::size_t match_index,
    const std::string& incoming,
    const std::string& incoming_checkpoint, json::object solver,
    const std::string& root,
    const std::string& plan_checkpoint,
    std::int32_t epsilon_min = 0) {
  return request(json::object{
      {"schema", 2}, {"op", "transport.consume_physical_value_hop"},
      {"session", session}, {"tile_plan", plan},
      {"tile_plan_checkpoint_identity", plan_checkpoint},
      {"arm", arm_name}, {"match", match_index},
      {"value_solver", std::move(solver)},
      {"incoming", incoming},
      {"incoming_checkpoint_identity", incoming_checkpoint},
      {"epsilon", json::object{{"min", epsilon_min}, {"max", 2},
                                  {"required_complete_max", 2}}},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
           {"root", root}}}});
}

json::object consume_basis_hop(
    const std::string& session, const std::string& plan,
    const std::string& plan_checkpoint, const std::string& arm_name,
    const std::string& incoming, const std::string& incoming_checkpoint,
    const std::vector<std::string>& basis, const std::string& root) {
  json::array encoded_basis;
  for (const auto& local : basis) encoded_basis.emplace_back(local);
  return request(json::object{
      {"schema", 2}, {"op", "transport.consume_hop"},
      {"session", session}, {"tile_plan", plan},
      {"tile_plan_checkpoint_identity", plan_checkpoint},
      {"arm", arm_name}, {"match", 0},
      {"receiving_basis", std::move(encoded_basis)},
      {"incoming", incoming},
      {"incoming_checkpoint_identity", incoming_checkpoint},
      {"epsilon", json::object{{"min", 0}, {"max", 2},
                                  {"required_complete_max", 2}}},
      {"refinement", refinement()},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
           {"root", root}}}});
}

json::object published_side(const std::string& anchor,
                            const json::object& hop) {
  const auto next = hop.at("next_local").as_object().at("local");
  return json::object{
      {"tile_sources", json::array{anchor, next}}};
}

json::object publish_consumed_states(
    const std::string& session, const std::string& plan,
    const std::string& anchor, const json::object& lower,
    const json::object& upper, const std::string& root,
    const std::string& plan_checkpoint = "streaming-state-plan",
    const std::string& anchor_checkpoint =
        "streaming-state-anchor") {
  return request(json::object{
      {"schema", 2}, {"op", "transport.publish_consumed_states"},
      {"session", session}, {"tile_plan", plan},
      {"tile_plan_checkpoint_identity", plan_checkpoint},
      {"anchor", anchor},
      {"anchor_checkpoint_identity", anchor_checkpoint},
      {"epsilon", epsilon_contract()}, {"refinement", refinement()},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
           {"root", root}}},
      {"lower", published_side(anchor, lower)},
      {"upper", published_side(anchor, upper)}});
}

json::object run_single(const std::string& session,
                        const std::string& plan,
                        const std::string& anchor,
                        const std::string& arm_name,
                        const std::string& basis,
                        const std::string& root) {
  auto receiving = receiving_basis(basis);
  return request(json::object{
      {"schema", 2}, {"op", "transport.run_arm"},
      {"session", session}, {"tile_plan", plan}, {"anchor", anchor},
      {"tile_plan_checkpoint_identity", "paired-state-plan"},
      {"anchor_checkpoint_identity", "paired-state-anchor"},
      {"arm", arm_name},
      {"receiving_basis", std::move(receiving.at("receiving_basis"))},
      {"epsilon", epsilon_contract()}, {"refinement", refinement()},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
           {"root", root}}}});
}

json::object integrand_row(const std::string& identity,
                           std::size_t taylor_width = 5) {
  if (taylor_width == 0)
    throw std::invalid_argument("integrand row Taylor width is zero");
  json::array one;
  json::array zero;
  one.reserve(taylor_width);
  zero.reserve(taylor_width);
  for (std::size_t n = 0; n < taylor_width; ++n) {
    one.emplace_back(n == 0 ? "1" : "0");
    zero.emplace_back("0");
  }
  return json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", 1}, {"exact_identity", identity},
      {"entries", json::array{json::object{
           {"column", 0},
           {"multiplier", json::object{
                {"epsilon_shift", 0}, {"center_pole_order", 0},
                {"kernels", json::array{one, zero, zero}},
                {"exact_identity", identity},
                {"proven_zero", false}}}}}}};
}

json::object cancellation_integrand_row(
    const std::string& identity,
    std::size_t taylor_width = 5) {
  if (taylor_width == 0)
    throw std::invalid_argument(
        "cancellation integrand row Taylor width is zero");
  const auto kernels = [&](const std::string& constant) {
    json::array leading;
    json::array zero;
    leading.reserve(taylor_width);
    zero.reserve(taylor_width);
    for (std::size_t n = 0; n < taylor_width; ++n) {
      leading.emplace_back(n == 0 ? constant : "0");
      zero.emplace_back("0");
    }
    return json::array{std::move(leading), zero, zero};
  };
  const auto analytic = [&](const std::string& constant) {
    return json::array{
        json::object{{"numerator", json::array{constant}},
                     {"denominator", json::array{"1"}}},
        json::object{{"numerator", json::array{"0"}},
                     {"denominator", json::array{"1"}}},
        json::object{{"numerator", json::array{"0"}},
                     {"denominator", json::array{"1"}}}};
  };
  return json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", 2}, {"exact_identity", identity},
      {"entries", json::array{
           json::object{
               {"column", 0},
               {"multiplier", json::object{
                    {"epsilon_shift", 0},
                    {"center_pole_order", 0},
                    {"kernels", kernels("1")},
                    {"analytic_coefficients", analytic("1")},
                    {"exact_identity", identity + ":plus"},
                    {"proven_zero", false}}}},
           json::object{
               {"column", 1},
               {"multiplier", json::object{
                    {"epsilon_shift", 0},
                    {"center_pole_order", 0},
                    {"kernels", kernels("-1")},
                    {"analytic_coefficients", analytic("-1")},
                    {"exact_identity", identity + ":minus"},
                    {"proven_zero", false}}}}}}};
}

json::value contract_and_export(const std::string& session,
                                const json::object& state,
                                const std::string& root,
                                std::size_t taylor_width = 5) {
  const auto contracted = request(json::object{
      {"schema", 2}, {"op", "transport.contract"},
      {"session", session},
      {"transport_state", state.at("transport_state")},
      {"transport_state_checkpoint_identity",
       state.at("checkpoint_identity")},
      {"transport_state_provenance_identity",
       state.at("provenance_identity")},
      {"checkpoint_policy", json::object{
           {"schema",
            "diffexp2-deterministic-transport-contraction-checkpoints-v1"},
           {"root", root}}},
      {"observables", json::array{json::object{
           {"identity", root + ":observable"},
           {"checkpoint_identity", root + ":line"},
           {"integrand_rows", json::array{
                integrand_row(root + ":row:1", taylor_width),
                integrand_row(root + ":row:2", taylor_width)}},
           {"epsilon", json::object{{"min", 0}, {"max", 2},
                                      {"required_complete_max", 1}}},
           {"tail_policy", "stored"}}}}});
  require_ok(contracted, "transport.contract");
  const auto& line = contracted.at("lines").as_array().front().as_object();
  const auto exported = request(json::object{
      {"schema", 2}, {"op", "integration.export"},
      {"session", session}, {"line", line.at("line")},
      {"checkpoint_identity", line.at("checkpoint_identity")},
      {"output_digits", 40}});
  require_ok(exported, "integration.export");
  const auto value = exported.at("value");
  require_ok(request(json::object{
      {"schema", 2}, {"op", "integration.release"},
      {"session", session}, {"line", line.at("line")}}),
      "integration.release");
  return value;
}

json::value contract_cancellation_and_export(
    const std::string& session, const json::object& state,
    const std::string& root, std::size_t taylor_width = 5,
    bool release_line = true) {
  const auto contracted = request(json::object{
      {"schema", 2}, {"op", "transport.contract"},
      {"session", session},
      {"transport_state", state.at("transport_state")},
      {"transport_state_checkpoint_identity",
       state.at("checkpoint_identity")},
      {"transport_state_provenance_identity",
       state.at("provenance_identity")},
      {"checkpoint_policy", json::object{
           {"schema",
            "diffexp2-deterministic-transport-contraction-checkpoints-v1"},
           {"root", root}}},
      {"observables", json::array{json::object{
           {"identity", root + ":observable"},
           {"checkpoint_identity", root + ":line"},
           {"integrand_rows", json::array{
                cancellation_integrand_row(
                    root + ":row:1", taylor_width),
                cancellation_integrand_row(
                    root + ":row:2", taylor_width)}},
           {"epsilon", json::object{
                {"min", 0}, {"max", 2},
                {"required_complete_max", 1}}},
           {"tail_policy", "stored"},
           {"divergent_cancellation", json::object{
                {"mode", "bounded-relative-acb"},
                {"relative_tolerance", "1e-40"},
                {"provenance",
                 "terminal-single-contraction-policy-regression"}}}}}}});
  require_ok(contracted, "terminal cancellation transport.contract");
  const auto& line =
      contracted.at("lines").as_array().front().as_object();
  const auto exported = request(json::object{
      {"schema", 2}, {"op", "integration.export"},
      {"session", session}, {"line", line.at("line")},
      {"checkpoint_identity", line.at("checkpoint_identity")},
      {"output_digits", 40}});
  require_ok(exported, "terminal cancellation integration.export");
  const auto value = exported.at("value");
  if (release_line)
    require_ok(request(json::object{
        {"schema", 2}, {"op", "integration.release"},
        {"session", session}, {"line", line.at("line")}}),
        "terminal cancellation integration.release");
  return value;
}

json::object session_stats(const std::string& session) {
  return request(json::object{{"schema", 2}, {"op", "session.stats"},
                              {"session", session}});
}

void release_local(const std::string& session, const std::string& local) {
  require_ok(request(json::object{{"schema", 2}, {"op", "local.release"},
                                  {"session", session}, {"local", local}}),
             "local.release");
}

void release_state(const std::string& session, const std::string& state) {
  require_ok(request(json::object{
      {"schema", 2}, {"op", "transport.release"},
      {"session", session}, {"transport_state", state}}),
      "transport.release");
}

void test_regular_value_hop_checkpoint() {
  const std::string checkpoint =
      "/tmp/diffexp2-regular-value-hop-roundtrip.de2cp";
  const std::string tampered_checkpoint =
      "/tmp/diffexp2-regular-value-hop-retained-value-tampered.de2cp";
  std::remove(checkpoint.c_str());
  std::remove(tampered_checkpoint.c_str());
  std::string session;
  std::string restored;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "acb"}, {"precision_bits", 256},
        {"output_digits", 40},
        {"chart_capacity", 6}, {"local_capacity", 12},
        {"match_capacity", 4}, {"tile_plan_capacity", 2},
        {"transport_state_capacity", 4}, {"line_result_capacity", 4}});
    require_ok(created, "value-hop session.create");
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_chart(
        session, "value-hop-anchor-chart", "0", true, "acb", 4, -1,
        "1/100", true);
    const auto lower_chart = prepare_chart(
        session, "value-hop-lower-chart", "-2/3", true, "acb", 4);
    const auto upper_chart = prepare_chart(
        session, "value-hop-upper-chart", "2/3", true, "acb", 4);
    const auto anchor = solve_local(
        session, anchor_chart, "0", "streaming-state-anchor", "2", 30,
        true);
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "streaming-state-plan"},
        {"division_order", 3},
        {"lower", arm("-2/3", anchor_chart, lower_chart)},
        {"upper", arm("2/3", anchor_chart, upper_chart)}});
    require_ok(planned, "value-hop tile.plan");
    const auto plan = std::string(planned.at("tile_plan").as_string());
    const auto legacy_chart = prepare_chart(
        session, "value-hop-legacy-chart", "-1/3", true, "acb", 4,
        -1, "");
    const auto legacy_local = solve_local(
        session, legacy_chart, "-1/3", "value-hop-legacy-local", "1", 0,
        true);
    release_local(session, legacy_local);
    const auto legacy_planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "value-hop-legacy-plan"},
        {"division_order", 3},
        {"lower", arm("-1/3", anchor_chart, legacy_chart)},
        {"upper", arm("2/3", anchor_chart, upper_chart)}});
    require_ok(legacy_planned, "legacy-contract value-hop tile.plan");
    const auto legacy_plan = std::string(
        legacy_planned.at("tile_plan").as_string());
    const auto before_legacy_hop = session_stats(session);
    const auto legacy_hop = consume_value_hop(
        session, legacy_plan, "lower", anchor, "streaming-state-anchor",
        value_solver("-1/3", true), "value-hop-legacy-contract",
        "value-hop-legacy-plan");
    require_ok(legacy_hop, "legacy-contract value hop");
    const auto after_legacy_hop = session_stats(session);
    if (legacy_hop.at("used") != false ||
        legacy_hop.at("reason") !=
            "receiving-owner-has-no-relative-accuracy-contract" ||
        after_legacy_hop.at("locals") != before_legacy_hop.at("locals") ||
        counter(after_legacy_hop, "local_solves") !=
            counter(before_legacy_hop, "local_solves") ||
        after_legacy_hop.at("pending_local_solves") != 0)
      throw std::runtime_error(
          "legacy chart did not remain usable while failing closed for value transport: " +
          json::serialize(legacy_hop));
    require_ok(request(json::object{
        {"schema", 2}, {"op", "tile.release"}, {"session", session},
        {"tile_plan", legacy_plan}}), "legacy-contract tile.release");

    const auto low_order_anchor = solve_local(
        session, anchor_chart, "0", "value-tail-order-probe-anchor", "2",
        0, true);
    const auto before = session_stats(session);
    const auto unsafe_tail = consume_value_hop(
        session, plan, "lower", low_order_anchor,
        "value-tail-order-probe-anchor",
        value_solver("-2/3", false, 4, "1/100", "1/10000"),
        "value-hop-success");
    require_ok(unsafe_tail, "unsafe-tail value hop");
    const auto after_unsafe_tail = session_stats(session);
    if (unsafe_tail.at("used") != false ||
        unsafe_tail.at("reason") !=
            "inflated-center-evaluation-fails-relative-accuracy-contract" ||
        counter(after_unsafe_tail, "local_solves") !=
            counter(before, "local_solves"))
      throw std::runtime_error(
          "low-order producer tail did not fail closed before solving: " +
          json::serialize(unsafe_tail));
    release_local(session, low_order_anchor);

    const auto expect_rejected_template = [&](json::object solver,
                                               const char* label,
                                               const char* detail) {
      const auto before_rejection = session_stats(session);
      const auto rejected = consume_value_hop(
          session, plan, "lower", anchor, "streaming-state-anchor",
          std::move(solver), "value-hop-invalid-template");
      const auto after_rejection = session_stats(session);
      if (rejected.at("status") != "error" ||
          std::string(rejected.at("detail").as_string()).find(detail) ==
              std::string::npos ||
          after_rejection.at("locals") != before_rejection.at("locals") ||
          counter(after_rejection, "local_solves") !=
              counter(before_rejection, "local_solves") ||
          after_rejection.at("pending_local_solves") != 0)
        throw std::runtime_error(
            std::string("regular value-hop accepted or retained state for ") +
            label + ": " + json::serialize(rejected) + " / " +
            json::serialize(after_rejection));
    };
    {
      auto solver = value_solver("-2/3", false, 4, "1/100", "1/10000");
      solver.at("run").as_object()["adaptive_probe"] = true;
      expect_rejected_template(std::move(solver), "adaptive probe",
                               "not one homogeneous (0,0,0) value run");
    }
    {
      auto solver = value_solver("-2/3", false, 4, "1/100", "1/10000");
      solver.at("run").as_object()["a_shift_min"] = -1;
      expect_rejected_template(std::move(solver), "nonzero shift origin",
                               "not one homogeneous (0,0,0) value run");
    }
    {
      auto solver = value_solver("-2/3", false, 4, "1/100", "1/10000");
      solver.at("run").as_object().at("a_shifts").as_array().emplace_back(
          "5");
      expect_rejected_template(std::move(solver), "extra a-shift",
                               "not one homogeneous (0,0,0) value run");
    }
    {
      auto solver = value_solver("-2/3", false, 4, "1/100", "1/10000");
      solver.at("run").as_object().at("a_shifts").as_array()[2] = "7";
      expect_rejected_template(std::move(solver), "non-Taylor a-shift",
                               "a-shifts are not the exact Taylor indices");
    }
    {
      auto solver = value_solver("-2/3", false, 4, "1/100", "1/10000");
      solver.at("run").as_object().at("schedule").as_array()[2]
          .as_array()[0].as_object()["da"] = "7";
      expect_rejected_template(std::move(solver), "non-Taylor da",
                               "Taylor by exact index");
    }
    {
      auto solver = value_solver("-2/3", false, 4, "1/100", "1/10000");
      solver.at("run").as_object().at("schedule").as_array()[2]
          .as_array()[0].as_object()["db"] = "1";
      expect_rejected_template(std::move(solver), "nonzero db",
                               "Taylor by exact index");
    }

    const auto logarithmic = solve_log_local(
        session, anchor_chart, "0", "value-hop-logarithmic-source");
    const auto before_log_rejection = session_stats(session);
    const auto log_rejection = consume_value_hop(
        session, plan, "lower", logarithmic,
        "value-hop-logarithmic-source",
        value_solver("-2/3", true, 4, "1/100", "1/10000"),
        "value-hop-nonsingle-valued");
    const auto after_log_rejection = session_stats(session);
    if (log_rejection.at("status") != "ok" ||
        log_rejection.at("used") != false ||
        log_rejection.at("reason") !=
            "incoming-local-has-no-certified-regular-tail-model" ||
        counter(after_log_rejection, "local_solves") !=
            counter(before_log_rejection, "local_solves") ||
        after_log_rejection.at("pending_local_solves") != 0)
      throw std::runtime_error(
          "non-single-valued prescribed source entered the value-hop path: " +
          json::serialize(log_rejection));
    release_local(session, logarithmic);

    const auto before_hops = session_stats(session);
    const auto lower = consume_value_hop(
        session, plan, "lower", anchor, "streaming-state-anchor",
        value_solver("-2/3", true, 30, "1/100", "1/10000"),
        "value-hop-success");
    const auto upper = consume_value_hop(
        session, plan, "upper", anchor, "streaming-state-anchor",
        value_solver("2/3", true, 30, "1/100", "1/10000"),
        "value-hop-success");
    require_ok(lower, "lower regular value hop");
    require_ok(upper, "upper regular value hop");
    const auto lower_local_stats = request(json::object{
        {"schema", 2}, {"op", "local.stats"}, {"session", session},
        {"local", lower.at("next_local").as_object().at("local")}});
    require_ok(lower_local_stats, "lower regular value local.stats");
    const auto& compact_derivation =
        lower_local_stats.at("retained_derivation").as_object();
    const auto lower_stats_json = json::serialize(lower_local_stats);
    if (compact_derivation.size() != 5 ||
        compact_derivation.at("schema") !=
            "diffexp2-retained-plan-value-handoff-v2" ||
        compact_derivation.at("capability") !=
            "retained-native-regular-value-handoff-v2" ||
        compact_derivation.at("checkpoint_identity") !=
            "value-hop-success:lower:local:1" ||
        compact_derivation.at("ownership") !=
            "sealed-plan-match-lineage" ||
        compact_derivation.at("provenance").as_object()
                .at("algorithm") != "fnv1a64-v1" ||
        std::string(compact_derivation.at("provenance").as_object()
                        .at("fingerprint").as_string()).empty() ||
        lower_stats_json.size() >= 16384 ||
        lower_stats_json.find("\"source_model\"") != std::string::npos ||
        lower_stats_json.find("\"inflation\"") != std::string::npos)
      throw std::runtime_error(
          "value-hop local.stats exported a non-compact derivation: " +
          lower_stats_json);
    const auto after_hops = session_stats(session);
    if (!lower.at("used").as_bool() || !upper.at("used").as_bool() ||
        lower.at("value_hops") != 1 || lower.at("basis_matches") != 0 ||
        upper.at("value_hops") != 1 || upper.at("basis_matches") != 0 ||
        counter(after_hops, "local_matches") !=
            counter(before_hops, "local_matches") ||
        counter(after_hops, "local_solves") !=
            counter(before_hops, "local_solves") + 2)
      throw std::runtime_error(
          "eligible regular hops did not use exactly one value solve each: " +
          json::serialize(lower) + " / " + json::serialize(upper));

    // Compare the optimized local directly with the established one-column
    // basis/match path before either becomes dependency-only state.
    const auto comparison_anchor = solve_local(
        session, anchor_chart, "0", "value-basis-comparison-anchor", "2",
        4, true);
    const auto comparison_basis = solve_local(
        session, lower_chart, "-2/3", "value-basis-comparison-column", "1",
        4, true);
    const auto comparison_planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "value-basis-comparison-plan"},
        {"division_order", 3},
        {"lower", arm("-2/3", anchor_chart, lower_chart)},
        {"upper", arm("2/3", anchor_chart, upper_chart)}});
    require_ok(comparison_planned, "value/basis comparison tile.plan");
    const auto comparison_plan = std::string(
        comparison_planned.at("tile_plan").as_string());
    const auto basis_hop = consume_hop(
        session, comparison_plan, "lower", comparison_anchor,
        "value-basis-comparison-anchor", comparison_basis,
        "value-basis-comparison", "value-basis-comparison-plan");
    require_ok(basis_hop, "value/basis comparison hop");
    const auto evaluate_center = [&](const std::string& session_handle,
                                     const json::value& local) {
      const auto evaluated = request(json::object{
          {"schema", 2}, {"op", "local.evaluate"},
          {"session", session_handle},
          {"local", local},
          {"point", json::object{{"exact", "0"}}},
          {"options", json::object{{"tail_estimate", false}}},
          {"output_digits", 40}});
      require_ok(evaluated, "value/basis center evaluation");
      return evaluated.at("value");
    };
    if (!epsilon_vectors_overlap(
            evaluate_center(
                session,
                lower.at("next_local").as_object().at("local")),
            evaluate_center(
                session,
                basis_hop.at("next_local").as_object().at("local"))))
      throw std::runtime_error(
          "eligible value hop differs from the ordinary basis/match result");
    release_local(session, comparison_anchor);
    release_local(session, std::string(
        basis_hop.at("next_local").as_object().at("local").as_string()));
    require_ok(request(json::object{
        {"schema", 2}, {"op", "tile.release"}, {"session", session},
        {"tile_plan", comparison_plan}}),
        "value/basis comparison tile.release");

    const auto published = publish_consumed_states(
        session, plan, anchor, lower, upper, "value-hop-success");
    require_ok(published, "value-hop transport.publish_consumed_states");
    const auto& states = published.at("states").as_object();
    const auto lower_state = states.at("lower").as_object();
    const auto upper_state = states.at("upper").as_object();
    if (lower_state.at("value_hops") != 1 ||
        lower_state.at("basis_matches") != 0 ||
        upper_state.at("value_hops") != 1 ||
        upper_state.at("basis_matches") != 0)
      throw std::runtime_error(
          "published value-hop instrumentation is inconsistent");
    const auto lower_value = contract_and_export(
        session, lower_state, "value-hop-before-checkpoint", 31);
    require_ok(request(json::object{
        {"schema", 2}, {"op", "tile.release"}, {"session", session},
        {"tile_plan", plan}}), "value-hop tile.release");
    release_local(session, anchor);
    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", session}, {"path", checkpoint},
        {"checkpoint_identity", "value-hop-roundtrip"}});
    require_ok(saved, "value-hop checkpoint.save");
    const auto container = diffexp2::checkpoint::read(checkpoint);
    auto tampered_payload =
        json::parse(container.payload_json).as_object();
    bool tampered_retained_value = false;
    for (auto& raw_local :
         tampered_payload.at("retained_locals").as_array()) {
      auto& local = raw_local.as_object();
      if (local.at("retained_derivation").is_null()) continue;
      auto& derivation = local.at("retained_derivation").as_object();
      if (derivation.at("schema") !=
          "diffexp2-retained-plan-value-handoff-v2")
        continue;
      auto& coefficients = derivation.at("tail_contract").as_object()
          .at("inflation").as_object()
          .at("retained_value").as_object()
          .at("coefficients").as_array();
      if (coefficients.size() < 2 || coefficients[0] == coefficients[1])
        throw std::runtime_error(
            "value-hop checkpoint has no distinct retained coefficient for tamper test");
      coefficients[0] = coefficients[1];
      tampered_retained_value = true;
      break;
    }
    if (!tampered_retained_value)
      throw std::runtime_error(
          "value-hop checkpoint did not contain a v2 retained value");
    diffexp2::checkpoint::write_atomic(
        tampered_checkpoint, container.header_json,
        json::serialize(tampered_payload));
    const auto tampered_restore = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", tampered_checkpoint},
        {"expected_identity", "value-hop-roundtrip"}});
    if (tampered_restore.at("status") == "ok") {
      (void)request(json::object{
          {"schema", 2}, {"op", "session.close"},
          {"session", tampered_restore.at("session")}});
      throw std::runtime_error(
          "checkpoint restore accepted a tampered certified retained value");
    }
    if (std::string(tampered_restore.at("detail").as_string()).find(
            "differs from its restored incoming local") ==
        std::string::npos)
      throw std::runtime_error(
          "retained-value tamper was rejected for the wrong reason: " +
          json::serialize(tampered_restore));
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"}, {"session", session}}),
        "value-hop session.close");
    session.clear();
    const auto restored_record = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint},
        {"expected_identity", "value-hop-roundtrip"}});
    require_ok(restored_record, "value-hop checkpoint.restore");
    restored = std::string(restored_record.at("session").as_string());
    const auto restored_lower = request(json::object{
        {"schema", 2}, {"op", "transport.stats"},
        {"session", restored},
        {"transport_state", lower_state.at("transport_state")}});
    require_ok(restored_lower, "restored value-hop transport.stats");
    if (restored_record.at("transport_states").as_array().size() != 2 ||
        restored_record.at("locals").as_array().size() != 0 ||
        restored_lower.at("value_hops") != 1 ||
        restored_lower.at("basis_matches") != 0 ||
        !epsilon_vectors_overlap(
            contract_and_export(restored, lower_state,
                                "value-hop-after-checkpoint", 31),
            lower_value))
      throw std::runtime_error(
          "regular value-hop checkpoint roundtrip changed its sealed state");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"}, {"session", restored}}),
        "restored value-hop session.close");
    restored.clear();
    std::remove(checkpoint.c_str());
    std::remove(tampered_checkpoint.c_str());
  } catch (...) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    if (!restored.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored}});
    std::remove(checkpoint.c_str());
    std::remove(tampered_checkpoint.c_str());
    throw;
  }
}

void test_acb_prescribed_value_hops_match_basis() {
  for (const std::int32_t sign : {-1, 1}) {
    std::string session;
    try {
      const auto suffix = sign < 0 ? std::string("minus")
                                   : std::string("plus");
      const auto created = request(json::object{
          {"schema", 2}, {"op", "session.create"},
          {"domain", "acb"}, {"precision_bits", 256},
          {"output_digits", 40}, {"chart_capacity", 6},
          {"local_capacity", 20}, {"match_capacity", 6},
          {"tile_plan_capacity", 2}, {"transport_state_capacity", 2},
          {"line_result_capacity", 2}});
      require_ok(created, "prescribed Acb session.create");
      session = std::string(created.at("session").as_string());
      const auto anchor_chart = prepare_chart(
          session, "prescribed-acb-anchor-" + suffix, "0", true,
          "acb", 4, sign);
      const auto lower_chart = prepare_chart(
          session, "prescribed-acb-lower-" + suffix, "-2/3", true,
          "acb", 4, sign);
      const auto upper_chart = prepare_chart(
          session, "prescribed-acb-upper-" + suffix, "2/3", true,
          "acb", 4, sign);
      const auto value_anchor = solve_local(
          session, anchor_chart, "0",
          "prescribed-acb-value-anchor-" + suffix, "2", 10, true, sign);
      const auto value_plan_record = request(json::object{
          {"schema", 2}, {"op", "tile.plan"}, {"session", session},
          {"checkpoint_identity", "prescribed-acb-value-plan-" + suffix},
          {"division_order", 3},
          {"lower", arm("-2/3", anchor_chart, lower_chart, true, sign)},
          {"upper", arm("2/3", anchor_chart, upper_chart, true, sign)}});
      require_ok(value_plan_record, "prescribed Acb value tile.plan");
      const auto value_plan = std::string(
          value_plan_record.at("tile_plan").as_string());
      if (sign < 0) {
        const auto before_loose_threshold = session_stats(session);
        const auto loose_threshold = consume_value_hop(
            session, value_plan, "lower", value_anchor,
            "prescribed-acb-value-anchor-" + suffix,
            value_solver("-2/3", true, 10, "9/10", "1/10000", sign),
            "prescribed-acb-loose-threshold",
            "prescribed-acb-value-plan-" + suffix);
        const auto after_loose_threshold = session_stats(session);
        if (loose_threshold.at("status") != "error" ||
            std::string(loose_threshold.at("detail").as_string()).find(
                "relative-accuracy threshold differs from its prepared owner") ==
                std::string::npos ||
            after_loose_threshold.at("locals") !=
                before_loose_threshold.at("locals") ||
            after_loose_threshold.at("matches") !=
                before_loose_threshold.at("matches") ||
            after_loose_threshold.at("transport_states") !=
                before_loose_threshold.at("transport_states") ||
            counter(after_loose_threshold, "local_solves") !=
                counter(before_loose_threshold, "local_solves") ||
            after_loose_threshold.at("pending_local_solves") != 0)
          throw std::runtime_error(
              "looser caller-supplied Acb accuracy threshold mutated retained state: " +
              json::serialize(loose_threshold) + " / " +
              json::serialize(after_loose_threshold));
      }
      const auto lower_value = consume_value_hop(
          session, value_plan, "lower", value_anchor,
          "prescribed-acb-value-anchor-" + suffix,
          value_solver("-2/3", true, 10, "1/100", "1/10000", sign),
          "prescribed-acb-value-" + suffix,
          "prescribed-acb-value-plan-" + suffix);
      const auto upper_value = consume_value_hop(
          session, value_plan, "upper", value_anchor,
          "prescribed-acb-value-anchor-" + suffix,
          value_solver("2/3", true, 10, "1/100", "1/10000", sign),
          "prescribed-acb-value-" + suffix,
          "prescribed-acb-value-plan-" + suffix);
      require_ok(lower_value, "prescribed Acb lower value hop");
      require_ok(upper_value, "prescribed Acb upper value hop");
      if (lower_value.at("used") != true || upper_value.at("used") != true)
        throw std::runtime_error(
            "prescribed Acb value hop did not accept both local signs");

      const auto basis_anchor = solve_local(
          session, anchor_chart, "0",
          "prescribed-acb-basis-anchor-" + suffix, "2", 10, true, sign);
      const auto lower_basis = solve_local(
          session, lower_chart, "-2/3",
          "prescribed-acb-lower-basis-" + suffix, "1", 10, true, sign);
      const auto upper_basis = solve_local(
          session, upper_chart, "2/3",
          "prescribed-acb-upper-basis-" + suffix, "1", 10, true, sign);
      const auto basis_plan_record = request(json::object{
          {"schema", 2}, {"op", "tile.plan"}, {"session", session},
          {"checkpoint_identity", "prescribed-acb-basis-plan-" + suffix},
          {"division_order", 3},
          {"lower", arm("-2/3", anchor_chart, lower_chart, true, sign)},
          {"upper", arm("2/3", anchor_chart, upper_chart, true, sign)}});
      require_ok(basis_plan_record, "prescribed Acb basis tile.plan");
      const auto basis_plan = std::string(
          basis_plan_record.at("tile_plan").as_string());
      const auto lower_basis_hop = consume_basis_hop(
          session, basis_plan, "prescribed-acb-basis-plan-" + suffix,
          "lower", basis_anchor,
          "prescribed-acb-basis-anchor-" + suffix, {lower_basis},
          "prescribed-acb-basis-" + suffix);
      const auto upper_basis_hop = consume_basis_hop(
          session, basis_plan, "prescribed-acb-basis-plan-" + suffix,
          "upper", basis_anchor,
          "prescribed-acb-basis-anchor-" + suffix, {upper_basis},
          "prescribed-acb-basis-" + suffix);
      require_ok(lower_basis_hop, "prescribed Acb lower basis hop");
      require_ok(upper_basis_hop, "prescribed Acb upper basis hop");

      const auto evaluate_center = [&](const json::value& local) {
        const auto evaluated = request(json::object{
            {"schema", 2}, {"op", "local.evaluate"},
            {"session", session}, {"local", local},
            {"point", json::object{{"exact", "0"}}},
            {"options", json::object{{"tail_estimate", false}}},
            {"output_digits", 40}});
        require_ok(evaluated, "prescribed Acb center evaluation");
        return evaluated.at("value");
      };
      const auto& lower_value_local =
          lower_value.at("next_local").as_object().at("local");
      const auto& upper_value_local =
          upper_value.at("next_local").as_object().at("local");
      const auto& lower_basis_local =
          lower_basis_hop.at("next_local").as_object().at("local");
      const auto& upper_basis_local =
          upper_basis_hop.at("next_local").as_object().at("local");
      if (!epsilon_vectors_overlap(evaluate_center(lower_value_local),
                                   evaluate_center(lower_basis_local)) ||
          !epsilon_vectors_overlap(evaluate_center(upper_value_local),
                                   evaluate_center(upper_basis_local)))
        throw std::runtime_error(
            "prescribed Acb value and ordinary basis hops disagree");
      for (const auto* local : {&lower_value_local, &upper_value_local}) {
        const auto stats = request(json::object{
            {"schema", 2}, {"op", "local.stats"},
            {"session", session}, {"local", *local}});
        require_ok(stats, "prescribed Acb value local.stats");
        if (stats.at("metadata").as_object().at("prescriptions") !=
            prescriptions(true, sign))
          throw std::runtime_error(
              "prescribed Acb value hop did not preserve receiver prescriptions");
      }
      require_ok(request(json::object{
          {"schema", 2}, {"op", "session.close"}, {"session", session}}),
          "prescribed Acb session.close");
      session.clear();
    } catch (...) {
      if (!session.empty())
        (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                   {"session", session}});
      throw;
    }
  }
}

void test_acb_value_handoff_significance_gate() {
  std::string session;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "acb"}, {"precision_bits", 256},
        {"output_digits", 40}, {"chart_capacity", 6},
        {"local_capacity", 12}, {"match_capacity", 4},
        {"tile_plan_capacity", 2}, {"transport_state_capacity", 2},
        {"line_result_capacity", 2}});
    require_ok(created, "Acb value-gate session.create");
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_chart(
        session, "acb-value-anchor-chart", "0", false, "acb", 4);
    const auto lower_chart = prepare_chart(
        session, "acb-value-lower-chart", "-2/3", false, "acb", 4);
    const auto upper_chart = prepare_chart(
        session, "acb-value-upper-chart", "2/3", false, "acb", 4);
    const auto tiny_crossing_zero = solve_local(
        session, anchor_chart, "0", "streaming-state-anchor",
        "[1e-6 +/- 1e-3]", 10, false);
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "streaming-state-plan"},
        {"division_order", 3},
        {"lower", arm("-2/3", anchor_chart, lower_chart, false)},
        {"upper", arm("2/3", anchor_chart, upper_chart, false)}});
    require_ok(planned, "Acb value-gate tile.plan");
    const auto plan = std::string(planned.at("tile_plan").as_string());
    const auto accepted = consume_value_hop(
        session, plan, "lower", tiny_crossing_zero,
        "streaming-state-anchor",
        value_solver("-2/3", false, 10, "1/100", "1/10000"),
        "acb-value-gate");
    require_ok(accepted, "Acb tiny value hop");
    if (accepted.at("used") != true)
      throw std::runtime_error(
          "scale-aware Acb gate rejected a tiny, absolutely accurate ball: " +
          json::serialize(accepted));

    const auto inaccurate = solve_local(
        session, anchor_chart, "0", "acb-inaccurate-anchor",
        "[1e-6 +/- 1e-1]", 10, false);
    const auto before_rejection = session_stats(session);
    const auto rejected = consume_value_hop(
        session, plan, "upper", inaccurate, "acb-inaccurate-anchor",
        value_solver("2/3", false, 10, "1/100", "1/10000"),
        "acb-value-gate");
    require_ok(rejected, "Acb inaccurate value preflight");
    const auto after_rejection = session_stats(session);
    if (rejected.at("used") != false ||
        rejected.at("reason") !=
            "inflated-center-evaluation-fails-relative-accuracy-contract" ||
        counter(after_rejection, "local_solves") !=
            counter(before_rejection, "local_solves"))
      throw std::runtime_error(
          "Acb significance rejection did not fail closed before solving: " +
          json::serialize(rejected));
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"}, {"session", session}}),
        "Acb value-gate session.close");
    session.clear();
  } catch (...) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    throw;
  }
}

void test_acb_consuming_hop_reservoir_retry_lifecycle() {
  std::string session;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "acb"}, {"precision_bits", 256},
        {"output_digits", 40}, {"chart_capacity", 3},
        {"local_capacity", 4}, {"match_capacity", 2},
        {"tile_plan_capacity", 1}, {"transport_state_capacity", 1},
        {"line_result_capacity", 1}});
    require_ok(created, "reservoir-retry session.create");
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_chart(
        session, "reservoir-retry-anchor-chart", "0", false, "acb");
    const auto lower_chart = prepare_chart(
        session, "reservoir-retry-lower-chart", "-2/3", false, "acb");
    const auto upper_chart = prepare_chart(
        session, "reservoir-retry-upper-chart", "2/3", false, "acb");
    const auto anchor = solve_local(
        session, anchor_chart, "0", "reservoir-retry-anchor", "2", 0,
        false);
    const auto short_basis = solve_local(
        session, lower_chart, "-2/3", "reservoir-retry-short", "1", 0,
        false, -1, 0);
    const auto full_basis = solve_local(
        session, lower_chart, "-2/3", "reservoir-retry-full", "1", 0,
        false);
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "streaming-state-plan"},
        {"division_order", 3},
        {"lower", arm("-2/3", anchor_chart, lower_chart, false)},
        {"upper", arm("2/3", anchor_chart, upper_chart, false)}});
    require_ok(planned, "reservoir-retry tile.plan");
    const auto plan = std::string(planned.at("tile_plan").as_string());

    const auto before = session_stats(session);
    const auto retry = consume_hop(
        session, plan, "lower", anchor, "reservoir-retry-anchor",
        short_basis, "reservoir-retry");
    const auto after = session_stats(session);
    if (retry.at("status") != "error" ||
        retry.at("reason") != "acb_match_residual_inconclusive" ||
        retry.at("retryable_epsilon_reservoir") != true ||
        retry.at("retryable_matching_clearance") != false ||
        retry.at("required_additional_epsilon_orders") != 2 ||
        after.at("locals") != before.at("locals") ||
        after.at("pending_local_solves") != 0 ||
        counter(after, "local_solves") != counter(before, "local_solves"))
      throw std::runtime_error(
          "Acb reservoir retry was not structured and transactional: " +
          json::serialize(retry) + " / " + json::serialize(after));

    release_local(session, short_basis);
    const auto recovered = consume_hop(
        session, plan, "lower", anchor, "reservoir-retry-anchor",
        full_basis, "reservoir-retry-recovered");
    require_ok(recovered, "reservoir-retry recovered consume_hop");
    if (session_stats(session).at("pending_local_solves") != 0)
      throw std::runtime_error(
          "Acb reservoir retry did not leave reusable local capacity");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"}, {"session", session}}),
        "reservoir-retry session.close");
    session.clear();
  } catch (...) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    throw;
  }
}

void test_multiblock_regular_value_fallback_owner() {
  std::string session;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "rational"}, {"output_digits", 40},
        {"chart_capacity", 6}, {"local_capacity", 16},
        {"match_capacity", 6}, {"tile_plan_capacity", 2},
        {"transport_state_capacity", 2}, {"line_result_capacity", 2}});
    require_ok(created, "multiblock fallback session.create");
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_multiblock_chart(
        session, "multiblock-anchor", "0");
    const auto lower_chart = prepare_multiblock_chart(
        session, "multiblock-lower", "-2/3");
    const auto upper_chart = prepare_multiblock_chart(
        session, "multiblock-upper", "2/3");
    const auto incoming = solve_multiblock_local(
        session, anchor_chart, "0", "multiblock-anchor-local", {"2", "3"},
        0);
    const auto basis_0 = solve_multiblock_local(
        session, lower_chart, "-2/3", "multiblock-basis-0", {"1", "0"});
    const auto basis_1 = solve_multiblock_local(
        session, lower_chart, "-2/3", "multiblock-basis-1", {"0", "1"});
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "multiblock-plan"},
        {"division_order", 3},
        {"lower", arm("-2/3", anchor_chart, lower_chart)},
        {"upper", arm("2/3", anchor_chart, upper_chart)}});
    require_ok(planned, "multiblock fallback tile.plan");
    const auto plan = std::string(planned.at("tile_plan").as_string());
    const auto before_preflight = session_stats(session);
    const auto ineligible = consume_value_hop(
        session, plan, "lower", incoming, "multiblock-anchor-local",
        value_solver("-2/3", true, 0), "multiblock-fallback",
        "multiblock-plan");
    require_ok(ineligible, "multiblock value preflight");
    const auto after_preflight = session_stats(session);
    if (ineligible.at("used") != false ||
        ineligible.at("reason") !=
            "rational-value-handoff-has-no-exact-polynomial-tail-zero-certificate" ||
        counter(after_preflight, "local_solves") !=
            counter(before_preflight, "local_solves"))
      throw std::runtime_error(
          "multiblock regular value hop did not become an ineligible fallback");
    const auto fallback = consume_basis_hop(
        session, plan, "multiblock-plan", "lower", incoming,
        "multiblock-anchor-local", {basis_0, basis_1},
        "multiblock-fallback");
    require_ok(fallback, "multiblock primitive-owner fallback");
    const auto& next = fallback.at("next_local").as_object();
    const auto& operator_reference =
        next.at("source_operator_reference").as_object();
    if (std::string(next.at("chart").as_string()) != lower_chart ||
        next.if_contains("source_operator_identity") != nullptr ||
        operator_reference.at("algorithm") != "fnv1a64-v1" ||
        counter(operator_reference, "identity_bytes") !=
            std::string("de2-operator-multiblock-lower").size())
      throw std::runtime_error(
          "multiblock fallback did not retain the primitive receiving owner: " +
          json::serialize(fallback));
    const auto evaluate = [&](const std::string& local,
                              const std::string& point) {
      const auto result = request(json::object{
          {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
          {"local", local}, {"point", json::object{{"exact", point}}},
          {"options", json::object{{"tail_estimate", false}}},
          {"output_digits", 40}});
      require_ok(result, "multiblock fallback evaluation");
      return result.at("value");
    };
    if (evaluate(incoming, "-2/3") !=
        evaluate(std::string(next.at("local").as_string()), "0"))
      throw std::runtime_error(
          "multiblock primitive-owner fallback changed the ordinary basis result");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"}, {"session", session}}),
        "multiblock fallback session.close");
    session.clear();
  } catch (...) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    throw;
  }
}

void test_streaming_consumed_transport() {
  const std::string checkpoint =
      "/tmp/diffexp2-streaming-consumed-roundtrip.de2cp";
  std::remove(checkpoint.c_str());
  std::string session;
  std::string restored;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "rational"}, {"output_digits", 40},
        {"chart_capacity", 6}, {"local_capacity", 12},
        {"match_capacity", 4}, {"tile_plan_capacity", 2},
        {"transport_state_capacity", 4}, {"line_result_capacity", 4}});
    require_ok(created, "streaming session.create");
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_chart(
        session, "streaming-anchor-chart", "0");
    const auto lower_chart = prepare_chart(
        session, "streaming-lower-chart", "-2/3");
    const auto upper_chart = prepare_chart(
        session, "streaming-upper-chart", "2/3");
    const auto anchor = solve_local(
        session, anchor_chart, "0", "streaming-state-anchor", "2");
    const auto lower_basis = solve_local(
        session, lower_chart, "-2/3", "streaming-lower-basis", "1");
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "streaming-state-plan"},
        {"division_order", 3},
        {"lower", arm("-2/3", anchor_chart, lower_chart)},
        {"upper", arm("2/3", anchor_chart, upper_chart)}});
    require_ok(planned, "streaming tile.plan");
    const auto plan = std::string(planned.at("tile_plan").as_string());

    const auto lower_before = session_stats(session);
    const auto ineligible_value = consume_value_hop(
        session, plan, "lower", anchor, "streaming-state-anchor",
        value_solver("-2/3", true, 0), "streaming-state-success");
    require_ok(ineligible_value, "unsafe-tail value-hop preflight");
    const auto after_ineligible_value = session_stats(session);
    if (ineligible_value.at("used") != false ||
        ineligible_value.at("reason") !=
            "rational-value-handoff-has-no-exact-polynomial-tail-zero-certificate" ||
        after_ineligible_value.at("locals") != lower_before.at("locals") ||
        counter(after_ineligible_value, "local_solves") !=
            counter(lower_before, "local_solves"))
      throw std::runtime_error(
          "unsafe-tail value hop did not fail closed before solving: " +
          json::serialize(ineligible_value));
    const auto lower = consume_hop(
        session, plan, "lower", anchor, "streaming-state-anchor",
        lower_basis, "streaming-state-success");
    require_ok(lower, "streaming lower consume_hop");
    const auto lower_after = session_stats(session);
    if (lower_after.at("locals") != 2 ||
        lower_after.at("matches") != 0 ||
        lower.at("consumed_basis_handles").as_array().size() != 1 ||
        counter(lower_after, "local_coefficient_count") >
            counter(lower_before, "local_coefficient_count"))
      throw std::runtime_error(
          "streaming lower hop retained its consumed basis slab: " +
          json::serialize(lower_after));
    if (request(json::object{
            {"schema", 2}, {"op", "local.stats"}, {"session", session},
            {"local", lower_basis}}).at("status") != "error")
      throw std::runtime_error(
          "streaming lower basis remained publicly visible");

    // The upper basis is not solved until the lower basis has been consumed.
    const auto upper_basis = solve_local(
        session, upper_chart, "2/3", "streaming-upper-basis", "1");
    const auto upper_ready = session_stats(session);
    if (upper_ready.at("locals") != 3 ||
        counter(upper_ready, "local_coefficient_count") <
            counter(lower_after, "local_coefficient_count"))
      throw std::runtime_error(
          "streaming fixture did not admit exactly one next receiving basis");
    const auto upper = consume_hop(
        session, plan, "upper", anchor, "streaming-state-anchor",
        upper_basis, "streaming-state-success");
    require_ok(upper, "streaming upper consume_hop");
    const auto upper_after = session_stats(session);
    if (upper_after.at("locals") != 3 ||
        upper_after.at("matches") != 0 ||
        request(json::object{
            {"schema", 2}, {"op", "local.stats"}, {"session", session},
            {"local", upper_basis}}).at("status") != "error")
      throw std::runtime_error(
          "streaming upper hop retained its consumed basis slab");

    const auto published = publish_consumed_states(
        session, plan, anchor, lower, upper, "streaming-state-success");
    require_ok(published, "transport.publish_consumed_states");
    const auto after_publish = session_stats(session);
    if (after_publish.at("locals") != 1 ||
        after_publish.at("transport_states") != 2 ||
        counter(after_publish, "transport_arm_marches") !=
            counter(upper_after, "transport_arm_marches") + 2 ||
        published.at("consumed_tile_local_handles").as_array().size() != 2)
      throw std::runtime_error(
          "streaming state publication did not consume its tile-local tokens: " +
          json::serialize(after_publish));
    const auto& states = published.at("states").as_object();
    const auto lower_state = states.at("lower").as_object();
    const auto upper_state = states.at("upper").as_object();
    if (lower_state.at("value_hops") != 0 ||
        lower_state.at("basis_matches") != 1 ||
        upper_state.at("value_hops") != 0 ||
        upper_state.at("basis_matches") != 1)
      throw std::runtime_error(
          "basis-fallback instrumentation is inconsistent");
    const auto lower_value = contract_and_export(
        session, lower_state, "streaming-before-checkpoint");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "tile.release"}, {"session", session},
        {"tile_plan", plan}}), "streaming tile.release");
    release_local(session, anchor);
    require_ok(request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
        {"path", checkpoint},
        {"checkpoint_identity", "streaming-roundtrip"}}),
        "streaming checkpoint.save");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"}, {"session", session}}),
        "streaming session.close");
    session.clear();
    const auto restored_record = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"}, {"path", checkpoint},
        {"expected_identity", "streaming-roundtrip"}});
    require_ok(restored_record, "streaming checkpoint.restore");
    restored = std::string(restored_record.at("session").as_string());
    if (restored_record.at("transport_states").as_array().size() != 2 ||
        restored_record.at("planned_match_hops").as_array().size() != 0 ||
        contract_and_export(restored, lower_state,
                            "streaming-after-checkpoint") != lower_value)
      throw std::runtime_error(
          "streaming compact checkpoint roundtrip changed contraction output");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"}, {"session", restored}}),
        "streaming restored session.close");
    restored.clear();
    std::remove(checkpoint.c_str());
  } catch (...) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    if (!restored.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored}});
    std::remove(checkpoint.c_str());
    throw;
  }
}

void test_consuming_transport() {
  const std::string checkpoint =
      "/tmp/diffexp2-consuming-transport-roundtrip.de2cp";
  std::remove(checkpoint.c_str());
  std::string session;
  std::string restored;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "rational"}, {"output_digits", 40},
        {"chart_capacity", 6}, {"local_capacity", 24},
        {"match_capacity", 12}, {"tile_plan_capacity", 2},
        {"transport_state_capacity", 8}, {"line_result_capacity", 8}});
    require_ok(created, "consuming session.create");
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_chart(
        session, "consuming-anchor-chart", "0");
    const auto lower_chart = prepare_chart(
        session, "consuming-lower-chart", "-2/3");
    const auto upper_chart = prepare_chart(
        session, "consuming-upper-chart", "2/3");
    const auto anchor = solve_local(
        session, anchor_chart, "0", "consuming-state-anchor", "2");
    const auto lower_basis = solve_local(
        session, lower_chart, "-2/3", "consuming-lower-basis", "1");
    const auto upper_basis = solve_local(
        session, upper_chart, "2/3", "consuming-upper-basis", "1");
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "consuming-state-plan"},
        {"division_order", 3},
        {"lower", arm("-2/3", anchor_chart, lower_chart)},
        {"upper", arm("2/3", anchor_chart, upper_chart)}});
    require_ok(planned, "consuming tile.plan");
    const auto plan = std::string(planned.at("tile_plan").as_string());
    const auto before = session_stats(session);
    const auto consumed = run_consuming_pair(
        session, plan, anchor, lower_basis, upper_basis,
        "consuming-state-success");
    require_ok(consumed, "transport.run_arms_consuming");
    const auto after = session_stats(session);
    const auto& diagnostics = consumed.at("consumption").as_array();
    if (diagnostics.size() != 2 ||
        diagnostics[0].as_object().at("session_locals_after_hop") != 2 ||
        diagnostics[1].as_object().at("session_locals_after_hop") != 1 ||
        counter(diagnostics[0].as_object(),
                "session_local_coefficient_count_after_hop") >=
            counter(before, "local_coefficient_count") ||
        counter(diagnostics[1].as_object(),
                "session_local_coefficient_count_after_hop") >=
            counter(diagnostics[0].as_object(),
                    "session_local_coefficient_count_after_hop") ||
        after.at("locals") != 1 ||
        consumed.at("consuming_basis_handles") != true ||
        consumed.at("workers") != 1)
      throw std::runtime_error(
          "consuming transport did not monotonically release basis slabs: " +
          json::serialize(consumed) + " / " + json::serialize(after));
    for (const auto& handle : {lower_basis, upper_basis}) {
      const auto released = request(json::object{
          {"schema", 2}, {"op", "local.stats"}, {"session", session},
          {"local", handle}});
      if (released.at("status") != "error")
        throw std::runtime_error(
            "consumed basis handle remained publicly visible");
    }
    const auto& states = consumed.at("states").as_object();
    const auto lower_state = states.at("lower").as_object();
    const auto lower_value = contract_and_export(
        session, lower_state, "consuming-before-checkpoint");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "tile.release"}, {"session", session},
        {"tile_plan", plan}}), "consuming tile.release");
    release_local(session, anchor);
    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", session}, {"path", checkpoint},
        {"checkpoint_identity", "consuming-roundtrip"}});
    require_ok(saved, "consuming checkpoint.save");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"}, {"session", session}}),
        "consuming session.close");
    session.clear();
    const auto restored_record = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint},
        {"expected_identity", "consuming-roundtrip"}});
    require_ok(restored_record, "consuming checkpoint.restore");
    restored = std::string(restored_record.at("session").as_string());
    if (restored_record.at("transport_states").as_array().size() != 2 ||
        restored_record.at("planned_match_hops").as_array().size() != 0 ||
        contract_and_export(restored, lower_state,
                            "consuming-after-checkpoint") != lower_value)
      throw std::runtime_error(
          "consuming compact checkpoint roundtrip changed retained state");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"}, {"session", restored}}),
        "consuming restored session.close");
    restored.clear();
    std::remove(checkpoint.c_str());
  } catch (...) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    if (!restored.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored}});
    std::remove(checkpoint.c_str());
    throw;
  }
}

void test_frame_independent_physical_value_wiring() {
  const std::string checkpoint =
      "/tmp/diffexp2-frame-independent-physical-hop.de2cp";
  const std::string resaved_checkpoint =
      "/tmp/diffexp2-frame-independent-physical-hop-resaved.de2cp";
  std::remove(checkpoint.c_str());
  std::remove(resaved_checkpoint.c_str());
  std::string session;
  std::string restored_session;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "acb"}, {"precision_bits", 256},
        {"output_digits", 40}, {"chart_capacity", 8},
        {"local_capacity", 16}, {"match_capacity", 8},
        {"tile_plan_capacity", 4}});
    require_ok(created, "physical wiring session.create");
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_chart(
        session, "physical-wiring-anchor", "0", true, "acb", 0);
    const auto first_equation = prepare_regular_equation_owner(
        session, "physical-wiring-first", "-2/3");
    const auto second_equation = prepare_regular_equation_owner(
        session, "physical-wiring-second", "-4/3");
    const auto second_framed_chart = prepare_chart(
        session, "physical-wiring-second", "-4/3", true, "acb", 0);
    const auto anchor = solve_local(
        session, anchor_chart, "0", "physical-wiring-anchor-local", "2",
        30, true);

    const std::string plan_checkpoint = "physical-wiring-plan";
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan_arm"}, {"session", session},
        {"checkpoint_identity", plan_checkpoint}, {"division_order", 3},
        {"arm", json::object{
             {"from_exact", "0"}, {"to_exact", "-4/3"},
             {"charts", json::array{
                  anchor_chart, first_equation, second_equation}},
             {"topology", topology(true)}}}});
    require_ok(planned, "physical wiring tile.plan_arm");
    const auto plan = std::string(planned.at("tile_plan").as_string());

    const auto oversized = consume_physical_value_hop(
        session, plan, "lower", 0, anchor,
        "physical-wiring-anchor-local",
        physical_value_solver("-2/3", 12), "physical-wiring-oversized",
        plan_checkpoint, std::numeric_limits<std::int32_t>::min());
    if (oversized.at("status") != "error" ||
        std::string(oversized.at("detail").as_string()).find(
            "unreasonably large") == std::string::npos)
      throw std::runtime_error(
          "ordinary physical hop did not reject an oversized zero-padding request before allocation: " +
          json::serialize(oversized));

    const auto direct = consume_physical_value_hop(
        session, plan, "lower", 0, anchor,
        "physical-wiring-anchor-local",
        physical_value_solver("-2/3", 12), "physical-wiring-hop",
        plan_checkpoint, -1);
    require_ok(direct, "physical wiring direct hop");
    if (direct.at("used") != true ||
        direct.at("execution_mode") !=
            "causal-ordinary-physical-evolution" ||
        direct.at("output_tail_status") != "unsupported" ||
        direct.at("next_hop_policy") != "exact-framed-fallback")
      throw std::runtime_error(
          "ordinary physical hop did not advertise its exact tail/fallback contract: " +
          json::serialize(direct));
    const auto& direct_local = direct.at("next_local").as_object();
    if (std::string(direct_local.at("chart").as_string()) != first_equation)
      throw std::runtime_error(
          "ordinary physical hop did not retain its equation owner");
    const auto direct_stats = request(json::object{
        {"schema", 2}, {"op", "local.stats"}, {"session", session},
        {"local", direct_local.at("local")}});
    require_ok(direct_stats, "physical wiring direct local.stats");
    if (direct_stats.at("epsilon_min") != -1 ||
        direct_stats.at("epsilon_max") != 2 ||
        direct_stats.at("top_valid") != 2 ||
        direct_stats.at("tail_majorant").as_object().at("status") !=
            "unsupported")
      throw std::runtime_error(
          "ordinary physical hop did not zero-pad its structural lower row or claimed a tail theorem: " +
          json::serialize(direct_stats));
    const auto canonical_zero_row_trimmed = [&](const char* point) {
      const auto evaluated = request(json::object{
          {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
          {"local", direct_local.at("local")},
          {"point", json::object{{"exact", point}}},
          {"options", json::object{{"tail_estimate", false}}},
          {"output_digits", 40}});
      require_ok(evaluated, "physical wiring padded local.evaluate");
      const auto& value = evaluated.at("value").as_object();
      const auto& coefficients = value.at("coefficients").as_array();
      // Point evaluation intentionally canonicalizes away exact leading
      // epsilon-zero rows.  The retained frame assertion above proves the
      // structural row exists; this verifies it cannot reappear numerically.
      return value.at("min") == 0 && value.at("max") == 2 &&
          coefficients.size() == 3;
    };
    if (!canonical_zero_row_trimmed("0") ||
        !canonical_zero_row_trimmed("1/3"))
      throw std::runtime_error(
          "ordinary physical hop did not canonically trim its exact structural lower row during local evaluation");

    const auto unavailable_tail = consume_physical_value_hop(
        session, plan, "lower", 1,
        std::string(direct_local.at("local").as_string()),
        std::string(direct_local.at("checkpoint_identity").as_string()),
        physical_value_solver("-4/3", 12), "physical-wiring-hop",
        plan_checkpoint);
    require_ok(unavailable_tail, "physical wiring unavailable-tail hop");
    if (unavailable_tail.at("used") != false ||
        unavailable_tail.at("execution_mode") !=
            "framed-fallback-required" ||
        unavailable_tail.at("reason") !=
            "incoming-local-has-no-certified-regular-tail-model")
      throw std::runtime_error(
          "unsupported physical tail did not fail side-effect-free into the framed path: " +
          json::serialize(unavailable_tail));

    const auto fallback_basis = solve_local(
        session, second_framed_chart, "-4/3",
        "physical-wiring-framed-fallback", "1", 12, true, -1, 2,
        second_equation);
    const auto fallback = request(json::object{
        {"schema", 2}, {"op", "transport.consume_hop"},
        {"session", session}, {"tile_plan", plan},
        {"tile_plan_checkpoint_identity", plan_checkpoint},
        {"arm", "lower"}, {"match", 1},
        {"receiving_basis", json::array{fallback_basis}},
        {"incoming", direct_local.at("local")},
        {"incoming_checkpoint_identity",
         direct_local.at("checkpoint_identity")},
        {"epsilon", json::object{{"min", 0}, {"max", 2},
                                    {"required_complete_max", 2}}},
        {"refinement", json::object{{"relative_tolerance", "1e-6"},
                                      {"max_steps", 2}}},
        {"checkpoint_policy", json::object{
             {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
             {"root", "physical-wiring-framed"}}}});
    require_ok(fallback, "physical wiring framed fallback");
    if (std::string(fallback.at("next_local").as_object()
                        .at("chart").as_string()) != second_equation ||
        fallback.at("consumed_basis_handles").as_array().size() != 1)
      throw std::runtime_error(
          "framed fallback did not materialize under the exact equation owner: " +
          json::serialize(fallback));

    const auto counters = request(json::object{
        {"schema", 2}, {"op", "session.counters"},
        {"session", session}});
    require_ok(counters, "physical wiring session.counters");
    if (counter(counters, "transport_physical_value_hop_attempts") != 3 ||
        counter(counters, "transport_physical_value_hop_successes") != 1 ||
        counter(counters, "transport_physical_value_hop_ineligible") != 1 ||
        counter(counters, "transport_framed_basis_hops") != 1)
      throw std::runtime_error(
          "physical/framed wiring diagnostics do not expose the alternating path: " +
          json::serialize(counters));

    for (const auto& owner : {first_equation, second_equation})
      require_ok(request(json::object{
          {"schema", 2}, {"op", "regular_equation.release"},
          {"session", session}, {"equation_owner", owner}}),
          "physical wiring regular_equation.release");
    const auto saved_checkpoint = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", session}, {"path", checkpoint},
        {"checkpoint_identity", "physical-wiring-hidden-owner"}});
    require_ok(saved_checkpoint,
               "physical wiring hidden-owner checkpoint.save");
    if (saved_checkpoint.at("regular_equation_owners") != 0 ||
        !std::filesystem::exists(checkpoint))
      throw std::runtime_error(
          "checkpoint did not distinguish hidden regular equation owners from public registry visibility: " +
          json::serialize(saved_checkpoint));
    const auto restored_record = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint},
        {"expected_identity", "physical-wiring-hidden-owner"}});
    require_ok(restored_record,
               "physical wiring hidden-owner checkpoint.restore");
    restored_session =
        std::string(restored_record.at("session").as_string());
    if (!restored_record.at("regular_equation_owners").as_array().empty())
      throw std::runtime_error(
          "hidden regular equation owner became publicly visible after restore");
    const auto resaved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", restored_session}, {"path", resaved_checkpoint},
        {"checkpoint_identity", "physical-wiring-hidden-owner-resaved"}});
    require_ok(resaved,
               "physical wiring restored hidden-owner checkpoint.save");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"},
        {"session", restored_session}}),
        "physical wiring restored session.close");
    restored_session.clear();
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"}, {"session", session}}),
        "physical wiring session.close");
    session.clear();
  } catch (...) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    if (!restored_session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored_session}});
    std::remove(checkpoint.c_str());
    std::remove(resaved_checkpoint.c_str());
    throw;
  }
  std::remove(checkpoint.c_str());
  std::remove(resaved_checkpoint.c_str());
}

void test_acb_terminal_factorized_consumed_checkpoint() {
  const std::string checkpoint =
      "/tmp/diffexp2-terminal-factorized-consumed.de2cp";
  const std::string corrupted =
      "/tmp/diffexp2-terminal-factorized-consumed-corrupt.de2cp";
  const std::string prepublication_checkpoint =
      "/tmp/diffexp2-terminal-factorized-prepublication.de2cp";
  std::remove(checkpoint.c_str());
  std::remove(corrupted.c_str());
  std::remove(prepublication_checkpoint.c_str());
  std::string session;
  std::string restored;
  std::string restored_single;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "acb"}, {"precision_bits", 256},
        {"output_digits", 40}, {"chart_capacity", 6},
        {"local_capacity", 24}, {"match_capacity", 8},
        {"tile_plan_capacity", 2},
        {"transport_state_capacity", 4},
        {"line_result_capacity", 4}, {"endpoint_capacity", 4}});
    require_ok(created, "terminal factorized session.create");
    session = std::string(created.at("session").as_string());

    const auto anchor_chart = prepare_multiblock_chart(
        session, "terminal-factorized-anchor", "0", "acb");
    const auto lower_chart = prepare_multiblock_chart(
        session, "terminal-factorized-lower", "-2/3", "acb");
    const auto upper_chart = prepare_multiblock_chart(
        session, "terminal-factorized-upper", "2/3", "acb");
    const std::string anchor_checkpoint =
        "terminal-factorized-anchor-local";
    const std::string plan_checkpoint =
        "terminal-factorized-plan";
    const auto anchor = solve_multiblock_local(
        session, anchor_chart, "0", anchor_checkpoint,
        {"1", "1"}, 8);
    const auto lower_basis_0 = solve_multiblock_local(
        session, lower_chart, "-2/3",
        "terminal-factorized-lower-basis-0", {"1", "1"}, 8);
    const auto lower_basis_1 = solve_multiblock_local(
        session, lower_chart, "-2/3",
        "terminal-factorized-lower-basis-1", {"0", "1"}, 8);
    const auto upper_basis_0 = solve_multiblock_local(
        session, upper_chart, "2/3",
        "terminal-factorized-upper-basis-0", {"1", "1"}, 8);
    const auto upper_basis_1 = solve_multiblock_local(
        session, upper_chart, "2/3",
        "terminal-factorized-upper-basis-1", {"0", "1"}, 8);
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", plan_checkpoint},
        {"division_order", 3},
        {"lower", arm("-2/3", anchor_chart, lower_chart)},
        {"upper", arm("2/3", anchor_chart, upper_chart)}});
    require_ok(planned, "terminal factorized tile.plan");
    const auto plan =
        std::string(planned.at("tile_plan").as_string());
    const auto lower = consume_basis_hop(
        session, plan, plan_checkpoint, "lower", anchor,
        anchor_checkpoint, {lower_basis_0, lower_basis_1},
        "terminal-factorized");
    const auto upper = consume_basis_hop(
        session, plan, plan_checkpoint, "upper", anchor,
        anchor_checkpoint, {upper_basis_0, upper_basis_1},
        "terminal-factorized");
    require_ok(lower, "terminal factorized lower consume_hop");
    require_ok(upper, "terminal factorized upper consume_hop");

    const auto saved_prepublication = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", session}, {"path", prepublication_checkpoint},
        {"checkpoint_identity",
         "terminal-factorized-prepublication"}});
    require_ok(saved_prepublication,
               "terminal factorized prepublication checkpoint.save");
    const auto restored_prepublication = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", prepublication_checkpoint},
        {"expected_identity",
         "terminal-factorized-prepublication"}});
    require_ok(restored_prepublication,
               "terminal factorized prepublication checkpoint.restore");
    restored_single = std::string(
        restored_prepublication.at("session").as_string());
    const auto single_published = request(json::object{
        {"schema", 2},
        {"op", "transport.publish_consumed_state"},
        {"session", restored_single},
        {"tile_plan", plan},
        {"tile_plan_checkpoint_identity", plan_checkpoint},
        {"anchor", anchor},
        {"anchor_checkpoint_identity", anchor_checkpoint},
        {"epsilon", epsilon_contract()},
        {"refinement", refinement()},
        {"checkpoint_policy", json::object{
             {"schema",
              "diffexp2-deterministic-arm-checkpoints-v1"},
             {"root", "terminal-factorized-single"}}},
        {"arm", "lower"},
        {"tile_sources", json::array{
             anchor,
             lower.at("next_local").as_object().at("local")}}});
    require_ok(single_published,
               "terminal factorized single-state publication");
    const auto& single_states =
        single_published.at("states").as_object();
    if (single_states.size() != 1 ||
        single_states.if_contains("lower") == nullptr ||
        single_states.if_contains("upper") != nullptr)
      throw std::runtime_error(
          "single consumed-state publication changed its requested arm");
    const auto single_value = contract_cancellation_and_export(
        restored_single, single_states.at("lower").as_object(),
        "terminal-factorized-single", 9);
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"},
        {"session", restored_single}}),
        "terminal factorized single-state session.close");
    restored_single.clear();

    const auto published = publish_consumed_states(
        session, plan, anchor, lower, upper, "terminal-factorized",
        plan_checkpoint, anchor_checkpoint);
    require_ok(published, "terminal factorized state publication");
    const auto& states = published.at("states").as_object();
    const auto lower_state = states.at("lower").as_object();
    const auto upper_state = states.at("upper").as_object();
    if (lower_state.at("terminal_factorized_match") != true ||
        upper_state.at("terminal_factorized_match") != true)
      throw std::runtime_error(
          "Acb terminal state publication dropped its factorized match owner: " +
          json::serialize(published));
    for (const auto* state : {&lower_state, &upper_state}) {
      const auto& contracts =
          state->at("tile_consumer_epsilon").as_array();
      if (contracts.size() != counter(*state, "tiles") ||
          contracts.empty())
        throw std::runtime_error(
            "terminal state did not publish one consumer epsilon contract per tile");
      const auto& terminal = contracts.back().as_object();
      if (terminal.at("sources").as_array().size() != 2 ||
          !terminal.at("output_min_shift").is_int64() ||
          !terminal.at("basis_min_power").is_int64() ||
          !terminal.at("projection_weight_min_power").is_int64())
        throw std::runtime_error(
            "terminal state did not publish its physical consumer sources, factorized-basis minimum, and transformed-weight minimum");
    }

    const auto require_small_zero = [](const json::value& raw,
                                       const char* label) {
      const auto& value = raw.as_object();
      if (value.at("dimension") != 1)
        throw std::runtime_error(
            std::string(label) + " did not remain scalar");
      for (const auto& coefficient :
           value.at("coefficients").as_array()) {
        const auto ball = decode_encoded_ball(coefficient);
        if (!ball.contains_zero() ||
            diffexp2::Magnitude::upper_abs(ball)
                    .approximate_upper() > 1e-45)
          throw std::runtime_error(
              std::string(label) +
              " did not preserve the terminal cancellation");
      }
    };
    require_small_zero(single_value,
                       "single-state terminal contraction");
    // Keep one bounded-cancellation observable public through the checkpoint
    // round trip.  This exercises the compact single-arm recipe schema and
    // binds its restored diagnostics back to the explicit cancellation
    // policy, rather than only restoring the exact-singleton pair form.
    const auto lower_value = contract_cancellation_and_export(
        session, lower_state, "terminal-factorized-lower", 9, false);
    const auto upper_value = contract_cancellation_and_export(
        session, upper_state, "terminal-factorized-upper", 9);
    require_small_zero(lower_value, "lower terminal contraction");
    require_small_zero(upper_value, "upper terminal contraction");
    if (setenv("DE2_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE",
               "factorized", 1) != 0)
      throw std::runtime_error(
          "could not select explicit factorized terminal contraction");
    json::value explicit_factorized_value;
    try {
      explicit_factorized_value = contract_cancellation_and_export(
          session, lower_state, "terminal-explicit-factorized", 9);
    } catch (...) {
      unsetenv("DE2_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE");
      throw;
    }
    unsetenv("DE2_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE");
    require_small_zero(
        explicit_factorized_value,
        "explicit factorized terminal contraction");
    if (setenv("DE2_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE",
               "compare", 1) != 0)
      throw std::runtime_error(
          "could not select terminal direct/adjoint comparison");
    if (setenv("DE2_DIAGNOSTIC_TERMINAL_COMPOSED_ADJOINT",
               "require", 1) != 0) {
      unsetenv("DE2_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE");
      throw std::runtime_error(
          "could not require the certified composed terminal adjoint");
    }
    json::value compared_value;
    try {
      compared_value = contract_cancellation_and_export(
          session, lower_state, "terminal-compare-factorized", 9);
    } catch (...) {
      unsetenv("DE2_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE");
      unsetenv("DE2_DIAGNOSTIC_TERMINAL_COMPOSED_ADJOINT");
      throw;
    }
    unsetenv("DE2_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE");
    unsetenv("DE2_DIAGNOSTIC_TERMINAL_COMPOSED_ADJOINT");
    require_small_zero(
        compared_value,
        "compared direct/adjoint terminal contraction");
    // The direct physical route contracts L(F) with the certified physical
    // weights w.  It must remain a genuine independent route: in particular,
    // it must not pass through the epsilon-shifted factorization (F T) P and
    // accidentally inherit that frame's much deeper convolution window.
    if (setenv("DE2_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE",
               "physical", 1) != 0)
      throw std::runtime_error(
          "could not enable direct physical terminal contraction");
    json::value physical_value;
    try {
      physical_value = contract_cancellation_and_export(
          session, lower_state, "terminal-direct-physical", 9);
    } catch (...) {
      unsetenv("DE2_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE");
      throw;
    }
    unsetenv("DE2_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE");
    require_small_zero(
        physical_value, "direct physical terminal contraction");
    if (setenv("DE2_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE",
               "invalid-route", 1) != 0)
      throw std::runtime_error(
          "could not select invalid terminal contraction route");
    bool invalid_route_rejected = false;
    try {
      (void)contract_cancellation_and_export(
          session, lower_state, "terminal-invalid-route", 9);
    } catch (...) {
      invalid_route_rejected = true;
    }
    unsetenv("DE2_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE");
    if (!invalid_route_rejected)
      throw std::runtime_error(
          "invalid terminal contraction route was not rejected");

    const auto endpoints = request(json::object{
        {"schema", 2}, {"op", "transport.endpoint_batch"},
        {"session", session},
        {"transport_state", lower_state.at("transport_state")},
        {"transport_state_checkpoint_identity",
         lower_state.at("checkpoint_identity")},
        {"transport_state_provenance_identity",
         lower_state.at("provenance_identity")},
        {"checkpoint_policy", json::object{
             {"schema",
              "diffexp2-deterministic-transport-endpoint-checkpoints-v1"},
             {"root", "terminal-factorized-endpoint"}}},
        {"observables", json::array{json::object{
             {"identity", "terminal-factorized-endpoint-observable"},
             {"checkpoint_identity",
              "terminal-factorized-endpoint-result"},
             {"integrand_row", cancellation_integrand_row(
                  "terminal-factorized-endpoint-row", 9)},
             {"epsilon", json::object{
                  {"min", 0}, {"max", 2},
                  {"required_complete_max", 1}}}}}}});
    require_ok(endpoints, "terminal factorized endpoint_batch");
    if (endpoints.at("endpoints").as_array().size() != 1 ||
        endpoints.at("no_projected_local_publication") != true)
      throw std::runtime_error(
          "terminal factorized endpoint batch published scratch state");
    const auto& endpoint =
        endpoints.at("endpoints").as_array().front().as_object();
    const auto endpoint_export = request(json::object{
        {"schema", 2}, {"op", "endpoint.export"}, {"session", session},
        {"endpoint", endpoint.at("endpoint")},
        {"checkpoint_identity", endpoint.at("checkpoint_identity")},
        {"output_digits", 40}});
    require_ok(endpoint_export, "terminal factorized endpoint.export");
    require_small_zero(endpoint_export.at("value"),
                       "terminal factorized endpoint");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "endpoint.release"},
        {"session", session}, {"endpoint", endpoint.at("endpoint")}}),
        "terminal factorized endpoint.release");

    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", session}, {"path", checkpoint},
        {"checkpoint_identity", "terminal-factorized-roundtrip"}});
    require_ok(saved, "terminal factorized checkpoint.save");
    const auto container = diffexp2::checkpoint::read(checkpoint);
    const auto payload =
        json::parse(container.payload_json).as_object();
    const auto& retained_states =
        payload.at("retained_transport_states").as_array();
    const auto& retained_hops =
        payload.at("retained_planned_match_hops").as_array();
    if (retained_states.size() != 2 || retained_hops.size() != 2)
      throw std::runtime_error(
          "terminal checkpoint did not retain exactly two state/match closures");
    for (const auto& raw_state : retained_states) {
      const auto& state = raw_state.as_object();
      if (state.at("schema") !=
              "diffexp2-retained-transport-arm-state-v6" ||
          state.at("provenance").as_object().if_contains(
              "terminal_match") == nullptr)
        throw std::runtime_error(
            "terminal checkpoint lost its v6 match reference");
    }

    auto corrupted_payload = payload;
    corrupted_payload.at("retained_transport_states")
        .as_array().front().as_object()
        .at("provenance").as_object()
        .at("terminal_match").as_object()
        .at("checkpoint_identity") =
            "tampered-terminal-match-checkpoint";
    diffexp2::checkpoint::write_atomic(
        corrupted, container.header_json,
        json::serialize(corrupted_payload));
    const auto rejected = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", corrupted},
        {"expected_identity", "terminal-factorized-roundtrip"}});
    if (rejected.at("status") != "error")
      throw std::runtime_error(
          "checkpoint restore accepted a stale terminal match reference");

    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"},
        {"session", session}}),
        "terminal factorized session.close");
    session.clear();
    const auto restored_record = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint},
        {"expected_identity", "terminal-factorized-roundtrip"}});
    require_ok(restored_record,
               "terminal factorized checkpoint.restore");
    restored =
        std::string(restored_record.at("session").as_string());
    if (restored_record.at("transport_states").as_array().size() != 2 ||
        !restored_record.at("planned_match_hops").as_array().empty())
      throw std::runtime_error(
          "terminal factorized restore lost its public state/match closure: " +
          json::serialize(restored_record));
    const auto restored_stats = request(json::object{
        {"schema", 2}, {"op", "transport.stats"},
        {"session", restored},
        {"transport_state", lower_state.at("transport_state")}});
    require_ok(restored_stats,
               "restored terminal factorized transport.stats");
    if (restored_stats.at("terminal_factorized_match") != true ||
        restored_stats.at("strong_ownership").as_object().at(
            "terminal_factorized_match") != true)
      throw std::runtime_error(
          "restored terminal state did not privately retain its match owner");
    const auto restored_lower = contract_cancellation_and_export(
        restored, lower_state, "terminal-factorized-restored", 9);
    require_small_zero(restored_lower,
                       "restored terminal contraction");
    if (!epsilon_vectors_overlap(lower_value, restored_lower))
      throw std::runtime_error(
          "terminal factorized checkpoint roundtrip changed its result");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"},
        {"session", restored}}),
        "restored terminal factorized session.close");
    restored.clear();
  } catch (...) {
    if (!session.empty())
      (void)request(json::object{
          {"schema", 2}, {"op", "session.close"},
          {"session", session}});
    if (!restored.empty())
      (void)request(json::object{
          {"schema", 2}, {"op", "session.close"},
          {"session", restored}});
    if (!restored_single.empty())
      (void)request(json::object{
          {"schema", 2}, {"op", "session.close"},
          {"session", restored_single}});
    std::remove(checkpoint.c_str());
    std::remove(corrupted.c_str());
    std::remove(prepublication_checkpoint.c_str());
    throw;
  }
  std::remove(checkpoint.c_str());
  std::remove(corrupted.c_str());
  std::remove(prepublication_checkpoint.c_str());
}

}  // namespace

int main() {
  const std::string first_checkpoint =
      "/tmp/diffexp2-transport-run-arms-first.de2cp";
  const std::string second_checkpoint =
      "/tmp/diffexp2-transport-run-arms-second.de2cp";
  std::remove(first_checkpoint.c_str());
  std::remove(second_checkpoint.c_str());
  std::string session;
  std::string restored_first;
  std::string restored_second;
  try {
    test_regular_value_hop_checkpoint();
    test_acb_prescribed_value_hops_match_basis();
    test_acb_value_handoff_significance_gate();
    test_acb_consuming_hop_reservoir_retry_lifecycle();
    test_multiblock_regular_value_fallback_owner();
    test_streaming_consumed_transport();
    test_consuming_transport();
    test_frame_independent_physical_value_wiring();
    test_acb_terminal_factorized_consumed_checkpoint();
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "rational"}, {"output_digits", 40},
        {"chart_capacity", 6}, {"local_capacity", 24},
        {"match_capacity", 12}, {"tile_plan_capacity", 2},
        {"transport_state_capacity", 8}, {"line_result_capacity", 8}});
    require_ok(created, "session.create");
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_chart(
        session, "paired-state-anchor-chart", "0");
    const auto lower_chart = prepare_chart(
        session, "paired-state-lower-chart", "-2/3");
    const auto upper_chart = prepare_chart(
        session, "paired-state-upper-chart", "2/3");
    const auto anchor = solve_local(
        session, anchor_chart, "0", "paired-state-anchor", "2");
    const auto lower_basis = solve_local(
        session, lower_chart, "-2/3", "paired-state-lower-basis", "1");
    const auto upper_basis = solve_local(
        session, upper_chart, "2/3", "paired-state-upper-basis", "1");
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "paired-state-plan"},
        {"division_order", 3},
        {"lower", arm("-2/3", anchor_chart, lower_chart)},
        {"upper", arm("2/3", anchor_chart, upper_chart)}});
    require_ok(planned, "tile.plan");
    const auto plan = std::string(planned.at("tile_plan").as_string());
    if (planned.at("lower_matches") != 1 ||
        planned.at("upper_matches") != 1)
      throw std::runtime_error(
          "fixture did not produce one match per arm: " +
          json::serialize(planned));

    const auto before_failure = session_stats(session);
    const auto failed = run_pair(
        session, plan, anchor, lower_basis, lower_basis,
        "paired-state-failure");
    const auto after_failure = session_stats(session);
    if (failed.at("status") != "error" ||
        before_failure.at("locals") != after_failure.at("locals") ||
        before_failure.at("matches") != after_failure.at("matches") ||
        before_failure.at("transport_states") !=
            after_failure.at("transport_states") ||
        before_failure.at("local_matches") !=
            after_failure.at("local_matches") ||
        before_failure.at("transport_arm_marches") !=
            after_failure.at("transport_arm_marches") ||
        after_failure.at("pending_local_solves") != 0 ||
        after_failure.at("pending_matches") != 0 ||
        after_failure.at("pending_transport_states") != 0)
      throw std::runtime_error(
          "one-arm failure was not a complete rollback: " +
          json::serialize(failed) + " / " + json::serialize(after_failure));

    const auto before_pair = session_stats(session);
    const auto paired = run_pair(
        session, plan, anchor, lower_basis, upper_basis,
        "paired-state-success");
    require_ok(paired, "transport.run_arms");
    const auto after_pair = session_stats(session);
    const auto& states = paired.at("states").as_object();
    const auto lower_state = states.at("lower").as_object();
    const auto upper_state = states.at("upper").as_object();
    if (paired.at("workers") != 2 ||
        paired.at("max_parallel_arms") != 2 ||
        paired.at("worker_overlap") != true ||
        paired.at("atomic_publication") != true ||
        paired.at("public_result_tokens") != "transport_states_only" ||
        paired.at("dependency_only_final_locals") != true ||
        paired.at("marches") != 2 ||
        paired.at("matches").as_object().at("total") != 2 ||
        lower_state.at("matches") != 1 || lower_state.at("tiles") != 2 ||
        upper_state.at("matches") != 1 || upper_state.at("tiles") != 2 ||
        lower_state.at("final_local").as_object().at("public_token") != false ||
        upper_state.at("final_local").as_object().at("public_token") != false ||
        before_pair.at("locals") != after_pair.at("locals") ||
        before_pair.at("matches") != after_pair.at("matches") ||
        before_pair.at("line_results") != after_pair.at("line_results") ||
        counter(after_pair, "transport_states") !=
            counter(before_pair, "transport_states") + 2 ||
        counter(after_pair, "local_matches") !=
            counter(before_pair, "local_matches") + 2 ||
        counter(after_pair, "transport_arm_marches") !=
            counter(before_pair, "transport_arm_marches") + 2 ||
        paired.at("plan_stats").as_object().at("lower_match_advances") != 1 ||
        paired.at("plan_stats").as_object().at("upper_match_advances") != 1)
      throw std::runtime_error(
          "two-arm result violated publication/overlap accounting: " +
          json::serialize(paired) + " / " + json::serialize(after_pair));

    const auto hidden_lower = std::string(
        lower_state.at("final_local").as_object().at("local").as_string());
    const auto hidden_release = request(json::object{
        {"schema", 2}, {"op", "local.release"}, {"session", session},
        {"local", hidden_lower}});
    if (hidden_release.at("status") != "error")
      throw std::runtime_error(
          "dependency-only final local became a public token: " +
          json::serialize(hidden_release));

    const auto single_lower = run_single(
        session, plan, anchor, "lower", lower_basis,
        "paired-state-single-lower");
    const auto single_upper = run_single(
        session, plan, anchor, "upper", upper_basis,
        "paired-state-single-upper");
    require_ok(single_lower, "single lower transport");
    require_ok(single_upper, "single upper transport");
    if (single_lower.at("matches") != lower_state.at("matches") ||
        single_lower.at("tiles") != lower_state.at("tiles") ||
        single_lower.at("epsilon") != lower_state.at("epsilon") ||
        single_upper.at("matches") != upper_state.at("matches") ||
        single_upper.at("tiles") != upper_state.at("tiles") ||
        single_upper.at("epsilon") != upper_state.at("epsilon"))
      throw std::runtime_error(
          "paired states differ topologically from single-arm marches");
    const auto paired_lower_value = contract_and_export(
        session, lower_state, "paired-value-lower");
    const auto single_lower_value = contract_and_export(
        session, single_lower, "single-value-lower");
    const auto paired_upper_value = contract_and_export(
        session, upper_state, "paired-value-upper");
    const auto single_upper_value = contract_and_export(
        session, single_upper, "single-value-upper");
    if (paired_lower_value != single_lower_value ||
        paired_upper_value != single_upper_value)
      throw std::runtime_error(
          "paired transport values differ from single-arm values");

    release_state(session, std::string(
        single_lower.at("transport_state").as_string()));
    release_state(session, std::string(
        single_upper.at("transport_state").as_string()));
    release_local(session, std::string(single_lower.at("final_local")
        .as_object().at("local").as_string()));
    release_local(session, std::string(single_upper.at("final_local")
        .as_object().at("local").as_string()));
    require_ok(request(json::object{
        {"schema", 2}, {"op", "tile.release"}, {"session", session},
        {"tile_plan", plan}}), "tile.release");
    release_local(session, anchor);
    release_local(session, lower_basis);
    release_local(session, upper_basis);
    const auto hidden_only = session_stats(session);
    if (hidden_only.at("locals") != 0 || hidden_only.at("matches") != 0 ||
        hidden_only.at("tile_plans") != 0 ||
        hidden_only.at("transport_states") != 2 ||
        hidden_only.at("line_results") != 0)
      throw std::runtime_error(
          "public source release did not leave exactly two state tokens: " +
          json::serialize(hidden_only));

    const auto saved_first = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
        {"path", first_checkpoint},
        {"checkpoint_identity", "paired-state-roundtrip-1"}});
    require_ok(saved_first, "first checkpoint.save");
    release_state(session, std::string(
        lower_state.at("transport_state").as_string()));
    release_state(session, std::string(
        upper_state.at("transport_state").as_string()));
    const auto released_stats = request(json::object{
        {"schema", 2}, {"op", "transport.stats"}, {"session", session},
        {"transport_state", lower_state.at("transport_state")}});
    if (released_stats.at("status") != "error")
      throw std::runtime_error("transport release did not hide the state");
    require_ok(request(json::object{{"schema", 2}, {"op", "session.close"},
                                    {"session", session}}),
               "session.close");
    session.clear();

    const auto first_restore = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", first_checkpoint},
        {"expected_identity", "paired-state-roundtrip-1"}});
    require_ok(first_restore, "first checkpoint.restore");
    restored_first = std::string(first_restore.at("session").as_string());
    if (first_restore.at("locals").as_array().size() != 0 ||
        first_restore.at("tile_plans").as_array().size() != 0 ||
        first_restore.at("planned_match_hops").as_array().size() != 0 ||
        first_restore.at("transport_states").as_array().size() != 2)
      throw std::runtime_error(
          "first restore exposed hidden closure objects: " +
          json::serialize(first_restore));
    const auto restored_lower_stats = request(json::object{
        {"schema", 2}, {"op", "transport.stats"},
        {"session", restored_first},
        {"transport_state", lower_state.at("transport_state")}});
    require_ok(restored_lower_stats, "restored transport.stats");
    if (contract_and_export(restored_first, lower_state,
                            "restored-value-lower") != paired_lower_value)
      throw std::runtime_error(
          "first restored hidden closure changed the lower value");
    const auto saved_second = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", restored_first}, {"path", second_checkpoint},
        {"checkpoint_identity", "paired-state-roundtrip-2"}});
    require_ok(saved_second, "second checkpoint.save");
    require_ok(request(json::object{{"schema", 2}, {"op", "session.close"},
                                    {"session", restored_first}}),
               "first restored session.close");
    restored_first.clear();

    const auto second_restore = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", second_checkpoint},
        {"expected_identity", "paired-state-roundtrip-2"}});
    require_ok(second_restore, "second checkpoint.restore");
    restored_second = std::string(second_restore.at("session").as_string());
    if (second_restore.at("locals").as_array().size() != 0 ||
        second_restore.at("tile_plans").as_array().size() != 0 ||
        second_restore.at("planned_match_hops").as_array().size() != 0 ||
        second_restore.at("transport_states").as_array().size() != 2 ||
        contract_and_export(restored_second, upper_state,
                            "restored-value-upper") != paired_upper_value)
      throw std::runtime_error(
          "second restore did not preserve the two hidden state closures");
    release_state(restored_second, std::string(
        lower_state.at("transport_state").as_string()));
    release_state(restored_second, std::string(
        upper_state.at("transport_state").as_string()));
    const auto second_release = request(json::object{
        {"schema", 2}, {"op", "transport.release"},
        {"session", restored_second},
        {"transport_state", lower_state.at("transport_state")}});
    if (second_release.at("status") != "error")
      throw std::runtime_error("second state release was not loud");
    require_ok(request(json::object{{"schema", 2}, {"op", "session.close"},
                                    {"session", restored_second}}),
               "second restored session.close");
    restored_second.clear();
    std::remove(first_checkpoint.c_str());
    std::remove(second_checkpoint.c_str());
    std::cout << "persistent atomic transport.run_arms passed\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    if (!restored_first.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored_first}});
    if (!restored_second.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored_second}});
    std::remove(first_checkpoint.c_str());
    std::remove(second_checkpoint.c_str());
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
