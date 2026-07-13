#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <array>
#include <cstdint>
#include <cstdio>
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

json::array prescriptions(bool branch_sensitive = true) {
  if (!branch_sensitive) return json::array{};
  return json::array{json::object{
      {"factor_exact", "paired-state-f"}, {"sign", -1},
      {"multiplicity", 1}, {"leading_coefficient_sign", 1}}};
}

json::object epsilon_rational_one() {
  return json::object{{"zero", false}, {"valuation", 0},
                      {"numerator", json::array{"1"}},
                      {"denominator", json::array{"1"}}};
}

std::string prepare_chart(const std::string& session,
                          const std::string& name,
                          const std::string& center,
                          bool branch_sensitive = true,
                          const std::string& domain = "rational",
                          int chop_digits = 0) {
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
  json::array blocks;
  blocks.push_back(std::move(component));
  const auto owner = "de2-operator-" + name;
  const auto one = epsilon_rational_one();
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
                {"q", json::array{one}},
                {"c", std::move(physical_c)}}},
           {"chop_digits", chop_digits}};
  if (domain == "acb") problem["precision_bits"] = 256;
  const auto prepared = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", name}, {"identity", owner},
      {"analytic", json::object{
           {"geometry", json::object{
                {"center_exact", center}, {"scale_exact", "1"},
                {"radius_exact", "2"}, {"infinite_radius", false},
                {"prescriptions", prescriptions(branch_sensitive)}}},
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
                        bool branch_sensitive = true) {
  json::array schedule;
  json::array shifts;
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
           {"nmax", nmax}, {"p", 0}, {"has_initial", true},
           {"adaptive_probe", false}, {"a_target", "0"},
           {"b_target", "0"}, {"a_shift_min", 0},
           {"a_shifts", std::move(shifts)},
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
           {"prescriptions", prescriptions(branch_sensitive)},
           {"checkpoint_identity", checkpoint}}}});
  require_ok(solved, "local.solve");
  return std::string(solved.at("local").as_string());
}

std::string prepare_multiblock_chart(const std::string& session,
                                     const std::string& name,
                                     const std::string& center) {
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
  const auto prepared = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", name}, {"identity", owner},
      {"analytic", json::object{
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
      {"problem", json::object{
           {"domain", "rational"}, {"d", dimension}, {"fb", 0},
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
           {"chop_digits", 0}}}});
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

json::object topology(bool branch_sensitive = true) {
  return json::object{
      {"singular_points", json::array{}},
      {"boundary_points", json::array{}},
      {"complex_projections", json::array{}},
      {"branch_sheets", branch_sensitive
           ? json::array{json::object{
                 {"factor_exact", "paired-state-f"}, {"sign", -1}}}
           : json::array{}}};
}

json::object arm(const std::string& endpoint,
                 const std::string& anchor_chart,
                 const std::string& receiving_chart,
                 bool branch_sensitive = true) {
  return json::object{
      {"from_exact", "0"}, {"to_exact", endpoint},
      {"charts", json::array{anchor_chart, receiving_chart}},
      {"topology", topology(branch_sensitive)}};
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
                          const std::string& tail_proxy_max = "1/100") {
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
           {"prescriptions", prescriptions(branch_sensitive)},
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
    const json::object& upper, const std::string& root) {
  return request(json::object{
      {"schema", 2}, {"op", "transport.publish_consumed_states"},
      {"session", session}, {"tile_plan", plan},
      {"tile_plan_checkpoint_identity", "streaming-state-plan"},
      {"anchor", anchor},
      {"anchor_checkpoint_identity", "streaming-state-anchor"},
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

json::object integrand_row(const std::string& identity) {
  const json::array one{"1", "0", "0", "0", "0"};
  const json::array zero{"0", "0", "0", "0", "0"};
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

json::value contract_and_export(const std::string& session,
                                const json::object& state,
                                const std::string& root) {
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
                integrand_row(root + ":row:1"),
                integrand_row(root + ":row:2")}},
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
    require_ok(created, "value-hop session.create");
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_chart(
        session, "value-hop-anchor-chart", "0", false);
    const auto lower_chart = prepare_chart(
        session, "value-hop-lower-chart", "-2/3", false);
    const auto upper_chart = prepare_chart(
        session, "value-hop-upper-chart", "2/3", false);
    const auto anchor = solve_local(
        session, anchor_chart, "0", "streaming-state-anchor", "2", 4,
        false);
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "streaming-state-plan"},
        {"division_order", 3},
        {"lower", arm("-2/3", anchor_chart, lower_chart, false)},
        {"upper", arm("2/3", anchor_chart, upper_chart, false)}});
    require_ok(planned, "value-hop tile.plan");
    const auto plan = std::string(planned.at("tile_plan").as_string());
    const auto before = session_stats(session);
    const auto unsafe_tail = consume_value_hop(
        session, plan, "lower", anchor, "streaming-state-anchor",
        value_solver("-2/3", false, 0), "value-hop-success");
    require_ok(unsafe_tail, "unsafe-tail value hop");
    const auto after_unsafe_tail = session_stats(session);
    if (unsafe_tail.at("used") != false ||
        unsafe_tail.at("reason") !=
            "receiver-center-fails-exact-truncation-tail-contract" ||
        counter(after_unsafe_tail, "local_solves") !=
            counter(before, "local_solves"))
      throw std::runtime_error(
          "unsafe exact center tail did not fail closed before solving: " +
          json::serialize(unsafe_tail));

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
      auto solver = value_solver("-2/3");
      solver.at("run").as_object()["adaptive_probe"] = true;
      expect_rejected_template(std::move(solver), "adaptive probe",
                               "not one homogeneous (0,0,0) value run");
    }
    {
      auto solver = value_solver("-2/3");
      solver.at("run").as_object()["a_shift_min"] = -1;
      expect_rejected_template(std::move(solver), "nonzero shift origin",
                               "not one homogeneous (0,0,0) value run");
    }
    {
      auto solver = value_solver("-2/3");
      solver.at("run").as_object().at("a_shifts").as_array().emplace_back(
          "5");
      expect_rejected_template(std::move(solver), "extra a-shift",
                               "not one homogeneous (0,0,0) value run");
    }
    {
      auto solver = value_solver("-2/3");
      solver.at("run").as_object().at("a_shifts").as_array()[2] = "7";
      expect_rejected_template(std::move(solver), "non-Taylor a-shift",
                               "a-shifts are not the exact Taylor indices");
    }
    {
      auto solver = value_solver("-2/3");
      solver.at("run").as_object().at("schedule").as_array()[2]
          .as_array()[0].as_object()["da"] = "7";
      expect_rejected_template(std::move(solver), "non-Taylor da",
                               "Taylor by exact index");
    }
    {
      auto solver = value_solver("-2/3");
      solver.at("run").as_object().at("schedule").as_array()[2]
          .as_array()[0].as_object()["db"] = "1";
      expect_rejected_template(std::move(solver), "nonzero db",
                               "Taylor by exact index");
    }

    const auto lower = consume_value_hop(
        session, plan, "lower", anchor, "streaming-state-anchor",
        value_solver("-2/3"), "value-hop-success");
    const auto upper = consume_value_hop(
        session, plan, "upper", anchor, "streaming-state-anchor",
        value_solver("2/3"), "value-hop-success");
    require_ok(lower, "lower regular value hop");
    require_ok(upper, "upper regular value hop");
    const auto after_hops = session_stats(session);
    if (!lower.at("used").as_bool() || !upper.at("used").as_bool() ||
        lower.at("value_hops") != 1 || lower.at("basis_matches") != 0 ||
        upper.at("value_hops") != 1 || upper.at("basis_matches") != 0 ||
        counter(after_hops, "local_matches") !=
            counter(before, "local_matches") ||
        counter(after_hops, "local_solves") !=
            counter(before, "local_solves") + 2)
      throw std::runtime_error(
          "eligible regular hops did not use exactly one value solve each: " +
          json::serialize(lower) + " / " + json::serialize(upper));

    // Compare the optimized local directly with the established one-column
    // basis/match path before either becomes dependency-only state.
    const auto comparison_anchor = solve_local(
        session, anchor_chart, "0", "value-basis-comparison-anchor", "2",
        4, false);
    const auto comparison_basis = solve_local(
        session, lower_chart, "-2/3", "value-basis-comparison-column", "1",
        4, false);
    const auto comparison_planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "value-basis-comparison-plan"},
        {"division_order", 3},
        {"lower", arm("-2/3", anchor_chart, lower_chart, false)},
        {"upper", arm("2/3", anchor_chart, upper_chart, false)}});
    require_ok(comparison_planned, "value/basis comparison tile.plan");
    const auto comparison_plan = std::string(
        comparison_planned.at("tile_plan").as_string());
    const auto basis_hop = consume_hop(
        session, comparison_plan, "lower", comparison_anchor,
        "value-basis-comparison-anchor", comparison_basis,
        "value-basis-comparison", "value-basis-comparison-plan");
    require_ok(basis_hop, "value/basis comparison hop");
    const auto evaluate_center = [&](const json::value& local) {
      const auto evaluated = request(json::object{
          {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
          {"local", local},
          {"point", json::object{{"exact", "0"}}},
          {"options", json::object{{"tail_estimate", false}}},
          {"output_digits", 40}});
      require_ok(evaluated, "value/basis center evaluation");
      return evaluated.at("value");
    };
    if (evaluate_center(lower.at("next_local").as_object().at("local")) !=
        evaluate_center(basis_hop.at("next_local").as_object().at("local")))
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
        session, lower_state, "value-hop-before-checkpoint");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "tile.release"}, {"session", session},
        {"tile_plan", plan}}), "value-hop tile.release");
    release_local(session, anchor);
    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", session}, {"path", checkpoint},
        {"checkpoint_identity", "value-hop-roundtrip"}});
    require_ok(saved, "value-hop checkpoint.save");
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
        contract_and_export(restored, lower_state,
                            "value-hop-after-checkpoint") != lower_value)
      throw std::runtime_error(
          "regular value-hop checkpoint roundtrip changed its sealed state");
    require_ok(request(json::object{
        {"schema", 2}, {"op", "session.close"}, {"session", restored}}),
        "restored value-hop session.close");
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
            "center-evaluation-fails-relative-accuracy-contract" ||
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
        session, anchor_chart, "0", "multiblock-anchor-local", {"2", "3"});
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
        value_solver("-2/3", true), "multiblock-fallback",
        "multiblock-plan");
    require_ok(ineligible, "multiblock value preflight");
    const auto after_preflight = session_stats(session);
    if (ineligible.at("used") != false ||
        ineligible.at("reason") != "branch-sensitive-crossing" ||
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
    if (std::string(next.at("chart").as_string()) != lower_chart ||
        next.at("source_operator_identity") !=
            "de2-operator-multiblock-lower")
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
        value_solver("-2/3", true), "streaming-state-success");
    require_ok(ineligible_value, "branch-sensitive value-hop preflight");
    const auto after_ineligible_value = session_stats(session);
    if (ineligible_value.at("used") != false ||
        ineligible_value.at("reason") != "branch-sensitive-crossing" ||
        after_ineligible_value.at("locals") != lower_before.at("locals") ||
        counter(after_ineligible_value, "local_solves") !=
            counter(lower_before, "local_solves"))
      throw std::runtime_error(
          "branch-sensitive value hop did not fail closed before solving: " +
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
    test_acb_value_handoff_significance_gate();
    test_multiblock_regular_value_fallback_owner();
    test_streaming_consumed_transport();
    test_consuming_transport();
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
