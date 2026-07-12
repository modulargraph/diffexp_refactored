#include "diffexp2/checkpoint.hpp"
#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace json = boost::json;

namespace {

constexpr const char* kSingularProofSchema =
    "diffexp2-native-acb-singular-scc-valuation-zero-saturation-proof-v1";
constexpr const char* kOrdinaryProofSchema =
    "diffexp2-native-acb-unit-leading-saturation-proof-v1";

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
      .as_object();
}

std::uint64_t unsigned_value(const json::object& object, const char* key) {
  const auto& value = object.at(key);
  if (value.is_uint64()) return value.as_uint64();
  if (value.is_int64() && value.as_int64() >= 0)
    return static_cast<std::uint64_t>(value.as_int64());
  throw std::runtime_error(std::string("expected nonnegative counter: ") + key);
}

json::object geometry(const std::string& center) {
  return json::object{{"center_exact", center}, {"scale_exact", "1"},
                      {"radius_exact", "2"},
                      {"infinite_radius", false},
                      {"prescriptions", json::array{}}};
}

json::array singleton_components() {
  json::array component;
  component.push_back(0);
  json::array components;
  components.push_back(std::move(component));
  return components;
}

json::array unit_d_lags() {
  json::array lag;
  lag.push_back(json::object{{"s", 0}, {"v", "1"}});
  json::array lags;
  lags.push_back(std::move(lag));
  return lags;
}

json::array self_edges() {
  json::array edge;
  edge.push_back(0);
  edge.push_back(0);
  json::array edges;
  edges.push_back(std::move(edge));
  return edges;
}

json::array scalar_entry(const std::string& value) {
  json::array entry;
  entry.push_back(0);
  entry.push_back(0);
  entry.emplace_back(value);
  json::array entries;
  entries.push_back(std::move(entry));
  return entries;
}

std::string prepare_regular_chart(const std::string& session,
                                  const std::string& key,
                                  const std::string& center) {
  json::array principal_row;
  principal_row.push_back(
      json::object{{"exact", "0"}, {"proven_zero", true}});
  json::array principal_matrix;
  principal_matrix.push_back(std::move(principal_row));
  const auto prepared = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", key}, {"identity", key + "-identity"},
      {"analytic", json::object{
           {"geometry", geometry(center)},
           {"principal_matrix", std::move(principal_matrix)},
           {"native_scc_capabilities", json::object{
                {"regular", true}, {"identity_gauge", true},
                {"identity_v", true}, {"no_pseudo", true}}}}},
      {"scc", json::object{
           {"components", singleton_components()},
           {"structural_edges", json::array{}},
           {"condensation_edges", json::array{}},
           {"topological_order", json::array{0}},
           {"coupling_depth", 0}}},
      {"problem", json::object{
           {"domain", "acb"}, {"precision_bits", 256},
           {"d", 1}, {"fb", -1}, {"w", 4},
           {"d_lags", unit_d_lags()},
           {"denominators", json::array{}},
           {"nhat_lags", json::array{json::object{
                {"poly", json::array{}}, {"rat", json::array{}},
                {"val", json::array{nullptr}}}}},
           {"d0_inverse", "1"},
           {"blocks", singleton_components()},
           {"assembly", json::object{
                {"identity", true}, {"poly", json::array{}},
                {"rat", json::array{}}, {"val", json::array{0}}}},
           {"chop_digits", 0}}}});
  if (prepared.at("status") != "ok")
    throw std::runtime_error("regular chart.prepare: " +
                             json::serialize(prepared));
  return std::string(prepared.at("chart").as_string());
}

std::string prepare_singular_chart(const std::string& session) {
  json::array principal_row;
  principal_row.push_back(json::object{
      {"exact", "((1/2)+(eps/3))/t"}, {"proven_zero", false}});
  json::array principal_matrix;
  principal_matrix.push_back(std::move(principal_row));
  const auto prepared = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", "singular-arm-block"},
      {"identity", "singular-arm-block-identity"},
      {"analytic", json::object{
           {"geometry", geometry("-2/3")},
           {"principal_matrix", std::move(principal_matrix)},
           {"native_scc_capabilities", json::object{
                {"regular", false}, {"identity_gauge", true},
                {"identity_v", true}, {"no_pseudo", true}}}}},
      {"scc", json::object{
           {"components", singleton_components()},
           {"structural_edges", self_edges()},
           {"condensation_edges", json::array{}},
           {"topological_order", json::array{0}},
           {"coupling_depth", 0}}},
      {"problem", json::object{
           {"domain", "acb"}, {"precision_bits", 256},
           {"d", 1}, {"fb", -1}, {"w", 4},
           {"d_lags", unit_d_lags()},
           {"denominators", json::array{}},
           {"nhat_lags", json::array{json::object{
                {"poly", json::array{
                     json::object{{"s", 0},
                                  {"e", scalar_entry("1/2")}},
                     json::object{{"s", 1},
                                  {"e", scalar_entry("1/3")}}}},
                {"rat", json::array{}}, {"val", json::array{0}}}}},
           {"d0_inverse", "1"},
           {"blocks", singleton_components()},
           {"assembly", json::object{
                {"identity", true}, {"poly", json::array{}},
                {"rat", json::array{}}, {"val", json::array{0}}}},
           {"chop_digits", 0}}}});
  if (prepared.at("status") != "ok")
    throw std::runtime_error("singular chart.prepare: " +
                             json::serialize(prepared));
  return std::string(prepared.at("chart").as_string());
}

json::object regular_run(const std::string& value) {
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= 4; ++n) {
    shifts.emplace_back(std::to_string(n));
    schedule.push_back(json::array{json::object{
        {"case", n == 0 ? "R" : "T"}, {"da", std::to_string(n)},
        {"db", "0"}}});
  }
  return json::object{
      {"nmax", 4}, {"p", 0}, {"has_initial", true},
      {"adaptive_probe", false}, {"a_target", "0"},
      {"b_target", "0"}, {"a_shift_min", 0},
      {"a_shifts", std::move(shifts)}, {"schedule", std::move(schedule)},
      {"initial", json::array{"0", value, "0", "0"}},
      {"initial_validity", json::array{2}}, {"source", nullptr},
      {"return_u", false}};
}

json::object singular_run() {
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= 4; ++n) {
    shifts.emplace_back(n == 0 ? "1/2"
                                : std::to_string(2 * n + 1) + "/2");
    schedule.push_back(json::array{json::object{
        {"case", n == 0 ? "R" : "T"}, {"da", std::to_string(n)},
        {"db", "0"}}});
  }
  return json::object{
      {"nmax", 4}, {"p", 0}, {"has_initial", true},
      {"adaptive_probe", false}, {"a_target", "1/2"},
      {"b_target", "1/3"}, {"a_shift_min", 0},
      {"a_shifts", std::move(shifts)}, {"schedule", std::move(schedule)},
      {"initial", json::array{"0", "1", "0", "0"}},
      {"initial_validity", json::array{2}}, {"source", nullptr},
      {"return_u", false}};
}

json::object metadata(const std::string& center, const std::string& checkpoint,
                      bool singular) {
  json::object tag{
      {"a", json::object{{"domain", "rational"},
                           {"canonical", singular ? "1/2" : "0"}}},
      {"b", json::object{{"domain", "rational"},
                           {"canonical", singular ? "1/3" : "0"}}}};
  if (singular)
    tag["p"] = json::object{{"domain", "integer"}, {"canonical", "0"}};
  return json::object{
      {"chart", json::object{{"center_exact", center},
                              {"scale_exact", "1"},
                              {"radius", "2"},
                              {"infinite_radius", false}}},
      {"tag", std::move(tag)},
      {"prescriptions", json::array{}},
      {"checkpoint_identity", checkpoint}};
}

std::string solve_regular(const std::string& session,
                          const std::string& chart,
                          const std::string& center,
                          const std::string& checkpoint,
                          const std::string& value) {
  const auto solved = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart}, {"run", regular_run(value)},
      {"metadata", metadata(center, checkpoint, false)}});
  if (solved.at("status") != "ok")
    throw std::runtime_error("regular local.solve: " +
                             json::serialize(solved));
  return std::string(solved.at("local").as_string());
}

std::string prepare_singular_scc(const std::string& session,
                                 const std::string& block_chart,
                                 const std::string& key,
                                 const std::string& identity) {
  json::array exact_system_row;
  exact_system_row.push_back(json::object{
      {"exact", "((1/2)+(eps/3))/t"}, {"proven_zero", false}});
  json::array exact_system;
  exact_system.push_back(std::move(exact_system_row));
  json::array exact_theta_row;
  exact_theta_row.push_back(json::object{
      {"exact", "(1/2)+(eps/3)"}, {"proven_zero", false}});
  json::array exact_theta;
  exact_theta.push_back(std::move(exact_theta_row));
  const auto prepared = request(json::object{
      {"schema", 2}, {"op", "scc.prepare"}, {"session", session},
      {"key", key}, {"identity", identity},
      {"parent", json::object{
           {"dimension", 1},
           {"exact_system_record", std::move(exact_system)},
           {"exact_theta_record", std::move(exact_theta)},
           {"chart", geometry("-2/3")},
           {"scc", json::object{
                {"components", singleton_components()},
                {"structural_edges", self_edges()},
                {"condensation_edges", json::array{}},
                {"topological_order", json::array{0}},
                {"coupling_depth", 0}}},
           {"execution", json::object{
                {"mode", "BlockSequentialStrict"},
                {"work_t_order", 4}}},
           {"work_contract", json::object{
                {"work_min", -1}, {"requested_min", -1},
                {"requested_max", 2}, {"work_complete_max", 2},
                {"public_t_order", 0},
                {"wolfram_coupling_depth", 1}}}}},
      {"blocks", json::array{json::object{
           {"block", 0}, {"vertices", json::array{0}},
           {"chart", block_chart},
           {"principal_identity", "singular-arm-block-identity"},
           {"regular", false}, {"identity_gauge", true},
           {"identity_v", true}, {"no_pseudo", true},
           {"exact_affine_jordan_indicial", json::object{
                {"schema", "diffexp2-exact-affine-jordan-indicial-v1"},
                {"dimension", 1},
                {"blocks", json::array{json::object{
                     {"block", 0}, {"columns", json::array{0}},
                     {"a", "1/2"}, {"b", "1/3"}}}}}}}}},
      {"couplings", json::array{}}});
  if (prepared.at("status") != "ok")
    throw std::runtime_error("scc.prepare: " + json::serialize(prepared));
  return std::string(prepared.at("scc").as_string());
}

std::string solve_singular_column(const std::string& session,
                                  const std::string& scc,
                                  const std::string& checkpoint) {
  const auto solved = request(json::object{
      {"schema", 2}, {"op", "scc.solve_column"}, {"session", session},
      {"scc", scc}, {"checkpoint_identity", checkpoint},
      {"seed", json::object{
           {"block", 0}, {"run", singular_run()},
           {"metadata", metadata("-2/3", checkpoint + ":seed", true)}}},
      {"targets", json::array{}}});
  if (solved.at("status") != "ok")
    throw std::runtime_error("scc.solve_column: " +
                             json::serialize(solved));
  if (solved.at("execution_capability") !=
      "acb-regular-singular-scalar-block-dag-column-v1")
    throw std::runtime_error("singular SCC lacks scalar capability: " +
                             json::serialize(solved));
  return std::string(solved.at("local").as_string());
}

json::object topology(bool singular) {
  return json::object{
      {"singular_points",
       singular ? json::array{"-2/3"} : json::array{}},
      {"boundary_points", json::array{}},
      {"complex_projections", json::array{}},
      {"branch_sheets", json::array{}}};
}

json::object arm(const std::string& endpoint, const std::string& anchor,
                 const std::string& receiver, bool singular) {
  return json::object{{"from_exact", "0"}, {"to_exact", endpoint},
                      {"charts", json::array{anchor, receiver}},
                      {"topology", topology(singular)}};
}

json::object zero_match_arm(const std::string& endpoint,
                            const std::string& anchor) {
  return json::object{{"from_exact", "0"}, {"to_exact", endpoint},
                      {"charts", json::array{anchor}},
                      {"topology", topology(false)}};
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
                {"kernels", json::array{one, zero, zero, zero}},
                {"exact_identity", identity},
                {"proven_zero", false}}}}}}};
}

json::object execution(const std::string& basis,
                       const std::string& prefix) {
  json::array basis_set;
  basis_set.emplace_back(basis);
  json::array receiving_basis;
  receiving_basis.push_back(std::move(basis_set));
  return json::object{
      {"receiving_basis", std::move(receiving_basis)},
      {"integrand_rows", json::array{
           integrand_row(prefix + ":row:0"),
           integrand_row(prefix + ":row:1")}}};
}

json::object run_arms(const std::string& session, const std::string& plan,
                      const std::string& anchor,
                      const std::string& lower_basis,
                      const std::string& upper_basis,
                      const std::string& checkpoint_root) {
  return request(json::object{
      {"schema", 2}, {"op", "integration.run_arms"},
      {"session", session}, {"tile_plan", plan}, {"anchor", anchor},
      {"tile_plan_checkpoint_identity", "singular-arm-plan"},
      {"anchor_checkpoint_identity", "singular-arm-anchor"},
      {"epsilon", json::object{
           {"min", -1}, {"max", 2}, {"required_complete_max", 1},
           {"match_required_complete_max", 2}}},
      {"refinement", json::object{
           {"relative_tolerance", "1e-30"}, {"max_steps", 2}}},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
           {"root", checkpoint_root}}},
      {"lower", execution(lower_basis, checkpoint_root + ":lower")},
      {"upper", execution(upper_basis, checkpoint_root + ":upper")}});
}

json::object run_transport_arm(
    const std::string& session, const std::string& plan,
    const std::string& plan_checkpoint, const std::string& anchor,
    const std::string& basis, const std::string& checkpoint_root) {
  json::array basis_set;
  basis_set.emplace_back(basis);
  json::array receiving_basis;
  receiving_basis.push_back(std::move(basis_set));
  return request(json::object{
      {"schema", 2}, {"op", "transport.run_arm"},
      {"session", session}, {"tile_plan", plan}, {"anchor", anchor},
      {"tile_plan_checkpoint_identity", plan_checkpoint},
      {"anchor_checkpoint_identity", "singular-arm-anchor"},
      {"arm", "lower"},
      {"receiving_basis", std::move(receiving_basis)},
      {"epsilon", json::object{
           {"min", -1}, {"max", 2}, {"required_complete_max", 1},
           {"match_required_complete_max", 2}}},
      {"refinement", json::object{
           {"relative_tolerance", "1e-30"}, {"max_steps", 2}}},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
           {"root", checkpoint_root}}}});
}

json::object run_zero_match_transport_arm(
    const std::string& session, const std::string& plan,
    const std::string& anchor, const std::string& checkpoint_root) {
  return request(json::object{
      {"schema", 2}, {"op", "transport.run_arm"},
      {"session", session}, {"tile_plan", plan}, {"anchor", anchor},
      {"tile_plan_checkpoint_identity", "zero-state-plan"},
      {"anchor_checkpoint_identity", "singular-arm-anchor"},
      {"arm", "lower"}, {"receiving_basis", json::array{}},
      {"epsilon", json::object{
           {"min", -1}, {"max", 2}, {"required_complete_max", 1},
           {"match_required_complete_max", 2}}},
      {"refinement", json::object{
           {"relative_tolerance", "1e-30"}, {"max_steps", 2}}},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
           {"root", checkpoint_root}}}});
}

json::object observable(const std::string& identity,
                        const std::string& checkpoint_identity,
                        const std::string& tail_policy = "stored",
                        bool malformed_second_row = false,
                        bool shifted_output_window = false) {
  json::array rows{
      integrand_row(identity + ":tile:0"),
      integrand_row(identity + ":tile:1")};
  if (malformed_second_row)
    rows[1].as_object()["columns"] = 2;
  if (shifted_output_window)
    for (auto& raw_row : rows)
      raw_row.as_object().at("entries").as_array().front().as_object()
          .at("multiplier").as_object()["epsilon_shift"] = 1;
  return json::object{
      {"identity", identity},
      {"checkpoint_identity", checkpoint_identity},
      {"integrand_rows", std::move(rows)},
      {"epsilon", json::object{{"min", shifted_output_window ? 0 : -1},
                                {"max", shifted_output_window ? 3 : 2},
                                {"required_complete_max", 1}}},
      {"tail_policy", tail_policy}};
}

json::object contract_observables(
    const std::string& session, const std::string& state,
    const std::string& state_checkpoint,
    const std::string& state_provenance, json::array observables,
    const std::string& checkpoint_root) {
  return request(json::object{
      {"schema", 2}, {"op", "transport.contract"},
      {"session", session}, {"transport_state", state},
      {"transport_state_checkpoint_identity", state_checkpoint},
      {"transport_state_provenance_identity", state_provenance},
      {"checkpoint_policy", json::object{
           {"schema",
            "diffexp2-deterministic-transport-contraction-checkpoints-v1"},
           {"root", checkpoint_root}}},
      {"observables", std::move(observables)}});
}

json::object line_stats(const std::string& session,
                        const std::string& line) {
  return request(json::object{{"schema", 2}, {"op", "integration.stats"},
                              {"session", session}, {"line", line}});
}

json::object line_export(const std::string& session, const std::string& line,
                         const std::string& checkpoint_identity) {
  return request(json::object{
      {"schema", 2}, {"op", "integration.export"},
      {"session", session}, {"line", line},
      {"checkpoint_identity", checkpoint_identity},
      {"output_digits", 30}});
}

json::object direct_singular_match_request(
    const std::string& session, const std::string& scc,
    const std::string& basis, const std::string& incoming) {
  return request(json::object{
      {"schema", 2}, {"op", "local.match_acb"}, {"session", session},
      {"basis", json::array{basis}}, {"incoming", incoming},
      {"basis_chart", scc}, {"incoming_chart", scc},
      {"basis_point", json::object{{"exact", "1/2"}}},
      {"incoming_point", json::object{{"exact", "1/2"}}},
      {"epsilon", json::object{{"min", -1}, {"max", 2},
                                {"required_complete_max", 2}}},
      {"basis_checkpoint_identities",
       json::array{"singular-arm-good-column"}},
      {"incoming_checkpoint_identity", "singular-arm-second-column"},
      {"checkpoint_identity", "singular-arm-direct-match"},
      {"native_singular_scc_saturation", json::object{
           {"schema",
            "diffexp2-native-acb-singular-scc-valuation-zero-saturation-request-v1"},
           {"session_configuration_identity", "caller-config"},
           {"tile_plan", "tile:caller-string"},
           {"tile_plan_checkpoint_identity", "caller-string"},
           {"tile_plan_provenance_identity", "caller-string"},
           {"arm", "lower"}, {"match", 0},
           {"match_checkpoint_identity", "singular-arm-direct-match"},
           {"receiving_scc", scc},
           {"receiving_scc_exact_identity", "singular-arm-scc-identity"},
           {"receiving_execution_capability",
            "acb-regular-singular-scalar-block-dag-column-v1"},
           {"receiving_basis_point_exact", "1/2"},
           {"physical_match_point_exact", "-1/2"},
           {"receiving_rim", nullptr}}},
      {"refinement", json::object{
           {"relative_tolerance", "1e-30"}, {"max_steps", 2}}}});
}

bool tamper_singular_proof_session(const std::string& source,
                                   const std::string& target) {
  auto container = diffexp2::checkpoint::read(source);
  auto payload = json::parse(container.payload_json).as_object();
  bool changed = false;
  for (auto& raw_hop :
       payload.at("retained_planned_match_hops").as_array()) {
    auto& hop = raw_hop.as_object();
    auto& native = hop.at("native_match").as_object();
    auto witness = json::parse(std::string(
        native.at("exact_lattice_canonical_witness").as_string()));
    if (witness.as_object().at("schema") != kSingularProofSchema)
      continue;
    witness.as_object().at("native_request").as_object()[
        "session_configuration_identity"] = "tampered-config";
    native["exact_lattice_canonical_witness"] = json::serialize(witness);
    changed = true;
    break;
  }
  if (!changed) return false;
  diffexp2::checkpoint::write_atomic(
      target, container.header_json, json::serialize(payload));
  return true;
}

std::pair<bool, bool> proof_schemas(const std::string& path) {
  const auto payload = json::parse(
      diffexp2::checkpoint::read(path).payload_json).as_object();
  bool singular = false;
  bool ordinary = false;
  for (const auto& raw_hop :
       payload.at("retained_planned_match_hops").as_array()) {
    const auto& native = raw_hop.as_object().at("native_match").as_object();
    const auto witness = json::parse(std::string(
        native.at("exact_lattice_canonical_witness").as_string()));
    const auto schema = std::string(
        witness.as_object().at("schema").as_string());
    singular = singular || schema == kSingularProofSchema;
    ordinary = ordinary || schema == kOrdinaryProofSchema;
  }
  return {singular, ordinary};
}

}  // namespace

int main() {
  const std::string checkpoint =
      "/tmp/diffexp2-singular-scc-arms.de2cp";
  const std::string checkpoint_second =
      "/tmp/diffexp2-singular-scc-arms-second.de2cp";
  const std::string tampered =
      "/tmp/diffexp2-singular-scc-arms-tampered.de2cp";
  std::remove(checkpoint.c_str());
  std::remove(checkpoint_second.c_str());
  std::remove(tampered.c_str());
  std::string session;
  std::string restored_session;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"}, {"domain", "acb"},
        {"precision_bits", 256}, {"output_digits", 30},
        {"chart_capacity", 6}, {"scc_capacity", 2},
        {"local_capacity", 32}, {"match_capacity", 8},
        {"tile_plan_capacity", 2}});
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_regular_chart(
        session, "singular-arm-anchor-chart", "0");
    const auto singular_block = prepare_singular_chart(session);
    const auto upper_chart = prepare_regular_chart(
        session, "singular-arm-upper-chart", "2/3");
    const auto singular_scc = prepare_singular_scc(
        session, singular_block, "singular-arm-scc",
        "singular-arm-scc-identity");
    const auto foreign_scc = prepare_singular_scc(
        session, singular_block, "singular-arm-foreign-scc",
        "singular-arm-foreign-scc-identity");
    const auto anchor = solve_regular(
        session, anchor_chart, "0", "singular-arm-anchor", "2");
    const auto good_singular = solve_singular_column(
        session, singular_scc, "singular-arm-good-column");
    const auto second_singular = solve_singular_column(
        session, singular_scc, "singular-arm-second-column");
    const auto foreign_singular = solve_singular_column(
        session, foreign_scc, "singular-arm-foreign-column");
    const auto upper_basis = solve_regular(
        session, upper_chart, "2/3", "singular-arm-upper-basis", "1");
    const auto direct_rejected = direct_singular_match_request(
        session, singular_scc, good_singular, second_singular);
    if (direct_rejected.at("status") != "error" ||
        std::string(direct_rejected.at("detail").as_string()).find(
            "only through a retained planned match") == std::string::npos)
      throw std::runtime_error(
          "direct singular-SCC proof admission was not rejected: " +
          json::serialize(direct_rejected));

    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", "singular-arm-plan"},
        {"division_order", 3},
        {"lower", arm("-2/3", anchor_chart, singular_scc, true)},
        {"upper", arm("2/3", anchor_chart, upper_chart, false)}});
    if (planned.at("status") != "ok")
      throw std::runtime_error("tile.plan: " + json::serialize(planned));
    const auto plan = std::string(planned.at("tile_plan").as_string());
    if (planned.at("lower").as_object().at("matches").as_array().size() != 1 ||
        planned.at("upper").as_object().at("matches").as_array().size() != 1)
      throw std::runtime_error(
          "fixture did not produce one match per arm: " +
          json::serialize(planned));

    const auto before_rejection = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto rejected = run_arms(
        session, plan, anchor, foreign_singular, upper_basis,
        "singular-arm-rejected");
    const auto after_rejection = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const bool rejected_atomically = rejected.at("status") == "error" &&
        std::string(rejected.at("detail").as_string()).find(
            "basis chart provenance mismatch") != std::string::npos &&
        before_rejection.at("locals") == after_rejection.at("locals") &&
        before_rejection.at("matches") == after_rejection.at("matches") &&
        before_rejection.at("line_results") ==
            after_rejection.at("line_results") &&
        after_rejection.at("pending_local_solves") == 0 &&
        after_rejection.at("pending_matches") == 0 &&
        after_rejection.at("pending_line_integrations") == 0;
    if (!rejected_atomically)
      throw std::runtime_error(
          "foreign singular-SCC basis was not rejected atomically: " +
          json::serialize(rejected) + " / " +
          json::serialize(after_rejection));

    const auto marched = run_arms(
        session, plan, anchor, good_singular, upper_basis,
        "singular-arm-success");
    if (marched.at("status") != "ok")
      throw std::runtime_error("integration.run_arms: " +
                               json::serialize(marched));

    const auto single_planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan_arm"}, {"session", session},
        {"checkpoint_identity", "singular-state-plan"},
        {"division_order", 3},
        {"arm", arm("-2/3", anchor_chart, singular_scc, true)}});
    if (single_planned.at("status") != "ok" ||
        single_planned.at("arm_name") != "lower" ||
        single_planned.at("matches") != 1)
      throw std::runtime_error(
          "single-arm tile.plan_arm: " + json::serialize(single_planned));
    const auto single_plan =
        std::string(single_planned.at("tile_plan").as_string());

    const auto before_transport_rejection = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto transport_rejected = run_transport_arm(
        session, single_plan, "singular-state-plan", anchor,
        foreign_singular, "singular-state-rejected");
    const auto after_transport_rejection = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const bool transport_rejected_atomically =
        transport_rejected.at("status") == "error" &&
        std::string(transport_rejected.at("detail").as_string()).find(
            "basis chart provenance mismatch") != std::string::npos &&
        before_transport_rejection.at("locals") ==
            after_transport_rejection.at("locals") &&
        before_transport_rejection.at("matches") ==
            after_transport_rejection.at("matches") &&
        before_transport_rejection.at("transport_states") ==
            after_transport_rejection.at("transport_states") &&
        after_transport_rejection.at("pending_local_solves") == 0 &&
        after_transport_rejection.at("pending_matches") == 0 &&
        after_transport_rejection.at("pending_transport_states") == 0;
    if (!transport_rejected_atomically)
      throw std::runtime_error(
          "foreign singular transport basis was not rejected atomically: " +
          json::serialize(transport_rejected) + " / " +
          json::serialize(after_transport_rejection));

    const auto transported = run_transport_arm(
        session, single_plan, "singular-state-plan", anchor,
        good_singular, "singular-state-success");
    if (transported.at("status") != "ok" ||
        transported.at("matches") != 1 || transported.at("tiles") != 2 ||
        transported.at("atomic_publication") != true)
      throw std::runtime_error(
          "transport.run_arm: " + json::serialize(transported));
    const auto transport_state =
        std::string(transported.at("transport_state").as_string());
    const auto transport_checkpoint = std::string(
        transported.at("checkpoint_identity").as_string());
    const auto transport_provenance = std::string(
        transported.at("provenance_identity").as_string());
    const auto transport_final = std::string(
        transported.at("final_local").as_object().at("local").as_string());
    const auto transport_stats = request(json::object{
        {"schema", 2}, {"op", "transport.stats"}, {"session", session},
        {"transport_state", transport_state}});
    if (transport_stats.at("status") != "ok" ||
        transport_stats.at("final_local").as_object().at("local").as_string() !=
            transport_final)
      throw std::runtime_error(
          "transport.stats: " + json::serialize(transport_stats));

    std::vector<std::pair<std::string, std::string>> contracted_lines;
    json::value compact_export_value;
    const auto before_empty_contract = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto before_empty_plan = request(json::object{
        {"schema", 2}, {"op", "tile.stats"}, {"session", session},
        {"tile_plan", single_plan}});
    const auto empty_contract = contract_observables(
        session, transport_state, transport_checkpoint,
        transport_provenance, json::array{}, "contract-empty");
    const auto after_empty_contract = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto after_empty_plan = request(json::object{
        {"schema", 2}, {"op", "tile.stats"}, {"session", session},
        {"tile_plan", single_plan}});
    const auto after_empty_state = request(json::object{
        {"schema", 2}, {"op", "transport.stats"}, {"session", session},
        {"transport_state", transport_state}});
    if (empty_contract.at("status") != "ok" ||
        empty_contract.at("observables") != 0 ||
        !empty_contract.at("lines").as_array().empty() ||
        before_empty_contract.at("locals") != after_empty_contract.at("locals") ||
        before_empty_contract.at("matches") != after_empty_contract.at("matches") ||
        before_empty_contract.at("line_results") !=
            after_empty_contract.at("line_results") ||
        before_empty_contract.at("line_integrations") !=
            after_empty_contract.at("line_integrations") ||
        unsigned_value(after_empty_contract, "transport_contractions") !=
            unsigned_value(before_empty_contract, "transport_contractions") + 1 ||
        before_empty_contract.at("transport_observables") !=
            after_empty_contract.at("transport_observables") ||
        before_empty_plan.at("integrations") != after_empty_plan.at("integrations") ||
        unsigned_value(after_empty_state, "contraction_operations") != 1 ||
        unsigned_value(after_empty_state, "contracted_observables") != 0)
      throw std::runtime_error(
          "zero-observable contraction was not a handle-free batch no-op: " +
          json::serialize(empty_contract) + " / " +
          json::serialize(after_empty_contract) + " / " +
          json::serialize(after_empty_state));

    json::array one_observable;
    one_observable.push_back(observable(
        "observable-one", "observable-one-checkpoint"));
    const auto before_one_contract = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto one_contract = contract_observables(
        session, transport_state, transport_checkpoint,
        transport_provenance, std::move(one_observable), "contract-one");
    const auto after_one_contract = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    if (one_contract.at("status") != "ok" ||
        one_contract.at("observables") != 1 ||
        one_contract.at("json_coefficients") != 0 ||
        one_contract.at("no_rematching") != true ||
        one_contract.at("lines").as_array().size() != 1 ||
        before_one_contract.at("locals") != after_one_contract.at("locals") ||
        before_one_contract.at("matches") != after_one_contract.at("matches") ||
        unsigned_value(after_one_contract, "line_results") !=
            unsigned_value(before_one_contract, "line_results") + 1 ||
        unsigned_value(after_one_contract, "line_integrations") !=
            unsigned_value(before_one_contract, "line_integrations") + 2 ||
        before_one_contract.at("local_matches") !=
            after_one_contract.at("local_matches"))
      throw std::runtime_error(
          "one-observable contraction did not publish one compact line: " +
          json::serialize(one_contract) + " / " +
          json::serialize(after_one_contract));
    const auto& one_line = one_contract.at("lines").as_array().front().as_object();
    if (one_line.at("request_index") != 0 ||
        one_line.at("observable_identity") != "observable-one" ||
        one_line.at("session").as_string() != session)
      throw std::runtime_error(
          "one-observable contraction lost ordered public identity");
    contracted_lines.emplace_back(
        std::string(one_line.at("line").as_string()),
        std::string(one_line.at("checkpoint_identity").as_string()));

    json::array three_observables;
    for (int index = 0; index < 3; ++index)
      three_observables.push_back(observable(
          "observable-three-" + std::to_string(index),
          "observable-three-checkpoint-" + std::to_string(index),
          index == 1 ? "attempt" : "stored", false, index == 2));
    const auto before_three_contract = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto three_contract = contract_observables(
        session, transport_state, transport_checkpoint,
        transport_provenance, std::move(three_observables), "contract-three");
    const auto after_three_contract = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    if (three_contract.at("status") != "ok" ||
        three_contract.at("observables") != 3 ||
        three_contract.at("lines").as_array().size() != 3 ||
        before_three_contract.at("locals") != after_three_contract.at("locals") ||
        before_three_contract.at("matches") != after_three_contract.at("matches") ||
        unsigned_value(after_three_contract, "line_results") !=
            unsigned_value(before_three_contract, "line_results") + 3 ||
        unsigned_value(after_three_contract, "line_integrations") !=
            unsigned_value(before_three_contract, "line_integrations") + 6 ||
        before_three_contract.at("local_matches") !=
            after_three_contract.at("local_matches"))
      throw std::runtime_error(
          "three-observable contraction did not publish three compact lines: " +
          json::serialize(three_contract) + " / " +
          json::serialize(after_three_contract));
    for (std::size_t index = 0;
         index < three_contract.at("lines").as_array().size(); ++index) {
      const auto& line =
          three_contract.at("lines").as_array()[index].as_object();
      if (unsigned_value(line, "request_index") != index ||
          line.at("observable_identity").as_string() !=
              "observable-three-" + std::to_string(index) ||
          line.at("session").as_string() != session)
        throw std::runtime_error(
            "three-observable contraction reordered its outputs");
      contracted_lines.emplace_back(
          std::string(line.at("line").as_string()),
          std::string(line.at("checkpoint_identity").as_string()));
    }
    if (three_contract.at("lines").as_array()[1].as_object().at("scope") !=
            "stored_truncation" ||
        unsigned_value(
            three_contract.at("lines").as_array()[2].as_object()
                .at("epsilon").as_object(),
            "max") != 2 ||
        three_contract.at("streaming_tile_contraction") != true)
      throw std::runtime_error(
          "attempt-tail downgrade or shifted output epsilon contract was not preserved: " +
          json::serialize(three_contract));
    const auto after_success_state = request(json::object{
        {"schema", 2}, {"op", "transport.stats"}, {"session", session},
        {"transport_state", transport_state}});
    if (unsigned_value(after_success_state, "contraction_operations") != 3 ||
        unsigned_value(after_success_state, "contracted_observables") != 4)
      throw std::runtime_error(
          "successful contraction counters are inconsistent: " +
          json::serialize(after_success_state));

    json::array malformed_observables;
    malformed_observables.push_back(observable(
        "malformed-prefix-good", "malformed-prefix-good-checkpoint"));
    malformed_observables.push_back(observable(
        "malformed-row", "malformed-row-checkpoint", "stored", true));
    const auto before_malformed = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto before_malformed_state = request(json::object{
        {"schema", 2}, {"op", "transport.stats"}, {"session", session},
        {"transport_state", transport_state}});
    const auto before_malformed_plan = request(json::object{
        {"schema", 2}, {"op", "tile.stats"}, {"session", session},
        {"tile_plan", single_plan}});
    const auto malformed_contract = contract_observables(
        session, transport_state, transport_checkpoint,
        transport_provenance, std::move(malformed_observables),
        "contract-malformed");
    const auto after_malformed = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto after_malformed_state = request(json::object{
        {"schema", 2}, {"op", "transport.stats"}, {"session", session},
        {"transport_state", transport_state}});
    const auto after_malformed_plan = request(json::object{
        {"schema", 2}, {"op", "tile.stats"}, {"session", session},
        {"tile_plan", single_plan}});
    if (malformed_contract.at("status") != "error" ||
        std::string(malformed_contract.at("detail").as_string()).find(
            "dimension differs") == std::string::npos ||
        before_malformed.at("locals") != after_malformed.at("locals") ||
        before_malformed.at("matches") != after_malformed.at("matches") ||
        before_malformed.at("line_results") != after_malformed.at("line_results") ||
        before_malformed.at("line_integrations") !=
            after_malformed.at("line_integrations") ||
        before_malformed.at("transport_contractions") !=
            after_malformed.at("transport_contractions") ||
        before_malformed.at("transport_observables") !=
            after_malformed.at("transport_observables") ||
        before_malformed_state.at("contraction_operations") !=
            after_malformed_state.at("contraction_operations") ||
        before_malformed_state.at("contracted_observables") !=
            after_malformed_state.at("contracted_observables") ||
        before_malformed_plan.at("integrations") !=
            after_malformed_plan.at("integrations") ||
        after_malformed.at("pending_line_integrations") != 0)
      throw std::runtime_error(
          "malformed observable contraction was not rolled back atomically: " +
          json::serialize(malformed_contract) + " / " +
          json::serialize(after_malformed));

    json::array require_observable;
    require_observable.push_back(observable(
        "require-tail", "require-tail-checkpoint", "require"));
    const auto before_require = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto before_require_state = request(json::object{
        {"schema", 2}, {"op", "transport.stats"}, {"session", session},
        {"transport_state", transport_state}});
    const auto before_require_plan = request(json::object{
        {"schema", 2}, {"op", "tile.stats"}, {"session", session},
        {"tile_plan", single_plan}});
    const auto require_contract = contract_observables(
        session, transport_state, transport_checkpoint,
        transport_provenance, std::move(require_observable),
        "contract-require-tail");
    const auto after_require = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto after_require_state = request(json::object{
        {"schema", 2}, {"op", "transport.stats"}, {"session", session},
        {"transport_state", transport_state}});
    const auto after_require_plan = request(json::object{
        {"schema", 2}, {"op", "tile.stats"}, {"session", session},
        {"tile_plan", single_plan}});
    if (require_contract.at("status") != "error" ||
        std::string(require_contract.at("detail").as_string()).find(
            "requires a certified full-local tail") == std::string::npos ||
        before_require.at("locals") != after_require.at("locals") ||
        before_require.at("line_results") != after_require.at("line_results") ||
        before_require.at("line_integrations") !=
            after_require.at("line_integrations") ||
        before_require.at("transport_contractions") !=
            after_require.at("transport_contractions") ||
        before_require.at("transport_observables") !=
            after_require.at("transport_observables") ||
        before_require_state.at("contraction_operations") !=
            after_require_state.at("contraction_operations") ||
        before_require_state.at("contracted_observables") !=
            after_require_state.at("contracted_observables") ||
        before_require_plan.at("integrations") !=
            after_require_plan.at("integrations") ||
        after_require.at("pending_line_integrations") != 0)
      throw std::runtime_error(
          "required-tail contraction did not fail atomically and honestly: " +
          json::serialize(require_contract) + " / " +
          json::serialize(after_require));

    const auto released_plan = request(json::object{
        {"schema", 2}, {"op", "tile.release"}, {"session", session},
        {"tile_plan", single_plan}});
    const auto released_final = request(json::object{
        {"schema", 2}, {"op", "local.release"}, {"session", session},
        {"local", transport_final}});
    if (released_plan.at("status") != "ok" ||
        released_final.at("status") != "ok")
      throw std::runtime_error(
          "transport closure dependency release failed: " +
          json::serialize(released_plan) + " / " +
          json::serialize(released_final));
    const auto hidden_plan = request(json::object{
        {"schema", 2}, {"op", "tile.stats"}, {"session", session},
        {"tile_plan", single_plan}});
    if (hidden_plan.at("status") != "error")
      throw std::runtime_error(
          "released transport-owned tile plan remained public: " +
          json::serialize(hidden_plan));

    const auto zero_planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan_arm"}, {"session", session},
        {"checkpoint_identity", "zero-state-plan"}, {"division_order", 3},
        {"arm", zero_match_arm("-1/3", anchor_chart)}});
    if (zero_planned.at("status") != "ok" ||
        zero_planned.at("arm_name") != "lower" ||
        zero_planned.at("matches") != 0 || zero_planned.at("tiles") != 1)
      throw std::runtime_error(
          "zero-match tile.plan_arm: " + json::serialize(zero_planned));
    const auto zero_plan =
        std::string(zero_planned.at("tile_plan").as_string());
    const auto before_zero_transport = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    const auto zero_transport = run_zero_match_transport_arm(
        session, zero_plan, anchor, "zero-state-success");
    const auto after_zero_transport = request(json::object{
        {"schema", 2}, {"op", "session.stats"}, {"session", session}});
    if (zero_transport.at("status") != "ok" ||
        zero_transport.at("matches") != 0 || zero_transport.at("tiles") != 1 ||
        zero_transport.at("final_local").as_object().at("local").as_string() !=
            anchor ||
        before_zero_transport.at("locals") != after_zero_transport.at("locals") ||
        before_zero_transport.at("matches") != after_zero_transport.at("matches") ||
        after_zero_transport.at("transport_states").as_int64() !=
            before_zero_transport.at("transport_states").as_int64() + 1)
      throw std::runtime_error(
          "zero-match transport did not retain the anchor without publishing a local/match: " +
          json::serialize(zero_transport) + " / " +
          json::serialize(after_zero_transport));
    const auto zero_state =
        std::string(zero_transport.at("transport_state").as_string());
    const auto released_zero_state = request(json::object{
        {"schema", 2}, {"op", "transport.release"}, {"session", session},
        {"transport_state", zero_state}});
    const auto released_zero_plan = request(json::object{
        {"schema", 2}, {"op", "tile.release"}, {"session", session},
        {"tile_plan", zero_plan}});
    if (released_zero_state.at("status") != "ok" ||
        released_zero_plan.at("status") != "ok")
      throw std::runtime_error(
          "zero-match transport release failed: " +
          json::serialize(released_zero_state) + " / " +
          json::serialize(released_zero_plan));

    const auto released_transport = request(json::object{
        {"schema", 2}, {"op", "transport.release"}, {"session", session},
        {"transport_state", transport_state}});
    const auto hidden_transport = request(json::object{
        {"schema", 2}, {"op", "transport.stats"}, {"session", session},
        {"transport_state", transport_state}});
    const auto compact_stats = line_stats(
        session, contracted_lines.front().first);
    const auto compact_export = line_export(
        session, contracted_lines.front().first,
        contracted_lines.front().second);
    if (released_transport.at("status") != "ok" ||
        hidden_transport.at("status") != "error" ||
        compact_stats.at("status") != "ok" ||
        compact_export.at("status") != "ok" ||
        compact_export.at("json_coefficients") == 0)
      throw std::runtime_error(
          "compact observable line did not survive transport-state release: " +
          json::serialize(released_transport) + " / " +
          json::serialize(compact_stats) + " / " +
          json::serialize(compact_export));
    compact_export_value = compact_export.at("value");

    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
        {"path", checkpoint},
        {"checkpoint_identity", "singular-arm-roundtrip"}});
    if (saved.at("status") != "ok")
      throw std::runtime_error("checkpoint.save: " +
                               json::serialize(saved));
    const auto saved_payload = json::parse(
        diffexp2::checkpoint::read(checkpoint).payload_json).as_object();
    const auto& saved_visibility = saved_payload.at("session").as_object()
        .at("registry_visibility").as_object();
    const auto contains_string = [](const json::array& values,
                                    const std::string& expected) {
      for (const auto& value : values)
        if (value.is_string() && value.as_string() == expected) return true;
      return false;
    };
    const auto contains_record = [](const json::array& values,
                                    const char* key,
                                    const std::string& expected) {
      for (const auto& value : values)
        if (value.is_object() && value.as_object().at(key).as_string() ==
                                     expected)
          return true;
      return false;
    };
    if (saved_payload.at("retained_transport_states").as_array().size() != 1)
      throw std::runtime_error(
          "checkpoint did not retain exactly one hidden transport state");
    const auto& retained_state_record =
        saved_payload.at("retained_transport_states").as_array().front()
            .as_object();
    const auto& retained_state_provenance =
        retained_state_record.at("provenance").as_object();
    const auto& retained_state_plan =
        retained_state_provenance.at("tile_plan").as_object();
    bool recursive_match_provenance = false;
    for (const auto& raw_match :
         retained_state_provenance.at("matches").as_array())
      recursive_match_provenance |=
          raw_match.as_object().if_contains("provenance_identity") != nullptr;
    if (!contains_record(
            saved_payload.at("retained_transport_states").as_array(),
            "handle", transport_state) ||
        contains_string(saved_visibility.at("transport_states").as_array(),
                        transport_state) ||
        contains_string(saved_visibility.at("tile_plans").as_array(),
                        single_plan) ||
        contains_string(saved_visibility.at("locals").as_array(),
                        transport_final) ||
        !contains_record(saved_payload.at("retained_tile_plans").as_array(),
                         "handle", single_plan) ||
        !contains_record(saved_payload.at("retained_locals").as_array(),
                         "handle", transport_final) ||
        retained_state_record.at("schema") !=
            "diffexp2-retained-transport-arm-state-v3" ||
        retained_state_plan.if_contains("provenance_identity") != nullptr ||
        recursive_match_provenance ||
        unsigned_value(retained_state_record.at("runtime_stats").as_object(),
                       "contraction_operations") != 3 ||
        unsigned_value(retained_state_record.at("runtime_stats").as_object(),
                       "contracted_observables") != 4)
      throw std::runtime_error(
          "checkpoint did not separate transport-owned closure from registry visibility");
    for (const auto& [line_handle, ignored] : contracted_lines) {
      bool found_compact = false;
      for (const auto& raw_line :
           saved_payload.at("retained_line_results").as_array()) {
        const auto& line = raw_line.as_object();
        if (line.at("handle").as_string() != line_handle) continue;
        found_compact =
            line.at("schema") ==
                "diffexp2-retained-transport-observable-line-v2" &&
            line.at("provenance").as_object().at("aggregate").as_object()
                    .if_contains("components") == nullptr;
      }
      if (!found_compact)
        throw std::runtime_error(
            "transport observable checkpoint retained projected component state");
    }
    for (const auto& raw_local :
         saved_payload.at("retained_locals").as_array())
      if (std::string(raw_local.as_object().at("handle").as_string())
              .starts_with("private:"))
        throw std::runtime_error(
            "transport observable checkpoint retained a private projected local");
    const auto [singular_proof, ordinary_proof] = proof_schemas(checkpoint);

    const auto restored = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint},
        {"expected_identity", "singular-arm-roundtrip"}});
    if (restored.at("status") == "ok")
      restored_session = std::string(restored.at("session").as_string());
    if (restored.at("status") != "ok")
      throw std::runtime_error("checkpoint.restore: " +
                               json::serialize(restored));
    const auto restored_transport_stats = request(json::object{
        {"schema", 2}, {"op", "transport.stats"},
        {"session", restored_session},
        {"transport_state", transport_state}});
    const auto restored_hidden_plan = request(json::object{
        {"schema", 2}, {"op", "tile.stats"},
        {"session", restored_session}, {"tile_plan", single_plan}});
    const auto restored_compact_stats = line_stats(
        restored_session, contracted_lines.front().first);
    const auto restored_compact_export = line_export(
        restored_session, contracted_lines.front().first,
        contracted_lines.front().second);
    if (restored_transport_stats.at("status") != "error" ||
        restored_hidden_plan.at("status") != "error" ||
        restored_compact_stats.at("status") != "ok" ||
        restored_compact_export.at("status") != "ok" ||
        restored_compact_export.at("value") != compact_export_value)
      throw std::runtime_error(
          "hidden transport-owned compact line did not restore exactly: " +
          json::serialize(restored_transport_stats) + " / " +
          json::serialize(restored_hidden_plan) + " / " +
          json::serialize(restored_compact_export));
    const auto saved_second = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", restored_session}, {"path", checkpoint_second},
        {"checkpoint_identity", "singular-arm-roundtrip-second"}});
    if (saved_second.at("status") != "ok")
      throw std::runtime_error("second checkpoint.save: " +
                               json::serialize(saved_second));
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", restored_session}});
    restored_session.clear();

    const auto restored_second = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint_second},
        {"expected_identity", "singular-arm-roundtrip-second"}});
    if (restored_second.at("status") == "ok")
      restored_session =
          std::string(restored_second.at("session").as_string());
    if (restored_second.at("status") != "ok")
      throw std::runtime_error("second checkpoint.restore: " +
                               json::serialize(restored_second));
    const auto second_transport_stats = request(json::object{
        {"schema", 2}, {"op", "transport.stats"},
        {"session", restored_session},
        {"transport_state", transport_state}});
    const auto second_compact_stats = line_stats(
        restored_session, contracted_lines.front().first);
    const auto second_compact_export = line_export(
        restored_session, contracted_lines.front().first,
        contracted_lines.front().second);
    const auto after_transport_release = request(json::object{
        {"schema", 2}, {"op", "session.stats"},
        {"session", restored_session}});
    if (second_transport_stats.at("status") != "error" ||
        second_compact_stats.at("status") != "ok" ||
        second_compact_export.at("status") != "ok" ||
        second_compact_export.at("value") != compact_export_value ||
        after_transport_release.at("transport_states") != 0)
      throw std::runtime_error(
          "second-generation compact transport line was not stable: " +
          json::serialize(second_transport_stats) + " / " +
          json::serialize(second_compact_export) + " / " +
          json::serialize(after_transport_release));
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", restored_session}});
    restored_session.clear();

    if (!tamper_singular_proof_session(checkpoint_second, tampered))
      throw std::runtime_error("checkpoint contains no singular-SCC proof");
    const auto tamper_rejected = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", tampered},
        {"expected_identity", "singular-arm-roundtrip-second"}});
    const auto& arms = marched.at("arms").as_object();
    const bool ok = singular_proof && ordinary_proof &&
        arms.at("lower").as_object().at("matches") == 1 &&
        arms.at("upper").as_object().at("matches") == 1 &&
        tamper_rejected.at("status") == "error" &&
        std::string(tamper_rejected.at("detail").as_string()).find(
            "session-configuration binding") != std::string::npos;
    if (!ok)
      std::cerr << "marched: " << json::serialize(marched) << '\n'
                << "schemas: singular=" << singular_proof
                << " ordinary=" << ordinary_proof << '\n'
                << "tamper: " << json::serialize(tamper_rejected) << '\n';
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    session.clear();
    std::remove(checkpoint.c_str());
    std::remove(checkpoint_second.c_str());
    std::remove(tampered.c_str());
    std::cout << (ok ? "PASS" : "FAIL")
              << ": Acb singular-SCC valuation-zero whole-arm proof\n";
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    if (!restored_session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored_session}});
    std::remove(checkpoint.c_str());
    std::remove(checkpoint_second.c_str());
    std::remove(tampered.c_str());
    std::cerr << "FAIL: Acb singular-SCC whole-arm proof: "
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
