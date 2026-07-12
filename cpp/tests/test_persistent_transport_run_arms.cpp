#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

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

json::array prescriptions() {
  return json::array{json::object{
      {"factor_exact", "paired-state-f"}, {"sign", -1},
      {"multiplicity", 1}, {"leading_coefficient_sign", 1}}};
}

std::string prepare_chart(const std::string& session,
                          const std::string& name,
                          const std::string& center) {
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
  const auto prepared = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", name}, {"identity", name + ":identity"},
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
           {"domain", "rational"}, {"precision_bits", 256},
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
           {"chop_digits", 0}}}});
  require_ok(prepared, "chart.prepare");
  return std::string(prepared.at("chart").as_string());
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
  const auto solved = request(json::object{
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
  require_ok(solved, "local.solve");
  return std::string(solved.at("local").as_string());
}

json::object topology() {
  return json::object{
      {"singular_points", json::array{}},
      {"boundary_points", json::array{}},
      {"complex_projections", json::array{}},
      {"branch_sheets", json::array{json::object{
           {"factor_exact", "paired-state-f"}, {"sign", -1}}}}};
}

json::object arm(const std::string& endpoint,
                 const std::string& anchor_chart,
                 const std::string& receiving_chart) {
  return json::object{
      {"from_exact", "0"}, {"to_exact", endpoint},
      {"charts", json::array{anchor_chart, receiving_chart}},
      {"topology", topology()}};
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
