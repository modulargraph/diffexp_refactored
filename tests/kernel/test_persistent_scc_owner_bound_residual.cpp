#include "diffexp/kernel/checkpoint.hpp"
#include "diffexp/kernel/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <string>

#include <unistd.h>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp::kernel::run_recurrence_json(json::serialize(value)))
      .as_object();
}

bool is_error(const json::object& response) {
  return response.if_contains("status") != nullptr &&
         response.at("status") == "error";
}

void require_ok(const json::object& response, const char* label) {
  if (response.if_contains("status") == nullptr ||
      response.at("status") != "ok")
    throw std::runtime_error(std::string(label) + ": " +
                             json::serialize(response));
}

json::value nested(json::array value) {
  return json::value(std::move(value));
}

json::object epsilon_rational(std::int32_t valuation,
                              std::initializer_list<const char*> numerator,
                              std::initializer_list<const char*> denominator) {
  json::array p;
  json::array q;
  for (const auto* value : numerator) p.emplace_back(value);
  for (const auto* value : denominator) q.emplace_back(value);
  return json::object{{"zero", false},
                      {"valuation", valuation},
                      {"numerator", std::move(p)},
                      {"denominator", std::move(q)}};
}

json::object physical_entry(std::uint32_t row, std::uint32_t column,
                            json::object value) {
  return json::object{{"r", row}, {"c", column}, {"v", std::move(value)}};
}

json::object parent_physical_ode(const std::string& owner) {
  const auto one = epsilon_rational(0, {"1"}, {"1"});
  const auto diagonal_one = epsilon_rational(0, {"1"}, {"1", "-1"});
  const auto diagonal_two = epsilon_rational(0, {"2"}, {"1", "-1"});
  const auto cross = epsilon_rational(1, {"1"}, {"1", "-1"});
  return json::object{
      {"schema", "diffexp3-physical-cleared-ode-v1"},
      {"basis", "physical-original-master"},
      {"theta_coordinate", "local-t"},
      {"owner_signature_identity", owner},
      {"payload_identity", "de2-physical-ode-scc-full-parent-v1"},
      {"q", json::array{one}},
      {"c", json::array{
          nested(json::array{}),
          nested(json::array{
              physical_entry(0, 0, diagonal_one),
              physical_entry(1, 0, cross),
              physical_entry(1, 1, diagonal_two)})}}};
}

json::object exact_cell(const char* exact, bool zero) {
  return json::object{{"exact", exact}, {"proven_zero", zero}};
}

json::object scalar_problem(const char* diagonal,
                            bool nontrivial_assembly) {
  json::object assembly;
  if (nontrivial_assembly) {
    assembly = json::object{
        {"identity", false},
        {"poly", json::array{json::object{
            {"s", 0},
            {"e", json::array{nested(json::array{0, 0, "2"})}}}}},
        {"rat", json::array{}}, {"val", json::array{0}}};
  } else {
    assembly = json::object{{"identity", true},
                            {"poly", json::array{}},
                            {"rat", json::array{}},
                            {"val", json::array{0}}};
  }
  return json::object{
      {"domain", "rational"}, {"d", 1}, {"fb", 0}, {"w", 8},
      {"d_lags", json::array{nested(json::array{
          json::object{{"s", 0}, {"v", "1"}}})}},
      {"denominators", json::array{nested(json::array{"1", "-1"})}},
      {"nhat_lags", json::array{
          json::object{{"poly", json::array{}},
                       {"rat", json::array{}},
                       {"val", json::array{nullptr}}},
          json::object{{"poly", json::array{}},
                       {"rat", json::array{json::object{
                           {"q", 0},
                           {"num", json::array{json::object{
                               {"s", 0},
                               {"e", json::array{nested(
                                   json::array{0, 0, diagonal})}}}}}}}},
                       {"val", json::array{0}}}}},
      {"d0_inverse", "1"},
      {"blocks", json::array{nested(json::array{0})}},
      {"assembly", std::move(assembly)}, {"chop_digits", 0}};
}

json::object prepare_scalar_chart(const std::string& session,
                                  const std::string& key,
                                  const std::string& identity,
                                  const char* diagonal,
                                  bool nontrivial_assembly = false) {
  return request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", key}, {"identity", identity},
      {"analytic", json::object{
          {"geometry", json::object{
              {"center_exact", "0"}, {"scale_exact", "1"},
              {"radius_exact", "2"}, {"infinite_radius", false},
              {"prescriptions", json::array{}}}},
          {"principal_matrix", json::array{nested(json::array{
              exact_cell(diagonal[0] == '1' ? "1/(1-eps)"
                                               : "2/(1-eps)", false)})}},
          {"native_scc_capabilities", json::object{
              {"regular", true}, {"identity_gauge", true},
              {"identity_v", !nontrivial_assembly},
              {"no_pseudo", true}}}}},
      {"scc", json::object{
          {"components", json::array{nested(json::array{0})}},
          {"structural_edges", json::array{nested(json::array{0, 0})}},
          {"condensation_edges", json::array{}},
          {"topological_order", json::array{0}}, {"coupling_depth", 0}}},
      {"problem", scalar_problem(diagonal, nontrivial_assembly)}});
}

json::array coupling_kernels() {
  json::array kernels;
  for (std::size_t epsilon = 0; epsilon < 8; ++epsilon) {
    json::array taylor;
    for (std::size_t n = 0; n <= 24; ++n)
      taylor.emplace_back(n == 1 ? "1" : "0");
    kernels.push_back(nested(std::move(taylor)));
  }
  return kernels;
}

json::object scc_manifest(const std::string& first,
                          const std::string& second,
                          const std::string& identity) {
  return json::object{
      {"identity", identity},
      {"parent", json::object{
          {"dimension", 2},
          {"exact_system_record", json::array{
              nested(json::array{exact_cell("1/(1-eps)", false),
                                  exact_cell("0", true)}),
              nested(json::array{exact_cell("eps/(1-eps)", false),
                                  exact_cell("2/(1-eps)", false)})}},
          {"exact_theta_record", json::array{
              nested(json::array{exact_cell("t/(1-eps)", false),
                                  exact_cell("0", true)}),
              nested(json::array{exact_cell("eps*t/(1-eps)", false),
                                  exact_cell("2*t/(1-eps)", false)})}},
          {"chart", json::object{
              {"center_exact", "0"}, {"scale_exact", "1"},
              {"radius_exact", "2"}, {"infinite_radius", false},
              {"prescriptions", json::array{}}}},
          {"scc", json::object{
              {"components", json::array{nested(json::array{0}),
                                          nested(json::array{1})}},
              {"structural_edges", json::array{
                  nested(json::array{0, 0}), nested(json::array{0, 1}),
                  nested(json::array{1, 1})}},
              {"condensation_edges", json::array{
                  nested(json::array{0, 1})}},
              {"topological_order", json::array{0, 1}},
              {"coupling_depth", 1}}},
          {"execution", json::object{
              {"mode", "BlockSequentialStrict"}, {"work_t_order", 24}}},
          {"work_contract", json::object{
              {"work_min", 0}, {"requested_min", 0},
              {"requested_max", 3}, {"work_complete_max", 7},
              {"public_t_order", 18}, {"wolfram_coupling_depth", 2}}}}},
      {"blocks", json::array{
          json::object{{"block", 0}, {"vertices", json::array{0}},
                       {"chart", first},
                       {"principal_identity", "scc-rational-diagonal-one"},
                       {"regular", true}, {"identity_gauge", true},
                       {"identity_v", true}, {"no_pseudo", true}},
          json::object{{"block", 1}, {"vertices", json::array{1}},
                       {"chart", second},
                       {"principal_identity", "scc-rational-diagonal-two"},
                       {"regular", true}, {"identity_gauge", true},
                       {"identity_v", true}, {"no_pseudo", true}}}},
      {"couplings", json::array{json::object{
          {"source_block", 0}, {"target_block", 1},
          {"source_vertices", json::array{0}},
          {"target_vertices", json::array{1}},
          {"rows", 1}, {"columns", 1},
          {"exact_identity", "scc-full-parent-cross-v1"},
          {"domain", "rational"}, {"symbols", json::array{}},
          {"entries", json::array{json::object{
              {"row", 0}, {"column", 0},
              {"source_vertex", 0}, {"target_vertex", 1},
              {"exact_original_entry", "eps/(1-eps)"},
              {"exact_theta_entry", "eps*t/(1-eps)"},
              {"multiplier", json::object{
                  {"epsilon_shift", 1}, {"center_pole_order", 0},
                  {"kernels", coupling_kernels()},
                  {"exact_identity", "eps*t/(1-eps)"},
                  {"proven_zero", false}}}}}}}}},
      {"physical_ode", parent_physical_ode(identity)}};
}

json::object prepare_scc(const std::string& session, const std::string& key,
                         json::object manifest) {
  manifest["schema"] = 2;
  manifest["op"] = "scc.prepare";
  manifest["session"] = session;
  manifest["key"] = key;
  return request(std::move(manifest));
}

json::object regular_run(bool seed) {
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= 24; ++n) {
    shifts.emplace_back(std::to_string(n));
    schedule.push_back(nested(json::array{json::object{
        {"case", n == 0 ? "R" : "T"},
        {"da", std::to_string(n)}, {"db", "0"}}}));
  }
  json::array initial;
  for (std::size_t epsilon = 0; epsilon < 8; ++epsilon)
    initial.emplace_back(seed && epsilon == 0 ? "1" : "0");
  return json::object{
      {"nmax", 24}, {"p", 0}, {"has_initial", seed},
      {"adaptive_probe", false}, {"a_target", "0"}, {"b_target", "0"},
      {"a_shift_min", 0}, {"a_shifts", std::move(shifts)},
      {"schedule", std::move(schedule)}, {"initial", std::move(initial)},
      {"initial_validity", json::array{seed ? json::value(7)
                                             : json::value(nullptr)}},
      {"source", nullptr}, {"return_u", false}};
}

json::object metadata(const std::string& checkpoint) {
  return json::object{
      {"chart", json::object{{"center_exact", "0"},
                              {"scale_exact", "1"}, {"radius", "2"},
                              {"infinite_radius", false}}},
      {"tag", json::object{
          {"a", json::object{{"domain", "rational"}, {"canonical", "0"}}},
          {"b", json::object{{"domain", "rational"}, {"canonical", "0"}}}}},
      {"prescriptions", json::array{}},
      {"checkpoint_identity", checkpoint}};
}

json::object residual_request(const std::string& session,
                              const std::string& local,
                              const json::object& binding) {
  return json::object{
      {"schema", 2}, {"op", "local.certify_residual"},
      {"session", session}, {"local", local},
      {"point", json::object{{"exact", "1/10"}}},
      {"options", json::object{{"tail_estimate", false}}},
      {"relative_tolerance", "2e-24"},
      {"scope", "stored_truncation"}, {"include_residual", true},
      {"operator_identity", binding.at("operator_identity")},
      {"source_identity", binding.at("source_identity")},
      {"checkpoint_identity", binding.at("local_checkpoint_identity")},
      {"analytic_metadata", binding.at("analytic_metadata")},
      {"owner_signature_identity", binding.at("owner_signature_identity")},
      {"physical_payload_identity", binding.at("physical_payload_identity")},
      {"provenance_identity", binding.at("provenance_identity")}};
}

json::object checkpoint_payload(const std::filesystem::path& path) {
  return json::parse(diffexp::kernel::checkpoint::read(path.string()).payload_json)
      .as_object();
}

const json::object& record_by_handle(const json::object& payload,
                                     const char* collection,
                                     const char* key,
                                     const std::string& handle) {
  for (const auto& raw : payload.at(collection).as_array()) {
    const auto& record = raw.as_object();
    if (std::string(record.at(key).as_string()) == handle) return record;
  }
  throw std::runtime_error(std::string("checkpoint lost ") + collection +
                           " handle " + handle);
}

json::object& record_by_handle(json::object& payload,
                               const char* collection,
                               const char* key,
                               const std::string& handle) {
  for (auto& raw : payload.at(collection).as_array()) {
    auto& record = raw.as_object();
    if (std::string(record.at(key).as_string()) == handle) return record;
  }
  throw std::runtime_error(std::string("checkpoint lost ") + collection +
                           " handle " + handle);
}

}  // namespace

int main() {
  const auto stem = std::filesystem::temp_directory_path() /
      ("diffexp2_scc_physical_owner_" +
       std::to_string(static_cast<long long>(::getpid())));
  const auto saved_path = stem.string() + "_saved.de2cp";
  const auto tampered_path = stem.string() + "_tampered.de2cp";
  const auto resaved_path = stem.string() + "_resaved.de2cp";

  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", "rational"},
      {"precision_bits", 256}, {"output_digits", 50},
      {"chart_capacity", 4}, {"scc_capacity", 4}, {"local_capacity", 5}});
  require_ok(created, "session.create");
  const auto session = std::string(created.at("session").as_string());

  const auto first_prepared = prepare_scalar_chart(
      session, "scc-rational-diagonal-one-key",
      "scc-rational-diagonal-one", "1");
  const auto second_prepared = prepare_scalar_chart(
      session, "scc-rational-diagonal-two-key",
      "scc-rational-diagonal-two", "2");
  const auto nontrivial_prepared = prepare_scalar_chart(
      session, "scc-rational-nontrivial-v-key",
      "scc-rational-nontrivial-v", "1", true);
  require_ok(first_prepared, "first chart.prepare");
  require_ok(second_prepared, "second chart.prepare");
  require_ok(nontrivial_prepared, "V!=I chart.prepare");
  const auto first = std::string(first_prepared.at("chart").as_string());
  const auto second = std::string(second_prepared.at("chart").as_string());
  const auto nontrivial =
      std::string(nontrivial_prepared.at("chart").as_string());

  const std::string owner = "scc-full-parent-physical-owner-v1";
  const auto prepared = prepare_scc(
      session, "scc-full-parent-owner-key",
      scc_manifest(first, second, owner));
  require_ok(prepared, "scc.prepare");
  const auto scc = std::string(prepared.at("scc").as_string());

  auto wrong_owner_manifest = scc_manifest(
      first, second, "scc-full-parent-wrong-owner-v1");
  wrong_owner_manifest.at("physical_ode").as_object()
      ["owner_signature_identity"] = "attacker-owner";
  const auto wrong_owner = prepare_scc(
      session, "scc-full-parent-wrong-owner-key",
      std::move(wrong_owner_manifest));

  auto nontrivial_manifest = scc_manifest(
      nontrivial, second, "scc-full-parent-nontrivial-v-v1");
  auto& nontrivial_block =
      nontrivial_manifest.at("blocks").as_array().front().as_object();
  nontrivial_block["principal_identity"] = "scc-rational-nontrivial-v";
  nontrivial_block["identity_v"] = false;
  const auto nontrivial_rejected = prepare_scc(
      session, "scc-full-parent-nontrivial-v-key",
      std::move(nontrivial_manifest));

  const auto intermediate = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", first}, {"run", regular_run(true)},
      {"metadata", metadata("scc-intermediate-diagonal-v1")}});
  require_ok(intermediate, "standalone diagonal local.solve");

  const auto solved = request(json::object{
      {"schema", 2}, {"op", "scc.solve_column"}, {"session", session},
      {"scc", scc}, {"checkpoint_identity", "scc-full-parent-local-v1"},
      {"seed", json::object{{"block", 0}, {"run", regular_run(true)},
                              {"metadata", metadata("scc-seed-work-v1")}}},
      {"targets", json::array{json::object{
          {"block", 1}, {"run", regular_run(false)},
          {"metadata", metadata("scc-target-work-v1")}}}}});
  require_ok(solved, "scc.solve_column");
  const auto local = std::string(solved.at("local").as_string());
  const auto local_stats = request(json::object{
      {"schema", 2}, {"op", "local.stats"}, {"session", session},
      {"local", local}});
  require_ok(local_stats, "full-parent local.stats");
  const auto binding = local_stats.at("residual_binding").as_object()
      .at("binding").as_object();
  const auto certified = request(residual_request(session, local, binding));
  require_ok(certified, "full-parent residual");

  json::array batch_columns;
  batch_columns.push_back(json::object{
      {"checkpoint_identity", "scc-full-parent-batch-local-v1"},
      {"seed", json::object{{"block", 0}, {"run", regular_run(true)},
                              {"metadata", metadata("scc-batch-seed-v1")}}},
      {"targets", json::array{json::object{
          {"block", 1}, {"run", regular_run(false)},
          {"metadata", metadata("scc-batch-target-v1")}}}}});
  const auto solved_batch = request(json::object{
      {"schema", 2}, {"op", "scc.solve_columns"},
      {"session", session}, {"scc", scc},
      {"columns", std::move(batch_columns)}, {"threads", 1}});
  require_ok(solved_batch, "scc.solve_columns");
  const auto& batch_column =
      solved_batch.at("results").as_array().front().as_object();
  const auto& batch_residual =
      batch_column.at("residual_binding").as_object();
  const auto& batch_reference =
      batch_residual.at("binding_reference").as_object();
  const auto batch_local =
      std::string(batch_column.at("local").as_string());

  auto owner_tamper = residual_request(session, local, binding);
  owner_tamper["owner_signature_identity"] = "attacker-owner";
  const auto owner_tamper_rejected = request(std::move(owner_tamper));
  auto payload_tamper = residual_request(session, local, binding);
  payload_tamper["physical_payload_identity"] = "attacker-payload";
  const auto payload_tamper_rejected = request(std::move(payload_tamper));

  const auto intermediate_local =
      std::string(intermediate.at("local").as_string());
  require_ok(request(json::object{{"schema", 2}, {"op", "local.release"},
                                  {"session", session},
                                  {"local", intermediate_local}}),
             "intermediate local.release");
  require_ok(request(json::object{{"schema", 2}, {"op", "local.release"},
                                  {"session", session},
                                  {"local", batch_local}}),
             "batch local.release");
  require_ok(request(json::object{{"schema", 2}, {"op", "scc.release"},
                                  {"session", session}, {"scc", scc}}),
             "scc.release");
  for (const auto& chart : {first, second, nontrivial})
    require_ok(request(json::object{{"schema", 2}, {"op", "chart.release"},
                                    {"session", session}, {"chart", chart}}),
               "chart.release");
  const auto after_release = request(residual_request(session, local, binding));
  require_ok(after_release, "residual after hidden SCC release");

  const auto saved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
      {"path", saved_path},
      {"checkpoint_identity", "scc-full-parent-session-v1"}});
  require_ok(saved, "checkpoint.save");
  const auto original_payload = checkpoint_payload(saved_path);
  const auto& original_local = record_by_handle(
      original_payload, "retained_locals", "handle", local);
  const auto& original_scc = record_by_handle(
      original_payload, "prepared_scc", "handle", scc);
  const auto& original_physical = original_scc.at("request").as_object()
      .at("physical_ode");

  auto tampered_payload = original_payload;
  auto& tampered_scc = record_by_handle(
      tampered_payload, "prepared_scc", "handle", scc);
  tampered_scc.at("request").as_object().at("physical_ode").as_object()
      .at("c").as_array()[1].as_array()[0].as_object()
      .at("v").as_object().at("numerator").as_array()[0] = "9";
  const auto saved_container = diffexp::kernel::checkpoint::read(saved_path);
  diffexp::kernel::checkpoint::write_atomic(
      tampered_path, saved_container.header_json,
      json::serialize(tampered_payload));

  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  const auto tampered_restore = request(json::object{
      {"schema", 2}, {"op", "checkpoint.restore"},
      {"path", tampered_path},
      {"expected_identity", "scc-full-parent-session-v1"}});
  const auto restored = request(json::object{
      {"schema", 2}, {"op", "checkpoint.restore"}, {"path", saved_path},
      {"expected_identity", "scc-full-parent-session-v1"}});
  require_ok(restored, "checkpoint.restore");
  const auto restored_session =
      std::string(restored.at("session").as_string());
  const auto restored_certified = request(
      residual_request(restored_session, local, binding));
  require_ok(restored_certified, "residual after hidden SCC restore");
  const auto hidden_scc = request(json::object{
      {"schema", 2}, {"op", "scc.stats"}, {"session", restored_session},
      {"scc", scc}});
  const auto resaved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"},
      {"session", restored_session}, {"path", resaved_path},
      {"checkpoint_identity", "scc-full-parent-resaved-v1"}});
  require_ok(resaved, "checkpoint resave");
  const auto restored_payload = checkpoint_payload(resaved_path);
  const auto& restored_local = record_by_handle(
      restored_payload, "retained_locals", "handle", local);
  const auto& restored_scc = record_by_handle(
      restored_payload, "prepared_scc", "handle", scc);
  const auto& restored_physical = restored_scc.at("request").as_object()
      .at("physical_ode");

  const auto provenance = json::parse(
      std::string(binding.at("provenance_identity").as_string())).as_object();
  const auto& diagnostics = solved.at("block_diagnostics").as_array();
  const auto& owner_restore =
      original_local.at("equation_owner_restore").as_object();
  const auto registry_visibility = original_payload.at("session").as_object()
      .at("registry_visibility").as_object();
  const bool ok = is_error(wrong_owner) && is_error(nontrivial_rejected) &&
      std::string(nontrivial_rejected.at("detail").as_string())
              .find(
                  "exact t-independent Laurent-unimodular spectral frame") !=
          std::string::npos &&
      intermediate.at("residual_binding").as_object().at("status") ==
          "unsupported" &&
      prepared.at("physical_ode_owner") == "full-parent" &&
      prepared.at("physical_payload_identity") ==
          "de2-physical-ode-scc-full-parent-v1" &&
      solved.at("residual_binding").as_object().at("status") ==
          "available" &&
      solved.at("residual_binding").as_object().if_contains("binding") ==
          nullptr &&
      solved.at("residual_binding").as_object()
          .at("binding_reference").as_object().at("authority") ==
          "retained-native-exact-residual-owner" &&
      local_stats.at("residual_binding").as_object().at("status") ==
          "available" &&
      batch_residual.at("status") == "available" &&
      batch_residual.if_contains("binding") == nullptr &&
      batch_column.if_contains("tail_majorant") == nullptr &&
      batch_column.if_contains("metadata") == nullptr &&
      batch_column.if_contains("source_operator_identity") == nullptr &&
      batch_reference.at("schema") ==
          "diffexp3-owner-bound-residual-reference-v1" &&
      batch_reference.at("authority") ==
          "retained-native-exact-residual-owner" &&
      batch_reference.at("local_checkpoint_identity") ==
          "scc-full-parent-batch-local-v1" &&
      batch_reference.at("identity_diagnostics").as_object()
          .at("physical_payload_identity").as_object()
          .at("identity_bytes") ==
          std::string("de2-physical-ode-scc-full-parent-v1").size() &&
      diagnostics.size() == 2 &&
      diagnostics[1].as_object().at("role") == "particular" &&
      certified.at("verdict") == "pass" && certified.at("dimension") == 2 &&
      certified.at("epsilon_min") == 0 && certified.at("epsilon_max") == 7 &&
      provenance.at("owner_kind") == "composite-scc" &&
      std::string(provenance.at("owner_handle").as_string()) == scc &&
      is_error(owner_tamper_rejected) && is_error(payload_tamper_rejected) &&
      after_release.at("verdict") == "pass" &&
      original_payload.at("schema") == 9 &&
      registry_visibility.at("sccs").as_array().empty() &&
      owner_restore.at("owner_kind") == "composite-scc" &&
      std::string(owner_restore.at("owner_handle").as_string()) == scc &&
      original_physical == restored_physical &&
      original_local.at("equation_owner_restore") ==
          restored_local.at("equation_owner_restore") &&
      is_error(tampered_restore) && restored.at("sccs").as_array().empty() &&
      is_error(hidden_scc) && restored_certified.at("verdict") == "pass";

  if (!ok) {
    std::cerr << "wrong owner: " << json::serialize(wrong_owner) << '\n'
              << "V!=I: " << json::serialize(nontrivial_rejected) << '\n'
              << "intermediate: " << json::serialize(intermediate) << '\n'
              << "solved: " << json::serialize(solved) << '\n'
              << "solved batch: " << json::serialize(solved_batch) << '\n'
              << "certified: " << json::serialize(certified) << '\n'
              << "after release: " << json::serialize(after_release) << '\n'
              << "tampered restore: " << json::serialize(tampered_restore)
              << '\n' << "restored: " << json::serialize(restored) << '\n'
              << "restored certified: "
              << json::serialize(restored_certified) << '\n';
  }

  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", restored_session}});
  std::error_code error;
  std::filesystem::remove(saved_path, error);
  std::filesystem::remove(tampered_path, error);
  std::filesystem::remove(resaved_path, error);
  std::cout << (ok ? "PASS" : "FAIL")
            << ": CompositeSCC full-parent physical residual owner\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
