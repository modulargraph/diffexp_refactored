#include "diffexp2/checkpoint.hpp"
#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

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

json::object metadata(const std::string& checkpoint_identity,
                      const std::string& a = "0") {
  return json::object{
      {"chart", json::object{{"center_exact", "0"},
                              {"scale_exact", "1"},
                              {"radius", "2"},
                              {"infinite_radius", false}}},
      {"tag", json::object{
          {"a", json::object{{"domain", "rational"},
                              {"canonical", a}}},
          {"b", json::object{{"domain", "rational"},
                              {"canonical", "0"}}}}},
      {"prescriptions", json::array{}},
      {"checkpoint_identity", checkpoint_identity}};
}

json::object regular_run(std::uint32_t nmax, const std::string& a = "0") {
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    const auto affine_a = n == 0 ? a : std::to_string(n);
    shifts.emplace_back(affine_a);
    schedule.push_back(json::array{json::object{
        {"case", n == 0 && a == "0" ? "R" : "T"},
        {"da", affine_a},
        {"db", "0"}}});
  }
  return json::object{
      {"nmax", nmax}, {"p", 0}, {"has_initial", true},
      {"adaptive_probe", false}, {"a_target", a},
      {"b_target", "0"}, {"a_shift_min", 0},
      {"a_shifts", std::move(shifts)},
      {"schedule", std::move(schedule)},
      {"initial", json::array{"1"}},
      {"initial_validity", json::array{0}},
      {"source", nullptr}, {"return_u", false}};
}

json::object residual_request(const std::string& session,
                              const std::string& local,
                              const json::object& binding) {
  return json::object{
      {"schema", 2}, {"op", "local.certify_residual"},
      {"session", session}, {"local", local},
      {"point", json::object{{"exact", "1/10"}}},
      {"options", json::object{{"tail_estimate", false}}},
      {"relative_tolerance", "1e-18"},
      {"scope", "stored_truncation"}, {"include_residual", false},
      {"operator_identity", binding.at("operator_identity")},
      {"source_identity", binding.at("source_identity")},
      {"checkpoint_identity", binding.at("local_checkpoint_identity")},
      {"analytic_metadata", binding.at("analytic_metadata")},
      {"provenance_identity", binding.at("provenance_identity")}};
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

json::object checkpoint_payload(const std::filesystem::path& path) {
  return json::parse(diffexp2::checkpoint::read(path.string()).payload_json)
      .as_object();
}

const json::object& local_record(const json::object& payload,
                                 const std::string& handle) {
  for (const auto& raw : payload.at("retained_locals").as_array()) {
    const auto& record = raw.as_object();
    if (std::string(record.at("handle").as_string()) == handle)
      return record;
  }
  throw std::runtime_error("checkpoint lost retained local " + handle);
}

}  // namespace

int main() {
  const auto stem = std::filesystem::temp_directory_path() /
      ("diffexp2_owner_bound_residual_" +
       std::to_string(static_cast<long long>(::getpid())));
  const auto saved_path = stem.string() + "_saved.de2cp";
  const auto resaved_path = stem.string() + "_resaved.de2cp";

  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"},
      {"domain", "rational"}, {"precision_bits", 256},
      {"output_digits", 40}, {"chart_capacity", 2},
      {"local_capacity", 4}});
  require_ok(created, "session.create");
  const auto session = std::string(created.at("session").as_string());
  const auto prepared = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", "owner-bound-exp@0"},
      {"identity", "owner-bound-exp-operator-v1"},
      {"analytic", json::object{{"prescription", "+i0"}}},
      {"scc", json::object{
          {"components", json::array{json::value(json::array{0})}},
          {"structural_edges", json::array{}},
          {"condensation_edges", json::array{}},
          {"topological_order", json::array{0}},
          {"coupling_depth", 0}}},
      {"problem", json::object{
          {"domain", "rational"}, {"d", 1}, {"fb", 0}, {"w", 1},
          {"d_lags", json::array{json::value(json::array{
              json::object{{"s", 0}, {"v", "1"}}})}},
          {"denominators", json::array{}},
          {"nhat_lags", json::array{
              json::object{{"poly", json::array{}},
                           {"rat", json::array{}},
                           {"val", json::array{nullptr}}},
              json::object{{"poly", json::array{json::object{
                               {"s", 0},
                               {"e", json::array{json::value(
                                   json::array{0, 0, "1"})}}}}},
                           {"rat", json::array{}},
                           {"val", json::array{0}}}}},
          {"d0_inverse", "1"},
          {"blocks", json::array{json::value(json::array{0})}},
          {"assembly", json::object{{"identity", true},
                                      {"poly", json::array{}},
                                      {"rat", json::array{}},
                                      {"val", json::array{0}}}},
          {"chop_digits", 0}}}});
  require_ok(prepared, "chart.prepare");
  const auto chart = std::string(prepared.at("chart").as_string());

  const auto solved = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart}, {"run", regular_run(12)},
      {"metadata", metadata("owner-bound-exp-local-v1")}});
  require_ok(solved, "regular local.solve");
  const auto local = std::string(solved.at("local").as_string());
  const auto& residual_status = solved.at("residual_binding").as_object();
  const auto binding = residual_status.at("binding").as_object();

  // A fractional-power local deliberately lies outside this first theorem
  // slice.  Its unsupported reason must survive a checkpoint round trip.
  const auto unsupported_solved = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart}, {"run", regular_run(0, "1/2")},
      {"metadata", metadata("owner-bound-unsupported-local-v1", "1/2")}});
  require_ok(unsupported_solved, "unsupported local.solve");
  const auto unsupported =
      std::string(unsupported_solved.at("local").as_string());

  const auto released = request(json::object{
      {"schema", 2}, {"op", "chart.release"}, {"session", session},
      {"chart", chart}});
  require_ok(released, "chart.release");
  const auto valid = request(residual_request(session, local, binding));
  require_ok(valid, "owner-bound residual after chart release");
  auto full_scope_request = residual_request(session, local, binding);
  full_scope_request["scope"] = "full_local_solution";
  const auto full_scope = request(full_scope_request);
  require_ok(full_scope, "full-local owner-bound residual");

  auto operator_payload_tamper = residual_request(session, local, binding);
  operator_payload_tamper["theta_operator"] = json::object{
      {"min", 0}, {"max", 0}, {"dimension", 1},
      {"coefficients", json::array{"0"}}};
  const auto operator_payload_rejected = request(operator_payload_tamper);

  auto source_payload_tamper = residual_request(session, local, binding);
  source_payload_tamper["source"] = json::object{
      {"epsilon", json::object{{"min", 0}, {"max", 0}}},
      {"dimension", 1}, {"coefficients", json::array{"0"}},
      {"error", json::object{}}};
  const auto source_payload_rejected = request(source_payload_tamper);

  auto operator_identity_tamper = residual_request(session, local, binding);
  operator_identity_tamper["operator_identity"] = "attacker-operator";
  const auto operator_identity_rejected = request(operator_identity_tamper);

  auto source_identity_tamper = residual_request(session, local, binding);
  source_identity_tamper["source_identity"] = "attacker-source";
  const auto source_identity_rejected = request(source_identity_tamper);

  auto provenance_tamper = residual_request(session, local, binding);
  provenance_tamper["provenance_identity"] = "attacker-provenance";
  const auto provenance_rejected = request(provenance_tamper);

  auto analytic_tamper = residual_request(session, local, binding);
  analytic_tamper.at("analytic_metadata").as_object()
      .at("chart").as_object()["center_exact"] = "1";
  const auto analytic_rejected = request(analytic_tamper);

  auto checkpoint_tamper = residual_request(session, local, binding);
  checkpoint_tamper["checkpoint_identity"] = "attacker-checkpoint";
  const auto checkpoint_rejected = request(checkpoint_tamper);

  const auto saved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
      {"path", saved_path},
      {"checkpoint_identity", "owner-bound-residual-session-v1"}});
  require_ok(saved, "checkpoint.save");
  const auto original_payload = checkpoint_payload(saved_path);
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});

  const auto restored = request(json::object{
      {"schema", 2}, {"op", "checkpoint.restore"},
      {"path", saved_path},
      {"expected_identity", "owner-bound-residual-session-v1"}});
  require_ok(restored, "checkpoint.restore");
  const auto restored_session =
      std::string(restored.at("session").as_string());
  const auto restored_valid = request(
      residual_request(restored_session, local, binding));
  require_ok(restored_valid, "owner-bound residual after restore");
  const auto resaved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"},
      {"session", restored_session}, {"path", resaved_path},
      {"checkpoint_identity", "owner-bound-residual-resaved-v1"}});
  require_ok(resaved, "checkpoint resave");
  const auto restored_payload = checkpoint_payload(resaved_path);

  const auto& original_bound = local_record(original_payload, local);
  const auto& restored_bound = local_record(restored_payload, local);
  const auto& original_unsupported =
      local_record(original_payload, unsupported);
  const auto& restored_unsupported =
      local_record(restored_payload, unsupported);
  const auto& saved_binding = original_bound.at(
      "residual_operator_restore").as_object();
  const auto& saved_unsupported = original_unsupported.at(
      "residual_operator_restore").as_object();

  const bool ok = created.at("status") == "ok" &&
      prepared.at("status") == "ok" && solved.at("status") == "ok" &&
      residual_status.at("status") == "available" &&
      unsupported_solved.at("status") == "ok" &&
      unsupported_solved.at("residual_binding").as_object().at("status") ==
          "unsupported" &&
      released.at("status") == "ok" && valid.at("status") == "ok" &&
      valid.at("verdict") == "pass" &&
      full_scope.at("status") == "ok" &&
      full_scope.at("scope") == "full_local_solution" &&
      full_scope.at("verdict") == "inconclusive" &&
      valid.at("operator_identity") == binding.at("operator_identity") &&
      valid.at("source_identity") == binding.at("source_identity") &&
      valid.at("provenance_identity") == binding.at("provenance_identity") &&
      valid.at("analytic_metadata") == binding.at("analytic_metadata") &&
      is_error(operator_payload_rejected) &&
      is_error(source_payload_rejected) &&
      is_error(operator_identity_rejected) &&
      is_error(source_identity_rejected) && is_error(provenance_rejected) &&
      is_error(analytic_rejected) && is_error(checkpoint_rejected) &&
      saved.at("status") == "ok" && restored.at("status") == "ok" &&
      restored_valid.at("status") == "ok" &&
      restored_valid.at("verdict") == "pass" &&
      resaved.at("status") == "ok" &&
      original_bound.at("schema") == "diffexp2-retained-local-v3" &&
      saved_binding.at("serialized") == true &&
      saved_binding.at("status") == "available" &&
      saved_binding.at("binding") == residual_status.at("binding") &&
      original_bound.at("residual_operator_restore") ==
          restored_bound.at("residual_operator_restore") &&
      saved_unsupported.at("serialized") == false &&
      saved_unsupported.at("status") == "unsupported" &&
      original_unsupported.at("residual_operator_restore") ==
          restored_unsupported.at("residual_operator_restore");

  if (!ok) {
    std::cerr << "solve: " << json::serialize(solved) << '\n'
              << "unsupported solve: "
              << json::serialize(unsupported_solved) << '\n'
              << "valid after release: " << json::serialize(valid) << '\n'
              << "full-local scope: " << json::serialize(full_scope) << '\n'
              << "operator payload tamper: "
              << json::serialize(operator_payload_rejected) << '\n'
              << "source payload tamper: "
              << json::serialize(source_payload_rejected) << '\n'
              << "operator identity tamper: "
              << json::serialize(operator_identity_rejected) << '\n'
              << "source identity tamper: "
              << json::serialize(source_identity_rejected) << '\n'
              << "provenance tamper: "
              << json::serialize(provenance_rejected) << '\n'
              << "analytic tamper: " << json::serialize(analytic_rejected)
              << '\n' << "checkpoint tamper: "
              << json::serialize(checkpoint_rejected) << '\n'
              << "restore: " << json::serialize(restored) << '\n'
              << "valid after restore: " << json::serialize(restored_valid)
              << '\n';
  }

  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", restored_session}});
  std::error_code error;
  std::filesystem::remove(saved_path, error);
  std::filesystem::remove(resaved_path, error);
  std::cout << (ok ? "PASS" : "FAIL")
            << ": owner-bound persistent residual and checkpoint restore\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
