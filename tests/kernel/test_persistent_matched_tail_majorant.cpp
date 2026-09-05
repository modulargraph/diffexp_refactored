#include "diffexp/kernel/checkpoint.hpp"
#include "diffexp/kernel/json_codec.hpp"

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
  return json::parse(diffexp::kernel::run_recurrence_json(json::serialize(value)))
      .as_object();
}

json::array nested(json::array value) {
  json::array result;
  result.push_back(std::move(value));
  return result;
}

json::object prepared_problem(const std::string& domain) {
  json::object result{
      {"domain", domain}, {"d", 1}, {"fb", 0}, {"w", 1},
      {"d_lags", nested(json::array{
           json::object{{"s", 0}, {"v", "1"}}})},
      {"denominators", json::array{}},
      {"nhat_lags", json::array{
           json::object{{"poly", json::array{}}, {"rat", json::array{}},
                        {"val", json::array{nullptr}}},
           json::object{{"poly", json::array{json::object{
                                {"s", 0},
                                {"e", nested(json::array{0, 0, "1"})}}}},
                        {"rat", json::array{}}, {"val", json::array{0}}}}},
      {"d0_inverse", "1"}, {"blocks", nested(json::array{0})},
      {"assembly", json::object{
           {"identity", true}, {"poly", json::array{}},
           {"rat", json::array{}}, {"val", json::array{0}}}},
      {"chop_digits", 50}};
  if (domain == "acb") result["precision_bits"] = 256;
  return result;
}

std::string prepare_chart(const std::string& session,
                          const std::string& domain,
                          const std::string& key,
                          const std::string& center) {
  const auto identity = key + "-operator-v1";
  const auto response = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", key}, {"identity", identity},
      {"analytic", json::object{
           {"geometry", json::object{
                {"center_exact", center}, {"scale_exact", "1"},
                {"radius_exact", "2"}, {"infinite_radius", false},
                {"prescriptions", json::array{}}}},
           {"principal_matrix", nested(json::array{
                json::object{{"exact", "0"}, {"proven_zero", true}}})},
           {"native_scc_capabilities", json::object{
                {"regular", true}, {"identity_gauge", true},
                {"identity_v", true}, {"no_pseudo", true}}}}},
      {"scc", json::object{
           {"components", nested(json::array{0})},
           {"structural_edges", json::array{}},
           {"condensation_edges", json::array{}},
           {"topological_order", json::array{0}}, {"coupling_depth", 0}}},
      {"problem", prepared_problem(domain)}});
  if (response.at("status") != "ok")
    throw std::runtime_error("prepare: " + json::serialize(response));
  return std::string(response.at("chart").as_string());
}

json::object regular_run(bool sourced) {
  constexpr std::uint32_t nmax = 4;
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    shifts.push_back(json::string(std::to_string(n)));
    schedule.push_back(json::array{json::object{
        {"case", n == 0 ? "R" : "T"},
        {"da", json::string(std::to_string(n))}, {"db", "0"}}});
  }
  json::value source = nullptr;
  if (sourced) {
    json::array frames;
    json::array validity;
    json::array present;
    for (std::uint32_t n = 0; n <= nmax; ++n) {
      frames.push_back(n == 1 ? "1" : "0");
      validity.push_back(0);
      present.push_back(true);
    }
    source = json::object{{"frames", std::move(frames)},
                          {"validity", std::move(validity)},
                          {"present", std::move(present)}};
  }
  return json::object{
      {"nmax", nmax}, {"p", 0}, {"has_initial", true},
      {"adaptive_probe", false}, {"a_target", "0"},
      {"b_target", "0"}, {"a_shift_min", 0},
      {"a_shifts", std::move(shifts)}, {"schedule", std::move(schedule)},
      {"initial", json::array{"1"}}, {"initial_validity", json::array{0}},
      {"source", std::move(source)}, {"return_u", false}};
}

json::object metadata(const std::string& center,
                      const std::string& checkpoint) {
  return json::object{
      {"chart", json::object{{"center_exact", center},
                              {"scale_exact", "1"}, {"radius", "2"},
                              {"infinite_radius", false}}},
      {"tag", json::object{
           {"a", json::object{{"domain", "rational"},
                               {"canonical", "0"}}},
           {"b", json::object{{"domain", "rational"},
                               {"canonical", "0"}}}}},
      {"prescriptions", json::array{}},
      {"checkpoint_identity", checkpoint}};
}

json::object solve_local(const std::string& session,
                         const std::string& chart,
                         const std::string& center,
                         const std::string& checkpoint,
                         bool sourced = false) {
  const auto response = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart}, {"run", regular_run(sourced)},
      {"metadata", metadata(center, checkpoint)}});
  if (response.at("status") != "ok")
    throw std::runtime_error("solve: " + json::serialize(response));
  return response;
}

json::object topology() {
  return json::object{{"singular_points", json::array{}},
                      {"boundary_points", json::array{}},
                      {"complex_projections", json::array{}},
                      {"branch_sheets", json::array{}}};
}

json::object arm(const std::string& endpoint,
                 const std::string& anchor,
                 const std::string& endpoint_chart) {
  return json::object{{"from_exact", "0"}, {"to_exact", endpoint},
                      {"charts", json::array{anchor, endpoint_chart}},
                      {"topology", topology()}};
}

json::object exact_lattice() {
  return json::object{
      {"schema", "diffexp3-exact-evaluated-epsilon-lattice-v1"},
      {"identity", "matched-tail-exp-lattice-v1"},
      {"evaluated_basis", nested(json::array{
           json::object{{"min", 0}, {"max", 0},
                        {"coefficients", json::array{"1393/1944"}}}})}};
}

json::object advance(const std::string& domain,
                     const std::string& session,
                     const std::string& plan,
                     const std::string& basis,
                     const std::string& incoming,
                     const std::string& checkpoint) {
  json::object match{
      {"schema", 2}, {"op", "tile.match_advance"},
      {"session", session}, {"tile_plan", plan}, {"arm", "upper"},
      {"match", 0}, {"basis", json::array{basis}},
      {"incoming", incoming},
      {"epsilon", json::object{{"min", 0}, {"max", 0},
                                {"required_complete_max", 0}}},
      {"checkpoint_identity", checkpoint}};
  if (domain == "acb") {
    match["exact_lattice"] = exact_lattice();
    match["refinement"] = json::object{
        {"relative_tolerance", "1e-40"}, {"max_steps", 2}};
  }
  const auto response = request(std::move(match));
  if (response.at("status") != "ok")
    throw std::runtime_error("match: " + json::serialize(response));
  return response;
}

json::object materialize(const std::string& session,
                         const std::string& match,
                         const std::string& checkpoint) {
  const auto response = request(json::object{
      {"schema", 2}, {"op", "match.materialize_local"},
      {"session", session}, {"match", match},
      {"checkpoint_identity", checkpoint}});
  if (response.at("status") != "ok")
    throw std::runtime_error("materialize: " + json::serialize(response));
  return response;
}

json::object evaluate(const std::string& session,
                      const std::string& local) {
  return request(json::object{
      {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
      {"local", local}, {"point", json::object{{"exact", "1/4"}}},
      {"options", json::object{{"tail_estimate", false},
                                {"certified_tail_radius_exact", "1"}}}});
}

json::object integrate(const std::string& session,
                       const std::string& plan,
                       const std::string& local,
                       const std::string& source_checkpoint,
                       const std::string& result_checkpoint) {
  return request(json::object{
      {"schema", 2}, {"op", "integration.line"},
      {"session", session}, {"tile_plan", plan}, {"local", local},
      {"arm", "upper"}, {"tile", 1},
      {"epsilon", json::object{{"min", 0}, {"max", 0}}},
      {"source_checkpoint_identity", source_checkpoint},
      {"tile_plan_checkpoint_identity", "matched-tail-plan-v1"},
      {"checkpoint_identity", result_checkpoint}, {"certify_tail", true}});
}

const json::object& retained_local_record(const json::object& payload,
                                          const std::string& handle) {
  for (const auto& raw : payload.at("retained_locals").as_array()) {
    const auto& item = raw.as_object();
    if (std::string(item.at("handle").as_string()) == handle)
      return item;
  }
  throw std::runtime_error("checkpoint omitted retained local " + handle);
}

json::object& retained_local_record(json::object& payload,
                                    const std::string& handle) {
  for (auto& raw : payload.at("retained_locals").as_array()) {
    auto& item = raw.as_object();
    if (std::string(item.at("handle").as_string()) == handle)
      return item;
  }
  throw std::runtime_error("checkpoint omitted retained local " + handle);
}

json::object checkpoint_payload(const std::string& path) {
  return json::parse(diffexp::kernel::checkpoint::read(path).payload_json)
      .as_object();
}

void write_corrupt_tail_checkpoint(const std::string& source_path,
                                   const std::string& corrupt_path,
                                   const std::string& local_handle) {
  const auto container = diffexp::kernel::checkpoint::read(source_path);
  auto header = json::parse(container.header_json).as_object();
  auto payload = json::parse(container.payload_json).as_object();
  auto& tail = retained_local_record(payload, local_handle)
      .at("tail_model_restore").as_object();
  if (tail.if_contains("model") == nullptr)
    throw std::runtime_error(
        "matched tail checkpoint has no corruptible model: " +
        json::serialize(tail));
  auto& model = tail.at("model").as_object();
  auto& n_coefficients = model.at("n_coefficients").as_array();
  auto& n_norms = model.at("n_row_sum_upper_exact").as_array();
  n_coefficients.at(1).as_array().front() =
      n_coefficients.front().as_array().front();
  n_norms.at(1) = n_norms.front();
  diffexp::kernel::checkpoint::write_atomic(
      corrupt_path, json::serialize(header), json::serialize(payload));
}

bool certified_domain(const std::string& domain, bool incompatible_case) {
  json::object create{{"schema", 2}, {"op", "session.create"},
                      {"domain", domain}, {"output_digits", 40},
                      {"local_capacity", 12}, {"match_capacity", 6},
                      {"tile_plan_capacity", 2},
                      {"line_result_capacity", 4}};
  if (domain == "acb") create["precision_bits"] = 256;
  const auto created = request(std::move(create));
  const auto session = std::string(created.at("session").as_string());
  const auto anchor = prepare_chart(session, domain, domain + "-anchor", "0");
  const auto receiver = prepare_chart(
      session, domain, domain + "-receiver", "2/3");
  const auto lower = prepare_chart(
      session, domain, domain + "-lower", "-2/3");
  const auto incoming = solve_local(
      session, anchor, "0", domain + "-incoming-v1");
  const auto basis = solve_local(
      session, receiver, "2/3", domain + "-basis-v1");
  const auto basis_handle = std::string(basis.at("local").as_string());
  const auto planned = request(json::object{
      {"schema", 2}, {"op", "tile.plan"}, {"session", session},
      {"checkpoint_identity", "matched-tail-plan-v1"},
      {"division_order", 3}, {"lower", arm("-2/3", anchor, lower)},
      {"upper", arm("2/3", anchor, receiver)}});
  const auto plan = std::string(planned.at("tile_plan").as_string());
  const auto hop = advance(
      domain, session, plan,
      basis_handle,
      std::string(incoming.at("local").as_string()),
      domain + "-matched-hop-v1");
  const auto matched = materialize(
      session, std::string(hop.at("match").as_string()),
      domain + "-matched-local-v1");
  const auto matched_handle = std::string(matched.at("local").as_string());
  const auto evaluated = evaluate(session, matched_handle);
  const auto direct_evaluated = evaluate(session, basis_handle);
  const auto line = integrate(
      session, plan, matched_handle, domain + "-matched-local-v1",
      domain + "-matched-line-v1");
  const auto direct_line = integrate(
      session, plan, basis_handle, domain + "-basis-v1",
      domain + "-direct-line-v1");

  const auto& tail = matched.at("tail_majorant").as_object();
  const auto& evaluation_tail = evaluated.at("tail_certificate").as_object();
  const auto& line_diagnostics = line.at("diagnostics").as_object();
  bool ok = tail.at("status") == "certified" &&
      tail.at("attached") == true &&
      tail.at("checkpoint_serialized") == true &&
      std::string(tail.at("operator_identity").as_string()) ==
          domain + "-receiver-operator-v1" &&
      std::string(tail.at("local_checkpoint_identity").as_string()) ==
          domain + "-matched-local-v1" &&
      std::string(tail.at("provenance").as_string()).find(
          "match materialization") != std::string::npos &&
      evaluation_tail.at("status") == "certified" &&
      evaluated.at("value").as_object().at("error").as_object().at(
          "guarantee") == "certified" &&
      line.at("scope") == "full_local_with_certified_tail" &&
      line.at("error").as_object().at("guarantee") == "certified" &&
      line_diagnostics.at("tail_certificate_status") == "certified" &&
      basis.at("tail_majorant").as_object().at("status") == "certified" &&
      basis.at("tail_majorant").as_object().at("checkpoint_serialized") ==
          true &&
      direct_evaluated.at("tail_certificate").as_object().at("status") ==
          "certified" &&
      direct_line.at("scope") == "full_local_with_certified_tail" &&
      direct_line.at("error").as_object().at("guarantee") == "certified";

  json::object incompatible;
  json::object incompatible_evaluation;
  json::object incompatible_line;
  std::string incompatible_handle;
  if (incompatible_case) {
    const auto sourced_basis = solve_local(
        session, receiver, "2/3", "sourced-basis-v1", true);
    const auto sourced_hop = advance(
        domain, session, plan,
        std::string(sourced_basis.at("local").as_string()),
        std::string(incoming.at("local").as_string()),
        "sourced-hop-v1");
    incompatible = materialize(
        session, std::string(sourced_hop.at("match").as_string()),
        "sourced-matched-local-v1");
    incompatible_handle = std::string(incompatible.at("local").as_string());
    incompatible_evaluation = evaluate(session, incompatible_handle);
    incompatible_line = integrate(
        session, plan, incompatible_handle, "sourced-matched-local-v1",
        "sourced-matched-line-v1");
    ok = ok &&
        incompatible.at("tail_majorant").as_object().at("status") ==
            "unsupported" &&
        incompatible.at("tail_majorant").as_object().at("attached") ==
            false &&
        incompatible_evaluation.at("tail_certificate").as_object().at(
            "status") == "unsupported" &&
        incompatible_evaluation.at("value").as_object().if_contains(
            "error") == nullptr &&
        incompatible_line.at("scope") == "stored_truncation" &&
        incompatible_line.at("error").as_object().at("guarantee") ==
            "none" &&
        incompatible_line.at("diagnostics").as_object().at(
            "tail_certificate_status") == "unsupported";
  }

  const auto release_line = [&](const json::object& result) {
    if (result.if_contains("line") == nullptr) return;
    (void)request(json::object{
        {"schema", 2}, {"op", "integration.release"},
        {"session", session}, {"line", result.at("line")}});
  };
  release_line(line);
  release_line(direct_line);
  release_line(incompatible_line);
  const auto base = std::filesystem::temp_directory_path() /
      ("diffexp3-matched-tail-checkpoint-" + domain + "-" +
       std::to_string(::getpid()));
  const auto path = base.string() + ".de2cp";
  const auto corrupt_path = base.string() + "-corrupt.de2cp";
  const auto resaved_path = base.string() + "-resaved.de2cp";
  const auto checkpoint_identity = domain + "-tail-checkpoint-v1";
  const auto saved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
      {"path", path}, {"checkpoint_identity", checkpoint_identity}});
  const auto saved_payload = checkpoint_payload(path);
  const auto& saved_direct_tail = retained_local_record(
      saved_payload, basis_handle).at("tail_model_restore").as_object();
  const auto& saved_matched_tail = retained_local_record(
      saved_payload, matched_handle).at("tail_model_restore").as_object();
  if (saved_matched_tail.if_contains("model") == nullptr)
    throw std::runtime_error(
        domain + " matched tail was not retained: " +
        json::serialize(matched.at("tail_majorant")));
  write_corrupt_tail_checkpoint(path, corrupt_path, matched_handle);
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  const auto corruption = request(json::object{
      {"schema", 2}, {"op", "checkpoint.restore"},
      {"path", corrupt_path}, {"expected_identity", checkpoint_identity}});
  const auto restored = request(json::object{
      {"schema", 2}, {"op", "checkpoint.restore"}, {"path", path},
      {"expected_identity", checkpoint_identity}});
  if (restored.at("status") != "ok")
    throw std::runtime_error(
        domain + " matched-tail checkpoint restore: " +
        json::serialize(restored));
  const auto restored_session =
      std::string(restored.at("session").as_string());
  const auto resaved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"},
      {"session", restored_session}, {"path", resaved_path},
      {"checkpoint_identity", domain + "-tail-resaved-v1"}});
  const auto resaved_payload = checkpoint_payload(resaved_path);
  const auto restored_direct_evaluated = evaluate(
      restored_session, basis_handle);
  const auto restored_matched_evaluated = evaluate(
      restored_session, matched_handle);
  const auto restored_direct_line = integrate(
      restored_session, plan, basis_handle, domain + "-basis-v1",
      domain + "-restored-direct-line-v1");
  const auto restored_matched_line = integrate(
      restored_session, plan, matched_handle,
      domain + "-matched-local-v1",
      domain + "-restored-matched-line-v1");
  json::object restored_incompatible;
  if (incompatible_case)
    restored_incompatible = evaluate(
        restored_session, incompatible_handle);

  ok = ok && saved.at("status") == "ok" &&
      saved_direct_tail.at("serialized") == true &&
      saved_direct_tail.at("status") == "certified" &&
      saved_direct_tail.at("attached_before_save") == true &&
      saved_direct_tail.if_contains("model") != nullptr &&
      saved_matched_tail.at("serialized") == true &&
      saved_matched_tail.at("status") == "certified" &&
      saved_matched_tail.at("attached_before_save") == true &&
      saved_matched_tail.if_contains("model") != nullptr &&
      corruption.at("status") == "error" &&
      std::string(corruption.at("detail").as_string()).find(
          "q/N payload") != std::string::npos &&
      restored.at("status") == "ok" &&
      resaved.at("status") == "ok" &&
      retained_local_record(saved_payload, basis_handle) ==
          retained_local_record(resaved_payload, basis_handle) &&
      retained_local_record(saved_payload, matched_handle) ==
          retained_local_record(resaved_payload, matched_handle) &&
      restored_direct_evaluated.at("tail_certificate").as_object().at(
          "status") == "certified" &&
      restored_matched_evaluated.at("tail_certificate").as_object().at(
          "status") == "certified" &&
      restored_direct_line.at("scope") ==
          "full_local_with_certified_tail" &&
      restored_direct_line.at("error").as_object().at("guarantee") ==
          "certified" &&
      restored_matched_line.at("scope") ==
          "full_local_with_certified_tail" &&
      restored_matched_line.at("error").as_object().at("guarantee") ==
          "certified";
  if (incompatible_case) {
    const auto& marker = retained_local_record(
        saved_payload, incompatible_handle).at(
            "tail_model_restore").as_object();
    ok = ok && marker.at("serialized") == false &&
        marker.at("status") == "unsupported" &&
        marker.if_contains("model") == nullptr &&
        restored_incompatible.at("tail_certificate").as_object().at(
            "status") == "unsupported";
  }

  if (!ok) {
    std::cerr << domain << " matched: " << json::serialize(matched) << '\n'
              << domain << " evaluate: " << json::serialize(evaluated)
              << '\n' << domain << " line: " << json::serialize(line)
              << '\n';
    if (incompatible_case)
      std::cerr << "incompatible: " << json::serialize(incompatible) << '\n'
                << "incompatible evaluate: "
                << json::serialize(incompatible_evaluation) << '\n'
                << "incompatible line: "
                << json::serialize(incompatible_line) << '\n';
    std::cerr << domain << " saved: " << json::serialize(saved) << '\n'
              << domain << " corruption: " << json::serialize(corruption)
              << '\n' << domain << " restored: "
              << json::serialize(restored) << '\n'
              << domain << " resaved: " << json::serialize(resaved) << '\n'
              << domain << " restored direct evaluate: "
              << json::serialize(restored_direct_evaluated) << '\n'
              << domain << " restored matched evaluate: "
              << json::serialize(restored_matched_evaluated) << '\n'
              << domain << " restored direct line: "
              << json::serialize(restored_direct_line) << '\n'
              << domain << " restored matched line: "
              << json::serialize(restored_matched_line) << '\n';
  }
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", restored_session}});
  std::error_code ignored;
  std::filesystem::remove(path, ignored);
  std::filesystem::remove(corrupt_path, ignored);
  std::filesystem::remove(resaved_path, ignored);
  return ok;
}

}  // namespace

int main() {
  const bool rational = certified_domain("rational", true);
  const bool acb = certified_domain("acb", false);
  const bool ok = rational && acb;
  std::cout << (ok ? "PASS" : "FAIL")
            << ": certified regular tails through retained match materialization\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
