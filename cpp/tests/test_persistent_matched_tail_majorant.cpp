#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
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
      {"schema", "diffexp2-exact-evaluated-epsilon-lattice-v1"},
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
  const auto planned = request(json::object{
      {"schema", 2}, {"op", "tile.plan"}, {"session", session},
      {"checkpoint_identity", "matched-tail-plan-v1"},
      {"division_order", 3}, {"lower", arm("-2/3", anchor, lower)},
      {"upper", arm("2/3", anchor, receiver)}});
  const auto plan = std::string(planned.at("tile_plan").as_string());
  const auto hop = advance(
      domain, session, plan,
      std::string(basis.at("local").as_string()),
      std::string(incoming.at("local").as_string()),
      domain + "-matched-hop-v1");
  const auto matched = materialize(
      session, std::string(hop.at("match").as_string()),
      domain + "-matched-local-v1");
  const auto matched_handle = std::string(matched.at("local").as_string());
  const auto evaluated = evaluate(session, matched_handle);
  const auto line = integrate(
      session, plan, matched_handle, domain + "-matched-local-v1",
      domain + "-matched-line-v1");

  const auto& tail = matched.at("tail_majorant").as_object();
  const auto& evaluation_tail = evaluated.at("tail_certificate").as_object();
  const auto& line_diagnostics = line.at("diagnostics").as_object();
  bool ok = tail.at("status") == "certified" &&
      tail.at("attached") == true &&
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
      line_diagnostics.at("tail_certificate_status") == "certified";

  json::object incompatible;
  json::object incompatible_evaluation;
  json::object incompatible_line;
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
    const auto incompatible_handle =
        std::string(incompatible.at("local").as_string());
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
  }
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
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
