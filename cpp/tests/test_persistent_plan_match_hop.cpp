#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdlib>
#include <future>
#include <iostream>
#include <stdexcept>
#include <string>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
      .as_object();
}

json::array prescriptions() {
  return json::array{json::object{
      {"factor_exact", "physical-f"}, {"sign", -1},
      {"multiplicity", 1}, {"leading_coefficient_sign", 1}}};
}

std::string prepare_chart(const std::string& session,
                          const std::string& domain,
                          const std::string& name,
                          const std::string& center) {
  json::array d_lag;
  d_lag.push_back(json::object{{"s", 0}, {"v", "1"}});
  json::array d_lags;
  d_lags.push_back(std::move(d_lag));
  json::array block;
  block.emplace_back(0);
  json::array blocks;
  blocks.push_back(block);
  json::array principal_row;
  principal_row.push_back(
      json::object{{"exact", "0"}, {"proven_zero", true}});
  json::array principal_matrix;
  principal_matrix.push_back(std::move(principal_row));
  json::array components;
  components.push_back(std::move(block));
  json::object problem{
      {"domain", domain}, {"d", 1}, {"fb", 0}, {"w", 3},
      {"d_lags", std::move(d_lags)},
      {"denominators", json::array{}},
      {"nhat_lags", json::array{json::object{
           {"poly", json::array{}}, {"rat", json::array{}},
           {"val", json::array{nullptr}}}}},
      {"d0_inverse", "1"}, {"blocks", std::move(blocks)},
      {"assembly", json::object{
           {"identity", true}, {"poly", json::array{}},
           {"rat", json::array{}}, {"val", json::array{0}}}},
      {"chop_digits", 0}};
  if (domain == "acb") problem["precision_bits"] = 256;
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
      {"problem", std::move(problem)}});
  if (response.at("status") != "ok")
    throw std::runtime_error("prepare: " + json::serialize(response));
  return std::string(response.at("chart").as_string());
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
  const auto response = request(json::object{
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
           {"factor_exact", "physical-f"}, {"sign", -1}}}}};
}

json::object arm(const std::string& endpoint,
                 const std::string& anchor,
                 const std::string& endpoint_chart) {
  return json::object{
      {"from_exact", "0"}, {"to_exact", endpoint},
      {"charts", json::array{anchor, endpoint_chart}},
      {"topology", topology()}};
}

json::object exact_lattice() {
  json::array row;
  row.push_back(json::object{
      {"min", 0}, {"max", 2},
      {"coefficients", json::array{"1", "0", "0"}}});
  json::array matrix;
  matrix.push_back(std::move(row));
  return json::object{
      {"schema", "diffexp2-exact-evaluated-epsilon-lattice-v1"},
      {"identity", "planned-hop-unit-lattice"},
      {"evaluated_basis", std::move(matrix)}};
}

json::object hop_request(const std::string& domain,
                         const std::string& session,
                         const std::string& plan,
                         const std::string& selected_arm,
                         const std::string& basis,
                         const std::string& incoming,
                         const std::string& checkpoint) {
  json::object value{
      {"schema", 2}, {"op", "tile.match_advance"},
      {"session", session}, {"tile_plan", plan},
      {"arm", selected_arm}, {"match", 0},
      {"basis", json::array{basis}}, {"incoming", incoming},
      {"epsilon", json::object{{"min", 0}, {"max", 2},
                                {"required_complete_max", 1}}},
      {"checkpoint_identity", checkpoint}};
  if (domain == "acb") {
    value["exact_lattice"] = exact_lattice();
    value["refinement"] = json::object{
        {"relative_tolerance", "1e-40"}, {"max_steps", 1}};
  }
  return value;
}

bool run_domain(const std::string& domain) {
  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", domain},
      {"precision_bits", 256}, {"output_digits", 30},
      {"chart_capacity", 3}, {"local_capacity", 3},
      {"match_capacity", 2}, {"tile_plan_capacity", 1}});
  const auto session = std::string(created.at("session").as_string());
  const auto anchor = prepare_chart(session, domain, domain + "-anchor", "0");
  const auto lower_chart = prepare_chart(
      session, domain, domain + "-lower", "-2/3");
  const auto upper_chart = prepare_chart(
      session, domain, domain + "-upper", "2/3");
  const auto incoming = solve_local(
      session, anchor, "0", domain + "-incoming", "2");
  const auto lower_basis = solve_local(
      session, lower_chart, "-2/3", domain + "-lower-basis", "1");
  const auto upper_basis = solve_local(
      session, upper_chart, "2/3", domain + "-upper-basis", "1");

  const auto planned = request(json::object{
      {"schema", 2}, {"op", "tile.plan"}, {"session", session},
      {"checkpoint_identity", domain + "-tile-plan"},
      {"division_order", 3},
      {"lower", arm("-2/3", anchor, lower_chart)},
      {"upper", arm("2/3", anchor, upper_chart)}});
  const auto plan = std::string(planned.at("tile_plan").as_string());

  auto submit = [&](const std::string& selected_arm,
                    const std::string& basis) {
    return request(hop_request(
        domain, session, plan, selected_arm, basis, incoming,
        domain + "-" + selected_arm + "-hop"));
  };
  auto lower_future = std::async(
      std::launch::async, submit, "lower", lower_basis);
  auto upper_future = std::async(
      std::launch::async, submit, "upper", upper_basis);
  const auto lower = lower_future.get();
  const auto upper = upper_future.get();
  const auto lower_match = std::string(lower.at("match").as_string());
  const auto upper_match = std::string(upper.at("match").as_string());
  const auto plan_stats = request(json::object{
      {"schema", 2}, {"op", "tile.stats"}, {"session", session},
      {"tile_plan", plan}});

  for (const auto& local : {incoming, lower_basis, upper_basis})
    (void)request(json::object{
        {"schema", 2}, {"op", "local.release"}, {"session", session},
        {"local", local}});
  (void)request(json::object{
      {"schema", 2}, {"op", "tile.release"}, {"session", session},
      {"tile_plan", plan}});
  const auto retained_lower = request(json::object{
      {"schema", 2}, {"op", "match.stats"}, {"session", session},
      {"match", lower_match}});
  const auto checkpoint_rejection = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
      {"path", "/tmp/diffexp2-plan-match-hop.chk"},
      {"checkpoint_identity", domain + "-must-retain-hop"}});
  const auto stats = request(json::object{
      {"schema", 2}, {"op", "session.stats"}, {"session", session}});

  const auto& lower_hop = lower.at("planned_hop").as_object();
  const auto& upper_hop = upper.at("planned_hop").as_object();
  const auto& lower_geometry = lower_hop.at("geometry").as_object();
  const auto& upper_geometry = upper_hop.at("geometry").as_object();
  const auto& lower_advance = lower_hop.at("advance").as_object();
  const auto& lower_producing = lower_hop.at("producing").as_object();
  const auto& lower_receiving = lower_hop.at("receiving").as_object();
  const bool domain_match_ok = domain == "rational"
      ? lower.at("residual").as_object().at("status") == "exact-zero"
      : lower.at("residual").as_object().at("verdict") == "pass";
  const bool ok =
      created.at("planned_match_hop_capability") ==
          "retained-exact-plan-driven-local-match-hop-v1" &&
      lower.at("status") == "ok" && upper.at("status") == "ok" &&
      lower.at("plan_driven") == true &&
      lower.at("planned_hop_capability") ==
          "retained-exact-plan-driven-local-match-hop-v1" &&
      lower_geometry.at("physical_exact") == "-1/3" &&
      lower_geometry.at("producing_local_exact") == "-1/3" &&
      lower_geometry.at("receiving_local_exact") == "1/3" &&
      upper_geometry.at("physical_exact") == "1/3" &&
      upper_geometry.at("producing_local_exact") == "1/3" &&
      upper_geometry.at("receiving_local_exact") == "-1/3" &&
      lower_geometry.at("branch_sheets").as_array().front().as_object()
          .at("sign") == -1 &&
      std::string(lower_hop.at("tile_plan_checkpoint_identity").as_string()) ==
          domain + "-tile-plan" &&
      lower_producing.at("effective_rim") == -1 &&
      lower_receiving.at("effective_rim") == -1 &&
      std::string(lower_producing.at("incoming").as_object()
                      .at("checkpoint_identity").as_string()) ==
          domain + "-incoming" &&
      std::string(lower_receiving.at("basis").as_array().front().as_object()
                      .at("checkpoint_identity").as_string()) ==
          domain + "-lower-basis" &&
      lower_advance.at("scope") == "single-match-handoff" &&
      lower_advance.at("state") ==
          "retained-receiving-basis-weights" &&
      lower_advance.at("whole_arm_complete") == false &&
      lower.at("strong_ownership").as_object().at("tile_plan") == true &&
      plan_stats.at("lower_match_advances") == 1 &&
      plan_stats.at("upper_match_advances") == 1 && domain_match_ok &&
      retained_lower.at("status") == "ok" &&
      retained_lower.at("planned_hop_provenance_identity") ==
          lower.at("planned_hop_provenance_identity") &&
      checkpoint_rejection.at("status") == "error" &&
      std::string(checkpoint_rejection.at("detail").as_string()).find(
          "plan-driven") != std::string::npos &&
      stats.at("locals") == 0 && stats.at("tile_plans") == 0 &&
      stats.at("matches") == 2 && stats.at("local_matches") == 2 &&
      stats.at("pending_matches") == 0;

  if (!ok)
    std::cerr << domain << " lower: " << json::serialize(lower) << '\n'
              << domain << " upper: " << json::serialize(upper) << '\n'
              << domain << " plan stats: " << json::serialize(plan_stats)
              << '\n' << domain << " retained: "
              << json::serialize(retained_lower) << '\n'
              << domain << " checkpoint: "
              << json::serialize(checkpoint_rejection) << '\n'
              << domain << " stats: " << json::serialize(stats) << '\n';

  for (const auto& match : {lower_match, upper_match})
    (void)request(json::object{
        {"schema", 2}, {"op", "match.release"}, {"session", session},
        {"match", match}});
  (void)request(json::object{
      {"schema", 2}, {"op", "session.close"}, {"session", session}});
  return ok;
}

}  // namespace

int main() {
  const bool rational = run_domain("rational");
  const bool acb = run_domain("acb");
  const bool ok = rational && acb;
  std::cout << (ok ? "PASS" : "FAIL")
            << ": retained exact plan-driven rational/Acb match hops\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
