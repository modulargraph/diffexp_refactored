#include "diffexp2/checkpoint.hpp"
#include "diffexp2/json_codec.hpp"
#include "diffexp2/local_algebra.hpp"

#include <boost/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <future>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <unistd.h>

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

json::object epsilon_one() {
  return json::object{{"zero", false}, {"valuation", 0},
                      {"numerator", json::array{"1"}},
                      {"denominator", json::array{"1"}}};
}

json::object physical_zero_ode(const std::string& operator_identity,
                               const std::string& payload_identity) {
  json::array q;
  q.push_back(epsilon_one());
  json::array c;
  c.push_back(json::array{});
  return json::object{
      {"schema", "diffexp2-physical-cleared-ode-v1"},
      {"basis", "physical-original-master"},
      {"theta_coordinate", "local-t"},
      {"owner_signature_identity", operator_identity},
      {"payload_identity", payload_identity},
      {"q", std::move(q)},
      {"c", std::move(c)}};
}

std::string prepare_chart(const std::string& session,
                          const std::string& domain,
                          const std::string& name,
                          const std::string& center) {
  const auto operator_identity = "de2-operator-" + name;
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
      {"physical_ode", physical_zero_ode(
           operator_identity, "de2-physical-ode-" + name)},
      {"chop_digits", 0}};
  if (domain == "acb") problem["precision_bits"] = 256;
  const auto response = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", name}, {"identity", operator_identity},
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
                        const std::string& value,
                        const std::string& epsilon_value = "0") {
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
           {"initial", json::array{value, epsilon_value, "0"}},
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

json::object residual_request(const std::string& session,
                              const std::string& local,
                              const json::object& binding,
                              const std::string& point) {
  return json::object{
      {"schema", 2}, {"op", "local.certify_residual"},
      {"session", session}, {"local", local},
      {"point", json::object{{"exact", point}}},
      {"options", json::object{{"tail_estimate", false}}},
      {"relative_tolerance", "1e-25"},
      {"scope", "stored_truncation"}, {"include_residual", false},
      {"operator_identity", binding.at("operator_identity")},
      {"source_identity", binding.at("source_identity")},
      {"checkpoint_identity", binding.at("local_checkpoint_identity")},
      {"analytic_metadata", binding.at("analytic_metadata")},
      {"owner_signature_identity", binding.at("owner_signature_identity")},
      {"physical_payload_identity", binding.at("physical_payload_identity")},
      {"provenance_identity", binding.at("provenance_identity")}};
}

double value_midpoint(const json::object& evaluation,
                      std::int32_t epsilon_power) {
  const auto& value = evaluation.at("value").as_object();
  const auto minimum = static_cast<std::int32_t>(
      value.at("min").as_int64());
  const auto& coefficient = value.at("coefficients").as_array().at(
      static_cast<std::size_t>(epsilon_power - minimum));
  return coefficient.is_string()
      ? std::stod(std::string(coefficient.as_string()))
      : std::stod(std::string(
            coefficient.as_array().front().as_string()));
}

bool contains_string(const json::array& values, const std::string& expected) {
  return std::any_of(values.begin(), values.end(), [&](const auto& value) {
    return value.is_string() &&
        std::string(value.as_string()) == expected;
  });
}

bool finite_laurent_pole_materialization() {
  diffexp2::LocalSolution<diffexp2::Rational> column;
  column.chart.center_exact = "0";
  column.chart.scale_exact = "1";
  column.chart.radius = diffexp2::ComplexBall(2);
  column.epsilon = {0, 2};
  column.taylor_complete_max = 0;
  column.dimension = 1;
  column.checkpoint_identity = "pole-basis";
  diffexp2::LocalSector<diffexp2::Rational> sector;
  sector.a = diffexp2::ExactScalarDescriptor::rational("0");
  sector.b = diffexp2::ExactScalarDescriptor::rational("0");
  sector.coefficients = {
      diffexp2::Rational(0), diffexp2::Rational(1),
      diffexp2::Rational(0)};
  column.sectors.push_back(std::move(sector));
  const diffexp2::EpsilonFrame<diffexp2::Rational> pole_weight(
      {-1, 1}, {diffexp2::Rational(2), diffexp2::Rational(0),
                diffexp2::Rational(0)});
  const auto materialized = diffexp2::materialize_local_basis_weights(
      std::vector<const diffexp2::LocalSolution<diffexp2::Rational>*>{
          &column},
      diffexp2::FiniteLaurentVector<diffexp2::Rational>{pole_weight},
      "pole-materialized");
  const auto& coefficients = materialized.sectors.front().coefficients;
  return materialized.epsilon.min_power == -1 &&
      materialized.epsilon.complete_max == 1 &&
      coefficients.size() == 3 && coefficients[0].is_zero() &&
      coefficients[1] == diffexp2::Rational(2) &&
      coefficients[2].is_zero();
}

bool run_domain(const std::string& domain) {
  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", domain},
      {"precision_bits", 256}, {"output_digits", 30},
      {"chart_capacity", 6}, {"local_capacity", 8},
      {"match_capacity", 4}, {"tile_plan_capacity", 2}});
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
  if (lower.at("status") != "ok" || upper.at("status") != "ok")
    throw std::runtime_error(
        "match: " + json::serialize(lower) + " / " +
        json::serialize(upper));
  const auto lower_match = std::string(lower.at("match").as_string());
  const auto upper_match = std::string(upper.at("match").as_string());

  auto materialize = [&](const std::string& match,
                         const std::string& checkpoint) {
    return request(json::object{
        {"schema", 2}, {"op", "match.materialize_local"},
        {"session", session}, {"match", match},
        {"checkpoint_identity", checkpoint}});
  };
  auto lower_materialized_future = std::async(
      std::launch::async, materialize, lower_match,
      domain + "-lower-receiving-local");
  auto upper_materialized_future = std::async(
      std::launch::async, materialize, upper_match,
      domain + "-upper-receiving-local");
  const auto lower_materialized = lower_materialized_future.get();
  const auto upper_materialized = upper_materialized_future.get();
  if (lower_materialized.at("status") != "ok" ||
      upper_materialized.at("status") != "ok")
    throw std::runtime_error(
        "materialize: " + json::serialize(lower_materialized) + " / " +
        json::serialize(upper_materialized));
  const auto lower_local =
      std::string(lower_materialized.at("local").as_string());
  const auto upper_local =
      std::string(upper_materialized.at("local").as_string());
  const auto plan_stats = request(json::object{
      {"schema", 2}, {"op", "tile.stats"}, {"session", session},
      {"tile_plan", plan}});

  const auto retained_lower = request(json::object{
      {"schema", 2}, {"op", "match.stats"}, {"session", session},
      {"match", lower_match}});
  for (const auto& local : {incoming, lower_basis, upper_basis})
    (void)request(json::object{
        {"schema", 2}, {"op", "local.release"}, {"session", session},
        {"local", local}});
  (void)request(json::object{
      {"schema", 2}, {"op", "tile.release"}, {"session", session},
      {"tile_plan", plan}});
  for (const auto& match : {lower_match, upper_match})
    (void)request(json::object{
        {"schema", 2}, {"op", "match.release"}, {"session", session},
        {"match", match}});
  for (const auto& chart : {anchor, lower_chart, upper_chart})
    (void)request(json::object{
        {"schema", 2}, {"op", "chart.release"}, {"session", session},
        {"chart", chart}});

  const auto lower_evaluation = request(json::object{
      {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
      {"local", lower_local},
      {"point", json::object{{"exact", "1/5"}}},
      {"options", json::object{{"tail_estimate", false}}}});
  const auto upper_evaluation = request(json::object{
      {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
      {"local", upper_local},
      {"point", json::object{{"exact", "-1/5"}}},
      {"options", json::object{{"tail_estimate", false}}}});
  const auto& lower_residual_binding =
      lower_materialized.at("residual_binding").as_object()
          .at("binding").as_object();
  const auto lower_residual = request(residual_request(
      session, lower_local, lower_residual_binding, "1/5"));

  const auto base = std::filesystem::temp_directory_path() /
      ("diffexp2-plan-match-hop-" + domain + "-" +
       std::to_string(::getpid()));
  const auto path = base.string() + ".de2cp";
  const auto resaved_path = base.string() + "-resaved.de2cp";
  const auto corrupt_path = base.string() + "-corrupt.de2cp";
  const auto checkpoint_identity = domain + "-planned-hop-session-v2";
  const auto saved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
      {"path", path}, {"checkpoint_identity", checkpoint_identity}});
  const auto container = diffexp2::checkpoint::read(path);
  const auto saved_header = json::parse(container.header_json).as_object();
  const auto saved_payload = json::parse(container.payload_json).as_object();
  const auto& saved_visibility = saved_payload.at("session").as_object()
      .at("registry_visibility").as_object();
  const auto& planned_records =
      saved_payload.at("retained_planned_match_hops").as_array();
  const json::object* lower_planned = nullptr;
  for (const auto& record : planned_records) {
    const auto& object = record.as_object();
    if (std::string(object.at("handle").as_string()) == lower_match)
      lower_planned = &object;
  }
  if (lower_planned == nullptr)
    throw std::runtime_error("checkpoint omitted the lower planned match hop");
  const auto& embedded_native = lower_planned->at("native_match").as_object();
  const auto expected_native_schema = domain == "rational"
      ? "diffexp2-retained-exact-rational-match-v2"
      : "diffexp2-retained-acb-match-v2";

  json::object corruption{{"status", "ok"}};
  if (domain == "rational") {
    auto corrupt_header = saved_header;
    auto corrupt_payload = saved_payload;
    auto& corrupt_native = corrupt_payload
        .at("retained_planned_match_hops").as_array().front().as_object()
        .at("native_match").as_object();
    corrupt_native["schema"] = "diffexp2-retained-unknown-match-v9";
    corrupt_header.at("planned_match_identities").as_array().front()
        .as_object()["native_match_schema"] = corrupt_native.at("schema");
    diffexp2::checkpoint::write_atomic(
        corrupt_path, json::serialize(corrupt_header),
        json::serialize(corrupt_payload));
    corruption = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", corrupt_path}, {"expected_identity", checkpoint_identity}});
  }

  (void)request(json::object{
      {"schema", 2}, {"op", "session.close"}, {"session", session}});
  const auto restored = request(json::object{
      {"schema", 2}, {"op", "checkpoint.restore"}, {"path", path},
      {"expected_identity", checkpoint_identity}});
  const auto restored_session =
      std::string(restored.at("session").as_string());
  const auto restored_lower_evaluation = request(json::object{
      {"schema", 2}, {"op", "local.evaluate"},
      {"session", restored_session}, {"local", lower_local},
      {"point", json::object{{"exact", "1/5"}}},
      {"options", json::object{{"tail_estimate", false}}}});
  const auto restored_upper_evaluation = request(json::object{
      {"schema", 2}, {"op", "local.evaluate"},
      {"session", restored_session}, {"local", upper_local},
      {"point", json::object{{"exact", "-1/5"}}},
      {"options", json::object{{"tail_estimate", false}}}});
  const auto restored_lower_residual = request(residual_request(
      restored_session, lower_local, lower_residual_binding, "1/5"));
  const auto restored_local_stats = request(json::object{
      {"schema", 2}, {"op", "local.stats"},
      {"session", restored_session}, {"local", lower_local}});
  const auto hidden_match = request(json::object{
      {"schema", 2}, {"op", "match.stats"},
      {"session", restored_session}, {"match", lower_match}});
  const auto hidden_plan = request(json::object{
      {"schema", 2}, {"op", "tile.stats"},
      {"session", restored_session}, {"tile_plan", plan}});
  const auto hidden_source = request(json::object{
      {"schema", 2}, {"op", "local.stats"},
      {"session", restored_session}, {"local", incoming}});
  const auto resaved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"},
      {"session", restored_session}, {"path", resaved_path},
      {"checkpoint_identity", domain + "-planned-hop-resaved-v2"}});
  const auto resaved_payload = json::parse(
      diffexp2::checkpoint::read(resaved_path).payload_json).as_object();

  // The saved next-handle counters include every released strong owner.
  const auto next_anchor = prepare_chart(
      restored_session, domain, domain + "-next-anchor", "0");
  const auto next_lower_chart = prepare_chart(
      restored_session, domain, domain + "-next-lower", "-2/3");
  const auto next_upper_chart = prepare_chart(
      restored_session, domain, domain + "-next-upper", "2/3");
  const auto next_incoming = solve_local(
      restored_session, next_anchor, "0", domain + "-next-incoming", "2");
  const auto next_basis = solve_local(
      restored_session, next_lower_chart, "-2/3",
      domain + "-next-lower-basis", "1");
  const auto next_planned = request(json::object{
      {"schema", 2}, {"op", "tile.plan"}, {"session", restored_session},
      {"checkpoint_identity", domain + "-next-tile-plan"},
      {"division_order", 3},
      {"lower", arm("-2/3", next_anchor, next_lower_chart)},
      {"upper", arm("2/3", next_anchor, next_upper_chart)}});
  const auto next_plan =
      std::string(next_planned.at("tile_plan").as_string());
  const auto next_match_result = request(hop_request(
      domain, restored_session, next_plan, "lower", next_basis,
      next_incoming, domain + "-next-hop"));
  const auto stats = request(json::object{
      {"schema", 2}, {"op", "session.stats"},
      {"session", restored_session}});

  const auto& lower_hop = lower.at("planned_hop").as_object();
  const auto& upper_hop = upper.at("planned_hop").as_object();
  const auto& lower_geometry = lower_hop.at("geometry").as_object();
  const auto& upper_geometry = upper_hop.at("geometry").as_object();
  const auto& lower_advance = lower_hop.at("advance").as_object();
  const auto& lower_producing = lower_hop.at("producing").as_object();
  const auto& lower_receiving = lower_hop.at("receiving").as_object();
  const auto& materialized_derivation =
      lower_materialized.at("retained_derivation").as_object();
  const auto& materialized_output =
      materialized_derivation.at("output").as_object();
  const bool domain_match_ok = domain == "rational"
      ? lower.at("residual").as_object().at("status") == "exact-zero"
      : lower.at("residual").as_object().at("verdict") == "pass";
  const bool ok =
      created.at("planned_match_hop_capability") ==
          "retained-exact-plan-driven-local-match-hop-v1" &&
      created.at("planned_match_materialization_capability") ==
          "retained-native-plan-match-local-materialization-v1" &&
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
      retained_lower.at("materializations") == 1 &&
      lower_materialized.at("status") == "ok" &&
      upper_materialized.at("status") == "ok" &&
      lower_materialized.at("materialization_capability") ==
          "retained-native-plan-match-local-materialization-v1" &&
      lower_materialized.at("native_retained") == true &&
      lower_materialized.at("json_coefficients") == 0 &&
      std::string(lower_materialized.at("chart").as_string()) ==
          lower_chart &&
      std::string(upper_materialized.at("chart").as_string()) ==
          upper_chart &&
      lower_materialized.at("epsilon_min") == 0 &&
      lower_materialized.at("epsilon_max") ==
          (domain == "rational" ? 2 : 1) &&
      lower_materialized.at("strong_derivation_ownership") == true &&
      lower_materialized.at("residual_binding").as_object().at("status") ==
          "available" &&
      std::string(materialized_derivation.at("source_match").as_string()) ==
          lower_match &&
      materialized_derivation.at("planned_hop_provenance_identity") ==
          lower.at("planned_hop_provenance_identity") &&
      materialized_derivation.at("coefficient_transport") ==
          "native-retained-only" &&
      materialized_derivation.at("match_certified_complete_max") ==
          (domain == "rational" ? 2 : 1) &&
      materialized_derivation.at("whole_arm_complete") == false &&
      std::string(materialized_output.at("chart").as_string()) ==
          lower_chart &&
      std::string(materialized_output.at("checkpoint_identity").as_string()) ==
          domain + "-lower-receiving-local" &&
      lower_evaluation.at("status") == "ok" &&
      upper_evaluation.at("status") == "ok" &&
      lower_residual.at("status") == "ok" &&
      lower_residual.at("verdict") == "pass" &&
      std::abs(value_midpoint(lower_evaluation, 0) - 2.0) < 1e-25 &&
      std::abs(value_midpoint(upper_evaluation, 0) - 2.0) < 1e-25 &&
      saved.at("status") == "ok" && saved.at("locals") == 2 &&
      saved.at("planned_match_hops") == 2 &&
      saved.at("exact_matches") == 0 && saved.at("acb_matches") == 0 &&
      saved.at("tile_plans") == 0 &&
      saved.at("deferred_handle_kinds").as_array() ==
          json::array{"symbolic-local"} &&
      saved_payload.at("prepared_charts").as_array().size() == 3 &&
      saved_payload.at("retained_locals").as_array().size() == 5 &&
      planned_records.size() == 2 &&
      saved_payload.at("retained_tile_plans").as_array().size() == 1 &&
      saved_visibility.at("charts").as_array().empty() &&
      saved_visibility.at("locals").as_array().size() == 2 &&
      contains_string(saved_visibility.at("locals").as_array(), lower_local) &&
      contains_string(saved_visibility.at("locals").as_array(), upper_local) &&
      saved_visibility.at("matches").as_array().empty() &&
      saved_visibility.at("tile_plans").as_array().empty() &&
      lower_planned->at("schema") ==
          "diffexp2-retained-planned-match-hop-v2" &&
      embedded_native.at("schema") == expected_native_schema &&
      lower_planned->at("handoff") == lower.at("planned_hop") &&
      lower_planned->at("provenance_identity") ==
          lower.at("planned_hop_provenance_identity") &&
      (domain != "rational" ||
       (corruption.at("status") == "error" &&
        std::string(corruption.at("detail").as_string()).find(
            "unsupported native match kind") != std::string::npos)) &&
      restored.at("status") == "ok" &&
      restored.at("charts").as_array().empty() &&
      restored.at("locals").as_array().size() == 2 &&
      restored.at("exact_matches").as_array().empty() &&
      restored.at("acb_matches").as_array().empty() &&
      restored.at("planned_match_hops").as_array().empty() &&
      restored.at("tile_plans").as_array().empty() &&
      restored_lower_evaluation.at("status") == "ok" &&
      restored_upper_evaluation.at("status") == "ok" &&
      restored_lower_residual.at("status") == "ok" &&
      restored_lower_residual.at("verdict") == "pass" &&
      std::abs(value_midpoint(restored_lower_evaluation, 0) - 2.0) < 1e-25 &&
      std::abs(value_midpoint(restored_upper_evaluation, 0) - 2.0) < 1e-25 &&
      restored_local_stats.at("strong_derivation_ownership") == true &&
      std::string(restored_local_stats.at("retained_derivation").as_object()
                      .at("source_match").as_string()) == lower_match &&
      hidden_match.at("status") == "error" &&
      hidden_plan.at("status") == "error" &&
      hidden_source.at("status") == "error" &&
      resaved.at("status") == "ok" &&
      resaved.at("planned_match_hops") == 2 &&
      resaved_payload.at("retained_locals").as_array().size() == 5 &&
      resaved_payload.at("retained_planned_match_hops").as_array().size() == 2 &&
      resaved_payload.at("retained_tile_plans").as_array().size() == 1 &&
      next_anchor == "c:4" && next_lower_chart == "c:5" &&
      next_upper_chart == "c:6" && next_incoming == "l:6" &&
      next_basis == "l:7" && next_plan == "tile:2" &&
      next_match_result.at("match") == "m:3" &&
      stats.at("locals") == 4 && stats.at("tile_plans") == 1 &&
      stats.at("matches") == 1 && stats.at("local_matches") == 3 &&
      stats.at("pending_matches") == 0;

  if (!ok)
    std::cerr << domain << " lower: " << json::serialize(lower) << '\n'
              << domain << " upper: " << json::serialize(upper) << '\n'
              << domain << " lower materialized: "
              << json::serialize(lower_materialized) << '\n'
              << domain << " upper materialized: "
              << json::serialize(upper_materialized) << '\n'
              << domain << " lower evaluation: "
              << json::serialize(lower_evaluation) << '\n'
              << domain << " upper evaluation: "
              << json::serialize(upper_evaluation) << '\n'
              << domain << " plan stats: " << json::serialize(plan_stats)
              << '\n' << domain << " retained: "
              << json::serialize(retained_lower) << '\n'
              << domain << " saved: " << json::serialize(saved) << '\n'
              << domain << " restored: " << json::serialize(restored) << '\n'
              << domain << " restored local: "
              << json::serialize(restored_local_stats) << '\n'
              << domain << " corruption: " << json::serialize(corruption)
              << '\n'
              << domain << " next match: "
              << json::serialize(next_match_result) << '\n'
              << domain << " stats: " << json::serialize(stats) << '\n';

  (void)request(json::object{
      {"schema", 2}, {"op", "session.close"},
      {"session", restored_session}});
  std::error_code ignored;
  std::filesystem::remove(path, ignored);
  std::filesystem::remove(resaved_path, ignored);
  std::filesystem::remove(corrupt_path, ignored);
  return ok;
}

}  // namespace

int main() {
  const bool finite_pole = finite_laurent_pole_materialization();
  const bool rational = run_domain("rational");
  const bool acb = run_domain("acb");
  const bool ok = finite_pole && rational && acb;
  std::cout << (ok ? "PASS" : "FAIL")
            << ": plan-driven rational/Acb checkpoint ownership roundtrip\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
