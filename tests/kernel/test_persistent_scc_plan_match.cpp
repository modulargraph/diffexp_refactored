#include "diffexp/kernel/checkpoint.hpp"
#include "diffexp/kernel/json_codec.hpp"

#include <boost/json.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp::kernel::run_recurrence_json(json::serialize(value)))
      .as_object();
}

json::array prescriptions() {
  return json::array{json::object{
      {"factor_exact", "mixed-scc-f"}, {"sign", -1},
      {"multiplicity", 1}, {"leading_coefficient_sign", 1}}};
}

json::object geometry(const std::string& center) {
  return json::object{{"center_exact", center}, {"scale_exact", "1"},
                      {"radius_exact", "2"},
                      {"infinite_radius", false},
                      {"prescriptions", prescriptions()}};
}

json::array zero_matrix(std::uint32_t dimension) {
  json::array matrix;
  for (std::uint32_t row = 0; row < dimension; ++row) {
    json::array encoded_row;
    for (std::uint32_t column = 0; column < dimension; ++column)
      encoded_row.push_back(
          json::object{{"exact", "0"}, {"proven_zero", true}});
    matrix.push_back(std::move(encoded_row));
  }
  return matrix;
}

json::array singleton_components(std::uint32_t dimension) {
  json::array components;
  for (std::uint32_t component = 0; component < dimension; ++component)
    components.push_back(json::array{component});
  return components;
}

std::string prepare_chart(const std::string& session,
                          const std::string& key,
                          const std::string& identity,
                          const std::string& center,
                          std::uint32_t dimension) {
  json::array d_lags;
  d_lags.push_back(json::array{json::object{{"s", 0}, {"v", "1"}}});
  json::array null_matrix;
  json::array assembly_matrix;
  for (std::uint32_t row = 0; row < dimension; ++row) {
    for (std::uint32_t column = 0; column < dimension; ++column) {
      null_matrix.push_back(nullptr);
      assembly_matrix.push_back(row == column ? json::value(0)
                                               : json::value(nullptr));
    }
  }
  const auto components = singleton_components(dimension);
  const auto response = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", key}, {"identity", identity},
      {"analytic",
       json::object{
           {"geometry", geometry(center)},
           {"principal_matrix", zero_matrix(dimension)},
           {"native_scc_capabilities",
            json::object{{"regular", true}, {"identity_gauge", true},
                         {"identity_v", true}, {"no_pseudo", true}}}}},
      {"scc",
       json::object{{"components", components},
                    {"structural_edges", json::array{}},
                    {"condensation_edges", json::array{}},
                    {"topological_order", [&] {
                       json::array order;
                       for (std::uint32_t index = 0; index < dimension;
                            ++index)
                         order.push_back(index);
                       return order;
                     }()},
                    {"coupling_depth", 0}}},
      {"problem",
       json::object{
           {"domain", "rational"}, {"d", dimension}, {"fb", 0},
           {"w", 3}, {"d_lags", std::move(d_lags)},
           {"denominators", json::array{}},
           {"nhat_lags",
            json::array{json::object{{"poly", json::array{}},
                                     {"rat", json::array{}},
                                     {"val", std::move(null_matrix)}}}},
           {"d0_inverse", "1"},
           {"blocks", singleton_components(dimension)},
           {"assembly",
            json::object{{"identity", true}, {"poly", json::array{}},
                         {"rat", json::array{}},
                         {"val", std::move(assembly_matrix)}}},
           {"chop_digits", 0}}}});
  if (response.at("status") != "ok")
    throw std::runtime_error("chart.prepare: " + json::serialize(response));
  return std::string(response.at("chart").as_string());
}

json::object regular_scalar_run() {
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
      {"initial", json::array{"1", "0", "0"}},
      {"initial_validity", json::array{2}}, {"source", nullptr},
      {"return_u", false}};
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
      {"prescriptions", prescriptions()},
      {"checkpoint_identity", checkpoint}};
}

std::string prepare_scc(const std::string& session,
                        const std::string& first,
                        const std::string& second) {
  json::array blocks;
  blocks.push_back(json::object{
      {"block", 0}, {"vertices", json::array{0}}, {"chart", first},
      {"principal_identity", "mixed-scc-block-0"}, {"regular", true},
      {"identity_gauge", true}, {"identity_v", true},
      {"no_pseudo", true}});
  blocks.push_back(json::object{
      {"block", 1}, {"vertices", json::array{1}}, {"chart", second},
      {"principal_identity", "mixed-scc-block-1"}, {"regular", true},
      {"identity_gauge", true}, {"identity_v", true},
      {"no_pseudo", true}});
  const auto response = request(json::object{
      {"schema", 2}, {"op", "scc.prepare"}, {"session", session},
      {"key", "mixed-scc-owner"}, {"identity", "mixed-scc-parent"},
      {"parent",
       json::object{
           {"dimension", 2}, {"exact_system_record", zero_matrix(2)},
           {"exact_theta_record", zero_matrix(2)},
           {"chart", geometry("-3/5")},
           {"scc", json::object{{"components", singleton_components(2)},
                                 {"structural_edges", json::array{}},
                                 {"condensation_edges", json::array{}},
                                 {"topological_order", json::array{0, 1}},
                                 {"coupling_depth", 0}}},
           {"execution",
            json::object{{"mode", "BlockSequentialStrict"},
                         {"work_t_order", 4}}},
           {"work_contract",
            json::object{{"work_min", 0}, {"requested_min", 0},
                         {"requested_max", 2}, {"work_complete_max", 2},
                         {"public_t_order", 0},
                         {"wolfram_coupling_depth", 1}}}}},
      {"blocks", std::move(blocks)}, {"couplings", json::array{}}});
  if (response.at("status") != "ok")
    throw std::runtime_error("scc.prepare: " + json::serialize(response));
  return std::string(response.at("scc").as_string());
}

std::vector<std::string> solve_scc_basis(const std::string& session,
                                         const std::string& scc) {
  json::array columns;
  for (std::uint32_t block = 0; block < 2; ++block) {
    const auto checkpoint = "mixed-scc-basis-" + std::to_string(block);
    columns.push_back(json::object{
        {"checkpoint_identity", checkpoint},
        {"seed", json::object{{"block", block},
                               {"run", regular_scalar_run()},
                               {"metadata", metadata("-3/5", checkpoint)}}},
        {"targets", json::array{}}});
  }
  const auto response = request(json::object{
      {"schema", 2}, {"op", "scc.solve_columns"}, {"session", session},
      {"scc", scc}, {"columns", std::move(columns)}, {"threads", 2}});
  if (response.at("status") != "ok")
    throw std::runtime_error("scc.solve_columns: " +
                             json::serialize(response));
  std::vector<std::string> result;
  for (const auto& raw : response.at("results").as_array())
    result.emplace_back(raw.as_object().at("local").as_string());
  return result;
}

std::string solve_incoming(const std::string& session,
                           const std::string& chart) {
  json::array schedule_row;
  schedule_row.push_back(json::object{{"case", "R"}, {"da", "0"},
                                      {"db", "0"}});
  schedule_row.push_back(json::object{{"case", "R"}, {"da", "0"},
                                      {"db", "0"}});
  json::array schedule;
  schedule.push_back(std::move(schedule_row));
  const auto response = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart},
      {"run", json::object{{"nmax", 0}, {"p", 0},
                            {"has_initial", true},
                            {"adaptive_probe", false},
                            {"a_target", "0"}, {"b_target", "0"},
                            {"a_shift_min", 0},
                            {"a_shifts", json::array{"0"}},
                            {"schedule", std::move(schedule)},
                            {"initial", json::array{"2", "0", "0",
                                                     "3", "0", "0"}},
                            {"initial_validity", json::array{2, 2}},
                            {"source", nullptr}, {"return_u", false}}},
      {"metadata", metadata("0", "mixed-scc-incoming")}});
  if (response.at("status") != "ok")
    throw std::runtime_error("local.solve: " + json::serialize(response));
  return std::string(response.at("local").as_string());
}

json::object topology(bool singular_endpoint) {
  return json::object{
      {"singular_points",
       singular_endpoint ? json::array{"-3/5"} : json::array{}},
      {"boundary_points", json::array{}},
      {"complex_projections", json::array{}},
      {"branch_sheets", json::array{json::object{
                            {"factor_exact", "mixed-scc-f"},
                            {"sign", -1}}}}};
}

json::object arm(const std::string& endpoint, const std::string& anchor,
                 const std::string& owner, bool singular_endpoint) {
  return json::object{{"from_exact", "0"}, {"to_exact", endpoint},
                      {"charts", json::array{anchor, owner}},
                      {"topology", topology(singular_endpoint)}};
}

double midpoint(const json::value& coefficient) {
  return std::stod(
      std::string(coefficient.as_array().front().as_string()));
}

}  // namespace

int main() {
  const std::string path = "/tmp/diffexp3-mixed-scc-plan-owner.chk";
  std::remove(path.c_str());
  std::string original_session;
  std::string restored_session;
  std::string stage = "create";
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"}, {"domain", "rational"},
        {"precision_bits", 256}, {"output_digits", 30},
        {"chart_capacity", 4}, {"scc_capacity", 1},
        {"local_capacity", 5}, {"match_capacity", 1},
        {"tile_plan_capacity", 1}});
    original_session = std::string(created.at("session").as_string());

    stage = "prepare";
    const auto anchor = prepare_chart(original_session, "mixed-anchor",
                                      "mixed-anchor-identity", "0", 2);
    const auto block0 = prepare_chart(original_session, "mixed-block-0",
                                      "mixed-scc-block-0", "-3/5", 1);
    const auto block1 = prepare_chart(original_session, "mixed-block-1",
                                      "mixed-scc-block-1", "-3/5", 1);
    const auto upper = prepare_chart(original_session, "mixed-upper",
                                     "mixed-upper-identity", "2/3", 2);
    const auto scc = prepare_scc(original_session, block0, block1);
    const auto basis = solve_scc_basis(original_session, scc);
    const auto incoming = solve_incoming(original_session, anchor);

    stage = "plan";
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", original_session},
        {"checkpoint_identity", "mixed-scc-plan"}, {"division_order", 3},
        {"lower", arm("-3/5", anchor, scc, true)},
        {"upper", arm("2/3", anchor, upper, false)}});
    if (planned.at("status") != "ok")
      throw std::runtime_error("tile.plan: " + json::serialize(planned));
    const auto plan = std::string(planned.at("tile_plan").as_string());
    const auto plan_stats = request(json::object{
        {"schema", 2}, {"op", "tile.stats"},
        {"session", original_session}, {"tile_plan", plan}});
    if (plan_stats.at("status") != "ok")
      throw std::runtime_error("tile.stats: " + json::serialize(plan_stats));
    const auto planned_scc = plan_stats.at("lower").as_object()
                                  .at("charts").as_array()[1].as_object();

    stage = "release-owners";
    (void)request(json::object{{"schema", 2}, {"op", "scc.release"},
                               {"session", original_session}, {"scc", scc}});
    for (const auto& chart : {anchor, block0, block1, upper})
      (void)request(json::object{{"schema", 2}, {"op", "chart.release"},
                                 {"session", original_session},
                                 {"chart", chart}});

    stage = "checkpoint";
    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", original_session}, {"path", path},
        {"checkpoint_identity", "mixed-scc-owner-checkpoint"}});
    if (saved.at("status") != "ok")
      throw std::runtime_error("checkpoint.save: " + json::serialize(saved));
    const auto payload = json::parse(
        diffexp::kernel::checkpoint::read(path).payload_json).as_object();
    const auto& visibility = payload.at("session").as_object()
                                 .at("registry_visibility").as_object();

    stage = "restore";
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", original_session}});
    original_session.clear();
    const auto restored = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"}, {"path", path},
        {"expected_identity", "mixed-scc-owner-checkpoint"}});
    if (restored.at("status") != "ok")
      throw std::runtime_error("checkpoint.restore: " +
                               json::serialize(restored));
    restored_session = std::string(restored.at("session").as_string());
    const auto stats = request(json::object{
        {"schema", 2}, {"op", "session.stats"},
        {"session", restored_session}});
    const auto retained_plan = request(json::object{
        {"schema", 2}, {"op", "tile.stats"}, {"session", restored_session},
        {"tile_plan", plan}});

    stage = "match";
    const auto matched = request(json::object{
        {"schema", 2}, {"op", "tile.match_advance"},
        {"session", restored_session}, {"tile_plan", plan},
        {"arm", "lower"}, {"match", 0},
        {"basis", json::array{basis[0], basis[1]}},
        {"incoming", incoming},
        {"epsilon", json::object{{"min", 0}, {"max", 2},
                                  {"required_complete_max", 2}}},
        {"checkpoint_identity", "mixed-scc-hop"}});
    if (matched.at("status") != "ok")
      throw std::runtime_error("tile.match_advance: " +
                               json::serialize(matched));
    const auto match = std::string(matched.at("match").as_string());

    stage = "materialize";
    const auto materialized = request(json::object{
        {"schema", 2}, {"op", "match.materialize_local"},
        {"session", restored_session}, {"match", match},
        {"checkpoint_identity", "mixed-scc-materialized"}});
    if (materialized.at("status") != "ok")
      throw std::runtime_error("match.materialize_local: " +
                               json::serialize(materialized));
    const auto local = std::string(materialized.at("local").as_string());
    const auto evaluated = request(json::object{
        {"schema", 2}, {"op", "local.evaluate"},
        {"session", restored_session}, {"local", local},
        {"point", json::object{{"exact", "1/5"}}},
        {"options", json::object{{"tail_estimate", false}}}});

    const auto& restored_lower = retained_plan.at("lower").as_object();
    const auto& restored_scc =
        restored_lower.at("charts").as_array()[1].as_object();
    const auto& hop = matched.at("planned_hop").as_object();
    const auto& receiving = hop.at("receiving").as_object();
    const auto& derivation =
        materialized.at("retained_derivation").as_object();
    const auto& coefficients =
        evaluated.at("value").as_object().at("coefficients").as_array();
    const bool ok =
        std::string(planned_scc.at("chart").as_string()) == scc &&
        planned_scc.at("identity") == "mixed-scc-parent" &&
        planned_scc.at("singular_center") == true &&
        planned_scc.at("prescriptions").as_array().size() == 1 &&
        saved.at("charts") == 0 && saved.at("sccs") == 0 &&
        payload.at("prepared_charts").as_array().size() == 4 &&
        payload.at("prepared_scc").as_array().size() == 1 &&
        visibility.at("charts").as_array().empty() &&
        visibility.at("sccs").as_array().empty() &&
        restored.at("charts").as_array().empty() &&
        restored.at("sccs").as_array().empty() &&
        restored.at("locals").as_array().size() == 3 &&
        restored.at("tile_plans").as_array().size() == 1 &&
        stats.at("charts") == 0 && stats.at("scc_charts") == 0 &&
        std::string(restored_scc.at("chart").as_string()) == scc &&
        restored_scc.at("identity") == "mixed-scc-parent" &&
        restored_scc.at("singular_center") == true &&
        std::string(receiving.at("chart").as_string()) == scc &&
        receiving.at("chart_identity") == "mixed-scc-parent" &&
        receiving.at("basis").as_array().size() == 2 &&
        std::string(materialized.at("chart").as_string()) == scc &&
        derivation.size() == 5 &&
        derivation.at("schema") ==
            "diffexp3-retained-plan-match-local-materialization-v1" &&
        derivation.at("capability") ==
            "retained-native-plan-match-local-materialization-v1" &&
        derivation.at("checkpoint_identity") == "mixed-scc-materialized" &&
        derivation.at("ownership") == "strong" &&
        evaluated.at("status") == "ok" && coefficients.size() == 6 &&
        std::abs(midpoint(coefficients[0]) - 2.0) < 1e-25 &&
        std::abs(midpoint(coefficients[1]) - 3.0) < 1e-25;

    if (!ok)
      std::cerr << "planned: " << json::serialize(planned) << '\n'
                << "saved: " << json::serialize(saved) << '\n'
                << "payload: " << json::serialize(payload) << '\n'
                << "restored: " << json::serialize(restored) << '\n'
                << "stats: " << json::serialize(stats) << '\n'
                << "retained plan: " << json::serialize(retained_plan)
                << '\n' << "matched: " << json::serialize(matched) << '\n'
                << "materialized: " << json::serialize(materialized) << '\n'
                << "evaluated: " << json::serialize(evaluated) << '\n';

    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", restored_session}});
    restored_session.clear();
    std::remove(path.c_str());
    std::cout << (ok ? "PASS" : "FAIL")
              << ": mixed prepared/SCC tile owner plan-match materialization\n";
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    if (!original_session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", original_session}});
    if (!restored_session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored_session}});
    std::remove(path.c_str());
    std::cerr << "FAIL at " << stage << ": " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
