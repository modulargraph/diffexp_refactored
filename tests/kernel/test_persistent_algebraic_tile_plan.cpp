#include "diffexp/kernel/checkpoint.hpp"
#include "diffexp/kernel/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp::kernel::run_recurrence_json(json::serialize(value)))
      .as_object();
}

json::array ball(const char* real) { return json::array{real, "0"}; }

json::object certified(const char* exact, const char* real, int sign) {
  return json::object{{"exact", exact}, {"value", ball(real)},
                      {"sign", sign}};
}

std::string prepare_chart(const std::string& session,
                          const std::string& key,
                          json::object geometry) {
  auto value = json::parse(R"json({
    "schema":2,"op":"chart.prepare","session":"placeholder",
    "key":"placeholder","identity":"placeholder",
    "analytic":{
      "geometry":{},
      "principal_matrix":[[{"exact":"0","proven_zero":true}]],
      "native_scc_capabilities":{"regular":true,"identity_gauge":true,
        "identity_v":true,"no_pseudo":true}},
    "scc":{"components":[[0]],"structural_edges":[],
      "condensation_edges":[],"topological_order":[0],
      "coupling_depth":0},
    "problem":{"domain":"acb","precision_bits":256,
      "d":1,"fb":0,"w":3,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],"val":[null]}],
      "d0_inverse":"1","blocks":[[0]],
      "assembly":{"identity":true,"poly":[],"rat":[],"val":[0]},
      "chop_digits":0}
  })json").as_object();
  value["session"] = session;
  value["key"] = key;
  value["identity"] = key + "-identity";
  value.at("analytic").as_object()["geometry"] = std::move(geometry);
  const auto response = request(std::move(value));
  if (response.at("status") != "ok")
    throw std::runtime_error("chart.prepare: " + json::serialize(response));
  return std::string(response.at("chart").as_string());
}

std::string solve_local(const std::string& session,
                        const std::string& chart,
                        const std::string& center,
                        const std::string& scale,
                        const std::string& radius_exact,
                        const std::string& checkpoint,
                        const std::string& value,
                        const std::string& radius_numeric = {}) {
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
                {"center_exact", center}, {"scale_exact", scale},
                {"radius_exact", radius_exact},
                {"radius", radius_numeric.empty()
                     ? json::value(radius_exact)
                     : json::value(radius_numeric)},
                {"infinite_radius", false}}},
           {"tag", json::object{
                {"a", json::object{{"domain", "rational"},
                                     {"canonical", "0"}}},
                {"b", json::object{{"domain", "rational"},
                                     {"canonical", "0"}}}}},
           {"prescriptions", json::array{}},
           {"checkpoint_identity", checkpoint}}}});
  if (response.at("status") != "ok")
    throw std::runtime_error("local.solve: " + json::serialize(response));
  return std::string(response.at("local").as_string());
}

json::object exact_lattice() {
  json::array row;
  row.push_back(json::object{
      {"min", 0}, {"max", 2},
      {"coefficients", json::array{"1", "0", "0"}}});
  json::array matrix;
  matrix.push_back(std::move(row));
  return json::object{
      {"schema", "diffexp3-exact-evaluated-epsilon-lattice-v1"},
      {"identity", "algebraic-transition-unit-lattice-v1"},
      {"evaluated_basis", std::move(matrix)}};
}

json::object geometry(const char* center_exact, const char* center_numeric,
                      const char* scale_exact, const char* scale_numeric,
                      const char* radius_exact, const char* radius_numeric) {
  return json::object{
      {"center_exact", center_exact},
      {"center_numeric", ball(center_numeric)},
      {"scale_exact", scale_exact},
      {"scale_numeric", ball(scale_numeric)},
      {"radius_exact", radius_exact},
      {"radius_numeric", ball(radius_numeric)},
      {"infinite_radius", false}, {"prescriptions", json::array{}}};
}

json::object topology() {
  return json::object{{"singular_points", json::array{}},
                      {"boundary_points", json::array{}},
                      {"complex_projections", json::array{}},
                      {"branch_sheets", json::array{}}};
}

json::object algebraic_arm(const std::string& anchor,
                           const std::string& receiver) {
  constexpr auto center_exact = "Sqrt[2]/2";
  constexpr auto scale_exact = "Sqrt[2]";
  constexpr auto center_value =
      "[0.7071067811865475244008443621 +/- 1e-27]";
  constexpr auto scale_value =
      "[1.4142135623730950488016887242 +/- 1e-27]";
  constexpr auto physical_match_exact = "1/2+Sqrt[2]/10";
  constexpr auto physical_match_value =
      "[0.6414213562373095048801688724 +/- 1e-27]";
  constexpr auto producing_match_exact = "1/4+Sqrt[2]/20";
  constexpr auto producing_match_value =
      "[0.3207106781186547524400844362 +/- 1e-27]";
  constexpr auto receiving_match_exact = "1/(2*Sqrt[2])-2/5";
  constexpr auto receiving_match_value =
      "[-0.04644660940672623779957781895 +/- 1e-27]";
  constexpr auto receiving_end_exact = "1/Sqrt[2]-1/2";
  constexpr auto receiving_end_value =
      "[0.2071067811865475244008443621 +/- 1e-27]";

  json::array planning{
      json::object{{"center_exact", "0"}, {"scale_exact", "2"},
                   {"radius_exact", "2"},
                   {"certificate_identity", "rational-anchor-v1"}},
      json::object{{"center_exact", "13/16"},
                   {"scale_exact", "11/8"},
                   {"radius_exact", "3/4"},
                   {"certificate_identity", "algebraic-surrogate-v1"}}};
  json::array charts{
      json::object{{"index", 0}, {"center", certified("0", "0", 0)},
                   {"scale", certified("2", "2", 1)},
                   {"radius", certified("2", "2", 1)}},
      json::object{{"index", 1},
                   {"center", certified(center_exact, center_value, 1)},
                   {"scale", certified(scale_exact, scale_value, 1)},
                   {"radius", certified(
                        "3/4", "[0.75 +/- 1e-27]", 1)}}};
  json::array matches{json::object{
      {"index", 0},
      {"physical", certified(physical_match_exact, physical_match_value, 1)},
      {"producing_local", certified(producing_match_exact,
                                      producing_match_value, 1)},
      {"receiving_local", certified(receiving_match_exact,
                                     receiving_match_value, -1)}}};
  json::array tiles{
      json::object{{"index", 0}, {"chart_index", 0},
                   {"physical_begin", certified("0", "0", 0)},
                   {"physical_end", certified(physical_match_exact,
                                                physical_match_value, 1)},
                   {"local_begin", certified("0", "0", 0)},
                   {"local_end", certified(producing_match_exact,
                                             producing_match_value, 1)}},
      json::object{{"index", 1}, {"chart_index", 1},
                   {"physical_begin", certified(physical_match_exact,
                                                  physical_match_value, 1)},
                   {"physical_end", certified("1", "1", 1)},
                   {"local_begin", certified(receiving_match_exact,
                                               receiving_match_value, -1)},
                   {"local_end", certified(receiving_end_exact,
                                             receiving_end_value, 1)}}};
  json::object certified_geometry{
      {"schema", "diffexp3-wolfram-certified-algebraic-arm-v1"},
      {"exact_identity", "algebraic-arm-exact-v1"},
      {"charts", std::move(charts)}, {"matches", std::move(matches)},
      {"tiles", std::move(tiles)}};
  return json::object{{"from_exact", "0"}, {"to_exact", "1"},
                      {"charts", json::array{anchor, receiver}},
                      {"planning_charts", std::move(planning)},
                      {"certified_geometry", std::move(certified_geometry)},
                      {"topology", topology()}};
}

json::object center_match_arm(const std::string& anchor,
                              const std::string& receiver) {
  auto arm = algebraic_arm(anchor, receiver);
  auto& certified_geometry = arm.at("certified_geometry").as_object();
  certified_geometry["exact_identity"] = "center-match-arm-exact-v1";
  auto& match = certified_geometry.at("matches").as_array().front().as_object();
  match["physical"] = certified("0", "0", 0);
  match["producing_local"] = certified("0", "0", 0);
  match["receiving_local"] = certified("-1/2", "-1/2", -1);
  auto& tiles = certified_geometry.at("tiles").as_array();
  auto& producing = tiles.front().as_object();
  producing["physical_end"] = certified("0", "0", 0);
  producing["local_end"] = certified("0", "0", 0);
  auto& receiving = tiles.back().as_object();
  receiving["physical_begin"] = certified("0", "0", 0);
  receiving["local_begin"] = certified("-1/2", "-1/2", -1);
  return arm;
}

}  // namespace

int main() {
  const std::string path =
      "/tmp/diffexp3-persistent-algebraic-tile-plan.de2cp";
  std::remove(path.c_str());
  std::string session;
  std::string restored_session;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "acb"}, {"precision_bits", 256},
        {"output_digits", 50}, {"chart_capacity", 2},
        {"local_capacity", 3}, {"match_capacity", 2},
        {"tile_plan_capacity", 2}});
    session = std::string(created.at("session").as_string());
    const auto anchor = prepare_chart(
        session, "algebraic-anchor",
        geometry("0", "0", "2", "2", "2", "2"));
    const auto receiver = prepare_chart(
        session, "algebraic-receiver",
        geometry("Sqrt[2]/2",
                 "[0.7071067811865475244008443621 +/- 1e-27]",
                 "Sqrt[2]",
                 "[1.4142135623730950488016887242 +/- 1e-27]",
                 "3/4", "[0.75 +/- 1e-30]"));
    const auto incoming = solve_local(
        session, anchor, "0", "2", "2", "algebraic-anchor-local", "2");
    const auto basis = solve_local(
        session, receiver, "Sqrt[2]/2", "Sqrt[2]", "3/4",
        "algebraic-receiver-local", "1", "[0.75 +/- 1e-25]");
    auto arm = algebraic_arm(anchor, receiver);
    const auto plan = request(json::object{
        {"schema", 2}, {"op", "tile.plan_arm"}, {"session", session},
        {"checkpoint_identity", "algebraic-plan-checkpoint-v1"},
        {"division_order", 3}, {"arm", arm}});
    if (plan.at("status") != "ok")
      throw std::runtime_error("tile.plan_arm: " + json::serialize(plan));
    const auto plan_handle = std::string(plan.at("tile_plan").as_string());
    const auto& published_arm = plan.at("arm").as_object();
    const auto& provenance_reference =
        plan.at("provenance_reference").as_object();
    const bool bounded_publication =
        plan.if_contains("provenance_identity") == nullptr &&
        provenance_reference.at("algorithm") == "fnv1a64-v1" &&
        provenance_reference.at("identity_bytes").as_int64() > 0 &&
        published_arm.if_contains("planning_charts") == nullptr &&
        published_arm.if_contains("certified_geometry") == nullptr &&
        published_arm.if_contains("charts") == nullptr &&
        published_arm.if_contains("matches") == nullptr &&
        published_arm.if_contains("tiles") == nullptr &&
        published_arm.at("chart_count") == 2 &&
        published_arm.at("match_count") == 1 &&
        published_arm.at("tile_count") == 2 &&
        published_arm.at("exact_intervals_published") == false &&
        json::serialize(plan).size() < 16384;
    const auto match = request(json::object{
        {"schema", 2}, {"op", "tile.match_interval"},
        {"session", session}, {"tile_plan", plan_handle},
        {"arm", "upper"}, {"match", 0}});
    const auto hop = request(json::object{
        {"schema", 2}, {"op", "tile.match_advance"},
        {"session", session}, {"tile_plan", plan_handle},
        {"arm", "upper"}, {"match", 0},
        {"basis", json::array{basis}}, {"incoming", incoming},
        {"epsilon", json::object{{"min", 0}, {"max", 2},
                                   {"required_complete_max", 1}}},
        {"checkpoint_identity", "algebraic-transition-hop-v1"},
        {"exact_lattice", exact_lattice()},
        {"refinement", json::object{
             {"relative_tolerance", "1e-40"}, {"max_steps", 1}}}});
    if (hop.at("status") != "ok")
      throw std::runtime_error("tile.match_advance: " +
                               json::serialize(hop));

    auto bad_arm = arm;
    bad_arm.at("certified_geometry").as_object().at("charts").as_array()[1]
        .as_object().at("center").as_object()["exact"] = "Sqrt[3]/2";
    const auto rejected = request(json::object{
        {"schema", 2}, {"op", "tile.plan_arm"}, {"session", session},
        {"checkpoint_identity", "bad-algebraic-plan"},
        {"division_order", 3}, {"arm", std::move(bad_arm)}});
    if (rejected.at("status") != "error")
      throw std::runtime_error(
          "mismatched algebraic owner identity was accepted");

    auto bad_radius_arm = arm;
    bad_radius_arm.at("certified_geometry").as_object()
        .at("charts").as_array()[1].as_object()
        .at("radius").as_object()["exact"] = "4/5";
    const auto radius_rejected = request(json::object{
        {"schema", 2}, {"op", "tile.plan_arm"}, {"session", session},
        {"checkpoint_identity", "bad-algebraic-radius-plan"},
        {"division_order", 3}, {"arm", std::move(bad_radius_arm)}});
    if (radius_rejected.at("status") != "error")
      throw std::runtime_error(
          "mismatched exact radius identity was accepted");

    const auto center_plan = request(json::object{
        {"schema", 2}, {"op", "tile.plan_arm"}, {"session", session},
        {"checkpoint_identity", "center-match-plan-v1"},
        {"division_order", 3},
        {"arm", center_match_arm(anchor, receiver)}});
    if (center_plan.at("status") != "ok")
      throw std::runtime_error("center tile.plan_arm: " +
                               json::serialize(center_plan));
    const auto center_hop = request(json::object{
        {"schema", 2}, {"op", "tile.match_advance"},
        {"session", session},
        {"tile_plan", center_plan.at("tile_plan")},
        {"arm", "upper"}, {"match", 0},
        {"basis", json::array{basis}}, {"incoming", incoming},
        {"epsilon", json::object{{"min", 0}, {"max", 2},
                                   {"required_complete_max", 1}}},
        {"checkpoint_identity", "regular-center-transition-hop-v1"},
        {"exact_lattice", exact_lattice()},
        {"refinement", json::object{
             {"relative_tolerance", "1e-40"}, {"max_steps", 1}}}});
    if (center_hop.at("status") != "ok")
      throw std::runtime_error("center tile.match_advance: " +
                               json::serialize(center_hop));

    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
        {"path", path}, {"checkpoint_identity", "algebraic-roundtrip-v1"}});
    if (saved.at("status") != "ok")
      throw std::runtime_error("checkpoint.save: " + json::serialize(saved));
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    session.clear();
    const auto restored = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"}, {"path", path},
        {"expected_identity", "algebraic-roundtrip-v1"}});
    if (restored.at("status") != "ok")
      throw std::runtime_error("checkpoint.restore: " +
                               json::serialize(restored));
    restored_session = std::string(restored.at("session").as_string());
    const auto stats = request(json::object{
        {"schema", 2}, {"op", "tile.stats"},
        {"session", restored_session}, {"tile_plan", plan_handle}});
    const auto& hop_geometry =
        hop.at("planned_hop").as_object().at("geometry").as_object();
    const bool ok = bounded_publication && match.at("status") == "ok" &&
        match.at("physical_exact") != "1/2" &&
        hop.at("physical_match_point_exact") == "1/2+Sqrt[2]/10" &&
        hop.at("basis_point_exact") == "1/(2*Sqrt[2])-2/5" &&
        hop.at("incoming_point_exact") == "1/4+Sqrt[2]/20" &&
        hop_geometry.at("physical_exact") != "1/2+Sqrt[2]/10" &&
        center_hop.at("physical_match_point_exact") == "0" &&
        center_hop.at("incoming_point_exact") == "0" &&
        center_hop.at("basis_point_exact") == "-1/2" &&
        stats.at("status") == "ok" && stats.at("matches") == 1 &&
        stats.at("tiles") == 2;
    // The public interval remains the rational planner's combinatorial
    // witness; the certified physical identity is consumed only by matching
    // and line integration, so it must not replace that topology record.
    if (!ok)
      std::cerr << "plan=" << json::serialize(plan)
                << " match=" << json::serialize(match)
                << " hop=" << json::serialize(hop)
                << " center_hop=" << json::serialize(center_hop)
                << " stats=" << json::serialize(stats) << '\n';
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", restored_session}});
    restored_session.clear();
    std::remove(path.c_str());
    std::cout << (ok ? "PASS" : "FAIL")
              << ": exact radius identity survives unequal rigorous specializations\n";
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    if (!restored_session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored_session}});
    std::remove(path.c_str());
    std::cerr << "FAIL: persistent algebraic tile plan: "
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
