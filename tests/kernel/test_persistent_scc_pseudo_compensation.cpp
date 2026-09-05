#include "diffexp/kernel/json_codec.hpp"
#include "diffexp/kernel/scalar.hpp"

#include <boost/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace json = boost::json;

namespace {

json::object request(const std::string& text) {
  return json::parse(diffexp::kernel::run_recurrence_json(text)).as_object();
}

json::object request(json::object value) {
  return request(json::serialize(value));
}

json::object metadata(const std::string& checkpoint,
                      const std::string& a,
                      const std::string& b) {
  return json::object{
      {"chart", json::object{{"center_exact", "0"},
                              {"scale_exact", "1"},
                              {"radius", "2"},
                              {"infinite_radius", false}}},
      {"tag", json::object{
          {"a", json::object{{"domain", "rational"},
                              {"canonical", a}}},
          {"b", json::object{{"domain", "rational"},
                              {"canonical", b}}},
          {"p", json::object{{"domain", "integer"},
                              {"canonical", "0"}}}}},
      {"prescriptions", json::array{}},
      {"checkpoint_identity", checkpoint}};
}

json::object scalar_run(bool seed, const std::string& a,
                        const std::string& b,
                        const std::string& root_a,
                        const std::string& root_b) {
  constexpr std::uint32_t nmax = 8;
  constexpr std::uint32_t width = 9;
  const diffexp::kernel::Rational aa(a), bb(b), ra(root_a), rb(root_b);
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    const auto an = aa + diffexp::kernel::Rational(std::to_string(n));
    const auto da = an - ra;
    const auto db = bb - rb;
    const auto kind = !da.is_zero() ? "T" : !db.is_zero() ? "P" : "R";
    shifts.push_back(json::string(an.str()));
    schedule.push_back(json::array{json::object{
        {"case", kind}, {"da", da.str()}, {"db", db.str()}}});
  }
  json::array initial;
  for (std::uint32_t epsilon = 0; epsilon < width; ++epsilon)
    initial.push_back(seed && epsilon == 4 ? "1" : "0");
  return json::object{
      {"nmax", nmax}, {"p", 0}, {"has_initial", seed},
      {"adaptive_probe", false}, {"a_target", a}, {"b_target", b},
      {"a_shift_min", 0}, {"a_shifts", std::move(shifts)},
      {"schedule", std::move(schedule)}, {"initial", std::move(initial)},
      {"initial_validity", json::array{
          seed ? json::value(4) : json::value(nullptr)}},
      {"source", nullptr}, {"return_u", false}};
}

json::object jordan_particular_run(bool adaptive_probe) {
  constexpr std::uint32_t nmax = 8;
  constexpr std::uint32_t width = 9;
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    shifts.push_back(json::string(
        n == 0 ? "1" : std::to_string(n + 1)));
    schedule.push_back(json::array{json::object{
        {"case", n == 0 ? "P" : "T"},
        {"da", std::to_string(n)}, {"db", "-1"}}});
  }
  json::array initial;
  for (std::uint32_t component = 0; component < 2; ++component)
    for (std::uint32_t epsilon = 0; epsilon < width; ++epsilon)
      initial.push_back("0");
  return json::object{
      {"nmax", nmax}, {"p", 0}, {"has_initial", false},
      {"adaptive_probe", adaptive_probe}, {"a_target", "1"},
      {"b_target", "0"}, {"a_shift_min", 0},
      {"a_shifts", std::move(shifts)}, {"schedule", std::move(schedule)},
      {"initial", std::move(initial)},
      {"initial_validity", json::array{nullptr, nullptr}},
      {"source", nullptr}, {"return_u", false}};
}

json::object casep_column(const std::string& checkpoint,
                          bool adaptive_probe) {
  json::array targets;
  targets.push_back(json::object{
      {"block", 1}, {"run", jordan_particular_run(adaptive_probe)},
      {"metadata", metadata(
          checkpoint + ":jordan-particular", "1", "0")}});
  targets.push_back(json::object{
      {"block", 2},
      {"run", scalar_run(false, "1", "0", "3/2", "3")},
      {"metadata", metadata(
          checkpoint + ":tail-particular", "1", "0")}});
  return json::object{
      {"checkpoint_identity", checkpoint},
      {"seed", json::object{
          {"block", 0}, {"run", scalar_run(true, "1", "0", "1", "0")},
          {"metadata", metadata(checkpoint + ":seed", "1", "0")}}},
      {"targets", std::move(targets)}};
}

std::string prepare_scalar_chart(const std::string& session,
                                 const std::string& key,
                                 const std::string& identity,
                                 const std::string& principal,
                                 const std::string& a,
                                 const std::string& b,
                                 bool no_pseudo) {
  const auto response = request(
      std::string("{\"schema\":2,\"op\":\"chart.prepare\",\"session\":\"") +
      session + "\",\"key\":\"" + key + "\",\"identity\":\"" + identity +
      "\",\"analytic\":{\"geometry\":{\"center_exact\":\"0\","
      "\"scale_exact\":\"1\",\"radius_exact\":\"2\","
      "\"infinite_radius\":false,\"prescriptions\":[]},"
      "\"principal_matrix\":[[{\"exact\":\"" + principal +
      "\",\"proven_zero\":false}]],\"native_scc_capabilities\":{"
      "\"regular\":false,\"identity_gauge\":true,\"identity_v\":true,"
      "\"no_pseudo\":" + (no_pseudo ? "true" : "false") +
      "}},\"scc\":{\"components\":[[0]],\"structural_edges\":[[0,0]],"
      "\"condensation_edges\":[],\"topological_order\":[0],"
      "\"coupling_depth\":0},\"problem\":{\"domain\":\"rational\","
      "\"d\":1,\"fb\":-4,\"w\":9,\"d_lags\":[[{\"s\":0,\"v\":\"1\"}]],"
      "\"denominators\":[],\"nhat_lags\":[{\"poly\":[{\"s\":0,"
      "\"e\":[[0,0,\"" + a + "\"]]},{\"s\":1,\"e\":[[0,0,\"" + b +
      "\"]]}],\"rat\":[],\"val\":[0]}],\"d0_inverse\":\"1\","
      "\"blocks\":[[0]],\"assembly\":{\"identity\":true,\"poly\":[],"
      "\"rat\":[],\"val\":[0]},\"chop_digits\":0}}");
  if (response.at("status") != "ok")
    throw std::runtime_error(json::serialize(response));
  return std::string(response.at("chart").as_string());
}

std::string prepare_jordan_chart(const std::string& session) {
  const auto response = request(std::string(R"json({
    "schema":2,"op":"chart.prepare","session":")json") + session +
    R"json(","key":"casep-jordan-key","identity":"casep-jordan",
    "analytic":{
      "geometry":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,"prescriptions":[]},
      "principal_matrix":[
        [{"exact":"(1+eps)/t","proven_zero":false},
         {"exact":"1/t","proven_zero":false}],
        [{"exact":"1","proven_zero":false},
         {"exact":"(1+eps)/t","proven_zero":false}]],
      "native_scc_capabilities":{"regular":false,"identity_gauge":true,
        "identity_v":true,"no_pseudo":false}},
    "scc":{"components":[[0,1]],
      "structural_edges":[[0,0],[0,1],[1,0],[1,1]],
      "condensation_edges":[],"topological_order":[0],"coupling_depth":0},
    "problem":{"domain":"rational","d":2,"fb":-4,"w":9,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[
        {"poly":[
          {"s":0,"e":[[0,0,"1"],[0,1,"1"],[1,1,"1"]]},
          {"s":1,"e":[[0,0,"1"],[1,1,"1"]]}],
         "rat":[],"val":[0,0,null,0]},
        {"poly":[{"s":0,"e":[[1,0,"1"]]}],
         "rat":[],"val":[null,null,0,null]}],
      "d0_inverse":"1","blocks":[[0,1]],
      "assembly":{"identity":true,"poly":[],"rat":[],
        "val":[0,null,null,0]},"chop_digits":0}
  })json");
  if (response.at("status") != "ok")
    throw std::runtime_error(json::serialize(response));
  return std::string(response.at("chart").as_string());
}

double real_midpoint(const json::value& coefficient) {
  return std::stod(
      std::string(coefficient.as_array().front().as_string()));
}

}  // namespace

int main() {
  constexpr std::size_t column_count = 4;
  const auto created = request(R"json({
    "schema":2,"op":"session.create","domain":"rational",
    "precision_bits":384,"output_digits":50,"scc_capacity":1,
    "local_capacity":4
  })json");
  const auto session = std::string(created.at("session").as_string());
  const auto seed_chart = prepare_scalar_chart(
      session, "casep-seed-key", "casep-seed", "1/t", "1", "0", true);
  const auto jordan_chart = prepare_jordan_chart(session);
  const auto tail_chart = prepare_scalar_chart(
      session, "casep-tail-key", "casep-tail",
      "((3/2)+(3*eps))/t", "3/2", "3", true);

  const auto prepared = request(std::string(R"json({
    "schema":2,"op":"scc.prepare","session":")json") + session +
    R"json(","key":"casep-composite-key","identity":"casep-parent-v1",
    "rational_shadow_identity":"casep-shadow-v1",
    "parent":{"dimension":4,
      "exact_system_record":[
        [{"exact":"1/t","proven_zero":false},
         {"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true}],
        [{"exact":"0","proven_zero":true},
         {"exact":"(1+eps)/t","proven_zero":false},
         {"exact":"1/t","proven_zero":false},
         {"exact":"0","proven_zero":true}],
        [{"exact":"1/t","proven_zero":false},
         {"exact":"1","proven_zero":false},
         {"exact":"(1+eps)/t","proven_zero":false},
         {"exact":"0","proven_zero":true}],
        [{"exact":"0","proven_zero":true},
         {"exact":"1/t","proven_zero":false},
         {"exact":"0","proven_zero":true},
         {"exact":"((3/2)+(3*eps))/t","proven_zero":false}]],
      "exact_theta_record":[
        [{"exact":"1","proven_zero":false},
         {"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true}],
        [{"exact":"0","proven_zero":true},
         {"exact":"1+eps","proven_zero":false},
         {"exact":"1","proven_zero":false},
         {"exact":"0","proven_zero":true}],
        [{"exact":"1","proven_zero":false},
         {"exact":"t","proven_zero":false},
         {"exact":"1+eps","proven_zero":false},
         {"exact":"0","proven_zero":true}],
        [{"exact":"0","proven_zero":true},
         {"exact":"1","proven_zero":false},
         {"exact":"0","proven_zero":true},
         {"exact":"(3/2)+(3*eps)","proven_zero":false}]],
      "chart":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,"prescriptions":[]},
      "scc":{"components":[[0],[1,2],[3]],
        "structural_edges":[[0,0],[0,2],[1,1],[1,2],[1,3],[2,1],[2,2],[3,3]],
        "condensation_edges":[[0,1],[1,2]],"topological_order":[0,1,2],
        "coupling_depth":2},
      "execution":{"mode":"BlockSequentialStrict","work_t_order":8},
      "work_contract":{"work_min":-4,"requested_min":-2,
        "requested_max":0,"work_complete_max":4,"public_t_order":0,
        "wolfram_coupling_depth":3}},
    "blocks":[
      {"block":0,"vertices":[0],"chart":")json" + seed_chart +
    R"json(","principal_identity":"casep-seed","regular":false,
       "identity_gauge":true,"identity_v":true,"no_pseudo":true},
      {"block":1,"vertices":[1,2],"chart":")json" + jordan_chart +
    R"json(","principal_identity":"casep-jordan","regular":false,
       "identity_gauge":true,"identity_v":true,"no_pseudo":false},
      {"block":2,"vertices":[3],"chart":")json" + tail_chart +
    R"json(","principal_identity":"casep-tail","regular":false,
       "identity_gauge":true,"identity_v":true,"no_pseudo":true}],
    "couplings":[
      {"source_block":0,"target_block":1,"source_vertices":[0],
       "target_vertices":[1,2],"rows":2,"columns":1,
       "exact_identity":"casep-seed-to-jordan","domain":"rational",
       "symbols":[],"entries":[
        {"row":1,"column":0,"source_vertex":0,"target_vertex":2,
         "exact_original_entry":"1/t","exact_theta_entry":"1",
         "multiplier":{"epsilon_shift":0,"center_pole_order":0,
          "kernels":[["1","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"]],
          "exact_identity":"1","proven_zero":false}}]},
      {"source_block":1,"target_block":2,"source_vertices":[1,2],
       "target_vertices":[3],"rows":1,"columns":2,
       "exact_identity":"casep-jordan-to-tail","domain":"rational",
       "symbols":[],"entries":[
        {"row":0,"column":0,"source_vertex":1,"target_vertex":3,
         "exact_original_entry":"1/t","exact_theta_entry":"1",
         "multiplier":{"epsilon_shift":0,"center_pole_order":0,
          "kernels":[["1","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0","0","0"]],
          "exact_identity":"1","proven_zero":false}}]}]
  })json");
  if (prepared.at("status") != "ok") {
    std::cerr << "prepare: " << json::serialize(prepared) << '\n';
    return EXIT_FAILURE;
  }
  const auto scc = std::string(prepared.at("scc").as_string());

  json::array columns;
  std::vector<std::string> checkpoints;
  for (std::size_t index = 0; index < column_count; ++index) {
    checkpoints.push_back("casep-column-" + std::to_string(index));
    // Reproduce the production double-box path: a cached homogeneous target
    // is first discovered through a submitted seed-style template with the
    // lower-frame probe disabled, then through derived particulars with it
    // enabled. The target cache must canonicalize this diagnostic policy.
    columns.push_back(casep_column(checkpoints.back(), index != 0));
  }
  const auto solved_batch = request(json::object{
      {"schema", 2}, {"op", "scc.solve_columns"}, {"session", session},
      {"scc", scc}, {"columns", std::move(columns)}, {"threads", 2}});
  const auto& solved_results = solved_batch.at("results").as_array();

  std::vector<json::object> evaluated;
  if (solved_batch.at("status") == "ok" &&
      solved_results.size() == column_count) {
    for (const auto& raw_solved : solved_results) {
      const auto& solved = raw_solved.as_object();
      evaluated.push_back(request(json::object{
          {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
          {"local", solved.at("local")},
          {"point", json::object{{"exact", "1/2"}}},
          {"options", json::object{{"tail_estimate", false}}}}));
    }
  }
  const auto stats = request(json::object{
      {"schema", 2}, {"op", "scc.stats"}, {"session", session},
      {"scc", scc}});

  bool collision_ok = solved_results.size() == column_count;
  bool checkpoints_ok = solved_results.size() == column_count;
  std::set<std::string> all_tail_tags;
  for (std::size_t solved_index = 0;
       solved_index < solved_results.size(); ++solved_index) {
    const auto& solved = solved_results[solved_index].as_object();
    checkpoints_ok = checkpoints_ok &&
        std::string(solved.at("checkpoint_identity").as_string()) ==
            checkpoints[solved_index];
    bool column_collision = false;
    std::set<std::string> tail_tags;
    for (const auto& raw : solved.at("block_diagnostics").as_array()) {
      const auto& diagnostic = raw.as_object();
      const auto block = diagnostic.at("block").as_int64();
      if (block == 1 && diagnostic.at("source_b") == "0") {
        column_collision = diagnostic.at("pseudo_hit_count") == 1 &&
            diagnostic.at("pseudo_compensation_count") == 2 &&
            diagnostic.at("max_pseudo_depth") == 2 &&
            diagnostic.at("pseudo_value_certified") == true &&
            diagnostic.at("uncompensated_pseudo_hit_count") == 0;
      }
      if (block == 2)
        tail_tags.insert(std::string(diagnostic.at("source_b").as_string()));
    }
    collision_ok = collision_ok && column_collision &&
        tail_tags == std::set<std::string>{"0", "1"};
    all_tail_tags.insert(tail_tags.begin(), tail_tags.end());
  }

  bool evaluated_finite = evaluated.size() == column_count;
  bool evaluated_identical = evaluated.size() == column_count;
  std::string reference_value;
  for (const auto& item : evaluated) {
    evaluated_finite = evaluated_finite && item.at("status") == "ok";
    if (!evaluated_finite) break;
    const auto& value = item.at("value").as_object();
    const auto serialized_value = json::serialize(value);
    if (reference_value.empty())
      reference_value = serialized_value;
    else
      evaluated_identical = evaluated_identical &&
          serialized_value == reference_value;
    const auto& coefficients = value.at("coefficients").as_array();
    evaluated_finite = evaluated_finite && value.at("min") == -2 &&
        coefficients.size() >= 8;
    for (std::size_t index = 0; index < 8 && evaluated_finite; ++index)
      evaluated_finite = std::abs(real_midpoint(coefficients[index])) < 1e-30;
  }

  const bool ok = solved_batch.at("status") == "ok" &&
      solved_batch.at("columns") == column_count &&
      solved_batch.at("worker_threads") == 2 &&
      std::all_of(solved_results.begin(), solved_results.end(),
          [](const auto& raw) {
            return raw.as_object().at("pseudo_hit_count") == 0;
          }) &&
      checkpoints_ok && collision_ok &&
      all_tail_tags == std::set<std::string>{"0", "1"} &&
      evaluated_finite && evaluated_identical &&
      stats.at("execution_scope") ==
          "exact-rational-regular-singular-jordan-block-dag-column-v2" &&
      stats.at("rational_shadow_identity") == "casep-shadow-v1" &&
      stats.at("casep_homogeneous_target_cache_scope") ==
          "immutable-composite" &&
      stats.at("casep_homogeneous_target_cache_entries") == 2 &&
      stats.at("casep_homogeneous_target_cache_builds") == 2 &&
      stats.at("casep_homogeneous_target_cache_hits") ==
          2 * (column_count - 1) &&
      stats.at("capability_evidence").as_object()
          .at("pseudo_schedule_execution") ==
          "exact-rational-joint-compensation-and-formal-overlap-certificate";
  if (!ok) {
    std::cerr << "solved: " << json::serialize(solved_batch) << '\n'
              << "evaluated:";
    for (const auto& item : evaluated)
      std::cerr << ' ' << json::serialize(item);
    std::cerr << '\n'
              << "stats: " << json::serialize(stats) << '\n';
  }
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                              {"session", session}});
  std::cout << (ok ? "PASS" : "FAIL")
            << ": mixed-probe CASE-P columns share two immutable targets\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
