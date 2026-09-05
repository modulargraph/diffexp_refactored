#include "diffexp/kernel/json_codec.hpp"

#include <boost/json.hpp>

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

std::string prepare_regular_chart(const std::string& session) {
  auto response = request(json::parse(std::string(R"json({
    "schema":2,"op":"chart.prepare","session":")json") + session +
      R"json(","key":"exact-regular-local-match-chart",
    "identity":"exact-regular-local-match-chart-v1",
    "analytic":{
      "geometry":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,"prescriptions":[]},
      "principal_matrix":[
        [{"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true}],
        [{"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true}]],
      "native_scc_capabilities":{"regular":true,"identity_gauge":true,
        "identity_v":true,"no_pseudo":true}},
    "scc":{"components":[[0],[1]],"structural_edges":[],
      "condensation_edges":[],"topological_order":[0,1],
      "coupling_depth":0},
    "problem":{"domain":"rational","d":2,"fb":0,"w":3,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],
        "val":[null,null,null,null]}],
      "d0_inverse":"1","blocks":[[0],[1]],
      "assembly":{"identity":true,"poly":[],"rat":[],
        "val":[0,null,null,0]},"chop_digits":0}
  })json").as_object());
  if (response.at("status") != "ok")
    throw std::runtime_error(json::serialize(response));
  return std::string(response.at("chart").as_string());
}

json::object regular_run(std::vector<std::string> initial) {
  json::array encoded_initial;
  for (auto& value : initial)
    encoded_initial.push_back(json::string(std::move(value)));
  json::array schedule_row;
  schedule_row.push_back(
      json::object{{"case", "R"}, {"da", "0"}, {"db", "0"}});
  schedule_row.push_back(
      json::object{{"case", "R"}, {"da", "0"}, {"db", "0"}});
  json::array schedule;
  schedule.push_back(std::move(schedule_row));
  return json::object{
      {"nmax", 0}, {"p", 0}, {"has_initial", true},
      {"adaptive_probe", false}, {"a_target", "0"},
      {"b_target", "0"}, {"a_shift_min", 0},
      {"a_shifts", json::array{"0"}},
      {"schedule", std::move(schedule)},
      {"initial", std::move(encoded_initial)},
      {"initial_validity", json::array{2, 2}},
      {"source", nullptr}, {"return_u", false}};
}

json::object metadata(const std::string& checkpoint) {
  return json::object{
      {"chart", json::object{{"center_exact", "0"}, {"scale_exact", "1"},
                              {"radius", "2"},
                              {"infinite_radius", false}}},
      {"tag", json::object{
          {"a", json::object{{"domain", "rational"}, {"canonical", "0"}}},
          {"b", json::object{{"domain", "rational"}, {"canonical", "0"}}}}},
      {"prescriptions", json::array{}},
      {"checkpoint_identity", checkpoint}};
}

std::string solve_local(const std::string& session, const std::string& chart,
                        std::vector<std::string> initial,
                        const std::string& checkpoint) {
  auto response = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart}, {"run", regular_run(std::move(initial))},
      {"metadata", metadata(checkpoint)}});
  if (response.at("status") != "ok")
    throw std::runtime_error(json::serialize(response));
  return std::string(response.at("local").as_string());
}

}  // namespace

int main() {
  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", "rational"},
      {"precision_bits", 256}, {"output_digits", 30},
      {"local_capacity", 3}, {"match_capacity", 1}});
  const auto session = std::string(created.at("session").as_string());
  const auto chart = prepare_regular_chart(session);

  // F(eps) = {{1,1},{0,eps}} has one epsilon-lattice saturation action.
  const auto first = solve_local(
      session, chart, {"1", "0", "0", "0", "0", "0"}, "basis-0");
  const auto second = solve_local(
      session, chart, {"1", "0", "0", "0", "1", "0"}, "basis-1");
  // v = first + 2 second = {3,2 eps}; the exact original-frame weights are
  // {1,2} after the saturation transformation is mapped back.
  const auto incoming = solve_local(
      session, chart, {"3", "0", "0", "0", "2", "0"}, "incoming");

  const auto matched = request(json::object{
      {"schema", 2}, {"op", "local.match"}, {"session", session},
      {"basis", json::array{first, second}}, {"incoming", incoming},
      {"basis_chart", chart}, {"incoming_chart", chart},
      {"basis_point", json::object{{"exact", "1/2"}}},
      {"incoming_point", json::object{{"exact", "1/2"}}},
      {"epsilon", json::object{{"min", 0}, {"max", 2},
                                {"required_complete_max", 1}}},
      {"basis_checkpoint_identities", json::array{"basis-0", "basis-1"}},
      {"incoming_checkpoint_identity", "incoming"},
      {"checkpoint_identity", "exact-regular-match-state-v1"}});

  const auto match = std::string(matched.at("match").as_string());
  for (const auto& local : {first, second, incoming})
    (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                               {"session", session}, {"local", local}});
  const auto retained = request(json::object{
      {"schema", 2}, {"op", "match.stats"}, {"session", session},
      {"match", match}});
  const auto stats = request(json::object{
      {"schema", 2}, {"op", "session.stats"}, {"session", session}});

  const auto& residual = matched.at("residual").as_object();
  const bool ok = created.at("status") == "ok" &&
      created.at("local_match_capability") ==
          "exact-rational-regular-local-match-v1" &&
      matched.at("status") == "ok" && matched.at("native_retained") == true &&
      matched.at("json_coefficients") == 0 && matched.at("dimension") == 2 &&
      matched.at("normalized_determinant_valuation") == 1 &&
      matched.at("initial_leading_rank") == 1 &&
      matched.at("final_leading_rank") == 2 &&
      matched.at("saturation_actions").as_array().size() == 1 &&
      residual.at("status") == "exact-zero" &&
      residual.at("scope") == "stored-taylor-truncation" &&
      retained.at("status") == "ok" &&
      std::string(retained.at("match").as_string()) == match &&
      retained.at("retained_state") ==
          "exact-lattice-transformation-and-laurent-weights" &&
      stats.at("locals") == 0 && stats.at("matches") == 1 &&
      stats.at("local_matches") == 1 && stats.at("pending_matches") == 0 &&
      stats.at("local_match_capability") ==
          "exact-rational-regular-local-match-v1";

  if (!ok) {
    std::cerr << "local.match: " << json::serialize(matched) << '\n'
              << "match.stats: " << json::serialize(retained) << '\n'
              << "session.stats: " << json::serialize(stats) << '\n';
  }

  (void)request(json::object{{"schema", 2}, {"op", "match.release"},
                             {"session", session}, {"match", match}});
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  std::cout << (ok ? "PASS" : "FAIL")
            << ": persistent exact regular local match smoke\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
