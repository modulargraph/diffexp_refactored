#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
      .as_object();
}

std::string prepare_regular_chart(const std::string& session) {
  auto response = request(json::parse(std::string(R"json({
    "schema":2,"op":"chart.prepare","session":")json") + session +
      R"json(","key":"acb-local-match-chart",
    "identity":"acb-local-match-chart-v1",
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
    "problem":{"domain":"acb","precision_bits":256,
      "d":2,"fb":0,"w":6,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],
        "val":[null,null,null,null]}],
      "d0_inverse":"1","blocks":[[0],[1]],
      "assembly":{"identity":true,"poly":[],"rat":[],
        "val":[0,null,null,0]},"chop_digits":0}
  })json").as_object());
  if (response.at("status") != "ok")
    throw std::runtime_error("prepare: " + json::serialize(response));
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
      {"initial_validity", json::array{5, 5}},
      {"source", nullptr}, {"return_u", false}};
}

json::object metadata(const std::string& checkpoint) {
  return json::object{
      {"chart", json::object{{"center_exact", "0"},
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

std::string solve_local(const std::string& session, const std::string& chart,
                        std::vector<std::string> initial,
                        const std::string& checkpoint) {
  auto response = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart}, {"run", regular_run(std::move(initial))},
      {"metadata", metadata(checkpoint)}});
  if (response.at("status") != "ok")
    throw std::runtime_error("solve " + checkpoint + ": " +
                             json::serialize(response));
  return std::string(response.at("local").as_string());
}

json::object exact_frame(std::vector<std::string> coefficients) {
  json::array encoded;
  for (auto& coefficient : coefficients)
    encoded.emplace_back(std::move(coefficient));
  return json::object{{"min", 0}, {"max", 5},
                      {"coefficients", std::move(encoded)}};
}

json::object exact_lattice() {
  const std::vector<std::string> one{"1", "0", "0", "0", "0", "0"};
  const std::vector<std::string> zero{"0", "0", "0", "0", "0", "0"};
  const std::vector<std::string> epsilon{"0", "1", "0", "0", "0", "0"};
  json::array first_row{exact_frame(one), exact_frame(one)};
  json::array second_row{exact_frame(zero), exact_frame(epsilon)};
  return json::object{
      {"schema", "diffexp2-exact-evaluated-epsilon-lattice-v1"},
      {"identity", "nontrivial-eps-pole-lattice-v1"},
      {"evaluated_basis",
       json::array{std::move(first_row), std::move(second_row)}}};
}

json::object match_request(const std::string& session,
                           const std::string& chart,
                           const std::string& first,
                           const std::string& second,
                           const std::string& incoming,
                           const std::string& first_checkpoint,
                           const std::string& checkpoint) {
  return json::object{
      {"schema", 2}, {"op", "local.match_acb"}, {"session", session},
      {"basis", json::array{first, second}}, {"incoming", incoming},
      {"basis_chart", chart}, {"incoming_chart", chart},
      {"basis_point", json::object{{"exact", "1/2"}}},
      {"incoming_point", json::object{{"exact", "1/2"}}},
      {"epsilon", json::object{{"min", 0}, {"max", 5},
                                {"required_complete_max", 3}}},
      {"basis_checkpoint_identities",
       json::array{first_checkpoint, "basis-1"}},
      {"incoming_checkpoint_identity", "incoming"},
      {"checkpoint_identity", checkpoint},
      {"exact_lattice", exact_lattice()},
      {"refinement", json::object{{"relative_tolerance", "1e-50"},
                                    {"max_steps", 2}}}};
}

std::vector<std::string> column(std::string first_component,
                                std::string second_epsilon = "0") {
  return {std::move(first_component), "0", "0", "0", "0", "0",
          "0", std::move(second_epsilon), "0", "0", "0", "0"};
}

}  // namespace

int main() {
  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", "acb"},
      {"precision_bits", 256}, {"output_digits", 30},
      {"local_capacity", 4}, {"match_capacity", 2}});
  const auto session = std::string(created.at("session").as_string());
  const auto chart = prepare_regular_chart(session);

  const auto first = solve_local(session, chart, column("1"), "basis-0");
  const auto second = solve_local(session, chart, column("1", "1"),
                                  "basis-1");
  const auto incoming = solve_local(session, chart, column("3", "2"),
                                    "incoming");
  const auto ambiguous = solve_local(
      session, chart, column("[0 +/- 1]"), "basis-ambiguous");

  const auto matched = request(match_request(
      session, chart, first, second, incoming, "basis-0",
      "acb-match-state-v1"));
  const auto match = std::string(matched.at("match").as_string());
  const auto retained = request(json::object{
      {"schema", 2}, {"op", "match.stats"}, {"session", session},
      {"match", match}});

  // The exact witness still supplies T, but it can never license a numeric
  // pivot whose enclosure overlaps zero.
  const auto ambiguous_result = request(match_request(
      session, chart, ambiguous, second, incoming, "basis-ambiguous",
      "acb-ambiguous-match-v1"));

  for (const auto& local : {first, second, incoming, ambiguous})
    (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                               {"session", session}, {"local", local}});

  const auto stamp = std::chrono::steady_clock::now()
                         .time_since_epoch().count();
  const auto checkpoint_path =
      std::filesystem::temp_directory_path() /
      ("diffexp2-live-acb-match-" + std::to_string(stamp) + ".chk");
  std::error_code ignored;
  std::filesystem::remove(checkpoint_path, ignored);
  const auto live_checkpoint = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
      {"path", checkpoint_path.string()},
      {"checkpoint_identity", "live-acb-match-checkpoint-v1"}});

  const auto& exact = matched.at("exact_lattice").as_object();
  const auto& residual = matched.at("residual").as_object();
  const auto& refinement = matched.at("refinement").as_object();
  const auto& retained_basis = matched.at("basis").as_array();
  const auto& retained_first = retained_basis.front().as_object();
  const auto encoded_match = json::serialize(matched);
  const bool ambiguous_rejected =
      ambiguous_result.at("status") == "error" &&
      std::string(ambiguous_result.at("detail").as_string())
          .find("overlaps zero") != std::string::npos;
  const bool checkpoint_rejected =
      live_checkpoint.at("status") == "error" &&
      std::string(live_checkpoint.at("detail").as_string())
          .find("match") != std::string::npos;

  const auto released = request(json::object{
      {"schema", 2}, {"op", "match.release"}, {"session", session},
      {"match", match}});
  const auto released_again = request(json::object{
      {"schema", 2}, {"op", "match.release"}, {"session", session},
      {"match", match}});
  const auto saved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
      {"path", checkpoint_path.string()},
      {"checkpoint_identity", "released-acb-match-checkpoint-v1"}});
  const auto stats = request(json::object{
      {"schema", 2}, {"op", "session.stats"}, {"session", session}});

  const bool ok =
      created.at("status") == "ok" &&
      created.at("local_match_capability") == "unsupported" &&
      created.at("acb_local_match_capability") ==
          "exact-lattice-guided-acb-local-match-v1" &&
      matched.at("status") == "ok" && matched.at("native_retained") == true &&
      matched.at("json_coefficients") == 0 && matched.at("dimension") == 2 &&
      matched.at("capability") ==
          "exact-lattice-guided-acb-local-match-v1" &&
      std::string(retained_first.at("chart").as_string()) == chart &&
      retained_first.at("checkpoint_identity") == "basis-0" &&
      retained_first.at("requested_imaginary_sign").is_null() &&
      retained_first.at("effective_imaginary_sign").is_null() &&
      retained_first.if_contains("analytic_metadata") != nullptr &&
      exact.at("identity") == "nontrivial-eps-pole-lattice-v1" &&
      exact.at("canonical_witness_retained") == true &&
      exact.at("transformation_min_power") == -1 &&
      exact.at("normalized_determinant_valuation") == 1 &&
      exact.at("saturation_actions").as_array().size() == 1 &&
      refinement.at("factorizations") == 1 &&
      residual.at("verdict") == "pass" &&
      residual.at("complete_through_required") == true &&
      residual.at("history").as_array().size() >= 1 &&
      encoded_match.find("evaluated_basis") == std::string::npos &&
      encoded_match.find("\"coefficients\"") == std::string::npos &&
      retained.at("status") == "ok" &&
      std::string(retained.at("match").as_string()) == match &&
      retained.at("provenance_identity") == matched.at("provenance_identity") &&
      ambiguous_rejected && checkpoint_rejected &&
      released.at("status") == "ok" &&
      released_again.at("status") == "error" &&
      saved.at("status") == "ok" && stats.at("locals") == 0 &&
      stats.at("matches") == 0 && stats.at("pending_matches") == 0 &&
      stats.at("local_matches") == 1 &&
      stats.at("acb_local_match_capability") ==
          "exact-lattice-guided-acb-local-match-v1";

  if (!ok) {
    std::cerr << "local.match_acb: " << json::serialize(matched) << '\n'
              << "match.stats: " << json::serialize(retained) << '\n'
              << "ambiguous: " << json::serialize(ambiguous_result) << '\n'
              << "live checkpoint: " << json::serialize(live_checkpoint)
              << '\n' << "released checkpoint: " << json::serialize(saved)
              << '\n' << "session.stats: " << json::serialize(stats) << '\n';
  }

  std::filesystem::remove(checkpoint_path, ignored);
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  std::cout << (ok ? "PASS" : "FAIL")
            << ": persistent exact-lattice-guided Acb local match smoke\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
