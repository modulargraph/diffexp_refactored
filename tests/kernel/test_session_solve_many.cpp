#include "diffexp/kernel/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdlib>
#include <iostream>
#include <string>

namespace json = boost::json;

namespace {

json::object request(const std::string& text) {
  return json::parse(diffexp::kernel::run_recurrence_json(text)).as_object();
}

std::string prepare_chart(const std::string& session,
                          const std::string& domain,
                          const std::string& key,
                          const std::string& coefficient,
                          bool symbolic) {
  const auto symbols = symbolic ? R"json(,"symbols":["rho"])json" : "";
  const auto prepared = request(
      std::string(R"json({
        "schema":2,"op":"chart.prepare","session":")json") + session +
      R"json(","key":")json" + key + R"json(","identity":")json" + key +
      R"json(-identity","analytic":{"prescription":"none"},
        "scc":{"components":[[0]],"structural_edges":[],
          "condensation_edges":[],"topological_order":[0],
          "coupling_depth":0},
        "problem":{"domain":")json" + domain + "\"" + symbols + R"json(,
          "d":1,"fb":-2,"w":8,
          "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
          "nhat_lags":[{"poly":[],"rat":[],"val":[null]},
            {"poly":[{"s":0,"e":[[0,0,")json" + coefficient +
      R"json("]]}],"rat":[],"val":[0]}],
          "d0_inverse":"1","blocks":[[0]],"assembly":null,
          "chop_digits":10}})json");
  return std::string(prepared.at("chart").as_string());
}

std::string complete_run(const std::string& initial) {
  return std::string(R"json({
    "nmax":1,"p":0,"has_initial":true,"adaptive_probe":false,
    "a_target":"0","b_target":"0","a_shift_min":0,
    "a_shifts":["0","1"],
    "schedule":[[{"case":"R","da":"0","db":"0"}],
      [{"case":"T","da":"1","db":"0"}]],
    "initial":["0","0",")json") + initial + R"json(","0","0","0","0","0"],
    "initial_validity":[5],"source":null,"return_u":true})json";
}

bool rational_smoke() {
  const auto created = request(R"json({
    "schema":2,"op":"session.create","domain":"rational",
    "output_digits":30
  })json");
  const auto session = std::string(created.at("session").as_string());
  const auto chart_one = prepare_chart(
      session, "rational", "cross-one", "1", false);
  const auto chart_two = prepare_chart(
      session, "rational", "cross-two", "2", false);
  const auto run = complete_run("1");

  // The valid leading job must not run when a later handle fails validation.
  const auto rejected = request(
      std::string(R"json({
        "schema":2,"op":"session.solve_many","session":")json") + session +
      R"json(","threads":8,"jobs":[{"chart":")json" + chart_one +
      R"json(","run":)json" + run +
      R"json(},{"chart":"c:missing","run":)json" + run + R"json(}]})json");
  const auto before = request(
      std::string(R"json({"schema":2,"op":"session.stats","session":")json") +
      session + R"json("})json");

  const auto solved = request(
      std::string(R"json({
        "schema":2,"op":"session.solve_many","session":")json") + session +
      R"json(","threads":64,"jobs":[{"chart":")json" + chart_two +
      R"json(","output_digits":17,"run":)json" + run +
      R"json(},{"chart":")json" + chart_one + R"json(","run":)json" + run +
      R"json(},{"chart":")json" + chart_two +
      R"json(","run":{"nmax":"malformed"}}]})json");
  const auto after = request(
      std::string(R"json({"schema":2,"op":"session.stats","session":")json") +
      session + R"json("})json");
  const auto counters = request(
      std::string(R"json({"schema":2,"op":"session.counters","session":")json") +
      session + R"json("})json");
  const auto rejected_counters = request(
      std::string(R"json({"schema":2,"op":"session.counters","session":")json") +
      session + R"json(","detail":true})json");

  const auto& results = solved.at("results").as_array();
  const auto& first_u = results[0].as_object().at("u").as_array();
  const auto& second_u = results[1].as_object().at("u").as_array();
  const bool ok = rejected.at("status") == "error" &&
      rejected.at("id") == "CPP" && before.at("runs") == 0 &&
      solved.at("status") == "ok" && solved.at("attempted") == 3 &&
      solved.at("succeeded") == 2 && solved.at("failed") == 1 &&
      solved.at("worker_threads") == 3 &&
      solved.at("symbolic_serialized") == false &&
      std::string(results[0].as_object().at("chart").as_string()) ==
          chart_two &&
      std::string(results[1].as_object().at("chart").as_string()) ==
          chart_one &&
      std::string(results[2].as_object().at("chart").as_string()) ==
          chart_two &&
      first_u[18] == "2" && second_u[18] == "1" &&
      results[2].as_object().at("status") == "error" &&
      results[2].as_object().at("id") == "CPP" &&
      after.at("runs") == 2 && after.at("static_tensor_copies") == 0 &&
      counters.at("status") == "ok" &&
      counters.at("scope") == "fixed-session-counters" &&
      std::string(counters.at("session").as_string()) == session &&
      counters.at("charts") == 2 &&
      counters.at("locals") == 0 && counters.at("local_solves") == 0 &&
      counters.if_contains("chart_stats") == nullptr &&
      counters.if_contains("local_stats") == nullptr &&
      counters.if_contains("retained_derivation") == nullptr &&
      json::serialize(counters).size() < 4096 &&
      rejected_counters.at("status") == "error" &&
      rejected_counters.at("id") == "CPP";

  if (!ok) {
    std::cerr << "rejected: " << json::serialize(rejected) << '\n'
              << "stats before: " << json::serialize(before) << '\n'
              << "solved: " << json::serialize(solved) << '\n'
              << "stats after: " << json::serialize(after) << '\n'
              << "counters: " << json::serialize(counters) << '\n'
              << "rejected counters: " << json::serialize(rejected_counters)
              << '\n';
  }
  (void)request(
      std::string(R"json({"schema":2,"op":"session.close","session":")json") +
      session + R"json("})json");
  return ok;
}

bool symbolic_smoke() {
  const auto created = request(R"json({
    "schema":2,"op":"session.create","domain":"symbolic",
    "symbols":["rho"],"output_digits":30
  })json");
  const auto session = std::string(created.at("session").as_string());
  const auto chart_one = prepare_chart(
      session, "symbolic", "symbolic-one", "1+rho", true);
  const auto chart_two = prepare_chart(
      session, "symbolic", "symbolic-two", "2-rho", true);
  const auto run = complete_run("1");
  const auto solved = request(
      std::string(R"json({
        "schema":2,"op":"session.solve_many","session":")json") + session +
      R"json(","threads":8,"jobs":[{"chart":")json" + chart_one +
      R"json(","run":)json" + run + R"json(},{"chart":")json" + chart_two +
      R"json(","run":)json" + run + R"json(}]})json");
  const auto stats = request(
      std::string(R"json({"schema":2,"op":"session.stats","session":")json") +
      session + R"json("})json");
  const auto& results = solved.at("results").as_array();
  const bool ok = solved.at("status") == "ok" &&
      solved.at("attempted") == 2 && solved.at("succeeded") == 2 &&
      solved.at("worker_threads") == 1 &&
      solved.at("symbolic_serialized") == true && stats.at("runs") == 2 &&
      std::string(results[0].as_object().at("chart").as_string()) ==
          chart_one &&
      std::string(results[1].as_object().at("chart").as_string()) ==
          chart_two;
  if (!ok) {
    std::cerr << "symbolic solved: " << json::serialize(solved) << '\n'
              << "symbolic stats: " << json::serialize(stats) << '\n';
  }
  (void)request(
      std::string(R"json({"schema":2,"op":"session.close","session":")json") +
      session + R"json("})json");
  return ok;
}

}  // namespace

int main() {
  const bool rational = rational_smoke();
  const bool symbolic = symbolic_smoke();
  const bool ok = rational && symbolic;
  std::cout << (ok ? "PASS" : "FAIL")
            << ": persistent cross-chart solve_many smoke\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
