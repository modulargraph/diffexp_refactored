#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdlib>
#include <iostream>
#include <string>

namespace json = boost::json;

namespace {

json::object request(const std::string& text) {
  return json::parse(diffexp2::run_recurrence_json(text)).as_object();
}

}  // namespace

int main() {
  const auto created = request(R"json({
    "schema":2,"op":"session.create","domain":"rational",
    "output_digits":30,"local_capacity":2
  })json");
  const auto session = std::string(created.at("session").as_string());
  const auto prepared = request(std::string(R"json({
    "schema":2,"op":"chart.prepare","session":")json") + session + R"json(",
    "key":"sqrt@0[-2,8]","identity":"persistent-local-smoke-v1",
    "analytic":{"prescription":"+i0"},
    "scc":{"components":[[0]],"structural_edges":[],
      "condensation_edges":[],"topological_order":[0],"coupling_depth":0},
    "problem":{"domain":"rational","d":1,"fb":-2,"w":8,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],"val":[null]}],
      "d0_inverse":"1","blocks":[[0]],
      "assembly":{"identity":false,
        "poly":[{"s":0,"e":[[0,0,"1"]]}],"rat":[],"val":[0]},
      "chop_digits":10}
  })json");
  const auto chart = std::string(prepared.at("chart").as_string());
  const auto local_result = request(std::string(R"json({
    "schema":2,"op":"local.solve","session":")json") + session +
    R"json(","chart":")json" + chart + R"json(",
    "run":{"nmax":0,"p":0,"has_initial":true,"adaptive_probe":false,
      "a_target":"1/2","b_target":"0","a_shift_min":0,
      "a_shifts":["1/2"],
      "schedule":[[{"case":"R","da":"0","db":"0"}]],
      "initial":["0","0","1","0","0","0","0","0"],
      "initial_validity":[5],"source":null,"return_u":false},
    "metadata":{"chart":{"center_exact":"0","scale_exact":"1",
        "radius":"2","infinite_radius":false},
      "tag":{"a":{"domain":"rational","canonical":"1/2"},
        "b":{"domain":"rational","canonical":"0"}},
      "prescriptions":[{"factor_exact":"t","sign":1,
        "multiplicity":1,"leading_coefficient_sign":1}],
      "checkpoint_identity":"sqrt-local-v1"}
  })json");
  const auto local = std::string(local_result.at("local").as_string());

  const auto plus = request(std::string(R"json({
    "schema":2,"op":"local.evaluate","session":")json") + session +
    R"json(","local":")json" + local + R"json(",
    "point":{"exact":"-1/2"},"options":{"tail_estimate":false}
  })json");
  (void)request(std::string(R"json({
    "schema":2,"op":"chart.release","session":")json") + session +
    R"json(","chart":")json" + chart + R"json("})json");
  const auto minus = request(std::string(R"json({
    "schema":2,"op":"local.evaluate","session":")json") + session +
    R"json(","local":")json" + local + R"json(",
    "point":{"exact":"-1/2"},
    "options":{"imaginary_sign":-1,"tail_estimate":false}
  })json");
  const auto stats = request(std::string(R"json({
    "schema":2,"op":"session.stats","session":")json") + session + "\"}");
  const auto local_stats = request(std::string(R"json({
    "schema":2,"op":"local.stats","session":")json") + session +
    R"json(","local":")json" + local + R"json("})json");

  const auto& plus_coefficient = plus.at("value").as_object()
      .at("coefficients").as_array().front().as_array();
  const auto& minus_coefficient = minus.at("value").as_object()
      .at("coefficients").as_array().front().as_array();
  const bool ok = created.at("status") == "ok" &&
      prepared.at("status") == "ok" &&
      local_result.at("status") == "ok" &&
      local_result.at("native_retained") == true &&
      local_result.at("json_coefficients") == 0 &&
      local_result.at("epsilon_min") == 0 &&
      local_result.if_contains("metadata") != nullptr &&
      plus.at("status") == "ok" && plus.at("imaginary_sign") == 1 &&
      minus.at("status") == "ok" && minus.at("imaginary_sign") == -1 &&
      !std::string(plus_coefficient[1].as_string()).starts_with("-") &&
      std::string(minus_coefficient[1].as_string()).starts_with("-") &&
      plus.if_contains("theta") != nullptr &&
      stats.at("charts") == 0 && stats.at("locals") == 1 &&
      stats.at("local_solves") == 1 &&
      stats.at("local_evaluations") == 2 &&
      local_stats.at("status") == "ok" &&
      local_stats.at("evaluations") == 2;

  if (!ok) {
    std::cerr << "local.solve: " << json::serialize(local_result) << '\n'
              << "+i0 evaluate: " << json::serialize(plus) << '\n'
              << "-i0 evaluate: " << json::serialize(minus) << '\n'
              << "stats: " << json::serialize(stats) << '\n';
  }

  (void)request(std::string(R"json({
    "schema":2,"op":"local.release","session":")json") + session +
    R"json(","local":")json" + local + R"json("})json");
  (void)request(std::string(R"json({
    "schema":2,"op":"session.close","session":")json") + session + "\"}");
  std::cout << (ok ? "PASS" : "FAIL")
            << ": persistent native local handle smoke\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
