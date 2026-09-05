#include "diffexp/kernel/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace json = boost::json;

namespace {

json::object request(const std::string& text) {
  return json::parse(diffexp::kernel::run_recurrence_json(text)).as_object();
}

std::string prepare_scalar_chart(const std::string& session,
                                 const std::string& key,
                                 const std::string& identity) {
  const auto response = request(std::string(R"json({
    "schema":2,"op":"chart.prepare","session":")json") + session +
    R"json(","key":")json" + key + R"json(","identity":")json" +
    identity + R"json(",
    "analytic":{
      "geometry":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,"prescriptions":[]},
      "principal_matrix":[[{"exact":"0","proven_zero":true}]],
      "native_scc_capabilities":{"regular":true,"identity_gauge":true,
        "identity_v":true,"no_pseudo":true}},
    "scc":{"components":[[0]],"structural_edges":[],
      "condensation_edges":[],"topological_order":[0],"coupling_depth":0},
    "problem":{"domain":"rational","d":1,"fb":0,"w":2,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],"val":[null]}],
      "d0_inverse":"1","blocks":[[0]],
      "assembly":{"identity":true,"poly":[],"rat":[],"val":[0]},
      "chop_digits":0}
  })json");
  if (response.at("status") != "ok")
    throw std::runtime_error(json::serialize(response));
  return std::string(response.at("chart").as_string());
}

}  // namespace

int main() {
  const auto created = request(R"json({
    "schema":2,"op":"session.create","domain":"rational",
    "output_digits":30,"scc_capacity":2
  })json");
  const auto session = std::string(created.at("session").as_string());
  const auto first = prepare_scalar_chart(
      session, "principal-0-key", "principal-0");
  const auto second = prepare_scalar_chart(
      session, "principal-1-key", "principal-1");

  const auto prepared = request(std::string(R"json({
    "schema":2,"op":"scc.prepare","session":")json") + session +
    R"json(","key":"two-block-key","identity":"two-block-parent-v1",
    "parent":{
      "dimension":2,
      "exact_system_record":[
        [{"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true}],
        [{"exact":"g","proven_zero":false},
         {"exact":"0","proven_zero":true}]],
      "exact_theta_record":[
        [{"exact":"0","proven_zero":true},
         {"exact":"0","proven_zero":true}],
        [{"exact":"theta-g","proven_zero":false},
         {"exact":"0","proven_zero":true}]],
      "chart":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,"prescriptions":[]},
      "scc":{"components":[[0],[1]],"structural_edges":[[0,1]],
        "condensation_edges":[[0,1]],"topological_order":[0,1],
        "coupling_depth":1},
      "execution":{"mode":"BlockSequentialStrict","work_t_order":10},
      "work_contract":{"work_min":0,"requested_min":0,
        "requested_max":1,"work_complete_max":1,"public_t_order":0,
        "singular_tail_t_order_halo":4,
        "wolfram_coupling_depth":2}},
    "blocks":[
      {"block":0,"vertices":[0],"chart":")json" + first +
    R"json(","principal_identity":"principal-0","regular":true,
       "identity_gauge":true,"identity_v":true,"no_pseudo":true},
      {"block":1,"vertices":[1],"chart":")json" + second +
    R"json(","principal_identity":"principal-1","regular":true,
       "identity_gauge":true,"identity_v":true,"no_pseudo":true}],
    "couplings":[{
      "source_block":0,"target_block":1,
      "source_vertices":[0],"target_vertices":[1],
      "rows":1,"columns":1,"exact_identity":"block-0-to-1",
      "domain":"rational","symbols":[],
      "entries":[{"row":0,"column":0,
        "source_vertex":0,"target_vertex":1,
        "exact_original_entry":"g","exact_theta_entry":"theta-g",
        "multiplier":{"epsilon_shift":-1,"center_pole_order":0,
          "kernels":[
            ["0","1","0","0","0","0","0","0","0","0","0"],
            ["0","0","0","0","0","0","0","0","0","0","0"]],
          "exact_identity":"theta-g","proven_zero":false}}]}],
    "rational_shadow_singular_tail":{
      "schema":"diffexp3-scc-singular-tail-frame-v1",
      "equation":{
        "schema":"diffexp3-physical-cleared-ode-v1",
        "basis":"physical-original-master","theta_coordinate":"local-t",
        "owner_signature_identity":"two-block-parent-v1",
        "payload_identity":"de2-physical-ode-tail-fixture",
        "q":[{"zero":false,"valuation":0,
          "numerator":["1"],"denominator":["1"]}],
        "c":[[{"r":1,"c":0,"v":{"zero":false,"valuation":0,
          "numerator":["1"],"denominator":["1"]}}]]},
      "epsilon_shifts":[1,0]}
  })json");
  if (prepared.at("status") != "ok") {
    std::cerr << "scc.prepare failed: "
              << json::serialize(prepared) << '\n';
    return EXIT_FAILURE;
  }
  const auto scc = std::string(prepared.at("scc").as_string());
  const auto stats = request(std::string(R"json({
    "schema":2,"op":"scc.stats","session":")json") + session +
    R"json(","scc":")json" + scc + R"json("})json");
  const auto first_stats_page = request(std::string(R"json({
    "schema":2,"op":"scc.stats","session":")json") + session +
    R"json(","scc":")json" + scc +
    R"json(","page_offset":0,"page_size":1})json");
  const auto last_stats_page = request(std::string(R"json({
    "schema":2,"op":"scc.stats","session":")json") + session +
    R"json(","scc":")json" + scc +
    R"json(","page_offset":1,"page_size":1})json");
  const auto released = request(std::string(R"json({
    "schema":2,"op":"scc.release","session":")json") + session +
    R"json(","scc":")json" + scc + R"json("})json");
  const auto session_stats = request(std::string(R"json({
    "schema":2,"op":"session.stats","session":")json") + session +
    R"json("})json");

  const bool ok = prepared.at("status") == "ok" &&
      prepared.at("diagnostic_detail") == "bounded-scc-publication" &&
      prepared.if_contains("block_charts") == nullptr &&
      prepared.at("blocks") == 2 &&
      prepared.at("active_coupling_entries") == 1 &&
      prepared.at("execution_implemented") == true &&
      prepared.at("scalar_block_dag_column_execution") == true &&
      prepared.at("min_coupling_shift") == -1 &&
      prepared.at("rational_shadow_singular_tail") ==
          "fuchsian-epsilon-sheared" &&
      prepared.at(
          "rational_shadow_singular_tail_max_epsilon_shift") == 1 &&
      stats.at("status") == "ok" &&
      stats.at("block_charts").as_array().size() == 2 &&
      std::string(stats.at("scc").as_string()) == scc &&
      first_stats_page.at("block_charts").as_array().size() == 1 &&
      first_stats_page.at("block_charts").as_array().front()
          .as_object().at("block") == 0 &&
      first_stats_page.at("diagnostics_page").as_object()
          .at("total") == 2 &&
      first_stats_page.at("diagnostics_page").as_object()
          .at("has_more") == true &&
      first_stats_page.at("diagnostics_page").as_object()
          .at("next_offset") == 1 &&
      last_stats_page.at("block_charts").as_array().size() == 1 &&
      last_stats_page.at("block_charts").as_array().front()
          .as_object().at("block") == 1 &&
      last_stats_page.at("diagnostics_page").as_object()
          .at("has_more") == false &&
      stats.at("frame_base") == 0 && stats.at("frame_width") == 2 &&
      stats.at("public_t_order") == 0 &&
      stats.at("work_t_order") == 10 &&
      stats.at("singular_tail_t_order_halo") == 4 &&
      released.at("status") == "ok" &&
      session_stats.at("scc_charts") == 0;
  if (!ok) {
    std::cerr << "scc.prepare: " << json::serialize(prepared) << '\n'
              << "scc.stats: " << json::serialize(stats) << '\n'
              << "scc.stats first page: "
              << json::serialize(first_stats_page) << '\n'
              << "scc.stats last page: "
              << json::serialize(last_stats_page) << '\n'
              << "scc.release: " << json::serialize(released) << '\n'
              << "session.stats: " << json::serialize(session_stats) << '\n';
  }

  (void)request(std::string(R"json({
    "schema":2,"op":"session.close","session":")json") + session +
    R"json("})json");
  std::cout << (ok ? "PASS" : "FAIL")
            << ": persistent two-block SCC prepare/release smoke\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
