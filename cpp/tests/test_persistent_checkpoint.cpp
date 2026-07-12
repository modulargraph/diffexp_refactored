#include "diffexp2/checkpoint.hpp"
#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include <unistd.h>

namespace json = boost::json;

namespace {

json::object request(const std::string& text) {
  return json::parse(diffexp2::run_recurrence_json(text)).as_object();
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

std::string prepare_two_block_scc(const std::string& session,
                                  const std::string& first,
                                  const std::string& second) {
  const auto prepared = request(std::string(R"json({
    "schema":2,"op":"scc.prepare","session":")json") + session +
    R"json(","key":"checkpoint-two-block-key",
    "identity":"checkpoint-two-block-parent-v1",
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
      "execution":{"mode":"BlockSequentialStrict","work_t_order":6},
      "work_contract":{"work_min":0,"requested_min":0,
        "requested_max":1,"work_complete_max":1,"public_t_order":0,
        "wolfram_coupling_depth":2}},
    "blocks":[
      {"block":0,"vertices":[0],"chart":")json" + first +
    R"json(","principal_identity":"checkpoint-principal-0","regular":true,
       "identity_gauge":true,"identity_v":true,"no_pseudo":true},
      {"block":1,"vertices":[1],"chart":")json" + second +
    R"json(","principal_identity":"checkpoint-principal-1","regular":true,
       "identity_gauge":true,"identity_v":true,"no_pseudo":true}],
    "couplings":[{
      "source_block":0,"target_block":1,
      "source_vertices":[0],"target_vertices":[1],
      "rows":1,"columns":1,"exact_identity":"checkpoint-block-0-to-1",
      "domain":"rational","symbols":[],
      "entries":[{"row":0,"column":0,
        "source_vertex":0,"target_vertex":1,
        "exact_original_entry":"g","exact_theta_entry":"theta-g",
        "multiplier":{"epsilon_shift":-1,"center_pole_order":0,
          "kernels":[["0","1","0","0","0","0","0"],
                     ["0","0","0","0","0","0","0"]],
          "exact_identity":"theta-g","proven_zero":false}}]}]
  })json");
  if (prepared.at("status") != "ok")
    throw std::runtime_error(json::serialize(prepared));
  return std::string(prepared.at("scc").as_string());
}

}  // namespace

int main() {
  diffexp2::ComplexBall::set_precision(384);
  const auto original_ball = diffexp2::ComplexBall::from_strings(
      "[1.234567890123456789 +/- 1e-90]",
      "[-2.500000000000000001 +/- 1e-110]");
  const auto ball_dump =
      diffexp2::checkpoint::dump_complex_ball_exact(original_ball);
  const auto restored_ball =
      diffexp2::checkpoint::load_complex_ball_exact(ball_dump);
  const bool exact_ball_roundtrip =
      acb_equal(original_ball.raw(), restored_ball.raw());

  const std::string checkpoint_path =
      "/tmp/diffexp2_persistent_checkpoint_" +
      std::to_string(static_cast<long long>(::getpid())) + ".de2cp";
  std::filesystem::remove(checkpoint_path);

  const auto created = request(R"json({
    "schema":2,"op":"session.create","domain":"rational",
    "output_digits":30,"chart_capacity":4,"scc_capacity":2,
    "analytic":{"regulators":["rho"],"delta_sign":-1,
      "branch_policy":"prescription-specialized"}
  })json");
  const auto session = std::string(created.at("session").as_string());
  const auto first = prepare_scalar_chart(
      session, "checkpoint-principal-0-key", "checkpoint-principal-0");
  const auto second = prepare_scalar_chart(
      session, "checkpoint-principal-1-key", "checkpoint-principal-1");
  const auto scc = prepare_two_block_scc(session, first, second);

  const auto saved = request(std::string(R"json({
    "schema":2,"op":"checkpoint.save","session":")json") + session +
    R"json(","path":")json" + checkpoint_path +
    R"json(","checkpoint_identity":"native-checkpoint-smoke-v1"})json");
  const auto closed = request(std::string(R"json({
    "schema":2,"op":"session.close","session":")json") + session +
    R"json("})json");
  const auto wrong_identity = request(std::string(R"json({
    "schema":2,"op":"checkpoint.restore","path":")json") +
    checkpoint_path +
    R"json(","expected_identity":"wrong-checkpoint"})json");
  const auto restored = request(std::string(R"json({
    "schema":2,"op":"checkpoint.restore","path":")json") +
    checkpoint_path +
    R"json(","expected_identity":"native-checkpoint-smoke-v1"})json");
  const auto restored_session =
      std::string(restored.at("session").as_string());
  const auto stats = request(std::string(R"json({
    "schema":2,"op":"session.stats","session":")json") +
    restored_session + R"json("})json");
  const auto restored_scc_stats = request(std::string(R"json({
    "schema":2,"op":"scc.stats","session":")json") +
    restored_session + R"json(","scc":")json" + scc + R"json("})json");
  const auto solved = request(std::string(R"json({
    "schema":2,"op":"chart.solve","session":")json") +
    restored_session + R"json(","chart":")json" + first + R"json(","run":{
      "nmax":0,"p":0,"has_initial":true,"adaptive_probe":false,
      "a_target":"0","b_target":"0","a_shift_min":0,
      "a_shifts":["0"],
      "schedule":[[{"case":"R","da":"0","db":"0"}]],
      "initial":["1","0"],"initial_validity":[1],
      "source":null,"return_u":true}}
  )json");

  const bool ok = exact_ball_roundtrip && saved.at("status") == "ok" &&
      saved.at("atomic") == true && saved.at("charts") == 2 &&
      saved.at("sccs") == 1 && std::filesystem::exists(checkpoint_path) &&
      closed.at("status") == "ok" && wrong_identity.at("status") == "error" &&
      restored.at("status") == "ok" && restored_session != session &&
      restored.at("replayed_wolfram_preprocessing") == false &&
      restored.at("charts").as_array().size() == 2 &&
      restored.at("sccs").as_array().size() == 1 &&
      restored.at("analytic_identity").as_object().at("delta_sign") == -1 &&
      stats.at("charts") == 2 && stats.at("scc_charts") == 1 &&
      stats.at("checkpoint_generation") == 1 &&
      stats.at("checkpoint_restore_count") == 1 &&
      stats.at("restored_from_checkpoint_identity") ==
          "native-checkpoint-smoke-v1" &&
      restored_scc_stats.at("identity") ==
          "checkpoint-two-block-parent-v1" &&
      restored_scc_stats.at("native_coupling_depth") == 1 &&
      solved.at("status") == "ok" && solved.at("u").as_array()[0] == "1";
  if (!ok) {
    std::cerr << "saved: " << json::serialize(saved) << '\n'
              << "wrong identity: " << json::serialize(wrong_identity) << '\n'
              << "restored: " << json::serialize(restored) << '\n'
              << "stats: " << json::serialize(stats) << '\n'
              << "scc stats: " << json::serialize(restored_scc_stats) << '\n'
              << "solve: " << json::serialize(solved) << '\n';
  }
  (void)request(std::string(R"json({
    "schema":2,"op":"session.close","session":")json") +
    restored_session + R"json("})json");
  std::filesystem::remove(checkpoint_path);
  std::cout << (ok ? "PASS" : "FAIL")
            << ": exact Arb state and persistent prepared chart/SCC checkpoint restore\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
