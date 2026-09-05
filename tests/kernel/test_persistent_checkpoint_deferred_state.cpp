#include "diffexp/kernel/checkpoint.hpp"
#include "diffexp/kernel/json_codec.hpp"

#include <boost/json.hpp>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <unistd.h>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp::kernel::run_recurrence_json(json::serialize(value)))
      .as_object();
}

std::string prepare_chart(const std::string& session,
                          const std::string& key,
                          const std::string& identity,
                          const std::string& center) {
  auto value = json::parse(R"json({
    "schema":2,"op":"chart.prepare","session":"placeholder",
    "key":"placeholder","identity":"placeholder",
    "analytic":{
      "geometry":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,
        "prescriptions":[{"factor_exact":"t","sign":-1,
          "multiplicity":1,"leading_coefficient_sign":1}]},
      "principal_matrix":[[{"exact":"0","proven_zero":true}]],
      "native_scc_capabilities":{"regular":true,"identity_gauge":true,
        "identity_v":true,"no_pseudo":true}},
    "scc":{"components":[[0]],"structural_edges":[],
      "condensation_edges":[],"topological_order":[0],
      "coupling_depth":0},
    "problem":{"domain":"rational","d":1,"fb":0,"w":3,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],"val":[null]}],
      "d0_inverse":"1","blocks":[[0]],
      "assembly":{"identity":true,"poly":[],"rat":[],"val":[0]},
      "chop_digits":0}
  })json").as_object();
  value["session"] = session;
  value["key"] = key;
  value["identity"] = identity;
  value.at("analytic").as_object().at("geometry").as_object()
      ["center_exact"] = center;
  const auto response = request(std::move(value));
  if (response.at("status") != "ok")
    throw std::runtime_error("chart.prepare: " + json::serialize(response));
  return std::string(response.at("chart").as_string());
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
      {"prescriptions", json::array{
          json::object{{"factor_exact", "t"}, {"sign", -1},
                       {"multiplicity", 1},
                       {"leading_coefficient_sign", 1}}}},
      {"checkpoint_identity", checkpoint}};
}

std::string solve_local(const std::string& session,
                        const std::string& chart,
                        const std::string& checkpoint) {
  json::array schedule_row;
  schedule_row.push_back(
      json::object{{"case", "R"}, {"da", "0"}, {"db", "0"}});
  json::array schedule;
  schedule.push_back(std::move(schedule_row));
  auto response = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart},
      {"run", json::object{
          {"nmax", 0}, {"p", 0}, {"has_initial", true},
          {"adaptive_probe", false}, {"a_target", "0"},
          {"b_target", "0"}, {"a_shift_min", 0},
          {"a_shifts", json::array{"0"}},
          {"schedule", std::move(schedule)},
          {"initial", json::array{"3", "0", "0"}},
          {"initial_validity", json::array{2}},
          {"source", nullptr}, {"return_u", false}}},
      {"metadata", metadata(checkpoint)}});
  if (response.at("status") != "ok")
    throw std::runtime_error("local.solve: " + json::serialize(response));
  return std::string(response.at("local").as_string());
}

json::object exact_match(const std::string& session,
                         const std::string& chart,
                         const std::string& basis,
                         const std::string& basis_checkpoint,
                         const std::string& incoming,
                         const std::string& incoming_checkpoint,
                         const std::string& match_checkpoint) {
  return request(json::object{
      {"schema", 2}, {"op", "local.match"}, {"session", session},
      {"basis", json::array{basis}}, {"incoming", incoming},
      {"basis_chart", chart}, {"incoming_chart", chart},
      {"basis_point", json::object{{"exact", "1/2"}}},
      {"incoming_point", json::object{{"exact", "1/2"}}},
      {"epsilon", json::object{{"min", 0}, {"max", 2},
                                {"required_complete_max", 1}}},
      {"basis_checkpoint_identities", json::array{basis_checkpoint}},
      {"incoming_checkpoint_identity", incoming_checkpoint},
      {"checkpoint_identity", match_checkpoint}});
}

json::object topology() {
  return json::object{
      {"singular_points", json::array{}},
      {"boundary_points", json::array{}},
      {"complex_projections", json::array{}},
      {"branch_sheets", json::array{
          json::object{{"factor_exact", "t"}, {"sign", -1}}}}};
}

json::object arm(const std::string& endpoint,
                 const std::string& anchor,
                 const std::string& endpoint_chart) {
  return json::object{{"from_exact", "0"}, {"to_exact", endpoint},
                      {"charts", json::array{anchor, endpoint_chart}},
                      {"topology", topology()}};
}

json::object tile_plan(const std::string& session,
                       const std::string& anchor,
                       const std::string& lower,
                       const std::string& upper,
                       const std::string& checkpoint) {
  return request(json::object{
      {"schema", 2}, {"op", "tile.plan"}, {"session", session},
      {"checkpoint_identity", checkpoint}, {"division_order", 3},
      {"lower", arm("-2/3", anchor, lower)},
      {"upper", arm("2/3", anchor, upper)}});
}

json::object integrate(const std::string& session,
                       const std::string& plan,
                       const std::string& plan_checkpoint,
                       const std::string& local,
                       const std::string& local_checkpoint,
                       const std::string& checkpoint) {
  return request(json::object{
      {"schema", 2}, {"op", "integration.line"}, {"session", session},
      {"tile_plan", plan}, {"local", local}, {"arm", "upper"},
      {"tile", 0},
      {"epsilon", json::object{{"min", 0}, {"max", 0}}},
      {"source_checkpoint_identity", local_checkpoint},
      {"tile_plan_checkpoint_identity", plan_checkpoint},
      {"checkpoint_identity", checkpoint}});
}

double exported_value(const json::object& response) {
  return std::stod(std::string(response.at("value").as_object()
      .at("coefficients").as_array().front().as_array().front().as_string()));
}

bool corrupt_last_payload_byte(const std::filesystem::path& source,
                               const std::filesystem::path& target) {
  std::error_code error;
  std::filesystem::copy_file(
      source, target, std::filesystem::copy_options::overwrite_existing,
      error);
  if (error) return false;
  std::fstream file(target, std::ios::in | std::ios::out | std::ios::binary);
  if (!file) return false;
  file.seekg(-1, std::ios::end);
  char byte = 0;
  file.read(&byte, 1);
  if (!file) return false;
  byte ^= 0x5a;
  file.seekp(-1, std::ios::end);
  file.write(&byte, 1);
  file.flush();
  return static_cast<bool>(file);
}

}  // namespace

int main() {
  const auto base = std::filesystem::temp_directory_path() /
      ("diffexp3-checkpoint-deferred-" + std::to_string(::getpid()));
  const auto path = base.string() + ".de2cp";
  const auto corrupt_path = base.string() + "-corrupt.de2cp";
  bool ok = false;
  std::string stage = "create";
  std::string original_session;
  std::string restored_session;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "rational"}, {"precision_bits", 256},
        {"output_digits", 40}, {"chart_capacity", 8},
        {"local_capacity", 4}, {"match_capacity", 4}});
    original_session = std::string(created.at("session").as_string());
    const auto anchor = prepare_chart(
        original_session, "checkpoint-anchor", "checkpoint-anchor-v1", "0");
    const auto lower = prepare_chart(
        original_session, "checkpoint-lower", "checkpoint-lower-v1", "-2/3");
    const auto upper = prepare_chart(
        original_session, "checkpoint-upper", "checkpoint-upper-v1", "2/3");
    const auto local = solve_local(
        original_session, anchor, "deferred-basis-v1");
    const auto incoming = solve_local(
        original_session, anchor, "deferred-incoming-v1");
    stage = "first-match";
    const auto matched = exact_match(
        original_session, anchor, local, "deferred-basis-v1", incoming,
        "deferred-incoming-v1", "deferred-match-v1");
    if (matched.at("status") != "ok")
      throw std::runtime_error(json::serialize(matched));
    const auto match = std::string(matched.at("match").as_string());
    stage = "first-plan";
    const auto planned = tile_plan(
        original_session, anchor, lower, upper, "deferred-plan-v1");
    if (planned.at("status") != "ok")
      throw std::runtime_error(json::serialize(planned));
    const auto plan = std::string(planned.at("tile_plan").as_string());
    (void)request(json::object{
        {"schema", 2}, {"op", "tile.match_interval"},
        {"session", original_session}, {"tile_plan", plan},
        {"arm", "upper"}, {"match", 0}});
    stage = "first-line";
    const auto integrated = integrate(
        original_session, plan, "deferred-plan-v1", local,
        "deferred-basis-v1", "deferred-line-v1");
    if (integrated.at("status") != "ok")
      throw std::runtime_error(json::serialize(integrated));
    const auto line = std::string(integrated.at("line").as_string());
    stage = "first-export";
    const auto initial_export = request(json::object{
        {"schema", 2}, {"op", "integration.export"},
        {"session", original_session}, {"line", line},
        {"checkpoint_identity", "deferred-line-v1"},
        {"output_digits", 40}});

    // The line owns the plan/local and the plan owns both released endpoint
    // charts.  The exact match independently owns the same local snapshot.
    for (const auto& source : {local, incoming})
      (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                                 {"session", original_session},
                                 {"local", source}});
    (void)request(json::object{{"schema", 2}, {"op", "tile.release"},
                               {"session", original_session},
                               {"tile_plan", plan}});
    for (const auto& chart : {lower, upper})
      (void)request(json::object{{"schema", 2}, {"op", "chart.release"},
                                 {"session", original_session},
                                 {"chart", chart}});

    stage = "save";
    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", original_session}, {"path", path},
        {"checkpoint_identity", "deferred-session-v2"}});
    stage = "payload";
    const auto payload = json::parse(
        diffexp::kernel::checkpoint::read(path).payload_json).as_object();
    const auto& saved_session = payload.at("session").as_object();
    const auto& visibility =
        saved_session.at("registry_visibility").as_object();
    const auto& line_record = payload.at("retained_line_results")
        .as_array().front().as_object();
    const auto& ball = line_record.at("result").as_object()
        .at("value").as_object().at("coefficients").as_array()
        .front().as_object();

    stage = "close-original";
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", original_session}});
    original_session.clear();
    stage = "restore";
    const auto restored = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"}, {"path", path},
        {"expected_identity", "deferred-session-v2"}});
    restored_session = std::string(restored.at("session").as_string());
    stage = "inspect-restored";
    const auto retained_match = request(json::object{
        {"schema", 2}, {"op", "match.stats"},
        {"session", restored_session}, {"match", match}});
    const auto retained_line = request(json::object{
        {"schema", 2}, {"op", "integration.stats"},
        {"session", restored_session}, {"line", line}});
    const auto restored_export = request(json::object{
        {"schema", 2}, {"op", "integration.export"},
        {"session", restored_session}, {"line", line},
        {"checkpoint_identity", "deferred-line-v1"},
        {"output_digits", 40}});
    const auto hidden_local = request(json::object{
        {"schema", 2}, {"op", "local.stats"},
        {"session", restored_session}, {"local", local}});
    const auto hidden_plan = request(json::object{
        {"schema", 2}, {"op", "tile.stats"},
        {"session", restored_session}, {"tile_plan", plan}});

    // Creation continues from the exact saved counters even though c:2,
    // c:3, l:1 and tile:1 are dependency-only after restoration.
    stage = "continue-counters";
    const auto lower2 = prepare_chart(
        restored_session, "checkpoint-lower-2", "checkpoint-lower-v2", "-2/3");
    const auto upper2 = prepare_chart(
        restored_session, "checkpoint-upper-2", "checkpoint-upper-v2", "2/3");
    const auto local2 = solve_local(
        restored_session, anchor, "deferred-basis-v2");
    const auto incoming2 = solve_local(
        restored_session, anchor, "deferred-incoming-v2");
    const auto matched2 = exact_match(
        restored_session, anchor, local2, "deferred-basis-v2", incoming2,
        "deferred-incoming-v2", "deferred-match-v2");
    const auto planned2 = tile_plan(
        restored_session, anchor, lower2, upper2, "deferred-plan-v2");
    const auto line2 = integrate(
        restored_session, std::string(planned2.at("tile_plan").as_string()),
        "deferred-plan-v2", local2, "deferred-basis-v2",
        "deferred-line-v2");
    const auto stats = request(json::object{
        {"schema", 2}, {"op", "session.stats"},
        {"session", restored_session}});

    stage = "corruption";
    const bool corrupted = corrupt_last_payload_byte(path, corrupt_path);
    const auto corruption = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", corrupt_path},
        {"expected_identity", "deferred-session-v2"}});

    stage = "assertions";
    ok = saved.at("status") == "ok" && saved.at("exact_matches") == 1 &&
        saved.at("tile_plans") == 0 && saved.at("line_results") == 1 &&
        saved.at("deferred_handle_kinds").as_array() ==
            json::array{"symbolic-local"} &&
        payload.at("retained_locals").as_array().size() == 2 &&
        payload.at("retained_locals").as_array().front().as_object()
            .at("runtime_stats").as_object().at("line_integrations") == 1 &&
        payload.at("retained_exact_matches").as_array().size() == 1 &&
        payload.at("retained_tile_plans").as_array().size() == 1 &&
        payload.at("retained_line_results").as_array().size() == 1 &&
        payload.at("retained_tile_plans").as_array().front().as_object()
            .at("runtime_stats").as_object().at("integrations") == 1 &&
        line_record.at("runtime_stats").as_object().at("exports") == 1 &&
        visibility.at("charts").as_array() == json::array{anchor} &&
        visibility.at("locals").as_array().empty() &&
        visibility.at("tile_plans").as_array().empty() &&
        ball.if_contains("real") != nullptr &&
        ball.if_contains("imaginary") != nullptr &&
        restored.at("status") == "ok" &&
        restored.at("charts").as_array().size() == 1 &&
        restored.at("locals").as_array().empty() &&
        restored.at("exact_matches").as_array().size() == 1 &&
        restored.at("tile_plans").as_array().empty() &&
        restored.at("line_results").as_array().size() == 1 &&
        retained_match.at("status") == "ok" &&
        std::string(retained_match.at("match").as_string()) == match &&
        retained_line.at("status") == "ok" &&
        retained_line.at("exports") == 1 &&
        retained_line.at("effective_rim") == -1 &&
        std::abs(exported_value(initial_export) - 1.0) < 1e-30 &&
        std::abs(exported_value(restored_export) - 1.0) < 1e-30 &&
        hidden_local.at("status") == "error" &&
        hidden_plan.at("status") == "error" &&
        lower2 == "c:4" && upper2 == "c:5" && local2 == "l:3" &&
        incoming2 == "l:4" &&
        matched2.at("match") == "m:2" &&
        planned2.at("tile_plan") == "tile:2" &&
        line2.at("line") == "line:2" &&
        stats.at("local_matches") == 2 &&
        stats.at("tile_plans_created") == 2 &&
        stats.at("line_integrations") == 2 &&
        stats.at("line_exports") == 2 && corrupted &&
        corruption.at("status") == "error" &&
        std::string(corruption.at("detail").as_string()).find("checksum") !=
            std::string::npos;
    if (!ok) {
      std::cerr << "saved: " << json::serialize(saved) << '\n'
                << "payload: " << json::serialize(payload) << '\n'
                << "restored: " << json::serialize(restored) << '\n'
                << "match: " << json::serialize(retained_match) << '\n'
                << "line: " << json::serialize(retained_line) << '\n'
                << "stats: " << json::serialize(stats) << '\n'
                << "corruption: " << json::serialize(corruption) << '\n';
    }
  } catch (const std::exception& error) {
    std::cerr << "exception at " << stage << ": " << error.what() << '\n';
  }
  if (!original_session.empty())
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", original_session}});
  if (!restored_session.empty())
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", restored_session}});
  std::error_code ignored;
  std::filesystem::remove(path, ignored);
  std::filesystem::remove(corrupt_path, ignored);
  std::cout << (ok ? "PASS" : "FAIL")
            << ": checkpoint deferred retained state roundtrip and CRC\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
