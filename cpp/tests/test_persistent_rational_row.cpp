#include "diffexp2/checkpoint.hpp"
#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#include <unistd.h>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
      .as_object();
}

json::array prescriptions() {
  return json::array{json::object{
      {"factor_exact", "x"}, {"sign", -1}, {"multiplicity", 1},
      {"leading_coefficient_sign", 1}}};
}

std::string prepare_chart(const std::string& session,
                          const std::string& domain,
                          const std::string& suffix,
                          const std::string& center) {
  auto value = json::parse(R"json({
    "schema":2,"op":"chart.prepare","session":"placeholder",
    "key":"placeholder","identity":"placeholder",
    "analytic":{
      "geometry":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,
        "prescriptions":[]},
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
    "problem":{"domain":"rational","precision_bits":256,
      "d":2,"fb":0,"w":3,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],
        "val":[null,null,null,null]}],
      "d0_inverse":"1","blocks":[[0],[1]],
      "assembly":{"identity":true,"poly":[],"rat":[],
        "val":[0,null,null,0]},"chop_digits":0}
  })json").as_object();
  value["session"] = session;
  value["key"] = "rational-row-" + domain + "-" + suffix;
  value["identity"] = "rational-row-chart-" + domain + "-" + suffix;
  value.at("problem").as_object()["domain"] = domain;
  auto& geometry = value.at("analytic").as_object()
      .at("geometry").as_object();
  geometry["center_exact"] = center;
  geometry["prescriptions"] = prescriptions();
  const auto response = request(std::move(value));
  if (response.at("status") != "ok")
    throw std::runtime_error("prepare: " + json::serialize(response));
  return std::string(response.at("chart").as_string());
}

json::object run() {
  json::array schedule_row{
      json::object{{"case", "R"}, {"da", "0"}, {"db", "0"}},
      json::object{{"case", "R"}, {"da", "0"}, {"db", "0"}}};
  json::array schedule;
  schedule.push_back(std::move(schedule_row));
  return json::object{
      {"nmax", 0}, {"p", 0}, {"has_initial", true},
      {"adaptive_probe", false}, {"a_target", "1/2"},
      {"b_target", "0"}, {"a_shift_min", 0},
      {"a_shifts", json::array{"1/2"}},
      {"schedule", std::move(schedule)},
      {"initial", json::array{"1", "0", "0", "2", "0", "0"}},
      {"initial_validity", json::array{2, 2}},
      {"source", nullptr}, {"return_u", false}};
}

json::object metadata(const std::string& checkpoint) {
  return json::object{
      {"chart", json::object{{"center_exact", "0"},
                              {"scale_exact", "1"}, {"radius", "2"},
                              {"infinite_radius", false}}},
      {"tag", json::object{
          {"a", json::object{{"domain", "rational"},
                              {"canonical", "1/2"}}},
          {"b", json::object{{"domain", "rational"},
                              {"canonical", "0"}}}}},
      {"prescriptions", prescriptions()},
      {"checkpoint_identity", checkpoint}};
}

json::object solve_local(const std::string& session,
                         const std::string& chart,
                         const std::string& checkpoint) {
  const auto response = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart}, {"run", run()},
      {"metadata", metadata(checkpoint)}});
  if (response.at("status") != "ok")
    throw std::runtime_error("solve: " + json::serialize(response));
  return response;
}

json::array zero_kernel() {
  return json::array{"0"};
}

json::object multiplier(std::int32_t shift, std::uint32_t pole,
                        json::array leading,
                        std::string identity) {
  return json::object{
      {"epsilon_shift", shift}, {"center_pole_order", pole},
      {"kernels", json::array{
           std::move(leading), zero_kernel(), zero_kernel()}},
      {"exact_identity", std::move(identity)},
      {"proven_zero", false}};
}

json::object rational_row() {
  return json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", 2}, {"exact_identity", "[1/t,3eps]"},
      {"entries", json::array{
           json::object{{"column", 0},
                        {"multiplier", multiplier(
                             0, 1, json::array{"1"},
                             "1/t")}},
           json::object{{"column", 1},
                        {"multiplier", multiplier(
                             1, 0, json::array{"3"},
                             "3eps")}}}}};
}

json::object zero_row() {
  return json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", 2}, {"exact_identity", "[0,0]"},
      {"entries", json::array{}}};
}

json::object scalar_identity_row() {
  return json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", 1}, {"exact_identity", "[1]"},
      {"entries", json::array{json::object{
           {"column", 0},
           {"multiplier", multiplier(
                0, 0, json::array{"1"}, "scalar-unit")}}}}};
}

json::object apply_row(const std::string& session,
                       const std::string& local,
                       json::object row,
                       const std::string& source_checkpoint,
                       const std::string& checkpoint) {
  return request(json::object{
      {"schema", 2}, {"op", "local.apply_rational_row"},
      {"session", session}, {"local", local}, {"row", std::move(row)},
      {"source_checkpoint_identity", source_checkpoint},
      {"checkpoint_identity", checkpoint}});
}

json::object topology() {
  return json::object{
      {"singular_points", json::array{}},
      {"boundary_points", json::array{}},
      {"complex_projections", json::array{}},
      {"branch_sheets", json::array{
           json::object{{"factor_exact", "x"}, {"sign", -1}}}}};
}

json::object arm(const std::string& endpoint,
                 const std::string& anchor,
                 const std::string& endpoint_chart) {
  return json::object{
      {"from_exact", "0"}, {"to_exact", endpoint},
      {"charts", json::array{anchor, endpoint_chart}},
      {"topology", topology()}};
}

double coefficient_midpoint(const json::object& response,
                            std::size_t index) {
  const auto& encoded = response.at("value").as_object()
      .at("coefficients").as_array().at(index).as_array();
  return std::stod(std::string(encoded.at(0).as_string()));
}

const json::object& retained_local_record(const json::object& payload,
                                          const std::string& handle) {
  for (const auto& raw : payload.at("retained_locals").as_array()) {
    const auto& item = raw.as_object();
    if (std::string(item.at("handle").as_string()) == handle)
      return item;
  }
  throw std::runtime_error("checkpoint omitted retained local " + handle);
}

json::object& retained_local_record(json::object& payload,
                                    const std::string& handle) {
  for (auto& raw : payload.at("retained_locals").as_array()) {
    auto& item = raw.as_object();
    if (std::string(item.at("handle").as_string()) == handle)
      return item;
  }
  throw std::runtime_error("checkpoint omitted retained local " + handle);
}

json::object checkpoint_payload(const std::string& path) {
  return json::parse(diffexp2::checkpoint::read(path).payload_json)
      .as_object();
}

json::object write_corrupt_row_checkpoint(
    const std::string& source_path, const std::string& corrupt_path,
    const std::string& derived_handle) {
  const auto container = diffexp2::checkpoint::read(source_path);
  auto header = json::parse(container.header_json).as_object();
  auto payload = json::parse(container.payload_json).as_object();
  auto& item = retained_local_record(payload, derived_handle);
  item.at("retained_owner_lineage").as_object()["row_exact_identity"] =
      "forged-row-identity";
  diffexp2::checkpoint::write_atomic(
      corrupt_path, json::serialize(header), json::serialize(payload));
  return payload;
}

bool rational_protocol() {
  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"},
      {"domain", "rational"}, {"output_digits", 50},
      {"local_capacity", 8}});
  const auto session = std::string(created.at("session").as_string());
  const auto anchor = prepare_chart(session, "rational", "anchor", "0");
  const auto lower = prepare_chart(session, "rational", "lower", "-2/3");
  const auto upper = prepare_chart(session, "rational", "upper", "2/3");
  const auto source = solve_local(session, anchor, "row-source-rational");
  const auto source_handle = std::string(source.at("local").as_string());

  auto malformed_row = rational_row();
  malformed_row.at("entries").as_array().front().as_object()
      .at("multiplier").as_object().at("kernels").as_array()
      .front() = json::array{};
  const auto malformed = apply_row(
      session, source_handle, std::move(malformed_row),
      "row-source-rational", "malformed-row");
  const auto projected = apply_row(
      session, source_handle, rational_row(), "row-source-rational",
      "row-result-rational");
  const auto projected_handle =
      std::string(projected.at("local").as_string());
  const auto zero = apply_row(
      session, source_handle, zero_row(), "row-source-rational",
      "zero-row-result");
  const auto zero_handle = std::string(zero.at("local").as_string());
  const auto chained = apply_row(
      session, projected_handle, scalar_identity_row(),
      "row-result-rational", "row-chain-rational");
  const auto chained_handle =
      std::string(chained.at("local").as_string());

  const auto plan = request(json::object{
      {"schema", 2}, {"op", "tile.plan"}, {"session", session},
      {"checkpoint_identity", "row-plan"}, {"division_order", 3},
      {"lower", arm("-2/3", anchor, lower)},
      {"upper", arm("2/3", anchor, upper)}});
  const auto plan_handle = std::string(plan.at("tile_plan").as_string());

  // The projected local strongly owns its source and remains fully usable
  // after the public vector-local token is released.
  (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                             {"session", session},
                             {"local", source_handle}});
  (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                             {"session", session},
                             {"local", projected_handle}});
  const auto positive = request(json::object{
      {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
      {"local", chained_handle},
      {"point", json::object{{"exact", "1/4"}}},
      {"options", json::object{{"tail_estimate", false}}}});
  const auto negative = request(json::object{
      {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
      {"local", chained_handle},
      {"point", json::object{{"exact", "-1/4"}}},
      {"options", json::object{{"tail_estimate", false}}}});
  const auto zero_value = request(json::object{
      {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
      {"local", zero_handle},
      {"point", json::object{{"exact", "1/4"}}},
      {"options", json::object{{"tail_estimate", false}}}});

  const auto line = request(json::object{
      {"schema", 2}, {"op", "integration.line"},
      {"session", session}, {"tile_plan", plan_handle},
      {"local", chained_handle}, {"arm", "upper"}, {"tile", 0},
      {"epsilon", json::object{{"min", 0}, {"max", 2}}},
      {"source_checkpoint_identity", "row-chain-rational"},
      {"tile_plan_checkpoint_identity", "row-plan"},
      {"checkpoint_identity", "row-line"}});
  const auto line_handle = std::string(line.at("line").as_string());
  const auto exported = request(json::object{
      {"schema", 2}, {"op", "integration.export"},
      {"session", session}, {"line", line_handle},
      {"checkpoint_identity", "row-line"}, {"output_digits", 50}});
  const auto stats = request(json::object{
      {"schema", 2}, {"op", "session.stats"}, {"session", session}});
  (void)request(json::object{{"schema", 2}, {"op", "integration.release"},
                             {"session", session}, {"line", line_handle}});

  const auto base = std::filesystem::temp_directory_path() /
      ("diffexp2-rational-row-rational-" + std::to_string(::getpid()));
  const auto path = base.string() + ".de2cp";
  const auto corrupt_path = base.string() + "-corrupt.de2cp";
  const auto resaved_path = base.string() + "-resaved.de2cp";
  const auto checkpoint_identity = "rational-row-roundtrip-v1";
  const auto saved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
      {"path", path}, {"checkpoint_identity", checkpoint_identity}});
  const auto saved_payload = checkpoint_payload(path);
  const auto& saved_visibility = saved_payload.at("session").as_object()
      .at("registry_visibility").as_object();
  (void)write_corrupt_row_checkpoint(
      path, corrupt_path, chained_handle);
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  const auto corruption = request(json::object{
      {"schema", 2}, {"op", "checkpoint.restore"},
      {"path", corrupt_path}, {"expected_identity", checkpoint_identity}});
  const auto restored = request(json::object{
      {"schema", 2}, {"op", "checkpoint.restore"}, {"path", path},
      {"expected_identity", checkpoint_identity}});
  const auto restored_session =
      std::string(restored.at("session").as_string());
  const auto resaved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"},
      {"session", restored_session}, {"path", resaved_path},
      {"checkpoint_identity", "rational-row-resaved-v1"}});
  const auto resaved_payload = checkpoint_payload(resaved_path);
  const auto hidden_source = request(json::object{
      {"schema", 2}, {"op", "local.stats"},
      {"session", restored_session}, {"local", source_handle}});
  const auto hidden_intermediate = request(json::object{
      {"schema", 2}, {"op", "local.stats"},
      {"session", restored_session}, {"local", projected_handle}});
  const auto restored_stats = request(json::object{
      {"schema", 2}, {"op", "session.stats"},
      {"session", restored_session}});
  const auto restored_positive = request(json::object{
      {"schema", 2}, {"op", "local.evaluate"},
      {"session", restored_session}, {"local", chained_handle},
      {"point", json::object{{"exact", "1/4"}}},
      {"options", json::object{{"tail_estimate", false}}}});
  const auto restored_line = request(json::object{
      {"schema", 2}, {"op", "integration.line"},
      {"session", restored_session}, {"tile_plan", plan_handle},
      {"local", chained_handle}, {"arm", "upper"}, {"tile", 0},
      {"epsilon", json::object{{"min", 0}, {"max", 2}}},
      {"source_checkpoint_identity", "row-chain-rational"},
      {"tile_plan_checkpoint_identity", "row-plan"},
      {"checkpoint_identity", "row-line-restored"}});
  const auto next_local = solve_local(
      restored_session, anchor, "row-next-rational");
  const auto final_stats = request(json::object{
      {"schema", 2}, {"op", "session.stats"},
      {"session", restored_session}});

  const auto projected_json = json::serialize(projected);
  const auto expected_eps0 = 2.0 / std::sqrt(3.0);
  const auto expected_eps1 = 4.0 / (3.0 * std::sqrt(3.0));
  const bool ok =
      created.at("rational_row_application_capability") ==
          "retained-native-rational-row-local-application-v1" &&
      malformed.at("status") == "error" &&
      std::string(malformed.at("detail").as_string()).find(
          "exact source Taylor width") != std::string::npos &&
      projected.at("status") == "ok" &&
      projected.at("application_capability") ==
          "retained-native-rational-row-local-application-v1" &&
      projected.at("native_retained") == true &&
      projected.at("json_coefficients") == 0 &&
      projected.at("dimension") == 1 &&
      projected.at("epsilon_min") == 0 &&
      projected.at("epsilon_max") == 2 &&
      projected.at("checkpoint_identity") == "row-result-rational" &&
      projected.at("strong_derivation_ownership") == true &&
      chained.at("status") == "ok" &&
      chained.at("checkpoint_identity") == "row-chain-rational" &&
      chained.at("strong_derivation_ownership") == true &&
      chained.at("tail_majorant").as_object().at("status") ==
          "unsupported" &&
      projected_json.find("\"kernels\"") == std::string::npos &&
      projected_json.find("\"coefficients\"") == std::string::npos &&
      zero.at("dimension") == 1 &&
      zero.at("epsilon_min") == 0 && zero.at("epsilon_max") == 2 &&
      positive.at("status") == "ok" && positive.at("imaginary_sign").is_null() &&
      std::abs(coefficient_midpoint(positive, 0) - 2.0) < 1e-30 &&
      std::abs(coefficient_midpoint(positive, 1) - 3.0) < 1e-30 &&
      negative.at("status") == "ok" && negative.at("imaginary_sign") == -1 &&
      std::abs(coefficient_midpoint(zero_value, 0)) < 1e-30 &&
      line.at("status") == "ok" && line.at("dimension") == 1 &&
      line.at("json_coefficients") == 0 &&
      exported.at("json_coefficients") == 3 &&
      std::abs(coefficient_midpoint(exported, 0) - expected_eps0) < 1e-14 &&
      std::abs(coefficient_midpoint(exported, 1) - expected_eps1) < 1e-14 &&
      stats.at("locals") == 2 && stats.at("pending_local_solves") == 0 &&
      saved.at("status") == "ok" &&
      saved_payload.at("retained_locals").as_array().size() == 4 &&
      saved_visibility.at("locals").as_array().size() == 2 &&
      std::find(saved_visibility.at("locals").as_array().begin(),
                saved_visibility.at("locals").as_array().end(),
                json::value(source_handle)) ==
          saved_visibility.at("locals").as_array().end() &&
      std::find(saved_visibility.at("locals").as_array().begin(),
                saved_visibility.at("locals").as_array().end(),
                json::value(projected_handle)) ==
          saved_visibility.at("locals").as_array().end() &&
      corruption.at("status") == "error" &&
      std::string(corruption.at("detail").as_string()).find(
          "rational-row owner lineage") != std::string::npos &&
      restored.at("status") == "ok" &&
      restored.at("locals").as_array().size() == 2 &&
      resaved.at("status") == "ok" &&
      retained_local_record(saved_payload, source_handle) ==
          retained_local_record(resaved_payload, source_handle) &&
      retained_local_record(saved_payload, projected_handle) ==
          retained_local_record(resaved_payload, projected_handle) &&
      retained_local_record(saved_payload, chained_handle) ==
          retained_local_record(resaved_payload, chained_handle) &&
      saved_visibility == resaved_payload.at("session").as_object()
                              .at("registry_visibility").as_object() &&
      hidden_source.at("status") == "error" &&
      hidden_intermediate.at("status") == "error" &&
      restored_stats.at("locals") == 2 &&
      restored_stats.at("line_results") == 0 &&
      restored_positive.at("status") == "ok" &&
      std::abs(coefficient_midpoint(restored_positive, 0) - 2.0) < 1e-30 &&
      restored_line.at("status") == "ok" &&
      restored_line.at("line") == "line:2" &&
      next_local.at("local") == "l:6" &&
      final_stats.at("locals") == 3 &&
      final_stats.at("line_results") == 1;

  if (!ok) {
    std::cerr << "malformed: " << json::serialize(malformed) << '\n'
              << "projected: " << projected_json << '\n'
              << "positive: " << json::serialize(positive) << '\n'
              << "negative: " << json::serialize(negative) << '\n'
              << "zero: " << json::serialize(zero_value) << '\n'
              << "line: " << json::serialize(line) << '\n'
              << "exported: " << json::serialize(exported) << '\n'
              << "stats: " << json::serialize(stats) << '\n'
              << "saved: " << json::serialize(saved) << '\n'
              << "corruption: " << json::serialize(corruption) << '\n'
              << "restored: " << json::serialize(restored) << '\n'
              << "resaved: " << json::serialize(resaved) << '\n'
              << "hidden source: " << json::serialize(hidden_source) << '\n'
              << "hidden intermediate: "
              << json::serialize(hidden_intermediate) << '\n'
              << "restored positive: "
              << json::serialize(restored_positive) << '\n'
              << "restored line: " << json::serialize(restored_line) << '\n'
              << "next local: " << json::serialize(next_local) << '\n'
              << "final stats: " << json::serialize(final_stats) << '\n';
  }

  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", restored_session}});
  std::error_code ignored;
  std::filesystem::remove(path, ignored);
  std::filesystem::remove(corrupt_path, ignored);
  std::filesystem::remove(resaved_path, ignored);
  return ok;
}

bool acb_protocol() {
  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", "acb"},
      {"precision_bits", 256}, {"output_digits", 50},
      {"local_capacity", 4}});
  const auto session = std::string(created.at("session").as_string());
  const auto chart = prepare_chart(session, "acb", "anchor", "0");
  const auto lower = prepare_chart(session, "acb", "lower", "-2/3");
  const auto upper = prepare_chart(session, "acb", "upper", "2/3");
  const auto source = solve_local(session, chart, "row-source-acb");
  const auto source_handle = std::string(source.at("local").as_string());
  json::array acb_leading;
  acb_leading.push_back(json::array{"[1 +/- 1e-40]", "0"});
  auto uncertain = multiplier(
      0, 0, std::move(acb_leading),
      "uncertain-unit");
  const auto row = json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", 2}, {"exact_identity", "[uncertain-one,0]"},
      {"entries", json::array{json::object{
           {"column", 0}, {"multiplier", std::move(uncertain)}}}}};
  const auto projected = apply_row(
      session, source_handle, row, "row-source-acb", "row-result-acb");
  if (projected.at("status") != "ok") {
    std::cerr << "Acb apply failed: " << json::serialize(projected) << '\n';
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    return false;
  }
  const auto projected_handle =
      std::string(projected.at("local").as_string());
  const auto chained = apply_row(
      session, projected_handle, scalar_identity_row(),
      "row-result-acb", "row-chain-acb");
  if (chained.at("status") != "ok") {
    std::cerr << "Acb chained apply failed: " << json::serialize(chained)
              << '\n';
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    return false;
  }
  const auto chained_handle =
      std::string(chained.at("local").as_string());
  const auto plan = request(json::object{
      {"schema", 2}, {"op", "tile.plan"}, {"session", session},
      {"checkpoint_identity", "row-plan-acb"}, {"division_order", 3},
      {"lower", arm("-2/3", chart, lower)},
      {"upper", arm("2/3", chart, upper)}});
  const auto plan_handle = std::string(plan.at("tile_plan").as_string());
  (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                             {"session", session},
                             {"local", source_handle}});
  (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                             {"session", session},
                             {"local", projected_handle}});
  const auto evaluated = request(json::object{
      {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
      {"local", chained_handle},
      {"point", json::object{{"exact", "1/4"}}},
      {"options", json::object{{"tail_estimate", false}}}});
  const auto line = request(json::object{
      {"schema", 2}, {"op", "integration.line"},
      {"session", session}, {"tile_plan", plan_handle},
      {"local", chained_handle}, {"arm", "upper"}, {"tile", 0},
      {"epsilon", json::object{{"min", 0}, {"max", 2}}},
      {"source_checkpoint_identity", "row-chain-acb"},
      {"tile_plan_checkpoint_identity", "row-plan-acb"},
      {"checkpoint_identity", "row-line-acb"}});
  const auto line_handle = std::string(line.at("line").as_string());
  (void)request(json::object{{"schema", 2}, {"op", "integration.release"},
                             {"session", session}, {"line", line_handle}});

  const auto base = std::filesystem::temp_directory_path() /
      ("diffexp2-rational-row-acb-" + std::to_string(::getpid()));
  const auto path = base.string() + ".de2cp";
  const auto corrupt_path = base.string() + "-corrupt.de2cp";
  const auto resaved_path = base.string() + "-resaved.de2cp";
  const auto checkpoint_identity = "acb-rational-row-roundtrip-v1";
  const auto saved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
      {"path", path}, {"checkpoint_identity", checkpoint_identity}});
  const auto saved_payload = checkpoint_payload(path);
  const auto& saved_visibility = saved_payload.at("session").as_object()
      .at("registry_visibility").as_object();
  (void)write_corrupt_row_checkpoint(
      path, corrupt_path, chained_handle);
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  const auto corruption = request(json::object{
      {"schema", 2}, {"op", "checkpoint.restore"},
      {"path", corrupt_path}, {"expected_identity", checkpoint_identity}});
  const auto restored = request(json::object{
      {"schema", 2}, {"op", "checkpoint.restore"}, {"path", path},
      {"expected_identity", checkpoint_identity}});
  const auto restored_session =
      std::string(restored.at("session").as_string());
  const auto resaved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"},
      {"session", restored_session}, {"path", resaved_path},
      {"checkpoint_identity", "acb-rational-row-resaved-v1"}});
  const auto resaved_payload = checkpoint_payload(resaved_path);
  const auto hidden_source = request(json::object{
      {"schema", 2}, {"op", "local.stats"},
      {"session", restored_session}, {"local", source_handle}});
  const auto hidden_intermediate = request(json::object{
      {"schema", 2}, {"op", "local.stats"},
      {"session", restored_session}, {"local", projected_handle}});
  const auto restored_stats = request(json::object{
      {"schema", 2}, {"op", "session.stats"},
      {"session", restored_session}});
  const auto restored_evaluated = request(json::object{
      {"schema", 2}, {"op", "local.evaluate"},
      {"session", restored_session}, {"local", chained_handle},
      {"point", json::object{{"exact", "1/4"}}},
      {"options", json::object{{"tail_estimate", false}}}});
  const auto restored_line = request(json::object{
      {"schema", 2}, {"op", "integration.line"},
      {"session", restored_session}, {"tile_plan", plan_handle},
      {"local", chained_handle}, {"arm", "upper"}, {"tile", 0},
      {"epsilon", json::object{{"min", 0}, {"max", 2}}},
      {"source_checkpoint_identity", "row-chain-acb"},
      {"tile_plan_checkpoint_identity", "row-plan-acb"},
      {"checkpoint_identity", "row-line-acb-restored"}});
  const auto next_local = solve_local(
      restored_session, chart, "row-next-acb");
  const auto final_stats = request(json::object{
      {"schema", 2}, {"op", "session.stats"},
      {"session", restored_session}});
  const auto& first = evaluated.at("value").as_object()
      .at("coefficients").as_array().front().as_array();
  const auto& restored_first = restored_evaluated.at("value").as_object()
      .at("coefficients").as_array().front().as_array();
  const bool ok =
      created.at("rational_row_application_capability") ==
          "retained-native-rational-row-local-application-v1" &&
      projected.at("status") == "ok" && projected.at("dimension") == 1 &&
      projected.at("json_coefficients") == 0 &&
      projected.at("metadata").as_object()
          .at("prescriptions").as_array().size() == 1 &&
      projected.at("tail_majorant").as_object().at("status") ==
          "unsupported" &&
      chained.at("status") == "ok" &&
      chained.at("strong_derivation_ownership") == true &&
      chained.at("tail_majorant").as_object().at("status") ==
          "unsupported" &&
      evaluated.at("status") == "ok" && evaluated.at("arithmetic_enclosed") == true &&
      first.size() == 4 && first.at(2) != "zero" &&
      line.at("status") == "ok" &&
      saved.at("status") == "ok" &&
      saved_payload.at("retained_locals").as_array().size() == 3 &&
      saved_visibility.at("locals").as_array().size() == 1 &&
      saved_visibility.at("locals").as_array().front() ==
          json::value(chained_handle) &&
      corruption.at("status") == "error" &&
      std::string(corruption.at("detail").as_string()).find(
          "rational-row owner lineage") != std::string::npos &&
      restored.at("status") == "ok" &&
      restored.at("locals").as_array().size() == 1 &&
      resaved.at("status") == "ok" &&
      retained_local_record(saved_payload, source_handle) ==
          retained_local_record(resaved_payload, source_handle) &&
      retained_local_record(saved_payload, projected_handle) ==
          retained_local_record(resaved_payload, projected_handle) &&
      retained_local_record(saved_payload, chained_handle) ==
          retained_local_record(resaved_payload, chained_handle) &&
      saved_visibility == resaved_payload.at("session").as_object()
                              .at("registry_visibility").as_object() &&
      hidden_source.at("status") == "error" &&
      hidden_intermediate.at("status") == "error" &&
      restored_stats.at("locals") == 1 &&
      restored_stats.at("line_results") == 0 &&
      restored_evaluated.at("status") == "ok" &&
      restored_evaluated.at("arithmetic_enclosed") == true &&
      restored_first == first &&
      restored_line.at("status") == "ok" &&
      restored_line.at("line") == "line:2" &&
      next_local.at("local") == "l:4" &&
      final_stats.at("locals") == 2 &&
      final_stats.at("line_results") == 1;
  if (!ok) {
    std::cerr << "Acb projected: " << json::serialize(projected) << '\n'
              << "Acb evaluated: " << json::serialize(evaluated) << '\n'
              << "Acb line: " << json::serialize(line) << '\n'
              << "Acb saved: " << json::serialize(saved) << '\n'
              << "Acb corruption: " << json::serialize(corruption) << '\n'
              << "Acb restored: " << json::serialize(restored) << '\n'
              << "Acb resaved: " << json::serialize(resaved) << '\n'
              << "Acb hidden source: " << json::serialize(hidden_source)
              << '\n' << "Acb hidden intermediate: "
              << json::serialize(hidden_intermediate) << '\n'
              << "Acb restored evaluation: "
              << json::serialize(restored_evaluated) << '\n'
              << "Acb restored line: " << json::serialize(restored_line)
              << '\n' << "Acb next local: " << json::serialize(next_local)
              << '\n' << "Acb final stats: " << json::serialize(final_stats)
              << '\n';
  }
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", restored_session}});
  std::error_code ignored;
  std::filesystem::remove(path, ignored);
  std::filesystem::remove(corrupt_path, ignored);
  std::filesystem::remove(resaved_path, ignored);
  return ok;
}

}  // namespace

int main() {
  const bool ok = rational_protocol() && acb_protocol();
  std::cout << (ok ? "PASS" : "FAIL")
            << ": retained Rational/Acb rational-row projection and tile integration\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
