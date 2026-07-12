#include "diffexp2/checkpoint.hpp"
#include "diffexp2/json_codec.hpp"

#include <boost/crc.hpp>
#include <boost/json.hpp>

#include <array>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <unistd.h>

namespace json = boost::json;

namespace {

struct RawContainer {
  std::string header;
  std::string payload;
};

void append_u32(std::vector<unsigned char>& output, std::uint32_t value) {
  for (int shift = 24; shift >= 0; shift -= 8)
    output.push_back(static_cast<unsigned char>((value >> shift) & 0xffu));
}

void append_u64(std::vector<unsigned char>& output, std::uint64_t value) {
  for (int shift = 56; shift >= 0; shift -= 8)
    output.push_back(static_cast<unsigned char>((value >> shift) & 0xffu));
}

std::uint32_t read_u32(const unsigned char* bytes) {
  std::uint32_t value = 0;
  for (int index = 0; index < 4; ++index)
    value = (value << 8) | bytes[index];
  return value;
}

std::uint32_t checksum(std::string_view bytes) {
  boost::crc_32_type crc;
  crc.process_bytes(bytes.data(), bytes.size());
  return crc.checksum();
}

void write_container_fixture(const std::filesystem::path& path,
                             std::uint32_t schema,
                             std::string_view header,
                             std::string_view payload) {
  std::vector<unsigned char> bytes{
      'D', 'E', '2', 'C', 'P', '0', '0', '1'};
  append_u32(bytes, schema);
  append_u32(bytes, static_cast<std::uint32_t>(header.size()));
  append_u64(bytes, payload.size());
  append_u32(bytes, checksum(header));
  append_u32(bytes, checksum(payload));
  bytes.insert(bytes.end(), header.begin(), header.end());
  bytes.insert(bytes.end(), payload.begin(), payload.end());
  std::ofstream file(path, std::ios::binary | std::ios::trunc);
  file.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
  if (!file) throw std::runtime_error("could not write checkpoint fixture");
}

RawContainer read_raw_container_fixture(const std::filesystem::path& path) {
  std::ifstream file(path, std::ios::binary);
  std::array<unsigned char, 32> fixed{};
  file.read(reinterpret_cast<char*>(fixed.data()), fixed.size());
  if (!file) throw std::runtime_error("could not read checkpoint fixture");
  const auto header_size = read_u32(fixed.data() + 12);
  std::uint64_t payload_size = 0;
  for (int index = 0; index < 8; ++index)
    payload_size = (payload_size << 8) | fixed[16 + index];
  RawContainer output;
  output.header.resize(header_size);
  output.payload.resize(static_cast<std::size_t>(payload_size));
  file.read(output.header.data(), output.header.size());
  file.read(output.payload.data(), output.payload.size());
  if (!file) throw std::runtime_error("checkpoint fixture is truncated");
  return output;
}

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
      .as_object();
}

std::string prepare_rational_chart(const std::string& session) {
  auto response = request(json::parse(std::string(R"json({
    "schema":2,"op":"chart.prepare","session":")json") + session +
      R"json(","key":"checkpoint-rational-chart",
    "identity":"checkpoint-rational-chart-v1",
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
  })json").as_object());
  if (response.at("status") != "ok")
    throw std::runtime_error("rational prepare: " +
                             json::serialize(response));
  return std::string(response.at("chart").as_string());
}

std::string prepare_acb_chart(const std::string& session) {
  auto response = request(json::parse(std::string(R"json({
    "schema":2,"op":"chart.prepare","session":")json") + session +
      R"json(","key":"checkpoint-acb-chart",
    "identity":"checkpoint-acb-chart-v1",
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
    throw std::runtime_error("Acb prepare: " + json::serialize(response));
  return std::string(response.at("chart").as_string());
}

json::object metadata(const std::string& checkpoint_identity) {
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
      {"checkpoint_identity", checkpoint_identity}};
}

json::object rational_run(const std::string& leading) {
  json::array schedule_row;
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
      {"initial", json::array{leading, "0"}},
      {"initial_validity", json::array{1}},
      {"source", nullptr}, {"return_u", false}};
}

json::object acb_run(std::vector<std::string> initial) {
  json::array encoded_initial;
  for (auto& value : initial)
    encoded_initial.emplace_back(std::move(value));
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

std::string solve_local(const std::string& session, const std::string& chart,
                        json::object run,
                        const std::string& checkpoint_identity) {
  auto response = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart}, {"run", std::move(run)},
      {"metadata", metadata(checkpoint_identity)}});
  if (response.at("status") != "ok")
    throw std::runtime_error("local solve: " + json::serialize(response));
  return std::string(response.at("local").as_string());
}

std::vector<std::string> acb_column(std::string first_component,
                                    std::string second_epsilon = "0") {
  return {std::move(first_component), "0", "0", "0", "0", "0",
          "0", std::move(second_epsilon), "0", "0", "0", "0"};
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
  return json::object{
      {"schema", "diffexp2-exact-evaluated-epsilon-lattice-v1"},
      {"identity", "checkpoint-nontrivial-lattice-v1"},
      {"evaluated_basis", json::array{
          json::array{exact_frame(one), exact_frame(one)},
          json::array{exact_frame(zero), exact_frame(epsilon)}}}};
}

json::object match_request(const std::string& session,
                           const std::string& chart,
                           const std::string& first,
                           const std::string& second,
                           const std::string& incoming) {
  return json::object{
      {"schema", 2}, {"op", "local.match_acb"}, {"session", session},
      {"basis", json::array{first, second}}, {"incoming", incoming},
      {"basis_chart", chart}, {"incoming_chart", chart},
      {"basis_point", json::object{{"exact", "1/2"}}},
      {"incoming_point", json::object{{"exact", "1/2"}}},
      {"epsilon", json::object{{"min", 0}, {"max", 5},
                                {"required_complete_max", 3}}},
      {"basis_checkpoint_identities", json::array{"basis-0", "basis-1"}},
      {"incoming_checkpoint_identity", "incoming"},
      {"checkpoint_identity", "checkpoint-acb-match-v1"},
      {"exact_lattice", exact_lattice()},
      {"refinement", json::object{{"relative_tolerance", "1e-50"},
                                    {"max_steps", 2}}}};
}

json::object endpoint_request(const std::string& session,
                              const std::string& local,
                              const std::string& source_checkpoint,
                              const std::string& checkpoint_identity) {
  return json::object{
      {"schema", 2}, {"op", "local.endpoint_limit"},
      {"session", session}, {"local", local},
      {"source_checkpoint_identity", source_checkpoint},
      {"checkpoint_identity", checkpoint_identity},
      {"approach_direction", 1},
      {"cancellation", json::object{{"mode", "exact-or-acb-singleton"}}}};
}

json::object payload_at(const std::filesystem::path& path) {
  return json::parse(diffexp2::checkpoint::read(path.string()).payload_json)
      .as_object();
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

bool canonical_json_dag_container_roundtrip(
    const std::filesystem::path& path,
    const std::filesystem::path& legacy_path,
    const std::filesystem::path& unused_path) {
  const std::string blob(4096, 'x');
  const json::object leaf{{"schema", "nested-leaf-v1"}, {"blob", blob}};
  const auto leaf_identity = json::serialize(leaf);
  json::array leaf_references;
  for (int index = 0; index < 8; ++index)
    leaf_references.emplace_back(leaf_identity);
  const json::object provenance{
      {"schema", "nested-provenance-v1"},
      {"leaf", leaf},
      {"leaf_identity", leaf_identity},
      {"leaf_references", std::move(leaf_references)}};
  const auto provenance_identity = json::serialize(provenance);
  json::array provenance_references;
  for (int index = 0; index < 8; ++index)
    provenance_references.emplace_back(provenance_identity);
  const json::object closure{
      {"schema", "nested-closure-v1"},
      {"provenance", provenance},
      {"provenance_identity", provenance_identity},
      {"provenance_references", std::move(provenance_references)}};
  const auto closure_identity = json::serialize(closure);
  std::string reserved_literal;
  reserved_literal.push_back('\0');
  reserved_literal += "R7";
  const json::object header{
      {"schema", 8},
      {"closure_identity", closure_identity},
      {"closure_identity_copy", closure_identity},
      {"reserved_literal", reserved_literal}};
  const json::object payload{
      {"schema", 8},
      {"closure", closure},
      {"closure_identity", closure_identity},
      {"closure_identity_copy", closure_identity},
      {"reserved_literal", reserved_literal}};

  diffexp2::checkpoint::write_atomic(path.string(), header, payload);
  const auto decoded = diffexp2::checkpoint::read_json(path.string());
  const auto text = diffexp2::checkpoint::read(path.string());
  write_container_fixture(legacy_path, 1, text.header_json,
                          text.payload_json);
  const auto legacy = diffexp2::checkpoint::read_json(legacy_path.string());
  auto raw = read_raw_container_fixture(path);
  auto raw_payload = json::parse(raw.payload).as_object();
  raw_payload.at("canonical_json_strings").as_array().push_back(
      json::object{{"unused", true}});
  write_container_fixture(unused_path, 2, raw.header,
                          json::serialize(raw_payload));
  bool unused_rejected = false;
  try {
    (void)diffexp2::checkpoint::read_json(unused_path.string());
  } catch (const std::invalid_argument& error) {
    unused_rejected =
        std::string(error.what()).find("unused entries") != std::string::npos;
  }
  std::array<unsigned char, 12> fixed{};
  std::ifstream file(path, std::ios::binary);
  file.read(reinterpret_cast<char*>(fixed.data()), fixed.size());
  const auto semantic_bytes = json::serialize(header).size() +
      json::serialize(payload).size() + 32;
  const bool schema_v2 = file && fixed[8] == 0 && fixed[9] == 0 &&
      fixed[10] == 0 && fixed[11] == 2;
  const bool compact = std::filesystem::file_size(path) < semantic_bytes / 8;
  const bool exact = decoded.header == header && decoded.payload == payload &&
      legacy.header == header && legacy.payload == payload &&
      json::parse(text.header_json) == header &&
      json::parse(text.payload_json) == payload;
  if (!schema_v2 || !compact || !exact || !unused_rejected)
    std::cerr << "canonical JSON DAG codec: schema_v2=" << schema_v2
              << " compact=" << compact << " exact=" << exact
              << " unused_rejected=" << unused_rejected
              << " file_bytes=" << std::filesystem::file_size(path)
              << " semantic_bytes=" << semantic_bytes << '\n';
  return schema_v2 && compact && exact && unused_rejected;
}

bool rational_local_roundtrip(const std::filesystem::path& path) {
  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", "rational"},
      {"precision_bits", 256}, {"output_digits", 30},
      {"chart_capacity", 2}, {"local_capacity", 4}});
  const auto session = std::string(created.at("session").as_string());
  const auto chart = prepare_rational_chart(session);
  const auto local = solve_local(session, chart, rational_run("7/3"),
                                 "checkpoint-rational-local-v1");
  const auto saved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
      {"path", path.string()},
      {"checkpoint_identity", "checkpoint-rational-session-v2"}});
  const auto original_payload = payload_at(path);
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});
  const auto restored = request(json::object{
      {"schema", 2}, {"op", "checkpoint.restore"},
      {"path", path.string()},
      {"expected_identity", "checkpoint-rational-session-v2"}});
  if (restored.at("status") != "ok") {
    std::cerr << "rational restore: " << json::serialize(restored) << '\n';
    return false;
  }
  const auto restored_session =
      std::string(restored.at("session").as_string());
  const auto retained = request(json::object{
      {"schema", 2}, {"op", "local.stats"},
      {"session", restored_session}, {"local", local}});
  const auto next = solve_local(restored_session, chart, rational_run("5/2"),
                                "checkpoint-rational-postrestore-v1");
  const auto& local_items = original_payload.at("retained_locals").as_array();
  const auto& coefficient = local_items.front().as_object()
      .at("solution").as_object().at("sectors").as_array().front().as_object()
      .at("coefficients").as_array().front();
  const bool ok = saved.at("status") == "ok" && saved.at("locals") == 1 &&
      restored.at("locals").as_array().size() == 1 &&
      retained.at("status") == "ok" &&
      std::string(retained.at("local").as_string()) == local &&
      retained.at("checkpoint_identity") ==
          "checkpoint-rational-local-v1" && coefficient == "7/3" &&
      next == "l:2";
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", restored_session}});
  if (!ok)
    std::cerr << "rational saved: " << json::serialize(saved) << '\n'
              << "rational retained: " << json::serialize(retained) << '\n'
              << "rational restored: " << json::serialize(restored) << '\n';
  return ok;
}

bool acb_match_roundtrip(const std::filesystem::path& path,
                         const std::filesystem::path& resaved_path,
                         const std::filesystem::path& corrupted_path) {
  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", "acb"},
      {"precision_bits", 256}, {"output_digits", 30},
      {"chart_capacity", 2}, {"local_capacity", 8},
      {"match_capacity", 4}});
  const auto session = std::string(created.at("session").as_string());
  const auto chart = prepare_acb_chart(session);
  const auto first = solve_local(session, chart, acb_run(acb_column("1")),
                                 "basis-0");
  const auto second = solve_local(
      session, chart, acb_run(acb_column("1", "1")), "basis-1");
  const auto incoming = solve_local(
      session, chart, acb_run(acb_column("3", "2")), "incoming");
  const auto matched = request(match_request(session, chart, first, second,
                                             incoming));
  if (matched.at("status") != "ok")
    throw std::runtime_error("Acb match: " + json::serialize(matched));
  const auto match = std::string(matched.at("match").as_string());
  const auto endpoint_result = request(endpoint_request(
      session, first, "basis-0", "checkpoint-endpoint-v1"));
  if (endpoint_result.at("status") != "ok")
    throw std::runtime_error("endpoint: " +
                             json::serialize(endpoint_result));
  const auto endpoint =
      std::string(endpoint_result.at("endpoint").as_string());
  (void)request(json::object{
      {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
      {"local", first}, {"point", json::object{{"exact", "1/2"}}},
      {"options", json::object{{"tail_estimate", false}}}});
  const auto saved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
      {"path", path.string()},
      {"checkpoint_identity", "checkpoint-acb-session-v2"}});
  const auto original_payload = payload_at(path);
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});

  const auto restored = request(json::object{
      {"schema", 2}, {"op", "checkpoint.restore"},
      {"path", path.string()},
      {"expected_identity", "checkpoint-acb-session-v2"}});
  if (restored.at("status") != "ok") {
    std::cerr << "Acb restore: " << json::serialize(restored) << '\n';
    return false;
  }
  const auto restored_session =
      std::string(restored.at("session").as_string());
  const auto resaved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"},
      {"session", restored_session}, {"path", resaved_path.string()},
      {"checkpoint_identity", "checkpoint-acb-resaved-v2"}});
  const auto restored_payload = payload_at(resaved_path);
  const auto retained_match = request(json::object{
      {"schema", 2}, {"op", "match.stats"},
      {"session", restored_session}, {"match", match}});
  const auto retained_local = request(json::object{
      {"schema", 2}, {"op", "local.stats"},
      {"session", restored_session}, {"local", first}});
  const auto retained_endpoint = request(json::object{
      {"schema", 2}, {"op", "endpoint.stats"},
      {"session", restored_session}, {"endpoint", endpoint}});
  const auto next = solve_local(restored_session, chart,
      acb_run(acb_column("4")), "postrestore");
  const auto next_endpoint_result = request(endpoint_request(
      restored_session, first, "basis-0", "postrestore-endpoint"));

  const bool corrupted = corrupt_last_payload_byte(path, corrupted_path);
  const auto corruption_result = request(json::object{
      {"schema", 2}, {"op", "checkpoint.restore"},
      {"path", corrupted_path.string()},
      {"expected_identity", "checkpoint-acb-session-v2"}});

  const auto& original_locals = original_payload.at("retained_locals");
  const auto& original_matches = original_payload.at("retained_acb_matches");
  const auto& original_endpoints = original_payload.at("retained_endpoints");
  const auto serialized_live_state = json::serialize(json::object{
      {"locals", original_locals}, {"matches", original_matches},
      {"endpoints", original_endpoints}});
  const auto& first_ball = original_locals.as_array().front().as_object()
      .at("solution").as_object().at("sectors").as_array().front().as_object()
      .at("coefficients").as_array().front().as_object();
  const auto& match_state = original_matches.as_array().front().as_object();
  const auto& refined = match_state.at("refined").as_object();
  const auto& first_weight_ball = refined.at("weights").as_array().front()
      .as_object().at("coefficients").as_array().front().as_object();
  const auto& first_bound = refined.at("residual_history").as_array().front()
      .as_object().at("coefficients").as_array().front().as_object()
      .at("residual_upper_exact");
  const auto& endpoint_ball = original_endpoints.as_array().front().as_object()
      .at("result").as_object().at("values").as_array().front().as_object()
      .at("coefficients").as_array().front().as_object();
  const auto corruption_detail = corruption_result.at("status") == "error"
      ? std::string(corruption_result.at("detail").as_string())
      : std::string();
  const bool ok = saved.at("status") == "ok" && saved.at("locals") == 3 &&
      saved.at("acb_matches") == 1 && saved.at("endpoints") == 1 &&
      restored.at("status") == "ok" &&
      restored.at("locals").as_array().size() == 3 &&
      restored.at("acb_matches").as_array().size() == 1 &&
      restored.at("endpoints").as_array().size() == 1 &&
      resaved.at("status") == "ok" &&
      original_locals == restored_payload.at("retained_locals") &&
      original_matches == restored_payload.at("retained_acb_matches") &&
      original_endpoints == restored_payload.at("retained_endpoints") &&
      retained_match.at("status") == "ok" &&
      retained_match.at("provenance_identity") ==
          matched.at("provenance_identity") &&
      retained_match.at("residual").as_object().at("history").as_array()
          .size() >= 1 && retained_local.at("evaluations") == 1 &&
      retained_endpoint.at("status") == "ok" &&
      retained_endpoint.at("checkpoint_identity") ==
          "checkpoint-endpoint-v1" &&
      next == "l:4" &&
      next_endpoint_result.at("endpoint") == "e:2" &&
      first_ball.at("real").is_string() &&
      first_ball.at("imaginary").is_string() &&
      first_weight_ball.at("real").is_string() && first_bound.is_string() &&
      endpoint_ball.at("real").is_string() &&
      serialized_live_state.find("\"radius_ball\"") == std::string::npos &&
      serialized_live_state.find("radius_exact_ball") != std::string::npos &&
      corrupted && corruption_result.at("status") == "error" &&
      corruption_detail.find("checksum mismatch") != std::string::npos;
  if (!ok) {
    std::cerr << "Acb saved: " << json::serialize(saved) << '\n'
              << "Acb restored: " << json::serialize(restored) << '\n'
              << "Acb resaved: " << json::serialize(resaved) << '\n'
              << "Acb match: " << json::serialize(retained_match) << '\n'
              << "Acb local: " << json::serialize(retained_local) << '\n'
              << "endpoint: " << json::serialize(retained_endpoint) << '\n'
              << "corruption: " << json::serialize(corruption_result) << '\n';
  }
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", restored_session}});
  return ok;
}

}  // namespace

int main() {
  const auto stem = std::filesystem::temp_directory_path() /
      ("diffexp2_checkpoint_v2_" +
       std::to_string(static_cast<long long>(::getpid())));
  const auto rational_path = stem.string() + "_rational.de2cp";
  const auto acb_path = stem.string() + "_acb.de2cp";
  const auto resaved_path = stem.string() + "_acb_resaved.de2cp";
  const auto corrupted_path = stem.string() + "_acb_corrupt.de2cp";
  const auto codec_path = stem.string() + "_codec.de2cp";
  const auto legacy_path = stem.string() + "_legacy.de2cp";
  const auto unused_path = stem.string() + "_unused.de2cp";
  for (const auto& path : {rational_path, acb_path, resaved_path,
                           corrupted_path, codec_path, legacy_path,
                           unused_path})
    std::filesystem::remove(path);

  bool ok = false;
  try {
    ok = canonical_json_dag_container_roundtrip(
             codec_path, legacy_path, unused_path) &&
         rational_local_roundtrip(rational_path) &&
         acb_match_roundtrip(acb_path, resaved_path, corrupted_path);
  } catch (const std::exception& error) {
    std::cerr << "checkpoint v2 exception: " << error.what() << '\n';
  }
  for (const auto& path : {rational_path, acb_path, resaved_path,
                           corrupted_path, codec_path, legacy_path,
                           unused_path})
    std::filesystem::remove(path);
  std::cout << (ok ? "PASS" : "FAIL")
            << ": compact canonical-JSON DAG, exact live checkpoint roundtrip, and CRC corruption\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
