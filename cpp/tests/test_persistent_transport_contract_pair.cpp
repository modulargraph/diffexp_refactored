#include "diffexp2/checkpoint.hpp"
#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
      .as_object();
}

std::uint64_t counter(const json::object& object, const char* key) {
  const auto& value = object.at(key);
  if (value.is_uint64()) return value.as_uint64();
  if (value.is_int64() && value.as_int64() >= 0)
    return static_cast<std::uint64_t>(value.as_int64());
  throw std::runtime_error(std::string("invalid counter: ") + key);
}

json::array prescriptions() {
  return json::array{json::object{
      {"factor_exact", "t"}, {"sign", -1}, {"multiplicity", 1},
      {"leading_coefficient_sign", 1}}};
}

std::string prepare_anchor_chart(const std::string& session) {
  auto value = json::parse(R"json({
    "schema":2,"op":"chart.prepare","session":"placeholder",
    "key":"pair-anchor-chart","identity":"pair-anchor-operator",
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
  const auto prepared = request(std::move(value));
  if (prepared.at("status") != "ok")
    throw std::runtime_error("chart.prepare: " + json::serialize(prepared));
  return std::string(prepared.at("chart").as_string());
}

std::string solve_anchor(const std::string& session,
                         const std::string& chart) {
  auto value = json::parse(R"json({
    "schema":2,"op":"local.solve","session":"placeholder",
    "chart":"placeholder",
    "run":{"nmax":0,"p":0,"has_initial":true,
      "adaptive_probe":false,"a_target":"0","b_target":"0",
      "a_shift_min":0,"a_shifts":["0"],
      "schedule":[[{"case":"R","da":"0","db":"0"}]],
      "initial":["1","0","0"],"initial_validity":[2],
      "source":null,"return_u":false},
    "metadata":{"chart":{"center_exact":"0","scale_exact":"1",
        "radius":"2","infinite_radius":false},
      "tag":{"a":{"domain":"rational","canonical":"0"},
        "b":{"domain":"rational","canonical":"0"}},
      "prescriptions":[{"factor_exact":"t","sign":-1,
        "multiplicity":1,"leading_coefficient_sign":1}],
      "checkpoint_identity":"pair-common-anchor"}
  })json").as_object();
  value["session"] = session;
  value["chart"] = chart;
  const auto solved = request(std::move(value));
  if (solved.at("status") != "ok")
    throw std::runtime_error("local.solve: " + json::serialize(solved));
  return std::string(solved.at("local").as_string());
}

std::string solve_regulated_anchor(const std::string& session,
                                   const std::string& chart) {
  auto value = json::parse(R"json({
    "schema":2,"op":"local.solve","session":"placeholder",
    "chart":"placeholder",
    "run":{"nmax":0,"p":0,"has_initial":true,
      "adaptive_probe":false,"a_target":"-1","b_target":"1",
      "a_shift_min":0,"a_shifts":["-1"],
      "schedule":[[{"case":"R","da":"0","db":"0"}]],
      "initial":["1","0","0"],"initial_validity":[2],
      "source":null,"return_u":false},
    "metadata":{"chart":{"center_exact":"0","scale_exact":"1",
        "radius":"2","infinite_radius":false},
      "tag":{"a":{"domain":"rational","canonical":"-1"},
        "b":{"domain":"rational","canonical":"1"}},
      "prescriptions":[{"factor_exact":"t","sign":-1,
        "multiplicity":1,"leading_coefficient_sign":1}],
      "checkpoint_identity":"pair-common-anchor"}
  })json").as_object();
  value["session"] = session;
  value["chart"] = chart;
  const auto solved = request(std::move(value));
  if (solved.at("status") != "ok")
    throw std::runtime_error(
        "regulated local.solve: " + json::serialize(solved));
  return std::string(solved.at("local").as_string());
}

json::object topology() {
  return json::object{
      {"singular_points", json::array{}},
      {"boundary_points", json::array{}},
      {"complex_projections", json::array{}},
      {"branch_sheets", json::array{
           json::object{{"factor_exact", "t"}, {"sign", -1}}}}};
}

json::object one_chart_arm(const std::string& endpoint,
                           const std::string& chart) {
  return json::object{
      {"from_exact", "0"}, {"to_exact", endpoint},
      {"charts", json::array{chart}}, {"topology", topology()}};
}

json::object prepare_plan(const std::string& session,
                          const std::string& checkpoint,
                          const std::string& endpoint,
                          const std::string& chart) {
  return request(json::object{
      {"schema", 2}, {"op", "tile.plan_arm"}, {"session", session},
      {"checkpoint_identity", checkpoint}, {"division_order", 3},
      {"arm", one_chart_arm(endpoint, chart)}});
}

json::object prepare_pair_plan(const std::string& session,
                               const std::string& checkpoint,
                               const std::string& chart) {
  return request(json::object{
      {"schema", 2}, {"op", "tile.plan"}, {"session", session},
      {"checkpoint_identity", checkpoint}, {"division_order", 3},
      {"lower", one_chart_arm("-1/3", chart)},
      {"upper", one_chart_arm("1/2", chart)}});
}

json::object run_state(const std::string& session,
                       const json::object& plan,
                       const std::string& plan_checkpoint,
                       const std::string& anchor,
                       const std::string& arm,
                       const std::string& checkpoint_root) {
  return request(json::object{
      {"schema", 2}, {"op", "transport.run_arm"},
      {"session", session},
      {"tile_plan", plan.at("tile_plan")}, {"anchor", anchor},
      {"tile_plan_checkpoint_identity", plan_checkpoint},
      {"anchor_checkpoint_identity", "pair-common-anchor"},
      {"arm", arm}, {"receiving_basis", json::array{}},
      {"epsilon", json::object{
           {"min", 0}, {"max", 2}, {"required_complete_max", 0},
           {"match_required_complete_max", 2}}},
      {"refinement", json::object{
           {"relative_tolerance", "1e-30"}, {"max_steps", 0}}},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp2-deterministic-arm-checkpoints-v1"},
           {"root", checkpoint_root}}}});
}

json::object multiplier(const std::string& identity) {
  return json::object{
      {"epsilon_shift", 0}, {"center_pole_order", 0},
      {"kernels", json::array{
           json::array{"1"}, json::array{"0"}, json::array{"0"}}},
      {"exact_identity", identity}, {"proven_zero", false}};
}

json::object row(const std::string& identity, bool malformed = false) {
  return json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", malformed ? 2 : 1}, {"exact_identity", identity},
      {"entries", json::array{json::object{
           {"column", 0}, {"multiplier", multiplier(identity)}}}}};
}

json::object public_prefix_row(const std::string& identity) {
  json::array analytic;
  analytic.push_back(json::object{
      {"numerator", json::array{"1"}},
      {"denominator", json::array{"1"}}});
  analytic.push_back(json::object{
      {"numerator", json::array{"0"}},
      {"denominator", json::array{"1"}}});
  return json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", 1}, {"exact_identity", identity},
      {"entries", json::array{json::object{
           {"column", 0},
           {"multiplier", json::object{
                {"epsilon_shift", 0}, {"center_pole_order", 0},
                {"analytic_coefficients", std::move(analytic)},
                {"exact_identity", identity},
                {"proven_zero", false}}}}}}};
}

json::object regulated_row(const std::string& identity,
                           const std::string& coefficient) {
  return json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", 1}, {"exact_identity", identity},
      {"entries", json::array{json::object{
           {"column", 0},
           {"multiplier", json::object{
                {"epsilon_shift", 0}, {"center_pole_order", 0},
                {"kernels", json::array{
                     json::array{coefficient}, json::array{"0"},
                     json::array{"0"}}},
                {"exact_identity", coefficient},
                {"proven_zero", false}}}}}}};
}

json::object observable(const std::string& identity,
                        const std::string& checkpoint,
                        const std::string& tail = "stored",
                        bool malformed_upper = false) {
  return json::object{
      {"identity", identity}, {"checkpoint_identity", checkpoint},
      {"lower_integrand_rows", json::array{row(identity + ":lower")}},
      {"upper_integrand_rows",
       json::array{row(identity + ":upper", malformed_upper)}},
      {"epsilon", json::object{
           {"min", 0}, {"max", 0}, {"required_complete_max", 0}}},
      {"tail_policy", tail}};
}

json::object public_prefix_observable(const std::string& identity,
                                      const std::string& checkpoint) {
  return json::object{
      {"identity", identity}, {"checkpoint_identity", checkpoint},
      {"lower_integrand_rows",
       json::array{public_prefix_row(identity + ":lower")}},
      {"upper_integrand_rows",
       json::array{public_prefix_row(identity + ":upper")}},
      {"epsilon", json::object{
           {"min", 0}, {"max", 0}, {"required_complete_max", 0}}},
      {"tail_policy", "stored"}};
}

json::object regulated_public_prefix_row(
    const std::string& identity, const std::string& coefficient,
    bool omit_primitive_halo = false) {
  json::array analytic;
  analytic.push_back(json::object{
      {"numerator", json::array{coefficient}},
      {"denominator", json::array{"1"}}});
  if (!omit_primitive_halo)
    analytic.push_back(json::object{
        {"numerator", json::array{"0"}},
        {"denominator", json::array{"1"}}});
  return json::object{
      {"schema", "diffexp2-prepared-rational-local-row-v1"},
      {"columns", 1}, {"exact_identity", identity},
      {"entries", json::array{json::object{
           {"column", 0},
           {"multiplier", json::object{
                {"epsilon_shift", 0}, {"center_pole_order", 0},
                {"analytic_coefficients", std::move(analytic)},
                {"exact_identity", coefficient},
                {"proven_zero", false}}}}}}};
}

json::object regulated_public_prefix_observable(
    const std::string& identity, const std::string& checkpoint,
    bool omit_primitive_halo = false) {
  return json::object{
      {"identity", identity}, {"checkpoint_identity", checkpoint},
      {"lower_integrand_rows",
       json::array{
           regulated_public_prefix_row(
               identity + ":lower", "1", omit_primitive_halo)}},
      {"upper_integrand_rows",
       json::array{
           regulated_public_prefix_row(
               identity + ":upper", "2", omit_primitive_halo)}},
      {"epsilon", json::object{
           {"min", -1}, {"max", 0}, {"required_complete_max", 0}}},
      {"tail_policy", "stored"}};
}

json::object regulated_observable(const std::string& identity,
                                  const std::string& checkpoint) {
  return json::object{
      {"identity", identity}, {"checkpoint_identity", checkpoint},
      {"lower_integrand_rows",
       json::array{regulated_row(identity + ":lower", "1")}},
      {"upper_integrand_rows",
       json::array{regulated_row(identity + ":upper", "2")}},
      {"epsilon", json::object{
           {"min", -1}, {"max", 2}, {"required_complete_max", 0}}},
      {"tail_policy", "stored"}};
}

json::object state_reference(const json::object& state) {
  return json::object{
      {"transport_state", state.at("transport_state")},
      {"checkpoint_identity", state.at("checkpoint_identity")},
      {"provenance_identity", state.at("provenance_identity")}};
}

json::object contract_pair(const std::string& session,
                           const json::object& lower,
                           const json::object& upper,
                           json::array observables,
                           const std::string& checkpoint_root,
                           std::uint32_t threads = 0) {
  auto payload = json::object{
      {"schema", 2}, {"op", "transport.contract_pair"},
      {"session", session}, {"lower", state_reference(lower)},
      {"upper", state_reference(upper)},
      {"checkpoint_policy", json::object{
           {"schema",
            "diffexp2-deterministic-transport-pair-contraction-checkpoints-v1"},
           {"root", checkpoint_root}}},
      {"observables", std::move(observables)}};
  if (threads != 0) payload["threads"] = threads;
  return request(std::move(payload));
}

json::object stream_observable(const std::string& identity,
                               const std::string& checkpoint,
                               std::int32_t minimum = 0,
                               std::int32_t maximum = 0) {
  return json::object{
      {"identity", identity}, {"checkpoint_identity", checkpoint},
      {"epsilon", json::object{
           {"min", minimum}, {"max", maximum},
           {"required_complete_max", 0}}},
      {"tail_policy", "stored"}};
}

json::object begin_stream(const std::string& session,
                          const json::object& lower,
                          const json::object& upper,
                          const std::string& identity,
                          const std::string& checkpoint,
                          const std::string& checkpoint_root,
                          std::int32_t minimum = 0,
                          std::int32_t maximum = 0) {
  return request(json::object{
      {"schema", 2}, {"op", "transport.contract_pair_stream_begin"},
      {"session", session}, {"lower", state_reference(lower)},
      {"upper", state_reference(upper)},
      {"checkpoint_policy", json::object{
           {"schema",
            "diffexp2-deterministic-transport-pair-contraction-checkpoints-v1"},
           {"root", checkpoint_root}}},
      {"observable",
       stream_observable(identity, checkpoint, minimum, maximum)}});
}

json::object begin_rolling_stream(const std::string& session,
                                  const json::object& plan,
                                  const std::string& anchor,
                                  const std::string& identity,
                                  const std::string& checkpoint,
                                  const std::string& checkpoint_root,
                                  std::int32_t minimum = 0,
                                  std::int32_t maximum = 0) {
  return request(json::object{
      {"schema", 2},
      {"op", "transport.contract_pair_rolling_stream_begin"},
      {"session", session}, {"tile_plan", plan.at("tile_plan")},
      {"tile_plan_checkpoint_identity", plan.at("checkpoint_identity")},
      {"anchor", anchor},
      {"anchor_checkpoint_identity", "pair-common-anchor"},
      {"checkpoint_policy", json::object{
           {"schema",
            "diffexp2-deterministic-transport-pair-contraction-checkpoints-v1"},
           {"root", checkpoint_root}}},
      {"observable",
       stream_observable(identity, checkpoint, minimum, maximum)}});
}

json::object add_stream_tile(const std::string& session,
                             const json::object& stream,
                             const std::string& side,
                             std::size_t tile,
                             json::object prepared_row) {
  return request(json::object{
      {"schema", 2}, {"op", "transport.contract_pair_stream_add_tile"},
      {"session", session}, {"stream", stream.at("stream")},
      {"stream_checkpoint_identity",
       stream.at("stream_checkpoint_identity")},
      {"side", side}, {"tile", tile}, {"row", std::move(prepared_row)}});
}

json::object add_rolling_stream_tile(
    const std::string& session, const json::object& stream,
    const std::string& side, std::size_t tile,
    const std::string& local, json::object prepared_row) {
  return request(json::object{
      {"schema", 2},
      {"op", "transport.contract_pair_rolling_stream_add_tile"},
      {"session", session}, {"stream", stream.at("stream")},
      {"stream_checkpoint_identity",
       stream.at("stream_checkpoint_identity")},
      {"side", side}, {"tile", tile}, {"local", local},
      {"local_checkpoint_identity", "pair-common-anchor"},
      {"row", std::move(prepared_row)}});
}

json::object finish_stream(const std::string& session,
                           const json::object& stream) {
  return request(json::object{
      {"schema", 2}, {"op", "transport.contract_pair_stream_finish"},
      {"session", session}, {"stream", stream.at("stream")},
      {"stream_checkpoint_identity",
       stream.at("stream_checkpoint_identity")}});
}

json::object abort_stream(const std::string& session,
                          const json::object& stream) {
  return request(json::object{
      {"schema", 2}, {"op", "transport.contract_pair_stream_abort"},
      {"session", session}, {"stream", stream.at("stream")},
      {"stream_checkpoint_identity",
       stream.at("stream_checkpoint_identity")}});
}

json::object session_stats(const std::string& session) {
  return request(json::object{{"schema", 2}, {"op", "session.stats"},
                              {"session", session}});
}

json::object state_stats(const std::string& session,
                         const json::object& state) {
  return request(json::object{
      {"schema", 2}, {"op", "transport.stats"}, {"session", session},
      {"transport_state", state.at("transport_state")}});
}

json::object export_line(const std::string& session,
                         const json::object& line) {
  return request(json::object{
      {"schema", 2}, {"op", "integration.export"},
      {"session", session}, {"line", line.at("line")},
      {"checkpoint_identity", line.at("checkpoint_identity")},
      {"output_digits", 50}});
}

json::object rolling_endpoint_batch(
    const std::string& session, const json::object& plan,
    const std::string& arm, const std::string& local,
    const std::string& identity, const std::string& checkpoint,
    const std::string& checkpoint_root) {
  return request(json::object{
      {"schema", 2}, {"op", "transport.rolling_endpoint_batch"},
      {"session", session}, {"tile_plan", plan.at("tile_plan")},
      {"tile_plan_checkpoint_identity", plan.at("checkpoint_identity")},
      {"arm", arm}, {"local", local},
      {"local_checkpoint_identity", "pair-common-anchor"},
      {"checkpoint_policy", json::object{
           {"schema",
            "diffexp2-deterministic-transport-endpoint-checkpoints-v1"},
           {"root", checkpoint_root}}},
      {"observables", json::array{json::object{
           {"identity", identity}, {"checkpoint_identity", checkpoint},
           {"integrand_row", row(identity + ":row")},
           {"publication_relative_tolerance", "1e-30"},
           {"epsilon", json::object{
                {"min", 0}, {"max", 0},
                {"required_complete_max", 0}}}}}},
      {"threads", 1}});
}

json::object export_endpoint(const std::string& session,
                             const json::object& endpoint) {
  return request(json::object{
      {"schema", 2}, {"op", "endpoint.export"}, {"session", session},
      {"endpoint", endpoint.at("endpoint")},
      {"checkpoint_identity", endpoint.at("checkpoint_identity")},
      {"output_digits", 50}});
}

double midpoint(const json::object& exported) {
  const auto& coefficient = exported.at("value").as_object()
      .at("coefficients").as_array().front().as_array();
  return std::stod(std::string(coefficient.front().as_string()));
}

bool valid_pair_conditioning(const json::object& line,
                             std::int32_t epsilon_min,
                             std::int32_t epsilon_max) {
  const auto* raw = line.if_contains("conditioning");
  if (raw == nullptr || !raw->is_object()) return false;
  const auto& diagnostics = raw->as_object();
  if (diagnostics.at("schema") !=
          "diffexp2-transport-pair-conditioning-diagnostics-v1" ||
      diagnostics.at("dimension") != 1)
    return false;
  const auto& epsilon = diagnostics.at("epsilon").as_object();
  if (epsilon.at("min") != epsilon_min ||
      epsilon.at("max") != epsilon_max)
    return false;
  const auto& entries = diagnostics.at("entries").as_array();
  if (entries.size() != static_cast<std::size_t>(
                            epsilon_max - epsilon_min + 1))
    return false;
  for (std::size_t index = 0; index < entries.size(); ++index) {
    const auto& entry = entries[index].as_object();
    if (entry.at("power") !=
            epsilon_min + static_cast<std::int32_t>(index) ||
        entry.at("component") != 0 ||
        !entry.at("lower_midpoint").is_array() ||
        !entry.at("upper_midpoint").is_array() ||
        !entry.at("combined_midpoint").is_array() ||
        !entry.at("lower_radius2exp").is_array() ||
        !entry.at("upper_radius2exp").is_array() ||
        !entry.at("combined_radius2exp").is_array() ||
        !entry.at("combined_contains_zero").is_bool())
      return false;
  }
  for (const auto* key : {"lower_tiles", "upper_tiles"}) {
    const auto& tiles = diagnostics.at(key).as_array();
    if (tiles.empty()) return false;
    for (const auto& raw_tile : tiles) {
      const auto& tile = raw_tile.as_object();
      const auto& value = tile.at("value").as_object();
      for (const auto& raw_entry : value.at("entries").as_array())
        if (!raw_entry.as_object().at("radius2exp").is_array())
          return false;
    }
  }
  return true;
}

bool payload_has_handle(const json::array& records,
                        const std::string& handle) {
  for (const auto& raw : records)
    if (raw.as_object().at("handle").as_string() == handle) return true;
  return false;
}

void regulated_center_primitive_preserves_pair_lower_halo() {
  std::string session;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "rational"}, {"output_digits", 50},
        {"chart_capacity", 1}, {"local_capacity", 1},
        {"match_capacity", 1}, {"tile_plan_capacity", 2},
        {"transport_state_capacity", 2}, {"line_result_capacity", 2}});
    if (created.at("status") != "ok")
      throw std::runtime_error(
          "regulated session.create: " + json::serialize(created));
    session = std::string(created.at("session").as_string());
    const auto chart = prepare_anchor_chart(session);
    const auto anchor = solve_regulated_anchor(session, chart);
    const auto lower_plan = prepare_plan(
        session, "regulated-lower-plan", "-1/3", chart);
    const auto upper_plan = prepare_plan(
        session, "regulated-upper-plan", "1/2", chart);
    const auto lower = run_state(
        session, lower_plan, "regulated-lower-plan", anchor, "lower",
        "regulated-lower-state");
    const auto upper = run_state(
        session, upper_plan, "regulated-upper-plan", anchor, "upper",
        "regulated-upper-state");
    if (lower.at("status") != "ok" || upper.at("status") != "ok")
      throw std::runtime_error(
          "regulated retained states failed: " + json::serialize(lower) +
          " / " + json::serialize(upper));

    const auto paired = contract_pair(
        session, lower, upper,
        json::array{regulated_observable(
            "pair-regulated-halo", "pair-regulated-halo-checkpoint")},
        "pair-regulated-halo-root");
    if (paired.at("status") != "ok" ||
        paired.at("lines").as_array().size() != 1)
      throw std::runtime_error(
          "regulated pair contraction failed: " + json::serialize(paired));
    const auto& paired_line =
        paired.at("lines").as_array().front().as_object();
    const auto paired_export = export_line(session, paired_line);
    if (paired_line.at("epsilon").as_object().at("min") != -1 ||
        paired_line.at("epsilon").as_object().at("max") != 1 ||
        !valid_pair_conditioning(paired_line, -1, 1) ||
        paired_export.at("status") != "ok" ||
        paired_export.at("value").as_object().at("min") != -1 ||
        paired_export.at("value").as_object().at("max") != 1 ||
        std::abs(midpoint(paired_export) - 1.0) > 1e-40)
      throw std::runtime_error(
          "one-shot pair clipped the regulated primitive lower halo: " +
          json::serialize(paired_export));

    // The source owns three epsilon rows, while a public epsilon^0 result
    // behind a genuine 1/epsilon primitive needs exactly two multiplier
    // rows.  This exercises the consumer halo, not merely an ordinary
    // primitive whose epsilon valuation is zero.
    const auto public_prefix = contract_pair(
        session, lower, upper,
        json::array{regulated_public_prefix_observable(
            "pair-regulated-public-prefix",
            "pair-regulated-public-prefix-checkpoint")},
        "pair-regulated-public-prefix-root");
    if (public_prefix.at("status") != "ok" ||
        public_prefix.at("lines").as_array().size() != 1)
      throw std::runtime_error(
          "regulated public-prefix pair contraction failed: " +
          json::serialize(public_prefix));
    const auto& public_prefix_line =
        public_prefix.at("lines").as_array().front().as_object();
    const auto public_prefix_export =
        export_line(session, public_prefix_line);
    const auto& full_coefficients =
        paired_export.at("value").as_object().at("coefficients").as_array();
    const auto& prefix_coefficients =
        public_prefix_export.at("value").as_object()
            .at("coefficients").as_array();
    if (public_prefix_line.at("epsilon").as_object().at("min") != -1 ||
        public_prefix_line.at("epsilon").as_object().at("max") != 0 ||
        public_prefix_export.at("status") != "ok" ||
        public_prefix_export.at("value").as_object().at("min") != -1 ||
        public_prefix_export.at("value").as_object().at("max") != 0 ||
        prefix_coefficients.size() != 2 ||
        prefix_coefficients.at(0) != full_coefficients.at(0) ||
        prefix_coefficients.at(1) != full_coefficients.at(1))
      throw std::runtime_error(
          "regulated public-prefix row lost its exact primitive halo: " +
          json::serialize(public_prefix_export));

    const auto undercovered = contract_pair(
        session, lower, upper,
        json::array{regulated_public_prefix_observable(
            "pair-regulated-undercovered",
            "pair-regulated-undercovered-checkpoint", true)},
        "pair-regulated-undercovered-root");
    if (undercovered.at("status") != "error" ||
        std::string(undercovered.at("detail").as_string()).find(
            "consumer prefix is too short") == std::string::npos)
      throw std::runtime_error(
          "regulated row without the primitive halo was not rejected: " +
          json::serialize(undercovered));

    const auto stream = begin_stream(
        session, lower, upper, "pair-stream-regulated-halo",
        "pair-stream-regulated-halo-checkpoint",
        "pair-stream-regulated-halo-root", -1, 2);
    const auto lower_added = add_stream_tile(
        session, stream, "lower", 0,
        regulated_row("pair-stream-regulated-halo:lower", "1"));
    const auto upper_added = add_stream_tile(
        session, stream, "upper", 0,
        regulated_row("pair-stream-regulated-halo:upper", "2"));
    const auto streamed = finish_stream(session, stream);
    if (stream.at("status") != "ok" ||
        lower_added.at("status") != "ok" ||
        upper_added.at("status") != "ok" ||
        streamed.at("status") != "ok" ||
        streamed.at("lines").as_array().size() != 1)
      throw std::runtime_error(
          "regulated pair stream failed: " + json::serialize(streamed));
    const auto& streamed_line =
        streamed.at("lines").as_array().front().as_object();
    const auto streamed_export = export_line(session, streamed_line);
    if (streamed_line.at("epsilon").as_object().at("min") != -1 ||
        streamed_line.at("epsilon").as_object().at("max") != 1 ||
        !valid_pair_conditioning(streamed_line, -1, 1) ||
        streamed_line.at("conditioning") !=
            paired_line.at("conditioning") ||
        streamed_export.at("status") != "ok" ||
        streamed_export.at("value").as_object().at("min") != -1 ||
        streamed_export.at("value").as_object().at("max") != 1 ||
        std::abs(midpoint(streamed_export) - 1.0) > 1e-40 ||
        streamed_export.at("value") != paired_export.at("value"))
      throw std::runtime_error(
          "streamed pair clipped or changed the regulated primitive lower halo: " +
          json::serialize(streamed_export));

    for (const auto& state : {lower, upper})
      (void)request(json::object{
          {"schema", 2}, {"op", "transport.release"},
          {"session", session},
          {"transport_state", state.at("transport_state")}});
    for (const auto& plan : {lower_plan, upper_plan})
      (void)request(json::object{
          {"schema", 2}, {"op", "tile.release"}, {"session", session},
          {"tile_plan", plan.at("tile_plan")}});
    (void)request(json::object{{"schema", 2}, {"op", "local.release"},
                               {"session", session}, {"local", anchor}});
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    session.clear();
  } catch (...) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    throw;
  }
}

void rolling_stream_contracts_before_state_publication() {
  const std::string checkpoint =
      "/tmp/diffexp2-rolling-transport-contract-pair.de2cp";
  std::remove(checkpoint.c_str());
  std::string session;
  std::string restored_session;
  std::string stage = "session.create";
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "rational"}, {"output_digits", 50},
        {"chart_capacity", 1}, {"local_capacity", 1},
        {"match_capacity", 1}, {"tile_plan_capacity", 1},
        {"transport_state_capacity", 1}, {"line_result_capacity", 2},
        {"endpoint_capacity", 2}});
    if (created.at("status") != "ok")
      throw std::runtime_error(
          "rolling session.create: " + json::serialize(created));
    session = std::string(created.at("session").as_string());
    stage = "prepare owners";
    const auto chart = prepare_anchor_chart(session);
    const auto anchor = solve_anchor(session, chart);
    const auto plan = prepare_pair_plan(
        session, "rolling-pair-plan", chart);
    const auto stream = begin_rolling_stream(
        session, plan, anchor, "rolling-pair-observable",
        "rolling-pair-observable-checkpoint", "rolling-pair-root");
    stage = "validate begin";
    if (plan.at("status") != "ok" || stream.at("status") != "ok" ||
        stream.at("capability") !=
            "rolling-transport-pair-observable-stream-v1" ||
        stream.at("lower_tiles") != 1 || stream.at("upper_tiles") != 1 ||
        stream.at("coefficient_retention") !=
            "current-and-terminal-only")
      throw std::runtime_error(
          "rolling stream did not bind the paired plan exactly: " +
          json::serialize(stream));

    stage = "add tiles";
    const auto lower = add_rolling_stream_tile(
        session, stream, "lower", 0, anchor,
        row("rolling-pair:lower"));
    const auto upper = add_rolling_stream_tile(
        session, stream, "upper", 0, anchor,
        row("rolling-pair:upper"));
    if (counter(lower, "retained_row_record_bytes") == 0 ||
        counter(lower, "retained_tile_diagnostic_bytes") == 0 ||
        counter(upper, "retained_row_record_bytes") <=
            counter(lower, "retained_row_record_bytes") ||
        counter(upper, "retained_tile_diagnostic_bytes") <=
            counter(lower, "retained_tile_diagnostic_bytes"))
      throw std::runtime_error(
          "rolling stream did not report monotonically retained compact byte classes");
    stage = "finish";
    const auto finished = finish_stream(session, stream);
    stage = "validate finish";
    if (lower.at("status") != "ok" || upper.at("status") != "ok" ||
        upper.at("next_side") != nullptr ||
        finished.at("status") != "ok" ||
        finished.at("capability") !=
            "rolling-transport-pair-observable-stream-v1" ||
        finished.at("tile_integrations") != 2 ||
        finished.at("lines").as_array().size() != 1 ||
        finished.at("lower").as_object().at("coefficient_retention") !=
            "terminal-local-only" ||
        finished.at("upper").as_object().at("coefficient_retention") !=
            "terminal-local-only")
      throw std::runtime_error(
          "rolling stream failed to publish a compact pair: " +
          json::serialize(finished));
    const auto line =
        finished.at("lines").as_array().front().as_object();
    const auto& conditioning = line.at("conditioning").as_object();
    for (const auto* side : {"lower_tiles", "upper_tiles"}) {
      const auto& compact = conditioning.at(side).as_array().front()
          .as_object().at("value").as_object();
      if (compact.at("schema") !=
              "diffexp2-compact-rolling-tile-diagnostic-v1" ||
          compact.if_contains("entries") != nullptr ||
          compact.at("coefficients") != 1)
        throw std::runtime_error(
            "rolling line retained expanded per-coefficient tile diagnostics");
    }
    stage = "export";
    const auto exported = export_line(session, line);
    if (exported.at("status") != "ok" ||
        std::abs(midpoint(exported) - 5.0 / 6.0) > 1e-40)
      throw std::runtime_error(
          "rolling stream changed the exact paired value: " +
          json::serialize(exported));

    stage = "rolling endpoint";
    const auto endpoint_batch = rolling_endpoint_batch(
        session, plan, "lower", anchor, "rolling-lower-endpoint",
        "rolling-lower-endpoint-checkpoint",
        "rolling-lower-endpoint-root");
    if (endpoint_batch.at("status") != "ok" ||
        endpoint_batch.at("capability") !=
            "rolling-transport-endpoint-batch-v1" ||
        endpoint_batch.at("coefficient_retention") !=
            "terminal-local-only" ||
        endpoint_batch.at("endpoints").as_array().size() != 1)
      throw std::runtime_error(
          "rolling terminal-local endpoint batch failed: " +
          json::serialize(endpoint_batch));
    const auto endpoint =
        endpoint_batch.at("endpoints").as_array().front().as_object();
    const auto endpoint_export = export_endpoint(session, endpoint);
    if (endpoint_export.at("status") != "ok" ||
        std::abs(midpoint(endpoint_export) - 1.0) > 1e-40)
      throw std::runtime_error(
          "rolling terminal-local endpoint changed its exact value: " +
          json::serialize(endpoint_export));

    stage = "release owners";
    (void)request(json::object{
        {"schema", 2}, {"op", "tile.release"}, {"session", session},
        {"tile_plan", plan.at("tile_plan")}});
    (void)request(json::object{
        {"schema", 2}, {"op", "local.release"}, {"session", session},
        {"local", anchor}});
    if (export_line(session, line).at("status") != "ok")
      throw std::runtime_error(
          "rolling line did not retain its compact plan/local closure");
    if (export_endpoint(session, endpoint).at("status") != "ok")
      throw std::runtime_error(
          "rolling endpoint did not retain its compact plan/local closure");

    stage = "checkpoint save";
    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
        {"path", checkpoint},
        {"checkpoint_identity", "rolling-pair-roundtrip"}});
    if (saved.at("status") != "ok")
      throw std::runtime_error(
          "rolling checkpoint.save: " + json::serialize(saved));
    (void)request(json::object{
        {"schema", 2}, {"op", "session.close"}, {"session", session}});
    session.clear();

    stage = "checkpoint restore";
    const auto restored = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint},
        {"expected_identity", "rolling-pair-roundtrip"}});
    if (restored.at("status") != "ok")
      throw std::runtime_error(
          "rolling checkpoint.restore: " + json::serialize(restored));
    restored_session = std::string(restored.at("session").as_string());
    const auto restored_export = export_line(restored_session, line);
    if (restored_export.at("status") != "ok" ||
        restored_export.at("value") != exported.at("value"))
      throw std::runtime_error(
          "rolling line changed across checkpoint restore: " +
          json::serialize(restored_export));
    const auto restored_endpoint_export =
        export_endpoint(restored_session, endpoint);
    if (restored_endpoint_export.at("status") != "ok" ||
        restored_endpoint_export.at("value") !=
            endpoint_export.at("value"))
      throw std::runtime_error(
          "rolling endpoint changed across checkpoint restore: " +
          json::serialize(restored_endpoint_export));
    (void)request(json::object{
        {"schema", 2}, {"op", "session.close"},
        {"session", restored_session}});
    restored_session.clear();
    std::remove(checkpoint.c_str());
  } catch (const std::exception& error) {
    if (!session.empty())
      (void)request(json::object{
          {"schema", 2}, {"op", "session.close"}, {"session", session}});
    if (!restored_session.empty())
      (void)request(json::object{
          {"schema", 2}, {"op", "session.close"},
          {"session", restored_session}});
    std::remove(checkpoint.c_str());
    throw std::runtime_error(
        "rolling stream stage " + stage + ": " + error.what());
  } catch (...) {
    if (!session.empty())
      (void)request(json::object{
          {"schema", 2}, {"op", "session.close"}, {"session", session}});
    if (!restored_session.empty())
      (void)request(json::object{
          {"schema", 2}, {"op", "session.close"},
          {"session", restored_session}});
    std::remove(checkpoint.c_str());
    throw;
  }
}

}  // namespace

int main() {
  const std::string checkpoint =
      "/tmp/diffexp2-transport-contract-pair.de2cp";
  const std::string checkpoint_second =
      "/tmp/diffexp2-transport-contract-pair-second.de2cp";
  const std::string checkpoint_tampered =
      "/tmp/diffexp2-transport-contract-pair-tampered.de2cp";
  const std::string checkpoint_count_tampered =
      "/tmp/diffexp2-transport-contract-pair-count-tampered.de2cp";
  const std::string checkpoint_legacy =
      "/tmp/diffexp2-transport-contract-pair-legacy.de2cp";
  std::remove(checkpoint.c_str());
  std::remove(checkpoint_second.c_str());
  std::remove(checkpoint_tampered.c_str());
  std::remove(checkpoint_count_tampered.c_str());
  std::remove(checkpoint_legacy.c_str());
  std::string session;
  std::string restored_session;
  try {
    regulated_center_primitive_preserves_pair_lower_halo();
    rolling_stream_contracts_before_state_publication();
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "rational"}, {"output_digits", 50},
        {"chart_capacity", 2}, {"local_capacity", 4},
        {"match_capacity", 2}, {"tile_plan_capacity", 2},
        {"transport_state_capacity", 2}, {"line_result_capacity", 8}});
    session = std::string(created.at("session").as_string());
    const auto chart = prepare_anchor_chart(session);
    const auto anchor = solve_anchor(session, chart);
    const auto lower_plan = prepare_plan(
        session, "pair-lower-plan", "-1/3", chart);
    const auto upper_plan = prepare_plan(
        session, "pair-upper-plan", "1/2", chart);
    if (lower_plan.at("status") != "ok" ||
        upper_plan.at("status") != "ok" ||
        lower_plan.at("tile_plan") == upper_plan.at("tile_plan") ||
        lower_plan.at("arm_name") != "lower" ||
        upper_plan.at("arm_name") != "upper")
      throw std::runtime_error(
          "independent single-arm plans were not retained exactly");
    const auto lower = run_state(
        session, lower_plan, "pair-lower-plan", anchor, "lower",
        "pair-lower-state");
    const auto upper = run_state(
        session, upper_plan, "pair-upper-plan", anchor, "upper",
        "pair-upper-state");
    if (lower.at("status") != "ok" || upper.at("status") != "ok" ||
        lower.at("tiles") != 1 || upper.at("tiles") != 1 ||
        lower.at("tile_source_epsilon") !=
            json::array{json::object{{"min", 0}, {"max", 2}}} ||
        upper.at("tile_source_epsilon") !=
            json::array{json::object{{"min", 0}, {"max", 2}}})
      throw std::runtime_error(
          "paired retained states did not expose exact live tile-source epsilon windows: " +
          json::serialize(lower) + " / " + json::serialize(upper));

    const auto before = session_stats(session);
    const auto zero = contract_pair(
        session, lower, upper, json::array{}, "pair-zero");
    const auto after_zero = session_stats(session);
    const auto lower_after_zero = state_stats(session, lower);
    const auto upper_after_zero = state_stats(session, upper);
    if (zero.at("status") != "ok" || zero.at("observables") != 0 ||
        !zero.at("lines").as_array().empty() ||
        zero.at("max_parallel_arms") != 0 ||
        before.at("line_results") != after_zero.at("line_results") ||
        before.at("line_integrations") != after_zero.at("line_integrations") ||
        before.at("locals") != after_zero.at("locals") ||
        before.at("matches") != after_zero.at("matches") ||
        counter(after_zero, "transport_contractions") !=
            counter(before, "transport_contractions") + 2 ||
        counter(after_zero, "transport_pair_contractions") !=
            counter(before, "transport_pair_contractions") + 1 ||
        counter(lower_after_zero, "contraction_operations") != 1 ||
        counter(upper_after_zero, "contraction_operations") != 1)
      throw std::runtime_error(
          "zero-observable pair contraction changed allocation or counters incorrectly");

    const auto one = contract_pair(
        session, lower, upper,
        json::array{observable("pair-one", "pair-one-checkpoint")},
        "pair-one-root");
    if (one.at("status") != "ok" || one.at("observables") != 1 ||
        one.at("lines").as_array().size() != 1 ||
        one.at("combination") != "negative-lower-plus-upper" ||
        one.at("max_parallel_arms") != 1 ||
        one.at("concurrent_arms") != false ||
        one.at("streaming_tile_contraction") != true ||
        one.at("no_remarching") != true || one.at("no_rematching") != true)
      throw std::runtime_error(
          "one-observable pair contraction failed: " + json::serialize(one));
    if (!valid_pair_conditioning(
            one.at("lines").as_array().front().as_object(), 0, 0))
      throw std::runtime_error(
          "one-observable pair contraction omitted its bounded conditioning diagnostics");
    const auto one_export = export_line(
        session, one.at("lines").as_array().front().as_object());
    if (one_export.at("status") != "ok" ||
        std::abs(midpoint(one_export) - 5.0 / 6.0) > 1e-40)
      throw std::runtime_error(
          "fixed -lower+upper value is incorrect: " +
          json::serialize(one_export));

    auto direct_publication_observable = observable(
        "pair-direct-publication",
        "pair-direct-publication-checkpoint");
    direct_publication_observable[
        "publication_relative_tolerance"] = "1e-20";
    const auto direct_publication = contract_pair(
        session, lower, upper,
        json::array{std::move(direct_publication_observable)},
        "pair-direct-publication-root");
    if (direct_publication.at("status") != "ok" ||
        direct_publication.at("ordinary_factorization_retries") != 0 ||
        direct_publication.at("lines").as_array().size() != 1)
      throw std::runtime_error(
          "a rigorous direct paired line that met publication was "
          "unnecessarily factorized: " +
          json::serialize(direct_publication));

    // The retained arm source owns three epsilon rows, but a public
    // epsilon^0 observable needs only the first multiplier row.  Private
    // source reservoirs used by matching must not force an otherwise
    // sufficient public integrand row to reproduce those proof-only orders.
    const auto public_prefix = contract_pair(
        session, lower, upper,
        json::array{public_prefix_observable(
            "pair-public-prefix", "pair-public-prefix-checkpoint")},
        "pair-public-prefix-root");
    if (public_prefix.at("status") != "ok" ||
        public_prefix.at("lines").as_array().size() != 1)
      throw std::runtime_error(
          "public-prefix row did not ignore the private source reservoir: " +
          json::serialize(public_prefix));
    const auto public_prefix_export = export_line(
        session,
        public_prefix.at("lines").as_array().front().as_object());
    if (public_prefix_export.at("status") != "ok" ||
        public_prefix_export.at("value").as_object().at("min") != 0 ||
        public_prefix_export.at("value").as_object().at("max") != 0 ||
        std::abs(midpoint(public_prefix_export) - 5.0 / 6.0) > 1e-40)
      throw std::runtime_error(
          "public-prefix row changed the contracted public coefficient: " +
          json::serialize(public_prefix_export));

    json::array many_observables;
    for (int index = 0; index < 3; ++index)
      many_observables.push_back(observable(
          "pair-many-" + std::to_string(index),
          "pair-many-checkpoint-" + std::to_string(index)));
    const auto many = contract_pair(
        session, lower, upper, std::move(many_observables),
        "pair-many-root", 2);
    std::set<std::string> many_handles;
    if (many.at("status") != "ok" ||
        many.at("lines").as_array().size() != 3 ||
        many.at("requested_observable_threads") != 2 ||
        many.at("observable_worker_threads") != 2)
      throw std::runtime_error(
          "many-observable pair contraction failed: " +
          json::serialize(many));
    for (std::size_t index = 0; index < 3; ++index) {
      const auto& line = many.at("lines").as_array()[index].as_object();
      if (counter(line, "request_index") != index ||
          std::string(line.at("observable_identity").as_string()) !=
              "pair-many-" + std::to_string(index) ||
          !many_handles.insert(std::string(line.at("line").as_string())).second)
        throw std::runtime_error(
            "many-observable pair order or handles changed");
    }

    const auto before_malformed = session_stats(session);
    const auto lower_before_malformed = state_stats(session, lower);
    const auto upper_before_malformed = state_stats(session, upper);
    const auto malformed = contract_pair(
        session, lower, upper,
        json::array{
            observable("pair-valid-prefix", "pair-valid-prefix-checkpoint"),
            observable("pair-malformed", "pair-malformed-checkpoint",
                       "stored", true)},
        "pair-malformed-root");
    const auto after_malformed = session_stats(session);
    const auto lower_after_malformed = state_stats(session, lower);
    const auto upper_after_malformed = state_stats(session, upper);
    if (malformed.at("status") != "error" ||
        before_malformed.at("line_results") !=
            after_malformed.at("line_results") ||
        before_malformed.at("line_integrations") !=
            after_malformed.at("line_integrations") ||
        before_malformed.at("transport_pair_contractions") !=
            after_malformed.at("transport_pair_contractions") ||
        lower_before_malformed.at("contraction_operations") !=
            lower_after_malformed.at("contraction_operations") ||
        upper_before_malformed.at("contraction_operations") !=
            upper_after_malformed.at("contraction_operations") ||
        after_malformed.at("pending_line_integrations") != 0)
      throw std::runtime_error(
          "malformed upper arm did not roll back atomically: " +
          json::serialize(malformed));

    const auto before_require = session_stats(session);
    const auto require = contract_pair(
        session, lower, upper,
        json::array{observable("pair-require", "pair-require-checkpoint",
                               "require")},
        "pair-require-root");
    const auto after_require = session_stats(session);
    if (require.at("status") != "error" ||
        before_require.at("line_results") != after_require.at("line_results") ||
        before_require.at("line_integrations") !=
            after_require.at("line_integrations") ||
        before_require.at("transport_pair_contractions") !=
            after_require.at("transport_pair_contractions") ||
        after_require.at("pending_line_integrations") != 0)
      throw std::runtime_error(
          "required pair tail did not fail atomically: " +
          json::serialize(require));

    const auto before_stream_abort = session_stats(session);
    const auto aborted_stream = begin_stream(
        session, lower, upper, "pair-stream-abort",
        "pair-stream-abort-checkpoint", "pair-stream-abort-root");
    const auto open_stream_stats = session_stats(session);
    const auto blocked_checkpoint = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
        {"path", checkpoint},
        {"checkpoint_identity", "pair-open-stream-must-not-save"}});
    const auto upper_first = add_stream_tile(
        session, aborted_stream, "upper", 0,
        row("pair-stream-abort:upper"));
    const auto aborted = abort_stream(session, aborted_stream);
    const auto aborted_again = abort_stream(session, aborted_stream);
    const auto after_stream_abort = session_stats(session);
    if (aborted_stream.at("status") != "ok" ||
        open_stream_stats.at("transport_pair_streams") != 1 ||
        open_stream_stats.at("pending_line_integrations") != 1 ||
        blocked_checkpoint.at("status") != "error" ||
        upper_first.at("status") != "error" ||
        aborted.at("status") != "ok" || aborted.at("aborted") != true ||
        aborted_again.at("status") != "ok" ||
        aborted_again.at("already_absent") != true ||
        after_stream_abort.at("transport_pair_streams") != 0 ||
        after_stream_abort.at("pending_line_integrations") != 0 ||
        before_stream_abort.at("line_results") !=
            after_stream_abort.at("line_results") ||
        before_stream_abort.at("transport_pair_contractions") !=
            after_stream_abort.at("transport_pair_contractions") ||
        before_stream_abort.at("line_integrations") !=
            after_stream_abort.at("line_integrations"))
      throw std::runtime_error(
          "stream abort/checkpoint transaction was not atomic");

    const auto before_poison = session_stats(session);
    const auto poisoned_stream = begin_stream(
        session, lower, upper, "pair-stream-poison",
        "pair-stream-poison-checkpoint", "pair-stream-poison-root");
    const auto malformed_tile = add_stream_tile(
        session, poisoned_stream, "lower", 0,
        row("pair-stream-poison:malformed", true));
    const auto poisoned_retry = add_stream_tile(
        session, poisoned_stream, "lower", 0,
        row("pair-stream-poison:retry"));
    const auto poison_abort = abort_stream(session, poisoned_stream);
    const auto after_poison = session_stats(session);
    if (malformed_tile.at("status") != "error" ||
        poisoned_retry.at("status") != "error" ||
        poison_abort.at("status") != "ok" ||
        after_poison.at("transport_pair_streams") != 0 ||
        after_poison.at("pending_line_integrations") != 0 ||
        before_poison.at("line_results") != after_poison.at("line_results") ||
        before_poison.at("transport_pair_contractions") !=
            after_poison.at("transport_pair_contractions") ||
        before_poison.at("line_integrations") !=
            after_poison.at("line_integrations"))
      throw std::runtime_error(
          "malformed streamed row did not poison and abort atomically");

    const auto before_premature = session_stats(session);
    const auto premature_stream = begin_stream(
        session, lower, upper, "pair-stream-premature",
        "pair-stream-premature-checkpoint", "pair-stream-premature-root");
    const auto premature_finish = finish_stream(session, premature_stream);
    const auto after_premature = session_stats(session);
    if (premature_finish.at("status") != "error" ||
        after_premature.at("transport_pair_streams") != 0 ||
        after_premature.at("pending_line_integrations") != 0 ||
        before_premature.at("line_results") !=
            after_premature.at("line_results") ||
        before_premature.at("transport_pair_contractions") !=
            after_premature.at("transport_pair_contractions") ||
        before_premature.at("line_integrations") !=
            after_premature.at("line_integrations"))
      throw std::runtime_error(
          "premature streamed finish did not terminate without publication");

    const auto streamed_begin = begin_stream(
        session, lower, upper, "pair-streamed",
        "pair-streamed-checkpoint", "pair-streamed-root");
    const auto streamed_lower = add_stream_tile(
        session, streamed_begin, "lower", 0,
        row("pair-streamed:lower"));
    const auto streamed_upper = add_stream_tile(
        session, streamed_begin, "upper", 0,
        row("pair-streamed:upper"));
    const auto streamed = finish_stream(session, streamed_begin);
    if (streamed_begin.at("status") != "ok" ||
        streamed_lower.at("status") != "ok" ||
        streamed_upper.at("status") != "ok" ||
        streamed_upper.at("next_side") != nullptr ||
        streamed.at("status") != "ok" ||
        streamed.at("persistent_tile_stream") != true ||
        streamed.at("lines").as_array().size() != 1 ||
        streamed.at("tile_integrations") != 2)
      throw std::runtime_error(
          "one-row transport-pair stream did not finish exactly: " +
          json::serialize(streamed));
    if (!valid_pair_conditioning(
            streamed.at("lines").as_array().front().as_object(), 0, 0) ||
        streamed.at("lines").as_array().front().as_object().at(
            "conditioning") !=
            one.at("lines").as_array().front().as_object().at(
                "conditioning"))
      throw std::runtime_error(
          "streamed pair conditioning diagnostics differ from one-shot contraction");
    const auto streamed_export = export_line(
        session, streamed.at("lines").as_array().front().as_object());
    if (streamed_export.at("status") != "ok" ||
        std::abs(midpoint(streamed_export) - 5.0 / 6.0) > 1e-40)
      throw std::runtime_error(
          "streamed pair value differs from the one-shot contraction");

    const auto after_success = session_stats(session);
    const auto lower_after_success = state_stats(session, lower);
    const auto upper_after_success = state_stats(session, upper);
    if (counter(after_success, "transport_pair_contractions") != 6 ||
        counter(after_success, "transport_pair_observables") != 7 ||
        counter(after_success, "transport_contractions") != 12 ||
        counter(after_success, "transport_observables") != 14 ||
        counter(after_success, "line_integrations") != 14 ||
        counter(lower_after_success, "contraction_operations") != 6 ||
        counter(upper_after_success, "contraction_operations") != 6 ||
        counter(lower_after_success, "contracted_observables") != 7 ||
        counter(upper_after_success, "contracted_observables") != 7)
      throw std::runtime_error(
          "successful pair counters are not additive and honest");

    std::vector<json::object> lines;
    lines.push_back(one.at("lines").as_array().front().as_object());
    lines.push_back(
        direct_publication.at("lines").as_array().front().as_object());
    for (const auto& raw : many.at("lines").as_array())
      lines.push_back(raw.as_object());
    lines.push_back(streamed.at("lines").as_array().front().as_object());
    const auto first_export_value = one_export.at("value");
    for (const auto& state : {lower, upper})
      if (request(json::object{
              {"schema", 2}, {"op", "transport.release"},
              {"session", session},
              {"transport_state", state.at("transport_state")}})
              .at("status") != "ok")
        throw std::runtime_error("transport state release failed");
    for (const auto& plan : {lower_plan, upper_plan})
      if (request(json::object{
              {"schema", 2}, {"op", "tile.release"},
              {"session", session}, {"tile_plan", plan.at("tile_plan")}})
              .at("status") != "ok")
        throw std::runtime_error("tile plan release failed");
    if (request(json::object{
            {"schema", 2}, {"op", "local.release"},
            {"session", session}, {"local", anchor}}).at("status") != "ok")
      throw std::runtime_error("common anchor release failed");
    for (const auto& line : lines)
      if (export_line(session, line).at("status") != "ok")
        throw std::runtime_error(
            "paired line did not survive public owner release");

    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
        {"path", checkpoint}, {"checkpoint_identity", "pair-roundtrip"}});
    if (saved.at("status") != "ok")
      throw std::runtime_error("checkpoint.save: " + json::serialize(saved));
    const auto container = diffexp2::checkpoint::read(checkpoint);
    const auto payload = json::parse(container.payload_json).as_object();
    if (payload.at("schema") != 9 ||
        payload.at("retained_transport_states").as_array().size() != 2 ||
        !payload_has_handle(payload.at("retained_transport_states").as_array(),
                            std::string(lower.at("transport_state").as_string())) ||
        !payload_has_handle(payload.at("retained_transport_states").as_array(),
                            std::string(upper.at("transport_state").as_string())) ||
        container.payload_json.find("private:") != std::string::npos)
      throw std::runtime_error(
          "paired checkpoint retained scratch objects or lost its two-state closure");
    for (const auto& raw : payload.at("retained_line_results").as_array()) {
      const auto& record = raw.as_object();
      if (record.at("schema") !=
              "diffexp2-retained-transport-pair-observable-line-v2" ||
          record.at("provenance").as_object().at("aggregate").as_object()
                  .if_contains("components") != nullptr)
        throw std::runtime_error(
            "paired checkpoint emitted a noncompact line record");
      const auto& provenance = record.at("provenance").as_object();
      const auto& source = provenance.at("source").as_object();
      const auto& aggregate = provenance.at("aggregate").as_object();
      if (aggregate.at("projection_mode") !=
          "fused-stored-hash-monomial-stream-v2")
        throw std::runtime_error(
            "paired checkpoint lost its bounded-memory projection algorithm");
      for (const auto* side_name : {"lower", "upper"}) {
        const auto& side = source.at(side_name).as_object();
        const auto& state = side.at("transport_state").as_object();
        if (side.if_contains("tile_plan_provenance_identity") != nullptr ||
            state.if_contains("provenance_identity") != nullptr)
          throw std::runtime_error(
              "paired checkpoint recursively embedded owner provenance");
        const auto& recipe_side = aggregate.at(side_name).as_object();
        const auto& rows = recipe_side.at("rows").as_array();
        if (rows.size() != 1)
          throw std::runtime_error(
              "paired checkpoint lost its compact row order");
        const auto& row = rows.front().as_object();
        if (row.size() != 3 || row.if_contains("tile") == nullptr ||
            row.if_contains("exact_identity") == nullptr ||
            row.if_contains("entries") == nullptr ||
            row.if_contains("prepared_row") != nullptr)
          throw std::runtime_error(
              "paired checkpoint retained a deep prepared-row copy");
        const auto& entry_facts = row.at("entries").as_array();
        if (entry_facts.size() != 1 ||
            entry_facts.front().as_object().size() != 4)
          throw std::runtime_error(
              "paired checkpoint lost its compact exact entry facts");
      }
    }

    auto tampered_payload = payload;
    auto& tampered_diagnostics = tampered_payload
        .at("retained_line_results").as_array().front().as_object()
        .at("result").as_object().at("diagnostics").as_object();
    if (tampered_diagnostics.at("divergent_cancellation_mode") !=
            "exact-singleton" ||
        tampered_diagnostics.at("divergent_relative_tolerance") != "" ||
        tampered_diagnostics.at("divergent_cancellation_provenance") != "")
      throw std::runtime_error(
          "paired checkpoint did not expose an exact-singleton aggregate");
    tampered_diagnostics["divergent_relative_tolerance"] = "1e-20";
    diffexp2::checkpoint::write_atomic(
        checkpoint_tampered, container.header_json,
        json::serialize(tampered_payload));
    const auto tampered_restore = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint_tampered},
        {"expected_identity", "pair-roundtrip"}});
    if (tampered_restore.at("status") == "ok") {
      (void)request(json::object{
          {"schema", 2}, {"op", "session.close"},
          {"session", tampered_restore.at("session")}});
      throw std::runtime_error(
          "exact-singleton aggregate accepted a stray bounded tolerance");
    }
    if (std::string(tampered_restore.at("detail").as_string()).find(
            "aggregate divergent-cancellation policy is incomplete") ==
        std::string::npos)
      throw std::runtime_error(
          "malformed aggregate cancellation failed for the wrong reason: " +
          json::serialize(tampered_restore));

    auto count_tampered_payload = payload;
    count_tampered_payload.at("retained_line_results").as_array()
        .front().as_object().at("result").as_object()
        .at("diagnostics").as_object()
        ["bounded_cancelled_divergent_coefficients"] = 1;
    diffexp2::checkpoint::write_atomic(
        checkpoint_count_tampered, container.header_json,
        json::serialize(count_tampered_payload));
    const auto count_tampered_restore = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint_count_tampered},
        {"expected_identity", "pair-roundtrip"}});
    if (count_tampered_restore.at("status") == "ok") {
      (void)request(json::object{
          {"schema", 2}, {"op", "session.close"},
          {"session", count_tampered_restore.at("session")}});
      throw std::runtime_error(
          "exact-singleton aggregate accepted a bounded-cancelled count");
    }
    if (std::string(count_tampered_restore.at("detail").as_string()).find(
            "exact-singleton aggregate cannot report bounded-cancelled") ==
        std::string::npos)
      throw std::runtime_error(
          "aggregate bounded-count tamper failed for the wrong reason: " +
          json::serialize(count_tampered_restore));

    auto legacy_payload = payload;
    auto& legacy_diagnostics = legacy_payload
        .at("retained_line_results").as_array().front().as_object()
        .at("result").as_object().at("diagnostics").as_object();
    for (const auto* key : {
             "bounded_cancelled_divergent_coefficients",
             "divergent_cancellation_mode",
             "divergent_relative_tolerance",
             "divergent_cancellation_provenance"})
      if (legacy_diagnostics.erase(key) != 1)
        throw std::runtime_error(
            "paired checkpoint omitted a current cancellation diagnostic");
    diffexp2::checkpoint::write_atomic(
        checkpoint_legacy, container.header_json,
        json::serialize(legacy_payload));
    const auto legacy_restore = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint_legacy},
        {"expected_identity", "pair-roundtrip"}});
    if (legacy_restore.at("status") != "ok")
      throw std::runtime_error(
          "legacy exact-singleton aggregate failed to restore: " +
          json::serialize(legacy_restore));
    (void)request(json::object{
        {"schema", 2}, {"op", "session.close"},
        {"session", legacy_restore.at("session")}});

    const auto restored = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint}, {"expected_identity", "pair-roundtrip"}});
    restored_session = std::string(restored.at("session").as_string());
    const auto restored_stats = session_stats(restored_session);
    const auto restored_export = export_line(restored_session, lines.front());
    if (restored.at("status") != "ok" ||
        restored_stats.at("transport_states") != 0 ||
        restored_stats.at("tile_plans") != 0 ||
        restored_stats.at("locals") != 0 ||
        restored_stats.at("line_results") != 7 ||
        restored_export.at("value") != first_export_value)
      throw std::runtime_error(
          "first paired hidden-owner restore changed visibility or value");
    const auto saved_second = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", restored_session}, {"path", checkpoint_second},
        {"checkpoint_identity", "pair-roundtrip-second"}});
    if (saved_second.at("status") != "ok")
      throw std::runtime_error("second checkpoint save failed");
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", restored_session}});
    restored_session.clear();
    const auto restored_second = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint_second},
        {"expected_identity", "pair-roundtrip-second"}});
    restored_session =
        std::string(restored_second.at("session").as_string());
    const auto second_stats = session_stats(restored_session);
    const auto second_export = export_line(restored_session, lines.front());
    if (restored_second.at("status") != "ok" ||
        second_stats.at("transport_states") != 0 ||
        second_stats.at("tile_plans") != 0 ||
        second_export.at("value") != first_export_value)
      throw std::runtime_error(
          "second paired hidden-owner restore changed visibility or value");

    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", restored_session}});
    restored_session.clear();
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    session.clear();
    std::remove(checkpoint.c_str());
    std::remove(checkpoint_second.c_str());
    std::remove(checkpoint_tampered.c_str());
    std::remove(checkpoint_count_tampered.c_str());
    std::remove(checkpoint_legacy.c_str());
    std::cout << "PASS: retained transport-pair contraction\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    if (!restored_session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored_session}});
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    std::remove(checkpoint.c_str());
    std::remove(checkpoint_second.c_str());
    std::remove(checkpoint_tampered.c_str());
    std::remove(checkpoint_count_tampered.c_str());
    std::remove(checkpoint_legacy.c_str());
    std::cerr << "FAIL: retained transport-pair contraction: "
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
