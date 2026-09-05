#include "diffexp/kernel/checkpoint.hpp"
#include "diffexp/kernel/json_codec.hpp"

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
  return json::parse(diffexp::kernel::run_recurrence_json(json::serialize(value)))
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

json::object topology() {
  return json::object{
      {"singular_points", json::array{}},
      {"boundary_points", json::array{}},
      {"complex_projections", json::array{}},
      {"branch_sheets", json::array{
           json::object{{"factor_exact", "t"}, {"sign", -1}}}}};
}

std::string prepare_chart(const std::string& session,
                          const std::string& key,
                          const std::string& center,
                          const std::string& radius) {
  auto value = json::parse(R"json({
    "schema":2,"op":"chart.prepare","session":"placeholder",
    "key":"endpoint-batch-chart","identity":"endpoint-batch-operator",
    "analytic":{
      "geometry":{"center_exact":"1","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,
        "prescriptions":[{"factor_exact":"t","sign":-1,
          "multiplicity":1,"leading_coefficient_sign":1}]},
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
  })json").as_object();
  value["session"] = session;
  value["key"] = key;
  auto& geometry = value.at("analytic").as_object()
      .at("geometry").as_object();
  geometry["center_exact"] = center;
  geometry["scale_exact"] = "2";
  geometry["radius_exact"] = radius;
  const auto prepared = request(std::move(value));
  if (prepared.at("status") != "ok")
    throw std::runtime_error("chart.prepare: " + json::serialize(prepared));
  return std::string(prepared.at("chart").as_string());
}

json::object solve_local(const std::string& session,
                         const std::string& chart,
                         const std::string& center,
                         const std::string& radius,
                         const std::string& checkpoint,
                         json::array initial) {
  json::array schedule_row{
      json::object{{"case", "R"}, {"da", "0"}, {"db", "0"}},
      json::object{{"case", "R"}, {"da", "0"}, {"db", "0"}}};
  json::array schedule;
  schedule.push_back(std::move(schedule_row));
  const auto solved = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart},
      {"run", json::object{
           {"nmax", 0}, {"p", 0}, {"has_initial", true},
           {"adaptive_probe", false}, {"a_target", "0"},
           {"b_target", "0"}, {"a_shift_min", 0},
           {"a_shifts", json::array{"0"}},
           {"schedule", std::move(schedule)},
           {"initial", std::move(initial)},
           {"initial_validity", json::array{2, 2}},
           {"source", nullptr}, {"return_u", false}}},
      {"metadata", json::object{
           {"chart", json::object{
                {"center_exact", center}, {"scale_exact", "2"},
                {"radius", radius}, {"infinite_radius", false}}},
           {"tag", json::object{
                {"a", json::object{{"domain", "rational"},
                                     {"canonical", "0"}}},
                {"b", json::object{{"domain", "rational"},
                                     {"canonical", "0"}}}}},
           {"prescriptions", prescriptions()},
           {"checkpoint_identity", checkpoint}}}});
  if (solved.at("status") != "ok")
    throw std::runtime_error("local.solve: " + json::serialize(solved));
  return solved;
}

json::object plan_arm(const std::string& endpoint,
                      json::array charts) {
  return json::object{
      {"from_exact", "0"}, {"to_exact", endpoint},
      {"charts", std::move(charts)},
      {"topology", topology()}};
}

json::object prepare_plan(const std::string& session,
                          const std::string& checkpoint,
                          const std::string& endpoint,
                          json::array charts) {
  return request(json::object{
      {"schema", 2}, {"op", "tile.plan_arm"}, {"session", session},
      {"checkpoint_identity", checkpoint}, {"division_order", 3},
      {"arm", plan_arm(endpoint, std::move(charts))}});
}

json::object run_state(const std::string& session,
                       const json::object& plan,
                       const std::string& plan_checkpoint,
                       const json::object& anchor,
                       const std::vector<json::object>& basis,
                       const std::string& checkpoint_root) {
  json::array receiving_columns;
  for (const auto& column : basis)
    receiving_columns.push_back(column.at("local"));
  json::array receiving_basis;
  if (!receiving_columns.empty())
    receiving_basis.push_back(std::move(receiving_columns));
  return request(json::object{
      {"schema", 2}, {"op", "transport.run_arm"},
      {"session", session}, {"tile_plan", plan.at("tile_plan")},
      {"tile_plan_checkpoint_identity", plan_checkpoint},
      {"anchor", anchor.at("local")},
      {"anchor_checkpoint_identity", anchor.at("checkpoint_identity")},
      {"arm", "upper"},
      {"receiving_basis", std::move(receiving_basis)},
      {"epsilon", json::object{
           {"min", 0}, {"max", 2}, {"required_complete_max", 2},
           {"match_required_complete_max", 2}}},
      {"refinement", json::object{
           {"relative_tolerance", "1e-30"}, {"max_steps", 0}}},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp3-deterministic-arm-checkpoints-v1"},
           {"root", checkpoint_root}}}});
}

json::object multiplier(const std::string& scale,
                        std::int32_t epsilon_shift,
                        const std::string& identity) {
  return json::object{
      {"epsilon_shift", epsilon_shift}, {"center_pole_order", 0},
      {"kernels", json::array{
           json::array{scale}, json::array{"0"}, json::array{"0"}}},
      {"exact_identity", identity}, {"proven_zero", false}};
}

json::object cancellation_row(const std::string& identity,
                              std::int32_t shift = 0) {
  return json::object{
      {"schema", "diffexp3-prepared-rational-local-row-v1"},
      {"columns", 2}, {"exact_identity", identity},
      {"entries", json::array{
           json::object{{"column", 0},
                        {"multiplier", multiplier("1", shift,
                                                  identity + ":c0")}},
           json::object{{"column", 1},
                        {"multiplier", multiplier("1", shift,
                                                  identity + ":c1")}}}}};
}

json::object first_component_row(const std::string& identity,
                                 std::int32_t shift = 0,
                                 bool malformed = false) {
  return json::object{
      {"schema", "diffexp3-prepared-rational-local-row-v1"},
      {"columns", 2}, {"exact_identity", identity},
      {"entries", json::array{json::object{
           {"column", 0},
           {"multiplier", multiplier(
                malformed ? "not-a-rational" : "1", shift,
                identity + ":c0")}}}}};
}

json::object public_prefix_first_component_row(
    const std::string& identity, std::int32_t shift) {
  return json::object{
      {"schema", "diffexp3-prepared-rational-local-row-v1"},
      {"columns", 2}, {"exact_identity", identity},
      {"entries", json::array{json::object{
           {"column", 0},
           {"multiplier", json::object{
                {"epsilon_shift", shift}, {"center_pole_order", 0},
                {"analytic_coefficients", json::array{
                     json::object{
                         {"numerator", json::array{"1"}},
                         {"denominator", json::array{"1"}}},
                     json::object{
                         {"numerator", json::array{"0"}},
                         {"denominator", json::array{"1"}}}}},
                {"exact_identity", identity + ":c0"},
                {"proven_zero", false}}}}}}};
}

json::object observable(const std::string& identity,
                        const std::string& checkpoint,
                        json::object row,
                        std::int32_t epsilon_min = 0,
                        std::int32_t epsilon_max = 1,
                        std::int32_t required_max = 1) {
  return json::object{
      {"identity", identity}, {"checkpoint_identity", checkpoint},
      {"integrand_row", std::move(row)},
      {"publication_relative_tolerance", "1e-8"},
      {"epsilon", json::object{
           {"min", epsilon_min}, {"max", epsilon_max},
           {"required_complete_max", required_max}}}};
}

json::object endpoint_batch(const std::string& session,
                            const json::object& state,
                            json::array observables,
                            const std::string& checkpoint_root,
                            std::uint32_t threads = 1) {
  return request(json::object{
      {"schema", 2}, {"op", "transport.endpoint_batch"},
      {"session", session},
      {"transport_state", state.at("transport_state")},
      {"transport_state_checkpoint_identity",
       state.at("checkpoint_identity")},
      {"transport_state_provenance_identity",
       state.at("provenance_identity")},
      {"checkpoint_policy", json::object{
           {"schema",
            "diffexp3-deterministic-transport-endpoint-checkpoints-v1"},
           {"root", checkpoint_root}}},
      {"observables", std::move(observables)}, {"threads", threads}});
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

json::object export_endpoint(const std::string& session,
                             const json::object& endpoint) {
  return request(json::object{
      {"schema", 2}, {"op", "endpoint.export"}, {"session", session},
      {"endpoint", endpoint.at("endpoint")},
      {"checkpoint_identity", endpoint.at("checkpoint_identity")},
      {"output_digits", 50}});
}

double midpoint(const json::object& exported, std::size_t index) {
  const auto& coefficient = exported.at("value").as_object()
      .at("coefficients").as_array().at(index).as_array();
  return std::stod(std::string(coefficient.front().as_string()));
}

bool payload_has_handle(const json::array& records,
                        const std::string& handle) {
  for (const auto& raw : records)
    if (raw.as_object().at("handle").as_string() == handle) return true;
  return false;
}

}  // namespace

int main() {
  const std::string checkpoint =
      "/tmp/diffexp3-transport-endpoint-batch.de2cp";
  const std::string checkpoint_second =
      "/tmp/diffexp3-transport-endpoint-batch-second.de2cp";
  std::remove(checkpoint.c_str());
  std::remove(checkpoint_second.c_str());
  std::string session;
  std::string restored_session;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "rational"}, {"output_digits", 50},
        {"chart_capacity", 2}, {"local_capacity", 4},
        {"match_capacity", 2}, {"tile_plan_capacity", 2},
        {"transport_state_capacity", 2}, {"endpoint_capacity", 5}});
    session = std::string(created.at("session").as_string());
    const auto anchor_chart = prepare_chart(
        session, "endpoint-batch-anchor-chart", "0", "2");
    const auto endpoint_chart = prepare_chart(
        session, "endpoint-batch-final-chart", "1", "2");
    const auto anchor = solve_local(
        session, anchor_chart, "0", "2", "endpoint-batch-anchor",
        json::array{"1", "0", "0", "-1", "0", "0"});
    const auto basis_first = solve_local(
        session, endpoint_chart, "1", "2", "endpoint-batch-basis-1",
        json::array{"1", "0", "0", "0", "0", "0"});
    const auto basis_second = solve_local(
        session, endpoint_chart, "1", "2", "endpoint-batch-basis-2",
        json::array{"0", "0", "0", "1", "0", "0"});
    const std::vector<json::object> receiving_basis{
        basis_first, basis_second};
    const auto centered_plan = prepare_plan(
        session, "endpoint-centered-plan", "1",
        json::array{anchor_chart, endpoint_chart});
    const auto regular_plan = prepare_plan(
        session, "endpoint-regular-plan", "1/2",
        json::array{anchor_chart});
    if (centered_plan.at("status") != "ok" ||
        regular_plan.at("status") != "ok" ||
        centered_plan.at("tile_plan") == regular_plan.at("tile_plan"))
      throw std::runtime_error(
          "independent endpoint plans failed: " +
          json::serialize(centered_plan) + " / " +
          json::serialize(regular_plan));
    const auto centered_state = run_state(
        session, centered_plan, "endpoint-centered-plan", anchor,
        receiving_basis, "endpoint-centered-state");
    const auto regular_state = run_state(
        session, regular_plan, "endpoint-regular-plan", anchor,
        {}, "endpoint-regular-state");
    if (centered_state.at("status") != "ok" ||
        regular_state.at("status") != "ok")
      throw std::runtime_error(
          "endpoint transport states failed: " +
          json::serialize(centered_state) + " / " +
          json::serialize(regular_state));

    const auto before_zero = session_stats(session);
    const auto zero = endpoint_batch(
        session, centered_state, json::array{}, "endpoint-zero");
    const auto after_zero = session_stats(session);
    const auto centered_after_zero = state_stats(session, centered_state);
    if (zero.at("status") != "ok" || zero.at("observables") != 0 ||
        !zero.at("endpoints").as_array().empty() ||
        before_zero.at("endpoints") != after_zero.at("endpoints") ||
        before_zero.at("locals") != after_zero.at("locals") ||
        before_zero.at("matches") != after_zero.at("matches") ||
        counter(after_zero, "transport_endpoint_batches") !=
            counter(before_zero, "transport_endpoint_batches") + 1 ||
        counter(after_zero, "transport_endpoint_rows") !=
            counter(before_zero, "transport_endpoint_rows") ||
        counter(centered_after_zero, "endpoint_batch_operations") != 1 ||
        counter(centered_after_zero, "endpoint_rows") != 0)
      throw std::runtime_error(
          "zero-row endpoint batch changed handles or counters incorrectly");

    const auto one = endpoint_batch(
        session, centered_state,
        json::array{observable(
            "center-cancel", "center-cancel-checkpoint",
            cancellation_row("center-cancel-row"))},
        "endpoint-one");
    if (one.at("status") != "ok" || one.at("centered") != true ||
        one.at("local_endpoint_exact") != "0" ||
        one.at("endpoints").as_array().size() != 1 ||
        one.at("no_projected_local_publication") != true)
      throw std::runtime_error(
          "centered cancellation endpoint failed: " + json::serialize(one));
    const auto one_export = export_endpoint(
        session, one.at("endpoints").as_array().front().as_object());
    if (one_export.at("status") != "ok" ||
        std::abs(midpoint(one_export, 0)) > 1e-45)
      throw std::runtime_error(
          "centered row cancellation was not applied before endpoint limiting");

    const auto before_malformed = session_stats(session);
    const auto centered_before_malformed = state_stats(
        session, centered_state);
    const auto malformed = endpoint_batch(
        session, centered_state,
        json::array{
            observable("valid-prefix", "valid-prefix-checkpoint",
                       cancellation_row("valid-prefix-row")),
            observable("malformed", "malformed-checkpoint",
                       first_component_row("malformed-row", 0, true))},
        "endpoint-malformed");
    const auto after_malformed = session_stats(session);
    const auto centered_after_malformed = state_stats(
        session, centered_state);
    if (malformed.at("status") != "error" ||
        before_malformed.at("endpoints") != after_malformed.at("endpoints") ||
        before_malformed.at("endpoint_limits") !=
            after_malformed.at("endpoint_limits") ||
        before_malformed.at("transport_endpoint_batches") !=
            after_malformed.at("transport_endpoint_batches") ||
        centered_before_malformed.at("endpoint_batch_operations") !=
            centered_after_malformed.at("endpoint_batch_operations") ||
        after_malformed.at("pending_endpoint_limits") != 0)
      throw std::runtime_error(
          "malformed later row did not roll back atomically");

    const auto before_bad_tolerance = session_stats(session);
    auto bad_tolerance_observable = observable(
        "bad-publication-tolerance",
        "bad-publication-tolerance-checkpoint",
        cancellation_row("bad-publication-tolerance-row"));
    bad_tolerance_observable["publication_relative_tolerance"] = "0";
    const auto bad_tolerance = endpoint_batch(
        session, centered_state,
        json::array{std::move(bad_tolerance_observable)},
        "endpoint-bad-publication-tolerance");
    const auto after_bad_tolerance = session_stats(session);
    if (bad_tolerance.at("status") != "error" ||
        before_bad_tolerance.at("endpoints") !=
            after_bad_tolerance.at("endpoints") ||
        before_bad_tolerance.at("transport_endpoint_batches") !=
            after_bad_tolerance.at("transport_endpoint_batches") ||
        after_bad_tolerance.at("pending_endpoint_limits") != 0)
      throw std::runtime_error(
          "invalid endpoint publication tolerance changed retained state");

    json::array many_observables;
    for (int index = 0; index < 3; ++index)
      many_observables.push_back(observable(
          "center-many-" + std::to_string(index),
          "center-many-checkpoint-" + std::to_string(index),
          cancellation_row("center-many-row-" + std::to_string(index))));
    const auto many = endpoint_batch(
        session, centered_state, std::move(many_observables),
        "endpoint-many", 3);
    std::set<std::string> endpoint_handles;
    if (many.at("status") != "ok" ||
        many.at("endpoints").as_array().size() != 3 ||
        many.at("requested_observable_threads") != 3 ||
        many.at("observable_worker_threads") != 3)
      throw std::runtime_error(
          "three-row endpoint batch failed: " + json::serialize(many));
    for (std::size_t index = 0; index < 3; ++index) {
      const auto& endpoint = many.at("endpoints").as_array()[index].as_object();
      if (counter(endpoint, "request_index") != index ||
          std::string(endpoint.at("observable_identity").as_string()) !=
              "center-many-" + std::to_string(index) ||
          !endpoint_handles.insert(
              std::string(endpoint.at("endpoint").as_string())).second)
        throw std::runtime_error(
            "three-row endpoint ordering or identity changed");
    }

    const auto regular = endpoint_batch(
        session, regular_state,
        json::array{observable(
            "regular-shift", "regular-shift-checkpoint",
            public_prefix_first_component_row("regular-shift-row", 1),
            0, 2, 2)},
        "endpoint-regular");
    if (regular.at("status") != "ok" || regular.at("centered") != false ||
        regular.at("local_endpoint_exact") != "1/4" ||
        regular.at("derived_rim") != -1 ||
        regular.at("endpoints").as_array().size() != 1)
      throw std::runtime_error(
          "noncentered regular-point endpoint failed: " +
          json::serialize(regular));
    const auto regular_export = export_endpoint(
        session, regular.at("endpoints").as_array().front().as_object());
    if (regular_export.at("status") != "ok" ||
        regular_export.at("value").as_object().at("min") != 0 ||
        regular_export.at("value").as_object().at("max") != 2 ||
        std::abs(midpoint(regular_export, 0)) > 1e-45 ||
        std::abs(midpoint(regular_export, 1) - 1.0) > 1e-40)
      throw std::runtime_error(
          "noncentered public-prefix evaluation or shifted epsilon window is wrong: " +
          json::serialize(regular_export));

    const auto before_capacity = session_stats(session);
    const auto capacity = endpoint_batch(
        session, centered_state,
        json::array{observable(
            "over-capacity", "over-capacity-checkpoint",
            cancellation_row("over-capacity-row"))},
        "endpoint-capacity");
    const auto after_capacity = session_stats(session);
    if (capacity.at("status") != "error" ||
        before_capacity.at("endpoints") != after_capacity.at("endpoints") ||
        before_capacity.at("transport_endpoint_batches") !=
            after_capacity.at("transport_endpoint_batches") ||
        after_capacity.at("pending_endpoint_limits") != 0)
      throw std::runtime_error(
          "endpoint capacity rejection changed retained state");

    const auto final_stats = session_stats(session);
    const auto centered_stats = state_stats(session, centered_state);
    const auto regular_stats = state_stats(session, regular_state);
    if (counter(final_stats, "transport_endpoint_batches") != 4 ||
        counter(final_stats, "transport_endpoint_rows") != 5 ||
        counter(final_stats, "endpoint_limits") != 5 ||
        counter(centered_stats, "endpoint_batch_operations") != 3 ||
        counter(centered_stats, "endpoint_rows") != 4 ||
        counter(regular_stats, "endpoint_batch_operations") != 1 ||
        counter(regular_stats, "endpoint_rows") != 1 ||
        final_stats.at("locals") != 4 || final_stats.at("matches") != 0 ||
        final_stats.at("line_results") != 0)
      throw std::runtime_error(
          "endpoint batch counters or scratch allocation are dishonest");

    std::vector<json::object> retained_endpoints;
    retained_endpoints.push_back(
        one.at("endpoints").as_array().front().as_object());
    for (const auto& raw : many.at("endpoints").as_array())
      retained_endpoints.push_back(raw.as_object());
    retained_endpoints.push_back(
        regular.at("endpoints").as_array().front().as_object());
    const auto regular_export_value = regular_export.at("value");
    for (const auto& state : {centered_state, regular_state})
      if (request(json::object{
              {"schema", 2}, {"op", "transport.release"},
              {"session", session},
              {"transport_state", state.at("transport_state")}})
              .at("status") != "ok")
        throw std::runtime_error("transport state release failed");
    for (const auto& plan : {centered_plan, regular_plan})
      if (request(json::object{
              {"schema", 2}, {"op", "tile.release"},
              {"session", session}, {"tile_plan", plan.at("tile_plan")}})
              .at("status") != "ok")
        throw std::runtime_error("endpoint plan release failed");
    std::vector<json::object> public_locals{
        anchor, basis_first, basis_second};
    const auto centered_final =
        centered_state.at("final_local").as_object();
    if (centered_final.at("local") != anchor.at("local"))
      public_locals.push_back(centered_final);
    for (const auto& local : public_locals)
      if (request(json::object{
              {"schema", 2}, {"op", "local.release"},
              {"session", session}, {"local", local.at("local")}})
              .at("status") != "ok")
        throw std::runtime_error("endpoint fixture local release failed");
    for (const auto& endpoint : retained_endpoints)
      if (export_endpoint(session, endpoint).at("status") != "ok")
        throw std::runtime_error(
            "compact endpoint did not survive state/plan/local release");

    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
        {"path", checkpoint},
        {"checkpoint_identity", "endpoint-batch-roundtrip"}});
    if (saved.at("status") != "ok")
      throw std::runtime_error("checkpoint.save: " + json::serialize(saved));
    const auto container = diffexp::kernel::checkpoint::read(checkpoint);
    const auto payload = json::parse(container.payload_json).as_object();
    if (payload.at("schema") != 9 ||
        payload.at("retained_transport_states").as_array().size() != 2 ||
        !payload_has_handle(payload.at("retained_transport_states").as_array(),
                            std::string(centered_state.at("transport_state").as_string())) ||
        !payload_has_handle(payload.at("retained_transport_states").as_array(),
                            std::string(regular_state.at("transport_state").as_string())) ||
        container.payload_json.find("private:") != std::string::npos)
      throw std::runtime_error(
          "compact endpoint checkpoint lost its state closure or retained scratch locals");
    for (const auto& raw : payload.at("retained_endpoints").as_array()) {
      const auto& record = raw.as_object();
      if (record.at("schema") !=
              "diffexp3-retained-transport-endpoint-result-v2" ||
          record.at("source").as_object().if_contains("row") == nullptr ||
          record.if_contains("prepared_row") == nullptr ||
          record.at("source").as_object().at("row").as_object()
                  .if_contains("prepared_row") != nullptr)
        throw std::runtime_error(
            "endpoint checkpoint is not compact and row-bound");
    }

    const auto restored = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint},
        {"expected_identity", "endpoint-batch-roundtrip"}});
    if (restored.at("status") != "ok")
      throw std::runtime_error(
          "checkpoint.restore: " + json::serialize(restored));
    restored_session = std::string(restored.at("session").as_string());
    const auto restored_stats = session_stats(restored_session);
    const auto restored_export = export_endpoint(
        restored_session, retained_endpoints.back());
    if (restored_stats.at("transport_states") != 0 ||
        restored_stats.at("tile_plans") != 0 ||
        restored_stats.at("locals") != 0 ||
        restored_stats.at("endpoints") != 5 ||
        restored_export.at("value") != regular_export_value)
      throw std::runtime_error(
          "first hidden transport-endpoint restore changed visibility or value");
    const auto saved_second = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", restored_session}, {"path", checkpoint_second},
        {"checkpoint_identity", "endpoint-batch-roundtrip-second"}});
    if (saved_second.at("status") != "ok")
      throw std::runtime_error("second endpoint checkpoint save failed");
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", restored_session}});
    restored_session.clear();
    const auto restored_second = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint_second},
        {"expected_identity", "endpoint-batch-roundtrip-second"}});
    if (restored_second.at("status") != "ok")
      throw std::runtime_error(
          "second endpoint checkpoint restore failed: " +
          json::serialize(restored_second));
    restored_session =
        std::string(restored_second.at("session").as_string());
    const auto second_stats = session_stats(restored_session);
    const auto second_export = export_endpoint(
        restored_session, retained_endpoints.back());
    if (second_stats.at("transport_states") != 0 ||
        second_stats.at("tile_plans") != 0 ||
        second_export.at("value") != regular_export_value)
      throw std::runtime_error(
          "second hidden transport-endpoint restore changed visibility or value");

    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", restored_session}});
    restored_session.clear();
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    session.clear();
    std::remove(checkpoint.c_str());
    std::remove(checkpoint_second.c_str());
    std::cout << "PASS: retained transport endpoint batch\n";
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
    std::cerr << "FAIL: retained transport endpoint batch: "
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
