#include "diffexp2/checkpoint.hpp"
#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <initializer_list>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
      .as_object();
}

bool is_error(const json::object& response) {
  return response.if_contains("status") != nullptr &&
         response.at("status") == "error";
}

void require_ok(const json::object& response, const char* label) {
  if (response.if_contains("status") == nullptr ||
      response.at("status") != "ok")
    throw std::runtime_error(std::string(label) + ": " +
                             json::serialize(response));
}

json::value nested(json::array value) {
  return json::value(std::move(value));
}

json::object epsilon_rational(std::int32_t valuation,
                              std::initializer_list<const char*> numerator,
                              std::initializer_list<const char*> denominator) {
  json::array p;
  json::array q;
  for (const auto* value : numerator) p.emplace_back(value);
  for (const auto* value : denominator) q.emplace_back(value);
  return json::object{{"zero", false},
                      {"valuation", valuation},
                      {"numerator", std::move(p)},
                      {"denominator", std::move(q)}};
}

json::object physical_entry(std::uint32_t row, std::uint32_t column,
                            json::object value) {
  return json::object{{"r", row}, {"c", column}, {"v", std::move(value)}};
}

json::object physical_ode(const std::string& owner,
                          const std::string& payload,
                          json::array q, json::array c) {
  return json::object{
      {"schema", "diffexp2-physical-cleared-ode-v1"},
      {"basis", "physical-original-master"},
      {"theta_coordinate", "local-t"},
      {"owner_signature_identity", owner},
      {"payload_identity", payload},
      {"q", std::move(q)},
      {"c", std::move(c)}};
}

json::object scalar_shift(std::int32_t shift, const std::string& value) {
  return json::object{{"s", shift}, {"v", value}};
}

json::object matrix_shift(
    std::int32_t shift,
    std::initializer_list<json::array> entries) {
  json::array encoded;
  for (const auto& entry : entries) encoded.push_back(entry);
  return json::object{{"s", shift}, {"e", std::move(encoded)}};
}

json::array null_valuations(std::uint32_t dimension) {
  json::array values;
  for (std::uint32_t i = 0; i < dimension * dimension; ++i)
    values.push_back(nullptr);
  return values;
}

json::object scc_certificate(std::uint32_t dimension) {
  if (dimension == 1)
    return json::object{
        {"components", json::array{nested(json::array{0})}},
        {"structural_edges", json::array{}},
        {"condensation_edges", json::array{}},
        {"topological_order", json::array{0}},
        {"coupling_depth", 0}};
  return json::object{
      {"components", json::array{nested(json::array{0, 1})}},
      {"structural_edges", json::array{nested(json::array{0, 1}),
                                        nested(json::array{1, 0})}},
      {"condensation_edges", json::array{}},
      {"topological_order", json::array{0}},
      {"coupling_depth", 0}};
}

json::object regular_problem(const std::string& domain,
                             const std::string& owner,
                             const std::string& payload,
                             bool include_physical = true) {
  const auto multiplier = epsilon_rational(-1, {"1", "-1"},
                                            {"1", "1"});
  json::object problem{
      {"domain", domain},
      {"d", 1}, {"fb", 0}, {"w", 5},
      {"d_lags", json::array{nested(json::array{scalar_shift(0, "1")})}},
      {"denominators", json::array{}},
      {"nhat_lags", json::array{
          json::object{{"poly", json::array{}},
                       {"rat", json::array{}},
                       {"val", json::array{nullptr}}},
          json::object{{"poly", json::array{matrix_shift(
                                0, {json::array{0, 0, "1"}})}},
                       {"rat", json::array{}},
                       {"val", json::array{0}}}}},
      {"d0_inverse", "1"},
      {"blocks", json::array{nested(json::array{0})}},
      {"assembly", json::object{{"identity", true},
                                  {"poly", json::array{}},
                                  {"rat", json::array{}},
                                  {"val", json::array{0}}}},
      {"chop_digits", 0}};
  if (domain == "acb") problem["precision_bits"] = 256;
  if (include_physical) {
    json::array c{
        nested(json::array{}),
        nested(json::array{physical_entry(0, 0, multiplier)})};
    problem["physical_ode"] = physical_ode(
        owner, payload, json::array{multiplier}, std::move(c));
  }
  return problem;
}

json::object cancellation_problem(const std::string& owner,
                                  const std::string& payload) {
  auto problem = regular_problem("rational", owner, payload, false);
  const auto one = epsilon_rational(0, {"1"}, {"1"});
  const auto minus_two = epsilon_rational(0, {"-2"}, {"1"});
  problem["physical_ode"] = physical_ode(
      owner, payload, json::array{one, minus_two},
      json::array{nested(json::array{}),
                  nested(json::array{physical_entry(0, 0, one)}),
                  nested(json::array{physical_entry(0, 0, minus_two)})});
  return problem;
}

json::object spectral_problem(const std::string& owner,
                              const std::string& payload) {
  const auto one = epsilon_rational(0, {"1"}, {"1"});
  const auto two = epsilon_rational(0, {"2"}, {"1"});
  const auto eps = epsilon_rational(1, {"1"}, {"1"});
  return json::object{
      {"domain", "acb"}, {"precision_bits", 256},
      {"d", 2}, {"fb", 0}, {"w", 5},
      {"d_lags", json::array{nested(json::array{scalar_shift(0, "1")})}},
      {"denominators", json::array{}},
      {"nhat_lags", json::array{
          json::object{{"poly", json::array{}},
                       {"rat", json::array{}},
                       {"val", null_valuations(2)}},
          json::object{{"poly", json::array{matrix_shift(
              0, {json::array{0, 0, "1"},
                  json::array{1, 1, "2"}})}},
                       {"rat", json::array{}},
                       {"val", json::array{0, nullptr, nullptr, 0}}}}},
      {"d0_inverse", "1"},
      {"blocks", json::array{nested(json::array{0}),
                               nested(json::array{1})}},
      {"assembly", json::object{
          {"identity", false},
          {"poly", json::array{
              matrix_shift(0, {json::array{0, 0, "1"},
                               json::array{1, 1, "1"}}),
              matrix_shift(1, {json::array{0, 1, "1"}})}},
          {"rat", json::array{}},
          {"val", json::array{0, 1, nullptr, 0}}}},
      {"physical_ode", physical_ode(
          owner, payload, json::array{one},
          json::array{nested(json::array{}),
                      nested(json::array{physical_entry(0, 0, one),
                                         physical_entry(0, 1, eps),
                                         physical_entry(1, 1, two)})})},
      {"chop_digits", 50}};
}

json::object logarithmic_problem(const std::string& owner,
                                 const std::string& payload) {
  const auto one = epsilon_rational(0, {"1"}, {"1"});
  const auto diagonal = epsilon_rational(0, {"1/2", "1/3"}, {"1"});
  return json::object{
      {"domain", "rational"},
      {"d", 2}, {"fb", 0}, {"w", 4},
      {"d_lags", json::array{nested(json::array{scalar_shift(0, "1")})}},
      {"denominators", json::array{}},
      {"nhat_lags", json::array{json::object{
          {"poly", json::array{
              matrix_shift(0, {json::array{0, 0, "1/2"},
                               json::array{0, 1, "1"},
                               json::array{1, 1, "1/2"}}),
              matrix_shift(1, {json::array{0, 0, "1/3"},
                               json::array{1, 1, "1/3"}})}},
          {"rat", json::array{}},
          {"val", json::array{0, 0, nullptr, 0}}}}},
      {"d0_inverse", "1"},
      {"blocks", json::array{nested(json::array{0, 1})}},
      {"assembly", json::object{{"identity", true},
                                  {"poly", json::array{}},
                                  {"rat", json::array{}},
                                  {"val", json::array{0, nullptr,
                                                       nullptr, 0}}}},
      {"physical_ode", physical_ode(
          owner, payload, json::array{one},
          json::array{nested(json::array{
              physical_entry(0, 0, diagonal),
              physical_entry(0, 1, one),
              physical_entry(1, 1, diagonal)})})},
      {"chop_digits", 0}};
}

json::object prepare_chart(const std::string& session,
                           const std::string& key,
                           const std::string& identity,
                           json::object problem) {
  return request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", key}, {"identity", identity},
      {"analytic", json::object{{"residual_fixture", key}}},
      {"scc", scc_certificate(
          static_cast<std::uint32_t>(problem.at("d").as_int64()))},
      {"problem", std::move(problem)}});
}

json::object regular_run(std::uint32_t dimension, std::uint32_t block_count,
                         std::uint32_t width, std::uint32_t nmax,
                         std::uint32_t seed_component,
                         bool explicit_empty_source = false) {
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    shifts.emplace_back(std::to_string(n));
    json::array row;
    for (std::uint32_t block = 0; block < block_count; ++block)
      row.push_back(json::object{{"case", n == 0 ? "R" : "T"},
                                 {"da", std::to_string(n)},
                                 {"db", "0"}});
    schedule.push_back(std::move(row));
  }
  json::array initial;
  json::array validity;
  for (std::uint32_t component = 0; component < dimension; ++component) {
    for (std::uint32_t epsilon = 0; epsilon < width; ++epsilon)
      initial.emplace_back(component == seed_component && epsilon == 0
                               ? "1" : "0");
    validity.push_back(static_cast<std::int64_t>(width - 1));
  }
  json::value source = nullptr;
  if (explicit_empty_source) {
    const auto points = static_cast<std::size_t>(nmax + 1);
    json::array frames;
    json::array source_validity;
    json::array present;
    for (std::size_t point = 0; point < points; ++point) {
      present.push_back(false);
      for (std::uint32_t component = 0; component < dimension; ++component) {
        source_validity.push_back(nullptr);
        for (std::uint32_t epsilon = 0; epsilon < width; ++epsilon)
          frames.emplace_back("0");
      }
    }
    source = json::object{{"frames", std::move(frames)},
                          {"validity", std::move(source_validity)},
                          {"present", std::move(present)}};
  }
  return json::object{
      {"nmax", nmax}, {"p", 0}, {"has_initial", true},
      {"adaptive_probe", false}, {"a_target", "0"},
      {"b_target", "0"}, {"a_shift_min", 0},
      {"a_shifts", std::move(shifts)}, {"schedule", std::move(schedule)},
      {"initial", std::move(initial)},
      {"initial_validity", std::move(validity)},
      {"source", std::move(source)}, {"return_u", false}};
}

json::object logarithmic_run() {
  // The native log convention is (eps log(t))^p/p!.  These boundary rows
  // encode f_2=eps t^(1/2+eps/3) and
  // f_1=eps log(t) t^(1/2+eps/3), which solve the retained Jordan equation.
  return json::object{
      {"nmax", 0}, {"p", 1}, {"has_initial", true},
      {"adaptive_probe", false}, {"a_target", "1/2"},
      {"b_target", "1/3"}, {"a_shift_min", 0},
      {"a_shifts", json::array{"1/2"}},
      {"schedule", json::array{nested(json::array{json::object{
          {"case", "R"}, {"da", "0"}, {"db", "0"}}})}},
      {"initial", json::array{
          // log 0: component 0, component 1
          "0", "0", "0", "0",  "0", "1", "0", "0",
          // log 1: component 0, component 1
          "1", "0", "0", "0",  "0", "0", "0", "0"}},
      {"initial_validity", json::array{3, 3, 3, 3}},
      {"source", nullptr}, {"return_u", false}};
}

json::array prescriptions(std::int32_t sign) {
  if (sign == 0) return {};
  return json::array{json::object{
      {"factor_exact", "t"}, {"sign", sign}, {"multiplicity", 1},
      {"leading_coefficient_sign", 1}}};
}

json::object local_metadata(const std::string& checkpoint,
                            const std::string& a = "0",
                            const std::string& b = "0",
                            std::int32_t prescription_sign = 0) {
  return json::object{
      {"chart", json::object{{"center_exact", "0"},
                              {"scale_exact", "1"},
                              {"radius", "2"},
                              {"infinite_radius", false}}},
      {"tag", json::object{
          {"a", json::object{{"domain", "rational"},
                              {"canonical", a}}},
          {"b", json::object{{"domain", "rational"},
                              {"canonical", b}}}}},
      {"prescriptions", prescriptions(prescription_sign)},
      {"checkpoint_identity", checkpoint}};
}

json::object solve_local(const std::string& session, const std::string& chart,
                         json::object run, json::object metadata) {
  return request(json::object{{"schema", 2}, {"op", "local.solve"},
                              {"session", session}, {"chart", chart},
                              {"run", std::move(run)},
                              {"metadata", std::move(metadata)}});
}

json::object residual_request(const std::string& session,
                              const std::string& local,
                              const json::object& binding,
                              const std::string& point = "1/10",
                              std::optional<std::int32_t> rim = std::nullopt) {
  json::object options{{"tail_estimate", false}};
  if (rim.has_value()) options["imaginary_sign"] = *rim;
  return json::object{
      {"schema", 2}, {"op", "local.certify_residual"},
      {"session", session}, {"local", local},
      {"point", json::object{{"exact", point}}},
      {"options", std::move(options)},
      {"relative_tolerance", "1e-25"},
      {"scope", "stored_truncation"}, {"include_residual", false},
      {"operator_identity", binding.at("operator_identity")},
      {"source_identity", binding.at("source_identity")},
      {"checkpoint_identity", binding.at("local_checkpoint_identity")},
      {"analytic_metadata", binding.at("analytic_metadata")},
      {"owner_signature_identity", binding.at("owner_signature_identity")},
      {"physical_payload_identity", binding.at("physical_payload_identity")},
      {"provenance_identity", binding.at("provenance_identity")}};
}

json::object checkpoint_payload(const std::filesystem::path& path) {
  return json::parse(diffexp2::checkpoint::read(path.string()).payload_json)
      .as_object();
}

const json::object& record_by_handle(const json::object& payload,
                                     const char* collection,
                                     const std::string& key,
                                     const std::string& handle) {
  for (const auto& raw : payload.at(collection).as_array()) {
    const auto& record = raw.as_object();
    if (std::string(record.at(key).as_string()) == handle) return record;
  }
  throw std::runtime_error(std::string("checkpoint lost ") + collection +
                           " handle " + handle);
}

}  // namespace

int main() {
  const auto stem = std::filesystem::temp_directory_path() /
      ("diffexp2_physical_owner_residual_" +
       std::to_string(static_cast<long long>(::getpid())));
  const auto saved_path = stem.string() + "_saved.de2cp";
  const auto resaved_path = stem.string() + "_resaved.de2cp";

  const auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", "rational"},
      {"precision_bits", 256}, {"output_digits", 50},
      {"chart_capacity", 8}, {"local_capacity", 10}});
  require_ok(created, "rational session.create");
  const auto session = std::string(created.at("session").as_string());

  const std::string regular_owner = "de2-operator-regular-rational-v2";
  const auto regular_prepared = prepare_chart(
      session, "physical-regular", regular_owner,
      regular_problem("rational", regular_owner,
                      "de2-physical-ode-regular-rational-v2"));
  require_ok(regular_prepared, "regular chart.prepare");
  const auto regular_chart =
      std::string(regular_prepared.at("chart").as_string());
  const auto regular_solved = solve_local(
      session, regular_chart, regular_run(1, 1, 5, 24, 0),
      local_metadata("physical-regular-local-v2"));
  require_ok(regular_solved, "regular local.solve");
  const auto regular_local =
      std::string(regular_solved.at("local").as_string());
  const auto regular_binding = regular_solved.at("residual_binding")
                                   .as_object().at("binding").as_object();

  // An explicitly sourced primitive remains unsupported even if every source
  // row is absent: source ownership is a provenance fact, not a numeric-zero
  // guess made from this finite run.
  const auto sourced_solved = solve_local(
      session, regular_chart, regular_run(1, 1, 5, 3, 0, true),
      local_metadata("physical-sourced-local-v1"));
  require_ok(sourced_solved, "sourced local.solve");

  // A legacy hand-authored chart is still solvable, but cannot silently gain
  // a residual equation it never retained.
  const auto legacy_prepared = prepare_chart(
      session, "legacy-no-physical", "legacy-no-physical-identity",
      regular_problem("rational", "unused", "unused", false));
  require_ok(legacy_prepared, "legacy chart.prepare");
  const auto legacy_solved = solve_local(
      session, std::string(legacy_prepared.at("chart").as_string()),
      regular_run(1, 1, 5, 4, 0),
      local_metadata("legacy-no-physical-local-v1"));
  require_ok(legacy_solved, "legacy local.solve");

  // Q(0)=1 is part of the unique exact payload contract.
  auto bad_q_problem = regular_problem(
      "rational", "de2-operator-bad-q",
      "de2-physical-ode-bad-q");
  bad_q_problem.at("physical_ode").as_object().at("q").as_array()
      .front().as_object().at("denominator").as_array().front() = "2";
  const auto bad_q = prepare_chart(
      session, "bad-q", "de2-operator-bad-q", std::move(bad_q_problem));

  // The payload owner token and native chart exact identity are one binding.
  const auto wrong_owner = prepare_chart(
      session, "wrong-owner", "de2-operator-different-owner",
      regular_problem("rational", "de2-operator-payload-owner",
                      "de2-physical-ode-wrong-owner"));

  // q(1/2,eps)=0 makes the cleared equation vacuous.  It must reject rather
  // than certify 0=0.
  const std::string cancellation_owner = "de2-operator-q-cancellation";
  const auto cancellation_prepared = prepare_chart(
      session, "q-cancellation", cancellation_owner,
      cancellation_problem(cancellation_owner,
                           "de2-physical-ode-q-cancellation"));
  require_ok(cancellation_prepared, "q-cancellation chart.prepare");
  const auto cancellation_solved = solve_local(
      session, std::string(cancellation_prepared.at("chart").as_string()),
      regular_run(1, 1, 5, 12, 0),
      local_metadata("q-cancellation-local-v1"));
  require_ok(cancellation_solved, "q-cancellation local.solve");
  const auto cancellation_binding = cancellation_solved.at("residual_binding")
      .as_object().at("binding").as_object();
  const auto cancellation_rejected = request(residual_request(
      session, std::string(cancellation_solved.at("local").as_string()),
      cancellation_binding, "1/2"));

  // Fractional power + regulator exponent + a genuine p=1 log sector on the
  // negative axis, bound to the retained -i0 prescription.
  const std::string log_owner = "de2-operator-fractional-log-minus-rim";
  const auto log_prepared = prepare_chart(
      session, "fractional-log", log_owner,
      logarithmic_problem(log_owner,
                          "de2-physical-ode-fractional-log-minus-rim"));
  require_ok(log_prepared, "fractional/log chart.prepare");
  const auto log_solved = solve_local(
      session, std::string(log_prepared.at("chart").as_string()),
      logarithmic_run(),
      local_metadata("fractional-log-local-v1", "1/2", "1/3", -1));
  require_ok(log_solved, "fractional/log local.solve");
  const auto log_binding = log_solved.at("residual_binding")
      .as_object().at("binding").as_object();
  const auto log_minus = request(residual_request(
      session, std::string(log_solved.at("local").as_string()),
      log_binding, "-1/2", -1));
  auto wrong_rim_request = residual_request(
      session, std::string(log_solved.at("local").as_string()),
      log_binding, "-1/2", 1);
  const auto wrong_rim = request(std::move(wrong_rim_request));

  const auto released = request(json::object{
      {"schema", 2}, {"op", "chart.release"}, {"session", session},
      {"chart", regular_chart}});
  require_ok(released, "regular chart.release");
  const auto valid = request(residual_request(
      session, regular_local, regular_binding));
  require_ok(valid, "physical residual after chart release");
  auto full_scope_request = residual_request(
      session, regular_local, regular_binding);
  full_scope_request["scope"] = "full_local_solution";
  const auto full_scope = request(std::move(full_scope_request));
  require_ok(full_scope, "full-local residual");

  auto operator_payload_tamper = residual_request(
      session, regular_local, regular_binding);
  operator_payload_tamper["theta_operator"] = json::object{};
  const auto operator_payload_rejected = request(operator_payload_tamper);
  auto source_payload_tamper = residual_request(
      session, regular_local, regular_binding);
  source_payload_tamper["source"] = json::object{};
  const auto source_payload_rejected = request(source_payload_tamper);
  auto operator_identity_tamper = residual_request(
      session, regular_local, regular_binding);
  operator_identity_tamper["operator_identity"] = "attacker-operator";
  const auto operator_identity_rejected = request(operator_identity_tamper);
  auto source_identity_tamper = residual_request(
      session, regular_local, regular_binding);
  source_identity_tamper["source_identity"] = "attacker-source";
  const auto source_identity_rejected = request(source_identity_tamper);
  auto owner_identity_tamper = residual_request(
      session, regular_local, regular_binding);
  owner_identity_tamper["owner_signature_identity"] = "attacker-owner";
  const auto owner_identity_rejected = request(owner_identity_tamper);
  auto physical_identity_tamper = residual_request(
      session, regular_local, regular_binding);
  physical_identity_tamper["physical_payload_identity"] = "attacker-qc";
  const auto physical_identity_rejected = request(physical_identity_tamper);
  auto provenance_tamper = residual_request(
      session, regular_local, regular_binding);
  provenance_tamper["provenance_identity"] = "attacker-provenance";
  const auto provenance_rejected = request(provenance_tamper);
  auto analytic_tamper = residual_request(
      session, regular_local, regular_binding);
  analytic_tamper.at("analytic_metadata").as_object()
      .at("chart").as_object()["center_exact"] = "1";
  const auto analytic_rejected = request(analytic_tamper);
  auto checkpoint_tamper = residual_request(
      session, regular_local, regular_binding);
  checkpoint_tamper["checkpoint_identity"] = "attacker-checkpoint";
  const auto checkpoint_rejected = request(checkpoint_tamper);
  auto tail_tamper = residual_request(session, regular_local, regular_binding);
  tail_tamper.at("options").as_object()["tail_estimate"] = true;
  const auto tail_rejected = request(tail_tamper);

  const auto saved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
      {"path", saved_path},
      {"checkpoint_identity", "physical-owner-session-v1"}});
  require_ok(saved, "checkpoint.save");
  const auto original_payload = checkpoint_payload(saved_path);
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", session}});

  const auto restored = request(json::object{
      {"schema", 2}, {"op", "checkpoint.restore"}, {"path", saved_path},
      {"expected_identity", "physical-owner-session-v1"}});
  require_ok(restored, "checkpoint.restore");
  const auto restored_session =
      std::string(restored.at("session").as_string());
  const auto restored_valid = request(residual_request(
      restored_session, regular_local, regular_binding));
  require_ok(restored_valid, "physical residual after hidden-owner restore");
  const auto resaved = request(json::object{
      {"schema", 2}, {"op", "checkpoint.save"},
      {"session", restored_session}, {"path", resaved_path},
      {"checkpoint_identity", "physical-owner-resaved-v1"}});
  require_ok(resaved, "checkpoint resave");
  const auto restored_payload = checkpoint_payload(resaved_path);

  const auto& original_local = record_by_handle(
      original_payload, "retained_locals", "handle", regular_local);
  const auto& restored_local = record_by_handle(
      restored_payload, "retained_locals", "handle", regular_local);
  const auto& original_chart = record_by_handle(
      original_payload, "prepared_charts", "handle", regular_chart);
  const auto& restored_chart = record_by_handle(
      restored_payload, "prepared_charts", "handle", regular_chart);
  const auto& original_physical = original_chart.at("request").as_object()
      .at("problem").as_object().at("physical_ode");
  const auto& restored_physical = restored_chart.at("request").as_object()
      .at("problem").as_object().at("physical_ode");

  // Acb specialization, with V=[[1,eps],[0,1]].  The physical C has an
  // off-diagonal eps term and is not Nhat; using the spectral matrix directly
  // would fail this residual.
  const auto acb_created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", "acb"},
      {"precision_bits", 256}, {"output_digits", 50},
      {"chart_capacity", 2}, {"local_capacity", 2}});
  require_ok(acb_created, "Acb session.create");
  const auto acb_session =
      std::string(acb_created.at("session").as_string());
  auto acb_bad_q_problem = spectral_problem(
      "de2-operator-acb-bad-q", "de2-physical-ode-acb-bad-q");
  acb_bad_q_problem.at("physical_ode").as_object().at("q").as_array()
      .front().as_object().at("denominator").as_array().front() =
          "[1 +/- 1e-40]";
  const auto acb_bad_q = prepare_chart(
      acb_session, "acb-bad-q", "de2-operator-acb-bad-q",
      std::move(acb_bad_q_problem));
  const std::string spectral_owner = "de2-operator-acb-spectral-v";
  const auto spectral_prepared = prepare_chart(
      acb_session, "acb-spectral-v", spectral_owner,
      spectral_problem(spectral_owner,
                       "de2-physical-ode-acb-spectral-v"));
  require_ok(spectral_prepared, "Acb spectral chart.prepare");
  const auto spectral_solved = solve_local(
      acb_session, std::string(spectral_prepared.at("chart").as_string()),
      regular_run(2, 2, 5, 24, 1),
      local_metadata("acb-spectral-local-v1"));
  require_ok(spectral_solved, "Acb spectral local.solve");
  const auto spectral_binding = spectral_solved.at("residual_binding")
      .as_object().at("binding").as_object();
  auto spectral_residual_request = residual_request(
      acb_session, std::string(spectral_solved.at("local").as_string()),
      spectral_binding);
  spectral_residual_request["include_residual"] = true;
  const auto spectral_valid = request(std::move(spectral_residual_request));
  require_ok(spectral_valid, "Acb spectral physical residual");

  const bool ok =
      regular_solved.at("residual_binding").as_object().at("status") ==
          "available" &&
      valid.at("verdict") == "pass" && valid.at("epsilon_min") == -1 &&
      valid.at("epsilon_max") == 3 &&
      valid.at("basis") == "physical-original-master" &&
      valid.at("owner_signature_identity") ==
          regular_binding.at("owner_signature_identity") &&
      valid.at("physical_payload_identity") ==
          regular_binding.at("physical_payload_identity") &&
      full_scope.at("verdict") == "inconclusive" &&
      sourced_solved.at("residual_binding").as_object().at("status") ==
          "unsupported" &&
      legacy_solved.at("residual_binding").as_object().at("status") ==
          "unsupported" &&
      is_error(bad_q) && is_error(wrong_owner) &&
      is_error(cancellation_rejected) &&
      log_solved.at("residual_binding").as_object().at("status") ==
          "available" &&
      log_minus.at("status") == "ok" && log_minus.at("verdict") == "pass" &&
      log_minus.at("imaginary_sign") == -1 && is_error(wrong_rim) &&
      is_error(operator_payload_rejected) &&
      is_error(source_payload_rejected) &&
      is_error(operator_identity_rejected) &&
      is_error(source_identity_rejected) &&
      is_error(owner_identity_rejected) &&
      is_error(physical_identity_rejected) && is_error(provenance_rejected) &&
      is_error(analytic_rejected) && is_error(checkpoint_rejected) &&
      is_error(tail_rejected) &&
      restored_valid.at("verdict") == "pass" &&
      original_payload.at("schema") == 9 &&
      original_local.at("schema") == "diffexp2-retained-local-v4" &&
      !original_local.at("equation_owner_restore").is_null() &&
      original_local.at("equation_owner_restore") ==
          restored_local.at("equation_owner_restore") &&
      original_local.at("residual_operator_restore") ==
          restored_local.at("residual_operator_restore") &&
      original_physical == restored_physical &&
      is_error(acb_bad_q) &&
      spectral_solved.at("residual_binding").as_object().at("status") ==
          "available" &&
      spectral_valid.at("verdict") == "pass" &&
      spectral_valid.at("arithmetic_enclosed") == true;

  if (!ok) {
    std::cerr << "regular solve: " << json::serialize(regular_solved) << '\n'
              << "valid: " << json::serialize(valid) << '\n'
              << "full scope: " << json::serialize(full_scope) << '\n'
              << "sourced: " << json::serialize(sourced_solved) << '\n'
              << "legacy: " << json::serialize(legacy_solved) << '\n'
              << "bad Q: " << json::serialize(bad_q) << '\n'
              << "wrong owner: " << json::serialize(wrong_owner) << '\n'
              << "q cancellation: " << json::serialize(cancellation_rejected)
              << '\n' << "log -i0: " << json::serialize(log_minus) << '\n'
              << "wrong rim: " << json::serialize(wrong_rim) << '\n'
              << "restored valid: " << json::serialize(restored_valid) << '\n'
              << "Acb bad Q: " << json::serialize(acb_bad_q) << '\n'
              << "spectral solve: " << json::serialize(spectral_solved) << '\n'
              << "spectral residual: " << json::serialize(spectral_valid)
              << '\n';
  }

  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", restored_session}});
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                             {"session", acb_session}});
  std::error_code error;
  std::filesystem::remove(saved_path, error);
  std::filesystem::remove(resaved_path, error);
  std::cout << (ok ? "PASS" : "FAIL")
            << ": physical equation-owned persistent residual\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
