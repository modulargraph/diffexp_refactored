#include "diffexp/kernel/checkpoint.hpp"
#include "diffexp/kernel/json_codec.hpp"

#include <boost/json.hpp>

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp::kernel::run_recurrence_json(json::serialize(value)))
      .as_object();
}

json::array prescriptions(std::int32_t sign) {
  return json::array{json::object{
      {"factor_exact", "x-2/3"}, {"sign", sign},
      {"multiplicity", 1}, {"leading_coefficient_sign", 1}}};
}

std::string prepare_fractional_chart(const std::string& session,
                                     const std::string& prefix,
                                     std::int32_t sign) {
  auto value = json::parse(R"json({
    "schema":2,"op":"chart.prepare","session":"placeholder",
    "key":"placeholder","identity":"placeholder",
    "analytic":{
      "geometry":{"center_exact":"2/3","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,
        "prescriptions":[]},
      "principal_matrix":[[{"exact":"((1/2)+eps)/t",
        "proven_zero":false}]],
      "native_scc_capabilities":{"regular":false,"identity_gauge":true,
        "identity_v":true,"no_pseudo":true}},
    "scc":{"components":[[0]],"structural_edges":[[0,0]],
      "condensation_edges":[],"topological_order":[0],
      "coupling_depth":0},
    "problem":{"domain":"acb","precision_bits":256,
      "d":1,"fb":0,"w":3,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[
          {"s":0,"e":[[0,0,"1/2"]]},
          {"s":1,"e":[[0,0,"1"]]}],
        "rat":[],"val":[0]}],
      "d0_inverse":"1","blocks":[[0]],
      "assembly":{"identity":true,"poly":[],"rat":[],"val":[0]},
      "chop_digits":0}
  })json").as_object();
  value["session"] = session;
  value["key"] = prefix + "-fractional-chart";
  value["identity"] = prefix + "-fractional-chart-identity";
  value.at("analytic").as_object().at("geometry").as_object()
      ["prescriptions"] = prescriptions(sign);
  const auto prepared = request(std::move(value));
  if (prepared.at("status") != "ok")
    throw std::runtime_error(
        "fractional chart.prepare: " + json::serialize(prepared));
  return std::string(prepared.at("chart").as_string());
}

std::string prepare_regular_chart(const std::string& session,
                                  const std::string& prefix,
                                  const std::string& role,
                                  const std::string& center,
                                  std::int32_t sign) {
  auto value = json::parse(R"json({
    "schema":2,"op":"chart.prepare","session":"placeholder",
    "key":"placeholder","identity":"placeholder",
    "analytic":{
      "geometry":{"center_exact":"0","scale_exact":"1",
        "radius_exact":"2","infinite_radius":false,
        "prescriptions":[]},
      "principal_matrix":[[{"exact":"0","proven_zero":true}]],
      "native_scc_capabilities":{"regular":true,"identity_gauge":true,
        "identity_v":true,"no_pseudo":true}},
    "scc":{"components":[[0]],"structural_edges":[],
      "condensation_edges":[],"topological_order":[0],
      "coupling_depth":0},
    "problem":{"domain":"acb","precision_bits":256,
      "d":1,"fb":0,"w":3,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],"val":[null]}],
      "d0_inverse":"1","blocks":[[0]],
      "assembly":{"identity":true,"poly":[],"rat":[],"val":[0]},
      "chop_digits":0}
  })json").as_object();
  value["session"] = session;
  value["key"] = prefix + "-" + role + "-chart";
  value["identity"] = prefix + "-" + role + "-chart-identity";
  auto& geometry = value.at("analytic").as_object()
                       .at("geometry").as_object();
  geometry["center_exact"] = center;
  geometry["prescriptions"] = prescriptions(sign);
  const auto prepared = request(std::move(value));
  if (prepared.at("status") != "ok")
    throw std::runtime_error(
        "regular chart.prepare: " + json::serialize(prepared));
  return std::string(prepared.at("chart").as_string());
}

json::object fixed_run(const std::string& a, const std::string& b) {
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= 4; ++n) {
    shifts.emplace_back(
        n == 0 ? a : (a == "0" ? std::to_string(n)
                                : std::to_string(2 * n + 1) + "/2"));
    schedule.push_back(json::array{json::object{
        {"case", n == 0 ? "R" : "T"}, {"da", std::to_string(n)},
        {"db", "0"}}});
  }
  return json::object{
      {"nmax", 4}, {"p", 0}, {"has_initial", true},
      {"adaptive_probe", false}, {"a_target", a}, {"b_target", b},
      {"a_shift_min", 0}, {"a_shifts", std::move(shifts)},
      {"schedule", std::move(schedule)},
      {"initial", json::array{"1", "0", "0"}},
      {"initial_validity", json::array{2}}, {"source", nullptr},
      {"return_u", false}};
}

std::string solve_local(const std::string& session,
                        const std::string& chart,
                        const std::string& prefix,
                        const std::string& role,
                        const std::string& center,
                        const std::string& a,
                        const std::string& b,
                        std::int32_t sign) {
  const auto checkpoint = prefix + "-" + role + "-local";
  const auto solved = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart}, {"run", fixed_run(a, b)},
      {"metadata", json::object{
           {"chart", json::object{{"center_exact", center},
                                   {"scale_exact", "1"},
                                   {"radius", "2"},
                                   {"infinite_radius", false}}},
           {"tag", json::object{
                {"a", json::object{{"domain", "rational"},
                                     {"canonical", a}}},
                {"b", json::object{{"domain", "rational"},
                                     {"canonical", b}}}}},
           {"prescriptions", prescriptions(sign)},
           {"checkpoint_identity", checkpoint}}}});
  if (solved.at("status") != "ok")
    throw std::runtime_error("local.solve: " + json::serialize(solved));
  return std::string(solved.at("local").as_string());
}

json::object topology(std::int32_t sign) {
  return json::object{
      {"singular_points", json::array{"2/3"}},
      {"boundary_points", json::array{}},
      {"complex_projections", json::array{}},
      {"branch_sheets", json::array{json::object{
           {"factor_exact", "x-2/3"}, {"sign", sign}}}}};
}

json::object arm(const std::string& endpoint,
                 const std::string& producing,
                 const std::string& receiving,
                 std::int32_t sign) {
  return json::object{{"from_exact", "2/3"}, {"to_exact", endpoint},
                      {"charts", json::array{producing, receiving}},
                      {"topology", topology(sign)}};
}

json::object integrand_row(const std::string& identity) {
  const json::array one{"1", "0", "0", "0", "0"};
  const json::array zero{"0", "0", "0", "0", "0"};
  return json::object{
      {"schema", "diffexp3-prepared-rational-local-row-v1"},
      {"columns", 1}, {"exact_identity", identity},
      {"entries", json::array{json::object{
           {"column", 0},
           {"multiplier", json::object{
                {"epsilon_shift", 0}, {"center_pole_order", 0},
                {"kernels", json::array{one, zero, zero}},
                {"exact_identity", identity},
                {"proven_zero", false}}}}}}};
}

json::object arm_execution(const std::string& basis,
                           const std::string& prefix,
                           const std::string& name) {
  json::array basis_set;
  basis_set.emplace_back(basis);
  json::array receiving_basis;
  receiving_basis.push_back(std::move(basis_set));
  return json::object{
      {"receiving_basis", std::move(receiving_basis)},
      {"integrand_rows", json::array{
           integrand_row(prefix + ":" + name + ":row:0"),
           integrand_row(prefix + ":" + name + ":row:1")}}};
}

json::object run_arms_request(const std::string& session,
                              const std::string& plan,
                              const std::string& anchor,
                              const std::string& lower_basis,
                              const std::string& upper_basis,
                              const std::string& prefix) {
  return json::object{
      {"schema", 2}, {"op", "integration.run_arms"},
      {"session", session}, {"tile_plan", plan}, {"anchor", anchor},
      {"tile_plan_checkpoint_identity", prefix + "-plan"},
      {"anchor_checkpoint_identity", prefix + "-anchor-local"},
      {"epsilon", json::object{{"min", 0}, {"max", 2},
                                {"required_complete_max", 1},
                                {"match_required_complete_max", 1}}},
      {"refinement", json::object{
           {"relative_tolerance", "1e-40"}, {"max_steps", 2}}},
      {"checkpoint_policy", json::object{
           {"schema", "diffexp3-deterministic-arm-checkpoints-v1"},
           {"root", prefix + "-march"}}},
      {"lower", arm_execution(lower_basis, prefix, "lower")},
      {"upper", arm_execution(upper_basis, prefix, "upper")}};
}

std::complex<double> coefficient(const json::object& evaluated,
                                 std::size_t epsilon_index) {
  const auto& raw = evaluated.at("value").as_object()
                        .at("coefficients").as_array().at(epsilon_index)
                        .as_array();
  return {std::stod(std::string(raw.at(0).as_string())),
          std::stod(std::string(raw.at(1).as_string()))};
}

json::object evaluate_local(const std::string& session,
                            const std::string& local) {
  return request(json::object{
      {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
      {"local", local}, {"point", json::object{{"exact", "0"}}},
      {"options", json::object{{"tail_estimate", false}}},
      {"output_digits", 50}});
}

bool match_provenance_ok(const std::string& path, std::int32_t sign) {
  const auto payload = json::parse(
      diffexp::kernel::checkpoint::read(path).payload_json).as_object();
  for (const auto& raw_hop :
       payload.at("retained_planned_match_hops").as_array()) {
    const auto& hop = raw_hop.as_object();
    const auto& handoff = hop.at("handoff").as_object();
    if (handoff.at("arm") != "lower" || handoff.at("match") != 0)
      continue;
    const auto& native = hop.at("native_match").as_object();
    const auto& basis = native.at("basis_sources").as_array().front()
                            .as_object();
    const auto& incoming = native.at("incoming_source").as_object();
    const auto witness = json::parse(std::string(
        native.at("exact_lattice_canonical_witness").as_string()));
    return native.at("basis_point_exact") == "1/3" &&
        native.at("incoming_point_exact") == "-1/3" &&
        native.at("physical_match_point_exact") == "1/3" &&
        basis.at("requested_imaginary_sign") == sign &&
        basis.at("effective_imaginary_sign").is_null() &&
        incoming.at("requested_imaginary_sign") == sign &&
        incoming.at("effective_imaginary_sign") == sign &&
        witness.as_object().at("schema") ==
            "diffexp3-native-acb-unit-leading-saturation-proof-v1";
  }
  return false;
}

struct PhaseCase {
  std::complex<double> epsilon0;
  std::complex<double> epsilon1;
  std::complex<double> restored_epsilon0;
  std::complex<double> restored_epsilon1;
  bool provenance = false;
  bool caller_proof_free = false;
};

PhaseCase run_phase_case(std::int32_t sign, const std::string& label) {
  const auto prefix = "phase-" + label;
  const auto checkpoint_path =
      (std::filesystem::temp_directory_path() /
       (prefix + "-planned-match.de2cp")).string();
  std::remove(checkpoint_path.c_str());
  std::string session;
  std::string restored_session;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"}, {"domain", "acb"},
        {"precision_bits", 256}, {"output_digits", 50},
        {"chart_capacity", 4}, {"local_capacity", 20},
        {"match_capacity", 4}, {"tile_plan_capacity", 2},
        {"line_result_capacity", 8}});
    session = std::string(created.at("session").as_string());
    const auto producing = prepare_fractional_chart(
        session, prefix, sign);
    const auto receiving = prepare_regular_chart(
        session, prefix, "receiving", "0", sign);
    const auto upper = prepare_regular_chart(
        session, prefix, "upper", "4/3", sign);
    const auto anchor = solve_local(
        session, producing, prefix, "anchor", "2/3", "1/2", "1",
        sign);
    const auto lower_basis = solve_local(
        session, receiving, prefix, "lower-basis", "0", "0", "0",
        sign);
    const auto upper_basis = solve_local(
        session, upper, prefix, "upper-basis", "4/3", "0", "0",
        sign);
    const auto planned = request(json::object{
        {"schema", 2}, {"op", "tile.plan"}, {"session", session},
        {"checkpoint_identity", prefix + "-plan"},
        {"division_order", 3},
        {"lower", arm("0", producing, receiving, sign)},
        {"upper", arm("4/3", producing, upper, sign)}});
    if (planned.at("status") != "ok")
      throw std::runtime_error("tile.plan: " + json::serialize(planned));
    auto march_request = run_arms_request(
        session, std::string(planned.at("tile_plan").as_string()), anchor,
        lower_basis, upper_basis, prefix);
    const auto serialized_request = json::serialize(march_request);
    const bool caller_proof_free =
        serialized_request.find("exact_lattice") == std::string::npos &&
        serialized_request.find("native_unit_saturation") ==
            std::string::npos &&
        serialized_request.find("basis_imaginary_sign") ==
            std::string::npos &&
        serialized_request.find("incoming_imaginary_sign") ==
            std::string::npos &&
        serialized_request.find("receiving_rim") == std::string::npos;
    const auto marched = request(std::move(march_request));
    if (marched.at("status") != "ok")
      throw std::runtime_error(
          "integration.run_arms: " + json::serialize(marched));
    const auto& final_local = marched.at("arms").as_object()
                                  .at("lower").as_object()
                                  .at("final_local").as_object();
    const auto final_handle = std::string(
        final_local.at("local").as_string());
    const auto evaluated = evaluate_local(session, final_handle);
    if (evaluated.at("status") != "ok")
      throw std::runtime_error(
          "local.evaluate: " + json::serialize(evaluated));
    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
        {"path", checkpoint_path},
        {"checkpoint_identity", prefix + "-checkpoint"}});
    if (saved.at("status") != "ok")
      throw std::runtime_error(
          "checkpoint.save: " + json::serialize(saved));
    const bool provenance = match_provenance_ok(checkpoint_path, sign);
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    session.clear();
    const auto restored = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint_path},
        {"expected_identity", prefix + "-checkpoint"}});
    if (restored.at("status") != "ok")
      throw std::runtime_error(
          "checkpoint.restore: " + json::serialize(restored));
    restored_session = std::string(restored.at("session").as_string());
    const auto restored_evaluated = evaluate_local(
        restored_session, final_handle);
    if (restored_evaluated.at("status") != "ok")
      throw std::runtime_error(
          "restored local.evaluate: " +
          json::serialize(restored_evaluated));
    PhaseCase result{
        coefficient(evaluated, 0), coefficient(evaluated, 1),
        coefficient(restored_evaluated, 0),
        coefficient(restored_evaluated, 1), provenance,
        caller_proof_free};
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", restored_session}});
    restored_session.clear();
    std::remove(checkpoint_path.c_str());
    return result;
  } catch (...) {
    if (!session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", session}});
    if (!restored_session.empty())
      (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                                 {"session", restored_session}});
    std::remove(checkpoint_path.c_str());
    throw;
  }
}

bool close(std::complex<double> left, std::complex<double> right,
           double tolerance = 1e-12) {
  return std::abs(left - right) < tolerance;
}

}  // namespace

int main() {
  try {
    const auto plus = run_phase_case(1, "plus");
    const auto minus = run_phase_case(-1, "minus");
    const std::complex<double> i(0.0, 1.0);
    const auto plus_leading = i / std::sqrt(3.0);
    const auto minus_leading = -i / std::sqrt(3.0);
    const auto plus_ratio = -std::log(3.0) + i * std::acos(-1.0);
    const auto minus_ratio = -std::log(3.0) - i * std::acos(-1.0);
    const bool ok = plus.provenance && minus.provenance &&
        plus.caller_proof_free && minus.caller_proof_free &&
        close(plus.epsilon0, plus_leading) &&
        close(minus.epsilon0, minus_leading) &&
        close(plus.epsilon1 / plus.epsilon0, plus_ratio) &&
        close(minus.epsilon1 / minus.epsilon0, minus_ratio) &&
        close(minus.epsilon0, std::conj(plus.epsilon0)) &&
        close(minus.epsilon1, std::conj(plus.epsilon1)) &&
        close(plus.restored_epsilon0, plus.epsilon0) &&
        close(plus.restored_epsilon1, plus.epsilon1) &&
        close(minus.restored_epsilon0, minus.epsilon0) &&
        close(minus.restored_epsilon1, minus.epsilon1);
    if (!ok)
      std::cerr << "plus eps0=" << plus.epsilon0
                << " eps1=" << plus.epsilon1
                << " ratio=" << plus.epsilon1 / plus.epsilon0
                << " restored=" << plus.restored_epsilon0 << ','
                << plus.restored_epsilon1
                << " provenance=" << plus.provenance
                << " caller-proof-free=" << plus.caller_proof_free << '\n'
                << "minus eps0=" << minus.epsilon0
                << " eps1=" << minus.epsilon1
                << " ratio=" << minus.epsilon1 / minus.epsilon0
                << " restored=" << minus.restored_epsilon0 << ','
                << minus.restored_epsilon1
                << " provenance=" << minus.provenance
                << " caller-proof-free=" << minus.caller_proof_free << '\n';
    std::cout << (ok ? "PASS" : "FAIL")
              << ": phase-sensitive planned production match\n";
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << "FAIL: phase-sensitive planned production match: "
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
