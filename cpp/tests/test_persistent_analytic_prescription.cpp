#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

#include <unistd.h>

namespace json = boost::json;

namespace {

json::object request(json::object value) {
  return json::parse(diffexp2::run_recurrence_json(json::serialize(value)))
      .as_object();
}

json::array prescriptions(std::optional<std::int32_t> sign,
                          std::uint32_t multiplicity = 1) {
  if (!sign.has_value()) return {};
  return json::array{json::object{
      {"factor_exact", "x"}, {"sign", *sign},
      {"multiplicity", multiplicity},
      {"leading_coefficient_sign", 1}}};
}

std::string prepare_chart(const std::string& session,
                          const std::string& identity,
                          const std::string& scale,
                          std::optional<std::int32_t> sign,
                          std::uint32_t multiplicity = 1) {
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
    "problem":{"domain":"rational","d":1,"fb":0,"w":3,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],"val":[null]}],
      "d0_inverse":"1","blocks":[[0]],
      "assembly":{"identity":true,"poly":[],"rat":[],"val":[0]},
      "chop_digits":0}
  })json").as_object();
  value["session"] = session;
  value["key"] = identity;
  value["identity"] = identity;
  auto& geometry = value.at("analytic").as_object()
      .at("geometry").as_object();
  geometry["scale_exact"] = scale;
  geometry["prescriptions"] = prescriptions(sign, multiplicity);
  const auto response = request(std::move(value));
  if (response.at("status") != "ok")
    throw std::runtime_error("chart.prepare: " + json::serialize(response));
  return std::string(response.at("chart").as_string());
}

std::string solve_local(const std::string& session,
                        const std::string& chart,
                        const std::string& checkpoint,
                        const std::string& scale,
                        const std::string& a,
                        std::optional<std::int32_t> sign) {
  json::array schedule_row;
  schedule_row.push_back(json::object{{"case", "R"}, {"da", "0"},
                                      {"db", "0"}});
  json::array schedule;
  schedule.push_back(std::move(schedule_row));
  const auto response = request(json::object{
      {"schema", 2}, {"op", "local.solve"}, {"session", session},
      {"chart", chart},
      {"run", json::object{
          {"nmax", 0}, {"p", 0}, {"has_initial", true},
          {"adaptive_probe", false}, {"a_target", a},
          {"b_target", "0"}, {"a_shift_min", 0},
          {"a_shifts", json::array{a}},
          {"schedule", std::move(schedule)},
          {"initial", json::array{"1", "0", "0"}},
          {"initial_validity", json::array{2}},
          {"source", nullptr}, {"return_u", false}}},
      {"metadata", json::object{
          {"chart", json::object{{"center_exact", "0"},
                                  {"scale_exact", scale},
                                  {"radius", "2"},
                                  {"infinite_radius", false}}},
          {"tag", json::object{
              {"a", json::object{{"domain", "rational"},
                                  {"canonical", a}}},
              {"b", json::object{{"domain", "rational"},
                                  {"canonical", "0"}}}}},
          {"prescriptions", prescriptions(sign)},
          {"checkpoint_identity", checkpoint}}}});
  if (response.at("status") != "ok")
    throw std::runtime_error("local.solve: " + json::serialize(response));
  return std::string(response.at("local").as_string());
}

json::object topology(std::optional<std::int32_t> sign) {
  json::array sheets;
  if (sign.has_value())
    sheets.push_back(json::object{{"factor_exact", "x"},
                                  {"sign", *sign}});
  return json::object{{"singular_points", json::array{}},
                      {"boundary_points", json::array{}},
                      {"complex_projections", json::array{}},
                      {"branch_sheets", std::move(sheets)}};
}

json::object arm(const std::string& endpoint,
                 const std::string& chart,
                 std::optional<std::int32_t> sign) {
  return json::object{{"from_exact", "0"}, {"to_exact", endpoint},
                      {"charts", json::array{chart}},
                      {"topology", topology(sign)}};
}

json::object make_plan(const std::string& session,
                       const std::string& chart,
                       const std::string& checkpoint,
                       std::optional<std::int32_t> sign) {
  return request(json::object{
      {"schema", 2}, {"op", "tile.plan"}, {"session", session},
      {"checkpoint_identity", checkpoint}, {"division_order", 3},
      {"lower", arm("-1/2", chart, sign)},
      {"upper", arm("1/2", chart, sign)}});
}

json::object integrate(const std::string& session,
                       const std::string& plan,
                       const std::string& plan_checkpoint,
                       const std::string& local,
                       const std::string& local_checkpoint,
                       const std::string& arm_name,
                       const std::string& checkpoint) {
  return request(json::object{
      {"schema", 2}, {"op", "integration.line"},
      {"session", session}, {"tile_plan", plan}, {"local", local},
      {"arm", arm_name}, {"tile", 0},
      {"epsilon", json::object{{"min", 0}, {"max", 0}}},
      {"source_checkpoint_identity", local_checkpoint},
      {"tile_plan_checkpoint_identity", plan_checkpoint},
      {"checkpoint_identity", checkpoint}});
}

json::object export_line(const std::string& session,
                         const json::object& integrated,
                         const std::string& checkpoint) {
  return request(json::object{
      {"schema", 2}, {"op", "integration.export"},
      {"session", session},
      {"line", std::string(integrated.at("line").as_string())},
      {"checkpoint_identity", checkpoint}, {"output_digits", 40}});
}

double imaginary_midpoint(const json::object& exported) {
  const auto& encoded = exported.at("value").as_object()
      .at("coefficients").as_array().front().as_array();
  return std::stod(std::string(encoded.at(1).as_string()));
}

}  // namespace

int main() {
  const auto checkpoint_path =
      (std::filesystem::temp_directory_path() /
       ("diffexp2-analytic-prescription-" +
        std::to_string(::getpid()) + ".de2cp")).string();
  std::string session;
  std::string restored_session;
  bool ok = false;
  try {
    const auto created = request(json::object{
        {"schema", 2}, {"op", "session.create"},
        {"domain", "rational"}, {"precision_bits", 256},
        {"output_digits", 40}, {"chart_capacity", 8},
        {"local_capacity", 8}, {"tile_plan_capacity", 8},
        {"line_result_capacity", 8}});
    session = std::string(created.at("session").as_string());

    const auto plus_chart = prepare_chart(
        session, "analytic-plus-chart", "1", 1);
    const auto minus_chart = prepare_chart(
        session, "analytic-minus-chart", "1", -1);
    const auto plain_chart = prepare_chart(
        session, "analytic-plain-chart", "1", std::nullopt);
    const auto reversed_chart = prepare_chart(
        session, "analytic-reversed-chart", "-1", 1);
    const auto tangential_chart = prepare_chart(
        session, "analytic-tangential-chart", "1", 1, 2);

    const auto plus_local = solve_local(
        session, plus_chart, "analytic-plus-local", "1", "1/2", 1);
    const auto minus_local = solve_local(
        session, minus_chart, "analytic-minus-local", "1", "1/2", -1);
    const auto missing_local = solve_local(
        session, plain_chart, "analytic-missing-local", "1", "1/2",
        std::nullopt);
    const auto single_valued_local = solve_local(
        session, plain_chart, "analytic-single-valued-local", "1", "0",
        std::nullopt);
    const auto reversed_local = solve_local(
        session, reversed_chart, "analytic-reversed-local", "-1", "1/2",
        1);

    const auto plus_plan = make_plan(
        session, plus_chart, "analytic-plus-plan", 1);
    const auto minus_plan = make_plan(
        session, minus_chart, "analytic-minus-plan", -1);
    const auto plain_plan = make_plan(
        session, plain_chart, "analytic-plain-plan", std::nullopt);
    const auto reversed_plan = make_plan(
        session, reversed_chart, "analytic-reversed-plan", 1);
    const auto tangential_plan = make_plan(
        session, tangential_chart, "analytic-tangential-plan", 1);
    for (const auto* plan : {&plus_plan, &minus_plan, &plain_plan,
                             &reversed_plan})
      if (plan->at("status") != "ok")
        throw std::runtime_error("tile.plan: " + json::serialize(*plan));

    const auto plus = integrate(
        session, std::string(plus_plan.at("tile_plan").as_string()),
        "analytic-plus-plan", plus_local, "analytic-plus-local", "lower",
        "analytic-plus-line");
    const auto minus = integrate(
        session, std::string(minus_plan.at("tile_plan").as_string()),
        "analytic-minus-plan", minus_local, "analytic-minus-local", "lower",
        "analytic-minus-line");
    const auto missing = integrate(
        session, std::string(plain_plan.at("tile_plan").as_string()),
        "analytic-plain-plan", missing_local, "analytic-missing-local",
        "lower", "analytic-missing-line");
    const auto single_valued = integrate(
        session, std::string(plain_plan.at("tile_plan").as_string()),
        "analytic-plain-plan", single_valued_local,
        "analytic-single-valued-local", "lower",
        "analytic-single-valued-line");
    // Physical +i0 with x=-t becomes the lower local rim.  Integrating the
    // upper physical arm therefore exercises a negative local coordinate.
    const auto reversed = integrate(
        session, std::string(reversed_plan.at("tile_plan").as_string()),
        "analytic-reversed-plan", reversed_local, "analytic-reversed-local",
        "upper", "analytic-reversed-line");

    const auto plus_export = export_line(
        session, plus, "analytic-plus-line");
    const auto minus_export = export_line(
        session, minus, "analytic-minus-line");
    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"}, {"session", session},
        {"path", checkpoint_path},
        {"checkpoint_identity", "analytic-prescription-session-v1"}});
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
    session.clear();

    const auto restored = request(json::object{
        {"schema", 2}, {"op", "checkpoint.restore"},
        {"path", checkpoint_path},
        {"expected_identity", "analytic-prescription-session-v1"}});
    if (restored.at("status") == "ok")
      restored_session = std::string(restored.at("session").as_string());
    const auto restored_single_valued = restored.at("status") == "ok"
        ? request(json::object{
              {"schema", 2}, {"op", "integration.stats"},
              {"session", restored_session},
              {"line", std::string(single_valued.at("line").as_string())}})
        : json::object{{"status", "error"}};

    const auto plus_imaginary = imaginary_midpoint(plus_export);
    const auto minus_imaginary = imaginary_midpoint(minus_export);
    ok = plus.at("status") == "ok" && plus.at("effective_rim") == 1 &&
        minus.at("status") == "ok" && minus.at("effective_rim") == -1 &&
        plus_imaginary < 0.0 && minus_imaginary > 0.0 &&
        std::abs(plus_imaginary + minus_imaginary) < 1e-30 &&
        missing.at("status") == "error" && missing.at("id") == "E3" &&
        std::string(missing.at("detail").as_string()).find(
            "no derivable imaginary sign") != std::string::npos &&
        single_valued.at("status") == "ok" &&
        single_valued.at("effective_rim").is_null() &&
        reversed.at("status") == "ok" &&
        reversed.at("effective_rim") == -1 &&
        tangential_plan.at("status") == "error" &&
        std::string(tangential_plan.at("detail").as_string()).find(
            "even-multiplicity tangential") != std::string::npos &&
        saved.at("status") == "ok" && restored.at("status") == "ok" &&
        restored_single_valued.at("status") == "ok" &&
        restored_single_valued.at("effective_rim").is_null();

    if (!ok) {
      std::cerr << "+i0: " << json::serialize(plus) << '\n'
                << "-i0: " << json::serialize(minus) << '\n'
                << "missing: " << json::serialize(missing) << '\n'
                << "single-valued: " << json::serialize(single_valued)
                << '\n' << "reversed: " << json::serialize(reversed) << '\n'
                << "tangential: " << json::serialize(tangential_plan) << '\n'
                << "save: " << json::serialize(saved) << '\n'
                << "restore: " << json::serialize(restored) << '\n'
                << "restored single-valued: "
                << json::serialize(restored_single_valued) << '\n';
    }
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << '\n';
  }

  if (!session.empty())
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", session}});
  if (!restored_session.empty())
    (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                               {"session", restored_session}});
  std::remove(checkpoint_path.c_str());
  std::cout << (ok ? "PASS" : "FAIL")
            << ": persistent analytic-prescription line protocol\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
