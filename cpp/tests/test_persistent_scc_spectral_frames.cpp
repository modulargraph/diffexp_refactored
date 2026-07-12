#include "diffexp2/checkpoint.hpp"
#include "diffexp2/json_codec.hpp"

#include <boost/json.hpp>

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace json = boost::json;

namespace {

constexpr std::uint32_t kFrameWidth = 4;
constexpr std::uint32_t kWorkTaylorOrder = 7;

using Matrix2 = std::array<std::array<std::string, 2>, 2>;

json::object regular_run(bool seed);

json::object request(json::object value) {
  return json::parse(
      diffexp2::run_recurrence_json(json::serialize(value))).as_object();
}

json::value nested(json::array value) {
  return json::value(std::move(value));
}

json::object exact_cell(const std::string& exact, bool zero) {
  return json::object{{"exact", exact}, {"proven_zero", zero}};
}

json::object epsilon_rational(const std::string& value) {
  return json::object{{"zero", value == "0"},
                      {"valuation", value == "0" ? 0 : 0},
                      {"numerator", json::array{value}},
                      {"denominator", json::array{"1"}}};
}

json::object physical_entry(std::uint32_t row, std::uint32_t column,
                            const std::string& value) {
  return json::object{{"r", row}, {"c", column},
                      {"v", epsilon_rational(value)}};
}

json::array constant_kernels(const std::string& value) {
  json::array kernels;
  for (std::uint32_t epsilon = 0; epsilon < kFrameWidth; ++epsilon) {
    json::array taylor;
    for (std::uint32_t n = 0; n <= kWorkTaylorOrder; ++n)
      taylor.emplace_back(epsilon == 0 && n == 0 ? value : "0");
    kernels.push_back(nested(std::move(taylor)));
  }
  return kernels;
}

json::array linear_kernels(const std::string& value) {
  json::array kernels;
  for (std::uint32_t epsilon = 0; epsilon < kFrameWidth; ++epsilon) {
    json::array taylor;
    for (std::uint32_t n = 0; n <= kWorkTaylorOrder; ++n)
      taylor.emplace_back(epsilon == 0 && n == 1 ? value : "0");
    kernels.push_back(nested(std::move(taylor)));
  }
  return kernels;
}

json::object source_transform(const std::string& role,
                              const std::string& v,
                              const std::string& vinv,
                              const std::string& domain,
                              bool identity = false) {
  const auto v_identity = role + ":V=" + v;
  const auto vinv_identity = role + ":VInv=" + vinv;
  const auto det_identity = role + ":det=" + v;
  const auto entry_identity = role + ":entry=" + vinv;
  const auto producer_identity = json::serialize(json::object{
      {"schema",
       "diffexp2-scc-spectral-source-transform-identity-v1"},
      {"state_basis", "reduced-g-after-spectral-assembly"},
      {"target_recurrence_basis", "spectral-u"},
      {"dimension", 1}, {"identity", identity},
      {"epsilon_unimodular", true}, {"det_epsilon_valuation", 0},
      {"v_exact_identity", v_identity},
      {"vinv_exact_identity", vinv_identity},
      {"det_exact_identity", det_identity},
      {"source_window", json::object{
          {"epsilon_min", 0}, {"epsilon_complete_max", 3},
          {"taylor_complete_max", kWorkTaylorOrder}}},
      {"serialization", json::object{
          {"domain", domain}, {"symbols", json::array{}}}},
      {"entries", json::array{json::object{
          {"row", 0}, {"column", 0},
          {"exact_entry", entry_identity}, {"epsilon_shift", 0},
          {"center_pole_order", 0}}}}});
  return json::object{
      {"schema", "diffexp2-scc-spectral-source-transform-v1"},
      {"rows", 1}, {"columns", 1}, {"identity", identity},
      {"epsilon_unimodular", true}, {"det_epsilon_valuation", 0},
      {"v_exact_identity", v_identity},
      {"vinv_exact_identity", vinv_identity},
      {"det_exact_identity", det_identity},
      {"exact_identity", producer_identity},
      {"domain", domain}, {"symbols", json::array{}},
      {"entries", json::array{json::object{
          {"row", 0}, {"column", 0},
          {"exact_entry", entry_identity},
          {"multiplier", json::object{
              {"epsilon_shift", 0}, {"center_pole_order", 0},
              {"kernels", constant_kernels(vinv)},
              {"exact_identity", entry_identity},
              {"proven_zero", false}}}}}}};
}

json::object gauge_transform(const std::string& frame,
                             const std::string& role,
                             const std::string& gauge,
                             const std::string& inverse,
                             const std::string& multiplier,
                             const std::string& domain) {
  const auto gauge_identity = frame + ":Gauge=" + gauge;
  const auto inverse_identity = frame + ":GaugeInverse=" + inverse;
  const auto determinant_identity = frame + ":det=" + gauge;
  const auto entry_identity = frame + ":" + role + "=" + multiplier;
  const auto proof = json::serialize(json::object{
      {"schema", "diffexp2-scc-gauge-transform-identity-v1"},
      {"role", role}, {"dimension", 1}, {"identity", false},
      {"gauge_exact_identity", gauge_identity},
      {"gauge_inverse_exact_identity", inverse_identity},
      {"gauge_det_exact_identity", determinant_identity},
      {"source_window", json::object{{"epsilon_min", 0},
          {"epsilon_complete_max", 3},
          {"taylor_complete_max", kWorkTaylorOrder}}},
      {"entries", json::array{json::object{{"row", 0}, {"column", 0},
          {"exact_entry", entry_identity}, {"epsilon_shift", 0},
          {"center_pole_order", 0}}}}});
  return json::object{
      {"schema", "diffexp2-scc-gauge-transform-v1"}, {"role", role},
      {"rows", 1}, {"columns", 1}, {"identity", false},
      {"gauge_exact_identity", gauge_identity},
      {"gauge_inverse_exact_identity", inverse_identity},
      {"gauge_det_exact_identity", determinant_identity},
      {"exact_identity", proof}, {"domain", domain},
      {"symbols", json::array{}},
      {"entries", json::array{json::object{{"row", 0}, {"column", 0},
          {"exact_entry", entry_identity},
          {"multiplier", json::object{{"epsilon_shift", 0},
              {"center_pole_order", 0},
              {"kernels", constant_kernels(multiplier)},
              {"exact_identity", entry_identity},
              {"proven_zero", false}}}}}}};
}

json::object dense_source_transform(const std::string& role,
                                    const Matrix2& v,
                                    const Matrix2& vinv,
                                    const std::string& determinant,
                                    const std::string& domain) {
  const auto v_identity = role + ":dense-V";
  const auto vinv_identity = role + ":dense-VInv";
  json::array entries;
  json::array identity_entries;
  for (std::uint32_t row = 0; row < 2; ++row) {
    for (std::uint32_t column = 0; column < 2; ++column) {
      if (vinv[row][column] == "0") continue;
      const auto exact_entry = role + ":VInv:" +
          std::to_string(row) + ":" + std::to_string(column) + "=" +
          vinv[row][column];
      entries.push_back(json::object{
          {"row", row}, {"column", column},
          {"exact_entry", exact_entry},
          {"multiplier", json::object{
              {"epsilon_shift", 0}, {"center_pole_order", 0},
              {"kernels", constant_kernels(vinv[row][column])},
              {"exact_identity", exact_entry},
              {"proven_zero", false}}}});
      identity_entries.push_back(json::object{
          {"row", row}, {"column", column},
          {"exact_entry", exact_entry}, {"epsilon_shift", 0},
          {"center_pole_order", 0}});
    }
  }
  const auto producer_identity = json::serialize(json::object{
      {"schema",
       "diffexp2-scc-spectral-source-transform-identity-v1"},
      {"state_basis", "reduced-g-after-spectral-assembly"},
      {"target_recurrence_basis", "spectral-u"},
      {"dimension", 2}, {"identity", false},
      {"epsilon_unimodular", true}, {"det_epsilon_valuation", 0},
      {"v_exact_identity", v_identity},
      {"vinv_exact_identity", vinv_identity},
      {"det_exact_identity", role + ":det=" + determinant},
      {"source_window", json::object{
          {"epsilon_min", 0}, {"epsilon_complete_max", 3},
          {"taylor_complete_max", kWorkTaylorOrder}}},
      {"serialization", json::object{
          {"domain", domain}, {"symbols", json::array{}}}},
      {"entries", identity_entries}});
  return json::object{
      {"schema", "diffexp2-scc-spectral-source-transform-v1"},
      {"rows", 2}, {"columns", 2}, {"identity", false},
      {"epsilon_unimodular", true}, {"det_epsilon_valuation", 0},
      {"v_exact_identity", v_identity},
      {"vinv_exact_identity", vinv_identity},
      {"det_exact_identity", role + ":det=" + determinant},
      {"exact_identity", producer_identity},
      {"domain", domain}, {"symbols", json::array{}},
      {"entries", std::move(entries)}};
}

json::object matrix_shift(const Matrix2& matrix) {
  json::array entries;
  for (std::uint32_t row = 0; row < 2; ++row)
    for (std::uint32_t column = 0; column < 2; ++column)
      if (matrix[row][column] != "0")
        entries.push_back(nested(json::array{
            row, column, matrix[row][column]}));
  return json::object{{"s", 0}, {"e", std::move(entries)}};
}

json::array matrix_valuations(const Matrix2& matrix) {
  json::array values;
  for (const auto& row : matrix)
    for (const auto& value : row)
      values.push_back(value == "0" ? json::value(nullptr) : json::value(0));
  return values;
}

json::object dense_problem(const std::string& domain,
                           const std::string& role,
                           const Matrix2& recurrence,
                           const Matrix2& assembly) {
  json::object problem{
      {"domain", domain}, {"d", 2}, {"fb", 0}, {"w", kFrameWidth},
      {"d_lags", json::array{nested(json::array{
          json::object{{"s", 0}, {"v", "1"}}})}},
      {"denominators", json::array{}},
      {"nhat_lags", json::array{
          json::object{{"poly", json::array{}},
                       {"rat", json::array{}},
                       {"val", json::array{
                           nullptr, nullptr, nullptr, nullptr}}},
          json::object{{"poly", json::array{matrix_shift(recurrence)}},
                       {"rat", json::array{}},
                       {"val", matrix_valuations(recurrence)}}}},
      {"d0_inverse", "1"},
      {"blocks", json::array{nested(json::array{0}),
                               nested(json::array{1})}},
      {"assembly", json::object{
          {"identity", false},
          {"exact_identity", role + ":dense-V"},
          {"poly", json::array{matrix_shift(assembly)}},
          {"rat", json::array{}},
          {"val", matrix_valuations(assembly)}}},
      {"chop_digits", 0}};
  if (domain == "acb") problem["precision_bits"] = 256;
  return problem;
}

std::string prepare_dense_chart(const std::string& session,
                                const std::string& domain,
                                const std::string& role,
                                const Matrix2& recurrence,
                                const Matrix2& assembly) {
  const Matrix2 physical{{{{"0", "1"}}, {{"1", "0"}}}};
  json::array principal;
  for (const auto& row : physical)
    principal.push_back(nested(json::array{
        exact_cell(row[0], row[0] == "0"),
        exact_cell(row[1], row[1] == "0")}));
  auto response = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", domain + ":" + role + ":dense-key"},
      {"identity", domain + ":" + role + ":dense-principal"},
      {"analytic", json::object{
          {"geometry", json::object{
              {"center_exact", "0"}, {"scale_exact", "1"},
              {"radius_exact", "2"}, {"infinite_radius", false},
              {"prescriptions", json::array{}}}},
          {"principal_matrix", std::move(principal)},
          {"native_scc_capabilities", json::object{
              {"regular", true}, {"identity_gauge", true},
              {"identity_v", false}, {"epsilon_unimodular_v", true},
              {"no_pseudo", true}}}}},
      {"scc", json::object{
          {"components", json::array{nested(json::array{0, 1})}},
          {"structural_edges", json::array{
              nested(json::array{0, 1}), nested(json::array{1, 0})}},
          {"condensation_edges", json::array{}},
          {"topological_order", json::array{0}},
          {"coupling_depth", 0}}},
      {"problem", dense_problem(domain, role, recurrence, assembly)}});
  if (response.at("status") != "ok")
    throw std::runtime_error(
        "dense chart.prepare: " + json::serialize(response));
  return std::string(response.at("chart").as_string());
}

json::object dense_run(bool seed) {
  auto run = regular_run(seed);
  json::array initial;
  for (std::uint32_t component = 0; component < 2; ++component)
    for (std::uint32_t epsilon = 0; epsilon < kFrameWidth; ++epsilon)
      initial.emplace_back(seed && component == 0 && epsilon == 0
                               ? "1" : "0");
  run["initial"] = std::move(initial);
  run["initial_validity"] = seed ? json::array{3, 3}
                                  : json::array{nullptr, nullptr};
  auto& schedule = run.at("schedule").as_array();
  for (auto& raw_row : schedule) {
    const auto step = raw_row.as_array().front();
    raw_row = json::array{step, step};
  }
  return run;
}

json::object scalar_problem(const std::string& domain,
                            const std::string& role,
                            const std::string& v,
                            bool identity) {
  json::object assembly{
      {"identity", identity}, {"exact_identity", role + ":V=" + v},
      {"poly", identity ? json::array{} : json::array{json::object{
          {"s", 0},
          {"e", json::array{nested(json::array{0, 0, v})}}}}},
      {"rat", json::array{}}, {"val", json::array{0}}};
  json::object problem{
      {"domain", domain}, {"d", 1}, {"fb", 0}, {"w", kFrameWidth},
      {"d_lags", json::array{nested(json::array{
          json::object{{"s", 0}, {"v", "1"}}})}},
      {"denominators", json::array{}},
      {"nhat_lags", json::array{json::object{
          {"poly", json::array{}}, {"rat", json::array{}},
          {"val", json::array{nullptr}}}}},
      {"d0_inverse", "1"},
      {"blocks", json::array{nested(json::array{0})}},
      {"assembly", std::move(assembly)},
      {"chop_digits", 0}};
  if (domain == "acb") problem["precision_bits"] = 256;
  return problem;
}

std::string prepare_chart(const std::string& session,
                          const std::string& domain,
                          const std::string& role,
                          const std::string& v,
                          bool identity = false,
                          bool identity_gauge = true) {
  auto response = request(json::object{
      {"schema", 2}, {"op", "chart.prepare"}, {"session", session},
      {"key", domain + ":" + role + ":key"},
      {"identity", domain + ":" + role + ":principal"},
      {"analytic", json::object{
          {"geometry", json::object{
              {"center_exact", "0"}, {"scale_exact", "1"},
              {"radius_exact", "2"}, {"infinite_radius", false},
              {"prescriptions", json::array{}}}},
          {"principal_matrix", json::array{
              nested(json::array{exact_cell("0", true)})}},
          {"native_scc_capabilities", json::object{
              {"regular", true}, {"identity_gauge", identity_gauge},
              {"exact_gauge", true},
              {"identity_v", identity},
              {"epsilon_unimodular_v", true},
              {"no_pseudo", true}}}}},
      {"scc", json::object{
          {"components", json::array{nested(json::array{0})}},
          {"structural_edges", json::array{}},
          {"condensation_edges", json::array{}},
          {"topological_order", json::array{0}},
          {"coupling_depth", 0}}},
      {"problem", scalar_problem(domain, role, v, identity)}});
  if (response.at("status") != "ok")
    throw std::runtime_error("chart.prepare: " + json::serialize(response));
  return std::string(response.at("chart").as_string());
}

json::object regular_run(bool seed) {
  json::array shifts;
  json::array schedule;
  for (std::uint32_t n = 0; n <= kWorkTaylorOrder; ++n) {
    shifts.emplace_back(std::to_string(n));
    schedule.push_back(json::array{json::object{
        {"case", n == 0 ? "R" : "T"},
        {"da", std::to_string(n)}, {"db", "0"}}});
  }
  json::array initial;
  for (std::uint32_t epsilon = 0; epsilon < kFrameWidth; ++epsilon)
    initial.emplace_back(seed && epsilon == 0 ? "1" : "0");
  return json::object{
      {"nmax", kWorkTaylorOrder}, {"p", 0}, {"has_initial", seed},
      {"adaptive_probe", false}, {"a_target", "0"},
      {"b_target", "0"}, {"a_shift_min", 0},
      {"a_shifts", std::move(shifts)}, {"schedule", std::move(schedule)},
      {"initial", std::move(initial)},
      {"initial_validity", seed ? json::array{3}
                                  : json::array{nullptr}},
      {"source", nullptr}, {"return_u", false}};
}

json::object metadata(const std::string& checkpoint) {
  return json::object{
      {"chart", json::object{
          {"center_exact", "0"}, {"scale_exact", "1"},
          {"radius", "2"}, {"infinite_radius", false}}},
      {"tag", json::object{
          {"a", json::object{{"domain", "rational"},
                               {"canonical", "0"}}},
          {"b", json::object{{"domain", "rational"},
                               {"canonical", "0"}}}}},
      {"prescriptions", json::array{}},
      {"checkpoint_identity", checkpoint}};
}

json::object parent_physical_ode(const std::string& identity) {
  return json::object{
      {"schema", "diffexp2-physical-cleared-ode-v1"},
      {"basis", "physical-original-master"},
      {"theta_coordinate", "local-t"},
      {"owner_signature_identity", identity},
      {"payload_identity", "de2-physical-ode-" + identity},
      {"q", json::array{epsilon_rational("1")}},
      {"c", json::array{
          nested(json::array{}),
          nested(json::array{physical_entry(1, 0, "5")})}}};
}

json::object manifest(const std::string& domain,
                      const std::string& first,
                      const std::string& second,
                      const std::string& identity) {
  auto source = source_transform("source", "2", "1/2", domain);
  auto target = source_transform("target", "3", "1/3", domain);
  json::object coupling{
      {"source_block", 0}, {"target_block", 1},
      {"source_vertices", json::array{0}},
      {"target_vertices", json::array{1}},
      {"rows", 1}, {"columns", 1},
      {"exact_identity", "physical-C=5*t"},
      {"domain", domain}, {"symbols", json::array{}},
      {"entries", json::array{json::object{
          {"row", 0}, {"column", 0},
          {"source_vertex", 0}, {"target_vertex", 1},
          {"exact_original_entry", "5"},
          {"exact_theta_entry", "5*t"},
          {"multiplier", json::object{
              {"epsilon_shift", 0}, {"center_pole_order", 0},
              {"kernels", linear_kernels("5")},
              {"exact_identity", "5*t"},
              {"proven_zero", false}}}}}}};
  return json::object{
      {"identity", identity},
      {"parent", json::object{
          {"dimension", 2},
          {"exact_system_record", json::array{
              nested(json::array{exact_cell("0", true),
                                  exact_cell("0", true)}),
              nested(json::array{exact_cell("5", false),
                                  exact_cell("0", true)})}},
          {"exact_theta_record", json::array{
              nested(json::array{exact_cell("0", true),
                                  exact_cell("0", true)}),
              nested(json::array{exact_cell("5*t", false),
                                  exact_cell("0", true)})}},
          {"chart", json::object{
              {"center_exact", "0"}, {"scale_exact", "1"},
              {"radius_exact", "2"}, {"infinite_radius", false},
              {"prescriptions", json::array{}}}},
          {"scc", json::object{
              {"components", json::array{nested(json::array{0}),
                                           nested(json::array{1})}},
              {"structural_edges", json::array{
                  nested(json::array{0, 1})}},
              {"condensation_edges", json::array{
                  nested(json::array{0, 1})}},
              {"topological_order", json::array{0, 1}},
              {"coupling_depth", 1}}},
          {"execution", json::object{
              {"mode", "BlockSequentialStrict"},
              {"work_t_order", kWorkTaylorOrder}}},
          {"work_contract", json::object{
              {"work_min", 0}, {"requested_min", 0},
              {"requested_max", 2}, {"work_complete_max", 3},
              {"public_t_order", 1}, {"wolfram_coupling_depth", 2}}}}},
      {"blocks", json::array{
          json::object{
              {"block", 0}, {"vertices", json::array{0}},
              {"chart", first},
              {"principal_identity", domain + ":source:principal"},
              {"regular", true}, {"identity_gauge", true},
              {"identity_v", false}, {"epsilon_unimodular_v", true},
              {"source_transform", std::move(source)},
              {"no_pseudo", true}},
          json::object{
              {"block", 1}, {"vertices", json::array{1}},
              {"chart", second},
              {"principal_identity", domain + ":target:principal"},
              {"regular", true}, {"identity_gauge", true},
              {"identity_v", false}, {"epsilon_unimodular_v", true},
              {"source_transform", std::move(target)},
              {"no_pseudo", true}}}},
      {"couplings", json::array{std::move(coupling)}},
      {"physical_ode", parent_physical_ode(identity)}};
}

std::string theta_entry(const std::string& value) {
  if (value == "0") return "0";
  if (value == "1") return "t";
  if (value == "-1") return "-t";
  return value + "*t";
}

json::array exact_matrix_record(
    const std::array<std::array<std::string, 4>, 4>& matrix,
    bool theta) {
  json::array rows;
  for (const auto& row : matrix) {
    json::array columns;
    for (const auto& value : row) {
      const auto exact = theta ? theta_entry(value) : value;
      columns.push_back(exact_cell(exact, value == "0"));
    }
    rows.push_back(nested(std::move(columns)));
  }
  return rows;
}

json::object dense_parent_physical_ode(
    const std::string& identity,
    const std::array<std::array<std::string, 4>, 4>& matrix) {
  json::array entries;
  for (std::uint32_t row = 0; row < 4; ++row)
    for (std::uint32_t column = 0; column < 4; ++column)
      if (matrix[row][column] != "0")
        entries.push_back(physical_entry(
            row, column, matrix[row][column]));
  return json::object{
      {"schema", "diffexp2-physical-cleared-ode-v1"},
      {"basis", "physical-original-master"},
      {"theta_coordinate", "local-t"},
      {"owner_signature_identity", identity},
      {"payload_identity", "de2-physical-ode-" + identity},
      {"q", json::array{epsilon_rational("1")}},
      {"c", json::array{nested(json::array{}),
                          nested(std::move(entries))}}};
}

json::object dense_manifest(const std::string& domain,
                            const std::string& source_chart,
                            const std::string& target_chart,
                            const std::string& identity) {
  const Matrix2 source_v{{{{"1", "1"}}, {{"1", "2"}}}};
  const Matrix2 source_vinv{{{{"2", "-1"}}, {{"-1", "1"}}}};
  const Matrix2 target_v{{{{"1", "2"}}, {{"3", "5"}}}};
  const Matrix2 target_vinv{{{{"-5", "2"}}, {{"3", "-1"}}}};
  const Matrix2 coupling{{{{"2", "1"}}, {{"-1", "3"}}}};
  std::array<std::array<std::string, 4>, 4> physical{};
  for (auto& row : physical) row.fill("0");
  physical[0][1] = "1";
  physical[1][0] = "1";
  physical[2][3] = "1";
  physical[3][2] = "1";
  for (std::uint32_t row = 0; row < 2; ++row)
    for (std::uint32_t column = 0; column < 2; ++column)
      physical[row + 2][column] = coupling[row][column];

  json::array coupling_entries;
  for (std::uint32_t row = 0; row < 2; ++row) {
    for (std::uint32_t column = 0; column < 2; ++column) {
      coupling_entries.push_back(json::object{
          {"row", row}, {"column", column},
          {"source_vertex", column}, {"target_vertex", row + 2},
          {"exact_original_entry", coupling[row][column]},
          {"exact_theta_entry", theta_entry(coupling[row][column])},
          {"multiplier", json::object{
              {"epsilon_shift", 0}, {"center_pole_order", 0},
              {"kernels", linear_kernels(coupling[row][column])},
              {"exact_identity", theta_entry(coupling[row][column])},
              {"proven_zero", false}}}});
    }
  }
  return json::object{
      {"identity", identity},
      {"parent", json::object{
          {"dimension", 4},
          {"exact_system_record", exact_matrix_record(physical, false)},
          {"exact_theta_record", exact_matrix_record(physical, true)},
          {"chart", json::object{
              {"center_exact", "0"}, {"scale_exact", "1"},
              {"radius_exact", "2"}, {"infinite_radius", false},
              {"prescriptions", json::array{}}}},
          {"scc", json::object{
              {"components", json::array{
                  nested(json::array{0, 1}),
                  nested(json::array{2, 3})}},
              {"structural_edges", json::array{
                  nested(json::array{0, 1}),
                  nested(json::array{0, 2}),
                  nested(json::array{0, 3}),
                  nested(json::array{1, 0}),
                  nested(json::array{1, 2}),
                  nested(json::array{1, 3}),
                  nested(json::array{2, 3}),
                  nested(json::array{3, 2})}},
              {"condensation_edges", json::array{
                  nested(json::array{0, 1})}},
              {"topological_order", json::array{0, 1}},
              {"coupling_depth", 1}}},
          {"execution", json::object{
              {"mode", "BlockSequentialStrict"},
              {"work_t_order", kWorkTaylorOrder}}},
          {"work_contract", json::object{
              {"work_min", 0}, {"requested_min", 0},
              {"requested_max", 2}, {"work_complete_max", 3},
              {"public_t_order", 1},
              {"wolfram_coupling_depth", 2}}}}},
      {"blocks", json::array{
          json::object{
              {"block", 0}, {"vertices", json::array{0, 1}},
              {"chart", source_chart},
              {"principal_identity", domain + ":dense-source:dense-principal"},
              {"regular", true}, {"identity_gauge", true},
              {"identity_v", false}, {"epsilon_unimodular_v", true},
              {"source_transform", dense_source_transform(
                  "dense-source", source_v, source_vinv, "1", domain)},
              {"no_pseudo", true}},
          json::object{
              {"block", 1}, {"vertices", json::array{2, 3}},
              {"chart", target_chart},
              {"principal_identity", domain + ":dense-target:dense-principal"},
              {"regular", true}, {"identity_gauge", true},
              {"identity_v", false}, {"epsilon_unimodular_v", true},
              {"source_transform", dense_source_transform(
                  "dense-target", target_v, target_vinv, "-1", domain)},
              {"no_pseudo", true}}}},
      {"couplings", json::array{json::object{
          {"source_block", 0}, {"target_block", 1},
          {"source_vertices", json::array{0, 1}},
          {"target_vertices", json::array{2, 3}},
          {"rows", 2}, {"columns", 2},
          {"exact_identity", "dense-physical-cross-C"},
          {"domain", domain}, {"symbols", json::array{}},
          {"entries", std::move(coupling_entries)}}}},
      {"physical_ode", dense_parent_physical_ode(identity, physical)}};
}

json::object prepare_scc(const std::string& session,
                         const std::string& key,
                         json::object value) {
  value["schema"] = 2;
  value["op"] = "scc.prepare";
  value["session"] = session;
  value["key"] = key;
  return request(std::move(value));
}

double real_midpoint(const json::value& value) {
  return std::stod(std::string(value.as_array().front().as_string()));
}

bool run_dense_domain(const std::string& domain) {
  const Matrix2 source_v{{{{"1", "1"}}, {{"1", "2"}}}};
  const Matrix2 target_v{{{{"1", "2"}}, {{"3", "5"}}}};
  // V^-1 A V for A={{0,1},{1,0}}.  These non-symmetric matrices also make
  // a transpose/right-action bug visible independently of source assembly.
  const Matrix2 source_recurrence{{{{"1", "3"}}, {{"0", "-1"}}}};
  const Matrix2 target_recurrence{{{{"-13", "-21"}}, {{"8", "13"}}}};
  auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", domain},
      {"precision_bits", 256}, {"output_digits", 40},
      {"chart_capacity", 2}, {"scc_capacity", 1},
      {"local_capacity", 1}});
  if (created.at("status") != "ok") return false;
  const auto session = std::string(created.at("session").as_string());
  const auto source = prepare_dense_chart(
      session, domain, "dense-source", source_recurrence, source_v);
  const auto target = prepare_dense_chart(
      session, domain, "dense-target", target_recurrence, target_v);
  const auto identity = domain + ":dense-spectral-parent-v1";
  const auto prepared = prepare_scc(
      session, domain + ":dense-spectral-key",
      dense_manifest(domain, source, target, identity));
  json::object solved{{"status", "not-run"}};
  json::object evaluated{{"status", "not-run"}};
  json::object stats{{"status", "not-run"}};
  if (prepared.at("status") == "ok") {
    const auto scc = std::string(prepared.at("scc").as_string());
    solved = request(json::object{
        {"schema", 2}, {"op", "scc.solve_column"},
        {"session", session}, {"scc", scc},
        {"checkpoint_identity", domain + ":dense-spectral-column"},
        {"seed", json::object{
            {"block", 0}, {"run", dense_run(true)},
            {"metadata", metadata(domain + ":dense-seed")}}},
        {"targets", json::array{json::object{
            {"block", 1}, {"run", dense_run(false)},
            {"metadata", metadata(domain + ":dense-target")}}}}});
    if (solved.at("status") == "ok")
      evaluated = request(json::object{
          {"schema", 2}, {"op", "local.evaluate"},
          {"session", session}, {"local", solved.at("local")},
          {"point", json::object{{"exact", "1/2"}}},
          {"options", json::object{{"tail_estimate", false}}}});
    stats = request(json::object{
        {"schema", 2}, {"op", "scc.stats"},
        {"session", session}, {"scc", scc}});
  }

  bool values_ok = false;
  if (evaluated.at("status") == "ok") {
    const auto& coefficients = evaluated.at("value").as_object()
                                   .at("coefficients").as_array();
    // u_S(0)=e0, so g_S(0)=V_S e0=(1,1).  To public order one,
    // g_S(1/2)=(3/2,3/2).  The noncommuting physical cross block
    // C={{2,1},{-1,3}} gives C g_S(0)=(3,2), hence
    // g_T(1/2)=(3/2,1).  Omitting V_T^-1 gives (7/2,19/2), applying
    // VInv^T/right-oriented gives another vector, and an extra V_S gives
    // (7/2,7/2), so every orientation/double-transform bug fails here.
    const std::array<double, 4> expected{1.5, 1.5, 1.5, 1.0};
    values_ok = coefficients.size() == 12;
    for (std::size_t index = 0; index < expected.size() && values_ok;
         ++index)
      values_ok = std::abs(real_midpoint(coefficients[index]) -
                           expected[index]) < 1e-25;
    for (std::size_t index = expected.size();
         index < coefficients.size() && values_ok; ++index)
      values_ok = std::abs(real_midpoint(coefficients[index])) < 1e-25;
  }
  const bool ok = prepared.at("status") == "ok" &&
      solved.at("status") == "ok" && values_ok &&
      stats.at("execution_implemented") == true &&
      stats.at("execution_scope") ==
          (domain == "acb" ? "acb-regular-block-dag-column-v2"
                           : "exact-rational-regular-block-dag-column-v2") &&
      stats.at("scalar_block_dag_column_execution") == false &&
      stats.at("active_coupling_entries") == 4;
  if (!ok) {
    std::cerr << domain << " dense prepared: " << json::serialize(prepared)
              << '\n' << domain << " dense solved: "
              << json::serialize(solved) << '\n'
              << domain << " dense evaluated: "
              << json::serialize(evaluated) << '\n'
              << domain << " dense stats: " << json::serialize(stats)
              << '\n';
  }
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                              {"session", session}});
  return ok;
}

bool run_domain(const std::string& domain, bool test_rejections) {
  auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", domain},
      {"precision_bits", 256}, {"output_digits", 40},
      {"chart_capacity", 4}, {"scc_capacity", 4},
      {"local_capacity", 2}});
  if (created.at("status") != "ok") return false;
  const auto session = std::string(created.at("session").as_string());
  const auto first = prepare_chart(session, domain, "source", "2");
  const auto second = prepare_chart(session, domain, "target", "3");
  const auto identity = domain + ":spectral-parent-v1";
  const auto prepared = prepare_scc(
      session, domain + ":spectral-key",
      manifest(domain, first, second, identity));

  json::object determinant_rejected{{"status", "skipped"}};
  json::object laurent_accepted{{"status", "skipped"}};
  json::object binding_rejected{{"status", "skipped"}};
  json::object identity_rejected{{"status", "skipped"}};
  if (test_rejections) {
    auto nonunimodular = manifest(
        domain, first, second, domain + ":nonunimodular-parent");
    nonunimodular.at("blocks").as_array()[1].as_object()
        .at("source_transform").as_object()["det_epsilon_valuation"] = 1;
    determinant_rejected = prepare_scc(
        session, domain + ":nonunimodular-key", std::move(nonunimodular));

    auto laurent = manifest(
        domain, first, second, domain + ":laurent-parent");
    auto& laurent_transform = laurent.at("blocks").as_array()[1]
        .as_object().at("source_transform").as_object();
    laurent_transform["det_epsilon_valuation"] = 1;
    auto laurent_identity = json::parse(
        laurent_transform.at("exact_identity").as_string()).as_object();
    laurent_identity["det_epsilon_valuation"] = 1;
    laurent_transform["exact_identity"] = json::serialize(laurent_identity);
    laurent_accepted = prepare_scc(
        session, domain + ":laurent-key", std::move(laurent));

    auto mismatched = manifest(
        domain, first, second, domain + ":mismatched-v-parent");
    mismatched.at("blocks").as_array()[1].as_object()
        .at("source_transform").as_object()["v_exact_identity"] =
        "attacker-V";
    binding_rejected = prepare_scc(
        session, domain + ":mismatched-v-key", std::move(mismatched));

    const auto identity_first = prepare_chart(
        session, domain, "identity-source", "1", true);
    const auto identity_second = prepare_chart(
        session, domain, "identity-target", "1", true);
    auto malformed_identity = manifest(
        domain, identity_first, identity_second,
        domain + ":malformed-identity-parent");
    auto& identity_blocks = malformed_identity.at("blocks").as_array();
    identity_blocks[0].as_object()["principal_identity"] =
        domain + ":identity-source:principal";
    identity_blocks[1].as_object()["principal_identity"] =
        domain + ":identity-target:principal";
    identity_blocks[0].as_object()["identity_v"] = true;
    identity_blocks[1].as_object()["identity_v"] = true;
    identity_blocks[0].as_object()["source_transform"] =
        source_transform("identity-source", "1", "1", domain, true);
    identity_blocks[1].as_object()["source_transform"] =
        source_transform("identity-target", "1", "1", domain, true);
    identity_blocks[0].as_object().at("source_transform").as_object()
        .at("entries").as_array()[0].as_object()
        .at("multiplier").as_object().at("kernels").as_array()[0]
        .as_array()[0] = "2";
    identity_rejected = prepare_scc(
        session, domain + ":malformed-identity-key",
        std::move(malformed_identity));
  }

  json::object solved{{"status", "not-run"}};
  json::object evaluated{{"status", "not-run"}};
  json::object stats{{"status", "not-run"}};
  if (prepared.at("status") == "ok") {
    const auto scc = std::string(prepared.at("scc").as_string());
    solved = request(json::object{
        {"schema", 2}, {"op", "scc.solve_column"},
        {"session", session}, {"scc", scc},
        {"checkpoint_identity", domain + ":spectral-column"},
        {"seed", json::object{
            {"block", 0}, {"run", regular_run(true)},
            {"metadata", metadata(domain + ":seed")}}},
        {"targets", json::array{json::object{
            {"block", 1}, {"run", regular_run(false)},
            {"metadata", metadata(domain + ":target")}}}}});
    if (solved.at("status") == "ok")
      evaluated = request(json::object{
          {"schema", 2}, {"op", "local.evaluate"},
          {"session", session}, {"local", solved.at("local")},
          {"point", json::object{{"exact", "1/2"}}},
          {"options", json::object{{"tail_estimate", false}}}});
    stats = request(json::object{
        {"schema", 2}, {"op", "scc.stats"},
        {"session", session}, {"scc", scc}});
  }

  bool checkpoint_ok = true;
  if (test_rejections && solved.at("status") == "ok") {
    const auto checkpoint =
        (std::filesystem::temp_directory_path() /
         "diffexp2-spectral-frame-checkpoint.bin").string();
    std::filesystem::remove(checkpoint);
    const auto saved = request(json::object{
        {"schema", 2}, {"op", "checkpoint.save"},
        {"session", session}, {"path", checkpoint},
        {"checkpoint_identity", "spectral-frame-checkpoint"}});
    if (saved.at("status") == "ok") {
      const auto container = diffexp2::checkpoint::read(checkpoint);
      const auto payload = json::parse(container.payload_json).as_object();
      const auto& retained = payload.at("prepared_scc").as_array();
      const json::object* spectral_record = nullptr;
      for (const auto& raw_record : retained) {
        const auto& record = raw_record.as_object();
        const auto& request = record.at("request").as_object();
        if (std::string(request.at("identity").as_string()) ==
            domain + ":spectral-parent-v1") {
          spectral_record = &record;
          break;
        }
      }
      if (spectral_record == nullptr) {
        checkpoint_ok = false;
      } else {
        const auto& blocks = spectral_record->at("request").as_object()
                                 .at("blocks").as_array();
        checkpoint_ok = blocks.size() == 2 &&
            blocks[0].as_object().at("source_transform").as_object()
                    .at("v_exact_identity") == "source:V=2" &&
            blocks[0].as_object().at("source_transform").as_object()
                    .at("vinv_exact_identity") == "source:VInv=1/2" &&
            blocks[1].as_object().at("source_transform").as_object()
                    .at("v_exact_identity") == "target:V=3" &&
            blocks[1].as_object().at("source_transform").as_object()
                    .at("vinv_exact_identity") == "target:VInv=1/3";
      }
    } else {
      checkpoint_ok = false;
    }
    std::filesystem::remove(checkpoint);
  }

  bool values_ok = false;
  if (evaluated.at("status") == "ok") {
    const auto& coefficients = evaluated.at("value").as_object()
                                   .at("coefficients").as_array();
    // g_source = V_S*1 = 2.  The physical cross source is
    // 5 t g_source = 10 t.  V_T^-1 gives (10/3)t to the target
    // recurrence and retained assembly V_T restores 10t, hence 5 at t=1/2.
    // Omitting V_T^-1 would give 15; multiplying V_S twice would give 10.
    values_ok = coefficients.size() == 6 &&
        std::abs(real_midpoint(coefficients[0]) - 2.0) < 1e-25 &&
        std::abs(real_midpoint(coefficients[1]) - 5.0) < 1e-25;
    for (std::size_t index = 2; index < coefficients.size(); ++index)
      values_ok = values_ok &&
          std::abs(real_midpoint(coefficients[index])) < 1e-25;
  }
  const bool rejection_ok = !test_rejections ||
      (determinant_rejected.at("status") == "error" &&
       std::string(determinant_rejected.at("detail").as_string())
               .find("exact structural identity") != std::string::npos &&
       laurent_accepted.at("status") == "ok" &&
       binding_rejected.at("status") == "error" &&
       std::string(binding_rejected.at("detail").as_string())
               .find("retained assembly operator") != std::string::npos &&
       identity_rejected.at("status") == "error" &&
       std::string(identity_rejected.at("detail").as_string())
               .find("nonidentity retained coefficient") !=
           std::string::npos);
  const bool ok = prepared.at("status") == "ok" &&
      solved.at("status") == "ok" && values_ok && rejection_ok &&
      checkpoint_ok &&
      stats.at("execution_implemented") == true &&
      stats.at("capability_evidence").as_object().at("identity_v") ==
          "native-retained-spectral-assembly-and-target-inverse";
  if (!ok) {
    std::cerr << domain << " prepared: " << json::serialize(prepared)
              << '\n' << domain << " solved: " << json::serialize(solved)
              << '\n' << domain << " evaluated: "
              << json::serialize(evaluated) << '\n'
              << domain << " determinant rejection: "
              << json::serialize(determinant_rejected) << '\n'
              << domain << " Laurent admission: "
              << json::serialize(laurent_accepted) << '\n'
              << domain << " V binding rejection: "
              << json::serialize(binding_rejected) << '\n'
              << domain << " identity fast-path rejection: "
              << json::serialize(identity_rejected) << '\n'
              << domain << " stats: " << json::serialize(stats) << '\n';
  }
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                              {"session", session}});
  return ok;
}

bool run_gauge_domain(const std::string& domain) {
  auto created = request(json::object{
      {"schema", 2}, {"op", "session.create"}, {"domain", domain},
      {"precision_bits", 256}, {"output_digits", 40},
      {"chart_capacity", 2}, {"scc_capacity", 1}, {"local_capacity", 1}});
  if (created.at("status") != "ok") return false;
  const auto session = std::string(created.at("session").as_string());
  const auto first = prepare_chart(session, domain, "gauge-source", "2",
                                   false, false);
  const auto second = prepare_chart(session, domain, "gauge-target", "3",
                                    false, false);
  const auto identity = domain + ":gauge-parent-v1";
  auto payload = manifest(domain, first, second, identity);
  auto& blocks = payload.at("blocks").as_array();
  blocks[0].as_object()["principal_identity"] =
      domain + ":gauge-source:principal";
  blocks[1].as_object()["principal_identity"] =
      domain + ":gauge-target:principal";
  blocks[0].as_object()["source_transform"] = source_transform(
      "gauge-source", "2", "1/2", domain);
  blocks[1].as_object()["source_transform"] = source_transform(
      "gauge-target", "3", "1/3", domain);
  for (auto& value : blocks) {
    value.as_object()["identity_gauge"] = false;
    value.as_object()["exact_gauge"] = true;
  }
  blocks[0].as_object()["to_physical"] = gauge_transform(
      "gauge-source", "to_physical", "7", "1/7", "7", domain);
  blocks[0].as_object()["to_reduced"] = gauge_transform(
      "gauge-source", "to_reduced", "7", "1/7", "1/7", domain);
  blocks[1].as_object()["to_physical"] = gauge_transform(
      "gauge-target", "to_physical", "11", "1/11", "11", domain);
  blocks[1].as_object()["to_reduced"] = gauge_transform(
      "gauge-target", "to_reduced", "11", "1/11", "1/11", domain);
  const auto prepared = prepare_scc(
      session, domain + ":gauge-key", std::move(payload));
  json::object solved{{"status", "not-run"}}, evaluated{{"status", "not-run"}};
  if (prepared.at("status") == "ok") {
    solved = request(json::object{
        {"schema", 2}, {"op", "scc.solve_column"}, {"session", session},
        {"scc", prepared.at("scc")}, {"checkpoint_identity", "gauge-column"},
        {"seed", json::object{{"block", 0}, {"run", regular_run(true)},
            {"metadata", metadata("gauge-seed")}}},
        {"targets", json::array{json::object{{"block", 1},
            {"run", regular_run(false)},
            {"metadata", metadata("gauge-target")}}}}});
    if (solved.at("status") == "ok")
      evaluated = request(json::object{
          {"schema", 2}, {"op", "local.evaluate"}, {"session", session},
          {"local", solved.at("local")},
          {"point", json::object{{"exact", "1/2"}}},
          {"options", json::object{{"tail_estimate", false}}}});
  }
  bool ok = false;
  if (evaluated.at("status") == "ok") {
    const auto& coefficients = evaluated.at("value").as_object()
                                   .at("coefficients").as_array();
    ok = coefficients.size() >= 2 &&
        std::abs(real_midpoint(coefficients[0]) - 14.0) < 1e-25 &&
        std::abs(real_midpoint(coefficients[1]) - 35.0) < 1e-25;
  }
  if (!ok)
    std::cerr << domain << " gauge prepared: " << json::serialize(prepared)
              << '\n' << domain << " gauge solved: " << json::serialize(solved)
              << '\n' << domain << " gauge evaluated: "
              << json::serialize(evaluated) << '\n';
  (void)request(json::object{{"schema", 2}, {"op", "session.close"},
                              {"session", session}});
  return ok;
}

}  // namespace

int main() {
  const bool rational = run_domain("rational", true);
  const bool acb = run_domain("acb", false);
  const bool dense_rational = run_dense_domain("rational");
  const bool dense_acb = run_dense_domain("acb");
  const bool gauge_rational = run_gauge_domain("rational");
  const bool gauge_acb = run_gauge_domain("acb");
  const bool ok = rational && acb && dense_rational && dense_acb &&
      gauge_rational && gauge_acb;
  std::cout << (ok ? "PASS" : "FAIL")
            << ": nonidentity Laurent-unimodular CompositeSCC spectral frames"
            << '\n';
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
