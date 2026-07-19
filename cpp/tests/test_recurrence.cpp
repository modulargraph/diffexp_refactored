#include "diffexp2/json_codec.hpp"
#include "diffexp2/recurrence.hpp"

#include <boost/json.hpp>

#include <cstdlib>
#include <iostream>
#include <string>

using diffexp2::BlockStep;
using diffexp2::JordanBlock;
using diffexp2::MatrixEntry;
using diffexp2::MatrixShift;
using diffexp2::PreparedLag;
using diffexp2::PreparedMatrix;
using diffexp2::Rational;
using diffexp2::RecurrenceProblem;
using diffexp2::RecurrenceSolver;
using diffexp2::ScalarShift;
using diffexp2::StepCase;
using diffexp2::SymbolicRational;

namespace {

int passed = 0;
int failed = 0;

void check(const std::string& label, bool condition) {
  if (condition) {
    ++passed;
    std::cout << "  PASS: " << label << '\n';
  } else {
    ++failed;
    std::cout << "  FAIL: " << label << '\n';
  }
}

RecurrenceProblem<Rational> exponential_problem(std::uint32_t nmax) {
  RecurrenceProblem<Rational> p;
  p.dimension = 1;
  p.nmax = nmax;
  p.log_max = 0;
  p.frame_base = -2;
  p.frame_width = 8;
  p.a_target = Rational(0);
  p.b_target = Rational(0);
  p.a_shift_min = 0;
  for (std::uint32_t n = 0; n <= nmax; ++n) p.a_shifts.emplace_back(n);
  p.d_lags = {{{0, Rational(1)}}};
  p.nhat_lags.resize(2);
  p.nhat_lags[0].valuations = {diffexp2::kCompleteInfinity};
  MatrixShift<Rational> identity;
  identity.shift = 0;
  identity.entries.push_back(MatrixEntry<Rational>{0, 0, Rational(1)});
  p.nhat_lags[1].polynomial.push_back(identity);
  p.nhat_lags[1].valuations = {0};
  p.d0_inverse_scalar = Rational(1);
  p.blocks = {JordanBlock{{0}}};
  p.schedule.resize(nmax + 1);
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    p.schedule[n] = {{n == 0 ? StepCase::Resonant : StepCase::Taylor,
                      Rational(n), Rational(0)}};
  }
  p.initial.assign(p.frame_width, Rational(0));
  p.initial[-p.frame_base] = Rational(1);
  p.initial_validity = {5};
  return p;
}

void test_exponential() {
  auto result = RecurrenceSolver<Rational>(exponential_problem(8)).run();
  Rational factorial(1);
  bool equal = true;
  for (std::uint32_t n = 0; n <= 8; ++n) {
    if (n > 0) factorial *= Rational(n);
    const auto index = (((static_cast<std::size_t>(n) * 2) * 1) * 8) + 2;
    equal = equal && result.u[index] == Rational(1) / factorial;
  }
  check("regular exponential exact coefficients", equal);
  check("regular exponential honest validity", result.top_valid == 5);
}

void test_epsilon_denominator() {
  // (1+eps) theta g = t g gives u_n = u_{n-1}/(n(1+eps)).
  auto p = exponential_problem(4);
  p.d_lags[0] = {{0, Rational(1)}, {1, Rational(1)}};
  p.d0_inverse_scalar.reset();
  auto result = RecurrenceSolver<Rational>(p).run();
  const auto at = [&](std::uint32_t n, std::uint32_t eps_index) -> const Rational& {
    return result.u[(((static_cast<std::size_t>(n) * 2) * 1) * 8) + eps_index];
  };
  check("epsilon denominator leading coefficient", at(1, 2) == Rational(1));
  check("epsilon denominator alternating coefficients",
        at(1, 3) == Rational(-1) && at(1, 4) == Rational(1) &&
        at(1, 5) == Rational(-1));
  check("epsilon denominator second Taylor coefficient",
        at(2, 2) == Rational("1/2") && at(2, 3) == Rational(-1));
}

void test_reusable_dense_factorization() {
  const std::vector<std::vector<Rational>> matrix{
      {Rational(0), Rational(2), Rational(1)},
      {Rational(1), Rational(1), Rational(0)},
      {Rational(2), Rational(0), Rational(1)}};
  const auto factorization =
      diffexp2::detail::factor_dense_system(matrix);
  const auto first = diffexp2::detail::solve_factored_dense_system(
      factorization, {Rational(7), Rational(3), Rational(5)});
  const auto second = diffexp2::detail::solve_factored_dense_system(
      factorization, {Rational(10), Rational(3), Rational(0)});
  check("one exact dense factorization solves multiple epsilon coefficients",
        first == std::vector<Rational>{
            Rational(1), Rational(2), Rational(3)} &&
            second == std::vector<Rational>{
                Rational(-1), Rational(4), Rational(2)});
}

void test_lower_frame_guard() {
  auto p = exponential_problem(1);
  p.nhat_lags[1].polynomial[0].shift = -3;
  bool loud = false;
  try {
    (void)RecurrenceSolver<Rational>(p).run();
  } catch (const diffexp2::RecurrenceError& error) {
    loud = error.id == "E4";
  }
  check("negative epsilon shift underflow is loud", loud);
}

RecurrenceProblem<Rational> resonant_jordan_source_problem(
    const Rational& top_first_component) {
  RecurrenceProblem<Rational> p;
  p.dimension = 2;
  p.nmax = 0;
  p.log_max = 1;
  p.frame_base = -1;
  p.frame_width = 4;
  p.has_initial = false;
  p.a_target = Rational(0);
  p.b_target = Rational(0);
  p.a_shift_min = 0;
  p.a_shifts = {Rational(0)};
  p.d_lags = {{{0, Rational(1)}}};
  p.nhat_lags.resize(1);
  p.nhat_lags[0].valuations.assign(4, diffexp2::kCompleteInfinity);
  p.d0_inverse_scalar = Rational(1);
  p.blocks = {JordanBlock{{0, 1}}};
  p.schedule = {{{StepCase::Resonant, Rational(0), Rational(0)}}};

  diffexp2::SourceData<Rational> source;
  source.present = {1, 1};
  source.validity.assign(4, 2);
  source.frames.assign(16, Rational(0));
  const auto at = [&](std::uint32_t log, std::uint32_t component,
                      std::int32_t power) -> Rational& {
    const auto point = static_cast<std::size_t>(log);
    const auto epsilon = static_cast<std::size_t>(power - p.frame_base);
    return source.frames[
        (point * p.dimension + component) * p.frame_width + epsilon];
  };
  // The last-row log-zero equation fixes y[1,1] = 2.  Consequently the
  // highest-log first-row equation requires Rtilde[1,0] = -2 exactly.
  at(0, 0, 0) = Rational(3);
  at(0, 1, 1) = Rational(2);
  at(1, 0, 0) = top_first_component;
  p.source = std::move(source);
  return p;
}

void test_resonant_jordan_log_ceiling_compatibility() {
  const auto compatible =
      RecurrenceSolver<Rational>(resonant_jordan_source_problem(Rational(-2)))
          .run();
  const auto coefficient = [&](std::uint32_t log, std::uint32_t component,
                               std::int32_t power) -> const Rational& {
    const auto epsilon = static_cast<std::size_t>(power + 1);
    return compatible.u[
        ((static_cast<std::size_t>(log) * 2 + component) * 4) + epsilon];
  };
  check("resonant Jordan ladder accepts an exactly compatible log ceiling",
        coefficient(1, 1, 0) == Rational(2) &&
            coefficient(0, 1, 0) == Rational(-3));

  bool loud = false;
  try {
    (void)RecurrenceSolver<Rational>(
        resonant_jordan_source_problem(Rational(-1))).run();
  } catch (const diffexp2::RecurrenceError& error) {
    loud = error.id == "E5" &&
        std::string(error.what()).find("captured log ceiling") !=
            std::string::npos;
  }
  check("resonant Jordan ladder rejects a hidden cross-log inconsistency",
        loud);
}

PreparedMatrix<Rational> identity_prepared_matrix(std::uint32_t dimension) {
  PreparedMatrix<Rational> matrix;
  matrix.identity = true;
  matrix.valuations.assign(
      static_cast<std::size_t>(dimension) * dimension,
      diffexp2::kCompleteInfinity);
  MatrixShift<Rational> identity;
  identity.shift = 0;
  for (std::uint32_t component = 0; component < dimension; ++component) {
    identity.entries.push_back(
        MatrixEntry<Rational>{component, component, Rational(1)});
    matrix.valuations[
        static_cast<std::size_t>(component) * dimension + component] = 0;
  }
  matrix.polynomial.push_back(std::move(identity));
  return matrix;
}

RecurrenceProblem<Rational> pseudo_jordan_problem(
    bool epsilon_regular_principal) {
  RecurrenceProblem<Rational> problem;
  problem.dimension = 2;
  problem.nmax = 1;
  problem.log_max = 0;
  problem.frame_base = -3;
  problem.frame_width = 7;
  problem.has_initial = true;
  problem.a_target = Rational(-1);
  problem.b_target = Rational(0);
  problem.a_shift_min = 0;
  problem.a_shifts = {Rational(-1), Rational(0)};
  problem.d_lags = {{{0, Rational(1)}}};
  problem.d0_inverse_scalar = Rational(1);
  problem.blocks = {JordanBlock{{0, 1}}};
  problem.schedule = {
      {{StepCase::Resonant, Rational(0), Rational(0)}},
      {{StepCase::Pseudo, Rational(0), Rational(1)}}};

  problem.nhat_lags.resize(2);
  problem.nhat_lags[0].valuations.assign(
      4, diffexp2::kCompleteInfinity);
  problem.nhat_lags[1] = identity_prepared_matrix(2);

  problem.initial.assign(
      static_cast<std::size_t>(problem.dimension) * problem.frame_width,
      Rational(0));
  problem.initial[
      static_cast<std::size_t>(-problem.frame_base)] = Rational(1);
  problem.initial_validity.assign(problem.dimension, 3);

  problem.epsilon_regular_principal = epsilon_regular_principal;
  if (epsilon_regular_principal) {
    PreparedLag<Rational> spectral_principal;
    spectral_principal.valuations.assign(
        4, diffexp2::kCompleteInfinity);
    problem.spectral_principal_lag = std::move(spectral_principal);
    problem.spectral_source_matrix = identity_prepared_matrix(2);
    problem.assembly_matrix = identity_prepared_matrix(2);
  }
  return problem;
}

void test_epsilon_regular_pseudo_jordan_transaction() {
  const auto ordinary =
      RecurrenceSolver<Rational>(pseudo_jordan_problem(false)).run();
  const auto epsilon_regular =
      RecurrenceSolver<Rational>(pseudo_jordan_problem(true)).run();
  const bool same_hits =
      ordinary.hits.size() == 1 && epsilon_regular.hits.size() == 1 &&
      ordinary.hits.front().n == epsilon_regular.hits.front().n &&
      ordinary.hits.front().columns ==
          epsilon_regular.hits.front().columns &&
      ordinary.hits.front().delta_b ==
          epsilon_regular.hits.front().delta_b &&
      ordinary.hits.front().gamma_frames ==
          epsilon_regular.hits.front().gamma_frames &&
      ordinary.hits.front().gamma_validity ==
          epsilon_regular.hits.front().gamma_validity;
  check("epsilon-regular exceptional layer preserves a size-two CASE-P Jordan transaction",
        ordinary.u == epsilon_regular.u &&
            ordinary.validity == epsilon_regular.validity && same_hits);
}

void test_json_error_contract() {
  const auto value = boost::json::parse(diffexp2::run_recurrence_json("{}"));
  check("malformed JSON request returns typed error",
        value.as_object().at("status") == "error" &&
        value.as_object().at("id") == "CPP");
}

void test_malformed_tensor_is_typed_error() {
  // This is otherwise a valid 1x1 request, but the assembly valuation tensor
  // is empty. It must be rejected before any unchecked native indexing; a
  // malformed public request must never be able to terminate WolframKernel.
  const std::string request = R"json({
    "schema":1,"domain":"rational","output_digits":30,
    "d":1,"nmax":0,"p":0,"fb":-1,"w":3,
    "has_initial":true,"adaptive_probe":false,
    "a_target":"0","b_target":"0","a_shift_min":0,
    "a_shifts":["0"],
    "d_lags":[[{"s":0,"v":"1"}]],
    "denominators":[],
    "nhat_lags":[{"poly":[],"rat":[],"val":[null]}],
    "d0_inverse":"1","blocks":[[0]],
    "schedule":[[{"case":"R","da":"0","db":"0"}]],
    "initial":["0","1","0"],"initial_validity":[1],
    "source":null,
    "assembly":{"identity":true,"poly":[],"rat":[],"val":[]},
    "chop_digits":10,"return_u":false
  })json";
  const auto value = boost::json::parse(diffexp2::run_recurrence_json(request));
  check("malformed native tensor returns typed error without crashing",
        value.as_object().at("status") == "error" &&
        value.as_object().at("id") == "E5");
}

void test_symbolic_rational_field() {
  SymbolicRational::configure({"rho"});
  const SymbolicRational rho("rho");
  const SymbolicRational rate("(1+rho)/(2-rho)");
  const auto identity = rate * SymbolicRational("2-rho") -
                        SymbolicRational("1+rho");
  check("symbolic rational field canonical cancellation", identity.is_zero());
  check("symbolic rational field retains regulator",
        rate.str().find("rho") != std::string::npos);
}

boost::json::object json_request(const std::string& request) {
  return boost::json::parse(diffexp2::run_recurrence_json(request)).as_object();
}

void test_persistent_operator_session() {
  const auto created = json_request(R"json({
    "schema":2,"op":"session.create","domain":"rational",
    "output_digits":30,"analytic":{"regulators":[],"branch":"euclidean"}
  })json");
  const auto session = std::string(created.at("session").as_string());
  const auto prepared = json_request(std::string(R"json({
    "schema":2,"op":"chart.prepare","session":")json") + session + R"json(",
    "key":"exp@0[-2,8]","identity":"exact-exp-v1",
    "analytic":{"prescription":"none"},
    "scc":{"components":[[0]],"structural_edges":[],
      "condensation_edges":[],"topological_order":[0],"coupling_depth":0},
    "problem":{"domain":"rational","d":1,"fb":-2,"w":8,
      "d_lags":[[{"s":0,"v":"1"}]],"denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],"val":[null]},
        {"poly":[{"s":0,"e":[[0,0,"1"]]}],"rat":[],"val":[0]}],
      "d0_inverse":"1","blocks":[[0]],"assembly":null,
      "chop_digits":10}
  })json");
  const auto chart = std::string(prepared.at("chart").as_string());
  const auto solved = json_request(std::string(R"json({
    "schema":2,"op":"chart.solve","session":")json") + session +
    R"json(","chart":")json" + chart + R"json(","run":{
      "nmax":2,"p":0,"has_initial":true,"adaptive_probe":false,
      "a_target":"0","b_target":"0","a_shift_min":0,
      "a_shifts":["0","1","2"],
      "schedule":[[{"case":"R","da":"0","db":"0"}],
        [{"case":"T","da":"1","db":"0"}],
        [{"case":"T","da":"2","db":"0"}]],
      "initial":["0","0","1","0","0","0","0","0"],
      "initial_validity":[5],"source":null,"return_u":true}
  })json");
  const auto batched = json_request(std::string(R"json({
    "schema":2,"op":"chart.solve_batch","session":")json") + session +
    R"json(","chart":")json" + chart + R"json(","threads":32,"runs":[{
      "nmax":0,"p":0,"has_initial":true,"adaptive_probe":false,
      "a_target":"0","b_target":"0","a_shift_min":0,
      "a_shifts":["0"],
      "schedule":[[{"case":"R","da":"0","db":"0"}]],
      "initial":["0","0","1","0","0","0","0","0"],
      "initial_validity":[5],"source":null,"return_u":true
    },{
      "nmax":"malformed"
    },{
      "nmax":0,"p":0,"has_initial":true,"adaptive_probe":false,
      "a_target":"0","b_target":"0","a_shift_min":0,
      "a_shifts":["0"],
      "schedule":[[{"case":"R","da":"0","db":"0"}]],
      "initial":["0","0","3","0","0","0","0","0"],
      "initial_validity":[5],"source":null,"return_u":true
    }]})json");
  const auto stats = json_request(std::string(R"json({
    "schema":2,"op":"session.stats","session":")json") + session + "\"}");
  check("persistent typed operator prepares and solves without static copies",
        prepared.at("status") == "ok" && solved.at("status") == "ok" &&
        prepared.at("scc_components") == 1 &&
        prepared.at("scc_structural_edges") == 0 &&
        prepared.at("scc_condensation_edges") == 0 &&
        prepared.at("scc_topological_order").as_array().size() == 1 &&
        prepared.at("scc_topological_order").as_array()[0] == 0 &&
        solved.at("persistent").as_object().at("static_tensor_copies") == 0 &&
        stats.at("runs") == 3 && stats.at("static_tensor_copies") == 0);
  const auto& batch_results = batched.at("results").as_array();
  check("persistent batch preserves order, typed errors, and run statistics",
        batched.at("status") == "ok" && batched.at("attempted") == 3 &&
        batched.at("succeeded") == 2 && batched.at("failed") == 1 &&
        batched.at("worker_threads") == 3 &&
        batch_results[0].as_object().at("status") == "ok" &&
        batch_results[0].as_object().at("u").as_array()[2] == "1" &&
        batch_results[1].as_object().at("status") == "error" &&
        batch_results[1].as_object().at("id") == "CPP" &&
        batch_results[2].as_object().at("status") == "ok" &&
        batch_results[2].as_object().at("u").as_array()[2] == "3");
  const auto& retained_scc =
      stats.at("chart_stats").as_array()[0].as_object();
  check("persistent stats retain typed SCC graph order",
        retained_scc.at("scc_components") == 1 &&
        retained_scc.at("scc_structural_edges") == 0 &&
        retained_scc.at("scc_condensation_edges") == 0 &&
        retained_scc.at("scc_topological_order").as_array()[0] == 0 &&
        retained_scc.at("scc_coupling_depth") == 0);
  (void)json_request(std::string(R"json({
    "schema":2,"op":"session.close","session":")json") + session + "\"}");

  const auto acb = json_request(R"json({
    "schema":2,"op":"session.create","domain":"acb",
    "precision_bits":128,"output_digits":30})json");
  const auto other_precision = json_request(R"json({
    "schema":2,"op":"session.create","domain":"acb",
    "precision_bits":256,"output_digits":30})json");
  check("persistent Acb sessions retain independent idle precisions",
        acb.at("status") == "ok" && other_precision.at("status") == "ok");
  (void)json_request(std::string(R"json({
    "schema":2,"op":"session.close","session":")json") +
    std::string(acb.at("session").as_string()) + "\"}");
  (void)json_request(std::string(R"json({
    "schema":2,"op":"session.close","session":")json") +
    std::string(other_precision.at("session").as_string()) + "\"}");
}

void test_persistent_framed_d0_inverse() {
  const auto created = json_request(R"json({
    "schema":2,"op":"session.create","domain":"rational",
    "output_digits":30,"analytic":{"regulators":[],"branch":"euclidean"}
  })json");
  const auto session = std::string(created.at("session").as_string());
  const auto prepared = json_request(std::string(R"json({
    "schema":2,"op":"chart.prepare","session":")json") + session + R"json(",
    "key":"exp-eps@0[-2,8]","identity":"exact-exp-eps-v1",
    "analytic":{"prescription":"none"},
    "scc":{"components":[[0]],"structural_edges":[],
      "condensation_edges":[],"topological_order":[0],"coupling_depth":0},
    "problem":{"domain":"rational","d":1,"fb":-2,"w":8,
      "d_lags":[[{"s":0,"v":"1"},{"s":1,"v":"1"}]],
      "denominators":[],
      "nhat_lags":[{"poly":[],"rat":[],"val":[null]},
        {"poly":[{"s":0,"e":[[0,0,"1"]]}],"rat":[],"val":[0]}],
      "d0_inverse":null,"blocks":[[0]],"assembly":null,
      "chop_digits":10}
  })json");
  const auto chart = std::string(prepared.at("chart").as_string());
  const auto run = std::string(R"json({"schema":2,"op":"chart.solve","session":")json") +
      session + R"json(","chart":")json" + chart + R"json(","run":{
        "nmax":1,"p":0,"has_initial":true,"adaptive_probe":false,
        "a_target":"0","b_target":"0","a_shift_min":0,
        "a_shifts":["0","1"],
        "schedule":[[{"case":"R","da":"0","db":"0"}],
          [{"case":"T","da":"1","db":"0"}]],
        "initial":["0","0","1","0","0","0","0","0"],
        "initial_validity":[5],"source":null,"return_u":true}})json";
  const auto first = json_request(run);
  const auto second = json_request(run);
  const auto stats = json_request(std::string(R"json({
    "schema":2,"op":"session.stats","session":")json") + session + "\"}");
  const auto& coefficients = first.at("u").as_array();
  const auto& chart_stats = stats.at("chart_stats").as_array().front().as_object();
  check("persistent chart retains one framed d0 inverse across runs",
        prepared.at("d0_inverse_mode") == "retained-frame" &&
            first.at("status") == "ok" && second.at("status") == "ok" &&
            coefficients[18] == "1" && coefficients[19] == "-1" &&
            chart_stats.at("runs") == 2 &&
            chart_stats.at("d0_inverse_mode") == "retained-frame");
  (void)json_request(std::string(R"json({
    "schema":2,"op":"session.close","session":")json") + session + "\"}");
}

}  // namespace

int main() {
  test_exponential();
  test_epsilon_denominator();
  test_reusable_dense_factorization();
  test_lower_frame_guard();
  test_resonant_jordan_log_ceiling_compatibility();
  test_epsilon_regular_pseudo_jordan_transaction();
  test_json_error_contract();
  test_malformed_tensor_is_typed_error();
  test_symbolic_rational_field();
  test_persistent_operator_session();
  test_persistent_framed_d0_inverse();
  std::cout << "Results: " << passed << " / " << (passed + failed)
            << " tests passed\n";
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
