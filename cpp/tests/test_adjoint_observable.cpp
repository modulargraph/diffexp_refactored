#include "diffexp2/adjoint_observable.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace diffexp2;

namespace {

ComplexBall ball(const std::string& value) {
  return ComplexBall::from_strings(value);
}

EpsilonFrame<ComplexBall> series(
    std::initializer_list<const char*> coefficients,
    std::int32_t min_power = 0) {
  std::vector<ComplexBall> values;
  values.reserve(coefficients.size());
  for (const auto* coefficient : coefficients)
    values.push_back(ball(coefficient));
  return EpsilonFrame<ComplexBall>(min_power, std::move(values));
}

FiniteLaurentMatrix<ComplexBall> constant_matrix(
    const std::vector<std::vector<std::string>>& values,
    std::int32_t epsilon_complete_max = 0) {
  FiniteLaurentMatrix<ComplexBall> result;
  result.reserve(values.size());
  for (const auto& source_row : values) {
    FiniteLaurentVector<ComplexBall> row;
    row.reserve(source_row.size());
    for (const auto& value : source_row) {
      std::vector<ComplexBall> coefficients(
          static_cast<std::size_t>(epsilon_complete_max) + 1,
          ComplexBall(0));
      coefficients.front() = ball(value);
      row.emplace_back(0, std::move(coefficients));
    }
    result.push_back(std::move(row));
  }
  return result;
}

FiniteLaurentVector<ComplexBall> constant_vector(
    const std::vector<std::string>& values,
    std::int32_t epsilon_complete_max = 0) {
  FiniteLaurentVector<ComplexBall> result;
  result.reserve(values.size());
  for (const auto& value : values) {
    std::vector<ComplexBall> coefficients(
        static_cast<std::size_t>(epsilon_complete_max) + 1,
        ComplexBall(0));
    coefficients.front() = ball(value);
    result.emplace_back(0, std::move(coefficients));
  }
  return result;
}

bool contains_exact(const ComplexBall& value, const std::string& expected) {
  return (value - ball(expected)).contains_zero();
}

BackwardAdjointTaylorProblem scalar_problem(const std::string& alpha) {
  BackwardAdjointTaylorProblem problem;
  problem.dimension = 1;
  problem.taylor_complete_max = 4;
  problem.required_epsilon_complete_max = 0;
  problem.q_lags = {series({"1"})};
  problem.c_transpose_lags = {constant_matrix({{alpha}})};
  // q theta(lambda) + alpha lambda = -t.
  problem.forcing = {constant_vector({"-1"})};
  return problem;
}

bool scalar_regular_and_fractional() {
  const auto regular = solve_backward_adjoint_taylor(
      scalar_problem("0"), "scalar regular adjoint");
  const auto regular_value = evaluate_backward_adjoint_taylor(
      regular, ball("1/2"));
  if (!contains_exact(regular.coefficients.front().front().coefficient(0),
                      "-1") ||
      !contains_exact(regular_value.front().coefficient(0), "-1/2"))
    return false;

  // A physical mode f=t^(1/2) has an analytic composed adjoint even though
  // the separately constructed receiving basis is fractional-power.
  const auto fractional = solve_backward_adjoint_taylor(
      scalar_problem("1/2"), "scalar fractional-mode adjoint");
  const auto fractional_value = evaluate_backward_adjoint_taylor(
      fractional, ball("1/4"));
  return contains_exact(
             fractional.coefficients.front().front().coefficient(0),
             "-2/3") &&
         contains_exact(fractional_value.front().coefficient(0), "-1/6");
}

bool jordan_block_without_basis_logs() {
  BackwardAdjointTaylorProblem problem;
  problem.dimension = 2;
  problem.taylor_complete_max = 3;
  problem.required_epsilon_complete_max = 0;
  problem.q_lags = {series({"1"})};
  // Physical C=[[0,1],[0,0]], hence C^T below.  The physical fundamental
  // matrix has a logarithmic Jordan mode, while the composed adjoint for
  // this row remains the ordinary vector (-t,t).
  problem.c_transpose_lags = {
      constant_matrix({{"0", "0"}, {"1", "0"}})};
  problem.forcing = {constant_vector({"-1", "0"})};
  const auto solved = solve_backward_adjoint_taylor(
      problem, "Jordan composed adjoint");
  return contains_exact(solved.coefficients[0][0].coefficient(0), "-1") &&
         contains_exact(solved.coefficients[0][1].coefficient(0), "1");
}

bool epsilon_coupled_triangular() {
  BackwardAdjointTaylorProblem problem;
  problem.dimension = 2;
  problem.taylor_complete_max = 3;
  problem.required_epsilon_complete_max = 1;
  problem.q_lags = {series({"1", "0"})};
  auto c0 = constant_matrix({{"0", "0"}, {"0", "0"}}, 1);
  c0[1][0] = series({"1"}, 1);  // epsilon below the diagonal in C^T.
  problem.c_transpose_lags = {std::move(c0)};
  problem.forcing = {constant_vector({"-1", "0"}, 1)};
  const auto solved = solve_backward_adjoint_taylor(
      problem, "epsilon-coupled composed adjoint");
  const auto& first = solved.coefficients.front();
  return contains_exact(first[0].coefficient(0), "-1") &&
         contains_exact(first[1].coefficient(0), "0") &&
         contains_exact(first[1].coefficient(1), "1") &&
         solved.common_epsilon_complete_max >= 1;
}

bool direct_integral_equivalence() {
  // f=t satisfies theta(f)=f.  For the row r=1 and integration a -> 0,
  // lambda=-t/2 and lambda(a)f(a)=-a^2/2, exactly the direct integral.
  const auto solved = solve_backward_adjoint_taylor(
      scalar_problem("1"), "direct/composed scalar comparison");
  const auto adjoint = evaluate_backward_adjoint_taylor(
      solved, ball("1/4"));
  const auto physical = constant_vector({"1/4"});
  const auto composed = contract_backward_adjoint(adjoint, physical);
  return contains_exact(composed.coefficient(0), "-1/32");
}

bool resonance_fails_closed() {
  try {
    (void)solve_backward_adjoint_taylor(
        scalar_problem("-1"), "resonant composed adjoint");
  } catch (const MatchingArithmeticError& error) {
    return std::string(error.what()).find(
               "logarithmic/resonant Fuchsian completion is required") !=
           std::string::npos;
  }
  return false;
}

bool rational_shadow_pseudo_resonance() {
  BackwardAdjointTaylorProblem problem;
  problem.dimension = 1;
  problem.taylor_complete_max = 2;
  problem.required_epsilon_complete_max = 0;
  problem.q_lags = {series({"1", "0", "0"})};
  auto uncertain_minus_one = ball("-1");
  arb_add_error_2exp_si(acb_realref(uncertain_minus_one.raw()), -160);
  FiniteLaurentMatrix<ComplexBall> c0(1);
  c0.front().emplace_back(
      0, std::vector<ComplexBall>{
             std::move(uncertain_minus_one), ball("1"), ball("0")});
  problem.c_transpose_lags.push_back(std::move(c0));
  problem.forcing = {constant_vector({"-1"}, 2)};

  bool unguided_failed = false;
  try {
    (void)solve_backward_adjoint_taylor(
        problem, "unguided pseudo-resonant adjoint");
  } catch (const MatchingArithmeticError&) {
    unguided_failed = true;
  }
  if (!unguided_failed) return false;

  const auto exact_value = [](std::initializer_list<const char*> numerator) {
    ExactEpsilonRational<Rational> result;
    result.zero = false;
    for (const auto* coefficient : numerator)
      result.numerator.emplace_back(coefficient);
    result.denominator = {Rational(1)};
    return result;
  };
  PreparedPhysicalClearedODE<Rational> exact;
  exact.dimension = 1;
  exact.q_lags = {exact_value({"1"})};
  exact.c_lags = {{{0, 0, exact_value({"-1", "1"})}}};
  exact.owner_signature_identity = "pseudo-resonant-exact-owner";
  exact.payload_identity = "pseudo-resonant-exact-payload";
  exact.exact_payload_record = "pseudo-resonant-exact-record";

  const auto solved = solve_backward_adjoint_taylor(
      problem, &exact, "Rational-shadow pseudo-resonant adjoint");
  const auto& coefficient = solved.coefficients.front().front();
  return coefficient.min_power() == -1 &&
         contains_exact(coefficient.coefficient(-1), "-1") &&
         solved.common_epsilon_complete_max >= 0;
}

bool rational_shadow_true_scalar_resonance() {
  auto problem = scalar_problem("-1");
  problem.required_epsilon_complete_max = -1;
  problem.q_lags = {series({"1", "0", "0"})};
  problem.c_transpose_lags = {constant_matrix({{"-1"}}, 2)};
  problem.forcing = {constant_vector({"-1"}, 2)};
  const auto exact_value = [](const std::string& value) {
    ExactEpsilonRational<Rational> result;
    result.zero = value == "0";
    if (!result.zero) {
      result.numerator = {Rational(value)};
      result.denominator = {Rational(1)};
    }
    return result;
  };
  PreparedPhysicalClearedODE<Rational> exact;
  exact.dimension = 1;
  exact.q_lags = {exact_value("1")};
  exact.c_lags = {{{0, 0, exact_value("-1")}}};
  exact.owner_signature_identity = "true-resonant-exact-owner";
  exact.payload_identity = "true-resonant-exact-payload";
  exact.exact_payload_record = "true-resonant-exact-record";

  const auto solved = solve_backward_adjoint_taylor(
      problem, &exact, "Rational-shadow true scalar resonance");
  if (solved.max_log_power != 1 ||
      solved.higher_log_coefficients.front().size() != 1 ||
      !contains_exact(
          solved.coefficients.front().front().coefficient(0), "0") ||
      !contains_exact(
          solved.higher_log_coefficients.front().front().front()
              .coefficient(-1),
          "-1"))
    return false;
  const auto point = ball("1/2");
  const auto evaluated = evaluate_backward_adjoint_taylor(
      solved, point, local_detail::cb_log(point),
      "true scalar resonant evaluation");
  const auto expected = -point * local_detail::cb_log(point);
  return (evaluated.front().coefficient(0) - expected).contains_zero();
}

ExactEpsilonRational<ComplexBall> exact_constant(const std::string& value) {
  ExactEpsilonRational<ComplexBall> result;
  result.zero = value == "0";
  if (!result.zero) {
    result.numerator = {ball(value)};
    result.denominator = {ball("1")};
  }
  return result;
}

PreparedPhysicalClearedODE<ComplexBall> scalar_ode(
    const std::string& alpha) {
  PreparedPhysicalClearedODE<ComplexBall> ode;
  ode.dimension = 1;
  ode.q_lags = {exact_constant("1")};
  ode.c_lags = {{{0, 0, exact_constant(alpha)}}};
  if (alpha == "0") ode.c_lags.front().clear();
  ode.owner_signature_identity = "adjoint-test-owner";
  ode.payload_identity = "adjoint-test-payload";
  ode.exact_payload_record = "adjoint-test-record";
  return ode;
}

PreparedSparseLocalMultiplierMatrix<ComplexBall> scalar_row(
    std::uint32_t center_pole_order = 0, bool geometric = false,
    bool doubled_geometric = false) {
  PreparedSparseLocalMultiplierMatrix<ComplexBall> row;
  row.rows = 1;
  row.columns = 1;
  row.exact_identity = "adjoint-test-row";
  PreparedRationalTaylorMultiplier<ComplexBall> multiplier;
  multiplier.center_pole_order = center_pole_order;
  multiplier.kernels = {{ball("1"),
                         ball(geometric ?
                                  (doubled_geometric ? "2" : "1") : "0"),
                         ball(geometric ?
                                  (doubled_geometric ? "4" : "1") : "0"),
                         ball(geometric ?
                                  (doubled_geometric ? "8" : "1") : "0")}};
  PreparedRationalAnalyticCoefficient<ComplexBall> analytic;
  analytic.numerator = {ball("1")};
  analytic.denominator = geometric
      ? std::vector<ComplexBall>{
            ball("1"), ball(doubled_geometric ? "-2" : "-1")}
      : std::vector<ComplexBall>{ball("1")};
  multiplier.analytic_coefficients =
      std::vector<PreparedRationalAnalyticCoefficient<ComplexBall>>{
          std::move(analytic)};
  multiplier.exact_identity = "adjoint-test-multiplier";
  row.entries.push_back({0, 0, std::move(multiplier)});
  return row;
}

bool physical_ode_row_adapter() {
  const auto prepared = prepare_backward_adjoint_taylor_problem(
      scalar_ode("1"), scalar_row(), 4, 0, 0, ball("1"),
      "physical ODE/row composed adapter");
  const auto solved = solve_backward_adjoint_taylor(
      prepared, "physical ODE/row composed solve");
  if (!contains_exact(solved.coefficients[0][0].coefficient(0), "-1/2"))
    return false;
  try {
    (void)prepare_backward_adjoint_taylor_problem(
        scalar_ode("1"), scalar_row(1), 4, 0, 0, ball("1"),
        "center-pole composed adapter");
  } catch (const std::domain_error& error) {
    return std::string(error.what()).find(
               "Laurent/log backward Fuchsian completion is required") !=
           std::string::npos;
  }
  return false;
}

bool certified_geometric_tail() {
  const auto row = scalar_row(0, true);
  const auto problem = prepare_backward_adjoint_taylor_problem(
      scalar_ode("0"), row, 4, 0, 0, ball("1"),
      "geometric composed-tail adapter");
  const auto solved = solve_backward_adjoint_taylor(
      problem, "geometric composed-tail solve");
  const auto evaluation_point = ball("1/4");
  const auto witness_radius = ball("1/2");
  const auto forcing_bound =
      backward_adjoint_forcing_cauchy_numerator_upper(
          problem, row, ball("1"), witness_radius,
          "geometric composed forcing bound");
  const auto certificate = certify_backward_adjoint_taylor_tail(
      problem, solved, evaluation_point, witness_radius, forcing_bound,
      "geometric composed-tail certificate");
  const auto truncated = evaluate_backward_adjoint_taylor(
      solved, evaluation_point).front().coefficient(0);
  // lambda(t)=log(1-t) for -lambda'=1/(1-t), lambda(0)=0.
  const auto exact = local_detail::cb_log(ball("3/4"));
  return Magnitude::upper_abs(exact - truncated) <=
             certificate.absolute_vector_tail_upper &&
         certificate.recurrence_contraction_upper.is_zero() &&
         certificate.certified_after_taylor_order == 4;
}

bool tail_certificate_fails_closed() {
  const auto row = scalar_row(0, true);
  auto problem = prepare_backward_adjoint_taylor_problem(
      scalar_ode("0"), row, 4, 0, 0, ball("1"),
      "fail-closed composed-tail adapter");
  const auto solved = solve_backward_adjoint_taylor(
      problem, "fail-closed composed-tail solve");
  try {
    (void)backward_adjoint_forcing_cauchy_numerator_upper(
        problem, row, ball("1"), ball("2"),
        "crossed-denominator forcing bound");
    return false;
  } catch (const std::domain_error&) {
  }

  // A negative epsilon operator shift makes a finite causal stack depend on
  // unknown higher epsilon coefficients and therefore cannot be certified by
  // this theorem.
  problem.c_transpose_lags.front().front().front() = series({"1"}, -1);
  try {
    (void)certify_backward_adjoint_taylor_tail(
        problem, solved, ball("1/4"), ball("1/2"), Magnitude::one(),
        "negative-epsilon fail-closed tail");
  } catch (const std::domain_error& error) {
    return std::string(error.what()).find(
               "negative-epsilon C coupling") != std::string::npos;
  }
  return false;
}

bool adaptive_witness_respects_row_pole() {
  const auto row = scalar_row(0, true, true);
  const auto problem = prepare_backward_adjoint_taylor_problem(
      scalar_ode("0"), row, 4, 0, 0, ball("1"),
      "adaptive-witness composed-tail adapter");
  const auto solved = solve_backward_adjoint_taylor(
      problem, "adaptive-witness composed-tail solve");
  const auto certificate =
      certify_backward_adjoint_taylor_tail_adaptive_witness(
          problem, solved, row, ball("1"), ball("1/4"),
          Rational("1/4"), Rational(1), 8,
          "adaptive-witness row-pole regression");
  return Rational("1/4") < certificate.witness_radius_exact &&
         certificate.witness_radius_exact < Rational("1/2") &&
         certificate.tail.certified_after_taylor_order == 4;
}

}  // namespace

int main() {
  ComplexBall::set_precision(256);
  const bool ok = scalar_regular_and_fractional() &&
                  jordan_block_without_basis_logs() &&
                  epsilon_coupled_triangular() &&
                  direct_integral_equivalence() &&
                  resonance_fails_closed() &&
                  rational_shadow_pseudo_resonance() &&
                  rational_shadow_true_scalar_resonance() &&
                  physical_ode_row_adapter() &&
                  certified_geometric_tail() &&
                  adaptive_witness_respects_row_pole() &&
                  tail_certificate_fails_closed();
  std::cout << (ok ? "PASS" : "FAIL")
            << ": composed backward Fuchsian adjoint Taylor recurrence\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
