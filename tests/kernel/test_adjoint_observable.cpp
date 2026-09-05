#include "diffexp/kernel/adjoint_observable.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace diffexp::kernel;

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
  if (!contains_exact(composed.coefficient(0), "-1/32")) return false;

  // The forcing convention depends on the requested a -> 0 orientation,
  // not on the sign of a or on whether a direct primitive implementation
  // sorts its endpoints first.
  const auto negative_adjoint = evaluate_backward_adjoint_taylor(
      solved, ball("-1/4"), ball("0"),
      "negative-point direct/composed scalar comparison");
  const auto negative_physical = constant_vector({"-1/4"});
  const auto negative_composed = contract_backward_adjoint(
      negative_adjoint, negative_physical);
  if (!contains_exact(negative_composed.coefficient(0), "-1/32"))
    return false;

  // With x=beta*t and beta=-2, b=-beta*t*q*r=2t.  The composed
  // contraction already includes the physical Jacobian and must equal
  // beta*Integral_a^0 t dt = 1/16 without another orientation factor.
  auto negative_scale_problem = scalar_problem("1");
  negative_scale_problem.forcing = {constant_vector({"2"})};
  const auto negative_scale = solve_backward_adjoint_taylor(
      negative_scale_problem, "negative-scale direct/composed comparison");
  const auto negative_scale_adjoint = evaluate_backward_adjoint_taylor(
      negative_scale, ball("1/4"));
  const auto negative_scale_composed = contract_backward_adjoint(
      negative_scale_adjoint, physical);
  return contains_exact(negative_scale_composed.coefficient(0), "1/16");
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

ExactEpsilonRational<Rational> exact_rational_constant(
    const std::string& value) {
  ExactEpsilonRational<Rational> result;
  result.zero = value == "0";
  if (!result.zero) {
    result.numerator = {Rational(value)};
    result.denominator = {Rational(1)};
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
  } catch (const BackwardAdjointCenterAnchoringError& error) {
    return error.forcing_power == 0 &&
           std::string(error.what()).find("endpoint pairing are required") !=
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
          problem, row, ball("1"), witness_radius, std::nullopt,
          nullptr, nullptr,
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
        problem, row, ball("1"), ball("2"), std::nullopt,
        nullptr, nullptr,
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

bool adaptive_witness_handles_narrow_denominator_clearance() {
  const Rational evaluation_exact("1/4");
  const Rational pole_exact =
      evaluation_exact + Rational(1) / Rational(1048576);
  const auto evaluation = ball(evaluation_exact.str());
  const auto inverse_pole = ball("1") / ball(pole_exact.str());
  auto row = scalar_row();
  auto& multiplier = row.entries.front().multiplier;
  auto power = ball("1");
  for (auto& coefficient : multiplier.kernels.front()) {
    coefficient = power;
    power *= inverse_pole;
  }
  multiplier.analytic_coefficients->front().denominator = {
      ball("1"), -inverse_pole};
  const auto problem = prepare_backward_adjoint_taylor_problem(
      scalar_ode("0"), row, 4, 0, 0, ball("1"),
      "narrow-clearance composed-tail adapter");
  const auto solved = solve_backward_adjoint_taylor(
      problem, "narrow-clearance composed-tail solve");
  bool shallow_failed = false;
  try {
    (void)certify_backward_adjoint_taylor_tail_adaptive_witness(
        problem, solved, row, ball("1"), evaluation,
        evaluation_exact, Rational(1), 16,
        "narrow-clearance shallow witness");
  } catch (const std::domain_error&) {
    shallow_failed = true;
  }
  if (!shallow_failed) return false;
  const auto deep =
      certify_backward_adjoint_taylor_tail_adaptive_witness(
          problem, solved, row, ball("1"), evaluation,
          evaluation_exact, Rational(1), 64,
          "narrow-clearance deep witness");
  return evaluation_exact < deep.witness_radius_exact &&
         deep.witness_radius_exact < pole_exact;
}

bool exact_combined_forcing_cancels_cleared_row_pole() {
  PreparedPhysicalClearedODE<ComplexBall> numeric_ode;
  numeric_ode.dimension = 1;
  numeric_ode.q_lags = {exact_constant("1"), exact_constant("-2")};
  numeric_ode.c_lags = {{}, {}};
  numeric_ode.owner_signature_identity = "cleared-pole-numeric-owner";
  numeric_ode.payload_identity = "cleared-pole-numeric-payload";
  numeric_ode.exact_payload_record = "cleared-pole-numeric-record";

  auto numeric_row = scalar_row();
  auto& numeric_multiplier = numeric_row.entries.front().multiplier;
  numeric_multiplier.kernels.front() = {
      ball("1"), ball("2"), ball("4"), ball("8")};
  numeric_multiplier.analytic_coefficients->front().denominator = {
      ball("1"), ball("-2")};
  const auto problem = prepare_backward_adjoint_taylor_problem(
      numeric_ode, numeric_row, 4, 0, 0, ball("1"),
      "cleared-row-pole composed adapter");
  const auto solved = solve_backward_adjoint_taylor(
      problem, "cleared-row-pole composed solve");

  bool separate_bound_failed = false;
  try {
    (void)certify_backward_adjoint_taylor_tail_adaptive_witness(
        problem, solved, numeric_row, ball("1"), ball("1/4"),
        Rational("1/4"), Rational(1), 1,
        "separate cleared-row-pole witness", 0);
  } catch (const std::domain_error&) {
    separate_bound_failed = true;
  }
  if (!separate_bound_failed) return false;

  PreparedPhysicalClearedODE<Rational> exact_ode;
  exact_ode.dimension = 1;
  exact_ode.q_lags = {
      exact_rational_constant("1"), exact_rational_constant("-2")};
  exact_ode.c_lags = {{}, {}};
  exact_ode.owner_signature_identity = "cleared-pole-exact-owner";
  exact_ode.payload_identity = "cleared-pole-exact-payload";
  exact_ode.exact_payload_record = "cleared-pole-exact-record";

  PreparedSparseLocalMultiplierMatrix<Rational> exact_row;
  exact_row.rows = 1;
  exact_row.columns = 1;
  exact_row.exact_identity = "cleared-pole-exact-row";
  PreparedRationalTaylorMultiplier<Rational> exact_multiplier;
  exact_multiplier.kernels = {{Rational(1), Rational(2), Rational(4),
                               Rational(8)}};
  exact_multiplier.analytic_coefficients =
      std::vector<PreparedRationalAnalyticCoefficient<Rational>>{{
          {Rational(1)}, {Rational(1), Rational(-2)}}};
  exact_multiplier.exact_identity = "cleared-pole-exact-multiplier";
  exact_row.entries.push_back({0, 0, std::move(exact_multiplier)});

  const auto crossed_combined_bound =
      exact_combined_backward_adjoint_forcing_disk_upper(
          problem, exact_row, exact_ode, ball("1"), ball("5/8"), 0,
          "crossed combined cleared-row-pole bound");
  const auto combined =
      certify_backward_adjoint_taylor_tail_adaptive_witness(
          problem, solved, numeric_row, ball("1"), ball("1/4"),
          Rational("1/4"), Rational(1), 2,
          "combined cleared-row-pole witness", 0,
          &exact_row, &exact_ode);
  return crossed_combined_bound.is_finite() &&
         combined.witness_radius_exact == Rational("7/16") &&
         combined.tail.absolute_vector_tail_upper.is_finite();
}

bool exact_q_normalization_removes_clearing_contraction() {
  PreparedPhysicalClearedODE<Rational> cleared;
  cleared.dimension = 1;
  cleared.q_lags = {
      exact_rational_constant("1"), exact_rational_constant("-2")};
  cleared.c_lags = {
      {{0, 0, exact_rational_constant("3")}},
      {{0, 0, exact_rational_constant("-6")}}};
  cleared.owner_signature_identity = "q-normalization-owner";
  cleared.payload_identity = "q-normalization-payload";
  cleared.exact_payload_record = "q-normalization-record";

  const auto normalized =
      adjoint_observable_detail::normalize_backward_adjoint_exact_ode_by_q(
          cleared, 4, 0, "q-normalization regression");
  if (normalized.truncated_ode.q_lags.size() != 1 ||
      normalized.truncated_ode.c_lags.size() != 5 ||
      normalized.truncated_ode.c_lags.front().size() != 1)
    return false;
  const auto normalized_c0 =
      adjoint_observable_detail::expand_exact_epsilon_rational(
          normalized.truncated_ode.c_lags.front().front().value, 0,
          "q-normalization C0");
  if (!(normalized_c0.coefficient(0) == Rational(3))) return false;
  for (std::size_t lag = 1;
       lag < normalized.truncated_ode.c_lags.size(); ++lag)
    if (!normalized.truncated_ode.c_lags[lag].empty()) return false;

  const auto row = scalar_row();
  const auto normalized_problem = prepare_backward_adjoint_taylor_problem(
      normalized.truncated_ode, row, 4, 0, 0, ball("1"),
      "q-normalized composed adapter");
  const auto normalized_solution = solve_backward_adjoint_taylor(
      normalized_problem, &normalized.truncated_ode,
      "q-normalized composed solve");
  if (!contains_exact(
          normalized_solution.coefficients.front().front().coefficient(0),
          "-1/4"))
    return false;
  for (std::size_t order = 1;
       order < normalized_solution.coefficients.size(); ++order)
    if (!contains_exact(
            normalized_solution.coefficients[order].front().coefficient(0),
            "0"))
      return false;

  const auto witness = ball("2/5");
  const auto forcing = backward_adjoint_forcing_cauchy_numerator_upper(
      normalized_problem, row, ball("1"), witness, 0,
      nullptr, nullptr, "q-normalized forcing bound");
  const auto c_tail =
      adjoint_observable_detail::
          normalized_backward_adjoint_c_tail_disk_upper(
              normalized, witness, 0, "q-normalized C-tail bound");
  const auto certificate = certify_backward_adjoint_taylor_tail(
      normalized_problem, normalized_solution, ball("1/4"), witness,
      forcing, "q-normalized tail certificate", c_tail);
  if (!certificate.recurrence_contraction_upper.is_zero()) return false;

  const auto cleared_problem = prepare_backward_adjoint_taylor_problem(
      cleared, row, 4, 0, 0, ball("1"),
      "cleared composed adapter control");
  const auto cleared_solution = solve_backward_adjoint_taylor(
      cleared_problem, &cleared, "cleared composed solve control");
  try {
    (void)certify_backward_adjoint_taylor_tail(
        cleared_problem, cleared_solution, ball("1/4"), witness,
        Magnitude::one(), "cleared tail control");
  } catch (const std::domain_error& error) {
    return std::string(error.what()).find("not contractive") !=
           std::string::npos;
  }
  return false;
}

bool nonunit_q_normalization_is_typed_applicability() {
  const auto make_ode = [] {
    PreparedPhysicalClearedODE<Rational> ode;
    ode.dimension = 1;
    ode.c_lags = {{{0, 0, exact_rational_constant("1")}}};
    ode.owner_signature_identity = "nonunit-q-owner";
    ode.payload_identity = "nonunit-q-payload";
    ode.exact_payload_record = "nonunit-q-record";
    return ode;
  };

  auto t_nonunit = make_ode();
  t_nonunit.q_lags = {
      exact_rational_constant("0"), exact_rational_constant("1")};
  try {
    (void)adjoint_observable_detail::
        normalize_backward_adjoint_exact_ode_by_q(
            t_nonunit, 4, 0, "t-nonunit q regression");
    return false;
  } catch (const BackwardAdjointCenterUnitError& error) {
    if (!error.t_valuation.has_value() || *error.t_valuation != 1)
      return false;
  }

  auto epsilon_nonunit = make_ode();
  ExactEpsilonRational<Rational> epsilon;
  epsilon.zero = false;
  epsilon.valuation = 1;
  epsilon.numerator = {Rational(1)};
  epsilon.denominator = {Rational(1)};
  epsilon_nonunit.q_lags = {std::move(epsilon)};
  try {
    (void)adjoint_observable_detail::
        normalize_backward_adjoint_exact_ode_by_q(
            epsilon_nonunit, 4, 1, "epsilon-nonunit q regression");
  } catch (const BackwardAdjointCenterUnitError& error) {
    return !error.t_valuation.has_value();
  }
  return false;
}

bool normalized_a_posteriori_defect_encloses_exact_tail() {
  PreparedPhysicalClearedODE<Rational> exact_ode;
  exact_ode.dimension = 1;
  exact_ode.q_lags = {exact_rational_constant("1")};
  exact_ode.c_lags = {{}};
  exact_ode.owner_signature_identity = "defect-tail-owner";
  exact_ode.payload_identity = "defect-tail-payload";
  exact_ode.exact_payload_record = "defect-tail-record";
  const auto normalized =
      adjoint_observable_detail::normalize_backward_adjoint_exact_ode_by_q(
          exact_ode, 4, 0, "a-posteriori defect normalization");

  PreparedSparseLocalMultiplierMatrix<Rational> exact_row;
  exact_row.rows = 1;
  exact_row.columns = 1;
  exact_row.exact_identity = "defect-tail-row";
  PreparedRationalTaylorMultiplier<Rational> multiplier;
  multiplier.kernels = {{Rational(1), Rational(1), Rational(1),
                         Rational(1), Rational(1)}};
  multiplier.analytic_coefficients =
      std::vector<PreparedRationalAnalyticCoefficient<Rational>>{{
          {Rational(1)}, {Rational(1), Rational(-1)}}};
  multiplier.exact_identity = "defect-tail-multiplier";
  exact_row.entries.push_back({0, 0, std::move(multiplier)});
  const auto row =
      adjoint_observable_detail::specialize_exact_backward_adjoint_row_taylor(
          exact_row, 4, "a-posteriori defect row");
  const auto problem = prepare_backward_adjoint_taylor_problem(
      normalized.truncated_ode, row, 4, 0, 0, ball("1"),
      "a-posteriori defect adapter");
  const auto solved = solve_backward_adjoint_taylor(
      problem, &normalized.truncated_ode,
      "a-posteriori defect solve");
  const auto certificate =
      certify_backward_adjoint_taylor_tail_adaptive_witness(
          problem, solved, row, ball("1"), ball("1/4"),
          Rational("1/4"), Rational("3/4"), 4,
          "a-posteriori defect certificate", 0,
          &exact_row, &normalized.truncated_ode, &normalized);
  const auto truncated = evaluate_backward_adjoint_taylor(
      solved, ball("1/4")).front().coefficient(0);
  const auto exact = local_detail::cb_log(ball("3/4"));
  return Magnitude::upper_abs(exact - truncated) <=
             certificate.tail.absolute_vector_tail_upper &&
         certificate.tail.recurrence_contraction_upper.is_zero() &&
         certificate.tail.absolute_vector_tail_upper.approximate_upper() <
             1e-3;
}

bool batched_real_interval_remainders_match_individual_quotients() {
  const adjoint_observable_detail::ExactTRationalFunction exact{
      {Rational(3), Rational(-5), Rational(7), Rational(2)},
      {Rational(2), Rational(-3), Rational(1)}};
  const auto interval = adjoint_observable_detail::exact_real_interval_ball(
      Rational("-1/10"), Rational("1/8"));
  const auto batched = adjoint_observable_detail::
      exact_t_rational_all_taylor_remainder_quotient_real_interval_uppers(
          exact, 12, interval, "batched real-interval remainder");
  if (batched.size() != 13) return false;
  for (std::uint32_t retained = 0; retained <= 12; ++retained) {
    const auto quotient = adjoint_observable_detail::
        exact_t_rational_taylor_remainder_quotient(
            exact, retained, "individual real-interval remainder");
    const auto individual = adjoint_observable_detail::
        exact_t_rational_real_interval_upper(
            quotient, interval, "individual real-interval upper");
    const auto left = batched[retained].approximate_upper();
    const auto right = individual.approximate_upper();
    if (std::abs(left - right) >
        1e-12 * std::max({1.0, std::abs(left), std::abs(right)}))
      return false;
  }
  return true;
}

bool real_ray_operator_bounds_are_cached_by_exact_interval() {
  PreparedPhysicalClearedODE<Rational> exact_ode;
  exact_ode.dimension = 1;
  exact_ode.q_lags = {exact_rational_constant("1")};
  exact_ode.c_lags = {
      {{0, 0, exact_rational_constant("3")}},
      {{0, 0, exact_rational_constant("1")}}};
  exact_ode.owner_signature_identity = "real-ray-cache-owner";
  exact_ode.payload_identity = "real-ray-cache-payload";
  exact_ode.exact_payload_record = "real-ray-cache-record";
  const auto normalized =
      adjoint_observable_detail::normalize_backward_adjoint_exact_ode_by_q(
          exact_ode, 6, 0, "real-ray cache normalization");
  const auto normalized_extended =
      adjoint_observable_detail::normalize_backward_adjoint_exact_ode_by_q(
          exact_ode, 6, 1, "real-ray cache extended normalization");
  adjoint_observable_detail::BackwardAdjointRealRayOperatorCache cache;
  const Rational left(0);
  const Rational right("1/4");
  const auto interval = adjoint_observable_detail::exact_real_interval_ball(
      left, right);
  const auto first = cache.get_or_build(
      normalized, left, right, interval, 6, 0,
      "real-ray cache first lookup");
  const auto second = cache.get_or_build(
      normalized_extended, left, right, interval, 6, 0,
      "real-ray cache prefix-compatible extended lookup");
  const auto first_stats = cache.stats();
  if (first.get() != second.get() || first_stats.entries != 1 ||
      first_stats.builds != 1 || first_stats.hits != 1)
    return false;

  const Rational split("1/8");
  const auto split_interval =
      adjoint_observable_detail::exact_real_interval_ball(left, split);
  (void)cache.get_or_build(
      normalized, left, split, split_interval, 6, 0,
      "real-ray cache split lookup");
  const auto split_stats = cache.stats();
  return split_stats.entries == 2 && split_stats.builds == 2 &&
         split_stats.hits == 1;
}

bool real_ray_tail_ignores_pole_beyond_endpoint() {
  PreparedPhysicalClearedODE<Rational> exact_ode;
  exact_ode.dimension = 1;
  exact_ode.q_lags = {exact_rational_constant("1")};
  exact_ode.c_lags = {{}};
  exact_ode.owner_signature_identity = "real-ray-owner";
  exact_ode.payload_identity = "real-ray-payload";
  exact_ode.exact_payload_record = "real-ray-record";
  const auto normalized =
      adjoint_observable_detail::normalize_backward_adjoint_exact_ode_by_q(
          exact_ode, 4, 0, "real-ray normalization");

  const Rational endpoint("1/4");
  const Rational pole = endpoint + Rational(1) / Rational(1048576);
  const Rational inverse_pole = Rational(1) / pole;
  PreparedSparseLocalMultiplierMatrix<Rational> exact_row;
  exact_row.rows = 1;
  exact_row.columns = 1;
  exact_row.exact_identity = "real-ray-row";
  PreparedRationalTaylorMultiplier<Rational> multiplier;
  multiplier.analytic_coefficients =
      std::vector<PreparedRationalAnalyticCoefficient<Rational>>{{
          {Rational(1)}, {Rational(1), -inverse_pole}}};
  multiplier.kernels = {{Rational(1)}};
  multiplier.exact_identity = "real-ray-multiplier";
  exact_row.entries.push_back({0, 0, std::move(multiplier)});
  const auto row =
      adjoint_observable_detail::specialize_exact_backward_adjoint_row_taylor(
          exact_row, 4, "real-ray exact row");
  const auto problem = prepare_backward_adjoint_taylor_problem(
      normalized.truncated_ode, row, 4, 0, 0, ball("1"),
      "real-ray adapter");
  const auto solved = solve_backward_adjoint_taylor(
      problem, &normalized.truncated_ode, "real-ray solve");
  const auto certificate = certify_backward_adjoint_real_ray_tail(
      normalized, exact_row, solved, ball("1"), endpoint, 0, 128,
      "real-ray beyond-endpoint-pole certificate");
  if (!certificate.absolute_vector_tail_upper.is_finite() ||
      !certificate.stability_margin_lower.is_finite() ||
      certificate.stability_margin_lower.is_zero() ||
      certificate.accepted_intervals == 0)
    return false;

  auto crossed_row = exact_row;
  const Rational crossed_pole = endpoint - Rational(1) / Rational(1048576);
  crossed_row.entries[0]
      .multiplier.analytic_coefficients->at(0).denominator =
      {Rational(1), -(Rational(1) / crossed_pole)};
  const auto crossed_taylor =
      adjoint_observable_detail::specialize_exact_backward_adjoint_row_taylor(
          crossed_row, 4, "real-ray crossed-pole row");
  const auto crossed_problem = prepare_backward_adjoint_taylor_problem(
      normalized.truncated_ode, crossed_taylor, 4, 0, 0, ball("1"),
      "real-ray crossed-pole adapter");
  const auto crossed_solution = solve_backward_adjoint_taylor(
      crossed_problem, &normalized.truncated_ode,
      "real-ray crossed-pole solve");
  try {
    (void)certify_backward_adjoint_real_ray_tail(
        normalized, crossed_row, crossed_solution, ball("1"), endpoint, 0,
        24, "real-ray crossed-pole certificate");
  } catch (const std::domain_error&) {
    return true;
  }
  return false;
}

bool coefficientwise_tail_ignores_irrelevant_high_epsilon_input() {
  FiniteLaurentVector<ComplexBall> adjoint{
      EpsilonFrame<ComplexBall>(0, {ball("1")})};
  std::vector<ComplexBall> incoming_coefficients(11, ball("0"));
  incoming_coefficients.front() = ball("2");
  incoming_coefficients.back() = ball("1e100");
  FiniteLaurentVector<ComplexBall> incoming{
      EpsilonFrame<ComplexBall>(0, std::move(incoming_coefficients))};
  const auto tails = backward_adjoint_contracted_tail_by_output(
      Magnitude::one(), adjoint, incoming, {0, 0},
      "coefficientwise irrelevant-input regression");
  return tails.size() == 1 &&
         tails.front().approximate_upper() < 3.0 &&
         tails.front().approximate_upper() >= 2.0;
}

bool forcing_bound_ignores_epsilon_above_private_cap() {
  auto row = scalar_row(0, true);
  auto& multiplier = row.entries.front().multiplier;
  multiplier.kernels.push_back(
      std::vector<ComplexBall>(4, ball("1e100")));
  PreparedRationalAnalyticCoefficient<ComplexBall> huge;
  huge.numerator = {ball("1e100")};
  huge.denominator = {ball("1"), ball("-1")};
  multiplier.analytic_coefficients->push_back(std::move(huge));
  const auto problem = prepare_backward_adjoint_taylor_problem(
      scalar_ode("0"), row, 4, 0, 0, ball("1"),
      "private-cap forcing-bound adapter");
  const auto capped = backward_adjoint_forcing_cauchy_numerator_upper(
      problem, row, ball("1"), ball("1/2"), 0, nullptr, nullptr,
      "private-cap forcing bound");
  const auto uncapped = backward_adjoint_forcing_cauchy_numerator_upper(
      problem, row, ball("1"), ball("1/2"), std::nullopt,
      nullptr, nullptr,
      "uncapped forcing bound");
  return capped.approximate_upper() < 2.0 &&
         uncapped.approximate_upper() > 1e90;
}

bool private_epsilon_reservoir_covers_inverse_loss() {
  const auto factory = [](std::int32_t input_complete_max) {
    const auto extended = [input_complete_max](
        std::int32_t minimum,
        std::initializer_list<std::pair<std::int32_t, const char*>> terms) {
      std::vector<ComplexBall> coefficients(
          static_cast<std::size_t>(input_complete_max - minimum) + 1,
          ball("0"));
      for (const auto& [power, value] : terms) {
        if (power > input_complete_max) continue;
        coefficients.at(static_cast<std::size_t>(power - minimum)) =
            ball(value);
      }
      return EpsilonFrame<ComplexBall>(minimum, std::move(coefficients));
    };
    BackwardAdjointTaylorProblem problem;
    problem.dimension = 2;
    problem.taylor_complete_max = 1;
    problem.required_epsilon_complete_max = 0;
    problem.q_lags = {extended(0, {{0, "1"}})};
    FiniteLaurentMatrix<ComplexBall> c(
        2, FiniteLaurentVector<ComplexBall>());
    c[0] = {extended(-1, {{-1, "1"}, {0, "-1"}}),
            extended(0, {})};
    c[1] = {extended(0, {}),
            extended(0, {{0, "-1"}, {1, "1"}})};
    problem.c_transpose_lags = {std::move(c)};
    std::vector<ComplexBall> rhs(
        static_cast<std::size_t>(input_complete_max) + 1,
        ball("0"));
    rhs.front() = ball("1");
    problem.forcing = {FiniteLaurentVector<ComplexBall>{
        EpsilonFrame<ComplexBall>(0, rhs),
        EpsilonFrame<ComplexBall>(0, std::move(rhs))}};
    return problem;
  };
  const auto solved =
      solve_backward_adjoint_taylor_with_epsilon_reservoir(
          factory, 0, 2, nullptr,
          "private epsilon inverse-loss regression");
  const auto larger_factory = [&](std::int32_t input_complete_max) {
    auto problem = factory(input_complete_max);
    problem.required_epsilon_complete_max = 1;
    return problem;
  };
  const auto larger =
      solve_backward_adjoint_taylor_with_epsilon_reservoir(
          larger_factory, 1, 3, nullptr,
          "larger private epsilon inverse-loss regression");
  if (solved.input_epsilon_complete_max != 2 ||
      solved.result.common_epsilon_complete_max < 0 ||
      larger.result.common_epsilon_complete_max < 1 ||
      backward_adjoint_prefix_input_complete_max(
          larger, 0, "private epsilon sliced-prefix regression") !=
          solved.input_epsilon_complete_max ||
      !contains_exact(
          solved.result.coefficients.front()[0].coefficient(1), "1") ||
      !contains_exact(
          solved.result.coefficients.front()[1].coefficient(-1), "1"))
    return false;
  for (std::size_t component = 0; component < 2; ++component) {
    const auto& smaller = solved.result.coefficients.front()[component];
    const auto& extended = larger.result.coefficients.front()[component];
    const auto common_min = std::max(
        smaller.min_power(), extended.min_power());
    for (std::int32_t power = common_min; power <= 0; ++power)
      if (!(smaller.coefficient(power) -
            extended.coefficient(power)).contains_zero())
        return false;
  }
  return true;
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
                  adaptive_witness_handles_narrow_denominator_clearance() &&
                  exact_combined_forcing_cancels_cleared_row_pole() &&
                  exact_q_normalization_removes_clearing_contraction() &&
                  nonunit_q_normalization_is_typed_applicability() &&
                  normalized_a_posteriori_defect_encloses_exact_tail() &&
                  batched_real_interval_remainders_match_individual_quotients() &&
                  real_ray_operator_bounds_are_cached_by_exact_interval() &&
                  real_ray_tail_ignores_pole_beyond_endpoint() &&
                  coefficientwise_tail_ignores_irrelevant_high_epsilon_input() &&
                  forcing_bound_ignores_epsilon_above_private_cap() &&
                  private_epsilon_reservoir_covers_inverse_loss() &&
                  tail_certificate_fails_closed();
  std::cout << (ok ? "PASS" : "FAIL")
            << ": composed backward Fuchsian adjoint Taylor recurrence\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
