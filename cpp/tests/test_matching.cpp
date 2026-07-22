#include "diffexp2/matching.hpp"

#include <chrono>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

using diffexp2::ComplexBall;
using diffexp2::EpsilonFrame;
using diffexp2::ExactLaurentMatrix;
using diffexp2::ExactLaurentPolynomial;
using diffexp2::FiniteLaurentMatrix;
using diffexp2::FiniteLaurentVector;
using diffexp2::AcbLaurentRefinementOptions;
using diffexp2::AcbMatchingResidualVerdict;
using diffexp2::MatchingArithmeticError;
using diffexp2::MatchingArithmeticErrorCode;
using diffexp2::Magnitude;
using diffexp2::Rational;

namespace {

int failed = 0;
int checked = 0;

void check(const std::string& label, bool condition) {
  ++checked;
  std::cout << "  " << (condition ? "PASS: " : "FAIL: ") << label << '\n';
  if (!condition) ++failed;
}

EpsilonFrame<Rational> frame(std::int32_t minimum,
                             std::initializer_list<long> values) {
  std::vector<Rational> coefficients;
  coefficients.reserve(values.size());
  for (const auto value : values) coefficients.emplace_back(value);
  return EpsilonFrame<Rational>(minimum, std::move(coefficients));
}

EpsilonFrame<Rational> rational_constant_frame(const std::string& value,
                                                std::size_t width = 6) {
  std::vector<Rational> coefficients;
  coefficients.reserve(width);
  coefficients.emplace_back(value);
  for (std::size_t i = 1; i < width; ++i) coefficients.emplace_back(0);
  return EpsilonFrame<Rational>(0, std::move(coefficients));
}

EpsilonFrame<ComplexBall> ball_constant_frame(const std::string& value,
                                               std::size_t width = 6) {
  std::vector<ComplexBall> coefficients;
  coefficients.reserve(width);
  coefficients.push_back(ComplexBall::from_strings(value));
  for (std::size_t i = 1; i < width; ++i) coefficients.emplace_back(0);
  return EpsilonFrame<ComplexBall>(0, std::move(coefficients));
}

EpsilonFrame<Rational> rational_epsilon_frame(const std::string& value,
                                               std::size_t width = 6) {
  std::vector<Rational> coefficients(width, Rational(0));
  coefficients[1] = Rational(value);
  return EpsilonFrame<Rational>(0, std::move(coefficients));
}

EpsilonFrame<ComplexBall> ball_epsilon_frame(const std::string& value,
                                              std::size_t width = 6) {
  std::vector<ComplexBall> coefficients(width, ComplexBall(0));
  coefficients[1] = ComplexBall::from_strings(value);
  return EpsilonFrame<ComplexBall>(0, std::move(coefficients));
}

EpsilonFrame<ComplexBall> ball_value_frame(const ComplexBall& value,
                                            std::size_t width = 6) {
  std::vector<ComplexBall> coefficients;
  coefficients.reserve(width);
  coefficients.push_back(value);
  for (std::size_t i = 1; i < width; ++i) coefficients.emplace_back(0);
  return EpsilonFrame<ComplexBall>(0, std::move(coefficients));
}

ComplexBall real_ball_with_error(long midpoint, slong error_exponent) {
  ComplexBall result(midpoint);
  arb_add_error_2exp_si(acb_realref(result.raw()), error_exponent);
  return result;
}

template <typename Scalar>
std::optional<std::vector<std::pair<std::size_t, std::size_t>>>
reference_exhaustive_full_rank_plan(
    diffexp2::matching_detail::DenseScalarMatrix<Scalar> matrix,
    std::size_t position = 0) {
  if (position == matrix.size())
    return std::vector<std::pair<std::size_t, std::size_t>>{};
  for (std::size_t row = position; row < matrix.size(); ++row) {
    for (std::size_t column = position; column < matrix.size(); ++column) {
      if (diffexp2::matching_detail::zero_decision(matrix[row][column]) !=
          diffexp2::matching_detail::ZeroDecision::Nonzero)
        continue;
      auto next = matrix;
      diffexp2::matching_detail::apply_certified_pivot(
          next, position, {row, column});
      auto tail = reference_exhaustive_full_rank_plan(
          std::move(next), position + 1);
      if (!tail.has_value()) continue;
      tail->insert(tail->begin(), {row, column});
      return tail;
    }
  }
  return std::nullopt;
}

void transformation_support_smoke() {
  const auto input = frame(0, {1, 2, 3});
  const auto constant = ExactLaurentPolynomial<Rational>::monomial(
      0, Rational(2));
  const auto pole = ExactLaurentPolynomial<Rational>::monomial(
      -1, Rational(1));
  const auto constant_result =
      diffexp2::exact_laurent_times_frame(constant, input);
  const auto pole_result = diffexp2::exact_laurent_times_frame(pole, input);
  check("exact constant does not truncate the finite frame",
        constant_result.min_power() == 0 &&
            constant_result.complete_max() == 2 &&
            constant_result.coefficient(2) == Rational(6));
  check("exact eps^-1 monomial shifts without truncating",
        pole_result.min_power() == -1 &&
            pole_result.complete_max() == 1 &&
            pole_result.coefficient(-1) == Rational(1) &&
            pole_result.coefficient(1) == Rational(3));

  ExactLaurentMatrix<Rational> transformation = {
      {pole, ExactLaurentPolynomial<Rational>::one()},
      {ExactLaurentPolynomial<Rational>::zero(),
       ExactLaurentPolynomial<Rational>::one()}};
  const auto transformed = diffexp2::apply_exact_laurent_matrix(
      transformation, {input, frame(0, {4, 5, 6})});
  check("formal Laurent matrix applies to a finite vector honestly",
        transformed[0].min_power() == -1 &&
            transformed[0].complete_max() == 1 &&
            transformed[0].coefficient(-1) == Rational(1) &&
            transformed[0].coefficient(0) == Rational(6) &&
            transformed[0].coefficient(1) == Rational(8) &&
            transformed[1].complete_max() == 2);

  const auto zero = frame(0, {0, 0, 0, 0, 0, 0});
  const auto one = frame(0, {1, 0, 0, 0, 0, 0});
  const auto one_plus_epsilon = frame(0, {1, 1, 0, 0, 0, 0});
  const auto epsilon = frame(0, {0, 1, 0, 0, 0, 0});
  const FiniteLaurentMatrix<Rational> finite_basis{
      {frame(0, {2, 3, 1, 0, 0, 0}), one},
      {epsilon, frame(0, {1, 4, 0, 0, 0, 0})}};
  const FiniteLaurentMatrix<Rational> right_frame{
      {one_plus_epsilon, epsilon}, {zero, one}};
  const ExactLaurentMatrix<Rational> exact_right{
      {ExactLaurentPolynomial<Rational>::one(),
       ExactLaurentPolynomial<Rational>::monomial(-1, Rational(1))},
      {ExactLaurentPolynomial<Rational>::zero(),
       ExactLaurentPolynomial<Rational>::one()}};
  const auto sequential = diffexp2::right_multiply_finite_by_exact_laurent(
      diffexp2::right_multiply_finite_laurent_matrices(
          finite_basis, right_frame),
      exact_right);
  const auto fused = diffexp2::right_multiply_finite_laurent_matrices(
      finite_basis,
      diffexp2::right_multiply_finite_by_exact_laurent(
          right_frame, exact_right));
  bool associative_prefix = true;
  for (std::size_t row = 0; row < 2; ++row)
    for (std::size_t column = 0; column < 2; ++column) {
      const auto minimum = std::max(
          sequential[row][column].min_power(),
          fused[row][column].min_power());
      const auto maximum = std::min(
          sequential[row][column].complete_max(),
          fused[row][column].complete_max());
      associative_prefix = associative_prefix && minimum <= maximum;
      for (std::int32_t power = minimum;
           power <= maximum; ++power)
        associative_prefix =
            associative_prefix &&
            sequential[row][column].coefficient(power) ==
                fused[row][column].coefficient(power);
    }
  check("finite right frame composes before an exact right map honestly",
        associative_prefix);
}

void quotient_and_solve_smoke() {
  const auto padded_incoming =
      diffexp2::matching_detail::
          trim_leading_exact_zeros_preserve_ambiguous(
              EpsilonFrame<ComplexBall>(
                  -3, std::vector<ComplexBall>{
                          ComplexBall(0), ComplexBall(0), ComplexBall(2),
                          ComplexBall(0)}));
  check("aligned incoming trims only its exact-zero padding",
        padded_incoming.min_power() == -1 &&
            padded_incoming.complete_max() == 0 &&
            padded_incoming.coefficient(-1).is_zero() == false);

  const auto ambiguous_leading = real_ball_with_error(0, -80);
  const auto ambiguous_incoming =
      diffexp2::matching_detail::
          trim_leading_exact_zeros_preserve_ambiguous(
              EpsilonFrame<ComplexBall>(
                  -3, std::vector<ComplexBall>{
                          ComplexBall(0), ambiguous_leading,
                          ComplexBall(2), ComplexBall(0)}));
  check("aligned incoming retains its first zero-ambiguous coefficient",
        ambiguous_incoming.min_power() == -2 &&
            !ambiguous_incoming.coefficient(-2).is_zero() &&
            ambiguous_incoming.coefficient(-2).contains_zero());

  const auto quotient = diffexp2::finite_laurent_quotient(
      frame(0, {2, 5, 8}), frame(-1, {0, 2, 1, 0}), "quotient smoke");
  check("finite quotient trims only certified denominator zeros",
        quotient.min_power() == 0 && quotient.complete_max() == 2 &&
            quotient.coefficient(0) == Rational(1) &&
            quotient.coefficient(1) == Rational(2) &&
            quotient.coefficient(2) == Rational(3));

  const auto inverse_epsilon = frame(-1, {1, 0, 0, 0, 0, 0});
  const auto one = frame(0, {1, 0, 0, 0, 0});
  const auto inverse_epsilon_plus_one = frame(-1, {1, 1, 0, 0, 0, 0});
  FiniteLaurentMatrix<Rational> matrix = {
      {inverse_epsilon, one}, {inverse_epsilon_plus_one, one}};
  const std::vector<EpsilonFrame<Rational>> rhs = {
      frame(0, {2, 0, 0, 0, 0}), frame(0, {2, 1, 0, 0, 0})};
  const auto weights = diffexp2::solve_finite_laurent_system(
      matrix, rhs, "two-by-two Laurent smoke");
  const auto reconstructed =
      diffexp2::apply_finite_laurent_matrix(matrix, weights);
  check("full-pivot finite Laurent solve preserves original column order",
        weights[0].min_power() <= 1 &&
            weights[0].complete_max() >= 1 &&
            weights[0].coefficient(1) == Rational(1) &&
            weights[1].min_power() <= 0 &&
            weights[1].complete_max() >= 0 &&
            weights[1].coefficient(0) == Rational(1) &&
            reconstructed[0].coefficient(0) == Rational(2) &&
            reconstructed[1].coefficient(1) == Rational(1));

  const auto wide_one = frame(0, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
  const auto wide_zero = frame(0, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
  const auto short_zero = frame(0, {0, 0, 0, 0, 0});
  const auto triangular = diffexp2::solve_finite_laurent_system<Rational>(
      {{wide_one, wide_zero}, {short_zero, wide_one}},
      {frame(0, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}),
       frame(0, {2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0})},
      "certified-zero triangular smoke");
  check("finite-zero below pivot shrinks RHS completeness honestly",
        triangular[0].complete_max() == 10 &&
            triangular[1].complete_max() == 4 &&
            triangular[0].coefficient(0) == Rational(1) &&
            triangular[1].coefficient(0) == Rational(2));

  const auto epsilon = frame(0, {0, 1, 0, 0, 0, 0});
  const auto epsilon_squared = frame(0, {0, 0, 1, 0, 0, 0});
  const auto epsilon_cubed = frame(0, {0, 0, 0, 1, 0, 0});
  const auto zero_wide = frame(0, {0, 0, 0, 0, 0, 0});
  const auto valuation_factorization =
      diffexp2::factor_finite_laurent_system<Rational>(
          {{epsilon, zero_wide}, {epsilon_squared, epsilon_cubed}},
          "certified below-pivot valuation window");
  const auto& valuation_factor =
      valuation_factorization.steps.front().eliminations.front().factor;
  check("certified below-pivot zeros do not consume quotient completeness",
        valuation_factor.min_power() == 1 &&
            valuation_factor.complete_max() == 4 &&
            valuation_factor.coefficient(1) == Rational(1));

  const auto twice_epsilon_cubed = frame(0, {0, 0, 0, 2, 0, 0});
  const auto schur_factorization =
      diffexp2::factor_finite_laurent_system<Rational>(
          {{epsilon, epsilon_squared},
           {epsilon_squared, twice_epsilon_cubed}},
          "certified pivot-row valuation window");
  check("certified pivot-row zeros do not consume Schur completeness",
        schur_factorization.upper[1][1].min_power() == 3 &&
            schur_factorization.upper[1][1].complete_max() == 5 &&
            schur_factorization.upper[1][1].coefficient(3) == Rational(1));

  const auto upper_window_matrix = FiniteLaurentMatrix<Rational>{
      {rational_constant_frame("1"), epsilon_squared},
      {zero_wide, rational_constant_frame("1")}};
  const auto upper_window_factorization =
      diffexp2::factor_exact_nonnegative_finite_laurent_system(
          upper_window_matrix,
          "certified stored-upper back-substitution window");
  const auto short_polar_weight = frame(-2, {1, 0, 0, 0});
  const auto wide_regular_weight = frame(0, {2, 0, 0, 0});
  const auto canonical_epsilon_squared = frame(2, {1, 0, 0, 0});
  const std::vector<EpsilonFrame<Rational>> upper_window_rhs = {
      wide_regular_weight + canonical_epsilon_squared * short_polar_weight,
      short_polar_weight};
  const auto upper_window_solution =
      diffexp2::solve_factorized_finite_laurent_system(
          upper_window_factorization, upper_window_rhs,
          "certified stored-upper back-substitution solve");
  check("certified stored-upper zeros preserve solved CompleteMax",
        upper_window_factorization.upper[0][1].min_power() == 2 &&
            upper_window_solution[0].complete_max() == 3 &&
            upper_window_solution[0].coefficient(0) == Rational(2) &&
            upper_window_solution[1].complete_max() == 1);

  const auto finite_zero_through_two =
      diffexp2::matching_detail::
          canonicalize_certified_or_preserve_ambiguous(
              frame(0, {0, 0, 0}), "finite-zero contract");
  const auto finite_zero_through_zero =
      diffexp2::matching_detail::
          canonicalize_certified_or_preserve_ambiguous(
              frame(0, {0}), "short finite-zero contract");
  const auto negative_power_zero_factor = diffexp2::finite_laurent_quotient(
      finite_zero_through_zero, epsilon,
      "finite zero divided by positive-valuation pivot");
  const auto finite_zero_schur_update =
      epsilon_squared -
      negative_power_zero_factor * finite_zero_through_two;
  check("finite zero times a negative-power factor loses completeness honestly",
        negative_power_zero_factor.complete_max() == -1 &&
            finite_zero_schur_update.complete_max() == 1);
}

void coefficientwise_power_series_solve_smoke() {
  constexpr std::size_t dimension = 8;
  constexpr std::size_t matrix_width = 7;
  FiniteLaurentMatrix<Rational> matrix(
      dimension, FiniteLaurentVector<Rational>());
  for (std::size_t row = 0; row < dimension; ++row) {
    matrix[row].reserve(dimension);
    for (std::size_t column = 0; column < dimension; ++column) {
      std::vector<Rational> coefficients(matrix_width, Rational(0));
      coefficients[0] = Rational(
          static_cast<long>(row == column ? 3 : row < column ? 1 : 0));
      if (row == column) coefficients[1] = Rational(1);
      matrix[row].emplace_back(0, std::move(coefficients));
    }
  }
  FiniteLaurentVector<Rational> expected;
  expected.reserve(dimension);
  for (std::size_t column = 0; column < dimension; ++column) {
    std::vector<Rational> coefficients(8, Rational(0));
    coefficients[0] = Rational(static_cast<long>(column + 1));
    coefficients[1] = Rational(static_cast<long>(2 * column + 1));
    coefficients[4] = Rational(static_cast<long>(column + 3));
    expected.emplace_back(-1, std::move(coefficients));
  }
  const auto rhs = diffexp2::apply_finite_laurent_matrix(matrix, expected);
  const auto factorization =
      diffexp2::factor_exact_nonnegative_power_series_system(
          matrix, "dimension-independent power-series solve");
  const auto solved =
      diffexp2::solve_factorized_exact_nonnegative_power_series_system(
          factorization, rhs,
          "dimension-independent coefficient recurrence");
  bool same_prefix = true;
  for (std::size_t column = 0; column < dimension; ++column) {
    same_prefix = same_prefix && solved[column].min_power() == -1 &&
        solved[column].complete_max() == 5;
    for (std::int32_t power = -1; power <= 5; ++power)
      same_prefix = same_prefix &&
          solved[column].coefficient(power) ==
              expected[column].coefficient(power);
  }
  const auto reconstructed =
      diffexp2::apply_finite_laurent_matrix(matrix, solved);
  bool exact_residual = true;
  for (std::size_t row = 0; row < dimension; ++row) {
    const auto residual = reconstructed[row] - rhs[row];
    exact_residual = exact_residual &&
        residual.complete_max() == 5 &&
        !diffexp2::finite_laurent_leading_power(
             residual, "coefficient recurrence residual").has_value();
  }
  check("coefficient recurrence preserves the honest top independently of dimension",
        same_prefix && exact_residual);

  const auto polar_rhs = FiniteLaurentVector<Rational>{
      EpsilonFrame<Rational>(-2, std::vector<Rational>(8, Rational(1)))};
  const auto scalar_factorization =
      diffexp2::factor_exact_nonnegative_power_series_system<Rational>(
          {{EpsilonFrame<Rational>(
              0, std::vector<Rational>(6, Rational(1)))}},
          "negative-floor power-series bound");
  const auto polar_solution =
      diffexp2::solve_factorized_exact_nonnegative_power_series_system(
          scalar_factorization, polar_rhs,
          "negative-floor coefficient recurrence");
  check("coefficient recurrence loses only the true unknown matrix tail",
        polar_solution.front().min_power() == -2 &&
            polar_solution.front().complete_max() == 3);

  // A transformed CASE-P basis can carry an Acb enclosure around a negative
  // coefficient which the exact saturation witness proves is identically
  // zero.  Both the recurrence and its residual certificate must consume the
  // same witness-projected matrix; replaying the unprojected numerical
  // remnant would test a different linear system.
  const FiniteLaurentMatrix<Rational> numerical_remnant = {
      {frame(-1, {7, 1, 0, 0})}};
  const auto witness_projected =
      diffexp2::exact_nonnegative_finite_laurent_matrix(
          numerical_remnant, "witness-projected residual matrix");
  const FiniteLaurentVector<Rational> witness_rhs = {
      frame(-2, {3, 5, 0, 0})};
  const auto witness_factorization =
      diffexp2::factor_exact_nonnegative_power_series_system(
          witness_projected, "witness-projected residual factorization");
  const auto witness_solution =
      diffexp2::solve_factorized_exact_nonnegative_power_series_system(
          witness_factorization, witness_rhs,
          "witness-projected residual solve");
  const auto projected_reconstruction =
      diffexp2::apply_finite_laurent_matrix(
          witness_projected, witness_solution);
  const auto unprojected_reconstruction =
      diffexp2::apply_finite_laurent_matrix(
          numerical_remnant, witness_solution);
  const auto projected_residual =
      projected_reconstruction.front() - witness_rhs.front();
  const auto unprojected_residual =
      unprojected_reconstruction.front() - witness_rhs.front();
  check("exact-witness residual reuses the projected nonnegative matrix",
        !diffexp2::finite_laurent_leading_power(
             projected_residual,
             "witness-projected coefficient residual").has_value() &&
            diffexp2::finite_laurent_leading_power(
                unprojected_residual,
                "unprojected numerical-remnant residual").has_value());
}

void epsilon_lattice_saturation_smoke() {
  const auto one = frame(0, {1, 0, 0});
  const auto zero = frame(0, {0, 0, 0});
  const auto epsilon = frame(0, {0, 1, 0});
  const FiniteLaurentMatrix<Rational> basis = {
      {one, one}, {zero, epsilon}};
  const auto saturated = diffexp2::saturate_finite_laurent_basis(
      basis, "one-action saturation witness");

  check("epsilon-lattice witness has one certified saturation action",
        saturated.diagnostics.normalized_determinant_valuation == 1 &&
            saturated.diagnostics.initial_leading_rank == 1 &&
            saturated.diagnostics.final_leading_rank == 2 &&
            saturated.diagnostics.actions.size() == 1 &&
            saturated.diagnostics.actions.front().target_column == 1);
  check("saturation action records the deterministic null relation",
        saturated.diagnostics.actions.front().null_relation[0] ==
                Rational(-1) &&
            saturated.diagnostics.actions.front().null_relation[1] ==
                Rational(1));

  const auto& product = saturated.basis_times_transformation;
  check("saturation makes the epsilon-zero leading frame invertible",
        product[0][0].coefficient(0) == Rational(1) &&
            product[0][1].coefficient(0) == Rational(0) &&
            product[1][0].coefficient(0) == Rational(0) &&
            product[1][1].coefficient(0) == Rational(1));

  const auto& transformation = saturated.transformation;
  check("exact-support T retains both shift and column-combination semantics",
        transformation[0][0].coefficient(0) == Rational(1) &&
            transformation[0][1].coefficient(-1) == Rational(-1) &&
            transformation[1][0].is_zero() &&
            transformation[1][1].coefficient(-1) == Rational(1));

  const auto independently_multiplied =
      diffexp2::right_multiply_finite_by_exact_laurent(basis,
                                                        transformation);
  check("returned transformed basis agrees coefficientwise with F*T",
        independently_multiplied[0][0].coefficient(0) ==
                product[0][0].coefficient(0) &&
            independently_multiplied[0][1].coefficient(0) ==
                product[0][1].coefficient(0) &&
            independently_multiplied[1][0].coefficient(0) ==
                product[1][0].coefficient(0) &&
            independently_multiplied[1][1].coefficient(0) ==
                product[1][1].coefficient(0));

  const auto shifted = diffexp2::saturate_finite_laurent_basis<Rational>(
      {{frame(2, {1, 0, 0})}}, "initial shift witness");
  check("initial column valuation is accumulated as the same monomial in T",
        shifted.diagnostics.initial_column_valuations ==
                std::vector<std::int32_t>{2} &&
            shifted.diagnostics.initial_column_shifts ==
                std::vector<std::int32_t>{-2} &&
            shifted.diagnostics.actions.empty() &&
            shifted.transformation[0][0].coefficient(-2) == Rational(1) &&
            shifted.basis_times_transformation[0][0].coefficient(0) ==
                Rational(1));

  // The determinant of this confluent basis is 2 eps^3, although every
  // entry is known only through eps^2.  Sequential exact saturation has
  // enough local information to prove three divisions; expanding the
  // original determinant in the common entrywise rectangle does not.
  const FiniteLaurentMatrix<Rational> confluent = {
      {frame(0, {1, 0, 0}), frame(0, {1, 0, 0}),
       frame(0, {1, 0, 0})},
      {frame(0, {0, 0, 0}), frame(0, {0, 1, 0}),
       frame(0, {0, 2, 0})},
      {frame(0, {0, 0, 0}), frame(0, {0, 0, 1}),
       frame(0, {0, 0, 4})}};
  const auto confluent_saturated =
      diffexp2::saturate_finite_laurent_basis(
          confluent, "confluent action-derived determinant witness");
  check("exact saturation certifies determinant valuation beyond the common entry top",
        confluent_saturated.diagnostics.actions.size() == 3 &&
            confluent_saturated.diagnostics
                    .normalized_determinant_valuation == 3 &&
            confluent_saturated.diagnostics.normalized_determinant
                    .min_power() == 3 &&
            confluent_saturated.diagnostics.normalized_determinant
                    .complete_max() == 3 &&
            confluent_saturated.diagnostics.normalized_determinant
                    .coefficient(3) == Rational(2) &&
            confluent_saturated.diagnostics.final_leading_rank == 3);
}

void acb_saturation_candidate_chop_smoke() {
  ComplexBall::set_precision(1024);
  const auto tenth = ComplexBall::from_strings("0.1");
  const auto three_tenths = ComplexBall(3) * tenth;
  std::vector<ComplexBall> epsilon_perturbed{
      three_tenths, ComplexBall(1), ComplexBall(0)};
  const FiniteLaurentMatrix<ComplexBall> basis = {
      {ball_constant_frame("1", 3), ball_value_frame(tenth, 3)},
      {ball_constant_frame("3", 3),
       EpsilonFrame<ComplexBall>(0, std::move(epsilon_perturbed))}};

  bool strict_rejected = false;
  try {
    (void)diffexp2::saturate_finite_laurent_basis(
        basis, "strict Acb dependency-loss saturation");
  } catch (const MatchingArithmeticError& error) {
    strict_rejected =
        error.code == MatchingArithmeticErrorCode::AmbiguousZero;
  }

  const auto saturated = diffexp2::saturate_finite_laurent_basis(
      basis, "candidate-chopped Acb dependency-loss saturation", 100);
  const auto& product = saturated.basis_times_transformation;
  const bool normalized =
      saturated.diagnostics.normalized_determinant_valuation == 1 &&
      saturated.diagnostics.initial_leading_rank == 1 &&
      saturated.diagnostics.final_leading_rank == 2 &&
      saturated.diagnostics.actions.size() == 1 &&
      product[0][0].coefficient(0).contains_zero() == false &&
      product[0][1].coefficient(0).is_zero() &&
      product[1][0].coefficient(0).contains_zero() == false &&
      (product[1][1].coefficient(0) - ComplexBall(1)).contains_zero();
  check("strict Acb saturation rejects dependency-loss zero overlaps",
        strict_rejected);
  check("certified candidate chop recovers the nontrivial Acb lattice",
        normalized);
  ComplexBall::set_precision(256);
}

void ambiguous_acb_pivot_smoke() {
  ComplexBall::set_precision(256);
  ComplexBall ambiguous;
  arb_add_error_2exp_si(acb_realref(ambiguous.raw()), 0);
  bool rejected = false;
  try {
    (void)diffexp2::finite_laurent_quotient(
        EpsilonFrame<ComplexBall>(0, {ComplexBall(1)}),
        EpsilonFrame<ComplexBall>(0, {ambiguous}), "Acb ambiguity smoke");
  } catch (const MatchingArithmeticError& error) {
    rejected = error.code == MatchingArithmeticErrorCode::AmbiguousZero &&
               error.epsilon_power == 0;
  }
  check("Acb zero-enclosing Laurent divisor is loud", rejected);

  bool saturation_rejected = false;
  try {
    (void)diffexp2::saturate_finite_laurent_basis<ComplexBall>(
        {{EpsilonFrame<ComplexBall>(
            0, {ambiguous, ComplexBall(1), ComplexBall(0)})}},
        "Acb saturation ambiguity smoke");
  } catch (const MatchingArithmeticError& error) {
    saturation_rejected =
        error.code == MatchingArithmeticErrorCode::AmbiguousZero &&
        error.epsilon_power == 0;
  }
  check("Acb zero-overlap in saturation valuation is loud",
        saturation_rejected);
}

void refined_acb_match_smoke() {
  ComplexBall::set_precision(256);
  const FiniteLaurentMatrix<Rational> exact_basis = {
      {rational_constant_frame("1"), rational_constant_frame("1")},
      {rational_constant_frame("0"), rational_epsilon_frame("1")}};
  const auto exact_record = diffexp2::saturate_finite_laurent_basis(
      exact_basis, "exact nontrivial-T Acb match record");
  AcbLaurentRefinementOptions options;
  options.relative_tolerance = Magnitude::decimal("1e-50");
  options.required_complete_max = 3;
  options.max_refinement_steps = 2;
  const auto matched = diffexp2::refine_acb_finite_laurent_match(
      {{ball_constant_frame("1"), ball_constant_frame("1")},
       {ball_constant_frame("0"), ball_epsilon_frame("1")}},
      {ball_constant_frame("3"), ball_epsilon_frame("2")}, exact_record,
      options, "nontrivial-T refined Acb match");
  check("exact saturation T drives a certified Acb match",
        exact_record.diagnostics.actions.size() == 1 &&
        matched.refinement_steps == 0 &&
            matched.residual_history.size() == 1 &&
            matched.residual_history.back().verdict ==
                AcbMatchingResidualVerdict::Pass &&
            matched.residual_history.back().complete_through_required &&
            matched.weights[0].complete_max() >= 3 &&
            matched.weights[1].complete_max() >= 3 &&
            (matched.weights[0].coefficient(0) - ComplexBall(1)).is_zero() &&
            (matched.weights[1].coefficient(0) - ComplexBall(2)).is_zero());
}

void refined_acb_match_without_public_upper_slack_smoke() {
  ComplexBall::set_precision(256);
  constexpr std::size_t width = 4;
  const FiniteLaurentMatrix<Rational> exact_basis = {
      {rational_constant_frame("1", width),
       rational_constant_frame("1", width)},
      {rational_constant_frame("0", width),
       rational_epsilon_frame("1", width)}};
  const auto exact_record = diffexp2::saturate_finite_laurent_basis(
      exact_basis, "exact no-upper-slack Acb match record");
  AcbLaurentRefinementOptions options;
  options.relative_tolerance = Magnitude::decimal("1e-50");
  options.required_complete_max = 3;
  options.max_refinement_steps = 2;
  const auto matched = diffexp2::refine_acb_finite_laurent_match(
      {{ball_constant_frame("1", width),
        ball_constant_frame("1", width)},
       {ball_constant_frame("0", width),
        ball_epsilon_frame("1", width)}},
      {ball_constant_frame("3", width),
       ball_epsilon_frame("2", width)},
      exact_record, options, "no-upper-slack refined Acb match");
  check("private candidate zeros preserve a certified public match edge",
        matched.residual_history.size() == 1 &&
            matched.residual_history.back().verdict ==
                AcbMatchingResidualVerdict::Pass &&
            matched.residual_history.back().complete_through_required &&
            matched.weights[0].complete_max() == 3 &&
            matched.weights[1].complete_max() == 3 &&
            (matched.weights[0].coefficient(0) - ComplexBall(1)).is_zero() &&
            (matched.weights[1].coefficient(0) - ComplexBall(2)).is_zero());
}

void exact_right_frame_residual_smoke() {
  ComplexBall::set_precision(256);
  constexpr std::size_t width = 6;
  const FiniteLaurentMatrix<Rational> exact_basis = {
      {rational_constant_frame("1", width),
       rational_constant_frame("1", width)},
      {rational_constant_frame("0", width),
       rational_epsilon_frame("1", width)}};
  const auto exact_record = diffexp2::saturate_finite_laurent_basis(
      exact_basis, "exact right-frame residual record");
  const FiniteLaurentMatrix<ComplexBall> physical_basis = {
      {ball_constant_frame("1", width),
       ball_constant_frame("1", width)},
      {ball_constant_frame("0", width),
       ball_epsilon_frame("1", width)}};
  std::vector<ComplexBall> second_rhs(width, ComplexBall(0));
  second_rhs[1] = ComplexBall::from_strings("1/3");
  const FiniteLaurentVector<ComplexBall> rhs = {
      ball_constant_frame("1", width),
      EpsilonFrame<ComplexBall>(0, std::move(second_rhs))};
  AcbLaurentRefinementOptions options;
  options.relative_tolerance = Magnitude::zero();
  options.required_min_power = 0;
  options.required_complete_max = 0;
  options.max_refinement_steps = 0;

  const auto matched = diffexp2::refine_acb_finite_laurent_match(
      physical_basis, rhs, exact_record, options,
      "exact right-frame physical residual");
  const auto direct =
      diffexp2::matching_detail::evaluate_acb_matching_residual(
          physical_basis, matched.weights, rhs, options,
          "dependency-lost physical right-frame residual");
  check("exact right lattice certifies a physical residual without replaying confluent columns",
        exact_record.diagnostics.actions.size() == 1 &&
            direct.diagnostics.verdict ==
                AcbMatchingResidualVerdict::Inconclusive &&
            matched.residual_history.size() == 1 &&
            matched.residual_history.back().verdict ==
                AcbMatchingResidualVerdict::Pass &&
            matched.residual_history.back().complete_through_required);
}

void retained_transformed_basis_authority_smoke() {
  ComplexBall::set_precision(256);
  constexpr std::size_t width = 6;
  const FiniteLaurentMatrix<ComplexBall> unrelated_base_basis = {
      {ball_constant_frame("0", width),
       ball_constant_frame("0", width)},
      {ball_constant_frame("0", width),
       ball_constant_frame("0", width)}};
  const FiniteLaurentMatrix<ComplexBall> retained_proposal_basis = {
      {ball_constant_frame("1", width),
       ball_constant_frame("0", width)},
      {ball_constant_frame("0", width),
       ball_constant_frame("1", width)}};
  const FiniteLaurentMatrix<ComplexBall> authoritative_physical_basis = {
      {ball_constant_frame("2", width),
       ball_constant_frame("0", width)},
      {ball_constant_frame("0", width),
       ball_constant_frame("3", width)}};
  const ExactLaurentMatrix<ComplexBall> coordinate_transformation = {
      {ExactLaurentPolynomial<ComplexBall>::one(),
       ExactLaurentPolynomial<ComplexBall>::zero()},
      {ExactLaurentPolynomial<ComplexBall>::zero(),
       ExactLaurentPolynomial<ComplexBall>::one()}};
  const FiniteLaurentVector<ComplexBall> physical_rhs = {
      ball_constant_frame("4", width),
      ball_constant_frame("9", width)};
  const FiniteLaurentVector<ComplexBall> deliberately_rough_proposal_rhs = {
      ball_constant_frame("1", width),
      ball_constant_frame("1", width)};
  const std::function<FiniteLaurentVector<ComplexBall>(
      const FiniteLaurentVector<ComplexBall>&)>
      physical_residual_to_proposal_rhs =
          [](const FiniteLaurentVector<ComplexBall>& residual) {
            FiniteLaurentVector<ComplexBall> normalized;
            normalized.reserve(2);
            normalized.push_back(
                residual[0].scaled(
                    ComplexBall::from_strings("1/2")));
            normalized.push_back(
                residual[1].scaled(
                    ComplexBall::from_strings("1/3")));
            return normalized;
          };
  AcbLaurentRefinementOptions options;
  options.relative_tolerance = Magnitude::decimal("1e-50");
  options.required_min_power = 0;
  options.required_complete_max = 3;
  options.max_refinement_steps = 1;

  const auto matched = diffexp2::refine_acb_finite_laurent_match(
      unrelated_base_basis, deliberately_rough_proposal_rhs,
      coordinate_transformation,
      options, "retained transformed physical basis authority",
      false, true, &retained_proposal_basis,
      &authoritative_physical_basis, &physical_rhs,
      &physical_residual_to_proposal_rhs, true, 512);
  check("normalized proposal corrections can be certified in authoritative physical coordinates",
        matched.refinement_steps == 1 &&
            matched.residual_history.size() == 2 &&
            matched.residual_history.back().verdict ==
                AcbMatchingResidualVerdict::Pass &&
            matched.residual_history.back().complete_through_required &&
            (matched.transformed_weights[0].coefficient(0) -
             ComplexBall(2)).contains_zero() &&
            (matched.transformed_weights[1].coefficient(0) -
             ComplexBall(3)).contains_zero());

  bool missing_authority_rejected = false;
  try {
    (void)diffexp2::refine_acb_finite_laurent_match(
        retained_proposal_basis, physical_rhs,
        coordinate_transformation, options,
        "missing retained transformed basis authority",
        false, true, &retained_proposal_basis,
        &authoritative_physical_basis, nullptr,
        &physical_residual_to_proposal_rhs, true, 512);
  } catch (const MatchingArithmeticError& error) {
    missing_authority_rejected =
        error.code ==
        MatchingArithmeticErrorCode::InvalidSaturationLattice;
  }
  check("transformed-basis authority cannot be enabled without its retained matrix",
        missing_authority_rejected &&
            ComplexBall::precision() == 256);

  constexpr const char* denominator =
      "1532495540865888858358347027150309183618739122183602176";
  const std::string one_plus_delta =
      std::string(
          "1532495540865888858358347027150309183618739122183602177/") +
      denominator;
  const std::string two_plus_delta =
      std::string(
          "3064991081731777716716694054300618367237478244367204353/") +
      denominator;
  const FiniteLaurentMatrix<ComplexBall> ill_conditioned_basis = {
      {ball_constant_frame("1", width),
       ball_constant_frame("1", width)},
      {ball_constant_frame("1", width),
       ball_constant_frame(one_plus_delta, width)}};
  auto uncertain_two = ComplexBall(2);
  auto uncertain_two_plus_delta =
      ComplexBall::from_strings(two_plus_delta);
  arb_add_error_2exp_si(
      acb_realref(uncertain_two.raw()), -150);
  arb_add_error_2exp_si(
      acb_realref(uncertain_two_plus_delta.raw()), -150);
  const FiniteLaurentVector<ComplexBall> uncertain_rhs = {
      ball_value_frame(uncertain_two, width),
      ball_value_frame(uncertain_two_plus_delta, width)};
  AcbLaurentRefinementOptions ill_conditioned_options;
  ill_conditioned_options.relative_tolerance =
      Magnitude::decimal("1e-20");
  ill_conditioned_options.required_min_power = 0;
  ill_conditioned_options.required_complete_max = 3;
  ill_conditioned_options.max_refinement_steps = 0;
  bool enclosing_inverse_not_certified = false;
  try {
    const auto enclosing =
        diffexp2::refine_acb_finite_laurent_match(
            ill_conditioned_basis, uncertain_rhs,
            coordinate_transformation, ill_conditioned_options,
            "ill-conditioned enclosing inverse proposal");
    enclosing_inverse_not_certified =
        enclosing.residual_history.back().verdict !=
        AcbMatchingResidualVerdict::Pass;
  } catch (const MatchingArithmeticError&) {
    enclosing_inverse_not_certified = true;
  }
  const std::function<FiniteLaurentVector<ComplexBall>(
      const FiniteLaurentVector<ComplexBall>&)>
      identity_correction =
          [](const FiniteLaurentVector<ComplexBall>& residual) {
            return residual;
          };
  const auto midpoint = diffexp2::refine_acb_finite_laurent_match(
      ill_conditioned_basis, uncertain_rhs,
      coordinate_transformation, ill_conditioned_options,
      "ill-conditioned midpoint proposal with full-ball certificate",
      false, true, &ill_conditioned_basis,
      &ill_conditioned_basis, &uncertain_rhs,
      &identity_correction, true, 512);
  check("midpoint proposal avoids inverse-radius blow-up while the full-ball forward residual remains authoritative",
        enclosing_inverse_not_certified &&
            midpoint.residual_history.back().verdict ==
                AcbMatchingResidualVerdict::Pass &&
            (midpoint.transformed_weights[0].coefficient(0) -
             ComplexBall(1)).contains_zero() &&
            (midpoint.transformed_weights[1].coefficient(0) -
             ComplexBall(1)).contains_zero() &&
            ComplexBall::precision() == 256);

  ComplexBall::set_precision(1024);
  ComplexBall tiny_delta;
  arf_set_ui_2exp_si(
      arb_midref(acb_realref(tiny_delta.raw())), 1, -600);
  const auto one_plus_tiny_delta =
      ComplexBall(1) + tiny_delta;
  const auto two_plus_tiny_delta =
      ComplexBall(2) + tiny_delta;
  const FiniteLaurentMatrix<ComplexBall>
      retained_high_precision_basis = {
          {ball_constant_frame("1", width),
           ball_constant_frame("1", width)},
          {ball_constant_frame("1", width),
           ball_value_frame(one_plus_tiny_delta, width)}};
  auto high_precision_rhs_first = ComplexBall(2);
  auto high_precision_rhs_second = two_plus_tiny_delta;
  arb_add_error_2exp_si(
      acb_realref(high_precision_rhs_first.raw()), -700);
  arb_add_error_2exp_si(
      acb_realref(high_precision_rhs_second.raw()), -700);
  const FiniteLaurentVector<ComplexBall>
      retained_high_precision_rhs = {
          ball_value_frame(high_precision_rhs_first, width),
          ball_value_frame(high_precision_rhs_second, width)};
  ComplexBall::set_precision(256);
  bool base_precision_midpoint_not_certified = false;
  try {
    const auto base_precision_midpoint =
        diffexp2::refine_acb_finite_laurent_match(
            retained_high_precision_basis,
            retained_high_precision_rhs,
            coordinate_transformation, ill_conditioned_options,
            "retained high-precision basis at base proposal precision",
            false, true, &retained_high_precision_basis,
            &retained_high_precision_basis,
            &retained_high_precision_rhs,
            &identity_correction, true);
    base_precision_midpoint_not_certified =
        base_precision_midpoint.residual_history.back().verdict !=
        AcbMatchingResidualVerdict::Pass;
  } catch (const MatchingArithmeticError&) {
    base_precision_midpoint_not_certified = true;
  }
  const auto elevated_precision_midpoint =
      diffexp2::refine_acb_finite_laurent_match(
          retained_high_precision_basis,
          retained_high_precision_rhs,
          coordinate_transformation, ill_conditioned_options,
          "retained high-precision basis at elevated proposal precision",
          false, true, &retained_high_precision_basis,
          &retained_high_precision_basis,
          &retained_high_precision_rhs,
          &identity_correction, true, 1024);
  check("scoped proposal precision recovers retained input bits and restores the session precision",
        base_precision_midpoint_not_certified &&
            elevated_precision_midpoint.residual_history.back().verdict ==
                AcbMatchingResidualVerdict::Pass &&
            ComplexBall::precision() == 256);
}

void certified_tiny_physical_weight_publication_smoke() {
  ComplexBall::set_precision(1024);
  constexpr std::size_t width = 3;
  const FiniteLaurentMatrix<ComplexBall> identity_basis = {
      {ball_constant_frame("1", width)}};
  const FiniteLaurentVector<ComplexBall> tiny_rhs = {
      ball_constant_frame("1e-150", width)};
  const ExactLaurentMatrix<ComplexBall> identity_transformation = {
      {ExactLaurentPolynomial<ComplexBall>::one()}};
  const std::function<FiniteLaurentVector<ComplexBall>(
      const FiniteLaurentVector<ComplexBall>&)>
      identity_correction =
          [](const FiniteLaurentVector<ComplexBall>& residual) {
            return residual;
          };
  AcbLaurentRefinementOptions options;
  options.relative_tolerance = Magnitude::decimal("1e-180");
  options.required_min_power = 0;
  options.required_complete_max = 0;
  options.max_refinement_steps = 0;

  const auto matched = diffexp2::refine_acb_finite_laurent_match(
      identity_basis, tiny_rhs, identity_transformation, options,
      "certified tiny physical publication weight", false, true,
      &identity_basis, &identity_basis, &tiny_rhs, &identity_correction,
      true);
  const auto expected = ComplexBall::from_strings("1e-150");
  check("a certified tiny endpoint-sensitive weight is published unchopped",
        matched.residual_history.back().verdict ==
                AcbMatchingResidualVerdict::Pass &&
            !matched.weights[0].coefficient(0).is_zero() &&
            (matched.weights[0].coefficient(0) - expected).contains_zero() &&
            !matched.transformed_weights[0].coefficient(0).is_zero());
  ComplexBall::set_precision(256);
}

void structural_transformed_weight_floor_smoke() {
  ComplexBall::set_precision(256);
  constexpr std::size_t width = 4;
  std::vector<ComplexBall> identity_coefficients(width, ComplexBall(0));
  identity_coefficients[1] = ComplexBall(1);
  const auto identity_frame = EpsilonFrame<ComplexBall>(
      -1, std::move(identity_coefficients));
  std::vector<ComplexBall> rhs_coefficients(width, ComplexBall(0));
  rhs_coefficients[0] = ComplexBall::from_strings("1e-14");
  rhs_coefficients[1] = ComplexBall(1);
  const FiniteLaurentMatrix<ComplexBall> identity_basis = {
      {identity_frame}};
  const FiniteLaurentVector<ComplexBall> rhs = {
      EpsilonFrame<ComplexBall>(-1, std::move(rhs_coefficients))};
  const ExactLaurentMatrix<ComplexBall> identity_transformation = {
      {ExactLaurentPolynomial<ComplexBall>::one()}};
  AcbLaurentRefinementOptions options;
  options.relative_tolerance = Magnitude::decimal("1e-50");
  options.required_min_power = -1;
  options.required_complete_max = 0;
  options.max_refinement_steps = 0;

  const auto unconstrained = diffexp2::refine_acb_finite_laurent_match(
      identity_basis, rhs, identity_transformation, options,
      "unconstrained transverse Taylor defect", false, true,
      &identity_basis);
  const auto floored = diffexp2::refine_acb_finite_laurent_match(
      identity_basis, rhs, identity_transformation, options,
      "structurally floored transverse Taylor defect", false, true,
      &identity_basis, nullptr, nullptr, nullptr, false, 0,
      std::int32_t{0});
  check("an exact transformed-weight floor cannot be fitted away by a negative Laurent coefficient",
        unconstrained.residual_history.back().verdict ==
                AcbMatchingResidualVerdict::Pass &&
            unconstrained.transformed_weights[0].min_power() == -1 &&
            floored.residual_history.back().verdict !=
                AcbMatchingResidualVerdict::Pass &&
            floored.transformed_weights[0].min_power() == 0 &&
            (floored.transformed_weights[0].coefficient(0) -
             ComplexBall(1)).contains_zero());
}

void ill_scaled_refinement_smoke() {
  ComplexBall::set_precision(256);
  const std::string big =
      "1000000000000000000000000000000000000000000000000000000000000";
  const std::string twice_over_big = "2/" + big;
  const FiniteLaurentMatrix<Rational> exact_basis = {
      {rational_constant_frame(big), rational_constant_frame("1")},
      {rational_constant_frame("1"),
       rational_constant_frame(twice_over_big)}};
  const auto exact_record = diffexp2::saturate_finite_laurent_basis(
      exact_basis, "ill-scaled exact saturation record");

  AcbLaurentRefinementOptions options;
  // A zero tolerance deliberately cannot accept a nonzero-radius enclosure.
  // This exercises bounded correction proposals and must remain honest by
  // retaining the best Inconclusive prefix rather than manufacturing a
  // numerical zero or publishing a non-improving iterate.
  options.relative_tolerance = Magnitude::zero();
  options.required_complete_max = 4;
  options.max_refinement_steps = 2;
  const auto matched = diffexp2::refine_acb_finite_laurent_match(
      {{ball_constant_frame(big), ball_constant_frame("1")},
       {ball_constant_frame("1"), ball_constant_frame(twice_over_big)}},
      {ball_constant_frame("0"), ball_constant_frame("-1")}, exact_record,
      options, "ill-scaled refined Acb match");
  const auto expected_large = ComplexBall::from_strings("-" + big);
  check("ill-scaled Acb solve retains the best bounded refinement honestly",
        matched.refinement_steps == 0 &&
            matched.residual_history.size() == 1 &&
            matched.residual_history.back().verdict ==
                AcbMatchingResidualVerdict::Inconclusive &&
            matched.residual_history.back().complete_through_required &&
            (matched.weights[0].coefficient(0) - ComplexBall(1))
                .contains_zero() &&
            (matched.weights[1].coefficient(0) - expected_large)
                .contains_zero());
}

void incomplete_refinement_rollback_smoke() {
  ComplexBall::set_precision(256);
  std::vector<Rational> exact_basis_coefficients(19, Rational(0));
  exact_basis_coefficients[3] = Rational(3);
  std::vector<ComplexBall> ball_basis_coefficients(19, ComplexBall(0));
  ball_basis_coefficients[3] = ComplexBall(3);
  std::vector<ComplexBall> rhs_coefficients(17, ComplexBall(0));
  rhs_coefficients[1] = ComplexBall(1);

  const auto exact_record = diffexp2::saturate_finite_laurent_basis(
      FiniteLaurentMatrix<Rational>{{EpsilonFrame<Rational>(
          -3, std::move(exact_basis_coefficients))}},
      "rollback exact saturation record");
  AcbLaurentRefinementOptions options;
  // Division by three creates a nonzero-radius zero enclosure.  With a zero
  // tolerance the first residual is deliberately inconclusive but complete
  // through epsilon^12; feeding that residual back consumes three upper
  // orders and would previously replace it with an incomplete epsilon^9
  // iterate, causing an endless +3 reservoir retry.
  options.relative_tolerance = Magnitude::zero();
  options.required_min_power = -1;
  options.required_complete_max = 12;
  options.max_refinement_steps = 2;
  const auto matched = diffexp2::refine_acb_finite_laurent_match(
      {{EpsilonFrame<ComplexBall>(
          -3, std::move(ball_basis_coefficients))}},
      {EpsilonFrame<ComplexBall>(-1, std::move(rhs_coefficients))},
      exact_record, options, "incomplete correction rollback", true);
  check("incomplete Acb correction retains the last complete prefix",
        matched.refinement_steps == 0 &&
            matched.residual_history.size() == 1 &&
            matched.residual_history.back().complete_through_required &&
            matched.residual_history.back().complete_window.complete_max >=
                12 &&
            matched.residual_history.back().verdict ==
                AcbMatchingResidualVerdict::Inconclusive);
}

void refined_acb_ambiguous_pivot_smoke() {
  ComplexBall::set_precision(256);
  const auto exact_record = diffexp2::saturate_finite_laurent_basis<Rational>(
      {{rational_constant_frame("1")}},
      "exact ambiguous-pivot control record");
  ComplexBall ambiguous;
  arb_add_error_2exp_si(acb_realref(ambiguous.raw()), 0);
  std::vector<ComplexBall> coefficients;
  coefficients.push_back(std::move(ambiguous));
  for (std::size_t i = 1; i < 6; ++i) coefficients.emplace_back(0);
  bool rejected = false;
  try {
    AcbLaurentRefinementOptions options;
    options.required_complete_max = 4;
    (void)diffexp2::refine_acb_finite_laurent_match(
        {{EpsilonFrame<ComplexBall>(0, std::move(coefficients))}},
        {ball_constant_frame("1")}, exact_record, options,
        "ambiguous refined Acb pivot");
  } catch (const MatchingArithmeticError& error) {
    rejected = error.code == MatchingArithmeticErrorCode::AmbiguousZero &&
               error.epsilon_power == 0;
  }
  check("exact T never licenses an Acb pivot whose enclosure overlaps zero",
        rejected);
}

void laurent_off_pivot_ambiguity_smoke() {
  ComplexBall::set_precision(256);
  ComplexBall ambiguous;
  arb_add_error_2exp_si(acb_realref(ambiguous.raw()), -80);
  const auto one = ball_constant_frame("1");
  const auto zero = ball_constant_frame("0");
  const auto epsilon = ball_epsilon_frame("1");
  const auto overlap = ball_value_frame(ambiguous);
  const FiniteLaurentMatrix<ComplexBall> matrix = {
      {one, overlap}, {zero, epsilon}};
  const std::vector<EpsilonFrame<ComplexBall>> expected = {
      ball_constant_frame("2"), ball_constant_frame("3")};
  const auto right_hand_side =
      diffexp2::apply_finite_laurent_matrix(matrix, expected);
  const auto solved = diffexp2::solve_finite_laurent_system(
      matrix, right_hand_side, "ambiguous off-pivot Laurent solve");
  check("Laurent factorization ignores an ambiguous off-pivot valuation",
        (solved[0].coefficient(0) - ComplexBall(2)).contains_zero() &&
            (solved[1].coefficient(0) - ComplexBall(3)).contains_zero());

  bool rejected = false;
  try {
    (void)diffexp2::factor_finite_laurent_system<ComplexBall>(
        {{overlap}}, "ambiguous only Laurent pivot");
  } catch (const MatchingArithmeticError& error) {
    rejected = error.code == MatchingArithmeticErrorCode::AmbiguousZero &&
               error.row == 0 && error.column == 0 &&
               error.epsilon_power == 0;
  }
  check("Laurent factorization rejects an ambiguous required pivot loudly",
        rejected);
}

void refined_acb_ambiguous_off_pivot_smoke() {
  ComplexBall::set_precision(256);
  ComplexBall two(2);
  ComplexBall a;
  acb_sqrt(a.raw(), two.raw(), 256);
  const auto exact_record = diffexp2::saturate_finite_laurent_basis<Rational>(
      {{rational_constant_frame("1"), rational_constant_frame("0"),
        rational_constant_frame("0")},
       {rational_constant_frame("0"), rational_constant_frame("1"),
        rational_constant_frame("0")},
       {rational_constant_frame("0"), rational_constant_frame("0"),
        rational_constant_frame("1")}},
      "exact ambiguous-off-pivot control record");
  const auto a_plus_two = a + ComplexBall(2);
  const auto a_plus_five = a + ComplexBall(5);
  AcbLaurentRefinementOptions options;
  options.relative_tolerance = Magnitude::decimal("1e-50");
  options.required_complete_max = 4;
  const auto matched = diffexp2::refine_acb_finite_laurent_match(
      {{ball_value_frame(a), ball_constant_frame("1"),
        ball_constant_frame("0")},
       {ball_value_frame(a), ball_constant_frame("1"),
        ball_constant_frame("1")},
       {ball_constant_frame("0"), ball_constant_frame("1"),
        ball_constant_frame("1")}},
      {ball_value_frame(a_plus_two), ball_value_frame(a_plus_five),
       ball_constant_frame("5")},
      exact_record, options, "ambiguous off-pivot Acb match", true);
  check("certified pivot plan tolerates an ambiguous exact Schur zero",
        matched.residual_history.back().verdict ==
                AcbMatchingResidualVerdict::Pass &&
            (matched.weights[0].coefficient(0) - ComplexBall(1))
                .contains_zero() &&
            (matched.weights[1].coefficient(0) - ComplexBall(2))
                .contains_zero() &&
            (matched.weights[2].coefficient(0) - ComplexBall(3))
                .contains_zero());
}

void empty_residual_window_retry_metadata_smoke() {
  ComplexBall::set_precision(256);
  const auto short_frame = []() {
    std::vector<ComplexBall> coefficients;
    coefficients.emplace_back(1);
    coefficients.emplace_back(0);
    coefficients.emplace_back(0);
    return EpsilonFrame<ComplexBall>(-4, std::move(coefficients));
  };
  AcbLaurentRefinementOptions options;
  options.required_min_power = -1;
  options.required_complete_max = 4;
  bool reported = false;
  try {
    (void)diffexp2::matching_detail::evaluate_acb_matching_residual(
        FiniteLaurentMatrix<ComplexBall>{{short_frame()}},
        std::vector<EpsilonFrame<ComplexBall>>{
            ball_constant_frame("1")},
        std::vector<EpsilonFrame<ComplexBall>>{short_frame()}, options,
        "empty residual retry metadata");
  } catch (const MatchingArithmeticError& error) {
    reported =
        error.code == MatchingArithmeticErrorCode::InsufficientCompleteWindow &&
        error.epsilon_power.has_value() && *error.epsilon_power == -2 &&
        options.required_complete_max - *error.epsilon_power == 6;
  }
  check("empty Acb residual window reports the exact retry reservoir loss",
        reported);
}

void precomputed_physical_residual_certification_smoke() {
  ComplexBall::set_precision(256);
  const auto uncertain = real_ball_with_error(1, -20);
  const FiniteLaurentMatrix<ComplexBall> physical_basis{
      {ball_value_frame(uncertain)}};
  const FiniteLaurentVector<ComplexBall> physical_weights{
      ball_constant_frame("1")};
  const FiniteLaurentVector<ComplexBall> physical_incoming{
      ball_value_frame(uncertain)};
  AcbLaurentRefinementOptions options;
  options.relative_tolerance = Magnitude::decimal("1e-8");
  options.required_min_power = 0;
  options.required_complete_max = 0;

  const auto direct =
      diffexp2::matching_detail::evaluate_acb_matching_residual(
          physical_basis, physical_weights, physical_incoming, options,
          "dependency-heavy physical subtraction");
  const auto pushed =
      diffexp2::matching_detail::certify_precomputed_acb_matching_residual(
          physical_basis, physical_weights, physical_incoming,
          FiniteLaurentVector<ComplexBall>{ball_constant_frame("0")},
          options, "exact normal-frame residual pushforward");
  const auto pushed_fail =
      diffexp2::matching_detail::certify_precomputed_acb_matching_residual(
          physical_basis, physical_weights, physical_incoming,
          FiniteLaurentVector<ComplexBall>{ball_constant_frame("1e-2")},
          options, "incorrect normal-frame residual pushforward");
  const auto cancellation_scaled_fail =
      diffexp2::matching_detail::certify_precomputed_acb_matching_residual(
          FiniteLaurentMatrix<ComplexBall>{
              {ball_constant_frame("1e40")}},
          FiniteLaurentVector<ComplexBall>{
              ball_constant_frame("1")},
          FiniteLaurentVector<ComplexBall>{
              ball_constant_frame("1")},
          FiniteLaurentVector<ComplexBall>{
              ball_constant_frame("1e20")},
          options,
          "ill-conditioned cancellation cannot weaken publication accuracy");
  auto short_options = options;
  short_options.required_complete_max = 3;
  const auto short_scale =
      diffexp2::matching_detail::certify_precomputed_acb_matching_residual(
          FiniteLaurentMatrix<ComplexBall>{
              {ball_constant_frame("1", 2)}},
          FiniteLaurentVector<ComplexBall>{
              ball_constant_frame("1", 2)},
          FiniteLaurentVector<ComplexBall>{
              ball_constant_frame("1", 2)},
          FiniteLaurentVector<ComplexBall>{
              ball_constant_frame("0", 6)},
          short_options,
          "pushed residual cannot outlive physical scale frames");

  check("direct dependency-heavy physical residual is inconclusive",
        direct.diagnostics.verdict ==
            AcbMatchingResidualVerdict::Inconclusive);
  check("exactly pushed physical residual certifies against physical scale",
        pushed.diagnostics.verdict ==
                AcbMatchingResidualVerdict::Pass &&
            pushed.diagnostics.complete_through_required);
  check("a pushed residual with a certified physical discrepancy is rejected",
        pushed_fail.diagnostics.verdict ==
            AcbMatchingResidualVerdict::Fail);
  check("large individual contributions cannot authorize a large physical residual",
        cancellation_scaled_fail.diagnostics.verdict ==
            AcbMatchingResidualVerdict::Fail);
  check("a pushed residual cannot manufacture physical scale completeness",
        !short_scale.diagnostics.complete_through_required &&
            short_scale.diagnostics.complete_window.complete_max == 1 &&
            short_scale.diagnostics.verdict ==
                AcbMatchingResidualVerdict::Inconclusive);
}

struct NormalFramePrefixRun {
  std::vector<EpsilonFrame<ComplexBall>> physical_weights;
  diffexp2::AcbMatchingResidualDiagnostics normalized_residual;
  diffexp2::AcbMatchingResidualDiagnostics physical_residual;
};

NormalFramePrefixRun run_laurent_normal_frame_prefix(std::size_t width) {
  const auto zero = ball_constant_frame("0", width);
  const auto one = ball_constant_frame("1", width);
  const FiniteLaurentMatrix<ComplexBall> physical_basis = {
      {one, zero}, {zero, one}};
  std::vector<ComplexBall> first_coefficients;
  std::vector<ComplexBall> second_coefficients;
  first_coefficients.reserve(width);
  second_coefficients.reserve(width);
  for (std::size_t power = 0; power < width; ++power) {
    first_coefficients.emplace_back(static_cast<long>(power + 2));
    second_coefficients.emplace_back(static_cast<long>(2 * power + 3));
  }
  const std::vector<EpsilonFrame<ComplexBall>> physical_incoming{
      EpsilonFrame<ComplexBall>(0, std::move(first_coefficients)),
      EpsilonFrame<ComplexBall>(0, std::move(second_coefficients))};

  const auto exact_zero = ExactLaurentPolynomial<ComplexBall>::zero();
  const auto exact_one = ExactLaurentPolynomial<ComplexBall>::one();
  const auto inverse_epsilon =
      ExactLaurentPolynomial<ComplexBall>::monomial(
          -1, ComplexBall(1));
  const auto negative_inverse_epsilon =
      ExactLaurentPolynomial<ComplexBall>::monomial(
          -1, ComplexBall(-1));
  const ExactLaurentMatrix<ComplexBall> frame = {
      {exact_one, inverse_epsilon}, {exact_zero, exact_one}};
  const ExactLaurentMatrix<ComplexBall> inverse_frame = {
      {exact_one, negative_inverse_epsilon}, {exact_zero, exact_one}};

  FiniteLaurentMatrix<ComplexBall> left_normalized(
      2, FiniteLaurentVector<ComplexBall>());
  for (auto& row : left_normalized) row.reserve(2);
  for (std::size_t column = 0; column < 2; ++column) {
    const auto normalized_column = diffexp2::apply_exact_laurent_matrix(
        inverse_frame,
        {physical_basis[0][column], physical_basis[1][column]});
    for (std::size_t row = 0; row < 2; ++row)
      left_normalized[row].push_back(normalized_column[row]);
  }
  const auto normalized_basis =
      diffexp2::right_multiply_finite_by_exact_laurent(
          left_normalized, frame);
  const auto normalized_incoming = diffexp2::apply_exact_laurent_matrix(
      inverse_frame, physical_incoming);

  AcbLaurentRefinementOptions options;
  options.relative_tolerance = Magnitude::decimal("1e-50");
  options.required_min_power = -1;
  options.required_complete_max = 1;
  options.max_refinement_steps = 2;
  const auto matched = diffexp2::refine_acb_finite_laurent_match(
      normalized_basis, normalized_incoming,
      diffexp2::identity_exact_laurent_matrix<ComplexBall>(2), options,
      "two-width SCC Laurent normal frame", true);
  auto physical_weights = diffexp2::apply_exact_laurent_matrix(
      frame, matched.weights);
  auto physical_options = options;
  physical_options.required_min_power = 0;
  const auto physical =
      diffexp2::matching_detail::evaluate_acb_matching_residual(
          physical_basis, physical_weights, physical_incoming,
          physical_options, "two-width physical prefix check");
  return {std::move(physical_weights),
          matched.residual_history.back(), physical.diagnostics};
}

void normal_frame_prefix_monotonicity_smoke() {
  ComplexBall::set_precision(256);
  const auto narrow = run_laurent_normal_frame_prefix(7);
  const auto wide = run_laurent_normal_frame_prefix(11);
  bool same_old_prefix = true;
  for (std::size_t component = 0; component < 2; ++component) {
    const auto old_top = narrow.physical_weights[component].complete_max();
    if (wide.physical_weights[component].complete_max() < old_top) {
      same_old_prefix = false;
      continue;
    }
    for (std::int32_t power = 0; power <= old_top; ++power)
      same_old_prefix = same_old_prefix &&
          (narrow.physical_weights[component].coefficient(power) -
           wide.physical_weights[component].coefficient(power)).is_zero();
  }
  check("widened Laurent normal-frame weights preserve the old physical prefix",
        same_old_prefix);
  check("widened Laurent normal frame cannot lower residual completeness",
        wide.normalized_residual.complete_window.complete_max >=
                narrow.normalized_residual.complete_window.complete_max &&
            wide.physical_residual.complete_window.complete_max >=
                narrow.physical_residual.complete_window.complete_max);
  check("physical-prefix gate rejects a normal frame beyond its honest edge",
        narrow.physical_residual.complete_through_required &&
            narrow.physical_residual.complete_window.complete_max < 6);
}

void verified_midpoint_preconditioner_smoke() {
  // A 13x13 Hilbert matrix at the minimum production precision is a compact
  // deterministic wrapping stress case.  Direct interval elimination cannot
  // find a complete certified pivot path, although the interval matrix is
  // regular.  A fixed dyadic midpoint-inverse row transform exposes that
  // regularity without classifying any zero-containing coefficient.
  ComplexBall::set_precision(64);
  constexpr std::size_t size = 13;
  FiniteLaurentMatrix<ComplexBall> matrix(
      size, std::vector<EpsilonFrame<ComplexBall>>());
  for (std::size_t row = 0; row < size; ++row) {
    matrix[row].reserve(size);
    for (std::size_t column = 0; column < size; ++column)
      matrix[row].push_back(ball_constant_frame(
          "1/" + std::to_string(row + column + 1), 40));
  }
  const auto direct = diffexp2::matching_detail::certified_full_rank_plan(
      diffexp2::matching_detail::epsilon_zero_matrix(
          matrix, "Hilbert direct leading matrix"),
      0);
  const auto factorization =
      diffexp2::factor_exact_nonnegative_finite_laurent_system(
          matrix, "verified Hilbert midpoint preconditioning");
  std::vector<EpsilonFrame<ComplexBall>> ones(
      size, ball_constant_frame("1", 40));
  const auto right_hand_side =
      diffexp2::apply_finite_laurent_matrix(matrix, ones);
  const auto solved = diffexp2::solve_factorized_finite_laurent_system(
      factorization, right_hand_side,
      "verified Hilbert preconditioned solve");
  const bool encloses_expected = std::all_of(
      solved.begin(), solved.end(), [](const auto& value) {
        return (value.coefficient(0) - ComplexBall(1)).contains_zero();
      });
  check("verified midpoint inverse rescues an Acb wrapping-only rank failure",
        !direct.has_value() &&
            factorization.left_preconditioner.has_value() &&
            encloses_expected);
  ComplexBall::set_precision(256);
}

void certified_pivot_quality_and_parity_smoke() {
  ComplexBall::set_precision(256);
  const auto weak = real_ball_with_error(1, -12);
  const auto strong = real_ball_with_error(1, -100);
  const auto ranked = diffexp2::matching_detail::certified_full_rank_plan(
      diffexp2::matching_detail::DenseScalarMatrix<ComplexBall>{
          {weak, ComplexBall(0)}, {ComplexBall(0), strong}},
      0);
  check("certified planner ranks tighter relative-radius pivot first",
        ranked.has_value() && ranked->front() == std::pair{1UL, 1UL});

  // Property comparison against the former exhaustive row/column DFS.  The
  // matrices deliberately mix exact zeros, zero-overlapping balls, and
  // certified nonzero balls; all are small enough for exhaustive reference.
  std::uint64_t state = 0x9e3779b97f4a7c15ULL;
  bool parity = true;
  for (std::size_t sample = 0; sample < 40 && parity; ++sample) {
    diffexp2::matching_detail::DenseScalarMatrix<ComplexBall> matrix(
        3, std::vector<ComplexBall>(3, ComplexBall(0)));
    for (std::size_t row = 0; row < 3; ++row) {
      for (std::size_t column = 0; column < 3; ++column) {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        const auto midpoint = static_cast<long>((state >> 32) % 7) - 3;
        if (midpoint == 0 && ((state >> 8) & 1U) != 0)
          matrix[row][column] = real_ball_with_error(0, -80);
        else if (midpoint != 0 && ((state >> 9) & 1U) != 0)
          matrix[row][column] = real_ball_with_error(midpoint, -80);
        else
          matrix[row][column] = ComplexBall(midpoint);
      }
    }
    const auto reference = reference_exhaustive_full_rank_plan(matrix);
    const auto planned =
        diffexp2::matching_detail::certified_full_rank_plan(matrix, 0);
    parity = reference.has_value() == planned.has_value();
  }
  check("small Acb property sweep matches exhaustive planner existence",
        parity);
}

void dense_eleven_by_eleven_pivot_budget_smoke() {
  ComplexBall::set_precision(256);
  constexpr std::size_t size = 11;
  diffexp2::matching_detail::DenseScalarMatrix<ComplexBall> matrix(
      size, std::vector<ComplexBall>(size, ComplexBall(0)));
  for (std::size_t row = 0; row + 1 < size; ++row) {
    for (std::size_t column = 0; column < size; ++column) {
      const auto midpoint = static_cast<long>(
          1 + (row == column ? 1 : 0));
      matrix[row][column] = real_ball_with_error(midpoint, -180);
    }
  }
  // A duplicated final row with independent tiny radii is the adversarial
  // interval analogue of a dense rank-10 matrix.  It offers certified pivots
  // almost everywhere but no certifiable last pivot, which made the old DFS
  // enumerate pivot permutations.
  for (std::size_t column = 0; column < size; ++column)
    matrix.back()[column] = matrix[size - 2][column];

  const auto start = std::chrono::steady_clock::now();
  const auto plan =
      diffexp2::matching_detail::certified_full_rank_plan(matrix, 0);
  const auto proof =
      diffexp2::matching_detail::certified_full_rank_search(matrix, 0);
  const auto elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
  std::cout << "  INFO: dense 11x11 bounded pivot planning " << elapsed
            << " s\n";
  check("dense 11x11 ambiguous-rank planner is globally bounded",
        !plan.has_value() &&
            proof == diffexp2::matching_detail::FullRankProofResult::
                         SearchBudgetExhausted &&
            elapsed < 2.0);
}

}  // namespace

int main() {
  transformation_support_smoke();
  quotient_and_solve_smoke();
  coefficientwise_power_series_solve_smoke();
  epsilon_lattice_saturation_smoke();
  acb_saturation_candidate_chop_smoke();
  ambiguous_acb_pivot_smoke();
  refined_acb_match_smoke();
  refined_acb_match_without_public_upper_slack_smoke();
  exact_right_frame_residual_smoke();
  retained_transformed_basis_authority_smoke();
  certified_tiny_physical_weight_publication_smoke();
  structural_transformed_weight_floor_smoke();
  ill_scaled_refinement_smoke();
  incomplete_refinement_rollback_smoke();
  refined_acb_ambiguous_pivot_smoke();
  laurent_off_pivot_ambiguity_smoke();
  refined_acb_ambiguous_off_pivot_smoke();
  empty_residual_window_retry_metadata_smoke();
  precomputed_physical_residual_certification_smoke();
  normal_frame_prefix_monotonicity_smoke();
  verified_midpoint_preconditioner_smoke();
  certified_pivot_quality_and_parity_smoke();
  dense_eleven_by_eleven_pivot_budget_smoke();
  std::cout << "Results: " << (checked - failed) << " / " << checked
            << " tests passed\n";
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
