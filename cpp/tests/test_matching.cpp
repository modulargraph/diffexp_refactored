#include "diffexp2/matching.hpp"

#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <string>
#include <vector>

using diffexp2::ComplexBall;
using diffexp2::EpsilonFrame;
using diffexp2::ExactLaurentMatrix;
using diffexp2::ExactLaurentPolynomial;
using diffexp2::FiniteLaurentMatrix;
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
}

void quotient_and_solve_smoke() {
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
  check("certified-zero below pivot does not shorten completeness",
        triangular[0].complete_max() == 10 &&
            triangular[1].complete_max() == 10 &&
            triangular[0].coefficient(0) == Rational(1) &&
            triangular[1].coefficient(0) == Rational(2));
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
  // This exercises the bounded correction replay and must remain honest by
  // returning Inconclusive rather than manufacturing a numerical zero.
  options.relative_tolerance = Magnitude::zero();
  options.required_complete_max = 4;
  options.max_refinement_steps = 2;
  const auto matched = diffexp2::refine_acb_finite_laurent_match(
      {{ball_constant_frame(big), ball_constant_frame("1")},
       {ball_constant_frame("1"), ball_constant_frame(twice_over_big)}},
      {ball_constant_frame("0"), ball_constant_frame("-1")}, exact_record,
      options, "ill-scaled refined Acb match");
  const auto expected_large = ComplexBall::from_strings("-" + big);
  check("ill-scaled Acb solve reuses bounded refinement honestly",
        matched.refinement_steps == 2 &&
            matched.residual_history.size() == 3 &&
            matched.residual_history.back().verdict ==
                AcbMatchingResidualVerdict::Inconclusive &&
            matched.residual_history.back().complete_through_required &&
            (matched.weights[0].coefficient(0) - ComplexBall(1))
                .contains_zero() &&
            (matched.weights[1].coefficient(0) - expected_large)
                .contains_zero());
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

}  // namespace

int main() {
  transformation_support_smoke();
  quotient_and_solve_smoke();
  epsilon_lattice_saturation_smoke();
  ambiguous_acb_pivot_smoke();
  refined_acb_match_smoke();
  ill_scaled_refinement_smoke();
  refined_acb_ambiguous_pivot_smoke();
  refined_acb_ambiguous_off_pivot_smoke();
  std::cout << "Results: " << (checked - failed) << " / " << checked
            << " tests passed\n";
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
