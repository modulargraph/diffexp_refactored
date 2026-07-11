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
using diffexp2::MatchingArithmeticError;
using diffexp2::MatchingArithmeticErrorCode;
using diffexp2::Rational;

namespace {

int failed = 0;

void check(const std::string& label, bool condition) {
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
}

}  // namespace

int main() {
  transformation_support_smoke();
  quotient_and_solve_smoke();
  ambiguous_acb_pivot_smoke();
  std::cout << "Results: " << (6 - failed) << " / 6 tests passed\n";
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
