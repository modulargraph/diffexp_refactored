#pragma once

#include "diffexp2/local_algebra.hpp"
#include "diffexp2/physical_ode.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace diffexp2 {

// Finite-epsilon Taylor problem for the backward observable equation
//
//             q(t,eps) theta(lambda) + C(t,eps)^T lambda = b(t,eps).
//
// The endpoint condition lambda(0)=0 is structural: coefficients start at
// t^1.  This is the composed alternative to separately truncating a terminal
// basis A, a line row J, and the ill-conditioned connection solve A^-1 u.
// Each q/C/forcing entry is an honest finite Laurent frame.  A resonant layer
// is rejected by the certified Laurent factorization; logarithmic completion
// must be supplied by the later singular extension, never guessed here.
struct BackwardAdjointTaylorProblem {
  std::uint32_t dimension = 0;
  std::uint32_t taylor_complete_max = 0;
  std::int32_t required_epsilon_complete_max = 0;
  std::vector<EpsilonFrame<ComplexBall>> q_lags;
  std::vector<FiniteLaurentMatrix<ComplexBall>> c_transpose_lags;
  // forcing[n-1] is the coefficient of t^n.  A short forcing is completed by
  // structural zeros, while coefficients beyond taylor_complete_max are not
  // consumed by this finite recurrence.
  std::vector<FiniteLaurentVector<ComplexBall>> forcing;
};

struct BackwardAdjointTaylorResult {
  std::uint32_t dimension = 0;
  std::uint32_t taylor_complete_max = 0;
  std::int32_t common_epsilon_complete_max = 0;
  // coefficients[n-1][component] is lambda_{component,n}, n >= 1.
  std::vector<FiniteLaurentVector<ComplexBall>> coefficients;
  // higher_log_coefficients[n-1][k-1][component] is the coefficient of
  // t^n (epsilon Log[t])^k/k!, k >= 1.  Keeping the normalized logarithm
  // convention identical to LocalSolution is essential: every log raise
  // consumes one lower epsilon order and remains visible to completeness.
  std::vector<std::vector<FiniteLaurentVector<ComplexBall>>>
      higher_log_coefficients;
  std::uint32_t max_log_power = 0;

  // This is only a planning signal.  It is not a tail certificate and must
  // never be used to publish a full-local result.
  std::optional<double> last_coefficient_ratio_upper;
};

struct BackwardAdjointTailCertificate {
  Magnitude absolute_vector_tail_upper = Magnitude::zero();
  Magnitude coefficient_majorant_upper = Magnitude::zero();
  Magnitude recurrence_contraction_upper = Magnitude::zero();
  Magnitude evaluation_ratio_upper = Magnitude::zero();
  std::uint32_t certified_after_taylor_order = 0;
};

struct BackwardAdjointAdaptiveTailCertificate {
  BackwardAdjointTailCertificate tail;
  Rational witness_radius_exact = Rational(0);
};

namespace adjoint_observable_detail {

inline EpsilonFrame<ComplexBall> structural_zero_frame(
    std::int32_t complete_max) {
  return EpsilonFrame<ComplexBall>::zero(complete_max);
}

inline FiniteLaurentVector<ComplexBall> structural_zero_vector(
    std::uint32_t dimension, std::int32_t complete_max) {
  FiniteLaurentVector<ComplexBall> result;
  result.reserve(dimension);
  for (std::uint32_t component = 0; component < dimension; ++component)
    result.push_back(structural_zero_frame(complete_max));
  return result;
}

inline void require_matrix_shape(
    const FiniteLaurentMatrix<ComplexBall>& matrix,
    std::uint32_t dimension, const char* label) {
  if (matrix.size() != dimension)
    throw std::invalid_argument(std::string(label) +
                                " has the wrong row count");
  for (const auto& row : matrix)
    if (row.size() != dimension)
      throw std::invalid_argument(std::string(label) +
                                  " has the wrong column count");
}

inline EpsilonFrame<ComplexBall> scaled_integer(
    const EpsilonFrame<ComplexBall>& value, std::int64_t factor) {
  return value.scaled(ComplexBall::from_strings(std::to_string(factor)));
}

inline FiniteLaurentVector<ComplexBall> apply_lag_operator(
    const EpsilonFrame<ComplexBall>& q,
    const FiniteLaurentMatrix<ComplexBall>* c_transpose,
    std::int64_t theta_factor,
    const FiniteLaurentVector<ComplexBall>& input,
    std::int32_t structural_complete_max) {
  const auto dimension = static_cast<std::uint32_t>(input.size());
  auto output = structural_zero_vector(dimension, structural_complete_max);
  for (std::uint32_t row = 0; row < dimension; ++row) {
    if (theta_factor != 0)
      output[row] = output[row] +
          scaled_integer(q, theta_factor) * input[row];
    if (c_transpose == nullptr) continue;
    for (std::uint32_t column = 0; column < dimension; ++column)
      output[row] = output[row] +
          (*c_transpose)[row][column] * input[column];
  }
  return output;
}

inline FiniteLaurentVector<ComplexBall> epsilon_shift_vector(
    const FiniteLaurentVector<ComplexBall>& input) {
  FiniteLaurentVector<ComplexBall> result;
  result.reserve(input.size());
  for (const auto& component : input)
    result.push_back(component.shifted(1));
  return result;
}

inline FiniteLaurentVector<ComplexBall> apply_log_lag_operator(
    const EpsilonFrame<ComplexBall>& q,
    const FiniteLaurentMatrix<ComplexBall>* c_transpose,
    std::int64_t theta_factor,
    const FiniteLaurentVector<ComplexBall>& current_log,
    const FiniteLaurentVector<ComplexBall>* next_log,
    std::int32_t structural_complete_max) {
  auto output = apply_lag_operator(
      q, c_transpose, theta_factor, current_log,
      structural_complete_max);
  if (next_log == nullptr) return output;
  const auto shifted = epsilon_shift_vector(*next_log);
  const auto log_derivative = apply_lag_operator(
      q, nullptr, 1, shifted, structural_complete_max);
  for (std::size_t component = 0; component < output.size(); ++component)
    output[component] = output[component] + log_derivative[component];
  return output;
}

inline void subtract_vector_in_place(
    FiniteLaurentVector<ComplexBall>& target,
    const FiniteLaurentVector<ComplexBall>& value) {
  if (target.size() != value.size())
    throw std::invalid_argument(
        "backward adjoint vector dimensions disagree");
  for (std::size_t component = 0; component < target.size(); ++component)
    target[component] = target[component] - value[component];
}

inline void add_vector_in_place(
    FiniteLaurentVector<ComplexBall>& target,
    const FiniteLaurentVector<ComplexBall>& value) {
  if (target.size() != value.size())
    throw std::invalid_argument(
        "backward adjoint vector dimensions disagree");
  for (std::size_t component = 0; component < target.size(); ++component)
    target[component] = target[component] + value[component];
}

inline double vector_upper_norm(
    const FiniteLaurentVector<ComplexBall>& value) {
  double result = 0.0;
  for (const auto& component : value)
    for (std::int64_t raw_power = component.min_power();
         raw_power <= component.complete_max(); ++raw_power)
      result = std::max(
          result,
          Magnitude::upper_abs(component.coefficient(
              static_cast<std::int32_t>(raw_power))).approximate_upper());
  return result;
}

template <typename Scalar>
EpsilonFrame<ComplexBall> expand_epsilon_rational(
    const ExactEpsilonRational<Scalar>& value,
    std::int32_t complete_max, const std::string& context) {
  physical_ode_detail::validate_rational(
      value, "backward adjoint epsilon-rational coefficient");
  if (value.zero || value.valuation > complete_max)
    return structural_zero_frame(complete_max);
  const auto width = EpsilonWindow{
      value.valuation, complete_max}.width();
  std::vector<ComplexBall> coefficients(width, ComplexBall(0));
  const auto denominator0 = local_detail::to_ball(value.denominator.front());
  if (denominator0.contains_zero())
    throw std::domain_error(
        context + ": epsilon-rational denominator contains zero");
  for (std::size_t offset = 0; offset < width; ++offset) {
    ComplexBall coefficient(0);
    if (offset < value.numerator.size())
      coefficient = local_detail::to_ball(value.numerator[offset]);
    for (std::size_t degree = 1;
         degree < value.denominator.size() && degree <= offset; ++degree)
      coefficient -= local_detail::to_ball(value.denominator[degree]) *
                     coefficients[offset - degree];
    coefficients[offset] = coefficient / denominator0;
  }
  return EpsilonFrame<ComplexBall>(value.valuation,
                                   std::move(coefficients));
}

inline EpsilonFrame<Rational> expand_exact_epsilon_rational(
    const ExactEpsilonRational<Rational>& value,
    std::int32_t complete_max, const std::string& context) {
  physical_ode_detail::validate_rational(
      value, "exact backward adjoint epsilon-rational coefficient");
  if (value.zero || value.valuation > complete_max)
    return EpsilonFrame<Rational>::zero(complete_max);
  const auto width = EpsilonWindow{
      value.valuation, complete_max}.width();
  std::vector<Rational> coefficients(width, Rational(0));
  const auto& denominator0 = value.denominator.front();
  if (denominator0.is_zero())
    throw std::domain_error(
        context + ": exact epsilon-rational denominator is zero");
  for (std::size_t offset = 0; offset < width; ++offset) {
    Rational coefficient(0);
    if (offset < value.numerator.size())
      coefficient = value.numerator[offset];
    for (std::size_t degree = 1;
         degree < value.denominator.size() && degree <= offset; ++degree)
      coefficient -= value.denominator[degree] *
                     coefficients[offset - degree];
    coefficients[offset] = coefficient / denominator0;
  }
  return EpsilonFrame<Rational>(value.valuation,
                                std::move(coefficients));
}

inline FiniteLaurentMatrix<Rational> exact_backward_adjoint_layer(
    const PreparedPhysicalClearedODE<Rational>& ode,
    std::uint32_t taylor_order, std::int32_t complete_max,
    const std::string& context) {
  physical_ode_detail::validate_ode(ode);
  if (ode.q_lags.empty() || ode.c_lags.empty())
    throw std::invalid_argument(
        context + ": exact physical equation has no indicial layer");
  auto q0 = expand_exact_epsilon_rational(
      ode.q_lags.front(), complete_max, context + ": exact q0");
  FiniteLaurentMatrix<Rational> layer(
      ode.dimension, FiniteLaurentVector<Rational>());
  for (auto& row : layer)
    for (std::uint32_t column = 0; column < ode.dimension; ++column)
      row.push_back(EpsilonFrame<Rational>::zero(complete_max));
  for (const auto& entry : ode.c_lags.front())
    layer[entry.column][entry.row] = expand_exact_epsilon_rational(
        entry.value, complete_max, context + ": exact C0 transpose");
  const Rational order(std::to_string(taylor_order));
  for (std::uint32_t diagonal = 0; diagonal < ode.dimension; ++diagonal)
    layer[diagonal][diagonal] =
        layer[diagonal][diagonal] + q0.scaled(order);
  return layer;
}

inline void apply_exact_backward_adjoint_layer_zeros(
    FiniteLaurentMatrix<ComplexBall>& numeric,
    const FiniteLaurentMatrix<Rational>& exact,
    const std::string& context) {
  if (numeric.size() != exact.size())
    throw std::invalid_argument(
        context + ": numeric and exact layer dimensions differ");
  for (std::size_t row = 0; row < numeric.size(); ++row) {
    if (numeric[row].size() != exact[row].size())
      throw std::invalid_argument(
          context + ": numeric and exact layer shapes differ");
    for (std::size_t column = 0; column < numeric[row].size(); ++column) {
      const auto leading = matching_detail::certified_laurent_leading_power(
          exact[row][column]);
      if (leading.first_ambiguous_power.has_value())
        throw std::logic_error(
            context + ": exact Rational layer has an ambiguous zero");
      if (!leading.power.has_value()) {
        numeric[row][column] = EpsilonFrame<ComplexBall>::zero(
            numeric[row][column].complete_max());
        continue;
      }
      const auto power = *leading.power;
      if (power < numeric[row][column].min_power() ||
          power > numeric[row][column].complete_max())
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            context +
                ": numeric layer does not contain its exact leading power",
            row, column, power);
      std::vector<ComplexBall> coefficients;
      coefficients.reserve(EpsilonWindow{
          power, numeric[row][column].complete_max()}.width());
      for (std::int64_t raw_power = power;
           raw_power <= numeric[row][column].complete_max(); ++raw_power)
        coefficients.push_back(numeric[row][column].coefficient(
            static_cast<std::int32_t>(raw_power)));
      numeric[row][column] = EpsilonFrame<ComplexBall>(
          power, std::move(coefficients));
    }
  }
}

struct ExactSimpleResonanceFrame {
  FiniteLaurentVector<Rational> right_kernel;
  FiniteLaurentVector<Rational> left_kernel;
};

inline FiniteLaurentMatrix<Rational> transpose_exact_matrix(
    const FiniteLaurentMatrix<Rational>& matrix) {
  const auto dimension = matrix.size();
  FiniteLaurentMatrix<Rational> result(
      dimension, FiniteLaurentVector<Rational>());
  for (auto& row : result) row.reserve(dimension);
  for (std::size_t row = 0; row < dimension; ++row) {
    if (matrix[row].size() != dimension)
      throw std::invalid_argument(
          "exact backward adjoint layer is not square");
    for (std::size_t column = 0; column < dimension; ++column)
      result[row].push_back(matrix[column][row]);
  }
  return result;
}

inline matching_detail::DenseScalarMatrix<Rational>
epsilon_zero_matrix(const FiniteLaurentMatrix<Rational>& matrix,
                    const std::string& context) {
  matching_detail::DenseScalarMatrix<Rational> result(
      matrix.size(), std::vector<Rational>(matrix.size(), Rational(0)));
  for (std::size_t row = 0; row < matrix.size(); ++row) {
    if (matrix[row].size() != matrix.size())
      throw std::invalid_argument(context + ": matrix is not square");
    for (std::size_t column = 0; column < matrix.size(); ++column) {
      if (matrix[row][column].complete_max() < 0)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            context + ": exact layer does not contain epsilon^0",
            row, column, matrix[row][column].complete_max());
      result[row][column] = matrix[row][column].coefficient(0);
    }
  }
  return result;
}

inline void require_exact_zero_vector(
    const FiniteLaurentVector<Rational>& value,
    const std::string& context) {
  for (std::size_t component = 0; component < value.size(); ++component)
    for (std::int64_t raw_power = value[component].min_power();
         raw_power <= value[component].complete_max(); ++raw_power)
      if (!value[component]
               .coefficient(static_cast<std::int32_t>(raw_power))
               .is_zero())
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::SaturationFailure,
            context + ": exact formal null-vector residual is nonzero",
            component, std::nullopt,
            static_cast<std::int32_t>(raw_power));
}

inline FiniteLaurentVector<Rational> exact_formal_null_vector(
    const FiniteLaurentMatrix<Rational>& matrix,
    const matching_detail::LeadingNullRelation<Rational>& leading,
    const std::vector<Rational>& dual_normalizer,
    const std::string& context) {
  const auto dimension = matrix.size();
  if (leading.rank + 1 != dimension || !leading.vector.has_value() ||
      leading.pivot_rows.size() != leading.rank ||
      dual_normalizer.size() != dimension)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::SingularOrIncompleteSystem,
        context +
            ": only a one-dimensional semisimple resonant kernel is supported");
  auto complete_max = matrix.front().front().complete_max();
  for (const auto& row : matrix)
    for (const auto& entry : row)
      complete_max = std::min(complete_max, entry.complete_max());
  FiniteLaurentMatrix<Rational> bordered;
  bordered.reserve(dimension);
  for (const auto pivot_row : leading.pivot_rows)
    bordered.push_back(matrix.at(pivot_row));
  FiniteLaurentVector<Rational> normalization;
  normalization.reserve(dimension);
  for (const auto& coefficient : dual_normalizer) {
    std::vector<Rational> coefficients(
        static_cast<std::size_t>(complete_max) + 1, Rational(0));
    coefficients.front() = coefficient;
    normalization.emplace_back(0, std::move(coefficients));
  }
  bordered.push_back(std::move(normalization));
  FiniteLaurentVector<Rational> exact_rhs;
  exact_rhs.reserve(dimension);
  for (std::size_t row = 0; row < dimension; ++row) {
    std::vector<Rational> coefficients(
        static_cast<std::size_t>(complete_max) + 1, Rational(0));
    if (row + 1 == dimension) coefficients.front() = Rational(1);
    exact_rhs.emplace_back(0, std::move(coefficients));
  }
  auto factorization = factor_finite_laurent_system(
      std::move(bordered), context + ": normalized null-vector solve");
  auto result = solve_factorized_finite_laurent_system(
      factorization, std::move(exact_rhs),
      context + ": normalized null-vector coefficients");
  require_exact_zero_vector(
      apply_finite_laurent_matrix(matrix, result),
      context + ": full null-vector audit");
  return result;
}

inline std::optional<ExactSimpleResonanceFrame>
exact_simple_resonance_frame(
    const FiniteLaurentMatrix<Rational>& layer,
    const std::string& context) {
  try {
    (void)factor_finite_laurent_system(
        layer, context + ": exact full-rank classification");
    return std::nullopt;
  } catch (const MatchingArithmeticError& error) {
    if (error.code !=
        MatchingArithmeticErrorCode::SingularOrIncompleteSystem)
      throw;
  }
  const auto leading_matrix = epsilon_zero_matrix(layer, context);
  const auto right = matching_detail::leading_null_relation(
      leading_matrix, context + ": exact right kernel");
  const auto transposed_layer = transpose_exact_matrix(layer);
  const auto left = matching_detail::leading_null_relation(
      epsilon_zero_matrix(transposed_layer, context),
      context + ": exact left kernel");
  if (right.rank + 1 != layer.size() ||
      left.rank + 1 != layer.size() || !right.vector.has_value() ||
      !left.vector.has_value())
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::SingularOrIncompleteSystem,
        context +
            ": resonant layer is not a simple one-dimensional kernel");
  Rational pairing(0);
  for (std::size_t component = 0; component < layer.size(); ++component)
    pairing += (*left.vector)[component] * (*right.vector)[component];
  if (pairing.is_zero())
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::SingularOrIncompleteSystem,
        context +
            ": resonant zero eigenvalue is nonsemisimple and needs the exact Jordan ladder");
  return ExactSimpleResonanceFrame{
      exact_formal_null_vector(
          layer, right, *left.vector,
          context + ": exact right formal kernel"),
      exact_formal_null_vector(
          transposed_layer, left, *right.vector,
          context + ": exact left formal kernel")};
}

inline EpsilonFrame<ComplexBall> exact_frame_to_ball(
    const EpsilonFrame<Rational>& exact) {
  std::vector<ComplexBall> coefficients;
  coefficients.reserve(exact.coefficients().size());
  for (const auto& coefficient : exact.coefficients())
    coefficients.push_back(ComplexBall::from_strings(coefficient.str()));
  return EpsilonFrame<ComplexBall>(exact.window(),
                                   std::move(coefficients));
}

inline FiniteLaurentVector<ComplexBall> exact_vector_to_ball(
    const FiniteLaurentVector<Rational>& exact) {
  FiniteLaurentVector<ComplexBall> result;
  result.reserve(exact.size());
  for (const auto& component : exact)
    result.push_back(exact_frame_to_ball(component));
  return result;
}

inline bool material_frame(const EpsilonFrame<ComplexBall>& value) {
  return std::any_of(
      value.coefficients().begin(), value.coefficients().end(),
      [](const auto& coefficient) { return !coefficient.is_zero(); });
}

inline Magnitude frame_l1_upper(
    const EpsilonFrame<ComplexBall>& value) {
  auto result = Magnitude::zero();
  for (const auto& coefficient : value.coefficients())
    result += Magnitude::upper_abs(coefficient);
  return result;
}

inline Magnitude vector_infinity_upper(
    const FiniteLaurentVector<ComplexBall>& value) {
  auto result = Magnitude::zero();
  for (const auto& component : value)
    for (const auto& coefficient : component.coefficients())
      result = Magnitude::maximum(
          result, Magnitude::upper_abs(coefficient));
  return result;
}

inline Magnitude matrix_infinity_upper(
    const FiniteLaurentMatrix<ComplexBall>& matrix) {
  auto result = Magnitude::zero();
  for (const auto& row : matrix) {
    auto row_sum = Magnitude::zero();
    for (const auto& entry : row) row_sum += frame_l1_upper(entry);
    result = Magnitude::maximum(result, row_sum);
  }
  return result;
}

inline Magnitude polynomial_disk_upper(
    const std::vector<ComplexBall>& coefficients,
    const Magnitude& radius_upper) {
  auto result = Magnitude::zero();
  auto power = Magnitude::one();
  for (const auto& coefficient : coefficients) {
    result += Magnitude::upper_abs(coefficient) * power;
    power = power * radius_upper;
  }
  return result;
}

inline Magnitude rational_disk_upper(
    const PreparedRationalAnalyticCoefficient<ComplexBall>& rational,
    const Magnitude& radius_upper, const std::string& context) {
  if (rational.numerator.empty() || rational.denominator.empty())
    throw std::invalid_argument(
        context + ": analytic rational coefficient is empty");
  const auto numerator =
      polynomial_disk_upper(rational.numerator, radius_upper);
  auto denominator_tail = Magnitude::zero();
  auto power = radius_upper;
  for (std::size_t degree = 1; degree < rational.denominator.size();
       ++degree) {
    denominator_tail +=
        Magnitude::upper_abs(rational.denominator[degree]) * power;
    power = power * radius_upper;
  }
  const auto denominator_lower = Magnitude::positive_difference_lower(
      Magnitude::lower_abs(rational.denominator.front()),
      denominator_tail);
  if (denominator_lower.is_zero())
    throw std::domain_error(
        context +
        ": witness disk does not prove the analytic denominator nonzero");
  return numerator / denominator_lower;
}

inline void trim_exact_t_polynomial(std::vector<Rational>& value) {
  while (value.size() > 1 && value.back().is_zero()) value.pop_back();
  if (value.empty()) value.emplace_back(0);
}

inline bool exact_t_polynomial_zero(const std::vector<Rational>& value) {
  return std::all_of(value.begin(), value.end(),
                     [](const Rational& coefficient) {
                       return coefficient.is_zero();
                     });
}

inline std::vector<Rational> multiply_exact_t_polynomials(
    const std::vector<Rational>& left,
    const std::vector<Rational>& right) {
  std::vector<Rational> result(
      left.size() + right.size() - 1, Rational(0));
  for (std::size_t i = 0; i < left.size(); ++i)
    for (std::size_t j = 0; j < right.size(); ++j)
      result[i + j] += left[i] * right[j];
  trim_exact_t_polynomial(result);
  return result;
}

inline std::vector<Rational> add_exact_t_polynomials(
    std::vector<Rational> left, const std::vector<Rational>& right) {
  left.resize(std::max(left.size(), right.size()), Rational(0));
  for (std::size_t i = 0; i < right.size(); ++i) left[i] += right[i];
  trim_exact_t_polynomial(left);
  return left;
}

inline std::pair<std::vector<Rational>, std::vector<Rational>>
divide_exact_t_polynomials(std::vector<Rational> numerator,
                           std::vector<Rational> denominator) {
  trim_exact_t_polynomial(numerator);
  trim_exact_t_polynomial(denominator);
  if (exact_t_polynomial_zero(denominator))
    throw std::domain_error("exact t-polynomial division by zero");
  if (numerator.size() < denominator.size())
    return {{Rational(0)}, std::move(numerator)};
  std::vector<Rational> quotient(
      numerator.size() - denominator.size() + 1, Rational(0));
  while (!exact_t_polynomial_zero(numerator) &&
         numerator.size() >= denominator.size()) {
    const auto shift = numerator.size() - denominator.size();
    const auto factor = numerator.back() / denominator.back();
    quotient[shift] += factor;
    for (std::size_t i = 0; i < denominator.size(); ++i)
      numerator[i + shift] -= factor * denominator[i];
    trim_exact_t_polynomial(numerator);
  }
  trim_exact_t_polynomial(quotient);
  return {std::move(quotient), std::move(numerator)};
}

inline std::vector<Rational> gcd_exact_t_polynomials(
    std::vector<Rational> left, std::vector<Rational> right) {
  trim_exact_t_polynomial(left);
  trim_exact_t_polynomial(right);
  while (!exact_t_polynomial_zero(right)) {
    auto remainder = divide_exact_t_polynomials(left, right).second;
    left = std::move(right);
    right = std::move(remainder);
  }
  if (exact_t_polynomial_zero(left)) return {Rational(1)};
  const auto leading = left.back();
  for (auto& coefficient : left) coefficient = coefficient / leading;
  trim_exact_t_polynomial(left);
  return left;
}

struct ExactTRationalFunction {
  std::vector<Rational> numerator{Rational(0)};
  std::vector<Rational> denominator{Rational(1)};
};

inline Magnitude exact_t_rational_disk_upper(
    const ExactTRationalFunction& exact,
    const Magnitude& radius_upper, const std::string& context);

inline ExactTRationalFunction canonical_exact_t_rational(
    ExactTRationalFunction value) {
  trim_exact_t_polynomial(value.numerator);
  trim_exact_t_polynomial(value.denominator);
  if (exact_t_polynomial_zero(value.denominator))
    throw std::domain_error("exact t-rational denominator is zero");
  if (exact_t_polynomial_zero(value.numerator)) return {};
  const auto divisor = gcd_exact_t_polynomials(
      value.numerator, value.denominator);
  if (!(divisor.size() == 1 && divisor.front() == Rational(1))) {
    auto numerator_division = divide_exact_t_polynomials(
        value.numerator, divisor);
    auto denominator_division = divide_exact_t_polynomials(
        value.denominator, divisor);
    if (!exact_t_polynomial_zero(numerator_division.second) ||
        !exact_t_polynomial_zero(denominator_division.second))
      throw std::logic_error("exact t-rational gcd did not divide exactly");
    value.numerator = std::move(numerator_division.first);
    value.denominator = std::move(denominator_division.first);
  }
  return value;
}

inline ExactTRationalFunction add_exact_t_rationals(
    const ExactTRationalFunction& left,
    const ExactTRationalFunction& right) {
  return canonical_exact_t_rational({
      add_exact_t_polynomials(
          multiply_exact_t_polynomials(left.numerator, right.denominator),
          multiply_exact_t_polynomials(right.numerator, left.denominator)),
      multiply_exact_t_polynomials(left.denominator, right.denominator)});
}

inline ExactTRationalFunction negate_exact_t_rational(
    ExactTRationalFunction value) {
  for (auto& coefficient : value.numerator) coefficient = -coefficient;
  return value;
}

inline ExactTRationalFunction subtract_exact_t_rationals(
    const ExactTRationalFunction& left,
    const ExactTRationalFunction& right) {
  return add_exact_t_rationals(left, negate_exact_t_rational(right));
}

inline ExactTRationalFunction multiply_exact_t_rationals(
    const ExactTRationalFunction& left,
    const ExactTRationalFunction& right) {
  return canonical_exact_t_rational({
      multiply_exact_t_polynomials(left.numerator, right.numerator),
      multiply_exact_t_polynomials(left.denominator, right.denominator)});
}

inline ExactTRationalFunction divide_exact_t_rationals(
    const ExactTRationalFunction& numerator,
    const ExactTRationalFunction& denominator) {
  if (exact_t_polynomial_zero(denominator.numerator))
    throw std::domain_error("exact t-rational division by zero");
  return canonical_exact_t_rational({
      multiply_exact_t_polynomials(
          numerator.numerator, denominator.denominator),
      multiply_exact_t_polynomials(
          numerator.denominator, denominator.numerator)});
}

inline std::vector<Rational> expand_exact_t_rational_taylor(
    const ExactTRationalFunction& exact, std::uint32_t complete_max,
    const std::string& context) {
  if (exact.denominator.empty() || exact.denominator.front().is_zero())
    throw std::domain_error(
        context + ": t-rational denominator vanishes at the center");
  std::vector<Rational> result(
      static_cast<std::size_t>(complete_max) + 1, Rational(0));
  for (std::uint32_t order = 0; order <= complete_max; ++order) {
    auto coefficient = order < exact.numerator.size()
        ? exact.numerator[order] : Rational(0);
    const auto maximum_lag = std::min<std::size_t>(
        order, exact.denominator.size() - 1);
    for (std::size_t lag = 1; lag <= maximum_lag; ++lag)
      coefficient -= exact.denominator[lag] * result[order - lag];
    result[order] = coefficient / exact.denominator.front();
  }
  return result;
}

struct NormalizedBackwardAdjointExactODE {
  PreparedPhysicalClearedODE<Rational> truncated_ode;
  // Physical orientation: c_by_entry[row][column][epsilon_power].  Every
  // entry is the exact analytic t-rational coefficient of C/q.
  std::vector<std::vector<std::vector<ExactTRationalFunction>>> c_by_entry;
  std::int32_t epsilon_complete_max = 0;
};

inline ExactEpsilonRational<Rational> exact_epsilon_polynomial(
    const std::vector<Rational>& coefficients) {
  auto first = coefficients.size();
  while (first > 0 && coefficients[first - 1].is_zero()) --first;
  if (first == 0) return {};
  first = 0;
  while (first < coefficients.size() && coefficients[first].is_zero())
    ++first;
  ExactEpsilonRational<Rational> result;
  result.zero = false;
  result.valuation = matching_detail::checked_power(
      static_cast<std::int64_t>(first),
      "normalized backward adjoint epsilon valuation");
  result.numerator.assign(
      coefficients.begin() + static_cast<std::ptrdiff_t>(first),
      coefficients.end());
  result.denominator = {Rational(1)};
  return result;
}

inline NormalizedBackwardAdjointExactODE
normalize_backward_adjoint_exact_ode_by_q(
    const PreparedPhysicalClearedODE<Rational>& exact_ode,
    std::uint32_t taylor_complete_max,
    std::int32_t epsilon_complete_max,
    const std::string& context) {
  physical_ode_detail::validate_ode(exact_ode);
  if (epsilon_complete_max < 0)
    throw std::invalid_argument(
        context + ": normalized finite epsilon stack must begin at zero");
  const auto epsilon_width = static_cast<std::size_t>(
      static_cast<std::int64_t>(epsilon_complete_max) + 1);
  const auto dimension = exact_ode.dimension;

  std::vector<std::vector<Rational>> q_by_epsilon(
      epsilon_width,
      std::vector<Rational>(exact_ode.q_lags.size(), Rational(0)));
  for (std::size_t lag = 0; lag < exact_ode.q_lags.size(); ++lag) {
    const auto frame = expand_exact_epsilon_rational(
        exact_ode.q_lags[lag], epsilon_complete_max,
        context + ": q lag " + std::to_string(lag));
    if (frame.min_power() < 0)
      throw std::domain_error(
          context + ": q has negative epsilon coupling");
    for (std::int64_t raw_power = frame.min_power();
         raw_power <= frame.complete_max(); ++raw_power)
      q_by_epsilon[static_cast<std::size_t>(raw_power)][lag] =
          frame.coefficient(static_cast<std::int32_t>(raw_power));
  }
  for (auto& polynomial : q_by_epsilon) trim_exact_t_polynomial(polynomial);
  const ExactTRationalFunction q0{q_by_epsilon.front(), {Rational(1)}};
  if (exact_t_polynomial_zero(q0.numerator) ||
      q0.numerator.front().is_zero())
    throw std::domain_error(
        context + ": q(t,epsilon=0) is not a center unit");

  std::vector<std::vector<std::vector<std::vector<Rational>>>>
      c_polynomials(
          dimension,
          std::vector<std::vector<std::vector<Rational>>>(
              dimension,
              std::vector<std::vector<Rational>>(
                  epsilon_width, std::vector<Rational>(1, Rational(0)))));
  for (std::size_t lag = 0; lag < exact_ode.c_lags.size(); ++lag) {
    for (const auto& entry : exact_ode.c_lags[lag]) {
      const auto frame = expand_exact_epsilon_rational(
          entry.value, epsilon_complete_max,
          context + ": C lag " + std::to_string(lag));
      if (frame.min_power() < 0)
        throw std::domain_error(
            context + ": C has negative epsilon coupling");
      for (std::int64_t raw_power = frame.min_power();
           raw_power <= frame.complete_max(); ++raw_power) {
        auto& polynomial = c_polynomials[entry.row][entry.column]
            [static_cast<std::size_t>(raw_power)];
        if (polynomial.size() <= lag)
          polynomial.resize(lag + 1, Rational(0));
        polynomial[lag] +=
            frame.coefficient(static_cast<std::int32_t>(raw_power));
      }
    }
  }

  NormalizedBackwardAdjointExactODE result;
  result.epsilon_complete_max = epsilon_complete_max;
  result.c_by_entry.assign(
      dimension,
      std::vector<std::vector<ExactTRationalFunction>>(
          dimension,
          std::vector<ExactTRationalFunction>(epsilon_width)));
  for (std::uint32_t row = 0; row < dimension; ++row) {
    for (std::uint32_t column = 0; column < dimension; ++column) {
      for (std::size_t epsilon = 0; epsilon < epsilon_width; ++epsilon) {
        trim_exact_t_polynomial(c_polynomials[row][column][epsilon]);
        auto rhs = ExactTRationalFunction{
            c_polynomials[row][column][epsilon], {Rational(1)}};
        for (std::size_t q_power = 1; q_power <= epsilon; ++q_power) {
          const ExactTRationalFunction q_term{
              q_by_epsilon[q_power], {Rational(1)}};
          rhs = subtract_exact_t_rationals(
              rhs, multiply_exact_t_rationals(
                       q_term,
                       result.c_by_entry[row][column]
                           [epsilon - q_power]));
        }
        result.c_by_entry[row][column][epsilon] =
            divide_exact_t_rationals(rhs, q0);
      }
    }
  }

  auto& normalized = result.truncated_ode;
  normalized.dimension = dimension;
  normalized.owner_signature_identity =
      exact_ode.owner_signature_identity + ":adjoint-q-normalized";
  normalized.payload_identity =
      exact_ode.payload_identity + ":adjoint-q-normalized:t" +
      std::to_string(taylor_complete_max) + ":e" +
      std::to_string(epsilon_complete_max);
  normalized.exact_payload_record =
      "diffexp2-backward-adjoint-q-normalized-v1:" +
      normalized.payload_identity;
  ExactEpsilonRational<Rational> one;
  one.zero = false;
  one.valuation = 0;
  one.numerator = {Rational(1)};
  one.denominator = {Rational(1)};
  normalized.q_lags = {one};
  normalized.c_lags.resize(
      static_cast<std::size_t>(taylor_complete_max) + 1);
  for (std::uint32_t row = 0; row < dimension; ++row) {
    for (std::uint32_t column = 0; column < dimension; ++column) {
      std::vector<std::vector<Rational>> by_epsilon;
      by_epsilon.reserve(epsilon_width);
      for (std::size_t epsilon = 0; epsilon < epsilon_width; ++epsilon)
        by_epsilon.push_back(expand_exact_t_rational_taylor(
            result.c_by_entry[row][column][epsilon],
            taylor_complete_max,
            context + ": normalized C(" + std::to_string(row) + "," +
                std::to_string(column) + ") epsilon " +
                std::to_string(epsilon)));
      for (std::uint32_t lag = 0; lag <= taylor_complete_max; ++lag) {
        std::vector<Rational> epsilon_coefficients(epsilon_width, Rational(0));
        for (std::size_t epsilon = 0; epsilon < epsilon_width; ++epsilon)
          epsilon_coefficients[epsilon] = by_epsilon[epsilon][lag];
        auto value = exact_epsilon_polynomial(epsilon_coefficients);
        if (value.zero) continue;
        normalized.c_lags[lag].push_back(
            {row, column, std::move(value)});
      }
    }
  }
  physical_ode_detail::validate_ode(normalized);
  return result;
}

inline Magnitude normalized_backward_adjoint_c_tail_disk_upper(
    const NormalizedBackwardAdjointExactODE& normalized,
    const ComplexBall& witness_radius,
    std::int32_t epsilon_complete_max,
    const std::string& context) {
  if (epsilon_complete_max < 0 ||
      epsilon_complete_max > normalized.epsilon_complete_max)
    throw std::invalid_argument(
        context + ": normalized C tail epsilon cap is out of range");
  const auto radius_upper = Magnitude::upper_abs(witness_radius);
  auto matrix_upper = Magnitude::zero();
  // The recurrence contains C^T, hence a row sum here is a physical-column
  // sum in c_by_entry.
  for (std::size_t transpose_row = 0;
       transpose_row < normalized.c_by_entry.size(); ++transpose_row) {
    auto row_upper = Magnitude::zero();
    for (std::size_t physical_row = 0;
         physical_row < normalized.c_by_entry.size(); ++physical_row) {
      const auto& by_epsilon =
          normalized.c_by_entry[physical_row][transpose_row];
      for (std::int32_t epsilon = 0;
           epsilon <= epsilon_complete_max; ++epsilon) {
        const auto& coefficient =
            by_epsilon[static_cast<std::size_t>(epsilon)];
        if (coefficient.denominator.empty() ||
            coefficient.denominator.front().is_zero())
          throw std::domain_error(
              context + ": normalized C denominator vanishes at center");
        const auto constant = coefficient.numerator.front() /
                              coefficient.denominator.front();
        const auto positive_tail = subtract_exact_t_rationals(
            coefficient,
            ExactTRationalFunction{{constant}, {Rational(1)}});
        row_upper += exact_t_rational_disk_upper(
            positive_tail, radius_upper,
            context + ": transpose row " +
                std::to_string(transpose_row) + ": physical row " +
                std::to_string(physical_row) + ": epsilon " +
                std::to_string(epsilon));
      }
    }
    matrix_upper = Magnitude::maximum(matrix_upper, row_upper);
  }
  return matrix_upper;
}

inline Magnitude exact_t_rational_disk_upper(
    const ExactTRationalFunction& exact,
    const Magnitude& radius_upper, const std::string& context) {
  PreparedRationalAnalyticCoefficient<ComplexBall> numeric;
  for (const auto& coefficient : exact.numerator)
    numeric.numerator.push_back(
        ComplexBall::from_strings(coefficient.str()));
  for (const auto& coefficient : exact.denominator)
    numeric.denominator.push_back(
        ComplexBall::from_strings(coefficient.str()));
  return rational_disk_upper(numeric, radius_upper, context);
}

inline Magnitude exact_t_rational_taylor_tail_l1_upper(
    const ExactTRationalFunction& exact,
    std::uint32_t retained_complete_max,
    const Magnitude& radius_upper, const std::string& context) {
  const auto prefix = expand_exact_t_rational_taylor(
      exact, retained_complete_max, context + ": Taylor prefix");
  auto residual = add_exact_t_polynomials(
      exact.numerator,
      negate_exact_t_rational(ExactTRationalFunction{
          multiply_exact_t_polynomials(exact.denominator, prefix),
          {Rational(1)}}).numerator);
  const auto first_unseen =
      static_cast<std::size_t>(retained_complete_max) + 1;
  if (residual.size() <= first_unseen) return Magnitude::zero();
  for (std::size_t order = 0; order < first_unseen; ++order)
    if (!residual[order].is_zero())
      throw std::logic_error(
          context + ": exact rational Taylor remainder did not factor");
  residual.erase(
      residual.begin(),
      residual.begin() + static_cast<std::ptrdiff_t>(first_unseen));
  trim_exact_t_polynomial(residual);
  if (exact_t_polynomial_zero(residual)) return Magnitude::zero();

  std::vector<ComplexBall> residual_numeric;
  residual_numeric.reserve(residual.size());
  for (const auto& coefficient : residual)
    residual_numeric.push_back(
        ComplexBall::from_strings(coefficient.str()));
  std::vector<ComplexBall> denominator_numeric;
  denominator_numeric.reserve(exact.denominator.size());
  for (const auto& coefficient : exact.denominator)
    denominator_numeric.push_back(
        ComplexBall::from_strings(coefficient.str()));
  auto denominator_tail = Magnitude::zero();
  auto power = radius_upper;
  for (std::size_t degree = 1; degree < denominator_numeric.size(); ++degree) {
    denominator_tail +=
        Magnitude::upper_abs(denominator_numeric[degree]) * power;
    power = power * radius_upper;
  }
  const auto denominator_lower = Magnitude::positive_difference_lower(
      Magnitude::lower_abs(denominator_numeric.front()), denominator_tail);
  if (denominator_lower.is_zero())
    throw std::domain_error(
        context + ": witness disk does not separate the exact denominator");
  return radius_upper.power_upper(
             static_cast<ulong>(first_unseen)) *
         polynomial_disk_upper(residual_numeric, radius_upper) /
         denominator_lower;
}

inline Magnitude normalized_backward_adjoint_c_taylor_tail_l1_upper(
    const NormalizedBackwardAdjointExactODE& normalized,
    const ComplexBall& witness_radius,
    std::int32_t epsilon_complete_max,
    std::uint32_t retained_complete_max,
    const std::string& context) {
  if (epsilon_complete_max < 0 ||
      epsilon_complete_max > normalized.epsilon_complete_max)
    throw std::invalid_argument(
        context + ": normalized C Taylor-tail epsilon cap is out of range");
  const auto radius_upper = Magnitude::upper_abs(witness_radius);
  auto matrix_upper = Magnitude::zero();
  for (std::size_t transpose_row = 0;
       transpose_row < normalized.c_by_entry.size(); ++transpose_row) {
    auto row_upper = Magnitude::zero();
    for (std::size_t physical_row = 0;
         physical_row < normalized.c_by_entry.size(); ++physical_row) {
      const auto& by_epsilon =
          normalized.c_by_entry[physical_row][transpose_row];
      for (std::int32_t epsilon = 0;
           epsilon <= epsilon_complete_max; ++epsilon)
        row_upper += exact_t_rational_taylor_tail_l1_upper(
            by_epsilon[static_cast<std::size_t>(epsilon)],
            retained_complete_max, radius_upper,
            context + ": transpose row " +
                std::to_string(transpose_row) + ": physical row " +
                std::to_string(physical_row) + ": epsilon " +
                std::to_string(epsilon));
    }
    matrix_upper = Magnitude::maximum(matrix_upper, row_upper);
  }
  return matrix_upper;
}

inline Magnitude vector_infinity_upper_through(
    const FiniteLaurentVector<ComplexBall>& value,
    std::int32_t epsilon_complete_max) {
  auto result = Magnitude::zero();
  for (const auto& component : value)
    for (std::int64_t raw_power = component.min_power();
         raw_power <= std::min<std::int64_t>(
                          component.complete_max(),
                          epsilon_complete_max);
         ++raw_power)
      result = Magnitude::maximum(
          result, Magnitude::upper_abs(component.coefficient(
                      static_cast<std::int32_t>(raw_power))));
  return result;
}

inline Magnitude normalized_backward_adjoint_defect_l1_upper(
    const NormalizedBackwardAdjointExactODE& normalized,
    const PreparedSparseLocalMultiplierMatrix<Rational>& exact_row,
    const BackwardAdjointTaylorResult& solution,
    const ComplexBall& oriented_physical_jacobian,
    const ComplexBall& witness_radius,
    std::int32_t epsilon_complete_max,
    const std::string& context) {
  if (solution.taylor_complete_max == 0 ||
      solution.coefficients.size() != solution.taylor_complete_max ||
      exact_row.rows != 1 ||
      exact_row.columns != normalized.truncated_ode.dimension)
    throw std::invalid_argument(
        context + ": malformed normalized adjoint defect payload");
  const auto radius_upper = Magnitude::upper_abs(witness_radius);
  const auto order = solution.taylor_complete_max;
  auto defect = Magnitude::zero();

  // g=-beta*t*r.  Sum the exact rational coefficient tails; this is an l1
  // bound on the enlarged finite epsilon/component vector.
  for (const auto& entry : exact_row.entries) {
    const auto& multiplier = entry.multiplier;
    if (multiplier.center_pole_order != 0 ||
        !multiplier.analytic_coefficients.has_value())
      throw std::domain_error(
          context + ": normalized defect needs an ordinary exact row");
    for (std::size_t epsilon = 0;
         epsilon < multiplier.analytic_coefficients->size(); ++epsilon) {
      const auto raw_power = matching_detail::checked_power(
          static_cast<std::int64_t>(multiplier.epsilon_shift) +
              static_cast<std::int64_t>(epsilon),
          "normalized defect row epsilon power");
      if (raw_power > epsilon_complete_max) continue;
      const auto& rational = (*multiplier.analytic_coefficients)[epsilon];
      auto numerator = rational.numerator;
      numerator.insert(numerator.begin(), Rational(0));
      const auto forcing_tail = exact_t_rational_taylor_tail_l1_upper(
          canonical_exact_t_rational(
              {std::move(numerator), rational.denominator}),
          order, radius_upper,
          context + ": forcing component " +
              std::to_string(entry.column) + ": epsilon " +
              std::to_string(raw_power));
      defect += Magnitude::upper_abs(oriented_physical_jacobian) *
                forcing_tail;
    }
  }

  // For p_N=sum_{k=1}^N p_k t^k, only A_j with j>N-k can enter the
  // coefficient defect above N.  Bound each exact rational A tail before
  // multiplying by the corresponding certified prefix coefficient.
  auto radius_power = radius_upper;
  for (std::uint32_t prefix_order = 1; prefix_order <= order;
       ++prefix_order) {
    const auto prefix_upper = vector_infinity_upper_through(
        solution.coefficients[prefix_order - 1],
        epsilon_complete_max);
    if (!prefix_upper.is_zero()) {
      const auto operator_tail =
          normalized_backward_adjoint_c_taylor_tail_l1_upper(
              normalized, witness_radius, epsilon_complete_max,
              order - prefix_order,
              context + ": C tail for prefix order " +
                  std::to_string(prefix_order));
      defect += prefix_upper * radius_power * operator_tail;
    }
    radius_power = radius_power * radius_upper;
  }
  return defect;
}

inline PreparedSparseLocalMultiplierMatrix<ComplexBall>
specialize_exact_backward_adjoint_row_taylor(
    const PreparedSparseLocalMultiplierMatrix<Rational>& exact_row,
    std::uint32_t taylor_complete_max,
    const std::string& context) {
  PreparedSparseLocalMultiplierMatrix<ComplexBall> result;
  result.rows = exact_row.rows;
  result.columns = exact_row.columns;
  result.exact_identity = exact_row.exact_identity;
  result.entries.reserve(exact_row.entries.size());
  for (const auto& entry : exact_row.entries) {
    const auto& source = entry.multiplier;
    if (!source.analytic_coefficients.has_value() ||
        source.analytic_coefficients->size() != source.kernels.size())
      throw std::domain_error(
          context + ": exact row has no complete analytic coefficients");
    PreparedRationalTaylorMultiplier<ComplexBall> target;
    target.epsilon_shift = source.epsilon_shift;
    target.center_pole_order = source.center_pole_order;
    target.exact_identity = source.exact_identity;
    target.proven_zero = source.proven_zero;
    target.kernels.reserve(source.analytic_coefficients->size());
    std::vector<PreparedRationalAnalyticCoefficient<ComplexBall>> analytic;
    analytic.reserve(source.analytic_coefficients->size());
    for (std::size_t epsilon = 0;
         epsilon < source.analytic_coefficients->size(); ++epsilon) {
      const auto& coefficient = (*source.analytic_coefficients)[epsilon];
      const auto exact = canonical_exact_t_rational(
          {coefficient.numerator, coefficient.denominator});
      const auto expanded = expand_exact_t_rational_taylor(
          exact, taylor_complete_max,
          context + ": entry " + std::to_string(entry.column) +
              ": epsilon " + std::to_string(epsilon));
      std::vector<ComplexBall> numeric;
      numeric.reserve(expanded.size());
      for (const auto& value : expanded)
        numeric.push_back(ComplexBall::from_strings(value.str()));
      target.kernels.push_back(std::move(numeric));
      PreparedRationalAnalyticCoefficient<ComplexBall> numeric_exact;
      for (const auto& value : exact.numerator)
        numeric_exact.numerator.push_back(
            ComplexBall::from_strings(value.str()));
      for (const auto& value : exact.denominator)
        numeric_exact.denominator.push_back(
            ComplexBall::from_strings(value.str()));
      analytic.push_back(std::move(numeric_exact));
    }
    target.analytic_coefficients = std::move(analytic);
    result.entries.push_back(
        {entry.row, entry.column, std::move(target)});
  }
  return result;
}

}  // namespace adjoint_observable_detail

inline Magnitude exact_combined_backward_adjoint_forcing_disk_upper(
    const BackwardAdjointTaylorProblem& problem,
    const PreparedSparseLocalMultiplierMatrix<Rational>& exact_row,
    const PreparedPhysicalClearedODE<Rational>& exact_ode,
    const ComplexBall& oriented_physical_jacobian,
    const ComplexBall& witness_radius,
    std::int32_t epsilon_complete_max,
    const std::string& context) {
  using namespace adjoint_observable_detail;
  physical_ode_detail::validate_ode(exact_ode);
  if (exact_row.rows != 1 || exact_row.columns != exact_ode.dimension ||
      exact_row.columns != problem.dimension ||
      problem.q_lags.size() > exact_ode.q_lags.size())
    throw std::invalid_argument(
        context + ": exact q/row forcing payload has the wrong shape");
  const auto radius_upper = Magnitude::upper_abs(witness_radius);
  std::vector<EpsilonFrame<Rational>> exact_q;
  exact_q.reserve(problem.q_lags.size());
  for (std::size_t lag = 0; lag < problem.q_lags.size(); ++lag)
    exact_q.push_back(expand_exact_epsilon_rational(
        exact_ode.q_lags[lag], problem.q_lags[lag].complete_max(),
        context + ": exact q lag " + std::to_string(lag)));

  auto row_disk_upper = Magnitude::zero();
  for (const auto& entry : exact_row.entries) {
    const auto& multiplier = entry.multiplier;
    if (multiplier.center_pole_order != 0 ||
        !multiplier.analytic_coefficients.has_value() ||
        multiplier.analytic_coefficients->size() !=
            multiplier.kernels.size())
      throw std::domain_error(
          context + ": exact combined forcing requires an ordinary complete row");
    std::map<std::int32_t, ExactTRationalFunction> combined;
    auto q_min = std::numeric_limits<std::int32_t>::max();
    auto q_max = std::numeric_limits<std::int32_t>::min();
    for (const auto& q : exact_q) {
      q_min = std::min(q_min, q.min_power());
      q_max = std::max(q_max, q.complete_max());
    }
    for (std::int64_t raw_q = q_min; raw_q <= q_max; ++raw_q) {
      const auto q_power = static_cast<std::int32_t>(raw_q);
      std::vector<Rational> q_polynomial(
          exact_q.size(), Rational(0));
      for (std::size_t lag = 0; lag < exact_q.size(); ++lag)
        if (q_power >= exact_q[lag].min_power() &&
            q_power <= exact_q[lag].complete_max())
          q_polynomial[lag] = exact_q[lag].coefficient(q_power);
      trim_exact_t_polynomial(q_polynomial);
      if (exact_t_polynomial_zero(q_polynomial)) continue;
      for (std::size_t offset = 0;
           offset < multiplier.analytic_coefficients->size(); ++offset) {
        const auto row_power = matching_detail::checked_power(
            static_cast<std::int64_t>(multiplier.epsilon_shift) +
                static_cast<std::int64_t>(offset),
            "exact combined row epsilon power");
        const auto output_power = matching_detail::checked_power(
            raw_q + static_cast<std::int64_t>(row_power),
            "exact combined forcing epsilon power");
        if (output_power > epsilon_complete_max) continue;
        const auto& row_rational =
            (*multiplier.analytic_coefficients)[offset];
        ExactTRationalFunction term{
            multiply_exact_t_polynomials(
                q_polynomial, row_rational.numerator),
            row_rational.denominator};
        term = canonical_exact_t_rational(std::move(term));
        const auto found = combined.find(output_power);
        if (found == combined.end())
          combined.emplace(output_power, std::move(term));
        else
          found->second = add_exact_t_rationals(found->second, term);
      }
    }
    auto component_upper = Magnitude::zero();
    for (const auto& [epsilon_power, rational] : combined)
      component_upper += exact_t_rational_disk_upper(
          rational, radius_upper,
          context + ": component " + std::to_string(entry.column) +
              ": epsilon " + std::to_string(epsilon_power));
    row_disk_upper = Magnitude::maximum(
        row_disk_upper, component_upper);
  }
  return Magnitude::upper_abs(oriented_physical_jacobian) *
         radius_upper * row_disk_upper;
}

template <typename Scalar>
BackwardAdjointTaylorProblem prepare_backward_adjoint_taylor_problem(
    const PreparedPhysicalClearedODE<Scalar>& ode,
    const PreparedSparseLocalMultiplierMatrix<ComplexBall>& row,
    std::uint32_t taylor_complete_max,
    std::int32_t coefficient_epsilon_complete_max,
    std::int32_t required_epsilon_complete_max,
    const ComplexBall& oriented_physical_jacobian,
    const std::string& context = "prepare backward Fuchsian adjoint") {
  using namespace adjoint_observable_detail;
  physical_ode_detail::validate_ode(ode);
  if (taylor_complete_max == 0 || row.rows != 1 ||
      row.columns != ode.dimension || row.entries.empty() ||
      coefficient_epsilon_complete_max < required_epsilon_complete_max)
    throw std::invalid_argument(
        context + ": inconsistent ODE, row, Taylor, or epsilon contract");

  BackwardAdjointTaylorProblem problem;
  problem.dimension = ode.dimension;
  problem.taylor_complete_max = taylor_complete_max;
  problem.required_epsilon_complete_max =
      required_epsilon_complete_max;
  problem.q_lags.reserve(
      std::min<std::size_t>(ode.q_lags.size(), taylor_complete_max + 1));
  for (std::size_t lag = 0;
       lag < ode.q_lags.size() && lag <= taylor_complete_max; ++lag)
    problem.q_lags.push_back(expand_epsilon_rational(
        ode.q_lags[lag], coefficient_epsilon_complete_max,
        context + ": q lag " + std::to_string(lag)));

  const auto c_lag_count = std::min<std::size_t>(
      ode.c_lags.size(), taylor_complete_max + 1);
  problem.c_transpose_lags.reserve(c_lag_count);
  for (std::size_t lag = 0; lag < c_lag_count; ++lag) {
    FiniteLaurentMatrix<ComplexBall> transposed(
        ode.dimension, FiniteLaurentVector<ComplexBall>());
    for (auto& target_row : transposed)
      for (std::uint32_t column = 0; column < ode.dimension; ++column)
        target_row.push_back(
            structural_zero_frame(coefficient_epsilon_complete_max));
    for (const auto& entry : ode.c_lags[lag])
      transposed[entry.column][entry.row] = expand_epsilon_rational(
          entry.value, coefficient_epsilon_complete_max,
          context + ": C lag " + std::to_string(lag));
    problem.c_transpose_lags.push_back(std::move(transposed));
  }

  problem.forcing.assign(
      taylor_complete_max,
      structural_zero_vector(ode.dimension,
                             coefficient_epsilon_complete_max));
  for (const auto& entry : row.entries) {
    if (entry.row != 0 || entry.column >= ode.dimension ||
        entry.multiplier.structurally_zero())
      throw std::invalid_argument(
          context + ": prepared row contains an invalid active entry");
    const auto& multiplier = entry.multiplier;
    if (multiplier.kernels.empty())
      throw std::invalid_argument(
          context + ": prepared row entry has no Taylor kernels");
    const auto taylor_width = multiplier.kernels.front().size();
    if (taylor_width == 0 ||
        std::any_of(multiplier.kernels.begin(), multiplier.kernels.end(),
                    [&](const auto& kernel) {
                      return kernel.size() != taylor_width;
                    }))
      throw std::invalid_argument(
          context + ": prepared row entry has a ragged Taylor kernel");

    for (std::size_t row_taylor = 0; row_taylor < taylor_width;
         ++row_taylor) {
      std::vector<ComplexBall> row_epsilon;
      row_epsilon.reserve(multiplier.kernels.size());
      for (const auto& epsilon_kernel : multiplier.kernels)
        row_epsilon.push_back(epsilon_kernel[row_taylor]);
      EpsilonFrame<ComplexBall> row_coefficient(
          multiplier.epsilon_shift, std::move(row_epsilon));
      if (!material_frame(row_coefficient)) continue;
      const auto row_power =
          static_cast<std::int64_t>(row_taylor) -
          multiplier.center_pole_order;
      for (std::size_t q_lag = 0; q_lag < problem.q_lags.size(); ++q_lag) {
        if (!material_frame(problem.q_lags[q_lag])) continue;
        const auto forcing_power = row_power +
            static_cast<std::int64_t>(q_lag) + 1;
        auto forcing_coefficient =
            problem.q_lags[q_lag] * row_coefficient;
        forcing_coefficient = forcing_coefficient.scaled(
            -oriented_physical_jacobian);
        if (!material_frame(forcing_coefficient)) continue;
        if (forcing_power <= 0)
          throw std::domain_error(
              context +
              ": center-pole forcing reaches t^0 or below; a Laurent/log "
              "backward Fuchsian completion is required");
        if (forcing_power > taylor_complete_max) continue;
        auto& target = problem.forcing[
            static_cast<std::size_t>(forcing_power - 1)][entry.column];
        target = target + forcing_coefficient;
      }
    }
  }
  return problem;
}

inline Magnitude backward_adjoint_forcing_cauchy_numerator_upper(
    const BackwardAdjointTaylorProblem& problem,
    const PreparedSparseLocalMultiplierMatrix<ComplexBall>& row,
    const ComplexBall& oriented_physical_jacobian,
    const ComplexBall& witness_radius,
    std::optional<std::int32_t> epsilon_complete_max,
    const PreparedSparseLocalMultiplierMatrix<Rational>* exact_row,
    const PreparedPhysicalClearedODE<Rational>* exact_ode,
    const std::string& context =
        "backward Fuchsian adjoint forcing Cauchy bound") {
  using namespace adjoint_observable_detail;
  const auto radius_lower = Magnitude::lower_abs(witness_radius);
  const auto radius_upper = Magnitude::upper_abs(witness_radius);
  if (radius_lower.is_zero() || !radius_upper.is_finite() ||
      row.rows != 1 || row.columns != problem.dimension)
    throw std::invalid_argument(
        context + ": invalid witness radius or row dimension");
  if ((exact_row == nullptr) != (exact_ode == nullptr))
    throw std::invalid_argument(
        context + ": exact row and exact physical equation must be paired");
  if (exact_row != nullptr) {
    if (!epsilon_complete_max.has_value())
      throw std::invalid_argument(
          context + ": exact combined forcing needs an epsilon cap");
    return exact_combined_backward_adjoint_forcing_disk_upper(
        problem, *exact_row, *exact_ode, oriented_physical_jacobian,
        witness_radius, *epsilon_complete_max, context);
  }

  auto q_disk_upper = Magnitude::zero();
  auto t_power = Magnitude::one();
  for (const auto& q : problem.q_lags) {
    q_disk_upper += frame_l1_upper(q) * t_power;
    t_power = t_power * radius_upper;
  }

  auto row_disk_upper = Magnitude::zero();
  for (const auto& entry : row.entries) {
    const auto& multiplier = entry.multiplier;
    if (multiplier.center_pole_order != 0)
      throw std::domain_error(
          context +
          ": ordinary composed-tail theorem does not admit a center pole");
    if (!multiplier.analytic_coefficients.has_value() ||
        multiplier.analytic_coefficients->size() !=
            multiplier.kernels.size())
      throw std::domain_error(
          context +
          ": complete analytic row coefficients are unavailable");
    auto component_upper = Magnitude::zero();
    for (std::size_t epsilon = 0;
         epsilon < multiplier.analytic_coefficients->size(); ++epsilon) {
      const auto raw_power = matching_detail::checked_power(
          static_cast<std::int64_t>(multiplier.epsilon_shift) +
              static_cast<std::int64_t>(epsilon),
          "backward adjoint analytic row epsilon power");
      if (epsilon_complete_max.has_value() &&
          raw_power > *epsilon_complete_max)
        continue;
      component_upper += rational_disk_upper(
          (*multiplier.analytic_coefficients)[epsilon], radius_upper,
          context + ": row component " +
              std::to_string(entry.column) + ": epsilon offset " +
              std::to_string(epsilon));
    }
    row_disk_upper = Magnitude::maximum(
        row_disk_upper, component_upper);
  }
  if (row_disk_upper.is_zero()) return Magnitude::zero();
  return Magnitude::upper_abs(oriented_physical_jacobian) *
         radius_upper * q_disk_upper * row_disk_upper;
}

inline BackwardAdjointTailCertificate
certify_backward_adjoint_taylor_tail(
    const BackwardAdjointTaylorProblem& problem,
    const BackwardAdjointTaylorResult& solution,
    const ComplexBall& evaluation_point,
    const ComplexBall& witness_radius,
    const Magnitude& forcing_cauchy_numerator_upper,
    const std::string& context =
        "backward Fuchsian adjoint Taylor-tail certificate",
    std::optional<Magnitude> analytic_c_tail_sum_upper = std::nullopt,
    std::optional<Magnitude> a_posteriori_defect_l1_upper = std::nullopt) {
  using namespace adjoint_observable_detail;
  if (solution.dimension != problem.dimension ||
      solution.taylor_complete_max != problem.taylor_complete_max ||
      solution.coefficients.size() != problem.taylor_complete_max ||
      problem.q_lags.empty() || problem.c_transpose_lags.empty())
    throw std::invalid_argument(
        context + ": problem and Taylor solution disagree");
  if (solution.max_log_power != 0)
    throw std::domain_error(
        context +
        ": logarithmic composed adjoint requires the augmented log-stack tail theorem");
  const auto evaluation_upper = Magnitude::upper_abs(evaluation_point);
  const auto radius_lower = Magnitude::lower_abs(witness_radius);
  const auto radius_upper = Magnitude::upper_abs(witness_radius);
  if (radius_lower.is_zero() || !radius_upper.is_finite())
    throw std::invalid_argument(context + ": invalid witness radius");
  const auto evaluation_ratio = evaluation_upper / radius_lower;
  if (Magnitude::one() <= evaluation_ratio)
    throw std::domain_error(
        context + ": evaluation point is not strictly inside witness disk");

  const auto& q0 = problem.q_lags.front();
  if (q0.min_power() != 0 || q0.coefficient(0).contains_zero())
    throw std::domain_error(
        context + ": q(0,epsilon) is not a certified formal unit");
  for (const auto& q : problem.q_lags)
    if (q.min_power() < 0)
      throw std::domain_error(
          context +
          ": negative-epsilon q coupling is not closed in the finite causal stack");
  for (const auto& matrix : problem.c_transpose_lags)
    for (const auto& row : matrix)
      for (const auto& entry : row)
        if (entry.min_power() < 0)
          throw std::domain_error(
              context +
              ": negative-epsilon C coupling is not closed in the finite causal stack");
  std::vector<ComplexBall> one_coefficients(q0.window().width(),
                                            ComplexBall(0));
  one_coefficients.front() = ComplexBall(1);
  const auto q0_inverse = finite_laurent_quotient(
      EpsilonFrame<ComplexBall>(q0.window(), std::move(one_coefficients)),
      q0, context + ": q0 inverse");
  const auto q_inverse_upper = frame_l1_upper(q0_inverse);
  const auto c0_upper =
      matrix_infinity_upper(problem.c_transpose_lags.front());
  const auto indicial_upper = q_inverse_upper * c0_upper;

  const auto next_order =
      static_cast<ulong>(problem.taylor_complete_max) + 1;
  const auto next_order_magnitude = Magnitude::from_ui(next_order);
  const auto inverse_denominator_lower =
      Magnitude::positive_difference_lower(
          next_order_magnitude, indicial_upper);
  if (inverse_denominator_lower.is_zero())
    throw std::domain_error(
        context +
        ": retained Taylor order does not dominate the indicial norm");

  auto q_tail_sum = Magnitude::zero();
  for (std::size_t lag = 1; lag < problem.q_lags.size(); ++lag)
    q_tail_sum += frame_l1_upper(problem.q_lags[lag]) *
                  radius_upper.power_upper(static_cast<ulong>(lag));
  auto c_tail_sum = analytic_c_tail_sum_upper.value_or(Magnitude::zero());
  if (!analytic_c_tail_sum_upper.has_value())
    for (std::size_t lag = 1;
         lag < problem.c_transpose_lags.size(); ++lag)
      c_tail_sum += matrix_infinity_upper(
                        problem.c_transpose_lags[lag]) *
                    radius_upper.power_upper(static_cast<ulong>(lag));
  const auto recurrence_contraction =
      q_inverse_upper *
      (next_order_magnitude * q_tail_sum + c_tail_sum) /
      inverse_denominator_lower;
  if (Magnitude::one() <= recurrence_contraction)
    throw std::domain_error(
        context +
        ": recurrence majorant is not contractive on the witness disk");
  const auto contraction_gap_lower =
      Magnitude::positive_difference_lower(
          Magnitude::one(), recurrence_contraction);
  if (contraction_gap_lower.is_zero())
    throw std::domain_error(
        context + ": recurrence contraction has no positive gap");

  auto coefficient_majorant = Magnitude::zero();
  if (a_posteriori_defect_l1_upper.has_value()) {
    // Weighted-l1 coefficient theorem for the normalized equation.  If
    // L_n=nI+A0 and M bounds every ||L_n^-1|| for n>N, then the unseen
    // coefficient l1 norm on |t|=R is at most M*D/(1-M*S).  D is the exact
    // high-order defect of the retained prefix and S is the positive-lag
    // operator l1 norm already used in recurrence_contraction.
    coefficient_majorant =
        q_inverse_upper * *a_posteriori_defect_l1_upper /
        inverse_denominator_lower / contraction_gap_lower;
  } else {
    coefficient_majorant =
        q_inverse_upper * forcing_cauchy_numerator_upper /
        inverse_denominator_lower / contraction_gap_lower;
    const auto lag_memory = std::max(
        problem.q_lags.size(), problem.c_transpose_lags.size());
    const auto first_known = problem.taylor_complete_max + 1 > lag_memory
        ? problem.taylor_complete_max + 1 - lag_memory
        : 1;
    for (std::uint32_t order = first_known;
         order <= problem.taylor_complete_max; ++order) {
      const auto known = vector_infinity_upper(
          solution.coefficients.at(order - 1)) *
          radius_upper.power_upper(order);
      coefficient_majorant = Magnitude::maximum(
          coefficient_majorant, known);
    }
  }

  const auto ratio_gap_lower = Magnitude::positive_difference_lower(
      Magnitude::one(), evaluation_ratio);
  if (ratio_gap_lower.is_zero())
    throw std::domain_error(
        context + ": evaluation/witness ratio has no positive gap");
  const auto tail = a_posteriori_defect_l1_upper.has_value()
      ? coefficient_majorant * evaluation_ratio.power_upper(next_order)
      : coefficient_majorant * evaluation_ratio.power_upper(next_order) /
            ratio_gap_lower;
  return {tail, coefficient_majorant, recurrence_contraction,
          evaluation_ratio, problem.taylor_complete_max};
}

inline BackwardAdjointAdaptiveTailCertificate
certify_backward_adjoint_taylor_tail_adaptive_witness(
    const BackwardAdjointTaylorProblem& problem,
    const BackwardAdjointTaylorResult& solution,
    const PreparedSparseLocalMultiplierMatrix<ComplexBall>& row,
    const ComplexBall& oriented_physical_jacobian,
    const ComplexBall& evaluation_point,
    const Rational& evaluation_modulus_exact,
    const Rational& chart_radius_exact,
    std::uint32_t maximum_attempts = 16,
    const std::string& context =
        "backward Fuchsian adjoint adaptive Taylor-tail certificate",
    std::optional<std::int32_t> epsilon_complete_max = std::nullopt,
    const PreparedSparseLocalMultiplierMatrix<Rational>* exact_row = nullptr,
    const PreparedPhysicalClearedODE<Rational>* exact_ode = nullptr,
    const adjoint_observable_detail::NormalizedBackwardAdjointExactODE*
        normalized_exact_ode = nullptr) {
  if (evaluation_modulus_exact.sign() < 0 ||
      !(evaluation_modulus_exact < chart_radius_exact) ||
      maximum_attempts == 0)
    throw std::invalid_argument(
        context + ": invalid exact evaluation/chart witness interval");
  const auto gap = chart_radius_exact - evaluation_modulus_exact;
  std::string last_rejection;
  Rational denominator(2);
  for (std::uint32_t attempt = 0; attempt < maximum_attempts; ++attempt) {
    const auto witness_exact =
        evaluation_modulus_exact + gap / denominator;
    const auto witness = ComplexBall::from_strings(witness_exact.str());
    try {
      const auto forcing =
          backward_adjoint_forcing_cauchy_numerator_upper(
              problem, row, oriented_physical_jacobian, witness,
              epsilon_complete_max, exact_row, exact_ode,
              context + ": forcing attempt " +
                  std::to_string(attempt + 1));
      std::optional<Magnitude> normalized_c_tail;
      std::optional<Magnitude> normalized_defect;
      if (normalized_exact_ode != nullptr) {
        if (!epsilon_complete_max.has_value())
          throw std::invalid_argument(
              context + ": normalized C tail needs an epsilon cap");
        normalized_c_tail =
            adjoint_observable_detail::
                normalized_backward_adjoint_c_tail_disk_upper(
                    *normalized_exact_ode, witness,
                    *epsilon_complete_max,
                    context + ": normalized C tail attempt " +
                        std::to_string(attempt + 1));
        if (exact_row == nullptr)
          throw std::invalid_argument(
              context + ": normalized defect needs the exact row");
        normalized_defect =
            adjoint_observable_detail::
                normalized_backward_adjoint_defect_l1_upper(
                    *normalized_exact_ode, *exact_row, solution,
                    oriented_physical_jacobian, witness,
                    *epsilon_complete_max,
                    context + ": normalized defect attempt " +
                        std::to_string(attempt + 1));
      }
      auto tail = certify_backward_adjoint_taylor_tail(
          problem, solution, evaluation_point, witness, forcing,
          context + ": tail attempt " + std::to_string(attempt + 1),
          normalized_c_tail, normalized_defect);
      return {std::move(tail), witness_exact};
    } catch (const std::domain_error& error) {
      last_rejection = error.what();
    }
    denominator *= Rational(2);
  }
  throw std::domain_error(
      context + ": no certified witness radius in " +
      std::to_string(maximum_attempts) +
      " deterministic dyadic attempts; last rejection=" +
      last_rejection);
}

inline BackwardAdjointTaylorResult solve_backward_adjoint_taylor(
    const BackwardAdjointTaylorProblem& problem,
    const PreparedPhysicalClearedODE<Rational>* exact_physical_ode,
    const std::string& context) {
  using namespace adjoint_observable_detail;
  if (problem.dimension == 0 || problem.taylor_complete_max == 0 ||
      problem.q_lags.empty() || problem.c_transpose_lags.empty())
    throw std::invalid_argument(
        context + ": incomplete dimension, Taylor order, q, or C payload");
  if (problem.forcing.size() > problem.taylor_complete_max)
    throw std::invalid_argument(
        context + ": forcing extends beyond the requested Taylor order");
  for (const auto& matrix : problem.c_transpose_lags)
    require_matrix_shape(matrix, problem.dimension,
                         "backward adjoint C lag");
  for (const auto& rhs : problem.forcing)
    if (rhs.size() != problem.dimension)
      throw std::invalid_argument(
          context + ": forcing vector has the wrong dimension");
  if (exact_physical_ode != nullptr) {
    physical_ode_detail::validate_ode(*exact_physical_ode);
    if (exact_physical_ode->dimension != problem.dimension)
      throw std::invalid_argument(
          context + ": exact Rational shadow has the wrong dimension");
  }
  auto structural_complete_max = problem.q_lags.front().complete_max();
  for (const auto& q : problem.q_lags)
    structural_complete_max = std::min(
        structural_complete_max, q.complete_max());
  for (const auto& matrix : problem.c_transpose_lags)
    for (const auto& row : matrix)
      for (const auto& entry : row)
        structural_complete_max = std::min(
            structural_complete_max, entry.complete_max());
  if (structural_complete_max < problem.required_epsilon_complete_max)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InsufficientCompleteWindow,
        context + ": operator does not cover the required epsilon output",
        std::nullopt, std::nullopt, structural_complete_max);

  BackwardAdjointTaylorResult result;
  result.dimension = problem.dimension;
  result.taylor_complete_max = problem.taylor_complete_max;
  result.common_epsilon_complete_max =
      std::numeric_limits<std::int32_t>::max();
  result.coefficients.reserve(problem.taylor_complete_max);
  result.higher_log_coefficients.reserve(problem.taylor_complete_max);

  for (std::uint32_t n = 1; n <= problem.taylor_complete_max; ++n) {
    const auto rhs_log_max = result.max_log_power;
    std::vector<FiniteLaurentVector<ComplexBall>> rhs_logs;
    rhs_logs.reserve(static_cast<std::size_t>(rhs_log_max) + 1);
    rhs_logs.push_back(n <= problem.forcing.size()
                           ? problem.forcing[n - 1]
                           : structural_zero_vector(
                                 problem.dimension,
                                 structural_complete_max));
    for (std::uint32_t log = 1; log <= rhs_log_max; ++log)
      rhs_logs.push_back(structural_zero_vector(
          problem.dimension, structural_complete_max));

    const auto previous_log = [&](std::uint32_t order,
                                  std::uint32_t log)
        -> const FiniteLaurentVector<ComplexBall>* {
      if (log == 0) return &result.coefficients.at(order - 1);
      const auto& higher = result.higher_log_coefficients.at(order - 1);
      return log - 1 < higher.size() ? &higher[log - 1] : nullptr;
    };

    const auto maximum_lag = std::min<std::size_t>(
        n - 1,
        std::max(problem.q_lags.size(),
                 problem.c_transpose_lags.size()) - 1);
    for (std::size_t lag = 1; lag <= maximum_lag; ++lag) {
      const auto previous_n = n - static_cast<std::uint32_t>(lag);
      const auto& q = lag < problem.q_lags.size()
          ? problem.q_lags[lag]
          : structural_zero_frame(structural_complete_max);
      const auto* c = lag < problem.c_transpose_lags.size()
          ? &problem.c_transpose_lags[lag]
          : nullptr;
      for (std::uint32_t log = 0; log <= rhs_log_max; ++log) {
        const auto* current = previous_log(previous_n, log);
        if (current == nullptr) continue;
        const auto* above = previous_log(previous_n, log + 1);
        const auto contribution = apply_log_lag_operator(
            q, c, previous_n, *current, above,
            structural_complete_max);
        subtract_vector_in_place(rhs_logs[log], contribution);
      }
    }

    FiniteLaurentMatrix<ComplexBall> layer(
        problem.dimension, FiniteLaurentVector<ComplexBall>());
    for (std::uint32_t row = 0; row < problem.dimension; ++row) {
      layer[row].reserve(problem.dimension);
      for (std::uint32_t column = 0; column < problem.dimension; ++column) {
        auto value = problem.c_transpose_lags.front()[row][column];
        if (row == column)
          value = value + scaled_integer(problem.q_lags.front(), n);
        layer[row].push_back(std::move(value));
      }
    }
    std::optional<FiniteLaurentMatrix<Rational>> exact_layer;
    std::optional<ExactSimpleResonanceFrame> resonance;
    if (exact_physical_ode != nullptr) {
      auto exact_complete_max = layer.front().front().complete_max();
      for (const auto& row : layer)
        for (const auto& entry : row)
          exact_complete_max = std::min(
              exact_complete_max, entry.complete_max());
      exact_layer = exact_backward_adjoint_layer(
          *exact_physical_ode, n, exact_complete_max,
          context + ": Taylor layer " + std::to_string(n));
      apply_exact_backward_adjoint_layer_zeros(
          layer, *exact_layer,
          context + ": Taylor layer " + std::to_string(n) +
              ": Rational-shadow structural projection");
      resonance = exact_simple_resonance_frame(
          *exact_layer,
          context + ": Taylor layer " + std::to_string(n));
    }

    std::vector<FiniteLaurentVector<ComplexBall>> solved_logs;
    if (!resonance.has_value()) {
      try {
        auto factorization =
            factor_preconditioned_acb_finite_laurent_system(
                layer, context + ": Taylor layer " + std::to_string(n));
        solved_logs.resize(rhs_logs.size());
        for (std::size_t reverse = rhs_logs.size(); reverse-- > 0;) {
          auto adjusted = rhs_logs[reverse];
          if (reverse + 1 < solved_logs.size()) {
            const auto shifted =
                epsilon_shift_vector(solved_logs[reverse + 1]);
            const auto derivative = apply_lag_operator(
                problem.q_lags.front(), nullptr, 1, shifted,
                structural_complete_max);
            subtract_vector_in_place(adjusted, derivative);
          }
          solved_logs[reverse] = solve_factorized_finite_laurent_system(
              factorization, std::move(adjusted),
              context + ": Taylor coefficient " + std::to_string(n) +
                  ", log " + std::to_string(reverse));
        }
      } catch (const MatchingArithmeticError& error) {
        throw MatchingArithmeticError(
            error.code,
            context + ": Taylor layer " + std::to_string(n) +
                " is singular or epsilon-incomplete; a logarithmic/resonant "
                "Fuchsian completion is required; detail=" + error.what(),
            error.row, error.column, error.epsilon_power);
      }
    } else {
      const auto right_kernel =
          exact_vector_to_ball(resonance->right_kernel);
      const auto left_kernel =
          exact_vector_to_ball(resonance->left_kernel);
      if (right_kernel.size() != problem.dimension ||
          left_kernel.size() != problem.dimension)
        throw std::logic_error(
            context + ": exact resonant frame changed dimension");

      FiniteLaurentMatrix<ComplexBall> bordered;
      bordered.reserve(static_cast<std::size_t>(problem.dimension) + 1);
      for (std::uint32_t row = 0; row < problem.dimension; ++row) {
        auto bordered_row = layer[row];
        bordered_row.push_back(
            (problem.q_lags.front() * right_kernel[row]).shifted(1));
        bordered.push_back(std::move(bordered_row));
      }
      FiniteLaurentVector<ComplexBall> gauge = left_kernel;
      gauge.push_back(structural_zero_frame(
          structural_complete_max));
      bordered.push_back(std::move(gauge));

      const auto exact_complete_max = exact_layer->front().front()
                                              .complete_max();
      const auto exact_q0 = expand_exact_epsilon_rational(
          exact_physical_ode->q_lags.front(), exact_complete_max,
          context + ": exact resonant q0");
      FiniteLaurentMatrix<Rational> exact_bordered;
      exact_bordered.reserve(
          static_cast<std::size_t>(problem.dimension) + 1);
      for (std::uint32_t row = 0; row < problem.dimension; ++row) {
        auto bordered_row = (*exact_layer)[row];
        bordered_row.push_back(
            (exact_q0 * resonance->right_kernel[row]).shifted(1));
        exact_bordered.push_back(std::move(bordered_row));
      }
      auto exact_gauge = resonance->left_kernel;
      exact_gauge.push_back(EpsilonFrame<Rational>::zero(
          exact_complete_max));
      exact_bordered.push_back(std::move(exact_gauge));
      apply_exact_backward_adjoint_layer_zeros(
          bordered, exact_bordered,
          context + ": Taylor layer " + std::to_string(n) +
              ": exact resonant bordered projection");
      auto factorization =
          factor_preconditioned_acb_finite_laurent_system(
              bordered, context + ": Taylor layer " +
                            std::to_string(n) +
                            ": simple resonant bordered solve");

      solved_logs.assign(
          rhs_logs.size() + 1,
          structural_zero_vector(
              problem.dimension, structural_complete_max));
      for (std::size_t log = 0; log < rhs_logs.size(); ++log) {
        auto adjusted = rhs_logs[log];
        subtract_vector_in_place(
            adjusted, apply_finite_laurent_matrix(layer, solved_logs[log]));
        auto bordered_rhs = std::move(adjusted);
        bordered_rhs.push_back(structural_zero_frame(
            structural_complete_max));
        auto solved = solve_factorized_finite_laurent_system(
            factorization, std::move(bordered_rhs),
            context + ": Taylor coefficient " + std::to_string(n) +
                ", resonant log " + std::to_string(log));
        FiniteLaurentVector<ComplexBall> correction;
        correction.reserve(problem.dimension);
        for (std::uint32_t component = 0;
             component < problem.dimension; ++component)
          correction.push_back(std::move(solved[component]));
        add_vector_in_place(solved_logs[log], correction);
        const auto obstruction = std::move(solved.back());
        FiniteLaurentVector<ComplexBall> raised;
        raised.reserve(problem.dimension);
        for (const auto& component : right_kernel)
          raised.push_back(component * obstruction);
        add_vector_in_place(solved_logs[log + 1], raised);
      }
      const auto top_residual = apply_finite_laurent_matrix(
          layer, solved_logs.back());
      for (std::size_t component = 0;
           component < top_residual.size(); ++component)
        for (std::int64_t raw_power = top_residual[component].min_power();
             raw_power <= top_residual[component].complete_max();
             ++raw_power)
          if (!top_residual[component]
                   .coefficient(static_cast<std::int32_t>(raw_power))
                   .contains_zero())
            throw MatchingArithmeticError(
                MatchingArithmeticErrorCode::SaturationFailure,
                context + ": resonant top-log residual excludes zero",
                component, std::nullopt,
                static_cast<std::int32_t>(raw_power));
    }

    while (solved_logs.size() > 1) {
      const auto& last = solved_logs.back();
      bool material = false;
      for (const auto& component : last)
        material = material || material_frame(component);
      if (material) break;
      solved_logs.pop_back();
    }
    for (std::size_t log_power = 0; log_power < solved_logs.size();
         ++log_power)
      for (const auto& component : solved_logs[log_power]) {
        if (!material_frame(component)) continue;
        const auto evaluated_complete = matching_detail::checked_power(
            static_cast<std::int64_t>(component.complete_max()) +
                static_cast<std::int64_t>(log_power),
            "backward adjoint evaluated log completeness");
        result.common_epsilon_complete_max = std::min(
            result.common_epsilon_complete_max, evaluated_complete);
        if (evaluated_complete <
            problem.required_epsilon_complete_max)
          throw MatchingArithmeticError(
              MatchingArithmeticErrorCode::InsufficientCompleteWindow,
              context + ": Taylor layer " + std::to_string(n) +
                  " does not cover the required epsilon output",
              std::nullopt, std::nullopt,
              evaluated_complete);
      }
    result.max_log_power = std::max<std::uint32_t>(
        result.max_log_power,
        static_cast<std::uint32_t>(solved_logs.size() - 1));
    result.coefficients.push_back(std::move(solved_logs.front()));
    solved_logs.erase(solved_logs.begin());
    result.higher_log_coefficients.push_back(std::move(solved_logs));
  }

  if (result.coefficients.size() >= 2) {
    const auto previous = vector_upper_norm(
        result.coefficients[result.coefficients.size() - 2]);
    const auto last = vector_upper_norm(result.coefficients.back());
    if (previous > 0.0 && std::isfinite(previous) &&
        std::isfinite(last))
      result.last_coefficient_ratio_upper = last / previous;
  }
  return result;
}

inline BackwardAdjointTaylorResult solve_backward_adjoint_taylor(
    const BackwardAdjointTaylorProblem& problem,
    const std::string& context = "backward Fuchsian adjoint") {
  return solve_backward_adjoint_taylor(problem, nullptr, context);
}

struct BackwardAdjointReservoirSolve {
  BackwardAdjointTaylorProblem problem;
  BackwardAdjointTaylorResult result;
  std::int32_t input_epsilon_complete_max = 0;
};

// A Laurent matrix inverse can consume high epsilon coefficients even when
// only a low output prefix is requested.  Grow that private input reservoir
// inside the recurrence instead of forcing the public observable to request
// extra epsilon orders.  The factory must return the same exact problem
// truncated honestly through `input_complete_max`; `maximum_complete_max` is
// the last coefficient actually retained by the prepared row.
template <typename ProblemFactory>
BackwardAdjointReservoirSolve
solve_backward_adjoint_taylor_with_epsilon_reservoir(
    ProblemFactory&& factory,
    std::int32_t initial_complete_max,
    std::int32_t maximum_complete_max,
    const PreparedPhysicalClearedODE<Rational>* exact_physical_ode,
    const std::string& context =
        "backward Fuchsian adjoint private epsilon reservoir") {
  if (initial_complete_max > maximum_complete_max)
    throw std::invalid_argument(
        context + ": initial epsilon reservoir exceeds its honest maximum");
  auto input_complete_max = initial_complete_max;
  while (true) {
    auto problem = factory(input_complete_max);
    try {
      auto result = solve_backward_adjoint_taylor(
          problem, exact_physical_ode,
          context + ": input_complete_max=" +
              std::to_string(input_complete_max));
      return {std::move(problem), std::move(result), input_complete_max};
    } catch (const MatchingArithmeticError& error) {
      const bool retryable_window_failure =
          error.code ==
              MatchingArithmeticErrorCode::InsufficientCompleteWindow ||
          error.code == MatchingArithmeticErrorCode::AmbiguousZero ||
          error.code == MatchingArithmeticErrorCode::ZeroDivisor ||
          error.code ==
              MatchingArithmeticErrorCode::SingularOrIncompleteSystem ||
          error.code ==
              MatchingArithmeticErrorCode::UnresolvedDeterminantTail;
      if (!retryable_window_failure ||
          input_complete_max >= maximum_complete_max)
        throw;
      std::int64_t additional = 1;
      if (error.epsilon_power.has_value() &&
          *error.epsilon_power < problem.required_epsilon_complete_max)
        additional = static_cast<std::int64_t>(
            problem.required_epsilon_complete_max) -
            *error.epsilon_power;
      const auto next = std::min<std::int64_t>(
          maximum_complete_max,
          static_cast<std::int64_t>(input_complete_max) +
              std::max<std::int64_t>(1, additional));
      if (next <= input_complete_max) throw;
      input_complete_max = matching_detail::checked_power(
          next, "backward adjoint private epsilon reservoir maximum");
    }
  }
}

inline FiniteLaurentVector<ComplexBall> evaluate_backward_adjoint_taylor(
    const BackwardAdjointTaylorResult& solution,
    const ComplexBall& point,
    const ComplexBall& logarithm,
    const std::string& context) {
  using namespace adjoint_observable_detail;
  if (solution.dimension == 0 || solution.coefficients.empty() ||
      solution.coefficients.size() != solution.taylor_complete_max ||
      solution.higher_log_coefficients.size() !=
          solution.taylor_complete_max)
    throw std::invalid_argument(context + ": malformed Taylor solution");
  auto value = structural_zero_vector(
      solution.dimension, solution.common_epsilon_complete_max);
  auto log_scale = ComplexBall(1);
  for (std::uint32_t log = 0; log <= solution.max_log_power; ++log) {
    auto polynomial = structural_zero_vector(
        solution.dimension, solution.common_epsilon_complete_max);
    for (std::size_t reverse = solution.coefficients.size(); reverse-- > 0;) {
      const FiniteLaurentVector<ComplexBall>* coefficient = nullptr;
      if (log == 0) {
        coefficient = &solution.coefficients[reverse];
      } else {
        const auto& higher = solution.higher_log_coefficients[reverse];
        if (log - 1 < higher.size()) coefficient = &higher[log - 1];
      }
      for (std::uint32_t component = 0;
           component < solution.dimension; ++component) {
        polynomial[component] = polynomial[component].scaled(point);
        if (coefficient != nullptr)
          polynomial[component] = polynomial[component] +
                                  (*coefficient)[component];
      }
    }
    for (std::uint32_t component = 0;
         component < solution.dimension; ++component) {
      auto term = polynomial[component].scaled(point);
      if (log > 0) term = term.shifted(static_cast<std::int32_t>(log));
      term = term.scaled(log_scale);
      value[component] = value[component] + term;
    }
    if (log < solution.max_log_power)
      log_scale = local_detail::cb_div_ui(
          log_scale * logarithm, log + 1);
  }
  return value;
}

inline FiniteLaurentVector<ComplexBall> evaluate_backward_adjoint_taylor(
    const BackwardAdjointTaylorResult& solution,
    const ComplexBall& point,
    const std::string& context = "backward Fuchsian adjoint evaluation") {
  return evaluate_backward_adjoint_taylor(
      solution, point, local_detail::cb_log(point), context);
}

inline EpsilonFrame<ComplexBall> contract_backward_adjoint(
    const FiniteLaurentVector<ComplexBall>& adjoint,
    const FiniteLaurentVector<ComplexBall>& physical_value,
    const std::string& context = "backward Fuchsian adjoint contraction") {
  if (adjoint.empty() || adjoint.size() != physical_value.size())
    throw std::invalid_argument(context + ": vector dimensions disagree");
  std::optional<EpsilonFrame<ComplexBall>> result;
  for (std::size_t component = 0; component < adjoint.size(); ++component) {
    auto term = adjoint[component] * physical_value[component];
    result = result.has_value() ? *result + term : std::move(term);
  }
  if (!result.has_value())
    throw std::logic_error(context + ": empty contraction");
  return *result;
}

inline std::vector<Magnitude>
backward_adjoint_contracted_tail_by_output(
    const Magnitude& adjoint_coefficient_tail_upper,
    const FiniteLaurentVector<ComplexBall>& adjoint_shape,
    const FiniteLaurentVector<ComplexBall>& physical_value,
    EpsilonWindow output,
    const std::string& context =
        "backward Fuchsian adjoint coefficientwise tail contraction") {
  if (adjoint_shape.empty() ||
      adjoint_shape.size() != physical_value.size())
    throw std::invalid_argument(context + ": vector dimensions disagree");
  std::vector<Magnitude> result(
      output.width(), Magnitude::zero());
  for (std::int64_t raw_output = output.min_power;
       raw_output <= output.complete_max; ++raw_output) {
    auto incoming_l1 = Magnitude::zero();
    for (std::size_t component = 0;
         component < adjoint_shape.size(); ++component) {
      const auto& adjoint_component = adjoint_shape[component];
      const auto& incoming_component = physical_value[component];
      for (std::int64_t raw_incoming =
               incoming_component.min_power();
           raw_incoming <= incoming_component.complete_max();
           ++raw_incoming) {
        const auto raw_adjoint = raw_output - raw_incoming;
        if (raw_adjoint < adjoint_component.min_power() ||
            raw_adjoint > adjoint_component.complete_max())
          continue;
        incoming_l1 += Magnitude::upper_abs(
            incoming_component.coefficient(
                static_cast<std::int32_t>(raw_incoming)));
      }
    }
    result[static_cast<std::size_t>(
        raw_output - output.min_power)] =
        adjoint_coefficient_tail_upper * incoming_l1;
  }
  return result;
}

}  // namespace diffexp2
