#pragma once

#include "diffexp2/integration.hpp"

#include <algorithm>
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

enum class MatchingArithmeticErrorCode : std::uint8_t {
  DimensionMismatch,
  AmbiguousZero,
  ZeroDivisor,
  SingularOrIncompleteSystem,
  StructurallySingularTransformation,
  ExponentOverflow
};

class MatchingArithmeticError : public std::runtime_error {
 public:
  MatchingArithmeticError(MatchingArithmeticErrorCode error_code,
                          std::string detail,
                          std::optional<std::size_t> error_row = std::nullopt,
                          std::optional<std::size_t> error_column = std::nullopt,
                          std::optional<std::int32_t> error_power = std::nullopt)
      : std::runtime_error(std::move(detail)),
        code(error_code),
        row(error_row),
        column(error_column),
        epsilon_power(error_power) {}

  MatchingArithmeticErrorCode code;
  std::optional<std::size_t> row;
  std::optional<std::size_t> column;
  std::optional<std::int32_t> epsilon_power;
};

namespace matching_detail {

inline std::int32_t checked_power(std::int64_t value,
                                  const char* description) {
  if (value < std::numeric_limits<std::int32_t>::min() ||
      value > std::numeric_limits<std::int32_t>::max())
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::ExponentOverflow,
        std::string(description) + " exceeds the int32 epsilon lattice");
  return static_cast<std::int32_t>(value);
}

enum class ZeroDecision : std::uint8_t { Zero, Nonzero, Ambiguous };

template <typename Scalar>
ZeroDecision zero_decision(const Scalar& value) {
  return ScalarTraits<Scalar>::is_zero(value) ? ZeroDecision::Zero
                                               : ZeroDecision::Nonzero;
}

template <>
inline ZeroDecision zero_decision<ComplexBall>(const ComplexBall& value) {
  if (value.is_zero()) return ZeroDecision::Zero;
  return value.contains_zero() ? ZeroDecision::Ambiguous
                               : ZeroDecision::Nonzero;
}

template <typename Scalar>
std::optional<EpsilonFrame<Scalar>> trim_leading_certified_zeros(
    const EpsilonFrame<Scalar>& input, const std::string& context,
    std::optional<std::size_t> row = std::nullopt,
    std::optional<std::size_t> column = std::nullopt) {
  const auto& coefficients = input.coefficients();
  for (std::size_t i = 0; i < coefficients.size(); ++i) {
    const auto power = checked_power(
        static_cast<std::int64_t>(input.min_power()) +
            static_cast<std::int64_t>(i),
        "finite Laurent leading power");
    switch (zero_decision(coefficients[i])) {
      case ZeroDecision::Zero:
        continue;
      case ZeroDecision::Ambiguous:
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::AmbiguousZero,
            context +
                ": Acb enclosure overlaps zero at a required Laurent "
                "leading-coefficient decision",
            row, column, power);
      case ZeroDecision::Nonzero: {
        if (i == 0) return input;
        return EpsilonFrame<Scalar>(
            {power, input.complete_max()},
            std::vector<Scalar>(coefficients.begin() +
                                    static_cast<std::ptrdiff_t>(i),
                                coefficients.end()));
      }
    }
  }
  return std::nullopt;
}

template <typename Scalar>
EpsilonFrame<Scalar> canonical_leading_frame(
    const EpsilonFrame<Scalar>& input, const std::string& context,
    std::optional<std::size_t> row = std::nullopt,
    std::optional<std::size_t> column = std::nullopt) {
  auto trimmed = trim_leading_certified_zeros(input, context, row, column);
  if (trimmed.has_value()) return std::move(*trimmed);
  // Every coefficient through CompleteMax was proved exactly zero.  Moving
  // MinPower to CompleteMax records precisely that fact; coefficients above
  // CompleteMax remain unknown, as required by EpsilonFrame.
  return EpsilonFrame<Scalar>::zero(input.complete_max());
}

template <typename Matrix>
std::size_t rectangular_columns(const Matrix& matrix,
                                const char* description) {
  if (matrix.empty())
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        std::string(description) + " cannot be empty");
  const auto columns = matrix.front().size();
  if (columns == 0 ||
      !std::all_of(matrix.begin(), matrix.end(),
                   [columns](const auto& row) {
                     return row.size() == columns;
                   }))
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        std::string(description) + " must be a nonempty rectangular matrix");
  return columns;
}

}  // namespace matching_detail

// Return the first coefficient which is proved nonzero.  An exact zero is
// skipped, while an Acb enclosure containing zero is deliberately not
// guessed away.  A missing result means that the whole finite complete
// window is exactly zero; it does not claim anything above CompleteMax.
template <typename Scalar>
std::optional<std::int32_t> finite_laurent_leading_power(
    const EpsilonFrame<Scalar>& input, const std::string& context =
                                              "finite Laurent valuation") {
  const auto trimmed =
      matching_detail::trim_leading_certified_zeros(input, context);
  if (!trimmed.has_value()) return std::nullopt;
  return trimmed->min_power();
}

// Solve denominator * quotient == numerator coefficient by coefficient.
// If the denominator starts L rows above its declared MinPower, the exact
// zero rows are trimmed first; the output then has
//
//   MinPower = numerator.MinPower - denominator.actual_leading_power
//   width    = min(numerator.width, trimmed_denominator.width).
//
// This is the finite honest analogue of Laurent-series division.  It never
// samples epsilon or a symbolic regulator.
template <typename Scalar>
EpsilonFrame<Scalar> finite_laurent_quotient(
    const EpsilonFrame<Scalar>& numerator,
    const EpsilonFrame<Scalar>& denominator,
    const std::string& context = "finite Laurent quotient") {
  const auto divisor = matching_detail::trim_leading_certified_zeros(
      denominator, context + ": denominator");
  if (!divisor.has_value())
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::ZeroDivisor,
        context +
            ": denominator is zero throughout its complete epsilon window");

  const auto width =
      std::min(numerator.coefficients().size(),
               divisor->coefficients().size());
  const auto minimum = matching_detail::checked_power(
      static_cast<std::int64_t>(numerator.min_power()) -
          divisor->min_power(),
      "finite Laurent quotient minimum");
  std::vector<Scalar> quotient(width, ScalarTraits<Scalar>::zero());
  const auto& a = numerator.coefficients();
  const auto& b = divisor->coefficients();
  for (std::size_t m = 0; m < width; ++m) {
    Scalar value = a[m];
    const auto last = std::min(m, b.size() - 1);
    for (std::size_t j = 1; j <= last; ++j)
      value -= b[j] * quotient[m - j];
    quotient[m] = value / b.front();
  }
  return EpsilonFrame<Scalar>(minimum, std::move(quotient));
}

// A transformation entry is an exact finite Laurent polynomial: powers not
// present in `terms` are structural zero at every order, including beyond
// the largest stored power.  This support contract is intentionally distinct
// from EpsilonFrame, whose coefficients above CompleteMax are unknown.
// `Scalar` remains generic: for Acb the support is exact even though a stored
// coefficient is an enclosure and therefore must never be omitted merely
// because it happens to contain zero.
template <typename Scalar>
class ExactLaurentPolynomial {
 public:
  using Terms = std::map<std::int32_t, Scalar>;

  ExactLaurentPolynomial() = default;

  explicit ExactLaurentPolynomial(
      std::vector<std::pair<std::int32_t, Scalar>> terms) {
    for (auto& [power, coefficient] : terms)
      add_term(power, std::move(coefficient));
  }

  static ExactLaurentPolynomial zero() { return {}; }

  static ExactLaurentPolynomial one() {
    return monomial(0, ScalarTraits<Scalar>::one());
  }

  static ExactLaurentPolynomial monomial(std::int32_t power,
                                         Scalar coefficient) {
    ExactLaurentPolynomial result;
    result.add_term(power, std::move(coefficient));
    return result;
  }

  [[nodiscard]] bool is_zero() const { return terms_.empty(); }
  [[nodiscard]] const Terms& terms() const { return terms_; }

  [[nodiscard]] std::optional<std::int32_t> minimum_power() const {
    if (terms_.empty()) return std::nullopt;
    return terms_.begin()->first;
  }

  [[nodiscard]] std::optional<std::int32_t> maximum_power() const {
    if (terms_.empty()) return std::nullopt;
    return terms_.rbegin()->first;
  }

  [[nodiscard]] Scalar coefficient(std::int32_t power) const {
    const auto found = terms_.find(power);
    return found == terms_.end() ? ScalarTraits<Scalar>::zero()
                                : found->second;
  }

  void add_term(std::int32_t power, Scalar coefficient) {
    // `is_zero` is exact for all native scalar domains.  In particular it is
    // acb_is_zero, not acb_contains_zero, for ComplexBall.
    if (ScalarTraits<Scalar>::is_zero(coefficient)) return;
    const auto found = terms_.find(power);
    if (found == terms_.end()) {
      terms_.emplace(power, std::move(coefficient));
      return;
    }
    found->second += coefficient;
    if (ScalarTraits<Scalar>::is_zero(found->second)) terms_.erase(found);
  }

  [[nodiscard]] ExactLaurentPolynomial shifted(std::int32_t amount) const {
    ExactLaurentPolynomial result;
    for (const auto& [power, coefficient] : terms_)
      result.add_term(matching_detail::checked_power(
                          static_cast<std::int64_t>(power) + amount,
                          "exact Laurent-polynomial shift"),
                      coefficient);
    return result;
  }

  [[nodiscard]] ExactLaurentPolynomial scaled(const Scalar& factor) const {
    ExactLaurentPolynomial result;
    for (const auto& [power, coefficient] : terms_)
      result.add_term(power, coefficient * factor);
    return result;
  }

  friend ExactLaurentPolynomial operator+(const ExactLaurentPolynomial& left,
                                           const ExactLaurentPolynomial& right) {
    auto result = left;
    for (const auto& [power, coefficient] : right.terms_)
      result.add_term(power, coefficient);
    return result;
  }

  friend ExactLaurentPolynomial operator-(const ExactLaurentPolynomial& value) {
    ExactLaurentPolynomial result;
    for (const auto& [power, coefficient] : value.terms_)
      result.add_term(power, -coefficient);
    return result;
  }

  friend ExactLaurentPolynomial operator-(const ExactLaurentPolynomial& left,
                                           const ExactLaurentPolynomial& right) {
    return left + (-right);
  }

  friend ExactLaurentPolynomial operator*(const ExactLaurentPolynomial& left,
                                           const ExactLaurentPolynomial& right) {
    ExactLaurentPolynomial result;
    for (const auto& [left_power, left_coefficient] : left.terms_)
      for (const auto& [right_power, right_coefficient] : right.terms_)
        result.add_term(
            matching_detail::checked_power(
                static_cast<std::int64_t>(left_power) + right_power,
                "exact Laurent-polynomial product"),
            left_coefficient * right_coefficient);
    return result;
  }

 private:
  Terms terms_;
};

template <typename Scalar>
using FiniteLaurentVector = std::vector<EpsilonFrame<Scalar>>;

template <typename Scalar>
using FiniteLaurentMatrix = std::vector<FiniteLaurentVector<Scalar>>;

template <typename Scalar>
using ExactLaurentMatrix =
    std::vector<std::vector<ExactLaurentPolynomial<Scalar>>>;

// Exact-polynomial multiplication preserves the whole finite input width.
// For a polynomial with lowest power p, the result is complete through
// frame.CompleteMax+p: a higher output coefficient could depend on the first
// unknown input coefficient through the p monomial.  This is why a constant
// or eps^-1 transformation must not be represented as a one-row frame.
template <typename Scalar>
EpsilonFrame<Scalar> exact_laurent_times_frame(
    const ExactLaurentPolynomial<Scalar>& polynomial,
    const EpsilonFrame<Scalar>& frame) {
  if (polynomial.is_zero())
    return EpsilonFrame<Scalar>(
        frame.window(),
        std::vector<Scalar>(frame.coefficients().size(),
                            ScalarTraits<Scalar>::zero()));

  const auto polynomial_minimum = *polynomial.minimum_power();
  const auto minimum = matching_detail::checked_power(
      static_cast<std::int64_t>(frame.min_power()) + polynomial_minimum,
      "exact-polynomial frame-product minimum");
  const auto maximum = matching_detail::checked_power(
      static_cast<std::int64_t>(frame.complete_max()) + polynomial_minimum,
      "exact-polynomial frame-product complete maximum");
  std::vector<Scalar> result(frame.coefficients().size(),
                             ScalarTraits<Scalar>::zero());
  for (std::int64_t k = minimum; k <= maximum; ++k) {
    Scalar value = ScalarTraits<Scalar>::zero();
    for (const auto& [power, coefficient] : polynomial.terms()) {
      const auto frame_power = k - power;
      if (frame_power < frame.min_power()) continue;
      // k <= CompleteMax + polynomial_minimum and power >= minimum imply
      // frame_power <= CompleteMax, so no unknown coefficient is read.
      value += coefficient *
               frame.coefficient(static_cast<std::int32_t>(frame_power));
    }
    result[static_cast<std::size_t>(k - minimum)] = std::move(value);
  }
  return EpsilonFrame<Scalar>({minimum, maximum}, std::move(result));
}

template <typename Scalar>
ExactLaurentMatrix<Scalar> identity_exact_laurent_matrix(std::size_t size) {
  if (size == 0)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        "an exact Laurent transformation cannot have zero dimension");
  ExactLaurentMatrix<Scalar> result(
      size, std::vector<ExactLaurentPolynomial<Scalar>>(size));
  for (std::size_t i = 0; i < size; ++i)
    result[i][i] = ExactLaurentPolynomial<Scalar>::one();
  return result;
}

template <typename Scalar>
ExactLaurentMatrix<Scalar> multiply_exact_laurent_matrices(
    const ExactLaurentMatrix<Scalar>& left,
    const ExactLaurentMatrix<Scalar>& right) {
  const auto inner = matching_detail::rectangular_columns(
      left, "left exact Laurent transformation");
  const auto output_columns = matching_detail::rectangular_columns(
      right, "right exact Laurent transformation");
  if (right.size() != inner)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        "exact Laurent transformation dimensions do not compose");
  ExactLaurentMatrix<Scalar> result(
      left.size(),
      std::vector<ExactLaurentPolynomial<Scalar>>(output_columns));
  for (std::size_t row = 0; row < left.size(); ++row)
    for (std::size_t column = 0; column < output_columns; ++column)
      for (std::size_t k = 0; k < inner; ++k)
        result[row][column] =
            result[row][column] + left[row][k] * right[k][column];
  return result;
}

template <typename Scalar>
FiniteLaurentVector<Scalar> apply_exact_laurent_matrix(
    const ExactLaurentMatrix<Scalar>& transformation,
    const FiniteLaurentVector<Scalar>& input) {
  const auto columns = matching_detail::rectangular_columns(
      transformation, "exact Laurent transformation");
  if (columns != input.size())
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        "exact Laurent transformation and vector dimensions disagree");

  FiniteLaurentVector<Scalar> result;
  result.reserve(transformation.size());
  for (std::size_t row = 0; row < transformation.size(); ++row) {
    std::optional<EpsilonFrame<Scalar>> value;
    for (std::size_t column = 0; column < columns; ++column) {
      if (transformation[row][column].is_zero()) continue;
      auto term = exact_laurent_times_frame(
          transformation[row][column], input[column]);
      value = value.has_value() ? *value + term : std::move(term);
    }
    if (!value.has_value())
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::StructurallySingularTransformation,
          "exact Laurent transformation has a structurally zero row", row);
    result.push_back(std::move(*value));
  }
  return result;
}

// This is the arithmetic needed for the two future matching identities:
// form F*T once from the raw evaluated basis, and separately form T*w before
// checking F*(T*w).  No matching policy is encoded here.
template <typename Scalar>
FiniteLaurentMatrix<Scalar> right_multiply_finite_by_exact_laurent(
    const FiniteLaurentMatrix<Scalar>& finite,
    const ExactLaurentMatrix<Scalar>& transformation) {
  const auto inner = matching_detail::rectangular_columns(
      finite, "finite Laurent matrix");
  const auto output_columns = matching_detail::rectangular_columns(
      transformation, "exact Laurent transformation");
  if (transformation.size() != inner)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        "finite and exact Laurent matrix dimensions do not compose");

  FiniteLaurentMatrix<Scalar> result(finite.size());
  for (std::size_t row = 0; row < finite.size(); ++row) {
    result[row].reserve(output_columns);
    for (std::size_t column = 0; column < output_columns; ++column) {
      std::optional<EpsilonFrame<Scalar>> value;
      for (std::size_t k = 0; k < inner; ++k) {
        if (transformation[k][column].is_zero()) continue;
        auto term = exact_laurent_times_frame(transformation[k][column],
                                              finite[row][k]);
        value = value.has_value() ? *value + term : std::move(term);
      }
      if (!value.has_value())
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::StructurallySingularTransformation,
            "exact Laurent transformation has a structurally zero column",
            std::nullopt, column);
      result[row].push_back(std::move(*value));
    }
  }
  return result;
}

template <typename Scalar>
FiniteLaurentVector<Scalar> apply_finite_laurent_matrix(
    const FiniteLaurentMatrix<Scalar>& matrix,
    const FiniteLaurentVector<Scalar>& input) {
  const auto columns = matching_detail::rectangular_columns(
      matrix, "finite Laurent matrix");
  if (columns != input.size())
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        "finite Laurent matrix and vector dimensions disagree");
  FiniteLaurentVector<Scalar> result;
  result.reserve(matrix.size());
  for (const auto& row : matrix) {
    auto value = row.front() * input.front();
    for (std::size_t column = 1; column < columns; ++column)
      value = value + row[column] * input[column];
    result.push_back(std::move(value));
  }
  return result;
}

// Deterministic full-pivot Gaussian elimination over finite Laurent frames.
// It is intentionally only an arithmetic primitive: epsilon-lattice
// saturation, tolerance policy, refinement and residual certification belong
// to the future matching orchestration.  Pivots are ordered by exact Laurent
// valuation; ties use row/column order.  Any Acb zero ambiguity encountered
// while deciding a valuation is loud.
template <typename Scalar>
FiniteLaurentVector<Scalar> solve_finite_laurent_system(
    FiniteLaurentMatrix<Scalar> matrix,
    FiniteLaurentVector<Scalar> right_hand_side,
    const std::string& context = "finite Laurent linear solve") {
  const auto size = matching_detail::rectangular_columns(
      matrix, "finite Laurent coefficient matrix");
  if (matrix.size() != size || right_hand_side.size() != size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": coefficient matrix must be square and match the rhs");

  for (std::size_t row = 0; row < size; ++row)
    for (std::size_t column = 0; column < size; ++column)
      matrix[row][column] = matching_detail::canonical_leading_frame(
          matrix[row][column], context + ": input valuation", row, column);

  std::vector<std::size_t> column_permutation(size);
  for (std::size_t i = 0; i < size; ++i) column_permutation[i] = i;

  for (std::size_t position = 0; position < size; ++position) {
    struct Pivot {
      std::size_t row;
      std::size_t column;
      std::int32_t power;
    };
    std::optional<Pivot> pivot;
    for (std::size_t row = position; row < size; ++row) {
      for (std::size_t column = position; column < size; ++column) {
        const auto power = finite_laurent_leading_power(
            matrix[row][column], context + ": pivot search");
        if (power.has_value() &&
            (!pivot.has_value() || *power < pivot->power))
          pivot = Pivot{row, column, *power};
      }
    }
    if (!pivot.has_value())
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::SingularOrIncompleteSystem,
          context +
              ": no proved nonzero pivot remains in the complete epsilon "
              "windows",
          position, position);

    if (pivot->row != position) {
      std::swap(matrix[position], matrix[pivot->row]);
      std::swap(right_hand_side[position], right_hand_side[pivot->row]);
    }
    if (pivot->column != position) {
      for (auto& row : matrix)
        std::swap(row[position], row[pivot->column]);
      std::swap(column_permutation[position],
                column_permutation[pivot->column]);
    }

    for (std::size_t row = position + 1; row < size; ++row) {
      const auto factor = finite_laurent_quotient(
          matrix[row][position], matrix[position][position],
          context + ": elimination quotient");
      for (std::size_t column = position; column < size; ++column) {
        matrix[row][column] = matching_detail::canonical_leading_frame(
            matrix[row][column] - factor * matrix[position][column],
            context + ": Schur cancellation", row, column);
      }
      right_hand_side[row] =
          right_hand_side[row] - factor * right_hand_side[position];
    }
  }

  FiniteLaurentVector<Scalar> permuted_solution(size,
                                                right_hand_side.front());
  for (std::size_t reverse = size; reverse-- > 0;) {
    auto residual = right_hand_side[reverse];
    for (std::size_t column = reverse + 1; column < size; ++column)
      residual = residual - matrix[reverse][column] *
                                permuted_solution[column];
    permuted_solution[reverse] = finite_laurent_quotient(
        residual, matrix[reverse][reverse],
        context + ": back-substitution quotient");
  }

  FiniteLaurentVector<Scalar> solution(size, permuted_solution.front());
  for (std::size_t current = 0; current < size; ++current)
    solution[column_permutation[current]] =
        std::move(permuted_solution[current]);
  return solution;
}

}  // namespace diffexp2
