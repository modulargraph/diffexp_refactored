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
  ExponentOverflow,
  InsufficientCompleteWindow,
  InvalidSaturationLattice,
  SaturationFailure
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
struct EpsilonLatticeSaturationAction {
  std::size_t leading_rank_before = 0;
  std::size_t target_column = 0;
  std::vector<Scalar> null_relation;
};

template <typename Scalar>
struct EpsilonLatticeSaturationDiagnostics {
  std::vector<std::int32_t> initial_column_valuations;
  std::vector<std::int32_t> initial_column_shifts;
  EpsilonFrame<Scalar> normalized_determinant;
  std::int32_t normalized_determinant_valuation = 0;
  std::size_t initial_leading_rank = 0;
  std::size_t final_leading_rank = 0;
  std::vector<EpsilonLatticeSaturationAction<Scalar>> actions;
};

template <typename Scalar>
struct EpsilonLatticeSaturationResult {
  // This is the sequentially evaluated product F*T.  Its windows include
  // every honest coefficient loss caused by an epsilon division.
  FiniteLaurentMatrix<Scalar> basis_times_transformation;
  // Unlike EpsilonFrame, every omitted coefficient in T is structural zero.
  ExactLaurentMatrix<Scalar> transformation;
  EpsilonLatticeSaturationDiagnostics<Scalar> diagnostics;
};

namespace matching_detail {

template <typename Scalar>
using DenseScalarMatrix = std::vector<std::vector<Scalar>>;

template <typename Scalar>
using TruncatedPolynomial = std::vector<Scalar>;

template <typename Scalar>
using TruncatedPolynomialMatrix =
    std::vector<std::vector<TruncatedPolynomial<Scalar>>>;

template <typename Scalar>
TruncatedPolynomial<Scalar> zero_polynomial(std::size_t width) {
  return TruncatedPolynomial<Scalar>(width, ScalarTraits<Scalar>::zero());
}

template <typename Scalar>
TruncatedPolynomial<Scalar> add_polynomials(
    TruncatedPolynomial<Scalar> left,
    const TruncatedPolynomial<Scalar>& right) {
  for (std::size_t i = 0; i < left.size(); ++i) left[i] += right[i];
  return left;
}

template <typename Scalar>
TruncatedPolynomial<Scalar> multiply_polynomials(
    const TruncatedPolynomial<Scalar>& left,
    const TruncatedPolynomial<Scalar>& right) {
  auto result = zero_polynomial<Scalar>(left.size());
  for (std::size_t i = 0; i < left.size(); ++i) {
    if (ScalarTraits<Scalar>::is_zero(left[i])) continue;
    for (std::size_t j = 0; j + i < right.size(); ++j) {
      if (ScalarTraits<Scalar>::is_zero(right[j])) continue;
      result[i + j] += left[i] * right[j];
    }
  }
  return result;
}

template <typename Scalar>
TruncatedPolynomialMatrix<Scalar> multiply_polynomial_matrices(
    const TruncatedPolynomialMatrix<Scalar>& left,
    const TruncatedPolynomialMatrix<Scalar>& right) {
  const auto size = left.size();
  const auto width = left.front().front().size();
  TruncatedPolynomialMatrix<Scalar> result(
      size, std::vector<TruncatedPolynomial<Scalar>>(
                size, zero_polynomial<Scalar>(width)));
  for (std::size_t row = 0; row < size; ++row)
    for (std::size_t column = 0; column < size; ++column)
      for (std::size_t k = 0; k < size; ++k)
        result[row][column] = add_polynomials(
            std::move(result[row][column]),
            multiply_polynomials(left[row][k], right[k][column]));
  return result;
}

// Newton's identities evaluate the determinant in polynomial time while
// staying inside one finite, honest coefficient rectangle.  All normalized
// entries are known through `common_top`, so every determinant coefficient
// through that same order is complete.
template <typename Scalar>
EpsilonFrame<Scalar> nonnegative_determinant_frame(
    const FiniteLaurentMatrix<Scalar>& matrix,
    const std::string& context) {
  const auto size = rectangular_columns(matrix, "saturation basis matrix");
  if (matrix.size() != size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": determinant requires a square basis matrix");

  std::int32_t common_top = std::numeric_limits<std::int32_t>::max();
  for (std::size_t row = 0; row < size; ++row) {
    for (std::size_t column = 0; column < size; ++column) {
      const auto& entry = matrix[row][column];
      if (entry.complete_max() < 0)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            context +
                ": normalized basis is not complete through epsilon^0",
            row, column, entry.complete_max());
      const auto leading = finite_laurent_leading_power(
          entry, context + ": determinant input valuation");
      if (leading.has_value() && *leading < 0)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InvalidSaturationLattice,
            context +
                ": column normalization retained a negative valuation",
            row, column, *leading);
      common_top = std::min(common_top, entry.complete_max());
    }
  }

  const auto width = static_cast<std::size_t>(
      static_cast<std::int64_t>(common_top) + 1);
  TruncatedPolynomialMatrix<Scalar> value(
      size, std::vector<TruncatedPolynomial<Scalar>>(
                size, zero_polynomial<Scalar>(width)));
  for (std::size_t row = 0; row < size; ++row)
    for (std::size_t column = 0; column < size; ++column)
      for (std::int64_t power = 0; power <= common_top; ++power)
        value[row][column][static_cast<std::size_t>(power)] =
            matrix[row][column].coefficient(
                static_cast<std::int32_t>(power));

  auto one = zero_polynomial<Scalar>(width);
  one.front() = ScalarTraits<Scalar>::one();
  std::vector<TruncatedPolynomial<Scalar>> traces;
  traces.reserve(size);
  std::vector<TruncatedPolynomial<Scalar>> elementary;
  elementary.reserve(size + 1);
  elementary.push_back(one);
  auto power_matrix = value;

  for (std::size_t degree = 1; degree <= size; ++degree) {
    auto trace = zero_polynomial<Scalar>(width);
    for (std::size_t i = 0; i < size; ++i)
      trace = add_polynomials(std::move(trace), power_matrix[i][i]);
    traces.push_back(std::move(trace));

    auto sum = zero_polynomial<Scalar>(width);
    for (std::size_t i = 1; i <= degree; ++i) {
      auto term = multiply_polynomials(elementary[degree - i],
                                       traces[i - 1]);
      if (i % 2 == 0)
        for (auto& coefficient : term) coefficient = -coefficient;
      sum = add_polynomials(std::move(sum), term);
    }
    if (degree > static_cast<std::size_t>(
                     std::numeric_limits<long>::max()))
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::DimensionMismatch,
          context + ": determinant dimension exceeds scalar integer range");
    const auto divisor =
        ScalarTraits<Scalar>::integer(static_cast<long>(degree));
    for (auto& coefficient : sum) coefficient = coefficient / divisor;
    elementary.push_back(std::move(sum));
    if (degree < size)
      power_matrix = multiply_polynomial_matrices(power_matrix, value);
  }

  return canonical_leading_frame(
      EpsilonFrame<Scalar>(0, std::move(elementary.back())),
      context + ": determinant valuation");
}

template <typename Scalar>
struct LeadingNullRelation {
  std::size_t rank = 0;
  std::optional<std::vector<Scalar>> vector;
  std::optional<std::size_t> target_column;
};

// Full-pivot elimination at epsilon^0.  The first proved nonzero entry in
// row-major order is the deterministic pivot; no magnitude sampling or
// tolerance enters the formal rank.  If the rank is deficient, the first
// free permuted column gives a right-null relation whose target coefficient
// is exactly one, making the later elementary column action invertible.
template <typename Scalar>
LeadingNullRelation<Scalar> leading_null_relation(
    DenseScalarMatrix<Scalar> matrix, const std::string& context) {
  const auto size = rectangular_columns(matrix, "leading coefficient matrix");
  if (matrix.size() != size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": leading coefficient matrix must be square");
  const auto original = matrix;
  std::vector<std::size_t> column_permutation(size);
  for (std::size_t i = 0; i < size; ++i) column_permutation[i] = i;

  std::size_t rank = 0;
  for (std::size_t position = 0; position < size; ++position) {
    std::optional<std::pair<std::size_t, std::size_t>> pivot;
    for (std::size_t row = position; row < size; ++row) {
      for (std::size_t column = position; column < size; ++column) {
        const auto decision = zero_decision(matrix[row][column]);
        if (decision == ZeroDecision::Ambiguous)
          throw MatchingArithmeticError(
              MatchingArithmeticErrorCode::AmbiguousZero,
              context +
                  ": Acb enclosure overlaps zero during saturation rank",
              row, column, 0);
        if (decision == ZeroDecision::Nonzero && !pivot.has_value())
          pivot = std::pair{row, column};
      }
    }
    if (!pivot.has_value()) break;

    if (pivot->first != position)
      std::swap(matrix[position], matrix[pivot->first]);
    if (pivot->second != position) {
      for (auto& row : matrix)
        std::swap(row[position], row[pivot->second]);
      std::swap(column_permutation[position],
                column_permutation[pivot->second]);
    }

    for (std::size_t row = position + 1; row < size; ++row) {
      const auto decision = zero_decision(matrix[row][position]);
      if (decision == ZeroDecision::Ambiguous)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::AmbiguousZero,
            context +
                ": Acb enclosure overlaps zero below a saturation pivot",
            row, position, 0);
      if (decision == ZeroDecision::Zero) continue;
      const auto factor =
          matrix[row][position] / matrix[position][position];
      matrix[row][position] = ScalarTraits<Scalar>::zero();
      for (std::size_t column = position + 1; column < size; ++column)
        matrix[row][column] -= factor * matrix[position][column];
    }
    ++rank;
  }

  if (rank == size) return LeadingNullRelation<Scalar>{rank, {}, {}};

  std::vector<Scalar> permuted(size, ScalarTraits<Scalar>::zero());
  permuted[rank] = ScalarTraits<Scalar>::one();
  for (std::size_t reverse = rank; reverse-- > 0;) {
    auto value = ScalarTraits<Scalar>::zero();
    for (std::size_t column = reverse + 1; column < size; ++column)
      value += matrix[reverse][column] * permuted[column];
    permuted[reverse] = -value / matrix[reverse][reverse];
  }

  std::vector<Scalar> relation(size, ScalarTraits<Scalar>::zero());
  for (std::size_t column = 0; column < size; ++column)
    relation[column_permutation[column]] = std::move(permuted[column]);
  const auto target = column_permutation[rank];

  for (std::size_t row = 0; row < size; ++row) {
    auto residual = ScalarTraits<Scalar>::zero();
    for (std::size_t column = 0; column < size; ++column)
      residual += original[row][column] * relation[column];
    const auto decision = zero_decision(residual);
    if (decision == ZeroDecision::Ambiguous)
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::AmbiguousZero,
          context +
              ": Acb null-relation residual overlaps zero and is not exact",
          row, std::nullopt, 0);
    if (decision == ZeroDecision::Nonzero)
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::SaturationFailure,
          context + ": computed saturation relation is not in the nullspace",
          row, std::nullopt, 0);
  }
  return LeadingNullRelation<Scalar>{rank, std::move(relation), target};
}

template <typename Scalar>
DenseScalarMatrix<Scalar> epsilon_zero_matrix(
    const FiniteLaurentMatrix<Scalar>& matrix,
    const std::string& context) {
  const auto size = matrix.size();
  DenseScalarMatrix<Scalar> result(
      size, std::vector<Scalar>(size, ScalarTraits<Scalar>::zero()));
  for (std::size_t row = 0; row < size; ++row) {
    for (std::size_t column = 0; column < size; ++column) {
      const auto& entry = matrix[row][column];
      if (entry.complete_max() < 0)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            context + ": basis entry is not complete through epsilon^0",
            row, column, entry.complete_max());
      result[row][column] = entry.coefficient(0);
    }
  }
  return result;
}

template <typename Scalar>
EpsilonFrame<Scalar> saturation_divide(
    const FiniteLaurentVector<Scalar>& row,
    const std::vector<Scalar>& relation,
    const std::string& context) {
  std::optional<EpsilonFrame<Scalar>> combination;
  for (std::size_t column = 0; column < relation.size(); ++column) {
    if (ScalarTraits<Scalar>::is_zero(relation[column])) continue;
    auto term = row[column].scaled(relation[column]);
    combination = combination.has_value()
                      ? *combination + term
                      : std::move(term);
  }
  if (!combination.has_value())
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::SaturationFailure,
        context + ": saturation null relation has empty support");
  if (combination->min_power() < 0)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InvalidSaturationLattice,
        context + ": saturation combination retained a negative valuation",
        std::nullopt, std::nullopt, combination->min_power());

  if (combination->min_power() == 0) {
    const auto decision = zero_decision(combination->coefficient(0));
    if (decision == ZeroDecision::Ambiguous)
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::AmbiguousZero,
          context +
              ": Acb enclosure overlaps zero at epsilon divisibility check",
          std::nullopt, std::nullopt, 0);
    if (decision == ZeroDecision::Nonzero)
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::SaturationFailure,
          context +
              ": null-column combination is not divisible by epsilon",
          std::nullopt, std::nullopt, 0);
    if (combination->complete_max() == 0)
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::InsufficientCompleteWindow,
          context +
              ": no complete coefficient remains after epsilon division");
    combination = EpsilonFrame<Scalar>(
        {1, combination->complete_max()},
        std::vector<Scalar>(combination->coefficients().begin() + 1,
                            combination->coefficients().end()));
  }

  auto canonical = canonical_leading_frame(
      *combination, context + ": saturation quotient");
  const auto shifted_min = checked_power(
      static_cast<std::int64_t>(canonical.min_power()) - 1,
      "saturation quotient minimum");
  const auto shifted_max = checked_power(
      static_cast<std::int64_t>(canonical.complete_max()) - 1,
      "saturation quotient complete maximum");
  if (shifted_max < 0)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InsufficientCompleteWindow,
        context +
            ": epsilon division lost the complete leading coefficient");
  return EpsilonFrame<Scalar>({shifted_min, shifted_max},
                              canonical.coefficients());
}

}  // namespace matching_detail

// Deterministically saturate the epsilon lattice generated by the columns
// of a finite Laurent basis F.  Initial column monomials normalize their
// valuations.  Every deficient epsilon^0 frame then supplies one certified
// null relation a and replaces a target column by (Sum_i a_i F_i)/epsilon.
// The same right actions are accumulated in the exact-support Laurent matrix
// T, maintaining the invariant `basis_times_transformation == F*T`.
template <typename Scalar>
EpsilonLatticeSaturationResult<Scalar> saturate_finite_laurent_basis(
    const FiniteLaurentMatrix<Scalar>& basis,
    const std::string& context = "finite Laurent basis saturation") {
  const auto size = matching_detail::rectangular_columns(
      basis, "finite Laurent saturation basis");
  if (basis.size() != size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": saturation basis must be square");

  auto transformed = basis;
  for (std::size_t row = 0; row < size; ++row)
    for (std::size_t column = 0; column < size; ++column)
      transformed[row][column] = matching_detail::canonical_leading_frame(
          transformed[row][column], context + ": input valuation", row,
          column);

  std::vector<std::int32_t> valuations(size);
  std::vector<std::int32_t> shifts(size);
  auto transformation = identity_exact_laurent_matrix<Scalar>(size);
  for (std::size_t column = 0; column < size; ++column) {
    std::optional<std::int32_t> valuation;
    for (std::size_t row = 0; row < size; ++row) {
      const auto candidate = finite_laurent_leading_power(
          transformed[row][column], context + ": column valuation");
      if (candidate.has_value() &&
          (!valuation.has_value() || *candidate < *valuation))
        valuation = candidate;
    }
    if (!valuation.has_value())
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::SingularOrIncompleteSystem,
          context +
              ": a basis column is zero throughout its complete window",
          std::nullopt, column);
    valuations[column] = *valuation;
    shifts[column] = matching_detail::checked_power(
        -static_cast<std::int64_t>(*valuation),
        "saturation initial column shift");
    transformation[column][column] =
        ExactLaurentPolynomial<Scalar>::monomial(
            shifts[column], ScalarTraits<Scalar>::one());
    for (std::size_t row = 0; row < size; ++row) {
      const auto shifted_min = matching_detail::checked_power(
          static_cast<std::int64_t>(
              transformed[row][column].min_power()) + shifts[column],
          "normalized column minimum");
      const auto shifted_max = matching_detail::checked_power(
          static_cast<std::int64_t>(
              transformed[row][column].complete_max()) + shifts[column],
          "normalized column complete maximum");
      transformed[row][column] = EpsilonFrame<Scalar>(
          {shifted_min, shifted_max},
          transformed[row][column].coefficients());
      transformed[row][column] = matching_detail::canonical_leading_frame(
          transformed[row][column], context + ": normalized column", row,
          column);
      if (transformed[row][column].complete_max() < 0)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            context +
                ": column normalization lost the epsilon^0 coefficient",
            row, column, transformed[row][column].complete_max());
    }
  }

  auto determinant = matching_detail::nonnegative_determinant_frame(
      transformed, context);
  const auto determinant_valuation = finite_laurent_leading_power(
      determinant, context + ": normalized determinant");
  if (!determinant_valuation.has_value())
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::SingularOrIncompleteSystem,
        context +
            ": determinant is unresolved in the complete epsilon window");
  if (*determinant_valuation < 0)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InvalidSaturationLattice,
        context + ": normalized determinant retained an epsilon pole",
        std::nullopt, std::nullopt, *determinant_valuation);

  std::vector<EpsilonLatticeSaturationAction<Scalar>> actions;
  auto relation = matching_detail::leading_null_relation(
      matching_detail::epsilon_zero_matrix(transformed, context),
      context + "#sat0");
  const auto initial_rank = relation.rank;
  std::size_t steps = 0;
  while (relation.rank != size) {
    if (steps >= static_cast<std::size_t>(*determinant_valuation))
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::SaturationFailure,
          context +
              ": determinant valuation was consumed before full rank");
    if (!relation.vector.has_value() ||
        !relation.target_column.has_value())
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::SaturationFailure,
          context + ": deficient leading frame has no null relation");

    const auto& null_vector = *relation.vector;
    const auto target = *relation.target_column;
    for (std::size_t row = 0; row < size; ++row)
      transformed[row][target] = matching_detail::saturation_divide(
          transformed[row], null_vector,
          context + "#sat" + std::to_string(steps + 1) + "/row" +
              std::to_string(row));

    auto elementary = identity_exact_laurent_matrix<Scalar>(size);
    for (std::size_t row = 0; row < size; ++row)
      elementary[row][target] =
          ExactLaurentPolynomial<Scalar>::monomial(-1, null_vector[row]);
    transformation = multiply_exact_laurent_matrices(transformation,
                                                       elementary);
    actions.push_back(
        {relation.rank, target, std::move(*relation.vector)});
    ++steps;
    relation = matching_detail::leading_null_relation(
        matching_detail::epsilon_zero_matrix(transformed, context),
        context + "#sat" + std::to_string(steps));
  }

  if (steps != static_cast<std::size_t>(*determinant_valuation))
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::SaturationFailure,
        context +
            ": full rank did not consume the determinant valuation exactly");

  EpsilonLatticeSaturationDiagnostics<Scalar> diagnostics{
      std::move(valuations),
      std::move(shifts),
      std::move(determinant),
      *determinant_valuation,
      initial_rank,
      relation.rank,
      std::move(actions)};
  return {std::move(transformed), std::move(transformation),
          std::move(diagnostics)};
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
      const auto below_pivot = finite_laurent_leading_power(
          matrix[row][position],
          context + ": below-pivot elimination entry");
      if (!below_pivot.has_value()) continue;
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
