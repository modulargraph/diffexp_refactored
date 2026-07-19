#pragma once

#include "diffexp2/integration.hpp"

#include <flint/acb_mat.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <tuple>
#include <utility>
#include <vector>

namespace diffexp2 {

enum class MatchingArithmeticErrorCode : std::uint8_t {
  DimensionMismatch,
  AmbiguousZero,
  ZeroDivisor,
  SingularOrIncompleteSystem,
  UnresolvedDeterminantTail,
  StructurallySingularTransformation,
  ExponentOverflow,
  InsufficientCompleteWindow,
  InvalidSaturationLattice,
  SaturationFailure,
  SearchBudgetExhausted
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

struct CertifiedLaurentLeadingPower {
  std::optional<std::int32_t> power;
  std::optional<std::int32_t> first_ambiguous_power;
};

// Probe whether a frame is usable as a Laurent pivot without forcing a
// decision on a zero-overlapping Acb coefficient.  Once an ambiguous
// coefficient is encountered, a later proved-nonzero coefficient cannot be
// called the leading one, so this candidate is skipped.  Exact-zero frames
// remain distinguishable from ambiguous candidates.
template <typename Scalar>
CertifiedLaurentLeadingPower certified_laurent_leading_power(
    const EpsilonFrame<Scalar>& input) {
  for (std::size_t index = 0; index < input.coefficients().size(); ++index) {
    const auto power = checked_power(
        static_cast<std::int64_t>(input.min_power()) +
            static_cast<std::int64_t>(index),
        "certified Laurent leading power");
    switch (zero_decision(input.coefficients()[index])) {
      case ZeroDecision::Zero:
        continue;
      case ZeroDecision::Ambiguous:
        return {std::nullopt, power};
      case ZeroDecision::Nonzero:
        return {power, std::nullopt};
    }
  }
  return {};
}

// Canonicalize exactly as far as the enclosure proves.  A certified nonzero
// leading coefficient permits removal of the exact-zero prefix; a wholly
// certified-zero frame becomes structural zero through CompleteMax.  If the
// first undecided coefficient overlaps zero, preserve the original frame and
// its full enclosure instead of guessing its valuation.
template <typename Scalar>
EpsilonFrame<Scalar> canonicalize_certified_or_preserve_ambiguous(
    const EpsilonFrame<Scalar>& input, const std::string& context,
    std::optional<std::size_t> row = std::nullopt,
    std::optional<std::size_t> column = std::nullopt) {
  const auto leading = certified_laurent_leading_power(input);
  if (leading.first_ambiguous_power.has_value()) return input;
  if (!leading.power.has_value())
    return EpsilonFrame<Scalar>::zero(input.complete_max());
  return canonical_leading_frame(input, context, row, column);
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

// Specialize an exact-rational Laurent transformation coefficientwise while
// retaining its exact support.  In particular, absence from `terms()` remains
// structural zero; an inexact Acb coefficient is never used to decide whether
// a monomial exists.
inline ExactLaurentMatrix<ComplexBall>
specialize_exact_rational_laurent_matrix_to_acb(
    const ExactLaurentMatrix<Rational>& exact) {
  const auto columns = matching_detail::rectangular_columns(
      exact, "exact-rational Laurent transformation");
  ExactLaurentMatrix<ComplexBall> result(
      exact.size(),
      std::vector<ExactLaurentPolynomial<ComplexBall>>(columns));
  for (std::size_t row = 0; row < exact.size(); ++row)
    for (std::size_t column = 0; column < columns; ++column)
      for (const auto& [power, coefficient] :
           exact[row][column].terms())
        result[row][column].add_term(
            power, ComplexBall::from_strings(coefficient.str()));
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
    const std::string& context,
    std::optional<int> certified_candidate_chop_digits = std::nullopt) {
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

  auto coefficients = std::move(elementary.back());
  if constexpr (std::is_same_v<Scalar, ComplexBall>) {
    if (certified_candidate_chop_digits.has_value())
      for (auto& coefficient : coefficients)
        coefficient = ScalarTraits<ComplexBall>::canonicalized(
            coefficient, *certified_candidate_chop_digits);
  }
  return canonical_leading_frame(
      EpsilonFrame<Scalar>(0, std::move(coefficients)),
      context + ": determinant valuation");
}

template <typename Scalar>
struct LeadingNullRelation {
  std::size_t rank = 0;
  std::optional<std::vector<Scalar>> vector;
  std::optional<std::size_t> target_column;
  std::vector<std::size_t> pivot_rows;
  // When rank is full, this is the determinant of the selected leading
  // minor with rows ordered as pivot_rows and columns in their original
  // order.  The elimination already needed for the rank/null-relation proof
  // supplies it essentially for free.
  std::optional<Scalar> full_rank_leading_determinant;
};

// Full-pivot elimination at epsilon^0.  The first proved nonzero entry in
// row-major order is the deterministic pivot; no magnitude sampling or
// tolerance enters the formal rank.  If the rank is deficient, the first
// free permuted column gives a right-null relation whose target coefficient
// is exactly one, making the later elementary column action invertible.
template <typename Scalar>
LeadingNullRelation<Scalar> leading_null_relation(
    DenseScalarMatrix<Scalar> matrix, const std::string& context,
    std::optional<int> certified_candidate_chop_digits = std::nullopt) {
  const auto size = rectangular_columns(matrix, "leading coefficient matrix");
  const auto row_count = matrix.size();
  if (row_count < size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
      context +
            ": leading coefficient matrix must have at least as many rows as columns");
  if constexpr (std::is_same_v<Scalar, ComplexBall>) {
    if (certified_candidate_chop_digits.has_value())
      for (auto& row : matrix)
        for (auto& coefficient : row)
          coefficient = ScalarTraits<ComplexBall>::canonicalized(
              coefficient, *certified_candidate_chop_digits);
  }
  const auto original = matrix;
  std::vector<std::size_t> column_permutation(size);
  for (std::size_t i = 0; i < size; ++i) column_permutation[i] = i;
  std::vector<std::size_t> row_permutation(row_count);
  for (std::size_t i = 0; i < row_count; ++i) row_permutation[i] = i;

  std::size_t rank = 0;
  auto pivot_product = ScalarTraits<Scalar>::one();
  bool odd_column_permutation = false;
  for (std::size_t position = 0; position < size; ++position) {
    std::optional<std::pair<std::size_t, std::size_t>> pivot;
    for (std::size_t row = position; row < row_count; ++row) {
      for (std::size_t column = position; column < size; ++column) {
        const auto decision = zero_decision(matrix[row][column]);
        if (decision == ZeroDecision::Ambiguous) {
          std::string ball_detail;
          if constexpr (std::is_same_v<Scalar, ComplexBall>) {
            const auto& coefficient = matrix[row][column];
            ball_detail = "; row=" + std::to_string(row) +
                "; column=" + std::to_string(column) +
                "; midpoint=(" + coefficient.real_midpoint(16) + "," +
                coefficient.imag_midpoint(16) + ")" +
                "; radius2exp=(" +
                coefficient.real_radius_exponent() + "," +
                coefficient.imag_radius_exponent() + ")";
          }
          throw MatchingArithmeticError(
              MatchingArithmeticErrorCode::AmbiguousZero,
              context +
                  ": Acb enclosure overlaps zero during saturation rank" +
                  ball_detail,
              row, column, 0);
        }
        if (decision == ZeroDecision::Nonzero && !pivot.has_value())
          pivot = std::pair{row, column};
      }
    }
    if (!pivot.has_value()) break;

    if (pivot->first != position) {
      std::swap(matrix[position], matrix[pivot->first]);
      std::swap(row_permutation[position], row_permutation[pivot->first]);
    }
    if (pivot->second != position) {
      for (auto& row : matrix)
        std::swap(row[position], row[pivot->second]);
      std::swap(column_permutation[position],
                column_permutation[pivot->second]);
      odd_column_permutation = !odd_column_permutation;
    }

    pivot_product *= matrix[position][position];
    for (std::size_t row = position + 1; row < row_count; ++row) {
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
      for (std::size_t column = position + 1; column < size; ++column) {
        matrix[row][column] -= factor * matrix[position][column];
        if constexpr (std::is_same_v<Scalar, ComplexBall>) {
          if (certified_candidate_chop_digits.has_value())
            matrix[row][column] =
                ScalarTraits<ComplexBall>::canonicalized(
                    matrix[row][column],
                    *certified_candidate_chop_digits);
        }
      }
    }
    ++rank;
  }

  std::vector<std::size_t> pivot_rows(row_permutation.begin(),
                                      row_permutation.begin() + rank);
  if (rank == size) {
    if (odd_column_permutation) pivot_product = -pivot_product;
    return LeadingNullRelation<Scalar>{rank, {}, {},
                                       std::move(pivot_rows),
                                       std::move(pivot_product)};
  }

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

  for (std::size_t row = 0; row < row_count; ++row) {
    auto residual = ScalarTraits<Scalar>::zero();
    for (std::size_t column = 0; column < size; ++column)
      residual += original[row][column] * relation[column];
    if constexpr (std::is_same_v<Scalar, ComplexBall>) {
      if (certified_candidate_chop_digits.has_value())
        residual = ScalarTraits<ComplexBall>::canonicalized(
            residual, *certified_candidate_chop_digits);
    }
    const auto decision = zero_decision(residual);
    if (decision == ZeroDecision::Ambiguous) {
      std::string ball_detail;
      if constexpr (std::is_same_v<Scalar, ComplexBall>)
        ball_detail =
            "; midpoint=(" + residual.real_midpoint(16) + "," +
            residual.imag_midpoint(16) + ")" +
            "; radius2exp=(" + residual.real_radius_exponent() + "," +
            residual.imag_radius_exponent() + ")";
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::AmbiguousZero,
          context +
              ": Acb null-relation residual overlaps zero and is not exact" +
              ball_detail,
          row, std::nullopt, 0);
    }
    if (decision == ZeroDecision::Nonzero)
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::SaturationFailure,
          context + ": computed saturation relation is not in the nullspace",
          row, std::nullopt, 0);
  }
  return LeadingNullRelation<Scalar>{rank, std::move(relation), target,
                                     std::move(pivot_rows), std::nullopt};
}

enum class FullRankProofResult : std::uint8_t {
  Proved,
  Ambiguous,
  Deficient,
  SearchBudgetExhausted
};

template <typename Scalar>
std::optional<std::vector<std::pair<std::size_t, std::size_t>>>
certified_full_rank_plan_bounded_search(
    DenseScalarMatrix<Scalar> matrix, std::size_t position,
    std::size_t& remaining_nodes, std::size_t branch_width);

struct CertifiedPivotCandidate {
  std::size_t row = 0;
  std::size_t column = 0;
};

// Prefer the ball which is separated from zero by the largest number of
// relative-accuracy bits.  A lower absolute-modulus bound breaks the (quite
// common) accuracy tie, implementing full rather than row-only pivoting.
// Exact scalar fields deliberately retain row/column order: there is no
// enclosure quality to rank and the stable order minimizes behavior churn.
template <typename Scalar>
int compare_certified_pivot_strength(const Scalar&, const Scalar&) {
  return 0;
}

template <>
inline int compare_certified_pivot_strength<ComplexBall>(
    const ComplexBall& left, const ComplexBall& right) {
  const auto left_accuracy = acb_rel_accuracy_bits(left.raw());
  const auto right_accuracy = acb_rel_accuracy_bits(right.raw());
  if (left_accuracy != right_accuracy)
    return left_accuracy > right_accuracy ? 1 : -1;

  arf_t left_lower;
  arf_t right_lower;
  arf_init(left_lower);
  arf_init(right_lower);
  acb_get_abs_lbound_arf(left_lower, left.raw(), ComplexBall::precision());
  acb_get_abs_lbound_arf(right_lower, right.raw(), ComplexBall::precision());
  const auto comparison = arf_cmp(left_lower, right_lower);
  arf_clear(left_lower);
  arf_clear(right_lower);
  return comparison > 0 ? 1 : comparison < 0 ? -1 : 0;
}

template <typename Scalar>
std::vector<CertifiedPivotCandidate> ranked_certified_pivots(
    const DenseScalarMatrix<Scalar>& matrix, std::size_t position) {
  const auto size = matrix.size();
  std::vector<CertifiedPivotCandidate> candidates;
  candidates.reserve((size - position) * (size - position));
  for (std::size_t pivot_row = position; pivot_row < size; ++pivot_row) {
    for (std::size_t pivot_column = position; pivot_column < size;
         ++pivot_column) {
      if (zero_decision(matrix[pivot_row][pivot_column]) !=
          ZeroDecision::Nonzero)
        continue;
      candidates.push_back({pivot_row, pivot_column});
    }
  }
  std::stable_sort(
      candidates.begin(), candidates.end(),
      [&matrix](const auto& left, const auto& right) {
        const auto strength = compare_certified_pivot_strength(
            matrix[left.row][left.column], matrix[right.row][right.column]);
        if (strength != 0) return strength > 0;
        if (left.row != right.row) return left.row < right.row;
        return left.column < right.column;
      });
  return candidates;
}

template <typename Scalar>
void apply_certified_pivot(DenseScalarMatrix<Scalar>& matrix,
                           std::size_t position,
                           const CertifiedPivotCandidate& pivot) {
  const auto size = matrix.size();
  if (pivot.row != position)
    std::swap(matrix[position], matrix[pivot.row]);
  if (pivot.column != position)
    for (auto& row : matrix)
      std::swap(row[position], row[pivot.column]);

  for (std::size_t row = position + 1; row < size; ++row) {
    if (zero_decision(matrix[row][position]) == ZeroDecision::Zero) {
      matrix[row][position] = ScalarTraits<Scalar>::zero();
      continue;
    }
    const auto factor = matrix[row][position] /
                        matrix[position][position];
    matrix[row][position] = ScalarTraits<Scalar>::zero();
    for (std::size_t column = position + 1; column < size; ++column)
      matrix[row][column] -= factor * matrix[position][column];
  }
}

// The normal route is a single in-place rank-revealing elimination: O(n^3)
// arithmetic and no matrix copy per candidate.  Every chosen pivot is still
// individually certified not to contain zero, so ordering changes efficiency
// only, never the proof obligation.
template <typename Scalar>
std::optional<std::vector<std::pair<std::size_t, std::size_t>>>
certified_full_rank_greedy_plan(DenseScalarMatrix<Scalar> matrix,
                                std::size_t position) {
  std::vector<std::pair<std::size_t, std::size_t>> plan;
  plan.reserve(matrix.size() - position);
  for (; position < matrix.size(); ++position) {
    auto candidates = ranked_certified_pivots(matrix, position);
    if (candidates.empty()) return std::nullopt;
    const auto pivot = candidates.front();
    plan.emplace_back(pivot.row, pivot.column);
    apply_certified_pivot(matrix, position, pivot);
  }
  return plan;
}

// Interval widening can exceptionally make the strongest greedy route
// inconclusive even though another certified pivot order succeeds.  Explore
// alternatives in the same deterministic quality order, but never permit an
// unbounded factorial search at production dimensions.  Small matrices retain
// exhaustive existence parity with the former planner.
template <typename Scalar>
std::optional<std::vector<std::pair<std::size_t, std::size_t>>>
certified_full_rank_plan_bounded_search(
    DenseScalarMatrix<Scalar> matrix, std::size_t position,
    std::size_t& remaining_nodes, std::size_t branch_width) {
  if (position == matrix.size())
    return std::vector<std::pair<std::size_t, std::size_t>>{};
  auto candidates = ranked_certified_pivots(matrix, position);
  if (candidates.size() > branch_width) candidates.resize(branch_width);
  for (const auto& pivot : candidates) {
    if (remaining_nodes == 0) return std::nullopt;
    --remaining_nodes;
    auto next = matrix;
    apply_certified_pivot(next, position, pivot);
    auto tail = certified_full_rank_plan_bounded_search(
        std::move(next), position + 1, remaining_nodes, branch_width);
    if (!tail.has_value()) continue;
    tail->insert(tail->begin(), {pivot.row, pivot.column});
    return tail;
  }
  return std::nullopt;
}

template <typename Scalar>
std::optional<std::vector<std::pair<std::size_t, std::size_t>>>
certified_full_rank_plan(DenseScalarMatrix<Scalar> matrix,
                         std::size_t position) {
  auto greedy = certified_full_rank_greedy_plan(matrix, position);
  if (greedy.has_value()) return greedy;

  constexpr std::size_t exhaustive_size = 5;
  constexpr std::size_t production_node_budget = 32;
  constexpr std::size_t production_branch_width = 3;
  const auto remaining_dimension = matrix.size() - position;
  auto remaining_nodes = remaining_dimension <= exhaustive_size
      ? std::numeric_limits<std::size_t>::max()
      : production_node_budget;
  const auto branch_width = remaining_dimension <= exhaustive_size
      ? std::numeric_limits<std::size_t>::max()
      : production_branch_width;
  return certified_full_rank_plan_bounded_search(
      std::move(matrix), position, remaining_nodes, branch_width);
}

template <typename Scalar>
bool active_submatrix_has_ambiguous_entry(
    const DenseScalarMatrix<Scalar>& matrix, std::size_t position) {
  for (std::size_t row = position; row < matrix.size(); ++row)
    for (std::size_t column = position; column < matrix.size(); ++column)
      if (zero_decision(matrix[row][column]) == ZeroDecision::Ambiguous)
        return true;
  return false;
}

template <typename Scalar>
FullRankProofResult certified_full_rank_bounded_search(
    DenseScalarMatrix<Scalar> matrix, std::size_t position,
    std::size_t& remaining_nodes, std::size_t branch_width) {
  if (position == matrix.size()) return FullRankProofResult::Proved;
  auto candidates = ranked_certified_pivots(matrix, position);
  const auto has_ambiguous_candidate =
      active_submatrix_has_ambiguous_entry(matrix, position);
  if (candidates.empty())
    return has_ambiguous_candidate ? FullRankProofResult::Ambiguous
                                   : FullRankProofResult::Deficient;
  const auto candidates_truncated = candidates.size() > branch_width;
  if (candidates_truncated) candidates.resize(branch_width);

  bool saw_ambiguous = has_ambiguous_candidate;
  // Unvisited certified candidates are logically an exhausted search budget,
  // never evidence of rank deficiency.
  bool exhausted_budget = candidates_truncated;
  for (const auto& pivot : candidates) {
    if (remaining_nodes == 0) {
      exhausted_budget = true;
      break;
    }
    --remaining_nodes;
    auto next = matrix;
    apply_certified_pivot(next, position, pivot);
    const auto continuation = certified_full_rank_bounded_search(
        std::move(next), position + 1, remaining_nodes, branch_width);
    if (continuation == FullRankProofResult::Proved) return continuation;
    saw_ambiguous |= continuation == FullRankProofResult::Ambiguous;
    exhausted_budget |=
        continuation == FullRankProofResult::SearchBudgetExhausted;
  }
  if (exhausted_budget) return FullRankProofResult::SearchBudgetExhausted;
  return saw_ambiguous ? FullRankProofResult::Ambiguous
                       : FullRankProofResult::Deficient;
}

// Full-rank certification has weaker zero-decision requirements than lattice
// saturation: arbitrary off-pivot balls are harmless when a sequence of
// certified nonzero pivots proves invertibility.  Use the same O(n^3) greedy
// proof first as the plan-producing path.  Only interval-inconclusive greedy
// routes backtrack, exhaustively for small matrices and under a global node
// budget at production dimensions.
template <typename Scalar>
FullRankProofResult certified_full_rank_search(
    DenseScalarMatrix<Scalar> matrix, std::size_t position) {
  auto greedy_matrix = matrix;
  for (auto current = position; current < greedy_matrix.size(); ++current) {
    auto candidates = ranked_certified_pivots(greedy_matrix, current);
    if (candidates.empty()) {
      if (!active_submatrix_has_ambiguous_entry(greedy_matrix, current))
        return FullRankProofResult::Deficient;
      break;
    }
    apply_certified_pivot(greedy_matrix, current, candidates.front());
    if (current + 1 == greedy_matrix.size())
      return FullRankProofResult::Proved;
  }
  if (position == matrix.size()) return FullRankProofResult::Proved;

  constexpr std::size_t exhaustive_size = 5;
  constexpr std::size_t production_node_budget = 32;
  constexpr std::size_t production_branch_width = 3;
  const auto remaining_dimension = matrix.size() - position;
  auto remaining_nodes = remaining_dimension <= exhaustive_size
      ? std::numeric_limits<std::size_t>::max()
      : production_node_budget;
  const auto branch_width = remaining_dimension <= exhaustive_size
      ? std::numeric_limits<std::size_t>::max()
      : production_branch_width;
  return certified_full_rank_bounded_search(
      std::move(matrix), position, remaining_nodes, branch_width);
}

template <typename Scalar>
std::size_t certify_full_rank_by_nonzero_pivots(
    DenseScalarMatrix<Scalar> matrix, const std::string& context) {
  const auto size = rectangular_columns(
      matrix, "certified full-rank matrix");
  if (size == 0 || matrix.size() != size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": full-rank certificate requires a nonempty square matrix");
  const auto result = certified_full_rank_search(std::move(matrix), 0);
  if (result == FullRankProofResult::Proved) return size;
  if (result == FullRankProofResult::Ambiguous)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::AmbiguousZero,
        context +
            ": no certified nonzero pivot sequence proves full rank; remaining Acb enclosures overlap zero");
  if (result == FullRankProofResult::SearchBudgetExhausted)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::SearchBudgetExhausted,
        context +
            ": bounded certified pivot search exhausted its global node budget");
  throw MatchingArithmeticError(
      MatchingArithmeticErrorCode::SingularOrIncompleteSystem,
      context + ": the leading matrix is certifiably rank deficient");
}

template <typename Scalar>
DenseScalarMatrix<Scalar> epsilon_zero_matrix(
    const FiniteLaurentMatrix<Scalar>& matrix,
    const std::string& context) {
  const auto columns = rectangular_columns(matrix, "epsilon-zero matrix");
  DenseScalarMatrix<Scalar> result(
      matrix.size(),
      std::vector<Scalar>(columns, ScalarTraits<Scalar>::zero()));
  for (std::size_t row = 0; row < matrix.size(); ++row) {
    for (std::size_t column = 0; column < columns; ++column) {
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
FiniteLaurentVector<Scalar> left_multiply_dense_by_finite_vector(
    const DenseScalarMatrix<Scalar>& left,
    const FiniteLaurentVector<Scalar>& right,
    const std::string& context) {
  const auto inner = rectangular_columns(left, "dense left multiplier");
  if (inner != right.size())
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": dense multiplier and finite vector dimensions disagree");
  FiniteLaurentVector<Scalar> result;
  result.reserve(left.size());
  for (std::size_t row = 0; row < left.size(); ++row) {
    std::optional<EpsilonFrame<Scalar>> value;
    for (std::size_t column = 0; column < inner; ++column) {
      if (ScalarTraits<Scalar>::is_zero(left[row][column])) continue;
      auto term = right[column].scaled(left[row][column]);
      value = value.has_value() ? *value + term : std::move(term);
    }
    if (!value.has_value())
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::StructurallySingularTransformation,
          context + ": dense multiplier has a structurally zero row", row);
    result.push_back(std::move(*value));
  }
  return result;
}

template <typename Scalar>
DenseScalarMatrix<Scalar> multiply_dense_matrices(
    const DenseScalarMatrix<Scalar>& left,
    const DenseScalarMatrix<Scalar>& right,
    const std::string& context) {
  const auto inner = rectangular_columns(left, "dense left matrix");
  const auto columns = rectangular_columns(right, "dense right matrix");
  if (inner != right.size())
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": dense matrix dimensions disagree");
  DenseScalarMatrix<Scalar> result(
      left.size(),
      std::vector<Scalar>(columns, ScalarTraits<Scalar>::zero()));
  for (std::size_t row = 0; row < left.size(); ++row)
    for (std::size_t column = 0; column < columns; ++column)
      for (std::size_t k = 0; k < inner; ++k)
        result[row][column] += left[row][k] * right[k][column];
  return result;
}

template <typename Scalar>
FiniteLaurentMatrix<Scalar> left_multiply_dense_by_finite_matrix(
    const DenseScalarMatrix<Scalar>& left,
    const FiniteLaurentMatrix<Scalar>& right,
    const std::string& context) {
  const auto inner = rectangular_columns(left, "dense left multiplier");
  const auto columns = rectangular_columns(right, "finite right matrix");
  if (inner != right.size())
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": dense and finite matrix dimensions disagree");
  FiniteLaurentMatrix<Scalar> result(
      left.size(), FiniteLaurentVector<Scalar>());
  for (std::size_t row = 0; row < left.size(); ++row) {
    result[row].reserve(columns);
    for (std::size_t column = 0; column < columns; ++column) {
      std::optional<EpsilonFrame<Scalar>> value;
      for (std::size_t k = 0; k < inner; ++k) {
        if (ScalarTraits<Scalar>::is_zero(left[row][k])) continue;
        auto term = right[k][column].scaled(left[row][k]);
        value = value.has_value() ? *value + term : std::move(term);
      }
      if (!value.has_value())
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::StructurallySingularTransformation,
            context + ": dense multiplier has a structurally zero row", row,
            column);
      result[row].push_back(std::move(*value));
    }
  }
  return result;
}

// Return a fixed zero-radius dyadic approximation to the inverse of the
// midpoint matrix.  The inverse computation itself may be rounded: taking
// its midpoint simply chooses a concrete preconditioner.  Subsequent Acb
// multiplication encloses R*A and R*b, and factorization still requires
// every actual pivot to exclude zero.  Therefore this improves wrapping
// without licensing any midpoint zero/nonzero decision.
inline std::optional<DenseScalarMatrix<ComplexBall>>
acb_midpoint_inverse_preconditioner(
    const DenseScalarMatrix<ComplexBall>& matrix) {
  const auto size = rectangular_columns(matrix, "Acb midpoint matrix");
  if (size == 0 || matrix.size() != size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        "Acb midpoint preconditioner requires a nonempty square matrix");

  struct AcbMatrixOwner {
    explicit AcbMatrixOwner(slong size) { acb_mat_init(value, size, size); }
    ~AcbMatrixOwner() { acb_mat_clear(value); }
    acb_mat_t value;
  } midpoint(static_cast<slong>(size)), inverse(static_cast<slong>(size));
  for (std::size_t row = 0; row < size; ++row)
    for (std::size_t column = 0; column < size; ++column)
      acb_get_mid(
          acb_mat_entry(midpoint.value, static_cast<slong>(row),
                        static_cast<slong>(column)),
          matrix[row][column].raw());
  if (acb_mat_inv(inverse.value, midpoint.value,
                  ComplexBall::precision()) == 0)
    return std::nullopt;

  DenseScalarMatrix<ComplexBall> result(
      size, std::vector<ComplexBall>(size));
  for (std::size_t row = 0; row < size; ++row)
    for (std::size_t column = 0; column < size; ++column)
      acb_get_mid(
          result[row][column].raw(),
          acb_mat_entry(inverse.value, static_cast<slong>(row),
                        static_cast<slong>(column)));
  return result;
}

template <typename Scalar>
EpsilonFrame<Scalar> saturation_divide(
    const FiniteLaurentVector<Scalar>& row,
    const std::vector<Scalar>& relation,
    const std::string& context,
    std::optional<int> certified_candidate_chop_digits = std::nullopt) {
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
  if constexpr (std::is_same_v<Scalar, ComplexBall>) {
    if (certified_candidate_chop_digits.has_value()) {
      auto coefficients = combination->coefficients();
      for (auto& coefficient : coefficients)
        coefficient = ScalarTraits<ComplexBall>::canonicalized(
            coefficient, *certified_candidate_chop_digits);
      combination = EpsilonFrame<ComplexBall>(
          combination->window(), std::move(coefficients));
    }
  }
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
EpsilonLatticeSaturationResult<Scalar>
saturate_rectangular_finite_laurent_basis(
    const FiniteLaurentMatrix<Scalar>& basis,
    const std::string& context =
        "rectangular finite Laurent basis saturation",
    std::optional<int> certified_candidate_chop_digits = std::nullopt) {
  const auto size = matching_detail::rectangular_columns(
      basis, "finite Laurent saturation basis");
  const auto row_count = basis.size();
  if (row_count < size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context +
            ": saturation basis must have at least as many rows as columns");

  auto transformed = basis;
  for (std::size_t row = 0; row < row_count; ++row)
    for (std::size_t column = 0; column < size; ++column)
      transformed[row][column] = matching_detail::canonical_leading_frame(
          transformed[row][column], context + ": input valuation", row,
          column);

  std::vector<std::int32_t> valuations(size);
  std::vector<std::int32_t> shifts(size);
  auto transformation = identity_exact_laurent_matrix<Scalar>(size);
  for (std::size_t column = 0; column < size; ++column) {
    std::optional<std::int32_t> valuation;
    for (std::size_t row = 0; row < row_count; ++row) {
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
    for (std::size_t row = 0; row < row_count; ++row) {
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

  const auto initial_normalized = transformed;

  // Every elementary action contains the old target column with coefficient
  // one and then divides it by epsilon.  Consequently the target's common
  // complete top drops by at least one.  The sum below is therefore a strict,
  // finite progress measure independent of any determinant-minor guess.
  std::size_t maximum_steps = 0;
  for (std::size_t column = 0; column < size; ++column) {
    std::int32_t common_top = std::numeric_limits<std::int32_t>::max();
    for (std::size_t row = 0; row < row_count; ++row)
      common_top = std::min(common_top,
                            transformed[row][column].complete_max());
    if (common_top < 0)
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::InsufficientCompleteWindow,
          context + ": normalized column is incomplete through epsilon^0",
          std::nullopt, column, common_top);
    maximum_steps += static_cast<std::size_t>(common_top) + 1;
  }

  const auto certified_saturation_relation =
      [&](const matching_detail::DenseScalarMatrix<Scalar>& leading,
          const std::string& relation_context) {
        if constexpr (std::is_same_v<Scalar, ComplexBall>) {
          // A square Acb leading frame can be proved full rank without
          // deciding arbitrary off-pivot balls.  Saturation needs a null
          // relation only when the frame is actually deficient; requiring
          // every Schur remainder to exclude or equal zero rejects regular
          // interval matrices through harmless dependency loss.
          bool full_rank = row_count == size &&
              matching_detail::certified_full_rank_plan(leading, 0)
                  .has_value();
          if (!full_rank && row_count == size) {
            const auto preconditioner =
                matching_detail::acb_midpoint_inverse_preconditioner(
                    leading);
            if (preconditioner.has_value())
              full_rank = matching_detail::certified_full_rank_plan(
                  matching_detail::multiply_dense_matrices(
                      *preconditioner, leading,
                      relation_context +
                          ": midpoint leading normalization"),
                  0).has_value();
          }
          if (full_rank) {
            std::vector<std::size_t> pivot_rows(size);
            for (std::size_t row = 0; row < size; ++row)
              pivot_rows[row] = row;
            return matching_detail::LeadingNullRelation<Scalar>{
                size, std::nullopt, std::nullopt, std::move(pivot_rows)};
          }
        }
        return matching_detail::leading_null_relation(
            leading, relation_context,
            certified_candidate_chop_digits);
      };
  std::vector<EpsilonLatticeSaturationAction<Scalar>> actions;
  auto relation = certified_saturation_relation(
      matching_detail::epsilon_zero_matrix(transformed, context),
      context + "#sat0");
  const auto initial_rank = relation.rank;
  std::size_t steps = 0;
  while (relation.rank != size) {
    if (steps >= maximum_steps)
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::SaturationFailure,
          context +
              ": complete epsilon windows were consumed before full rank");
    if (!relation.vector.has_value() ||
        !relation.target_column.has_value())
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::SaturationFailure,
          context + ": deficient leading frame has no null relation");

    const auto& null_vector = *relation.vector;
    const auto target = *relation.target_column;
    for (std::size_t row = 0; row < row_count; ++row)
      transformed[row][target] = matching_detail::saturation_divide(
          transformed[row], null_vector,
          context + "#sat" + std::to_string(steps + 1) + "/row" +
              std::to_string(row),
          certified_candidate_chop_digits);

    auto elementary = identity_exact_laurent_matrix<Scalar>(size);
    for (std::size_t row = 0; row < size; ++row)
      elementary[row][target] =
          ExactLaurentPolynomial<Scalar>::monomial(-1, null_vector[row]);
    transformation = multiply_exact_laurent_matrices(transformation,
                                                       elementary);
    actions.push_back(
        {relation.rank, target, std::move(*relation.vector)});
    ++steps;
    relation = certified_saturation_relation(
        matching_detail::epsilon_zero_matrix(transformed, context),
        context + "#sat" + std::to_string(steps));
  }

  if (relation.pivot_rows.size() != size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InvalidSaturationLattice,
        context + ": final full-rank proof did not retain all pivot rows");

  FiniteLaurentMatrix<Scalar> determinant_minor;
  determinant_minor.reserve(size);
  for (const auto row : relation.pivot_rows)
    determinant_minor.push_back(initial_normalized[row]);
  const auto determinant = [&]() {
    if constexpr (std::is_same_v<Scalar, Rational>) {
      const auto action_valuation = matching_detail::checked_power(
          static_cast<std::int64_t>(steps),
          "exact saturation action determinant valuation");
      if (!relation.full_rank_leading_determinant.has_value() ||
          relation.full_rank_leading_determinant->is_zero())
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InvalidSaturationLattice,
            context +
                ": exact saturation reached full leading rank without a nonzero leading determinant");

      // Every recorded elementary action replaces one target column by a
      // null combination divided by epsilon, with exact target coefficient
      // one.  Its determinant is therefore epsilon^-1.  The already-proved
      // full-rank epsilon^0 leading minor and the action count prove that the
      // same minor of the initial normalized basis has valuation `steps`;
      // its leading coefficient is the pivot product retained by that exact
      // rank proof.  Expanding any positive-epsilon determinant coefficient
      // is irrelevant to this certificate and can cause catastrophic
      // rational coefficient swell.
      return EpsilonFrame<Scalar>(
          {action_valuation, action_valuation},
          {*relation.full_rank_leading_determinant});
    } else {
      return matching_detail::nonnegative_determinant_frame(
          determinant_minor, context + ": pivot-row determinant witness",
          certified_candidate_chop_digits);
    }
  }();
  const auto determinant_valuation = finite_laurent_leading_power(
      determinant, context + ": normalized pivot-row determinant");
  if (!determinant_valuation.has_value() || *determinant_valuation < 0)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::UnresolvedDeterminantTail,
        context +
            ": pivot-row determinant is unresolved in the complete epsilon window",
        std::nullopt, std::nullopt, determinant.complete_max());
  if (steps != static_cast<std::size_t>(*determinant_valuation))
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::SaturationFailure,
        context +
            ": exact pivot-row determinant valuation disagrees with the number of epsilon divisions");

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
EpsilonLatticeSaturationResult<Scalar> saturate_finite_laurent_basis(
    const FiniteLaurentMatrix<Scalar>& basis,
    const std::string& context = "finite Laurent basis saturation",
    std::optional<int> certified_candidate_chop_digits = std::nullopt) {
  const auto size = matching_detail::rectangular_columns(
      basis, "finite Laurent saturation basis");
  if (basis.size() != size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": saturation basis must be square");
  return saturate_rectangular_finite_laurent_basis(
      basis, context, certified_candidate_chop_digits);
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

template <typename Scalar>
struct FiniteLaurentFactorization {
  struct Elimination {
    std::size_t row = 0;
    EpsilonFrame<Scalar> factor;
  };

  struct Step {
    std::size_t row_swap = 0;
    std::size_t column_swap = 0;
    std::vector<Elimination> eliminations;
  };

  FiniteLaurentMatrix<Scalar> upper;
  std::vector<std::size_t> column_permutation;
  std::vector<Step> steps;
  // A fixed constant row transformation chosen before elimination.  It is
  // applied to every right-hand side before replaying `steps`.  For Acb this
  // may be a zero-radius dyadic approximate inverse used solely to reduce
  // interval wrapping; no entry of it is treated as an exact inverse.
  std::optional<matching_detail::DenseScalarMatrix<Scalar>>
      left_preconditioner;
  std::string preconditioner_kind = "none";
};

// Deterministic full-pivot Gaussian elimination over finite Laurent frames.
// Pivots are ordered by certified Laurent valuation; ties use row/column
// order.  For Acb, a candidate leading coefficient must exclude zero.  An
// enclosure overlapping zero is never classified from its midpoint and is
// rejected loudly.  The recorded row operations can be replayed for several
// right-hand sides, which is the factorization contract needed by iterative
// refinement.
template <typename Scalar>
FiniteLaurentFactorization<Scalar> factor_finite_laurent_system(
    FiniteLaurentMatrix<Scalar> matrix,
    const std::string& context = "finite Laurent factorization",
    std::optional<int> certified_candidate_chop_digits = std::nullopt) {
  const auto size = matching_detail::rectangular_columns(
      matrix, "finite Laurent coefficient matrix");
  if (matrix.size() != size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": coefficient matrix must be square");

  std::vector<std::size_t> column_permutation(size);
  for (std::size_t i = 0; i < size; ++i) column_permutation[i] = i;
  std::vector<typename FiniteLaurentFactorization<Scalar>::Step> steps;
  steps.reserve(size);

  for (std::size_t position = 0; position < size; ++position) {
    // A numeric match is ultimately certified by a residual against the
    // untouched original basis.  It is therefore safe to discard an Acb
    // enclosure that lies wholly below a much stricter candidate floor while
    // choosing the factorization: the discarded ball is still present in the
    // final residual.  This prevents exact cancellations represented by
    // independent balls from manufacturing 10^-900-scale pseudo-pivots.
    if constexpr (std::is_same_v<Scalar, ComplexBall>) {
      if (certified_candidate_chop_digits.has_value()) {
        for (std::size_t row = position; row < size; ++row)
          for (std::size_t column = position; column < size; ++column) {
            auto coefficients = matrix[row][column].coefficients();
            for (auto& coefficient : coefficients)
              coefficient = ScalarTraits<ComplexBall>::canonicalized(
                  coefficient, *certified_candidate_chop_digits);
            matrix[row][column] = EpsilonFrame<ComplexBall>(
                matrix[row][column].window(), std::move(coefficients));
          }
      }
    }
    struct Pivot {
      std::size_t row;
      std::size_t column;
      std::int32_t power;
    };
    std::optional<Pivot> pivot;
    std::optional<std::tuple<std::size_t, std::size_t, std::int32_t>>
        first_ambiguous;
    for (std::size_t row = position; row < size; ++row) {
      for (std::size_t column = position; column < size; ++column) {
        const auto leading =
            matching_detail::certified_laurent_leading_power(
                matrix[row][column]);
        if (leading.first_ambiguous_power.has_value() &&
            !first_ambiguous.has_value())
          first_ambiguous =
              std::tuple{row, column, *leading.first_ambiguous_power};
        if (leading.power.has_value()) {
          const bool lower_valuation =
              !pivot.has_value() || *leading.power < pivot->power;
          const bool stronger_same_valuation =
              pivot.has_value() && *leading.power == pivot->power &&
              matching_detail::compare_certified_pivot_strength(
                  matrix[row][column].coefficient(*leading.power),
                  matrix[pivot->row][pivot->column].coefficient(
                      pivot->power)) > 0;
          if (lower_valuation || stronger_same_valuation)
            pivot = Pivot{row, column, *leading.power};
        }
      }
    }
    if (!pivot.has_value() && first_ambiguous.has_value()) {
      const auto [row, column, power] = *first_ambiguous;
      std::string ball_detail;
      if constexpr (std::is_same_v<Scalar, ComplexBall>) {
        const auto& coefficient = matrix[row][column].coefficient(power);
        ball_detail = " at elimination position " +
                      std::to_string(position) + ", row " +
                      std::to_string(row) + ", column " +
                      std::to_string(column) + ", epsilon^" +
                      std::to_string(power) + ", midpoint=(" +
                      coefficient.real_midpoint(12) + "," +
                      coefficient.imag_midpoint(12) + "), radius2exp=(" +
                      coefficient.real_radius_exponent() + "," +
                      coefficient.imag_radius_exponent() + ")";
      }
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::AmbiguousZero,
          context +
              ": no certified Laurent pivot remains; an Acb candidate "
              "overlaps zero at its required leading coefficient" +
              ball_detail,
          row, column, power);
    }
    if (!pivot.has_value())
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::SingularOrIncompleteSystem,
          context +
              ": no proved nonzero pivot remains in the complete epsilon "
              "windows",
          position, position);

    typename FiniteLaurentFactorization<Scalar>::Step step;
    step.row_swap = pivot->row;
    step.column_swap = pivot->column;
    if (pivot->row != position) {
      std::swap(matrix[position], matrix[pivot->row]);
    }
    if (pivot->column != position) {
      for (auto& row : matrix)
        std::swap(row[position], row[pivot->column]);
      std::swap(column_permutation[position],
                column_permutation[pivot->column]);
    }
    // Only the selected pivot needs a canonical leading frame.  Its probe
    // above certified every discarded coefficient as exact zero.
    matrix[position][position] = matching_detail::canonical_leading_frame(
        matrix[position][position], context + ": selected pivot", position,
        position);
    for (std::size_t column = position + 1; column < size; ++column) {
      const auto leading =
          matching_detail::certified_laurent_leading_power(
              matrix[position][column]);
      if (leading.power.has_value()) {
        // These entries are about to be multiplied by every elimination
        // factor.  Retaining certified leading zeros would unnecessarily
        // lower the honest product CompleteMax, just as it does for the
        // below-pivot quotient numerator.
        matrix[position][column] =
            matching_detail::canonical_leading_frame(
                matrix[position][column],
                context + ": certified pivot-row entry", position, column);
      } else if (!leading.first_ambiguous_power.has_value()) {
        matrix[position][column] = EpsilonFrame<Scalar>::zero(
            matrix[position][column].complete_max());
      }
    }

    for (std::size_t row = position + 1; row < size; ++row) {
      const auto below_pivot =
          matching_detail::certified_laurent_leading_power(
              matrix[row][position]);
      if (!below_pivot.power.has_value() &&
          !below_pivot.first_ambiguous_power.has_value()) {
        matrix[row][position] = EpsilonFrame<Scalar>::zero(
            matrix[row][position].complete_max());
      }
      if (below_pivot.power.has_value())
        matrix[row][position] =
            matching_detail::canonical_leading_frame(
                matrix[row][position],
                context + ": certified below-pivot numerator", row,
                position);
      const auto factor = finite_laurent_quotient(
          matrix[row][position], matrix[position][position],
          context + ": elimination quotient");
      // This entry is zero by the defining quotient recurrence, coefficient
      // by coefficient.  Re-evaluating A - (A/pivot)*pivot with independent
      // Acb balls only destroys that dependency and manufactures an
      // ambiguous enclosure around zero.  Retain exactly the complete top of
      // the finite cancellation, then record the proved eliminated entry as
      // structural zero.  Other Schur entries do not have this identity and
      // still pass through certified leading-coefficient decisions below.
      const auto eliminated =
          matrix[row][position] - factor * matrix[position][position];
      matrix[row][position] =
          EpsilonFrame<Scalar>::zero(eliminated.complete_max());
      for (std::size_t column = position + 1; column < size; ++column) {
        matrix[row][column] =
            matrix[row][column] - factor * matrix[position][column];
      }
      step.eliminations.push_back({row, factor});
    }
    steps.push_back(std::move(step));
  }

  for (std::size_t row = 0; row < size; ++row)
    for (std::size_t column = row; column < size; ++column)
      matrix[row][column] =
          matching_detail::canonicalize_certified_or_preserve_ambiguous(
              matrix[row][column], context + ": stored upper triangle", row,
              column);

  return {std::move(matrix), std::move(column_permutation),
          std::move(steps)};
}

// A constant midpoint-derived row preconditioner is safe for Acb Laurent
// systems: it is applied to the complete matrix and every right-hand side,
// while all subsequent pivots and the final residual remain ball-certified.
// Exact lattice saturation normally makes the epsilon^0 leading matrix
// nonsingular.  Normalize that matrix first: this is the finite-Laurent
// analogue of an epsilon-Fuchsian leading frame and prevents the formal
// coefficient solve from losing the condition number again at every order.
// General callers whose leading matrix is singular retain the older dyadic
// sample fallbacks.
inline FiniteLaurentFactorization<ComplexBall>
factor_preconditioned_acb_finite_laurent_system(
    const FiniteLaurentMatrix<ComplexBall>& matrix,
    const std::string& context =
        "preconditioned Acb finite Laurent factorization") {
  const auto size = matching_detail::rectangular_columns(
      matrix, "Acb Laurent preconditioner matrix");
  if (matrix.size() != size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": coefficient matrix must be square");

  std::string candidate_failures;
  bool complete_at_zero = true;
  matching_detail::DenseScalarMatrix<ComplexBall> leading(
      size, std::vector<ComplexBall>(size));
  for (std::size_t row = 0; row < size; ++row)
    for (std::size_t column = 0; column < size; ++column) {
      const auto& frame = matrix[row][column];
      if (frame.complete_max() < 0) {
        complete_at_zero = false;
        continue;
      }
      leading[row][column] = frame.coefficient(0);
    }
  if (complete_at_zero) {
    auto preconditioner =
        matching_detail::acb_midpoint_inverse_preconditioner(leading);
    if (preconditioner.has_value()) {
      try {
        auto factorization = factor_finite_laurent_system(
            matching_detail::left_multiply_dense_by_finite_matrix(
                *preconditioner, matrix,
                context + ":epsilon-zero-leading-normalization"),
            context + ":certified-epsilon-zero-leading-normalization", 100);
        if (factorization.left_preconditioner.has_value())
          throw std::logic_error(
              "general Laurent factorization unexpectedly retained a preconditioner");
        factorization.left_preconditioner = std::move(preconditioner);
        factorization.preconditioner_kind = "epsilon-zero-leading";
        return factorization;
      } catch (const MatchingArithmeticError& error) {
        candidate_failures +=
            " [epsilon-zero leading normalization: " +
            std::string(error.what()) + "]";
      }
    }
  }
  for (const auto dyadic_exponent : {1, 2, 3, 4}) {
    const auto point = ComplexBall::from_strings(
        "1/" + std::to_string(1U << dyadic_exponent));
    matching_detail::DenseScalarMatrix<ComplexBall> evaluated(
        size, std::vector<ComplexBall>(size));
    for (std::size_t row = 0; row < size; ++row)
      for (std::size_t column = 0; column < size; ++column) {
        const auto& frame = matrix[row][column];
        auto total = ScalarTraits<ComplexBall>::zero();
        for (std::int64_t power = frame.min_power();
             power <= frame.complete_max(); ++power) {
          auto monomial = ScalarTraits<ComplexBall>::one();
          if (power >= 0) {
            for (std::int64_t step = 0; step < power; ++step)
              monomial *= point;
          } else {
            for (std::int64_t step = 0; step > power; --step)
              monomial = monomial / point;
          }
          total += frame.coefficient(static_cast<std::int32_t>(power)) *
                   monomial;
        }
        evaluated[row][column] = std::move(total);
      }
    auto preconditioner =
        matching_detail::acb_midpoint_inverse_preconditioner(evaluated);
    if (!preconditioner.has_value()) continue;
    try {
      auto factorization = factor_finite_laurent_system(
          matching_detail::left_multiply_dense_by_finite_matrix(
              *preconditioner, matrix,
              context + ":dyadic-" + std::to_string(dyadic_exponent)),
          context + ":certified-dyadic-" +
              std::to_string(dyadic_exponent),
          100);
      if (factorization.left_preconditioner.has_value())
        throw std::logic_error(
            "general Laurent factorization unexpectedly retained a preconditioner");
      factorization.left_preconditioner = std::move(preconditioner);
      factorization.preconditioner_kind =
          "dyadic-" + std::to_string(dyadic_exponent);
      return factorization;
    } catch (const MatchingArithmeticError& error) {
      // A proposed midpoint preconditioner carries no proof obligation.
      // Try the next dyadic sample; only the returned certified factorization
      // is authoritative.
      candidate_failures += " [dyadic " +
                            std::to_string(dyadic_exponent) + ": " +
                            error.what() + "]";
    }
  }
  try {
    return factor_finite_laurent_system(
        matrix, context + ":unpreconditioned", 100);
  } catch (const MatchingArithmeticError& error) {
    throw MatchingArithmeticError(
        error.code, std::string(error.what()) +
                        "; rejected midpoint preconditioners:" +
                        candidate_failures,
        error.row, error.column, error.epsilon_power);
  }
}

// Factor a basis whose exact provenance proves that every negative epsilon
// coefficient vanishes.  Rebuild every entry at epsilon^0 before invoking
// the ordinary certified Laurent factorization.  The latter deliberately
// remains free to choose positive-valuation pivots: singular perturbations
// need not have a full-rank epsilon^0 value at a fixed matching point.
template <typename Scalar>
FiniteLaurentMatrix<Scalar>
exact_nonnegative_finite_laurent_matrix(
    FiniteLaurentMatrix<Scalar> matrix,
    const std::string& context =
        "exact nonnegative finite Laurent matrix") {
  const auto size = matching_detail::rectangular_columns(
      matrix, "exact valuation-zero coefficient matrix");
  if (matrix.size() != size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": coefficient matrix must be square");

  for (std::size_t row = 0; row < size; ++row) {
    for (std::size_t column = 0; column < size; ++column) {
      const auto& frame = matrix[row][column];
      if (frame.complete_max() < 0)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            context + ": transformed basis is incomplete through epsilon^0",
            row, column, frame.complete_max());
      std::vector<Scalar> nonnegative;
      nonnegative.reserve(
          static_cast<std::size_t>(frame.complete_max()) + 1);
      for (std::int64_t power = 0; power <= frame.complete_max(); ++power)
        nonnegative.push_back(frame.coefficient(
            static_cast<std::int32_t>(power)));
      matrix[row][column] = EpsilonFrame<Scalar>(0,
                                                  std::move(nonnegative));
    }
  }
  return matrix;
}

template <typename Scalar>
FiniteLaurentFactorization<Scalar>
factor_exact_nonnegative_finite_laurent_system(
    FiniteLaurentMatrix<Scalar> matrix,
    const std::string& context =
        "exact valuation-zero finite Laurent factorization") {
  matrix = exact_nonnegative_finite_laurent_matrix(
      std::move(matrix), context + ": certified nonnegative matrix");
  const auto size = matrix.size();

  // When the epsilon^0 value already has a certified full-pivot plan, replay
  // that plan on the complete Laurent frames.  The dense certificate permits
  // ambiguous off-pivot balls: they are arbitrary entries, not zero claims.
  // This matters for an exact Rational-shadow CASE-P basis, where Schur
  // subtraction can enclose an algebraically exact zero even at arbitrarily
  // high Acb precision.  Re-canonicalizing every Schur entry would demand an
  // unnecessary zero decision and reject a matrix whose nonzero pivot path
  // already proves invertibility.
  auto pivot_plan = matching_detail::certified_full_rank_plan(
      matching_detail::epsilon_zero_matrix(matrix, context), 0);
  std::optional<matching_detail::DenseScalarMatrix<Scalar>>
      left_preconditioner;
  if (!pivot_plan.has_value()) {
    if constexpr (std::is_same_v<Scalar, ComplexBall>) {
      // Exact CASE-P provenance proves that the transformed formal basis has
      // no negative epsilon coefficients and is a full-rank epsilon lattice.
      // Numerical evaluation can nevertheless make ordinary interval
      // elimination wrap badly.  A fixed midpoint-inverse row
      // preconditioner is safe: it is applied to A and every rhs, while the
      // resulting pivots remain subject to the usual Acb exclusion test.
      auto candidate = matching_detail::acb_midpoint_inverse_preconditioner(
          matching_detail::epsilon_zero_matrix(matrix, context));
      if (candidate.has_value()) {
        auto preconditioned =
            matching_detail::left_multiply_dense_by_finite_matrix(
                *candidate, matrix, context + ": midpoint preconditioning");
        auto candidate_plan = matching_detail::certified_full_rank_plan(
            matching_detail::epsilon_zero_matrix(preconditioned, context),
            0);
        if (candidate_plan.has_value()) {
          matrix = std::move(preconditioned);
          pivot_plan = std::move(candidate_plan);
          left_preconditioner = std::move(candidate);
        }
      }
    }
  }
  if (pivot_plan.has_value()) {
    std::vector<std::size_t> column_permutation(size);
    for (std::size_t index = 0; index < size; ++index)
      column_permutation[index] = index;
    std::vector<typename FiniteLaurentFactorization<Scalar>::Step> steps;
    steps.reserve(size);
    for (std::size_t position = 0; position < size; ++position) {
      const auto [pivot_row, pivot_column] = (*pivot_plan)[position];
      typename FiniteLaurentFactorization<Scalar>::Step step;
      step.row_swap = pivot_row;
      step.column_swap = pivot_column;
      if (pivot_row != position)
        std::swap(matrix[position], matrix[pivot_row]);
      if (pivot_column != position) {
        for (auto& row : matrix)
          std::swap(row[position], row[pivot_column]);
        std::swap(column_permutation[position],
                  column_permutation[pivot_column]);
      }
      // The epsilon^0 dense plan proves this pivot excludes zero.  The
      // finite quotient therefore has a certified denominator even when a
      // below-pivot numerator overlaps zero.
      (void)matching_detail::trim_leading_certified_zeros(
          matrix[position][position],
          context + ": planned pivot", position, position);
      for (std::size_t row = position + 1; row < size; ++row) {
        matrix[row][position] =
            matching_detail::canonicalize_certified_or_preserve_ambiguous(
                matrix[row][position],
                context + ": planned below-pivot numerator", row,
                position);
        const auto factor = finite_laurent_quotient(
            matrix[row][position], matrix[position][position],
            context + ": planned elimination quotient");
        const auto eliminated =
            matrix[row][position] - factor * matrix[position][position];
        matrix[row][position] =
            EpsilonFrame<Scalar>::zero(eliminated.complete_max());
        for (std::size_t column = position + 1; column < size; ++column)
          matrix[row][column] =
              matrix[row][column] - factor * matrix[position][column];
        step.eliminations.push_back({row, factor});
      }
      steps.push_back(std::move(step));
    }
    for (std::size_t row = 0; row < size; ++row)
      for (std::size_t column = row; column < size; ++column)
        matrix[row][column] =
            matching_detail::canonicalize_certified_or_preserve_ambiguous(
                matrix[row][column], context + ": stored upper triangle",
                row, column);
    const auto preconditioner_kind = left_preconditioner.has_value()
        ? std::string("epsilon-zero-leading")
        : std::string("none");
    return {std::move(matrix), std::move(column_permutation),
            std::move(steps), std::move(left_preconditioner),
            preconditioner_kind};
  }

  return factor_finite_laurent_system(
      std::move(matrix), context + ": certified Laurent pivots");
}

// A saturated basis with a proved invertible epsilon^0 coefficient is a
// matrix power series, not a generic Laurent matrix.  Eliminating whole
// finite series loses an upper coefficient at every pivot because every
// Schur update multiplies two independently truncated frames.  The formal
// coefficient recurrence below factors the constant matrix once and loses
// only the genuinely unavailable matrix tail, independent of dimension.
template <typename Scalar>
struct ExactNonnegativePowerSeriesFactorization {
  FiniteLaurentMatrix<Scalar> matrix;
  FiniteLaurentFactorization<Scalar> leading;
  std::int32_t common_complete_max = 0;
};

template <typename Scalar>
ExactNonnegativePowerSeriesFactorization<Scalar>
factor_exact_nonnegative_power_series_system(
    FiniteLaurentMatrix<Scalar> matrix,
    const std::string& context =
        "exact nonnegative power-series factorization") {
  matrix = exact_nonnegative_finite_laurent_matrix(
      std::move(matrix), context + ": certified nonnegative matrix");
  const auto size = matrix.size();
  std::int32_t common_complete_max =
      std::numeric_limits<std::int32_t>::max();
  FiniteLaurentMatrix<Scalar> leading(
      size, FiniteLaurentVector<Scalar>());
  for (std::size_t row = 0; row < size; ++row) {
    leading[row].reserve(size);
    for (std::size_t column = 0; column < size; ++column) {
      common_complete_max = std::min(
          common_complete_max, matrix[row][column].complete_max());
      leading[row].emplace_back(
          0, std::vector<Scalar>{matrix[row][column].coefficient(0)});
    }
  }
  auto leading_factorization = [&]() {
    if constexpr (std::is_same_v<Scalar, ComplexBall>)
      return factor_preconditioned_acb_finite_laurent_system(
          leading, context + ": certified epsilon-zero factorization");
    else
      return factor_finite_laurent_system(
          std::move(leading),
          context + ": exact epsilon-zero factorization");
  }();
  leading_factorization.preconditioner_kind =
      leading_factorization.preconditioner_kind == "none"
      ? "epsilon-coefficient-recurrence"
      : "epsilon-coefficient-recurrence/" +
            leading_factorization.preconditioner_kind;
  return {std::move(matrix), std::move(leading_factorization),
          common_complete_max};
}

template <typename Scalar>
FiniteLaurentVector<Scalar> solve_factorized_finite_laurent_system(
    const FiniteLaurentFactorization<Scalar>& factorization,
    FiniteLaurentVector<Scalar> right_hand_side,
    const std::string& context = "factorized finite Laurent solve") {
  const auto size = factorization.upper.size();
  if (size == 0 || right_hand_side.size() != size ||
      factorization.column_permutation.size() != size ||
      factorization.steps.size() != size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": factorization and rhs dimensions disagree");

  if (factorization.left_preconditioner.has_value())
    right_hand_side = matching_detail::left_multiply_dense_by_finite_vector(
        *factorization.left_preconditioner, right_hand_side,
        context + ": left preconditioning");

  for (std::size_t position = 0; position < size; ++position) {
    const auto& step = factorization.steps[position];
    if (step.row_swap >= size || step.column_swap >= size)
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::DimensionMismatch,
          context + ": malformed factorization permutation");
    if (step.row_swap != position)
      std::swap(right_hand_side[position], right_hand_side[step.row_swap]);
    right_hand_side[position] =
        matching_detail::canonicalize_certified_or_preserve_ambiguous(
            right_hand_side[position],
            context + ": forward rhs pivot row", position);
    for (const auto& elimination : step.eliminations) {
      if (elimination.row <= position || elimination.row >= size)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::DimensionMismatch,
            context + ": malformed factorization elimination row");
      right_hand_side[elimination.row] =
          right_hand_side[elimination.row] -
          elimination.factor * right_hand_side[position];
    }
  }

  FiniteLaurentVector<Scalar> permuted_solution(size,
                                                right_hand_side.front());
  for (std::size_t reverse = size; reverse-- > 0;) {
    auto residual = right_hand_side[reverse];
    for (std::size_t column = reverse + 1; column < size; ++column)
      residual = residual - factorization.upper[reverse][column] *
                                permuted_solution[column];
    residual =
        matching_detail::canonicalize_certified_or_preserve_ambiguous(
            residual, context + ": back-substitution numerator", reverse);
    permuted_solution[reverse] = finite_laurent_quotient(
        residual, factorization.upper[reverse][reverse],
        context + ": back-substitution quotient");
  }

  FiniteLaurentVector<Scalar> solution(size, permuted_solution.front());
  for (std::size_t current = 0; current < size; ++current)
    solution[factorization.column_permutation[current]] =
        std::move(permuted_solution[current]);
  return solution;
}

template <typename Scalar>
FiniteLaurentVector<Scalar>
solve_factorized_exact_nonnegative_power_series_system(
    const ExactNonnegativePowerSeriesFactorization<Scalar>& factorization,
    const FiniteLaurentVector<Scalar>& right_hand_side,
    const std::string& context =
        "factorized exact nonnegative power-series solve") {
  const auto size = factorization.matrix.size();
  if (size == 0 || right_hand_side.size() != size ||
      factorization.leading.upper.size() != size)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": factorization and rhs dimensions disagree");

  auto solution_min = right_hand_side.front().min_power();
  auto rhs_complete_max = right_hand_side.front().complete_max();
  for (const auto& row : right_hand_side) {
    solution_min = std::min(solution_min, row.min_power());
    rhs_complete_max = std::min(rhs_complete_max, row.complete_max());
  }
  const auto matrix_limited_max = matching_detail::checked_power(
      static_cast<std::int64_t>(factorization.common_complete_max) +
          solution_min,
      "power-series solution complete maximum");
  const auto solution_complete_max =
      std::min(rhs_complete_max, matrix_limited_max);
  if (solution_complete_max < solution_min)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InsufficientCompleteWindow,
        context +
            ": the matrix tail does not cover one complete solution coefficient",
        std::nullopt, std::nullopt, solution_complete_max);

  const auto width = EpsilonWindow{
      solution_min, solution_complete_max}.width();
  std::vector<std::vector<Scalar>> coefficients(
      size, std::vector<Scalar>(width, ScalarTraits<Scalar>::zero()));
  for (std::int64_t raw_power = solution_min;
       raw_power <= solution_complete_max; ++raw_power) {
    const auto power = static_cast<std::int32_t>(raw_power);
    const auto output_index = static_cast<std::size_t>(
        raw_power - solution_min);
    FiniteLaurentVector<Scalar> coefficient_rhs;
    coefficient_rhs.reserve(size);
    for (std::size_t row = 0; row < size; ++row) {
      auto value = right_hand_side[row].coefficient(power);
      for (std::int64_t shift = 1;
           shift <= raw_power - solution_min; ++shift) {
        const auto matrix_power = static_cast<std::int32_t>(shift);
        const auto previous = static_cast<std::size_t>(
            raw_power - shift - solution_min);
        for (std::size_t column = 0; column < size; ++column)
          value -= factorization.matrix[row][column].coefficient(
                       matrix_power) *
                   coefficients[column][previous];
      }
      coefficient_rhs.emplace_back(
          0, std::vector<Scalar>{std::move(value)});
    }
    const auto solved = solve_factorized_finite_laurent_system(
        factorization.leading, std::move(coefficient_rhs),
        context + ": coefficient epsilon^" + std::to_string(power));
    for (std::size_t column = 0; column < size; ++column)
      coefficients[column][output_index] = solved[column].coefficient(0);
  }

  FiniteLaurentVector<Scalar> solution;
  solution.reserve(size);
  for (auto& column : coefficients)
    solution.emplace_back(
        EpsilonWindow{solution_min, solution_complete_max},
        std::move(column));
  return solution;
}

// One-shot convenience wrapper.  Refinement callers should retain the
// factorization and call solve_factorized_finite_laurent_system repeatedly.
template <typename Scalar>
FiniteLaurentVector<Scalar> solve_finite_laurent_system(
    FiniteLaurentMatrix<Scalar> matrix,
    FiniteLaurentVector<Scalar> right_hand_side,
    const std::string& context = "finite Laurent linear solve") {
  auto factorization = factor_finite_laurent_system(
      std::move(matrix), context + ": factor");
  return solve_factorized_finite_laurent_system(
      factorization, std::move(right_hand_side), context + ": solve");
}

enum class AcbMatchingResidualVerdict : std::uint8_t {
  Pass,
  Fail,
  Inconclusive
};

struct AcbMatchingCoefficientResidual {
  std::size_t row = 0;
  std::int32_t epsilon_power = 0;
  Magnitude residual_lower;
  Magnitude residual_upper;
  Magnitude scale_lower;
  Magnitude scale_upper;
  AcbMatchingResidualVerdict verdict =
      AcbMatchingResidualVerdict::Inconclusive;
};

struct AcbMatchingResidualDiagnostics {
  AcbMatchingResidualVerdict verdict =
      AcbMatchingResidualVerdict::Inconclusive;
  EpsilonWindow complete_window;
  std::int32_t required_complete_max = 0;
  bool complete_through_required = false;
  std::vector<AcbMatchingCoefficientResidual> coefficients;
  std::string detail;
};

struct AcbLaurentRefinementOptions {
  Magnitude relative_tolerance = Magnitude::decimal("1e-30");
  // Matching lives in a finite Laurent quotient.  Coefficients below this
  // floor are deliberately outside the requested calculation and must not be
  // chased by residual refinement.  When absent, retain the historical
  // behavior and certify every structurally represented lower coefficient.
  std::optional<std::int32_t> required_min_power;
  std::int32_t required_complete_max = 0;
  std::size_t max_refinement_steps = 2;
};

struct RefinedAcbLaurentMatch {
  // Coordinates in the exactly saturated basis F*T.
  FiniteLaurentVector<ComplexBall> transformed_weights;
  // Coordinates in the caller's original basis F, equal to T times the
  // transformed weights with honest finite-window loss.
  FiniteLaurentVector<ComplexBall> weights;
  FiniteLaurentVector<ComplexBall> residual;
  std::vector<AcbMatchingResidualDiagnostics> residual_history;
  std::size_t refinement_steps = 0;
  std::string factorization_preconditioner = "none";
};

namespace matching_detail {

struct AcbResidualEvaluation {
  FiniteLaurentVector<ComplexBall> residual;
  AcbMatchingResidualDiagnostics diagnostics;
};

using AcbResidualContributions =
    std::vector<std::vector<EpsilonFrame<ComplexBall>>>;

inline AcbResidualEvaluation
certify_acb_matching_residual_with_contributions(
    const FiniteLaurentMatrix<ComplexBall>& matrix,
    const FiniteLaurentVector<ComplexBall>& weights,
    const FiniteLaurentVector<ComplexBall>& right_hand_side,
    FiniteLaurentVector<ComplexBall> residual,
    AcbResidualContributions contributions,
    const AcbLaurentRefinementOptions& options,
    const std::string& context) {
  const auto columns = rectangular_columns(matrix,
                                            "Acb matching residual matrix");
  if (matrix.size() != right_hand_side.size() ||
      matrix.size() != residual.size() ||
      matrix.size() != contributions.size() ||
      columns != weights.size() ||
      std::any_of(contributions.begin(), contributions.end(),
          [columns](const auto& row) { return row.size() != columns; }))
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": residual dimensions disagree");
  if (!options.relative_tolerance.is_finite())
    throw std::invalid_argument(
        context + ": residual tolerance must be finite");

  auto common_min = residual.front().min_power();
  auto common_max = residual.front().complete_max();
  for (const auto& row : residual) {
    common_min = std::min(common_min, row.min_power());
    common_max = std::min(common_max, row.complete_max());
  }
  for (std::size_t row = 0; row < right_hand_side.size(); ++row) {
    common_max = std::min(
        common_max, right_hand_side[row].complete_max());
    for (const auto& contribution : contributions[row])
      common_max = std::min(common_max, contribution.complete_max());
  }
  const auto certified_min = options.required_min_power.has_value()
      ? std::max(common_min, *options.required_min_power)
      : common_min;
  if (options.required_min_power.has_value() &&
      *options.required_min_power > options.required_complete_max)
    throw std::invalid_argument(
        context + ": residual Laurent floor exceeds its required maximum");
  if (common_max < certified_min)
  {
    std::string windows = "[";
    for (std::size_t row = 0; row < residual.size(); ++row) {
      if (row != 0) windows += ",";
      windows += "[" + std::to_string(residual[row].min_power()) + "," +
          std::to_string(residual[row].complete_max()) + "]";
    }
    windows += "]";
    std::string matrix_windows = "[";
    for (std::size_t row = 0; row < matrix.size(); ++row) {
      if (row != 0) matrix_windows += ",";
      matrix_windows += "[";
      for (std::size_t column = 0; column < columns; ++column) {
        if (column != 0) matrix_windows += ",";
        matrix_windows += "[" +
            std::to_string(matrix[row][column].min_power()) + "," +
            std::to_string(matrix[row][column].complete_max()) + "]";
      }
      matrix_windows += "]";
    }
    matrix_windows += "]";
    std::string weight_windows = "[";
    for (std::size_t column = 0; column < weights.size(); ++column) {
      if (column != 0) weight_windows += ",";
      weight_windows += "[" +
          std::to_string(weights[column].min_power()) + "," +
          std::to_string(weights[column].complete_max()) + "]";
    }
    weight_windows += "]";
    std::string rhs_windows = "[";
    for (std::size_t row = 0; row < right_hand_side.size(); ++row) {
      if (row != 0) rhs_windows += ",";
      rhs_windows += "[" +
          std::to_string(right_hand_side[row].min_power()) + "," +
          std::to_string(right_hand_side[row].complete_max()) + "]";
    }
    rhs_windows += "]";
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InsufficientCompleteWindow,
        context + ": residual rows have no common complete window; " +
            "certified_min=" + std::to_string(certified_min) +
            "; common_complete_max=" + std::to_string(common_max) +
            "; required_complete_max=" +
            std::to_string(options.required_complete_max) +
            "; row_windows=" + windows +
            "; matrix_windows=" + matrix_windows +
            "; weight_windows=" + weight_windows +
            "; rhs_windows=" + rhs_windows,
        std::nullopt, std::nullopt, common_max);
  }

  AcbMatchingResidualDiagnostics diagnostics;
  diagnostics.complete_window = {certified_min, common_max};
  diagnostics.required_complete_max = options.required_complete_max;
  diagnostics.complete_through_required =
      common_max >= options.required_complete_max;
  diagnostics.coefficients.reserve(
      diagnostics.complete_window.width() * matrix.size());

  bool any_fail = false;
  bool all_pass = true;
  for (std::size_t row = 0; row < matrix.size(); ++row) {
    for (std::int64_t power64 = certified_min; power64 <= common_max;
         ++power64) {
      const auto power = static_cast<std::int32_t>(power64);
      const auto residual_value = residual[row].coefficient(power);
      auto scale_lower = Magnitude::one();
      auto scale_upper = Magnitude::one();
      const auto rhs_value = right_hand_side[row].coefficient(power);
      scale_lower = Magnitude::maximum(
          scale_lower, Magnitude::lower_abs(rhs_value));
      scale_upper = Magnitude::maximum(
          scale_upper, Magnitude::upper_abs(rhs_value));

      // Publication is an accuracy claim about the reconstructed physical
      // right-hand side, not merely a backward-error claim for A w = b.
      // Scaling by max_j |A_j w_j| lets an arbitrarily ill-conditioned
      // cancellation authorize an arbitrarily large physical residual.  In
      // particular, O(10^70) individual terms which should sum to O(1)
      // would label an O(10^62) residual an "8 digit" pass.  The exact
      // right/normal-frame routes exist precisely to form that cancellation
      // tightly, so certify it against max(1, |b|).  Contribution bounds are
      // still retained above to intersect every honest finite-Laurent
      // completeness window; they no longer weaken the accuracy contract.

      const auto residual_lower = Magnitude::lower_abs(residual_value);
      const auto residual_upper = Magnitude::upper_abs(residual_value);
      const auto pass_bound =
          scale_lower * options.relative_tolerance;
      const auto fail_bound =
          scale_upper * options.relative_tolerance;
      AcbMatchingResidualVerdict verdict =
          AcbMatchingResidualVerdict::Inconclusive;
      if (residual_upper <= pass_bound)
        verdict = AcbMatchingResidualVerdict::Pass;
      else if (residual_lower > fail_bound)
        verdict = AcbMatchingResidualVerdict::Fail;
      if (power <= options.required_complete_max) {
        any_fail = any_fail || verdict == AcbMatchingResidualVerdict::Fail;
        all_pass = all_pass && verdict == AcbMatchingResidualVerdict::Pass;
      }
      diagnostics.coefficients.push_back(
          {row, power, residual_lower, residual_upper, scale_lower,
           scale_upper, verdict});
    }
  }

  if (any_fail) {
    diagnostics.verdict = AcbMatchingResidualVerdict::Fail;
    diagnostics.detail =
        "a certified residual lower bound exceeds the tolerance-scaled "
        "contribution upper bound";
  } else if (all_pass && diagnostics.complete_through_required) {
    diagnostics.verdict = AcbMatchingResidualVerdict::Pass;
    diagnostics.detail =
        "every coefficient in the required complete window satisfies the "
        "certified residual bound";
  } else {
    diagnostics.verdict = AcbMatchingResidualVerdict::Inconclusive;
    diagnostics.detail = diagnostics.complete_through_required
        ? "at least one residual enclosure overlaps the requested tolerance"
        : "the honest residual window does not reach the required complete "
          "epsilon power";
  }
  return {std::move(residual), std::move(diagnostics)};
}

inline AcbResidualEvaluation evaluate_acb_matching_residual(
    const FiniteLaurentMatrix<ComplexBall>& matrix,
    const FiniteLaurentVector<ComplexBall>& weights,
    const FiniteLaurentVector<ComplexBall>& right_hand_side,
    const AcbLaurentRefinementOptions& options,
    const std::string& context) {
  const auto columns = rectangular_columns(matrix,
                                            "Acb matching residual matrix");
  if (matrix.size() != right_hand_side.size() || columns != weights.size())
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": residual dimensions disagree");

  AcbResidualContributions contributions(matrix.size());
  FiniteLaurentVector<ComplexBall> residual;
  residual.reserve(matrix.size());
  for (std::size_t row = 0; row < matrix.size(); ++row) {
    auto& row_terms = contributions[row];
    row_terms.reserve(columns);
    for (std::size_t column = 0; column < columns; ++column)
      row_terms.push_back(matrix[row][column] * weights[column]);
    auto reconstructed = row_terms.front();
    for (std::size_t column = 1; column < columns; ++column)
      reconstructed = reconstructed + row_terms[column];
    residual.push_back(right_hand_side[row] - reconstructed);
  }
  return certify_acb_matching_residual_with_contributions(
      matrix, weights, right_hand_side, std::move(residual),
      std::move(contributions), options, context);
}

// Certify a residual formed by an independently certified linear route.  This
// is essential when a well-conditioned exact normal frame produces a tight
// residual but subtracting the corresponding large physical-frame quantities
// would lose that information to interval dependency.  The physical matrix,
// weights, and right-hand side still determine the relative scale; only the
// cancellation-heavy subtraction is replaced.
inline AcbResidualEvaluation certify_precomputed_acb_matching_residual(
    const FiniteLaurentMatrix<ComplexBall>& matrix,
    const FiniteLaurentVector<ComplexBall>& weights,
    const FiniteLaurentVector<ComplexBall>& right_hand_side,
    FiniteLaurentVector<ComplexBall> residual,
    const AcbLaurentRefinementOptions& options,
    const std::string& context) {
  const auto columns = rectangular_columns(matrix,
                                            "Acb matching residual matrix");
  if (matrix.size() != right_hand_side.size() ||
      matrix.size() != residual.size() || columns != weights.size())
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": residual dimensions disagree");

  AcbResidualContributions contributions(matrix.size());
  for (std::size_t row = 0; row < matrix.size(); ++row) {
    contributions[row].reserve(columns);
    for (std::size_t column = 0; column < columns; ++column)
      contributions[row].push_back(matrix[row][column] * weights[column]);
  }
  return certify_acb_matching_residual_with_contributions(
      matrix, weights, right_hand_side, std::move(residual),
      std::move(contributions), options, context);
}

inline EpsilonFrame<ComplexBall> project_acb_frame_to_laurent_floor(
    const EpsilonFrame<ComplexBall>& frame, std::int32_t floor,
    const std::string& context) {
  if (frame.min_power() >= floor) return frame;
  if (frame.complete_max() < floor)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InsufficientCompleteWindow,
        context + ": residual correction does not reach its Laurent floor");
  std::vector<ComplexBall> coefficients;
  coefficients.reserve(EpsilonWindow{floor, frame.complete_max()}.width());
  for (std::int64_t power = floor; power <= frame.complete_max(); ++power)
    coefficients.push_back(
        frame.coefficient(static_cast<std::int32_t>(power)));
  return EpsilonFrame<ComplexBall>(
      {floor, frame.complete_max()}, std::move(coefficients));
}

inline FiniteLaurentVector<ComplexBall>
project_acb_residual_to_laurent_floor(
    const FiniteLaurentVector<ComplexBall>& residual, std::int32_t floor,
    const std::string& context) {
  FiniteLaurentVector<ComplexBall> projected;
  projected.reserve(residual.size());
  for (const auto& row : residual)
    projected.push_back(project_acb_frame_to_laurent_floor(
        row, floor, context));
  return projected;
}

// Matching weights are candidates whose authority comes from the residual,
// not from the particular numerical solve that proposed them.  A singular
// Laurent change of coordinates can consume several *upper* coefficients
// while forming F*T even though those coefficients cancel again in the
// physical weights T*y.  Extend the candidate problem with zero-radius zeros
// so the proposal solve can carry those cancellations to completion.  This
// extension must never be passed to residual certification: coefficients
// above the input CompleteMax remain mathematically unknown.
inline EpsilonFrame<ComplexBall> zero_extend_acb_match_candidate(
    const EpsilonFrame<ComplexBall>& frame, std::size_t additional,
    const std::string& context) {
  if (additional == 0) return frame;
  const auto extension_context =
      context + ": candidate upper extension";
  const auto maximum = checked_power(
      static_cast<std::int64_t>(frame.complete_max()) +
          static_cast<std::int64_t>(additional),
      extension_context.c_str());
  auto coefficients = frame.coefficients();
  coefficients.resize(
      EpsilonWindow{frame.min_power(), maximum}.width(), ComplexBall(0));
  return EpsilonFrame<ComplexBall>(
      {frame.min_power(), maximum}, std::move(coefficients));
}

inline FiniteLaurentVector<ComplexBall> zero_extend_acb_match_candidate(
    const FiniteLaurentVector<ComplexBall>& vector, std::size_t additional,
    const std::string& context) {
  FiniteLaurentVector<ComplexBall> result;
  result.reserve(vector.size());
  for (const auto& frame : vector)
    result.push_back(zero_extend_acb_match_candidate(
        frame, additional, context));
  return result;
}

inline FiniteLaurentMatrix<ComplexBall> zero_extend_acb_match_candidate(
    const FiniteLaurentMatrix<ComplexBall>& matrix, std::size_t additional,
    const std::string& context) {
  FiniteLaurentMatrix<ComplexBall> result(matrix.size());
  for (std::size_t row = 0; row < matrix.size(); ++row) {
    result[row].reserve(matrix[row].size());
    for (const auto& frame : matrix[row])
      result[row].push_back(zero_extend_acb_match_candidate(
          frame, additional, context));
  }
  return result;
}

inline FiniteLaurentVector<ComplexBall>
canonicalize_acb_match_candidate(
    const FiniteLaurentVector<ComplexBall>& vector, int chop_digits,
    const std::string& context) {
  FiniteLaurentVector<ComplexBall> result;
  result.reserve(vector.size());
  for (std::size_t row = 0; row < vector.size(); ++row) {
    auto coefficients = vector[row].coefficients();
    for (auto& coefficient : coefficients)
      coefficient =
          ScalarTraits<ComplexBall>::canonicalized(coefficient, chop_digits);
    result.push_back(canonicalize_certified_or_preserve_ambiguous(
        EpsilonFrame<ComplexBall>(
            vector[row].window(), std::move(coefficients)),
        context, row));
  }
  return result;
}

inline EpsilonFrame<ComplexBall> truncate_acb_match_candidate(
    const EpsilonFrame<ComplexBall>& frame, std::int32_t requested_top) {
  if (frame.complete_max() <= requested_top) return frame;
  if (frame.min_power() > requested_top)
    return EpsilonFrame<ComplexBall>::zero(requested_top);
  return frame.truncated(requested_top);
}

inline FiniteLaurentVector<ComplexBall> truncate_acb_match_candidate(
    FiniteLaurentVector<ComplexBall> vector,
    std::int32_t requested_top) {
  for (auto& frame : vector)
    frame = truncate_acb_match_candidate(frame, requested_top);
  return vector;
}

inline void require_complete_exact_saturation_record(
    const EpsilonLatticeSaturationResult<Rational>& record,
    std::size_t dimension, const std::string& context) {
  const auto transformation_columns = rectangular_columns(
      record.transformation, "exact saturation transformation");
  const auto transformed_columns = rectangular_columns(
      record.basis_times_transformation,
      "exact saturation transformed basis");
  const auto& diagnostics = record.diagnostics;
  if (record.transformation.size() != dimension ||
      transformation_columns != dimension ||
      record.basis_times_transformation.size() < dimension ||
      transformed_columns != dimension ||
      diagnostics.initial_column_valuations.size() != dimension ||
      diagnostics.initial_column_shifts.size() != dimension ||
      diagnostics.final_leading_rank != dimension ||
      diagnostics.normalized_determinant_valuation < 0 ||
      diagnostics.actions.size() != static_cast<std::size_t>(
          diagnostics.normalized_determinant_valuation))
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InvalidSaturationLattice,
        context + ": exact saturation record is incomplete or inconsistent");
  for (const auto& action : diagnostics.actions)
    if (action.target_column >= dimension ||
        action.null_relation.size() != dimension)
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::InvalidSaturationLattice,
          context + ": exact saturation action dimensions are inconsistent");
}

}  // namespace matching_detail

// Certified finite Acb matching after exact epsilon-lattice saturation.
//
// The structural phase is deliberately not repeated numerically: `record`
// must have been produced by exact Rational saturation and supplies T plus
// the completed rank/action proof.  Acb is used only to specialize the exact
// coefficients, apply F*T, and certify that every chosen solve pivot excludes
// zero.  One Laurent factorization is replayed for all residual corrections.
// The returned windows remain the honest EpsilonFrame windows; refinement
// never pads a coefficient lost to finite Laurent arithmetic.  This arithmetic
// primitive intentionally does not bind provenance between the exact and Acb
// bases; the retaining session/checkpoint layer must establish that identity.
inline RefinedAcbLaurentMatch refine_acb_finite_laurent_match(
    const FiniteLaurentMatrix<ComplexBall>& basis,
    const FiniteLaurentVector<ComplexBall>& right_hand_side,
    const ExactLaurentMatrix<ComplexBall>& transformation,
    const AcbLaurentRefinementOptions& options,
    const std::string& context = "refined Acb finite Laurent match",
    bool exact_formal_negative_coefficients_are_zero = false,
    bool exact_right_transformation_residual = false) {
  const auto dimension = matching_detail::rectangular_columns(
      basis, "Acb matching basis");
  if (basis.size() != dimension || right_hand_side.size() != dimension)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::DimensionMismatch,
        context + ": basis must be square and match the rhs");
  if (transformation.size() != dimension ||
      matching_detail::rectangular_columns(
          transformation, "Acb saturation transformation") != dimension)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InvalidSaturationLattice,
        context + ": Acb saturation transformation is not square");
  const auto transformed_basis =
      right_multiply_finite_by_exact_laurent(basis, transformation);

  std::int32_t transformation_minimum = 0;
  bool identity_transformation = true;
  for (const auto& row : transformation)
    for (const auto& polynomial : row)
      if (const auto minimum = polynomial.minimum_power();
          minimum.has_value())
        transformation_minimum =
            std::min(transformation_minimum, *minimum);
  for (std::size_t row = 0; row < dimension; ++row)
    for (std::size_t column = 0; column < dimension; ++column) {
      const auto& polynomial = transformation[row][column];
      if (row != column) {
        identity_transformation =
            identity_transformation && polynomial.is_zero();
        continue;
      }
      identity_transformation = identity_transformation &&
          polynomial.terms().size() == 1 &&
          polynomial.terms().begin()->first == 0 &&
          (polynomial.terms().begin()->second - ComplexBall(1)).is_zero();
    }
  const auto transformation_depth = static_cast<std::size_t>(
      -static_cast<std::int64_t>(transformation_minimum));
  bool negative_input_support = false;
  for (const auto& row : basis)
    for (const auto& frame : row)
      negative_input_support =
          negative_input_support || frame.min_power() < 0;
  // A nonnegative identity problem keeps the ordinary fast path byte-for-byte
  // finite.  An identity problem with a represented Laurent tail needs one
  // proposal coefficient to avoid losing that tail in the candidate solve.
  // A nontrivial T gets enough private proposal width to carry its polar
  // shift (if any) and a small guard for the quotient solve.  If this is ever
  // insufficient, the untouched physical residual remains incomplete and
  // requests real input width; no padded coefficient can authorize success.
  const auto candidate_padding = identity_transformation
      ? (negative_input_support ? std::size_t{1} : std::size_t{0})
      : transformation_depth + std::size_t{3};
  auto candidate_basis =
      matching_detail::zero_extend_acb_match_candidate(
          transformed_basis, candidate_padding,
          context + ": speculative transformed basis");
  const auto candidate_right_hand_side =
      matching_detail::zero_extend_acb_match_candidate(
          right_hand_side, candidate_padding,
          context + ": speculative transformed rhs");

  std::optional<ExactNonnegativePowerSeriesFactorization<ComplexBall>>
      power_series_factorization;
  std::optional<FiniteLaurentFactorization<ComplexBall>>
      laurent_factorization;
  if (exact_formal_negative_coefficients_are_zero)
    power_series_factorization =
        factor_exact_nonnegative_power_series_system(
            std::move(candidate_basis),
            context +
                ": certified transformed power-series factorization");
  else
    laurent_factorization =
        factor_preconditioned_acb_finite_laurent_system(
            candidate_basis,
            context + ": certified transformed factorization");
  const auto solve = [&](const FiniteLaurentVector<ComplexBall>& rhs,
                         const std::string& solve_context) {
    return power_series_factorization.has_value()
        ? solve_factorized_exact_nonnegative_power_series_system(
              *power_series_factorization, rhs, solve_context)
        : solve_factorized_finite_laurent_system(
              *laurent_factorization, rhs, solve_context);
  };
  auto transformed_weights = solve(
      candidate_right_hand_side, context + ": initial candidate solve");

  RefinedAcbLaurentMatch result;
  result.factorization_preconditioner =
      power_series_factorization.has_value()
      ? power_series_factorization->leading.preconditioner_kind
      : laurent_factorization->preconditioner_kind;
  std::string factorization_diagnostics =
      "; factorization=" + result.factorization_preconditioner;
  if (laurent_factorization.has_value()) {
    factorization_diagnostics += "; diagonal_pivots=[";
    for (std::size_t row = 0;
         row < laurent_factorization->upper.size(); ++row) {
      if (row != 0) factorization_diagnostics += ",";
      const auto& pivot = laurent_factorization->upper[row][row];
      const auto leading =
          matching_detail::certified_laurent_leading_power(pivot);
      factorization_diagnostics += "{min=" +
          std::to_string(pivot.min_power()) + ",max=" +
          std::to_string(pivot.complete_max()) + ",leading=" +
          (leading.power.has_value()
               ? std::to_string(*leading.power)
               : leading.first_ambiguous_power.has_value()
               ? "ambiguous@" +
                     std::to_string(*leading.first_ambiguous_power)
               : std::string("zero")) +
          "}";
    }
    factorization_diagnostics += "]";
  }
  std::optional<FiniteLaurentVector<ComplexBall>> last_complete_transformed;
  std::optional<FiniteLaurentVector<ComplexBall>> last_complete_weights;
  std::optional<FiniteLaurentVector<ComplexBall>> last_complete_residual;
  std::optional<std::int64_t> last_complete_pass_prefix;
  std::size_t last_complete_refinement_steps = 0;
  std::size_t last_complete_history_size = 0;
  for (;;) {
    auto weights = matching_detail::canonicalize_acb_match_candidate(
        apply_exact_laurent_matrix(transformation, transformed_weights),
        100, context + ": physical candidate weights");
    const auto residual_context =
        context + ": residual#" +
        std::to_string(result.residual_history.size()) +
        factorization_diagnostics;
    auto residual = [&] {
      if (!exact_right_transformation_residual)
        return matching_detail::evaluate_acb_matching_residual(
            basis, weights, right_hand_side, options, residual_context);

      // T is supplied by the retained exact Rational lattice proof.  The
      // two products F*(T*y) and (F*T)*y are therefore the same physical
      // vector, but evaluating the former with independent balls can destroy
      // the exact column dependencies which made F*T epsilon-regular.  Form
      // the residual through the latter route and certify that enclosure
      // against the original physical basis, weights, rhs, and scale.  Honest
      // finite-Laurent windows from both products are still intersected by
      // certify_precomputed_acb_matching_residual, so the exact identity
      // cannot manufacture an unknown upper coefficient.
      auto transformed =
          matching_detail::evaluate_acb_matching_residual(
              transformed_basis, transformed_weights, right_hand_side,
              options, residual_context + ":exact-right-frame");
      return matching_detail::certify_precomputed_acb_matching_residual(
          basis, weights, right_hand_side, std::move(transformed.residual),
          options, residual_context + ":physical-via-exact-right-frame");
    }();
    result.residual_history.push_back(residual.diagnostics);
    result.weights = std::move(weights);
    result.residual = std::move(residual.residual);

    const auto& latest = result.residual_history.back();
    const auto contiguous_pass_prefix = [&latest]() {
      const auto upper = std::min(
          latest.complete_window.complete_max,
          latest.required_complete_max);
      std::int64_t prefix =
          static_cast<std::int64_t>(latest.complete_window.min_power) - 1;
      for (std::int64_t power = latest.complete_window.min_power;
           power <= upper; ++power) {
        bool saw_coefficient = false;
        bool all_rows_pass = true;
        for (const auto& coefficient : latest.coefficients) {
          if (coefficient.epsilon_power != power) continue;
          saw_coefficient = true;
          all_rows_pass = all_rows_pass &&
              coefficient.verdict == AcbMatchingResidualVerdict::Pass;
        }
        if (!saw_coefficient || !all_rows_pass) break;
        prefix = power;
      }
      return prefix;
    }();

    // A finite-Laurent correction can consume upper coefficients or worsen
    // the contiguous certified residual prefix even when an earlier iterate
    // already covered the caller's requirement.  Keep corrections
    // speculative and retain the best complete prefix transactionally.  A
    // later correction may recover and win; if refinement stops first, the
    // worse proposal must not replace the best usable iterate.  Otherwise a
    // wider private reservoir merely permits one more destructive correction
    // and can publish the same reservoir retry forever.
    if (latest.complete_through_required &&
        (!last_complete_pass_prefix.has_value() ||
         contiguous_pass_prefix > *last_complete_pass_prefix)) {
      last_complete_transformed = transformed_weights;
      last_complete_weights = result.weights;
      last_complete_residual = result.residual;
      last_complete_pass_prefix = contiguous_pass_prefix;
      last_complete_refinement_steps = result.refinement_steps;
      last_complete_history_size = result.residual_history.size();
    }
    if (latest.verdict == AcbMatchingResidualVerdict::Pass ||
        !latest.complete_through_required ||
        result.refinement_steps >= options.max_refinement_steps) {
      if (last_complete_transformed.has_value() &&
          last_complete_history_size != result.residual_history.size()) {
        transformed_weights = std::move(*last_complete_transformed);
        result.weights = std::move(*last_complete_weights);
        result.residual = std::move(*last_complete_residual);
        result.refinement_steps = last_complete_refinement_steps;
        result.residual_history.resize(last_complete_history_size);
      }
      break;
    }

    // Only correct the declared finite Laurent quotient.  Roundoff remnants
    // below its lower floor are retained in `result.residual` for diagnosis,
    // but feeding them back would manufacture an ever-deeper negative tail
    // that the caller neither requested nor can use.
    auto correction_rhs = matching_detail::project_acb_residual_to_laurent_floor(
        result.residual, latest.complete_window.min_power,
        context + ": residual correction projection");
    correction_rhs = matching_detail::zero_extend_acb_match_candidate(
        correction_rhs, candidate_padding,
        context + ": speculative correction rhs");
    auto correction = solve(
        correction_rhs,
        context + ": correction#" +
            std::to_string(result.refinement_steps + 1));
    for (std::size_t column = 0; column < dimension; ++column)
      transformed_weights[column] =
          transformed_weights[column] + correction[column];
    ++result.refinement_steps;
  }
  auto published_weight_top = right_hand_side.front().complete_max();
  for (const auto& row : right_hand_side)
    published_weight_top =
        std::min(published_weight_top, row.complete_max());
  result.weights = matching_detail::truncate_acb_match_candidate(
      std::move(result.weights), published_weight_top);
  const auto transformed_published_top = matching_detail::checked_power(
      static_cast<std::int64_t>(published_weight_top) -
          transformation_minimum,
      "published transformed matching-weight maximum");
  result.transformed_weights =
      matching_detail::truncate_acb_match_candidate(
          std::move(transformed_weights), transformed_published_top);
  return result;
}

inline RefinedAcbLaurentMatch refine_acb_finite_laurent_match(
    const FiniteLaurentMatrix<ComplexBall>& basis,
    const FiniteLaurentVector<ComplexBall>& right_hand_side,
    const EpsilonLatticeSaturationResult<Rational>& exact_saturation_record,
    const AcbLaurentRefinementOptions& options,
    const std::string& context = "refined Acb finite Laurent match",
    bool exact_formal_negative_coefficients_are_zero = false) {
  const auto dimension = matching_detail::rectangular_columns(
      basis, "Acb matching basis");
  matching_detail::require_complete_exact_saturation_record(
      exact_saturation_record, dimension, context);
  return refine_acb_finite_laurent_match(
      basis, right_hand_side,
      specialize_exact_rational_laurent_matrix_to_acb(
          exact_saturation_record.transformation),
      options, context, exact_formal_negative_coefficients_are_zero, true);
}

}  // namespace diffexp2
