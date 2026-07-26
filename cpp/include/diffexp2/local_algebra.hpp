#pragma once

#include "diffexp2/local_solution.hpp"
#include "diffexp2/matching.hpp"

#include <flint/acb_poly.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace diffexp2 {

// Full analytic completion of one epsilon coefficient of a prepared rational
// multiplier after its common center pole has been removed.  The finite
// Taylor kernels below are sufficient for stored algebra, but never determine
// an unseen tail.  A caller which wants a rigorous projected-line certificate
// may additionally retain these numerator/denominator polynomials.  The tail
// layer independently replays polynomial division through the complete stored
// Taylor order before trusting this completion.
template <typename Scalar>
struct PreparedRationalAnalyticCoefficient {
  std::vector<Scalar> numerator;    // ascending powers of the chart variable
  std::vector<Scalar> denominator;  // ascending powers; constant is nonzero
};

// Prepared exact rational multiplication for a finite local-sector slab.
// Wolfram still owns exact pole classification in the first SCC milestone:
// it expands c(t,eps) as
//
//   eps^epsilon_shift * t^-center_pole_order
//     * sum_j eps^j sum_n kernels[j][n] t^n.
//
// The native operation below owns only the finite triangular convolutions.
// Every kernel row must contain at least the input Taylor width.  `kernels`
// may be a consumer prefix shorter than the input epsilon width; in that
// case the complete product stops at the shorter of the two finite prefixes.
// This matches the finite convolution used by SectorSeries`MultiplyRational
// without requiring private source coefficients that the caller's capped
// output can never consume.
template <typename Scalar>
struct PreparedRationalTaylorMultiplier {
  std::int32_t epsilon_shift = 0;
  std::uint32_t center_pole_order = 0;
  std::vector<std::vector<Scalar>> kernels;  // [epsilon][Taylor]
  std::optional<std::vector<PreparedRationalAnalyticCoefficient<Scalar>>>
      analytic_coefficients;
  std::string exact_identity;
  // Exact preparation provenance, never a conclusion drawn from a numeric
  // slab.  An Acb tensor whose entries currently enclose exact zero remains
  // active unless the rational multiplier was proved structurally zero
  // before specialization.
  bool proven_zero = false;

  [[nodiscard]] bool structurally_zero() const {
    return proven_zero;
  }
};

template <typename Scalar>
struct PreparedSparseLocalMultiplierMatrix {
  struct Entry {
    std::uint32_t row = 0;
    std::uint32_t column = 0;
    PreparedRationalTaylorMultiplier<Scalar> multiplier;
  };

  std::uint32_t rows = 0;
  std::uint32_t columns = 0;
  std::vector<Entry> entries;
  std::string exact_identity;
};

inline PreparedRationalTaylorMultiplier<ComplexBall>
specialize_prepared_rational_multiplier_to_acb(
    const PreparedRationalTaylorMultiplier<Rational>& exact) {
  PreparedRationalTaylorMultiplier<ComplexBall> numeric;
  numeric.epsilon_shift = exact.epsilon_shift;
  numeric.center_pole_order = exact.center_pole_order;
  numeric.exact_identity = exact.exact_identity;
  numeric.proven_zero = exact.proven_zero;
  numeric.kernels.reserve(exact.kernels.size());
  for (const auto& exact_kernel : exact.kernels) {
    std::vector<ComplexBall> kernel;
    kernel.reserve(exact_kernel.size());
    for (const auto& coefficient : exact_kernel)
      kernel.push_back(ComplexBall::from_strings(coefficient.str()));
    numeric.kernels.push_back(std::move(kernel));
  }
  if (exact.analytic_coefficients.has_value()) {
    numeric.analytic_coefficients.emplace();
    numeric.analytic_coefficients->reserve(
        exact.analytic_coefficients->size());
    for (const auto& exact_rational :
         *exact.analytic_coefficients) {
      PreparedRationalAnalyticCoefficient<ComplexBall> rational;
      rational.numerator.reserve(exact_rational.numerator.size());
      for (const auto& coefficient : exact_rational.numerator)
        rational.numerator.push_back(
            ComplexBall::from_strings(coefficient.str()));
      rational.denominator.reserve(
          exact_rational.denominator.size());
      for (const auto& coefficient : exact_rational.denominator)
        rational.denominator.push_back(
            ComplexBall::from_strings(coefficient.str()));
      numeric.analytic_coefficients->push_back(std::move(rational));
    }
  }
  return numeric;
}

inline PreparedSparseLocalMultiplierMatrix<ComplexBall>
specialize_prepared_rational_matrix_to_acb(
    const PreparedSparseLocalMultiplierMatrix<Rational>& exact) {
  PreparedSparseLocalMultiplierMatrix<ComplexBall> numeric;
  numeric.rows = exact.rows;
  numeric.columns = exact.columns;
  numeric.exact_identity = exact.exact_identity;
  numeric.entries.reserve(exact.entries.size());
  for (const auto& exact_entry : exact.entries)
    numeric.entries.push_back(
        typename PreparedSparseLocalMultiplierMatrix<
            ComplexBall>::Entry{
            exact_entry.row, exact_entry.column,
            specialize_prepared_rational_multiplier_to_acb(
                exact_entry.multiplier)});
  return numeric;
}

namespace local_algebra_detail {

inline std::optional<ComplexBall> signed_real_evaluation_ball(
    const RealEvaluationPoint& point) {
  if (point.exact_coordinate.empty() || point.sign < -1 || point.sign > 1 ||
      !arb_is_zero(acb_imagref(point.modulus.raw())) ||
      !arb_is_nonnegative(acb_realref(point.modulus.raw())))
    return std::nullopt;
  if (point.sign == 0)
    return point.modulus.is_zero()
        ? std::optional<ComplexBall>(point.modulus)
        : std::nullopt;
  if (point.modulus.contains_zero()) return std::nullopt;
  return point.sign < 0
      ? std::optional<ComplexBall>(-point.modulus)
      : std::optional<ComplexBall>(point.modulus);
}

inline std::optional<ComplexBall> evaluate_acb_polynomial_at(
    const std::vector<ComplexBall>& coefficients,
    const ComplexBall& point) {
  if (coefficients.empty()) return std::nullopt;
  auto value = coefficients.back();
  for (std::size_t index = coefficients.size() - 1; index > 0; --index)
    value = value * point + coefficients[index - 1];
  return value;
}

// Specialize the complete retained rational function at an exact real chart
// coordinate.  This deliberately does not replay `kernels`: those are only a
// finite Taylor prefix around the chart center and cannot certify the gauge
// at an off-center match point.
inline std::optional<EpsilonFrame<ComplexBall>>
specialize_prepared_rational_multiplier_at_real_point(
    const PreparedRationalTaylorMultiplier<ComplexBall>& multiplier,
    const RealEvaluationPoint& point) {
  if (multiplier.proven_zero ||
      !multiplier.analytic_coefficients.has_value() ||
      multiplier.analytic_coefficients->empty() ||
      multiplier.analytic_coefficients->size() != multiplier.kernels.size())
    return std::nullopt;
  const auto signed_point = signed_real_evaluation_ball(point);
  if (!signed_point.has_value() ||
      (multiplier.center_pole_order != 0 &&
       signed_point->contains_zero()))
    return std::nullopt;

  ComplexBall center_pole_factor(1);
  if (multiplier.center_pole_order != 0) {
    acb_pow_ui(center_pole_factor.raw(), signed_point->raw(),
               multiplier.center_pole_order, ComplexBall::precision());
    if (center_pole_factor.contains_zero()) return std::nullopt;
  }

  std::vector<ComplexBall> specialized;
  specialized.reserve(multiplier.analytic_coefficients->size());
  for (const auto& rational : *multiplier.analytic_coefficients) {
    auto numerator =
        evaluate_acb_polynomial_at(rational.numerator, *signed_point);
    auto denominator =
        evaluate_acb_polynomial_at(rational.denominator, *signed_point);
    if (!numerator.has_value() || !denominator.has_value() ||
        denominator->contains_zero())
      return std::nullopt;
    auto value = *numerator / *denominator;
    if (multiplier.center_pole_order != 0)
      value = value / center_pole_factor;
    specialized.push_back(std::move(value));
  }
  return EpsilonFrame<ComplexBall>(
      multiplier.epsilon_shift, std::move(specialized));
}

inline std::optional<std::string>
diagnose_prepared_rational_specialization_at_real_point(
    const PreparedRationalTaylorMultiplier<ComplexBall>& multiplier,
    const RealEvaluationPoint& point) {
  if (multiplier.proven_zero)
    return "structural-zero multiplier was retained as an active entry";
  if (!multiplier.analytic_coefficients.has_value())
    return "complete analytic rational coefficients are unavailable";
  if (multiplier.analytic_coefficients->empty())
    return "complete analytic rational coefficient list is empty";
  if (multiplier.analytic_coefficients->size() != multiplier.kernels.size())
    return "analytic rational and epsilon-kernel counts differ";
  const auto signed_point = signed_real_evaluation_ball(point);
  if (!signed_point.has_value())
    return "receiving point has no certified signed real specialization";
  if (multiplier.center_pole_order != 0 &&
      signed_point->contains_zero())
    return "receiving point intersects the gauge center pole";

  if (multiplier.center_pole_order != 0) {
    ComplexBall center_pole_factor(1);
    acb_pow_ui(center_pole_factor.raw(), signed_point->raw(),
               multiplier.center_pole_order, ComplexBall::precision());
    if (center_pole_factor.contains_zero())
      return "evaluated gauge center-pole factor contains zero";
  }
  for (std::size_t epsilon = 0;
       epsilon < multiplier.analytic_coefficients->size(); ++epsilon) {
    const auto& rational =
        (*multiplier.analytic_coefficients)[epsilon];
    auto numerator =
        evaluate_acb_polynomial_at(rational.numerator, *signed_point);
    if (!numerator.has_value())
      return "analytic numerator is empty at epsilon coefficient " +
          std::to_string(epsilon);
    auto denominator =
        evaluate_acb_polynomial_at(rational.denominator, *signed_point);
    if (!denominator.has_value())
      return "analytic denominator is empty at epsilon coefficient " +
          std::to_string(epsilon);
    if (denominator->contains_zero())
      return std::string(denominator->is_zero()
              ? "analytic denominator is exactly zero"
              : "analytic denominator enclosure contains zero") +
          " at epsilon coefficient " + std::to_string(epsilon);
  }
  return std::nullopt;
}

inline std::optional<FiniteLaurentVector<ComplexBall>>
apply_prepared_sparse_epsilon_matrix_at_real_point(
    const PreparedSparseLocalMultiplierMatrix<ComplexBall>& matrix,
    const FiniteLaurentVector<ComplexBall>& input,
    const RealEvaluationPoint& point) {
  if (matrix.rows != input.size() || matrix.columns != input.size())
    throw std::invalid_argument(
        "point-specialized epsilon matrix dimensions differ from its vector");
  std::vector<std::optional<EpsilonFrame<ComplexBall>>> rows(matrix.rows);
  for (const auto& entry : matrix.entries) {
    if (entry.row >= matrix.rows || entry.column >= matrix.columns)
      throw std::invalid_argument(
          "point-specialized epsilon matrix entry is outside its shape");
    if (entry.multiplier.proven_zero) continue;
    auto multiplier =
        specialize_prepared_rational_multiplier_at_real_point(
            entry.multiplier, point);
    if (!multiplier.has_value()) return std::nullopt;
    auto term = *multiplier * input[entry.column];
    rows[entry.row] = rows[entry.row].has_value()
        ? *rows[entry.row] + term
        : std::move(term);
  }

  FiniteLaurentVector<ComplexBall> output;
  output.reserve(rows.size());
  for (auto& row : rows) {
    if (!row.has_value())
      throw std::invalid_argument(
          "point-specialized invertible epsilon matrix has an empty row");
    output.push_back(std::move(*row));
  }
  return output;
}

// Move-only owner for a FLINT polynomial.  Acb local algebra uses these to
// replace the scalar O(N^2) Taylor convolution with Arb's polynomial kernel
// while preserving the same finite epsilon/Taylor rectangle.
class AcbPolynomial {
 public:
  AcbPolynomial() { acb_poly_init(value_); }
  AcbPolynomial(const AcbPolynomial&) = delete;
  AcbPolynomial& operator=(const AcbPolynomial&) = delete;
  AcbPolynomial(AcbPolynomial&& other) noexcept : AcbPolynomial() {
    acb_poly_swap(value_, other.value_);
  }
  AcbPolynomial& operator=(AcbPolynomial&& other) noexcept {
    acb_poly_swap(value_, other.value_);
    return *this;
  }
  ~AcbPolynomial() { acb_poly_clear(value_); }

  acb_poly_struct* raw() { return value_; }
  const acb_poly_struct* raw() const { return value_; }

 private:
  acb_poly_t value_;
};

// Materialize the finite Taylor kernel of a compact rational coefficient.
// Transport rows arrive with the complete low-degree numerator/denominator
// polynomials, so expanding hundreds of coefficients in Wolfram and shipping
// them as JSON is redundant.  Acb uses FLINT's native power-series division;
// exact coefficient fields retain the same elementary division recurrence.
template <typename Scalar>
std::vector<Scalar> expand_rational_taylor(
    const PreparedRationalAnalyticCoefficient<Scalar>& rational,
    std::size_t taylor_width) {
  if (rational.numerator.empty() || rational.denominator.empty())
    throw std::invalid_argument(
        "compact rational Taylor source has an empty polynomial");
  if (taylor_width == 0)
    throw std::invalid_argument(
        "compact rational Taylor source requested an empty kernel");

  if constexpr (std::is_same_v<Scalar, ComplexBall>) {
    if (rational.denominator.front().contains_zero())
      throw std::domain_error(
          "compact rational Taylor denominator contains zero at the chart center");
    if (taylor_width >
        static_cast<std::size_t>(std::numeric_limits<slong>::max()))
      throw std::overflow_error(
          "compact rational Taylor width exceeds FLINT range");
    AcbPolynomial numerator;
    AcbPolynomial denominator;
    AcbPolynomial quotient;
    for (std::size_t i = 0; i < rational.numerator.size(); ++i)
      if (!rational.numerator[i].is_zero())
        acb_poly_set_coeff_acb(numerator.raw(), static_cast<slong>(i),
                               rational.numerator[i].raw());
    for (std::size_t i = 0; i < rational.denominator.size(); ++i)
      if (!rational.denominator[i].is_zero())
        acb_poly_set_coeff_acb(denominator.raw(), static_cast<slong>(i),
                               rational.denominator[i].raw());
    acb_poly_div_series(quotient.raw(), numerator.raw(), denominator.raw(),
                        static_cast<slong>(taylor_width),
                        ComplexBall::precision());
    std::vector<Scalar> kernel(taylor_width);
    for (std::size_t i = 0; i < taylor_width; ++i) {
      acb_poly_get_coeff_acb(kernel[i].raw(), quotient.raw(),
                             static_cast<slong>(i));
      if (!kernel[i].is_finite())
        throw std::domain_error(
            "compact rational Taylor expansion produced a nonfinite Acb coefficient");
    }
    return kernel;
  } else {
    if (ScalarTraits<Scalar>::is_zero(rational.denominator.front()))
      throw std::domain_error(
          "compact rational Taylor denominator is zero at the chart center");
    std::vector<Scalar> kernel(taylor_width, ScalarTraits<Scalar>::zero());
    for (std::size_t m = 0; m < taylor_width; ++m) {
      auto value = m < rational.numerator.size()
          ? rational.numerator[m]
          : ScalarTraits<Scalar>::zero();
      const auto denominator_degree = rational.denominator.size() - 1;
      for (std::size_t l = 1; l <= std::min(m, denominator_degree); ++l)
        value -= rational.denominator[l] * kernel[m - l];
      kernel[m] = value / rational.denominator.front();
    }
    return kernel;
  }
}

inline std::int32_t checked_i32(std::int64_t value, const char* label) {
  if (value < std::numeric_limits<std::int32_t>::min() ||
      value > std::numeric_limits<std::int32_t>::max())
    throw std::overflow_error(std::string(label) + " exceeds int32 range");
  return static_cast<std::int32_t>(value);
}

inline bool same_descriptor(const ExactScalarDescriptor& left,
                            const ExactScalarDescriptor& right) {
  if (left.domain != right.domain || left.canonical != right.canonical ||
      left.symbols != right.symbols || left.is_zero != right.is_zero ||
      left.is_integer != right.is_integer || left.sign != right.sign ||
      left.specialization.has_value() != right.specialization.has_value())
    return false;
  return !left.specialization.has_value() ||
         acb_equal(left.specialization->raw(), right.specialization->raw());
}

inline bool same_prescriptions(const std::vector<Prescription>& left,
                               const std::vector<Prescription>& right) {
  if (left.size() != right.size()) return false;
  for (std::size_t i = 0; i < left.size(); ++i) {
    if (left[i].factor_exact != right[i].factor_exact ||
        left[i].sign != right[i].sign ||
        left[i].multiplicity != right[i].multiplicity ||
        left[i].leading_coefficient_sign !=
            right[i].leading_coefficient_sign)
      return false;
  }
  return true;
}

inline bool same_chart(const ChartGeometry& left,
                       const ChartGeometry& right) {
  if (left.center_exact != right.center_exact ||
      left.scale_exact != right.scale_exact ||
      left.infinite_radius != right.infinite_radius)
    return false;
  if (left.infinite_radius) return true;
  if (!left.radius_exact.empty() && !right.radius_exact.empty())
    return left.radius_exact == right.radius_exact;
  return acb_equal(left.radius.raw(), right.radius.raw());
}

template <typename Scalar>
void require_factorized_local_basis_contract(
    const std::vector<LocalSolution<Scalar>>& basis,
    std::uint32_t dimension, EpsilonWindow requested_window,
    std::int32_t required_complete_max,
    const std::string& context) {
  if (basis.empty() || basis.size() != dimension)
    throw std::logic_error(
        context +
        ": factorized local basis has the wrong column count");
  const auto& reference = basis.front();
  validate_local_solution(reference, false);
  if (reference.dimension != dimension ||
      reference.sectors.empty() ||
      reference.epsilon.min_power >
          requested_window.complete_max ||
      reference.epsilon.complete_max <
          required_complete_max)
    throw std::logic_error(
        context +
        ": factorized reference column has invalid shape or epsilon "
        "coverage");
  for (const auto& column : basis) {
    validate_local_solution(column, false);
    if (column.dimension != dimension ||
        column.sectors.empty() ||
        !same_chart(column.chart, reference.chart) ||
        !same_prescriptions(
            column.prescriptions, reference.prescriptions) ||
        column.taylor_complete_max !=
            reference.taylor_complete_max ||
        column.epsilon.min_power >
            requested_window.complete_max ||
        column.epsilon.complete_max <
            required_complete_max)
      throw std::logic_error(
          context +
          ": factorized column has incompatible chart, prescriptions, "
          "shape, Taylor width, or epsilon coverage");
  }
}

template <typename Scalar>
void require_factorized_receiving_local_compatibility(
    const LocalSolution<Scalar>& factorized,
    const LocalSolution<Scalar>& receiving,
    const std::string& context) {
  validate_local_solution(factorized, false);
  validate_local_solution(receiving, false);
  if (!same_chart(factorized.chart, receiving.chart) ||
      !same_prescriptions(
          factorized.prescriptions, receiving.prescriptions) ||
      factorized.dimension != receiving.dimension ||
      factorized.taylor_complete_max !=
          receiving.taylor_complete_max)
    throw std::logic_error(
        context +
        ": factorized local is incompatible with its receiving chart, "
        "prescriptions, dimension, or Taylor width");
}

template <typename Scalar>
void require_same_local_space(const LocalSolution<Scalar>& left,
                              const LocalSolution<Scalar>& right) {
  if (!same_chart(left.chart, right.chart) ||
      !same_prescriptions(left.prescriptions, right.prescriptions) ||
      left.dimension != right.dimension)
    throw std::invalid_argument(
        "local algebra requires identical chart, prescription and component spaces");
}

inline std::size_t flat_index(std::size_t epsilon_index,
                              std::size_t taylor_index,
                              std::uint32_t component,
                              std::size_t taylor_width,
                              std::uint32_t dimension) {
  return ((epsilon_index * taylor_width + taylor_index) * dimension) +
         component;
}

inline ExactScalarDescriptor subtract_nonnegative_integer(
    const ExactScalarDescriptor& input, std::uint32_t amount) {
  if (amount == 0) return input;
  const auto amount_string = std::to_string(amount);
  if (input.domain == ExactDomain::Rational) {
    const auto value = Rational(input.canonical) - Rational(amount_string);
    return ExactScalarDescriptor::rational(value.str());
  }

  auto output = input;
  output.canonical = "(" + input.canonical + ")-" + amount_string;
  // Subtracting an integer preserves integer/noninteger status, but zero and
  // sign need an exact comparison which this descriptor layer must not guess.
  output.is_zero = TruthValue::Unknown;
  output.sign = ExactSign::Unknown;
  if (input.specialization.has_value())
    output.specialization = *input.specialization -
        ComplexBall::from_strings(amount_string);
  return output;
}

template <typename Scalar>
bool same_sector_tag(const LocalSector<Scalar>& left,
                     const LocalSector<Scalar>& right) {
  return left.log_power == right.log_power &&
         same_descriptor(left.a, right.a) &&
         same_descriptor(left.b, right.b);
}

template <typename Scalar>
bool sector_less(const LocalSector<Scalar>& left,
                 const LocalSector<Scalar>& right) {
  return std::tie(left.a.domain, left.a.canonical, left.b.domain,
                  left.b.canonical, left.log_power) <
         std::tie(right.a.domain, right.a.canonical, right.b.domain,
                  right.b.canonical, right.log_power);
}

template <typename Scalar>
LocalSolution<Scalar> with_selected_component(
    const LocalSolution<Scalar>& input, std::uint32_t component) {
  if (component >= input.dimension)
    throw std::out_of_range("selected local component is outside dimension");
  LocalSolution<Scalar> output;
  output.chart = input.chart;
  output.epsilon = input.epsilon;
  output.taylor_complete_max = input.taylor_complete_max;
  output.dimension = 1;
  output.prescriptions = input.prescriptions;
  output.checkpoint_identity = input.checkpoint_identity + ":component:" +
      std::to_string(component);
  output.sectors.reserve(input.sectors.size());
  for (const auto& sector : input.sectors) {
    LocalSector<Scalar> selected;
    selected.a = sector.a;
    selected.b = sector.b;
    selected.log_power = sector.log_power;
    selected.coefficients.assign(output.sector_size(),
                                 ScalarTraits<Scalar>::zero());
    for (std::size_t ei = 0; ei < input.epsilon.width(); ++ei)
      for (std::size_t n = 0; n < input.taylor_width(); ++n)
        selected.coefficients[flat_index(ei, n, 0, output.taylor_width(), 1)] =
            sector.coefficients[flat_index(
                ei, n, component, input.taylor_width(), input.dimension)];
    output.sectors.push_back(std::move(selected));
  }
  return output;
}

template <typename Scalar>
LocalSolution<Scalar> with_selected_components(
    const LocalSolution<Scalar>& input,
    const std::vector<std::uint32_t>& components) {
  if (components.empty())
    throw std::invalid_argument("selected local component list is empty");
  std::vector<std::uint8_t> seen(input.dimension, 0);
  for (const auto component : components) {
    if (component >= input.dimension || seen[component])
      throw std::out_of_range(
          "selected local components are outside dimension or duplicated");
    seen[component] = 1;
  }
  LocalSolution<Scalar> output;
  output.chart = input.chart;
  output.epsilon = input.epsilon;
  output.taylor_complete_max = input.taylor_complete_max;
  output.dimension = static_cast<std::uint32_t>(components.size());
  output.prescriptions = input.prescriptions;
  output.checkpoint_identity =
      input.checkpoint_identity + ":components";
  for (const auto component : components)
    output.checkpoint_identity += ":" + std::to_string(component);
  output.sectors.reserve(input.sectors.size());
  for (const auto& sector : input.sectors) {
    LocalSector<Scalar> selected;
    selected.a = sector.a;
    selected.b = sector.b;
    selected.log_power = sector.log_power;
    selected.coefficients.assign(output.sector_size(),
                                 ScalarTraits<Scalar>::zero());
    for (std::size_t ei = 0; ei < input.epsilon.width(); ++ei)
      for (std::size_t n = 0; n < input.taylor_width(); ++n)
        for (std::uint32_t local = 0; local < output.dimension; ++local)
          selected.coefficients[flat_index(
              ei, n, local, output.taylor_width(), output.dimension)] =
              sector.coefficients[flat_index(
                  ei, n, components[local], input.taylor_width(),
                  input.dimension)];
    output.sectors.push_back(std::move(selected));
  }
  return output;
}

template <typename Scalar>
LocalSolution<Scalar> embedded_components(
    const LocalSolution<Scalar>& input,
    const std::vector<std::uint32_t>& components,
    std::uint32_t dimension) {
  if (input.dimension == 0 || input.dimension != components.size() ||
      dimension == 0)
    throw std::invalid_argument("invalid local component embedding");
  std::vector<std::uint8_t> seen(dimension, 0);
  for (const auto component : components) {
    if (component >= dimension || seen[component])
      throw std::invalid_argument("invalid local component embedding");
    seen[component] = 1;
  }
  LocalSolution<Scalar> output;
  output.chart = input.chart;
  output.epsilon = input.epsilon;
  output.taylor_complete_max = input.taylor_complete_max;
  output.dimension = dimension;
  output.prescriptions = input.prescriptions;
  output.checkpoint_identity = input.checkpoint_identity + ":embedded-block:" +
      std::to_string(input.dimension) + ":" + std::to_string(dimension);
  for (const auto component : components)
    output.checkpoint_identity += ":" + std::to_string(component);
  output.sectors.reserve(input.sectors.size());
  for (const auto& sector : input.sectors) {
    LocalSector<Scalar> embedded;
    embedded.a = sector.a;
    embedded.b = sector.b;
    embedded.log_power = sector.log_power;
    embedded.coefficients.assign(output.sector_size(),
                                 ScalarTraits<Scalar>::zero());
    for (std::size_t ei = 0; ei < input.epsilon.width(); ++ei)
      for (std::size_t n = 0; n < input.taylor_width(); ++n)
        for (std::uint32_t local = 0; local < input.dimension; ++local)
          embedded.coefficients[flat_index(
              ei, n, components[local], output.taylor_width(), dimension)] =
              sector.coefficients[flat_index(
                  ei, n, local, input.taylor_width(), input.dimension)];
    output.sectors.push_back(std::move(embedded));
  }
  return output;
}

template <typename Scalar>
LocalSolution<Scalar> embedded_component(const LocalSolution<Scalar>& input,
                                         std::uint32_t component,
                                         std::uint32_t dimension) {
  auto output = embedded_components(
      input, std::vector<std::uint32_t>{component}, dimension);
  output.checkpoint_identity = input.checkpoint_identity + ":embedded:" +
      std::to_string(component) + ":" + std::to_string(dimension);
  return output;
}

}  // namespace local_algebra_detail

// Merge only exactly identical tags.  Integer-spaced tower compaction is a
// representation optimization, not a semantic requirement, and is deferred
// until the native exact field can prove the integer difference and the
// finite Taylor slab can prove that shifting discards no nonzero tail.
template <typename Scalar>
LocalSolution<Scalar> canonicalize_identical_local_sectors(
    LocalSolution<Scalar> input) {
  validate_local_solution(input, false);
  std::vector<LocalSector<Scalar>> merged;
  for (auto& sector : input.sectors) {
    const auto found = std::find_if(merged.begin(), merged.end(),
        [&](const auto& candidate) {
          return local_algebra_detail::same_sector_tag(candidate, sector);
        });
    if (found == merged.end()) {
      merged.push_back(std::move(sector));
    } else {
      if (found->coefficients.size() != sector.coefficients.size())
        throw std::invalid_argument("identical local tags have unequal slabs");
      for (std::size_t i = 0; i < found->coefficients.size(); ++i)
        found->coefficients[i] += sector.coefficients[i];
    }
  }
  std::stable_sort(merged.begin(), merged.end(),
      local_algebra_detail::sector_less<Scalar>);
  input.sectors = std::move(merged);
  validate_local_solution(input, false);
  return input;
}

template <typename Scalar>
LocalSolution<Scalar> multiply_prepared_rational(
    const LocalSolution<Scalar>& input,
    const PreparedRationalTaylorMultiplier<Scalar>& multiplier,
    std::string checkpoint_identity = {}) {
  validate_local_solution(input, false);
  if (!input.error.empty())
    throw std::invalid_argument(
        "native rational multiplication needs explicit error-envelope propagation");
  const auto epsilon_width = input.epsilon.width();
  const auto taylor_width = input.taylor_width();
  if (multiplier.kernels.size() < epsilon_width)
    throw std::invalid_argument(
        "prepared rational multiplier has too few epsilon kernels");
  for (std::size_t j = 0; j < epsilon_width; ++j)
    if (multiplier.kernels[j].size() < taylor_width)
      throw std::invalid_argument(
          "prepared rational multiplier has too few Taylor coefficients");

  LocalSolution<Scalar> output;
  output.chart = input.chart;
  output.epsilon = {
      local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(input.epsilon.min_power) +
              multiplier.epsilon_shift,
          "rational-product epsilon minimum"),
      local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(input.epsilon.complete_max) +
              multiplier.epsilon_shift,
          "rational-product epsilon maximum")};
  output.taylor_complete_max = input.taylor_complete_max;
  output.dimension = input.dimension;
  output.prescriptions = input.prescriptions;
  output.checkpoint_identity = checkpoint_identity.empty()
      ? input.checkpoint_identity + ":mul:" + multiplier.exact_identity
      : std::move(checkpoint_identity);
  output.sectors.reserve(input.sectors.size());

  std::vector<local_algebra_detail::AcbPolynomial>
      acb_kernel_polynomials;
  std::vector<bool> acb_kernel_material;
  bool rational_taylor_constant = false;
  if constexpr (std::is_same_v<Scalar, ComplexBall>) {
    acb_kernel_polynomials.reserve(epsilon_width);
    acb_kernel_material.assign(epsilon_width, false);
    for (std::size_t epsilon = 0; epsilon < epsilon_width; ++epsilon) {
      acb_kernel_polynomials.emplace_back();
      for (std::size_t taylor = 0; taylor < taylor_width; ++taylor) {
        const auto& coefficient = multiplier.kernels[epsilon][taylor];
        if (coefficient.is_zero()) continue;
        acb_poly_set_coeff_acb(acb_kernel_polynomials.back().raw(),
                               static_cast<slong>(taylor),
                               coefficient.raw());
        acb_kernel_material[epsilon] = true;
      }
    }
  } else if constexpr (std::is_same_v<Scalar, Rational>) {
    rational_taylor_constant = std::all_of(
        multiplier.kernels.begin(),
        multiplier.kernels.begin() +
            static_cast<std::ptrdiff_t>(epsilon_width),
        [taylor_width](const auto& kernel) {
          return std::all_of(
              kernel.begin() + 1,
              kernel.begin() + static_cast<std::ptrdiff_t>(taylor_width),
              [](const auto& coefficient) { return coefficient.is_zero(); });
        });
  }

  for (const auto& sector : input.sectors) {
    LocalSector<Scalar> product;
    product.a = local_algebra_detail::subtract_nonnegative_integer(
        sector.a, multiplier.center_pole_order);
    product.b = sector.b;
    product.log_power = sector.log_power;
    product.coefficients.assign(output.sector_size(),
                                ScalarTraits<Scalar>::zero());
    if constexpr (std::is_same_v<Scalar, ComplexBall>) {
      // The old loop below performed a fresh scalar triangular Taylor
      // convolution for every output epsilon row.  At N=600 this dominates
      // every SCC coupling and gauge transform.  Build one source polynomial
      // per epsilon row/component and use Arb's capped product for the same
      // finite convolution rectangle.
      for (std::uint32_t component = 0;
           component < input.dimension; ++component) {
        std::vector<local_algebra_detail::AcbPolynomial>
            source_polynomials;
        std::vector<bool> source_material(epsilon_width, false);
        source_polynomials.reserve(epsilon_width);
        for (std::size_t epsilon = 0; epsilon < epsilon_width; ++epsilon) {
          source_polynomials.emplace_back();
          for (std::size_t taylor = 0; taylor < taylor_width; ++taylor) {
            const auto& coefficient = sector.coefficients[
                local_algebra_detail::flat_index(
                    epsilon, taylor, component, taylor_width,
                    input.dimension)];
            if (coefficient.is_zero()) continue;
            acb_poly_set_coeff_acb(source_polynomials.back().raw(),
                                   static_cast<slong>(taylor),
                                   coefficient.raw());
            source_material[epsilon] = true;
          }
        }
        local_algebra_detail::AcbPolynomial polynomial_product;
        for (std::size_t kernel_epsilon = 0;
             kernel_epsilon < epsilon_width; ++kernel_epsilon) {
          if (!acb_kernel_material[kernel_epsilon]) continue;
          for (std::size_t input_epsilon = 0;
               input_epsilon + kernel_epsilon < epsilon_width;
               ++input_epsilon) {
            if (!source_material[input_epsilon]) continue;
            acb_poly_mullow(
                polynomial_product.raw(),
                acb_kernel_polynomials[kernel_epsilon].raw(),
                source_polynomials[input_epsilon].raw(),
                static_cast<slong>(taylor_width),
                ComplexBall::precision());
            const auto output_epsilon = input_epsilon + kernel_epsilon;
            const auto product_length = std::min<std::size_t>(
                taylor_width, static_cast<std::size_t>(
                    acb_poly_length(polynomial_product.raw())));
            for (std::size_t taylor = 0;
                 taylor < product_length; ++taylor) {
              ComplexBall coefficient;
              acb_poly_get_coeff_acb(
                  coefficient.raw(), polynomial_product.raw(),
                  static_cast<slong>(taylor));
              product.coefficients[local_algebra_detail::flat_index(
                  output_epsilon, taylor, component, taylor_width,
                  input.dimension)] += coefficient;
            }
          }
        }
      }
    } else if constexpr (std::is_same_v<Scalar, Rational>) {
      if (rational_taylor_constant) {
        // CASE-P polar weights are constant in the local Taylor variable.
        // Keep the generic loop's epsilon-kernel accumulation order while
        // avoiding its triangular scan over coefficients proven to be zero.
        for (std::size_t out_ei = 0; out_ei < epsilon_width; ++out_ei) {
          for (std::size_t kernel_ei = 0; kernel_ei <= out_ei; ++kernel_ei) {
            const auto& coefficient = multiplier.kernels[kernel_ei][0];
            if (coefficient.is_zero()) continue;
            const auto input_ei = out_ei - kernel_ei;
            for (std::size_t n = 0; n < taylor_width; ++n) {
              for (std::uint32_t component = 0;
                   component < input.dimension; ++component) {
                product.coefficients[local_algebra_detail::flat_index(
                    out_ei, n, component, taylor_width, input.dimension)] +=
                    coefficient * sector.coefficients[
                        local_algebra_detail::flat_index(
                            input_ei, n, component, taylor_width,
                            input.dimension)];
              }
            }
          }
        }
      } else {
        for (std::size_t out_ei = 0; out_ei < epsilon_width; ++out_ei) {
          for (std::size_t kernel_ei = 0; kernel_ei <= out_ei; ++kernel_ei) {
            const auto input_ei = out_ei - kernel_ei;
            const auto& kernel = multiplier.kernels[kernel_ei];
            for (std::size_t n = 0; n < taylor_width; ++n) {
              for (std::size_t m = 0; m <= n; ++m) {
                if (kernel[m].is_zero()) continue;
                for (std::uint32_t component = 0;
                     component < input.dimension; ++component) {
                  product.coefficients[local_algebra_detail::flat_index(
                      out_ei, n, component, taylor_width, input.dimension)] +=
                      kernel[m] * sector.coefficients[
                          local_algebra_detail::flat_index(
                              input_ei, n - m, component, taylor_width,
                              input.dimension)];
                }
              }
            }
          }
        }
      }
    } else {
      for (std::size_t out_ei = 0; out_ei < epsilon_width; ++out_ei) {
        for (std::size_t kernel_ei = 0; kernel_ei <= out_ei; ++kernel_ei) {
          const auto input_ei = out_ei - kernel_ei;
          const auto& kernel = multiplier.kernels[kernel_ei];
          for (std::size_t n = 0; n < taylor_width; ++n) {
            for (std::size_t m = 0; m <= n; ++m) {
              if (ScalarTraits<Scalar>::is_zero(kernel[m])) continue;
              for (std::uint32_t component = 0;
                   component < input.dimension; ++component) {
                product.coefficients[local_algebra_detail::flat_index(
                    out_ei, n, component, taylor_width, input.dimension)] +=
                    kernel[m] * sector.coefficients[
                        local_algebra_detail::flat_index(
                            input_ei, n - m, component, taylor_width,
                            input.dimension)];
              }
            }
          }
        }
      }
    }
    output.sectors.push_back(std::move(product));
  }
  return canonicalize_identical_local_sectors(std::move(output));
}

namespace local_algebra_detail {

template <typename Scalar>
void add_product(Scalar& accumulator, const Scalar& left,
                 const Scalar& right) {
  accumulator += left * right;
}

template <>
inline void add_product<ComplexBall>(ComplexBall& accumulator,
                                     const ComplexBall& left,
                                     const ComplexBall& right) {
  acb_addmul(accumulator.raw(), left.raw(), right.raw(),
             ComplexBall::precision());
}

}  // namespace local_algebra_detail

// Direct specialization of a sparse one-row local matrix.  Observable
// contraction needs only one scalar local, while the generic matrix path
// selects, multiplies, embeds, and combines one complete LocalSolution per
// active column.  Accumulate the same finite convolutions directly into the
// scalar exact-tag groups instead.  A caller may cap the honest upper frame;
// coefficients above it are deliberately not computed, while the structural
// lower edge and the least complete active-entry edge remain unchanged.
template <typename Scalar>
std::optional<LocalSolution<Scalar>> apply_prepared_scalar_row_window(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& matrix,
    const LocalSolution<Scalar>& input,
    std::int32_t target_complete_max,
    std::string checkpoint_identity = {}) {
  validate_local_solution(input, false);
  if (!input.error.empty())
    throw std::invalid_argument(
        "native scalar-row application needs explicit error-envelope propagation");
  if (matrix.rows != 1 || matrix.columns != input.dimension)
    throw std::invalid_argument(
        "prepared scalar row dimensions disagree with its local");

  const auto epsilon_width = input.epsilon.width();
  const auto taylor_width = input.taylor_width();
  std::int32_t natural_min = std::numeric_limits<std::int32_t>::max();
  std::int32_t natural_complete = std::numeric_limits<std::int32_t>::max();
  bool active = false;
  for (const auto& entry : matrix.entries) {
    if (entry.row != 0 || entry.column >= matrix.columns)
      throw std::invalid_argument(
          "prepared scalar-row entry is out of range");
    if (entry.multiplier.structurally_zero()) continue;
    active = true;
    const auto multiplier_width = entry.multiplier.kernels.size();
    if (multiplier_width == 0)
      throw std::invalid_argument(
          "prepared scalar-row multiplier has no epsilon kernels");
    for (std::size_t j = 0; j < multiplier_width; ++j)
      if (entry.multiplier.kernels[j].size() < taylor_width)
        throw std::invalid_argument(
            "prepared scalar-row multiplier has too few Taylor coefficients");
    const auto term_min = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(input.epsilon.min_power) +
            entry.multiplier.epsilon_shift,
        "scalar-row epsilon minimum");
    const auto term_complete = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(term_min) +
            static_cast<std::int64_t>(
                std::min(epsilon_width, multiplier_width)) -
            1,
        "scalar-row epsilon complete maximum");
    natural_min = std::min(natural_min, term_min);
    natural_complete = std::min(
        natural_complete, term_complete);
  }
  if (!active) return std::nullopt;

  const auto complete_max = std::min(
      natural_complete, target_complete_max);
  if (complete_max < natural_min)
    throw std::invalid_argument(
        "prepared scalar row has no coefficient in its requested upper window");

  LocalSolution<Scalar> output;
  output.chart = input.chart;
  output.epsilon = {natural_min, complete_max};
  output.taylor_complete_max = input.taylor_complete_max;
  output.dimension = 1;
  output.prescriptions = input.prescriptions;
  output.checkpoint_identity = checkpoint_identity.empty()
      ? input.checkpoint_identity + ":scalar-row:" + matrix.exact_identity
      : std::move(checkpoint_identity);
  // One output tag per input sector is the common row-projection shape.
  // Reserve only that metadata scale; varying pole orders may grow beyond it,
  // but reserving the full sector-entry cross product would recreate the
  // memory problem this specialization is intended to avoid.
  output.sectors.reserve(input.sectors.size());

  const auto find_group = [&](const LocalSector<Scalar>& source,
                              std::uint32_t pole) -> LocalSector<Scalar>& {
    const auto shifted_a =
        local_algebra_detail::subtract_nonnegative_integer(source.a, pole);
    const auto found = std::find_if(
        output.sectors.begin(), output.sectors.end(),
        [&](const auto& candidate) {
          return candidate.log_power == source.log_power &&
                 local_algebra_detail::same_descriptor(candidate.a,
                                                       shifted_a) &&
                 local_algebra_detail::same_descriptor(candidate.b,
                                                       source.b);
        });
    if (found != output.sectors.end()) return *found;
    LocalSector<Scalar> group;
    group.a = shifted_a;
    group.b = source.b;
    group.log_power = source.log_power;
    group.coefficients.assign(output.sector_size(),
                              ScalarTraits<Scalar>::zero());
    output.sectors.push_back(std::move(group));
    return output.sectors.back();
  };

  for (const auto& entry : matrix.entries) {
    if (entry.multiplier.structurally_zero()) continue;
    const auto term_min = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(input.epsilon.min_power) +
            entry.multiplier.epsilon_shift,
        "scalar-row term epsilon minimum");
    const auto multiplier_width = entry.multiplier.kernels.size();
    const bool multiplier_can_contribute = [&]() {
      for (std::size_t kernel_epsilon = 0;
           kernel_epsilon < multiplier_width; ++kernel_epsilon) {
        if (static_cast<std::int64_t>(term_min) +
                static_cast<std::int64_t>(kernel_epsilon) >
            complete_max)
          break;
        const auto& kernel = entry.multiplier.kernels[kernel_epsilon];
        if (std::any_of(
                kernel.begin(), kernel.begin() + taylor_width,
                [](const auto& value) {
                  return !ScalarTraits<Scalar>::is_zero(value);
                }))
          return true;
      }
      return false;
    }();
    if (!multiplier_can_contribute) continue;
    for (const auto& sector : input.sectors) {
      const bool selected_component_can_contribute = [&]() {
        for (std::size_t input_epsilon = 0;
             input_epsilon < epsilon_width; ++input_epsilon)
          for (std::size_t input_taylor = 0;
               input_taylor < taylor_width; ++input_taylor)
            if (!ScalarTraits<Scalar>::is_zero(
                    sector.coefficients[
                        local_algebra_detail::flat_index(
                            input_epsilon, input_taylor, entry.column,
                            taylor_width, input.dimension)]))
              return true;
        return false;
      }();
      if (!selected_component_can_contribute) continue;
      LocalSector<Scalar>* group = nullptr;
      for (std::size_t kernel_epsilon = 0;
           kernel_epsilon < multiplier_width; ++kernel_epsilon) {
        const auto& kernel = entry.multiplier.kernels[kernel_epsilon];
        for (std::size_t input_epsilon = 0;
             input_epsilon + kernel_epsilon < epsilon_width;
             ++input_epsilon) {
          const auto output_power = local_algebra_detail::checked_i32(
              static_cast<std::int64_t>(term_min) +
                  static_cast<std::int64_t>(input_epsilon) +
                  static_cast<std::int64_t>(kernel_epsilon),
              "scalar-row output epsilon power");
          if (output_power > complete_max) break;
          const auto output_epsilon = static_cast<std::size_t>(
              output_power - output.epsilon.min_power);
          for (std::size_t kernel_taylor = 0;
               kernel_taylor < taylor_width; ++kernel_taylor) {
            const auto& kernel_coefficient = kernel[kernel_taylor];
            if (ScalarTraits<Scalar>::is_zero(kernel_coefficient)) continue;
            for (std::size_t input_taylor = 0;
                 input_taylor + kernel_taylor < taylor_width;
                 ++input_taylor) {
              const auto& source_coefficient = sector.coefficients[
                  local_algebra_detail::flat_index(
                      input_epsilon, input_taylor, entry.column,
                      taylor_width, input.dimension)];
              if (ScalarTraits<Scalar>::is_zero(source_coefficient))
                continue;
              // Allocate the potentially large scalar sector slab only
              // after both exact operands prove that this selected component
              // can contribute.  For Acb, is_zero means the exact singleton
              // zero; an enclosure merely containing zero remains material.
              if (group == nullptr)
                group = &find_group(
                    sector, entry.multiplier.center_pole_order);
              auto& destination = group->coefficients[
                  local_algebra_detail::flat_index(
                      output_epsilon, input_taylor + kernel_taylor, 0,
                      taylor_width, 1)];
              local_algebra_detail::add_product(
                  destination, kernel_coefficient, source_coefficient);
            }
          }
        }
      }
    }
  }

  if (output.sectors.empty()) {
    // Every active exact entry still owns the frame computed above even when
    // all of its finite kernels or selected source cells are exact zero.
    // LocalSolution requires at least one sector, so retain one genuine
    // product tag as a zero representative without materializing the full
    // input-sector x active-entry cross product.
    const auto active_entry = std::find_if(
        matrix.entries.begin(), matrix.entries.end(), [](const auto& entry) {
          return !entry.multiplier.structurally_zero();
        });
    if (active_entry == matrix.entries.end() || input.sectors.empty())
      throw std::logic_error(
          "active scalar row lost its zero-sector representative");
    (void)find_group(input.sectors.front(),
                     active_entry->multiplier.center_pole_order);
  }

  std::stable_sort(output.sectors.begin(), output.sectors.end(),
                   local_algebra_detail::sector_less<Scalar>);
  validate_local_solution(output, false);
  return output;
}

// Retain exactly the Taylor prefix used by a certified match.  A matching
// retry may intentionally reject a numerically poisoned stored suffix; that
// suffix must not silently re-enter when the receiving local is materialized.
template <typename Scalar>
LocalSolution<Scalar> restrict_local_taylor_prefix(
    const LocalSolution<Scalar>& input,
    std::uint32_t target_complete_max,
    std::string checkpoint_identity = {}) {
  validate_local_solution(input, false);
  if (!input.error.empty())
    throw std::invalid_argument(
        "native Taylor-prefix restriction needs explicit error-envelope propagation");
  if (target_complete_max > input.taylor_complete_max)
    throw std::invalid_argument(
        "requested Taylor prefix exceeds the retained local order; "
        "requested_complete_max=" +
        std::to_string(target_complete_max) +
        ";retained_complete_max=" +
        std::to_string(input.taylor_complete_max) +
        ";source_checkpoint_identity=" + input.checkpoint_identity +
        ";requested_checkpoint_identity=" + checkpoint_identity);
  if (target_complete_max == input.taylor_complete_max)
    return input;

  LocalSolution<Scalar> output;
  output.chart = input.chart;
  output.epsilon = input.epsilon;
  output.taylor_complete_max = target_complete_max;
  output.dimension = input.dimension;
  output.prescriptions = input.prescriptions;
  output.checkpoint_identity = checkpoint_identity.empty()
      ? input.checkpoint_identity + ":taylor-prefix:" +
            std::to_string(target_complete_max)
      : std::move(checkpoint_identity);
  output.sectors.reserve(input.sectors.size());
  for (const auto& sector : input.sectors) {
    LocalSector<Scalar> restricted;
    restricted.a = sector.a;
    restricted.b = sector.b;
    restricted.log_power = sector.log_power;
    restricted.coefficients.assign(output.sector_size(),
                                    ScalarTraits<Scalar>::zero());
    for (std::size_t epsilon = 0;
         epsilon < input.epsilon.width(); ++epsilon)
      for (std::size_t taylor = 0;
           taylor < output.taylor_width(); ++taylor)
        for (std::uint32_t component = 0;
             component < input.dimension; ++component)
          restricted.coefficients[local_algebra_detail::flat_index(
              epsilon, taylor, component, output.taylor_width(),
              output.dimension)] =
              sector.coefficients[local_algebra_detail::flat_index(
                  epsilon, taylor, component, input.taylor_width(),
                  input.dimension)];
    output.sectors.push_back(std::move(restricted));
  }
  validate_local_solution(output, false);
  return output;
}

// Publish factorized matching columns only through the Taylor prefix which
// the receiving basis actually used.  A singular endpoint may retain a wider
// private reservoir for its tail proof; that is safe input to this boundary,
// but it must not become part of the public retained-match contract.  A
// narrower factorized column is never repairable by padding unknown terms.
template <typename Scalar>
void restrict_factorized_local_basis_taylor_prefix(
    std::vector<LocalSolution<Scalar>>& basis,
    std::uint32_t target_complete_max,
    const std::string& context) {
  if (basis.empty())
    throw std::logic_error(
        context + ": factorized local basis is empty");
  for (std::size_t column = 0; column < basis.size(); ++column) {
    validate_local_solution(basis[column], false);
    if (basis[column].taylor_complete_max < target_complete_max)
      throw std::logic_error(
          context + ": factorized local column " +
          std::to_string(column) +
          " is narrower than the receiving Taylor prefix; "
          "factorized_complete_max=" +
          std::to_string(basis[column].taylor_complete_max) +
          ";receiving_complete_max=" +
          std::to_string(target_complete_max));
    if (basis[column].taylor_complete_max == target_complete_max)
      continue;
    basis[column] = restrict_local_taylor_prefix(
        basis[column], target_complete_max,
        basis[column].checkpoint_identity +
            ":receiving-taylor-prefix:" +
            std::to_string(target_complete_max));
  }
}

// Restrict a finite local slab to a target recurrence frame.  Discarding
// known upper coefficients is safe because the target never consumes them.
// Discarding lower coefficients is only safe when every discarded exact
// coefficient is zero: epsilon-dependent Frobenius exponents can mix a pole
// row into finite and positive orders during evaluation.
template <typename Scalar>
LocalSolution<Scalar> restrict_local_epsilon_frame_strict_lower(
    const LocalSolution<Scalar>& input, std::int32_t target_min,
    std::int32_t target_complete_max,
    std::string checkpoint_identity = {}) {
  validate_local_solution(input, false);
  if (!input.error.empty())
    throw std::invalid_argument(
        "native epsilon-frame restriction needs explicit error-envelope propagation");
  if (target_complete_max < target_min)
    throw std::invalid_argument("empty target epsilon frame");

  if (input.epsilon.min_power < target_min) {
    const auto discarded_max = std::min<std::int32_t>(
        input.epsilon.complete_max,
        local_algebra_detail::checked_i32(
            static_cast<std::int64_t>(target_min) - 1,
            "lower epsilon-frame boundary"));
    for (const auto& sector : input.sectors)
      for (std::int64_t power = input.epsilon.min_power;
           power <= discarded_max; ++power) {
        const auto epsilon_index = static_cast<std::size_t>(
            power - input.epsilon.min_power);
        for (std::size_t n = 0; n < input.taylor_width(); ++n)
          for (std::uint32_t component = 0;
               component < input.dimension; ++component)
            if (!ScalarTraits<Scalar>::is_zero(
                    sector.coefficients[local_algebra_detail::flat_index(
                        epsilon_index, n, component, input.taylor_width(),
                        input.dimension)]))
              throw std::invalid_argument(
                  "signed epsilon shift has nonzero content below the target frame; insufficient lower halo");
      }
  }

  if (input.epsilon.complete_max < target_min)
    throw std::invalid_argument(
        "signed epsilon shift leaves no complete coefficient in the target frame");

  const bool structurally_zero_on_target =
      input.epsilon.min_power > target_complete_max;
  const EpsilonWindow output_window = structurally_zero_on_target
      ? EpsilonWindow{target_min, target_complete_max}
      : EpsilonWindow{
            std::max(input.epsilon.min_power, target_min),
            std::min(input.epsilon.complete_max, target_complete_max)};

  LocalSolution<Scalar> output;
  output.chart = input.chart;
  output.epsilon = output_window;
  output.taylor_complete_max = input.taylor_complete_max;
  output.dimension = input.dimension;
  output.prescriptions = input.prescriptions;
  output.checkpoint_identity = checkpoint_identity.empty()
      ? input.checkpoint_identity + ":epsilon-frame"
      : std::move(checkpoint_identity);
  output.sectors.reserve(input.sectors.size());
  for (const auto& sector : input.sectors) {
    LocalSector<Scalar> restricted;
    restricted.a = sector.a;
    restricted.b = sector.b;
    restricted.log_power = sector.log_power;
    restricted.coefficients.assign(output.sector_size(),
                                    ScalarTraits<Scalar>::zero());
    if (!structurally_zero_on_target) {
      for (std::int64_t power = output.epsilon.min_power;
           power <= output.epsilon.complete_max; ++power) {
        const auto input_epsilon = static_cast<std::size_t>(
            power - input.epsilon.min_power);
        const auto output_epsilon = static_cast<std::size_t>(
            power - output.epsilon.min_power);
        for (std::size_t n = 0; n < output.taylor_width(); ++n)
          for (std::uint32_t component = 0;
               component < output.dimension; ++component)
            restricted.coefficients[local_algebra_detail::flat_index(
                output_epsilon, n, component, output.taylor_width(),
                output.dimension)] = sector.coefficients[
                    local_algebra_detail::flat_index(
                        input_epsilon, n, component, input.taylor_width(),
                        input.dimension)];
      }
    }
    output.sectors.push_back(std::move(restricted));
  }
  validate_local_solution(output, false);
  return output;
}

template <typename Scalar>
LocalSolution<Scalar> combine_local_solutions(
    const std::vector<LocalSolution<Scalar>>& inputs,
    std::string checkpoint_identity = {}) {
  if (inputs.empty())
    throw std::invalid_argument("cannot combine an empty local-solution list");
  for (const auto& input : inputs) {
    validate_local_solution(input, false);
    if (!input.error.empty())
      throw std::invalid_argument(
          "native local combination needs explicit error-envelope propagation");
    local_algebra_detail::require_same_local_space(inputs.front(), input);
  }

  const auto min_power = std::min_element(inputs.begin(), inputs.end(),
      [](const auto& left, const auto& right) {
        return left.epsilon.min_power < right.epsilon.min_power;
      })->epsilon.min_power;
  const auto complete_max = std::min_element(inputs.begin(), inputs.end(),
      [](const auto& left, const auto& right) {
        return left.epsilon.complete_max < right.epsilon.complete_max;
      })->epsilon.complete_max;
  if (complete_max < min_power)
    throw std::invalid_argument("combined local solutions have an empty epsilon window");
  const auto taylor_complete_max = std::min_element(
      inputs.begin(), inputs.end(), [](const auto& left, const auto& right) {
        return left.taylor_complete_max < right.taylor_complete_max;
      })->taylor_complete_max;

  LocalSolution<Scalar> output;
  output.chart = inputs.front().chart;
  output.epsilon = {min_power, complete_max};
  output.taylor_complete_max = taylor_complete_max;
  output.dimension = inputs.front().dimension;
  output.prescriptions = inputs.front().prescriptions;
  output.checkpoint_identity = checkpoint_identity.empty()
      ? inputs.front().checkpoint_identity + ":combined"
      : std::move(checkpoint_identity);

  const auto out_taylor_width = output.taylor_width();
  for (const auto& input : inputs) {
    for (const auto& sector : input.sectors) {
      LocalSector<Scalar> aligned;
      aligned.a = sector.a;
      aligned.b = sector.b;
      aligned.log_power = sector.log_power;
      aligned.coefficients.assign(output.sector_size(),
                                  ScalarTraits<Scalar>::zero());
      const auto copy_min = std::max(min_power, input.epsilon.min_power);
      const auto copy_max = std::min(complete_max, input.epsilon.complete_max);
      for (std::int32_t power = copy_min; power <= copy_max; ++power) {
        const auto input_ei = static_cast<std::size_t>(
            power - input.epsilon.min_power);
        const auto output_ei = static_cast<std::size_t>(power - min_power);
        for (std::size_t n = 0; n < out_taylor_width; ++n)
          for (std::uint32_t component = 0;
               component < output.dimension; ++component)
            aligned.coefficients[local_algebra_detail::flat_index(
                output_ei, n, component, out_taylor_width,
                output.dimension)] = sector.coefficients[
                    local_algebra_detail::flat_index(
                        input_ei, n, component, input.taylor_width(),
                        input.dimension)];
      }
      output.sectors.push_back(std::move(aligned));
    }
  }
  return canonicalize_identical_local_sectors(std::move(output));
}

// Apply one exact epsilon-Laurent right transformation to a retained local
// basis without first expanding the corresponding physical weights.
//
// Matching in a saturated lattice solves
//
//                         (F T) y = b.
//
// Reconstructing w = T y and then materializing F w is algebraically
// equivalent, but for Acb it recreates every confluent-column cancellation
// independently in every Taylor coefficient.  Keeping the association
// (F T) y preserves the exact structural dependencies used by the matching
// solve.  Omitted terms of ExactLaurentPolynomial are structural zero, so
// unlike a finite EpsilonFrame they do not impose an upper completeness edge.
template <typename Scalar>
std::vector<LocalSolution<Scalar>> right_transform_local_basis_exact(
    const std::vector<const LocalSolution<Scalar>*>& basis,
    const ExactLaurentMatrix<Scalar>& transformation,
    const std::string& checkpoint_identity_root) {
  const auto dimension = basis.size();
  if (dimension == 0 || transformation.size() != dimension ||
      checkpoint_identity_root.empty())
    throw std::invalid_argument(
        "exact right local-basis transformation requires a nonempty square basis and checkpoint root");
  for (const auto& row : transformation)
    if (row.size() != dimension)
      throw std::invalid_argument(
          "exact right local-basis transformation is not square");
  for (const auto* column : basis) {
    if (column == nullptr)
      throw std::invalid_argument(
          "exact right local-basis transformation received a null column");
    validate_local_solution(*column, false);
    if (!column->error.empty())
      throw std::invalid_argument(
          "exact right local-basis transformation cannot discard an error envelope");
    local_algebra_detail::require_same_local_space(*basis.front(), *column);
  }

  std::vector<LocalSolution<Scalar>> transformed;
  transformed.reserve(dimension);
  for (std::size_t output_column = 0; output_column < dimension;
       ++output_column) {
    std::optional<std::int32_t> minimum;
    std::optional<std::int32_t> complete_max;
    std::uint32_t taylor_complete_max =
        std::numeric_limits<std::uint32_t>::max();
    for (std::size_t input_column = 0; input_column < dimension;
         ++input_column) {
      const auto& polynomial =
          transformation[input_column][output_column];
      const auto polynomial_minimum = polynomial.minimum_power();
      if (!polynomial_minimum.has_value()) continue;
      const auto term_minimum = local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(
              basis[input_column]->epsilon.min_power) +
              *polynomial_minimum,
          "exact right transformed local epsilon minimum");
      const auto term_complete = local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(
              basis[input_column]->epsilon.complete_max) +
              *polynomial_minimum,
          "exact right transformed local complete maximum");
      minimum = !minimum.has_value()
          ? term_minimum : std::min(*minimum, term_minimum);
      complete_max = !complete_max.has_value()
          ? term_complete : std::min(*complete_max, term_complete);
      taylor_complete_max = std::min(
          taylor_complete_max,
          basis[input_column]->taylor_complete_max);
    }
    if (!minimum.has_value() || !complete_max.has_value() ||
        *complete_max < *minimum ||
        taylor_complete_max ==
            std::numeric_limits<std::uint32_t>::max())
      throw std::invalid_argument(
          "exact right local-basis transformation has an empty or invalid output column");

    LocalSolution<Scalar> output;
    output.chart = basis.front()->chart;
    output.epsilon = {*minimum, *complete_max};
    output.taylor_complete_max = taylor_complete_max;
    output.dimension = basis.front()->dimension;
    output.prescriptions = basis.front()->prescriptions;
    output.checkpoint_identity =
        checkpoint_identity_root + ":column:" +
        std::to_string(output_column);

    const auto output_taylor_width = output.taylor_width();
    for (std::size_t input_column = 0; input_column < dimension;
         ++input_column) {
      const auto& polynomial =
          transformation[input_column][output_column];
      if (polynomial.is_zero()) continue;
      const auto& source = *basis[input_column];
      for (const auto& sector : source.sectors) {
        LocalSector<Scalar> materialized;
        materialized.a = sector.a;
        materialized.b = sector.b;
        materialized.log_power = sector.log_power;
        materialized.coefficients.assign(
            output.sector_size(), ScalarTraits<Scalar>::zero());
        for (std::int64_t output_power = output.epsilon.min_power;
             output_power <= output.epsilon.complete_max;
             ++output_power) {
          const auto output_epsilon = static_cast<std::size_t>(
              output_power - output.epsilon.min_power);
          for (const auto& [transformation_power,
                            transformation_coefficient] :
               polynomial.terms()) {
            const auto source_power =
                output_power - transformation_power;
            if (source_power < source.epsilon.min_power ||
                source_power > source.epsilon.complete_max)
              continue;
            const auto source_epsilon = static_cast<std::size_t>(
                source_power - source.epsilon.min_power);
            for (std::size_t n = 0; n < output_taylor_width; ++n)
              for (std::uint32_t component = 0;
                   component < output.dimension; ++component)
                materialized.coefficients[
                    local_algebra_detail::flat_index(
                        output_epsilon, n, component,
                        output_taylor_width, output.dimension)] +=
                    transformation_coefficient *
                    sector.coefficients[
                        local_algebra_detail::flat_index(
                            source_epsilon, n, component,
                            source.taylor_width(), source.dimension)];
          }
        }
        output.sectors.push_back(std::move(materialized));
      }
    }
    output = canonicalize_identical_local_sectors(std::move(output));
    std::optional<std::int32_t> first_material_power;
    for (std::int64_t power = output.epsilon.min_power;
         power <= output.epsilon.complete_max &&
             !first_material_power.has_value();
         ++power) {
      const auto epsilon_index = static_cast<std::size_t>(
          power - output.epsilon.min_power);
      for (const auto& sector : output.sectors) {
        const auto begin = epsilon_index * output.taylor_width() *
            output.dimension;
        const auto end = begin + output.taylor_width() *
            output.dimension;
        if (std::any_of(
                sector.coefficients.begin() +
                    static_cast<std::ptrdiff_t>(begin),
                sector.coefficients.begin() +
                    static_cast<std::ptrdiff_t>(end),
                [](const Scalar& value) {
                  return !ScalarTraits<Scalar>::is_zero(value);
                })) {
          first_material_power = static_cast<std::int32_t>(power);
          break;
        }
      }
    }
    if (first_material_power.has_value() &&
        *first_material_power > output.epsilon.min_power)
      output = restrict_local_epsilon_frame_strict_lower(
          output, *first_material_power, output.epsilon.complete_max,
          output.checkpoint_identity);
    transformed.push_back(std::move(output));
  }
  return transformed;
}

// Materialize one retained matching state in its receiving chart without
// exporting either the basis slabs or the Laurent weights.  Each basis column
// is a full local vector S_j(t,eps), while weights[j] is a finite Laurent
// frame independent of t.  The complete upper edge is the first edge at which
// an unseen coefficient of either factor could contribute:
//
//   CompleteMax(S_j w_j) =
//     min(CompleteMax(S_j) + Min(w_j), Min(S_j) + CompleteMax(w_j)).
//
// Lower frame edges are exact structural bounds in both representations.
// No zero padding above a complete edge, midpoint zero test, or coefficient
// export is involved.
template <typename Scalar>
LocalSolution<Scalar> materialize_local_basis_weights(
    const std::vector<const LocalSolution<Scalar>*>& basis,
    const FiniteLaurentVector<Scalar>& weights,
    std::string checkpoint_identity) {
  if (basis.empty() || basis.size() != weights.size())
    throw std::invalid_argument(
        "local basis materialization requires one weight per nonempty basis column");
  if (checkpoint_identity.empty())
    throw std::invalid_argument(
        "local basis materialization checkpoint identity is empty");
  for (const auto* column : basis) {
    if (column == nullptr)
      throw std::invalid_argument(
          "local basis materialization received a null basis column");
    validate_local_solution(*column, false);
    if (!column->error.empty())
      throw std::invalid_argument(
          "local basis materialization cannot discard an error envelope");
    local_algebra_detail::require_same_local_space(*basis.front(), *column);
  }

  std::vector<EpsilonWindow> product_windows;
  product_windows.reserve(basis.size());
  for (std::size_t column = 0; column < basis.size(); ++column) {
    const auto minimum = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(basis[column]->epsilon.min_power) +
            weights[column].min_power(),
        "materialized local epsilon minimum");
    const auto basis_complete = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(basis[column]->epsilon.complete_max) +
            weights[column].min_power(),
        "materialized local basis-complete edge");
    const auto weight_complete = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(basis[column]->epsilon.min_power) +
            weights[column].complete_max(),
        "materialized local weight-complete edge");
    product_windows.push_back(
        {minimum, std::min(basis_complete, weight_complete)});
  }
  const auto minimum = std::min_element(
      product_windows.begin(), product_windows.end(),
      [](const auto& left, const auto& right) {
        return left.min_power < right.min_power;
      })->min_power;
  const auto complete_max = std::min_element(
      product_windows.begin(), product_windows.end(),
      [](const auto& left, const auto& right) {
        return left.complete_max < right.complete_max;
      })->complete_max;
  if (complete_max < minimum)
    throw std::invalid_argument(
        "materialized local basis has no common complete epsilon window");

  LocalSolution<Scalar> output;
  output.chart = basis.front()->chart;
  output.epsilon = {minimum, complete_max};
  output.taylor_complete_max = (*std::min_element(
      basis.begin(), basis.end(), [](const auto* left, const auto* right) {
        return left->taylor_complete_max < right->taylor_complete_max;
      }))->taylor_complete_max;
  output.dimension = basis.front()->dimension;
  output.prescriptions = basis.front()->prescriptions;
  output.checkpoint_identity = std::move(checkpoint_identity);

  const auto output_taylor_width = output.taylor_width();
  for (std::size_t column = 0; column < basis.size(); ++column) {
    const auto& source = *basis[column];
    const auto& weight = weights[column];
    for (const auto& sector : source.sectors) {
      LocalSector<Scalar> materialized;
      materialized.a = sector.a;
      materialized.b = sector.b;
      materialized.log_power = sector.log_power;
      materialized.coefficients.assign(
          output.sector_size(), ScalarTraits<Scalar>::zero());
      for (std::int64_t power = output.epsilon.min_power;
           power <= output.epsilon.complete_max; ++power) {
        const auto output_epsilon = static_cast<std::size_t>(
            power - output.epsilon.min_power);
        for (std::int64_t weight_power = weight.min_power();
             weight_power <= weight.complete_max(); ++weight_power) {
          const auto source_power = power - weight_power;
          if (source_power < source.epsilon.min_power ||
              source_power > source.epsilon.complete_max)
            continue;
          const auto source_epsilon = static_cast<std::size_t>(
              source_power - source.epsilon.min_power);
          const auto& scalar_weight = weight.coefficient(
              static_cast<std::int32_t>(weight_power));
          if (ScalarTraits<Scalar>::is_zero(scalar_weight)) continue;
          for (std::size_t n = 0; n < output_taylor_width; ++n)
            for (std::uint32_t component = 0;
                 component < output.dimension; ++component)
              materialized.coefficients[local_algebra_detail::flat_index(
                  output_epsilon, n, component, output_taylor_width,
                  output.dimension)] += scalar_weight * sector.coefficients[
                      local_algebra_detail::flat_index(
                          source_epsilon, n, component,
                          source.taylor_width(), source.dimension)];
        }
      }
      output.sectors.push_back(std::move(materialized));
    }
  }
  return canonicalize_identical_local_sectors(std::move(output));
}

template <typename Scalar>
std::optional<LocalSolution<Scalar>> apply_prepared_sparse_local_matrix(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& matrix,
    const LocalSolution<Scalar>& input,
    std::string checkpoint_identity = {}) {
  validate_local_solution(input, false);
  if (!input.error.empty())
    throw std::invalid_argument(
        "native local matrix application needs explicit error-envelope propagation");
  if (matrix.rows == 0 || matrix.columns != input.dimension)
    throw std::invalid_argument("prepared local matrix dimensions disagree");
  std::vector<LocalSolution<Scalar>> terms;
  terms.reserve(matrix.entries.size());
  for (const auto& entry : matrix.entries) {
    if (entry.row >= matrix.rows || entry.column >= matrix.columns)
      throw std::invalid_argument("prepared local matrix entry is out of range");
    if (entry.multiplier.structurally_zero()) continue;
    auto selected = local_algebra_detail::with_selected_component(
        input, entry.column);
    auto product = multiply_prepared_rational(
        selected, entry.multiplier,
        input.checkpoint_identity + ":matrix-entry:" +
            std::to_string(entry.row) + ":" +
            std::to_string(entry.column));
    // Keep a zero finite slab from an exact nonzero matrix entry: in the
    // Wolfram algebra it still constrains the honest intersection window
    // when another entry of the same matrix product is nonzero.  Only the
    // fully combined matrix result may be discarded as a structural source.
    terms.push_back(local_algebra_detail::embedded_component(
        product, entry.row, matrix.rows));
  }
  if (terms.empty()) return std::nullopt;
  return combine_local_solutions(terms,
      checkpoint_identity.empty()
          ? input.checkpoint_identity + ":matrix:" + matrix.exact_identity
          : std::move(checkpoint_identity));
}

}  // namespace diffexp2
