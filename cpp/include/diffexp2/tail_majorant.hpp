#pragma once

#include "diffexp2/line_integration.hpp"
#include "diffexp2/local_algebra.hpp"
#include "diffexp2/matching.hpp"
#include "diffexp2/physical_ode.hpp"
#include "diffexp2/recurrence.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace diffexp2 {

// A full-function certificate is intentionally a different object from the
// stored-truncation residual.  This first rigorous slice covers an ordinary
// homogeneous regular chart
//
//        q(t) theta f(t) = N(t) f(t),       N(0) = 0,
//
// with an identity assembly and an epsilon-decoupled retained operator.  On a
// witness disk |t| <= R we prove q != 0 by a triangle lower bound, use
// Gronwall to bound f on the circle, and then use Cauchy's coefficient bound
// for the unseen Taylor tail.  Nothing in this layer promotes singular,
// logarithmic, sourced, or unresolved epsilon-coupled data to a certificate.

enum class TailMajorantStatus : std::uint8_t {
  Certified,
  Inconclusive,
  Unsupported
};

struct RegularTaylorTailModel {
  EpsilonWindow epsilon;
  std::uint32_t dimension = 0;
  std::uint32_t taylor_complete_max = 0;
  std::vector<ComplexBall> q_coefficients;  // exact/Acb enclosure by t lag
  // Retain the actual prepared numerator matrices, not only their norms.
  // Match materialization may share a tail theorem only after proving that
  // every receiving column was certified against this same immutable
  // operator payload.  Each lag is flattened in row-major order.
  std::vector<std::vector<ComplexBall>> n_coefficients;
  std::vector<Magnitude> n_row_sum_upper;   // infinity norm by t lag
  std::vector<Magnitude> initial_row_upper; // one bound per epsilon row
  ChartGeometry chart;
  std::vector<Prescription> prescriptions;
  std::string operator_identity;
  std::string local_checkpoint_identity;
  std::string provenance;
};

struct RegularTaylorTailModelResult {
  TailMajorantStatus status = TailMajorantStatus::Unsupported;
  std::optional<RegularTaylorTailModel> model;
  std::string detail;
};

struct RegularTaylorDiskCertificate {
  TailMajorantStatus status = TailMajorantStatus::Inconclusive;
  std::string witness_radius_exact;
  Magnitude q_lower = Magnitude::zero();
  Magnitude ode_norm_upper = Magnitude::zero();
  std::vector<Magnitude> cauchy_circle_upper;
  std::string detail;
};

struct RegularTaylorPointTailCertificate {
  TailMajorantStatus status = TailMajorantStatus::Inconclusive;
  ErrorEnvelope value;
  ErrorEnvelope theta;
  RegularTaylorDiskCertificate disk;
  std::string detail;
};

struct RegularTaylorLineTailCertificate {
  TailMajorantStatus status = TailMajorantStatus::Inconclusive;
  ErrorEnvelope integral;
  RegularTaylorDiskCertificate disk;
  std::string detail;
};

struct CertifiedLocalEvaluation {
  LocalEvaluation evaluation;
  RegularTaylorPointTailCertificate tail;
};

struct CertifiedStoredLineIntegral {
  StoredLineIntegral integral;
  RegularTaylorLineTailCertificate tail;
};

// Composite SCC locals satisfy the retained physical equation
//
//     Q(t, eps) theta F(t, eps) = C(t, eps) F(t, eps)
//
// on a finite epsilon window.  Unlike RegularTaylorTailModel, the epsilon
// rows need not decouple.  This transient terminal model treats the whole
// finite epsilon stack as one enlarged ordinary ODE, proves the scalar causal
// Q(t,epsilon) invertible on a witness disk by interval evaluation and finite
// triangular back-substitution, and reconstructs the stored Taylor prefix
// from the authoritative physical equation.  It is deliberately not a
// checkpoint payload: all fields are recomputed from the strongly held
// physical owner whenever a terminal observable needs the certificate.
struct PhysicalRegularTaylorTailModel {
  EpsilonWindow epsilon;
  std::uint32_t dimension = 0;
  std::uint32_t taylor_complete_max = 0;
  // Entry k bounds the invariant causal epsilon prefix
  // [epsilon.min_power, epsilon.min_power + k].  A single norm over the
  // complete private reservoir can be arbitrarily larger and must not be
  // assigned to lower public coefficients that cannot depend on it.
  std::vector<Magnitude> q0_inverse_prefix_norm_upper;
  std::vector<std::vector<Magnitude>>
      q_operator_prefix_norm_upper;
  // Scalar multiplication by q(t,epsilon) is lower-triangular Toeplitz on
  // every finite causal epsilon frame.  Entry [t_lag][epsilon_lag] is the
  // corresponding first-column coefficient.  Retaining the coefficients,
  // rather than only a norm of Q(t)-Q(0), lets the disk proof respect the
  // nilpotent epsilon structure instead of rejecting invertible operators
  // through a coarse Neumann condition.
  std::vector<std::vector<ComplexBall>>
      q_causal_coefficients;
  std::vector<std::vector<Magnitude>>
      c_operator_prefix_norm_upper;
  // Entry [t_lag] is a causal epsilon series of physical matrices, flattened
  // as [epsilon_lag][row][column].  It supports interval evaluation of
  // Q(t,epsilon)^-1 C(t,epsilon)/t before taking norms, preserving the large
  // scalar-clearing cancellations that separate Q/C norm bounds destroy.
  std::vector<std::vector<ComplexBall>>
      c_causal_coefficients;
  std::vector<Magnitude> initial_row_upper;
  ChartGeometry chart;
  std::vector<Prescription> prescriptions;
  std::string physical_payload_identity;
  std::string local_checkpoint_identity;
  std::string provenance;
  LocalSolution<ComplexBall> reconstructed;
};

struct PhysicalRegularTaylorTailModelResult {
  TailMajorantStatus status = TailMajorantStatus::Unsupported;
  std::optional<PhysicalRegularTaylorTailModel> model;
  std::string detail;
};

// A rational-row projection is generally not a solution of the source
// homogeneous q/C system, so it must not inherit that equation owner.  This
// separate model retains only what is needed to certify the scalar line tail:
// a source full-local Cauchy model plus independently verified analytic
// rational completions for every active multiplier epsilon coefficient.
template <typename Scalar>
struct RationalRowLineTailEntryModel {
  std::uint32_t column = 0;
  std::int32_t epsilon_shift = 0;
  std::uint32_t center_pole_order = 0;
  std::vector<PreparedRationalAnalyticCoefficient<Scalar>>
      analytic_coefficients;
  std::string exact_identity;
};

template <typename Scalar>
struct RationalRowLineTailModel {
  EpsilonWindow epsilon;
  std::uint32_t taylor_complete_max = 0;
  ChartGeometry chart;
  std::vector<Prescription> prescriptions;
  RegularTaylorTailModel source;
  std::vector<RationalRowLineTailEntryModel<Scalar>> entries;
  std::string row_exact_identity;
  std::string source_checkpoint_identity;
  std::string local_checkpoint_identity;
  std::string provenance;
};

template <typename Scalar>
struct RationalRowLineTailModelResult {
  TailMajorantStatus status = TailMajorantStatus::Unsupported;
  std::optional<RationalRowLineTailModel<Scalar>> model;
  std::string detail;
};

namespace tail_majorant_detail {

inline Rational abs_rational(const Rational& value) {
  return value.sign() < 0 ? -value : value;
}

template <typename Scalar>
bool exact_scalar_zero(const Scalar& value) {
  return ScalarTraits<Scalar>::is_zero(value);
}

inline bool finite_scalar(const Rational&) { return true; }

inline bool finite_scalar(const ComplexBall& value) {
  return value.is_finite();
}

inline std::string status_prefix(TailMajorantStatus status) {
  switch (status) {
    case TailMajorantStatus::Certified:
      return "certified";
    case TailMajorantStatus::Inconclusive:
      return "inconclusive";
    case TailMajorantStatus::Unsupported:
      return "unsupported";
  }
  throw std::logic_error("unknown tail-majorant status");
}

inline RegularTaylorTailModelResult unsupported_model(std::string detail) {
  return {TailMajorantStatus::Unsupported, std::nullopt, std::move(detail)};
}

inline RegularTaylorTailModelResult inconclusive_model(std::string detail) {
  return {TailMajorantStatus::Inconclusive, std::nullopt, std::move(detail)};
}

template <typename Scalar>
RationalRowLineTailModelResult<Scalar> unsupported_rational_row_model(
    std::string detail) {
  return {TailMajorantStatus::Unsupported, std::nullopt, std::move(detail)};
}

template <typename Scalar>
RationalRowLineTailModelResult<Scalar> inconclusive_rational_row_model(
    std::string detail) {
  return {TailMajorantStatus::Inconclusive, std::nullopt, std::move(detail)};
}

inline bool same_epsilon_window(const EpsilonWindow& left,
                                const EpsilonWindow& right) {
  return left.min_power == right.min_power &&
         left.complete_max == right.complete_max;
}

inline bool same_chart_geometry(const ChartGeometry& left,
                                const ChartGeometry& right) {
  const bool exact_radius_available =
      !left.radius_exact.empty() && !right.radius_exact.empty();
  return left.center_exact == right.center_exact &&
      left.scale_exact == right.scale_exact &&
      left.infinite_radius == right.infinite_radius &&
      (left.infinite_radius ||
       (exact_radius_available
            ? left.radius_exact == right.radius_exact
            : acb_equal(left.radius.raw(), right.radius.raw())));
}

inline bool same_prescriptions(const std::vector<Prescription>& left,
                               const std::vector<Prescription>& right) {
  if (left.size() != right.size()) return false;
  for (std::size_t index = 0; index < left.size(); ++index)
    if (left[index].factor_exact != right[index].factor_exact ||
        left[index].sign != right[index].sign ||
        left[index].multiplicity != right[index].multiplicity ||
        left[index].leading_coefficient_sign !=
            right[index].leading_coefficient_sign)
      return false;
  return true;
}

inline bool same_ball(const ComplexBall& left, const ComplexBall& right) {
  return acb_equal(left.raw(), right.raw());
}

inline bool same_operator_payload(const RegularTaylorTailModel& left,
                                  const RegularTaylorTailModel& right) {
  if (left.operator_identity != right.operator_identity ||
      left.dimension != right.dimension ||
      left.q_coefficients.size() != right.q_coefficients.size() ||
      left.n_coefficients.size() != right.n_coefficients.size() ||
      !same_chart_geometry(left.chart, right.chart) ||
      !same_prescriptions(left.prescriptions, right.prescriptions))
    return false;
  for (std::size_t lag = 0; lag < left.q_coefficients.size(); ++lag)
    if (!same_ball(left.q_coefficients[lag], right.q_coefficients[lag]))
      return false;
  for (std::size_t lag = 0; lag < left.n_coefficients.size(); ++lag) {
    if (left.n_coefficients[lag].size() !=
        right.n_coefficients[lag].size())
      return false;
    for (std::size_t entry = 0;
         entry < left.n_coefficients[lag].size(); ++entry)
      if (!same_ball(left.n_coefficients[lag][entry],
                     right.n_coefficients[lag][entry]))
        return false;
  }
  return true;
}

template <typename Scalar>
std::vector<Scalar> aggregate_nhat_lag(const PreparedLag<Scalar>& lag,
                                       std::uint32_t dimension,
                                       bool* supported,
                                       std::string* detail) {
  std::vector<Scalar> matrix(
      static_cast<std::size_t>(dimension) * dimension,
      ScalarTraits<Scalar>::zero());
  if (!lag.rational.empty()) {
    *supported = false;
    *detail =
        "regular tail majorant does not yet prove epsilon-rational matrix "
        "groups are globally decoupled";
    return matrix;
  }
  for (const auto& shifted : lag.polynomial) {
    if (shifted.shift != 0) {
      *supported = false;
      *detail =
          "regular tail majorant requires an exactly epsilon-decoupled "
          "operator (every retained matrix shift must be zero)";
      return matrix;
    }
    for (const auto& entry : shifted.entries) {
      if (entry.row >= dimension || entry.col >= dimension)
        throw std::invalid_argument(
            "tail-majorant matrix entry lies outside its dimension");
      matrix[static_cast<std::size_t>(entry.row) * dimension + entry.col] +=
          entry.value;
    }
  }
  return matrix;
}

template <typename Scalar>
Magnitude matrix_infinity_norm_upper(const std::vector<Scalar>& matrix,
                                      std::uint32_t dimension) {
  auto maximum = Magnitude::zero();
  for (std::uint32_t row = 0; row < dimension; ++row) {
    auto sum = Magnitude::zero();
    for (std::uint32_t column = 0; column < dimension; ++column)
      sum += Magnitude::upper_abs(local_detail::to_ball(
          matrix[static_cast<std::size_t>(row) * dimension + column]));
    maximum = Magnitude::maximum(maximum, sum);
  }
  return maximum;
}

inline bool encloses_recurrence_value(const Rational& stored,
                                      const Rational& expected) {
  return stored == expected;
}

inline bool encloses_recurrence_value(const ComplexBall& stored,
                                      const ComplexBall& expected) {
  return acb_contains(stored.raw(), expected.raw());
}

inline bool overlaps_recurrence_value(const Rational& stored,
                                      const Rational& expected) {
  return stored == expected;
}

inline bool overlaps_recurrence_value(const ComplexBall& stored,
                                      const ComplexBall& expected) {
  return acb_overlaps(stored.raw(), expected.raw());
}

inline bool divisor_contains_zero(const Rational& value) {
  return value.is_zero();
}

inline bool divisor_contains_zero(const ComplexBall& value) {
  return value.contains_zero();
}

inline RegularTaylorDiskCertificate inconclusive_disk(
    const std::string& radius, std::string detail) {
  RegularTaylorDiskCertificate result;
  result.status = TailMajorantStatus::Inconclusive;
  result.witness_radius_exact = radius;
  result.detail = std::move(detail);
  return result;
}

inline Magnitude exact_rational_upper(const Rational& value) {
  return Magnitude::upper_abs(ComplexBall::from_strings(value.str()));
}

inline Magnitude exact_rational_lower(const Rational& value) {
  return Magnitude::lower_abs(ComplexBall::from_strings(value.str()));
}

inline void require_model_binding(const RegularTaylorTailModel& model,
                                  const auto& solution) {
  if (model.local_checkpoint_identity != solution.checkpoint_identity ||
      model.dimension != solution.dimension ||
      model.taylor_complete_max != solution.taylor_complete_max ||
      !same_epsilon_window(model.epsilon, solution.epsilon) ||
      !same_chart_geometry(model.chart, solution.chart) ||
      !same_prescriptions(model.prescriptions, solution.prescriptions))
    throw std::invalid_argument(
        "tail-majorant model is not bound to this retained local solution");
}

template <typename Scalar>
void validate_restored_regular_taylor_tail_model(
    const RegularTaylorTailModel& model,
    const LocalSolution<Scalar>& solution,
    const std::string& expected_operator_identity) {
  static_assert(std::is_same_v<Scalar, Rational> ||
                    std::is_same_v<Scalar, ComplexBall>,
                "restored regular tail models support Rational or Acb "
                "locals only");
  validate_local_solution(solution, false);
  require_model_binding(model, solution);
  if (expected_operator_identity.empty() ||
      model.operator_identity != expected_operator_identity ||
      model.provenance.empty())
    throw std::invalid_argument(
        "restored tail model lost its retained operator/provenance binding");
  if (!solution.error.empty() || solution.sectors.size() != 1 ||
      solution.sectors.front().log_power != 0 ||
      solution.sectors.front().a.is_zero != TruthValue::Yes ||
      solution.sectors.front().b.is_zero != TruthValue::Yes)
    throw std::invalid_argument(
        "restored tail model is attached to a non-ordinary local solution");
  if (model.q_coefficients.empty() || model.n_coefficients.empty() ||
      model.n_row_sum_upper.size() != model.n_coefficients.size() ||
      model.initial_row_upper.size() != solution.epsilon.width())
    throw std::invalid_argument(
        "restored tail model has incomplete q/N or magnitude payloads");
  if (std::any_of(model.q_coefficients.begin(),
                  model.q_coefficients.end(),
                  [](const ComplexBall& value) {
                    return !value.is_finite();
                  }) ||
      model.q_coefficients.front().contains_zero())
    throw std::invalid_argument(
        "restored tail denominator payload is nonfinite or contains zero at the center");

  const auto matrix_size = static_cast<std::size_t>(solution.dimension) *
                           solution.dimension;
  for (std::size_t lag = 0; lag < model.n_coefficients.size(); ++lag) {
    const auto& matrix = model.n_coefficients[lag];
    if (matrix.size() != matrix_size ||
        std::any_of(matrix.begin(), matrix.end(),
                    [](const ComplexBall& value) {
                      return !value.is_finite();
                    }))
      throw std::invalid_argument(
          "restored tail numerator matrix payload is malformed");
    if (lag == 0 &&
        std::any_of(matrix.begin(), matrix.end(),
                    [](const ComplexBall& value) {
                      return !value.is_zero();
                    }))
      throw std::invalid_argument(
          "restored ordinary tail model has nonzero N(0)");
    const auto recomputed = matrix_infinity_norm_upper(
        matrix, solution.dimension);
    if (recomputed.dump_exact() !=
        model.n_row_sum_upper[lag].dump_exact())
      throw std::invalid_argument(
          "restored tail numerator norm does not equal its exact matrix norm");
  }

  const auto& sector = solution.sectors.front();
  for (std::size_t epsilon = 0; epsilon < solution.epsilon.width();
       ++epsilon) {
    auto recomputed_initial = Magnitude::zero();
    for (std::uint32_t component = 0; component < solution.dimension;
         ++component)
      recomputed_initial = Magnitude::maximum(
          recomputed_initial,
          Magnitude::upper_abs(local_detail::to_ball(
              sector.coefficients[local_detail::sector_index(
                  solution, epsilon, 0, component)])));
    if (recomputed_initial.dump_exact() !=
        model.initial_row_upper[epsilon].dump_exact())
      throw std::invalid_argument(
          "restored tail initial-row magnitude disagrees with its local tensor");
  }

  // Recheck the complete stored Taylor slab directly against the serialized
  // epsilon-decoupled ODE payload.  A zero-containing interval residual is
  // required at every retained coefficient; metadata and magnitudes alone
  // can never reactivate a certificate.
  for (std::size_t epsilon = 0; epsilon < solution.epsilon.width();
       ++epsilon)
    for (std::uint32_t n = 0; n <= solution.taylor_complete_max; ++n)
      for (std::uint32_t row = 0; row < solution.dimension; ++row) {
        auto residual = ComplexBall(0);
        const auto max_q = std::min<std::size_t>(
            n, model.q_coefficients.size() - 1);
        for (std::size_t lag = 0; lag <= max_q; ++lag) {
          const auto degree = n - static_cast<std::uint32_t>(lag);
          residual += model.q_coefficients[lag] *
              ComplexBall(static_cast<long>(degree)) *
              local_detail::to_ball(sector.coefficients[
                  local_detail::sector_index(
                      solution, epsilon, degree, row)]);
        }
        const auto max_n = std::min<std::size_t>(
            n, model.n_coefficients.size() - 1);
        for (std::size_t lag = 0; lag <= max_n; ++lag) {
          const auto degree = n - static_cast<std::uint32_t>(lag);
          for (std::uint32_t column = 0;
               column < solution.dimension; ++column)
            residual -= model.n_coefficients[lag][
                static_cast<std::size_t>(row) * solution.dimension +
                column] *
                local_detail::to_ball(sector.coefficients[
                    local_detail::sector_index(
                        solution, epsilon, degree, column)]);
        }
        if (!residual.contains_zero())
          throw std::invalid_argument(
              "restored tail q/N payload does not enclose the retained Taylor recurrence at n=" +
              std::to_string(n) + ", row=" + std::to_string(row) +
              ", epsilon_index=" + std::to_string(epsilon));
      }
}

inline Magnitude geometric_tail_factor(const Magnitude& ratio_upper,
                                       ulong first_power,
                                       Magnitude* gap_lower_out = nullptr) {
  const auto gap_lower = Magnitude::positive_difference_lower(
      Magnitude::one(), ratio_upper);
  if (gap_lower.is_zero()) return Magnitude::zero();
  if (gap_lower_out != nullptr) *gap_lower_out = gap_lower;
  return ratio_upper.power_upper(first_power) /
         gap_lower;
}

inline std::string certificate_provenance(
    const RegularTaylorTailModel& model, const std::string& radius,
    const std::string& target) {
  return "certified full-local " + target +
      " tail; regular homogeneous finite-width ODE; q nonvanishing by "
      "triangle lower bound on |t|<=" + radius +
      "; Gronwall circle bound plus Cauchy coefficients; "
      "operator_identity=" + model.operator_identity +
      "; local_checkpoint_identity=" + model.local_checkpoint_identity +
      "; analytic prescriptions retained (ordinary Taylor tail is "
      "rim-independent in absolute value)";
}

}  // namespace tail_majorant_detail

namespace tail_majorant_detail {

inline std::vector<ComplexBall> finite_causal_multiplier_matrix(
    const physical_ode_detail::PreparedCausalEpsilonMultiplier& multiplier,
    EpsilonWindow epsilon) {
  const auto width = epsilon.width();
  if (width > std::numeric_limits<std::size_t>::max() / width)
    throw std::overflow_error(
        "finite causal multiplier matrix size overflows");
  std::vector<ComplexBall> matrix(width * width, ComplexBall(0));
  for (std::size_t column = 0; column < width; ++column) {
    auto source = physical_ode_detail::zero_epsilon_vector(epsilon, 1);
    source.coefficients[column] = ComplexBall(1);
    auto target = physical_ode_detail::zero_epsilon_vector(epsilon, 1);
    physical_ode_detail::accumulate_causal_multiplier(
        multiplier, source, 0, target, 0, ComplexBall(1));
    for (std::size_t row = 0; row < width; ++row)
      matrix[row * width + column] = target.coefficients[row];
  }
  return matrix;
}

inline Magnitude finite_matrix_infinity_norm_upper(
    const std::vector<ComplexBall>& matrix, std::size_t dimension) {
  if (dimension == 0 ||
      dimension > std::numeric_limits<std::size_t>::max() / dimension ||
      matrix.size() != dimension * dimension)
    throw std::invalid_argument(
        "finite matrix infinity norm received an invalid square matrix");
  auto maximum = Magnitude::zero();
  for (std::size_t row = 0; row < dimension; ++row) {
    auto sum = Magnitude::zero();
    for (std::size_t column = 0; column < dimension; ++column)
      sum += Magnitude::upper_abs(matrix[row * dimension + column]);
    maximum = Magnitude::maximum(maximum, sum);
  }
  return maximum;
}

inline std::vector<Magnitude> finite_causal_matrix_prefix_norm_upper(
    const std::vector<ComplexBall>& matrix, std::size_t width) {
  if (width == 0 ||
      width > std::numeric_limits<std::size_t>::max() / width ||
      matrix.size() != width * width)
    throw std::invalid_argument(
        "finite causal prefix norm received an invalid square matrix");
  std::vector<Magnitude> prefix;
  prefix.reserve(width);
  auto maximum = Magnitude::zero();
  for (std::size_t row = 0; row < width; ++row) {
    auto sum = Magnitude::zero();
    // Causality makes the finite epsilon matrix lower triangular.  Summing
    // only through the current row therefore gives the exact induced norm of
    // the invariant prefix, without allowing later private rows to enter it.
    for (std::size_t column = 0; column <= row; ++column)
      sum += Magnitude::upper_abs(matrix[row * width + column]);
    maximum = Magnitude::maximum(maximum, sum);
    prefix.push_back(maximum);
  }
  return prefix;
}

inline Magnitude finite_causal_multiplier_norm_upper(
    const physical_ode_detail::PreparedCausalEpsilonMultiplier& multiplier,
    EpsilonWindow epsilon) {
  const auto matrix = finite_causal_multiplier_matrix(multiplier, epsilon);
  return finite_matrix_infinity_norm_upper(matrix, epsilon.width());
}

inline std::vector<Magnitude> finite_causal_multiplier_prefix_norm_upper(
    const physical_ode_detail::PreparedCausalEpsilonMultiplier& multiplier,
    EpsilonWindow epsilon) {
  const auto matrix = finite_causal_multiplier_matrix(multiplier, epsilon);
  return finite_causal_matrix_prefix_norm_upper(
      matrix, epsilon.width());
}

inline Magnitude finite_causal_unit_inverse_norm_upper(
    const physical_ode_detail::PreparedCausalEpsilonMultiplier& q0,
    EpsilonWindow epsilon) {
  const auto width = epsilon.width();
  if (width > std::numeric_limits<std::size_t>::max() / width)
    throw std::overflow_error(
        "finite causal inverse matrix size overflows");
  std::vector<ComplexBall> inverse(width * width, ComplexBall(0));
  for (std::size_t column = 0; column < width; ++column) {
    auto rhs = physical_ode_detail::zero_epsilon_vector(epsilon, 1);
    rhs.coefficients[column] = ComplexBall(1);
    const auto solved =
        physical_ode_detail::solve_formal_unit_q0(q0, rhs, 1);
    for (std::size_t row = 0; row < width; ++row)
      inverse[row * width + column] = solved.coefficients[row];
  }
  return finite_matrix_infinity_norm_upper(inverse, width);
}

inline std::vector<Magnitude>
finite_causal_unit_inverse_prefix_norm_upper(
    const physical_ode_detail::PreparedCausalEpsilonMultiplier& q0,
    EpsilonWindow epsilon) {
  const auto width = epsilon.width();
  if (width > std::numeric_limits<std::size_t>::max() / width)
    throw std::overflow_error(
        "finite causal inverse prefix matrix size overflows");
  std::vector<ComplexBall> inverse(width * width, ComplexBall(0));
  for (std::size_t column = 0; column < width; ++column) {
    auto rhs = physical_ode_detail::zero_epsilon_vector(epsilon, 1);
    rhs.coefficients[column] = ComplexBall(1);
    const auto solved =
        physical_ode_detail::solve_formal_unit_q0(q0, rhs, 1);
    for (std::size_t row = 0; row < width; ++row)
      inverse[row * width + column] = solved.coefficients[row];
  }
  return finite_causal_matrix_prefix_norm_upper(inverse, width);
}

struct FiniteCausalQDiskInverseCertificate {
  TailMajorantStatus status = TailMajorantStatus::Inconclusive;
  std::vector<Magnitude> prefix_norm_upper;
  Magnitude diagonal_lower = Magnitude::zero();
  std::size_t accepted_tiles = 0;
  std::size_t subdivided_tiles = 0;
  std::string detail;
};

struct ExactComplexDiskTile {
  Rational real_center;
  Rational imag_center;
  Rational half_width;
  std::uint32_t depth = 0;
};

inline bool tile_is_strictly_outside_disk(
    const ExactComplexDiskTile& tile, const Rational& radius) {
  auto real_distance =
      abs_rational(tile.real_center) - tile.half_width;
  auto imag_distance =
      abs_rational(tile.imag_center) - tile.half_width;
  if (real_distance.sign() < 0) real_distance = Rational(0);
  if (imag_distance.sign() < 0) imag_distance = Rational(0);
  return radius * radius <
      real_distance * real_distance + imag_distance * imag_distance;
}

inline ComplexBall exact_complex_tile_ball(
    const ExactComplexDiskTile& tile) {
  auto interval = ComplexBall::from_strings(
      tile.real_center.str(), tile.imag_center.str());
  const auto error =
      ComplexBall::from_strings(tile.half_width.str());
  arb_add_error(
      acb_realref(interval.raw()), acb_realref(error.raw()));
  arb_add_error(
      acb_imagref(interval.raw()), acb_realref(error.raw()));
  return interval;
}

inline std::vector<ComplexBall>
evaluate_causal_q_on_tile(
    const std::vector<std::vector<ComplexBall>>& coefficients,
    const ComplexBall& tile, std::size_t width) {
  if (coefficients.empty() ||
      std::any_of(coefficients.begin(), coefficients.end(),
                  [width](const auto& lag) {
                    return lag.size() != width;
                  }))
    throw std::invalid_argument(
        "causal Q disk evaluation received incomplete coefficients");
  std::vector<ComplexBall> value(width, ComplexBall(0));
  for (std::size_t reverse = coefficients.size();
       reverse-- > 0;) {
    for (std::size_t epsilon_lag = 0;
         epsilon_lag < width; ++epsilon_lag)
      value[epsilon_lag] =
          value[epsilon_lag] * tile +
          coefficients[reverse][epsilon_lag];
  }
  return value;
}

inline std::vector<Magnitude>
finite_causal_series_inverse_prefix_norm_upper(
    const std::vector<ComplexBall>& coefficients) {
  if (coefficients.empty() || coefficients.front().contains_zero())
    throw std::domain_error(
        "finite causal series inverse requires a separated diagonal");
  const auto width = coefficients.size();
  std::vector<ComplexBall> inverse(width, ComplexBall(0));
  inverse.front() = ComplexBall(1) / coefficients.front();
  for (std::size_t epsilon_lag = 1;
       epsilon_lag < width; ++epsilon_lag) {
    ComplexBall convolution(0);
    for (std::size_t degree = 1;
         degree <= epsilon_lag; ++degree)
      convolution += coefficients[degree] *
          inverse[epsilon_lag - degree];
    inverse[epsilon_lag] =
        -convolution / coefficients.front();
  }
  std::vector<Magnitude> prefix;
  prefix.reserve(width);
  auto row_sum = Magnitude::zero();
  for (const auto& coefficient : inverse) {
    row_sum += Magnitude::upper_abs(coefficient);
    prefix.push_back(row_sum);
  }
  return prefix;
}

inline std::vector<ComplexBall>
finite_causal_series_inverse_coefficients(
    const std::vector<ComplexBall>& coefficients) {
  if (coefficients.empty() || coefficients.front().contains_zero())
    throw std::domain_error(
        "finite causal series inverse requires a separated diagonal");
  const auto width = coefficients.size();
  std::vector<ComplexBall> inverse(width, ComplexBall(0));
  inverse.front() = ComplexBall(1) / coefficients.front();
  for (std::size_t epsilon_lag = 1;
       epsilon_lag < width; ++epsilon_lag) {
    ComplexBall convolution(0);
    for (std::size_t degree = 1;
         degree <= epsilon_lag; ++degree)
      convolution += coefficients[degree] *
          inverse[epsilon_lag - degree];
    inverse[epsilon_lag] =
        -convolution / coefficients.front();
  }
  return inverse;
}

struct PhysicalNormalizedODEDiskBounds {
  TailMajorantStatus status = TailMajorantStatus::Inconclusive;
  // [epsilon_lag][physical_row] bounds the row sum of the corresponding
  // causal coefficient of Q^-1 C/t over the complete disk.
  std::vector<std::vector<Magnitude>> degree_row_upper;
  Magnitude q_diagonal_lower = Magnitude::zero();
  std::size_t accepted_tiles = 0;
  std::size_t subdivided_tiles = 0;
  std::string detail;
};

inline std::vector<ComplexBall>
evaluate_causal_matrix_polynomial_on_tile(
    const std::vector<std::vector<ComplexBall>>& coefficients,
    const ComplexBall& tile, std::size_t width,
    std::size_t matrix_size, std::size_t first_lag) {
  if (coefficients.size() <= first_lag ||
      std::any_of(coefficients.begin(), coefficients.end(),
                  [width, matrix_size](const auto& lag) {
                    return lag.size() != width * matrix_size;
                  }))
    throw std::invalid_argument(
        "causal matrix disk evaluation received incomplete coefficients");
  std::vector<ComplexBall> value(
      width * matrix_size, ComplexBall(0));
  for (std::size_t reverse = coefficients.size();
       reverse-- > first_lag;) {
    for (std::size_t entry = 0;
         entry < value.size(); ++entry)
      value[entry] =
          value[entry] * tile + coefficients[reverse][entry];
  }
  return value;
}

inline PhysicalNormalizedODEDiskBounds
certify_physical_normalized_ode_disk_bounds(
    const std::vector<std::vector<ComplexBall>>& q_coefficients,
    const std::vector<std::vector<ComplexBall>>& c_coefficients,
    std::uint32_t dimension, const Rational& radius,
    std::uint32_t minimum_subdivision_depth = 3) {
  PhysicalNormalizedODEDiskBounds result;
  if (dimension == 0 || radius.sign() <= 0 ||
      q_coefficients.empty() || q_coefficients.front().empty() ||
      c_coefficients.empty()) {
    result.detail =
        "normalized physical ODE disk bound requires a positive radius, "
        "dimension, Q, and a structural C lag";
    return result;
  }
  const auto width = q_coefficients.front().size();
  const auto matrix_size =
      static_cast<std::size_t>(dimension) * dimension;
  if (std::any_of(q_coefficients.begin(), q_coefficients.end(),
                  [width](const auto& lag) {
                    return lag.size() != width;
                  }) ||
      std::any_of(c_coefficients.begin(), c_coefficients.end(),
                  [width, matrix_size](const auto& lag) {
                    return lag.size() != width * matrix_size;
                  })) {
    result.detail =
        "normalized physical ODE disk bound has inconsistent coefficient "
        "shape";
    return result;
  }
  if (std::any_of(c_coefficients.front().begin(),
                  c_coefficients.front().end(),
                  [](const ComplexBall& value) {
                    return !value.is_zero();
                  })) {
    result.detail =
        "normalized physical ODE disk bound requires structurally exact "
        "C(0)=0 before forming C/t";
    return result;
  }

  constexpr std::uint32_t kMaximumDepth = 18;
  constexpr std::size_t kMaximumVisitedTiles = 262144;
  std::vector<ExactComplexDiskTile> pending;
  pending.push_back(
      {Rational(0), Rational(0), radius, 0});
  result.degree_row_upper.assign(
      width, std::vector<Magnitude>(
                 dimension, Magnitude::zero()));
  std::optional<Magnitude> minimum_diagonal;
  std::size_t visited_tiles = 0;
  while (!pending.empty()) {
    if (++visited_tiles > kMaximumVisitedTiles) {
      result.detail =
          "normalized physical ODE disk cover exceeded its tile budget";
      return result;
    }
    auto tile = std::move(pending.back());
    pending.pop_back();
    if (tile_is_strictly_outside_disk(tile, radius))
      continue;
    const auto tile_ball = exact_complex_tile_ball(tile);
    const auto q =
        evaluate_causal_q_on_tile(q_coefficients, tile_ball, width);
    const auto diagonal_lower = Magnitude::lower_abs(q.front());
    if (tile.depth < minimum_subdivision_depth ||
        diagonal_lower.is_zero()) {
      if (tile.depth >= kMaximumDepth) {
        result.detail =
            "normalized physical ODE disk cover cannot separate the scalar "
            "Q diagonal from zero";
        return result;
      }
      const auto child_half = tile.half_width / Rational(2);
      for (const int real_sign : {-1, 1})
        for (const int imag_sign : {-1, 1})
          pending.push_back({
              tile.real_center +
                  Rational(real_sign) * child_half,
              tile.imag_center +
                  Rational(imag_sign) * child_half,
              child_half, tile.depth + 1});
      ++result.subdivided_tiles;
      continue;
    }

    const auto q_inverse =
        finite_causal_series_inverse_coefficients(q);
    const auto c_over_t =
        c_coefficients.size() == 1
            ? std::vector<ComplexBall>(
                  width * matrix_size, ComplexBall(0))
            : evaluate_causal_matrix_polynomial_on_tile(
                  c_coefficients, tile_ball, width,
                  matrix_size, 1);
    std::vector<ComplexBall> normalized(
        width * matrix_size, ComplexBall(0));
    for (std::size_t epsilon_degree = 0;
         epsilon_degree < width; ++epsilon_degree)
      for (std::size_t input_degree = 0;
           input_degree <= epsilon_degree; ++input_degree) {
        const auto& scalar = q_inverse[input_degree];
        for (std::size_t entry = 0;
             entry < matrix_size; ++entry)
          normalized[epsilon_degree * matrix_size + entry] +=
              scalar *
              c_over_t[(epsilon_degree - input_degree) *
                           matrix_size +
                       entry];
      }
    for (std::size_t epsilon_degree = 0;
         epsilon_degree < width; ++epsilon_degree)
      for (std::uint32_t row = 0; row < dimension; ++row) {
        auto row_sum = Magnitude::zero();
        for (std::uint32_t column = 0;
             column < dimension; ++column)
          row_sum += Magnitude::upper_abs(
              normalized[epsilon_degree * matrix_size +
                         static_cast<std::size_t>(row) * dimension +
                         column]);
        result.degree_row_upper[epsilon_degree][row] =
            Magnitude::maximum(
                result.degree_row_upper[epsilon_degree][row],
                row_sum);
      }
    if (!minimum_diagonal.has_value() ||
        diagonal_lower <= *minimum_diagonal)
      minimum_diagonal = diagonal_lower;
    ++result.accepted_tiles;
  }
  if (result.accepted_tiles == 0 ||
      !minimum_diagonal.has_value()) {
    result.detail =
        "normalized physical ODE disk cover produced no accepted tiles";
    return result;
  }
  result.status = TailMajorantStatus::Certified;
  result.q_diagonal_lower = *minimum_diagonal;
  result.detail =
      "Q^-1 C/t was interval-evaluated before norms on an adaptive "
      "complex-disk cover";
  return result;
}

inline std::vector<Magnitude>
weighted_physical_ode_prefix_norm_upper(
    const PhysicalNormalizedODEDiskBounds& bounds,
    ulong epsilon_weight_base) {
  if (epsilon_weight_base == 0 ||
      (epsilon_weight_base & (epsilon_weight_base - 1UL)) != 0 ||
      bounds.status != TailMajorantStatus::Certified ||
      bounds.degree_row_upper.empty())
    throw std::invalid_argument(
        "weighted physical ODE prefix norm requires certified bounds and "
        "an exact power-of-two weight");
  const auto dimension = bounds.degree_row_upper.front().size();
  if (dimension == 0 ||
      std::any_of(bounds.degree_row_upper.begin(),
                  bounds.degree_row_upper.end(),
                  [dimension](const auto& rows) {
                    return rows.size() != dimension;
                  }))
    throw std::invalid_argument(
        "weighted physical ODE prefix norm has inconsistent dimension");
  const auto epsilon_weight =
      Magnitude::from_ui(epsilon_weight_base);
  std::vector<Magnitude> prefix;
  prefix.reserve(bounds.degree_row_upper.size());
  auto maximum = Magnitude::zero();
  std::vector<Magnitude> row_sums(
      dimension, Magnitude::zero());
  for (std::size_t epsilon_degree = 0;
       epsilon_degree < bounds.degree_row_upper.size();
       ++epsilon_degree) {
    const auto scale = epsilon_weight.power_upper(
        static_cast<ulong>(epsilon_degree));
    for (std::size_t row = 0; row < dimension; ++row) {
      row_sums[row] +=
          bounds.degree_row_upper[epsilon_degree][row] / scale;
      maximum = Magnitude::maximum(maximum, row_sums[row]);
    }
    prefix.push_back(maximum);
  }
  return prefix;
}

inline FiniteCausalQDiskInverseCertificate
certify_finite_causal_q_disk_inverse(
    const std::vector<std::vector<ComplexBall>>& coefficients,
    const Rational& radius) {
  FiniteCausalQDiskInverseCertificate result;
  if (radius.sign() <= 0 || coefficients.empty() ||
      coefficients.front().empty()) {
    result.detail =
        "finite causal Q disk inverse requires positive radius and "
        "nonempty coefficients";
    return result;
  }
  const auto width = coefficients.front().size();
  if (std::any_of(coefficients.begin(), coefficients.end(),
                  [width](const auto& lag) {
                    return lag.size() != width;
                  })) {
    result.detail =
        "finite causal Q disk inverse has inconsistent epsilon width";
    return result;
  }

  // A square cover is deliberately refined only where interval dependency
  // prevents separation of the scalar diagonal q_0(t).  Squares wholly
  // outside the closed disk are discarded exactly.  Thus an exterior root
  // near the Cauchy circle cannot poison the proof merely because it lies in
  // the initial bounding square.
  constexpr std::uint32_t kMaximumDepth = 18;
  constexpr std::size_t kMaximumVisitedTiles = 262144;
  std::vector<ExactComplexDiskTile> pending;
  pending.push_back(
      {Rational(0), Rational(0), radius, 0});
  result.prefix_norm_upper.assign(width, Magnitude::zero());
  std::optional<Magnitude> minimum_diagonal;
  std::size_t visited_tiles = 0;
  while (!pending.empty()) {
    if (++visited_tiles > kMaximumVisitedTiles) {
      result.detail =
          "adaptive finite-causal Q disk cover exceeded its tile budget";
      return result;
    }
    auto tile = std::move(pending.back());
    pending.pop_back();
    if (tile_is_strictly_outside_disk(tile, radius))
      continue;

    const auto tile_ball = exact_complex_tile_ball(tile);
    const auto q =
        evaluate_causal_q_on_tile(coefficients, tile_ball, width);
    const auto diagonal_lower = Magnitude::lower_abs(q.front());
    if (diagonal_lower.is_zero()) {
      if (tile.depth >= kMaximumDepth) {
        result.detail =
            "adaptive finite-causal Q disk cover cannot separate the "
            "scalar diagonal from zero";
        return result;
      }
      const auto child_half = tile.half_width / Rational(2);
      for (const int real_sign : {-1, 1})
        for (const int imag_sign : {-1, 1})
          pending.push_back({
              tile.real_center +
                  Rational(real_sign) * child_half,
              tile.imag_center +
                  Rational(imag_sign) * child_half,
              child_half, tile.depth + 1});
      ++result.subdivided_tiles;
      continue;
    }

    const auto inverse_prefix =
        finite_causal_series_inverse_prefix_norm_upper(q);
    for (std::size_t epsilon_index = 0;
         epsilon_index < width; ++epsilon_index)
      result.prefix_norm_upper[epsilon_index] =
          Magnitude::maximum(
              result.prefix_norm_upper[epsilon_index],
              inverse_prefix[epsilon_index]);
    if (!minimum_diagonal.has_value() ||
        diagonal_lower <= *minimum_diagonal)
      minimum_diagonal = diagonal_lower;
    ++result.accepted_tiles;
  }
  if (result.accepted_tiles == 0 ||
      !minimum_diagonal.has_value() ||
      std::any_of(result.prefix_norm_upper.begin(),
                  result.prefix_norm_upper.end(),
                  [](const Magnitude& value) {
                    return value.is_zero() || !value.is_finite();
                  })) {
    result.detail =
        "adaptive finite-causal Q disk cover produced no finite inverse "
        "bound";
    return result;
  }
  result.status = TailMajorantStatus::Certified;
  result.diagonal_lower = *minimum_diagonal;
  result.detail =
      "scalar causal Q is separated from zero on an adaptive complex-disk "
      "cover and every finite epsilon prefix is bounded by triangular "
      "back-substitution";
  return result;
}

inline Magnitude finite_physical_c_lag_norm_upper(
    const std::vector<PhysicalODEMatrixEntry<ComplexBall>>& entries,
    EpsilonWindow epsilon, std::uint32_t dimension) {
  const auto width = epsilon.width();
  std::vector<Magnitude> row_sums(
      width * static_cast<std::size_t>(dimension), Magnitude::zero());
  for (const auto& entry : entries) {
    const auto multiplier =
        physical_ode_detail::prepare_causal_multiplier(entry.value);
    const auto matrix = finite_causal_multiplier_matrix(
        multiplier, epsilon);
    for (std::size_t target_epsilon = 0;
         target_epsilon < width; ++target_epsilon)
      for (std::size_t source_epsilon = 0;
           source_epsilon < width; ++source_epsilon)
        row_sums[target_epsilon * dimension + entry.row] +=
            Magnitude::upper_abs(
                matrix[target_epsilon * width + source_epsilon]);
  }
  auto maximum = Magnitude::zero();
  for (const auto& sum : row_sums)
    maximum = Magnitude::maximum(maximum, sum);
  return maximum;
}

inline std::vector<Magnitude> finite_physical_c_lag_prefix_norm_upper(
    const std::vector<PhysicalODEMatrixEntry<ComplexBall>>& entries,
    EpsilonWindow epsilon, std::uint32_t dimension) {
  const auto width = epsilon.width();
  std::vector<Magnitude> row_sums(
      width * static_cast<std::size_t>(dimension), Magnitude::zero());
  for (const auto& entry : entries) {
    const auto multiplier =
        physical_ode_detail::prepare_causal_multiplier(entry.value);
    const auto matrix = finite_causal_multiplier_matrix(
        multiplier, epsilon);
    for (std::size_t target_epsilon = 0;
         target_epsilon < width; ++target_epsilon)
      for (std::size_t source_epsilon = 0;
           source_epsilon <= target_epsilon; ++source_epsilon)
        row_sums[target_epsilon * dimension + entry.row] +=
            Magnitude::upper_abs(
                matrix[target_epsilon * width + source_epsilon]);
  }
  std::vector<Magnitude> prefix;
  prefix.reserve(width);
  auto maximum = Magnitude::zero();
  for (std::size_t target_epsilon = 0;
       target_epsilon < width; ++target_epsilon) {
    for (std::uint32_t row = 0; row < dimension; ++row)
      maximum = Magnitude::maximum(
          maximum,
          row_sums[target_epsilon * dimension + row]);
    prefix.push_back(maximum);
  }
  return prefix;
}

}  // namespace tail_majorant_detail

template <typename Scalar>
RegularTaylorTailModelResult prepare_regular_homogeneous_tail_model(
    const PreparedRecurrenceOperator<Scalar>& prepared,
    const RecurrenceProblem<Scalar>& problem,
    const LocalSolution<Scalar>& solution,
    std::string operator_identity) {
  static_assert(std::is_same_v<Scalar, Rational> ||
                    std::is_same_v<Scalar, ComplexBall>,
                "certified regular tail models support Rational or Acb "
                "recurrences only");
  using namespace tail_majorant_detail;

  validate_local_solution(solution, false);
  if (operator_identity.empty() || solution.checkpoint_identity.empty())
    throw std::invalid_argument(
        "tail-majorant provenance identities must be nonempty");
  if (prepared.dimension != problem.dimension ||
      prepared.dimension != solution.dimension ||
      prepared.frame_base != problem.frame_base ||
      prepared.frame_width != problem.frame_width ||
      problem.nmax != solution.taylor_complete_max)
    throw std::invalid_argument(
        "tail-majorant prepared/run/local dimensions or windows disagree");
  if (!solution.error.empty())
    return unsupported_model(
        "regular tail model cannot absorb an existing local error envelope");
  if (!problem.has_initial)
    return unsupported_model(
        "regular tail model requires a retained homogeneous initial value");
  if (problem.source.has_value())
    return unsupported_model(
        "sourced recurrences need a certified unseen source tail before a "
        "full local majorant can be proved");
  if (problem.log_max != 0 || solution.sectors.size() != 1 ||
      solution.sectors.front().log_power != 0)
    return unsupported_model(
        "logarithmic local sectors are outside the ordinary regular-tail "
        "certificate");
  const auto& sector = solution.sectors.front();
  if (sector.a.is_zero != TruthValue::Yes ||
      sector.b.is_zero != TruthValue::Yes ||
      !exact_scalar_zero(problem.a_target) ||
      !exact_scalar_zero(problem.b_target))
    return unsupported_model(
        "fractional or regulator-dependent local powers require a singular "
        "sector majorant; they are never treated as ordinary Taylor tails");
  if (!prepared.assembly_matrix.has_value() ||
      !prepared.assembly_matrix->identity)
    return unsupported_model(
        "regular tail model currently requires exact identity assembly so "
        "the bounded ODE and retained coefficients share one basis");
  if (prepared.d_lags.empty() || prepared.nhat_lags.empty())
    throw std::invalid_argument(
        "tail-majorant prepared operator has no recurrence lags");
  std::vector<std::uint8_t> block_seen(prepared.dimension, 0);
  bool singleton_partition = prepared.blocks.size() == prepared.dimension;
  for (const auto& block : prepared.blocks) {
    if (block.columns.size() != 1 ||
        block.columns.front() >= prepared.dimension ||
        block_seen[block.columns.front()]) {
      singleton_partition = false;
      break;
    }
    block_seen[block.columns.front()] = 1;
  }
  if (!singleton_partition)
    return unsupported_model(
        "ordinary regular tail model requires the exact zero indicial "
        "matrix to be represented by singleton Jordan blocks");

  RegularTaylorTailModel model;
  model.epsilon = solution.epsilon;
  model.dimension = solution.dimension;
  model.taylor_complete_max = solution.taylor_complete_max;
  model.chart = solution.chart;
  model.prescriptions = solution.prescriptions;
  model.operator_identity = std::move(operator_identity);
  model.local_checkpoint_identity = solution.checkpoint_identity;

  std::vector<Scalar> exact_q_coefficients;
  exact_q_coefficients.reserve(prepared.d_lags.size());
  model.q_coefficients.reserve(prepared.d_lags.size());
  for (const auto& lag : prepared.d_lags) {
    auto coefficient = ScalarTraits<Scalar>::zero();
    for (const auto& shifted : lag) {
      if (shifted.shift != 0)
        return unsupported_model(
            "regular tail majorant requires an exactly epsilon-decoupled "
            "scalar denominator (every retained q shift must be zero)");
      coefficient += shifted.value;
    }
    exact_q_coefficients.push_back(coefficient);
    model.q_coefficients.push_back(local_detail::to_ball(coefficient));
  }
  if (std::any_of(model.q_coefficients.begin(),
                  model.q_coefficients.end(),
                  [](const ComplexBall& value) {
                    return !value.is_finite();
                  }))
    return inconclusive_model(
        "nonfinite q coefficient enclosure cannot support a tail proof");
  if (model.q_coefficients.front().is_zero())
    return unsupported_model(
        "regular tail equation has structurally zero q(0)");
  if (divisor_contains_zero(exact_q_coefficients.front()))
    return inconclusive_model(
        "q(0) enclosure contains zero, so the retained coefficient "
        "recurrence cannot be certified at this precision");
  if (!prepared.d0_inverse_scalar.has_value())
    return unsupported_model(
        "regular tail slice requires a retained scalar d0 inverse which "
        "certifiably inverts q(0)");
  const auto recomputed_d0_inverse =
      ScalarTraits<Scalar>::one() / exact_q_coefficients.front();
  if (!encloses_recurrence_value(
          *prepared.d0_inverse_scalar, recomputed_d0_inverse))
    return inconclusive_model(
        "retained scalar d0 inverse does not enclose a direct certified "
        "inverse of q(0)");

  std::vector<std::vector<Scalar>> exact_nhat_coefficients;
  exact_nhat_coefficients.reserve(prepared.nhat_lags.size());
  model.n_coefficients.reserve(prepared.nhat_lags.size());
  model.n_row_sum_upper.reserve(prepared.nhat_lags.size());
  for (std::size_t lag_index = 0;
       lag_index < prepared.nhat_lags.size(); ++lag_index) {
    bool supported = true;
    std::string detail;
    auto matrix = aggregate_nhat_lag(
        prepared.nhat_lags[lag_index], prepared.dimension,
        &supported, &detail);
    if (!supported) return unsupported_model(std::move(detail));
    if (lag_index == 0 &&
        std::any_of(matrix.begin(), matrix.end(),
                    [](const Scalar& value) {
                      return !exact_scalar_zero(value);
                    }))
      return unsupported_model(
          "N(0) is not the exact zero matrix, so theta f=N f is a "
          "regular-singular rather than ordinary regular equation");
    std::vector<ComplexBall> enclosed_matrix;
    enclosed_matrix.reserve(matrix.size());
    for (const auto& entry : matrix)
      enclosed_matrix.push_back(local_detail::to_ball(entry));
    model.n_coefficients.push_back(std::move(enclosed_matrix));
    exact_nhat_coefficients.push_back(matrix);
    auto norm = matrix_infinity_norm_upper(matrix, prepared.dimension);
    if (!norm.is_finite())
      return inconclusive_model(
          "nonfinite N coefficient enclosure cannot support a tail proof");
    model.n_row_sum_upper.push_back(std::move(norm));
  }

  // Bind the theorem to the actual retained polynomial.  Starting from the
  // parsed initial frame, propagate the ordinary coefficient recurrence
  // forward using enclosures.  Every retained Acb coefficient must contain
  // that forward enclosure (Rational coefficients must agree exactly).
  // Therefore the stored polynomial encloses the same unique analytic
  // solution whose unseen tail is bounded by the certificate constructed
  // here; metadata alone is never a sufficient certificate.
  const auto expected_initial_size =
      static_cast<std::size_t>(prepared.dimension) * prepared.frame_width;
  if (problem.initial.size() != expected_initial_size)
    throw std::invalid_argument(
        "tail-majorant initial tensor has the wrong ordinary-sector size");
  const auto coefficient_at = [&](std::size_t epsilon,
                                  std::uint32_t n,
                                  std::uint32_t component)
      -> const Scalar& {
    return sector.coefficients[local_detail::sector_index(
        solution, epsilon, n, component)];
  };
  for (std::size_t epsilon = 0; epsilon < solution.epsilon.width();
       ++epsilon) {
    const auto power = static_cast<std::int64_t>(
        solution.epsilon.min_power) + static_cast<std::int64_t>(epsilon);
    const auto work_index = power - prepared.frame_base;
    if (work_index < 0 ||
        work_index >= static_cast<std::int64_t>(prepared.frame_width))
      throw std::invalid_argument(
          "tail-majorant local epsilon window lies outside the prepared "
          "work frame");
    for (std::uint32_t component = 0; component < prepared.dimension;
         ++component) {
      const auto& initial = problem.initial[
          static_cast<std::size_t>(component) * prepared.frame_width +
          static_cast<std::size_t>(work_index)];
      if (!encloses_recurrence_value(
              coefficient_at(epsilon, 0, component), initial))
        return inconclusive_model(
            "retained Taylor constant does not enclose its parsed initial "
            "recurrence value");
    }
    for (std::uint32_t n = 1; n <= solution.taylor_complete_max; ++n) {
      std::vector<Scalar> rhs(
          prepared.dimension, ScalarTraits<Scalar>::zero());
      const auto max_nhat = std::min<std::size_t>(
          n, exact_nhat_coefficients.size() - 1);
      for (std::size_t lag = 1; lag <= max_nhat; ++lag)
        for (std::uint32_t row = 0; row < prepared.dimension; ++row)
          for (std::uint32_t column = 0; column < prepared.dimension;
               ++column)
            rhs[row] += exact_nhat_coefficients[lag][
                static_cast<std::size_t>(row) * prepared.dimension +
                column] *
                coefficient_at(epsilon, n - static_cast<std::uint32_t>(lag),
                               column);
      const auto max_q = std::min<std::size_t>(
          n, exact_q_coefficients.size() - 1);
      for (std::size_t lag = 1; lag <= max_q; ++lag) {
        const auto degree = n - static_cast<std::uint32_t>(lag);
        const auto scale = exact_q_coefficients[lag] *
                           ScalarTraits<Scalar>::integer(degree);
        for (std::uint32_t row = 0; row < prepared.dimension; ++row)
          rhs[row] -= scale * coefficient_at(epsilon, degree, row);
      }
      for (std::uint32_t row = 0; row < prepared.dimension; ++row) {
        // Match the certified recurrence kernel's operation order: first the
        // affine n inverse, then the retained d0 inverse.  This avoids asking
        // a valid Acb coefficient to contain a differently rounded but
        // algebraically equivalent enclosure.
        const auto expected =
            (rhs[row] / ScalarTraits<Scalar>::integer(n)) *
            *prepared.d0_inverse_scalar;
        if (!encloses_recurrence_value(
                coefficient_at(epsilon, n, row), expected))
          return inconclusive_model(
              "retained Taylor tensor does not forward-enclose its exact "
              "regular recurrence through the claimed complete order "
              "(n=" + std::to_string(n) + ", row=" +
              std::to_string(row) + ", epsilon_index=" +
              std::to_string(epsilon) + ")");
      }
    }
  }

  model.initial_row_upper.assign(
      solution.epsilon.width(), Magnitude::zero());
  for (std::size_t epsilon = 0; epsilon < solution.epsilon.width();
       ++epsilon)
    for (std::uint32_t component = 0; component < solution.dimension;
         ++component)
      model.initial_row_upper[epsilon] = Magnitude::maximum(
          model.initial_row_upper[epsilon], Magnitude::upper_abs(
              local_detail::to_ball(sector.coefficients[
                  local_detail::sector_index(
                      solution, epsilon, 0, component)])));
  if (std::any_of(model.initial_row_upper.begin(),
                  model.initial_row_upper.end(),
                  [](const Magnitude& value) {
                    return !value.is_finite();
                  }))
    return inconclusive_model(
        "nonfinite initial-value enclosure cannot support a tail proof");

  model.provenance =
      "regular homogeneous finite-width equation captured from retained "
      "prepared operator; exact ordinary sector; exact epsilon-decoupled "
      "q/N shifts; identity assembly; no source; stored coefficients "
      "forward-enclose the recurrence; operator_identity=" +
      model.operator_identity + "; local_checkpoint_identity=" +
      model.local_checkpoint_identity;
  return {TailMajorantStatus::Certified, std::move(model),
          "regular homogeneous tail model prepared"};
}

inline PhysicalRegularTaylorTailModelResult
prepare_physical_regular_homogeneous_tail_model(
    const PreparedPhysicalClearedODE<ComplexBall>& equation,
    const LocalSolution<ComplexBall>& solution,
    std::optional<std::uint32_t> reconstruction_complete_max =
        std::nullopt) {
  using namespace tail_majorant_detail;
  const auto unsupported = [](std::string detail) {
    return PhysicalRegularTaylorTailModelResult{
        TailMajorantStatus::Unsupported, std::nullopt, std::move(detail)};
  };
  const auto inconclusive = [](std::string detail) {
    return PhysicalRegularTaylorTailModelResult{
        TailMajorantStatus::Inconclusive, std::nullopt, std::move(detail)};
  };

  physical_ode_detail::validate_ode(equation);
  validate_local_solution(solution, false);
  if (solution.dimension != equation.dimension ||
      solution.checkpoint_identity.empty() ||
      equation.payload_identity.empty())
    throw std::invalid_argument(
        "physical tail model lost its dimension or provenance binding");
  if (!solution.error.empty())
    return unsupported(
        "physical regular tail model cannot absorb a local error envelope");
  if (solution.sectors.size() != 1 ||
      solution.sectors.front().log_power != 0 ||
      solution.sectors.front().a.is_zero != TruthValue::Yes ||
      solution.sectors.front().b.is_zero != TruthValue::Yes)
    return unsupported(
        "physical regular tail model requires one ordinary a=b=0, log=0 sector");
  if (!equation.c_lags.front().empty())
    return unsupported(
        "physical regular tail model requires a structurally ordinary center C(0)=0");
  const auto reconstructed_complete_max =
      reconstruction_complete_max.value_or(
          solution.taylor_complete_max);
  if (reconstructed_complete_max < solution.taylor_complete_max)
    throw std::invalid_argument(
        "physical tail reconstruction cannot discard retained Taylor coefficients");

  EpsilonVector initial;
  initial.epsilon = solution.epsilon;
  initial.dimension = solution.dimension;
  initial.coefficients.assign(
      initial.epsilon.width() * initial.dimension, ComplexBall(0));
  const auto& source_sector = solution.sectors.front();
  for (std::int64_t raw_power = solution.epsilon.min_power;
       raw_power <= solution.epsilon.complete_max; ++raw_power) {
    const auto power = static_cast<std::int32_t>(raw_power);
    const auto epsilon_index = static_cast<std::size_t>(
        raw_power - solution.epsilon.min_power);
    for (std::uint32_t component = 0;
         component < solution.dimension; ++component)
      initial.at(power, component) =
          source_sector.coefficients[local_detail::sector_index(
              solution, epsilon_index, 0, component)];
  }
  const auto evolution = evolve_ordinary_center_value(
      equation, initial, reconstructed_complete_max);
  if (!evolution.eligible)
    return unsupported(
        "physical regular Taylor reconstruction is ineligible: " +
        evolution.reason);

  PhysicalRegularTaylorTailModel model;
  model.epsilon = solution.epsilon;
  model.dimension = solution.dimension;
  model.taylor_complete_max = reconstructed_complete_max;
  model.chart = solution.chart;
  model.prescriptions = solution.prescriptions;
  model.physical_payload_identity = equation.payload_identity;
  model.local_checkpoint_identity = solution.checkpoint_identity;
  model.reconstructed = solution;
  model.reconstructed.taylor_complete_max =
      reconstructed_complete_max;
  auto& reconstructed_sector = model.reconstructed.sectors.front();
  reconstructed_sector.coefficients.assign(
      model.reconstructed.sector_size(), ComplexBall(0));
  for (std::uint32_t taylor = 0;
       taylor <= reconstructed_complete_max; ++taylor)
    for (std::int64_t raw_power = solution.epsilon.min_power;
         raw_power <= solution.epsilon.complete_max; ++raw_power) {
      const auto power = static_cast<std::int32_t>(raw_power);
      const auto epsilon_index = static_cast<std::size_t>(
          raw_power - solution.epsilon.min_power);
      for (std::uint32_t component = 0;
           component < solution.dimension; ++component) {
        const auto coefficient = local_detail::sector_index(
            model.reconstructed, epsilon_index, taylor, component);
        const auto& expected = evolution.at(taylor).at(power, component);
        // The retained recurrence and this physical q/C replay are two
        // independently rounded rigorous evaluations of the same exact
        // coefficient.  Neither Arb ball is required to contain the other.
        // A nonempty intersection is the fail-closed consistency condition;
        // the freshly reconstructed ball remains the authoritative enclosure.
        if (taylor <= solution.taylor_complete_max) {
          const auto retained_coefficient =
              local_detail::sector_index(
                  solution, epsilon_index, taylor, component);
          if (!overlaps_recurrence_value(
                  source_sector.coefficients[retained_coefficient],
                  expected))
            return inconclusive(
                "retained Taylor tensor does not overlap its physical "
                "q/C recurrence through the claimed complete order "
                "(n=" + std::to_string(taylor) + ", component=" +
                std::to_string(component) + ", epsilon_index=" +
                std::to_string(epsilon_index) + ")");
        }
        reconstructed_sector.coefficients[coefficient] = expected;
      }
    }

  model.q_operator_prefix_norm_upper.reserve(equation.q_lags.size());
  model.q_causal_coefficients.reserve(equation.q_lags.size());
  for (const auto& coefficient : equation.q_lags) {
    if (coefficient.zero) {
      model.q_operator_prefix_norm_upper.push_back(
          std::vector<Magnitude>(
              solution.epsilon.width(), Magnitude::zero()));
      model.q_causal_coefficients.push_back(
          std::vector<ComplexBall>(
              solution.epsilon.width(), ComplexBall(0)));
      continue;
    }
    if (coefficient.valuation < 0)
      return unsupported(
          "physical tail model requires causal nonnegative q epsilon valuations");
    const auto multiplier =
        physical_ode_detail::prepare_causal_multiplier(coefficient);
    const auto matrix =
        finite_causal_multiplier_matrix(
            multiplier, solution.epsilon);
    model.q_operator_prefix_norm_upper.push_back(
        finite_causal_matrix_prefix_norm_upper(
            matrix, solution.epsilon.width()));
    std::vector<ComplexBall> causal_coefficients(
        solution.epsilon.width(), ComplexBall(0));
    for (std::size_t epsilon_lag = 0;
         epsilon_lag < solution.epsilon.width(); ++epsilon_lag)
      causal_coefficients[epsilon_lag] =
          matrix[epsilon_lag * solution.epsilon.width()];
    model.q_causal_coefficients.push_back(
        std::move(causal_coefficients));
  }
  const auto q0 = physical_ode_detail::prepare_causal_multiplier(
      equation.q_lags.front());
  if (q0.valuation != 0)
    return unsupported(
        "physical tail model requires q(0,eps) to be a formal epsilon unit");
  model.q0_inverse_prefix_norm_upper =
      finite_causal_unit_inverse_prefix_norm_upper(
          q0, solution.epsilon);

  model.c_operator_prefix_norm_upper.reserve(equation.c_lags.size());
  model.c_causal_coefficients.reserve(equation.c_lags.size());
  for (const auto& entries : equation.c_lags) {
    if (std::any_of(entries.begin(), entries.end(),
                    [](const auto& entry) {
                      return entry.value.valuation < 0;
                    }))
      return unsupported(
          "physical tail model requires causal nonnegative C epsilon valuations");
    model.c_operator_prefix_norm_upper.push_back(
        finite_physical_c_lag_prefix_norm_upper(
            entries, solution.epsilon, solution.dimension));
    const auto matrix_size =
        static_cast<std::size_t>(solution.dimension) *
        solution.dimension;
    std::vector<ComplexBall> causal_coefficients(
        solution.epsilon.width() * matrix_size,
        ComplexBall(0));
    for (const auto& entry : entries) {
      const auto multiplier =
          finite_causal_multiplier_matrix(
              physical_ode_detail::prepare_causal_multiplier(
                  entry.value),
              solution.epsilon);
      for (std::size_t epsilon_lag = 0;
           epsilon_lag < solution.epsilon.width(); ++epsilon_lag)
        causal_coefficients[
            epsilon_lag * matrix_size +
            static_cast<std::size_t>(entry.row) *
                solution.dimension +
            entry.column] +=
            multiplier[
                epsilon_lag * solution.epsilon.width()];
    }
    model.c_causal_coefficients.push_back(
        std::move(causal_coefficients));
  }
  const auto valid_prefix = [&](const auto& prefix) {
    return prefix.size() == solution.epsilon.width() &&
        std::all_of(prefix.begin(), prefix.end(),
                    [](const Magnitude& value) {
                      return value.is_finite();
                    });
  };
  if (model.q_operator_prefix_norm_upper.empty() ||
      model.q_causal_coefficients.empty() ||
      model.c_operator_prefix_norm_upper.empty() ||
      model.c_causal_coefficients.empty() ||
      !valid_prefix(model.q0_inverse_prefix_norm_upper) ||
      std::any_of(model.q0_inverse_prefix_norm_upper.begin(),
                  model.q0_inverse_prefix_norm_upper.end(),
                  [](const Magnitude& value) {
                    return value.is_zero();
                  }) ||
      !std::all_of(model.q_operator_prefix_norm_upper.begin(),
                   model.q_operator_prefix_norm_upper.end(),
                   valid_prefix) ||
      !std::all_of(model.q_causal_coefficients.begin(),
                   model.q_causal_coefficients.end(),
                   [width = solution.epsilon.width()](
                       const auto& coefficients) {
                     return coefficients.size() == width &&
                         std::all_of(
                             coefficients.begin(), coefficients.end(),
                             [](const ComplexBall& value) {
                               return value.is_finite();
                             });
                   }) ||
      !std::all_of(model.c_operator_prefix_norm_upper.begin(),
                   model.c_operator_prefix_norm_upper.end(),
                   valid_prefix) ||
      !std::all_of(
          model.c_causal_coefficients.begin(),
          model.c_causal_coefficients.end(),
          [expected = solution.epsilon.width() *
                      static_cast<std::size_t>(solution.dimension) *
                      solution.dimension](
              const auto& coefficients) {
            return coefficients.size() == expected &&
                std::all_of(
                    coefficients.begin(), coefficients.end(),
                    [](const ComplexBall& value) {
                      return value.is_finite();
                    });
          }))
    return inconclusive(
        "physical finite-epsilon prefix operator norm is zero, incomplete, or nonfinite");

  model.initial_row_upper.assign(
      solution.epsilon.width(), Magnitude::zero());
  for (std::size_t epsilon_index = 0;
       epsilon_index < solution.epsilon.width(); ++epsilon_index)
    for (std::uint32_t component = 0;
         component < solution.dimension; ++component)
      model.initial_row_upper[epsilon_index] = Magnitude::maximum(
          model.initial_row_upper[epsilon_index],
          Magnitude::upper_abs(initial.coefficients[
              epsilon_index * solution.dimension + component]));
  model.provenance =
      "transient finite-epsilon physical regular tail model; retained "
      "q(t,eps)/C(t,eps) payload reconstructed the Taylor prefix; the "
      "epsilon stack is bounded by its nested invariant causal prefixes; "
      "retained_taylor_complete_max=" +
      std::to_string(solution.taylor_complete_max) +
      "; reconstructed_taylor_complete_max=" +
      std::to_string(reconstructed_complete_max) + "; "
      "physical_payload_identity=" + model.physical_payload_identity +
      "; local_checkpoint_identity=" +
      model.local_checkpoint_identity;
  return {TailMajorantStatus::Certified, std::move(model),
          "physical finite-epsilon regular tail model prepared"};
}

// Certify the equation owned by an immutable regular-tail payload.  The
// payload is not caller data: it was captured from PreparedRecurrenceOperator
// and RecurrenceProblem, checked against every retained Taylor coefficient,
// and (on restore) replayed before it can be attached to a local again.
// Keeping the denominator cleared avoids introducing a ball division at a
// point where q may be small:
//
//                    q(t) theta f(t) - N(t) f(t) = 0.
//
// This slice is deliberately homogeneous.  A sourced local can only acquire
// an owner-bound residual payload after its source function and unseen source
// tail have equally strong retained provenance.
template <typename Scalar>
ResidualCertificate certify_regular_owner_bound_residual(
    const LocalSolution<Scalar>& solution,
    const RegularTaylorTailModel& model,
    const LocalEvaluation& evaluation,
    const RealEvaluationPoint& point,
    const Magnitude& relative_tolerance,
    ResidualScope scope = ResidualScope::StoredTruncation) {
  static_assert(std::is_same_v<Scalar, Rational> ||
                    std::is_same_v<Scalar, ComplexBall>,
                "owner-bound regular residuals support Rational or Acb "
                "locals only");
  using namespace tail_majorant_detail;

  validate_restored_regular_taylor_tail_model(
      model, solution, model.operator_identity);
  if (evaluation.value.dimension != model.dimension ||
      evaluation.theta_value.dimension != model.dimension ||
      evaluation.value.epsilon.min_power != model.epsilon.min_power ||
      evaluation.value.epsilon.complete_max != model.epsilon.complete_max ||
      evaluation.theta_value.epsilon.min_power != model.epsilon.min_power ||
      evaluation.theta_value.epsilon.complete_max !=
          model.epsilon.complete_max)
    throw std::invalid_argument(
        "owner-bound residual evaluation differs from its retained operator window");

  auto t = point.modulus;
  if (point.sign < 0) t = -t;

  auto q_at_point = ComplexBall(0);
  for (auto lag = model.q_coefficients.rbegin();
       lag != model.q_coefficients.rend(); ++lag)
    q_at_point = q_at_point * t + *lag;

  const auto matrix_size = static_cast<std::size_t>(model.dimension) *
                           model.dimension;
  std::vector<ComplexBall> n_at_point(matrix_size, ComplexBall(0));
  for (auto lag = model.n_coefficients.rbegin();
       lag != model.n_coefficients.rend(); ++lag)
    for (std::size_t entry = 0; entry < matrix_size; ++entry)
      n_at_point[entry] = n_at_point[entry] * t + (*lag)[entry];

  LocalEvaluation cleared = evaluation;
  for (auto& coefficient : cleared.theta_value.coefficients)
    coefficient *= q_at_point;
  EpsilonMatrix numerator;
  numerator.epsilon = {0, 0};
  numerator.dimension = model.dimension;
  numerator.coefficients = std::move(n_at_point);
  auto certificate = certify_theta_residual(
      cleared, numerator, std::nullopt, relative_tolerance, scope);
  if (scope == ResidualScope::StoredTruncation)
    certificate.detail +=
        "; equation derived from retained owner payload "
        "q(t) theta(f)=N(t) f (homogeneous zero source)";
  else
    certificate.detail +=
        "; retained q/N ownership is proven, but the requested full-local "
        "residual still requires propagation of certified value/theta tails";
  return certificate;
}

// Transfer the ordinary homogeneous theorem through one retained match
// materialization.  The Laurent weights are independent of t, so every
// complete epsilon row of their linear combination solves the same ODE.
// This is sound only when every receiving basis column has already been
// forward-certified against byte-for-byte equal prepared q/N enclosures.
// We also independently replay the finite convolution into `materialized`;
// merely presenting compatible provenance strings is never sufficient.
template <typename Scalar>
RegularTaylorTailModelResult
derive_materialized_regular_homogeneous_tail_model(
    const std::vector<const LocalSolution<Scalar>*>& basis,
    const std::vector<const RegularTaylorTailModelResult*>& basis_models,
    const FiniteLaurentVector<Scalar>& weights,
    const LocalSolution<Scalar>& materialized,
    const std::string& expected_operator_identity,
    const std::string& match_checkpoint_identity) {
  static_assert(std::is_same_v<Scalar, Rational> ||
                    std::is_same_v<Scalar, ComplexBall>,
                "matched regular tail models support Rational or Acb "
                "locals only");
  using namespace tail_majorant_detail;

  validate_local_solution(materialized, false);
  if (basis.empty() || basis.size() != weights.size() ||
      basis.size() != basis_models.size())
    throw std::invalid_argument(
        "matched tail derivation requires one model and weight per basis column");
  if (expected_operator_identity.empty() ||
      match_checkpoint_identity.empty() ||
      materialized.checkpoint_identity.empty())
    throw std::invalid_argument(
        "matched tail derivation identities must be nonempty");
  if (!materialized.error.empty())
    return unsupported_model(
        "matched regular tail model cannot absorb an existing local error envelope");
  if (materialized.sectors.size() != 1 ||
      materialized.sectors.front().log_power != 0 ||
      materialized.sectors.front().a.is_zero != TruthValue::Yes ||
      materialized.sectors.front().b.is_zero != TruthValue::Yes)
    return unsupported_model(
        "matched tail propagation requires one ordinary a=b=0, log=0 sector");

  const RegularTaylorTailModel* reference = nullptr;
  for (std::size_t column = 0; column < basis.size(); ++column) {
    if (basis[column] == nullptr || basis_models[column] == nullptr)
      throw std::invalid_argument(
          "matched tail derivation received a null retained basis object");
    validate_local_solution(*basis[column], false);
    const auto& candidate = *basis_models[column];
    if (candidate.status != TailMajorantStatus::Certified ||
        !candidate.model.has_value()) {
      const auto detail =
          "receiving basis column " + std::to_string(column) +
          " has no certified ordinary homogeneous tail model: " +
          candidate.detail;
      return candidate.status == TailMajorantStatus::Inconclusive
          ? inconclusive_model(detail) : unsupported_model(detail);
    }
    const auto& model = *candidate.model;
    try {
      require_model_binding(model, *basis[column]);
    } catch (const std::invalid_argument&) {
      return inconclusive_model(
          "receiving basis tail model is not immutably bound to column " +
          std::to_string(column));
    }
    if (model.operator_identity != expected_operator_identity)
      return unsupported_model(
          "receiving basis columns do not share the requested prepared operator identity");
    if (basis[column]->sectors.size() != 1 ||
        basis[column]->sectors.front().log_power != 0 ||
        basis[column]->sectors.front().a.is_zero != TruthValue::Yes ||
        basis[column]->sectors.front().b.is_zero != TruthValue::Yes)
      return unsupported_model(
          "receiving basis contains a singular, fractional, or logarithmic local sector");
    if (basis[column]->taylor_complete_max <
        materialized.taylor_complete_max)
      return inconclusive_model(
          "materialized Taylor order exceeds a certified receiving basis order");
    if (reference == nullptr) {
      reference = &model;
    } else if (!same_operator_payload(*reference, model)) {
      return unsupported_model(
          "receiving basis tail models do not share one identical prepared operator payload");
    }
  }
  if (reference == nullptr)
    throw std::logic_error("matched tail derivation lost its reference model");
  if (!same_chart_geometry(reference->chart, materialized.chart) ||
      !same_prescriptions(reference->prescriptions,
                          materialized.prescriptions) ||
      reference->dimension != materialized.dimension)
    return unsupported_model(
        "materialized local left the certified receiving chart/operator space");
  if (reference->q_coefficients.empty() ||
      reference->n_coefficients.empty())
    return inconclusive_model(
        "certified receiving operator payload is incomplete");
  const auto matrix_size = static_cast<std::size_t>(reference->dimension) *
                           reference->dimension;
  for (const auto& matrix : reference->n_coefficients)
    if (matrix.size() != matrix_size)
      return inconclusive_model(
          "certified receiving operator matrix payload is malformed");

  // Check the honest product window before looking at any coefficient.
  std::int64_t expected_min = std::numeric_limits<std::int64_t>::max();
  std::int64_t expected_complete =
      std::numeric_limits<std::int64_t>::max();
  for (std::size_t column = 0; column < basis.size(); ++column) {
    const auto product_min =
        static_cast<std::int64_t>(basis[column]->epsilon.min_power) +
        weights[column].min_power();
    const auto product_complete = std::min(
        static_cast<std::int64_t>(basis[column]->epsilon.complete_max) +
            weights[column].min_power(),
        static_cast<std::int64_t>(basis[column]->epsilon.min_power) +
            weights[column].complete_max());
    expected_min = std::min(expected_min, product_min);
    expected_complete = std::min(expected_complete, product_complete);
  }
  if (materialized.epsilon.min_power != expected_min ||
      materialized.epsilon.complete_max > expected_complete)
    return inconclusive_model(
        "materialized epsilon frame is not a lower-exact, upper-complete subset of its Laurent products");

  // Replay the exact operation grouping used by materialize_local_basis_weights:
  // sum each column's weight convolution first, then add the columns.  For
  // Acb this avoids mistaking harmless reassociation widening for a mismatch.
  const auto& output_sector = materialized.sectors.front();
  for (std::int64_t power = materialized.epsilon.min_power;
       power <= materialized.epsilon.complete_max; ++power) {
    const auto output_epsilon = static_cast<std::size_t>(
        power - materialized.epsilon.min_power);
    for (std::uint32_t n = 0; n <= materialized.taylor_complete_max; ++n)
      for (std::uint32_t component = 0;
           component < materialized.dimension; ++component) {
        auto expected = ScalarTraits<Scalar>::zero();
        for (std::size_t column = 0; column < basis.size(); ++column) {
          auto column_sum = ScalarTraits<Scalar>::zero();
          const auto& source = *basis[column];
          const auto& source_sector = source.sectors.front();
          for (std::int64_t weight_power = weights[column].min_power();
               weight_power <= weights[column].complete_max();
               ++weight_power) {
            const auto source_power = power - weight_power;
            if (source_power < source.epsilon.min_power ||
                source_power > source.epsilon.complete_max)
              continue;
            const auto scalar_weight = weights[column].coefficient(
                static_cast<std::int32_t>(weight_power));
            if (ScalarTraits<Scalar>::is_zero(scalar_weight)) continue;
            const auto source_epsilon = static_cast<std::size_t>(
                source_power - source.epsilon.min_power);
            column_sum += scalar_weight * source_sector.coefficients[
                local_detail::sector_index(
                    source, source_epsilon, n, component)];
          }
          expected += column_sum;
        }
        const auto& stored = output_sector.coefficients[
            local_detail::sector_index(
                materialized, output_epsilon, n, component)];
        if (!encloses_recurrence_value(stored, expected))
          return inconclusive_model(
              "materialized local does not enclose the retained basis/weight convolution at epsilon=" +
              std::to_string(power) + ", n=" + std::to_string(n) +
              ", component=" + std::to_string(component));
      }
  }

  RegularTaylorTailModel model = *reference;
  model.epsilon = materialized.epsilon;
  model.taylor_complete_max = materialized.taylor_complete_max;
  model.chart = materialized.chart;
  model.prescriptions = materialized.prescriptions;
  model.local_checkpoint_identity = materialized.checkpoint_identity;
  model.n_row_sum_upper.clear();
  model.n_row_sum_upper.reserve(model.n_coefficients.size());
  for (const auto& matrix : model.n_coefficients)
    model.n_row_sum_upper.push_back(
        matrix_infinity_norm_upper(matrix, model.dimension));
  model.initial_row_upper.assign(
      materialized.epsilon.width(), Magnitude::zero());
  for (std::size_t epsilon = 0;
       epsilon < materialized.epsilon.width(); ++epsilon)
    for (std::uint32_t component = 0;
         component < materialized.dimension; ++component)
      model.initial_row_upper[epsilon] = Magnitude::maximum(
          model.initial_row_upper[epsilon],
          Magnitude::upper_abs(local_detail::to_ball(
              output_sector.coefficients[local_detail::sector_index(
                  materialized, epsilon, 0, component)])));
  if (std::any_of(model.initial_row_upper.begin(),
                  model.initial_row_upper.end(),
                  [](const Magnitude& value) {
                    return !value.is_finite();
                  }))
    return inconclusive_model(
        "materialized initial-value enclosure is nonfinite");
  model.provenance =
      "regular homogeneous tail model derived from a native retained "
      "match materialization; every receiving basis column was bound to "
      "the identical prepared q/N payload and the finite Laurent "
      "convolution was replayed; this certifies the outgoing local tail, "
      "not a separate whole-path matching-error bound; "
      "operator_identity=" + model.operator_identity +
      "; local_checkpoint_identity=" +
      model.local_checkpoint_identity +
      "; match_checkpoint_identity=" + match_checkpoint_identity;
  return {TailMajorantStatus::Certified, std::move(model),
          "regular homogeneous tail model propagated through retained match materialization"};
}

// Bind a prepared rational row to the full analytic source theorem without
// pretending that the projected scalar satisfies the source vector ODE.  The
// model is attached only at the private direct-output seam of the finite row
// application, and each supplied rational completion must independently
// reproduce every stored Taylor-kernel coefficient by polynomial division.
// Missing completions leave stored integration fully available but make
// full-local certification explicitly unsupported.
template <typename Scalar>
RationalRowLineTailModelResult<Scalar>
derive_rational_row_line_tail_model(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& matrix,
    const LocalSolution<Scalar>& source,
    const RegularTaylorTailModelResult& source_tail,
    const LocalSolution<Scalar>& projected) {
  static_assert(std::is_same_v<Scalar, Rational> ||
                    std::is_same_v<Scalar, ComplexBall>,
                "rational-row tail models support Rational or Acb locals only");
  using namespace tail_majorant_detail;

  validate_local_solution(source, false);
  validate_local_solution(projected, false);
  if (matrix.rows != 1 || matrix.columns != source.dimension ||
      matrix.exact_identity.empty())
    throw std::invalid_argument(
        "rational-row tail derivation received a malformed prepared row");
  if (source_tail.status != TailMajorantStatus::Certified ||
      !source_tail.model.has_value()) {
    const auto detail =
        "rational-row line tail requires a certified ordinary source local: " +
        source_tail.detail;
    return source_tail.status == TailMajorantStatus::Inconclusive
        ? inconclusive_rational_row_model<Scalar>(detail)
        : unsupported_rational_row_model<Scalar>(detail);
  }
  try {
    require_model_binding(*source_tail.model, source);
  } catch (const std::invalid_argument&) {
    return inconclusive_rational_row_model<Scalar>(
        "rational-row source tail model is not immutably bound to its retained local");
  }
  if (!source.error.empty() || !projected.error.empty())
    return unsupported_rational_row_model<Scalar>(
        "rational-row tail derivation cannot discard an existing local error envelope");
  if (projected.dimension != 1 ||
      projected.taylor_complete_max != source.taylor_complete_max ||
      !same_chart_geometry(projected.chart, source.chart) ||
      !same_prescriptions(projected.prescriptions, source.prescriptions))
    return unsupported_rational_row_model<Scalar>(
        "rational-row projection left its source chart or complete Taylor order");

  // build_rational_row_local calls this immediately on the direct output of
  // apply_prepared_sparse_local_matrix, before either object is moved.  Do not
  // repeat that potentially dominant convolution merely for certification;
  // instead recheck its exact frame contract here, while the private StoredLocal
  // attachment gate prevents this model from being assigned to any other
  // derivation kind.  The analytic completion itself is still independently
  // replayed coefficient by coefficient below.
  std::int32_t expected_min = source.epsilon.min_power;
  std::int32_t expected_complete = source.epsilon.complete_max;
  if (!matrix.entries.empty()) {
    expected_min = std::numeric_limits<std::int32_t>::max();
    expected_complete = std::numeric_limits<std::int32_t>::max();
    for (const auto& entry : matrix.entries) {
      expected_min = std::min(expected_min,
          local_algebra_detail::checked_i32(
              static_cast<std::int64_t>(source.epsilon.min_power) +
                  entry.multiplier.epsilon_shift,
              "rational-row tail minimum"));
      expected_complete = std::min(expected_complete,
          local_algebra_detail::checked_i32(
              static_cast<std::int64_t>(source.epsilon.complete_max) +
                  entry.multiplier.epsilon_shift,
              "rational-row tail complete maximum"));
    }
  }
  if (projected.epsilon.min_power != expected_min ||
      projected.epsilon.complete_max > expected_complete)
    return inconclusive_rational_row_model<Scalar>(
        "retained rational-row epsilon frame is not a complete restriction of its direct finite product");

  RationalRowLineTailModel<Scalar> model;
  model.epsilon = projected.epsilon;
  model.taylor_complete_max = projected.taylor_complete_max;
  model.chart = projected.chart;
  model.prescriptions = projected.prescriptions;
  model.source = *source_tail.model;
  model.row_exact_identity = matrix.exact_identity;
  model.source_checkpoint_identity = source.checkpoint_identity;
  model.local_checkpoint_identity = projected.checkpoint_identity;
  model.entries.reserve(matrix.entries.size());

  for (const auto& entry : matrix.entries) {
    if (entry.row != 0 || entry.column >= source.dimension)
      throw std::invalid_argument(
          "rational-row tail entry lies outside its scalar row");
    const auto& multiplier = entry.multiplier;
    if (!multiplier.analytic_coefficients.has_value())
      return unsupported_rational_row_model<Scalar>(
          "prepared rational-row entry has only a finite Taylor kernel; no analytic numerator/denominator completion was retained");
    const auto& analytic = *multiplier.analytic_coefficients;
    if (analytic.size() != multiplier.kernels.size())
      return inconclusive_rational_row_model<Scalar>(
          "prepared rational-row analytic completion count differs from its epsilon kernels");
    for (std::size_t epsilon = 0; epsilon < analytic.size(); ++epsilon) {
      const auto& rational = analytic[epsilon];
      if (rational.numerator.empty() || rational.denominator.empty() ||
          std::any_of(rational.numerator.begin(), rational.numerator.end(),
                      [](const Scalar& value) {
                        return !finite_scalar(value);
                      }) ||
          std::any_of(rational.denominator.begin(), rational.denominator.end(),
                      [](const Scalar& value) {
                        return !finite_scalar(value);
                      }))
        return inconclusive_rational_row_model<Scalar>(
            "prepared rational-row analytic completion is empty or nonfinite");
      if (divisor_contains_zero(rational.denominator.front()))
        return inconclusive_rational_row_model<Scalar>(
            "prepared rational-row analytic denominator contains zero at the chart center");
      std::vector<Scalar> replayed_kernel;
      replayed_kernel.reserve(source.taylor_width());
      for (std::size_t n = 0; n < source.taylor_width(); ++n) {
        auto rhs = n < rational.numerator.size()
            ? rational.numerator[n] : ScalarTraits<Scalar>::zero();
        const auto max_lag = std::min(n, rational.denominator.size() - 1);
        for (std::size_t lag = 1; lag <= max_lag; ++lag)
          rhs -= rational.denominator[lag] * replayed_kernel[n - lag];
        auto coefficient = rhs / rational.denominator.front();
        if (n >= multiplier.kernels[epsilon].size() ||
            !encloses_recurrence_value(
                multiplier.kernels[epsilon][n], coefficient))
          return inconclusive_rational_row_model<Scalar>(
              "prepared rational-row analytic completion does not reproduce every retained Taylor-kernel coefficient");
        replayed_kernel.push_back(std::move(coefficient));
      }
    }
    model.entries.push_back(RationalRowLineTailEntryModel<Scalar>{
        entry.column, multiplier.epsilon_shift,
        multiplier.center_pole_order, analytic,
        multiplier.exact_identity});
  }
  model.provenance =
      "certified rational-row line-tail model; source ordinary homogeneous "
      "q/N Cauchy theorem retained; full analytic numerator/denominator "
      "payload replayed through every finite multiplier kernel; projected "
      "scalar receives no homogeneous q/C equation ownership; "
      "source_operator_identity=" + model.source.operator_identity +
      "; source_checkpoint_identity=" + model.source_checkpoint_identity +
      "; row_exact_identity=" + model.row_exact_identity +
      "; local_checkpoint_identity=" + model.local_checkpoint_identity;
  return {TailMajorantStatus::Certified, std::move(model),
          "rational-row analytic completion is bound to its certified source and finite projection"};
}

inline RegularTaylorDiskCertificate certify_regular_taylor_disk(
    const RegularTaylorTailModel& model,
    const std::string& witness_radius_exact) {
  using namespace tail_majorant_detail;
  Rational radius(0);
  try {
    radius = Rational(witness_radius_exact);
  } catch (const std::invalid_argument&) {
    throw std::invalid_argument(
        "tail witness radius must be an exact rational string");
  }
  if (radius.sign() <= 0)
    throw std::invalid_argument("tail witness radius must be positive");
  const auto radius_ball = ComplexBall::from_strings(radius.str());
  if (!model.chart.infinite_radius &&
      !arb_lt(acb_realref(radius_ball.raw()),
              acb_realref(model.chart.radius.raw())))
    return inconclusive_disk(
        radius.str(),
        "witness disk is not provably strictly inside the retained chart "
        "radius");

  const auto radius_upper = exact_rational_upper(radius);
  auto q_remainder_upper = Magnitude::zero();
  for (std::size_t lag = 1; lag < model.q_coefficients.size(); ++lag)
    q_remainder_upper +=
        Magnitude::upper_abs(model.q_coefficients[lag]) *
        radius_upper.power_upper(static_cast<ulong>(lag));
  const auto q_lower = Magnitude::positive_difference_lower(
      Magnitude::lower_abs(model.q_coefficients.front()),
      q_remainder_upper);
  if (q_lower.is_zero() || !q_lower.is_finite())
    return inconclusive_disk(
        radius.str(),
        "triangle inequality does not prove q(t) nonzero on the witness "
        "disk; no full-function claim is made");

  auto numerator_upper = Magnitude::zero();
  for (std::size_t lag = 1; lag < model.n_row_sum_upper.size(); ++lag)
    numerator_upper += model.n_row_sum_upper[lag] *
        radius_upper.power_upper(static_cast<ulong>(lag - 1));
  const auto ode_norm_upper = numerator_upper / q_lower;
  const auto gronwall =
      (ode_norm_upper * radius_upper).exponential_upper();
  if (!ode_norm_upper.is_finite() || !gronwall.is_finite())
    return inconclusive_disk(
        radius.str(),
        "ODE/Gronwall magnitude overflowed; no finite full-function bound "
        "is available at this precision/radius");

  RegularTaylorDiskCertificate result;
  result.status = TailMajorantStatus::Certified;
  result.witness_radius_exact = radius.str();
  result.q_lower = q_lower;
  result.ode_norm_upper = ode_norm_upper;
  result.cauchy_circle_upper.reserve(model.initial_row_upper.size());
  for (const auto& initial : model.initial_row_upper)
    result.cauchy_circle_upper.push_back(initial * gronwall);
  result.detail =
      "q is separated from zero on the witness disk and Gronwall bounds "
      "the regular solution on its Cauchy circle";
  return result;
}

inline RegularTaylorPointTailCertificate certify_regular_taylor_point_tail(
    const RegularTaylorTailModel& model,
    const RealEvaluationPoint& input_point,
    const std::string& witness_radius_exact,
    std::optional<std::uint32_t> retained_complete_max = std::nullopt) {
  using namespace tail_majorant_detail;
  const auto point = line_integration_detail::require_exact_rational_point(
      input_point, "tail-evaluation");
  const Rational exact_point(point.exact_coordinate);
  const auto modulus = abs_rational(exact_point);
  const Rational radius(witness_radius_exact);
  const auto retained = retained_complete_max.value_or(
      model.taylor_complete_max);
  if (retained > model.taylor_complete_max)
    throw std::invalid_argument(
        "tail retained order exceeds the model's complete Taylor order");

  RegularTaylorPointTailCertificate result;
  result.disk = certify_regular_taylor_disk(model, radius.str());
  if (result.disk.status != TailMajorantStatus::Certified) {
    result.status = result.disk.status;
    result.detail = result.disk.detail;
    return result;
  }
  if (!(modulus < radius)) {
    result.status = TailMajorantStatus::Inconclusive;
    result.detail =
        "evaluation point is not strictly inside the tail witness disk";
    return result;
  }

  const auto radius_lower = exact_rational_lower(radius);
  if (radius_lower.is_zero()) {
    result.status = TailMajorantStatus::Inconclusive;
    result.detail =
        "finite-precision radius enclosure has no positive lower bound";
    return result;
  }
  const auto ratio_upper = exact_rational_upper(modulus) / radius_lower;
  Magnitude gap_lower;
  const auto first = static_cast<ulong>(retained) + 1UL;
  const auto value_factor = geometric_tail_factor(
      ratio_upper, first, &gap_lower);
  if (gap_lower.is_zero()) {
    result.status = TailMajorantStatus::Inconclusive;
    result.detail =
        "rounded |t|/R bound does not remain strictly below one";
    return result;
  }
  const auto inverse_gap = gap_lower.reciprocal_upper();
  const auto theta_factor = ratio_upper.power_upper(first) *
      (Magnitude::from_ui(first) * inverse_gap +
       ratio_upper * inverse_gap * inverse_gap);

  result.value.frame = model.epsilon;
  result.value.guarantee = ErrorGuarantee::Certified;
  result.theta.frame = model.epsilon;
  result.theta.guarantee = ErrorGuarantee::Certified;
  result.value.absolute.reserve(result.disk.cauchy_circle_upper.size());
  result.theta.absolute.reserve(result.disk.cauchy_circle_upper.size());
  for (const auto& circle : result.disk.cauchy_circle_upper) {
    result.value.absolute.push_back(circle * value_factor);
    result.theta.absolute.push_back(circle * theta_factor);
  }
  result.value.provenance = certificate_provenance(
      model, radius.str(), "point-value");
  result.theta.provenance = certificate_provenance(
      model, radius.str(), "point-theta");
  result.status = TailMajorantStatus::Certified;
  result.detail =
      "Cauchy coefficient bounds certify every Taylor coefficient beyond "
      "the retained order";
  return result;
}

inline RegularTaylorDiskCertificate certify_physical_regular_taylor_disk(
    const PhysicalRegularTaylorTailModel& model,
    const std::string& witness_radius_exact) {
  using namespace tail_majorant_detail;
  Rational radius(0);
  try {
    radius = Rational(witness_radius_exact);
  } catch (const std::invalid_argument&) {
    throw std::invalid_argument(
        "physical tail witness radius must be an exact rational string");
  }
  if (radius.sign() <= 0)
    throw std::invalid_argument(
        "physical tail witness radius must be positive");
  const auto radius_ball = ComplexBall::from_strings(radius.str());
  if (!model.chart.infinite_radius &&
      !arb_lt(acb_realref(radius_ball.raw()),
              acb_realref(model.chart.radius.raw())))
    return inconclusive_disk(
        radius.str(),
        "physical witness disk is not provably strictly inside the retained chart radius");
  const auto width = model.epsilon.width();
  const auto complete_prefix_family = [width](const auto& family) {
    return !family.empty() &&
        std::all_of(family.begin(), family.end(),
                    [width](const auto& prefix) {
                      return prefix.size() == width;
                    });
  };
  if (!complete_prefix_family(
          model.q_operator_prefix_norm_upper) ||
      !complete_prefix_family(
          model.q_causal_coefficients) ||
      !complete_prefix_family(
          model.c_operator_prefix_norm_upper) ||
      model.c_causal_coefficients.empty() ||
      !std::all_of(
          model.c_causal_coefficients.begin(),
          model.c_causal_coefficients.end(),
          [width, dimension = model.dimension](
              const auto& lag) {
            return lag.size() ==
                width * static_cast<std::size_t>(dimension) *
                    dimension;
          }) ||
      model.q0_inverse_prefix_norm_upper.size() != width ||
      model.initial_row_upper.size() != width)
    throw std::invalid_argument(
        "physical tail model has incomplete operator or initial bounds");

  const auto radius_upper = exact_rational_upper(radius);
  const auto normalized =
      certify_physical_normalized_ode_disk_bounds(
          model.q_causal_coefficients,
          model.c_causal_coefficients, model.dimension, radius);
  if (normalized.status != TailMajorantStatus::Certified)
    return inconclusive_disk(
        radius.str(),
        "finite-epsilon normalized physical ODE disk bound is "
        "inconclusive: " +
            normalized.detail);
  const std::vector<ulong> epsilon_weights{
      1UL, 2UL, 4UL, 8UL, 16UL, 32UL, 64UL, 128UL,
      256UL, 512UL, 1024UL};
  std::vector<std::vector<Magnitude>> weighted_ode_norms;
  weighted_ode_norms.reserve(epsilon_weights.size());
  for (const auto weight : epsilon_weights)
    weighted_ode_norms.push_back(
        weighted_physical_ode_prefix_norm_upper(
            normalized, weight));
  RegularTaylorDiskCertificate result;
  result.status = TailMajorantStatus::Certified;
  result.witness_radius_exact = radius.str();
  auto maximum_ode_norm = Magnitude::zero();
  ulong minimum_selected_weight =
      std::numeric_limits<ulong>::max();
  ulong maximum_selected_weight = 0;
  result.cauchy_circle_upper.reserve(width);
  for (std::size_t epsilon_index = 0;
       epsilon_index < width; ++epsilon_index) {
    std::optional<Magnitude> best_circle;
    Magnitude best_ode_norm = Magnitude::zero();
    ulong best_weight = 0;
    for (std::size_t candidate = 0;
         candidate < epsilon_weights.size(); ++candidate) {
      const auto weight = epsilon_weights[candidate];
      const auto epsilon_weight = Magnitude::from_ui(weight);
      auto weighted_initial = Magnitude::zero();
      for (std::size_t source_index = 0;
           source_index <= epsilon_index; ++source_index)
        weighted_initial = Magnitude::maximum(
            weighted_initial,
            model.initial_row_upper[source_index] /
                epsilon_weight.power_upper(
                    static_cast<ulong>(source_index)));
      const auto& ode_norm =
          weighted_ode_norms[candidate][epsilon_index];
      const auto gronwall =
          (ode_norm * radius_upper).exponential_upper();
      const auto circle =
          weighted_initial * gronwall *
          epsilon_weight.power_upper(
              static_cast<ulong>(epsilon_index));
      if (!circle.is_finite())
        continue;
      if (!best_circle.has_value() || circle <= *best_circle) {
        best_circle = circle;
        best_ode_norm = ode_norm;
        best_weight = weight;
      }
    }
    if (!best_circle.has_value())
      return inconclusive_disk(
          radius.str(),
          "physical normalized finite-epsilon ODE/Gronwall bound "
          "overflowed for every tested geometric epsilon weight on the "
          "causal prefix ending at epsilon power " +
              std::to_string(
                  static_cast<std::int64_t>(
                      model.epsilon.min_power) +
                  static_cast<std::int64_t>(epsilon_index)));
    maximum_ode_norm =
        Magnitude::maximum(maximum_ode_norm, best_ode_norm);
    minimum_selected_weight =
        std::min(minimum_selected_weight, best_weight);
    maximum_selected_weight =
        std::max(maximum_selected_weight, best_weight);
    result.cauchy_circle_upper.push_back(*best_circle);
  }
  // These scalar fields remain conservative summaries for diagnostics and
  // checkpoint replay.  The actual coefficient bounds above use the
  // corresponding invariant prefix, not these extrema.
  result.q_lower = normalized.q_diagonal_lower;
  result.ode_norm_upper = maximum_ode_norm;
  result.detail =
      "Q(t,epsilon)^-1 C(t,epsilon)/t is interval-evaluated before norms on "
      "an adaptive disk cover; independently optimized geometric epsilon "
      "weights and prefixwise Gronwall bound every invariant causal prefix; "
      "accepted_tiles=" +
      std::to_string(normalized.accepted_tiles) +
      "; subdivided_tiles=" +
      std::to_string(normalized.subdivided_tiles) +
      "; selected_epsilon_weight_range=" +
      std::to_string(minimum_selected_weight) + ".." +
      std::to_string(maximum_selected_weight);
  return result;
}

inline RegularTaylorPointTailCertificate
certify_physical_regular_taylor_point_tail(
    const PhysicalRegularTaylorTailModel& model,
    const RealEvaluationPoint& input_point,
    const std::string& witness_radius_exact,
    std::optional<std::uint32_t> retained_complete_max = std::nullopt) {
  using namespace tail_majorant_detail;
  const auto point = line_integration_detail::require_exact_rational_point(
      input_point, "physical-tail-evaluation");
  const Rational exact_point(point.exact_coordinate);
  const auto modulus = abs_rational(exact_point);
  const Rational radius(witness_radius_exact);
  const auto retained = retained_complete_max.value_or(
      model.taylor_complete_max);
  if (retained > model.taylor_complete_max)
    throw std::invalid_argument(
        "physical tail retained order exceeds the model order");

  RegularTaylorPointTailCertificate result;
  result.disk = certify_physical_regular_taylor_disk(
      model, radius.str());
  if (result.disk.status != TailMajorantStatus::Certified) {
    result.status = result.disk.status;
    result.detail = result.disk.detail;
    return result;
  }
  if (!(modulus < radius)) {
    result.status = TailMajorantStatus::Inconclusive;
    result.detail =
        "physical evaluation point is not strictly inside the witness disk";
    return result;
  }
  const auto radius_lower = exact_rational_lower(radius);
  if (radius_lower.is_zero()) {
    result.status = TailMajorantStatus::Inconclusive;
    result.detail =
        "physical witness-radius enclosure has no positive lower bound";
    return result;
  }
  const auto ratio_upper = exact_rational_upper(modulus) / radius_lower;
  Magnitude gap_lower;
  const auto first = static_cast<ulong>(retained) + 1UL;
  const auto value_factor = geometric_tail_factor(
      ratio_upper, first, &gap_lower);
  if (gap_lower.is_zero()) {
    result.status = TailMajorantStatus::Inconclusive;
    result.detail =
        "physical rounded |t|/R bound does not remain below one";
    return result;
  }
  const auto inverse_gap = gap_lower.reciprocal_upper();
  const auto theta_factor = ratio_upper.power_upper(first) *
      (Magnitude::from_ui(first) * inverse_gap +
       ratio_upper * inverse_gap * inverse_gap);
  result.value.frame = model.epsilon;
  result.value.guarantee = ErrorGuarantee::Certified;
  result.theta.frame = model.epsilon;
  result.theta.guarantee = ErrorGuarantee::Certified;
  for (const auto& circle : result.disk.cauchy_circle_upper) {
    result.value.absolute.push_back(circle * value_factor);
    result.theta.absolute.push_back(circle * theta_factor);
  }
  const auto provenance =
      "certified full-local point tail; augmented finite-epsilon physical "
      "ODE; witness_radius_exact=" + radius.str() + "; " +
      model.provenance;
  result.value.provenance = provenance;
  result.theta.provenance = provenance;
  result.status = TailMajorantStatus::Certified;
  result.detail =
      "Cauchy bounds certify the unseen Taylor tail of the finite-epsilon physical solution";
  return result;
}

inline CertifiedLocalEvaluation
evaluate_physical_local_solution_with_certified_tail(
    const PhysicalRegularTaylorTailModel& model,
    const RealEvaluationPoint& point,
    const std::string& witness_radius_exact,
    EvaluationOptions options = {}) {
  if (model.reconstructed.checkpoint_identity !=
          model.local_checkpoint_identity ||
      model.reconstructed.dimension != model.dimension ||
      !tail_majorant_detail::same_epsilon_window(
          model.reconstructed.epsilon, model.epsilon))
    throw std::invalid_argument(
        "physical tail model is not bound to its reconstructed local");
  if (options.t_order_reduction > model.taylor_complete_max)
    throw std::invalid_argument(
        "physical tail-aware evaluation reduction exceeds retained order");
  const auto retained = model.taylor_complete_max -
                        options.t_order_reduction;
  options.compute_tail_estimate = false;
  CertifiedLocalEvaluation result;
  result.evaluation = evaluate_local_solution(
      model.reconstructed, point, options);
  result.tail = certify_physical_regular_taylor_point_tail(
      model, point, witness_radius_exact, retained);
  if (result.tail.status != TailMajorantStatus::Certified) {
    result.evaluation.value.error.provenance = result.tail.detail;
    result.evaluation.theta_value.error.provenance = result.tail.detail;
    return result;
  }
  result.evaluation.value.error = result.tail.value;
  ErrorEnvelope theta_error;
  theta_error.frame = result.evaluation.theta_value.epsilon;
  theta_error.guarantee = ErrorGuarantee::Certified;
  theta_error.provenance = result.tail.theta.provenance;
  for (std::int64_t raw_power = theta_error.frame.min_power;
       raw_power <= theta_error.frame.complete_max; ++raw_power)
    theta_error.absolute.push_back(result.tail.theta.absolute.at(
        static_cast<std::size_t>(
            raw_power - result.tail.theta.frame.min_power)));
  result.evaluation.theta_value.error = std::move(theta_error);
  return result;
}

template <typename Scalar>
CertifiedLocalEvaluation evaluate_local_solution_with_certified_tail(
    const LocalSolution<Scalar>& solution,
    const RegularTaylorTailModel& model,
    const RealEvaluationPoint& point,
    const std::string& witness_radius_exact,
    EvaluationOptions options = {}) {
  using namespace tail_majorant_detail;
  require_model_binding(model, solution);
  if (options.t_order_reduction > solution.taylor_complete_max)
    throw std::invalid_argument(
        "tail-aware evaluation Taylor reduction exceeds retained order");
  const auto retained = solution.taylor_complete_max -
                        options.t_order_reduction;
  options.compute_tail_estimate = false;
  CertifiedLocalEvaluation result;
  result.evaluation = evaluate_local_solution(solution, point, options);
  result.tail = certify_regular_taylor_point_tail(
      model, point, witness_radius_exact, retained);
  if (result.tail.status == TailMajorantStatus::Certified) {
    if (!same_epsilon_window(result.evaluation.value.epsilon,
                             result.tail.value.frame) ||
        result.evaluation.theta_value.epsilon.min_power <
            result.tail.theta.frame.min_power ||
        result.evaluation.theta_value.epsilon.complete_max >
            result.tail.theta.frame.complete_max)
      throw std::logic_error(
          "ordinary tail certificate and evaluation epsilon windows differ");
    result.evaluation.value.error = result.tail.value;
    // A structurally zero theta value may canonically trim leading epsilon
    // rows even though the theorem certifies the full retained frame.  Attach
    // the exact certified sub-envelope that corresponds to the evaluated
    // theta vector; the full theorem remains available in result.tail.
    ErrorEnvelope theta_error;
    theta_error.frame = result.evaluation.theta_value.epsilon;
    theta_error.guarantee = ErrorGuarantee::Certified;
    theta_error.provenance = result.tail.theta.provenance;
    theta_error.absolute.reserve(theta_error.frame.width());
    for (std::int64_t raw_power = theta_error.frame.min_power;
         raw_power <= theta_error.frame.complete_max; ++raw_power)
      theta_error.absolute.push_back(result.tail.theta.absolute.at(
          static_cast<std::size_t>(
              raw_power - result.tail.theta.frame.min_power)));
    result.evaluation.theta_value.error = std::move(theta_error);
  } else {
    result.evaluation.value.error.provenance =
        status_prefix(result.tail.status) + ": " + result.tail.detail;
    result.evaluation.theta_value.error.provenance =
        result.evaluation.value.error.provenance;
  }
  return result;
}

inline RegularTaylorLineTailCertificate certify_regular_taylor_line_tail(
    const RegularTaylorTailModel& model,
    const RealEvaluationPoint& lower_input,
    const RealEvaluationPoint& upper_input,
    const std::string& witness_radius_exact,
    const EpsilonWindow& delivered_epsilon) {
  using namespace tail_majorant_detail;
  const auto lower = line_integration_detail::require_exact_rational_point(
      lower_input, "tail-line lower");
  const auto upper = line_integration_detail::require_exact_rational_point(
      upper_input, "tail-line upper");
  integration_detail::validate_interval(lower, upper);
  (void)delivered_epsilon.width();
  if (delivered_epsilon.min_power < model.epsilon.min_power) {
    RegularTaylorLineTailCertificate result;
    result.status = TailMajorantStatus::Unsupported;
    result.detail =
        "full-local line certificate does not infer structural epsilon "
        "zeros below the retained model window";
    return result;
  }
  if (delivered_epsilon.complete_max > model.epsilon.complete_max)
    throw std::invalid_argument(
        "tail line requests epsilon rows beyond the certified model window");

  const Rational radius(witness_radius_exact);
  RegularTaylorLineTailCertificate result;
  result.disk = certify_regular_taylor_disk(model, radius.str());
  if (result.disk.status != TailMajorantStatus::Certified) {
    result.status = result.disk.status;
    result.detail = result.disk.detail;
    return result;
  }

  const auto lower_modulus = abs_rational(Rational(lower.exact_coordinate));
  const auto upper_modulus = abs_rational(Rational(upper.exact_coordinate));
  if (!(lower_modulus < radius) || !(upper_modulus < radius)) {
    result.status = TailMajorantStatus::Inconclusive;
    result.detail =
        "line endpoint is not strictly inside the tail witness disk";
    return result;
  }
  const auto radius_lower = exact_rational_lower(radius);
  const auto radius_upper = exact_rational_upper(radius);
  if (radius_lower.is_zero()) {
    result.status = TailMajorantStatus::Inconclusive;
    result.detail =
        "finite-precision radius enclosure has no positive lower bound";
    return result;
  }

  const auto arm_factor = [&](const Rational& arm_radius,
                              bool* certified) {
    if (arm_radius.is_zero()) return Magnitude::zero();
    const auto ratio = exact_rational_upper(arm_radius) / radius_lower;
    Magnitude gap;
    const auto geometric = geometric_tail_factor(
        ratio, static_cast<ulong>(model.taylor_complete_max) + 2UL, &gap);
    if (gap.is_zero()) {
      *certified = false;
      return Magnitude::zero();
    }
    return radius_upper * geometric /
           Magnitude::from_ui(
               static_cast<ulong>(model.taylor_complete_max) + 2UL);
  };

  bool factor_certified = true;
  Magnitude integral_factor = Magnitude::zero();
  if (lower.sign < 0 && upper.sign > 0) {
    integral_factor += arm_factor(lower_modulus, &factor_certified);
    integral_factor += arm_factor(upper_modulus, &factor_certified);
  } else {
    const auto outer = lower_modulus < upper_modulus
        ? upper_modulus : lower_modulus;
    integral_factor = arm_factor(outer, &factor_certified);
  }
  if (!factor_certified) {
    result.status = TailMajorantStatus::Inconclusive;
    result.detail =
        "rounded line endpoint ratio does not remain strictly below one";
    return result;
  }

  result.integral.frame = delivered_epsilon;
  result.integral.guarantee = ErrorGuarantee::Certified;
  result.integral.absolute.reserve(delivered_epsilon.width());
  for (std::int64_t power64 = delivered_epsilon.min_power;
       power64 <= delivered_epsilon.complete_max; ++power64) {
    const auto power = static_cast<std::int32_t>(power64);
    if (power < model.epsilon.min_power) {
      result.integral.absolute.push_back(Magnitude::zero());
      continue;
    }
    const auto index = static_cast<std::size_t>(
        static_cast<std::int64_t>(power) - model.epsilon.min_power);
    result.integral.absolute.push_back(
        result.disk.cauchy_circle_upper.at(index) * integral_factor);
  }
  result.integral.provenance = certificate_provenance(
      model, radius.str(), "line-integral");
  result.status = TailMajorantStatus::Certified;
  result.detail =
      "integrated Cauchy coefficient majorant certifies the unseen Taylor "
      "tail on the complete line interval";
  return result;
}

template <typename Scalar>
RegularTaylorLineTailCertificate certify_rational_row_line_tail(
    const RationalRowLineTailModel<Scalar>& model,
    const LocalSolution<Scalar>& projected,
    const RealEvaluationPoint& lower_input,
    const RealEvaluationPoint& upper_input,
    const std::string& witness_radius_exact,
    const EpsilonWindow& delivered_epsilon) {
  static_assert(std::is_same_v<Scalar, Rational> ||
                    std::is_same_v<Scalar, ComplexBall>,
                "rational-row line tails support Rational or Acb locals only");
  using namespace tail_majorant_detail;

  validate_local_solution(projected, false);
  if (projected.checkpoint_identity != model.local_checkpoint_identity ||
      projected.dimension != 1 ||
      projected.epsilon.min_power != model.epsilon.min_power ||
      projected.epsilon.complete_max != model.epsilon.complete_max ||
      projected.taylor_complete_max != model.taylor_complete_max ||
      !same_chart_geometry(projected.chart, model.chart) ||
      !same_prescriptions(projected.prescriptions, model.prescriptions))
    throw std::invalid_argument(
        "rational-row line-tail model is not bound to this retained scalar local");
  (void)delivered_epsilon.width();
  // The fused primitive contract conservatively reserves one lower epsilon
  // row for a possible regulated-center pole.  For the unseen regular Taylor
  // tail, first_exponent > 0 is certified below before any epsilon
  // convolution is accumulated, so rows below the projected model window are
  // rigorously zero.  The accumulation loop already represents them by an
  // empty source sum.  Upper rows, in contrast, would require unavailable
  // source coefficients and remain a hard error.
  if (delivered_epsilon.complete_max > model.epsilon.complete_max)
    throw std::invalid_argument(
        "rational-row line tail requests epsilon rows outside its complete projected window");

  const auto lower = line_integration_detail::require_exact_rational_point(
      lower_input, "rational-row tail-line lower");
  const auto upper = line_integration_detail::require_exact_rational_point(
      upper_input, "rational-row tail-line upper");
  integration_detail::validate_interval(lower, upper);
  const Rational radius(witness_radius_exact);
  RegularTaylorLineTailCertificate result;
  result.disk = certify_regular_taylor_disk(model.source, radius.str());
  if (result.disk.status != TailMajorantStatus::Certified) {
    result.status = result.disk.status;
    result.detail = result.disk.detail;
    return result;
  }

  const auto lower_modulus = abs_rational(Rational(lower.exact_coordinate));
  const auto upper_modulus = abs_rational(Rational(upper.exact_coordinate));
  if (!(lower_modulus < radius) || !(upper_modulus < radius)) {
    result.status = TailMajorantStatus::Inconclusive;
    result.detail =
        "rational-row line endpoint is not strictly inside the tail witness disk";
    return result;
  }
  const auto radius_lower = exact_rational_lower(radius);
  const auto radius_upper = exact_rational_upper(radius);
  if (radius_lower.is_zero()) {
    result.status = TailMajorantStatus::Inconclusive;
    result.detail =
        "finite-precision rational-row witness radius has no positive lower bound";
    return result;
  }

  const auto rational_circle_upper = [&](const auto& rational,
                                         bool* certified) {
    auto numerator_upper = Magnitude::zero();
    for (std::size_t degree = 0; degree < rational.numerator.size(); ++degree)
      numerator_upper += Magnitude::upper_abs(
          local_detail::to_ball(rational.numerator[degree])) *
          radius_upper.power_upper(static_cast<ulong>(degree));
    auto denominator_remainder = Magnitude::zero();
    for (std::size_t degree = 1; degree < rational.denominator.size(); ++degree)
      denominator_remainder += Magnitude::upper_abs(
          local_detail::to_ball(rational.denominator[degree])) *
          radius_upper.power_upper(static_cast<ulong>(degree));
    const auto denominator_lower = Magnitude::positive_difference_lower(
        Magnitude::lower_abs(
            local_detail::to_ball(rational.denominator.front())),
        denominator_remainder);
    if (denominator_lower.is_zero() || !denominator_lower.is_finite()) {
      *certified = false;
      return Magnitude::zero();
    }
    const auto bound = numerator_upper / denominator_lower;
    if (!bound.is_finite()) {
      *certified = false;
      return Magnitude::zero();
    }
    return bound;
  };

  const auto arm_factor = [&](const Rational& arm_radius,
                              std::uint32_t pole_order,
                              bool* certified) {
    if (arm_radius.is_zero()) return Magnitude::zero();
    const auto first_exponent =
        static_cast<std::int64_t>(model.taylor_complete_max) + 2 -
        static_cast<std::int64_t>(pole_order);
    if (first_exponent <= 0) {
      *certified = false;
      return Magnitude::zero();
    }
    const auto ratio = exact_rational_upper(arm_radius) / radius_lower;
    Magnitude gap;
    const auto geometric = geometric_tail_factor(
        ratio, static_cast<ulong>(first_exponent), &gap);
    if (gap.is_zero()) {
      *certified = false;
      return Magnitude::zero();
    }
    Magnitude radius_scale;
    if (pole_order == 0) {
      radius_scale = radius_upper;
    } else if (pole_order == 1) {
      radius_scale = Magnitude::one();
    } else {
      Rational positive_power(1);
      for (std::uint32_t count = 1; count < pole_order; ++count)
        positive_power *= radius;
      const auto lower_power = exact_rational_lower(positive_power);
      if (lower_power.is_zero()) {
        *certified = false;
        return Magnitude::zero();
      }
      radius_scale = lower_power.reciprocal_upper();
    }
    return radius_scale * geometric /
           Magnitude::from_ui(static_cast<ulong>(first_exponent));
  };

  struct EntryBounds {
    const RationalRowLineTailEntryModel<Scalar>* entry = nullptr;
    std::vector<Magnitude> rational_circle;
    Magnitude integral_factor = Magnitude::zero();
  };
  std::vector<EntryBounds> entry_bounds;
  entry_bounds.reserve(model.entries.size());
  bool factors_certified = true;
  for (const auto& entry : model.entries) {
    EntryBounds bounds;
    bounds.entry = &entry;
    bounds.rational_circle.reserve(entry.analytic_coefficients.size());
    for (const auto& rational : entry.analytic_coefficients)
      bounds.rational_circle.push_back(
          rational_circle_upper(rational, &factors_certified));
    if (lower.sign < 0 && upper.sign > 0) {
      bounds.integral_factor += arm_factor(
          lower_modulus, entry.center_pole_order, &factors_certified);
      bounds.integral_factor += arm_factor(
          upper_modulus, entry.center_pole_order, &factors_certified);
    } else {
      const auto outer = lower_modulus < upper_modulus
          ? upper_modulus : lower_modulus;
      bounds.integral_factor = arm_factor(
          outer, entry.center_pole_order, &factors_certified);
    }
    entry_bounds.push_back(std::move(bounds));
  }
  if (!factors_certified) {
    result.status = TailMajorantStatus::Inconclusive;
    result.detail =
        "rational-row denominator separation, endpoint ratio, or first-unseen integrability bound was inconclusive";
    return result;
  }

  result.integral.frame = delivered_epsilon;
  result.integral.guarantee = ErrorGuarantee::Certified;
  result.integral.absolute.reserve(delivered_epsilon.width());
  for (std::int64_t power64 = delivered_epsilon.min_power;
       power64 <= delivered_epsilon.complete_max; ++power64) {
    auto output_bound = Magnitude::zero();
    for (const auto& bounds : entry_bounds) {
      const auto& entry = *bounds.entry;
      for (std::size_t epsilon = 0;
           epsilon < bounds.rational_circle.size(); ++epsilon) {
        const auto source_power = power64 -
            static_cast<std::int64_t>(entry.epsilon_shift) -
            static_cast<std::int64_t>(epsilon);
        if (source_power < model.source.epsilon.min_power ||
            source_power > model.source.epsilon.complete_max)
          continue;
        const auto source_index = static_cast<std::size_t>(
            source_power - model.source.epsilon.min_power);
        output_bound += bounds.rational_circle[epsilon] *
            result.disk.cauchy_circle_upper.at(source_index) *
            bounds.integral_factor;
      }
    }
    result.integral.absolute.push_back(std::move(output_bound));
  }
  result.integral.provenance =
      "certified full-local rational-row line tail; source Gronwall/Cauchy "
      "circle bound times verified analytic rational multiplier bounds; "
      "entry-wise center-pole integration majorant; strict epsilon "
      "convolution; analytic prescriptions retained in absolute value; " +
      model.provenance + "; witness_radius_exact=" + radius.str();
  result.status = TailMajorantStatus::Certified;
  result.detail =
      "verified rational multiplier and source Cauchy majorants certify every unseen projected Taylor contribution to the line integral";
  return result;
}

template <typename Scalar>
CertifiedStoredLineIntegral
integrate_rational_row_local_line_with_certified_tail(
    const LocalSolution<Scalar>& projected,
    const RationalRowLineTailModel<Scalar>& model,
    const RealEvaluationPoint& lower,
    const RealEvaluationPoint& upper,
    const StoredLineIntegrationOptions& options,
    const std::string& witness_radius_exact) {
  CertifiedStoredLineIntegral result;
  result.integral = integrate_stored_local_line(
      projected, lower, upper, options);
  result.tail = certify_rational_row_line_tail(
      model, projected, lower, upper, witness_radius_exact,
      options.delivered_epsilon);
  result.integral.diagnostics.tail_certificate_requested = true;
  result.integral.diagnostics.tail_certificate_status =
      tail_majorant_detail::status_prefix(result.tail.status);
  result.integral.diagnostics.tail_witness_radius_exact =
      witness_radius_exact;
  if (result.tail.status == TailMajorantStatus::Certified) {
    result.integral.value.error = result.tail.integral;
    result.integral.scope =
        LineIntegrationScope::FullLocalWithCertifiedTail;
    result.integral.diagnostics.detail = result.tail.detail;
  } else {
    result.integral.value.error.provenance =
        tail_majorant_detail::status_prefix(result.tail.status) + ": " +
        result.tail.detail +
        "; returned value remains stored Taylor truncation only";
    result.integral.diagnostics.detail =
        result.integral.value.error.provenance;
  }
  return result;
}

template <typename Scalar>
CertifiedStoredLineIntegral integrate_regular_local_line_with_certified_tail(
    const LocalSolution<Scalar>& solution,
    const RegularTaylorTailModel& model,
    const RealEvaluationPoint& lower,
    const RealEvaluationPoint& upper,
    const StoredLineIntegrationOptions& options,
    const std::string& witness_radius_exact) {
  using namespace tail_majorant_detail;
  require_model_binding(model, solution);
  CertifiedStoredLineIntegral result;
  result.integral = integrate_stored_local_line(
      solution, lower, upper, options);
  result.tail = certify_regular_taylor_line_tail(
      model, lower, upper, witness_radius_exact,
      options.delivered_epsilon);
  result.integral.diagnostics.tail_certificate_requested = true;
  result.integral.diagnostics.tail_certificate_status =
      status_prefix(result.tail.status);
  result.integral.diagnostics.tail_witness_radius_exact =
      witness_radius_exact;
  if (result.tail.status == TailMajorantStatus::Certified) {
    result.integral.value.error = result.tail.integral;
    result.integral.scope =
        LineIntegrationScope::FullLocalWithCertifiedTail;
    result.integral.diagnostics.detail = result.tail.detail;
  } else {
    result.integral.value.error.provenance =
        status_prefix(result.tail.status) + ": " + result.tail.detail +
        "; returned value remains stored Taylor truncation only";
    result.integral.diagnostics.detail =
        result.integral.value.error.provenance;
  }
  return result;
}

}  // namespace diffexp2
