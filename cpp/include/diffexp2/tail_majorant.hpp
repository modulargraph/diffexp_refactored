#pragma once

#include "diffexp2/line_integration.hpp"
#include "diffexp2/matching.hpp"
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

namespace tail_majorant_detail {

inline Rational abs_rational(const Rational& value) {
  return value.sign() < 0 ? -value : value;
}

template <typename Scalar>
bool exact_scalar_zero(const Scalar& value) {
  return ScalarTraits<Scalar>::is_zero(value);
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

inline bool same_epsilon_window(const EpsilonWindow& left,
                                const EpsilonWindow& right) {
  return left.min_power == right.min_power &&
         left.complete_max == right.complete_max;
}

inline bool same_chart_geometry(const ChartGeometry& left,
                                const ChartGeometry& right) {
  return left.center_exact == right.center_exact &&
      left.scale_exact == right.scale_exact &&
      left.infinite_radius == right.infinite_radius &&
      (left.infinite_radius ||
       acb_equal(left.radius.raw(), right.radius.raw()));
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
        !same_epsilon_window(result.evaluation.theta_value.epsilon,
                             result.tail.theta.frame))
      throw std::logic_error(
          "ordinary tail certificate and evaluation epsilon windows differ");
    result.evaluation.value.error = result.tail.value;
    result.evaluation.theta_value.error = result.tail.theta;
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
