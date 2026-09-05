#pragma once

#include "diffexp/kernel/scc_completeness.hpp"
#include "diffexp/kernel/tail_majorant.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace diffexp::kernel {

// A finite Frobenius prefix is not a certified value merely because its
// cleared-ODE residual is exactly zero through the retained order.  This
// transient model additionally proves an all-future contraction for every
// canonical (a,b) tower and every invariant causal epsilon prefix.
struct SingularFrobeniusTowerTailModel {
  std::string a_exact;
  std::string b_exact;
  std::uint32_t log_complete_max = 0;
  ulong epsilon_weight_base = 1;
  // Entry k applies to the causal prefix ending at
  // epsilon.min_power+k.  M is the induction bound
  // max ||u_n|| R^n over the last polynomial-lag window.
  std::vector<Magnitude> induction_m_upper;
  std::vector<Magnitude> contraction_upper;
};

struct SingularFrobeniusTailModel {
  EpsilonWindow epsilon;
  std::uint32_t dimension = 0;
  std::uint32_t taylor_complete_max = 0;
  std::uint32_t polynomial_lag_degree = 0;
  std::string witness_radius_exact;
  ChartGeometry chart;
  std::vector<Prescription> prescriptions;
  std::string physical_payload_identity;
  std::string local_checkpoint_identity;
  std::string provenance;
  std::vector<SingularFrobeniusTowerTailModel> towers;
};

struct SingularFrobeniusTailModelResult {
  TailMajorantStatus status = TailMajorantStatus::Unsupported;
  std::optional<SingularFrobeniusTailModel> model;
  std::string detail;
};

struct IntervalTransferContractionResult {
  TailMajorantStatus status = TailMajorantStatus::Inconclusive;
  std::vector<Magnitude> contraction_upper;
  ulong epsilon_weight_base = 1;
  std::string detail;
};

struct IntervalTransferDegreeBoundsResult {
  TailMajorantStatus status = TailMajorantStatus::Inconclusive;
  // Entry [epsilon_degree][physical_row] is the radius-weighted row sum of
  // that Toeplitz coefficient, including every polynomial lag and finite
  // logarithmic shift. Geometric epsilon weights are applied only after this
  // enclosure has been constructed.
  std::vector<std::vector<Magnitude>> degree_row_upper;
  std::string detail;
};

struct IntervalTransferLagDegreeBoundsResult {
  TailMajorantStatus status = TailMajorantStatus::Inconclusive;
  std::size_t epsilon_width = 0;
  std::size_t physical_dimension = 0;
  // Entry [polynomial_lag][epsilon_degree][physical_row] excludes the witness
  // radius. A certified exact prefix can therefore retain this expensive
  // enclosure while bridge trials vary only their candidate radius.
  std::vector<std::vector<std::vector<Magnitude>>> lag_degree_row_upper;
  std::string detail;
};

struct SingularFrobeniusTowerTransferCertificate {
  std::string a_exact;
  std::string b_exact;
  std::uint32_t log_complete_max = 0;
  IntervalTransferLagDegreeBoundsResult interval_transfer;
};

struct SingularFrobeniusPrefixCertificate {
  EpsilonWindow epsilon;
  std::uint32_t taylor_complete_max = 0;
  std::string physical_payload_identity;
  std::string local_checkpoint_identity;
  SCCFormalResidualCertificate exact_residual;
  // Radius-independent recurrence transfers are bound to the same exact
  // equation/local/window as the residual. Bridge seed/radius searches must
  // reuse them rather than reconstructing identical interval matrices.
  std::vector<SingularFrobeniusTowerTransferCertificate> tower_transfers;
  std::string tower_transfer_detail;
};

struct SingularFrobeniusPointTailCertificate {
  TailMajorantStatus status = TailMajorantStatus::Inconclusive;
  ErrorEnvelope value;
  std::string detail;
};

struct CertifiedSingularSeedEvaluation {
  LocalEvaluation evaluation;
  SingularFrobeniusPointTailCertificate tail;
};

struct CertifiedSingularOrdinaryBridge {
  TailMajorantStatus status = TailMajorantStatus::Inconclusive;
  LocalEvaluation evaluation;
  SingularFrobeniusPointTailCertificate singular_tail;
  RegularTaylorPointTailCertificate ordinary_tail;
  std::string seed_exact;
  std::string target_delta_exact;
  std::string singular_witness_radius_exact;
  std::string ordinary_witness_radius_exact;
  std::uint32_t singular_taylor_complete_max = 0;
  std::uint32_t ordinary_taylor_complete_max = 0;
  std::string detail;
};

namespace singular_tail_majorant_detail {

inline ExactEpsilonRational<Rational>
shifted_valuation(const ExactEpsilonRational<Rational> &value,
                  std::int32_t common_valuation) {
  if (value.zero)
    return value;
  auto shifted = value;
  shifted.valuation = local_detail::checked_i32(
      static_cast<std::int64_t>(value.valuation) - common_valuation,
      "singular tail normalized epsilon valuation");
  return shifted;
}

inline std::int32_t common_equation_valuation(
    const PreparedPhysicalClearedODE<Rational> &equation) {
  auto common = std::numeric_limits<std::int32_t>::max();
  for (const auto &value : equation.q_lags)
    if (!value.zero)
      common = std::min(common, value.valuation);
  for (const auto &lag : equation.c_lags)
    for (const auto &entry : lag)
      common = std::min(common, entry.value.valuation);
  if (common == std::numeric_limits<std::int32_t>::max())
    throw std::invalid_argument(
        "singular tail equation has no nonzero epsilon coefficient");
  return common;
}

inline std::uint32_t
polynomial_lag_degree(const PreparedPhysicalClearedODE<Rational> &equation) {
  std::size_t last = 0;
  for (std::size_t lag = 0; lag < equation.q_lags.size(); ++lag)
    if (!equation.q_lags[lag].zero)
      last = std::max(last, lag);
  for (std::size_t lag = 0; lag < equation.c_lags.size(); ++lag)
    if (!equation.c_lags[lag].empty())
      last = std::max(last, lag);
  if (last > std::numeric_limits<std::uint32_t>::max())
    throw std::overflow_error(
        "singular tail polynomial lag degree exceeds uint32");
  return static_cast<std::uint32_t>(last);
}

inline std::vector<Magnitude> rational_causal_c_lag_prefix_norm_upper(
    const std::vector<PhysicalODEMatrixEntry<Rational>> &entries,
    EpsilonWindow epsilon, std::uint32_t dimension,
    std::int32_t common_valuation) {
  const auto width = epsilon.width();
  std::vector<Magnitude> row_sums(width * static_cast<std::size_t>(dimension),
                                  Magnitude::zero());
  for (const auto &entry : entries) {
    const auto adjusted = shifted_valuation(entry.value, common_valuation);
    if (adjusted.valuation < 0)
      throw std::logic_error(
          "common epsilon monomial cancellation left a negative C valuation");
    const auto multiplier =
        physical_ode_detail::prepare_causal_multiplier(adjusted);
    const auto matrix = tail_majorant_detail::finite_causal_multiplier_matrix(
        multiplier, epsilon);
    for (std::size_t target_epsilon = 0; target_epsilon < width;
         ++target_epsilon)
      for (std::size_t source_epsilon = 0; source_epsilon <= target_epsilon;
           ++source_epsilon)
        row_sums[target_epsilon * dimension + entry.row] +=
            Magnitude::upper_abs(
                matrix[target_epsilon * width + source_epsilon]);
  }
  std::vector<Magnitude> prefix;
  prefix.reserve(width);
  auto maximum = Magnitude::zero();
  for (std::size_t target_epsilon = 0; target_epsilon < width;
       ++target_epsilon) {
    for (std::uint32_t row = 0; row < dimension; ++row)
      maximum = Magnitude::maximum(maximum,
                                   row_sums[target_epsilon * dimension + row]);
    prefix.push_back(maximum);
  }
  return prefix;
}

struct PreparedEntryMatrix {
  std::uint32_t row = 0;
  std::uint32_t column = 0;
  std::vector<ComplexBall> epsilon_matrix;
};

// Bound Q_k(A-kI)-C_k on every epsilon prefix, retaining cancellations
// between the indicial carrier A and C_k.  A acts on canonical tower slabs
// by
//
//   A = a I + b eps I + d/d(log),
//
// where d/d(log) maps log level p+1 to p under the stored log^p/p!
// convention.
inline std::vector<Magnitude> shifted_indicial_operator_prefix_norm_upper(
    const ExactEpsilonRational<Rational> *q,
    const std::vector<PhysicalODEMatrixEntry<Rational>> &c_entries,
    EpsilonWindow epsilon, std::uint32_t dimension, const Rational &a,
    const Rational &b, std::uint32_t log_complete_max, std::uint32_t lag,
    std::int32_t common_valuation) {
  const auto width = epsilon.width();
  std::vector<ComplexBall> q_matrix(width * width, ComplexBall(0));
  if (q != nullptr && !q->zero) {
    const auto adjusted = shifted_valuation(*q, common_valuation);
    if (adjusted.valuation < 0)
      throw std::logic_error(
          "common epsilon monomial cancellation left a negative q valuation");
    q_matrix = tail_majorant_detail::finite_causal_multiplier_matrix(
        physical_ode_detail::prepare_causal_multiplier(adjusted), epsilon);
  }
  std::vector<PreparedEntryMatrix> c_matrices;
  c_matrices.reserve(c_entries.size());
  for (const auto &entry : c_entries) {
    const auto adjusted = shifted_valuation(entry.value, common_valuation);
    if (adjusted.valuation < 0)
      throw std::logic_error(
          "common epsilon monomial cancellation left a negative C valuation");
    c_matrices.push_back(PreparedEntryMatrix{
        entry.row, entry.column,
        tail_majorant_detail::finite_causal_multiplier_matrix(
            physical_ode_detail::prepare_causal_multiplier(adjusted),
            epsilon)});
  }

  const auto log_width = static_cast<std::size_t>(log_complete_max) + 1;
  const auto vector_dimension =
      width * log_width * static_cast<std::size_t>(dimension);
  if (dimension == 0 || vector_dimension == 0)
    throw std::invalid_argument(
        "singular indicial norm received a zero dimension");
  std::vector<ComplexBall> row(vector_dimension, ComplexBall(0));
  const auto source_index = [width, dimension](std::size_t epsilon_index,
                                               std::uint32_t log_power,
                                               std::uint32_t component) {
    return ((static_cast<std::size_t>(log_power) * width + epsilon_index) *
            dimension) +
           component;
  };
  const auto a_minus_lag =
      local_detail::to_ball(a - Rational(std::to_string(lag)));
  const auto b_ball = local_detail::to_ball(b);
  std::vector<Magnitude> prefix;
  prefix.reserve(width);
  auto maximum = Magnitude::zero();
  for (std::size_t target_epsilon = 0; target_epsilon < width;
       ++target_epsilon) {
    for (std::uint32_t target_log = 0; target_log <= log_complete_max;
         ++target_log)
      for (std::uint32_t target_component = 0; target_component < dimension;
           ++target_component) {
        std::fill(row.begin(), row.end(), ComplexBall(0));
        for (std::size_t source_epsilon = 0; source_epsilon <= target_epsilon;
             ++source_epsilon) {
          const auto q_weight =
              q_matrix[target_epsilon * width + source_epsilon];
          row[source_index(source_epsilon, target_log, target_component)] +=
              q_weight * a_minus_lag;
          if (target_log < log_complete_max)
            row[source_index(source_epsilon, target_log + 1,
                             target_component)] += q_weight;
        }
        // Multiplication by b*eps happens before the causal Q_k
        // multiplier.  If u is the intermediate epsilon row, its source is
        // u-1.
        for (std::size_t intermediate_epsilon = 1;
             intermediate_epsilon <= target_epsilon; ++intermediate_epsilon)
          row[source_index(intermediate_epsilon - 1, target_log,
                           target_component)] +=
              q_matrix[target_epsilon * width + intermediate_epsilon] * b_ball;
        for (const auto &entry : c_matrices) {
          if (entry.row != target_component)
            continue;
          for (std::size_t source_epsilon = 0; source_epsilon <= target_epsilon;
               ++source_epsilon)
            row[source_index(source_epsilon, target_log, entry.column)] -=
                entry.epsilon_matrix[target_epsilon * width + source_epsilon];
        }
        auto row_sum = Magnitude::zero();
        for (const auto &coefficient : row)
          row_sum += Magnitude::upper_abs(coefficient);
        maximum = Magnitude::maximum(maximum, row_sum);
      }
    prefix.push_back(maximum);
  }
  return prefix;
}

struct StructuredFutureInverseCertificate {
  TailMajorantStatus status = TailMajorantStatus::Inconclusive;
  std::vector<Magnitude> normalized_inverse_n_times_prefix_upper;
  std::vector<Magnitude> inverse_n_times_prefix_upper;
  std::uint64_t exact_resolvent_checks = 0;
  std::uint64_t asymptotic_start = 0;
  std::string detail;
};

inline Rational exact_abs(const Rational &value) {
  return value.sign() < 0 ? -value : value;
}

inline Rational exact_matrix_infinity_norm(const std::vector<Rational> &matrix,
                                           std::size_t dimension) {
  if (dimension == 0 ||
      dimension > std::numeric_limits<std::size_t>::max() / dimension ||
      matrix.size() != dimension * dimension)
    throw std::invalid_argument(
        "exact matrix infinity norm received an invalid square matrix");
  Rational maximum(0);
  for (std::size_t row = 0; row < dimension; ++row) {
    Rational sum(0);
    for (std::size_t column = 0; column < dimension; ++column)
      sum += exact_abs(matrix[row * dimension + column]);
    if (maximum < sum)
      maximum = std::move(sum);
  }
  return maximum;
}

inline std::optional<Rational>
exact_inverse_infinity_norm(std::vector<Rational> matrix,
                            std::size_t dimension) {
  if (dimension == 0 ||
      dimension > std::numeric_limits<std::size_t>::max() / dimension ||
      matrix.size() != dimension * dimension)
    throw std::invalid_argument(
        "exact inverse norm received an invalid square matrix");
  std::vector<Rational> inverse(dimension * dimension, Rational(0));
  for (std::size_t diagonal = 0; diagonal < dimension; ++diagonal)
    inverse[diagonal * dimension + diagonal] = Rational(1);
  for (std::size_t column = 0; column < dimension; ++column) {
    auto pivot = column;
    while (pivot < dimension && matrix[pivot * dimension + column].is_zero())
      ++pivot;
    if (pivot == dimension)
      return std::nullopt;
    if (pivot != column)
      for (std::size_t entry = 0; entry < dimension; ++entry) {
        std::swap(matrix[column * dimension + entry],
                  matrix[pivot * dimension + entry]);
        std::swap(inverse[column * dimension + entry],
                  inverse[pivot * dimension + entry]);
      }
    const auto pivot_value = matrix[column * dimension + column];
    for (std::size_t entry = 0; entry < dimension; ++entry) {
      matrix[column * dimension + entry] =
          matrix[column * dimension + entry] / pivot_value;
      inverse[column * dimension + entry] =
          inverse[column * dimension + entry] / pivot_value;
    }
    for (std::size_t row = 0; row < dimension; ++row) {
      if (row == column)
        continue;
      const auto factor = matrix[row * dimension + column];
      if (factor.is_zero())
        continue;
      for (std::size_t entry = 0; entry < dimension; ++entry) {
        matrix[row * dimension + entry] -=
            factor * matrix[column * dimension + entry];
        inverse[row * dimension + entry] -=
            factor * inverse[column * dimension + entry];
      }
    }
  }
  return exact_matrix_infinity_norm(inverse, dimension);
}

inline std::vector<ComplexBall> finite_causal_unit_inverse_matrix(
    const physical_ode_detail::PreparedCausalEpsilonMultiplier &q0,
    EpsilonWindow epsilon) {
  const auto width = epsilon.width();
  if (width > std::numeric_limits<std::size_t>::max() / width)
    throw std::overflow_error(
        "structured singular inverse epsilon matrix size overflows");
  std::vector<ComplexBall> inverse(width * width, ComplexBall(0));
  for (std::size_t column = 0; column < width; ++column) {
    auto rhs = physical_ode_detail::zero_epsilon_vector(epsilon, 1);
    rhs.coefficients[column] = ComplexBall(1);
    const auto solved = physical_ode_detail::solve_formal_unit_q0(q0, rhs, 1);
    for (std::size_t row = 0; row < width; ++row)
      inverse[row * width + column] = solved.coefficients[row];
  }
  return inverse;
}

// Prove and bound every future Frobenius diagonal without confusing
// nilpotent epsilon/log couplings with spectral resonances.  On each fixed
// epsilon/log layer the exact diagonal block is
//
//        (n+a) I - q00^-1 C00.
//
// We invert that small Rational block for the finite range in which a
// resonance can occur.  Once its exact infinity norm is below n, a Neumann
// estimate proves the entire remaining half-line.  Positive-epsilon and log
// derivative terms are strictly triangular in the ordering
//
//        epsilon ascending, log power descending,
//
// so finite forward substitution gives a uniform inverse bound without any
// small-norm assumption on those couplings.
inline StructuredFutureInverseCertificate certify_structured_future_inverse(
    const PreparedPhysicalClearedODE<Rational> &equation, EpsilonWindow epsilon,
    const Rational &a, const Rational &b, std::uint32_t log_complete_max,
    std::int32_t common_valuation, std::uint64_t n0) {
  StructuredFutureInverseCertificate result;
  const auto dimension = static_cast<std::size_t>(equation.dimension);
  const auto width = epsilon.width();
  const auto log_width = static_cast<std::size_t>(log_complete_max) + 1;
  if (dimension == 0 || width == 0 || log_width == 0 || n0 == 0)
    throw std::invalid_argument(
        "structured singular inverse received a zero dimension");
  if (n0 > std::numeric_limits<ulong>::max())
    throw std::overflow_error(
        "structured singular inverse Taylor index exceeds FLINT ulong");

  const auto adjusted_q0 =
      shifted_valuation(equation.q_lags.front(), common_valuation);
  if (adjusted_q0.zero || adjusted_q0.valuation != 0)
    throw std::invalid_argument(
        "structured singular inverse requires a causal formal-unit q0");
  const auto q00 =
      adjusted_q0.numerator.front() / adjusted_q0.denominator.front();
  if (q00.is_zero())
    throw std::logic_error(
        "structured singular inverse lost the exact q0 leading unit");

  std::vector<Rational> residue(dimension * dimension, Rational(0));
  const auto &c0 = equation.c_lags.front();
  for (const auto &entry : c0) {
    const auto adjusted = shifted_valuation(entry.value, common_valuation);
    if (adjusted.valuation != 0)
      continue;
    const auto leading =
        adjusted.numerator.front() / adjusted.denominator.front();
    residue[static_cast<std::size_t>(entry.row) * dimension + entry.column] +=
        leading / q00;
  }
  std::vector<Rational> base_operator = residue;
  for (auto &entry : base_operator)
    entry = -entry;
  for (std::size_t diagonal = 0; diagonal < dimension; ++diagonal)
    base_operator[diagonal * dimension + diagonal] += a;

  const auto base_norm_exact =
      exact_matrix_infinity_norm(base_operator, dimension);
  const auto base_norm_upper =
      tail_majorant_detail::exact_rational_upper(base_norm_exact);
  constexpr std::uint64_t kMaximumFiniteResolventChecks = 4096;
  auto asymptotic_start = n0;
  for (;;) {
    if (asymptotic_start > std::numeric_limits<ulong>::max())
      throw std::overflow_error(
          "structured singular inverse asymptotic index exceeds FLINT ulong");
    const auto n_magnitude =
        Magnitude::from_ui(static_cast<ulong>(asymptotic_start));
    const auto gap = Magnitude::positive_difference_lower(
        Magnitude::one(), base_norm_upper / n_magnitude);
    if (!gap.is_zero())
      break;
    if (asymptotic_start - n0 >= kMaximumFiniteResolventChecks) {
      result.detail = "structured Frobenius finite resolvent range exceeds " +
                      std::to_string(kMaximumFiniteResolventChecks) +
                      " Taylor layers";
      return result;
    }
    ++asymptotic_start;
  }

  const auto q0_multiplier =
      physical_ode_detail::prepare_causal_multiplier(adjusted_q0);
  const auto q0_inverse =
      finite_causal_unit_inverse_matrix(q0_multiplier, epsilon);
  const auto q0_inverse_prefix =
      tail_majorant_detail::finite_causal_matrix_prefix_norm_upper(q0_inverse,
                                                                   width);

  struct PreparedC0Entry {
    std::uint32_t row = 0;
    std::uint32_t column = 0;
    std::vector<ComplexBall> matrix;
  };
  std::vector<PreparedC0Entry> prepared_c0;
  prepared_c0.reserve(c0.size());
  for (const auto &entry : c0) {
    const auto adjusted = shifted_valuation(entry.value, common_valuation);
    prepared_c0.push_back(PreparedC0Entry{
        entry.row, entry.column,
        tail_majorant_detail::finite_causal_multiplier_matrix(
            physical_ode_detail::prepare_causal_multiplier(adjusted),
            epsilon)});
  }

  std::vector<Magnitude> epsilon_block_norms(width * width, Magnitude::zero());
  const auto b_ball = local_detail::to_ball(b);
  for (std::size_t target_epsilon = 1; target_epsilon < width;
       ++target_epsilon) {
    for (std::size_t source_epsilon = 0; source_epsilon < target_epsilon;
         ++source_epsilon) {
      std::vector<ComplexBall> block(dimension * dimension, ComplexBall(0));
      for (const auto &entry : prepared_c0) {
        ComplexBall value(0);
        for (std::size_t intermediate = source_epsilon;
             intermediate <= target_epsilon; ++intermediate)
          value += q0_inverse[target_epsilon * width + intermediate] *
                   entry.matrix[intermediate * width + source_epsilon];
        block[static_cast<std::size_t>(entry.row) * dimension + entry.column] -=
            value;
      }
      if (target_epsilon == source_epsilon + 1)
        for (std::size_t diagonal = 0; diagonal < dimension; ++diagonal)
          block[diagonal * dimension + diagonal] += b_ball;
      epsilon_block_norms[target_epsilon * width + source_epsilon] =
          tail_majorant_detail::finite_matrix_infinity_norm_upper(block,
                                                                  dimension);
    }
  }

  struct Predecessor {
    std::size_t block = 0;
    Magnitude norm = Magnitude::zero();
  };
  const auto block_count = width * log_width;
  std::vector<std::vector<Predecessor>> predecessors(block_count);
  const auto block_index = [log_width](std::size_t epsilon_index,
                                       std::size_t descending_log_index) {
    return epsilon_index * log_width + descending_log_index;
  };
  for (std::size_t target_epsilon = 0; target_epsilon < width;
       ++target_epsilon) {
    for (std::size_t log_index = 0; log_index < log_width; ++log_index) {
      auto &incoming = predecessors[block_index(target_epsilon, log_index)];
      for (std::size_t source_epsilon = 0; source_epsilon < target_epsilon;
           ++source_epsilon) {
        const auto &norm =
            epsilon_block_norms[target_epsilon * width + source_epsilon];
        if (!norm.is_zero())
          incoming.push_back(
              Predecessor{block_index(source_epsilon, log_index), norm});
      }
      if (log_index != 0)
        incoming.push_back(Predecessor{
            block_index(target_epsilon, log_index - 1), Magnitude::one()});
    }
  }

  result.normalized_inverse_n_times_prefix_upper.assign(width,
                                                        Magnitude::zero());
  result.inverse_n_times_prefix_upper.assign(width, Magnitude::zero());
  const auto accumulate_bounds =
      [&](std::vector<Magnitude> &output,
          const std::vector<Magnitude> &rhs_prefix,
          const Magnitude &base_inverse_norm,
          const Magnitude &base_n_times_inverse_norm) {
        std::vector<Magnitude> z(block_count, Magnitude::zero());
        std::vector<Magnitude> prefix(width, Magnitude::zero());
        auto cumulative = Magnitude::zero();
        for (std::size_t target = 0; target < block_count; ++target) {
          auto propagated = Magnitude::zero();
          for (const auto &incoming : predecessors[target])
            propagated += incoming.norm * z[incoming.block];
          const auto epsilon_index = target / log_width;
          z[target] = base_n_times_inverse_norm * rhs_prefix[epsilon_index] +
                      base_inverse_norm * propagated;
          cumulative = Magnitude::maximum(cumulative, z[target]);
          if (target % log_width == log_width - 1)
            prefix[epsilon_index] = cumulative;
        }
        for (std::size_t index = 0; index < width; ++index)
          output[index] = Magnitude::maximum(output[index], prefix[index]);
      };
  const std::vector<Magnitude> identity_rhs_prefix(width, Magnitude::one());

  for (auto n = n0; n < asymptotic_start; ++n) {
    auto diagonal = base_operator;
    const Rational n_exact(std::to_string(n));
    for (std::size_t component = 0; component < dimension; ++component)
      diagonal[component * dimension + component] += n_exact;
    const auto inverse_norm =
        exact_inverse_infinity_norm(std::move(diagonal), dimension);
    if (!inverse_norm.has_value()) {
      result.detail = "exact future Frobenius resonance at Taylor index " +
                      std::to_string(n);
      return result;
    }
    const auto inverse_norm_upper =
        tail_majorant_detail::exact_rational_upper(*inverse_norm);
    const auto n_magnitude = Magnitude::from_ui(static_cast<ulong>(n));
    accumulate_bounds(result.normalized_inverse_n_times_prefix_upper,
                      identity_rhs_prefix, inverse_norm_upper,
                      n_magnitude * inverse_norm_upper);
    accumulate_bounds(result.inverse_n_times_prefix_upper, q0_inverse_prefix,
                      inverse_norm_upper, n_magnitude * inverse_norm_upper);
    ++result.exact_resolvent_checks;
  }

  const auto asymptotic_magnitude =
      Magnitude::from_ui(static_cast<ulong>(asymptotic_start));
  const auto asymptotic_gap = Magnitude::positive_difference_lower(
      Magnitude::one(), base_norm_upper / asymptotic_magnitude);
  if (asymptotic_gap.is_zero())
    throw std::logic_error(
        "structured singular inverse lost its asymptotic Neumann gap");
  const auto asymptotic_base_n_times_inverse =
      asymptotic_gap.reciprocal_upper();
  const auto asymptotic_base_inverse =
      asymptotic_base_n_times_inverse / asymptotic_magnitude;
  accumulate_bounds(result.normalized_inverse_n_times_prefix_upper,
                    identity_rhs_prefix, asymptotic_base_inverse,
                    asymptotic_base_n_times_inverse);
  accumulate_bounds(result.inverse_n_times_prefix_upper, q0_inverse_prefix,
                    asymptotic_base_inverse, asymptotic_base_n_times_inverse);
  result.asymptotic_start = asymptotic_start;
  result.status = TailMajorantStatus::Certified;
  result.detail = "exact finite base resolvents plus causal epsilon/log "
                  "back-substitution; exact_checks=" +
                  std::to_string(result.exact_resolvent_checks) +
                  "; asymptotic_start=" + std::to_string(asymptotic_start);
  return result;
}

struct NormalizedRecurrenceLagPrefixBounds {
  std::vector<Magnitude> scalar_quotient_upper;
  std::vector<Magnitude> correction_upper;
};

// Let S_k=Q_0^-1 Q_k, R_k=Q_0^-1 C_k and
//
//   D_n = Q_0 (n I + G),       G=A-R_0.
//
// Causal scalar epsilon multipliers commute with G, hence the exact lag-k
// recurrence transfer satisfies
//
//   D_n^-1 [Q_k(n-k+A)-C_k]
//     = S_k + (n I+G)^-1[-k S_k-R_k+S_k R_0].
//
// Taking this identity before norms removes the arbitrary scalar clearing
// polynomial from the dominant n term.
inline NormalizedRecurrenceLagPrefixBounds
normalized_recurrence_lag_prefix_bounds(
    const PreparedPhysicalClearedODE<Rational> &equation, EpsilonWindow epsilon,
    std::uint32_t lag, std::int32_t common_valuation) {
  if (lag == 0)
    throw std::invalid_argument(
        "normalized recurrence lag bound requires a positive lag");
  const auto dimension = static_cast<std::size_t>(equation.dimension);
  const auto width = epsilon.width();
  const auto adjusted_q0 =
      shifted_valuation(equation.q_lags.front(), common_valuation);
  const auto q0_inverse = finite_causal_unit_inverse_matrix(
      physical_ode_detail::prepare_causal_multiplier(adjusted_q0), epsilon);

  std::vector<ComplexBall> qk(width * width, ComplexBall(0));
  if (lag < equation.q_lags.size() && !equation.q_lags[lag].zero) {
    const auto adjusted =
        shifted_valuation(equation.q_lags[lag], common_valuation);
    qk = tail_majorant_detail::finite_causal_multiplier_matrix(
        physical_ode_detail::prepare_causal_multiplier(adjusted), epsilon);
  }
  std::vector<ComplexBall> scalar_quotient(width * width, ComplexBall(0));
  for (std::size_t target = 0; target < width; ++target)
    for (std::size_t source = 0; source <= target; ++source)
      for (std::size_t intermediate = source; intermediate <= target;
           ++intermediate)
        scalar_quotient[target * width + source] +=
            q0_inverse[target * width + intermediate] *
            qk[intermediate * width + source];

  const auto matrix_size = width * width * dimension * dimension;
  std::vector<ComplexBall> r0(matrix_size, ComplexBall(0));
  std::vector<ComplexBall> rk(matrix_size, ComplexBall(0));
  const auto flat_index = [width,
                           dimension](std::size_t target, std::size_t source,
                                      std::size_t row, std::size_t column) {
    return (((target * width + source) * dimension + row) * dimension) + column;
  };
  const auto accumulate_normalized_c =
      [&](const std::vector<PhysicalODEMatrixEntry<Rational>> &entries,
          std::vector<ComplexBall> &output) {
        for (const auto &entry : entries) {
          const auto adjusted =
              shifted_valuation(entry.value, common_valuation);
          const auto multiplier =
              tail_majorant_detail::finite_causal_multiplier_matrix(
                  physical_ode_detail::prepare_causal_multiplier(adjusted),
                  epsilon);
          for (std::size_t target = 0; target < width; ++target)
            for (std::size_t source = 0; source <= target; ++source)
              for (std::size_t intermediate = source; intermediate <= target;
                   ++intermediate)
                output[flat_index(target, source, entry.row, entry.column)] +=
                    q0_inverse[target * width + intermediate] *
                    multiplier[intermediate * width + source];
        }
      };
  accumulate_normalized_c(equation.c_lags.front(), r0);
  if (lag < equation.c_lags.size())
    accumulate_normalized_c(equation.c_lags[lag], rk);

  std::vector<ComplexBall> correction(matrix_size, ComplexBall(0));
  const auto lag_ball = ComplexBall(static_cast<long>(lag));
  for (std::size_t target = 0; target < width; ++target)
    for (std::size_t source = 0; source <= target; ++source)
      for (std::size_t row = 0; row < dimension; ++row)
        for (std::size_t column = 0; column < dimension; ++column) {
          auto &value = correction[flat_index(target, source, row, column)];
          value -= rk[flat_index(target, source, row, column)];
          for (std::size_t intermediate = source; intermediate <= target;
               ++intermediate)
            value += scalar_quotient[target * width + intermediate] *
                     r0[flat_index(intermediate, source, row, column)];
          if (row == column)
            value -= lag_ball * scalar_quotient[target * width + source];
        }

  NormalizedRecurrenceLagPrefixBounds result;
  result.scalar_quotient_upper =
      tail_majorant_detail::finite_causal_matrix_prefix_norm_upper(
          scalar_quotient, width);
  result.correction_upper.reserve(width);
  auto maximum = Magnitude::zero();
  for (std::size_t target = 0; target < width; ++target) {
    for (std::size_t row = 0; row < dimension; ++row) {
      auto row_sum = Magnitude::zero();
      for (std::size_t source = 0; source <= target; ++source)
        for (std::size_t column = 0; column < dimension; ++column)
          row_sum += Magnitude::upper_abs(
              correction[flat_index(target, source, row, column)]);
      maximum = Magnitude::maximum(maximum, row_sum);
    }
    result.correction_upper.push_back(maximum);
  }
  return result;
}

using BallMatrix = std::vector<ComplexBall>;
using BallMatrixSeries = std::vector<BallMatrix>;

inline BallMatrix multiply_ball_matrices(const BallMatrix &left,
                                         const BallMatrix &right,
                                         std::size_t dimension) {
  BallMatrix output(dimension * dimension, ComplexBall(0));
  for (std::size_t row = 0; row < dimension; ++row)
    for (std::size_t intermediate = 0; intermediate < dimension;
         ++intermediate) {
      const auto &weight = left[row * dimension + intermediate];
      if (weight.is_zero())
        continue;
      for (std::size_t column = 0; column < dimension; ++column)
        output[row * dimension + column] +=
            weight * right[intermediate * dimension + column];
    }
  return output;
}

inline std::optional<BallMatrix> invert_ball_matrix(BallMatrix matrix,
                                                    std::size_t dimension) {
  BallMatrix inverse(dimension * dimension, ComplexBall(0));
  for (std::size_t diagonal = 0; diagonal < dimension; ++diagonal)
    inverse[diagonal * dimension + diagonal] = ComplexBall(1);
  for (std::size_t column = 0; column < dimension; ++column) {
    auto pivot = column;
    while (pivot < dimension &&
           matrix[pivot * dimension + column].contains_zero())
      ++pivot;
    if (pivot == dimension)
      return std::nullopt;
    if (pivot != column)
      for (std::size_t entry = 0; entry < dimension; ++entry) {
        std::swap(matrix[column * dimension + entry],
                  matrix[pivot * dimension + entry]);
        std::swap(inverse[column * dimension + entry],
                  inverse[pivot * dimension + entry]);
      }
    const auto pivot_value = matrix[column * dimension + column];
    for (std::size_t entry = 0; entry < dimension; ++entry) {
      matrix[column * dimension + entry] =
          matrix[column * dimension + entry] / pivot_value;
      inverse[column * dimension + entry] =
          inverse[column * dimension + entry] / pivot_value;
    }
    for (std::size_t row = 0; row < dimension; ++row) {
      if (row == column)
        continue;
      const auto factor = matrix[row * dimension + column];
      if (factor.is_zero())
        continue;
      for (std::size_t entry = 0; entry < dimension; ++entry) {
        matrix[row * dimension + entry] -=
            factor * matrix[column * dimension + entry];
        inverse[row * dimension + entry] -=
            factor * inverse[column * dimension + entry];
      }
    }
  }
  return inverse;
}

inline BallMatrixSeries
multiply_ball_matrix_series(const BallMatrixSeries &left,
                            const BallMatrixSeries &right, std::size_t width,
                            std::size_t dimension) {
  BallMatrixSeries output(width,
                          BallMatrix(dimension * dimension, ComplexBall(0)));
  for (std::size_t degree = 0; degree < width; ++degree)
    for (std::size_t left_degree = 0; left_degree <= degree; ++left_degree) {
      const auto product = multiply_ball_matrices(
          left[left_degree], right[degree - left_degree], dimension);
      for (std::size_t entry = 0; entry < product.size(); ++entry)
        output[degree][entry] += product[entry];
    }
  return output;
}

inline std::optional<BallMatrixSeries>
invert_ball_matrix_series(const BallMatrixSeries &input, std::size_t width,
                          std::size_t dimension) {
  auto constant = invert_ball_matrix(input.front(), dimension);
  if (!constant.has_value())
    return std::nullopt;
  BallMatrixSeries inverse(width,
                           BallMatrix(dimension * dimension, ComplexBall(0)));
  inverse.front() = std::move(*constant);
  for (std::size_t degree = 1; degree < width; ++degree) {
    BallMatrix sum(dimension * dimension, ComplexBall(0));
    for (std::size_t input_degree = 1; input_degree <= degree; ++input_degree) {
      const auto product = multiply_ball_matrices(
          input[input_degree], inverse[degree - input_degree], dimension);
      for (std::size_t entry = 0; entry < product.size(); ++entry)
        sum[entry] += product[entry];
    }
    inverse[degree] = multiply_ball_matrices(inverse.front(), sum, dimension);
    for (auto &entry : inverse[degree])
      entry = -entry;
  }
  return inverse;
}

inline ComplexBall exact_real_interval_ball(const Rational &raw_left,
                                            const Rational &raw_right) {
  const auto left = raw_left < raw_right ? raw_left : raw_right;
  const auto right = raw_left < raw_right ? raw_right : raw_left;
  const auto midpoint = (left + right) / Rational(2);
  const auto radius = (right - left) / Rational(2);
  auto interval = ComplexBall::from_strings(midpoint.str());
  if (!radius.is_zero()) {
    const auto error = ComplexBall::from_strings(radius.str());
    arb_add_error(acb_realref(interval.raw()), acb_realref(error.raw()));
  }
  return interval;
}

// Enclose the exact transfer operator simultaneously for every real
// n>=n0 by evaluating it on x=1/n in the interval [0,1/n0].  The operator is
// a lower-triangular Toeplitz series in epsilon.  Its logarithmic ladder
// commutes with that series, so a finite nilpotent expansion retains all
// cancellations in
//
//       S_k + x (I+xG)^-1 F_k
//
// before absolute values are taken.
inline IntervalTransferLagDegreeBoundsResult
certify_interval_transfer_lag_degree_bounds_on_x_interval(
    const PreparedPhysicalClearedODE<Rational> &equation, EpsilonWindow epsilon,
    const Rational &a, const Rational &b, std::uint32_t log_complete_max,
    std::int32_t common_valuation, const Rational &x_left,
    const Rational &x_right) {
  IntervalTransferLagDegreeBoundsResult result;
  const auto width = epsilon.width();
  const auto dimension = static_cast<std::size_t>(equation.dimension);
  result.epsilon_width = width;
  result.physical_dimension = dimension;
  const auto x = exact_real_interval_ball(x_left, x_right);

  const auto adjusted_q0 =
      shifted_valuation(equation.q_lags.front(), common_valuation);
  const auto q0_inverse = finite_causal_unit_inverse_matrix(
      physical_ode_detail::prepare_causal_multiplier(adjusted_q0), epsilon);
  const auto normalized_c_series =
      [&](const std::vector<PhysicalODEMatrixEntry<Rational>> &entries) {
        BallMatrixSeries output(
            width, BallMatrix(dimension * dimension, ComplexBall(0)));
        for (const auto &entry : entries) {
          const auto adjusted =
              shifted_valuation(entry.value, common_valuation);
          const auto multiplier =
              tail_majorant_detail::finite_causal_multiplier_matrix(
                  physical_ode_detail::prepare_causal_multiplier(adjusted),
                  epsilon);
          for (std::size_t degree = 0; degree < width; ++degree)
            for (std::size_t intermediate = 0; intermediate <= degree;
                 ++intermediate)
              output[degree][static_cast<std::size_t>(entry.row) * dimension +
                             entry.column] +=
                  q0_inverse[degree * width + intermediate] *
                  multiplier[intermediate * width];
        }
        return output;
      };
  const auto r0 = normalized_c_series(equation.c_lags.front());

  BallMatrixSeries g = r0;
  for (auto &coefficient : g)
    for (auto &entry : coefficient)
      entry = -entry;
  for (std::size_t diagonal = 0; diagonal < dimension; ++diagonal)
    g.front()[diagonal * dimension + diagonal] += local_detail::to_ball(a);
  if (width > 1)
    for (std::size_t diagonal = 0; diagonal < dimension; ++diagonal)
      g[1][diagonal * dimension + diagonal] += local_detail::to_ball(b);

  BallMatrixSeries scaled_diagonal(
      width, BallMatrix(dimension * dimension, ComplexBall(0)));
  for (std::size_t degree = 0; degree < width; ++degree)
    for (std::size_t entry = 0; entry < dimension * dimension; ++entry)
      scaled_diagonal[degree][entry] = x * g[degree][entry];
  for (std::size_t diagonal = 0; diagonal < dimension; ++diagonal)
    scaled_diagonal.front()[diagonal * dimension + diagonal] += ComplexBall(1);
  auto inverse = invert_ball_matrix_series(scaled_diagonal, width, dimension);
  if (!inverse.has_value()) {
    result.detail =
        "x-interval base inversion did not separate every future diagonal";
    return result;
  }

  const auto degree = polynomial_lag_degree(equation);
  std::vector<std::vector<std::vector<Magnitude>>> lag_degree_row_upper(
      static_cast<std::size_t>(degree) + 1);
  for (std::uint32_t lag = 1; lag <= degree; ++lag) {
    BallMatrixSeries qk_series(
        width, BallMatrix(dimension * dimension, ComplexBall(0)));
    if (lag < equation.q_lags.size() && !equation.q_lags[lag].zero) {
      const auto adjusted =
          shifted_valuation(equation.q_lags[lag], common_valuation);
      const auto qk = tail_majorant_detail::finite_causal_multiplier_matrix(
          physical_ode_detail::prepare_causal_multiplier(adjusted), epsilon);
      std::vector<ComplexBall> scalar(width, ComplexBall(0));
      for (std::size_t output_degree = 0; output_degree < width;
           ++output_degree)
        for (std::size_t intermediate = 0; intermediate <= output_degree;
             ++intermediate)
          scalar[output_degree] +=
              q0_inverse[output_degree * width + intermediate] *
              qk[intermediate * width];
      for (std::size_t epsilon_degree = 0; epsilon_degree < width;
           ++epsilon_degree)
        for (std::size_t diagonal = 0; diagonal < dimension; ++diagonal)
          qk_series[epsilon_degree][diagonal * dimension + diagonal] =
              scalar[epsilon_degree];
    }
    const auto rk =
        lag < equation.c_lags.size()
            ? normalized_c_series(equation.c_lags[lag])
            : BallMatrixSeries(
                  width, BallMatrix(dimension * dimension, ComplexBall(0)));
    auto correction =
        multiply_ball_matrix_series(qk_series, r0, width, dimension);
    const auto lag_ball = ComplexBall(static_cast<long>(lag));
    for (std::size_t epsilon_degree = 0; epsilon_degree < width;
         ++epsilon_degree)
      for (std::size_t entry = 0; entry < dimension * dimension; ++entry) {
        correction[epsilon_degree][entry] -= rk[epsilon_degree][entry];
        correction[epsilon_degree][entry] -=
            lag_ball * qk_series[epsilon_degree][entry];
      }

    auto inverse_power = *inverse;
    auto x_log_power = x;
    std::vector<std::vector<Magnitude>> lag_rows(
        width, std::vector<Magnitude>(dimension, Magnitude::zero()));
    for (std::uint32_t log_shift = 0; log_shift <= log_complete_max;
         ++log_shift) {
      auto transfer = multiply_ball_matrix_series(inverse_power, correction,
                                                  width, dimension);
      for (auto &coefficient : transfer)
        for (auto &entry : coefficient)
          entry *= x_log_power;
      if (log_shift == 0)
        for (std::size_t epsilon_degree = 0; epsilon_degree < width;
             ++epsilon_degree)
          for (std::size_t entry = 0; entry < dimension * dimension; ++entry)
            transfer[epsilon_degree][entry] += qk_series[epsilon_degree][entry];
      for (std::size_t epsilon_degree = 0; epsilon_degree < width;
           ++epsilon_degree)
        for (std::size_t row = 0; row < dimension; ++row) {
          auto row_sum = Magnitude::zero();
          for (std::size_t column = 0; column < dimension; ++column)
            row_sum += Magnitude::upper_abs(
                transfer[epsilon_degree][row * dimension + column]);
          lag_rows[epsilon_degree][row] += row_sum;
        }
      if (log_shift != log_complete_max) {
        inverse_power = multiply_ball_matrix_series(inverse_power, *inverse,
                                                    width, dimension);
        x_log_power *= -x;
      }
    }
    lag_degree_row_upper[lag] = std::move(lag_rows);
  }

  result.lag_degree_row_upper = std::move(lag_degree_row_upper);
  result.status = TailMajorantStatus::Certified;
  result.detail =
      "all-future x=1/n interval transfer with causal epsilon Toeplitz and "
      "finite log-nilpotent expansion";
  return result;
}

inline IntervalTransferDegreeBoundsResult
radius_weighted_interval_transfer_degree_bounds(
    const IntervalTransferLagDegreeBoundsResult &lag_bounds,
    const Rational &radius) {
  IntervalTransferDegreeBoundsResult result;
  result.status = lag_bounds.status;
  result.detail = lag_bounds.detail;
  if (lag_bounds.status != TailMajorantStatus::Certified)
    return result;
  if (lag_bounds.epsilon_width == 0 ||
      lag_bounds.physical_dimension == 0 ||
      lag_bounds.lag_degree_row_upper.empty()) {
    result.detail = "interval transfer lag enclosure lost its shape";
    result.status = TailMajorantStatus::Inconclusive;
    return result;
  }
  if (lag_bounds.lag_degree_row_upper.size() == 1) {
    result.degree_row_upper.assign(
        lag_bounds.epsilon_width,
        std::vector<Magnitude>(lag_bounds.physical_dimension,
                               Magnitude::zero()));
    return result;
  }
  const auto first_material = std::find_if(
      std::next(lag_bounds.lag_degree_row_upper.begin()),
      lag_bounds.lag_degree_row_upper.end(),
      [](const auto &lag) { return !lag.empty(); });
  if (first_material == lag_bounds.lag_degree_row_upper.end() ||
      first_material->front().empty()) {
    result.detail = "interval transfer lag enclosure has no material rows";
    result.status = TailMajorantStatus::Inconclusive;
    return result;
  }
  const auto width = lag_bounds.epsilon_width;
  const auto dimension = lag_bounds.physical_dimension;
  if (first_material->size() != width ||
      first_material->front().size() != dimension)
    throw std::logic_error(
        "interval transfer lag enclosure disagrees with its retained shape");
  result.degree_row_upper.assign(
      width, std::vector<Magnitude>(dimension, Magnitude::zero()));
  const auto radius_upper =
      tail_majorant_detail::exact_rational_upper(radius);
  for (std::size_t lag = 1;
       lag < lag_bounds.lag_degree_row_upper.size(); ++lag) {
    const auto &lag_rows = lag_bounds.lag_degree_row_upper[lag];
    if (lag_rows.empty())
      continue;
    if (lag_rows.size() != width)
      throw std::logic_error(
          "interval transfer lag enclosure changed epsilon width");
    const auto radius_factor =
        radius_upper.power_upper(static_cast<ulong>(lag));
    for (std::size_t epsilon_degree = 0; epsilon_degree < width;
         ++epsilon_degree) {
      if (lag_rows[epsilon_degree].size() != dimension)
        throw std::logic_error(
            "interval transfer lag enclosure changed physical dimension");
      for (std::size_t row = 0; row < dimension; ++row)
        result.degree_row_upper[epsilon_degree][row] +=
            radius_factor * lag_rows[epsilon_degree][row];
    }
  }
  return result;
}

inline IntervalTransferDegreeBoundsResult
certify_interval_transfer_degree_bounds_on_x_interval(
    const PreparedPhysicalClearedODE<Rational> &equation, EpsilonWindow epsilon,
    const Rational &a, const Rational &b, std::uint32_t log_complete_max,
    std::int32_t common_valuation, const Rational &x_left,
    const Rational &x_right, const Rational &radius) {
  return radius_weighted_interval_transfer_degree_bounds(
      certify_interval_transfer_lag_degree_bounds_on_x_interval(
          equation, epsilon, a, b, log_complete_max, common_valuation, x_left,
          x_right),
      radius);
}

inline IntervalTransferContractionResult weighted_interval_transfer_contraction(
    const IntervalTransferDegreeBoundsResult &degree_bounds,
    ulong epsilon_weight_base) {
  if (epsilon_weight_base == 0)
    throw std::invalid_argument(
        "interval transfer epsilon weight base must be positive");
  IntervalTransferContractionResult result;
  result.status = degree_bounds.status;
  result.epsilon_weight_base = epsilon_weight_base;
  result.detail = degree_bounds.detail;
  if (degree_bounds.status != TailMajorantStatus::Certified)
    return result;
  if (degree_bounds.degree_row_upper.empty()) {
    result.detail = "interval transfer degree enclosure is empty";
    result.status = TailMajorantStatus::Inconclusive;
    return result;
  }
  const auto dimension = degree_bounds.degree_row_upper.front().size();
  const auto epsilon_weight = Magnitude::from_ui(epsilon_weight_base);
  auto maximum = Magnitude::zero();
  result.contraction_upper.reserve(degree_bounds.degree_row_upper.size());
  for (std::size_t target_epsilon = 0;
       target_epsilon < degree_bounds.degree_row_upper.size();
       ++target_epsilon) {
    for (std::size_t row = 0; row < dimension; ++row) {
      auto row_sum = Magnitude::zero();
      for (std::size_t epsilon_degree = 0; epsilon_degree <= target_epsilon;
           ++epsilon_degree)
        row_sum +=
            degree_bounds.degree_row_upper[epsilon_degree][row] /
            epsilon_weight.power_upper(static_cast<ulong>(epsilon_degree));
      maximum = Magnitude::maximum(maximum, row_sum);
    }
    result.contraction_upper.push_back(maximum);
  }
  return result;
}

inline IntervalTransferContractionResult
certify_interval_transfer_contraction_on_x_interval(
    const PreparedPhysicalClearedODE<Rational> &equation, EpsilonWindow epsilon,
    const Rational &a, const Rational &b, std::uint32_t log_complete_max,
    std::int32_t common_valuation, const Rational &x_left,
    const Rational &x_right, const Rational &radius,
    ulong epsilon_weight_base = 1) {
  return weighted_interval_transfer_contraction(
      certify_interval_transfer_degree_bounds_on_x_interval(
          equation, epsilon, a, b, log_complete_max, common_valuation, x_left,
          x_right, radius),
      epsilon_weight_base);
}

inline IntervalTransferContractionResult
certify_weighted_interval_transfer_contraction(
    const IntervalTransferDegreeBoundsResult &degree_bounds) {
  if (degree_bounds.status != TailMajorantStatus::Certified)
    return weighted_interval_transfer_contraction(degree_bounds, 1);
  IntervalTransferContractionResult last;
  for (const auto epsilon_weight_base :
       {1UL,       2UL,       4UL,       8UL,       16UL,      32UL,
        64UL,      128UL,     256UL,     512UL,     1024UL,    2048UL,
        4096UL,    8192UL,    16384UL,   32768UL,   65536UL,   131072UL,
        262144UL,  524288UL,  1048576UL}) {
    last = weighted_interval_transfer_contraction(degree_bounds,
                                                  epsilon_weight_base);
    if (last.status == TailMajorantStatus::Certified &&
        !last.contraction_upper.empty() &&
        !Magnitude::positive_difference_lower(Magnitude::one(),
                                              last.contraction_upper.back())
             .is_zero()) {
      last.detail +=
          "; epsilon_weight_base=" + std::to_string(epsilon_weight_base);
      return last;
    }
  }
  last.detail +=
      "; no tested geometric epsilon weight proves contraction below one";
  return last;
}

inline IntervalTransferContractionResult certify_interval_transfer_contraction(
    const PreparedPhysicalClearedODE<Rational> &equation, EpsilonWindow epsilon,
    const Rational &a, const Rational &b, std::uint32_t log_complete_max,
    std::int32_t common_valuation, std::uint64_t n0, const Rational &radius) {
  if (n0 == 0 || n0 > std::numeric_limits<ulong>::max())
    throw std::invalid_argument(
        "interval transfer contraction received an invalid Taylor index");
  const auto degree_bounds =
      certify_interval_transfer_degree_bounds_on_x_interval(
          equation, epsilon, a, b, log_complete_max, common_valuation,
          Rational(0), Rational(1) / Rational(std::to_string(n0)), radius);
  return certify_weighted_interval_transfer_contraction(degree_bounds);
}

inline Magnitude tower_coefficient_prefix_norm_upper(
    const scc_completeness_detail::FormalSlabs &slabs,
    const std::string &a_exact, const std::string &b_exact,
    std::uint32_t log_complete_max, EpsilonWindow model_epsilon,
    std::uint32_t taylor_index, std::size_t prefix_index,
    std::uint32_t dimension, ulong epsilon_weight_base = 1) {
  if (epsilon_weight_base == 0)
    throw std::invalid_argument(
        "tower coefficient epsilon weight base must be positive");
  const auto epsilon_weight = Magnitude::from_ui(epsilon_weight_base);
  auto maximum = Magnitude::zero();
  const auto prefix_power = local_detail::checked_i32(
      static_cast<std::int64_t>(model_epsilon.min_power) +
          static_cast<std::int64_t>(prefix_index),
      "singular tail prefix epsilon power");
  for (std::uint32_t log = 0; log <= log_complete_max; ++log) {
    const scc_completeness_detail::FormalTag tag{a_exact, b_exact, log};
    const auto found = slabs.find(tag);
    if (found == slabs.end())
      continue;
    if (taylor_index >= found->second.size())
      throw std::invalid_argument(
          "singular tail tower coefficient lies outside its formal slab");
    for (std::int32_t power = model_epsilon.min_power; power <= prefix_power;
         ++power)
      for (std::uint32_t component = 0; component < dimension; ++component)
        maximum = Magnitude::maximum(
            maximum, tail_majorant_detail::exact_rational_upper(
                         found->second[taylor_index].at(power, component)) /
                         epsilon_weight.power_upper(static_cast<ulong>(
                             power - model_epsilon.min_power)));
  }
  return maximum * epsilon_weight.power_upper(prefix_index);
}

inline std::string tower_key(const std::string &a, const std::string &b) {
  return a + "\x1f" + b;
}

struct FrobeniusTowerRecord {
  std::string a;
  std::string b;
  std::uint32_t max_log = 0;
};

inline std::map<std::string, FrobeniusTowerRecord>
collect_frobenius_tower_records(const LocalSolution<Rational> &solution) {
  const auto slabs = scc_completeness_detail::collect_formal_slabs(solution);
  std::map<std::string, FrobeniusTowerRecord> records;
  for (const auto &[tag, slab] : slabs) {
    (void)slab;
    const auto key = tower_key(tag.a, tag.b);
    auto [found, inserted] = records.try_emplace(
        key, FrobeniusTowerRecord{tag.a, tag.b, tag.log_power});
    if (!inserted)
      found->second.max_log =
          std::max(found->second.max_log, tag.log_power);
  }
  return records;
}

inline LocalSolution<Rational>
restrict_exact_local_epsilon(const LocalSolution<Rational> &source,
                             EpsilonWindow epsilon) {
  if (epsilon.min_power < source.epsilon.min_power ||
      epsilon.complete_max > source.epsilon.complete_max)
    throw std::invalid_argument(
        "restricted singular local epsilon window lies outside its source");
  auto restricted = source;
  restricted.epsilon = epsilon;
  restricted.error = {};
  for (std::size_t sector_index = 0; sector_index < source.sectors.size();
       ++sector_index) {
    auto &target_sector = restricted.sectors[sector_index];
    const auto &source_sector = source.sectors[sector_index];
    target_sector.coefficients.assign(restricted.sector_size(), Rational(0));
    for (std::int32_t power = epsilon.min_power; power <= epsilon.complete_max;
         ++power) {
      const auto source_epsilon_index =
          static_cast<std::size_t>(power - source.epsilon.min_power);
      const auto target_epsilon_index =
          static_cast<std::size_t>(power - epsilon.min_power);
      for (std::uint32_t taylor = 0; taylor <= source.taylor_complete_max;
           ++taylor)
        for (std::uint32_t component = 0; component < source.dimension;
             ++component)
          target_sector.coefficients[local_detail::sector_index(
              restricted, target_epsilon_index, taylor, component)] =
              source_sector.coefficients[local_detail::sector_index(
                  source, source_epsilon_index, taylor, component)];
    }
  }
  return restricted;
}

} // namespace singular_tail_majorant_detail

inline SingularFrobeniusPrefixCertificate
certify_singular_rational_shadow_prefix(
    const PreparedPhysicalClearedODE<Rational> &equation,
    const LocalSolution<Rational> &solution, EpsilonWindow claimed_epsilon,
    std::uint32_t retained_taylor_complete_max) {
  physical_ode_detail::validate_ode(equation);
  local_detail::validate_local_solution(solution, true);
  if (solution.dimension != equation.dimension ||
      solution.checkpoint_identity.empty() || equation.payload_identity.empty())
    throw std::invalid_argument(
        "singular prefix certificate lost its dimension or provenance binding");
  SingularFrobeniusPrefixCertificate result;
  result.epsilon = claimed_epsilon;
  result.taylor_complete_max = retained_taylor_complete_max;
  result.physical_payload_identity = equation.payload_identity;
  result.local_checkpoint_identity = solution.checkpoint_identity;
  result.exact_residual = certify_scc_parent_exact_formal_residual(
      equation, solution, claimed_epsilon, retained_taylor_complete_max);
  try {
    const auto common_valuation =
        singular_tail_majorant_detail::common_equation_valuation(equation);
    const auto records =
        singular_tail_majorant_detail::collect_frobenius_tower_records(
            solution);
    const auto n0 =
        static_cast<std::uint64_t>(retained_taylor_complete_max) + 1;
    if (n0 > std::numeric_limits<ulong>::max())
      throw std::overflow_error(
          "singular tail first unseen order exceeds FLINT ulong");
    result.tower_transfers.reserve(records.size());
    for (const auto &[key, record] : records) {
      (void)key;
      result.tower_transfers.push_back(
          SingularFrobeniusTowerTransferCertificate{
              record.a, record.b, record.max_log,
              singular_tail_majorant_detail::
                  certify_interval_transfer_lag_degree_bounds_on_x_interval(
                      equation, claimed_epsilon, Rational(record.a),
                      Rational(record.b), record.max_log, common_valuation,
                      Rational(0),
                      Rational(1) / Rational(std::to_string(n0)))});
    }
    result.tower_transfer_detail =
        "radius-independent interval transfers prepared for " +
        std::to_string(result.tower_transfers.size()) + " Frobenius towers";
  } catch (const std::exception &error) {
    result.tower_transfers.clear();
    result.tower_transfer_detail =
        std::string("radius-independent interval transfer cache unavailable: ") +
        error.what();
  }
  return result;
}

inline SingularFrobeniusTailModelResult
prepare_singular_rational_shadow_tail_model(
    const PreparedPhysicalClearedODE<Rational> &equation,
    const LocalSolution<Rational> &solution, EpsilonWindow claimed_epsilon,
    std::uint32_t retained_taylor_complete_max,
    const std::string &witness_radius_exact,
    const SingularFrobeniusPrefixCertificate *prefix_certificate = nullptr) {
  using namespace singular_tail_majorant_detail;
  const auto unsupported = [](std::string detail) {
    return SingularFrobeniusTailModelResult{TailMajorantStatus::Unsupported,
                                            std::nullopt, std::move(detail)};
  };
  const auto inconclusive = [](std::string detail) {
    return SingularFrobeniusTailModelResult{TailMajorantStatus::Inconclusive,
                                            std::nullopt, std::move(detail)};
  };

  physical_ode_detail::validate_ode(equation);
  local_detail::validate_local_solution(solution, true);
  if (solution.dimension != equation.dimension ||
      solution.checkpoint_identity.empty() || equation.payload_identity.empty())
    throw std::invalid_argument(
        "singular tail model lost its dimension or provenance binding");
  if (!solution.error.empty())
    return unsupported(
        "singular tail model cannot absorb an existing error envelope");
  if (claimed_epsilon.min_power < solution.epsilon.min_power ||
      claimed_epsilon.complete_max > solution.epsilon.complete_max ||
      retained_taylor_complete_max > solution.taylor_complete_max)
    throw std::invalid_argument(
        "singular tail claim lies outside the retained exact local");
  Rational radius(0);
  try {
    radius = Rational(witness_radius_exact);
  } catch (const std::invalid_argument &) {
    throw std::invalid_argument(
        "singular tail witness radius must be exact rational");
  }
  if (radius.sign() <= 0)
    throw std::invalid_argument(
        "singular tail witness radius must be positive");

  if (prefix_certificate != nullptr) {
    if (prefix_certificate->epsilon.min_power != claimed_epsilon.min_power ||
        prefix_certificate->epsilon.complete_max !=
            claimed_epsilon.complete_max ||
        prefix_certificate->taylor_complete_max !=
            retained_taylor_complete_max ||
        prefix_certificate->physical_payload_identity !=
            equation.payload_identity ||
        prefix_certificate->local_checkpoint_identity !=
            solution.checkpoint_identity ||
        prefix_certificate->exact_residual.epsilon.min_power !=
            claimed_epsilon.min_power ||
        prefix_certificate->exact_residual.epsilon.complete_max !=
            claimed_epsilon.complete_max ||
        prefix_certificate->exact_residual.taylor_complete_max !=
            retained_taylor_complete_max)
      throw std::invalid_argument("singular tail prefix certificate does not "
                                  "bind this equation/local/window");
  } else {
    try {
      (void)certify_scc_parent_exact_formal_residual(
          equation, solution, claimed_epsilon,
          retained_taylor_complete_max);
    } catch (const std::exception &error) {
      return inconclusive(std::string("exact Frobenius prefix does not satisfy "
                                      "its physical q/C equation: ") +
                          error.what());
    }
  }

  const auto common_valuation = common_equation_valuation(equation);
  const auto adjusted_q0 =
      shifted_valuation(equation.q_lags.front(), common_valuation);
  if (adjusted_q0.zero || adjusted_q0.valuation != 0)
    return unsupported("common epsilon monomial cancellation does not leave "
                       "q(0,eps) as a formal unit");
  for (const auto &value : equation.q_lags)
    if (!value.zero && shifted_valuation(value, common_valuation).valuation < 0)
      return unsupported(
          "common epsilon monomial cancellation leaves a negative q valuation");
  for (const auto &lag : equation.c_lags)
    for (const auto &entry : lag)
      if (shifted_valuation(entry.value, common_valuation).valuation < 0)
        return unsupported("common epsilon monomial cancellation leaves a "
                           "negative C valuation");

  const auto slabs = scc_completeness_detail::collect_formal_slabs(solution);
  const auto tower_records = collect_frobenius_tower_records(solution);
  if (tower_records.empty())
    return unsupported(
        "singular tail model has no canonical Rational Frobenius towers");

  SingularFrobeniusTailModel model;
  model.epsilon = claimed_epsilon;
  model.dimension = solution.dimension;
  model.taylor_complete_max = retained_taylor_complete_max;
  model.polynomial_lag_degree = polynomial_lag_degree(equation);
  model.witness_radius_exact = radius.str();
  model.chart = solution.chart;
  model.prescriptions = solution.prescriptions;
  model.physical_payload_identity = equation.payload_identity;
  model.local_checkpoint_identity = solution.checkpoint_identity;
  if (model.polynomial_lag_degree != 0 &&
      retained_taylor_complete_max + 1 < model.polynomial_lag_degree)
    return unsupported("retained Frobenius prefix is shorter than the physical "
                       "polynomial lag window");
  const auto radius_upper = tail_majorant_detail::exact_rational_upper(radius);
  const auto n0 = static_cast<std::uint64_t>(retained_taylor_complete_max) + 1;
  if (n0 > std::numeric_limits<ulong>::max())
    throw std::overflow_error(
        "singular tail first unseen order exceeds FLINT ulong");
  for (const auto &[key, record] : tower_records) {
    (void)key;
    const Rational a(record.a);
    const Rational b(record.b);
    const auto future_inverse = certify_structured_future_inverse(
        equation, claimed_epsilon, a, b, record.max_log, common_valuation, n0);
    if (future_inverse.status != TailMajorantStatus::Certified)
      return inconclusive(
          "structured all-future Frobenius inverse is inconclusive for "
          "tower (a=" +
          record.a + ",b=" + record.b + "): " + future_inverse.detail);
    IntervalTransferContractionResult interval_contraction;
    const SingularFrobeniusTowerTransferCertificate *cached_transfer =
        nullptr;
    if (prefix_certificate != nullptr) {
      const auto found = std::find_if(
          prefix_certificate->tower_transfers.begin(),
          prefix_certificate->tower_transfers.end(),
          [&](const auto &candidate) {
            return candidate.a_exact == record.a &&
                   candidate.b_exact == record.b &&
                   candidate.log_complete_max == record.max_log;
          });
      if (found != prefix_certificate->tower_transfers.end())
        cached_transfer = &*found;
    }
    if (cached_transfer != nullptr &&
        cached_transfer->interval_transfer.status ==
            TailMajorantStatus::Certified &&
        cached_transfer->interval_transfer.epsilon_width ==
            claimed_epsilon.width() &&
        cached_transfer->interval_transfer.physical_dimension ==
            static_cast<std::size_t>(equation.dimension) &&
        cached_transfer->interval_transfer.lag_degree_row_upper.size() ==
            static_cast<std::size_t>(model.polynomial_lag_degree) + 1) {
      interval_contraction =
          certify_weighted_interval_transfer_contraction(
              radius_weighted_interval_transfer_degree_bounds(
                  cached_transfer->interval_transfer, radius));
    } else {
      interval_contraction = certify_interval_transfer_contraction(
          equation, claimed_epsilon, a, b, record.max_log, common_valuation,
          n0, radius);
    }
    if (interval_contraction.status != TailMajorantStatus::Certified)
      return inconclusive(
          "structured all-future Frobenius contraction is inconclusive for "
          "tower (a=" +
          record.a + ",b=" + record.b + "): " + interval_contraction.detail);
    SingularFrobeniusTowerTailModel tower;
    tower.a_exact = record.a;
    tower.b_exact = record.b;
    tower.log_complete_max = record.max_log;
    tower.epsilon_weight_base = interval_contraction.epsilon_weight_base;
    tower.induction_m_upper.reserve(claimed_epsilon.width());
    tower.contraction_upper.reserve(claimed_epsilon.width());

    for (std::size_t prefix_index = 0; prefix_index < claimed_epsilon.width();
         ++prefix_index) {
      const auto rho = interval_contraction.contraction_upper.at(prefix_index);
      const auto contraction_gap =
          Magnitude::positive_difference_lower(Magnitude::one(), rho);
      if (contraction_gap.is_zero())
        return inconclusive(
            "Frobenius all-future contraction rho is not below one for tower "
            "(a=" +
            record.a + ",b=" + record.b + ") on epsilon prefix ending at " +
            std::to_string(
                static_cast<std::int64_t>(claimed_epsilon.min_power) +
                static_cast<std::int64_t>(prefix_index)) +
            "; rho_upper=" + std::to_string(rho.approximate_upper()));

      auto induction_m = Magnitude::zero();
      if (model.polynomial_lag_degree != 0) {
        const auto first =
            retained_taylor_complete_max + 1 - model.polynomial_lag_degree;
        for (std::uint32_t n = first; n <= retained_taylor_complete_max; ++n) {
          const auto coefficient = tower_coefficient_prefix_norm_upper(
              slabs, record.a, record.b, record.max_log, claimed_epsilon, n,
              prefix_index, solution.dimension,
              interval_contraction.epsilon_weight_base);
          induction_m = Magnitude::maximum(
              induction_m,
              coefficient * radius_upper.power_upper(static_cast<ulong>(n)));
        }
      }
      tower.induction_m_upper.push_back(induction_m);
      tower.contraction_upper.push_back(rho);
    }
    model.towers.push_back(std::move(tower));
  }

  model.provenance =
      "certified exact Rational Frobenius all-future recurrence contraction; "
      "canonical integer-shift towers; nested invariant causal epsilon "
      "prefixes; "
      "retained_taylor_complete_max=" +
      std::to_string(retained_taylor_complete_max) +
      "; witness_radius_exact=" + radius.str() +
      "; physical_payload_identity=" + model.physical_payload_identity +
      "; local_checkpoint_identity=" + model.local_checkpoint_identity;
  return {TailMajorantStatus::Certified, std::move(model),
          "singular Rational-shadow Frobenius tail model prepared"};
}

inline SingularFrobeniusPointTailCertificate
certify_singular_frobenius_point_tail(const SingularFrobeniusTailModel &model,
                                      const LocalSolution<Rational> &solution,
                                      const RealEvaluationPoint &input_point,
                                      EvaluationOptions options = {}) {
  using namespace singular_tail_majorant_detail;
  if (solution.checkpoint_identity != model.local_checkpoint_identity ||
      solution.dimension != model.dimension ||
      model.taylor_complete_max > solution.taylor_complete_max ||
      !tail_majorant_detail::same_chart_geometry(solution.chart, model.chart) ||
      !tail_majorant_detail::same_prescriptions(solution.prescriptions,
                                                model.prescriptions))
    throw std::invalid_argument(
        "singular tail model is not bound to this exact local solution");
  const auto point = line_integration_detail::require_exact_rational_point(
      input_point, "singular-frobenius-tail-evaluation");
  const Rational exact_point(point.exact_coordinate);
  const auto modulus = tail_majorant_detail::abs_rational(exact_point);
  const Rational radius(model.witness_radius_exact);
  SingularFrobeniusPointTailCertificate result;
  if (!(modulus < radius)) {
    result.detail = "Frobenius seed point is not strictly inside its "
                    "recurrence witness radius";
    return result;
  }
  const auto radius_lower = tail_majorant_detail::exact_rational_lower(radius);
  if (radius_lower.is_zero()) {
    result.detail =
        "Frobenius witness radius has no positive finite-precision lower bound";
    return result;
  }
  const auto ratio_upper =
      tail_majorant_detail::exact_rational_upper(modulus) / radius_lower;
  Magnitude ratio_gap;
  const auto geometric = tail_majorant_detail::geometric_tail_factor(
      ratio_upper, static_cast<ulong>(model.taylor_complete_max) + 1UL,
      &ratio_gap);
  if (ratio_gap.is_zero()) {
    result.detail =
        "rounded Frobenius seed modulus/witness ratio is not below one";
    return result;
  }

  options.t_order_reduction =
      solution.taylor_complete_max - model.taylor_complete_max;
  options.compute_tail_estimate = false;
  const auto restricted = restrict_exact_local_epsilon(solution, model.epsilon);
  const auto finite = evaluate_local_solution(restricted, point, options);
  result.value.frame = finite.value.epsilon;
  result.value.guarantee = ErrorGuarantee::Certified;
  result.value.absolute.assign(finite.value.epsilon.width(), Magnitude::zero());

  auto sigma = options.imaginary_sign.has_value()
                   ? options.imaginary_sign
                   : local_detail::chart_imaginary_sign(solution);
  if (point.sign < 0 && local_detail::branch_sensitive(solution) &&
      !sigma.has_value())
    throw std::domain_error(
        "negative Frobenius seed has no analytic-continuation rim");
  ComplexBall logarithm(0);
  if (point.sign != 0) {
    logarithm = local_detail::cb_log(point.modulus);
    if (point.sign < 0 && sigma.has_value())
      logarithm += local_detail::imaginary_pi(*sigma);
  }

  for (const auto &tower : model.towers) {
    const auto a = ComplexBall::from_strings(tower.a_exact);
    const auto b = ComplexBall::from_strings(tower.b_exact);
    auto t_to_a = local_detail::cb_pow(point.modulus, a);
    if (point.sign < 0)
      t_to_a *= local_detail::cb_exp(
          local_detail::imaginary_pi(sigma.value_or(1)) * a);
    const auto t_to_a_upper = Magnitude::upper_abs(t_to_a);
    const auto b_log = b * logarithm;
    for (std::uint32_t log_power = 0; log_power <= tower.log_complete_max;
         ++log_power) {
      auto outer = local_detail::cb_pow_ui(logarithm, log_power);
      for (std::uint32_t divisor = 2; divisor <= log_power; ++divisor)
        outer = local_detail::cb_div_ui(outer, divisor);
      const auto outer_upper = t_to_a_upper * Magnitude::upper_abs(outer);
      for (std::int32_t output_power = result.value.frame.min_power;
           output_power <= result.value.frame.complete_max; ++output_power) {
        ComplexBall exponential_coefficient(1);
        const auto jmax64 = static_cast<std::int64_t>(output_power) -
                            static_cast<std::int64_t>(log_power) -
                            model.epsilon.min_power;
        if (jmax64 < 0)
          continue;
        for (std::int64_t j = 0; j <= jmax64; ++j) {
          const auto input_power64 = static_cast<std::int64_t>(output_power) -
                                     static_cast<std::int64_t>(log_power) - j;
          if (input_power64 >= model.epsilon.min_power &&
              input_power64 <= model.epsilon.complete_max) {
            const auto prefix_index = static_cast<std::size_t>(
                input_power64 - model.epsilon.min_power);
            const auto analytic_tail =
                tower.induction_m_upper.at(prefix_index) * geometric;
            result.value.absolute.at(static_cast<std::size_t>(
                output_power - result.value.frame.min_power)) +=
                outer_upper * Magnitude::upper_abs(exponential_coefficient) *
                analytic_tail;
          }
          exponential_coefficient =
              local_detail::cb_div_ui(exponential_coefficient * b_log,
                                      static_cast<std::uint32_t>(j + 1));
        }
      }
    }
  }
  result.value.provenance =
      "certified singular Frobenius seed-value tail; " + model.provenance;
  result.status = TailMajorantStatus::Certified;
  result.detail =
      "all omitted Frobenius coefficients are bounded at the ordinary seed";
  return result;
}

inline CertifiedSingularSeedEvaluation
evaluate_singular_rational_shadow_with_certified_tail(
    const SingularFrobeniusTailModel &model,
    const LocalSolution<Rational> &solution, const RealEvaluationPoint &point,
    EvaluationOptions options = {}) {
  options.t_order_reduction =
      solution.taylor_complete_max - model.taylor_complete_max;
  options.compute_tail_estimate = false;
  CertifiedSingularSeedEvaluation result;
  const auto restricted =
      singular_tail_majorant_detail::restrict_exact_local_epsilon(
          solution, model.epsilon);
  result.evaluation = evaluate_local_solution(restricted, point, options);
  result.tail =
      certify_singular_frobenius_point_tail(model, solution, point, options);
  if (result.tail.status == TailMajorantStatus::Certified) {
    if (!tail_majorant_detail::same_epsilon_window(
            result.evaluation.value.epsilon, result.tail.value.frame))
      throw std::logic_error(
          "singular seed finite value and tail certificate frames disagree");
    result.evaluation.value.error = result.tail.value;
  } else {
    result.evaluation.value.error.provenance = result.tail.detail;
  }
  return result;
}

// The ordinary leg of a singular-to-ordinary bridge is still governed by an
// exact Rational q/C recurrence.  A disk-wide Gronwall norm can lose the
// cancellations in that recurrence by many orders of magnitude after a
// singular gauge transformation.  Reuse the all-future transfer theorem
// above with the ordinary tower (a=b=0, log=0), and take the induction seed
// from the rigorously reconstructed Acb coefficients.  This proves the tail
// of every solution in the propagated seed balls; it does not replace an Acb
// enclosure by a midpoint approximation.
inline RegularTaylorPointTailCertificate
certify_physical_ordinary_recurrence_point_tail(
    const PreparedPhysicalClearedODE<Rational> &equation,
    const PhysicalRegularTaylorTailModel &model,
    const RealEvaluationPoint &input_point,
    const std::string &witness_radius_exact,
    std::optional<std::uint32_t> retained_complete_max = std::nullopt) {
  using namespace singular_tail_majorant_detail;
  RegularTaylorPointTailCertificate result;
  const auto fail = [&](TailMajorantStatus status, std::string detail) {
    result.status = status;
    result.detail = std::move(detail);
    result.disk.status = status;
    result.disk.detail = result.detail;
    return result;
  };

  physical_ode_detail::validate_ode(equation);
  local_detail::validate_local_solution(model.reconstructed, false);
  const auto &solution = model.reconstructed;
  if (solution.dimension != equation.dimension ||
      solution.dimension != model.dimension ||
      equation.payload_identity != model.physical_payload_identity ||
      solution.checkpoint_identity != model.local_checkpoint_identity ||
      !tail_majorant_detail::same_epsilon_window(solution.epsilon,
                                                 model.epsilon))
    throw std::invalid_argument(
        "ordinary recurrence tail lost its equation/local provenance binding");
  if (!solution.error.empty() || solution.sectors.size() != 1 ||
      solution.sectors.front().log_power != 0 ||
      solution.sectors.front().a.is_zero != TruthValue::Yes ||
      solution.sectors.front().b.is_zero != TruthValue::Yes)
    return fail(TailMajorantStatus::Unsupported,
                "ordinary recurrence tail requires one error-free "
                "a=b=0, log=0 Taylor sector");
  if (!equation.c_lags.front().empty())
    return fail(TailMajorantStatus::Unsupported,
                "ordinary recurrence tail requires C(0,epsilon)=0");

  const auto retained =
      retained_complete_max.value_or(model.taylor_complete_max);
  if (retained > model.taylor_complete_max ||
      retained > solution.taylor_complete_max)
    throw std::invalid_argument(
        "ordinary recurrence tail retained order exceeds its reconstructed "
        "local");
  const auto degree = polynomial_lag_degree(equation);
  if (degree != 0 && retained + 1 < degree)
    return fail(TailMajorantStatus::Unsupported,
                "ordinary recurrence prefix is shorter than its polynomial "
                "lag window");

  const auto point = line_integration_detail::require_exact_rational_point(
      input_point, "ordinary-recurrence-tail-evaluation");
  const Rational exact_point(point.exact_coordinate);
  const auto modulus = tail_majorant_detail::abs_rational(exact_point);
  Rational radius(0);
  try {
    radius = Rational(witness_radius_exact);
  } catch (const std::invalid_argument &) {
    throw std::invalid_argument(
        "ordinary recurrence witness radius must be exact rational");
  }
  if (radius.sign() <= 0)
    throw std::invalid_argument(
        "ordinary recurrence witness radius must be positive");
  result.disk.witness_radius_exact = radius.str();
  if (!(modulus < radius))
    return fail(TailMajorantStatus::Inconclusive,
                "ordinary recurrence evaluation point is not strictly "
                "inside its witness disk");
  if (!solution.chart.infinite_radius &&
      !arb_lt(acb_realref(ComplexBall::from_strings(radius.str()).raw()),
              acb_realref(solution.chart.radius.raw())))
    return fail(TailMajorantStatus::Inconclusive,
                "ordinary recurrence witness disk is not strictly inside "
                "the retained chart radius");

  const auto common_valuation = common_equation_valuation(equation);
  const auto adjusted_q0 =
      shifted_valuation(equation.q_lags.front(), common_valuation);
  if (adjusted_q0.zero || adjusted_q0.valuation != 0)
    return fail(TailMajorantStatus::Unsupported,
                "ordinary recurrence common cancellation does not leave "
                "q(0,epsilon) as a formal unit");
  for (const auto &value : equation.q_lags)
    if (!value.zero &&
        shifted_valuation(value, common_valuation).valuation < 0)
      return fail(TailMajorantStatus::Unsupported,
                  "ordinary recurrence has a negative q epsilon valuation");
  for (const auto &lag : equation.c_lags)
    for (const auto &entry : lag)
      if (shifted_valuation(entry.value, common_valuation).valuation < 0)
        return fail(TailMajorantStatus::Unsupported,
                    "ordinary recurrence has a negative C epsilon valuation");

  const auto n0 = static_cast<std::uint64_t>(retained) + 1;
  const auto future_inverse = certify_structured_future_inverse(
      equation, model.epsilon, Rational(0), Rational(0), 0,
      common_valuation, n0);
  if (future_inverse.status != TailMajorantStatus::Certified)
    return fail(TailMajorantStatus::Inconclusive,
                "ordinary all-future inverse is inconclusive: " +
                    future_inverse.detail);
  const auto contraction = certify_interval_transfer_contraction(
      equation, model.epsilon, Rational(0), Rational(0), 0,
      common_valuation, n0, radius);
  if (contraction.status != TailMajorantStatus::Certified)
    return fail(TailMajorantStatus::Inconclusive,
                "ordinary all-future contraction is inconclusive: " +
                    contraction.detail);

  const auto radius_upper =
      tail_majorant_detail::exact_rational_upper(radius);
  const auto epsilon_weight =
      Magnitude::from_ui(contraction.epsilon_weight_base);
  std::vector<Magnitude> induction_m(model.epsilon.width(),
                                     Magnitude::zero());
  const auto &sector = solution.sectors.front();
  for (std::size_t prefix_index = 0;
       prefix_index < model.epsilon.width(); ++prefix_index) {
    const auto rho = contraction.contraction_upper.at(prefix_index);
    if (Magnitude::positive_difference_lower(Magnitude::one(), rho)
            .is_zero())
      return fail(
          TailMajorantStatus::Inconclusive,
          "ordinary all-future contraction rho is not below one on epsilon "
          "prefix ending at " +
              std::to_string(
                  static_cast<std::int64_t>(model.epsilon.min_power) +
                  static_cast<std::int64_t>(prefix_index)) +
              "; rho_upper=" + std::to_string(rho.approximate_upper()));
    if (degree == 0)
      continue;
    const auto first = retained + 1 - degree;
    for (std::uint32_t n = first; n <= retained; ++n) {
      auto weighted_prefix = Magnitude::zero();
      for (std::size_t epsilon_index = 0;
           epsilon_index <= prefix_index; ++epsilon_index)
        for (std::uint32_t component = 0;
             component < solution.dimension; ++component)
          weighted_prefix = Magnitude::maximum(
              weighted_prefix,
              Magnitude::upper_abs(
                  sector.coefficients[local_detail::sector_index(
                      solution, epsilon_index, n, component)]) /
                  epsilon_weight.power_upper(
                      static_cast<ulong>(epsilon_index)));
      const auto coefficient_bound =
          weighted_prefix *
          epsilon_weight.power_upper(
              static_cast<ulong>(prefix_index));
      induction_m[prefix_index] = Magnitude::maximum(
          induction_m[prefix_index],
          coefficient_bound *
              radius_upper.power_upper(static_cast<ulong>(n)));
    }
  }

  const auto radius_lower =
      tail_majorant_detail::exact_rational_lower(radius);
  if (radius_lower.is_zero())
    return fail(TailMajorantStatus::Inconclusive,
                "ordinary recurrence witness radius has no positive lower "
                "bound");
  const auto ratio_upper =
      tail_majorant_detail::exact_rational_upper(modulus) / radius_lower;
  Magnitude gap_lower;
  const auto first_unseen = static_cast<ulong>(retained) + 1UL;
  const auto value_factor = tail_majorant_detail::geometric_tail_factor(
      ratio_upper, first_unseen, &gap_lower);
  if (gap_lower.is_zero())
    return fail(TailMajorantStatus::Inconclusive,
                "ordinary recurrence evaluation ratio is not below one");
  const auto inverse_gap = gap_lower.reciprocal_upper();
  const auto theta_factor =
      ratio_upper.power_upper(first_unseen) *
      (Magnitude::from_ui(first_unseen) * inverse_gap +
       ratio_upper * inverse_gap * inverse_gap);

  result.value.frame = model.epsilon;
  result.value.guarantee = ErrorGuarantee::Certified;
  result.theta.frame = model.epsilon;
  result.theta.guarantee = ErrorGuarantee::Certified;
  result.value.absolute.reserve(induction_m.size());
  result.theta.absolute.reserve(induction_m.size());
  for (const auto &bound : induction_m) {
    result.value.absolute.push_back(bound * value_factor);
    result.theta.absolute.push_back(bound * theta_factor);
  }
  const auto provenance =
      "certified ordinary all-future q/C recurrence tail; "
      "epsilon_weight_base=" +
      std::to_string(contraction.epsilon_weight_base) +
      "; retained_taylor_complete_max=" + std::to_string(retained) +
      "; witness_radius_exact=" + radius.str() +
      "; physical_payload_identity=" + equation.payload_identity;
  result.value.provenance = provenance;
  result.theta.provenance = provenance;
  result.disk.status = TailMajorantStatus::Certified;
  result.disk.cauchy_circle_upper = induction_m;
  result.disk.detail =
      "all-future exact q/C recurrence contraction; " +
      contraction.detail;
  result.status = TailMajorantStatus::Certified;
  result.detail =
      "all omitted ordinary coefficients are bounded by the exact q/C "
      "recurrence";
  return result;
}

inline CertifiedLocalEvaluation
evaluate_physical_local_solution_with_recurrence_tail(
    const PreparedPhysicalClearedODE<Rational> &equation,
    const PhysicalRegularTaylorTailModel &model,
    const RealEvaluationPoint &point,
    const std::string &witness_radius_exact,
    EvaluationOptions options = {}) {
  if (options.t_order_reduction > model.taylor_complete_max)
    throw std::invalid_argument(
        "ordinary recurrence tail evaluation reduction exceeds retained "
        "order");
  const auto retained =
      model.taylor_complete_max - options.t_order_reduction;
  options.compute_tail_estimate = false;
  CertifiedLocalEvaluation result;
  result.evaluation =
      evaluate_local_solution(model.reconstructed, point, options);
  result.tail = certify_physical_ordinary_recurrence_point_tail(
      equation, model, point, witness_radius_exact, retained);
  if (result.tail.status != TailMajorantStatus::Certified) {
    result.evaluation.value.error.provenance = result.tail.detail;
    result.evaluation.theta_value.error.provenance = result.tail.detail;
    return result;
  }
  if (result.evaluation.value.epsilon.min_power <
          result.tail.value.frame.min_power ||
      result.evaluation.value.epsilon.complete_max >
          result.tail.value.frame.complete_max ||
      result.evaluation.theta_value.epsilon.min_power <
          result.tail.theta.frame.min_power ||
      result.evaluation.theta_value.epsilon.complete_max >
          result.tail.theta.frame.complete_max)
    throw std::logic_error(
        "ordinary recurrence tail does not cover the evaluated epsilon "
        "windows");
  ErrorEnvelope value_error;
  value_error.frame = result.evaluation.value.epsilon;
  value_error.guarantee = ErrorGuarantee::Certified;
  value_error.provenance = result.tail.value.provenance;
  for (std::int64_t power = value_error.frame.min_power;
       power <= value_error.frame.complete_max; ++power)
    value_error.absolute.push_back(result.tail.value.absolute.at(
        static_cast<std::size_t>(
            power - result.tail.value.frame.min_power)));
  result.evaluation.value.error = std::move(value_error);
  ErrorEnvelope theta_error;
  theta_error.frame = result.evaluation.theta_value.epsilon;
  theta_error.guarantee = ErrorGuarantee::Certified;
  theta_error.provenance = result.tail.theta.provenance;
  for (std::int64_t power = theta_error.frame.min_power;
       power <= theta_error.frame.complete_max; ++power)
    theta_error.absolute.push_back(result.tail.theta.absolute.at(
        static_cast<std::size_t>(
            power - result.tail.theta.frame.min_power)));
  result.evaluation.theta_value.error = std::move(theta_error);
  return result;
}

inline EpsilonVector inflate_certified_singular_seed_value(
    const CertifiedSingularSeedEvaluation &certified) {
  if (certified.tail.status != TailMajorantStatus::Certified ||
      certified.evaluation.value.error.guarantee != ErrorGuarantee::Certified ||
      !tail_majorant_detail::same_epsilon_window(
          certified.evaluation.value.epsilon,
          certified.evaluation.value.error.frame) ||
      certified.evaluation.value.error.absolute.size() !=
          certified.evaluation.value.epsilon.width())
    throw std::invalid_argument(
        "singular seed value has no complete certified tail envelope");
  auto inflated = certified.evaluation.value;
  inflated.error = ErrorEnvelope{};
  for (std::size_t epsilon_index = 0; epsilon_index < inflated.epsilon.width();
       ++epsilon_index)
    for (std::uint32_t component = 0; component < inflated.dimension;
         ++component)
      certified.tail.value.absolute[epsilon_index].add_error_to(
          inflated
              .coefficients[epsilon_index * inflated.dimension + component]);
  return inflated;
}

template <typename PhysicalizeSeed>
inline CertifiedSingularOrdinaryBridge
certify_singular_rational_shadow_ordinary_bridge_with_seed_map(
    const PreparedPhysicalClearedODE<Rational> &tail_equation,
    const LocalSolution<Rational> &tail_solution,
    const PreparedPhysicalClearedODE<Rational> &physical_equation,
    PhysicalizeSeed &&physicalize_seed, EpsilonWindow claimed_epsilon,
    std::uint32_t singular_taylor_complete_max,
    const RealEvaluationPoint &seed_point,
    const std::string &singular_witness_radius_exact,
    const RealEvaluationPoint &target_point,
    std::uint32_t ordinary_taylor_complete_max,
    const std::string &ordinary_witness_radius_exact,
    EvaluationOptions singular_options = {},
    const SingularFrobeniusPrefixCertificate *prefix_certificate = nullptr,
    bool allow_disk_gronwall_fallback = true) {
  CertifiedSingularOrdinaryBridge result;
  result.seed_exact = seed_point.exact_coordinate;
  result.singular_witness_radius_exact = singular_witness_radius_exact;
  result.ordinary_witness_radius_exact = ordinary_witness_radius_exact;
  result.singular_taylor_complete_max = singular_taylor_complete_max;
  result.ordinary_taylor_complete_max = ordinary_taylor_complete_max;
  const auto fail = [&](TailMajorantStatus status, std::string detail) {
    result.status = status;
    result.detail = std::move(detail);
    return result;
  };
  if (seed_point.exact_coordinate.empty() ||
      target_point.exact_coordinate.empty())
    throw std::invalid_argument("singular ordinary bridge requires exact "
                                "rational seed and target points");
  const Rational seed(seed_point.exact_coordinate);
  const Rational target(target_point.exact_coordinate);
  if (seed.is_zero())
    return fail(TailMajorantStatus::Unsupported,
                "singular ordinary bridge seed is the singular center");
  const auto delta = target - seed;
  result.target_delta_exact = delta.str();
  const auto normalized_tail_equation =
      cancel_physical_cleared_ode_common_epsilon_monomial(tail_equation);
  const auto normalized_physical_equation =
      cancel_physical_cleared_ode_common_epsilon_monomial(physical_equation);

  auto singular_model = prepare_singular_rational_shadow_tail_model(
      normalized_tail_equation, tail_solution, claimed_epsilon,
      singular_taylor_complete_max, singular_witness_radius_exact,
      prefix_certificate);
  if (singular_model.status != TailMajorantStatus::Certified ||
      !singular_model.model.has_value())
    return fail(singular_model.status,
                "singular seed model: " + singular_model.detail);
  auto singular_evaluation =
      evaluate_singular_rational_shadow_with_certified_tail(
          *singular_model.model, tail_solution, seed_point, singular_options);
  result.singular_tail = singular_evaluation.tail;
  if (singular_evaluation.tail.status != TailMajorantStatus::Certified)
    return fail(singular_evaluation.tail.status,
                "singular seed evaluation: " + singular_evaluation.tail.detail);
  auto tail_seed_value =
      inflate_certified_singular_seed_value(singular_evaluation);
  auto seed_value = std::forward<PhysicalizeSeed>(physicalize_seed)(
      tail_seed_value, seed_point);
  if (!seed_value.has_value())
    return fail(TailMajorantStatus::Unsupported,
                "singular seed frame conversion was unavailable");

  auto recentered = recenter_singular_physical_cleared_ode_to_ordinary(
      normalized_physical_equation, seed);
  if (!recentered.eligible || !recentered.equation.has_value())
    return fail(TailMajorantStatus::Unsupported,
                "singular-to-ordinary equation translation: " +
                    recentered.reason);
  const auto evolution = evolve_ordinary_center_value(
      *recentered.equation, *seed_value, ordinary_taylor_complete_max);
  if (!evolution.eligible)
    return fail(TailMajorantStatus::Unsupported,
                "ordinary bridge evolution: " + evolution.reason);

  const auto seed_modulus = tail_majorant_detail::abs_rational(seed);
  const auto delta_modulus = tail_majorant_detail::abs_rational(delta);
  if (!(delta_modulus < seed_modulus))
    return fail(TailMajorantStatus::Unsupported,
                "ordinary target is not inside the disk excluding the original "
                "singular center");
  ChartGeometry ordinary_chart;
  ordinary_chart.center_exact = seed.str();
  ordinary_chart.scale_exact = "1";
  ordinary_chart.radius_exact = seed_modulus.str();
  ordinary_chart.radius = ComplexBall::from_strings(seed_modulus.str());
  ordinary_chart.infinite_radius = false;
  auto ordinary_local = ordinary_evolution_local_solution(
      evolution, std::move(ordinary_chart), {},
      tail_solution.checkpoint_identity +
          ":certified-singular-ordinary-bridge:" + seed.str());
  const auto numeric_equation =
      specialize_rational_physical_cleared_ode_to_acb(*recentered.equation);
  auto ordinary_model = prepare_physical_regular_homogeneous_tail_model(
      numeric_equation, ordinary_local);
  if (ordinary_model.status != TailMajorantStatus::Certified ||
      !ordinary_model.model.has_value())
    return fail(ordinary_model.status,
                "ordinary bridge tail model: " + ordinary_model.detail);
  EvaluationOptions ordinary_options;
  ordinary_options.compute_tail_estimate = false;
  auto ordinary_evaluation =
      evaluate_physical_local_solution_with_recurrence_tail(
          *recentered.equation, *ordinary_model.model,
          RealEvaluationPoint::rational(delta.str()),
          ordinary_witness_radius_exact, ordinary_options);
  std::string ordinary_tail_method = "all-future-q/C-recurrence";
  std::string recurrence_tail_detail;
  if (ordinary_evaluation.tail.status != TailMajorantStatus::Certified) {
    recurrence_tail_detail = ordinary_evaluation.tail.detail;
    if (!allow_disk_gronwall_fallback)
      return fail(
          ordinary_evaluation.tail.status,
          "ordinary recurrence bridge target evaluation: " +
              recurrence_tail_detail);
    ordinary_evaluation =
        evaluate_physical_local_solution_with_certified_tail(
          *ordinary_model.model, RealEvaluationPoint::rational(delta.str()),
          ordinary_witness_radius_exact, ordinary_options);
    ordinary_tail_method = "disk-gronwall-fallback";
  }
  result.ordinary_tail = ordinary_evaluation.tail;
  if (ordinary_evaluation.tail.status != TailMajorantStatus::Certified)
    return fail(ordinary_evaluation.tail.status,
                "ordinary bridge target evaluation: " +
                    ordinary_evaluation.tail.detail +
                    (recurrence_tail_detail.empty()
                         ? std::string()
                         : "; recurrence fallback reason: " +
                               recurrence_tail_detail));

  auto final_value = ordinary_evaluation.evaluation.value;
  const auto certified_error = final_value.error;
  final_value.error = ErrorEnvelope{};
  if (!tail_majorant_detail::same_epsilon_window(
          final_value.epsilon, certified_error.frame) ||
      certified_error.guarantee != ErrorGuarantee::Certified ||
      certified_error.absolute.size() !=
          final_value.epsilon.width())
    throw std::logic_error(
        "ordinary bridge finite value and tail frames disagree");
  for (std::size_t epsilon_index = 0;
       epsilon_index < final_value.epsilon.width(); ++epsilon_index)
    for (std::uint32_t component = 0; component < final_value.dimension;
         ++component)
      certified_error.absolute[epsilon_index].add_error_to(
          final_value
              .coefficients[epsilon_index * final_value.dimension + component]);
  result.evaluation = std::move(ordinary_evaluation.evaluation);
  result.evaluation.value = std::move(final_value);
  result.evaluation.value.error = ErrorEnvelope{};
  result.status = TailMajorantStatus::Certified;
  result.detail = "certified Frobenius seed translated and propagated through "
                  "an ordinary physical q/C chart; ordinary_tail_method=" +
                  ordinary_tail_method +
                  (recurrence_tail_detail.empty()
                       ? std::string()
                       : "; recurrence_fallback_reason=" +
                             recurrence_tail_detail);
  return result;
}

inline CertifiedSingularOrdinaryBridge
certify_singular_rational_shadow_ordinary_bridge(
    const PreparedPhysicalClearedODE<Rational> &equation,
    const LocalSolution<Rational> &solution, EpsilonWindow claimed_epsilon,
    std::uint32_t singular_taylor_complete_max,
    const RealEvaluationPoint &seed_point,
    const std::string &singular_witness_radius_exact,
    const RealEvaluationPoint &target_point,
    std::uint32_t ordinary_taylor_complete_max,
    const std::string &ordinary_witness_radius_exact,
    EvaluationOptions singular_options = {},
    const SingularFrobeniusPrefixCertificate *prefix_certificate = nullptr,
    bool allow_disk_gronwall_fallback = true) {
  return certify_singular_rational_shadow_ordinary_bridge_with_seed_map(
      equation, solution, equation,
      [](const EpsilonVector &seed, const RealEvaluationPoint &)
          -> std::optional<EpsilonVector> { return seed; },
      claimed_epsilon, singular_taylor_complete_max, seed_point,
      singular_witness_radius_exact, target_point, ordinary_taylor_complete_max,
      ordinary_witness_radius_exact, singular_options, prefix_certificate,
      allow_disk_gronwall_fallback);
}

} // namespace diffexp::kernel
