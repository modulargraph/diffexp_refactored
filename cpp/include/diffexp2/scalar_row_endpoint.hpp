#pragma once

#include "diffexp2/local_algebra.hpp"
#include "diffexp2/recurrence.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace diffexp2 {

// The centered endpoint only consumes projected Taylor cells whose absolute
// power is nonpositive.  This versioned mode records whether the row was
// contracted directly through that exact liveness ceiling or whether the
// existing full materialization was used as a strict semantic fallback.
enum class PreparedScalarRowEndpointMode : std::uint8_t {
  CenteredB0RationalTaylorLivenessV1,
  MaterializedFallbackV1
};

inline const char* prepared_scalar_row_endpoint_mode_name(
    PreparedScalarRowEndpointMode mode) {
  switch (mode) {
    case PreparedScalarRowEndpointMode::
        CenteredB0RationalTaylorLivenessV1:
      return "centered-b0-rational-taylor-liveness-v1";
    case PreparedScalarRowEndpointMode::MaterializedFallbackV1:
      return "materialized-fallback-v1";
  }
  throw std::logic_error("unknown prepared scalar-row endpoint mode");
}

template <typename Scalar>
struct PreparedScalarRowEndpointResult {
  EndpointLimitResult endpoint;
  // The endpoint provenance needs the exact projected sector tags.  In the
  // liveness mode this slab contains only the nonpositive-power Taylor cells,
  // while its tag set and epsilon frame equal the full materialized oracle.
  LocalSolution<Scalar> projected;
  std::int32_t projected_top_valid = kCompleteInfinity;
  PreparedScalarRowEndpointMode mode =
      PreparedScalarRowEndpointMode::MaterializedFallbackV1;
  std::uint32_t source_taylor_complete_max = 0;
  std::uint32_t projected_taylor_complete_max = 0;
  std::string fallback_reason;
};

struct PreparedScalarRowEndpointPlan {
  std::int32_t projected_top_valid = kCompleteInfinity;
  PreparedScalarRowEndpointMode mode =
      PreparedScalarRowEndpointMode::MaterializedFallbackV1;
  std::optional<std::uint32_t> taylor_liveness_cap;
  std::uint32_t source_taylor_complete_max = 0;
  std::uint32_t projected_taylor_complete_max = 0;
  std::string fallback_reason;
};

namespace scalar_row_endpoint_detail {

inline std::int32_t shifted_validity(std::int32_t validity,
                                     std::int32_t shift) {
  if (validity == kCompleteInfinity) return kCompleteInfinity;
  return local_algebra_detail::checked_i32(
      static_cast<std::int64_t>(validity) + shift,
      "rational-row output validity");
}

template <typename Scalar>
LocalSolution<Scalar> zero_scalar_like(
    const LocalSolution<Scalar>& source,
    const std::string& checkpoint_identity) {
  auto result = local_algebra_detail::with_selected_component(source, 0);
  for (auto& sector : result.sectors)
    std::fill(sector.coefficients.begin(), sector.coefficients.end(),
              ScalarTraits<Scalar>::zero());
  result.checkpoint_identity = checkpoint_identity;
  validate_local_solution(result, false);
  return result;
}

template <typename Scalar>
std::int32_t row_top_valid(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& matrix,
    std::int32_t source_top_valid) {
  std::int32_t output_top_valid = matrix.entries.empty()
      ? source_top_valid : kCompleteInfinity;
  // This deliberately includes entries proved structurally zero.  The
  // retained materialized-row pipeline binds validity before its algebraic
  // zero pruning, and endpoint behavior must remain byte-for-byte compatible.
  for (const auto& entry : matrix.entries)
    output_top_valid = std::min(
        output_top_valid,
        shifted_validity(source_top_valid,
                         entry.multiplier.epsilon_shift));
  return output_top_valid;
}

template <typename Scalar>
LocalSolution<Scalar> restrict_to_top_valid(
    LocalSolution<Scalar> solution, std::int32_t output_top_valid,
    const std::string& checkpoint_identity) {
  output_top_valid = std::min(
      output_top_valid, solution.epsilon.complete_max);
  if (output_top_valid < solution.epsilon.min_power)
    throw std::domain_error(
        "rational-row application has no valid output epsilon coefficient");
  if (output_top_valid < solution.epsilon.complete_max)
    solution = restrict_local_epsilon_frame_strict_lower(
        solution, solution.epsilon.min_power, output_top_valid,
        checkpoint_identity);
  if (solution.dimension != 1)
    throw std::logic_error(
        "rational-row application did not produce a scalar local solution");
  return solution;
}

template <typename Scalar>
LocalSolution<Scalar> materialized_projection(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& matrix,
    const LocalSolution<Scalar>& input,
    std::int32_t output_top_valid,
    const std::string& checkpoint_identity) {
  auto applied = apply_prepared_sparse_local_matrix(
      matrix, input, checkpoint_identity);
  auto solution = applied.has_value()
      ? std::move(*applied)
      : zero_scalar_like(input, checkpoint_identity);
  return restrict_to_top_valid(
      std::move(solution), output_top_valid, checkpoint_identity);
}

// Contract one prepared row through a Taylor ceiling while retaining the
// exact arithmetic order of the generic materialized pipeline:
//
//   entry -> source sector -> epsilon convolution -> Taylor convolution,
//
// followed by one whole-slab addition when identical projected tags merge.
// In particular, Acb does not use an addmul reassociation here.  A temporary
// slab for one entry/sector preserves the old product-then-combine rounding
// order without allocating the full Taylor tower for every term.
template <typename Scalar>
std::optional<LocalSolution<Scalar>> capped_projection_materialized_order(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& matrix,
    const LocalSolution<Scalar>& input,
    std::uint32_t target_taylor_complete_max,
    const std::string& checkpoint_identity) {
  validate_local_solution(input, false);
  if (!input.error.empty())
    throw std::invalid_argument(
        "native local matrix application needs explicit error-envelope propagation");
  if (matrix.rows != 1 || matrix.columns != input.dimension)
    throw std::invalid_argument(
        "prepared scalar row dimensions disagree with its local");
  if (target_taylor_complete_max > input.taylor_complete_max)
    throw std::invalid_argument(
        "prepared scalar-row Taylor liveness ceiling exceeds its source window");

  const auto epsilon_width = input.epsilon.width();
  const auto input_taylor_width = input.taylor_width();
  std::int32_t natural_min = std::numeric_limits<std::int32_t>::max();
  std::int32_t natural_complete = std::numeric_limits<std::int32_t>::max();
  bool active = false;
  for (const auto& entry : matrix.entries) {
    if (entry.row != 0 || entry.column >= matrix.columns)
      throw std::invalid_argument(
          "prepared local matrix entry is out of range");
    if (entry.multiplier.structurally_zero()) continue;
    active = true;
    if (entry.multiplier.kernels.size() < epsilon_width)
      throw std::invalid_argument(
          "prepared rational multiplier has too few kernels");
    for (std::size_t epsilon = 0; epsilon < epsilon_width; ++epsilon)
      if (entry.multiplier.kernels[epsilon].size() <
          input_taylor_width)
        throw std::invalid_argument(
            "prepared rational multiplier has too few Taylor coefficients");
    natural_min = std::min(
        natural_min, local_algebra_detail::checked_i32(
                         static_cast<std::int64_t>(
                             input.epsilon.min_power) +
                             entry.multiplier.epsilon_shift,
                         "rational-product epsilon minimum"));
    natural_complete = std::min(
        natural_complete, local_algebra_detail::checked_i32(
                              static_cast<std::int64_t>(
                                  input.epsilon.complete_max) +
                                  entry.multiplier.epsilon_shift,
                              "rational-product epsilon maximum"));
  }
  if (!active) return std::nullopt;

  LocalSolution<Scalar> output;
  output.chart = input.chart;
  output.epsilon = {natural_min, natural_complete};
  output.taylor_complete_max = target_taylor_complete_max;
  output.dimension = 1;
  output.prescriptions = input.prescriptions;
  output.checkpoint_identity = checkpoint_identity;
  output.sectors.reserve(input.sectors.size());
  const auto output_taylor_width = output.taylor_width();

  for (const auto& entry : matrix.entries) {
    if (entry.multiplier.structurally_zero()) continue;
    const auto term_min = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(input.epsilon.min_power) +
            entry.multiplier.epsilon_shift,
        "rational-product epsilon minimum");
    for (const auto& source_sector : input.sectors) {
      LocalSector<Scalar> term;
      term.a = local_algebra_detail::subtract_nonnegative_integer(
          source_sector.a, entry.multiplier.center_pole_order);
      term.b = source_sector.b;
      term.log_power = source_sector.log_power;
      term.coefficients.assign(output.sector_size(),
                               ScalarTraits<Scalar>::zero());

      for (std::size_t out_epsilon = 0;
           out_epsilon < epsilon_width; ++out_epsilon) {
        const auto output_power = local_algebra_detail::checked_i32(
            static_cast<std::int64_t>(term_min) + out_epsilon,
            "rational-product epsilon power");
        if (output_power > natural_complete) break;
        const auto aligned_epsilon = static_cast<std::size_t>(
            output_power - natural_min);
        for (std::size_t kernel_epsilon = 0;
             kernel_epsilon <= out_epsilon; ++kernel_epsilon) {
          const auto input_epsilon = out_epsilon - kernel_epsilon;
          const auto& kernel =
              entry.multiplier.kernels[kernel_epsilon];
          for (std::size_t n = 0; n < output_taylor_width; ++n) {
            for (std::size_t m = 0; m <= n; ++m) {
              if (ScalarTraits<Scalar>::is_zero(kernel[m])) continue;
              term.coefficients[local_algebra_detail::flat_index(
                  aligned_epsilon, n, 0, output_taylor_width, 1)] +=
                  kernel[m] * source_sector.coefficients[
                      local_algebra_detail::flat_index(
                          input_epsilon, n - m, entry.column,
                          input_taylor_width, input.dimension)];
            }
          }
        }
      }

      const auto found = std::find_if(
          output.sectors.begin(), output.sectors.end(),
          [&](const auto& candidate) {
            return local_algebra_detail::same_sector_tag(candidate, term);
          });
      if (found == output.sectors.end()) {
        output.sectors.push_back(std::move(term));
      } else {
        for (std::size_t coefficient = 0;
             coefficient < found->coefficients.size(); ++coefficient)
          found->coefficients[coefficient] += term.coefficients[coefficient];
      }
    }
  }

  std::stable_sort(output.sectors.begin(), output.sectors.end(),
                   local_algebra_detail::sector_less<Scalar>);
  validate_local_solution(output, false);
  return output;
}

template <typename Scalar>
struct LivenessDecision {
  std::optional<std::uint32_t> cap;
  std::string fallback_reason;
};

template <typename Scalar>
LivenessDecision<Scalar> centered_liveness_decision(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& matrix,
    const LocalSolution<Scalar>& input) {
  if (matrix.rows != 1 || matrix.columns != input.dimension)
    return {std::nullopt, "non-scalar-or-mismatched-row"};
  if (matrix.entries.empty())
    return {std::nullopt, "empty-row"};

  bool active = false;
  std::uint32_t cap = 0;
  bool has_live_cell = false;
  const auto first_unseen = input.taylor_width();
  for (const auto& entry : matrix.entries) {
    if (entry.row != 0 || entry.column >= matrix.columns)
      return {std::nullopt, "invalid-row-entry"};
    if (entry.multiplier.structurally_zero()) continue;
    active = true;
    for (const auto& sector : input.sectors) {
      if (sector.b.domain != ExactDomain::Rational ||
          sector.b.is_zero != TruthValue::Yes ||
          sector.b.canonical != "0")
        return {std::nullopt, "nonzero-or-unsupported-regulator"};
      if (sector.a.domain != ExactDomain::Rational)
        return {std::nullopt, "non-rational-local-power"};
      const auto shifted = Rational(sector.a.canonical) - Rational(
          std::to_string(entry.multiplier.center_pole_order));
      std::optional<std::size_t> positive;
      for (std::size_t n = 0; n <= first_unseen; ++n) {
        if ((shifted + Rational(std::to_string(n))).sign() > 0) {
          positive = n;
          break;
        }
      }
      // A leading positive power has no live endpoint cell, and a tower which
      // never reaches a positive first-unseen cell needs the materialized
      // oracle to preserve its exact empty/incomplete behavior.
      if (!positive.has_value())
        return {std::nullopt, "incomplete-taylor-tower"};
      if (*positive == 0)
        return {std::nullopt, "positive-leading-power"};
      if (*positive == first_unseen)
        return {std::nullopt, "no-taylor-reduction"};
      has_live_cell = true;
      cap = std::max(
          cap, static_cast<std::uint32_t>(*positive - 1));
    }
  }
  if (!active) return {std::nullopt, "structurally-zero-row"};
  if (!has_live_cell || cap >= input.taylor_complete_max)
    return {std::nullopt, "no-taylor-reduction"};
  return {cap, {}};
}

template <typename Scalar>
std::int32_t projected_complete_max(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& matrix,
    const LocalSolution<Scalar>& input) {
  std::int32_t complete_max = input.epsilon.complete_max;
  bool active = false;
  for (const auto& entry : matrix.entries) {
    if (entry.multiplier.structurally_zero()) continue;
    const auto term_complete = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(input.epsilon.complete_max) +
            entry.multiplier.epsilon_shift,
        "rational-product epsilon maximum");
    complete_max = active ? std::min(complete_max, term_complete)
                          : term_complete;
    active = true;
  }
  // With no active term the materialized compatibility path retains the
  // selected first source component as an exact-zero representative.
  return complete_max;
}

}  // namespace scalar_row_endpoint_detail

template <typename Scalar>
PreparedScalarRowEndpointPlan centered_prepared_scalar_row_endpoint_plan(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& matrix,
    const LocalSolution<Scalar>& input,
    std::int32_t source_top_valid) {
  validate_local_solution(input, false);
  const auto output_top_valid =
      scalar_row_endpoint_detail::row_top_valid(matrix, source_top_valid);
  const auto decision =
      scalar_row_endpoint_detail::centered_liveness_decision(matrix, input);

  PreparedScalarRowEndpointPlan plan;
  plan.projected_top_valid = std::min(
      output_top_valid,
      scalar_row_endpoint_detail::projected_complete_max(matrix, input));
  plan.source_taylor_complete_max = input.taylor_complete_max;
  if (decision.cap.has_value()) {
    plan.mode = PreparedScalarRowEndpointMode::
        CenteredB0RationalTaylorLivenessV1;
    plan.taylor_liveness_cap = decision.cap;
    plan.projected_taylor_complete_max = *decision.cap;
  } else {
    plan.mode = PreparedScalarRowEndpointMode::MaterializedFallbackV1;
    plan.projected_taylor_complete_max = input.taylor_complete_max;
    plan.fallback_reason = decision.fallback_reason;
  }
  return plan;
}

// Reconstruct only the exact projected chart/tag/prescription shape.  This is
// used to validate compact endpoint checkpoints without repeating either the
// finite convolution or the endpoint limit during restore.
template <typename Scalar>
LocalSolution<Scalar> prepared_scalar_row_endpoint_analytic_shape(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& matrix,
    const LocalSolution<Scalar>& input,
    std::string checkpoint_identity = {}) {
  validate_local_solution(input, false);
  if (matrix.rows != 1 || matrix.columns != input.dimension)
    throw std::invalid_argument(
        "prepared scalar row dimensions disagree with its local");

  LocalSolution<Scalar> shape;
  shape.chart = input.chart;
  shape.epsilon = {0, 0};
  shape.taylor_complete_max = 0;
  shape.dimension = 1;
  shape.prescriptions = input.prescriptions;
  shape.checkpoint_identity = checkpoint_identity.empty()
      ? input.checkpoint_identity +
            ":centered-endpoint-analytic-shape:" + matrix.exact_identity
      : std::move(checkpoint_identity);

  bool active = false;
  for (const auto& entry : matrix.entries) {
    if (entry.row != 0 || entry.column >= matrix.columns)
      throw std::invalid_argument(
          "prepared local matrix entry is out of range");
    if (entry.multiplier.structurally_zero()) continue;
    active = true;
    for (const auto& source_sector : input.sectors) {
      LocalSector<Scalar> sector;
      sector.a = local_algebra_detail::subtract_nonnegative_integer(
          source_sector.a, entry.multiplier.center_pole_order);
      sector.b = source_sector.b;
      sector.log_power = source_sector.log_power;
      sector.coefficients.assign(1, ScalarTraits<Scalar>::zero());
      shape.sectors.push_back(std::move(sector));
    }
  }
  if (!active) {
    for (const auto& source_sector : input.sectors) {
      LocalSector<Scalar> sector;
      sector.a = source_sector.a;
      sector.b = source_sector.b;
      sector.log_power = source_sector.log_power;
      sector.coefficients.assign(1, ScalarTraits<Scalar>::zero());
      shape.sectors.push_back(std::move(sector));
    }
  }
  return canonicalize_identical_local_sectors(std::move(shape));
}

template <typename Scalar>
PreparedScalarRowEndpointResult<Scalar>
centered_prepared_scalar_row_endpoint_limit(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& matrix,
    const LocalSolution<Scalar>& input,
    std::int32_t source_top_valid,
    const EndpointLimitOptions& options,
    std::string checkpoint_identity = {}) {
  static_assert(std::is_same_v<Scalar, Rational> ||
                    std::is_same_v<Scalar, ComplexBall>,
                "prepared scalar-row endpoint limits require Rational or Acb coefficients");
  validate_local_solution(input, false);
  if (!input.error.empty())
    throw std::domain_error(
        "native rational-row application requires explicit source error-envelope propagation");
  if (checkpoint_identity.empty())
    checkpoint_identity = input.checkpoint_identity +
        ":centered-endpoint-row:" + matrix.exact_identity;

  const auto plan = centered_prepared_scalar_row_endpoint_plan(
      matrix, input, source_top_valid);

  PreparedScalarRowEndpointResult<Scalar> result;
  result.source_taylor_complete_max = plan.source_taylor_complete_max;
  if (plan.taylor_liveness_cap.has_value()) {
    auto projected =
        scalar_row_endpoint_detail::capped_projection_materialized_order(
            matrix, input, *plan.taylor_liveness_cap,
            checkpoint_identity);
    if (!projected.has_value())
      throw std::logic_error(
          "eligible centered scalar row lost every active entry");
    result.projected =
        scalar_row_endpoint_detail::restrict_to_top_valid(
            std::move(*projected), plan.projected_top_valid,
            checkpoint_identity);
  } else {
    result.projected =
        scalar_row_endpoint_detail::materialized_projection(
            matrix, input, plan.projected_top_valid,
            checkpoint_identity);
  }
  result.mode = plan.mode;
  result.fallback_reason = plan.fallback_reason;
  result.projected_top_valid = plan.projected_top_valid;
  result.projected_taylor_complete_max =
      result.projected.taylor_complete_max;
  result.endpoint = endpoint_sector_limit(result.projected, options);
  return result;
}

}  // namespace diffexp2
