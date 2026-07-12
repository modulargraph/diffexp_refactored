#pragma once

#include "diffexp2/integration.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace diffexp2 {

// This layer integrates exactly the Taylor coefficients retained by a typed
// LocalSolution.  It deliberately makes no statement about the unseen Taylor
// tail: a future line/tile driver may attach a separate tail majorant without
// weakening the finite-frame contract implemented here.
enum class LineIntegrationScope : std::uint8_t { StoredTruncation };

struct StoredLineIntegrationOptions {
  // The caller, rather than the primitive, owns the epsilon delivery contract.
  // Rows below min_power are requested structural zeros; every row through
  // complete_max must be computable without reading an unknown coefficient.
  EpsilonWindow delivered_epsilon;
  std::optional<std::int32_t> imaginary_sign;
};

struct StoredLineIntegrationDiagnostics {
  std::size_t input_monomial_cells = 0;
  std::size_t grouped_monomials = 0;
  std::size_t zero_groups_skipped = 0;
  std::size_t cancelled_divergent_groups = 0;
  std::size_t primitive_evaluations = 0;
  std::size_t primitive_component_applications = 0;
  std::size_t primitive_component_reuses = 0;
  bool has_center_endpoint = false;
  std::string detail;
};

struct StoredLineIntegral {
  EpsilonVector value;
  LineIntegrationScope scope = LineIntegrationScope::StoredTruncation;
  std::optional<std::int32_t> imaginary_sign;
  StoredLineIntegrationDiagnostics diagnostics;
};

namespace line_integration_detail {

inline RealEvaluationPoint require_exact_rational_point(
    const RealEvaluationPoint& point, const char* label) {
  if (point.exact_coordinate.empty())
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::InvalidInterval, "E9",
        std::string(label) +
            " endpoint has no retained exact local-coordinate identity");

  Rational exact(0);
  try {
    exact = Rational(point.exact_coordinate);
  } catch (const std::invalid_argument&) {
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::InvalidInterval, "E9",
        std::string(label) +
            " endpoint is not an exact rational local coordinate");
  }
  const auto normalized = RealEvaluationPoint::rational(exact.str());
  if (normalized.sign != point.sign ||
      !acb_equal(normalized.modulus.raw(), point.modulus.raw()))
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::InvalidInterval, "E9",
        std::string(label) +
            " endpoint ball/sign disagrees with its retained exact "
            "local-coordinate identity");
  return normalized;
}

template <typename Scalar>
void require_inside_chart(const LocalSolution<Scalar>& solution,
                          const RealEvaluationPoint& point,
                          const char* label) {
  if (solution.chart.infinite_radius) return;
  if (!arb_lt(acb_realref(point.modulus.raw()),
              acb_realref(solution.chart.radius.raw())))
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::InvalidInterval, "E9",
        std::string(label) +
            " endpoint is not provably inside the finite chart radius");
}

inline bool exact_singleton_zero(const Rational& value) {
  return value.is_zero();
}

inline bool exact_singleton_zero(const ComplexBall& value) {
  return value.is_zero();
}

template <typename Scalar>
ComplexBall as_ball(const Scalar& value) {
  return local_detail::to_ball(value);
}

using MonomialKey =
    std::tuple<std::string, std::uint8_t, std::string, std::uint32_t>;

template <typename Scalar>
struct MonomialGroup {
  SectorMonomialTag tag;
  // Flat [input epsilon][component], component fastest.
  std::vector<Scalar> coefficients;
  std::size_t cell_count = 0;
  bool had_material_input = false;
};

inline MonomialKey monomial_key(const SectorMonomialTag& tag) {
  return {tag.m.canonical, static_cast<std::uint8_t>(tag.b.domain),
          tag.b.canonical, tag.log_power};
}

inline bool same_exact_descriptor(const ExactScalarDescriptor& left,
                                  const ExactScalarDescriptor& right) {
  if (left.domain != right.domain || left.canonical != right.canonical ||
      left.symbols != right.symbols || left.is_zero != right.is_zero ||
      left.is_integer != right.is_integer || left.sign != right.sign ||
      left.specialization.has_value() != right.specialization.has_value())
    return false;
  return !left.specialization.has_value() ||
         acb_equal(left.specialization->raw(), right.specialization->raw());
}

template <typename Scalar>
std::map<MonomialKey, MonomialGroup<Scalar>> group_monomials(
    const LocalSolution<Scalar>& solution,
    StoredLineIntegrationDiagnostics* diagnostics) {
  std::map<MonomialKey, MonomialGroup<Scalar>> groups;
  const auto coefficient_count =
      solution.epsilon.width() * solution.dimension;
  for (const auto& sector : solution.sectors) {
    for (std::size_t n = 0; n < solution.taylor_width(); ++n) {
      auto tag = sector_monomial_tag(
          sector, static_cast<std::uint32_t>(n));
      integration_detail::validate_tag(tag);
      const auto key = monomial_key(tag);
      auto found = groups.find(key);
      if (found == groups.end()) {
        MonomialGroup<Scalar> group;
        group.tag = std::move(tag);
        group.coefficients.assign(coefficient_count,
                                  ScalarTraits<Scalar>::zero());
        found = groups.emplace(key, std::move(group)).first;
      } else if (!same_exact_descriptor(found->second.tag.b, tag.b)) {
        throw NativeIntegrationError(
            NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
            "equal exact monomial keys carry inconsistent regulator "
            "descriptor facts or numeric specializations");
      }
      auto& group = found->second;
      ++group.cell_count;
      ++diagnostics->input_monomial_cells;
      for (std::size_t ei = 0; ei < solution.epsilon.width(); ++ei) {
        for (std::uint32_t component = 0; component < solution.dimension;
             ++component) {
          const auto& value =
              sector.coefficients[local_detail::sector_index(
                  solution, ei, n, component)];
          if (!exact_singleton_zero(value)) group.had_material_input = true;
          group.coefficients[ei * solution.dimension + component] += value;
        }
      }
    }
  }
  diagnostics->grouped_monomials = groups.size();
  return groups;
}

template <typename Scalar>
void require_integrable_first_unseen_taylor(
    const LocalSolution<Scalar>& solution, bool has_center_endpoint) {
  if (!has_center_endpoint) return;
  // The first unseen monomial is t^(a + N + 1), and its primitive is
  // integrable exactly when a + N + 2 is positive.  This is a structural
  // exponent check, not a tail-size claim.
  const Rational unseen_shift(
      std::to_string(solution.taylor_width() + 1));
  for (const auto& sector : solution.sectors) {
    if (sector.b.is_zero != TruthValue::Yes) continue;
    if (sector.a.domain != ExactDomain::Rational)
      throw NativeIntegrationError(
          NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
          "center-touching line integration requires a rational local "
          "power for the first-unseen Taylor check");
    const auto first_unseen_alpha =
        Rational(sector.a.canonical) + unseen_shift;
    if (ExactScalarDescriptor::rational(first_unseen_alpha.str()).sign !=
        ExactSign::Positive)
      throw NativeIntegrationError(
          NativeIntegrationErrorCode::IncompleteTaylorWindow, "E10",
          "center-touching stored Taylor window is incomplete: the first "
          "unseen b=0 monomial is not provably integrable");
  }
}

inline bool center_divergent(const SectorMonomialTag& tag,
                             bool has_center_endpoint) {
  return has_center_endpoint && tag.b.is_zero == TruthValue::Yes &&
         tag.alpha0.sign != ExactSign::Positive;
}

inline std::int32_t primitive_min_power(const SectorMonomialTag& tag,
                                        bool has_center_endpoint) {
  // Only a regulated t^(-1+b eps) primitive at a center endpoint exposes the
  // genuine 1/eps row.  Same-side and crossing primitives pair their two
  // boundaries before returning and therefore start at eps^p.
  if (has_center_endpoint && tag.alpha0.is_zero == TruthValue::Yes &&
      tag.b.is_zero == TruthValue::No)
    return -1;
  return local_detail::checked_i32(tag.log_power,
                                   "line primitive epsilon minimum");
}

inline std::int32_t checked_difference(std::int32_t left,
                                       std::int32_t right,
                                       const char* label) {
  return local_detail::checked_i32(
      static_cast<std::int64_t>(left) - right, label);
}

template <typename Scalar>
std::optional<std::int32_t> first_nonzero_power(
    const LocalSolution<Scalar>& solution,
    const MonomialGroup<Scalar>& group, std::uint32_t component) {
  for (std::int64_t power64 = solution.epsilon.min_power;
       power64 <= solution.epsilon.complete_max; ++power64) {
    const auto power = static_cast<std::int32_t>(power64);
    const auto ei = static_cast<std::size_t>(
        static_cast<std::int64_t>(power) - solution.epsilon.min_power);
    if (!exact_singleton_zero(
            group.coefficients[ei * solution.dimension + component]))
      return power;
  }
  return std::nullopt;
}

template <typename Scalar>
[[noreturn]] void throw_divergent_group(
    const LocalSolution<Scalar>& solution,
    const MonomialGroup<Scalar>& group,
    const std::vector<std::optional<std::int32_t>>& first_nonzero) {
  for (std::uint32_t component = 0; component < solution.dimension;
       ++component) {
    if (!first_nonzero[component].has_value()) continue;
    const auto power = *first_nonzero[component];
    const auto ei = static_cast<std::size_t>(
        static_cast<std::int64_t>(power) - solution.epsilon.min_power);
    const auto& coefficient =
        group.coefficients[ei * solution.dimension + component];
    bool uncertified = false;
    if constexpr (std::is_same_v<Scalar, ComplexBall>)
      uncertified = coefficient.contains_zero();
    NativeIntegrationError error(
        uncertified
            ? NativeIntegrationErrorCode::UncertifiedCancellation
            : NativeIntegrationErrorCode::DivergentEndpoint,
        uncertified ? "E10" : "E2",
        uncertified
            ? "grouped Acb coefficient of a divergent center monomial "
              "contains zero but is not the exact singleton zero"
            : "grouped coefficient of a divergent center monomial is "
              "nonzero");
    error.absolute_power = group.tag.m.canonical;
    error.log_power = group.tag.log_power;
    error.epsilon_power = power;
    error.component = component;
    throw error;
  }
  throw std::logic_error("divergent group has no nonzero coefficient");
}

}  // namespace line_integration_detail

template <typename Scalar>
StoredLineIntegral integrate_stored_local_line(
    const LocalSolution<Scalar>& solution,
    const RealEvaluationPoint& lower_input,
    const RealEvaluationPoint& upper_input,
    const StoredLineIntegrationOptions& options) {
  static_assert(std::is_same_v<Scalar, Rational> ||
                    std::is_same_v<Scalar, ComplexBall>,
                "stored line integration supports Rational or ComplexBall "
                "coefficient frames only");
  using namespace line_integration_detail;

  validate_local_solution(solution, false);
  if (!solution.error.empty())
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
        "stored line integration cannot discard an input error envelope; "
        "native tail/error propagation is not implemented yet");
  (void)options.delivered_epsilon.width();
  if (options.imaginary_sign.has_value() &&
      *options.imaginary_sign != 1 && *options.imaginary_sign != -1)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
        "line imaginary sign must be exactly +1 or -1");

  const auto lower = require_exact_rational_point(lower_input, "lower");
  const auto upper = require_exact_rational_point(upper_input, "upper");
  integration_detail::validate_interval(lower, upper);
  require_inside_chart(solution, lower, "lower");
  require_inside_chart(solution, upper, "upper");

  std::optional<std::int32_t> chart_sign;
  try {
    chart_sign = derive_chart_imaginary_sign(solution);
  } catch (const std::domain_error& error) {
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
        std::string("invalid prepared chart branch prescription: ") +
            error.what());
  }
  if (chart_sign.has_value() && options.imaginary_sign.has_value() &&
      *chart_sign != *options.imaginary_sign)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
        "explicit line branch sign conflicts with the prepared chart");

  StoredLineIntegral result;
  result.value.epsilon = options.delivered_epsilon;
  result.value.dimension = solution.dimension;
  result.value.coefficients.assign(
      options.delivered_epsilon.width() * solution.dimension,
      ComplexBall(0));
  result.diagnostics.has_center_endpoint =
      lower.sign == 0 || upper.sign == 0;
  require_integrable_first_unseen_taylor(
      solution, result.diagnostics.has_center_endpoint);

  auto groups = group_monomials(solution, &result.diagnostics);
  const bool has_negative_arm = lower.sign < 0 || upper.sign < 0;
  auto effective_sign = options.imaginary_sign.has_value()
                            ? options.imaginary_sign
                            : chart_sign;
  if (has_negative_arm && !effective_sign.has_value())
    effective_sign = 1;  // Mathematica-compatible principal +i0 rim.
  result.imaginary_sign = effective_sign;

  for (const auto& [key, group] : groups) {
    (void)key;
    std::vector<std::optional<std::int32_t>> first_nonzero(
        solution.dimension);
    std::size_t material_components = 0;
    std::int32_t earliest_coefficient =
        std::numeric_limits<std::int32_t>::max();
    for (std::uint32_t component = 0; component < solution.dimension;
         ++component) {
      first_nonzero[component] =
          first_nonzero_power(solution, group, component);
      if (first_nonzero[component].has_value()) {
        ++material_components;
        earliest_coefficient =
            std::min(earliest_coefficient, *first_nonzero[component]);
      }
    }
    const bool divergent = center_divergent(
        group.tag, result.diagnostics.has_center_endpoint);
    if (material_components != 0 && divergent)
      throw_divergent_group(solution, group, first_nonzero);

    const auto primitive_min = primitive_min_power(
        group.tag, result.diagnostics.has_center_endpoint);
    const auto deliverable_max = local_detail::checked_i32(
        static_cast<std::int64_t>(solution.epsilon.complete_max) +
            primitive_min,
        "line deliverable epsilon maximum");
    if (deliverable_max < options.delivered_epsilon.complete_max) {
      std::ostringstream detail;
      detail << "delivered epsilon complete_max "
             << options.delivered_epsilon.complete_max
             << " requires coefficient complete_max at least "
             << (static_cast<std::int64_t>(
                     options.delivered_epsilon.complete_max) -
                 primitive_min)
             << " because the grouped primitive begins at eps^"
             << primitive_min << "; retained complete_max is "
             << solution.epsilon.complete_max;
      throw NativeIntegrationError(
          NativeIntegrationErrorCode::IncompleteEpsilonWindow, "E10",
          detail.str());
    }

    // Even an all-zero retained group cannot bypass the halo gate: its next
    // unknown epsilon coefficient can feed a requested row through a
    // negative-power primitive.
    if (material_components == 0) {
      ++result.diagnostics.zero_groups_skipped;
      if (divergent && group.had_material_input && group.cell_count > 1)
        ++result.diagnostics.cancelled_divergent_groups;
      continue;
    }

    MonomialIntegrationOptions primitive_options;
    primitive_options.imaginary_sign = effective_sign;
    primitive_options.complete_max = std::max(
        primitive_min,
        checked_difference(options.delivered_epsilon.complete_max,
                           earliest_coefficient,
                           "required primitive complete maximum"));
    const auto primitive = integrate_sector_monomial(
        group.tag, lower, upper, primitive_options);
    if (primitive.min_power() != primitive_min)
      throw std::logic_error(
          "line primitive minimum disagrees with its exact preflight");
    ++result.diagnostics.primitive_evaluations;
    result.diagnostics.primitive_component_applications +=
        material_components;
    result.diagnostics.primitive_component_reuses +=
        material_components - 1;

    for (std::uint32_t component = 0; component < solution.dimension;
         ++component) {
      if (!first_nonzero[component].has_value()) continue;
      const auto first = *first_nonzero[component];
      std::vector<ComplexBall> coefficient_values;
      coefficient_values.reserve(static_cast<std::size_t>(
          static_cast<std::int64_t>(solution.epsilon.complete_max) - first +
          1));
      for (std::int64_t power64 = first;
           power64 <= solution.epsilon.complete_max; ++power64) {
        const auto power = static_cast<std::int32_t>(power64);
        const auto ei = static_cast<std::size_t>(
            static_cast<std::int64_t>(power) -
            solution.epsilon.min_power);
        coefficient_values.push_back(as_ball(
            group.coefficients[ei * solution.dimension + component]));
      }
      const EpsilonFrame<ComplexBall> coefficient_frame(
          first, std::move(coefficient_values));
      const auto contribution = coefficient_frame * primitive;
      if (contribution.complete_max() <
          options.delivered_epsilon.complete_max)
        throw std::logic_error(
            "line epsilon preflight did not deliver its promised window");
      for (std::int64_t power64 = options.delivered_epsilon.min_power;
           power64 <= options.delivered_epsilon.complete_max; ++power64) {
        const auto power = static_cast<std::int32_t>(power64);
        result.value.at(power, component) += contribution.coefficient(power);
      }
    }
  }

  result.value.error.guarantee = ErrorGuarantee::None;
  result.value.error.provenance =
      "stored Taylor truncation only; no unseen-tail majorant or full-local "
      "certificate";
  result.diagnostics.detail = result.value.error.provenance;
  return result;
}

}  // namespace diffexp2
