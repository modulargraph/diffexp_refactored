#pragma once

#include "diffexp2/integration.hpp"
#include "diffexp2/local_algebra.hpp"
#include "diffexp2/physical_ode.hpp"

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
enum class LineIntegrationScope : std::uint8_t {
  StoredTruncation,
  FullLocalWithCertifiedTail
};

struct StoredLineIntegrationOptions {
  // The caller, rather than the primitive, owns the epsilon delivery contract.
  // Rows below min_power are requested structural zeros; every row through
  // complete_max must be computable without reading an unknown coefficient.
  EpsilonWindow delivered_epsilon;
  // Most callers require the complete delivery window above exactly.  A
  // transport observable instead has a requested ceiling and a smaller
  // public completeness floor.  When this floor is present, fused row
  // integration may lower the returned complete_max to the exact maximum
  // supported by all retained monomial primitives, but never below this
  // value.  This charges the eps^-1 halo only to regulated centre primitives
  // without discarding an upper row from ordinary tiles.
  std::optional<std::int32_t> required_complete_max;
  std::optional<std::int32_t> imaginary_sign;
  // A retained tile plan may bind a genuinely algebraic chart scale.  Its
  // exact string remains opaque; this independently certified sign is the
  // only scale fact needed to validate the prepared branch rim.
  std::optional<std::int32_t> certified_chart_scale_sign;
  // Generic native line integration remains exact-singleton strict.  The
  // only caller allowed to populate this policy is the explicitly declared
  // FT transport-pair observable path.  Its tolerance and producer
  // provenance are therefore part of the retained request, never a hidden
  // global threshold.
  struct BoundedDivergentCancellation {
    Magnitude relative_tolerance;
    std::string relative_tolerance_text;
    std::string provenance;
  };
  std::optional<BoundedDivergentCancellation> divergent_cancellation;
  // Terminal factorized transport evaluates the physical basis columns
  // independently and contracts them only after exact monomial grouping.
  // In that one internal mode a divergent group must be returned to the
  // caller instead of being rejected column by column.
  bool defer_divergent_groups = false;
};

struct StoredLineIntegrationDiagnostics {
  std::size_t input_monomial_cells = 0;
  std::size_t grouped_monomials = 0;
  std::size_t zero_groups_skipped = 0;
  std::size_t cancelled_divergent_groups = 0;
  std::size_t bounded_cancelled_divergent_coefficients = 0;
  std::size_t primitive_evaluations = 0;
  std::size_t primitive_component_applications = 0;
  std::size_t primitive_component_reuses = 0;
  bool has_center_endpoint = false;
  bool tail_certificate_requested = false;
  std::string tail_certificate_status = "not-requested";
  std::string tail_witness_radius_exact;
  std::string divergent_cancellation_mode = "exact-singleton";
  std::string divergent_relative_tolerance;
  std::string divergent_cancellation_provenance;
  std::string detail;
};

struct StoredLineIntegral {
  struct DeferredDivergentGroup {
    SectorMonomialTag tag;
    EpsilonFrame<ComplexBall> coefficients;
    std::vector<Magnitude> contribution_scale_uppers;
    std::size_t cell_count = 0;
    bool had_material_input = false;
  };

  EpsilonVector value;
  LineIntegrationScope scope = LineIntegrationScope::StoredTruncation;
  std::optional<std::int32_t> imaginary_sign;
  StoredLineIntegrationDiagnostics diagnostics;
  std::vector<DeferredDivergentGroup> deferred_divergent_groups;
};

inline StoredLineIntegral certified_zero_physical_line(
    const EpsilonWindow& epsilon, std::uint32_t dimension,
    std::optional<std::int32_t> imaginary_sign,
    bool has_center_endpoint, bool tail_certificate_requested) {
  StoredLineIntegral zero;
  zero.value.epsilon = epsilon;
  zero.value.dimension = dimension;
  zero.value.coefficients.assign(
      epsilon.width() * dimension, ComplexBall(0));
  zero.value.error.frame = epsilon;
  zero.value.error.guarantee = ErrorGuarantee::Certified;
  zero.value.error.absolute.assign(epsilon.width(), Magnitude::zero());
  zero.value.error.provenance =
      "exact zero certified physical tile; no primitive or unseen-tail contribution";
  zero.scope = LineIntegrationScope::FullLocalWithCertifiedTail;
  zero.imaginary_sign = imaginary_sign;
  zero.diagnostics.has_center_endpoint = has_center_endpoint;
  zero.diagnostics.tail_certificate_requested =
      tail_certificate_requested;
  zero.diagnostics.tail_certificate_status =
      "exact-zero-certified-interval";
  zero.diagnostics.detail = zero.value.error.provenance;
  return zero;
}

namespace line_integration_detail {

inline RealEvaluationPoint require_exact_rational_point(
    const RealEvaluationPoint& point, const char* label) {
  if (point.exact_coordinate.empty())
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::InvalidInterval, "E9",
        std::string(label) +
            " endpoint has no retained exact local-coordinate identity");

  if (point.certified_algebraic) {
    if (!arb_is_zero(acb_imagref(point.modulus.raw())) ||
        !arb_is_nonnegative(acb_realref(point.modulus.raw())) ||
        (point.sign == 0 && !point.modulus.is_zero()) ||
        (point.sign != 0 && point.modulus.contains_zero()))
      throw NativeIntegrationError(
          NativeIntegrationErrorCode::InvalidInterval, "E9",
          std::string(label) +
              " certified algebraic endpoint has an invalid real modulus/sign specialization");
    return point;
  }

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
  fmpq_t exact_modulus;
  fmpq_init(exact_modulus);
  const auto parse_status =
      fmpq_set_str(exact_modulus, exact.str().c_str(), 10);
  if (parse_status == 0) {
    fmpq_canonicalise(exact_modulus);
    fmpq_abs(exact_modulus, exact_modulus);
  }
  const bool specialization_contains_identity =
      parse_status == 0 &&
      acb_contains_fmpq(point.modulus.raw(), exact_modulus);
  fmpq_clear(exact_modulus);
  // The exact coordinate is the identity authority.  Its Acb specialization
  // may have been created under a different precision lease, so bitwise
  // equality with a freshly rounded ball is neither necessary nor stable.
  // Require the retained enclosure to contain the exact modulus instead.
  if (normalized.sign != point.sign ||
      !specialization_contains_identity)
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
  // Per [input epsilon][component], this is the rigorous l1 upper magnitude
  // of the ungrouped contributions.  A cancellation residual accumulates
  // the uncertainty of every addend, so a maximum would make the admissible
  // relative error spuriously tighter by the number of terms.  The scale is
  // deliberately kept separate from the grouped coefficient so cancellation
  // cannot shrink its own reference norm.
  std::vector<Magnitude> contribution_scale_uppers;
  std::size_t cell_count = 0;
  bool had_material_input = false;
};

inline MonomialKey monomial_key(const SectorMonomialTag& tag) {
  return {tag.m.canonical, static_cast<std::uint8_t>(tag.b.domain),
          tag.b.canonical, tag.log_power};
}

// Compact bounded-memory index for the fused transport path.  Exact tag
// strings remain owned once by OutputSector; each projected sector/Taylor
// cell contributes only this 16-byte POD record.
struct FusedMonomialWorkItem {
  std::uint64_t stable_key_hash = 0;
  std::uint32_t sector_ordinal = 0;
  std::uint32_t taylor_order = 0;
};

static_assert(std::is_trivially_copyable_v<FusedMonomialWorkItem>);
static_assert(sizeof(FusedMonomialWorkItem) == 16);

inline void fnv1a64_append_byte(std::uint64_t& hash, std::uint8_t byte) {
  hash ^= byte;
  hash *= UINT64_C(1099511628211);
}

inline void fnv1a64_append_u64(std::uint64_t& hash, std::uint64_t value) {
  for (unsigned shift = 0; shift < 64; shift += 8)
    fnv1a64_append_byte(
        hash, static_cast<std::uint8_t>((value >> shift) & UINT64_C(0xff)));
}

inline void fnv1a64_append_string(std::uint64_t& hash,
                                  const std::string& value) {
  fnv1a64_append_u64(hash, value.size());
  for (const auto byte : value)
    fnv1a64_append_byte(hash, static_cast<std::uint8_t>(byte));
}

inline std::uint64_t stable_monomial_key_hash(
    const SectorMonomialTag& tag) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  fnv1a64_append_string(hash, tag.m.canonical);
  fnv1a64_append_byte(hash, static_cast<std::uint8_t>(tag.b.domain));
  fnv1a64_append_string(hash, tag.b.canonical);
  fnv1a64_append_u64(hash, tag.log_power);
  return hash;
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
        group.contribution_scale_uppers.assign(coefficient_count,
                                               Magnitude::zero());
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
          const auto cell = ei * solution.dimension + component;
          group.coefficients[cell] += value;
          group.contribution_scale_uppers[cell] =
              group.contribution_scale_uppers[cell] +
              Magnitude::upper_abs(as_ball(value));
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
bool material_group(const MonomialGroup<Scalar>& group) {
  return std::any_of(group.coefficients.begin(), group.coefficients.end(),
                     [](const Scalar& value) {
                       return !exact_singleton_zero(value);
                     });
}

template <typename Scalar>
bool accept_or_throw_divergent_group(
    const LocalSolution<Scalar>& solution,
    const MonomialGroup<Scalar>& group,
    const StoredLineIntegrationOptions& options,
    StoredLineIntegrationDiagnostics* diagnostics) {
  bool had_material_coefficient = false;
  for (std::uint32_t component = 0; component < solution.dimension;
       ++component) {
    for (std::int64_t power64 = solution.epsilon.min_power;
         power64 <= solution.epsilon.complete_max; ++power64) {
      const auto power = static_cast<std::int32_t>(power64);
      const auto ei = static_cast<std::size_t>(
          power64 - solution.epsilon.min_power);
      const auto cell = ei * solution.dimension + component;
      const auto& coefficient = group.coefficients[cell];
      if (exact_singleton_zero(coefficient)) continue;
      had_material_coefficient = true;

      bool bounded = false;
      Magnitude scale = Magnitude::one();
      Magnitude bound = Magnitude::zero();
      Magnitude coefficient_upper = Magnitude::zero();
      if constexpr (std::is_same_v<Scalar, ComplexBall>) {
        if (options.divergent_cancellation.has_value()) {
          scale = Magnitude::maximum(
              Magnitude::one(), group.contribution_scale_uppers[cell]);
          bound = scale * options.divergent_cancellation->relative_tolerance;
          coefficient_upper = Magnitude::upper_abs(coefficient);
          bounded = coefficient_upper <= bound;
        }
      }
      if (bounded) {
        ++diagnostics->bounded_cancelled_divergent_coefficients;
        continue;
      }

      bool uncertified = false;
      if constexpr (std::is_same_v<Scalar, ComplexBall>)
        uncertified = coefficient.contains_zero();
      std::ostringstream detail;
      detail << (uncertified
          ? "grouped Acb coefficient of a divergent center monomial "
            "contains zero but is not certified by the active cancellation policy"
          : "grouped coefficient of a divergent center monomial is nonzero");
      if constexpr (std::is_same_v<Scalar, ComplexBall>) {
        if (options.divergent_cancellation.has_value())
          detail << "; coefficient_upper="
                 << coefficient_upper.dump_exact()
                 << "; contribution_scale_upper=" << scale.dump_exact()
                 << "; relative_tolerance="
                 << options.divergent_cancellation->relative_tolerance_text
                 << "; tolerance_bound_upper=" << bound.dump_exact();
      }
      NativeIntegrationError error(
          uncertified
              ? NativeIntegrationErrorCode::UncertifiedCancellation
              : NativeIntegrationErrorCode::DivergentEndpoint,
          uncertified ? "E10" : "E2", detail.str());
      error.absolute_power = group.tag.m.canonical;
      error.log_power = group.tag.log_power;
      error.epsilon_power = power;
      error.component = component;
      throw error;
    }
  }
  return had_material_coefficient;
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
    chart_sign = options.certified_chart_scale_sign.has_value()
        ? derive_chart_imaginary_sign(
              solution, *options.certified_chart_scale_sign)
        : derive_chart_imaginary_sign(solution);
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
  if (options.divergent_cancellation.has_value()) {
    if (!options.divergent_cancellation->relative_tolerance.is_finite() ||
        options.divergent_cancellation->relative_tolerance.is_zero() ||
        Magnitude::one() <=
            options.divergent_cancellation->relative_tolerance ||
        options.divergent_cancellation->relative_tolerance_text.empty() ||
        options.divergent_cancellation->provenance.empty())
      throw NativeIntegrationError(
          NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
          "bounded divergent-cancellation policy is malformed");
    result.diagnostics.divergent_cancellation_mode =
        "bounded-relative-acb";
    result.diagnostics.divergent_relative_tolerance =
        options.divergent_cancellation->relative_tolerance_text;
    result.diagnostics.divergent_cancellation_provenance =
        options.divergent_cancellation->provenance;
  }
  if (options.defer_divergent_groups &&
      !std::is_same_v<Scalar, ComplexBall>)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
        "deferred divergent groups are restricted to Acb terminal transport");
  require_integrable_first_unseen_taylor(
      solution, result.diagnostics.has_center_endpoint);

  auto groups = group_monomials(solution, &result.diagnostics);
  const bool has_negative_arm = lower.sign < 0 || upper.sign < 0;
  auto effective_sign = options.imaginary_sign.has_value()
                            ? options.imaginary_sign
                            : chart_sign;
  const bool material_branch_sensitive_negative_arm =
      has_negative_arm &&
      std::any_of(groups.begin(), groups.end(), [](const auto& entry) {
        const auto& group = entry.second;
        return material_group(group) &&
               integration_detail::branch_sensitive_on_negative_arm(
                   group.tag);
      });
  if (material_branch_sensitive_negative_arm && !effective_sign.has_value())
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
        "negative branch-sensitive line integration has no derivable "
        "imaginary sign");
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
    const bool bounded_divergent_cancellation =
        material_components != 0 && divergent;
    const bool deferred_divergent =
        bounded_divergent_cancellation &&
        options.defer_divergent_groups;
    if (bounded_divergent_cancellation && !deferred_divergent)
      (void)accept_or_throw_divergent_group(
          solution, group, options, &result.diagnostics);

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
    if (deferred_divergent) {
      if (solution.dimension != 1)
        throw std::logic_error(
            "deferred terminal divergence requires one scalar physical column");
      std::vector<ComplexBall> coefficients;
      coefficients.reserve(solution.epsilon.width());
      std::vector<Magnitude> scales;
      scales.reserve(solution.epsilon.width());
      for (std::size_t epsilon = 0;
           epsilon < solution.epsilon.width(); ++epsilon) {
        coefficients.push_back(
            as_ball(group.coefficients[epsilon]));
        scales.push_back(
            group.contribution_scale_uppers[epsilon]);
      }
      result.deferred_divergent_groups.push_back(
          StoredLineIntegral::DeferredDivergentGroup{
              group.tag,
              EpsilonFrame<ComplexBall>(
                  solution.epsilon, std::move(coefficients)),
              std::move(scales), group.cell_count,
              group.had_material_input});
      ++result.diagnostics.zero_groups_skipped;
      continue;
    }
    if (material_components == 0 || bounded_divergent_cancellation) {
      ++result.diagnostics.zero_groups_skipped;
      if (bounded_divergent_cancellation ||
          (divergent && group.had_material_input && group.cell_count > 1))
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

// Transport-only stored-truncation seam.  It reproduces
//
//   apply_prepared_sparse_local_matrix(row, source)
//   -> canonicalize identical output sectors
//   -> group_monomials
//   -> integrate_stored_local_line
//
// without ever owning the complete projected scalar LocalSolution.  Only
// exact-tag/contributor metadata is retained globally.  One final monomial
// group and one canonical output-sector/Taylor cell epsilon vector are live
// at a time, while divergent cells still meet in the same global monomial
// group before the cancellation gate runs.
template <typename Scalar>
StoredLineIntegral integrate_prepared_scalar_row_stored(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& matrix,
    const LocalSolution<Scalar>& source,
    std::int32_t projected_complete_cap,
    const RealEvaluationPoint& lower_input,
    const RealEvaluationPoint& upper_input,
    const StoredLineIntegrationOptions& options) {
  static_assert(std::is_same_v<Scalar, Rational> ||
                    std::is_same_v<Scalar, ComplexBall>,
                "fused row integration supports Rational or Acb only");
  using namespace line_integration_detail;
  validate_local_solution(source, false);
  if (!source.error.empty())
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
        "fused scalar-row integration cannot discard a source error envelope");
  if (matrix.rows != 1 || matrix.columns != source.dimension)
    throw std::invalid_argument(
        "fused scalar-row dimensions disagree with the source local");

  const auto epsilon_width = source.epsilon.width();
  const auto taylor_width = source.taylor_width();
  auto projected_min = source.epsilon.min_power;
  auto projected_complete = source.epsilon.complete_max;
  bool active = false;
  for (const auto& entry : matrix.entries) {
    if (entry.row != 0 || entry.column >= matrix.columns)
      throw std::invalid_argument(
          "fused scalar-row entry is out of range");
    if (entry.multiplier.structurally_zero()) continue;
    const auto multiplier_width = entry.multiplier.kernels.size();
    if (multiplier_width == 0)
      throw std::invalid_argument(
          "fused scalar-row multiplier has no epsilon kernels");
    for (std::size_t epsilon = 0; epsilon < multiplier_width; ++epsilon)
      if (entry.multiplier.kernels[epsilon].size() < taylor_width)
        throw std::invalid_argument(
            "fused scalar-row multiplier has too few Taylor coefficients");
    const auto term_min = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(source.epsilon.min_power) +
            entry.multiplier.epsilon_shift,
        "fused scalar-row epsilon minimum");
    const auto term_complete = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(term_min) +
            static_cast<std::int64_t>(
                std::min(epsilon_width, multiplier_width)) -
            1,
        "fused scalar-row epsilon complete maximum");
    if (!active) {
      projected_min = term_min;
      projected_complete = term_complete;
      active = true;
    } else {
      projected_min = std::min(projected_min, term_min);
      projected_complete = std::min(projected_complete, term_complete);
    }
  }
  projected_complete = std::min(
      projected_complete, projected_complete_cap);
  if (projected_complete < projected_min)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::IncompleteEpsilonWindow, "E10",
        "fused scalar row has no coefficient in its requested upper window");
  const EpsilonWindow projected_epsilon{
      projected_min, projected_complete};

  const auto lower = require_exact_rational_point(lower_input, "lower");
  const auto upper = require_exact_rational_point(upper_input, "upper");
  integration_detail::validate_interval(lower, upper);
  require_inside_chart(source, lower, "lower");
  require_inside_chart(source, upper, "upper");
  (void)options.delivered_epsilon.width();
  if (options.required_complete_max.has_value() &&
      (*options.required_complete_max <
           options.delivered_epsilon.min_power ||
       *options.required_complete_max >
           options.delivered_epsilon.complete_max))
    throw std::invalid_argument(
        "fused line required epsilon maximum must lie in its requested delivery window");
  if (options.imaginary_sign.has_value() &&
      *options.imaginary_sign != 1 && *options.imaginary_sign != -1)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
        "fused line imaginary sign must be exactly +1 or -1");

  std::optional<std::int32_t> chart_sign;
  try {
    chart_sign = options.certified_chart_scale_sign.has_value()
        ? derive_chart_imaginary_sign(
              source, *options.certified_chart_scale_sign)
        : derive_chart_imaginary_sign(source);
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
        "explicit fused-line branch sign conflicts with the prepared chart");
  const auto effective_sign = options.imaginary_sign.has_value()
      ? options.imaginary_sign : chart_sign;
  const bool has_negative_arm = lower.sign < 0 || upper.sign < 0;
  const bool has_center_endpoint = lower.sign == 0 || upper.sign == 0;

  StoredLineIntegral result;
  result.imaginary_sign = effective_sign;
  result.diagnostics.has_center_endpoint = has_center_endpoint;
  if (options.divergent_cancellation.has_value()) {
    if (!options.divergent_cancellation->relative_tolerance.is_finite() ||
        options.divergent_cancellation->relative_tolerance.is_zero() ||
        Magnitude::one() <=
            options.divergent_cancellation->relative_tolerance ||
        options.divergent_cancellation->relative_tolerance_text.empty() ||
        options.divergent_cancellation->provenance.empty())
      throw NativeIntegrationError(
          NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
          "bounded divergent-cancellation policy is malformed");
    result.diagnostics.divergent_cancellation_mode =
        "bounded-relative-acb";
    result.diagnostics.divergent_relative_tolerance =
        options.divergent_cancellation->relative_tolerance_text;
    result.diagnostics.divergent_cancellation_provenance =
        options.divergent_cancellation->provenance;
  }
  if (options.defer_divergent_groups &&
      !std::is_same_v<Scalar, ComplexBall>)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
        "deferred divergent groups are restricted to Acb terminal transport");

  using OutputSectorKey =
      std::tuple<std::uint8_t, std::string, std::uint8_t, std::string,
                 std::uint32_t>;
  using MatrixEntry =
      typename PreparedSparseLocalMultiplierMatrix<Scalar>::Entry;
  struct Contributor {
    const MatrixEntry* entry = nullptr;
    const LocalSector<Scalar>* sector = nullptr;
  };
  struct OutputSector {
    LocalSector<Scalar> tag;
    std::vector<Contributor> contributors;
  };
  std::map<OutputSectorKey, OutputSector> output_sectors;
  const auto admit_output_sector = [&](const LocalSector<Scalar>& sector,
                                       const MatrixEntry* entry,
                                       bool add_contributor = true) {
    LocalSector<Scalar> tag;
    tag.a = entry == nullptr
        ? sector.a
        : local_algebra_detail::subtract_nonnegative_integer(
              sector.a, entry->multiplier.center_pole_order);
    tag.b = sector.b;
    tag.log_power = sector.log_power;
    const OutputSectorKey key{
        static_cast<std::uint8_t>(tag.a.domain), tag.a.canonical,
        static_cast<std::uint8_t>(tag.b.domain), tag.b.canonical,
        tag.log_power};
    auto [found, inserted] = output_sectors.try_emplace(key);
    if (inserted) {
      found->second.tag = std::move(tag);
    } else if (!same_exact_descriptor(found->second.tag.a, tag.a) ||
               !same_exact_descriptor(found->second.tag.b, tag.b)) {
      throw NativeIntegrationError(
          NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
          "equal fused output-sector keys carry inconsistent exact descriptors");
    }
    if (entry != nullptr && add_contributor)
      found->second.contributors.push_back(Contributor{entry, &sector});
  };
  if (active) {
    for (const auto& entry : matrix.entries) {
      if (entry.multiplier.structurally_zero()) continue;
      const auto term_min = local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(source.epsilon.min_power) +
              entry.multiplier.epsilon_shift,
          "fused scalar-row admission epsilon minimum");
      for (const auto& sector : source.sectors) {
        // Match apply_prepared_scalar_row_window exactly: an active matrix
        // entry does not by itself materialize every selected source tag.
        // Admit this pair only when at least one raw finite convolution
        // product in the capped projected frame is exact-nonzero.  For Acb,
        // is_zero means the singleton zero; an enclosure containing zero is
        // deliberately retained as material.
        bool pair_can_contribute = false;
        for (std::size_t kernel_epsilon = 0;
             kernel_epsilon < entry.multiplier.kernels.size() &&
                 !pair_can_contribute;
             ++kernel_epsilon) {
          const auto output_base = static_cast<std::int64_t>(term_min) +
              static_cast<std::int64_t>(kernel_epsilon);
          if (output_base > projected_epsilon.complete_max) break;
          const auto& kernel = entry.multiplier.kernels[kernel_epsilon];
          for (std::size_t input_epsilon = 0;
               input_epsilon < epsilon_width && !pair_can_contribute;
               ++input_epsilon) {
            const auto output_power = output_base +
                static_cast<std::int64_t>(input_epsilon);
            if (output_power > projected_epsilon.complete_max) break;
            if (output_power < projected_epsilon.min_power) continue;
            for (std::size_t kernel_taylor = 0;
                 kernel_taylor < taylor_width && !pair_can_contribute;
                 ++kernel_taylor) {
              const auto& multiplier = kernel[kernel_taylor];
              if (exact_singleton_zero(multiplier)) continue;
              for (std::size_t input_taylor = 0;
                   input_taylor + kernel_taylor < taylor_width;
                   ++input_taylor) {
                const auto& coefficient = sector.coefficients[
                    local_detail::sector_index(
                        source, input_epsilon, input_taylor,
                        entry.column)];
                if (!exact_singleton_zero(coefficient)) {
                  pair_can_contribute = true;
                  break;
                }
              }
            }
          }
        }
        if (pair_can_contribute)
          admit_output_sector(sector, &entry);
      }
    }
    if (output_sectors.empty()) {
      // Preserve the same one-tag zero representative as the sparse scalar
      // projector when every active pair is exact zero in the capped frame.
      const auto active_entry = std::find_if(
          matrix.entries.begin(), matrix.entries.end(),
          [](const auto& entry) {
            return !entry.multiplier.structurally_zero();
          });
      if (active_entry == matrix.entries.end() || source.sectors.empty())
        throw std::logic_error(
            "active fused scalar row lost its zero-sector representative");
      admit_output_sector(source.sectors.front(), &*active_entry, false);
    }
  } else {
    // A structurally zero row follows exact_zero_scalar_local_like: its
    // frame and sector tags are inherited from one scalar source component.
    for (const auto& sector : source.sectors)
      admit_output_sector(sector, nullptr);
  }
  if (output_sectors.empty())
    throw std::logic_error("fused scalar row produced no output-sector metadata");

  if (has_center_endpoint) {
    const Rational unseen_shift(std::to_string(taylor_width + 1));
    for (const auto& [key, output_sector] : output_sectors) {
      (void)key;
      const auto& sector = output_sector.tag;
      if (sector.b.is_zero != TruthValue::Yes) continue;
      if (sector.a.domain != ExactDomain::Rational)
        throw NativeIntegrationError(
            NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
            "center-touching fused row requires rational local powers");
      const auto first_unseen =
          Rational(sector.a.canonical) + unseen_shift;
      if (ExactScalarDescriptor::rational(first_unseen.str()).sign !=
          ExactSign::Positive)
        throw NativeIntegrationError(
            NativeIntegrationErrorCode::IncompleteTaylorWindow, "E10",
            "center-touching fused Taylor window has a nonintegrable first unseen monomial");
    }
  }

  std::vector<const OutputSector*> canonical_output_sectors;
  canonical_output_sectors.reserve(output_sectors.size());
  for (const auto& [sector_key, output_sector] : output_sectors) {
    (void)sector_key;
    canonical_output_sectors.push_back(&output_sector);
  }
  if (canonical_output_sectors.size() >
      std::numeric_limits<std::uint32_t>::max())
    throw std::overflow_error(
        "fused output-sector count exceeds its compact work index");
  if (taylor_width != 0 && canonical_output_sectors.size() >
          std::numeric_limits<std::size_t>::max() / taylor_width)
    throw std::overflow_error("fused monomial work-item count overflows");

  std::vector<FusedMonomialWorkItem> work_items;
  work_items.reserve(canonical_output_sectors.size() * taylor_width);
  for (std::size_t ordinal = 0;
       ordinal < canonical_output_sectors.size(); ++ordinal) {
    const auto& output_sector = *canonical_output_sectors[ordinal];
    for (std::size_t n = 0; n < taylor_width; ++n) {
      auto tag = sector_monomial_tag(
          output_sector.tag, static_cast<std::uint32_t>(n));
      integration_detail::validate_tag(tag);
      work_items.push_back(FusedMonomialWorkItem{
          stable_monomial_key_hash(tag), static_cast<std::uint32_t>(ordinal),
          static_cast<std::uint32_t>(n)});
    }
  }
  std::sort(work_items.begin(), work_items.end(),
            [](const auto& left, const auto& right) {
              return std::tie(left.stable_key_hash, left.sector_ordinal,
                              left.taylor_order) <
                     std::tie(right.stable_key_hash, right.sector_ordinal,
                              right.taylor_order);
            });

  const auto tag_for_item = [&](const FusedMonomialWorkItem& item) {
    return sector_monomial_tag(
        canonical_output_sectors.at(item.sector_ordinal)->tag,
        item.taylor_order);
  };
  const auto key_for_item = [&](const FusedMonomialWorkItem& item) {
    return monomial_key(tag_for_item(item));
  };

  auto delivered_epsilon = options.delivered_epsilon;
  if (options.required_complete_max.has_value()) {
    auto exact_complete_max = delivered_epsilon.complete_max;
    for (const auto& item : work_items) {
      const auto primitive_min = primitive_min_power(
          tag_for_item(item), has_center_endpoint);
      exact_complete_max = std::min(
          exact_complete_max,
          local_detail::checked_i32(
              static_cast<std::int64_t>(projected_epsilon.complete_max) +
                  primitive_min,
              "fused exact deliverable epsilon maximum"));
    }
    if (exact_complete_max >= *options.required_complete_max)
      delivered_epsilon.complete_max = exact_complete_max;
  }
  result.value.epsilon = delivered_epsilon;
  result.value.dimension = 1;
  result.value.coefficients.assign(
      delivered_epsilon.width(), ComplexBall(0));

  // A 64-bit hash owns no exact strings and makes the normal sort entirely
  // integer-only.  Exact equality remains authoritative: detect a true hash
  // collision and sort just that bucket by the full key, with the canonical
  // sector/Taylor ordinal as the deterministic Acb addition tie-break.
  std::size_t grouped_monomials = 0;
  std::vector<std::uint64_t> collision_hashes;
  for (std::size_t bucket_begin = 0; bucket_begin < work_items.size();) {
    std::size_t bucket_end = bucket_begin + 1;
    while (bucket_end < work_items.size() &&
           work_items[bucket_end].stable_key_hash ==
               work_items[bucket_begin].stable_key_hash)
      ++bucket_end;
    const auto first_key = key_for_item(work_items[bucket_begin]);
    bool hash_collision = false;
    for (std::size_t index = bucket_begin + 1;
         index < bucket_end; ++index)
      if (key_for_item(work_items[index]) != first_key) {
        hash_collision = true;
        break;
      }
    if (hash_collision) {
      collision_hashes.push_back(
          work_items[bucket_begin].stable_key_hash);
      std::sort(
          work_items.begin() + static_cast<std::ptrdiff_t>(bucket_begin),
          work_items.begin() + static_cast<std::ptrdiff_t>(bucket_end),
          [&](const auto& left, const auto& right) {
            const auto left_key = key_for_item(left);
            const auto right_key = key_for_item(right);
            if (left_key != right_key) return left_key < right_key;
            return std::tie(left.sector_ordinal, left.taylor_order) <
                   std::tie(right.sector_ordinal, right.taylor_order);
          });
      auto previous_key = key_for_item(work_items[bucket_begin]);
      ++grouped_monomials;
      for (std::size_t index = bucket_begin + 1;
           index < bucket_end; ++index) {
        auto current_key = key_for_item(work_items[index]);
        if (current_key != previous_key) {
          ++grouped_monomials;
          previous_key = std::move(current_key);
        }
      }
    } else {
      ++grouped_monomials;
    }
    bucket_begin = bucket_end;
  }
  result.diagnostics.grouped_monomials = grouped_monomials;

  LocalSolution<Scalar> projected_shape;
  projected_shape.epsilon = projected_epsilon;
  projected_shape.dimension = 1;

  // The scalar implementation below used to recompute the complete
  // two-dimensional epsilon/Taylor convolution independently for every
  // output Taylor cell.  At the production banana settings (roughly twenty
  // epsilon rows and six hundred Taylor rows) that means tens of millions of
  // 4000-bit acb_mul calls per tile.  For Acb, form the same finite Taylor
  // convolutions with FLINT's polynomial kernel once per epsilon pair and
  // retain only the capped projected rectangle.  Epsilon ordering, Taylor
  // truncation, structural-zero decisions, and the later exact monomial
  // grouping are unchanged.
  std::vector<std::vector<ComplexBall>> acb_projected_sector_cells;
  if constexpr (std::is_same_v<Scalar, ComplexBall>) {
    const auto projected_width = projected_epsilon.width();
    acb_projected_sector_cells.resize(canonical_output_sectors.size());
    for (std::size_t ordinal = 0;
         ordinal < canonical_output_sectors.size(); ++ordinal) {
      const auto& output_sector = *canonical_output_sectors[ordinal];
      auto& cells = acb_projected_sector_cells[ordinal];
      cells.assign(projected_width * taylor_width, ComplexBall(0));
      for (const auto& contributor : output_sector.contributors) {
        const auto& entry = *contributor.entry;
        const auto& sector = *contributor.sector;
        const auto term_min = local_algebra_detail::checked_i32(
            static_cast<std::int64_t>(source.epsilon.min_power) +
                entry.multiplier.epsilon_shift,
            "fused Acb polynomial row-cell epsilon minimum");

        std::vector<local_algebra_detail::AcbPolynomial>
            kernel_polynomials;
        std::vector<local_algebra_detail::AcbPolynomial>
            source_polynomials;
        const auto multiplier_width = entry.multiplier.kernels.size();
        std::vector<bool> kernel_material(multiplier_width, false);
        std::vector<bool> source_material(epsilon_width, false);
        kernel_polynomials.reserve(multiplier_width);
        source_polynomials.reserve(epsilon_width);
        for (std::size_t epsilon = 0; epsilon < multiplier_width;
             ++epsilon) {
          kernel_polynomials.emplace_back();
          const auto& kernel = entry.multiplier.kernels[epsilon];
          for (std::size_t taylor = 0; taylor < taylor_width; ++taylor) {
            const auto& multiplier = kernel[taylor];
            if (!exact_singleton_zero(multiplier)) {
              acb_poly_set_coeff_acb(kernel_polynomials.back().raw(),
                                     static_cast<slong>(taylor),
                                     multiplier.raw());
              kernel_material[epsilon] = true;
            }
          }
        }
        for (std::size_t epsilon = 0; epsilon < epsilon_width; ++epsilon) {
          source_polynomials.emplace_back();
          for (std::size_t taylor = 0; taylor < taylor_width; ++taylor) {
            const auto& coefficient = sector.coefficients[
                local_detail::sector_index(
                    source, epsilon, taylor, entry.column)];
            if (!exact_singleton_zero(coefficient)) {
              acb_poly_set_coeff_acb(source_polynomials.back().raw(),
                                     static_cast<slong>(taylor),
                                     coefficient.raw());
              source_material[epsilon] = true;
            }
          }
        }

        local_algebra_detail::AcbPolynomial product;
        for (std::size_t kernel_epsilon = 0;
             kernel_epsilon < multiplier_width; ++kernel_epsilon) {
          if (!kernel_material[kernel_epsilon]) continue;
          for (std::size_t input_epsilon = 0;
               input_epsilon < epsilon_width; ++input_epsilon) {
            if (!source_material[input_epsilon]) continue;
            const auto raw_power = static_cast<std::int64_t>(term_min) +
                static_cast<std::int64_t>(kernel_epsilon) +
                static_cast<std::int64_t>(input_epsilon);
            if (raw_power < projected_epsilon.min_power ||
                raw_power > projected_epsilon.complete_max)
              continue;
            acb_poly_mullow(product.raw(),
                            kernel_polynomials[kernel_epsilon].raw(),
                            source_polynomials[input_epsilon].raw(),
                            static_cast<slong>(taylor_width),
                            ComplexBall::precision());
            const auto output_epsilon = static_cast<std::size_t>(
                raw_power - projected_epsilon.min_power);
            const auto product_length = std::min<std::size_t>(
                taylor_width,
                static_cast<std::size_t>(acb_poly_length(product.raw())));
            for (std::size_t taylor = 0; taylor < product_length; ++taylor) {
              ComplexBall coefficient;
              acb_poly_get_coeff_acb(coefficient.raw(), product.raw(),
                                     static_cast<slong>(taylor));
              cells[output_epsilon * taylor_width + taylor] += coefficient;
            }
          }
        }
      }
    }
  }

  // Recover the rigorous l1 norm before the fast Acb polynomial projection
  // combines its matrix/Taylor products.  Once those products have been
  // summed into one ball, |sum| no longer records the cancellation scale and
  // can make an otherwise relative policy term-count dependent.  This is
  // evaluated only for the few center-divergent monomial cells below, so the
  // production polynomial fast path remains O(N log N) for ordinary cells.
  const auto fused_cell_contribution_scale = [&source, &projected_epsilon,
      epsilon_width](const OutputSector& output_sector,
                     std::size_t output_taylor,
                     std::size_t output_epsilon) {
    Magnitude scale = Magnitude::zero();
    const auto raw_power = static_cast<std::int64_t>(
        projected_epsilon.min_power) +
        static_cast<std::int64_t>(output_epsilon);
    for (const auto& contributor : output_sector.contributors) {
      const auto& entry = *contributor.entry;
      const auto& sector = *contributor.sector;
      const auto term_min = local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(source.epsilon.min_power) +
              entry.multiplier.epsilon_shift,
          "fused cancellation-scale epsilon minimum");
      const auto term_index = raw_power - term_min;
      if (term_index < 0) continue;
      for (std::int64_t kernel_epsilon = 0;
           kernel_epsilon <= term_index; ++kernel_epsilon) {
        const auto input_epsilon = term_index - kernel_epsilon;
        if (kernel_epsilon >= static_cast<std::int64_t>(
                                  entry.multiplier.kernels.size()) ||
            input_epsilon >= static_cast<std::int64_t>(epsilon_width))
          continue;
        const auto& kernel = entry.multiplier.kernels[
            static_cast<std::size_t>(kernel_epsilon)];
        for (std::size_t kernel_taylor = 0;
             kernel_taylor <= output_taylor; ++kernel_taylor) {
          const auto& multiplier = kernel[kernel_taylor];
          if (exact_singleton_zero(multiplier)) continue;
          const auto input_taylor = output_taylor - kernel_taylor;
          const auto& coefficient = sector.coefficients[
              local_detail::sector_index(
                  source, static_cast<std::size_t>(input_epsilon),
                  input_taylor, entry.column)];
          if (exact_singleton_zero(coefficient)) continue;
          scale = scale + Magnitude::upper_abs(as_ball(multiplier)) *
              Magnitude::upper_abs(as_ball(coefficient));
        }
      }
    }
    return scale;
  };

  const auto build_group = [&](std::size_t group_begin,
                               std::size_t group_end,
                               const SectorMonomialTag& group_tag) {
    MonomialGroup<Scalar> group;
    group.tag = group_tag;
    group.coefficients.assign(projected_epsilon.width(),
                              ScalarTraits<Scalar>::zero());
    group.contribution_scale_uppers.assign(projected_epsilon.width(),
                                           Magnitude::zero());
    for (std::size_t item_index = group_begin;
         item_index < group_end; ++item_index) {
      const auto& item = work_items[item_index];
      const auto* output_sector =
          canonical_output_sectors.at(item.sector_ordinal);
      const auto taylor = item.taylor_order;
      if (!same_exact_descriptor(group.tag.b, output_sector->tag.b))
        throw NativeIntegrationError(
            NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
            "equal fused monomial keys carry inconsistent regulator descriptors");
      std::vector<Scalar> cell(projected_epsilon.width(),
                               ScalarTraits<Scalar>::zero());
      if constexpr (std::is_same_v<Scalar, ComplexBall>) {
        const auto& projected =
            acb_projected_sector_cells.at(item.sector_ordinal);
        for (std::size_t epsilon = 0;
             epsilon < projected_epsilon.width(); ++epsilon)
          cell[epsilon] = projected[epsilon * taylor_width + taylor];
      } else {
        for (const auto& contributor : output_sector->contributors) {
          const auto& entry = *contributor.entry;
          const auto& sector = *contributor.sector;
          const auto term_min = local_algebra_detail::checked_i32(
              static_cast<std::int64_t>(source.epsilon.min_power) +
                  entry.multiplier.epsilon_shift,
              "fused row-cell epsilon minimum");
          for (std::int64_t raw_power = projected_epsilon.min_power;
               raw_power <= projected_epsilon.complete_max; ++raw_power) {
            const auto term_index = raw_power - term_min;
            if (term_index < 0) continue;
            const auto output_epsilon = static_cast<std::size_t>(
                raw_power - projected_epsilon.min_power);
            for (std::int64_t kernel_epsilon = 0;
                 kernel_epsilon <= term_index; ++kernel_epsilon) {
              const auto input_epsilon = term_index - kernel_epsilon;
              if (kernel_epsilon >= static_cast<std::int64_t>(
                                        entry.multiplier.kernels.size()) ||
                  input_epsilon >=
                      static_cast<std::int64_t>(epsilon_width))
                continue;
              const auto& kernel = entry.multiplier.kernels[
                  static_cast<std::size_t>(kernel_epsilon)];
              for (std::size_t kernel_taylor = 0;
                   kernel_taylor <= taylor; ++kernel_taylor) {
                const auto& multiplier = kernel[kernel_taylor];
                if (exact_singleton_zero(multiplier)) continue;
                const auto input_taylor = taylor - kernel_taylor;
                const auto& coefficient = sector.coefficients[
                    local_detail::sector_index(
                        source, static_cast<std::size_t>(input_epsilon),
                        input_taylor, entry.column)];
                if (exact_singleton_zero(coefficient)) continue;
                local_algebra_detail::add_product(
                    cell[output_epsilon], multiplier, coefficient);
              }
            }
          }
        }
      }
      ++group.cell_count;
      for (std::size_t epsilon = 0; epsilon < cell.size(); ++epsilon) {
        const auto& value = cell[epsilon];
        if (!exact_singleton_zero(value)) group.had_material_input = true;
        group.coefficients[epsilon] += value;
        Magnitude contribution_scale = Magnitude::upper_abs(as_ball(value));
        if constexpr (std::is_same_v<Scalar, ComplexBall>) {
          if (center_divergent(group.tag, has_center_endpoint))
            contribution_scale = fused_cell_contribution_scale(
                *output_sector, taylor, epsilon);
        }
        group.contribution_scale_uppers[epsilon] =
            group.contribution_scale_uppers[epsilon] +
            contribution_scale;
      }
    }
    return group;
  };

  const auto for_each_exact_group = [&](auto&& visitor) {
    for (std::size_t group_begin = 0; group_begin < work_items.size();) {
      const auto group_tag = tag_for_item(work_items[group_begin]);
      std::size_t hash_end = group_begin + 1;
      while (hash_end < work_items.size() &&
             work_items[hash_end].stable_key_hash ==
                 work_items[group_begin].stable_key_hash)
        ++hash_end;
      if (!std::binary_search(
              collision_hashes.begin(), collision_hashes.end(),
              work_items[group_begin].stable_key_hash)) {
        visitor(group_begin, hash_end, group_tag);
        group_begin = hash_end;
        continue;
      }
      const auto group_key = monomial_key(group_tag);
      std::size_t group_end = group_begin + 1;
      while (group_end < hash_end &&
             key_for_item(work_items[group_end]) == group_key)
        ++group_end;
      visitor(group_begin, group_end, group_tag);
      group_begin = group_end;
    }
  };

  // The materialized path performs this as one global preflight before any
  // divergence or halo error.  Preserve that ordering without retaining all
  // projected coefficients: only the exceptional missing-rim path needs a
  // read-only streaming pass over the monomial groups.
  if (has_negative_arm && !effective_sign.has_value()) {
    for_each_exact_group(
        [&](std::size_t group_begin, std::size_t group_end,
            const SectorMonomialTag& group_tag) {
          if (!integration_detail::branch_sensitive_on_negative_arm(
                  group_tag))
            return;
          if (material_group(
                  build_group(group_begin, group_end, group_tag)))
            throw NativeIntegrationError(
                NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
                "negative branch-sensitive fused line has no derivable imaginary sign");
        });
  }

  for_each_exact_group(
      [&](std::size_t group_begin, std::size_t group_end,
          const SectorMonomialTag& group_tag) {
    auto group = build_group(group_begin, group_end, group_tag);
    result.diagnostics.input_monomial_cells += group.cell_count;

    const auto first_nonzero =
        first_nonzero_power(projected_shape, group, 0);
    const auto material_components = first_nonzero.has_value() ? 1U : 0U;
    const bool divergent = center_divergent(
        group.tag, has_center_endpoint);
    const bool cancelled_divergent = material_components != 0 && divergent;
    const bool deferred_divergent =
        cancelled_divergent && options.defer_divergent_groups;
    if (cancelled_divergent && !deferred_divergent)
      (void)accept_or_throw_divergent_group(
          projected_shape, group, options, &result.diagnostics);

    const auto primitive_min = primitive_min_power(
        group.tag, has_center_endpoint);
    const auto deliverable_max = local_detail::checked_i32(
        static_cast<std::int64_t>(projected_epsilon.complete_max) +
            primitive_min,
        "fused line deliverable epsilon maximum");
    if (deliverable_max < delivered_epsilon.complete_max)
      throw NativeIntegrationError(
          NativeIntegrationErrorCode::IncompleteEpsilonWindow, "E10",
          "fused row coefficient frame does not cover the requested primitive output");

    if (deferred_divergent) {
      std::vector<ComplexBall> coefficients;
      coefficients.reserve(group.coefficients.size());
      for (const auto& coefficient : group.coefficients)
        coefficients.push_back(as_ball(coefficient));
      result.deferred_divergent_groups.push_back(
          StoredLineIntegral::DeferredDivergentGroup{
              group.tag,
              EpsilonFrame<ComplexBall>(
                  projected_epsilon, std::move(coefficients)),
              group.contribution_scale_uppers,
              group.cell_count,
              group.had_material_input});
      ++result.diagnostics.zero_groups_skipped;
      return;
    }
    if (material_components == 0 || cancelled_divergent) {
      ++result.diagnostics.zero_groups_skipped;
      if (cancelled_divergent ||
          (divergent && group.had_material_input && group.cell_count > 1))
        ++result.diagnostics.cancelled_divergent_groups;
      return;
    }
    MonomialIntegrationOptions primitive_options;
    primitive_options.imaginary_sign = effective_sign;
    primitive_options.complete_max = std::max(
        primitive_min,
        checked_difference(delivered_epsilon.complete_max,
                           *first_nonzero,
                           "fused required primitive complete maximum"));
    const auto primitive = integrate_sector_monomial(
        group.tag, lower, upper, primitive_options);
    if (primitive.min_power() != primitive_min)
      throw std::logic_error(
          "fused line primitive minimum disagrees with its exact preflight");
    ++result.diagnostics.primitive_evaluations;
    ++result.diagnostics.primitive_component_applications;

    std::vector<ComplexBall> coefficient_values;
    coefficient_values.reserve(static_cast<std::size_t>(
        static_cast<std::int64_t>(projected_epsilon.complete_max) -
        *first_nonzero + 1));
    for (std::int64_t raw_power = *first_nonzero;
         raw_power <= projected_epsilon.complete_max; ++raw_power)
      coefficient_values.push_back(as_ball(group.coefficients[
          static_cast<std::size_t>(raw_power - projected_epsilon.min_power)]));
    const EpsilonFrame<ComplexBall> coefficient_frame(
        *first_nonzero, std::move(coefficient_values));
    const auto contribution = coefficient_frame * primitive;
    if (contribution.complete_max() <
        delivered_epsilon.complete_max)
      throw std::logic_error(
          "fused line epsilon preflight did not deliver its promised window");
    for (std::int64_t raw_power = delivered_epsilon.min_power;
         raw_power <= delivered_epsilon.complete_max; ++raw_power)
      result.value.at(static_cast<std::int32_t>(raw_power), 0) +=
          contribution.coefficient(static_cast<std::int32_t>(raw_power));
  });

  result.value.error.guarantee = ErrorGuarantee::None;
  result.value.error.provenance =
      "stored Taylor truncation only; fused transport row projection; no "
      "unseen-tail majorant or full-local certificate";
  result.diagnostics.detail = result.value.error.provenance;
  return result;
}

struct FactorizedOrdinaryStoredRowIntegral {
  bool eligible = false;
  std::string reason;
  StoredLineIntegral integral;
  std::size_t operator_columns = 0;
};

// Integrate a finite ordinary Taylor recurrence as a transfer operator from
// its retained center value.  This is the line-integral analogue of
// evaluate_ordinary_center_value_factorized: exact unit center impulses are
// evolved and integrated first, then each original center-value ball is
// applied once.  The caller remains responsible for certifying that `source`
// is an authoritative replay of the retained physical recurrence and for
// checking overlap with the historical direct integral.
inline FactorizedOrdinaryStoredRowIntegral
integrate_ordinary_center_stored_row_factorized(
    const PreparedPhysicalClearedODE<ComplexBall>& equation,
    const LocalSolution<ComplexBall>& source,
    const PreparedSparseLocalMultiplierMatrix<ComplexBall>& matrix,
    std::int32_t projected_complete_cap,
    const RealEvaluationPoint& lower,
    const RealEvaluationPoint& upper,
    const StoredLineIntegrationOptions& options,
    const StoredLineIntegral& direct,
    std::size_t maximum_operator_columns =
        std::numeric_limits<std::size_t>::max()) {
  physical_ode_detail::validate_ode(equation);
  validate_local_solution(source, true);
  FactorizedOrdinaryStoredRowIntegral result;
  const auto ineligible = [&](std::string reason) {
    result.reason = std::move(reason);
    result.integral = StoredLineIntegral{};
    result.operator_columns = 0;
    return result;
  };
  if (source.dimension != equation.dimension ||
      source.sectors.size() != 1 ||
      source.sectors.front().a.is_zero != TruthValue::Yes ||
      source.sectors.front().b.is_zero != TruthValue::Yes ||
      source.sectors.front().log_power != 0)
    return ineligible(
        "factorized ordinary row integration requires one ordinary sector "
        "bound to the same physical equation dimension");
  if (!direct.deferred_divergent_groups.empty())
    return ineligible(
        "factorized ordinary row integration cannot reassociate deferred "
        "divergent groups");
  if (source.epsilon.width() >
      std::numeric_limits<std::size_t>::max() / source.dimension)
    throw std::overflow_error(
        "factorized ordinary row integration column count overflows size_t");
  const auto column_count =
      source.epsilon.width() * source.dimension;
  if (column_count > maximum_operator_columns)
    return ineligible(
        "factorized ordinary row integration exceeds its operator-column "
        "cap; required=" + std::to_string(column_count) +
        "; cap=" + std::to_string(maximum_operator_columns));

  EpsilonVector initial =
      physical_ode_detail::zero_epsilon_vector(
          source.epsilon, source.dimension);
  const auto& sector = source.sectors.front();
  for (std::size_t epsilon_index = 0;
       epsilon_index < source.epsilon.width(); ++epsilon_index)
    for (std::uint32_t component = 0;
         component < source.dimension; ++component)
      initial.coefficients[
          epsilon_index * source.dimension + component] =
          sector.coefficients[local_detail::sector_index(
              source, epsilon_index, 0, component)];

  auto factorized = direct;
  std::fill(factorized.value.coefficients.begin(),
            factorized.value.coefficients.end(), ComplexBall(0));
  std::size_t column = 0;
  for (std::int64_t raw_input_power = source.epsilon.min_power;
       raw_input_power <= source.epsilon.complete_max;
       ++raw_input_power) {
    const auto input_power =
        static_cast<std::int32_t>(raw_input_power);
    for (std::uint32_t input_component = 0;
         input_component < source.dimension; ++input_component) {
      auto impulse = physical_ode_detail::zero_epsilon_vector(
          source.epsilon, source.dimension);
      impulse.at(input_power, input_component) = ComplexBall(1);
      const auto evolved = evolve_ordinary_center_value<ComplexBall>(
          equation, impulse, source.taylor_complete_max);
      if (!evolved.eligible)
        return ineligible(
            "factorized row impulse evolution is ineligible: " +
            evolved.reason);
      auto impulse_local = ordinary_evolution_local_solution(
          evolved, source.chart, source.prescriptions,
          source.checkpoint_identity +
              ":factorized-ordinary-row-column:" +
              std::to_string(column));
      StoredLineIntegral response;
      try {
        response = integrate_prepared_scalar_row_stored(
            matrix, impulse_local, projected_complete_cap,
            lower, upper, options);
      } catch (const NativeIntegrationError& error) {
        return ineligible(
            "factorized row impulse integration is unsupported: " +
            std::string(error.what()));
      }
      if (!response.deferred_divergent_groups.empty() ||
          response.diagnostics.cancelled_divergent_groups != 0 ||
          response.diagnostics
                  .bounded_cancelled_divergent_coefficients != 0)
        return ineligible(
            "factorized row impulse integration encountered a divergent "
            "group whose cancellation is not column-separable");
      if (response.value.dimension != factorized.value.dimension ||
          response.value.epsilon.min_power !=
              factorized.value.epsilon.min_power ||
          response.value.epsilon.complete_max !=
              factorized.value.epsilon.complete_max ||
          response.scope != factorized.scope ||
          response.imaginary_sign != factorized.imaginary_sign)
        throw std::logic_error(
            "factorized ordinary row response changed the direct integral "
            "contract");
      const auto& amplitude =
          initial.at(input_power, input_component);
      for (std::int64_t raw_output_power =
               factorized.value.epsilon.min_power;
           raw_output_power <=
               factorized.value.epsilon.complete_max;
           ++raw_output_power) {
        const auto output_power =
            static_cast<std::int32_t>(raw_output_power);
        for (std::uint32_t output_component = 0;
             output_component < factorized.value.dimension;
             ++output_component)
          factorized.value.at(output_power, output_component) +=
              local_detail::coefficient_or_zero(
                  response.value, output_power,
                  output_component) *
              amplitude;
      }
      ++column;
    }
  }
  result.eligible = true;
  result.reason =
      "finite ordinary stored-row transfer operator applied once to "
      "retained center-value balls";
  result.operator_columns = column_count;
  result.integral = std::move(factorized);
  return result;
}

}  // namespace diffexp2
