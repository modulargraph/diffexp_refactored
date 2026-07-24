#pragma once

#include "diffexp2/scalar.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace diffexp2 {

// The planner is deliberately an exact, protocol-independent geometry
// layer.  It receives affine charts whose physical coordinate is
//
//                     x = center + scale * t
//
// and returns the exact physical/local intervals consumed by matching and
// stored line integration.  JSON/session ownership belongs to the persistent
// driver, not to this header.

enum class ExactPathPlanningErrorCode : std::uint8_t {
  InvalidInput,
  InvalidTopology,
  UnsafeGeometry,
  ForbiddenMatch,
  InvalidPlan
};

class ExactPathPlanningError : public std::runtime_error {
 public:
  ExactPathPlanningError(ExactPathPlanningErrorCode code,
                         const std::string& detail)
      : std::runtime_error(detail), code(code) {}

  ExactPathPlanningErrorCode code;
};

// A branch sheet is input topology.  The planner never guesses a rim from a
// floating imaginary part.  Multiple factor identities may pin the same rim;
// conflicting duplicates are rejected exactly.
struct ExactBranchSheet {
  std::string factor_exact;
  std::int32_t imaginary_sign = 0;

  friend bool operator==(const ExactBranchSheet&,
                         const ExactBranchSheet&) = default;
};

// Faithful metadata for the classic real projection of one complex
// singularity.  Retention/suppression is decided by the symbolic producer and
// carried here as topology; the native planner does not redo a numeric
// occupancy test.  The three possible waypoints are Re(z), Re(z)-|Im(z)| and
// Re(z)+|Im(z)|.
struct ExactComplexProjection {
  std::string source_identity;
  Rational real_part;
  Rational imaginary_magnitude;
  bool retain_minus_imaginary = true;
  bool retain_real_part = true;
  bool retain_plus_imaginary = true;

  [[nodiscard]] Rational minus_imaginary() const {
    return real_part - imaginary_magnitude;
  }
  [[nodiscard]] Rational plus_imaginary() const {
    return real_part + imaginary_magnitude;
  }
  [[nodiscard]] std::vector<Rational> retained_waypoints() const {
    std::vector<Rational> out;
    if (retain_minus_imaginary) out.push_back(minus_imaginary());
    if (retain_real_part) out.push_back(real_part);
    if (retain_plus_imaginary) out.push_back(plus_imaginary());
    return out;
  }
};

struct ExactPathTopology {
  // True real singularities.  Every one in the open arm must be represented
  // by an explicitly singular chart; none may be used as a handoff point.
  std::vector<Rational> singular_points;
  // Additional exact points that may be chart centers or integration split
  // points but must never be selected as matches.
  std::vector<Rational> boundary_points;
  std::vector<ExactComplexProjection> complex_projections;
  std::vector<ExactBranchSheet> branch_sheets;
};

struct ExactAffineChart {
  std::string identity;
  Rational center;
  Rational scale = Rational(1);
  // True convergence radius in the local coordinate t.  The conditioning
  // scale is |scale|: an ordinary classic handoff targets |t|=1/k and is
  // additionally capped at half this true radius.
  Rational radius = Rational(1);
  bool singular_center = false;
};

struct ExactArmRequest {
  Rational from;
  Rational to;
  std::vector<ExactAffineChart> charts;
  ExactPathTopology topology;
};

struct ExactPathPlanOptions {
  std::uint32_t division_order = 3;
};

enum class ExactMatchKind : std::uint8_t {
  // The legacy GetCPL/GetCPR invariant: exactly +1/k in the producing chart
  // and -1/k in the receiving chart, with signs reversed on a reversed arm.
  SymmetricDivisionPoint,
  // Unequal/clipped geometry: the common point balances the two exact safe
  // reaches while retaining the same <=1/k and <=radius/2 envelope.
  BalancedSafeOverlap,
  // A singular receiving chart uses the exact FixWithin balance retained by
  // the Wolfram planner.  Both charts remain in their affine 1/k/half-radius
  // envelopes, and the receiver is approached from the physical near side.
  SingularBalancedApproach,
  // The balanced/symmetric candidate was a declared singular/boundary point;
  // a deterministic exact interior rational was selected instead.
  ForbiddenPointAvoidance
};

struct ExactMatchPoint {
  std::size_t producing_chart = 0;
  std::size_t receiving_chart = 0;
  Rational physical;
  Rational producing_local;
  Rational receiving_local;
  ExactMatchKind kind = ExactMatchKind::BalancedSafeOverlap;
  std::vector<ExactBranchSheet> branch_sheets;
};

struct ExactPathTile {
  std::size_t chart = 0;
  Rational physical_begin;
  Rational physical_end;
  Rational local_begin;
  Rational local_end;
  bool crosses_singular_center = false;
  std::vector<ExactBranchSheet> branch_sheets;
};

struct ExactArmPlan {
  Rational from;
  Rational to;
  std::int32_t direction = 0;
  std::uint32_t division_order = 3;
  std::vector<ExactAffineChart> charts;
  std::vector<ExactMatchPoint> matches;
  std::vector<ExactPathTile> tiles;
  ExactPathTopology topology;
};

// Two immutable value plans sharing one exact anchor chart.  They contain no
// mutable solver/session state and can therefore be submitted to independent
// workers by the persistent driver.
struct ExactIndependentArmPlans {
  ExactAffineChart anchor_chart;
  ExactArmPlan lower;
  ExactArmPlan upper;
};

namespace exact_path_detail {

inline Rational abs(const Rational& value) {
  return value.sign() < 0 ? -value : value;
}

inline Rational minimum(const Rational& left, const Rational& right) {
  return left < right ? left : right;
}

inline Rational maximum(const Rational& left, const Rational& right) {
  return left > right ? left : right;
}

inline Rational directed(const Rational& value, std::int32_t direction) {
  return direction < 0 ? -value : value;
}

inline Rational local_coordinate(const ExactAffineChart& chart,
                                 const Rational& physical) {
  return (physical - chart.center) / chart.scale;
}

inline Rational physical_coordinate(const ExactAffineChart& chart,
                                    const Rational& local) {
  return chart.center + chart.scale * local;
}

inline bool same_chart_geometry(const ExactAffineChart& left,
                                const ExactAffineChart& right) {
  return left.identity == right.identity && left.center == right.center &&
         left.scale == right.scale && left.radius == right.radius &&
         left.singular_center == right.singular_center;
}

inline bool contains(const std::vector<Rational>& points,
                     const Rational& point) {
  return std::any_of(points.begin(), points.end(),
                     [&](const Rational& candidate) {
                       return candidate == point;
                     });
}

inline bool same_branch_sheets(const std::vector<ExactBranchSheet>& left,
                               const std::vector<ExactBranchSheet>& right) {
  return left == right;
}

inline void validate_topology(const ExactPathTopology& topology) {
  std::unordered_set<std::string> projection_ids;
  for (const auto& projection : topology.complex_projections) {
    if (projection.source_identity.empty() ||
        projection.imaginary_magnitude.sign() <= 0 ||
        (!projection.retain_minus_imaginary &&
         !projection.retain_real_part &&
         !projection.retain_plus_imaginary))
      throw ExactPathPlanningError(
          ExactPathPlanningErrorCode::InvalidTopology,
          "complex projection requires a unique identity, a positive exact "
          "imaginary magnitude, and at least one retained waypoint");
    if (!projection_ids.insert(projection.source_identity).second)
      throw ExactPathPlanningError(
          ExactPathPlanningErrorCode::InvalidTopology,
          "duplicate complex-projection identity");
  }

  std::vector<ExactBranchSheet> seen;
  for (const auto& sheet : topology.branch_sheets) {
    if (sheet.factor_exact.empty() ||
        (sheet.imaginary_sign != -1 && sheet.imaginary_sign != 1))
      throw ExactPathPlanningError(
          ExactPathPlanningErrorCode::InvalidTopology,
          "branch sheet requires an exact factor identity and sign +/-1");
    const auto found = std::find_if(
        seen.begin(), seen.end(), [&](const ExactBranchSheet& candidate) {
          return candidate.factor_exact == sheet.factor_exact;
        });
    if (found != seen.end()) {
      if (found->imaginary_sign != sheet.imaginary_sign)
        throw ExactPathPlanningError(
            ExactPathPlanningErrorCode::InvalidTopology,
            "conflicting exact branch signs for one factor identity");
    } else {
      seen.push_back(sheet);
    }
  }
}

inline void validate_request_geometry(const ExactArmRequest& request,
                                      std::int32_t direction) {
  if (!(request.charts.front().center == request.from))
    throw ExactPathPlanningError(
        ExactPathPlanningErrorCode::InvalidInput,
        "the first chart must be the exact boundary anchor");
  const auto arm_begin = directed(request.from, direction);
  const auto arm_end = directed(request.to, direction);
  for (std::size_t index = 0; index < request.charts.size(); ++index) {
    const auto& chart = request.charts[index];
    const auto center = directed(chart.center, direction);
    if (chart.identity.empty() || chart.scale.is_zero() ||
        chart.radius.sign() <= 0 || center < arm_begin || center > arm_end ||
        (index > 0 &&
         !(center > directed(request.charts[index - 1].center,
                             direction))))
      throw ExactPathPlanningError(
          ExactPathPlanningErrorCode::InvalidInput,
          "chart identity/scale/radius/order is invalid");
    if (contains(request.topology.singular_points, chart.center) !=
        chart.singular_center)
      throw ExactPathPlanningError(
          ExactPathPlanningErrorCode::InvalidTopology,
          "chart singular flag disagrees with exact path topology");
  }
  for (const auto& singularity : request.topology.singular_points) {
    const auto coordinate = directed(singularity, direction);
    if (coordinate < arm_begin || coordinate > arm_end) continue;
    if (std::none_of(request.charts.begin(), request.charts.end(),
                     [&](const ExactAffineChart& chart) {
                       return chart.singular_center &&
                              chart.center == singularity;
                     }))
      throw ExactPathPlanningError(
          ExactPathPlanningErrorCode::InvalidTopology,
          "an on-arm real singularity has no explicit singular chart");
  }
}

inline Rational safe_local_limit(const ExactAffineChart& chart,
                                 std::uint32_t division_order) {
  const Rational k(std::to_string(division_order));
  return minimum(Rational(1) / k, chart.radius / Rational(2));
}

inline Rational safe_physical_reach(const ExactAffineChart& chart,
                                    std::uint32_t division_order) {
  return abs(chart.scale) * safe_local_limit(chart, division_order);
}

inline Rational physical_radius(const ExactAffineChart& chart) {
  return abs(chart.scale) * chart.radius;
}

// Exact rational counterpart of Transport`singularMatchPoint.  Balance the
// handoff using AFFINE safe reaches, not the charts' full analytic physical
// radii.  A singular chart can have a large true radius but a tiny transport
// scale; using the former here placed endpoint reads near |t|=1 and amplified
// otherwise harmless Taylor tails by hundreds of orders of magnitude.
//
// The 9/10 producer margin leaves exact headroom for an algebraic chart whose
// scale is subsequently replaced by a strict inward rational lower bound.
inline Rational singular_balanced_approach(
    const ExactAffineChart& left, const ExactAffineChart& right,
    std::uint32_t division_order, std::int32_t direction) {
  const auto y_left = directed(left.center, direction);
  const auto y_right = directed(right.center, direction);
  const auto gap = y_right - y_left;
  const Rational margin("9/10");
  const auto left_reach =
      margin * safe_physical_reach(left, division_order);
  const auto right_reach = safe_physical_reach(right, division_order);
  const auto denominator = left_reach + right_reach;
  if (!(gap > Rational(0)) || !(gap < denominator))
    throw ExactPathPlanningError(
        ExactPathPlanningErrorCode::UnsafeGeometry,
        "singular receiving chart has no conditioned affine overlap");
  return directed(y_left + gap * left_reach / denominator, direction);
}

inline bool at_or_inside_safe_envelope(const ExactAffineChart& chart,
                                       const Rational& local,
                                       std::uint32_t division_order);

inline std::vector<Rational> forbidden_matches(const ExactArmRequest& request) {
  std::vector<Rational> out = request.topology.singular_points;
  out.insert(out.end(), request.topology.boundary_points.begin(),
             request.topology.boundary_points.end());
  out.push_back(request.from);
  out.push_back(request.to);
  return out;
}

inline Rational choose_nonforbidden_interior(
    const Rational& lower, const Rational& upper,
    const std::vector<Rational>& forbidden) {
  if (!(lower < upper))
    throw ExactPathPlanningError(
        ExactPathPlanningErrorCode::ForbiddenMatch,
        "the only common safe handoff is a forbidden singular/boundary point");

  // There are F forbidden values and F+1 distinct interior candidates, so
  // one must survive.  This is deterministic and remains entirely in Q.
  const auto candidate_count = forbidden.size() + 1;
  const Rational denominator(std::to_string(candidate_count + 1));
  for (std::size_t index = 1; index <= candidate_count; ++index) {
    const auto candidate =
        lower + (upper - lower) * Rational(std::to_string(index)) /
                    denominator;
    if (!contains(forbidden, candidate)) return candidate;
  }
  throw ExactPathPlanningError(
      ExactPathPlanningErrorCode::ForbiddenMatch,
      "could not select a nonforbidden exact handoff from a nonempty interval");
}

inline ExactMatchPoint plan_match(const ExactArmRequest& request,
                                  std::size_t left_index,
                                  std::uint32_t division_order,
                                  std::int32_t direction) {
  const auto& left = request.charts[left_index];
  const auto& right = request.charts[left_index + 1];
  if (left.scale.sign() != right.scale.sign())
    throw ExactPathPlanningError(
        ExactPathPlanningErrorCode::UnsafeGeometry,
        "adjacent affine chart scales have opposite orientation, so an "
        "inter-center match cannot have opposite local signs");

  const Rational k(std::to_string(division_order));
  const Rational local_direction(direction * left.scale.sign());
  const auto preferred_left_local = local_direction / k;
  const auto preferred_right_local = -local_direction / k;
  const auto preferred_from_left =
      physical_coordinate(left, preferred_left_local);
  const auto preferred_from_right =
      physical_coordinate(right, preferred_right_local);
  const auto forbidden = forbidden_matches(request);

  ExactMatchPoint result;
  result.producing_chart = left_index;
  result.receiving_chart = left_index + 1;
  result.branch_sheets = request.topology.branch_sheets;
  if (right.singular_center) {
    if (left.singular_center)
      throw ExactPathPlanningError(
          ExactPathPlanningErrorCode::UnsafeGeometry,
          "a singular receiving chart requires a regular producing chart");
    result.physical = singular_balanced_approach(
        left, right, division_order, direction);
    result.producing_local = local_coordinate(left, result.physical);
    result.receiving_local = local_coordinate(right, result.physical);
    if (!at_or_inside_safe_envelope(left, result.producing_local,
                                    division_order) ||
        !at_or_inside_safe_envelope(right, result.receiving_local,
                                    division_order) ||
        contains(forbidden, result.physical))
      throw ExactPathPlanningError(
          ExactPathPlanningErrorCode::UnsafeGeometry,
          "canonical singular approach lies outside an affine 1/k or "
          "half-radius envelope");
    result.kind = ExactMatchKind::SingularBalancedApproach;
    return result;
  }
  if (preferred_from_left == preferred_from_right &&
      abs(preferred_left_local) <= left.radius / Rational(2) &&
      abs(preferred_right_local) <= right.radius / Rational(2) &&
      !contains(forbidden, preferred_from_left)) {
    result.physical = preferred_from_left;
    result.producing_local = preferred_left_local;
    result.receiving_local = preferred_right_local;
    result.kind = ExactMatchKind::SymmetricDivisionPoint;
    return result;
  }

  const auto y_left = directed(left.center, direction);
  const auto y_right = directed(right.center, direction);
  // SegmentLine gives the receiving chart ownership of an ordinary
  // handoff: its canonical point is -1/k in the receiving coordinate, while
  // the producer only has to remain inside the certified half-radius error
  // envelope.  Requiring 1/k on both sides rejects valid asymmetric chains
  // after exact rational scale bridging and makes increasing DivisionOrder
  // paradoxically reduce planner clearance.
  const auto left_reach =
      abs(left.scale) * left.radius / Rational(2);
  const auto right_reach = safe_physical_reach(right, division_order);
  const auto preferred_right_y = directed(preferred_from_right, direction);
  const auto preferred_left_at_right =
      local_coordinate(left, preferred_from_right);
  if (preferred_right_y > y_left && preferred_right_y < y_right &&
      abs(preferred_left_at_right) <= left.radius / Rational(2) &&
      abs(preferred_right_local) <=
          safe_local_limit(right, division_order) &&
      !contains(forbidden, preferred_from_right)) {
    result.physical = preferred_from_right;
    result.producing_local = preferred_left_at_right;
    result.receiving_local = preferred_right_local;
    result.kind = ExactMatchKind::BalancedSafeOverlap;
    return result;
  }
  const auto gap = y_right - y_left;
  if (!(gap > Rational(0)) || gap > left_reach + right_reach)
    throw ExactPathPlanningError(
        ExactPathPlanningErrorCode::UnsafeGeometry,
        "adjacent charts have no common exact handoff inside the producer "
        "half-radius and receiving 1/k envelopes; pair=" +
        std::to_string(left_index) +
        "; left=" + left.identity + "; right=" + right.identity +
        "; gap=" + gap.str() + "; left_reach=" + left_reach.str() +
        "; right_reach=" + right_reach.str() +
        "; division_order=" + std::to_string(division_order));

  auto y_match = y_left + gap * left_reach / (left_reach + right_reach);
  auto physical = directed(y_match, direction);
  result.kind = ExactMatchKind::BalancedSafeOverlap;
  if (contains(forbidden, physical)) {
    const auto lower = maximum(y_left, y_right - right_reach);
    const auto upper = minimum(y_right, y_left + left_reach);
    const auto y_alternative = choose_nonforbidden_interior(
        lower, upper,
        [&]() {
          std::vector<Rational> directed_forbidden;
          directed_forbidden.reserve(forbidden.size());
          for (const auto& point : forbidden)
            directed_forbidden.push_back(directed(point, direction));
          return directed_forbidden;
        }());
    physical = directed(y_alternative, direction);
    result.kind = ExactMatchKind::ForbiddenPointAvoidance;
  }
  result.physical = physical;
  result.producing_local = local_coordinate(left, physical);
  result.receiving_local = local_coordinate(right, physical);
  return result;
}

inline bool strictly_inside(const ExactAffineChart& chart,
                            const Rational& local) {
  return abs(local) < chart.radius;
}

inline bool at_or_inside_safe_envelope(const ExactAffineChart& chart,
                                       const Rational& local,
                                       std::uint32_t division_order) {
  return abs(local) <= safe_local_limit(chart, division_order);
}

inline bool valid_match_envelopes(const ExactAffineChart& left,
                                  const ExactAffineChart& right,
                                  const ExactMatchPoint& match,
                                  std::uint32_t division_order,
                                  std::int32_t direction) {
  const bool singular_approach =
      match.kind == ExactMatchKind::SingularBalancedApproach;
  if (right.singular_center != singular_approach) return false;
  if (!singular_approach)
    return abs(match.producing_local) <= left.radius / Rational(2) &&
           at_or_inside_safe_envelope(right, match.receiving_local,
                                      division_order);
  if (left.singular_center || !right.singular_center ||
      !at_or_inside_safe_envelope(left, match.producing_local,
                                  division_order) ||
      !at_or_inside_safe_envelope(right, match.receiving_local,
                                  division_order))
    return false;
  try {
    return match.physical ==
           singular_balanced_approach(left, right, division_order, direction);
  } catch (const ExactPathPlanningError&) {
    return false;
  }
}

}  // namespace exact_path_detail

inline void validate_exact_arm_plan(const ExactArmPlan& plan) {
  using namespace exact_path_detail;
  if (plan.direction != -1 && plan.direction != 1)
    throw ExactPathPlanningError(ExactPathPlanningErrorCode::InvalidPlan,
                                 "plan direction is not +/-1");
  if (plan.division_order < 2 || plan.charts.empty() ||
      plan.matches.size() + 1 != plan.charts.size() ||
      plan.tiles.size() != plan.charts.size())
    throw ExactPathPlanningError(
        ExactPathPlanningErrorCode::InvalidPlan,
        "plan has an invalid division order or inconsistent chart/match/tile counts");
  if (!(directed(plan.to, plan.direction) >
        directed(plan.from, plan.direction)) ||
      !(plan.charts.front().center == plan.from))
    throw ExactPathPlanningError(
        ExactPathPlanningErrorCode::InvalidPlan,
        "plan is empty/reversed relative to its direction or lacks the exact anchor chart");
  validate_topology(plan.topology);

  const auto arm_begin = directed(plan.from, plan.direction);
  const auto arm_end = directed(plan.to, plan.direction);
  for (std::size_t index = 0; index < plan.charts.size(); ++index) {
    const auto& chart = plan.charts[index];
    const auto center = directed(chart.center, plan.direction);
    if (chart.identity.empty() || chart.scale.is_zero() ||
        chart.radius.sign() <= 0 || center < arm_begin || center > arm_end ||
        (index > 0 &&
         !(center > directed(plan.charts[index - 1].center,
                             plan.direction))))
      throw ExactPathPlanningError(
          ExactPathPlanningErrorCode::InvalidPlan,
          "chart identity/scale/radius/order is invalid");
    const bool declared_singular = contains(plan.topology.singular_points,
                                            chart.center);
    if (declared_singular != chart.singular_center)
      throw ExactPathPlanningError(
          ExactPathPlanningErrorCode::InvalidPlan,
          "chart singular flag disagrees with exact path topology");
  }
  for (const auto& singularity : plan.topology.singular_points) {
    const auto y = directed(singularity, plan.direction);
    if (y < arm_begin || y > arm_end) continue;
    if (std::none_of(plan.charts.begin(), plan.charts.end(),
                     [&](const ExactAffineChart& chart) {
                       return chart.singular_center &&
                              chart.center == singularity;
                     }))
      throw ExactPathPlanningError(
          ExactPathPlanningErrorCode::InvalidPlan,
          "an on-arm real singularity has no explicit singular chart");
  }

  const std::vector<Rational> forbidden = [&]() {
    std::vector<Rational> points = plan.topology.singular_points;
    points.insert(points.end(), plan.topology.boundary_points.begin(),
                  plan.topology.boundary_points.end());
    points.push_back(plan.from);
    points.push_back(plan.to);
    return points;
  }();
  Rational previous_physical = plan.from;
  for (std::size_t index = 0; index < plan.matches.size(); ++index) {
    const auto& match = plan.matches[index];
    const auto& left = plan.charts[index];
    const auto& right = plan.charts[index + 1];
    const auto y = directed(match.physical, plan.direction);
    if (match.producing_chart != index ||
        match.receiving_chart != index + 1 ||
        !(y > directed(left.center, plan.direction)) ||
        !(y < directed(right.center, plan.direction)) ||
        !(directed(match.physical, plan.direction) >
          directed(previous_physical, plan.direction)) ||
        contains(forbidden, match.physical) ||
        !(physical_coordinate(left, match.producing_local) ==
          match.physical) ||
        !(physical_coordinate(right, match.receiving_local) ==
          match.physical) ||
        match.producing_local.sign() * match.receiving_local.sign() >= 0 ||
        !strictly_inside(left, match.producing_local) ||
        !strictly_inside(right, match.receiving_local) ||
        !valid_match_envelopes(left, right, match, plan.division_order,
                               plan.direction) ||
        !same_branch_sheets(match.branch_sheets,
                            plan.topology.branch_sheets))
      throw ExactPathPlanningError(
          ExactPathPlanningErrorCode::InvalidPlan,
          "match violates exact identity, ordering, topology, opposite-sign, or safe-radius invariants");
    previous_physical = match.physical;
  }

  for (std::size_t index = 0; index < plan.tiles.size(); ++index) {
    const auto& tile = plan.tiles[index];
    const auto& chart = plan.charts[index];
    const auto expected_begin =
        index == 0 ? plan.from : plan.matches[index - 1].physical;
    const auto expected_end = index + 1 == plan.charts.size()
                                  ? plan.to
                                  : plan.matches[index].physical;
    const bool expected_crossing = chart.singular_center &&
        directed(expected_begin, plan.direction) <
            directed(chart.center, plan.direction) &&
        directed(expected_end, plan.direction) >
            directed(chart.center, plan.direction);
    if (tile.chart != index || !(tile.physical_begin == expected_begin) ||
        !(tile.physical_end == expected_end) ||
        !(directed(tile.physical_end, plan.direction) >=
          directed(tile.physical_begin, plan.direction)) ||
        !(physical_coordinate(chart, tile.local_begin) ==
          tile.physical_begin) ||
        !(physical_coordinate(chart, tile.local_end) == tile.physical_end) ||
        !strictly_inside(chart, tile.local_begin) ||
        !strictly_inside(chart, tile.local_end) ||
        tile.crosses_singular_center != expected_crossing ||
        !same_branch_sheets(tile.branch_sheets,
                            plan.topology.branch_sheets))
      throw ExactPathPlanningError(
          ExactPathPlanningErrorCode::InvalidPlan,
          "tile chain has a gap/loop, unstable local identity, invalid radius, or lost topology metadata");
  }
  if (!(plan.tiles.front().physical_begin == plan.from) ||
      !(plan.tiles.back().physical_end == plan.to))
    throw ExactPathPlanningError(
        ExactPathPlanningErrorCode::InvalidPlan,
        "tile chain does not cover the arm endpoints exactly");
}

inline ExactArmPlan plan_exact_arm(
    const ExactArmRequest& request,
    const ExactPathPlanOptions& options = ExactPathPlanOptions{}) {
  using namespace exact_path_detail;
  if (options.division_order < 2 || request.from == request.to ||
      request.charts.empty())
    throw ExactPathPlanningError(
        ExactPathPlanningErrorCode::InvalidInput,
        "exact arm requires k>=2, distinct endpoints, and at least one chart");
  validate_topology(request.topology);
  const auto direction = request.to > request.from ? 1 : -1;
  validate_request_geometry(request, direction);

  ExactArmPlan plan;
  plan.from = request.from;
  plan.to = request.to;
  plan.direction = direction;
  plan.division_order = options.division_order;
  plan.charts = request.charts;
  plan.topology = request.topology;

  plan.matches.reserve(plan.charts.size() - 1);
  for (std::size_t index = 0; index + 1 < plan.charts.size(); ++index)
    plan.matches.push_back(
        plan_match(request, index, options.division_order, direction));

  plan.tiles.reserve(plan.charts.size());
  for (std::size_t index = 0; index < plan.charts.size(); ++index) {
    const auto begin = index == 0 ? request.from
                                  : plan.matches[index - 1].physical;
    const auto end = index + 1 == plan.charts.size()
                         ? request.to
                         : plan.matches[index].physical;
    const auto& chart = plan.charts[index];
    ExactPathTile tile;
    tile.chart = index;
    tile.physical_begin = begin;
    tile.physical_end = end;
    tile.local_begin = local_coordinate(chart, begin);
    tile.local_end = local_coordinate(chart, end);
    tile.crosses_singular_center = chart.singular_center &&
        directed(begin, direction) < directed(chart.center, direction) &&
        directed(end, direction) > directed(chart.center, direction);
    tile.branch_sheets = request.topology.branch_sheets;
    plan.tiles.push_back(std::move(tile));
  }
  validate_exact_arm_plan(plan);
  return plan;
}

inline ExactIndependentArmPlans plan_exact_independent_arms(
    const ExactArmRequest& lower, const ExactArmRequest& upper,
    const ExactPathPlanOptions& options = ExactPathPlanOptions{}) {
  if (!(lower.from == upper.from) || lower.charts.empty() ||
      upper.charts.empty() ||
      !exact_path_detail::same_chart_geometry(lower.charts.front(),
                                              upper.charts.front()))
    throw ExactPathPlanningError(
        ExactPathPlanningErrorCode::InvalidInput,
        "independent arms require one stable exact anchor chart");
  const auto lower_direction = lower.to > lower.from ? 1 : -1;
  const auto upper_direction = upper.to > upper.from ? 1 : -1;
  if (lower_direction == upper_direction)
    throw ExactPathPlanningError(
        ExactPathPlanningErrorCode::InvalidInput,
        "lower and upper arm endpoints must lie on opposite sides of the anchor");

  ExactIndependentArmPlans result;
  result.anchor_chart = lower.charts.front();
  result.lower = plan_exact_arm(lower, options);
  result.upper = plan_exact_arm(upper, options);
  return result;
}

}  // namespace diffexp2
