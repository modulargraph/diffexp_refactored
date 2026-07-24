#include "diffexp2/path_planner.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using diffexp2::ExactAffineChart;
using diffexp2::ExactArmRequest;
using diffexp2::ExactBranchSheet;
using diffexp2::ExactComplexProjection;
using diffexp2::ExactMatchKind;
using diffexp2::ExactPathPlanOptions;
using diffexp2::ExactPathPlanningError;
using diffexp2::ExactPathTopology;
using diffexp2::Rational;

namespace {

int passed = 0;
int failed = 0;

void check(const std::string& label, bool condition) {
  if (condition) {
    ++passed;
    std::cout << "  PASS: " << label << '\n';
  } else {
    ++failed;
    std::cout << "  FAIL: " << label << '\n';
  }
}

ExactPathTopology topology() {
  ExactPathTopology result;
  result.singular_points = {Rational(1)};
  result.boundary_points = {Rational("7/4")};
  result.complex_projections.push_back(ExactComplexProjection{
      "complex-root-pair", Rational("5/2"), Rational("1/2"), true, true,
      true});
  result.branch_sheets.push_back(
      ExactBranchSheet{"threshold-factor", -1});
  return result;
}

ExactAffineChart left_chart() {
  return ExactAffineChart{"left", Rational(0), Rational(1), Rational(2),
                          false};
}

ExactAffineChart middle_chart() {
  return ExactAffineChart{"middle-singular", Rational(1), Rational(2),
                          Rational(2), true};
}

ExactAffineChart right_chart() {
  return ExactAffineChart{"right", Rational(2), Rational(1),
                          Rational("4/3"), false};
}

void classic_exact_chain_and_reverse() {
  ExactPathPlanOptions options;
  check("DivisionOrder defaults to three", options.division_order == 3);

  ExactAffineChart forward_clearance{"forward-clearance", Rational("1/2"),
                                     Rational(1), Rational(2), false};
  ExactArmRequest forward{
      Rational(0), Rational(2),
      {left_chart(), forward_clearance, middle_chart(), right_chart()},
      topology()};
  const auto plan = diffexp2::plan_exact_arm(forward, options);
  check("singular receiving and ordinary matches retain distinct exact contracts",
        plan.matches.size() == 3 &&
            plan.matches[1].kind ==
                ExactMatchKind::SingularBalancedApproach &&
            plan.matches[1].physical == Rational("19/29") &&
            plan.matches[1].producing_local == Rational("9/58") &&
            plan.matches[1].receiving_local == Rational("-5/29") &&
            plan.matches[2].kind ==
                ExactMatchKind::SymmetricDivisionPoint &&
            plan.matches[2].physical == Rational("5/3") &&
            plan.matches[2].producing_local == Rational("1/3") &&
            plan.matches[2].receiving_local == Rational("-1/3"));
  check("tiles cover exactly once with no gap and flag the singular crossing",
        plan.tiles.size() == 4 &&
            plan.tiles[0].physical_begin == Rational(0) &&
            plan.tiles[0].physical_end == plan.tiles[1].physical_begin &&
            plan.tiles[1].physical_begin == plan.tiles[0].physical_end &&
            plan.tiles[1].physical_end == plan.tiles[2].physical_begin &&
            plan.tiles[2].physical_end == plan.tiles[3].physical_begin &&
            plan.tiles[3].physical_end == Rational(2) &&
            plan.tiles[2].crosses_singular_center);
  check("explicit -i0 topology propagates to every match and tile",
        std::all_of(plan.matches.begin(), plan.matches.end(), [](const auto& m) {
          return m.branch_sheets.size() == 1 &&
                 m.branch_sheets.front().imaginary_sign == -1;
        }) &&
            std::all_of(plan.tiles.begin(), plan.tiles.end(), [](const auto& t) {
              return t.branch_sheets.size() == 1 &&
                     t.branch_sheets.front().imaginary_sign == -1;
            }));
  const auto waypoints =
      plan.topology.complex_projections.front().retained_waypoints();
  check("Re and Re+/-Im projection topology remains exact and uninferred",
        waypoints == std::vector<Rational>({Rational(2), Rational("5/2"),
                                            Rational(3)}));

  ExactAffineChart reverse_clearance{"reverse-clearance", Rational("3/2"),
                                     Rational(1), Rational(2), false};
  ExactArmRequest reverse{
      Rational(2), Rational(0),
      {right_chart(), reverse_clearance, middle_chart(), left_chart()},
      topology()};
  const auto reversed = diffexp2::plan_exact_arm(reverse, options);
  check("reversed arm mirrors physical order and local match signs",
        reversed.direction == -1 && reversed.matches.size() == 3 &&
            reversed.matches[1].kind ==
                ExactMatchKind::SingularBalancedApproach &&
            reversed.matches[1].physical == Rational("39/29") &&
            reversed.matches[1].producing_local == Rational("-9/58") &&
            reversed.matches[1].receiving_local == Rational("5/29") &&
            reversed.matches[2].kind ==
                ExactMatchKind::SymmetricDivisionPoint &&
            reversed.matches[2].physical == Rational("1/3") &&
            reversed.matches[2].producing_local == Rational("-1/3") &&
            reversed.matches[2].receiving_local == Rational("1/3"));
}

void box_bubble_singular_endpoint_handoffs() {
  ExactAffineChart anchor{"box-bubble-anchor", Rational("11/23"),
                          Rational("11/23"), Rational(1), false};
  ExactAffineChart lower_endpoint{"box-bubble-lower-singular", Rational(0),
                                  Rational("22/23"), Rational(1), true};
  ExactAffineChart upper_endpoint{"box-bubble-upper-singular", Rational(1),
                                  Rational(1), Rational(1), true};

  ExactPathTopology lower_topology;
  lower_topology.singular_points = {Rational(0)};
  lower_topology.boundary_points = {Rational(0)};
  lower_topology.branch_sheets = {ExactBranchSheet{"xx1", 1}};
  ExactPathTopology upper_topology;
  upper_topology.singular_points = {Rational(1)};
  upper_topology.boundary_points = {Rational(1)};
  upper_topology.branch_sheets = {ExactBranchSheet{"-1 + xx1", 1}};

  ExactAffineChart lower_clearance{"box-bubble-lower-clearance",
                                   Rational("1/3"), Rational(1),
                                   Rational(1), false};
  ExactAffineChart upper_clearance{"box-bubble-upper-clearance",
                                   Rational("2/3"), Rational(1),
                                   Rational(1), false};
  ExactArmRequest lower{Rational("11/23"), Rational(0),
                        {anchor, lower_clearance, lower_endpoint},
                        lower_topology};
  ExactArmRequest upper{Rational("11/23"), Rational(1),
                        {anchor, upper_clearance, upper_endpoint},
                        upper_topology};
  const auto arms = diffexp2::plan_exact_independent_arms(lower, upper);
  const auto& lower_match = arms.lower.matches.back();
  const auto& upper_match = arms.upper.matches.back();

  check("box-bubble lower singular endpoint uses conditioned affine balance",
        lower_match.kind == ExactMatchKind::SingularBalancedApproach &&
            lower_match.physical == Rational("220/1281") &&
            lower_match.producing_local == Rational("-69/427") &&
            lower_match.receiving_local == Rational("230/1281"));
  check("box-bubble upper singular endpoint stays well inside both affine frames",
        upper_match.kind == ExactMatchKind::SingularBalancedApproach &&
            upper_match.physical == Rational("47/57") &&
            upper_match.producing_local == Rational("3/19") &&
            upper_match.receiving_local == Rational("-10/57"));
  check("both singular handoff coordinates stay within 1/k",
        diffexp2::exact_path_detail::abs(lower_match.producing_local) <=
                Rational("1/3") &&
            diffexp2::exact_path_detail::abs(lower_match.receiving_local) <=
                Rational("1/3") &&
            diffexp2::exact_path_detail::abs(upper_match.producing_local) <=
                Rational("1/3") &&
            diffexp2::exact_path_detail::abs(upper_match.receiving_local) <=
                Rational("1/3"));

  ExactPathTopology unsafe_topology;
  unsafe_topology.singular_points = {Rational(1)};
  ExactArmRequest beyond_half{
      Rational(0), Rational(1),
      {ExactAffineChart{"half-left", Rational(0), Rational(2),
                        Rational("1/2"), false},
       ExactAffineChart{"half-right-singular", Rational(1), Rational(1),
                        Rational(1), true}},
      unsafe_topology};
  bool half_rejected = false;
  try {
    (void)diffexp2::plan_exact_arm(beyond_half);
  } catch (const ExactPathPlanningError&) {
    half_rejected = true;
  }
  check("singular balanced approach beyond the conditioned overlap is rejected",
        half_rejected);

  auto rewrite_match = [](diffexp2::ExactArmPlan& plan,
                          const Rational& physical) {
    auto& match = plan.matches.back();
    const auto& left = plan.charts[plan.charts.size() - 2];
    const auto& right = plan.charts.back();
    match.physical = physical;
    match.producing_local =
        diffexp2::exact_path_detail::local_coordinate(left, physical);
    match.receiving_local =
        diffexp2::exact_path_detail::local_coordinate(right, physical);
    plan.tiles[plan.tiles.size() - 2].physical_end = physical;
    plan.tiles[plan.tiles.size() - 2].local_end = match.producing_local;
    plan.tiles.back().physical_begin = physical;
    plan.tiles.back().local_begin = match.receiving_local;
  };
  auto rejects_plan = [](const diffexp2::ExactArmPlan& plan) {
    try {
      diffexp2::validate_exact_arm_plan(plan);
      return false;
    } catch (const ExactPathPlanningError&) {
      return true;
    }
  };
  auto wrong_side = arms.upper;
  rewrite_match(wrong_side, Rational("101/100"));
  check("singular balanced handoff on the far side is rejected",
        rejects_plan(wrong_side));
  auto noncanonical = arms.upper;
  rewrite_match(noncanonical, Rational("55/91"));
  check("an otherwise safe singular handoff cannot forge the canonical balance",
        rejects_plan(noncanonical));
}

void singular_true_radius_does_not_override_affine_conditioning() {
  ExactPathPlanOptions options;
  options.division_order = 4;
  ExactPathTopology endpoint_topology;
  endpoint_topology.singular_points = {Rational(1)};
  const ExactAffineChart endpoint{"endpoint", Rational(1), Rational("1/101"),
                                  Rational(100), true};
  ExactArmRequest unsafe{
      Rational("100/101"), Rational(1),
      {ExactAffineChart{"pre-endpoint", Rational("100/101"),
                        Rational("1/101"), Rational(100), false},
       endpoint},
      endpoint_topology};
  bool rejected = false;
  try {
    (void)diffexp2::plan_exact_arm(unsafe, options);
  } catch (const ExactPathPlanningError&) {
    rejected = true;
  }
  check("large true radii cannot authorize a near-unit singular local coordinate",
        rejected);

  ExactArmRequest conditioned{
      Rational("999/1000"), Rational(1),
      {ExactAffineChart{"conditioned-pre-endpoint", Rational("999/1000"),
                        Rational("1/101"), Rational(100), false},
       endpoint},
      endpoint_topology};
  const auto plan = diffexp2::plan_exact_arm(conditioned, options);
  const auto& match = plan.matches.front();
  check("a closer regular chart gives a singular handoff within both 1/k frames",
        diffexp2::exact_path_detail::abs(match.producing_local) <=
                Rational("1/4") &&
            diffexp2::exact_path_detail::abs(match.receiving_local) <=
                Rational("1/4"));
}

void forbidden_match_and_independent_arms() {
  ExactPathTopology forbidden_topology;
  // A symbolic producer has declared the receiver-owned -1/k point x=1/3
  // unavailable as a handoff (for example a singular/branch split), so the
  // planner must choose another exact point in the overlap.
  forbidden_topology.boundary_points = {Rational("1/3"), Rational("3/5")};
  forbidden_topology.branch_sheets = {
      ExactBranchSheet{"forbidden-threshold", -1}};
  ExactArmRequest avoided{
      Rational(0), Rational(1),
      {ExactAffineChart{"avoid-left", Rational(0), Rational(2), Rational(1),
                        false},
       ExactAffineChart{"avoid-right", Rational(1), Rational(2), Rational(1),
                        false}},
      forbidden_topology};
  const auto avoided_plan = diffexp2::plan_exact_arm(avoided);
  check("forbidden singular/boundary candidate is never used as a match",
        avoided_plan.matches.size() == 1 &&
            avoided_plan.matches.front().kind ==
                ExactMatchKind::ForbiddenPointAvoidance &&
            !(avoided_plan.matches.front().physical == Rational("1/3")) &&
            avoided_plan.matches.front().producing_local.sign() > 0 &&
            avoided_plan.matches.front().receiving_local.sign() < 0);

  const auto common_topology = topology();
  ExactArmRequest lower{Rational(1), Rational(0),
                        {middle_chart(), left_chart()}, common_topology};
  ExactArmRequest upper{Rational(1), Rational(2),
                        {middle_chart(), right_chart()}, common_topology};
  const auto arms = diffexp2::plan_exact_independent_arms(lower, upper);
  check("lower and upper arms are independent plans sharing one stable anchor",
        arms.anchor_chart.identity == "middle-singular" &&
            arms.lower.direction == -1 && arms.upper.direction == 1 &&
            arms.lower.matches.front().physical == Rational("1/3") &&
            arms.upper.matches.front().physical == Rational("5/3") &&
            arms.lower.topology.branch_sheets.front().imaginary_sign == -1 &&
            arms.upper.topology.branch_sheets.front().imaginary_sign == -1);

  ExactPathTopology missed_singularity;
  missed_singularity.singular_points = {Rational("1/2")};
  avoided.topology = missed_singularity;
  bool rejected = false;
  try {
    (void)diffexp2::plan_exact_arm(avoided);
  } catch (const ExactPathPlanningError&) {
    rejected = true;
  }
  check("an undeclared-chart real singularity cannot be hidden by match avoidance",
        rejected);
}

void asymmetric_receiver_clearance() {
  ExactPathPlanOptions options;
  options.division_order = 4;
  ExactArmRequest request{
      Rational(0), Rational("3/46"),
      {ExactAffineChart{"asymmetric-left", Rational(0), Rational("2/23"),
                        Rational(1), false},
       ExactAffineChart{"asymmetric-right", Rational("1/23"),
                        Rational("1/23"), Rational(1), false}},
      ExactPathTopology{}};
  const auto plan = diffexp2::plan_exact_arm(request, options);
  const auto& match = plan.matches.front();
  check("ordinary exact planning honors receiver -1/k with producer half-radius clearance",
        match.kind == ExactMatchKind::BalancedSafeOverlap &&
            match.physical == Rational("3/92") &&
            match.producing_local == Rational("3/8") &&
            match.receiving_local == Rational("-1/4"));
}

}  // namespace

int main() {
  try {
    classic_exact_chain_and_reverse();
    box_bubble_singular_endpoint_handoffs();
    singular_true_radius_does_not_override_affine_conditioning();
    forbidden_match_and_independent_arms();
    asymmetric_receiver_clearance();
  } catch (const std::exception& error) {
    ++failed;
    std::cout << "  FAIL: unexpected exception: " << error.what() << '\n';
  }
  std::cout << "Exact path planner tests: " << passed << " passed, " << failed
            << " failed\n";
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
