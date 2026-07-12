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
                          Rational("3/2"), true};
}

ExactAffineChart right_chart() {
  return ExactAffineChart{"right", Rational(2), Rational(1),
                          Rational("4/3"), false};
}

void classic_exact_chain_and_reverse() {
  ExactPathPlanOptions options;
  check("DivisionOrder defaults to three", options.division_order == 3);

  ExactArmRequest forward{Rational(0), Rational(2),
                          {left_chart(), middle_chart(), right_chart()},
                          topology()};
  const auto plan = diffexp2::plan_exact_arm(forward, options);
  check("unequal scales and radii retain classic coupled +/-1/3 matches",
        plan.matches.size() == 2 &&
            plan.matches[0].kind ==
                ExactMatchKind::SymmetricDivisionPoint &&
            plan.matches[0].physical == Rational("1/3") &&
            plan.matches[0].producing_local == Rational("1/3") &&
            plan.matches[0].receiving_local == Rational("-1/3") &&
            plan.matches[1].physical == Rational("5/3") &&
            plan.matches[1].producing_local == Rational("1/3") &&
            plan.matches[1].receiving_local == Rational("-1/3"));
  check("tiles cover exactly once with no gap and flag the singular crossing",
        plan.tiles.size() == 3 &&
            plan.tiles[0].physical_begin == Rational(0) &&
            plan.tiles[0].physical_end == Rational("1/3") &&
            plan.tiles[1].physical_begin == plan.tiles[0].physical_end &&
            plan.tiles[1].physical_end == plan.tiles[2].physical_begin &&
            plan.tiles[2].physical_end == Rational(2) &&
            plan.tiles[1].crosses_singular_center);
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

  ExactArmRequest reverse{Rational(2), Rational(0),
                          {right_chart(), middle_chart(), left_chart()},
                          topology()};
  const auto reversed = diffexp2::plan_exact_arm(reverse, options);
  check("reversed arm mirrors physical order and local match signs",
        reversed.direction == -1 && reversed.matches.size() == 2 &&
            reversed.matches[0].physical == Rational("5/3") &&
            reversed.matches[0].producing_local == Rational("-1/3") &&
            reversed.matches[0].receiving_local == Rational("1/3") &&
            reversed.matches[1].physical == Rational("1/3") &&
            reversed.matches[1].producing_local == Rational("-1/3") &&
            reversed.matches[1].receiving_local == Rational("1/3"));
}

void forbidden_match_and_independent_arms() {
  ExactPathTopology forbidden_topology;
  // A symbolic producer has declared x=1/2 unavailable as a handoff (for
  // example a singular/branch split).  The naive balanced point is exactly
  // 1/2, so the planner must choose another exact point in the overlap.
  forbidden_topology.boundary_points = {Rational("1/2")};
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
            !(avoided_plan.matches.front().physical == Rational("1/2")) &&
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

}  // namespace

int main() {
  try {
    classic_exact_chain_and_reverse();
    forbidden_match_and_independent_arms();
  } catch (const std::exception& error) {
    ++failed;
    std::cout << "  FAIL: unexpected exception: " << error.what() << '\n';
  }
  std::cout << "Exact path planner tests: " << passed << " passed, " << failed
            << " failed\n";
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
