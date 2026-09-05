#include "diffexp/kernel/local_solution.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

using diffexp::kernel::ComplexBall;
using diffexp::kernel::EpsilonMatrix;
using diffexp::kernel::EpsilonWindow;
using diffexp::kernel::EvaluationOptions;
using diffexp::kernel::ExactScalarDescriptor;
using diffexp::kernel::LocalSector;
using diffexp::kernel::LocalSolution;
using diffexp::kernel::Magnitude;
using diffexp::kernel::Prescription;
using diffexp::kernel::RealEvaluationPoint;
using diffexp::kernel::ResidualVerdict;

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

LocalSolution<ComplexBall> regular_linear_solution() {
  LocalSolution<ComplexBall> solution;
  solution.chart.center_exact = "0";
  solution.chart.scale_exact = "1";
  solution.chart.radius = ComplexBall(2);
  solution.epsilon = EpsilonWindow{0, 0};
  solution.taylor_complete_max = 1;
  solution.dimension = 1;
  LocalSector<ComplexBall> sector;
  sector.a = ExactScalarDescriptor::rational("0");
  sector.b = ExactScalarDescriptor::rational("0");
  sector.log_power = 0;
  sector.coefficients = {ComplexBall(1), ComplexBall(2)};
  solution.sectors.push_back(std::move(sector));
  return solution;
}

void test_direct_value_theta_and_residual() {
  const auto evaluation = diffexp::kernel::evaluate_local_solution(
      regular_linear_solution(), RealEvaluationPoint::rational("1/4"));
  const auto expected_value = ComplexBall::from_strings("3/2");
  const auto expected_theta = ComplexBall::from_strings("1/2");
  check("direct local value enclosure",
        (evaluation.value.at(0, 0) - expected_value).contains_zero());
  check("direct theta-value enclosure",
        (evaluation.theta_value.at(0, 0) - expected_theta).contains_zero());

  EpsilonMatrix op;
  op.epsilon = {0, 0};
  op.dimension = 1;
  op.coefficients = {ComplexBall::from_strings("1/3")};
  const auto certificate = diffexp::kernel::certify_theta_residual(
      evaluation, op, std::nullopt, Magnitude::decimal("1e-60"));
  check("Acb stored-polynomial residual certificate",
        certificate.verdict == ResidualVerdict::Pass);
}

LocalSolution<ComplexBall> half_power_solution() {
  LocalSolution<ComplexBall> solution;
  solution.chart.center_exact = "0";
  solution.chart.radius = ComplexBall(2);
  solution.epsilon = {0, 0};
  solution.taylor_complete_max = 0;
  solution.dimension = 1;
  LocalSector<ComplexBall> sector;
  sector.a = ExactScalarDescriptor::rational("1/2");
  sector.b = ExactScalarDescriptor::rational("0");
  sector.coefficients = {ComplexBall(1)};
  solution.sectors.push_back(std::move(sector));
  solution.prescriptions.push_back(Prescription{"x", 1, 1, 1});
  return solution;
}

void test_two_branch_rims() {
  const auto solution = half_power_solution();
  const auto upper = diffexp::kernel::evaluate_local_solution(
      solution, RealEvaluationPoint::rational("-1"));
  EvaluationOptions lower_options;
  lower_options.imaginary_sign = -1;
  const auto lower = diffexp::kernel::evaluate_local_solution(
      solution, RealEvaluationPoint::rational("-1"), lower_options);
  check("fractional power upper rim",
        arb_is_positive(acb_imagref(upper.value.at(0, 0).raw())));
  check("fractional power lower rim",
        arb_is_negative(acb_imagref(lower.value.at(0, 0).raw())));
  check("opposite rims are conjugate",
        (upper.value.at(0, 0) + lower.value.at(0, 0)).contains_zero());
}

void test_epsilon_exponent_and_log_shift() {
  LocalSolution<ComplexBall> solution;
  solution.chart.center_exact = "0";
  solution.chart.radius = ComplexBall(2);
  solution.epsilon = {0, 2};
  solution.taylor_complete_max = 0;
  solution.dimension = 1;
  LocalSector<ComplexBall> sector;
  sector.a = ExactScalarDescriptor::rational("0");
  sector.b = ExactScalarDescriptor::rational("1");
  sector.coefficients = {ComplexBall(1), ComplexBall(0), ComplexBall(0)};
  solution.sectors.push_back(std::move(sector));
  const auto evaluation = diffexp::kernel::evaluate_local_solution(
      solution, RealEvaluationPoint::rational("1/2"));
  ComplexBall log_half;
  acb_log(log_half.raw(), ComplexBall::from_strings("1/2").raw(),
          ComplexBall::precision());
  check("epsilon exponent first coefficient",
        (evaluation.value.at(1, 0) - log_half).contains_zero());

  solution.epsilon = {-1, 1};
  solution.sectors[0].b = ExactScalarDescriptor::rational("0");
  solution.sectors[0].log_power = 1;
  solution.sectors[0].coefficients = {ComplexBall(2), ComplexBall(0),
                                      ComplexBall(0)};
  const auto log_evaluation = diffexp::kernel::evaluate_local_solution(
      solution, RealEvaluationPoint::rational("1/2"));
  check("explicit log sector shifts epsilon minimum",
        log_evaluation.value.epsilon.min_power == 0 &&
        (log_evaluation.value.at(0, 0) - ComplexBall(2) * log_half)
            .contains_zero());
}

void test_pole_rows_feed_finite_orders() {
  LocalSolution<ComplexBall> solution;
  solution.chart.center_exact = "0";
  solution.chart.scale_exact = "1";
  solution.chart.radius = ComplexBall(2);
  solution.epsilon = {-1, 1};
  solution.taylor_complete_max = 0;
  solution.dimension = 1;

  LocalSector<ComplexBall> fixed;
  fixed.a = ExactScalarDescriptor::rational("0");
  fixed.b = ExactScalarDescriptor::rational("0");
  fixed.log_power = 0;
  fixed.coefficients = {
      ComplexBall(-1), ComplexBall(0), ComplexBall(0)};
  solution.sectors.push_back(std::move(fixed));

  LocalSector<ComplexBall> moving;
  moving.a = ExactScalarDescriptor::rational("0");
  moving.b = ExactScalarDescriptor::rational("1");
  moving.log_power = 0;
  moving.coefficients = {
      ComplexBall(1), ComplexBall(0), ComplexBall(0)};
  solution.sectors.push_back(std::move(moving));

  const auto evaluation = diffexp::kernel::evaluate_local_solution(
      solution, RealEvaluationPoint::rational("1/2"));
  ComplexBall log_half;
  acb_log(log_half.raw(), ComplexBall::from_strings("1/2").raw(),
          ComplexBall::precision());
  check("cancelled pole rows still determine the finite coefficient",
        evaluation.value.at(-1, 0).contains_zero() &&
        (evaluation.value.at(0, 0) - log_half).contains_zero());
}

void test_loud_branch_failures() {
  auto missing = half_power_solution();
  missing.prescriptions.clear();
  bool missing_loud = false;
  try {
    (void)diffexp::kernel::evaluate_local_solution(
        missing, RealEvaluationPoint::rational("-1"));
  } catch (const std::domain_error&) {
    missing_loud = true;
  }
  check("missing branch prescription is loud", missing_loud);

  auto conflict = half_power_solution();
  conflict.prescriptions.push_back(Prescription{"1-x", 1, 1, -1});
  bool conflict_loud = false;
  try {
    (void)diffexp::kernel::evaluate_local_solution(
        conflict, RealEvaluationPoint::rational("-1"));
  } catch (const std::domain_error&) {
    conflict_loud = true;
  }
  check("conflicting branch prescriptions are loud", conflict_loud);
}

}  // namespace

int main() {
  ComplexBall::set_precision(256);
  test_direct_value_theta_and_residual();
  test_two_branch_rims();
  test_epsilon_exponent_and_log_shift();
  test_pole_rows_feed_finite_orders();
  test_loud_branch_failures();
  std::cout << "Results: " << passed << " / " << (passed + failed)
            << " tests passed\n";
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
