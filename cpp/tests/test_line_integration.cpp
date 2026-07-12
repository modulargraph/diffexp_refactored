#include "diffexp2/line_integration.hpp"

#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

using diffexp2::ComplexBall;
using diffexp2::ExactScalarDescriptor;
using diffexp2::LineIntegrationScope;
using diffexp2::LocalSector;
using diffexp2::LocalSolution;
using diffexp2::MonomialIntegrationOptions;
using diffexp2::NativeIntegrationError;
using diffexp2::NativeIntegrationErrorCode;
using diffexp2::Prescription;
using diffexp2::Rational;
using diffexp2::RealEvaluationPoint;
using diffexp2::SectorMonomialTag;
using diffexp2::StoredLineIntegrationOptions;

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

bool overlaps(const ComplexBall& value, const ComplexBall& expected) {
  return (value - expected).contains_zero();
}

template <typename Scalar>
LocalSolution<Scalar> base_solution(std::uint32_t dimension) {
  LocalSolution<Scalar> solution;
  solution.chart.center_exact = "0";
  solution.chart.scale_exact = "1";
  solution.chart.radius = ComplexBall(2);
  solution.epsilon = {0, 0};
  solution.taylor_complete_max = 0;
  solution.dimension = dimension;
  return solution;
}

template <typename Scalar>
LocalSector<Scalar> sector(const std::string& a, const std::string& b,
                           std::uint32_t p,
                           std::vector<Scalar> coefficients) {
  LocalSector<Scalar> result;
  result.a = ExactScalarDescriptor::rational(a);
  result.b = ExactScalarDescriptor::rational(b);
  result.log_power = p;
  result.coefficients = std::move(coefficients);
  return result;
}

void grouped_center_and_component_reuse() {
  auto solution = base_solution<Rational>(2);
  solution.sectors = {
      sector<Rational>("-1", "0", 0, {Rational(1), Rational(2)}),
      sector<Rational>("-1", "0", 0, {Rational(-1), Rational(-2)}),
      sector<Rational>("0", "0", 0, {Rational(3), Rational(5)})};

  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {0, 0};
  const auto result = diffexp2::integrate_stored_local_line(
      solution, RealEvaluationPoint::rational("0"),
      RealEvaluationPoint::rational("1/2"), options);

  check("identical divergent center monomials cancel exactly before the primitive",
        result.diagnostics.cancelled_divergent_groups == 1 &&
            result.diagnostics.zero_groups_skipped == 1 &&
            overlaps(result.value.at(0, 0),
                     ComplexBall::from_strings("3/2")) &&
            overlaps(result.value.at(0, 1),
                     ComplexBall::from_strings("5/2")));
  check("one grouped primitive is reused across both components",
        result.diagnostics.input_monomial_cells == 3 &&
            result.diagnostics.grouped_monomials == 2 &&
            result.diagnostics.primitive_evaluations == 1 &&
            result.diagnostics.primitive_component_applications == 2 &&
            result.diagnostics.primitive_component_reuses == 1);
  check("line result states its stored-truncation-only scope",
        result.scope == LineIntegrationScope::StoredTruncation &&
            result.value.error.guarantee == diffexp2::ErrorGuarantee::None &&
            result.diagnostics.detail.find("no unseen-tail") !=
                std::string::npos);
}

void acb_cancellation_is_exact_only() {
  auto solution = base_solution<ComplexBall>(1);
  ComplexBall uncertain(1);
  arb_add_error_2exp_si(acb_realref(uncertain.raw()), -100);
  solution.sectors = {
      sector<ComplexBall>("-1", "0", 0, {uncertain}),
      sector<ComplexBall>("-1", "0", 0, {ComplexBall(-1)})};
  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {0, 0};

  bool enclosure_rejected = false;
  try {
    (void)diffexp2::integrate_stored_local_line(
        solution, RealEvaluationPoint::rational("0"),
        RealEvaluationPoint::rational("1/2"), options);
  } catch (const NativeIntegrationError& error) {
    enclosure_rejected =
        error.code == NativeIntegrationErrorCode::UncertifiedCancellation;
  }
  solution.sectors[0].coefficients[0] = ComplexBall(1);
  const auto exact = diffexp2::integrate_stored_local_line(
      solution, RealEvaluationPoint::rational("0"),
      RealEvaluationPoint::rational("1/2"), options);
  check("Acb divergent cancellation requires the exact singleton zero",
        enclosure_rejected && exact.value.at(0, 0).is_zero() &&
            exact.diagnostics.cancelled_divergent_groups == 1);
}

void stored_branch_selection() {
  auto solution = base_solution<Rational>(1);
  solution.sectors = {
      sector<Rational>("-1/2", "0", 0, {Rational(1)})};
  solution.prescriptions.push_back(Prescription{"x", -1, 1, 1});
  const auto lower = RealEvaluationPoint::rational("-1/2");
  const auto upper = RealEvaluationPoint::rational("-1/4");
  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {0, 0};
  const auto result = diffexp2::integrate_stored_local_line(
      solution, lower, upper, options);

  MonomialIntegrationOptions expected_options;
  expected_options.complete_max = 0;
  expected_options.imaginary_sign = -1;
  const auto expected = diffexp2::integrate_sector_monomial(
      SectorMonomialTag::rational("-1/2", "0", 0), lower, upper,
      expected_options);
  bool conflict_rejected = false;
  options.imaginary_sign = 1;
  try {
    (void)diffexp2::integrate_stored_local_line(solution, lower, upper,
                                                options);
  } catch (const NativeIntegrationError& error) {
    conflict_rejected =
        error.code == NativeIntegrationErrorCode::MissingBranchPrescription;
  }
  check("stored -i0 branch is selected and a conflicting override is rejected",
        result.imaginary_sign == std::optional<std::int32_t>(-1) &&
            overlaps(result.value.at(0, 0), expected.coefficient(0)) &&
            conflict_rejected);
}

void upper_halo_and_first_unseen_fail_loudly() {
  auto halo = base_solution<Rational>(1);
  halo.sectors = {
      sector<Rational>("-1", "2", 1, {Rational(1)})};
  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {0, 0};
  bool halo_rejected = false;
  try {
    (void)diffexp2::integrate_stored_local_line(
        halo, RealEvaluationPoint::rational("0"),
        RealEvaluationPoint::rational("1/2"), options);
  } catch (const NativeIntegrationError& error) {
    halo_rejected =
        error.code == NativeIntegrationErrorCode::IncompleteEpsilonWindow &&
        std::string(error.what()).find("eps^-1") != std::string::npos;
  }

  auto unseen = base_solution<Rational>(1);
  unseen.sectors = {
      sector<Rational>("-3", "0", 0, {Rational(0)})};
  bool unseen_rejected = false;
  try {
    (void)diffexp2::integrate_stored_local_line(
        unseen, RealEvaluationPoint::rational("0"),
        RealEvaluationPoint::rational("1/2"), options);
  } catch (const NativeIntegrationError& error) {
    unseen_rejected =
        error.code == NativeIntegrationErrorCode::IncompleteTaylorWindow;
  }
  check("1/eps primitive requires an honest coefficient upper halo",
        halo_rejected);
  check("nonintegrable first unseen b=0 Taylor monomial is rejected",
        unseen_rejected);
}

}  // namespace

int main() {
  ComplexBall::set_precision(512);
  grouped_center_and_component_reuse();
  acb_cancellation_is_exact_only();
  stored_branch_selection();
  upper_halo_and_first_unseen_fail_loudly();
  std::cout << "Results: " << passed << " / " << (passed + failed)
            << " tests passed\n";
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
