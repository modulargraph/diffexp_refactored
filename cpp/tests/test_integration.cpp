#include "diffexp2/integration.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using diffexp2::ComplexBall;
using diffexp2::EndpointLimitOptions;
using diffexp2::EpsilonFrame;
using diffexp2::ExactScalarDescriptor;
using diffexp2::LocalSector;
using diffexp2::LocalSolution;
using diffexp2::MonomialIntegrationOptions;
using diffexp2::NativeIntegrationError;
using diffexp2::NativeIntegrationErrorCode;
using diffexp2::Rational;
using diffexp2::RealEvaluationPoint;
using diffexp2::SectorMonomialTag;

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

ComplexBall log_ball(const std::string& value) {
  ComplexBall out;
  const auto input = ComplexBall::from_strings(value);
  acb_log(out.raw(), input.raw(), ComplexBall::precision());
  return out;
}

ComplexBall imaginary_pi() {
  ComplexBall out;
  acb_const_pi(out.raw(), ComplexBall::precision());
  acb_mul_onei(out.raw(), out.raw());
  return out;
}

void frame_smoke() {
  const EpsilonFrame<Rational> a(0, {Rational(1), Rational(2), Rational(3)});
  const EpsilonFrame<Rational> b(-1, {Rational(4), Rational(5), Rational(6)});
  const auto sum = a + b;
  const auto product = a * b;
  check("honest frame add intersects complete maxima",
        sum.min_power() == -1 && sum.complete_max() == 1 &&
            sum.coefficient(-1) == Rational(4) &&
            sum.coefficient(0) == Rational(6) &&
            sum.coefficient(1) == Rational(8));
  check("honest frame product uses two-sided Cauchy window",
        product.min_power() == -1 && product.complete_max() == 1 &&
            product.coefficient(-1) == Rational(4) &&
            product.coefficient(0) == Rational(13));
}

void endpoint_limit_smoke() {
  LocalSolution<Rational> solution;
  solution.chart.center_exact = "0";
  solution.chart.scale_exact = "1";
  solution.chart.radius = ComplexBall(2);
  solution.epsilon = {0, 1};
  solution.taylor_complete_max = 1;
  solution.dimension = 1;

  LocalSector<Rational> first;
  first.a = ExactScalarDescriptor::rational("-1");
  first.b = ExactScalarDescriptor::rational("0");
  first.log_power = 0;
  first.coefficients = {Rational(3), Rational(7), Rational(5), Rational(11)};
  LocalSector<Rational> cancelling = first;
  cancelling.coefficients =
      {Rational(-3), Rational(0), Rational(-5), Rational(0)};
  LocalSector<Rational> regulated = first;
  regulated.a = ExactScalarDescriptor::rational("-3");
  regulated.b = ExactScalarDescriptor::rational("2");
  regulated.coefficients =
      {Rational(13), Rational(17), Rational(19), Rational(23)};
  solution.sectors = {first, cancelling, regulated};

  EndpointLimitOptions options;
  options.imaginary_sign = 1;
  const auto result = diffexp2::endpoint_sector_limit(solution, options);
  check("endpoint gate merges absolute monomials and reads buried t^0",
        result.values.size() == 1 &&
            result.values.front().min_power() == 0 &&
            result.values.front().complete_max() == 1 &&
            overlaps(result.values.front().coefficient(0), ComplexBall(7)) &&
            overlaps(result.values.front().coefficient(1), ComplexBall(11)) &&
            result.dropped_regulated_sectors == 1);
}

void primitive_smoke() {
  const auto zero = RealEvaluationPoint::rational("0");
  MonomialIntegrationOptions options;
  options.complete_max = 2;

  const auto regularized = diffexp2::integrate_sector_monomial(
      SectorMonomialTag::rational("-2", "2", 0), zero,
      RealEvaluationPoint::rational("1/2"), options);
  check("b!=0 endpoint analytically regularizes a negative power",
        regularized.min_power() == 0 &&
            overlaps(regularized.coefficient(0), ComplexBall(-2)));

  const auto pole_log = diffexp2::integrate_sector_monomial(
      SectorMonomialTag::rational("-1", "2", 1), zero,
      RealEvaluationPoint::rational("1/2"), options);
  check("endpoint pole cell retains its genuine Laurent row at log depth",
        pole_log.min_power() == -1 &&
            overlaps(pole_log.coefficient(-1),
                     ComplexBall::from_strings("-1/4")));

  options.imaginary_sign = 1;
  const auto paired = diffexp2::integrate_sector_monomial(
      SectorMonomialTag::rational("-1", "2", 0),
      RealEvaluationPoint::rational("-1/5"),
      RealEvaluationPoint::rational("1/4"), options);
  const auto expected_paired = log_ball("5/4") - imaginary_pi();
  check("m=-1 crossing is paired before epsilon-window arithmetic",
        paired.min_power() == 0 &&
            overlaps(paired.coefficient(0), expected_paired));

  const auto pv = diffexp2::integrate_sector_monomial(
      SectorMonomialTag::rational("-1", "0", 0),
      RealEvaluationPoint::rational("-1/5"),
      RealEvaluationPoint::rational("1/4"), options);
  check("integer b=0 crossing uses the real principal value",
        pv.min_power() == 0 &&
            arb_is_zero(acb_imagref(pv.coefficient(0).raw())) &&
            overlaps(pv.coefficient(0), log_ball("5/4")));

  auto symbolic = SectorMonomialTag::rational("-1", "2", 0);
  symbolic.b = ExactScalarDescriptor::symbolic(
      "rho", {"rho"}, diffexp2::TruthValue::No,
      diffexp2::TruthValue::Unknown, diffexp2::ExactSign::Positive,
      ComplexBall(2));
  bool rejected = false;
  try {
    (void)diffexp2::integrate_sector_monomial(
        symbolic, zero, RealEvaluationPoint::rational("1/2"), options);
  } catch (const NativeIntegrationError& error) {
    rejected =
        error.code == NativeIntegrationErrorCode::UnsupportedSymbolicRegulator;
  }
  check("symbolic regulator is rejected instead of numerically sampled",
        rejected);
}

}  // namespace

int main() {
  ComplexBall::set_precision(512);
  frame_smoke();
  endpoint_limit_smoke();
  primitive_smoke();
  std::cout << "Results: " << passed << " / " << (passed + failed)
            << " tests passed\n";
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
