#include "diffexp2/tail_majorant.hpp"

#include <cstdlib>
#include <iostream>
#include <optional>
#include <string>
#include <utility>

using diffexp2::ComplexBall;
using diffexp2::EpsilonWindow;
using diffexp2::ErrorGuarantee;
using diffexp2::ExactScalarDescriptor;
using diffexp2::JordanBlock;
using diffexp2::LineIntegrationScope;
using diffexp2::LocalSector;
using diffexp2::LocalSolution;
using diffexp2::Magnitude;
using diffexp2::MatrixEntry;
using diffexp2::MatrixShift;
using diffexp2::PreparedMatrix;
using diffexp2::PreparedRecurrenceOperator;
using diffexp2::Rational;
using diffexp2::RealEvaluationPoint;
using diffexp2::RecurrenceProblem;
using diffexp2::ScalarTraits;
using diffexp2::StepCase;
using diffexp2::StoredLineIntegrationOptions;
using diffexp2::TailMajorantStatus;

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

template <typename Scalar>
struct Fixture {
  PreparedRecurrenceOperator<Scalar> prepared;
  RecurrenceProblem<Scalar> problem;
  LocalSolution<Scalar> solution;
};

template <typename Scalar>
Fixture<Scalar> exponential_fixture(
    std::optional<long> q_linear = std::nullopt) {
  Fixture<Scalar> fixture;
  auto& prepared = fixture.prepared;
  prepared.dimension = 1;
  prepared.frame_base = 0;
  prepared.frame_width = 1;
  prepared.d_lags = {{{0, ScalarTraits<Scalar>::one()}}};
  if (q_linear.has_value())
    prepared.d_lags.push_back(
        {{0, ScalarTraits<Scalar>::integer(*q_linear)}});
  prepared.nhat_lags.resize(2);
  prepared.nhat_lags[0].valuations = {diffexp2::kCompleteInfinity};
  MatrixShift<Scalar> linear;
  linear.shift = 0;
  linear.entries.push_back(
      MatrixEntry<Scalar>{0, 0, ScalarTraits<Scalar>::one()});
  prepared.nhat_lags[1].polynomial.push_back(std::move(linear));
  prepared.nhat_lags[1].valuations = {0};
  prepared.d0_inverse_scalar = ScalarTraits<Scalar>::one();
  prepared.chop_digits = 50;
  prepared.blocks = {JordanBlock{{0}}};
  PreparedMatrix<Scalar> assembly;
  assembly.identity = true;
  assembly.valuations = {0};
  prepared.assembly_matrix = std::move(assembly);

  auto& problem = fixture.problem;
  problem.dimension = 1;
  problem.nmax = 4;
  problem.log_max = 0;
  problem.frame_base = 0;
  problem.frame_width = 1;
  problem.has_initial = true;
  problem.a_target = ScalarTraits<Scalar>::zero();
  problem.b_target = ScalarTraits<Scalar>::zero();
  problem.a_shift_min = 0;
  for (long n = 0; n <= 4; ++n)
    problem.a_shifts.push_back(ScalarTraits<Scalar>::integer(n));
  problem.schedule.resize(5);
  for (long n = 0; n <= 4; ++n)
    problem.schedule[static_cast<std::size_t>(n)] = {
        {StepCase::Taylor, ScalarTraits<Scalar>::integer(n),
         ScalarTraits<Scalar>::zero()}};
  problem.initial = {ScalarTraits<Scalar>::one()};
  problem.initial_validity = {diffexp2::kCompleteInfinity};
  problem.return_u = true;

  auto recurrence = diffexp2::RecurrenceSolver<Scalar>(
      problem, prepared).run();
  auto assembled = diffexp2::assemble_recurrence(
      prepared, problem, recurrence);
  auto& solution = fixture.solution;
  solution.chart.center_exact = "0";
  solution.chart.scale_exact = "1";
  solution.chart.radius = ComplexBall::from_strings("2");
  solution.epsilon = {assembled.min_power, assembled.complete_max};
  solution.taylor_complete_max = problem.nmax;
  solution.dimension = 1;
  solution.checkpoint_identity = "local:exp:eo4";
  LocalSector<Scalar> sector;
  sector.a = ExactScalarDescriptor::rational("0");
  sector.b = ExactScalarDescriptor::rational("0");
  sector.log_power = 0;
  sector.coefficients = std::move(assembled.coefficients);
  solution.sectors.push_back(std::move(sector));
  return fixture;
}

ComplexBall exponential(const std::string& value) {
  const auto input = ComplexBall::from_strings(value);
  ComplexBall output;
  acb_exp(output.raw(), input.raw(), ComplexBall::precision());
  return output;
}

void test_rational_vertical_slice() {
  auto fixture = exponential_fixture<Rational>();
  auto built = diffexp2::prepare_regular_homogeneous_tail_model(
      fixture.prepared, fixture.problem, fixture.solution,
      "operator:qtheta-minus-t");
  check("rational regular model is structurally certified",
        built.status == TailMajorantStatus::Certified &&
        built.model.has_value());
  if (!built.model.has_value()) return;

  const auto point = RealEvaluationPoint::rational("1/2");
  auto evaluated = diffexp2::evaluate_local_solution_with_certified_tail(
      fixture.solution, *built.model, point, "1");
  check("point value and theta carry certified unseen-tail envelopes",
        evaluated.tail.status == TailMajorantStatus::Certified &&
        evaluated.evaluation.value.error.guarantee ==
            ErrorGuarantee::Certified &&
        evaluated.evaluation.theta_value.error.guarantee ==
            ErrorGuarantee::Certified);

  const auto exact = exponential("1/2");
  const auto value_remainder =
      exact - evaluated.evaluation.value.at(0, 0);
  const auto theta_remainder =
      exact - evaluated.evaluation.theta_value.at(0, 0);
  check("certified point envelopes dominate the true e^t remainders",
        Magnitude::upper_abs(value_remainder) <=
            evaluated.evaluation.value.error.absolute.at(0) &&
        Magnitude::upper_abs(theta_remainder) <=
            evaluated.evaluation.theta_value.error.absolute.at(0));

  StoredLineIntegrationOptions options;
  options.delivered_epsilon = EpsilonWindow{0, 0};
  auto integrated =
      diffexp2::integrate_regular_local_line_with_certified_tail(
          fixture.solution, *built.model,
          RealEvaluationPoint::rational("-1/4"),
          RealEvaluationPoint::rational("1/2"), options, "1");
  check("line result is promoted only with a certified full-local tail",
        integrated.tail.status == TailMajorantStatus::Certified &&
        integrated.integral.scope ==
            LineIntegrationScope::FullLocalWithCertifiedTail &&
        integrated.integral.value.error.guarantee ==
            ErrorGuarantee::Certified);
  const auto exact_integral = exponential("1/2") - exponential("-1/4");
  const auto line_remainder =
      exact_integral - integrated.integral.value.at(0, 0);
  check("integrated tail envelope dominates the true line remainder",
        Magnitude::upper_abs(line_remainder) <=
            integrated.integral.value.error.absolute.at(0));
}

void test_acb_and_loud_noncertification() {
  auto acb = exponential_fixture<ComplexBall>();
  auto acb_model = diffexp2::prepare_regular_homogeneous_tail_model(
      acb.prepared, acb.problem, acb.solution,
      "operator:acb:qtheta-minus-t");
  bool acb_certified = acb_model.status == TailMajorantStatus::Certified &&
                       acb_model.model.has_value();
  if (!acb_certified)
    std::cout << "    Acb model detail: " << acb_model.detail << '\n';
  if (acb_certified) {
    const auto evaluated =
        diffexp2::evaluate_local_solution_with_certified_tail(
            acb.solution, *acb_model.model,
            RealEvaluationPoint::rational("1/3"), "1");
    acb_certified =
        evaluated.tail.status == TailMajorantStatus::Certified &&
        evaluated.evaluation.value.error.guarantee ==
            ErrorGuarantee::Certified;
  }
  check("Acb recurrence coefficients use the same rigorous majorant path",
        acb_certified);

  auto zero_separation = exponential_fixture<Rational>(-2);
  auto zero_model = diffexp2::prepare_regular_homogeneous_tail_model(
      zero_separation.prepared, zero_separation.problem,
      zero_separation.solution, "operator:q-equals-1-minus-2t");
  const auto disk = zero_model.model.has_value()
      ? diffexp2::certify_regular_taylor_disk(*zero_model.model, "1")
      : diffexp2::RegularTaylorDiskCertificate{};
  check("failure to separate q from zero is explicitly inconclusive",
        zero_model.status == TailMajorantStatus::Certified &&
        disk.status == TailMajorantStatus::Inconclusive);

  auto regulated = exponential_fixture<Rational>();
  regulated.problem.b_target = Rational(1);
  regulated.solution.sectors.front().b =
      ExactScalarDescriptor::rational("1");
  auto regulated_model =
      diffexp2::prepare_regular_homogeneous_tail_model(
          regulated.prepared, regulated.problem, regulated.solution,
          "operator:regulated");
  check("regulated power is explicit unsupported, never a Taylor proof",
        regulated_model.status == TailMajorantStatus::Unsupported &&
        !regulated_model.model.has_value());

  auto sourced = exponential_fixture<Rational>();
  sourced.problem.source = diffexp2::SourceData<Rational>{};
  auto sourced_model = diffexp2::prepare_regular_homogeneous_tail_model(
      sourced.prepared, sourced.problem, sourced.solution,
      "operator:sourced");
  check("unseen source tail is explicit unsupported",
        sourced_model.status == TailMajorantStatus::Unsupported &&
        !sourced_model.model.has_value());

  auto corrupted = exponential_fixture<Rational>();
  corrupted.solution.sectors.front().coefficients.back() += Rational(1);
  auto corrupted_model =
      diffexp2::prepare_regular_homogeneous_tail_model(
          corrupted.prepared, corrupted.problem, corrupted.solution,
          "operator:corrupted-local");
  check("unrelated retained coefficients cannot inherit a tail proof",
        corrupted_model.status == TailMajorantStatus::Inconclusive &&
        !corrupted_model.model.has_value());
}

}  // namespace

int main() {
  ComplexBall::set_precision(256);
  test_rational_vertical_slice();
  test_acb_and_loud_noncertification();
  std::cout << "Results: " << passed << " / " << (passed + failed)
            << " tests passed\n";
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
