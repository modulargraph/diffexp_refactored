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

void test_finite_epsilon_physical_tail_model() {
  diffexp2::PreparedPhysicalClearedODE<ComplexBall> equation;
  equation.dimension = 1;
  diffexp2::ExactEpsilonRational<ComplexBall> q0;
  q0.zero = false;
  q0.valuation = 0;
  q0.numerator = {ComplexBall(1)};
  q0.denominator = {ComplexBall(1)};
  equation.q_lags = {q0};
  equation.c_lags.resize(2);
  diffexp2::ExactEpsilonRational<ComplexBall> one_plus_epsilon;
  one_plus_epsilon.zero = false;
  one_plus_epsilon.valuation = 0;
  one_plus_epsilon.numerator = {ComplexBall(1), ComplexBall(1)};
  one_plus_epsilon.denominator = {ComplexBall(1)};
  equation.c_lags[1].push_back(
      {0, 0, std::move(one_plus_epsilon)});
  equation.owner_signature_identity = "physical-owner:coupled-exp";
  equation.payload_identity = "physical-payload:coupled-exp";
  equation.exact_payload_record = "physical-record:coupled-exp";

  diffexp2::EpsilonVector initial;
  initial.epsilon = {0, 2};
  initial.dimension = 1;
  initial.coefficients = {
      ComplexBall(1), ComplexBall(0),
      ComplexBall::from_strings("1e100")};
  const auto evolution = diffexp2::evolve_ordinary_center_value(
      equation, initial, 4);
  LocalSolution<ComplexBall> solution;
  solution.chart.center_exact = "0";
  solution.chart.scale_exact = "1";
  solution.chart.radius = ComplexBall::from_strings("2");
  solution.epsilon = initial.epsilon;
  solution.taylor_complete_max = 4;
  solution.dimension = 1;
  solution.checkpoint_identity = "local:physical-coupled-exp:eo4";
  LocalSector<ComplexBall> sector;
  sector.a = ExactScalarDescriptor::rational("0");
  sector.b = ExactScalarDescriptor::rational("0");
  sector.log_power = 0;
  sector.coefficients.assign(solution.sector_size(), ComplexBall(0));
  for (std::uint32_t taylor = 0; taylor <= 4; ++taylor)
    for (std::int32_t power = 0; power <= 2; ++power)
      sector.coefficients[diffexp2::local_detail::sector_index(
          solution, static_cast<std::size_t>(power), taylor, 0)] =
          evolution.at(taylor).at(power, 0);
  solution.sectors.push_back(std::move(sector));

  auto built =
      diffexp2::prepare_physical_regular_homogeneous_tail_model(
          equation, solution);
  check("epsilon-coupled physical ODE builds an augmented tail model",
        built.status == TailMajorantStatus::Certified &&
        built.model.has_value());
  if (!built.model.has_value()) {
    std::cout << "    physical model detail: " << built.detail << '\n';
    return;
  }
  const auto evaluated =
      diffexp2::evaluate_physical_local_solution_with_certified_tail(
          *built.model, RealEvaluationPoint::rational("1/2"), "1");
  check("epsilon-coupled physical point evaluation carries a certified tail",
        evaluated.tail.status == TailMajorantStatus::Certified &&
        evaluated.evaluation.value.error.guarantee ==
            ErrorGuarantee::Certified);
  const auto midpoint_witness =
      diffexp2::evaluate_physical_local_solution_with_certified_tail(
          *built.model, RealEvaluationPoint::rational("1/2"), "5/4");
  const auto outward_witness =
      diffexp2::evaluate_physical_local_solution_with_certified_tail(
          *built.model, RealEvaluationPoint::rational("1/2"),
          "262141/131072");
  check("outward certified witness can satisfy a contract that the midpoint misses",
        midpoint_witness.tail.status ==
                TailMajorantStatus::Certified &&
            outward_witness.tail.status ==
                TailMajorantStatus::Certified &&
            2.0 * outward_witness.tail.value.absolute.at(0)
                      .approximate_upper() <
                midpoint_witness.tail.value.absolute.at(0)
                    .approximate_upper());
  const auto exact0 = exponential("1/2");
  const auto exact1 = ComplexBall::from_strings("1/2") * exact0;
  const auto remainder0 = exact0 - evaluated.evaluation.value.at(0, 0);
  const auto remainder1 = exact1 - evaluated.evaluation.value.at(1, 0);
  check("augmented tail envelope covers coupled epsilon coefficients",
        Magnitude::upper_abs(remainder0) <=
                evaluated.evaluation.value.error.absolute.at(0) &&
            Magnitude::upper_abs(remainder1) <=
                evaluated.evaluation.value.error.absolute.at(1));
  check("private high epsilon data cannot inflate a causal low-prefix tail",
        evaluated.evaluation.value.error.absolute.at(0)
                    .approximate_upper() < 1.0 &&
            evaluated.evaluation.value.error.absolute.at(2)
                    .approximate_upper() > 1e50);

  const auto extended =
      diffexp2::prepare_physical_regular_homogeneous_tail_model(
          equation, solution, 29);
  check("physical handoff can reconstruct a private Taylor certification prefix",
        extended.status == TailMajorantStatus::Certified &&
            extended.model.has_value() &&
            extended.model->taylor_complete_max == 29 &&
            extended.model->reconstructed.taylor_complete_max == 29 &&
            solution.taylor_complete_max == 4);
  if (extended.model.has_value()) {
    const auto extended_evaluated =
        diffexp2::evaluate_physical_local_solution_with_certified_tail(
            *extended.model, RealEvaluationPoint::rational("1/2"), "1");
    check("private certification order tightens the public tail without mutating retained order",
          extended_evaluated.tail.status ==
                  TailMajorantStatus::Certified &&
              extended_evaluated.evaluation.value.error.absolute.at(0)
                      .approximate_upper() <
                  evaluated.evaluation.value.error.absolute.at(0)
                      .approximate_upper() &&
              extended_evaluated.evaluation.value.error.absolute.at(0)
                      .approximate_upper() <
                  1e-6);
  }

  auto independently_rounded = solution;
  const auto independently_rounded_coefficient =
      diffexp2::local_detail::sector_index(
          independently_rounded, 0, 3, 0);
  const auto& recurrence_ball = evolution.at(3).at(0, 0);
  auto overlapping_non_enclosure = recurrence_ball;
  mag_set_ui_2exp_si(
      arb_radref(acb_realref(overlapping_non_enclosure.raw())), 1, -200);
  arf_t shift;
  arf_init(shift);
  arf_set_ui_2exp_si(shift, 1, -200);
  arf_add(arb_midref(acb_realref(overlapping_non_enclosure.raw())),
          arb_midref(acb_realref(overlapping_non_enclosure.raw())),
          shift, ARF_PREC_EXACT, ARF_RND_NEAR);
  arf_clear(shift);
  check("fixture constructs an overlapping non-enclosing Arb replay",
        acb_overlaps(overlapping_non_enclosure.raw(),
                     recurrence_ball.raw()) &&
            !acb_contains(overlapping_non_enclosure.raw(),
                          recurrence_ball.raw()));
  independently_rounded.sectors.front().coefficients[
      independently_rounded_coefficient] =
      std::move(overlapping_non_enclosure);
  const auto independently_rounded_model =
      diffexp2::prepare_physical_regular_homogeneous_tail_model(
          equation, independently_rounded);
  check("physical replay accepts overlapping independent rounding",
        independently_rounded_model.status ==
                TailMajorantStatus::Certified &&
            independently_rounded_model.model.has_value());

  auto corrupted = solution;
  corrupted.sectors.front().coefficients[
      diffexp2::local_detail::sector_index(
          corrupted, 0, corrupted.taylor_complete_max, 0)] +=
      ComplexBall(1);
  const auto corrupted_model =
      diffexp2::prepare_physical_regular_homogeneous_tail_model(
          equation, corrupted, 29);
  check("extended physical tail proof still rejects an unrelated retained Taylor tensor",
        corrupted_model.status == TailMajorantStatus::Inconclusive &&
            !corrupted_model.model.has_value());
}

void test_finite_causal_q_disk_inverse_certificate() {
  // q(t,epsilon) = 1 + t + t^2 + 100 epsilon t.  The diagonal
  // polynomial has both roots on |t|=1, so it is invertible on |t|<=9/10.
  // The former Neumann test rejected even the one-row prefix because
  // |t|+|t|^2 > 1, while the large nilpotent epsilon term made wider
  // prefixes still less likely to pass despite not changing invertibility.
  const std::vector<std::vector<ComplexBall>> safe = {
      {ComplexBall(1), ComplexBall(0)},
      {ComplexBall(1), ComplexBall(100)},
      {ComplexBall(1), ComplexBall(0)}};
  const auto certified =
      diffexp2::tail_majorant_detail::
          certify_finite_causal_q_disk_inverse(
              safe, Rational("9/10"));
  check("adaptive Q disk proof certifies a polynomial rejected by the coarse Neumann norm",
        certified.status == TailMajorantStatus::Certified &&
            certified.prefix_norm_upper.size() == 2 &&
            !certified.diagonal_lower.is_zero() &&
            certified.accepted_tiles != 0);
  check("finite back-substitution retains rather than ignores nilpotent epsilon coupling",
        certified.status == TailMajorantStatus::Certified &&
            certified.prefix_norm_upper.at(0).approximate_upper() <
                certified.prefix_norm_upper.at(1).approximate_upper());

  // q(t) = 1 + 2t + 2t^2 has roots (-1 +/- i)/2, both strictly inside
  // |t|<=4/5.  Adaptive refinement must fail closed rather than manufacture
  // an inverse bound around a genuine zero.
  const std::vector<std::vector<ComplexBall>> unsafe = {
      {ComplexBall(1)},
      {ComplexBall(2)},
      {ComplexBall(2)}};
  const auto rejected =
      diffexp2::tail_majorant_detail::
          certify_finite_causal_q_disk_inverse(
              unsafe, Rational("4/5"));
  check("adaptive Q disk proof remains inconclusive when the disk contains a true zero",
        rejected.status == TailMajorantStatus::Inconclusive);
  const std::vector<std::vector<ComplexBall>> boundary_zero = {
      {ComplexBall(1)},
      {ComplexBall(1)},
      {ComplexBall(1)}};
  const auto boundary_rejected =
      diffexp2::tail_majorant_detail::
          certify_finite_causal_q_disk_inverse(
              boundary_zero, Rational("1"));
  check("Q root on the closed Cauchy circle remains inconclusive",
        boundary_rejected.status ==
            TailMajorantStatus::Inconclusive);

  // A large scalar clearing polynomial is algebraically irrelevant to the
  // physical ODE.  Form C=t*q*(2+100 epsilon), then require the normalized
  // disk path to cancel q before norms and the geometric epsilon weight to
  // recognize the strictly causal 100 epsilon coupling.
  const auto scale = ComplexBall::from_strings("1e20");
  const std::vector<std::vector<ComplexBall>> cleared_q = {
      {scale, ComplexBall(0)},
      {scale, ComplexBall(0)},
      {scale, ComplexBall(0)}};
  const auto c0 = std::vector<ComplexBall>{
      ComplexBall(0), ComplexBall(0)};
  const auto c = std::vector<ComplexBall>{
      scale * ComplexBall(2), scale * ComplexBall(100)};
  const std::vector<std::vector<ComplexBall>> cleared_c = {
      c0, c, c, c};
  const auto normalized =
      diffexp2::tail_majorant_detail::
          certify_physical_normalized_ode_disk_bounds(
              cleared_q, cleared_c, 1, Rational("1/2"), 3);
  const auto unweighted =
      normalized.status == TailMajorantStatus::Certified
          ? diffexp2::tail_majorant_detail::
                weighted_physical_ode_prefix_norm_upper(
                    normalized, 1)
          : std::vector<Magnitude>{};
  const auto weighted =
      normalized.status == TailMajorantStatus::Certified
          ? diffexp2::tail_majorant_detail::
                weighted_physical_ode_prefix_norm_upper(
                    normalized, 16)
          : std::vector<Magnitude>{};
  check("normalized disk proof cancels a large common clearing polynomial before norms",
        normalized.status == TailMajorantStatus::Certified &&
            !unweighted.empty() &&
            unweighted.back().approximate_upper() < 1e6);
  check("geometric epsilon norm suppresses strictly causal growth without dropping it",
        !weighted.empty() &&
            weighted.back().approximate_upper() <
                unweighted.back().approximate_upper());

  const auto zero_ode =
      diffexp2::tail_majorant_detail::
          certify_physical_normalized_ode_disk_bounds(
              std::vector<std::vector<ComplexBall>>{
                  {ComplexBall(1)}},
              std::vector<std::vector<ComplexBall>>{
                  {ComplexBall(0)}},
              1, Rational("1"), 1);
  const auto zero_norm =
      zero_ode.status == TailMajorantStatus::Certified
          ? diffexp2::tail_majorant_detail::
                weighted_physical_ode_prefix_norm_upper(
                    zero_ode, 1)
          : std::vector<Magnitude>{};
  check("identically zero physical C has an exact zero normalized growth bound",
        zero_ode.status == TailMajorantStatus::Certified &&
            zero_norm.size() == 1 && zero_norm.front().is_zero());
  const auto corrupt_c0 =
      diffexp2::tail_majorant_detail::
          certify_physical_normalized_ode_disk_bounds(
              std::vector<std::vector<ComplexBall>>{
                  {ComplexBall(1)}},
              std::vector<std::vector<ComplexBall>>{
                  {ComplexBall(1)}, {ComplexBall(0)}},
              1, Rational("1"), 1);
  check("normalized physical disk proof rejects nonzero structural C(0)",
        corrupt_c0.status == TailMajorantStatus::Inconclusive);
  bool rejected_inexact_weight = false;
  try {
    (void)diffexp2::tail_majorant_detail::
        weighted_physical_ode_prefix_norm_upper(
            normalized, 3);
  } catch (const std::invalid_argument&) {
    rejected_inexact_weight = true;
  }
  check("weighted physical norm rejects non-power-of-two magnitude scaling",
        rejected_inexact_weight);

  diffexp2::PhysicalRegularTaylorTailModel causal_solution;
  causal_solution.epsilon = {0, 1};
  causal_solution.dimension = 1;
  causal_solution.chart.center_exact = "0";
  causal_solution.chart.scale_exact = "1";
  causal_solution.chart.radius_exact = "2";
  causal_solution.chart.radius =
      ComplexBall::from_strings("2");
  causal_solution.q0_inverse_prefix_norm_upper = {
      Magnitude::one(), Magnitude::one()};
  causal_solution.q_operator_prefix_norm_upper = {{
      Magnitude::one(), Magnitude::one()}};
  causal_solution.q_causal_coefficients = {{
      ComplexBall(1), ComplexBall(0)}};
  causal_solution.c_operator_prefix_norm_upper = {
      {Magnitude::zero(), Magnitude::zero()},
      {Magnitude::zero(), Magnitude::from_ui(100)}};
  causal_solution.c_causal_coefficients = {
      {ComplexBall(0), ComplexBall(0)},
      {ComplexBall(0), ComplexBall(100)}};
  causal_solution.initial_row_upper = {
      Magnitude::one(), Magnitude::zero()};
  const auto causal_disk =
      diffexp2::certify_physical_regular_taylor_disk(
          causal_solution, "1/2");
  // df/dt = 100 epsilon f with f(0)=1 has
  // [epsilon^1]f(t)=100t.  Any circle bound at R=1/2 must therefore
  // unweight back to at least 50.
  check("weighted physical Gronwall bound unweights the delivered epsilon coefficient",
        causal_disk.status == TailMajorantStatus::Certified &&
            causal_disk.cauchy_circle_upper.size() == 2 &&
            causal_disk.cauchy_circle_upper[1]
                    .approximate_upper() >= 50.0);
}

}  // namespace

int main() {
  ComplexBall::set_precision(256);
  test_rational_vertical_slice();
  test_acb_and_loud_noncertification();
  test_finite_epsilon_physical_tail_model();
  test_finite_causal_q_disk_inverse_certificate();
  std::cout << "Results: " << passed << " / " << (passed + failed)
            << " tests passed\n";
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
