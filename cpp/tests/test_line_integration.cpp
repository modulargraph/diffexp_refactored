#include "diffexp2/line_integration.hpp"

#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using diffexp2::ComplexBall;
using diffexp2::EpsilonVector;
using diffexp2::ExactEpsilonRational;
using diffexp2::ExactScalarDescriptor;
using diffexp2::LineIntegrationScope;
using diffexp2::LocalSector;
using diffexp2::LocalSolution;
using diffexp2::MonomialIntegrationOptions;
using diffexp2::NativeIntegrationError;
using diffexp2::NativeIntegrationErrorCode;
using diffexp2::PreparedRationalTaylorMultiplier;
using diffexp2::PreparedPhysicalClearedODE;
using diffexp2::PreparedSparseLocalMultiplierMatrix;
using diffexp2::PhysicalODEMatrixEntry;
using diffexp2::Prescription;
using diffexp2::Rational;
using diffexp2::RealEvaluationPoint;
using diffexp2::SectorMonomialTag;
using diffexp2::StoredLineIntegrationOptions;

using FusedWorkItem =
    diffexp2::line_integration_detail::FusedMonomialWorkItem;
static_assert(std::is_trivially_copyable_v<FusedWorkItem>);
static_assert(sizeof(FusedWorkItem) == 16);

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

std::size_t coefficient_index(std::size_t epsilon, std::size_t taylor,
                              std::uint32_t component,
                              std::size_t taylor_width,
                              std::uint32_t dimension) {
  return ((epsilon * taylor_width + taylor) * dimension) + component;
}

template <typename Scalar>
PreparedRationalTaylorMultiplier<Scalar> prepared_multiplier(
    std::int32_t shift, std::uint32_t pole,
    std::vector<std::vector<Scalar>> kernels,
    std::string identity) {
  PreparedRationalTaylorMultiplier<Scalar> result;
  result.epsilon_shift = shift;
  result.center_pole_order = pole;
  result.kernels = std::move(kernels);
  result.exact_identity = std::move(identity);
  return result;
}

template <typename Scalar>
bool same_integral_values(const diffexp2::StoredLineIntegral& left,
                          const diffexp2::StoredLineIntegral& right) {
  (void)sizeof(Scalar);
  if (left.value.epsilon.min_power != right.value.epsilon.min_power ||
      left.value.epsilon.complete_max !=
          right.value.epsilon.complete_max ||
      left.value.dimension != right.value.dimension)
    return false;
  for (std::int64_t power = left.value.epsilon.min_power;
       power <= left.value.epsilon.complete_max; ++power)
    if (!overlaps(left.value.at(static_cast<std::int32_t>(power), 0),
                  right.value.at(static_cast<std::int32_t>(power), 0)))
      return false;
  return true;
}

template <typename Scalar>
std::pair<diffexp2::StoredLineIntegral, diffexp2::StoredLineIntegral>
materialized_and_fused(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& row,
    const LocalSolution<Scalar>& source, std::int32_t projected_complete_cap,
    const RealEvaluationPoint& lower, const RealEvaluationPoint& upper,
    const StoredLineIntegrationOptions& options) {
  const auto projected = diffexp2::apply_prepared_scalar_row_window(
      row, source, projected_complete_cap, "line-test-materialized-row");
  if (!projected.has_value())
    throw std::logic_error("line-test row unexpectedly had no active entry");
  return {diffexp2::integrate_stored_local_line(
              *projected, lower, upper, options),
          diffexp2::integrate_prepared_scalar_row_stored(
              row, source, projected_complete_cap, lower, upper, options)};
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

void acb_bounded_divergent_cancellation_is_explicit_and_relative() {
  auto solution = base_solution<ComplexBall>(1);
  ComplexBall within(1);
  // 2^-66 is larger than 1e-20 times either individual O(1) addend,
  // but smaller than 1e-20 times their rigorous l1 contribution norm.
  // This pins the relative cancellation policy to the accumulated scale,
  // rather than a term-count-dependent maximum norm.
  arb_add_error_2exp_si(acb_realref(within.raw()), -66);
  solution.sectors = {
      sector<ComplexBall>("-1", "0", 0, {within}),
      sector<ComplexBall>("-1", "0", 0, {ComplexBall(-1)})};
  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {0, 0};
  options.divergent_cancellation =
      StoredLineIntegrationOptions::BoundedDivergentCancellation{
          diffexp2::Magnitude::decimal("1e-20"), "1e-20",
          "focused-ft-test"};

  const auto accepted = diffexp2::integrate_stored_local_line(
      solution, RealEvaluationPoint::rational("0"),
      RealEvaluationPoint::rational("1/2"), options);
  check("explicit bounded Acb policy accepts only a relative upper-bound proof",
        accepted.value.at(0, 0).is_zero() &&
            accepted.diagnostics.cancelled_divergent_groups == 1 &&
            accepted.diagnostics.bounded_cancelled_divergent_coefficients ==
                1 &&
            accepted.diagnostics.divergent_cancellation_mode ==
                "bounded-relative-acb" &&
            accepted.diagnostics.divergent_relative_tolerance == "1e-20" &&
            accepted.diagnostics.divergent_cancellation_provenance ==
                "focused-ft-test");

  ComplexBall outside(1);
  arb_add_error_2exp_si(acb_realref(outside.raw()), -10);
  solution.sectors[0].coefficients[0] = outside;
  bool above_bound_rejected = false;
  try {
    (void)diffexp2::integrate_stored_local_line(
        solution, RealEvaluationPoint::rational("0"),
        RealEvaluationPoint::rational("1/2"), options);
  } catch (const NativeIntegrationError& error) {
    above_bound_rejected =
        error.code == NativeIntegrationErrorCode::UncertifiedCancellation &&
        std::string(error.what()).find("relative_tolerance=1e-20") !=
            std::string::npos;
  }
  check("bounded Acb policy rejects a zero-containing coefficient above its bound",
        above_bound_rejected);

  auto rational = base_solution<Rational>(1);
  rational.sectors = {
      sector<Rational>("-1", "0", 0, {Rational(1)}),
      sector<Rational>("-1", "0", 0, {Rational("-999/1000")})};
  bool exact_nonzero_rejected = false;
  try {
    (void)diffexp2::integrate_stored_local_line(
        rational, RealEvaluationPoint::rational("0"),
        RealEvaluationPoint::rational("1/2"), options);
  } catch (const NativeIntegrationError& error) {
    exact_nonzero_rejected =
        error.code == NativeIntegrationErrorCode::DivergentEndpoint;
  }
  check("bounded Acb policy never weakens exact Rational divergence",
        exact_nonzero_rejected);
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

void input_error_envelope_is_not_discarded() {
  auto solution = base_solution<Rational>(1);
  solution.sectors = {
      sector<Rational>("0", "0", 0, {Rational(1)})};
  solution.error.frame = {0, 0};
  solution.error.guarantee = diffexp2::ErrorGuarantee::Advisory;
  solution.error.absolute = {diffexp2::Magnitude::one()};
  solution.error.provenance = "test input tail";
  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {0, 0};
  bool rejected = false;
  try {
    (void)diffexp2::integrate_stored_local_line(
        solution, RealEvaluationPoint::rational("1/4"),
        RealEvaluationPoint::rational("1/2"), options);
  } catch (const NativeIntegrationError& error) {
    rejected = error.code == NativeIntegrationErrorCode::UnsupportedExactTag &&
        std::string(error.what()).find("error envelope") != std::string::npos;
  }
  check("stored line integration never silently drops an input error envelope",
        rejected);
}

void fused_rational_row_matches_materialized_projection() {
  auto source = base_solution<Rational>(2);
  source.epsilon = {-1, 1};
  source.taylor_complete_max = 2;
  source.checkpoint_identity = "fused-rational-parity-source";
  source.prescriptions.push_back(Prescription{"x", 1, 1, 1});

  auto first = sector<Rational>(
      "0", "0", 0,
      std::vector<Rational>(source.sector_size(), Rational(0)));
  auto second = sector<Rational>(
      "1", "1/3", 1,
      std::vector<Rational>(source.sector_size(), Rational(0)));
  for (std::size_t epsilon = 0; epsilon < source.epsilon.width(); ++epsilon) {
    for (std::size_t taylor = 0; taylor < source.taylor_width(); ++taylor) {
      for (std::uint32_t component = 0; component < source.dimension;
           ++component) {
        const auto index = coefficient_index(
            epsilon, taylor, component, source.taylor_width(),
            source.dimension);
        first.coefficients[index] = Rational(std::to_string(
            1 + 17 * epsilon + 5 * taylor + 2 * component));
        second.coefficients[index] = Rational(std::to_string(
            -3 - 11 * epsilon + 7 * taylor - component));
      }
    }
  }
  source.sectors = {std::move(first), std::move(second)};

  PreparedSparseLocalMultiplierMatrix<Rational> row;
  row.rows = 1;
  row.columns = 2;
  row.exact_identity = "fused-rational-parity-row";
  row.entries.push_back(
      {0, 0,
       prepared_multiplier<Rational>(
           0, 0,
           {{Rational(1), Rational("1/2"), Rational(0)},
            {Rational(2), Rational(0), Rational("-1/3")},
            {Rational(0), Rational(1), Rational(0)}},
           "rational-parity-left")});
  row.entries.push_back(
      {0, 1,
       prepared_multiplier<Rational>(
           -1, 1,
           {{Rational(-1), Rational(1), Rational("2/3")},
            {Rational("3/2"), Rational(-2), Rational(0)},
            {Rational(1), Rational(0), Rational(0)}},
           "rational-parity-right")});

  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {-2, 0};
  const auto [materialized, fused] = materialized_and_fused(
      row, source, 0, RealEvaluationPoint::rational("1/4"),
      RealEvaluationPoint::rational("1/2"), options);
  check("fused rational row integration matches a materialized multi-entry, multi-sector projection",
        same_integral_values<Rational>(materialized, fused) &&
            materialized.diagnostics.input_monomial_cells ==
                fused.diagnostics.input_monomial_cells &&
            materialized.diagnostics.grouped_monomials ==
                fused.diagnostics.grouped_monomials &&
            materialized.diagnostics.primitive_evaluations ==
                fused.diagnostics.primitive_evaluations);
}

void fused_row_uses_exact_rational_endpoint_order() {
  auto source = base_solution<Rational>(1);
  source.checkpoint_identity = "fused-exact-interval-source";
  source.sectors = {
      sector<Rational>("0", "0", 0, {Rational(1)})};
  PreparedSparseLocalMultiplierMatrix<Rational> row;
  row.rows = 1;
  row.columns = 1;
  row.exact_identity = "fused-exact-interval-row";
  row.entries.push_back(
      {0, 0,
       prepared_multiplier<Rational>(
           0, 0, {{Rational(1)}}, "1")});
  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {0, 0};

  const std::string pentagon_endpoint =
      "10376293541461622784/130481795891926181249";
  const auto pentagon = diffexp2::integrate_prepared_scalar_row_stored(
      row, source, 0, RealEvaluationPoint::rational("0"),
      RealEvaluationPoint::rational(pentagon_endpoint), options);
  check("fused row accepts the exact rational pentagon tile interval",
        overlaps(pentagon.value.at(0, 0),
                 ComplexBall::from_strings(pentagon_endpoint)));

  const auto lower = RealEvaluationPoint::rational(pentagon_endpoint);
  const auto close_exact = Rational(pentagon_endpoint) + Rational(
      "1/1" + std::string(1000, '0'));
  const auto upper = RealEvaluationPoint::rational(close_exact.str());
  const bool enclosures_overlap = arb_overlaps(
      acb_realref(lower.modulus.raw()), acb_realref(upper.modulus.raw()));
  bool exact_order_accepted = false;
  try {
    (void)diffexp2::integrate_prepared_scalar_row_stored(
        row, source, 0, lower, upper, options);
    exact_order_accepted = true;
  } catch (const NativeIntegrationError&) {
  }
  check("fused row orders distinct rational endpoints exactly when their Acb enclosures overlap",
        enclosures_overlap && exact_order_accepted);
}

ExactEpsilonRational<ComplexBall> exact_ball_rational(
    std::int32_t valuation,
    std::initializer_list<const char*> numerator) {
  ExactEpsilonRational<ComplexBall> value;
  value.zero = false;
  value.valuation = valuation;
  for (const auto* coefficient : numerator)
    value.numerator.push_back(
        ComplexBall::from_strings(coefficient));
  value.denominator.push_back(ComplexBall(1));
  return value;
}

void factorized_ordinary_row_applies_center_ball_once() {
  PreparedPhysicalClearedODE<ComplexBall> equation;
  equation.dimension = 1;
  equation.owner_signature_identity =
      "factorized-line-equation-owner";
  equation.payload_identity = "factorized-line-equation";
  equation.exact_payload_record =
      "factorized-line-equation-record";
  equation.q_lags = {exact_ball_rational(0, {"1"})};
  equation.c_lags.resize(2);
  equation.c_lags[1].push_back(
      PhysicalODEMatrixEntry<ComplexBall>{
          0, 0, exact_ball_rational(0, {"-1"})});

  EpsilonVector initial;
  initial.epsilon = {0, 0};
  initial.dimension = 1;
  initial.coefficients = {
      ComplexBall::from_strings("[1 +/- 0.1]")};
  const auto evolution =
      diffexp2::evolve_ordinary_center_value(
          equation, initial, 12);
  if (!evolution.eligible)
    throw std::runtime_error(
        "factorized line fixture ordinary evolution is ineligible");
  diffexp2::ChartGeometry chart;
  chart.center_exact = "0";
  chart.scale_exact = "1";
  chart.radius_exact = "1";
  chart.radius = ComplexBall(1);
  const auto source =
      diffexp2::ordinary_evolution_local_solution(
          evolution, std::move(chart), {},
          "factorized-line-source");

  std::vector<ComplexBall> identity_kernel(
      source.taylor_width(), ComplexBall(0));
  identity_kernel.front() = ComplexBall(1);
  PreparedSparseLocalMultiplierMatrix<ComplexBall> row;
  row.rows = 1;
  row.columns = 1;
  row.exact_identity = "factorized-line-identity-row";
  row.entries.push_back(
      {0, 0,
       prepared_multiplier<ComplexBall>(
           0, 0, {std::move(identity_kernel)}, "1")});
  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {0, 0};
  const auto lower = RealEvaluationPoint::rational("0");
  const auto upper = RealEvaluationPoint::rational("1/2");
  const auto direct =
      diffexp2::integrate_prepared_scalar_row_stored(
          row, source, 0, lower, upper, options);
  const auto factorized =
      diffexp2::integrate_ordinary_center_stored_row_factorized(
          equation, source, row, 0, lower, upper,
          options, direct);
  const auto capped =
      diffexp2::integrate_ordinary_center_stored_row_factorized(
          equation, source, row, 0, lower, upper,
          options, direct, 0);
  const bool same_value =
      factorized.eligible &&
      overlaps(direct.value.at(0, 0),
               factorized.integral.value.at(0, 0));
  const auto direct_radius = std::stoi(
      direct.value.at(0, 0).real_radius_exponent());
  const auto factorized_radius = factorized.eligible
      ? std::stoi(
            factorized.integral.value.at(0, 0)
                .real_radius_exponent())
      : direct_radius;
  if (!(same_value && factorized.operator_columns == 1 &&
        factorized_radius < direct_radius &&
        factorized.integral.diagnostics.detail ==
            direct.diagnostics.detail))
    std::cerr
        << "factorized-line fixture: eligible="
        << factorized.eligible << "; reason="
        << factorized.reason << "; columns="
        << factorized.operator_columns << "; direct_radius="
        << direct_radius << "; factorized_radius="
        << factorized_radius << "; overlap=" << same_value
        << '\n';
  check("factorized ordinary row integration overlaps the direct finite integral and applies the uncertain center once",
        same_value && factorized.operator_columns == 1 &&
            factorized_radius < direct_radius &&
            factorized.integral.diagnostics.detail ==
                direct.diagnostics.detail);
  check("factorized ordinary row integration obeys its resource cap",
        !capped.eligible &&
            capped.reason.find("column cap") !=
                std::string::npos);
}

void fused_row_charges_upper_halo_only_to_regulated_center_primitive() {
  auto source = base_solution<Rational>(1);
  source.epsilon = {0, 2};
  source.checkpoint_identity = "fused-regulated-halo-source";
  source.sectors = {
      sector<Rational>("-1", "1", 0,
                       {Rational(1), Rational(0), Rational(0)})};
  PreparedSparseLocalMultiplierMatrix<Rational> row;
  row.rows = 1;
  row.columns = 1;
  row.exact_identity = "fused-regulated-halo-row";
  row.entries.push_back(
      {0, 0,
       prepared_multiplier<Rational>(
           0, 0,
           {{Rational(1)}, {Rational(0)}, {Rational(0)}}, "1")});

  StoredLineIntegrationOptions strict;
  strict.delivered_epsilon = {-1, 2};
  bool strict_rejected = false;
  try {
    (void)diffexp2::integrate_prepared_scalar_row_stored(
        row, source, 2, RealEvaluationPoint::rational("0"),
        RealEvaluationPoint::rational("1/2"), strict);
  } catch (const NativeIntegrationError& error) {
    strict_rejected =
        error.code == NativeIntegrationErrorCode::IncompleteEpsilonWindow;
  }

  auto observable = strict;
  observable.required_complete_max = 0;
  const auto clipped = diffexp2::integrate_prepared_scalar_row_stored(
      row, source, 2, RealEvaluationPoint::rational("0"),
      RealEvaluationPoint::rational("1/2"), observable);
  check("fused observable returns the exact complete window after a regulated centre primitive consumes one upper row",
        strict_rejected && clipped.value.epsilon.min_power == -1 &&
            clipped.value.epsilon.complete_max == 1);

  source.sectors = {
      sector<Rational>("0", "0", 0,
                       {Rational(1), Rational(0), Rational(0)})};
  const auto ordinary = diffexp2::integrate_prepared_scalar_row_stored(
      row, source, 2, RealEvaluationPoint::rational("0"),
      RealEvaluationPoint::rational("1/2"), observable);
  check("fused observable preserves the projected upper row for an ordinary centre primitive",
        ordinary.value.epsilon.complete_max == 2);
}

template <typename Scalar>
LocalSolution<Scalar> cross_cell_cancellation_source(const Scalar& positive,
                                                      const Scalar& negative) {
  auto source = base_solution<Scalar>(1);
  source.taylor_complete_max = 1;
  source.checkpoint_identity = "fused-cross-cell-source";
  source.sectors = {
      sector<Scalar>("-1", "0", 0,
                     {positive, Scalar(0)}),
      sector<Scalar>("-2", "0", 0,
                     {Scalar(0), negative})};
  return source;
}

template <typename Scalar>
PreparedSparseLocalMultiplierMatrix<Scalar> cross_cell_identity_row() {
  PreparedSparseLocalMultiplierMatrix<Scalar> row;
  row.rows = 1;
  row.columns = 1;
  row.exact_identity = "fused-cross-cell-identity";
  row.entries.push_back(
      {0, 0,
       prepared_multiplier<Scalar>(
           0, 0, {{Scalar(1), Scalar(0)}}, "1")});
  return row;
}

void fused_row_cancels_divergence_across_sector_taylor_cells() {
  const auto source = cross_cell_cancellation_source<Rational>(
      Rational(1), Rational(-1));
  const auto row = cross_cell_identity_row<Rational>();
  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {0, 0};
  const auto [materialized, fused] = materialized_and_fused(
      row, source, 0, RealEvaluationPoint::rational("0"),
      RealEvaluationPoint::rational("1/2"), options);
  check("fused row performs the divergent cancellation only after distinct sector/Taylor cells meet globally",
        same_integral_values<Rational>(materialized, fused) &&
            fused.value.at(0, 0).is_zero() &&
            materialized.diagnostics.cancelled_divergent_groups == 1 &&
            fused.diagnostics.cancelled_divergent_groups == 1 &&
            materialized.diagnostics.zero_groups_skipped ==
                fused.diagnostics.zero_groups_skipped);
}

void fused_acb_row_matches_bounded_cross_cell_cancellation() {
  ComplexBall uncertain(1);
  // Exercise the same l1-scale distinction through the production fused
  // row path: the residual is above tolerance times one contribution but
  // below tolerance times the two-contribution cancellation norm.
  arb_add_error_2exp_si(acb_realref(uncertain.raw()), -66);
  const auto source = cross_cell_cancellation_source<ComplexBall>(
      uncertain, ComplexBall(-1));
  const auto row = cross_cell_identity_row<ComplexBall>();
  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {0, 0};
  options.divergent_cancellation =
      StoredLineIntegrationOptions::BoundedDivergentCancellation{
          diffexp2::Magnitude::decimal("1e-20"), "1e-20",
          "fused-cross-cell-test"};
  const auto [materialized, fused] = materialized_and_fused(
      row, source, 0, RealEvaluationPoint::rational("0"),
      RealEvaluationPoint::rational("1/2"), options);
  check("fused Acb row preserves bounded relative cancellation across cells",
        same_integral_values<ComplexBall>(materialized, fused) &&
            materialized.diagnostics.bounded_cancelled_divergent_coefficients ==
                1 &&
            fused.diagnostics.bounded_cancelled_divergent_coefficients == 1 &&
            materialized.diagnostics.divergent_cancellation_mode ==
                fused.diagnostics.divergent_cancellation_mode &&
            materialized.diagnostics.divergent_cancellation_provenance ==
                fused.diagnostics.divergent_cancellation_provenance);
}

void fused_acb_row_keeps_preprojection_cancellation_scale() {
  auto source = base_solution<ComplexBall>(2);
  ComplexBall uncertain(1);
  arb_add_error_2exp_si(acb_realref(uncertain.raw()), -66);
  source.sectors = {sector<ComplexBall>(
      "-1", "0", 0, {uncertain, ComplexBall(-1)})};

  PreparedSparseLocalMultiplierMatrix<ComplexBall> row;
  row.rows = 1;
  row.columns = 2;
  row.exact_identity = "fused-preprojection-l1-row";
  row.entries.push_back(
      {0, 0, prepared_multiplier<ComplexBall>(
          0, 0, {{ComplexBall(1)}}, "column-zero")});
  row.entries.push_back(
      {0, 1, prepared_multiplier<ComplexBall>(
          0, 0, {{ComplexBall(1)}}, "column-one")});

  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {0, 0};
  options.divergent_cancellation =
      StoredLineIntegrationOptions::BoundedDivergentCancellation{
          diffexp2::Magnitude::decimal("1e-20"), "1e-20",
          "fused-preprojection-l1-test"};
  const auto fused = diffexp2::integrate_prepared_scalar_row_stored(
      row, source, 0, RealEvaluationPoint::rational("0"),
      RealEvaluationPoint::rational("1/2"), options);
  check("fused Acb row retains the l1 scale before matrix projection cancellation",
        fused.value.at(0, 0).is_zero() &&
            fused.diagnostics.bounded_cancelled_divergent_coefficients == 1 &&
            fused.diagnostics.cancelled_divergent_groups == 1);
}

void fused_sparse_row_ignores_unselected_sector_payloads() {
  auto source = base_solution<Rational>(2);
  source.checkpoint_identity = "fused-sparse-stress-source";
  constexpr std::size_t sector_count = 256;
  source.sectors.reserve(sector_count);
  for (std::size_t sector_number = 0; sector_number < sector_count;
       ++sector_number) {
    auto current = sector<Rational>(
        Rational(std::to_string(sector_number + 1) + "/2").str(), "0", 0,
        std::vector<Rational>(source.sector_size(), Rational(0)));
    current.coefficients[0] = Rational(1);
    if (sector_number + 1 == sector_count)
      current.coefficients[1] = Rational(7);
    source.sectors.push_back(std::move(current));
  }
  const auto row = [&] {
    PreparedSparseLocalMultiplierMatrix<Rational> result;
    result.rows = 1;
    result.columns = 2;
    result.exact_identity = "fused-sparse-stress-row";
    result.entries.push_back(
        {0, 1,
         prepared_multiplier<Rational>(
             0, 0, {{Rational(1)}}, "1")});
    return result;
  }();
  const auto projected = diffexp2::apply_prepared_scalar_row_window(
      row, source, 0, "fused-sparse-stress-materialized");
  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {0, 0};
  const auto [materialized, fused] = materialized_and_fused(
      row, source, 0, RealEvaluationPoint::rational("1/4"),
      RealEvaluationPoint::rational("1/2"), options);
  check("fused sparse row matches materialization while only one of 256 source sectors is selected",
        projected.has_value() && projected->sectors.size() == 1 &&
            same_integral_values<Rational>(materialized, fused) &&
            fused.diagnostics.input_monomial_cells == 1);
  check("fused work metadata is a 16-byte POD with no owned exact-tag strings",
        std::is_trivially_copyable_v<FusedWorkItem> &&
            sizeof(FusedWorkItem) == 16);
}

void fused_missing_rim_preflight_matches_materialized_error_order() {
  auto source = base_solution<Rational>(1);
  source.checkpoint_identity = "fused-missing-rim-source";
  source.sectors = {
      sector<Rational>("-1", "0", 0, {Rational(1)}),
      sector<Rational>("-1/2", "0", 0, {Rational(1)})};
  PreparedSparseLocalMultiplierMatrix<Rational> row;
  row.rows = 1;
  row.columns = 1;
  row.exact_identity = "fused-missing-rim-row";
  row.entries.push_back(
      {0, 0,
       prepared_multiplier<Rational>(
           0, 0, {{Rational(1)}}, "1")});
  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {0, 0};
  const auto projected = diffexp2::apply_prepared_scalar_row_window(
      row, source, 0, "fused-missing-rim-materialized");
  bool materialized_missing_rim = false;
  bool fused_missing_rim = false;
  try {
    (void)diffexp2::integrate_stored_local_line(
        *projected, RealEvaluationPoint::rational("-1/2"),
        RealEvaluationPoint::rational("0"), options);
  } catch (const NativeIntegrationError& error) {
    materialized_missing_rim =
        error.code == NativeIntegrationErrorCode::MissingBranchPrescription;
  }
  try {
    (void)diffexp2::integrate_prepared_scalar_row_stored(
        row, source, 0, RealEvaluationPoint::rational("-1/2"),
        RealEvaluationPoint::rational("0"), options);
  } catch (const NativeIntegrationError& error) {
    fused_missing_rim =
        error.code == NativeIntegrationErrorCode::MissingBranchPrescription;
  }
  check("fused row keeps the global missing-rim preflight ahead of a divergent center group",
        projected.has_value() && materialized_missing_rim &&
            fused_missing_rim);
}

void rational_endpoint_identity_survives_precision_change() {
  auto solution = base_solution<Rational>(1);
  solution.sectors = {
      sector<Rational>("0", "0", 0, {Rational(1)})};
  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {0, 0};

  ComplexBall::set_precision(128);
  const auto lower = RealEvaluationPoint::rational("1/101");
  const auto upper = RealEvaluationPoint::rational("2/101");
  auto inconsistent = lower;
  inconsistent.modulus = ComplexBall::from_strings("2/101");

  ComplexBall::set_precision(1024);
  const auto accepted = diffexp2::integrate_stored_local_line(
      solution, lower, upper, options);
  bool inconsistent_rejected = false;
  try {
    (void)diffexp2::integrate_stored_local_line(
        solution, inconsistent, upper, options);
  } catch (const NativeIntegrationError& error) {
    inconsistent_rejected =
        error.code == NativeIntegrationErrorCode::InvalidInterval;
  }
  check("exact rational endpoint identity survives a higher-precision consumer lease",
        overlaps(accepted.value.at(0, 0),
                 ComplexBall::from_strings("1/101")) &&
            inconsistent_rejected);
  ComplexBall::set_precision(512);
}

}  // namespace

int main() {
  ComplexBall::set_precision(512);
  grouped_center_and_component_reuse();
  acb_cancellation_is_exact_only();
  acb_bounded_divergent_cancellation_is_explicit_and_relative();
  stored_branch_selection();
  upper_halo_and_first_unseen_fail_loudly();
  input_error_envelope_is_not_discarded();
  fused_rational_row_matches_materialized_projection();
  fused_row_uses_exact_rational_endpoint_order();
  factorized_ordinary_row_applies_center_ball_once();
  fused_row_charges_upper_halo_only_to_regulated_center_primitive();
  fused_row_cancels_divergence_across_sector_taylor_cells();
  fused_acb_row_matches_bounded_cross_cell_cancellation();
  fused_acb_row_keeps_preprojection_cancellation_scale();
  fused_sparse_row_ignores_unselected_sector_payloads();
  fused_missing_rim_preflight_matches_materialized_error_order();
  rational_endpoint_identity_survives_precision_change();
  std::cout << "Results: " << passed << " / " << (passed + failed)
            << " tests passed\n";
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
