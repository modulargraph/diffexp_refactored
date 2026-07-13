#include "diffexp2/scalar_row_endpoint.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using diffexp2::ComplexBall;
using diffexp2::EndpointLimitOptions;
using diffexp2::EndpointLimitResult;
using diffexp2::ExactScalarDescriptor;
using diffexp2::LocalSector;
using diffexp2::LocalSolution;
using diffexp2::NativeIntegrationError;
using diffexp2::NativeIntegrationErrorCode;
using diffexp2::PreparedRationalTaylorMultiplier;
using diffexp2::PreparedScalarRowEndpointMode;
using diffexp2::PreparedSparseLocalMultiplierMatrix;
using diffexp2::Rational;

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

std::size_t index(const auto& solution, std::size_t epsilon,
                  std::size_t taylor, std::uint32_t component) {
  return ((epsilon * solution.taylor_width() + taylor) *
          solution.dimension) + component;
}

template <typename Scalar>
Scalar scalar(const std::string& value) {
  if constexpr (std::is_same_v<Scalar, Rational>)
    return Rational(value);
  else
    return ComplexBall::from_strings(value);
}

template <typename Scalar>
LocalSolution<Scalar> local(std::int32_t epsilon_min,
                            std::int32_t epsilon_max,
                            std::uint32_t taylor_max,
                            std::uint32_t dimension,
                            const std::vector<std::pair<
                                std::string, std::uint32_t>>& tags,
                            const std::string& b = "0") {
  LocalSolution<Scalar> result;
  result.chart.center_exact = "0";
  result.chart.scale_exact = "1";
  result.chart.radius = ComplexBall::from_strings("2");
  result.epsilon = {epsilon_min, epsilon_max};
  result.taylor_complete_max = taylor_max;
  result.dimension = dimension;
  result.checkpoint_identity = "scalar-row-endpoint-source";
  result.prescriptions.push_back({"t", -1, 1, 1});
  for (const auto& [a, log_power] : tags) {
    LocalSector<Scalar> sector;
    sector.a = ExactScalarDescriptor::rational(a);
    sector.b = ExactScalarDescriptor::rational(b);
    sector.log_power = log_power;
    sector.coefficients.assign(result.sector_size(), scalar<Scalar>("0"));
    result.sectors.push_back(std::move(sector));
  }
  diffexp2::validate_local_solution(result, false);
  return result;
}

template <typename Scalar>
PreparedRationalTaylorMultiplier<Scalar> multiplier(
    std::size_t epsilon_width, std::size_t taylor_width,
    std::int32_t epsilon_shift, std::uint32_t pole,
    const std::string& constant, const std::string& linear_taylor = "0") {
  PreparedRationalTaylorMultiplier<Scalar> result;
  result.epsilon_shift = epsilon_shift;
  result.center_pole_order = pole;
  result.exact_identity = "test-multiplier";
  result.kernels.assign(
      epsilon_width,
      std::vector<Scalar>(taylor_width, scalar<Scalar>("0")));
  result.kernels[0][0] = scalar<Scalar>(constant);
  if (taylor_width > 1)
    result.kernels[0][1] = scalar<Scalar>(linear_taylor);
  return result;
}

template <typename Scalar>
PreparedSparseLocalMultiplierMatrix<Scalar> row(
    const LocalSolution<Scalar>& source,
    const std::vector<std::tuple<std::uint32_t, std::int32_t,
                                 std::uint32_t, std::string,
                                 std::string>>& entries) {
  PreparedSparseLocalMultiplierMatrix<Scalar> result;
  result.rows = 1;
  result.columns = source.dimension;
  result.exact_identity = "scalar-row-endpoint-test-row";
  for (const auto& [column, shift, pole, constant, linear] : entries)
    result.entries.push_back({
        0, column,
        multiplier<Scalar>(source.epsilon.width(), source.taylor_width(),
                           shift, pole, constant, linear)});
  return result;
}

template <typename Scalar>
LocalSolution<Scalar> exact_zero_like(
    const LocalSolution<Scalar>& source,
    const std::string& checkpoint) {
  auto result = diffexp2::local_algebra_detail::with_selected_component(
      source, 0);
  for (auto& sector : result.sectors)
    std::fill(sector.coefficients.begin(), sector.coefficients.end(),
              scalar<Scalar>("0"));
  result.checkpoint_identity = checkpoint;
  return result;
}

template <typename Scalar>
struct OracleResult {
  EndpointLimitResult endpoint;
  LocalSolution<Scalar> projected;
};

template <typename Scalar>
OracleResult<Scalar> materialized_oracle(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& matrix,
    const LocalSolution<Scalar>& source,
    std::int32_t source_top_valid,
    const EndpointLimitOptions& options) {
  std::int32_t top_valid = matrix.entries.empty()
      ? source_top_valid : diffexp2::kCompleteInfinity;
  for (const auto& entry : matrix.entries) {
    const auto shifted = source_top_valid == diffexp2::kCompleteInfinity
        ? diffexp2::kCompleteInfinity
        : diffexp2::local_algebra_detail::checked_i32(
              static_cast<std::int64_t>(source_top_valid) +
                  entry.multiplier.epsilon_shift,
              "oracle shifted validity");
    top_valid = std::min(top_valid, shifted);
  }
  auto applied = diffexp2::apply_prepared_sparse_local_matrix(
      matrix, source, "materialized-oracle");
  auto projected = applied.has_value()
      ? std::move(*applied)
      : exact_zero_like(source, "materialized-oracle");
  top_valid = std::min(top_valid, projected.epsilon.complete_max);
  if (top_valid < projected.epsilon.min_power)
    throw std::domain_error(
        "rational-row application has no valid output epsilon coefficient");
  if (top_valid < projected.epsilon.complete_max)
    projected = diffexp2::restrict_local_epsilon_frame_strict_lower(
        projected, projected.epsilon.min_power, top_valid,
        "materialized-oracle");
  return {diffexp2::endpoint_sector_limit(projected, options),
          std::move(projected)};
}

bool same_ball(const ComplexBall& left, const ComplexBall& right) {
  return acb_equal(left.raw(), right.raw());
}

bool same_endpoint(const EndpointLimitResult& left,
                   const EndpointLimitResult& right) {
  if (left.values.size() != right.values.size() ||
      left.dropped_regulated_sectors != right.dropped_regulated_sectors ||
      left.cancelled_divergent_coefficients !=
          right.cancelled_divergent_coefficients ||
      left.imaginary_sign != right.imaginary_sign)
    return false;
  for (std::size_t component = 0; component < left.values.size();
       ++component) {
    const auto& a = left.values[component];
    const auto& b = right.values[component];
    if (a.min_power() != b.min_power() ||
        a.complete_max() != b.complete_max() ||
        a.coefficients().size() != b.coefficients().size())
      return false;
    for (std::size_t coefficient = 0;
         coefficient < a.coefficients().size(); ++coefficient)
      if (!same_ball(a.coefficients()[coefficient],
                     b.coefficients()[coefficient]))
        return false;
  }
  return true;
}

template <typename Scalar>
bool same_tags(const LocalSolution<Scalar>& left,
               const LocalSolution<Scalar>& right) {
  if (left.sectors.size() != right.sectors.size() ||
      left.epsilon.min_power != right.epsilon.min_power ||
      left.epsilon.complete_max != right.epsilon.complete_max)
    return false;
  for (std::size_t sector = 0; sector < left.sectors.size(); ++sector) {
    const auto& a = left.sectors[sector];
    const auto& b = right.sectors[sector];
    if (a.a.domain != b.a.domain || a.a.canonical != b.a.canonical ||
        a.b.domain != b.b.domain || a.b.canonical != b.b.canonical ||
        a.log_power != b.log_power)
      return false;
  }
  return true;
}

template <typename Scalar>
bool same_live_coefficients(const LocalSolution<Scalar>& capped,
                            const LocalSolution<Scalar>& full) {
  if (!same_tags(capped, full) ||
      capped.taylor_complete_max > full.taylor_complete_max)
    return false;
  for (std::size_t sector = 0; sector < capped.sectors.size(); ++sector)
    for (std::size_t epsilon = 0;
         epsilon < capped.epsilon.width(); ++epsilon)
      for (std::size_t taylor = 0;
           taylor < capped.taylor_width(); ++taylor) {
        const auto& a = capped.sectors[sector].coefficients[
            index(capped, epsilon, taylor, 0)];
        const auto& b = full.sectors[sector].coefficients[
            index(full, epsilon, taylor, 0)];
        if constexpr (std::is_same_v<Scalar, Rational>) {
          if (!(a == b)) return false;
        } else {
          if (!same_ball(a, b)) return false;
        }
      }
  return true;
}

void cross_component_and_distinct_cell_tests() {
  EndpointLimitOptions options;
  options.approach_direction = 1;
  options.imaginary_sign = -1;

  auto source = local<Rational>(0, 1, 4, 2, {{"-1", 0}});
  auto& sector = source.sectors.front();
  sector.coefficients[index(source, 0, 0, 0)] = Rational(3);
  sector.coefficients[index(source, 0, 0, 1)] = Rational(-3);
  sector.coefficients[index(source, 0, 1, 0)] = Rational(5);
  sector.coefficients[index(source, 0, 1, 1)] = Rational(7);
  sector.coefficients[index(source, 1, 0, 0)] = Rational(2);
  sector.coefficients[index(source, 1, 0, 1)] = Rational(-2);
  sector.coefficients[index(source, 1, 1, 0)] = Rational(11);
  sector.coefficients[index(source, 1, 1, 1)] = Rational(13);
  const auto matrix = row<Rational>(
      source, {{0, 0, 0, "1", "0"}, {1, 0, 0, "1", "0"}});
  const auto oracle = materialized_oracle(matrix, source, 1, options);
  const auto fused = diffexp2::centered_prepared_scalar_row_endpoint_limit(
      matrix, source, 1, options, "cross-component-fused");
  check("cross-component divergent cancellation preserves finite remainder",
        fused.mode == PreparedScalarRowEndpointMode::
                          CenteredB0RationalTaylorLivenessV1 &&
            fused.projected.taylor_complete_max == 1 &&
            same_endpoint(fused.endpoint, oracle.endpoint) &&
            same_live_coefficients(fused.projected, oracle.projected) &&
            same_ball(fused.endpoint.values.front().coefficient(0),
                      ComplexBall(12)) &&
            same_ball(fused.endpoint.values.front().coefficient(1),
                      ComplexBall(24)));

  auto distinct = local<Rational>(0, 0, 5, 1,
                                  {{"-2", 0}, {"-1", 0}});
  distinct.sectors[0].coefficients[index(distinct, 0, 1, 0)] = Rational(4);
  distinct.sectors[0].coefficients[index(distinct, 0, 2, 0)] = Rational(6);
  distinct.sectors[1].coefficients[index(distinct, 0, 0, 0)] = Rational(-4);
  distinct.sectors[1].coefficients[index(distinct, 0, 1, 0)] = Rational(8);
  const auto distinct_row = row<Rational>(
      distinct, {{0, 0, 0, "1", "0"}});
  const auto distinct_oracle = materialized_oracle(
      distinct_row, distinct, 0, options);
  const auto distinct_fused =
      diffexp2::centered_prepared_scalar_row_endpoint_limit(
          distinct_row, distinct, 0, options, "distinct-cell-fused");
  check("distinct sector/Taylor cells cancel by absolute-power class",
        distinct_fused.projected.taylor_complete_max == 2 &&
            same_endpoint(distinct_fused.endpoint,
                          distinct_oracle.endpoint) &&
            same_live_coefficients(distinct_fused.projected,
                                   distinct_oracle.projected) &&
            same_ball(distinct_fused.endpoint.values.front().coefficient(0),
                      ComplexBall(14)));
}

void mixed_shift_pole_and_fallback_tests() {
  EndpointLimitOptions options;
  options.imaginary_sign = -1;
  auto source = local<Rational>(-2, 3, 5, 2, {{"0", 0}});
  auto& sector = source.sectors.front();
  const auto ep = [&](std::int32_t power) {
    return static_cast<std::size_t>(power - source.epsilon.min_power);
  };
  sector.coefficients[index(source, ep(0), 0, 0)] = Rational(5);
  sector.coefficients[index(source, ep(-2), 0, 1)] = Rational(-5);
  sector.coefficients[index(source, ep(1), 0, 0)] = Rational(7);
  sector.coefficients[index(source, ep(-1), 0, 1)] = Rational(-7);
  sector.coefficients[index(source, ep(0), 1, 0)] = Rational(2);
  sector.coefficients[index(source, ep(-2), 1, 1)] = Rational(3);
  const auto matrix = row<Rational>(
      source, {{0, -1, 1, "1", "1"}, {1, 1, 1, "1", "0"}});
  const auto oracle = materialized_oracle(matrix, source, 2, options);
  const auto fused = diffexp2::centered_prepared_scalar_row_endpoint_limit(
      matrix, source, 2, options, "mixed-shift-pole-fused");
  check("nonzero pole and mixed epsilon shifts preserve natural frame/order",
        fused.mode == PreparedScalarRowEndpointMode::
                          CenteredB0RationalTaylorLivenessV1 &&
            fused.projected.epsilon.min_power == -3 &&
            fused.projected.epsilon.complete_max == 1 &&
            fused.projected_top_valid == 1 &&
            same_endpoint(fused.endpoint, oracle.endpoint) &&
            same_live_coefficients(fused.projected, oracle.projected));

  auto regulated = local<Rational>(0, 1, 4, 1, {{"-2", 0}}, "1");
  regulated.sectors.front().coefficients[index(regulated, 0, 0, 0)] =
      Rational(9);
  const auto regulated_row = row<Rational>(
      regulated, {{0, 0, 0, "1", "0"}});
  const auto regulated_oracle = materialized_oracle(
      regulated_row, regulated, 1, options);
  const auto regulated_fused =
      diffexp2::centered_prepared_scalar_row_endpoint_limit(
          regulated_row, regulated, 1, options, "regulated-fallback");
  check("b!=0 uses strict materialized fallback and preserves drop count",
        regulated_fused.mode ==
                PreparedScalarRowEndpointMode::MaterializedFallbackV1 &&
            regulated_fused.fallback_reason ==
                "nonzero-or-unsupported-regulator" &&
            same_endpoint(regulated_fused.endpoint,
                          regulated_oracle.endpoint) &&
            regulated_fused.endpoint.dropped_regulated_sectors == 1);

  auto incomplete = local<Rational>(0, 0, 1, 1, {{"-4", 0}});
  incomplete.sectors.front().coefficients[index(incomplete, 0, 0, 0)] =
      Rational(1);
  const auto incomplete_row = row<Rational>(
      incomplete, {{0, 0, 0, "1", "0"}});
  std::optional<NativeIntegrationErrorCode> oracle_error;
  std::optional<NativeIntegrationErrorCode> fused_error;
  try {
    (void)materialized_oracle(incomplete_row, incomplete, 0, options);
  } catch (const NativeIntegrationError& error) {
    oracle_error = error.code;
  }
  try {
    (void)diffexp2::centered_prepared_scalar_row_endpoint_limit(
        incomplete_row, incomplete, 0, options, "incomplete-fallback");
  } catch (const NativeIntegrationError& error) {
    fused_error = error.code;
  }
  check("incomplete tower preserves full-oracle failure ordering/code",
        oracle_error == NativeIntegrationErrorCode::IncompleteTaylorWindow &&
            fused_error == oracle_error);
}

void zero_and_acb_tests() {
  EndpointLimitOptions options;
  options.imaginary_sign = -1;
  auto zero_source = local<Rational>(0, 1, 4, 1, {{"-1", 0}});
  auto active_zero = row<Rational>(
      zero_source, {{0, 0, 0, "0", "0"}});
  const auto zero_oracle = materialized_oracle(
      active_zero, zero_source, 1, options);
  const auto zero_fused =
      diffexp2::centered_prepared_scalar_row_endpoint_limit(
          active_zero, zero_source, 1, options, "active-zero-fused");
  check("active exact-zero row retains natural frame and fast tag support",
        zero_fused.mode == PreparedScalarRowEndpointMode::
                               CenteredB0RationalTaylorLivenessV1 &&
            same_endpoint(zero_fused.endpoint, zero_oracle.endpoint) &&
            same_tags(zero_fused.projected, zero_oracle.projected));

  PreparedSparseLocalMultiplierMatrix<Rational> empty;
  empty.rows = 1;
  empty.columns = 1;
  empty.exact_identity = "empty-row";
  const auto empty_oracle = materialized_oracle(
      empty, zero_source, 1, options);
  const auto empty_fused =
      diffexp2::centered_prepared_scalar_row_endpoint_limit(
          empty, zero_source, 1, options, "empty-fallback");
  check("empty row preserves materialized zero representative",
        empty_fused.mode ==
                PreparedScalarRowEndpointMode::MaterializedFallbackV1 &&
            empty_fused.fallback_reason == "empty-row" &&
            same_endpoint(empty_fused.endpoint, empty_oracle.endpoint) &&
            same_tags(empty_fused.projected, empty_oracle.projected));

  auto acb = local<ComplexBall>(0, 0, 4, 2, {{"-1", 0}});
  acb.sectors.front().coefficients[index(acb, 0, 0, 0)] = ComplexBall(3);
  acb.sectors.front().coefficients[index(acb, 0, 0, 1)] = ComplexBall(-3);
  acb.sectors.front().coefficients[index(acb, 0, 1, 0)] = ComplexBall(5);
  acb.sectors.front().coefficients[index(acb, 0, 1, 1)] = ComplexBall(7);
  const auto acb_row = row<ComplexBall>(
      acb, {{0, 0, 0, "1", "0"}, {1, 0, 0, "1", "0"}});
  const auto acb_oracle = materialized_oracle(acb_row, acb, 0, options);
  const auto acb_fused =
      diffexp2::centered_prepared_scalar_row_endpoint_limit(
          acb_row, acb, 0, options, "acb-singleton-fused");
  check("Acb singleton-zero cancellation preserves materialized order",
        acb_fused.mode == PreparedScalarRowEndpointMode::
                               CenteredB0RationalTaylorLivenessV1 &&
            same_endpoint(acb_fused.endpoint, acb_oracle.endpoint) &&
            same_live_coefficients(acb_fused.projected,
                                   acb_oracle.projected));

  auto uncertain = acb;
  auto ambiguous = ComplexBall(-3);
  arb_add_error_2exp_si(acb_realref(ambiguous.raw()), -100);
  uncertain.sectors.front().coefficients[index(uncertain, 0, 0, 1)] =
      ambiguous;
  std::optional<NativeIntegrationErrorCode> oracle_error;
  std::optional<NativeIntegrationErrorCode> fused_error;
  try {
    (void)materialized_oracle(acb_row, uncertain, 0, options);
  } catch (const NativeIntegrationError& error) {
    oracle_error = error.code;
  }
  try {
    (void)diffexp2::centered_prepared_scalar_row_endpoint_limit(
        acb_row, uncertain, 0, options, "acb-uncertain-fused");
  } catch (const NativeIntegrationError& error) {
    fused_error = error.code;
  }
  check("Acb enclosure containing zero is never accepted as cancellation",
        oracle_error ==
                NativeIntegrationErrorCode::UncertifiedCancellation &&
            fused_error == oracle_error);
}

}  // namespace

int main() {
  ComplexBall::set_precision(512);
  cross_component_and_distinct_cell_tests();
  mixed_shift_pole_and_fallback_tests();
  zero_and_acb_tests();
  std::cout << "Results: " << passed << " / " << (passed + failed)
            << " tests passed\n";
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
