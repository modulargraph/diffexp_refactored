#include "diffexp2/local_algebra.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {

using diffexp2::ExactScalarDescriptor;
using diffexp2::ComplexBall;
using diffexp2::LocalSector;
using diffexp2::LocalSolution;
using diffexp2::PreparedEpsilonMultiplier;
using diffexp2::PreparedRationalTaylorMultiplier;
using diffexp2::PreparedSparseLocalMultiplierMatrix;
using diffexp2::Rational;

std::size_t index(std::size_t ei, std::size_t n, std::uint32_t c,
                  std::size_t tw, std::uint32_t d) {
  return ((ei * tw + n) * d) + c;
}

LocalSolution<Rational> sample() {
  LocalSolution<Rational> out;
  out.chart.center_exact = "0";
  out.chart.scale_exact = "1";
  out.chart.radius = diffexp2::ComplexBall::from_strings("2");
  out.epsilon = {-1, 1};
  out.taylor_complete_max = 2;
  out.dimension = 2;
  out.checkpoint_identity = "local-algebra-input";
  out.prescriptions.push_back({"x", 1, 1, 1});
  LocalSector<Rational> sector;
  sector.a = ExactScalarDescriptor::rational("1/2");
  sector.b = ExactScalarDescriptor::rational("0");
  sector.log_power = 0;
  sector.coefficients.assign(out.sector_size(), Rational(0));
  for (std::size_t ei = 0; ei < out.epsilon.width(); ++ei)
    for (std::size_t n = 0; n < out.taylor_width(); ++n)
      for (std::uint32_t c = 0; c < out.dimension; ++c)
        sector.coefficients[index(ei, n, c, out.taylor_width(), out.dimension)] =
            Rational(std::to_string(100 * ei + 10 * n + c + 1));
  out.sectors.push_back(std::move(sector));
  diffexp2::validate_local_solution(out, false);
  return out;
}

PreparedRationalTaylorMultiplier<Rational> multiplier(
    std::int32_t shift, std::uint32_t pole) {
  PreparedRationalTaylorMultiplier<Rational> out;
  out.epsilon_shift = shift;
  out.center_pole_order = pole;
  out.exact_identity = "(1+2t)+3eps";
  out.kernels = {
      {Rational(1), Rational(2), Rational(0)},
      {Rational(3), Rational(0), Rational(0)},
      {Rational(0), Rational(0), Rational(0)}};
  return out;
}

bool equal(const Rational& value, const std::string& expected) {
  return value == Rational(expected);
}

bool same_rational_local(const LocalSolution<Rational>& left,
                         const LocalSolution<Rational>& right) {
  if (left.chart.center_exact != right.chart.center_exact ||
      left.chart.scale_exact != right.chart.scale_exact ||
      left.epsilon.min_power != right.epsilon.min_power ||
      left.epsilon.complete_max != right.epsilon.complete_max ||
      left.taylor_complete_max != right.taylor_complete_max ||
      left.dimension != right.dimension)
    return false;
  std::vector<const LocalSector<Rational>*> left_material;
  std::vector<const LocalSector<Rational>*> right_material;
  for (const auto& sector : left.sectors)
    if (std::any_of(sector.coefficients.begin(), sector.coefficients.end(),
                    [](const auto& value) { return !value.is_zero(); }))
      left_material.push_back(&sector);
  for (const auto& sector : right.sectors)
    if (std::any_of(sector.coefficients.begin(), sector.coefficients.end(),
                    [](const auto& value) { return !value.is_zero(); }))
      right_material.push_back(&sector);
  if (left_material.size() != right_material.size()) return false;
  for (std::size_t sector = 0; sector < left_material.size(); ++sector) {
    const auto& a = *left_material[sector];
    const auto& b = *right_material[sector];
    if (a.a.canonical != b.a.canonical ||
        a.b.canonical != b.b.canonical ||
        a.log_power != b.log_power ||
        a.coefficients != b.coefficients)
      return false;
  }
  return true;
}

LocalSolution<ComplexBall> acb_sample() {
  const auto exact = sample();
  LocalSolution<ComplexBall> out;
  out.chart = exact.chart;
  out.epsilon = exact.epsilon;
  out.taylor_complete_max = exact.taylor_complete_max;
  out.dimension = exact.dimension;
  out.prescriptions = exact.prescriptions;
  out.checkpoint_identity = "local-algebra-acb-input";
  for (const auto& source : exact.sectors) {
    LocalSector<ComplexBall> sector;
    sector.a = source.a;
    sector.b = source.b;
    sector.log_power = source.log_power;
    for (const auto& coefficient : source.coefficients)
      sector.coefficients.push_back(
          ComplexBall::from_strings(coefficient.str()));
    out.sectors.push_back(std::move(sector));
  }
  diffexp2::validate_local_solution(out, false);
  return out;
}

PreparedRationalTaylorMultiplier<ComplexBall> acb_multiplier(
    std::int32_t shift, std::uint32_t pole) {
  const auto exact = multiplier(shift, pole);
  PreparedRationalTaylorMultiplier<ComplexBall> out;
  out.epsilon_shift = exact.epsilon_shift;
  out.center_pole_order = exact.center_pole_order;
  out.exact_identity = exact.exact_identity;
  for (const auto& exact_kernel : exact.kernels) {
    std::vector<ComplexBall> kernel;
    for (const auto& coefficient : exact_kernel)
      kernel.push_back(ComplexBall::from_strings(coefficient.str()));
    out.kernels.push_back(std::move(kernel));
  }
  return out;
}

bool overlapping_acb_local(const LocalSolution<ComplexBall>& left,
                           const LocalSolution<ComplexBall>& right) {
  if (left.epsilon.min_power != right.epsilon.min_power ||
      left.epsilon.complete_max != right.epsilon.complete_max ||
      left.taylor_complete_max != right.taylor_complete_max ||
      left.dimension != right.dimension)
    return false;
  std::vector<const LocalSector<ComplexBall>*> left_material;
  std::vector<const LocalSector<ComplexBall>*> right_material;
  for (const auto& sector : left.sectors)
    if (std::any_of(sector.coefficients.begin(), sector.coefficients.end(),
                    [](const auto& value) { return !value.is_zero(); }))
      left_material.push_back(&sector);
  for (const auto& sector : right.sectors)
    if (std::any_of(sector.coefficients.begin(), sector.coefficients.end(),
                    [](const auto& value) { return !value.is_zero(); }))
      right_material.push_back(&sector);
  if (left_material.size() != right_material.size()) return false;
  for (std::size_t sector = 0; sector < left_material.size(); ++sector) {
    const auto& a = *left_material[sector];
    const auto& b = *right_material[sector];
    if (a.a.canonical != b.a.canonical ||
        a.b.canonical != b.b.canonical ||
        a.log_power != b.log_power ||
        a.coefficients.size() != b.coefficients.size())
      return false;
    for (std::size_t coefficient = 0;
         coefficient < a.coefficients.size(); ++coefficient)
      if (!(a.coefficients[coefficient] - b.coefficients[coefficient])
               .contains_zero())
        return false;
  }
  return true;
}

PreparedRationalTaylorMultiplier<Rational> identity_multiplier() {
  auto out = multiplier(0, 0);
  for (auto& kernel : out.kernels)
    std::fill(kernel.begin(), kernel.end(), Rational(0));
  out.kernels.front().front() = Rational(1);
  out.exact_identity = "1";
  return out;
}

LocalSolution<Rational> component_sparse_sample(std::size_t sector_count) {
  auto out = sample();
  out.checkpoint_identity = "component-sparse-input";
  out.sectors.clear();
  out.sectors.reserve(sector_count);
  for (std::size_t sector_index = 0; sector_index < sector_count;
       ++sector_index) {
    LocalSector<Rational> sector;
    sector.a = ExactScalarDescriptor::rational(
        Rational(std::to_string(2 * sector_index + 1) + "/2").str());
    sector.b = ExactScalarDescriptor::rational("0");
    sector.log_power = 0;
    sector.coefficients.assign(out.sector_size(), Rational(0));
    // Component zero is material everywhere, but the projected component
    // one occurs in only the final sector.
    sector.coefficients[index(0, 0, 0, out.taylor_width(), out.dimension)] =
        Rational(1);
    if (sector_index + 1 == sector_count)
      sector.coefficients[index(
          0, 0, 1, out.taylor_width(), out.dimension)] = Rational(7);
    out.sectors.push_back(std::move(sector));
  }
  diffexp2::validate_local_solution(out, false);
  return out;
}

}  // namespace

int main() {
  diffexp2::ComplexBall::set_precision(128);
  const auto input = sample();
  const auto product = diffexp2::multiply_prepared_rational(
      input, multiplier(-1, 1));

  bool ok = product.epsilon.min_power == -2 &&
            product.epsilon.complete_max == 0 &&
            product.taylor_complete_max == 2 &&
            product.dimension == 2 && product.sectors.size() == 1 &&
            product.sectors.front().a.canonical == "-1/2" &&
            product.prescriptions.size() == 1;
  const auto& coefficients = product.sectors.front().coefficients;
  // out(eps index 0,t^1,c0) = in(0,t^1,c0) + 2 in(0,t^0,c0)
  ok = ok && equal(coefficients[index(0, 1, 0, 3, 2)], "13");
  // out(eps index 1,t^0,c0) = in(1,t^0,c0) + 3 in(0,t^0,c0)
  ok = ok && equal(coefficients[index(1, 0, 0, 3, 2)], "104");
  // out(eps index 2,t^2,c1) = in(2,t^2,c1) + 2 in(2,t^1,c1)
  //                              + 3 in(1,t^2,c1)
  ok = ok && equal(coefficients[index(2, 2, 1, 3, 2)], "1012");

  auto epsilon_generic = multiplier(-2, 0);
  for (auto& kernel : epsilon_generic.kernels) {
    const auto constant = kernel.front();
    std::fill(kernel.begin(), kernel.end(), Rational(0));
    kernel.front() = constant;
  }
  PreparedEpsilonMultiplier<Rational> epsilon_specialized;
  epsilon_specialized.epsilon_shift = epsilon_generic.epsilon_shift;
  epsilon_specialized.exact_identity = "1+3eps";
  epsilon_specialized.kernels = {
      Rational(1), Rational(3), Rational(0)};
  const auto epsilon_reference = diffexp2::multiply_prepared_rational(
      input, epsilon_generic, "epsilon-reference");
  const auto epsilon_product = diffexp2::multiply_prepared_epsilon(
      input, epsilon_specialized, "epsilon-specialized");
  ok = ok && same_rational_local(epsilon_reference, epsilon_product) &&
       epsilon_product.epsilon.min_power == -3 &&
       epsilon_product.epsilon.complete_max == -1 &&
       epsilon_product.sectors.front().a.canonical == "1/2";

  PreparedSparseLocalMultiplierMatrix<Rational> matrix;
  matrix.rows = 2;
  matrix.columns = 2;
  matrix.exact_identity = "one-coupling";
  matrix.entries.push_back({1, 0, multiplier(-1, 1)});
  const auto applied = diffexp2::apply_prepared_sparse_local_matrix(
      matrix, input);
  ok = ok && applied.has_value() && applied->dimension == 2 &&
       applied->epsilon.min_power == -2 &&
       applied->epsilon.complete_max == 0;
  if (applied.has_value()) {
    const auto& slab = applied->sectors.front().coefficients;
    ok = ok && equal(slab[index(1, 0, 0, 3, 2)], "0") &&
         equal(slab[index(1, 0, 1, 3, 2)], "104");
  }

  const auto unshifted = diffexp2::multiply_prepared_rational(
      input, multiplier(0, 0));
  const auto combined = diffexp2::combine_local_solutions(
      std::vector<LocalSolution<Rational>>{product, unshifted});
  ok = ok && combined.epsilon.min_power == -2 &&
       combined.epsilon.complete_max == 0;

  PreparedSparseLocalMultiplierMatrix<Rational> empty;
  empty.rows = 2;
  empty.columns = 2;
  ok = ok && !diffexp2::apply_prepared_sparse_local_matrix(
      empty, input).has_value();

  PreparedSparseLocalMultiplierMatrix<Rational> scalar_row;
  scalar_row.rows = 1;
  scalar_row.columns = 2;
  scalar_row.exact_identity = "direct-scalar-row";
  scalar_row.entries.push_back({0, 0, multiplier(-1, 1)});
  scalar_row.entries.push_back({0, 1, multiplier(0, 0)});
  const auto scalar_full = diffexp2::apply_prepared_sparse_local_matrix(
      scalar_row, input, "scalar-full");
  const auto scalar_direct = diffexp2::apply_prepared_scalar_row_window(
      scalar_row, input, -1, "scalar-direct");
  ok = ok && scalar_full.has_value() && scalar_direct.has_value();
  if (scalar_full.has_value() && scalar_direct.has_value()) {
    const auto scalar_restricted =
        diffexp2::restrict_local_epsilon_frame_strict_lower(
            *scalar_full, scalar_full->epsilon.min_power, -1,
            "scalar-restricted");
    ok = ok && same_rational_local(
                   scalar_restricted, *scalar_direct) &&
         scalar_direct->epsilon.min_power == -2 &&
         scalar_direct->epsilon.complete_max == -1 &&
         scalar_direct->dimension == 1;
  }

  PreparedSparseLocalMultiplierMatrix<Rational> zero_guard_row;
  zero_guard_row.rows = 1;
  zero_guard_row.columns = 2;
  zero_guard_row.exact_identity = "direct-zero-frame-guard";
  zero_guard_row.entries.push_back({0, 0, multiplier(-1, 1)});
  auto zero_active = multiplier(-2, 0);
  for (auto& kernel : zero_active.kernels)
    std::fill(kernel.begin(), kernel.end(), Rational(0));
  zero_guard_row.entries.push_back({0, 1, std::move(zero_active)});
  const auto zero_guard_full =
      diffexp2::apply_prepared_sparse_local_matrix(zero_guard_row, input);
  const auto zero_guard_direct =
      diffexp2::apply_prepared_scalar_row_window(
          zero_guard_row, input, 0);
  ok = ok && zero_guard_full.has_value() &&
       zero_guard_direct.has_value() &&
       same_rational_local(*zero_guard_full, *zero_guard_direct) &&
       zero_guard_direct->epsilon.min_power == -3 &&
       zero_guard_direct->epsilon.complete_max == -1;

  PreparedSparseLocalMultiplierMatrix<Rational> all_zero_row;
  all_zero_row.rows = 1;
  all_zero_row.columns = 2;
  all_zero_row.exact_identity = "active-all-zero-row";
  auto all_zero_multiplier = multiplier(-2, 1);
  for (auto& kernel : all_zero_multiplier.kernels)
    std::fill(kernel.begin(), kernel.end(), Rational(0));
  all_zero_row.entries.push_back(
      {0, 1, std::move(all_zero_multiplier)});
  const auto all_zero_direct =
      diffexp2::apply_prepared_scalar_row_window(
          all_zero_row, input, 0, "active-all-zero-direct");
  ok = ok && all_zero_direct.has_value() &&
       all_zero_direct->epsilon.min_power == -3 &&
       all_zero_direct->epsilon.complete_max == -1 &&
       all_zero_direct->sectors.size() == 1 &&
       all_zero_direct->sectors.front().a.canonical == "-1/2" &&
       std::all_of(
           all_zero_direct->sectors.front().coefficients.begin(),
           all_zero_direct->sectors.front().coefficients.end(),
           [](const auto& value) { return value.is_zero(); });

  const auto sparse_input = component_sparse_sample(128);
  PreparedSparseLocalMultiplierMatrix<Rational> sparse_row;
  sparse_row.rows = 1;
  sparse_row.columns = 2;
  sparse_row.exact_identity = "component-sparse-row";
  sparse_row.entries.push_back({0, 1, identity_multiplier()});
  const auto sparse_full = diffexp2::apply_prepared_sparse_local_matrix(
      sparse_row, sparse_input, "component-sparse-full");
  const auto sparse_direct =
      diffexp2::apply_prepared_scalar_row_window(
          sparse_row, sparse_input, sparse_input.epsilon.complete_max,
          "component-sparse-direct");
  ok = ok && sparse_full.has_value() && sparse_direct.has_value() &&
       same_rational_local(*sparse_full, *sparse_direct) &&
       sparse_full->sectors.size() == 128 &&
       sparse_direct->sectors.size() == 1 &&
       sparse_direct->sectors.front().a.canonical == "255/2";

  PreparedSparseLocalMultiplierMatrix<ComplexBall> acb_row;
  acb_row.rows = 1;
  acb_row.columns = 2;
  acb_row.exact_identity = "direct-acb-row";
  acb_row.entries.push_back({0, 0, acb_multiplier(-1, 1)});
  const auto acb_input = acb_sample();
  const auto acb_full = diffexp2::apply_prepared_sparse_local_matrix(
      acb_row, acb_input, "acb-full");
  const auto acb_direct = diffexp2::apply_prepared_scalar_row_window(
      acb_row, acb_input, -1, "acb-direct");
  ok = ok && acb_full.has_value() && acb_direct.has_value();
  if (acb_full.has_value() && acb_direct.has_value()) {
    const auto acb_restricted =
        diffexp2::restrict_local_epsilon_frame_strict_lower(
            *acb_full, acb_full->epsilon.min_power, -1,
            "acb-restricted");
    ok = ok && overlapping_acb_local(acb_restricted, *acb_direct);
  }

  auto ambiguous_input = acb_sample();
  for (auto& sector : ambiguous_input.sectors)
    for (auto& coefficient : sector.coefficients)
      coefficient = ComplexBall(0);
  auto ambiguous = ComplexBall(0);
  arb_add_error_2exp_si(acb_realref(ambiguous.raw()), -80);
  ambiguous_input.sectors.front().coefficients[index(
      0, 0, 0, ambiguous_input.taylor_width(),
      ambiguous_input.dimension)] = ambiguous;
  PreparedSparseLocalMultiplierMatrix<ComplexBall> ambiguous_row;
  ambiguous_row.rows = 1;
  ambiguous_row.columns = 2;
  ambiguous_row.exact_identity = "ambiguous-acb-row";
  auto ambiguous_multiplier = acb_multiplier(0, 0);
  for (auto& kernel : ambiguous_multiplier.kernels)
    for (auto& coefficient : kernel)
      coefficient = ComplexBall(0);
  ambiguous_multiplier.kernels.front().front() = ComplexBall(1);
  ambiguous_row.entries.push_back(
      {0, 0, std::move(ambiguous_multiplier)});
  const auto ambiguous_direct =
      diffexp2::apply_prepared_scalar_row_window(
          ambiguous_row, ambiguous_input,
          ambiguous_input.epsilon.complete_max, "ambiguous-acb-direct");
  ok = ok && ambiguous_direct.has_value() &&
       ambiguous_direct->sectors.size() == 1 &&
       !ambiguous_direct->sectors.front().coefficients.front().is_zero() &&
       ambiguous_direct->sectors.front().coefficients.front().contains_zero();

  PreparedSparseLocalMultiplierMatrix<Rational> empty_scalar;
  empty_scalar.rows = 1;
  empty_scalar.columns = 2;
  ok = ok && !diffexp2::apply_prepared_scalar_row_window(
      empty_scalar, input, 0).has_value();

  std::cout << (ok ? "PASS" : "FAIL")
            << ": native local rational/SCC and direct scalar-row algebra\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
