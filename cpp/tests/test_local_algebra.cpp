#include "diffexp2/local_algebra.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {

using diffexp2::ExactScalarDescriptor;
using diffexp2::LocalSector;
using diffexp2::LocalSolution;
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

  std::cout << (ok ? "PASS" : "FAIL")
            << ": native local rational/SCC source algebra\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
