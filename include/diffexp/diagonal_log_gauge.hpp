#pragma once
#include "diffexp/exact.hpp"
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly_factor.h>

namespace diffexp {
struct DiagonalLogGaugeOptions {
  unsigned max_degree = 32;
  unsigned max_coefficient_bits = 4096;
  unsigned max_factors = 32;
  bool half_integer_only = true;
};
struct DiagonalLogFactor {
  Exact polynomial;
  Rational exponent;
};
struct DiagonalLogGaugeResult {
  bool supported = false;
  std::vector<DiagonalLogFactor> factors;
  std::string rejection_reason;
};
namespace diagonal_log_gauge_detail {
struct Polynomials {
  fmpz_poly_t numerator, denominator, one;
  fmpz_poly_factor_t factors;
  fmpq_poly_t n, d, p, quotient, remainder, derivative, product, rn, rd;
  fmpq_t a, b;
  Polynomials() {
    fmpz_poly_init(numerator);
    fmpz_poly_init(denominator);
    fmpz_poly_init(one);
    fmpz_poly_one(one);
    fmpz_poly_factor_init(factors);
    for (auto v : {n, d, p, quotient, remainder, derivative, product, rn, rd})
      fmpq_poly_init(v);
    fmpq_init(a);
    fmpq_init(b);
  }
  ~Polynomials() {
    fmpq_clear(a);
    fmpq_clear(b);
    for (auto v : {n, d, p, quotient, remainder, derivative, product, rn, rd})
      fmpq_poly_clear(v);
    fmpz_poly_factor_clear(factors);
    fmpz_poly_clear(one);
    fmpz_poly_clear(denominator);
    fmpz_poly_clear(numerator);
  }
  Polynomials(const Polynomials &) = delete;
};
} // namespace diagonal_log_gauge_detail
// q = sum exponent * polynomial'/polynomial, verified over the owned exact
// field. This only identifies a gauge; callers own path branches and boundary
// transforms.
inline DiagonalLogGaugeResult
diagonal_log_gauge(const Exact &q, std::size_t variable,
                   const DiagonalLogGaugeOptions &options = {}) {
  auto reject = [](const char *reason) {
    return DiagonalLogGaugeResult{false, {}, reason};
  };
  if (variable >= q.variable_count())
    return reject("invalid polynomial variable");
  if (!options.max_degree || !options.max_coefficient_bits ||
      !options.max_factors)
    return reject("invalid diagonal gauge budget");
  if (q.is_zero())
    return {true, {}, {}};
  // Check sparse exponents before allocating dense univariate polynomials.
  for (const auto &terms : {q.numerator_terms(), q.denominator_terms()})
    for (const auto &term : terms) {
      for (std::size_t i = 0; i < term.powers.size(); ++i) {
        if (i != variable && term.powers[i])
          return reject("coefficient is not univariate");
        if (term.powers[i] > options.max_degree)
          return reject("polynomial degree budget exceeded");
      }
    }
  diagonal_log_gauge_detail::Polynomials z;
  q.univariate_polynomials(z.numerator, z.denominator, variable);
  for (auto v : {z.numerator, z.denominator})
    if (fmpz_poly_max_bits(v) >
            static_cast<slong>(options.max_coefficient_bits) ||
        fmpz_poly_max_bits(v) <
            -static_cast<slong>(options.max_coefficient_bits))
      return reject("polynomial coefficient bit budget exceeded");
  if (fmpz_poly_degree(z.numerator) >= fmpz_poly_degree(z.denominator))
    return reject("nonzero polynomial part");
  fmpz_poly_factor(z.factors, z.denominator);
  if (z.factors->num > static_cast<slong>(options.max_factors))
    return reject("factor count budget exceeded");
  for (slong i = 0; i < z.factors->num; ++i)
    if (z.factors->exp[i] != 1)
      return reject("higher order pole");
  fmpq_poly_set_fmpz_poly(z.n, z.numerator);
  fmpq_poly_set_fmpz_poly(z.d, z.denominator);
  DiagonalLogGaugeResult result{true, {}, {}};
  auto reconstructed = q.constant(0);
  for (slong i = 0; i < z.factors->num; ++i) {
    auto factor = z.factors->p + i;
    fmpq_poly_set_fmpz_poly(z.p, factor);
    fmpq_poly_divrem(z.quotient, z.remainder, z.d, z.p);
    if (!fmpq_poly_is_zero(z.remainder))
      return reject("factor division failed");
    fmpq_poly_derivative(z.derivative, z.p);
    fmpq_poly_mul(z.product, z.derivative, z.quotient);
    fmpq_poly_rem(z.rn, z.n, z.p);
    fmpq_poly_rem(z.rd, z.product, z.p);
    auto degree = fmpq_poly_degree(z.rd);
    if (degree < 0 || fmpq_poly_degree(z.rn) != degree)
      return reject("residue is not constant on an irreducible factor");
    fmpq_poly_get_coeff_fmpq(z.a, z.rn, degree);
    fmpq_poly_get_coeff_fmpq(z.b, z.rd, degree);
    fmpq_div(z.a, z.a, z.b);
    if (options.half_integer_only && !fmpz_is_one(fmpq_denref(z.a)) &&
        !fmpz_equal_ui(fmpq_denref(z.a), 2))
      return reject("exponent is neither integer nor half integer");
    char *raw = fmpq_get_str(nullptr, 10, z.a);
    if (!raw)
      throw std::bad_alloc();
    std::string text(raw);
    flint_free(raw);
    Rational exponent(text);
    auto polynomial = q.from_univariate_polynomials(factor, z.one, variable);
    reconstructed = reconstructed + q.constant(exponent) *
                                        polynomial.derivative(variable) /
                                        polynomial;
    result.factors.push_back({std::move(polynomial), std::move(exponent)});
  }
  if (reconstructed != q)
    return reject("exact logarithmic derivative reconstruction failed");
  return result;
}
} // namespace diffexp
