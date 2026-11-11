#pragma once
// Irreducible denominator factors of an exact connection, classified by what
// they depend on. A factor of the merge parameter alone is a fixed singular
// point and a factor of epsilon alone a pole in the dimension. A factor that
// mixes both, e.g. 5x-1+2eps, is a singular point that moves with the
// dimension. Feynman integrals have no such singularities (Landau points do
// not depend on d): a mixed factor is an apparent singularity of the master
// basis, a point where the chosen masters become linearly dependent.
#include "diffexp/affine_frobenius.hpp"
#include <map>
#include <set>

namespace diffexp::denominator_factors {
enum class Kind { parameter, epsilon, mixed, constant };
struct Factor {
  Exact polynomial;
  Kind kind;
};
using Powers = std::vector<std::pair<std::size_t, unsigned>>;

// Distinct irreducible factors, collected incrementally. A new polynomial is
// first divided by every known factor; only what is left is factored, so a
// scan over many candidates of one level factors each new factor once.
class Catalogue {
  std::size_t xi_, ei_;
  std::vector<Factor> factors_;
  std::map<std::string, Powers> known_;

  static bool divides(const Exact &factor, Exact &value) {
    auto quotient = value / factor;
    if (!quotient.denominator().is_rational())
      return false;
    // The quotient of two integer polynomials can carry an integer content
    // in its denominator; only the polynomial part matters here.
    value = quotient.numerator();
    return true;
  }
  Kind classify(const Exact &p) const {
    const bool x = !p.derivative(xi_).is_zero(), e = !p.derivative(ei_).is_zero();
    return x && e ? Kind::mixed : x ? Kind::parameter : e ? Kind::epsilon : Kind::constant;
  }

public:
  Catalogue(std::size_t parameter, std::size_t epsilon) : xi_(parameter), ei_(epsilon) {}
  const std::vector<Factor> &factors() const { return factors_; }
  // Factor indices and multiplicities of a polynomial (the numerator of the
  // argument is read; pass value.denominator() for a pole structure).
  const Powers &split(const Exact &polynomial) {
    static const Powers none;
    if (polynomial.is_zero() || polynomial.is_rational())
      return none;
    auto [found, inserted] = known_.try_emplace(polynomial.str());
    if (!inserted)
      return found->second;
    Powers powers;
    auto rest = polynomial.numerator();
    for (std::size_t k = 0; k < factors_.size() && !rest.is_rational(); ++k) {
      unsigned power = 0;
      while (divides(factors_[k].polynomial, rest))
        ++power;
      if (power)
        powers.emplace_back(k, power);
    }
    if (!rest.is_rational())
      for (auto &[factor, power] : affine_frobenius_detail::polynomial_factors(rest)) {
        const auto kind = classify(factor);
        if (kind == Kind::constant)
          continue;
        powers.emplace_back(factors_.size(), static_cast<unsigned>(power));
        factors_.push_back({std::move(factor), kind});
      }
    return found->second = std::move(powers);
  }
  // Order of a factor in a rational function: numerator minus denominator.
  long order(const Exact &value, std::size_t factor) {
    long result = 0;
    for (const auto &[k, power] : split(value.numerator()))
      if (k == factor) result += power;
    for (const auto &[k, power] : split(value.denominator()))
      if (k == factor) result -= power;
    return result;
  }
};

// Sum over the distinct mixed factors of a matrix of the highest power in
// which each occurs in a denominator: the total order of the moving poles.
// Zero means every denominator is a function of the parameter times a
// function of epsilon.
struct MixedWeight {
  std::size_t factors = 0, total_order = 0, entries = 0;
};
inline MixedWeight mixed_weight(const std::vector<std::vector<Exact>> &matrix, Catalogue &catalogue) {
  std::map<std::size_t, unsigned> highest;
  MixedWeight weight;
  for (const auto &row : matrix)
    for (const auto &value : row) {
      if (value.is_zero())
        continue;
      bool mixed = false;
      for (const auto &[k, power] : catalogue.split(value.denominator()))
        if (catalogue.factors()[k].kind == Kind::mixed) {
          mixed = true;
          auto &seen = highest[k];
          seen = std::max(seen, power);
        }
      weight.entries += mixed;
    }
  weight.factors = highest.size();
  for (const auto &[k, power] : highest)
    weight.total_order += power;
  return weight;
}
} // namespace diffexp::denominator_factors
