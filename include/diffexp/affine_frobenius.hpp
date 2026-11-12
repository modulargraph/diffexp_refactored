#pragma once
#include "diffexp/fuchsify.hpp"
#include "diffexp/jet.hpp"
#include "diffexp/univariate_rational.hpp"
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mpoly_factor.h>
#include <atomic>
#include <exception>
#include <functional>
#include <iterator>
#include <map>
#include <mutex>
#include <thread>
#include <tuple>

namespace diffexp {
namespace affine_frobenius_detail {
using Matrix = fuchsify::Matrix;
using Poly = std::vector<Exact>;
struct Exponent {
  Rational power, slope;
};
// Exponent offsets are tested on the canonical fraction itself; printing a
// rational only to look for '/' allocates in the innermost projection loops.
inline bool is_integer(const Rational &q) {
  return fmpz_is_one(fmpq_denref(q.raw()));
}
inline long integer_value(const Rational &q) {
  if (!is_integer(q) || !fmpz_fits_si(fmpq_numref(q.raw())))
    throw std::domain_error("affine exponent offset is not a machine integer");
  return fmpz_get_si(fmpq_numref(q.raw()));
}
template<class Coefficient> inline void trim(std::vector<Coefficient> &p) {
  while (!p.empty() && p.back().is_zero())
    p.pop_back();
}
template<class Coefficient> inline std::vector<Coefficient> solve(const std::vector<Coefficient> &rhs, const Coefficient &delta) {
  if (rhs.empty())
    return {};
  std::vector<Coefficient> out(rhs.size() + (delta.is_zero() ? 1 : 0), delta.constant(0));
  if (delta.is_zero()) {
    for (unsigned l = 0; l < rhs.size(); ++l)
      out[l + 1] = rhs[l] / delta.constant(l + 1);
  } else
    for (unsigned l = rhs.size(); l-- > 0;)
      out[l] =
          (rhs[l] - (l + 1 < out.size() ? delta.constant(l + 1) * out[l + 1]
                                        : delta.constant(0))) /
          delta;
  trim(out);
  return out;
}
inline std::vector<Exact> taylor(const Exact &e, std::size_t xi, unsigned top) {
  auto collect = [&](const auto &terms) {
    std::vector<Exact> out(top + 1, e.constant(0));
    for (auto &t : terms) {
      auto n = t.powers[xi];
      if (n > top)
        continue;
      auto c = e.constant(t.coefficient);
      for (unsigned j = 0; j < t.powers.size(); ++j)
        if (j != xi && t.powers[j])
          c = c * e.variable(j).pow(t.powers[j]);
      out[n] = out[n] + c;
    }
    return out;
  };
  auto numerator = collect(e.numerator_terms()),
       denominator = collect(e.denominator_terms());
  if (denominator[0].is_zero())
    throw std::domain_error(
        "affine Frobenius coefficient is not x-analytic over Q(eps)");
  std::vector<Exact> out(top + 1, e.constant(0));
  for (unsigned n = 0; n <= top; ++n) {
    auto v = numerator[n];
    for (unsigned m = 1; m <= n; ++m)
      if (!denominator[m].is_zero())
        v = v - denominator[m] * out[n - m];
    out[n] = v / denominator[0];
  }
  return out;
}
inline Exact characteristic(const Matrix &r, std::size_t xi) {
  using namespace fuchsify::detail;
  unsigned d = r.size();
  auto z = r[0][0].constant(0), lambda = z.variable(xi);
  auto b = identity(d, z);
  Exact p = lambda.pow(d);
  for (unsigned k = 1; k <= d; ++k) {
    b = multiply(r, b);
    auto c = z;
    for (unsigned i = 0; i < d; ++i)
      c = c - b[i][i];
    c = c / z.constant(k);
    p = p + c * lambda.pow(d - k);
    for (unsigned i = 0; i < d; ++i)
      b[i][i] = b[i][i] + c;
  }
  return p;
}
inline std::vector<std::pair<Exact,long>> polynomial_factors(const Exact& p) {
  const auto z=p.constant(0);
  auto names = p.variables();
  std::vector<const char *> symbols;
  for (auto &n : names)
    symbols.push_back(n.c_str());
  // The factor terms are rebuilt with throwing Exact arithmetic, so the FLINT
  // workspace is owned rather than cleared by hand on the success path only.
  struct Workspace {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t raw;
    fmpz_mpoly_factor_t factors;
    fmpz_t coefficient;
    explicit Workspace(slong variables) {
      fmpz_mpoly_ctx_init(ctx, variables, ORD_DEGLEX);
      fmpz_mpoly_init(raw, ctx);
      fmpz_mpoly_factor_init(factors, ctx);
      fmpz_init(coefficient);
    }
    ~Workspace() {
      fmpz_clear(coefficient);
      fmpz_mpoly_factor_clear(factors, ctx);
      fmpz_mpoly_clear(raw, ctx);
      fmpz_mpoly_ctx_clear(ctx);
    }
    Workspace(const Workspace &) = delete;
    Workspace &operator=(const Workspace &) = delete;
  } w(static_cast<slong>(names.size()));
  if (fmpz_mpoly_set_str_pretty(w.raw, p.numerator().str().c_str(),
                                symbols.data(), w.ctx) != 0 ||
      fmpz_mpoly_factor(w.factors, w.raw, w.ctx) == 0)
    throw std::runtime_error("exact affine polynomial factorization failed");
  std::vector<std::pair<Exact, long>> polynomials;
  for (slong f = 0; f < w.factors->num; ++f) {
    auto poly = z;
    for (slong t = 0; t < fmpz_mpoly_length(w.factors->poly + f, w.ctx); ++t) {
      fmpz_mpoly_get_term_coeff_fmpz(w.coefficient, w.factors->poly + f, t,
                                     w.ctx);
      auto c = z.constant(Rational::from_integer(w.coefficient));
      std::vector<unsigned long> powers(names.size());
      fmpz_mpoly_get_term_exp_ui(powers.data(), w.factors->poly + f, t, w.ctx);
      for (unsigned j = 0; j < powers.size(); ++j)
        if (powers[j])
          c = c * z.variable(j).pow(powers[j]);
      poly = poly + c;
    }
    polynomials.emplace_back(
        poly, fmpz_mpoly_factor_get_exp_si(w.factors, f, w.ctx));
  }
  return polynomials;
}
// Factor the exact characteristic polynomial; every spectral factor must be
// linear and its root exactly affine in epsilon. No sampled eigenvalues.
inline std::vector<Exponent> spectrum(const Matrix &r, std::size_t xi,
                                      std::size_t ei) {
  const auto p=characteristic(r,xi),z=p.constant(0);
  auto polynomials=polynomial_factors(p);
  std::vector<Exponent> out;
  Exact product = z.constant(1);
  unsigned count = 0;
  for (auto &[factor, multiplicity] : polynomials) {
    auto linear = factor.derivative(xi);
    if (linear.is_zero())
      continue;
    if (!linear.derivative(xi).is_zero())
      throw std::domain_error(
          "unsupported residue spectrum: nonlinear spectral factor");
    auto root = -(factor - linear * z.variable(xi)) / linear;
    if (!root.derivative(ei).derivative(ei).is_zero())
      throw std::domain_error(
          "unsupported residue spectrum: eigenvalue is nonlinear in epsilon");
    auto slope = root.derivative(ei), power = root - slope * z.variable(ei);
    if (!power.is_rational() || !slope.is_rational())
      throw std::domain_error(
          "unsupported residue spectrum: nonrational affine coefficients");
    Exponent exponent{power.rational(), slope.rational()};
    out.push_back(exponent);
    for (long k = 0; k < multiplicity; ++k) {
      product = product * (z.variable(xi) - root);
      ++count;
    }
  }
  if (count != r.size() || !(product == p))
    throw std::logic_error("affine characteristic factor certificate failed");
  std::sort(out.begin(), out.end(), [](const auto &a, const auto &b) {
    return std::tie(a.power, a.slope) < std::tie(b.power, b.slope);
  });
  return out;
}
struct Frame {
  Matrix transform;
  std::vector<Exponent> exponents;
  std::vector<int> successor;
};
inline Frame frame(const Matrix &r, std::size_t xi, std::size_t ei) {
  using namespace fuchsify::detail;
  unsigned d = r.size();
  auto z = r[0][0].constant(0);
  Frame out{zeros(d, d, z), {}, {}};
  Matrix selected;
  auto jordan = zeros(d, d, z);
  for (auto &root : spectrum(r, xi, ei)) {
    auto lambda =
        z.constant(root.power) + z.constant(root.slope) * z.variable(ei);
    auto m = r;
    for (unsigned i = 0; i < d; ++i)
      m[i][i] = m[i][i] - lambda;
    auto power = identity(d, z);
    std::vector<Matrix> kernels(1);
    for (unsigned l = 1; l <= d; ++l) {
      power = multiply(power, m);
      kernels.push_back(nullspace(power));
      // Kernels form an increasing chain. Equal consecutive dimensions
      // prove stabilization; a full kernel also ends every Jordan ladder.
      if (kernels[l].size() == kernels[l - 1].size()) {
        kernels.pop_back();
        break;
      }
      if (kernels[l].size() == d)
        break;
    }
    for (unsigned l = kernels.size() - 1; l > 0; --l) {
      Matrix span;
      for (auto &v : selected)
        extend(span, v);
      for (auto &v : kernels[l - 1])
        extend(span, v);
      for (auto &head : kernels[l])
        if (extend(span, head)) {
          std::vector<std::vector<Exact>> chain(l);
          chain.back() = head;
          for (unsigned j = l - 1; j > 0; --j)
            chain[j - 1] = action(m, chain[j]);
          for (unsigned j = 0; j < l; ++j) {
            unsigned column = out.exponents.size();
            if (column >= d)
              throw std::logic_error("affine Jordan frame overflow");
            for (unsigned i = 0; i < d; ++i)
              out.transform[i][column] = chain[j][i];
            jordan[column][column] = lambda;
            if (j + 1 < l)
              jordan[column][column + 1] = z.constant(1);
            out.exponents.push_back(root);
            out.successor.push_back(j + 1 < l ? static_cast<int>(column + 1)
                                              : -1);
            selected.push_back(chain[j]);
            extend(span, chain[j]);
          }
        }
    }
  }
  if (out.exponents.size() != d ||
      multiply(r, out.transform) != multiply(out.transform, jordan))
    throw std::logic_error("affine Jordan residue certificate failed");
  (void)inverse(out.transform);
  return out;
}
inline Exact determinant(Matrix a) {
  auto value = a[0][0].constant(1);
  const unsigned d = a.size();
  for (unsigned k = 0; k < d; ++k) {
    unsigned pivot = k;
    while (pivot < d && a[pivot][k].is_zero())
      ++pivot;
    if (pivot == d)
      return value.constant(0);
    if (pivot != k) {
      std::swap(a[pivot], a[k]);
      value = -value;
    }
    const auto diagonal = a[k][k];
    value = value * diagonal;
    for (unsigned i = k + 1; i < d; ++i)
      if (!a[i][k].is_zero()) {
        auto factor = a[i][k] / diagonal;
        for (unsigned j = k + 1; j < d; ++j)
          a[i][j] = a[i][j] - factor * a[k][j];
        a[i][k] = value.constant(0);
      }
  }
  return value;
}
inline kernel::ComplexBall ball(const Rational &q) {
  kernel::ComplexBall b;
  acb_set_fmpq(b.raw(), q.raw(), kernel::ComplexBall::precision());
  return b;
}
// Read the already exact univariate coefficient directly from FLINT. Building
// a nested expression tree from its printed polynomial makes long endpoint
// coefficients quadratic to parse and copies large subtrees repeatedly.
class EpsilonRational {
  struct Polynomial {
    fmpz_poly_t value;
    Polynomial(){fmpz_poly_init(value);}
    ~Polynomial(){fmpz_poly_clear(value);}
    Polynomial(const Polynomial&)=delete;
    Polynomial& operator=(const Polynomial&)=delete;
  } numerator_,denominator_;
  slong numerator_low_=0,denominator_low_=0;
  static slong lowest(const fmpz_poly_t p){
    slong k=0;while(k<fmpz_poly_length(p)&&fmpz_is_zero(fmpz_poly_get_coeff_ptr(p,k)))++k;return k;
  }
 public:
  EpsilonRational(const Exact& value,std::size_t variable){
    if(value.is_zero())throw std::invalid_argument("zero epsilon coefficient has no valuation");
    value.univariate_polynomials(numerator_.value,denominator_.value,variable);
    numerator_low_=lowest(numerator_.value);denominator_low_=lowest(denominator_.value);
  }
  long valuation()const{return numerator_low_-denominator_low_;}
  // Exact Taylor coefficients of value/eps^valuation() through eps^(length-1):
  // one power-series division in Q[[eps]].
  std::vector<Rational> regular_taylor(unsigned length)const{
    struct Series{fmpq_poly_t value;Series(){fmpq_poly_init(value);}~Series(){fmpq_poly_clear(value);}
      Series(const Series&)=delete;Series& operator=(const Series&)=delete;}n,d,quotient;
    // Only the first `length` coefficients above the valuation can matter;
    // the full numerator of a late x order is far longer than that.
    const auto leading=[&](fmpq_poly_t out,const fmpz_poly_t p,slong low){
      for(slong k=std::min<slong>(fmpz_poly_length(p)-low,length);k-->0;)fmpq_poly_set_coeff_fmpz(out,k,fmpz_poly_get_coeff_ptr(p,low+k));
    };
    leading(n.value,numerator_.value,numerator_low_);leading(d.value,denominator_.value,denominator_low_);
    fmpq_poly_div_series(quotient.value,n.value,d.value,length);
    fmpq_t coefficient;fmpq_init(coefficient);std::vector<Rational> out;out.reserve(length);
    for(unsigned k=0;k<length;++k){fmpq_poly_get_coeff_fmpq(coefficient,quotient.value,k);out.push_back(Rational::from_fmpq(coefficient));}
    fmpq_clear(coefficient);return out;
  }
  Jet regular_jet(unsigned length)const{
    using B=kernel::ComplexBall;
    auto polynomial=[&](const fmpz_poly_t p,slong low){
      Jet out(0,length,B::precision());
      for(slong k=low;k<fmpz_poly_length(p)&&k-low<length;++k){
        B coefficient;acb_set_fmpz(coefficient.raw(),fmpz_poly_get_coeff_ptr(p,k));out.set(k-low,coefficient);
      }return out;
    };
    return polynomial(numerator_.value,numerator_low_)/polynomial(denominator_.value,denominator_low_);
  }
};
} // namespace affine_frobenius_detail
// A reader asked a truncated endpoint series for epsilon orders it does not
// hold; `missing` more orders would have answered it. Never a failure of the
// method: the owner of the series builds it again with more orders.
struct EpsilonOrdersNeeded : std::runtime_error {
  long missing;
  explicit EpsilonOrdersNeeded(long orders)
      : std::runtime_error("truncated endpoint series: " + std::to_string(orders) +
                           " more epsilon orders are needed"),
        missing(orders) {}
};
namespace affine_frobenius_detail {
// A Laurent series in epsilon with exact rational coefficients, known modulo
// eps^order(): eps^shift*u(eps) + O(eps^order). An order of exact_order is a
// Laurent polynomial known exactly. Sums and products carry the order the way
// truncated power series do; nothing is ever read beyond it, and a value whose
// known part vanishes is "zero to its order", never the exact zero.
class EpsilonSeries {
public:
  static constexpr long exact_order = std::numeric_limits<long>::max() / 4;
  EpsilonSeries() { fmpq_poly_init(u_); }
  explicit EpsilonSeries(long n) : EpsilonSeries() { fmpq_poly_set_si(u_, n); }
  EpsilonSeries(const EpsilonSeries &o) : EpsilonSeries() {
    fmpq_poly_set(u_, o.u_);
    shift_ = o.shift_;
    order_ = o.order_;
  }
  EpsilonSeries(EpsilonSeries &&o) noexcept : EpsilonSeries() { swap(o); }
  EpsilonSeries &operator=(EpsilonSeries o) noexcept {
    swap(o);
    return *this;
  }
  ~EpsilonSeries() { fmpq_poly_clear(u_); }
  // The rational function `value` of one variable, modulo eps^order. With
  // exact_order the denominator must be a monomial.
  EpsilonSeries(const Exact &value, std::size_t variable, long order)
      : EpsilonSeries() {
    order_ = order;
    if (value.is_zero())
      return;
    Polynomial numerator, denominator;
    value.univariate_polynomials(numerator.value, denominator.value, variable);
    const slong low = lowest(numerator.value), pole = lowest(denominator.value);
    fmpz_poly_shift_right(numerator.value, numerator.value, low);
    fmpz_poly_shift_right(denominator.value, denominator.value, pole);
    shift_ = low - pole;
    fmpq_poly_set_fmpz_poly(u_, numerator.value);
    if (fmpz_poly_length(denominator.value) == 1)
      fmpq_poly_scalar_div_fmpz(u_, u_, denominator.value->coeffs);
    else {
      if (order >= exact_order)
        throw std::logic_error("an exact epsilon series needs a monomial denominator");
      const long length = order - shift_;
      if (length > 0) {
        fmpq_poly_t d;
        fmpq_poly_init(d);
        fmpq_poly_set_fmpz_poly(d, denominator.value);
        fmpq_poly_div_series(u_, u_, d, length);
        fmpq_poly_clear(d);
      }
    }
    normalize();
  }
  // `terms` coefficients from the lowest one on; exact if that is all of it.
  static EpsilonSeries expand(const Exact &value, std::size_t variable, long terms) {
    if (value.is_zero())
      return {};
    if (value.denominator_terms().size() == 1)
      return EpsilonSeries(value, variable, exact_order);
    return EpsilonSeries(value, variable, fuchsify::detail::valuation(value, variable) + terms);
  }
  // eps^shift*numerator(eps)/denominator + O(eps^order).
  static EpsilonSeries from_integers(const fmpz_poly_t numerator, const fmpz_t denominator,
                                     long shift, long order) {
    EpsilonSeries out;
    fmpq_poly_set_fmpz_poly(out.u_, numerator);
    fmpq_poly_scalar_div_fmpz(out.u_, out.u_, denominator);
    out.shift_ = shift;
    out.order_ = order;
    out.normalize();
    return out;
  }
  static EpsilonSeries linear(const Rational &constant, const Rational &slope) {
    EpsilonSeries out;
    fmpq_poly_set_coeff_fmpq(out.u_, 0, constant.raw());
    fmpq_poly_set_coeff_fmpq(out.u_, 1, slope.raw());
    out.normalize();
    return out;
  }
  EpsilonSeries constant(long n) const { return EpsilonSeries(n); }
  // Zero to the known order (which includes the exact zero).
  bool is_zero() const { return fmpq_poly_is_zero(u_); }
  bool is_exact() const { return order_ >= exact_order; }
  long order() const { return order_; }
  // Lowest known exponent; the order itself for a value that is zero to it.
  long valuation() const { return is_zero() ? order_ : shift_; }
  // Forget the orders from `order` on.
  void truncate(long order) {
    if (order < order_) {
      order_ = order;
      normalize();
    }
  }
  // A series in t = eps/scale, rewritten as the series in eps.
  void divide_variable(const fmpz_t scale) {
    if (is_zero() || fmpz_is_one(scale))
      return;
    fmpq_t a;
    fmpq_init(a);
    fmpq_set_fmpz(a, scale);
    fmpq_inv(a, a);
    fmpq_poly_rescale(u_, u_, a);
    fmpq_pow_si(a, a, shift_);
    fmpq_poly_scalar_mul_fmpq(u_, u_, a);
    fmpq_clear(a);
  }
  // The series divided by a power series with a nonzero constant term.
  void divide_unit(const fmpz_poly_t unit) {
    if (is_zero())
      return;
    if (fmpz_poly_length(unit) == 1) {
      fmpq_poly_scalar_div_fmpz(u_, u_, unit->coeffs);
      return;
    }
    if (is_exact())
      throw std::logic_error("an exact epsilon series needs a monomial denominator");
    fmpq_poly_t d;
    fmpq_poly_init(d);
    fmpq_poly_set_fmpz_poly(d, unit);
    fmpq_poly_div_series(u_, u_, d, order_ - shift_);
    fmpq_poly_clear(d);
    normalize();
  }
  void shift(long powers) {
    if (!is_zero())
      shift_ += powers;
    if (!is_exact())
      order_ += powers;
  }
  Rational coefficient(long exponent) const {
    fmpq_t c;
    fmpq_init(c);
    if (!is_zero() && exponent >= shift_ && exponent - shift_ < fmpq_poly_length(u_))
      fmpq_poly_get_coeff_fmpq(c, u_, exponent - shift_);
    auto out = Rational::from_fmpq(c);
    fmpq_clear(c);
    return out;
  }
  // The known part as a Laurent polynomial of the field of `sample`.
  Exact exact(const Exact &sample, std::size_t variable) const {
    if (is_zero())
      return sample.constant(0);
    Polynomial numerator, denominator;
    fmpq_poly_get_numerator(numerator.value, u_);
    fmpz_poly_set_fmpz(denominator.value, fmpq_poly_denref(u_));
    if (shift_ > 0)
      fmpz_poly_shift_left(numerator.value, numerator.value, shift_);
    else if (shift_ < 0)
      fmpz_poly_shift_left(denominator.value, denominator.value, -shift_);
    return sample.from_reduced_univariate_polynomials(numerator.value, denominator.value, variable);
  }
  EpsilonSeries &operator+=(const EpsilonSeries &b) { return add(b, false); }
  EpsilonSeries &operator-=(const EpsilonSeries &b) { return add(b, true); }
  friend EpsilonSeries operator+(EpsilonSeries a, const EpsilonSeries &b) { a += b; return a; }
  friend EpsilonSeries operator-(EpsilonSeries a, const EpsilonSeries &b) { a -= b; return a; }
  friend EpsilonSeries operator-(EpsilonSeries a) { fmpq_poly_neg(a.u_, a.u_); return a; }
  friend EpsilonSeries operator*(const EpsilonSeries &a, const EpsilonSeries &b) {
    EpsilonSeries out;
    out.order_ = std::min(sum(a.valuation(), b.order_), sum(b.valuation(), a.order_));
    if (a.is_zero() || b.is_zero())
      return out;
    out.shift_ = a.shift_ + b.shift_;
    if (out.is_exact())
      fmpq_poly_mul(out.u_, a.u_, b.u_);
    else if (out.order_ > out.shift_)
      fmpq_poly_mullow(out.u_, a.u_, b.u_, out.order_ - out.shift_);
    out.normalize();
    return out;
  }
  friend EpsilonSeries operator/(const EpsilonSeries &a, const EpsilonSeries &b) {
    if (b.is_zero())
      throw std::domain_error("division by an epsilon series that is zero to its known order");
    EpsilonSeries out;
    // 1/b is known through b.order - 2*valuation(b).
    out.order_ = std::min(sum(a.order_, -b.shift_),
                          sum(a.valuation(), sum(b.order_, -2 * b.shift_)));
    if (a.is_zero())
      return out;
    out.shift_ = a.shift_ - b.shift_;
    if (fmpq_poly_length(b.u_) == 1) {
      fmpq_t c;
      fmpq_init(c);
      fmpq_poly_get_coeff_fmpq(c, b.u_, 0);
      fmpq_poly_scalar_div_fmpq(out.u_, a.u_, c);
      fmpq_clear(c);
    } else if (out.is_exact())
      throw std::logic_error("an exact epsilon series division needs a monomial divisor");
    else if (out.order_ > out.shift_)
      fmpq_poly_div_series(out.u_, a.u_, b.u_, out.order_ - out.shift_);
    out.normalize();
    return out;
  }

private:
  struct Polynomial {
    fmpz_poly_t value;
    Polynomial() { fmpz_poly_init(value); }
    ~Polynomial() { fmpz_poly_clear(value); }
    Polynomial(const Polynomial &) = delete;
    Polynomial &operator=(const Polynomial &) = delete;
  };
  fmpq_poly_t u_;
  long shift_ = 0, order_ = exact_order;
  static slong lowest(const fmpz_poly_t p) {
    slong k = 0;
    while (k < fmpz_poly_length(p) && fmpz_is_zero(fmpz_poly_get_coeff_ptr(p, k)))
      ++k;
    return k;
  }
  static long sum(long a, long b) {
    return a >= exact_order || b >= exact_order ? exact_order : a + b;
  }
  void swap(EpsilonSeries &o) noexcept {
    fmpq_poly_swap(u_, o.u_);
    std::swap(shift_, o.shift_);
    std::swap(order_, o.order_);
  }
  // No coefficient at or above the order, a nonzero lowest coefficient.
  void normalize() {
    if (!is_exact()) {
      if (order_ - shift_ <= 0)
        fmpq_poly_zero(u_);
      else
        fmpq_poly_truncate(u_, order_ - shift_);
    }
    slong low = 0;
    const slong length = fmpq_poly_length(u_);
    while (low < length && fmpz_is_zero(fmpq_poly_numref(u_) + low))
      ++low;
    if (low == length) {
      fmpq_poly_zero(u_);
      shift_ = 0;
    } else if (low) {
      fmpq_poly_shift_right(u_, u_, low);
      shift_ += low;
    }
  }
  EpsilonSeries &add(const EpsilonSeries &b, bool subtract) {
    const long order = std::min(order_, b.order_);
    if (!b.is_zero()) {
      if (is_zero()) {
        fmpq_poly_set(u_, b.u_);
        if (subtract)
          fmpq_poly_neg(u_, u_);
        shift_ = b.shift_;
      } else {
        const long low = std::min(shift_, b.shift_);
        fmpq_poly_t other;
        fmpq_poly_init(other);
        fmpq_poly_shift_left(other, b.u_, b.shift_ - low);
        fmpq_poly_shift_left(u_, u_, shift_ - low);
        if (subtract)
          fmpq_poly_sub(u_, u_, other);
        else
          fmpq_poly_add(u_, u_, other);
        fmpq_poly_clear(other);
        shift_ = low;
      }
    }
    order_ = order;
    normalize();
    return *this;
  }
};
struct Integer {
  fmpz_t value;
  Integer() { fmpz_init(value); }
  explicit Integer(long n) { fmpz_init_set_si(value, n); }
  Integer(const Integer &other) { fmpz_init_set(value, other.value); }
  Integer &operator=(const Integer &other) {
    fmpz_set(value, other.value);
    return *this;
  }
  ~Integer() { fmpz_clear(value); }
};
// rho*(m + s*eps) in Z[eps] for the least positive integer rho.
inline IntegerPolynomial linear(const Rational &m, const Rational &s,
                                Integer &rho) {
  fmpz_lcm(rho.value, fmpq_denref(m.raw()), fmpq_denref(s.raw()));
  IntegerPolynomial out;
  Integer c;
  fmpz_divexact(c.value, rho.value, fmpq_denref(s.raw()));
  fmpz_mul(c.value, c.value, fmpq_numref(s.raw()));
  fmpz_poly_set_coeff_fmpz(out.raw(), 1, c.value);
  fmpz_divexact(c.value, rho.value, fmpq_denref(m.raw()));
  fmpz_mul(c.value, c.value, fmpq_numref(m.raw()));
  fmpz_poly_set_coeff_fmpz(out.raw(), 0, c.value);
  return out;
}
// Fraction-free Frobenius recurrence in the constant Jordan frame.
//
// All coefficients of one x order are numerators in Z[eps] over ONE shared
// denominator eps^pole*common; common is the product of the per-order steps
// (the resonance denominators delta and the row scales). Earlier orders are
// brought to the current denominator whenever it grows, so the dominant work,
// a connection coefficient times an earlier numerator, is a small-by-large
// integer polynomial product with no gcd. A fraction is reduced against the
// few small factors its own row introduces, and each emitted physical
// coefficient is put in lowest terms once by the caller.
class IntegerRecurrence {
public:
  using Z = IntegerPolynomial;
  using Value = std::vector<Z>; // coefficients of log(x)^l, no trailing zero
  IntegerRecurrence(const std::vector<Matrix> &g,
                    const std::vector<std::vector<Exact>> &derivative_lags,
                    const Matrix &frame, std::size_t ei) {
    const unsigned d = frame.size(), n = g.size() - 1;
    const auto one = frame[0][0].constant(1);
    const auto widen = [](Exact &scale, const Exact &value) {
      const auto denominator = value.denominator();
      if (!denominator.is_one())
        scale = scale.polynomial_lcm(denominator);
    };
    rows_.resize(d);
    frame_.resize(d);
    for (unsigned row = 0; row < d; ++row) {
      const auto &q = derivative_lags[row];
      auto scale = one;
      unsigned last = 0;
      for (unsigned lag = 1; lag <= n; ++lag) {
        for (unsigned j = 0; j < d; ++j)
          if (!g[lag][row][j].is_zero()) {
            widen(scale, g[lag][row][j]);
            last = lag;
          }
        if (lag < q.size() && !q[lag].is_zero()) {
          widen(scale, q[lag]);
          last = lag;
        }
      }
      auto &r = rows_[row];
      r.scale = integral(scale, ei);
      r.lags.resize(last + 1);
      if (!q.empty())
        r.derivative.resize(last + 1);
      for (unsigned lag = 1; lag <= last; ++lag) {
        for (unsigned j = 0; j < d; ++j)
          if (!g[lag][row][j].is_zero())
            r.lags[lag].emplace_back(j, integral(scale * g[lag][row][j], ei));
        if (lag < q.size() && !q[lag].is_zero())
          r.derivative[lag] = integral(scale * q[lag], ei);
      }
      window_ = std::max(window_, last);
      scale = one;
      for (unsigned j = 0; j < d; ++j)
        widen(scale, frame[row][j]);
      frame_[row].scale = integral(scale, ei);
      for (unsigned j = 0; j < d; ++j)
        if (!frame[row][j].is_zero())
          frame_[row].entries.emplace_back(j,
                                           integral(scale * frame[row][j], ei));
    }
  }
  // Solves one fundamental column through x order `top`. For every order and
  // physical row with a nonzero coefficient, emit(k, row, numerators, scale,
  // common, pole) receives the log coefficients numerators[l] over the
  // denominator scale*common*eps^pole (pole may be negative).
  template <class Emit>
  void column(const std::vector<Exponent> &eigen,
              const std::vector<int> &successor, unsigned column, unsigned top,
              std::size_t ladder_bound, Emit &&emit) const {
    const unsigned d = rows_.size();
    Integer rho, shift;
    // rho*alpha at order m is exponent + rho*m.
    const Z exponent = linear(eigen[column].power, eigen[column].slope, rho);
    const bool unit = fmpz_is_one(rho.value);
    std::vector<std::vector<Value>> orders(top + 1, std::vector<Value>(d));
    Z common(1), step, cofactor, weight;
    long pole = 0;
    for (unsigned k = 0; k <= top; ++k) {
      std::vector<Fraction> solved(d);
      for (unsigned row = d; row-- > 0;) {
        const Row &r = rows_[row];
        const bool cleared = !r.derivative.empty();
        // The row equation is multiplied by scale (and by rho when its
        // derivative lags bring the column exponent into the coefficients).
        Fraction rhs;
        rhs.denominator = r.scale;
        if (cleared && !unit)
          fmpz_poly_scalar_mul_fmpz(rhs.denominator.raw(),
                                    rhs.denominator.raw(), rho.value);
        const unsigned last =
            std::min<std::size_t>(k, r.lags.size() - 1);
        for (unsigned lag = 1; lag <= last; ++lag) {
          const auto &source = orders[k - lag];
          for (const auto &[j, coefficient] : r.lags[lag]) {
            if (source[j].empty())
              continue;
            if (cleared && !unit) {
              fmpz_poly_scalar_mul_fmpz(weight.raw(), coefficient.raw(),
                                        rho.value);
              accumulate(rhs.numerator, weight, source[j], 0);
            } else
              accumulate(rhs.numerator, coefficient, source[j], 0);
          }
          const auto &previous = source[row];
          if (!cleared || r.derivative[lag].zero() || previous.empty())
            continue;
          // -q~[lag]*((alpha+k-lag)*v + dv/dL) of the earlier coefficient.
          fmpz_poly_get_coeff_fmpz(shift.value, exponent.raw(), 0);
          fmpz_addmul_ui(shift.value, rho.value, k - lag);
          fmpz_poly_set(weight.raw(), exponent.raw());
          fmpz_poly_set_coeff_fmpz(weight.raw(), 0, shift.value);
          fmpz_poly_mul(weight.raw(), weight.raw(), r.derivative[lag].raw());
          fmpz_poly_neg(weight.raw(), weight.raw());
          accumulate(rhs.numerator, weight, previous, 0);
          fmpz_poly_scalar_mul_fmpz(weight.raw(), r.derivative[lag].raw(),
                                    rho.value);
          fmpz_poly_neg(weight.raw(), weight.raw());
          accumulate(rhs.numerator, weight, previous, 1);
        }
        trim(rhs.numerator);
        if (successor[row] >= 0 && !solved[successor[row]].numerator.empty())
          rhs = sum(rhs, solved[successor[row]]);
        solve(rhs,
              eigen[column].power + Rational(k) - eigen[row].power,
              eigen[column].slope - eigen[row].slope);
        if (k == 0 && row == column) {
          // Rows above the seed carry no source at order zero.
          if (!rhs.numerator.empty())
            throw std::logic_error("seeded Frobenius coefficient has a source");
          rhs = Fraction{};
          rhs.numerator.emplace_back(1);
        }
        if (rhs.numerator.size() > ladder_bound)
          throw std::logic_error("affine Frobenius log ladder bound exceeded");
        reduce(rhs);
        solved[row] = std::move(rhs);
      }
      // One denominator for the whole order: the least common multiple of
      // the row denominators and the deepest epsilon pole.
      long rise = 0;
      fmpz_poly_one(step.raw());
      for (const auto &value : solved)
        if (!value.numerator.empty()) {
          rise = std::max(rise, value.pole);
          fmpz_poly_lcm(step.raw(), step.raw(), value.denominator.raw());
        }
      for (unsigned row = 0; row < d; ++row) {
        auto &value = solved[row];
        if (value.numerator.empty())
          continue;
        fmpz_poly_div(cofactor.raw(), step.raw(), value.denominator.raw());
        for (auto &c : value.numerator)
          scale(c, cofactor, rise - value.pole);
        orders[k][row] = std::move(value.numerator);
      }
      if (rise || !step.one()) {
        for (unsigned old = k + 1 > window_ ? k + 1 - window_ : 0; old < k; ++old)
          for (auto &value : orders[old])
            for (auto &c : value)
              scale(c, step, rise);
        fmpz_poly_mul(common.raw(), common.raw(), step.raw());
        pole += rise;
      }
      Value physical;
      for (unsigned row = 0; row < d; ++row) {
        physical.clear();
        for (const auto &[j, coefficient] : frame_[row].entries)
          if (!orders[k][j].empty())
            accumulate(physical, coefficient, orders[k][j], 0);
        trim(physical);
        if (!physical.empty())
          emit(k, row, physical, frame_[row].scale, common, pole);
      }
      // No later order reaches back further than the longest lag.
      if (k >= window_)
        orders[k - window_] = std::vector<Value>(d);
    }
  }

  // ---- The same recurrence in truncated Laurent arithmetic ----
  //
  // Every coefficient is carried modulo a power of epsilon instead of as a
  // rational function. In Q[[eps]] the resonance denominators m+s*eps (m
  // nonzero) and the row scales are units, so nothing but an integer and a
  // power of epsilon is left in a denominator: all coefficients of one x order
  // are numerators in Z[t] over ONE shared integer and t^pole. The variable is
  // t = eps/L, with L the least common multiple of the integers whose inverse
  // series occur, so that these series have integer coefficients and the
  // t^j coefficient of a numerator grows like L^j instead of every coefficient
  // carrying L^(length).
  //
  // What is known is tracked, not assumed. A value carries the absolute
  // epsilon order through which it is known and a structural bound of its
  // logarithmic degree (from which sources reach it, whatever cancels):
  // dividing by s*eps at a resonance of different slopes costs one order per
  // logarithmic power, a row scale that vanishes at eps=0 its valuation. This
  // bookkeeping does not depend on the numbers, so a pass without arithmetic
  // (structure_only) tells how many orders a column loses, and the arithmetic
  // pass starts with exactly that many more than it must deliver.
  struct Meta {
    long known = EpsilonSeries::exact_order; // absolute epsilon order known through
    int bound = -1;                          // log-degree bound; negative: exactly zero
    void merge(const Meta &source) {
      known = std::min(known, source.known);
      bound = std::max(bound, source.bound);
    }
  };
  // t^valuation times a power series with a nonzero constant term.
  struct Unit {
    Z unit{1};
    slong valuation = 0;
  };
  void enable_truncation(const std::vector<Exponent> &eigen, unsigned top) {
    const unsigned d = rows_.size();
    auto t = std::make_unique<Truncation>();
    fmpz_one(t->scale.value);
    Integer part;
    // Every resonance denominator m+s*eps with m nonzero: (s/m)*L is integral.
    for (unsigned column = 0; column < d; ++column)
      for (unsigned row = 0; row < d; ++row) {
        const auto slope = eigen[column].slope - eigen[row].slope;
        if (slope.is_zero())
          continue;
        for (unsigned k = 0; k <= top; ++k) {
          const auto m = eigen[column].power + Rational(k) - eigen[row].power;
          if (m.is_zero())
            continue;
          const auto ratio = slope / m;
          fmpz_lcm(t->scale.value, t->scale.value, fmpq_denref(ratio.raw()));
        }
      }
    // The constant term of a scale, where that is a short integer.
    const auto constant = [&](const Z &scale) {
      const slong low = scale.valuation();
      fmpz_abs(part.value, scale.raw()->coeffs + low);
      if (fmpz_poly_length(scale.raw()) > low + 1 && fmpz_bits(part.value) <= 64)
        fmpz_lcm(t->scale.value, t->scale.value, part.value);
    };
    for (const auto &row : rows_)
      constant(row.scale);
    for (const auto &row : frame_)
      constant(row.scale);
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_fmpz(q, t->scale.value);
    t->scale_rational = Rational::from_fmpq(q);
    fmpq_clear(q);
    // p(eps) = p(L*t).
    const auto substitute = [&](Z value) {
      fmpz_one(part.value);
      for (slong e = 1; e < fmpz_poly_length(value.raw()); ++e) {
        fmpz_mul(part.value, part.value, t->scale.value);
        fmpz_mul(value.raw()->coeffs + e, value.raw()->coeffs + e, part.value);
      }
      return value;
    };
    const auto unit = [&](const Z &scale) {
      Unit out;
      out.unit = substitute(scale);
      out.valuation = out.unit.valuation();
      fmpz_poly_shift_right(out.unit.raw(), out.unit.raw(), out.valuation);
      return out;
    };
    t->rows = rows_;
    t->frame = frame_;
    for (auto &row : t->rows) {
      t->row_units.push_back(unit(row.scale));
      for (auto &lag : row.lags)
        for (auto &[j, coefficient] : lag)
          coefficient = substitute(coefficient);
      for (auto &q : row.derivative)
        q = substitute(q);
    }
    for (auto &row : t->frame) {
      t->frame_units.push_back(unit(row.scale));
      for (auto &[j, coefficient] : row.entries)
        coefficient = substitute(coefficient);
    }
    truncation_ = std::move(t);
  }
  // eps = scale()*t.
  const fmpz *scale() const { return truncation_->scale.value; }
  // One fundamental column through x order `top`, every value cut at the
  // absolute epsilon order `cap`. For every order and physical row that is not
  // structurally zero, emit(k, row, numerators, unit, denominator, pole, meta)
  // receives the log coefficients numerators[l](t) of
  //   numerators/(denominator * t^(pole+unit.valuation) * unit.unit(t)),
  // valid below t^(meta.known+pole); the physical coefficient is known through
  // the epsilon order meta.known-unit.valuation. With structure_only no number
  // is computed and the numerators are empty.
  template <class Emit>
  void truncated_column(const std::vector<Exponent> &eigen,
                        const std::vector<int> &successor, unsigned column,
                        unsigned top, std::size_t ladder_bound, long cap,
                        bool structure_only, Emit &&emit) const {
    const auto &T = *truncation_;
    const unsigned d = T.rows.size();
    Integer rho, shift;
    const Z exponent = linear(eigen[column].power,
                              eigen[column].slope * T.scale_rational, rho);
    const bool unit = fmpz_is_one(rho.value);
    std::vector<std::vector<Value>> orders(top + 1, std::vector<Value>(d));
    std::vector<std::vector<Meta>> meta(top + 1, std::vector<Meta>(d));
    std::vector<Inverse> inverses(d);
    Integer common(1), step, cofactor;
    Z weight;
    long pole = 0;
    for (unsigned k = 0; k <= top; ++k) {
      std::vector<Solved> solved(d);
      for (unsigned row = d; row-- > 0;) {
        const Row &r = T.rows[row];
        const bool cleared = !r.derivative.empty();
        Solved value;
        const unsigned last = std::min<std::size_t>(k, r.lags.size() - 1);
        for (unsigned lag = 1; lag <= last; ++lag) {
          for (const auto &[j, coefficient] : r.lags[lag])
            if (meta[k - lag][j].bound >= 0)
              value.meta.merge(meta[k - lag][j]);
          if (cleared && !r.derivative[lag].zero() && meta[k - lag][row].bound >= 0)
            value.meta.merge(meta[k - lag][row]);
        }
        if (value.meta.bound >= 0) {
          value.meta.known = std::min(value.meta.known, cap);
          const long length = value.meta.known + pole;
          if (!structure_only && length > 0) {
            for (unsigned lag = 1; lag <= last; ++lag) {
              const auto &source = orders[k - lag];
              for (const auto &[j, coefficient] : r.lags[lag]) {
                if (source[j].empty())
                  continue;
                if (cleared && !unit) {
                  fmpz_poly_scalar_mul_fmpz(weight.raw(), coefficient.raw(), rho.value);
                  accumulate_low(value.numerator, weight, source[j], 0, length);
                } else
                  accumulate_low(value.numerator, coefficient, source[j], 0, length);
              }
              const auto &previous = source[row];
              if (!cleared || r.derivative[lag].zero() || previous.empty())
                continue;
              fmpz_poly_get_coeff_fmpz(shift.value, exponent.raw(), 0);
              fmpz_addmul_ui(shift.value, rho.value, k - lag);
              fmpz_poly_set(weight.raw(), exponent.raw());
              fmpz_poly_set_coeff_fmpz(weight.raw(), 0, shift.value);
              fmpz_poly_mul(weight.raw(), weight.raw(), r.derivative[lag].raw());
              fmpz_poly_neg(weight.raw(), weight.raw());
              accumulate_low(value.numerator, weight, previous, 0, length);
              fmpz_poly_scalar_mul_fmpz(weight.raw(), r.derivative[lag].raw(), rho.value);
              fmpz_poly_neg(weight.raw(), weight.raw());
              accumulate_low(value.numerator, weight, previous, 1, length);
            }
            trim(value.numerator);
            // The row equation was multiplied by its scale (and by rho).
            if (!value.numerator.empty()) {
              if (cleared && !unit)
                fmpz_mul(value.denominator.value, value.denominator.value, rho.value);
              divide_unit(value, T.row_units[row].unit, inverses[row], length);
            }
          }
          value.rel = T.row_units[row].valuation;
          value.meta.known -= T.row_units[row].valuation;
        }
        if (successor[row] >= 0 && solved[successor[row]].meta.bound >= 0)
          value = sum(value, solved[successor[row]], pole, cap, structure_only);
        solve_low(value,
                  eigen[column].power + Rational(k) - eigen[row].power,
                  (eigen[column].slope - eigen[row].slope) * T.scale_rational,
                  pole, cap);
        if (k == 0 && row == column) {
          if (value.meta.bound >= 0)
            throw std::logic_error("seeded Frobenius coefficient has a source");
          value = Solved{};
          value.meta.bound = 0;
          if (!structure_only)
            value.numerator.emplace_back(1);
        }
        if (value.numerator.size() > ladder_bound)
          throw std::logic_error("affine Frobenius log ladder bound exceeded");
        reduce_low(value);
        solved[row] = std::move(value);
      }
      long rise = 0;
      fmpz_one(step.value);
      for (const auto &value : solved)
        if (value.meta.bound >= 0) {
          rise = std::max(rise, value.rel);
          if (!value.numerator.empty())
            fmpz_lcm(step.value, step.value, value.denominator.value);
        }
      for (unsigned row = 0; row < d; ++row) {
        auto &value = solved[row];
        meta[k][row] = value.meta;
        if (value.numerator.empty())
          continue;
        fmpz_divexact(cofactor.value, step.value, value.denominator.value);
        const long length = std::min(value.meta.known, cap) + pole + rise;
        for (auto &c : value.numerator) {
          if (!fmpz_is_one(cofactor.value))
            fmpz_poly_scalar_mul_fmpz(c.raw(), c.raw(), cofactor.value);
          if (rise > value.rel)
            fmpz_poly_shift_left(c.raw(), c.raw(), rise - value.rel);
          fmpz_poly_truncate(c.raw(), std::max(0l, length));
        }
        trim(value.numerator);
        orders[k][row] = std::move(value.numerator);
      }
      if (rise || !fmpz_is_one(step.value)) {
        for (unsigned old = k + 1 > window_ ? k + 1 - window_ : 0; old < k; ++old)
          for (auto &value : orders[old])
            for (auto &c : value) {
              if (!fmpz_is_one(step.value))
                fmpz_poly_scalar_mul_fmpz(c.raw(), c.raw(), step.value);
              if (rise)
                fmpz_poly_shift_left(c.raw(), c.raw(), rise);
            }
        fmpz_mul(common.value, common.value, step.value);
        pole += rise;
      }
      Value physical;
      for (unsigned row = 0; row < d; ++row) {
        Meta known;
        for (const auto &[j, coefficient] : T.frame[row].entries)
          if (meta[k][j].bound >= 0)
            known.merge(meta[k][j]);
        if (known.bound < 0)
          continue;
        known.known = std::min(known.known, cap);
        physical.clear();
        if (!structure_only && known.known + pole > 0)
          for (const auto &[j, coefficient] : T.frame[row].entries)
            if (!orders[k][j].empty())
              accumulate_low(physical, coefficient, orders[k][j], 0, known.known + pole);
        trim(physical);
        emit(k, row, physical, T.frame_units[row], common.value, pole, known);
      }
      if (k >= window_) {
        orders[k - window_] = std::vector<Value>(d);
      }
    }
  }

private:
  struct Row {
    Z scale; // scale*g[lag][row][*] and scale*q~[lag] are integral
    std::vector<std::vector<std::pair<unsigned, Z>>> lags; // by lag: (column, coefficient)
    std::vector<Z> derivative; // by lag; empty for a Taylor (uncleared) row
  };
  struct FrameRow {
    Z scale;
    std::vector<std::pair<unsigned, Z>> entries;
  };
  // numerator/(denominator*eps^pole), relative to the shared denominator.
  struct Fraction {
    Value numerator;
    Z denominator{1};
    long pole = 0;
  };
  std::vector<Row> rows_;
  std::vector<FrameRow> frame_;
  unsigned window_ = 0;
  // The recurrence in t = eps/scale: the coefficients of rows_ and frame_
  // with eps replaced by scale*t, and every scale as t^valuation times a unit.
  struct Truncation {
    Integer scale;
    Rational scale_rational;
    std::vector<Row> rows;
    std::vector<FrameRow> frame;
    std::vector<Unit> row_units, frame_units;
  };
  std::unique_ptr<Truncation> truncation_;
  // numerator/(denominator * t^rel) relative to the shared denominator; valid
  // below t^(meta.known + pole + rel).
  struct Solved {
    Value numerator;
    Integer denominator{1};
    long rel = 0;
    Meta meta;
  };
  // numerator/denominator = 1/unit modulo t^length.
  struct Inverse {
    Z numerator;
    Integer denominator;
    slong length = 0;
  };
  static void mul_low(fmpz_poly_struct *out, const fmpz_poly_struct *a,
                      const fmpz_poly_struct *b, long length) {
    const long full = a->length + b->length - 1;
    if (length <= 0 || !a->length || !b->length)
      fmpz_poly_zero(out);
    else if (b->length == 1)
      fmpz_poly_scalar_mul_fmpz(out, a, b->coeffs), fmpz_poly_truncate(out, length);
    else
      fmpz_poly_mullow(out, a, b, std::min(length, full));
  }
  // target += source*coefficient modulo t^length, for a short coefficient.
  static void addmul_low(fmpz_poly_struct *target, const fmpz_poly_struct *source,
                         const fmpz_poly_struct *coefficient, long length) {
    length = std::min<long>(length, source->length + coefficient->length - 1);
    if (length <= 0)
      return;
    fmpz_poly_fit_length(target, length);
    if (target->length < length) {
      _fmpz_vec_zero(target->coeffs + target->length, length - target->length);
      _fmpz_poly_set_length(target, length);
    }
    for (slong e = 0; e < coefficient->length && e < length; ++e)
      if (!fmpz_is_zero(coefficient->coeffs + e))
        _fmpz_vec_scalar_addmul_fmpz(target->coeffs + e, source->coeffs,
                                     std::min<slong>(source->length, length - e),
                                     coefficient->coeffs + e);
    _fmpz_poly_normalise(target);
  }
  static void accumulate_low(Value &out, const Z &coefficient, const Value &source,
                             unsigned drop, long length) {
    if (coefficient.zero() || source.size() <= drop || length <= 0)
      return;
    if (out.size() < source.size() - drop)
      out.resize(source.size() - drop);
    Z scaled;
    for (std::size_t l = drop; l < source.size(); ++l) {
      if (source[l].zero())
        continue;
      if (drop) {
        fmpz_poly_scalar_mul_ui(scaled.raw(), coefficient.raw(), l);
        addmul_low(out[l - drop].raw(), source[l].raw(), scaled.raw(), length);
      } else
        addmul_low(out[l].raw(), source[l].raw(), coefficient.raw(), length);
    }
  }
  static void divide_unit(Solved &value, const Z &unit, Inverse &inverse, long length) {
    if (fmpz_poly_length(unit.raw()) == 1) {
      const fmpz *c = unit.raw()->coeffs;
      if (fmpz_is_one(c))
        return;
      if (fmpz_sgn(c) < 0)
        for (auto &n : value.numerator)
          fmpz_poly_neg(n.raw(), n.raw());
      Integer magnitude;
      fmpz_abs(magnitude.value, c);
      fmpz_mul(value.denominator.value, value.denominator.value, magnitude.value);
      return;
    }
    if (inverse.length < length) {
      inverse.length = length + 8;
      fmpq_poly_t q, r;
      fmpq_poly_init(q);
      fmpq_poly_init(r);
      fmpq_poly_set_fmpz_poly(q, unit.raw());
      fmpq_poly_inv_series(r, q, inverse.length);
      fmpq_poly_get_numerator(inverse.numerator.raw(), r);
      fmpz_set(inverse.denominator.value, fmpq_poly_denref(r));
      fmpq_poly_clear(q);
      fmpq_poly_clear(r);
    }
    for (auto &n : value.numerator)
      mul_low(n.raw(), n.raw(), inverse.numerator.raw(), length);
    fmpz_mul(value.denominator.value, value.denominator.value, inverse.denominator.value);
    trim(value.numerator);
  }
  static Solved sum(const Solved &a, const Solved &b, long pole, long cap,
                    bool structure_only) {
    if (a.meta.bound < 0)
      return b;
    if (b.meta.bound < 0)
      return a;
    Solved out;
    out.meta = a.meta;
    out.meta.merge(b.meta);
    out.meta.known = std::min(out.meta.known, cap);
    out.rel = std::max(a.rel, b.rel);
    if (structure_only)
      return out;
    const long length = out.meta.known + pole + out.rel;
    out.numerator.resize(std::max(a.numerator.size(), b.numerator.size()));
    Z term;
    const auto add = [&](const Solved &value, const Solved &other) {
      for (std::size_t l = 0; l < value.numerator.size(); ++l) {
        fmpz_poly_scalar_mul_fmpz(term.raw(), value.numerator[l].raw(), other.denominator.value);
        fmpz_poly_shift_left(term.raw(), term.raw(), out.rel - value.rel);
        fmpz_poly_add(out.numerator[l].raw(), out.numerator[l].raw(), term.raw());
      }
    };
    add(a, b);
    add(b, a);
    for (auto &c : out.numerator)
      fmpz_poly_truncate(c.raw(), std::max(0l, length));
    fmpz_mul(out.denominator.value, a.denominator.value, b.denominator.value);
    trim(out.numerator);
    return out;
  }
  // (delta + d/dL) out = rhs for delta = m + s*t.
  static void solve_low(Solved &value, const Rational &m, const Rational &s,
                        long pole, long cap) {
    if (value.meta.bound < 0)
      return;
    value.meta.known = std::min(value.meta.known, cap);
    auto &rhs = value.numerator;
    Integer rho;
    const Z lambda = linear(m, s, rho); // rho*delta
    if (lambda.zero()) {
      // Exact resonance: integrate in L with zero integration constant.
      ++value.meta.bound;
      if (rhs.empty())
        return;
      const std::size_t top = rhs.size() - 1;
      Integer weights(1), factor;
      for (std::size_t l = 2; l <= top + 1; ++l) {
        fmpz_set_ui(factor.value, l);
        fmpz_lcm(weights.value, weights.value, factor.value);
      }
      rhs.emplace_back();
      for (std::size_t l = top + 1; l-- > 0;) {
        fmpz_poly_scalar_mul_fmpz(rhs[l + 1].raw(), rhs[l].raw(), weights.value);
        fmpz_poly_scalar_divexact_ui(rhs[l + 1].raw(), rhs[l + 1].raw(), l + 1);
      }
      fmpz_poly_zero(rhs[0].raw());
      fmpz_mul(value.denominator.value, value.denominator.value, weights.value);
      return;
    }
    // delta = t^a * unit; a = 1 where the powers resonate and the slopes differ.
    const long a = m.is_zero() ? 1 : 0;
    const long length = value.meta.known + pole + value.rel;
    value.meta.known -= a * (value.meta.bound + 1);
    if (rhs.empty())
      return;
    const std::size_t top = rhs.size() - 1;
    // 1/delta = inverse(t)/(divisor*t^a).
    Z inverse;
    Integer divisor, c0, c1, power;
    fmpz_poly_get_coeff_fmpz(c0.value, lambda.raw(), 0);
    fmpz_poly_get_coeff_fmpz(c1.value, lambda.raw(), 1);
    if (a || fmpz_is_zero(c1.value) || length <= 1) {
      fmpz_poly_set_fmpz(inverse.raw(), rho.value);
      fmpz_set(divisor.value, a ? c1.value : c0.value);
    } else if (fmpz_divisible(c1.value, c0.value)) {
      // (1/c0) * sum (-c1/c0)^j t^j.
      fmpz_divexact(c1.value, c1.value, c0.value);
      fmpz_neg(c1.value, c1.value);
      fmpz_set(power.value, rho.value);
      for (long j = 0; j < length; ++j) {
        fmpz_poly_set_coeff_fmpz(inverse.raw(), j, power.value);
        fmpz_mul(power.value, power.value, c1.value);
      }
      fmpz_set(divisor.value, c0.value);
    } else {
      // sum (-c1)^j c0^(length-1-j) t^j / c0^length.
      fmpz_neg(c1.value, c1.value);
      fmpz_set(power.value, rho.value);
      for (long j = length; j-- > 0;) {
        fmpz_poly_set_coeff_fmpz(inverse.raw(), j, power.value);
        fmpz_mul(power.value, power.value, c0.value);
      }
      fmpz_one(power.value);
      for (long j = 0; j < length; ++j) {
        fmpz_mul(inverse.raw()->coeffs + j, inverse.raw()->coeffs + j, power.value);
        fmpz_mul(power.value, power.value, c1.value);
      }
      _fmpz_poly_normalise(inverse.raw());
      fmpz_pow_ui(divisor.value, c0.value, length);
    }
    if (fmpz_sgn(divisor.value) < 0) {
      fmpz_neg(divisor.value, divisor.value);
      fmpz_poly_neg(inverse.raw(), inverse.raw());
    }
    // out[l] = T[l]/(divisor*t^a)^(top+1-l), T[top] = rhs[top]*inverse,
    // T[l] = (rhs[l]*(divisor*t^a)^(top-l) - (l+1)*T[l+1])*inverse.
    std::vector<Integer> powers(top + 2, Integer(1));
    for (std::size_t i = 1; i < powers.size(); ++i)
      fmpz_mul(powers[i].value, powers[i - 1].value, divisor.value);
    Z above, current;
    for (std::size_t l = top + 1; l-- > 0;) {
      fmpz_poly_scalar_mul_fmpz(current.raw(), rhs[l].raw(), powers[top - l].value);
      if (a && top > l)
        fmpz_poly_shift_left(current.raw(), current.raw(), a * (top - l));
      if (l < top) {
        fmpz_poly_scalar_mul_ui(above.raw(), above.raw(), l + 1);
        fmpz_poly_sub(current.raw(), current.raw(), above.raw());
      }
      mul_low(current.raw(), current.raw(), inverse.raw(), length);
      fmpz_poly_scalar_mul_fmpz(rhs[l].raw(), current.raw(), powers[l].value);
      if (a && l)
        fmpz_poly_shift_left(rhs[l].raw(), rhs[l].raw(), a * l);
      above = current;
    }
    fmpz_mul(value.denominator.value, value.denominator.value, powers[top + 1].value);
    value.rel += a * (top + 1);
    const long valid = value.meta.known + pole + value.rel;
    for (auto &c : rhs)
      fmpz_poly_truncate(c.raw(), std::max(0l, valid));
    trim(rhs);
  }
  // The integer denominator against the content of the numerators.
  static void reduce_low(Solved &value) {
    if (value.numerator.empty()) {
      fmpz_one(value.denominator.value);
      return;
    }
    Integer g(value.denominator), content;
    for (const auto &c : value.numerator) {
      if (fmpz_is_one(g.value))
        return;
      if (c.zero())
        continue;
      fmpz_poly_content(content.value, c.raw());
      fmpz_gcd(g.value, g.value, content.value);
    }
    if (fmpz_is_one(g.value))
      return;
    fmpz_divexact(value.denominator.value, value.denominator.value, g.value);
    for (auto &c : value.numerator)
      fmpz_poly_scalar_divexact_fmpz(c.raw(), c.raw(), g.value);
  }
  static Z integral(const Exact &value, std::size_t ei) {
    Z numerator, denominator;
    value.univariate_polynomials(numerator.raw(), denominator.raw(), ei);
    if (!denominator.one())
      throw std::logic_error(
          "fraction-free Frobenius coefficient kept a denominator");
    return numerator;
  }
  static void trim(Value &value) {
    while (!value.empty() && value.back().zero())
      value.pop_back();
  }
  static void scale(Z &value, const Z &factor, long shift) {
    if (value.zero())
      return;
    if (!factor.one())
      fmpz_poly_mul(value.raw(), value.raw(), factor.raw());
    if (shift)
      fmpz_poly_shift_left(value.raw(), value.raw(), shift);
  }
  // target += source*coefficient for a short coefficient: one scalar pass per
  // coefficient term, with no temporary product.
  static void addmul(fmpz_poly_struct *target, const fmpz_poly_struct *source,
                     const fmpz_poly_struct *coefficient) {
    const slong length = source->length + coefficient->length - 1;
    fmpz_poly_fit_length(target, length);
    if (target->length < length) {
      _fmpz_vec_zero(target->coeffs + target->length, length - target->length);
      _fmpz_poly_set_length(target, length);
    }
    for (slong e = 0; e < coefficient->length; ++e)
      if (!fmpz_is_zero(coefficient->coeffs + e))
        _fmpz_vec_scalar_addmul_fmpz(target->coeffs + e, source->coeffs,
                                     source->length, coefficient->coeffs + e);
    _fmpz_poly_normalise(target);
  }
  // out[l-drop] += coefficient*(drop ? l : 1)*source[l]; drop=1 is d/dL.
  static void accumulate(Value &out, const Z &coefficient, const Value &source,
                         unsigned drop) {
    if (coefficient.zero() || source.size() <= drop)
      return;
    if (out.size() < source.size() - drop)
      out.resize(source.size() - drop);
    Z scaled;
    for (std::size_t l = drop; l < source.size(); ++l) {
      if (source[l].zero())
        continue;
      if (drop) {
        fmpz_poly_scalar_mul_ui(scaled.raw(), coefficient.raw(), l);
        addmul(out[l - drop].raw(), source[l].raw(), scaled.raw());
      } else
        addmul(out[l].raw(), source[l].raw(), coefficient.raw());
    }
  }
  // a + b over the product of their denominators (Jordan successor only).
  static Fraction sum(const Fraction &a, const Fraction &b) {
    Fraction out;
    out.pole = std::max(a.pole, b.pole);
    out.numerator.resize(std::max(a.numerator.size(), b.numerator.size()));
    Z term;
    const auto add = [&](const Fraction &value, const Fraction &other) {
      for (std::size_t l = 0; l < value.numerator.size(); ++l) {
        term = value.numerator[l];
        scale(term, other.denominator, out.pole - value.pole);
        fmpz_poly_add(out.numerator[l].raw(), out.numerator[l].raw(),
                      term.raw());
      }
    };
    add(a, b);
    add(b, a);
    fmpz_poly_mul(out.denominator.raw(), a.denominator.raw(),
                  b.denominator.raw());
    trim(out.numerator);
    return out;
  }
  // (delta + d/dL) out = rhs for delta = m + s*eps, kept in integers.
  static void solve(Fraction &value, const Rational &m, const Rational &s) {
    auto &rhs = value.numerator;
    if (rhs.empty()) {
      value = Fraction{};
      return;
    }
    Integer rho;
    const Z lambda = linear(m, s, rho); // rho*delta
    const std::size_t top = rhs.size() - 1;
    auto *denominator = value.denominator.raw();
    if (lambda.zero()) {
      // Exact resonance: integrate in L with zero integration constant.
      Integer weights, factor;
      fmpz_one(weights.value);
      for (std::size_t l = 2; l <= top + 1; ++l) {
        fmpz_set_ui(factor.value, l);
        fmpz_lcm(weights.value, weights.value, factor.value);
      }
      rhs.emplace_back();
      for (std::size_t l = top + 1; l-- > 0;) {
        fmpz_poly_scalar_mul_fmpz(rhs[l + 1].raw(), rhs[l].raw(), weights.value);
        fmpz_poly_scalar_divexact_ui(rhs[l + 1].raw(), rhs[l + 1].raw(), l + 1);
      }
      fmpz_poly_zero(rhs[0].raw());
      fmpz_poly_scalar_mul_fmpz(denominator, denominator, weights.value);
      return;
    }
    if (!top) {
      fmpz_poly_scalar_mul_fmpz(rhs[0].raw(), rhs[0].raw(), rho.value);
      fmpz_poly_mul(denominator, denominator, lambda.raw());
      return;
    }
    // out[l] = T[l]/lambda^(top+1-l) with T[top] = rho*rhs[top] and
    // T[l] = rho*(rhs[l]*lambda^(top-l) - (l+1)*T[l+1]).
    std::vector<Z> power(top + 2, Z(1));
    for (std::size_t i = 1; i < power.size(); ++i)
      fmpz_poly_mul(power[i].raw(), power[i - 1].raw(), lambda.raw());
    Z above, current;
    for (std::size_t l = top + 1; l-- > 0;) {
      fmpz_poly_mul(current.raw(), rhs[l].raw(), power[top - l].raw());
      if (l < top) {
        fmpz_poly_scalar_mul_ui(above.raw(), above.raw(), l + 1);
        fmpz_poly_sub(current.raw(), current.raw(), above.raw());
      }
      fmpz_poly_scalar_mul_fmpz(current.raw(), current.raw(), rho.value);
      fmpz_poly_mul(rhs[l].raw(), current.raw(), power[l].raw());
      above = current;
    }
    fmpz_poly_mul(denominator, denominator, power[top + 1].raw());
    trim(rhs);
  }
  // Epsilon powers move into the pole; the rest of the denominator is small,
  // so its gcd with the large numerators is cheap and keeps the shared
  // denominator free of factors that this row never needed.
  static void reduce(Fraction &value) {
    auto &numerator = value.numerator;
    if (numerator.empty()) {
      value = Fraction{};
      return;
    }
    auto *denominator = value.denominator.raw();
    const slong below = value.denominator.valuation();
    if (below)
      fmpz_poly_shift_right(denominator, denominator, below);
    slong above = -1;
    for (const auto &c : numerator)
      if (!c.zero())
        above = above < 0 ? c.valuation() : std::min(above, c.valuation());
    if (above > 0)
      for (auto &c : numerator)
        fmpz_poly_shift_right(c.raw(), c.raw(), above);
    value.pole += below - above;
    Z g(value.denominator);
    for (const auto &c : numerator) {
      if (g.one())
        break;
      if (!c.zero())
        fmpz_poly_gcd(g.raw(), g.raw(), c.raw());
    }
    if (fmpz_sgn(fmpz_poly_lead(denominator)) < 0)
      fmpz_poly_neg(g.raw(), g.raw());
    if (g.one())
      return;
    fmpz_poly_div(denominator, denominator, g.raw());
    for (auto &c : numerator)
      if (fmpz_poly_length(g.raw()) == 1)
        fmpz_poly_scalar_divexact_fmpz(c.raw(), c.raw(), g.raw()->coeffs);
      else
        fmpz_poly_div(c.raw(), c.raw(), g.raw());
  }
};
} // namespace affine_frobenius_detail

class AffineFrobeniusSeries {
public:
  using Matrix = fuchsify::Matrix;
  using B = kernel::ComplexBall;
  using Exponent = affine_frobenius_detail::Exponent;
  struct Options {
    unsigned max_dimension = 12, max_x_order = 128, max_terms = 250000,
             max_epsilon_depth = 512;
    // Clear polynomial row denominators in the constant Jordan frame, so a
    // row reaches back only deg(denominator) orders. False, or a row beyond
    // max_clearing_degree, uses the Taylor coefficients of the connection at
    // every lag instead; both feed the same fraction-free recurrence.
    bool finite_lag_recurrence = true;
    unsigned max_clearing_degree = 32;
    std::size_t max_projection_products=100000000;
    // Called with (columns completed, columns) under a lock; the columns are
    // solved on `threads` threads (0: one per hardware thread).
    std::function<void(unsigned,unsigned)> column_progress;
    unsigned threads = 0;
    // Zero: every coefficient is the exact rational function of epsilon.
    // Otherwise the coefficients are carried as truncated Laurent series in
    // epsilon: column c is known through at least this many orders, counted
    // from the lowest epsilon order of column c of the residue frame, and
    // Expansion::epsilon_truncation states what is known of every entry.
    unsigned epsilon_orders = 0;
    // The system is the recursion at one numeric dimension (no epsilon left).
    // Its exponents a+b*shift then carry no slope, and the dimensionally
    // regulated sectors (b nonzero) are recognised by a non-integer power at a
    // generic shift. Such columns are labelled with a unit slope, so every
    // endpoint rule treats them as regulated; the recurrence itself uses the
    // true exponents. Only order zero in epsilon is meaningful afterwards.
    bool regulate_noninteger_powers = false;
  };
  struct Term {
    unsigned row, column, log_degree;
    Rational power, slope;
    Exact coefficient;
  };
  struct Expansion {
    unsigned rows, columns;
    std::vector<Term> terms;
    bool omitted_tail_certified = false;
    // Full-system Wronskian prefactor, never inferred from the finite series.
    std::optional<Exact> wronskian_prefactor;
    bool coherent_x_frontier = false;
    // Row i is complete in column c only through the absolute power
    // absolute_x_frontiers()[c] + row_x_valuations[i] (Valuations::zero_valuation
    // marks an exact zero row). Empty when no such frontier is known, e.g. for
    // contracted columns; beyond it coefficients are unknown, never zero.
    std::vector<long> row_x_valuations;
    // Empty: the coefficients are exact. Otherwise entry [row][column] is the
    // epsilon order T from which that entry is unknown: each coefficient is
    // the Laurent polynomial of the true one modulo eps^T, and a term that is
    // absent is O(eps^T), not zero. exact_order marks an exact entry.
    std::vector<std::vector<long>> epsilon_truncation;
    static constexpr long exact_order = affine_frobenius_detail::EpsilonSeries::exact_order;
    long truncation(unsigned row, unsigned column) const {
      return epsilon_truncation.empty() ? exact_order : epsilon_truncation[row][column];
    }
  };
  struct Valuations {
    std::vector<long> minimum_by_column;
    unsigned maximum_pole = 0;
    std::vector<long> required_source_top(long desired_top) const {
      std::vector<long> out;
      for (auto v : minimum_by_column)
        out.push_back(v == zero_valuation ? 0 : desired_top - v);
      return out;
    }
    static constexpr long zero_valuation = std::numeric_limits<long>::max() / 4;
  };
  struct LaurentMatrix {
    // An empty window, so a result returned on a failure path is well defined.
    int low = 0, high = -1;
    std::vector<std::vector<std::vector<B>>> coefficients;
  };
  static constexpr bool omitted_tail_certified = false;
  unsigned dimension() const { return d_; }
  unsigned x_order() const { return n_; }
  unsigned finite_lag_rows() const { return finite_lag_rows_; }
  // Coefficients are truncated Laurent series in epsilon (Options::epsilon_orders).
  bool truncated() const { return !expansion_.epsilon_truncation.empty(); }
  const Matrix &residue_frame() const { return frame_; }
  const std::vector<Exponent> &exponents() const { return eigen_; }
  const std::vector<int> &jordan_successors() const { return successor_; }
  static AffineFrobeniusSeries prepare(const Matrix &a, std::size_t xi,
                                       std::size_t ei, unsigned n) {
    return prepare(a, xi, ei, n, Options{});
  }
  static AffineFrobeniusSeries prepare(const Matrix &a, std::size_t xi,
                                       std::size_t ei, unsigned n,
                                       const Options &options) {
    return prepare_impl(a,xi,ei,n,options);
  }
private:
  static AffineFrobeniusSeries prepare_impl(const Matrix &a,std::size_t xi,
      std::size_t ei,unsigned n,const Options &options) {
    using namespace affine_frobenius_detail;
    using namespace fuchsify::detail;
    if (a.empty() || a[0].empty() || a.size() > options.max_dimension ||
        n > options.max_x_order || xi == ei || !options.max_terms ||
        !options.max_epsilon_depth)
      throw std::invalid_argument(
          "invalid affine Frobenius dimensions or finite budgets");
    AffineFrobeniusSeries out(a[0][0].constant(0));
    out.d_ = a.size();
    out.n_ = n;
    out.xi_ = xi;
    out.ei_ = ei;
    out.options_ = options;
    unsigned d = out.d_;
    auto z = out.zero_;
    auto x = z.variable(xi);
    (void)z.variable(ei);
    std::vector<Matrix> g(n + 1, zeros(d, d, z));
    auto connection = zeros(d, d, z);
    for (unsigned i = 0; i < d; ++i) {
      if (a[i].size() != d)
        throw std::invalid_argument("affine Frobenius matrix must be square");
      for (unsigned j = 0; j < d; ++j) {
        for (const auto &ts :
             {a[i][j].numerator_terms(), a[i][j].denominator_terms()})
          for (auto &t : ts)
            for (unsigned v = 0; v < t.powers.size(); ++v)
              if (v != xi && v != ei && t.powers[v])
                throw std::invalid_argument(
                    "unsubstituted affine Frobenius parameter");
        connection[i][j] = x * a[i][j];
        g[0][i][j] = taylor(connection[i][j], xi, 0)[0];
      }
    }
    auto prepared = affine_frobenius_detail::frame(g[0], xi, ei);
    out.frame_ = prepared.transform;
    out.eigen_ = prepared.exponents;
    out.successor_ = prepared.successor;
    out.slopes_.clear();
    for (const auto &exponent : out.eigen_)
      out.slopes_.push_back(options.regulate_noninteger_powers && exponent.slope.is_zero() &&
                            !is_integer(exponent.power) ? Rational(1) : exponent.slope);
    // Constant Laurent column relations can mix exactly those base powers
    // that differ by integers (epsilon expansion makes slopes logarithms).
    // A common absolute cutoff commutes with these relations; a separate
    // number of Taylor coefficients per column does not.
    out.frontiers_.resize(d);
    for (unsigned c = 0; c < d; ++c) {
      auto minimum = out.eigen_[c].power;
      for (unsigned j = 0; j < d; ++j)
        if (is_integer(out.eigen_[j].power - out.eigen_[c].power))
          minimum = std::min(minimum, out.eigen_[j].power);
      out.frontiers_[c] = minimum + Rational(n);
      if (out.frontiers_[c] < out.eigen_[c].power)
        throw std::domain_error("affine x order does not span the absolute "
                                "exponent frontier; increase x order");
    }
    auto inv = inverse(out.frame_);
    std::vector<std::vector<Exact>> derivative_lags(d);
    if (!options.finite_lag_recurrence) {
      // Retain the original rational Taylor convolution as a reference path.
      for (unsigned i = 0; i < d; ++i)
        for (unsigned j = 0; j < d; ++j) {
          auto coefficients = taylor(connection[i][j], xi, n);
          for (unsigned k = 1; k <= n; ++k)
            g[k][i][j] = coefficients[k];
        }
      for (auto &m : g)
        m = multiply(multiply(inv, m), out.frame_);
    } else {
      g[0] = multiply(multiply(inv, g[0]), out.frame_);
      connection = multiply(multiply(inv, connection), out.frame_);
      // In the constant Jordan frame, q_i(x)*theta(v_i)=sum_j p_ij(x)*v_j.
      // Its zero lag is still theta-J; all positive denominator lags involve
      // previously known coefficients and their logarithmic derivatives.
      for (unsigned row = 0; row < d; ++row) {
        bool cleared = n > 0;
        auto denominator = z.constant(1);
        auto degree = [&](const Exact &value) {
          unsigned result = 0;
          for (const auto &term : value.numerator_terms())
            result = std::max(result, static_cast<unsigned>(term.powers[xi]));
          return result;
        };
        unsigned lag_degree = 0;
        if (cleared)
          for (unsigned j = 0; j < d; ++j) {
            denominator = denominator.polynomial_lcm(connection[row][j].denominator());
            lag_degree = degree(denominator);
            if (lag_degree > options.max_clearing_degree ||
                denominator.numerator_terms().size() > options.max_terms) {
              cleared = false;
              break;
            }
          }
        std::vector<Exact> numerators;
        if (cleared)
          for (unsigned j = 0; j < d; ++j) {
            numerators.push_back(denominator * connection[row][j]);
            lag_degree = std::max(lag_degree, degree(numerators.back()));
            if (lag_degree > options.max_clearing_degree ||
                numerators.back().numerator_terms().size() > options.max_terms) {
              cleared = false;
              break;
            }
          }
        if (cleared) {
          derivative_lags[row] = taylor(denominator, xi, lag_degree);
          const auto leading = derivative_lags[row][0];
          if (leading.is_zero())
            throw std::logic_error("cleared Frobenius denominator is not x-analytic");
          for (auto &q : derivative_lags[row])
            q = q / leading;
          for (unsigned j = 0; j < d; ++j) {
            auto coefficients = taylor(numerators[j], xi, lag_degree);
            for (unsigned k = 1; k <= std::min(lag_degree,n); ++k)
              g[k][row][j] = coefficients[k] / leading;
            if (coefficients[0] / leading != g[0][row][j])
              throw std::logic_error("cleared Frobenius residue changed");
          }
          ++out.finite_lag_rows_;
        } else {
          for (unsigned j = 0; j < d; ++j) {
            auto coefficients = taylor(connection[row][j], xi, n);
            for (unsigned k = 1; k <= n; ++k)
              g[k][row][j] = coefficients[k];
          }
        }
      }
    }
    out.expansion_ = {d, d, {}};
    out.expansion_.coherent_x_frontier = true;
    out.expansion_.row_x_valuations.assign(d, 0);
    auto regular_trace = z;
    for (unsigned i = 0; i < d; ++i)
      regular_trace = regular_trace + a[i][i] - g[0][i][i] / x;
    // The regular Wronskian exponential is an epsilon unit only when its
    // integrand is jointly analytic at (x,eps)=(0,0). Generic epsilon
    // valuation is insufficient: 1/(x+eps) has a moving pole.
    std::vector<Exact> origin;
    for (unsigned i = 0; i < z.variable_count(); ++i)
      origin.push_back(i == xi || i == ei ? z : z.variable(i));
    if (!regular_trace.denominator().substitute(origin).is_zero())
      out.expansion_.wronskian_prefactor =
          affine_frobenius_detail::determinant(out.frame_);
    IntegerRecurrence recurrence(g, derivative_lags, out.frame_, ei);
    const bool truncated = options.epsilon_orders != 0;
    out.log_bounds_.assign(d, 0);
    if (truncated) {
      recurrence.enable_truncation(out.eigen_, n);
      out.expansion_.epsilon_truncation.assign(d, std::vector<long>(d, Expansion::exact_order));
    }
    // Columns are independent, so they are solved on separate threads and
    // their terms appended in column order: the result does not depend on the
    // thread count. Orders beyond a column's absolute frontier cannot
    // contribute to any retained coefficient.
    std::vector<unsigned> tops(d);
    for (unsigned column = 0; column < d; ++column) {
      const auto cutoff=out.frontiers_[column]-out.eigen_[column].power;
      if(!is_integer(cutoff)||cutoff.sign()<0)throw std::logic_error("nonintegral affine column frontier");
      if(cutoff>Rational(n))throw std::logic_error("affine column frontier exceeds prepared order");
      tops[column]=static_cast<unsigned>(integer_value(cutoff));
    }
    std::vector<std::vector<Term>> columns(d);
    std::atomic<unsigned> next{0};
    std::atomic<bool> failed{false};
    std::exception_ptr failure;
    std::mutex shared_state;
    unsigned completed = 0;
    const auto solve_columns = [&] {
      IntegerPolynomial numerator, denominator, shared, divisor;
      for (unsigned column; !failed && (column = next++) < d;) try {
        auto &terms = columns[column];
        if (truncated) {
          using S = EpsilonSeries;
          using R = IntegerRecurrence;
          long lead = S::exact_order;
          for (unsigned i = 0; i < d; ++i)
            if (!out.frame_[i][column].is_zero())
              lead = std::min(lead, fuchsify::detail::valuation(out.frame_[i][column], ei));
          // How many orders the column loses, from a pass without numbers.
          const long probe = 1l << 20;
          long lowest = S::exact_order;
          int logs = 0;
          recurrence.truncated_column(out.eigen_, out.successor_, column, tops[column],
            static_cast<std::size_t>(d) * (n + 1), probe, true,
            [&](unsigned, unsigned, const R::Value &, const R::Unit &unit, const fmpz *, long,
                const R::Meta &meta) {
              lowest = std::min(lowest, meta.known - static_cast<long>(unit.valuation));
              logs = std::max(logs, meta.bound);
            });
          out.log_bounds_[column] = logs;
          const long cap = lead + static_cast<long>(options.epsilon_orders) + (probe - lowest);
          struct Pending { unsigned k, row, log; S value; };
          std::vector<Pending> pending;
          std::vector<long> known(d, S::exact_order);
          recurrence.truncated_column(out.eigen_, out.successor_, column, tops[column],
            static_cast<std::size_t>(d) * (n + 1), cap, false,
            [&](unsigned k, unsigned row, const R::Value &values, const R::Unit &unit,
                const fmpz *common, long pole, const R::Meta &meta) {
              const long order = meta.known - static_cast<long>(unit.valuation);
              known[row] = std::min(known[row], order);
              for (unsigned l = 0; l < values.size(); ++l) {
                if (values[l].zero())
                  continue;
                auto value = S::from_integers(values[l].raw(), common,
                                              -(pole + static_cast<long>(unit.valuation)), order);
                value.divide_unit(unit.unit.raw());
                value.divide_variable(recurrence.scale());
                pending.push_back({k, row, l, std::move(value)});
              }
              if (pending.size() > options.max_terms)
                throw std::length_error("affine Frobenius term budget exhausted");
            });
          for (unsigned row = 0; row < d; ++row)
            out.expansion_.epsilon_truncation[row][column] = known[row];
          for (auto &entry : pending) {
            entry.value.truncate(known[entry.row]);
            if (!entry.value.is_zero())
              terms.push_back({entry.row, column, entry.log,
                               out.eigen_[column].power + Rational(entry.k),
                               out.slopes_[column], entry.value.exact(z, ei)});
          }
        } else
        recurrence.column(out.eigen_, out.successor_, column, tops[column],
          static_cast<std::size_t>(d) * (n + 1),
          [&](unsigned k, unsigned row, const IntegerRecurrence::Value &logs,
              const IntegerPolynomial &scale, const IntegerPolynomial &common,
              long pole) {
            fmpz_poly_mul(shared.raw(), scale.raw(), common.raw());
            if (pole > 0)
              fmpz_poly_shift_left(shared.raw(), shared.raw(), pole);
            for (unsigned l = 0; l < logs.size(); ++l) {
              if (logs[l].zero())
                continue;
              // The only gcd of a stored coefficient: lowest terms on output.
              fmpz_poly_shift_left(numerator.raw(), logs[l].raw(),
                                   pole < 0 ? -pole : 0);
              fmpz_poly_gcd(divisor.raw(), numerator.raw(), shared.raw());
              if (fmpz_sgn(fmpz_poly_lead(shared.raw())) < 0)
                fmpz_poly_neg(divisor.raw(), divisor.raw());
              fmpz_poly_div(numerator.raw(), numerator.raw(), divisor.raw());
              fmpz_poly_div(denominator.raw(), shared.raw(), divisor.raw());
              terms.push_back(
                  {row, column, l, out.eigen_[column].power + Rational(k),
                   out.slopes_[column],
                   z.from_reduced_univariate_polynomials(numerator.raw(),
                                                         denominator.raw(), ei)});
            }
            if (terms.size() > options.max_terms)
              throw std::length_error("affine Frobenius term budget exhausted");
          });
        const std::lock_guard lock(shared_state);
        ++completed;
        if(options.column_progress)options.column_progress(completed,d);
      } catch (...) {
        const std::lock_guard lock(shared_state);
        if (!failure) failure = std::current_exception();
        failed = true;
      }
    };
    const unsigned workers = std::min<unsigned>(d,
        options.threads ? options.threads : std::max(1u, std::thread::hardware_concurrency()));
    {
      std::vector<std::jthread> pool;
      for (unsigned i = 1; i < workers; ++i)
        pool.emplace_back([&] { solve_columns(); flint_cleanup(); });
      solve_columns();
    }
    if (failure) std::rethrow_exception(failure);
    std::size_t total = 0;
    for (const auto &terms : columns) total += terms.size();
    if (total > options.max_terms)
      throw std::length_error("affine Frobenius term budget exhausted");
    out.expansion_.terms.reserve(total);
    for (auto &terms : columns)
      std::move(terms.begin(), terms.end(), std::back_inserter(out.expansion_.terms));
    return out;
  }
public:
  const Expansion &terms() const { return expansion_; }
  // A rational change of basis can introduce apparent moving poles even when
  // the physical trace has an epsilon-regular Wronskian. Certify the determinant
  // of G*F from the exact gauge identity, independently of retained x terms.
  std::optional<Exact> gauged_wronskian_prefactor(const Matrix& gauge,
      const Matrix& source_connection,const Matrix& physical_connection)const {
    using namespace affine_frobenius_detail;
    for(const auto* matrix:{&gauge,&source_connection,&physical_connection}){
      if(matrix->size()!=d_)throw std::invalid_argument("Wronskian gauge matrix height");
      for(const auto& row:*matrix)if(row.size()!=d_)throw std::invalid_argument("Wronskian gauge matrix width");
    }
    auto det=determinant(gauge);
    if(det.is_zero())throw std::domain_error("singular Wronskian gauge");
    auto source_trace=zero_,physical_trace=zero_,residue_trace=zero_;
    for(unsigned i=0;i<d_;++i){
      source_trace=source_trace+source_connection[i][i];physical_trace=physical_trace+physical_connection[i][i];
      residue_trace=residue_trace+zero_.constant(eigen_[i].power)+zero_.constant(eigen_[i].slope)*zero_.variable(ei_);
    }
    if(source_trace+det.derivative(xi_)/det!=physical_trace)
      throw std::domain_error("exact Wronskian gauge trace identity failed");
    const auto k=fuchsify::detail::valuation(det,xi_);
    const auto x=zero_.variable(xi_);
    auto unit_trace=physical_trace-(residue_trace+zero_.constant(k))/x;
    std::vector<Exact> origin,at_zero;
    for(unsigned i=0;i<zero_.variable_count();++i){
      origin.push_back(i==xi_||i==ei_?zero_:zero_.variable(i));
      at_zero.push_back(i==xi_?zero_:zero_.variable(i));
    }
    auto rational_unit=zero_.constant(1);
    if(unit_trace.denominator().substitute(origin).is_zero()){
      // Remove simple linear moving poles only when their EXACT residues are
      // constant integers. Epsilon-dependent residues produce log(eps), and
      // higher-order poles can produce essential singularities; reject both.
      if(unit_trace.denominator_terms().size()>4096)return std::nullopt;
      auto factors=polynomial_factors(unit_trace.denominator());
      if(factors.size()>128)return std::nullopt;
      for(const auto& [factor,multiplicity]:factors){
        if(!factor.substitute(origin).is_zero())continue;
        auto derivative=factor.derivative(xi_);
        if(multiplicity!=1||derivative.is_zero()||!derivative.derivative(xi_).is_zero()||factor.substitute(at_zero).is_zero())return std::nullopt;
        auto root=-(factor-derivative*x)/derivative;
        auto at_root=at_zero;at_root[xi_]=root;
        auto residue=(unit_trace*factor/derivative).substitute(at_root);
        if(!residue.is_rational())return std::nullopt;
        const auto exponent=residue.rational();
        if(!is_integer(exponent)||exponent<Rational(-128)||exponent>Rational(128))return std::nullopt;
        const long power=integer_value(exponent);
        rational_unit=rational_unit*fuchsify::detail::power(factor,power);
        if(rational_unit.numerator_terms().size()>options_.max_terms||rational_unit.denominator_terms().size()>options_.max_terms)return std::nullopt;
      }
      unit_trace=unit_trace-rational_unit.derivative(xi_)/rational_unit;
      if(unit_trace.denominator().substitute(origin).is_zero())return std::nullopt;
      auto initial=rational_unit.substitute(at_zero);
      if(initial.is_zero())return std::nullopt;
      rational_unit=rational_unit/initial;
    }
    // det(G)=x^k*g, g(0,eps)!=0. Liouville's identity gives
    // det(GF)=det(frame)*x^(tr(residue)+k)*g(0,eps)*exp(int U).
    // U is jointly analytic, so its exponential is an epsilon unit. Keep x^k
    // in the prefactor because project() retains the source base exponents.
    auto regular=det/fuchsify::detail::power(x,k);
    auto leading=regular.substitute(at_zero);
    if(leading.is_zero())throw std::logic_error("Wronskian gauge has zero leading coefficient");
    return determinant(frame_)*fuchsify::detail::power(x,k)*leading*rational_unit;
  }
  const std::vector<Rational> &absolute_x_frontiers() const {
    return frontiers_;
  }
  std::optional<std::pair<long, B>>
  wronskian_valuation(const Expansion &expansion, const B &point) const {
    if (!expansion.wronskian_prefactor)
      return std::nullopt;
    auto determinant = *expansion.wronskian_prefactor;
    if (determinant.is_zero())
      throw std::domain_error(
          "physical Wronskian prefactor is identically zero");
    long valuation = fuchsify::detail::valuation(determinant, ei_);
    auto regular =
        determinant / fuchsify::detail::power(zero_.variable(ei_), valuation);
    std::vector<Exact> origin;
    for (unsigned i = 0; i < zero_.variable_count(); ++i)
      origin.push_back(i == ei_ ? zero_ : zero_.variable(i));
    auto leading = regular.substitute(origin);
    Jet context(0, 1, B::precision()), x = context;
    x.set(0, point);
    auto value = diffexp::evaluate(data::Reader(leading.str()).read(), context,
                                    {{zero_.variables()[xi_], x}})
                     .at(0);
    return std::pair<long, B>{valuation, value};
  }
  Expansion project(const std::vector<Exact> &row) const {
    return project(Matrix{row});
  }
  // Rational observable rows may have Laurent poles in x. The retained output
  // contains all convolutions through each source column's existing x order,
  // shifted by each row-entry valuation; no omitted-tail certificate is made.
  Expansion project(const Matrix &rows) const {
    return project_impl(rows,false);
  }
private:
  // Shared projection preamble. Every row width is checked before an entry is
  // indexed, and each row records the x valuation that shifts its frontier.
  Expansion projection_shape(const Matrix &rows,bool endpoint_only) const {
    for(const auto& row:rows)if(row.size()!=d_)throw std::invalid_argument("affine observable row dimension");
    Expansion out{static_cast<unsigned>(rows.size()),d_,{}};out.coherent_x_frontier=!endpoint_only;
    if(!endpoint_only&&rows.size()==d_&&expansion_.wronskian_prefactor)
      out.wronskian_prefactor=*expansion_.wronskian_prefactor*affine_frobenius_detail::determinant(rows);
    out.row_x_valuations.assign(rows.size(),Valuations::zero_valuation);
    for(std::size_t i=0;i<rows.size();++i)for(const auto& value:rows[i])if(!value.is_zero())
      out.row_x_valuations[i]=std::min(out.row_x_valuations[i],fuchsify::detail::valuation(value,xi_));
    return out;
  }
  // Every row is projected through its cleared form q(x)*result = sum_j
  // p_j(x)*y_j with polynomial q and p_j: reading polynomial coefficients
  // replaces a Taylor expansion of each rational entry, and dividing by q is a
  // recurrence over deg(q) lags. All arithmetic stays exact in Q(eps).
  Expansion project_impl(const Matrix& rows,bool endpoint_only) const {
    using C=UnivariateRational;using namespace affine_frobenius_detail;
    if(!options_.max_projection_products||rows.size()>options_.max_terms)throw std::length_error("affine projection receiving budget exhausted");
    auto out=projection_shape(rows,endpoint_only);
    const auto x=zero_.variable(xi_);std::size_t products=0;
    auto multiply=[&](const C& a,const C& b){if(products++>=options_.max_projection_products)throw std::length_error("affine projection product budget exhausted");return a*b;};
    auto coefficients=[&](const Exact& p,unsigned top){std::vector<C> result;for(const auto& c:taylor(p,xi_,top))result.emplace_back(c,ei_);trim(result);return result;};
    std::vector<std::vector<std::size_t>> columns(d_);
    for(std::size_t t=0;t<expansion_.terms.size();++t)columns[expansion_.terms[t].column].push_back(t);
    // Each stored coefficient is converted once for the series, not once per
    // projection: an endpoint projects every observable row by itself, and
    // each of those calls read the whole series again.
    std::call_once(sources_->once,[&]{
      sources_->values.reserve(expansion_.terms.size());
      for(const auto& t:expansion_.terms)sources_->values.emplace_back(t.coefficient,ei_);
    });
    const auto& sources=sources_->values;
    std::vector<long> cutoffs(d_);
    for(unsigned i=0;i<rows.size();++i){
      const long valuation=out.row_x_valuations[i];
      if(valuation==Valuations::zero_valuation)continue;
      if(std::abs(valuation)>options_.max_x_order)throw std::length_error("observable Laurent order exceeds affine budget");
      // Retained orders per column; an endpoint needs fixed powers <= 0 only.
      long top=-1;
      for(unsigned c=0;c<d_;++c){
        cutoffs[c]=integer_value(frontiers_[c]-eigen_[c].power);
        if(endpoint_only){
          const auto lowest=eigen_[c].power+Rational(valuation);
          if(!slopes_[c].is_zero()||lowest>Rational(0))cutoffs[c]=-1;
          else while(cutoffs[c]&&lowest+Rational(cutoffs[c])>Rational(0))--cutoffs[c];
        }
        top=std::max(top,cutoffs[c]);
      }
      if(top<0)continue;
      std::vector<Exact> analytic;auto q=zero_.constant(1);
      for(const auto& value:rows[i]){analytic.push_back(value/fuchsify::detail::power(x,valuation));q=q.polynomial_lcm(analytic.back().denominator());}
      auto denominator=coefficients(q,top);if(denominator.empty()||denominator.front().is_zero())throw std::domain_error("affine projection clearing denominator is singular at origin");
      std::vector<std::vector<C>> numerator;for(const auto& value:analytic)numerator.push_back(coefficients(q*value,top));
      for(unsigned c=0;c<d_;++c){
        if(cutoffs[c]<0)continue;
        const unsigned cutoff=cutoffs[c];
        unsigned log_high=0;for(const auto index:columns[c])log_high=std::max(log_high,expansion_.terms[index].log_degree);
        std::vector<std::vector<C>> rhs(log_high+1,std::vector<C>(cutoff+1));
        for(const auto index:columns[c]){
          const auto& t=expansion_.terms[index];
          const auto k=static_cast<unsigned long>(integer_value(t.power-eigen_[c].power));if(k>cutoff||numerator[t.row].empty())continue;
          const auto& value=sources[index];
          for(unsigned lag=0;lag<numerator[t.row].size()&&k+lag<=cutoff;++lag)if(!numerator[t.row][lag].is_zero())rhs[t.log_degree][k+lag]+=multiply(value,numerator[t.row][lag]);
        }
        for(unsigned l=0;l<=log_high;++l)for(unsigned k=0;k<=cutoff;++k){
          auto& value=rhs[l][k];for(unsigned lag=1;lag<denominator.size()&&lag<=k;++lag)if(!denominator[lag].is_zero()&&!rhs[l][k-lag].is_zero())value-=multiply(denominator[lag],rhs[l][k-lag]);
          value=value/denominator[0];
          if(!rhs[l][k].is_zero()){
            if(out.terms.size()>=options_.max_terms)throw std::length_error("affine projection accumulated coordinate budget exhausted");
            out.terms.push_back({i,c,l,eigen_[c].power+Rational(valuation)+Rational(k),slopes_[c],rhs[l][k].exact(zero_,ei_)});
          }
        }
      }
    }
    return normalize(std::move(out));
  }
public:
  // Exact contraction supports pole-bearing matching constants without dropping
  // colliding exponents. Their epsilon cancellations remain visible on
  // expansion.
  Expansion contract(const Expansion &input,
                     const std::vector<Exact> &constants) const {
    if (constants.size() != input.columns)
      throw std::invalid_argument("affine contraction dimension");
    auto out = input;
    out.columns = 1;
    out.wronskian_prefactor.reset();
    // One contracted column mixes several source frontiers.
    out.row_x_valuations.clear();
    for (auto &t : out.terms) {
      if (t.column >= constants.size())
        throw std::invalid_argument(
            "affine expansion coordinate outside shape");
      t.coefficient = t.coefficient * constants[t.column];
      t.column = 0;
    }
    return normalize(std::move(out));
  }
  Expansion primitive(const Expansion &input) const {
    using namespace affine_frobenius_detail;
    Expansion out{input.rows, input.columns, {}};
    for (auto &t : normalize(input).terms) {
      auto a = t.power + Rational(1);
      auto delta =
          zero_.constant(a) + zero_.constant(t.slope) * zero_.variable(ei_);
      if (delta.is_zero())
        throw std::domain_error(
            "unregulated fixed x^-1 logarithm has no meromorphic DR primitive");
      Poly p(t.log_degree + 1, zero_);
      p.back() = t.coefficient;
      p = solve(p, delta);
      for (unsigned l = 0; l < p.size(); ++l)
        if (!p[l].is_zero())
          out.terms.push_back({t.row, t.column, l, a, t.slope, p[l]});
    }
    return normalize(std::move(out));
  }
  // Formal DR sector projection, NOT an ordinary endpoint limit. Exponents
  // with nonzero epsilon slope remain distinct even if their base power is
  // zero.
  Expansion dr_analytic_sector(const Expansion &input) const {
    auto out = normalize(input);
    out.wronskian_prefactor.reset();
    std::erase_if(out.terms, [](const Term &t) {
      return !t.slope.is_zero() || t.log_degree || t.power.sign() < 0 ||
             !affine_frobenius_detail::is_integer(t.power);
    });
    return out;
  }
  Matrix dr_endpoint_constant(const Expansion &input) const {
    auto canonical = normalize(input);
    require_fixed_sector(canonical, false);
    auto out = fuchsify::detail::zeros(input.rows, input.columns, zero_);
    for (auto &t : canonical.terms)
      if (t.slope.is_zero()) {
        if (t.power.sign() < 0 || (t.power.is_zero() && t.log_degree))
          throw std::domain_error(
              "unregulated fixed divergent sector has no DR endpoint constant");
        if (t.power.is_zero() && !t.log_degree)
          out[t.row][t.column] = out[t.row][t.column] + t.coefficient;
      }
    return out;
  }
  struct ZeroConstraint {
    unsigned row;
    Rational power;
    unsigned log_degree;
    std::vector<Exact> coefficients;
  };
  struct DRDomain {
    Expansion admissible;
    std::vector<ZeroConstraint> zero_constraints;
  };
  // The DR functional is defined on a subspace of matching constants when
  // fixed divergent modes occur. Keep that domain explicit: callers must
  // establish every coefficient-row constraint on the physical solution.
  // The strict standalone endpoint/integral methods above/below still reject
  // these modes when no physical boundary has been supplied.
  DRDomain dr_domain(const Expansion& input,bool integral) const {
    auto canonical=normalize(input);require_fixed_sector(canonical,integral);DRDomain out{canonical,{}};
    out.admissible.terms.clear();out.admissible.wronskian_prefactor.reset();
    std::map<std::tuple<unsigned,Rational,unsigned>,std::vector<Exact>> constraints;
    for(const auto& t:canonical.terms) {
      const bool forbidden=t.slope.is_zero() && (integral?t.power<=Rational(-1):
          (t.power.sign()<0 || (t.power.is_zero() && t.log_degree)));
      if(!forbidden){out.admissible.terms.push_back(t);continue;}
      auto [at,inserted]=constraints.try_emplace({t.row,t.power,t.log_degree},input.columns,zero_);
      at->second[t.column]=at->second[t.column]+t.coefficient;
    }
    for(auto& [coordinate,row]:constraints){const auto& [r,p,l]=coordinate;out.zero_constraints.push_back({r,p,l,std::move(row)});}
    return out;
  }
  // A limit needs only fixed nonpositive-power sectors. Positive powers and
  // nonzero epsilon slopes are annihilated by the DR endpoint prescription.
  // Compute exactly the same constant/domain rows without their unused tails.
  DRDomain project_endpoint_domain(const Matrix& rows) const {
    return dr_domain(project_impl(rows,true),false);
  }
  // Validates fixed sectors before applying the meromorphic, termwise lower-end
  // DR prescription. Nonzero-slope sectors are analytically continued; no
  // common convergence domain or omitted-tail bound is inferred by this
  // operation.
  Expansion dr_integral_from_zero(const Expansion &input) const {
    auto canonical = normalize(input);
    require_fixed_sector(canonical, true);
    for (auto &t : canonical.terms)
      if (t.slope.is_zero() && t.power <= Rational(-1))
        throw std::domain_error(
            "unregulated fixed nonintegrable endpoint power in DR integral");
    // Every fixed power is now > -1, so its primitive has a positive power and
    // vanishes at the lower end; nonzero slopes vanish there by continuation.
    return primitive(canonical);
  }
  // Epsilon-adic valuation read from the stored exponents alone; no big
  // coefficient is materialized. A zero coefficient has Valuations::zero_valuation.
  long epsilon_valuation(const Exact &coefficient) const {
    if (coefficient.is_zero())
      return Valuations::zero_valuation;
    const auto low = coefficient.minimum_exponent(ei_),
               pole = coefficient.minimum_exponent(ei_, true);
    if (std::max(low, pole) >
        static_cast<unsigned long>(Valuations::zero_valuation))
      throw std::length_error("endpoint valuation exceeds finite integer range");
    return static_cast<long>(low) - static_cast<long>(pole);
  }
  Valuations valuation_metadata(const Expansion &input) const {
    Valuations out;
    out.minimum_by_column.assign(input.columns, Valuations::zero_valuation);
    for (auto &t : normalize(input).terms) {
      auto val = epsilon_valuation(t.coefficient);
      out.minimum_by_column[t.column] =
          std::min(out.minimum_by_column[t.column], val);
      if (val < 0)
        out.maximum_pole =
            std::max(out.maximum_pole, static_cast<unsigned>(-val));
    }
    return out;
  }
  std::vector<std::vector<B>> evaluate(const Expansion &input, const B &x,
                                       const B &epsilon) const {
    using namespace affine_frobenius_detail;
    if (x.contains_zero())
      throw std::domain_error("affine evaluation requires a nonzero x");
    B logx;
    acb_log(logx.raw(), x.raw(), B::precision());
    auto out = std::vector(input.rows, std::vector<B>(input.columns, B(0)));
    Jet context(0, 1, B::precision()), ep = context;
    ep.set(0, epsilon);
    for (auto &t : normalize(input).terms) {
      auto c = diffexp::evaluate(data::Reader(t.coefficient.str()).read(),
                                  context, {{zero_.variables()[ei_], ep}})
                   .at(0);
      auto exponent = (ball(t.power) + ball(t.slope) * epsilon) * logx;
      B factor;
      acb_exp(factor.raw(), exponent.raw(), B::precision());
      for (unsigned l = 0; l < t.log_degree; ++l)
        factor = factor * logx;
      out[t.row][t.column] = out[t.row][t.column] + c * factor;
    }
    return out;
  }
  LaurentMatrix evaluate_laurent(const Expansion &input, const B &x, int low,
                                 int high) const {
    using namespace affine_frobenius_detail;
    if (x.contains_zero() || low > high ||
        static_cast<long>(high) - low + 1 > options_.max_epsilon_depth)
      throw std::invalid_argument("invalid affine Laurent evaluation window");
    B logx;
    acb_log(logx.raw(), x.raw(), B::precision());
    LaurentMatrix out{
        low, high,
        std::vector(
            input.rows,
            std::vector(input.columns, std::vector<B>(high - low + 1, B(0))))};
    for (auto &t : normalize(input).terms) {
      EpsilonRational coefficient(t.coefficient,ei_);
      long valuation = coefficient.valuation();
      if (valuation > high)
        continue;
      long length = static_cast<long>(high) - valuation + 1;
      if (length > options_.max_epsilon_depth)
        throw std::length_error(
            "affine coefficient poles require additional epsilon-depth budget");
      Jet ep(0, length, B::precision());
      if (length > 1)
        ep.set(1, B(1));
      auto value = coefficient.regular_jet(length);
      auto slope = ep.constant(0);
      slope.set(0, ball(t.slope) * logx);
      value = value * (ep * slope).exp();
      B xp;
      auto exponent = ball(t.power) * logx;
      acb_exp(xp.raw(), exponent.raw(), B::precision());
      for (unsigned l = 0; l < t.log_degree; ++l)
        xp = xp * logx;
      for (long k = std::max<long>(low, valuation); k <= high; ++k)
        out.coefficients[t.row][t.column][k - low] =
            out.coefficients[t.row][t.column][k - low] +
            xp * value.at(k - valuation);
    }
    return out;
  }

private:
  Exact zero_;
  // The coefficients of expansion_ as univariate fractions (project_impl).
  struct Sources { std::once_flag once; std::vector<UnivariateRational> values; };
  std::shared_ptr<Sources> sources_ = std::make_shared<Sources>();
  unsigned d_ = 0, n_ = 0, finite_lag_rows_ = 0;
  std::size_t xi_ = 0, ei_ = 0;
  Options options_;
  Matrix frame_;
  std::vector<Exponent> eigen_;
  // Slope label of each column's terms: the exponent's own, or one for a
  // column regulated by Options::regulate_noninteger_powers.
  std::vector<Rational> slopes_;
  std::vector<Rational> frontiers_;
  std::vector<int> successor_;
  // Structural bound of the logarithmic degree of each column (truncated mode).
  std::vector<int> log_bounds_;
  Expansion expansion_{0, 0, {}};
  explicit AffineFrobeniusSeries(Exact z) : zero_(std::move(z)) {}
  // The fixed (slope-zero) sector is exact only inside the retained frontier:
  // a column known through x^F says nothing about its lattice powers above F.
  // An endpoint limit reads powers <= 0 and an integral rejects powers <= -1,
  // so every such lattice power F+1, F+2, ... must already exceed that range.
  void require_fixed_sector(const Expansion &input, bool integral) const {
    if (input.row_x_valuations.empty())
      return;
    if (input.row_x_valuations.size() != input.rows || input.columns != d_)
      throw std::invalid_argument(
          "affine expansion frontier metadata does not match its shape");
    const Rational bound(integral ? -2 : -1);
    for (unsigned i = 0; i < input.rows; ++i) {
      if (input.row_x_valuations[i] == Valuations::zero_valuation)
        continue;
      for (unsigned c = 0; c < d_; ++c) {
        if (!slopes_[c].is_zero())
          continue;
        const auto complete = frontiers_[c] + Rational(input.row_x_valuations[i]);
        if (!(complete > bound))
          throw std::domain_error(
              "affine x order does not reach the fixed endpoint sector: row " +
              std::to_string(i) + ", column " + std::to_string(c) +
              " is known only through x^(" + complete.str() +
              "); increase the endpoint x order");
      }
    }
  }
  Expansion normalize(Expansion out) const {
    if (out.rows > options_.max_terms || out.columns > d_ ||
        out.terms.size() > options_.max_terms)
      throw std::length_error("affine expansion shape or term budget exceeded: rows="+
        std::to_string(out.rows)+", columns="+std::to_string(out.columns)+", terms="+
        std::to_string(out.terms.size())+", term cap="+std::to_string(options_.max_terms));
    out.omitted_tail_certified = false;
    std::map<std::tuple<unsigned, unsigned, Rational, Rational, unsigned>,
             Exact>
        combined;
    for (auto &t : out.terms) {
      if (t.row >= out.rows || t.column >= out.columns)
        throw std::invalid_argument(
            "affine expansion coordinate outside shape");
      if (t.log_degree >
          static_cast<std::size_t>(options_.max_dimension) *
              (static_cast<std::size_t>(options_.max_x_order) + 1))
        throw std::length_error(
            "affine logarithmic degree exceeds finite budget");
      zero_.require_same_field(t.coefficient);
      if(!t.coefficient.is_univariate(ei_))
        throw std::invalid_argument("affine term coefficient must depend only on epsilon");
      auto key =
          std::make_tuple(t.row, t.column, t.power, t.slope, t.log_degree);
      auto [p, inserted] = combined.try_emplace(key, std::move(t.coefficient));
      if(!inserted)p->second = p->second + t.coefficient;
    }
    out.terms.clear();
    for (auto &[key, c] : combined)
      if (!c.is_zero()) {
        auto [i, j, a, b, l] = key;
        out.terms.push_back({i, j, l, a, b, std::move(c)});
      }
    if (out.terms.size() > options_.max_terms)
      throw std::length_error("affine expansion term budget exhausted");
    return out;
  }
};
} // namespace diffexp
