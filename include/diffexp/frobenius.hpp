#pragma once
#include "diffexp/exact.hpp"
#include <flint/acb_mat.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <limits>
#include <map>
#include <tuple>

namespace diffexp {
namespace frobenius_detail {
using Q = Rational;
using Matrix = std::vector<std::vector<Q>>;
using Poly = std::vector<Q>;
inline Matrix zeros(unsigned r, unsigned c) {
  return Matrix(r, std::vector<Q>(c));
}
inline Matrix multiply(const Matrix &a, const Matrix &b) {
  auto c = zeros(a.size(), b[0].size());
  for (unsigned i = 0; i < a.size(); ++i)
    for (unsigned k = 0; k < b.size(); ++k)
      if (!a[i][k].is_zero())
        for (unsigned j = 0; j < b[0].size(); ++j)
          c[i][j] += a[i][k] * b[k][j];
  return c;
}
inline std::vector<unsigned> rref(Matrix &a, unsigned columns) {
  std::vector<unsigned> pivots;
  unsigned r = 0;
  for (unsigned c = 0; c < columns && r < a.size(); ++c) {
    unsigned p = r;
    while (p < a.size() && a[p][c].is_zero())
      ++p;
    if (p == a.size())
      continue;
    std::swap(a[p], a[r]);
    auto q = a[r][c];
    for (auto &v : a[r])
      v = v / q;
    for (unsigned i = 0; i < a.size(); ++i)
      if (i != r) {
        q = a[i][c];
        for (unsigned j = 0; j < a[i].size(); ++j)
          a[i][j] -= q * a[r][j];
      }
    pivots.push_back(c);
    ++r;
  }
  return pivots;
}
inline Matrix inverse(Matrix a) {
  unsigned d = a.size();
  for (unsigned i = 0; i < d; ++i) {
    a[i].resize(2 * d);
    a[i][d + i] = Q(1);
  }
  if (rref(a, d).size() != d)
    throw std::domain_error("singular rational frame");
  auto b = zeros(d, d);
  for (unsigned i = 0; i < d; ++i)
    for (unsigned j = 0; j < d; ++j)
      b[i][j] = a[i][d + j];
  return b;
}
inline Q rational(const fmpq_t q) {
  char *s = fmpq_get_str(nullptr, 10, q);
  Q out(s);
  flint_free(s);
  return out;
}
inline std::vector<Q> spectrum(const Matrix &r) {
  unsigned d = r.size();
  fmpq_mat_t m;
  fmpq_mat_init(m, d, d);
  for (unsigned i = 0; i < d; ++i)
    for (unsigned j = 0; j < d; ++j)
      fmpq_set_str(fmpq_mat_entry(m, i, j), r[i][j].str().c_str(), 10);
  fmpq_poly_t p;
  fmpq_poly_init(p);
  fmpq_mat_charpoly(p, m);
  fmpz_poly_t z;
  fmpz_poly_init(z);
  fmpq_poly_get_numerator(z, p);
  fmpz_poly_factor_t f;
  fmpz_poly_factor_init(f);
  fmpz_poly_factor(f, z);
  std::vector<Q> roots;
  bool supported = true;
  fmpq_t q;
  fmpq_init(q);
  for (slong i = 0; i < f->num; ++i) {
    if (fmpz_poly_degree(f->p + i) != 1) {
      supported = false;
      break;
    }
    fmpz_neg(fmpq_numref(q), fmpz_poly_get_coeff_ptr(f->p + i, 0));
    fmpz_set(fmpq_denref(q), fmpz_poly_get_coeff_ptr(f->p + i, 1));
    fmpq_canonicalise(q);
    roots.push_back(rational(q));
  }
  fmpq_clear(q);
  fmpz_poly_factor_clear(f);
  fmpz_poly_clear(z);
  fmpq_poly_clear(p);
  fmpq_mat_clear(m);
  if (!supported)
    throw std::domain_error("Frobenius residue requires rational eigenvalues");
  std::sort(roots.begin(), roots.end());
  return roots;
}
// Nullspace vectors are stored as rows; multiplication elsewhere remains the
// ordinary matrix product. This filtration is finite: powers stop at dimension.
inline Matrix nullspace(Matrix m) {
  const unsigned d = m.size();
  auto pivots = rref(m, d);
  Matrix basis;
  for (unsigned free = 0; free < d; ++free)
    if (std::find(pivots.begin(), pivots.end(), free) == pivots.end()) {
      std::vector<Q> v(d);
      v[free] = Q(1);
      for (unsigned i = 0; i < pivots.size(); ++i)
        v[pivots[i]] = -m[i][free];
      basis.push_back(std::move(v));
    }
  return basis;
}
inline std::vector<Q> action(const Matrix &m, const std::vector<Q> &v) {
  std::vector<Q> result(v.size());
  for (unsigned i = 0; i < v.size(); ++i)
    for (unsigned j = 0; j < v.size(); ++j)
      result[i] += m[i][j] * v[j];
  return result;
}
inline bool extend_span(Matrix &span, const std::vector<Q> &v) {
  auto trial = span;
  trial.push_back(v);
  auto reduced = trial;
  if (rref(reduced, v.size()).size() == span.size())
    return false;
  span = std::move(reduced);
  while (!span.empty() && std::all_of(span.back().begin(), span.back().end(),
                                      [](const Q &q) { return q.is_zero(); }))
    span.pop_back();
  return true;
}
struct JordanFrame {
  Matrix transform;
  std::vector<Q> eigenvalues;
  std::vector<int> successor;
};
inline JordanFrame frame(const Matrix &r) {
  const unsigned d = r.size();
  JordanFrame result{zeros(d, d), {}, {}};
  Matrix selected;
  for (const auto &lambda : spectrum(r)) {
    auto m = r;
    for (unsigned i = 0; i < d; ++i)
      m[i][i] -= lambda;
    auto power = zeros(d, d);
    for (unsigned i = 0; i < d; ++i)
      power[i][i] = Q(1);
    std::vector<Matrix> kernels(d + 1);
    for (unsigned l = 1; l <= d; ++l) {
      power = multiply(power, m);
      kernels[l] = nullspace(power);
    }
    // At level l choose heads modulo ker(M^(l-1)) and already selected chains.
    // Descending lengths prevent shorter chains from consuming a long head.
    for (unsigned l = d; l > 0; --l) {
      Matrix span;
      for (const auto &v : selected)
        extend_span(span, v);
      for (const auto &v : kernels[l - 1])
        extend_span(span, v);
      for (const auto &head : kernels[l])
        if (extend_span(span, head)) {
          std::vector<std::vector<Q>> chain(l);
          chain[l - 1] = head;
          for (unsigned j = l - 1; j > 0; --j)
            chain[j - 1] = action(m, chain[j]);
          for (unsigned j = 0; j < l; ++j) {
            const unsigned column = result.eigenvalues.size();
            if (column >= d)
              throw std::logic_error("Jordan chain count exceeds dimension");
            for (unsigned row = 0; row < d; ++row)
              result.transform[row][column] = chain[j][row];
            result.eigenvalues.push_back(lambda);
            result.successor.push_back(j + 1 < l ? static_cast<int>(column + 1)
                                                 : -1);
            selected.push_back(chain[j]);
            extend_span(span, chain[j]);
          }
        }
    }
  }
  if (result.eigenvalues.size() != d)
    throw std::logic_error("incomplete exact Jordan frame");
  // Both assertions are exact certificates, independent of construction logic.
  (void)inverse(result.transform);
  auto jordan = zeros(d, d);
  for (unsigned i = 0; i < d; ++i) {
    jordan[i][i] = result.eigenvalues[i];
    if (result.successor[i] >= 0)
      jordan[i][result.successor[i]] = Q(1);
  }
  if (multiply(r, result.transform) != multiply(result.transform, jordan))
    throw std::logic_error("exact Jordan frame fails R*T=T*J");
  return result;
}
inline void trim(Poly &p) {
  while (!p.empty() && p.back().is_zero())
    p.pop_back();
}
inline void add_scaled(Poly &out, const Poly &p, const Q &q) {
  if (q.is_zero())
    return;
  if (out.size() < p.size())
    out.resize(p.size());
  for (unsigned l = 0; l < p.size(); ++l)
    out[l] += q * p[l];
  trim(out);
}
// Unique polynomial particular solution, with zero integration constant at
// resonance.
inline Poly solve(const Poly &rhs, const Q &a) {
  if (rhs.empty())
    return {};
  Poly out(rhs.size() + (a.is_zero() ? 1 : 0));
  if (a.is_zero()) {
    for (unsigned l = 0; l < rhs.size(); ++l)
      out[l + 1] = rhs[l] / Q(l + 1);
  } else {
    for (unsigned l = rhs.size(); l-- > 0;)
      out[l] =
          (rhs[l] - (l + 1 < out.size() ? Q(l + 1) * out[l + 1] : Q(0))) / a;
  }
  trim(out);
  return out;
}

// Rectangular, jointly analytic bivariate expansion. All other parameters must
// already have been substituted exactly.
inline std::vector<Q> taylor(const Exact &e, std::size_t xi, std::size_t ei,
                             unsigned n, unsigned k) {
  if (xi == ei || xi >= e.variable_count() || ei >= e.variable_count())
    throw std::invalid_argument("Taylor variable indices");
  std::vector<Q> num((n + 1) * (k + 1)), den(num.size()), out(num.size());
  auto collect = [&](const auto &terms, auto &v) {
    for (const auto &t : terms) {
      for (unsigned z = 0; z < t.powers.size(); ++z)
        if (z != xi && z != ei && t.powers[z])
          throw std::invalid_argument("unsubstituted Frobenius parameter");
      auto nx = t.powers[xi], ke = t.powers[ei];
      if (nx <= n && ke <= k)
        v[nx * (k + 1) + ke] += t.coefficient;
    }
  };
  collect(e.numerator_terms(), num);
  collect(e.denominator_terms(), den);
  if (den[0].is_zero())
    throw std::domain_error(
        "rational expression is not jointly analytic at (x,epsilon)=(0,0)");
  for (unsigned nx = 0; nx <= n; ++nx)
    for (unsigned ke = 0; ke <= k; ++ke) {
      auto v = num[nx * (k + 1) + ke];
      for (unsigned dx = 0; dx <= nx; ++dx)
        for (unsigned de = 0; de <= ke; ++de)
          if (dx || de)
            v -= den[dx * (k + 1) + de] * out[(nx - dx) * (k + 1) + ke - de];
      out[nx * (k + 1) + ke] = v / den[0];
    }
  return out;
}
inline kernel::ComplexBall ball(const Q &q) {
  kernel::ComplexBall b;
  fmpq_t v;
  fmpq_init(v);
  fmpq_set_str(v, q.str().c_str(), 10);
  acb_set_fmpq(b.raw(), v, kernel::ComplexBall::precision());
  fmpq_clear(v);
  return b;
}
} // namespace frobenius_detail

// Formal epsilon expansion about a residue with rational spectrum at epsilon=0.
// Acb radii enclose retained arithmetic only; no omitted-x-tail bound is
// claimed.
class FrobeniusSeries {
public:
  using B = kernel::ComplexBall;
  using Boundary = std::vector<std::vector<B>>;
  using MatrixSeries =
      std::vector<std::vector<std::vector<B>>>; // row,column,epsilon
  using ExactMatrix = std::vector<std::vector<Exact>>;
  using Poly = frobenius_detail::Poly;
  struct Monomial {
    unsigned row, column, epsilon, log_degree;
    Rational power, coefficient;
  };
  static constexpr bool omitted_tail_certified = false;
  unsigned dimension() const { return d_; }
  unsigned x_order() const { return n_; }
  unsigned epsilon_order() const { return k_; }
  const std::vector<Rational> &exponents() const { return eigen_; }
  static FrobeniusSeries prepare(const ExactMatrix &a, std::size_t xi,
                                 std::size_t ei, unsigned n, unsigned k) {
    using namespace frobenius_detail;
    if (a.empty() || xi == ei || n == std::numeric_limits<unsigned>::max() ||
        k == std::numeric_limits<unsigned>::max())
      throw std::invalid_argument("invalid Frobenius matrix or variables");
    // Check products before allocating. A causal recurrence path has at most
    // n+k lower layers, each of which can add at most d logarithms; thus the
    // exact log ladder also has a finite structural bound, without chopping.
    const std::size_t dsize = a.size();
    const std::size_t nsize = static_cast<std::size_t>(n) + 1;
    const std::size_t ksize = static_cast<std::size_t>(k) + 1;
    const auto max = std::numeric_limits<std::size_t>::max();
    if (dsize > std::numeric_limits<int>::max() || nsize > max / ksize ||
        nsize * ksize > std::numeric_limits<unsigned>::max() / dsize ||
        dsize > max / dsize || dsize * dsize > max / (nsize * ksize) ||
        dsize > max / (static_cast<std::size_t>(n) + k + 1))
      throw std::length_error(
          "Frobenius dimensions exceed addressable finite work");
    const std::size_t log_slots = dsize * (static_cast<std::size_t>(n) + k + 1);
    FrobeniusSeries out;
    out.d_ = a.size();
    out.n_ = n;
    out.k_ = k;
    unsigned d = out.d_;
    // Expand x A as a bivariate rational Taylor series. A denominator with a
    // zero constant is rejected, rather than silently choosing an iterated
    // limit.
    std::vector<Matrix> g(nsize * ksize, zeros(d, d));
    for (unsigned i = 0; i < d; ++i) {
      if (a[i].size() != d)
        throw std::invalid_argument("Frobenius matrix must be square");
      for (unsigned j = 0; j < d; ++j) {
        auto values = taylor(a[i][j] * a[i][j].variable(xi), xi, ei, n, k);
        for (unsigned nx = 0; nx <= n; ++nx)
          for (unsigned ke = 0; ke <= k; ++ke)
            g[nx * (k + 1) + ke][i][j] = values[nx * (k + 1) + ke];
      }
    }
    auto [t, eigen, successor] = frame(g[0]);
    out.eigen_ = eigen;
    auto ti = inverse(t);
    for (auto &m : g)
      m = multiply(multiply(ti, m), t);
    out.coefficients_.resize(dsize * dsize * nsize * ksize);
    for (unsigned column = 0; column < d; ++column) {
      std::vector<Poly> v(nsize * ksize * dsize);
      auto at = [&](unsigned nx, unsigned ke, unsigned row) -> Poly & {
        return v[(nx * (k + 1) + ke) * d + row];
      };
      for (unsigned nx = 0; nx <= n; ++nx)
        for (unsigned ke = 0; ke <= k; ++ke)
          for (unsigned row = d; row-- > 0;) {
            if (nx == 0 && ke == 0) {
              Poly rhs;
              if (successor[row] >= 0)
                rhs = at(nx, ke, successor[row]);
              at(nx, ke, row) = solve(rhs, eigen[column] - eigen[row]);
              if (row == column)
                add_scaled(at(nx, ke, row), Poly{Q(1)}, Q(1));
              continue;
            }
            Poly rhs;
            for (unsigned dx = 0; dx <= nx; ++dx)
              for (unsigned de = 0; de <= ke; ++de)
                if (dx || de)
                  for (unsigned j = 0; j < d; ++j)
                    add_scaled(rhs, at(nx - dx, ke - de, j),
                               g[dx * (k + 1) + de][row][j]);
            // Reverse chain order solves (D+delta-N)v=rhs. At delta=0
            // each row is integrated with zero particular constant; otherwise
            // scalar polynomial inversion is the finite (delta I-N)^-1 solve.
            if (successor[row] >= 0)
              add_scaled(rhs, at(nx, ke, successor[row]), Q(1));
            at(nx, ke, row) = solve(rhs, eigen[column] + Q(nx) - eigen[row]);
            if (at(nx, ke, row).size() > log_slots)
              throw std::logic_error(
                  "Frobenius log ladder exceeds causal bound");
          }
      for (unsigned nx = 0; nx <= n; ++nx)
        for (unsigned ke = 0; ke <= k; ++ke)
          for (unsigned row = 0; row < d; ++row)
            for (unsigned j = 0; j < d; ++j)
              add_scaled(out.at(row, column, nx, ke), at(nx, ke, j), t[row][j]);
    }
    return out;
  }
  std::vector<Monomial> monomials() const {
    std::vector<Monomial> out;
    for (unsigned i = 0; i < d_; ++i)
      for (unsigned c = 0; c < d_; ++c)
        for (unsigned n = 0; n <= n_; ++n)
          for (unsigned k = 0; k <= k_; ++k) {
            const auto &p = at(i, c, n, k);
            for (unsigned l = 0; l < p.size(); ++l)
              if (!p[l].is_zero())
                out.push_back({i, c, k, l, eigen_[c] + Rational(n), p[l]});
          }
    return out;
  }
  // Project an analytic rational observable row BEFORE numerical matching, so
  // exact cancellations of nonintegrable powers are preserved. The same Taylor
  // cutoff is applied after convolution, relative to each column exponent.
  std::vector<Monomial> project(const std::vector<Exact> &row, std::size_t xi,
                                std::size_t ei) const {
    using namespace frobenius_detail;
    if (row.size() != d_)
      throw std::invalid_argument("Frobenius observable row dimension");
    std::vector<std::vector<Q>> weights;
    for (const auto &e : row)
      weights.push_back(taylor(e, xi, ei, n_, k_));
    std::map<std::tuple<unsigned, unsigned, Q, unsigned>, Q> combined;
    for (const auto &m : monomials())
      for (unsigned n = 0; n <= n_; ++n) {
        auto power = m.power + Q(n);
        if (power > eigen_[m.column] + Q(n_))
          break;
        for (unsigned e = 0; e + m.epsilon <= k_; ++e) {
          const auto &w = weights[m.row][n * (k_ + 1) + e];
          if (!w.is_zero())
            combined[{m.column, m.epsilon + e, power, m.log_degree}] +=
                m.coefficient * w;
        }
      }
    std::vector<Monomial> result;
    for (const auto &[key, q] : combined)
      if (!q.is_zero()) {
        auto [c, e, p, l] = key;
        result.push_back({0, c, e, l, p, q});
      }
    return result;
  }
  // Integral from zero of an exactly projected, locally integrable finite
  // series. Rejects every retained power <= -1; no DR/finite-part shortcut.
  // Output balls contain arithmetic error only, not the omitted Taylor tail.
  std::vector<B> integrate_projected(const std::vector<Monomial> &terms,
                                     const B &x,
                                     const Boundary &constants) const {
    using namespace frobenius_detail;
    validate(constants);
    if (x.contains_zero())
      throw std::domain_error(
          "projected integration requires nonzero endpoint");
    B logx;
    acb_log(logx.raw(), x.raw(), B::precision());
    std::vector<B> out(k_ + 1, B(0));
    for (const auto &m : terms) {
      if (m.column >= d_ || m.epsilon > k_)
        throw std::invalid_argument("projected monomial context mismatch");
      if (m.coefficient.is_zero())
        continue;
      const auto power = m.power + Q(1);
      if (power.sign() <= 0)
        throw std::domain_error(
            "projected series has a nonintegrable endpoint power");
      Poly p(m.log_degree + 1);
      p.back() = m.coefficient;
      p = solve(p, power);
      B value(0);
      for (unsigned l = p.size(); l-- > 0;)
        value = value * logx + ball(p[l]);
      B xp;
      auto exponent = ball(power) * logx;
      acb_exp(xp.raw(), exponent.raw(), B::precision());
      value = value * xp;
      for (unsigned e = m.epsilon; e <= k_; ++e)
        out[e] = out[e] + value * constants[m.column][e - m.epsilon];
    }
    return out;
  }
  MatrixSeries evaluate(const B &x) const {
    return evaluate_impl(x, Rational(0), false);
  }
  Boundary solution(const B &x, const Boundary &constants) const {
    return apply(evaluate(x), constants);
  }
  Boundary primitive(const B &x, const Boundary &constants,
                     const Rational &weight_power = Rational(0)) const {
    return apply(evaluate_impl(x, weight_power, true), constants);
  }
  // Hadamard finite part of this canonical primitive at zero. All nonzero
  // powers and logarithms are discarded. It is NOT a convergent endpoint claim.
  Boundary
  finite_part_at_zero(const Boundary &constants,
                      const Rational &weight_power = Rational(0)) const {
    (void)weight_power;
    validate(constants);
    return Boundary(d_, std::vector<B>(k_ + 1, B(0)));
  }
  Boundary match(const B &x, const Boundary &boundary) const {
    validate(boundary);
    auto f = evaluate(x);
    acb_mat_t a, inv;
    acb_mat_init(a, d_, d_);
    acb_mat_init(inv, d_, d_);
    for (unsigned i = 0; i < d_; ++i)
      for (unsigned j = 0; j < d_; ++j)
        acb_set(acb_mat_entry(a, i, j), f[i][j][0].raw());
    int ok = acb_mat_inv(inv, a, B::precision());
    if (!ok) {
      acb_mat_clear(a);
      acb_mat_clear(inv);
      throw std::domain_error(
          "Frobenius matching matrix not invertible at working precision");
    }
    Boundary c(d_, std::vector<B>(k_ + 1, B(0)));
    for (unsigned k = 0; k <= k_; ++k) {
      std::vector<B> rhs(d_, B(0));
      for (unsigned i = 0; i < d_; ++i) {
        rhs[i] = boundary[i][k];
        for (unsigned j = 0; j < d_; ++j)
          for (unsigned e = 1; e <= k; ++e)
            rhs[i] = rhs[i] - f[i][j][e] * c[j][k - e];
      }
      for (unsigned i = 0; i < d_; ++i)
        for (unsigned j = 0; j < d_; ++j)
          acb_addmul(c[i][k].raw(), acb_mat_entry(inv, i, j), rhs[j].raw(),
                     B::precision());
    }
    acb_mat_clear(a);
    acb_mat_clear(inv);
    return c;
  }

private:
  unsigned d_ = 0, n_ = 0, k_ = 0;
  std::vector<Rational> eigen_;
  std::vector<Poly> coefficients_;
  Poly &at(unsigned i, unsigned c, unsigned n, unsigned k) {
    return coefficients_[((static_cast<std::size_t>(i) * d_ + c) * (n_ + 1) +
                          n) *
                             (k_ + 1) +
                         k];
  }
  const Poly &at(unsigned i, unsigned c, unsigned n, unsigned k) const {
    return coefficients_[((static_cast<std::size_t>(i) * d_ + c) * (n_ + 1) +
                          n) *
                             (k_ + 1) +
                         k];
  }
  void validate(const Boundary &b) const {
    if (b.size() != d_)
      throw std::invalid_argument("Frobenius boundary dimension");
    for (auto &r : b)
      if (r.size() != k_ + 1)
        throw std::invalid_argument("Frobenius boundary epsilon window");
  }
  Boundary apply(const MatrixSeries &f, const Boundary &c) const {
    validate(c);
    Boundary b(d_, std::vector<B>(k_ + 1, B(0)));
    for (unsigned i = 0; i < d_; ++i)
      for (unsigned j = 0; j < d_; ++j)
        for (unsigned k = 0; k <= k_; ++k)
          for (unsigned e = 0; e <= k; ++e)
            b[i][k] = b[i][k] + f[i][j][e] * c[j][k - e];
    return b;
  }
  MatrixSeries evaluate_impl(const B &x, const Rational &shift,
                             bool integrate) const {
    using namespace frobenius_detail;
    if (x.contains_zero())
      throw std::domain_error("Frobenius evaluation requires nonzero x");
    B logx;
    acb_log(logx.raw(), x.raw(), B::precision());
    MatrixSeries f(
        d_, std::vector<std::vector<B>>(d_, std::vector<B>(k_ + 1, B(0))));
    for (unsigned c = 0; c < d_; ++c)
      for (unsigned n = 0; n <= n_; ++n) {
        auto power = eigen_[c] + Rational(n) +
                     (integrate ? shift + Rational(1) : Rational(0));
        B xp;
        auto exponent = ball(power) * logx;
        acb_exp(xp.raw(), exponent.raw(), B::precision());
        for (unsigned i = 0; i < d_; ++i)
          for (unsigned k = 0; k <= k_; ++k) {
            auto p = integrate ? solve(at(i, c, n, k), power) : at(i, c, n, k);
            B sum(0);
            for (unsigned l = p.size(); l-- > 0;)
              sum = sum * logx + ball(p[l]);
            f[i][c][k] = f[i][c][k] + xp * sum;
          }
      }
    return f;
  }
};
} // namespace diffexp
