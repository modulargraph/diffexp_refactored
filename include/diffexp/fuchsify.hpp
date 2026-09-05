#pragma once
#include "diffexp/exact.hpp"
#include <limits>
#include <optional>

namespace diffexp::fuchsify {
using Matrix = std::vector<std::vector<Exact>>;
struct Options {
  unsigned max_steps = 24, max_dimension = 64, max_shear_candidates = 256,
           max_power = 128;
};
namespace detail {
inline Matrix zeros(unsigned r, unsigned c, const Exact &z) {
  return Matrix(r, std::vector<Exact>(c, z.constant(0)));
}
inline Matrix identity(unsigned d, const Exact &z) {
  auto m = zeros(d, d, z);
  for (unsigned i = 0; i < d; ++i)
    m[i][i] = z.constant(1);
  return m;
}
inline Matrix multiply(const Matrix &a, const Matrix &b) {
  auto c = zeros(a.size(), b[0].size(), a[0][0]);
  for (unsigned i = 0; i < a.size(); ++i)
    for (unsigned k = 0; k < b.size(); ++k)
      if (!a[i][k].is_zero())
        for (unsigned j = 0; j < b[0].size(); ++j)
          c[i][j] = c[i][j] + a[i][k] * b[k][j];
  return c;
}
inline Matrix subtract(Matrix a, const Matrix &b) {
  for (unsigned i = 0; i < a.size(); ++i)
    for (unsigned j = 0; j < a[i].size(); ++j)
      a[i][j] = a[i][j] - b[i][j];
  return a;
}
inline Matrix derivative(Matrix a, std::size_t xi) {
  for (auto &row : a)
    for (auto &e : row)
      e = e.derivative(xi);
  return a;
}
inline bool zero(const Matrix &a) {
  for (auto &row : a)
    for (auto &e : row)
      if (!e.is_zero())
        return false;
  return true;
}
inline std::vector<unsigned> rref(Matrix &a, unsigned columns) {
  std::vector<unsigned> pivots;
  unsigned row = 0;
  for (unsigned c = 0; c < columns && row < a.size(); ++c) {
    unsigned p = row;
    while (p < a.size() && a[p][c].is_zero())
      ++p;
    if (p == a.size())
      continue;
    std::swap(a[p], a[row]);
    auto q = a[row][c];
    for (auto &e : a[row])
      e = e / q;
    for (unsigned i = 0; i < a.size(); ++i)
      if (i != row) {
        q = a[i][c];
        if (!q.is_zero())
          for (unsigned j = 0; j < a[i].size(); ++j)
            a[i][j] = a[i][j] - q * a[row][j];
      }
    pivots.push_back(c);
    ++row;
  }
  return pivots;
}
inline Matrix inverse(Matrix a) {
  unsigned d = a.size();
  const auto z = a[0][0].constant(0);
  for (unsigned i = 0; i < d; ++i) {
    a[i].resize(2 * d, z);
    a[i][d + i] = z.constant(1);
  }
  if (rref(a, d).size() != d)
    throw std::domain_error("Fuchsification gauge is singular");
  auto inv = zeros(d, d, z);
  for (unsigned i = 0; i < d; ++i)
    for (unsigned j = 0; j < d; ++j)
      inv[i][j] = a[i][d + j];
  return inv;
}
inline long valuation(const Exact &e, std::size_t xi) {
  if (e.is_zero())
    return std::numeric_limits<long>::max() / 4;
  auto order = [&](const auto &terms) {
    unsigned long v = std::numeric_limits<unsigned long>::max();
    for (const auto &t : terms)
      v = std::min(v, t.powers.at(xi));
    if (v > static_cast<unsigned long>(std::numeric_limits<long>::max() / 4))
      throw std::length_error(
          "endpoint valuation exceeds finite integer range");
    return static_cast<long>(v);
  };
  return order(e.numerator_terms()) - order(e.denominator_terms());
}
inline unsigned pole_order(const Matrix &a, std::size_t xi) {
  long p = 0;
  for (auto &row : a)
    for (auto &e : row)
      p = std::max(p, -valuation(e, xi));
  if (p > std::numeric_limits<unsigned>::max())
    throw std::length_error("endpoint pole order exceeds finite range");
  return p;
}
inline Exact power(const Exact &x, long n) {
  return n >= 0 ? x.pow(n) : x.constant(1) / x.pow(-n);
}
inline Matrix leading(const Matrix &a, std::size_t xi, unsigned p) {
  auto out = a;
  auto z = a[0][0].constant(0), x = z.variable(xi);
  std::vector<Exact> replacement;
  for (unsigned i = 0; i < z.variable_count(); ++i)
    replacement.push_back(i == xi ? z : z.variable(i));
  for (auto &row : out)
    for (auto &e : row)
      e = (e * x.pow(p)).substitute(replacement);
  return out;
}
inline std::pair<unsigned, unsigned> score(const Matrix &a, std::size_t xi) {
  unsigned p = pole_order(a, xi);
  auto l = leading(a, xi, p);
  return {p, rref(l, l.size()).size()};
}
inline Matrix nullspace(Matrix m) {
  unsigned d = m.size();
  auto p = rref(m, d);
  Matrix out;
  auto z = m[0][0].constant(0);
  for (unsigned f = 0; f < d; ++f)
    if (std::find(p.begin(), p.end(), f) == p.end()) {
      std::vector<Exact> v(d, z);
      v[f] = z.constant(1);
      for (unsigned i = 0; i < p.size(); ++i)
        v[p[i]] = -m[i][f];
      out.push_back(v);
    }
  return out;
}
inline bool extend(Matrix &span, const std::vector<Exact> &v) {
  auto next = span;
  next.push_back(v);
  auto rank = rref(next, v.size()).size();
  if (rank == span.size())
    return false;
  next.resize(rank);
  span = std::move(next);
  return true;
}
inline std::vector<Exact> action(const Matrix &m, const std::vector<Exact> &v) {
  std::vector<Exact> out(v.size(), v[0].constant(0));
  for (unsigned i = 0; i < v.size(); ++i)
    for (unsigned j = 0; j < v.size(); ++j)
      out[i] = out[i] + m[i][j] * v[j];
  return out;
}
// Exact Jordan-chain frame of a nilpotent leading coefficient over
// Q(parameters).
inline Matrix nilpotent_frame(const Matrix &m) {
  unsigned d = m.size();
  auto z = m[0][0].constant(0);
  auto pow = identity(d, z);
  std::vector<Matrix> kernels(d + 1);
  for (unsigned l = 1; l <= d; ++l) {
    pow = multiply(pow, m);
    kernels[l] = nullspace(pow);
  }
  if (!zero(pow))
    throw std::domain_error(
        "irregular endpoint: highest-pole coefficient is not nilpotent");
  Matrix selected;
  auto t = zeros(d, d, z), jordan = zeros(d, d, z);
  unsigned column = 0;
  for (unsigned l = d; l > 0; --l) {
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
          if (column >= d)
            throw std::logic_error("nilpotent frame overflow");
          for (unsigned i = 0; i < d; ++i)
            t[i][column] = chain[j][i];
          if (j + 1 < l)
            jordan[column][column + 1] = z.constant(1);
          ++column;
          selected.push_back(chain[j]);
          extend(span, chain[j]);
        }
      }
  }
  if (column != d || multiply(m, t) != multiply(t, jordan))
    throw std::logic_error("nilpotent frame certificate failed");
  (void)inverse(t);
  return t;
}
inline Matrix gauge(const Matrix &a, const Matrix &t, std::size_t xi) {
  return multiply(inverse(t), subtract(multiply(a, t), derivative(t, xi)));
}
inline Matrix diagonal(const Exact &x, const std::vector<long> &s) {
  auto t = zeros(s.size(), s.size(), x);
  for (unsigned i = 0; i < s.size(); ++i)
    t[i][i] = power(x, s[i]);
  return t;
}
// Difference constraints: ord(Aij)+s[j]-s[i]>=-1. Bellman-Ford either
// constructs a complete monomial gauge or proves this coordinate ansatz fails.
inline std::optional<std::vector<long>> valuations(const Matrix &a,
                                                   std::size_t xi) {
  unsigned d = a.size();
  std::vector<long> s(d, 0);
  for (unsigned pass = 0; pass < d; ++pass) {
    bool changed = false;
    for (unsigned i = 0; i < d; ++i)
      for (unsigned j = 0; j < d; ++j)
        if (!a[i][j].is_zero()) {
          long bound = s[j] + valuation(a[i][j], xi) + 1;
          if (s[i] > bound) {
            s[i] = bound;
            changed = true;
          }
        }
    if (!changed)
      return s;
  }
  return std::nullopt;
}
} // namespace detail
struct Result {
  bool success = false, identity_verified = false;
  std::string reason;
  Matrix matrix, transform, inverse_transform;
  unsigned initial_pole_order = 0, final_pole_order = 0, steps = 0;
  // Y=T Z, so every physical observable row transforms by right multiplication.
  Matrix transform_rows(const Matrix &rows) const {
    if (!success)
      throw std::domain_error("Fuchsification was unsuccessful");
    if (rows.empty())
      return {};
    for (const auto &r : rows)
      if (r.size() != transform.size())
        throw std::invalid_argument("observable gauge dimension");
    return detail::multiply(rows, transform);
  }
};
inline unsigned pole_order(const Matrix &a, std::size_t xi) {
  return detail::pole_order(a, xi);
}
inline Result prepare(const Matrix &input, std::size_t xi,
                      const Options &options = {}) {
  using namespace detail;
  Result result;
  try {
    if (input.empty() || input.size() > options.max_dimension ||
        !options.max_steps || !options.max_shear_candidates ||
        !options.max_power)
      throw std::invalid_argument("Fuchsification requires a nonempty matrix "
                                  "and positive finite budgets");
    unsigned d = input.size();
    for (const auto &row : input) {
      if (row.size() != d)
        throw std::invalid_argument("Fuchsification matrix must be square");
      for (const auto &e : row)
        if (xi >= e.variable_count())
          throw std::invalid_argument("Fuchsification variable index");
    }
    auto current = input;
    auto x = input[0][0].variable(xi);
    auto total = identity(d, x);
    result.initial_pole_order = detail::pole_order(input, xi);
    if (result.initial_pole_order > options.max_power)
      throw std::length_error(
          "endpoint pole order exceeds Fuchsification power budget");
    auto apply = [&](const Matrix &t) {
      current = gauge(current, t, xi);
      total = multiply(total, t);
      ++result.steps;
    };
    while (detail::pole_order(current, xi) > 1) {
      if (result.steps >= options.max_steps)
        throw std::runtime_error(
            "Fuchsification step budget exhausted; provide a better rational "
            "master basis or increase max_steps");
      if (auto shifts = valuations(current, xi)) {
        for (long s : *shifts)
          if (s > options.max_power ||
              s < -static_cast<long>(options.max_power))
            throw std::length_error(
                "Fuchsification shear exceeds power budget");
        apply(diagonal(x, *shifts));
        continue;
      }
      auto before = score(current, xi);
      auto frame = nilpotent_frame(leading(current, xi, before.first));
      if (frame != identity(d, x)) {
        apply(frame);
        if (result.steps >= options.max_steps)
          continue;
      }
      if (auto shifts = valuations(current, xi)) {
        for (long s : *shifts)
          if (s > options.max_power ||
              s < -static_cast<long>(options.max_power))
            throw std::length_error(
                "Fuchsification shear exceeds power budget");
        apply(diagonal(x, *shifts));
        continue;
      }
      // Bounded block shears lower the pole/rank pair. Failure is unsupported,
      // not a claim that a nilpotent-leading system is mathematically
      // irregular.
      before = score(current, xi);
      std::optional<Matrix> best;
      auto best_score = before;
      unsigned candidates = 0;
      auto consider = [&](const std::vector<long> &s) {
        if (candidates++ >= options.max_shear_candidates)
          return;
        auto t = diagonal(x, s);
        auto candidate = gauge(current, t, xi);
        auto sc = score(candidate, xi);
        if (sc < best_score) {
          best_score = sc;
          best = std::move(t);
        }
      };
      for (unsigned cut = 1;
           cut < d && candidates < options.max_shear_candidates; ++cut)
        for (long direction : {-1L, 1L}) {
          std::vector<long> s(d, 0);
          for (unsigned j = cut; j < d; ++j)
            s[j] = direction;
          consider(s);
        }
      for (unsigned j = 0; j < d && candidates < options.max_shear_candidates;
           ++j)
        for (long direction : {-1L, 1L}) {
          std::vector<long> s(d, 0);
          s[j] = direction;
          consider(s);
        }
      if (!best)
        throw std::runtime_error(
            "unsupported endpoint gauge: bounded valuation/Jordan block shears "
            "did not reduce pole rank; a general Moser reduction or supplied "
            "rational gauge is required");
      apply(*best);
    }
    result.matrix = current;
    result.transform = total;
    result.inverse_transform = inverse(total);
    result.final_pole_order = detail::pole_order(current, xi);
    result.identity_verified =
        subtract(multiply(input, total), derivative(total, xi)) ==
        multiply(total, current);
    if (!result.identity_verified ||
        multiply(total, result.inverse_transform) != identity(d, x))
      throw std::logic_error("exact rational gauge identity failed");
    result.success = true;
  } catch (const std::exception &e) {
    result.reason = e.what();
  }
  return result;
}
} // namespace diffexp::fuchsify
