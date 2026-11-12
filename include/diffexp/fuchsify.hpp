#pragma once
#include "diffexp/exact.hpp"
#include <limits>
#include <optional>

namespace diffexp::fuchsify {
using Matrix = std::vector<std::vector<Exact>>;
struct Options {
  unsigned max_steps = 24, max_dimension = 64, max_shear_candidates = 256,
           max_power = 128;
  // A parameter at whose zero the gauge should stay regular (epsilon). A
  // Jordan frame over the field of the parameters can be singular there
  // although the matrix is not; the system then has poles in the parameter
  // that its solutions do not have. Where such a frame occurs, Jordan chains
  // that are regular with a regular inverse are used instead if they exist;
  // nothing changes where the frame is regular already.
  std::optional<std::size_t> regular_variable;
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
          if (!b[k][j].is_zero())
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
      if (!e.is_zero())
        e = e / q;
    for (unsigned i = 0; i < a.size(); ++i)
      if (i != row) {
        q = a[i][c];
        if (!q.is_zero())
          for (unsigned j = 0; j < a[i].size(); ++j)
            if (!a[row][j].is_zero())
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
      if (!e.is_zero())
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
// Incremental echelon basis of a row span.  Row i is normalised to one in
// column pivots[i] and vanishes in the pivot columns of all earlier rows, so
// reducing a candidate by the rows in insertion order clears every pivot
// column; the candidate lies in the span exactly when nothing is left.
// Span as a plain matrix, re-reduced at every insertion: used where spans are
// small (affine_frobenius.hpp).
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
struct Echelon {
  Matrix rows;
  std::vector<unsigned> pivots;
};
inline bool extend(Echelon &span, std::vector<Exact> v) {
  for (std::size_t i = 0; i < span.rows.size(); ++i) {
    const auto q = v[span.pivots[i]];
    if (q.is_zero())
      continue;
    const auto &row = span.rows[i];
    for (std::size_t j = 0; j < v.size(); ++j)
      if (!row[j].is_zero())
        v[j] = v[j] - q * row[j];
  }
  unsigned p = 0;
  while (p < v.size() && v[p].is_zero())
    ++p;
  if (p == v.size())
    return false;
  const auto q = v[p];
  for (auto &e : v)
    if (!e.is_zero())
      e = e / q;
  span.rows.push_back(std::move(v));
  span.pivots.push_back(p);
  return true;
}
inline std::vector<Exact> action(const Matrix &m, const std::vector<Exact> &v) {
  std::vector<Exact> out(v.size(), v[0].constant(0));
  for (unsigned i = 0; i < v.size(); ++i)
    for (unsigned j = 0; j < v.size(); ++j)
      out[i] = out[i] + m[i][j] * v[j];
  return out;
}
// Echelon basis over the local ring of one parameter at zero: every vector
// vanishes in the pivot rows of the earlier ones and is a unit (valuation
// zero) in its own, so the vectors span a saturated lattice: what they span
// over the field and is regular, they span over the ring.
struct LocalEchelon {
  std::size_t variable;
  Matrix rows;
  std::vector<unsigned> pivots;
  // v minus the combination of the basis that clears its pivot rows, divided
  // by the power of the parameter it then carries; empty inside the span.
  std::vector<Exact> reduce(std::vector<Exact> v) const {
    for (unsigned k = 0; k < rows.size(); ++k) {
      const unsigned p = pivots[k];
      if (v[p].is_zero())
        continue;
      const auto factor = v[p] / rows[k][p];
      for (unsigned j = 0; j < v.size(); ++j)
        if (!rows[k][j].is_zero())
          v[j] = v[j] - factor * rows[k][j];
    }
    std::optional<long> lowest;
    for (const auto &e : v)
      if (!e.is_zero()) {
        const long order = valuation(e, variable);
        if (!lowest || order < *lowest)
          lowest = order;
      }
    if (!lowest)
      return {};
    if (*lowest)
      for (auto &e : v)
        if (!e.is_zero())
          e = e / power(e.variable(variable), *lowest);
    return v;
  }
  void insert(const std::vector<Exact> &v) {
    auto reduced = reduce(v);
    for (unsigned i = 0; i < reduced.size(); ++i)
      if (!reduced[i].is_zero() && valuation(reduced[i], variable) == 0) {
        rows.push_back(std::move(reduced));
        pivots.push_back(i);
        return;
      }
  }
};
// Exact Jordan-chain frame of a nilpotent leading coefficient over
// Q(parameters). With `regular`, the head of every chain is taken regular and
// primitive at the zero of that parameter modulo the kernel below it and the
// chains before it (a head is free up to exactly those), which gives a frame
// that is regular with a regular inverse wherever the Jordan type at the zero
// is the generic one.
inline Matrix nilpotent_frame(const Matrix &m, std::optional<std::size_t> regular = std::nullopt) {
  unsigned d = m.size();
  auto z = m[0][0].constant(0);
  auto pow = identity(d, z);
  // kernels[l] is the kernel of m^l up to the nilpotency index: the first l
  // with m^l = 0.  Every higher power vanishes as well, its kernel is the
  // whole space, and a level above the index starts from a span that is
  // already everything and selects no chain head, so those levels are not
  // formed at all.
  std::vector<Matrix> kernels(1);
  unsigned index = 0;
  for (unsigned l = 1; l <= d && !index; ++l) {
    pow = multiply(pow, m);
    kernels.push_back(nullspace(pow));
    if (zero(pow))
      index = l;
  }
  if (!index)
    throw std::domain_error(
        "irregular endpoint: highest-pole coefficient is not nilpotent");
  Matrix selected;
  std::vector<unsigned> depths; // selected[i] is killed by m^depths[i]
  auto t = zeros(d, d, z), jordan = zeros(d, d, z);
  unsigned column = 0;
  for (unsigned l = index; l > 0; --l) {
    Echelon span;
    for (auto &v : selected)
      extend(span, v);
    for (auto &v : kernels[l - 1])
      extend(span, v);
    LocalEchelon local{regular.value_or(0), {}, {}};
    if (regular) {
      // A head of level l stays in the kernel of m^l: it may be changed by
      // that kernel's vectors only, which of the longer chains are the
      // members at most l steps above their ends.
      for (auto &v : kernels[l - 1])
        local.insert(v);
      for (unsigned i = 0; i < selected.size(); ++i)
        if (depths[i] <= l)
          local.insert(selected[i]);
    }
    for (auto &head : kernels[l])
      if (extend(span, head)) {
        std::vector<std::vector<Exact>> chain(l);
        chain.back() = regular ? local.reduce(head) : head;
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
          depths.push_back(j + 1);
          extend(span, chain[j]);
          if (regular)
            local.insert(chain[j]);
        }
      }
  }
  if (column != d || multiply(m, t) != multiply(t, jordan))
    throw std::logic_error("nilpotent frame certificate failed");
  (void)inverse(t);
  return t;
}
// No entry of the matrix or of its inverse has a pole at the zero of the parameter.
inline bool unimodular(const Matrix &t, std::size_t variable) {
  for (const auto *matrix : {&t})
    for (const auto &row : *matrix)
      for (const auto &e : row)
        if (!e.is_zero() && valuation(e, variable) < 0)
          return false;
  for (const auto &row : inverse(t))
    for (const auto &e : row)
      if (!e.is_zero() && valuation(e, variable) < 0)
        return false;
  return true;
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
// gauge(a, diagonal(x, s), xi) without the matrix products: entry (i,j) is
// multiplied by x^(s[j]-s[i]) and s[i]/x leaves the diagonal. The general
// formula inverts the shear and multiplies three dense matrices, d^3 exact
// operations for a result that has the entries of `a` and no others.
inline Matrix diagonal_gauge(const Matrix &a, const std::vector<long> &s,
                             const Exact &x) {
  auto out = a;
  const unsigned d = a.size();
  for (unsigned i = 0; i < d; ++i) {
    for (unsigned j = 0; j < d; ++j)
      if (s[j] != s[i] && !out[i][j].is_zero())
        out[i][j] = out[i][j] * power(x, s[j] - s[i]);
    if (s[i])
      out[i][i] = out[i][i] - x.constant(Rational(s[i])) / x;
  }
  return out;
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
    // The same step for the shear diagonal(x, s).
    auto apply_shear = [&](const std::vector<long> &s) {
      current = diagonal_gauge(current, s, x);
      for (auto &row : total)
        for (unsigned j = 0; j < d; ++j)
          if (s[j] && !row[j].is_zero())
            row[j] = row[j] * power(x, s[j]);
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
        apply_shear(*shifts);
        continue;
      }
      auto before = score(current, xi);
      auto frame = nilpotent_frame(leading(current, xi, before.first));
      if (options.regular_variable && !unimodular(frame, *options.regular_variable)) {
        auto regular = nilpotent_frame(leading(current, xi, before.first), options.regular_variable);
        if (unimodular(regular, *options.regular_variable))
          frame = std::move(regular);
      }
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
        apply_shear(*shifts);
        continue;
      }
      // Bounded block shears lower the pole/rank pair. Failure is unsupported,
      // not a claim that a nilpotent-leading system is mathematically
      // irregular.
      before = score(current, xi);
      std::optional<std::vector<long>> best;
      auto best_score = before;
      unsigned candidates = 0;
      // The score of a candidate is its pole order and the rank of its
      // leading coefficient. A shear moves the valuation of entry (i,j) by
      // s[j]-s[i] and leaves its leading coefficient alone, so both follow
      // from the valuations and leading coefficients of the current matrix,
      // taken once for all candidates; and a candidate with a higher pole
      // order than the best one needs no rank. Only a candidate of pole
      // order one or less, where the term s[i]/x of the shear enters the
      // leading coefficient, is formed and scored as a matrix.
      std::vector<std::vector<long>> order(d, std::vector<long>(d, 0));
      auto lead = zeros(d, d, x);
      {
        std::vector<Exact> origin;
        for (unsigned i = 0; i < x.variable_count(); ++i)
          origin.push_back(i == xi ? x.constant(0) : x.variable(i));
        for (unsigned i = 0; i < d; ++i)
          for (unsigned j = 0; j < d; ++j)
            if (!current[i][j].is_zero()) {
              order[i][j] = valuation(current[i][j], xi);
              lead[i][j] = (current[i][j] * power(x, -order[i][j])).substitute(origin);
            }
      }
      auto consider = [&](const std::vector<long> &s) {
        if (candidates++ >= options.max_shear_candidates)
          return;
        long pole = 0;
        for (unsigned i = 0; i < d; ++i)
          for (unsigned j = 0; j < d; ++j)
            if (!current[i][j].is_zero())
              pole = std::max(pole, -(order[i][j] + s[j] - s[i]));
        std::pair<unsigned, unsigned> sc;
        if (pole < 2)
          sc = score(diagonal_gauge(current, s, x), xi);
        else {
          if (static_cast<unsigned long>(pole) > best_score.first)
            return;
          auto top = zeros(d, d, x);
          for (unsigned i = 0; i < d; ++i)
            for (unsigned j = 0; j < d; ++j)
              if (!current[i][j].is_zero() && order[i][j] + s[j] - s[i] == -pole)
                top[i][j] = lead[i][j];
          sc = {static_cast<unsigned>(pole), static_cast<unsigned>(rref(top, d).size())};
        }
        if (sc < best_score) {
          best_score = sc;
          best = s;
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
      apply_shear(*best);
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
