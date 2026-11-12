#pragma once
// Removal of moving poles from the connection of a master basis.
//
// A denominator factor f(x,eps) of A that mixes the parameter with epsilon is
// a singular point whose position depends on the dimension. The integrals
// have none (Landau singularities do not move with d), and in every case met
// so far the pole is apparent: all solutions of J'=A J are regular there. It
// says that the masters satisfy a linear relation at f=0, i.e. that they stop
// being a basis at that point. The closure below ends exactly in that case.
//
// In terms of rows w (standing for the function w.J) with D w=w'+w A: the
// rows regular at f form a lattice M over the local ring at f, and A has a
// pole at f exactly when M is not closed under D. A basis whose connection is
// regular at f spans a D-closed lattice, and there is only one that contains
// M (two such bases differ by a ratio of fundamental matrices that are both
// regular and invertible at f). It is finite over M when the singularity is
// apparent, and is reached one generator at a time. If row i of A has a pole
// of order p at f, then f^(p-1) D e_i lies in the closure, and so does
//
//   q = r/f,   r = (f^p A_i) reduced modulo f,
//
// because the two differ by a row regular at f. Replacing one master J_k with
// r_k a unit by q.J enlarges the lattice by one step; the old masters stay
// polynomial combinations of the new ones (J_k = (f q.J - sum_j r_j J_j)/r_k),
// so observables acquire no pole at f. For a factor linear in x the reduction
// modulo f is the substitution x=x0(eps) and r depends on epsilon alone. In
// general r has coefficients in the residue field at f; scaled to r_k=1 it
// gives J_k back as a polynomial combination. No new singular point in x
// appears either way.
//
// The same step can be made with an integral instead of the abstract row: an
// integral whose reduction to the masters has a pole at f lies in the closure
// (it is regular at f) and outside M. Replacing a master by it is the known
// remedy "choose another dotted integral of the sector" (Smirnov and Smirnov
// 2002.08042, Usovitsch 2002.08173), and usually gives the smaller
// connection. Both kinds of step are tried and the smaller result is kept.
#include "diffexp/denominator_factors.hpp"
#include "diffexp/recursion_graph.hpp"
#include "diffexp/recursion_physical_basis.hpp"
#include "diffexp/artifact_store.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <unistd.h>

namespace diffexp::moving_poles {
using Matrix = ExactEpsilonMatrix;
struct Script;
struct Options {
  // Steps in all, and per factor: a relation of order n takes n steps.
  std::size_t max_steps = 64, max_factor_steps = 8;
  double seconds = 60;
  // Experimental: also treat the singular points in x alone that lie at
  // neither end of the interval. Those that are apparent are removed like
  // the moving ones; a genuine one exhausts its step bound and is kept.
  bool fixed_points = false;
  // Steps tried per factor: the closure row at up to two positions, and each
  // integral with a pole at the factor at its last max_positions positions.
  std::size_t max_positions = 2, max_integrals = 12;
  bool closure_rows = true;
  // A factor whose removal leaves the connection more than this many times
  // larger (in characters) is kept. Transport and the exact endpoint series
  // cost in proportion to the size of the connection, and a moving pole
  // costs a detour: a closure row in the residue field of a factor of high
  // degree in x can multiply the size by dozens. Zero: no bound.
  double max_growth = 4;
  // Optional restriction of the masters an integral may replace.
  std::function<bool(const ibp::Integral &, std::size_t)> admissible;
  // The outcome of an earlier search on the same input: only its steps are
  // tried and the factors it left in place are not attempted. The steps are
  // applied and checked like any others; see remove_cached.
  const Script *script = nullptr;
};
// A coordinate row in the ORIGINAL masters, e.g. the exact reduction of an
// integral (then named; empty for a row that is no single integral).
struct Row {
  std::vector<Exact> coordinates;
  std::string label;
  ibp::Integral integral;
};
struct Step {
  std::string factor, source;
  std::size_t replaced;
  unsigned order;
  std::size_t trials;
  // The integral that took the master's place; empty for a closure row.
  ibp::Integral integral;
};
struct Result {
  recursion::PhysicalBasisTransform transform;
  denominator_factors::MixedWeight before, after;
  std::size_t original_characters = 0, characters = 0;
  std::vector<Step> steps;
  // Target factors left in place: nothing tried cleared them.
  std::vector<std::string> remaining;
  // Of those, the factors whose removal was found and undone for its size
  // (Options::max_growth).
  std::size_t size_refusals = 0;
  bool budget_reached = false;
  double seconds = 0;
  bool changed() const { return !steps.empty(); }
};
// What a completed search decided: enough to repeat it without the search.
struct Script {
  std::vector<Step> steps;
  std::vector<std::string> remaining;
  std::size_t size_refusals = 0;
};
inline std::size_t characters(const Matrix &matrix) {
  std::size_t count = 0;
  for (const auto &row : matrix)
    for (const auto &value : row)
      if (!value.is_zero())
        count += value.str().size();
  return count;
}

// J' = H J with H the identity except for row pos, which is q: the rank-one
// update of physical_basis_scan, on a complete transform.
inline void replace_row(recursion::PhysicalBasisTransform &basis, std::size_t pos,
                        const std::vector<Exact> &q, std::size_t xi) {
  const auto d = basis.connection.size();
  const auto zero = q[pos].constant(0), one = q[pos].constant(1);
  if (q[pos].is_zero())
    throw std::invalid_argument("moving-pole row replacement is singular");
  const auto row_times = [&](const Matrix &m) {
    std::vector<Exact> out(d, zero);
    for (std::size_t k = 0; k < d; ++k)
      if (!q[k].is_zero())
        for (std::size_t j = 0; j < d; ++j)
          if (!m[k][j].is_zero())
            out[j] = out[j] + q[k] * m[k][j];
    return out;
  };
  std::vector<Exact> inverse_row(d, zero);
  for (std::size_t k = 0; k < d; ++k)
    inverse_row[k] = k == pos ? one / q[pos] : -q[k] / q[pos];
  const auto times_inverse = [&](Matrix &m) {
    for (std::size_t i = 0; i < d; ++i) {
      const Exact pivot = m[i][pos];
      if (pivot.is_zero())
        continue;
      for (std::size_t k = 0; k < d; ++k)
        if (k == pos)
          m[i][k] = pivot * inverse_row[k];
        else if (!inverse_row[k].is_zero())
          m[i][k] = m[i][k] + pivot * inverse_row[k];
    }
  };
  auto connection_row = row_times(basis.connection);
  for (std::size_t k = 0; k < d; ++k)
    connection_row[k] = connection_row[k] + q[k].derivative(xi);
  basis.to_basis[pos] = row_times(basis.to_basis);
  basis.connection[pos] = std::move(connection_row);
  times_inverse(basis.connection);
  times_inverse(basis.from_basis);
}

// The epsilon content of a row: a unit of the local ring at any mixed factor,
// so dividing it out changes neither the lattice nor any pole in x, and keeps
// the new master from carrying a spurious power of epsilon.
inline std::vector<Exact> primitive_in_epsilon(std::vector<Exact> row) {
  std::optional<Exact> denominator, content;
  for (const auto &value : row)
    if (!value.is_zero())
      denominator = denominator ? denominator->polynomial_lcm(value.denominator()) : value.denominator();
  if (!denominator)
    return row;
  for (auto &value : row)
    if (!value.is_zero()) {
      value = value * *denominator;
      // gcd(a,b) = a b / lcm(a,b), up to a rational factor.
      content = content ? *content * value / content->polynomial_lcm(value) : value;
    }
  for (auto &value : row)
    if (!value.is_zero())
      value = value / *content;
  return row;
}

// Arithmetic modulo an irreducible factor f, as a polynomial in the parameter
// over the field of rational functions of the other variables: the residue
// field at f. An element is a polynomial of degree below that of f, stored by
// its coefficients. For a linear factor this is the substitution x=x0.
class ResidueField {
public:
  using Polynomial = std::vector<Exact>;
  ResidueField(const Exact &f, std::size_t xi) : xi_(xi), zero_(f.constant(0)), one_(f.constant(1)) {
    modulus_ = coefficients(f);
    if (modulus_.size() < 2)
      throw std::invalid_argument("a residue field needs a factor that depends on the parameter");
    const auto lead = modulus_.back();
    for (auto &c : modulus_)
      c = c / lead;
  }
  std::size_t degree() const { return modulus_.size() - 1; }
  // An integer polynomial by powers of the parameter.
  Polynomial coefficients(const Exact &polynomial) const {
    Polynomial out;
    for (const auto &term : polynomial.numerator_terms()) {
      const auto power = term.powers.at(xi_);
      if (out.size() <= power)
        out.resize(power + 1, zero_);
      auto c = zero_.constant(term.coefficient);
      for (std::size_t v = 0; v < term.powers.size(); ++v)
        if (v != xi_ && term.powers[v])
          c = c * zero_.variable(v).pow(term.powers[v]);
      out[power] = out[power] + c;
    }
    return out;
  }
  Polynomial reduce(Polynomial p) const {
    const auto n = degree();
    while (p.size() > n) {
      const auto lead = p.back();
      p.pop_back();
      if (lead.is_zero())
        continue;
      const auto shift = p.size() - n;
      for (std::size_t i = 0; i < n; ++i)
        if (!modulus_[i].is_zero())
          p[shift + i] = p[shift + i] - lead * modulus_[i];
    }
    p.resize(n, zero_);
    return p;
  }
  Polynomial multiply(const Polynomial &a, const Polynomial &b) const { return reduce(product(a, b)); }
  // Extended Euclid against f; f is irreducible, so only zero has no inverse.
  Polynomial inverse(const Polynomial &a) const {
    auto r0 = modulus_, r1 = trimmed(a);
    Polynomial t0, t1{one_};
    if (r1.empty())
      throw std::domain_error("zero has no inverse in a residue field");
    while (!r1.empty()) {
      Polynomial quotient(r0.size() >= r1.size() ? r0.size() - r1.size() + 1 : 0, zero_);
      for (std::size_t i = quotient.size(); i-- > 0;) {
        const auto c = r0[i + r1.size() - 1] / r1.back();
        quotient[i] = c;
        if (!c.is_zero())
          for (std::size_t j = 0; j < r1.size(); ++j)
            r0[i + j] = r0[i + j] - c * r1[j];
      }
      auto t2 = product(quotient, t1);
      t2.resize(std::max(t2.size(), t0.size()), zero_);
      for (std::size_t i = 0; i < t2.size(); ++i)
        t2[i] = (i < t0.size() ? t0[i] : zero_) - t2[i];
      auto remainder = trimmed(std::move(r0));
      r0 = std::move(r1);
      r1 = std::move(remainder);
      t0 = std::move(t1);
      t1 = trimmed(std::move(t2));
    }
    if (r0.size() != 1)
      throw std::domain_error("a residue shares a factor with the modulus, which is then not irreducible");
    for (auto &c : t0)
      c = c / r0[0];
    return reduce(std::move(t0));
  }
  // The residue of a rational function whose denominator is prime to f.
  Polynomial residue(const Exact &value) const {
    if (value.is_zero())
      return Polynomial(degree(), zero_);
    return multiply(reduce(coefficients(value.numerator())), inverse(reduce(coefficients(value.denominator()))));
  }
  bool is_zero(const Polynomial &p) const {
    return std::all_of(p.begin(), p.end(), [](const Exact &c) { return c.is_zero(); });
  }
  Exact value(const Polynomial &p) const {
    auto out = zero_;
    for (std::size_t i = p.size(); i-- > 0;)
      out = out * zero_.variable(xi_) + p[i];
    return out;
  }

private:
  std::size_t xi_;
  Exact zero_, one_;
  Polynomial modulus_;
  Polynomial product(const Polynomial &a, const Polynomial &b) const {
    if (a.empty() || b.empty())
      return {};
    Polynomial out(a.size() + b.size() - 1, zero_);
    for (std::size_t i = 0; i < a.size(); ++i)
      if (!a[i].is_zero())
        for (std::size_t j = 0; j < b.size(); ++j)
          if (!b[j].is_zero())
            out[i + j] = out[i + j] + a[i] * b[j];
    return out;
  }
  static Polynomial trimmed(Polynomial p) {
    while (!p.empty() && p.back().is_zero())
      p.pop_back();
    return p;
  }
};

// The integrals a level can name without new IBP work: every cached identity
// that expresses one integral through the selected masters is its reduction.
inline std::vector<Row> cached_integrals(const level::Result &level) {
  std::vector<Row> rows;
  const auto d = level.ordered_basis.size();
  if (!d)
    return rows;
  const auto zero = level.matrix[0][0].constant(0);
  std::map<ibp::Integral, std::size_t> position;
  for (std::size_t i = 0; i < d; ++i)
    position.emplace(level.ordered_basis[i], i);
  std::set<ibp::Integral> seen;
  const auto label = [](const ibp::Integral &a) {
    std::string out = "integral [";
    for (std::size_t i = 0; i < a.size(); ++i)
      out += (i ? "," : "") + std::to_string(a[i]);
    return out + "]";
  };
  for (std::size_t i = 0; i < level.target_rows.size(); ++i)
    if (!position.count(level.requested_integrals[i]) && seen.insert(level.requested_integrals[i]).second)
      rows.push_back({level.target_rows[i], label(level.requested_integrals[i]), level.requested_integrals[i]});
  for (const auto &relation : level.source_identities) {
    const ibp::Integral *reduced = nullptr;
    bool usable = true;
    for (const auto &[integral, coefficient] : relation)
      if (!position.count(integral)) {
        usable = usable && !reduced && coefficient.is_one();
        reduced = &integral;
      }
    if (!usable || !reduced || !seen.insert(*reduced).second)
      continue;
    std::vector<Exact> row(d, zero);
    for (const auto &[integral, coefficient] : relation)
      if (&integral != reduced)
        row[position.at(integral)] = -coefficient;
    rows.push_back({std::move(row), label(*reduced), *reduced});
  }
  return rows;
}

inline Result remove(const Matrix &connection, std::size_t xi, std::size_t ei,
                     const std::vector<Row> &integrals = {}, const Options &options = {}) {
  namespace df = denominator_factors;
  const auto start = std::chrono::steady_clock::now();
  const auto elapsed = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(); };
  const auto d = connection.size();
  if (!d || connection[0].size() != d)
    throw std::invalid_argument("moving-pole removal needs a square connection");
  const auto &z = connection[0][0];
  Matrix unit(d, std::vector<Exact>(d, z.constant(0)));
  for (std::size_t i = 0; i < d; ++i)
    unit[i][i] = z.constant(1);
  Result result{{connection, unit, unit, connection}, {}, {}, 0, 0, {}, {}, false, 0};
  df::Catalogue catalogue(xi, ei);
  // The factors to remove: the mixed ones and, on request, those of x alone
  // that vanish at neither end of the integration interval. Whether one of
  // the latter is apparent is not known beforehand; the closure decides.
  std::vector<Exact> at_zero, at_one;
  for (std::size_t v = 0; v < z.variable_count(); ++v) {
    at_zero.push_back(v == xi ? z.constant(0) : z.variable(v));
    at_one.push_back(v == xi ? z.constant(1) : z.variable(v));
  }
  std::map<std::size_t, bool> targets;
  const auto target = [&](std::size_t k) {
    const auto &factor = catalogue.factors()[k];
    if (factor.kind == df::Kind::mixed)
      return true;
    if (!options.fixed_points || factor.kind != df::Kind::parameter)
      return false;
    auto [found, inserted] = targets.try_emplace(k, false);
    if (inserted)
      found->second = !factor.polynomial.substitute(at_zero).is_zero() && !factor.polynomial.substitute(at_one).is_zero();
    return found->second;
  };
  const auto weigh = [&](const Matrix &matrix) {
    std::map<std::size_t, unsigned> highest_power;
    df::MixedWeight weight;
    for (const auto &matrix_row : matrix)
      for (const auto &value : matrix_row) {
        if (value.is_zero())
          continue;
        bool open = false;
        for (const auto &[k, power] : catalogue.split(value.denominator()))
          if (target(k)) {
            open = true;
            auto &seen = highest_power[k];
            seen = std::max(seen, power);
          }
        weight.entries += open;
      }
    weight.factors = highest_power.size();
    for (const auto &[k, power] : highest_power)
      weight.total_order += power;
    return weight;
  };
  result.before = result.after = weigh(connection);
  result.original_characters = result.characters = characters(connection);
  // Only an integral whose reduction has such a pole can change a lattice.
  std::vector<const Row *> useful;
  if (result.before.factors)
    for (const auto &row : integrals) {
      if (row.coordinates.size() != d)
        throw std::invalid_argument("moving-pole integral row width");
      bool pole = false;
      for (const auto &value : row.coordinates)
        if (!value.is_zero())
          for (const auto &[k, power] : catalogue.split(value.denominator()))
            pole = pole || target(k);
      if (pole)
        useful.push_back(&row);
    }
  // Singular points of the original connection: a step may not add one.
  std::set<std::size_t> singular;
  for (const auto &matrix_row : connection)
    for (const auto &value : matrix_row)
      if (!value.is_zero())
        for (const auto &[k, power] : catalogue.split(value.denominator()))
          singular.insert(k);
  const auto new_singular_point = [&](const Matrix &matrix) {
    for (const auto &matrix_row : matrix)
      for (const auto &value : matrix_row)
        if (!value.is_zero())
          for (const auto &[k, power] : catalogue.split(value.denominator()))
            if (catalogue.factors()[k].kind != df::Kind::epsilon && !singular.count(k))
              return true;
    return false;
  };
  std::set<std::size_t> abandoned;
  auto &basis = result.transform;
  // The highest pole of one factor in the current connection, and a row
  // that carries it.
  const auto highest = [&](std::size_t factor) {
    std::pair<unsigned, std::size_t> found{0, 0};
    for (std::size_t i = 0; i < d; ++i)
      for (std::size_t j = 0; j < d; ++j)
        if (!basis.connection[i][j].is_zero())
          for (const auto &[k, power] : catalogue.split(basis.connection[i][j].denominator()))
            if (k == factor && power > found.first)
              found = {power, i};
    return found;
  };
  for (;;) {
    // The first target factor still open.
    std::optional<std::size_t> factor;
    for (const auto &matrix_row : basis.connection)
      for (const auto &value : matrix_row)
        if (!value.is_zero())
          for (const auto &[k, power] : catalogue.split(value.denominator()))
            if (target(k) && !abandoned.count(k) && (!factor || k < *factor))
              factor = k;
    if (!factor)
      break;
    if (result.steps.size() >= options.max_steps || elapsed() >= options.seconds) {
      result.budget_reached = true;
      break;
    }
    const auto f = catalogue.factors()[*factor].polynomial;
    const auto factor_name = f.str();
    if (options.script && std::find(options.script->remaining.begin(), options.script->remaining.end(), factor_name) !=
                              options.script->remaining.end()) {
      abandoned.insert(*factor);
      continue;
    }
    // A factor is cleared completely or not touched: steps that do not end
    // (a singular point that is not apparent, or a relation longer than the
    // step bound) are undone.
    const auto snapshot = basis;
    const auto steps_before = result.steps.size();
    const auto weight_before = result.after;
    const auto characters_before = result.characters;
    const auto too_large = [&] {
      return options.max_growth > 0 &&
             static_cast<double>(result.characters) > options.max_growth * static_cast<double>(characters_before);
    };
    for (std::size_t local = 0; local < options.max_factor_steps; ++local) {
      const auto [order, row] = highest(*factor);
      if (!order)
        break;
      struct Trial {
        recursion::PhysicalBasisTransform basis;
        df::MixedWeight weight;
        std::size_t characters, replaced;
        std::string source;
        ibp::Integral integral;
      };
      std::optional<Trial> best;
      std::size_t trials = 0;
      // A closure step always enlarges the lattice, so it is admitted even
      // when the pole order does not drop at once; an integral only if it
      // does, and if the other denominators of its reduction add no singular
      // point in x.
      const auto attempt = [&](const std::vector<Exact> &q, std::size_t replaced, const std::string &source, bool closure,
                               const ibp::Integral &named = {}) {
        ++trials;
        if (options.script) {
          // Only the step the earlier search took at this point.
          const auto index = result.steps.size();
          if (index >= options.script->steps.size())
            return;
          const auto &expected = options.script->steps[index];
          if (expected.factor != factor_name || expected.replaced != replaced || expected.source != source)
            return;
        }
        auto next = basis;
        replace_row(next, replaced, q, xi);
        const auto weight = weigh(next.connection);
        if (!closure && (weight.total_order >= result.after.total_order || new_singular_point(next.connection)))
          return;
        const auto size = characters(next.connection);
        if (!best || std::tie(weight.total_order, size) < std::tie(best->weight.total_order, best->characters))
          best = Trial{std::move(next), weight, size, replaced, source, named};
      };
      if (options.closure_rows) {
        const ResidueField field(f, xi);
        const auto scale = f.pow(order);
        std::vector<ResidueField::Polynomial> residues(d);
        std::vector<std::size_t> support;
        for (std::size_t j = 0; j < d; ++j)
          if (!basis.connection[row][j].is_zero()) {
            residues[j] = field.residue(basis.connection[row][j] * scale);
            if (!field.is_zero(residues[j]))
              support.push_back(j);
          }
        if (support.empty())
          throw std::logic_error("moving-pole residue row vanished");
        std::vector<Exact> r(d, z.constant(0));
        if (field.degree() == 1) {
          // The relation depends on epsilon alone. 1/r_k enters the old
          // masters and the new connection. The last master of the relation
          // keeps the new coordinate in the highest sector that carries it;
          // the last one whose coefficient does not vanish at epsilon = 0
          // (one exists, r being primitive) keeps the new basis from
          // degenerating there. They can differ, and both are tried.
          for (auto j : support)
            r[j] = residues[j][0];
          r = primitive_in_epsilon(std::move(r));
          std::optional<std::size_t> regular;
          for (std::size_t j = d; j-- > 0 && !regular;)
            if (!r[j].is_zero() && *exact_epsilon_valuation(r[j], ei) == 0)
              regular = j;
          for (auto j : support)
            r[j] = r[j] / f;
          attempt(r, support.back(), "closure row", true);
          if (regular && *regular != support.back())
            attempt(r, *regular, "closure row, coefficient regular at epsilon 0", true);
        } else {
          // The relation has coefficients in the residue field. Normalised to
          // r_k = 1 the replaced master returns as f times the new one minus
          // polynomials, so that no root of r_k becomes a singular point.
          for (std::size_t tried = 0; tried < options.max_positions && tried < support.size(); ++tried) {
            const auto k = support[support.size() - 1 - tried];
            const auto inverse = field.inverse(residues[k]);
            for (auto j : support)
              r[j] = (j == k ? z.constant(1) : field.value(field.multiply(residues[j], inverse))) / f;
            attempt(r, k, "closure row", true);
          }
        }
      }
      std::size_t used = 0;
      for (const auto *integral : useful) {
        if (used >= options.max_integrals || elapsed() >= options.seconds)
          break;
        // The integral in the current coordinates.
        std::vector<Exact> q(d, z.constant(0));
        for (std::size_t k = 0; k < d; ++k)
          if (!integral->coordinates[k].is_zero())
            for (std::size_t j = 0; j < d; ++j)
              if (!basis.from_basis[k][j].is_zero())
                q[j] = q[j] + integral->coordinates[k] * basis.from_basis[k][j];
        // A master can be exchanged for it only where its coefficient has
        // the pole: elsewhere the old master would leave the lattice.
        std::size_t positions = 0;
        for (std::size_t j = d; j-- > 0 && positions < options.max_positions;)
          if (!q[j].is_zero() && catalogue.order(q[j], *factor) < 0 &&
              (!options.admissible || options.admissible(integral->integral, j))) {
            attempt(q, j, integral->label, false, integral->integral);
            ++positions;
          }
        used += positions > 0;
      }
      if (!best)
        break;
      result.steps.push_back({factor_name, best->source, best->replaced, order, trials, best->integral});
      basis = std::move(best->basis);
      result.after = best->weight;
      result.characters = best->characters;
      if (too_large())
        break;
    }
    if (highest(*factor).first || too_large()) {
      if (!highest(*factor).first)
        ++result.size_refusals;
      basis = snapshot;
      result.steps.resize(steps_before);
      result.after = weight_before;
      result.characters = characters_before;
      abandoned.insert(*factor);
    }
  }
  std::map<std::size_t, bool> open;
  for (const auto &matrix_row : basis.connection)
    for (const auto &value : matrix_row)
      if (!value.is_zero())
        for (const auto &[k, power] : catalogue.split(value.denominator()))
          if (target(k) && open.emplace(k, true).second)
            result.remaining.push_back(catalogue.factors()[k].polynomial.str());
  result.seconds = elapsed();
  return result;
}

// remove(), with the outcome of the search kept in `directory` (nothing is
// kept if it is empty). The search tries a dozen exchanges per factor on the
// whole connection, and every process that prepares the same level (each
// epsilon sample of a family) would repeat it. A stored outcome names the
// steps that were taken and the factors that were left in place; they are
// applied again by the same code, so the result is the one of the search or,
// if anything about it differs, the search is run. `tag` names what the
// options cannot be hashed for (Options::admissible).
inline Result remove_cached(const std::filesystem::path &directory, const std::string &tag, const Matrix &connection,
                            std::size_t xi, std::size_t ei, const std::vector<Row> &integrals, Options options) {
  namespace json = boost::json;
  if (directory.empty() || options.script)
    return remove(connection, xi, ei, integrals, options);
  std::string text = "moving-poles-script-v1|" + tag + "|" + std::to_string(xi) + "," + std::to_string(ei) + "," +
                     std::to_string(options.max_steps) + "," + std::to_string(options.max_factor_steps) + "," +
                     std::to_string(options.fixed_points) + "," + std::to_string(options.max_positions) + "," +
                     std::to_string(options.max_integrals) + "," + std::to_string(options.closure_rows) + "," +
                     std::to_string(options.max_growth) + "|";
  for (const auto &name : connection.at(0).at(0).variables())
    text += name + ",";
  for (const auto &row : connection)
    for (const auto &value : row)
      text += value.is_zero() ? std::string(";") : value.str() + ";";
  for (const auto &row : integrals) {
    text += "|" + row.label + ":";
    for (const auto &value : row.coordinates)
      text += value.is_zero() ? std::string(";") : value.str() + ";";
  }
  const auto file = directory / (artifacts::detail::sha256(text) + ".json");
  const auto same = [](const Result &result, const Script &script) {
    if (result.steps.size() != script.steps.size() || result.remaining != script.remaining || result.budget_reached)
      return false;
    for (std::size_t i = 0; i < result.steps.size(); ++i)
      if (result.steps[i].factor != script.steps[i].factor || result.steps[i].replaced != script.steps[i].replaced ||
          result.steps[i].source != script.steps[i].source)
        return false;
    return true;
  };
  try {
    std::ifstream in(file);
    if (in) {
      const std::string bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      const auto record = json::parse(bytes);
      const auto &object = record.as_object();
      Script script;
      for (const auto &item : object.at("steps").as_array()) {
        const auto &step = item.as_object();
        Step value{std::string(step.at("factor").as_string()), std::string(step.at("source").as_string()),
                   static_cast<std::size_t>(step.at("replaced").to_number<std::uint64_t>()),
                   static_cast<unsigned>(step.at("order").to_number<std::uint64_t>()),
                   static_cast<std::size_t>(step.at("trials").to_number<std::uint64_t>()), {}};
        for (const auto &power : step.at("integral").as_array())
          value.integral.push_back(static_cast<int>(power.to_number<std::int64_t>()));
        script.steps.push_back(std::move(value));
      }
      for (const auto &factor : object.at("remaining").as_array())
        script.remaining.emplace_back(factor.as_string());
      script.size_refusals = static_cast<std::size_t>(object.at("size_refusals").to_number<std::uint64_t>());
      auto scripted = options;
      scripted.script = &script;
      auto result = remove(connection, xi, ei, integrals, scripted);
      if (same(result, script)) {
        for (std::size_t i = 0; i < result.steps.size(); ++i)
          result.steps[i].trials = script.steps[i].trials;
        result.size_refusals = script.size_refusals;
        return result;
      }
    }
  } catch (const std::exception &) {
    // An unreadable or stale record: the search decides.
  }
  auto result = remove(connection, xi, ei, integrals, options);
  if (result.budget_reached)
    return result;
  try {
    json::array steps, remaining;
    for (const auto &step : result.steps) {
      json::array integral;
      for (int power : step.integral)
        integral.push_back(power);
      steps.push_back(json::object{{"factor", step.factor}, {"source", step.source}, {"replaced", step.replaced},
                                   {"order", step.order}, {"trials", step.trials}, {"integral", std::move(integral)}});
    }
    for (const auto &factor : result.remaining)
      remaining.emplace_back(factor);
    const auto bytes = json::serialize(json::object{{"schema", "DiffExp.MovingPoleScript/v1"}, {"steps", std::move(steps)},
                                                    {"remaining", std::move(remaining)}, {"size_refusals", result.size_refusals}});
    std::filesystem::create_directories(directory);
    const auto temporary = directory / (file.filename().string() + "." + std::to_string(::getpid()) + ".tmp");
    {
      std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
      out << bytes;
      if (!out)
        throw std::runtime_error("moving-pole record not written");
    }
    std::filesystem::rename(temporary, file);
  } catch (const std::exception &) {
    // The record is a convenience; the result does not depend on it.
  }
  return result;
}

// The internal bases of a prepared recursion in which the moving poles of
// every level are removed, keyed by zero-based depth as
// recursion::NumericalOptions::physical_bases expects them; levels without
// moving poles, or where nothing lowered them, are absent. The graph must
// still hold epsilon: a recursion specialised to a numeric dimension shows
// its moving poles as ordinary ones. Shift the result with shift_epsilon.
using Report = std::function<void(std::size_t, const Result &)>;
inline std::map<std::size_t, recursion::PhysicalBasisTransform> bases(const recursion::Graph &graph,
    const Options &options = {}, const Report &report = {}, const std::filesystem::path &directory = {}) {
  if (graph.specialised)
    throw std::invalid_argument("moving poles are found before a recursion is specialised to a numeric dimension");
  const auto &names = graph.dimension.variables();
  const auto x = std::find(names.begin(), names.end(), "x"), eps = std::find(names.begin(), names.end(), "eps");
  if (x == names.end() || eps == names.end())
    throw std::invalid_argument("moving-pole removal needs the x and eps variables");
  std::map<std::size_t, recursion::PhysicalBasisTransform> found;
  for (std::size_t depth = 0; depth < graph.nodes.size(); ++depth) {
    const auto &node = graph.nodes[depth];
    if (node.scalar_leaf || node.closure.matrix.empty())
      continue;
    auto result = remove_cached(directory, "bases", node.closure.matrix, x - names.begin(), eps - names.begin(),
                                cached_integrals(node.closure), options);
    if (result.changed())
      recursion::verify_physical_basis(result.transform, node.closure.matrix, x - names.begin());
    if (report)
      report(depth, result);
    if (result.changed())
      found.emplace(depth, std::move(result.transform));
  }
  return found;
}
// The level in masters without moving poles, as far as exchanging a master
// for a cached integral of the same sector achieves it. Unlike an internal
// basis this changes what the level asks of its child, so the recursion below
// it is prepared for the new masters; every later stage (basis scan, epsilon
// shift, caches) then sees an ordinary level. The new closure is the certified
// one transformed by an exact, verified change of basis; the IBP witnesses of
// the old masters do not apply to it and are dropped.
inline level::Result exchange(level::Result level, std::size_t xi, std::size_t ei, Options options = {},
                              Result *report = nullptr, const std::filesystem::path &directory = {}) {
  if (!level.success || level.matrix.empty())
    return level;
  options.closure_rows = false;
  const auto physical = [](const ibp::Integral &a) {
    std::vector<bool> sector;
    for (int power : a)
      sector.push_back(power > 0);
    return sector;
  };
  const auto masters = level.ordered_basis;
  options.admissible = [&](const ibp::Integral &integral, std::size_t position) {
    return !integral.empty() && physical(integral) == physical(masters[position]);
  };
  // The admissible positions follow from the masters, which the record names.
  std::string tag = "exchange";
  for (const auto &master : masters) {
    tag += "[";
    for (int power : master)
      tag += std::to_string(power) + ",";
    tag += "]";
  }
  auto result = remove_cached(directory, tag, level.matrix, xi, ei, cached_integrals(level), options);
  if (result.changed()) {
    recursion::verify_physical_basis(result.transform, level.matrix, xi);
    if (!level.target_rows.empty())
      level.target_rows = fuchsify::detail::multiply(level.target_rows, result.transform.from_basis);
    level.matrix = result.transform.connection;
    for (const auto &step : result.steps)
      level.ordered_basis[step.replaced] = step.integral;
    level.differential_witnesses.clear();
    level.target_witnesses.clear();
  }
  if (report)
    *report = std::move(result);
  return level;
}
// The companion of recursion::shift_epsilon for an internal basis of a
// non-leaf level: eps -> eps+shift, or eps -> shift when the levels are
// specialised to the number.
inline recursion::PhysicalBasisTransform shift_epsilon(recursion::PhysicalBasisTransform basis,
    const Rational &shift, bool specialise = true) {
  const auto &z = basis.connection.at(0).at(0);
  const auto &names = z.variables();
  const auto eps = std::find(names.begin(), names.end(), "eps");
  if (eps == names.end())
    throw std::invalid_argument("epsilon shift needs the eps variable");
  std::vector<Exact> point;
  for (std::size_t i = 0; i < z.variable_count(); ++i)
    point.push_back(z.variable(i));
  point[eps - names.begin()] = specialise ? z.constant(shift) : point[eps - names.begin()] + z.constant(shift);
  for (auto *matrix : {&basis.original_connection, &basis.to_basis, &basis.from_basis, &basis.connection})
    for (auto &row : *matrix)
      for (auto &value : row)
        if (!value.is_zero())
          value = value.substitute(point);
  return basis;
}
} // namespace diffexp::moving_poles
