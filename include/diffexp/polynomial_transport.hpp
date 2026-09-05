#pragma once
#include "diffexp/laurent_transport.hpp"
#include <flint/acb_poly.h>

namespace diffexp::polynomial_transport {
using B = Jet::Ball;
struct Options {
  unsigned max_dimension = 5000, max_degree = 256, max_epsilon_degree = 256;
  std::size_t max_terms = 200000, max_cells = 20000000,
              max_operations = 200000000;
  bool allow_fallback = true;
};
struct Scalar {
  Rational real{0}, imaginary{0};
};
using Polynomial = std::vector<Scalar>;
struct Block {
  unsigned epsilon, polynomial;
};
struct Entry {
  unsigned column;
  std::vector<Block> coefficients;
};
struct Row {
  bool polynomial = true;
  std::vector<Block> denominator;
  std::vector<Entry> entries;
};
// Exact coefficients are interned once and reused for every ordinary chart,
// including equal rows in the block-diagonal adjoint system.
struct SharedFallback {
  unsigned epsilon;
  data::Expr coefficient;
  std::vector<std::pair<unsigned, unsigned>> positions;
};
struct Compiled {
  unsigned dimension = 0, epsilon_high = 0, expected_order = 0,
           polynomial_rows = 0, fallback_rows = 0;
  std::vector<Polynomial> polynomials;
  std::vector<Row> rows;
  std::vector<SharedFallback> fallback;
  Options options;
};
namespace detail {
inline B ball(const Rational &q) {
  auto s = q.str();
  auto slash = s.find('/');
  return slash == std::string::npos ? B::from_strings(s)
                                    : B::from_strings(s.substr(0, slash)) /
                                          B::from_strings(s.substr(slash + 1));
}
inline void validate_parameters(const Exact &e, std::size_t xi,
                                std::optional<std::size_t> ei) {
  for (const auto &list : {e.numerator_terms(), e.denominator_terms()})
    for (const auto &term : list)
      for (unsigned i = 0; i < term.powers.size(); ++i)
        if (term.powers[i] && i != xi && (!ei || i != *ei) &&
            e.variables()[i] != "I")
          throw std::invalid_argument(
              "polynomial transport has an unsubstituted parameter");
}
struct Compiler {
  Compiled output;
  std::map<std::string, unsigned> interned;
  std::size_t terms = 0;
  std::vector<Block> blocks(const Exact &e, std::size_t xi,
                            std::optional<std::size_t> ei,
                            unsigned offset = 0) {
    if (!e.denominator().is_rational())
      throw std::logic_error(
          "denominator clearing did not produce a polynomial");
    const auto divisor = e.denominator().rational();
    std::map<unsigned, std::map<unsigned, Scalar>> groups;
    for (const auto &t : e.numerator_terms()) {
      if (++terms > output.options.max_terms)
        throw std::length_error("polynomial clearing term budget exhausted");
      auto xpower = t.powers[xi];
      unsigned epower = offset + (ei ? t.powers[*ei] : 0);
      if (xpower > output.options.max_degree ||
          epower > output.options.max_epsilon_degree)
        throw std::length_error("polynomial clearing degree budget exhausted");
      auto q = t.coefficient / divisor;
      unsigned imaginary = 0;
      for (unsigned i = 0; i < t.powers.size(); ++i)
        if (e.variables()[i] == "I")
          imaginary = t.powers[i] % 4;
      if (imaginary >= 2)
        q = -q;
      auto &target = groups[epower][xpower];
      if (imaginary % 2)
        target.imaginary += q;
      else
        target.real += q;
    }
    std::vector<Block> result;
    for (const auto &[epsilon, coefficients] : groups) {
      Polynomial p(coefficients.rbegin()->first + 1);
      for (const auto &[n, c] : coefficients)
        p[n] = c;
      while (!p.empty() && p.back().real.is_zero() &&
             p.back().imaginary.is_zero())
        p.pop_back();
      if (p.empty())
        continue;
      std::string key;
      for (const auto &c : p)
        key += c.real.str() + "," + c.imaginary.str() + ";";
      auto [it, inserted] = interned.emplace(key, output.polynomials.size());
      if (inserted)
        output.polynomials.push_back(std::move(p));
      result.push_back({epsilon, it->second});
    }
    return result;
  }
};
inline std::size_t cells(std::initializer_list<std::size_t> factors,
                         std::size_t maximum) {
  std::size_t n = 1;
  for (auto k : factors) {
    if (k && n > maximum / k)
      throw std::length_error("polynomial transport storage budget exhausted");
    n *= k;
  }
  return n;
}
} // namespace detail
// General sparse input can retain rational epsilon dependence. For expanded
// entries, epsilon_offset carries their preexisting epsilon monomial.
inline Compiled compile(const std::vector<RationalLineEntry> &entries,
                        unsigned dimension, std::size_t xi,
                        std::optional<std::size_t> ei, unsigned epsilon_high,
                        unsigned expected_order, const Options &options = {}) {
  if (entries.size() > options.max_terms || !dimension ||
      dimension > options.max_dimension || !expected_order ||
      expected_order > 1000 || epsilon_high > options.max_epsilon_degree ||
      !options.max_terms || !options.max_cells || !options.max_operations)
    throw std::invalid_argument(
        "polynomial transport compilation dimensions or budgets");
  detail::Compiler compiler;
  auto &out = compiler.output;
  out.dimension = dimension;
  out.epsilon_high = epsilon_high;
  out.expected_order = expected_order;
  out.options = options;
  out.rows.resize(dimension);
  if (entries.empty()) {
    for (auto &row : out.rows) {
      row.denominator = {{0, 0}};
    }
    out.polynomials.push_back(Polynomial{{Rational(1), Rational(0)}});
    out.polynomial_rows = dimension;
    return out;
  }
  const auto zero = entries[0].coefficient.constant(0);
  std::vector<std::vector<const RationalLineEntry *>> byrow(dimension);
  for (const auto &entry : entries) {
    if (entry.row >= dimension || entry.column >= dimension)
      throw std::invalid_argument("polynomial transport sparse index");
    detail::validate_parameters(entry.coefficient, xi, ei);
    if (entry.epsilon <= epsilon_high && !entry.coefficient.is_zero())
      byrow[entry.row].push_back(&entry);
  }
  std::vector<RationalLineEntry> fallback;
  for (unsigned i = 0; i < dimension; ++i) {
    auto &row = out.rows[i];
    bool use = true;
    try {
      auto denominator = zero.constant(1);
      for (const auto *e : byrow[i]) {
        denominator = denominator.polynomial_lcm(e->coefficient.denominator());
        if (denominator.numerator_terms().size() > options.max_terms)
          throw std::length_error("common denominator term budget exhausted");
        for (const auto &t : denominator.numerator_terms())
          if (t.powers[xi] > options.max_degree ||
              (ei && t.powers[*ei] > options.max_epsilon_degree))
            throw std::length_error(
                "common denominator degree budget exhausted");
      }
      row.denominator = compiler.blocks(denominator, xi, ei);
      std::size_t finite_cost = 0;
      for (auto b : row.denominator)
        finite_cost += out.polynomials[b.polynomial].size();
      for (const auto *e : byrow[i]) {
        auto p = denominator * e->coefficient;
        auto blocks = compiler.blocks(p, xi, ei, e->epsilon);
        for (auto b : blocks)
          finite_cost += out.polynomials[b.polynomial].size();
        row.entries.push_back({e->column, std::move(blocks)});
      }
      // This deliberately conservative estimate does not assert a speedup.
      // Large clearing degrees retain the established rational convolution.
      std::size_t rational_cost = byrow[i].size() * ((expected_order + 1) / 2);
      if (!byrow[i].empty() && finite_cost > rational_cost)
        use = false;
    } catch (const std::length_error &) {
      use = false;
    }
    if (!use) {
      if (!options.allow_fallback)
        throw std::length_error("denominator clearing exceeds the finite-lag "
                                "cost or degree budget");
      row = {false, {}, {}};
      ++out.fallback_rows;
      for (const auto *e : byrow[i]) {
        auto expanded = feynman::scalar_functional_detail::epsilon_series(
            e->coefficient, ei, epsilon_high - e->epsilon);
        for (unsigned k = 0; k < expanded.size(); ++k)
          if (!expanded[k].is_zero()) {
            if (fallback.size() >= options.max_terms)
              throw std::length_error(
                  "polynomial fallback term budget exhausted");
            fallback.push_back({i, e->column, k + e->epsilon, expanded[k]});
          }
      }
    } else
      ++out.polynomial_rows;
  }
  std::map<std::string, unsigned> shared_fallback;
  for (const auto &entry : fallback) {
    auto key = std::to_string(entry.epsilon) + "|" + entry.coefficient.str();
    auto [it, inserted] = shared_fallback.emplace(key, out.fallback.size());
    if (inserted) {
      auto compiled = compile_rational_entries({entry});
      out.fallback.push_back(
          {entry.epsilon, std::move(compiled[0].coefficient), {}});
    }
    out.fallback[it->second].positions.emplace_back(entry.row, entry.column);
  }
  return std::move(out);
}
// Integration adapter for the existing epsilon-expanded sparse line entries.
inline Compiled compile(const std::vector<RationalLineEntry> &entries,
                        unsigned dimension, unsigned expected_order,
                        unsigned epsilon_high, const Options &options = {}) {
  if (entries.empty())
    return compile(entries, dimension, 0, std::nullopt, epsilon_high,
                   expected_order, options);
  const auto &names = entries[0].coefficient.variables();
  auto x = std::find(names.begin(), names.end(), "x"),
       eps = std::find(names.begin(), names.end(), "eps");
  if (x == names.end())
    throw std::invalid_argument("polynomial transport requires x variable");
  return compile(entries, dimension, x - names.begin(),
                 eps == names.end()
                     ? std::nullopt
                     : std::optional<std::size_t>(eps - names.begin()),
                 epsilon_high, expected_order, options);
}
// Preferred compact input: denominators remain polynomial in (x,epsilon),
// avoiding powers created by first expanding moving poles in epsilon.
inline Compiled compile(const ExactEpsilonMatrix &matrix, std::size_t xi,
                        std::size_t ei, unsigned expected_order,
                        unsigned epsilon_high, const Options &options = {}) {
  std::vector<RationalLineEntry> entries;
  for (unsigned i = 0; i < matrix.size(); ++i) {
    if (matrix[i].size() != matrix.size())
      throw std::invalid_argument("polynomial transport square matrix");
    for (unsigned j = 0; j < matrix.size(); ++j)
      if (!matrix[i][j].is_zero())
        entries.push_back({i, j, 0, matrix[i][j]});
  }
  return compile(entries, matrix.size(), xi, ei, epsilon_high, expected_order,
                 options);
}
// Computes the same retained Taylor polynomial as ordinary rational transport.
// Ball rounding is enclosed; no omitted-Taylor-tail certificate is inferred.
inline Boundary chart(const Compiled &compiled, const Boundary &boundary,
                      const B &center, const B &step, unsigned order) {
  const unsigned d = compiled.dimension;
  if (boundary.size() != d || boundary.empty() || boundary[0].empty() ||
      !order || order > 1000 || !center.is_finite() || !step.is_finite())
    throw std::invalid_argument("polynomial chart dimensions, order or point");
  unsigned high = boundary[0].size() - 1;
  if (high > compiled.epsilon_high)
    throw std::invalid_argument(
        "polynomial chart exceeds compiled epsilon window");
  for (const auto &row : boundary) {
    if (row.size() != high + 1)
      throw std::invalid_argument("polynomial chart boundary shape");
    for (const auto &b : row)
      if (!b.is_finite())
        throw std::invalid_argument("nonfinite polynomial chart boundary");
  }
  auto bits = B::precision();
  detail::cells({order + 1, d, high + 1}, compiled.options.max_cells);
  std::vector<std::vector<B>> polynomials;
  polynomials.reserve(compiled.polynomials.size());
  for (const auto &p : compiled.polynomials) {
    acb_poly_t a, b;
    acb_poly_init(a);
    acb_poly_init(b);
    for (unsigned j = 0; j < p.size(); ++j) {
      auto value = detail::ball(p[j].real);
      auto imaginary = detail::ball(p[j].imaginary);
      arb_set(acb_imagref(value.raw()), acb_realref(imaginary.raw()));
      acb_poly_set_coeff_acb(a, j, value.raw());
    }
    acb_poly_taylor_shift(b, a, center.raw(), bits);
    std::vector<B> coefficients(p.size());
    for (unsigned j = 0; j < p.size(); ++j)
      acb_poly_get_coeff_acb(coefficients[j].raw(), b, j);
    polynomials.push_back(std::move(coefficients));
    acb_poly_clear(a);
    acb_poly_clear(b);
  }
  std::vector<B> leading(d, B(1));
  for (unsigned i = 0; i < d; ++i)
    if (compiled.rows[i].polynomial) {
      leading[i] = B(0);
      for (auto block : compiled.rows[i].denominator)
        if (!block.epsilon)
          leading[i] += polynomials[block.polynomial][0];
      if (!leading[i].is_finite() || leading[i].contains_zero())
        throw std::domain_error(
            "polynomial chart denominator at (center,epsilon=0) is not "
            "separated from zero");
    }
  struct Fallback {
    unsigned epsilon;
    const std::vector<std::pair<unsigned, unsigned>> *positions;
    std::vector<B> coefficients;
  };
  std::vector<Fallback> fallback;
  Jet x(0, order + 1, bits);
  x.set(0, center);
  x.set(1, B(1));
  auto imaginary = x.constant(0);
  imaginary.set(0, B::from_strings("0", "1"));
  for (const auto &e : compiled.fallback) {
    auto jet = evaluate(e.coefficient, x, {{"x", x}, {"I", imaginary}});
    Fallback item{e.epsilon, &e.positions, {}};
    for (unsigned n = 0; n < order; ++n)
      item.coefficients.push_back(jet.at(n));
    fallback.push_back(std::move(item));
  }
  std::vector<B> values(static_cast<std::size_t>(order + 1) * d * (high + 1),
                        B(0));
  auto at = [&](unsigned n, unsigned i, unsigned k) -> B & {
    return values[(static_cast<std::size_t>(n) * d + i) * (high + 1) + k];
  };
  for (unsigned i = 0; i < d; ++i)
    for (unsigned k = 0; k <= high; ++k)
      at(0, i, k) = boundary[i][k];
  std::size_t work = 0;
  auto spend = [&] {
    if (++work > compiled.options.max_operations)
      throw std::length_error("polynomial chart operation budget exhausted");
  };
  for (unsigned n = 0; n < order; ++n) {
    for (const auto &e : fallback)
      for (const auto &[row, column] : *e.positions)
        for (unsigned k = e.epsilon; k <= high; ++k)
          for (unsigned m = 0; m <= n; ++m)
            if (!e.coefficients[m].is_zero() &&
                !at(n - m, column, k - e.epsilon).is_zero()) {
              spend();
              acb_addmul(at(n + 1, row, k).raw(), e.coefficients[m].raw(),
                         at(n - m, column, k - e.epsilon).raw(), bits);
            }
    // q[0,0]*(n+1)*y[n+1,k] is the sole unknown pivot. q[0,s>0]
    // couples only already computed epsilon orders at the same Taylor order.
    for (unsigned k = 0; k <= high; ++k)
      for (unsigned i = 0; i < d; ++i) {
        const auto &row = compiled.rows[i];
        auto &value = at(n + 1, i, k);
        if (row.polynomial) {
          for (const auto &e : row.entries)
            for (auto block : e.coefficients)
              if (block.epsilon <= k) {
                const auto &p = polynomials[block.polynomial];
                for (unsigned m = 0; m < p.size() && m <= n; ++m)
                  if (!p[m].is_zero() &&
                      !at(n - m, e.column, k - block.epsilon).is_zero()) {
                    spend();
                    acb_addmul(value.raw(), p[m].raw(),
                               at(n - m, e.column, k - block.epsilon).raw(),
                               bits);
                  }
              }
          for (auto block : row.denominator)
            if (block.epsilon <= k) {
              const auto &q = polynomials[block.polynomial];
              for (unsigned m = 0; m < q.size() && m <= n; ++m) {
                if (!m && !block.epsilon)
                  continue;
                if (!q[m].is_zero() &&
                    !at(n - m + 1, i, k - block.epsilon).is_zero()) {
                  spend();
                  B product;
                  acb_mul_ui(product.raw(), q[m].raw(), n - m + 1, bits);
                  acb_submul(value.raw(), product.raw(),
                             at(n - m + 1, i, k - block.epsilon).raw(), bits);
                }
              }
            }
          acb_div(value.raw(), value.raw(), leading[i].raw(), bits);
        }
        acb_div_ui(value.raw(), value.raw(), n + 1, bits);
      }
  }
  Boundary out(d, std::vector<B>(high + 1, B(0)));
  for (unsigned i = 0; i < d; ++i)
    for (unsigned k = 0; k <= high; ++k)
      for (unsigned n = order + 1; n-- > 0;)
        out[i][k] = out[i][k] * step + at(n, i, k);
  return out;
}
} // namespace diffexp::polynomial_transport
