#pragma once
#include "diffexp/fuchsify.hpp"
#include "diffexp/jet.hpp"
#include "diffexp/univariate_rational.hpp"
#include <flint/fmpz_mpoly_factor.h>
#include <map>
#include <tuple>
#include <functional>

namespace diffexp {
namespace cached_affine { struct StateAccess; }
namespace affine_frobenius_detail {
using Matrix = fuchsify::Matrix;
using Poly = std::vector<Exact>;
struct Exponent {
  Rational power, slope;
};
template<class Coefficient> inline void trim(std::vector<Coefficient> &p) {
  while (!p.empty() && p.back().is_zero())
    p.pop_back();
}
template<class Coefficient> inline void add(std::vector<Coefficient> &out, const std::vector<Coefficient> &p, const Coefficient &q) {
  if (q.is_zero() || p.empty())
    return;
  if (out.size() < p.size())
    out.resize(p.size(), q.constant(0));
  for (unsigned l = 0; l < p.size(); ++l)
    out[l] = out[l] + q * p[l];
  trim(out);
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
// Factor the exact characteristic polynomial; every spectral factor must be
// linear and its root exactly affine in epsilon. No sampled eigenvalues.
inline std::vector<Exponent> spectrum(const Matrix &r, std::size_t xi,
                                      std::size_t ei) {
  const auto p = characteristic(r, xi), z = p.constant(0);
  auto names = p.variables();
  std::vector<const char *> symbols;
  for (auto &n : names)
    symbols.push_back(n.c_str());
  fmpz_mpoly_ctx_t ctx;
  fmpz_mpoly_ctx_init(ctx, names.size(), ORD_DEGLEX);
  fmpz_mpoly_t raw;
  fmpz_mpoly_init(raw, ctx);
  fmpz_mpoly_factor_t factors;
  fmpz_mpoly_factor_init(factors, ctx);
  bool ok = fmpz_mpoly_set_str_pretty(raw, p.numerator().str().c_str(),
                                      symbols.data(), ctx) == 0;
  if (ok)
    ok = fmpz_mpoly_factor(factors, raw, ctx) != 0;
  std::vector<std::pair<Exact, long>> polynomials;
  if (ok)
    for (slong f = 0; f < factors->num; ++f) {
      auto poly = z;
      for (slong t = 0; t < fmpz_mpoly_length(factors->poly + f, ctx); ++t) {
        fmpz_t coefficient;
        fmpz_init(coefficient);
        fmpz_mpoly_get_term_coeff_fmpz(coefficient, factors->poly + f, t, ctx);
        char *s = fmpz_get_str(nullptr, 10, coefficient);
        auto c = z.constant(Rational(s));
        flint_free(s);
        fmpz_clear(coefficient);
        std::vector<unsigned long> powers(names.size());
        fmpz_mpoly_get_term_exp_ui(powers.data(), factors->poly + f, t, ctx);
        for (unsigned j = 0; j < powers.size(); ++j)
          if (powers[j])
            c = c * z.variable(j).pow(powers[j]);
        poly = poly + c;
      }
      polynomials.emplace_back(poly,
                               fmpz_mpoly_factor_get_exp_si(factors, f, ctx));
    }
  fmpz_mpoly_factor_clear(factors, ctx);
  fmpz_mpoly_clear(raw, ctx);
  fmpz_mpoly_ctx_clear(ctx);
  if (!ok)
    throw std::runtime_error("exact affine residue factorization failed");
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
  fmpq_t v;
  fmpq_init(v);
  fmpq_set_str(v, q.str().c_str(), 10);
  acb_set_fmpq(b.raw(), v, kernel::ComplexBall::precision());
  fmpq_clear(v);
  return b;
}
} // namespace affine_frobenius_detail

class AffineFrobeniusSeries {
public:
  using Matrix = fuchsify::Matrix;
  using B = kernel::ComplexBall;
  using Exponent = affine_frobenius_detail::Exponent;
  struct Options {
    unsigned max_dimension = 12, max_x_order = 128, max_terms = 250000,
             max_epsilon_depth = 512;
    // Clear polynomial row denominators in the constant Jordan frame.
    // False retains the original rational Taylor convolution for comparison.
    bool finite_lag_recurrence = true;
    unsigned max_clearing_degree = 32;
    // Disable only to compare the bounded cleared recurrence at small orders.
    bool finite_lag_cost_fallback = true;
    bool univariate_epsilon_recurrence = false;
    std::size_t max_projection_products=100000000;
    std::function<void(unsigned,unsigned)> column_progress;
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
    int low, high;
    std::vector<std::vector<std::vector<B>>> coefficients;
  };
  static constexpr bool omitted_tail_certified = false;
  unsigned dimension() const { return d_; }
  unsigned x_order() const { return n_; }
  unsigned finite_lag_rows() const { return finite_lag_rows_; }
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
    return prepare_impl(a,xi,ei,n,options,nullptr);
  }
private:
  struct ColumnCallbacks {
    std::function<std::optional<std::vector<Term>>(unsigned)> load;
    std::function<void(unsigned,std::span<const Term>)> save;
  };
  static AffineFrobeniusSeries prepare_impl(const Matrix &a,std::size_t xi,
      std::size_t ei,unsigned n,const Options &options,const ColumnCallbacks* columns) {
    using namespace affine_frobenius_detail;
    using namespace fuchsify::detail;
    if (a.empty() || a.size() > options.max_dimension ||
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
    // Constant Laurent column relations can mix exactly those base powers
    // that differ by integers (epsilon expansion makes slopes logarithms).
    // A common absolute cutoff commutes with these relations; a separate
    // number of Taylor coefficients per column does not.
    out.frontiers_.resize(d);
    for (unsigned c = 0; c < d; ++c) {
      auto minimum = out.eigen_[c].power;
      for (unsigned j = 0; j < d; ++j)
        if ((out.eigen_[j].power - out.eigen_[c].power).str().find('/') ==
            std::string::npos)
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
        // Denominator lags also contribute logarithmic derivatives. Require
        // at least a factor-two reduction in lag degree before selecting it.
        if (!n || (options.finite_lag_cost_fallback && lag_degree > (n-1)/2))
          cleared = false;
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
    auto generate = [&](const auto& coefficient_g,const auto& coefficient_lags,
        const auto& coefficient_frame,const auto& cz,const auto& cepsilon,auto to_exact) {
      using Coefficient=std::decay_t<decltype(cz)>;using CPoly=std::vector<Coefficient>;
    for (unsigned column = 0; column < d; ++column) {
      if(columns&&columns->load)if(auto restored=columns->load(column)) {
        if(restored->size()>options.max_terms-out.expansion_.terms.size())throw std::length_error("affine column checkpoint term budget exhausted");
        for(auto& term:*restored) {
          if(term.column!=column||term.row>=d||term.coefficient.is_zero())throw std::invalid_argument("affine column checkpoint coordinate mismatch");
          out.expansion_.terms.push_back(std::move(term));
        }
        if(options.column_progress)options.column_progress(column+1,d);
        continue;
      }
      const auto first_term=out.expansion_.terms.size();
      // Columns are independent. Orders beyond this column's absolute
      // frontier cannot contribute to any retained coefficient.
      const auto cutoff=(out.frontiers_[column]-out.eigen_[column].power).str();
      if(cutoff.find('/')!=std::string::npos)throw std::logic_error("nonintegral affine column frontier");
      const auto column_order=static_cast<unsigned>(std::stoul(cutoff));
      if(column_order>n)throw std::logic_error("affine column frontier exceeds prepared order");
      std::vector<std::vector<CPoly>> v(column_order + 1, std::vector<CPoly>(d));
      for (unsigned k = 0; k <= column_order; ++k)
        for (unsigned row = d; row-- > 0;) {
          CPoly rhs;
          if (k) {
            const auto &q = coefficient_lags[row];
            const unsigned last = q.empty() ? k : std::min(k, static_cast<unsigned>(q.size()-1));
            for (unsigned lag = 1; lag <= last; ++lag) {
              const bool derivative_lag = !q.empty() && !q[lag].is_zero();
              if (!derivative_lag) {
                for (unsigned j = 0; j < d; ++j)
                  add(rhs, v[k-lag][j], coefficient_g[lag][row][j]);
              } else {
                const auto alpha = cz.constant(out.eigen_[column].power + Rational(k-lag)) +
                    cz.constant(out.eigen_[column].slope) * cepsilon;
                // Cancel scalar weights before multiplying the large exact
                // previous-order coefficient, avoiding redundant gcd work.
                const auto diagonal = coefficient_g[lag][row][row] - q[lag] * alpha;
                for (unsigned j = 0; j < d; ++j)
                  add(rhs, v[k-lag][j], j == row ? diagonal : coefficient_g[lag][row][j]);
                CPoly derivative;
                const auto &previous = v[k-lag][row];
                for (unsigned l = 1; l < previous.size(); ++l)
                  derivative.push_back(previous[l] * cz.constant(l));
                add(rhs, derivative, -q[lag]);
              }
            }
          }
          if (out.successor_[row] >= 0)
            add(rhs, v[k][out.successor_[row]], cz.constant(1));
          auto delta =
              cz.constant(out.eigen_[column].power + Rational(k) -
                         out.eigen_[row].power) +
              cz.constant(out.eigen_[column].slope - out.eigen_[row].slope) *
                  cepsilon;
          v[k][row] = solve(rhs, delta);
          if (k == 0 && row == column)
            add(v[k][row], CPoly{cz.constant(1)}, cz.constant(1));
          if (v[k][row].size() > static_cast<std::size_t>(d) * (n + 1))
            throw std::logic_error(
                "affine Frobenius log ladder bound exceeded");
        }
      for (unsigned k = 0; k <= column_order; ++k)
        for (unsigned row = 0; row < d; ++row) {
          CPoly p;
          for (unsigned j = 0; j < d; ++j)
            add(p, v[k][j], coefficient_frame[row][j]);
          for (unsigned l = 0; l < p.size(); ++l)
            if (!p[l].is_zero() && out.eigen_[column].power + Rational(k) <=
                                       out.frontiers_[column])
              out.expansion_.terms.push_back(
                  {row, column, l, out.eigen_[column].power + Rational(k),
                   out.eigen_[column].slope, to_exact(p[l])});
          if (out.expansion_.terms.size() > options.max_terms)
            throw std::length_error("affine Frobenius term budget exhausted");
        }
      if(columns&&columns->save)columns->save(column,std::span<const Term>(out.expansion_.terms).subspan(first_term));
      if(options.column_progress)options.column_progress(column+1,d);
    }
    };
    if(options.univariate_epsilon_recurrence) {
      auto convert_matrix=[&](const auto& matrix){std::vector<std::vector<UnivariateRational>> result;
        for(const auto& row:matrix){result.emplace_back();for(const auto& value:row)result.back().emplace_back(value,ei);}return result;};
      std::vector<std::vector<std::vector<UnivariateRational>>> coefficient_g;
      for(const auto& matrix:g)coefficient_g.push_back(convert_matrix(matrix));
      generate(coefficient_g,convert_matrix(derivative_lags),convert_matrix(out.frame_),
        UnivariateRational(z,ei),UnivariateRational(z.variable(ei),ei),
        [&](const UnivariateRational& value){return value.exact(z,ei);});
    } else generate(g,derivative_lags,out.frame_,z,z.variable(ei),[](const Exact& value){return value;});
    return out;
  }
public:
  const Expansion &terms() const { return expansion_; }
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
  Expansion project_impl(const Matrix &rows,bool endpoint_only) const {
    using namespace affine_frobenius_detail;
    if(!options_.max_projection_products || rows.size()>options_.max_terms)
      throw std::length_error("affine projection receiving budget exhausted");
    Expansion out{static_cast<unsigned>(rows.size()), d_, {}};
    using Coordinate=std::tuple<unsigned,unsigned,Rational,Rational,unsigned>;
    std::map<Coordinate,Exact> accumulated;
    std::size_t products=0;
    out.coherent_x_frontier = !endpoint_only;
    if (!endpoint_only && rows.size() == d_ && expansion_.wronskian_prefactor)
      out.wronskian_prefactor = *expansion_.wronskian_prefactor *
                                affine_frobenius_detail::determinant(rows);
    for (unsigned i = 0; i < rows.size(); ++i) {
      if (rows[i].size() != d_)
        throw std::invalid_argument("affine observable row dimension");
      // Projection may have Laurent poles: only this shared row frontier
      // is determined by every input column through its known frontier.
      long row_valuation = std::numeric_limits<long>::max() / 4;
      for (const auto &coefficient : rows[i])
        if (!coefficient.is_zero())
          row_valuation = std::min(
              row_valuation, fuchsify::detail::valuation(coefficient, xi_));
      for (unsigned j = 0; j < d_; ++j) {
        if (rows[i][j].is_zero())
          continue;
        long val = fuchsify::detail::valuation(rows[i][j], xi_);
        if (std::abs(val) > options_.max_x_order)
          throw std::length_error(
              "observable Laurent order exceeds affine budget");
        auto regular =
            rows[i][j] / fuchsify::detail::power(zero_.variable(xi_), val);
        unsigned needed=n_;
        if(endpoint_only) {
          needed=0;
          for(const auto& term:expansion_.terms)if(term.row==j && term.slope.is_zero())
            for(unsigned k=0;k<=n_ && term.power+Rational(val)+Rational(k)<=Rational(0);++k)needed=std::max(needed,k);
        }
        auto coefficients = taylor(regular, xi_, needed);
        for (auto &t : expansion_.terms)
          if (t.row == j && (!endpoint_only || t.slope.is_zero()))
            for (unsigned n = 0; n <= needed; ++n) {
              if(endpoint_only && t.power+Rational(val)+Rational(n)>Rational(0))break;
              if (t.power + Rational(val) + Rational(n) >
                  frontiers_[t.column] + Rational(row_valuation))
                break;
              if (!coefficients[n].is_zero()) {
                if(products++>=options_.max_projection_products)
                  throw std::length_error("affine projection product budget exhausted");
                auto [at,inserted]=accumulated.try_emplace(
                  Coordinate{i,t.column,t.power+Rational(val)+Rational(n),t.slope,t.log_degree},zero_);
                at->second=at->second+t.coefficient*coefficients[n];
                if(accumulated.size()>options_.max_terms)
                  throw std::length_error("affine projection accumulated coordinate budget exhausted");
              }
            }
      }
    }
    for(auto& [coordinate,coefficient]:accumulated)if(!coefficient.is_zero()) {
      const auto& [row,column,power,slope,log]=coordinate;
      out.terms.push_back({row,column,log,power,slope,std::move(coefficient)});
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
    for (auto &t : out.terms) {
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
             t.power.str().find('/') != std::string::npos;
    });
    return out;
  }
  Matrix dr_endpoint_constant(const Expansion &input) const {
    auto canonical = normalize(input);
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
    auto canonical=normalize(input);DRDomain out{canonical,{}};
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
    for (auto &t : canonical.terms)
      if (t.slope.is_zero() && t.power <= Rational(-1))
        throw std::domain_error(
            "unregulated fixed nonintegrable endpoint power in DR integral");
    auto out = primitive(canonical);
    auto lower = dr_endpoint_constant(out);
    for (unsigned i = 0; i < out.rows; ++i)
      for (unsigned j = 0; j < out.columns; ++j)
        if (!lower[i][j].is_zero())
          out.terms.push_back(
              {i, j, 0, Rational(0), Rational(0), -lower[i][j]});
    return normalize(std::move(out));
  }
  Valuations valuation_metadata(const Expansion &input) const {
    Valuations out;
    out.minimum_by_column.assign(input.columns, Valuations::zero_valuation);
    for (auto &t : normalize(input).terms) {
      auto val = fuchsify::detail::valuation(t.coefficient, ei_);
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
      long valuation = fuchsify::detail::valuation(t.coefficient, ei_);
      if (valuation > high)
        continue;
      long length = static_cast<long>(high) - valuation + 1;
      if (length > options_.max_epsilon_depth)
        throw std::length_error(
            "affine coefficient poles require additional epsilon-depth budget");
      auto regular = t.coefficient /
                     fuchsify::detail::power(zero_.variable(ei_), valuation);
      Jet ep(0, length, B::precision());
      if (length > 1)
        ep.set(1, B(1));
      auto value = diffexp::evaluate(data::Reader(regular.str()).read(), ep,
                                      {{zero_.variables()[ei_], ep}});
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
  friend struct cached_affine::StateAccess;
  Exact zero_;
  unsigned d_ = 0, n_ = 0, finite_lag_rows_ = 0;
  std::size_t xi_ = 0, ei_ = 0;
  Options options_;
  Matrix frame_;
  std::vector<Exponent> eigen_;
  std::vector<Rational> frontiers_;
  std::vector<int> successor_;
  Expansion expansion_{0, 0, {}};
  explicit AffineFrobeniusSeries(Exact z) : zero_(std::move(z)) {}
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
      for (std::size_t v = 0; v < t.coefficient.variable_count(); ++v)
        if (v != ei_ && !t.coefficient.derivative(v).is_zero())
          throw std::invalid_argument(
              "affine term coefficient must depend only on epsilon");
      auto key =
          std::make_tuple(t.row, t.column, t.power, t.slope, t.log_degree);
      auto [p, inserted] = combined.try_emplace(key, zero_);
      p->second = p->second + t.coefficient;
    }
    out.terms.clear();
    for (auto &[key, c] : combined)
      if (!c.is_zero()) {
        auto [i, j, a, b, l] = key;
        out.terms.push_back({i, j, l, a, b, c});
      }
    if (out.terms.size() > options_.max_terms)
      throw std::length_error("affine expansion term budget exhausted");
    return out;
  }
};
} // namespace diffexp
