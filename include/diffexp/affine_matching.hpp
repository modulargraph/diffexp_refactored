#pragma once
#include "diffexp/affine_frobenius.hpp"
#include "diffexp/deepest_beta.hpp"
#include <bit>
#include <flint/acb_mat.h>

namespace diffexp::affine_matching {
using B = kernel::ComplexBall;
struct Window {
  int low, high;
};
// Coefficients below low are declared structural zeros. Missing coefficients
// above high are unknown and are never padded with zeros.
struct Boundary {
  int low = 0, high = -1;
  std::vector<std::vector<B>> coefficients;
};
enum class Status {
  Success,
  NeedMoreBoundary,
  NeedMoreSource,
  NeedMorePrecision,
  NeedMoreXOrder,
  Unsupported
};
struct Options {
  unsigned max_dimension = 6, max_determinant_dimension = 6,
           max_determinant_depth = 8, max_epsilon_depth = 64;
  unsigned max_row_normalization_steps = 8, max_symbolic_x_power = 4096;
  std::size_t max_exact_operations = 20000000, max_exact_monomials = 100000,
              max_series_cells = 8000000;
  // Operator-only alternative; exact normalization and K0 remain mandatory.
  bool numeric_operator_convolution = false;
};
struct Result {
  Status status = Status::Unsupported;
  std::string reason;
  Boundary value;
  int required_boundary_high = 0, required_source_high = 0;
  long determinant_valuation = 0, expected_determinant_valuation = 0;
  bool wronskian_verified = false;
  bool omitted_tail_certified = false;
  bool used_regular_leading = false;
  bool used_row_normalization = false;
  unsigned normalization_steps = 0, row_pole_loss = 0;
  bool success() const { return status == Status::Success; }
};
inline Boundary from_laurent(const std::vector<feynman::LaurentValue> &input) {
  if (input.empty())
    throw std::invalid_argument("empty Laurent boundary");
  int low = input[0].epsilon_low, high = std::numeric_limits<int>::max();
  for (const auto &v : input) {
    if (v.coefficients.empty())
      throw std::invalid_argument("empty Laurent boundary row");
    low = std::min(low, v.epsilon_low);
    const long row_high = static_cast<long>(v.epsilon_low) +
                          static_cast<long>(v.coefficients.size()) - 1;
    if (row_high > std::numeric_limits<int>::max())
      throw std::length_error("Laurent boundary upper order overflow");
    high = std::min(high, static_cast<int>(row_high));
  }
  if (high < low)
    throw std::invalid_argument("Laurent boundary has no common window");
  Boundary out{low, high,
               std::vector(input.size(), std::vector<B>(high - low + 1, B(0)))};
  for (unsigned i = 0; i < input.size(); ++i)
    for (int k = std::max(low, input[i].epsilon_low); k <= high; ++k)
      out.coefficients[i][k - low] = input[i].at(k);
  return out;
}
namespace detail {
using Q = Rational;
using Polynomial = std::map<std::pair<Q, unsigned>, Q>;
using Series = std::vector<Polynomial>;
using Matrix = std::vector<std::vector<Series>>;
struct Budget {
  Options options;
  std::size_t work = 0;
  void series_cells(std::size_t rows, std::size_t columns,
                    std::size_t depth) const {
    auto limit = options.max_series_cells;
    if (!limit || (columns && rows > limit / columns) ||
        (depth && rows * columns > limit / depth))
      throw std::length_error("affine matching series-cell budget exhausted");
  }
  void spend(std::size_t n = 1) {
    if (n > options.max_exact_operations - work)
      throw std::length_error(
          "affine matching exact-operation budget exhausted");
    work += n;
  }
};
inline void add(Polynomial &out, const Polynomial &p, const Q &factor = Q(1)) {
  for (auto &[key, q] : p) {
    auto [it, inserted] = out.try_emplace(key, Q(0));
    it->second += q * factor;
    if (it->second.is_zero())
      out.erase(it);
  }
}
inline Polynomial multiply(const Polynomial &a, const Polynomial &b,
                           Budget &budget) {
  Polynomial out;
  for (auto &[ka, qa] : a)
    for (auto &[kb, qb] : b) {
      budget.spend();
      auto key = std::make_pair(ka.first + kb.first, ka.second + kb.second);
      auto [it, inserted] = out.try_emplace(key, Q(0));
      it->second += qa * qb;
      if (it->second.is_zero())
        out.erase(it);
    }
  if (out.size() > budget.options.max_exact_monomials)
    throw std::length_error(
        "affine matching symbolic monomial budget exhausted");
  return out;
}
inline Series multiply(const Series &a, const Series &b, unsigned top,
                       Budget &budget) {
  Series out(top + 1);
  for (unsigned i = 0; i < a.size() && i <= top; ++i)
    if (!a[i].empty())
      for (unsigned j = 0; j < b.size() && i + j <= top; ++j)
        if (!b[j].empty())
          add(out[i + j], multiply(a[i], b[j], budget));
  return out;
}
inline void add(Series &out, const Series &p, int sign) {
  for (unsigned n = 0; n < p.size(); ++n)
    add(out[n], p[n], Q(sign));
}
inline Series one(unsigned top) {
  Series out(top + 1);
  out[0][{Q(0), 0}] = Q(1);
  return out;
}
inline Series determinant(const Matrix &a, unsigned top, Budget &budget) {
  unsigned d = a.size();
  if (!d)
    return one(top);
  if (d > 12)
    throw std::length_error(
        "affine determinant subset fallback is limited to dimension twelve");
  budget.series_cells(std::size_t(1) << d, 1,
                      static_cast<std::size_t>(top) + 1);
  std::vector<Series> dp(std::size_t(1) << d, Series(top + 1));
  dp[0] = one(top);
  for (std::size_t mask = 0; mask + 1 < dp.size(); ++mask) {
    unsigned row = std::popcount(mask);
    for (unsigned column = 0; column < d; ++column)
      if (!(mask & (std::size_t(1) << column))) {
        int sign = std::popcount(mask >> (column + 1)) % 2 ? -1 : 1;
        auto product = multiply(dp[mask], a[row][column], top, budget);
        add(dp[mask | (std::size_t(1) << column)], product, sign);
      }
  }
  return dp.back();
}
inline Matrix minor(const Matrix &a, unsigned row, unsigned column) {
  Matrix out;
  for (unsigned i = 0; i < a.size(); ++i)
    if (i != row) {
      std::vector<Series> r;
      for (unsigned j = 0; j < a.size(); ++j)
        if (j != column)
          r.push_back(a[i][j]);
      out.push_back(std::move(r));
    }
  return out;
}
inline std::pair<long, std::vector<Q>> rational_coefficients(const Exact &c,
                                                             unsigned top) {
  if (c.is_rational()) {
    std::vector<Q> out(top + 1);
    out[0] = c.rational();
    return {0, out};
  }
  std::optional<std::size_t> variable;
  for (auto &ts : {c.numerator_terms(), c.denominator_terms()})
    for (auto &t : ts)
      for (unsigned j = 0; j < t.powers.size(); ++j)
        if (t.powers[j]) {
          if (variable && *variable != j)
            throw std::invalid_argument(
                "matching coefficient depends on multiple parameters");
          variable = j;
        }
  if (!variable)
    throw std::logic_error("nonconstant matching coefficient has no variable");
  auto val = fuchsify::detail::valuation(c, *variable);
  auto regular = c / fuchsify::detail::power(c.variable(*variable), val);
  auto exact = affine_frobenius_detail::taylor(regular, *variable, top);
  std::vector<Q> out;
  for (auto &e : exact)
    out.push_back(e.rational());
  return {val, out};
}
// Expand epsilon in an exact ring of x^r log(x)^k. Cancellation is decided
// here, before evaluating at the matchpoint; no interval is classified as exact
// zero.
inline Matrix coefficients(const AffineFrobeniusSeries::Expansion &expansion,
                           const std::vector<long> &shifts, unsigned top,
                           Budget &budget) {
  unsigned d = expansion.rows;
  budget.series_cells(d, d, static_cast<std::size_t>(top) + 1);
  Matrix out(d, std::vector<Series>(d, Series(top + 1)));
  for (auto &t : expansion.terms) {
    auto [val, c] = rational_coefficients(t.coefficient, top);
    long start = val - shifts[t.column];
    if (start < 0)
      throw std::logic_error(
          "column normalization left a negative epsilon order");
    if (start > top)
      continue;
    for (unsigned k = 0; k + start <= top; ++k)
      if (!c[k].is_zero()) {
        Q slope_power(1);
        for (unsigned m = 0; k + m + start <= top; ++m) {
          budget.spend();
          auto key = std::make_pair(t.power, t.log_degree + m);
          auto &p = out[t.row][t.column][k + m + start];
          auto [it, inserted] = p.try_emplace(key, Q(0));
          it->second += c[k] * slope_power;
          if (it->second.is_zero())
            p.erase(it);
          slope_power = slope_power * t.slope / Q(m + 1);
          if (slope_power.is_zero())
            break;
        }
      }
  }
  return out;
}
inline B evaluate(const Polynomial &p, const B &x) {
  B logx;
  acb_log(logx.raw(), x.raw(), B::precision());
  B out(0);
  for (auto &[key, q] : p) {
    B xp;
    auto exponent = affine_frobenius_detail::ball(key.first) * logx;
    acb_exp(xp.raw(), exponent.raw(), B::precision());
    for (unsigned l = 0; l < key.second; ++l)
      xp = xp * logx;
    out = out + affine_frobenius_detail::ball(q) * xp;
  }
  return out;
}
inline void validate(const Boundary &b, unsigned rows) {
  if (b.coefficients.size() != rows || b.low > b.high)
    throw std::invalid_argument("Laurent boundary shape/window mismatch");
  for (auto &r : b.coefficients) {
    if (r.size() !=
        static_cast<std::size_t>(static_cast<long>(b.high) - b.low + 1))
      throw std::invalid_argument("Laurent boundary row/window mismatch");
    for (auto &v : r)
      if (!v.is_finite())
        throw std::invalid_argument("nonfinite Laurent boundary coefficient");
  }
}
inline std::optional<std::vector<std::vector<B>>>
verified_inverse(const Matrix &matrix, const B &point) {
  const unsigned d = matrix.size();
  acb_mat_t a, inv;
  acb_mat_init(a, d, d);
  acb_mat_init(inv, d, d);
  for (unsigned i = 0; i < d; ++i)
    for (unsigned j = 0; j < d; ++j) {
      auto value = evaluate(matrix[i][j][0], point);
      acb_set(acb_mat_entry(a, i, j), value.raw());
    }
  const bool success = acb_mat_inv(inv, a, B::precision()) != 0;
  std::optional<std::vector<std::vector<B>>> out;
  if (success) {
    out = std::vector(d, std::vector<B>(d, B(0)));
    for (unsigned i = 0; i < d; ++i)
      for (unsigned j = 0; j < d; ++j)
        acb_set((*out)[i][j].raw(), acb_mat_entry(inv, i, j));
  }
  acb_mat_clear(a);
  acb_mat_clear(inv);
  return out;
}
inline int integer(long n) {
  if (n < std::numeric_limits<int>::min() ||
      n > std::numeric_limits<int>::max())
    throw std::length_error("affine matching epsilon window overflow");
  return n;
}
// A verified invertible leading epsilon matrix avoids determinant expansion.
// This path scales to larger master spaces; it never infers a pivot from size.
inline Result regular_match(const AffineFrobeniusSeries::Expansion &physical,
                            const std::vector<long> &shifts, const B &point,
                            const Boundary &boundary, Window requested,
                            const std::vector<std::vector<B>> &inverse,
                            Budget &budget) {
  Result out;
  out.used_regular_leading = true;
  for (long v : shifts)
    out.determinant_valuation += v;
  unsigned d = physical.rows;
  long max_shift = *std::max_element(shifts.begin(), shifts.end());
  long required = static_cast<long>(requested.high) + max_shift;
  out.required_boundary_high = integer(required);
  if (required > boundary.high) {
    out.status = Status::NeedMoreBoundary;
    out.reason = "matching requires additional upper boundary coefficients; "
                 "missing coefficients are unknown";
    return out;
  }
  long depth = std::max(0L, required - boundary.low);
  if (depth >= budget.options.max_epsilon_depth)
    throw std::length_error(
        "regular affine matching epsilon-depth budget exhausted");
  auto exact = coefficients(physical, shifts, depth, budget);
  std::vector<std::vector<std::vector<B>>> h(
      depth + 1, std::vector(d, std::vector<B>(d, B(0))));
  for (unsigned n = 1; n <= depth; ++n)
    for (unsigned i = 0; i < d; ++i)
      for (unsigned j = 0; j < d; ++j)
        h[n][i][j] = evaluate(exact[i][j][n], point);
  std::vector<std::vector<B>> w(depth + 1, std::vector<B>(d, B(0)));
  for (unsigned n = 0; n <= depth; ++n) {
    std::vector<B> rhs(d, B(0));
    for (unsigned i = 0; i < d; ++i) {
      if (static_cast<long>(boundary.low) + n > boundary.high)
        throw std::logic_error(
            "regular matching attempted an unavailable boundary coefficient");
      rhs[i] = boundary.coefficients[i][n];
      for (unsigned lag = 1; lag <= n; ++lag)
        for (unsigned j = 0; j < d; ++j)
          rhs[i] = rhs[i] - h[lag][i][j] * w[n - lag][j];
    }
    for (unsigned i = 0; i < d; ++i)
      for (unsigned j = 0; j < d; ++j)
        w[n][i] = w[n][i] + inverse[i][j] * rhs[j];
  }
  int low = integer(std::min<long>(
      requested.low, static_cast<long>(boundary.low) - max_shift));
  out.value = {
      low, requested.high,
      std::vector(d, std::vector<B>(static_cast<long>(requested.high) - low + 1,
                                    B(0)))};
  for (unsigned i = 0; i < d; ++i)
    for (int k = low; k <= requested.high; ++k) {
      long order = static_cast<long>(k) + shifts[i] - boundary.low;
      if (order >= 0)
        out.value.coefficients[i][k - low] = w.at(order)[i];
    }
  out.status = Status::Success;
  return out;
}

struct PrecisionNeeded : std::runtime_error {
  using std::runtime_error::runtime_error;
};
using ExactMatrix = fuchsify::Matrix;
using ExactSeries = std::vector<ExactMatrix>;
// Native sparse polynomial evaluation avoids constructing/copying a deeply
// nested text-expression tree for each normalized matrix coefficient.
class ExactPointEvaluator {
  B point_, logarithm_;
  std::map<unsigned long, B> x_powers_{{0, B(1)}}, log_powers_{{0, B(1)}};
  const B &power(bool logarithm, unsigned long exponent) {
    auto &cache = logarithm ? log_powers_ : x_powers_;
    auto found = cache.find(exponent);
    if (found != cache.end())
      return found->second;
    B value;
    acb_pow_ui(value.raw(), logarithm ? logarithm_.raw() : point_.raw(),
               exponent, B::precision());
    return cache.emplace(exponent, std::move(value)).first->second;
  }

public:
  explicit ExactPointEvaluator(const B &point) : point_(point) {
    acb_log(logarithm_.raw(), point.raw(), B::precision());
  }
  B operator()(const Exact &value) {
    if (value.is_zero())
      return B(0);
    const auto polynomial = [&](const std::vector<Exact::Term> &terms) {
      B sum(0);
      for (const auto &term : terms) {
        auto coefficient = affine_frobenius_detail::ball(term.coefficient);
        for (unsigned i = 0; i < term.powers.size(); ++i)
          if (term.powers[i]) {
            const auto &name = value.variables()[i];
            if (name != "x" && name != "ell")
              throw std::invalid_argument(
                  "unsubstituted exact matching coefficient parameter");
            coefficient = coefficient * power(name == "ell", term.powers[i]);
          }
        sum = sum + coefficient;
      }
      return sum;
    };
    auto denominator = polynomial(value.denominator_terms());
    if (!denominator.is_finite() || denominator.contains_zero())
      throw PrecisionNeeded(
          "a pointwise exact row-normalization denominator contains zero; "
          "increase precision or change matchpoint");
    auto result = polynomial(value.numerator_terms()) / denominator;
    if (!result.is_finite())
      throw PrecisionNeeded(
          "a pointwise exact row-normalization coefficient is nonfinite");
    return result;
  }
};
inline B evaluate_exact(const Exact &value, const B &point) {
  return ExactPointEvaluator(point)(value);
}
inline Exact as_exact(const Polynomial &p, const Q &column_power,
                      const Exact &zero, Budget &budget) {
  auto out = zero, x = zero.variable(0), ell = zero.variable(1);
  for (auto &[key, q] : p) {
    budget.spend();
    auto exponent = key.first - column_power;
    auto text = exponent.str();
    if (text.find('/') != std::string::npos)
      throw std::domain_error("physical affine column does not have integer x "
                              "shifts after its base exponent is factored");
    long power = std::stol(text);
    if (std::abs(power) > budget.options.max_symbolic_x_power)
      throw std::length_error(
          "row normalization symbolic x-power budget exhausted");
    out = out + zero.constant(q) * fuchsify::detail::power(x, power) *
                    ell.pow(key.second);
  }
  return out;
}
inline ExactSeries
exact_coefficients(const AffineFrobeniusSeries &series,
                   const AffineFrobeniusSeries::Expansion &physical,
                   const std::vector<long> &shifts, unsigned top,
                   const Exact &zero, Budget &budget) {
  auto symbolic = coefficients(physical, shifts, top, budget);
  unsigned d = physical.rows;
  ExactSeries out(top + 1, fuchsify::detail::zeros(d, d, zero));
  for (unsigned n = 0; n <= top; ++n)
    for (unsigned i = 0; i < d; ++i)
      for (unsigned j = 0; j < d; ++j)
        out[n][i][j] = as_exact(symbolic[i][j][n], series.exponents()[j].power,
                                zero, budget);
  return out;
}
inline ExactMatrix exact_product(const ExactMatrix &a, const ExactMatrix &b,
                                 Budget &budget) {
  auto out = fuchsify::detail::zeros(a.size(), b[0].size(), a[0][0]);
  for (unsigned i = 0; i < a.size(); ++i)
    for (unsigned k = 0; k < b.size(); ++k)
      if (!a[i][k].is_zero())
        for (unsigned j = 0; j < b[0].size(); ++j)
          if (!b[k][j].is_zero()) {
            budget.spend();
            out[i][j] = out[i][j] + a[i][k] * b[k][j];
          }
  return out;
}
struct RowGauge {
  ExactSeries inverse_epsilon_powers;
  unsigned steps = 0, sheared_rows = 0;
};
template <class Coefficient>
inline RowGauge row_gauge_coefficients(unsigned d, const Exact &zero,
                                       const B &point, Budget &budget,
                                       Coefficient &&coefficient) {
  ExactPointEvaluator evaluate_at(point);
  RowGauge gauge{{fuchsify::detail::identity(d, zero)}};
  for (unsigned step = 0;; ++step) {
    // Reconstruct only the leading coefficient; propagating each Gaussian
    // pivot through all epsilon orders would cause unnecessary rational swell.
    auto leading = fuchsify::detail::zeros(d, d, zero);
    for (unsigned p = 0; p < gauge.inverse_epsilon_powers.size(); ++p) {
      auto product = exact_product(gauge.inverse_epsilon_powers[p],
                                   coefficient(p), budget);
      for (unsigned i = 0; i < d; ++i)
        for (unsigned j = 0; j < d; ++j)
          leading[i][j] = leading[i][j] + product[i][j];
    }
    unsigned rank = 0;
    std::vector<bool> used(d, false);
    for (;;) {
      std::optional<std::pair<unsigned, unsigned>> pivot;
      bool structural_nonzero = false;
      for (unsigned column = 0; column < d && !pivot; ++column)
        if (!used[column])
          for (unsigned row = rank; row < d; ++row)
            if (!leading[row][column].is_zero()) {
              structural_nonzero = true;
              try {
                auto value = evaluate_at(leading[row][column]);
                if (value.is_finite() && !value.contains_zero()) {
                  pivot = {{row, column}};
                  break;
                }
              } catch (const PrecisionNeeded &) {
              }
            }
      if (!pivot) {
        if (structural_nonzero)
          throw PrecisionNeeded(
              "structurally nonzero epsilon-leading row pivots cannot be "
              "separated from zero at the matchpoint; increase precision or "
              "change point");
        break;
      }
      auto [row, column] = *pivot;
      used[column] = true;
      if (row != rank) {
        std::swap(leading[row], leading[rank]);
        for (auto &m : gauge.inverse_epsilon_powers)
          std::swap(m[row], m[rank]);
      }
      const auto divisor = leading[rank][column];
      for (auto &v : leading[rank])
        if (!v.is_zero()) {
          budget.spend();
          v = v / divisor;
        }
      for (auto &m : gauge.inverse_epsilon_powers)
        for (auto &v : m[rank])
          if (!v.is_zero()) {
            budget.spend();
            v = v / divisor;
          }
      for (unsigned i = 0; i < d; ++i)
        if (i != rank && !leading[i][column].is_zero()) {
          const auto factor = leading[i][column];
          for (unsigned j = 0; j < d; ++j)
            if (!leading[rank][j].is_zero()) {
              budget.spend();
              leading[i][j] = leading[i][j] - factor * leading[rank][j];
            }
          for (auto &m : gauge.inverse_epsilon_powers)
            for (unsigned j = 0; j < d; ++j)
              if (!m[rank][j].is_zero()) {
                budget.spend();
                m[i][j] = m[i][j] - factor * m[rank][j];
              }
        }
      if (++rank == d)
        break;
    }
    if (rank == d)
      return gauge;
    if (step >= budget.options.max_row_normalization_steps)
      throw std::length_error(
          "singular-leading epsilon row-normalization step budget exhausted");
    // Only exact zero leading rows are shifted. Numerical smallness is never
    // a license to divide a nonzero row by epsilon.
    for (unsigned i = rank; i < d; ++i)
      for (auto &v : leading[i])
        if (!v.is_zero())
          throw std::logic_error(
              "epsilon shear attempted on a structurally nonzero leading row");
    ExactSeries next(gauge.inverse_epsilon_powers.size() + 1,
                     fuchsify::detail::zeros(d, d, zero));
    for (unsigned p = 0; p < gauge.inverse_epsilon_powers.size(); ++p)
      for (unsigned i = 0; i < d; ++i)
        next[p + (i >= rank)][i] = gauge.inverse_epsilon_powers[p][i];
    while (next.size() > 1 && fuchsify::detail::zero(next.back()))
      next.pop_back();
    gauge.inverse_epsilon_powers = std::move(next);
    ++gauge.steps;
    gauge.sheared_rows += d - rank;
  }
}
inline RowGauge row_gauge(const ExactSeries &h, const B &point,
                          Budget &budget) {
  return row_gauge_coefficients(
      h.at(0).size(), h.at(0).at(0).at(0).constant(0), point, budget,
      [&](unsigned n) -> const ExactMatrix & { return h.at(n); });
}
// Epsilon orders are generated only when an exact zero leading row requires
// another shear. The configured maximum is a bound, never an eager demand.
inline RowGauge
adaptive_row_gauge(const AffineFrobeniusSeries &series,
                   const AffineFrobeniusSeries::Expansion &physical,
                   const std::vector<long> &shifts, const B &point,
                   const Exact &zero, Budget &budget) {
  ExactSeries coefficients;
  return row_gauge_coefficients(
      physical.rows, zero, point, budget,
      [&](unsigned n) -> const ExactMatrix & {
        if (n >= coefficients.size())
          coefficients =
              exact_coefficients(series, physical, shifts, n, zero, budget);
        return coefficients[n];
      });
}
inline Result
row_normalized_match(const AffineFrobeniusSeries &series,
                     const AffineFrobeniusSeries::Expansion &physical,
                     const std::vector<long> &shifts, const B &point,
                     const Boundary &boundary, Window requested,
                     Budget &budget) {
  if (!budget.options.max_row_normalization_steps)
    throw std::length_error("singular-leading row normalization is disabled");
  if (budget.options.max_row_normalization_steps >=
      budget.options.max_epsilon_depth)
    throw std::length_error(
        "row-normalization probe exceeds the epsilon-depth budget");
  ExactField field({"x", "ell"});
  Exact zero(field, 0);
  unsigned d = physical.rows;
  auto gauge =
      adaptive_row_gauge(series, physical, shifts, point, zero, budget);
  unsigned loss = gauge.inverse_epsilon_powers.size() - 1;
  long max_shift = *std::max_element(shifts.begin(), shifts.end());
  Result out;
  out.used_row_normalization = true;
  out.normalization_steps = gauge.steps;
  out.row_pole_loss = loss;
  out.determinant_valuation = gauge.sheared_rows;
  for (long v : shifts)
    out.determinant_valuation += v;
  const long w_top = static_cast<long>(requested.high) + max_shift,
             required = w_top + loss;
  out.required_boundary_high = integer(required);
  if (required > boundary.high) {
    out.status = Status::NeedMoreBoundary;
    out.reason = "exact epsilon row normalization requires additional boundary "
                 "upper coefficients";
    return out;
  }
  const long transformed_low = static_cast<long>(boundary.low) - loss,
             depth = std::max(0L, w_top - transformed_low);
  if (depth + loss >= budget.options.max_epsilon_depth)
    throw std::length_error(
        "normalized affine matching epsilon-depth budget exhausted");
  auto original =
      exact_coefficients(series, physical, shifts, depth + loss, zero, budget);
  // Certify that the Laurent row gauge introduces no negative epsilon terms
  // in the normalized matrix, with identities checked exactly in Q(x,ell).
  for (unsigned q = 1; q <= loss; ++q) {
    auto residual = fuchsify::detail::zeros(d, d, zero);
    for (unsigned p = q; p <= loss; ++p) {
      auto product = exact_product(gauge.inverse_epsilon_powers[p],
                                   original[p - q], budget);
      for (unsigned i = 0; i < d; ++i)
        for (unsigned j = 0; j < d; ++j)
          residual[i][j] = residual[i][j] + product[i][j];
    }
    if (!fuchsify::detail::zero(residual))
      throw std::logic_error(
          "epsilon row-gauge negative-order certificate failed");
  }
  ExactSeries normalized(depth + 1, fuchsify::detail::zeros(d, d, zero));
  for (unsigned n = 0; n <= depth; ++n)
    for (unsigned p = 0; p <= loss; ++p) {
      auto product = exact_product(gauge.inverse_epsilon_powers[p],
                                   original[n + p], budget);
      for (unsigned i = 0; i < d; ++i)
        for (unsigned j = 0; j < d; ++j)
          normalized[n][i][j] = normalized[n][i][j] + product[i][j];
    }
  std::vector<std::vector<std::vector<B>>> k(
      depth + 1, std::vector(d, std::vector<B>(d, B(0)))),
      r(loss + 1, std::vector(d, std::vector<B>(d, B(0))));
  ExactPointEvaluator evaluate_at(point);
  for (unsigned n = 0; n <= depth; ++n)
    for (unsigned i = 0; i < d; ++i)
      for (unsigned j = 0; j < d; ++j)
        k[n][i][j] = evaluate_at(normalized[n][i][j]);
  for (unsigned p = 0; p <= loss; ++p)
    for (unsigned i = 0; i < d; ++i)
      for (unsigned j = 0; j < d; ++j)
        r[p][i][j] = evaluate_at(gauge.inverse_epsilon_powers[p][i][j]);
  acb_mat_t leading, inverse;
  acb_mat_init(leading, d, d);
  acb_mat_init(inverse, d, d);
  for (unsigned i = 0; i < d; ++i)
    for (unsigned j = 0; j < d; ++j)
      acb_set(acb_mat_entry(leading, i, j), k[0][i][j].raw());
  bool ok = acb_mat_inv(inverse, leading, B::precision());
  std::vector<std::vector<B>> inv(d, std::vector<B>(d, B(0)));
  if (ok)
    for (unsigned i = 0; i < d; ++i)
      for (unsigned j = 0; j < d; ++j)
        acb_set(inv[i][j].raw(), acb_mat_entry(inverse, i, j));
  acb_mat_clear(leading);
  acb_mat_clear(inverse);
  if (!ok)
    throw PrecisionNeeded("normalized unit-leading matrix is not invertible at "
                          "working precision");
  std::vector<std::vector<B>> w(depth + 1, std::vector<B>(d, B(0)));
  for (unsigned n = 0; n <= depth; ++n) {
    std::vector<B> rhs(d, B(0));
    for (unsigned p = 0; p <= loss; ++p) {
      long order = transformed_low + n + p;
      if (order < boundary.low)
        continue;
      if (order > boundary.high)
        throw std::logic_error("epsilon normalization attempted an unknown "
                               "upper boundary coefficient");
      for (unsigned i = 0; i < d; ++i)
        for (unsigned j = 0; j < d; ++j)
          rhs[i] = rhs[i] +
                   r[p][i][j] * boundary.coefficients[j][order - boundary.low];
    }
    for (unsigned lag = 1; lag <= n; ++lag)
      for (unsigned i = 0; i < d; ++i)
        for (unsigned j = 0; j < d; ++j)
          rhs[i] = rhs[i] - k[lag][i][j] * w[n - lag][j];
    for (unsigned i = 0; i < d; ++i)
      for (unsigned j = 0; j < d; ++j)
        w[n][i] = w[n][i] + inv[i][j] * rhs[j];
  }
  int low = integer(std::min<long>(requested.low, transformed_low - max_shift));
  out.value = {
      low, requested.high,
      std::vector(d, std::vector<B>(static_cast<long>(requested.high) - low + 1,
                                    B(0)))};
  B logarithm;
  acb_log(logarithm.raw(), point.raw(), B::precision());
  for (unsigned i = 0; i < d; ++i) {
    B factor;
    auto exponent =
        -affine_frobenius_detail::ball(series.exponents()[i].power) * logarithm;
    acb_exp(factor.raw(), exponent.raw(), B::precision());
    for (int order = low; order <= requested.high; ++order) {
      long n = static_cast<long>(order) + shifts[i] - transformed_low;
      if (n >= 0)
        out.value.coefficients[i][order - low] = factor * w.at(n)[i];
    }
  }
  out.status = Status::Success;
  return out;
}
} // namespace detail

inline Result match(const AffineFrobeniusSeries &series,
                    const AffineFrobeniusSeries::Expansion &physical,
                    const B &point, const Boundary &boundary, Window requested,
                    const Options &options = {}) {
  using namespace detail;
  Result result;
  try {
    unsigned d = physical.rows;
    if (!d || d != physical.columns || d != series.dimension() ||
        d > options.max_dimension || options.max_dimension > 256 ||
        !options.max_exact_operations || !options.max_epsilon_depth ||
        requested.low > requested.high ||
        static_cast<long>(requested.high) - requested.low + 1 >
            options.max_epsilon_depth ||
        point.contains_zero())
      throw std::invalid_argument(
          "invalid square affine matching request or finite budgets");
    validate(boundary, d);
    auto guard = series.wronskian_valuation(physical, point);
    if (!guard || !physical.coherent_x_frontier)
      throw std::domain_error(
          "affine matching requires coherent absolute x frontiers and an "
          "analytic-epsilon Wronskian guard");
    if (!guard->second.is_finite() || guard->second.contains_zero()) {
      result.status = Status::NeedMorePrecision;
      result.reason = "full-system Wronskian prefactor is not separated from "
                      "zero at the matchpoint";
      return result;
    }
    const auto audit = [&](Result r) {
      r.expected_determinant_valuation = guard->first;
      if (r.status == Status::Success || r.status == Status::NeedMoreBoundary) {
        if (r.determinant_valuation != guard->first) {
          r.status = Status::NeedMoreXOrder;
          r.reason = "retained determinant epsilon valuation " +
                     std::to_string(r.determinant_valuation) +
                     " disagrees with full-system Wronskian valuation " +
                     std::to_string(guard->first) +
                     "; increase coherent absolute x order";
          r.value = {};
        } else
          r.wronskian_verified = true;
      }
      return r;
    };
    Budget budget{options};
    auto valuations = series.valuation_metadata(physical);
    auto shifts = valuations.minimum_by_column;
    for (auto v : shifts)
      if (v == AffineFrobeniusSeries::Valuations::zero_valuation)
        throw std::domain_error(
            "affine matching frame has an exact zero column");
    auto leading_matrix = coefficients(physical, shifts, 0, budget);
    if (auto inverse = verified_inverse(leading_matrix, point))
      return audit(regular_match(physical, shifts, point, boundary, requested,
                                 *inverse, budget));
    if (d > options.max_determinant_dimension || d > 12)
      return audit(row_normalized_match(series, physical, shifts, point,
                                        boundary, requested, budget));
    unsigned probe =
        std::min(options.max_determinant_depth, options.max_epsilon_depth - 1);
    auto matrix = coefficients(physical, shifts, probe, budget);
    auto det = determinant(matrix, probe, budget);
    std::optional<unsigned> order;
    B leading;
    for (unsigned n = 0; n <= probe; ++n) {
      if (det[n].empty())
        continue;
      auto value = evaluate(det[n], point);
      if (value.is_zero())
        continue;
      if (!value.is_finite() || value.contains_zero()) {
        result.status = Status::NeedMorePrecision;
        result.reason =
            "the exact determinant coefficient cannot be separated from zero "
            "at this matchpoint; increase precision or change point";
        return result;
      }
      order = n;
      leading = value;
      break;
    }
    if (!order)
      throw std::domain_error(
          "normalized determinant valuation exceeds the bounded search, or the "
          "retained frame is singular at this point");
    const long m = *order,
               max_shift = *std::max_element(shifts.begin(), shifts.end());
    result.determinant_valuation = m;
    for (auto v : shifts)
      result.determinant_valuation += v;
    result.status = Status::Success;
    result = audit(std::move(result));
    if (result.status == Status::NeedMoreXOrder)
      return result;
    const long required = static_cast<long>(requested.high) + m + max_shift;
    result.required_boundary_high = integer(required);
    if (required > boundary.high) {
      result.status = Status::NeedMoreBoundary;
      result.reason = "matching requires additional upper boundary "
                      "coefficients; missing coefficients are unknown";
      return result;
    }
    const long depth = std::max(0L, required - boundary.low),
               matrix_depth = m + depth;
    if (matrix_depth >= options.max_epsilon_depth)
      throw std::length_error(
          "affine inverse requires a larger epsilon-depth budget");
    if (matrix_depth != probe) {
      matrix = coefficients(physical, shifts, matrix_depth, budget);
      det = determinant(matrix, matrix_depth, budget);
    }
    std::vector<B> inverse(depth + 1, B(0)), denominator(depth + 1, B(0));
    for (unsigned n = 0; n <= depth; ++n)
      denominator[n] = evaluate(det[n + m], point);
    if (denominator[0].contains_zero()) {
      result.status = Status::NeedMorePrecision;
      result.reason = "normalized determinant lost numerical separation";
      return result;
    }
    inverse[0] = B(1) / denominator[0];
    for (unsigned n = 1; n <= depth; ++n) {
      B v(0);
      for (unsigned k = 1; k <= n; ++k)
        v = v + denominator[k] * inverse[n - k];
      inverse[n] = -v / denominator[0];
    }
    const int low = integer(std::min<long>(
        requested.low, static_cast<long>(boundary.low) - m - max_shift));
    result.value = {
        low, requested.high,
        std::vector(d, std::vector<B>(
                           static_cast<long>(requested.high) - low + 1, B(0)))};
    for (unsigned i = 0; i < d; ++i) {
      std::vector<B> numerator(depth + 1, B(0));
      for (unsigned j = 0; j < d; ++j) {
        auto cofactor = determinant(minor(matrix, j, i), depth, budget);
        for (unsigned n = 0; n <= depth; ++n) {
          B value = evaluate(cofactor[n], point);
          if ((i + j) % 2)
            value = -value;
          for (long k = boundary.low;
               k <= boundary.high && n + k - boundary.low <= depth; ++k)
            numerator[n + k - boundary.low] =
                numerator[n + k - boundary.low] +
                value * boundary.coefficients[j][k - boundary.low];
        }
      }
      for (int k = low; k <= requested.high; ++k) {
        long degree = static_cast<long>(k) + m + shifts[i] - boundary.low;
        if (degree < 0)
          continue;
        if (degree > depth)
          throw std::logic_error("affine inverse depth accounting failed");
        B value(0);
        for (long n = 0; n <= degree; ++n)
          value = value + numerator[n] * inverse[degree - n];
        result.value.coefficients[i][k - low] = value;
      }
    }
    result.status = Status::Success;
  } catch (const detail::PrecisionNeeded &e) {
    result.status = Status::NeedMorePrecision;
    result.reason = e.what();
  } catch (const std::exception &e) {
    result.status = Status::Unsupported;
    result.reason = e.what();
  }
  return result;
}
inline Result match(const AffineFrobeniusSeries &series, const B &point,
                    const Boundary &boundary, Window requested,
                    const Options &options = {}) {
  return match(series, series.terms(), point, boundary, requested, options);
}
inline Result match(const AffineFrobeniusSeries &series, const B &point,
                    const std::vector<feynman::LaurentValue> &boundary,
                    Window requested, const Options &options = {}) {
  return match(series, point, from_laurent(boundary), requested, options);
}

// Apply a projected primitive (or any affine expansion) to matched constants.
// Its exact epsilon valuations determine the required source depth first.
inline Result apply(const AffineFrobeniusSeries &series,
                    const AffineFrobeniusSeries::Expansion &functional,
                    const B &point, const Boundary &source, Window requested,
                    const Options &options = {}) {
  Result result;
  try {
    detail::validate(source, functional.columns);
    if (requested.low > requested.high ||
        static_cast<long>(requested.high) - requested.low + 1 >
            options.max_epsilon_depth)
      throw std::invalid_argument("invalid affine functional output window");
    auto valuations = series.valuation_metadata(functional);
    long required = std::numeric_limits<long>::min(),
         minimum = AffineFrobeniusSeries::Valuations::zero_valuation;
    for (long v : valuations.minimum_by_column)
      if (v != AffineFrobeniusSeries::Valuations::zero_valuation) {
        required = std::max(required, static_cast<long>(requested.high) - v);
        minimum = std::min(minimum, v);
      }
    if (minimum == AffineFrobeniusSeries::Valuations::zero_valuation ||
        static_cast<long>(requested.high) <
            static_cast<long>(source.low) + minimum) {
      result.value = {
          requested.low, requested.high,
          std::vector(functional.rows,
                      std::vector<B>(static_cast<long>(requested.high) -
                                         requested.low + 1,
                                     B(0)))};
      result.status = Status::Success;
      return result;
    }
    result.required_source_high = detail::integer(required);
    if (required > source.high) {
      result.status = Status::NeedMoreSource;
      result.reason = "the projected functional's epsilon poles require "
                      "additional matching-constant coefficients";
      return result;
    }
    int low =
        detail::integer(std::min<long>(requested.low, source.low + minimum));
    auto matrix = series.evaluate_laurent(
        functional, point, detail::integer(minimum),
        detail::integer(static_cast<long>(requested.high) - source.low));
    result.value = {
        low, requested.high,
        std::vector(
            functional.rows,
            std::vector<B>(static_cast<long>(requested.high) - low + 1, B(0)))};
    for (unsigned i = 0; i < functional.rows; ++i)
      for (unsigned j = 0; j < functional.columns; ++j)
        for (int k = low; k <= requested.high; ++k)
          for (long e = minimum; e <= matrix.high; ++e) {
            long c = static_cast<long>(k) - e;
            if (c < source.low || c > source.high)
              continue;
            result.value.coefficients[i][k - low] =
                result.value.coefficients[i][k - low] +
                matrix.coefficients[i][j][e - matrix.low] *
                    source.coefficients[j][c - source.low];
          }
    result.status = Status::Success;
  } catch (const std::exception &e) {
    result.status = Status::Unsupported;
    result.reason = e.what();
  }
  return result;
}
} // namespace diffexp::affine_matching
