// Endpoint series in truncated Laurent arithmetic against the exact rational
// functions of epsilon: every coefficient, every projection and every DR
// functional must be the exact one modulo the stated epsilon order.
#include "diffexp/affine_operator.hpp"
#include <iostream>
#include <map>
#include <tuple>
using namespace diffexp;
using B = kernel::ComplexBall;
using Series = AffineFrobeniusSeries;
using S = affine_frobenius_detail::EpsilonSeries;
using Key = std::tuple<unsigned, unsigned, Rational, Rational, unsigned>;
static void check(bool v, const std::string &why) {
  if (!v)
    throw std::runtime_error(why);
}
// `truncated` is `exact` modulo its stated orders; at least `orders` of them
// are known in every entry that `require` names. Returns the terms compared.
static std::size_t same(const Series::Expansion &exact, const Series::Expansion &truncated,
                        const std::string &what, long minimum = std::numeric_limits<long>::min()) {
  check(exact.rows == truncated.rows && exact.columns == truncated.columns, what + ": shape");
  check(truncated.epsilon_truncation.size() == truncated.rows, what + ": no truncation orders");
  check(exact.row_x_valuations == truncated.row_x_valuations, what + ": x valuations");
  check(exact.coherent_x_frontier == truncated.coherent_x_frontier, what + ": frontier flag");
  std::map<Key, Exact> left, right;
  for (const auto &t : exact.terms)
    left.emplace(Key{t.row, t.column, t.power, t.slope, t.log_degree}, t.coefficient);
  for (const auto &t : truncated.terms) {
    check(!t.coefficient.is_zero(), what + ": a stored zero");
    right.emplace(Key{t.row, t.column, t.power, t.slope, t.log_degree}, t.coefficient);
  }
  std::size_t compared = 0;
  for (const auto &[key, value] : left) {
    const long order = truncated.truncation(std::get<0>(key), std::get<1>(key));
    check(order >= minimum, what + ": fewer epsilon orders than asked for");
    const auto expected = order >= Series::Expansion::exact_order
                              ? value : S(value, 1, order).exact(value, 1);
    const auto found = right.find(key);
    if (expected.is_zero())
      check(found == right.end(), what + ": a coefficient beyond the truncation order is stored");
    else {
      check(found != right.end(), what + ": a coefficient is missing");
      check(found->second == expected, what + ": a coefficient differs from the exact one");
      ++compared;
    }
  }
  for (const auto &[key, value] : right)
    check(left.count(key), what + ": a coefficient that the exact series does not have");
  return compared;
}
static Series::Matrix matrix(const ExactField &field, const std::vector<std::vector<const char *>> &text) {
  Series::Matrix out;
  for (const auto &row : text) {
    out.emplace_back();
    for (const auto *entry : row)
      out.back().emplace_back(field, entry);
  }
  return out;
}
int main() {
  try {
    B::set_precision(256);
    ExactField field({"x", "eps"});
    auto e = [&](const char *s) { return Exact(field, s); };
    {
      // The series type itself: orders of sums, products and quotients.
      const S a(e("(1+eps)/(eps^2*(1-2*eps))"), 1, 3), b(e("eps^3*(2+eps)/(3-eps)"), 1, 6);
      check(a.valuation() == -2 && a.order() == 3 && b.valuation() == 3 && b.order() == 6, "series orders");
      const auto product = a * b, quotient = a / b, sum = a + b;
      check(product.order() == 4 && product.valuation() == 1, "product order");
      check(product.exact(e("0"), 1) == S(e("(1+eps)/(eps^2*(1-2*eps))") * e("eps^3*(2+eps)/(3-eps)"), 1, 4).exact(e("0"), 1),
            "product coefficients");
      check(quotient.order() == -2 && quotient.valuation() == -5, "quotient order");
      check(quotient.exact(e("0"), 1) == S(e("(1+eps)/(eps^2*(1-2*eps))") / e("eps^3*(2+eps)/(3-eps)"), 1, -2).exact(e("0"), 1),
            "quotient coefficients");
      check(sum.order() == 3 && sum.exact(e("0"), 1) == S(e("(1+eps)/(eps^2*(1-2*eps))"), 1, 3).exact(e("0"), 1), "sum");
      const S zero(e("eps^5"), 1, 4);
      check(zero.is_zero() && !zero.is_exact() && zero.valuation() == 4 && (zero * a).order() == 2, "zero to its order");
      check(S().is_exact() && (S() * a).is_exact() && (S() * a).is_zero(), "exact zero");
      bool refused = false;
      try { (void)(a / zero); } catch (const std::domain_error &) { refused = true; }
      check(refused, "division by a series that is zero to its order");
      check(S::expand(e("3/eps^2+eps"), 1, 1).is_exact() && !S::expand(e("1/(1-eps)"), 1, 4).is_exact() &&
                S::expand(e("1/(eps*(1-eps))"), 1, 4).order() == 3, "expansion of a rational function");
    }
    struct Case { const char *name; Series::Matrix a; unsigned order; };
    const std::vector<Case> cases{
        {"triangular", matrix(field, {{"eps/x", "1/x"}, {"0", "eps/x+(x+eps)/(1-x)"}}), 8},
        {"fixed", matrix(field, {{"0", "0"}, {"1/(1-x)", "1/x"}}), 8},
        {"jordan", matrix(field, {{"(-1+2*eps)/x", "1/x"}, {"0", "(-1+2*eps)/x"}}), 4},
        {"colliding frame", matrix(field, {{"0", "1/x"}, {"0", "eps/x"}}), 4},
        // A fixed integer resonance: an exact logarithm.
        {"logarithm", matrix(field, {{"eps/x", "0"}, {"1", "(1+eps)/x"}}), 5},
        // A resonance of different slopes: a pole in epsilon.
        {"pole", matrix(field, {{"0", "0"}, {"1", "(1+eps)/x"}}), 5},
        {"rational jordan", matrix(field, {{"eps/x", "1/((1-eps)*x)"}, {"0", "eps/x"}}), 4},
        {"half integers", matrix(field, {{"(1/2+eps)/x", "1/(1-x)"}, {"x/(1+x)", "(-1/2-eps)/x+1/(2-x)"}}), 7},
        // Exponents eps, 1-2*eps, 2+eps and a fixed 1: poles at two orders,
        // logarithms where the fixed exponents meet, and moving row scales.
        {"mixed", matrix(field, {{"eps/x", "1/(1-x)", "x/(2-x)", "1"},
                                 {"(1+eps)/(1+x)", "(1-2*eps)/x", "eps/(1-x)", "0"},
                                 {"1/(3-eps)", "1", "(2+eps)/x", "1/(1-2*x)"},
                                 {"1/(1+x)", "eps", "1/(1-eps)", "1/x+1/(3+x)"}}), 9},
        // Three integer exponents without epsilon: logarithms up to degree two.
        {"logarithms", matrix(field, {{"0", "0", "0"}, {"1/(1-x)", "1/x", "0"}, {"eps", "1+x", "2/x+1/(1+x)"}}), 7}};
    std::size_t compared = 0, logarithmic = 0, poles = 0;
    for (const auto &item : cases)
      for (const bool finite_lag : {true, false})
        for (const unsigned orders : {1u, 4u, 9u}) {
          Series::Options exact_options, options;
          exact_options.finite_lag_recurrence = options.finite_lag_recurrence = finite_lag;
          options.epsilon_orders = orders;
          options.threads = exact_options.threads = 1;
          const auto exact = Series::prepare(item.a, 0, 1, item.order, exact_options);
          const auto series = Series::prepare(item.a, 0, 1, item.order, options);
          const std::string what = std::string(item.name) + (finite_lag ? "" : " (Taylor lags)") + " at " + std::to_string(orders) + " orders";
          check(series.truncated() && !exact.truncated(), what + ": mode");
          check(series.exponents().size() == exact.exponents().size(), what + ": exponents");
          // Every column is known through `orders` orders from the lowest
          // order of its frame column.
          for (unsigned c = 0; c < series.dimension(); ++c) {
            long lead = Series::Expansion::exact_order;
            for (unsigned i = 0; i < series.dimension(); ++i)
              if (!series.residue_frame()[i][c].is_zero())
                lead = std::min(lead, fuchsify::detail::valuation(series.residue_frame()[i][c], 1));
            for (unsigned i = 0; i < series.dimension(); ++i)
              check(series.terms().truncation(i, c) >= lead + static_cast<long>(orders), what + ": a column is short of its orders");
          }
          compared += same(exact.terms(), series.terms(), what);
          for (const auto &t : series.terms().terms) {
            logarithmic += t.log_degree != 0;
            poles += series.epsilon_valuation(t.coefficient) < 0;
          }
        }
    check(compared > 1000 && logarithmic > 50 && poles > 50, "the truncated cases lost their logarithms or poles");
    std::cout << "Truncated epsilon series: " << compared << " coefficients agree with the exact rational functions ("
              << logarithmic << " logarithmic, " << poles << " with poles)\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
