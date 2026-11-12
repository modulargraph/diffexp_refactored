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
// The constraints of a truncated domain are those of the exact one modulo
// their orders; one whose coefficients all vanish that far may be absent.
static void same_constraints(const Series::DRDomain &exact, const Series::DRDomain &truncated, const std::string &what) {
  std::size_t matched = 0;
  for (const auto &expected : exact.zero_constraints) {
    const Series::ZeroConstraint *found = nullptr;
    for (const auto &candidate : truncated.zero_constraints)
      if (candidate.row == expected.row && candidate.power == expected.power && candidate.log_degree == expected.log_degree)
        found = &candidate;
    if (found) {
      ++matched;
      check(found->truncation.size() == expected.coefficients.size(), what + ": constraint without orders");
    }
    for (unsigned j = 0; j < expected.coefficients.size(); ++j) {
      const long order = found ? found->truncation[j] : truncated.admissible.truncation(expected.row, j);
      const auto value = order >= Series::Expansion::exact_order || expected.coefficients[j].is_zero()
                             ? expected.coefficients[j] : S(expected.coefficients[j], 1, order).exact(expected.coefficients[j], 1);
      check(found ? found->coefficients[j] == value : value.is_zero(), what + ": constraint coefficient differs");
    }
  }
  check(matched == truncated.zero_constraints.size(), what + ": a constraint that the exact domain does not have");
}
static bool overlap(const Series::LaurentMatrix &a, const Series::LaurentMatrix &b, int high) {
  if (a.coefficients.size() != b.coefficients.size())
    return false;
  for (unsigned i = 0; i < a.coefficients.size(); ++i)
    for (unsigned j = 0; j < a.coefficients[i].size(); ++j)
      for (int k = std::min(a.low, b.low); k <= high; ++k) {
        const B left = k >= a.low && k <= a.high ? a.coefficients[i][j][k - a.low] : B(0);
        const B right = k >= b.low && k <= b.high ? b.coefficients[i][j][k - b.low] : B(0);
        if (!acb_overlaps(left.raw(), right.raw()))
          return false;
      }
  return true;
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
    // Projections, endpoint domains and DR primitives: rows with poles in x,
    // with poles in epsilon, and with a pole that moves with epsilon (whose
    // Taylor coefficients in x lose an epsilon order each).
    const std::vector<const char *> entries{"(1+eps)/(x^2*(1+x))", "(x+eps)/(3*x*(1-x))", "1/(2-eps)", "0",
                                            "1/(eps*(1+2*x))", "x/(x+eps)", "(1-x)/(1+eps*x)", "eps^2/x"};
    std::size_t projected = 0, primitives = 0, refused = 0, operators = 0, demands = 0;
    for (const auto &item : cases)
      for (const unsigned orders : {2u, 6u}) {
        Series::Options options;
        options.epsilon_orders = orders;
        const auto exact = Series::prepare(item.a, 0, 1, item.order);
        const auto series = Series::prepare(item.a, 0, 1, item.order, options);
        const unsigned d = series.dimension();
        const std::string what = std::string(item.name) + " at " + std::to_string(orders) + " orders";
        for (unsigned variant = 0; variant < 4; ++variant) {
          Series::Matrix rows(2 + variant % 2, std::vector<Exact>(d, e("0")));
          for (unsigned i = 0; i < rows.size(); ++i)
            for (unsigned j = 0; j < d; ++j)
              rows[i][j] = e(entries[(3 * variant + 5 * i + j + (variant > 1 ? 4 : 0)) % entries.size()]);
          const auto full = exact.project(rows), cut = series.project(rows);
          projected += same(full, cut, what + ": projection");
          // Limits: the fixed nonpositive sector, its constant and its domain.
          Series::DRDomain limit_exact, limit_cut;
          bool reaches = true;
          try { limit_exact = exact.project_endpoint_domain(rows); } catch (const std::domain_error &) { reaches = false; }
          if (reaches) {
            limit_cut = series.project_endpoint_domain(rows);
            same(limit_exact.admissible, limit_cut.admissible, what + ": endpoint domain");
            same_constraints(limit_exact, limit_cut, what + ": endpoint domain");
            same(exact.dr_endpoint_functional(limit_exact.admissible), series.dr_endpoint_functional(limit_cut.admissible), what + ": endpoint constant");
          } else {
            bool also = false;
            try { series.project_endpoint_domain(rows); } catch (const std::domain_error &) { also = true; }
            check(also, what + ": the truncated endpoint domain reached further than the exact one");
          }
          // Integrals: the admissible part and its DR primitive.
          bool integrable = true;
          Series::DRDomain domain_exact;
          Series::Expansion primitive_exact;
          try {
            domain_exact = exact.dr_domain(full, true);
            primitive_exact = exact.dr_integral_from_zero(domain_exact.admissible);
          } catch (const std::domain_error &) { integrable = false; }
          if (integrable) {
            const auto domain_cut = series.dr_domain(cut, true);
            same_constraints(domain_exact, domain_cut, what + ": integral domain");
            const auto primitive_cut = series.dr_integral_from_zero(domain_cut.admissible);
            // The primitive drops the frontier metadata on both paths.
            primitives += same(primitive_exact, primitive_cut, what + ": primitive");
            const auto point = B::from_strings("0.25");
            long lowest = 0, known = Series::Expansion::exact_order;
            for (const auto &t : primitive_cut.terms)
              lowest = std::min(lowest, series.epsilon_valuation(t.coefficient));
            for (const auto &row : primitive_cut.epsilon_truncation)
              for (const auto order : row)
                known = std::min(known, order);
            if (known < Series::Expansion::exact_order && known - 1 >= lowest) {
              check(overlap(exact.evaluate_laurent(primitive_exact, point, lowest, known - 1),
                            series.evaluate_laurent(primitive_cut, point, lowest, known - 1), known - 1),
                    what + ": Laurent evaluation of the primitive");
              bool asked = false;
              try { series.evaluate_laurent(primitive_cut, point, lowest, known + 2); }
              catch (const EpsilonOrdersNeeded &demand) { asked = demand.missing == 3; }
              check(asked, what + ": an evaluation beyond the known orders did not ask for them");
            }
          } else {
            bool also = false;
            try { series.dr_integral_from_zero(series.dr_domain(cut, true).admissible); } catch (const std::domain_error &) { also = true; }
            check(also, what + ": the truncated integral accepted what the exact one refuses");
            ++refused;
          }
        }
        // The endpoint operator: the inverse of the frame and a functional
        // composed with it. Too few orders must be asked for, never guessed.
        if (!exact.terms().wronskian_prefactor)
          continue;
        const auto point = B::from_strings("0.125");
        const int high = 2;
        const auto reference = affine_operator::prepare(exact, exact.terms(), point, high);
        if (!reference.success())
          continue;
        Series::Matrix row(1, std::vector<Exact>(d, e("0")));
        for (unsigned j = 0; j < d; ++j)
          row[0][j] = e(entries[(2 * j + 1) % 3]);
        Series::DRDomain reference_domain;
        try { reference_domain = exact.project_endpoint_domain(row); } catch (const std::domain_error &) { continue; }
        const auto reference_operator = affine_operator::compose(reference, exact, exact.dr_endpoint_functional(reference_domain.admissible), point, high);
        auto current = options;
        for (unsigned attempt = 0;; ++attempt) {
          check(attempt < 6, what + ": the demand for epsilon orders does not settle");
          const auto refined = Series::prepare(item.a, 0, 1, item.order, current);
          try {
            const auto prepared = affine_operator::prepare(refined, refined.terms(), point, high);
            check(prepared.success(), what + ": truncated inverse: " + prepared.reason);
            check(prepared.row_lower_bounds == reference.row_lower_bounds && prepared.determinant_valuation == reference.determinant_valuation &&
                      prepared.row_pole_loss == reference.row_pole_loss, what + ": the structure of the truncated inverse differs");
            check(overlap(reference.inverse, prepared.inverse, high), what + ": truncated inverse differs from the exact one");
            const auto functional = refined.dr_endpoint_functional(refined.project_endpoint_domain(row).admissible);
            const auto composed = affine_operator::compose(prepared, refined, functional, point, high);
            check(composed.success() == reference_operator.success(), what + ": truncated operator: " + composed.reason);
            if (composed.success()) {
              check(composed.row_lower_bounds == reference_operator.row_lower_bounds, what + ": operator bounds differ");
              check(overlap(reference_operator.matrix, composed.matrix, high), what + ": truncated operator differs from the exact one");
            }
            ++operators;
            break;
          } catch (const EpsilonOrdersNeeded &demand) {
            check(demand.missing > 0, what + ": an empty demand");
            current.epsilon_orders += static_cast<unsigned>(demand.missing);
            ++demands;
          }
        }
      }
    {
      // A fixed 1/x has no DR primitive, truncated or not.
      Series::Options options;
      options.epsilon_orders = 3;
      const auto constant = Series::prepare({{e("0")}}, 0, 1, 2, options);
      try { constant.dr_integral_from_zero(constant.project(std::vector<Exact>{e("1/x")})); }
      catch (const std::domain_error &) { ++refused; }
      const auto scalar = Series::prepare({{e("eps/x")}}, 0, 1, 3, options);
      const auto integral = scalar.dr_integral_from_zero(scalar.project(std::vector<Exact>{e("1/x")}));
      check(integral.terms.size() == 1 && integral.terms[0].coefficient == e("1/eps") && integral.truncation(0, 0) == 2,
            "the truncated DR x^(eps-1) integral lost its pole or its order");
    }
    check(projected > 2000 && primitives > 500 && refused > 0 && operators > 4 && demands > 2,
          "the truncated projections lost their coverage");
    {
      // Nilpotent coefficients whose kernels have echelon bases with poles at
      // eps=0 (eps*v2+v3=0 gives (0,-1/eps,1)): Jordan chains that are regular
      // with a regular inverse exist, and are found.
      using fuchsify::detail::multiply;
      for (const auto &n : {matrix(field, {{"0", "eps", "1"}, {"0", "0", "0"}, {"0", "0", "0"}}),
                            matrix(field, {{"0", "eps", "1", "0"}, {"0", "0", "0", "1+eps"}, {"0", "0", "0", "-eps-eps^2"}, {"0", "0", "0", "0"}}),
                            matrix(field, {{"0", "1", "eps", "2"}, {"0", "0", "1", "1"}, {"0", "0", "0", "1-eps"}, {"0", "0", "0", "0"}})}) {
        const auto plain = fuchsify::detail::nilpotent_frame(n), regular = fuchsify::detail::nilpotent_frame(n, 1);
        check(!fuchsify::detail::unimodular(plain, 1) || plain == regular, "the echelon Jordan frame was expected to be singular at eps=0");
        check(fuchsify::detail::unimodular(regular, 1), "no regular Jordan chains were found for a matrix of dimension " + std::to_string(n.size()) + " with entry " + n[0][1].str());
        const auto jordan = multiply(fuchsify::detail::inverse(regular), multiply(n, regular));
        for (unsigned i = 0; i < jordan.size(); ++i)
          for (unsigned j = 0; j < jordan.size(); ++j)
            check(jordan[i][j].is_zero() || (j == i + 1 && jordan[i][j] == e("1")), "regular chains are not Jordan chains");
      }
      // A connection with a pole at eps=0 is answered with exact coefficients.
      Series::Options options;
      options.epsilon_orders = 5;
      check(!Series::prepare(matrix(field, {{"eps/x", "1/(eps*(1-x))"}, {"0", "(1+eps)/x"}}), 0, 1, 4, options).truncated(),
            "a connection with poles at eps=0 was truncated");
    }
    std::cout << "Truncated projections: " << projected << " projected and " << primitives << " integrated coefficients agree, "
              << refused << " refusals shared, " << operators << " operators within the exact enclosures after " << demands << " demands\n";
    std::cout << "Truncated epsilon series: " << compared << " coefficients agree with the exact rational functions ("
              << logarithmic << " logarithmic, " << poles << " with poles)\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
