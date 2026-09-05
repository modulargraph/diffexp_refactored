#include "diffexp/affine_matching.hpp"
#include "diffexp/families.hpp"
#include "diffexp/frobenius.hpp"
#include "diffexp/level_cache.hpp"
#include "diffexp/level_preparation.hpp"
#include "diffexp/recursion_graph.hpp"
#include <iostream>
using namespace diffexp;
using B = kernel::ComplexBall;
namespace am = affine_matching;
void check(bool v, const std::string &why) {
  if (!v)
    throw std::runtime_error(why);
}
void near(const B &a, const B &b) {
  auto e = a - b;
  arf_t upper;
  arf_init(upper);
  acb_get_abs_ubound_arf(upper, e.raw(), B::precision());
  double error = arf_get_d(upper, ARF_RND_CEIL);
  arf_clear(upper);
  check(error < 1e-35, "Laurent matching analytic oracle mismatch");
}
B at(const am::Boundary &b, unsigned row, int n) {
  check(n >= b.low && n <= b.high, "test coefficient outside window");
  return b.coefficients[row][n - b.low];
}
int main(int argc, char **argv) {
  try {
    B::set_precision(256);
    ExactField field({"x", "eps"});
    auto e = [&](const char *s) { return Exact(field, s); };
    auto half = B::from_strings("0.5");
    B logx;
    acb_log(logx.raw(), half.raw(), B::precision());
    auto scalar = AffineFrobeniusSeries::prepare({{e("eps/x")}}, 0, 1, 2);
    auto basic = am::match(
        scalar, B(1), am::Boundary{-2, 2, {{B(3), B(2), B(1), B(4), B(5)}}},
        {-2, 2});
    check(basic.success(), basic.reason);
    near(at(basic.value, 0, -2), B(3));
    near(at(basic.value, 0, 2), B(5));
    auto collision = AffineFrobeniusSeries::prepare(
        {{e("0"), e("1/x")}, {e("0"), e("eps/x")}}, 0, 1, 2);
    auto missing = am::match(collision, B(1),
                             am::Boundary{0, 0, {{B(0)}, {B(1)}}}, {-1, 0});
    check(missing.status == am::Status::NeedMoreBoundary &&
              missing.required_boundary_high == 1,
          "missing boundary upper order must be demanded explicitly");
    am::Boundary analytic{0, 4, std::vector(2, std::vector<B>(5, B(0)))};
    B power(1);
    for (unsigned n = 0; n <= 4; ++n) {
      analytic.coefficients[1][n] = power;
      power = power * logx / B(n + 1);
      analytic.coefficients[0][n] = power;
    }
    auto matched = am::match(collision, half, analytic, {-1, 2});
    check(matched.success(), matched.reason);
    near(at(matched.value, 0, -1), B(-1));
    near(at(matched.value, 1, 0), B(1));
    for (int k = 0; k <= 2; ++k)
      near(at(matched.value, 0, k), B(0));
    auto pseudo = AffineFrobeniusSeries::prepare(
        {{e("0"), e("0")}, {e("1"), e("(1+eps)/x")}}, 0, 1, 2);
    auto pb = analytic;
    for (unsigned n = 0; n <= 4; ++n) {
      pb.coefficients[0][n] = B(n == 0 ? 1 : 0);
      pb.coefficients[1][n] = half * analytic.coefficients[0][n];
    }
    auto pm = am::match(pseudo, half, pb, {-1, 2});
    check(pm.success(), pm.reason);
    near(at(pm.value, 0, 0), B(1));
    near(at(pm.value, 1, -1), B(1));
    // Match the physical expansion directly, including an epsilon-singular
    // gauge.
    auto constant = AffineFrobeniusSeries::prepare(
        {{e("0"), e("0")}, {e("0"), e("0")}}, 0, 1, 2);
    auto physical = constant.project(
        AffineFrobeniusSeries::Matrix{{e("1"), e("x/eps")}, {e("0"), e("1")}});
    auto gm =
        am::match(constant, physical, half,
                  am::Boundary{0, 1, {{B(2), B(0)}, {B(1), B(0)}}}, {-1, 0});
    check(gm.success(), gm.reason);
    near(at(gm.value, 0, -1), -half);
    near(at(gm.value, 0, 0), B(2));
    near(at(gm.value, 1, 0), B(1));
    auto primitive = scalar.dr_integral_from_zero(
        scalar.project(std::vector<Exact>{e("1/x")}));
    auto need = am::apply(scalar, primitive, half, am::Boundary{0, 0, {{B(1)}}},
                          {-1, 0});
    check(need.status == am::Status::NeedMoreSource &&
              need.required_source_high == 1,
          "primitive pole must request matching-constant lookahead");
    auto integrated = am::apply(scalar, primitive, half,
                                am::Boundary{0, 1, {{B(1), B(2)}}}, {-1, 0});
    check(integrated.success(), integrated.reason);
    near(at(integrated.value, 0, -1), B(1));
    near(at(integrated.value, 0, 0), logx + B(2));
    auto vanishing = scalar.project(std::vector<Exact>{e("eps^2")});
    auto below_support = am::apply(scalar, vanishing, half,
                                   am::Boundary{0, 0, {{B(1)}}}, {0, 0});
    check(below_support.success(), below_support.reason);
    near(at(below_support.value, 0, 0), B(0));
    // An uncertain determinant coefficient must never become a nonzero pivot.
    auto trivial = AffineFrobeniusSeries::prepare({{e("0")}}, 0, 1, 2);
    auto zero_at_point = trivial.project(std::vector<Exact>{e("x-1/2")});
    auto ambiguous = am::match(trivial, zero_at_point, half,
                               am::Boundary{0, 0, {{B(1)}}}, {0, 0});
    check(ambiguous.status == am::Status::NeedMorePrecision ||
              ambiguous.status == am::Status::Unsupported,
          "zero-containing determinant coefficient must not be inverted");
    // Independent epsilon-regular Frobenius matcher on a resonant log system.
    AffineFrobeniusSeries::Matrix connection{{e("eps/x"), e("0")},
                                             {e("1"), e("(1+eps)/x")}};
    auto affine = AffineFrobeniusSeries::prepare(connection, 0, 1, 3);
    auto regular = FrobeniusSeries::prepare(connection, 0, 1, 3, 3);
    FrobeniusSeries::Boundary rb{{B(2), B(0), B(0), B(0)},
                                 {B(3), B(0), B(0), B(0)}};
    auto rc = regular.match(half, rb);
    auto ac = am::match(affine, half, am::Boundary{0, 3, rb}, {0, 3});
    check(ac.success(), ac.reason);
    for (unsigned i = 0; i < 2; ++i)
      for (int k = 0; k <= 3; ++k)
        near(at(ac.value, i, k), rc[i][k]);
    // The regular-leading path must support master spaces larger than the
    // bounded determinant fallback, including a physical matrix gauge.
    constexpr unsigned large_dimension = 15;
    AffineFrobeniusSeries::Matrix large_matrix(
        large_dimension, std::vector<Exact>(large_dimension, e("0"))),
        large_gauge = large_matrix;
    for (unsigned i = 0; i < large_dimension; ++i) {
      large_matrix[i][i] = e("eps/x");
      large_gauge[i][i] = e("1");
      if (i + 1 < large_dimension)
        large_gauge[i][i + 1] = e("x");
    }
    AffineFrobeniusSeries::Options large_series_options;
    large_series_options.max_dimension = 32;
    auto large = AffineFrobeniusSeries::prepare(large_matrix, 0, 1, 1,
                                                large_series_options);
    am::Options large_options;
    large_options.max_dimension = 32;
    am::Boundary large_boundary{
        0, 2, std::vector(large_dimension, std::vector<B>{B(1), B(0), B(0)})};
    auto large_match = am::match(large, large.project(large_gauge), half,
                                 large_boundary, {0, 2}, large_options);
    check(large_match.success() && large_match.used_regular_leading,
          large_match.reason);
    B geometric(0);
    for (unsigned i = large_dimension; i-- > 0;) {
      geometric = B(1) - half * geometric;
      near(at(large_match.value, i, 0), geometric);
      near(at(large_match.value, i, 1), -geometric * logx);
    }
    auto singular_gauge = large_gauge;
    singular_gauge[0][0] = e("1");
    singular_gauge[0][1] = e("1");
    singular_gauge[1][1] = e("eps");
    auto normalized = am::match(large, large.project(singular_gauge), half,
                                large_boundary, {0, 0}, large_options);
    check(normalized.success() && normalized.used_row_normalization,
          normalized.reason);
    check(normalized.row_pole_loss == 1, "single epsilon-rank loss metadata");
    std::vector<B> geometric_values(large_dimension, B(0));
    for (unsigned i = large_dimension; i-- > 0;)
      geometric_values[i] =
          B(1) -
          (i + 1 < large_dimension ? half * geometric_values[i + 1] : B(0));
    near(at(normalized.value, 1, -1), geometric_values[1]);
    near(at(normalized.value, 0, -1), -geometric_values[1]);
    near(at(normalized.value, 1, 0), -geometric_values[1] * logx);
    auto shallow = large_boundary;
    shallow.high = 0;
    for (auto &row : shallow.coefficients)
      row.resize(1);
    auto require_guard = am::match(large, large.project(singular_gauge), half,
                                   shallow, {0, 0}, large_options);
    check(require_guard.status == am::Status::NeedMoreBoundary &&
              require_guard.required_boundary_high == 1,
          "row normalization must demand source lookahead");
    singular_gauge[1][1] = e("eps^2");
    auto twice = am::match(large, large.project(singular_gauge), half,
                           large_boundary, {0, 0}, large_options);
    check(twice.success() && twice.row_pole_loss == 2, twice.reason);
    near(at(twice.value, 1, -2), geometric_values[1]);
    near(at(twice.value, 1, 0), geometric_values[1] * logx * logx / B(2));
    auto restricted = large_options;
    restricted.max_row_normalization_steps = 1;
    auto incomplete = am::match(large, large.project(singular_gauge), half,
                                large_boundary, {0, 0}, restricted);
    check(incomplete.status == am::Status::Unsupported &&
              incomplete.reason.find("step budget") != std::string::npos,
          "epsilon row-normalization budget must be enforced");
    auto rank_one = large_gauge;
    for (auto &row : rank_one)
      for (auto &entry : row)
        entry = e("0");
    rank_one[0][0] = e("1");
    for (unsigned j = 1; j < large_dimension; ++j) {
      rank_one[0][j] = e("1");
      rank_one[j][j] = e("eps");
    }
    auto collapse = am::match(large, large.project(rank_one), half,
                              large_boundary, {0, 0}, large_options);
    check(collapse.success() && collapse.row_pole_loss == 1 &&
              collapse.determinant_valuation == 14,
          collapse.reason);
    near(at(collapse.value, 0, -1), B(-14));
    near(at(collapse.value, 1, -1), B(1));
    // Factoring a rational base exponent must precede the Q(x,log x) lift.
    auto fractional = AffineFrobeniusSeries::prepare(
        {{e("eps/x"), e("0")}, {e("0"), e("(1/2+eps)/x")}}, 0, 1, 1);
    auto fractional_gauge = fractional.project(
        AffineFrobeniusSeries::Matrix{{e("1"), e("1")}, {e("0"), e("eps")}});
    am::Options force_rows;
    force_rows.max_determinant_dimension = 1;
    auto fm = am::match(fractional, fractional_gauge, half,
                        am::Boundary{0, 1, {{B(1), B(0)}, {B(1), B(0)}}},
                        {0, 0}, force_rows);
    check(fm.success() && fm.used_row_normalization, fm.reason);
    B sqrt_two;
    acb_sqrt(sqrt_two.raw(), B(2).raw(), B::precision());
    near(at(fm.value, 1, -1), sqrt_two);
    // Independent physical solution: an analytic x-dependent gauge of
    // h(x)*{1, x*(x^eps-1)/eps}, h=1/(1-x). Unequal column-relative
    // cutoffs spuriously destroy the exact constant epsilon-pole relation.
    using namespace fuchsify::detail;
    AffineFrobeniusSeries::Matrix gauge{{e("1"), e("x")}, {e("1"), e("1")}};
    AffineFrobeniusSeries::Matrix base{{e("1/(1-x)"), e("0")},
                                       {e("1"), e("(1+eps)/x+1/(1-x)")}};
    auto physical_connection =
        subtract(multiply(multiply(gauge, base), inverse(gauge)),
                 multiply(gauge, derivative(inverse(gauge), 0)));
    B small = B::from_strings("0.125"), small_log;
    acb_log(small_log.raw(), small.raw(), B::precision());
    am::Boundary exact_boundary{0, 4, std::vector(2, std::vector<B>(5, B(0)))};
    B log_coefficient = small_log, h = B(1) / (B(1) - small);
    for (unsigned k = 0; k <= 4; ++k) {
      exact_boundary.coefficients[0][k] =
          h * (B(k == 0 ? 1 : 0) + small * small * log_coefficient);
      exact_boundary.coefficients[1][k] =
          h * (B(k == 0 ? 1 : 0) + small * log_coefficient);
      log_coefficient = log_coefficient * small_log / B(k + 2);
    }
    auto magnitude = [](const B &v) {
      arf_t upper;
      arf_init(upper);
      acb_get_abs_ubound_arf(upper, v.raw(), B::precision());
      double result = arf_get_d(upper, ARF_RND_CEIL);
      arf_clear(upper);
      return result;
    };
    double previous = 1;
    for (unsigned order : {6, 12}) {
      auto coherent =
          AffineFrobeniusSeries::prepare(physical_connection, 0, 1, order);
      auto recovered = am::match(coherent, small, exact_boundary, {-1, 0});
      check(recovered.success() && recovered.wronskian_verified &&
                recovered.determinant_valuation == 0,
            recovered.reason);
      double error = std::max({magnitude(at(recovered.value, 0, 0) - B(1)),
                               magnitude(at(recovered.value, 1, -1) - B(1)),
                               magnitude(at(recovered.value, 1, 0))});
      check(error < previous * 0.01, "coherent matching did not converge to "
                                     "the independent physical solution");
      previous = error;
    }
    check(previous < 1e-8, "independent physical matching accuracy");
    auto longer = AffineFrobeniusSeries::prepare(physical_connection, 0, 1, 7);
    auto broken = longer.terms();
    std::erase_if(broken.terms, [](const auto &t) {
      return t.column == 0 && t.power > Rational(6);
    });
    auto prevented = am::match(longer, broken, small, exact_boundary, {-1, 0});
    check(prevented.status == am::Status::NeedMoreXOrder,
          "Wronskian guard accepted an artificial tiny determinant from "
          "unequal x cutoffs");
    auto moving = AffineFrobeniusSeries::prepare({{e("1/(x+eps)")}}, 0, 1, 3);
    check(!moving.terms().wronskian_prefactor,
          "moving pole was incorrectly assigned an epsilon-unit Wronskian");
    am::Boundary moving_boundary{0, 2, {{B(1), B(0), B(0)}}};
    check(!am::match(moving, small, moving_boundary, {0, 0}).success(),
          "matching accepted an unsupported moving pole guard");
    if (argc > 1 && std::string(argv[1]) != "--cache") {
      Exact x = e("x"), d = e("2-2*eps");
      auto family = ibp::merge(
          ibp::quadratic_family(
              feynman::banana(2, {Rational(1), Rational(1), Rational(1)}), x),
          0, 1, x);
      ibp::PropagatorBasis basis(family);
      level::Options options;
      options.provider.executable = argv[1];
      auto prepared =
          level::prepare(basis, d, field, 0, {{3, 1, 0, 0, 0}}, options);
      check(prepared.success, prepared.reason);
      auto fs = AffineFrobeniusSeries::prepare(prepared.matrix, 0, 1, 3);
      auto values = fs.evaluate_laurent(fs.terms(), half, -4, 4);
      am::Boundary boundary{
          -4, 4, std::vector(prepared.matrix.size(), std::vector<B>(9, B(0)))};
      for (unsigned i = 0; i < boundary.coefficients.size(); ++i)
        for (unsigned j = 0; j < boundary.coefficients.size(); ++j)
          for (unsigned k = 0; k < 9; ++k)
            boundary.coefficients[i][k] =
                boundary.coefficients[i][k] +
                B(j + 1) * values.coefficients[i][j][k];
      auto fire = am::match(fs, half, boundary, {-2, 0});
      check(fire.success(), fire.reason);
      for (unsigned i = 0; i < boundary.coefficients.size(); ++i)
        near(at(fire.value, i, 0), B(i + 1));
      std::cout
          << "Actual FIRE sunrise affine frame matched; determinant valuation="
          << fire.determinant_valuation << "\n";
    }
    if (argc > 2 && std::string(argv[1]) == "--cache") {
      const unsigned cached_order =
          argc > 3 ? static_cast<unsigned>(std::stoul(argv[3])) : 3;
      artifacts::Store store(argv[2]);
      recursion::Options opts;
      auto provider = [&](const auto &basis, const auto &dimension,
                          const auto &field, auto parameter,
                          const auto &sources, const auto &budget) {
        auto reject = [](const std::vector<ibp::Integral> &,
                         const fire::Options &) -> fire::Result {
          throw std::runtime_error(
              "cache-only affine validation refuses to start FIRE");
        };
        auto reused = cached_level::prepare(store, basis, dimension, field,
                                            parameter, sources, budget, reject);
        check(reused.cache_hit && reused.result.success,
              "required Banana4 exact level missing from cache");
        return reused.result;
      };
      auto graph = recursion::prepare(feynman::example_family("banana4"), {},
                                      opts, provider);
      B::set_precision(384);
      auto point = B::from_strings("0.0625");
      for (const auto &node : graph.nodes)
        if (!node.scalar_leaf) {
          std::cout << "Cached Banana4 matching: preparing "
                    << node.closure.matrix.size() << " masters\n"
                    << std::flush;
          auto fuchs = fuchsify::prepare(node.closure.matrix, 0);
          check(fuchs.success, fuchs.reason);
          AffineFrobeniusSeries::Options endpoint_options;
          endpoint_options.max_dimension = 32;
          auto endpoint = AffineFrobeniusSeries::prepare(
              fuchs.matrix, 0, 1, cached_order, endpoint_options);
          std::cout << "  affine basis ready at x order " << cached_order
                    << "\n"
                    << std::flush;
          auto physical = endpoint.project(fuchs.transform);
          auto valuations = endpoint.valuation_metadata(physical);
          int low = static_cast<int>(std::min(
                  0L, *std::min_element(valuations.minimum_by_column.begin(),
                                        valuations.minimum_by_column.end()))),
              high = 2;
          am::Options matching_options;
          matching_options.max_dimension = 32;
          am::Result recovered;
          for (unsigned attempt = 0; attempt < 3; ++attempt) {
            auto evaluation =
                endpoint.evaluate_laurent(physical, point, low, high);
            am::Boundary source{
                low, high,
                std::vector(physical.rows,
                            std::vector<B>(high - low + 1, B(0)))};
            for (unsigned i = 0; i < physical.rows; ++i)
              for (unsigned j = 0; j < physical.columns; ++j)
                for (int k = low; k <= high; ++k)
                  source.coefficients[i][k - low] =
                      source.coefficients[i][k - low] +
                      B(j + 1) * evaluation.coefficients[i][j][k - low];
            recovered = am::match(endpoint, physical, point, source, {-2, 0},
                                  matching_options);
            if (recovered.status != am::Status::NeedMoreBoundary)
              break;
            high = recovered.required_boundary_high;
          }
          check(recovered.success() && recovered.wronskian_verified &&
                    recovered.determinant_valuation ==
                        recovered.expected_determinant_valuation,
                recovered.reason);
          for (unsigned i = 0; i < physical.columns; ++i)
            near(at(recovered.value, i, 0), B(i + 1));
          std::cout << "Cached Banana4 " << physical.rows
                    << " masters: matched; row loss=" << recovered.row_pole_loss
                    << ", determinant valuation="
                    << recovered.determinant_valuation << "\n"
                    << std::flush;
        }
    }
    std::cout << "Affine Laurent matching, explicit source demands and "
                 "physical gauges passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
