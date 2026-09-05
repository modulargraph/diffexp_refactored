#include "diffexp/affine_operator.hpp"
#include "diffexp/epsilon_gauge.hpp"
#include "diffexp/level_cache.hpp"
#include "diffexp/recursion_graph.hpp"
#include <chrono>
#include <iostream>
using namespace diffexp;
namespace ao = affine_operator;
namespace am = affine_matching;
using B = kernel::ComplexBall;
void check(bool yes, const std::string &why) {
  if (!yes)
    throw std::runtime_error(why);
}
void near(const B &a, const B &b) {
  check(acb_overlaps(a.raw(), b.raw()),
        "operator analytic coefficient mismatch");
}
int main(int argc, char **argv) {
  try {
    B::set_precision(256);
    ExactField field({"x", "eps"});
    auto e = [&](const char *s) { return Exact(field, s); };
    B x = B::from_strings("0.125");
    ExactField symbolic({"x", "ell"});
    Exact sx(symbolic, "x"), sl(symbolic, "ell");
    auto dense = (sx + sl + sx.constant(1)).pow(30);
    am::detail::ExactPointEvaluator evaluator(x);
    auto dense_value = evaluator(dense);
    B logx, oracle;
    acb_log(logx.raw(), x.raw(), B::precision());
    auto base = B(1) + x + logx;
    acb_pow_ui(oracle.raw(), base.raw(), 30, B::precision());
    near(dense_value, oracle);
    bool denominator_rejected = false;
    try {
      (void)am::detail::ExactPointEvaluator(B(1))(sx.constant(1) /
                                                  (sx - sx.constant(1)));
    } catch (const am::detail::PrecisionNeeded &) {
      denominator_rejected = true;
    }
    check(denominator_rejected,
          "direct exact evaluator accepted zero denominator");

    auto scalar = AffineFrobeniusSeries::prepare({{e("eps/x")}}, 0, 1, 4);
    ao::Options lazy_options;
    lazy_options.max_row_normalization_steps = 32;
    lazy_options.max_exact_operations = 40;
    auto lazy = ao::prepare(scalar, x, 1, lazy_options);
    check(lazy.success(),
          "normalization eagerly spent the unused maximum probe depth: " +
              lazy.reason);
    auto prepared = ao::prepare(scalar, x, 8);
    check(prepared.success(), prepared.reason);
    auto primitive = scalar.dr_integral_from_zero(scalar.terms());
    auto op = ao::compose(prepared, scalar, primitive, x, 5);
    check(op.success(), op.reason);
    for (int k = 0; k <= 5; ++k)
      near(op.matrix.coefficients[0][0][k - op.matrix.low], k % 2 ? -x : x);
    am::Boundary boundary{
        -2, 5, {{B(3), B(2), B(1), B(0), B(0), B(0), B(0), B(0)}}};
    auto result = ao::apply_operator(op, boundary, {-2, 0});
    check(result.success(), result.reason);
    near(result.value.coefficients[0][-2 - result.value.low], B(3) * x);
    near(result.value.coefficients[0][-result.value.low], B(2) * x);
    auto short_op = ao::compose(prepared, scalar, primitive, x, 1);
    check(ao::apply_operator(short_op, boundary, {-2, 0}).status ==
              am::Status::NeedMoreSource,
          "missing operator high padded to zero");
    auto colliding = AffineFrobeniusSeries::prepare(
        {{e("eps/x"), e("1/x")}, {e("0"), e("0")}}, 0, 1, 4);
    auto collision_inverse = ao::prepare(colliding, x, 8);
    check(collision_inverse.success() && collision_inverse.wronskian_verified,
          collision_inverse.reason);
    ao::Options numeric_options;
    numeric_options.numeric_operator_convolution = true;
    auto numeric_collision = ao::prepare(colliding, x, 8, numeric_options);
    check(numeric_collision.success() && numeric_collision.used_numeric_convolution,
          numeric_collision.reason);
    for (unsigned i = 0; i < 2; ++i)
      for (unsigned j = 0; j < 2; ++j)
        for (int k = collision_inverse.inverse.low; k <= 8; ++k)
          near(numeric_collision.inverse.coefficients[i][j][k - numeric_collision.inverse.low],
               collision_inverse.inverse.coefficients[i][j][k - collision_inverse.inverse.low]);
    auto numeric_identity = ao::compose(numeric_collision, colliding, colliding.terms(), x, 3);
    check(numeric_identity.success(), numeric_identity.reason);
    for (unsigned i = 0; i < 2; ++i)
      for (unsigned j = 0; j < 2; ++j)
        for (int k = numeric_identity.matrix.low; k <= 3; ++k)
          near(numeric_identity.matrix.coefficients[i][j][k - numeric_identity.matrix.low],
               B(k == 0 && i == j ? 1 : 0));
    auto identity =
        ao::compose(collision_inverse, colliding, colliding.terms(), x, 3);
    check(identity.success(), identity.reason);
    for (unsigned i = 0; i < 2; ++i)
      for (unsigned j = 0; j < 2; ++j)
        for (int k = identity.matrix.low; k <= 3; ++k)
          near(identity.matrix.coefficients[i][j][k - identity.matrix.low],
               B(k == 0 && i == j ? 1 : 0));
    // A near-dependent polynomial frame makes the normalizer subtract
    // two huge independently rounded coefficients; the exact result is zero.
    auto constant_series = AffineFrobeniusSeries::prepare(
        {{e("0"), e("0")}, {e("0"), e("0")}}, 0, 1, 2);
    auto big = e("2^512+1");
    auto adversarial = constant_series.project(std::vector<std::vector<Exact>>{
        {e("1+eps+eps^2"), e("1")},
        {big * e("1+eps+eps^2"), big + e("eps")}});
    auto guarded = ao::prepare(constant_series, adversarial, x, 3, numeric_options);
    check(guarded.success(), guarded.reason);
    check(guarded.numeric_convolution_fallbacks > 0,
          "cancelling operator factors did not activate exact fallback");
    for (int k = guarded.inverse.low; k <= 3; ++k) {
      auto coefficient = [&](unsigned i, unsigned j) -> const B& {
        return guarded.inverse.coefficients[i][j][k - guarded.inverse.low];
      };
      B huge = am::detail::ExactPointEvaluator(x)(big);
      near(coefficient(1, 0), k == -1 ? -huge : B(0));
      near(coefficient(1, 1), B(k == -1 ? 1 : 0));
    }
    auto pole = colliding.project(std::vector<Exact>{e("1/eps^2"), e("0")});
    auto insufficient = ao::compose(collision_inverse, colliding, pole, x, 8);
    check(insufficient.status == am::Status::NeedMoreSource &&
              insufficient.required_inverse_high > 8,
          "inverse composition demand missing");
    auto projected = ao::compose(collision_inverse, colliding, pole, x, 3);
    check(projected.success(), projected.reason);
    am::Boundary uncertain{0, 0, {{B(1)}, {B(2)}}};
    check(ao::apply_operator(projected, uncertain, {0, 0}).status ==
              am::Status::NeedMoreBoundary,
          "operator poles did not demand boundary upper");
    // Fifteen columns with rank-one leading frame; one normalization serves
    // every inverse column and two unrelated physical observables.
    const unsigned d = 15;
    auto a = fuchsify::detail::zeros(d, d, e("0"));
    auto many = AffineFrobeniusSeries::prepare(a, 0, 1, 4, [] {
      AffineFrobeniusSeries::Options o;
      o.max_dimension = 32;
      return o;
    }());
    auto gauge = fuchsify::detail::zeros(d, d, e("0"));
    for (unsigned i = 0; i < d; ++i)
      for (unsigned j = 0; j < d; ++j)
        gauge[i][j] = i == j ? e("1+eps") : e("1");
    ao::Options opts;
    opts.max_dimension = 32;
    auto physical = many.project(gauge);
    auto inverse = ao::prepare(many, physical, x, 4, opts);
    check(inverse.success() && inverse.normalization_steps == 1,
          inverse.reason);
    auto plain = ao::compose(inverse, many, many.terms(), x, 3, opts);
    check(plain.success(), plain.reason);
    // (eps I + 11^T)^-1 = I/eps - 11^T/[eps(eps+d)].
    for (unsigned i = 0; i < d; ++i)
      for (unsigned j = 0; j < d; ++j) {
        near(plain.matrix.coefficients[i][j][-1 - plain.matrix.low],
             B(i == j ? 1 : 0) - B(1) / B(d));
        near(plain.matrix.coefficients[i][j][-plain.matrix.low],
             B(1) / B(d * d));
      }
    auto same = ao::compose(inverse, many, physical, x, 2, opts);
    check(same.success(), same.reason);
    for (unsigned i = 0; i < d; ++i)
      for (unsigned j = 0; j < d; ++j)
        near(same.matrix.coefficients[i][j][-same.matrix.low],
             B(i == j ? 1 : 0));
    // The prepared identity operator contracts shared boundary uncertainty only
    // after its exact algebraic cancellation is enclosed in matrix
    // coefficients.
    am::Boundary balls{0, 3, std::vector(d, std::vector<B>(4, B(0)))};
    for (auto &row : balls.coefficients) {
      row[0] = B(1);
      arb_add_error_2exp_si(acb_realref(row[0].raw()), -40);
    }
    auto contracted = ao::apply_operator(same, balls, {0, 0}, opts);
    check(contracted.success(), contracted.reason);
    for (unsigned i = 0; i < d; ++i) {
      auto &value = contracted.value.coefficients[i][-contracted.value.low];
      near(value, B(1));
      mag_t bound;
      mag_init(bound);
      mag_set_ui_2exp_si(bound, 1, -39);
      check(mag_cmp(arb_radref(acb_realref(value.raw())), bound) < 0,
            "shared boundary radius amplified before contraction");
      mag_clear(bound);
    }
    if (argc > 1) {
      artifacts::Store store(argv[1]);
      recursion::Options graph_options;
      auto provider = [&](const auto &basis, const auto &dim, const auto &f,
                          auto xi, const auto &sources, const auto &budget) {
        auto reject = [](const std::vector<ibp::Integral> &,
                         const fire::Options &) -> fire::Result {
          throw std::runtime_error("operator test cache miss; FIRE prohibited");
        };
        auto cached = cached_level::prepare(store, basis, dim, f, xi, sources,
                                            budget, reject);
        check(cached.cache_hit && cached.result.success,
              "operator cache unavailable");
        return cached.result;
      };
      auto graph = recursion::prepare(feynman::example_family("banana"), {},
                                      graph_options, provider);
      const auto &node = graph.nodes[0];
      auto epsg = epsilon_diagonal_gauge(node.closure.matrix, 1);
      auto fuchs = fuchsify::prepare(epsg.matrix, 0);
      check(fuchs.success, fuchs.reason);
      AffineFrobeniusSeries::Options endpoint_options;
      endpoint_options.max_dimension = 32;
      B::set_precision(384);
      const unsigned retained_order = argc > 2 ? std::stoul(argv[2]) : 8;
      auto endpoint = AffineFrobeniusSeries::prepare(
          fuchs.matrix, 0, 1, retained_order, endpoint_options);
      auto frame = endpoint.project(fuchs.transform);
      const B matchpoint = B::from_strings("1/128");
      auto started = std::chrono::steady_clock::now();
      auto cached_inverse = ao::prepare(endpoint, frame, matchpoint, 6, opts);
      check(cached_inverse.success(), cached_inverse.reason);
      auto midpoint = std::chrono::steady_clock::now();
      auto numerical_options = opts;
      numerical_options.numeric_operator_convolution = true;
      auto numeric_inverse =
          ao::prepare(endpoint, frame, matchpoint, 6, numerical_options);
      check(numeric_inverse.success(), numeric_inverse.reason);
      check(numeric_inverse.numeric_convolution_fallbacks == 0,
            "well-conditioned cached seven-master path unexpectedly fell back");
      auto finished = std::chrono::steady_clock::now();
      double old_radius = 0, new_radius = 0;
      for (unsigned i = 0; i < frame.rows; ++i)
        for (unsigned j = 0; j < frame.rows; ++j)
          for (int k = cached_inverse.inverse.low; k <= 6; ++k) {
            const auto &old =
                cached_inverse.inverse
                    .coefficients[i][j][k - cached_inverse.inverse.low];
            const auto &now =
                numeric_inverse.inverse
                    .coefficients[i][j][k - numeric_inverse.inverse.low];
            near(old, now);
            old_radius = std::max(
                old_radius, mag_get_d(arb_radref(acb_realref(old.raw()))));
            new_radius = std::max(
                new_radius, mag_get_d(arb_radref(acb_realref(now.raw()))));
          }
      std::cout << "Cached N" << retained_order << " bits384 inverse: exact="
                << std::chrono::duration<double>(midpoint - started).count()
                << "s, numeric="
                << std::chrono::duration<double>(finished - midpoint).count()
                << "s; fallbacks=" << numeric_inverse.numeric_convolution_fallbacks
                << "; maximum radii " << old_radius << " / " << new_radius
                << "\n";
      auto cached_identity =
          ao::compose(cached_inverse, endpoint, frame, matchpoint, 2, opts);
      check(cached_identity.success(), cached_identity.reason);
      for (unsigned i = 0; i < frame.rows; ++i)
        for (unsigned j = 0; j < frame.rows; ++j)
          for (int k = cached_identity.matrix.low; k <= 2; ++k)
            near(cached_identity.matrix
                     .coefficients[i][j][k - cached_identity.matrix.low],
                 B(i == j && k == 0 ? 1 : 0));
      std::cout
          << "Cached seven-master inverse prepared once; determinant valuation="
          << cached_inverse.determinant_valuation << "\n";
    }
    std::cout << "Prepared affine operator reuse, analytic inverse, pole "
                 "demand and matrix composition passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
