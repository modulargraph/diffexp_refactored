#include "diffexp/polynomial_transport.hpp"
#include <chrono>
#include <iostream>
using namespace diffexp;
namespace pt = polynomial_transport;
using B = Jet::Ball;
void check(bool b, const std::string &why) {
  if (!b)
    throw std::runtime_error(why);
}
void near(const B &a, const B &b) {
  check(acb_overlaps(a.raw(), b.raw()),
        "polynomial analytic coefficient mismatch");
}
void close(const B &a, const B &b) {
  B delta = a - b;
  mag_t bound;
  mag_init(bound);
  mag_set_ui_2exp_si(bound, 1, -140);
  mag_t upper;
  mag_init(upper);
  acb_get_mag(upper, delta.raw());
  bool ok = mag_cmp(upper, bound) < 0;
  mag_clear(upper);
  mag_clear(bound);
  check(ok, "polynomial/rational same-order comparison failed");
}
int main(int argc, char **argv) {
  try {
    B::set_precision(256);
    ExactField field({"x", "eps", "I"});
    auto e = [&](const char *s) { return Exact(field, s); };
    // y' = 1/(1+x+eps), represented by a constant augmented source.
    // y(x)=log(1+x+eps)-log(1+eps), analytic oracle at each epsilon order.
    ExactEpsilonMatrix a{{e("0"), e("1/(1+x+eps)")}, {e("0"), e("0")}};
    auto compact = pt::compile(a, 0, 1, 80, 4);
    check(compact.fallback_rows == 0,
          "compact moving denominator did not use polynomial recurrence");
    Boundary boundary(2, std::vector<B>(5, B(0)));
    boundary[1][0] = B(1);
    B step = B::from_strings("0.125");
    auto value = pt::chart(compact, boundary, B(0), step, 80);
    B expected;
    acb_log(expected.raw(), (B(1) + step).raw(), B::precision());
    close(value[0][0], expected);
    for (unsigned k = 1; k <= 4; ++k) {
      B power;
      acb_pow_ui(power.raw(), (B(1) + step).raw(), k, B::precision());
      close(value[0][k], (B(1) / power - B(1)) * B(k % 2 ? 1 : -1) / B(k));
    }
    std::vector<RationalLineEntry> entries;
    auto expanded = feynman::scalar_functional_detail::epsilon_series(
        e("1/(1+x+eps)"), 1, 4);
    for (unsigned k = 0; k <= 4; ++k)
      entries.push_back({0, 1, k, expanded[k]});
    auto old = rational_chart(entries, boundary, B(0), step, 80);
    for (unsigned i = 0; i < 2; ++i)
      for (unsigned k = 0; k <= 4; ++k)
        close(value[i][k], old[i][k]);
    auto compiled = pt::compile(entries, 2, 80, 4);
    auto expanded_value = pt::chart(compiled, boundary, B(0), step, 80);
    for (unsigned k = 0; k <= 4; ++k)
      close(value[0][k], expanded_value[0][k]);
    // A genuine coupled matrix including epsilon denominator coupling and
    // complex polynomial coefficients; compare equal retained N and precision.
    ExactEpsilonMatrix coupled{{e("eps/(2-x+eps)"), e("I/(1+x)")},
                               {e("1/(3-x)"), e("(1+eps)/(2+x)")}};
    auto system = pt::compile(coupled, 0, 1, 64, 3);
    Boundary data{{B(2), B(1), B(0), B(0)}, {B(3), B(0), B(1), B(0)}};
    std::vector<RationalLineEntry> rational;
    for (unsigned i = 0; i < 2; ++i)
      for (unsigned j = 0; j < 2; ++j) {
        auto cs = feynman::scalar_functional_detail::epsilon_series(
            coupled[i][j], 1, 3);
        for (unsigned k = 0; k <= 3; ++k)
          if (!cs[k].is_zero())
            rational.push_back({i, j, k, cs[k]});
      }
    B center = B::from_strings("0.1"),
      h = B::from_strings("0.03125", "0.015625");
    auto new_matrix = pt::chart(system, data, center, h, 64),
         old_matrix = rational_chart(rational, data, center, h, 64);
    for (unsigned i = 0; i < 2; ++i)
      for (unsigned k = 0; k <= 3; ++k)
        close(new_matrix[i][k], old_matrix[i][k]);
    // Interned exact polynomials are shared between repeated adjoint blocks.
    std::vector<RationalLineEntry> repeated{{0, 0, 0, e("1/(1+x)")},
                                            {1, 1, 0, e("1/(1+x)")},
                                            {2, 2, 0, e("1/(1+x)")}};
    auto reuse = pt::compile(repeated, 3, 64, 0);
    check(reuse.polynomials.size() == 2,
          "equal block rows duplicated polynomial storage");
    // High-degree clearing safely selects the existing rational method.
    pt::Options fallback_options;
    fallback_options.max_degree = 4;
    auto fallback =
        pt::compile(std::vector<RationalLineEntry>{{0, 0, 0, e("1/(1+x^5)")}},
                    1, 16, 0, fallback_options);
    check(fallback.fallback_rows == 1,
          "clearing degree cap did not select safe fallback");
    auto fallback_value = pt::chart(fallback, {{B(1)}}, B(0), step, 16),
         fallback_old = rational_chart(
             std::vector<RationalLineEntry>{{0, 0, 0, e("1/(1+x^5)")}},
             {{B(1)}}, B(0), step, 16);
    near(fallback_value[0][0], fallback_old[0][0]);
    ExactEpsilonMatrix mixed{{e("1/(1+x^5)"), e("1")}, {e("1/(1+x)"), e("0")}};
    auto mixed_system = pt::compile(mixed, 0, 1, 32, 0, fallback_options);
    check(mixed_system.polynomial_rows == 1 && mixed_system.fallback_rows == 1,
          "mixed per-row fallback selection failed");
    std::vector<RationalLineEntry> mixed_entries{
        {0, 0, 0, mixed[0][0]}, {0, 1, 0, mixed[0][1]}, {1, 0, 0, mixed[1][0]}};
    auto mixed_new = pt::chart(mixed_system, {{B(1)}, {B(2)}}, B(0), step, 32);
    auto mixed_old =
        rational_chart(mixed_entries, {{B(1)}, {B(2)}}, B(0), step, 32);
    for (unsigned i = 0; i < 2; ++i)
      close(mixed_new[i][0], mixed_old[i][0]);
    auto costly = pt::compile(
        std::vector<RationalLineEntry>{{0, 0, 0, e("1/(1+x^8)")}}, 1, 8, 0);
    check(costly.fallback_rows == 1,
          "clearing cost estimate did not retain rational fallback");
    bool rejected = false;
    try {
      pt::chart(compact, boundary, B(-1), step, 8);
    } catch (const std::domain_error &) {
      rejected = true;
    }
    check(rejected, "zero leading polynomial denominator accepted");
    auto limited = system;
    limited.options.max_operations = 1;
    rejected = false;
    try {
      pt::chart(limited, data, center, h, 64);
    } catch (const std::length_error &) {
      rejected = true;
    }
    check(rejected, "finite-lag operation budget ignored");
    if (argc > 1) {
      auto old_compiled = compile_rational_entries(rational);
      auto start = std::chrono::steady_clock::now();
      for (unsigned trial = 0; trial < 20; ++trial)
        (void)rational_chart(old_compiled, data, center, h, 64);
      auto middle = std::chrono::steady_clock::now();
      for (unsigned trial = 0; trial < 20; ++trial)
        (void)pt::chart(system, data, center, h, 64);
      auto finish = std::chrono::steady_clock::now();
      std::cout << "Same N64/bits256, 20 coupled charts: rational="
                << std::chrono::duration<double>(middle - start).count()
                << " s; polynomial="
                << std::chrono::duration<double>(finish - middle).count()
                << " s\n";
    }
    std::cout << "Finite-lag polynomial recurrence, compact epsilon "
                 "denominator, shared blocks and fallback passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
