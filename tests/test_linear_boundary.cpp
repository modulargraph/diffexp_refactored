#include "diffexp/linear_boundary.hpp"
#include <iostream>
using namespace diffexp;
namespace lb = linear_boundary;
using B = Jet::Ball;
void check(bool b, const char *why) {
  if (!b)
    throw std::runtime_error(why);
}
void near(const B &a, const B &b) {
  check(acb_overlaps(a.raw(), b.raw()),
        "linear expression analytic result mismatch");
}
LaurentRows matrix(const std::vector<std::vector<B>> &values, int power,
                   int high) {
  LaurentRows out{
      std::min(power, high), high,
      std::vector(
          values.size(),
          std::vector(values[0].size(),
                      std::vector<B>(high - std::min(power, high) + 1, B(0))))};
  if (power <= high)
    for (unsigned i = 0; i < values.size(); ++i)
      for (unsigned j = 0; j < values[i].size(); ++j)
        out.coefficients[i][j][power - out.low] = values[i][j];
  return out;
}
int main() {
  try {
    B::set_precision(256);
    B huge;
    acb_one(huge.raw());
    acb_mul_2exp_si(huge.raw(), huge.raw(), 40);
    auto source = std::make_shared<LaurentBoundary>(
        LaurentBoundary{-1,
                        {{B(2), B(3), B(0), B(0), B(0), B(0)},
                         {B(5), B(7), B(0), B(0), B(0), B(0)}},
                        true});
    for (auto &row : source->values)
      for (auto &c : row)
        arb_add_error_2exp_si(acb_realref(c.raw()), -50);
    auto a = matrix({{huge, B(1)}, {huge, B(0)}}, -1, 5);
    auto inv = matrix({{B(0), B(1) / huge}, {B(1), B(-1)}}, 1, 5);
    auto identity = lb::identity(source, 5);
    auto first = lb::compose(a, identity, 4);
    auto second = lb::compose(inv, first, 3);
    check(second.leaf_source.get() == source.get(),
          "leaf source identity was not preserved");
    auto final = lb::materialize(second, 0);
    near(final.values[0][-final.low], B(3));
    near(final.values[1][-final.low], B(7));
    auto independent1 = apply_laurent_rows(a, *source, 2);
    auto independent2 = apply_laurent_rows(inv, independent1, 1);
    const auto &good = final.values[1][-final.low];
    const auto &bad = independent2.values[1][-independent2.low];
    check(mag_cmp(arb_radref(acb_realref(good.raw())),
                  arb_radref(acb_realref(bad.raw()))) < 0,
          "global composition did not preserve shared correlation");
    mag_t bound;
    mag_init(bound);
    mag_set_ui_2exp_si(bound, 1, -49);
    check(mag_cmp(arb_radref(acb_realref(good.raw())), bound) < 0,
          "composed identity amplified source uncertainty");
    mag_clear(bound);
    // Third ill-conditioned map acts on both reconstructed rows. The shared
    // source remains the original one even through an additional Laurent pole.
    auto third = lb::compose(a, second, 2);
    auto three = lb::materialize(third, 0);
    auto analytic = apply_laurent_rows(a, *source, 0);
    auto independent3 = apply_laurent_rows(a, independent2, 0);
    check(mag_cmp(arb_radref(acb_realref(three.values[0][-three.low].raw())),
                  arb_radref(acb_realref(
                      independent3.values[0][-independent3.low].raw()))) < 0,
          "three-map expression lost the shared source correlation");
    for (unsigned i = 0; i < 2; ++i)
      for (int k = analytic.low; k <= 0; ++k)
        near(three.values[i][k - three.low],
             analytic.values[i][k - analytic.low]);
    // Explicitly distinguish insufficient outer, inner, transform and leaf
    // windows; never silently pad unknown upper coefficients with zero.
    bool demanded = false;
    try {
      lb::compose(matrix({{B(1)}}, -2, 0), matrix({{B(1)}}, -3, 0), 0);
    } catch (const lb::CompositionDemand &d) {
      demanded = d.required_outer_high == 3 && d.required_inner_high == 2;
    }
    check(demanded, "composition window demand incorrect");
    auto short_identity = lb::identity(source, 0);
    demanded = false;
    try {
      lb::materialize(short_identity, 0);
    } catch (const lb::OperatorDemand &d) {
      demanded = d.required_high == 1;
    }
    check(demanded, "transform lookahead for leaf poles missing");
    auto short_source =
        std::make_shared<LaurentBoundary>(LaurentBoundary{0, {{B(1)}}, true});
    lb::Expression pole{matrix({{B(1)}}, -2, 2), short_source};
    demanded = false;
    try {
      lb::materialize(pole, 0);
    } catch (const BoundaryDemand &d) {
      demanded = d.required_high == 2;
    }
    check(demanded, "leaf upper source demand missing");
    lb::Options budget;
    budget.max_operations = 1;
    bool limited = false;
    try {
      lb::compose(a, identity, 4, budget);
    } catch (const std::length_error &) {
      limited = true;
    }
    check(limited, "composition ignored finite work limit");
    check(second.transform.low == 0, "Laurent lower-bound addition incorrect");

    // A wide master basis still represents one uncertain scalar source.
    // Adjacent differences cancel large exact coefficients before that source
    // is materialized. This covers the Henn-sized intermediate dimension that
    // previously hit the unrelated 64-column leaf budget.
    B seed(2);
    arb_add_error_2exp_si(acb_realref(seed.raw()), -50);
    auto one_leaf = std::make_shared<const LaurentBoundary>(
        LaurentBoundary{0, {{seed}}, true});
    std::vector<std::vector<B>> master_values(109, std::vector<B>(1));
    for (unsigned j = 0; j < 109; ++j)
      master_values[j][0] = huge + B(j);
    lb::Expression wide{matrix(master_values, 0, 0), one_leaf};
    std::vector<std::vector<B>> observable_values(257,
                                                 std::vector<B>(109));
    for (unsigned i = 0; i < 257; ++i) {
      observable_values[i][i % 108] = B(1);
      observable_values[i][i % 108 + 1] = B(-1);
    }
    auto observables = matrix(observable_values, 0, 0);
    auto composed = lb::compose(observables, wide, 0);
    check(composed.leaf_source == one_leaf,
          "wide master composition replaced the shared leaf");
    auto materialized = lb::materialize(composed, 0);
    const auto expected = -seed;
    mag_init(bound);
    mag_set_ui_2exp_si(bound, 1, -49);
    for (unsigned i = 0; i < materialized.values.size(); ++i) {
      check(acb_equal_si(composed.transform.coefficients[i][0][0].raw(), -1),
            "wide master difference failed exact cancellation");
      // Multiplication by -1 may round an Arb magnitude outward by an ulp;
      // require enclosure with less than twice the original tiny radius.
      const auto &value = materialized.values[i][0];
      check(acb_contains(value.raw(), expected.raw()) &&
                mag_cmp(arb_radref(acb_realref(value.raw())), bound) < 0,
            "wide master cancellation amplified the scalar leaf radius");
    }
    mag_clear(bound);
    auto constrained = lb::Options{};
    constrained.max_columns = 108;
    limited = false;
    try { lb::compose(observables, wide, 0, constrained); }
    catch (const std::invalid_argument &) { limited = true; }
    check(limited, "wide composition ignored its explicit column budget");
    constrained.max_columns = 1;
    limited = false;
    try { lb::identity(source, 0, constrained); }
    catch (const std::length_error &) { limited = true; }
    check(limited, "identity ignored the explicit matrix column budget");
    auto too_many_leaves = std::make_shared<const LaurentBoundary>(
        LaurentBoundary{0, Boundary(65, std::vector<B>{B(1)}), true});
    limited = false;
    try { lb::identity(too_many_leaves, 0); }
    catch (const std::length_error &) { limited = true; }
    check(limited, "wide intermediate support removed the scalar leaf budget");
    std::cout << "Global shared-leaf linear composition, Laurent demands and "
                 "correlation preservation passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
