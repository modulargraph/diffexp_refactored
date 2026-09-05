#include "diffexp/families.hpp"
#include "diffexp/frobenius.hpp"
#include "diffexp/fuchsify.hpp"
#include "diffexp/jet.hpp"
#include "diffexp/level_preparation.hpp"
#include <iostream>
using namespace diffexp;
using B = kernel::ComplexBall;
void check(bool v, const char *why) {
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
  check(error < 1e-35, "transformed solution analytic oracle mismatch");
}
FrobeniusSeries::Boundary numerical(const fuchsify::Matrix &t, const B &x,
                                    const FrobeniusSeries::Boundary &v) {
  unsigned k = v[0].size();
  Jet eps(0, k, B::precision());
  if (k > 1)
    eps.set(1, B(1));
  auto xx = eps.constant(0);
  xx.set(0, x);
  FrobeniusSeries::Boundary out(t.size(), std::vector<B>(k, B(0)));
  for (unsigned i = 0; i < t.size(); ++i)
    for (unsigned j = 0; j < t[i].size(); ++j) {
      auto a = evaluate(data::Reader(t[i][j].str()).read(), eps,
                        {{"x", xx}, {"eps", eps}});
      for (unsigned e = 0; e < k; ++e)
        for (unsigned l = 0; l <= e; ++l)
          out[i][e] = out[i][e] + a.at(l) * v[j][e - l];
    }
  return out;
}
int main(int argc, char **argv) {
  try {
    B::set_precision(256);
    ExactField field({"x", "eps"});
    Exact x(field, "x");
    auto e = [&](const char *s) { return Exact(field, s); };
    using namespace fuchsify::detail;
    // Manufacture a rational gauge of an analytic diagonal system. The physical
    // fundamental columns are T*(1,x), giving a fully independent solution
    // oracle.
    fuchsify::Matrix base{{e("0"), e("0")}, {e("0"), e("1/x")}};
    fuchsify::Matrix original_gauge{{e("1"), e("1/x^2")},
                                    {e("1"), e("1+1/x^2")}};
    auto original = subtract(
        multiply(multiply(original_gauge, base), inverse(original_gauge)),
        multiply(original_gauge, derivative(inverse(original_gauge), 0)));
    auto prepared = fuchsify::prepare(original, 0);
    check(prepared.success, prepared.reason.c_str());
    check(prepared.identity_verified && prepared.final_pole_order <= 1 &&
              prepared.initial_pole_order > 1,
          "rational gauge reduction certificate");
    auto series = FrobeniusSeries::prepare(prepared.matrix, 0, 1, 40, 2);
    B at = B::from_strings("0.5");
    FrobeniusSeries::Boundary physical{{B(2), B(0), B(0)}, {B(3), B(0), B(0)}};
    auto constants = series.match(
        B(1), numerical(prepared.inverse_transform, B(1), physical));
    auto solution =
        numerical(prepared.transform, at, series.solution(at, constants));
    near(solution[0][0], B(1) + B(1) / at);
    near(solution[1][0], B(1) + B(1) / at + at);
    auto row = prepared.transform_rows({{e("1"), e("-1")}});
    auto observable = numerical(row, at, series.solution(at, constants));
    near(observable[0][0], -at);
    // Diagonal valuation shears with epsilon-dependent residue.
    auto direct = fuchsify::prepare(
        {{e("eps/x"), e("1/x^3")}, {e("0"), e("2*eps/x")}}, 0);
    check(direct.success && direct.steps == 1 && direct.identity_verified,
          "direct valuation shear");
    (void)FrobeniusSeries::prepare(direct.matrix, 0, 1, 3, 3);
    auto irregular = fuchsify::prepare({{e("1/x^2")}}, 0);
    check(!irregular.success &&
              irregular.reason.find("irregular") != std::string::npos,
          "scalar irregular endpoint rejection");
    auto irregular_matrix =
        fuchsify::prepare({{e("0"), e("1/x^2")}, {e("1/x^2"), e("0")}}, 0);
    check(!irregular_matrix.success &&
              irregular_matrix.reason.find("irregular") != std::string::npos,
          "nonnilpotent leading matrix rejection");
    auto cap = fuchsify::Options{};
    cap.max_power = 1;
    auto limited = fuchsify::prepare(original, 0, cap);
    check(!limited.success &&
              limited.reason.find("budget") != std::string::npos,
          "finite endpoint power budget");
    if (argc > 1) {
      level::Options opts;
      opts.provider.executable = argv[1];
      Exact d(field, "2-2*eps");
      for (unsigned loops : {2, 4}) {
        auto family = ibp::quadratic_family(
            feynman::banana(loops,
                            std::vector<Rational>(loops + 1, Rational(1))),
            x);
        if (loops == 4) {
          family = ibp::merge(family, 0, 1, e("1/3"));
          family = ibp::merge(family, 0, 1, e("2/5"));
        }
        family = ibp::merge(family, 0, 1, x);
        ibp::PropagatorBasis basis(family);
        ibp::Integral target(basis.denominators.size(), 0);
        target[0] = loops == 2 ? 3 : 4;
        target[1] = 1;
        auto level = level::prepare(basis, d, field, 0, {target}, opts);
        check(level.success, level.reason.c_str());
        for (unsigned end = 0; end < 2; ++end) {
          auto matrix = level.matrix;
          if (end)
            for (auto &row : matrix)
              for (auto &v : row)
                v = -v.substitute(std::vector<Exact>{e("1-x"), e("eps")});
          auto ff = fuchsify::prepare(matrix, 0);
          check(ff.success, ff.reason.c_str());
          check(ff.identity_verified, "actual FIRE endpoint certificate");
          (void)FrobeniusSeries::prepare(ff.matrix, 0, 1, 4, 2);
          std::cout << "Banana" << loops << " endpoint=" << end
                    << " masters=" << matrix.size()
                    << " poles=" << ff.initial_pole_order << "->"
                    << ff.final_pole_order << "\n";
        }
      }
    }
    std::cout << "Exact endpoint gauge, transformed solution and irregular "
                 "rejection passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
