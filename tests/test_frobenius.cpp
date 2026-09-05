#include "diffexp/frobenius.hpp"
#include "diffexp/rational_transport.hpp"
#include <iostream>
using namespace diffexp;
using B = kernel::ComplexBall;
void check(bool b, const char *m) {
  if (!b)
    throw std::runtime_error(m);
}
void near(const B &a, const B &b, double tolerance = 1e-40) {
  B e = a - b;
  arf_t bound;
  arf_init(bound);
  acb_get_abs_ubound_arf(bound, e.raw(), B::precision());
  double error = arf_get_d(bound, ARF_RND_CEIL);
  arf_clear(bound);
  check(error < tolerance, "numerical reference disagreement");
}
template <class F> void rejects(F f) {
  bool rejected = false;
  try {
    f();
  } catch (const std::domain_error &) {
    rejected = true;
  }
  check(rejected, "unsupported singular input accepted");
}
int main() {
  B::set_precision(256);
  ExactField field({"x", "eps"});
  auto e = [&](const char *s) { return Exact(field, s); };
  // eps/(x*(1-eps)) gives exp((eps+eps^2+...) log x).
  auto scalar = FrobeniusSeries::prepare({{e("eps/(x*(1-eps))")}}, 0, 1, 3, 3);
  auto x = B::from_strings("0.5");
  B l;
  acb_log(l.raw(), x.raw(), 256);
  auto f = scalar.evaluate(x);
  near(f[0][0][0], B(1));
  near(f[0][0][1], l);
  near(f[0][0][2], l + l * l / B(2));
  near(f[0][0][3], l + l * l + l * l * l / B(6));
  check(!FrobeniusSeries::omitted_tail_certified,
        "must not certify a missing tail");
  // Integer-separated exponents: resonant forcing creates x log x.
  auto resonant = FrobeniusSeries::prepare(
      {{e("0"), e("0")}, {e("1"), e("1/x")}}, 0, 1, 4, 0);
  bool has_log = false;
  for (auto &m : resonant.monomials())
    if (m.row == 1 && m.power == Rational(1) && m.log_degree == 1 &&
        m.coefficient == Rational(1))
      has_log = true;
  check(has_log, "integer resonance did not produce exact logarithm");
  FrobeniusSeries::Boundary boundary{{B(2)}, {B(3)}};
  auto c = resonant.match(B(1), boundary);
  auto y = resonant.solution(x, c);
  near(y[0][0], B(2));
  near(y[1][0], x * (B(3) + B(2) * l));
  auto primitive = resonant.primitive(x, c);
  near(primitive[0][0], B(2) * x);
  near(primitive[1][0], x * x * (B(1) + l));
  // A regular chart gives an independent ordinary-system continuation.
  std::vector<RationalLineEntry> entries{{1, 0, 0, e("1")},
                                         {1, 1, 0, e("1/x")}};
  auto ordinary =
      rational_chart(entries, boundary, B(1), B::from_strings("-0.25"), 90);
  auto singular = resonant.solution(B::from_strings("0.75"), c);
  for (unsigned i = 0; i < 2; ++i)
    near(ordinary[i][0], singular[i][0]);
  // Nontrivial constant eigenvector frame with colliding affine eps roots.
  auto collision = FrobeniusSeries::prepare(
      {{e("(-1+eps)/x"), e("0")}, {e("eps/x"), e("(-1+2*eps)/x")}}, 0, 1, 0, 3);
  auto cf = collision.evaluate(x);
  near(cf[0][0][0], B(1) / x);
  near(cf[1][0][1], l / x);
  near(cf[1][0][2], B(3) * l * l / (B(2) * x));
  auto analytic =
      FrobeniusSeries::prepare({{e("eps/x+1/(1-x)")}}, 0, 1, 100, 3);
  auto af = analytic.evaluate(B::from_strings("0.25"));
  B ll;
  acb_log(ll.raw(), B::from_strings("0.25").raw(), 256);
  near(af[0][0][0], B(4) / B(3));
  near(af[0][0][2], B(2) * ll * ll / B(3));
  // Fractional rational exponents and a resonant primitive logarithm.
  auto half = FrobeniusSeries::prepare({{e("1/(2*x)")}}, 0, 1, 0, 0);
  near(half.evaluate(B(4))[0][0][0], B(2));
  auto pole = FrobeniusSeries::prepare({{e("-1/x")}}, 0, 1, 0, 0);
  near(pole.primitive(x, {{B(1)}})[0][0], l);
  // Exact rational-row projection removes a divergent eigenmode before Acb
  // matching; a numerical near-zero coefficient is never chopped.
  auto projected = FrobeniusSeries::prepare(
      {{e("-1/x"), e("0")}, {e("-1/x"), e("0")}}, 0, 1, 4, 1);
  auto pc = projected.match(B(1), {{B(2), B(0)}, {B(5), B(0)}});
  auto terms = projected.project({e("1"), e("-1")}, 0, 1);
  for (const auto &m : terms)
    check(m.power >= Rational(0), "projection did not cancel divergent power");
  near(projected.integrate_projected(terms, x, pc)[0], -B(3) * x);
  auto weighted = projected.project({e("x"), e("0")}, 0, 1);
  near(projected.integrate_projected(weighted, x, pc)[0], B(2) * x);
  rejects([&] {
    projected.integrate_projected(projected.project({e("1"), e("0")}, 0, 1), x,
                                  pc);
  });
  // M^2=eps I: exp(M log x) gives an independent closed-form oracle
  // even though the nonzero-epsilon eigenvalues involve sqrt(eps).
  auto jordan2 = FrobeniusSeries::prepare(
      {{e("0"), e("1/x")}, {e("eps/x"), e("0")}}, 0, 1, 2, 3);
  FrobeniusSeries::Boundary seed2(2, std::vector<B>(4, B(0)));
  seed2[1][0] = B(1);
  auto jc2 = jordan2.match(B(1), seed2);
  auto jy2 = jordan2.solution(x, jc2);
  auto divided_power = [&](unsigned n) {
    B v(1);
    for (unsigned j = 1; j <= n; ++j)
      v = v * l / B(j);
    return v;
  };
  for (unsigned k = 0; k < 4; ++k) {
    near(jy2[0][k], divided_power(2 * k + 1));
    near(jy2[1][k], divided_power(2 * k));
  }
  // Length-three Jordan ladder, M^3=eps I, with exact log^(3k+2)
  // coefficients. Projection and matching retain these enhanced log degrees.
  auto jordan3 = FrobeniusSeries::prepare({{e("0"), e("1/x"), e("0")},
                                           {e("0"), e("0"), e("1/x")},
                                           {e("eps/x"), e("0"), e("0")}},
                                          0, 1, 2, 3);
  FrobeniusSeries::Boundary seed3(3, std::vector<B>(4, B(0)));
  seed3[2][0] = B(1);
  auto jc3 = jordan3.match(B(1), seed3);
  auto jy3 = jordan3.solution(x, jc3);
  for (unsigned k = 0; k < 4; ++k)
    for (unsigned i = 0; i < 3; ++i)
      near(jy3[i][k], divided_power(3 * k + 2 - i));
  bool exact_ladder = false;
  for (const auto &m : jordan3.monomials())
    if (m.row == 0 && m.column == 2 && m.epsilon == 0 && m.log_degree == 2)
      exact_ladder = m.coefficient == Rational("1/2");
  check(exact_ladder, "missing exact Jordan factorial");
  auto jp = jordan3.integrate_projected(
      jordan3.project({e("1"), e("0"), e("0")}, 0, 1), x, jc3);
  for (unsigned k = 0; k < 4; ++k) {
    unsigned degree = 3 * k + 2;
    B oracle(0);
    for (unsigned j = 0; j <= degree; ++j)
      oracle = oracle + B((degree - j) % 2 ? -1 : 1) * divided_power(j);
    near(jp[k], x * oracle);
  }
  // Nonresonant Taylor levels solve delta I-N, including nilpotent mixing.
  auto ordinary_jordan = FrobeniusSeries::prepare(
      {{e("1"), e("1/x")}, {e("0"), e("1")}}, 0, 1, 80, 0);
  auto oc = ordinary_jordan.match(B(1), {{B(0)}, {B(1)}});
  auto oy = ordinary_jordan.solution(x, oc);
  B expx;
  auto xm1 = x - B(1);
  acb_exp(expx.raw(), xm1.raw(), B::precision());
  near(oy[0][0], expx * l);
  near(oy[1][0], expx);
  // Mixed Jordan blocks after a nontrivial rational similarity: frame
  // completeness and R*T=T*J are checked inside preparation.
  using namespace frobenius_detail;
  Matrix jt = zeros(4, 4),
         change{{Rational(1), Rational(2), Rational(0), Rational(1)},
                {Rational(0), Rational(1), Rational(1), Rational(0)},
                {Rational(0), Rational(0), Rational(1), Rational(1)},
                {Rational(0), Rational(0), Rational(0), Rational(1)}};
  jt[0][1] = Q(1);
  jt[2][2] = Q(2);
  jt[3][3] = Q(2);
  auto residue = multiply(multiply(change, jt), inverse(change));
  FrobeniusSeries::ExactMatrix transformed(4, std::vector<Exact>(4, e("0")));
  for (unsigned i = 0; i < 4; ++i)
    for (unsigned j = 0; j < 4; ++j)
      transformed[i][j] = e("1/x") * e("0").constant(residue[i][j]);
  auto mixed = FrobeniusSeries::prepare(transformed, 0, 1, 0, 0);
  FrobeniusSeries::Boundary mb(4, std::vector<B>(1, B(0)));
  for (unsigned i = 0; i < 4; ++i)
    mb[i][0] = ball(change[i][1] + change[i][2] + change[i][3]);
  auto my = mixed.solution(x, mixed.match(B(1), mb));
  for (unsigned i = 0; i < 4; ++i)
    near(my[i][0], ball(change[i][0]) * l + ball(change[i][1]) +
                       x * x * ball(change[i][2] + change[i][3]));
  Matrix j32 = zeros(5, 5), similarity = zeros(5, 5);
  for (unsigned i = 0; i < 5; ++i) {
    similarity[i][i] = Q(1);
    if (i + 1 < 5)
      similarity[i][i + 1] = Q(i + 1);
  }
  j32[0][1] = Q(1);
  j32[1][2] = Q(1);
  j32[3][4] = Q(1);
  auto jf = frame(multiply(multiply(similarity, j32), inverse(similarity)));
  std::vector<unsigned> chain_lengths;
  unsigned length = 0;
  for (int next : jf.successor) {
    ++length;
    if (next < 0) {
      chain_lengths.push_back(length);
      length = 0;
    }
  }
  check(chain_lengths == std::vector<unsigned>({3, 2}),
        "mixed-length Jordan filtration");
  rejects([&] { FrobeniusSeries::prepare({{e("1/x^2")}}, 0, 1, 2, 0); });
  rejects([&] {
    FrobeniusSeries::prepare({{e("0"), e("2/x")}, {e("1/x"), e("0")}}, 0, 1, 2,
                             0);
  });
  std::cout << "Frobenius exact recurrence, epsilon collision, matching and "
               "primitives passed\n";
}
