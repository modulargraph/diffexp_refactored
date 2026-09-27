#include "diffexp/affine_frobenius.hpp"
#include "diffexp/families.hpp"
#include "diffexp/level_preparation.hpp"
#include <iostream>
using namespace diffexp;
using B = kernel::ComplexBall;
void check(bool v, const char *why) {
  if (!v)
    throw std::runtime_error(why);
}
void near(const B &a, const B &b) {
  auto error = a - b;
  arf_t upper;
  arf_init(upper);
  acb_get_abs_ubound_arf(upper, error.raw(), B::precision());
  double e = arf_get_d(upper, ARF_RND_CEIL);
  arf_clear(upper);
  check(e < 1e-40, "affine analytic oracle mismatch");
}
template <class F> void rejects(F f, const char *why) {
  bool bad = false;
  try {
    f();
  } catch (const std::domain_error &) {
    bad = true;
  }
  check(bad, why);
}
int main(int argc, char **argv) {
  try {
    B::set_precision(256);
    ExactField field({"x", "eps"});
    for(unsigned mode:{0,1,2,3}) {
      auto coefficient=[&](const char* text){return Exact(field,text);};
      AffineFrobeniusSeries::Options limited;limited.max_terms=40;limited.univariate_epsilon_projection=mode==1;limited.cleared_epsilon_projection=mode==2;limited.finite_lag_projection=mode==3;
      AffineFrobeniusSeries::Matrix diagonal{{coefficient("1/(1-x)"),coefficient("0")},
        {coefficient("0"),coefficient("1/(1-x)")}};
      auto f=AffineFrobeniusSeries::prepare(diagonal,0,1,8,limited);
      // Ninety convolution products combine into eighteen coordinates. A
      // bounded accumulator must fit although the old temporary list did not.
      auto projected=f.project(std::vector<Exact>{coefficient("1/(1-x)"),coefficient("1/(1-x)")});
      check(projected.terms.size()==18,"streamed projection failed receiving coordinate budget");
      for(const auto& term:projected.terms)check(term.coefficient==coefficient("1").constant(term.power+Rational(1)),
        "streamed projection differs from exact (1-x)^-2 coefficients");
      limited.max_projection_products=1;auto bounded=AffineFrobeniusSeries::prepare(diagonal,0,1,8,limited);
      bool rejected=false;try{bounded.project(std::vector<Exact>{coefficient("1/(1-x)"),coefficient("1/(1-x)")});}
      catch(const std::length_error&){rejected=true;}check(rejected,"projection work budget was not enforced");
      limited.max_projection_products=2;auto endpoint=AffineFrobeniusSeries::prepare(diagonal,0,1,8,limited);
      auto domain=endpoint.project_endpoint_domain({{coefficient("1/(1-x)"),coefficient("1/(1-x)")}});
      check(domain.zero_constraints.empty() && domain.admissible.terms.size()==2,
        "endpoint projection spent its work budget on positive powers");
      auto constants=endpoint.dr_endpoint_constant(domain.admissible);
      check(constants[0][0]==coefficient("1") && constants[0][1]==coefficient("1"),"endpoint-only projection changed exact limits");
    }
    auto e = [&](const char *s) { return Exact(field, s); };
    {
      // Compare direct FLINT coefficient ingestion with the independent text
      // evaluator, including numerator zeros, epsilon poles and large powers.
      for(const auto& value:std::vector<Exact>{e("7/3"),e("eps^4*(2+eps)/(3-eps)"),e("(1+eps)^60/(eps^3*(1-2*eps)^50)")}){
        affine_frobenius_detail::EpsilonRational native(value,1);
        auto valuation=fuchsify::detail::valuation(value,1);
        check(native.valuation()==valuation,"native epsilon valuation changed");
        Jet ep(0,10,B::precision());ep.set(1,B(1));
        auto regular=value/fuchsify::detail::power(e("eps"),valuation);
        auto parsed=diffexp::evaluate(data::Reader(regular.str()).read(),ep,{{"eps",ep}});
        auto direct=native.regular_jet(10);
        for(unsigned k=0;k<10;++k)near(direct.at(k),parsed.at(k));
      }
      bool rejected=false;try{affine_frobenius_detail::EpsilonRational nonunivariate(e("x+eps"),1);}
      catch(const std::domain_error&){rejected=true;}check(rejected,"unsubstituted coefficient parameter accepted by native epsilon evaluator");
    }
    {
      auto compare=[](const auto& a,const auto& b){
        check(a.rows==b.rows&&a.columns==b.columns&&a.terms.size()==b.terms.size()&&a.wronskian_prefactor==b.wronskian_prefactor&&a.coherent_x_frontier==b.coherent_x_frontier,"projection metadata mismatch");
        for(std::size_t i=0;i<a.terms.size();++i){const auto& x=a.terms[i];const auto& y=b.terms[i];check(x.row==y.row&&x.column==y.column&&x.log_degree==y.log_degree&&x.power==y.power&&x.slope==y.slope&&x.coefficient==y.coefficient,"univariate projection exact coefficient mismatch");}
      };
      for(unsigned mode:{1,2,3})for(const auto& a:std::vector<AffineFrobeniusSeries::Matrix>{
          {{e("eps/x"),e("1/x")},{e("0"),e("eps/x+(x+eps)/(1-x)")}},
          {{e("0"),e("0")},{e("1/(1-x)"),e("1/x")}}}) {
        AffineFrobeniusSeries::Options options;auto generic=AffineFrobeniusSeries::prepare(a,0,1,8,options);
        options.univariate_epsilon_projection=mode==1;options.cleared_epsilon_projection=mode==2;options.finite_lag_projection=mode==3;auto fast=AffineFrobeniusSeries::prepare(a,0,1,8,options);
        AffineFrobeniusSeries::Matrix rows{{e("(1+eps)/(x^2*(1+x))"),e("(x+eps)/(3*x*(1-x))")},{e("1/(2-eps)"),e("0")}};
        compare(generic.project(rows),fast.project(rows));
        auto g=generic.project_endpoint_domain(rows),f=fast.project_endpoint_domain(rows);compare(g.admissible,f.admissible);
        check(g.zero_constraints.size()==f.zero_constraints.size(),"univariate projection domain size mismatch");
        for(std::size_t i=0;i<g.zero_constraints.size();++i){const auto& x=g.zero_constraints[i];const auto& y=f.zero_constraints[i];check(x.row==y.row&&x.power==y.power&&x.log_degree==y.log_degree&&x.coefficients==y.coefficients,"univariate projection domain mismatch");}
      }
    }
    Exact x = e("x");
    auto half = B::from_strings("0.5");
    B logx;
    acb_log(logx.raw(), half.raw(), B::precision());
    auto scalar = AffineFrobeniusSeries::prepare({{e("eps/x")}}, 0, 1, 3);
    {ExactField foreign({"x","eps"});auto corrupted=scalar.terms();corrupted.terms.front().coefficient=Exact(foreign,"eps");bool caught=false;try{scalar.dr_domain(corrupted,false);}catch(const std::invalid_argument&){caught=true;}check(caught,"normalization accepted a coefficient from another exact field");}
    check(scalar.dr_endpoint_constant(scalar.terms())[0][0].is_zero(),
          "x^eps must not be confused with a constant endpoint sector");
    auto integral = scalar.dr_integral_from_zero(
        scalar.project(std::vector<Exact>{e("1/x")}));
    check(integral.terms.size() == 1 &&
              integral.terms[0].coefficient == e("1/eps"),
          "DR x^(eps-1) integral retains 1/eps");
    auto metadata = scalar.valuation_metadata(integral);
    check(metadata.maximum_pole == 1 &&
              metadata.required_source_top(4) == std::vector<long>{5},
          "primitive source lookahead");
    auto expansion = scalar.evaluate_laurent(integral, half, -1, 3);
    near(expansion.coefficients[0][0][0], B(1));
    near(expansion.coefficients[0][0][1], logx);
    near(expansion.coefficients[0][0][2], logx * logx / B(2));
    auto constant = AffineFrobeniusSeries::prepare({{e("0")}}, 0, 1, 2);
    check(constant.dr_endpoint_constant(constant.terms())[0][0] == e("1"),
          "constant DR endpoint sector preserved");
    rejects(
        [&] {
          constant.dr_integral_from_zero(
              constant.project(std::vector<Exact>{e("1/x")}));
        },
        "fixed unregulated 1/x must reject");
    rejects(
        [&] {
          constant.dr_integral_from_zero(
              constant.project(std::vector<Exact>{e("1/x^2")}));
        },
        "fixed unregulated power divergence must reject");
    // Equal-slope rational Jordan residue: x^(-1+2eps) log(x).
    auto jordan = AffineFrobeniusSeries::prepare(
        {{e("(-1+2*eps)/x"), e("1/x")}, {e("0"), e("(-1+2*eps)/x")}}, 0, 1, 2);
    auto logsector = jordan.contract(
        jordan.project(std::vector<Exact>{e("1"), e("0")}), {e("0"), e("1")});
    auto moment = jordan.dr_integral_from_zero(logsector);
    auto atone = jordan.evaluate_laurent(moment, B(1), -2, 2);
    near(atone.coefficients[0][0][0], -B(1) / B(4));
    for (unsigned k = 1; k < 5; ++k)
      near(atone.coefficients[0][0][k], B(0));
    check(jordan.valuation_metadata(moment).maximum_pole == 2,
          "log moment requires two epsilon lookahead slots");
    // Colliding eigenvectors create 1/eps coefficients, which must cancel only
    // after preserving both distinct affine exponents through evaluation.
    auto collision = AffineFrobeniusSeries::prepare(
        {{e("0"), e("1/x")}, {e("0"), e("eps/x")}}, 0, 1, 2);
    check(collision.residue_frame()[0][1] == e("1/eps"),
          "exact colliding affine frame");
    auto colliding =
        collision.contract(collision.terms(), {e("-1/eps"), e("1")});
    auto ce = collision.evaluate_laurent(colliding, half, -1, 3);
    near(ce.coefficients[0][0][0], B(0));
    near(ce.coefficients[0][0][1], logx);
    near(ce.coefficients[0][0][2], logx * logx / B(2));
    near(ce.coefficients[1][0][1], B(1));
    // Fixed integer resonance yields an exact log independently of epsilon
    // poles.
    auto resonant = AffineFrobeniusSeries::prepare(
        {{e("eps/x"), e("0")}, {e("1"), e("(1+eps)/x")}}, 0, 1, 3);
    auto rv =
        resonant.evaluate(resonant.terms(), half, B::from_strings("0.25"));
    B quarter;
    auto qlog = logx / B(4);
    acb_exp(quarter.raw(), qlog.raw(), B::precision());
    near(rv[1][0], half * quarter * logx);
    rejects(
        [&] {
          AffineFrobeniusSeries::prepare(
              {{e("0"), e("1/x")}, {e("eps/x"), e("0")}}, 0, 1, 2);
        },
        "nonlinear algebraic epsilon spectrum rejected");
    rejects([&] { AffineFrobeniusSeries::prepare({{e("eps^2/x")}}, 0, 1, 2); },
            "nonaffine rational epsilon spectrum rejected");
    // Near-resonance at positive x order introduces a genuine coefficient pole.
    auto pseudo = AffineFrobeniusSeries::prepare(
        {{e("0"), e("0")}, {e("1"), e("(1+eps)/x")}}, 0, 1, 3);
    check(pseudo.valuation_metadata(pseudo.terms()).maximum_pole == 1,
          "positive-order affine collision lost its epsilon pole");
    auto ps = pseudo.contract(pseudo.terms(), {e("1"), e("1/eps")});
    auto pe = pseudo.evaluate_laurent(ps, half, -1, 2);
    near(pe.coefficients[1][0][0], B(0));
    near(pe.coefficients[1][0][1], half * logx);
    // Exact cancellation precedes classification of a fixed divergent sector.
    auto pair = AffineFrobeniusSeries::prepare(
        {{e("0"), e("0")}, {e("0"), e("0")}}, 0, 1, 2);
    auto cancel =
        pair.contract(pair.project(std::vector<Exact>{e("1/x"), e("-1/x")}),
                      {e("1"), e("1")});
    check(pair.dr_integral_from_zero(cancel).terms.empty(),
          "exact cancelled fixed divergence was rejected");
    auto rational_jordan = AffineFrobeniusSeries::prepare(
        {{e("eps/x"), e("1/((1-eps)*x)")}, {e("0"), e("eps/x")}}, 0, 1, 1);
    auto rj = rational_jordan.evaluate(rational_jordan.terms(), half,
                                       B::from_strings("0.25"));
    near(rj[0][1], quarter * logx * B(4) / B(3));
    if (argc > 1) {
      level::Options options;
      options.provider.executable = argv[1];
      auto d = e("2-2*eps");
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
        auto prepared = level::prepare(basis, d, field, 0, {target}, options);
        check(prepared.success, prepared.reason.c_str());
        for (unsigned endpoint = 0; endpoint < 2; ++endpoint) {
          auto a = prepared.matrix;
          if (endpoint)
            for (auto &row : a)
              for (auto &v : row)
                v = -v.substitute(std::vector<Exact>{e("1-x"), e("eps")});
          auto f = fuchsify::prepare(a, 0);
          check(f.success, f.reason.c_str());
          auto affine = AffineFrobeniusSeries::prepare(f.matrix, 0, 1, 3);
          std::cout << "Banana" << loops << " endpoint=" << endpoint
                    << " affine roots:";
          for (auto &r : affine.exponents())
            std::cout << " (" << r.power.str() << "," << r.slope.str() << ")";
          std::cout << "\n";
        }
      }
    }
    std::cout << "Affine exact spectra, DR primitives, Laurent cancellation "
                 "and lookahead passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
