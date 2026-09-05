#pragma once
#include "diffexp/affine_matching.hpp"
#include "diffexp/reference_comparison.hpp"
#include <boost/json.hpp>

namespace diffexp {
// Native counterpart of Examples/Direct/SingularEndpointAndSegments.wl.
// A single exact monomial chart covers this solution; no subdivision is needed.
inline boost::json::object singular_endpoint_example() {
  namespace json=boost::json;using B=kernel::ComplexBall;
  struct PrecisionScope {slong before=B::precision();PrecisionScope(){B::set_precision(384);}~PrecisionScope(){B::set_precision(before);}} precision;
  ExactField field({"x","eps"});Exact x(field,"x"),eps(field,"eps");
  auto series=AffineFrobeniusSeries::prepare({{eps/x}},0,1,40);
  const auto& expansion=series.terms();
  if(expansion.terms.size()!=1 || expansion.terms[0].power!=Rational(0) ||
      expansion.terms[0].slope!=Rational(1) || expansion.terms[0].log_degree ||
      expansion.terms[0].coefficient!=x.constant(1))
    throw std::runtime_error("singular direct example is not the exact x^eps sector");
  // Verify the retained expression solves the full connection identically.
  if(!(x.constant(expansion.terms[0].power)+eps*x.constant(expansion.terms[0].slope)==eps))
    throw std::runtime_error("singular direct monomial differential identity");
  B anchor=B::from_strings("1/2"),log_anchor;acb_log(log_anchor.raw(),anchor.raw(),384);
  Jet epsilon(0,4,384);epsilon.set(1,B(1));auto log_jet=epsilon.constant(0);log_jet.set(0,log_anchor);
  auto initial=(epsilon*log_jet).exp();
  affine_matching::Boundary boundary{0,3,{std::vector<B>(4)}};
  for(unsigned k=0;k<4;++k)boundary.coefficients[0][k]=initial.at(k);
  auto matched=affine_matching::match(series,anchor,boundary,{0,3});
  if(!matched.success())throw std::runtime_error("singular direct matching: "+matched.reason);
  auto endpoint=series.dr_endpoint_constant(expansion);
  if(!endpoint[0][0].is_zero())throw std::runtime_error("x^eps dimensional endpoint must be zero");
  B log_lower,log_upper;auto lower=B::from_strings("1/1000000"),upper=B::from_strings("1/4");
  acb_log(log_lower.raw(),lower.raw(),384);acb_log(log_upper.raw(),upper.raw(),384);
  json::array samples;double max_error=0;
  const auto text=[](const B& value){char* raw=arb_get_str(acb_realref(value.raw()),45,0);std::string out(raw);flint_free(raw);return out;};
  for(unsigned i=0;i<=80;++i) {
    auto logarithm=log_lower+(log_upper-log_lower)*B(i)/B(80);B at;acb_exp(at.raw(),logarithm.raw(),384);
    auto value=affine_matching::apply(series,expansion,at,matched.value,{0,3});
    if(!value.success())throw std::runtime_error("singular direct sample: "+value.reason);
    max_error=std::max(max_error,finite_reference_error(value.value.coefficients[0][1],logarithm,384));
    samples.push_back(json::object{{"x",text(at)},{"epsilon_one",text(value.value.coefficients[0][1])},
      {"plot_x",arf_get_d(arb_midref(acb_realref(at.raw())),ARF_RND_NEAR)},
      {"plot_y",arf_get_d(arb_midref(acb_realref(value.value.coefficients[0][1].raw())),ARF_RND_NEAR)}});
  }
  if(!(max_error<1e-60))throw std::runtime_error("singular direct epsilon-one samples disagree with logarithm");
  return {{"schema","DiffExp3.SingularEndpointExample/v1"},{"epsilon_low",0},{"epsilon_high",3},
    {"segments",json::array{json::object{{"center","0"},{"scale","1"},{"incoming_match_point","1/2"},
      {"singular",true},{"radius",nullptr},{"radius_reason","exact monomial solution on the positive real branch"}}}},
    {"sectors",json::array{json::object{{"a","0"},{"b","1"},{"p",0}}}},
    {"dimreg_endpoint_constant","0"},{"ordinary_endpoint_limit_claimed",false},
    {"exact_monomial_differential_identity_verified",true},{"max_epsilon_one_error",max_error},
    {"samples",std::move(samples)}};
}
} // namespace diffexp
