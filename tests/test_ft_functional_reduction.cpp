#include "diffexp/ft_functional_reduction.hpp"
#include "diffexp/direct_adjoint_endpoint.hpp"
#include <iostream>
using namespace diffexp;
namespace fr=ft_functional_reduction;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
template<class F>void rejects(F f,const char* reason){try{f();}catch(const std::exception& e){require(std::string(e.what()).find(reason)!=std::string::npos,"unexpected certificate rejection");return;}throw std::runtime_error("malformed certificate accepted");}
int main(){try {
  ExactField field({"x","eps"});Exact x(field,"x"),e(field,"eps"),z(field,0),o(field,1),two(field,2);
  require(direct_adjoint_endpoint::exact_residual_zero({{x}},{{z}},{{o}},0,1),"exact primitive rejected");
  require(!direct_adjoint_endpoint::exact_residual_zero({{x+x.pow(2)/two}},{{z}},{{o/(o-x)}},0,1),
    "finite spatial primitive was promoted to exact endpoint data");
  fr::Functional first{{{e}},{{-o/e}},{{z}},{{o/e}},z,o};
  auto a=fr::search_counterterm(first,o);
  require(a.status==fr::SearchStatus::certified && a.certificate.has_value(),"first-order endpoint cancellation not found");
  require(a.certificate->rewritten.lower[0][0].is_zero() && a.certificate->rewritten.upper[0][0].is_zero() &&
    a.certificate->rewritten.weight[0][0]==o,"(Y(1)-Y(0))/eps must become integral Y");
  fr::verify_counterterm(first,*a.certificate);
  auto wrong=*a.certificate;wrong.rewritten.weight[0][0]=two;
  rejects([&]{fr::verify_counterterm(first,wrong);},"differential identity");
  wrong=*a.certificate;wrong.rewritten.upper_point=two;
  rejects([&]{fr::verify_counterterm(first,wrong);},"metadata mismatch");

  // Arbitrary-lookahead example d=2: a regular weight 1-x replaces epsilon^-2 endpoints.
  fr::Functional second{{{e}},{{-(o+e)/e.pow(2)}},{{z}},{{o/e.pow(2)}},z,o};
  auto b=fr::search_counterterm(second,o);
  require(b.status==fr::SearchStatus::certified,"second-order principal part not removed");
  require(b.certificate->rewritten.weight[0][0]==o-x,"second-order Taylor remainder weight");
  require(b.certificate->rewritten.lower[0][0].is_zero() && b.certificate->rewritten.upper[0][0].is_zero(),"second-order endpoints remain");
  fr::SearchOptions shallow;shallow.pole_order=1;
  require(fr::search_counterterm(second,o,shallow).status==fr::SearchStatus::no_certificate_in_ansatz,"bounded ansatz failure misreported");

  // Rational epsilon denominators are cleared exactly, rather than expanded at a guessed order.
  auto rational_epsilon=first;rational_epsilon.connection={{e/(o+e)}};
  auto re=fr::search_counterterm(rational_epsilon,o);
  require(re.status==fr::SearchStatus::certified &&
    re.certificate->rewritten.weight[0][0]==o/(o+e),"rational epsilon denominator certificate");
  auto multiple=first;multiple.lower.push_back(second.lower[0]);multiple.weight.push_back({z});multiple.upper.push_back(second.upper[0]);
  auto multi=fr::search_counterterm(multiple,o);
  require(multi.status==fr::SearchStatus::certified && multi.certificate->rewritten.weight[1][0]==o-x,"multiple independent functionals");

  // A supplied allowed spatial divisor enables a rational counterterm unavailable to a constant ansatz.
  fr::Functional rational{{{z}},{{-o/e}},{{o/(e*(o+x).pow(2))}},{{o/(two*e)}},z,o};
  fr::SearchOptions constant;constant.pole_order=1;constant.spatial_degree=0;
  require(fr::search_counterterm(rational,o,constant).status==fr::SearchStatus::no_certificate_in_ansatz,"polynomial ansatz falsely solved rational primitive");
  auto c=fr::search_counterterm(rational,o+x,constant);
  require(c.status==fr::SearchStatus::certified && c.certificate->rewritten.weight[0][0].is_zero(),"allowed rational denominator not used");

  // Regular bulk is insufficient: two incompatible endpoint residues cannot both be removed.
  fr::Functional incompatible{{{z}},{{z}},{{z}},{{o/e}},z,o};
  require(fr::search_counterterm(incompatible,o).status==fr::SearchStatus::no_certificate_in_ansatz,"incompatible endpoints accepted");
  rejects([&]{fr::certify_counterterm(incompatible,{{z}});},"epsilon poles");
  // Scalar epsilon-regular logarithmic connection; its endpoint information still has an epsilon pole.
  fr::Functional endpoint_pole{{{e/x}},{{z}},{{z}},{{o/e}},o,two};
  require(fr::search_counterterm(endpoint_pole,o).status==fr::SearchStatus::no_certificate_in_ansatz,"regular bulk falsely removed endpoint pole");
  auto singular_cut=rational;singular_cut.lower_point=-o;
  require(fr::search_counterterm(singular_cut,o+x,constant).status==fr::SearchStatus::unsupported,"undefined counterterm endpoint accepted");
  auto nonregular=first;nonregular.connection={{o/e}};
  require(fr::search_counterterm(nonregular,o).status==fr::SearchStatus::unsupported,"irregular connection accepted by search");
  auto tiny=constant;tiny.max_unknowns=1;tiny.spatial_degree=1;
  require(fr::search_counterterm(first,o,tiny).status==fr::SearchStatus::budget_exceeded,"unknown budget ignored");
  tiny=constant;tiny.max_terms=1;
  require(fr::search_counterterm(first,o,tiny).status==fr::SearchStatus::budget_exceeded,"term budget ignored");
  auto malformed=first;malformed.weight[0].push_back(z);
  rejects([&]{fr::search_counterterm(malformed,o);},"column count");
  rejects([&]{fr::search_counterterm(first,e);},"spatial denominator");

  // Rectangular U=(Y,(Y-S)/eps) drops an unused old state. Initialization is
  // simplified symbolically, so U(0)=(s,0), without computing s through H+1.
  fr::Matrix M{{e,z,z},{z,z,z},{z,z,z}};
  fr::Matrix R{{o,z,z},{o/e,z,-o/e}};
  fr::Matrix Ahat{{e,z},{o,z}};
  fr::Segment segment{M,R,Ahat,z,o};
  fr::RectangularRealization realization{{segment},{},{{o},{z},{o}},{{o},{z}},{{o/e,z,-o/e}},{{z,o}}};
  fr::verify_regular_realization(realization);
  // A nontrivial evaluated junction: reset Y and anchor S by the same factor.
  auto chain=realization;chain.segments.push_back({M,R,Ahat,o,two});
  chain.interfaces.push_back({{{two,z,z},{z,o,z},{z,z,two}},{{two,z},{z,two}}});
  fr::verify_regular_realization(chain);
  auto bad=chain;bad.interfaces[0].reduced[0][0]=o;
  rejects([&]{fr::verify_regular_realization(bad);},"interface identity");
  bad=realization;bad.segments[0].reduced_connection[1][0]=two;
  rejects([&]{fr::verify_regular_realization(bad);},"differential identity");
  bad=realization;bad.reduced_initialization[1][0]=o;
  rejects([&]{fr::verify_regular_realization(bad);},"initialization identity");
  bad=realization;bad.reduced_output[0][0]=o;
  rejects([&]{fr::verify_regular_realization(bad);},"output identity");
  bad=realization;bad.segments[0].projection[0].pop_back();
  rejects([&]{fr::verify_regular_realization(bad);},"column count");
  bad=realization;bad.interfaces.push_back(chain.interfaces[0]);
  rejects([&]{fr::verify_regular_realization(bad);},"segment/interface count");

  // Coordinate derivatives and evaluated interface cuts must not be dropped.
  fr::Segment varying{{{z}},{{x}},{{o/x}},o,two};
  fr::RectangularRealization varying_chain{{varying,{{{z}},{{x}},{{o/x}},z.constant(3),z.constant(4)}},
    {{{{o}},{{z.constant(Rational("3/2"))}}}},{{o}},{{o}},{{z.constant(4)}},{{o}}};
  fr::verify_regular_realization(varying_chain);
  auto varying_bad=varying_chain;varying_bad.segments[0].reduced_connection={{z}};
  rejects([&]{fr::verify_regular_realization(varying_bad);},"differential identity");
  varying_bad=varying_chain;varying_bad.interfaces[0].reduced={{o}};
  rejects([&]{fr::verify_regular_realization(varying_bad);},"interface identity");

  // Genuine information obstruction: (s1-s2)/eps cannot be initialized from
  // independent source jets of the same order, even though its bulk equation is zero.
  fr::RectangularRealization obstruction{{{{{z,z},{z,z}},{{o/e,-o/e}},{{z}},z,o}},
    {},{{o,z},{z,o}},{{o/e,-o/e}},{{o/e,-o/e}},{{o}}};
  rejects([&]{fr::verify_regular_realization(obstruction);},"initialization/output has epsilon pole");
  // Giving the divided difference as a NEW declared independent source changes
  // the information contract; only then is this initializer regular.
  obstruction.source_initialization={{o,e},{o,z}};obstruction.reduced_initialization={{z,o}};
  fr::verify_regular_realization(obstruction);

  // Exact rational field context and variable ordering are part of the certificate.
  ExactField other({"eps","x"});bad=realization;bad.output[0][0]=Exact(other,1);
  rejects([&]{fr::verify_regular_realization(bad);},"field");
  std::cout<<"Whole-functional exact counterterms, bounded rational search, rectangular regular realizations, endpoint obstructions and malformed certificates passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
