#include "diffexp/exact_jet_composition.hpp"
#include <iostream>
using namespace diffexp;
namespace jc=exact_jet_composition;
void require(bool value,const char* reason){if(!value)throw std::runtime_error(reason);}
template<class F>void rejects(F f,const char* reason){try{f();}catch(const std::exception& e){require(std::string(e.what()).find(reason)!=std::string::npos,"unexpected exact composition rejection");return;}throw std::runtime_error("invalid exact composition accepted");}
int main(){try {
  ExactField field({"x","eps"});Exact x(field,"x"),e(field,"eps"),z(field,0),o(field,1),two(field,2),three(field,3);
  // Unrelated pole directions: a sum of global minima falsely predicts loss 3.
  auto diagonal=jc::fuse({{o/e,z},{z,o}},{{o,z},{z,o/e.pow(2)}},1);
  require(diagonal.separate_global_lower_sum==jc::High(-3) && diagonal.pathwise_lower==jc::High(-2) &&
    diagonal.fused_lower==jc::High(-2),"unrelated pole directions were globally charged");
  require(diagonal.product==jc::Matrix{{o/e,z},{z,o/e.pow(2)}},"incorrect exact diagonal product");
  auto schedule=jc::schedule(diagonal.product,{0,0},0,1);
  require(schedule.input_highs==jc::Highs{1,2},"per-input pullback did not retain distinct highs");
  require(schedule.entries.size()==2 && schedule.map_coefficients==5,"sparse exact map expansion not scheduled");
  auto selected=jc::schedule(diagonal.product,{0,0},jc::Highs{std::nullopt,0},1);
  require(!selected.input_highs[0] && selected.input_highs[1]==jc::High(2),"unrequested output created demand");
  auto value=jc::apply(diagonal.product,{{0,-1,{}},{0,2,{o,two,three}}},jc::Highs{std::nullopt,0},1);
  require(!value.outputs[0] && value.outputs[1]->low==-2 && value.outputs[1]->coefficients==std::vector<Exact>{o,two,three},"componentwise output jet incorrect");

  // Sum cancellation must precede epsilon valuation and truncation.
  jc::Matrix outer{{o/e,-o/e}},inner{{o},{o+e}};
  auto cancellation=jc::fuse(outer,inner,1);
  require(cancellation.separate_global_lower_sum==jc::High(-1) && cancellation.pathwise_lower==jc::High(-1) &&
    cancellation.fused_lower==jc::High(0) && cancellation.product[0][0]==-o,"exact cancellation was lost");
  auto result=jc::apply(cancellation.product,{{0,0,{three}}},0,1);
  require(result.outputs[0]->coefficients==std::vector<Exact>{-three},"fused application requested hidden high coefficient");
  // Separately truncated inner rows really are incomplete; identical retained
  // constants must not be interpreted as an exact identity of their unknown tails.
  bool missing=false;
  try{(void)jc::apply(outer,{{0,0,{three}},{0,0,{three}}},0,1);}
  catch(const jc::MissingCoefficient& d){missing=d.input==0 && d.required_high==1 && d.available_high==0;}
  require(missing,"unknown higher jet coefficient silently used as zero");
  auto zero=jc::fuse(outer,{{o},{o}},1);
  require(!zero.fused_lower && zero.product[0][0].is_zero(),"exact zero product must have infinite valuation");
  auto empty=jc::apply(zero.product,{{0,-1,{}}},0,1);
  require(!empty.demand.input_highs[0] && empty.outputs[0]->coefficients[0].is_zero(),"exact zero product demanded unknown source");

  auto expanded=jc::expand_demanded({{o/(o-e)}},{0},jc::Highs{2},1);
  require(expanded.size()==1 && expanded[0].coefficients.coefficients==std::vector<Exact>{o,o,o},"exact rational epsilon expansion");
  auto rational=jc::apply({{o/(o-e)}},{{0,2,{two,three,z.constant(5)}}},2,1);
  require(rational.outputs[0]->coefficients==std::vector<Exact>{two,z.constant(5),z.constant(10)},"finite convolution is incorrect");
  auto symbolic=jc::apply({{(o+x*e)/(o-e)}},{{0,2,{o,z,z}}},2,1);
  require(symbolic.outputs[0]->coefficients==std::vector<Exact>{o,o+x,o+x},"non-epsilon coefficient field changed");
  auto signed_window=jc::apply({{e}},{{-2,0,{o,two,three}}},-1,1);
  require(signed_window.outputs[0]->low==-1 && signed_window.outputs[0]->coefficients==std::vector<Exact>{o},"signed Laurent lower bound lost");
  // A product beginning above the requested high is already known zero, even
  // when the first potentially nonzero source coefficient is unavailable.
  auto beyond=jc::apply({{e.pow(3)}},{{0,-1,{}}},1,1);
  require(!beyond.demand.input_highs[0] && beyond.outputs[0]->coefficients[0].is_zero(),"known lower support ignored");

  rejects([&]{jc::fuse_available(std::nullopt,inner,1);},"requires both exact maps");
  rejects([&]{jc::fuse({},inner,1);},"missing exact map");
  rejects([&]{jc::fuse({{o,z}},{{o}},1);},"adjacent map dimensions");
  rejects([&]{jc::fuse({{o},{o,z}},{{o}},1);},"ragged");
  rejects([&]{jc::schedule({{o}},{0},jc::Highs{},1);},"demand shape");
  rejects([&]{jc::apply({{o}},{{0,1,{o}}},0,1);},"input jet shape");
  rejects([&]{jc::apply({{o}},{{0,0,{e}}},0,1);},"depends on epsilon");
  jc::Limits small;small.max_products=1;
  rejects([&]{jc::fuse({{o,o}},{{o},{o}},1,small);},"product budget");
  small={};small.max_term_work=1;
  rejects([&]{jc::fuse({{o}},{{o}},1,small);},"term-work budget");
  small={};small.max_coefficients=1;
  rejects([&]{jc::schedule({{o}},{0},1,1,small);},"coefficient budget");
  small={};small.max_abs_order=1;
  rejects([&]{jc::schedule({{o/e.pow(2)}},{0},0,1,small);},"order budget");
  ExactField other({"eps","x"});
  rejects([&]{jc::fuse({{o}},{{Exact(other,1)}},1);},"field");
  std::cout<<"Exact adjacent fusion, cancellation-aware component schedules, sparse rational expansion and unknown-jet rejection passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
