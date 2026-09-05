#include "diffexp/recursion_leaf.hpp"
#include "diffexp/certified_bubble.hpp"
#include "diffexp/sunrise.hpp"
#include <iostream>
using namespace diffexp;
double size(const Jet::Ball& x) {arf_t u;arf_init(u);acb_get_abs_ubound_arf(u,x.raw(),384);auto d=arf_get_d(u,ARF_RND_CEIL);arf_clear(u);return d;}
int main(){try {
  using B=Jet::Ball;B::set_precision(384);
  auto graph=recursion::prepare(feynman::example_family("bubble"),{});
  feynman::CertifiedDeepestBetaOptions proof;proof.working_bits=384;proof.requested_digits=28;
  auto bubble=recursion::evaluate_leaf(graph.nodes[0],graph.dimension,2,4,proof);
  CertifiedBubbleOptions bubble_options;bubble_options.working_bits=384;bubble_options.requested_digits=28; bubble_options.taylor_order=160;
  auto reference=certified_feynman_bubble(bubble_options);
  for(unsigned k=0;k<=4;++k)if(size(bubble.values[0][k-bubble.epsilon_low]-reference.coefficients[k])>1e-27)
    throw std::runtime_error("generic scalar leaf disagrees with independently derived bubble system");
  ExactField field({"x","eps","I"});Exact x(field,"x"),dimension(field,"2-2*eps");
  const Rational anchor("1/3");
  auto raw=ibp::merge(ibp::quadratic_family(feynman::banana(2,{Rational(1),Rational(1),Rational(1)}),x),0,1,x.constant(anchor));
  recursion::Node leaf{ibp::PropagatorBasis(raw),ibp::PropagatorBasis(ibp::merge(raw,0,1,x)),0,1,Rational("1/3"),
    {{1,0,0,0,0},{1,1,0,0,0},{1,1,-1,0,0},{0,1,0,0,0}}, {}, {}, {}, {},true};
  for(const auto& target:leaf.requested)leaf.operations.push_back(pullback::plan(leaf.incoming,leaf.merged,0,1,target,x));
  auto actual=recursion::evaluate_leaf(leaf,dimension,2,4,proof);
  auto expected=feynman::sunrise_boundary(feynman::prepare_sunrise(),anchor,4);
  double maximum=0;
  for(unsigned i=0;i<3;++i)for(int k=std::max(actual.epsilon_low,expected.epsilon_low);k<=4;++k)
    maximum=std::max(maximum,size(actual.values[i][k-actual.epsilon_low]-expected.values[i][k-expected.epsilon_low]));
  if(maximum>1e-27)throw std::runtime_error("generic Gaussian numerator leaf disagrees with independent sunrise boundary coordinates");
  for(const auto& value:actual.values[3])if(!value.is_zero())throw std::runtime_error("proven free-loop endpoint is not zero");
  if(!actual.taylor_tail_certified)throw std::runtime_error("generic leaf lost scalar certificates");
  std::cout<<"Generic beta/numerator/endpoint scalar leaf discrepancy "<<maximum<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
