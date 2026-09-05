#include "diffexp/system.hpp"
#include "diffexp/scheduler.hpp"
#include "diffexp/linear_uncertainty.hpp"
#include "diffexp/feynman.hpp"
#include <iostream>

using namespace diffexp;
void check(bool b,const char* message) { if (!b) throw std::runtime_error(message); }
template<class F> void rejects(F f,const char* message) {
  bool rejected=false; try {f();} catch(const std::exception&) {rejected=true;} check(rejected,message);
}

void exact_tests() {
  ExactField field({"x","eps"}), other({"s","d"});
  Exact x(field,"x"),e(field,"eps"),s(other,"s");
  check(((x*x-x.constant(1))/(x-x.constant(1)))==x+x.constant(1),"rational normalization");
  check(Exact(field,"(x+eps)/(1-x)").derivative(0)==Exact(field,"(1+eps)/(1-x)^2"),"exact derivative");
  check(Exact(field,"x+eps").substitute(std::vector<Exact>{e,x})==x+e,"simultaneous substitution");
  check(Exact(field,"x^2-1").polynomial_lcm(Exact(field,"x-1"))==Exact(field,"x^2-1"),"denominator lcm");
  check((s+s)==Exact(other,"2*s"),"independent fields coexist");
  rejects([&]{(void)(x+s);},"mixed field rejected");
  rejects([&]{(void)(x/x.constant(0));},"zero division rejected");
  rejects([&]{(void)Exact(field,"x/(eps-eps)");},"input division by zero rejected before FLINT");
  auto retained=[] {ExactField f({"z"});return Exact(f,"z^2");}();
  check(retained.derivative(0).str()=="2*z","value owns field lifetime");
}

void recurrence_tests() {
  ExactField f({"x","eps"});
  // y'=y/(x-2), y(0)=1/2: exact y=1/2-x/4, no truncation ambiguity.
  RationalSystem linear({{Exact(f,"1/(x-2)")}},0,1);
  auto y=regular_series(linear,Rational(0),8,0,0,{{Rational("1/2")}});
  check(y.at(1,0,0)==Rational("-1/4"),"compiler produces regular recurrence");
  for(unsigned n=2;n<=8;++n) check(y.at(n,0,0).is_zero(),"linear exact tail vanishes");
  check(y.evaluate_polynomial(Rational(1))[0][0]==Rational("1/4"),"standalone minimal example");
  rejects([&]{regular_series(linear,Rational(2),4,0,0,{{Rational(1)}});},"singular center cannot masquerade as ordinary");
  // y' = eps*y/(1-x), y(0)=1. epsilon^1 is -log(1-x).
  RationalSystem logarithm({{Exact(f,"eps/(1-x)")}},0,1);
  auto log=regular_series(logarithm,Rational(0),12,0,3,{{Rational(1),Rational(0),Rational(0),Rational(0)}});
  for(unsigned n=1;n<=12;++n) check(log.at(n,0,1)==Rational(1)/Rational(n),"canonical logarithm recurrence");
  check(log.at(2,0,2)==Rational("1/2"),"epsilon triangular convolution");
  // Nonpolynomial epsilon dependence is retained exactly until recurrence.
  RationalSystem rational_epsilon({{Exact(f,"1/(1+eps)")}},0,1);
  auto exp=regular_series(rational_epsilon,Rational(0),3,0,3,{{Rational(1),Rational(0),Rational(0),Rational(0)}});
  check(exp.at(2,0,2)==Rational("3/2"),"rational epsilon denominator");
}

void h48_test() {
  using X=Executor<int>; X executor;
  Demand upstream{0,22,25,256,8},basis{0,21,25,256,8};
  for(int level=7;level>=4;--level) {
    auto id="level-"+std::to_string(level);
    std::vector<X::Dependency> parents;
    if(level<7) parents.push_back({"level-"+std::to_string(level+1),[=](const Demand&){return upstream;}});
    executor.add({id,id,parents,[=](const Demand& d,const auto&) {return X::Artifact{d,level,id,{}};}});
  }
  for(const std::string arm:{"lower","upper"}) {
    auto id="basis-"+arm;
    executor.add({id,id,{},[=](const Demand& d,const auto&) {return X::Artifact{d,d.epsilon_high,id,{}};}});
    auto match="match-"+arm;
    executor.add({match,match,{{"level-4",[=](const Demand&){return upstream;}},{id,[=](const Demand&){return basis;}}},
      [=](const Demand& d,const auto& parents) {
        if(parents[1]->guarantee.epsilon_high<22) throw Refinement(id,upstream,"H-48 receiving basis needs epsilon 22");
        return X::Artifact{d,22,match,{}};
      }});
  }
  executor.add({"output","output",{{"match-lower",[=](const Demand&){return upstream;}},{"match-upper",[=](const Demand&){return upstream;}}},
    [](const Demand& d,const auto&) {return X::Artifact{d,44,"output",{}};}});
  check(executor.run("output",upstream)->value==44,"H48 both matches complete");
  for(int level=4;level<=7;++level) check(executor.builds("level-"+std::to_string(level))==1,"H48 upstream levels must not replay");
  check(executor.builds("basis-lower")==2 && executor.builds("basis-upper")==2,"H48 only deficient bases extend");
  check(executor.builds("match-lower")==2,"upper refinement preserves completed lower match");
  X bad; bad.add({"bad","bad",{},[](const Demand& d,const auto&)->X::Artifact {throw Refinement("bad",d,"no progress");}});
  rejects([&]{bad.run("bad",upstream);},"identical failed work must not loop");
  X cycle;
  for(const std::string id:{"a","b"}) cycle.add({id,id,{{id=="a"?"b":"a",[](const Demand& d){return d;}}},
    [=](const Demand& d,const auto&){return X::Artifact{d,0,id,{}};}});
  rejects([&]{cycle.run("a",upstream);},"cycles rejected before execution");
}

void correlation_test() {
  using L=LinearUncertainty; using B=L::Ball;
  auto source=std::make_shared<const L::Source>(L::Source{"common-tail",{B::from_strings("[1 +/- 0.001]"),B::from_strings("[2 +/- 0.001]")}});
  L u(source);
  Rational big("100000000000000000000000000000000000000000000000000"),one(1),zero(0);
  auto same=u.transformed({{one,big},{zero,one}}).transformed({{one,-big},{zero,one}}).independent_hull();
  check(acb_equal(same[0].raw(),source->values[0].raw()),"ill-conditioned inverse maps preserve shared tail uncertainty");
  auto canceled=u.transformed({{one,one},{one,one}}).transformed({{one,-one}}).independent_hull();
  check(canceled[0].is_zero(),"shared source cancellation is exact before enclosure");
}

void feynman_tests() {
  using namespace feynman;
  ExactField field({"x1","x2","x3","x4","eps"});
  Exact x(field,"x1");
  auto bubble=banana(1,{Rational(1),Rational("3/2")});
  auto s=symanzik(bubble,{x,x.constant(1)-x});
  check(s.U==x.constant(1),"bubble U");
  check(s.F==Exact(field,"x1+(1-x1)*3/2+x1*(1-x1)"),"bubble F with Euclidean external Gram");
  std::vector<Exact> parameters;
  for(unsigned i=0;i<4;++i)parameters.push_back(x.variable(i));
  auto weights=merge_weights(x,5,{{0,1},{0,2},{0,3},{0,4}},parameters);
  auto sum=x.constant(0);for(const auto& w:weights)sum=sum+w;
  check(sum==x.constant(1),"recursive convex weights sum to one");
  auto b4=symanzik(banana(4,std::vector<Rational>(5,Rational(1))),weights);
  auto expectedU=x.constant(0),product=x.constant(1);
  for(unsigned i=0;i<5;++i) {
    auto term=x.constant(1);
    for(unsigned j=0;j<5;++j)if(i!=j)term=term*weights[j];
    expectedU=expectedU+term;product=product*weights[i];
  }
  check(b4.U==expectedU,"banana4 first Symanzik polynomial");
  check(b4.F==expectedU+product,"banana4 second Symanzik polynomial");
  auto beta=merge_request({2,3,1},0,1);
  check(beta.operation==Operation::BetaIntegral && beta.normalization==Rational(12) && beta.source_indices==std::vector<int>{5,0,1},"beta normalization and merged powers");
  check(merge_request({0,2},0,1).operation==Operation::LowerLimit,"absent left line is lower limit");
  check(merge_request({2,0},0,1).operation==Operation::UpperLimit,"absent right line is upper limit");
  check(merge_request({0,0},0,1).operation==Operation::Direct,"both absent is direct");
  check(source_epsilon_top(4,-2,3)==9,"log sector primitive demand includes full pole depth");
  rejects([&]{merge_weights(x,3,{{0,1},{1,2}},{x,x});},"eliminated merge slot rejected");
}

int main() {try {exact_tests();recurrence_tests();h48_test();correlation_test();feynman_tests();
  std::cout<<"Exact compiler, recurrence, H48 local refinement and shared-tail tests passed\n";
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}}
