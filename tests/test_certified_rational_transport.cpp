#include "diffexp/certified_rational_transport.hpp"
#include <iostream>
using namespace diffexp;
using B=kernel::ComplexBall;
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
template<class F>void fails(F f,const std::string& reason){
  try{f();}catch(const std::exception& e){require(std::string(e.what()).find(reason)!=std::string::npos,"unexpected certified rational failure reason");return;}
  throw std::runtime_error("uncertified rational transport was accepted");
}
int main(){try {
  B::set_precision(384);ExactField field({"x","I"});Exact x(field,"x");
  // Non-diagonal variable matrix with independent closed form:
  // y0=exp(x), y1=exp(x)*(1+log(1+x)); y1'=y1+y0/(1+x).
  std::vector<RationalLineEntry> entries{{0,0,0,Exact(field,1)},
    {1,1,0,Exact(field,1)},{1,0,0,Exact(field,"1/(1+x)")}};
  auto actual=certified_rational_line(entries,{{B(1)},{B(1)}});
  B exponential,logtwo;acb_exp(exponential.raw(),B(1).raw(),384);acb_log(logtwo.raw(),B(2).raw(),384);
  auto second=exponential*(B(1)+logtwo);
  require(acb_contains(actual.boundary[0][0].raw(),exponential.raw()),"non-diagonal exponential not contained");
  require(acb_contains(actual.boundary[1][0].raw(),second.raw()),"non-diagonal rational integral not contained");
  require(actual.conditional_on_initial_enclosure && native_enclosure_meets_tolerance(actual.boundary,actual.absolute_tolerance),"result must declare/check conditional enclosure guarantee");
  require(B::precision()==384,"rational solver must restore precision");
  require(actual.charts>0 && actual.charts<=64,"ordinary transport chart count unreasonable");

  // Large ordinary norm: a holomorphic disk can still yield a useless
  // exponential majorant. The solver must shrink it using the tail proof.
  CertifiedRationalOptions large_options;large_options.initial_radius="4";
  large_options.taylor_order=192;large_options.working_bits=384;
  auto large=certified_rational_line({{0,0,0,Exact(field,100)}},{{B(1)}},large_options);
  B exp100;acb_exp(exp100.raw(),B(100).raw(),384);
  require(acb_contains(large.boundary[0][0].raw(),exp100.raw()),"large-norm ordinary result not contained");
  require(large.disk_attempts>large.charts,"large exponential majorant did not trigger radius shrinking");

  // Fixed h/R=1/4 cannot meet the local target for N80 and this large
  // boundary. Reducing the ratio independently removes that artificial floor.
  CertifiedRationalOptions ratio_options;ratio_options.taylor_order=80;
  ratio_options.working_bits=384;ratio_options.requested_digits=28;
  auto ratio_case=certified_rational_line({{0,0,0,Exact(field,"1/(1+x)")}},{{B(1000000)}},ratio_options);
  require(acb_contains(ratio_case.boundary[0][0].raw(),B(2000000).raw()),"independent step/radius adaptation failed");

  // Epsilon coupling: exp(eps*log(2)+2 eps^2) at the endpoint.
  std::vector<RationalLineEntry> epsilon_matrix{{0,0,1,Exact(field,"1/(1+x)")},{0,0,2,Exact(field,2)}};
  Boundary seed{{B(1),B(0),B(0),B(0),B(0)}};
  auto coupled=certified_rational_line(epsilon_matrix,seed);
  Jet epsilon(0,5,384);epsilon.set(1,logtwo);epsilon.set(2,B(2));auto expected=epsilon.exp();
  for(unsigned k=0;k<5;++k)require(acb_contains(coupled.boundary[0][k].raw(),expected.at(k).raw()),"epsilon-polynomial coupled result not contained");

  // Complex contour pullback z=i*x of y'=y/(1+z), and genuinely complex
  // initial data. Exact solution is y=(1+2i)*(1+i*x).
  auto complex=certified_rational_line({{0,0,0,Exact(field,"I/(1+I*x)")}},{{B::from_strings("1","2")}});
  auto complex_expected=B::from_strings("-1","3");
  require(acb_contains(complex.boundary[0][0].raw(),complex_expected.raw()),"complex contour coefficient/I substitution not contained");
  // A rational logarithmic solution follows its analytic sheet on the leg.
  auto complex_log=certified_rational_line({{0,1,0,Exact(field,"I/(1+I*x)")}},{{B(0)},{B(1)}});
  B principal;acb_log(principal.raw(),B::from_strings("1","1").raw(),384);
  require(acb_contains(complex_log.boundary[0][0].raw(),principal.raw()),"complex rational logarithm continuation not contained");

  // Small inherited uncertainty is transported, not replaced by a midpoint.
  Boundary uncertain{{B(1)}};NativeTailMagnitude::decimal("1e-25").add_error_to(uncertain[0][0]);
  auto inherited=certified_rational_line({{0,0,0,Exact(field,1)}},uncertain);
  auto shift=B::from_strings("1/10000000000000000000000000");
  auto lower=(B(1)-shift)*exponential,upper=(B(1)+shift)*exponential;
  require(acb_contains(inherited.boundary[0][0].raw(),lower.raw()) && acb_contains(inherited.boundary[0][0].raw(),upper.raw()),"inherited uncertainty was narrowed");
  NativeTailMagnitude::decimal("0.01").add_error_to(uncertain[0][0]);
  fails([&]{certified_rational_line({{0,0,0,Exact(field,1)}},uncertain);},"final enclosure does not meet");

  CertifiedRationalOptions budget;budget.max_charts=1;
  fails([&]{certified_rational_line(entries,{{B(1)},{B(1)}},budget);},"chart budget exhausted");
  budget={};budget.max_disk_attempts=1;budget.initial_radius="4";
  fails([&]{certified_rational_line(entries,{{B(1)},{B(1)}},budget);},"disk/tail-proof budget exhausted");
  budget={};budget.working_bits=64;
  fails([&]{certified_rational_line(entries,{{B(1)},{B(1)}},budget);},"final enclosure does not meet");
  budget={};budget.max_charts=32;
  fails([&]{certified_rational_line({{0,0,0,Exact(field,"1/(1-2*x)")}},{{B(1)}},budget);},"budget exhausted");
  // A singular initial point cannot be patched over by sampling away from it.
  fails([&]{certified_rational_line({{0,0,0,Exact(field,"1/x")}},{{B(1)}});},"disk/tail-proof budget exhausted");
  ExactField extra({"x","mass"});
  fails([&]{certified_rational_line({{0,0,0,Exact(extra,"mass/(1+x)")}},{{B(1)}});},"unsubstituted parameter");
  require(B::precision()==384,"rational solver must restore precision after failure");
  std::cout<<"Certified rational systems: non-diagonal, epsilon coupling, complex contours, inherited uncertainty and proof/precision failure tests passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
