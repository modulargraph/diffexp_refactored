#include "diffexp/laurent_transport.hpp"
#include <iostream>
using namespace diffexp;
double size(const Jet::Ball& x){arf_t u;arf_init(u);acb_get_abs_ubound_arf(u,x.raw(),384);auto d=arf_get_d(u,ARF_RND_CEIL);arf_clear(u);return d;}
int main(){try {
  using B=Jet::Ball;B::set_precision(384);ExactField field({"eps","I","x"});Exact x(field,"x"),eps(field,"eps"),ii(field,"I");
  // A Laurent window is transported coefficient by coefficient along the
  // complex contour: Y'=eps*Y from Y(0)=1 gives exp(eps) at x=1.
  LaurentBoundary initial{0,{{B(1),B(0),B(0),B(0),B(0)}},true};
  std::vector<Exact> contour{x.constant(0),ii/x.constant(10),x.constant(1)+ii/x.constant(10),x.constant(1)};
  auto transported=transport_laurent({{eps}},initial,contour,60);
  Rational factorial(1);
  for(int k=0;k<=4;++k){if(k)factorial*=Rational(k);if(size(transported.values[0][k]-B::from_strings((Rational(1)/factorial).str()))>1e-60)
    throw std::runtime_error("complex-contour Laurent transport disagrees with exp(eps)");}
  if(transported.taylor_tail_certified)throw std::runtime_error("transported Laurent window kept a tail certificate");
  // Entire but stiff: no pole bounds the geometric step, so the retained
  // Taylor tail alone must shorten it. One order-80 chart of unit length
  // would be wrong at 8e-9 (rate 40) and 6e-3 (rate 60).
  double tail=0;(void)rational_chart(std::vector<RationalLineEntry>{{0,0,0,Exact(field,"40")}},{{B(1)}},B(0),B(1),80,&tail);
  if(taylor_tail_consistent(tail,80,384))throw std::runtime_error("unconverged unit chart reported a consistent Taylor tail");
  (void)rational_chart(std::vector<RationalLineEntry>{{0,0,0,Exact(field,"2")}},{{B(1)}},B(0),B(1),80,&tail);
  if(!taylor_tail_consistent(tail,80,384))throw std::runtime_error("converged chart reported an inconsistent Taylor tail");
  for(const char* rate:{"40","60"}) {
    auto stiff=rational_line({{0,0,0,Exact(field,rate)}},{{B(1)}},80);
    B expected;acb_exp(expected.raw(),B::from_strings(rate).raw(),384);
    if(size(stiff[0][0]/expected-B(1))>1e-30)throw std::runtime_error("stiff pole-free line was not subdivided to its Taylor accuracy");
  }
  // Two simple roots 2^-20 apart are not separable at the requested 16 bits;
  // an exact squarefree polynomial must be retried at higher precision.
  auto near_roots=polynomial_roots((x-x.constant(1))*(x.constant(1048576)*x-x.constant(1048577)),2,16);
  if(near_roots.size()!=2 || size(near_roots[0]-near_roots[1])>1e-5 || near_roots[0].contains_zero())
    throw std::runtime_error("close simple roots were not isolated after a precision retry");
  std::cout<<"Laurent contour transport, Taylor-tail step control and root isolation passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
