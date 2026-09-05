#include "diffexp/laurent_transport.hpp"
#include <iostream>
using namespace diffexp;
double size(const Jet::Ball& x){arf_t u;arf_init(u);acb_get_abs_ubound_arf(u,x.raw(),384);auto d=arf_get_d(u,ARF_RND_CEIL);arf_clear(u);return d;}
int main(){try {
  using B=Jet::Ball;B::set_precision(384);ExactField field({"eps","I","x"});Exact x(field,"x"),eps(field,"eps"),ii(field,"I");
  LaurentBoundary initial{0,{{B(1),B(0),B(0),B(0),B(0)}},true};
  std::vector<Exact> contour{x.constant(0),ii/x.constant(10),x.constant(1)+ii/x.constant(10),x.constant(1)};
  auto result=integrate_regular_rows({{eps}},{{x.constant(1)/eps.pow(2)}},initial,contour,2,60);
  if(result.integrals.low!=-2 || result.integrals.high()!=2)throw std::runtime_error("regular integral shifted Laurent window");
  Rational factorial(1);
  for(int k=-2;k<=2;++k){factorial*=Rational(k+3);if(size(result.integrals.values[0][k+2]-B::from_strings((Rational(1)/factorial).str()))>1e-60)
    throw std::runtime_error("complex-contour epsilon-pole integral disagrees with (exp(eps)-1)/eps^3");}
  auto scalar=apply_rational_rows({{x.constant(1)/(eps*(x.constant(1)-eps))}},x.constant(1),result.endpoint,2);
  if(scalar.low!=-1)throw std::runtime_error("rational projection lost epsilon valuation");
  bool needed=false;try{apply_rational_rows({{x.constant(1)/eps.pow(3)}},x.constant(1),result.endpoint,2);}
  catch(const BoundaryDemand& demand){needed=demand.required_high==5;}
  if(!needed)throw std::runtime_error("missing lookahead silently padded instead of requested");
  auto rational=integrate_regular_rows({{-x.constant(1)/(x.constant(2)-x)}},{{x.constant(1)/(eps*(x.constant(2)-x))}},
    {0,{{B(2),B(0),B(0)}},true},contour,1,60);
  if(size(rational.integrals.values[0][0]-B(1))>1e-60 || size(rational.endpoint.values[0][0]-B(1))>1e-60)
    throw std::runtime_error("rational contour pullback/accumulator Jacobian");
  for(unsigned k=1;k<rational.integrals.values[0].size();++k)if(size(rational.integrals.values[0][k])>1e-60)
    throw std::runtime_error("rational zero coefficients");
  std::cout<<"Laurent contour transport and explicit boundary demands passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
