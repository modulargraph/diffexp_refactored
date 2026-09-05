#include "diffexp/univariate_rational.hpp"
#include "diffexp/cached_affine.hpp"
#include <iostream>
using namespace diffexp;
int main(){try{
  ExactField field({"x","eps","unused"});Exact x(field,"x"),e(field,"eps"),z(field,0);
  const auto a=(e*e+e.constant(3))/(e-e.constant(2)),b=e.parse("(98765432109876543210987654321*eps+7)/(eps^3+11)");
  UnivariateRational ua(a,1),ub(b,1);
  if(ua.exact(z,1)!=a||(ua+ub).exact(z,1)!=a+b||(ua-ub).exact(z,1)!=a-b||(ua*ub).exact(z,1)!=a*b||(ua/ub).exact(z,1)!=a/b)throw std::runtime_error("univariate exact arithmetic differs from Q(x,eps)");
  if(ua.constant(Rational("-7/3")).exact(z,1)!=e.parse("-7/3"))throw std::runtime_error("rational constant embedding");
  bool rejected=false;try{UnivariateRational wrong(x+e,1);}catch(const std::domain_error&){rejected=true;}
  if(!rejected)throw std::runtime_error("multivariate coefficient silently specialized");
  rejected=false;try{auto invalid=ua/UnivariateRational(0);(void)invalid;}catch(const std::domain_error&){rejected=true;}
  if(!rejected)throw std::runtime_error("univariate zero divisor accepted");
  AffineFrobeniusSeries::Matrix matrix={{e/x,e.constant(1)/x,z},{z,e/x,e.constant(1)/(e.constant(1)-x)},{z,z,(e.constant(1)-e)/x+e/(e.constant(1)+x)}};
  for(bool cleared:{false,true}){
    AffineFrobeniusSeries::Options reference;reference.finite_lag_recurrence=cleared;
    auto optimized=reference;optimized.univariate_epsilon_recurrence=true;
    auto old=AffineFrobeniusSeries::prepare(matrix,0,1,14,reference),next=AffineFrobeniusSeries::prepare(matrix,0,1,14,optimized);
    if(cached_affine::detail::payload(old)!=cached_affine::detail::payload(next))throw std::runtime_error("univariate Jordan/resonance recurrence differs from existing exact coefficients");
    auto state=cached_affine::detail::decode(cached_affine::detail::payload(next),e,3,14,optimized);
    cached_affine::detail::verify(matrix,0,1,14,optimized,state);
  }
  std::cout<<"Exact univariate conversion/arithmetic, rejected specialization, full resonant Frobenius coefficients and independent polynomial ODE checks passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
