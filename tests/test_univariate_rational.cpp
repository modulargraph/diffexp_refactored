#include "diffexp/univariate_rational.hpp"
#include "affine_series_verifier.hpp"
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
  // The Taylor and the cleared row form must give the same exact series, and
  // each must pass the independent polynomial ODE/normalization verifier.
  std::optional<AffineFrobeniusSeries> first;
  const auto same=[](const AffineFrobeniusSeries& a,const AffineFrobeniusSeries& b) {
    const auto& s=a.terms().terms;const auto& t=b.terms().terms;
    bool equal=s.size()==t.size() && a.residue_frame()==b.residue_frame() && a.jordan_successors()==b.jordan_successors() &&
      a.absolute_x_frontiers()==b.absolute_x_frontiers() && a.terms().wronskian_prefactor==b.terms().wronskian_prefactor;
    for(std::size_t i=0;equal && i<s.size();++i)equal=s[i].row==t[i].row && s[i].column==t[i].column && s[i].log_degree==t[i].log_degree &&
      s[i].power==t[i].power && s[i].slope==t[i].slope && s[i].coefficient==t[i].coefficient;
    return equal;
  };
  for(bool cleared:{false,true}){
    AffineFrobeniusSeries::Options options;options.finite_lag_recurrence=cleared;
    auto series=AffineFrobeniusSeries::prepare(matrix,0,1,14,options);
    if(first&&!same(*first,series))throw std::runtime_error("cleared and Taylor row forms give different Jordan/resonance coefficients");
    affine_series_verifier::verify(matrix,0,1,series,options);
    first=std::move(series);
  }
  std::cout<<"Exact univariate conversion/arithmetic, rejected specialization, full resonant Frobenius coefficients and independent polynomial ODE checks passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
