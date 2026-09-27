#pragma once
#include "diffexp/exact.hpp"
#include <flint/fmpq_poly.h>

namespace diffexp::dense_polynomial_residual {
// Exact Q[epsilon] arithmetic for checking polynomial differential identities.
// This representation is independent of the Frobenius construction recurrence.
class Polynomial {
  fmpq_poly_t value_;
 public:
  Polynomial(){fmpq_poly_init(value_);}
  Polynomial(const Polynomial& p):Polynomial(){fmpq_poly_set(value_,p.value_);}
  Polynomial(Polynomial&& p)noexcept:Polynomial(){fmpq_poly_swap(value_,p.value_);}
  Polynomial& operator=(Polynomial p)noexcept{fmpq_poly_swap(value_,p.value_);return *this;}
  ~Polynomial(){fmpq_poly_clear(value_);}
  Polynomial(const Exact& p,std::size_t variable):Polynomial(){
    fmpz_poly_t n,d;fmpz_poly_init(n);fmpz_poly_init(d);
    try{p.univariate_polynomials(n,d,variable);if(fmpz_poly_degree(d)!=0)throw std::invalid_argument("residual coefficient is not polynomial");
      fmpq_poly_set_fmpz_poly(value_,n);fmpq_poly_scalar_div_fmpz(value_,value_,d->coeffs);
    }catch(...){fmpz_poly_clear(n);fmpz_poly_clear(d);throw;}
    fmpz_poly_clear(n);fmpz_poly_clear(d);
  }
  Exact exact(const Exact& sample,std::size_t variable)const{
    fmpz_poly_t n,d;fmpz_poly_init(n);fmpz_poly_init(d);
    fmpq_poly_get_numerator(n,value_);fmpz_poly_set_fmpz(d,value_->den);
    try{auto result=sample.from_univariate_polynomials(n,d,variable);fmpz_poly_clear(n);fmpz_poly_clear(d);return result;}
    catch(...){fmpz_poly_clear(n);fmpz_poly_clear(d);throw;}
  }
  bool zero()const{return fmpq_poly_is_zero(value_);}
  std::size_t terms()const{std::size_t n=0;for(slong i=0;i<fmpq_poly_length(value_);++i)n+=!fmpz_is_zero(value_->coeffs+i);return n;}
  void set(unsigned exponent,const Rational& coefficient){fmpq_t q;fmpq_init(q);if(fmpq_set_str(q,coefficient.str().c_str(),10)){fmpq_clear(q);throw std::invalid_argument("polynomial coefficient");}fmpq_canonicalise(q);fmpq_poly_set_coeff_fmpq(value_,exponent,q);fmpq_clear(q);}
  Polynomial& operator+=(const Polynomial& p){fmpq_poly_add(value_,value_,p.value_);return *this;}
  Polynomial& operator-=(const Polynomial& p){fmpq_poly_sub(value_,value_,p.value_);return *this;}
  Polynomial scaled(long n)const{Polynomial p;fmpq_poly_scalar_mul_si(p.value_,value_,n);return p;}
  friend Polynomial operator*(const Polynomial& a,const Polynomial& b){Polynomial p;fmpq_poly_mul(p.value_,a.value_,b.value_);return p;}
};
using Bivariate=std::vector<Polynomial>; // coefficients in x, low order first
inline Bivariate split(const Exact& value,std::size_t xi,std::size_t ei,unsigned cutoff){
  if(!value.denominator().is_rational())throw std::invalid_argument("residual expected polynomial");
  const auto denominator=value.denominator().rational();Bivariate out;
  for(const auto& term:value.numerator_terms()){
    for(std::size_t v=0;v<term.powers.size();++v)if(v!=xi&&v!=ei&&term.powers[v])throw std::invalid_argument("unassigned residual polynomial variable");
    if(term.powers[xi]>cutoff)continue;
    if(out.size()<=term.powers[xi])out.resize(term.powers[xi]+1);
    out[term.powers[xi]].set(term.powers[ei],term.coefficient/denominator);
  }
  return out;
}
template<class Check> inline Bivariate multiply(const Bivariate& a,const Bivariate& b,unsigned cutoff,Check&& check){
  if(a.empty()||b.empty())return {};
  Bivariate out(std::min<std::size_t>(cutoff+1,a.size()+b.size()-1));
  for(std::size_t i=0;i<a.size()&&i<out.size();++i)if(!a[i].zero())
    for(std::size_t j=0;j<b.size()&&i+j<out.size();++j)if(!b[j].zero()){
      check(a[i].terms(),b[j].terms());out[i+j]+=a[i]*b[j];
    }
  return out;
}
} // namespace diffexp::dense_polynomial_residual
