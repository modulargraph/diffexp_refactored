#pragma once
#include "diffexp/exact.hpp"
#include <flint/fmpz_poly_q.h>

namespace diffexp {
// A coefficient domain for exact recurrences after the path variable has been
// expanded. Dense univariate FLINT arithmetic avoids repeated multivariate
// deflation/inflation and uses the polynomial multiplication algorithms for Q(t).
class UnivariateRational {
 public:
  UnivariateRational(){fmpz_poly_q_init(value_);}
  explicit UnivariateRational(long n):UnivariateRational(){fmpz_poly_q_set_si(value_,n);}
  UnivariateRational(const Exact& value,std::size_t variable):UnivariateRational(){
    value.univariate_polynomials(value_->num,value_->den,variable);
  }
  UnivariateRational(const UnivariateRational& other):UnivariateRational(){fmpz_poly_q_set(value_,other.value_);}
  UnivariateRational(UnivariateRational&& other)noexcept:UnivariateRational(){fmpz_poly_q_swap(value_,other.value_);}
  UnivariateRational& operator=(UnivariateRational other)noexcept{fmpz_poly_q_swap(value_,other.value_);return *this;}
  ~UnivariateRational(){fmpz_poly_q_clear(value_);}
  bool is_zero()const{return fmpz_poly_q_is_zero(value_);}
  UnivariateRational constant(long n)const{return UnivariateRational(n);}
  UnivariateRational constant(const Rational& value)const {
    fmpq_t q;fmpq_init(q);fmpq_set_str(q,value.str().c_str(),10);fmpq_canonicalise(q);
    UnivariateRational out;fmpz_poly_set_fmpz(out.value_->num,fmpq_numref(q));fmpz_poly_set_fmpz(out.value_->den,fmpq_denref(q));fmpq_clear(q);return out;
  }
  Exact exact(const Exact& sample,std::size_t variable)const{return sample.from_univariate_polynomials(value_->num,value_->den,variable);}
  friend UnivariateRational operator+(const UnivariateRational& a,const UnivariateRational& b){UnivariateRational out;fmpz_poly_q_add(out.value_,a.value_,b.value_);return out;}
  friend UnivariateRational operator-(const UnivariateRational& a,const UnivariateRational& b){UnivariateRational out;fmpz_poly_q_sub(out.value_,a.value_,b.value_);return out;}
  friend UnivariateRational operator-(const UnivariateRational& a){UnivariateRational out;fmpz_poly_q_neg(out.value_,a.value_);return out;}
  friend UnivariateRational operator*(const UnivariateRational& a,const UnivariateRational& b){UnivariateRational out;fmpz_poly_q_mul(out.value_,a.value_,b.value_);return out;}
  friend UnivariateRational operator/(const UnivariateRational& a,const UnivariateRational& b){if(b.is_zero())throw std::domain_error("zero univariate divisor");UnivariateRational out;fmpz_poly_q_div(out.value_,a.value_,b.value_);return out;}
  friend bool operator==(const UnivariateRational& a,const UnivariateRational& b){return fmpz_poly_q_equal(a.value_,b.value_);}
 private:
  fmpz_poly_q_t value_;
};
} // namespace diffexp
