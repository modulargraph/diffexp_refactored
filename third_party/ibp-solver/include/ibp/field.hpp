#pragma once
#include <flint/nmod.h>
#include <flint/ulong_extras.h>
#include <flint/fmpq.h>
#include <cstdint>
#include <string>
#include <stdexcept>
namespace ibp {
using Word=std::uint64_t;
class Field {
 Word p_; nmod_t modulus_;
 public:
 explicit Field(Word p=2305843009213693951ULL):p_(p) {
  if(p<3 || p>2305843009213693951ULL || !n_is_prime(p))throw std::invalid_argument("prime must be odd, prime, and at most 2^61-1");nmod_init(&modulus_,p);
 }
 Word prime()const{return p_;}
 Word add(Word a,Word b)const {auto s=a+b;return s>=p_?s-p_:s;}
 Word sub(Word a,Word b)const {return a>=b?a-b:p_-(b-a);}
 Word mul(Word a,Word b)const {
  if(p_==2305843009213693951ULL){__uint128_t v=static_cast<__uint128_t>(a)*b;Word s=static_cast<Word>(v)&p_;s+=static_cast<Word>(v>>61);return s>=p_?s-p_:s;}
  if(p_==2147483647ULL){Word v=a*b,s=(v&p_)+(v>>31);return s>=p_?s-p_:s;}
  return nmod_mul(a,b,modulus_);
 }
 Word inv(Word a)const {if(!a)throw std::domain_error("zero pivot at this prime/sample");return n_invmod(a,p_);}
 Word integer(std::int64_t n)const {if(n>=0)return static_cast<Word>(n)%p_;return sub(0,(static_cast<Word>(-(n+1))+1)%p_);}
 Word rational(const std::string& text)const {
  fmpq_t q;fmpq_init(q);if(fmpq_set_str(q,text.c_str(),10)){fmpq_clear(q);throw std::invalid_argument("expected exact rational: "+text);}if(fmpz_is_zero(fmpq_denref(q))){fmpq_clear(q);throw std::invalid_argument("zero rational denominator");}fmpq_canonicalise(q);
  auto a=fmpz_fdiv_ui(fmpq_numref(q),p_),b=fmpz_fdiv_ui(fmpq_denref(q),p_);fmpq_clear(q);return mul(a,inv(b));
 }
};
}
