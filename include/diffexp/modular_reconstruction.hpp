#pragma once
#include "diffexp/exact.hpp"
#include <flint/nmod_mat.h>
#include <flint/nmod.h>
#include <flint/fmpq.h>
#include <optional>
#include <numeric>

namespace diffexp::modular {
using Word=mp_limb_t;
// FIRE7 sources/tools/primes.cpp, entries 1..16. The provider records both
// index and modulus; no platform-dependent signed arithmetic is used.
inline constexpr Word primes[]={18446744073709551557ULL,18446744073709551533ULL,
  18446744073709551521ULL,18446744073709551437ULL,18446744073709551427ULL,
  18446744073709551359ULL,18446744073709551337ULL,18446744073709551293ULL,
  18446744073709551263ULL,18446744073709551253ULL,18446744073709551191ULL,
  18446744073709551163ULL,18446744073709551113ULL,18446744073709550873ULL,
  18446744073709550791ULL,18446744073709550773ULL};
inline Word prime(unsigned index) {
  if(!index||index>std::size(primes))throw std::invalid_argument("unsupported FIRE reconstruction prime index");
  return primes[index-1];
}
struct Integer {
  fmpz_t value;
  Integer(){fmpz_init(value);} explicit Integer(Word n):Integer(){fmpz_set_ui(value,n);}
  Integer(const Integer& n):Integer(){fmpz_set(value,n.value);}
  Integer(Integer&& n)noexcept:Integer(){fmpz_swap(value,n.value);}
  Integer& operator=(Integer n){fmpz_swap(value,n.value);return *this;}
  ~Integer(){fmpz_clear(value);}
};
inline Word rational_mod(const Rational& q,Word p) {
  fmpq_t raw;fmpq_init(raw);fmpq_set_str(raw,q.str().c_str(),10);fmpq_canonicalise(raw);
  auto n=fmpz_fdiv_ui(fmpq_numref(raw),p),d=fmpz_fdiv_ui(fmpq_denref(raw),p);fmpq_clear(raw);
  if(!d)throw std::domain_error("rational coefficient has a pole modulo prime");
  nmod_t mod;nmod_init(&mod,p);return nmod_mul(n,nmod_inv(d,mod),mod);
}
inline Word evaluate(const Exact& value,const std::vector<Word>& point,Word p) {
  if(point.size()!=value.variable_count())throw std::invalid_argument("modular evaluation point arity");
  nmod_t mod;nmod_init(&mod,p);
  auto polynomial=[&](const auto& terms){Word sum=0;for(const auto& term:terms){auto c=rational_mod(term.coefficient,p);
    for(std::size_t i=0;i<point.size();++i)c=nmod_mul(c,nmod_pow_ui(point[i]%p,term.powers[i],mod),mod);
    sum=nmod_add(sum,c,mod);}return sum;};
  auto n=polynomial(value.numerator_terms()),d=polynomial(value.denominator_terms());
  if(!d)throw std::domain_error("rational function has a pole at modular point");
  return nmod_mul(n,nmod_inv(d,mod),mod);
}
using Power=std::vector<unsigned>;
inline std::vector<Power> monomials(std::size_t variables,unsigned degree,std::size_t limit=1024) {
  std::vector<Power> out;Power powers(variables);
  auto visit=[&](auto&& self,std::size_t i,unsigned left)->void {
    if(i==variables){out.push_back(powers);if(out.size()>limit)throw std::runtime_error("modular monomial budget exhausted");return;}
    for(unsigned n=0;n<=left;++n){powers[i]=n;self(self,i+1,left-n);}
  };visit(visit,0,degree);return out;
}
inline std::vector<Word> values(const std::vector<Power>& powers,const std::vector<Word>& point,nmod_t mod) {
  std::vector<Word> out;for(const auto& term:powers){Word v=1;for(std::size_t i=0;i<term.size();++i)v=nmod_mul(v,nmod_pow_ui(point.at(i)%mod.n,term[i],mod),mod);out.push_back(v);}return out;
}
struct Sample {std::vector<Word> point;std::vector<Word> coefficients;};
struct Ansatz {
  std::vector<Power> numerator,denominator;
  std::size_t pivot=0;
};
struct Image {Ansatz ansatz;std::vector<Word> coefficients;};
// Discover support at one prime, then fit only that support at later primes.
// This is a hypothesis: independent held-out primes must still validate it.
inline Image compact_support(const Image& image) {
  if(image.coefficients.size()!=image.ansatz.numerator.size()+image.ansatz.denominator.size()||image.ansatz.pivot>=image.ansatz.denominator.size())throw std::invalid_argument("invalid modular image shape");
  Image result;const auto n=image.ansatz.numerator.size();
  for(std::size_t i=0;i<n;++i)if(image.coefficients[i]){result.ansatz.numerator.push_back(image.ansatz.numerator[i]);result.coefficients.push_back(image.coefficients[i]);}
  // Keep a numerator slot for the identically zero function.
  if(result.ansatz.numerator.empty()){if(!n)throw std::invalid_argument("empty modular numerator");result.ansatz.numerator.push_back(image.ansatz.numerator[0]);result.coefficients.push_back(0);}
  for(std::size_t i=0;i<image.ansatz.denominator.size();++i)if(image.coefficients[n+i]){
    if(i==image.ansatz.pivot)result.ansatz.pivot=result.ansatz.denominator.size();
    result.ansatz.denominator.push_back(image.ansatz.denominator[i]);result.coefficients.push_back(image.coefficients[n+i]);}
  if(!image.coefficients[n+image.ansatz.pivot])throw std::invalid_argument("zero modular normalization pivot");
  return result;
}
inline std::optional<Word> evaluate(const Image& image,const std::vector<Word>& point,Word p) {
  nmod_t mod;nmod_init(&mod,p);auto n=values(image.ansatz.numerator,point,mod),d=values(image.ansatz.denominator,point,mod);
  Word numerator=0,denominator=0;
  for(std::size_t i=0;i<n.size();++i)numerator=nmod_add(numerator,nmod_mul(n[i],image.coefficients[i],mod),mod);
  for(std::size_t i=0;i<d.size();++i)denominator=nmod_add(denominator,nmod_mul(d[i],image.coefficients[n.size()+i],mod),mod);
  if(!denominator)return std::nullopt;return nmod_mul(numerator,nmod_inv(denominator,mod),mod);
}
// Homogeneous rational interpolation avoids assuming a nonzero constant
// denominator. A one-dimensional nullspace fixes the normalized coefficient
// vector. Degenerate points or an insufficient degree budget are explicit.
inline std::optional<Image> fit(const std::vector<Sample>& samples,std::size_t coefficient,
    const Ansatz& ansatz,Word p,bool fixed_pivot=false) {
  const auto columns=ansatz.numerator.size()+ansatz.denominator.size();
  if(samples.size()+1<columns||columns<2)return std::nullopt;
  nmod_mat_t matrix,kernel;nmod_mat_init(matrix,samples.size(),columns,p);nmod_mat_init(kernel,columns,columns,p);
  for(std::size_t i=0;i<samples.size();++i) {
    auto n=values(ansatz.numerator,samples[i].point,matrix->mod),d=values(ansatz.denominator,samples[i].point,matrix->mod);
    for(std::size_t j=0;j<n.size();++j)nmod_mat_entry(matrix,i,j)=n[j];
    for(std::size_t j=0;j<d.size();++j)nmod_mat_entry(matrix,i,n.size()+j)=nmod_neg(nmod_mul(samples[i].coefficients.at(coefficient),d[j],matrix->mod),matrix->mod);
  }
  auto rank=nmod_mat_nullspace(kernel,matrix);nmod_mat_clear(matrix);
  if(rank!=1){nmod_mat_clear(kernel);return std::nullopt;}
  auto pivot=ansatz.pivot;
  if(!fixed_pivot){pivot=0;while(pivot<ansatz.denominator.size()&&!nmod_mat_entry(kernel,ansatz.numerator.size()+pivot,0))++pivot;}
  if(pivot>=ansatz.denominator.size()||!nmod_mat_entry(kernel,ansatz.numerator.size()+pivot,0)){nmod_mat_clear(kernel);return std::nullopt;}
  auto scale=nmod_inv(nmod_mat_entry(kernel,ansatz.numerator.size()+pivot,0),kernel->mod);
  Image result{ansatz,{}};result.ansatz.pivot=pivot;
  for(std::size_t i=0;i<columns;++i)result.coefficients.push_back(nmod_mul(nmod_mat_entry(kernel,i,0),scale,kernel->mod));
  nmod_mat_clear(kernel);return result;
}
inline std::optional<Image> discover(const std::vector<Sample>& samples,std::size_t coefficient,
    std::size_t variables,unsigned max_degree,Word p) {
  if(samples.size()<4)return std::nullopt;
  const std::vector<Sample> training(samples.begin(),samples.end()-3);
  std::vector<std::vector<Power>> terms;for(unsigned d=0;d<=max_degree;++d){
    try{terms.push_back(monomials(variables,d));}catch(const std::runtime_error&){break;}}
  struct Choice{std::size_t count;unsigned n,d;};std::vector<Choice> choices;
  for(unsigned n=0;n<terms.size();++n)for(unsigned d=0;d<terms.size();++d)
    if(terms[n].size()+terms[d].size()<=training.size()+1)choices.push_back({terms[n].size()+terms[d].size(),n,d});
  std::sort(choices.begin(),choices.end(),[](auto a,auto b){return std::tie(a.count,a.n,a.d)<std::tie(b.count,b.n,b.d);});
  for(auto c:choices)if(auto result=fit(training,coefficient,{terms[c.n],terms[c.d],0},p)) {
    bool valid=true;for(std::size_t i=training.size();i<samples.size();++i)
      if(evaluate(*result,samples[i].point,p)!=std::optional<Word>(samples[i].coefficients.at(coefficient)))valid=false;
    if(valid)return result;
  }
  return std::nullopt;
}
struct Lift {
  Ansatz ansatz;std::vector<Integer> residues;Integer modulus;
  explicit Lift(const Image& image,Word p):ansatz(image.ansatz),modulus(p){for(auto c:image.coefficients)residues.emplace_back(c);}
  void append(const Image& image,Word p) {
    if(image.ansatz.numerator!=ansatz.numerator||image.ansatz.denominator!=ansatz.denominator||image.ansatz.pivot!=ansatz.pivot||image.coefficients.size()!=residues.size()||fmpz_fdiv_ui(modulus.value,p)==0)
      throw std::invalid_argument("incompatible or repeated modular image");
    for(std::size_t i=0;i<residues.size();++i)fmpz_CRT_ui(residues[i].value,residues[i].value,modulus.value,image.coefficients[i],p,0);
    fmpz_mul_ui(modulus.value,modulus.value,p);
  }
  std::optional<Exact> reconstruct(const Exact& sample,const std::vector<std::size_t>& active)const {
    auto n=sample.constant(0),d=sample.constant(0);
    for(std::size_t i=0;i<residues.size();++i) {
      fmpq_t q;fmpq_init(q);if(!fmpq_reconstruct_fmpz(q,residues[i].value,modulus.value)){fmpq_clear(q);return std::nullopt;}
      char* text=fmpq_get_str(nullptr,10,q);auto value=sample.parse(text);flint_free(text);fmpq_clear(q);
      auto& powers=i<ansatz.numerator.size()?ansatz.numerator[i]:ansatz.denominator[i-ansatz.numerator.size()];
      if(powers.size()!=active.size())throw std::invalid_argument("reconstruction active variable arity");
      for(std::size_t j=0;j<active.size();++j)if(powers[j])value=value*sample.variable(active[j]).pow(powers[j]);
      if(i<ansatz.numerator.size())n=n+value;else d=d+value;
    }
    if(d.is_zero())return std::nullopt;return n/d;
  }
};
inline std::uint64_t mix(std::uint64_t x) {
  x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);
}
inline std::vector<Word> point(std::size_t variables,unsigned prime_index,std::size_t ordinal) {
  std::vector<Word> out;for(std::size_t j=0;j<variables;++j)out.push_back(17+mix(mix(ordinal)^mix(j+123)^mix(prime_index+987))%1000003);return out;
}
} // namespace diffexp::modular
