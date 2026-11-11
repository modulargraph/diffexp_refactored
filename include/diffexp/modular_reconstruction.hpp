#pragma once
#include "diffexp/exact.hpp"
#include <flint/nmod_mat.h>
#include <flint/nmod.h>
#include <flint/fmpq.h>
#include <optional>
#include <functional>
#include <numeric>
#include <span>

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
  auto n=fmpz_fdiv_ui(fmpq_numref(q.raw()),p),d=fmpz_fdiv_ui(fmpq_denref(q.raw()),p);
  if(!d)throw std::domain_error("rational coefficient has a pole modulo prime");
  nmod_t mod;nmod_init(&mod,p);return nmod_mul(n,nmod_inv(d,mod),mod);
}
inline Word evaluate(const Exact& value,const std::vector<Word>& point,Word p) {
  nmod_t mod;nmod_init(&mod,p);auto image=value.modular_value(point,mod);
  if(!image)throw std::domain_error("rational function has a pole at modular point");
  return *image;
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
// Fits read samples in place: one row carries every coefficient of a probe, so
// copying a training prefix per coefficient would cost O(samples*coefficients).
using Samples=std::span<const Sample>;
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
inline std::optional<Image> fit(Samples samples,std::size_t coefficient,
    const Ansatz& ansatz,Word p,bool fixed_pivot=false) {
  const auto columns=ansatz.numerator.size()+ansatz.denominator.size();
  if(samples.size()+1<columns||columns<2)return std::nullopt;
  nmod_mat_t matrix,kernel;nmod_mat_init(matrix,samples.size(),columns,p);nmod_mat_init(kernel,columns,columns,p);
  // The powers of each coordinate are computed once per sample and shared by
  // all monomials of the ansatz.
  std::vector<unsigned> degree;
  for(const auto* terms:{&ansatz.numerator,&ansatz.denominator})for(const auto& term:*terms){
    if(term.size()>degree.size())degree.resize(term.size());for(std::size_t v=0;v<term.size();++v)degree[v]=std::max(degree[v],term[v]);}
  std::vector<std::vector<Word>> powers(degree.size());for(std::size_t v=0;v<degree.size();++v)powers[v].resize(degree[v]+1);
  const auto monomial=[&](const Power& term){Word value=1;for(std::size_t v=0;v<term.size();++v)if(term[v])value=nmod_mul(value,powers[v][term[v]],matrix->mod);return value;};
  const auto n=ansatz.numerator.size();
  for(std::size_t i=0;i<samples.size();++i) {
    for(std::size_t v=0;v<degree.size();++v){const Word x=samples[i].point.at(v)%p;powers[v][0]=1;for(unsigned e=1;e<=degree[v];++e)powers[v][e]=nmod_mul(powers[v][e-1],x,matrix->mod);}
    const Word negative=nmod_neg(samples[i].coefficients.at(coefficient)%p,matrix->mod);
    for(std::size_t j=0;j<n;++j)nmod_mat_entry(matrix,i,j)=monomial(ansatz.numerator[j]);
    for(std::size_t j=0;j<ansatz.denominator.size();++j)nmod_mat_entry(matrix,i,n+j)=nmod_mul(negative,monomial(ansatz.denominator[j]),matrix->mod);
  }
  // The kernel of the leading rows decides without eliminating all of them.
  // If it is one vector, the kernel of the complete system is that vector
  // when the remaining rows annihilate it and empty otherwise; if it is
  // empty, so is the complete kernel. Only a larger kernel of the leading
  // rows needs the complete system. The result is the same in every case.
  slong rank=-1;const slong rows=nmod_mat_nrows(matrix),lead=static_cast<slong>(columns)+2;
  if(rows>lead+2) {
    nmod_mat_t leading;nmod_mat_window_init(leading,matrix,0,0,lead,columns);rank=nmod_mat_nullspace(kernel,leading);nmod_mat_window_clear(leading);
    if(rank==1)for(slong i=lead;i<rows&&rank==1;++i){Word dot=0;
      for(std::size_t j=0;j<columns;++j)dot=nmod_add(dot,nmod_mul(nmod_mat_entry(matrix,i,j),nmod_mat_entry(kernel,j,0),matrix->mod),matrix->mod);
      if(dot)rank=0;}
    else if(rank>1)rank=-1;
  }
  if(rank<0)rank=nmod_mat_nullspace(kernel,matrix);
  nmod_mat_clear(matrix);
  if(rank!=1){nmod_mat_clear(kernel);return std::nullopt;}
  auto pivot=ansatz.pivot;
  if(!fixed_pivot){pivot=0;while(pivot<ansatz.denominator.size()&&!nmod_mat_entry(kernel,ansatz.numerator.size()+pivot,0))++pivot;}
  if(pivot>=ansatz.denominator.size()||!nmod_mat_entry(kernel,ansatz.numerator.size()+pivot,0)){nmod_mat_clear(kernel);return std::nullopt;}
  auto scale=nmod_inv(nmod_mat_entry(kernel,ansatz.numerator.size()+pivot,0),kernel->mod);
  Image result{ansatz,{}};result.ansatz.pivot=pivot;
  for(std::size_t i=0;i<columns;++i)result.coefficients.push_back(nmod_mul(nmod_mat_entry(kernel,i,0),scale,kernel->mod));
  nmod_mat_clear(kernel);return result;
}
// The last `holdouts` samples never enter a fit. A candidate is accepted only if
// it reproduces them: with rows = columns-1 a kernel vector exists for ANY data,
// so agreement with the training rows alone proves nothing.
inline constexpr std::size_t holdouts=3;
inline std::optional<Image> fit_with_holdouts(Samples samples,std::size_t coefficient,
    const Ansatz& ansatz,Word p) {
  if(samples.size()<=holdouts)return std::nullopt;
  auto image=fit(samples.first(samples.size()-holdouts),coefficient,ansatz,p);if(!image)return std::nullopt;
  for(const auto& held:samples.last(holdouts))
    if(evaluate(*image,held.point,p)!=std::optional<Word>(held.coefficients.at(coefficient)))return std::nullopt;
  return image;
}
// `checkpoint`, when supplied, is called before every dense ansatz so a caller's
// deadline is honoured inside the (n,d) scan, not only between coefficients.
inline std::optional<Image> discover(Samples samples,std::size_t coefficient,
    std::size_t variables,unsigned max_degree,Word p,const std::function<void()>& checkpoint={}) {
  if(samples.size()<=holdouts)return std::nullopt;
  const auto training=samples.size()-holdouts;
  std::vector<std::vector<Power>> terms;for(unsigned d=0;d<=max_degree;++d){
    try{terms.push_back(monomials(variables,d));}catch(const std::runtime_error&){break;}}
  struct Choice{std::size_t count;unsigned n,d;};std::vector<Choice> choices;
  for(unsigned n=0;n<terms.size();++n)for(unsigned d=0;d<terms.size();++d)
    if(terms[n].size()+terms[d].size()<=training+1)choices.push_back({terms[n].size()+terms[d].size(),n,d});
  std::sort(choices.begin(),choices.end(),[](auto a,auto b){return std::tie(a.count,a.n,a.d)<std::tie(b.count,b.n,b.d);});
  for(auto c:choices){if(checkpoint)checkpoint();if(auto result=fit_with_holdouts(samples,coefficient,{terms[c.n],terms[c.d],0},p))return result;}
  return std::nullopt;
}
// True once every (numerator,denominator) total-degree pair up to max_degree
// already fits the training rows: further samples cannot help the dense scan.
inline bool dense_scan_saturated(std::size_t samples,std::size_t variables,unsigned max_degree) {
  if(samples<=holdouts)return false;
  try{return 2*monomials(variables,max_degree).size()<=samples-holdouts+1;}
  catch(const std::runtime_error&){return false;}
}
// A generic univariate slice preserves each coordinate degree and valuation
// unless that specialization cancels a factor. Treat the resulting box only
// as a hypothesis: fit it at unrelated multivariate points and retain heldouts.
// The box is NOT cut by a total degree: a coefficient of bidegree (15,2) is
// far cheaper than its dense total-degree-17 hull, and a global total-degree
// cap would silently delete monomials such as x^15 eps^2 from the ansatz.
inline std::optional<Ansatz> degree_box(const std::vector<Image>& slices,unsigned max_variable_degree,
    std::size_t limit=1024) {
  Ansatz result;
  auto terms=[&](bool denominator)->std::optional<std::vector<Power>> {
    std::vector<std::pair<unsigned,unsigned>> ranges;
    for(const auto& image:slices){
      const auto& powers=denominator?image.ansatz.denominator:image.ansatz.numerator;
      const auto offset=denominator?image.ansatz.numerator.size():0;
      if(image.coefficients.size()!=image.ansatz.numerator.size()+image.ansatz.denominator.size())throw std::invalid_argument("slice image shape");
      unsigned lo=max_variable_degree+1,hi=0;
      for(std::size_t j=0;j<powers.size();++j){if(powers[j].size()!=1)throw std::invalid_argument("degree slice must be univariate");
        if(image.coefficients[offset+j]){lo=std::min(lo,powers[j][0]);hi=std::max(hi,powers[j][0]);}}
      if(lo>max_variable_degree)return std::nullopt;ranges.emplace_back(lo,hi);
    }
    std::vector<Power> out;Power power(slices.size());bool overflow=false;
    auto visit=[&](auto&& self,std::size_t i)->void {
      if(overflow)return;if(i==ranges.size()){out.push_back(power);overflow=out.size()>limit;return;}
      for(unsigned e=ranges[i].first;e<=ranges[i].second&&!overflow;++e){power[i]=e;self(self,i+1);}
    };visit(visit,0);if(overflow||out.empty())return std::nullopt;return out;
  };
  auto n=terms(false),d=terms(true);if(!n||!d)return std::nullopt;
  result.numerator=std::move(*n);result.denominator=std::move(*d);return result;
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
  // The rational coefficients are cleared of denominators together, and the
  // integer numerator and denominator are divided by their gcd once.
  std::optional<Exact> reconstruct(const Exact& sample,const std::vector<std::size_t>& active)const {
    struct Rationals {std::vector<fmpq> q;fmpz_t scale;explicit Rationals(std::size_t n):q(n){for(auto& v:q)fmpq_init(&v);fmpz_init_set_ui(scale,1);}
      ~Rationals(){for(auto& v:q)fmpq_clear(&v);fmpz_clear(scale);}Rationals(const Rationals&)=delete;}values(residues.size());
    for(std::size_t i=0;i<residues.size();++i) {
      if(!fmpq_reconstruct_fmpz(&values.q[i],residues[i].value,modulus.value))return std::nullopt;
      fmpz_lcm(values.scale,values.scale,fmpq_denref(&values.q[i]));
    }
    std::vector<std::vector<unsigned long>> powers(residues.size(),std::vector<unsigned long>(sample.variable_count()));
    std::vector<Exact::IntegerTerm> n,d;
    for(std::size_t i=0;i<residues.size();++i) {
      const auto& term=i<ansatz.numerator.size()?ansatz.numerator[i]:ansatz.denominator[i-ansatz.numerator.size()];
      if(term.size()!=active.size())throw std::invalid_argument("reconstruction active variable arity");
      for(std::size_t j=0;j<active.size();++j)powers[i].at(active[j])+=term[j];
      // q = a/b with b | scale: the integer coefficient is a*(scale/b).
      auto* a=fmpq_numref(&values.q[i]);fmpz_mul(a,a,values.scale);fmpz_divexact(a,a,fmpq_denref(&values.q[i]));
      if(!fmpz_is_zero(a))(i<ansatz.numerator.size()?n:d).push_back({a,&powers[i]});
    }
    if(d.empty())return std::nullopt;
    try{return sample.quotient(n,d);}catch(const std::domain_error&){return std::nullopt;}
  }
};
inline std::uint64_t mix(std::uint64_t x) {
  x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);
}
inline std::vector<Word> point(std::size_t variables,unsigned prime_index,std::size_t ordinal) {
  std::vector<Word> out;for(std::size_t j=0;j<variables;++j)out.push_back(17+mix(mix(ordinal)^mix(j+123)^mix(prime_index+987))%1000003);return out;
}
} // namespace diffexp::modular
