#pragma once
#include "diffexp/exact.hpp"
#include <flint/nmod_mat.h>
#include <flint/nmod.h>
#include <flint/fmpq.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>
#include <flint/nmod_mpoly.h>
#include <flint/nmod_mpoly_factor.h>
#include <algorithm>
#include <map>
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
// The denominators of the coefficient functions of a reduction are products
// of few irreducible polynomials with small rational coefficients. Once
// those are known over the rationals, the denominator of every function is
// known modulo any further prime, and a fit at such a prime has the
// numerator alone to find: as many points as the numerator has terms, not as
// numerator and denominator together.
//
// The denominators of one prime are factored modulo that prime. A factor is
// identified, at every prime, by the functions whose denominators it divides
// and how often: factors that divide the same denominators equally often
// are multiplied together (the pieces a polynomial irreducible over the
// rationals splits into modulo a prime are such, and so are factors that
// only occur together). The monic products of several primes are combined
// by the Chinese remainder theorem until their coefficients are rational
// numbers well within the modulus.
//
// Everything here is a hypothesis: the fits made with it carry spare points
// of their own (fit_numerator), and a complete fit can be compared with what
// it predicts.
class DenominatorFactors {
  struct Context {nmod_mpoly_ctx_t value;Context(std::size_t variables,Word p){nmod_mpoly_ctx_init(value,static_cast<slong>(variables),ORD_LEX,p);}~Context(){nmod_mpoly_ctx_clear(value);}Context(const Context&)=delete;};
  struct Polynomial {nmod_mpoly_t value;const Context* context;explicit Polynomial(const Context& c):context(&c){nmod_mpoly_init(value,c.value);}
    Polynomial(Polynomial&& other)noexcept:context(other.context){nmod_mpoly_init(value,context->value);nmod_mpoly_swap(value,other.value,context->value);}
    ~Polynomial(){nmod_mpoly_clear(value,context->value);}};
  // Where a factor occurs: (function, multiplicity), ascending.
  using Occurrence=std::vector<std::pair<std::uint32_t,std::uint32_t>>;
  struct Term {Power power;Word residue=0;};
  // The monic products at one prime, by occurrence.
  static std::optional<std::map<Occurrence,std::vector<Term>>> products(const std::vector<Image>& images,std::size_t variables,Word p,const std::function<void()>& checkpoint) {
    Context context(variables,p);std::vector<Polynomial> pool;std::vector<Occurrence> occurs;Polynomial denominator(context),quotient(context);std::vector<ulong> exponent(variables);
    for(std::size_t i=0;i<images.size();++i){const auto& image=images[i];
      if(checkpoint)checkpoint();
      nmod_mpoly_zero(denominator.value,context.value);const auto n=image.ansatz.numerator.size();
      for(std::size_t j=0;j<image.ansatz.denominator.size();++j){const auto& term=image.ansatz.denominator[j];if(term.size()!=variables)return std::nullopt;
        for(std::size_t v=0;v<variables;++v)exponent[v]=term[v];
        if(image.coefficients.at(n+j))nmod_mpoly_push_term_ui_ui(denominator.value,image.coefficients[n+j],exponent.data(),context.value);}
      nmod_mpoly_sort_terms(denominator.value,context.value);nmod_mpoly_combine_like_terms(denominator.value,context.value);
      if(nmod_mpoly_is_zero(denominator.value,context.value))return std::nullopt;
      // Most denominators are made of factors already found.
      for(std::size_t f=0;f<pool.size();++f){std::uint32_t multiplicity=0;
        while(nmod_mpoly_divides(quotient.value,denominator.value,pool[f].value,context.value)){nmod_mpoly_swap(quotient.value,denominator.value,context.value);++multiplicity;}
        if(multiplicity)occurs[f].emplace_back(i,multiplicity);}
      if(nmod_mpoly_is_ui(denominator.value,context.value))continue;
      nmod_mpoly_factor_t factors;nmod_mpoly_factor_init(factors,context.value);
      const bool factored=nmod_mpoly_factor(factors,denominator.value,context.value);
      if(factored)for(slong f=0;f<factors->num;++f){
        Polynomial monic(context);nmod_mpoly_make_monic(monic.value,factors->poly+f,context.value);
        occurs.push_back({{static_cast<std::uint32_t>(i),static_cast<std::uint32_t>(fmpz_get_ui(factors->exp+f))}});pool.push_back(std::move(monic));}
      nmod_mpoly_factor_clear(factors,context.value);
      if(!factored)return std::nullopt;
    }
    std::map<Occurrence,Polynomial> merged;
    for(std::size_t f=0;f<pool.size();++f){auto [it,inserted]=merged.try_emplace(occurs[f],context);
      if(inserted)nmod_mpoly_set(it->second.value,pool[f].value,context.value);else nmod_mpoly_mul(it->second.value,it->second.value,pool[f].value,context.value);}
    std::map<Occurrence,std::vector<Term>> out;
    for(const auto& [where,product]:merged){auto& terms=out[where];
      for(slong t=0;t<nmod_mpoly_length(product.value,context.value);++t){nmod_mpoly_get_term_exp_ui(exponent.data(),product.value,t,context.value);
        auto& term=terms.emplace_back();for(auto e:exponent)term.power.push_back(static_cast<unsigned>(e));term.residue=nmod_mpoly_get_term_coeff_ui(product.value,t,context.value);}}
    return out;
  }
 public:
  // From the denominators of `images` at prime p (their terms are powers of
  // `variables` variables). Nothing if a denominator is not factored.
  static std::optional<DenominatorFactors> find(const std::vector<Image>& images,std::size_t variables,Word p,const std::function<void()>& checkpoint={}) {
    if(!variables||images.empty())return std::nullopt;
    auto first=products(images,variables,p,checkpoint);if(!first)return std::nullopt;
    DenominatorFactors out;out.variables_=variables;out.powers_.resize(images.size());fmpz_set_ui(out.modulus_.value,p);out.primes_=1;
    for(auto& [where,terms]:*first){
      for(const auto& [function,multiplicity]:where)out.powers_[function].emplace_back(out.factors_.size(),multiplicity);
      auto& factor=out.factors_.emplace_back();factor.where=where;
      for(const auto& term:terms){factor.powers.push_back(term.power);factor.residues.emplace_back(term.residue);}}
    out.lift();return out;
  }
  // The same denominators at a further prime. False if their factors do not
  // occur as those of the first prime do, or have other terms.
  bool absorb(const std::vector<Image>& images,Word p,const std::function<void()>& checkpoint={}) {
    if(fmpz_fdiv_ui(modulus_.value,p)==0)return false;
    const auto next=products(images,variables_,p,checkpoint);if(!next||next->size()!=factors_.size())return false;
    for(const auto& factor:factors_){const auto it=next->find(factor.where);if(it==next->end()||it->second.size()!=factor.powers.size())return false;
      for(std::size_t t=0;t<factor.powers.size();++t)if(it->second[t].power!=factor.powers[t])return false;}
    for(auto& factor:factors_){const auto& terms=next->at(factor.where);
      for(std::size_t t=0;t<terms.size();++t)fmpz_CRT_ui(factor.residues[t].value,factor.residues[t].value,modulus_.value,terms[t].residue,p,0);}
    fmpz_mul_ui(modulus_.value,modulus_.value,p);++primes_;lift();return true;
  }
  // Whether every coefficient of every factor is a rational number with
  // numerator and denominator nine bits below the square root of half the
  // modulus: of all residues one in 2^18 is such a number by chance.
  bool lifted()const{return lifted_;}
  std::size_t factors()const{return factors_.size();}
  unsigned primes()const{return primes_;}
  // The denominators at prime p, each with the coefficients of its ansatz in
  // order and one at the pivot. Nothing if a factor has a pole modulo p, or
  // if the product does not have exactly the terms of the ansatz.
  class At {
   public:
    At(const DenominatorFactors& factors,Word p):owner_(factors),context_(factors.variables_,p),product_(context_),power_(context_){
      if(!factors.lifted_){valid_=false;return;}
      nmod_t mod;nmod_init(&mod,p);std::vector<ulong> exponent(factors.variables_);
      for(const auto& source:factors.factors_){auto& factor=factors_.emplace_back(context_);
        for(std::size_t t=0;t<source.powers.size();++t){const auto d=fmpz_fdiv_ui(source.denominators[t].value,p);if(!d){valid_=false;return;}
          for(std::size_t v=0;v<exponent.size();++v)exponent[v]=source.powers[t][v];
          if(const auto c=nmod_mul(fmpz_fdiv_ui(source.numerators[t].value,p),nmod_inv(d,mod),mod))nmod_mpoly_push_term_ui_ui(factor.value,c,exponent.data(),context_.value);}
        nmod_mpoly_sort_terms(factor.value,context_.value);nmod_mpoly_combine_like_terms(factor.value,context_.value);}
    }
    std::optional<std::vector<Word>> denominator(std::size_t index,const Ansatz& ansatz) {
      if(!valid_)return std::nullopt;
      nmod_mpoly_one(product_.value,context_.value);
      for(const auto& [factor,multiplicity]:owner_.powers_.at(index)){
        if(!nmod_mpoly_pow_ui(power_.value,factors_[factor].value,multiplicity,context_.value))return std::nullopt;
        nmod_mpoly_mul(product_.value,product_.value,power_.value,context_.value);}
      if(static_cast<std::size_t>(nmod_mpoly_length(product_.value,context_.value))!=ansatz.denominator.size())return std::nullopt;
      std::vector<Word> out;std::vector<ulong> exponent(owner_.variables_);
      for(const auto& term:ansatz.denominator){for(std::size_t v=0;v<exponent.size();++v)exponent[v]=term.at(v);
        out.push_back(nmod_mpoly_get_coeff_ui_ui(product_.value,exponent.data(),context_.value));if(!out.back())return std::nullopt;}
      nmod_t mod;nmod_init(&mod,context_.value->mod.n);const auto scale=nmod_inv(out.at(ansatz.pivot),mod);
      for(auto& c:out)c=nmod_mul(c,scale,mod);
      return out;
    }
   private:
    const DenominatorFactors& owner_;Context context_;std::vector<Polynomial> factors_;Polynomial product_,power_;bool valid_=true;
  };
 private:
  void lift() {
    Integer limit;fmpz_fdiv_q_2exp(limit.value,modulus_.value,1);fmpz_sqrt(limit.value,limit.value);fmpz_fdiv_q_2exp(limit.value,limit.value,9);
    fmpq_t q;fmpq_init(q);lifted_=!fmpz_is_zero(limit.value);
    for(auto& factor:factors_){factor.numerators.assign(factor.powers.size(),Integer());factor.denominators.assign(factor.powers.size(),Integer());
      for(std::size_t t=0;t<factor.powers.size()&&lifted_;++t){
        lifted_=fmpq_reconstruct_fmpz_2(q,factor.residues[t].value,modulus_.value,limit.value,limit.value);
        fmpz_set(factor.numerators[t].value,fmpq_numref(q));fmpz_set(factor.denominators[t].value,fmpq_denref(q));}
      if(!lifted_)break;}
    fmpq_clear(q);
  }
  // A factor: its terms, the residues of their coefficients modulo modulus_
  // (monic) and, once lifted, the rational coefficients.
  struct Factor {Occurrence where;std::vector<Power> powers;std::vector<Integer> residues,numerators,denominators;};
  std::size_t variables_=0;unsigned primes_=0;bool lifted_=false;Integer modulus_;
  std::vector<Factor> factors_;
  std::vector<std::vector<std::pair<std::size_t,std::uint32_t>>> powers_;    // of each denominator: (factor, multiplicity)
};
// The fit of a function whose denominator is known at this prime
// (`denominator`: the coefficients of the ansatz' denominator, one at the
// pivot): the numerator alone is unknown. With more rows than numerator
// terms the system is overdetermined, so a denominator that is not the
// function's leaves no solution.
inline std::optional<Image> fit_numerator(Samples samples,std::size_t coefficient,const Ansatz& ansatz,const std::vector<Word>& denominator,Word p) {
  const auto n=ansatz.numerator.size(),columns=n+1;
  if(denominator.size()!=ansatz.denominator.size()||samples.size()<columns+1)return std::nullopt;
  nmod_mat_t matrix,kernel;nmod_mat_init(matrix,samples.size(),columns,p);nmod_mat_init(kernel,columns,columns,p);
  std::vector<unsigned> degree;
  for(const auto* terms:{&ansatz.numerator,&ansatz.denominator})for(const auto& term:*terms){
    if(term.size()>degree.size())degree.resize(term.size());for(std::size_t v=0;v<term.size();++v)degree[v]=std::max(degree[v],term[v]);}
  std::vector<std::vector<Word>> powers(degree.size());for(std::size_t v=0;v<degree.size();++v)powers[v].resize(degree[v]+1);
  const auto monomial=[&](const Power& term){Word value=1;for(std::size_t v=0;v<term.size();++v)if(term[v])value=nmod_mul(value,powers[v][term[v]],matrix->mod);return value;};
  for(std::size_t i=0;i<samples.size();++i) {
    for(std::size_t v=0;v<degree.size();++v){const Word x=samples[i].point.at(v)%p;powers[v][0]=1;for(unsigned e=1;e<=degree[v];++e)powers[v][e]=nmod_mul(powers[v][e-1],x,matrix->mod);}
    for(std::size_t j=0;j<n;++j)nmod_mat_entry(matrix,i,j)=monomial(ansatz.numerator[j]);
    Word value=0;for(std::size_t j=0;j<denominator.size();++j)value=nmod_add(value,nmod_mul(denominator[j],monomial(ansatz.denominator[j]),matrix->mod),matrix->mod);
    nmod_mat_entry(matrix,i,n)=nmod_neg(nmod_mul(samples[i].coefficients.at(coefficient)%p,value,matrix->mod),matrix->mod);
  }
  const slong rank=nmod_mat_nullspace(kernel,matrix);nmod_mat_clear(matrix);
  if(rank!=1||!nmod_mat_entry(kernel,n,0)){nmod_mat_clear(kernel);return std::nullopt;}
  const auto scale=nmod_inv(nmod_mat_entry(kernel,n,0),kernel->mod);
  Image result{ansatz,{}};
  for(std::size_t j=0;j<n;++j)result.coefficients.push_back(nmod_mul(nmod_mat_entry(kernel,j,0),scale,kernel->mod));
  result.coefficients.insert(result.coefficients.end(),denominator.begin(),denominator.end());
  nmod_mat_clear(kernel);return result;
}
struct Lift {
  Ansatz ansatz;std::vector<Integer> residues;Integer modulus;
  // The coefficients of one rational function share a denominator d (the
  // integer its normalising coefficient stands for). Reconstructed one by
  // one, each needs a modulus beyond 2*|a|*d, twice the bits of the integers
  // themselves. Together they need little more than those bits: (d*r_1 mod
  // M, ..., d*r_k mod M, d) is a short vector of the lattice spanned by
  // M*e_1, ..., M*e_k and (r_1, ..., r_k, 1), of length about 2^H where the
  // lattice has determinant M^k, so that M beyond 2^(H(k+1)/k) singles it
  // out. The d found from k coefficients is accepted only if it makes every
  // coefficient of the function that small, 16 bits below what a lattice of
  // this determinant always contains; held-out primes validate the result as
  // they do any other.
  bool together(std::vector<fmpq>& out)const {
    constexpr std::size_t most=11;constexpr unsigned margin=16;
    Integer d(1),a,half,bound,factor;fmpz_fdiv_q_2exp(half.value,modulus.value,1);
    // M^(11/12) / 2^16: what the largest lattice used always contains, less 16 bits.
    fmpz_pow_ui(bound.value,modulus.value,most);fmpz_root(bound.value,bound.value,most+1);fmpz_fdiv_q_2exp(bound.value,bound.value,margin);
    std::vector<Integer> scaled(residues.size());
    const auto rescale=[&]{std::vector<std::size_t> large;
      for(std::size_t i=0;i<residues.size();++i){auto* v=scaled[i].value;fmpz_mul(v,residues[i].value,d.value);fmpz_mod(v,v,modulus.value);
        if(fmpz_cmp(v,half.value)>0)fmpz_sub(v,v,modulus.value);if(fmpz_cmpabs(v,bound.value)>=0)large.push_back(i);}
      return large;};
    // A denominator found from some coefficients need not clear the others:
    // those it leaves large give the factor that is missing.
    auto large=rescale();
    // The coefficients of a lattice are spread over the large ones. Those
    // next to each other in the ansatz (one power of a variable, successive
    // powers of the other) often are one large integer times small ones:
    // their lattice then holds vectors far shorter than the one sought, and
    // they share only a small part of the denominator. For the same reason
    // the factor is looked for in every vector of the reduced basis: the one
    // taken is the one that leaves the fewest coefficients large.
    // Where fewer coefficients are large than a lattice takes, it is filled
    // with others: the factor sought leaves those small as well, and the
    // lattice then singles it out as the largest one does.
    // Every round leaves fewer coefficients large: that ends the loop.
    while(!large.empty()) {
      std::vector<std::size_t> used;
      for(std::size_t j=0,count=std::min(most,large.size());j<count;++j)used.push_back(large[j*large.size()/count]);
      if(used.size()<most){std::vector<std::size_t> others;
        for(std::size_t i=0;i<residues.size();++i)if(!fmpz_is_zero(scaled[i].value)&&std::find(used.begin(),used.end(),i)==used.end())others.push_back(i);
        for(std::size_t j=0,count=std::min(most-used.size(),others.size());j<count;++j)used.push_back(others[j*others.size()/count]);}
      const auto k=static_cast<slong>(used.size());
      fmpz_mat_t basis;fmpz_mat_init(basis,k+1,k+1);
      for(slong j=0;j<k;++j){fmpz_set(fmpz_mat_entry(basis,j,j),modulus.value);fmpz_mod(fmpz_mat_entry(basis,k,j),scaled[used[j]].value,modulus.value);}
      fmpz_one(fmpz_mat_entry(basis,k,k));
      fmpz_lll_t context;fmpz_lll_context_init_default(context);fmpz_lll(basis,nullptr,context);
      Integer best,candidate,product,value;std::size_t fewest=large.size();
      for(slong row=0;row<=k;++row){
        fmpz_abs(candidate.value,fmpz_mat_entry(basis,row,k));if(fmpz_is_zero(candidate.value)||fmpz_is_one(candidate.value))continue;
        fmpz_mul(product.value,d.value,candidate.value);if(fmpz_cmpabs(product.value,bound.value)>=0)continue;
        // The count stops once it cannot fall below the best so far.
        std::size_t left=0;
        for(const auto& old:scaled){auto* v=value.value;fmpz_mul(v,old.value,candidate.value);fmpz_mod(v,v,modulus.value);
          if(fmpz_cmp(v,half.value)>0)fmpz_sub(v,v,modulus.value);if(fmpz_cmpabs(v,bound.value)>=0&&++left>=fewest)break;}
        if(left<fewest){fewest=left;fmpz_set(best.value,candidate.value);}
      }
      fmpz_mat_clear(basis);
      if(fmpz_is_zero(best.value))return false;
      fmpz_set(factor.value,best.value);fmpz_mul(d.value,d.value,factor.value);
      large=rescale();
      // A smaller lattice always contains a shorter vector: what it found
      // must be as far below that as well.
      if(static_cast<std::size_t>(k)<most){Integer small;fmpz_pow_ui(small.value,modulus.value,static_cast<ulong>(k));fmpz_root(small.value,small.value,k+1);fmpz_fdiv_q_2exp(small.value,small.value,margin);
        if(fmpz_cmpabs(factor.value,small.value)>=0)return false;
        for(auto i:used)if(fmpz_cmpabs(scaled[i].value,small.value)>=0)return false;}
    }
    for(std::size_t i=0;i<residues.size();++i)fmpq_set_fmpz_frac(&out[i],scaled[i].value,d.value);
    return true;
  }
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
    if(!together(values.q))for(std::size_t i=0;i<residues.size();++i)
      if(!fmpq_reconstruct_fmpz(&values.q[i],residues[i].value,modulus.value))return std::nullopt;
    for(std::size_t i=0;i<residues.size();++i)fmpz_lcm(values.scale,values.scale,fmpq_denref(&values.q[i]));
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
