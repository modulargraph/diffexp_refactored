#pragma once
#include "ibp/solver.hpp"
#include <flint/ulong_extras.h>
#include <flint/nmod.h>
#include <array>
#include <bit>
#include <cstdlib>
#include <cstring>
#include <new>
#if defined(__linux__)
#include <sys/mman.h>
#endif
#if defined(__AVX512F__)&&defined(__AVX512IFMA__)&&!defined(IBP_BATCH_PORTABLE)
#include <immintrin.h>
#define IBP_BATCH_IFMA 1
#endif

namespace ibp {
// The elimination of Solver at K points of one prime at once.
//
// What an elimination does besides arithmetic (which pivot a column has, the
// order of the columns of a row, where a row is stored, which columns a row
// occupies) is generically the same at every point. Here it is done once for
// K points: an entry of a row is a column and K values, one per point (a
// lane), and a row under reduction K values per column.
//
// A pivot found here is a column where the row is nonzero at some lane. At
// a lane where it is zero there, the elimination at that point alone would
// have taken another column: that lane is marked failed (failed()), what is
// computed for it is not used, and its point is left to Solver. At every
// other lane each step is the step Solver takes at that point, in the same
// order, so that the reduced targets are the same residues.
//
// Arithmetic. A row under reduction accumulates the products factor*entry of
// a lane unreduced and is reduced modulo the prime once, when its column is
// read: a product costs no division. The prime is any the Field takes
// (below 2^61). There are two forms of the accumulator, with the same
// results:
//  - two words a lane, an integer below prime*2^64 (one 64x64 multiplication,
//    an addition with carry and a conditional subtraction per product);
//  - with AVX-512 IFMA and K a multiple of 8, three words a lane holding
//    sums of 52-bit halves of the four partial products of the 52-bit and
//    9-bit parts of factor and entry: seven multiply-adds for eight lanes.
//    A word takes 1024 such sums; a row reduced by more pivots than that is
//    brought back to residues in between.
template<unsigned K> class BatchSolver {
 public:
  static constexpr unsigned lanes=K;
#ifdef IBP_BATCH_IFMA
  static constexpr bool wide=K%8==0;
#else
  static constexpr bool wide=false;
#endif
  struct Limits {std::size_t max_fill=std::numeric_limits<std::size_t>::max();double seconds=1e9;};
  using Values=std::array<Word,K>;
  struct Reduced {std::vector<Column> columns;std::vector<Word> values;}; // values[K*i+lane] of columns[i]
  BatchSolver(std::size_t columns,const Field& field):p_(field.prime()),lookup_(columns,-1),slot_(columns,none),occupied_((columns+63)/64,0),words_((occupied_.size()+63)/64,0){
    static_assert(K>=1&&K<=64);nmod_init(&mod_,p_);
    // 2^104 modulo the prime.
    carry_=nmod_pow_ui(2%p_,104,mod_);
    // The quotient of 2^106 by the prime, for primes of 61 bits (see residues).
    if(p_>>60)quotient_=static_cast<Word>((static_cast<unsigned __int128>(1)<<106)/p_);
  }
  // Eliminates the selected equations in order. `forms` holds the K values of
  // every form, forms[K*form+lane]. complete_rows as SolveOptions::complete_rows.
  void solve(const FormRows& rows,const Word* forms,const std::vector<std::uint32_t>& selected,bool complete_rows,const Limits& limits={}) {
    started_=std::chrono::steady_clock::now();limits_=limits;
    for(const auto& pivot:pivots_)lookup_[columns_[pivot.first]]=-1;
    pivots_.clear();columns_.clear();values_.size=0;stats_={};failed_=0;
    try{for(const auto source:selected){
      if((stats_.equations++&255)==0)budget();
      const auto& terms=rows.at(source);if(terms.empty())continue;
      Column low=terms.front().column,high=terms.back().column;
      acc_.room(stride*(used_+terms.size()));
      for(const auto& term:terms){
        if(term.column>=lookup_.size())throw std::invalid_argument("reduction column out of range");
        // The cell of a column is zero: its first word takes the residue.
        std::memcpy(cell(term.column),forms+K*std::size_t(term.form),K*sizeof(Word));
      }
      sums_=1;insert(low,high,complete_rows);
    }}catch(...){clear();throw;}
    stats_.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count();
  }
  // The reduction of a column by the pivots: the columns left, ascending, with
  // their K values. A column whose value is zero at a lane is not part of
  // that lane's reduction.
  void reduce(Column target,Reduced& out) {
    if(target>=lookup_.size())throw std::invalid_argument("reduction column out of range");
    out.columns.clear();out.values.clear();
    acc_.room(stride*(used_+1));Word* a=cell(target);for(unsigned k=0;k<K;++k)a[k]=1;
    Column high=target;Values value;std::size_t steps=0;sums_=1;
    try{for(Column column=target;(column=next(column,high))!=none;++column){
      const bool any=take(column,value);
      if(!any)continue;
      if((steps++&1023)==0)budget();
      if(lookup_[column]>=0){ahead(column,high);eliminate(pivots_[lookup_[column]],value,column,high);continue;}
      out.columns.push_back(column);out.values.insert(out.values.end(),value.begin(),value.end());
    }}catch(...){clear();throw;}
    used_=0;
  }
  // Whether a column has a pivot.
  bool pivot(Column column)const{return lookup_.at(column)>=0;}
  // Lanes whose point does not have the pivots found here, one bit each.
  std::uint64_t failed()const{return failed_;}
  const Statistics& statistics()const{return stats_;}
 private:
  static constexpr Column none=std::numeric_limits<Column>::max();
  // The words of a column of the row under reduction: `limbs` for each lane,
  // the first word of all lanes, then the second, ...
  static constexpr unsigned limbs=wide?3:2,stride=limbs*K;
  static constexpr Word low52=(Word(1)<<52)-1;
  // Words on 64-byte boundaries that keep their contents when they grow;
  // new words are zero.
  struct Words {
    Word* data=nullptr;std::size_t size=0,capacity=0;
    Words()=default;Words(const Words&)=delete;Words& operator=(const Words&)=delete;~Words(){std::free(data);}
    void room(std::size_t wanted){
      if(wanted<=capacity)return;const std::size_t grown=std::max<std::size_t>({wanted,2*capacity,1024});
      // Large blocks on huge pages where the system gives them on request:
      // pivot rows are read in no order, a row at a time.
      const std::size_t bytes=grown*sizeof(Word),huge=std::size_t(1)<<21,alignment=bytes>=2*huge?huge:64;
      auto* fresh=static_cast<Word*>(std::aligned_alloc(alignment,(bytes+alignment-1)/alignment*alignment));if(!fresh)throw std::bad_alloc();
#if defined(__linux__)&&defined(MADV_HUGEPAGE)&&!defined(IBP_NO_HUGE)
      if(alignment==huge)madvise(fresh,(bytes+alignment-1)/alignment*alignment,MADV_HUGEPAGE);
#endif
      if(capacity)std::memcpy(fresh,data,capacity*sizeof(Word));std::memset(fresh+capacity,0,(grown-capacity)*sizeof(Word));std::free(data);data=fresh;capacity=grown;
    }
  };
  // first: where the row starts in columns_ (and, times K, in values_); second: its length.
  using Pivot=std::pair<std::size_t,std::size_t>;
  void budget(){
    if(stats_.fill>limits_.max_fill)throw std::length_error("pivot fill budget exceeded");
    if(std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count()>limits_.seconds)throw std::runtime_error("elimination time budget exceeded");
  }
  // After a budget ended a reduction: nothing is left of its row.
  void clear(){if(acc_.capacity)std::memset(acc_.data,0,acc_.capacity*sizeof(Word));std::fill(occupied_.begin(),occupied_.end(),0);std::fill(words_.begin(),words_.end(),0);std::fill(slot_.begin(),slot_.end(),none);used_=0;}
  // The words of a column of the row under reduction; a column not yet
  // occupied gets the next free place, which holds zero. The caller has made
  // room for it.
  Word* cell(Column column){
    auto& slot=slot_[column];
    if(slot==none){slot=static_cast<Column>(used_++);occupied_[column>>6]|=std::uint64_t(1)<<(column&63);words_[column>>12]|=std::uint64_t(1)<<((column>>6)&63);}
    return acc_.data+stride*std::size_t(slot);
  }
  // The first occupied column in [from,high]. The columns of a row lie far
  // apart: words_ has one bit per nonzero word of occupied_.
  Column next(Column from,Column high)const {
    if(from>high)return none;std::size_t w=from>>6;
    if(const std::uint64_t word=occupied_[w]&(~std::uint64_t(0)<<(from&63))){const Column c=static_cast<Column>((w<<6)+std::countr_zero(word));return c<=high?c:none;}
    ++w;std::size_t g=w>>6;if(g>=words_.size())return none;std::uint64_t group=words_[g]&(~std::uint64_t(0)<<(w&63));
    for(;;){
      if(group){w=(g<<6)+std::countr_zero(group);const Column c=static_cast<Column>((w<<6)+std::countr_zero(occupied_[w]));return c<=high?c:none;}
      if(++g>(high>>12))return none;group=words_[g];
    }
  }
  // The residue the words of one lane stand for.
  Word residue(const Word* a,unsigned k)const {
    if constexpr(wide){
      // first + second*2^52 + third*2^104, the third below 2^29.
      const Word first=a[k],second=a[K+k],third=a[2*K+k];
      if(!(second|third))return first<p_?first:first%p_;
      const unsigned __int128 sum=static_cast<unsigned __int128>(first)+(static_cast<unsigned __int128>(second)<<52)+static_cast<unsigned __int128>(third)*carry_;
      Word hi=static_cast<Word>(sum>>64),r;const Word lo=static_cast<Word>(sum);if(hi>=p_)hi%=p_;
      NMOD_RED2(r,hi,lo,mod_);return r;
    }else{
      // The high word is below the prime.
      const Word lo=a[k],hi=a[K+k];Word r;
      if(hi==0&&lo<p_)return lo;
      NMOD_RED2(r,hi,lo,mod_);return r;
    }
  }
  // The pivot row of the next occupied column is asked for while the present
  // one is subtracted: rows are read in no order, and the next row operation
  // would otherwise begin by waiting for its row.
  void ahead(Column column,Column high)const {
#ifndef IBP_NO_PREFETCH
    const Column following=next(column+1,high);if(following==none)return;
    const int p=lookup_[following];if(p<0)return;
    const auto first=pivots_[p].first;const char* v=reinterpret_cast<const char*>(values_.data+K*first);
    __builtin_prefetch(columns_.data()+first);for(unsigned line=0;line<4;++line)__builtin_prefetch(v+64*line);
#else
    (void)column;(void)high;
#endif
  }
#ifdef IBP_BATCH_IFMA
  // The residues of eight lanes of three words (see residue), for a prime of
  // 61 bits. With the words brought to 52 bits each, the third is folded
  // into the others by 2^104 modulo the prime, which leaves an integer z
  // below 2^106. An estimate q of its quotient by the prime follows from the
  // top 47 bits of z and the quotient of 2^106 by the prime; it is at most
  // two too small, so z-q*prime, computed modulo 2^64, is below three times
  // the prime.
  __m512i residues(__m512i first,__m512i second,__m512i third)const {
    const __m512i mask=_mm512_set1_epi64(static_cast<long long>(low52)),prime=_mm512_set1_epi64(static_cast<long long>(p_));
    const __m512i fold0=_mm512_set1_epi64(static_cast<long long>(carry_&low52)),fold1=_mm512_set1_epi64(static_cast<long long>(carry_>>52));
    const __m512i c0=_mm512_and_si512(first,mask),t=_mm512_add_epi64(second,_mm512_srli_epi64(first,52));
    const __m512i c1=_mm512_and_si512(t,mask),c2=_mm512_add_epi64(third,_mm512_srli_epi64(t,52));
    const __m512i d0=_mm512_madd52lo_epu64(c0,c2,fold0);
    const __m512i d1=_mm512_madd52lo_epu64(_mm512_madd52hi_epu64(c1,c2,fold0),c2,fold1);
    const __m512i g0=_mm512_and_si512(d0,mask),u=_mm512_add_epi64(d1,_mm512_srli_epi64(d0,52));
    const __m512i g1=_mm512_and_si512(u,mask),g2=_mm512_srli_epi64(u,52);
    const __m512i top=_mm512_or_si512(_mm512_srli_epi64(g1,7),_mm512_slli_epi64(g2,45));
    const __m512i zero=_mm512_setzero_si512(),estimate=_mm512_set1_epi64(static_cast<long long>(quotient_));
    const __m512i q=_mm512_or_si512(_mm512_slli_epi64(_mm512_madd52hi_epu64(zero,top,estimate),5),_mm512_srli_epi64(_mm512_madd52lo_epu64(zero,top,estimate),47));
    __m512i r=_mm512_sub_epi64(_mm512_add_epi64(g0,_mm512_slli_epi64(g1,52)),_mm512_mullo_epi64(q,prime));
    r=_mm512_min_epu64(r,_mm512_sub_epi64(r,prime));return _mm512_min_epu64(r,_mm512_sub_epi64(r,prime));
  }
  // The products of eight lanes, as three words a lane.
  static __m512i lower(__m512i x){return _mm512_and_si512(x,_mm512_set1_epi64(static_cast<long long>(low52)));}
  static __m512i upper(__m512i x){return _mm512_srli_epi64(x,52);}
  __m512i product(__m512i a,__m512i b)const {
    const __m512i zero=_mm512_setzero_si512(),a0=lower(a),a1=upper(a),b0=lower(b),b1=upper(b);
    const __m512i first=_mm512_madd52lo_epu64(zero,a0,b0);
    const __m512i second=_mm512_madd52lo_epu64(_mm512_madd52lo_epu64(_mm512_madd52hi_epu64(zero,a0,b0),a0,b1),a1,b0);
    const __m512i third=_mm512_madd52lo_epu64(_mm512_madd52hi_epu64(_mm512_madd52hi_epu64(zero,a0,b1),a1,b0),a1,b1);
    return residues(first,second,third);
  }
#endif
  // The residues of all lanes of a column of the row under reduction.
  Word residues(const Word* a,Values& value)const {
#ifdef IBP_BATCH_IFMA
    if constexpr(wide)if(quotient_){__m512i any=_mm512_setzero_si512();
      for(unsigned b=0;b<K/8;++b){const __m512i r=residues(_mm512_load_si512(a+8*b),_mm512_load_si512(a+K+8*b),_mm512_load_si512(a+2*K+8*b));
        _mm512_storeu_si512(value.data()+8*b,r);any=_mm512_or_si512(any,r);}
      return _mm512_test_epi64_mask(any,any)?1:0;}
#endif
    Word any=0;for(unsigned k=0;k<K;++k){const Word r=residue(a,k);value[k]=r;any|=r;}
    return any;
  }
  // value*factor modulo the prime, lane by lane.
  void scale(Values& value,const Values& factor)const {
#ifdef IBP_BATCH_IFMA
    if constexpr(wide)if(quotient_){for(unsigned b=0;b<K/8;++b)_mm512_storeu_si512(value.data()+8*b,product(_mm512_loadu_si512(value.data()+8*b),_mm512_loadu_si512(factor.data()+8*b)));return;}
#endif
    for(unsigned k=0;k<K;++k)value[k]=nmod_mul(value[k],factor[k],mod_);
  }
  // Reads the K residues of a column of the row under reduction and leaves
  // the column empty. Whether any of them is nonzero.
  bool take(Column column,Values& value) {
    Word* a=acc_.data+stride*std::size_t(slot_[column]);const Word any=residues(a,value);
    std::memset(a,0,stride*sizeof(Word));
    if(!(occupied_[column>>6]&=~(std::uint64_t(1)<<(column&63))))words_[column>>12]&=~(std::uint64_t(1)<<((column>>6)&63));
    slot_[column]=none;
    return any!=0;
  }
  // Subtracts value*row from the row under reduction; the row is monic and
  // its own column, already taken, is left out. `current` is that column:
  // what is occupied lies after it.
  void eliminate(const Pivot& pivot,const Values& value,Column current,Column& high) {
    const Column* c=columns_.data()+pivot.first;const Word* v=values_.data+K*pivot.first;const std::size_t n=pivot.second;
    acc_.room(stride*(used_+n));
#ifdef IBP_BATCH_IFMA
    if constexpr(wide){
      const __m512i mask=_mm512_set1_epi64(static_cast<long long>(low52));
      __m512i f0[K/8],f1[K/8];
      for(unsigned b=0;b<K/8;++b){alignas(64) Word factor[8];for(unsigned k=0;k<8;++k)factor[k]=p_-value[8*b+k];
        const __m512i f=_mm512_load_si512(factor);f0[b]=_mm512_and_si512(f,mask);f1[b]=_mm512_srli_epi64(f,52);}
      for(std::size_t j=1;j<n;++j){
        Word* a=cell(c[j]);const Word* e=v+K*j;
        for(unsigned b=0;b<K/8;++b){
          const __m512i x=_mm512_load_si512(e+8*b),x0=_mm512_and_si512(x,mask),x1=_mm512_srli_epi64(x,52);
          __m512i first=_mm512_load_si512(a+8*b),second=_mm512_load_si512(a+K+8*b),third=_mm512_load_si512(a+2*K+8*b);
          first=_mm512_madd52lo_epu64(first,f0[b],x0);
          second=_mm512_madd52hi_epu64(second,f0[b],x0);second=_mm512_madd52lo_epu64(second,f0[b],x1);second=_mm512_madd52lo_epu64(second,f1[b],x0);
          third=_mm512_madd52hi_epu64(third,f0[b],x1);third=_mm512_madd52hi_epu64(third,f1[b],x0);third=_mm512_madd52lo_epu64(third,f1[b],x1);
          _mm512_store_si512(a+8*b,first);_mm512_store_si512(a+K+8*b,second);_mm512_store_si512(a+2*K+8*b,third);
        }
      }
      high=std::max(high,c[n-1]);++stats_.eliminations;stats_.multiplications+=K*(n-1);
      // A word has taken up to three 52-bit halves from this row.
      if(++sums_>=1024){Values rest;
        for(Column column=current+1;(column=next(column,high))!=none;++column){Word* a=acc_.data+stride*std::size_t(slot_[column]);
          residues(a,rest);std::memset(a,0,stride*sizeof(Word));std::memcpy(a,rest.data(),K*sizeof(Word));}
        sums_=1;}
      return;
    }
#endif
    (void)current;
    Values factor;for(unsigned k=0;k<K;++k)factor[k]=p_-value[k];
    const Word p=p_;
    for(std::size_t j=1;j<n;++j){
      Word* a=cell(c[j]);const Word* e=v+K*j;
      for(unsigned k=0;k<K;++k){
        const unsigned __int128 product=static_cast<unsigned __int128>(factor[k])*e[k];
        const unsigned __int128 sum=product+((static_cast<unsigned __int128>(a[K+k])<<64)|a[k]);
        Word hi=static_cast<Word>(sum>>64);hi=hi>=p?hi-p:hi;
        a[k]=static_cast<Word>(sum);a[K+k]=hi;
      }
    }
    high=std::max(high,c[n-1]);++stats_.eliminations;stats_.multiplications+=K*(n-1);
  }
  // The inverses of K values with one inversion; a zero is marked failed
  // and inverted as if it were one.
  void invert(Values& value) {
    Values prefix;Word product=1;
    for(unsigned k=0;k<K;++k){if(!value[k]){failed_|=std::uint64_t(1)<<k;value[k]=1;}prefix[k]=product;product=nmod_mul(product,value[k],mod_);}
    Word inverse=n_invmod(product,p_);
    for(unsigned k=K;k-->0;){const Word own=nmod_mul(inverse,prefix[k],mod_);inverse=nmod_mul(inverse,value[k],mod_);value[k]=own;}
  }
  void insert(Column low,Column high,bool complete_rows) {
    const std::size_t start=columns_.size();Values value,inverse;bool led=false;std::size_t steps=0;
    for(Column column=low;(column=next(column,high))!=none;++column){
      const bool any=take(column,value);
      if(!any)continue;
      if((!led||complete_rows)&&lookup_[column]>=0){ahead(column,high);eliminate(pivots_[lookup_[column]],value,column,high);if((steps++&1023)==0)budget();continue;}
      if(!led){led=true;inverse=value;invert(inverse);for(unsigned k=0;k<K;++k)value[k]=1;}
      else scale(value,inverse);
      columns_.push_back(column);values_.room(values_.size+K);std::memcpy(values_.data+values_.size,value.data(),K*sizeof(Word));values_.size+=K;
    }
    used_=0;
    const std::size_t length=columns_.size()-start;if(!length)return;
    stats_.peak_row=std::max(stats_.peak_row,length);stats_.fill+=length;++stats_.rank;
    lookup_[columns_[start]]=static_cast<int>(pivots_.size());pivots_.push_back({start,length});
  }
  Word p_,carry_=0,quotient_=0;nmod_t mod_;Limits limits_;std::chrono::steady_clock::time_point started_;
  std::vector<int> lookup_;std::vector<Pivot> pivots_;std::vector<Column> columns_;
  // The values of the pivot rows: K words for every entry of columns_.
  Words values_;
  // The row under reduction. Its columns take places in acc_ in the order
  // they become occupied (slot_), so that what a reduction touches lies
  // together whatever the number of columns of the system. One bit per
  // occupied column gives their order. Between reductions no column has a
  // place and all of acc_ is zero.
  Words acc_;std::vector<Column> slot_;std::vector<std::uint64_t> occupied_,words_;std::size_t used_=0;
  // Row operations since the words of the row were last residues.
  unsigned sums_=0;
  Statistics stats_;std::uint64_t failed_=0;
};
}
