#pragma once
#include "ibp/solver.hpp"
#include <flint/ulong_extras.h>
#include <flint/nmod.h>
#include <array>
#include <bit>

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
// Arithmetic: a row under reduction accumulates the products factor*entry of
// a lane as a two-word integer below prime*2^64 and is reduced modulo the
// prime once, when its column is read; a product is then one multiplication
// with no division. The prime is any the Field takes (below 2^61).
template<unsigned K> class BatchSolver {
 public:
  static constexpr unsigned lanes=K;
  struct Limits {std::size_t max_fill=std::numeric_limits<std::size_t>::max();double seconds=1e9;};
  using Values=std::array<Word,K>;
  struct Reduced {std::vector<Column> columns;std::vector<Word> values;}; // values[K*i+lane] of columns[i]
  BatchSolver(std::size_t columns,const Field& field):p_(field.prime()),lookup_(columns,-1),slot_(columns,none),occupied_((columns+63)/64,0),words_((occupied_.size()+63)/64,0){
    static_assert(K>=1&&K<=64);nmod_init(&mod_,p_);}
  // Eliminates the selected equations in order. `forms` holds the K values of
  // every form, forms[K*form+lane]. complete_rows as SolveOptions::complete_rows.
  void solve(const FormRows& rows,const Word* forms,const std::vector<std::uint32_t>& selected,bool complete_rows,const Limits& limits={}) {
    started_=std::chrono::steady_clock::now();limits_=limits;
    for(const auto& pivot:pivots_)lookup_[columns_[pivot.first]]=-1;
    pivots_.clear();columns_.clear();values_.clear();stats_={};failed_=0;
    try{for(const auto source:selected){
      if((stats_.equations++&255)==0)budget();
      const auto& terms=rows.at(source);if(terms.empty())continue;
      // An equation whose first column has no pivot is stored as it is.
      Column low=terms.front().column,high=terms.back().column;
      for(const auto& term:terms){
        if(term.column>=lookup_.size())throw std::invalid_argument("reduction column out of range");
        room(1);Word* a=cell(term.column);const Word* v=forms+K*std::size_t(term.form);
        for(unsigned k=0;k<K;++k){a[2*k]=v[k];a[2*k+1]=0;}
      }
      insert(low,high,complete_rows);
    }}catch(...){clear();throw;}
    stats_.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count();
  }
  // The reduction of a column by the pivots: the columns left, ascending, with
  // their K values. A column whose value is zero at a lane is not part of
  // that lane's reduction.
  void reduce(Column target,Reduced& out) {
    if(target>=lookup_.size())throw std::invalid_argument("reduction column out of range");
    out.columns.clear();out.values.clear();
    room(1);Word* a=cell(target);for(unsigned k=0;k<K;++k){a[2*k]=1;a[2*k+1]=0;}
    Column high=target;Values value;std::size_t steps=0;
    try{for(Column column=target;(column=next(column,high))!=none;++column){
      const bool any=take(column,value);
      if(!any)continue;
      if((steps++&1023)==0)budget();
      if(lookup_[column]>=0){eliminate(pivots_[lookup_[column]],value,high);continue;}
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
  // first: where the row starts in columns_ (and, times K, in values_); second: its length.
  using Pivot=std::pair<std::size_t,std::size_t>;
  void budget(){
    if(stats_.fill>limits_.max_fill)throw std::length_error("pivot fill budget exceeded");
    if(std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count()>limits_.seconds)throw std::runtime_error("elimination time budget exceeded");
  }
  // After a budget ended a reduction: nothing is left of its row.
  void clear(){std::fill(acc_.begin(),acc_.end(),0);std::fill(occupied_.begin(),occupied_.end(),0);std::fill(words_.begin(),words_.end(),0);std::fill(slot_.begin(),slot_.end(),none);used_=0;}
  // Room for `more` further columns in the row under reduction.
  void room(std::size_t more){if(2*K*(used_+more)>acc_.size())acc_.resize(std::max(2*acc_.size(),2*K*(used_+more)),0);}
  // The 2K words of a column of the row under reduction; a column not yet
  // occupied gets the next free place, which holds zero.
  Word* cell(Column column){
    auto& slot=slot_[column];
    if(slot==none){slot=static_cast<Column>(used_++);occupied_[column>>6]|=std::uint64_t(1)<<(column&63);words_[column>>12]|=std::uint64_t(1)<<((column>>6)&63);}
    return acc_.data()+2*K*std::size_t(slot);
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
  // Reads the K residues of a column of the row under reduction and leaves
  // the column empty. Whether any of them is nonzero.
  bool take(Column column,Values& value) {
    Word* a=acc_.data()+2*K*std::size_t(slot_[column]);Word any=0;
    for(unsigned k=0;k<K;++k){
      // The high word is below the prime.
      const Word hi=a[2*k+1],lo=a[2*k];Word r;
      if(hi==0&&lo<p_)r=lo;else NMOD_RED2(r,hi,lo,mod_);
      value[k]=r;any|=r;a[2*k]=0;a[2*k+1]=0;
    }
    if(!(occupied_[column>>6]&=~(std::uint64_t(1)<<(column&63))))words_[column>>12]&=~(std::uint64_t(1)<<((column>>6)&63));
    slot_[column]=none;
    return any!=0;
  }
  // Subtracts value*row from the row under reduction; the row is monic and
  // its own column, already taken, is left out.
  void eliminate(const Pivot& pivot,const Values& value,Column& high) {
    Values factor;for(unsigned k=0;k<K;++k)factor[k]=p_-value[k];
    const Column* c=columns_.data()+pivot.first;const Word* v=values_.data()+K*pivot.first;const std::size_t n=pivot.second;
    const Word p=p_;room(n);
    for(std::size_t j=1;j<n;++j){
      Word* a=cell(c[j]);const Word* e=v+K*j;
      for(unsigned k=0;k<K;++k){
        const unsigned __int128 product=static_cast<unsigned __int128>(factor[k])*e[k];
        const unsigned __int128 sum=product+((static_cast<unsigned __int128>(a[2*k+1])<<64)|a[2*k]);
        Word hi=static_cast<Word>(sum>>64);hi=hi>=p?hi-p:hi;
        a[2*k]=static_cast<Word>(sum);a[2*k+1]=hi;
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
      if((!led||complete_rows)&&lookup_[column]>=0){eliminate(pivots_[lookup_[column]],value,high);if((steps++&1023)==0)budget();continue;}
      if(!led){led=true;inverse=value;invert(inverse);for(unsigned k=0;k<K;++k)value[k]=1;}
      else for(unsigned k=0;k<K;++k)value[k]=nmod_mul(value[k],inverse[k],mod_);
      columns_.push_back(column);values_.insert(values_.end(),value.begin(),value.end());
    }
    used_=0;
    const std::size_t length=columns_.size()-start;if(!length)return;
    stats_.peak_row=std::max(stats_.peak_row,length);stats_.fill+=length;++stats_.rank;
    lookup_[columns_[start]]=static_cast<int>(pivots_.size());pivots_.push_back({start,length});
  }
  Word p_;nmod_t mod_;Limits limits_;std::chrono::steady_clock::time_point started_;
  std::vector<int> lookup_;std::vector<Pivot> pivots_;std::vector<Column> columns_;std::vector<Word> values_;
  // The row under reduction. Its columns take places in acc_ in the order
  // they become occupied (slot_), so that what a reduction touches lies
  // together whatever the number of columns of the system: for the column at
  // place s and lane k the two words acc_[2*(K*s+k)] (low) and
  // acc_[2*(K*s+k)+1] (high). One bit per occupied column gives their order.
  // Between reductions no column has a place and all of acc_ is zero.
  std::vector<Word> acc_;std::vector<Column> slot_;std::vector<std::uint64_t> occupied_,words_;std::size_t used_=0;
  Statistics stats_;std::uint64_t failed_=0;
};
}
