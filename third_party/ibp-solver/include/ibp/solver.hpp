#pragma once
#include "ibp/geometry.hpp"
#include <limits>
#include <unordered_map>
#include <chrono>
#include <unordered_set>
#include <functional>
namespace ibp {
using Column=std::uint32_t;
struct LinearEntry {Column column;Word constant,dimension;};
struct Entry {Column column;Word value;bool operator==(const Entry&)const=default;};
// A term of an equation whose coefficient is one of a table of values (the
// forms of a ParametricProgram at a point): eight bytes instead of the
// twenty-four of a LinearEntry, and no copy of the equations per point.
struct FormTerm {Column column;std::uint32_t form;};
using FormRows=std::vector<std::vector<FormTerm>>;
using Row=std::vector<Entry>;
struct SeedOptions {unsigned dots=1,numerators=1;std::size_t max_seeds=100000,max_terms=10000000,max_states=2000000;double seconds=120;};
struct Program {
 std::vector<Integral> integrals,targets;std::vector<std::vector<LinearEntry>> equations;
 std::vector<Column> target_columns;std::size_t seeds=0;
};
inline Program generate(const Geometry& g,const Field& f,const SeedOptions& options,std::vector<Integral> targets={}) {
 if(options.dots>Integral::max_seed_power || options.numerators>Integral::max_seed_power)throw std::invalid_argument("seed powers exceed what an integral stores");
 const auto generation_start=std::chrono::steady_clock::now();std::size_t states=0;
 Program out;std::unordered_map<Integral,Column,IntegralHash> ids;
 auto intern=[&](const Integral& a){auto [it,inserted]=ids.try_emplace(a,out.integrals.size());if(inserted)out.integrals.push_back(a);return it->second;};
 if(targets.empty()){Integral top;for(unsigned k=0;k<g.physical;++k)top.powers[k]=1;targets.push_back(top);for(unsigned k=0;k<g.n;++k){auto target=top;target.powers[k]+=k<g.physical?1:-1;targets.push_back(target);}}
 out.targets=targets;for(const auto& target:targets){for(unsigned k=g.physical;k<Integral::slots;++k)if(target.powers[k]>0 || (k>=g.n && target.powers[k]))throw std::invalid_argument("invalid target numerator slots");out.target_columns.push_back(intern(target));}
 for(const auto& target:targets)if(g.zero_sectors[sector(target,g.physical)])out.equations.push_back({{intern(target),1,0}});
 std::size_t terms=0;Integral seed;
 auto consume=[&]{if(g.zero_sectors[sector(seed,g.physical)])return;if(++out.seeds>options.max_seeds)throw std::length_error("seed budget exceeded");
  for(unsigned op=0;op<g.contractions.size();++op){std::vector<LinearEntry> raw;
   auto add=[&](const Integral& a,Word c,Word d=0){if((c||d)&&!g.zero_sectors[sector(a,g.physical)])raw.push_back({intern(a),c,d});};
   if(g.trace[op])add(seed,0,1);
   for(unsigned k=0;k<g.n;++k)if(seed.powers[k]){auto raised=seed;++raised.powers[k];Word multiplier=f.integer(-seed.powers[k]);const auto& contraction=g.contractions[op][k];add(raised,f.mul(multiplier,contraction[g.n]));for(unsigned j=0;j<g.n;++j)if(contraction[j]){auto lowered=raised;--lowered.powers[j];add(lowered,f.mul(multiplier,contraction[j]));}}
   std::sort(raw.begin(),raw.end(),[](auto a,auto b){return a.column<b.column;});std::vector<LinearEntry> row;
   for(const auto& v:raw){if(!row.empty()&&row.back().column==v.column){row.back().constant=f.add(row.back().constant,v.constant);row.back().dimension=f.add(row.back().dimension,v.dimension);}else row.push_back(v);}
   std::erase_if(row,[](auto v){return !v.constant&&!v.dimension;});if(!row.empty()){terms+=row.size();if(terms>options.max_terms)throw std::length_error("generated term budget exceeded");out.equations.push_back(std::move(row));}
  }
 };
 auto enumerate=[&](auto&& self,unsigned k,unsigned dots,unsigned nums)->void {
  if(++states>options.max_states)throw std::length_error("seed enumeration state budget exceeded");if((states&4095)==0 && std::chrono::duration<double>(std::chrono::steady_clock::now()-generation_start).count()>options.seconds)throw std::runtime_error("seed generation time budget exceeded");
  if(k==g.n){consume();return;}seed.powers[k]=0;self(self,k+1,dots,nums);
  for(unsigned s=1;s<=nums;++s){seed.powers[k]=-static_cast<int>(s);self(self,k+1,dots,nums-s);}
  if(k<g.physical)for(unsigned d=0;d<=dots;++d){seed.powers[k]=1+d;self(self,k+1,dots-d,nums);}
 };enumerate(enumerate,0,options.dots,options.numerators);
 std::vector<Column> permutation(out.integrals.size()),inverse(out.integrals.size());std::iota(permutation.begin(),permutation.end(),0);
 std::sort(permutation.begin(),permutation.end(),[&](auto a,auto b){return grade(out.integrals[a])>grade(out.integrals[b]);});auto old=out.integrals;for(unsigned i=0;i<permutation.size();++i){out.integrals[i]=old[permutation[i]];inverse[permutation[i]]=i;}
 for(auto& row:out.equations){for(auto& v:row)v.column=inverse[v.column];std::sort(row.begin(),row.end(),[](auto a,auto b){return a.column<b.column;});}
 for(auto& column:out.target_columns)column=inverse[column];
 // Easy sectors first; sparse equations first within one leading integral.
 std::stable_sort(out.equations.begin(),out.equations.end(),[](const auto& a,const auto& b){return a.front().column==b.front().column?a.size()<b.size():a.front().column>b.front().column;});
 return out;
}
struct Statistics {std::size_t equations=0,rank=0,eliminations=0,fill=0,peak_row=0;double seconds=0;};
// `progress`, if set, is told the state of a long elimination every
// `progress_seconds`: its statistics so far and the seconds it has run.
// Two ways to store a pivot row, with the same pivot columns, the same rank
// and the same reduced targets.
//  - Leading (complete_rows false): an equation is reduced until its leading
//    column has no pivot and stored with its tail as it is. `dependencies`
//    are the columns eliminated on the way. Little is computed that a target
//    does not need, which makes this the elimination of a few selected
//    equations and of an arithmetic trace; of a whole system it leaves long
//    tails that every later row and every target walks through again.
//  - Complete (complete_rows true): an equation is reduced by every pivot
//    before it, so that besides its own column a row holds only columns that
//    had no pivot then. Equations come with the simplest leading integral
//    first, so the rows it is reduced by are short already. This is the
//    elimination of a whole system: far less fill and time. Which equations
//    the targets need is then found by redoing, for them alone, what the
//    leading way would have done (Solver::selection).
struct SolveOptions {std::size_t max_fill=12000000;double seconds=120;std::function<void(const Statistics&,double)> progress={};double progress_seconds=30;bool complete_rows=false;};
struct Pivot {Row row;std::uint32_t source;std::vector<Column> dependencies;};
class Solver {
 const Program& program_;const Field& field_;SolveOptions options_;std::chrono::steady_clock::time_point started_;
 std::vector<int> lookup_;std::vector<Pivot> pivots_;Statistics stats_;
 // If set, the equations: rows of terms and the values their forms take.
 const FormRows* rows_=nullptr;const std::vector<Word>* forms_=nullptr;
 double told_=0;
 void budget(){if(stats_.fill>options_.max_fill)throw std::length_error("pivot fill budget exceeded");const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count();
  if(elapsed>options_.seconds)throw std::runtime_error("elimination time budget exceeded");
  if(options_.progress&&elapsed-told_>=options_.progress_seconds){told_=elapsed;options_.progress(stats_,elapsed);}}
 // A row under reduction lives in a dense vector with one bit per occupied
 // column. Subtracting a pivot row then costs the length of the pivot row,
 // not that of the row it is subtracted from, which a merge of two sorted
 // rows pays at every step. Between reductions both are all zero.
 static constexpr Column none=std::numeric_limits<Column>::max();
 std::vector<Word> dense_;std::vector<std::uint64_t> occupied_;
 void scatter(const Row& row){
  if(dense_.size()!=lookup_.size()){dense_.assign(lookup_.size(),0);occupied_.assign((lookup_.size()+63)/64,0);}
  for(const auto& e:row)if(e.column>=lookup_.size())throw std::invalid_argument("reduction column out of range");
  for(const auto& e:row){dense_[e.column]=e.value;occupied_[e.column>>6]|=std::uint64_t(1)<<(e.column&63);}
 }
 // The first occupied column in [from,high].
 Column next(Column from,Column high)const {
  if(from>high)return none;std::size_t w=from>>6;std::uint64_t word=occupied_[w]&(~std::uint64_t(0)<<(from&63));
  for(;;){if(word){const Column c=static_cast<Column>((w<<6)+std::countr_zero(word));return c<=high?c:none;}if(++w>(high>>6))return none;word=occupied_[w];}
 }
 void vacate(Column column){occupied_[column>>6]&=~(std::uint64_t(1)<<(column&63));}
 // Everything occupied belongs to the row under reduction.
 void clear(Column low,Column high){for(std::size_t w=low>>6;w<=(high>>6);++w){auto word=occupied_[w];while(word){dense_[(w<<6)+std::countr_zero(word)]=0;word&=word-1;}occupied_[w]=0;}}
 // Pivot rows are monic: this leaves zero in the pivot's own column.
 void eliminate(Column column,Word factor,Column& high){eliminate(pivots_[lookup_[column]].row,factor,high);}
 void eliminate(const Row& pivot,Word factor,Column& high){
  for(const auto& e:pivot){auto& value=dense_[e.column];value=field_.sub(value,field_.mul(factor,e.value));occupied_[e.column>>6]|=std::uint64_t(1)<<(e.column&63);}
  high=std::max(high,pivot.back().column);++stats_.eliminations;
 }
 // The rows of the leading way, for the pivots a selection asks about.
 struct Leading {Row row;std::vector<Column> dependencies;bool known=false;};
 std::vector<Leading> leading_;std::vector<std::uint32_t> mark_;std::uint32_t generation_=0;
 using Fetch=std::function<Row(std::uint32_t)>;
 Word dimension_=0;
 Row source_row(std::uint32_t source)const{Row row;
  if(rows_){const auto& terms=rows_->at(source);row.reserve(terms.size());for(const auto& a:terms)if(const auto value=(*forms_)[a.form])row.push_back({a.column,value});return row;}
  for(const auto& a:program_.equations.at(source)){const auto value=field_.add(a.constant,field_.mul(dimension_,a.dimension));if(value)row.push_back({a.column,value});}return row;}
 std::size_t equation_count()const{return rows_?rows_->size():program_.equations.size();}
 // The leading row of pivot `wanted`: its equation with the columns before
 // the pivot's own eliminated by the leading rows of their pivots, which are
 // found first. What can stand before the pivot column is known from the
 // equation and the rows it meets, so that nothing is asked for while the
 // dense vector is in use.
 void lead(int wanted,const Fetch& fetch){
  if(leading_.size()!=pivots_.size())leading_.resize(pivots_.size());
  if(mark_.size()!=lookup_.size())mark_.assign(lookup_.size(),0);
  std::vector<int> stack{wanted};std::vector<Column> work;std::size_t steps=0;
  while(!stack.empty()){
   const int p=stack.back();if(leading_[p].known){stack.pop_back();continue;}
   if((steps++&63)==0)budget();
   const Column own=pivots_[p].row.front().column;Row row=fetch(pivots_[p].source);bool missing=false;
   if(++generation_==0){std::fill(mark_.begin(),mark_.end(),0);generation_=1;}
   work.clear();for(const auto& e:row)if(e.column<own)work.push_back(e.column);
   while(!work.empty()){const Column c=work.back();work.pop_back();if(mark_[c]==generation_)continue;mark_[c]=generation_;
    const int q=lookup_[c];if(q<0||q>=p)continue;
    if(!leading_[q].known){stack.push_back(q);missing=true;continue;}
    for(const auto& e:leading_[q].row)if(e.column<own&&mark_[e.column]!=generation_)work.push_back(e.column);
   }
   if(missing)continue;
   Leading out;
   if(!row.empty()&&row.front().column<own){
    const Column low=row.front().column;Column high=row.back().column,column=low;scatter(row);row.clear();
    try{
     for(;(column=next(column,high))!=none&&column<own;++column){
      const Word value=dense_[column];dense_[column]=0;vacate(column);if(!value)continue;
      const int q=lookup_[column];if(q<0||q>=p||!leading_[q].known)throw std::logic_error("leading row: a column before the pivot has no earlier pivot");
      out.dependencies.push_back(column);dense_[column]=value;eliminate(leading_[q].row,value,high);dense_[column]=0;vacate(column);
     }
     for(;column!=none&&(column=next(column,high))!=none;++column){if(dense_[column])row.push_back({column,dense_[column]});dense_[column]=0;vacate(column);}
    }catch(...){clear(low,high);throw;}
   }
   if(row.empty()||row.front().column!=own)throw std::logic_error("leading row does not end at its pivot column");
   const auto inverse=field_.inv(row.front().value);for(auto& entry:row)entry.value=field_.mul(entry.value,inverse);
   out.row=std::move(row);out.known=true;leading_[p]=std::move(out);stack.pop_back();
  }
 }
 Row subtract(const Row& a,const Row& b,Word factor){Row out;out.reserve(a.size()+b.size());std::size_t i=0,j=0;while(i<a.size()||j<b.size()){
  if(j==b.size()||(i<a.size()&&a[i].column<b[j].column))out.push_back(a[i++]);else if(i==a.size()||b[j].column<a[i].column){auto v=field_.sub(0,field_.mul(factor,b[j].value));if(v)out.push_back({b[j].column,v});++j;}else{auto v=field_.sub(a[i].value,field_.mul(factor,b[j].value));if(v)out.push_back({a[i].column,v});++i;++j;}}
  stats_.peak_row=std::max(stats_.peak_row,out.size());return out;
 }
 public:
 Solver(const Program& program,const Field& field,SolveOptions options={}):program_(program),field_(field),options_(options),lookup_(program.integrals.size(),-1){}
 // The equations are `rows` with the coefficients `forms`; both must outlive
 // the solver, and `forms` keep their values while it is used.
 Solver(const Program& program,const Field& field,SolveOptions options,const FormRows& rows,const std::vector<Word>& forms):program_(program),field_(field),options_(options),lookup_(program.integrals.size(),-1),rows_(&rows),forms_(&forms){}
 // `kept` continues an elimination: pivots of a subsystem of this program,
 // in this program's columns and equation indices. With them, `selected`
 // names the equations that subsystem did not contain.
 void solve(Word dimension,const std::vector<std::uint32_t>* selected=nullptr,std::vector<Pivot> kept={}){started_=std::chrono::steady_clock::now();told_=0;pivots_=std::move(kept);leading_.clear();dimension_=dimension;std::fill(lookup_.begin(),lookup_.end(),-1);stats_={};
  for(std::size_t i=0;i<pivots_.size();++i){const auto& row=pivots_[i].row;if(row.empty()||row.front().column>=lookup_.size()||lookup_[row.front().column]>=0)throw std::invalid_argument("kept pivots are not an echelon basis of this program");lookup_[row.front().column]=i;stats_.fill+=row.size();++stats_.rank;}
  auto insert=[&](std::uint32_t source){if((stats_.equations++&127)==0)budget();Row row=source_row(source);std::vector<Column> dependencies;
   if(options_.complete_rows){
    bool reducible=false;for(const auto& entry:row)if(lookup_[entry.column]>=0){reducible=true;break;}
    if(reducible)row=reduce(std::move(row));
   }else if(!row.empty()&&lookup_[row.front().column]>=0){
    const Column low=row.front().column;Column high=row.back().column,column=low;scatter(row);row.clear();
    try{
     for(;(column=next(column,high))!=none;++column){
      const Word value=dense_[column];if(!value){vacate(column);continue;}
      if(lookup_[column]<0)break;
      dependencies.push_back(column);eliminate(column,value,high);vacate(column);if((stats_.eliminations&2047)==0)budget();
     }
     for(;column!=none&&(column=next(column,high))!=none;++column){if(dense_[column])row.push_back({column,dense_[column]});dense_[column]=0;vacate(column);}
    }catch(...){clear(low,high);throw;}
    stats_.peak_row=std::max(stats_.peak_row,row.size());
   }
   if(row.empty())return;auto inverse=field_.inv(row.front().value);for(auto& entry:row)entry.value=field_.mul(entry.value,inverse);auto column=row.front().column;row.shrink_to_fit();dependencies.shrink_to_fit();stats_.fill+=row.size();lookup_[column]=pivots_.size();pivots_.push_back({std::move(row),source,std::move(dependencies)});++stats_.rank;budget();
  };
  if(selected)for(auto row:*selected){if(row>=equation_count())throw std::invalid_argument("selected equation out of range");insert(row);}else for(std::uint32_t row=0;row<equation_count();++row)insert(row);
  stats_.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count();
 }
 Row reduce(Row row,std::vector<Column>* touched=nullptr){
  if(row.empty())return row;
  Column low=row.front().column,high=low;for(const auto& e:row){low=std::min(low,e.column);high=std::max(high,e.column);}
  scatter(row);row.clear();std::size_t steps=0;
  try{
   for(Column column=low;(column=next(column,high))!=none;++column){
    if((steps++&255)==0)budget();
    const Word value=dense_[column];
    if(value&&lookup_[column]>=0){if(touched)touched->push_back(column);eliminate(column,value,high);if((stats_.eliminations&2047)==0)budget();}
    else if(value)row.push_back({column,value});
    dense_[column]=0;vacate(column);
   }
  }catch(...){clear(low,high);throw;}
  stats_.peak_row=std::max(stats_.peak_row,row.size());return row;
 }
 std::vector<Row> targets(std::vector<Column>* touched=nullptr){return targets(program_.target_columns,touched);}
 std::vector<std::uint32_t> selection(){return selection(program_.target_columns);}
 // The same for columns other than the program's targets, and what an
 // elimination that is kept for later needs: a clock of its own for each use
 // and its pivots.
 std::vector<Row> targets(const std::vector<Column>& columns,std::vector<Column>* touched=nullptr){std::vector<Row> out;for(auto column:columns)out.push_back(reduce({{column,1}},touched));return out;}
 void restart(double seconds){options_.seconds=seconds;started_=std::chrono::steady_clock::now();told_=0;}
 std::vector<Pivot> release(){std::fill(lookup_.begin(),lookup_.end(),-1);stats_={};leading_.clear();return std::move(pivots_);}
 bool complete_rows()const{return options_.complete_rows;}
 // The equations that reproduce the reductions of the given columns when
 // eliminated alone, in order, the leading way: those of the pivots the
 // reductions meet and of the pivots their equations were reduced by. After
 // a complete elimination the leading rows are made for these pivots alone;
 // `fetch` then returns an equation at the point of the elimination (by
 // default from the program, if its coefficients are still those).
 std::vector<std::uint32_t> selection(const std::vector<Column>& columns,const Fetch& fetch={}){
  std::vector<bool> needed(pivots_.size());std::vector<int> stack;
  auto need=[&](Column column){const auto p=lookup_[column];if(p>=0&&!needed[p]){needed[p]=true;stack.push_back(p);}};
  if(!options_.complete_rows){
   std::vector<Column> roots;targets(columns,&roots);for(auto root:roots)need(root);
   while(!stack.empty()){const auto p=stack.back();stack.pop_back();for(auto column:pivots_[p].dependencies)need(column);}
  }else{
   const Fetch from=fetch?fetch:Fetch([this](std::uint32_t source){return source_row(source);});
   // A second dense vector for the targets: leading rows are made, in
   // dense_, while a target is being reduced.
   std::vector<Word> values(lookup_.size(),0);std::vector<std::uint64_t> bits((lookup_.size()+63)/64,0);
   for(auto target:columns){
    if(target>=lookup_.size())throw std::invalid_argument("reduction column out of range");
    Row reduced;Column high=target;values[target]=1;bits[target>>6]|=std::uint64_t(1)<<(target&63);
    for(Column column=target;column<=high;++column){
     std::size_t w=column>>6;std::uint64_t word=bits[w]&(~std::uint64_t(0)<<(column&63));
     while(!word&&++w<=(high>>6))word=bits[w];
     if(!word)break;column=static_cast<Column>((w<<6)+std::countr_zero(word));if(column>high)break;
     const Word value=values[column];values[column]=0;bits[column>>6]&=~(std::uint64_t(1)<<(column&63));if(!value)continue;
     const int p=lookup_[column];if(p<0){reduced.push_back({column,value});continue;}
     lead(p,from);need(column);
     for(const auto& e:leading_[p].row)if(e.column!=column){auto& v=values[e.column];v=field_.sub(v,field_.mul(value,e.value));bits[e.column>>6]|=std::uint64_t(1)<<(e.column&63);}
     high=std::max(high,leading_[p].row.back().column);
    }
    // Both ways reduce a target to the same combination of free columns.
    if(reduced!=reduce({{target,1}}))throw std::logic_error("the leading rows of a complete elimination reduce a target differently");
   }
   while(!stack.empty()){const auto p=stack.back();stack.pop_back();for(auto column:leading_[p].dependencies)need(column);}
  }
  std::vector<std::uint32_t> result;for(std::size_t p=0;p<pivots_.size();++p)if(needed[p])result.push_back(pivots_[p].source);std::sort(result.begin(),result.end());return result;
 }
 const Statistics& statistics()const{return stats_;}
 std::vector<Column> free_columns()const {std::vector<Column> result;for(unsigned i=0;i<lookup_.size();++i)if(lookup_[i]<0)result.push_back(i);return result;}
};
}
