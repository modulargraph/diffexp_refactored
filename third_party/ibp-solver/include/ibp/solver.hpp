#pragma once
#include "ibp/geometry.hpp"
#include <unordered_map>
#include <chrono>
#include <unordered_set>
#include <functional>
namespace ibp {
using Column=std::uint32_t;
struct LinearEntry {Column column;Word constant,dimension;};
struct Entry {Column column;Word value;bool operator==(const Entry&)const=default;};
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
struct SolveOptions {std::size_t max_fill=12000000;double seconds=120;std::function<void(const Statistics&,double)> progress={};double progress_seconds=30;};
struct Pivot {Row row;std::uint32_t source;std::vector<Column> dependencies;};
class Solver {
 const Program& program_;const Field& field_;SolveOptions options_;std::chrono::steady_clock::time_point started_;
 std::vector<int> lookup_;std::vector<Pivot> pivots_;Statistics stats_;
 double told_=0;
 void budget(){if(stats_.fill>options_.max_fill)throw std::length_error("pivot fill budget exceeded");const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count();
  if(elapsed>options_.seconds)throw std::runtime_error("elimination time budget exceeded");
  if(options_.progress&&elapsed-told_>=options_.progress_seconds){told_=elapsed;options_.progress(stats_,elapsed);}}
 Row subtract(const Row& a,const Row& b,Word factor){Row out;out.reserve(a.size()+b.size());std::size_t i=0,j=0;while(i<a.size()||j<b.size()){
  if(j==b.size()||(i<a.size()&&a[i].column<b[j].column))out.push_back(a[i++]);else if(i==a.size()||b[j].column<a[i].column){auto v=field_.sub(0,field_.mul(factor,b[j].value));if(v)out.push_back({b[j].column,v});++j;}else{auto v=field_.sub(a[i].value,field_.mul(factor,b[j].value));if(v)out.push_back({a[i].column,v});++i;++j;}}
  stats_.peak_row=std::max(stats_.peak_row,out.size());return out;
 }
 public:
 Solver(const Program& program,const Field& field,SolveOptions options={}):program_(program),field_(field),options_(options),lookup_(program.integrals.size(),-1){}
 // `kept` continues an elimination: pivots of a subsystem of this program,
 // in this program's columns and equation indices. With them, `selected`
 // names the equations that subsystem did not contain.
 void solve(Word dimension,const std::vector<std::uint32_t>* selected=nullptr,std::vector<Pivot> kept={}){started_=std::chrono::steady_clock::now();told_=0;pivots_=std::move(kept);std::fill(lookup_.begin(),lookup_.end(),-1);stats_={};
  for(std::size_t i=0;i<pivots_.size();++i){const auto& row=pivots_[i].row;if(row.empty()||row.front().column>=lookup_.size()||lookup_[row.front().column]>=0)throw std::invalid_argument("kept pivots are not an echelon basis of this program");lookup_[row.front().column]=i;stats_.fill+=row.size();++stats_.rank;}
  auto insert=[&](std::uint32_t source){if((stats_.equations++&127)==0)budget();Row row;for(const auto& a:program_.equations[source]){auto value=field_.add(a.constant,field_.mul(dimension,a.dimension));if(value)row.push_back({a.column,value});}std::vector<Column> dependencies;
   while(!row.empty()&&lookup_[row.front().column]>=0){auto column=row.front().column;dependencies.push_back(column);row=subtract(row,pivots_[lookup_[column]].row,row.front().value);++stats_.eliminations;if((stats_.eliminations&2047)==0)budget();}
   if(row.empty())return;auto inverse=field_.inv(row.front().value);for(auto& entry:row)entry.value=field_.mul(entry.value,inverse);auto column=row.front().column;stats_.fill+=row.size();lookup_[column]=pivots_.size();pivots_.push_back({std::move(row),source,std::move(dependencies)});++stats_.rank;budget();
  };
  if(selected)for(auto row:*selected){if(row>=program_.equations.size())throw std::invalid_argument("selected equation out of range");insert(row);}else for(std::uint32_t row=0;row<program_.equations.size();++row)insert(row);
  stats_.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count();
 }
 Row reduce(Row row,std::vector<Column>* touched=nullptr){std::size_t i=0,steps=0;while(i<row.size()){if((steps++&255)==0)budget();if(row[i].column>=lookup_.size())throw std::invalid_argument("reduction column out of range");auto column=row[i].column;if(lookup_[column]<0){++i;continue;}if(touched)touched->push_back(column);auto factor=row[i].value;row=subtract(row,pivots_[lookup_[column]].row,factor);}return row;}
 std::vector<Row> targets(std::vector<Column>* touched=nullptr){return targets(program_.target_columns,touched);}
 std::vector<std::uint32_t> selection(){return selection(program_.target_columns);}
 // The same for columns other than the program's targets, and what an
 // elimination that is kept for later needs: a clock of its own for each use
 // and its pivots.
 std::vector<Row> targets(const std::vector<Column>& columns,std::vector<Column>* touched=nullptr){std::vector<Row> out;for(auto column:columns)out.push_back(reduce({{column,1}},touched));return out;}
 void restart(double seconds){options_.seconds=seconds;started_=std::chrono::steady_clock::now();told_=0;}
 std::vector<Pivot> release(){std::fill(lookup_.begin(),lookup_.end(),-1);stats_={};return std::move(pivots_);}
 std::vector<std::uint32_t> selection(const std::vector<Column>& columns){std::vector<Column> roots;targets(columns,&roots);std::vector<bool> needed(pivots_.size());auto visit=[&](auto&& self,Column column)->void{auto p=lookup_[column];if(p<0||needed[p])return;needed[p]=true;for(auto dependency:pivots_[p].dependencies)self(self,dependency);};for(auto root:roots)visit(visit,root);std::vector<std::uint32_t> result;for(unsigned p=0;p<pivots_.size();++p)if(needed[p])result.push_back(pivots_[p].source);std::sort(result.begin(),result.end());return result;}
 const Statistics& statistics()const{return stats_;}
 std::vector<Column> free_columns()const {std::vector<Column> result;for(unsigned i=0;i<lookup_.size();++i)if(lookup_[i]<0)result.push_back(i);return result;}
};
}
