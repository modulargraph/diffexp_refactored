#pragma once
#include "ibp/solver.hpp"
#include <limits>

namespace ibp {
// A coefficient is linear in supplied contraction values and dimension.
// Geometry may depend arbitrarily on parameters; integral/row identities do not.
struct InputGeometry {
 static constexpr std::uint32_t zero=std::numeric_limits<std::uint32_t>::max();
 unsigned loops=0,externals=0,physical=0,n=0;
 std::size_t inputs=0;std::uint32_t dimension=0;
 std::vector<bool> trace,zero_sectors,constant_inputs;
 std::vector<std::vector<std::vector<std::uint32_t>>> contractions;
};
struct InputCoefficient {
 std::int64_t constant=0;
 std::vector<std::pair<std::uint32_t,std::int64_t>> terms;
 Word evaluate(const std::vector<Word>& inputs,const Field& f)const {
  auto v=f.integer(constant);for(auto [id,k]:terms)v=f.add(v,f.mul(f.integer(k),inputs.at(id)));return v;
 }
 void normalize(){
  std::sort(terms.begin(),terms.end());std::size_t count=0;
  for(auto term:terms){if(count&&terms[count-1].first==term.first)terms[count-1].second+=term.second;else terms[count++]=term;}
  terms.resize(count);std::erase_if(terms,[](auto t){return !t.second;});
 }
 bool empty()const{return !constant&&terms.empty();}
};
enum class ParametricOrdering {dots_first,total_degree};
// Treat numerator powers and extra denominator powers at the same degree.
// At equal degree eliminate numerators first, retaining positive-index masters.
inline auto total_degree_grade(const Integral& a){auto [sector,dots,nums,powers]=grade(a);return std::tuple{sector,dots+nums,nums,powers};}
class ParametricProgram {
 public:
 Program program;
 std::vector<std::vector<InputCoefficient>> coefficients;
 std::vector<bool> constant_inputs;
 std::size_t input_count=0;
 void bind(const std::vector<Word>& inputs,const Field& f,const std::vector<std::uint32_t>* selection=nullptr){
  if(inputs.size()!=input_count)throw std::invalid_argument("parametric input arity");
  for(auto v:inputs)if(v>=f.prime())throw std::invalid_argument("noncanonical parametric input");
  auto row=[&](std::size_t i){if(i>=program.equations.size())throw std::invalid_argument("parametric source index");
   for(std::size_t j=0;j<program.equations[i].size();++j)program.equations[i][j].constant=coefficients[i][j].evaluate(inputs,f);};
  if(selection)for(auto i:*selection)row(i);else for(std::size_t i=0;i<program.equations.size();++i)row(i);
 }
 static ParametricProgram compile(const InputGeometry& g,const SeedOptions& options,std::vector<Integral> targets={},ParametricOrdering ordering=ParametricOrdering::dots_first){
  if(!g.loops||g.loops>4||!g.physical||g.physical>12||g.n>16||g.n<g.physical||!g.inputs||g.inputs>=InputGeometry::zero||g.dimension>=g.inputs||g.trace.size()!=g.contractions.size()||g.zero_sectors.size()!=(std::size_t(1)<<g.physical)||options.dots>8||options.numerators>8)
   throw std::invalid_argument("parametric geometry or seed bounds");
  if(!g.constant_inputs.empty()&&g.constant_inputs.size()!=g.inputs)throw std::invalid_argument("constant input mask arity");
  for(const auto& op:g.contractions){if(op.size()!=g.n)throw std::invalid_argument("contraction slot arity");for(const auto& row:op){if(row.size()!=g.n+1)throw std::invalid_argument("contraction affine arity");for(auto id:row)if(id!=InputGeometry::zero&&id>=g.inputs)throw std::invalid_argument("contraction input id");}}
  const auto start=std::chrono::steady_clock::now();std::size_t states=0,total_terms=0;
  ParametricProgram out;out.input_count=g.inputs;out.constant_inputs=g.constant_inputs;if(out.constant_inputs.empty())out.constant_inputs.resize(g.inputs);
  auto& p=out.program;std::unordered_map<Integral,Column,IntegralHash> ids;
  auto intern=[&](const Integral& a){auto [it,inserted]=ids.try_emplace(a,p.integrals.size());if(inserted)p.integrals.push_back(a);return it->second;};
  struct Term {Column column;InputCoefficient value;};using Equation=std::vector<Term>;std::vector<Equation> equations;
  if(targets.empty()){Integral top;for(unsigned k=0;k<g.physical;++k)top.powers[k]=1;targets.push_back(top);for(unsigned k=0;k<g.n;++k){auto a=top;a.powers[k]+=k<g.physical?1:-1;targets.push_back(a);}}
  p.targets=targets;
  for(const auto& target:targets){for(unsigned k=g.physical;k<16;++k)if(target.powers[k]>0||(k>=g.n&&target.powers[k]))throw std::invalid_argument("invalid target auxiliary slots");p.target_columns.push_back(intern(target));
   if(g.zero_sectors[sector(target,g.physical)])equations.push_back({{intern(target),{1,{}}}});}
  Integral seed;
  auto consume=[&]{if(g.zero_sectors[sector(seed,g.physical)])return;if(++p.seeds>options.max_seeds)throw std::length_error("parametric seed budget");
   for(unsigned op=0;op<g.contractions.size();++op){Equation raw;
    auto add=[&](const Integral& a,std::uint32_t input,std::int64_t scale){if(input!=InputGeometry::zero&&scale&&!g.zero_sectors[sector(a,g.physical)])raw.push_back({intern(a),{0,{{input,scale}}}});};
    if(g.trace[op])add(seed,g.dimension,1);
    for(unsigned k=0;k<g.n;++k)if(seed.powers[k]){auto raised=seed;++raised.powers[k];auto scale=-std::int64_t(seed.powers[k]);const auto& c=g.contractions[op][k];add(raised,c[g.n],scale);
     for(unsigned j=0;j<g.n;++j)if(c[j]!=InputGeometry::zero){auto lowered=raised;--lowered.powers[j];add(lowered,c[j],scale);}}
    std::sort(raw.begin(),raw.end(),[](const auto& a,const auto& b){return a.column<b.column;});Equation row;
    for(auto& term:raw){if(!row.empty()&&row.back().column==term.column){auto& terms=row.back().value.terms;terms.insert(terms.end(),term.value.terms.begin(),term.value.terms.end());}else row.push_back(std::move(term));}
    for(auto& term:row)term.value.normalize();std::erase_if(row,[](const auto& t){return t.value.empty();});
    if(!row.empty()){total_terms+=row.size();if(total_terms>options.max_terms)throw std::length_error("parametric term budget");equations.push_back(std::move(row));}
   }
  };
  auto visit=[&](auto&& self,unsigned k,unsigned dots,unsigned nums)->void {
   if(++states>options.max_states)throw std::length_error("parametric seed state budget");
   if((states&4095)==0&&std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()>options.seconds)throw std::runtime_error("parametric generation time budget");
   if(k==g.n){consume();return;}seed.powers[k]=0;self(self,k+1,dots,nums);
   for(unsigned s=1;s<=nums;++s){seed.powers[k]=-static_cast<int>(s);self(self,k+1,dots,nums-s);}
   if(k<g.physical)for(unsigned d=0;d<=dots;++d){seed.powers[k]=1+d;self(self,k+1,dots-d,nums);}
  };visit(visit,0,options.dots,options.numerators);
  std::vector<Column> permutation(p.integrals.size()),inverse(p.integrals.size());std::iota(permutation.begin(),permutation.end(),0);
  std::sort(permutation.begin(),permutation.end(),[&](auto a,auto b){return ordering==ParametricOrdering::total_degree?total_degree_grade(p.integrals[a])>total_degree_grade(p.integrals[b]):grade(p.integrals[a])>grade(p.integrals[b]);});auto old=p.integrals;
  for(unsigned i=0;i<permutation.size();++i){p.integrals[i]=old[permutation[i]];inverse[permutation[i]]=i;}
  for(auto& row:equations){for(auto& t:row)t.column=inverse[t.column];std::sort(row.begin(),row.end(),[](const auto& a,const auto& b){return a.column<b.column;});}
  for(auto& c:p.target_columns)c=inverse[c];
  std::stable_sort(equations.begin(),equations.end(),[](const auto& a,const auto& b){return a.front().column==b.front().column?a.size()<b.size():a.front().column>b.front().column;});
  for(auto& row:equations){p.equations.emplace_back();out.coefficients.emplace_back();for(auto& t:row){p.equations.back().push_back({t.column,0,0});out.coefficients.back().push_back(std::move(t.value));}}
  return out;
 }
};
}
