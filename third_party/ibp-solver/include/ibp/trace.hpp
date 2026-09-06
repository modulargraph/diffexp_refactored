#pragma once
#include "ibp/parametric.hpp"
#include <limits>

namespace ibp {
// Traces are prime-specific. Dimension traces fix the geometry; parametric
// traces accept all changing contraction coefficients as explicit inputs.
class ArithmeticTrace {
 using Id=std::uint32_t;
 enum class Op:std::uint8_t {constant,dimension,add,sub,mul,inv,input};
 struct Node {Op op;Id a=0,b=0;Word value=0;};
 struct Key {Op op;Id a,b;bool operator==(const Key&)const=default;};
 struct Hash {std::size_t operator()(Key k)const{return (static_cast<std::size_t>(k.a)*0x9e3779b1u)^(static_cast<std::size_t>(k.b)*0x85ebca6bu)^static_cast<unsigned>(k.op);}};
 struct TraceEntry {Column column;Id id;};
 using TraceRow=std::vector<TraceEntry>;
 struct Instruction {Op op;Id dest,a,b;Word constant;};
 Word prime_=0;std::vector<Instruction> code_;std::vector<Id> guards_;
 std::vector<std::vector<std::pair<Column,Id>>> outputs_;
 std::size_t nodes_=0,registers_=0,inversions_=0,input_count_=0;
 std::vector<std::pair<Id,Word>> constant_inputs_;
 public:
 struct Options {std::size_t max_nodes=4000000,max_fill=12000000;double seconds=120;bool cse=true;};
 struct Workspace {std::vector<Word> values;};
 struct DeviceProgram {
  Word prime;std::uint32_t registers;
  std::vector<std::array<std::uint32_t,6>> instructions;
  std::vector<std::uint32_t> guards;
  std::vector<std::vector<std::pair<Column,std::uint32_t>>> outputs;
 };
 DeviceProgram device_program()const {
  if(input_count_)throw std::invalid_argument("device replay currently supports dimension traces only");
  DeviceProgram out{prime_,static_cast<std::uint32_t>(registers_),{},guards_,outputs_};
  for(auto i:code_)out.instructions.push_back({static_cast<std::uint32_t>(i.op),i.dest,i.a,i.b,static_cast<std::uint32_t>(i.constant),static_cast<std::uint32_t>(i.constant>>32)});
  return out;
 }
 private:
 static ArithmeticTrace build(const Program& program,const Field& f,Word dimension,
                              const std::vector<std::uint32_t>& selection,Options options,
                              const ParametricProgram* parametric=nullptr,const std::vector<Word>& inputs={}){
  const auto start=std::chrono::steady_clock::now();
  auto budget=[&]{if(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()>options.seconds)throw std::runtime_error("trace learning time budget exceeded");};
  std::vector<Node> nodes;std::unordered_map<Word,Id> constants;std::unordered_map<Key,Id,Hash> expressions;
  auto append=[&](Node n)->Id{if(nodes.size()>=options.max_nodes)throw std::length_error("trace node budget exceeded");if((nodes.size()&16383)==0)budget();nodes.push_back(n);return nodes.size()-1;};
  auto constant=[&](Word v)->Id{auto it=constants.find(v);if(it!=constants.end())return it->second;auto id=append({Op::constant,0,0,v});constants.emplace(v,id);return id;};
  const Id zero=constant(0),one=constant(1),dim=append({Op::dimension,0,0,dimension});
  auto operation=[&](Op op,Id a,Id b=0)->Id{
   if(op==Op::add){if(a==zero)return b;if(b==zero)return a;}
   if(op==Op::sub){if(b==zero)return a;if(a==b)return zero;}
   if(op==Op::mul){if(a==zero||b==zero)return zero;if(a==one)return b;if(b==one)return a;}
   if((op==Op::add||op==Op::mul)&&b<a)std::swap(a,b);
   Word v=0;switch(op){case Op::add:v=f.add(nodes[a].value,nodes[b].value);break;case Op::sub:v=f.sub(nodes[a].value,nodes[b].value);break;case Op::mul:v=f.mul(nodes[a].value,nodes[b].value);break;case Op::inv:v=f.inv(nodes[a].value);break;default:throw std::logic_error("invalid trace operation");}
   if(nodes[a].op==Op::constant&&(op==Op::inv||nodes[b].op==Op::constant))return constant(v);
   Key key{op,a,b};if(options.cse){auto it=expressions.find(key);if(it!=expressions.end())return it->second;}
   auto id=append({op,a,b,v});if(options.cse)expressions.emplace(key,id);return id;
  };
  std::vector<Id> input_ids;
  if(parametric)for(Id i=0;i<inputs.size();++i)input_ids.push_back(parametric->constant_inputs[i]?constant(inputs[i]):append({Op::input,i,0,inputs[i]}));
  std::vector<Id> zero_guards;
  auto nonzero=[&](Id id){if(nodes[id].value)return true;if(nodes[id].op!=Op::constant)zero_guards.push_back(id);return false;};
  auto subtract=[&](const TraceRow& a,const TraceRow& b,Id factor){TraceRow row;row.reserve(a.size()+b.size());std::size_t i=0,j=0;
   while(i<a.size()||j<b.size()){
    if(j==b.size()||(i<a.size()&&a[i].column<b[j].column))row.push_back(a[i++]);
    else if(i==a.size()||b[j].column<a[i].column){auto id=operation(Op::sub,zero,operation(Op::mul,factor,b[j].id));if(nonzero(id))row.push_back({b[j].column,id});++j;}
    else{auto id=operation(Op::sub,a[i].id,operation(Op::mul,factor,b[j].id));if(nonzero(id))row.push_back({a[i].column,id});++i;++j;}
   }return row;
  };
  std::vector<int> lookup(program.integrals.size(),-1);std::vector<TraceRow> pivots;std::size_t fill=0;
  for(auto source:selection){budget();if(source>=program.equations.size())throw std::invalid_argument("trace source out of range");TraceRow row;
   for(std::size_t j=0;j<program.equations[source].size();++j){auto e=program.equations[source][j];Id id;
    if(parametric){const auto& c=parametric->coefficients.at(source).at(j);id=constant(f.integer(c.constant));for(auto [slot,k]:c.terms)id=operation(Op::add,id,operation(Op::mul,constant(f.integer(k)),input_ids.at(slot)));}
    else id=operation(Op::add,constant(e.constant),operation(Op::mul,constant(e.dimension),dim));
    if(nonzero(id))row.push_back({e.column,id});}
   while(!row.empty()&&lookup[row.front().column]>=0)row=subtract(row,pivots[lookup[row.front().column]],row.front().id);
   if(row.empty())continue;auto inv=operation(Op::inv,row.front().id);for(auto& e:row)e.id=operation(Op::mul,e.id,inv);
   // Normalization is algebraically one, even when the builder cannot simplify
   // a*inv(a). Keeping it symbolic would impede structural cancellation.
   row.front().id=one;lookup[row.front().column]=pivots.size();fill+=row.size();if(fill>options.max_fill)throw std::length_error("trace pivot fill budget exceeded");pivots.push_back(std::move(row));
  }
  std::vector<TraceRow> results;for(auto column:program.target_columns){TraceRow row{{column,one}};std::size_t i=0;
   while(i<row.size()){budget();auto pivot=lookup[row[i].column];if(pivot<0){++i;continue;}row=subtract(row,pivots[pivot],row[i].id);}results.push_back(std::move(row));
  }
  // Retain outputs and zero guards. A sample-specific cancellation must never
  // silently become a symbolic identity. Nonzero pivot guards are inv checks.
  std::vector<bool> live(nodes.size());std::vector<Id> stack=zero_guards;for(const auto& row:results)for(auto e:row)stack.push_back(e.id);
  // Every learned inverse is a structural pivot guard, including pivots whose
  // numeric result would otherwise disappear during output dead-code removal.
  for(Id i=0;i<nodes.size();++i)if(nodes[i].op==Op::inv)stack.push_back(i);
  while(!stack.empty()){auto i=stack.back();stack.pop_back();if(live[i])continue;live[i]=true;auto op=nodes[i].op;if(op==Op::add||op==Op::sub||op==Op::mul){stack.push_back(nodes[i].a);stack.push_back(nodes[i].b);}else if(op==Op::inv)stack.push_back(nodes[i].a);}
  std::vector<std::size_t> uses(nodes.size());for(Id i=0;i<nodes.size();++i)if(live[i]){auto op=nodes[i].op;if(op==Op::add||op==Op::sub||op==Op::mul){++uses[nodes[i].a];++uses[nodes[i].b];}else if(op==Op::inv)++uses[nodes[i].a];}
  std::sort(zero_guards.begin(),zero_guards.end());zero_guards.erase(std::unique(zero_guards.begin(),zero_guards.end()),zero_guards.end());
  for(auto i:zero_guards)++uses[i];for(const auto& row:results)for(auto e:row)++uses[e.id];
  ArithmeticTrace trace;trace.prime_=f.prime();trace.nodes_=nodes.size();trace.input_count_=inputs.size();
  if(parametric)for(Id i=0;i<inputs.size();++i)if(parametric->constant_inputs[i])trace.constant_inputs_.emplace_back(i,inputs[i]);
  std::vector<Id> locations(nodes.size()),free;
  auto release=[&](Id id){if(!--uses[id])free.push_back(locations[id]);};
  for(Id i=0;i<nodes.size();++i)if(live[i]){auto n=nodes[i];Id a=n.op==Op::input?n.a:0,b=0;if(n.op==Op::add||n.op==Op::sub||n.op==Op::mul){a=locations[n.a];b=locations[n.b];release(n.a);release(n.b);}else if(n.op==Op::inv){a=locations[n.a];release(n.a);++trace.inversions_;}
   Id dest;if(free.empty())dest=trace.registers_++;else{dest=free.back();free.pop_back();}locations[i]=dest;trace.code_.push_back({n.op,dest,a,b,n.value});if(!uses[i])free.push_back(dest);
  }
  for(auto id:zero_guards)trace.guards_.push_back(locations[id]);for(const auto& row:results){trace.outputs_.emplace_back();for(auto e:row)trace.outputs_.back().emplace_back(e.column,locations[e.id]);}
  return trace;
 }
 public:
 static ArithmeticTrace learn(const Program& p,const Field& f,Word d,const std::vector<std::uint32_t>& s,Options o){return build(p,f,d,s,o);}
 static ArithmeticTrace learn(const Program& p,const Field& f,Word d,const std::vector<std::uint32_t>& s){return learn(p,f,d,s,Options{});}
 static ArithmeticTrace learn_parametric(const ParametricProgram& p,const Field& f,const std::vector<Word>& inputs,const std::vector<std::uint32_t>& s,Options o){
  if(!p.input_count||inputs.size()!=p.input_count||p.constant_inputs.size()!=p.input_count||p.coefficients.size()!=p.program.equations.size())throw std::invalid_argument("parametric trace input arity");
  for(auto v:inputs)if(v>=f.prime())throw std::invalid_argument("noncanonical parametric trace input");
  return build(p.program,f,0,s,o,&p,inputs);
 }
 static ArithmeticTrace learn_parametric(const ParametricProgram& p,const Field& f,const std::vector<Word>& inputs,const std::vector<std::uint32_t>& s){return learn_parametric(p,f,inputs,s,Options{});}
 std::vector<Row> evaluate_inputs(const std::vector<Word>& inputs,const Field& f,Workspace& workspace)const {
  if(!input_count_||inputs.size()!=input_count_)throw std::invalid_argument("parametric trace input arity");
  for(auto v:inputs)if(v>=f.prime())throw std::invalid_argument("noncanonical parametric trace input");
  for(auto [i,v]:constant_inputs_)if(inputs[i]!=v)throw std::invalid_argument("parametric trace constant input changed");
  return execute(0,inputs,f,workspace);
 }
 std::vector<Row> evaluate(Word dimension,const Field& f,Workspace& workspace)const{
  if(input_count_)throw std::invalid_argument("parametric trace requires coefficient inputs");return execute(dimension,{},f,workspace);
 }
 private:
 std::vector<Row> execute(Word dimension,const std::vector<Word>& inputs,const Field& f,Workspace& workspace)const{
  if(f.prime()!=prime_||dimension>=prime_)throw std::invalid_argument("trace prime or sample mismatch");workspace.values.resize(registers_);auto& v=workspace.values;
  for(auto ins:code_){Word value;switch(ins.op){case Op::constant:value=ins.constant;break;case Op::dimension:value=dimension;break;case Op::input:value=inputs.at(ins.a);break;case Op::add:value=f.add(v[ins.a],v[ins.b]);break;case Op::sub:value=f.sub(v[ins.a],v[ins.b]);break;case Op::mul:value=f.mul(v[ins.a],v[ins.b]);break;case Op::inv:value=f.inv(v[ins.a]);break;}v[ins.dest]=value;}
  for(auto i:guards_)if(v[i])throw std::domain_error("trace cancellation guard failed; relearn at a generic point");
  std::vector<Row> result;result.reserve(outputs_.size());for(const auto& row:outputs_){result.emplace_back();for(auto [column,id]:row)if(v[id])result.back().push_back({column,v[id]});}return result;
 }
 public:
 Word prime()const{return prime_;}
 std::size_t output_rows()const{return outputs_.size();}
 std::size_t output_terms()const{std::size_t n=0;for(const auto& row:outputs_)n+=row.size();return n;}
 std::size_t instructions()const{return code_.size();}std::size_t registers()const{return registers_;}std::size_t learned_nodes()const{return nodes_;}std::size_t inversions()const{return inversions_;}
};
}
