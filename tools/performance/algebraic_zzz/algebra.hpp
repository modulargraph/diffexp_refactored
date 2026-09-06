#pragma once
#include "diffexp/transport.hpp"
#include "diffexp/polynomial_transport.hpp"
#include <fstream>
#include <iostream>
using namespace diffexp;
using B=kernel::ComplexBall;
using Clock=std::chrono::steady_clock;
double seconds(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
std::map<unsigned,Exact> decompose(const transport::Compiled& c,Exact q){
 auto num=q.numerator(),den=q.denominator();
 for(unsigned r=c.squares.size();r-->0;){
  bool present=false;for(const auto& t:den.numerator_terms())if(t.powers[2+r])present=true;
  if(!present)continue;
  std::vector<Exact> subs;for(const auto& name:q.variables())subs.emplace_back(c.field,name);
  subs[2+r]=-subs[2+r];auto conjugate=den.substitute(subs);
  num=c.reduced(num*conjugate);den=c.reduced(den*conjugate);
  auto reduced=c.reduced(num/den);num=reduced.numerator();den=reduced.denominator();
 }
 for(const auto& t:den.numerator_terms())for(unsigned r=0;r<c.squares.size();++r)if(t.powers[2+r])throw std::runtime_error("unrationalized denominator");
 std::map<unsigned,Exact> out;
 for(const auto& t:num.numerator_terms()){
  unsigned mask=0;Exact term(c.field,t.coefficient.str());
  for(unsigned v=0;v<t.powers.size();++v){if(v<2)term=term*term.variable(v).pow(t.powers[v]);else{if(t.powers[v]>1)throw std::runtime_error("unreduced root power");if(t.powers[v])mask|=1u<<(v-2);}}
  auto [it,added]=out.try_emplace(mask,Exact(c.field,0));it->second=it->second+term/den;
 }
 Exact reconstructed(c.field,0);for(const auto& [m,a]:out){auto term=a;for(unsigned r=0;r<c.squares.size();++r)if(m&(1u<<r))term=term*term.variable(r+2);reconstructed=reconstructed+term;}
 if(!c.reduced(reconstructed-q).is_zero())throw std::runtime_error("exact algebraic reconstruction failed");
 return out;
}
