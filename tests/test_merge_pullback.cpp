#include "diffexp/merge_pullback.hpp"
#include "diffexp/families.hpp"
#include "diffexp/henn_observables.hpp"
#include <iostream>
#include <random>
using namespace diffexp;
void check(bool c,const char* message){if(!c)throw std::runtime_error(message);}
template<class F> void rejects(F f,const char* message){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,message);}
void numerator_identity(const ibp::PropagatorBasis& oldb,const ibp::PropagatorBasis& newb,
    const ibp::Integral& target,const pullback::Plan& plan,const Exact& x) {
  std::mt19937 random(1967);
  for(int trial=0;trial<4;++trial) {
    std::vector<Exact> products;for(std::size_t k=0;k<oldb.space.size();++k)products.push_back(x.constant(Rational(int(random()%19)-9)/Rational(7)));
    auto evaluate=[&](const ibp::Affine& d){auto value=d.constant;for(std::size_t k=0;k<products.size();++k)value=value+d.linear[k]*products[k];return value;};
    auto expected=x.constant(1);for(std::size_t k=0;k<target.size();++k)if(target[k]<0)expected=expected*evaluate(oldb.denominators[k]).pow(-target[k]);
    auto actual=x.constant(0);
    for(const auto& [source,c]:plan.source_row) {
      auto term=c;for(std::size_t k=0;k<source.size();++k) {
        auto power=plan.denominator_indices[k]-source[k];check(power>=0,"expanded numerator only lowers source indices");
        term=term*evaluate(newb.denominators[k]).pow(power);
      }
      actual=actual+term;
    }
    check(actual==expected,"random exact scalar-product numerator identity");
  }
}
int main(int argc,char** argv){try {
  ExactField field({"x"});Exact x(field,"x");
  auto raw=ibp::quadratic_family(feynman::banana(1,{Rational(1),Rational(2)}),x);
  ibp::PropagatorBasis oldb(raw),newb(ibp::merge(raw,0,1,x));
  auto beta=pullback::plan(oldb,newb,0,1,{2,3},x);
  check(beta.operation==feynman::Operation::BetaIntegral&&beta.normalization==Rational(12)&&beta.left_power==1&&beta.right_power==2,"normalized beta weights");
  check(beta.denominator_indices==ibp::Integral({5,0})&&beta.source_row.size()==1&&beta.source_row.at({5,0})==x.constant(1),"positive powers sum before numerator expansion");
  for(const auto& target:std::vector<ibp::Integral>{{-2,3},{3,-2},{-2,-3},{0,3},{3,0},{0,0}}) {
    auto plan=pullback::plan(oldb,newb,0,1,target,x);
    auto expected=target[0]>0?feynman::Operation::UpperLimit:target[1]>0?feynman::Operation::LowerLimit:feynman::Operation::Direct;
    check(plan.operation==expected,"typed endpoint or direct operation");
    check(plan.denominator_indices[0]==std::max(0,target[0])+std::max(0,target[1]),"negative merged powers stay numerator factors");
    numerator_identity(oldb,newb,target,plan,x);
  }
  auto endpoint=pullback::plan(oldb,newb,0,1,{3,-2},x);bool singular=false;
  for(const auto& [a,c]:endpoint.source_row)try{(void)c.substitute(std::vector<Exact>{x.constant(1)});}catch(const std::exception&){singular=true;}
  check(singular&&endpoint.operation==feynman::Operation::UpperLimit,"singular endpoint coefficient remains symbolic");
  auto sunrise=ibp::quadratic_family(feynman::banana(2,{Rational(1),Rational(1),Rational(1)}),x);
  ibp::PropagatorBasis so(sunrise);auto sm=ibp::merge(sunrise,0,1,x);ibp::PropagatorBasis sn(sm);
  auto mapped=pullback::plan(so,sn,0,1,{2,3,4,0,0},x);
  check(mapped.denominator_indices==ibp::Integral({5,4,0,0,0}),"nonmerged positive denominator slots preserved");
  sm.physical[1].constant=sm.physical[1].constant+x.constant(1);ibp::PropagatorBasis wrong(sm);
  rejects([&]{pullback::plan(so,wrong,0,1,{2,3,4,0,0},x);},"changed nonmerged propagator rejected");
  auto gram_changed=ibp::merge(raw,0,1,x);gram_changed.scalar_products.external_gram[0][0]=x.constant(-2);
  ibp::PropagatorBasis incompatible(gram_changed);
  rejects([&]{pullback::plan(oldb,incompatible,0,1,{1,1},x);},"incompatible external Gram rejected");
  // Replace all auxiliary coordinates genuinely, including affine constants.
  ibp::ScalarProducts space(1,{{x.constant(2),x.constant(1)},{x.constant(1),x.constant(3)}},x);
  auto d0=ibp::scaled(space.dot(0,0),x.constant(-1));d0.constant=x.constant(1);
  auto d1=ibp::plus(d0,ibp::scaled(space.dot(0,1),x.constant(2)));d1.constant=x.constant(4);
  auto numerator=ibp::plus(space.dot(0,2),ibp::scaled(space.dot(0,0),x.constant(2)));numerator.constant=x.constant(3);
  ibp::QuadraticFamily family{space,{d0,d1},{numerator}};
  ibp::PropagatorBasis old_changed(family);auto merged=ibp::merge(family,0,1,x);
  merged.auxiliary={ibp::plus(space.dot(0,1),space.dot(0,2)),space.dot(0,2)};
  ibp::PropagatorBasis new_changed(merged);
  for(const auto& target:std::vector<ibp::Integral>{{2,3,-2},{-2,3,-1},{3,-2,-1},{-1,-2,-2}}) {
    auto plan=pullback::plan(old_changed,new_changed,0,1,target,x);
    numerator_identity(old_changed,new_changed,target,plan,x);
  }
  rejects([&]{pullback::plan(old_changed,new_changed,0,1,{1,1,1},x);},"positive auxiliary rejected");
  rejects([&]{pullback::plan(old_changed,new_changed,0,1,{1,1,-3},x,{2,100,1000});},"numerator power budget");
  rejects([&]{pullback::plan(old_changed,new_changed,0,1,{1,1,-2},x,{3,1,1000});},"expanded term budget");
  rejects([&]{pullback::plan(oldb,newb,0,1,{INT_MIN,1},x);},"extreme negative power bounded safely");
  // A reversed merge has a different surviving slot and beta orientation.
  ibp::PropagatorBasis reverse(ibp::merge(raw,1,0,x));auto rp=pullback::plan(oldb,reverse,1,0,{3,2},x);
  check(rp.merged_slot==0&&rp.left_power==1&&rp.right_power==2,"reversed merge ordering");
  if(argc>1) {
    auto example=feynman::example_family("henn_double_pentagon_x0");
    auto hennraw=ibp::quadratic_family(example.momenta,x,example.physical_count);
    ibp::PropagatorBasis hb(hennraw),hm(ibp::merge(hennraw,0,1,x));
    auto imported=henn::read_x0(argv[1]);std::size_t numerators=0,terms=0;
    for(const auto& target:imported.scalar_targets) {
      auto plan=pullback::plan(hb,hm,0,1,target,x);terms+=plan.source_row.size();
      for(const auto& [source,c]:plan.source_row){check(!c.is_zero(),"zero expansion coefficients removed");for(std::size_t k=hm.physical_count;k<source.size();++k)check(source[k]<=0,"new Henn auxiliary slots remain numerators");}
      if(*std::min_element(target.begin(),target.end())<0){++numerators;numerator_identity(hb,hm,target,plan,x);}
    }
    check(imported.scalar_targets.size()==257&&numerators==174,"all published Henn numerator targets exercised");
    std::cout<<"Henn targets="<<imported.scalar_targets.size()<<" numerator targets="<<numerators<<" expanded terms="<<terms<<"\n";
  }
  std::cout<<"native merge pullback passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
