#pragma once
#include "diffexp/canonical.hpp"
#include "diffexp/geometry.hpp"
#include <functional>

namespace diffexp {
inline Jet::Ball multiple_polylogarithm(const std::vector<std::string>& word,
    const Rational& endpoint,unsigned order=72) {
  using B=Jet::Ball;
  if(endpoint<=Rational(0))throw std::invalid_argument("this MPL adapter requires a positive real endpoint");
  if(word.size()>100 || order<2 || order>500)throw std::invalid_argument("MPL resource limits");
  std::map<std::vector<std::string>,B> cache;
  B endpoint_ball=B::from_strings(endpoint.str()),logarithm;
  acb_log(logarithm.raw(),endpoint_ball.raw(),B::precision());
  std::function<B(const std::vector<std::string>&)> solve=[&](const std::vector<std::string>& indices)->B {
    if(indices.empty())return B(1);
    if(auto found=cache.find(indices);found!=cache.end())return found->second;
    ExactField field({"t","I"});Exact t(field,"t"),imag(field,"I");
    unsigned zeros=0;
    for(auto i=indices.rbegin();i!=indices.rend() && Exact(field,*i).is_zero();++i)++zeros;
    B value;
    if(zeros==indices.size()) {
      value=B(1);for(unsigned i=1;i<=zeros;++i)value=value*logarithm/B(i);
    } else if(zeros) {
      auto prefix=indices;prefix.pop_back();value=solve(prefix)*logarithm;
      // Shuffle with G(0;z). Exactly `zeros` insertions reproduce the original
      // word; every other word has fewer trailing zeros.
      for(unsigned position=0;position<indices.size()-zeros;++position) {
        auto shuffled=prefix;shuffled.insert(shuffled.begin()+position,"0");value-=solve(shuffled);
      }
      value=value/B(zeros);
    } else {
      unsigned weight=indices.size(),dimension=weight+1;
      auto z=t*t.constant(endpoint)-imag*t.constant(Rational("2/5"))*t*(t.constant(1)-t);
      std::vector<B> roots;
      for(const auto& letter:indices) {
        auto candidate=polynomial_norm(z-Exact(field,letter),1,t.constant(-1));
        auto found=polynomial_roots(candidate,0,B::precision());roots.insert(roots.end(),found.begin(),found.end());
      }
      std::vector<CanonicalEntry> entries;
      for(unsigned i=0;i<weight;++i)entries.push_back({i,i+1,i,Rational(1)});
      Boundary boundary(dimension,std::vector<B>(weight+1,B(0)));boundary.back()[0]=B(1);
      double center=0;unsigned charts=0;
      while(center<1) {
        if(++charts>20000)throw std::runtime_error("MPL chart budget exhausted");
        auto finite_roots=roots;
        if(center==0)finite_roots.erase(std::remove_if(finite_roots.begin(),finite_roots.end(),[](const B& r){return r.contains_zero();}),finite_roots.end());
        double next=clearance_endpoint(center,finite_roots);
        B c,end;acb_set_d(c.raw(),center);acb_set_d(end.raw(),next);
        Jet tau(0,order+4,B::precision());tau.set(0,c);tau.set(1,B(1));
        auto e=tau.constant(0),im=tau.constant(0);e.set(0,endpoint_ball);im.set(0,B::from_strings("0","-2/5"));
        auto point=e*tau+im*tau*(tau.constant(1)-tau);
        std::vector<Jet> letters;
        for(const auto& index:indices)letters.push_back(point-evaluate(data::Reader(index).read(),tau,{}));
        boundary=canonical_chart(entries,letters,boundary,B(0),end-c,order);center=next;
      }
      value=boundary[0][weight];
    }
    cache.emplace(indices,value);return value;
  };
  return solve(word);
}

inline int run_mpl(bool weight20=true) {
  using B=Jet::Ball;B::set_precision(384);
  struct Case {std::string name;std::vector<std::string> word;Rational endpoint;std::string real,imag;double tolerance;};
  std::vector<Case> cases{
    {"G[1,0,1;4]",{"1","0","1"},Rational(4),"-6.7782180257804207212554826775005988168291802221955692129682","0.9250147943833369547396749852220309435917997631163983727603",1e-34},
    {"G[1,-10,0;4]",{"1","-10","0"},Rational(4),"-0.0191508840720296721365611597236750922866172732200324064383","-0.3066358899483403657463434439014286049874538907865239005438",1e-34},
    {"G[10,-10+I,-1/2,-50;1]",{"10","-10+I","-1/2","-50"},Rational(1),"-0.0000098802442781507281548895360764863423574760704710022738","-0.0000009352314628872620198852585457725779647560856726857223",1e-34}};
  if(weight20) {
    std::vector<std::string> indices;for(unsigned i=1;i<=20;++i)indices.push_back(std::to_string(i));
    cases.push_back({"weight20",indices,Rational(21),"0.00000000000513066731719907533179589918813462949766948546641803107466616580097","0",1e-40});
  }
  bool passed=true;
  for(const auto& test:cases) {
    std::cout<<"MPL "<<test.name<<": starting"<<std::endl;
    auto start=std::chrono::steady_clock::now();
    auto value=multiple_polylogarithm(test.word,test.endpoint);
    auto difference=value-B::from_strings(test.real,test.imag);arf_t upper;arf_init(upper);
    acb_get_abs_ubound_arf(upper,difference.raw(),B::precision());auto error=arf_get_d(upper,ARF_RND_CEIL);arf_clear(upper);
    passed=passed && std::isfinite(error) && error<test.tolerance;
    std::cout<<"MPL "<<test.name<<": discrepancy = "<<error<<", seconds = "
      <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
  }
  std::cout<<"MPL numerical reference comparison; Taylor remainders not yet certified.\n";
  return passed?0:1;
}
} // namespace diffexp
