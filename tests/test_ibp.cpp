#include "diffexp/ibp.hpp"
#include "diffexp/families.hpp"
#include <iostream>

using namespace diffexp;
using namespace diffexp::ibp;
void check(bool value,const char* message) {if(!value) throw std::runtime_error(message);}
template<class F> void rejects(F f,const char* message) {
  bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,message);
}

void massive_tadpole() {
  ExactField field({"m2","d"});Exact m2(field,"m2"),d(field,"d");
  ScalarProducts space(1,{},m2);
  auto denominator=space.zero();denominator.constant=m2;denominator.linear[0]=-m2.constant(1);
  Generator generator(PropagatorBasis({space,{denominator}}),d);
  auto identities=generator.relations({3});
  check(identities.size()==1,"one vacuum IBP vector");
  check(identities[0].size()==2 && identities[0].at({3})==Exact(field,"d-6") &&
        identities[0].at({4})==Exact(field,"6*m2"),"massive tadpole IBP signs and divergence");
  ExactReducer reducer(m2);
  for(int a=1;a<=4;++a) for(auto& relation:generator.relations({a})) reducer.insert(std::move(relation));
  Relation input{{{5},m2.constant(1)}};
  auto reduction=reducer.reduce(input);
  check(reduction.remainder.size()==1 && reduction.remainder.at({1})==
    Exact(field,"(1-d/2)*(2-d/2)*(3-d/2)*(4-d/2)/(24*m2^4)"),"tadpole reduction agrees with independent gamma ratio");
  check(reducer.verify(input,reduction),"reduction witness reconstructs all original equations");
  reduction.remainder.begin()->second=m2.constant(0);
  check(!reducer.verify(input,reduction),"tampered reduction witness rejected");
  rejects([&]{ExactReducer limited(m2,0);limited.insert(identities[0]);},"finite equation budget");
  std::set<Integral> seeds;
  check(for_each_seed(1,2,{1,1,7},[&](const Integral& seed){seeds.insert(seed);})==7 && seeds.size()==7,
    "finite seed range includes physical numerators and no duplicate seeds");
  rejects([&]{for_each_seed(1,2,{1,1,6},[](const Integral&){});},"finite seed bound terminates enumeration");
}

void merged_bubble() {
  ExactField field({"x","d"});Exact x(field,"x"),d(field,"d");
  auto raw=quadratic_family(feynman::banana(1,{Rational(1),Rational(1)}),x);
  Generator generator(PropagatorBasis(merge(raw,0,1,x)),d);
  check(generator.basis().physical_count==1 && generator.basis().denominators.size()==2,
    "one merged physical propagator and one explicit ISP");
  ExactReducer reducer(x);
  for(int a=1;a<=4;++a) for(int isp=0;isp>=-2;--isp)
    for(auto& equation:generator.relations({a,isp})) reducer.insert(std::move(equation));
  // Completing the square gives I_2=Gamma(2-d/2)*F^(d/2-2), F=1+x-x^2.
  // The generated parameter derivative and IBP system must recover its exact
  // logarithmic derivative without supplying that analytic reduction to them.
  auto derivative=generator.derivative({2,0},0);
  add(derivative,Integral{2,0},-Exact(field,"(d/2-2)*(1-2*x)/(1+x-x^2)"));
  const auto result=reducer.reduce(derivative);
  check(result.remainder.empty(),"merged bubble differential equation from native IBP");
  check(reducer.verify(derivative,result),"parameter derivative exact IBP witness");
  auto closed=differential_system(generator,reducer,{{2,0}},0,x);
  check(closed.matrix.size()==1 && closed.matrix[0][0]==Exact(field,"(d/2-2)*(1-2*x)/(1+x-x^2)"),
    "requested bubble basis is closed without importing an analytic connection");
  rejects([&]{differential_system(generator,reducer,{{1,0},{2,0}},0,x);},"dependent requested masters rejected");
  ExactReducer empty(x);
  rejects([&]{differential_system(generator,empty,{{2,0}},0,x);},"incomplete reduction cannot claim derivative closure");
  rejects([&]{generator.derivative({2,0},1);},"dimensional differentiation is not a Feynman-parameter derivative");
  rejects([&]{generator.relations({2,1});},"ISP positive powers rejected");
  std::cout<<"Merged bubble: "<<reducer.equation_count()<<" exact identities, rank "<<reducer.rank()<<'\n';
}

void banana4_family() {
  ExactField field({"x","d"});Exact x(field,"x"),d(field,"d");
  auto raw=quadratic_family(feynman::banana(4,std::vector<Rational>(5,Rational(1))),x);
  auto merged=merge(raw,0,1,x);
  Generator generator(PropagatorBasis(raw),d);
  check(generator.basis().denominators.size()==14 && generator.basis().physical_count==5,
    "banana4 has 14 scalar products and nine auxiliary numerator slots");
  Integral target(14,0);std::fill_n(target.begin(),5,1);
  auto equations=generator.relations(target);
  check(equations.size()==20,"banana4 has 4 times 5 native IBP vector identities");
  check(generator.derivative(target,0).empty(),"unmerged family has no Feynman parameter dependence");
  const auto& basis=generator.basis();
  for(std::size_t i=0;i<basis.denominators.size();++i) {
    auto rewritten=basis.rewrite(basis.denominators[i]);
    check(rewritten.constant.is_zero(),"scalar product inverse cancels constants");
    for(std::size_t j=0;j<rewritten.linear.size();++j)
      check(rewritten.linear[j]==x.constant(i==j),"scalar product inverse reconstructs original basis");
  }
  Generator level(PropagatorBasis(std::move(merged)),d);
  target.assign(14,0);std::fill_n(target.begin(),4,1);target[0]=2;
  check(!level.derivative(target,0).empty(),"merged banana4 derivative is generated natively");
  raw.physical.push_back(raw.physical.front());
  rejects([&]{PropagatorBasis dependent(raw);},"dependent physical denominators require explicit partial fractions");
}

void example_registry() {
  ExactField field({"x","eps"});Exact x(field,"x");
  for(const auto& name:feynman::example_names()) {
    auto example=feynman::example_family(name);
    auto raw=quadratic_family(example.momenta,x,example.physical_count);
    Generator generator(PropagatorBasis(raw),x.constant(example.dimension_at_epsilon_zero)-x.constant(2)*x.variable(1));
    const auto& basis=generator.basis();
    check(basis.physical_count==example.physical_count,"native family physical-slot count");
    if(name=="henn_double_pentagon_x0") {
      check(basis.denominators.size()==11 && basis.physical_count==8,"Henn preserves eight physical and three named numerator slots");
      for(unsigned i=0;i<3;++i)
        check(basis.denominators[i+8].constant==raw.auxiliary[i].constant &&
              basis.denominators[i+8].linear==raw.auxiliary[i].linear,"Henn original numerator basis order");
    }
    // Every intermediate FT family must still have a complete scalar-product
    // basis. Use distinct exact interior anchors to avoid accidental diagonal
    // degeneracies. This checks preparation only, not numerical integration.
    while(raw.physical.size()>1) {
      raw=merge(raw,0,1,x.constant(Rational(static_cast<long>(raw.physical.size()))/Rational(17)));
      PropagatorBasis intermediate(raw);
      check(intermediate.denominators.size()==basis.space.size(),"complete scalar-product basis survives every merge");
    }
  }
  auto massive=feynman::example_family("pentagon_massive");
  const std::vector<Rational> masses{Rational(1),Rational("3/2"),Rational("4/3"),Rational("5/4"),Rational("6/5")};
  for(unsigned i=0;i<5;++i)check(massive.momenta.lines[i].mass_squared==masses[i],"massive pentagon reference masses");
  rejects([&]{feynman::example_family("missing");},"unknown native family rejected");
}

int main() {try {
  massive_tadpole();merged_bubble();banana4_family();example_registry();
  std::cout<<"Native quadratic compiler, IBP generator, exact reduction and witnesses passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
