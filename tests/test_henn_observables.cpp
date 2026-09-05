#include "diffexp/henn_observables.hpp"
#include "diffexp/jet.hpp"
#include <iostream>
using namespace diffexp;
void require(bool yes,const char* why){if(!yes)throw std::runtime_error(why);}
std::set<henn::Integral> source_targets(const data::Expr& e) {
  std::set<henn::Integral> result;
  if(e.head=="XB") {
    henn::Integral index;for(const auto& a:e.args)index.push_back(data::integer(a));result.insert(index);
  } else for(const auto& child:e.args) {
    auto targets=source_targets(child);result.insert(targets.begin(),targets.end());
  }
  return result;
}
// Independent reconstruction check: assign deterministic values to every
// scalar target, evaluate the original AST with the generic Acb evaluator,
// and compare with the imported linear combination for every component.
long scalar_weight(const henn::Integral& index) {
  long n=19;for(int i:index)n=(n*17+i+1024)%101;return n+1;
}
data::Expr scalar_assignment(data::Expr e) {
  if(e.head=="XB") {
    henn::Integral index;long sum=0;for(const auto& a:e.args){auto i=data::integer(a);index.push_back(i);sum+=i;}
    auto weight=scalar_weight(index);if(sum%2)weight=-weight;
    return data::Reader(std::to_string(weight)).read();
  }
  for(auto& child:e.args)child=scalar_assignment(std::move(child));return e;
}
int main(int argc,char** argv) {try {
  ExactField field({"eps"});henn::detail::Importer importer(field);
  const std::string xb="XB[1,1,0,1,1,0,1,0,0,0,0]";
  auto first=importer.observable(data::Reader("(s23-s45)*"+xb).read());
  require(first.size()==1 && first.begin()->second.rational.rational()==Rational(2) &&
      first.begin()->second.parity_odd.is_zero(),"first component X0 coefficient/sign must be +2");
  auto odd=importer.observable(data::Reader("(eps5+1/eps5)*(d-4)/eps*"+xb).read());
  require(odd.size()==1 && odd.begin()->second.rational.is_zero() &&
      odd.begin()->second.parity_odd.rational()==Rational("4/3"),"algebraic relation and d=4-2eps exact substitution");
  auto zero=importer.observable(data::Reader("(s23-s15)*"+xb).read());
  require(zero.empty() && importer.targets.size()==1,"source target retained when X0 coefficient vanishes");
  for(const auto& bad:std::vector<std::string>{xb+"*"+xb,xb+"^2","1/"+xb,"(0*"+xb+")*"+xb,
      "Sin[s12]*"+xb,"Sqrt[3]*"+xb,"XB[1,2]","XB[1,1,0,1,1,0,1,0,1,0,0]","1+"+xb,"0.5*"+xb,"unknown*"+xb}) {
    bool rejected=false;try {importer.observable(data::Reader(bad).read());}catch(const std::exception&){rejected=true;}
    require(rejected,"invalid or nonlinear observable accepted");
  }
  bool rejected=false;try{henn::import_x0(data::Reader("{"+xb+"}").read());}catch(const std::invalid_argument&){rejected=true;}
  require(rejected,"wrong basis dimension accepted");
  data::Expr fixture{"List",{}};for(unsigned i=0;i<108;++i)fixture.args.push_back(data::Reader("(s23-s45)*"+xb).read());
  auto small=henn::import_x0(fixture);require(small.components.size()==108 && small.scalar_targets.size()==1,"small basis fixture import");
  if(argc>1) {
    auto basis=henn::read_x0(argv[1]);auto source=data::read_file(argv[1]);auto targets=source_targets(source);
    require(basis.components.size()==108,"published basis has 108 components");
    require(basis.scalar_targets.size()==257,"published source has 257 distinct scalar targets");
    require(targets==std::set<henn::Integral>(basis.scalar_targets.begin(),basis.scalar_targets.end()),"all source targets preserved exactly");
    require(basis.components.front().size()==1 && basis.components.front().begin()->second.rational.rational()==Rational(2),"published component one coefficient");
    std::size_t numerator=0,odd_count=0;
    for(const auto& index:basis.scalar_targets) {
      require(index.size()==11 && index[8]<=0 && index[9]<=0 && index[10]<=0,"slot semantics preserved");
      if(*std::min_element(index.begin(),index.end())<0)++numerator;
    }
    for(const auto& observable:basis.components)for(const auto& [index,c]:observable) {
      require(!c.is_zero() && c.rational.variables()==std::vector<std::string>{"eps"},"coefficients canonical and exact X0 substitutions complete");
      if(!c.parity_odd.is_zero())++odd_count;
    }
    require(numerator==174 && odd_count>0,"published numerator and parity-odd terms preserved");
    Jet::Ball::set_precision(256);Jet context(0,1,256);
    auto root=evaluate(data::Reader("I*Sqrt[3]").read(),context,{});
    std::map<std::string,Jet> values{{"s12",context.constant(3)},{"s23",context.constant(-1)},
      {"s34",context.constant(1)},{"s45",context.constant(1)},{"s15",context.constant(-1)},
      {"eps5",root},{"eps",context.constant(1)/context.constant(7)}};
    for(std::size_t i=0;i<108;++i) {
      auto direct=evaluate(scalar_assignment(source.args[i]),context,values);
      auto imported=context.constant(0);
      for(const auto& [index,c]:basis.components[i])
        imported=imported+(evaluate(data::Reader(c.rational.str()).read(),context,values)+
          root*evaluate(data::Reader(c.parity_odd.str()).read(),context,values))*context.constant(scalar_weight(index));
      auto difference=(direct-imported).at(0);
      require(difference.contains_zero(),"full original-AST reconstruction disagrees with imported observable");
      mag_t bound;mag_init(bound);acb_get_mag(bound,difference.raw());
      const bool accurate=mag_get_d(bound)<1e-60;mag_clear(bound);
      require(accurate,"full original-AST reconstruction is not accurate");
    }
    std::cout<<"published components="<<basis.components.size()<<" source_targets="<<basis.scalar_targets.size()
        <<" nonzero_targets="<<basis.nonzero_targets.size()<<" numerator_targets="<<numerator<<" parity_odd_terms="<<odd_count<<"\n";
  }
  std::cout<<"Henn canonical observable import PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
