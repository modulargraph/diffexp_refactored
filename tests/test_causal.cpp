#include "diffexp/causal.hpp"
#include <iostream>
using namespace diffexp;
static void require(bool b){if(!b)throw std::runtime_error("causal assertion failed");}
template<class F>void rejects(F f){bool failed=false;try{f();}catch(const std::invalid_argument&){failed=true;}require(failed);}
int main(){try {
  rejects([]{causal::Prescription{}.validate(0);});
  auto henn=causal::current_example("henn_double_pentagon_x0",7,causal::henn_anchors(),true);
  require(henn.f_rim==-1 && henn.levels[0].x_detour_sign==1 && henn.assurance==causal::Assurance::SuppliedPrescription);
  for(unsigned i=1;i<7;++i)require(henn.levels[i].x_detour_sign==-1);
  rejects([&]{henn.validate(6);});
  rejects([]{causal::current_example("henn_double_pentagon_x0",7);});
  rejects([]{causal::current_example("henn_double_pentagon_x0",7,causal::henn_anchors(),false);});
  auto wrong=causal::henn_anchors();wrong[1]=wrong[0];rejects([&]{causal::current_example("henn_double_pentagon_x0",7,wrong,true);});
  require(causal::simple_root_orientation(-1,Rational(3),true)==-1);
  require(causal::simple_root_orientation(-1,Rational(-3),true)==1);
  rejects([]{causal::simple_root_orientation(-1,Rational(0),true);});
  rejects([]{causal::simple_root_orientation(-1,Rational(3),false);});
  for(const auto& name:feynman::example_names())if(name!="henn_double_pentagon_x0") {
    auto example=feynman::example_family(name);auto p=causal::current_example(name,example.physical_count-1);
    require(p.f_rim==-1 && p.assurance==causal::Assurance::EuclideanPhysicalDomain);
    std::cout<<name<<" exact physical U/F positivity PASS\n";
  }
  std::cout<<"causal prescription validation PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
