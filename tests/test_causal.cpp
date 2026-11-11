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
  // The automatic prescription.
  auto automatic=causal::derived(3);
  require(automatic.f_rim==-1 && automatic.derive_orientations && automatic.assurance==causal::Assurance::DerivedFromSymanzik && automatic.levels.size()==3);
  require(causal::Prescription{}.f_rim==-1);
  {
    // Leaf: F=1-9x(1-x) (equal-mass bubble at p^2=9m^2) has two simple zeros
    // of opposite slope; F - i0 passes the first above and the second below.
    ExactField field({"x","eps","I"});Exact x(field,"x");
    const auto zeros=causal::leaf_orientations(Exact(field,"1-9*x*(1-x)"),0);
    require(zeros.size()==2 && zeros[0].side==1 && zeros[1].side==-1 && zeros[0].high<Rational("1/2") && Rational("1/2")<zeros[1].low);
    require(causal::leaf_orientations(Exact(field,"1-9*x*(1-x)"),0,1)[0].side==-1);
    require(causal::leaf_orientations(Exact(field,"1-2*x"),0).size()==1 && causal::leaf_orientations(Exact(field,"1-2*x"),0)[0].side==1);
    require(causal::leaf_orientations(Exact(field,"1+x*(1-x)"),0).empty());
    bool repeated=false;try{causal::leaf_orientations(Exact(field,"(1-2*x)^2"),0);}catch(const std::domain_error&){repeated=true;}
    require(repeated);
    // One-loop level: the triangle D1=1-l^2, D2=1-(l+p)^2, D3=9-(l+q)^2 with
    // p^2=9, q^2=-1, p.q=0, first two lines merged. The two-particle cut of
    // the merged pair opens at x=(1 -/+ sqrt(5)/3)/2 with opposite slopes of
    // F, so no single side serves the level.
    feynman::Family triangle{1,{{Rational(9),Rational(0)},{Rational(0),Rational(-1)}},{
      {{Rational(1)},{Rational(0),Rational(0)},Rational(1)},{{Rational(1)},{Rational(1),Rational(0)},Rational(1)},
      {{Rational(1)},{Rational(0),Rational(1)},Rational(9)}}};
    const ibp::PropagatorBasis merged(ibp::merge(ibp::quadratic_family(triangle,x,3),0,1,x));
    const causal::Landau landau(merged);
    require(landau.loops()==1 && landau.propagators()==2);
    const auto signature=landau.signature(-1);
    std::cout<<"triangle level: "<<signature<<'\n';
    require(signature.find("0.12732200375003505060:1")!=std::string::npos && signature.find("0.87267799624996494940:-1")!=std::string::npos);
    require(landau.signature(1).find("0.12732200375003505060:-1")!=std::string::npos);
    const auto regular=landau.orient(Rational("1/2"),Rational("1/2"),-1);
    require(regular.kind==causal::Orientation::Kind::Free);
    const auto opened=landau.orient(Rational("12732200375003505/100000000000000000"),Rational("12732200375003506/100000000000000000"),-1);
    require(opened.kind==causal::Orientation::Kind::Oriented && opened.side==1);
    // The same triangle with spacelike p: nothing to orient.
    triangle.external_gram[0][0]=Rational(-9);
    const causal::Landau euclidean(ibp::PropagatorBasis(ibp::merge(ibp::quadratic_family(triangle,x,3),0,1,x)));
    require(euclidean.orient(Rational("1/8"),Rational("1/8"),-1).kind==causal::Orientation::Kind::Free);
    // Massless box, s=1, t=-1/3, first two lines merged:
    // F = a0 (-s x a1 - t (1-x) a2), whose x-derivative has one sign.
    feynman::Family box{1,{{Rational(0),Rational("1/2"),Rational("-1/3")},{Rational("1/2"),Rational(0),Rational("-1/6")},{Rational("-1/3"),Rational("-1/6"),Rational(0)}},{
      {{Rational(1)},{Rational(0),Rational(0),Rational(0)},Rational(0)},{{Rational(1)},{Rational(1),Rational(0),Rational(0)},Rational(0)},
      {{Rational(1)},{Rational(1),Rational(1),Rational(0)},Rational(0)},{{Rational(1)},{Rational(1),Rational(1),Rational(1)},Rational(0)}}};
    const causal::Landau massless(ibp::PropagatorBasis(ibp::merge(ibp::quadratic_family(box,x,4),0,1,x)));
    const auto adjacent=massless.orient(Rational("2/5"),Rational("2/5"),-1);
    std::cout<<"massless box, adjacent merge at 2/5: "<<adjacent.reason<<'\n';
    require(adjacent.kind!=causal::Orientation::Kind::Undetermined);
    // Two loops: the sunrise above its threshold has faces with vanishing U.
    // Nothing is guessed there; below it, in the Euclidean region, F is positive.
    auto sunrise=feynman::example_family("sunrise");sunrise.momenta.lines.resize(sunrise.physical_count);
    const causal::Landau two_loop_euclidean(ibp::PropagatorBasis(ibp::merge(ibp::quadratic_family(sunrise.momenta,x,3),0,1,x)));
    require(two_loop_euclidean.loops()==2 && two_loop_euclidean.orient(Rational("1/7"),Rational("1/7"),-1).kind==causal::Orientation::Kind::Free);
    sunrise.momenta.external_gram[0][0]=Rational(20);
    const causal::Landau two_loop(ibp::PropagatorBasis(ibp::merge(ibp::quadratic_family(sunrise.momenta,x,3),0,1,x)));
    const auto undecided=two_loop.orient(Rational("1/11"),Rational("1/11"),-1);
    require(undecided.kind==causal::Orientation::Kind::Undetermined && undecided.reason.find("one-loop")!=std::string::npos);
    std::cout<<"derived orientations: leaf zeros, one-loop faces, multi-loop refusal PASS\n";
  }
  std::cout<<"causal prescription validation PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
