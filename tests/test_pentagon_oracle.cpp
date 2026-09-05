#include "diffexp/pentagon_oracle.hpp"
#include <iostream>
using namespace diffexp;using B=kernel::ComplexBall;
void require(bool b,const char* why){if(!b)throw std::runtime_error(why);}
template<class F>void rejects(F f){bool failed=false;try{f();}catch(const std::exception&){failed=true;}require(failed,"unsupported kinematics accepted");}
int main(){try {
  B::set_precision(256);B gamma,pi;arb_const_euler(acb_realref(gamma.raw()),256);acb_const_pi(pi.raw(),256);
  auto symmetric=oracle::one_mass_box_reference(Rational(-1),Rational(-1),Rational(-1),256);
  require(acb_equal_si(symmetric.at(-2).raw(),2),"one-mass leading normalization");
  require(acb_overlaps(symmetric.at(-1).raw(),(-B(2)*gamma).raw()),"rGamma epsilon coefficient");
  require(acb_overlaps(symmetric.at(0).raw(),(gamma*gamma-pi*pi/B(2)).raw()),"one-mass symmetric finite normalization");
  // Independent Li2(-1)=-pi^2/12 closed case probes the nonzero dilogs.
  B logtwo(2);acb_log(logtwo.raw(),logtwo.raw(),256);
  auto minusone=oracle::one_mass_box_reference(Rational(-1),Rational(-1),Rational(-2),256);
  require(acb_overlaps(minusone.at(-1).raw(),(B(2)*(logtwo-gamma)).raw()),"one-mass log-mass pole");
  require(acb_overlaps(minusone.at(0).raw(),(gamma*gamma-B(2)*gamma*logtwo-logtwo*logtwo-pi*pi/B(6)).raw()),"independent dilogarithm special value");
  auto a=oracle::one_mass_box_reference(Rational(-2),Rational(-3),Rational(-5),256),b=oracle::one_mass_box_reference(Rational(-3),Rational(-2),Rational(-5),256);
  for(int k=-2;k<=0;++k)require(acb_overlaps(a.at(k).raw(),b.at(k).raw()),"one-mass s/t symmetry");
  auto family=feynman::example_family("pentagon").momenta;auto result=oracle::massless_pentagon_reference(family,256);
  for(unsigned i=0;i<5;++i)std::cout<<"pinch "<<i<<" weight="<<result.box_weights[i].str()<<" s,t,m2="<<result.box_kinematics[i][0].str()<<","<<result.box_kinematics[i][1].str()<<","<<result.box_kinematics[i][2].str()<<'\n';
  // Independent exact native F identity: sum_i b_i*dF/da_i=sum_i a_i.
  // This fixes both the reduction sign and the factor of two in S.
  ExactField field({"a","b","c","d","e"});Exact zero(field);std::vector<Exact> parameters;
  for(unsigned i=0;i<5;++i)parameters.push_back(zero.variable(i));
  auto geometry=feynman::symanzik(family,parameters);Exact contraction=zero,sum=zero;
  const std::array<Rational,5> expected_weights{Rational("31/70"),Rational("-17/28"),Rational("7/4"),Rational("-3/20"),Rational("13/70")};
  for(unsigned i=0;i<5;++i){require(result.box_weights[i]==expected_weights[i],"exact fixture reduction weights");contraction=contraction+geometry.F.derivative(i)*zero.constant(result.box_weights[i]);sum=sum+parameters[i];}
  require(contraction==sum,"native Symanzik reduction identity");
  require(acb_contains(result.at(-2).raw(),B::from_strings("13/42").raw()),"independent exact leading pentagon residue");
  // Cyclic relabelling and reflection test pinched-box identification.
  auto rotated=family;std::rotate(rotated.lines.begin(),rotated.lines.begin()+1,rotated.lines.end());auto reversed=family;std::reverse(reversed.lines.begin(),reversed.lines.end());
  auto cyclic=oracle::massless_pentagon_reference(rotated,256),reflection=oracle::massless_pentagon_reference(reversed,256);
  for(int k=-2;k<=0;++k){require(acb_overlaps(result.at(k).raw(),cyclic.at(k).raw()),"cyclic invariance");require(acb_overlaps(result.at(k).raw(),reflection.at(k).raw()),"reflection invariance");}
  // J5(lambda*s)=lambda^(-3-eps)J5(s), including pole-induced log terms.
  auto doubled=family;for(auto& row:doubled.external_gram)for(auto& c:row)c=c*Rational(2);auto scaled=oracle::massless_pentagon_reference(doubled,256);B log2(2);acb_log(log2.raw(),log2.raw(),256);
  require(acb_overlaps((scaled.at(-2)*B(8)).raw(),result.at(-2).raw()),"pentagon mass dimension");
  require(acb_overlaps((scaled.at(-1)*B(8)).raw(),(result.at(-1)-log2*result.at(-2)).raw()),"pentagon residue scaling");
  require(acb_overlaps((scaled.at(0)*B(8)).raw(),(result.at(0)-log2*result.at(-1)+log2*log2*result.at(-2)/B(2)).raw()),"pentagon finite scaling");
  rejects([]{oracle::one_mass_box_reference(Rational(1),Rational(-1),Rational(-1));});
  auto bad=family;bad.lines[0].mass_squared=Rational(1);rejects([&]{oracle::massless_pentagon_reference(bad);});
  rejects([]{oracle::massless_pentagon_reference(feynman::example_family("box").momenta);});
  for(int k=-2;k<=0;++k){std::cout<<"pentagon eps^"<<k<<" ";acb_printn(result.at(k).raw(),45,0);std::cout<<'\n';}
  std::cout<<"native independent pentagon oracle PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
