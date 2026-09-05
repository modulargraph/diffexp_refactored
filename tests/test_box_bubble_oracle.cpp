#include "diffexp/box_bubble_oracle.hpp"
#include "diffexp/families.hpp"
#include <iostream>
using namespace diffexp;using B=kernel::ComplexBall;
void require(bool b,const char* why){if(!b)throw std::runtime_error(why);}
int main(){try {
 B::set_precision(256);auto ref=oracle::box_bubble_reference(Rational(-1),0,256);
 B gamma,z2,z3;arb_const_euler(acb_realref(gamma.raw()),256);arb_zeta_ui(acb_realref(z2.raw()),2,256);arb_zeta_ui(acb_realref(z3.raw()),3,256);
 auto l1=B(2)-B(2)*gamma,l2=B(2)-z2,l3=B(8)/B(3)-B(32)*z3/B(3);
 std::vector<B> independent{B(1)/B(2),l1/B(2),(l2+l1*l1/B(2))/B(2),(l3+l1*l2+l1*l1*l1/B(6))/B(2)};
 for(int k=-3;k<=0;++k)require(acb_overlaps(ref.at(k).raw(),independent[k+3].raw()),"independent log-gamma zeta coefficients");
 auto f=feynman::example_family("box_bubble").momenta;
 require(f.lines.size()==5 && f.lines[0].loop_coefficients==std::vector<Rational>{Rational(1),Rational(0)} && f.lines[1].loop_coefficients==std::vector<Rational>{Rational(1),Rational(-1)},"bubble loop routing");
 require(f.lines[0].external_coefficients==f.lines[3].external_coefficients,"bubble raises remaining p12 denominator");
 feynman::Family triangle{1,f.external_gram,{}};
 for(unsigned i=2;i<5;++i){auto line=f.lines[i];line.loop_coefficients={Rational(1)};triangle.lines.push_back(line);}
 ExactField field({"a","b","c"});std::vector<Exact> parameters{Exact(field,"a"),Exact(field,"b"),Exact(field,"c")};
 auto uf=feynman::symanzik(triangle,parameters);require(uf.F==Exact(field,"a*b"),"one-offshell triangle F and powered slot");
 auto scaled=oracle::box_bubble_reference(Rational(-2),0,256);B log2(2);acb_log(log2.raw(),log2.raw(),256);
 for(int k=-3;k<=0;++k){B expected(0),power(1);for(int n=0;n<=k+3;++n){expected+=power*ref.at(k-n);power=power*(-B(2)*log2)/B(n+1);}require(acb_overlaps((scaled.at(k)*B(2)).raw(),expected.raw()),"two-loop scale normalization");}
 for(int k=-3;k<=0;++k){std::cout<<"box_bubble eps^"<<k<<" ";acb_printn(ref.at(k).raw(),45,0);std::cout<<'\n';}
 std::cout<<"native box-bubble independent gamma oracle PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
