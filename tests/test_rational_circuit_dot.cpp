#include "diffexp/rational_circuit.hpp"
#include <iostream>
using namespace diffexp;namespace rc=rational_circuit;using B=Jet::Ball;
void check(bool b,const char* why){if(!b)throw std::runtime_error(why);}
void close(const B&a,const B&b,slong bits){B delta=a-b;mag_t error,tol;mag_init(error);mag_init(tol);acb_get_mag(error,delta.raw());mag_set_ui_2exp_si(tol,1,-bits/2);bool ok=mag_cmp(error,tol)<0;mag_clear(error);mag_clear(tol);check(ok,"retained value/width comparison failed");}
int main(){try{
 ExactField field({"x","eps","I"});auto e=[&](const char*s){return Exact(field,s);};
 for(slong bits:{128,320}){
  B::set_precision(bits);B h=B::from_strings("1/16","1/32"),center=B::from_strings("1/8","1/16");
  std::vector<RationalLineEntry> entries{{0,0,0,e("eps/(2-x+eps)")},{0,1,0,e("I/(1+x)")},
    {1,0,0,e("(1+x^7)/(2-x+eps)")},{1,1,1,e("(1+eps)/(2+x)")},
    {0,2,0,e("3+x^2")},{1,2,0,e("1/(1+x+eps)")},{3,3,0,e("1/(1-x)")}};
  Boundary input{{B(2),B(1),B(0),B(-1)},{B(3),B(0),B(1),B(0)},
    {B(1),B(0),B(0),B(0)},{B(0),B(0),B(0),B(0)}};
  for(bool uncertain:{false,true}){
   if(uncertain)arb_add_error_2exp_si(acb_realref(input[0][0].raw()),-bits+20);
   rc::Options scalar,dot;scalar.grouped_dot=false;dot.grouped_dot=true;
   auto plain=rc::prepare(rc::compile(entries,4,0,1,3,20,scalar),center,20);
   auto grouped=rc::prepare(rc::compile(entries,4,0,1,3,20,dot),center,20);
   rc::Statistics a,b;auto expected=rc::chart(plain,input,h,&a),got=rc::chart(grouped,input,h,&b);
   check(a.addmul_operations==b.addmul_operations,"dot changed executed term count");
   check(b.dot_products>0 && a.dot_products==0,"dot call count missing");
   check(b.live_cells-a.live_cells==2*grouped.data().dot_capacity+1,"dot scratch accounting mismatch");
   for(unsigned i=0;i<input.size();++i)for(unsigned k=0;k<input[0].size();++k){
    check(acb_overlaps(got[i][k].raw(),expected[i][k].raw()),"dot/scalar disjoint enclosures");close(got[i][k],expected[i][k],bits);
   }
   for(auto&v:got.back())check(v.is_zero(),"dot activated unreachable source");
  }
  rc::Options dot;dot.grouped_dot=true;
  // Independent finite coefficient oracle: source forcing (1+x^7)/(1-x).
  auto forcing=rc::compile({{0,1,0,e("(1+x^7)/(1-x)")}},2,0,1,0,20,dot);
  auto value=rc::chart(forcing,{{B(0)},{B(1)}},B(0),h,20);B oracle(0),power=h;
  for(unsigned n=1;n<=20;++n){oracle+=power*B(n>=8?2:1)/B(n);power*=h;}
  close(value[0][0],oracle,bits);
  // Signed epsilon window lower=-3: e^(eps*h)*eps^-3, no omitted
  // coefficients contribute to these four powers at spatial order >=3.
  auto shifted=rc::compile({{0,0,0,e("eps")}},1,0,1,3,8,dot);
  auto shiftvalue=rc::chart(shifted,{{B(1),B(0),B(0),B(0)}},B(0),h,8);B coefficient(1);
  for(unsigned k=0;k<4;++k){close(shiftvalue[0][k],coefficient,bits);coefficient=coefficient*h/B(k+1);}
  // Q has only epsilon lag: dot subtraction must subtract the whole sum
  // from Y, including higher powers generated at the same Taylor degree.
  auto eps=rc::compile({{0,1,0,e("1/(1-eps-eps^2)")}},2,0,1,5,3,dot);
  auto epsvalue=rc::chart(eps,{{B(0),B(0),B(0),B(0),B(0),B(0)},
    {B(1),B(0),B(0),B(0),B(0),B(0)}},B(0),h,3);
  int fib[]={1,1,2,3,5,8};for(unsigned k=0;k<6;++k)close(epsvalue[0][k],h*B(fib[k]),bits);
 }
 std::cout<<"Grouped quotient/numerator dots, scalar enclosures, ring alias/lags, Laurent offsets and independent finite oracles passed\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
