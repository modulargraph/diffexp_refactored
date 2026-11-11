#include "diffexp/epsilon_sampling.hpp"
#include <iostream>
using namespace diffexp;
namespace es=epsilon_sampling;
using B=kernel::ComplexBall;
using Q=kernel::Rational;
static void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
static bool close(const B& value,const char* expected,slong bits) {
  const B difference=value-B::from_strings(expected);
  mag_t m;mag_init(m);acb_get_mag(m,difference.raw());const bool ok=mag_cmp_2exp_si(m,-bits)<0;mag_clear(m);return ok;
}
int main(){try {
  B::set_precision(1024);
  const auto x=es::nodes(Q("1/100"),10);
  check(x.size()==10 && x[0].str()=="1/100" && x[1].str()=="-1/100" && x[9].str()=="-1/20","symmetric nodes, innermost first");
  // f = 3/e^2 - 2/e + 5 + 7e - e^2 (+ i/e): a pole of order two and a
  // polynomial that ten nodes reproduce exactly.
  const auto f=[](const Q& q) {
    const B e=B::from_strings(q.str());
    return B(3)/(e*e)-B(2)/e+B(5)+B(7)*e-e*e+B::from_strings("0","1")/e;
  };
  std::vector<std::vector<B>> values(1);
  for(const auto& q:x)values[0].push_back(f(q));
  check(es::pole_order(x,values,4)==2,"pole order from the growth of the samples");
  const auto fitted=es::fit(x,values,2,2);
  check(fitted.poles==2 && fitted.coefficients.size()==1 && fitted.coefficients[0].size()==5,"fit shape");
  const char* expected[]{"3","-2","5","7","-1"};
  for(unsigned k=0;k<5;++k) {
    B real;acb_set_arb(real.raw(),acb_realref(fitted.coefficients[0][k].raw()));
    check(close(real,expected[k],600),"Laurent coefficient "+std::to_string(k));
    check(fitted.log2_change[0][k]<-600,"an exact polynomial does not change without its outermost samples");
  }
  B imaginary;acb_set_arb(imaginary.raw(),acb_imagref(fitted.coefficients[0][1].raw()));
  check(close(imaginary,"1",600),"imaginary part of the simple pole");
  // A function that is not a polynomial: the fit converges and its reported
  // change bounds the actual error here.
  std::vector<std::vector<B>> smooth(1);
  for(const auto& q:x){B e=B::from_strings(q.str()),value;acb_exp(value.raw(),e.raw(),B::precision());smooth[0].push_back(value/e);}
  const auto exponential=es::fit(x,smooth,1,1);
  check(close(exponential.coefficients[0][0],"1",40) && close(exponential.coefficients[0][1],"1",40) && close(exponential.coefficients[0][2],"1/2",30),"exp(e)/e");
  check(acb_contains(exponential.coefficients[0][2].raw(),B::from_strings("1/2").raw()),"the widened ball contains the coefficient");
  bool refused=false;
  try{(void)es::fit(es::nodes(Q("1/100"),4),{{B(1),B(1),B(1),B(1)}},2,1);}catch(const std::invalid_argument&){refused=true;}
  check(refused,"too few samples for the requested window are refused");
  refused=false;
  try{(void)es::nodes(Q("1/100"),5);}catch(const std::invalid_argument&){refused=true;}
  check(refused,"an odd number of samples is refused");
  std::cout<<"Epsilon sampling: nodes, pole order and Laurent fit\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
