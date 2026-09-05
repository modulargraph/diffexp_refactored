#include "diffexp/gaussian_boundary.hpp"
#include "diffexp/families.hpp"
#include <iostream>
using namespace diffexp;
void check(bool c,const char* message){if(!c)throw std::runtime_error(message);}
template<class F>void rejects(F f,const char* message){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,message);}
int main(){try {
  Jet::Ball::set_precision(192);ExactField field({"x","d"});Exact x(field,"x"),d(field,"d");
  ibp::ScalarProducts sp(1,{},d);auto denominator=ibp::scaled(sp.dot(0,0),d.constant(-1));denominator.constant=d.constant(3);
  ibp::PropagatorBasis basis({sp,{denominator},{}});
  auto scalar=gaussian::reduce(basis,{3},d);
  check(scalar.U==d.constant(1)&&scalar.F==d.constant(3)&&scalar.multiplier==d.constant(1),"scalar Gaussian geometry");
  auto q2=gaussian::reduce(basis,{3},d,{}, {sp.dot(0,0)});
  check(q2.multiplier==Exact(field,"3-6/(2-d/2)"),"one-loop q squared equals mass squared minus denominator");
  auto q4=gaussian::reduce(basis,{3},d,{}, {sp.dot(0,0),sp.dot(0,0)});
  check(q4.multiplier==Exact(field,"9*d*(d+2)/(4*(2-d/2)*(1-d/2))"),"quartic Wick Lorentz cycles");
  auto values=gaussian::evaluate(q2,1,4,2);check(values.epsilon_low==-1&&values.coefficients.size()==4,"numerator epsilon pole lookahead");
  auto scalar_values=gaussian::evaluate(scalar,1,4,2);auto oracle=feynman::tadpole(Rational(1),Rational(3),3,1,4,2);
  check(scalar_values.epsilon_low==oracle.epsilon_low,"scalar epsilon valuation");
  for(unsigned k=0;k<scalar_values.coefficients.size();++k)check(acb_overlaps(scalar_values.coefficients[k].raw(),oracle.coefficients[k].raw()),"scalar Laurent agreement");
  auto lower=feynman::tadpole(Rational(1),Rational(3),2,1,4,2);
  for(int order=-1;order<=2;++order) {
    auto expected=-lower.coefficients.at(order-lower.epsilon_low);
    if(order>=scalar_values.epsilon_low)expected+=Jet::Ball(3)*scalar_values.at(order);
    check(acb_overlaps(values.at(order).raw(),expected.raw()),"numerator Laurent convolution agrees with difference of scalar tadpoles");
  }
  check(gaussian::reduce(basis,{-2},d).multiplier.is_zero(),"nonpositive physical power is scaleless");
  rejects([&]{gaussian::reduce(basis,{3},d,{1,100,100,100},{sp.dot(0,0),sp.dot(0,0)});},"Gaussian degree budget");
  ibp::ScalarProducts s2(2,{{d.constant(5)}},d);auto D=s2.zero();
  // Q={{2,1},{1,3}}, R={{1},{-2}}, C=7; M=7+5*(R.Q^-1.R)=22.
  D.constant=d.constant(7);D.linear={d.constant(-2),d.constant(-2),d.constant(-2),d.constant(-3),d.constant(4)};
  ibp::PropagatorBasis shifted({s2,{D},{}});auto sh=gaussian::reduce(shifted,{4,0,0,0,0},d);
  check(sh.U==d.constant(5)&&sh.F==d.constant(110),"shifted nondiagonal two-loop complete square");
  auto linear=gaussian::reduce(shifted,{4,0,0,0,0},d,{}, {s2.dot(0,2)});
  check(linear.multiplier==d.constant(-5),"shifted external linear moment");
  auto cross=gaussian::reduce(shifted,{4,0,0,0,0},d,{}, {s2.dot(0,1)});
  check(cross.multiplier==Exact(field,"-5+11*d/(5*(3-d))"),"shifted cross-loop tensor moment");
  ibp::Generator generator(shifted,d);ibp::Integral ibp_seed{4,-1,0,-1,0};
  for(const auto& relation:generator.relations(ibp_seed)) {
    auto residual=d.constant(0);
    for(const auto& [indices,c]:relation) {
      auto value=gaussian::reduce(shifted,indices,d);auto ratio=d.constant(1);const auto mass=value.F/value.U;
      for(int p=4;p<indices[0];++p)ratio=ratio*(d.constant(p)-d)/(d.constant(p)*mass);
      for(int p=indices[0];p<4;++p)ratio=ratio*d.constant(p)*mass/(d.constant(p)-d);
      residual=residual+c*value.multiplier*ratio;
    }
    check(residual.is_zero(),"Gaussian moments satisfy independently generated exact two-loop IBPs");
  }
  auto singular=s2.zero();singular.constant=d.constant(1);singular.linear[0]=d.constant(-1);
  ibp::PropagatorBasis singular_basis({s2,{singular},{}});
  rejects([&]{gaussian::reduce(singular_basis,{3,0,0,0,0},d);},"singular loop Q rejected");
  // Symbolic current parameter survives through U,F and numerator coefficients.
  auto banana=ibp::quadratic_family(feynman::banana(4,{Rational(1),Rational(1),Rational(1),Rational(1),Rational(1)}),x);
  banana=ibp::merge(banana,0,1,x.constant(Rational("1/3")));banana=ibp::merge(banana,0,1,x.constant(Rational("2/5")));
  banana=ibp::merge(banana,0,1,x.constant(Rational("3/7")));banana=ibp::merge(banana,0,1,x);
  ibp::PropagatorBasis bb(banana);ibp::Integral target(14,0);target[0]=5;target[1]=-1;target[7]=-1;target[13]=-1;
  auto br=gaussian::reduce(bb,target,d);auto anchored=gaussian::specialize(br,std::vector<Exact>{x.constant(Rational("1/2")),d});
  auto bv=gaussian::evaluate(anchored,1,2,1);check(!bv.coefficients.empty(),"Banana4 14-slot cubic numerator boundary");
  auto henn=feynman::example_family("henn_double_pentagon_x0");auto hf=ibp::quadratic_family(henn.momenta,x,henn.physical_count);
  while(hf.physical.size()>1)hf=ibp::merge(hf,0,1,x.constant(Rational("1/3")));
  ibp::PropagatorBasis hb(hf);std::size_t tested=0;
  for(unsigned degree=1;degree<=3;++degree) {
    ibp::Integral indices(11,0);indices[0]=8;
    std::function<void(unsigned,unsigned)> visit=[&](unsigned slot,unsigned remaining) {
      if(!remaining){auto b=gaussian::reduce(hb,indices,d);check(b.U.is_rational()&&b.F.is_rational(),"Henn anchored numerator geometry");++tested;return;}
      for(unsigned k=slot;k<11;++k){--indices[k];visit(k,remaining-1);++indices[k];}
    };visit(1,degree);
  }
  check(tested==285,"all Henn mixed numerator monomials through degree three");
  std::cout<<"Gaussian boundary passed; Henn slot/degree cases="<<tested<<" Banana4 slots="<<bb.denominators.size()<<"\n";
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
