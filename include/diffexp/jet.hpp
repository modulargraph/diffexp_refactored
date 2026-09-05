#pragma once
#include "diffexp/kernel/scalar.hpp"
#include "diffexp/data.hpp"
#include <flint/acb_poly.h>
#include <map>
#include <limits>

namespace diffexp {
class Jet {
 public:
  using Ball=kernel::ComplexBall;
  explicit Jet(long value,unsigned length,slong bits):length_(length),bits_(bits) {
    if(!length || bits<64)throw std::invalid_argument("invalid jet context");
    acb_poly_init(p_);acb_poly_set_si(p_,value);
  }
  Jet(const Jet& b):Jet(0,b.length_,b.bits_){acb_poly_set(p_,b.p_);}
  Jet(Jet&& b) noexcept:Jet(0,b.length_,b.bits_){acb_poly_swap(p_,b.p_);}
  ~Jet(){acb_poly_clear(p_);}
  Jet& operator=(Jet b) noexcept {acb_poly_swap(p_,b.p_);std::swap(length_,b.length_);std::swap(bits_,b.bits_);return *this;}
  unsigned length() const{return length_;}
  slong bits() const{return bits_;}
  Ball at(unsigned n) const {Ball b;acb_poly_get_coeff_acb(b.raw(),p_,n);return b;}
  void set(unsigned n,const Ball& b){if(n>=length_)throw std::out_of_range("jet index");acb_poly_set_coeff_acb(p_,n,b.raw());}
  Jet constant(long n)const{return Jet(n,length_,bits_);}
  Jet decimal(std::string s)const {
    auto backtick=s.find('`');
    std::string precision;
    if(backtick!=std::string::npos) {
      auto exponent=s.find("*^",backtick);
      precision=s.substr(backtick+1,exponent==std::string::npos?std::string::npos:exponent-backtick-1);
      s=s.substr(0,backtick)+(exponent==std::string::npos?"":s.substr(exponent));
    }
    for(std::size_t i;(i=s.find("*^"))!=std::string::npos;)s.replace(i,2,"e");
    Jet j=constant(0);Ball b;
    if(arb_set_str(acb_realref(b.raw()),s.c_str(),bits_) || !b.is_finite())
      throw std::invalid_argument("invalid finite decimal literal");
    if(!precision.empty()) {
      // Published precision annotations describe significant decimal digits.
      const auto p=static_cast<long>(std::floor(std::stod(precision)));
      if(p<1 || p>1000000)throw std::invalid_argument("invalid input precision annotation");
      arb_t radius;arb_init(radius);
      arb_abs(radius,acb_realref(b.raw()));
      arb_t scale;arb_init(scale);arb_set_ui(scale,10);arb_pow_ui(scale,scale,p,bits_);
      arb_div(radius,radius,scale,bits_);arb_add_error(acb_realref(b.raw()),radius);
      arb_clear(scale);arb_clear(radius);
    }
    j.set(0,b);return j;
  }
  Jet derivative()const{Jet j=constant(0);acb_poly_derivative(j.p_,p_,bits_);return j;}
  Ball evaluate_polynomial(const Ball& point)const {
    Ball b;acb_poly_evaluate(b.raw(),p_,point.raw(),bits_);return b;
  }
  Jet log()const{
    if(at(0).contains_zero())throw std::domain_error("logarithmic branch point at chart center");
    Jet j=constant(0);acb_poly_log_series(j.p_,p_,length_,bits_);return j;
  }
  Jet exp()const {Jet j=constant(0);acb_poly_exp_series(j.p_,p_,length_,bits_);return j;}
  Jet gamma()const {
    Jet j=constant(0);acb_poly_gamma_series(j.p_,p_,length_,bits_);
    for(unsigned n=0;n<length_;++n)if(!j.at(n).is_finite())throw std::domain_error("gamma jet is not analytic at its constant term");
    return j;
  }
  Jet sqrt()const{
    if(at(0).contains_zero())throw std::domain_error("square-root branch point at chart center");
    Jet j=constant(0);acb_poly_sqrt_series(j.p_,p_,length_,bits_);return j;
  }
  Jet shifted_down(unsigned n)const {
    for(unsigned i=0;i<n;++i)if(!at(i).is_zero())throw std::domain_error("jet does not have an exact zero prefix");
    Jet out=constant(0);acb_poly_shift_right(out.p_,p_,n);return out;
  }
  Jet pow(long n)const{
    if(n==std::numeric_limits<long>::min())throw std::overflow_error("jet exponent");
    if(n<0)return constant(1)/pow(-n);
    Jet a=*this,out=constant(1);
    while(n){if(n&1)out=out*a;n>>=1;if(n)a=a*a;}return out;
  }
  friend Jet operator+(const Jet& a,const Jet& b){a.compatible(b);Jet c=a.constant(0);acb_poly_add_series(c.p_,a.p_,b.p_,a.length_,a.bits_);return c;}
  friend Jet operator-(const Jet& a,const Jet& b){a.compatible(b);Jet c=a.constant(0);acb_poly_sub_series(c.p_,a.p_,b.p_,a.length_,a.bits_);return c;}
  friend Jet operator-(const Jet& a){Jet c=a.constant(0);acb_poly_neg(c.p_,a.p_);return c;}
  friend Jet operator*(const Jet& a,const Jet& b){a.compatible(b);Jet c=a.constant(0);acb_poly_mullow(c.p_,a.p_,b.p_,a.length_,a.bits_);return c;}
  friend Jet operator/(const Jet& a,const Jet& b){a.compatible(b);if(b.at(0).contains_zero())throw std::domain_error("jet division by a zero-containing constant");Jet c=a.constant(0);acb_poly_div_series(c.p_,a.p_,b.p_,a.length_,a.bits_);return c;}
 private:
  unsigned length_;slong bits_;acb_poly_t p_;
  void compatible(const Jet& b)const{if(length_!=b.length_ || bits_!=b.bits_)throw std::invalid_argument("incompatible jet contexts");}
};

inline Jet evaluate(const data::Expr& e,const Jet& context,const std::map<std::string,Jet>& vars) {
  if(data::number(e))return context.decimal(e.head);
  if(e.atom()) {
    if(e.head=="Pi"){auto j=context.constant(0);Jet::Ball b;acb_const_pi(b.raw(),context.bits());j.set(0,b);return j;}
    if(e.head=="I"){auto j=context.constant(0);j.set(0,Jet::Ball::from_strings("0","1"));return j;}
    if(auto i=vars.find(e.head);i!=vars.end())return i->second;
    throw std::invalid_argument("unknown expression symbol: "+e.head);
  }
  const auto a=[&](unsigned i){return evaluate(e.args.at(i),context,vars);};
  if(e.head=="+" && e.args.size()==2)return a(0)+a(1);
  if(e.head=="-" && e.args.size()==2)return a(0)-a(1);
  if(e.head=="*" && e.args.size()==2)return a(0)*a(1);
  if(e.head=="/" && e.args.size()==2)return a(0)/a(1);
  if(e.head=="neg" && e.args.size()==1)return -a(0);
  if(e.head=="^" && e.args.size()==2)return a(0).pow(data::integer(e.args[1]));
  if(e.head=="Sqrt" && e.args.size()==1)return a(0).sqrt();
  if(e.head=="Log" && e.args.size()==1)return a(0).log();
  throw std::invalid_argument("unsupported ancillary expression: "+e.head);
}
}  // namespace diffexp
