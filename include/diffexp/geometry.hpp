#pragma once
#include "diffexp/exact.hpp"
#include "diffexp/jet.hpp"
#include <flint/acb_poly.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <set>
#include <map>

namespace diffexp {
inline Exact evaluate_exact(const data::Expr& e,const Exact& sample,const std::map<std::string,Exact>& vars) {
  if(data::number(e)) {
    if(!std::all_of(e.head.begin(),e.head.end(),[](unsigned char c){return std::isdigit(c);}))
      throw std::invalid_argument("exact polynomial data contains an approximate number");
    return sample.constant(Rational(e.head));
  }
  if(e.atom()) {
    auto p=vars.find(e.head);if(p==vars.end())throw std::invalid_argument("undeclared polynomial symbol: "+e.head);
    return p->second;
  }
  if(e.head=="neg")return -evaluate_exact(e.args.at(0),sample,vars);
  if(e.args.size()!=2)throw std::invalid_argument("unsupported polynomial expression: "+e.head);
  auto a=evaluate_exact(e.args[0],sample,vars);
  if(e.head=="^") {
    auto n=data::integer(e.args[1]);if(n>100000 || n< -100000)throw std::invalid_argument("polynomial exponent limit");
    return n<0?sample.constant(1)/a.pow(-n):a.pow(n);
  }
  auto b=evaluate_exact(e.args[1],sample,vars);
  if(e.head=="+")return a+b;if(e.head=="-")return a-b;
  if(e.head=="*")return a*b;if(e.head=="/")return a/b;
  throw std::invalid_argument("unsupported polynomial operation");
}

inline Exact reduce_square(const Exact& polynomial,unsigned variable,const Exact& square) {
  if(!polynomial.denominator().is_rational())throw std::invalid_argument("square reduction needs a polynomial");
  auto result=polynomial.constant(0);
  auto divisor=polynomial.denominator().rational();
  for(const auto& term:polynomial.numerator_terms()) {
    auto value=polynomial.constant(term.coefficient/divisor);
    for(unsigned i=0;i<term.powers.size();++i) {
      auto n=term.powers[i];
      if(i==variable) {value=value*square.pow(n/2);n%=2;}
      if(n)value=value*polynomial.variable(i).pow(n);
    }
    result=result+value;
  }
  return result;
}

inline Exact polynomial_norm(const Exact& polynomial,unsigned variable,const Exact& square) {
  bool depends=false;
  for(const auto& term:polynomial.numerator_terms())if(term.powers.at(variable))depends=true;
  if(!depends)return polynomial;
  std::vector<Exact> vars;
  for(unsigned i=0;i<polynomial.variable_count();++i)vars.push_back(polynomial.variable(i));
  vars[variable]=-vars[variable];
  return reduce_square(polynomial*polynomial.substitute(vars),variable,square);
}

// Isolate all complex roots of a real rational polynomial, stripping repeated
// factors first. Only a full-degree Acb isolation result is accepted.
inline std::vector<kernel::ComplexBall> polynomial_roots(const Exact& polynomial,unsigned variable,slong bits=256) {
  using B=kernel::ComplexBall;
  struct Polys {
    fmpq_poly_t q;fmpz_poly_t z;fmpz_poly_factor_t factors;acb_poly_t a;
    Polys(){fmpq_poly_init(q);fmpz_poly_init(z);fmpz_poly_factor_init(factors);acb_poly_init(a);}
    ~Polys(){acb_poly_clear(a);fmpz_poly_factor_clear(factors);fmpz_poly_clear(z);fmpq_poly_clear(q);}
  } p;
  if(polynomial.is_zero())throw std::domain_error("identically zero singularity polynomial");
  auto denominator=polynomial.denominator().rational();
  for(const auto& term:polynomial.numerator_terms()) {
    for(unsigned i=0;i<term.powers.size();++i)
      if(i!=variable && term.powers[i])throw std::invalid_argument("singularity polynomial not univariate");
    if(term.powers[variable]>10000)throw std::invalid_argument("root isolation degree limit");
    fmpq_t q;fmpq_init(q);fmpq_set_str(q,(term.coefficient/denominator).str().c_str(),10);fmpq_canonicalise(q);
    fmpq_poly_set_coeff_fmpq(p.q,term.powers[variable],q);fmpq_clear(q);
  }
  fmpq_poly_get_numerator(p.z,p.q);fmpz_poly_factor_squarefree(p.factors,p.z);
  std::vector<B> out;
  for(slong i=0;i<p.factors->num;++i) {
    auto degree=fmpz_poly_degree(p.factors->p+i);
    if(degree<1)continue;
    acb_poly_set_fmpz_poly(p.a,p.factors->p+i,bits);
    acb_ptr roots=_acb_vec_init(degree);
    auto found=acb_poly_find_roots(roots,p.a,nullptr,1000,bits);
    if(found!=degree) {_acb_vec_clear(roots,degree);throw std::runtime_error("complex root isolation incomplete");}
    for(slong j=0;j<degree;++j){B b;acb_set(b.raw(),roots+j);out.push_back(std::move(b));}
    _acb_vec_clear(roots,degree);
  }
  return out;
}

// Every step lies inside one quarter of the distance to every isolated
// singularity. Rounding is toward zero, so doubles only shorten certified steps.
inline double clearance_endpoint(double center,const std::vector<kernel::ComplexBall>& roots) {
  double clearance=4;
  kernel::ComplexBall point;acb_set_d(point.raw(),center);
  for(const auto& root:roots) {
    auto delta=point-root;arf_t lower;arf_init(lower);
    acb_get_abs_lbound_arf(lower,delta.raw(),kernel::ComplexBall::precision());
    clearance=std::min(clearance,arf_get_d(lower,ARF_RND_FLOOR));arf_clear(lower);
  }
  double h=std::min(1-center,clearance/4);
  if(!(h>0) || center+h==center)throw std::runtime_error("no positive certified chart clearance");
  double endpoint=std::min(1.0,center+h);
  // Certify the endpoint after floating-point addition, not just the proposed
  // step before rounding. Move inward by ulps if the rounded sum overshoots.
  for(unsigned attempt=0;attempt<32;++attempt) {
    kernel::ComplexBall end;acb_set_d(end.raw(),endpoint);auto distance=end-point;
    arf_t used;arf_init(used);acb_get_abs_ubound_arf(used,distance.raw(),kernel::ComplexBall::precision());arf_mul_2exp_si(used,used,2);
    bool safe=true;
    for(const auto& root:roots) {
      auto delta=point-root;arf_t lower;arf_init(lower);acb_get_abs_lbound_arf(lower,delta.raw(),kernel::ComplexBall::precision());
      safe=safe && arf_cmp(used,lower)<=0;arf_clear(lower);
    }
    arf_clear(used);
    if(safe && endpoint>center)return endpoint;
    endpoint=std::nextafter(endpoint,center);
  }
  throw std::runtime_error("rounded chart endpoint could not satisfy clearance");
}
inline double clearance_step(double center,const std::vector<kernel::ComplexBall>& roots) {
  return clearance_endpoint(center,roots)-center;
}

// Continue sqrt(p(t)) for a finite polynomial p along a real interval. Each
// subinterval must prove |p(t)/p(left)-1|<1, which gives one unambiguous analytic
// square root of the ratio. Splitting is local and bounded, never a global retry.
inline Jet::Ball continue_polynomial_sqrt(const Jet& p,const Jet::Ball& left,
    const Jet::Ball& right,const Jet::Ball& root_at_left,unsigned depth=0) {
  using B=Jet::Ball;
  if(depth>24)throw std::runtime_error("square-root continuation could not certify a nonzero image disk");
  if(!arb_is_zero(acb_imagref(left.raw())) || !arb_is_zero(acb_imagref(right.raw())))
    throw std::invalid_argument("root continuation parameter endpoints must be real");
  B interval;
  arb_union(acb_realref(interval.raw()),acb_realref(left.raw()),acb_realref(right.raw()),p.bits());
  auto start=p.evaluate_polynomial(left);
  if(start.contains_zero())throw std::domain_error("root continuation begins at a branch point");
  auto ratio=p.evaluate_polynomial(interval)/start-B(1);
  arf_t bound;arf_init(bound);acb_get_abs_ubound_arf(bound,ratio.raw(),p.bits());
  bool admissible=arf_cmp_ui(bound,1)<0;arf_clear(bound);
  if(admissible) {
    auto quotient=p.evaluate_polynomial(right)/start;B factor;
    acb_sqrt(factor.raw(),quotient.raw(),p.bits());return root_at_left*factor;
  }
  auto midpoint=(left+right)/B(2);
  auto middle_root=continue_polynomial_sqrt(p,left,midpoint,root_at_left,depth+1);
  return continue_polynomial_sqrt(p,midpoint,right,middle_root,depth+1);
}
} // namespace diffexp
