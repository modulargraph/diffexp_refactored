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

// Isolated complex roots of an exact real polynomial with their exact
// multiplicities, from its squarefree factorization. Only a full-degree Acb
// isolation result is accepted.
struct IsolatedRoot {kernel::ComplexBall value;unsigned multiplicity;};
inline std::vector<IsolatedRoot> isolated_roots(const fmpq_poly_t polynomial,slong bits) {
  struct Polys {
    fmpz_poly_t z;fmpz_poly_factor_t factors;acb_poly_t a;
    Polys(){fmpz_poly_init(z);fmpz_poly_factor_init(factors);acb_poly_init(a);}
    ~Polys(){acb_poly_clear(a);fmpz_poly_factor_clear(factors);fmpz_poly_clear(z);}
  } p;
  struct Roots {
    acb_ptr values;slong size;
    explicit Roots(slong n):values(_acb_vec_init(n)),size(n){}
    ~Roots(){_acb_vec_clear(values,size);}
    Roots(const Roots&)=delete;
  };
  if(fmpq_poly_is_zero(polynomial))throw std::domain_error("identically zero singularity polynomial");
  // Irreducible factors, not only squarefree ones: a singularity polynomial
  // of a large level has degree in the hundreds and is a product of small
  // factors, whose roots are found at once, while the simultaneous iteration
  // on the whole product does not converge in any reasonable time.
  fmpq_poly_get_numerator(p.z,polynomial);fmpz_poly_factor(p.factors,p.z);
  // Roots of different factors are distinct numbers, but their enclosures
  // can overlap at the requested precision: then all are found again at
  // twice the precision.
  std::vector<IsolatedRoot> out;
  for(unsigned separation=0;;++separation,bits*=2) {
  out.clear();
  for(slong i=0;i<p.factors->num;++i) {
    auto degree=fmpz_poly_degree(p.factors->p+i);
    if(degree<1)continue;
    Roots roots(degree);
    // The factor is exact and squarefree, so its roots are simple: failure to
    // isolate them is a precision shortfall of this attempt, not of the input.
    slong working=bits,found=0;
    for(unsigned attempt=0;;++attempt) {
      acb_poly_set_fmpz_poly(p.a,p.factors->p+i,working);
      found=acb_poly_find_roots(roots.values,p.a,nullptr,1000,working);
      if(found==degree || attempt==4)break;
      working*=2;
    }
    if(found!=degree)throw std::runtime_error("complex root isolation incomplete");
    for(slong j=0;j<degree;++j) {
      IsolatedRoot root{{},static_cast<unsigned>(p.factors->exp[i])};acb_set(root.value.raw(),roots.values+j);
      out.push_back(std::move(root));
    }
  }
  bool separated=true;
  for(std::size_t i=0;i<out.size() && separated;++i)for(std::size_t j=i+1;j<out.size();++j)
    if(acb_overlaps(out[i].value.raw(),out[j].value.raw())){separated=false;break;}
  if(separated || separation==4)break;
  }
  return out;
}
// Isolated roots, with multiplicities, of an exact polynomial in one variable.
inline std::vector<IsolatedRoot> isolated_roots(const Exact& polynomial,unsigned variable,slong bits=256) {
  struct Poly {fmpq_poly_t q;Poly(){fmpq_poly_init(q);}~Poly(){fmpq_poly_clear(q);}} p;
  if(polynomial.is_zero())throw std::domain_error("identically zero singularity polynomial");
  auto denominator=polynomial.denominator().rational();
  for(const auto& term:polynomial.numerator_terms()) {
    for(unsigned i=0;i<term.powers.size();++i)
      if(i!=variable && term.powers[i])throw std::invalid_argument("singularity polynomial not univariate");
    if(term.powers[variable]>10000)throw std::invalid_argument("root isolation degree limit");
    const auto coefficient=term.coefficient/denominator;
    fmpq_poly_set_coeff_fmpq(p.q,term.powers[variable],coefficient.raw());
  }
  return isolated_roots(p.q,bits);
}
// Isolate all complex roots of a real rational polynomial, stripping repeated
// factors first.
inline std::vector<kernel::ComplexBall> polynomial_roots(const Exact& polynomial,unsigned variable,slong bits=256) {
  std::vector<kernel::ComplexBall> out;
  for(auto& root:isolated_roots(polynomial,variable,bits))out.push_back(std::move(root.value));
  return out;
}
// Exact rational bound of a real ball.
inline Rational rational_bound(const arb_t value,bool upper) {
  arf_t side;arf_init(side);
  if(upper)arb_get_ubound_arf(side,value,ARF_PREC_EXACT);else arb_get_lbound_arf(side,value,ARF_PREC_EXACT);
  fmpq_t exact;fmpq_init(exact);arf_get_fmpq(exact,side);
  auto result=Rational::from_fmpq(exact);
  fmpq_clear(exact);arf_clear(side);return result;
}
// The real zeros of an exact real polynomial strictly between two rational
// points where it does not vanish, ascending, each with a rational enclosure
// and the sign of the derivative there. A zero is accepted as real only on an
// exact change of sign over its enclosure, which also proves it simple; a
// repeated real zero in the range is an error, because no contour passes it.
struct SimpleRealRoot {Rational low,high;int slope;};
inline std::vector<SimpleRealRoot> simple_real_roots(const Exact& polynomial,unsigned variable,
    const Rational& from,const Rational& to,slong bits=256) {
  std::vector<SimpleRealRoot> out;
  if(polynomial.is_rational())return out;
  std::vector<Exact> point;for(std::size_t i=0;i<polynomial.variable_count();++i)point.push_back(polynomial.variable(i));
  const auto sign_at=[&](const Rational& at){point[variable]=polynomial.constant(at);return polynomial.substitute(point).rational().sign();};
  if(!(from<to) || !sign_at(from) || !sign_at(to))throw std::invalid_argument("real root range must have regular ordered ends");
  for(const auto& root:isolated_roots(polynomial,variable,bits)) {
    if(!arb_contains_zero(acb_imagref(root.value.raw())))continue;
    const auto low=rational_bound(acb_realref(root.value.raw()),false),high=rational_bound(acb_realref(root.value.raw()),true);
    if(high<from || to<low)continue;
    if(!(from<low) || !(high<to))throw std::runtime_error("a zero near "+low.str()+" is not separated from the end of its range");
    if(root.multiplicity>1)throw std::domain_error("repeated zero at x = "+root.value.real_midpoint(20)+" (multiplicity "+std::to_string(root.multiplicity)+")");
    const int left=sign_at(low),right=sign_at(high);
    if(!left || !right) {
      // The enclosure ends on the zero itself: step off it.
      const auto exact=!left?low:high;const Rational width=(to-from)/Rational(1L<<20);
      Rational step=width;int before=0,after=0;
      for(unsigned attempt=0;attempt<64 && (!before || !after || before==after);++attempt,step=step/Rational(2)) {
        before=from<exact-step?sign_at(exact-step):0;after=exact+step<to?sign_at(exact+step):0;
        // Only this zero may lie between the two probes.
        bool alone=true;
        for(const auto& other:isolated_roots(polynomial,variable,bits)) {
          if(acb_overlaps(other.value.raw(),root.value.raw()) || !arb_contains_zero(acb_imagref(other.value.raw())))continue;
          const auto l=rational_bound(acb_realref(other.value.raw()),false),h=rational_bound(acb_realref(other.value.raw()),true);
          if(!(h<exact-step) && !(exact+step<l))alone=false;
        }
        if(!alone)before=after=0;
      }
      if(!before || !after || before==after)throw std::runtime_error("cannot orient the zero at x = "+exact.str());
      out.push_back({exact-step*Rational(2),exact+step*Rational(2),after});
      // step was halved once more by the loop increment.
      continue;
    }
    if(left==right)continue; // a conjugate pair next to the axis, not a real zero
    out.push_back({low,high,right});
  }
  std::sort(out.begin(),out.end(),[](const auto& a,const auto& b){return a.low<b.low;});
  return out;
}

// Every step lies inside `fraction` (one quarter unless the caller controls
// convergence by other means) of the distance to every isolated singularity.
// Rounding is toward zero, so doubles only shorten certified steps.
inline double clearance_endpoint(double center,const std::vector<kernel::ComplexBall>& roots,double fraction=0.25) {
  if(!(fraction>0) || !(fraction<1))throw std::invalid_argument("chart clearance fraction must lie in (0,1)");
  double clearance=1/fraction;
  kernel::ComplexBall point;acb_set_d(point.raw(),center);
  for(const auto& root:roots) {
    auto delta=point-root;arf_t lower;arf_init(lower);
    acb_get_abs_lbound_arf(lower,delta.raw(),kernel::ComplexBall::precision());
    clearance=std::min(clearance,arf_get_d(lower,ARF_RND_FLOOR));arf_clear(lower);
  }
  double h=std::min(1-center,clearance*fraction);
  if(!(h>0) || center+h==center)throw std::runtime_error("no positive certified chart clearance");
  double endpoint=std::min(1.0,center+h);
  // Certify the endpoint after floating-point addition, not just the proposed
  // step before rounding. Move inward by ulps if the rounded sum overshoots.
  arf_t scale;arf_init(scale);arf_set_d(scale,fraction);
  struct Clear {arf_struct* value;~Clear(){arf_clear(value);}} clear_scale{scale};
  for(unsigned attempt=0;attempt<32;++attempt) {
    kernel::ComplexBall end;acb_set_d(end.raw(),endpoint);auto distance=end-point;
    arf_t used;arf_init(used);acb_get_abs_ubound_arf(used,distance.raw(),kernel::ComplexBall::precision());
    bool safe=true;
    for(const auto& root:roots) {
      auto delta=point-root;arf_t lower;arf_init(lower);acb_get_abs_lbound_arf(lower,delta.raw(),kernel::ComplexBall::precision());
      arf_mul(lower,lower,scale,kernel::ComplexBall::precision(),ARF_RND_FLOOR);
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
