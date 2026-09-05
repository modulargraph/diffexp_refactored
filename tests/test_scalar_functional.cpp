#include "diffexp/scalar_functional.hpp"
#include <flint/acb_calc.h>
#include <flint/arb_calc.h>
#include <iostream>
using namespace diffexp;using B=kernel::ComplexBall;
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
struct Oracle {unsigned order,mode;std::vector<B> gamma;};
int integrand(acb_ptr output,const acb_t input,void* opaque,slong order,slong bits) {
  if(order>1){for(slong j=0;j<order;++j)acb_indeterminate(output+j);return 0;}
  const auto& p=*static_cast<Oracle*>(opaque);B x,f,l,d,sum(0);acb_set(x.raw(),input);
  f=B(1)+x;d=B(2)+x;acb_log_analytic(l.raw(),f.raw(),order!=0,bits);l=-l;
  if(f.contains_zero() || d.contains_zero() || !l.is_finite()){acb_indeterminate(output);return 0;}
  B denominator_power(1),weight;
  for(unsigned m=0;m<=p.order;++m) {
    if(p.mode==1){denominator_power=denominator_power/d;weight=m%2?-denominator_power:denominator_power;}
    else {if(m)break;weight=p.mode==0?B(1)+x+x*x:B(1);}
    B log_power(1);
    for(unsigned n=0;n+m<=p.order;++n) {
      sum+=weight*log_power*p.gamma[p.order-m-n];log_power=log_power*l/B(n+1);
    }
  }
  acb_set(output,(sum/f).raw());return 0;
}
B quadrature(unsigned power,unsigned mode,const std::vector<B>& gamma) {
  Oracle oracle{power,mode,gamma};B result;mag_t tolerance;mag_init(tolerance);mag_set_ui_2exp_si(tolerance,1,-260);
  acb_calc_integrate_opt_t options;acb_calc_integrate_opt_init(options);options->eval_limit=100000;options->depth_limit=64;
  auto status=acb_calc_integrate(result.raw(),integrand,&oracle,B(0).raw(),B(1).raw(),260,tolerance,options,384);
  mag_clear(tolerance);require(status==ARB_CALC_SUCCESS && result.is_finite(),"independent scalar functional quadrature failed");return result;
}
template<class F>void fails(F f,const std::string& reason){try{f();}catch(const std::exception& e){require(std::string(e.what()).find(reason)!=std::string::npos,"unexpected scalar functional rejection");return;}throw std::runtime_error("unsupported scalar functional accepted");}
int main(){try {
  B::set_precision(384);ExactField field({"x","eps","I"});Exact x(field,"x"),eps(field,"eps");
  feynman::Symanzik geometry{x.constant(1),x.constant(1)+x};
  std::vector<B> log_gamma(7),gamma(7);gamma[0]=B(1);
  arb_const_euler(acb_realref(log_gamma[1].raw()),384);log_gamma[1]=-log_gamma[1];
  for(unsigned k=2;k<=6;++k){arb_zeta_ui(acb_realref(log_gamma[k].raw()),k,384);log_gamma[k]=log_gamma[k]/B(k);if(k%2)log_gamma[k]=-log_gamma[k];}
  for(unsigned n=1;n<=6;++n){for(unsigned k=1;k<=n;++k)gamma[n]+=B(k)*log_gamma[k]*gamma[n-k];gamma[n]=gamma[n]/B(n);}
  auto zero=feynman::certified_scalar_functional(geometry,1,2,2,x.constant(0),4);
  require(zero.taylor_tail_certified && zero.coefficients.size()==5,"zero functional must remain exactly certified");
  for(const auto& coefficient:zero.coefficients)require(coefficient.is_zero(),"zero functional acquired numerical error");
  auto polynomial=feynman::certified_scalar_functional(geometry,1,2,2,Exact(field,"1+x+x^2"),4);
  auto rational=feynman::certified_scalar_functional(geometry,1,2,2,Exact(field,"1/(2+x+eps)"),4);
  // Exact cancellation must precede joint-analyticity/valuation tests.
  auto shifted=feynman::certified_scalar_functional(geometry,1,2,2,Exact(field,"(eps*x+eps^2)/(eps^3*(x+eps))"),4);
  require(shifted.epsilon_low==-2,"weight epsilon^-2 valuation must be preserved");
  for(int k=0;k<=4;++k) {
    auto p=quadrature(k,0,gamma),r=quadrature(k,1,gamma);
    require(acb_contains(polynomial.at(k).raw(),p.raw()),"polynomial weight must contain rigorous quadrature");
    require(acb_contains(rational.at(k).raw(),r.raw()),"rational joint weight must contain rigorous quadrature");
  }
  for(int k=-2;k<=4;++k) {
    auto oracle=quadrature(k+2,2,gamma);
    require(acb_contains(shifted.at(k).raw(),oracle.raw()),"epsilon^-2 weight requires gamma lookahead through epsilon six");
  }
  // Resonant DR primitive: Gamma(2eps) * integral x^(-1+eps) dx.
  auto resonant=feynman::certified_scalar_functional({x,x},2,2,2,x.constant(1),4);
  auto seed=feynman::tadpole(Rational(1),Rational(1),2,2,2,5);
  for(int k=-2;k<=4;++k)require(acb_contains(resonant.at(k).raw(),seed.coefficients[k+2].raw()),"scalar functional resonant pole/lookahead");
  // Existing beta functional is exactly the weight 2*x with v=3.
  auto beta=feynman::certified_deepest_beta(geometry,1,2,2,1,4);
  auto weighted_beta=feynman::certified_scalar_functional(geometry,1,2,3,Exact(field,"2*x"),4);
  for(int k=0;k<=4;++k)require(acb_overlaps(beta.at(k).raw(),weighted_beta.at(k).raw()),"general functional/beta equivalence");
  feynman::CertifiedDeepestBetaOptions refined;refined.taylor_order=80;refined.working_bits=384;refined.requested_digits=28;
  auto refined_functional=feynman::certified_scalar_functional(geometry,1,2,2,Exact(field,"1/(2+x+eps)"),4,refined);
  require(native_enclosure_meets_tolerance({refined_functional.coefficients},"1/10000000000000000000000000000"),"functional must adapt step/radius independently at N80");
  for(const auto* value:{&polynomial,&rational,&shifted,&resonant,&weighted_beta})
    require(value->taylor_tail_certified && native_enclosure_meets_tolerance({value->coefficients},"1/100000000000000000000"),"scalar functional final actual tolerance");
  // Independent Euler beta identities check both singular endpoint arms,
  // including two resonant DR primitives and their gamma lookahead.
  feynman::Symanzik two_ends{x.constant(1),x*(x.constant(1)-x)};
  Jet e(0,6,384);e.set(1,B(1));auto one=e.constant(1);
  auto regular_beta=(one+e).gamma()*(one-e).gamma().pow(2)/(e.constant(2)-e.constant(2)*e).gamma();
  auto dr_beta=-(one+e).gamma()*(one-e).gamma().pow(2)/(one-e.constant(2)*e).gamma();
  auto both=feynman::certified_scalar_functional(two_ends,1,4,2,x.constant(1),4);
  auto both_dr=feynman::certified_scalar_functional(two_ends,1,4,3,x.constant(1),4);
  feynman::CertifiedDeepestBetaOptions alternate;alternate.scalar_split_anchor="1/3";
  auto asymmetric=feynman::certified_scalar_functional(two_ends,1,4,3,x.constant(1),4,alternate);
  for(int k=-1;k<=4;++k) {
    require(acb_contains(both.at(k).raw(),regular_beta.at(k+1).raw()),"two endpoint Euler beta enclosure");
    require(acb_contains(both_dr.at(k).raw(),dr_beta.at(k+1).raw()),"two resonant endpoint DR enclosure");
    require(acb_contains(asymmetric.at(k).raw(),dr_beta.at(k+1).raw()),"asymmetric reversal Jacobian enclosure");
  }
  for(const auto* value:{&both,&both_dr,&asymmetric})
    require(value->taylor_tail_certified && native_enclosure_meets_tolerance({value->coefficients},"1/100000000000000000000"),"both endpoint tails and summed tolerance must be certified");
  feynman::Symanzik negative{x.constant(1),-x.constant(1)-x};
  feynman::CertifiedDeepestBetaOptions lower;lower.f_rim=-1;
  auto plus=feynman::certified_scalar_functional(negative,1,2,2,x.constant(1),4);
  auto minus=feynman::certified_scalar_functional(negative,1,2,2,x.constant(1),4,lower);
  auto deepest_plus=feynman::certified_deepest_beta(negative,1,2,1,1,4);
  auto deepest_minus=feynman::certified_deepest_beta(negative,1,2,1,1,4,lower);
  for(int k=0;k<=4;++k) {
    B conjugate;acb_conj(conjugate.raw(),plus.at(k).raw());
    require(acb_overlaps(conjugate.raw(),minus.at(k).raw()),"scalar negative-F rim conjugation");
    acb_conj(conjugate.raw(),deepest_plus.at(k).raw());
    require(acb_overlaps(conjugate.raw(),deepest_minus.at(k).raw()),"deepest negative-F rim conjugation");
  }
  require(arb_is_positive(acb_imagref(plus.at(1).raw())) && arb_is_negative(acb_imagref(minus.at(1).raw())),"F rim must select opposite nonzero imaginary coefficients");
  auto invalid=lower;invalid.f_rim=0;
  fails([&]{feynman::certified_scalar_functional(negative,1,2,2,x.constant(1),4,invalid);},"F rim");
  invalid=lower;invalid.scalar_split_anchor="0";
  fails([&]{feynman::certified_scalar_functional(two_ends,1,4,2,x.constant(1),4,invalid);},"split anchor");
  fails([&]{feynman::certified_scalar_functional(geometry,1,2,2,Exact(field,"1/(x+eps)"),4);},"jointly analytic");
  fails([&]{feynman::certified_scalar_functional(geometry,1,2,2,Exact(field,"1/(1-x+eps)"),4);},"jointly analytic");
  fails([&]{feynman::certified_scalar_functional(geometry,1,2,2,Exact(field,"1/x"),4);},"not regulated by epsilon");
  feynman::CertifiedDeepestBetaOptions budget;budget.max_charts=1;
  fails([&]{feynman::certified_scalar_functional(geometry,1,2,2,x.constant(1),4,budget);},"chart budget exhausted");
  fails([&]{feynman::certified_scalar_functional(two_ends,1,4,2,x.constant(1),4,budget);},"chart budget exhausted");
  require(B::precision()==384,"functional must restore precision");
  std::cout<<"Certified scalar gamma functional: polynomial/rational weights, exact cancellations, eps^-2 lookahead, DR poles, two endpoints, F rims, beta equivalence and failure contracts passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
