#include "diffexp/bubble.hpp"
#include <flint/acb_calc.h>
#include <flint/arb_calc.h>
#include <iostream>
using namespace diffexp;
using B=Jet::Ball;
void require(bool b,const char* why){if(!b)throw std::runtime_error(why);}
struct OracleParameters {unsigned power;std::vector<B> gamma;};
int oracle_integrand(acb_ptr output,const acb_t input,void* parameters,slong order,slong bits) {
  if(order>1){for(slong i=0;i<order;++i)acb_indeterminate(output+i);return 0;}
  const auto& p=*static_cast<OracleParameters*>(parameters);
  B f,log,term(1),sum(0),contribution;
  acb_mul(f.raw(),input,input,bits);acb_neg(f.raw(),f.raw());acb_add(f.raw(),f.raw(),input,bits);acb_add_ui(f.raw(),f.raw(),1,bits);
  acb_log_analytic(log.raw(),f.raw(),order!=0,bits);acb_neg(log.raw(),log.raw());
  if(!f.is_finite() || f.contains_zero() || !log.is_finite()){acb_indeterminate(output);return 0;}
  for(unsigned j=0;j<=p.power;++j) {
    acb_mul(contribution.raw(),p.gamma[p.power-j].raw(),term.raw(),bits);acb_add(sum.raw(),sum.raw(),contribution.raw(),bits);
    acb_mul(term.raw(),term.raw(),log.raw(),bits);acb_div_ui(term.raw(),term.raw(),j+1,bits);
  }
  acb_div(output,sum.raw(),f.raw(),bits);return 0;
}
int main(){try {
  B::set_precision(384);
  const Rational u("13/10"),f("7/5");
  auto negative=feynman::tadpole(u,f,1,2,4,8);
  auto expected=-f*f*f/(Rational(12)*u*u*u*u*u);
  require(negative.epsilon_low==-1 && (negative.coefficients[0]-B::from_strings(expected.str())).contains_zero(),"negative gamma-shift Laurent residue");
  auto zero=feynman::tadpole(u,f,2,1,4,8);
  require(zero.epsilon_low==-1 && (zero.coefficients[0]-B::from_strings((Rational(1)/(u*u)).str())).contains_zero(),"zero gamma-shift Laurent residue");
  auto positive=feynman::tadpole(u,f,4,1,4,8);
  require(positive.epsilon_low==0 && (positive.coefficients[0]-B::from_strings((Rational(1)/(Rational(6)*f*f)).str())).contains_zero(),"positive gamma-shift boundary");
  auto upper=feynman::tadpole(u,-f,1,2,4,8,1),lower=feynman::tadpole(u,-f,1,2,4,8,-1);
  for(unsigned i=0;i<upper.coefficients.size();++i) {
    B conjugate;acb_conj(conjugate.raw(),upper.coefficients[i].raw());
    require((conjugate-lower.coefficients[i]).contains_zero(),"explicit negative-F rims must conjugate");
  }
  auto actual=feynman_bubble(4);
  // Independent log-Gamma coefficient recurrence, rather than gamma_series.
  std::vector<B> log_gamma(5,B(0)),gamma(5,B(0));gamma[0]=B(1);
  arb_const_euler(acb_realref(log_gamma[1].raw()),384);log_gamma[1]=-log_gamma[1];
  for(unsigned k=2;k<=4;++k){arb_zeta_ui(acb_realref(log_gamma[k].raw()),k,384);log_gamma[k]=log_gamma[k]/B(k);if(k%2)log_gamma[k]=-log_gamma[k];}
  for(unsigned n=1;n<=4;++n){for(unsigned k=1;k<=n;++k)gamma[n]+=B(k)*log_gamma[k]*gamma[n-k];gamma[n]=gamma[n]/B(n);}
  mag_t tolerance;mag_init(tolerance);mag_set_ui_2exp_si(tolerance,1,-180);
  double maximum=0;
  for(unsigned k=0;k<=4;++k) {
    OracleParameters parameters{k,gamma};B oracle;
    acb_calc_integrate_opt_t options;acb_calc_integrate_opt_init(options);options->eval_limit=100000;options->depth_limit=64;
    auto status=acb_calc_integrate(oracle.raw(),oracle_integrand,&parameters,B(0).raw(),B(1).raw(),180,tolerance,options,384);
    require(status==ARB_CALC_SUCCESS && oracle.is_finite(),"independent bubble quadrature did not converge");
    auto difference=actual[k]-oracle;arf_t error;arf_init(error);acb_get_abs_ubound_arf(error,difference.raw(),384);
    maximum=std::max(maximum,arf_get_d(error,ARF_RND_CEIL));arf_clear(error);
  }
  mag_clear(tolerance);
  require(maximum<1e-38,"Feynman bubble epsilon coefficients disagree with rigorous quadrature");
  std::cout<<"Tadpole residues/rims and bubble epsilon 0..4 passed; maximum oracle discrepancy "<<maximum<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
