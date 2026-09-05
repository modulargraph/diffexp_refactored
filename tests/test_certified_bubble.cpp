#include "diffexp/certified_bubble.hpp"
#include <flint/acb_calc.h>
#include <flint/arb_calc.h>
#include <iostream>
using namespace diffexp;
using B=kernel::ComplexBall;
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
template<class F> void fails_with(F f,const std::string& expected) {
  try {f();}catch(const std::exception& e) {
    require(std::string(e.what()).find(expected)!=std::string::npos,"failure must identify the exhausted proof/precision budget");return;
  }
  throw std::runtime_error("an insufficient certification budget was accepted");
}
// Independent quadrature of Gamma(1+epsilon) (1+x-x^2)^(-1-epsilon).
// Gamma coefficients use log-Gamma constants, not the solver's gamma_series;
// no transport recurrence or connection matrix is shared with this oracle.
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
  B::set_precision(192);
  auto actual=certified_feynman_bubble();
  require(B::precision()==192,"solver must restore ambient precision");
  require(actual.coefficients.size()==5 && actual.requested_digits==20,"certified bubble output window/accuracy");
  require(actual.charts>0 && actual.charts<=16 && actual.disk_attempts<=32,"default bubble charts should remain economical");
  require(native_enclosure_meets_tolerance({actual.coefficients},actual.absolute_tolerance),"actual returned balls fail declared tolerance");
  B::set_precision(384);
  std::vector<B> log_gamma(5,B(0)),gamma(5,B(0));gamma[0]=B(1);
  arb_const_euler(acb_realref(log_gamma[1].raw()),384);log_gamma[1]=-log_gamma[1];
  for(unsigned k=2;k<=4;++k){arb_zeta_ui(acb_realref(log_gamma[k].raw()),k,384);log_gamma[k]=log_gamma[k]/B(k);if(k%2)log_gamma[k]=-log_gamma[k];}
  for(unsigned n=1;n<=4;++n){for(unsigned k=1;k<=n;++k)gamma[n]+=B(k)*log_gamma[k]*gamma[n-k];gamma[n]=gamma[n]/B(n);}
  mag_t tolerance;mag_init(tolerance);mag_set_ui_2exp_si(tolerance,1,-180);
  for(unsigned k=0;k<=4;++k) {
    OracleParameters parameters{k,gamma};B oracle;
    acb_calc_integrate_opt_t options;acb_calc_integrate_opt_init(options);options->eval_limit=100000;options->depth_limit=64;
    auto status=acb_calc_integrate(oracle.raw(),oracle_integrand,&parameters,B(0).raw(),B(1).raw(),180,tolerance,options,384);
    require(status==ARB_CALC_SUCCESS && oracle.is_finite(),"independent bubble quadrature did not converge");
    require(acb_contains(actual.coefficients[k].raw(),oracle.raw()),"certified bubble must CONTAIN rigorous quadrature, not merely agree in midpoints");
  }
  mag_clear(tolerance);
  CertifiedBubbleOptions budget;budget.max_charts=1;
  fails_with([&]{certified_feynman_bubble(budget);},"chart budget exhausted");
  require(B::precision()==384,"solver must restore precision after budget failure");
  budget=CertifiedBubbleOptions{};budget.max_disk_attempts=1;budget.initial_radius="4";
  fails_with([&]{certified_feynman_bubble(budget);},"disk-proof budget exhausted");
  budget=CertifiedBubbleOptions{};budget.taylor_order=8;
  fails_with([&]{certified_feynman_bubble(budget);},"final enclosure does not meet");
  budget=CertifiedBubbleOptions{};budget.working_bits=64;
  fails_with([&]{certified_feynman_bubble(budget);},"final enclosure does not meet");
  // A deliberately overlarge initial witness radius must shrink by proof,
  // with progress and the same final enclosure guarantee.
  budget=CertifiedBubbleOptions{};budget.initial_radius="1/2";
  auto adapted=certified_feynman_bubble(budget);
  require(adapted.disk_attempts>adapted.charts,"radius adaptation was not exercised");
  require(native_enclosure_meets_tolerance({adapted.coefficients},adapted.absolute_tolerance),"adapted certified result lost requested digits");
  std::cout<<"Certified native bubble epsilon 0..4 contains independent rigorous quadrature; 20 absolute digits, "
           <<actual.charts<<" charts; proof/precision budget failures passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
