#include "diffexp/certified_deepest_beta.hpp"
#include <flint/acb_calc.h>
#include <flint/arb_calc.h>
#include <iostream>
using namespace diffexp;
using B=Jet::Ball;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
double bound(const B& value) {
  arf_t upper;arf_init(upper);acb_get_abs_ubound_arf(upper,value.raw(),B::precision());
  auto result=arf_get_d(upper,ARF_RND_CEIL);arf_clear(upper);return result;
}
struct Oracle {data::Expr rational;};
int integrand(acb_ptr output,const acb_t input,void* opaque,slong order,slong bits) {
  if(order>1){for(slong i=0;i<order;++i)acb_indeterminate(output+i);return 0;}
  try {
    const auto& parameters=*static_cast<Oracle*>(opaque);B point;acb_set(point.raw(),input);
    Jet x(0,1,bits);x.set(0,point);auto value=evaluate(parameters.rational,x,{{"x",x}});
    acb_set(output,value.at(0).raw());
  }catch(const std::domain_error&){acb_indeterminate(output);}
  return 0;
}
B quadrature(const Exact& rational) {
  Oracle oracle{data::Reader(rational.str()).read()};B result;
  mag_t tolerance;mag_init(tolerance);mag_set_ui_2exp_si(tolerance,1,-260);
  acb_calc_integrate_opt_t options;acb_calc_integrate_opt_init(options);
  options->eval_limit=100000;options->depth_limit=64;
  auto status=acb_calc_integrate(result.raw(),integrand,&oracle,B(0).raw(),B(1).raw(),260,tolerance,options,B::precision());
  mag_clear(tolerance);require(status==ARB_CALC_SUCCESS && result.is_finite(),"deepest finite-part oracle did not converge");return result;
}
int main(){try {
  B::set_precision(384);ExactField field({"x","I"});Exact x(field,"x");
  // Exact meromorphic integrals test epsilon lookahead through a primitive
  // pole and dimensionally regularized power divergence, including eps^4.
  for(unsigned power:{1,2}) {
    auto result=feynman::certified_deepest_beta({x.pow(power),x.pow(power)},2,2,1,1,4);
    auto gamma=feynman::tadpole(Rational(1),Rational(1),2,2,2,5);
    if(power==1) {
      for(int k=-2;k<=4;++k)
        require(acb_contains(result.at(k).raw(),gamma.coefficients[k+2].raw()),"resonant primitive retains gamma lookahead through epsilon 4");
    } else {
      Jet denominator(-1,6,384);denominator.set(1,B(2));auto reciprocal=denominator.constant(1)/denominator;
      require(result.at(-2).contains_zero(),"power divergence does not create a spurious double pole");
      for(int k=-1;k<=4;++k) {
        B expected(0);for(int j=0;j<=k+1;++j)expected+=gamma.coefficients[j]*reciprocal.at(k+1-j);
        require(acb_contains(result.at(k).raw(),expected.raw()),"power-divergent primitive agrees with analytic dimensional continuation");
      }
    }
  }
  auto sunrise=feynman::symanzik(feynman::banana(2,{Rational(1),Rational(1),Rational(1)}),
    {x/x.constant(3),x*x.constant(Rational("2/3")),x.constant(1)-x});
  auto polar=feynman::certified_deepest_beta(sunrise,2,2,1,1,4);B euler,log;
  arb_const_euler(acb_realref(euler.raw()),384);acb_log(log.raw(),B::from_strings("2/9").raw(),384);
  require(polar.epsilon_low==-2 && acb_contains(polar.at(-2).raw(),B::from_strings("1/2").raw()),
    "sunrise lower endpoint exposes the physical double pole");
  auto residue=-euler-log/B(2);
  require(acb_contains(polar.at(-1).raw(),residue.raw()),"sunrise primitive residue agrees with independent endpoint subtraction");
  for(unsigned loops:{2,4}) {
    std::vector<Exact> parameters;std::vector<std::pair<unsigned,unsigned>> sequence;
    for(unsigned i=0;i<loops;++i) {
      sequence.emplace_back(0,i+1);
      parameters.push_back(i+1==loops?x:x.constant(Rational(i+1)/Rational(5)));
    }
    auto weights=feynman::merge_weights(x,loops+1,sequence,parameters);
    auto geometry=feynman::symanzik(feynman::banana(loops,std::vector<Rational>(loops+1,Rational(1))),weights);
    auto result=feynman::certified_deepest_beta(geometry,loops,2,loops,1,4);
    // Gamma(v-L)/[Gamma(L) Gamma(1)] = 1/(L-1)! at eps=0.
    // The beta factor L still multiplies the tadpole's 1/Gamma(L+1).
    Rational prefactor(1);for(unsigned k=1;k<loops;++k)prefactor=prefactor/Rational(k);
    auto oracle=quadrature(x.constant(prefactor)*x.pow(loops-1)/geometry.F);
    const auto error=bound(result.at(0)-oracle);
    require(acb_contains(result.at(0).raw(),oracle.raw()),"deepest sunrise/Banana4 finite boundary disagrees with independent quadrature");
    require(result.taylor_tail_certified && native_enclosure_meets_tolerance({result.coefficients},"1/100000000000000000000"),"certified deepest result must meet 20 actual absolute digits");
    std::cout<<"L="<<loops<<" deepest scalar beta, finite-part oracle discrepancy "<<error<<'\n';
  }
  require(polar.taylor_tail_certified && native_enclosure_meets_tolerance({polar.coefficients},"1/100000000000000000000"),"sunrise actual enclosure tolerance");
  auto anchor_fifth=feynman::symanzik(feynman::banana(2,{Rational(1),Rational(1),Rational(1)}),
    {x*x.constant(Rational("1/5")),x*x.constant(Rational("4/5")),x.constant(1)-x});
  feynman::CertifiedDeepestBetaOptions refined;refined.taylor_order=80;refined.working_bits=384;refined.requested_digits=28;
  for(unsigned left:{1,2}) {
    auto value=feynman::certified_deepest_beta(anchor_fifth,2,2,left,1,4,refined);
    require(value.taylor_tail_certified && native_enclosure_meets_tolerance({value.coefficients},"1/10000000000000000000000000000"),"N80 must reduce h/R to certify 28 digits at the fifth anchor");
  }
  const auto fails=[&](const feynman::Symanzik& g,const feynman::CertifiedDeepestBetaOptions& options,const std::string& expected) {
    bool rejected=false;
    try {feynman::certified_deepest_beta(g,2,2,1,1,4,options);}
    catch(const std::exception& error) {rejected=std::string(error.what()).find(expected)!=std::string::npos;}
    require(rejected,"insufficient certified deepest budget/domain must fail explicitly");
  };
  feynman::CertifiedDeepestBetaOptions options;options.max_charts=1;
  fails(sunrise,options,"chart budget exhausted");
  options={};options.max_disk_attempts=1;options.initial_radius="4";
  fails(sunrise,options,"endpoint disk-proof budget exhausted");
  options={};options.taylor_order=8;options.max_charts=32;
  fails(sunrise,options,"budget exhausted");
  options={};options.working_bits=64;
  fails(sunrise,options,"final enclosure does not meet");
  require(B::precision()==384,"certified deepest restores precision after failure");
  options={};fails({x,x*(x.constant(1)-x)},options,"singular upper-endpoint");
  std::cout<<"Certified affine-epsilon endpoint primitives, Laurent lookahead, sunrise poles, deepest Banana4 boundary and explicit failure contracts passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
