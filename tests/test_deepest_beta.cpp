#include "diffexp/deepest_beta.hpp"
#include "diffexp/bubble.hpp"
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
  mag_t tolerance;mag_init(tolerance);mag_set_ui_2exp_si(tolerance,1,-160);
  acb_calc_integrate_opt_t options;acb_calc_integrate_opt_init(options);
  options->eval_limit=100000;options->depth_limit=64;
  auto status=acb_calc_integrate(result.raw(),integrand,&oracle,B(0).raw(),B(1).raw(),160,tolerance,options,B::precision());
  mag_clear(tolerance);require(status==ARB_CALC_SUCCESS && result.is_finite(),"deepest finite-part oracle did not converge");return result;
}
int main(){try {
  B::set_precision(384);ExactField field({"x","I"});Exact x(field,"x");
  auto bubble=feynman::symanzik(feynman::banana(1,{Rational(1),Rational(1)}),{x,x.constant(1)-x});
  auto generic=feynman::deepest_beta(bubble,1,2,1,1,4);auto direct=feynman_bubble();
  double maximum=0;
  for(int k=0;k<=4;++k)maximum=std::max(maximum,bound(generic.at(k)-direct[k]));
  require(maximum<1e-38 && !generic.taylor_tail_certified,"generic beta series matches independently tested bubble window");
  // Exact meromorphic integrals test epsilon lookahead through a primitive
  // pole and dimensionally regularized power divergence, including eps^4.
  for(unsigned power:{1,2}) {
    auto result=feynman::deepest_beta({x.pow(power),x.pow(power)},2,2,1,1,4);
    auto gamma=feynman::tadpole(Rational(1),Rational(1),2,2,2,5);
    if(power==1) {
      for(int k=-2;k<=4;++k)
        require(bound(result.at(k)-gamma.coefficients[k+2])<1e-35,"resonant primitive retains gamma lookahead through epsilon 4");
    } else {
      Jet denominator(-1,6,384);denominator.set(1,B(2));auto reciprocal=denominator.constant(1)/denominator;
      require(bound(result.at(-2))<1e-35,"power divergence does not create a spurious double pole");
      for(int k=-1;k<=4;++k) {
        B expected(0);for(int j=0;j<=k+1;++j)expected+=gamma.coefficients[j]*reciprocal.at(k+1-j);
        require(bound(result.at(k)-expected)<1e-35,"power-divergent primitive agrees with analytic dimensional continuation");
      }
    }
  }
  auto sunrise=feynman::symanzik(feynman::banana(2,{Rational(1),Rational(1),Rational(1)}),
    {x/x.constant(3),x*x.constant(Rational("2/3")),x.constant(1)-x});
  auto polar=feynman::deepest_beta(sunrise,2,2,1,1,4);B euler,log;
  arb_const_euler(acb_realref(euler.raw()),384);acb_log(log.raw(),B::from_strings("2/9").raw(),384);
  require(polar.epsilon_low==-2 && bound(polar.at(-2)-B::from_strings("1/2"))<1e-35,
    "sunrise lower endpoint exposes the physical double pole");
  require(bound(polar.at(-1)+euler+log/B(2))<1e-35,"sunrise primitive residue agrees with independent endpoint subtraction");
  for(unsigned loops:{2,4}) {
    std::vector<Exact> parameters;std::vector<std::pair<unsigned,unsigned>> sequence;
    for(unsigned i=0;i<loops;++i) {
      sequence.emplace_back(0,i+1);
      parameters.push_back(i+1==loops?x:x.constant(Rational(i+1)/Rational(5)));
    }
    auto weights=feynman::merge_weights(x,loops+1,sequence,parameters);
    auto geometry=feynman::symanzik(feynman::banana(loops,std::vector<Rational>(loops+1,Rational(1))),weights);
    auto result=feynman::deepest_beta(geometry,loops,2,loops,1,4);
    // Gamma(v-L)/[Gamma(L) Gamma(1)] = 1/(L-1)! at eps=0.
    // The beta factor L still multiplies the tadpole's 1/Gamma(L+1).
    Rational prefactor(1);for(unsigned k=1;k<loops;++k)prefactor=prefactor/Rational(k);
    auto oracle=quadrature(x.constant(prefactor)*x.pow(loops-1)/geometry.F);
    const auto error=bound(result.at(0)-oracle);
    require(error<1e-35,"deepest sunrise/Banana4 finite boundary disagrees with independent quadrature");
    std::cout<<"L="<<loops<<" deepest scalar beta, finite-part oracle discrepancy "<<error<<'\n';
  }
  std::cout<<"Affine-epsilon endpoint primitives, Laurent lookahead and deepest beta boundaries passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
