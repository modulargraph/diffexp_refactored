#include "diffexp/ibp.hpp"
#include "diffexp/jet.hpp"
#include <flint/acb_calc.h>
#include <flint/arb_calc.h>
#include <iostream>

using namespace diffexp;
using B=Jet::Ball;
using namespace diffexp::ibp;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}

struct Parameters {unsigned a,b;bool derivative;B x,gamma;};
// An independent Feynman-parameter representation at d=1/2, where every
// selected integral is absolutely convergent. Set z=t^4 to remove its fractional
// endpoint power. This oracle neither reads nor uses the generated ODE.
int integrand(acb_ptr out,const acb_t in,void* opaque,slong order,slong bits) {
  if(order>1){for(slong i=0;i<order;++i)acb_indeterminate(out+i);return 0;}
  const auto& p=*static_cast<Parameters*>(opaque);const auto v=p.a+p.b;
  B z,u,f,lu,lf,power;
  acb_pow_ui(z.raw(),in,4,bits);
  auto q=p.x*(B(1)-p.x);u=B(1)-z+z*q;f=u+q*z*(B(1)-z);
  acb_log_analytic(lu.raw(),u.raw(),order!=0,bits);acb_log_analytic(lf.raw(),f.raw(),order!=0,bits);
  if(!lu.is_finite() || !lf.is_finite()){acb_indeterminate(out);return 0;}
  const auto alpha=B(v)-B::from_strings("3/4"),beta=B::from_strings("1/2")-B(v);
  power=alpha*lu+beta*lf;acb_exp(power.raw(),power.raw(),bits);
  B t_power,end_power;acb_pow_ui(t_power.raw(),in,4*p.a-2,bits);
  acb_pow_ui(end_power.raw(),(B(1)-z).raw(),p.b-1,bits);
  auto result=B(4)*p.gamma*t_power*end_power*power;
  if(p.derivative) {
    auto du=z*(B(1)-B(2)*p.x),df=du*(B(2)-z);
    result=result*(alpha*du/u+beta*df/f);
  }
  acb_set(out,result.raw());return 0;
}
B quadrature(unsigned a,unsigned b,bool derivative,const B& x) {
  Parameters p{a,b,derivative,x,B(0)};B ga,gb;
  acb_gamma(p.gamma.raw(),(B(a+b)-B::from_strings("1/2")).raw(),B::precision());
  acb_gamma(ga.raw(),B(a).raw(),B::precision());acb_gamma(gb.raw(),B(b).raw(),B::precision());
  p.gamma=p.gamma/(ga*gb);
  mag_t tolerance;mag_init(tolerance);mag_set_ui_2exp_si(tolerance,1,-150);
  acb_calc_integrate_opt_t options;acb_calc_integrate_opt_init(options);
  options->eval_limit=100000;options->depth_limit=64;B result;
  const auto status=acb_calc_integrate(result.raw(),integrand,&p,B(0).raw(),B(1).raw(),150,tolerance,options,B::precision());
  mag_clear(tolerance);
  require(status==ARB_CALC_SUCCESS && result.is_finite(),"sunrise parameter oracle failed to converge");
  return result;
}
int main(){try {
  B::set_precision(256);ExactField field({"x","d"});Exact x(field,"x"),d(field,"d");
  auto family=quadratic_family(feynman::banana(2,{Rational(1),Rational(1),Rational(1)}),x);
  Generator generator(PropagatorBasis(merge(family,0,1,x)),d);ExactReducer reducer(x,300);
  const auto seeds=for_each_seed(2,5,{1,1,50},[&](const Integral& seed){
    for(auto& equation:generator.relations(seed))reducer.insert(std::move(equation));
  });
  auto system=differential_system(generator,reducer,{{1,0,0,0,0},{1,1,0,0,0},{2,1,0,0,0}},0,x);
  require(system.matrix[0][0]==Exact(field,"-d*(1-2*x)/(2*x*(1-x))"),"sunrise single-propagator scaling equation");
  std::vector<Exact> point{x.constant(Rational("1/3")),d.constant(Rational("1/2"))};
  const auto anchor=B::from_strings("1/3"),u=anchor*(B(1)-anchor);B gamma,upower;
  acb_gamma(gamma.raw(),B::from_strings("1/2").raw(),256);
  acb_pow(upower.raw(),u.raw(),B::from_strings("-1/4").raw(),256);
  const auto tadpole=gamma*upower;
  std::vector<B> values{tadpole,quadrature(1,1,false,anchor),quadrature(2,1,false,anchor)};
  std::vector<B> derivatives{tadpole*B::from_strings("-1/4")*(B(1)-B(2)*anchor)/u,
    quadrature(1,1,true,anchor),quadrature(2,1,true,anchor)};
  double maximum=0;
  for(unsigned i=0;i<3;++i) {
    auto residual=-derivatives[i];
    for(unsigned j=0;j<3;++j)residual+=B::from_strings(system.matrix[i][j].substitute(point).rational().str())*values[j];
    arf_t error;arf_init(error);acb_get_abs_ubound_arf(error,residual.raw(),256);
    maximum=std::max(maximum,arf_get_d(error,ARF_RND_CEIL));arf_clear(error);
  }
  require(maximum<1e-35,"native sunrise IBP differential matrix disagrees with independent parameter integrals");
  std::cout<<"Sunrise merged level: "<<seeds<<" seeds, "<<reducer.equation_count()
    <<" exact identities, 3-integral derivative closure, oracle residual "<<maximum<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
