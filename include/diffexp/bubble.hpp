#pragma once
#include "diffexp/rational_transport.hpp"
#include "diffexp/tadpole.hpp"
#include "diffexp/feynman.hpp"
#include "diffexp/ibp.hpp"

namespace diffexp {
// One complete native FT level. The deepest one-denominator integral supplies
// the gamma boundary; the beta functional is an accumulator ODE row.
inline std::vector<Jet::Ball> feynman_bubble(unsigned epsilon_order=4) {
  using B=Jet::Ball;
  ExactField field({"x","I","eps"});Exact x(field,"x");
  auto family=feynman::banana(1,{Rational(1),Rational(1)});
  auto geometry=feynman::symanzik(family,{x,x.constant(1)-x});
  auto request=feynman::merge_request({1,1},0,1);
  auto quadratic=ibp::quadratic_family(family,x);
  ibp::Generator equations(ibp::PropagatorBasis(ibp::merge(quadratic,0,1,x)),Exact(field,"2-2*eps"));
  ibp::ExactReducer reduction(x);
  for(int power=1;power<=4;++power) for(int isp=0;isp>=-2;--isp)
    for(auto& identity:equations.relations({power,isp})) reduction.insert(std::move(identity));
  auto system=ibp::differential_system(equations,reduction,{{2,0}},0,x);
  auto substitutions=std::vector<Exact>{x,x.variable(1),x.constant(0)};
  const auto connection0=system.matrix[0][0].substitute(substitutions);
  const auto connection1=system.matrix[0][0].derivative(2).substitute(substitutions);
  if(!system.matrix[0][0].derivative(2).derivative(2).is_zero())
    throw std::logic_error("bubble connection unexpectedly nonaffine in epsilon");
  auto origin=substitutions;origin[0]=x.constant(0);
  auto seed=feynman::tadpole(geometry.U.substitute(origin).rational(),geometry.F.substitute(origin).rational(),
    request.source_indices[0],family.loops,2,epsilon_order);
  if(seed.epsilon_low!=0)throw std::logic_error("bubble deepest boundary unexpectedly polar");
  Boundary boundary{seed.coefficients,std::vector<B>(epsilon_order+1,B(0))};
  auto integrand=x.constant(request.normalization)*x.pow(request.left_power)*(x.constant(1)-x).pow(request.right_power);
  std::vector<RationalLineEntry> entries{{0,0,0,connection0},{0,0,1,connection1},{1,0,0,integrand}};
  return rational_line(entries,std::move(boundary),72)[1];
}
inline int run_bubble() {
  using B=Jet::Ball;B::set_precision(384);auto result=feynman_bubble();
  for(unsigned k=0;k<result.size();++k)std::cout<<"Bubble epsilon^"<<k<<" = "<<result[k].real_midpoint(40)<<'\n';
  // Independent elementary epsilon^0 integral: 4 atanh(1/sqrt(5))/sqrt(5).
  B root,argument,atanh;acb_sqrt(root.raw(),B(5).raw(),B::precision());argument=B(1)/root;
  acb_atanh(atanh.raw(),argument.raw(),B::precision());auto reference=B(4)*atanh/root;
  auto error=result[0]-reference;arf_t bound;arf_init(bound);acb_get_abs_ubound_arf(bound,error.raw(),B::precision());
  const auto discrepancy=arf_get_d(bound,ARF_RND_CEIL);arf_clear(bound);
  std::cout<<"Bubble finite-part discrepancy = "<<discrepancy<<"; the test suite independently checks epsilon 0 through 4. Taylor tails are not yet certified.\n";
  return std::isfinite(discrepancy) && discrepancy<1e-30?0:1;
}
} // namespace diffexp
