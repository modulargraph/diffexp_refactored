#include <diffexp/ibp_solver_provider.hpp>
int main(){
  diffexp::ExactField field({"d"});diffexp::Exact d(field,"d");
  diffexp::ibp::ScalarProducts products(1,{},d);auto denominator=products.zero();denominator.constant=d.constant(1);denominator.linear[0]=-d.constant(1);
  diffexp::ibp::PropagatorBasis basis({products,{denominator},{}});
  diffexp::ibp_solver::Sampler sampler(basis,d,{});auto result=sampler({{2}},{7},diffexp::ibp_solver::prime(1),{});
  if(!result.success)return 1;
  auto expected=diffexp::modular::rational_mod(diffexp::Rational("-5/2"),diffexp::ibp_solver::prime(1));
  return result.reductions.at({2}).at({1}).str()!=std::to_string(expected);
}
