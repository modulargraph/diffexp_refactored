#include "diffexp/epsilon_gauge.hpp"
#include "diffexp/certified_rational_transport.hpp"
#include <iostream>
using namespace diffexp;
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
template<class F>void fails(F f,const std::string& reason){try{f();}catch(const std::exception& e){require(std::string(e.what()).find(reason)!=std::string::npos,"unexpected epsilon gauge rejection");return;}throw std::runtime_error("invalid epsilon gauge accepted");}
int main(){try {
  ExactField field({"x","eps"});Exact x(field,"x"),eps(field,"eps"),zero(field,0),one(field,1);
  auto nilpotent=epsilon_diagonal_gauge({{zero,one/eps},{zero,zero}},1);
  require(nilpotent.shifts==std::vector<std::int64_t>{0,1},"nilpotent shear orientation/normalization");
  require(nilpotent.matrix[0][1]==one,"nilpotent epsilon pole not removed");
  // Y(0)=(0,1) gives Z(0)=(0,epsilon^-1). Solve the transformed
  // ordinary system and explicitly demand/recover physical coefficients.
  auto demand=gauged_epsilon_demand(nilpotent,1,{-1,2});
  require(demand.low==-2 && demand.high==1,"gauged demand must shift rather than pad epsilon rows");
  auto recovered=physical_epsilon_demand(nilpotent,1,demand);
  require(recovered.low==-1 && recovered.high==2,"epsilon demand inverse mapping");
  using B=kernel::ComplexBall;B::set_precision(256);
  // This supplied window is epsilon -1..1 in the gauged basis.
  auto transported=certified_rational_line({{0,1,0,nilpotent.matrix[0][1]}},
    {{B(0),B(0),B(0)},{B(1),B(0),B(0)}});
  require(acb_contains(transported.boundary[0][0].raw(),B(1).raw()),"gauged nilpotent transport must recover physical Y0=1/epsilon");
  require(transported.boundary[0][1].contains_zero(),"gauged nilpotent structural epsilon zero");

  ExactEpsilonMatrix cycle{{zero,one/eps.pow(2)},{eps.pow(2),zero}};
  auto mixed=epsilon_diagonal_gauge(cycle,1);
  require(mixed.shifts==std::vector<std::int64_t>{0,2},"mixed epsilon^-2/epsilon^2 cycle shear");
  require(mixed.matrix[0][1]==one && mixed.matrix[1][0]==one,"balanced valuation cycle must be holomorphic");
  // Verify D*Anew=Aold*D entrywise, with no numerical specialization.
  for(std::size_t i=0;i<2;++i)for(std::size_t j=0;j<2;++j)
    require(eps.pow(mixed.shifts[i])*mixed.matrix[i][j]==cycle[i][j]*eps.pow(mixed.shifts[j]),"exact diagonal gauge identity");
  ExactEpsilonMatrix observable{{one+x,(one-x)/eps}};
  auto transformed=gauge_observable_columns(mixed,observable);
  require(transformed[0][0]==one+x && transformed[0][1]==eps*(one-x),"observable column shear identity");

  auto unchanged=epsilon_diagonal_gauge({{one,eps/(one+eps*x)},{zero,one+x}},1);
  require(unchanged.shifts==std::vector<std::int64_t>{0,0},"holomorphic matrix should require no shear");
  require(!exact_epsilon_valuation(zero,1) && *exact_epsilon_valuation(eps.pow(3)/(one+eps),1)==3,"exact zero/rational epsilon valuation");
  fails([&]{epsilon_diagonal_gauge({{one/eps}},1);},"negative epsilon-valuation cycle");
  fails([&]{epsilon_diagonal_gauge({{zero,one/eps},{one,zero}},1);},"negative epsilon-valuation cycle");
  fails([&]{epsilon_diagonal_gauge({{one,zero}},1);},"square matrix");
  fails([&]{gauged_epsilon_demand(mixed,2,{0,4});},"row index");
  fails([&]{gauged_epsilon_demand(mixed,0,{4,0});},"inverted");
  fails([&]{gauged_epsilon_demand(mixed,1,{std::numeric_limits<std::int64_t>::min(),0});},"overflow");
  fails([&]{physical_epsilon_demand(mixed,1,{0,std::numeric_limits<std::int64_t>::max()});},"overflow");
  EpsilonGaugeOptions small;small.max_abs_shift=1;
  fails([&]{epsilon_diagonal_gauge(cycle,1,small);},"shift exceeds");
  small={};small.max_relaxations=1;
  fails([&]{epsilon_diagonal_gauge(cycle,1,small);},"relaxation budget");
  small={};small.max_valuation_degree=1;
  fails([&]{epsilon_diagonal_gauge(cycle,1,small);},"degree bound");
  std::cout<<"Exact epsilon diagonal gauges: nilpotent analytic transport, balanced cycles, demands, observables, essential-pole rejection and finite bounds passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
