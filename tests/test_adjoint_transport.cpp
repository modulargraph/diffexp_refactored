#include "diffexp/adjoint_transport.hpp"
#include <iostream>
using namespace diffexp;using B=Jet::Ball;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
double size(const B& b){return NativeTailMagnitude::upper_abs(b).approximate_upper();}
template<class F>void fails(F f,const std::string& text){try{f();}catch(const std::exception& e){require(std::string(e.what()).find(text)!=std::string::npos,"unexpected adjoint failure");return;}throw std::runtime_error("unsupported adjoint call accepted");}
int main(){try {
  B::set_precision(384);ExactField field({"x","eps","I"});Exact x(field,"x"),e(field,"eps"),zero(field,"0"),one(field,"1");
  ExactEpsilonMatrix a{{zero,one},{zero,zero}},forcing{{x,one/e},{zero,zero}};
  auto initial=exact_laurent_rows({{one,zero},{zero,one}},zero,4);
  auto output=transport_adjoint_rows(a,initial,forcing,{zero,Exact(field,"I/3"),Exact(field,"1+I/3"),one});
  require(output.low==-1,"forcing epsilon pole must extend operator lower window");
  require(size(output.coefficients[0][0][1]-B::from_strings("3/2"))<1e-80,"nondiagonal forced adjoint first component");
  require(size(output.coefficients[0][1][0]-B(1))<1e-80 && size(output.coefficients[0][1][1]+B::from_strings("7/6"))<1e-80,"nondiagonal adjoint sign, transpose and Laurent forcing");
  require(size(output.coefficients[1][1][1]-B(1))<1e-80,"batched independent row propagation");
  auto backward=transport_adjoint_rows(a,output,forcing,{one,zero});
  for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j)for(int k=-1;k<=4;++k)
    require(size(backward.coefficients[i][j][k+1]-(k>=0?initial.coefficients[i][j][k]:B(0)))<1e-75,"adjoint roundtrip exact analytic polynomial");
  auto epsilon=transport_adjoint_rows({{zero,e},{zero,zero}},initial,{{zero,zero},{zero,zero}},{zero,one});
  require(size(epsilon.coefficients[0][1][1]+B(1))<1e-80,"adjoint epsilon coupling");
  auto exponential=transport_adjoint_rows({{one+one}},exact_laurent_rows({{one}},zero,4),{{zero}},{zero,one});
  B expected;acb_exp(expected.raw(),B(-2).raw(),384);
  require(size(exponential.coefficients[0][0][0]-expected)<1e-75,"adjoint scalar exponential independent solution");
  // Combine shared-source coefficients before the one uncertain application.
  auto left=exact_laurent_rows({{Exact(field,"100000000000000000000+1")}},zero,0);
  auto right=exact_laurent_rows({{Exact(field,"-100000000000000000000")}},zero,0);
  B uncertain(1);arb_add_error_2exp_si(acb_realref(uncertain.raw()),-70);
  LaurentBoundary boundary{0,{{uncertain}},true};auto combined=apply_laurent_rows(add_laurent_rows(left,right),boundary,0);
  auto separate=apply_laurent_rows(left,boundary,0).values[0][0]+apply_laurent_rows(right,boundary,0).values[0][0];
  require(size(combined.values[0][0]-B(1))<1e-20 && size(separate-B(1))>0.1,"shared source cancellation must precede interval application");
  fails([&]{apply_laurent_rows(output,{0,{{B(1)},{B(1)}},true},0);},"additional source");
  fails([&]{apply_laurent_rows(initial,{-1,{{B(1),B(1)},{B(1),B(1)}},true},4);},"additional source");
  fails([&]{transport_adjoint_rows({{one/e}},exact_laurent_rows({{one}},zero,0),{{zero}},{zero,one});},"epsilon gauge");
  AdjointOptions budget;budget.max_charts_per_leg=1;
  fails([&]{transport_adjoint_rows({{one/(x+one)}},exact_laurent_rows({{one}},zero,0),{{zero}},{zero,one},budget);},"chart budget");
  fails([&]{apply_laurent_rows(initial,{-1,Boundary(2,std::vector<B>(7,B(1))),true},4);},"additional operator");
  auto many=exact_laurent_rows(ExactEpsilonMatrix(128,std::vector<Exact>{one}),zero,0);
  auto many_result=transport_adjoint_rows({{zero}},many,ExactEpsilonMatrix(128,std::vector<Exact>{zero}),{zero,one});
  require(many_result.coefficients.size()==128 && acb_equal_si(many_result.coefficients.back()[0][0].raw(),1),"observable row count is independent of master dimension");
  budget.max_taylor_cells=1;
  fails([&]{transport_adjoint_rows({{zero}},many,ExactEpsilonMatrix(128,std::vector<Exact>{zero}),{zero,one},budget);},"workspace budget");
  require(!combined.taylor_tail_certified,"adjoint retained Taylor rows must remain uncertified");
  std::cout<<"Adjoint rows: nondiagonal/forced/epsilon/complex solutions, shared-source cancellation, windows and finite budgets passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
