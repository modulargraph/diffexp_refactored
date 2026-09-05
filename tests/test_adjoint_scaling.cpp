#include "diffexp/adjoint_transport.hpp"
#include <iostream>
using namespace diffexp;
using B=Jet::Ball;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void same(const LaurentRows& a,const LaurentRows& b){
  require(a.low==b.low && a.high==b.high && a.columns()==b.columns() && a.coefficients.size()==b.coefficients.size(),"batch result shape");
  const auto columns=a.columns();
  for(unsigned i=0;i<a.coefficients.size();++i)for(unsigned j=0;j<columns;++j)
    for(unsigned k=0;k<a.coefficients[i][j].size();++k)
      require(acb_equal(a.coefficients[i][j][k].raw(),b.coefficients[i][j][k].raw()),"batch changed retained ball at identical precision/order");
}
int main(){try {
  ExactField field({"x","eps","I"});Exact zero(field,"0"),one(field,"1"),x(field,"x"),eps(field,"eps");
  // Complex legs, Laurent forcing, epsilon-dependent pole, and unequal rows.
  // Check both finite-lag and rational recurrence without loosening precision.
  for(auto bits:{128,256})for(bool polynomial:{false,true}) {
    B::set_precision(bits);
    ExactEpsilonMatrix matrix{{one/(x+one+eps),eps},{zero,one/(x+one)}};
    ExactEpsilonMatrix forcing(7,std::vector<Exact>(2,zero));
    ExactEpsilonMatrix rows=forcing;
    for(unsigned i=0;i<7;++i){rows[i][0]=one.constant(i+1);rows[i][1]=eps*one.constant(i+2);forcing[i][i%2]=x/eps*one.constant(i+1);}
    auto initial=exact_laurent_rows(rows,zero,2);
    std::vector<Exact> path{zero,Exact(field,"I/5"),Exact(field,"1/3+I/5"),Exact(field,"1/3")};
    AdjointOptions full;full.taylor_order=24;full.polynomial_recurrence=polynomial;
    auto expected=transport_adjoint_rows(matrix,initial,forcing,path,full);
    auto bounded=full;bounded.max_taylor_cells=(2*2+1)*4*25; // two observable rows, including common source
    same(expected,transport_adjoint_rows(matrix,initial,forcing,path,bounded));
    bounded=full;bounded.max_rows_per_batch=1;
    same(expected,transport_adjoint_rows(matrix,initial,forcing,path,bounded));
  }
  B::set_precision(128);
  constexpr unsigned masters=109,observables=257;
  ExactEpsilonMatrix matrix(masters,std::vector<Exact>(masters,zero));
  matrix[0][1]=one;matrix[1][2]=eps;
  ExactEpsilonMatrix forcing(observables,std::vector<Exact>(masters,zero));
  LaurentRows initial{0,1,std::vector(observables,std::vector(masters,std::vector<B>(2,B(0))))};
  for(unsigned i=0;i<observables;++i){initial.coefficients[i][0][0]=B(i+1);initial.coefficients[i][108][0]=B(i+2);forcing[i][3]=one.constant(i+3);}
  AdjointOptions full;full.taylor_order=8;
  auto expected=transport_adjoint_rows(matrix,initial,forcing,{zero,one},full);
  auto bounded=full;bounded.max_taylor_cells=(7*masters+1)*2*9;
  auto result=transport_adjoint_rows(matrix,initial,forcing,{zero,one},bounded);
  same(expected,result);
  require(acb_equal_si(result.coefficients[256][1][0].raw(),-257),"109-master transpose/sign");
  require(acb_equal_si(result.coefficients[256][3][0].raw(),259),"257th observable forcing");
  require(acb_equal_si(result.coefficients[256][108][0].raw(),258),"109th master retained");
  bounded.max_taylor_cells=(masters+1)*2*9-1;
  bool rejected=false;
  try{transport_adjoint_rows(matrix,initial,forcing,{zero,one},bounded);}catch(const std::invalid_argument& e){rejected=std::string(e.what()).find("workspace budget")!=std::string::npos;}
  require(rejected,"insufficient single-row workspace accepted");
  std::cout<<"Adjoint scaling: 109 masters, 257 observables; automatic and explicit batches exactly match retained balls at fixed order/precision; single-row memory floor enforced\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
