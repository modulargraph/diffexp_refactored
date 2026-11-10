#include "diffexp/adjoint_transport.hpp"
#include <iostream>
using namespace diffexp;
using B=Jet::Ball;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
// Row batches share every chart endpoint, but the endpoint a tail-limited
// chart settles on depends on which batch shortened it first: results agree
// to the truncation target (2^-48 at order 24), not bit for bit.
void same(const LaurentRows& a,const LaurentRows& b){
  require(a.low==b.low && a.high==b.high && a.columns()==b.columns() && a.coefficients.size()==b.coefficients.size(),"batch result shape");
  const auto columns=a.columns();
  for(unsigned i=0;i<a.coefficients.size();++i)for(unsigned j=0;j<columns;++j)
    for(unsigned k=0;k<a.coefficients[i][j].size();++k) {
      const auto delta=a.coefficients[i][j][k]-b.coefficients[i][j][k];
      mag_t error,scale;mag_init(error);mag_init(scale);acb_get_mag(error,delta.raw());acb_get_mag(scale,a.coefficients[i][j][k].raw());
      mag_add_ui(scale,scale,1);mag_mul_2exp_si(scale,scale,-40);
      const bool close=mag_cmp(error,scale)<0;mag_clear(error);mag_clear(scale);
      require(close,"batching changed a retained coefficient beyond the truncation target");
    }
}
int main(){try {
  ExactField field({"x","eps","I"});Exact zero(field,"0"),one(field,"1"),x(field,"x"),eps(field,"eps");
  // Complex legs, Laurent forcing, epsilon-dependent pole, and unequal rows.
  for(auto bits:{128,256}) {
    B::set_precision(bits);
    ExactEpsilonMatrix matrix{{one/(x+one+eps),eps},{zero,one/(x+one)}};
    ExactEpsilonMatrix forcing(7,std::vector<Exact>(2,zero));
    ExactEpsilonMatrix rows=forcing;
    for(unsigned i=0;i<7;++i){rows[i][0]=one.constant(i+1);rows[i][1]=eps*one.constant(i+2);forcing[i][i%2]=x/eps*one.constant(i+1);}
    auto initial=exact_laurent_rows(rows,zero,2);
    std::vector<Exact> path{zero,Exact(field,"I/5"),Exact(field,"1/3+I/5"),Exact(field,"1/3")};
    AdjointOptions full;full.taylor_order=24;
    auto expected=transport_adjoint_rows(matrix,initial,forcing,path,full);
    auto bounded=full;bounded.max_taylor_cells=(2*2+1)*4*25; // two observable rows, including common source
    same(expected,transport_adjoint_rows(matrix,initial,forcing,path,bounded));
    bounded=full;bounded.max_rows_per_batch=1;
    same(expected,transport_adjoint_rows(matrix,initial,forcing,path,bounded));
  }
  {
    // Exact-term and operation budgets are charged per replicated row: a
    // batch exceeding them must split like one exceeding the cell budget and
    // reproduce the explicitly batched retained balls exactly.
    B::set_precision(256);
    const ExactEpsilonMatrix scalar{{one/(one+x)}},none(4,std::vector<Exact>{zero});
    LaurentRows rows{0,0,{{{B(1)}},{{B(2)}},{{B(3)}},{{B(4)}}}};
    AdjointOptions single;single.max_rows_per_batch=1;
    const auto expected=transport_adjoint_rows(scalar,rows,none,{zero,one},single);
    const auto one_row=rational_circuit::compile({{0,0,0,Exact(field,"-1/(1+x)")}},2,0,1,0,80);
    for(unsigned kind=0;kind<2;++kind) {
      AdjointOptions tight;AdjointConditioningStats work;tight.conditioning_stats=&work;
      if(kind==0)tight.max_circuit_operations=rational_circuit::operation_upper_bound(one_row.data(),80,1);
      else tight.max_circuit_terms=one_row.data().visited_terms;
      const auto split=transport_adjoint_rows(scalar,rows,none,{zero,one},tight);
      require(work.circuit_batch_reductions==2,"circuit budget overrun did not split the row batch");
      for(unsigned i=0;i<4;++i)require(acb_equal(split.coefficients[i][0][0].raw(),expected.coefficients[i][0][0].raw()),
        "budget-split batches changed retained balls");
      if(kind==0)--tight.max_circuit_operations;else --tight.max_circuit_terms;
      bool rejected=false;
      try{transport_adjoint_rows(scalar,rows,none,{zero,one},tight);}catch(const std::length_error&){rejected=true;}
      require(rejected,"a single row beyond the circuit budget was accepted");
    }
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
  std::cout<<"Adjoint scaling: 109 masters, 257 observables; automatic, explicit and budget-split batches agree; single-row memory floor enforced\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
