#include "diffexp/affine_frobenius.hpp"
#include <iostream>
using namespace diffexp;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void equal(const AffineFrobeniusSeries& a,const AffineFrobeniusSeries& b){
  require(a.dimension()==b.dimension() && a.x_order()==b.x_order(),"dimensions/order changed");
  require(a.residue_frame()==b.residue_frame(),"residue frame changed");
  require(a.exponents().size()==b.exponents().size(),"spectrum size changed");
  for(unsigned i=0;i<a.exponents().size();++i)
    require(a.exponents()[i].power==b.exponents()[i].power && a.exponents()[i].slope==b.exponents()[i].slope,"affine spectrum changed");
  require(a.jordan_successors()==b.jordan_successors(),"Jordan successors changed");
  require(a.absolute_x_frontiers()==b.absolute_x_frontiers(),"absolute frontier changed");
  const auto &x=a.terms(),&y=b.terms();
  require(x.rows==y.rows && x.columns==y.columns && x.terms.size()==y.terms.size(),"expansion shape changed");
  require(x.coherent_x_frontier==y.coherent_x_frontier && !x.omitted_tail_certified && !y.omitted_tail_certified,"certificate status changed");
  require(x.wronskian_prefactor==y.wronskian_prefactor,"Wronskian certificate changed");
  for(unsigned i=0;i<x.terms.size();++i){const auto &u=x.terms[i],&v=y.terms[i];
    require(u.row==v.row && u.column==v.column && u.log_degree==v.log_degree && u.power==v.power && u.slope==v.slope && u.coefficient==v.coefficient,"exact retained coefficient changed");}
}
int main(){try{
  ExactField field({"x","eps"});auto e=[&](const char* s){return Exact(field,s);};
  std::vector<AffineFrobeniusSeries::Matrix> examples{
    {{e("eps/x+1/(1-x)")}},
    {{e("eps/x"),e("1/x")},{e("0"),e("eps/x+1/(1-x)")}},
    {{e("0"),e("0")},{e("1/(1-x)"),e("1/x")}},
    {{e("eps/x+1/(1-x-eps*x)"),e("1/(1+x)")},{e("eps/(1-x)"),e("2*eps/x")}},
    {{e("1/x+1/(x+eps)"),e("0")},{e("1/(1-x)"),e("0")}}
  };
  for(const auto& matrix:examples)for(unsigned n:{4,8}){
    AffineFrobeniusSeries::Options old;old.finite_lag_recurrence=false;
    auto reference=AffineFrobeniusSeries::prepare(matrix,0,1,n,old);
    auto enabled=old;enabled.finite_lag_recurrence=true;
    auto finite=AffineFrobeniusSeries::prepare(matrix,0,1,n,enabled);
    require(finite.finite_lag_rows()>0,"finite lag path not exercised");
    equal(reference,finite);
    auto fallback=old;fallback.finite_lag_recurrence=true;fallback.max_clearing_degree=0;
    equal(reference,AffineFrobeniusSeries::prepare(matrix,0,1,n,fallback));
  }
  AffineFrobeniusSeries::Matrix high_degree{{e("eps/x+1/(1-x^7)")}};
  AffineFrobeniusSeries::Options forced;forced.finite_lag_recurrence=true;forced.finite_lag_cost_fallback=false;
  auto finite=AffineFrobeniusSeries::prepare(high_degree,0,1,4,forced);
  require(finite.finite_lag_rows()==1,"bounded forced clearing was not exercised");
  auto reference=forced;reference.finite_lag_recurrence=false;
  equal(finite,AffineFrobeniusSeries::prepare(high_degree,0,1,4,reference));
  std::cout<<"Finite-lag Frobenius exact retained coefficients, Jordan resonance, moving epsilon poles, frames, frontiers, Wronskians and bounded fallback match old recurrence\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
