#include "diffexp/singular_transport.hpp"
#include <iostream>
using namespace diffexp;
using B=Jet::Ball;
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
void near(const B& a,const B& b){check(acb_overlaps(a.raw(),b.raw()),"singular continuation analytic mismatch");}
void residual(const FrobeniusCrossing& crossing,const std::vector<std::vector<std::string>>& source) {
  ExactField field({"x","eps","ell"});Exact x(field,"x"),eps(field,"eps"),ell(field,"ell");
  auto d=crossing.series.dimension();
  std::vector<std::vector<Exact>> f(d,std::vector<Exact>(d,x.constant(0)));
  for(const auto& term:crossing.series.monomials()) {
    check(term.power.str().find('/')==std::string::npos,"analytic residual test needs integer powers");
    f[term.row][term.column]=f[term.row][term.column]+x.constant(term.coefficient)*x.pow(std::stol(term.power.str()))*eps.pow(term.epsilon)*ell.pow(term.log_degree);
  }
  for(unsigned i=0;i<d;++i)for(unsigned j=0;j<d;++j) {
    auto r=f[i][j].derivative(0)+f[i][j].derivative(2)/x;
    for(unsigned k=0;k<d;++k)r=r-Exact(field,source[i][k])*f[k][j];
    for(const auto& term:r.numerator_terms())
      check(term.powers[1]>crossing.series.epsilon_order(),"exact retained singular-series ODE residual is nonzero");
  }
}
int main(int argc,char** argv){try{
  B::set_precision(192);
  bool failed=false;
  try{(void)original_banana_equal_real("/nonexistent-diffexp3-crossing-fixture");}
  catch(const std::exception&){failed=true;}
  check(failed && B::precision()==192,"failed fixture changed ambient precision");
  B invalid;acb_indeterminate(invalid.raw());failed=false;
  try{(void)singular_transport_detail::finite_discrepancy(invalid,B(0));}
  catch(const std::domain_error&){failed=true;}
  check(failed,"nonfinite endpoint was hidden by discrepancy aggregation");
  failed=false;
  try{(void)singular_transport_detail::finite_discrepancy(B(0),invalid);}
  catch(const std::domain_error&){failed=true;}
  check(failed,"nonfinite reference was hidden by discrepancy aggregation");
  B::set_precision(384);ExactField field({"x","eps"});auto e=[&](const char* s){return Exact(field,s);};
  B pi;acb_const_pi(pi.raw(),384);auto phase=-B::from_strings("0","1")*pi;
  auto scalar=FrobeniusCrossing::prepare({{e("eps/x")}},0,1,4,4);
  Boundary seed{{B(1),B(0),B(0),B(0),B(0)}};
  auto value=scalar.upper_bank(seed,Rational("1/2"),Rational("1/2"));
  B expected(1);for(unsigned k=0;k<=4;++k){near(value[0][k],expected);expected=expected*phase/B(k+1);}
  auto jordan=FrobeniusCrossing::prepare({{e("eps/x"),e("1/x")},{e("0"),e("eps/x")}},0,1,4,4);
  Boundary jordan_seed(2,std::vector<B>(5,B(0)));jordan_seed[1][0]=B(1);
  auto ladder=jordan.upper_bank(jordan_seed,Rational("1/2"),Rational("1/2"));
  expected=B(1);for(unsigned k=0;k<=4;++k){near(ladder[1][k],expected);near(ladder[0][k],phase*expected);expected=expected*phase/B(k+1);}
  auto resonance=FrobeniusCrossing::prepare({{e("0"),e("0")},{e("1"),e("1/x")}},0,1,4,0);
  auto resonant=resonance.upper_bank({{B(1)},{B(0)}},Rational("1/2"),Rational("1/2"));
  near(resonant[0][0],B(1));near(resonant[1][0],phase/B(2));
  residual(scalar,{{"eps/x"}});
  residual(jordan,{{"eps/x","1/x"},{"0","eps/x"}});
  residual(resonance,{{"0","0"},{"1","1/x"}});
  check(!FrobeniusCrossing::omitted_tail_certified,"crossing invented omitted-tail certificate");
  std::cout<<"Scalar epsilon branch, Jordan log ladder and positive resonance crossing passed\n";
  if(argc>1)return run_original_banana_equal_real(argv[1]);
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
