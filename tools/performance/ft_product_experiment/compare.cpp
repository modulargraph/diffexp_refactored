#include "product_recurrence.hpp"
#include <chrono>
#include <iostream>
using namespace diffexp;
using namespace experiment;
void check(bool yes,const char* why){if(!yes)throw std::runtime_error(why);}
bool contains(const B& b,const Rational& q) {
  fmpq_t a;fmpq_init(a);fmpq_set_str(a,q.str().c_str(),10);fmpq_canonicalise(a);
  bool yes=arb_contains_fmpq(acb_realref(b.raw()),a)&&arb_contains_zero(acb_imagref(b.raw()));
  fmpq_clear(a);return yes;
}
int main(int argc,char** argv) {
  if(argc!=2||std::string(argv[1])!="--run") {
    std::cerr<<"Use --run for bounded validation and timing.\n";return 2;
  }
  try {
    B::set_precision(256);
    ExactField field({"x","eps","I"});
    const unsigned d=10,high=3,N=128,small=12;
    std::vector<RationalLineEntry> entries;
    for(unsigned i=0;i<d;++i) {
      entries.push_back({i,i,0,Exact(field,"1/(2-x)")});
      for(unsigned j=0;j<i;++j)
        entries.push_back({i,j,j%2,Exact(field,"1/("+std::to_string(j+3)+"-x)")});
    }
    auto start=std::chrono::steady_clock::now();
    auto c=prepare(entries,d,N,high);
    double prep=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    check(c.supported,"manufactured input unsupported");
    std::vector<std::vector<Rational>> exact(d,std::vector<Rational>(high+1,Rational(0)));
    Boundary boundary(d,std::vector<B>(high+1,B(0)));
    for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=high;++k) {
      exact[i][k]=Rational(long(i+k+1));boundary[i][k]=B(i+k+1);
    }
    auto a=coefficients<Rational>(c,exact,small),b=coefficients<Rational>(c,exact,small,true);
    check(a==b,"exact coefficient mismatch against independent rational convolution");
    // Carried real uncertainty: exact interior boundary realizations must remain
    // enclosed by the retained-polynomial result, including epsilon mixing.
    for(auto& r:boundary)for(auto& v:r)arb_add_error_2exp_si(acb_realref(v.raw()),-80);
    B h=polynomial_transport::detail::ball(Rational("1/8"));
    auto uncertain=chart(c,boundary,B(0),h,small);
    for(int sign:{-1,0,1}) {
      auto sample=exact;
      for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=high;++k)
        sample[i][k]+=Rational(sign*((i+k)%2?1:-1))*Rational("1/2417851639229258349412352"); // 2^-81
      auto oracle=horner(coefficients<Rational>(c,sample,small,true),Rational("1/8"));
      for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=high;++k)
        check(contains(uncertain[i][k],oracle[i][k]),"lost carried uncertainty");
    }
    // Compare unchanged production row-LCM recurrence, with the same balls/N.
    start=std::chrono::steady_clock::now();
    auto baseline=polynomial_transport::chart(c.baseline,boundary,B(0),h,N);
    double base_first=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    start=std::chrono::steady_clock::now();
    auto candidate=chart(c,boundary,B(0),h,N);
    double candidate_first=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=high;++k)
      check(acb_overlaps(baseline[i][k].raw(),candidate[i][k].raw()),"baseline overlap");
    bool fallback=false;
    auto shifted=chart(c,boundary,h,h,small,&fallback);
    auto shifted_ref=polynomial_transport::chart(c.baseline,boundary,h,h,small);
    check(fallback,"nonzero-center fallback missing");
    for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=high;++k)
      check(acb_overlaps(shifted[i][k].raw(),shifted_ref[i][k].raw()),"fallback mismatch");
    auto moving=prepare({{0,0,0,Exact(field,"1/(2-x+eps)")}},1,small,high);
    check(!moving.supported,"unexpanded epsilon should use baseline");
    chart(moving,{boundary[0]},B(0),h,small,&fallback);check(fallback,"epsilon fallback missing");
    std::size_t product_ops=0,dense_ops=0;
    coefficients<B>(c,boundary,N,false,&product_ops);
    coefficients<B>(c,boundary,N,true,&dense_ops);
    auto timed=[&](auto fn) {
      auto t=std::chrono::steady_clock::now();
      for(unsigned repeat=0;repeat<5;++repeat)fn();
      return std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()/5;
    };
    double base_seconds=timed([&]{polynomial_transport::chart(c.baseline,boundary,B(0),h,N);});
    double candidate_seconds=timed([&]{chart(c,boundary,B(0),h,N);});
    std::cout<<"{\"precision_bits\":256,\"taylor_order\":128,\"epsilon_high\":3,\"dimension\":10,"
      <<"\"edges\":"<<entries.size()<<",\"compile_seconds\":"<<prep
      <<",\"first_full_order_baseline_seconds\":"<<base_first
      <<",\"first_full_order_candidate_seconds\":"<<candidate_first
      <<",\"warm_repetitions\":5,\"baseline_chart_seconds\":"<<base_seconds
      <<",\"candidate_chart_seconds\":"<<candidate_seconds
      <<",\"product_multiply_adds\":"<<product_ops<<",\"dense_multiply_adds\":"<<dense_ops
      <<",\"baseline_fallback_rows\":"<<c.baseline.fallback_rows
      <<",\"validation\":\"passed\",\"omitted_tails_certified\":false}\n";
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
