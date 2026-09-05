#include "diffexp/fire_session.hpp"
#include "diffexp/level_preparation.hpp"
#include "diffexp/families.hpp"
#include <iostream>
using namespace diffexp;
void check(bool value,const char* why){if(!value)throw std::runtime_error(why);}
int main(int argc,char** argv){try {
  std::cout.setf(std::ios::unitbuf);
  if(argc<2){std::cout<<"FIRE session integration requires an executable\n";return 0;}
  ExactField field({"x","eps","I"});Exact x(field,"x"),eps(field,"eps");auto d=eps.constant(4)-eps.constant(2)*eps;
  auto raw=ibp::quadratic_family(feynman::banana(1,{Rational(1),Rational(2)}),x);ibp::PropagatorBasis basis(raw);
  fire::Options options;options.executable=argv[1];fire::Session session(basis,d,field,options);
  auto first=session({{2,1}},options);check(first.success,first.reason.c_str());
  auto second=session({{2,1},{4,1}},options);check(second.success,second.reason.c_str());
  auto cold=fire::reduce(basis,d,field,{{2,1},{4,1}},options);check(cold.success,cold.reason.c_str());
  check(second.reductions==cold.reductions,"resumed exact rules differ from cold bubble");
  auto same=session({{2,1},{4,1}},options);check(same.success&&session.runs()==2&&session.cache_hits()==1,"same demand must reuse completed result");
  auto reduced_limits=options;reduced_limits.timeout_seconds=1;reduced_limits.memory_bytes=1024;
  check(session({{2,1},{4,1}},reduced_limits).success,"shrinking resource budgets do not invalidate exact results");
  auto shrink=session({{2,1}},options);check(!shrink.success,"shrinking demand must reject");
  auto wrong=options;wrong.zero_sectors={{1,-1}};check(!session({{2,1},{4,1}},wrong).success,"changed scientific zero sectors must reject");
  std::cout<<"FIRE session exact growing-demand bubble passed\n";
  if(argc>2) {
    auto banana=ibp::quadratic_family(feynman::banana(4,{Rational(1),Rational(1),Rational(1),Rational(1),Rational(1)}),x);
    ibp::PropagatorBasis b(ibp::merge(banana,0,1,x));ibp::Integral target(14,0);target[0]=2;target[1]=target[2]=target[3]=1;
    level::Options limits;limits.provider=options;limits.provider.timeout_seconds=60;limits.total_timeout_seconds=80;
    fire::Session shared(b,d,field,options);
    auto warm=level::prepare(b,d,field,0,{target},limits,[&](const auto& demand,const auto& resource){return shared(demand,resource);});
    std::cout<<"Banana4 rules-session success="<<warm.success<<" passes="<<warm.passes<<" masters="<<warm.ordered_basis.size()<<" seconds="<<warm.elapsed_seconds<<" reason="<<warm.reason<<"\n";
    for(const auto& path:warm.fire_directories)std::cout<<path<<"\n";
    check(warm.success,"session Banana4 closure failed");
    auto baseline=level::prepare(b,d,field,0,{target},limits,[&](const auto& demand,const auto& resource){return fire::reduce(b,d,field,demand,resource);});
    std::cout<<"Banana4 cold success="<<baseline.success<<" passes="<<baseline.passes<<" masters="<<baseline.ordered_basis.size()<<" seconds="<<baseline.elapsed_seconds<<" reason="<<baseline.reason<<"\n";
    check(baseline.success,"cold Banana4 closure failed");
    check(warm.ordered_basis==baseline.ordered_basis&&warm.matrix==baseline.matrix&&warm.target_rows==baseline.target_rows,"session and cold exact closure mismatch");
  }
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
