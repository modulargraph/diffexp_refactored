#include "diffexp/ibp_solver_provider.hpp"
#include "diffexp/families.hpp"
#include <iostream>
using namespace diffexp;
void require(bool v,const char* s){if(!v)throw std::runtime_error(s);}
int main(int argc,char** argv){try{
  ExactField field({"x","d"});Exact x(field,"x"),d(field,"d");
  auto raw=diffexp::ibp::quadratic_family(feynman::banana(1,{Rational(1),Rational(2)}),d);
  raw.physical[0].constant=x;diffexp::ibp::PropagatorBasis basis(raw);
  auto name=(std::filesystem::temp_directory_path()/"diffexp-ibp-provider-XXXXXX").string();std::vector<char> temp(name.begin(),name.end());temp.push_back(0);require(mkdtemp(temp.data()),"private test cache");
  ibp_solver::Options options;options.cache_directory=temp.data();options.max_degree=8;
  auto malformed=basis;malformed.denominators.pop_back();bool rejected=false;try{ibp_solver::Sampler invalid(malformed,d,options);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"incomplete denominator basis was accepted");
  ibp_solver::Session session(basis,d,field,options);fire::Options limits;limits.timeout_seconds=90;
  std::vector<diffexp::ibp::Integral> requests{{2,1},{1,2}};auto result=session(requests,limits);require(result.success,result.reason.c_str());require(session.statistics().probes>0,"upstream solver was not called");
  require(session.statistics().trace_replays>0&&session.statistics().full_solves<session.statistics().probes,"repeated probes did not use arithmetic replay");
  require(session.statistics().templates==1&&session.statistics().full_solves>=10,"template reuse and independent held-out full reductions");
  // A sample-specific massless cancellation must trigger guarded relearning.
  ibp_solver::Sampler guarded(basis,d,options),independent(basis,d,options);
  auto special=guarded(requests,{0,12345},ibp_solver::prime(1),limits);require(special.success,special.reason.c_str());
  auto general=guarded(requests,{7,23456},ibp_solver::prime(1),limits);
  auto reference=independent(requests,{7,23456},ibp_solver::prime(1),limits,true);
  require(general.success&&reference.success&&general.reductions==reference.reductions,"exceptional point relearning disagrees with full reduction");
  require(guarded.statistics().trace_fallbacks>0,"sample-specific cancellation was not guarded");
  diffexp::ibp::Generator generator(basis,d);diffexp::ibp::ExactReducer exact(d);
  for(int a=-1;a<=4;++a)for(int b=-1;b<=4;++b)for(auto row:generator.relations({a,b}))exact.insert(std::move(row));
  for(const auto& [a,row]:result.reductions){diffexp::ibp::Relation residual{{a,d.constant(1)}};diffexp::ibp::add_scaled(residual,row,d.constant(-1));require(exact.reduce(residual).remainder.empty(),"reconstruction differs from independent exact IBPs");}
  ibp_solver::Session resumed(basis,d,field,options);auto recovered=resumed(requests,limits);require(recovered.success&&recovered.reductions==result.reductions,"cache recovery");require(resumed.statistics().probes==0,"completed cache reran probes");
  auto completed=result.directory/"completed.json";auto record=fire_modular::detail::read(completed);auto changed=result.reductions;changed.at(requests[0]).begin()->second=changed.at(requests[0]).begin()->second+d.constant(1);record["tables"]=fire_modular::detail::table(changed,fire::SymbolMap(field.variables()));fire_modular::detail::save(completed,record);
  ibp_solver::Session invalid(basis,d,field,options);require(!invalid(requests,limits).success,"corrupt reconstruction was accepted");
  require(ibp_solver::prime(1)==2305843009213693951ULL&&ibp_solver::prime(2)<ibp_solver::prime(1),"61-bit default and distinct primes");
  std::filesystem::remove_all(options.cache_directory);
  if(argc>1){
    ExactField one_field({"d"});Exact dim(one_field,"d");auto example=feynman::example_family("double_box_planar");
    diffexp::ibp::PropagatorBasis db(diffexp::ibp::quadratic_family(example.momenta,dim));
    diffexp::ibp::Integral top(db.denominators.size(),0);std::fill_n(top.begin(),db.physical_count,1);
    std::vector<diffexp::ibp::Integral> targets{top};for(unsigned k=0;k<top.size();++k){auto a=top;a[k]+=k<db.physical_count?1:-1;targets.push_back(a);}
    ibp_solver::Sampler sampler(db,dim,{});auto first=sampler(targets,{1234567},ibp_solver::prime(1),limits);require(first.success,first.reason.c_str());
    std::set<diffexp::ibp::Integral> all(targets.begin(),targets.end());for(const auto& [a,row]:first.reductions){all.insert(a);for(const auto& [b,c]:row)all.insert(b);}
    auto fire_limits=limits;fire_limits.executable=argv[1];auto reference=fire::reduce(db,dim,one_field,{all.begin(),all.end()},fire_limits);require(reference.success,reference.reason.c_str());
    for(unsigned pi=1;pi<=2;++pi)for(modular::Word at:{1234567UL,98765431UL}){
      auto p=ibp_solver::prime(pi);::ibp::Field finite(p);auto native=sampler(targets,{at},p,limits);require(native.success,native.reason.c_str());
      for(const auto& a:targets){std::map<diffexp::ibp::Integral,modular::Word> residual;
        for(const auto& [b,c]:reference.reductions.at(a))residual[b]=modular::evaluate(c,{at},p);
        for(const auto& [b,c]:native.reductions.at(a))for(const auto& [m,q]:reference.reductions.at(b))
          residual[m]=finite.sub(residual[m],finite.mul(modular::evaluate(c,{at},p),modular::evaluate(q,{at},p)));
        for(const auto& [m,c]:residual)require(c==0,"double-box basis conversion disagrees with exact FIRE");
      }
    }
    std::cout<<"Double-box ten targets agree with exact FIRE after basis conversion at two dimensions and two 61-bit primes\n";
  }
  std::cout<<"IBP solver x/d reconstruction, exact identities, retained basis, cache recovery and corruption rejection passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
