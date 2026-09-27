#include "diffexp/ibp_solver_provider.hpp"
#include "diffexp/families.hpp"
#include "diffexp/level_preparation.hpp"
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
  auto sliced_options=options;sliced_options.degree_slices=true;
  ibp_solver::Session sliced(basis,d,field,sliced_options);auto sliced_result=sliced(requests,limits);
  require(sliced_result.success,sliced_result.reason.c_str());
  require(sliced_result.reductions==result.reductions,"coordinate slice reconstruction differs from dense reconstruction");
  require(sliced_result.directory!=result.directory,"slice sample convention shares dense cache identity");
  ibp_solver::Session sliced_resume(basis,d,field,sliced_options);auto sliced_recovered=sliced_resume(requests,limits);
  require(sliced_recovered.success&&sliced_recovered.reductions==result.reductions&&sliced_resume.statistics().probes==0,"slice reconstruction cache recovery");
  // Saving fewer ordinary probes changes only restart cost. Every independent
  // validation probe remains durable, and completion is checked on reopening.
  auto compact_options=sliced_options;compact_options.cache_directory=options.cache_directory/"compact";compact_options.sample_checkpoint_interval=16;
  ibp_solver::Session compact(basis,d,field,compact_options);auto compact_result=compact(requests,limits);
  require(compact_result.success&&compact_result.reductions==result.reductions,"compact checkpoints changed reconstruction");
  unsigned saved=0,validation=0;
  for(const auto& entry:std::filesystem::directory_iterator(compact_result.directory)){
    auto name=entry.path().filename().string();if(!name.starts_with("sample-"))continue;
    auto ordinal=std::stoul(name.substr(name.find_last_of('-')+1));++saved;
    if(ordinal>=100000)++validation;else require(ordinal%16==0,"checkpoint interval was ignored");
  }
  require(validation==6&&saved<compact.statistics().probes,"compact checkpoints omitted validation or retained every probe");
  ibp_solver::Session compact_resume(basis,d,field,compact_options);auto compact_recovered=compact_resume(requests,limits);
  require(compact_recovered.success&&compact_recovered.reductions==result.reductions&&compact_resume.statistics().probes==0,"compact completed cache recovery failed");
  std::filesystem::remove(compact_result.directory/"completed.json");
  ibp_solver::Session compact_interrupted(basis,d,field,compact_options);auto compact_rebuilt=compact_interrupted(requests,limits);
  require(compact_rebuilt.success&&compact_rebuilt.reductions==result.reductions,"compact interrupted reconstruction recovery failed");
  // Independent validation can eliminate a dependency subset of ORIGINAL
  // equations. It must still compare complete rows at fresh points/primes.
  auto selected_options=sliced_options;selected_options.cache_directory=options.cache_directory/"selected";selected_options.selected_source_rows=true;
  ibp_solver::Session selected_session(basis,d,field,selected_options);auto selected_result=selected_session(requests,limits);
  require(selected_result.success&&selected_result.reductions==result.reductions,"original-row subset validation changed exact reconstruction");
  require(selected_session.statistics().full_solves==1&&selected_session.statistics().selected_solves>=9,"subset validation repeated full elimination or skipped fresh primal solves");
  require(selected_result.directory!=sliced_result.directory,"different validation protocols share completion identity");
  unsigned selected_validation=0;for(const auto& entry:std::filesystem::directory_iterator(selected_result.directory)){
    auto name=entry.path().filename().string();if(name.starts_with("sample-")&&std::stoul(name.substr(name.find_last_of('-')+1))>=100000)++selected_validation;
  }require(selected_validation==6,"subset validation omitted fresh held-out probes");
  ibp_solver::Sampler selected_points(basis,d,selected_options),full_points(basis,d,options);
  for(unsigned pi:{1u,2u,1u})for(modular::Word value:{7UL,0UL,19UL}){
    auto point=std::vector<modular::Word>{value,12345+value};
    auto subset=selected_points(requests,point,ibp_solver::prime(pi),limits,true);
    auto complete=full_points(requests,point,ibp_solver::prime(pi),limits,true);
    require(subset.success&&complete.success&&subset.reductions==complete.reductions,"fresh subset/fallback retained stale point or prime bindings");
  }
  require(selected_points.statistics().selected_solves>0&&selected_points.statistics().selection_fallbacks>0,"exceptional subset support did not trigger complete fallback");
  // Plan derivatives with finite-field probes, then reconstruct one combined
  // batch. Exact closure and a fresh-process completion replay remain required.
  for(bool sector_limited:{false,true}){
    auto planned_options=sliced_options;planned_options.cache_directory=options.cache_directory/(sector_limited?"planned-sectors":"planned-global");planned_options.closure_parameter=0;planned_options.target_sector_seeds=sector_limited;planned_options.selected_source_rows=sector_limited;
    ibp_solver::Session planned(basis,d,field,planned_options);level::Options stage_options;stage_options.total_timeout_seconds=90;stage_options.provider=limits;
    auto stage=level::prepare(basis,d,field,0,requests,stage_options,[&](const auto& r,const auto& l){return planned(r,l);});
    require(stage.success,stage.reason.c_str());require(stage.passes==1&&planned.planning_statistics().probes>=2,"planned derivative closure reconstructed more than one batch");
    ibp_solver::Session reopened(basis,d,field,planned_options);
    auto recovered_stage=level::prepare(basis,d,field,0,requests,stage_options,[&](const auto& r,const auto& l){return reopened(r,l);});
    require(recovered_stage.success&&recovered_stage.matrix==stage.matrix&&recovered_stage.target_rows==stage.target_rows&&reopened.statistics().probes==0&&reopened.planning_statistics().cache_hits==1,"planned closure did not recover without fresh probes");
    for(const auto& a:requests){auto direct=std::find(stage.requested_integrals.begin(),stage.requested_integrals.end(),a)-stage.requested_integrals.begin();diffexp::ibp::Relation residual{{a,d.constant(1)}};for(std::size_t j=0;j<stage.ordered_basis.size();++j)diffexp::ibp::add(residual,stage.ordered_basis[j],-stage.target_rows[direct][j]);
      diffexp::ibp::ExactReducer check(d);for(const auto& [b,row]:result.reductions){diffexp::ibp::Relation equation{{b,d.constant(1)}};diffexp::ibp::add_scaled(equation,row,d.constant(-1));check.insert(equation);}require(check.reduce(residual).remainder.empty(),"planned closure differs from unplanned exact reductions");}
    auto manifest=*std::filesystem::directory_iterator(planned_options.cache_directory/"closure-plans");auto corrupted=fire_modular::detail::read(manifest.path());corrupted["demands"]=boost::json::array{};fire_modular::detail::save(manifest.path(),corrupted);
    ibp_solver::Session invalid_hint(basis,d,field,planned_options);auto invalid_stage=level::prepare(basis,d,field,0,requests,stage_options,[&](const auto& r,const auto& l){return invalid_hint(r,l);});require(!invalid_stage.success,"valid-hash closure plan omitted original requests");
  }
  // A failed specialized probe must only disable the degree heuristic. This
  // deliberately exercises fallback without weakening cache-integrity checks.
  ibp_solver::Sampler fallback_sampler(basis,d,options);bool failed_slice=false;
  fire_modular::Options fallback_options;fallback_options.cache_directory=options.cache_directory;
  fallback_options.modulus=ibp_solver::prime;fallback_options.degree_slices=true;fallback_options.sparse_lifting=true;
  fallback_options.provider_identity="slice-failure-regression";
  auto anchor=modular::point(2,1,50000);
  fallback_options.sample_provider=[&](const auto& r,const auto& point,auto p,const auto& l){
    if(point[1]==anchor[1]){failed_slice=true;fire::Result rejected;rejected.reason="deliberately unusable slice";return rejected;}
    return fallback_sampler(r,point,p,l);
  };
  fallback_options.validation_provider=[&](const auto& r,const auto& point,auto p,const auto& l){return fallback_sampler(r,point,p,l,true);};
  fire_modular::Session fallback(basis,d,field,fallback_options);auto fallback_result=fallback(requests,limits);
  require(failed_slice&&fallback_result.success&&fallback_result.reductions==result.reductions,"failed coordinate probe did not fall back safely");
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
