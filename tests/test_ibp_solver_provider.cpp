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
  // One complete elimination at a fresh prime and point anchors every validation
  // attempt; the remaining held-out points eliminate the pilot's dependency subset.
  auto anchored_options=selected_options;anchored_options.cache_directory=options.cache_directory/"anchored";anchored_options.complete_validation_anchor=true;
  ibp_solver::Session anchored_session(basis,d,field,anchored_options);auto anchored_result=anchored_session(requests,limits);
  require(anchored_result.success&&anchored_result.reductions==result.reductions,"anchored subset validation changed exact reconstruction");
  require(anchored_session.statistics().full_solves==2&&anchored_session.statistics().selected_solves>=8,"anchored validation must run the complete pilot and exactly one complete held-out elimination");
  require(anchored_result.directory!=selected_result.directory&&anchored_result.directory!=sliced_result.directory,"anchored validation shares another protocol's completion identity");
  ibp_solver::Sampler anchor_points(basis,d,selected_options);const auto first_prime=ibp_solver::prime(1),second_prime=ibp_solver::prime(2);
  auto pilot_rows=anchor_points(requests,{7,12352},first_prime,limits,true),subset_rows=anchor_points(requests,{9,12354},second_prime,limits,true);
  auto complete_rows=anchor_points(requests,{9,12354},second_prime,limits,true,true),later_rows=anchor_points(requests,{11,12356},first_prime,limits,true);
  require(pilot_rows.success&&subset_rows.success&&complete_rows.success&&later_rows.success&&subset_rows.reductions==complete_rows.reductions,"complete anchor disagrees with the dependency subset at a generic point");
  require(anchor_points.statistics().full_solves==2&&anchor_points.statistics().selected_solves==2&&anchor_points.statistics().selection_fallbacks==0,"complete anchor was not a full elimination or replaced the pilot's subset");
  // A trace that cannot fit its node budget costs only the replay shortcut: the
  // primal reduction of the same probe is valid and sampling continues with it.
  auto untraced_options=options;untraced_options.cache_directory=options.cache_directory/"untraced";untraced_options.max_trace_nodes=4;
  ibp_solver::Session untraced(basis,d,field,untraced_options);auto untraced_result=untraced(requests,limits);
  require(untraced_result.success,untraced_result.reason.c_str());require(untraced_result.reductions==result.reductions,"primal-only sampling changed the reconstruction");
  require(untraced.statistics().trace_disabled>=1&&untraced.statistics().trace_replays==0&&untraced.statistics().full_solves==untraced.statistics().probes,"trace budget failure was not confined to the replay shortcut");
  // Ordinary sampling never returns to an earlier prime, so its trace is released.
  ibp_solver::Sampler released(basis,d,options),unreleased(basis,d,options);
  for(modular::Word value:{21UL,22UL,23UL})require(released(requests,{value,400+value},first_prime,limits).success,"trace warm-up probe");
  const auto warm=released.statistics();require(warm.full_solves==2&&warm.trace_replays==2,"trace was not learned, checked once and then replayed");
  require(released(requests,{24,424},second_prime,limits).success,"second-prime probe");
  auto returned=released(requests,{25,425},first_prime,limits),expected_rows=unreleased(requests,{25,425},first_prime,limits,true);
  require(returned.success&&expected_rows.success&&returned.reductions==expected_rows.reductions,"relearned trace disagrees with fresh elimination");
  require(released.statistics().trace_replays==warm.trace_replays&&released.statistics().full_solves==warm.full_solves+2,"a finished prime kept its arithmetic trace");
  // The session owns its exact field; a caller's temporary handle may die first.
  {auto owned_options=sliced_options;owned_options.cache_directory=options.cache_directory/"owned-field";owned_options.closure_parameter=0;
    auto owned=std::make_unique<ibp_solver::Session>(basis,d,ExactField(field),owned_options);level::Options stage_options;stage_options.total_timeout_seconds=90;stage_options.provider=limits;
    auto stage=level::prepare(basis,d,field,0,requests,stage_options,[&](const auto& r,const auto& l){return (*owned)(r,l);});require(stage.success,stage.reason.c_str());}
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
  // An unfinished plan is still a valid hint: the level preparer closes the rest.
  {auto short_options=sliced_options;short_options.cache_directory=options.cache_directory/"short-plan";short_options.closure_parameter=0;short_options.max_planning_passes=1;
    ibp_solver::Session shortened(basis,d,field,short_options);level::Options stage_options;stage_options.total_timeout_seconds=90;stage_options.provider=limits;
    auto stage=level::prepare(basis,d,field,0,requests,stage_options,[&](const auto& r,const auto& l){return shortened(r,l);});
    require(stage.success,stage.reason.c_str());require(shortened.planning_statistics().passes>=1&&!std::filesystem::exists(short_options.cache_directory/"closure-plans"),"an unstabilized closure plan was saved as complete");}
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
  // An exceptional point is replaced by the next ordinal: a pole of an input
  // coefficient, or a point where the terminal integrals differ, including the
  // very first one. The choice is deterministic, bounded and recorded.
  fire_modular::Options skip_options;skip_options.cache_directory=options.cache_directory/"exceptional";skip_options.modulus=ibp_solver::prime;skip_options.sparse_lifting=true;
  skip_options.provider_identity="exceptional-point-regression";
  const auto moved_point=modular::point(2,1,0),pole_point=modular::point(2,1,4),second_pole=modular::point(2,2,1);
  skip_options.sample_provider=[&](const auto& r,const auto& point,auto p,const auto& l){
    if((p==first_prime&&point==pole_point)||(p==second_prime&&point==second_pole)){fire::Result pole;pole.exceptional_point=true;pole.reason="deliberate input pole";return pole;}
    auto reduced=fallback_sampler(r,point,p,l);
    if(reduced.success&&p==first_prime&&point==moved_point){diffexp::ibp::Relation moved{{diffexp::ibp::Integral{3,3},d.constant(1)}};reduced.reductions[diffexp::ibp::Integral{3,3}]=moved;reduced.reductions[requests[0]]=moved;}
    return reduced;
  };
  skip_options.validation_provider=[&](const auto& r,const auto& point,auto p,const auto& l){return fallback_sampler(r,point,p,l,true);};
  auto skip_run=[&](const fire_modular::Options& o){fire_modular::Session session(basis,d,field,o);return session(requests,limits);};
  auto skipped_result=skip_run(skip_options);require(skipped_result.success,skipped_result.reason.c_str());
  require(skipped_result.reductions==result.reductions,"exceptional sample points changed the reconstruction");
  auto skip_record=[&]{return fire_modular::detail::read(skipped_result.directory/"completed.json").at("skipped_samples").as_array();};
  const auto first_record=skip_record();
  auto skipped_at=[&](std::size_t i){const auto& entry=first_record.at(i).as_object();return std::pair{entry.at("prime_index").to_number<unsigned>(),entry.at("ordinal").to_number<std::size_t>()};};
  require(first_record.size()==3&&skipped_at(0)==std::pair{1u,std::size_t(0)}&&skipped_at(1)==std::pair{1u,std::size_t(4)}&&skipped_at(2)==std::pair{2u,std::size_t(1)},"exceptional samples were not skipped exactly once each");
  std::filesystem::remove(skipped_result.directory/"completed.json");auto skipped_again=skip_run(skip_options);
  require(skipped_again.success&&skipped_again.reductions==result.reductions&&skip_record()==first_record,"resumed reconstruction chose a different sample point set");
  auto bounded_skip=skip_options;bounded_skip.max_skipped_samples=1;bounded_skip.provider_identity="exceptional-point-budget";auto exhausted=skip_run(bounded_skip);
  require(!exhausted.success&&exhausted.reason.find("exceptional-sample budget")!=std::string::npos,"unbounded exceptional-sample skipping");
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
