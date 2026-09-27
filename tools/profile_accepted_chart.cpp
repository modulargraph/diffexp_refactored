// Replay one accepted chart from immutable local snapshots, retaining all guards.
#include "diffexp/adjoint_checkpoint.hpp"
#include <chrono>
#include <iostream>
#include <cstring>
#include <sys/mman.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <signal.h>
using namespace diffexp;namespace json=boost::json;using B=Jet::Ball;
using Clock=std::chrono::steady_clock;
struct Shared {
  AdjointConditioningStats stats;
  unsigned status=0,leg=0,observed_charts=0;double from=0,to=0,output_radius=-INFINITY;
  unsigned long long finite=1,disjoint=0;unsigned comparison=0;double midpoint_distance=0,baseline_width=0,result_width=0;char error[1024]{};
};
struct Accepted {};
static std::string bytes(const std::filesystem::path& p){
  if(std::filesystem::file_size(p)>64*1024*1024)throw std::length_error("snapshot file budget");
  std::ifstream f(p);if(!f)throw std::runtime_error("snapshot open failed");return {std::istreambuf_iterator<char>(f),{}};
}
static json::object counters(const AdjointConditioningStats& s){json::object out;
#define FIELD(name) out[#name]=s.name
  FIELD(polynomial_charts);FIELD(rational_cross_checks);FIELD(rational_compilations);
  FIELD(centered_charts);FIELD(homogeneous_chart_maps);FIELD(centered_budget_skips);FIELD(conditioning_subdivisions);
  FIELD(polynomial_homogeneous_columns);FIELD(rational_homogeneous_columns);FIELD(polynomial_midpoint_charts);FIELD(rational_midpoint_charts);
  FIELD(reduced_precision_homogeneous_maps);FIELD(full_precision_map_retries);FIELD(centered_before_rational_charts);FIELD(centered_first_fallbacks);
  FIELD(centered_only_charts);FIELD(centered_only_fallbacks);FIELD(exact_input_map_columns_skipped);
  FIELD(max_homogeneous_map_workers);FIELD(circuit_charts);FIELD(circuit_homogeneous_columns);FIELD(circuit_preparations);
  FIELD(circuit_addmul_operations);FIELD(circuit_scalar_operations);FIELD(circuit_peak_live_cells);
  FIELD(compilation_seconds);FIELD(preparation_seconds);FIELD(source_seconds);FIELD(homogeneous_seconds);FIELD(reference_seconds);
#undef FIELD
  return out;
}
int main(int argc,char** argv){try{
  if(argc<4)throw std::invalid_argument("usage: profile_accepted_chart SNAPSHOT_PROBLEM CHECKPOINT|initial SECONDS(1..60) [--seek-fallback PREFIX] [--compact-recovery true|false] [--output-prefix PREFIX] [--compare-after V6_CHECKPOINT] [--complete] [--centered-only true|false] [--grouped-dot true|false]");
  bool seek=false,recovery=true,centered_only=true,grouped_dot=true,complete=false;std::string output_prefix,compare_after;
  for(int i=4;i<argc;++i){std::string flag=argv[i];if(flag=="--complete"){complete=true;continue;}if(i+1==argc)throw std::invalid_argument("missing option value");std::string value=argv[++i];
    if(flag=="--seek-fallback"){seek=true;output_prefix=value;}else if(flag=="--output-prefix")output_prefix=value;else if(flag=="--compare-after")compare_after=value;
    else if(flag=="--centered-only"&&(value=="true"||value=="false"))centered_only=value=="true";
    else if(flag=="--grouped-dot"&&(value=="true"||value=="false"))grouped_dot=value=="true";
    else if(flag=="--compact-recovery"&&(value=="true"||value=="false"))recovery=value=="true";
    else if(flag=="--compact-only"&&(value=="true"||value=="false"))centered_only=value=="true";else throw std::invalid_argument("invalid profiler option");}
  if(complete&&seek)throw std::invalid_argument("complete and seek modes are exclusive");
  const bool from_initial=std::string(argv[2])=="initial";
  unsigned seconds=std::stoul(argv[3]);if(!seconds||seconds>60)throw std::invalid_argument("wall budget must be 1..60 seconds");
  auto problem_bytes=bytes(argv[1]),checkpoint_bytes=from_initial?std::string{}:bytes(argv[2]);auto object=json::parse(problem_bytes).as_object();
  if(object.at("schema")!="DiffExp.FTTransportProblem/v1")throw std::invalid_argument("capture schema");
  auto precision=artifacts::detail::integer(object.at("working_bits"));if(precision<64||precision>1000000)throw std::invalid_argument("precision budget");B::set_precision(precision);
  std::vector<std::string> names;for(const auto& n:object.at("variables").as_array())names.emplace_back(n.as_string());ExactField field(names);
  auto matrix=[&](const json::value& value){ExactEpsilonMatrix result;for(const auto& row:value.as_array()){result.emplace_back();for(const auto& e:row.as_array())result.back().emplace_back(field,std::string(e.as_string()));}return result;};
  auto connection=matrix(object.at("matrix")),forcing=matrix(object.at("forcing"));auto initial=numerical_rows_io::read_rows(object.at("initial"));
  std::vector<Exact> path;for(const auto& v:object.at("path").as_array())path.emplace_back(field,std::string(v.as_string()));
  const auto& settings=object.at("settings").as_object();AdjointOptions options;
  if(settings.contains("ordinary_order"))options.taylor_order=artifacts::detail::integer(settings.at("ordinary_order"));
  if(settings.contains("centered_before_rational"))options.centered_before_rational=settings.at("centered_before_rational").as_bool();
  if(settings.contains("centered_map_working_bits"))options.centered_map_working_bits=artifacts::detail::integer(settings.at("centered_map_working_bits"));
  if(settings.contains("centered_map_workers"))options.centered_map_workers=artifacts::detail::integer(settings.at("centered_map_workers"));
  if(settings.contains("rational_circuit_recurrence"))options.rational_circuit_recurrence=settings.at("rational_circuit_recurrence").as_bool();
  if(!options.rational_circuit_recurrence)throw std::invalid_argument("expected current circuit checkpoint identity");
  options.compact_centered_recovery=false;options.compact_centered_only=false;options.circuit_grouped_dot=false; // Explicitly verify imported v6 scientific input.
  auto key=adjoint_checkpoint::identity(connection,initial,forcing,path,options);
  auto saved=from_initial?AdjointContinuation{0,0,0,initial}:adjoint_checkpoint::detail::decode(json::parse(checkpoint_bytes),key);
  options.compact_centered_recovery=recovery;options.compact_centered_only=centered_only;options.circuit_grouped_dot=grouped_dot;
  const auto execution_key=adjoint_checkpoint::identity(connection,initial,forcing,path,options);
  std::optional<AdjointContinuation> baseline_after;
  if(!compare_after.empty()) {
    auto value=json::parse(bytes(compare_after));
    const auto comparison_identity=artifacts::detail::string(value.as_object().at("payload").as_object().at("identity"));
    bool matches=false;
    for(unsigned flags=0;flags<8;++flags){auto candidate=options;candidate.compact_centered_recovery=flags&1;candidate.compact_centered_only=flags&2;candidate.circuit_grouped_dot=flags&4;
      matches|=comparison_identity==adjoint_checkpoint::identity(connection,initial,forcing,path,candidate);}
    if(!matches)throw std::invalid_argument("comparison checkpoint does not match the same scientific inputs under an admitted policy");
    baseline_after=adjoint_checkpoint::detail::decode(value,comparison_identity);
  }
  auto start_leg=saved.leg;double start_parameter=saved.parameter;if(start_parameter==1){++start_leg;start_parameter=0;}
  if(start_leg+1>=path.size())throw std::invalid_argument("saved checkpoint already completed entire arm; cannot profile a next chart");
  json::object report{{"schema","DiffExp21.AcceptedChartProfile/v1"},{"status","starting"},
    {"problem_snapshot",argv[1]},{"checkpoint_snapshot",argv[2]},{"problem_sha256",artifacts::detail::sha256(problem_bytes)},
    {"checkpoint_sha256",artifacts::detail::sha256(checkpoint_bytes)},{"verified_checkpoint_identity",key},
    {"saved_leg",saved.leg},{"saved_parameter",saved.parameter},{"saved_accepted_charts",saved.accepted_charts},
    {"start_leg",start_leg},{"start_parameter",start_parameter},{"source_rows",initial.coefficients.size()},
    {"physical_dimension",initial.columns()},{"extended_dimension",initial.coefficients.size()*initial.columns()+1},
    {"epsilon_low",saved.rows.low},{"epsilon_high",saved.rows.high},{"epsilon_width",saved.rows.high-saved.rows.low+1},
    {"working_bits",precision},{"order",options.taylor_order},{"centered_before_rational",options.centered_before_rational},
    {"centered_map_working_bits",options.centered_map_working_bits},{"centered_map_workers",options.centered_map_workers},
    {"circuit_grouped_dot",grouped_dot},{"compact_centered_only",centered_only},{"complete_path",complete},{"compact_centered_recovery",recovery},{"execution_identity",execution_key},{"seek_first_whole_input_crosscheck",seek},{"output_prefix",output_prefix},{"wall_budget_seconds",seconds},{"omitted_tails_included",false},{"full_guards_preserved",true}};
  std::cout<<json::serialize(report)<<'\n'<<std::flush;
  void* memory=mmap(nullptr,sizeof(Shared),PROT_READ|PROT_WRITE,MAP_SHARED|MAP_ANON,-1,0);if(memory==MAP_FAILED)throw std::runtime_error("shared profiling state allocation failed");
  auto* shared=new(memory) Shared;options.conditioning_stats=&shared->stats;
  if(!from_initial)options.continuation=std::make_shared<const AdjointContinuation>(saved);
  std::optional<AdjointContinuation> previous=from_initial?std::nullopt:std::optional<AdjointContinuation>(saved),current;
  AdjointConditioningStats previous_stats;std::string previous_key=key;
  const auto publish=[&](const std::string& suffix,const AdjointContinuation& state,const std::string& provenance_key){
    auto payload=adjoint_checkpoint::detail::encode(state,provenance_key);auto out=artifacts::detail::canonical(json::object{{"payload",payload},{"sha256",artifacts::detail::sha256(artifacts::detail::canonical(payload))}});
    std::ofstream file(output_prefix+suffix);file<<out<<'\n';if(!file)throw std::runtime_error("cannot publish profiler continuation");
  };
  options.continuation_observer=[&](const AdjointContinuation& state){current=state;};
  options.chart_observer=[&](unsigned leg,double from,double to,const LaurentRows& rows){
    ++shared->observed_charts;
    if(!current)throw std::logic_error("chart observer missing accepted continuation");
    bool final_chart=to==1;for(unsigned i=leg+1;i+1<path.size();++i)final_chart&=path[i]==path[i+1];
    const bool target=complete?final_chart:(!seek||shared->stats.rational_cross_checks>0);
    if(!output_prefix.empty()){
      if(target){if(previous)publish("-before.json",*previous,previous_key);publish("-after.json",*current,execution_key);
        std::ofstream meta(output_prefix+"-phase-counters.json");meta<<json::serialize(json::object{{"before",counters(previous_stats)},{"after",counters(shared->stats)},{"leg",leg},{"from",from},{"to",to}})<<'\n';}
      publish("-last.json",*current,execution_key);
    }
    std::cout<<json::serialize(json::object{{"status","chart_observed"},{"leg",leg},{"from",from},{"to",to},{"conditioning",counters(shared->stats)}})<<'\n'<<std::flush;
    if(!target){previous=*current;previous_stats=shared->stats;previous_key=execution_key;return;}
    shared->leg=leg;shared->from=from;shared->to=to;
    for(const auto& row:rows.coefficients)for(const auto& series:row)for(const auto& b:series){shared->finite&=b.is_finite();shared->output_radius=std::max({shared->output_radius,mag_get_d_log2_approx(arb_radref(acb_realref(b.raw()))),mag_get_d_log2_approx(arb_radref(acb_imagref(b.raw())))});}
    if(baseline_after){
      const auto& baseline=*baseline_after;shared->comparison=2;
      if(baseline.leg==leg&&baseline.parameter==to&&baseline.rows.low==rows.low&&baseline.rows.high==rows.high&&baseline.rows.coefficients.size()==rows.coefficients.size()&&baseline.rows.columns()==rows.columns()){
        shared->comparison=1;
        for(unsigned i=0;i<rows.coefficients.size();++i)for(unsigned j=0;j<rows.columns();++j)for(unsigned k=0;k<rows.coefficients[i][j].size();++k){
          const auto& a=baseline.rows.coefficients[i][j][k];const auto& b=rows.coefficients[i][j][k];shared->disjoint+=!acb_overlaps(a.raw(),b.raw());
          for(unsigned part=0;part<2;++part){auto x=part?acb_imagref(a.raw()):acb_realref(a.raw());auto y=part?acb_imagref(b.raw()):acb_realref(b.raw());
            arf_t delta;arf_init(delta);arf_sub(delta,arb_midref(x),arb_midref(y),B::precision(),ARF_RND_NEAR);shared->midpoint_distance=std::max(shared->midpoint_distance,std::abs(arf_get_d(delta,ARF_RND_NEAR)));arf_clear(delta);
            shared->baseline_width=std::max(shared->baseline_width,2*mag_get_d(arb_radref(x)));shared->result_width=std::max(shared->result_width,2*mag_get_d(arb_radref(y)));}
        }
      }
    }
    shared->status=complete?4:1;throw Accepted{};
  };
  const auto began=Clock::now();pid_t child=fork();if(child<0)throw std::runtime_error("fork failed");
  if(!child){try{transport_adjoint_rows(connection,initial,forcing,path,options);shared->status=2;}
    catch(const Accepted&){}catch(const std::exception& e){std::strncpy(shared->error,e.what(),sizeof(shared->error)-1);shared->status=3;}_exit(0);}
  int child_status=0;bool timeout=false;struct rusage usage{};
  for(;;){auto got=wait4(child,&child_status,WNOHANG,&usage);if(got==child)break;if(got<0&&errno!=EINTR)throw std::runtime_error("wait4 failed");
    if(std::chrono::duration<double>(Clock::now()-began).count()>=seconds){timeout=true;kill(child,SIGTERM);while(wait4(child,&child_status,0,&usage)<0&&errno==EINTR){}break;}usleep(10000);}
  report["elapsed_seconds"]=std::chrono::duration<double>(Clock::now()-began).count();report["conditioning"]=counters(shared->stats);report["observed_charts"]=shared->observed_charts;
  report["child_cpu_seconds"]=usage.ru_utime.tv_sec+usage.ru_utime.tv_usec/1e6+usage.ru_stime.tv_sec+usage.ru_stime.tv_usec/1e6;
  report["status"]=timeout?"wall_budget_exhausted":shared->status==4?"completed_path":shared->status==1?"accepted_one_chart":shared->status==2?"completed_without_chart":"error";
  report["phase_timers_complete"]=!timeout;report["timer_note"]=timeout?"In-flight RAII phase timer is not committed on timeout; completed phases and entered counters remain available.":"Completed chart timers include nested phases; do not sum them as disjoint times.";
  if(shared->status==1||shared->status==4){report["accepted_chart"]={{"leg",shared->leg},{"from",shared->from},{"to",shared->to},{"finite",bool(shared->finite)},{"maximum_component_radius_log2",std::isfinite(shared->output_radius)?json::value(shared->output_radius):json::value("-inf")}};}
  if(shared->comparison)report["baseline_comparison"]={{"same_position_and_shape",shared->comparison==1},{"all_overlap",shared->comparison==1&&shared->disjoint==0},{"disjoint_coefficients",shared->disjoint},{"maximum_component_midpoint_distance",shared->midpoint_distance},{"baseline_maximum_component_width",shared->baseline_width},{"result_maximum_component_width",shared->result_width}};
  if(shared->error[0])report["error"]=shared->error;
  if(WIFSIGNALED(child_status))report["child_signal"]=WTERMSIG(child_status);
  const int result_code=timeout?2:(WIFEXITED(child_status)&&shared->status!=3?0:1);
  std::cout<<json::serialize(report)<<'\n';shared->~Shared();munmap(memory,sizeof(Shared));return result_code;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
