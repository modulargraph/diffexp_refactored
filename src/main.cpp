#include "diffexp/canonical.hpp"
#include "diffexp/cli.hpp"
#include "diffexp/transport.hpp"
#include "diffexp/planar.hpp"
#include "diffexp/mpl.hpp"
#include "diffexp/rational_transport.hpp"
#include "diffexp/singular_transport.hpp"
#include "diffexp/original_banana_unequal.hpp"
#include "diffexp/bubble.hpp"
#include "diffexp/certified_bubble.hpp"
#include "diffexp/sunrise.hpp"
#include "diffexp/banana_oracle.hpp"
#include "diffexp/direct_examples.hpp"
#if defined(__unix__) || defined(__APPLE__)
#include "diffexp/recursion_graph.hpp"
#include "diffexp/family_config.hpp"
#include "diffexp/recursion_pipeline.hpp"
#include "diffexp/level_cache.hpp"
#include "diffexp/fire_modular.hpp"
#include "diffexp/ibp_solver_provider.hpp"
#include "diffexp/henn_boundary.hpp"
#include "diffexp/critical_points.hpp"
#include "diffexp/epsilon_sampling.hpp"
#include "diffexp/moving_poles.hpp"
#endif
#include <iostream>
#include <string>
#include <charconv>
#include <optional>
#include <chrono>
#include <cstdlib>
#include <thread>
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char** argv) {
#ifdef _OPENMP
  // Leave the cores that other work occupies at start-up: an oversubscribed
  // chart expansion waits a scheduler slice at every barrier. OMP_NUM_THREADS
  // overrides this.
  if(!std::getenv("OMP_NUM_THREADS")) {
    double load=0;
    if(getloadavg(&load,1)==1)omp_set_num_threads(std::max(1,omp_get_num_procs()-static_cast<int>(load+0.5)));
  }
#endif
  bool json_mode=false;
  try {
    const std::string command = argc > 1 ? argv[1] : "--help";
    if (command == "--version") {
      std::cout << "DiffExp 2.2.1 (standalone C++20 / FLINT)\n";
    } else if(command == "transport") {
      json_mode=true;if(argc!=3)throw std::invalid_argument("transport requires a JSON file or -");
      return diffexp::transport::run_file(argv[2]);
    } else if(command == "family-template") {
#if defined(__unix__) || defined(__APPLE__)
      if(argc!=3)throw std::invalid_argument("family-template requires one example name");
      std::cout<<boost::json::serialize(diffexp::family_config::describe(diffexp::feynman::example_family(argv[2])))<<'\n';
#else
      throw std::runtime_error("Feynman preparation requires a POSIX platform");
#endif
    } else if(command == "prepare" || command == "ft") {
#if defined(__unix__) || defined(__APPLE__)
      for(int i=3;i<argc;++i)if(std::string(argv[i])=="--json")json_mode=true;
      auto& progress=json_mode?std::cerr:std::cout;
      if(argc<3)throw std::invalid_argument(command+" requires a family configuration JSON file");
      auto configuration=diffexp::family_config::load(argv[2]);
      bool epsilon_order_given=configuration.epsilon_order_given;
      auto options=configuration.preparation;
      auto numerical=configuration.numerical;
      // The command line accumulates top-down wherever that is covered; the
      // library default is one row per requested integral.
      numerical.accumulation=diffexp::recursion::Accumulation::weighted;
      const auto preparation_start=std::chrono::steady_clock::now();
      diffexp::AdjointConditioningStats conditioning;
      numerical.adjoint.conditioning_stats=&conditioning;
      unsigned epsilon_order=configuration.epsilon_order;
      bool numerical_cache=true;
      bool fire_given=false,fire_prime_given=false;
      auto cache=std::filesystem::temp_directory_path()/"diffexp-exact-cache";
      std::string henn_basis,epsilon_shift,moving_poles;
      auto epsilon_strategy=configuration.epsilon_strategy,sample_step=configuration.sample_step;
      unsigned ibp_max_degree=0;std::string merge_selection="auto",boundary="auto";
      unsigned samples=configuration.samples,sample_jobs=0;int sample_poles=-1;bool leaf_digits_given=false;
      // What a sample process is started with: this command line without the
      // options that select and shape the sampling itself.
      std::vector<std::string> sample_arguments;
      std::filesystem::path fire_prime;
      auto ibp_provider=configuration.ibp_provider;std::string ibp_validation="anchored";auto ibp_dots=configuration.ibp_dots,ibp_numerators=configuration.ibp_numerators;
      std::string ibp_seeds="demand";
      // As for the chart expansion, the built-in solver leaves the cores that
      // other work occupies at start-up; --ibp-threads overrides this.
      unsigned ibp_threads=std::max(1u,std::thread::hardware_concurrency());
      {double load=0;if(getloadavg(&load,1)==1)ibp_threads=static_cast<unsigned>(std::max(1,static_cast<int>(ibp_threads)-static_cast<int>(load+0.5)));}
      for(int i=3;i<argc;++i) {
        const std::string option=argv[i];
        if(option=="--json")continue;
        if(command=="ft" && option=="--no-numerical-cache"){numerical_cache=false;sample_arguments.push_back(option);continue;}
        if(option=="--leaf-digits")leaf_digits_given=true;
        if(option!="--epsilon-strategy" && option!="--samples" && option!="--sample-step" && option!="--sample-jobs" &&
            option!="--sample-poles" && option!="--epsilon-order" && option!="--epsilon-shift" && i+1<argc) {
          sample_arguments.push_back(option);sample_arguments.push_back(argv[i+1]);
        }
        if(option=="--ft-transport" || option=="--ft-transport-digits")
          throw std::invalid_argument(option+" was removed: Taylor transport is the only ordinary method");
        if(i+1>=argc)throw std::invalid_argument("missing "+command+" option value");
        const char* option_value=argv[++i];
        if(option=="--fire"){options.reduction.provider.executable=option_value;fire_given=true;}
        else if(option=="--ibp-provider")ibp_provider=option_value;
        else if(option=="--ibp-validation"){ibp_validation=option_value;if(ibp_validation!="anchored"&&ibp_validation!="strict")throw std::invalid_argument("IBP validation must be anchored or strict");}
        else if(option=="--ibp-dots"||option=="--ibp-numerators"){
          const std::string value=option_value;unsigned n=0;auto parsed=std::from_chars(value.data(),value.data()+value.size(),n);
          if(parsed.ec!=std::errc{}||parsed.ptr!=value.data()+value.size()||n>8)throw std::invalid_argument("IBP seed powers must be 0..8");
          (option=="--ibp-dots"?ibp_dots:ibp_numerators)=n;
        }
        else if(option=="--ibp-seeds"){ibp_seeds=option_value;if(ibp_seeds!="demand"&&ibp_seeds!="sector"&&ibp_seeds!="subsectors"&&ibp_seeds!="all")throw std::invalid_argument("IBP seeds must be demand, sector, subsectors or all");}
        else if(option=="--ibp-max-degree"){
          const std::string value=option_value;auto parsed=std::from_chars(value.data(),value.data()+value.size(),ibp_max_degree);
          if(parsed.ec!=std::errc{}||parsed.ptr!=value.data()+value.size()||ibp_max_degree<4||ibp_max_degree>256)throw std::invalid_argument("IBP maximal degree must be 4..256");
        }
        else if(option=="--merge-selection"){
          merge_selection=option_value;
          if(merge_selection!="auto"&&merge_selection!="first")throw std::invalid_argument("--merge-selection must be auto or first");
        }
        else if(option=="--boundary"){
          boundary=option_value;
          if(boundary!="anchor"&&boundary!="endpoint"&&boundary!="auto")throw std::invalid_argument("--boundary must be anchor, endpoint or auto");
        }
        else if(option=="--ibp-threads"){
          const std::string value=option_value;auto parsed=std::from_chars(value.data(),value.data()+value.size(),ibp_threads);
          if(parsed.ec!=std::errc{}||parsed.ptr!=value.data()+value.size()||!ibp_threads||ibp_threads>256)throw std::invalid_argument("IBP threads must be 1..256");
        }
        else if(option=="--fire-prime"){fire_prime=option_value;fire_prime_given=true;}
        else if(option=="--cache")cache=option_value;
        else if(option=="--henn-basis")henn_basis=option_value;
        else if(command=="ft" && option=="--epsilon-shift")epsilon_shift=option_value;
        else if(command=="ft" && option=="--epsilon-strategy") {
          epsilon_strategy=option_value;
          if(epsilon_strategy!="expand" && epsilon_strategy!="sample")throw std::invalid_argument("--epsilon-strategy must be expand or sample");
        }
        else if(command=="ft" && option=="--sample-step") {
          sample_step=option_value;
          if(diffexp::Rational(sample_step).sign()<=0)throw std::invalid_argument("--sample-step must be a positive rational");
        }
        else if(command=="ft" && (option=="--samples" || option=="--sample-jobs" || option=="--sample-poles")) {
          const std::string value=option_value;unsigned number=0;
          const auto parsed=std::from_chars(value.data(),value.data()+value.size(),number);
          if(parsed.ec!=std::errc{} || parsed.ptr!=value.data()+value.size())throw std::invalid_argument(option+" must be a nonnegative integer");
          if(option=="--samples"){if(number<4 || number>200 || number%2)throw std::invalid_argument("--samples must be even, 4..200");samples=number;}
          else if(option=="--sample-jobs"){if(!number || number>1024)throw std::invalid_argument("--sample-jobs must be in 1..1024");sample_jobs=number;}
          else {if(number>64)throw std::invalid_argument("--sample-poles must be in 0..64");sample_poles=static_cast<int>(number);}
        }
        else if(option=="--moving-poles") {
          moving_poles=option_value;
          if(moving_poles!="keep" && moving_poles!="remove" && moving_poles!="exchange")throw std::invalid_argument("--moving-poles must be keep, remove or exchange");
        }
        else if(option=="--fire-seconds" || option=="--level-seconds" || option=="--total-seconds" ||
            option=="--fire-threads" || option=="--fire-simplifier-threads" || option=="--fire-memory-mib") {
          const std::string value=option_value;unsigned number=0;
          const auto parsed=std::from_chars(value.data(),value.data()+value.size(),number);
          if(parsed.ec!=std::errc{} || parsed.ptr!=value.data()+value.size() || !number || number>86400)
            throw std::invalid_argument("preparation resource budgets must be integers in 1..86400");
          if(option=="--fire-seconds")options.reduction.provider.timeout_seconds=number;
          else if(option=="--level-seconds")options.reduction.total_timeout_seconds=number;
          else if(option=="--total-seconds")options.total_timeout_seconds=number;
          else if(option=="--fire-threads")options.reduction.provider.threads=number;
          else if(option=="--fire-simplifier-threads")options.reduction.provider.simplifier_threads=number;
          else options.reduction.provider.memory_bytes=static_cast<std::size_t>(number)*1024*1024;
        }
        else if(command=="ft" && option=="--uncertainty") {
          const std::string mode=option_value;
          if(mode!="probe" && mode!="map")throw std::invalid_argument("FT uncertainty propagation must be probe or map");
          numerical.adjoint.uncertainty=mode=="map"?diffexp::AdjointOptions::Uncertainty::map:diffexp::AdjointOptions::Uncertainty::probe;
        }
        else if(command=="ft" && option=="--accumulation") {
          const std::string mode=option_value;
          if(mode!="rows" && mode!="weighted")throw std::invalid_argument("FT accumulation must be rows or weighted");
          numerical.accumulation=mode=="rows"?diffexp::recursion::Accumulation::rows:diffexp::recursion::Accumulation::weighted;
        }
        else if(command=="ft" && option=="--transport-bits") {
          char* end=nullptr;const auto value=std::strtoul(option_value,&end,10);
          if(!end || *end || value<16 || value>100000)throw std::invalid_argument("transport bits must be in 16..100000");
          numerical.adjoint.truncation_bits=static_cast<unsigned>(value);
        }
        else if(command=="ft" && option=="--method") {
          const std::string method=option_value;
          if(method=="values")throw std::invalid_argument("--method values was removed: use adjoint (default), factored or auto");
          if(method!="adjoint" && method!="factored" && method!="auto")
            throw std::invalid_argument("FT method must be adjoint, factored or auto");
          numerical.linear_method=method=="factored"?diffexp::recursion::LinearMethod::factored:
            method=="auto"?diffexp::recursion::LinearMethod::automatic:diffexp::recursion::LinearMethod::adjoint;
        }
        else if(command=="ft") {
          const std::string value=option_value;unsigned number=0;
          const auto parsed=std::from_chars(value.data(),value.data()+value.size(),number);
          if(parsed.ec!=std::errc{} || parsed.ptr!=value.data()+value.size())throw std::invalid_argument("FT numerical option must be a nonnegative integer");
          if(option=="--epsilon-order") {
            if(number>100)throw std::invalid_argument("epsilon order must be in 0..100");
            epsilon_order=number;epsilon_order_given=true;
          }
          // Same ranges as the "numerical" block of a configuration file.
          else if(option=="--endpoint-order"){if(number<1 || number>1000)throw std::invalid_argument("endpoint order must be in 1..1000");numerical.endpoint_order=number;}
          else if(option=="--overlap-ratio"){if(number<4 || number>65536)throw std::invalid_argument("overlap ratio must be in 4..65536");numerical.overlap_ratio=number;}
          else if(option=="--ordinary-order"){if(number<8 || number>1000)throw std::invalid_argument("ordinary order must be in 8..1000");numerical.ordinary_order=number;}
          else if(option=="--precision-bits"){if(number<64 || number>1000000)throw std::invalid_argument("precision bits must be in 64..1000000");numerical.working_bits=number;}
          else if(option=="--leaf-digits"){if(number<1 || number>100000)throw std::invalid_argument("leaf digits must be in 1..100000");numerical.leaf_digits=number;}
          else if(option=="--basis-scan-seconds"){if(number>60)throw std::invalid_argument("basis scan seconds must be in 0..60");numerical.automatic_basis_scan=number!=0;if(number)numerical.basis_scan_seconds=number;}
          else throw std::invalid_argument("unknown FT option: "+option);
        } else throw std::invalid_argument("unknown prepare option: "+option);
      }
      if(ibp_provider=="auto")ibp_provider=!fire_prime.empty()?"fire-modular":!options.reduction.provider.executable.empty()?"fire":"ibp-solver";
      if(ibp_provider!="ibp-solver"&&ibp_provider!="fire"&&ibp_provider!="fire-modular")throw std::invalid_argument("IBP provider must be auto, ibp-solver, fire or fire-modular");
      if(ibp_provider=="fire-modular"&&fire_prime.empty())throw std::invalid_argument("fire-modular requires --fire-prime PATH");
      if(ibp_provider=="fire"&&options.reduction.provider.executable.empty())throw std::invalid_argument("fire provider requires --fire PATH");
      // An executable that the selected provider never runs is almost always a
      // mistake about which provider the configuration pins.
      if(fire_given && ibp_provider!="fire")
        std::cerr<<"diffexp: warning: --fire is ignored because the IBP provider is "<<ibp_provider
          <<(configuration.ibp_provider==ibp_provider?" (pinned by the configuration)":"")<<"; pass --ibp-provider fire to use it\n";
      if(fire_prime_given && ibp_provider!="fire-modular")
        std::cerr<<"diffexp: warning: --fire-prime is ignored because the IBP provider is "<<ibp_provider<<"; pass --ibp-provider fire-modular to use it\n";
      options.reduction.batch_cache_directory=cache/"fire-batches";
      options.reduction.pending_cache_directory=cache/"fire-pending";
      if(numerical_cache){numerical.basis_cache_directory=cache/"physical-bases";numerical.ordinary_cache_directory=cache/"ordinary-transport";}
      auto requested=configuration.integrals;
      std::optional<diffexp::henn::CanonicalBasis> canonical_basis;
      if(!henn_basis.empty()) {
        if(!diffexp::family_config::same_geometry(configuration.family,diffexp::feynman::example_family("henn_double_pentagon_x0")))
          throw std::invalid_argument("--henn-basis requires the published Henn X0 geometry and denominator order; a family label alone is insufficient");
        canonical_basis=diffexp::henn::read_x0(henn_basis);
        requested=canonical_basis->scalar_targets;
        if(!epsilon_order_given)epsilon_order=4;
      }
      if(command=="ft") {
        // Everything the numerical phase would reject is rejected here, before
        // the exact IBP preparation, which can take minutes to hours.
        diffexp::recursion::validate_numerical_options(numerical);
        if(canonical_basis && (numerical.working_bits<64 || numerical.working_bits>4096))
          throw std::invalid_argument("--henn-basis projects at 64..4096 working bits; --precision-bits "+std::to_string(numerical.working_bits)+" is outside that range");
        if(canonical_basis) {
          const auto needed=diffexp::henn::needed_scalar_high(*canonical_basis,epsilon_order);
          if(needed>static_cast<int>(numerical.max_epsilon))
            throw std::invalid_argument("--henn-basis at epsilon order "+std::to_string(epsilon_order)+" needs scalar order "+std::to_string(needed)+", beyond the epsilon budget "+std::to_string(numerical.max_epsilon));
        }
      }
      const auto target_seeds=ibp_seeds=="all"?diffexp::ibp_solver::Options::TargetSeeds::all:ibp_seeds=="subsectors"?diffexp::ibp_solver::Options::TargetSeeds::subsectors:ibp_seeds=="sector"?diffexp::ibp_solver::Options::TargetSeeds::sector:diffexp::ibp_solver::Options::TargetSeeds::none;
      diffexp::artifacts::Store store(cache);unsigned reused=0,built=0;
      diffexp::ibp_solver::Statistics ibp_statistics;
      auto provider=[&](const auto& basis,const auto& dimension,const auto& field,auto parameter,const auto& sources,const auto& budget) {
        diffexp::level::Provider fallback;
        std::shared_ptr<diffexp::ibp_solver::Session> native_session;
        if(ibp_provider=="ibp-solver"){
          // --fire-seconds bounds one finite-field elimination of the built-in
          // solver as it bounds one FIRE call.
          auto native=diffexp::ibp_solver::recursion_options(cache/"ibp-solver",parameter,
            options.reduction.provider.timeout_seconds,ibp_validation=="strict");
          native.dots=ibp_dots;native.numerators=ibp_numerators;native.target_seeds=target_seeds;native.threads=ibp_threads;
          // A larger total degree comes with larger degrees per variable and
          // the samples to see them.
          if(ibp_max_degree) {
            native.max_degree=ibp_max_degree;
            native.max_variable_degree=std::max(native.max_variable_degree,4*ibp_max_degree);
            native.max_samples_per_prime=std::max<std::size_t>(native.max_samples_per_prime,64u*ibp_max_degree);
          }
          native.progress=[](const std::string& message){std::cerr<<"IBP solver: "<<message<<'\n';};
          // A reduction that meets a degree or sample bound is repeated with
          // the bounds doubled, as often as it takes: the time and memory
          // budgets are what ends a reduction, not a degree chosen in advance.
          fallback=[&,current=std::make_shared<diffexp::ibp_solver::Options>(native)](const auto& batch,const auto& limits){
            for(;;) {
              if(!native_session)native_session=std::make_shared<diffexp::ibp_solver::Session>(basis,dimension,field,*current);
              auto result=(*native_session)(batch,limits);
              if(result.success)return result;
              auto wider=diffexp::ibp_solver::escalate(*current,result.reason);
              if(!wider)return result;
              std::cerr<<"IBP solver: "<<result.reason<<"; repeating with total degree "<<wider->max_degree<<", degree per variable "
                <<wider->max_variable_degree<<" and "<<wider->max_samples_per_prime<<" samples per prime\n";
              *current=*wider;native_session.reset();
            }
          };
        }else if(ibp_provider=="fire-modular") {
          diffexp::fire_modular::Options modular;modular.executable=fire_prime;modular.cache_directory=cache/"fire-modular";
          modular.progress=[](const std::string& message){std::cerr<<message<<'\n';};
          auto session=std::make_shared<diffexp::fire_modular::Session>(basis,dimension,field,modular,budget.batch_cache_directory);
          fallback=[session](const auto& batch,const auto& limits){return (*session)(batch,limits);};
        }
        auto result=diffexp::cached_level::prepare(store,basis,dimension,field,parameter,sources,budget,fallback);
        if(native_session){const auto& stats=native_session->statistics();ibp_statistics.full_solve_seconds+=stats.full_solve_seconds;ibp_statistics.trace_learning_seconds+=stats.trace_learning_seconds;ibp_statistics.trace_replay_seconds+=stats.trace_replay_seconds;ibp_statistics.templates+=stats.templates;ibp_statistics.full_solves+=stats.full_solves;ibp_statistics.trace_replays+=stats.trace_replays;ibp_statistics.trace_fallbacks+=stats.trace_fallbacks;ibp_statistics.probes+=stats.probes;ibp_statistics.equations+=stats.equations;ibp_statistics.generation_seconds+=stats.generation_seconds;ibp_statistics.elimination_seconds+=stats.elimination_seconds;
          ibp_statistics.selected_solves+=stats.selected_solves;ibp_statistics.selected_solve_seconds+=stats.selected_solve_seconds;ibp_statistics.selection_fallbacks+=stats.selection_fallbacks;ibp_statistics.trace_disabled+=stats.trace_disabled;
          ibp_statistics.batched_replays+=stats.batched_replays;ibp_statistics.early_solves+=stats.early_solves;ibp_statistics.early_solves_used+=stats.early_solves_used;ibp_statistics.seed_widenings+=stats.seed_widenings;
          ibp_statistics.continued_solves+=stats.continued_solves;ibp_statistics.continued_equations+=stats.continued_equations;ibp_statistics.reused_solves+=stats.reused_solves;ibp_statistics.continued_solve_seconds+=stats.continued_solve_seconds;}
        if(result.cache_hit)++reused;else if(result.result.success)++built;
        return std::move(result.result);
      };
      // Masters without moving poles (docs/moving-poles.md). An expanded run
      // removes them by default; a sample gains nothing from it.
      // A sample exchanges the masters themselves: its endpoint series are then
      // free of the moving poles too, which otherwise force a tiny overlap and
      // hundreds of charts next to an endpoint (docs/performance.md).
      if(moving_poles.empty())moving_poles=command!="ft"?"keep":epsilon_shift.empty() && epsilon_strategy!="sample"?"remove":"exchange";
      boost::json::array moving_pole_levels;double moving_pole_seconds=0;
      const auto describe_moving_poles=[&](std::size_t depth,const diffexp::moving_poles::Result& result) {
        moving_pole_seconds+=result.seconds;
        if(!result.before.factors)return;
        progress<<"Level "<<depth+1<<": moving poles: "<<result.before.factors<<" mixed factors of total order "<<result.before.total_order
          <<" -> "<<result.after.total_order<<"; connection characters "<<result.original_characters<<" -> "<<result.characters<<"; "<<result.seconds<<" s\n";
        boost::json::array steps;
        for(const auto& step:result.steps) {
          progress<<"  ("<<step.factor<<"): master "<<step.replaced<<" replaced by "<<step.source<<"\n";
          steps.push_back(boost::json::object{{"factor",step.factor},{"replaced_master",step.replaced},{"by",step.source}});
        }
        boost::json::array remaining;
        for(const auto& factor:result.remaining){progress<<"  left in place: ("<<factor<<")\n";remaining.emplace_back(factor);}
        if(result.budget_reached)progress<<"  stopped by its step or time budget\n";
        moving_pole_levels.push_back(boost::json::object{{"level",depth+1},{"mixed_factors",result.before.factors},{"mixed_order_before",result.before.total_order},
          {"mixed_order_after",result.after.total_order},{"characters_before",result.original_characters},{"characters_after",result.characters},
          {"steps",std::move(steps)},{"left_in_place",std::move(remaining)},{"budget_reached",result.budget_reached},{"seconds",result.seconds}});
        progress<<std::flush;
      };
      // Each level merges the pair of propagators with the fewest masters.
      // Explicit anchors and a causal prescription are stated level by level
      // for the ladder as written, so they keep it.
      if(merge_selection=="auto" && ibp_provider=="ibp-solver" && options.merges.empty() && options.anchors.empty() &&
          (!numerical.causal_prescription || (numerical.causal_prescription->derive_orientations && numerical.causal_prescription->points.empty()))) {
        options.merge_chosen=[&progress](unsigned level,unsigned left,unsigned right,std::size_t candidates,std::size_t masters,std::size_t first) {
          if(left==0 && right==1)return;
          progress<<"Level "<<level<<": merging propagators "<<left<<" and "<<right<<", the pair with the fewest masters of "
            <<candidates<<" ("<<masters<<"; the first pair has "<<first<<")\n"<<std::flush;
        };
        // From the graph where Singular is installed: every pair, in seconds.
        if(diffexp::critical_points::available()) {
          auto physical=configuration.family.momenta;physical.lines.resize(configuration.family.physical_count);
          options.count_critical_points=[physical,directory=cache/"critical-points"](const auto& lines,unsigned left,unsigned right) {
            return diffexp::critical_points::count(physical,lines,left,right,directory);
          };
        } else progress<<"Singular is not installed: merges are chosen by counting masters with the IBP solver, for large levels and pairs on one loop momentum only.\n";
        options.count_masters=[&](const auto& basis,const auto& dimension,const auto& field,auto parameter,const auto& sources,const auto& budget)
            ->std::optional<std::size_t> {
          auto native=diffexp::ibp_solver::recursion_options(cache/"ibp-solver",parameter,options.reduction.provider.timeout_seconds,false);
          native.dots=ibp_dots;native.numerators=ibp_numerators;native.target_seeds=target_seeds;native.threads=1;
          try{return diffexp::ibp_solver::Session(basis,dimension,field,native).count_masters(sources,budget.provider);}
          catch(const std::exception&){return std::nullopt;}
        };
      }
      // A level whose endpoint qualifies takes its boundary there
      // (endpoint_boundary.hpp). Explicit anchors, a stated causal
      // prescription and the Henn projection are written for the ladder of
      // interior anchors and keep it.
      const bool stated_ladder=!options.anchors.empty() || !henn_basis.empty() ||
        (numerical.causal_prescription && !(numerical.causal_prescription->derive_orientations && numerical.causal_prescription->points.empty()));
      if(boundary=="endpoint" && stated_ladder)
        throw std::invalid_argument("--boundary endpoint cannot be combined with explicit anchors, stated causal sides or --henn-basis");
      options.boundary=boundary=="anchor" || stated_ladder?diffexp::recursion::BoundaryPolicy::anchor:
        boundary=="endpoint"?diffexp::recursion::BoundaryPolicy::endpoint:diffexp::recursion::BoundaryPolicy::automatic;
      options.endpoint_refused=[&progress](unsigned level,bool upper,const std::string& reason) {
        progress<<"Level "<<level<<": boundary not taken at x = "<<(upper?1:0)<<": "<<reason<<"\n"<<std::flush;
      };
      // exchange: a level swaps masters for cached integrals of the same sector
      // before its masters become the child's requests.
      if(moving_poles=="exchange")options.exchange_masters=[&](diffexp::level::Result level,unsigned depth) {
        diffexp::moving_poles::Result result;
        level=diffexp::moving_poles::exchange(std::move(level),0,1,{},&result);
        describe_moving_poles(depth,result);return level;
      };
      // A spent time budget names the three knobs that bound the preparation.
      auto graph=[&] {
        try {
          return diffexp::recursion::prepare(configuration.family,std::move(requested),options,provider,
            [&](const auto& event) {
              progress<<"Level "<<event.depth<<": "<<event.physical_propagators<<" propagators; ";
              if(event.scalar_leaf)progress<<"analytic scalar boundary";
              else progress<<event.masters<<" masters; "<<event.sources<<" shared source integrals";
              progress<<"; elapsed "<<event.elapsed_seconds<<" s\n";
              if(event.boundary!=diffexp::recursion::BoundarySite::anchor)progress<<"Level "<<event.depth<<": boundary at the endpoint x = "<<event.anchor.str()
                <<": the next level is the family with propagator "<<(event.boundary==diffexp::recursion::BoundarySite::lower?event.left:event.right)<<" removed\n";
              else if(event.anchor!=event.requested_anchor)progress<<"Level "<<event.depth<<": requested anchor "<<event.requested_anchor.str()
                <<" is singular in the selected master basis; matching at "<<event.anchor.str()<<" instead\n";
              progress<<std::flush;
            });
        } catch(const std::runtime_error& error) {
          const std::string message=error.what();
          if(message.find("time budget")==std::string::npos && message.find("exhausted reconstruction budget")==std::string::npos)throw;
          throw std::runtime_error(message+". Preparation budgets in force: --fire-seconds "+std::to_string(options.reduction.provider.timeout_seconds)
            +" (one IBP elimination), --level-seconds "+std::to_string(options.reduction.total_timeout_seconds)
            +", --total-seconds "+std::to_string(options.total_timeout_seconds)
            +"; they can also be set as preparation.fire_seconds, level_seconds and total_seconds in the family file");
        }
      }();
      progress<<"Exact recursion prepared: "<<graph.nodes.size()<<" levels, "<<built<<" systems built, "<<reused<<" systems reused.\n"
        <<"Verified cache: "<<std::filesystem::absolute(cache)<<"\n";
      boost::json::object report{{"schema","DiffExp.FeynmanTrick/v1"},{"family",configuration.family.name},
        {"status",command=="prepare"?"prepared":"completed"},{"levels",graph.nodes.size()},
        {"systems_built",built},{"systems_reused",reused}};
      report["ibp_provider"]=ibp_provider;
      {
        // Per level: its masters and where it takes its boundary.
        boost::json::array ladder;
        for(const auto& event:graph.events) {
          if(event.scalar_leaf)continue;
          ladder.push_back(boost::json::object{{"level",event.depth},{"masters",event.masters},
            {"boundary",event.boundary==diffexp::recursion::BoundarySite::lower?"lower":event.boundary==diffexp::recursion::BoundarySite::upper?"upper":"anchor"},
            {"x",event.anchor.str()}});
        }
        report["ladder"]=std::move(ladder);
      }
      // Verified level closures are exact and provider independent, and their
      // cache key deliberately omits the provider. A reused closure therefore
      // keeps the master basis of whichever provider and seeds first built it.
      report["ibp_cache"]=boost::json::object{{"systems_reused",reused},{"systems_built",built},{"provider_options_ignored",reused>0},
        {"note",reused?"reused level closures were not recomputed: --ibp-provider, --ibp-dots and --ibp-numerators did not apply to them; use a fresh --cache to apply them":"all level closures were built with the reported IBP provider and seeds"}};
      if(reused)std::cerr<<"diffexp: note: "<<reused<<" level system"<<(reused==1?" was":"s were")<<" reused from the verified cache; the IBP provider options (--ibp-provider "
        <<ibp_provider<<", --ibp-dots, --ibp-numerators) did not apply to "<<(reused==1?"it":"them")<<". Use a fresh --cache to rebuild.\n";
      report["ibp_statistics"]=boost::json::object{{"full_solve_seconds",ibp_statistics.full_solve_seconds},{"trace_learning_seconds",ibp_statistics.trace_learning_seconds},{"trace_replay_seconds",ibp_statistics.trace_replay_seconds},{"compiled_templates",ibp_statistics.templates},{"full_solves",ibp_statistics.full_solves},{"trace_replays",ibp_statistics.trace_replays},{"trace_fallbacks",ibp_statistics.trace_fallbacks},{"fresh_probes",ibp_statistics.probes},{"generated_equations",ibp_statistics.equations},{"probe_generation_seconds",ibp_statistics.generation_seconds},{"probe_elimination_seconds",ibp_statistics.elimination_seconds},
        {"validation",ibp_validation},{"subset_solves",ibp_statistics.selected_solves},{"subset_solve_seconds",ibp_statistics.selected_solve_seconds},{"subset_fallbacks",ibp_statistics.selection_fallbacks},{"traces_disabled",ibp_statistics.trace_disabled},
        {"threads",ibp_threads},{"batched_replays",ibp_statistics.batched_replays},{"early_solves",ibp_statistics.early_solves},{"early_solves_used",ibp_statistics.early_solves_used},{"seed_widenings",ibp_statistics.seed_widenings},
        {"continued_solves",ibp_statistics.continued_solves},{"continued_equations",ibp_statistics.continued_equations},{"reused_solves",ibp_statistics.reused_solves},{"continued_solve_seconds",ibp_statistics.continued_solve_seconds},{"seeds",ibp_seeds}};
      report["timings"]=boost::json::object{{"preparation_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-preparation_start).count()}};
      const auto report_moving_poles=[&] {
        if(moving_poles!="keep")report["moving_poles"]=boost::json::object{{"mode",moving_poles},{"levels",moving_pole_levels},{"seconds",moving_pole_seconds}};
      };
      report_moving_poles();
      if(command=="ft" && epsilon_strategy=="sample") {
        namespace es=diffexp::epsilon_sampling;
        if(!epsilon_shift.empty())throw std::invalid_argument("--epsilon-strategy sample chooses its own shifts; do not combine it with --epsilon-shift");
        if(canonical_basis)throw std::invalid_argument("--epsilon-strategy sample cannot be combined with --henn-basis");
        const unsigned loops=configuration.family.momenta.loops,pole_bound=2*loops;
        if(!samples){samples=pole_bound+epsilon_order+5;samples+=samples%2;samples=std::max(samples,8u);}
        // The leaf is re-expanded from its Laurent series at zero, and a deep
        // recursion multiplies what that loses: ask it for fewer digits.
        if(!leaf_digits_given) {
          sample_arguments.push_back("--leaf-digits");
          sample_arguments.push_back(std::to_string(std::min<unsigned>(numerical.leaf_digits,graph.nodes.size()>=6?12:16)));
        }
        unsigned cores=std::max(1u,std::thread::hardware_concurrency());
        {double load=0;if(getloadavg(&load,1)==1)cores=static_cast<unsigned>(std::max(1,static_cast<int>(cores)-static_cast<int>(load+0.5)));}
        if(!sample_jobs)sample_jobs=std::min(cores,samples);
        std::error_code ignored;
        auto self=std::filesystem::read_symlink("/proc/self/exe",ignored);
        if(ignored || self.empty())self=argv[0];
        std::vector<std::string> child{self.string(),"ft",argv[2]};
        child.insert(child.end(),sample_arguments.begin(),sample_arguments.end());
        if(std::find(sample_arguments.begin(),sample_arguments.end(),"--cache")==sample_arguments.end()){child.push_back("--cache");child.push_back(cache.string());}
        // A family file that selects sampling itself would make each sample
        // sample again: the samples are told to expand.
        child.insert(child.end(),{"--json","--epsilon-order","0","--epsilon-strategy","expand"});
        const auto points=es::nodes(diffexp::Rational(sample_step),samples);
        progress<<"Epsilon sampling: "<<samples<<" samples at epsilon = +-k*"<<sample_step<<", "<<sample_jobs<<" at a time.\n"<<std::flush;
        const auto numerical_start=std::chrono::steady_clock::now();
        const auto done=es::run(child,points,cache/"epsilon-samples",sample_jobs,std::max(1u,cores/sample_jobs),
          [&](const es::Sample& sample,std::size_t finished,std::size_t total) {
            progress<<"  sample "<<finished<<"/"<<total<<" at epsilon = "<<sample.epsilon.str()<<": "<<sample.seconds<<" s; elapsed "
              <<std::chrono::duration<double>(std::chrono::steady_clock::now()-numerical_start).count()<<" s\n"<<std::flush;
          });
        // The cancellations of the fit are exact rational arithmetic on the
        // nodes; carry the sample values well beyond their own accuracy.
        diffexp::kernel::ComplexBall::set_precision(std::max<slong>(numerical.working_bits,1024)+64*samples);
        std::vector<std::vector<diffexp::kernel::ComplexBall>> values(done.front().values.size());
        double cpu=0;boost::json::array sample_report;
        for(const auto& sample:done) {
          cpu+=sample.seconds;sample_report.push_back(boost::json::object{{"epsilon",sample.epsilon.str()},{"numerical_seconds",sample.seconds}});
          for(std::size_t c=0;c<values.size();++c)values[c].push_back(sample.values.at(c));
        }
        const unsigned poles=sample_poles>=0?static_cast<unsigned>(sample_poles):es::pole_order(points,values,pole_bound);
        const auto fitted=es::fit(points,values,poles,epsilon_order);
        report["timings"].as_object()["numerical_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-numerical_start).count();
        report["timings"].as_object()["sample_cpu_seconds"]=cpu;
        report["epsilon_strategy"]="sample";
        report["epsilon_sampling"]=boost::json::object{{"samples",samples},{"step",sample_step},{"jobs",sample_jobs},{"pole_order",poles},
          {"pole_order_detected",sample_poles<0},{"points",std::move(sample_report)}};
        boost::json::array rows,changes;
        for(std::size_t c=0;c<fitted.coefficients.size();++c) {
          boost::json::array powers,change;
          for(std::size_t k=0;k<fitted.coefficients[c].size();++k) {
            powers.push_back(diffexp::ball_json(fitted.coefficients[c][k],45));
            change.push_back(std::isfinite(fitted.log2_change[c][k])?std::pow(2.0,fitted.log2_change[c][k]):0.0);
            if(!json_mode){progress<<"Scalar "<<c+1<<" epsilon^"<<static_cast<long>(k)-static_cast<long>(poles)<<" = ";acb_printn(fitted.coefficients[c][k].raw(),32,0);progress<<'\n';}
          }
          rows.push_back(std::move(powers));changes.push_back(std::move(change));
        }
        report["epsilon_low"]=-static_cast<long>(poles);report["epsilon_high"]=epsilon_order;
        report["coefficient_order"]="component,epsilon";
        report["coefficients"]=std::move(rows);
        report["fit_change_without_outermost_samples"]=std::move(changes);
        report["component_basis"]="requested_scalar";
        report["ball_encoding"]="Arb decimal interval strings for real and imaginary parts; radii include the fit change, an estimate";
        report["omitted_tails_certified"]=false;
        report["working_bits"]=numerical.working_bits;
        progress<<"Laurent coefficients fitted from "<<samples<<" samples (pole order "<<poles<<(sample_poles<0?", detected":"")
          <<"); radii are estimates from dropping the two outermost samples.\n";
      }
      else if(command=="prepare")progress<<"This prepares equations and recursive operations; it does not report a numerical integral.\n";
      else {
        numerical.progress=[&progress,start=std::chrono::steady_clock::now()](std::size_t depth,const std::string& phase,int high) {
          const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
          progress<<"Numerical level "<<depth+1<<": "<<phase<<"; epsilon through "<<high<<"; elapsed "<<elapsed<<" s\n"<<std::flush;
        };
        // A transport that runs for minutes reports its position once a minute.
        numerical.ordinary_fraction_progress=[&progress,start=std::chrono::steady_clock::now(),last=std::chrono::steady_clock::now()](double fraction) mutable {
          const auto now=std::chrono::steady_clock::now();
          if(fraction<=0){last=now;return;}
          if(now-last<std::chrono::seconds(60))return;
          last=now;
          progress<<"  adjoint transport at "<<static_cast<int>(100*fraction)<<"% of its path; elapsed "
            <<std::chrono::duration<double>(now-start).count()<<" s\n"<<std::flush;
        };
        // remove: internal master bases without moving poles. They are found
        // on the exact graph, before any epsilon shift: at a numeric dimension
        // a moving pole is one more ordinary pole.
        if(moving_poles=="remove") {
          for(auto& [depth,basis]:diffexp::moving_poles::bases(graph,{},describe_moving_poles))
            numerical.physical_bases.insert_or_assign(depth,epsilon_shift.empty()?std::move(basis):
              diffexp::moving_poles::shift_epsilon(std::move(basis),diffexp::Rational(epsilon_shift),epsilon_order==0));
          // The scan for a smaller connection continues from these bases.
          numerical.scan_physical_bases=true;
        }
        // Experimental: expand around epsilon=shift; coefficient k of the
        // result then multiplies (epsilon-shift)^k.
        if(!epsilon_shift.empty()) {
          if(canonical_basis)throw std::invalid_argument("--epsilon-shift cannot be combined with --henn-basis");
          graph=diffexp::recursion::shift_epsilon(std::move(graph),diffexp::Rational(epsilon_shift),epsilon_order==0);
          // The Taylor tail of a sample follows its leaf digits like any other
          // run: the fit cannot resolve more than the leaf delivers.
          report["epsilon_shift"]=epsilon_shift;
        }
        report_moving_poles();
        const auto numerical_start=std::chrono::steady_clock::now();
        const auto raw_high=canonical_basis?std::max(0,diffexp::henn::needed_scalar_high(*canonical_basis,epsilon_order)):static_cast<int>(epsilon_order);
        diffexp::recursion::Evaluator evaluator(graph,numerical);auto result=evaluator.evaluate(raw_high);
        report["timings"].as_object()["numerical_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-numerical_start).count();
        if(canonical_basis)result=diffexp::henn::project_boundary(*canonical_basis,graph.nodes.front().requested,result,epsilon_order,numerical.working_bits);
        // The Taylor tail admitted per chart is not in the transported radii.
        // An uncertified result is widened by it: 2^6 charts' worth of the
        // tail, relative to the largest coefficient of each integral or, when
        // larger, to the largest term that cancels into it. This is
        // an estimate like the radii themselves; without it a sample asked
        // for few leaf digits would report far more digits than it has.
        if(!result.taylor_tail_certified) {
          const unsigned tail=numerical.adjoint.truncation_bits?numerical.adjoint.truncation_bits:static_cast<unsigned>(std::min<double>(
            std::min<double>(2.0*numerical.ordinary_order,static_cast<double>(numerical.working_bits)-32),
            std::ceil((numerical.leaf_digits+10.0)*std::log2(10.0))));
          for(std::size_t index=0;index<result.values.size();++index) {
            auto& row=result.values[index];
            mag_t scale,size;mag_init(scale);mag_init(size);
            for(const auto& value:row){acb_get_mag(size,value.raw());mag_max(scale,scale,size);}
            // The terms of a row before they cancel set the scale of its tail.
            if(index<result.term_log2.size() && std::isfinite(result.term_log2[index])){mag_set_ui_2exp_si(size,1,static_cast<slong>(std::ceil(result.term_log2[index])));mag_max(scale,scale,size);}
            mag_mul_2exp_si(scale,scale,6-static_cast<slong>(tail));
            for(auto& value:row)if(!value.is_zero())acb_add_error_mag(value.raw(),scale);
            mag_clear(scale);mag_clear(size);
          }
        }
        if(json_mode) {
          // Display enough decimal digits for the working binary precision;
          // Arb adds outward rounding to the printed interval radius.
          const long output_digits=std::max<long>(45,(static_cast<long>(numerical.working_bits)*30103+99999)/100000+3);
          boost::json::array rows;
          for(const auto& row:result.values) {
            boost::json::array powers;
            for(int k=result.low;k<=static_cast<int>(epsilon_order);++k)
              powers.push_back(diffexp::ball_json(row.at(k-result.low),output_digits));
            rows.push_back(std::move(powers));
          }
          report["epsilon_low"]=result.low;report["epsilon_high"]=epsilon_order;
          report["coefficient_order"]="component,epsilon";
          report["coefficients"]=std::move(rows);
          report["component_basis"]=canonical_basis?"canonical":"requested_scalar";
          report["ball_encoding"]="Arb decimal interval strings for real and imaginary parts";
          report["omitted_tails_certified"]=result.taylor_tail_certified;
          report["working_bits"]=numerical.working_bits;
          report["persistent_numerical_cache"]=numerical_cache;
        } else {
          if(!result.taylor_tail_certified)progress<<"Retained series values (printed radii exclude omitted endpoint and transport tails):\n";
          for(unsigned row=0;row<result.values.size();++row)for(int k=result.low;k<=static_cast<int>(epsilon_order);++k) {
            progress<<(canonical_basis?"Canonical component ":"Scalar ")<<row+1<<" epsilon^"<<k<<" = ";acb_printn(result.values[row][k-result.low].raw(),32,0);progress<<'\n';
          }
        }
        {
          // How the Feynman i0 was applied: which contour sides were stated,
          // which followed from positivity, and which were derived.
          const auto& prescription=evaluator.causal_prescription();
          const char* mode=prescription.assurance==diffexp::causal::Assurance::DerivedFromSymanzik?"derived":
            prescription.assurance==diffexp::causal::Assurance::EuclideanPhysicalDomain?"euclidean":"supplied";
          boost::json::array oriented,released;
          if(prescription.derive_orientations)for(const auto& point:evaluator.causal_points()) {
            if(!point.side){released.push_back(boost::json::object{{"level",static_cast<std::int64_t>(point.depth+1)},{"x",point.x}});continue;}
            oriented.push_back(boost::json::object{{"level",static_cast<std::int64_t>(point.depth+1)},{"x",point.x},{"side",point.side}});
          }
          if(prescription.derive_orientations) {
            progress<<"Causal contour derived from the Symanzik F of each level (F - i0): "<<oriented.size()<<" oriented real singular point"
              <<(oriented.size()==1?"":"s")<<", "<<released.size()<<" free.\n";
            for(const auto& entry:oriented)progress<<"  level "<<entry.at("level").as_int64()<<": x = "<<entry.at("x").as_string().c_str()
              <<" passed "<<(entry.at("side").as_int64()>0?"above":"below")<<" the real axis\n";
          }
          report["causal"]=boost::json::object{{"mode",mode},{"f_rim",prescription.f_rim}};
          if(prescription.derive_orientations) {
            report["causal"].as_object()["oriented_points"]=std::move(oriented);
            report["causal"].as_object()["free_points"]=std::move(released);
          }
        }
        const auto& stats=evaluator.statistics();
        report["ft_transport"]=boost::json::object{{"method","taylor"},{"direct_endpoints",stats.direct_endpoints},
          {"uncertainty",numerical.adjoint.uncertainty==diffexp::AdjointOptions::Uncertainty::probe?"probe":"map"}};
        report["ft_transport"].as_object()["accumulation"]=stats.weighted_levels?"weighted":"rows";
        report["ft_transport"].as_object()["weighted_levels"]=stats.weighted_levels;
        report["timings"].as_object()["ordinary_seconds"]=stats.ordinary_seconds;
        // Phases of numerical_seconds that precede transport. Disjoint from each
        // other and from ordinary_seconds.
        report["timings"].as_object()["endpoint_series_seconds"]=stats.endpoint_series_seconds;
        report["timings"].as_object()["endpoint_operator_seconds"]=stats.endpoint_operator_seconds;
        report["timings"].as_object()["basis_selection_seconds"]=stats.basis_selection_seconds;
        report["timings"].as_object()["leaf_seconds"]=stats.leaf_seconds;
        report["basis_selection"]=boost::json::object{{"scans",stats.basis_scans},
          {"selections",stats.basis_selections},{"cache_hits",stats.basis_cache_hits},{"endpoint_series_built",stats.endpoint_series_built},
          {"child_window_preflights",stats.demand_preflights},
          {"child_refinements",stats.refinements},{"structural_zero_coefficients",stats.structural_zero_coefficients},
          {"structural_zero_conflicts",stats.structural_zero_conflicts}};
        progress<<"Exact endpoint plans "<<stats.exact_plans<<"; series built "<<stats.endpoint_series_built<<"; numerical refinements "<<stats.refinements<<"; numerical cache hits "<<stats.cache_hits<<".\n"
          <<"Linear transport selections: adjoint "<<stats.adjoint_selections<<", factored "<<stats.factored_selections<<"; conditioning method fallbacks "<<stats.conditioning_method_fallbacks<<".\n"
          <<"Accumulation: "<<(stats.weighted_levels?"weighted, ":"rows, ")<<stats.weighted_levels<<" levels transported as weighted rows.\n"
          <<"Ordinary continuation: "<<stats.ordinary_checkpoints.loaded<<" saved arms resumed, "<<stats.ordinary_checkpoints.completed_reused
          <<" already complete; "<<stats.ordinary_checkpoints.saved<<" accepted charts saved; "<<stats.ordinary_checkpoints.skipped_oversize<<" transports too large to checkpoint.\n"
          <<"Child-window preflights "<<stats.demand_preflights<<"; local operator reuses "<<stats.local_operator_reuses<<".\n"
          <<"Basis scans "<<stats.basis_scans<<"; selections "<<stats.basis_selections<<".\n"
          <<"Seconds: endpoint series "<<stats.endpoint_series_seconds<<", endpoint operators "<<stats.endpoint_operator_seconds
          <<", basis selection "<<stats.basis_selection_seconds<<", leaf "<<stats.leaf_seconds<<", ordinary transport "<<stats.ordinary_seconds<<".\n"
          <<"Centered chart actions "<<conditioning.centered_charts<<"; conditioning subdivisions "<<conditioning.conditioning_subdivisions<<".\n"
          <<"Taylor recurrences: "<<diffexp::rational_circuit::compact_expansions.load()<<" on compact numbers, "
          <<diffexp::rational_circuit::ball_expansions.load()<<" in ball arithmetic.\n"
          <<(result.taylor_tail_certified?"Certified scalar-leaf coefficient enclosures.\n":"Numerical recursive series result; omitted endpoint and transport tails are not certified.\n");
        if(canonical_basis) {
          auto poles=diffexp::henn::audit_negative_poles(result,diffexp::Rational("1/100000000000000000000"),numerical.working_bits);
          progress<<"Canonical negative-pole audit at absolute tolerance 1e-20: "<<(poles.pass?"passed":"FAILED")
            <<" ("<<poles.checked_coefficients<<" coefficients checked).\n";
          report["negative_poles_pass"]=poles.pass;
          if(!poles.pass)throw std::runtime_error("canonical negative-pole audit failed");
        }
      }
      report["timings"].as_object()["total_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-preparation_start).count();
      if(json_mode)std::cout<<boost::json::serialize(report)<<'\n';
#else
      throw std::invalid_argument("native FIRE preparation currently requires a POSIX platform");
#endif
    } else if(command == "singular-endpoint") {
      if(argc!=2)throw std::invalid_argument("singular-endpoint takes no arguments; its JSON includes chart metadata and plot samples");
      std::cout<<boost::json::serialize(diffexp::singular_endpoint_example())<<'\n';
    } else if (command == "sunrise") {
      if(argc>3)throw std::invalid_argument("sunrise accepts an optional epsilon order");
      unsigned order=0;
      if(argc==3) {
        const std::string value=argv[2];auto parsed=std::from_chars(value.data(),value.data()+value.size(),order);
        if(parsed.ec!=std::errc{} || parsed.ptr!=value.data()+value.size())throw std::invalid_argument("epsilon order must be nonnegative integer");
      }
      diffexp::Jet::Ball::set_precision(384);auto result=diffexp::feynman::sunrise(order);
      for(unsigned k=0;k<=order;++k)std::cout<<"Sunrise epsilon^"<<k<<" = "<<result.at(k).real_midpoint(30)<<'\n';
      std::cout<<"Native two-level Feynman recursion; numerical series result with omitted Taylor tails not yet certified.\n";
    } else if(command == "banana4-oracle") {
      if(argc>3 || (argc==3 && std::string(argv[2])!="equal" && std::string(argv[2])!="unequal"))
        throw std::invalid_argument("banana4-oracle accepts equal or unequal");
      auto result=argc==3 && std::string(argv[2])=="unequal"?diffexp::oracle::banana4_unequal_bessel():diffexp::oracle::banana4_equal_bessel();
      std::cout<<"Independent certified Banana4 finite-part reference: ";acb_printn(result.value.raw(),40,0);
      std::cout<<"\nBessel integral with rigorous endpoint tails; this is not the Feynman-trick solver.\n";
    } else if (command == "bubble") {
      if(argc!=2)throw std::invalid_argument("bubble takes no arguments");
      auto result=diffexp::certified_feynman_bubble();
      for(unsigned k=0;k<result.coefficients.size();++k) {
        std::cout<<"Bubble epsilon^"<<k<<" = ";acb_printn(result.coefficients[k].raw(),30,0);std::cout<<'\n';
      }
      std::cout<<"Certified absolute radius <= "<<result.absolute_tolerance<<"; "<<result.charts<<" series charts.\n";
    } else if (command == "banana-unequal") {
      if((argc!=3 && argc!=5) || (argc==5 && std::string(argv[3])!="--route"))
        throw std::invalid_argument("banana-unequal requires ORIGINAL_DIR [--route momentum-first|mass-first]");
      const auto route=diffexp::original_banana_unequal_route(argc==5?argv[4]:"momentum-first");
      return diffexp::run_original_banana_unequal_route(argv[2],route);
    } else if (command == "banana-equal") {
      if(argc<3 || argc>7 || argc%2==0)throw std::invalid_argument("banana-equal requires ORIGINAL_DIR [--route real|contour] [--to ENDPOINT]");
      std::string route="real";diffexp::Rational endpoint(20);
      for(int i=3;i<argc;i+=2) {
        const std::string option=argv[i];
        if(option=="--route")route=argv[i+1];
        else if(option=="--to")endpoint=diffexp::Rational(argv[i+1]);
        else throw std::invalid_argument("unknown equal banana option");
      }
      if(route!="real" && route!="contour")
        throw std::invalid_argument("equal banana route must be real or contour");
      if(route=="contour" && endpoint!=diffexp::Rational(20))
        throw std::invalid_argument("the saved contour comparison ends at20; later endpoints use the real route");
      return route=="real"?diffexp::run_original_banana_equal_real(argv[2],endpoint):diffexp::run_original_banana_equal(argv[2]);
    } else if (command == "mpl") {
      if(argc>3 || (argc==3 && std::string(argv[2])!="--short"))throw std::invalid_argument("mpl accepts only --short");
      return diffexp::run_mpl(argc==2);
    } else if (command == "planar") {
      if(argc!=4)throw std::invalid_argument("planar requires ancillary directory and family");
      return diffexp::run_planar(argv[2],argv[3]);
    } else if (command == "series") {
      if(argc!=3)throw std::invalid_argument("series requires one JSON request file");
      return diffexp::run_series_file(argv[2]);
    } else if (command == "henn-nonplanar") {
      if(argc<3 || argc>4)throw std::invalid_argument("henn-nonplanar requires ancillary directory and optional Taylor order");
      unsigned order=140;
      if(argc==4) {
        const std::string value=argv[3];
        auto parsed=std::from_chars(value.data(),value.data()+value.size(),order);
        if(parsed.ec!=std::errc{} || parsed.ptr!=value.data()+value.size())
          throw std::invalid_argument("Henn Taylor order must be a positive integer");
      }
      return diffexp::run_henn_nonplanar(argv[2],order);
    } else if (command == "backend-info") {
      boost::json::object information{{"standalone_cpp",true},{"librarylink",false}};
      information["package_version"]="2.2.1";information["native_transport"]=true;
      information["feynman_configuration_schema"]="DiffExp.FeynmanFamily/v1";
#if defined(__unix__) || defined(__APPLE__)
      information["ibp_solver"]=true;information["ibp_default_provider"]="ibp-solver";information["ibp_prime_bits"]=61;
#else
      information["ibp_solver"]=false;
#endif
      information["transport_schema"]="DiffExp.Transport/v1";
      std::cout<<boost::json::serialize(information)<<'\n';
    } else if (command == "--help") {
      std::cout << R"(Usage: diffexp COMMAND [ARGUMENTS]

General
  --version | --help | backend-info
  series FILE|-                   exact regular-series coefficients of a JSON request
  transport FILE|-                numerical transport of a JSON request (docs/transport-protocol.md)

Feynman-trick recursion
  family-template NAME            print an editable family configuration
  prepare CONFIGURATION [OPTIONS] exact IBP preparation only
  ft CONFIGURATION [OPTIONS]      full calculation
  Use --json with prepare/ft for one JSON response; progress is written to stderr.

  --cache DIRECTORY               exact level closures and numerical checkpoints
                                  (default: diffexp-exact-cache in the temporary directory)
  --epsilon-order N               0..100 (default: the configuration's epsilon_order)
  --epsilon-strategy expand|sample
                                  expand (default) carries the Laurent coefficients through the
                                  recursion. sample evaluates the integral at shifted dimensions
                                  in parallel processes and fits the coefficients: the route for
                                  families whose expanded recursion is too deep
                                  (docs/performance.md). Its radii are estimates.
  --samples N                     even, 4..200 (default 2*loops + epsilon order + 5, at least 8)
  --sample-step Q                 samples at epsilon = +-k*Q (default 1/400)
  --sample-jobs N                 samples run at a time (default: one per idle hardware thread)
  --sample-poles N                order of the pole at epsilon = 0 (default: detected, at most
                                  2*loops). Samples ask the leaf for min(leaf digits, 16), or 12
                                  from six levels on, unless --leaf-digits is given.
  --endpoint-order N              1..1000: order of the exact series at the endpoints
  --overlap-ratio R               4..65536: the series is summed at most 1/R of the way to the
                                  nearest other singular point, where the transport starts.
                                  Defaults: R=64, N=22 (R=16, N=32 with --epsilon-shift), a
                                  2^-128 tail; with only N given, R follows from it.
  --ordinary-order N              8..1000 (default 80)
  --precision-bits N              64..1000000 (default 384; 64..4096 with --henn-basis)
  --leaf-digits N                 1..100000 (default 28)
  --method adjoint|factored|auto  linear transport of boundary data (default adjoint)
  --accumulation weighted|rows    how the levels are combined (adjoint method). weighted
                                  (default): from the top, every level transports one row per
                                  requested integral of the family, weighted by what the level
                                  above passes down, and the leaf values enter last. rows: from
                                  the leaf, every level transports one row per integral its
                                  parent needs. More than 4 requested integrals (--henn-basis)
                                  and the other methods are always evaluated by rows
                                  (docs/performance.md)
  --uncertainty probe|map         how a transport carries the uncertainty it inherits. probe
                                  (default): midpoints, with two probe rows (one beside a single
                                  row) following the size of the uncertainty, an estimate; steps
                                  are judged on the midpoints. map: the chart map applied to the
                                  intervals, a rigorous enclosure, one expansion per master on
                                  every chart that needs it
  --transport-bits N              relative Taylor tail admitted per chart, as a power of two
                                  (default min(2*ordinary order, precision bits-32))
  --basis-scan-seconds N          0..60 (default 30): per level, the search for an internal master
                                  basis with a smaller exact connection (cheaper transport) that does
                                  not widen the child epsilon window; 0 disables it
  --no-numerical-cache keeps only the exact IBP level closures under --cache. It disables the
                                  persisted internal-basis choice (physical-bases) and the
                                  ordinary-transport checkpoints (ordinary-transport)
  --henn-basis FILE               project onto a published canonical basis (Henn X0 geometry only)
  --epsilon-shift Q               experimental: expand around epsilon=Q (a rational); coefficient k
                                  then multiplies (epsilon-Q)^k, and order 0 is the integral in
                                  d0-2Q dimensions (docs/performance.md)
  --moving-poles keep|remove|exchange
                                  per level, the denominator factors of the connection that mix
                                  x with epsilon (docs/moving-poles.md; default remove for ft,
                                  exchange for epsilon samples). remove (ft): in an internal basis, a
                                  master replaced by a cached integral or an exact combination;
                                  the basis scan of such a level starts from this basis.
                                  exchange (prepare/ft): the level's masters themselves, a master
                                  replaced by a cached integral of its sector; the levels below
                                  are prepared for the new masters
  These ranges are checked before the IBP preparation starts. Ordinary transport is by Taylor charts;
  --ft-transport, --ft-transport-digits and --method values were removed. Endpoint series are never
  cached: recomputing them is faster than verifying a stored copy.

IBP preparation (prepare/ft)
  --ibp-provider auto|ibp-solver|fire|fire-modular
                                  auto selects the built-in ibp-solver unless --fire names a FIRE
                                  executable
  --ibp-dots N --ibp-numerators N the seeds every sector of the built-in solver starts from
                                  (0..8, default 1 and 1). They are grown where needed: a requested
                                  integral asks for the seeds that generate it, and an integral
                                  left unreduced outside the seeds of its sector widens that sector
  --ibp-seeds demand|sector|subsectors|all
                                  the sectors that receive in advance the seeds that generate a
                                  requested integral: none, widening alone decides (demand, the
                                  default); its own; its own and its subsectors; or all
  --merge-selection auto|first    which two propagators each level merges, unless the family file
                                  lists preparation.merges or anchors: auto (default) takes the
                                  pair whose merged graph has the fewest master integrals,
                                  counted from its Symanzik polynomials with Singular for every
                                  pair; without Singular, by IBP eliminations for large levels.
                                  first merges the first two
  --boundary anchor|endpoint|auto where a level takes its boundary data. anchor: at an interior value of
                                  its merge parameter (1/3 unless singular), where the next level is
                                  the merged family at that value. auto (default): at x = 0 or 1
                                  where that endpoint qualifies; the next level is then the family
                                  with one of the two merged propagators removed, at its original
                                  kinematics. An endpoint qualifies if the family's U and F have
                                  positive monomials, the two propagators share their loop momentum,
                                  the one that remains is massive, and the level's equation shows
                                  that the masters' values at the endpoint fix its analytic solution
                                  (docs/performance.md); such a level uses the adjoint method.
                                  endpoint: as auto, an error if no level qualifies. Explicit
                                  anchors, stated causal sides and --henn-basis keep the anchors
  --ibp-max-degree N              4..256 (default 16, doubled automatically when a reduction meets it): largest total degree, in the merge
                                  parameter and the dimension together, of a reduction
                                  coefficient the built-in solver reconstructs
  --ibp-threads N                 threads of the built-in solver (default: the hardware threads
                                  idle at start-up). Trace replays of announced sample points run
                                  in batches, and complete eliminations at points fixed in advance
                                  start early on free threads. Results do not depend on N.
  --ibp-validation anchored|strict
                                  anchored (default): held-out points eliminate the pilot's
                                  dependency subset of original equations, with one complete
                                  elimination per attempt; strict eliminates every equation at
                                  every point
  --fire-seconds N                one IBP elimination or FIRE call (default 3600)
  --level-seconds N               one recursion level (default 86400)
  --total-seconds N               the whole preparation (default 86400)
  --fire PATH | --fire-prime PATH FIRE executable | FIRE7p for fire-modular (durable samples under
                                  --cache, no symbolic fallback)
  --fire-threads N --fire-simplifier-threads N
                                  FIRE workers (default 1 each)
  --fire-memory-mib N             memory budget of one reduction, either provider (default 4096)
  The same budgets can be set per family as preparation.fire_seconds, level_seconds and
  total_seconds. Verified level closures under --cache are reused whatever
  --ibp-provider/--ibp-dots/--ibp-numerators say; such reuse is reported on stderr and in the JSON
  field ibp_cache.

Built-in checks and original DiffExp examples
  bubble | sunrise [EPSILON_ORDER] | banana4-oracle [equal|unequal] | mpl [--short]
  singular-endpoint               JSON chart, sectors, endpoint and plot samples
  henn-nonplanar DATA [ORDER] | planar DATA 1loop|zmz|mzz|zzz
  banana-equal ORIGINAL_DIR [--route real|contour] [--to ENDPOINT]
                                  default route real (local +i0 singular matching) and endpoint 20;
                                  the real route supports endpoints >= 20, including the historical 32
  banana-unequal ORIGINAL_DIR [--route momentum-first|mass-first]
                                  default momentum-first
)";
    } else {
      std::cerr << "Unknown command: " << command << '\n';
      return 2;
    }
  } catch (const std::exception& e) {
    if(json_mode)std::cout<<boost::json::serialize(boost::json::object{{"schema","DiffExp.Error/v1"},{"status","error"},{"message",e.what()}})<<'\n';
    std::cerr << "diffexp: " << e.what() << '\n';
    return 1;
  }
}
