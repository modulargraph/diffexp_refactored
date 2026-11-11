#pragma once
#include "diffexp/fire_modular.hpp"
#include <ibp/trace.hpp>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

namespace diffexp::ibp_solver {
// The pinned upstream core works in the exact denominator coordinates owned by
// DiffExp. No momentum reconstruction, renamed family, or new ISP completion.
struct Options {
  std::filesystem::path cache_directory;
  unsigned dots=1,numerators=2,max_degree=16,max_variable_degree=64,max_primes=32,probe_timeout_seconds=20;
  std::size_t max_samples_per_prime=512;
  // Ordinary samples are recomputed from their deterministic points in a few
  // milliseconds, so only every n-th is made durable; held-out validation
  // samples, which the completion record replays, are always retained.
  std::size_t sample_checkpoint_interval=32;
  bool degree_slices=false,dots_first=false,target_sector_seeds=false;
  // Fresh primal eliminations certify a reconstruction (see Sampler). By
  // default each eliminates ALL generated equations. With selected_source_rows
  // only the complete pilot does; later ones eliminate the pilot's dependency
  // subset of ORIGINAL equations. complete_validation_anchor additionally
  // eliminates all equations at the first held-out point of every validation
  // attempt, i.e. at a fresh prime and point.
  bool selected_source_rows=false,complete_validation_anchor=false;
  // Zero derives the arithmetic-trace node budget from the memory budget.
  std::size_t max_trace_nodes=0;
  // Threads for eliminations whose points are known before their results are
  // needed (see Sampler::announce and Sampler::anticipate); 0 uses one per
  // hardware thread. Reductions, sample points and counts do not depend on it.
  unsigned threads=0;
  std::optional<std::size_t> closure_parameter;
  unsigned max_planning_passes=6;
  std::function<void(const std::string&)> progress;
};
struct Statistics {
  std::size_t probes=0,equations=0,templates=0,full_solves=0,trace_replays=0,trace_fallbacks=0,trace_disabled=0;
  std::size_t selected_solves=0,selected_equations=0,selection_fallbacks=0;
  double generation_seconds=0,elimination_seconds=0,full_solve_seconds=0,trace_learning_seconds=0,trace_replay_seconds=0;
  double selected_solve_seconds=0;
  // Work done on other threads: replays served from a batch, and complete
  // eliminations started ahead and those of them that were used. Seconds of
  // such work are summed over threads.
  std::size_t batched_replays=0,early_solves=0,early_solves_used=0;
};
// Descending primes below 2^61; two beyond the reconstruction primes validate.
inline constexpr unsigned prime_count=64;
inline modular::Word prime(unsigned index) {
  if(index<1||index>prime_count)throw std::invalid_argument("IBP solver prime index must be 1.."+std::to_string(prime_count));
  static const auto values=[] {std::array<modular::Word,prime_count> out{};out[0]=2305843009213693951ULL;
    for(unsigned i=1;i<out.size();++i){auto p=out[i-1]-2;while(!n_is_prime(p))p-=2;out[i]=p;}return out;}();
  return values[index-1];
}
class Sampler {
 public:
  Sampler(ibp::PropagatorBasis basis,Exact dimension,Options options)
    :basis_(std::move(basis)),dimension_(std::move(dimension)),options_(std::move(options)),
     threads_(options_.threads?options_.threads:std::max(1u,std::thread::hardware_concurrency())) {
    if(!basis_.space.loops||basis_.space.loops>4||basis_.space.size()>16||basis_.physical_count>12||!basis_.physical_count||
       basis_.denominators.size()!=basis_.space.size()||basis_.physical_count>basis_.denominators.size())
      throw std::invalid_argument("IBP solver supports at most 4 loops, 16 scalar products and 12 physical denominators; select FIRE for larger systems");
    if(options_.dots>8||options_.numerators>8)throw std::invalid_argument("IBP solver seed powers must be 0..8");
    if(threads_>256)throw std::invalid_argument("IBP solver threads must be 0..256");
    contractions_=ibp::contractions(basis_);
    for(unsigned i=0;i<basis_.space.loops;++i)for(unsigned v=0;v<basis_.space.loops+basis_.space.externals();++v)trace_.push_back(i==v);
    zero_sectors_.resize(std::size_t(1)<<basis_.physical_count);
    for(const auto& a:fire::free_loop_sectors(basis_))zero_sectors_[sector(a)]=true;
    input_geometry_.loops=basis_.space.loops;input_geometry_.externals=basis_.space.externals();
    input_geometry_.physical=basis_.physical_count;input_geometry_.n=basis_.space.size();
    input_geometry_.trace=trace_;input_geometry_.zero_sectors=zero_sectors_;
    std::map<std::string,std::uint32_t> ids;
    auto input=[&](const Exact& value){if(value.is_zero())return ::ibp::InputGeometry::zero;
      auto [it,inserted]=ids.try_emplace(value.str(),inputs_.size());if(inserted){inputs_.push_back(value);input_geometry_.constant_inputs.push_back(value.is_rational());}return it->second;};
    // Even a constant zero dimension needs a valid explicit input slot.
    inputs_.push_back(dimension_);input_geometry_.constant_inputs.push_back(dimension_.is_rational());ids.emplace(dimension_.str(),0);input_geometry_.dimension=0;
    for(const auto& row:contractions_){input_geometry_.contractions.emplace_back();for(const auto& a:row){std::vector<std::uint32_t> c;for(const auto& v:a.linear)c.push_back(input(v));c.push_back(input(a.constant));input_geometry_.contractions.back().push_back(std::move(c));}}
    input_geometry_.inputs=inputs_.size();
  }
  const Statistics& statistics()const{return stats_;}
  // Whether announcements have any effect: with one thread they are ignored.
  bool works_ahead()const{return threads_>1;}
  struct Point {modular::Word modulus=0;std::vector<modular::Word> coordinates;bool operator==(const Point&)const=default;};
  // Working ahead. Neither announcement changes a result, a statistic other
  // than the three above, or the state a later call sees: what another thread
  // computed is used only where this thread would have computed exactly the
  // same, and is dropped otherwise. Each replaces the previous one of its kind.
  //
  // announce: the points of the next ordinary samples of one prime, in order.
  // Once the trace of that prime is learned and checked, a replay depends on
  // the trace and the point alone, so the announced points are replayed
  // together when the first of them is requested.
  void announce(std::vector<std::vector<modular::Word>> points,modular::Word modulus) {
    ahead_.clear();ahead_modulus_=modulus;if(threads_<2)return;
    for(auto& point:points)ahead_.push_back({std::move(point),std::nullopt,0,0,0});
  }
  // anticipate: points at which a complete elimination is likely to be
  // requested, in order of likelihood, within the given seconds. It depends on
  // the generated equations and the point, not on the targets, the pilot or a
  // trace, so it starts on a free thread as soon as the equations exist and
  // is adopted if they are still the same. One already started for a point
  // that is announced again is kept.
  void anticipate(std::vector<Point> points,double seconds=0) {
    std::vector<Expected> next;const auto until=std::chrono::steady_clock::now()+std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(std::clamp(seconds,0.0,1e6)));
    if(threads_>1)for(auto& point:points){
      auto old=std::find_if(expected_.begin(),expected_.end(),[&](const Expected& expected){return expected.at==point;});
      next.push_back({std::move(point),old==expected_.end()?nullptr:std::move(old->task),until});
    }
    expected_=std::move(next);
  }
  // `independent` solves by fresh primal elimination of original rows, without
  // the arithmetic trace. `complete` forces that elimination to use every
  // generated equation even when a dependency subset is available.
  //
  // What a subset elimination establishes: its rows are combinations of original
  // IBP equations at this point, hence identities of the complete system too.
  // If the terminal integrals are free columns of the complete system at this
  // point, uniqueness of the normal form makes them equal to the complete
  // system's rows. What only a complete elimination establishes: that no
  // omitted equation reduces a terminal integral here (rank and minimality).
  fire::Result operator()(const std::vector<ibp::Integral>& requested,const std::vector<modular::Word>& point,
      modular::Word modulus,const fire::Options& limits,bool independent=false,bool complete=false) {
    fire::Result out;
    try {
      if(!limits.timeout_seconds||!limits.memory_bytes||requested.empty()||requested.size()>10000)
        throw std::invalid_argument("IBP solver requires bounded nonempty requests");
      const auto begin=std::chrono::steady_clock::now();
      auto remaining=[&]{auto s=limits.timeout_seconds-std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        if(s<=0)throw std::runtime_error("IBP solver probe time budget exceeded");return s;};
      ::ibp::Field field(modulus);auto geometry=input_geometry_;
      for(const auto& a:limits.zero_sectors){fire_batch::validate_index(a,basis_);geometry.zero_sectors[sector(a)]=true;}
      std::vector<::ibp::Word> values;for(const auto& input:inputs_)values.push_back(modular::evaluate(input,point,modulus));
      ::ibp::SeedOptions seeds;seeds.dots=options_.dots;seeds.numerators=options_.numerators;
      std::vector<::ibp::Integral> targets;
      for(const auto& a:requested){fire_batch::validate_index(a,basis_);::ibp::Integral target;unsigned dots=0,nums=0;
        for(unsigned i=0;i<a.size();++i){target.powers[i]=a[i];if(a[i]>0)dots+=a[i]-1;else nums-=a[i];}
        // A seed with one fewer dot can produce the requested raised integral.
        seeds.dots=std::max(seeds.dots,dots?dots-1:0);seeds.numerators=std::max(seeds.numerators,nums);targets.push_back(target);
      }
      std::vector<std::pair<unsigned,unsigned>> sector_bounds;
      if(options_.target_sector_seeds){
        sector_bounds.resize(std::size_t(1)<<basis_.physical_count,{options_.dots,options_.numerators});
        for(const auto& target:requested){unsigned dots=0,nums=0;for(int power:target){if(power>0)dots+=power-1;else nums-=power;}
          const auto mask=sector(target);for(std::size_t sub=mask;;sub=(sub-1)&mask){
            auto& bound=sector_bounds[sub];bound.first=std::max(bound.first,dots?dots-1:0);bound.second=std::max(bound.second,nums);if(!sub)break;}
}
      }
      seeds.max_terms=std::min<std::size_t>(10000000,limits.memory_bytes/64);seeds.seconds=remaining();
      if(!seeds.max_terms)throw std::invalid_argument("IBP solver memory budget too small");
      if(!plan_||plan_->program->program.targets!=targets||plan_->zero_sectors!=geometry.zero_sectors||plan_->dots!=seeds.dots||plan_->numerators!=seeds.numerators||plan_->max_terms!=seeds.max_terms){
        auto fresh=std::make_unique<Plan>();fresh->program=std::make_shared<::ibp::ParametricProgram>(::ibp::ParametricProgram::compile(geometry,seeds,targets,options_.dots_first?::ibp::ParametricOrdering::dots_first: ::ibp::ParametricOrdering::total_degree,sector_bounds));
        fresh->zero_sectors=geometry.zero_sectors;fresh->dots=seeds.dots;fresh->numerators=seeds.numerators;fresh->max_terms=seeds.max_terms;
        // New targets usually leave the generated equations as they were, and
        // with them every complete elimination started ahead.
        fresh->system=plan_&&same_equations(*plan_->program,*fresh->program)?plan_->system:++systems_;
        plan_=std::move(fresh);++state_;++stats_.templates;stats_.equations+=plan_->program->program.equations.size();
      }
      auto& parametric=*plan_->program;auto& program=parametric.program;
      const auto generated=std::chrono::steady_clock::now();
      ::ibp::SolveOptions solve;solve.seconds=remaining();solve.max_fill=std::min<std::size_t>(12000000,limits.memory_bytes/64);
      start_anticipated(solve,{modulus,point});double elsewhere=0;
      // Reconstruction visits primes in order and never returns to one, so a
      // trace is dropped once ordinary sampling has moved to another prime.
      if(!independent&&plan_->sampling_modulus!=modulus){for(auto& [other,state]:plan_->primes)if(other!=modulus){state.trace.reset();state.workspace={};}plan_->sampling_modulus=modulus;++state_;}
      std::vector<::ibp::Row> rows;auto& learned=plan_->primes[modulus];bool replayed=false,learning_overran=false;
      // A replay that neither relearns nor is checked leaves no state behind:
      // an announced point is replayed with the announced points after it.
      if(!independent&&learned.trace&&learned.checked)if(auto* slot=announced(point,modulus)){
        if(!slot->result||slot->state!=state_)replay_announced(requested,modulus,field,*learned.trace,program,limits);
        if(slot->result&&slot->state==state_){
          out=std::move(*slot->result);++stats_.trace_replays;++stats_.batched_replays;stats_.trace_replay_seconds+=slot->replay_seconds;
          ++stats_.probes;stats_.generation_seconds+=std::chrono::duration<double>(generated-begin).count();stats_.elimination_seconds+=slot->seconds;
          ahead_.erase(ahead_.begin());return out;
        }
      }
      if(!independent&&learned.trace){const auto replay_begin=std::chrono::steady_clock::now();try{rows=learned.trace->evaluate_inputs(values,field,learned.workspace);++stats_.trace_replays;replayed=true;}
        catch(const std::domain_error&){learned.trace.reset();++state_;++stats_.trace_fallbacks;}stats_.trace_replay_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-replay_begin).count();}
      // Held-out samples use fresh primal elimination of original IBP rows,
      // independent of the arithmetic trace. A dependency subset proves the
      // same rowspace identities, but does not certify full-system minimality.
      // Check the first changed point too, before using a learned trace alone.
      if(independent||!replayed||!learned.checked){
        std::unique_ptr<::ibp::Solver> owned;std::shared_ptr<Elimination> early;::ibp::Solver* solver=nullptr;
        std::vector<::ibp::Row> full;std::optional<std::vector<std::uint32_t>> dependencies;
        auto terminals=[](const auto& reduced){std::set<::ibp::Column> out;for(const auto& row:reduced)for(const auto& entry:row)out.insert(entry.column);return out;};
        auto primal=[&](bool selected){
          auto started=std::chrono::steady_clock::now();double seconds=0;
          const auto* sources=selected?&*plan_->selection:nullptr;dependencies.reset();early.reset();
          if(!selected&&(early=adopt({modulus,point},solve.max_fill,remaining()))){
            // Its seconds count once: not those this thread spent waiting for it.
            elsewhere+=std::max(0.0,early->seconds-std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count());
            started=std::chrono::steady_clock::now();
            if(early->target_columns==program.target_columns){full=std::move(early->rows);dependencies=std::move(early->selection);}
            // Other targets are reduced here. The solver's clock has run since
            // its elimination began on the other thread: if it expires only
            // now, the elimination is repeated with a budget of its own.
            else try{early->program.target_columns=program.target_columns;full=early->solver->targets();dependencies=early->solver->selection();}
              catch(const std::runtime_error&){early.reset();}
            if(early){solver=&*early->solver;seconds=early->seconds;++stats_.early_solves_used;}
          }
          if(!early){
            // Rebind ALL original rows on fallback, including those left at an
            // earlier point or prime by the preceding subset solve.
            started=std::chrono::steady_clock::now();parametric.bind(values,field,sources);solve.seconds=remaining();
            owned=std::make_unique<::ibp::Solver>(program,field,solve);owned->solve(0,sources);solver=owned.get();full=solver->targets();
          }
          seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
          if(selected){++stats_.selected_solves;stats_.selected_equations+=sources->size();stats_.selected_solve_seconds+=seconds;}
          else {++stats_.full_solves;stats_.full_solve_seconds+=seconds;}
        };
        auto selection=[&]()->const std::vector<std::uint32_t>&{if(!dependencies)dependencies=solver->selection();return *dependencies;};
        bool selected=options_.selected_source_rows&&plan_->selection.has_value()&&!complete;primal(selected);
        if(selected&&terminals(full)!=plan_->terminal_columns){++stats_.selection_fallbacks;selected=false;primal(false);}
        // A complete anchor compares against the pilot; it does not replace the
        // pilot's subset, which the reconstruction samples were drawn from.
        if(options_.selected_source_rows&&!selected&&!(complete&&plan_->selection)){
          plan_->selection=selection();plan_->terminal_columns=terminals(full);
          if(options_.progress)options_.progress("full IBP pilot: "+std::to_string(program.equations.size())+" equations, "+std::to_string(plan_->selection->size())+" original rows needed for target identities");
        }
        if(replayed&&rows!=full){learned.trace.reset();++stats_.trace_fallbacks;replayed=false;}
        rows=std::move(full);
        if(!independent){
          ++state_;
          if(replayed)learned.checked=true;
          else if(!learned.disabled){const auto learning_begin=std::chrono::steady_clock::now();::ibp::ArithmeticTrace::Options tracing;tracing.seconds=remaining();tracing.max_fill=solve.max_fill;
            tracing.max_nodes=options_.max_trace_nodes?options_.max_trace_nodes:std::min<std::size_t>(std::numeric_limits<std::uint32_t>::max(),limits.memory_bytes/128);
            // The primal rows above are already a valid reduction. A trace that
            // does not fit its node, fill or time budget only costs the replay
            // shortcut: this prime continues with primal elimination alone.
            try{learned.trace=::ibp::ArithmeticTrace::learn_parametric(parametric,field,values,selection(),tracing);}
            catch(const std::length_error&){learned.disabled=true;}
            catch(const std::runtime_error&){learned.disabled=true;learning_overran=true;}
            if(learned.disabled){learned.trace.reset();++stats_.trace_disabled;if(options_.progress)options_.progress("arithmetic trace exceeds its budget at this prime; continuing with primal elimination");}
            else if(learned.trace->evaluate_inputs(values,field,learned.workspace)!=rows)throw std::runtime_error("IBP parametric trace disagrees with primal learning reduction");
            learned.checked=false;stats_.trace_learning_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-learning_begin).count();}
        }
      }
      // Time spent on an abandoned trace is not charged to the valid reduction.
      if(!learning_overran)remaining();
      emit(out,rows,program,requested);
      ++stats_.probes;stats_.generation_seconds+=std::chrono::duration<double>(generated-begin).count();
      stats_.elimination_seconds+=elsewhere+std::chrono::duration<double>(std::chrono::steady_clock::now()-generated).count();
    }catch(const std::domain_error& e){
      // Raised only while evaluating the input coefficients: a pole at this
      // point (or modulo this prime). Elimination itself never inverts zero.
      out.success=false;out.exceptional_point=true;out.reason=e.what();
    }catch(const std::exception& e){out.success=false;out.reason=e.what();}
    return out;
  }
 private:
  std::size_t sector(const ibp::Integral& a)const {std::size_t mask=0;for(unsigned i=0;i<basis_.physical_count;++i)if(a[i]>0)mask|=std::size_t(1)<<i;return mask;}
  // The reduction a probe returns for the rows of its targets.
  void emit(fire::Result& out,const std::vector<::ibp::Row>& rows,const ::ibp::Program& program,const std::vector<ibp::Integral>& requested)const {
    auto indices=[&](const ::ibp::Integral& a){ibp::Integral r;for(unsigned i=0;i<input_geometry_.n;++i)r.push_back(a.powers[i]);return r;};
    std::set<ibp::Integral> terminal;
    for(unsigned i=0;i<rows.size();++i){ibp::Relation row;for(const auto& term:rows[i]){auto a=indices(program.integrals[term.column]);terminal.insert(a);ibp::add(row,a,dimension_.parse(std::to_string(term.value)));}out.reductions.emplace(requested[i],std::move(row));}
    for(const auto& a:terminal){auto [it,inserted]=out.reductions.emplace(a,ibp::Relation{{a,dimension_.constant(1)}});
      if(!inserted&&it->second!=ibp::Relation{{a,dimension_.constant(1)}})throw std::runtime_error("IBP solver returned a nonterminal output basis");}
    out.success=true;fire_batch::validate(out,basis_,dimension_,requested);
  }
  struct Learned {std::optional<::ibp::ArithmeticTrace> trace;::ibp::ArithmeticTrace::Workspace workspace;bool checked=false,disabled=false;};
  // `system` names the generated equations: plans that differ in targets alone share it.
  struct Plan {std::shared_ptr<::ibp::ParametricProgram> program;std::vector<bool> zero_sectors;unsigned dots=0,numerators=0;std::size_t max_terms=0,system=0;std::map<modular::Word,Learned> primes;modular::Word sampling_modulus=0;
    std::optional<std::vector<std::uint32_t>> selection;std::set<::ibp::Column> terminal_columns;};
  static bool same_equations(const ::ibp::ParametricProgram& a,const ::ibp::ParametricProgram& b) {
    if(a.input_count!=b.input_count||a.constant_inputs!=b.constant_inputs||a.program.integrals!=b.program.integrals||a.program.equations.size()!=b.program.equations.size())return false;
    for(std::size_t i=0;i<a.program.equations.size();++i){const auto& x=a.program.equations[i];const auto& y=b.program.equations[i];if(x.size()!=y.size())return false;
      for(std::size_t j=0;j<x.size();++j)if(x[j].column!=y[j].column||a.coefficients[i][j].constant!=b.coefficients[i][j].constant||a.coefficients[i][j].terms!=b.coefficients[i][j].terms)return false;}
    return true;
  }
  // An announced ordinary sample. A result is current while state_, which
  // counts every change of plan or trace, is the one it was computed under.
  struct Ahead {std::vector<modular::Word> point;std::optional<fire::Result> result;std::uint64_t state=0;double replay_seconds=0,seconds=0;};
  // The slot of a requested point; announced points before it were not used.
  Ahead* announced(const std::vector<modular::Word>& point,modular::Word modulus) {
    if(modulus!=ahead_modulus_)return nullptr;
    auto it=std::find_if(ahead_.begin(),ahead_.end(),[&](const Ahead& slot){return slot.point==point;});
    if(it==ahead_.end())return nullptr;ahead_.erase(ahead_.begin(),it);return &ahead_.front();
  }
  // Replays a bounded stretch of announced points, the requested one first,
  // each worker in its own workspace. A point that fails in any way keeps no
  // result: its own request then takes the serial path and reports the failure
  // exactly as before. The requesting probe's budget bounds the whole batch.
  void replay_announced(const std::vector<ibp::Integral>& requested,modular::Word modulus,const ::ibp::Field& field,const ::ibp::ArithmeticTrace& trace,
      const ::ibp::Program& program,const fire::Options& limits) {
    using clock=std::chrono::steady_clock;auto since=[](clock::time_point t){return std::chrono::duration<double>(clock::now()-t).count();};
    std::vector<Ahead*> batch;const std::size_t stretch=std::max<std::size_t>(32,8*threads_);
    for(auto& slot:ahead_){if(batch.size()==stretch)break;if(!slot.result||slot.state!=state_){slot.result.reset();batch.push_back(&slot);}}
    const auto began=clock::now();std::atomic<std::size_t> next{0};
    auto work=[&]{::ibp::ArithmeticTrace::Workspace workspace;
      for(;;){const auto i=next.fetch_add(1);if(i>=batch.size()||since(began)>=limits.timeout_seconds)break;auto& slot=*batch[i];const auto start=clock::now();
        try{std::vector<::ibp::Word> values;for(const auto& input:inputs_)values.push_back(modular::evaluate(input,slot.point,modulus));
          auto rows=trace.evaluate_inputs(values,field,workspace);slot.replay_seconds=since(start);if(slot.replay_seconds>=limits.timeout_seconds)continue;
          fire::Result reduced;emit(reduced,rows,program,requested);slot.seconds=since(start);slot.state=state_;slot.result=std::move(reduced);
        }catch(...){}
      }};
    const auto busy=std::min<std::size_t>(running(),threads_-1);
    std::vector<std::jthread> pool;
    try{for(std::size_t i=1;i<std::min<std::size_t>(batch.size(),threads_-busy);++i)pool.emplace_back([&]{work();flint_cleanup();});}catch(const std::system_error&){}
    work();
  }
  // A complete elimination on another thread. It owns everything it reads
  // while it runs except the immutable equation structure it was started with.
  struct Elimination {
    std::size_t system=0,max_fill=0;::ibp::Field field;::ibp::Program program;std::optional<::ibp::Solver> solver;
    // Reduced for the targets known when it started, on its own clock.
    std::vector<::ibp::Column> target_columns;std::vector<::ibp::Row> rows;std::vector<std::uint32_t> selection;double seconds=0;
    std::mutex mutex;std::condition_variable ready;bool finished=false,solved=false;
    explicit Elimination(modular::Word modulus):field(modulus){}
  };
  struct Expected {Point at;std::shared_ptr<Elimination> task;std::chrono::steady_clock::time_point until;};
  std::size_t running() {
    std::erase_if(tasks_,[](auto& task){const std::lock_guard lock(task.second->mutex);return task.second->finished;});return tasks_.size();
  }
  // Starts anticipated eliminations of the present equations, in the order
  // announced, while a thread is free; the point being requested now is left
  // to this thread. At most `lookahead` wait to be asked for, since the later
  // an announced point the less likely its use. An elimination that fails is
  // repeated by the request it was meant for, which then reports the failure.
  // The solver's clock covers the time a result may wait, so that reducing
  // other targets later does not overrun it; adopt() applies the requesting
  // probe's own budget to the elimination.
  static constexpr std::size_t lookahead=3;
  void start_anticipated(const ::ibp::SolveOptions& probe,const Point& present) {
    std::size_t waiting=0;running();
    for(auto& expected:expected_){
      if(expected.task&&expected.task->system==plan_->system){++waiting;continue;}
      if(expected.at==present)continue;
      if(waiting>=lookahead||running()+1>=threads_)break;
      std::vector<::ibp::Word> values;
      try{for(const auto& input:inputs_)values.push_back(modular::evaluate(input,expected.at.coordinates,expected.at.modulus));}catch(const std::exception&){continue;}
      auto solve=probe;solve.seconds=std::max(solve.seconds,std::chrono::duration<double>(expected.until-std::chrono::steady_clock::now()).count());
      auto task=std::make_shared<Elimination>(expected.at.modulus);task->system=plan_->system;task->max_fill=solve.max_fill;task->target_columns=plan_->program->program.target_columns;
      try{tasks_.emplace_back(std::jthread([task,source=std::shared_ptr<const ::ibp::ParametricProgram>(plan_->program),values=std::move(values),solve]{
          const auto started=std::chrono::steady_clock::now();bool solved=false;
          try{const auto& equations=source->program.equations;task->program.integrals.resize(source->program.integrals.size());task->program.equations.resize(equations.size());
            for(std::size_t i=0;i<equations.size();++i){auto& row=task->program.equations[i];row.reserve(equations[i].size());
              for(std::size_t j=0;j<equations[i].size();++j)row.push_back({equations[i][j].column,source->coefficients[i][j].evaluate(values,task->field),0});}
            task->program.target_columns=task->target_columns;task->solver.emplace(task->program,task->field,solve);task->solver->solve(0);
            std::vector<std::vector<::ibp::LinearEntry>>().swap(task->program.equations);
            task->rows=task->solver->targets();task->selection=task->solver->selection();solved=true;
          }catch(...){}
          task->seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();flint_cleanup();
          {const std::lock_guard lock(task->mutex);task->solved=solved;task->finished=true;}task->ready.notify_all();
        }),task);}catch(const std::system_error&){break;}
      expected.task=std::move(task);++stats_.early_solves;++waiting;
    }
  }
  // The anticipated elimination at a requested point, if it used the present
  // equations and fill budget and fits the time the request has left. It
  // started before the request: if the request's own budget passes before it
  // finishes, that budget is spent.
  std::shared_ptr<Elimination> adopt(const Point& at,std::size_t max_fill,double seconds) {
    auto it=std::find_if(expected_.begin(),expected_.end(),[&](const Expected& expected){return expected.at==at;});
    if(it==expected_.end())return nullptr;auto task=std::move(it->task);expected_.erase(it);
    if(!task||task->system!=plan_->system||task->max_fill!=max_fill)return nullptr;
    std::unique_lock lock(task->mutex);
    if(!task->ready.wait_for(lock,std::chrono::duration<double>(seconds),[&]{return task->finished;}))throw std::runtime_error("elimination time budget exceeded");
    return task->solved&&task->seconds<=seconds?task:nullptr;
  }
  ::ibp::InputGeometry input_geometry_;std::vector<Exact> inputs_;std::unique_ptr<Plan> plan_;
  ibp::PropagatorBasis basis_;Exact dimension_;Options options_;Statistics stats_;
  std::vector<std::vector<ibp::Affine>> contractions_;std::vector<bool> trace_,zero_sectors_;
  unsigned threads_=1;std::uint64_t state_=0;std::size_t systems_=0;
  std::vector<Ahead> ahead_;modular::Word ahead_modulus_=0;std::vector<Expected> expected_;
  // Declared last: joined before anything else of the sampler goes away.
  std::vector<std::pair<std::jthread,std::shared_ptr<Elimination>>> tasks_;
};
struct ClosurePlanningStatistics {std::size_t passes=0,probes=0,demands=0,cache_hits=0;double seconds=0;};
class Session {
 public:
  Session(const ibp::PropagatorBasis& basis,const Exact& dimension,const ExactField& field,const Options& options)
    :sampler_(std::make_shared<Sampler>(basis,dimension,options)),session_(basis,dimension,field,configuration(options,sampler_)),basis_(basis),dimension_(dimension),field_(field),options_(options){}
  fire::Result operator()(const std::vector<ibp::Integral>& requested,const fire::Options& limits){
    if(!options_.closure_parameter)return session_(requested,limits);
    const auto started=std::chrono::steady_clock::now();
    try{auto planned=closure_plan(requested,limits);auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
      if(elapsed>=limits.timeout_seconds)throw std::runtime_error("IBP closure planning exhausted reconstruction budget");
      auto receiving=limits;receiving.timeout_seconds=std::max(1u,static_cast<unsigned>(limits.timeout_seconds-elapsed));return session_(planned,receiving);
    }catch(const std::exception& e){fire::Result result;result.reason=e.what();return result;}
  }
  const ClosurePlanningStatistics& planning_statistics()const{return planning_;}
  const Statistics& statistics()const{return sampler_->statistics();}
  // The number of master integrals of the level the requested integrals
  // generate: its masters and their derivatives in the closure parameter are
  // added to the demands until the set of masters is stable, each pass one
  // finite-field elimination at a fixed point, with no reconstruction. It
  // ranks alternative families (recursion::Options::count_masters) and is
  // kept on disk, as the ranking is repeated by every run on the same cache.
  // Counting only what the requests reduce to is not enough: a pair of
  // propagators can look small that way and close on a far larger level.
  std::size_t count_masters(const std::vector<ibp::Integral>& requested,const fire::Options& limits) {
    if(!options_.closure_parameter)throw std::invalid_argument("IBP master count needs a closure parameter");
    const auto parameter=*options_.closure_parameter;
    auto science=fire_batch::identity(basis_,dimension_,field_,requested,limits).json_value();
    boost::json::object identity{{"science",science},{"parameter",parameter},{"counter","one-point-derivative-closure-v1"},{"dots",options_.dots},{"numerators",options_.numerators},{"dots_first",options_.dots_first},{"target_sector_seeds",options_.target_sector_seeds}};
    const auto key=artifacts::detail::sha256(artifacts::detail::canonical(identity));
    const auto path=options_.cache_directory/"master-counts"/(key+".json");
    if(std::filesystem::exists(path)){
      auto cached=fire_modular::detail::read(path);
      if(artifacts::detail::string(cached.at("identity"))!=key||cached.at("schema")!="DiffExp.IBPMasterCount/v1")throw std::runtime_error("IBP master count identity mismatch");
      return cached.at("masters").to_number<std::size_t>();
    }
    const auto started=std::chrono::steady_clock::now();
    ibp::Generator generator(basis_,dimension_);
    std::set<ibp::Integral,ibp::IntegralOrder> demanded(requested.begin(),requested.end());
    std::size_t masters=0;bool stable=false;
    for(unsigned pass=0;pass<options_.max_planning_passes&&!stable;++pass){
      const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
      if(elapsed>=limits.timeout_seconds||demanded.size()>2000)throw std::runtime_error("IBP master count budget");
      auto receiving=limits;receiving.timeout_seconds=std::min(options_.probe_timeout_seconds,std::max(1u,static_cast<unsigned>(limits.timeout_seconds-elapsed)));
      auto reduced=(*sampler_)({demanded.begin(),demanded.end()},modular::point(dimension_.variable_count(),1,fire_modular::ordinals::planning_base+2*pass),prime(1),receiving);
      if(!reduced.success)throw std::runtime_error("IBP master count failed: "+reduced.reason);
      std::set<ibp::Integral> found;for(const auto& [a,row]:reduced.reductions)for(const auto& [master,c]:row)if(!c.is_zero())found.insert(master);
      const auto before=demanded.size();
      for(const auto& master:found){demanded.insert(master);for(const auto& [a,c]:generator.derivative(master,parameter))if(!c.is_zero())demanded.insert(a);}
      masters=found.size();stable=demanded.size()==before;
    }
    if(!stable)throw std::runtime_error("IBP master count did not stabilize");
    // The stable demands are what closure planning would find: leave them as
    // its hint, so that preparing the chosen level does not plan again.
    {
      boost::json::object plan_identity{{"science",science},{"parameter",parameter},{"planner","two-point-terminal-derivative-closure-v2"},{"dots",options_.dots},{"numerators",options_.numerators},{"dots_first",options_.dots_first},{"target_sector_seeds",options_.target_sector_seeds}};
      const auto plan_key=artifacts::detail::sha256(artifacts::detail::canonical(plan_identity));
      const auto plan_path=options_.cache_directory/"closure-plans"/(plan_key+".json");
      if(!std::filesystem::exists(plan_path)){
        boost::json::array rows;for(const auto& a:demanded){boost::json::array row;for(int power:a)row.push_back(power);rows.push_back(row);}
        std::filesystem::create_directories(plan_path.parent_path());
        fire_modular::detail::save(plan_path,{{"schema","DiffExp.IBPClosurePlan/v1"},{"identity",plan_key},{"planning_hint_only",true},{"demands",rows}});
      }
    }
    std::filesystem::create_directories(path.parent_path());
    fire_modular::detail::save(path,{{"schema","DiffExp.IBPMasterCount/v1"},{"identity",key},{"masters",masters}});
    return masters;
  }
 private:
  std::vector<ibp::Integral> closure_plan(const std::vector<ibp::Integral>& requested,const fire::Options& limits) {
    using namespace std::chrono;const auto started=steady_clock::now();
    struct Timing{ClosurePlanningStatistics& stats;steady_clock::time_point start;~Timing(){stats.seconds+=duration<double>(steady_clock::now()-start).count();}}timing{planning_,started};
    const auto parameter=*options_.closure_parameter;
    if(parameter>=dimension_.variable_count()||!options_.max_planning_passes||options_.max_planning_passes>16)throw std::invalid_argument("invalid finite-field closure planning bounds");
    auto science=fire_batch::identity(basis_,dimension_,field_,requested,limits).json_value();
    boost::json::object identity{{"science",science},{"parameter",parameter},{"planner","two-point-terminal-derivative-closure-v2"},{"dots",options_.dots},{"numerators",options_.numerators},{"dots_first",options_.dots_first},{"target_sector_seeds",options_.target_sector_seeds}};
    const auto key=artifacts::detail::sha256(artifacts::detail::canonical(identity));
    const auto path=options_.cache_directory/"closure-plans"/(key+".json");
    std::set<ibp::Integral,ibp::IntegralOrder> demanded(requested.begin(),requested.end());
    if(std::filesystem::exists(path)){
      auto cached=fire_modular::detail::read(path);
      if(artifacts::detail::string(cached.at("identity"))!=key||cached.at("schema")!="DiffExp.IBPClosurePlan/v1"||!cached.at("planning_hint_only").as_bool())throw std::runtime_error("IBP closure plan identity mismatch");
      std::set<ibp::Integral,ibp::IntegralOrder> planned;
      if(cached.at("demands").as_array().size()>2000)throw std::length_error("IBP closure plan demand budget");
      for(const auto& entry:cached.at("demands").as_array()){
        ibp::Integral a;for(const auto& power:entry.as_array())a.push_back(power.to_number<int>());fire_batch::validate_index(a,basis_);
        if(!planned.insert(a).second)throw std::runtime_error("duplicate closure plan demand");
      }
      for(const auto& a:requested)if(!planned.contains(a))throw std::runtime_error("closure plan omitted an original request");
      ++planning_.cache_hits;planning_.demands=planned.size();if(options_.progress)options_.progress("reused derivative-closure planning hint; exact reconstruction and closure checks remain required");
      return {planned.begin(),planned.end()};
    }
    ibp::Generator generator(basis_,dimension_);
    // Every pass begins with a complete elimination at a point fixed in
    // advance and, while the demands leave the seed bounds alone, of the same
    // equations; so does the validation of the reconstruction that follows.
    // A pass therefore announces the eliminations of the pass after it, of
    // that validation and of one pass more, in this order of likelihood: a
    // plan that settles in three passes waits for one elimination, not four.
    // Without selected rows the check at the second point is complete too.
    // Once planning ends only the validation's remain expected.
    std::vector<Sampler::Point> held;
    if(sampler_->works_ahead())try{held=complete(options_.selected_source_rows,options_.complete_validation_anchor,session_.held_out(limits,2));}catch(const std::exception&){}
    auto left=[&]{return limits.timeout_seconds-duration<double>(steady_clock::now()-started).count();};
    struct Forget{Sampler& sampler;const std::vector<Sampler::Point>& held;double seconds;~Forget(){sampler.anticipate(held,seconds);}}forget{*sampler_,held,double(limits.timeout_seconds)};
    auto expect=[&](unsigned pass){std::vector<Sampler::Point> points;const auto next=std::max(pass,1u);if(!sampler_->works_ahead())return;
      for(auto p=pass;p<=next+1&&p<options_.max_planning_passes;++p){
        for(unsigned check=0;check<(options_.selected_source_rows?1u:2u);++check)points.push_back({prime(1),modular::point(dimension_.variable_count(),1,fire_modular::ordinals::planning_base+2*p+check)});
        if(p==next)points.insert(points.end(),held.begin(),held.end());
      }
      sampler_->anticipate(std::move(points),left());};
    for(unsigned pass=0;pass<options_.max_planning_passes;++pass){
      expect(pass);
      if(demanded.size()>2000)throw std::length_error("IBP closure planning demand budget");
      std::vector<ibp::Integral> batch(demanded.begin(),demanded.end());std::set<ibp::Integral> terminal;
      for(unsigned check=0;check<2;++check){
        auto elapsed=duration<double>(steady_clock::now()-started).count();if(elapsed>=limits.timeout_seconds)throw std::runtime_error("IBP closure planning time budget");
        auto receiving=limits;receiving.timeout_seconds=std::min(options_.probe_timeout_seconds,std::max(1u,static_cast<unsigned>(limits.timeout_seconds-elapsed)));
        auto point=modular::point(dimension_.variable_count(),1,fire_modular::ordinals::planning_base+2*pass+check);
        auto reduced=(*sampler_)(batch,point,prime(1),receiving);++planning_.probes;
        if(!reduced.success)throw std::runtime_error("IBP closure planning probe failed: "+reduced.reason);
        std::set<ibp::Integral> found;for(const auto& [a,row]:reduced.reductions)for(const auto& [master,c]:row)if(!c.is_zero())found.insert(master);
        if(check&&found!=terminal)throw std::runtime_error("IBP closure planning basis changed at a second point");terminal=std::move(found);
      }
      if(terminal.size()>256)throw std::length_error("IBP closure planning master budget");
      ++planning_.passes;const auto before=demanded.size();
      for(const auto& master:terminal){demanded.insert(master);for(const auto& [a,c]:generator.derivative(master,parameter))if(!c.is_zero())demanded.insert(a);}
      planning_.demands=demanded.size();
      if(options_.progress)options_.progress("finite-field closure plan "+std::to_string(pass+1)+": "+std::to_string(terminal.size())+" candidate masters, "+std::to_string(demanded.size())+" combined demands");
      if(demanded.size()==before){
        boost::json::array rows;for(const auto& a:demanded){boost::json::array row;for(int power:a)row.push_back(power);rows.push_back(row);}
        std::filesystem::create_directories(path.parent_path());
        fire_modular::detail::save(path,{{"schema","DiffExp.IBPClosurePlan/v1"},{"identity",key},{"planning_hint_only",true},{"demands",rows}});
        return {demanded.begin(),demanded.end()};
      }
    }
    // The plan is a hint. One that has not stabilized is still a superset of
    // the request made of valid demands; it is not saved, and the level
    // preparer's own closure loop continues from its reconstruction.
    if(options_.progress)options_.progress("finite-field closure plan did not stabilize within "+std::to_string(options_.max_planning_passes)+" passes; reconstructing the "+std::to_string(demanded.size())+" demands found so far");
    return {demanded.begin(),demanded.end()};
  }
  // Of the held-out points, those eliminated completely can start early.
  static std::vector<Sampler::Point> complete(bool selected_rows,bool anchored,const std::vector<fire_modular::Options::HeldOut>& points) {
    std::vector<Sampler::Point> out;
    for(const auto& held:points)if(!selected_rows||(anchored&&held.anchor))out.push_back({prime(held.prime_index),held.point});
    return out;
  }
  static fire_modular::Options configuration(const Options& o,const std::shared_ptr<Sampler>& sampler){
    fire_modular::Options out;out.cache_directory=o.cache_directory;out.max_degree=o.max_degree;out.max_variable_degree=o.max_variable_degree;out.max_primes=o.max_primes;
    out.max_samples_per_prime=o.max_samples_per_prime;out.probe_timeout_seconds=o.probe_timeout_seconds;out.progress=o.progress;out.modulus=prime;
    out.sample_checkpoint_interval=o.sample_checkpoint_interval;
    out.provider_identity=std::string("ibp-solver-parametric-")+(o.dots_first?"dots-first":"total-degree")+"-v1-dots-"+std::to_string(o.dots)+"-numerators-"+std::to_string(o.numerators);
    if(o.target_sector_seeds)out.provider_identity+="-subsector-seeds-v1";
    if(o.selected_source_rows)out.provider_identity+="-fresh-original-row-subset-validation-v1";
    const bool anchored=o.selected_source_rows&&o.complete_validation_anchor;
    if(anchored)out.provider_identity+="-complete-anchor-v1";
    out.prime_table="descending-61-bit-primes-v1";out.prime_count=prime_count;out.sparse_lifting=true;out.degree_slices=o.degree_slices;
    out.sample_provider=[sampler](const auto& r,const auto& p,auto m,const auto& limits){return (*sampler)(r,p,m,limits);};
    out.validation_provider=[sampler](const auto& r,const auto& p,auto m,const auto& limits){return (*sampler)(r,p,m,limits,true);};
    if(anchored)out.anchor_provider=[sampler](const auto& r,const auto& p,auto m,const auto& limits){return (*sampler)(r,p,m,limits,true,true);};
    if(sampler->works_ahead()){
      out.upcoming_samples=[sampler](const auto& points,auto m){sampler->announce(points,m);};
      out.upcoming_validation=[sampler,selected=o.selected_source_rows,anchor=o.complete_validation_anchor](const auto& points,double seconds){sampler->anticipate(complete(selected,anchor,points),seconds);};
    }
    return out;
  }
  std::shared_ptr<Sampler> sampler_;fire_modular::Session session_;
  ibp::PropagatorBasis basis_;Exact dimension_;ExactField field_;Options options_;ClosurePlanningStatistics planning_;
};
// The settings every recursive Feynman-trick level uses with the built-in
// solver. Coefficients are typically of low degree in epsilon and high degree
// in the merge parameter: per-variable degree slices fit that shape with far
// fewer points than the dense total-degree scan, whose cap they also avoid.
// The closure is planned once per level, not once per pass. Unless `strict`,
// held-out validation eliminates the pilot's dependency subset of the original
// equations with one complete elimination per attempt; strict eliminates every
// generated equation at every point.
inline Options recursion_options(const std::filesystem::path& cache_directory,std::size_t parameter,
    unsigned elimination_seconds,bool strict=false) {
  Options native;native.cache_directory=cache_directory;
  native.degree_slices=true;native.max_samples_per_prime=1024;
  native.closure_parameter=parameter;native.probe_timeout_seconds=elimination_seconds;
  if(!strict){native.selected_source_rows=true;native.complete_validation_anchor=true;}
  return native;
}
} // namespace diffexp::ibp_solver
