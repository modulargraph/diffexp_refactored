#pragma once
#include "diffexp/fire_modular.hpp"
#include "diffexp/scaleless_sectors.hpp"
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
  // The seeds every sector starts from. They are the least a reduction is
  // tried with, not what it is confined to: a requested integral asks for the
  // seeds that generate it and an integral left unreduced at the edge of the
  // seeds widens them (see target_seeds and max_seed_widenings).
  unsigned dots=1,numerators=1,max_degree=16,max_variable_degree=64,max_primes=32,probe_timeout_seconds=20;
  std::size_t max_samples_per_prime=512;
  // Ordinary samples are recomputed from their deterministic points in a few
  // milliseconds, so only every n-th is made durable; held-out validation
  // samples, which the completion record replays, are always retained.
  std::size_t sample_checkpoint_interval=32;
  bool degree_slices=false,dots_first=false;
  // The sectors that receive, in advance, the seeds that generate a requested
  // integral (as many numerators, one dot fewer): `none`; `sector`, its own;
  // `subsectors`, its own and every subsector; `all`, every sector. With
  // `none` a requested integral is treated as any other: the smaller seeds
  // often reduce it already, and where they do not it is left unreduced
  // outside the seeds of its sector, which widens that sector alone.
  enum class TargetSeeds {none,sector,subsectors,all};
  TargetSeeds target_seeds=TargetSeeds::none;
  // Sectors proved scaleless from the Symanzik polynomials are zero sectors
  // (scaleless_sectors.hpp); otherwise only those with a free loop momentum.
  bool scaleless_sectors=true;
  // An integral left unreduced that lies outside the seeds of its sector (more
  // dots or a higher numerator degree than any seed there) was never the seed
  // of an equation: it is not a master integral but the edge of the system.
  // The seeds of its sector alone are then widened to contain it and the
  // elimination is continued with the new equations, at most this many times
  // for one batch: what they leave unreduced outside the seeds of a subsector
  // widens that subsector in the next repetition.
  unsigned max_seed_widenings=8;
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
  // How the arithmetic trace of the selected equations eliminates them
  // (::ibp::SolveOptions::complete_rows): reducing every row completely
  // shares, between the targets of a level, what the leading way repeats for
  // each of them.
  bool trace_complete_rows=false;
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
  // Eliminations of the selected equations served from a batch.
  std::size_t batched_solves=0;
  // Eliminations repeated with wider seeds (Options::max_seed_widenings).
  std::size_t seed_widenings=0;
  // Complete eliminations that did not start from nothing (they are not among
  // full_solves): continued from the pivots of an earlier one at the same
  // point, with the equations it had not seen, and those needing no new
  // equation at all. continued_equations counts the equations inserted.
  std::size_t continued_solves=0,continued_equations=0,reused_solves=0;double continued_solve_seconds=0;
  // The largest system eliminated: its equations and integrals, and the most
  // entries the pivot rows of one elimination held.
  std::size_t largest_system=0,largest_integrals=0,peak_fill=0;
  // The longest arithmetic trace replayed and the most equations it came from.
  std::size_t trace_instructions=0,trace_rows=0;
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
    if(!basis_.space.loops||!basis_.physical_count||basis_.denominators.size()!=basis_.space.size()||basis_.physical_count>basis_.denominators.size())
      throw std::invalid_argument("IBP solver needs a complete basis of scalar products with at least one loop and one physical denominator");
    // What is left of a size limit is the storage of one integral (a power
    // per scalar product) and the dense table of sectors; neither bounds the
    // number of loops or legs as such.
    if(basis_.space.size()>::ibp::Integral::slots)
      throw std::invalid_argument("IBP solver: "+std::to_string(basis_.space.size())+" scalar products exceed the "+std::to_string(::ibp::Integral::slots)+" powers an integral stores; rebuild with a larger IBP_SOLVER_INTEGRAL_SLOTS");
    if(basis_.physical_count>::ibp::Integral::max_physical)
      throw std::invalid_argument("IBP solver: the table of the 2^"+std::to_string(basis_.physical_count)+" sectors does not fit");
    if(options_.dots>::ibp::Integral::max_seed_power||options_.numerators>::ibp::Integral::max_seed_power)throw std::invalid_argument("IBP solver seed powers exceed what an integral stores");
    if(threads_>256)throw std::invalid_argument("IBP solver threads must be 0..256");
    contractions_=ibp::contractions(basis_);
    for(unsigned i=0;i<basis_.space.loops;++i)for(unsigned v=0;v<basis_.space.loops+basis_.space.externals();++v)trace_.push_back(i==v);
    zero_sectors_.resize(std::size_t(1)<<basis_.physical_count);
    for(const auto& a:fire::free_loop_sectors(basis_))zero_sectors_[sector(a)]=true;
    if(options_.scaleless_sectors){const auto scaleless=scaleless::sectors(basis_);for(std::size_t mask=0;mask<scaleless.size();++mask)if(scaleless[mask])zero_sectors_[mask]=true;}
    input_geometry_.loops=basis_.space.loops;input_geometry_.externals=basis_.space.externals();
    input_geometry_.physical=basis_.physical_count;input_geometry_.n=basis_.space.size();
    input_geometry_.trace=trace_;input_geometry_.zero_sectors=zero_sectors_;
    // Equal coefficients share one input slot.
    auto input=[&](const Exact& value){if(value.is_zero())return ::ibp::InputGeometry::zero;
      auto it=std::find(inputs_.begin(),inputs_.end(),value);if(it==inputs_.end()){inputs_.push_back(value);input_geometry_.constant_inputs.push_back(value.is_rational());return std::uint32_t(inputs_.size()-1);}
      return std::uint32_t(it-inputs_.begin());};
    // Even a constant zero dimension needs a valid explicit input slot.
    inputs_.push_back(dimension_);input_geometry_.constant_inputs.push_back(dimension_.is_rational());input_geometry_.dimension=0;
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
  // The seed bounds widened for the present batch, as (sector, dots,
  // numerators), and their restoration for a batch whose closure plan recorded
  // them: a later run then reduces the batch with the equations the plan was
  // made with, instead of finding the widening again.
  std::vector<std::array<unsigned,3>> widened()const {
    std::vector<std::array<unsigned,3>> out;
    for(std::size_t mask=0;mask<widened_.size();++mask)for(const auto& corner:widened_[mask])out.push_back({unsigned(mask),corner.first,corner.second});
    return out;
  }
  void widen(const std::vector<ibp::Integral>& batch,const std::vector<std::array<unsigned,3>>& bounds) {
    std::vector<::ibp::Integral> targets;
    for(const auto& a:batch){fire_batch::validate_index(a,basis_);::ibp::Integral target;for(unsigned i=0;i<a.size();++i)target.powers[i]=a[i];targets.push_back(target);}
    std::vector<Corners> wider;if(!bounds.empty())wider.resize(std::size_t(1)<<basis_.physical_count);
    for(const auto& [mask,dots,nums]:bounds){if(mask>=wider.size()||dots>::ibp::Integral::max_seed_power||nums>::ibp::Integral::max_seed_power)throw std::invalid_argument("IBP solver widened seed bounds");wider[mask].push_back({dots,nums});}
    widened_=std::move(wider);widened_for_=targets;std::sort(targets.begin(),targets.end());widened_sorted_=std::move(targets);widenings_=0;
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
      if(!limits.timeout_seconds||!limits.memory_bytes||requested.empty())
        throw std::invalid_argument("IBP solver requires a time and a memory budget and a nonempty request");
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
        if(options_.target_seeds==Options::TargetSeeds::all){seeds.dots=std::max(seeds.dots,dots?dots-1:0);seeds.numerators=std::max(seeds.numerators,nums);}
        targets.push_back(target);
      }
      std::vector<std::pair<unsigned,unsigned>> sector_bounds;
      if(options_.target_seeds!=Options::TargetSeeds::all){
        sector_bounds.resize(std::size_t(1)<<basis_.physical_count,{options_.dots,options_.numerators});
        if(options_.target_seeds!=Options::TargetSeeds::none)for(const auto& target:requested){unsigned dots=0,nums=0;for(int power:target){if(power>0)dots+=power-1;else nums-=power;}
          const auto mask=sector(target);for(std::size_t sub=mask;;sub=(sub-1)&mask){
            auto& bound=sector_bounds[sub];bound.first=std::max(bound.first,dots?dots-1:0);bound.second=std::max(bound.second,nums);
            if(!sub||options_.target_seeds==Options::TargetSeeds::sector)break;}
        }
      }
      // Widened seeds belong to the batch they were found for and to every
      // batch that contains it: the batches of a closure plan only grow, and
      // so do their equations, which lets an elimination continue from the
      // one before it. Any other batch starts from its own seeds.
      if(targets!=widened_for_){
        auto sorted=targets;std::sort(sorted.begin(),sorted.end());
        if(!std::includes(sorted.begin(),sorted.end(),widened_sorted_.begin(),widened_sorted_.end()))widened_.clear();
        widened_for_=targets;widened_sorted_=std::move(sorted);widenings_=0;}
      if(!widened_.empty()&&sector_bounds.empty())sector_bounds.assign(std::size_t(1)<<basis_.physical_count,{seeds.dots,seeds.numerators});
      // A sector's seeds are those of its starting box and of every corner
      // it was widened to; sector_bounds becomes the box around them all.
      const auto base_bounds=sector_bounds;const auto base_seeds=std::pair{seeds.dots,seeds.numerators};
      const auto base_of=[&](std::size_t mask){return base_bounds.empty()?base_seeds:base_bounds[mask];};
      for(std::size_t mask=0;mask<sector_bounds.size();++mask){auto& bound=sector_bounds[mask];
        if(!widened_.empty())for(const auto& corner:widened_[mask]){bound.first=std::max(bound.first,corner.first);bound.second=std::max(bound.second,corner.second);}
        seeds.dots=std::max(seeds.dots,bound.first);seeds.numerators=std::max(seeds.numerators,bound.second);}
      std::function<bool(std::size_t,unsigned,unsigned)> admitted;
      if(!widened_.empty())admitted=[&](std::size_t mask,unsigned dots,unsigned nums){
        const auto base=base_bounds[mask];if(dots<=base.first&&nums<=base.second)return true;
        for(const auto& corner:widened_[mask])if(dots<=corner.first&&nums<=corner.second)return true;
        return false;};
      seeds.max_terms=limits.memory_bytes/64;seeds.max_seeds=std::size_t(1)<<30;seeds.max_states=std::size_t(1)<<34;seeds.seconds=remaining();
      if(!seeds.max_terms)throw std::invalid_argument("IBP solver memory budget too small");
      if(!plan_||plan_->program->program.targets!=targets||plan_->zero_sectors!=geometry.zero_sectors||plan_->dots!=seeds.dots||plan_->numerators!=seeds.numerators||plan_->max_terms!=seeds.max_terms||plan_->bounds!=sector_bounds||plan_->corners!=widened_){
        auto fresh=std::make_unique<Plan>();fresh->program=std::make_shared<::ibp::ParametricProgram>(::ibp::ParametricProgram::compile(geometry,seeds,targets,options_.dots_first?::ibp::ParametricOrdering::dots_first: ::ibp::ParametricOrdering::total_degree,sector_bounds,admitted));
        fresh->zero_sectors=geometry.zero_sectors;fresh->dots=seeds.dots;fresh->numerators=seeds.numerators;fresh->max_terms=seeds.max_terms;fresh->bounds=sector_bounds;fresh->corners=widened_;
        // New targets usually leave the generated equations as they were, and
        // with them every complete elimination started ahead.
        fresh->system=plan_&&same_equations(*plan_->program,*fresh->program)?plan_->system:++systems_;
        for(const auto& equation:fresh->program->program.equations)fresh->terms+=equation.size();
        plan_=std::move(fresh);++state_;++stats_.templates;stats_.equations+=plan_->program->program.equations.size();
      }
      auto& parametric=*plan_->program;auto& program=parametric.program;
      const auto generated=std::chrono::steady_clock::now();
      // The pivot rows of one elimination may take a third of the memory
      // budget (an entry is 16 bytes); nothing else bounds them. Eliminations
      // kept for later and those started ahead share what is left (room).
      ::ibp::SolveOptions solve;solve.seconds=remaining();solve.max_fill=limits.memory_bytes/48;memory_=limits.memory_bytes;
      stats_.largest_system=std::max(stats_.largest_system,program.equations.size());stats_.largest_integrals=std::max(stats_.largest_integrals,program.integrals.size());
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
      // Without a trace (it did not fit its budget at this prime) the
      // announced points are eliminated together as well: each the selected
      // equations alone, as this thread would have done one after another.
      if(!independent&&learned.disabled&&options_.selected_source_rows&&plan_->selection&&threads_>1)if(auto* slot=announced(point,modulus)){
        if(!slot->result||slot->state!=state_)eliminate_announced(requested,modulus,solve,limits);
        if(slot->result&&slot->state==state_){
          out=std::move(*slot->result);++stats_.selected_solves;++stats_.batched_solves;stats_.selected_equations+=plan_->selection->size();stats_.selected_solve_seconds+=slot->seconds;
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
        std::shared_ptr<Owned> owned;std::shared_ptr<Elimination> early;std::shared_ptr<void> held;::ibp::Solver* solver=nullptr;
        std::vector<::ibp::Row> full;std::optional<std::vector<std::uint32_t>> dependencies;
        auto terminals=[](const auto& reduced){std::set<::ibp::Column> out;for(const auto& row:reduced)for(const auto& entry:row)out.insert(entry.column);return out;};
        // A point is revisited by a wider batch or wider seeds; a held-out
        // point is not, once the seeds of its batch are known to suffice.
        const bool revisited=!independent||!plan_->contained;
        double primal_seconds=0;
        auto primal=[&](bool selected){
          auto started=std::chrono::steady_clock::now();double seconds=0;
          const auto* sources=selected?&*plan_->selection:nullptr;dependencies.reset();early.reset();owned.reset();
          // A complete elimination of these equations at this point is kept:
          // only the targets are reduced.
          Kept* kept=selected?nullptr:find_kept({modulus,point},solve.max_fill);
          if(kept&&kept->system==plan_->system){
            held=kept->owner;solver=kept->solver;solver->restart(remaining());full=solver->targets(program.target_columns);
            ++stats_.reused_solves;stats_.continued_solve_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();return;
          }
          if(!selected&&(early=adopt({modulus,point},solve.max_fill,remaining()))){
            // Its seconds count once: not those this thread spent waiting for it.
            elsewhere+=std::max(0.0,early->seconds-std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count());
            started=std::chrono::steady_clock::now();
            if(early->target_columns==program.target_columns){full=std::move(early->rows);dependencies=std::move(early->selection);}
            // Other targets are reduced here. The solver's clock has run since
            // its elimination began on the other thread: if it expires only
            // now, the elimination is repeated with a budget of its own.
            else try{early->program.target_columns=program.target_columns;full=early->solver->targets();}
              catch(const std::runtime_error&){early.reset();}
            if(early){solver=&*early->solver;seconds=early->seconds;++stats_.early_solves_used;noted(*solver);
              if(revisited)keep({{modulus,point},plan_->system,solve.max_fill,plan_->zero_sectors,plan_->program,solver,early});}
          }
          std::optional<Continued> continued;
          if(!early){
            // A kept elimination of fewer equations at this point continues
            // with the equations it has not seen.
            started=std::chrono::steady_clock::now();
            if(kept)continued=continue_from(*kept,parametric);
            if(kept)drop_kept({modulus,point});
            if(continued)sources=&continued->rows;
            // Rebind ALL original rows on fallback, including those left at an
            // earlier point or prime by the preceding subset solve.
            parametric.bind(values,field,sources);solve.seconds=remaining();
            // A long elimination says how far it is: equations inserted of
            // those it was given, and the entries its pivot rows hold.
            // All equations: reduced completely. A few selected ones: by
            // their leading columns alone (::ibp::SolveOptions).
            auto telling=solve;telling.complete_rows=!selected;
            if(options_.progress)telling.progress=[this,total=sources?sources->size():program.equations.size()](const ::ibp::Statistics& s,double seconds){
              options_.progress("elimination: "+std::to_string(s.equations)+" of "+std::to_string(total)+" equations, rank "+std::to_string(s.rank)+", pivot fill "+std::to_string(s.fill)+", longest row "+std::to_string(s.peak_row)+", "+std::to_string(static_cast<long>(seconds))+" s");};
            owned=std::make_shared<Owned>(modulus,plan_->program,telling);solver=&owned->solver;
            solver->solve(0,sources,continued?std::move(continued->pivots):std::vector<::ibp::Pivot>{});full=solver->targets();
            if(!selected)noted(*solver);
            if(!selected&&revisited)keep({{modulus,point},plan_->system,solve.max_fill,plan_->zero_sectors,plan_->program,solver,owned});
          }
          seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();primal_seconds=seconds;
          if(selected){++stats_.selected_solves;stats_.selected_equations+=sources->size();stats_.selected_solve_seconds+=seconds;}
          else if(continued){++stats_.continued_solves;stats_.continued_equations+=continued->rows.size();stats_.continued_solve_seconds+=seconds;}
          else {++stats_.full_solves;stats_.full_solve_seconds+=seconds;}
        };
        // After a complete elimination the selection reads the equations it
        // needs again, at this point: the program's own coefficients may
        // since have been bound to another one.
        auto selection=[&]()->const std::vector<std::uint32_t>&{
          if(!dependencies){std::optional<std::vector<::ibp::Word>> forms;
            dependencies=solver->selection(program.target_columns,[&](std::uint32_t source){
              if(!forms)forms=parametric.evaluate_forms(values,field);
              ::ibp::Row row;const auto& equation=program.equations.at(source);
              for(std::size_t j=0;j<equation.size();++j)if(const auto value=(*forms)[parametric.coefficients[source][j]])row.push_back({equation[j].column,value});
              return row;});}
          return *dependencies;};
        bool selected=options_.selected_source_rows&&plan_->selection.has_value()&&!complete;primal(selected);
        // The first complete elimination of a batch shows whether its seeds
        // contain every integral it leaves unreduced.
        if(!selected&&!plan_->contained){
          // An integral left unreduced at the edge of seeds that were
          // widened (as many numerators, or as many dots, as the widest seed
          // of its sector) is no master either: seeds one step wider reduce
          // it. Taken for a master, its derivative would lie outside the
          // seeds again, and every planning pass would climb one step
          // further without end.
          // A widening adds a corner (dots, numerators) to its sector: the
          // seeds of the sector are the union of the boxes below its corners,
          // not the box around them, which holds many times more seeds.
          std::vector<Corners> wider;std::string example;bool at_edge=false;
          if(widenings_<options_.max_seed_widenings)for(const auto column:terminals(full)){
            const auto& a=program.integrals[column];unsigned dots=0,nums=0;std::size_t mask=0;
            for(unsigned i=0;i<input_geometry_.n;++i){if(a.powers[i]>0){dots+=a.powers[i]-1;if(i<basis_.physical_count)mask|=std::size_t(1)<<i;}else nums-=a.powers[i];}
            const auto base=base_of(mask);static const Corners none;const Corners& corners=widened_.empty()?none:widened_[mask];
            // The corner the integral lies under, the highest in numerators.
            const std::pair<unsigned,unsigned>* cover=nullptr;
            for(const auto& corner:corners)if(dots<=corner.first&&nums<=corner.second&&(!cover||std::pair{corner.second,corner.first}>std::pair{cover->second,cover->first}))cover=&corner;
            const bool outside=!cover&&(dots>base.first||nums>base.second);
            bool numerator_edge=!outside&&!corners.empty()&&nums>options_.numerators,dot_edge=!outside&&!corners.empty()&&dots>options_.dots;
            for(const auto& corner:corners){if(dots<=corner.first&&nums<corner.second)numerator_edge=false;if(dots<corner.first&&nums<=corner.second)dot_edge=false;}
            const bool edge=numerator_edge||dot_edge;
            if((!outside&&!edge)||dots>=::ibp::Integral::max_seed_power||nums>=::ibp::Integral::max_seed_power)continue;
            if(wider.empty()){wider=widened_;wider.resize(std::size_t(1)<<basis_.physical_count);}
            if(example.empty()){at_edge=edge;example="[";for(unsigned i=0;i<input_geometry_.n;++i)example+=std::to_string(a.powers[i])+",";example.back()=']';}
            const std::pair<unsigned,unsigned> corner=outside?std::pair{std::max(dots,base.first),std::max(nums,base.second)}:std::pair{(cover?cover->first:dots)+(dot_edge?1u:0u),(cover?cover->second:nums)+(numerator_edge?1u:0u)};
            auto& w=wider[mask];
            if(std::none_of(w.begin(),w.end(),[&](const auto& c){return corner.first<=c.first&&corner.second<=c.second;})){
              std::erase_if(w,[&](const auto& c){return c.first<=corner.first&&c.second<=corner.second;});w.push_back(corner);}
          }
          if(!wider.empty()){
            widened_=std::move(wider);++widenings_;++stats_.seed_widenings;
            if(options_.progress)options_.progress("integral "+example+(at_edge?" is left unreduced at the edge of the widened seeds of its sector":" is left unreduced outside the seeds of its sector")+"; widening the seeds of "+std::to_string(std::count_if(widened_.begin(),widened_.end(),[](const auto& w){return !w.empty();}))+" sectors");
            auto again=limits;again.timeout_seconds=std::max(1u,static_cast<unsigned>(remaining()));
            return (*this)(requested,point,modulus,again,independent,complete);
          }
          plan_->contained=true;
        }
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
          else if(!learned.disabled){const auto learning_begin=std::chrono::steady_clock::now();::ibp::ArithmeticTrace::Options tracing;tracing.max_fill=solve.max_fill;
            // A trace pays where it is learned in about the time of an
            // elimination: a replay then beats eliminating the selected
            // equations again. Where learning takes much longer (it also
            // takes far more memory), the selected equations are eliminated
            // at every point instead, which is as fast as replaying a long
            // trace and runs in batches as well (eliminate_announced).
            tracing.seconds=std::min(remaining(),std::max(0.25,2*primal_seconds));
            tracing.max_nodes=options_.max_trace_nodes?options_.max_trace_nodes:std::min<std::size_t>(std::numeric_limits<std::uint32_t>::max(),limits.memory_bytes/128);
            // The primal rows above are already a valid reduction. A trace that
            // does not fit its node, fill or time budget only costs the replay
            // shortcut: this prime continues with primal elimination alone.
            tracing.complete_rows=options_.trace_complete_rows;
            try{learned.trace=::ibp::ArithmeticTrace::learn_parametric(parametric,field,values,selection(),tracing);
              if(options_.progress&&!plan_->trace_told){plan_->trace_told=true;options_.progress("arithmetic trace: "+std::to_string(learned.trace->instructions())+" instructions, "+std::to_string(learned.trace->registers())+" registers, from "+std::to_string(selection().size())+" equations and "+std::to_string(program.target_columns.size())+" targets");}}
            catch(const std::length_error&){learned.disabled=true;}
            catch(const std::runtime_error&){learned.disabled=true;learning_overran=true;}
            if(learned.disabled){learned.trace.reset();++stats_.trace_disabled;
              if(options_.progress&&!plan_->trace_told){plan_->trace_told=true;options_.progress("no arithmetic trace (learning one takes longer than eliminating the "+std::to_string(selection().size())+" selected equations a few times, or it does not fit its memory); samples eliminate them at every point");}}
            else if(stats_.trace_instructions=std::max(stats_.trace_instructions,learned.trace->instructions()),stats_.trace_rows=std::max(stats_.trace_rows,selection().size());learned.trace->evaluate_inputs(values,field,learned.workspace)!=rows)throw std::runtime_error("IBP parametric trace disagrees with primal learning reduction");
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
    // The residues are handed on as machine words: the fit reads them as they are.
    auto& residues=out.residues.emplace();std::map<::ibp::Column,ibp::Integral> terminal;
    for(unsigned i=0;i<rows.size();++i){auto& row=residues[requested[i]];if(!row.empty())continue;
      for(const auto& term:rows[i])if(term.value)row.emplace(terminal.try_emplace(term.column,indices(program.integrals[term.column])).first->second,term.value);}
    for(const auto& [column,a]:terminal){auto [it,inserted]=residues.try_emplace(a);
      if(inserted)it->second.emplace(a,1);else if(it->second.size()!=1||it->second.begin()->first!=a||it->second.begin()->second!=1)throw std::runtime_error("IBP solver returned a nonterminal output basis");}
    out.success=true;fire_modular::detail::validate(residues,basis_,requested);
  }
  using Corners=std::vector<std::pair<unsigned,unsigned>>;
  struct Learned {std::optional<::ibp::ArithmeticTrace> trace;::ibp::ArithmeticTrace::Workspace workspace;bool checked=false,disabled=false;};
  // `system` names the generated equations: plans that differ in targets alone share it.
  struct Plan {std::shared_ptr<::ibp::ParametricProgram> program;std::vector<bool> zero_sectors;unsigned dots=0,numerators=0;std::size_t max_terms=0,system=0,terms=0;std::vector<std::pair<unsigned,unsigned>> bounds;std::vector<Corners> corners;bool contained=false,trace_told=false;std::map<modular::Word,Learned> primes;modular::Word sampling_modulus=0;
    std::optional<std::vector<std::uint32_t>> selection;std::set<::ibp::Column> terminal_columns;};
  static bool same_equations(const ::ibp::ParametricProgram& a,const ::ibp::ParametricProgram& b) {
    if(a.input_count!=b.input_count||a.constant_inputs!=b.constant_inputs||a.program.integrals!=b.program.integrals||a.program.equations.size()!=b.program.equations.size())return false;
    for(std::size_t i=0;i<a.program.equations.size();++i){const auto& x=a.program.equations[i];const auto& y=b.program.equations[i];if(x.size()!=y.size())return false;
      for(std::size_t j=0;j<x.size();++j)if(x[j].column!=y[j].column||!(a.coefficient(i,j)==b.coefficient(i,j)))return false;}
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
  // Eliminates the selected equations at a bounded stretch of announced
  // points, the requested one first, each worker with its own equations and
  // solver. A point whose elimination fails, or leaves other integrals
  // unreduced than the pilot did, keeps no result: its own request then takes
  // the serial path, which falls back to all equations and reports failures.
  void eliminate_announced(const std::vector<ibp::Integral>& requested,modular::Word modulus,const ::ibp::SolveOptions& probe,const fire::Options& limits) {
    using clock=std::chrono::steady_clock;auto since=[](clock::time_point t){return std::chrono::duration<double>(clock::now()-t).count();};
    std::vector<Ahead*> batch;const std::size_t stretch=4*threads_;
    for(auto& slot:ahead_){if(batch.size()==stretch)break;if(!slot.result||slot.state!=state_){slot.result.reset();batch.push_back(&slot);}}
    const auto began=clock::now();std::atomic<std::size_t> next{0};
    const auto& parametric=*plan_->program;const auto& program=parametric.program;const auto& selected=*plan_->selection;const auto& terminal=plan_->terminal_columns;
    auto work=[&]{
      ::ibp::Field field(modulus);::ibp::Program local;std::optional<::ibp::Solver> solver;
      for(;;){const auto i=next.fetch_add(1);if(i>=batch.size()||since(began)>=limits.timeout_seconds)break;auto& slot=*batch[i];const auto start=clock::now();
        try{std::vector<::ibp::Word> values;for(const auto& input:inputs_)values.push_back(modular::evaluate(input,slot.point,modulus));
          const auto forms=parametric.evaluate_forms(values,field);
          if(!solver){local.integrals.resize(program.integrals.size());local.equations.resize(program.equations.size());local.target_columns=program.target_columns;
            for(auto source:selected){auto& row=local.equations[source];row.reserve(program.equations[source].size());for(const auto& entry:program.equations[source])row.push_back({entry.column,0,0});}
            auto solve=probe;solve.complete_rows=false;solve.progress={};solve.seconds=std::max(1.0,limits.timeout_seconds-since(began));solver.emplace(local,field,solve);}
          for(auto source:selected){auto& row=local.equations[source];for(std::size_t j=0;j<row.size();++j)row[j].constant=forms[parametric.coefficients[source][j]];}
          solver->restart(std::max(1.0,limits.timeout_seconds-since(began)));solver->solve(0,&selected);auto rows=solver->targets();
          std::set<::ibp::Column> found;for(const auto& row:rows)for(const auto& entry:row)found.insert(entry.column);
          if(found!=terminal)continue;
          fire::Result reduced;emit(reduced,rows,program,requested);slot.seconds=since(start);slot.state=state_;slot.result=std::move(reduced);
        }catch(...){}
      }};
    const auto busy=std::min<std::size_t>(running(),threads_-1);
    std::vector<std::jthread> pool;
    try{for(std::size_t i=1;i<std::min<std::size_t>(batch.size(),threads_-busy);++i)pool.emplace_back([&]{work();flint_cleanup();});}catch(const std::system_error&){}
    work();
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
  // A complete elimination this thread made, with what its solver refers to.
  struct Owned {
    ::ibp::Field field;std::shared_ptr<const ::ibp::ParametricProgram> source;::ibp::Solver solver;
    Owned(modular::Word modulus,std::shared_ptr<const ::ibp::ParametricProgram> program,const ::ibp::SolveOptions& options)
      :field(modulus),source(std::move(program)),solver(source->program,field,options){}
  };
  // The last complete eliminations, by point. The batches of a closure plan
  // are reduced at the same two points, each a superset of the one before:
  // the same equations then need no elimination, and more equations only the
  // elimination of the new ones (continue_from).
  struct Kept {Point at;std::size_t system=0,max_fill=0;std::vector<bool> zero_sectors;std::shared_ptr<const ::ibp::ParametricProgram> program;::ibp::Solver* solver=nullptr;std::shared_ptr<void> owner;};
  static constexpr std::size_t kept_points=2;
  Kept* find_kept(const Point& at,std::size_t max_fill) {
    for(auto& kept:kept_)if(kept.at==at&&kept.max_fill==max_fill&&kept.zero_sectors==plan_->zero_sectors)return &kept;
    return nullptr;
  }
  void drop_kept(const Point& at){std::erase_if(kept_,[&](const Kept& kept){return kept.at==at;});}
  void keep(Kept kept){drop_kept(kept.at);if(kept_.size()==kept_points)kept_.erase(kept_.begin());if(!room(kept_.size()+running()+1)&&!kept_.empty())kept_.erase(kept_.begin());kept_.push_back(std::move(kept));}
  // The fill of the last complete elimination predicts that of the next one
  // of the same system. `sets` eliminations, each with its pivot rows (16
  // bytes an entry) and, while it runs on another thread, its own copy of the
  // equations (24 bytes a term), must fit in half the budget.
  void noted(const ::ibp::Solver& solver){last_fill_=solver.statistics().fill;stats_.peak_fill=std::max(stats_.peak_fill,last_fill_);}
  bool room(std::size_t sets)const{return !last_fill_||!memory_||2*sets*(last_fill_*16+(plan_?plan_->terms*24:0))<=memory_;}
  // What continues a kept elimination under the present equations: its pivots
  // in the present columns and equation indices, and the equations it has not
  // seen. Nothing unless every equation it eliminated is, term for term, one
  // of the present equations and the columns keep their order: the pivots are
  // then an echelon basis of a subspace of the present row space, inserting
  // the remaining equations completes it, and the reduced targets, which are
  // normal forms modulo that row space, are those of an elimination from
  // nothing. The kept elimination gives up its pivots.
  struct Continued {std::vector<::ibp::Pivot> pivots;std::vector<std::uint32_t> rows;};
  struct OriginHash {std::size_t operator()(const ::ibp::ParametricProgram::Origin& o)const{return ::ibp::IntegralHash{}(o.seed)*1099511628211ULL^o.op;}};
  static std::optional<Continued> continue_from(Kept& kept,const ::ibp::ParametricProgram& now) {
    const auto& old=*kept.program;constexpr auto absent=std::numeric_limits<::ibp::Column>::max();
    if(old.input_count!=now.input_count||old.constant_inputs!=now.constant_inputs||old.origins.size()!=old.program.equations.size()||now.origins.size()!=now.program.equations.size())return std::nullopt;
    std::unordered_map<::ibp::Integral,::ibp::Column,::ibp::IntegralHash> ids;ids.reserve(now.program.integrals.size());
    for(::ibp::Column i=0;i<now.program.integrals.size();++i)ids.emplace(now.program.integrals[i],i);
    std::vector<::ibp::Column> column(old.program.integrals.size(),absent);::ibp::Column last=absent;
    for(std::size_t i=0;i<column.size();++i)if(auto it=ids.find(old.program.integrals[i]);it!=ids.end()){
      if(last!=absent&&it->second<=last)return std::nullopt;column[i]=last=it->second;}
    std::unordered_map<::ibp::ParametricProgram::Origin,std::uint32_t,OriginHash> origins;origins.reserve(now.origins.size());
    for(std::uint32_t i=0;i<now.origins.size();++i)if(!origins.emplace(now.origins[i],i).second)return std::nullopt;
    std::vector<std::uint32_t> equation(old.origins.size());std::vector<bool> seen(now.origins.size());
    for(std::size_t i=0;i<old.origins.size();++i){
      auto it=origins.find(old.origins[i]);if(it==origins.end())return std::nullopt;
      const auto j=it->second;const auto& a=old.program.equations[i];const auto& b=now.program.equations[j];if(a.size()!=b.size()||seen[j])return std::nullopt;
      for(std::size_t t=0;t<a.size();++t)if(column[a[t].column]!=b[t].column||!(old.coefficient(i,t)==now.coefficient(j,t)))return std::nullopt;
      equation[i]=j;seen[j]=true;
    }
    Continued out;for(std::uint32_t j=0;j<seen.size();++j)if(!seen[j])out.rows.push_back(j);
    out.pivots=kept.solver->release();
    for(auto& pivot:out.pivots){pivot.source=equation.at(pivot.source);for(auto& entry:pivot.row)entry.column=column.at(entry.column);for(auto& c:pivot.dependencies)c=column.at(c);}
    return out;
  }
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
      if(waiting>=lookahead||running()+1>=threads_||!room(kept_.size()+running()+2))break;
      std::vector<::ibp::Word> values;
      try{for(const auto& input:inputs_)values.push_back(modular::evaluate(input,expected.at.coordinates,expected.at.modulus));}catch(const std::exception&){continue;}
      auto solve=probe;solve.seconds=std::max(solve.seconds,std::chrono::duration<double>(expected.until-std::chrono::steady_clock::now()).count());solve.complete_rows=true;solve.progress={};
      auto task=std::make_shared<Elimination>(expected.at.modulus);task->system=plan_->system;task->max_fill=solve.max_fill;task->target_columns=plan_->program->program.target_columns;
      try{tasks_.emplace_back(std::jthread([task,source=std::shared_ptr<const ::ibp::ParametricProgram>(plan_->program),values=std::move(values),solve]{
          const auto started=std::chrono::steady_clock::now();bool solved=false;
          try{const auto& equations=source->program.equations;task->program.integrals.resize(source->program.integrals.size());task->program.equations.resize(equations.size());
            const auto forms=source->evaluate_forms(values,task->field);
            for(std::size_t i=0;i<equations.size();++i){auto& row=task->program.equations[i];row.reserve(equations[i].size());
              for(std::size_t j=0;j<equations[i].size();++j)row.push_back({equations[i][j].column,forms[source->coefficients[i][j]],0});}
            task->program.target_columns=task->target_columns;task->solver.emplace(task->program,task->field,solve);task->solver->solve(0);
            task->rows=task->solver->targets();task->selection=task->solver->selection();
            std::vector<std::vector<::ibp::LinearEntry>>().swap(task->program.equations);solved=true;
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
  // The corners each sector's seeds were widened to for the present batch of
  // targets: none at all, or a list for every sector.
  std::vector<Corners> widened_;std::vector<::ibp::Integral> widened_for_,widened_sorted_;unsigned widenings_=0;
  std::vector<Kept> kept_;std::size_t memory_=0,last_fill_=0;
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
    boost::json::object identity{{"science",science},{"parameter",parameter},{"counter","one-point-derivative-closure-v4"},{"dots",options_.dots},{"numerators",options_.numerators},{"dots_first",options_.dots_first},{"target_seeds",target_seeds_name(options_.target_seeds)},{"scaleless_sectors",options_.scaleless_sectors},{"seed_widenings",options_.max_seed_widenings}};
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
      if(elapsed>=limits.timeout_seconds)throw std::runtime_error("IBP master count time budget");
      auto receiving=limits;receiving.timeout_seconds=std::min(options_.probe_timeout_seconds,std::max(1u,static_cast<unsigned>(limits.timeout_seconds-elapsed)));
      auto reduced=(*sampler_)({demanded.begin(),demanded.end()},modular::point(dimension_.variable_count(),1,fire_modular::ordinals::planning_base),prime(1),receiving);
      if(!reduced.success)throw std::runtime_error("IBP master count failed: "+reduced.reason);
      std::set<ibp::Integral> found;for(const auto& [a,row]:reduced.residues.value())for(const auto& [master,c]:row)found.insert(master);
      const auto before=demanded.size();
      for(const auto& master:found){demanded.insert(master);for(const auto& [a,c]:generator.derivative(master,parameter))if(!c.is_zero())demanded.insert(a);}
      masters=found.size();stable=demanded.size()==before;
    }
    if(!stable)throw std::runtime_error("IBP master count did not stabilize");
    // The stable demands are what closure planning would find: leave them as
    // its hint, so that preparing the chosen level does not plan again.
    {
      boost::json::object plan_identity{{"science",science},{"parameter",parameter},{"planner","two-point-terminal-derivative-closure-v5"},{"dots",options_.dots},{"numerators",options_.numerators},{"dots_first",options_.dots_first},{"target_seeds",target_seeds_name(options_.target_seeds)},{"scaleless_sectors",options_.scaleless_sectors},{"seed_widenings",options_.max_seed_widenings}};
      const auto plan_key=artifacts::detail::sha256(artifacts::detail::canonical(plan_identity));
      const auto plan_path=options_.cache_directory/"closure-plans"/(plan_key+".json");
      if(!std::filesystem::exists(plan_path)){
        boost::json::array rows;for(const auto& a:demanded){boost::json::array row;for(int power:a)row.push_back(power);rows.push_back(row);}
        std::filesystem::create_directories(plan_path.parent_path());
        fire_modular::detail::save(plan_path,{{"schema","DiffExp.IBPClosurePlan/v1"},{"identity",plan_key},{"planning_hint_only",true},{"demands",rows},{"widened_seeds",widened_json()}});
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
    if(parameter>=dimension_.variable_count()||!options_.max_planning_passes)throw std::invalid_argument("invalid finite-field closure planning bounds");
    auto science=fire_batch::identity(basis_,dimension_,field_,requested,limits).json_value();
    boost::json::object identity{{"science",science},{"parameter",parameter},{"planner","two-point-terminal-derivative-closure-v5"},{"dots",options_.dots},{"numerators",options_.numerators},{"dots_first",options_.dots_first},{"target_seeds",target_seeds_name(options_.target_seeds)},{"scaleless_sectors",options_.scaleless_sectors},{"seed_widenings",options_.max_seed_widenings}};
    const auto key=artifacts::detail::sha256(artifacts::detail::canonical(identity));
    const auto path=options_.cache_directory/"closure-plans"/(key+".json");
    std::set<ibp::Integral,ibp::IntegralOrder> demanded(requested.begin(),requested.end());
    if(std::filesystem::exists(path)){
      auto cached=fire_modular::detail::read(path);
      if(artifacts::detail::string(cached.at("identity"))!=key||cached.at("schema")!="DiffExp.IBPClosurePlan/v1"||!cached.at("planning_hint_only").as_bool())throw std::runtime_error("IBP closure plan identity mismatch");
      std::set<ibp::Integral,ibp::IntegralOrder> planned;
      for(const auto& entry:cached.at("demands").as_array()){
        ibp::Integral a;for(const auto& power:entry.as_array())a.push_back(power.to_number<int>());fire_batch::validate_index(a,basis_);
        if(!planned.insert(a).second)throw std::runtime_error("duplicate closure plan demand");
      }
      for(const auto& a:requested)if(!planned.contains(a))throw std::runtime_error("closure plan omitted an original request");
      if(auto bounds=cached.if_contains("widened_seeds")){std::vector<std::array<unsigned,3>> wider;
        for(const auto& entry:bounds->as_array()){const auto& b=entry.as_array();if(b.size()!=3)throw std::runtime_error("closure plan seed bounds");wider.push_back({b[0].to_number<unsigned>(),b[1].to_number<unsigned>(),b[2].to_number<unsigned>()});}
        sampler_->widen({planned.begin(),planned.end()},wider);}
      ++planning_.cache_hits;planning_.demands=planned.size();if(options_.progress)options_.progress("reused derivative-closure planning hint; exact reconstruction and closure checks remain required");
      return {planned.begin(),planned.end()};
    }
    ibp::Generator generator(basis_,dimension_);
    // Every pass reduces its batch at the same two points, and the batches
    // only grow: the complete elimination of a pass continues from that of
    // the pass before it with the equations that are new, and needs none
    // while the seeds stay as they are (Sampler::continue_from). What can be
    // started ahead is the validation of the reconstruction that follows.
    std::vector<Sampler::Point> held;
    if(sampler_->works_ahead())try{held=complete(options_.selected_source_rows,options_.complete_validation_anchor,session_.held_out(limits,2));}catch(const std::exception&){}
    sampler_->anticipate(held,limits.timeout_seconds);
    for(unsigned pass=0;pass<options_.max_planning_passes;++pass){
      std::vector<ibp::Integral> batch(demanded.begin(),demanded.end());std::set<ibp::Integral> terminal;
      for(unsigned check=0;check<2;++check){
        auto elapsed=duration<double>(steady_clock::now()-started).count();if(elapsed>=limits.timeout_seconds)throw std::runtime_error("IBP closure planning time budget");
        auto receiving=limits;receiving.timeout_seconds=std::min(options_.probe_timeout_seconds,std::max(1u,static_cast<unsigned>(limits.timeout_seconds-elapsed)));
        auto point=modular::point(dimension_.variable_count(),1,fire_modular::ordinals::planning_base+check);
        auto reduced=(*sampler_)(batch,point,prime(1),receiving);++planning_.probes;
        if(!reduced.success)throw std::runtime_error("IBP closure planning probe failed: "+reduced.reason);
        std::set<ibp::Integral> found;for(const auto& [a,row]:reduced.residues.value())for(const auto& [master,c]:row)found.insert(master);
        if(check&&found!=terminal)throw std::runtime_error("IBP closure planning basis changed at a second point");terminal=std::move(found);
      }
      ++planning_.passes;const auto before=demanded.size();
      for(const auto& master:terminal){demanded.insert(master);for(const auto& [a,c]:generator.derivative(master,parameter))if(!c.is_zero())demanded.insert(a);}
      planning_.demands=demanded.size();
      if(options_.progress)options_.progress("finite-field closure plan "+std::to_string(pass+1)+": "+std::to_string(terminal.size())+" candidate masters, "+std::to_string(demanded.size())+" combined demands");
      if(demanded.size()==before){
        boost::json::array rows;for(const auto& a:demanded){boost::json::array row;for(int power:a)row.push_back(power);rows.push_back(row);}
        std::filesystem::create_directories(path.parent_path());
        fire_modular::detail::save(path,{{"schema","DiffExp.IBPClosurePlan/v1"},{"identity",key},{"planning_hint_only",true},{"demands",rows},{"widened_seeds",widened_json()}});
        return {demanded.begin(),demanded.end()};
      }
    }
    // The plan is a hint. One that has not stabilized is still a superset of
    // the request made of valid demands; it is not saved, and the level
    // preparer's own closure loop continues from its reconstruction.
    if(options_.progress)options_.progress("finite-field closure plan did not stabilize within "+std::to_string(options_.max_planning_passes)+" passes; reconstructing the "+std::to_string(demanded.size())+" demands found so far");
    return {demanded.begin(),demanded.end()};
  }
  static const char* target_seeds_name(Options::TargetSeeds seeds){return seeds==Options::TargetSeeds::none?"none":seeds==Options::TargetSeeds::sector?"sector":seeds==Options::TargetSeeds::subsectors?"subsectors":"all";}
  boost::json::array widened_json()const {boost::json::array out;for(const auto& [mask,dots,nums]:sampler_->widened())out.push_back(boost::json::array{mask,dots,nums});return out;}
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
    out.provider_identity=std::string("ibp-solver-parametric-")+(o.dots_first?"dots-first":"total-degree")+"-v2-dots-"+std::to_string(o.dots)+"-numerators-"+std::to_string(o.numerators);
    if(!o.scaleless_sectors)out.provider_identity+="-free-loop-zero-sectors-only";
    if(o.max_seed_widenings!=4)out.provider_identity+="-seed-widenings-"+std::to_string(o.max_seed_widenings);
    if(o.target_seeds==Options::TargetSeeds::subsectors)out.provider_identity+="-subsector-seeds-v1";
    if(o.target_seeds==Options::TargetSeeds::sector)out.provider_identity+="-sector-seeds-v1";
    if(o.target_seeds==Options::TargetSeeds::none)out.provider_identity+="-widened-seeds-v1";
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
  // The degree bounds also set how many samples a reconstruction takes, so
  // they start small; a caller that meets them raises them (see escalate).
  native.degree_slices=true;native.max_samples_per_prime=1024;
  native.closure_parameter=parameter;native.probe_timeout_seconds=elimination_seconds;
  if(!strict){native.selected_source_rows=true;native.complete_validation_anchor=true;}
  return native;
}
// Whether a failed reduction stopped at a degree or sample bound, and if so
// the options with those bounds doubled. Progress made under the old bounds
// is kept by the reconstruction itself.
inline std::optional<Options> escalate(const Options& options,const std::string& reason) {
  if(reason.find("degree budget")==std::string::npos && reason.find("sample/degree budget")==std::string::npos &&
     reason.find("sample budget")==std::string::npos)return std::nullopt;
  // No degree is the last one: the time and memory budgets end a reduction.
  if(options.max_degree>(1u<<30) || options.max_variable_degree>(1u<<30))return std::nullopt;
  auto next=options;
  next.max_degree=2*options.max_degree;
  next.max_variable_degree=2*options.max_variable_degree;
  next.max_samples_per_prime=2*options.max_samples_per_prime;
  return next;
}
} // namespace diffexp::ibp_solver
