#pragma once
#include "diffexp/fire_modular.hpp"
#include <ibp/trace.hpp>

namespace diffexp::ibp_solver {
// The pinned upstream core works in the exact denominator coordinates owned by
// DiffExp. No momentum reconstruction, renamed family, or new ISP completion.
struct Options {
  std::filesystem::path cache_directory;
  unsigned dots=1,numerators=2,max_degree=16,max_primes=8,probe_timeout_seconds=20;
  std::size_t max_samples_per_prime=512;
  std::size_t sample_checkpoint_interval=1;
  bool degree_slices=false,dots_first=false,target_sector_seeds=false,selected_source_rows=false;
  std::optional<std::size_t> closure_parameter;
  unsigned max_planning_passes=6;
  std::function<void(const std::string&)> progress;
};
struct Statistics {
  std::size_t probes=0,equations=0,templates=0,full_solves=0,trace_replays=0,trace_fallbacks=0;
  std::size_t selected_solves=0,selected_equations=0,selection_fallbacks=0;
  double generation_seconds=0,elimination_seconds=0,full_solve_seconds=0,trace_learning_seconds=0,trace_replay_seconds=0;
  double selected_solve_seconds=0;
};
inline modular::Word prime(unsigned index) {
  if(index<1||index>16)throw std::invalid_argument("IBP solver prime index must be 1..16");
  static const auto values=[] {std::array<modular::Word,16> out{};out[0]=2305843009213693951ULL;
    for(unsigned i=1;i<out.size();++i){auto p=out[i-1]-2;while(!n_is_prime(p))p-=2;out[i]=p;}return out;}();
  return values[index-1];
}
class Sampler {
 public:
  Sampler(ibp::PropagatorBasis basis,Exact dimension,Options options)
    :basis_(std::move(basis)),dimension_(std::move(dimension)),options_(std::move(options)) {
    if(!basis_.space.loops||basis_.space.loops>4||basis_.space.size()>16||basis_.physical_count>12||!basis_.physical_count||
       basis_.denominators.size()!=basis_.space.size()||basis_.physical_count>basis_.denominators.size())
      throw std::invalid_argument("IBP solver supports at most 4 loops, 16 scalar products and 12 physical denominators; select FIRE for larger systems");
    if(options_.dots>8||options_.numerators>8)throw std::invalid_argument("IBP solver seed powers must be 0..8");
    for(unsigned i=0;i<basis_.space.loops;++i)for(unsigned v=0;v<basis_.space.loops+basis_.space.externals();++v){
      contractions_.emplace_back();trace_.push_back(i==v);
      for(const auto& d:basis_.denominators)contractions_.back().push_back(basis_.rewrite(basis_.space.contraction(d,i,v)));
    }
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
  fire::Result operator()(const std::vector<ibp::Integral>& requested,const std::vector<modular::Word>& point,
      modular::Word modulus,const fire::Options& limits,bool independent=false) {
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
      if(!plan_||plan_->program.program.targets!=targets||plan_->zero_sectors!=geometry.zero_sectors||plan_->dots!=seeds.dots||plan_->numerators!=seeds.numerators||plan_->max_terms!=seeds.max_terms){
        plan_=std::make_unique<Plan>();plan_->program=::ibp::ParametricProgram::compile(geometry,seeds,targets,options_.dots_first?::ibp::ParametricOrdering::dots_first: ::ibp::ParametricOrdering::total_degree,sector_bounds);
        plan_->zero_sectors=geometry.zero_sectors;plan_->dots=seeds.dots;plan_->numerators=seeds.numerators;plan_->max_terms=seeds.max_terms;++stats_.templates;stats_.equations+=plan_->program.program.equations.size();
      }
      auto& parametric=plan_->program;auto& program=parametric.program;
      const auto generated=std::chrono::steady_clock::now();
      ::ibp::SolveOptions solve;solve.seconds=remaining();solve.max_fill=std::min<std::size_t>(12000000,limits.memory_bytes/64);
      std::vector<::ibp::Row> rows;auto& learned=plan_->primes[modulus];bool replayed=false;
      if(!independent&&learned.trace){const auto replay_begin=std::chrono::steady_clock::now();try{rows=learned.trace->evaluate_inputs(values,field,learned.workspace);++stats_.trace_replays;replayed=true;}
        catch(const std::domain_error&){learned.trace.reset();++stats_.trace_fallbacks;}stats_.trace_replay_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-replay_begin).count();}
      // Held-out samples use fresh primal elimination of original IBP rows,
      // independent of the arithmetic trace. A dependency subset proves the
      // same rowspace identities, but does not certify full-system minimality.
      // Check the first changed point too, before using a learned trace alone.
      if(independent||!replayed||!learned.checked){
        std::unique_ptr<::ibp::Solver> solver;std::vector<::ibp::Row> full;
        auto terminals=[](const auto& reduced){std::set<::ibp::Column> out;for(const auto& row:reduced)for(const auto& entry:row)out.insert(entry.column);return out;};
        auto primal=[&](bool selected){
          const auto started=std::chrono::steady_clock::now();
          const auto* sources=selected?&*plan_->selection:nullptr;
          // Rebind ALL original rows on fallback, including those left at an
          // earlier point or prime by the preceding subset solve.
          parametric.bind(values,field,sources);solve.seconds=remaining();
          solver=std::make_unique<::ibp::Solver>(program,field,solve);solver->solve(0,sources);full=solver->targets();
          auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
          if(selected){++stats_.selected_solves;stats_.selected_equations+=sources->size();stats_.selected_solve_seconds+=seconds;}
          else {++stats_.full_solves;stats_.full_solve_seconds+=seconds;}
        };
        bool selected=options_.selected_source_rows&&plan_->selection.has_value();primal(selected);
        if(selected&&terminals(full)!=plan_->terminal_columns){++stats_.selection_fallbacks;selected=false;primal(false);}
        if(options_.selected_source_rows&&!selected){
          plan_->selection=solver->selection();plan_->terminal_columns=terminals(full);
          if(options_.progress)options_.progress("full IBP pilot: "+std::to_string(program.equations.size())+" equations, "+std::to_string(plan_->selection->size())+" original rows needed for target identities");
        }
        if(replayed&&rows!=full){learned.trace.reset();++stats_.trace_fallbacks;replayed=false;}
        rows=std::move(full);
        if(!independent){
          if(replayed)learned.checked=true;
          else {const auto learning_begin=std::chrono::steady_clock::now();::ibp::ArithmeticTrace::Options tracing;tracing.seconds=remaining();tracing.max_fill=solve.max_fill;
            tracing.max_nodes=std::min<std::size_t>(std::numeric_limits<std::uint32_t>::max(),limits.memory_bytes/128);
            learned.trace=::ibp::ArithmeticTrace::learn_parametric(parametric,field,values,solver->selection(),tracing);
            if(learned.trace->evaluate_inputs(values,field,learned.workspace)!=rows)throw std::runtime_error("IBP parametric trace disagrees with primal learning reduction");
            learned.checked=false;stats_.trace_learning_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-learning_begin).count();}
        }
      }
      remaining();
      auto indices=[&](const ::ibp::Integral& a){ibp::Integral r;for(unsigned i=0;i<geometry.n;++i)r.push_back(a.powers[i]);return r;};
      std::set<ibp::Integral> terminal;
      for(unsigned i=0;i<rows.size();++i){ibp::Relation row;for(const auto& term:rows[i]){auto a=indices(program.integrals[term.column]);terminal.insert(a);ibp::add(row,a,dimension_.parse(std::to_string(term.value)));}out.reductions.emplace(requested[i],std::move(row));}
      for(const auto& a:terminal){auto [it,inserted]=out.reductions.emplace(a,ibp::Relation{{a,dimension_.constant(1)}});
        if(!inserted&&it->second!=ibp::Relation{{a,dimension_.constant(1)}})throw std::runtime_error("IBP solver returned a nonterminal output basis");}
      out.success=true;fire_batch::validate(out,basis_,dimension_,requested);
      ++stats_.probes;stats_.generation_seconds+=std::chrono::duration<double>(generated-begin).count();
      stats_.elimination_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-generated).count();
    }catch(const std::exception& e){out.success=false;out.reason=e.what();}
    return out;
  }
 private:
  std::size_t sector(const ibp::Integral& a)const {std::size_t mask=0;for(unsigned i=0;i<basis_.physical_count;++i)if(a[i]>0)mask|=std::size_t(1)<<i;return mask;}
  struct Learned {std::optional<::ibp::ArithmeticTrace> trace;::ibp::ArithmeticTrace::Workspace workspace;bool checked=false;};
  struct Plan {::ibp::ParametricProgram program;std::vector<bool> zero_sectors;unsigned dots=0,numerators=0;std::size_t max_terms=0;std::map<modular::Word,Learned> primes;
    std::optional<std::vector<std::uint32_t>> selection;std::set<::ibp::Column> terminal_columns;};
  ::ibp::InputGeometry input_geometry_;std::vector<Exact> inputs_;std::unique_ptr<Plan> plan_;
  ibp::PropagatorBasis basis_;Exact dimension_;Options options_;Statistics stats_;
  std::vector<std::vector<ibp::Affine>> contractions_;std::vector<bool> trace_,zero_sectors_;
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
    for(unsigned pass=0;pass<options_.max_planning_passes;++pass){
      if(demanded.size()>2000)throw std::length_error("IBP closure planning demand budget");
      std::vector<ibp::Integral> batch(demanded.begin(),demanded.end());std::set<ibp::Integral> terminal;
      for(unsigned check=0;check<2;++check){
        auto elapsed=duration<double>(steady_clock::now()-started).count();if(elapsed>=limits.timeout_seconds)throw std::runtime_error("IBP closure planning time budget");
        auto receiving=limits;receiving.timeout_seconds=std::min(options_.probe_timeout_seconds,std::max(1u,static_cast<unsigned>(limits.timeout_seconds-elapsed)));
        auto point=modular::point(dimension_.variable_count(),1,60000+2*pass+check);
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
    throw std::runtime_error("IBP derivative closure planning did not stabilize within its pass budget");
  }
  static fire_modular::Options configuration(const Options& o,const std::shared_ptr<Sampler>& sampler){
    fire_modular::Options out;out.cache_directory=o.cache_directory;out.max_degree=o.max_degree;out.max_primes=o.max_primes;
    out.max_samples_per_prime=o.max_samples_per_prime;out.probe_timeout_seconds=o.probe_timeout_seconds;out.progress=o.progress;out.modulus=prime;
    out.sample_checkpoint_interval=o.sample_checkpoint_interval;
    out.provider_identity=std::string("ibp-solver-parametric-")+(o.dots_first?"dots-first":"total-degree")+"-v1-dots-"+std::to_string(o.dots)+"-numerators-"+std::to_string(o.numerators);
    if(o.target_sector_seeds)out.provider_identity+="-subsector-seeds-v1";
    if(o.selected_source_rows)out.provider_identity+="-fresh-original-row-subset-validation-v1";
    out.prime_table="descending-61-bit-primes-v1";out.sparse_lifting=true;out.degree_slices=o.degree_slices;
    out.sample_provider=[sampler](const auto& r,const auto& p,auto m,const auto& limits){return (*sampler)(r,p,m,limits);};
    out.validation_provider=[sampler](const auto& r,const auto& p,auto m,const auto& limits){return (*sampler)(r,p,m,limits,true);};return out;
  }
  std::shared_ptr<Sampler> sampler_;fire_modular::Session session_;
  ibp::PropagatorBasis basis_;Exact dimension_;const ExactField& field_;Options options_;ClosurePlanningStatistics planning_;
};
} // namespace diffexp::ibp_solver
