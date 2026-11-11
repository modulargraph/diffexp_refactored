// Opt-in diagnostic, not a test: the denominator structure of the master basis
// of every recursion level. For each level it lists the masters, the distinct
// irreducible denominator factors of the connection A(x,eps) by kind (x only,
// eps only, mixed), the highest pole in eps, and what the bounded internal
// basis scan makes of them. Optional experiments try to remove the mixed
// factors: by diagonal rescaling, by replacing one master with another cached
// integral, and by the lattice closure of moving_poles.hpp. Preparation is
// that of `diffexp prepare`, through the same verified cache.
#include "diffexp/denominator_factors.hpp"
#include "diffexp/family_config.hpp"
#include "diffexp/ibp_solver_provider.hpp"
#include "diffexp/level_cache.hpp"
#include "diffexp/moving_poles.hpp"
#include "diffexp/recursion_pipeline.hpp"
#include <iostream>
using namespace diffexp;
namespace df=denominator_factors;
using Matrix=ExactEpsilonMatrix;
using Integral=diffexp::ibp::Integral;

static std::string name(const Integral& a) {
  std::string out="[";for(std::size_t i=0;i<a.size();++i)out+=(i?",":"")+std::to_string(a[i]);return out+"]";
}
struct Use {unsigned power=0;std::size_t entries=0;std::set<std::size_t> rows,columns;};
struct Summary {
  std::size_t nonzero=0,characters=0;long epsilon_pole=0;
  std::map<std::size_t,Use> use;
  std::size_t count(const df::Catalogue& c,df::Kind kind)const {
    std::size_t n=0;for(const auto& [k,u]:use)n+=c.factors()[k].kind==kind;return n;
  }
  std::size_t order(const df::Catalogue& c,df::Kind kind)const {
    std::size_t n=0;for(const auto& [k,u]:use)if(c.factors()[k].kind==kind)n+=u.power;return n;
  }
};
static Summary summarise(const Matrix& matrix,df::Catalogue& catalogue,std::size_t ei) {
  Summary s;
  for(std::size_t i=0;i<matrix.size();++i)for(std::size_t j=0;j<matrix[i].size();++j) {
    const auto& value=matrix[i][j];if(value.is_zero())continue;
    ++s.nonzero;s.characters+=value.str().size();
    s.epsilon_pole=std::max<long>(s.epsilon_pole,-*exact_epsilon_valuation(value,ei));
    for(const auto& [k,power]:catalogue.split(value.denominator())) {
      auto& u=s.use[k];u.power=std::max(u.power,power);++u.entries;u.rows.insert(i);u.columns.insert(j);
    }
  }
  return s;
}
static std::string indices(const std::set<std::size_t>& values) {
  std::string out="{";bool first=true;for(auto v:values){out+=(first?"":",")+std::to_string(v);first=false;}return out+"}";
}
static void print(const std::string& label,const Summary& s,const df::Catalogue& c,bool list) {
  std::cout<<"  "<<label<<": nonzero "<<s.nonzero<<", characters "<<s.characters<<", eps pole "<<s.epsilon_pole
    <<"; factors x-only "<<s.count(c,df::Kind::parameter)<<", eps-only "<<s.count(c,df::Kind::epsilon)
    <<", mixed "<<s.count(c,df::Kind::mixed)<<" (total order "<<s.order(c,df::Kind::mixed)<<")\n";
  if(!list)return;
  for(auto kind:{df::Kind::parameter,df::Kind::epsilon,df::Kind::mixed}) {
    std::cout<<"    "<<(kind==df::Kind::parameter?"x-only:":kind==df::Kind::epsilon?"eps-only:":"mixed:");
    if(kind!=df::Kind::mixed){for(const auto& [k,u]:s.use)if(c.factors()[k].kind==kind)std::cout<<"  ("<<c.factors()[k].polynomial.str()<<")^"<<u.power;std::cout<<'\n';continue;}
    std::cout<<'\n';
    for(const auto& [k,u]:s.use)if(c.factors()[k].kind==kind)
      std::cout<<"      ("<<c.factors()[k].polynomial.str()<<")^"<<u.power<<"  in "<<u.entries<<" entries, rows "<<indices(u.rows)<<" columns "<<indices(u.columns)<<'\n';
  }
}
static recursion::PhysicalBasisTransform unchanged(const Matrix& matrix) {
  const auto unit=fuchsify::detail::identity(matrix.size(),matrix[0][0]);
  return {matrix,unit,unit,matrix};
}

// (a) Diagonal rescaling J_i -> f^(-s_i) J_i by a mixed factor f. Off the
// diagonal the order of f changes by s_j-s_i, so the shifts solve the same
// difference constraints as the epsilon gauge; the diagonal acquires
// -s_i f'/f, which only an existing simple pole of that residue cancels.
static void rescaling(const Matrix& matrix,df::Catalogue& catalogue,const Summary& summary,std::size_t xi,std::size_t ei) {
  const auto d=matrix.size();const auto& z=matrix[0][0];
  for(const auto& [k,use]:summary.use) {
    if(catalogue.factors()[k].kind!=df::Kind::mixed)continue;
    const auto f=catalogue.factors()[k].polynomial;
    Matrix orders(d,std::vector<Exact>(d,z.constant(0)));
    for(std::size_t i=0;i<d;++i)for(std::size_t j=0;j<d;++j)if(i!=j && !matrix[i][j].is_zero())
      orders[i][j]=epsilon_gauge_detail::multiply_power(z.constant(1),ei,catalogue.order(matrix[i][j],k));
    std::cout<<"  rescaling by ("<<f.str()<<"): ";
    try {
      auto gauge=epsilon_diagonal_gauge(orders,ei);
      auto basis=unchanged(matrix);
      std::string shifts;
      for(std::size_t i=0;i<d;++i) {
        shifts+=(i?",":"")+std::to_string(gauge.shifts[i]);
        if(!gauge.shifts[i])continue;
        std::vector<Exact> q(d,z.constant(0));q[i]=z.constant(1)/f.pow(static_cast<unsigned long>(gauge.shifts[i]));
        moving_poles::replace_row(basis,i,q,xi);
      }
      recursion::verify_physical_basis(basis,matrix,xi);
      auto after=summarise(basis.connection,catalogue,ei);
      const auto found=after.use.find(k);
      std::cout<<"shifts {"<<shifts<<"}; order of this factor "<<use.power<<" -> "<<(found==after.use.end()?0u:found->second.power)
        <<" (in "<<(found==after.use.end()?std::size_t(0):found->second.entries)<<" entries); all mixed: "<<after.count(catalogue,df::Kind::mixed)
        <<" factors, total order "<<after.order(catalogue,df::Kind::mixed)<<'\n';
    } catch(const UnsupportedEpsilonGauge&) {
      std::cout<<"no diagonal rescaling clears the off-diagonal entries (the pole orders contain a cycle)\n";
    }
  }
}
// (b) One master replaced by another integral whose exact reduction the level
// already holds (requested integrals and the integrals of the cached IBP
// identities): every integral at every position, ranked by the total order of
// the mixed factors that remain.
static void replacements(const recursion::Node& node,df::Catalogue& catalogue,std::size_t xi) {
  const auto& level=node.closure;const auto d=level.ordered_basis.size();const auto& z=level.matrix[0][0];
  std::vector<std::pair<Integral,std::vector<Exact>>> rows;
  std::set<Integral> seen(level.ordered_basis.begin(),level.ordered_basis.end());
  for(std::size_t i=0;i<level.target_rows.size();++i)if(seen.insert(level.requested_integrals[i]).second)rows.emplace_back(level.requested_integrals[i],level.target_rows[i]);
  diffexp::ibp::ExactReducer reducer(z,level.source_identities.size());
  std::set<Integral> possible;
  for(const auto& relation:level.source_identities){reducer.insert(relation);for(const auto& [a,c]:relation)possible.insert(a);}
  diffexp::ibp::BasisReduction coordinates(reducer,level.ordered_basis,z);
  for(const auto& a:possible)if(seen.insert(a).second)try{rows.emplace_back(a,coordinates.resolve({{a,z.constant(1)}}));}catch(const std::runtime_error&){}
  const auto original=df::mixed_weight(level.matrix,catalogue);
  struct Trial {std::size_t order,characters,position;Integral integral;bool pole;};
  std::vector<Trial> trials;std::size_t with_pole=0;
  for(const auto& [integral,row]:rows) {
    bool pole=false;
    for(const auto& value:row)if(!value.is_zero())for(const auto& [k,power]:catalogue.split(value.denominator()))pole=pole||catalogue.factors()[k].kind==df::Kind::mixed;
    with_pole+=pole;
    // Large levels: only the integrals that can change the lattice at a mixed
    // factor (see moving_poles.hpp); on small ones every integral is tried.
    if(!pole && d>16)continue;
    for(std::size_t position=0;position<d;++position) {
      if(row[position].is_zero())continue;
      auto basis=unchanged(level.matrix);
      moving_poles::replace_row(basis,position,row,xi);
      const auto weight=df::mixed_weight(basis.connection,catalogue);
      trials.push_back({weight.total_order,physical_basis_scan::expression_characters(basis.connection),position,integral,pole});
    }
  }
  std::sort(trials.begin(),trials.end(),[](const Trial& a,const Trial& b){return std::tie(a.order,a.characters)<std::tie(b.order,b.characters);});
  std::size_t better=0,clean=0;for(const auto& trial:trials){better+=trial.order<original.total_order;clean+=!trial.order;}
  std::cout<<"  single replacements: "<<rows.size()<<" integrals with a cached reduction ("<<with_pole<<" with a mixed factor in it), "<<trials.size()
    <<" trials; original mixed order "<<original.total_order<<"; "<<better<<" trials lower it, "<<clean<<" remove it\n";
  for(std::size_t i=0;i<trials.size() && i<6;++i)
    std::cout<<"    master "<<trials[i].position<<" "<<name(level.ordered_basis[trials[i].position])<<" -> "<<name(trials[i].integral)
      <<(trials[i].pole?" (reduction has a mixed factor)":"")<<": mixed order "<<trials[i].order<<", characters "<<trials[i].characters<<'\n';
}

int main(int argc,char** argv) {
  try {
    if(argc<3)throw std::invalid_argument("usage: master_basis_report CONFIGURATION --cache DIRECTORY [--list 0|1] [--scan SECONDS] [--scan-trials N] [--rescale 0|1] [--replace 0|1] [--remove 0|1] [--closure 0|1] [--integrals 0|1] [--fixed 0|1] [--profile 0|1] [--profile-level N] [--levels N] [--ibp-dots N] [--ibp-numerators N] [--dots-first 0|1]");
    auto configuration=family_config::load(argv[1]);
    std::filesystem::path cache;double scan_seconds=0;bool list=true,dots_first=false,remove=false,closure=true,use_integrals=true,fixed=false,rescale=false,replace=false,profile=false;
    unsigned dots=configuration.ibp_dots,numerators=configuration.ibp_numerators;std::size_t levels=100,profile_level=0,scan_trials=0;
    for(int i=2;i+1<argc;i+=2) {
      const std::string option=argv[i],value=argv[i+1];
      if(option=="--cache")cache=value;
      else if(option=="--scan")scan_seconds=std::stod(value);
      else if(option=="--scan-trials")scan_trials=std::stoul(value);
      else if(option=="--list")list=value!="0";
      else if(option=="--remove")remove=value!="0";
      else if(option=="--closure")closure=value!="0";
      else if(option=="--integrals")use_integrals=value!="0";
      else if(option=="--fixed")fixed=value!="0";
      else if(option=="--rescale")rescale=value!="0";
      else if(option=="--replace")replace=value!="0";
      else if(option=="--profile")profile=value!="0";
      else if(option=="--levels")levels=std::stoul(value);
      else if(option=="--profile-level")profile_level=std::stoul(value);
      else if(option=="--ibp-dots")dots=std::stoul(value);
      else if(option=="--ibp-numerators")numerators=std::stoul(value);
      else if(option=="--dots-first")dots_first=value!="0";
      else throw std::invalid_argument("unknown option "+option);
    }
    if(cache.empty())throw std::invalid_argument("--cache is required");
    auto options=configuration.preparation;
    options.reduction.batch_cache_directory=cache/"fire-batches";options.reduction.pending_cache_directory=cache/"fire-pending";
    artifacts::Store store(cache);
    auto provider=[&](const auto& basis,const auto& dimension,const auto& field,auto parameter,const auto& sources,const auto& budget) {
      std::shared_ptr<ibp_solver::Session> session;
      auto native=ibp_solver::recursion_options(cache/"ibp-solver",parameter,options.reduction.provider.timeout_seconds,false);
      native.dots=dots;native.numerators=numerators;native.dots_first=dots_first;
      level::Provider fallback=[&,native](const auto& batch,const auto& limits){
        if(!session)session=std::make_shared<ibp_solver::Session>(basis,dimension,field,native);return (*session)(batch,limits);};
      return std::move(cached_level::prepare(store,basis,dimension,field,parameter,sources,budget,fallback).result);
    };
    auto graph=recursion::prepare(configuration.family,configuration.integrals,options,provider);
    const std::size_t xi=0,ei=1;
    std::map<std::size_t,recursion::PhysicalBasisTransform> bases;
    std::cout<<"family "<<graph.family_name<<": "<<graph.nodes.size()<<" levels\n";
    for(std::size_t depth=0;depth<graph.nodes.size() && depth<levels;++depth) {
      const auto& node=graph.nodes[depth];
      if(node.scalar_leaf||node.closure.matrix.empty()){std::cout<<"level "<<depth+1<<": scalar leaf\n";continue;}
      df::Catalogue catalogue(xi,ei);
      const auto d=node.closure.ordered_basis.size();
      std::cout<<"level "<<depth+1<<": "<<node.merged.physical_count<<" propagators, "<<d<<" masters, "
        <<node.operations.size()<<" observables, anchor "<<node.anchor.str()<<'\n';
      if(list){std::cout<<"  masters:";for(std::size_t i=0;i<d;++i)std::cout<<' '<<i<<':'<<name(node.closure.ordered_basis[i]);std::cout<<'\n';}
      const auto connection=summarise(node.closure.matrix,catalogue,ei);
      print("connection",connection,catalogue,list);
      print("observable rows",summarise(node.observable_rows,catalogue,ei),catalogue,false);
      try {
        auto gauge=epsilon_diagonal_gauge(node.closure.matrix,ei);
        std::cout<<"  eps gauge: largest shift "<<*std::max_element(gauge.shifts.begin(),gauge.shifts.end())<<'\n';
      } catch(const UnsupportedEpsilonGauge&) {std::cout<<"  eps gauge: none (negative valuation cycle)\n";}
      if(scan_seconds>0) {
        physical_basis_scan::Options scan;scan.seconds=scan_seconds;scan.epsilon_variable=ei;scan.max_candidates=8;
        if(scan_trials)scan.max_trials=scan_trials;
        auto frontier=physical_basis_scan::scan_candidates(node.closure,node.merged.physical_count,xi,scan,
          [](const physical_basis_scan::Candidate& candidate)->std::optional<physical_basis_scan::DemandCost> {
            return physical_basis_scan::DemandCost{candidate.metadata.connection_characters};});
        std::cout<<"  scan (connection size, as the pipeline ranks): trials "<<frontier.stats.trials<<", accepted "<<frontier.stats.accepted
          <<", new poles "<<frontier.stats.new_poles<<", rows "<<frontier.stats.candidates<<", "<<frontier.stats.seconds<<" s"
          <<(frontier.stats.time_expired?" (time budget reached)":"")<<'\n';
        for(const auto& candidate:frontier.candidates) {
          if(candidate.metadata.provenance=="original")continue;
          print("  candidate ["+candidate.metadata.provenance+"]",summarise(candidate.transform.connection,catalogue,ei),catalogue,false);
        }
      }
      if(rescale)rescaling(node.closure.matrix,catalogue,connection,xi,ei);
      if(replace && connection.count(catalogue,df::Kind::mixed))replacements(node,catalogue,xi);
      if(remove) {
        moving_poles::Options removal;removal.closure_rows=closure;removal.fixed_points=fixed;
        const auto integrals=use_integrals?moving_poles::cached_integrals(node.closure):std::vector<moving_poles::Row>{};
        auto removed=moving_poles::remove(node.closure.matrix,xi,ei,integrals,removal);
        std::cout<<"  moving-pole removal ("<<integrals.size()<<" cached integrals"<<(closure?", closure rows":"")<<(fixed?", fixed points in x too":"")<<"): "<<removed.steps.size()<<" steps, "<<removed.seconds<<" s, "<<(fixed?"target":"mixed")<<" order "
          <<removed.before.total_order<<" -> "<<removed.after.total_order<<", characters "<<removed.original_characters<<" -> "<<removed.characters<<(removed.budget_reached?" (budget reached)":"")<<'\n';
        for(const auto& step:removed.steps)
          std::cout<<"    ("<<step.factor<<")^"<<step.order<<": master "<<step.replaced<<" "<<name(node.closure.ordered_basis[step.replaced])<<" replaced by "<<step.source<<" ("<<step.trials<<" trials)\n";
        for(const auto& factor:removed.remaining)std::cout<<"    left: ("<<factor<<")\n";
        if(removed.changed()) {
          recursion::verify_physical_basis(removed.transform,node.closure.matrix,xi);
          print("  new connection",summarise(removed.transform.connection,catalogue,ei),catalogue,list);
          print("  new observable rows",summarise(fuchsify::detail::multiply(node.observable_rows,removed.transform.from_basis),catalogue,ei),catalogue,false);
          print("  to_basis",summarise(removed.transform.to_basis,catalogue,ei),catalogue,false);
          print("  from_basis",summarise(removed.transform.from_basis,catalogue,ei),catalogue,false);
          bases.emplace(depth,std::move(removed.transform));
        }
      }
    }
    if(profile) {
      // Structural epsilon losses of every level: endpoint operator lows,
      // forcing low, anchor gauge loss, and the child window they add up to.
      for(int variant=0;variant<(bases.empty()?1:2);++variant) {
        auto numerical=configuration.numerical;numerical.automatic_basis_scan=false;
        if(variant)numerical.physical_bases=bases;
        numerical.progress=[](std::size_t depth,const std::string& phase,int high){std::cerr<<"  level "<<depth+1<<": "<<phase<<"; epsilon through "<<high<<'\n';};
        recursion::Evaluator evaluator(graph,numerical);
        std::cout<<(variant?"epsilon demand, moving poles removed:\n":"epsilon demand, original bases:\n");
        for(std::size_t depth=0;depth<graph.nodes.size() && depth<levels;++depth) {
          if(graph.nodes[depth].scalar_leaf||graph.nodes[depth].closure.matrix.empty()||(profile_level && depth+1!=profile_level))continue;
          const auto began=std::chrono::steady_clock::now();
          auto p=evaluator.epsilon_demand_profile(depth);
          std::cout<<"  level "<<depth+1<<": lower endpoint low "<<p.lower_endpoint_low<<", upper "<<p.upper_endpoint_low<<", forcing low "<<p.beta_low
            <<", anchor loss "<<p.anchor_inverse_pole_loss<<"; child loss common "<<p.adjoint_child_loss<<", by component "<<p.factored_child_loss
            <<" ("<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<" s)\n"<<std::flush;
        }
      }
    }
    return 0;
  } catch(const std::exception& error) {std::cerr<<"master_basis_report: "<<error.what()<<'\n';return 1;}
}
