#pragma once
#include <atomic>
#include <thread>
#include "diffexp/scaleless_sectors.hpp"
#include "diffexp/families.hpp"
#include "diffexp/level_preparation.hpp"
#include "diffexp/merge_pullback.hpp"
#include "diffexp/geometry.hpp"
#include "diffexp/causal.hpp"
#include "diffexp/endpoint_boundary.hpp"
#include <algorithm>
#include <map>
#include <numeric>
#include <limits>

namespace diffexp::recursion {
// Where a level takes its boundary data: at an interior anchor of the merge
// parameter, where the child is the merged family at that value, or at an
// endpoint, where the child is the incoming family with one of the two merged
// lines removed (endpoint_boundary.hpp).
enum class BoundarySite { anchor, lower, upper };
// anchor: every level matches at its interior anchor. automatic: a level
// whose endpoint qualifies matches there. endpoint: as automatic, and a
// preparation in which no level qualifies is an error.
enum class BoundaryPolicy { anchor, automatic, endpoint };
struct Options {
  level::Options reduction;
  pullback::Budget numerator_expansion;
  std::vector<Rational> anchors;
  std::vector<std::pair<unsigned,unsigned>> merges;
  unsigned max_levels=16,total_timeout_seconds=86400;
  std::size_t max_sources_per_level=5000;
  // Optional exchange of masters of a freshly closed level (zero-based depth)
  // for other integrals, before its masters become the child's requests. The
  // returned closure must describe the same requested integrals in its own
  // ordered basis; see moving_poles::exchange.
  std::function<level::Result(level::Result,unsigned)> exchange_masters;
  // Optional choice of the masters of a freshly closed level that may ask
  // for further reductions (finite_masters.hpp). It is given the closure for
  // the level's sources, the merged family, the field, the sources, the zero-based depth
  // and a function that closes the level for a longer list of requests, the
  // sources first, within a number of seconds; it returns a closure of the
  // same size with rows for the sources only. It runs before exchange_masters.
  std::function<level::Result(level::Result,const ibp::PropagatorBasis&,const ExactField&,const std::vector<ibp::Integral>&,unsigned,
    const std::function<level::Result(const std::vector<ibp::Integral>&,double)>&)> choose_masters;
  // Which two propagators a level merges. Without explicit `merges` and with
  // this counter, the pairs with a common loop momentum are tried at every
  // level and the one whose merged family has the fewest master integrals is
  // taken (the first pair unless
  // another has strictly fewer): the masters of a level are the rows its
  // child transports, and the child of the top level is where a run spends
  // its time. The counter returns nothing if it cannot count a candidate.
  // Without it, each level merges its first two propagators.
  std::function<std::optional<std::size_t>(const ibp::PropagatorBasis&,const Exact&,const ExactField&,std::size_t,
    const std::vector<ibp::Integral>&,const level::Options&)> count_masters;
  // The same choice from the graph alone: the number of masters of the family
  // that merges lines `left` and `right`, given each present line as a
  // combination of the original propagators (critical_points.hpp). It is
  // cheap, so every pair of every level is counted; where it answers, it
  // decides and count_masters is not asked.
  // The last arguments are the seconds one count may take and whether it was
  // given up for lack of time (now or, remembered, by an earlier run) rather
  // than failing for another reason. The first pair is
  // counted alone within selection_seconds. If it does not finish, the level
  // merges the first pair at once: the count only ranks the candidates and
  // is not worth more than that. If it does, the other pairs are counted on
  // selection_threads threads (0: one per hardware thread), each within
  // selection_seconds plus four times what the first pair took.
  std::function<std::optional<std::size_t>(const std::vector<std::vector<Rational>>&,unsigned,unsigned,unsigned,bool&)> count_critical_points;
  unsigned selection_seconds=10,selection_threads=0;
  // Told what the selection of a level is doing and how it ended.
  std::function<void(const std::string&)> selection_progress;
  // Told the choice of a level (one-based depth, pair, candidates, masters of
  // the pair and of the first pair) before that level is reduced.
  std::function<void(unsigned,unsigned,unsigned,std::size_t,std::size_t,std::size_t)> merge_chosen;
  // IBP counting only: alternatives are counted only when the first pair
  // gives a level of at least this many masters.
  std::size_t merge_selection_threshold=24;
  // Explicit anchors keep every level at its anchor whatever the policy.
  BoundaryPolicy boundary=BoundaryPolicy::anchor;
  endpoint_boundary::Options endpoint;
  // Told why an endpoint of a level (one-based depth) that passed the
  // kinematic rule was refused by the level's equation or its numerators.
  std::function<void(unsigned,bool,const std::string&)> endpoint_refused;
};
struct Event {
  unsigned depth,physical_propagators;
  std::size_t requests,sources,masters,provider_passes;
  double elapsed_seconds;
  bool scalar_leaf;
  Rational requested_anchor,anchor;
  // The merged pair, the number of pairs considered (1: no selection), and
  // the master counts of the chosen pair and of the first pair.
  unsigned left=0,right=1;std::size_t candidates=1,counted_masters=0,first_pair_masters=0;
  BoundarySite boundary=BoundarySite::anchor;
};
using PrepareLevel=std::function<level::Result(const ibp::PropagatorBasis&,const Exact&,
  const ExactField&,std::size_t,const std::vector<ibp::Integral>&,const level::Options&)>;
struct Node {
  ibp::PropagatorBasis incoming,merged;
  unsigned left,right;
  Rational anchor;
  std::vector<ibp::Integral> requested,source_indices;
  std::vector<pullback::Plan> operations;
  level::Result closure;
  std::vector<std::vector<Exact>> observable_rows;
  bool scalar_leaf=false;
  // At an endpoint the anchor is 0 or 1 and the child computes the masters'
  // values there, which are integrals of the pinched family.
  BoundarySite boundary=BoundarySite::anchor;
  // Dimension of the space of analytic solutions at that endpoint, from the
  // level's equation with symbolic epsilon.
  std::size_t analytic_solutions=0;
};
struct Graph {
  ExactField field;
  Exact dimension;
  std::string family_name;
  std::vector<Node> nodes;
  std::vector<Event> events;
  std::optional<feynman::ExampleFamily> definition;
  // The epsilon variable of every exact coefficient is epsilon_total minus
  // this shift: the dimension is d0-2*(epsilon_shift+eps). See shift_epsilon.
  Rational epsilon_shift=Rational(0);
  // The non-leaf levels hold no epsilon: they are the recursion at the one
  // dimension d0-2*epsilon_shift.
  bool specialised=false;
};

// The same prepared recursion at epsilon_total=shift instead of zero. At a
// generic shift nothing is singular there: no level has epsilon poles, and
// the integral in d0-2*shift dimensions follows from the equations at that
// dimension alone. Every level above the scalar leaf is therefore specialised
// to the number (no epsilon is left in its connection or observables, so its
// endpoint series, matching and transport work with rational functions of x
// only); the leaf keeps eps=shift+eps', which its closed form needs. With
// `specialise` false every level keeps the residual epsilon and higher
// orders in (epsilon-shift) can be requested.
// Sampling several shifts and fitting the Laurent coefficients replaces the
// widening epsilon windows of the expanded recursion.
inline Graph shift_epsilon(Graph graph,const Rational& shift,bool specialise=true) {
  const auto names=graph.dimension.variables();
  const auto found=std::find(names.begin(),names.end(),"eps");
  if(found==names.end())throw std::invalid_argument("epsilon shift needs the eps variable");
  const std::size_t epsilon=found-names.begin();
  std::vector<Exact> moved,fixed;
  for(std::size_t i=0;i<graph.dimension.variable_count();++i)moved.push_back(graph.dimension.variable(i));
  fixed=moved;
  moved[epsilon]=moved[epsilon]+graph.dimension.constant(shift);
  fixed[epsilon]=graph.dimension.constant(shift);
  graph.dimension=graph.dimension.substitute(moved);
  for(auto& node:graph.nodes) {
    const auto& point=specialise&&!node.scalar_leaf?fixed:moved;
    const auto move=[&](Exact& value){if(!value.is_zero())value=value.substitute(point);};
    for(auto* matrix:{&node.closure.matrix,&node.closure.target_rows,&node.observable_rows})
      for(auto& row:*matrix)for(auto& value:row)move(value);
    for(auto& operation:node.operations)for(auto& [integral,value]:operation.source_row)move(value);
    for(auto& identity:node.closure.source_identities)for(auto& [integral,value]:identity)move(value);
  }
  graph.epsilon_shift=graph.epsilon_shift+shift;graph.specialised=specialise;
  return graph;
}

// Test the denominator's leading epsilon coefficient, before specializing the
// rational function. For example eps/(x-a+eps) still forbids an anchor at a:
// its epsilon-zero value alone would hide the pole of higher coefficients.
inline bool regular_anchor(const std::vector<std::vector<Exact>>& matrix,
    std::size_t parameter,std::size_t epsilon,const Rational& anchor) {
  for(const auto& row:matrix)for(const auto& coefficient:row)if(!coefficient.is_zero()) {
    auto denominator=coefficient.denominator();
    auto valuation=std::numeric_limits<unsigned long>::max();
    for(const auto& term:denominator.numerator_terms())valuation=std::min(valuation,static_cast<unsigned long>(term.powers.at(epsilon)));
    if(valuation)denominator=denominator/denominator.variable(epsilon).pow(valuation);
    std::vector<Exact> point;for(std::size_t i=0;i<coefficient.variable_count();++i)point.push_back(coefficient.variable(i));
    point.at(parameter)=coefficient.constant(anchor);point.at(epsilon)=coefficient.constant(0);
    if(denominator.substitute(point).is_zero())return false;
  }
  return true;
}

// Real singular points of the epsilon-leading denominators, as isolating
// balls; the other roots are appended to `complex` if it is given.
inline std::vector<kernel::ComplexBall> real_singular_points(const std::vector<std::vector<Exact>>& matrix,
    std::size_t parameter,std::size_t epsilon,std::vector<kernel::ComplexBall>* complex=nullptr) {
  std::map<std::string,Exact> polynomials;
  for(const auto& row:matrix)for(const auto& coefficient:row)if(!coefficient.is_zero()) {
    auto denominator=coefficient.denominator();
    auto valuation=std::numeric_limits<unsigned long>::max();
    for(const auto& term:denominator.numerator_terms())valuation=std::min(valuation,static_cast<unsigned long>(term.powers.at(epsilon)));
    if(valuation)denominator=denominator/denominator.variable(epsilon).pow(valuation);
    std::vector<Exact> point;for(std::size_t k=0;k<coefficient.variable_count();++k)point.push_back(coefficient.variable(k));
    point.at(epsilon)=coefficient.constant(0);
    auto leading=denominator.substitute(point);
    if(!leading.is_rational())polynomials.emplace(leading.str(),std::move(leading));
  }
  std::vector<kernel::ComplexBall> out;
  for(const auto& [text,polynomial]:polynomials)
    for(auto& root:polynomial_roots(polynomial,static_cast<unsigned>(parameter))) {
      if(arb_contains_zero(acb_imagref(root.raw())))out.push_back(std::move(root));
      else if(complex)complex->push_back(std::move(root));
    }
  return out;
}

// A requested anchor that meets a pole of the closed system is moved: which
// points are singular depends on the master basis the reduction selected,
// which the request cannot know in advance. The replacement is a regular
// rational within 3/20 with no other real singular point of the system
// between it and the request, so the child family stays in the same real
// cell of the parameter and a supplied causal prescription keeps its meaning.
// Small denominators come first, because the anchor enters the kinematics of
// every deeper level and its height multiplies through their reductions; the
// other requested anchors are left to their own levels.
inline Rational select_regular_anchor(const std::vector<std::vector<Exact>>& connection,
    const std::vector<std::vector<Exact>>& observables,std::size_t parameter,std::size_t epsilon,
    const Rational& preferred,const std::vector<Rational>& requested={}) {
  const auto regular=[&](const Rational& at){return regular_anchor(connection,parameter,epsilon,at)&&regular_anchor(observables,parameter,epsilon,at);};
  if(regular(preferred))return preferred;
  if(!requested.empty()) {
    auto singular=real_singular_points(connection,parameter,epsilon);
    for(auto& root:real_singular_points(observables,parameter,epsilon))singular.push_back(std::move(root));
    const auto bits=kernel::ComplexBall::precision();
    const auto same_cell=[&](const Rational& candidate) {
      struct Real {arb_t value;Real(){arb_init(value);}~Real(){arb_clear(value);}} low,high,at;
      arb_set_fmpq(low.value,(candidate<preferred?candidate:preferred).raw(),bits);
      arb_set_fmpq(high.value,(candidate<preferred?preferred:candidate).raw(),bits);
      arb_set_fmpq(at.value,preferred.raw(),bits);
      for(const auto& root:singular) {
        const auto* real=acb_realref(root.raw());
        if(arb_overlaps(real,at.value))continue; // the pole at the request itself
        if(!arb_lt(real,low.value)&&!arb_gt(real,high.value))return false;
      }
      return true;
    };
    const Rational reach("3/20");
    for(unsigned denominator=2;denominator<=20;++denominator) {
      std::vector<std::pair<Rational,Rational>> row;
      for(unsigned numerator=1;numerator<denominator;++numerator) {
        if(std::gcd(numerator,denominator)!=1)continue;
        const auto candidate=Rational(numerator)/Rational(denominator);
        const auto distance=candidate<preferred?preferred-candidate:candidate-preferred;
        if(distance>reach||std::find(requested.begin(),requested.end(),candidate)!=requested.end())continue;
        row.emplace_back(distance,candidate);
      }
      std::sort(row.begin(),row.end(),[](const auto& a,const auto& b){return a.first<b.first;});
      for(const auto& [distance,candidate]:row)if(regular(candidate)&&same_cell(candidate))return candidate;
    }
    throw std::domain_error("requested recursive anchor "+preferred.str()+" is singular in a connection or observable epsilon coefficient, and no simple regular rational within 3/20 shares its real cell");
  }
  unsigned attempts=0;
  for(unsigned denominator=2;attempts<64;++denominator)for(unsigned numerator=1;numerator<denominator&&attempts<64;++numerator) {
    if(std::gcd(numerator,denominator)!=1)continue;
    ++attempts;Rational candidate=Rational(numerator)/Rational(denominator);
    if(regular(candidate))return candidate;
  }
  throw std::runtime_error("no regular recursive anchor found within the finite rational-candidate budget");
}

// Freeze the ACTUAL completed numerator coordinates at the matching anchor.
// Re-completing an earlier raw family here would silently change the meaning of
// master indices when a rank pivot changes between recursive levels.
inline ibp::QuadraticFamily at_anchor(const ibp::PropagatorBasis& basis,
    std::size_t parameter,const Rational& anchor) {
  auto sample=basis.space.zero().constant;
  std::vector<Exact> point;for(std::size_t i=0;i<sample.variable_count();++i)point.push_back(sample.variable(i));
  point.at(parameter)=sample.constant(anchor);
  auto gram=basis.space.external_gram;
  for(auto& row:gram)for(auto& c:row)c=c.substitute(point);
  ibp::QuadraticFamily out{ibp::ScalarProducts(basis.space.loops,std::move(gram),sample),{}, {}};
  for(std::size_t i=0;i<basis.denominators.size();++i) {
    auto d=basis.denominators[i];d.constant=d.constant.substitute(point);
    for(auto& c:d.linear)c=c.substitute(point);
    (i<basis.physical_count?out.physical:out.auxiliary).push_back(std::move(d));
  }
  return out;
}

// The family of an endpoint child: the merged family at x=0 or x=1. There the
// merged line is one of the two original lines, and a numerator slot that was
// independent of the propagators at generic x can depend on them (the slot
// k^2 next to a line 1-k^2). Such a slot is replaced by a scalar product that
// completes the basis, which changes nothing if no master carries a power in
// it; a slot some master uses must stay independent, or there is no family.
// Whether F of `family` contains every monomial of (parameter of `line`)*U with
// a coefficient of the sign all of F has: F >= c*b*U on the simplex for some
// c>0. A massive line gives m^2*b*U, but on-shell external momenta can cancel
// monomials of it (a massive line next to a leg of the same mass); the
// integrals are then not analytic where the other merged line is pinched.
inline bool endpoint_mass_bound(const ibp::QuadraticFamily& family,std::size_t line) {
  const auto polynomials=scaleless::symanzik(ibp::PropagatorBasis(family));
  if(polynomials.U.empty() || polynomials.F.empty())return false;
  int sign=0;
  for(const auto& [exponent,coefficient]:polynomials.F) {
    if(!coefficient.is_rational())return false;
    const int here=coefficient.rational().sign();
    if(!sign)sign=here;
    if(here!=sign)return false;
  }
  for(const auto& [exponent,coefficient]:polynomials.U) {
    auto shifted=exponent;
    if(line>=shifted.size())return false;
    ++shifted[line];
    if(!polynomials.F.count(shifted))return false;
  }
  return true;
}
inline std::optional<ibp::QuadraticFamily> endpoint_family(const ibp::PropagatorBasis& merged,
    std::size_t parameter,const Rational& endpoint,const std::vector<ibp::Integral>& masters) {
  auto family=at_anchor(merged,parameter,endpoint);
  const auto zero=merged.space.zero().constant;
  const auto size=merged.space.size();
  std::map<std::size_t,std::vector<Exact>> rank;
  const auto reduced=[&](std::vector<Exact> row) {
    for(const auto& [pivot,basis]:rank) {
      const auto c=row[pivot];
      if(!c.is_zero())for(std::size_t j=pivot;j<row.size();++j)row[j]=row[j]-c*basis[j];
    }
    return row;
  };
  const auto insert=[&](const std::vector<Exact>& candidate) {
    auto row=reduced(candidate);
    const auto p=std::find_if(row.begin(),row.end(),[](const Exact& value){return !value.is_zero();});
    if(p==row.end())return false;
    const auto coefficient=*p;for(auto& value:row)value=value/coefficient;
    rank.emplace(p-row.begin(),std::move(row));return true;
  };
  for(const auto& line:family.physical)if(!insert(line.linear))return std::nullopt;
  const auto physical=family.physical.size();
  std::vector<bool> used(family.auxiliary.size(),false),kept(family.auxiliary.size(),false);
  for(const auto& master:masters)for(std::size_t k=0;k<family.auxiliary.size();++k)if(master.at(physical+k))used[k]=true;
  for(std::size_t k=0;k<family.auxiliary.size();++k)if(used[k]) {
    if(!insert(family.auxiliary[k].linear))return std::nullopt;
    kept[k]=true;
  }
  for(std::size_t k=0;k<family.auxiliary.size();++k)if(!used[k])kept[k]=insert(family.auxiliary[k].linear);
  std::size_t next=0;
  for(std::size_t k=0;k<family.auxiliary.size();++k) {
    if(kept[k])continue;
    for(;;++next) {
      if(next==size)return std::nullopt;
      auto slot=merged.space.zero();slot.linear[next]=zero.constant(1);
      if(insert(slot.linear)){family.auxiliary[k]=std::move(slot);++next;break;}
    }
  }
  return family;
}

// One shared closed master vector feeds every operation at a level. A child
// computes precisely that ordered vector at the parent's anchor. Epsilon and
// Taylor demands belong to numerical execution, never to this exact graph.
inline Graph prepare(const feynman::ExampleFamily& family,
    std::vector<ibp::Integral> requested,const Options& options={},
    PrepareLevel provider={},std::function<void(const Event&)> progress={}) {
  if(!options.max_levels||!options.total_timeout_seconds||!options.max_sources_per_level)
    throw std::invalid_argument("recursion requires positive finite budgets");
  if(family.physical_count<2 || family.physical_count>family.momenta.lines.size())
    throw std::invalid_argument("recursion needs at least two physical propagators");
  const auto levels=family.physical_count-1;
  if(levels>options.max_levels || (!options.anchors.empty() && options.anchors.size()!=levels) ||
     (!options.merges.empty() && options.merges.size()!=levels))
    throw std::invalid_argument("recursive merge/anchor count or level budget");
  for(const auto& anchor:options.anchors)if(anchor<=Rational(0)||anchor>=Rational(1))
    throw std::invalid_argument("recursive matching anchors must lie strictly between zero and one");
  // Family names are labels. Special ladders belong in explicit configuration.
  const auto anchors=options.anchors.empty()?std::vector<Rational>(levels,Rational("1/3")):options.anchors;
  if(!provider)provider=[](const auto& b,const auto& d,const auto& f,auto p,const auto& r,const auto& o){
    return level::prepare(b,d,f,p,r,o);
  };
  ExactField field({"x","eps","I"});Exact x(field,"x"),eps(field,"eps");
  auto dimension=x.constant(family.dimension_at_epsilon_zero)-x.constant(2)*eps;
  Graph graph{field,dimension,family.name,{}, {},family};
  auto raw=ibp::quadratic_family(family.momenta,x,family.physical_count);
  const auto began=std::chrono::steady_clock::now();
  if(requested.empty()) {
    ibp::Integral scalar(raw.scalar_products.size(),0);
    std::fill_n(scalar.begin(),family.physical_count,1);requested.push_back(std::move(scalar));
  }
  // Each present line as a combination of the original propagators.
  std::vector<std::vector<Rational>> composition(family.physical_count,std::vector<Rational>(family.physical_count,Rational(0)));
  for(unsigned i=0;i<family.physical_count;++i)composition[i][i]=Rational(1);
  std::optional<bool> euclidean;
  for(unsigned depth=0;depth<levels;++depth) {
    const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
    if(elapsed>=options.total_timeout_seconds)throw std::runtime_error("recursive preparation total time budget exhausted");
    auto [left,right]=options.merges.empty()?std::pair{0u,1u}:options.merges[depth];
    std::size_t candidates=1,counted_masters=0,first_pair_masters=0;
    const unsigned lines=static_cast<unsigned>(ibp::PropagatorBasis(raw).physical_count);
    // Two lines with the same loop momentum merge into an ordinary propagator
    // with a shifted momentum and mass. A pair with different loop momenta
    // makes a new topology: where its master count is small, its reduction
    // has still been far more expensive, so such pairs are not candidates.
    const std::size_t quadratic=family.momenta.loops*(family.momenta.loops+1)/2;
    const auto same_routing=[&](unsigned l,unsigned r) {
      for(std::size_t k=0;k<quadratic && k<raw.physical[l].linear.size();++k)
        if(!(raw.physical[l].linear[k]-raw.physical[r].linear[k]).is_zero())return false;
      return true;
    };
    bool chosen_from_graph=false;
    if(options.merges.empty() && options.count_critical_points && lines>2) {
      std::vector<std::pair<unsigned,unsigned>> pairs;
      for(unsigned l=0;l<lines;++l)for(unsigned r=l+1;r<lines;++r)if((l==0 && r==1) || same_routing(l,r))pairs.emplace_back(l,r);
      std::vector<std::optional<std::size_t>> counts(pairs.size());
      const auto say=[&](const std::string& text){if(options.selection_progress)options.selection_progress("Level "+std::to_string(depth+1)+": "+text);};
      if(pairs.size()>1) {
        const auto since=[](std::chrono::steady_clock::time_point t){return std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();};
        const auto started=std::chrono::steady_clock::now();
        bool gave_up=false;
        try{counts[0]=options.count_critical_points(composition,pairs[0].first,pairs[0].second,options.selection_seconds,gave_up);}catch(const std::exception&){}
        const double first=since(started);
        // A count that fails for another reason is not a matter of time: as
        // before, the eliminations below are then asked.
        if(!counts[0] && gave_up) {
          say("the masters of the first merge were not counted within "+std::to_string(options.selection_seconds)+" s ("+std::to_string(lines)+
            " propagators); merging the first pair without comparing the other "+std::to_string(pairs.size()-1));
          // Counting by eliminations costs far more than what just failed.
          chosen_from_graph=true;
        } else if(counts[0]) {
          const unsigned each=options.selection_seconds+static_cast<unsigned>(std::min(4*first,86400.0));
          const unsigned threads=std::min<unsigned>(options.selection_threads?options.selection_threads:std::max(1u,std::thread::hardware_concurrency()),static_cast<unsigned>(pairs.size()-1));
          say("counting the masters of "+std::to_string(pairs.size()-1)+" other merges on "+std::to_string(threads)+" threads, at most "+std::to_string(each)+
            " s each (the first pair: "+std::to_string(*counts[0])+" masters)");
          std::atomic<std::size_t> next{1};
          const auto work=[&]{for(std::size_t i;(i=next++)<pairs.size();)
            try{bool late=false;counts[i]=options.count_critical_points(composition,pairs[i].first,pairs[i].second,each,late);}catch(const std::exception&){}};
          {
            std::vector<std::jthread> workers;
            for(unsigned t=1;t<threads;++t)workers.emplace_back(work);
            work();
          }
          std::size_t answered=0;for(const auto& count:counts)answered+=count.has_value();
          say(std::to_string(answered)+" of "+std::to_string(pairs.size())+" merges counted in "+std::to_string(since(started))+" s");
        }
      }
      // The first pair is kept unless another has strictly fewer masters.
      if(pairs.size()>1 && counts.front()) {
        chosen_from_graph=false;
        std::size_t best=0;candidates=0;
        for(std::size_t i=0;i<pairs.size();++i)if(counts[i]){++candidates;if(*counts[i]<*counts[best])best=i;}
        left=pairs[best].first;right=pairs[best].second;counted_masters=*counts[best];first_pair_masters=*counts.front();
        chosen_from_graph=true;
      }
    }
    if(!chosen_from_graph && options.merges.empty() && options.count_masters && lines>2) {
      // (masters, not the first pair, sources, pair): the first pair is kept
      // unless another has strictly fewer masters.
      std::optional<std::tuple<std::size_t,bool,std::size_t,unsigned,unsigned>> best;
      candidates=0;
      // The first pair is counted first and sets the time allowed for each
      // other pair, its own plus a second: a family that takes longer to
      // count is a larger one.
      auto allowance=options.reduction;
      // Candidates are the first pair and the pairs with the same loop
      // momentum, whose merged line is again an ordinary propagator with a
      // shifted momentum and mass. A pair with different loop momenta makes a
      // new topology; where counted, such pairs closed on larger levels and
      // cost the most to count.
      for(unsigned l=0;l<lines;++l)for(unsigned r=l+1;r<lines;++r) {
        if(!(l==0 && r==1) && !same_routing(l,r))continue;
        if(l==0 && r==1) {
          bool alternatives=false;
          for(unsigned a=0;a<lines && !alternatives;++a)for(unsigned b=a+1;b<lines;++b)if(!(a==0 && b==1) && same_routing(a,b)){alternatives=true;break;}
          if(!alternatives){candidates=1;l=lines;break;}
        }
        const auto candidate_start=std::chrono::steady_clock::now();
        if(std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()>=options.total_timeout_seconds)break;
        std::optional<std::size_t> masters;std::size_t source_count=0;
        try {
          const ibp::PropagatorBasis incoming(raw),merged(ibp::merge(raw,l,r,x));
          if(merged.physical_count==1){candidates=1;best.reset();l=lines;break;}
          std::set<ibp::Integral,ibp::IntegralOrder> candidate_sources;
          for(const auto& target:requested)
            for(const auto& [a,c]:pullback::plan(incoming,merged,l,r,target,x,options.numerator_expansion).source_row)if(!c.is_zero())candidate_sources.insert(a);
          if(candidate_sources.empty() || candidate_sources.size()>options.max_sources_per_level)continue;
          source_count=candidate_sources.size();
          masters=options.count_masters(merged,dimension,field,0,{candidate_sources.begin(),candidate_sources.end()},allowance);
        } catch(const std::exception&) {}
        if(l==0 && r==1) {
          // A small level is not worth the search: counting the alternatives
          // would cost more than transporting it.
          if(!masters || *masters<options.merge_selection_threshold){candidates=1;best.reset();l=lines;break;}
          const double first=std::chrono::duration<double>(std::chrono::steady_clock::now()-candidate_start).count();
          allowance.provider.timeout_seconds=std::min(allowance.provider.timeout_seconds,static_cast<unsigned>(first)+1u);
        }
        if(!masters)continue;
        ++candidates;
        if(l==0 && r==1)first_pair_masters=*masters;
        const std::tuple<std::size_t,bool,std::size_t,unsigned,unsigned> score{*masters,!(l==0 && r==1),source_count,l,r};
        if(!best || score<*best)best=score;
      }
      if(best){left=std::get<3>(*best);right=std::get<4>(*best);counted_masters=std::get<0>(*best);}
      candidates=std::max<std::size_t>(candidates,1);
    }
    if(options.merge_chosen && candidates>1)options.merge_chosen(depth+1,left,right,candidates,counted_masters,first_pair_masters);
    auto anchor=anchors[depth];
    auto merged_raw=ibp::merge(raw,left,right,x);
    Node node{ibp::PropagatorBasis(raw),ibp::PropagatorBasis(std::move(merged_raw)),
      left,right,anchor,std::move(requested),{}, {},{}, {},false};
    std::set<ibp::Integral,ibp::IntegralOrder> sources;
    for(const auto& target:node.requested) {
      auto operation=pullback::plan(node.incoming,node.merged,left,right,target,x,options.numerator_expansion);
      for(const auto& [a,c]:operation.source_row)if(!c.is_zero())sources.insert(a);
      if(sources.size()>options.max_sources_per_level)throw std::runtime_error("recursive source union budget exhausted");
      node.operations.push_back(std::move(operation));
    }
    node.source_indices.assign(sources.begin(),sources.end());
    node.scalar_leaf=node.merged.physical_count==1;
    if(!node.scalar_leaf && !node.source_indices.empty()) {
      auto budget=options.reduction;
      budget.total_timeout_seconds=std::min(budget.total_timeout_seconds,
        std::max(1u,static_cast<unsigned>(options.total_timeout_seconds-elapsed)));
      node.closure=provider(node.merged,dimension,field,0,node.source_indices,budget);
      if(!node.closure.success)throw std::runtime_error("recursive level "+std::to_string(depth+1)+": "+node.closure.reason);
      if(options.choose_masters) {
        const auto masters=node.closure.ordered_basis.size();
        node.closure=options.choose_masters(std::move(node.closure),node.merged,field,node.source_indices,depth,
          [&](const std::vector<ibp::Integral>& requests,double seconds) {
            auto limited=budget;
            if(seconds<static_cast<double>(limited.total_timeout_seconds))limited.total_timeout_seconds=std::max(1u,static_cast<unsigned>(seconds));
            limited.provider.timeout_seconds=std::min(limited.provider.timeout_seconds,limited.total_timeout_seconds);
            return provider(node.merged,dimension,field,0,requests,limited);
          });
        if(!node.closure.success || node.closure.ordered_basis.size()!=masters || node.closure.matrix.size()!=masters)
          throw std::logic_error("recursive master choice changed the size of a level");
      }
      if(options.exchange_masters) {
        const auto masters=node.closure.ordered_basis.size();
        node.closure=options.exchange_masters(std::move(node.closure),depth);
        if(!node.closure.success || node.closure.ordered_basis.size()!=masters || node.closure.matrix.size()!=masters)
          throw std::logic_error("recursive master exchange changed the size of a level");
      }
      if(node.closure.target_rows.size()!=node.source_indices.size())throw std::logic_error("recursive closure omitted source coordinates");
      std::map<ibp::Integral,std::size_t,ibp::IntegralOrder> source_position;
      for(std::size_t i=0;i<node.source_indices.size();++i)source_position.emplace(node.source_indices[i],i);
      for(const auto& operation:node.operations) {
        std::vector<Exact> row(node.closure.ordered_basis.size(),x.constant(0));
        for(const auto& [a,c]:operation.source_row) {
          const auto& source=node.closure.target_rows[source_position.at(a)];
          if(source.size()!=row.size())throw std::logic_error("recursive source-coordinate shape mismatch");
          for(std::size_t j=0;j<row.size();++j)row[j]=row[j]+c*source[j];
        }
        node.observable_rows.push_back(std::move(row));
      }
      // An endpoint is taken where the integrals are analytic there (the
      // kinematic rule), the level's numerator slots stay independent at it
      // (the masters keep their meaning in the pinched family) and the
      // level's equation confirms that the endpoint values fix the solution.
      // A level made of exact identities alone has no equation to solve.
      std::optional<ibp::QuadraticFamily> child;
      const bool transported=std::any_of(node.operations.begin(),node.operations.end(),
        [](const auto& operation){return operation.operation!=feynman::Operation::Direct;});
      if(options.boundary!=BoundaryPolicy::anchor && options.anchors.empty() && transported && !node.closure.ordered_basis.empty()) {
        if(!euclidean)try{euclidean=causal::euclidean(family);}catch(const std::exception&){euclidean=false;}
        for(const bool upper:{false,true}) {
          if(!endpoint_boundary::kinematics(family,*euclidean,composition,left,right,upper?left:right).qualifies)continue;
          std::string refusal;
          auto pinched=endpoint_family(node.merged,0,Rational(upper?1:0),node.closure.ordered_basis);
          if(!pinched)refusal="a numerator of the level's masters depends on the propagators at the endpoint";
          else if(left>=right || !endpoint_mass_bound(*pinched,left))
            refusal="F of the pinched family is not bounded below by the kept line's parameter times U (its mass is cancelled by on-shell external momenta)";
          else {
            const auto verdict=endpoint_boundary::examine(node.closure.matrix,0,1,upper,options.endpoint);
            if(verdict.qualifies) {
              node.boundary=upper?BoundarySite::upper:BoundarySite::lower;node.analytic_solutions=verdict.analytic_dimension;
              child=std::move(pinched);break;
            }
            refusal=verdict.reason;
          }
          if(options.endpoint_refused)options.endpoint_refused(depth+1,upper,refusal);
        }
      }
      if(node.boundary==BoundarySite::anchor)anchor=select_regular_anchor(node.closure.matrix,node.observable_rows,0,1,anchor,options.anchors);
      else anchor=Rational(node.boundary==BoundarySite::upper?1:0);
      node.anchor=anchor;
      requested=node.closure.ordered_basis;
      // This is the only numerical anchor substitution in exact preparation.
      raw=child?std::move(*child):at_anchor(node.merged,0,anchor);
      for(std::size_t i=0;i<composition[left].size();++i)
        composition[left][i]=anchor*composition[left][i]+(Rational(1)-anchor)*composition[right][i];
      composition.erase(composition.begin()+right);
    }
    const auto total=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
    Event event{depth+1,static_cast<unsigned>(node.merged.physical_count),node.requested.size(),
      node.source_indices.size(),node.closure.ordered_basis.size(),node.closure.passes,total,node.scalar_leaf,anchors[depth],anchor,
      left,right,candidates,counted_masters,first_pair_masters,node.boundary};
    const bool done=node.scalar_leaf || node.closure.ordered_basis.empty();
    graph.nodes.push_back(std::move(node));graph.events.push_back(event);if(progress)progress(event);
    // The total budget gates the START of a level (top of this loop). A level
    // that has just completed is never discarded for finishing late.
    if(done)break;
  }
  if(options.boundary==BoundaryPolicy::endpoint && options.anchors.empty() &&
      std::none_of(graph.nodes.begin(),graph.nodes.end(),[](const Node& node){return node.boundary!=BoundarySite::anchor;}))
    throw std::domain_error("an endpoint boundary was requested, but no level of this family qualifies for one");
  return graph;
}
} // namespace diffexp::recursion
