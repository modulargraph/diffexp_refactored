#pragma once
#include "diffexp/families.hpp"
#include "diffexp/level_preparation.hpp"
#include "diffexp/merge_pullback.hpp"
#include <numeric>
#include <limits>

namespace diffexp::recursion {
struct Options {
  level::Options reduction;
  pullback::Budget numerator_expansion;
  std::vector<Rational> anchors;
  std::vector<std::pair<unsigned,unsigned>> merges;
  unsigned max_levels=16,total_timeout_seconds=600;
  std::size_t max_sources_per_level=5000;
};
struct Event {
  unsigned depth,physical_propagators;
  std::size_t requests,sources,masters,provider_passes;
  double elapsed_seconds;
  bool scalar_leaf;
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
};
struct Graph {
  ExactField field;
  Exact dimension;
  std::string family_name;
  std::vector<Node> nodes;
  std::vector<Event> events;
  std::optional<feynman::ExampleFamily> definition;
};

inline std::vector<Rational> default_anchors(const std::string& family,unsigned levels) {
  // Published Henn ladder data: equal anchors hit the x2=x1 singularity.
  // See Examples/FeynmanTrick/HennDoublePentagonBoundary.wl fixedParameterValues.
  if(family=="henn_double_pentagon_x0") {
    if(levels!=7)throw std::invalid_argument("Henn example requires seven recursive levels");
    return {Rational("1/5"),Rational("3/10"),Rational("2/5"),Rational("1/2"),
      Rational("3/5"),Rational("7/10"),Rational("4/5")};
  }
  return std::vector<Rational>(levels,Rational("1/3"));
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

inline Rational select_regular_anchor(const std::vector<std::vector<Exact>>& connection,
    const std::vector<std::vector<Exact>>& observables,std::size_t parameter,std::size_t epsilon,
    const Rational& preferred,bool allow_alternative) {
  const auto regular=[&](const Rational& at){return regular_anchor(connection,parameter,epsilon,at)&&regular_anchor(observables,parameter,epsilon,at);};
  if(regular(preferred))return preferred;
  if(!allow_alternative)throw std::domain_error("requested recursive anchor "+preferred.str()+" is singular in a connection or observable epsilon coefficient");
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
  for(unsigned depth=0;depth<levels;++depth) {
    const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
    if(elapsed>=options.total_timeout_seconds)throw std::runtime_error("recursive preparation total time budget exhausted");
    auto [left,right]=options.merges.empty()?std::pair{0u,1u}:options.merges[depth];
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
      anchor=select_regular_anchor(node.closure.matrix,node.observable_rows,0,1,anchor,options.anchors.empty());
      node.anchor=anchor;
      requested=node.closure.ordered_basis;
      // This is the only numerical anchor substitution in exact preparation.
      raw=at_anchor(node.merged,0,anchor);
    }
    const auto total=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
    Event event{depth+1,static_cast<unsigned>(node.merged.physical_count),node.requested.size(),
      node.source_indices.size(),node.closure.ordered_basis.size(),node.closure.passes,total,node.scalar_leaf};
    const bool done=node.scalar_leaf || node.closure.ordered_basis.empty();
    graph.nodes.push_back(std::move(node));graph.events.push_back(event);if(progress)progress(event);
    if(total>options.total_timeout_seconds)throw std::runtime_error("recursive preparation total time budget exhausted");
    if(done)break;
  }
  return graph;
}
} // namespace diffexp::recursion
