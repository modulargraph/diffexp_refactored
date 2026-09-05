#pragma once

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace diffexp {
// A resource demand is not part of mathematical identity. An epsilon extension
// changes a guarantee, never the equation, branch, master ordering or parents.
struct Demand {
  int epsilon_low = 0;
  int epsilon_high = 0;
  unsigned taylor_order = 0;
  unsigned precision_bits = 64;
  unsigned publication_digits = 0;
  bool dominates(const Demand& b) const {
    return epsilon_low <= b.epsilon_low && epsilon_high >= b.epsilon_high &&
           taylor_order >= b.taylor_order && precision_bits >= b.precision_bits &&
           publication_digits >= b.publication_digits;
  }
  Demand joined(const Demand& b) const {
    return {std::min(epsilon_low,b.epsilon_low),std::max(epsilon_high,b.epsilon_high),
            std::max(taylor_order,b.taylor_order),std::max(precision_bits,b.precision_bits),
            std::max(publication_digits,b.publication_digits)};
  }
  void validate() const {
    if (epsilon_low > epsilon_high || precision_bits < 64)
      throw std::invalid_argument("invalid resource demand");
  }
  friend bool operator==(const Demand&, const Demand&) = default;
};

struct Refinement : std::runtime_error {
  std::string artifact;
  Demand needed;
  Refinement(std::string id, Demand demand, const std::string& reason)
      : std::runtime_error(reason), artifact(std::move(id)), needed(demand) {}
};

template<class Value> class Executor {
 public:
  struct Artifact {
    Demand guarantee;
    Value value;
    std::string semantic_identity;
    std::vector<std::shared_ptr<const Artifact>> parents;
  };
  using Ref = std::shared_ptr<const Artifact>;
  struct Dependency {
    std::string id;
    std::function<Demand(const Demand&)> demand;
  };
  struct Node {
    std::string id;
    // Exact equation/basis/branch/geometry identity, excluding resource knobs.
    std::string identity;
    std::vector<Dependency> parents;
    std::function<Artifact(const Demand&,const std::vector<Ref>&)> build;
  };
  struct Event {
    std::string kind, artifact, reason;
    Demand previous, required;
    std::vector<std::string> invalidated;
  };

  void add(Node node) {
    if (node.id.empty() || node.identity.empty() || !node.build || nodes_.contains(node.id))
      throw std::invalid_argument("invalid or duplicate operation node");
    std::set<std::string> ids;
    for (const auto& p : node.parents)
      if (!p.demand || !ids.insert(p.id).second)
        throw std::invalid_argument("invalid or duplicate dependency");
    nodes_.emplace(node.id, std::move(node));
  }
  Ref run(const std::string& target, Demand demand, unsigned max_refinements=16) {
    demand.validate();
    validate_graph(target);
    for (unsigned attempt=0;;++attempt) {
      try { return obtain(target,demand); }
      catch (const Refinement& r) {
        r.needed.validate();
        if (attempt >= max_refinements) throw std::runtime_error("local refinement budget exhausted: "+r.artifact);
        auto ancestors=closure(target);
        if (!ancestors.contains(r.artifact))
          throw std::logic_error("refinement targets an unrelated artifact: "+r.artifact);
        Demand previous = last_demand_.at(r.artifact);
        Demand next=previous.joined(r.needed);
        if (next==previous)
          throw std::logic_error("non-progressing refinement rejected: "+r.artifact);
        floors_.insert_or_assign(r.artifact,next);
        auto affected=descendants(r.artifact);
        // Invalidation removes eligibility in this execution, not the old
        // immutable artifact. Unchanged siblings and ancestors stay reusable.
        for (const auto& id : affected) cache_.erase(id);
        events_.push_back({"refine",r.artifact,r.what(),previous,next,affected});
      }
    }
  }
  const std::vector<Event>& events() const { return events_; }
  unsigned builds(const std::string& id) const {
    auto p=build_counts_.find(id); return p==build_counts_.end()?0:p->second;
  }
 private:
  std::map<std::string,Node> nodes_;
  std::map<std::string,std::vector<Ref>> cache_;
  std::map<std::string,Demand> floors_,last_demand_;
  std::map<std::string,unsigned> build_counts_;
  std::vector<Event> events_;
  void validate_graph(const std::string& target) const {
    std::set<std::string> visiting,done;
    std::function<void(const std::string&)> walk=[&](const std::string& id) {
      if (done.contains(id)) return;
      if (!nodes_.contains(id)) throw std::invalid_argument("missing dependency: "+id);
      if (!visiting.insert(id).second) throw std::invalid_argument("operation graph has a cycle");
      for (const auto& p : nodes_.at(id).parents) walk(p.id);
      visiting.erase(id); done.insert(id);
    };
    walk(target);
  }
  std::set<std::string> closure(const std::string& id) const {
    std::set<std::string> out{id};
    for (const auto& p : nodes_.at(id).parents) {
      auto parents=closure(p.id); out.insert(parents.begin(),parents.end());
    }
    return out;
  }
  std::vector<std::string> descendants(const std::string& id) const {
    std::set<std::string> affected{id};
    bool changed=true;
    while (changed) {
      changed=false;
      for (const auto& [name,n] : nodes_) {
        if (affected.contains(name)) continue;
        if (std::any_of(n.parents.begin(),n.parents.end(),[&](const auto& p) {return affected.contains(p.id);}))
          changed=affected.insert(name).second || changed;
      }
    }
    affected.erase(id); return {affected.begin(),affected.end()};
  }
  Ref obtain(const std::string& id,Demand demand) {
    if (auto f=floors_.find(id); f!=floors_.end()) demand=demand.joined(f->second);
    demand.validate(); last_demand_.insert_or_assign(id,demand);
    const auto& node=nodes_.at(id);
    // Checking the semantic node before visiting parents is essential: a
    // receiving-basis extension must not demand wider upstream boundaries.
    for (const auto& a : cache_[id])
      if (a->semantic_identity==node.identity && a->guarantee.dominates(demand)) {
        events_.push_back({"reuse",id,"guarantee dominates demand",a->guarantee,demand,{}});
        return a;
      }
    std::vector<Ref> parents;
    for (const auto& p : node.parents) parents.push_back(obtain(p.id,p.demand(demand)));
    ++build_counts_[id];
    auto value=node.build(demand,parents);
    value.guarantee.validate();
    if (!value.guarantee.dominates(demand))
      throw std::logic_error("producer did not satisfy its declared demand: "+id);
    if (value.semantic_identity!=node.identity)
      throw std::logic_error("producer changed its semantic identity: "+id);
    value.parents=std::move(parents);
    Ref out=std::make_shared<const Artifact>(std::move(value));
    cache_[id].push_back(out);
    events_.push_back({"build",id,"contract satisfied",{},demand,{}});
    return out;
  }
};
}  // namespace diffexp
