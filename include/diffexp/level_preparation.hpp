#pragma once
#include "diffexp/fire_session.hpp"

namespace diffexp::level {
struct Options {
  fire::Options provider;
  std::filesystem::path batch_cache_directory,pending_cache_directory;
  std::size_t max_passes=8,max_masters=256,max_demands=2000,max_equations=20000;
  unsigned total_timeout_seconds=180;
};
struct Result {
  bool success=false;
  std::string reason;
  std::vector<ibp::Integral> ordered_basis,requested_integrals;
  std::vector<std::vector<Exact>> matrix,target_rows;
  std::vector<std::filesystem::path> fire_directories;
  std::vector<ibp::Relation> source_identities;
  std::vector<ibp::ExactReducer::Witness> differential_witnesses,target_witnesses;
  std::size_t passes=0,demands=0,equations=0;
  double elapsed_seconds=0;
};
using Provider=std::function<fire::Result(const std::vector<ibp::Integral>&,const fire::Options&)>;

// This is a closure proof for the requested observables, conditional on the
// imported exact IBP identities. It makes no claim about a global master count.
// Demands and source identities only grow; each pass sends one shared batch.
inline Result prepare(const ibp::PropagatorBasis& basis,const Exact& dimension,
    const ExactField& field,std::size_t parameter,std::vector<ibp::Integral> requested,
    const Options& options,Provider provider={}) {
  Result result;result.requested_integrals=requested;const auto began=std::chrono::steady_clock::now();
  try {
    if(requested.empty()||!options.max_passes||!options.max_masters||!options.max_demands||
       !options.max_equations||!options.total_timeout_seconds)throw std::invalid_argument("level preparation requires positive finite budgets and requests");
    if(parameter>=field.variables().size())throw std::invalid_argument("level parameter outside exact field");
    (void)(dimension+Exact(field));ibp::Generator generator(basis,dimension);
    if(!provider) {
      auto session=std::make_shared<fire::Session>(basis,dimension,field,options.provider,options.batch_cache_directory,options.pending_cache_directory);
      provider=[session](const auto& batch,const auto& opts){return (*session)(batch,opts);};
    }
    auto one=dimension.constant(1);
    std::set<ibp::Integral,ibp::IntegralOrder> demand(requested.begin(),requested.end()),discovered;
    ibp::ExactReducer identities(dimension,options.max_equations);
    for(std::size_t pass=0;pass<options.max_passes;++pass) {
      const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
      if(elapsed>=options.total_timeout_seconds)throw std::runtime_error("level preparation total time budget exhausted");
      if(demand.size()>options.max_demands)throw std::runtime_error("level preparation demand budget exhausted");
      auto provider_options=options.provider;
      provider_options.timeout_seconds=std::min(provider_options.timeout_seconds,
        std::max(1u,static_cast<unsigned>(options.total_timeout_seconds-elapsed)));
      std::vector<ibp::Integral> batch(demand.begin(),demand.end());
      auto reduced=provider(batch,provider_options);++result.passes;result.demands=demand.size();
      if(!reduced.directory.empty())result.fire_directories.push_back(reduced.directory);
      if(!reduced.success)throw std::runtime_error("level FIRE pass failed: "+reduced.reason);
      for(const auto& a:batch)if(!reduced.reductions.count(a))throw std::runtime_error("level provider omitted a demanded integral");
      for(const auto& [a,row]:reduced.reductions) {
        // A reduction is a complete relation, including the zero case. Identity
        // master rows contribute no equation but remain candidate coordinates.
        ibp::Relation relation{{a,one}};ibp::add_scaled(relation,row,-one);
        if(!relation.empty()){result.source_identities.push_back(relation);identities.insert(std::move(relation));}
        for(const auto& [master,coefficient]:row)if(!coefficient.is_zero())discovered.insert(master);
      }
      result.equations=identities.equation_count();
      if(discovered.size()>options.max_masters)throw std::runtime_error("level preparation candidate-master budget exhausted: "+std::to_string(discovered.size())+" candidates exceed "+std::to_string(options.max_masters));
      // Drop dependencies revealed by later batches without changing the exact
      // span. Deterministic native integral ordering fixes the coordinate order.
      std::map<ibp::Integral,ibp::Relation,ibp::IntegralOrder> echelon;
      std::vector<ibp::Integral> selected;
      for(const auto& master:discovered) {
        auto row=identities.reduce({{master,one}}).remainder;
        for(auto it=echelon.rbegin();it!=echelon.rend();++it) {
          auto p=row.find(it->first);if(p!=row.end()){auto c=p->second;ibp::add_scaled(row,it->second,-c);}
        }
        if(row.empty())continue;
        auto pivot=std::prev(row.end());auto key=pivot->first;auto c=pivot->second;
        for(auto& [a,v]:row)v=v/c;echelon.emplace(key,std::move(row));selected.push_back(master);
      }
      result.ordered_basis=selected;
      if(selected.empty()) {
        bool all_zero=true;for(const auto& a:requested)all_zero=all_zero&&identities.reduce({{a,one}}).remainder.empty();
        if(!all_zero)throw std::runtime_error("level provider returned an empty unresolved span");
        result.target_rows.assign(requested.size(),{});
        for(const auto& a:requested)result.target_witnesses.push_back(identities.reduce({{a,one}}).witness);
        result.matrix.clear();result.success=true;break;
      }
      const auto before=demand.size();
      std::vector<ibp::Relation> derivatives;
      for(const auto& master:selected) {
        demand.insert(master);auto row=generator.derivative(master,parameter);
        for(const auto& [a,c]:row)if(!c.is_zero())demand.insert(a);
        derivatives.push_back(std::move(row));
      }
      ibp::BasisReduction coordinates(identities,selected,dimension);
      bool closed=true;std::vector<std::vector<Exact>> matrix,targets;
      try {
        for(const auto& row:derivatives)matrix.push_back(coordinates.resolve(row));
        for(const auto& a:requested)targets.push_back(coordinates.resolve({{a,one}}));
      } catch(const std::runtime_error&) {closed=false;}
      if(closed){
        auto certificate=[&](ibp::Relation row,const std::vector<Exact>& coordinates) {
          for(std::size_t j=0;j<selected.size();++j)ibp::add(row,selected[j],-coordinates[j]);
          auto proof=identities.reduce(row);
          if(!proof.remainder.empty()||!identities.verify(row,proof))throw std::logic_error("level closure certificate failed");
          return proof.witness;
        };
        for(std::size_t i=0;i<derivatives.size();++i)result.differential_witnesses.push_back(certificate(derivatives[i],matrix[i]));
        for(std::size_t i=0;i<requested.size();++i)result.target_witnesses.push_back(certificate({{requested[i],one}},targets[i]));
        result.matrix=std::move(matrix);result.target_rows=std::move(targets);result.success=true;break;
      }
      if(demand.size()==before)throw std::runtime_error("level span unresolved without a new demand; provider made no progress");
    }
    if(!result.success)throw std::runtime_error("level preparation pass budget exhausted");
  }catch(const std::exception& e){result.reason=e.what();}
  result.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
  if(result.success&&result.elapsed_seconds>options.total_timeout_seconds){result.success=false;result.reason="level preparation total time budget exhausted";}
  return result;
}
} // namespace diffexp::level
