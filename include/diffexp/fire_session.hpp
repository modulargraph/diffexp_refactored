#pragma once
#include "diffexp/fire_checkpoint.hpp"

namespace diffexp::fire {
// FIRE #storage completion flags describe a fixed demand, so replaying them for
// new integrals silently skips needed sectors. A session instead reloads exact
// solved #rules into fresh demand bookkeeping (parser.cpp add_rules). Every
// rule is a completed reduction from this session's immutable scientific input.
class Session {
 public:
  Session(ibp::PropagatorBasis basis,Exact dimension,ExactField field,const Options& options,std::filesystem::path batch_cache_directory={},std::filesystem::path pending_cache_directory={})
      : basis_(std::move(basis)),dimension_(std::move(dimension)),field_(std::move(field)),
        semantic_(semantic(options)) {if(!batch_cache_directory.empty())batch_cache_=std::make_unique<fire_batch::Cache>(std::move(batch_cache_directory));if(!pending_cache_directory.empty())pending_cache_=std::make_unique<fire_checkpoint::Store>(std::move(pending_cache_directory));}
  Session(const Session&)=delete;
  Session& operator=(const Session&)=delete;
  Result operator()(const std::vector<ibp::Integral>& requests,const Options& limits) {
    Result failure;
    try {
      if(!limits.timeout_seconds||!limits.memory_bytes)throw std::invalid_argument("FIRE session requires positive resource budgets");
      if(!limits.threads||limits.threads>64||!limits.simplifier_threads||limits.simplifier_threads>64||limits.threads*limits.simplifier_threads>256)throw std::invalid_argument("FIRE session worker counts outside bounded range");
      if(semantic(limits)!=semantic_)throw std::invalid_argument("FIRE session scientific context or executable changed");
      std::set<ibp::Integral> next(requests.begin(),requests.end());
      if(!std::includes(next.begin(),next.end(),demand_.begin(),demand_.end()))throw std::invalid_argument("FIRE session demands must be monotone supersets");
      if(complete_&&next==demand_){++cache_hits_;return last_;}
      if(batch_cache_)if(auto cached=batch_cache_->lookup(basis_,dimension_,field_,requests,limits)){demand_=std::move(next);last_=*cached;complete_=true;++cache_hits_;++durable_hits_;return *cached;}
      const auto* prior=complete_?&last_.reductions:nullptr;
      std::optional<boost::json::object> pending_identity;std::filesystem::path resume;
      if(pending_cache_){pending_identity=fire_checkpoint::identity(basis_,dimension_,field_,requests,limits,prior);if(auto found=pending_cache_->lookup(*pending_identity)){resume=*found;++checkpoint_hits_;}}
      auto result=pending_cache_?fire_checkpoint::reduce(basis_,dimension_,field_,requests,limits,prior,resume):fire::reduce(basis_,dimension_,field_,requests,limits,prior);++runs_;
      if(result.success){
        if(batch_cache_)batch_cache_->put(basis_,dimension_,field_,requests,limits,result);
        // Retire pending algebra only after its exact completed batch is durable.
        if(pending_cache_&&batch_cache_)pending_cache_->retire(*pending_identity);
        demand_=std::move(next);last_=result;complete_=true;
      }else if(pending_cache_&&!result.directory.empty()&&std::filesystem::exists(result.directory/"checkpoint.json"))pending_cache_->put(result.directory);
      return result;
    }catch(const std::exception& e){failure.reason=e.what();return failure;}
  }
  std::size_t runs()const{return runs_;}
  std::size_t cache_hits()const{return cache_hits_;}
  std::size_t durable_hits()const{return durable_hits_;}
  std::size_t checkpoint_hits()const{return checkpoint_hits_;}
 private:
  ibp::PropagatorBasis basis_;
  Exact dimension_;
  ExactField field_;
  std::string semantic_;
  std::set<ibp::Integral> demand_;
  Result last_;
  bool complete_=false;
  std::size_t runs_=0,cache_hits_=0,durable_hits_=0,checkpoint_hits_=0;
  std::unique_ptr<fire_batch::Cache> batch_cache_;
  std::unique_ptr<fire_checkpoint::Store> pending_cache_;
  std::string semantic(const Options& options)const {
    std::ostringstream out;
    out<<start(basis_,dimension_,options.zero_sectors)<<"\nphysical="<<basis_.physical_count<<"\n";
    out<<"loops="<<basis_.space.loops<<"\n";
    for(const auto& row:basis_.space.external_gram){for(const auto& x:row)out<<x.str()<<";";out<<"\n";}
    for(const auto& a:basis_.denominators){out<<a.constant.str()<<";";for(const auto& x:a.linear)out<<x.str()<<";";out<<"\n";}
    for(const auto& name:field_.variables())out<<name.size()<<":"<<name<<"\n";
    if(options.executable.empty()||access(options.executable.c_str(),X_OK)){out<<"executable-unavailable";return out.str();}
    auto binary=std::filesystem::canonical(options.executable);out<<binary.string()<<"\n"<<std::filesystem::file_size(binary)<<"\n";auto ticks=std::filesystem::last_write_time(binary).time_since_epoch().count();out<<static_cast<long long>(ticks/1000000000)<<":"<<static_cast<long long>(ticks%1000000000);
    return out.str();
  }
};
} // namespace diffexp::fire
