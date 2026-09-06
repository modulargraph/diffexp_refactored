#pragma once
#include "diffexp/fire_modular.hpp"
#include <ibp/solver.hpp>

namespace diffexp::ibp_solver {
// The pinned upstream core works in the exact denominator coordinates owned by
// DiffExp. No momentum reconstruction, renamed family, or new ISP completion.
struct Options {
  std::filesystem::path cache_directory;
  unsigned dots=1,numerators=2,max_degree=16,max_primes=8;
  std::size_t max_samples_per_prime=512;
  std::function<void(const std::string&)> progress;
};
struct Statistics {
  std::size_t probes=0,equations=0;
  double generation_seconds=0,elimination_seconds=0;
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
  }
  const Statistics& statistics()const{return stats_;}
  fire::Result operator()(const std::vector<ibp::Integral>& requested,const std::vector<modular::Word>& point,
      modular::Word modulus,const fire::Options& limits) {
    fire::Result out;
    try {
      if(!limits.timeout_seconds||!limits.memory_bytes||requested.empty()||requested.size()>10000)
        throw std::invalid_argument("IBP solver requires bounded nonempty requests");
      const auto begin=std::chrono::steady_clock::now();
      auto remaining=[&]{auto s=limits.timeout_seconds-std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        if(s<=0)throw std::runtime_error("IBP solver probe time budget exceeded");return s;};
      ::ibp::Field field(modulus);::ibp::Geometry geometry{basis_.space.loops,basis_.space.externals(),static_cast<unsigned>(basis_.physical_count),static_cast<unsigned>(basis_.space.size())};
      geometry.trace=trace_;geometry.zero_sectors=zero_sectors_;
      for(const auto& a:limits.zero_sectors){fire_batch::validate_index(a,basis_);geometry.zero_sectors[sector(a)]=true;}
      auto affine=[&](const ibp::Affine& a){::ibp::Affine r;for(const auto& c:a.linear)r.push_back(modular::evaluate(c,point,modulus));r.push_back(modular::evaluate(a.constant,point,modulus));return r;};
      for(const auto& row:contractions_){geometry.contractions.emplace_back();for(const auto& a:row)geometry.contractions.back().push_back(affine(a));}
      ::ibp::SeedOptions seeds;seeds.dots=options_.dots;seeds.numerators=options_.numerators;
      std::vector<::ibp::Integral> targets;
      for(const auto& a:requested){fire_batch::validate_index(a,basis_);::ibp::Integral target;unsigned dots=0,nums=0;
        for(unsigned i=0;i<a.size();++i){target.powers[i]=a[i];if(a[i]>0)dots+=a[i]-1;else nums-=a[i];}
        // A seed with one fewer dot can produce the requested raised integral.
        seeds.dots=std::max(seeds.dots,dots?dots-1:0);seeds.numerators=std::max(seeds.numerators,nums);targets.push_back(target);
      }
      seeds.max_terms=std::min<std::size_t>(10000000,limits.memory_bytes/64);seeds.seconds=remaining();
      if(!seeds.max_terms)throw std::invalid_argument("IBP solver memory budget too small");
      auto program=::ibp::generate(geometry,field,seeds,std::move(targets));
      const auto generated=std::chrono::steady_clock::now();
      ::ibp::SolveOptions solve;solve.seconds=remaining();solve.max_fill=std::min<std::size_t>(12000000,limits.memory_bytes/64);
      ::ibp::Solver solver(program,field,solve);solver.solve(modular::evaluate(dimension_,point,modulus));auto rows=solver.targets();remaining();
      auto indices=[&](const ::ibp::Integral& a){ibp::Integral r;for(unsigned i=0;i<geometry.n;++i)r.push_back(a.powers[i]);return r;};
      std::set<ibp::Integral> terminal;
      for(unsigned i=0;i<rows.size();++i){ibp::Relation row;for(const auto& term:rows[i]){auto a=indices(program.integrals[term.column]);terminal.insert(a);ibp::add(row,a,dimension_.parse(std::to_string(term.value)));}out.reductions.emplace(requested[i],std::move(row));}
      for(const auto& a:terminal){auto [it,inserted]=out.reductions.emplace(a,ibp::Relation{{a,dimension_.constant(1)}});
        if(!inserted&&it->second!=ibp::Relation{{a,dimension_.constant(1)}})throw std::runtime_error("IBP solver returned a nonterminal output basis");}
      out.success=true;fire_batch::validate(out,basis_,dimension_,requested);
      ++stats_.probes;stats_.equations+=program.equations.size();stats_.generation_seconds+=std::chrono::duration<double>(generated-begin).count();
      stats_.elimination_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-generated).count();
    }catch(const std::exception& e){out.success=false;out.reason=e.what();}
    return out;
  }
 private:
  std::size_t sector(const ibp::Integral& a)const {std::size_t mask=0;for(unsigned i=0;i<basis_.physical_count;++i)if(a[i]>0)mask|=std::size_t(1)<<i;return mask;}
  ibp::PropagatorBasis basis_;Exact dimension_;Options options_;Statistics stats_;
  std::vector<std::vector<ibp::Affine>> contractions_;std::vector<bool> trace_,zero_sectors_;
};
class Session {
 public:
  Session(const ibp::PropagatorBasis& basis,const Exact& dimension,const ExactField& field,const Options& options)
    :sampler_(std::make_shared<Sampler>(basis,dimension,options)),session_(basis,dimension,field,configuration(options,sampler_)){}
  fire::Result operator()(const std::vector<ibp::Integral>& requested,const fire::Options& limits){return session_(requested,limits);}
  const Statistics& statistics()const{return sampler_->statistics();}
 private:
  static fire_modular::Options configuration(const Options& o,const std::shared_ptr<Sampler>& sampler){
    fire_modular::Options out;out.cache_directory=o.cache_directory;out.max_degree=o.max_degree;out.max_primes=o.max_primes;
    out.max_samples_per_prime=o.max_samples_per_prime;out.progress=o.progress;out.modulus=prime;
    out.provider_identity="ibp-solver-574c5bd-adapter-v1-dots-"+std::to_string(o.dots)+"-numerators-"+std::to_string(o.numerators);
    out.prime_table="descending-61-bit-primes-v1";
    out.sample_provider=[sampler](const auto& r,const auto& p,auto m,const auto& limits){return (*sampler)(r,p,m,limits);};return out;
  }
  std::shared_ptr<Sampler> sampler_;fire_modular::Session session_;
};
} // namespace diffexp::ibp_solver
