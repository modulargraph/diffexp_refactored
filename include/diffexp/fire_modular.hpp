#pragma once
#include "diffexp/fire_batch_cache.hpp"
#include "diffexp/modular_reconstruction.hpp"

namespace diffexp::fire_modular {
namespace json=boost::json;
namespace fs=std::filesystem;
struct Options {
  fs::path executable,cache_directory;
  unsigned max_degree=16,max_primes=8,probe_timeout_seconds=20;
  std::size_t max_samples_per_prime=512;
  std::function<void(const std::string&)> progress;
};
namespace detail {
inline void publish(const fs::path& target,const std::string& bytes) {
  fs::create_directories(target.parent_path());auto name=target.string()+".new.XXXXXX";
  std::vector<char> tmp(name.begin(),name.end());tmp.push_back(0);int fd=mkstemp(tmp.data());
  if(fd<0)throw std::runtime_error("cannot create modular checkpoint");
  try{std::size_t offset=0;while(offset<bytes.size()){auto n=::write(fd,bytes.data()+offset,bytes.size()-offset);
      if(n<0&&errno==EINTR)continue;if(n<=0)throw std::runtime_error("cannot write modular checkpoint");offset+=n;}
    if(fsync(fd))throw std::runtime_error("cannot sync modular checkpoint");auto status=close(fd);fd=-1;
    if(status)throw std::runtime_error("cannot close modular checkpoint");fs::rename(tmp.data(),target);artifacts::detail::sync_directory(target.parent_path());
  }catch(...){if(fd>=0)close(fd);unlink(tmp.data());throw;}
}
inline void save(const fs::path& path,const json::object& payload) {
  auto bytes=artifacts::detail::canonical(payload);
  publish(path,json::serialize(json::object{{"sha256",artifacts::detail::sha256(bytes)},{"payload",payload}}));
}
inline json::object read(const fs::path& path) {
  auto record=json::parse(fire::read_text(path)).as_object();artifacts::detail::keys(record,{"sha256","payload"});
  if(artifacts::detail::sha256(artifacts::detail::canonical(record.at("payload")))!=artifacts::detail::string(record.at("sha256")))throw std::runtime_error("modular checkpoint checksum mismatch");
  return record.at("payload").as_object();
}
inline std::string table(const std::map<ibp::Integral,ibp::Relation>& reductions,const fire::SymbolMap& symbols) {
  std::map<ibp::Integral,std::size_t> ids;for(const auto& [a,row]:reductions){ids.try_emplace(a,ids.size()+1);for(const auto& [b,c]:row)ids.try_emplace(b,ids.size()+1);}
  std::ostringstream out;out<<"{{";bool first=true;
  for(const auto& [a,row]:reductions){if(!first)out<<",";first=false;out<<"{"<<ids.at(a)<<",{";bool term=true;
    for(const auto& [b,c]:row){if(!term)out<<",";term=false;out<<"{"<<ids.at(b)<<",\""<<fire::rename_symbols(c.str(),symbols.forward)<<"\"}";}out<<"}}";}
  out<<"},{";first=true;for(const auto& [a,id]:ids){if(!first)out<<",";first=false;out<<"{"<<id<<",{1,"<<fire::indices(a)<<"}}";}out<<"}}\n";return out.str();
}
inline std::vector<std::size_t> active_variables(const std::string& source,const fire::SymbolMap& symbols) {
  std::set<std::string> tokens;for(std::size_t i=0;i<source.size();) {
    if(std::isalpha(static_cast<unsigned char>(source[i]))||source[i]=='_'){auto begin=i++;while(i<source.size()&&(std::isalnum(static_cast<unsigned char>(source[i]))||source[i]=='_'))++i;tokens.insert(source.substr(begin,i-begin));}else ++i;
  }
  std::vector<std::size_t> out;for(std::size_t i=0;i<symbols.exported.size();++i)if(tokens.contains(symbols.exported[i]))out.push_back(i);return out;
}
using Key=std::pair<ibp::Integral,ibp::Integral>;
struct Probe {std::vector<modular::Word> point;std::map<ibp::Integral,ibp::Relation> rows;};
inline std::set<ibp::Integral> masters(const Probe& probe) {
  std::set<ibp::Integral> out;for(const auto& [a,row]:probe.rows)for(const auto& [b,c]:row)out.insert(b);return out;
}
inline modular::Word coefficient(const Probe& sample,const Key& key,modular::Word prime) {
  const auto& row=sample.rows.at(key.first);auto it=row.find(key.second);return it==row.end()?0:modular::evaluate(it->second,sample.point,prime);
}
} // namespace detail

// Completed symbolic batches remain reusable. An uncached demand uses only
// FIRE7p, native FLINT interpolation and CRT; there is no symbolic fallback.
// The generic level preparer subsequently proves its derivative/target closure
// conditional on these imported identities, exactly as for symbolic FIRE.
class Session {
 public:
  Session(ibp::PropagatorBasis basis,Exact dimension,ExactField field,Options options,fs::path exact_batches={})
    :basis_(std::move(basis)),dimension_(std::move(dimension)),field_(std::move(field)),options_(std::move(options)) {
    if(options_.cache_directory.empty()||options_.executable.empty()||!options_.max_degree||options_.max_degree>32||options_.max_primes<2||options_.max_primes>14||options_.max_samples_per_prime<8||options_.max_samples_per_prime>4096||!options_.probe_timeout_seconds)
      throw std::invalid_argument("invalid finite modular reconstruction budgets or paths");
    if(!exact_batches.empty())exact_batches_=std::make_unique<fire_batch::Cache>(std::move(exact_batches));
  }
  fire::Result operator()(const std::vector<ibp::Integral>& requested,const fire::Options& limits) {
    fire::Result result;auto began=std::chrono::steady_clock::now();
    try {
      if(!limits.timeout_seconds||!limits.memory_bytes)throw std::invalid_argument("modular provider requires finite positive budgets");
      if(exact_batches_)if(auto cached=exact_batches_->lookup(basis_,dimension_,field_,requested,limits))return *cached;
      const auto semantic=fire_batch::identity(basis_,dimension_,field_,requested,limits);
      fire::SymbolMap symbols(field_.variables());auto source=fire::start(basis_,dimension_,limits.zero_sectors,&symbols);
      auto active=detail::active_variables(source,symbols);
      json::array active_json;for(auto i:active)active_json.push_back(i);
      const auto identity=artifacts::detail::sha256(artifacts::detail::canonical(json::object{{"schema","DiffExp3.NativeModularProvider/v1"},{"science",semantic.json_value()},{"active",active_json},{"prime_table","FIRE7-primes-1-through-16"}}));
      result.directory=fs::absolute(options_.cache_directory)/identity;fs::create_directories(result.directory);
      detail::publish(result.directory/"family.start",source);
      auto status=[&](std::string message){if(options_.progress)options_.progress(message);};
      auto remaining=[&](){const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();if(elapsed>=limits.timeout_seconds)throw std::runtime_error("modular preparation time budget exhausted; samples retained");return std::max(1u,static_cast<unsigned>(limits.timeout_seconds-elapsed));};
      auto load_rows=[&](const json::object& payload){return fire::import_table_text(artifacts::detail::string(payload.at("tables")),field_,basis_.denominators.size(),&symbols);};
      auto verify_native=[&](const auto& rows) {
        ibp::Generator generator(basis_,dimension_);std::size_t checked=0;
        for(const auto& [seed,reduction]:rows)for(const auto& equation:generator.relations(seed)) {
          bool applicable=!equation.empty();for(const auto& [a,c]:equation)if(!rows.contains(a))applicable=false;
          if(!applicable)continue;ibp::Relation residual;
          for(const auto& [a,c]:equation)ibp::add_scaled(residual,rows.at(a),c);
          if(!residual.empty())throw std::runtime_error("reconstructed reductions fail an exact native IBP identity");++checked;
        }return checked;
      };
      if(fs::exists(result.directory/"completed.json")) {
        auto saved=detail::read(result.directory/"completed.json");
        if(artifacts::detail::string(saved.at("identity"))!=identity||artifacts::detail::string(saved.at("verification"))!="fresh-two-primes-three-points-and-applicable-native-IBPs")throw std::runtime_error("modular completion identity or verification mismatch");
        result.reductions=load_rows(saved);result.success=true;fire_batch::validate(result,basis_,dimension_,requested);
        const auto used=saved.at("reconstruction_primes").to_number<unsigned>();
        if(used<2||used>14||saved.at("validation_samples").as_array().size()!=6)throw std::runtime_error("modular completion lacks independent validation samples");
        std::size_t check_index=0;
        for(unsigned pi=used+1;pi<=used+2;++pi)for(std::size_t j=0;j<3;++j) {
          const auto ordinal=100000+100*used+j;
          const auto expected="sample-"+std::to_string(pi)+"-"+std::to_string(ordinal)+".json";
          if(artifacts::detail::string(saved.at("validation_samples").as_array()[check_index++])!=expected)throw std::runtime_error("modular held-out sample provenance mismatch");
          auto payload=detail::read(result.directory/expected);auto point=modular::point(field_.variables().size(),pi,ordinal);
          for(std::size_t i=0;i<point.size();++i)if(std::find(active.begin(),active.end(),i)==active.end())point[i]=3;
          json::array coordinates;for(auto value:point)coordinates.emplace_back(std::to_string(value));
          if(artifacts::detail::string(payload.at("identity"))!=identity||payload.at("point")!=coordinates||payload.at("prime_index").to_number<unsigned>()!=pi||artifacts::detail::string(payload.at("modulus"))!=std::to_string(modular::prime(pi)))throw std::runtime_error("modular held-out sample identity mismatch");
          auto check=load_rows(payload);std::map<ibp::Integral,ibp::Relation> image;
          for(const auto& [a,row]:result.reductions){image.emplace(a,ibp::Relation{});for(const auto& [b,c]:row){auto value=modular::evaluate(c,point,modular::prime(pi));if(value)image.at(a).emplace(b,dimension_.parse(std::to_string(value)));}}
          if(image!=check)throw std::runtime_error("cached modular reconstruction fails independent sample replay");
        }
        verify_native(result.reductions);
        status("reused verified modular reconstruction");return result;
      }
      std::size_t fresh=0,hits=0;
      auto probe=[&](unsigned prime_index,std::size_t ordinal)->detail::Probe {
        remaining();auto point=modular::point(field_.variables().size(),prime_index,ordinal);
        // An unused symbol is absent from the complete native IBP input, not
        // merely absent at one specialization. Keep its coordinate explicit.
        for(std::size_t i=0;i<point.size();++i)if(std::find(active.begin(),active.end(),i)==active.end())point[i]=3;
        fire::ModularPoint input;input.prime_index=prime_index;for(auto value:point)input.values.push_back(value);
        const auto path=result.directory/("sample-"+std::to_string(prime_index)+"-"+std::to_string(ordinal)+".json");
        json::array coordinates;for(auto value:point)coordinates.emplace_back(std::to_string(value));
        json::object payload;
        if(fs::exists(path)) {payload=detail::read(path);++hits;}
        else {
          auto opts=limits;opts.executable=options_.executable;opts.threads=1;opts.simplifier_threads=1;
          opts.max_completed_forward_sectors=opts.max_completed_backward_sectors=0;
          opts.timeout_seconds=std::min(options_.probe_timeout_seconds,remaining());
          auto reduced=fire::reduce(basis_,dimension_,field_,requested,opts,nullptr,nullptr,false,&input);
          if(!reduced.success)throw std::runtime_error("finite-field probe failed: "+reduced.reason);
          fire_batch::validate(reduced,basis_,dimension_,requested);
          auto tables=fire::read_text(reduced.directory/("result_"+input.suffix()+".tables"));
          payload={{"schema","DiffExp3.ModularSample/v1"},{"identity",identity},{"prime_index",prime_index},{"modulus",std::to_string(modular::prime(prime_index))},{"point",coordinates},{"tables",tables},{"directory",reduced.directory.string()}};
          detail::save(path,payload);++fresh;
        }
        if(artifacts::detail::string(payload.at("identity"))!=identity||payload.at("point")!=coordinates||artifacts::detail::string(payload.at("modulus"))!=std::to_string(modular::prime(prime_index))||payload.at("prime_index").to_number<unsigned>()!=prime_index)
          throw std::runtime_error("modular sample identity, point or prime mismatch");
        auto rows=load_rows(payload);fire::Result check;check.success=true;check.reductions=rows;fire_batch::validate(check,basis_,dimension_,requested);
        for(const auto& [a,row]:rows)for(const auto& [b,c]:row)if(!c.is_rational())throw std::runtime_error("nonconstant finite-field output");
        return {point,std::move(rows)};
      };
      std::vector<detail::Probe> probes;std::vector<modular::Image> images;std::vector<detail::Key> keys;
      std::set<ibp::Integral> terminal;
      auto flatten=[&](const auto& all,unsigned pi) {
        std::vector<modular::Sample> out;for(const auto& p:all) {
          if(detail::masters(p)!=terminal)throw std::runtime_error("modular master basis changed at a sample; no reconstruction accepted");
          modular::Sample row;for(auto i:active)row.point.push_back(p.point[i]);
          for(const auto& key:keys)row.coefficients.push_back(detail::coefficient(p,key,modular::prime(pi)));out.push_back(std::move(row));
        }return out;
      };
      std::size_t count=8;
      for(;;) {
        while(probes.size()<count)probes.push_back(probe(1,probes.size()));
        terminal=detail::masters(probes.front());std::set<detail::Key> all_keys;
        for(const auto& p:probes)for(const auto& [a,row]:p.rows)for(const auto& [b,c]:row)all_keys.emplace(a,b);
        keys.assign(all_keys.begin(),all_keys.end());auto samples=flatten(probes,1);images.clear();
        for(std::size_t i=0;i<keys.size();++i){remaining();auto image=modular::discover(samples,i,active.size(),options_.max_degree,modular::prime(1));if(!image)break;images.push_back(std::move(*image));}
        detail::save(result.directory/"progress.json",{{"schema","DiffExp3.ModularProgress/v1"},{"identity",identity},{"stage","degree-discovery"},{"samples_per_prime",probes.size()},{"coefficients",keys.size()},{"discovered",images.size()},{"fresh_samples",fresh},{"reused_samples",hits}});
        status("modular samples "+std::to_string(probes.size())+", reconstructed shapes "+std::to_string(images.size())+"/"+std::to_string(keys.size()));
        if(images.size()==keys.size())break;
        if(count==options_.max_samples_per_prime)throw std::runtime_error("modular sample/degree budget exhausted; progress retained");
        count=std::min(count*2,options_.max_samples_per_prime);
      }
      std::vector<modular::Lift> lifts;for(const auto& image:images)lifts.emplace_back(image,modular::prime(1));
      for(unsigned pi=2;pi<=options_.max_primes;++pi) {
        std::vector<detail::Probe> next;for(std::size_t i=0;i<probes.size();++i)next.push_back(probe(pi,i));
        auto samples=flatten(next,pi);std::vector<modular::Image> next_images;
        for(std::size_t i=0;i<keys.size();++i){remaining();auto image=modular::fit(samples,i,images[i].ansatz,modular::prime(pi),true);if(!image)throw std::runtime_error("modular rank/degree changed on independent reconstruction prime");next_images.push_back(*image);}
        for(std::size_t i=0;i<keys.size();++i)lifts[i].append(next_images[i],modular::prime(pi));
        std::vector<Exact> coefficients;for(const auto& lift:lifts){auto c=lift.reconstruct(dimension_,active);if(!c)break;coefficients.push_back(*c);}
        bool valid=coefficients.size()==keys.size();std::vector<std::string> validation_files;
        if(valid)for(unsigned check_prime=pi+1;check_prime<=pi+2&&valid;++check_prime)for(std::size_t j=0;j<3&&valid;++j) {
          const auto ordinal=100000+100*pi+j;auto check=probe(check_prime,ordinal);
          validation_files.push_back("sample-"+std::to_string(check_prime)+"-"+std::to_string(ordinal)+".json");
          if(detail::masters(check)!=terminal){valid=false;break;}
          std::map<ibp::Integral,ibp::Relation> candidate;for(const auto& [a,row]:probes.front().rows)candidate.emplace(a,ibp::Relation{});
          try{for(std::size_t k=0;k<keys.size();++k){const auto actual=modular::evaluate(coefficients[k],check.point,modular::prime(check_prime));
              if(actual!=detail::coefficient(check,keys[k],modular::prime(check_prime))){valid=false;break;}
              if(actual)candidate.at(keys[k].first).emplace(keys[k].second,dimension_.parse(std::to_string(actual)));}
            if(valid&&candidate!=check.rows)valid=false;
          }catch(const std::domain_error&){valid=false;}
        }
        json::array raw_lifts;for(const auto& lift:lifts){json::array residues;for(const auto& r:lift.residues){char* s=fmpz_get_str(nullptr,10,r.value);residues.emplace_back(s);flint_free(s);}raw_lifts.push_back(residues);}
        detail::save(result.directory/"progress.json",{{"schema","DiffExp3.ModularProgress/v1"},{"identity",identity},{"stage","CRT-and-validation"},{"reconstruction_primes",pi},{"samples_per_prime",probes.size()},{"residues",raw_lifts},{"fresh_samples",fresh},{"reused_samples",hits},{"candidate_validated",valid}});
        status("modular reconstruction primes "+std::to_string(pi)+", independent validation "+(valid?"passed":"requires more primes"));
        if(!valid)continue;
        for(const auto& [a,row]:probes.front().rows)result.reductions.emplace(a,ibp::Relation{});
        for(std::size_t i=0;i<keys.size();++i)ibp::add(result.reductions.at(keys[i].first),keys[i].second,coefficients[i]);
        result.success=true;fire_batch::validate(result,basis_,dimension_,requested);auto verified=verify_native(result.reductions);
        json::array checks;for(const auto& file:validation_files)checks.emplace_back(file);
        auto tables=detail::table(result.reductions,symbols);
        detail::save(result.directory/"completed.json",{{"schema","DiffExp3.ModularCompleted/v1"},{"identity",identity},{"tables",tables},{"verification","fresh-two-primes-three-points-and-applicable-native-IBPs"},{"validation_samples",checks},{"reconstruction_primes",pi},{"samples_per_prime",probes.size()},{"native_ibp_identities_checked",verified},{"assumptions","FIRE finite-field IBP reductions and declared zero sectors; rational reconstruction verified probabilistically, with exact checks of applicable native identities. Level closure is checked separately."}});
        status("accepted modular reduction; exact native IBP checks "+std::to_string(verified));return result;
      }
      throw std::runtime_error("modular rational coefficient prime budget exhausted; all samples retained");
    }catch(const std::exception& e){result.success=false;result.reason=e.what();return result;}
  }
 private:
  ibp::PropagatorBasis basis_;Exact dimension_;ExactField field_;Options options_;
  std::unique_ptr<fire_batch::Cache> exact_batches_;
};
} // namespace diffexp::fire_modular
