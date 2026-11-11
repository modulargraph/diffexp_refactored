#pragma once
#include "diffexp/fire_batch_cache.hpp"
#include "diffexp/modular_reconstruction.hpp"

namespace diffexp::fire_modular {
namespace json=boost::json;
namespace fs=std::filesystem;
struct Options {
  fs::path executable,cache_directory;
  using SampleProvider=std::function<fire::Result(const std::vector<ibp::Integral>&,const std::vector<modular::Word>&,modular::Word,const fire::Options&)>;
  // validation_provider serves held-out points. anchor_provider, when set, serves
  // the first held-out point of every validation attempt: a provider whose
  // ordinary held-out reductions are cheaper than its strongest one keeps at
  // least one strongest reduction per accepted candidate.
  SampleProvider sample_provider,validation_provider,anchor_provider;
  // Optional hints for a provider that works ahead on other threads; they must
  // not change what it returns. upcoming_samples lists the points of the next
  // sample_provider calls at one prime, in order. upcoming_validation lists the
  // held-out points of the next validation attempt, which may never be made,
  // and the seconds the reconstruction has left to make it.
  // Each call replaces the previous one of its kind; an empty list ends it.
  struct HeldOut {unsigned prime_index=0;std::size_t ordinal=0;std::vector<modular::Word> point;bool anchor=false;};
  std::function<void(const std::vector<std::vector<modular::Word>>&,modular::Word)> upcoming_samples;
  std::function<void(const std::vector<HeldOut>&,double)> upcoming_validation;
  bool sparse_lifting=false,degree_slices=false;
  std::function<modular::Word(unsigned)> modulus=modular::prime;
  std::string provider_identity="FIRE7p",prime_table="FIRE7-primes-1-through-16";
  // max_degree caps the dense total-degree scan. Coordinate slices and their
  // per-variable boxes have a separate cap: a coefficient of low degree in
  // epsilon may need a high degree in the merge parameter alone.
  unsigned max_degree=16,max_variable_degree=64,max_primes=8,probe_timeout_seconds=20;
  // Primes modulus() can supply; two beyond the reconstruction primes validate.
  unsigned prime_count=16;
  std::size_t max_samples_per_prime=512;
  std::size_t sample_checkpoint_interval=1;
  // A sample point may be exceptional for reasons unrelated to the reduction: a
  // pole of an input coefficient, or a coefficient vanishing there so that the
  // terminal integrals change. Such points are replaced by the next ordinal, at
  // most this many times per reconstruction.
  std::size_t max_skipped_samples=4;
  std::function<void(const std::string&)> progress;
};
// A point is identified by (prime index, ordinal) alone. The ordinal ranges are
// disjoint, so every purpose owns its points and a resumed run revisits them.
namespace ordinals {
inline constexpr std::size_t slice_base=50000,slice_stride=1000,slice_axes=4; // one-coordinate degree slices, per active axis
inline constexpr std::size_t planning_base=60000;                             // closure-planning probes of native providers
inline constexpr std::size_t validation_base=100000,validation_stride=100;    // held-out points, indexed by reconstruction primes used
inline constexpr std::size_t validation_points=3;
inline constexpr unsigned validation_primes=2;
inline bool slice(std::size_t ordinal){return ordinal>=slice_base&&ordinal<slice_base+slice_stride*slice_axes;}
inline bool validation(std::size_t ordinal){return ordinal>=validation_base;}
inline std::size_t validation_point(unsigned reconstruction_primes,std::size_t j){return validation_base+validation_stride*reconstruction_primes+j;}
}
namespace detail {
struct SliceFailure:std::runtime_error {using std::runtime_error::runtime_error;};
struct ExceptionalPoint:std::runtime_error {using std::runtime_error::runtime_error;};
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
// A finite-field table holds residues. Each is read once into a machine word;
// interpolation never returns to the exact-field representation.
using Rows=std::map<ibp::Integral,std::map<ibp::Integral,modular::Word>>;
struct Probe {std::vector<modular::Word> point;Rows rows;std::size_t ordinal=0;};
inline Rows residues(const std::map<ibp::Integral,ibp::Relation>& rows,const std::vector<modular::Word>& point,modular::Word prime) {
  Rows out;for(const auto& [a,row]:rows){auto& target=out[a];for(const auto& [b,c]:row){
    if(!c.is_rational())throw std::runtime_error("nonconstant finite-field output");target.emplace(b,modular::evaluate(c,point,prime));}}
  return out;
}
inline std::set<ibp::Integral> masters(const Probe& probe) {
  std::set<ibp::Integral> out;for(const auto& [a,row]:probe.rows)for(const auto& [b,c]:row)out.insert(b);return out;
}
// The point of a sample: its ordinal's, on the line of its axis if it belongs
// to a degree slice. Every cached sample records and verifies this complete
// actual point, including the specialised coordinates of a slice.
inline std::vector<modular::Word> sample_point(std::size_t variables,const std::vector<std::size_t>& active,bool slices,unsigned prime_index,std::size_t ordinal) {
  auto point=modular::point(variables,prime_index,ordinal);
  if(slices&&ordinals::slice(ordinal)){
    const auto axis=(ordinal-ordinals::slice_base)/ordinals::slice_stride;
    if(axis>=active.size())throw std::invalid_argument("degree slice axis");
    auto line=modular::point(variables,prime_index,ordinals::slice_base+ordinals::slice_stride*axis);
    line[active[axis]]=point[active[axis]];point=std::move(line);
  }
  // An unused symbol is absent from the complete native IBP input, not
  // merely absent at one specialization. Keep its coordinate explicit.
  for(std::size_t i=0;i<point.size();++i)if(std::find(active.begin(),active.end(),i)==active.end())point[i]=3;
  return point;
}
inline modular::Word coefficient(const Probe& sample,const Key& key) {
  const auto& row=sample.rows.at(key.first);auto it=row.find(key.second);return it==row.end()?0:it->second;
}
} // namespace detail

// Completed symbolic batches remain reusable. An uncached demand uses only
// the selected finite-field provider, native FLINT interpolation and CRT; there is no symbolic fallback.
// The generic level preparer subsequently proves its derivative/target closure
// conditional on these imported identities, exactly as for symbolic FIRE.
class Session {
 public:
  Session(ibp::PropagatorBasis basis,Exact dimension,ExactField field,Options options,fs::path exact_batches={})
    :basis_(std::move(basis)),dimension_(std::move(dimension)),field_(std::move(field)),options_(std::move(options)) {
    if(options_.cache_directory.empty()||(options_.executable.empty()&&!options_.sample_provider)||!options_.modulus||options_.provider_identity.empty()||options_.prime_table.empty()||!options_.max_degree||options_.max_degree>256||!options_.max_variable_degree||options_.max_variable_degree>4096||options_.max_primes<2||options_.prime_count<4||options_.max_primes>options_.prime_count-ordinals::validation_primes||options_.max_samples_per_prime<8||options_.max_samples_per_prime>(std::size_t(1)<<20)||!options_.probe_timeout_seconds)
      throw std::invalid_argument("invalid finite modular reconstruction budgets or paths");
    if(!options_.sample_checkpoint_interval||options_.sample_checkpoint_interval>4096)throw std::invalid_argument("modular checkpoint interval must be 1..4096");
    if(options_.max_skipped_samples>64)throw std::invalid_argument("modular exceptional-sample budget must be 0..64");
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
      json::object identity_input{{"schema","DiffExp3.NativeModularProvider/v1"},{"science",semantic.json_value()},{"active",active_json},{"prime_table",options_.prime_table}};
      // Preserve existing FIRE sample/checkpoint identities. Native providers
      // additionally bind their implementation and seed policy.
      if(options_.sample_provider||options_.provider_identity!="FIRE7p")identity_input["provider"]=options_.provider_identity;
      if(options_.degree_slices)identity_input["degree_discovery"]="coordinate-slices-v2";
      const auto identity=artifacts::detail::sha256(artifacts::detail::canonical(identity_input));
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
        if(used<2||used>options_.prime_count-ordinals::validation_primes||saved.at("validation_samples").as_array().size()!=ordinals::validation_primes*ordinals::validation_points)throw std::runtime_error("modular completion lacks independent validation samples");
        std::size_t check_index=0;
        for(unsigned pi=used+1;pi<=used+ordinals::validation_primes;++pi)for(std::size_t j=0;j<ordinals::validation_points;++j) {
          const auto ordinal=ordinals::validation_point(used,j);
          const auto expected="sample-"+std::to_string(pi)+"-"+std::to_string(ordinal)+".json";
          if(artifacts::detail::string(saved.at("validation_samples").as_array()[check_index++])!=expected)throw std::runtime_error("modular held-out sample provenance mismatch");
          auto payload=detail::read(result.directory/expected);auto point=modular::point(field_.variables().size(),pi,ordinal);
          for(std::size_t i=0;i<point.size();++i)if(std::find(active.begin(),active.end(),i)==active.end())point[i]=3;
          json::array coordinates;for(auto value:point)coordinates.emplace_back(std::to_string(value));
          if(artifacts::detail::string(payload.at("identity"))!=identity||payload.at("point")!=coordinates||payload.at("prime_index").to_number<unsigned>()!=pi||artifacts::detail::string(payload.at("modulus"))!=std::to_string(options_.modulus(pi)))throw std::runtime_error("modular held-out sample identity mismatch");
          auto check=load_rows(payload);std::map<ibp::Integral,ibp::Relation> image;
          for(const auto& [a,row]:result.reductions){image.emplace(a,ibp::Relation{});for(const auto& [b,c]:row){auto value=modular::evaluate(c,point,options_.modulus(pi));if(value)image.at(a).emplace(b,dimension_.parse(std::to_string(value)));}}
          if(image!=check)throw std::runtime_error("cached modular reconstruction fails independent sample replay");
        }
        verify_native(result.reductions);
        status("reused verified modular reconstruction");return result;
      }
      std::size_t fresh=0,hits=0;
      auto locate=[&](unsigned prime_index,std::size_t ordinal){return detail::sample_point(field_.variables().size(),active,options_.degree_slices,prime_index,ordinal);};
      auto record=[&](unsigned prime_index,std::size_t ordinal){return result.directory/("sample-"+std::to_string(prime_index)+"-"+std::to_string(ordinal)+".json");};
      // Points a provider may work on ahead: those of the given ordinals that
      // have no record yet and will therefore be requested from it.
      auto announce=[&](unsigned prime_index,std::size_t first,std::size_t count) {
        if(!options_.upcoming_samples)return;std::vector<std::vector<modular::Word>> points;
        for(auto ordinal=first;ordinal<first+count;++ordinal)if(!fs::exists(record(prime_index,ordinal)))points.push_back(locate(prime_index,ordinal));
        options_.upcoming_samples(points,options_.modulus(prime_index));
      };
      // The held-out points of a candidate lifted from `primes` primes that
      // have no record yet. They do not depend on the candidate.
      auto foresee=[&](unsigned primes) {
        if(!options_.upcoming_validation||primes>options_.max_primes)return;auto points=held_out(active,primes);
        std::erase_if(points,[&](const auto& held){return fs::exists(record(held.prime_index,held.ordinal));});
        options_.upcoming_validation(points,std::max(0.0,limits.timeout_seconds-std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()));
      };
      // No hint outlives the reconstruction, however it ends.
      struct Silence{const Options& options;~Silence(){try{if(options.upcoming_samples)options.upcoming_samples({},0);if(options.upcoming_validation)options.upcoming_validation({},0);}catch(...){}}}silence{options_};
      auto probe=[&](unsigned prime_index,std::size_t ordinal,bool anchor=false)->detail::Probe {
        remaining();auto point=locate(prime_index,ordinal);
        const bool slice=options_.degree_slices&&ordinals::slice(ordinal);
        fire::ModularPoint input;input.prime_index=prime_index;for(auto value:point)input.values.push_back(value);
        const auto path=record(prime_index,ordinal);
        const bool retain=ordinals::validation(ordinal)||ordinal%options_.sample_checkpoint_interval==0;
        json::array coordinates;for(auto value:point)coordinates.emplace_back(std::to_string(value));
        json::object payload;std::optional<detail::Rows> rows;
        if(fs::exists(path)) {payload=detail::read(path);++hits;}
        else {
          auto opts=limits;opts.executable=options_.executable;opts.threads=1;opts.simplifier_threads=1;
          opts.max_completed_forward_sectors=opts.max_completed_backward_sectors=0;
          opts.timeout_seconds=std::min(options_.probe_timeout_seconds,remaining());
          const auto& held_out=anchor&&options_.anchor_provider?options_.anchor_provider:options_.validation_provider;
          const auto& provider=ordinals::validation(ordinal)&&held_out?held_out:options_.sample_provider;
          auto reduced=provider?provider(requested,point,options_.modulus(prime_index),opts):fire::reduce(basis_,dimension_,field_,requested,opts,nullptr,nullptr,false,&input);
          if(!reduced.success){
            if(slice)throw detail::SliceFailure("coordinate probe failed: "+reduced.reason);
            if(reduced.exceptional_point)throw detail::ExceptionalPoint(reduced.reason);
            throw std::runtime_error("finite-field probe failed: "+reduced.reason);
          }
          fire_batch::validate(reduced,basis_,dimension_,requested);
          auto tables=provider?(retain?detail::table(reduced.reductions,symbols):std::string{}):fire::read_text(reduced.directory/("result_"+input.suffix()+".tables"));
          payload={{"schema","DiffExp3.ModularSample/v1"},{"identity",identity},{"prime_index",prime_index},{"modulus",std::to_string(options_.modulus(prime_index))},{"point",coordinates},{"tables",tables},{"directory",reduced.directory.string()}};
          if(retain)detail::save(path,payload);if(provider)rows=detail::residues(reduced.reductions,point,options_.modulus(prime_index));++fresh;
        }
        if(artifacts::detail::string(payload.at("identity"))!=identity||payload.at("point")!=coordinates||artifacts::detail::string(payload.at("modulus"))!=std::to_string(options_.modulus(prime_index))||payload.at("prime_index").to_number<unsigned>()!=prime_index)
          throw std::runtime_error("modular sample identity, point or prime mismatch");
        if(!rows){fire::Result check;check.success=true;check.reductions=load_rows(payload);fire_batch::validate(check,basis_,dimension_,requested);rows=detail::residues(check.reductions,point,options_.modulus(prime_index));}
        return {point,std::move(*rows),ordinal};
      };
      std::vector<detail::Probe> probes;std::vector<modular::Image> images;std::vector<detail::Key> keys;
      std::set<ibp::Integral> terminal;bool voted=false;
      auto flatten=[&](const auto& all) {
        std::vector<modular::Sample> out;for(const auto& p:all) {
          if(detail::masters(p)!=terminal)throw std::runtime_error("modular master basis changed at a sample; no reconstruction accepted");
          modular::Sample row;for(auto i:active)row.point.push_back(p.point[i]);
          for(const auto& key:keys)row.coefficients.push_back(detail::coefficient(p,key));out.push_back(std::move(row));
        }return out;
      };
      // Ordinary samples of a prime use consecutive ordinals; an exceptional
      // point is skipped. Points and the two skip criteria are deterministic, so
      // a resumed run repeats the same decisions; the record is provenance.
      std::map<unsigned,std::size_t> cursor;json::array skipped;
      auto skip=[&](unsigned prime_index,std::size_t ordinal,const std::string& why) {
        if(skipped.size()>=options_.max_skipped_samples)throw std::runtime_error("modular exceptional-sample budget ("+std::to_string(options_.max_skipped_samples)+") exhausted: "+why);
        skipped.push_back(json::object{{"prime_index",prime_index},{"ordinal",ordinal},{"reason",why}});
        status("skipped exceptional sample "+std::to_string(ordinal)+" at prime "+std::to_string(prime_index)+": "+why);
      };
      const std::string moved="terminal integrals differ from the reference samples";
      auto acquire=[&](unsigned prime_index,std::size_t wanted,std::vector<detail::Probe>& into) {
        std::size_t announced=0;
        while(into.size()<wanted){const auto ordinal=cursor[prime_index]++;
          // Every missing sample takes at least one further ordinal.
          if(ordinal>=announced){announce(prime_index,ordinal,wanted-into.size());announced=ordinal+wanted-into.size();}
          try{auto sample=probe(prime_index,ordinal);
            if(voted&&detail::masters(sample)!=terminal)skip(prime_index,ordinal,moved);else into.push_back(std::move(sample));
          }catch(const detail::ExceptionalPoint& e){skip(prime_index,ordinal,e.what());}
        }
      };
      // One exceptional point must not define the terminal integrals: the
      // reference is the set shared by at least two of the first three samples.
      foresee(2);acquire(1,3,probes);
      {std::size_t reference=probes.size();
        for(std::size_t i=0;i<probes.size()&&reference==probes.size();++i)for(std::size_t j=i+1;j<probes.size();++j)if(detail::masters(probes[i])==detail::masters(probes[j])){reference=i;break;}
        if(reference==probes.size())throw std::runtime_error("modular master basis differs at each of the first samples; no reconstruction accepted");
        terminal=detail::masters(probes[reference]);voted=true;
        for(std::size_t i=probes.size();i-->0;)if(detail::masters(probes[i])!=terminal){skip(1,probes[i].ordinal,moved);probes.erase(probes.begin()+i);}
      }
      std::size_t count=8;
      bool sliced=false;
      std::map<detail::Key,modular::Ansatz> boxes;
      std::map<detail::Key,modular::Image> discovered;
      const auto describe=[](const detail::Key& key){return "coefficient of "+fire::indices(key.second)+" in "+fire::indices(key.first);};
      const std::function<void()> deadline=[&]{remaining();};
      for(;;) {
        acquire(1,count,probes);std::set<detail::Key> all_keys;
        for(const auto& p:probes)for(const auto& [a,row]:p.rows)for(const auto& [b,c]:row)all_keys.emplace(a,b);
        keys.assign(all_keys.begin(),all_keys.end());auto samples=flatten(probes);images.clear();
        bool complete=true,dense_pending=false;std::size_t box_samples=0;std::optional<detail::Key> saturated;
        for(std::size_t i=0;i<keys.size();++i){remaining();std::optional<modular::Image> image;
          if(auto old=discovered.find(keys[i]);old!=discovered.end()){
            bool valid=true;for(const auto& sample:samples)if(modular::evaluate(old->second,sample.point,options_.modulus(1))!=std::optional<modular::Word>(sample.coefficients[i])){valid=false;break;}
            if(valid)image=old->second;else discovered.erase(old);
          }
          if(!image)if(auto box=boxes.find(keys[i]);box!=boxes.end()) {
            // A box states exactly how many generic points it needs. Wait for
            // them instead of repeating the dense scan at every sample count;
            // a box refuted with enough points returns to the dense scan.
            const auto required=box->second.numerator.size()+box->second.denominator.size()+3;
            if(samples.size()<required)box_samples=std::max(box_samples,required);
            else if(!(image=modular::fit_with_holdouts(samples,i,box->second,options_.modulus(1))))boxes.erase(box);
          }
          if(!image&&!boxes.contains(keys[i])) {
            image=modular::discover(samples,i,active.size(),options_.max_degree,options_.modulus(1),deadline);
            if(!image){if(modular::dense_scan_saturated(samples.size(),active.size(),options_.max_degree)){if(!saturated)saturated=keys[i];}else dense_pending=true;}
          }
          if(image){discovered.insert_or_assign(keys[i],modular::compact_support(*image));images.push_back(std::move(*image));}else complete=false;
        }
        detail::save(result.directory/"progress.json",{{"schema","DiffExp3.ModularProgress/v1"},{"identity",identity},{"stage","degree-discovery"},{"samples_per_prime",probes.size()},{"coefficients",keys.size()},{"discovered",images.size()},{"fresh_samples",fresh},{"reused_samples",hits},{"skipped_samples",skipped}});
        status("modular samples "+std::to_string(probes.size())+", reconstructed shapes "+std::to_string(images.size())+"/"+std::to_string(keys.size()));
        if(complete)break;
        if(options_.degree_slices&&!sliced&&active.size()>1&&active.size()<=ordinals::slice_axes){
          sliced=true;std::vector<std::vector<std::optional<modular::Image>>> axes;bool valid=true;
          for(std::size_t axis=0;axis<active.size()&&valid;++axis){
            std::vector<detail::Probe> line;std::size_t wanted=8;std::vector<std::optional<modular::Image>> shapes(keys.size());
            const auto maximum=std::min<std::size_t>(2*options_.max_variable_degree+4,options_.max_samples_per_prime);
            for(;;){
              announce(1,ordinals::slice_base+ordinals::slice_stride*axis+line.size(),wanted-line.size());
              try{while(line.size()<wanted){auto sample=probe(1,ordinals::slice_base+ordinals::slice_stride*axis+line.size());if(detail::masters(sample)!=terminal){valid=false;break;}line.push_back(std::move(sample));}}
              catch(const detail::SliceFailure& e){status(std::string(e.what())+"; using multivariate discovery");valid=false;}
              if(!valid)break;auto slice=flatten(line);for(auto& s:slice)s.point={s.point[axis]};std::size_t found=0;
              for(std::size_t i=0;i<keys.size();++i){remaining();
                if(shapes[i])for(const auto& sample:slice)if(modular::evaluate(*shapes[i],sample.point,options_.modulus(1))!=std::optional<modular::Word>(sample.coefficients[i])){shapes[i].reset();break;}
                if(!shapes[i])shapes[i]=modular::discover(slice,i,1,options_.max_variable_degree,options_.modulus(1),deadline);
                found+=shapes[i].has_value();
              }
              if(found==keys.size()||wanted==maximum){status("coordinate "+std::to_string(axis)+" degrees "+std::to_string(found)+"/"+std::to_string(keys.size()));break;}wanted=std::min(2*wanted,maximum);
            }
            if(valid)axes.push_back(std::move(shapes));
          }
          if(valid)for(std::size_t i=0;i<keys.size();++i){if(discovered.contains(keys[i]))continue;std::vector<modular::Image> slice;for(const auto& axis:axes){if(!axis[i])break;slice.push_back(*axis[i]);}
            if(slice.size()==active.size())if(auto box=modular::degree_box(slice,options_.max_variable_degree)){
              box_samples=std::max(box_samples,box->numerator.size()+box->denominator.size()+3);boxes.emplace(keys[i],std::move(*box));saturated.reset();}}
          status("coordinate degree hypotheses "+std::to_string(boxes.size())+"/"+std::to_string(keys.size()-discovered.size()));
          // Keys without a box still depend on the dense scan below.
          for(const auto& key:keys)if(!discovered.contains(key)&&!boxes.contains(key)){
            if(modular::dense_scan_saturated(samples.size(),active.size(),options_.max_degree)){if(!saturated)saturated=key;}else dense_pending=true;}
        }
        // Name the exhausted budget. A saturated dense scan cannot be helped
        // by more points; a box may simply need more than the sample cap.
        if(saturated&&!boxes.contains(*saturated))
          throw std::runtime_error("modular total-degree budget ("+std::to_string(options_.max_degree)+") exhausted for "+describe(*saturated)+
            "; enable coordinate degree slices or raise max_degree; progress retained");
        if(box_samples>options_.max_samples_per_prime)
          throw std::runtime_error("modular sample budget exhausted: a coefficient needs "+std::to_string(box_samples)+" points, the limit is "+
            std::to_string(options_.max_samples_per_prime)+"; progress retained");
        auto next=std::max(count,box_samples);if(dense_pending)next=std::max(next,count*2);
        next=std::min(next,options_.max_samples_per_prime);
        if(next==count)throw std::runtime_error("modular sample/degree budget exhausted; progress retained");
        count=next;
      }
      if(options_.sparse_lifting)for(auto& image:images)image=modular::compact_support(image);
      std::size_t lifting_samples=probes.size();
      if(options_.sparse_lifting){lifting_samples=8;for(const auto& image:images)lifting_samples=std::max(lifting_samples,image.coefficients.size()+2);lifting_samples=std::min(lifting_samples,probes.size());}
      std::vector<modular::Lift> lifts;for(const auto& image:images)lifts.emplace_back(image,options_.modulus(1));
      for(unsigned pi=2;pi<=options_.max_primes;++pi) {
        std::vector<detail::Probe> next;acquire(pi,lifting_samples,next);
        auto samples=flatten(next);std::vector<modular::Image> next_images;
        for(std::size_t i=0;i<keys.size();++i){remaining();auto image=modular::fit(samples,i,images[i].ansatz,options_.modulus(pi),true);if(!image)throw std::runtime_error("modular rank/degree changed on independent reconstruction prime");next_images.push_back(*image);}
        for(std::size_t i=0;i<keys.size();++i)lifts[i].append(next_images[i],options_.modulus(pi));
        std::vector<Exact> coefficients;for(const auto& lift:lifts){auto c=lift.reconstruct(dimension_,active);if(!c)break;coefficients.push_back(*c);}
        bool valid=coefficients.size()==keys.size();std::vector<std::string> validation_files;
        if(valid)for(unsigned check_prime=pi+1;check_prime<=pi+ordinals::validation_primes&&valid;++check_prime)for(std::size_t j=0;j<ordinals::validation_points&&valid;++j) {
          const auto ordinal=ordinals::validation_point(pi,j);std::optional<detail::Probe> held;
          // Held-out ordinals are fixed by the completion record. An exceptional
          // one rejects this attempt; the next prime brings fresh points.
          try{held=probe(check_prime,ordinal,anchor(pi,check_prime,j));}
          catch(const detail::ExceptionalPoint& e){status(std::string("exceptional held-out point: ")+e.what());valid=false;break;}
          const auto& check=*held;
          validation_files.push_back("sample-"+std::to_string(check_prime)+"-"+std::to_string(ordinal)+".json");
          if(detail::masters(check)!=terminal){valid=false;break;}
          detail::Rows candidate;for(const auto& [a,row]:probes.front().rows)candidate[a];
          try{for(std::size_t k=0;k<keys.size();++k){const auto actual=modular::evaluate(coefficients[k],check.point,options_.modulus(check_prime));
              if(actual!=detail::coefficient(check,keys[k])){valid=false;break;}
              if(actual)candidate.at(keys[k].first).emplace(keys[k].second,actual);}
            if(valid&&candidate!=check.rows)valid=false;
          }catch(const std::domain_error&){valid=false;}
        }
        json::array raw_lifts;for(const auto& lift:lifts){json::array residues;for(const auto& r:lift.residues){char* s=fmpz_get_str(nullptr,10,r.value);residues.emplace_back(s);flint_free(s);}raw_lifts.push_back(residues);}
        detail::save(result.directory/"progress.json",{{"schema","DiffExp3.ModularProgress/v1"},{"identity",identity},{"stage","CRT-and-validation"},{"reconstruction_primes",pi},{"samples_per_prime",probes.size()},{"residues",raw_lifts},{"fresh_samples",fresh},{"reused_samples",hits},{"skipped_samples",skipped},{"candidate_validated",valid}});
        status("modular reconstruction primes "+std::to_string(pi)+", independent validation "+(valid?"passed":"requires more primes"));
        if(!valid){foresee(pi+1);continue;}
        for(const auto& [a,row]:probes.front().rows)result.reductions.emplace(a,ibp::Relation{});
        for(std::size_t i=0;i<keys.size();++i)ibp::add(result.reductions.at(keys[i].first),keys[i].second,coefficients[i]);
        result.success=true;fire_batch::validate(result,basis_,dimension_,requested);auto verified=verify_native(result.reductions);
        json::array checks;for(const auto& file:validation_files)checks.emplace_back(file);
        auto tables=detail::table(result.reductions,symbols);
        detail::save(result.directory/"completed.json",{{"schema","DiffExp3.ModularCompleted/v1"},{"identity",identity},{"tables",tables},{"verification","fresh-two-primes-three-points-and-applicable-native-IBPs"},{"validation_samples",checks},{"reconstruction_primes",pi},{"samples_per_prime",probes.size()},{"skipped_samples",skipped},{"native_ibp_identities_checked",verified},{"assumptions",options_.provider_identity+" finite-field IBP reductions and declared zero sectors; rational reconstruction verified probabilistically, with exact checks of applicable native identities. Level closure is checked separately."}});
        status("accepted modular reduction; exact native IBP checks "+std::to_string(verified));return result;
      }
      throw std::runtime_error("modular rational coefficient prime budget exhausted; all samples retained");
    }catch(const std::exception& e){result.success=false;result.reason=e.what();return result;}
  }
  // The held-out points that validate a candidate lifted from `primes` primes,
  // in the order they are probed. They depend on the family and that number
  // alone, so a provider can be told of them before the batch is known.
  std::vector<Options::HeldOut> held_out(const fire::Options& limits,unsigned primes)const {
    fire::SymbolMap symbols(field_.variables());return held_out(detail::active_variables(fire::start(basis_,dimension_,limits.zero_sectors,&symbols),symbols),primes);
  }
 private:
  static bool anchor(unsigned primes,unsigned check_prime,std::size_t j){return check_prime==primes+1&&j==0;}
  std::vector<Options::HeldOut> held_out(const std::vector<std::size_t>& active,unsigned primes)const {
    std::vector<Options::HeldOut> out;if(primes+ordinals::validation_primes>options_.prime_count)return out;
    for(unsigned check_prime=primes+1;check_prime<=primes+ordinals::validation_primes;++check_prime)for(std::size_t j=0;j<ordinals::validation_points;++j)
      out.push_back({check_prime,ordinals::validation_point(primes,j),detail::sample_point(field_.variables().size(),active,options_.degree_slices,check_prime,ordinals::validation_point(primes,j)),anchor(primes,check_prime,j)});
    return out;
  }
  ibp::PropagatorBasis basis_;Exact dimension_;ExactField field_;Options options_;
  std::unique_ptr<fire_batch::Cache> exact_batches_;
};
} // namespace diffexp::fire_modular
