#pragma once
#include "diffexp/fire_batch_cache.hpp"
#include <regex>
#include <iomanip>

namespace diffexp::fire_checkpoint {
namespace json=boost::json;
inline bool sector_file(const std::string& name) {
  return name.size()==9&&name.substr(5)==".tmp"&&std::all_of(name.begin(),name.begin()+5,[](unsigned char c){return std::isdigit(c);});
}
inline json::object identity(const ibp::PropagatorBasis& b,const Exact& d,const ExactField& field,
    const std::vector<ibp::Integral>& requests,const fire::Options& options,
    const std::map<ibp::Integral,ibp::Relation>* prior) {
  fire::SymbolMap symbols(field.variables());auto scientific=fire_batch::identity(b,d,field,requests,options).json_value();
  return {{"schema","DiffExp3.FixedFIRECheckpoint/v1"},{"scientific",scientific},
    {"binary_sha256",artifacts::detail::sha256(fire::read_text(options.executable))},
    {"worker_binary_sha256",artifacts::detail::sha256(fire::read_text(options.executable.parent_path()/"FLAME7"))},
    {"threads",options.threads},{"simplifier_threads",options.simplifier_threads},
    {"prior_rules",prior?fire::rules_text(*prior,d,symbols):""},
    {"storage_policy","atomic forward snapshots only; discard all level flags; rebuild request propagation"}};
}
inline std::string validated_snapshot(const std::filesystem::path& path) {
  auto bytes=fire::read_text(path);if(bytes.size()<6||bytes.substr(0,6)!=std::string("KCSS\n\0",6))throw std::runtime_error("invalid FIRE snapshot magic");
  std::size_t pos=6,records=0;
  auto number=[&](){std::size_t value=0;unsigned parts=0;while(true){if(pos==bytes.size()||++parts>9)throw std::runtime_error("truncated snapshot integer");auto c=static_cast<unsigned char>(bytes[pos++]);if(value>bytes.size()/128)throw std::runtime_error("snapshot record exceeds byte budget");value=value*128+(c&127);if(c<128)return value;}};
  while(pos<bytes.size()){if(bytes[pos++]!=0||++records>20000000)throw std::runtime_error("invalid FIRE snapshot record");auto k=number(),v=number();if(k>bytes.size()-pos||v>bytes.size()-pos-k)throw std::runtime_error("truncated FIRE snapshot record");pos+=k+v;}
  return bytes;
}
inline json::array files(const std::filesystem::path& storage) {
  json::array out;std::vector<std::filesystem::path> paths;
  for(const auto& entry:std::filesystem::directory_iterator(storage))if(sector_file(entry.path().filename().string())) {
    if(entry.is_symlink()||!entry.is_regular_file())throw std::runtime_error("checkpoint has invalid sector file");paths.push_back(entry.path());
  }
  if(paths.size()>4096)throw std::runtime_error("checkpoint sector count budget");std::sort(paths.begin(),paths.end());std::uintmax_t total=0;
  for(const auto& path:paths){total+=std::filesystem::file_size(path);if(total>1024ull*1024*1024)throw std::runtime_error("checkpoint snapshot byte budget");out.push_back(json::object{{"name",path.filename().string()},{"sha256",artifacts::detail::sha256(validated_snapshot(path))}});}
  return out;
}
inline std::set<unsigned> completed_backward(const std::string& log) {
  std::set<unsigned> sectors;std::regex pattern(R"(FLAME time \(-([0-9]+)\):)");
  for(auto i=std::sregex_iterator(log.begin(),log.end(),pattern);i!=std::sregex_iterator();++i){auto n=std::stoul((*i)[1]);if(n<=1||n>99999||!sectors.insert(n).second)throw std::runtime_error("ambiguous backward checkpoint finalization evidence");}
  return sectors;
}
inline std::string sector_name(unsigned n){std::ostringstream out;out<<std::setw(5)<<std::setfill('0')<<n<<".tmp";return out.str();}
inline void copy_snapshot_list(const std::filesystem::path& source,const std::filesystem::path& destination,
    const json::array& list,bool replace=false) {
  std::filesystem::create_directories(destination);
  for(const auto& row:list){auto name=artifacts::detail::string(row.as_object().at("name"));if(!sector_file(name))throw std::runtime_error("invalid checkpoint snapshot name");std::filesystem::copy_file(source/name,destination/name,replace?std::filesystem::copy_options::overwrite_existing:std::filesystem::copy_options::none);}
}
// The identity's storage_policy describes the baseline forward-storage protocol.
// This explicit, separately checked extension adds finalized backward rows without
// changing scientific inputs, demands, or the old forward-only checkpoint keys.
inline void attach_backward(const std::filesystem::path& directory,json::object& manifest,const std::filesystem::path& completed_run={}) {
  const auto& source=completed_run.empty()?directory:completed_run;
  if(!std::filesystem::exists(source/"run.log"))return;
  auto log=fire::read_text(source/"run.log");auto closed=completed_backward(log);if(closed.empty())return;
  json::array retained;auto destination=directory/"backward";std::filesystem::create_directory(destination);
  for(auto n:closed){auto name=sector_name(n);auto original=source/"database"/name;if(!std::filesystem::exists(original))continue;
    auto bytes=validated_snapshot(original);std::ofstream out(destination/name,std::ios::binary);out<<bytes;if(!out)throw std::runtime_error("cannot persist finalized backward snapshot");retained.push_back(n);}
  if(retained.empty())return;
  {std::ofstream out(directory/"checkpoint-run.log");out<<log;if(!out)throw std::runtime_error("cannot persist backward finalization provenance");}
  manifest["backward_overlay"]=json::object{{"files",files(destination)},{"closed_sectors",retained},{"log_sha256",artifacts::detail::sha256(log)}};
}
inline void validate_manifest(const std::filesystem::path& directory,const json::object& manifest) {
  if(manifest.contains("backward_overlay"))artifacts::detail::keys(manifest,{"identity","files","state","backward_overlay"});
  else artifacts::detail::keys(manifest,{"identity","files","state"});
  if(artifacts::detail::string(manifest.at("state"))!="pending-not-an-exact-result"||artifacts::detail::canonical(files(directory/"storage"))!=artifacts::detail::canonical(manifest.at("files")))throw std::runtime_error("invalid pending FIRE state or forward snapshot checksum");
  if(manifest.contains("backward_overlay")) {
    const auto& overlay=manifest.at("backward_overlay").as_object();artifacts::detail::keys(overlay,{"files","closed_sectors","log_sha256"});
    if(artifacts::detail::canonical(files(directory/"backward"))!=artifacts::detail::canonical(overlay.at("files")))throw std::runtime_error("backward FIRE snapshot checksum mismatch");
    auto log=fire::read_text(directory/"checkpoint-run.log");if(artifacts::detail::sha256(log)!=artifacts::detail::string(overlay.at("log_sha256")))throw std::runtime_error("backward FIRE finalization log checksum mismatch");
    auto closed=completed_backward(log);std::set<unsigned> listed;
    for(const auto& sector:overlay.at("closed_sectors").as_array()){std::uint64_t n=0;if(sector.is_uint64())n=sector.as_uint64();else if(sector.is_int64()&&sector.as_int64()>0)n=static_cast<std::uint64_t>(sector.as_int64());if(n<=1||n>99999||!closed.count(n)||!listed.insert(n).second)throw std::runtime_error("backward FIRE sector lacks unique completed-worker provenance");}
    if(listed.size()!=overlay.at("files").as_array().size())throw std::runtime_error("backward FIRE provenance/file count mismatch");
    std::size_t k=0;for(auto n:listed)if(artifacts::detail::string(overlay.at("files").as_array().at(k++).as_object().at("name"))!=sector_name(n))throw std::runtime_error("backward FIRE provenance/file mapping mismatch");
  }
}
inline json::object read_manifest(const std::filesystem::path& directory) {
  auto text=fire::read_text(directory/"checkpoint.json");auto manifest=json::parse(text).as_object();if(artifacts::detail::canonical(manifest)!=text)throw std::runtime_error("noncanonical FIRE checkpoint manifest");validate_manifest(directory,manifest);return manifest;
}
// Upgrade an older forward-only producer after it stops, without mutating its
// checkpoint or any durable cache generation. No new FIRE work is performed.
inline std::filesystem::path upgrade_backward(const std::filesystem::path& source) {
  auto manifest=read_manifest(source);if(manifest.contains("backward_overlay")||!std::filesystem::exists(source/"run.log")||completed_backward(fire::read_text(source/"run.log")).empty())return source;
  auto original_hash=artifacts::detail::sha256(artifacts::detail::canonical(manifest));
  auto pattern=(std::filesystem::temp_directory_path()/"diffexp3-fire-backward-upgrade-XXXXXX").string();std::vector<char> name(pattern.begin(),pattern.end());name.push_back(0);if(!mkdtemp(name.data()))throw std::runtime_error("cannot create backward checkpoint upgrade");std::filesystem::path destination=name.data();
  copy_snapshot_list(source/"storage",destination/"storage",manifest.at("files").as_array());attach_backward(destination,manifest,source);
  {std::ofstream out(destination/"checkpoint.json");out<<artifacts::detail::canonical(manifest);if(!out)throw std::runtime_error("cannot persist backward checkpoint upgrade");}
  {std::ofstream out(destination/"recovery.json");out<<artifacts::detail::canonical(json::object{{"original_directory",source.string()},{"original_checkpoint_sha256",original_hash},{"operation","verified finalized-backward overlay; original checkpoint unchanged"}});}
  validate_manifest(destination,manifest);return destination;
}
// Migrate only forward sectors whose worker finished after closing its snapshot.
// The original raw database is never modified; active/request-only files are omitted.
inline std::filesystem::path recover_forward_run(const ibp::PropagatorBasis& b,const Exact& d,
    const ExactField& field,const std::vector<ibp::Integral>& requests,const fire::Options& options,
    const std::map<ibp::Integral,ibp::Relation>* prior,const std::filesystem::path& original) {
  auto bound=identity(b,d,field,requests,options,prior);fire::SymbolMap symbols(field.variables());
  if(fire::read_text(original/"family.start")!=fire::start(b,d,options.zero_sectors,&symbols))throw std::runtime_error("raw checkpoint scientific input mismatch");
  std::string demand="{";for(std::size_t i=0;i<requests.size();++i){if(i)demand+=", ";demand+="{1, "+fire::indices(requests[i])+"}";}demand+="}\n";
  if(fire::read_text(original/"integrals.m")!=demand)throw std::runtime_error("raw checkpoint ordered demand mismatch");
  auto rules=prior?fire::rules_text(*prior,d,symbols):"";if(!rules.empty()&&fire::read_text(original/"prior.rules")!=rules)throw std::runtime_error("raw checkpoint prior rules mismatch");
  std::ostringstream config;config<<"#threads "<<options.threads<<"\n#fthreads "<<options.simplifier_threads<<"\n#compressor none\n#variables ";
  for(std::size_t i=0;i<symbols.exported.size();++i){if(i)config<<",";config<<symbols.exported[i];}
  config<<"\n#start\n#folder "<<original.string()<<"/\n#problem 1 |"<<b.physical_count<<"|family.start\n";
  if(!rules.empty())config<<"#rules prior.rules\n";config<<"#integrals integrals.m\n#output "<<(original/"result.tables").string()<<"\n";
  if(config.str()!=fire::read_text(original/"job.config"))throw std::runtime_error("raw checkpoint native config mismatch");
  for(const auto& executable:{options.executable,options.executable.parent_path()/"FLAME7"})
    if(std::filesystem::last_write_time(executable)>std::filesystem::last_write_time(original/"job.config"))throw std::runtime_error("FIRE binary changed since raw checkpoint creation");
  auto log=fire::read_text(original/"run.log");
  if(log.find("Substituting in")!=std::string::npos)throw std::runtime_error("raw checkpoint recovery requires an interrupted forward pass");
  if(log.find("Path: "+options.executable.parent_path().string()+"/") ==std::string::npos)throw std::runtime_error("raw checkpoint executable provenance mismatch");
  std::set<unsigned> completed;std::regex finished(R"(FLAME time \(([0-9]+)\):)");
  for(auto i=std::sregex_iterator(log.begin(),log.end(),finished);i!=std::sregex_iterator();++i){auto sector=std::stoul((*i)[1]);if(sector<=1||sector>99999||!completed.insert(sector).second)throw std::runtime_error("ambiguous raw completed sector provenance");}
  if(completed.empty())throw std::runtime_error("raw checkpoint has no finalized forward sectors");
  auto pattern=(std::filesystem::temp_directory_path()/"diffexp3-fire-checkpoint-XXXXXX").string();std::vector<char> name(pattern.begin(),pattern.end());name.push_back(0);
  if(!mkdtemp(name.data()))throw std::runtime_error("cannot create private checkpoint recovery directory");std::filesystem::path recovered=name.data();std::filesystem::create_directory(recovered/"storage");
  for(auto sector:completed){std::ostringstream name;name<<std::setw(5)<<std::setfill('0')<<sector<<".tmp";auto bytes=validated_snapshot(original/"database"/name.str());std::ofstream out(recovered/"storage"/name.str(),std::ios::binary);out<<bytes;if(!out)throw std::runtime_error("cannot copy finalized sector snapshot");}
  json::object manifest{{"identity",bound},{"files",files(recovered/"storage")},{"state","pending-not-an-exact-result"}};
  {std::ofstream out(recovered/"checkpoint.json");out<<artifacts::detail::canonical(manifest);if(!out)throw std::runtime_error("cannot persist recovered checkpoint");}
  {std::ofstream out(recovered/"recovery.json");out<<artifacts::detail::canonical(json::object{{"original_directory",original.string()},{"log_sha256",artifacts::detail::sha256(log)},{"completed_forward_sectors",static_cast<std::int64_t>(completed.size())},{"evidence","exact native input/config/demand/rules; binaries older than invocation; completed FLAME worker; structurally complete snapshot"}});}
  return recovered;
}
// Pending snapshots are deliberately separate from certified completed batches.
// Resource limits may change; exact inputs, ordered demands, prior rules, binary,
// and reduction worker configuration must match the saved checkpoint exactly.
inline fire::Result reduce(const ibp::PropagatorBasis& b,const Exact& d,const ExactField& field,
    const std::vector<ibp::Integral>& requests,const fire::Options& options,
    const std::map<ibp::Integral,ibp::Relation>* prior=nullptr,
    const std::filesystem::path& resume_directory={}) {
  fire::Result result;std::filesystem::path overlay_directory;json::object resume_info;
  try {
    auto expected=identity(b,d,field,requests,options,prior);std::filesystem::path storage;
    if(!resume_directory.empty()) {
      auto saved=read_manifest(resume_directory);
      if(artifacts::detail::canonical(saved.at("identity"))!=artifacts::detail::canonical(expected))throw std::runtime_error("FIRE checkpoint exact demand/config identity mismatch");
      storage=resume_directory/"storage";resume_info={{"source",resume_directory.string()},{"manifest_sha256",artifacts::detail::sha256(artifacts::detail::canonical(saved))},{"backward_sector_snapshots",0}};
      if(saved.contains("backward_overlay")) {
        auto pattern=(std::filesystem::temp_directory_path()/"diffexp3-fire-overlay-XXXXXX").string();std::vector<char> name(pattern.begin(),pattern.end());name.push_back(0);if(!mkdtemp(name.data()))throw std::runtime_error("cannot create private backward overlay");overlay_directory=name.data();
        copy_snapshot_list(storage,overlay_directory,saved.at("files").as_array());
        copy_snapshot_list(resume_directory/"backward",overlay_directory,saved.at("backward_overlay").as_object().at("files").as_array(),true);storage=overlay_directory;resume_info["backward_sector_snapshots"]=saved.at("backward_overlay").as_object().at("files").as_array().size();
      }
    }
    result=fire::reduce(b,d,field,requests,options,prior,storage.empty()?nullptr:&storage,true);
    if(!overlay_directory.empty()){std::filesystem::remove_all(overlay_directory);overlay_directory.clear();}
    if(!result.directory.empty()&&!resume_info.empty()){std::ofstream proof(result.directory/"checkpoint-resume.json");proof<<artifacts::detail::canonical(resume_info);if(!proof)throw std::runtime_error("cannot persist checkpoint replay provenance");}
    // A successful solve goes directly to the exact batch store. Capturing more
    // pending work first would waste I/O and could hide a valid completed result.
    if(!result.success&&!result.directory.empty()&&std::filesystem::is_directory(result.directory/"storage")) {
      json::object manifest{{"identity",expected},{"files",files(result.directory/"storage")},{"state","pending-not-an-exact-result"}};
      try{attach_backward(result.directory,manifest);}catch(const std::exception& e){result.reason+="; backward checkpoint unavailable: "+std::string(e.what());}
      std::ofstream out(result.directory/"checkpoint.json");out<<artifacts::detail::canonical(manifest);if(!out)throw std::runtime_error("cannot persist pending FIRE checkpoint manifest");
    }
    return result;
  }catch(const std::exception& e){if(!overlay_directory.empty()){std::error_code ignored;std::filesystem::remove_all(overlay_directory,ignored);}if(!result.reason.empty())result.reason+="; ";result.reason+="checkpoint handling: "+std::string(e.what());return result;}
}
// Durable pending work uses its own directory tree, never the exact artifact store.
class Store {
 public:
  explicit Store(std::filesystem::path root):root_(std::move(root)){std::filesystem::create_directories(root_);}
  std::optional<std::filesystem::path> lookup(const json::object& expected)const {
    const auto semantic=artifacts::detail::sha256(artifacts::detail::canonical(expected));auto parent=root_/semantic;
    if(!std::filesystem::exists(parent))return std::nullopt;
    std::optional<std::filesystem::path> best;unsigned count=0;
    for(const auto& entry:std::filesystem::directory_iterator(parent)) {
      if(entry.path().filename().string().starts_with(".pending-"))continue;
      if(++count>64)throw std::runtime_error("pending FIRE generation budget");
      if(entry.is_symlink()||!entry.is_directory())throw std::runtime_error("invalid pending FIRE generation");
      auto text=fire::read_text(entry.path()/"checkpoint.json");auto manifest=read_manifest(entry.path());
      if(artifacts::detail::sha256(text)!=entry.path().filename().string()||artifacts::detail::canonical(manifest.at("identity"))!=artifacts::detail::canonical(expected))throw std::runtime_error("pending FIRE checkpoint identity or snapshot mismatch");
      if(!best||std::filesystem::last_write_time(entry.path())>std::filesystem::last_write_time(*best))best=entry.path();
    }
    return best;
  }
  std::filesystem::path put(const std::filesystem::path& source) {
    auto text=fire::read_text(source/"checkpoint.json");auto manifest=read_manifest(source);
    auto parent=root_/artifacts::detail::sha256(artifacts::detail::canonical(manifest.at("identity")));std::filesystem::create_directories(parent);
    auto target=parent/artifacts::detail::sha256(text);
    if(std::filesystem::exists(target)){read_manifest(target);if(fire::read_text(target/"checkpoint.json")!=text)throw std::runtime_error("pending FIRE immutable collision");return target;}
    auto pattern=(parent/".pending-XXXXXX").string();std::vector<char> name(pattern.begin(),pattern.end());name.push_back(0);if(!mkdtemp(name.data()))throw std::runtime_error("pending FIRE staging directory");
    std::filesystem::path stage=name.data();std::filesystem::create_directory(stage/"storage");
    copy_snapshot_list(source/"storage",stage/"storage",manifest.at("files").as_array());
    if(manifest.contains("backward_overlay")){copy_snapshot_list(source/"backward",stage/"backward",manifest.at("backward_overlay").as_object().at("files").as_array());std::filesystem::copy_file(source/"checkpoint-run.log",stage/"checkpoint-run.log");}
    {std::ofstream out(stage/"checkpoint.json");out<<text;if(!out)throw std::runtime_error("pending FIRE staging write");}
    if(std::filesystem::exists(source/"recovery.json"))std::filesystem::copy_file(source/"recovery.json",stage/"recovery.json");
    std::error_code error;std::filesystem::rename(stage,target,error);
    if(error){if(!std::filesystem::exists(target))throw std::runtime_error("pending FIRE publication failed: "+error.message());std::filesystem::remove_all(stage);}
    return target;
  }
  void retire(const json::object& expected) {
    auto parent=root_/artifacts::detail::sha256(artifacts::detail::canonical(expected));if(!std::filesystem::exists(parent))return;
    if(std::filesystem::is_symlink(parent)||!std::filesystem::is_directory(parent))throw std::runtime_error("invalid pending FIRE retirement identity");
    auto pattern=(root_/".retired-XXXXXX").string();std::vector<char> name(pattern.begin(),pattern.end());name.push_back(0);if(!mkdtemp(name.data()))throw std::runtime_error("pending FIRE retirement staging");
    std::filesystem::path retired=name.data();std::filesystem::remove(retired);std::error_code error;std::filesystem::rename(parent,retired,error);
    if(error){if(!std::filesystem::exists(parent))return;throw std::runtime_error("pending FIRE retirement failed: "+error.message());}
    std::filesystem::remove_all(retired);
  }
 private:std::filesystem::path root_;
};
} // namespace diffexp::fire_checkpoint
