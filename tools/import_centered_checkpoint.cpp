// Explicit migration of imported numerical boundary states, never implicit reuse.
#include "diffexp/adjoint_checkpoint.hpp"
#include <iostream>
#include <map>
#include <fcntl.h>
using namespace diffexp;namespace json=boost::json;namespace fs=std::filesystem;
using B=Jet::Ball;
struct Input {fs::path path;std::string bytes,sha;};
static Input read(const fs::path& path){
  if(fs::is_symlink(path)||!fs::is_regular_file(path)||fs::file_size(path)>64*1024*1024)throw std::invalid_argument("invalid snapshot file or size: "+path.string());
  std::ifstream in(path,std::ios::binary);if(!in)throw std::runtime_error("cannot read "+path.string());
  std::string bytes((std::istreambuf_iterator<char>(in)),{});return {fs::absolute(path),bytes,artifacts::detail::sha256(bytes)};
}
static bool checkpoint_name(const std::string& name){return name.size()==69&&name.substr(64)==".json"&&std::all_of(name.begin(),name.begin()+64,[](char c){return (c>='0'&&c<='9')||(c>='a'&&c<='f');});}
static void create_file(const fs::path& path,const std::string& bytes){
  int fd=::open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600);if(fd<0)throw std::runtime_error("cannot exclusively create "+path.string());
  try{std::size_t at=0;while(at<bytes.size()){auto count=::write(fd,bytes.data()+at,bytes.size()-at);if(count<0&&errno==EINTR)continue;if(count<=0)throw std::runtime_error("import write failed");at+=count;}
    if(::fsync(fd))throw std::runtime_error("import sync failed");if(::close(fd)){fd=-1;throw std::runtime_error("import close failed");}fd=-1;
  }catch(...){if(fd>=0)::close(fd);throw;}
}
int main(int argc,char** argv){try{
  if(argc<4||argc>5||(argc==5&&std::string(argv[4])!="--dry-run"))throw std::invalid_argument("usage: import_centered_checkpoint CAPTURE_FILE_OR_DIRECTORY OLD_CHECKPOINT_DIRECTORY NEW_EMPTY_DIRECTORY [--dry-run]; inputs must be immutable snapshots");
  bool dry=argc==5;fs::path captures=argv[1],source=argv[2],destination=fs::absolute(argv[3]);
  if(fs::is_symlink(captures)||fs::is_symlink(source)||fs::is_symlink(destination))throw std::invalid_argument("symlink inputs/destination forbidden");
  if(!fs::is_directory(source))throw std::invalid_argument("old checkpoint directory missing");
  if(fs::exists(destination)&&(!fs::is_directory(destination)||!fs::is_empty(destination)))throw std::invalid_argument("destination must be absent or empty; refusing overwrite");
  std::vector<Input> capture_files;std::map<std::string,Input> old_files;
  if(fs::is_regular_file(captures))capture_files.push_back(read(captures));else if(fs::is_directory(captures)){
    std::vector<fs::path> paths;for(const auto& entry:fs::directory_iterator(captures)){const auto name=entry.path().filename().string();if(name.starts_with("problem-")&&name.ends_with(".json"))paths.push_back(entry.path());}
    std::sort(paths.begin(),paths.end());for(const auto& path:paths)capture_files.push_back(read(path));
  }else throw std::invalid_argument("capture file or directory missing");
  for(const auto& entry:fs::directory_iterator(source)){auto name=entry.path().filename().string();
    if(name.starts_with("._"))continue; // ExFAT filesystem metadata, never records.
    if(!checkpoint_name(name))throw std::invalid_argument("unexpected source checkpoint filename: "+name);
    old_files.emplace(name.substr(0,64),read(entry.path()));}
  if(capture_files.empty()||old_files.empty())throw std::invalid_argument("no captures/checkpoints found");
  std::map<std::string,std::string> outputs;std::set<std::string> matched;json::array imports;
  for(const auto& captured:capture_files){
    auto object=json::parse(captured.bytes).as_object();
    if(object.at("schema")!="DiffExp.FTTransportProblem/v1")throw std::invalid_argument("capture schema mismatch");
    auto precision=artifacts::detail::integer(object.at("working_bits"));if(precision<128||precision>4096)throw std::invalid_argument("captured precision budget");B::set_precision(precision);
    std::vector<std::string> names;for(const auto& n:object.at("variables").as_array())names.emplace_back(n.as_string());if(names.empty()||names.size()>16)throw std::invalid_argument("field symbol budget");ExactField field(names);
    auto matrix=[&](const json::value& value){ExactEpsilonMatrix out;for(const auto& row:value.as_array()){out.emplace_back();for(const auto& item:row.as_array())out.back().emplace_back(field,std::string(item.as_string()));}return out;};
    auto connection=matrix(object.at("matrix")),forcing=matrix(object.at("forcing"));auto initial=numerical_rows_io::read_rows(object.at("initial"));
    auto d=initial.columns(),r=initial.coefficients.size();if(connection.size()!=d||forcing.size()!=r)throw std::invalid_argument("capture matrix/source dimension mismatch");
    for(const auto* m:{&connection,&forcing})for(const auto& row:*m)if(row.size()!=d)throw std::invalid_argument("capture matrix width mismatch");
    std::vector<Exact> path;for(const auto& point:object.at("path").as_array())path.emplace_back(field,std::string(point.as_string()));if(path.size()<2||path.size()>1000)throw std::invalid_argument("capture path size");
    auto [xi,ei]=path_epsilon_variables(path.front());(void)xi;
    const auto& settings=object.at("settings").as_object();AdjointOptions old_options;old_options.compact_centered_recovery=false;old_options.compact_centered_only=false;old_options.circuit_grouped_dot=false;
    const auto number=[&](const char* key,std::int64_t fallback,std::int64_t low,std::int64_t high){auto n=settings.contains(key)?artifacts::detail::integer(settings.at(key)):fallback;if(n<low||n>high)throw std::invalid_argument(std::string("captured option range: ")+key);return n;};
    if(number("working_bits",precision,128,4096)!=precision)throw std::invalid_argument("settings/captured precision disagree");
    old_options.taylor_order=number("ordinary_order",80,8,1000);
    old_options.centered_before_rational=settings.contains("centered_before_rational")&&settings.at("centered_before_rational").as_bool();
    old_options.centered_map_working_bits=number("centered_map_working_bits",0,0,4096);
    if(old_options.centered_map_working_bits&&old_options.centered_map_working_bits<64)throw std::invalid_argument("map precision range");
    old_options.centered_map_workers=number("centered_map_workers",1,1,4);
    old_options.rational_circuit_recurrence=!settings.contains("rational_circuit_recurrence")||settings.at("rational_circuit_recurrence").as_bool();
    if(!old_options.rational_circuit_recurrence)throw std::invalid_argument("capture is not an original v6 circuit problem");
    for(const char* flag:{"compact_centered_recovery","compact_centered_only","circuit_grouped_dot"})
      if(settings.contains(flag)&&settings.at(flag).as_bool())throw std::invalid_argument(std::string("capture enables a post-v6 policy: ")+flag);
    auto old_key=adjoint_checkpoint::identity(connection,initial,forcing,path,old_options);auto found=old_files.find(old_key);
    if(found==old_files.end())continue; // An unfinished capture may have no accepted checkpoint yet.
    if(!matched.insert(old_key).second)throw std::invalid_argument("multiple captures match one source checkpoint; provide an unambiguous snapshot");
    auto envelope=json::parse(found->second.bytes);auto saved=adjoint_checkpoint::detail::decode(envelope,old_key);
    int low=std::min(0,initial.low);for(const auto& row:forcing)for(const auto& e:row)if(!e.is_zero()){auto valuation=exact_epsilon_valuation(e,ei);if(!valuation)throw std::invalid_argument("forcing epsilon valuation unavailable");low=std::min(low,static_cast<int>(*valuation));}
    if(saved.rows.low!=low||saved.rows.high!=initial.high||saved.rows.columns()!=d||saved.rows.coefficients.size()!=r||saved.leg+1>=path.size()||path[saved.leg]==path[saved.leg+1]||!std::isfinite(saved.parameter)||saved.parameter<=0||saved.parameter>1||!saved.accepted_charts||saved.accepted_charts>old_options.max_charts_per_leg)throw std::invalid_argument("source continuation shape, window, geometry or chart budget mismatch");
    // Includes the original IEEE parameter bits, counts and every exact Arb dump.
    if(adjoint_checkpoint::detail::encode(saved,old_key)!=envelope.at("payload"))throw std::invalid_argument("source continuation does not round-trip unchanged");
    auto new_options=old_options;new_options.compact_centered_recovery=true;new_options.compact_centered_only=true;new_options.circuit_grouped_dot=true;
    auto new_key=adjoint_checkpoint::identity(connection,initial,forcing,path,new_options);if(old_key==new_key)throw std::logic_error("old/new identities failed to separate policies");
    auto payload=adjoint_checkpoint::detail::encode(saved,new_key);auto numerical=payload;numerical["identity"]=old_key;
    if(numerical!=envelope.at("payload"))throw std::logic_error("import modified numeric continuation data");
    auto new_bytes=artifacts::detail::canonical(json::object{{"payload",payload},{"sha256",artifacts::detail::sha256(artifacts::detail::canonical(payload))}})+"\n";
    if(!outputs.emplace(new_key,new_bytes).second)throw std::invalid_argument("duplicate destination identity");
    imports.push_back(json::object{{"capture",captured.path.string()},{"capture_sha256",captured.sha},{"captured_settings_sha256",artifacts::detail::sha256(artifacts::detail::canonical(settings))},
      {"source_checkpoint",found->second.path.string()},{"source_checkpoint_sha256",found->second.sha},{"old_identity",old_key},{"new_identity",new_key},{"destination_filename",new_key+".json"},{"destination_sha256",artifacts::detail::sha256(new_bytes)},
      {"unchanged_rows_sha256",artifacts::detail::sha256(artifacts::detail::canonical(payload.at("rows")))},{"leg",saved.leg},{"parameter_ieee_bits",payload.at("parameter_ieee_bits")},{"accepted_charts",saved.accepted_charts},{"working_bits",precision},{"order",old_options.taylor_order},{"epsilon_low",saved.rows.low},{"epsilon_high",saved.rows.high},{"source_rows",r},{"physical_dimension",d}});
  }
  if(matched.size()!=old_files.size()){std::string missing;for(const auto& [key,file]:old_files)if(!matched.contains(key))missing+=" "+key;throw std::invalid_argument("unmatched source checkpoints; capture/settings mismatch:"+missing);}
  // Detect changing source snapshots before publishing anything.
  for(const auto& captured:capture_files)if(read(captured.path).sha!=captured.sha)throw std::runtime_error("capture changed during import");
  for(const auto& [key,file]:old_files)if(read(file.path).sha!=file.sha)throw std::runtime_error("checkpoint changed during import");
  json::object manifest{{"schema","DiffExp21.ExplicitCenteredCheckpointImport/v1"},{"source_algorithm","DiffExp.AdjointOriginalPath/v6-circuit"},{"destination_algorithm","DiffExp.AdjointOriginalPath/v8-centered-dot"},
    {"changed_policy_flags",json::object{{"compact_centered_recovery",json::object{{"old",false},{"new",true}}},{"compact_centered_only",json::object{{"old",false},{"new",true}}},{"circuit_grouped_dot",json::object{{"old",false},{"new",true}}}}},{"all_other_options_unchanged",true},{"numeric_states_unchanged",true},{"old_charts_recomputed",false},{"analytic_tail_certificate",false},
    {"meaning","Explicitly imported retained numerical boundary states; old charts were computed with v6, not v8."},{"flint_version",FLINT_VERSION},{"destination",destination.string()},{"imports",imports},{"dry_run",dry}};
  if(!dry){
    if(!fs::is_directory(destination.parent_path()))throw std::invalid_argument("destination parent must already exist");
    auto pattern=destination.string()+".staging.XXXXXX";std::vector<char> name(pattern.begin(),pattern.end());name.push_back(0);if(!::mkdtemp(name.data()))throw std::runtime_error("cannot create import staging directory");fs::path stage=name.data();
    try{for(const auto& [key,content]:outputs)create_file(stage/(key+".json"),content);create_file(stage/"import-manifest.json",artifacts::detail::canonical(manifest)+"\n");artifacts::detail::sync_directory(stage);
      if(fs::is_symlink(destination)||(fs::exists(destination)&&(!fs::is_directory(destination)||!fs::is_empty(destination))))throw std::invalid_argument("destination became nonempty; refusing overwrite");
      fs::rename(stage,destination);artifacts::detail::sync_directory(destination.parent_path());
    }catch(...){std::error_code ignored;fs::remove_all(stage,ignored);throw;}
  }
  std::cout<<json::serialize(manifest)<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
