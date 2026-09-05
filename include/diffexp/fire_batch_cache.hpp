#pragma once
#include "diffexp/fire.hpp"
#include "diffexp/artifact_store.hpp"

namespace diffexp::fire_batch {
namespace json=boost::json;
inline json::array index_json(const ibp::Integral& a){json::array out;for(int n:a)out.push_back(n);return out;}
inline void validate_index(const ibp::Integral& a,const ibp::PropagatorBasis& b) {
  if(a.size()!=b.denominators.size())throw std::invalid_argument("cached FIRE index arity");
  for(std::size_t k=0;k<a.size();++k)if(a[k]<-100||a[k]>100||(k>=b.physical_count&&a[k]>0))throw std::invalid_argument("cached FIRE index range or numerator semantics");
}
inline json::array mapping_json(const std::map<ibp::Integral,ibp::Relation>& reductions) {
  json::array rows;
  for(const auto& [a,row]:reductions){json::array terms;for(const auto& [master,c]:row)terms.push_back(json::object{{"integral",index_json(master)},{"coefficient",c.str()}});rows.push_back(json::object{{"integral",index_json(a)},{"terms",terms}});}
  return rows;
}
inline void validate(const fire::Result& result,const ibp::PropagatorBasis& basis,const Exact& sample,
    const std::vector<ibp::Integral>& requested) {
  if(!result.success||requested.empty()||result.reductions.size()>20000)throw std::invalid_argument("invalid completed FIRE batch");
  for(const auto& a:requested){validate_index(a,basis);if(!result.reductions.count(a))throw std::invalid_argument("cached FIRE batch omits demand");}
  for(const auto& [a,row]:result.reductions) {
    validate_index(a,basis);
    for(const auto& [master,c]:row) {
      validate_index(master,basis);if(c.is_zero())throw std::invalid_argument("zero cached FIRE mapping term");
      auto found=result.reductions.find(master);
      if(found==result.reductions.end()||found->second.size()!=1||found->second.begin()->first!=master||found->second.begin()->second!=sample.constant(1))
        throw std::invalid_argument("cached FIRE mapping is not a complete terminal-master reduction");
    }
  }
}
inline artifacts::Identity identity(const ibp::PropagatorBasis& b,const Exact& dimension,const ExactField& field,
    const std::vector<ibp::Integral>& requested,const fire::Options& options) {
  (void)(dimension+Exact(field));if(requested.empty()||requested.size()>10000)throw std::invalid_argument("invalid FIRE batch demand");
  json::array gram,pairs,denominators,rewrites,symbols,demands,zeros;
  for(const auto& row:b.space.external_gram){json::array values;for(const auto& c:row)values.emplace_back(c.str());gram.push_back(values);}
  for(const auto& [a,c]:b.space.pairs)pairs.push_back(json::array{a,c});
  auto affine=[&](const ibp::Affine& a){json::array linear;for(const auto& c:a.linear)linear.emplace_back(c.str());return json::object{{"constant",a.constant.str()},{"linear",linear}};};
  for(const auto& a:b.denominators)denominators.push_back(affine(a));for(const auto& a:b.scalar_products_in_denominators)rewrites.push_back(affine(a));
  for(const auto& name:field.variables())symbols.emplace_back(name);
  for(const auto& a:requested){validate_index(a,b);demands.push_back(index_json(a));}
  auto sectors=options.zero_sectors;std::sort(sectors.begin(),sectors.end());sectors.erase(std::unique(sectors.begin(),sectors.end()),sectors.end());for(const auto& a:sectors)zeros.push_back(index_json(a));
  artifacts::Identity id;id.kind="exact_equation";id.algorithm_version="native-fire-completed-batch-v1";
  id.family={{"loops",b.space.loops},{"physical_count",static_cast<std::int64_t>(b.physical_count)},{"denominators",denominators},{"rewrites",rewrites}};
  id.ordered_basis={json::object{{"ordered_demands",demands}}};
  id.normalization={{"measure","native scalar integral"},{"propagators","supplied exact affine polynomials"}};
  id.branch={{"status","formal rational IBP identities"}};id.boundary={{"status","not-applicable to exact reduction"}};
  id.geometry={{"external_gram",gram},{"scalar_product_pairs",pairs}};
  id.scientific_inputs={{"dimension",dimension.str()},{"ordered_field_symbols",symbols},{"proven_zero_sectors",zeros},{"native_start",fire::start(b,dimension,sectors)}};
  return id;
}
class Cache {
 public:
  explicit Cache(std::filesystem::path directory):store_(std::move(directory)) {}
  std::optional<fire::Result> lookup(const ibp::PropagatorBasis& b,const Exact& d,const ExactField& field,
      const std::vector<ibp::Integral>& requested,const fire::Options& options)const {
    auto id=identity(b,d,field,requested,options);auto record=store_.lookup(id,resource(),{"exact",verifier(),"completed_provider_batch"});if(!record)return std::nullopt;
    const auto& payload=record->payload.as_object();artifacts::detail::keys(payload,{"schema","start","tables","mappings","directory","completion"});
    if(artifacts::detail::string(payload.at("schema"))!="DiffExp3.CompletedFIREBatch/v1"||artifacts::detail::string(payload.at("completion"))!="validated-success")throw std::runtime_error("cached FIRE completion schema");
    fire::SymbolMap symbols(field.variables());
    if(artifacts::detail::string(payload.at("start"))!=fire::start(b,d,options.zero_sectors,&symbols))throw std::runtime_error("cached FIRE native input mismatch");
    fire::Result out;out.success=true;out.directory=artifacts::detail::string(payload.at("directory"));
    out.reductions=fire::import_table_text(artifacts::detail::string(payload.at("tables")),field,b.denominators.size(),&symbols);
    validate(out,b,d,requested);
    if(artifacts::detail::canonical(mapping_json(out.reductions))!=artifacts::detail::canonical(payload.at("mappings")))throw std::runtime_error("cached FIRE exact mappings disagree with original table identities");
    return out;
  }
  void put(const ibp::PropagatorBasis& b,const Exact& d,const ExactField& field,
      const std::vector<ibp::Integral>& requested,const fire::Options& options,const fire::Result& result) {
    validate(result,b,d,requested);
    if(fire::read_text(result.directory/"success.receipt",128)!="DiffExp3.NativeFIRECompleted/v1\n")throw std::runtime_error("FIRE batch has no validated success receipt");
    auto text=fire::read_text(result.directory/"result.tables"),start=fire::read_text(result.directory/"family.start");fire::SymbolMap symbols(field.variables());
    if(start!=fire::start(b,d,options.zero_sectors,&symbols)||fire::import_table_text(text,field,b.denominators.size(),&symbols)!=result.reductions)throw std::runtime_error("FIRE receipt input/table mismatch");
    auto id=identity(b,d,field,requested,options);
    json::object payload{{"schema","DiffExp3.CompletedFIREBatch/v1"},{"start",start},{"tables",text},{"mappings",mapping_json(result.reductions)},{"directory",result.directory.string()},{"completion","validated-success"}};
    artifacts::Certificate certificate{"exact",verifier(),"completed_provider_batch",{{"assumptions","FIRE-produced IBP reductions and declared zero sectors"},{"check","successful exit and import receipt; input equality; exact table re-import; complete terminal-master mappings"}}};
    if(std::filesystem::exists(result.directory/"recovery.json"))certificate.evidence["replay_provenance"]=json::parse(fire::read_text(result.directory/"recovery.json",4096));
    store_.put(id,resource(),payload,certificate);
  }
  // Explicit migration of pre-receipt runs. A later native session could only
  // emit prior.rules after accepting the predecessor as a successful batch.
  // Bare output tables, even well-formed ones, never constitute this evidence.
  void recover_replayed(const ibp::PropagatorBasis& b,const Exact& d,const ExactField& field,
      const fire::Options& options,const std::filesystem::path& predecessor,
      const std::filesystem::path& successor) {
    fire::SymbolMap symbols(field.variables());auto expected=fire::start(b,d,options.zero_sectors,&symbols);
    if(fire::read_text(predecessor/"family.start")!=expected||fire::read_text(successor/"family.start")!=expected)
      throw std::runtime_error("replayed FIRE recovery scientific input mismatch");
    fire::Result result;result.success=true;
    auto tables=fire::read_text(predecessor/"result.tables");
    result.reductions=fire::import_table_text(tables,field,b.denominators.size(),&symbols);
    auto rules=fire::rules_text(result.reductions,d,symbols);
    if(rules.empty()||rules!=fire::read_text(successor/"prior.rules"))throw std::runtime_error("FIRE recovery has no exact successful replay evidence");
    auto parse_demands=[&](const std::filesystem::path& path){
      auto parsed=data::Reader(fire::read_text(path)).read();std::vector<ibp::Integral> requests;
      for(const auto& item:fire::list(parsed)){const auto& pair=fire::list(item,2);if(fire::table_integer(pair[0])!=1)throw std::runtime_error("FIRE recovery problem number");ibp::Integral a;for(const auto& power:fire::list(pair[1],b.denominators.size())){auto n=fire::table_integer(power);if(n< -100||n>100)throw std::runtime_error("FIRE recovery index range");a.push_back(static_cast<int>(n));}validate_index(a,b);requests.push_back(std::move(a));}return requests;
    };
    auto requests=parse_demands(predecessor/"integrals.m"),next=parse_demands(successor/"integrals.m");
    std::set<ibp::Integral> later(next.begin(),next.end());for(const auto& a:requests)if(!later.count(a))throw std::runtime_error("FIRE recovery demands are not monotone");
    validate(result,b,d,requests);
    auto pattern=(std::filesystem::temp_directory_path()/"diffexp3-fire-recovered-XXXXXX").string();std::vector<char> name(pattern.begin(),pattern.end());name.push_back(0);
    if(!mkdtemp(name.data()))throw std::runtime_error("FIRE recovery private directory");result.directory=name.data();
    auto write=[&](const char* file,const std::string& content){std::ofstream out(result.directory/file);out<<content;if(!out)throw std::runtime_error("FIRE recovery write failed");};
    write("family.start",expected);write("result.tables",tables);write("success.receipt","DiffExp3.NativeFIRECompleted/v1\n");
    write("recovery.json",artifacts::detail::canonical(json::object{{"predecessor",predecessor.string()},{"successor",successor.string()},{"evidence","exact accepted reduction replay through prior.rules"}}));
    put(b,d,field,requests,options,result);
  }
 private:
  artifacts::Store store_;
  static Demand resource(){return {0,0,0,64,0};}
  static const char* verifier(){return "native-fire-batch-reimport-v1";}
};
} // namespace diffexp::fire_batch
