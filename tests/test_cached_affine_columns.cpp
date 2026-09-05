#include "diffexp/cached_affine.hpp"
#include <iostream>
using namespace diffexp;
namespace ca=diffexp::cached_affine;
namespace ar=diffexp::artifacts;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* why){try{f();}catch(const std::exception&){return;}throw std::runtime_error(why);}
int main(){try{
  auto path=ar::fs::temp_directory_path()/("de3-affine-columns-"+std::to_string(::getpid()));
  struct Cleanup{ar::fs::path p;~Cleanup(){std::error_code e;ar::fs::remove_all(p,e);}}cleanup{path};
  ExactField field({"x","eps"});auto e=[&](const char* s){return Exact(field,s);};
  ca::Matrix a{{e("eps/x+1/(1-x)"),e("1/x"),e("0")},
    {e("0"),e("eps/x+1/(1-x)"),e("1/x")},{e("0"),e("0"),e("eps/x+1/(1-x)")}};
  ca::Series::Options options;options.univariate_epsilon_recurrence=true;const unsigned n=14;
  const auto base=ca::identity(a,0,1,n,options),key=ca::manifest_identity(base,3);const Demand demand{0,0,n,64,0};
  ar::Store store(path/"main");
  options.column_progress=[](unsigned completed,unsigned){if(completed==2)throw std::runtime_error("test interruption after two durable columns");};
  rejects([&]{ca::prepare(a,0,1,n,options,store);},"column interruption was not delivered");
  require(store.lookup(ca::column_identity(base,0),demand).has_value() && store.lookup(ca::column_identity(base,1),demand).has_value(),"completed columns not durable");
  require(!store.lookup(ca::column_identity(base,2),demand) && !store.lookup(key,demand),"incomplete work published");
  options.column_progress={};auto result=ca::prepare(a,0,1,n,options,store);
  require(result.columns_reused==2 && result.columns_prepared==1 && !result.cache_hit,"resume did not reuse exactly completed columns");
  auto direct=ca::Series::prepare(a,0,1,n,options);
  require(ca::detail::payload(result.series)==ca::detail::payload(direct),"resumed recurrence differs from direct recurrence");
  auto manifest=store.read(key,result.content_id);
  require(manifest.certificate.type=="exact" && !manifest.payload.as_object().contains("terms"),"manifest certificate/layout mismatch");
  unsigned verified_columns=0;std::size_t previous_products=0;ca::VerificationLimits progress;
  progress.column_progress=[&](unsigned c,unsigned d,std::size_t products){
    require(c==++verified_columns && d==3 && products>=previous_products,"verification progress is not cumulative/in column order");previous_products=products;
  };
  auto hit=ca::prepare(a,0,1,n,options,store,progress);
  require(verified_columns==3 && previous_products>0,"cached series did not report full verification work");
  require(hit.cache_hit && hit.columns_reused==3 && hit.columns_prepared==0 && hit.content_id==result.content_id,"verified manifest not reused");
  ExactField second({"x","eps"});ca::Matrix a2;for(const auto& row:a){a2.emplace_back();for(const auto& c:row)a2.back().emplace_back(second,c.str());}
  auto again=ca::prepare(a2,0,1,n,options,store);require(again.cache_hit,"fresh-field columns not reused");
  require(again.series.terms().terms.front().coefficient+Exact(second,1)==Exact(second,2),"column coefficients use wrong field");
  auto small=options;small.max_terms=2;rejects([&]{ca::prepare(a,0,1,n,small,store);},"receiving term cap bypassed");
  small=options;small.max_terms=direct.terms().terms.size()-1;
  rejects([&]{ca::prepare(a,0,1,n,small,store);},"cumulative column term cap bypassed");
  ca::VerificationLimits limited;limited.max_term_products=1;
  rejects([&]{ca::prepare(a,0,1,n,options,store,limited);},"verified cache bypassed independent verifier");
  ar::Store interrupted(path/"verify");rejects([&]{ca::prepare(a,0,1,n,options,interrupted,limited);},"fresh verification work cap bypassed");
  auto pending=interrupted.lookup(key,demand);require(pending && pending->certificate.type=="uncertified_numerical","unverified manifest not saved before verifier");
  auto recovered=ca::prepare(a,0,1,n,options,interrupted);require(recovered.cache_hit && recovered.columns_reused==3,"interrupted verifier recomputed columns");
  // Restrict receiving bytes to the largest column/manifest envelope: the old
  // monolithic payload cannot fit, while every chunk fits without raising caps.
  std::uintmax_t max_bytes=ar::fs::file_size(store.path(key.key(),result.content_id));
  for(unsigned c=0;c<3;++c)max_bytes=std::max(max_bytes,ar::fs::file_size(store.path(ca::column_identity(base,c).key(),manifest.parent_content_ids[c])));
  ar::Store bounded(path/"bounded",max_bytes+128);
  rejects([&]{ca::prepare_legacy(a,0,1,n,options,bounded);},"fixture does not exceed monolithic receiving byte limit");
  require(ca::prepare(a,0,1,n,options,bounded).columns_prepared==3,"column storage did not fit unchanged byte limit");
  unsigned test=0;
  auto tamper=[&](auto alter,bool publish_manifest,const char* why){
    ar::Store corrupt(path/("tamper"+std::to_string(++test)));std::vector<std::string> contents;
    for(unsigned c=0;c<3;++c){auto id=ca::column_identity(base,c);auto record=store.read(id,manifest.parent_content_ids[c]);if(c==1)alter(record.payload.as_object());contents.push_back(corrupt.put(id,demand,record.payload,record.certificate).content_id);}
    if(publish_manifest)corrupt.put(key,demand,manifest.payload,manifest.certificate,contents);
    rejects([&]{ca::prepare(a,0,1,n,options,corrupt);},why);
  };
  tamper([](auto& p){p.at("terms").as_array().back().as_object()["coefficient"]="123";},true,"valid-hash coefficient corruption passed manifest verification");
  tamper([](auto& p){p.at("terms").as_array().back().as_object()["coefficient"]="123";},false,"valid-hash unverified column corruption passed resumed verification");
  tamper([](auto& p){p["column"]=0;},true,"wrong column payload accepted");
  tamper([](auto& p){p.at("terms").as_array().front().as_object()["column"]=0;},false,"wrong column term accepted");
  tamper([](auto& p){p.at("terms").as_array().push_back(p.at("terms").as_array().front());},true,"duplicate column term accepted");
  auto first_column=ca::column_identity(base,0);auto file=store.path(first_column.key(),manifest.parent_content_ids[0]);
  ar::fs::remove(file);rejects([&]{ca::prepare(a,0,1,n,options,store);},"missing manifest parent accepted");
  // Existing complete v1 artifacts remain compatible and independently checked.
  ar::Store legacy(path/"legacy");auto old=ca::prepare_legacy(a,0,1,n,options,legacy);auto reused=ca::prepare(a,0,1,n,options,legacy);
  require(reused.cache_hit && reused.content_id==old.content_id,"legacy verified series was not reused");
  std::cout<<"Affine column cache: interrupted preparation/verifier recovery, bounded records, exact coefficient equivalence, fresh fields, legacy compatibility, missing parents and valid-hash tampering rejection passed\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
