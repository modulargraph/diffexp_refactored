#include "diffexp/cached_affine.hpp"
#include <iostream>
using namespace diffexp;
namespace ca=diffexp::cached_affine;
namespace ar=diffexp::artifacts;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* why){try{f();}catch(const std::exception&){return;}throw std::runtime_error(why);}
template<class F>void rejects_with(F f,const std::string& text,const char* why){try{f();}catch(const std::exception& e){require(std::string(e.what()).find(text)!=std::string::npos,why);return;}throw std::runtime_error(why);}
void same(const AffineFrobeniusSeries& a,const AffineFrobeniusSeries& b){
  require(a.dimension()==b.dimension() && a.x_order()==b.x_order(),"restored shape mismatch");
  require(ca::detail::payload(a)==ca::detail::payload(b),"restored exact payload mismatch");
}
int main(){try{
  auto path=ar::fs::temp_directory_path()/("de3-cached-affine-"+std::to_string(::getpid()));
  struct Cleanup{ar::fs::path p;~Cleanup(){std::error_code e;ar::fs::remove_all(p,e);}}cleanup{path};
  ar::Store store(path/"main");ExactField field({"x","eps"});auto e=[&](const char* s){return Exact(field,s);};
  AffineFrobeniusSeries::Options options;
  AffineFrobeniusSeries::Matrix a{{e("eps/x"),e("1/x")},{e("0"),e("eps/x+1/(1-x)")}};
  auto fresh=ca::prepare_legacy(a,0,1,8,options,store);
  require(!fresh.cache_hit && !fresh.series.terms().omitted_tail_certified,"first preparation or tail claim");
  auto hit=ca::prepare_legacy(a,0,1,8,options,store);require(hit.cache_hit && hit.content_id==fresh.content_id,"verified cache not reused");same(fresh.series,hit.series);
  auto key=ca::identity(a,0,1,8,options);auto record=store.read(key,fresh.content_id);
  require(record.certificate.type=="exact" && record.certificate.scope==ca::detail::scope,"incorrect retained-only certificate");
  ExactField second({"x","eps"});auto make2=[&](const char* s){return Exact(second,s);};
  AffineFrobeniusSeries::Matrix a2{{make2("eps/x"),make2("1/x")},{make2("0"),make2("eps/x+1/(1-x)")}};
  ar::Store reopened(path/"main");auto changed_options=options;changed_options.max_epsilon_depth=800;changed_options.finite_lag_recurrence=false;
  Jet::Ball::set_precision(512);auto again=ca::prepare_legacy(a2,0,1,8,changed_options,reopened);
  require(again.cache_hit && again.semantic_id==fresh.semantic_id,"fresh exact field or receiving resources prevented reuse");
  require(again.series.residue_frame()[0][0]+make2("1")==make2("2"),"restored coefficients do not belong to caller field");
  same(fresh.series,again.series);
  auto small=options;small.max_dimension=1;rejects([&]{ca::prepare_legacy(a,0,1,8,small,store);},"receiving dimension cap bypassed");
  small=options;small.max_terms=1;rejects([&]{ca::prepare_legacy(a,0,1,8,small,store);},"receiving term cap bypassed");
  ca::VerificationLimits work;work.max_term_products=1;rejects([&]{ca::prepare_legacy(a,0,1,8,options,store,work);},"receiving verification work cap bypassed");
  ca::VerificationLimits storage;storage.max_total_column_terms=1;
  rejects([&]{ca::prepare_legacy(a,0,1,8,options,store,storage);},"receiving column storage cap bypassed");
  storage=ca::VerificationLimits{};storage.max_shared_terms=1;
  rejects([&]{ca::prepare_legacy(a,0,1,8,options,store,storage);},"receiving shared storage cap bypassed");
  auto changed=a;changed[1][1]=e("eps/x+2/(1-x)");require(!ca::prepare_legacy(changed,0,1,8,options,store).cache_hit,"changed matrix reused");
  require(!ca::prepare_legacy(a,0,1,4,options,store).cache_hit,"changed retained order reused");
  ExactField reversed({"eps","x"});auto r=[&](const char* s){return Exact(reversed,s);};
  AffineFrobeniusSeries::Matrix reverse{{r("eps/x"),r("1/x")},{r("0"),r("eps/x+1/(1-x)")}};
  require(!ca::prepare_legacy(reverse,1,0,8,options,store).cache_hit,"ordered variables/indices reused");
  unsigned id=0;auto tamper=[&](auto change,const char* why){
    ar::Store corrupt(path/("bad-"+std::to_string(++id)));auto payload=record.payload;change(payload.as_object());
    corrupt.put(key,record.guarantee,payload,record.certificate);
    rejects([&]{ca::prepare_legacy(a,0,1,8,options,corrupt);},why);
  };
  tamper([](auto& p){p.at("terms").as_array().back().as_object()["coefficient"]="123";},"valid hashes hid incorrect retained coefficient");
  tamper([](auto& p){p.at("terms").as_array().erase(p.at("terms").as_array().begin());},"missing seeded solution accepted");
  tamper([](auto& p){p.at("terms").as_array().push_back(p.at("terms").as_array().front());},"duplicate term accepted");
  tamper([](auto& p){p.at("frame").as_array()[0].as_array()[0]="2";},"noncanonical frame accepted");
  tamper([](auto& p){p.at("frontiers").as_array()[0]="9";},"wrong frontier accepted");
  tamper([](auto& p){p.at("wronskian_prefactor")="2";},"wrong Wronskian accepted");
  tamper([](auto& p){p.at("omitted_tail_certified")=true;},"omitted tail claim accepted");
  tamper([](auto& p){p.at("terms").as_array()[0].as_object()["coefficient"]="x";},"non-epsilon coefficient accepted");
  // An altered homogeneous normalization still solves the ODE. The separate
  // resonance checks must reject it even when its polynomial residual vanishes.
  ar::Store scalar_store(path/"scalar");AffineFrobeniusSeries::Matrix scalar{{e("0")}};
  auto scalar_first=ca::prepare_legacy(scalar,0,1,4,options,scalar_store);auto scalar_key=ca::identity(scalar,0,1,4,options);
  auto scalar_record=scalar_store.read(scalar_key,scalar_first.content_id);auto scalar_payload=scalar_record.payload;
  scalar_payload.as_object().at("terms").as_array()[0].as_object()["coefficient"]="2";
  ar::Store wrong_constant(path/"constant");wrong_constant.put(scalar_key,scalar_record.guarantee,scalar_payload,scalar_record.certificate);
  rejects_with([&]{ca::prepare_legacy(scalar,0,1,4,options,wrong_constant);},"resonant integration constant","ODE-compatible wrong integration constant was not rejected by normalization");
  // Explicitly unverified checkpoint recovery still performs the full check.
  ar::Store checkpoint(path/"checkpoint");checkpoint.put(key,record.guarantee,record.payload);
  auto recovered=ca::prepare_legacy(a,0,1,8,options,checkpoint);require(recovered.cache_hit,"durable unverified checkpoint not recovered");same(fresh.series,recovered.series);
  require(checkpoint.read(key,recovered.content_id).certificate.type=="exact","checkpoint not promoted after verification");
  ar::Store interrupted(path/"interrupted-verifier");
  rejects([&]{ca::prepare_legacy(a,0,1,8,options,interrupted,work);},"producer verification budget did not fail");
  auto pending=interrupted.lookup(key,record.guarantee);
  require(pending && pending->certificate.type=="uncertified_numerical","expensive preparation was not durably checkpointed before verification");
  require(ca::prepare_legacy(a,0,1,8,options,interrupted).cache_hit,"interrupted verification recomputed exact preparation");
  // Integer resonance and moving epsilon poles exercise normalization/frontier
  // and the absence of a jointly analytic Wronskian unit certificate.
  AffineFrobeniusSeries::Matrix resonance{{e("0"),e("0")},{e("1/(1-x)"),e("1/x")}};
  auto resonant=ca::prepare_legacy(resonance,0,1,8,options,store);
  same(resonant.series,ca::prepare_legacy(resonance,0,1,8,options,store).series);
  auto resonance_key=ca::identity(resonance,0,1,8,options);
  auto resonance_record=store.read(resonance_key,resonant.content_id);auto resonance_payload=resonance_record.payload;
  resonance_payload.as_object().at("terms").as_array().push_back(boost::json::object{
      {"row",1},{"column",0},{"log_degree",0},{"power","1"},{"slope","0"},{"coefficient","1"}});
  ar::Store wrong_resonance(path/"positive-resonance");wrong_resonance.put(resonance_key,resonance_record.guarantee,resonance_payload,resonance_record.certificate);
  rejects_with([&]{ca::prepare_legacy(resonance,0,1,8,options,wrong_resonance);},"resonant integration constant","ODE-compatible positive-order resonance was not rejected by normalization");
  AffineFrobeniusSeries::Matrix moving{{e("eps/x+1/(x+eps)")}};
  auto moving_first=ca::prepare_legacy(moving,0,1,4,options,store);require(!moving_first.series.terms().wronskian_prefactor,"moving-pole Wronskian falsely certified");
  require(ca::prepare_legacy(moving,0,1,4,options,store).cache_hit,"moving epsilon pole cache reuse failed");
  std::cout<<"Cached affine preparation: durable reuse, fresh fields, receiving caps, exact polynomial residuals, resonance normalization, metadata and valid-hash tampering rejection passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
