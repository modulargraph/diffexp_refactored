#include "diffexp/level_cache.hpp"
#include <iostream>
using namespace diffexp;
namespace a=diffexp::artifacts;
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* why){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}require(rejected,why);}
ibp::PropagatorBasis tadpole(const ExactField& field,long gram=1) {
  Exact x(field,"x");ibp::ScalarProducts space(1,{{x.constant(gram)}},x);
  auto denominator=space.zero();denominator.constant=x;denominator.linear[0]=x.constant(-1);
  return ibp::PropagatorBasis({space,{denominator},{}});
}
int main(){try {
  auto root=a::fs::temp_directory_path()/("de3-level-cache-"+std::to_string(::getpid()));
  struct Cleanup{a::fs::path p;~Cleanup(){std::error_code ignored;a::fs::remove_all(p,ignored);}}cleanup{root};
  a::fs::remove_all(root);a::Store store(root/"main");ExactField field({"x","d"});Exact d(field,"d");auto basis=tadpole(field);
  level::Options options;std::size_t calls=0;
  auto provider=[&](const std::vector<ibp::Integral>& batch,const fire::Options&) {
    ++calls;fire::Result result;result.success=true;
    for(const auto& integral:batch) {
      if(integral==ibp::Integral{1,0})result.reductions.emplace(integral,ibp::Relation{{{1,0},Exact(field,1)}});
      else if(integral==ibp::Integral{2,0})result.reductions.emplace(integral,ibp::Relation{{{1,0},Exact(field,"(1-d/2)/x")}});
      else throw std::runtime_error("unexpected hermetic provider demand");
    }
    return result;
  };
  std::vector<ibp::Integral> requested{{1,0},{2,0}};
  auto first=cached_level::prepare(store,basis,d,field,0,requested,options,provider);
  require(first.result.success && !first.cache_hit && calls==1,"first exact closure builds once");
  require(first.result.matrix[0][0]==Exact(field,"(d/2-1)/x"),"tadpole differential matrix");
  auto scientific=cached_level::identity(basis,d,field,0,requested,options);
  auto record=store.read(scientific,first.content_id);
  require(record.certificate.type=="exact" && record.certificate.scope=="exact_closure","conditional exact closure certificate metadata");
  // New process-equivalent store and exact context; consumer precision and
  // receiving epsilon/Taylor refinements do not change this exact dependency.
  a::Store reopened(root/"main");ExactField context2({"x","d"});Exact d2(context2,"d");auto basis2=tadpole(context2);
  auto tiny=options;tiny.max_passes=0;tiny.provider.executable="/missing/FIRE";tiny.provider.timeout_seconds=0;tiny.provider.memory_bytes=1;
  kernel::ComplexBall::set_precision(512);
  for(const Demand receiving:std::vector<Demand>{{-2,4,48,128,8},{-2,8,80,512,20}}) {
    receiving.validate(); // Deliberately not passed into exact preparation.
    auto again=cached_level::prepare(reopened,basis2,d2,context2,0,requested,tiny);
    require(again.cache_hit && again.content_id==first.content_id && again.result.success,"refinement or missing FIRE prevented cache hit");
    require(again.result.matrix[0][0].variables()==context2.variables(),"cache decoded into wrong exact field");
  }
  require(calls==1,"receiving refinement invoked provider again");
  auto reverse=requested;std::reverse(reverse.begin(),reverse.end());
  auto changed=cached_level::prepare(store,basis,d,field,0,reverse,options,provider);
  require(!changed.cache_hit && calls==2 && changed.result.target_rows[0][0]==Exact(field,"(1-d/2)/x"),"requested ordering ignored");
  auto gram2=tadpole(field,2);changed=cached_level::prepare(store,gram2,d,field,0,requested,options,provider);
  require(!changed.cache_hit && calls==3,"changed exact Gram reused old closure");
  cached_level::Conventions conventions;conventions.normalization["measure"]="different named scalar normalization";
  changed=cached_level::prepare(store,basis,d,field,0,requested,options,provider,conventions);
  require(!changed.cache_hit && calls==4,"normalization identity omitted");
  auto zero_options=options;zero_options.provider.zero_sectors={{0,0}};
  require(cached_level::identity(basis,d,field,0,requested,zero_options).key()!=scientific.key(),"zero-sector assumptions omitted");
  ExactField wider({"x","d","unused"});auto wide_basis=tadpole(wider);Exact wide_d(wider,"d");
  require(cached_level::identity(wide_basis,wide_d,wider,0,requested,options).key()!=scientific.key(),"ordered context symbols omitted");
  auto no_provider=[](const auto&,const auto&){fire::Result result;result.reason="new context requires provider";return result;};
  auto wide=cached_level::prepare(store,wide_basis,wide_d,wider,0,requested,options,no_provider);
  require(!wide.cache_hit && !wide.result.success,"new context was reused");
  rejects([&]{cached_level::identity(basis,d,wider,0,requested,options);},"mismatched exact contexts accepted");
  auto branch=conventions;branch.branch["sheet"]="opposite";
  require(cached_level::identity(basis,d,field,0,requested,options,branch).key()!=
      cached_level::identity(basis,d,field,0,requested,options,conventions).key(),"branch convention omitted");
  require(cached_level::identity(basis,Exact(field,"d+2"),field,0,requested,options).key()!=scientific.key(),"dimension omitted");
  // Mathematical tampering is re-published with VALID payload/envelope hashes.
  // The level verifier must reject it rather than trusting checksum or claim.
  unsigned mutation=0;
  auto reject_math=[&](auto change,const char* why) {
    auto payload=record.payload;change(payload.as_object());a::Store corrupt(root/("math-"+std::to_string(++mutation)));
    corrupt.put(scientific,record.guarantee,payload,record.certificate);
    rejects([&]{cached_level::prepare(corrupt,basis,d,field,0,requested,tiny);},why);
  };
  reject_math([](auto& o){o.at("matrix").as_array()[0].as_array()[0]="0";},"tampered differential equation accepted");
  reject_math([](auto& o){o.at("target_rows").as_array()[0].as_array()[0]="2";},"tampered target coordinate accepted");
  reject_math([](auto& o){o.at("differential_witnesses").as_array()[0].as_array().clear();},"missing differential witness accepted");
  reject_math([](auto& o){o.at("source_identities").as_array()[0].as_array()[0].as_object()["coefficient"]="999";},"tampered source/witness reconstruction accepted");
  reject_math([](auto& o){o.at("matrix").as_array().clear();},"wrong matrix dimensions accepted");
  reject_math([](auto& o){o.at("matrix").as_array()[0].as_array()[0]="unknown";},"unknown coefficient symbol accepted");
  reject_math([](auto& o){std::swap(o.at("requested_integrals").as_array()[0],o.at("requested_integrals").as_array()[1]);},"payload request ordering mismatch accepted");
  reject_math([](auto& o){o.at("differential_witnesses").as_array()[0].as_array()[0].as_object()["source"]=99999;},"out-of-range witness source accepted");
  // Empty spans still need target witnesses; no phantom master is introduced.
  a::Store zeros(root/"zero-span");
  auto zero_provider=[](const auto& batch,const auto&){fire::Result r;r.success=true;for(const auto& a:batch)r.reductions.emplace(a,ibp::Relation{});return r;};
  auto zero=cached_level::prepare(zeros,basis,d,field,0,requested,options,zero_provider);
  auto zero_again=cached_level::prepare(zeros,basis,d,field,0,requested,tiny);
  require(zero.result.success && zero_again.cache_hit && zero_again.result.ordered_basis.empty() &&
      zero_again.result.target_witnesses.size()==requested.size(),"zero-span cache lost exact target witnesses");
  // Ordinary on-disk checksum corruption also fails closed, without FIRE.
  auto disk=store.path(scientific.key(),first.content_id);a::fs::permissions(disk,a::fs::perms::owner_write,a::fs::perm_options::add);
  std::ofstream(disk,std::ios::trunc)<<"{broken";
  rejects([&]{cached_level::prepare(store,basis,d,field,0,requested,tiny);},"corrupted cache silently rebuilt");
  std::cout<<"cached exact level preparation PASS: reuse, semantic inputs, independent witness reconstruction\n";
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
