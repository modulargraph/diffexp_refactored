#include "diffexp/artifact_store.hpp"
#include <iostream>
#include <thread>
using namespace diffexp;
namespace a=diffexp::artifacts;
void require(bool yes,const char* why){if(!yes)throw std::runtime_error(why);}
template<class F> void rejects(F f,const char* why){bool failed=false;try{f();}catch(const std::exception&){failed=true;}require(failed,why);}
a::Identity identity() {
  a::Identity i;i.kind="exact_equation";i.algorithm_version="native-v1";
  i.family={{"name","henn"},{"dimension","4-2*eps"}};
  i.ordered_basis={"I[1,0,1]","I[1,1,0]"};i.normalization={{"scalar","-q^2"},{"prefactor","eps^4*Exp[2*eps*EulerGamma]"}};
  i.branch={{"rim","lower"},{"eps5","+I*Sqrt[3]"}};i.geometry={{"x0","{3,-1,1,1,-1}"}};
  i.boundary={{"status","not-applicable"}};return i;
}
a::json::object read_raw(const a::fs::path& p){std::ifstream in(p);std::string text((std::istreambuf_iterator<char>(in)),{});return a::json::parse(text).as_object();}
void write_raw(const a::fs::path& p,const a::json::object& o) {
  a::fs::permissions(p,a::fs::perms::owner_write,a::fs::perm_options::add);
  std::ofstream out(p,std::ios::binary|std::ios::trunc);out<<a::detail::canonical(o);
}
int main(){try {
  require(a::detail::sha256("")=="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855","SHA256 empty vector");
  require(a::detail::sha256("abc")=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad","SHA256 abc vector");
  require(a::detail::sha256(std::string(1000000,'a'))=="cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0","SHA256 multiblock million-a vector");
  auto root=a::fs::temp_directory_path()/("de3-artifact-store-"+std::to_string(::getpid()));
  struct Cleanup{a::fs::path p;~Cleanup(){std::error_code ec;a::fs::remove_all(p,ec);}}cleanup{root};
  a::fs::remove_all(root);a::Store store(root);auto upstream=identity();
  auto reordered=upstream;reordered.family={{"dimension","4-2*eps"},{"name","henn"}};
  require(upstream.key()==reordered.key(),"object insertion order must not change scientific identity");
  auto wrong_basis=upstream;std::swap(wrong_basis.ordered_basis[0],wrong_basis.ordered_basis[1]);
  require(upstream.key()!=wrong_basis.key(),"ordered basis must affect identity");
  auto wrong_branch=upstream;wrong_branch.branch["rim"]="upper";
  require(upstream.key()!=wrong_branch.key(),"branch must affect identity");
  auto incomplete=upstream;incomplete.boundary.clear();rejects([&]{incomplete.key();},"missing identity axis accepted");
  Demand fixed{0,0,0,64,0};a::Certificate exact{"exact","native-rational-v1","full-equation",{{"method","exact canonical rational equality"}}};
  auto original=store.put(upstream,fixed,a::json::object{{"equation","(eps+1)/(x-1)"}},exact);
  auto same=store.put(reordered,fixed,a::json::object{{"equation","(eps+1)/(x-1)"}},exact);
  require(original.content_id==same.content_id,"idempotent immutable publication");
  a::Store reopened(root);auto found=reopened.lookup(upstream,fixed,{"exact","native-rational-v1","full-equation"});
  require(found && found->content_id==original.content_id,"durable exact roundtrip");
  require(!reopened.lookup(wrong_basis,fixed) && !reopened.lookup(wrong_branch,fixed),"scientific mismatch reused artifact");
  rejects([&]{reopened.read(wrong_basis,original.content_id);},"wrong identity direct load accepted");
  // Receiving-only extension across fresh store instances. Upstream semantic
  // identity/demand stay fixed; it is published once and reused at both widths.
  auto receiving=upstream;receiving.kind="receiving_basis";receiving.parents={{"equation",upstream.key()}};
  receiving.boundary={{"incoming","fixed-exact-boundary-v1"}};
  unsigned upstream_builds=1,receiving_builds=0;
  auto obtain=[&](Demand demand) {
    a::Store disk(root);if(auto cached=disk.lookup(receiving,demand))return *cached;
    auto parent=disk.lookup(upstream,fixed,{"exact","native-rational-v1","full-equation"});
    if(!parent){++upstream_builds;throw std::runtime_error("upstream unexpectedly rebuilt");}
    ++receiving_builds;
    return disk.put(receiving,demand,a::json::object{{"epsilon_high",demand.epsilon_high},{"taylor_order",demand.taylor_order}}, {},{parent->content_id});
  };
  Demand small{-2,4,48,128,8},large{-2,8,80,256,16};
  auto r48=obtain(small),r80=obtain(large),again=obtain(small);
  require(upstream_builds==1 && receiving_builds==2 && again.content_id==r48.content_id,"local refinement repeated upstream or chose excessive record");
  require(r48.semantic_id==r80.semantic_id && r48.content_id!=r80.content_id,"resource demand leaked into scientific identity");
  require(a::fs::exists(store.path(r48.semantic_id,r48.content_id)) && a::fs::exists(store.path(r80.semantic_id,r80.content_id)),"refinement overwrote an immutable artifact");
  require(!store.lookup(receiving,Demand{-3,8,80,256,16}),"epsilon-low dominance ignored");
  require(!store.lookup(receiving,Demand{-2,9,80,256,16}),"epsilon-high dominance ignored");
  require(!store.lookup(receiving,Demand{-2,8,81,256,16}),"Taylor dominance ignored");
  require(!store.lookup(receiving,Demand{-2,8,80,257,16}),"bit dominance ignored");
  require(!store.lookup(receiving,Demand{-2,8,80,256,17}),"publication resource dominance ignored");
  require(!store.lookup(receiving,small,{"certified_enclosure",{}, {}}),"requested digits misrepresented as numerical certification");
  rejects([&]{store.put(receiving,small,a::json::object{}, {},{});},"missing concrete parent accepted");
  rejects([&]{store.put(receiving,small,a::json::object{}, {},{r48.content_id});},"parent content/semantic mismatch accepted");
  rejects([&]{store.put(upstream,fixed,a::json::object{{"float",0.1}},exact);},"inexact JSON float accepted");
  rejects([&]{store.put(upstream,fixed,a::json::object{},a::Certificate{"certified_enclosure",{},{},{}});},"unbound certification claim accepted");
  // Crash leftovers are invisible. Parallel identical writers publish one
  // immutable content record, never a partially written destination.
  std::ofstream(root/upstream.key()/".publish-crashed")<<"partial";
  std::array<std::string,4> ids;std::array<std::exception_ptr,4> errors;std::vector<std::thread> threads;
  for(unsigned i=0;i<4;++i)threads.emplace_back([&,i]{try{a::Store local(root);ids[i]=local.put(upstream,fixed,a::json::object{{"equation","new exact artifact"}},exact).content_id;}catch(...){errors[i]=std::current_exception();}});
  for(auto& t:threads)t.join();for(unsigned i=0;i<4;++i){if(errors[i])std::rethrow_exception(errors[i]);require(ids[i]==ids[0],"parallel atomic publication differs");}
  // Tampering must fail closed, including when the file retains valid JSON.
  auto target=store.path(original.semantic_id,original.content_id);auto raw=read_raw(target);raw["payload"].as_object()["equation"]="CORRUPTED";write_raw(target,raw);
  rejects([&]{store.read(upstream,original.content_id);},"payload corruption accepted");
  rejects([&]{store.lookup(upstream,fixed);},"lookup silently ignored corruption");
  // Recompute only the envelope identity to exercise the independent payload
  // checksum and certificate-to-payload binding checks.
  raw.erase("content_id");auto forged=a::detail::sha256(a::detail::canonical(raw));raw["content_id"]=forged;
  auto forged_path=store.path(original.semantic_id,forged);std::ofstream(forged_path)<<a::detail::canonical(raw);
  rejects([&]{store.read(upstream,forged);},"independent payload checksum not checked");
  raw["payload_sha256"]=a::detail::sha256(a::detail::canonical(raw.at("payload")));raw.erase("content_id");
  forged=a::detail::sha256(a::detail::canonical(raw));raw["content_id"]=forged;
  std::ofstream(store.path(original.semantic_id,forged))<<a::detail::canonical(raw);
  rejects([&]{store.read(upstream,forged);},"certificate payload binding not checked");
  // A structurally invalid guarantee is rejected even when the outer digest
  // is recomputed, and schema migrations never silently reuse old records.
  auto valid=read_raw(store.path(original.semantic_id,ids[0]));
  auto publish_forged=[&](a::json::object object) {
    object.erase("content_id");auto id=a::detail::sha256(a::detail::canonical(object));object["content_id"]=id;
    std::ofstream(store.path(original.semantic_id,id))<<a::detail::canonical(object);return id;
  };
  auto bad_guarantee=valid;bad_guarantee["guarantee"].as_object()["epsilon_low"]=5;
  auto invalid_id=publish_forged(bad_guarantee);
  rejects([&]{store.read(upstream,invalid_id);},"invalid guarantee interval accepted");
  auto bad_schema=valid;bad_schema["schema"]="DiffExp3.ArtifactRecord/v999";invalid_id=publish_forged(bad_schema);
  rejects([&]{store.read(upstream,invalid_id);},"unknown record schema accepted");
  auto extra=valid;extra["unexpected_field"]=0;invalid_id=publish_forged(extra);
  rejects([&]{store.read(upstream,invalid_id);},"unknown schema field accepted");
  auto good_path=store.path(original.semantic_id,ids[0]);
  a::fs::permissions(good_path,a::fs::perms::owner_write,a::fs::perm_options::add);
  std::ofstream(good_path,std::ios::trunc)<<"{\"truncated\":";
  rejects([&]{store.read(upstream,ids[0]);},"truncated artifact accepted");
  std::cout<<"artifact store PASS: durable local refinement, identity, dominance, atomic writes, corruption rejection\n";
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
