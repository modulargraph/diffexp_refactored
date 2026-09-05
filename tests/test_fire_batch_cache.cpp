#include "diffexp/fire_session.hpp"
#include "diffexp/families.hpp"
#include <iostream>
using namespace diffexp;
void check(bool value,const char* why){if(!value)throw std::runtime_error(why);}
int main(int argc,char** argv){try {
  if(argc<2){std::cout<<"durable FIRE batch test requires executable\n";return 0;}
  auto pattern=(std::filesystem::temp_directory_path()/"de3-batch-cache-XXXXXX").string();std::vector<char> name(pattern.begin(),pattern.end());name.push_back(0);check(mkdtemp(name.data()),"temporary cache");std::filesystem::path cache=name.data();
  ExactField field({"x","eps","I"});Exact eps(field,"eps");auto d=eps.constant(4)-eps.constant(2)*eps;
  ibp::PropagatorBasis basis(ibp::quadratic_family(feynman::banana(1,{Rational(1),Rational(2)}),eps));fire::Options options;options.executable=argv[1];
  fire::Result first;
  {fire::Session producer(basis,d,field,options,cache);first=producer({{2,1}},options);check(first.success,first.reason.c_str());check(producer.runs()==1,"one cold completed batch");}
  {fire::Options offline=options;offline.executable.clear();offline.threads=4;fire::Session reader(basis,d,field,offline,cache);auto hit=reader({{2,1}},offline);check(hit.success,hit.reason.c_str());check(hit.reductions==first.reductions&&reader.runs()==0&&reader.durable_hits()==1,"durable re-import with no executable");}
  fire::Result grown;
  {fire::Session continuation(basis,d,field,options,cache);auto hit=continuation({{2,1}},options);check(hit.success&&continuation.runs()==0,"reuse predecessor before growth");grown=continuation({{2,1},{4,1}},options);check(grown.success,grown.reason.c_str());check(continuation.runs()==1,"only new demand executes");}
  fire_batch::Cache recovered(cache/"recovered");recovered.recover_replayed(basis,d,field,options,first.directory,grown.directory);
  check(recovered.lookup(basis,d,field,{{2,1}},options)->reductions==first.reductions,"old accepted replay recovers exact successful batch");
  auto cold=fire::reduce(basis,d,field,{{2,1},{4,1}},options);check(cold.success&&cold.reductions==grown.reductions,"durable continuation exact match with cold");
  fire_batch::Cache storage(cache);auto changed=options;changed.zero_sectors={{1,-1}};check(!storage.lookup(basis,d,field,{{2,1}},changed),"zero sector change invalidates durable batch");
  check(!storage.lookup(basis,d,field,{{4,1},{2,1}},options),"ordered demand identity is exact");
  auto failed=first;failed.success=false;bool rejected=false;try{storage.put(basis,d,field,{{2,1}},options,failed);}catch(const std::exception&){rejected=true;}check(rejected,"failed provider result cannot be cached");
  auto id=fire_batch::identity(basis,d,field,{{2,1}},options);artifacts::Store store(cache);auto record=store.lookup(id,{0,0,0,64,0});check(bool(record),"batch artifact published");
  auto damaged=record->payload.as_object();damaged["mappings"]=boost::json::array{};
  artifacts::Store malformed_store(cache/"malformed");malformed_store.put(id,{0,0,0,64,0},damaged,record->certificate);
  fire_batch::Cache malformed(cache/"malformed");rejected=false;
  try{malformed.lookup(basis,d,field,{{2,1}},options);}catch(const std::exception&){rejected=true;}
  check(rejected,"self-consistent artifact hash cannot hide incorrect exact mappings");
  // Independently reject corrupted bytes before mathematical payload decoding.
  std::filesystem::permissions(store.path(record->semantic_id,record->content_id),std::filesystem::perms::owner_write,std::filesystem::perm_options::add);
  {std::ofstream out(store.path(record->semantic_id,record->content_id),std::ios::app);out<<"broken";check(bool(out),"corruption fixture writes bytes");}
  rejected=false;try{storage.lookup(basis,d,field,{{2,1}},options);}catch(const std::exception&){rejected=true;}check(rejected,"corrupt durable batch fails closed");
  std::filesystem::remove_all(cache);std::cout<<"durable FIRE batch cache passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
