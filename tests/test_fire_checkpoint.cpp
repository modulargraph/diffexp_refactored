#include "diffexp/fire_checkpoint.hpp"
#include "diffexp/fire_session.hpp"
#include "diffexp/families.hpp"
#include <iostream>
using namespace diffexp;
void check(bool condition,const std::string& message){if(!condition)throw std::runtime_error(message);}
int main(int argc,char**argv){try {
 if(argc<2){std::cout<<"checkpoint integration requires FIRE executable\n";return 0;}
 ExactField field({"x","eps","I"});Exact x(field,"x"),eps(field,"eps");auto d=x.constant(4)-x.constant(2)*eps;
 auto f=feynman::example_family("kite");auto raw=ibp::quadratic_family(f.momenta,x,f.physical_count);ibp::PropagatorBasis b(ibp::merge(raw,0,1,x));
 ibp::Integral target(b.denominators.size(),0);std::fill_n(target.begin(),b.physical_count,1);target[0]=3;fire::Options options;options.executable=argv[1];options.timeout_seconds=60;options.max_completed_forward_sectors=1;
 auto first=fire_checkpoint::reduce(b,d,field,{target},options);check(!first.success&&first.reason.find("completed-sector budget")!=std::string::npos,"fixture must interrupt a real reduction");
 check(!std::filesystem::exists(first.directory/"success.receipt"),"interrupted table cannot be marked successful");
 fire_checkpoint::Store store(first.directory/"durable-pending");auto published=store.put(first.directory);
 auto bound=fire_checkpoint::identity(b,d,field,{target},options,nullptr);check(store.lookup(bound)==published,"durable pending fixed identity lookup");
 options.timeout_seconds=60;options.max_completed_forward_sectors=0;auto changed=target;++changed[0];
 auto wrong=fire_checkpoint::reduce(b,d,field,{changed},options,nullptr,first.directory);check(!wrong.success&&wrong.directory.empty(),"changed demand rejected before execution");
 auto workers=options;workers.threads=4;wrong=fire_checkpoint::reduce(b,d,field,{target},workers,nullptr,first.directory);check(!wrong.success&&wrong.directory.empty(),"changed checkpoint config rejected before execution");
 // Incomplete copies and level flags are deliberately excluded from replay.
 {std::ofstream out(published/"storage"/"00099.tmp.copying");out<<"partial";}
 {std::ofstream out(published/"storage"/"completed_in_storage");out<<"5 1\n4 1\n3 1\n2 1\n1 1\n";}
 fire::Session session(b,d,field,options,first.directory/"completed-batches",first.directory/"durable-pending");
 auto resumed=session({target},options);check(resumed.success,resumed.reason);
 check(session.checkpoint_hits()==1&&!store.lookup(bound),"Session resumes pending algebra and retires it after durable successful batch");
 auto log=fire::read_text(resumed.directory/"run.log");check(log.find("nothing to do")!=std::string::npos&&log.find("SKIPPING LEVEL")==std::string::npos,"reuse solved sectors while rebuilding all request propagation");
 auto cold=fire::reduce(b,d,field,{target},options);check(cold.success&&cold.reductions==resumed.reductions,"resumed exact output must equal cold output");
 options.timeout_seconds=60;options.max_completed_forward_sectors=1;auto raw_interrupted=fire::reduce(b,d,field,{target},options);check(!raw_interrupted.success,"raw recovery fixture interrupts");
 auto recovered=fire_checkpoint::recover_forward_run(b,d,field,{target},options,nullptr,raw_interrupted.directory);
 options.timeout_seconds=60;options.max_completed_forward_sectors=0;auto recovered_result=fire_checkpoint::reduce(b,d,field,{target},options,nullptr,recovered);
 check(recovered_result.success&&recovered_result.reductions==cold.reductions,"finalized raw forward snapshots resume to exact cold result");
 auto entries=fire_checkpoint::files(first.directory/"storage");check(!entries.empty(),"interrupted fixture retains a finalized snapshot");
 auto snapshot=first.directory/"storage"/artifacts::detail::string(entries.front().as_object().at("name"));{std::ofstream out(snapshot,std::ios::app);out<<"bad";}
 wrong=fire_checkpoint::reduce(b,d,field,{target},options,nullptr,first.directory);check(!wrong.success&&wrong.directory.empty(),"corrupt snapshot rejected before execution");
 std::cout<<"fixed-demand interrupted FIRE checkpoint verified against cold reduction\n";
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
