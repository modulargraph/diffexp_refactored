#include "diffexp/fire_checkpoint.hpp"
#include "diffexp/fire_session.hpp"
#include "diffexp/families.hpp"
#include <iostream>
using namespace diffexp;
void check(bool value,const std::string& why){if(!value)throw std::runtime_error(why);}
int main(int argc,char**argv){try{
 if(argc<2){std::cout<<"backward checkpoint integration requires FIRE executable\n";return 0;}
 ExactField field({"x","eps","I"});Exact x(field,"x"),eps(field,"eps");auto d=x.constant(4)-x.constant(2)*eps;
 auto family=feynman::example_family("box_triangle");auto raw=ibp::quadratic_family(family.momenta,x,family.physical_count);ibp::PropagatorBasis b(ibp::merge(raw,0,1,x));ibp::Integral target(b.denominators.size(),0);std::fill_n(target.begin(),b.physical_count,1);target[0]=3;
 auto dotted=target,numerator=target;dotted[0]=2;dotted[1]=2;numerator[0]=2;numerator.back()=-1;std::vector<ibp::Integral> requests{target,dotted,numerator};
 fire::Options options;options.executable=argv[1];options.timeout_seconds=90;options.max_completed_backward_sectors=3;
 auto interrupted=fire_checkpoint::reduce(b,d,field,requests,options);check(!interrupted.success&&interrupted.reason.find("completed-sector budget")!=std::string::npos,"interrupt after finalized backward sectors");
 auto manifest=boost::json::parse(fire::read_text(interrupted.directory/"checkpoint.json")).as_object();
 check(artifacts::detail::canonical(manifest.at("identity"))==artifacts::detail::canonical(fire_checkpoint::identity(b,d,field,requests,options,nullptr)),"overlay science and exact ordered demand binding");
 check(manifest.contains("backward_overlay"),"explicit finalized-backward checkpoint metadata");
 const auto& overlay=manifest.at("backward_overlay").as_object();auto count=overlay.at("files").as_array().size();check(count>=3,"completed backward worker provenance");
 fire_checkpoint::Store store(interrupted.directory/"pending");auto pending=store.put(interrupted.directory);
 options.max_completed_backward_sectors=0;fire::Session session(b,d,field,options,interrupted.directory/"batches",interrupted.directory/"pending");auto resumed=session(requests,options);check(resumed.success,resumed.reason);
 auto proof=boost::json::parse(fire::read_text(resumed.directory/"checkpoint-resume.json")).as_object();check(proof.at("backward_sector_snapshots").as_int64()==count&&session.checkpoint_hits()==1,"default Session replays all finalized backward snapshots");
 auto log=fire::read_text(resumed.directory/"run.log");check(log.find("SKIPPING LEVEL")==std::string::npos&&log.find("nothing to do")!=std::string::npos,"overlay rebuilds all request propagation without completion flags");
 auto cold=fire::reduce(b,d,field,requests,options);check(cold.success&&cold.reductions==resumed.reductions,"backward overlay final exact reductions differ from cold");
 auto changed=requests;changed[0][0]++;auto rejected=fire_checkpoint::reduce(b,d,field,changed,options,nullptr,interrupted.directory);check(!rejected.success&&rejected.directory.empty(),"backward overlay rejects changed demand before execution");
 {std::ofstream corrupt(interrupted.directory/"checkpoint-run.log",std::ios::app);corrupt<<"changed";}
 rejected=fire_checkpoint::reduce(b,d,field,requests,options,nullptr,interrupted.directory);check(!rejected.success&&rejected.directory.empty(),"backward finalization provenance corruption rejected before execution");
 std::cout<<"finalized backward overlay verified: "<<count<<" sectors; "<<resumed.reductions.size()<<" exact output identities; "<<resumed.directory<<std::endl;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
