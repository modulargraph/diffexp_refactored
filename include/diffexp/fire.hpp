#pragma once
#include "diffexp/ibp.hpp"
#include <filesystem>
#include <chrono>
#include <thread>
#include <fstream>
#include <sys/wait.h>
#include <sys/resource.h>
#include <fcntl.h>
#include <unistd.h>
#include <signal.h>
#include <cerrno>
#ifdef __APPLE__
#include <libproc.h>
#endif

namespace diffexp::fire {
// Native FIRE7 interchange. Source contract: sources/parser.cpp add_problem and
// split_coeff. No Wolfram evaluator or generated executable source is involved.
// FIRE restricts symbols to lowercase names. Rename complete identifier tokens
// simultaneously, so prefixes and names that coincide with temporary symbols
// cannot corrupt coefficients. The caller's exact field remains unchanged.
struct SymbolMap {
  std::vector<std::string> exported;
  std::map<std::string,std::string> forward,reverse;
  explicit SymbolMap(const std::vector<std::string>& names) {
    for(std::size_t i=0;i<names.size();++i) {
      auto number=i;std::string suffix;
      do{suffix.push_back(char('a'+number%26));number/=26;}while(number);
      auto temporary="devar"+suffix;exported.push_back(temporary);
      forward.emplace(names[i],temporary);reverse.emplace(temporary,names[i]);
    }
  }
};
inline std::string rename_symbols(const std::string& source,const std::map<std::string,std::string>& names) {
  std::string out;
  for(std::size_t i=0;i<source.size();) {
    if(std::isalpha(static_cast<unsigned char>(source[i]))||source[i]=='_') {
      auto begin=i++;while(i<source.size()&&(std::isalnum(static_cast<unsigned char>(source[i]))||source[i]=='_'))++i;
      auto symbol=source.substr(begin,i-begin);auto found=names.find(symbol);out+=found==names.end()?symbol:found->second;
    }else out.push_back(source[i++]);
  }
  return out;
}
inline std::string indices(const ibp::Integral& a) {
  std::string s="{";for(std::size_t i=0;i<a.size();++i){if(i)s+=", ";s+=std::to_string(a[i]);}return s+"}";
}
// A common null direction of every present quadratic leaves a free loop
// integration and therefore proves the sector scaleless. This criterion is
// conservative; it does not guess massless scaling or discrete symmetries.
inline std::vector<ibp::Integral> free_loop_sectors(const ibp::PropagatorBasis& b) {
  if(b.physical_count>12)throw std::invalid_argument("FIRE zero-sector enumeration exceeds 4096-sector budget");
  std::vector<ibp::Integral> result;auto zero=b.space.zero().constant;
  for(std::size_t mask=0;mask<(std::size_t(1)<<b.physical_count);++mask) {
    std::map<std::size_t,std::vector<Exact>> echelon;
    for(std::size_t k=0;k<b.physical_count;++k)if(mask&(std::size_t(1)<<k)) {
      std::vector<std::vector<Exact>> matrix(b.space.loops+b.space.externals(),std::vector<Exact>(b.space.loops,zero));
      for(std::size_t q=0;q<b.space.pairs.size();++q) {
        auto [a,c]=b.space.pairs[q];auto value=b.denominators[k].linear[q];
        matrix[c][a]=matrix[c][a]+value;
        if(c<b.space.loops)matrix[a][c]=matrix[a][c]+value;
      }
      for(auto row:matrix) {
        for(const auto& [pivot,r]:echelon){auto factor=row[pivot];for(std::size_t j=pivot;j<row.size();++j)row[j]=row[j]-factor*r[j];}
        auto p=std::find_if(row.begin(),row.end(),[](const Exact& c){return !c.is_zero();});
        if(p!=row.end()){auto pivot=p-row.begin();auto factor=*p;for(auto& c:row)c=c/factor;echelon.emplace(pivot,std::move(row));}
      }
    }
    if(echelon.size()<b.space.loops){ibp::Integral sector(b.denominators.size(),-1);for(std::size_t k=0;k<b.physical_count;++k)if(mask&(std::size_t(1)<<k))sector[k]=1;result.push_back(std::move(sector));}
  }
  return result;
}
inline std::string start(const ibp::PropagatorBasis& b,const Exact& dimension,
    const std::vector<ibp::Integral>& zero_sectors={},const SymbolMap* symbols=nullptr) {
  auto n=b.denominators.size();
  if(n>20)throw std::invalid_argument("native FIRE preparation limited to 20 slots");
  std::ostringstream out;out<<"ExampleDimension[1] = "<<n<<"\n\nSBasis0L[1] = "
    <<b.space.loops*(b.space.loops+b.space.externals())<<"\n\n";
  std::size_t r=0;
  for(unsigned i=0;i<b.space.loops;++i)for(unsigned v=0;v<b.space.loops+b.space.externals();++v) {
    ++r;std::map<ibp::Integral,std::map<unsigned,Exact>> coefficients;
    auto add=[&](const ibp::Integral& shift,unsigned index,const Exact& c){ibp::add(coefficients[shift],index,c);};
    ibp::Integral shift(n,0);if(i==v)add(shift,0,dimension);
    for(std::size_t k=0;k<n;++k) {
      auto contraction=b.rewrite(b.space.contraction(b.denominators[k],i,v));
      shift.assign(n,0);++shift[k];add(shift,k+1,-contraction.constant);
      for(std::size_t j=0;j<n;++j){--shift[j];add(shift,k+1,-contraction.linear[j]);++shift[j];}
    }
    for(const auto& [s,terms]:coefficients)if(!terms.empty()) {
      out<<"SBasis0C[1, "<<r<<", "<<indices(s)<<"] = {";bool comma=false;
      for(const auto& [index,c]:terms){if(comma)out<<", ";comma=true;out<<"{"<<(symbols?rename_symbols(c.str(),symbols->forward):c.str())<<", "<<index<<"}";}out<<"}\n\n";
    }
  }
  ibp::Integral identity(n),ones(n,1),zeros(n,0);std::iota(identity.begin(),identity.end(),1);
  out<<"SBasisS[1] = {{"<<indices(identity)<<", "<<indices(ones)<<", "<<indices(zeros)<<"}}\n\n";
  auto all_zero=free_loop_sectors(b);all_zero.insert(all_zero.end(),zero_sectors.begin(),zero_sectors.end());
  std::sort(all_zero.begin(),all_zero.end());all_zero.erase(std::unique(all_zero.begin(),all_zero.end()),all_zero.end());
  for(const auto& s:all_zero) {
    if(s.size()!=n || std::any_of(s.begin(),s.end(),[](int a){return a!=-1&&a!=1;}))
      throw std::invalid_argument("FIRE zero sector must contain one sign per slot");
    for(std::size_t k=b.physical_count;k<n;++k)if(s[k]>0)throw std::invalid_argument("positive numerator sector");
    out<<"SBasisR[1, "<<indices(s)<<"] = True\n\n";
  }
  return out.str();
}
struct Options {
  std::filesystem::path executable;
  std::filesystem::path scratch_root=std::filesystem::temp_directory_path();
  unsigned timeout_seconds=60;
  unsigned threads=1,simplifier_threads=1;
  unsigned max_completed_forward_sectors=0,max_completed_backward_sectors=0; // Optional deterministic progress slices; zero disables.
  std::size_t memory_bytes=std::size_t(4)*1024*1024*1024;
  std::vector<ibp::Integral> zero_sectors; // Caller supplies proven scaleless sectors only.
};
struct Result {
  bool success=false;
  std::string reason;
  std::filesystem::path directory;
  std::map<ibp::Integral,ibp::Relation> reductions;
};
// A probe is data for reconstruction, never a completed exact reduction.
struct ModularPoint {
  std::vector<std::uint64_t> values;
  unsigned prime_index=1;
  std::string suffix() const {
    std::string out;for(auto value:values){if(!out.empty())out+="_";out+=std::to_string(value);}
    return out+"_"+std::to_string(prime_index);
  }
};
inline std::string expression(const data::Expr& e) {
  if(e.args.empty())return e.head;
  if(e.head=="neg"&&e.args.size()==1)return "-("+expression(e.args[0])+")";
  if(e.args.size()==2&&(e.head=="+"||e.head=="-"||e.head=="*"||e.head=="/"||e.head=="^"))
    return "("+expression(e.args[0])+")"+e.head+"("+expression(e.args[1])+")";
  throw std::invalid_argument("unexpected FIRE coefficient expression");
}
inline std::int64_t table_integer(const data::Expr& e) {
  auto s=expression(e);if(e.head=="neg"&&e.args.size()==1)s="-"+e.args[0].head;
  std::size_t p=0;auto n=std::stoll(s,&p);if(p!=s.size())throw std::invalid_argument("FIRE integer expected");return n;
}
inline std::string table_id(const data::Expr& e) {
  if(!e.args.empty()||e.head.empty()||!std::all_of(e.head.begin(),e.head.end(),[](unsigned char c){return std::isdigit(c);}))throw std::invalid_argument("FIRE table id expected");return e.head;
}
inline const std::vector<data::Expr>& list(const data::Expr& e,std::size_t size=std::size_t(-1)) {
  if(e.head!="List"||(size!=std::size_t(-1)&&e.args.size()!=size))throw std::invalid_argument("FIRE table list shape");return e.args;
}
inline std::map<ibp::Integral,ibp::Relation> import_table_text(std::string text,const ExactField& field,std::size_t slots,const SymbolMap* symbols=nullptr) {
  if(text.size()>64*1024*1024)throw std::runtime_error("FIRE tables exceed 64 MiB import budget");
  // FIRE encloses rational-function coefficients in double quotes. They are
  // still parsed as data, and never interpreted by a shell or evaluator.
  text.erase(std::remove(text.begin(),text.end(),'"'),text.end());
  auto parsed=data::Reader(text).read();const auto& table=list(parsed,2);
  std::map<std::string,ibp::Integral> ids;
  for(const auto& item:list(table[1])) {
    const auto& pair=list(item,2);const auto& integral=list(pair[1],2);
    if(table_integer(integral[0])!=1)throw std::invalid_argument("unexpected FIRE problem number");
    ibp::Integral a;for(const auto& x:list(integral[1],slots)){auto value=table_integer(x);if(value< -100||value>100)throw std::invalid_argument("FIRE output index range");a.push_back(static_cast<int>(value));}
    if(!ids.emplace(table_id(pair[0]),a).second)throw std::invalid_argument("duplicate FIRE integral id");
  }
  std::map<ibp::Integral,ibp::Relation> result;
  for(const auto& item:list(table[0])) {
    const auto& pair=list(item,2);ibp::Relation row;
    for(const auto& term:list(pair[1])) {
      const auto& t=list(term,2);ibp::add(row,ids.at(table_id(t[0])),Exact(field,symbols?rename_symbols(expression(t[1]),symbols->reverse):expression(t[1])));
    }
    auto key=ids.at(table_id(pair[0]));auto it=result.find(key);
    if(it!=result.end() && it->second!=row)throw std::invalid_argument("inconsistent duplicate FIRE reduction");
    result.insert_or_assign(std::move(key),std::move(row));
  }
  return result;
}
inline std::string read_text(const std::filesystem::path& path,std::size_t limit=64*1024*1024) {
  std::ifstream in(path);if(!in)throw std::runtime_error("cannot read FIRE data: "+path.string());
  if(std::filesystem::file_size(path)>limit)throw std::runtime_error("FIRE data exceeds size budget");
  return std::string((std::istreambuf_iterator<char>(in)),{});
}
inline std::map<ibp::Integral,ibp::Relation> import_tables(const std::filesystem::path& path,const ExactField& field,std::size_t slots,const SymbolMap* symbols=nullptr) {
  return import_table_text(read_text(path),field,slots,symbols);
}
inline std::string rules_text(const std::map<ibp::Integral,ibp::Relation>& prior,const Exact& sample,const SymbolMap& symbols) {
  std::ostringstream rules;
  for(const auto& [a,row]:prior) {
    if(row.size()==1&&row.begin()->first==a&&row.begin()->second==sample.constant(1))continue;
    rules<<"G[1, "<<indices(a)<<"] -> {";bool comma=false;
    for(const auto& [master,c]:row){if(comma)rules<<", ";comma=true;rules<<"{"<<rename_symbols(c.str(),symbols.forward)<<", G[1, "<<indices(master)<<"]}";}
    rules<<"};\n\n";
  }
  return rules.str();
}
inline Result reduce(const ibp::PropagatorBasis& basis,const Exact& dimension,
    const ExactField& field,const std::vector<ibp::Integral>& requests,const Options& options,
    const std::map<ibp::Integral,ibp::Relation>* prior_rules=nullptr,
    const std::filesystem::path* checkpoint_storage=nullptr,bool preserve_checkpoint=false,
    const ModularPoint* modular_point=nullptr) {
  Result result;pid_t owned_pid=-1;
  try {
    (void)(dimension+basis.space.zero().constant);
    (void)(dimension+Exact(field)); // Reject mixed coefficient fields before launching.
    if(requests.empty())throw std::invalid_argument("FIRE requires at least one requested integral");
    if(modular_point&&(modular_point->values.size()!=field.variables().size()||!modular_point->prime_index||modular_point->prime_index>127||prior_rules||checkpoint_storage||preserve_checkpoint))
      throw std::invalid_argument("invalid modular FIRE point or symbolic checkpoint reuse");
    if(!options.timeout_seconds||!options.memory_bytes)throw std::invalid_argument("FIRE requires finite positive resource budgets");
    if(!options.threads||options.threads>64||!options.simplifier_threads||options.simplifier_threads>64||options.threads*options.simplifier_threads>256)throw std::invalid_argument("FIRE worker counts must be 1..64 with total concurrency at most 256");
    if(options.executable.empty()||access(options.executable.c_str(),X_OK))throw std::runtime_error("FIRE executable unavailable");
    auto executable=std::filesystem::absolute(options.executable);
    for(const auto& a:requests) {
      if(a.size()!=basis.denominators.size())throw std::invalid_argument("FIRE request shape");
      for(std::size_t k=0;k<a.size();++k)if(a[k]<-100||a[k]>100||(k>=basis.physical_count&&a[k]>0))throw std::invalid_argument("FIRE index outside supported range");
    }
    std::filesystem::create_directories(options.scratch_root);
    auto pattern=(std::filesystem::absolute(options.scratch_root)/"diffexp3-fire-XXXXXX").string();
    std::vector<char> name(pattern.begin(),pattern.end());name.push_back(0);
    if(!mkdtemp(name.data()))throw std::runtime_error("cannot create unique FIRE directory");result.directory=name.data();
    auto write=[&](const char* file,const std::string& data){std::ofstream out(result.directory/file);out<<data;if(!out)throw std::runtime_error("cannot write FIRE input");};
    if(preserve_checkpoint||checkpoint_storage) {
      const auto storage=result.directory/"storage";std::filesystem::create_directory(storage);
      if(checkpoint_storage)for(const auto& entry:std::filesystem::directory_iterator(*checkpoint_storage)) {
        const auto name=entry.path().filename().string();
        // Level flags are unsafe across interruption before request propagation;
        // unfinished .copying files are not committed sector snapshots either.
        if(name.size()!=9||name.substr(5)!=".tmp"||!std::all_of(name.begin(),name.begin()+5,[](unsigned char c){return std::isdigit(c);}))continue;
        if(entry.is_symlink()||!entry.is_regular_file())throw std::runtime_error("invalid FIRE checkpoint sector file");
        std::filesystem::copy_file(entry.path(),storage/name);
      }
    }
    SymbolMap symbols(field.variables());
    write("family.start",start(basis,dimension,options.zero_sectors,&symbols));
    std::ostringstream config;config<<"#threads "<<options.threads<<"\n#fthreads "<<options.simplifier_threads<<"\n#compressor none\n#variables ";
    for(std::size_t i=0;i<symbols.exported.size();++i){if(i)config<<",";config<<symbols.exported[i];}
    if(preserve_checkpoint||checkpoint_storage)config<<"\n#storage "<<(result.directory/"storage").string();
    config<<"\n#start\n#folder "<<result.directory.string()<<"/\n#problem 1 |"<<basis.physical_count<<"|family.start\n";
    if(prior_rules&&!prior_rules->empty()) {
      auto rules=rules_text(*prior_rules,dimension,symbols);
      if(!rules.empty()){write("prior.rules",rules);config<<"#rules prior.rules\n";}
    }
    config<<"#integrals integrals.m\n#output "<<(result.directory/"result.tables").string()<<"\n";
    write("job.config",config.str());std::string integrals="{";for(std::size_t i=0;i<requests.size();++i){if(i)integrals+=", ";integrals+="{1, "+indices(requests[i])+"}";}write("integrals.m",integrals+"}\n");
    auto job=(result.directory/"job").string(),db=(result.directory/"database").string(),log=(result.directory/"run.log").string();
    const auto probe_suffix=modular_point?modular_point->suffix():std::string();
    auto pid=fork();if(pid<0)throw std::runtime_error("cannot fork FIRE process");
    if(pid==0) {
      setpgid(0,0);
      // CPU time sums concurrent threads, whereas the parent enforces elapsed
      // wall time and aggregate resident memory over the whole process group.
      const rlim_t cpu_seconds=rlim_t(options.timeout_seconds)*options.threads*options.simplifier_threads;
      struct rlimit cpu{cpu_seconds+1,cpu_seconds+2};
      if(setrlimit(RLIMIT_CPU,&cpu))_exit(125);
#ifndef __APPLE__
      struct rlimit memory{options.memory_bytes,options.memory_bytes};
      if(setrlimit(RLIMIT_AS,&memory))_exit(125);
#endif
      int fd=open(log.c_str(),O_CREAT|O_WRONLY|O_TRUNC,0600);if(fd<0)_exit(126);dup2(fd,1);dup2(fd,2);close(fd);
      if(modular_point)execl(executable.c_str(),executable.c_str(),"-c",job.c_str(),"-d",db.c_str(),"--variables",probe_suffix.c_str(),static_cast<char*>(nullptr));
      else execl(executable.c_str(),executable.c_str(),"-c",job.c_str(),"-d",db.c_str(),static_cast<char*>(nullptr));_exit(127);
    }
    owned_pid=pid;setpgid(pid,pid);int status=0;auto deadline=std::chrono::steady_clock::now()+std::chrono::seconds(options.timeout_seconds);
    while(true) {
      auto waited=waitpid(pid,&status,WNOHANG);if(waited==pid)break;
      if(waited<0 && errno!=EINTR){throw std::runtime_error("cannot wait for FIRE process");}
#ifdef __APPLE__
      // Darwin does not implement RLIMIT_AS. Bound the resident memory of the
      // entire FIRE process group instead, including any FLAME workers.
      std::vector<pid_t> members(4096);
      int bytes=proc_listpids(PROC_PGRP_ONLY,pid,members.data(),members.size()*sizeof(pid_t));
      std::size_t resident=0;
      if(bytes<=0){if(waitpid(pid,&status,WNOHANG)==pid)break;throw std::runtime_error("cannot inspect FIRE memory usage");}
      for(std::size_t k=0;k<std::size_t(bytes)/sizeof(pid_t);++k){struct proc_taskinfo info{};
        if(proc_pidinfo(members[k],PROC_PIDTASKINFO,0,&info,sizeof(info))==sizeof(info))resident+=info.pti_resident_size;}
      if(resident>options.memory_bytes){throw std::runtime_error("FIRE memory budget exceeded; process group terminated");}
#endif
      if((options.max_completed_forward_sectors||options.max_completed_backward_sectors)&&std::filesystem::exists(log)) {
        auto progress=read_text(log);std::size_t position=0;unsigned completed=0,backward_completed=0;
        while((position=progress.find("FLAME time (",position))!=std::string::npos){position+=12;if(position<progress.size()&&progress.find("):",position)!=std::string::npos){if(std::isdigit(static_cast<unsigned char>(progress[position])))++completed;else if(progress[position]=='-')++backward_completed;}}
        if((options.max_completed_forward_sectors&&completed>=options.max_completed_forward_sectors)||(options.max_completed_backward_sectors&&backward_completed>=options.max_completed_backward_sectors)){throw std::runtime_error("FIRE completed-sector budget reached; process group terminated");}
      }
      if(std::chrono::steady_clock::now()>=deadline){throw std::runtime_error("FIRE timeout; process group terminated");}
      std::this_thread::sleep_for(std::chrono::milliseconds(25));
    }
    kill(-pid,SIGKILL);owned_pid=-1; // Terminate any orphaned simplifier/worker processes in our private group.
    if(!WIFEXITED(status)||WEXITSTATUS(status)!=0)throw std::runtime_error("FIRE exited unsuccessfully (status "+std::to_string(status)+"); inspect run.log");
    result.reductions=import_tables(result.directory/(modular_point?"result_"+probe_suffix+".tables":"result.tables"),field,basis.denominators.size(),&symbols);
    for(const auto& a:requests)if(!result.reductions.count(a))throw std::runtime_error("FIRE omitted requested reduction");
    write("success.receipt",modular_point?"DiffExp3.NativeFIREModularSample/v1\n":"DiffExp3.NativeFIRECompleted/v1\n");
    result.success=true;
  } catch(const std::exception& e){if(owned_pid>0){kill(-owned_pid,SIGKILL);int status=0;while(waitpid(owned_pid,&status,0)<0&&errno==EINTR){}}result.reason=e.what();if(!result.directory.empty())result.reason+="; log: "+(result.directory/"run.log").string();}
  return result;
}
} // namespace diffexp::fire
