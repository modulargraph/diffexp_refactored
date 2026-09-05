#include "diffexp/fire.hpp"
#include "diffexp/families.hpp"
#include <iostream>
using namespace diffexp;
int main(int argc,char** argv) {
  if(argc>1&&std::string(argv[1])=="-c") {std::this_thread::sleep_for(std::chrono::seconds(5));return 0;}
  try {
    fire::SymbolMap symbol_map({"I","I2","devara","Q_mass"});
    const auto original="I+I2*devara-Q_mass";
    if(fire::rename_symbols(fire::rename_symbols(original,symbol_map.forward),symbol_map.reverse)!=original)
      throw std::runtime_error("FIRE temporary symbols must round-trip whole tokens without collisions");
    ExactField field({"d"});Exact d(field,"d");
    ibp::PropagatorBasis basis(ibp::quadratic_family(feynman::banana(1,{Rational(1),Rational(2)}),d));
    auto text=fire::start(basis,d);if(text.find("SBasis0C")==std::string::npos)throw std::runtime_error("missing IBPs");
    fire::Options options;
    auto missing=fire::reduce(basis,d,field,{{2,1}},options);
    if(missing.success||missing.reason.empty()||!missing.directory.empty())throw std::runtime_error("missing executable must fail locally");
    options.executable=std::filesystem::absolute(argv[0]);options.timeout_seconds=1;
    auto timeout=fire::reduce(basis,d,field,{{2,1}},options);
    if(timeout.success||timeout.reason.find("timeout")==std::string::npos)throw std::runtime_error("bounded process timeout failed");
    std::filesystem::remove_all(timeout.directory);
    options.timeout_seconds=60;
    if(argc<2){std::cout<<"native preparation passed; FIRE executable argument enables integration\n";return 0;}options.executable=argv[1];
    auto result=fire::reduce(basis,d,field,{{2,1},{1,2}},options);
    std::cout<<result.directory<<"\n";if(!result.success)throw std::runtime_error(result.reason);
    auto parallel_options=options;parallel_options.threads=4;parallel_options.simplifier_threads=1;
    auto parallel=fire::reduce(basis,d,field,{{2,1},{1,2}},parallel_options);
    if(!parallel.success||parallel.reductions!=result.reductions)throw std::runtime_error("parallel multi-sector bubble differs from serial exact reductions: "+parallel.reason);
    if(fire::read_text(parallel.directory/"job.config").find("#threads 4\n#fthreads 1\n")==std::string::npos)throw std::runtime_error("parallel FIRE config omitted worker settings");
    ibp::Generator generator(basis,d);ibp::ExactReducer reducer(d);
    for(int a=-1;a<=4;++a)for(int b=-1;b<=4;++b)for(auto& row:generator.relations({a,b}))reducer.insert(std::move(row));
    for(const auto& [a,row]:result.reductions) {
      ibp::Relation residual{{a,d.constant(1)}};ibp::add_scaled(residual,row,d.constant(-1));
      if(!reducer.reduce(residual).remainder.empty())throw std::runtime_error("FIRE disagrees with native exact IBPs");
    }
    if(result.reductions.at({2,1}).count({2,1}))throw std::runtime_error("dotted target was not reduced");
    std::cout<<"native FIRE cold massive bubble reduction verified\n";
    ExactField sf({"x","d"});Exact x(sf,"x"),sd(sf,"d");
    ibp::PropagatorBasis sb(ibp::merge(ibp::quadratic_family(feynman::banana(2,{Rational(1),Rational(1),Rational(1)}),x),0,1,x));
    auto zeros=fire::free_loop_sectors(sb);
    if(zeros.size()!=2 || std::find(zeros.begin(),zeros.end(),ibp::Integral{-1,1,-1,-1,-1})==zeros.end())throw std::runtime_error("free loop zero-sector proof");
    auto sr=fire::reduce(sb,sd,sf,{{3,1,0,0,0}},options);
    std::cout<<sr.directory<<"\n";if(!sr.success)throw std::runtime_error(sr.reason);
    ibp::Generator sg(sb,sd);ibp::ExactReducer se(x,300);
    ibp::for_each_seed(2,5,{1,1,50},[&](const ibp::Integral& seed){for(auto& row:sg.relations(seed))se.insert(std::move(row));});
    ibp::Relation residual{{{3,1,0,0,0},x.constant(1)}};ibp::add_scaled(residual,sr.reductions.at({3,1,0,0,0}),x.constant(-1));
    if(!se.reduce(residual).remainder.empty())throw std::runtime_error("sunrise FIRE disagrees with exact native IBPs");
    if(sr.reductions.at({3,1,0,0,0}).count({3,1,0,0,0}))throw std::runtime_error("sunrise target was not reduced");
    std::cout<<"native FIRE cold merged sunrise reduction verified\n";
    ExactField capital_field({"x","eps","I","Q_mass"});Exact eps(capital_field,"eps"),mass(capital_field,"Q_mass");
    auto dimension=eps.constant(4)-eps.constant(2)*eps;
    auto cf=ibp::quadratic_family(feynman::banana(1,{Rational(1),Rational(2)}),eps);cf.physical[0].constant=mass;
    ibp::PropagatorBasis cb(cf);
    auto cr=fire::reduce(cb,dimension,capital_field,{{2,1}},options);
    if(!cr.success)throw std::runtime_error(cr.reason);
    ibp::Generator cg(cb,dimension);ibp::ExactReducer ce(eps);
    for(int a=-1;a<=3;++a)for(int b=-1;b<=3;++b)for(auto& equation:cg.relations({a,b}))ce.insert(std::move(equation));
    ibp::Relation capital_residual{{{2,1},eps.constant(1)}};ibp::add_scaled(capital_residual,cr.reductions.at({2,1}),eps.constant(-1));
    if(!ce.reduce(capital_residual).remainder.empty()||cr.reductions.at({2,1}).count({2,1}))
      throw std::runtime_error("uppercase used and unused symbols must import in the original field");
    std::cout<<"native FIRE arbitrary exact-field symbols verified\n";
  }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
