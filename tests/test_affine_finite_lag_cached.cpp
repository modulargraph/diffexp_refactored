#include "diffexp/recursion_pipeline.hpp"
#include "diffexp/level_cache.hpp"
#include <chrono>
#include <iostream>
using namespace diffexp;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void equal(const AffineFrobeniusSeries& a,const AffineFrobeniusSeries& b){
  require(a.dimension()==b.dimension() && a.x_order()==b.x_order(),"dimensions/order changed");
  require(a.residue_frame()==b.residue_frame(),"residue frame changed");
  require(a.exponents().size()==b.exponents().size(),"spectrum size changed");
  for(unsigned i=0;i<a.exponents().size();++i)
    require(a.exponents()[i].power==b.exponents()[i].power && a.exponents()[i].slope==b.exponents()[i].slope,"affine spectrum changed");
  require(a.jordan_successors()==b.jordan_successors(),"Jordan successors changed");
  require(a.absolute_x_frontiers()==b.absolute_x_frontiers(),"absolute frontier changed");
  const auto &x=a.terms(),&y=b.terms();
  require(x.rows==y.rows && x.columns==y.columns && x.terms.size()==y.terms.size(),"expansion shape changed");
  require(x.coherent_x_frontier==y.coherent_x_frontier && !x.omitted_tail_certified && !y.omitted_tail_certified,"certificate status changed");
  require(x.wronskian_prefactor==y.wronskian_prefactor,"Wronskian certificate changed");
  for(unsigned i=0;i<x.terms.size();++i){const auto &u=x.terms[i],&v=y.terms[i];
    require(u.row==v.row && u.column==v.column && u.log_degree==v.log_degree && u.power==v.power && u.slope==v.slope && u.coefficient==v.coefficient,"exact retained coefficient changed");}
}

int main(int argc,char** argv){try{
 if(argc<2)throw std::invalid_argument("supply verified exact cache directory, optional family and order");
 artifacts::Store store(argv[1]);recursion::Options graph_options;
 auto provider=[&](const auto& basis,const auto& dim,const auto& field,auto xi,const auto& sources,const auto& budget){auto reject=[](const std::vector<ibp::Integral>&,const fire::Options&)->fire::Result{throw std::runtime_error("cache miss; FIRE prohibited");};auto hit=cached_level::prepare(store,basis,dim,field,xi,sources,budget,reject);if(!hit.cache_hit||!hit.result.success)throw std::runtime_error("cache unavailable");return hit.result;};
 auto graph=recursion::prepare(feynman::example_family(argc>2?argv[2]:"banana4"),{},graph_options,provider);
 unsigned order=argc>3?std::stoi(argv[3]):8;
 for(unsigned node_index=0;node_index<graph.nodes.size();++node_index){
 auto& node=graph.nodes[node_index];if(node.closure.matrix.empty())continue;
 auto epsg=epsilon_diagonal_gauge(node.closure.matrix,1);auto fuchs=fuchsify::prepare(epsg.matrix,0);if(!fuchs.success)throw std::runtime_error(fuchs.reason);
 AffineFrobeniusSeries::Options options;options.max_dimension=256;options.finite_lag_recurrence=true;options.finite_lag_cost_fallback=order>8;
 auto begin=std::chrono::steady_clock::now();
 auto finite=AffineFrobeniusSeries::prepare(fuchs.matrix,0,1,order,options);
 auto middle=std::chrono::steady_clock::now();
 std::cerr<<"node="<<node_index<<" d="<<finite.dimension()<<" N="<<order<<" finite_rows="<<finite.finite_lag_rows()<<" finite_seconds="<<std::chrono::duration<double>(middle-begin).count()<<std::endl;
 if(argc<=4 || std::string(argv[4])!="finite-only"){
 options.finite_lag_recurrence=false;
 auto old=AffineFrobeniusSeries::prepare(fuchs.matrix,0,1,order,options);
 auto end=std::chrono::steady_clock::now();equal(old,finite);
 std::cerr<<"exact equality; old_seconds="<<std::chrono::duration<double>(end-middle).count()<<std::endl;
 }
 }
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
