// Bounded comparison of source-kernel arithmetic on one captured chart.
// This measures retained-polynomial work, not acceptance or an analytic tail.
#include "diffexp/adjoint_checkpoint.hpp"
#include <chrono>
#include <cmath>
#include <csignal>
#include <iostream>
#include <unistd.h>
using namespace diffexp;namespace rc=rational_circuit;namespace json=boost::json;
using B=Jet::Ball;using Clock=std::chrono::steady_clock;
std::string read(const char* path){std::ifstream f(path);if(!f)throw std::runtime_error("input open failed");return {std::istreambuf_iterator<char>(f),{}};}
double radius(const Boundary& b){double r=0;for(auto&row:b)for(auto&a:row){if(!a.is_finite())throw std::runtime_error("nonfinite result");r=std::max({r,mag_get_d(arb_radref(acb_realref(a.raw()))),mag_get_d(arb_radref(acb_imagref(a.raw())))});}return r;}
int main(int argc,char**argv){try{
 if(argc!=4)throw std::invalid_argument("usage: benchmark_circuit_source PROBLEM CHECKPOINT NEXT_PARAMETER");
 ::alarm(60);auto bytes=read(argv[1]),savedbytes=read(argv[2]);auto problem=json::parse(bytes).as_object(),envelope=json::parse(savedbytes).as_object();
 B::set_precision(artifacts::detail::integer(problem.at("working_bits")));
 auto saved=adjoint_checkpoint::detail::decode(envelope,artifacts::detail::string(envelope.at("payload").at("identity")));
 std::vector<std::string>names;for(auto&v:problem.at("variables").as_array())names.emplace_back(v.as_string());ExactField field(names);Exact zero(field);auto[xi,ei]=path_epsilon_variables(zero);
 std::vector<Exact>path;for(auto&v:problem.at("path").as_array())path.emplace_back(field,std::string(v.as_string()));
 unsigned leg=saved.leg;if(saved.parameter==1){leg++;saved.parameter=0;}if(leg+1>=path.size())throw std::invalid_argument("completed checkpoint");
 B center;acb_set_d(center.raw(),saved.parameter);B end;acb_set_d(end.raw(),std::stod(argv[3]));B step=end-center;
 if(!arb_is_positive(acb_realref(step.raw()))||!arb_le(acb_realref(end.raw()),acb_realref(B(1).raw())))throw std::invalid_argument("invalid next parameter");
 auto& rows=saved.rows;unsigned d=rows.columns(),R=rows.coefficients.size(),W=rows.high-rows.low+1,N=80,dimension=R*d+1;
 if(problem.at("matrix").as_array().size()!=d)throw std::invalid_argument("checkpoint dimension mismatch");
 if(auto p=problem.at("settings").as_object().if_contains("ordinary_order"))N=artifacts::detail::integer(*p);
 auto scale=path[leg+1]-path[leg];auto point=exact_point(zero,xi,path[leg]+scale*zero.variable(xi));std::vector<RationalLineEntry>entries;
 for(unsigned i=0;i<d;++i)for(unsigned j=0;j<d;++j){Exact a(field,std::string(problem.at("matrix").as_array()[j].as_array()[i].as_string()));if(a.is_zero())continue;a=-scale*a.substitute(point);for(unsigned r=0;r<R;++r)entries.push_back({r*d+i,r*d+j,0,a});}
 for(auto&r:problem.at("forcing").as_array())for(auto&a:r.as_array())if(!Exact(field,std::string(a.as_string())).is_zero())throw std::invalid_argument("nonzero forcing not supported by benchmark");
 Boundary input(dimension,std::vector<B>(W));input.back()[0]=B(1);for(unsigned r=0;r<R;++r)for(unsigned i=0;i<d;++i)input[r*d+i]=rows.coefficients[r][i];Boundary midpoint=input;for(auto&r:midpoint)for(auto&a:r)acb_get_mid(a.raw(),a.raw());
 std::vector<Boundary>baseline;json::array variants;
 for(bool grouped:{false,true}){
  rc::Options options;options.max_terms=2000000;options.max_operations=1000000000;options.grouped_dot=grouped;
  auto compiled=rc::compile(entries,dimension,xi,ei,W-1,N,options);auto prepared=rc::prepare(compiled,center,N,W);
  double elapsed=0;json::array passes;
  for(unsigned p=0;p<2;++p){auto start=Clock::now();rc::Statistics stats;auto value=rc::chart(prepared,p?midpoint:input,step,&stats);auto seconds=std::chrono::duration<double>(Clock::now()-start).count();elapsed+=seconds;
   json::object pass{{"input",p?"exact_midpoint":"saved_uncertain"},{"seconds",seconds},{"addmul_terms",stats.addmul_operations},{"scalar_operations",stats.scalar_operations},{"dot_products",stats.dot_products},{"live_cells",stats.live_cells},{"maximum_radius",radius(value)},{"quality_log2",std::log2(adjoint_detail::enclosure_quality(value).approximate_upper())}};
   if(!grouped)baseline.push_back(value);else{std::size_t disjoint=0;double distance=0,max_ratio=0;for(unsigned i=0;i<dimension;++i)for(unsigned k=0;k<W;++k){auto&a=value[i][k];auto&b=baseline[p][i][k];disjoint+=!acb_overlaps(a.raw(),b.raw());for(unsigned component=0;component<2;++component){auto aa=component?acb_imagref(a.raw()):acb_realref(a.raw());auto bb=component?acb_imagref(b.raw()):acb_realref(b.raw());arf_t delta;arf_init(delta);arf_sub(delta,arb_midref(aa),arb_midref(bb),B::precision(),ARF_RND_NEAR);distance=std::max(distance,std::abs(arf_get_d(delta,ARF_RND_NEAR)));arf_clear(delta);double rb=mag_get_d(arb_radref(bb)),ra=mag_get_d(arb_radref(aa));if(rb>0)max_ratio=std::max(max_ratio,ra/rb);else if(ra>0)max_ratio=INFINITY;}}
    pass["disjoint_coefficients"]=disjoint;pass["maximum_midpoint_distance"]=distance;pass["maximum_per_component_radius_ratio"]=std::isfinite(max_ratio)?json::value(max_ratio):json::value("inf");pass["maximum_radius_ratio"]=radius(value)/radius(baseline[p]);if(disjoint)throw std::runtime_error("grouped/scalar disjoint enclosures");}
   passes.push_back(std::move(pass));
  }
  variants.push_back(json::object{{"grouped_dot",grouped},{"total_seconds",elapsed},{"passes",passes},{"dot_capacity",prepared.data().dot_capacity}});
 }
 std::cout<<json::serialize(json::object{{"schema","DiffExp21.CircuitSourceBenchmark/v1"},{"problem_sha256",artifacts::detail::sha256(bytes)},{"checkpoint_sha256",artifacts::detail::sha256(savedbytes)},{"bits",B::precision()},{"order",N},{"width",W},{"dimension",dimension},{"from",saved.parameter},{"to",std::stod(argv[3])},{"variants",variants},{"omitted_tail_certified",false}})<<'\n';
 ::alarm(0);
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
