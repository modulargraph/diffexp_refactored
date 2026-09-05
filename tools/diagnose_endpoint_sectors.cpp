#include "diffexp/level_cache.hpp"
#include "diffexp/recursion_graph.hpp"
#include "diffexp/cached_affine.hpp"
#include "diffexp/epsilon_gauge.hpp"
#include "diffexp/causal.hpp"
#include "diffexp/henn_boundary.hpp"
#include <iostream>
#include <chrono>
// Read-only diagnostic of fixed divergent sectors in a saved endpoint. It
// reports coefficients, never a numerical integral or a new exact certificate.
int main(int argc,char** argv){using namespace diffexp;namespace json=boost::json;try{
 if(argc!=3 && argc!=4)throw std::invalid_argument("usage: diffexp_diagnose_endpoint_sectors CACHE HENN_BASIS_FILE [full-projection]");
 const bool full=argc==4 && std::string(argv[3])=="full-projection";
 if(argc==4&&!full)throw std::invalid_argument("unknown endpoint diagnostic mode");
 auto canonical=henn::read_x0(argv[2]);recursion::Options options;options.anchors=causal::henn_anchors();options.anchors.front()=Rational("1/3");
 artifacts::Store graph_store(argv[1]);
 auto provider=[&](const auto& basis,const auto& dimension,const auto& field,auto parameter,const auto& sources,const auto& budget){
   level::Provider no_reduction=[](const auto&,const auto&)->fire::Result{throw std::runtime_error("diagnostic refuses uncached reduction");};
   return cached_level::prepare(graph_store,basis,dimension,field,parameter,sources,budget,no_reduction).result;
 };
 auto graph=recursion::prepare(feynman::example_family("henn_double_pentagon_x0"),canonical.scalar_targets,options,provider);
 const auto& node=graph.nodes.front();auto gauge=epsilon_diagonal_gauge(node.closure.matrix,1);
 fuchsify::Options fo;fo.max_dimension=256;auto regular=fuchsify::prepare(gauge.matrix,0,fo);if(!regular.success)throw std::runtime_error(regular.reason);
 auto diagonal=fuchsify::detail::identity(gauge.matrix.size(),graph.dimension);
 for(unsigned j=0;j<diagonal.size();++j)diagonal[j][j]=graph.dimension.variable(1).pow(gauge.shifts[j]);
 auto physical=fuchsify::detail::multiply(diagonal,regular.transform);
 AffineFrobeniusSeries::Options so;so.max_dimension=256;const unsigned n=32,d=physical.size();
 auto base=cached_affine::identity(regular.matrix,0,1,n,so),key=cached_affine::manifest_identity(base,d);
 artifacts::Store endpoints(std::filesystem::path(argv[1])/"affine-series");
 auto record=endpoints.lookup(key,Demand{0,0,n,64,0},{"exact",cached_affine::detail::verifier,cached_affine::detail::scope});
 if(!record)throw std::runtime_error("verified column manifest missing; diagnostic refuses recurrence");
 if(full) {
   const auto start=std::chrono::steady_clock::now();unsigned loaded=0;
   auto load=[&](unsigned c)->std::optional<std::vector<AffineFrobeniusSeries::Term>>{
     auto column=endpoints.read(cached_affine::column_identity(base,c),record->parent_content_ids.at(c));++loaded;
     return cached_affine::detail::decode_column(column.payload,graph.dimension,d,n,c,so.max_terms);
   };
   auto series=cached_affine::StateAccess::build_unverified_columns(regular.matrix,0,1,n,so,load,
     [](unsigned,std::span<const AffineFrobeniusSeries::Term>){throw std::runtime_error("diagnostic refuses any newly prepared column");});
   // This read-only diagnostic deliberately does not repeat the ODE verifier
   // and must never supply a published integral or exact certificate.
   auto frame=series.project(regular.transform);json::array operations;
   const auto x=graph.dimension.variable(0);unsigned conditions=0;
   for(unsigned i=0;i<node.operations.size();++i){const auto& op=node.operations[i];const bool integral=op.operation==feynman::Operation::BetaIntegral;
     if(!integral&&op.operation!=feynman::Operation::LowerLimit)continue;
     const auto begin=std::chrono::steady_clock::now();auto row=node.observable_rows[i];
     if(integral)for(auto& v:row)v=v*x.constant(op.normalization)*x.pow(op.left_power)*(x.constant(1)-x).pow(op.right_power);
     row=fuchsify::detail::multiply({row},physical)[0];
     auto domain=integral?series.dr_domain(series.project(row),true):series.project_endpoint_domain({row});
     if(integral)(void)series.dr_integral_from_zero(domain.admissible);else (void)series.dr_endpoint_constant(domain.admissible);
     conditions+=domain.zero_constraints.size();const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
     std::size_t terms=domain.admissible.terms.size();for(const auto& c:domain.zero_constraints)for(const auto& v:c.coefficients)terms+=!v.is_zero();
     operations.push_back(json::object{{"index",i},{"scope",integral?"full series":"endpoint limit"},{"terms",terms},{"constraints",domain.zero_constraints.size()},{"seconds",seconds}});
     std::cerr<<"projected operation "<<i<<" terms "<<terms<<" constraints "<<domain.zero_constraints.size()<<" seconds "<<seconds<<'\n';
   }
   std::cout<<json::serialize(json::object{{"status","diagnostic pass"},{"scope","saved Henn lower endpoint observable projection and domain construction"},
     {"columns_loaded",loaded},{"columns_prepared",0},{"ode_verification_repeated",false},{"numerical_result_computed",false},
     {"exact_certificate_published",false},{"seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()},
     {"constraints",conditions},{"operations",operations}})<<'\n';return 0;
 }
 std::vector<AffineFrobeniusSeries::Term> terms;
 for(unsigned c=0;c<d;++c){auto column=endpoints.read(cached_affine::column_identity(base,c),record->parent_content_ids[c]);
   auto decoded=cached_affine::detail::decode_column(column.payload,graph.dimension,d,n,c,so.max_terms);
   for(auto& t:decoded)if(t.slope.is_zero())terms.push_back(std::move(t));}
 const auto z=graph.dimension.constant(0),x=z.variable(0);json::array failures;unsigned lower=0,beta=0;
 for(unsigned i=0;i<node.operations.size();++i){const auto& op=node.operations[i];const bool integral=op.operation==feynman::Operation::BetaIntegral;
   if(!integral&&op.operation!=feynman::Operation::LowerLimit)continue;
   integral?++beta:++lower;auto row=node.observable_rows[i];
   if(integral)for(auto& v:row)v=v*x.constant(op.normalization)*x.pow(op.left_power)*(x.constant(1)-x).pow(op.right_power);
   row=fuchsify::detail::multiply({row},physical)[0];
   using Coordinate=std::tuple<unsigned,Rational,unsigned>;std::map<Coordinate,Exact> combined;
   for(unsigned j=0;j<d;++j)if(!row[j].is_zero()){
     auto val=fuchsify::detail::valuation(row[j],0);long needed=0;
     for(const auto& t:terms)if(t.row==j){auto cutoff=-(t.power+Rational(val));if(cutoff.sign()>=0&&cutoff.str().find('/')==std::string::npos)needed=std::max(needed,std::stol(cutoff.str()));}
     if(needed>n)throw std::runtime_error("diagnostic projection exceeds retained order");
     auto coefficients=affine_frobenius_detail::taylor(row[j]/fuchsify::detail::power(x,val),0,needed);
     for(const auto& t:terms)if(t.row==j)for(unsigned k=0;k<coefficients.size();++k){auto power=t.power+Rational(val)+Rational(k);
       if(power>Rational(0))break;
       if((integral?power<=Rational(-1):(power<Rational(0)||(power.is_zero()&&t.log_degree)))&&!coefficients[k].is_zero()){
         auto at=Coordinate{t.column,power,t.log_degree};auto [where,inserted]=combined.try_emplace(at,z);where->second=where->second+t.coefficient*coefficients[k];}
     }
   }
   json::array bad;for(const auto& [coordinate,value]:combined)if(!value.is_zero()){
     const auto& [c,power,log]=coordinate;bad.push_back(json::object{{"column",c},{"power",power.str()},{"log_degree",log},{"coefficient",value.str()}});}
   if(!bad.empty()){json::array target;for(auto index:node.requested[i])target.push_back(index);failures.push_back(json::object{{"operation_index",i},{"kind",integral?"beta":"lower limit"},{"target",target},{"fixed_sectors",bad}});std::cerr<<"operation "<<i<<" has "<<bad.size()<<" fixed divergent terms\n";}
 }
 std::cout<<json::serialize(json::object{{"status","diagnostic only"},{"endpoint_manifest_semantic_id",key.key()},{"endpoint_manifest_content_id",record->content_id},{"dimension",d},{"order",n},{"lower_limit_operations",lower},{"beta_operations",beta},{"offending_operations",failures},{"numerical_result_computed",false}})<<'\n';
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
