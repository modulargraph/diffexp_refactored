#include "diffexp/cached_affine.hpp"
#include "diffexp/epsilon_gauge.hpp"
#include <chrono>
#include <iostream>
// Offline verifier of an existing column manifest. Never starts a recurrence.
int main(int argc,char** argv){using namespace diffexp;try{
 if(argc!=5)throw std::invalid_argument("usage: diffexp_verify_affine_cache MATRIX_JSON CACHE ORDER PRODUCT_BUDGET");
 std::ifstream in(argv[1]);std::string bytes((std::istreambuf_iterator<char>(in)),{});auto input=boost::json::parse(bytes).as_object();
 ExactField field({"x","eps","I"});AffineFrobeniusSeries::Matrix matrix;
 for(const auto& row:input.at("matrix").as_array()){matrix.emplace_back();for(const auto& value:row.as_array())matrix.back().emplace_back(field,std::string(value.as_string()));}
 unsigned order=std::stoul(argv[3]);if(order>128)throw std::invalid_argument("verification order exceeds128");
 auto gauge=epsilon_diagonal_gauge(matrix,1);fuchsify::Options fuchs;fuchs.max_dimension=256;
 auto regular=fuchsify::prepare(gauge.matrix,0,fuchs);if(!regular.success)throw std::runtime_error(regular.reason);
 AffineFrobeniusSeries::Options options;options.max_dimension=256;
 artifacts::Store store(argv[2]);auto base=cached_affine::identity(regular.matrix,0,1,order,options);
 auto key=cached_affine::manifest_identity(base,matrix.size());const Demand demand{0,0,order,64,0};
 if(!store.lookup(key,demand))throw std::runtime_error("required existing column manifest missing; refusing recurrence");
 cached_affine::VerificationLimits limits;limits.max_term_products=std::stoull(argv[4]);
 if(!limits.max_term_products || limits.max_term_products>200000000000ULL)throw std::invalid_argument("verification product budget outside1..200000000000");
 const auto start=std::chrono::steady_clock::now();std::size_t products=0;
 limits.column_progress=[&](unsigned c,unsigned d,std::size_t p){products=p;std::cerr<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"s verified column "<<c<<'/'<<d<<" cumulative product estimate "<<p<<'\n';};
 auto result=cached_affine::prepare(regular.matrix,0,1,order,options,store,limits);
 if(!result.cache_hit || result.columns_prepared)throw std::runtime_error("offline verifier unexpectedly prepared columns");
 std::cout<<boost::json::serialize(boost::json::object{{"status","pass"},{"dimension",matrix.size()},{"order",order},{"seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()},{"product_estimate",products},{"product_budget",limits.max_term_products},{"columns_reused",result.columns_reused},{"columns_prepared",result.columns_prepared},{"semantic_id",result.semantic_id},{"content_id",result.content_id}})<<'\n';
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
