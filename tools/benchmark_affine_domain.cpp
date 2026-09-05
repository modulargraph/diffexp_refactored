#include "diffexp/cached_affine.hpp"
#include "diffexp/epsilon_gauge.hpp"
#include <chrono>
#include <iostream>
int main(int argc,char** argv){using namespace diffexp;try{
 if(argc!=3)throw std::invalid_argument("usage: diffexp_benchmark_affine_domain MATRIX_JSON ORDER");
 std::ifstream in(argv[1]);std::string bytes((std::istreambuf_iterator<char>(in)),{});auto input=boost::json::parse(bytes).as_object();
 ExactField field({"x","eps","I"});Exact zero(field,0);AffineFrobeniusSeries::Matrix matrix;
 for(const auto& row:input.at("matrix").as_array()){matrix.emplace_back();for(const auto& value:row.as_array())matrix.back().emplace_back(field,std::string(value.as_string()));}
 const unsigned order=std::stoul(argv[2]);if(order>32)throw std::invalid_argument("bounded benchmark order exceeds32");
 auto gauge=epsilon_diagonal_gauge(matrix,1);fuchsify::Options fuchs;fuchs.max_dimension=256;
 auto regular=fuchsify::prepare(gauge.matrix,0,fuchs);if(!regular.success)throw std::runtime_error(regular.reason);
 AffineFrobeniusSeries::Options opts;opts.max_dimension=256;opts.univariate_epsilon_recurrence=true;
 auto begin=std::chrono::steady_clock::now();std::cerr<<"univariate begin\n";
 auto fast=AffineFrobeniusSeries::prepare(regular.matrix,0,1,order,opts);
 auto fast_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();std::cerr<<"univariate completed in "<<fast_seconds<<" seconds\n";
 begin=std::chrono::steady_clock::now();opts.univariate_epsilon_recurrence=false;std::cerr<<"multivariate reference begin\n";
 auto original=AffineFrobeniusSeries::prepare(regular.matrix,0,1,order,opts);auto original_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
 auto payload=cached_affine::detail::payload(fast);const bool equal=payload==cached_affine::detail::payload(original);
 if(!equal)throw std::runtime_error("univariate and multivariate retained coefficients differ");
 auto state=cached_affine::detail::decode(payload,zero,matrix.size(),order,opts);cached_affine::detail::verify(regular.matrix,0,1,order,opts,state);
 std::cout<<boost::json::serialize(boost::json::object{{"status","pass"},{"dimension",matrix.size()},{"order",order},{"terms",fast.terms().terms.size()},{"payload_bytes",boost::json::serialize(payload).size()},{"univariate_seconds",fast_seconds},{"multivariate_seconds",original_seconds},{"speed_ratio",original_seconds/fast_seconds},{"exact_payload_equal",equal},{"independent_polynomial_ode_check",true}})<<'\n';
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
