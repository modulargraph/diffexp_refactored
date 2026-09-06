#include "product_recurrence.hpp"
#include "stage_probe_io.hpp"
#include <fstream>
#include <chrono>
#include <iostream>
using namespace diffexp;
using namespace experiment;
namespace js=boost::json;
js::object read(const char* file) {
  std::ifstream in(file);if(!in)throw std::runtime_error("missing cache input");
  return js::parse(std::string(std::istreambuf_iterator<char>(in),{})).as_object();
}
void overlap(const Boundary& a,const Boundary& b) {
  if(a.size()!=b.size())throw std::runtime_error("shape mismatch");
  for(unsigned i=0;i<a.size();++i)for(unsigned k=0;k<a[i].size();++k)
    if(!acb_overlaps(a[i][k].raw(),b[i][k].raw()))throw std::runtime_error("retained-polynomial ball mismatch");
}
int main(int argc,char** argv) {
  if(argc!=4 || std::string(argv[1])!="--run") {
    std::cerr<<"Usage: cached_chart --run CLOSURE.json ROW_CHECKPOINT.json\n";return 2;
  }
  try {
    B::set_precision(384);const unsigned N=80;
    auto record=read(argv[2]),checkpoint=read(argv[3]);
    const auto& payload=record.at("payload");
    if(artifacts::detail::sha256(artifacts::detail::canonical(payload))!=
       artifacts::detail::string(record.at("certificate").as_object().at("payload_sha256")))
      throw std::runtime_error("closure payload checksum");
    if(artifacts::detail::sha256(artifacts::detail::canonical(checkpoint.at("payload")))!=
       artifacts::detail::string(checkpoint.at("sha256")))throw std::runtime_error("row checksum");
    auto rows=numerical_rows_io::read_rows(checkpoint.at("payload").as_object().at("rows"));
    if(rows.low!=0 || rows.coefficients.empty())throw std::runtime_error("unsupported boundary window");
    // A cached ball vector is used as common test data for the physical matrix.
    // It is NOT claimed to be the physical boundary or original chart replay:
    // the checkpoint does not encode the prepared stage gauge/path identity.
    Boundary boundary=rows.coefficients.front();
    const unsigned d=boundary.size(),high=rows.high;
    std::vector<std::string> names;
    for(const auto& v:record.at("identity").as_object().at("scientific_inputs").as_object().at("ordered_field_symbols").as_array())
      names.push_back(artifacts::detail::string(v));
    ExactField field(names);
    auto matrix=stage_probe::read_matrix(payload.as_object().at("matrix"),field,d,d);
    auto [xi,ei]=path_epsilon_variables(matrix[0][0]);
    std::vector<Exact> substitutions;for(const auto& name:names)substitutions.emplace_back(field,name);
    substitutions[xi]=substitutions[xi]+Exact(field,"1/2");
    std::vector<RationalLineEntry> translated,original;
    for(unsigned i=0;i<d;++i)for(unsigned j=0;j<d;++j) {
      auto cs=feynman::scalar_functional_detail::epsilon_series(matrix[i][j],ei,high);
      for(unsigned k=0;k<=high;++k)if(!cs[k].is_zero()) {
        original.push_back({i,j,k,cs[k]});
        translated.push_back({i,j,k,cs[k].substitute(substitutions)});
      }
    }
    auto t=std::chrono::steady_clock::now();auto c=prepare(translated,d,N,high);
    double prep=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
    if(!c.supported)throw std::runtime_error("cached matrix outside prototype scope");
    B h=polynomial_transport::detail::ball(Rational("1/64"));
    auto a=polynomial_transport::chart(c.baseline,boundary,B(0),h,N);
    auto b=chart(c,boundary,B(0),h,N);
    auto dense=horner(coefficients<B>(c,boundary,N,true),h);
    overlap(a,b);overlap(b,dense);
    // Cross-check the exact coordinate translation against production evaluation
    // at the original ordinary center, without replacing or rounding boundaries.
    auto unshifted=polynomial_transport::compile(original,d,N,high);
    overlap(b,polynomial_transport::chart(unshifted,boundary,B::from_strings("0.5"),h,N));
    std::size_t count=0,dense_count=0;
    coefficients<B>(c,boundary,N,false,&count);coefficients<B>(c,boundary,N,true,&dense_count);
    auto timed=[&](auto fn){auto start=std::chrono::steady_clock::now();for(unsigned i=0;i<5;++i)fn();return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/5;};
    double baseline=timed([&]{polynomial_transport::chart(c.baseline,boundary,B(0),h,N);});
    double candidate=timed([&]{chart(c,boundary,B(0),h,N);});
    std::cout<<js::serialize(js::object{{"case","cached sunrise matrix microbenchmark"},{"original_chart_replay",false},{"dimension",d},{"epsilon_high",high},{"working_bits",384},{"taylor_order",N},{"center","1/2"},{"step","1/64"},{"preparation_seconds",prep},{"warm_repetitions",5},{"baseline_seconds",baseline},{"candidate_seconds",candidate},{"product_operations",count},{"dense_operations",dense_count},{"validation","passed"},{"omitted_tails_certified",false}})<<'\n';
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
