#include "diffexp/family_config.hpp"
#include "diffexp/cli.hpp"
#include "diffexp/banana_oracle.hpp"
#include "diffexp/box_triangle_oracle.hpp"
#include "diffexp/double_box_oracle.hpp"
#include <iostream>
using namespace diffexp;
namespace j=boost::json;
using B=kernel::ComplexBall;
int main(int argc,char**argv){try {
 if(argc!=4)throw std::invalid_argument("usage: check REFERENCE_CASE CONFIG RESULT");
 const auto started=std::chrono::steady_clock::now();std::string name=argv[1];
 auto cfg=family_config::load(argv[2]);
 if(!family_config::same_geometry(cfg.family,feynman::example_family(name)) || !cfg.integrals.empty())
  throw std::invalid_argument("reference geometry or target does not match configuration");
 std::ifstream f(argv[3]);auto result=j::parse(std::string(std::istreambuf_iterator<char>(f),{})).as_object();
 if(result.at("coefficients").as_array().size()!=1)throw std::invalid_argument("one scalar target required");
 B::set_precision(384);std::map<int,B> refs;std::string kind,threshold="1e-20";int first=0;bool reference_certified=false;
 if(name=="pentagon_massive") {
  refs.emplace(0,B::from_strings("1813378668630195764/100000000000000000000"));
  refs.emplace(1,B::from_strings("7613115416144053565/1000000000000000000000"));
  refs.emplace(2,B::from_strings("5214475578477681142/1000000000000000000000"));
  kind="Independent Feynman-parameter integration pins; comparison limited to 18 absolute digits";threshold="1e-18";
 } else if(name=="box_triangle") {
  auto ref=oracle::box_triangle_reference(256);for(int k=-4;k<=0;++k)refs.emplace(k,ref.at(k));first=-4;
  kind="Independent original-IBP and Mellin-Barnes/Cauchy reference; estimated reference error 1e-30";
 } else if(name=="double_box_planar") {
  auto ref=oracle::double_box_planar_reference(256);for(int k=-4;k<=0;++k)refs.emplace(k,ref.at(k));first=-4;
  kind="Independent analytic Smirnov planar double-box Laurent coefficients";reference_certified=true;
 } else if(name=="banana_unequal") {
  std::vector<Rational> masses;for(auto& line:cfg.family.momenta.lines)masses.push_back(line.mass_squared);
  oracle::BananaOptions opt;opt.target_bits=100;opt.working_bits=256;auto ref=oracle::banana_bessel(masses,opt);refs.emplace(0,ref.value);
  kind="Independent certified coordinate-space Bessel integral at epsilon zero";reference_certified=true;
 } else throw std::invalid_argument("unsupported reference case");
 B::set_precision(384);const auto limit=NativeTailMagnitude::lower_abs(B::from_strings(threshold));
 int low=result.at("epsilon_low").as_int64();auto& row=result.at("coefficients").as_array()[0].as_array();
 auto coefficient=[&](int k){if(k<low)return B(0);auto& v=row.at(k-low).as_object();return B::from_strings(std::string(v.at("real").as_string()),std::string(v.at("imaginary").as_string()));};
 bool pass=true;j::array comparisons,poles;
 for(auto& [k,ref]:refs){auto value=coefficient(k),delta=value-ref;bool good=delta.is_finite()&&NativeTailMagnitude::upper_abs(delta)<=limit;pass&=good;comparisons.push_back(j::object{{"epsilon_order",k},{"actual",ball_json(value,80)},{"reference",ball_json(ref,80)},{"difference",ball_json(delta,50)},{"pass",good}});}
 for(int k=low;k<first;++k){auto v=coefficient(k);bool good=v.is_finite()&&NativeTailMagnitude::upper_abs(v)<=NativeTailMagnitude::lower_abs(B::from_strings("1e-20"));pass&=good;poles.push_back(j::object{{"epsilon_order",k},{"coefficient",ball_json(v,40)},{"pass",good}});}
 std::cout<<j::serialize(j::object{{"status",pass?"pass":"fail"},{"reference_kind",kind},{"reference_certified",reference_certified},{"absolute_threshold",threshold},{"comparisons",comparisons},{"forbidden_poles",poles},{"omitted_tails_certified",result.at("omitted_tails_certified")},{"reference_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()}})<<'\n';return pass?0:1;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
