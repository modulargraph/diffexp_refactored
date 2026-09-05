// Opt-in: test_original_banana_unequal_routes ORIGINAL_EXAMPLE_DIRECTORY
// Fixed original native settings N64/bits384/epsilon0..7; no accuracy retries.
#include "diffexp/original_banana_unequal.hpp"
#include <boost/json.hpp>
#include <iostream>
int main(int argc,char** argv) {
 using namespace diffexp;using B=Jet::Ball;B::set_precision(384);
 try {
  if(argc!=2)throw std::invalid_argument("requires original-example directory");
  bool rejected=false;try{original_banana_unequal_route("unsupported");}catch(const std::invalid_argument&){rejected=true;}
  if(!rejected)throw std::runtime_error("route parser accepted unsupported route");
  const auto start=std::chrono::steady_clock::now();auto elapsed=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
  auto progress=[&](const std::string& s){std::cerr<<elapsed()<<"s "<<s<<'\n'<<std::flush;};
  auto first=original_banana_unequal(argv[1],original_banana_unequal_route("momentum-first"),progress);auto first_seconds=elapsed();
  auto second=original_banana_unequal(argv[1],original_banana_unequal_route("mass-first"),progress);
  auto reference=read_boundary(data::read_file(std::string(argv[1])+"/Reference/BananaUnequalMassAt50.m"),15,4,384);
  auto upper=[](const B& a){arf_t u;arf_init(u);acb_get_abs_ubound_arf(u,a.raw(),384);double r=arf_get_d(u,ARF_RND_CEIL);arf_clear(u);return r;};
  boost::json::array route_errors,paper_first,paper_second;bool pass=true;
  for(unsigned k=0;k<8;++k){double error=0;for(unsigned i=0;i<15;++i)error=std::max(error,upper(first.values[i][k]-second.values[i][k]));route_errors.push_back(error);pass=pass&&std::isfinite(error)&&error<1e-20;}
  for(unsigned k=0;k<5;++k){double a=0,b=0;for(unsigned i=0;i<15;++i){a=std::max(a,upper(first.values[i][k]-reference[i][k]));b=std::max(b,upper(second.values[i][k]-reference[i][k]));}paper_first.push_back(a);paper_second.push_back(b);pass=pass&&std::isfinite(a)&&std::isfinite(b)&&a<1e-10&&b<1e-10;}
  boost::json::object out{{"schema","DiffExp.OriginalUnequalBananaRoutes/v1"},{"status",pass?"pass":"fail"},{"masters",15},{"evolved_epsilon_low",0},{"evolved_epsilon_high",7},{"route_agreement_absolute_threshold",1e-20},{"route_max_difference_by_epsilon",route_errors},{"published_reference_window",boost::json::array{0,4}},{"published_reference_absolute_threshold",1e-10},{"published_max_difference_momentum_first",paper_first},{"published_max_difference_mass_first",paper_second},{"taylor_order",64},{"working_bits",384},{"momentum_first_seconds",first_seconds},{"mass_first_seconds",elapsed()-first_seconds},{"total_seconds",elapsed()},{"full_tails_certified",false},{"agreement_is_independent_reference_for_epsilon5_to7",false},{"shared_external_seed","Reference/BananaBoundaryAtMinusOneEps7.m"},{"primary_route_provenance","Docs/OriginalDiffExpExamples.md: mass deformation at p_squared=1/2; momentum1/2->25+20i->50"},{"seed_connector","-1 -> -1+i ->1/2+i ->1/2 at equal masses; explicitly supplied upper-side connector"}};
  std::cout<<boost::json::serialize(out)<<'\n';return pass?0:1;
 }catch(const std::exception& e){std::cout<<boost::json::serialize(boost::json::object{{"schema","DiffExp.OriginalUnequalBananaRoutes/v1"},{"status","error"},{"error",e.what()}})<<'\n';return 1;}
}
