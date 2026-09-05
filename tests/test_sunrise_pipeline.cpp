#include "diffexp/sunrise.hpp"
#include "diffexp/banana_oracle.hpp"
#include <iostream>
using namespace diffexp;
int main(){try {
  using B=Jet::Ball;B::set_precision(384);auto plan=feynman::prepare_sunrise();
  const Rational from("1/5"),to("1/3");
  auto source=feynman::sunrise_boundary(plan,from,4),target=feynman::sunrise_boundary(plan,to,4);
  if(source.epsilon_low!=target.epsilon_low)throw std::runtime_error("sunrise boundary windows disagree");
  auto transported=rational_line(feynman::sunrise_line(plan,from,to),source.values,80);
  double maximum=0;
  for(unsigned i=0;i<3;++i)for(unsigned k=0;k<transported[i].size();++k) {
    auto difference=transported[i][k]-target.values[i][k];arf_t upper;arf_init(upper);
    acb_get_abs_ubound_arf(upper,difference.raw(),384);maximum=std::max(maximum,arf_get_d(upper,ARF_RND_CEIL));arf_clear(upper);
  }
  if(!std::isfinite(maximum) || maximum>1e-32)throw std::runtime_error("native recursive sunrise boundaries disagree after IBP series transport");
  std::cout<<"Sunrise complete merged-level boundary/transport consistency, epsilon -2..4, discrepancy "<<maximum<<'\n';
  auto full=feynman::sunrise();
  oracle::BananaOptions options;options.target_bits=112;options.working_bits=224;
  auto reference=oracle::equal_banana_bessel(2,options);
  auto difference=full.at(0)-reference.value;arf_t error;arf_init(error);acb_get_abs_ubound_arf(error,difference.raw(),384);
  auto discrepancy=arf_get_d(error,ARF_RND_CEIL);arf_clear(error);
  if(!std::isfinite(discrepancy) || discrepancy>1e-20)throw std::runtime_error("full native sunrise recursion disagrees with independent Bessel oracle");
  for(int k=full.epsilon_low;k<0;++k) {
    arf_t pole;arf_init(pole);acb_get_abs_ubound_arf(pole,full.at(k).raw(),384);
    auto size=arf_get_d(pole,ARF_RND_CEIL);arf_clear(pole);
    if(size>1e-20)throw std::runtime_error("full sunrise has an uncancelled Laurent pole");
  }
  std::cout<<"Full native sunrise finite-part Bessel discrepancy "<<discrepancy<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
