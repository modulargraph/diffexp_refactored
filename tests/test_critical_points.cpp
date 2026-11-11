#include "diffexp/critical_points.hpp"
#include "diffexp/families.hpp"
#include <iostream>
using namespace diffexp;
static void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
int main(){try {
  auto example=feynman::example_family("kite");
  auto family=example.momenta;family.lines.resize(example.physical_count);
  std::vector<std::vector<Rational>> lines(example.physical_count,std::vector<Rational>(example.physical_count,Rational(0)));
  for(unsigned i=0;i<example.physical_count;++i)lines[i][i]=Rational(1);
  const auto text=critical_points::script(family,lines,0,1);
  check(text.find("ring r=32003,(t,b0,b1,b2,b3),dp;")==0 && text.find("vdim(J)")!=std::string::npos,"likelihood ideal script");
  bool refused=false;
  try{(void)critical_points::script(family,lines,1,1);}catch(const std::invalid_argument&){refused=true;}
  check(refused,"a line cannot be merged with itself");
  if(!critical_points::available()){std::cout<<"Singular is not installed: critical point counts skipped\n";return 0;}
  const auto cache=std::filesystem::temp_directory_path()/("diffexp-critical-points-test-"+std::to_string(getpid()));
  const auto first=critical_points::count(family,lines,0,1,cache),again=critical_points::count(family,lines,0,1,cache);
  check(first && *first==16 && again==first,"the kite merged on its first pair has 16 masters");
  const auto other=critical_points::count(family,lines,0,2,cache);
  check(other && *other>*first,"a pair on different loop momenta gives a larger kite level");
  std::filesystem::remove_all(cache);
  std::cout<<"Critical point counts of merged families\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
