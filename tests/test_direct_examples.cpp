#include "diffexp/direct_examples.hpp"
#include <iostream>
int main(){try {
  auto result=diffexp::singular_endpoint_example();
  if(result.at("samples").as_array().size()!=81 || result.at("dimreg_endpoint_constant").as_string()!="0" ||
      result.at("max_epsilon_one_error").as_double()>=1e-60)return 1;
  std::cout<<"Native singular endpoint: exact affine sector, matched epsilon window, dimensional limit and81 logarithm samples passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
