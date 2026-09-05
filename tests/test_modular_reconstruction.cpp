#include "diffexp/modular_reconstruction.hpp"
#include <iostream>
using namespace diffexp;
int main(){try{
  ExactField f({"x","eps","unused"});Exact x(f,"x"),e(f,"eps");
  std::vector<Exact> functions={x.constant(0),x.constant(1),(x*x+e*x+x.constant(3))/(x*e+x.constant(7)),(x.constant(2)*e+x.constant(1))/(x*x),x.parse("(98765432109876543210987654321*x+eps)/(1234567890123456789*eps+x)")};
  std::vector<std::vector<modular::Sample>> samples(3);
  for(unsigned pi=1;pi<=3;++pi)for(unsigned i=0;i<40;++i){auto p=modular::point(3,pi,i);modular::Sample s{{p[0],p[1]}, {}};for(const auto& fn:functions)s.coefficients.push_back(modular::evaluate(fn,p,modular::prime(pi)));samples[pi-1].push_back(s);}
  for(std::size_t k=0;k<functions.size();++k){auto image=modular::discover(samples[0],k,2,6,modular::prime(1));if(!image)throw std::runtime_error("shape discovery failed");modular::Lift lift(*image,modular::prime(1));
    for(unsigned pi=2;pi<=3;++pi){auto next=modular::fit(samples[pi-1],k,image->ansatz,modular::prime(pi),true);if(!next)throw std::runtime_error("independent image failed");lift.append(*next,modular::prime(pi));}
    auto reconstructed=lift.reconstruct(x,{0,1});if(!reconstructed||*reconstructed!=functions[k])throw std::runtime_error("rational reconstruction differs from independent exact function");}
  auto corrupt=samples[0];corrupt.back().coefficients[2]^=1;
  if(modular::discover(corrupt,2,2,2,modular::prime(1)))throw std::runtime_error("corrupted held-out point accepted");
  bool pole=false;try{modular::evaluate(x/e,{1,0,3},modular::prime(1));}catch(const std::domain_error&){pole=true;}
  if(!pole)throw std::runtime_error("modular pole not rejected");
  std::cout<<"native multivariate modular interpolation, multi-prime CRT, large rational coefficients, zero-constant denominator and held-out corruption checks passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
