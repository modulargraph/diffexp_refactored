#include "diffexp/modular_reconstruction.hpp"
#include <iostream>
using namespace diffexp;
int main(){try{
  ExactField f({"x","eps","unused"});Exact x(f,"x"),e(f,"eps");
  std::vector<Exact> functions={x.constant(0),x.constant(1),(x*x+e*x+x.constant(3))/(x*e+x.constant(7)),(x.constant(2)*e+x.constant(1))/(x*x),x.parse("(98765432109876543210987654321*x+eps)/(1234567890123456789*eps+x)")};
  std::vector<std::vector<modular::Sample>> samples(3);
  for(unsigned pi=1;pi<=3;++pi)for(unsigned i=0;i<40;++i){auto p=modular::point(3,pi,i);modular::Sample s{{p[0],p[1]}, {}};for(const auto& fn:functions)s.coefficients.push_back(modular::evaluate(fn,p,modular::prime(pi)));samples[pi-1].push_back(s);}
  for(bool compact:{false,true})for(std::size_t k=0;k<functions.size();++k){auto image=modular::discover(samples[0],k,2,6,modular::prime(1));if(!image)throw std::runtime_error("shape discovery failed");if(compact)image=modular::compact_support(*image);modular::Lift lift(*image,modular::prime(1));
    for(unsigned pi=2;pi<=3;++pi){auto next=modular::fit(samples[pi-1],k,image->ansatz,modular::prime(pi),true);if(!next)throw std::runtime_error("independent image failed");lift.append(*next,modular::prime(pi));}
    auto reconstructed=lift.reconstruct(x,{0,1});if(!reconstructed||*reconstructed!=functions[k])throw std::runtime_error("rational reconstruction differs from independent exact function");}
  auto corrupt=samples[0];corrupt.back().coefficients[2]^=1;
  if(modular::discover(corrupt,2,2,2,modular::prime(1)))throw std::runtime_error("corrupted held-out point accepted");
  bool pole=false;try{modular::evaluate(x/e,{1,0,3},modular::prime(1));}catch(const std::domain_error&){pole=true;}
  if(!pole)throw std::runtime_error("modular pole not rejected");
  // Coordinate slices find an anisotropic rational function without a dense
  // total-degree search. Nonzero valuations matter for x*eps denominators.
  auto anisotropic=x.parse("(x^6+x+eps^2+1)/(x^7*eps+x*eps)");
  auto prime=modular::prime(1);std::vector<modular::Image> slices;
  for(unsigned axis=0;axis<2;++axis){std::vector<modular::Sample> line;
    for(unsigned i=0;i<24;++i){auto p=modular::point(3,1,100+i);auto t=p[axis];p[1-axis]=37;line.push_back({{t},{modular::evaluate(anisotropic,p,prime)}});}
    auto shape=modular::discover(line,0,1,8,prime);if(!shape)throw std::runtime_error("univariate degree discovery failed");slices.push_back(*shape);
  }
  auto box=modular::degree_box(slices,8);if(!box||box->numerator.size()+box->denominator.size()>29)throw std::runtime_error("anisotropic degree box failed");
  std::vector<modular::Sample> grid;
  for(unsigned i=0;i<40;++i){auto p=modular::point(3,1,i);grid.push_back({{p[0],p[1]},{modular::evaluate(anisotropic,p,prime)}});}
  auto image=modular::fit_with_holdouts(grid,0,*box,prime);if(!image)throw std::runtime_error("degree box fit failed");
  for(unsigned i=0;i<10;++i){auto p=modular::point(3,1,600+i);if(modular::evaluate(*image,{p[0],p[1]},prime)!=std::optional<modular::Word>(modular::evaluate(anisotropic,p,prime)))throw std::runtime_error("degree box independent evaluation failed");}
  auto corrupt_box=grid;corrupt_box.back().coefficients[0]^=1;
  if(modular::fit_with_holdouts(corrupt_box,0,*box,prime))throw std::runtime_error("degree box accepted corrupted heldout");
  // A nongeneric y=0 slice erases x dependence. Its candidate must fail on
  // independent multivariate samples; the caller then uses dense discovery.
  auto cancellation=x.parse("(x*eps+1)/(x*eps+eps+1)");std::vector<modular::Sample> cancelled,other,checks;
  for(unsigned i=0;i<16;++i){auto p=modular::point(3,1,900+i);cancelled.push_back({{p[0]},{modular::evaluate(cancellation,{p[0],0,3},prime)}});other.push_back({{p[1]},{modular::evaluate(cancellation,{37,p[1],3},prime)}});checks.push_back({{p[0],p[1]},{modular::evaluate(cancellation,p,prime)}});}
  auto sx=modular::discover(cancelled,0,1,4,prime),sy=modular::discover(other,0,1,4,prime);
  if(!sx||!sy)throw std::runtime_error("cancellation setup failed");auto bad_box=modular::degree_box({*sx,*sy},4);
  if(!bad_box||modular::fit_with_holdouts(checks,0,*bad_box,prime))throw std::runtime_error("nongeneric degree slice was accepted");
  if(!modular::discover(checks,0,2,4,prime))throw std::runtime_error("dense fallback failed");
  std::cout<<"native multivariate modular interpolation, multi-prime CRT, large rational coefficients, zero-constant denominator and held-out corruption checks passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
