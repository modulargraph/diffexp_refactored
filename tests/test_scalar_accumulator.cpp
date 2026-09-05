#include "diffexp/scalar_functional.hpp"
#include <iostream>
using namespace diffexp;using B=Jet::Ball;
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
int main(){try {
  B::set_precision(384);ExactField field({"x","eps","I"});Exact x(field,"x"),eps(field,"eps");
  feynman::CertifiedDeepestBetaOptions options;options.working_bits=384;options.requested_digits=28;options.taylor_order=80;options.f_rim=-1;
  for(const auto& f:{"x^2*(99+3959*x)/236196","x^2*(6+326*x)/19683"}) {
    feynman::Symanzik geometry{Exact(field,"x*(81-7*x)/729"),Exact(field,f)};
    // This exact gamma ratio equals J_2/J_1 for L=2,d=4-2eps. The
    // first representation has a large one-way forcing; the second does not.
    auto weight=(eps.constant(2)*eps-eps.constant(3))*geometry.U/geometry.F;
    auto weighted=feynman::certified_scalar_functional(geometry,2,4,1,weight,2,options);
    auto unit=feynman::certified_scalar_functional(geometry,2,4,2,x.constant(1),2,options);
    require(weighted.epsilon_low==unit.epsilon_low && weighted.epsilon_low==-1,"scalar gamma ratio Laurent window");
    for(int k=-1;k<=2;++k)require(acb_overlaps(weighted.at(k).raw(),unit.at(k).raw()),"large-forcing representation must agree with exact shifted-gamma identity");
    require(weighted.taylor_tail_certified && unit.taylor_tail_certified &&
      native_enclosure_meets_tolerance({weighted.coefficients,unit.coefficients},"1/10000000000000000000000000000"),"doublebox/boxtriangle scalar leaves must certify28digits without increasing original resource settings");
  }
  std::cout<<"Scalar accumulator: doublebox/boxtriangle large rational forcing versus exact gamma-ratio identity certified at original N80/384bits/28digits\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
