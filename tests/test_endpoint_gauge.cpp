#include "diffexp/affine_operator.hpp"
#include <iostream>
using namespace diffexp;
using B=kernel::ComplexBall;
void check(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
int main(){try {
  B::set_precision(256);
  ExactField field({"x","eps"});
  auto q=[&](const char* value){return Exact(field,value);};
  auto series=AffineFrobeniusSeries::prepare({{q("eps/x"),q("1/x")},{q("0"),q("0")}},0,1,4);
  std::vector<std::vector<Exact>> transform{{q("1"),q("1")},{q("1"),q("2")}};
  std::vector<std::vector<Exact>> diagonal{{q("eps^3"),q("0")},{q("0"),q("eps")}};
  auto physical_gauge=fuchsify::detail::multiply(diagonal,transform);
  auto physical=series.project(physical_gauge),gauged=series.project(transform);
  auto row=fuchsify::detail::multiply(std::vector<std::vector<Exact>>{{q("1"),q("1/eps")}},physical_gauge);
  auto functional=series.dr_integral_from_zero(series.project(row));
  const auto point=B::from_strings("1/8");
  affine_operator::Options settings;settings.max_epsilon_depth=64;
  auto original=affine_operator::prepare(series,physical,point,10,settings);
  auto cancelled=affine_operator::prepare(series,gauged,point,10,settings);
  check(original.success(),original.reason);check(cancelled.success(),cancelled.reason);
  check(original.determinant_valuation==cancelled.determinant_valuation+4,"exact diagonal determinant shift");
  check(cancelled.row_pole_loss<original.row_pole_loss,"cancelling the diagonal must reduce artificial row normalization loss");
  auto a=affine_operator::compose(original,series,functional,point,4,settings);
  auto b=affine_operator::compose(cancelled,series,functional,point,4,settings);
  check(a.success(),a.reason);check(b.success(),b.reason);
  const int shifts[]{3,1};
  const auto coefficient=[](const auto& matrix,unsigned j,int k){
    if(k<matrix.low)return B(0);
    if(k>matrix.high)throw std::runtime_error("unknown Laurent coefficient in test");
    return matrix.coefficients[0][j][k-matrix.low];
  };
  for(unsigned j=0;j<2;++j)for(int k=std::min(a.matrix.low+1,b.matrix.low);k<=4;++k)
    check(acb_overlaps(coefficient(a.matrix,j,k-shifts[j]).raw(),coefficient(b.matrix,j,k).raw()),"P*(D*T*F)^-1*D must equal P*(T*F)^-1 coefficient by coefficient");
  std::cout<<"Exact endpoint epsilon-gauge cancellation preserves Laurent operator and reduces normalization loss\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
