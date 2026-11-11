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
  // A frame whose leading epsilon coefficient is invertible at the point needs
  // no row gauge: the certified ball inverse replaces the exact elimination
  // and both give the same inverse, order by order.
  auto regular=AffineFrobeniusSeries::prepare({{q("eps/x"),q("1")},{q("0"),q("2*eps/x")}},0,1,6);
  auto frame=regular.project(transform);
  auto exact_settings=settings;exact_settings.numeric_operator_convolution=true;
  auto numeric_settings=exact_settings;numeric_settings.numeric_leading_gauge=true;
  auto by_elimination=affine_operator::prepare(regular,frame,point,3,exact_settings);
  auto by_inverse=affine_operator::prepare(regular,frame,point,3,numeric_settings);
  check(by_elimination.success(),by_elimination.reason);check(by_inverse.success(),by_inverse.reason);
  check(!by_elimination.recurrence->numeric_gauge && by_inverse.recurrence->numeric_gauge,"numeric row gauge selection");
  check(by_elimination.row_pole_loss==0 && by_inverse.row_pole_loss==0,"a regular leading frame needs no shear");
  check(by_elimination.inverse.low==by_inverse.inverse.low && by_elimination.row_lower_bounds==by_inverse.row_lower_bounds,"numeric row gauge window");
  for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j)for(int k=by_inverse.inverse.low;k<=3;++k) {
    const auto& u=by_elimination.inverse.coefficients[i][j][k-by_elimination.inverse.low];
    const auto& v=by_inverse.inverse.coefficients[i][j][k-by_inverse.inverse.low];
    check(acb_overlaps(u.raw(),v.raw()) && v.is_finite(),"numeric row gauge inverse");
    const bool tight=mag_cmp_2exp_si(arb_radref(acb_realref(v.raw())),-200)<0 &&
      mag_cmp_2exp_si(arb_radref(acb_imagref(v.raw())),-200)<0;
    check(tight,"numeric row gauge inverse is not accurate");
  }
  // The exact zeros of the inverse survive: rows built on it must not
  // carry columns that are not theirs.
  std::vector<std::vector<Exact>> unit{{q("1"),q("0")},{q("0"),q("1")}};
  auto triangular=regular.project(unit);
  auto exact_triangular=affine_operator::prepare(regular,triangular,point,2,exact_settings);
  auto numeric_triangular=affine_operator::prepare(regular,triangular,point,2,numeric_settings);
  check(exact_triangular.success() && numeric_triangular.success() && numeric_triangular.recurrence->numeric_gauge,"numeric row gauge on a triangular frame");
  unsigned zeros=0;
  for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j)for(int k=numeric_triangular.inverse.low;k<=2;++k) {
    const auto& u=exact_triangular.inverse.coefficients[i][j][k-exact_triangular.inverse.low];
    const auto& v=numeric_triangular.inverse.coefficients[i][j][k-numeric_triangular.inverse.low];
    check(acb_overlaps(u.raw(),v.raw()),"numeric row gauge triangular inverse");
    if(k==0 && u.is_zero()){++zeros;check(v.is_zero(),"numeric row gauge lost an exact zero of the inverse");}
  }
  check(zeros==1,"the triangular frame has one exact zero in its leading inverse");
  // The pole-bearing frame above is singular at leading order: the numeric
  // path must decline and leave the exact gauge in charge.
  auto declined=affine_operator::prepare(series,physical,point,10,numeric_settings);
  check(declined.success() && !declined.recurrence->numeric_gauge && declined.row_pole_loss==original.row_pole_loss,"numeric row gauge must not replace a shear");
  std::cout<<"Exact endpoint epsilon-gauge cancellation preserves Laurent operator and reduces normalization loss\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
