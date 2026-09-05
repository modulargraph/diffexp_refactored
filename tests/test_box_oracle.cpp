#include "diffexp/box_oracle.hpp"
#include "diffexp/families.hpp"
#include "diffexp/jet.hpp"
#include <flint/acb_poly.h>
#include <iostream>
using namespace diffexp;
void require(bool yes,const char* why){if(!yes)throw std::runtime_error(why);}
int main(){try {
  using B=kernel::ComplexBall;B::set_precision(256);auto reference=oracle::massless_box_reference(Rational(-1),Rational("-1/3"),224);
  require(acb_equal_si(reference.at(-2).raw(),12),"box leading coefficient normalization");
  // Verify the native fixture's exact F polynomial, propagator ordering and
  // massless external kinematics before identifying it with the paper's box.
  ExactField field({"a","b","c","d"});std::vector<Exact> parameters;for(unsigned i=0;i<4;++i)parameters.push_back(Exact(field).variable(i));
  auto family=feynman::example_family("box");auto geometry=feynman::symanzik(family.momenta,parameters);
  require(geometry.F==Exact(field,"a*c+b*d/3"),"native box invariant assignment does not match analytic reference");
  require(family.dimension_at_epsilon_zero==4 && family.physical_count==4,"native box dimension/arity");
  // Independent gamma-series normalization check directly expands the quoted
  // paper expression rather than using the closed coefficient formulas.
  acb_poly_t e,gp,gm,g2,rg,ex,bracket,temp;for(auto p:{e,gp,gm,g2,rg,ex,bracket,temp})acb_poly_init(p);
  acb_poly_set_coeff_si(e,0,1);acb_poly_set_coeff_si(e,1,1);acb_poly_gamma_series(gp,e,3,256);
  acb_poly_set_coeff_si(e,1,-1);acb_poly_gamma_series(gm,e,3,256);
  acb_poly_set_coeff_si(e,1,-2);acb_poly_gamma_series(g2,e,3,256);
  acb_poly_mullow(temp,gm,gm,3,256);acb_poly_mullow(rg,temp,gp,3,256);acb_poly_div_series(temp,rg,g2,3,256);acb_poly_set(rg,temp);
  B log3(3),pi;acb_log(log3.raw(),log3.raw(),256);arb_const_pi(acb_realref(pi.raw()),256);
  acb_poly_zero(e);acb_poly_set_coeff_acb(e,1,log3.raw());acb_poly_exp_series(ex,e,3,256);
  B two(2),finite;acb_poly_scalar_mul(bracket,ex,two.raw(),256);
  acb_poly_get_coeff_acb(finite.raw(),bracket,0);acb_add(finite.raw(),finite.raw(),two.raw(),256);acb_poly_set_coeff_acb(bracket,0,finite.raw());
  acb_mul(finite.raw(),log3.raw(),log3.raw(),256);acb_mul(pi.raw(),pi.raw(),pi.raw(),256);acb_add(finite.raw(),finite.raw(),pi.raw(),256);
  B coefficient;acb_poly_get_coeff_acb(coefficient.raw(),bracket,2);acb_sub(coefficient.raw(),coefficient.raw(),finite.raw(),256);acb_poly_set_coeff_acb(bracket,2,coefficient.raw());
  acb_poly_mullow(temp,rg,bracket,3,256);acb_poly_scalar_mul(temp,temp,B(3).raw(),256);
  for(unsigned i=0;i<3;++i){acb_poly_get_coeff_acb(coefficient.raw(),temp,i);require(arb_contains_zero(acb_imagref(coefficient.raw())),"gamma series imaginary enclosure");arb_zero(acb_imagref(coefficient.raw()));require(acb_overlaps(reference.coefficients[i].raw(),coefficient.raw()),"restored rGamma normalization mismatch");}
  for(auto p:{e,gp,gm,g2,rg,ex,bracket,temp})acb_poly_clear(p);
  auto swapped=oracle::massless_box_reference(Rational("-1/3"),Rational(-1),224);
  for(unsigned i=0;i<3;++i)require(acb_overlaps(reference.coefficients[i].raw(),swapped.coefficients[i].raw()),"s-t symmetry");
  auto scaled=oracle::massless_box_reference(Rational(-2),Rational("-2/3"),224);
  B log2(2);acb_log(log2.raw(),log2.raw(),256);
  require(acb_overlaps((scaled.at(-1)*B(4)).raw(),(reference.at(-1)-log2*reference.at(-2)).raw()),"dimensional scaling of box residue");
  require(acb_overlaps((scaled.at(0)*B(4)).raw(),(reference.at(0)-log2*reference.at(-1)+log2*log2*reference.at(-2)/B(2)).raw()),"dimensional scaling of box finite part");
  bool rejected=false;try{oracle::massless_box_reference(Rational(1),Rational(-1));}catch(const std::invalid_argument&){rejected=true;}require(rejected,"unsupported timelike continuation accepted");
  for(int k=-2;k<=0;++k){std::cout<<"box eps^"<<k<<" ";acb_printn(reference.at(k).raw(),45,0);std::cout<<"\n";}
  std::cout<<"independent analytic box oracle PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
