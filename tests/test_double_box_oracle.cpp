#include "diffexp/double_box_oracle.hpp"
#include <iostream>
using namespace diffexp;using B=kernel::ComplexBall;
int main(){try {
 B::set_precision(256);auto ref=oracle::double_box_planar_reference(256);
 // Independent decimal pins from Scripts/verify_double_box_planar.m (data only).
 const std::array<const char*,5> pins{"12","2.6260083723848567163723883918602272247005343047831","-77.808155752086180303874072147668656652688133264584","-201.73675884525049357623701386093884968969770043649","-247.00367822059256700943770828449562375927814454946"};
 for(unsigned i=0;i<5;++i){auto difference=ref.coefficients[i]-B::from_strings(pins[i]);B magnitude;acb_abs(acb_realref(magnitude.raw()),difference.raw(),256);if(!arb_lt(acb_realref(magnitude.raw()),acb_realref(B::from_strings("1e-45").raw())))throw std::runtime_error("Smirnov pin mismatch");}
 if(!acb_equal_si(ref.at(-4).raw(),12))throw std::runtime_error("doublebox leading sign normalization");
 B gamma,log3(3);arb_const_euler(acb_realref(gamma.raw()),256);acb_log(log3.raw(),log3.raw(),256);
 if(!acb_overlaps(ref.at(-3).raw(),(B(15)*log3-B(24)*gamma).raw()))throw std::runtime_error("exp(-2gammaeps) normalization");
 auto low=oracle::double_box_planar_reference(128);for(int k=-4;k<=0;++k)if(!acb_overlaps(low.at(k).raw(),ref.at(k).raw()))throw std::runtime_error("bounded series precision consistency");
 for(int k=-4;k<=0;++k){std::cout<<"doublebox eps^"<<k<<" ";acb_printn(ref.at(k).raw(),48,0);std::cout<<'\n';}
 std::cout<<"native independent Smirnov double-box oracle PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
