#include "diffexp/box_triangle_oracle.hpp"
#include <iostream>
#include "diffexp/families.hpp"
using B=diffexp::kernel::ComplexBall;
int main(){try{B::set_precision(256);
 auto full=diffexp::feynman::example_family("double_box_planar").momenta;
 auto diagonal=full;diagonal.lines.clear();for(unsigned i:{1u,2u,3u,4u,6u})diagonal.lines.push_back(full.lines[i]);
 diffexp::ExactField field({"a","b","c","d","e"});std::vector<diffexp::Exact> p;for(unsigned i=0;i<5;++i)p.push_back(diffexp::Exact(field).variable(i));
 auto uf=diffexp::feynman::symanzik(diagonal,p);
 if(!(uf.U==diffexp::Exact(field,"(a+b)*(d+e)+c*(a+b+d+e)")) || !(uf.F==diffexp::Exact(field,"c*(b*d+a*e/3)")))throw std::runtime_error("native diagonal MB geometry mismatch");
 auto r=diffexp::oracle::box_triangle_reference();for(int k=-4;k<=0;++k){std::cout<<"boxtriangle eps^"<<k<<" ";acb_printn(r.at(k).raw(),40,0);std::cout<<" resolution difference ";acb_printn(r.resolution_differences[k+4].raw(),8,0);std::cout<<'\n';B a;acb_abs(acb_realref(a.raw()),r.resolution_differences[k+4].raw(),256);if(!arb_lt(acb_realref(a.raw()),acb_realref(B::from_strings("1e-30").raw())))throw std::runtime_error("MB/Cauchy reference resolutions disagree");}B leading;auto leading_error=r.at(-4)-B::from_strings("9/4");acb_abs(acb_realref(leading.raw()),leading_error.raw(),256);if(!arb_lt(acb_realref(leading.raw()),acb_realref(B::from_strings("1e-30").raw())))throw std::runtime_error("boxtriangle leading pole normalization");std::cout<<"boxtriangle independent numerical reference PASS (not certified)\n";}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
