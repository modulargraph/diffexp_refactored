#include "diffexp/diagonal_log_gauge.hpp"
#include <iostream>
using namespace diffexp;
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
int main(){try {
 ExactField field({"x","eps"});Exact x(field,"x"),e(field,"eps"),one(field,1);
 auto check=[&](const Exact& q,std::size_t count){auto r=diagonal_log_gauge(q,0);require(r.supported,r.rejection_reason.c_str());require(r.factors.size()==count,"factor count");auto reconstructed=q.constant(0);for(const auto& f:r.factors)reconstructed=reconstructed+q.constant(f.exponent)*f.polynomial.derivative(0)/f.polynomial;require(reconstructed==q,"exact reconstruction");return r;};
 check(x.constant(0),0);check(one/x+one/(x-one),2);
 check(x.constant(3)/(x-one)-x.constant(5)/(x+x.constant(2)),2);
 auto half=check((x+x)/(x*x+one)/x.constant(2)-one/x,2);
 bool has_half=false;for(const auto& f:half.factors)has_half|=f.exponent.str()=="1/2";require(has_half,"half integer factor");
 check(x.constant(6)/(x.constant(3)*x+x.constant(2)),1);
 require(!diagonal_log_gauge(one/(x*x),0).supported,"double pole accepted");
 require(!diagonal_log_gauge(x+one/x,0).supported,"polynomial part accepted");
 require(!diagonal_log_gauge(one/(x*x+one),0).supported,"nonconstant algebraic residues accepted");
 require(!diagonal_log_gauge((x+one)/(x*x+one),0).supported,"nonproportional residue accepted");
 require(!diagonal_log_gauge(e/x,0).supported,"other variable accepted");
 require(!diagonal_log_gauge(one/(x.constant(3)*x),0).supported,"third exponent accepted by default");
 DiagonalLogGaugeOptions all;all.half_integer_only=false;require(diagonal_log_gauge(one/(x.constant(3)*x),0,all).supported,"rational exponent option");
 DiagonalLogGaugeOptions small;small.max_degree=1;require(!diagonal_log_gauge(x/(x*x+one),0,small).supported,"degree budget ignored");
 small={};small.max_coefficient_bits=2;require(!diagonal_log_gauge(one/(x+x.constant(100)),0,small).supported,"coefficient budget ignored");
 small={};small.max_factors=1;require(!diagonal_log_gauge(one/x+one/(x-one),0,small).supported,"factor budget ignored");
 require(!diagonal_log_gauge(one/x,2).supported,"invalid variable accepted");
 std::cout<<"diagonal log gauge tests passed\n";return 0;
 }catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}}
