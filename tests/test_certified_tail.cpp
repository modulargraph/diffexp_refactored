#include "diffexp/certified_tail.hpp"
#include "diffexp/rational_transport.hpp"
#include <iostream>
using namespace diffexp;
using B=kernel::ComplexBall;
void require(bool condition,const char* message) {if(!condition)throw std::runtime_error(message);}
template<class F> void rejects(F f,const char* message) {
  bool threw=false;try{f();}catch(const std::exception&){threw=true;}require(threw,message);
}
NativeAnalyticEntry entry(unsigned row,unsigned column,unsigned epsilon,const char* expression) {
  return {row,column,epsilon,data::Reader(expression).read()};
}
NativeTailBoundary enclose(NativeTailBoundary polynomial,const NativeTaylorTail& tail) {
  for(auto& row:polynomial)for(unsigned k=0;k<row.size();++k)tail.absolute[k].add_error_to(row[k]);
  return polynomial;
}
int main() {try {
  B::set_precision(256);
  const B zero(0),step=B::from_strings("1/8");
  ExactField field({"x"});
  // y=(1-x)^(-epsilon): exact logarithmic epsilon coefficients are independent
  // of both the Taylor recurrence and the Cauchy bound under test.
  auto disk=NativeAnalyticDisk::certify({entry(0,0,1,"1/(1-x)")},1,4,zero,"1/2");
  NativeTailBoundary initial(1,std::vector<B>(5));initial[0][0]=B(1);
  std::vector<RationalLineEntry> matrix{{0,0,1,Exact(field,"1/(1-x)")}};
  auto polynomial=rational_chart(matrix,initial,zero,step,24);
  auto tail=certify_native_taylor_tail(disk,initial,step,24);
  auto full=enclose(polynomial,tail);
  B log;auto argument=B(1)-step;acb_log(log.raw(),argument.raw(),256);log=-log;
  B reference(1);
  for(unsigned k=0;k<=4;++k) {
    if(k)reference=reference*log/B(k);
    require(acb_contains(full[0][k].raw(),reference.raw()),"logarithmic coefficient not enclosed");
    require(tail.absolute[k].approximate_upper()<1e-12,"canonical bound loses useful precision");
  }
  require(tail.absolute[0].is_zero(),"constant epsilon-zero coefficient has no tail");
  require(native_enclosure_meets_tolerance(full,"1/1000000000000"),"certified logarithm accuracy");
  require(!native_enclosure_meets_tolerance(full,"1/1000000000000000000000000000000"),"false precision must be rejected");
  require(!acb_contains(polynomial[0][1].raw(),log.raw()),"finite polynomial alone must not certify log");

  // Ordinary scalar exponential, with a nonzero epsilon-zero connection.
  auto ordinary=NativeAnalyticDisk::certify({entry(0,0,0,"2")},1,0,zero,"1");
  NativeTailBoundary one{{B(1)}};B h=B::from_strings("1/4");
  auto exp_tail=certify_native_taylor_tail(ordinary,one,h,24);
  auto exp_full=enclose(rational_chart(std::vector<RationalLineEntry>{{0,0,0,Exact(field,2)}},one,zero,h,24),exp_tail);
  B exponential;auto exponent=B(2)*h;acb_exp(exponential.raw(),exponent.raw(),256);
  require(acb_contains(exp_full[0][0].raw(),exponential.raw()),"ordinary exponential not enclosed");
  require(exp_tail.absolute[0].approximate_upper()<1e-13,"ordinary bound unreasonably loose");

  // Nilpotent coupled canonical system: y2=1, y1=2 eps L,
  // y0=eps^2 L^2, L=-log(1-x). This exercises both component and epsilon flow.
  auto coupled=NativeAnalyticDisk::certify({entry(0,1,1,"1/(1-x)"),entry(1,2,1,"2/(1-x)")},3,2,zero,"1/2");
  NativeTailBoundary coupled_initial(3,std::vector<B>(3));coupled_initial[2][0]=B(1);
  std::vector<RationalLineEntry> coupled_matrix{{0,1,1,Exact(field,"1/(1-x)")},{1,2,1,Exact(field,"2/(1-x)")}};
  auto coupled_tail=certify_native_taylor_tail(coupled,coupled_initial,step,24);
  auto coupled_full=enclose(rational_chart(coupled_matrix,coupled_initial,zero,step,24),coupled_tail);
  auto square=log*log;auto twice=B(2)*log;
  require(acb_contains(coupled_full[0][2].raw(),square.raw()),"coupled nilpotent square not enclosed");
  require(acb_contains(coupled_full[1][1].raw(),twice.raw()),"coupled nilpotent log not enclosed");

  // A1=1000 shows why finite epsilon depth must not be charged exp(1000R).
  auto large=NativeAnalyticDisk::certify({entry(0,0,1,"1000")},1,4,zero,"1");
  auto finite=certify_native_taylor_tail(large,initial,B::from_strings("1/16"),64);
  require(finite.circle_upper[4].approximate_upper()<5e10,"epsilon majorant improperly exponentiates full matrix norm");
  require(finite.absolute[4].approximate_upper()<1e-60,"finite epsilon depth should give useful tail");
  // Epsilon-quadratic connection tests the general finite convolution.
  auto mixed=NativeAnalyticDisk::certify({entry(0,0,1,"1"),entry(0,0,2,"2")},1,4,zero,"1/2");
  auto mixed_tail=certify_native_taylor_tail(mixed,initial,step,3);
  std::vector<RationalLineEntry> mixed_matrix{{0,0,1,Exact(field,1)},{0,0,2,Exact(field,2)}};
  auto mixed_full=enclose(rational_chart(mixed_matrix,initial,zero,step,3),mixed_tail);
  auto exact_fourth=B(2)*step*step+step*step*step+step*step*step*step/B(24);
  require(acb_contains(mixed_full[0][4].raw(),exact_fourth.raw()),"epsilon-quadratic convolution not enclosed");

  // A large one-way integral forcing must enter linearly, not in exp(M0 R).
  auto accumulated_disk=NativeAnalyticDisk::certify({entry(0,0,0,"1"),entry(0,0,1,"1"),
    entry(1,0,0,"1000000000000"),entry(1,0,2,"1000000000000")},2,4,zero,"1/2");
  require(accumulated_disk.scalar_accumulator(),"scalar integral triangular shape not recognized");
  NativeTailBoundary accumulated_initial(2,std::vector<B>(5,B(0)));accumulated_initial[0][0]=B(1);
  auto accumulated_tail=certify_native_taylor_tail(accumulated_disk,accumulated_initial,step,40);
  require(accumulated_tail.circle_upper[0].approximate_upper()<1e13 && accumulated_tail.absolute[0].approximate_upper()<1e-10,"large integral forcing must not cause exponential norm inflation");
  auto accumulated_full=enclose(rational_chart(std::vector<RationalLineEntry>{{0,0,0,Exact(field,"1")},{0,0,1,Exact(field,"1")},
    {1,0,0,Exact(field,"1000000000000")},{1,0,2,Exact(field,"1000000000000")}},accumulated_initial,zero,step,40),accumulated_tail);
  Jet epsilon(0,5,256);epsilon.set(1,B(1));auto unit=epsilon.constant(1),hjet=epsilon.constant(0);hjet.set(0,step);
  auto source_reference=((unit+epsilon)*hjet).exp();
  auto integral_reference=epsilon.decimal("1000000000000")*(unit+epsilon*epsilon)*(source_reference-unit)/(unit+epsilon);
  for(unsigned k=0;k<=4;++k)require(acb_contains(accumulated_full[0][k].raw(),source_reference.at(k).raw()) && acb_contains(accumulated_full[1][k].raw(),integral_reference.at(k).raw()),"triangular epsilon-polynomial integral must enclose analytic exponential solution");
  auto large_constant=accumulated_initial;large_constant[1][0]=B::from_strings("1000000000000000000000000000000");
  NativeTailMagnitude::decimal("0.01").add_error_to(large_constant[1][0]);
  auto constant_tail=certify_native_taylor_tail(accumulated_disk,large_constant,step,40);
  require(constant_tail.absolute[0].approximate_upper()<1e-10 && constant_tail.circle_upper[0].approximate_upper()>1e29,"constant accumulator must remain in full circle bound but cannot inflate omitted Taylor terms");
  require(!native_enclosure_meets_tolerance(large_constant,"1/1000"),"preexisting accumulator uncertainty must not be relabeled as a small Taylor tail");
  auto feedback=NativeAnalyticDisk::certify({entry(0,1,0,"1"),entry(1,0,0,"1000000000000")},2,0,zero,"1/2");
  require(!feedback.scalar_accumulator(),"feedback system must not receive a one-way forcing majorant");

  // Input uncertainty is propagated as an enclosure, separately from the tail.
  NativeTailBoundary uncertain{{B(1)}};NativeTailMagnitude::decimal("0.01").add_error_to(uncertain[0][0]);
  auto uncertain_tail=certify_native_taylor_tail(ordinary,uncertain,h,24);
  auto uncertain_full=enclose(rational_chart(std::vector<RationalLineEntry>{{0,0,0,Exact(field,2)}},uncertain,zero,h,24),uncertain_tail);
  auto low=B::from_strings("99/100")*exponential,high=B::from_strings("101/100")*exponential;
  require(acb_contains(uncertain_full[0][0].raw(),low.raw()) && acb_contains(uncertain_full[0][0].raw(),high.raw()),"input uncertainty lost");
  require(!native_enclosure_meets_tolerance(uncertain_full,"1/1000"),"tail size must not hide boundary uncertainty");

  // Reject a pole inside the witness disk even when the evaluated step is safe.
  rejects([&]{NativeAnalyticDisk::certify({entry(0,0,1,"1/(1-3*x)")},1,1,zero,"1/2");},"hidden outer-disk pole accepted");
  rejects([&]{certify_native_taylor_tail(disk,initial,B::from_strings("1/2"),24);},"disk boundary accepted");
  rejects([&]{NativeAnalyticDisk::certify({},1,1,zero,"-1");},"negative radius accepted");
  rejects([&]{certify_native_taylor_tail(disk,{{B(1)}},step,24);},"incomplete epsilon boundary accepted");
  rejects([&]{NativeAnalyticDisk::certify({entry(0,0,0,"Sqrt[x]")},1,0,zero,"1/2");},"branch point accepted");
  rejects([&]{NativeAnalyticDisk::certify({entry(0,0,0,"Sqrt[-2+x]")},1,0,zero,"1/2");},"principal cut accepted");
  auto algebraic_entry=entry(0,0,1,"Sqrt[2+x]");algebraic_entry.dlog=true;
  auto algebraic=NativeAnalyticDisk::certify({algebraic_entry},1,1,zero,"1/2");
  require(algebraic.matrix_norms()[1].is_finite(),"valid algebraic sheet not certified");
  NativeTailBoundary algebraic_initial{{B(1),B(0)}};
  auto algebraic_tail=certify_native_taylor_tail(algebraic,algebraic_initial,step,24);
  auto algebraic_full=enclose(rational_chart(std::vector<RationalLineEntry>{{0,0,1,Exact(field,"1/(2*(2+x))")}},algebraic_initial,zero,step,24),algebraic_tail);
  B algebraic_reference;auto algebraic_ratio=(B(2)+step)/B(2);
  acb_log(algebraic_reference.raw(),algebraic_ratio.raw(),256);algebraic_reference=algebraic_reference/B(2);
  require(acb_contains(algebraic_full[0][1].raw(),algebraic_reference.raw()),"algebraic dlog not enclosed");
  auto vanishing_letter=entry(0,0,1,"x");vanishing_letter.dlog=true;
  rejects([&]{NativeAnalyticDisk::certify({vanishing_letter},1,1,zero,"1/2");},"vanishing dlog letter accepted");
  auto small_tail=certify_native_taylor_tail(ordinary,one,h,64);
  require(small_tail.absolute[0].approximate_upper()<1e-35,"precision test requires negligible Taylor error");
  auto high_precision=enclose(rational_chart(std::vector<RationalLineEntry>{{0,0,0,Exact(field,2)}},one,zero,h,64),small_tail);
  require(native_enclosure_meets_tolerance(high_precision,"1/1000000000000000000000000000000"),"high arithmetic precision should meet tolerance");
  B::set_precision(64);
  auto low_precision=enclose(rational_chart(std::vector<RationalLineEntry>{{0,0,0,Exact(field,2)}},one,zero,h,64),small_tail);
  require(!native_enclosure_meets_tolerance(low_precision,"1/1000000000000000000000000000000"),"finite arithmetic false precision accepted");
  std::cout<<"certified native Taylor tails: logarithm, exponential, nilpotent/epsilon-polynomial coupling, domain and precision checks passed\n";
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;} }
