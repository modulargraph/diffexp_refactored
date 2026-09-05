// Bounded opt-in chart experiment. Production headers remain unchanged.
#include "dot_kernels.hpp"

#include <iomanip>
int main(int argc,char** argv) {try {
 require(argc==4,"usage: benchmark_charts METHOD DENSE REPEATS");
 const std::string method=argv[1];const unsigned input_kind=std::stoi(argv[2]);const bool dense=input_kind!=0;const unsigned repeats=std::stoul(argv[3]);
 require(repeats>=1 && repeats<=8,"repeat budget 1..8");
 constexpr unsigned order=80,high=5,masters=10,observables=4,d=41;
 B::set_precision(384);ExactField field({"x","eps","I"});Exact x(field,"x"),one(field,"1");
 std::vector<RationalLineEntry> sparse;
 for(unsigned block=0;block<observables;++block)for(unsigned row=0;row<masters;++row)for(unsigned edge=0;edge<3;++edge)
  sparse.push_back({block*masters+row,block*masters+(row+edge)%masters,edge,(edge%2?-one:one)/(one.constant(2+edge)+x)});
 for(unsigned block=0;block<observables;++block)sparse.push_back({block*masters+3,d-1,0,one/(one+x*x)});
 auto compiled=adjoint_detail::compile(sparse);
 const auto center=B::from_strings("1/8","1/16"),step=B::from_strings("1/32","1/64");
 Boundary initial(d,std::vector<B>(high+1,B(0)));initial.back()[0]=B(1);
 for(unsigned i=0;i<d-1;++i)for(unsigned k=0;k<=high;++k)if(dense || (k==0 && i%masters==0))
  initial[i][k]=B::from_strings(std::to_string(1+(i+3*k)%7)+"/16",std::to_string(static_cast<int>(i%3)-1)+"/32");
 if(input_kind==2)for(unsigned i=0;i<d-1;++i)for(auto& value:initial[i])arb_add_error_2exp_si(acb_realref(value.raw()),-190);
 auto run=[&] {
  if(method=="baseline")return adjoint_detail::scalar_chart(compiled,initial,center,step,order);
  if(method=="guarded_dot")return adjoint_detail::chart(compiled,initial,center,step,order);
  if(method=="raw_scalar")return experimental_chart(compiled,initial,center,step,order,0);
  if(method=="dot")return experimental_chart(compiled,initial,center,step,order,1);
  if(method=="precise_dot")return experimental_chart(compiled,initial,center,step,order,2);
  throw std::invalid_argument("unknown method");
 };
 std::cout<<std::setprecision(10);
 // First kernel invocation follows input preparation but no chart warm-up.
 auto first=time_charts(1,run);auto warm=time_charts(repeats,run);
 auto reference=adjoint_detail::scalar_chart(compiled,initial,center,step,order);auto result=run();
 auto candidate=adjoint_detail::dot_chart(compiled,initial,center,step,order);
 const bool fallback=adjoint_detail::needs_rational_cross_check(initial,candidate);
 if(method=="guarded_dot") {
   const auto& expected=fallback?reference:candidate;
   for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=high;++k)
     require(acb_equal(expected[i][k].raw(),result[i][k].raw()),"guard chose wrong candidate");
 }
 std::cerr<<"input_kind="<<input_kind<<" fallback="<<fallback<<"\n";
 unsigned equal=0,overlap=0,contained=0;slong ref_accuracy=1000000,new_accuracy=1000000;
 double max_radius_ratio=0;
 for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=high;++k) {
  auto a=reference[i][k].raw();auto b=result[i][k].raw();
  require(acb_is_finite(b),"nonfinite result");equal+=acb_equal(a,b);overlap+=acb_overlaps(a,b);contained+=acb_contains(a,b);
  double ar=std::max(mag_get_d(arb_radref(acb_realref(a))),mag_get_d(arb_radref(acb_imagref(a))));
  double br=std::max(mag_get_d(arb_radref(acb_realref(b))),mag_get_d(arb_radref(acb_imagref(b))));
  if(ar)max_radius_ratio=std::max(max_radius_ratio,br/ar);
  if(!acb_is_zero(a)){ref_accuracy=std::min(ref_accuracy,acb_rel_accuracy_bits(a));new_accuracy=std::min(new_accuracy,acb_rel_accuracy_bits(b));}
 }
 require(overlap==d*(high+1),"non-overlapping retained enclosure");
 if(method=="inplace_horner" || method=="raw_scalar")require(equal==d*(high+1),"scalar path is not bitwise equal");
 std::cout<<"{\"method\":\""<<method<<"\",\"dense\":"<<(dense?"true":"false")<<",\"repeats\":"<<repeats
 <<",\"first_wall_seconds\":"<<first.wall<<",\"first_cpu_seconds\":"<<first.cpu
 <<",\"warm_wall_seconds\":"<<warm.wall<<",\"warm_cpu_seconds\":"<<warm.cpu
 <<",\"equal\":"<<equal<<",\"overlap\":"<<overlap<<",\"contained\":"<<contained
 <<",\"cells\":246,\"max_radius_ratio\":"<<max_radius_ratio
 <<",\"reference_min_accuracy_bits\":"<<ref_accuracy<<",\"result_min_accuracy_bits\":"<<new_accuracy<<"}\n";
 return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
