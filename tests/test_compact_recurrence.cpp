// The compact recurrence core against the ball recurrence on the same charts.
#include "diffexp/rational_circuit.hpp"
#include <iostream>
using namespace diffexp;
namespace rc=rational_circuit;
using B=Jet::Ball;
static void check(bool value,const char* message){if(!value)throw std::runtime_error(message);}
static double log2_radius(const B& value) {
  mag_t r;mag_init(r);mag_hypot(r,arb_radref(acb_realref(value.raw())),arb_radref(acb_imagref(value.raw())));
  const double result=mag_is_zero(r)?-1e300:mag_get_d_log2_approx(r);mag_clear(r);return result;
}
static double log2_size(const B& value) {
  mag_t r;mag_init(r);acb_get_mag(r,value.raw());
  const double result=mag_is_zero(r)?-1e300:mag_get_d_log2_approx(r);mag_clear(r);return result;
}
// Both recurrences on one chart: enclosures must overlap everywhere, and the
// compact radii must stay within 2^slack of the ball radii or of 2^-(bits-40)
// relative to the coefficient.
static void compare(const std::vector<RationalLineEntry>& entries,unsigned dimension,unsigned epsilon_high,
    const B& center,const B& step,const Boundary& boundary,unsigned order,double slack,const char* name) {
  const auto compiled=rc::compile(entries,dimension,0,std::size_t(1),epsilon_high,order);
  const auto prepared=rc::prepare(compiled,center,order);
  rc::Statistics ball_stats,compact_stats;
  const auto reference=rc::expand_balls(prepared,boundary,&ball_stats);
  auto compact=rc::expand_compact(prepared,boundary,&compact_stats);
  check(compact.has_value(),"compact core declined an ordinary chart");
  check(compact->rows==reference.rows && compact->width==reference.width && compact->order==reference.order,"compact series shape");
  check(compact_stats.addmul_operations>0 && compact_stats.dot_products>0,"compact core reported no work");
  double step_log2=0;double tail_ball=0,tail_compact=0;
  const auto ball_values=rc::evaluate(reference,step,&tail_ball);
  const auto compact_values=rc::evaluate(*compact,step,&tail_compact);
  check(std::abs(tail_ball-tail_compact)<0.5 || (tail_ball< -1e200 && tail_compact< -1e200),"compact tail estimate differs");
  (void)step_log2;
  for(unsigned i=0;i<dimension;++i)for(unsigned k=0;k<boundary[0].size();++k) {
    check(acb_overlaps(ball_values[i][k].raw(),compact_values[i][k].raw()),"compact chart value excluded by the ball recurrence");
    check(log2_radius(compact_values[i][k])<=std::max(log2_radius(ball_values[i][k])+slack,log2_size(ball_values[i][k])-(B::precision()-40)),
      "compact chart value is far wider than the ball result");
  }
  rc::materialize(*compact);
  for(unsigned n=0;n<=order;++n)for(std::size_t r=0;r<reference.rows.size();++r)for(unsigned k=0;k<reference.width;++k) {
    const auto& a=reference.at(n,r,k);const auto& b=compact->at(n,r,k);
    if(!acb_overlaps(a.raw(),b.raw())){std::cerr<<name<<": order "<<n<<" row "<<r<<" epsilon "<<k<<"\n";check(false,"compact Taylor coefficient excluded by the ball recurrence");}
    check(a.is_zero()==b.is_zero() || !a.is_zero(),"compact core lost an exact structural zero");
    check(log2_radius(b)<=std::max(log2_radius(a)+slack,log2_size(a)-(B::precision()-40)),"compact Taylor coefficient is far wider than the ball one");
  }
}
// Both recurrences with defect radii (the midpoint transport of probe mode)
// on one chart. The radii are then the chart's own rounding, not enclosures
// of a common exact polynomial: the midpoints of the two cores may differ by
// that rounding, so they must agree within 2^slack times the larger radius,
// and the radii must be of the same order wherever they are above the
// rounding floor of the coefficient.
static double log2_distance(const B& a,const B& b) {
  B difference;acb_sub(difference.raw(),a.raw(),b.raw(),2*B::precision());
  acb_get_mid(difference.raw(),difference.raw());return log2_size(difference);
}
static double worst_defect_gap=-1e300,worst_defect_ratio=0;
static void compare_defect(const std::vector<RationalLineEntry>& entries,unsigned dimension,unsigned epsilon_high,
    const B& center,const B& step,Boundary boundary,unsigned order,double slack,const char* name) {
  // The transport folds every radius into its probes before such a chart.
  for(auto& row:boundary)for(auto& value:row)acb_get_mid(value.raw(),value.raw());
  const auto compiled=rc::compile(entries,dimension,0,std::size_t(1),epsilon_high,order);
  const auto prepared=rc::prepare(compiled,center,order);
  const auto reference=rc::expand_balls(prepared,boundary,nullptr,true);
  auto compact=rc::expand_compact(prepared,boundary,nullptr,true);
  check(compact.has_value(),"compact core declined an ordinary chart with defect radii");
  const auto enclosure=rc::expand_balls(prepared,boundary,nullptr,false);
  const auto judge=[&](const B& a,const B& b,const B& enclosed,const char* what,unsigned n,unsigned r,unsigned k) {
    const double floor=log2_size(a)-(B::precision()-40),ra=log2_radius(a),rb=log2_radius(b),widest=std::max(ra,rb);
    const double distance=log2_distance(a,b);
    const auto fail=[&](const char* why){std::cerr<<name<<": "<<what<<" order "<<n<<" row "<<r<<" epsilon "<<k<<": radii 2^"<<ra<<" (balls) 2^"<<rb
      <<" (compact), midpoints apart 2^"<<distance<<", size 2^"<<log2_size(a)<<"\n";check(false,why);};
    if(a.is_zero()!=b.is_zero())fail("the cores disagree on an exact zero with defect radii");
    if(distance>std::max(widest+slack,floor))fail("midpoints of the two cores differ by more than their defect radii allow");
    if(std::max(ra,rb)>floor && std::abs(ra-rb)>slack)fail("defect radii of the two cores are not of the same order");
    // Neither core may be farther from the enclosure than its radius allows.
    if(!acb_overlaps(a.raw(),enclosed.raw()) && log2_distance(a,enclosed)>std::max(ra+slack,floor))fail("ball defect result left the enclosure");
    if(!acb_overlaps(b.raw(),enclosed.raw()) && log2_distance(b,enclosed)>std::max(rb+slack,floor))fail("compact defect result left the enclosure");
    if(distance>-1e200 && widest>-1e200)worst_defect_gap=std::max(worst_defect_gap,distance-widest);
    if(ra>-1e200 && rb>-1e200 && std::max(ra,rb)>floor)worst_defect_ratio=std::max(worst_defect_ratio,std::abs(ra-rb));
  };
  double tail=0;
  const auto ball_values=rc::evaluate(reference,step,&tail),compact_values=rc::evaluate(*compact,step,&tail),
    enclosed_values=rc::evaluate(enclosure,step,&tail);
  for(unsigned i=0;i<dimension;++i)for(unsigned k=0;k<boundary[0].size();++k)
    judge(ball_values[i][k],compact_values[i][k],enclosed_values[i][k],"value",order,i,k);
  rc::materialize(*compact);
  for(unsigned n=0;n<=order;++n)for(std::size_t r=0;r<reference.rows.size();++r)for(unsigned k=0;k<reference.width;++k)
    judge(reference.at(n,r,k),compact->at(n,r,k),enclosure.at(n,r,k),"coefficient",n,unsigned(r),k);
}
int main(){try {
  ExactField field({"x","eps","I"});
  auto e=[&](const char* s){return Exact(field,s);};
  for(slong bits:{192,320,384,640}) {
    B::set_precision(bits);
    // Mixed epsilon blocks, quotient streams with epsilon in the denominator,
    // a constant source, complex entries.
    std::vector<RationalLineEntry> entries{
      {0,0,0,e("eps/(2-x+eps)")},{0,1,0,e("I/(1+x)")},{1,0,0,e("(3+x^2)/(5-2*x+x^2)")},{1,1,0,e("eps*(1-x)/(3+x)^2")},
      {2,0,0,e("(1+eps*x)/(7+x)")},{2,2,0,e("-2/(4-x)")},{0,3,0,e("1/(9+x^3)")},{1,3,0,e("x/(2-x+eps)")},{2,1,0,e("(x^5-3*x+1)/(11-x)^3")}};
    B a=B::from_strings("1.25","0.5"),b=B::from_strings("-0.75","0.125"),c=B::from_strings("0.3","-2");
    arb_add_error_2exp_si(acb_realref(a.raw()),-bits+30);arb_add_error_2exp_si(acb_imagref(b.raw()),-bits+20);
    Boundary boundary{{a,B(0),b,B(2)},{b,c,B(0),B(0)},{B(0),B(0),B(0),B(0)},{B(1),B(0),B(0),B(0)}};
    compare(entries,4,3,B::from_strings("0.125","0"),B::from_strings("0.03125","0"),boundary,60,12,"real centre");
    compare_defect(entries,4,3,B::from_strings("0.125","0"),B::from_strings("0.03125","0"),boundary,60,12,"real centre");
    compare(entries,4,3,B::from_strings("0.125","0.0625"),B::from_strings("0.046875","0"),boundary,80,12,"complex centre");
    compare_defect(entries,4,3,B::from_strings("0.125","0.0625"),B::from_strings("0.046875","0"),boundary,80,12,"complex centre");
    // A fine centre: coefficients of several limbs.
    compare(entries,4,3,B::from_strings("0.12548828125","-0.03173828125"),B::from_strings("0.01953125","0"),boundary,40,12,"fine centre");
    compare_defect(entries,4,3,B::from_strings("0.12548828125","-0.03173828125"),B::from_strings("0.01953125","0"),boundary,40,12,"fine centre");
    // Sparse rows: structural zeros must stay exact zeros.
    Boundary sparse{{B(0),B(0),B(0),B(0)},{B(0),B(0),B(0),B(0)},{B(1),b,B(0),B(0)},{B(0),B(0),B(0),B(0)}};
    compare(entries,4,3,B::from_strings("0.25","0"),B::from_strings("0.0625","0"),sparse,48,12,"sparse boundary");
    compare_defect(entries,4,3,B::from_strings("0.25","0"),B::from_strings("0.0625","0"),sparse,48,12,"sparse boundary");
  }
  std::cout<<"Compact recurrence agreed with the ball recurrence on every chart\n"
    <<"Defect radii: midpoints of the two cores at most 2^"<<worst_defect_gap<<" of the larger radius apart, radii within 2^"<<worst_defect_ratio<<"\n";
} catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
