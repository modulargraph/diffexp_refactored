#include "diffexp/rational_circuit.hpp"
#include "diffexp/polynomial_transport.hpp"
#include <iostream>
using namespace diffexp;
namespace rc=rational_circuit;
using B=Jet::Ball;
void check(bool value,const char *message) {
  if (!value) throw std::runtime_error(message);
}
void close(const B &a,const B &b) {
  auto delta=a-b;
  mag_t upper,bound;
  mag_init(upper); mag_init(bound);
  acb_get_mag(upper,delta.raw());
  mag_set_ui_2exp_si(bound,1,-170);
  bool good=mag_cmp(upper,bound)<0;
  mag_clear(upper); mag_clear(bound);
  check(good,"retained polynomial comparison or enclosure width failed");
}
void close(const Boundary &a,const Boundary &b) {
  check(a.size()==b.size(),"result dimension");
  for (unsigned i=0;i<a.size();++i) {
    check(a[i].size()==b[i].size(),"result width");
    for (unsigned k=0;k<a[i].size();++k) close(a[i][k],b[i][k]);
  }
}
template<class Exception,class Function> void rejects(Function f,const char *message) {
  bool rejected=false;
  try { f(); } catch (const Exception &) { rejected=true; }
  check(rejected,message);
}
int main() {
  try {
    B::set_precision(320);
    ExactField field({"x","eps","I","unused"});
    auto e=[&](const char *s) {return Exact(field,s);};
    const B center=B::from_strings("0.125","0.0625"), step=B::from_strings("0.03125","-0.015625");
    // One matrix plan serves independent RHS payloads; complex center, mixing,
    // augmented constant forcing, epsilon denominator and large consumer lag.
    std::vector<RationalLineEntry> entries{
      {0,0,0,e("eps/(2-x+eps)")}, {0,1,0,e("I/(1+x)")},
      {1,0,0,e("(1+x^7)/(2-x+eps)")}, {1,1,1,e("(1+eps)/(2+x)")},
      {0,2,0,e("3+x^2")}, {1,2,0,e("1/(1+x+eps)")}};
    auto compiled=rc::compile(entries,3,0,1,4,32);
    auto prepared=rc::prepare(compiled,center,32);
    check(prepared.data().auxiliary_cells<compiled.data().streams.size()*32*5,"ring storage is not bounded");
    bool lagfound=false,alias=false;
    for (unsigned id=0;id<compiled.data().streams.size();++id) {
      if (compiled.data().streams[id].lag==7) {
        lagfound=true; check(prepared.data().ring_lengths[id]==8,"consumer lag not retained");
      }
      if (!prepared.data().ring_lengths[id]) alias=true;
    }
    check(lagfound && alias,"missing consumer-lag or denominator-one alias case");
    std::vector<RationalLineEntry> expanded;
    for (const auto &entry:entries) {
      auto cs=feynman::scalar_functional_detail::epsilon_series(entry.coefficient,1,4-entry.epsilon);
      for (unsigned k=0;k<cs.size();++k)
        if (!cs[k].is_zero()) expanded.push_back({entry.row,entry.column,k+entry.epsilon,cs[k]});
    }
    for (const Boundary &boundary : {
         Boundary{{B(2),B(1),B(0),B(-1),B(3)},{B(3),B(0),B(1),B(0),B(0)},{B(1),B(0),B(0),B(0),B(0)}},
         Boundary{{B(-3),B(2),B(4),B(0),B(1)},{B(1),B(1),B(-1),B(2),B(0)},{B(0),B(1),B(0),B(0),B(0)}}}) {
      rc::Statistics stats;
      auto value=rc::chart(prepared,boundary,step,&stats);
      close(value,rational_chart(expanded,boundary,center,step,32));
      close(value,polynomial_transport::chart(polynomial_transport::compile(entries,3,0,1,4,32),boundary,center,step,32));
      check(stats.live_cells==prepared.data().live_cells,"live cell schedule differs from execution");
      check(stats.addmul_operations+stats.scalar_operations<=prepared.data().operation_upper_bound,"operation prediction is too small");
      check(stats.addmul_operations>0 && stats.scalar_operations>0,"work not counted");
    }
    // Smaller known windows are prefixes, including a plan prepared for that
    // narrower width. No coefficient beyond the supplied window is invented.
    Boundary wide{{B(2),B(1),B(3),B(4),B(5)},
                  {B(3),B(0),B(1),B(7),B(-3)},
                  {B(1),B(0),B(0),B(0),B(0)}};
    auto wide_value=rc::chart(prepared,wide,step);
    for (auto &row:wide) row.resize(2);
    auto narrow_value=rc::chart(rc::prepare(compiled,center,32,2),wide,step);
    for (auto &row:wide_value) row.resize(2);
    close(wide_value,narrow_value);
    rejects<std::invalid_argument>([&]{rc::prepare(compiled,center,32,6);},"unknown epsilon window admitted");
    // A homogeneous map column cannot reach the disconnected third component.
    auto sparse=rc::compile({{0,0,0,e("1/(1-x)")},{1,0,0,e("1+x^3")},
                              {2,2,0,e("1/(1+x)")}},3,0,1,0,12);
    auto sparse_prepared=rc::prepare(sparse,B(0),12);
    rc::Statistics sparse_stats;
    auto sparse_value=rc::chart(sparse_prepared,{{B(1)},{B(0)},{B(0)}},step,&sparse_stats);
    check(sparse_value[2][0].is_zero(),"unreachable row changed");
    check(sparse_stats.live_cells<sparse_prepared.data().live_cells,"reachability did not compact workspace");
    close(sparse_value,rational_chart({{0,0,0,e("1/(1-x)")},{1,0,0,e("1+x^3")},
                                      {2,2,0,e("1/(1+x)")}},{{B(1)},{B(0)},{B(0)}},B(0),step,12));
    rc::Statistics zero_stats;
    auto allzero=rc::chart(sparse_prepared,{{B(0)},{B(0)},{B(0)}},step,&zero_stats);
    check(zero_stats.addmul_operations==0 && zero_stats.scalar_operations==0 &&
          zero_stats.auxiliary_cells==0,"zero payload performed arithmetic");
    // Storage index zero denotes the certified physical lower bound -3 here.
    // y'=eps*y, initial y=eps^-3 gives exp(eps*h)*eps^-3.
    auto negative=rc::compile({{0,0,0,e("eps")}},1,0,1,4,10);
    auto neg=rc::chart(negative,{{B(1),B(0),B(0),B(0),B(0)}},B(0),step,10);
    B expected(1);
    for (unsigned k=0;k<5;++k) {close(neg[0][k],expected);expected=expected*step/B(k+1);}
    // Exact cancellation, unsigned offset cancels denominator epsilon valuation.
    auto normalized=rc::compile({{0,0,2,e("(1+x)/(eps*(1+x))")}},1,0,1,4,10);
    close(neg,rc::chart(normalized,{{B(1),B(0),B(0),B(0),B(0)}},B(0),step,10));
    // Gaussian coefficient reduction affects exact zero and denominator identity.
    auto gaussian=rc::compile({{0,0,0,e("(1+I^2)/(1+x)")},{0,0,0,e("I^3/(I*(1+x))")},
                               {1,0,0,e("-2/(2+2*x)")}},2,0,1,0,12);
    check(gaussian.data().entries.size()==2 && gaussian.data().denominators.size()==1 &&
          gaussian.data().streams.size()==1,"Gaussian normalization or shared stream identity failed");
    close(rc::chart(gaussian,{{B(1)},{B(2)}},B(0),step,12),
          rational_chart({{0,0,0,e("-1/(1+x)")},{1,0,0,e("-1/(1+x)")}},{{B(1)},{B(2)}},B(0),step,12));
    // Independent finite polynomial oracle: y'=(1+x^7)/(1-x), source=1.
    auto lag=rc::compile({{0,1,0,e("(1+x^7)/(1-x)")}},2,0,1,0,20);
    auto got=rc::chart(lag,{{B(0)},{B(1)}},B(0),step,20);
    B oracle(0),power=step;
    for (unsigned n=1;n<=20;++n) {oracle+=power*B(n>=8?2:1)/B(n);power*=step;}
    close(got[0][0],oracle);
    // epsilon denominator closed form at finite spatial order: forcing 1/(1-eps).
    auto epsq=rc::compile({{0,1,0,e("1/(1-eps)")}},2,0,1,4,3);
    auto epsvalue=rc::chart(epsq,{{B(0),B(0),B(0),B(0),B(0)},{B(1),B(0),B(0),B(0),B(0)}},B(0),step,3);
    for (auto &a:epsvalue[0]) close(a,step);
    auto zero=rc::compile({},2,0,std::nullopt,1,4);
    Boundary constants{{B(4),B(-2)},{B(0),B(7)}};
    close(rc::chart(zero,constants,B(0),step,4),constants);
    rejects<rc::UnsupportedDomain>([&]{rc::compile({{0,0,0,e("1/eps")}},1,0,1,2,4);},"negative epsilon valuation accepted");
    rejects<rc::UnsupportedDomain>([&]{rc::compile({{0,0,0,e("unused")}},1,0,1,2,4);},"unresolved parameter accepted");
    rejects<rc::UnsupportedDomain>([&]{rc::compile({{0,0,0,e("1/(1+I^2)")}},1,0,1,2,4);},"zero specialized denominator accepted");
    rejects<rc::InvalidPivot>([&]{rc::prepare(rc::compile({{0,0,0,e("1/(x+eps)")}},1,0,1,2,4),B(0),4);},"nonunit pivot accepted");
    B uncertain(0); arb_add_error_2exp_si(acb_realref(uncertain.raw()),0);
    rejects<rc::InvalidPivot>([&]{rc::prepare(rc::compile({{0,0,0,e("1/x")}},1,0,1,0,4),uncertain,4);},"uncertain pivot accepted");
    rejects<std::invalid_argument>([&]{rc::chart(prepared,{{B(0)}},step);},"wrong boundary dimension accepted");
    B::set_precision(256);
    rejects<std::invalid_argument>([&]{rc::chart(prepared,{{B(0)},{B(0)},{B(1)}},step);},"prepared precision mismatch accepted");
    B::set_precision(320);
    rc::Options options;options.max_cells=10;
    rejects<std::length_error>([&]{rc::prepare(rc::compile(entries,3,0,1,4,32,options),center,32);},"storage budget ignored");
    options=rc::Options{};options.max_operations=10;
    rejects<std::length_error>([&]{rc::prepare(rc::compile(entries,3,0,1,4,32,options),center,32);},"operation budget ignored");
    std::cout << "Rational circuit reference, epsilon, complex, ring, alias, independent RHS, budgets and domain checks passed\n";
  } catch (const std::exception &error) {std::cerr<<error.what()<<'\n';return 1;}
}
