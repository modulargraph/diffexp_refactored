#include "diffexp/transport.hpp"
#include <iostream>
using namespace diffexp;
using B=kernel::ComplexBall;
namespace j=boost::json;
static void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
static void close(const B& a,const B& b,const char* message) {
  auto difference=transport::magnitude(a-b);auto tolerance=B::from_strings("1e-60");
  check(arb_lt(acb_realref(difference.raw()),acb_realref(tolerance.raw())),message);
}
int main(){try {
  B::set_precision(384);
  auto request=j::parse(R"({"paths":{"t":"x"},"entries":[{"row":0,"column":1,"epsilon":1,"variable":"dlog","expression":"t+Sqrt[1+t^2]"}]})").as_object();
  auto c=transport::compile(request,2,2,384);auto plan=transport::finite_lag_plan(c);
  check(plan && plan->gauge.size()==2 && plan->gauge[0]!=plan->gauge[1],"root gauge was not selected");
  // Translation must preserve a polynomial at a nonzero complex chart center.
  std::vector<B> p{B(1),B(2),B(3),B(4)};B at=B::from_strings("1/4","1/8"),step=B::from_strings("1/32");
  auto q=transport::finite_lag_detail::shifted(p,at);B left,right;
  for(unsigned n=p.size();n-->0;){left=left*(at+step)+p[n];right=right*step+q[n];}
  check(acb_overlaps(left.raw(),right.raw()),"polynomial translation");
  at=B::from_strings("1/4");
  for(int sign:{-1,1}) {
    auto roots=transport::principal_roots_at(c,at);roots[0]=roots[0]*B(sign);
    Boundary initial{{B(0),B(0),B(0)},{B(1),B(0),B(0)}};
    auto result=transport::chart(c,transport::numerical_entries(c),initial,roots,at,step,80,false,110);
    check(result.finite_lag,"default did not use finite-lag recurrence");
    B a,b;acb_asinh(a.raw(),at.raw(),384);acb_asinh(b.raw(),(at+step).raw(),384);
    close(result.values[0][1],B(sign)*(b-a),"nonzero-center root sheet solution");
    auto series=transport::chart(c,transport::numerical_entries(c),initial,roots,at,step,80,false,110,true,10000000,true,false);
    close(result.values[0][1],series.values[0][1],"series comparison");
    auto saved=transport::chart(c,transport::numerical_entries(c),initial,roots,at,step,80,true,110);
    check(!saved.finite_lag && saved.saved.contains("coefficients"),"saved physical series fallback");
    // Nonzero incoming uncertainty is retained, including at zero midpoint.
    arb_add_error_2exp_si(acb_realref(initial[1][0].raw()),-100);
    auto uncertain=transport::chart(c,transport::numerical_entries(c),initial,roots,at,step,80,false,110);
    B rho(1);acb_mul_2exp_si(rho.raw(),rho.raw(),-101);
    for(int corner:{-1,0,1}) {
      auto expected=B(sign)*(b-a)*(B(1)+B(corner)*rho);
      check(acb_contains(uncertain.values[0][1].raw(),expected.raw()),"incoming uncertainty lost");
    }
  }
  auto conditioned_fallback=transport::chart(c,transport::numerical_entries(c),{{B(0),B(0),B(0)},{B(1),B(0),B(0)}},
      {B(1)},B(0),B::from_strings("99/100"),80,false,110);
  check(!conditioned_fallback.finite_lag,"denominator feedback guard was bypassed");
  // An odd-root self edge cannot be removed by a diagonal root gauge.
  request["entries"].as_array()[0].as_object()["column"]=0;
  auto unsupported=transport::compile(request,1,1,384);
  check(!transport::finite_lag_plan(unsupported),"inconsistent gauge accepted");
  auto fallback=transport::chart(unsupported,transport::numerical_entries(unsupported),{{B(1),B(0)}},{B(1)},B(0),step,64,false,110);
  check(!fallback.finite_lag,"unsupported connection failed to fall back");
  // A whole complex path checks repeated translations and branch continuation.
  auto run_request=j::parse(R"({"dimension":2,"epsilon_order":2,"taylor_order":80,"working_bits":384,"accuracy_goal":55,"paths":{"t":"x/2+I*x*(1-x)/4"},"entries":[{"row":0,"column":1,"epsilon":1,"variable":"dlog","expression":"t+Sqrt[1+t^2]"}],"boundary":[["0","0","0"],["1","0","0"]]})");
  run_request.as_object()["boundary_errors"]=j::array{j::array{"0","0","0"},j::array{"10^-65","0","0"}};
  auto result=transport::run(run_request).as_object();
  check(result.at("recurrence").as_object().at("finite_lag_charts").as_uint64()>0,"full run did not use default");
  auto number=result.at("values").as_array()[0].as_array()[1].as_object();
  B actual=B::from_strings(transport::string(number.at("real_midpoint")),transport::string(number.at("imaginary_midpoint"))),expected;
  acb_asinh(expected.raw(),B::from_strings("1/2").raw(),384);close(actual,expected,"full path analytic reference");
  std::cout<<"Finite-lag default, shifted charts, sheets, uncertainty and fallback passed\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
