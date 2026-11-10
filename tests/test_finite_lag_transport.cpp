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
  // A supplied one-form uses its derivative directly, with a separate exact weight.
  auto forms=j::parse(R"({"paths":{"t":"x"},"entries":[{"row":0,"column":1,"epsilon":1,"variable":"form","expression":"1/Sqrt[1+x^2]","coefficient":"2"}]})").as_object();
  auto form_system=transport::compile(forms,2,2,384);
  check(form_system.canonical && transport::finite_lag_plan(form_system),"supplied form lost sparse finite-lag support");
  Boundary form_initial{{B(0),B(0),B(0)},{B(1),B(0),B(0)}};
  B form_center=B::from_strings("1/4"),form_step=B::from_strings("1/32");
  auto form_roots=transport::principal_roots_at(form_system,form_center);
  auto form_result=transport::chart(form_system,transport::numerical_entries(form_system),form_initial,form_roots,form_center,form_step,80,false,110);
  auto form_series=transport::chart(form_system,transport::numerical_entries(form_system),form_initial,form_roots,form_center,form_step,80,false,110,true,10000000,true,false);
  B form_a,form_b;acb_asinh(form_a.raw(),form_center.raw(),384);acb_asinh(form_b.raw(),(form_center+form_step).raw(),384);
  check(form_result.finite_lag,"supplied form did not select finite-lag recurrence");
  close(form_result.values[0][1],B(2)*(form_b-form_a),"supplied form analytic result");
  close(form_series.values[0][1],form_result.values[0][1],"supplied form series fallback");
  // Identical syntax has distinct meaning as a letter and as a one-form.
  forms["entries"].as_array().push_back(j::parse(R"({"row":0,"column":1,"epsilon":1,"variable":"dlog","expression":"1/Sqrt[1+x^2]"})"));
  auto mixed=transport::compile(forms,2,2,384);
  check(mixed.letters.size()==2 && mixed.letter_forms[0] && !mixed.letter_forms[1],"mixed letter/form cache collision");
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
  // A truncation refinement shortens the step of the SAME retained polynomial.
  // It is independent of the recurrence: the finite-lag chart must neither be
  // blamed for it nor recomputed, and must agree with a fresh shorter chart.
  auto dilog=j::parse(R"({"dimension":3,"epsilon_order":2,"taylor_order":50,"working_bits":256,"accuracy_goal":30,"paths":{"t":"3*x/4"},"entries":[{"row":1,"column":0,"epsilon":1,"variable":"dlog","expression":"1-t","coefficient":"-1"},{"row":2,"column":1,"epsilon":1,"variable":"dlog","expression":"t","coefficient":"1"}],"boundary":[["1","0","0"],["0","0","0"],["0","0","0"]]})");
  auto refined=transport::run(dilog).as_object();const auto& recurrence=refined.at("recurrence").as_object();
  check(recurrence.at("finite_lag_charts").as_uint64()>0 && recurrence.at("finite_lag_rejections").as_uint64()==0 &&
      recurrence.at("series_charts").as_uint64()==0,"truncation refinement was charged to the finite-lag recurrence");
  B::set_precision(256);
  auto cell=refined.at("values").as_array()[2].as_array()[2].as_object();
  B li2=B::from_strings(transport::string(cell.at("real_midpoint"))),li2_expected;
  acb_polylog_si(li2_expected.raw(),2,B::from_strings("3/4").raw(),256);
  check(arb_lt(acb_realref(transport::magnitude(li2-li2_expected).raw()),acb_realref(B::from_strings("1e-30").raw())),"refined dilogarithm reference");
  B::set_precision(384);
  {
    const B wide=B::from_strings("1/4"),narrow=B::from_strings("1/9");transport::Chart first;
    auto roots=transport::principal_roots_at(c,at);Boundary initial{{B(0),B(0),B(0)},{B(1),B(0),B(0)}};
    arb_add_error_2exp_si(acb_realref(initial[1][0].raw()),-100);
    auto expansion=transport::expand(c,transport::numerical_entries(c),initial,roots,at,wide,80,false,110,first);
    check(expansion.finite_lag() && first.finite_lag,"expansion did not retain finite-lag coefficients");
    auto shorter=transport::evaluate_expansion(c,expansion,narrow,false,110);
    auto fresh=transport::chart(c,transport::numerical_entries(c),initial,roots,at,narrow,80,false,110);
    check(shorter && shorter->finite_lag && fresh.finite_lag,"finite-lag re-evaluation was not used");
    transport::Chart series_first;
    auto series=transport::expand(c,transport::numerical_entries(c),initial,roots,at,wide,80,false,110,series_first,true,10000000,true,false);
    auto series_shorter=transport::evaluate_expansion(c,series,narrow,false,110);
    auto series_fresh=transport::chart(c,transport::numerical_entries(c),initial,roots,at,narrow,80,false,110,true,10000000,true,false);
    check(!series.finite_lag() && series_shorter,"series expansion was not retained");
    for(unsigned i=0;i<2;++i)for(unsigned e=0;e<3;++e) {
      check(acb_equal(shorter->values[i][e].raw(),fresh.values[i][e].raw()) &&
          acb_equal(shorter->truncation_errors[i][e].raw(),fresh.truncation_errors[i][e].raw()),"finite-lag re-evaluation changed retained balls");
      check(acb_equal(series_shorter->values[i][e].raw(),series_fresh.values[i][e].raw()) &&
          acb_equal(series_shorter->errors[i][e].raw(),series_fresh.errors[i][e].raw()),"series re-evaluation changed retained balls");
    }
  }
  std::cout<<"Finite-lag default, shifted charts, sheets, uncertainty and fallback passed\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
