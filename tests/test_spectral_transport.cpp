#include "diffexp/transport.hpp"
#include <iostream>
using namespace diffexp;
using B=kernel::ComplexBall;
namespace j=boost::json;
static void check(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
static void near(const B& actual,const B& expected,const char* tolerance,const char* message){auto difference=transport::magnitude(actual-expected);auto limit=B::from_strings(tolerance);check(arb_lt(acb_realref(difference.raw()),acb_realref(limit.raw())),message);}
static B enclosed(const transport::SpectralResult& result,unsigned row,unsigned epsilon){auto value=result.values[row][epsilon];const auto& error=result.errors[row][epsilon];arb_add_error(acb_realref(value.raw()),acb_realref(error.raw()));arb_add_error(acb_imagref(value.raw()),acb_realref(error.raw()));return value;}
int main(){try{
  B::set_precision(256);
  auto tiny_factor=transport::spectral_detail::tail_factor(1e-300,4,32);
  check(tiny_factor && *tiny_factor>0,"positive geometric tail factor underflowed");
  transport::SpectralOptions options;options.accuracy_goal=40;options.seconds_budget=20;
  transport::SpectralDiagnostics diagnostics;
  auto scalar=j::parse(R"JSON({"paths":{"t":"x"},"entries":[{"row":0,"column":0,"epsilon":1,"variable":"form","expression":"1"}]})JSON").as_object();
  auto exponential=transport::compile(scalar,1,6,256);
  Boundary exp_initial{{B(1),B(0),B(0),B(0),B(0),B(0),B(0)}};
  auto exp_result=transport::spectral_try(exponential,exp_initial,options,diagnostics);
  if(!exp_result)throw std::runtime_error(diagnostics.rejection_reason);
  check(diagnostics.nodes_tried.size()>=3 && diagnostics.selected_nodes>=16,"spectral refinement skipped decay evidence");
  B factorial(1);for(unsigned epsilon=0;epsilon<=6;++epsilon){if(epsilon)factorial=factorial*B(epsilon);near(exp_result->values[0][epsilon],B(1)/factorial,"1e-60","spectral exponential coefficient");check(acb_contains(enclosed(*exp_result,0,epsilon).raw(),(B(1)/factorial).raw()),"spectral exponential estimated enclosure");}
  B::set_precision(384);
  auto cancellation_request=j::parse(R"JSON({"paths":{"t":"x"},"entries":[{"row":0,"column":1,"epsilon":1,"variable":"form","expression":"10^60*(1+2*x+3*x^2-12*x^3)+10^(-45)"}]})JSON").as_object();
  auto cancellation=transport::compile(cancellation_request,2,1,384);
  auto cancellation_result=transport::spectral_try(cancellation,{{B(0),B(0)},{B(1),B(0)}},options,diagnostics);
  if(!cancellation_result)throw std::runtime_error(diagnostics.rejection_reason);
  near(cancellation_result->values[0][1],B::from_strings("1e-45"),"1e-40","large internal cancellation with tiny endpoint");
  check(acb_contains(enclosed(*cancellation_result,0,1).raw(),B::from_strings("1e-45").raw()),"cancellation uncertainty lost");
  B::set_precision(256);
  auto log_request=j::parse(R"JSON({"paths":{"t":"x"},"entries":[{"row":0,"column":1,"epsilon":1,"variable":"dlog","expression":"1+t","coefficient":"2"}]})JSON").as_object();
  auto logarithm=transport::compile(log_request,2,2,256);
  Boundary log_initial{{B(0),B(0),B(0)},{B(1),B(0),B(0)}};
  options.accuracy_goal=25;options.max_nodes=12;
  check(!transport::spectral_try(logarithm,log_initial,options,diagnostics),"insufficient spectral node budget accepted");
  check(diagnostics.rejection_reason.find("node budget")!=std::string::npos,"missing node-budget rejection reason");
  options.max_nodes=128;options.matrix_bits=288;
  auto log_result=transport::spectral_try(logarithm,log_initial,options,diagnostics);
  if(!log_result)throw std::runtime_error(diagnostics.rejection_reason);
  B log2;acb_log(log2.raw(),B(2).raw(),256);
  near(log_result->values[0][1],B(2)*log2,"1e-25","spectral triangular logarithm");
  check(acb_contains(enclosed(*log_result,0,1).raw(),(B(2)*log2).raw()),"logarithm estimated tail enclosure");
  check(diagnostics.nodes_tried.size()>3,"nonpolynomial spectral tail was not refined");
  check(B::precision()==256,"spectral matrix precision escaped its scope");
  // Input boxes are carried independently of midpoint convergence. The source
  // interval here is much larger than the finite arithmetic rounding floor.
  arb_add_error_2exp_si(acb_realref(log_initial[1][0].raw()),-150);
  auto uncertain=transport::spectral_try(logarithm,log_initial,options,diagnostics);
  if(!uncertain)throw std::runtime_error(diagnostics.rejection_reason);
  B radius(1);acb_mul_2exp_si(radius.raw(),radius.raw(),-151);
  for(int sign:{-1,0,1}){auto expected=B(2)*log2*(B(1)+B(sign)*radius);check(acb_contains(enclosed(*uncertain,0,1).raw(),expected.raw()),"spectral input uncertainty lost");}
  check(!transport::arithmetic_error(uncertain->values[0][1]).is_zero(),"spectral endpoint box was replaced by its midpoint");
  auto contour_request=j::parse(R"JSON({"paths":{"t":"x/2+I*x*(1-x)/4"},"entries":[{"row":0,"column":1,"epsilon":1,"variable":"dlog","expression":"t+Sqrt[1+t^2]"}]})JSON").as_object();
  auto contour=transport::compile(contour_request,2,1,256);
  Boundary root_initial{{B(0),B(0)},{B(1),B(0)}};
  auto root_result=transport::spectral_try(contour,root_initial,options,diagnostics);
  if(!root_result)throw std::runtime_error(diagnostics.rejection_reason);
  B asinh;acb_asinh(asinh.raw(),B::from_strings("1/2").raw(),256);
  near(root_result->values[0][1],asinh,"1e-25","spectral algebraic complex contour");
  check(acb_contains(enclosed(*root_result,0,1).raw(),asinh.raw()),"complex contour estimated tail enclosure");
  // The radicand winds once around zero. Continuing its derivative integrates
  // to -2, whereas resampling principal roots would switch sheets midway.
  auto winding_request=j::parse(R"JSON({"paths":{"t":"x"},"entries":[{"row":0,"column":1,"epsilon":1,"variable":"form","expression":"(-8+16*x+I*(8-48*x+48*x^2))/(2*Sqrt[1-8*x*(1-x)+I*8*x*(1-x)*(1-2*x)])"}]})JSON").as_object();
  auto winding=transport::compile(winding_request,2,1,256);options.accuracy_goal=8;
  auto winding_result=transport::spectral_try(winding,root_initial,options,diagnostics);
  if(!winding_result)throw std::runtime_error(diagnostics.rejection_reason);
  near(winding_result->values[0][1],B(-2),"1e-8","spectral root winding continuation");
  check(acb_contains(enclosed(*winding_result,0,1).raw(),B(-2).raw()),"winding estimated tail enclosure");
  near(winding_result->endpoint_roots[0],B(-1),"1e-60","continued endpoint root was not retained");
  scalar["entries"].as_array()[0].as_object()["epsilon"]=0;
  auto noncanonical=transport::compile(scalar,1,1,256);
  check(!transport::spectral_try(noncanonical,{{B(1),B(0)}},options,diagnostics),"noncanonical spectral system accepted");
  auto pole_request=j::parse(R"JSON({"paths":{"t":"x"},"entries":[{"row":0,"column":0,"epsilon":1,"variable":"dlog","expression":"1-2*t"}]})JSON").as_object();
  auto pole=transport::compile(pole_request,1,1,256);
  check(!transport::spectral_try(pole,{{B(1),B(0)}},options,diagnostics),"spectral path singularity accepted");
  pole_request["entries"].as_array()[0].as_object()["expression"]="1-2*t+I/1000000";
  auto narrow=transport::compile(pole_request,1,1,256);
  check(!transport::spectral_try(narrow,{{B(1),B(0)}},options,diagnostics) && diagnostics.rejection_reason.find("nearby singularities")!=std::string::npos,"unresolved narrow feature accepted by spectral sampling");
  options.work_budget=1;
  check(!transport::spectral_try(exponential,exp_initial,options,diagnostics) && diagnostics.nodes_tried.empty(),"spectral work budget ignored");
  options.work_budget=100000000;options.seconds_budget=0;
  check(!transport::spectral_try(exponential,exp_initial,options,diagnostics),"invalid spectral time budget accepted");
  options.seconds_budget=20;options.accuracy_goal=40;options.matrix_bits=0;B::set_precision(64);
  check(!transport::spectral_try(exponential,exp_initial,options,diagnostics),"insufficient spectral precision accepted");
  B::set_precision(256);
  // Exercise selection through the public request interface, including the
  // native pullback of ordinary dt input and forced non-spectral choices.
  auto integration=j::parse(R"JSON({"dimension":1,"epsilon_order":3,"taylor_order":24,"working_bits":256,"accuracy_goal":20,"recurrence":"spectral","paths":{"t":"x/2"},"entries":[{"row":0,"column":0,"epsilon":1,"variable":"t","expression":"1"}],"boundary":[["1","0","0","0"]]})JSON").as_object();
  const auto endpoint=[](const j::object& output,unsigned row,unsigned epsilon){const auto& value=output.at("values").as_array()[row].as_array()[epsilon].as_object();return B::from_strings(transport::string(value.at("real_midpoint")),transport::string(value.at("imaginary_midpoint")));};
  auto explicit_result=transport::run(integration).as_object();
  check(explicit_result.at("spectral").as_object().at("accepted").as_bool(),"explicit spectral ordinary entry was not selected");
  near(endpoint(explicit_result,0,1),B::from_strings("1/2"),"1e-20","ordinary differential pullback in spectral integration");
  near(endpoint(explicit_result,0,3),B::from_strings("1/48"),"1e-20","ordinary differential epsilon hierarchy in spectral integration");
  check(!explicit_result.at("omitted_tails_certified").as_bool(),"spectral tails were mislabeled certified");
  auto named_form=integration;
  named_form["paths"]=j::object{{"form","2*x"}};
  named_form["entries"].as_array()[0].as_object()["variable"]="form";
  for(const auto* mode:{"spectral","series"}) {
    named_form["recurrence"]=mode;
    auto named_result=transport::run(named_form).as_object();
    near(endpoint(named_result,0,1),B(2),"1e-20","kinematic form lost its path derivative");
  }
  for(const auto* mode:{"taylor","series"}){
    integration["recurrence"]=mode;auto forced=transport::run(integration).as_object();
    check(!forced.at("spectral").as_object().at("attempted").as_bool(),"forced local recurrence attempted spectral transport");
    check(forced.at("recurrence").as_object().at("requested").as_string()==mode,"forced recurrence label changed");
    near(endpoint(forced,0,3),B::from_strings("1/48"),"1e-20","forced local recurrence endpoint");
  }
  integration["recurrence"]="auto";integration["dimension"]=64;
  j::array padded_boundary;for(unsigned row=0;row<64;++row)padded_boundary.push_back(j::array{row==0?"1":"0","0","0","0"});
  integration["boundary"]=std::move(padded_boundary);
  auto automatic=transport::run(integration).as_object();
  check(automatic.at("spectral").as_object().at("accepted").as_bool(),"eligible large auto request did not select spectral transport");
  near(endpoint(automatic,0,3),B::from_strings("1/48"),"1e-20","automatic spectral padded endpoint");
  integration["entries"].as_array()[0].as_object()["epsilon"]=0;
  auto local_fallback=transport::run(integration).as_object();
  check(!local_fallback.at("spectral").as_object().at("accepted").as_bool(),"noncanonical auto request accepted spectral transport");
  B exp_half;acb_exp(exp_half.raw(),B::from_strings("1/2").raw(),256);
  near(endpoint(local_fallback,0,0),exp_half,"1e-20","noncanonical automatic fallback endpoint");
  std::cout<<"Adaptive spectral precision, tails, budgets, boundaries and continued roots passed\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
