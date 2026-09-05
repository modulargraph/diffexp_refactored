#include "../tools/stage_probe_io.hpp"
#include <iostream>
using namespace diffexp;using namespace diffexp::stage_probe;
namespace json=boost::json;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* why){try{f();}catch(const std::exception&){return;}throw std::runtime_error(why);}
int main(){try {
  B::set_precision(384);ExactField field({"x","eps","I"});Exact zero(field,"0"),one(field,"1"),eps(field,"eps"),x(field,"x");
  B input=B::from_strings("1/3","-7/13");arb_add_error_2exp_si(acb_realref(input.raw()),-200);
  arb_add_error_2exp_si(acb_imagref(input.raw()),-250);
  LaurentRows rows{-2,3,std::vector(2,std::vector(2,std::vector<B>(6,B(0))))};
  rows.coefficients[0][0][0]=input;rows.coefficients[1][1][5]=B::from_strings("19/23");
  recursion::PreparedAdjointStage stage{{{one/(one+x),eps},{zero,zero}},{{x/eps,zero},{zero,one}},
    {{-x/eps,zero},{zero,-one}},rows,rows,{zero,Exact(field,"1/3")},{one,Exact(field,"1/3")},{0,1}};
  auto payload=stage_payload(stage,one.constant(2)-one.constant(2)*eps,{{1,0},{1,1}});
  json::object report{{"prepared_stage",payload},{"prepared_stage_sha256",artifacts::detail::sha256(artifacts::detail::canonical(payload))}};
  // Load at another ambient precision. No midpoint/radius may be rounded
  // during serialization, including tiny imaginary radii and exact zeros.
  B::set_precision(64);auto restored=read_stage(json::parse(json::serialize(report)).as_object());
  for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j)for(unsigned k=0;k<6;++k)
    require(acb_equal(rows.coefficients[i][j][k].raw(),restored.lower_endpoint.coefficients[i][j][k].raw()),"checkpoint changed stored ball bits");
  require(B::precision()==64,"checkpoint read changed ambient precision");
  require(stage.connection[0][0].str()==restored.connection[0][0].str() &&
    stage.upper_path[0].str()==restored.upper_path[0].str(),"checkpoint changed exact connection/path");
  B::set_precision(384);AdjointOptions options;options.taylor_order=24;
  auto fresh=transport_adjoint_rows(stage.connection,stage.lower_endpoint,stage.lower_forcing,stage.lower_path,options);
  auto replay=transport_adjoint_rows(restored.connection,restored.lower_endpoint,restored.lower_forcing,restored.lower_path,options);
  for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j)for(unsigned k=0;k<6;++k)
    require(acb_equal(fresh.coefficients[i][j][k].raw(),replay.coefficients[i][j][k].raw()),"replay differs from original retained transport");
  auto corrupt=report;corrupt["prepared_stage"].as_object()["dimension"]="4-2*eps";
  rejects([&]{read_stage(corrupt);},"checkpoint payload corruption accepted");
  auto malformed=exact_rows(rows);malformed["high"]=4;
  rejects([&]{read_rows(malformed);},"checkpoint accepted missing epsilon coefficients");
  malformed=exact_rows(rows);
  malformed["coefficients"].as_array()[0].as_array()[0].as_array()[0].as_array()[0]="garbage";
  rejects([&]{read_rows(malformed);},"malformed ball accepted");
  auto invalid=rows;acb_indeterminate(invalid.coefficients[0][0][0].raw());
  rejects([&]{read_rows(exact_rows(invalid));},"nonfinite checkpoint accepted");
  require(one.parse("x/(1+eps)")==x/(one+eps),"existing-field parse changed context or value");
  const std::vector<Exact> vertices{zero,one.constant(3),one.constant(4)};
  auto clipped=remaining_path(vertices,0,0.1);
  require(clipped.front()==one.parse("10808639105689191/36028797018963968"),"chart parameter was not preserved as an exact dyadic");
  require(remaining_path(vertices,0,1)==std::vector<Exact>({one.constant(3),one.constant(4)}),"completed leg was not removed");
  auto tiny=remaining_path(vertices,0,std::numeric_limits<double>::denorm_min());
  require(tiny.front()*one.constant(2).pow(1074)==one.constant(3),"subnormal chart parameter was rounded");
  rejects([&]{remaining_path(vertices,0,0);},"nonprogressing chart parameter accepted");
  // Interrupt after one accepted chart, serialize the state, and resume on
  // the remaining exact path. Reparameterization can change chart boundaries;
  // compare against y(x)=y(0)/(1+x), not bit equality of different truncations.
  ExactEpsilonMatrix scalar_connection{{one/(one+x)}},scalar_forcing{{zero}};
  B uncertain(1);arb_add_error_2exp_si(acb_realref(uncertain.raw()),-200);
  LaurentRows scalar_rows{0,0,{{{uncertain}}}},saved_rows=scalar_rows;
  std::vector<Exact> rest;unsigned completed=0;
  AdjointOptions interrupt_options;interrupt_options.taylor_order=80;
  interrupt_options.chart_observer=[&](unsigned leg,double from,double to,const LaurentRows& accepted) {
    require(to>from,"observer received an unaccepted chart");
    ++completed;saved_rows=read_rows(exact_rows(accepted));rest=remaining_path(vertices,leg,to);
    throw std::runtime_error("simulated interruption after saved chart");
  };
  rejects([&]{transport_adjoint_rows(scalar_connection,scalar_rows,scalar_forcing,vertices,interrupt_options);},"simulated interruption did not propagate");
  require(completed==1 && rest.size()==3,"completed chart snapshot missing");
  options=AdjointOptions{};
  auto resumed=transport_adjoint_rows(scalar_connection,saved_rows,scalar_forcing,rest,options);
  auto continuous=transport_adjoint_rows(scalar_connection,scalar_rows,scalar_forcing,vertices,options);
  const auto expected=uncertain/B(5);
  require(NativeTailMagnitude::upper_abs(resumed.coefficients[0][0][0]-expected).approximate_upper()<1e-45,"resumed transport disagrees with analytic solution");
  require(NativeTailMagnitude::upper_abs(continuous.coefficients[0][0][0]-expected).approximate_upper()<1e-45,"continuous transport disagrees with analytic solution");
  auto finished=transport_adjoint_rows(scalar_connection,resumed,scalar_forcing,{vertices.back()},options);
  require(acb_equal(resumed.coefficients[0][0][0].raw(),finished.coefficients[0][0][0].raw()),"completed-arm resume changed saved result");
  // A shared uncertain leaf must stay shared after persistence. Subtracting
  // equal large maps cancels exactly before the uncertain source is applied.
  B shared_value(1);arb_add_error_2exp_si(acb_realref(shared_value.raw()),-100);
  B large=B::from_strings("1000000000000000000000000000000000000000000000000000000000000");
  linear_boundary::Expression expression{LaurentRows{0,0,{{{large}},{{large}}}},
    std::make_shared<const LaurentBoundary>(LaurentBoundary{0,{{shared_value}},true})};
  auto encoded_expression=exact_expression(expression);
  auto expression_hash=artifacts::detail::sha256(artifacts::detail::canonical(encoded_expression));
  B::set_precision(64);auto decoded_expression=read_expression(encoded_expression,expression_hash);
  require(acb_equal(decoded_expression.leaf_source->values[0][0].raw(),shared_value.raw()),"expression source balls changed during loading");
  require(!decoded_expression.leaf_source->taylor_tail_certified,"uncertified checkpoint promoted a producer's certification claim");
  B::set_precision(384);
  LaurentRows subtraction{0,0,{{{B(1)},{B(-1)}}}};
  auto difference=linear_boundary::materialize(linear_boundary::compose(subtraction,decoded_expression,0),0);
  require(difference.values[0][0].is_zero(),"serialized expression lost common-source cancellation");
  auto separate=linear_boundary::materialize(decoded_expression,0);
  require(NativeTailMagnitude::upper_abs(separate.values[0][0]-separate.values[1][0]).approximate_upper()>1e20,"common-source test did not distinguish premature materialization");
  auto altered=encoded_expression;altered["producer_leaf_tail_certified"]=false;
  rejects([&]{read_expression(altered,expression_hash);},"expression checksum corruption accepted");
  LaurentRows before_shift{-1,2,{{{B(1),B(2),B(3),B(4)},{B(5),B(6),B(7),B(8)}}}};
  auto shifted=shift_laurent_columns(before_shift,{0,2},true);
  require(shifted.low==-3 && shifted.high==0 && shifted.coefficients[0][0][0].is_zero() &&
    shifted.coefficients[0][0][1].is_zero() && acb_equal_si(shifted.coefficients[0][0][2].raw(),1) &&
    acb_equal_si(shifted.coefficients[0][1][3].raw(),8),"inverse gauge did not preserve the known Laurent window");
  rejects([&]{shift_laurent_columns(before_shift,{0,std::numeric_limits<std::int64_t>::min()},true);},"invalid inverse shift overflow accepted");
  std::cout<<"Prepared stage: lossless balls, exact data, bit-identical replay, corruption rejection and interrupted/resumed analytic comparison passed\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
