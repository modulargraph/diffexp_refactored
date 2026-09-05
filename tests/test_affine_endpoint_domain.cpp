#include "diffexp/recursion_pipeline.hpp"
#include <iostream>
using namespace diffexp;using B=Jet::Ball;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* why){try{f();}catch(const std::domain_error&){return;}throw std::runtime_error(why);}
template<class F>void rejects_physical(F f){try{f();}catch(const std::domain_error& error){
  require(std::string(error.what()).find("physical endpoint constraint failed")!=std::string::npos,"invalid physical boundary failed for an unrelated reason");return;
}throw std::runtime_error("nonzero physical divergent mode was silently projected away");}
double upper(const B& value){return NativeTailMagnitude::upper_abs(value).approximate_upper();}
recursion::Graph model(feynman::Operation operation,bool logarithm=false,bool wrong=false) {
  auto graph=recursion::prepare(feynman::example_family("bubble"),
    logarithm?std::vector<ibp::Integral>{{1,0},wrong?ibp::Integral{1,0}:ibp::Integral{0,0}}:
      std::vector<ibp::Integral>{{1,0},wrong?ibp::Integral{2,0}:ibp::Integral{0,1}});
  auto parent=graph.nodes.front();parent.scalar_leaf=false;parent.anchor=Rational("1/2");
  parent.closure.success=true;parent.closure.ordered_basis=graph.nodes.front().requested;
  auto z=graph.dimension.constant(0),one=z.constant(1),x=z.variable(0);
  if(logarithm)parent.closure.matrix={{z,one/x},{z,z}};
  else if(operation==feynman::Operation::UpperLimit)parent.closure.matrix={{one/(one-x),-one/(one-x)},{z,z}};
  else parent.closure.matrix={{-one/x,one/x},{z,z}};
  parent.requested={{1,1}};parent.observable_rows={{one,z}};
  parent.operations={pullback::Plan{}};parent.operations[0].operation=operation;
  graph.nodes.insert(graph.nodes.begin(),std::move(parent));return graph;
}
int main(){try{
  B::set_precision(384);ExactField field({"x","eps"});Exact x(field,"x"),zero(field,0),one(field,1);
  auto series=AffineFrobeniusSeries::prepare({{-one/x,one/x},{zero,zero}},0,1,8);
  auto row=series.project({one,zero});
  rejects([&]{series.dr_endpoint_constant(row);},"unconstrained strict endpoint silently accepted divergence");
  rejects([&]{series.dr_integral_from_zero(row);},"unconstrained strict integral silently accepted divergence");
  auto domain=series.dr_domain(row,false);
  require(domain.zero_constraints.size()==1 && domain.zero_constraints.front().power==Rational(-1),"fixed-power endpoint domain not extracted");
  require(series.dr_endpoint_constant(domain.admissible)[0][1]==one,"conditional finite endpoint coefficient changed");
  // An exact cross-column cancellation must be accepted by the strict API
  // after contraction. The domain rows retain those cross-column relations.
  AffineFrobeniusSeries::Expansion mixed{1,2,{{0,0,0,Rational(-1),Rational(0),one},{0,1,0,Rational(-1),Rational(0),-one},
      {0,0,0,Rational(0),Rational(0),one}}};
  auto md=series.dr_domain(mixed,false);require(md.zero_constraints.size()==1 && md.zero_constraints[0].coefficients==std::vector<Exact>{one,-one},"domain split lost cancellation between fundamental columns");
  require(series.dr_endpoint_constant(series.contract(mixed,{one,one}))[0][0]==one,"exact physical contraction lost finite limit");
  for(auto operation:{feynman::Operation::LowerLimit,feynman::Operation::UpperLimit,feynman::Operation::BetaIntegral})
    for(int method=0;method<3;++method) {
      auto graph=model(operation);recursion::NumericalOptions options;options.endpoint_order=12;options.ordinary_order=32;
      options.observable_adjoint=method!=0;options.linear_method=method==2?recursion::LinearMethod::factored:recursion::LinearMethod::adjoint;
      recursion::Evaluator evaluator(graph,options);auto actual=evaluator.evaluate(0);auto leaf=evaluator.evaluate(1,0);
      require(actual.values.size()==1 && !actual.taylor_tail_certified,"domain auxiliaries leaked into physical outputs or certification");
      for(int k=std::max(actual.low,leaf.low);k<=0;++k)require(upper(actual.values[0][k-actual.low]-leaf.values[0][k-leaf.low])<1e-20,"constrained endpoint/integral differs from exact constant physical solution");
      require(evaluator.statistics().endpoint_constraint_rows>0 && evaluator.statistics().endpoint_constraint_coefficients>0,"physical endpoint domain was not checked");
      if(options.observable_adjoint)require(evaluator.linear_expression().transform.coefficients.size()==1,"domain rows leaked into parent expression");
      auto invalid_graph=model(operation,false,true);recursion::Evaluator invalid(invalid_graph,options);
      rejects_physical([&]{invalid.evaluate(0);});
    }
  for(bool wrong:{false,true}) {
    auto graph=model(feynman::Operation::LowerLimit,true,wrong);recursion::NumericalOptions options;options.endpoint_order=12;options.ordinary_order=32;
    recursion::Evaluator evaluator(graph,options);
    if(wrong)rejects_physical([&]{evaluator.evaluate(0);});
    else {auto actual=evaluator.evaluate(0);auto leaf=evaluator.evaluate(1,0);require(upper(actual.values[0][-actual.low]-leaf.values[0][-leaf.low])<1e-20,"zero physical logarithm did not cancel");}
  }
  // Equal and opposite divergences at different endpoints do not define an
  // ordinary DR integral. Their constraints must never cancel each other.
  auto opposite=model(feynman::Operation::BetaIntegral);auto z=opposite.dimension.constant(0),o=z.constant(1),t=z.variable(0);
  opposite.nodes[0].closure.matrix={{-o/t,z},{z,o/(o-t)}};opposite.nodes[0].observable_rows={{o,-o}};
  for(auto method:{recursion::LinearMethod::adjoint,recursion::LinearMethod::factored}){
    recursion::NumericalOptions settings;settings.linear_method=method;settings.endpoint_order=12;settings.ordinary_order=32;
    recursion::Evaluator invalid(opposite,settings);rejects_physical([&]{invalid.evaluate(0);});
  }
  auto graph=model(feynman::Operation::LowerLimit);recursion::NumericalOptions options;options.endpoint_order=12;
  recursion::Evaluator stage(graph,options);rejects([&]{stage.prepare_adjoint_stage(0,0);},"standalone stage exported undisclosed endpoint constraints");
  std::cout<<"Physical endpoint domains: exact cross-column constraints, power/log cancellation, lower/upper/integral operations in values/adjoint/factored routes, invalid physical boundary rejection, auxiliary row isolation passed\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
