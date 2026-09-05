#include "diffexp/recursion_pipeline.hpp"
#include <iostream>
using namespace diffexp;using B=Jet::Ball;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
double upper(const B& value){return NativeTailMagnitude::upper_abs(value).approximate_upper();}
level::Result native_sunrise(const ibp::PropagatorBasis& basis,const Exact& dimension,
    const ExactField&,std::size_t parameter,const std::vector<ibp::Integral>& requested,const level::Options&) {
  if(basis.physical_count!=2 || basis.denominators.size()!=5)throw std::invalid_argument("test closure only supports merged sunrise");
  ibp::Generator generator(basis,dimension);ibp::ExactReducer reducer(dimension,300);
  ibp::for_each_seed(2,5,{1,1,50},[&](const ibp::Integral& seed){for(auto& row:generator.relations(seed))reducer.insert(std::move(row));});
  ibp::BasisReduction coordinates(reducer,{{1,0,0,0,0},{1,1,0,0,0},{1,1,-1,0,0}},dimension);
  level::Result result;result.success=true;result.ordered_basis=coordinates.ordered_basis();
  for(const auto& master:result.ordered_basis)result.matrix.push_back(coordinates.resolve(generator.derivative(master,parameter)));
  for(const auto& target:requested)result.target_rows.push_back(coordinates.resolve({{target,dimension.constant(1)}}));
  return result;
}
int main(){try{
  auto graph=recursion::prepare(feynman::example_family("sunrise"),{{1,1,1,0,0},{2,1,1,0,0},{0,1,1,0,0}},{},native_sunrise);
  recursion::NumericalOptions settings;settings.working_bits=256;settings.endpoint_order=16;
  settings.ordinary_order=48;settings.leaf_digits=18;
  auto zero=graph.dimension.constant(0),one=graph.dimension.constant(1),eps=graph.dimension.variable(1);
  // A direct parent with a pole used to finish child output0 first, discover
  // missing orders only afterward, and repeat the child's entire transport.
  auto parent=graph.nodes[0];parent.closure.ordered_basis=parent.requested;
  parent.requested.resize(1);parent.operations.resize(1);parent.operations[0].operation=feynman::Operation::Direct;
  parent.observable_rows={{one/eps.pow(2),zero,zero}};
  graph.nodes.insert(graph.nodes.begin(),std::move(parent));
  for(auto method:{recursion::LinearMethod::adjoint,recursion::LinearMethod::factored}) {
    settings.linear_method=method;unsigned lower=0,upper_count=0,leaves=0;
    settings.progress=[&](std::size_t depth,const std::string& phase,int high){
      if(depth==1 && (phase=="lower adjoint transport" || phase=="lower factored transport"))++lower;
      if(depth==1 && (phase=="upper adjoint transport" || phase=="middle factored transport"))++upper_count;
      if(phase=="certified scalar leaf")++leaves;
    };
    recursion::Evaluator solver(graph,settings);
    bool unavailable=false;try{solver.linear_expression();}catch(const std::logic_error&){unavailable=true;}
    require(unavailable,"unevaluated linear expression was exposed");
    const auto ambient=B::precision();auto prepared=solver.prepare_adjoint_stage(1,5);
    require(B::precision()==ambient,"boundary-independent preparation leaked working precision");
    require(solver.statistics().leaf_evaluations==0 && solver.statistics().numeric_evaluations==0,"endpoint preparation evaluated a child");
    require(prepared.connection.size()==3 && prepared.lower_endpoint.columns()==3 &&
      prepared.lower_endpoint.high>=5 && prepared.upper_endpoint.high>=5,"prepared local stage shape/window");
    for(unsigned i=0;i<prepared.lower_forcing.size();++i)for(unsigned j=0;j<3;++j)
      require(prepared.lower_forcing[i][j]==-prepared.upper_forcing[i][j],"prepared upper arm lost forcing sign");
    auto value=solver.evaluate(0);
    auto expression=solver.linear_expression();
    const auto saved_precision=B::precision();B::set_precision(settings.working_bits);
    auto rematerialized=linear_boundary::materialize(expression,0);B::set_precision(saved_precision);
    require(rematerialized.low==value.low && acb_equal(rematerialized.values[0][-value.low].raw(),value.values[0][-value.low].raw()),"linear snapshot changed completed result");
    const auto saved_coefficient=expression.transform.coefficients[0][0][0];
    expression.transform.coefficients[0][0][0]=B(123);
    require(acb_equal(solver.linear_expression().transform.coefficients[0][0][0].raw(),saved_coefficient.raw()),"caller edit changed the evaluator's saved transform");
    require(lower==1 && upper_count==1 && leaves==1,"preflight repeated a transport or leaf evaluation");
    require(solver.statistics().demand_preflights==1 && solver.statistics().exact_plans==1,"structural demand plan not cached");
    auto child=solver.evaluate(1,2);
    require(upper(value.values[0][-value.low]-child.values[0][2-child.low])<1e-12,"pole preflight changed output");
    const auto before=solver.statistics();(void)solver.evaluate(0);
    require(solver.statistics().numeric_evaluations==before.numeric_evaluations,"preflight high-watermark cache missed");
  }
  // All-positive and zero rows must retain structural bounds even when the
  // order-zero operator probe has no nonzero retained coefficients.
  auto positive=graph;positive.nodes[0].observable_rows={{eps.pow(2),zero,zero}};
  settings.linear_method=recursion::LinearMethod::adjoint;settings.progress={};
  recursion::Evaluator plus(positive,settings);auto p=plus.evaluate(0);
  require(p.values[0][-p.low].is_finite(),"positive-order direct operator failed preflight");
  auto empty=graph;empty.nodes[0].observable_rows={{zero,zero,zero}};
  recursion::Evaluator zeros(empty,settings);auto z=zeros.evaluate(0);
  for(const auto& c:z.values[0])require(c.is_zero(),"structurally zero preflight result changed");
  std::cout<<"Laurent-demand preflight: one child transport per forced method, direct poles, positive/zero rows and high-watermark reuse passed\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
