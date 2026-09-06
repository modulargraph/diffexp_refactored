#include "diffexp/transport.hpp"
#include <iostream>
using namespace diffexp;
namespace j=boost::json;
static double real(const j::value& result,unsigned i=0,unsigned e=0) {
 return std::stod(std::string(result.as_object().at("values").as_array().at(i).as_array().at(e).as_object().at("real_midpoint").as_string()));
}
int main() {
 try {
  // Mathematica uncertain zeros use absolute accuracy (two backticks).
  Jet decimal_context(0,1,256);
  auto uncertain_zero=decimal_context.decimal("0``40.5").at(0);
  auto positive=Jet::Ball::from_strings("1e-41"),negative=Jet::Ball::from_strings("-1e-41");
  if(!acb_contains(uncertain_zero.raw(),positive.raw()) || !acb_contains(uncertain_zero.raw(),negative.raw()) || uncertain_zero.is_zero())
    throw std::runtime_error("absolute-accuracy zero lost its uncertainty");
  auto coarse_zero=evaluate(data::Reader("0``-2").read(),decimal_context,{}).at(0);
  if(!acb_contains(coarse_zero.raw(),Jet::Ball(10).raw()))throw std::runtime_error("negative absolute accuracy rejected");
  auto uncertain_value=decimal_context.decimal("100``4").at(0);
  auto inside=Jet::Ball::from_strings("100.00001"),outside=Jet::Ball::from_strings("100.001");
  if(!acb_contains(uncertain_value.raw(),inside.raw()) || acb_overlaps(uncertain_value.raw(),outside.raw()))
    throw std::runtime_error("absolute accuracy was treated as relative precision");
  auto request=j::parse(R"JSON({"dimension":1,"epsilon_order":0,"taylor_order":40,"working_bits":256,"division_order":4,"paths":{"t":"x/2"},"entries":[{"row":0,"column":0,"epsilon":0,"variable":"t","expression":"1/(1-t)"}],"boundary":[["1"]]})JSON");
  auto result=transport::run(request);if(std::abs(real(result)-2)>1e-15)throw std::runtime_error("ordinary physical-basis transport");
  auto& r=request.as_object();r["paths"]=j::object{{"t","x"}};
  r["entries"]=j::array{j::object{{"row",0},{"column",0},{"epsilon",0},{"variable","dlog"},{"expression","1+Sqrt[1+t]"}}};
  result=transport::run(request);if(std::abs(real(result)-(1+std::sqrt(2.0))/2)>1e-15)throw std::runtime_error("spurious conjugate-sheet root");
  r["paths"]=j::object{{"t","1-8*x*(1-x)+8*I*x*(1-x)*(1-2*x)"}};
  r["entries"].as_array()[0].as_object()["expression"]="2+Sqrt[t]";
  result=transport::run(request);if(std::abs(real(result)-1.0/3)>1e-12)throw std::runtime_error("continuous square-root monodromy");
  request=j::parse(R"JSON({"dimension":2,"epsilon_order":0,"taylor_order":32,"working_bits":256,"paths":{"t":"x"},"entries":[{"row":0,"column":1,"epsilon":0,"variable":"t","expression":"1/t"},{"row":1,"column":0,"epsilon":0,"variable":"t","expression":"1/t"}],"initial_only":true,"asymptotic":{"constraints":[{"row":0,"epsilon":0,"power":"1","log_degree":0,"value":"1"}],"cutoffs":[{"row":0,"power":"2"}]}})JSON");
  result=transport::run(request);auto parameter=std::stod(std::string(result.as_object().at("parameter").as_string()));
  if(std::abs(real(result)-parameter)>1e-15 || std::abs(real(result,1)-parameter)>1e-15)throw std::runtime_error("coupled partial asymptotic matching");
  request.as_object()["asymptotic"].as_object()["cutoffs"]=j::array{};
  bool rejected=false;try{transport::run(request);}catch(const std::invalid_argument&){rejected=true;}if(!rejected)throw std::runtime_error("underdetermined asymptotic constants accepted");
  auto simple=j::parse(R"JSON({"dimension":1,"epsilon_order":0,"taylor_order":8,"working_bits":256,"paths":{"t":"x"},"entries":[{"row":0,"column":0,"epsilon":0,"variable":"t","expression":"1"}],"initial_only":true,"asymptotic":{"constraints":[{"row":0,"epsilon":0,"power":"0","log_degree":0,"value":"1"}],"cutoffs":[{"row":0,"power":"2"}]}})JSON");
  rejected=false;try{transport::run(simple);}catch(const std::invalid_argument&){rejected=true;}if(!rejected)throw std::runtime_error("cutoff-implied zero ignored");
  simple.as_object()["asymptotic"].as_object()["cutoffs"].as_array()[0].as_object()["power"]="1";
  simple.as_object()["accuracy_goal"]=30;result=transport::run(simple);
  if(std::stod(std::string(result.as_object().at("values").as_array()[0].as_array()[0].as_object().at("radius").as_string()))>1e-30)throw std::runtime_error("initial-only accuracy goal ignored");
  auto malformed=simple;malformed.as_object()["asymptotic"].as_object()["constraints"].as_array()[0].as_object().erase("row");
  rejected=false;try{transport::run(malformed);}catch(const std::invalid_argument&){rejected=true;}if(!rejected)throw std::runtime_error("missing boundary row accepted");
  simple.as_object()["epsilon_order"]=1;simple.as_object()["accuracy_goal"]=0;simple.as_object()["entries"].as_array()[0].as_object()["epsilon"]=1;
  simple.as_object()["asymptotic"].as_object()["constraints"].as_array().push_back(j::object{{"row",0},{"epsilon",1},{"power","1"},{"log_degree",0},{"value","1"}});
  result=transport::run(simple);parameter=std::stod(std::string(result.as_object().at("parameter").as_string()));
  if(std::abs(real(result,0,1)-parameter)>1e-15)throw std::runtime_error("explicit constraint above zero-fill cutoff lost");
  auto canonical=j::parse(R"JSON({"dimension":3,"epsilon_order":2,"taylor_order":20,"working_bits":256,"division_order":2,"paths":{"t":"x"},"entries":[{"row":0,"column":1,"epsilon":1,"variable":"dlog","expression":"t"},{"row":1,"column":2,"epsilon":1,"variable":"dlog","expression":"100+t"}],"boundary":[["0","0","0"],["0","0","0"],["1","0","0"]]})JSON");
  result=transport::run(canonical);if(std::abs(real(result,1,1)-std::log(1.01))>1e-15)throw std::runtime_error("singular initial chart overshot endpoint");
  auto saved=j::parse(R"JSON({"dimension":1,"epsilon_order":0,"taylor_order":1000,"working_bits":100000,"save_segments":true,"paths":{"t":"x"},"entries":[],"boundary":[["1"]]})JSON");
  rejected=false;try{transport::run(saved);}catch(const std::length_error&){rejected=true;}if(!rejected)throw std::runtime_error("saved segment output budget missing");
  // A^2=0 resolves exactly before independent input radii are applied.
  auto nilpotent=j::parse(R"JSON({"paths":{"t":"x"},"entries":[{"row":0,"column":0,"epsilon":0,"variable":"t","expression":"1000000000000"},{"row":0,"column":1,"epsilon":0,"variable":"t","expression":"-1000000000000"},{"row":1,"column":0,"epsilon":0,"variable":"t","expression":"1000000000000"},{"row":1,"column":1,"epsilon":0,"variable":"t","expression":"-1000000000000"},{"row":0,"column":0,"epsilon":1,"variable":"t","expression":"1000000000000"},{"row":0,"column":1,"epsilon":1,"variable":"t","expression":"-1000000000000"},{"row":1,"column":0,"epsilon":1,"variable":"t","expression":"1000000000000"},{"row":1,"column":1,"epsilon":1,"variable":"t","expression":"-1000000000000"}]})JSON");
  using B=transport::B;B::set_precision(384);
  auto compiled=transport::compile(nilpotent.as_object(),2,1,384);auto entries=transport::numerical_entries(compiled);
  Boundary uncertain{{B(1),B(2)},{B(1),B(2)}};
  for(auto& row:uncertain)for(auto& v:row)arb_add_error_2exp_si(acb_realref(v.raw()),-100);
  const B h=B::from_strings("1/16"),scale=B::from_strings("1000000000000"),rho=B::from_strings("1/1267650600228229401496703205376");
  auto centered=transport::chart(compiled,entries,uncertain,{},B(0),h,16,true,120);
  auto naive=transport::chart(compiled,entries,uncertain,{},B(0),h,16,false,120,false);
  auto capped=transport::chart(compiled,entries,uncertain,{},B(0),h,16,false,120,true,1);
  if(!centered.centered || capped.centered)throw std::runtime_error("centered chart selector or map workspace cap");
  for(unsigned i=0;i<2;++i)for(unsigned e=0;e<2;++e) {
    if(!acb_equal(capped.values[i][e].raw(),naive.values[i][e].raw()))throw std::runtime_error("workspace fallback changed input radii");
    if(mag_get_d(arb_radref(acb_realref(centered.values[i][e].raw())))>1e-17 ||
       mag_get_d(arb_radref(acb_realref(naive.values[i][e].raw())))<1)
      throw std::runtime_error("nilpotent uncertainty wrapping was not removed");
  }
  for(unsigned corner=0;corner<16;++corner) {
    Boundary point{{B(1),B(2)},{B(1),B(2)}};
    for(unsigned i=0;i<2;++i)for(unsigned e=0;e<2;++e)point[i][e]+=rho*B(corner&(1u<<(2*i+e))?1:-1);
    for(unsigned i=0;i<2;++i)for(unsigned e=0;e<2;++e) {
      auto expected=point[i][e]+h*scale*(point[0][e]-point[1][e]+(e?point[0][0]-point[1][0]:B(0)));
      if(!acb_contains(centered.values[i][e].raw(),expected.raw()))throw std::runtime_error("centered epsilon map dropped an input corner");
    }
  }
  auto triangular=j::parse(R"JSON({"paths":{"t":"x"},"entries":[{"row":0,"column":1,"epsilon":0,"variable":"t","expression":"7"}]})JSON");
  auto triangular_compiled=transport::compile(triangular.as_object(),2,1,384);
  auto triangular_chart=transport::chart(triangular_compiled,transport::numerical_entries(triangular_compiled),uncertain,{},B(0),h,16,false,120);
  if(triangular_chart.centered)throw std::runtime_error("triangular recurrence built an unnecessary homogeneous map");
  for(unsigned e=0;e<2;++e) {
    auto expected=uncertain[0][e]+B(7)*h*uncertain[1][e];
    if(!acb_overlaps(triangular_chart.values[0][e].raw(),expected.raw()))throw std::runtime_error("triangular recurrence lost uncertainty");
  }
  // A zero midpoint still has a nonzero uncertainty contribution to the tail estimate.
  auto scalar=j::parse(R"JSON({"paths":{"t":"x"},"entries":[{"row":0,"column":0,"epsilon":0,"variable":"t","expression":"1"}]})JSON");
  auto scalar_compiled=transport::compile(scalar.as_object(),1,0,384);B zero;arb_add_error_2exp_si(acb_realref(zero.raw()),-100);
  auto uncertain_tail=transport::chart(scalar_compiled,transport::numerical_entries(scalar_compiled),{{zero}},{},B(0),B::from_strings("1/4"),12,true,120);
  B tail_bound(0),term=rho;for(unsigned n=1;n<=12;++n){term=term/B(4)/B(n);if(n>=9)tail_bound+=term;}
  if(!uncertain_tail.centered || !arb_contains(acb_realref(uncertain_tail.truncation_errors[0][0].raw()),acb_realref(tail_bound.raw())))
    throw std::runtime_error("centered chart dropped uncertainty from truncation estimate");
  if(std::stod(std::string(uncertain_tail.saved.at("coefficients").as_array()[0].as_array()[0].as_array()[0].as_object().at("radius").as_string()))<=0)
    throw std::runtime_error("saved centered coefficient dropped uncertainty");
  auto canonical_map_request=j::parse(R"JSON({"paths":{"t":"x"},"entries":[{"row":0,"column":0,"epsilon":1,"variable":"dlog","expression":"1+t"}]})JSON");
  auto canonical_map=transport::compile(canonical_map_request.as_object(),1,1,384);
  auto canonical_uncertain=transport::chart(canonical_map,transport::numerical_entries(canonical_map),{{uncertain[0][0],B(0)}},{},B(0),h,12,false,120);
  auto canonical_original=transport::chart(canonical_map,transport::numerical_entries(canonical_map),{{uncertain[0][0],B(0)}},{},B(0),h,12,false,120,false);
  if(canonical_uncertain.centered)throw std::runtime_error("canonical fast path expanded a homogeneous map");
  for(unsigned e=0;e<2;++e)if(!acb_equal(canonical_uncertain.values[0][e].raw(),canonical_original.values[0][e].raw()))
    throw std::runtime_error("canonical fast path changed retained balls");
  // Repeated source/letter pairs use the same dot result without changing
  // row accumulation order, including nonunit and negative rational weights.
  auto repeated=j::parse(R"JSON({"paths":{"t":"x"},"entries":[{"row":0,"column":0,"epsilon":1,"variable":"dlog","expression":"1+t","coefficient":"2/3"},{"row":1,"column":0,"epsilon":1,"variable":"dlog","expression":"1+t","coefficient":"-7/2"},{"row":1,"column":1,"epsilon":1,"variable":"dlog","expression":"1+t"},{"row":0,"column":1,"epsilon":1,"variable":"dlog","expression":"2+t","coefficient":"-1"}]})JSON");
  auto repeated_compiled=transport::compile(repeated.as_object(),2,2,384);
  Boundary repeat_boundary{{uncertain[0][0],B(2),B(3)},{B(0),B::from_strings("1/3"),B(0)}};
  auto repeated_entries=transport::numerical_entries(repeated_compiled);
  auto shared=transport::chart(repeated_compiled,repeated_entries,repeat_boundary,{},B(0),h,40,true,120);
  auto unshared=transport::chart(repeated_compiled,repeated_entries,repeat_boundary,{},B(0),h,40,true,120,true,10000000,false);
  for(unsigned i=0;i<2;++i)for(unsigned e=0;e<3;++e)
    if(!acb_equal(shared.values[i][e].raw(),unshared.values[i][e].raw()) || !acb_equal(shared.errors[i][e].raw(),unshared.errors[i][e].raw()))
      throw std::runtime_error("shared canonical convolution changed retained balls or errors");
  if(shared.saved!=unshared.saved)throw std::runtime_error("shared canonical convolution changed saved coefficients");
  auto direct=j::parse(R"JSON({"dimension":1,"epsilon_order":0,"taylor_order":8,"working_bits":256,"accuracy_goal":12,"division_order":4,"paths":{"t":"x"},"entries":[{"row":0,"column":0,"epsilon":0,"variable":"t","expression":"1/(2*(1+t))"}],"asymptotic":{"constraints":[{"row":0,"epsilon":0,"power":"0","log_degree":0,"value":"1"}],"cutoffs":[{"row":0,"power":"1"}]}})JSON");
  result=transport::run(direct);const auto& cell=result.as_object().at("values").as_array()[0].as_array()[0].as_object();
  B midpoint=B::from_strings(std::string(cell.at("real_midpoint").as_string())),radius=B::from_strings(std::string(cell.at("radius").as_string())),oracle;
  acb_sqrt(oracle.raw(),B(2).raw(),256);arb_add_error(acb_realref(midpoint.raw()),acb_realref(radius.raw()));
  if(!acb_overlaps(midpoint.raw(),oracle.raw()))throw std::runtime_error("direct asymptotic continuation lost initial uncertainty");
  auto unknown=direct;unknown.as_object()["accuracy_gola"]=12;rejected=false;
  try{transport::run(unknown);}catch(const std::invalid_argument&){rejected=true;}if(!rejected)throw std::runtime_error("unknown transport key accepted");
  auto principal_basis=j::parse(R"JSON({"dimension":1,"epsilon_order":0,"taylor_order":40,"working_bits":256,"paths":{"t":"1-8*x*(1-x)+8*I*x*(1-x)*(1-2*x)"},"entries":[{"row":0,"column":0,"epsilon":0,"variable":"dlog","expression":"2+Sqrt[t]"}],"boundary":[["1"]],"basis_prefactors":["2+Sqrt[t]"]})JSON");
  result=transport::run(principal_basis);
  if(std::abs(real(result)-1)>1e-12 || result.as_object().at("basis_convention")!="principal_endpoint")
    throw std::runtime_error("supplied principal-basis conversion lost algebraic loop normalization");
  // A root occurring only in the supplied basis data must also be continued.
  principal_basis.as_object()["entries"]=j::array{};
  principal_basis.as_object()["epsilon_order"]=1;
  principal_basis.as_object()["boundary"]=j::array{j::array{"1","2"}};
  principal_basis.as_object()["boundary_errors"]=j::array{j::array{"1/100000000000000000000","1/50000000000000000000"}};
  principal_basis.as_object()["basis_prefactors"]=j::array{"Sqrt[t]"};
  result=transport::run(principal_basis);
  if(std::abs(real(result)+1)>1e-12 || std::abs(real(result,0,1)+2)>1e-12)
    throw std::runtime_error("prefactor-only root or epsilon conversion lost monodromy");
  for(unsigned e=0;e<2;++e) {
    const auto& v=result.as_object().at("values").as_array()[0].as_array()[e].as_object();
    const auto& error=result.as_object().at("errors").as_array()[0].as_array()[e].as_object();
    if(std::stod(std::string(v.at("radius").as_string()))<1e-20 || std::stod(std::string(error.at("real_midpoint").as_string()))<1e-20)
      throw std::runtime_error("principal-basis conversion dropped carried errors");
  }
  // Reuse the exact helper independently of any differential-equation run.
  principal_basis.as_object()["basis_prefactors"]=j::array{"2+Sqrt[t]"};
  auto prefactor_compiled=transport::compile(principal_basis.as_object(),1,1,256);
  auto continued=transport::principal_roots_at(prefactor_compiled,B(0));
  for(unsigned i=0;i<continued.size();++i)continued[i]=continue_polynomial_sqrt(prefactor_compiled.square_polynomials[i],B(0),B(1),continued[i]);
  Boundary completed{{B(2),B(4)}},completed_errors{{B::from_strings("1/100"),B::from_strings("1/50")}};
  for(auto& v:completed[0])arb_add_error_2exp_si(acb_realref(v.raw()),-40);
  transport::apply_basis_prefactors(prefactor_compiled,B(1),continued,completed,completed_errors);
  if(!acb_contains(completed[0][0].raw(),B(6).raw()) || !acb_contains(completed[0][1].raw(),B(12).raw()) ||
     !arb_ge(acb_realref(completed_errors[0][0].raw()),acb_realref(B::from_strings("3/100").raw())))
    throw std::runtime_error("reusable principal-basis helper failed to scale values/errors");
  auto invalid_prefactor=principal_basis;invalid_prefactor.as_object()["basis_prefactors"]=j::array{"0"};
  rejected=false;try{transport::run(invalid_prefactor);}catch(const std::invalid_argument&){rejected=true;}
  if(!rejected)throw std::runtime_error("zero endpoint basis divisor accepted");
  invalid_prefactor.as_object()["basis_prefactors"]=j::array{};
  rejected=false;try{transport::run(invalid_prefactor);}catch(const std::invalid_argument&){rejected=true;}
  if(!rejected)throw std::runtime_error("basis prefactor shape mismatch accepted");
  // Exact endpoint cancellation must precede principal-root evaluation on its cut.
  auto cancel_cut=j::parse(R"JSON({"paths":{"t":"x"},"entries":[],"basis_prefactors":["Sqrt[-1+I*(t/3+t/7-10/21)]"]})JSON");
  auto cut_compiled=transport::compile(cancel_cut.as_object(),1,0,256);
  auto cut_roots=transport::principal_roots_at(cut_compiled,B(1));
  B plus_i=B::from_strings("0","1");
  if(!acb_equal(cut_roots[0].raw(),plus_i.raw()))throw std::runtime_error("exact endpoint cancellation straddled principal root cut");
  cancel_cut.as_object()["basis_prefactors"]=j::array{"Sqrt[-1+I*(t/3+t/7-5/21)]"};
  auto dyadic_compiled=transport::compile(cancel_cut.as_object(),1,0,256);
  auto dyadic_roots=transport::principal_roots_at(dyadic_compiled,B::from_strings("1/2"));
  if(!acb_equal(dyadic_roots[0].raw(),plus_i.raw()))throw std::runtime_error("exact dyadic endpoint cancellation lost principal branch");
  B interval_point(1);arb_add_error_2exp_si(acb_realref(interval_point.raw()),-80);
  auto interval_roots=transport::principal_roots_at(cut_compiled,interval_point);
  auto polynomial_value=cut_compiled.square_polynomials[0].evaluate_polynomial(interval_point);B enclosing_root;
  acb_sqrt(enclosing_root.raw(),polynomial_value.raw(),B::precision());
  if(!acb_equal(interval_roots[0].raw(),enclosing_root.raw()))throw std::runtime_error("non-exact endpoint lost enclosing root fallback");
  auto endpoint_ghost=j::parse(R"JSON({"dimension":1,"epsilon_order":0,"taylor_order":40,"working_bits":299,"paths":{"t":"1+3*x"},"entries":[{"row":0,"column":0,"epsilon":0,"variable":"t","expression":"1/(2*Sqrt[t]*(2+Sqrt[t]))"}],"boundary":[["1"]],"basis_prefactors":["2+Sqrt[t]"]})JSON");
  result=transport::run(endpoint_ghost);
  if(std::abs(real(result)-4.0/3)>1e-15)throw std::runtime_error("endpoint conjugate-sheet root blocked normalized chaining");
  endpoint_ghost.as_object()["entries"].as_array()[0].as_object()["expression"]="1/(2*Sqrt[t]*(2-Sqrt[t]))";
  auto endpoint_pole=transport::compile(endpoint_ghost.as_object(),1,0,299);
  if(!std::any_of(endpoint_pole.singularities.begin(),endpoint_pole.singularities.end(),[](const B& root){return acb_contains(root.raw(),B(1).raw());}))
    throw std::runtime_error("domain clipping removed a genuine endpoint pole");
  rejected=false;try{transport::run(endpoint_ghost);}catch(const std::runtime_error&){rejected=true;}
  if(!rejected)throw std::runtime_error("genuine singular endpoint accepted after ghost filtering");
  std::cout<<"Generic physical/canonical/algebraic transport and asymptotic matching passed\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
