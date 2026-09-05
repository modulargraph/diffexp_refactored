// Validation utility for reusing a completed ordinary transport when adding
// explicit algebraic basis data. No reference values enter this conversion.
// Build with the same include/library flags as the diffexp executable.
#include "diffexp/transport.hpp"
#include <fstream>
#include <iostream>
using namespace diffexp;
namespace j=boost::json;
using B=kernel::ComplexBall;
j::object read(const char* path) {
  std::ifstream stream(path);if(!stream)throw std::runtime_error("cannot open input");
  return j::parse(std::string((std::istreambuf_iterator<char>(stream)),{})).as_object();
}
Boundary decode(const j::value& input,bool error_bounds) {
  Boundary result;
  for(const auto& row:input.as_array()) {
    std::vector<B> values;
    for(const auto& cell:row.as_array()) {
      const auto& o=cell.as_object();
      auto value=B::from_strings(transport::string(o.at("real_midpoint")),transport::string(o.at("imaginary_midpoint")));
      auto radius=B::from_strings(transport::string(o.at("radius")));
      if(!radius.is_finite() || !arb_is_nonnegative(acb_realref(radius.raw())))throw std::runtime_error("invalid radius");
      arb_add_error(acb_realref(value.raw()),acb_realref(radius.raw()));
      if(!error_bounds)arb_add_error(acb_imagref(value.raw()),acb_realref(radius.raw()));
      values.push_back(std::move(value));
    }
    result.push_back(std::move(values));
  }return result;
}
int main(int argc,char** argv) {
  try {
    if(argc!=5)throw std::invalid_argument("usage: rebase_completed_transport OLD_REQUEST NEW_REQUEST RAW_RESPONSE OUTPUT");
    auto old=read(argv[1]),request=read(argv[2]),response=read(argv[3]);
    auto without_prefactors=request;without_prefactors.erase("basis_prefactors");
    if(old!=without_prefactors || !request.if_contains("basis_prefactors"))throw std::invalid_argument("requests differ beyond explicit basis prefactors");
    if(old.if_contains("asymptotic") || old.if_contains("initial_only") || response.if_contains("basis_convention") ||
       response.at("schema")!="DiffExp.TransportResult/v1" || response.at("parameter")!="1")
      throw std::invalid_argument("requires an ordinary completed continuous-basis response");
    auto started=std::chrono::steady_clock::now();
    const auto bits=transport::integer(request,"working_bits",384,64,100000);
    const auto d=transport::integer(request,"dimension",-1,1,1000),k=transport::integer(request,"epsilon_order",4,0,100);
    B::set_precision(bits);auto compiled=transport::compile(request,d,k,bits);
    auto roots=transport::principal_roots_at(compiled,B(0));
    for(unsigned i=0;i<roots.size();++i)roots[i]=continue_polynomial_sqrt(compiled.square_polynomials[i],B(0),B(1),roots[i]);
    auto values=decode(response.at("values"),false),errors=decode(response.at("errors"),true);
    transport::apply_basis_prefactors(compiled,B(1),roots,values,errors);
    const auto goal=transport::integer(request,"accuracy_goal",0,0,20000);
    if(goal>0)for(const auto& row:values)for(const auto& value:row) {
      auto bound=(transport::magnitude(value)+B(1))*B::from_strings("1e-"+std::to_string(goal));
      if(!arb_le(acb_realref(transport::arithmetic_error(value).raw()),acb_realref(bound.raw())))
        throw std::runtime_error("converted estimate exceeds AccuracyGoal: value="+j::serialize(transport::number(value,30))+", bound="+j::serialize(transport::number(bound,30)));
    }
    const auto digits=bits*30103/100000+5;
    response["values"]=transport::matrix_json(values,digits);response["errors"]=transport::matrix_json(errors,digits);
    response["basis_convention"]="principal_endpoint";response["segments_basis_convention"]="continued";
    response["endpoint_roots_continued"]=transport::matrix_json(Boundary{roots},digits).at(0);
    response["endpoint_roots_principal"]=transport::matrix_json(Boundary{transport::principal_roots_at(compiled,B(1))},digits).at(0);
    response["basis_rebase"]=j::object{{"requests_identical_except_basis_prefactors",true},{"ode_repeated",false},
      {"seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()},
      {"ratios",transport::matrix_json(Boundary{transport::basis_prefactor_ratios(compiled,B(1),roots)},digits).at(0)}};
    std::ofstream output(argv[4]);if(!output)throw std::runtime_error("cannot open output");
    output<<j::serialize(response)<<'\n';if(!output)throw std::runtime_error("cannot write output");
    std::cout<<j::serialize(response.at("basis_rebase"))<<'\n';
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
