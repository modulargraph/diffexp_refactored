#include <diffexp/ibp.hpp>
#include <diffexp/families.hpp>
#include <diffexp/system.hpp>
#include <diffexp/deepest_beta.hpp>
#include <diffexp/recursion_leaf.hpp>
#include <diffexp/cached_affine.hpp>
#include <diffexp/recursion_pipeline.hpp>
#include <iostream>
int main() {
  using namespace diffexp;
  ExactField field({"x","eps"});Exact x(field,"x");
  const auto example=feynman::example_family("banana4");
  ibp::PropagatorBasis family(ibp::quadratic_family(example.momenta,x));
  RationalSystem system({{Exact(field,"1/(x-2)")}},0,1);
  auto local=regular_series(system,Rational(0),8,0,0,{{Rational("1/2")}});
  if(local.evaluate_polynomial(Rational(1))[0][0]!=Rational("1/4") || family.denominators.size()!=14) return 1;
  auto bubble=feynman::example_family("bubble");
  ibp::Generator generator(ibp::PropagatorBasis(ibp::merge(ibp::quadratic_family(bubble.momenta,x),0,1,x)),
    Exact(field,"2-2*eps"));
  ibp::ExactReducer reducer(x,1000);
  ibp::for_each_seed(1,2,{2,2,100},[&](const ibp::Integral& seed) {
    for(auto& identity:generator.relations(seed))reducer.insert(std::move(identity));
  });
  auto equation=ibp::differential_system(generator,reducer,{{2,0}},0,x);
  if(!(equation.matrix[0][0]==Exact(field,"-(1+eps)*(1-2*x)/(1+x-x*x)")))return 1;
  auto graph=recursion::prepare(bubble,{});
  auto value=recursion::evaluate_leaf(graph.nodes[0],graph.dimension,2,0);
  if(!value.taylor_tail_certified || value.values.size()!=1)return 1;
  recursion::Evaluator evaluator(graph);auto evaluated=evaluator.evaluate(0);
  auto snapshot=evaluator.linear_expression();Jet::Ball::set_precision(384);
  auto restored=linear_boundary::materialize(snapshot,0);
  if(restored.low!=evaluated.low || !acb_overlaps(restored.values[0][0].raw(),evaluated.values[0][0].raw()))return 1;
  const auto directory=std::filesystem::temp_directory_path()/("de3-installed-cache-"+std::to_string(::getpid()));
  struct Cleanup {std::filesystem::path path;~Cleanup(){std::error_code error;std::filesystem::remove_all(path,error);}} cleanup{directory};
  artifacts::Store store(directory);
  AffineFrobeniusSeries::Matrix endpoint{{Exact(field,"eps/x")}};
  auto fresh=cached_affine::prepare(endpoint,0,1,4,{},store);
  auto hit=cached_affine::prepare(endpoint,0,1,4,{},store);
  if(fresh.cache_hit || !hit.cache_hit || fresh.content_id!=hit.content_id)return 1;
  const auto zero=x.constant(0),one=x.constant(1);
  const auto rows=exact_laurent_rows({{one}},zero,0);
  adjoint_checkpoint::Statistics ordinary_stats;
  adjoint_checkpoint::Options storage{directory/"ordinary",64*1024*1024,&ordinary_stats};
  const auto transported=adjoint_checkpoint::transport({{zero}},rows,{{one}},{zero,one},{},storage);
  const auto reused=adjoint_checkpoint::transport({{zero}},rows,{{one}},{zero,one},{},storage);
  if(ordinary_stats.completed_reused!=1 || !acb_equal(transported.coefficients[0][0][0].raw(),reused.coefficients[0][0][0].raw()) ||
      !acb_equal_si(reused.coefficients[0][0][0].raw(),2))return 1;
  recursion::NumericalOptions invalid;invalid.adjoint.continuation=std::make_shared<const AdjointContinuation>();
  bool rejected=false;try{recursion::Evaluator unsafe(graph,invalid);}catch(const std::invalid_argument&){rejected=true;}
  if(!rejected)return 1;
  std::cout<<"Installed C++ package: exact recurrence, Banana4 coordinates, certified scalar leaf, shared-source snapshot, verified endpoint cache and ordinary continuation passed\n";
}
