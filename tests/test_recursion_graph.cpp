#include "diffexp/recursion_graph.hpp"
#include <iostream>
#include <thread>
using namespace diffexp;
void require(bool test,const char* message){if(!test)throw std::runtime_error(message);}
int main(int argc,char** argv){try {
  {
    ExactField field({"x","eps"});
    std::vector<std::vector<Exact>> matrix{{Exact(field,"eps/(x-1/3+eps)")}};
    require(!recursion::regular_anchor(matrix,0,1,Rational("1/3")),"higher epsilon coefficients forbid a hidden singular anchor");
    require(recursion::regular_anchor(matrix,0,1,Rational("1/2")),"regular anchor remains available");
    auto chosen=recursion::select_regular_anchor(matrix,{},0,1,Rational("1/3"));
    require(chosen!=Rational("1/3")&&recursion::regular_anchor(matrix,0,1,chosen),"automatic anchor moves before child preparation");
    // A requested ladder keeps its other rungs free: 1/2 is skipped, and the
    // simplest regular rational within reach of 1/3 is 1/4.
    require(recursion::select_regular_anchor(matrix,{},0,1,Rational("1/3"),{Rational("1/3"),Rational("1/2")})==Rational("1/4"),
      "a singular requested anchor moves to the simplest regular neighbour");
    require(recursion::select_regular_anchor(matrix,{},0,1,Rational("2/5"),{Rational("2/5")})==Rational("2/5"),"a regular requested anchor is kept");
    // A second real pole at 7/24 separates 1/3 from 1/4: the replacement must
    // stay on the request's side of it.
    std::vector<std::vector<Exact>> fenced{{Exact(field,"eps/(x-1/3+eps)"),Exact(field,"1/(24*x-7)")}};
    const auto kept=recursion::select_regular_anchor(fenced,{},0,1,Rational("1/3"),{Rational("1/3"),Rational("1/2")});
    require(kept>Rational("7/24")&&kept!=Rational("1/3")&&recursion::regular_anchor(fenced,0,1,kept),
      "a replacement anchor does not cross another real singular point");
  }
  recursion::Options options;
  auto bubble=recursion::prepare(feynman::example_family("bubble"),{},options);
  require(bubble.nodes.size()==1&&bubble.nodes[0].scalar_leaf,"bubble reaches scalar leaf without a reduction provider");
  require(bubble.nodes[0].source_indices==std::vector<ibp::Integral>{{2,0}},"bubble scalar target follows actual merge coordinates");
  unsigned calls=0;
  auto native=[&](const ibp::PropagatorBasis& b,const Exact& d,const ExactField&,std::size_t p,
      const std::vector<ibp::Integral>& requested,const level::Options&) {
    ++calls;ibp::Generator generator(b,d);ibp::ExactReducer reducer(d,300);
    ibp::for_each_seed(b.physical_count,b.denominators.size(),{1,1,50},[&](const auto& seed){
      for(auto& relation:generator.relations(seed))reducer.insert(std::move(relation));
    });
    ibp::BasisReduction coordinates(reducer,{{1,0,0,0,0},{1,1,0,0,0},{1,1,-1,0,0}},d);
    level::Result result;result.ordered_basis=coordinates.ordered_basis();
    for(const auto& a:result.ordered_basis)result.matrix.push_back(coordinates.resolve(generator.derivative(a,p)));
    for(const auto& a:requested)result.target_rows.push_back(coordinates.resolve({{a,d.constant(1)}}));
    result.success=true;return result;
  };
  auto sunrise=recursion::prepare(feynman::example_family("sunrise"),{},options,native);
  require(calls==1&&sunrise.nodes.size()==2&&sunrise.nodes.back().scalar_leaf,"shared sunrise master vector descends once");
  const auto& outer=sunrise.nodes[0];const auto& inner=sunrise.nodes[1];
  require(inner.requested==outer.closure.ordered_basis,"child evaluates parent's exact ordered boundary vector");
  auto frozen=recursion::at_anchor(outer.merged,0,outer.anchor);ibp::PropagatorBasis frozen_basis(frozen);
  for(std::size_t i=0;i<inner.incoming.denominators.size();++i)
    require(pullback::equal(inner.incoming.denominators[i],frozen_basis.denominators[i]),"recursive numerator coordinate preservation");
  require(inner.operations[0].operation==feynman::Operation::UpperLimit,"single-denominator master is typed endpoint operation");
  require(inner.operations[2].source_row!=inner.operations[1].source_row,
    "recursive irreducible numerator remains distinct from its scalar integral");
  bool retained_numerator=false;
  for(const auto& [a,c]:inner.operations[2].source_row)
    for(std::size_t j=0;j<a.size();++j)retained_numerator=retained_numerator || a[j]<inner.operations[2].denominator_indices[j];
  require(retained_numerator,"recursive numerator factors lower source integral indices");
  bool rejected=false;options.anchors={Rational(0),Rational("1/3")};
  try{recursion::prepare(feynman::example_family("sunrise"),{},options,native);}catch(const std::invalid_argument&){rejected=true;}
  require(rejected,"singular matching anchor rejected before expensive work");options.anchors.clear();
  {
    // The total budget gates the START of a level. A level that completes
    // after the deadline (here the terminal, exactly zero one) is kept.
    auto hurried=options;hurried.total_timeout_seconds=1;
    auto slow_zero=[](const ibp::PropagatorBasis&,const Exact&,const ExactField&,std::size_t,
        const std::vector<ibp::Integral>& requested,const level::Options&) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1100));
      level::Result result;result.success=true;result.target_rows.assign(requested.size(),{});return result;
    };
    auto late=recursion::prepare(feynman::example_family("sunrise"),{},hurried,slow_zero);
    require(late.nodes.size()==1 && late.nodes[0].closure.ordered_basis.empty() && late.events.back().elapsed_seconds>1,
      "completed terminal level discarded for finishing after the total budget");
  }
  if(argc>1) {
    options.reduction.provider.executable=argv[1];options.reduction.provider.timeout_seconds=90;
    options.reduction.total_timeout_seconds=180;options.total_timeout_seconds=300;
    const std::string family=argc>2?argv[2]:"sunrise";
    auto graph=recursion::prepare(feynman::example_family(family),{},options,{},[](const auto& e){
      std::cout<<"level="<<e.depth<<" physical="<<e.physical_propagators<<" sources="<<e.sources
        <<" masters="<<e.masters<<" batches="<<e.provider_passes<<" seconds="<<e.elapsed_seconds<<std::endl;
    });
    require(graph.nodes.back().scalar_leaf,"real family descends to native scalar leaf");
  }
  std::cout<<"native recursive preparation passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
