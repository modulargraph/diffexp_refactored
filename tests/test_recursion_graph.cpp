#include "diffexp/recursion_graph.hpp"
#include <iostream>
using namespace diffexp;
void require(bool test,const char* message){if(!test)throw std::runtime_error(message);}
int main(int argc,char** argv){try {
  {
    ExactField field({"x","eps"});
    std::vector<std::vector<Exact>> matrix{{Exact(field,"eps/(x-1/3+eps)")}};
    require(!recursion::regular_anchor(matrix,0,1,Rational("1/3")),"higher epsilon coefficients forbid a hidden singular anchor");
    require(recursion::regular_anchor(matrix,0,1,Rational("1/2")),"regular anchor remains available");
    auto chosen=recursion::select_regular_anchor(matrix,{},0,1,Rational("1/3"),true);
    require(chosen!=Rational("1/3")&&recursion::regular_anchor(matrix,0,1,chosen),"automatic anchor moves before child preparation");
    bool rejected=false;
    try{recursion::select_regular_anchor(matrix,{},0,1,Rational("1/3"),false);}catch(const std::domain_error&){rejected=true;}
    require(rejected,"explicit singular anchor fails with its mathematical reason");
    auto henn=recursion::default_anchors("henn_double_pentagon_x0",7);
    require(henn.front()==Rational("1/5")&&henn.back()==Rational("4/5")&&std::set<Rational>(henn.begin(),henn.end()).size()==7,"Henn preserves its distinct prescribed anchors");
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
