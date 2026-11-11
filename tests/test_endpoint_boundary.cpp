// The endpoint boundary of a Feynman-trick level: the kinematic rule, the
// exact determination of the analytic solution by the endpoint values, the
// child family it prepares, and the agreement of an endpoint run with the
// interior-anchor run of the same family.
#include "diffexp/recursion_pipeline.hpp"
#include "diffexp/level_cache.hpp"
#include "diffexp/ibp_solver_provider.hpp"
#include <iostream>
#include <unistd.h>
using namespace diffexp;
namespace eb=endpoint_boundary;
static void check(bool test,const std::string& message){if(!test)throw std::runtime_error(message);}
int main(){try {
  // The kinematic rule on the kite: lines 0 and 1 carry the first loop
  // momentum, line 2 the second, and all are massive.
  {
    auto kite=feynman::example_family("kite");
    const bool euclidean=causal::euclidean(kite);
    check(euclidean,"the kite has positive Symanzik polynomials");
    std::vector<std::vector<Rational>> lines(5,std::vector<Rational>(5,Rational(0)));
    for(unsigned i=0;i<5;++i)lines[i][i]=Rational(1);
    check(eb::kinematics(kite,euclidean,lines,0,1,1).qualifies && eb::kinematics(kite,euclidean,lines,0,1,0).qualifies,
      "a massive pair on one loop momentum qualifies at both ends");
    check(!eb::kinematics(kite,false,lines,0,1,1).qualifies,"without proved positivity no endpoint qualifies");
    check(!eb::kinematics(kite,euclidean,lines,0,2,2).qualifies,"lines on different loop momenta do not qualify");
    auto massless=kite;massless.momenta.lines[1].mass_squared=Rational(0);
    check(!eb::kinematics(massless,euclidean,lines,0,1,1).qualifies,"a massless kept line does not qualify");
    check(eb::kinematics(massless,euclidean,lines,0,1,0).qualifies,"the massive line of the pair may be the one kept");
    // A line merged earlier from lines 0 and 1, now paired with line 4.
    auto merged=lines;merged[0]={Rational("1/3"),Rational("2/3"),Rational(0),Rational(0),Rational(0)};
    check(!eb::kinematics(kite,euclidean,merged,0,4,0).qualifies,"an earlier merge keeps its loop momentum");
    // Its mass is the weighted mass of its propagators.
    check(eb::kinematics(massless,euclidean,merged,0,1,0).qualifies,"a merged line with a massive part is massive");
  }
  // The equation half on systems with known solutions.
  ExactField field({"x","eps"});
  const auto matrix=[&](std::initializer_list<std::initializer_list<const char*>> rows) {
    eb::Matrix out;for(const auto& row:rows){out.emplace_back();for(const auto* entry:row)out.back().emplace_back(field,entry);}return out;
  };
  {
    // f'=0, g'=f+2*eps*g/x: g=c*x^(2 eps)+f*x/(1-2 eps). The analytic solution
    // is fixed by f(0).
    const auto system=matrix({{"0","0"},{"1","2*eps/x"}});
    const auto verdict=eb::examine(system,0,1,false);
    check(verdict.qualifies && verdict.analytic_dimension==1,"one analytic solution fixed by its endpoint value: "+verdict.reason);
    const auto fuchs=fuchsify::prepare(system,0);
    const auto series=AffineFrobeniusSeries::prepare(fuchs.matrix,0,1,6);
    const auto frame=series.project(fuchs.transform);
    const auto determination=eb::determine(series,frame,1,fuchsify::detail::pole_order(fuchs.inverse_transform,0));
    check(determination.determined && determination.columns.size()==1 && determination.analytic_dimension==1,"determination of the analytic column");
    const auto map=eb::taylor_map(frame,determination,Rational("1/4"));
    check(map[0][0]==Exact(field,"1") && map[1][0]==Exact(field,"1/(4*(1-2*eps))") && map[0][1].is_zero() && map[1][1].is_zero(),
      "the boundary map is the Taylor series of the analytic solution");
    // The same system seen from x=1.
    const auto mirrored=matrix({{"0","0"},{"1","2*eps/(x-1)"}});
    check(eb::examine(mirrored,0,1,true).qualifies && eb::examine(mirrored,0,1,false).qualifies,
      "the upper endpoint is examined in 1-x; a regular point qualifies trivially");
  }
  {
    // g=c*sqrt(x) is not analytic: its constant is zero like that of a
    // regulated power, unless the strict rule is asked for.
    const auto root=eb::examine(matrix({{"0","0"},{"0","1/(2*x)"}}),0,1,false);
    check(root.qualifies && root.analytic_dimension==1,"a fixed square root carries no analytic solution: "+root.reason);
    eb::Options strict;strict.fixed_noninteger_exponents=false;
    const auto refused=eb::examine(matrix({{"0","0"},{"0","1/(2*x)"}}),0,1,false,strict);
    check(!refused.qualifies && refused.reason.find("non-integer")!=std::string::npos,"the strict rule refuses a fixed square root: "+refused.reason);
    // f=c, g=c*log(x)+b: the analytic solutions are f=0, g=b.
    const auto logarithm=eb::examine(matrix({{"0","0"},{"1/x","0"}}),0,1,false);
    check(logarithm.qualifies && logarithm.analytic_dimension==1,"a logarithmic solution is excluded, the constant one kept: "+logarithm.reason);
    {
      const auto system=matrix({{"0","0"},{"1/x","0"}});
      const auto fuchs=fuchsify::prepare(system,0);
      const auto series=AffineFrobeniusSeries::prepare(fuchs.matrix,0,1,4);
      const auto frame=series.project(fuchs.transform);
      const auto map=eb::taylor_map(frame,eb::determine(series,frame,1,fuchsify::detail::pole_order(fuchs.inverse_transform,0)),Rational("1/8"));
      check(map[0][0].is_zero() && map[0][1].is_zero() && map[1][0].is_zero() && map[1][1]==Exact(field,"1"),
        "the boundary map of a logarithmic block keeps its analytic solution only");
    }
    // f=a/x+..., an analytic solution only for a=0: h'=-h/x, g'=h.
    const auto pole=eb::examine(matrix({{"-1/x","0"},{"1","0"}}),0,1,false);
    check(pole.qualifies && pole.analytic_dimension==1,"a solution with a pole is excluded: "+pole.reason);
    // g=c*x vanishes at the endpoint: its constant would need g'(0).
    const auto slope=eb::examine(matrix({{"0","0"},{"0","1/x"}}),0,1,false);
    check(!slope.qualifies && slope.reason.find("Taylor")!=std::string::npos,"a solution that vanishes at the endpoint is refused: "+slope.reason);
  }
  // The one-loop massive family: every level merges two massive lines on the
  // loop momentum, so every level matches at x=0.
  const auto cache=std::filesystem::temp_directory_path()/("diffexp-endpoint-boundary-test-"+std::to_string(::getpid()));
  struct Cleanup {std::filesystem::path path;~Cleanup(){std::error_code ignored;std::filesystem::remove_all(path,ignored);}} cleanup{cache};
  artifacts::Store store(cache);
  const auto provider=[&](const auto& basis,const auto& dimension,const auto& field,auto parameter,const auto& sources,const auto& budget) {
    std::shared_ptr<ibp_solver::Session> native;
    level::Provider fallback=[&](const auto& batch,const auto& limits) {
      if(!native)native=std::make_shared<ibp_solver::Session>(basis,dimension,field,ibp_solver::recursion_options(cache/"ibp-solver",parameter,600));
      return (*native)(batch,limits);
    };
    return cached_level::prepare(store,basis,dimension,field,parameter,sources,budget,fallback).result;
  };
  const auto family=feynman::example_family("pentagon_massive");
  recursion::Options anchored,automatic;automatic.boundary=recursion::BoundaryPolicy::automatic;
  check(anchored.boundary==recursion::BoundaryPolicy::anchor,"the library default keeps interior anchors");
  const auto reference=recursion::prepare(family,{},anchored,provider);
  const auto graph=recursion::prepare(family,{},automatic,provider);
  check(graph.nodes.size()==reference.nodes.size(),"the endpoint ladder has the same depth");
  for(const auto& node:reference.nodes)check(node.boundary==recursion::BoundarySite::anchor,"the anchor policy takes no endpoint");
  for(std::size_t depth=0;depth+1<graph.nodes.size();++depth) {
    const auto& node=graph.nodes[depth];const auto& child=graph.nodes[depth+1];
    check(node.boundary==recursion::BoundarySite::lower && node.anchor.is_zero(),"every level of the massive one-loop family matches at x=0");
    // The child is the incoming family without the left line of the pair:
    // the merged slot holds the right one, the others are unchanged.
    check(child.incoming.physical_count+1==node.incoming.physical_count,"the child has one propagator fewer");
    const unsigned slot=node.left-(node.left>node.right);
    for(std::size_t k=0;k<node.incoming.physical_count;++k) {
      if(k==node.left)continue;
      const auto position=k==node.right?slot:k-(k>node.right);
      check(pullback::equal(child.incoming.denominators[position],node.incoming.denominators[k]),"the child keeps the incoming propagators");
    }
    // What the level asks of it: its masters at the endpoint, which are
    // integrals of that family in the level's own numerator slots.
    check(child.requested==node.closure.ordered_basis,"the child is asked for the level's masters at the endpoint");
    for(std::size_t k=node.merged.physical_count;k<node.merged.denominators.size();++k)
      check(pullback::equal(child.incoming.denominators[k],recursion::at_anchor(node.merged,0,Rational(0)).auxiliary[k-node.merged.physical_count]),
        "the numerator slots of the level are those of the child");
  }
  {
    // A stated ladder of anchors is kept, and an endpoint that was asked for
    // but is not there is an error.
    auto stated=automatic;stated.anchors.assign(family.physical_count-1,Rational("1/3"));
    for(const auto& node:recursion::prepare(family,{},stated,provider).nodes)
      check(node.boundary==recursion::BoundarySite::anchor,"explicit anchors keep every level at its anchor");
    recursion::Options demanded;demanded.boundary=recursion::BoundaryPolicy::endpoint;
    bool refused=false;
    // The sunrise has no two lines on one loop momentum.
    try{(void)recursion::prepare(feynman::example_family("sunrise"),{},demanded,provider);}catch(const std::domain_error&){refused=true;}
    check(refused,"a demanded endpoint that no level offers is an error");
  }
  recursion::NumericalOptions numerical;numerical.automatic_basis_scan=false;
  recursion::Evaluator at_anchor(reference,numerical),at_endpoint(graph,numerical);
  const auto expected=at_anchor.evaluate(1),value=at_endpoint.evaluate(1);
  for(int order=std::min(expected.low,value.low);order<=1;++order) {
    using B=kernel::ComplexBall;
    const auto at=[&](const LaurentBoundary& result){return order<result.low?B(0):result.values.at(0).at(order-result.low);};
    const auto difference=NativeTailMagnitude::upper_abs(at(expected)-at(value)).approximate_upper();
    check(difference<1e-25,"endpoint and anchor runs agree at epsilon order "+std::to_string(order)+": difference "+std::to_string(difference));
  }
  bool refused=false;
  try{(void)at_endpoint.prepare_adjoint_stage(0,0);}catch(const std::domain_error&){refused=true;}
  check(refused,"the standalone stage export names an endpoint level instead of exporting an anchor");
  std::cout<<"endpoint boundary tests passed\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
