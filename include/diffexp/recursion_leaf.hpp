#pragma once
#include "diffexp/recursion_graph.hpp"
#include "diffexp/gaussian_boundary.hpp"
#include "diffexp/scalar_functional.hpp"

namespace diffexp::recursion {
struct LeafBoundary {
  int epsilon_low;
  Boundary values;
  bool taylor_tail_certified=true;
};
inline feynman::LaurentValue evaluate_anchored_gaussian(const gaussian::Boundary& source,
    std::size_t epsilon,int d0,unsigned top,int rim=1) {
  using B=Jet::Ball;
  if(top>100 || !(source.dimension==source.dimension.constant(d0)-source.dimension.constant(2)*source.dimension.variable(epsilon)))
    throw std::invalid_argument("recursive Gaussian dimension or epsilon budget");
  if(source.multiplier.is_zero())return {0,std::vector<B>(top+1,B(0)),true};
  const auto factor=feynman::factor_origin(source.multiplier,epsilon);
  if(factor.power< -100 || factor.power>100)throw std::runtime_error("recursive Gaussian epsilon valuation budget");
  auto scalar=feynman::tadpole(source.U.rational(),source.F.rational(),source.power,source.loops,d0,
    top+static_cast<unsigned>(std::max(0L,-factor.power)),rim);
  const int low=scalar.epsilon_low+factor.power;
  if(low>static_cast<int>(top))return {0,std::vector<B>(top+1,B(0)),true};
  Jet eps(0,top-low+1,B::precision());if(eps.length()>1)eps.set(1,B(1));
  auto regular=evaluate(data::Reader(factor.regular.str()).read(),eps,
    {{source.dimension.variables()[epsilon],eps}});
  std::vector<B> coefficients(top-low+1,B(0));
  for(unsigned k=0;k<coefficients.size();++k)for(unsigned j=0;j<=k;++j)
    coefficients[k]+=regular.at(j)*scalar.coefficients.at(k-j);
  return {low,std::move(coefficients),true};
}

// A leaf endpoint with only one positive denominator is itself a Gaussian
// integral in the incoming coordinates. Evaluating it directly preserves the
// complete numerator and avoids separately taking singular coefficient limits.
inline feynman::LaurentValue leaf_endpoint(const Node& node,const ibp::Integral& target,
    const Exact& dimension,int d0,unsigned top,const gaussian::Budget& budget={},int f_rim=1) {
  using B=Jet::Ball;
  std::vector<unsigned> positive;for(unsigned i=0;i<node.incoming.physical_count;++i)if(target[i]>0)positive.push_back(i);
  if(positive.empty())return {0,std::vector<B>(top+1,B(0)),true};
  if(positive.size()!=1)throw std::logic_error("leaf endpoint has more than one positive denominator");
  ibp::Integral signs(target.size(),-1);signs[positive[0]]=1;
  const auto zero_sectors=fire::free_loop_sectors(node.incoming);
  if(std::find(zero_sectors.begin(),zero_sectors.end(),signs)!=zero_sectors.end())
    return {0,std::vector<B>(top+1,B(0)),true};
  ibp::QuadraticFamily raw{node.incoming.space,{node.incoming.denominators[positive[0]]},{}};
  ibp::Integral reordered{target[positive[0]]};
  for(unsigned i=0;i<target.size();++i)if(i!=positive[0]) {
    raw.auxiliary.push_back(node.incoming.denominators[i]);reordered.push_back(target[i]);
  }
  auto result=gaussian::reduce(ibp::PropagatorBasis(std::move(raw)),reordered,dimension,budget);
  const auto& names=dimension.variables();auto e=std::find(names.begin(),names.end(),"eps");
  if(e==names.end())throw std::invalid_argument("recursive leaf requires epsilon variable eps");
  return evaluate_anchored_gaussian(result,e-names.begin(),d0,top,f_rim);
}

inline LeafBoundary evaluate_leaf(const Node& node,const Exact& dimension,int d0,unsigned top,
    const feynman::CertifiedDeepestBetaOptions& options={},const gaussian::Budget& budget={}) {
  using B=Jet::Ball;
  if(!node.scalar_leaf || node.merged.physical_count!=1 || node.incoming.physical_count!=2 ||
     node.requested.size()!=node.operations.size() || top>100)
    throw std::invalid_argument("recursive scalar leaf contract");
  struct PrecisionScope {slong previous;explicit PrecisionScope(slong bits):previous(B::precision()){B::set_precision(bits);}~PrecisionScope(){B::set_precision(previous);}} precision(options.working_bits);
  const auto& names=dimension.variables();auto xi=std::find(names.begin(),names.end(),"x");
  if(xi==names.end())throw std::invalid_argument("recursive leaf requires path variable x");
  const auto x=dimension.variable(xi-names.begin());
  std::map<ibp::Integral,gaussian::Boundary,ibp::IntegralOrder> sources;
  std::vector<feynman::LaurentValue> values;
  int low=0;
  for(std::size_t i=0;i<node.operations.size();++i) {
    const auto& operation=node.operations[i];
    if(operation.operation!=feynman::Operation::BetaIntegral) {
      values.push_back(leaf_endpoint(node,node.requested[i],dimension,d0,top,budget,options.f_rim));
    } else {
      auto weight=x.constant(0);std::optional<feynman::Symanzik> geometry;
      for(const auto& [a,c]:operation.source_row) {
        auto found=sources.find(a);
        if(found==sources.end())found=sources.emplace(a,gaussian::reduce(node.merged,a,dimension,budget)).first;
        const auto& source=found->second;if(source.multiplier.is_zero())continue;
        if(!geometry)geometry.emplace(feynman::Symanzik{source.U,source.F});
        else if(!(geometry->U==source.U) || !(geometry->F==source.F))throw std::logic_error("Gaussian leaf geometry depends on integral indices");
        // J_v / J_1 = product_{j=1}^{v-1}(j-L*d/2)/(j*M).
        // Combine all numerator/power reductions EXACTLY on this common gamma
        // source before either endpoint deflation or ball evaluation.
        auto ratio=x.constant(1),mass=source.F/source.U;
        for(unsigned j=1;j<source.power;++j)
          ratio=ratio*(x.constant(j)-x.constant(source.loops)*dimension/x.constant(2))/(x.constant(j)*mass);
        weight=weight+c*source.multiplier*ratio;
      }
      if(weight.is_zero())values.push_back({0,std::vector<B>(top+1,B(0)),true});
      else {
        weight=weight*x.constant(operation.normalization)*x.pow(operation.left_power)*(x.constant(1)-x).pow(operation.right_power);
        values.push_back(feynman::certified_scalar_functional(*geometry,node.merged.space.loops,d0,1,weight,top,options));
      }
    }
    low=std::min(low,values.back().epsilon_low);
  }
  Boundary boundary(values.size(),std::vector<B>(top-low+1,B(0)));
  for(unsigned i=0;i<values.size();++i)for(unsigned k=0;k<values[i].coefficients.size();++k)
    boundary[i][values[i].epsilon_low+k-low]=values[i].coefficients[k];
  if(!native_enclosure_meets_tolerance(boundary,"1/1"+std::string(options.requested_digits,'0')))
    throw std::runtime_error("recursive Gaussian leaf enclosure does not meet requested absolute digits");
  return {low,std::move(boundary),true};
}
} // namespace diffexp::recursion
