#pragma once
#include "diffexp/linear_boundary.hpp"

namespace diffexp::factored_transport {
using Expression=linear_boundary::Expression;
using B=Jet::Ball;
struct Options {
  AdjointOptions transport;
  unsigned max_augmented_dimension=512,max_leaf_width=64;
  // Optional orchestration adapter; the default mathematical path has no I/O.
  std::function<LaurentRows(const ExactEpsilonMatrix&,LaurentRows,const ExactEpsilonMatrix&,
      const std::vector<Exact>&,const AdjointOptions&)> transport_dispatch;
};
struct MapDemand : std::runtime_error {
  int required_high;
  explicit MapDemand(int high):std::runtime_error("factored transport requires additional initial map coefficients"),required_high(high){}
};
struct Result {
  Expression physical,integrated;
  int required_initial_high=0,observable_lower_bound=0;
  bool omitted_tail_certified=false;
};
// output_high is a coefficient-MAP order, not a materialized-source order.
// To materialize through T with source lower bound L, request maps through T-L.
// Evolve Y=M*s and Z'=B*M*s, Z(from)=0, without reading any source balls.
// The returned maps retain the very same immutable source. Exact epsilon gauges
// belong on A/M/B before this call; A must be epsilon regular. A negative B
// valuation is handled by an exact monomial rescaling of the accumulated map.
inline Result evolve(const ExactEpsilonMatrix& a,const Expression& initial,
    const ExactEpsilonMatrix& observable,const std::vector<Exact>& vertices,
    int output_high,const Options& options={}) {
  const auto d=a.size(),r=observable.size(),s=initial.transform.columns();
  if(!d || !r || !initial.leaf_source || initial.leaf_source->values.size()!=s ||
      initial.transform.coefficients.size()!=d || vertices.empty() ||
      d+r>options.max_augmented_dimension || s>options.max_leaf_width ||
      options.max_augmented_dimension>5000 || options.max_leaf_width>5000 ||
      output_high<initial.transform.low || output_high>1000)
    throw std::invalid_argument("factored transport dimensions, source or window budget");
  const auto [xi,ei]=path_epsilon_variables(vertices[0]);
  const auto zero=vertices[0].constant(0),eps=zero.variable(ei);
  int q=0;
  for(const auto& row:a) {
    if(row.size()!=d)throw std::invalid_argument("factored connection must be square");
    for(const auto& entry:row)if(!entry.is_zero() && *exact_epsilon_valuation(entry,ei)<0)
      throw std::domain_error("factored connection requires an exact epsilon gauge");
  }
  for(const auto& row:observable) {
    if(row.size()!=d)throw std::invalid_argument("factored observable width mismatch");
    for(const auto& entry:row)if(!entry.is_zero()) {
      const auto valuation=*exact_epsilon_valuation(entry,ei);
      if(valuation< -1000)throw std::length_error("factored observable epsilon pole budget");
      q=std::min(q,static_cast<int>(valuation));
    }
  }
  const long high=static_cast<long>(output_high)-q;
  const long integral_low=static_cast<long>(initial.transform.low)+q;
  if(high>1000 || integral_low< -1000)
    throw std::length_error("factored transport epsilon window budget");
  if(initial.transform.high<high)throw MapDemand(static_cast<int>(high));
  const auto augmented=d+r;
  // Using the existing shared-coefficient backend: g=M^T satisfies
  // g'=-g*(-A_aug^T). Each independent row is one leaf-source component.
  ExactEpsilonMatrix negative_transpose(augmented,std::vector<Exact>(augmented,zero));
  for(unsigned i=0;i<d;++i)for(unsigned j=0;j<d;++j)negative_transpose[j][i]=-a[i][j];
  const auto scale=eps.pow(-q);
  for(unsigned i=0;i<r;++i)for(unsigned j=0;j<d;++j)
    negative_transpose[j][d+i]=-scale*observable[i][j];
  const auto low=initial.transform.low;
  const auto width=static_cast<std::size_t>(high-low+1);
  if(s*augmented>20000000/width)
    throw std::length_error("factored persistent map storage budget");
  LaurentRows transposed{low,static_cast<int>(high),
      std::vector(s,std::vector(augmented,std::vector<B>(width,B(0))))};
  for(unsigned i=0;i<d;++i)for(unsigned j=0;j<s;++j)
    for(int k=low;k<=high;++k)transposed.coefficients[j][i][k-low]=initial.transform.coefficients[i][j][k-low];
  const auto dispatch=options.transport_dispatch?options.transport_dispatch:transport_adjoint_rows;
  auto transported=dispatch(negative_transpose,std::move(transposed),
      ExactEpsilonMatrix(s,std::vector<Exact>(augmented,zero)),vertices,options.transport);
  // The backend can conservatively lower its common bound to zero; retain the
  // sharper structural input lower bound because this augmented system is
  // homogeneous and epsilon regular, and its integrated initial map is zero.
  LaurentRows physical{low,output_high,std::vector(d,std::vector(s,std::vector<B>(output_high-low+1,B(0))))};
  LaurentRows integrated{static_cast<int>(integral_low),output_high,
      std::vector(r,std::vector(s,std::vector<B>(output_high-integral_low+1,B(0))))};
  for(unsigned i=0;i<d;++i)for(unsigned j=0;j<s;++j)for(int k=low;k<=output_high;++k)
    physical.coefficients[i][j][k-low]=transported.coefficients[j][i][k-transported.low];
  for(unsigned i=0;i<r;++i)for(unsigned j=0;j<s;++j)for(int k=integral_low;k<=output_high;++k)
    integrated.coefficients[i][j][k-integral_low]=transported.coefficients[j][d+i][k-q-transported.low];
  return {{std::move(physical),initial.leaf_source},{std::move(integrated),initial.leaf_source},
      static_cast<int>(high),q,false};
}
} // namespace diffexp::factored_transport
