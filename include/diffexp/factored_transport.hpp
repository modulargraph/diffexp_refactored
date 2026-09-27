#pragma once
#include "diffexp/linear_boundary.hpp"
#include "diffexp/epsilon_demand.hpp"

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
// Requested highs are coefficient-MAP orders, not materialized-source orders.
// To materialize through T with source lower bound L, request maps through T-L.
// Evolve Y=M*s and Z'=B*M*s, Z(from)=0, without reading any source balls.
// The returned maps retain the very same immutable source. Exact epsilon gauges
// belong on A/M/B before this call; A must be epsilon regular. A negative B
// valuation is handled by an exact monomial rescaling of the accumulated map.
// Options are explicit here so legacy evolve(..., high, {}) stays unambiguous.
inline Result evolve(const ExactEpsilonMatrix& a,const Expression& initial,
    const ExactEpsilonMatrix& observable,const std::vector<Exact>& vertices,
    int physical_high,int integrated_high,const Options& options) {
  const auto d=a.size(),r=observable.size(),s=initial.transform.columns();
  if(!d || !r || !initial.leaf_source || initial.leaf_source->values.size()!=s ||
      initial.transform.coefficients.size()!=d || vertices.empty() ||
      d+r>options.max_augmented_dimension || s>options.max_leaf_width ||
      options.max_augmented_dimension>5000 || options.max_leaf_width>5000 ||
      physical_high<initial.transform.low || physical_high>1000 || integrated_high>1000)
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
  // Endpoint composition and observable integration are independent consumers:
  // the observable pole raises only the integrated map requirement.
  const long high=std::max(static_cast<long>(physical_high),static_cast<long>(integrated_high)-q);
  const long integral_low=static_cast<long>(initial.transform.low)+q;
  if(integrated_high<integral_low)
    throw std::invalid_argument("factored integrated map window is below its lower bound");
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
  LaurentRows physical{low,physical_high,std::vector(d,std::vector(s,std::vector<B>(physical_high-low+1,B(0))))};
  LaurentRows integrated{static_cast<int>(integral_low),integrated_high,
      std::vector(r,std::vector(s,std::vector<B>(integrated_high-integral_low+1,B(0))))};
  for(unsigned i=0;i<d;++i)for(unsigned j=0;j<s;++j)for(int k=low;k<=physical_high;++k)
    physical.coefficients[i][j][k-low]=transported.coefficients[j][i][k-transported.low];
  for(unsigned i=0;i<r;++i)for(unsigned j=0;j<s;++j)for(int k=integral_low;k<=integrated_high;++k)
    integrated.coefficients[i][j][k-integral_low]=transported.coefficients[j][d+i][k-q-transported.low];
  return {{std::move(physical),initial.leaf_source},{std::move(integrated),initial.leaf_source},
      static_cast<int>(high),q,false};
}
// Compatibility entry point: retain both returned maps through one common high.
inline Result evolve(const ExactEpsilonMatrix& a,const Expression& initial,
    const ExactEpsilonMatrix& observable,const std::vector<Exact>& vertices,
    int output_high,const Options& options={}) {
  return evolve(a,initial,observable,vertices,output_high,output_high,options);
}

struct ComponentResult {
  epsilon_demand::Expression physical,integrated;
  epsilon_demand::Highs initial_high;
  std::size_t retained_slots=0,rectangular_slots=0,active_states=0;
};
// Compile a nonrectangular epsilon schedule into exact diagonal coordinates.
// The closure inequality h_j >= h_i-val(A_ij) makes the transformed connection
// epsilon-regular. Unknown coefficients are shifted ABOVE the retained top;
// the explicit returned contracts prevent reading them as zeros on extraction.
inline ComponentResult evolve_components(const ExactEpsilonMatrix& a,
    const epsilon_demand::Expression& initial,const ExactEpsilonMatrix& observable,
    const std::vector<Exact>& vertices,const epsilon_demand::Highs& physical_high,
    const epsilon_demand::Highs& integrated_high,const Options& options={}) {
  namespace ed=epsilon_demand;
  ed::validate(initial);
  const auto d=a.size(),r=observable.size(),s=initial.value.transform.columns();
  if(!d || vertices.empty() || physical_high.size()!=d || integrated_high.size()!=r ||
     initial.high.size()!=d || !initial.value.leaf_source || initial.value.leaf_source->values.size()!=s ||
     d+r>options.max_augmented_dimension || s>options.max_leaf_width)
    throw std::invalid_argument("component factored dimensions");
  const auto [xi,ei]=path_epsilon_variables(vertices.front());
  const auto z=vertices.front().constant(0),eps=z.variable(ei);
  ExactEpsilonMatrix forward(d+r,std::vector<Exact>(d+r,z));
  for(std::size_t i=0;i<d;++i) {
    if(a[i].size()!=d)throw std::invalid_argument("component connection width");
    for(std::size_t j=0;j<d;++j)forward[i][j]=a[i][j];
  }
  for(std::size_t i=0;i<r;++i) {
    if(observable[i].size()!=d)throw std::invalid_argument("component observable width");
    for(std::size_t j=0;j<d;++j)forward[d+i][j]=observable[i][j];
  }
  auto needs=physical_high;needs.insert(needs.end(),integrated_high.begin(),integrated_high.end());
  needs=ed::close(forward,std::move(needs),ei);
  const int top=ed::maximum(needs,initial.value.transform.low);
  const int input_low=initial.value.transform.low;
  ed::Highs input_needs(needs.begin(),needs.begin()+d);
  for(std::size_t i=0;i<d;++i)if(needs[i] && *needs[i]>=input_low &&
      (!initial.high[i] || *initial.high[i]<*needs[i]))throw MapDemand(*needs[i]);
  std::vector<std::size_t> active;std::vector<int> shifts(d+r);
  for(std::size_t i=0;i<needs.size();++i)if(needs[i]) {
    active.push_back(i);shifts[i]=ed::checked(static_cast<long>(top)-*needs[i]);
  }
  // Preserve certified lower zeros per input row. A ball containing zero is
  // still nonzero here; all incoming rounding uncertainty remains enclosed.
  int low=top;
  for(std::size_t i=0;i<d;++i)if(needs[i])for(std::size_t j=0;j<s;++j)
    for(int k=input_low;k<=*needs[i];++k)
      if(!initial.value.transform.coefficients[i][j][k-input_low].is_zero()) {
        low=std::min(low,ed::checked(static_cast<long>(k)+shifts[i]));break;
      }
  const auto width=static_cast<std::size_t>(top-low+1);
  linear_boundary::detail::product({s,active.size(),width},20000000);
  LaurentRows packed{low,top,std::vector(s,std::vector(active.size(),std::vector<B>(width,B(0))))};
  ExactEpsilonMatrix negative_transpose(active.size(),std::vector<Exact>(active.size(),z));
  for(std::size_t i=0;i<active.size();++i) {
    const auto original=active[i];
    if(original<d)for(std::size_t j=0;j<s;++j)for(int k=low;k<=top;++k) {
      const int original_order=k-shifts[original];
      if(original_order>=input_low && original_order<=*needs[original])
        packed.coefficients[j][i][k-low]=initial.value.transform.coefficients[original][j][original_order-input_low];
    }
    for(std::size_t j=0;j<active.size();++j) {
      const auto source=active[j];
      if(forward[original][source].is_zero())continue;
      auto coefficient=epsilon_gauge_detail::multiply_power(forward[original][source],ei,
        shifts[original]-shifts[source]);
      if(*exact_epsilon_valuation(coefficient,ei)<0)throw std::logic_error("component closure did not regularize connection");
      negative_transpose[j][i]=-coefficient;
    }
  }
  LaurentRows transported=packed;
  if(!active.empty()) {
    const auto dispatch=options.transport_dispatch?options.transport_dispatch:transport_adjoint_rows;
    transported=dispatch(negative_transpose,std::move(packed),
      ExactEpsilonMatrix(s,std::vector<Exact>(active.size(),z)),vertices,options.transport);
  }
  const auto extract=[&](std::size_t begin,const ed::Highs& requested) {
    const int high=ed::maximum(requested,input_low);
    int out_low=std::min(input_low,high);
    for(std::size_t i=0;i<requested.size();++i)if(requested[i])
      out_low=std::min(out_low,ed::checked(static_cast<long>(low)-shifts[begin+i]));
    LaurentRows result{out_low,high,std::vector(requested.size(),std::vector(s,std::vector<B>(high-out_low+1,B(0))))};
    for(std::size_t aidx=0;aidx<active.size();++aidx) {
      const auto original=active[aidx];
      if(original<begin || original>=begin+requested.size() || !requested[original-begin])continue;
      for(int k=out_low;k<=*requested[original-begin];++k) {
        const long packed_order=static_cast<long>(k)+shifts[original];
        if(packed_order<low)continue;
        if(packed_order>transported.high)throw std::logic_error("component extraction exceeds computed jet");
        for(std::size_t j=0;j<s;++j)
          result.coefficients[original-begin][j][k-out_low]=transported.coefficients[j][aidx][packed_order-transported.low];
      }
    }
    return ed::Expression{{std::move(result),initial.value.leaf_source},requested};
  };
  // Charge the actual packed storage, including shifted negative-order
  // accumulators. Logical output widths would undercount these registers.
  const auto slots=s*active.size()*width;
  return {extract(0,physical_high),extract(d,integrated_high),std::move(input_needs),
    slots,s*(d+r)*static_cast<std::size_t>(std::max(0,top-input_low+1)),active.size()};
}
} // namespace diffexp::factored_transport
