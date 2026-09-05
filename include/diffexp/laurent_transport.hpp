#pragma once
#include "diffexp/scalar_functional.hpp"
#include "diffexp/epsilon_gauge.hpp"

namespace diffexp {
struct LaurentBoundary {
  int low;
  Boundary values;
  bool taylor_tail_certified=false;
  int high() const {
    if(values.empty() || values[0].empty())throw std::invalid_argument("empty Laurent boundary");
    const auto n=values[0].size();
    for(const auto& row:values)if(row.size()!=n)throw std::invalid_argument("Laurent boundary window mismatch");
    if(n>1001 || low< -1000 || low>1000)throw std::invalid_argument("Laurent boundary window budget");
    return low+static_cast<int>(n)-1;
  }
};
struct BoundaryDemand : std::runtime_error {
  int required_high;
  explicit BoundaryDemand(int high,const std::string& reason):std::runtime_error(reason),required_high(high){}
};
inline std::pair<std::size_t,std::size_t> path_epsilon_variables(const Exact& sample) {
  const auto& names=sample.variables();auto x=std::find(names.begin(),names.end(),"x"),e=std::find(names.begin(),names.end(),"eps");
  if(x==names.end()||e==names.end())throw std::invalid_argument("Laurent transport needs x and eps variables");
  return {x-names.begin(),e-names.begin()};
}
inline std::vector<Exact> exact_point(const Exact& sample,std::size_t variable,const Exact& value) {
  std::vector<Exact> result;for(std::size_t i=0;i<sample.variable_count();++i)result.push_back(sample.variable(i));
  result.at(variable)=value;return result;
}

// Compile only the epsilon coefficients needed by this numerical window. The
// immutable rational matrix stays independent of Taylor/epsilon resource knobs.
inline std::vector<RationalLineEntry> mapped_rational_connection(
    const std::vector<std::vector<Exact>>& matrix,const Exact& from,const Exact& to,unsigned epsilon_top) {
  if(matrix.empty() || epsilon_top>1000)throw std::invalid_argument("mapped rational connection budget");
  auto [xi,ei]=path_epsilon_variables(from);auto width=to-from;
  const auto point=exact_point(from,xi,from+width*from.variable(xi));
  std::vector<RationalLineEntry> entries;
  for(unsigned i=0;i<matrix.size();++i) {
    if(matrix[i].size()!=matrix.size())throw std::invalid_argument("mapped connection must be square");
    for(unsigned j=0;j<matrix.size();++j) {
      const auto& coefficient=matrix[i][j];if(coefficient.is_zero())continue;
      auto valuation=*exact_epsilon_valuation(coefficient,ei);
      if(valuation<0)throw std::domain_error("ordinary connection requires an epsilon gauge before transport");
      auto expanded=feynman::scalar_functional_detail::epsilon_series(coefficient,ei,epsilon_top);
      for(unsigned e=0;e<expanded.size();++e)if(!expanded[e].is_zero())
        entries.push_back({i,j,e,width*expanded[e].substitute(point)});
    }
  }
  return entries;
}
inline LaurentBoundary transport_laurent(const std::vector<std::vector<Exact>>& matrix,
    LaurentBoundary initial,const std::vector<Exact>& vertices,unsigned taylor_order=80) {
  const auto high=initial.high();
  if(matrix.size()!=initial.values.size() || vertices.empty())throw std::invalid_argument("Laurent transport path or dimension");
  for(unsigned i=0;i+1<vertices.size();++i) {
    if(vertices[i]==vertices[i+1])continue;
    auto entries=mapped_rational_connection(matrix,vertices[i],vertices[i+1],high-initial.low);
    initial.values=rational_line(entries,std::move(initial.values),taylor_order);
    initial.taylor_tail_certified=false;
  }
  return initial;
}

inline LaurentBoundary apply_rational_rows(const std::vector<std::vector<Exact>>& rows,
    const Exact& at,const LaurentBoundary& source,int output_high) {
  using B=Jet::Ball;
  const auto high=source.high();auto [xi,ei]=path_epsilon_variables(at);
  if(rows.empty() || output_high< -1000 || output_high>1000)throw std::invalid_argument("rational observable output window");
  const auto point=exact_point(at,xi,at);auto zero=at.constant(0);
  std::vector<std::vector<Exact>> evaluated=rows;
  long minimum=1001;
  for(auto& row:evaluated) {
    if(row.size()!=source.values.size())throw std::invalid_argument("rational observable dimension");
    for(auto& c:row){c=c.substitute(point);if(!c.is_zero())minimum=std::min(minimum,static_cast<long>(*exact_epsilon_valuation(c,ei)));}
  }
  if(minimum==1001)return {std::min(0,output_high),Boundary(rows.size(),std::vector<B>(output_high-std::min(0,output_high)+1,B(0))),source.taylor_tail_certified};
  if(minimum< -100 || minimum>100)throw std::invalid_argument("rational observable epsilon valuation budget");
  const auto needed=output_high-minimum;
  if(needed>high)throw BoundaryDemand(needed,"rational observable needs additional source epsilon coefficients");
  const int low=std::min(output_high,static_cast<int>(source.low+minimum));
  Boundary result(rows.size(),std::vector<B>(output_high-low+1,B(0)));
  for(unsigned i=0;i<rows.size();++i)for(unsigned j=0;j<source.values.size();++j) {
    const auto& c=evaluated[i][j];if(c.is_zero())continue;
    const auto valuation=*exact_epsilon_valuation(c,ei);
    const int depth=output_high-source.low-valuation;
    if(depth<0)continue;
    if(depth>1000)throw std::invalid_argument("rational observable coefficient-depth budget");
    auto regular=valuation<0?c*c.variable(ei).pow(-valuation):c/c.variable(ei).pow(valuation);
    Jet epsilon(0,depth+1,B::precision());if(depth)epsilon.set(1,B(1));
    auto coefficients=evaluate(data::Reader(regular.str()).read(),epsilon,{{c.variables()[ei],epsilon}});
    for(int k=low;k<=output_high;++k)for(int n=0;n<=depth;++n) {
      const int source_order=k-valuation-n;
      if(source_order>=source.low && source_order<=high)
        result[i][k-low]+=coefficients.at(n)*source.values[j][source_order-source.low];
    }
  }
  return {low,std::move(result),source.taylor_tail_certified};
}

struct RegularIntegrals {LaurentBoundary endpoint,integrals;};
inline RegularIntegrals integrate_regular_rows(const std::vector<std::vector<Exact>>& matrix,
    const std::vector<std::vector<Exact>>& rows,const LaurentBoundary& initial,
    const std::vector<Exact>& vertices,int output_high,unsigned taylor_order=80) {
  using B=Jet::Ball;
  const auto high=initial.high();const auto d=matrix.size();
  if(!d || initial.values.size()!=d || rows.empty() || vertices.empty())throw std::invalid_argument("regular integral shape/path");
  for(const auto& row:matrix)if(row.size()!=d)throw std::invalid_argument("regular integral connection must be square");
  const auto& sample=vertices[0];auto [xi,ei]=path_epsilon_variables(sample);auto zero=sample.constant(0);
  auto augmented=matrix;for(auto& row:augmented)row.resize(d+rows.size(),zero);
  auto state=initial;
  std::vector<long> shifts;long lowest=100;
  for(const auto& row:rows) {
    if(row.size()!=d)throw std::invalid_argument("regular integral row dimension");
    long valuation=1001;for(const auto& c:row)if(!c.is_zero())valuation=std::min(valuation,static_cast<long>(*exact_epsilon_valuation(c,ei)));
    if(valuation==1001)valuation=0;
    if(valuation< -100 || valuation>100)throw std::invalid_argument("regular integral epsilon pole budget");
    if(output_high-valuation>high)throw BoundaryDemand(output_high-valuation,"regular beta accumulator needs additional source epsilon coefficients");
    shifts.push_back(valuation);lowest=std::min(lowest,valuation);
    auto normalized=row;
    for(auto& c:normalized)c=valuation<0?c*sample.variable(ei).pow(-valuation):c/sample.variable(ei).pow(valuation);
    normalized.resize(d+rows.size(),zero);augmented.push_back(std::move(normalized));
    state.values.push_back(std::vector<B>(high-initial.low+1,B(0)));
  }
  state=transport_laurent(augmented,std::move(state),vertices,taylor_order);
  const int low=std::min(output_high,static_cast<int>(initial.low+lowest));
  LaurentBoundary result{low,Boundary(rows.size(),std::vector<B>(output_high-low+1,B(0))),false};
  for(unsigned i=0;i<rows.size();++i)for(int k=low;k<=output_high;++k) {
    const int source=k-shifts[i]-initial.low;
    if(source>=0 && static_cast<std::size_t>(source)<state.values[d+i].size())result.values[i][k-low]=state.values[d+i][source];
  }
  state.values.resize(d);return {std::move(state),std::move(result)};
}
} // namespace diffexp
