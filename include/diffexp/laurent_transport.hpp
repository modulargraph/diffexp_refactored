#pragma once
#include "diffexp/scalar_functional.hpp"
#include "diffexp/epsilon_gauge.hpp"

namespace diffexp {
struct LaurentBoundary {
  int low;
  Boundary values;
  bool taylor_tail_certified=false;
  // By row and coefficient, like `values`: the estimated size of the error
  // that the transports behind a value committed (the Taylor tails they
  // admit and their rounding, AdjointOptions::error). It is NOT in the radii
  // of `values`, which remain those of the retained polynomials, exactly
  // what the linear expression of the result materializes to; a caller that
  // reports the values widens them (widen_by_error). Empty: not estimated.
  std::vector<std::vector<double>> error={};
  // Adds `factor` times the estimated error to the radius of every value
  // that is not an exact zero. The estimate is one realisation of an error
  // whose contributions can cancel: hence the factor.
  void widen_by_error(double factor=8) {
    mag_t size;mag_init(size);
    for(std::size_t i=0;i<values.size() && i<error.size();++i)for(std::size_t k=0;k<values[i].size() && k<error[i].size();++k) {
      if(!(error[i][k]>0) || values[i][k].is_zero())continue;
      mag_set_d(size,factor*error[i][k]);
      acb_add_error_mag(values[i][k].raw(),size);
    }
    mag_clear(size);
  }
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
} // namespace diffexp
