#pragma once
// Coefficient contracts for a directed linear computation. Unknown high
// coefficients are represented by bounds, never promoted to exact zeros.
#include "diffexp/linear_boundary.hpp"
#include <optional>

namespace diffexp::epsilon_demand {
using High = std::optional<int>;
using Highs = std::vector<High>;
using Matrix = ExactEpsilonMatrix;
inline int checked(long value) {
  if(value < -1000 || value > 1000)throw std::length_error("component epsilon demand budget");
  return static_cast<int>(value);
}
inline void raise(High& destination,int value) {
  if(!destination || *destination<value)destination=value;
}
inline Highs merge(Highs a,const Highs& b) {
  if(a.size()!=b.size())throw std::invalid_argument("component demand merge shape");
  for(std::size_t i=0;i<a.size();++i)if(b[i])raise(a[i],*b[i]);
  return a;
}
inline int maximum(const Highs& values,int empty=0) {
  High result;for(const auto& value:values)if(value)raise(result,*value);
  return result.value_or(empty);
}
// A_ij couples source j to destination i. This is backward reachability with
// epsilon valuations as weights, not merely an unweighted component graph.
inline Highs close(const Matrix& a,Highs highs,std::size_t ei) {
  const auto n=a.size();if(highs.size()!=n)throw std::invalid_argument("component demand shape");
  for(const auto& row:a)if(row.size()!=n)throw std::invalid_argument("component connection shape");
  for(std::size_t pass=0;pass<n;++pass) {
    bool changed=false;
    for(std::size_t i=0;i<n;++i)if(highs[i])for(std::size_t j=0;j<n;++j)
      if(!a[i][j].is_zero()) {
        const int need=checked(static_cast<long>(*highs[i])-*exact_epsilon_valuation(a[i][j],ei));
        if(!highs[j] || *highs[j]<need){highs[j]=need;changed=true;}
      }
    if(!changed)return highs;
    if(pass+1==n) {
      for(std::size_t i=0;i<n;++i)if(highs[i])for(std::size_t j=0;j<n;++j)
        if(!a[i][j].is_zero() && (!highs[j] || *highs[j]<*highs[i]-*exact_epsilon_valuation(a[i][j],ei)))
          throw std::domain_error("negative epsilon demand cycle requires another basis");
    }
  }
  return highs;
}
inline Highs pullback(const Matrix& map,const Highs& outputs,std::size_t ei) {
  if(map.size()!=outputs.size() || map.empty())throw std::invalid_argument("component exact map shape");
  Highs result(map.front().size());
  for(std::size_t i=0;i<map.size();++i) {
    if(map[i].size()!=result.size())throw std::invalid_argument("component exact map width");
    if(outputs[i])for(std::size_t j=0;j<result.size();++j)if(!map[i][j].is_zero())
      raise(result[j],checked(static_cast<long>(*outputs[i])-*exact_epsilon_valuation(map[i][j],ei)));
  }
  return result;
}
// The unknown tail starts at high+1 even if every retained coefficient is zero.
// Only exactly zero balls may improve the bound; zero-containing balls may not.
inline Highs pullback(const LaurentRows& map,const Highs& outputs) {
  if(map.coefficients.size()!=outputs.size())throw std::invalid_argument("component retained map shape");
  Highs result(map.columns());
  for(std::size_t i=0;i<outputs.size();++i)if(outputs[i])
    for(std::size_t j=0;j<result.size();++j) {
      int valuation=map.high+1;
      for(int k=map.low;k<=map.high;++k)if(!map.coefficients[i][j][k-map.low].is_zero()){valuation=k;break;}
      raise(result[j],checked(static_cast<long>(*outputs[i])-valuation));
    }
  return result;
}
struct Expression {
  linear_boundary::Expression value;
  Highs high;
};
inline Expression rectangular(const linear_boundary::Expression& value) {
  return {value,Highs(value.transform.coefficients.size(),value.transform.high)};
}
inline void validate(const Expression& expression) {
  linear_boundary::detail::validate(expression.value.transform,{});
  if(expression.high.size()!=expression.value.transform.coefficients.size())throw std::invalid_argument("component expression shape");
  for(const auto& h:expression.high)if(h && *h>expression.value.transform.high)
    throw std::invalid_argument("component known high exceeds stored range");
}
// Compute a complete output jet from individually truncated inputs. Every
// omitted input coefficient is checked against its actual consumer first.
inline linear_boundary::Expression compose(const LaurentRows& outer,const Expression& inner,int high) {
  validate(inner);const auto& in=inner.value.transform;
  linear_boundary::detail::validate(outer,{});
  if(outer.columns()!=inner.high.size())throw std::invalid_argument("component composition shape");
  const int outer_need=checked(static_cast<long>(high)-in.low);
  if(outer.high<outer_need)throw linear_boundary::CompositionDemand(outer_need,in.high);
  auto demand=pullback(outer,Highs(outer.coefficients.size(),high));
  for(std::size_t j=0;j<demand.size();++j)if(demand[j] && *demand[j]>=in.low &&
      (!inner.high[j] || *inner.high[j]<*demand[j]))
    throw linear_boundary::CompositionDemand(outer_need,*demand[j]);
  const int low=std::min(high,checked(static_cast<long>(outer.low)+in.low));
  LaurentRows result{low,high,std::vector(outer.coefficients.size(),std::vector(in.columns(),std::vector<Jet::Ball>(high-low+1,Jet::Ball(0))))};
  std::size_t work=0;
  for(std::size_t i=0;i<outer.coefficients.size();++i)for(std::size_t j=0;j<inner.high.size();++j)
    for(int p=outer.low;p<=outer_need;++p)if(!outer.coefficients[i][j][p-outer.low].is_zero())
      for(int q=in.low;q<=high-p;++q)for(std::size_t s=0;s<in.columns();++s) {
        if(++work>50000000)throw std::length_error("component composition operation budget");
        result.coefficients[i][s][p+q-low]+=outer.coefficients[i][j][p-outer.low]*in.coefficients[j][s][q-in.low];
      }
  return {std::move(result),inner.value.leaf_source};
}
inline std::size_t coefficient_slots(const Highs& highs,int low,std::size_t columns=1) {
  std::size_t count=0;for(const auto& h:highs)if(h && *h>=low)count+=static_cast<std::size_t>(*h-low+1)*columns;
  return count;
}
inline Expression compose(const LaurentRows& outer,const Expression& inner,const Highs& highs) {
  if(highs.size()!=outer.coefficients.size())throw std::invalid_argument("component output request shape");
  const int top=maximum(highs,inner.value.transform.low);
  int low=std::min(top,checked(static_cast<long>(outer.low)+inner.value.transform.low));
  LaurentRows result{low,top,std::vector(highs.size(),std::vector(inner.value.transform.columns(),std::vector<Jet::Ball>(top-low+1,Jet::Ball(0))))};
  for(std::size_t i=0;i<highs.size();++i)if(highs[i]) {
    LaurentRows row{outer.low,outer.high,{outer.coefficients[i]}};
    auto value=compose(row,inner,*highs[i]);
    for(std::size_t j=0;j<result.columns();++j)for(int k=value.transform.low;k<=*highs[i];++k)
      result.coefficients[i][j][k-low]=value.transform.coefficients[0][j][k-value.transform.low];
  }
  return {{std::move(result),inner.value.leaf_source},highs};
}
} // namespace diffexp::epsilon_demand
