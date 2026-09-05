#pragma once
#include "diffexp/exact.hpp"
#include <cstdint>
#include <limits>
#include <optional>

namespace diffexp {
using ExactEpsilonMatrix=std::vector<std::vector<Exact>>;
struct EpsilonGaugeOptions {
  std::size_t max_dimension=2048;
  std::size_t max_entries=1000000;
  std::uint64_t max_relaxations=20000000;
  std::int64_t max_abs_shift=1000000;
  unsigned long max_valuation_degree=1000000;
};
class UnsupportedEpsilonGauge:public std::domain_error {
 public:
  explicit UnsupportedEpsilonGauge(const std::string& message):std::domain_error(message){}
};

// nullopt is the valuation of an identically zero rational function (+infinity).
// This is the exact generic valuation in epsilon over the remaining rational
// function field, not a valuation after any numerical parameter specialization.
inline std::optional<std::int64_t> exact_epsilon_valuation(const Exact& value,
    std::size_t epsilon_variable,unsigned long max_degree=1000000) {
  if(epsilon_variable>=value.variable_count())throw std::out_of_range("epsilon valuation variable index");
  if(value.is_zero())return std::nullopt;
  auto minimum=[&](const std::vector<Exact::Term>& terms) {
    auto degree=std::numeric_limits<unsigned long>::max();
    for(const auto& term:terms)degree=std::min(degree,term.powers.at(epsilon_variable));
    if(degree>max_degree || degree>static_cast<unsigned long>(std::numeric_limits<std::int64_t>::max()))
      throw std::overflow_error("epsilon valuation exceeds the finite degree bound");
    return static_cast<std::int64_t>(degree);
  };
  return minimum(value.numerator_terms())-minimum(value.denominator_terms());
}

namespace epsilon_gauge_detail {
inline std::int64_t add(std::int64_t a,std::int64_t b) {
  if((b>0 && a>std::numeric_limits<std::int64_t>::max()-b) ||
      (b<0 && a<std::numeric_limits<std::int64_t>::min()-b))
    throw std::overflow_error("epsilon index arithmetic overflow");
  return a+b;
}
inline std::int64_t subtract(std::int64_t a,std::int64_t b) {
  if(b==std::numeric_limits<std::int64_t>::min()) {
    if(a>=0)throw std::overflow_error("epsilon index arithmetic overflow");
    return std::numeric_limits<std::int64_t>::max()+(a+1);
  }
  return add(a,-b);
}
inline Exact multiply_power(const Exact& value,std::size_t variable,std::int64_t exponent) {
  if(exponent==std::numeric_limits<std::int64_t>::min())throw std::overflow_error("epsilon gauge exponent overflow");
  auto power=static_cast<std::uint64_t>(exponent<0?-exponent:exponent);
  if(power>std::numeric_limits<unsigned long>::max())throw std::overflow_error("epsilon gauge exponent does not fit exact backend");
  auto monomial=value.variable(variable).pow(static_cast<unsigned long>(power));
  return exponent<0?value/monomial:value*monomial;
}
} // namespace epsilon_gauge_detail

struct EpsilonGaugeResult {
  // Y_i=epsilon^(shifts[i])*Z_i. The original row/column ordering is retained.
  ExactEpsilonMatrix matrix;
  std::vector<std::int64_t> shifts;
  std::size_t epsilon_variable=0;
};

// Exact diagonal, x-independent epsilon shearing. A'_ij has valuation
// val(A_ij)+s_j-s_i, so holomorphy requires s_i-s_j<=val(A_ij).
// All vertices start at distance zero (an implicit super-source); shortest
// path constraints give a feasible shear unless a negative cycle exists.
// A negative cycle excludes diagonal shearing, not every possible basis change.
// In particular exp(x/epsilon) cannot be represented by this holomorphic
// ordinary transport contract and is rejected explicitly.
inline EpsilonGaugeResult epsilon_diagonal_gauge(const ExactEpsilonMatrix& matrix,
    std::size_t epsilon_variable,const EpsilonGaugeOptions& options={}) {
  const auto n=matrix.size();
  if(!n || !options.max_dimension || !options.max_entries || !options.max_relaxations ||
      options.max_abs_shift<0 || !options.max_valuation_degree || n>options.max_dimension ||
      n>options.max_entries/n)
    throw std::invalid_argument("epsilon gauge matrix/options exceed finite shape budget");
  const auto& sample=matrix[0].empty()?throw std::invalid_argument("empty epsilon gauge matrix row"):matrix[0][0];
  if(epsilon_variable>=sample.variable_count())throw std::out_of_range("epsilon gauge variable index");
  struct Constraint {std::size_t from,to;std::int64_t weight;};
  std::vector<Constraint> constraints;
  for(std::size_t i=0;i<n;++i) {
    if(matrix[i].size()!=n)throw std::invalid_argument("epsilon gauge requires a square matrix");
    for(std::size_t j=0;j<n;++j) {
      (void)(sample.constant(0)+matrix[i][j]); // Enforce one owned exact field and ordering.
      auto valuation=exact_epsilon_valuation(matrix[i][j],epsilon_variable,options.max_valuation_degree);
      if(valuation)constraints.push_back({j,i,*valuation});
    }
  }
  std::vector<std::int64_t> distance(n,0);std::uint64_t relaxations=0;
  for(std::size_t pass=0;pass<n;++pass) {
    bool changed=false;
    for(const auto& edge:constraints) {
      if(relaxations++>=options.max_relaxations)throw std::runtime_error("epsilon gauge finite relaxation budget exhausted");
      auto proposed=epsilon_gauge_detail::add(distance[edge.from],edge.weight);
      if(proposed<distance[edge.to]) {
        distance[edge.to]=proposed;changed=true;
        if(pass+1==n)
          throw UnsupportedEpsilonGauge("negative epsilon-valuation cycle: no diagonal epsilon shear makes this connection holomorphic; a more general basis or essential-singularity treatment is required");
      }
    }
    if(!changed)break;
  }
  const auto minimum=*std::min_element(distance.begin(),distance.end());
  EpsilonGaugeResult result{matrix,{},epsilon_variable};result.shifts.reserve(n);
  for(auto value:distance) {
    auto shift=epsilon_gauge_detail::subtract(value,minimum);
    if(shift>options.max_abs_shift)throw std::overflow_error("epsilon gauge normalized shift exceeds finite index bound");
    result.shifts.push_back(shift);
  }
  for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j) {
    auto power=epsilon_gauge_detail::subtract(result.shifts[j],result.shifts[i]);
    result.matrix[i][j]=epsilon_gauge_detail::multiply_power(matrix[i][j],epsilon_variable,power);
    auto valuation=exact_epsilon_valuation(result.matrix[i][j],epsilon_variable,
      std::numeric_limits<unsigned long>::max());
    if(valuation && *valuation<0)throw std::logic_error("exact epsilon shear verification failed");
  }
  return result;
}

struct InclusiveEpsilonDemand {std::int64_t low,high;};
// Requests actual coefficients at shifted indices; it never pads missing rows.
inline InclusiveEpsilonDemand gauged_epsilon_demand(const EpsilonGaugeResult& gauge,
    std::size_t row,const InclusiveEpsilonDemand& physical) {
  if(row>=gauge.shifts.size())throw std::out_of_range("epsilon demand row index");
  if(physical.low>physical.high)throw std::invalid_argument("inverted physical epsilon demand");
  return {epsilon_gauge_detail::subtract(physical.low,gauge.shifts[row]),
    epsilon_gauge_detail::subtract(physical.high,gauge.shifts[row])};
}
inline InclusiveEpsilonDemand physical_epsilon_demand(const EpsilonGaugeResult& gauge,
    std::size_t row,const InclusiveEpsilonDemand& gauged) {
  if(row>=gauge.shifts.size())throw std::out_of_range("epsilon demand row index");
  if(gauged.low>gauged.high)throw std::invalid_argument("inverted gauged epsilon demand");
  return {epsilon_gauge_detail::add(gauged.low,gauge.shifts[row]),
    epsilon_gauge_detail::add(gauged.high,gauge.shifts[row])};
}
// For physical observable O=P Y, return P D so O=(P D) Z. This performs
// exact column multiplication only, preserving cancellations and ordering.
inline ExactEpsilonMatrix gauge_observable_columns(const EpsilonGaugeResult& gauge,
    const ExactEpsilonMatrix& physical) {
  ExactEpsilonMatrix result=physical;
  for(auto& row:result) {
    if(row.size()!=gauge.shifts.size())throw std::invalid_argument("observable/gauge column count mismatch");
    for(std::size_t j=0;j<row.size();++j) {
      if(gauge.matrix.empty() || gauge.matrix[0].empty())throw std::invalid_argument("empty epsilon gauge result");
      (void)(row[j]+gauge.matrix[0][0].constant(0));
      row[j]=epsilon_gauge_detail::multiply_power(row[j],gauge.epsilon_variable,gauge.shifts[j]);
    }
  }
  return result;
}
} // namespace diffexp
