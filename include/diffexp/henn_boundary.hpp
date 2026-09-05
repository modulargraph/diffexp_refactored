#pragma once
#include "diffexp/henn_observables.hpp"
#include "diffexp/laurent_transport.hpp"

namespace diffexp::henn {
namespace boundary_detail {
using B=kernel::ComplexBall;
inline std::optional<int> minimum_valuation(const CanonicalBasis& basis) {
  if(basis.field.variables()!=std::vector<std::string>{"eps"} || basis.components.empty() || basis.components.size()>10000)
    throw std::invalid_argument("invalid Henn canonical coefficient field or rows");
  std::optional<int> lowest;
  for(const auto& row:basis.components)for(const auto& [index,c]:row) {
    if(index.size()!=11)throw std::invalid_argument("Henn projection requires all 11 scalar slots");
    for(const auto* part:{&c.rational,&c.parity_odd}) {
      if(part->variables()!=basis.field.variables())throw std::invalid_argument("Henn projection coefficient context mismatch");
      if(auto v=exact_epsilon_valuation(*part,0)) {
        if(*v< -100 || *v>100)throw std::invalid_argument("Henn projection epsilon valuation exceeds limit");
        lowest=lowest?std::min(*lowest,static_cast<int>(*v)):static_cast<int>(*v);
      }
    }
  }
  return lowest;
}
struct Series {int low=0;std::vector<Rational> values;};
inline Series series(const Exact& coefficient,int top) {
  auto valuation=exact_epsilon_valuation(coefficient,0);
  if(!valuation || top<*valuation)return {};
  const int depth=top-*valuation;
  if(depth>1000)throw std::invalid_argument("Henn coefficient expansion exceeds limit");
  auto regular=*valuation<0?coefficient*coefficient.variable(0).pow(-*valuation):
      coefficient/coefficient.variable(0).pow(*valuation);
  auto values=feynman::scalar_functional_detail::epsilon_series(regular,std::size_t(0),depth);
  Series result;result.low=*valuation;for(const auto& value:values)result.values.push_back(value.rational());return result;
}
inline B rational_ball(const Rational& value,slong bits) {
  fmpq_t q;fmpq_init(q);fmpq_set_str(q,value.str().c_str(),10);
  B result;arb_set_fmpq(acb_realref(result.raw()),q,bits);fmpq_clear(q);return result;
}
inline B polynomial(const std::vector<Rational>& coefficients,const B& gamma,slong bits) {
  B value;
  for(auto it=coefficients.rbegin();it!=coefficients.rend();++it) {
    acb_mul(value.raw(),value.raw(),gamma.raw(),bits);auto c=rational_ball(*it,bits);acb_add(value.raw(),value.raw(),c.raw(),bits);
  }
  return value;
}
} // namespace boundary_detail

// Global normalization shifts epsilon by +4. This is the needed upper order
// of already-reduced scalar targets, NOT of masters before epsilon-polar FIRE
// reduction coefficients. The all-zero map needs no input; its returned bound
// is merely a conventional conservative value.
inline int needed_scalar_high(const CanonicalBasis& basis,int output_high=4) {
  if(output_high<0 || output_high>100)throw std::invalid_argument("Henn output epsilon order must lie in 0..100");
  return output_high-CanonicalBasis::normalization_epsilon_power-boundary_detail::minimum_valuation(basis).value_or(0);
}

inline LaurentBoundary project_boundary(const CanonicalBasis& basis,
    const std::vector<Integral>& ordered_scalar_targets,const LaurentBoundary& source,
    int output_high=4,slong bits=kernel::ComplexBall::precision()) {
  using B=kernel::ComplexBall;
  const int needed=needed_scalar_high(basis,output_high),source_high=source.high();
  const auto minimum=boundary_detail::minimum_valuation(basis);
  if(bits<64 || bits>4096 || source.values.size()!=ordered_scalar_targets.size())
    throw std::invalid_argument("Henn source shape or precision contract");
  std::map<Integral,std::size_t> positions;
  for(std::size_t i=0;i<ordered_scalar_targets.size();++i) {
    if(ordered_scalar_targets[i].size()!=11 || !positions.emplace(ordered_scalar_targets[i],i).second)
      throw std::invalid_argument("Henn source targets must be unique complete index vectors");
    for(const auto& value:source.values[i])if(!value.is_finite())throw std::invalid_argument("Henn source has a nonfinite coefficient");
  }
  if(minimum && needed>source_high)
    throw BoundaryDemand(needed,"Henn canonical projection needs additional raw scalar epsilon coefficients");
  // Structural zeros below source.low are part of LaurentBoundary's contract.
  // Every possible negative output order is retained for an explicit audit.
  const int low=minimum?std::min(0,source.low+4+*minimum):0;
  if(low< -1000 || output_high-low>1000)throw std::invalid_argument("Henn projected Laurent window exceeds limit");
  Boundary result(basis.components.size(),std::vector<B>(output_high-low+1));
  if(!minimum)return {low,std::move(result),source.taylor_tail_certified};
  B gamma,odd;arb_const_euler(acb_realref(gamma.raw()),bits);
  arb_set_ui(acb_imagref(odd.raw()),3);arb_sqrt(acb_imagref(odd.raw()),acb_imagref(odd.raw()),bits);
  std::vector<Rational> exponential(output_high-low+1,Rational(1));
  for(unsigned n=1;n<exponential.size();++n)exponential[n]=exponential[n-1]*Rational(2)/Rational(n);
  for(std::size_t row=0;row<basis.components.size();++row)for(const auto& [integral,c]:basis.components[row]) {
    if(c.is_zero())continue;
    auto found=positions.find(integral);if(found==positions.end())throw std::invalid_argument("Henn projection is missing a required scalar target");
    const auto real=boundary_detail::series(c.rational,output_high-4-source.low);
    const auto imaginary=boundary_detail::series(c.parity_odd,output_high-4-source.low);
    for(int target=low;target<=output_high;++target)for(int raw=source.low;raw<=source_high;++raw) {
      const int power=target-4-raw;
      std::vector<Rational> a(exponential.size()),b(exponential.size());bool nonzero=false;
      auto accumulate=[&](const boundary_detail::Series& series,std::vector<Rational>& polynomial) {
        for(unsigned j=0;j<series.values.size();++j) {
          const int gamma_power=power-series.low-static_cast<int>(j);
          if(gamma_power<0 || series.values[j].is_zero())continue;
          if(static_cast<std::size_t>(gamma_power)>=exponential.size())throw std::logic_error("Henn normalization window invariant");
          polynomial[gamma_power]+=series.values[j]*exponential[gamma_power];nonzero=true;
        }
      };
      accumulate(real,a);accumulate(imaginary,b);if(!nonzero)continue;
      // Combine exact rational coefficients of gamma^p and the algebraic odd
      // part BEFORE one multiplication by this source coefficient. The input
      // type carries independent balls; correlations between distinct source
      // coordinates cannot be reconstructed here and are never assumed.
      auto weight=boundary_detail::polynomial(a,gamma,bits),odd_weight=boundary_detail::polynomial(b,gamma,bits);
      acb_mul(odd_weight.raw(),odd_weight.raw(),odd.raw(),bits);acb_add(weight.raw(),weight.raw(),odd_weight.raw(),bits);
      if(weight.is_zero())continue;
      const auto& input=source.values[found->second][raw-source.low];B term;
      if(acb_is_one(weight.raw()))acb_set(term.raw(),input.raw());
      else acb_mul(term.raw(),weight.raw(),input.raw(),bits);
      if(result[row][target-low].is_zero())acb_set(result[row][target-low].raw(),term.raw());
      else acb_add(result[row][target-low].raw(),result[row][target-low].raw(),term.raw(),bits);
    }
  }
  return {low,std::move(result),source.taylor_tail_certified};
}

struct ForbiddenPole {
  std::size_t component; // 1-based canonical row
  int epsilon_order;
  kernel::ComplexBall coefficient;
};
struct PoleAudit {
  bool pass=true;
  std::size_t checked_coefficients=0;
  std::vector<ForbiddenPole> failures;
};
// Check the full absolute magnitude, including ball uncertainty. Containing
// zero alone is insufficient. This operation never removes negative orders;
// a successful tolerance audit is not a proof of exact analytic pole absence.
inline PoleAudit audit_negative_poles(const LaurentBoundary& canonical,
    const Rational& absolute_tolerance,slong bits=kernel::ComplexBall::precision()) {
  if(absolute_tolerance.sign()<=0 || bits<64 || bits>4096)throw std::invalid_argument("invalid Henn pole audit tolerance/precision");
  const auto high=canonical.high();auto tolerance=boundary_detail::rational_ball(absolute_tolerance,bits);
  mag_t threshold,bound;mag_init(threshold);mag_init(bound);
  arb_get_mag_lower(threshold,acb_realref(tolerance.raw()));PoleAudit result;
  for(std::size_t row=0;row<canonical.values.size();++row)for(int order=canonical.low;order<=std::min(-1,high);++order) {
    const auto& coefficient=canonical.values[row][order-canonical.low];++result.checked_coefficients;
    acb_get_mag(bound,coefficient.raw());
    if(!coefficient.is_finite() || mag_cmp(bound,threshold)>0) {
      result.pass=false;result.failures.push_back({row+1,order,coefficient});
    }
  }
  mag_clear(threshold);mag_clear(bound);return result;
}
} // namespace diffexp::henn
