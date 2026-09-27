#pragma once
// Fuse exact adjacent maps BEFORE assigning epsilon demands or expanding jets.
// Numerical/retained endpoint rows are deliberately not accepted as exact maps.
#include "diffexp/epsilon_gauge.hpp"
#include <map>

namespace diffexp::exact_jet_composition {
using Matrix=ExactEpsilonMatrix;
using High=std::optional<int>;
using Highs=std::vector<High>;
struct Limits {
  std::size_t max_rows=256,max_columns=256,max_products=1000000,
    max_term_work=2000000,max_coefficients=200000;
  int max_abs_order=1000;
};
namespace detail {
inline void limits(const Limits& l) {
  if(!l.max_rows || !l.max_columns || !l.max_products || !l.max_term_work ||
     !l.max_coefficients || l.max_abs_order<0 || l.max_abs_order>10000)
    throw std::invalid_argument("exact jet composition limits");
}
inline int order(std::int64_t n,const Limits& l) {
  if(n < -static_cast<std::int64_t>(l.max_abs_order) || n>l.max_abs_order)
    throw std::length_error("exact jet composition order budget");
  return static_cast<int>(n);
}
inline std::size_t shape(const Matrix& a,std::size_t epsilon,const Limits& l,const Exact* sample=nullptr) {
  limits(l);
  if(a.empty() || a.front().empty() || a.size()>l.max_rows || a.front().size()>l.max_columns)
    throw std::invalid_argument("missing exact map or matrix shape budget");
  const auto& z=sample?*sample:a.front().front();
  if(epsilon>=z.variable_count())throw std::invalid_argument("exact jet epsilon index");
  for(const auto& row:a) {
    if(row.size()!=a.front().size())throw std::invalid_argument("exact jet ragged map");
    for(const auto& v:row)z.require_same_field(v);
  }
  return a.front().size();
}
inline void raise(High& a,int b) {if(!a || *a<b)a=b;}
inline void lower(High& a,int b) {if(!a || *a>b)a=b;}
struct Budget {
  const Limits& limits;
  std::size_t products=0,term_work=0,coefficients=0;
  static void add(std::size_t& n,std::size_t delta,std::size_t limit,const char* message) {
    if(delta>limit-n)throw std::length_error(message);n+=delta;
  }
  void multiply_count(std::size_t a,std::size_t b) {
    if(a && b>limits.max_term_work/a)throw std::length_error("exact jet term-work budget");
    add(term_work,a*b,limits.max_term_work,"exact jet term-work budget");
  }
  std::pair<std::size_t,std::size_t> terms(const Exact& v) {
    auto n=v.numerator_terms().size(),d=v.denominator_terms().size();
    add(term_work,n,limits.max_term_work,"exact jet term-work budget");
    add(term_work,d,limits.max_term_work,"exact jet term-work budget");return {n,d};
  }
  void product(const Exact& a,const Exact& b) {
    add(products,1,limits.max_products,"exact jet product budget");
    auto [an,ad]=terms(a);auto [bn,bd]=terms(b);multiply_count(an,bn);multiply_count(ad,bd);
  }
  void sum(const Exact& a,const Exact& b) {
    auto [an,ad]=terms(a);auto [bn,bd]=terms(b);
    multiply_count(an,bd);multiply_count(bn,ad);multiply_count(ad,bd);
  }
  void cells(std::size_t n) {add(coefficients,n,limits.max_coefficients,"exact jet coefficient budget");}
};
inline High minimum(const Matrix& a,std::size_t epsilon,const Limits& l) {
  High result;for(const auto& row:a)for(const auto& v:row)
    if(auto value=exact_epsilon_valuation(v,epsilon))lower(result,order(*value,l));
  return result; // nullopt means identically zero, not an unknown bound.
}
} // namespace detail

struct Fusion {
  Matrix product;
  High separate_global_lower_sum,pathwise_lower,fused_lower;
  std::size_t nonzero_products=0,estimated_term_work=0;
};
inline Fusion fuse(const Matrix& outer,const Matrix& inner,std::size_t epsilon,const Limits& l={}) {
  const auto shared=detail::shape(outer,epsilon,l);
  const auto columns=detail::shape(inner,epsilon,l,&outer.front().front());
  if(shared!=inner.size())throw std::invalid_argument("exact adjacent map dimensions");
  const auto& z=outer.front().front();detail::Budget budget{l};
  // Count all nonzero scalar products before performing any matrix product.
  std::size_t preflight=0;
  for(std::size_t i=0;i<outer.size();++i)for(std::size_t k=0;k<shared;++k)if(!outer[i][k].is_zero())
    for(std::size_t j=0;j<columns;++j)if(!inner[k][j].is_zero())
      detail::Budget::add(preflight,1,l.max_products,"exact jet product budget");
  if(outer.size()>l.max_coefficients/columns)throw std::length_error("exact fused matrix storage budget");
  Fusion result;result.product.assign(outer.size(),std::vector<Exact>(columns,z.constant(0)));
  auto a=detail::minimum(outer,epsilon,l),b=detail::minimum(inner,epsilon,l);
  if(a && b)result.separate_global_lower_sum=detail::order(static_cast<std::int64_t>(*a)+*b,l);
  for(std::size_t i=0;i<outer.size();++i)for(std::size_t k=0;k<shared;++k)if(!outer[i][k].is_zero())
    for(std::size_t j=0;j<columns;++j)if(!inner[k][j].is_zero()) {
      auto valuation=*exact_epsilon_valuation(outer[i][k],epsilon)+*exact_epsilon_valuation(inner[k][j],epsilon);
      detail::lower(result.pathwise_lower,detail::order(valuation,l));
      budget.product(outer[i][k],inner[k][j]);auto term=outer[i][k]*inner[k][j];
      budget.sum(result.product[i][j],term);result.product[i][j]=result.product[i][j]+term;
      (void)budget.terms(result.product[i][j]);
    }
  result.fused_lower=detail::minimum(result.product,epsilon,l);
  result.nonzero_products=budget.products;result.estimated_term_work=budget.term_work;return result;
}
inline Fusion fuse_available(const std::optional<Matrix>& outer,const std::optional<Matrix>& inner,
                             std::size_t epsilon,const Limits& l={}) {
  if(!outer || !inner)throw std::invalid_argument("adjacent fusion requires both exact maps; numerical endpoint is not exact");
  return fuse(*outer,*inner,epsilon,l);
}

struct EntryDemand {std::size_t output,input;int low,high;};
struct Schedule {
  Highs input_highs,output_highs;
  std::vector<int> input_lows,output_lows;
  std::vector<EntryDemand> entries;
  std::size_t map_coefficients=0,output_coefficients=0;
};
// Valuations here are generic over the other rational-function variables.
// Specialize exact maps at the intended ordinary point before scheduling when
// that specialization can change epsilon valuation; numerical sampling is not valid.
// Inputs are identically zero below their independently certified lower bounds.
// nullopt output high means unrequested, not a zero output or a known prefix.
inline Schedule schedule(const Matrix& map,const std::vector<int>& input_lows,const Highs& highs,
                         std::size_t epsilon,const Limits& l={}) {
  const auto columns=detail::shape(map,epsilon,l);
  if(input_lows.size()!=columns || highs.size()!=map.size())throw std::invalid_argument("exact jet demand shape");
  for(auto low:input_lows)detail::order(low,l);
  Schedule s;s.input_highs.resize(columns);s.output_highs=highs;s.input_lows=input_lows;s.output_lows.assign(map.size(),0);
  for(std::size_t i=0;i<map.size();++i)if(highs[i]) {
    const auto high=detail::order(*highs[i],l);s.output_lows[i]=high;
    for(std::size_t j=0;j<columns;++j)if(auto valuation=exact_epsilon_valuation(map[i][j],epsilon)) {
      const auto low=detail::order(*valuation,l);
      const auto first=static_cast<std::int64_t>(low)+input_lows[j];
      if(first>high)continue; // Every requested contribution is structurally zero.
      s.output_lows[i]=std::min(s.output_lows[i],detail::order(first,l));
      const auto map_high=detail::order(static_cast<std::int64_t>(high)-input_lows[j],l);
      detail::raise(s.input_highs[j],detail::order(static_cast<std::int64_t>(high)-low,l));
      s.entries.push_back({i,j,low,map_high});
      detail::Budget::add(s.map_coefficients,static_cast<std::size_t>(map_high-low+1),l.max_coefficients,"exact jet coefficient budget");
    }
    detail::Budget::add(s.output_coefficients,static_cast<std::size_t>(high-s.output_lows[i]+1),l.max_coefficients,"exact jet coefficient budget");
  }
  if(s.output_coefficients>l.max_coefficients-s.map_coefficients)throw std::length_error("exact jet coefficient budget");
  return s;
}
inline Schedule schedule(const Matrix& map,const std::vector<int>& lows,int high,std::size_t epsilon,const Limits& l={}) {
  return schedule(map,lows,Highs(map.size(),high),epsilon,l);
}

struct Jet {
  int low,high;
  std::vector<Exact> coefficients;
  // Below low is known zero; above high is UNKNOWN, even when all stored values are zero.
};
struct ExpandedEntry {EntryDemand demand;Jet coefficients;};
namespace detail {
inline Jet expand(const Exact& f,int low,int high,std::size_t epsilon,Budget& budget) {
  const auto width=static_cast<std::size_t>(high-low+1);budget.cells(width);
  auto normalized=epsilon_gauge_detail::multiply_power(f,epsilon,-static_cast<std::int64_t>(low));
  // Check scratch before allocating it. This cap is separate from retained cells.
  if(width>budget.limits.max_coefficients/3)throw std::length_error("exact jet expansion scratch budget");
  std::vector<Exact> p(width,f.constant(0)),q(width,f.constant(0)),values(width,f.constant(0));
  auto collect=[&](const std::vector<Exact::Term>& terms,std::vector<Exact>& out) {
    Budget::add(budget.term_work,terms.size(),budget.limits.max_term_work,"exact jet term-work budget");
    for(const auto& t:terms)if(t.powers[epsilon]<width) {
      auto value=f.constant(t.coefficient);
      for(std::size_t v=0;v<t.powers.size();++v)if(v!=epsilon && t.powers[v]) {
        auto monomial=f.variable(v).pow(t.powers[v]);budget.product(value,monomial);value=value*monomial;
      }
      auto k=t.powers[epsilon];budget.sum(out[k],value);out[k]=out[k]+value;
    }
  };
  collect(normalized.numerator_terms(),p);collect(normalized.denominator_terms(),q);
  if(q[0].is_zero())throw std::logic_error("exact jet nonunit normalized denominator");
  for(std::size_t n=0;n<width;++n) {
    auto value=p[n];
    for(std::size_t j=1;j<=n;++j)if(!q[j].is_zero() && !values[n-j].is_zero()) {
      budget.product(q[j],values[n-j]);auto term=q[j]*values[n-j];budget.sum(value,term);value=value-term;
    }
    budget.product(value,q[0]); // Charge the scalar division as bounded arithmetic too.
    values[n]=value/q[0];(void)budget.terms(values[n]);
  }
  return {low,high,std::move(values)};
}
inline void validate_jet(const Jet& jet,const Exact& z,std::size_t epsilon,const Limits& l) {
  order(jet.low,l);order(jet.high,l);
  const auto width=jet.high>=jet.low?static_cast<std::size_t>(jet.high-jet.low+1):0;
  if(width!=jet.coefficients.size() || width>l.max_coefficients)throw std::invalid_argument("exact input jet shape/budget");
  for(const auto& v:jet.coefficients) {
    z.require_same_field(v);
    if(!v.derivative(epsilon).is_zero())throw std::invalid_argument("jet coefficient still depends on epsilon");
  }
}
} // namespace detail
// The public expansion recomputes the schedule, so a caller cannot forge its low/high metadata.
inline std::vector<ExpandedEntry> expand_demanded(const Matrix& map,const std::vector<int>& lows,
                                                const Highs& highs,std::size_t epsilon,const Limits& l={}) {
  auto plan=schedule(map,lows,highs,epsilon,l);detail::Budget budget{l};
  std::vector<ExpandedEntry> entries;
  for(const auto& e:plan.entries)entries.push_back({e,detail::expand(map[e.output][e.input],e.low,e.high,epsilon,budget)});
  return entries;
}
struct MissingCoefficient : std::runtime_error {
  std::size_t input;int required_high,available_high;
  MissingCoefficient(std::size_t i,int required,int available):std::runtime_error("exact composition requires unknown higher input coefficients"),
    input(i),required_high(required),available_high(available){}
};
struct Application {Schedule demand;std::vector<std::optional<Jet>> outputs;};
inline Application apply(const Matrix& map,const std::vector<Jet>& input,const Highs& highs,
                         std::size_t epsilon,const Limits& l={}) {
  const auto columns=detail::shape(map,epsilon,l);const auto& z=map.front().front();
  if(input.size()!=columns)throw std::invalid_argument("exact composition input count");
  std::vector<int> lows;for(const auto& jet:input){detail::validate_jet(jet,z,epsilon,l);lows.push_back(jet.low);}
  auto plan=schedule(map,lows,highs,epsilon,l);
  // Check all completeness before expanding any map coefficient or publishing any output.
  for(std::size_t j=0;j<columns;++j)if(plan.input_highs[j] && *plan.input_highs[j]>input[j].high)
    throw MissingCoefficient(j,*plan.input_highs[j],input[j].high);
  Application out{plan,std::vector<std::optional<Jet>>(map.size())};detail::Budget budget{l};
  for(std::size_t i=0;i<map.size();++i)if(highs[i]) {
    auto width=static_cast<std::size_t>(*highs[i]-plan.output_lows[i]+1);budget.cells(width);
    out.outputs[i]=Jet{plan.output_lows[i],*highs[i],std::vector<Exact>(width,z.constant(0))};
  }
  for(const auto& entry:plan.entries) {
    auto expanded=detail::expand(map[entry.output][entry.input],entry.low,entry.high,epsilon,budget);
    auto& dest=*out.outputs[entry.output];const auto& source=input[entry.input];
    for(int p=expanded.low;p<=expanded.high;++p)if(!expanded.coefficients[p-expanded.low].is_zero())
      for(int q=source.low;q<=dest.high-p;++q) {
        if(q>source.high)throw std::logic_error("exact composition completeness invariant");
        const auto& v=source.coefficients[q-source.low];if(v.is_zero())continue;
        const auto& a=expanded.coefficients[p-expanded.low];budget.product(a,v);auto term=a*v;
        auto& cell=dest.coefficients[p+q-dest.low];budget.sum(cell,term);cell=cell+term;
      }
  }
  return out;
}
inline Application apply(const Matrix& map,const std::vector<Jet>& input,int high,std::size_t epsilon,const Limits& l={}) {
  return apply(map,input,Highs(map.size(),high),epsilon,l);
}
} // namespace diffexp::exact_jet_composition
