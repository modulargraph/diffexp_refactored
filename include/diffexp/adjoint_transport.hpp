#pragma once
#include "diffexp/laurent_transport.hpp"
#include "diffexp/polynomial_transport.hpp"
#include <cmath>
#include <functional>
#include <memory>

namespace diffexp {
// A bounded numerical method exhausted its conditioning control. Distinct
// from invalid input, missing coefficients and failed exact preparation.
struct ArithmeticConditioningFailure : std::domain_error {
  using std::domain_error::domain_error;
  std::optional<std::size_t> recursion_depth;
};
// Coefficients of linear observables, independent of their eventual boundary.
// Orders below low vanish structurally; orders above high remain unknown.
struct LaurentRows {
  int low=0,high=0;
  std::vector<std::vector<std::vector<Jet::Ball>>> coefficients;
  std::size_t columns()const {
    if(low< -1000 || high>1000 || low>high || coefficients.empty() || coefficients.size()>5000 || coefficients[0].empty() || coefficients[0].size()>5000)
      throw std::invalid_argument("Laurent operator shape/window budget");
    const auto d=coefficients[0].size();
    if(coefficients.size()*d*static_cast<std::size_t>(high-low+1)>20000000)throw std::invalid_argument("Laurent operator coefficient storage budget");
    for(const auto& row:coefficients) {if(row.size()!=d)throw std::invalid_argument("Laurent operator row mismatch");
      for(const auto& c:row)if(c.size()!=static_cast<std::size_t>(high-low+1))throw std::invalid_argument("Laurent operator coefficient window mismatch");}
    return d;
  }
};
struct AdjointConditioningStats {
  std::size_t polynomial_charts=0, rational_cross_checks=0, rational_compilations=0;
  std::size_t centered_charts=0, homogeneous_chart_maps=0, centered_budget_skips=0, conditioning_subdivisions=0;
  std::size_t polynomial_homogeneous_columns=0, rational_homogeneous_columns=0;
};
// State after an accepted chart, expressed on the ORIGINAL path leg. Keeping
// its parameterization preserves the subsequent finite Taylor polynomials.
// This is numerical input, not a certificate of any omitted series tail.
struct AdjointContinuation {
  unsigned leg=0,accepted_charts=0;
  double parameter=0;
  LaurentRows rows;
};
struct AdjointOptions {
  unsigned taylor_order=80,max_charts_per_leg=20000;
  // Taylor-array cells per batch; persistent coefficients retain their own cap.
  std::size_t max_taylor_cells=20000000;
  bool polynomial_recurrence=true;
  // Apply inherited interval uncertainty once through the local homogeneous
  // map when both direct recurrences suffer arithmetic wrapping. This keeps
  // the same retained Taylor polynomial, epsilon window and working precision.
  bool centered_input=true;
  std::size_t max_centered_map_cells=20000000;
  // A remaining large loss of arithmetic accuracy rejects this step and
  // reduces its length. Zero preserves the original geometric chart sequence.
  unsigned max_conditioning_halvings=12;
  // Optional bounded counters; no output or polling inside transport.
  AdjointConditioningStats* conditioning_stats=nullptr;
  // Zero selects the largest independent-row batch fitting max_taylor_cells.
  std::size_t max_rows_per_batch=0;
  // Called only after a complete accepted chart. Parameters are the leg index
  // and its normalized start/end positions; rows exclude the constant source.
  std::function<void(unsigned,double,double,const LaurentRows&)> chart_observer;
  std::shared_ptr<const AdjointContinuation> continuation;
  std::function<void(const AdjointContinuation&)> continuation_observer;
};

inline LaurentRows exact_laurent_rows(const ExactEpsilonMatrix& rows,const Exact& at,int high) {
  using B=Jet::Ball;auto [xi,ei]=path_epsilon_variables(at);int low=0;
  auto evaluated=rows;for(auto& row:evaluated)for(auto& c:row){c=c.substitute(exact_point(at,xi,at));if(!c.is_zero())low=std::min(low,static_cast<int>(*exact_epsilon_valuation(c,ei)));}
  if(low< -1000 || high>1000 || high<low)throw std::invalid_argument("exact Laurent operator window budget");
  LaurentRows out{low,high,{}};
  for(const auto& row:evaluated) {out.coefficients.emplace_back();for(const auto& c:row) {
    std::vector<B> values(high-low+1,B(0));
    if(!c.is_zero()) {
      const int valuation=*exact_epsilon_valuation(c,ei);
      if(valuation<=high) {
        auto regular=valuation<0?c*c.variable(ei).pow(-valuation):c/c.variable(ei).pow(valuation);Jet epsilon(0,high-valuation+1,B::precision());if(epsilon.length()>1)epsilon.set(1,B(1));
        auto jet=evaluate(data::Reader(regular.str()).read(),epsilon,{{c.variables()[ei],epsilon}});
        for(int k=valuation;k<=high;++k)values[k-low]=jet.at(k-valuation);
      }
    }
    out.coefficients.back().push_back(std::move(values));
  }}out.columns();return out;
}

inline LaurentRows add_laurent_rows(const LaurentRows& a,const LaurentRows& b) {
  using B=Jet::Ball;auto d=a.columns();if(b.columns()!=d || a.coefficients.size()!=b.coefficients.size())throw std::invalid_argument("Laurent operator sum dimensions");
  const int low=std::min(a.low,b.low),high=std::min(a.high,b.high);
  if(high<low)throw std::invalid_argument("Laurent operator sum has no common window");
  LaurentRows out{low,high,std::vector(a.coefficients.size(),std::vector(d,std::vector<B>(high-low+1,B(0))))};
  for(const auto* input:{&a,&b})for(unsigned i=0;i<out.coefficients.size();++i)for(unsigned j=0;j<d;++j)
    for(int k=input->low;k<=high;++k)out.coefficients[i][j][k-low]+=input->coefficients[i][j][k-input->low];
  return out;
}

// Multiply column j by epsilon^shifts[j] (or its inverse), retaining the
// common known coefficient window. Shifted-out high orders remain unknown.
inline LaurentRows shift_laurent_columns(const LaurentRows& rows,
    const std::vector<std::int64_t>& shifts,bool inverse=false) {
  using B=Jet::Ball;const auto d=rows.columns();
  if(shifts.size()!=d)throw std::invalid_argument("Laurent operator gauge dimensions");
  int minimum=0;
  for(auto shift:shifts) {
    if(shift< -1000 || shift>1000)throw std::invalid_argument("Laurent column shift budget");
    minimum=std::min(minimum,static_cast<int>(inverse?-shift:shift));
  }
  const int low=rows.low+minimum,high=rows.high+minimum;
  if(low< -1000 || high>1000)throw std::invalid_argument("shifted Laurent operator window budget");
  LaurentRows out{low,high,std::vector(rows.coefficients.size(),std::vector(d,std::vector<B>(high-low+1,B(0))))};
  for(unsigned i=0;i<rows.coefficients.size();++i)for(unsigned j=0;j<d;++j) {
    const auto shift=inverse?-shifts[j]:shifts[j];
    for(int k=low;k<=high;++k)if(k-shift>=rows.low && k-shift<=rows.high)
      out.coefficients[i][j][k-low]=rows.coefficients[i][j][k-shift-rows.low];
  }
  return out;
}

inline LaurentBoundary apply_laurent_rows(const LaurentRows& rows,const LaurentBoundary& source,int output_high) {
  using B=Jet::Ball;const auto d=rows.columns();const int source_high=source.high();
  if(d!=source.values.size() || output_high< -1000 || output_high>1000)throw std::invalid_argument("Laurent operator application shape/window");
  if(output_high-rows.low>source_high)throw BoundaryDemand(output_high-rows.low,"Laurent operator requires additional source coefficients");
  if(output_high-source.low>rows.high)throw std::invalid_argument("Laurent operator requires additional operator coefficients");
  const int low=std::min(output_high,rows.low+source.low);
  if(low< -1000)throw std::invalid_argument("Laurent operator product lower-window budget");
  LaurentBoundary out{low,Boundary(rows.coefficients.size(),std::vector<B>(output_high-low+1,B(0))),false};
  for(unsigned i=0;i<rows.coefficients.size();++i)for(unsigned j=0;j<d;++j)
    for(int k=rows.low;k<=rows.high;++k)for(int n=source.low;n<=source_high && n+k<=output_high;++n)
      if(n+k>=low)out.values[i][n+k-low]+=rows.coefficients[i][j][k-rows.low]*source.values[j][n-source.low];
  return out;
}

namespace adjoint_detail {
// Absolute accuracy near zero and relative accuracy for large coefficients.
// This is an arithmetic conditioning guard, not an omitted-tail estimator.
inline NativeTailMagnitude enclosure_quality(const Jet::Ball& value) {
  using B=Jet::Ball;using M=NativeTailMagnitude;
  if(!value.is_finite())return M::upper_abs(value);
  B midpoint,error;acb_get_mid(midpoint.raw(),value.raw());
  arf_set_mag(arb_midref(acb_realref(error.raw())),arb_radref(acb_realref(value.raw())));
  arf_set_mag(arb_midref(acb_imagref(error.raw())),arb_radref(acb_imagref(value.raw())));
  return M::upper_abs(error)/M::maximum(M::one(),M::lower_abs(midpoint));
}
inline NativeTailMagnitude enclosure_quality(const Boundary& values) {
  NativeTailMagnitude result;
  for(const auto& row:values)for(const auto& value:row)
    result=NativeTailMagnitude::maximum(result,enclosure_quality(value));
  return result;
}
inline bool needs_rational_cross_check(const Boundary& input,const Boundary& output) {
  using B=Jet::Ball;using M=NativeTailMagnitude;
  B reserve(1),floor(1);acb_mul_2exp_si(reserve.raw(),reserve.raw(),-B::precision()/2);
  acb_mul_2exp_si(floor.raw(),floor.raw(),-B::precision()+16);
  const auto reserve_bound=M::upper_abs(reserve),floor_bound=M::upper_abs(floor);
  if(input.size()!=output.size())throw std::logic_error("adjoint conditioning shape mismatch");
  for(unsigned i=0;i<input.size();++i) {
    if(input[i].size()!=output[i].size())throw std::logic_error("adjoint conditioning window mismatch");
    for(unsigned k=0;k<input[i].size();++k) {
      const auto quality=enclosure_quality(output[i][k]);
      if(!output[i][k].is_finite() || !quality.is_finite() || quality>reserve_bound ||
         quality>M::from_ui(256)*M::maximum(enclosure_quality(input[i][k]),floor_bound))return true;
    }
  }
  return false;
}
inline Boundary intersect_retained_enclosures(Boundary polynomial,const Boundary& rational) {
  if(polynomial.size()!=rational.size())throw std::logic_error("adjoint enclosure intersection shape mismatch");
  for(unsigned i=0;i<polynomial.size();++i) {
    if(polynomial[i].size()!=rational[i].size())throw std::logic_error("adjoint enclosure intersection window mismatch");
    for(unsigned k=0;k<polynomial[i].size();++k) {
      auto& a=polynomial[i][k];const auto& b=rational[i][k];Jet::Ball intersection;
      if(!a.is_finite()) {
        if(!b.is_finite())throw std::domain_error("both adjoint chart recurrences produced a nonfinite coefficient");
        a=b;continue;
      }
      if(!b.is_finite())continue;
      if(!arb_intersection(acb_realref(intersection.raw()),acb_realref(a.raw()),acb_realref(b.raw()),Jet::Ball::precision()) ||
         !arb_intersection(acb_imagref(intersection.raw()),acb_imagref(a.raw()),acb_imagref(b.raw()),Jet::Ball::precision()))
        throw std::logic_error("adjoint retained-polynomial enclosures have empty intersection");
      a=std::move(intersection);
    }
  }
  return polynomial;
}
template<class Fast,class Reference>
inline Boundary conditioned_chart(const Boundary& input,Fast fast,Reference reference,AdjointConditioningStats* stats) {
  auto output=fast(input);if(stats)++stats->polynomial_charts;
  // The fixed half-working-precision reserve is independent of chart count:
  // repeated smaller losses cannot silently consume all working precision.
  if(needs_rational_cross_check(input,output)) {
    if(stats)++stats->rational_cross_checks;
    output=intersect_retained_enclosures(std::move(output),reference(input));
  }
  return output;
}
struct Entry {unsigned epsilon;data::Expr coefficient;std::vector<std::pair<unsigned,unsigned>> positions;};
inline std::vector<Entry> compile(const std::vector<RationalLineEntry>& entries) {
  std::vector<Entry> out;std::map<std::string,unsigned> common;
  for(const auto& entry:entries) {
    const auto key=std::to_string(entry.epsilon)+"|"+entry.coefficient.str();
    auto found=common.find(key);
    if(found==common.end()) {
      auto compiled=compile_rational_entries({entry});
      found=common.emplace(key,out.size()).first;out.push_back({entry.epsilon,std::move(compiled[0].coefficient),{}});
    }
    out[found->second].positions.emplace_back(entry.row,entry.column);
  }
  return out;
}
// Repeated observable rows share each exact coefficient's Taylor expansion.
inline Boundary scalar_chart(const std::vector<Entry>& entries,const Boundary& initial,const Jet::Ball& center,
    const Jet::Ball& step,unsigned order) {
  using B=Jet::Ball;const unsigned d=initial.size(),kmax=initial[0].size()-1;
  Jet x(0,order+1,B::precision());x.set(0,center);x.set(1,B(1));
  auto imaginary=x.constant(0);imaginary.set(0,B::from_strings("0","1"));
  std::vector<std::vector<B>> coefficients;
  for(const auto& entry:entries) {
    auto jet=evaluate(entry.coefficient,x,{{"x",x},{"I",imaginary}});coefficients.emplace_back();
    for(unsigned n=0;n<order;++n)coefficients.back().push_back(jet.at(n));
  }
  std::vector<B> values(static_cast<std::size_t>(order+1)*d*(kmax+1),B(0));
  const auto at=[&](unsigned n,unsigned i,unsigned k)->B&{return values[(static_cast<std::size_t>(n)*d+i)*(kmax+1)+k];};
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=kmax;++k)at(0,i,k)=initial[i][k];
  for(unsigned n=0;n<order;++n) {
    for(unsigned e=0;e<entries.size();++e)for(const auto& [row,column]:entries[e].positions)
      for(unsigned k=entries[e].epsilon;k<=kmax;++k)for(unsigned m=0;m<=n;++m)
        if(!coefficients[e][m].is_zero() && !at(n-m,column,k-entries[e].epsilon).is_zero())
          acb_addmul(at(n+1,row,k).raw(),coefficients[e][m].raw(),at(n-m,column,k-entries[e].epsilon).raw(),B::precision());
    for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=kmax;++k)acb_div_ui(at(n+1,i,k).raw(),at(n+1,i,k).raw(),n+1,B::precision());
  }
  Boundary out(d,std::vector<B>(kmax+1,B(0)));
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<=kmax;++k)for(unsigned n=order+1;n-->0;)out[i][k]=out[i][k]*step+at(n,i,k);
  return out;
}
struct ChartAcbArray {
  slong size;acb_ptr values;
  explicit ChartAcbArray(slong n):size(n),values(_acb_vec_init(n)){}
  ~ChartAcbArray(){_acb_vec_clear(values,size);}
  ChartAcbArray(const ChartAcbArray&)=delete;
};
inline Boundary dot_chart(const std::vector<adjoint_detail::Entry>& entries,const Boundary& initial,
    const Jet::Ball& center,const Jet::Ball& step,unsigned order) {
  using B=Jet::Ball;
  const unsigned d=initial.size(),width=initial[0].size();const auto bits=B::precision();
  Jet x(0,order+1,bits);x.set(0,center);x.set(1,B(1));
  auto imaginary=x.constant(0);imaginary.set(0,B::from_strings("0","1"));
  ChartAcbArray coefficients(static_cast<slong>(entries.size())*order);
  for(unsigned e=0;e<entries.size();++e){auto jet=evaluate(entries[e].coefficient,x,{{"x",x},{"I",imaginary}});
    for(unsigned n=0;n<order;++n){auto value=jet.at(n);acb_swap(coefficients.values+e*order+n,value.raw());}}
  const slong stride=static_cast<slong>(d)*width;
  ChartAcbArray values(static_cast<slong>(order+1)*stride);
  const auto at=[&](unsigned n,unsigned i,unsigned k){return values.values+static_cast<slong>(n)*stride+i*width+k;};
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<width;++k)acb_set(at(0,i,k),initial[i][k].raw());
  B sum;
  for(unsigned n=0;n<order;++n){
    for(unsigned e=0;e<entries.size();++e)for(const auto& [row,column]:entries[e].positions)
      for(unsigned k=entries[e].epsilon;k<width;++k){
        auto target=at(n+1,row,k);auto coefficient=coefficients.values+e*order;
        acb_dot(sum.raw(),target,0,coefficient,1,
            at(n,column,k-entries[e].epsilon),-stride,n+1,bits);
        acb_swap(target,sum.raw());
      }
    for(unsigned i=0;i<d;++i)for(unsigned k=0;k<width;++k)acb_div_ui(at(n+1,i,k),at(n+1,i,k),n+1,bits);
  }
  Boundary out(d,std::vector<B>(width,B(0)));
  for(unsigned i=0;i<d;++i)for(unsigned k=0;k<width;++k){
    acb_zero(out[i][k].raw());
    for(unsigned n=order+1;n-->0;){
      acb_mul(out[i][k].raw(),out[i][k].raw(),step.raw(),bits);
      acb_add(out[i][k].raw(),out[i][k].raw(),at(n,i,k),bits);
    }
  }
  return out;
}
// Apply the same existing arithmetic reserve and per-chart radius guard.
// If grouped rounding consumes the reserve, retain the scalar reference path.
inline Boundary chart(const std::vector<Entry>& entries,const Boundary& initial,
    const Jet::Ball& center,const Jet::Ball& step,unsigned order) {
  auto output=dot_chart(entries,initial,center,step,order);
  if(needs_rational_cross_check(initial,output))
    return scalar_chart(entries,initial,center,step,order);
  return output;
}

inline bool needs_centered_action(const Boundary& input,const Boundary& output) {
  using B=Jet::Ball;using M=NativeTailMagnitude;
  if(!needs_rational_cross_check(input,output))return false;
  bool uncertain=false,growing=false;
  B floor(1);acb_mul_2exp_si(floor.raw(),floor.raw(),-B::precision()+16);
  const auto floor_bound=M::upper_abs(floor);
  for(unsigned i=0;i<input.size();++i)for(unsigned k=0;k<input[i].size();++k) {
    if(!input[i][k].is_finite())throw std::domain_error("nonfinite centered chart input");
    uncertain|=!arb_is_exact(acb_realref(input[i][k].raw())) || !arb_is_exact(acb_imagref(input[i][k].raw()));
    growing|=!output[i][k].is_finite() || enclosure_quality(output[i][k])>
      M::from_ui(2)*M::maximum(enclosure_quality(input[i][k]),floor_bound);
  }
  return uncertain && growing;
}
inline bool needs_conditioning_subdivision(const Boundary& input,const Boundary& output) {
  using B=Jet::Ball;using M=NativeTailMagnitude;
  B reserve(1);acb_mul_2exp_si(reserve.raw(),reserve.raw(),-B::precision()/2);
  const auto reserve_bound=M::upper_abs(reserve);
  if(input.size()!=output.size())throw std::logic_error("conditioning subdivision shape mismatch");
  for(unsigned i=0;i<input.size();++i) {
    if(input[i].size()!=output[i].size())throw std::logic_error("conditioning subdivision window mismatch");
    for(unsigned k=0;k<input[i].size();++k)
      if(!output[i][k].is_finite() || enclosure_quality(output[i][k])>
          M::maximum(reserve_bound,M::from_ui(256)*enclosure_quality(input[i][k])))return true;
  }
  return false;
}
// A Taylor chart is affine in each observable's initial vector. Its constant
// forcing source is exact. Split y=m+delta, evolve m including forcing, and
// apply the homogeneous chart polynomial to delta only once. In particular,
// cancellations such as A^2=0 are resolved before input intervals enter.
// Every term remains a ball enclosure of the SAME retained polynomial; this
// operation neither discards radii nor supplies an omitted-series-tail bound.
template<class Reference>
inline Boundary centered_action(const Boundary& input,Reference reference,
    const LaurentRows& homogeneous,std::size_t block_dimension) {
  using B=Jet::Ball;const auto d=homogeneous.columns();
  if(!block_dimension || d!=block_dimension || homogeneous.coefficients.size()!=d ||
      homogeneous.low!=0 || input.size()<d+1 || (input.size()-1)%d ||
      input[0].size()!=static_cast<std::size_t>(homogeneous.high+1))
    throw std::invalid_argument("centered chart homogeneous map/input shape");
  const auto width=input[0].size();
  Boundary midpoint=input,delta=input;
  for(unsigned i=0;i<input.size();++i) {
    if(input[i].size()!=width)throw std::invalid_argument("centered chart epsilon window mismatch");
    for(unsigned k=0;k<width;++k) {
      if(!input[i][k].is_finite())throw std::domain_error("nonfinite centered chart input");
      acb_get_mid(midpoint[i][k].raw(),input[i][k].raw());
      acb_sub(delta[i][k].raw(),input[i][k].raw(),midpoint[i][k].raw(),B::precision());
    }
  }
  for(const auto& value:delta.back())if(!value.is_zero())
    throw std::invalid_argument("centered chart requires an exact constant forcing source");
  auto output=reference(midpoint);
  if(output.size()!=input.size())throw std::logic_error("centered chart reference shape mismatch");
  for(const auto& row:output)if(row.size()!=width)throw std::logic_error("centered chart reference window mismatch");
  for(std::size_t first=0;first+1<input.size();first+=d)
    for(unsigned i=0;i<d;++i)for(unsigned j=0;j<d;++j)
      for(unsigned k=0;k<width;++k)for(unsigned e=0;e<=k;++e)
        if(!homogeneous.coefficients[i][j][e].is_zero() && !delta[first+j][k-e].is_zero())
          acb_addmul(output[first+i][k].raw(),homogeneous.coefficients[i][j][e].raw(),
            delta[first+j][k-e].raw(),B::precision());
  return output;
}
} // namespace adjoint_detail

// g' = forcing - g*A. All rows share one compiled sparse recurrence and one
// geometric chart sequence. A must have nonnegative epsilon valuation. No
// uncertain physical boundary enters this propagation, and no omitted Taylor
// tail certificate is claimed. A common exact constant source carries forcing.
inline LaurentRows transport_adjoint_rows(const ExactEpsilonMatrix& matrix,LaurentRows initial,
    const ExactEpsilonMatrix& forcing,const std::vector<Exact>& vertices,const AdjointOptions& options={}) {
  using B=Jet::Ball;const auto d=initial.columns(),r=initial.coefficients.size();
  if(matrix.size()!=d || forcing.size()!=r || vertices.empty() || options.taylor_order<8 || options.taylor_order>1000 || !options.max_charts_per_leg ||
      !options.max_centered_map_cells || options.max_centered_map_cells>20000000 || options.max_conditioning_halvings>32)
    throw std::invalid_argument("adjoint transport dimensions/path/budget");
  for(const auto& row:matrix)if(row.size()!=d)throw std::invalid_argument("adjoint connection must be square");
  auto [xi,ei]=path_epsilon_variables(vertices[0]);int low=std::min(0,initial.low);
  for(const auto& row:forcing) {if(row.size()!=d)throw std::invalid_argument("adjoint forcing dimensions");for(const auto& c:row)if(!c.is_zero())low=std::min(low,static_cast<int>(*exact_epsilon_valuation(c,ei)));}
  if(low< -1000 || initial.high-low>1000)throw std::invalid_argument("adjoint epsilon depth budget");
  const auto epsilon_count=static_cast<std::size_t>(initial.high-low+1);
  if((r*d+1)*epsilon_count>20000000)throw std::invalid_argument("adjoint persistent coefficient storage budget");
  // max_taylor_cells bounds the large Taylor array, separately from the
  // persistent result and compiled sparse coefficients. Never lower order or
  // precision to fit: only independent observable rows may be separated.
  const auto available_dimension=options.max_taylor_cells/epsilon_count/(options.taylor_order+1);
  if(available_dimension<=d)throw std::invalid_argument("adjoint Taylor workspace budget cannot fit one observable row");
  auto batch_rows=std::min(r,(available_dimension-1)/d);
  if(options.max_rows_per_batch)batch_rows=std::min(batch_rows,options.max_rows_per_batch);
  Boundary state(r*d+1,std::vector<B>(initial.high-low+1,B(0)));state.back()[0]=B(1);
  for(unsigned i=0;i<r;++i)for(unsigned j=0;j<d;++j)for(int k=initial.low;k<=initial.high;++k)state[i*d+j][k-low]=initial.coefficients[i][j][k-initial.low];
  unsigned first_leg=0,first_charts=0;double first_parameter=0;
  if(options.continuation) {
    const auto& saved=*options.continuation;
    if(vertices.size()<2 || saved.leg>=vertices.size()-1 || !std::isfinite(saved.parameter) ||
        saved.parameter<0 || saved.parameter>1 || saved.accepted_charts>options.max_charts_per_leg ||
        (saved.parameter>0 && !saved.accepted_charts) || (saved.parameter==0 && saved.accepted_charts) ||
        vertices[saved.leg]==vertices[saved.leg+1] || saved.rows.low!=low || saved.rows.high!=initial.high ||
        saved.rows.columns()!=d || saved.rows.coefficients.size()!=r)
      throw std::invalid_argument("adjoint continuation geometry, budget or shape mismatch");
    first_leg=saved.leg;first_charts=saved.accepted_charts;first_parameter=saved.parameter;
    for(unsigned i=0;i<r;++i)for(unsigned j=0;j<d;++j) {
      for(const auto& value:saved.rows.coefficients[i][j])if(!value.is_finite())
        throw std::invalid_argument("adjoint continuation has a nonfinite coefficient");
      state[i*d+j]=saved.rows.coefficients[i][j];
    }
  }
  ExactEpsilonMatrix transpose(d,std::vector<Exact>(d,vertices[0].constant(0)));
  for(unsigned i=0;i<d;++i)for(unsigned j=0;j<d;++j)transpose[i][j]=-matrix[j][i];
  for(unsigned leg=first_leg;leg+1<vertices.size();++leg) {
    const auto& from=vertices[leg];const auto& to=vertices[leg+1];if(from==to)continue;
    double center=leg==first_leg?first_parameter:0;
    unsigned charts=leg==first_leg?first_charts:0;
    if(center==1)continue;
    // Retain rational epsilon denominators while clearing polynomial rows.
    // All replicated adjoint blocks share interned exact coefficient data.
    std::vector<RationalLineEntry> compact;
    const auto width=to-from;const auto point=exact_point(from,xi,from+width*from.variable(xi));
    for(unsigned j=0;j<d;++j)for(unsigned k=0;k<d;++k)if(!transpose[j][k].is_zero()) {
      if(*exact_epsilon_valuation(transpose[j][k],ei)<0)throw std::domain_error("ordinary connection requires an epsilon gauge before transport");
      auto coefficient=width*transpose[j][k].substitute(point);
      for(unsigned i=0;i<r;++i)compact.push_back({static_cast<unsigned>(i*d+j),static_cast<unsigned>(i*d+k),0,coefficient});
    }
    for(unsigned i=0;i<r;++i)for(unsigned j=0;j<d;++j)if(!forcing[i][j].is_zero()) {
      auto normalized=low<0?forcing[i][j]*from.variable(ei).pow(-low):forcing[i][j]/from.variable(ei).pow(low);
      compact.push_back({static_cast<unsigned>(i*d+j),static_cast<unsigned>(r*d),0,width*normalized.substitute(point)});
    }
    std::optional<polynomial_transport::Compiled> polynomial;
    std::vector<adjoint_detail::Entry> compiled;
    std::optional<std::vector<adjoint_detail::Entry>> homogeneous_compiled;
    std::optional<polynomial_transport::Compiled> homogeneous_polynomial;
    bool rational_compiled=false;
    const auto ensure_rational=[&] {
      if(rational_compiled)return;
      std::vector<RationalLineEntry> expanded;
      for(const auto& entry:compact) {
        auto coefficients=feynman::scalar_functional_detail::epsilon_series(entry.coefficient,ei,initial.high-low);
        for(unsigned k=0;k<coefficients.size();++k)if(!coefficients[k].is_zero())expanded.push_back({entry.row,entry.column,k,coefficients[k]});
      }
      compiled=adjoint_detail::compile(expanded);rational_compiled=true;
      if(options.conditioning_stats)++options.conditioning_stats->rational_compilations;
    };
    if(options.polynomial_recurrence) {
      polynomial_transport::Options limits;limits.max_cells=options.max_taylor_cells;
      limits.max_dimension=static_cast<unsigned>(r*d+1);
      limits.max_epsilon_degree=std::max(limits.max_epsilon_degree,static_cast<unsigned>(initial.high-low));
      polynomial=polynomial_transport::compile(compact,r*d+1,xi,ei,initial.high-low,options.taylor_order,limits);
    } else ensure_rational();
    std::set<std::string> seen;std::vector<B> roots;
    for(const auto& e:compact) {
      // Epsilon-zero denominator geometry includes poles that cancel only in
      // the leading numerator. Such poles still affect higher coefficients.
      auto p=e.coefficient.denominator();auto valuation=*exact_epsilon_valuation(p,ei);
      if(valuation)p=p/p.variable(ei).pow(valuation);
      p=p.substitute(exact_point(p,ei,p.constant(0)));
      const auto& names=p.variables();auto ii=std::find(names.begin(),names.end(),"I");
      if(ii!=names.end())p=polynomial_norm(p,ii-names.begin(),p.constant(-1));
      if(p.is_zero())throw std::domain_error("adjoint denominator vanishes on the selected complex sheet");
      if(p.is_rational())continue;p=p/p.constant(p.numerator_terms()[0].coefficient);
      if(seen.insert(p.str()).second){auto found=polynomial_roots(p,xi,B::precision());roots.insert(roots.end(),found.begin(),found.end());}
    }
    while(center<1) {
      if(++charts>options.max_charts_per_leg)throw std::runtime_error("adjoint transport chart budget exhausted");
      double next=clearance_endpoint(center,roots);
      if(!(next>center) || next>1)throw std::runtime_error("adjoint transport geometric chart made no progress");
      for(unsigned halvings=0;;++halvings) {
      // A rejected chart must not leak partially updated observable batches
      // into its retry. No precision, order or epsilon window is changed.
      Boundary before;
      if(options.max_conditioning_halvings)before=state;
      B c,end;acb_set_d(c.raw(),center);acb_set_d(end.raw(),next);
      // All observable blocks share this homogeneous map, including across
      // independent-row batches. Construct its columns sequentially, so the
      // Taylor workspace never exceeds the already checked single-row budget.
      std::optional<LaurentRows> homogeneous_map;
      const auto centered_guard=[&](const Boundary& input,Boundary output,const auto& reference) {
        if(!options.centered_input || !adjoint_detail::needs_centered_action(input,output))return output;
        if(d>options.max_centered_map_cells/epsilon_count/d) {
          if(options.conditioning_stats)++options.conditioning_stats->centered_budget_skips;
          return output;
        }
        if(!homogeneous_map) {
          ensure_rational();
          if(!homogeneous_compiled) {
            homogeneous_compiled.emplace();
            for(const auto& entry:compiled) {
              adjoint_detail::Entry selected{entry.epsilon,entry.coefficient,{}};
              for(const auto& [row,column]:entry.positions)if(row<d && column<d)
                selected.positions.emplace_back(row,column);
              if(!selected.positions.empty())homogeneous_compiled->push_back(std::move(selected));
            }
            if(polynomial) {
              // Compile the physical homogeneous equation independently of
              // forcing denominators. Its finite-lag recurrence describes the
              // same retained map as the rational coefficient recurrence.
              std::vector<RationalLineEntry> homogeneous_entries;
              for(const auto& entry:compact)if(entry.row<d && entry.column<d)
                homogeneous_entries.push_back(entry);
              auto limits=polynomial->options;limits.max_dimension=d;
              homogeneous_polynomial=polynomial_transport::compile(homogeneous_entries,d,xi,ei,
                  epsilon_count-1,options.taylor_order,limits);
            }
          }
          homogeneous_map=LaurentRows{0,static_cast<int>(epsilon_count)-1,
            std::vector(d,std::vector(d,std::vector<B>(epsilon_count,B(0))))};
          std::vector<bool> active_column(d,false);
          for(const auto& entry:*homogeneous_compiled)for(const auto& [row,column]:entry.positions)
            active_column[column]=true;
          for(unsigned column=0;column<d;++column) {
            // Zero connection columns have exactly constant unit solutions.
            // In particular, factored integral accumulators need no solves.
            if(!active_column[column]) {homogeneous_map->coefficients[column][column][0]=B(1);continue;}
            Boundary basis(d,std::vector<B>(epsilon_count,B(0)));basis[column][0]=B(1);
            Boundary mapped;bool rational_needed=true;
            if(homogeneous_polynomial) {
              mapped=polynomial_transport::chart(*homogeneous_polynomial,basis,c,end-c,options.taylor_order);
              if(options.conditioning_stats)++options.conditioning_stats->polynomial_homogeneous_columns;
              B reserve(1);acb_mul_2exp_si(reserve.raw(),reserve.raw(),-B::precision()/2);
              const auto quality=adjoint_detail::enclosure_quality(mapped);
              rational_needed=!quality.is_finite() || quality>NativeTailMagnitude::upper_abs(reserve);
            }
            if(rational_needed) {
              auto reference=adjoint_detail::chart(*homogeneous_compiled,basis,c,end-c,options.taylor_order);
              if(options.conditioning_stats)++options.conditioning_stats->rational_homogeneous_columns;
              mapped=homogeneous_polynomial?adjoint_detail::intersect_retained_enclosures(std::move(mapped),reference):std::move(reference);
            }
            for(unsigned row=0;row<d;++row)homogeneous_map->coefficients[row][column]=std::move(mapped[row]);
          }
          if(options.conditioning_stats)++options.conditioning_stats->homogeneous_chart_maps;
        }
        auto centered=adjoint_detail::centered_action(input,reference,*homogeneous_map,d);
        if(options.conditioning_stats)++options.conditioning_stats->centered_charts;
        return adjoint_detail::intersect_retained_enclosures(std::move(output),centered);
      };
      if(batch_rows==r) {
        const auto reference=[&](const Boundary& input){ensure_rational();return adjoint_detail::chart(compiled,input,c,end-c,options.taylor_order);};
        auto output=polynomial?adjoint_detail::conditioned_chart(state,
          [&](const Boundary& input){return polynomial_transport::chart(*polynomial,input,c,end-c,options.taylor_order);},
          reference,options.conditioning_stats):reference(state);
        state=centered_guard(state,std::move(output),reference);
      } else {
        // One exact compilation and pole geometry per leg. Each batch borrows
        // the interned polynomial pool instead of recompiling the connection.
        // Updating earlier rows is safe: blocks have no cross-observable edges.
        for(std::size_t first=0;first<r;first+=batch_rows) {
          const auto begin=first*d,count=std::min(batch_rows,r-first)*d,source=r*d;
          const auto remap=[&](unsigned column)->unsigned {
            if(column==source)return static_cast<unsigned>(count);
            if(column<begin || column>=begin+count)throw std::logic_error("adjoint batch has cross-observable coupling");
            return static_cast<unsigned>(column-begin);
          };
          Boundary batch(state.begin()+begin,state.begin()+begin+count);batch.push_back(state.back());
          const auto rational_batch=[&](const Boundary& input) {
            ensure_rational();std::vector<adjoint_detail::Entry> view;
            for(const auto& entry:compiled) {
              adjoint_detail::Entry selected{entry.epsilon,entry.coefficient,{}};
              for(const auto& [row,column]:entry.positions)if(row>=begin && row<begin+count)
                selected.positions.emplace_back(row-begin,remap(column));
              if(!selected.positions.empty())view.push_back(std::move(selected));
            }
            return adjoint_detail::chart(view,input,c,end-c,options.taylor_order);
          };
          if(polynomial) {
            polynomial_transport::Compiled view;
            view.dimension=count+1;view.epsilon_high=polynomial->epsilon_high;
            view.expected_order=polynomial->expected_order;view.options=polynomial->options;
            view.rows.assign(polynomial->rows.begin()+begin,polynomial->rows.begin()+begin+count);
            view.rows.push_back(polynomial->rows.back());
            for(auto& row:view.rows)for(auto& entry:row.entries)entry.column=remap(entry.column);
            for(const auto& entry:polynomial->fallback) {
              polynomial_transport::SharedFallback selected{entry.epsilon,entry.coefficient,{}};
              for(const auto& [row,column]:entry.positions)if(row>=begin && row<begin+count)
                selected.positions.emplace_back(row-begin,remap(column));
              if(!selected.positions.empty())view.fallback.push_back(std::move(selected));
            }
            view.polynomials=std::move(polynomial->polynomials);
            auto output=adjoint_detail::conditioned_chart(batch,
              [&](const Boundary& input){return polynomial_transport::chart(view,input,c,end-c,options.taylor_order);},
              rational_batch,options.conditioning_stats);
            polynomial->polynomials=std::move(view.polynomials);
            batch=centered_guard(batch,std::move(output),rational_batch);
          } else batch=centered_guard(batch,rational_batch(batch),rational_batch);
          for(std::size_t j=0;j<count;++j)state[begin+j]=std::move(batch[j]);
        }
      }
      if(options.max_conditioning_halvings && adjoint_detail::needs_conditioning_subdivision(before,state)) {
        if(halvings>=options.max_conditioning_halvings)
          throw ArithmeticConditioningFailure("adjoint arithmetic-conditioning step budget exhausted");
        const auto smaller=center+(next-center)/2;
        if(!(smaller>center) || !(smaller<next))throw ArithmeticConditioningFailure("adjoint conditioning step reached floating-point geometry floor");
        state=std::move(before);next=smaller;
        if(options.conditioning_stats)++options.conditioning_stats->conditioning_subdivisions;
        continue;
      }
      for(const auto& row:state)for(const auto& value:row)if(!value.is_finite())throw std::domain_error("adjoint transport produced a nonfinite coefficient");
      break;
      }
      if(options.chart_observer || options.continuation_observer) {
        LaurentRows snapshot{low,initial.high,std::vector(r,std::vector(d,std::vector<B>()))};
        for(unsigned i=0;i<r;++i)for(unsigned j=0;j<d;++j)snapshot.coefficients[i][j]=state[i*d+j];
        if(options.continuation_observer)options.continuation_observer(AdjointContinuation{leg,charts,next,snapshot});
        if(options.chart_observer)options.chart_observer(leg,center,next,snapshot);
      }
      center=next;
    }
  }
  LaurentRows out{low,initial.high,std::vector(r,std::vector(d,std::vector<B>(initial.high-low+1)))};
  for(unsigned i=0;i<r;++i)for(unsigned j=0;j<d;++j)out.coefficients[i][j]=std::move(state[i*d+j]);
  return out;
}
} // namespace diffexp
