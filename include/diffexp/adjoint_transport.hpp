#pragma once
#include "diffexp/laurent_transport.hpp"
#include "diffexp/rational_circuit.hpp"
#include <chrono>
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
struct AdjointPhaseTimer {
  using Clock=std::chrono::steady_clock;
  double* destination=nullptr;Clock::time_point start=Clock::now();
  ~AdjointPhaseTimer(){if(destination)*destination+=std::chrono::duration<double>(Clock::now()-start).count();}
};
// Bounded counters; no output or polling inside transport.
struct AdjointConditioningStats {
  // Accepted charts, and those whose endpoint the retained Taylor tail (not
  // the pole geometry) decided.
  std::size_t charts=0,tail_limited_charts=0;
  // Charts evaluated as midpoint + homogeneous map * uncertainty, and steps
  // halved because an enclosure still lost its arithmetic reserve.
  std::size_t centered_charts=0,conditioning_subdivisions=0;
  // Row batches split to satisfy the term/operation/cell budgets.
  std::size_t circuit_batch_reductions=0;
  std::size_t circuit_preparations=0,recurrences=0,circuit_terms=0;
  // Products spent on the probe rows, at their lower order and precision.
  std::size_t probe_terms=0;
  double compilation_seconds=0,preparation_seconds=0,recurrence_seconds=0;
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
  // Relative size, as a power of two, that the retained Taylor tail of every
  // chart must stay below; each chart is as long as this and the pole
  // geometry allow. Zero selects min(2*taylor_order, working bits - 32),
  // the accuracy of a quarter-radius step at this order.
  unsigned truncation_bits=0;
  // Live Taylor-array cells per row batch; persistent coefficients retain
  // their own cap. Never lower order or precision to fit: only independent
  // observable rows may be separated.
  std::size_t max_taylor_cells=20000000;
  // Exact-term and predicted-operation budgets of one row batch. They are
  // charged per replicated observable row, so an overrun splits the batch
  // exactly as a live-cell overrun does.
  std::size_t max_circuit_terms=2000000,max_circuit_operations=200000000;
  // Zero selects the largest independent-row batch fitting the budgets.
  std::size_t max_rows_per_batch=0;
  // Evaluate a chart whose input uncertainty would be wrapped by the
  // recurrence as midpoint chart + homogeneous chart map * uncertainty: the
  // map is computed from exact unit inputs, so cancellations such as A^2=0
  // are resolved before the input intervals enter, once.
  bool centered_input=true;
  // How the uncertainty a chart inherits is carried on. `map` is the centred
  // evaluation above: a rigorous enclosure that costs one expansion per
  // master whenever a chart would wrap its input. `probe` transports
  // midpoints only and follows the size of the discarded radii along two
  // extra homogeneous rows, into which every chart's radii are added without
  // cancellation; the final radii are eight times the largest probe
  // coefficient of each epsilon order. That is an estimate along two
  // directions, not an enclosure, at the cost of two rows. A transport of a
  // single row carries both directions in one probe row, the second one in
  // its imaginary part.
  // Once a leg transports midpoints, the radii of a chart are the defect
  // radii of its recurrence (rational_circuit::expand): the rounding the
  // chart itself commits, not a bound that the absolute values of the
  // connection inflate. They are what the probes receive next.
  enum class Uncertainty { map, probe } uncertainty=Uncertainty::map;
  // The probes also carry the Taylor tail every chart admits. What a chart
  // cuts off follows from its last retained terms (rational_circuit::
  // tail_estimate), and that estimate, a vector with its signs, is added to
  // the probes after each accepted chart: summed over the rows into the
  // first probe, with alternating row signs into the second. The
  // flow then carries that estimate of the truncation error to the end of the
  // path like any other row, cancellations included, where a bound on its
  // size per chart would lose them.
  bool tail_probe=true;
  // With probes, the estimated error of the initial rows and of the result.
  // On entry: at most as many rows as there are probes (one for a single
  // transported row, else two) in the window of the initial rows, which the
  // probes start from instead of zero; no rows mean zero. On return: the
  // probe rows at the end of the path, in the window of the result, whose
  // radii are then NOT widened by the probes. A caller that applies further
  // linear maps to the result applies them to these rows as well and takes
  // their size at the very end: an error estimate passed on as a radius
  // loses its direction, and the next cancellation multiplies it by the
  // condition number instead of cancelling it too.
  LaurentRows* error=nullptr;
  // A remaining large loss of arithmetic accuracy rejects the step and
  // halves it. Zero disables this last resort. With probes it is judged on
  // the defect radii of the midpoint chart: a step is not shortened because
  // an interval bound grew, only if the chart's own rounding lost the reserve.
  unsigned max_conditioning_halvings=12;
  AdjointConditioningStats* conditioning_stats=nullptr;
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
inline NativeTailMagnitude enclosure_uncertainty(const Jet::Ball& value) {
  using B=Jet::Ball;using M=NativeTailMagnitude;
  if(!value.is_finite())return M::upper_abs(value);
  B error;
  arf_set_mag(arb_midref(acb_realref(error.raw())),arb_radref(acb_realref(value.raw())));
  arf_set_mag(arb_midref(acb_imagref(error.raw())),arb_radref(acb_imagref(value.raw())));
  return M::upper_abs(error);
}
// Absolute accuracy near zero and relative accuracy for large coefficients.
// This is an arithmetic conditioning guard, not an omitted-tail estimator.
inline NativeTailMagnitude enclosure_quality(const Jet::Ball& value) {
  using B=Jet::Ball;using M=NativeTailMagnitude;
  if(!value.is_finite())return M::upper_abs(value);
  B midpoint;acb_get_mid(midpoint.raw(),value.raw());
  return enclosure_uncertainty(value)/M::maximum(M::one(),M::lower_abs(midpoint));
}
// Of the first `count` components (all by default).
inline NativeTailMagnitude enclosure_quality(const Boundary& values,std::size_t count=static_cast<std::size_t>(-1)) {
  NativeTailMagnitude result;
  for(std::size_t i=0;i<values.size() && i<count;++i)for(const auto& value:values[i])
    result=NativeTailMagnitude::maximum(result,enclosure_quality(value));
  return result;
}
// True when some coefficient has lost the fixed half-working-precision
// reserve. The reserve is independent of chart count and of step length:
// repeated smaller losses cannot silently consume all working precision.
inline bool beyond_reserve(const Boundary& values,std::size_t count=static_cast<std::size_t>(-1)) {
  using B=Jet::Ball;
  B reserve(1);acb_mul_2exp_si(reserve.raw(),reserve.raw(),-B::precision()/2);
  const auto quality=enclosure_quality(values,count);
  return !quality.is_finite() || quality>NativeTailMagnitude::upper_abs(reserve);
}
inline bool needs_conditioning_subdivision(const Boundary& input,const Boundary& output,std::size_t count=static_cast<std::size_t>(-1)) {
  using B=Jet::Ball;using M=NativeTailMagnitude;
  B reserve(1);acb_mul_2exp_si(reserve.raw(),reserve.raw(),-B::precision()/2);
  const auto reserve_bound=M::upper_abs(reserve);
  if(input.size()!=output.size())throw std::logic_error("conditioning subdivision shape mismatch");
  // Coupling transfers existing uncertainty into initially exact zero entries
  // (especially integral accumulators). That transfer is not newly lost
  // arithmetic precision and cannot be repaired by shrinking a step to zero.
  // Keep all output radii and the 256-fold growth guard, but include the
  // incoming absolute radius in the scale of the receiving component.
  M incoming;
  for(std::size_t i=0;i<input.size() && i<count;++i)for(const auto& value:input[i])
    incoming=M::maximum(incoming,enclosure_uncertainty(value));
  for(std::size_t i=0;i<input.size() && i<count;++i) {
    if(input[i].size()!=output[i].size())throw std::logic_error("conditioning subdivision window mismatch");
    for(unsigned k=0;k<input[i].size();++k) {
      if(!output[i][k].is_finite())return true;
      B midpoint;acb_get_mid(midpoint.raw(),output[i][k].raw());
      const auto inherited=incoming/M::maximum(M::one(),M::lower_abs(midpoint));
      if(enclosure_quality(output[i][k])>
          M::maximum(reserve_bound,M::from_ui(256)*M::maximum(enclosure_quality(input[i][k]),inherited)))return true;
    }
  }
  return false;
}
} // namespace adjoint_detail

// g' = forcing - g*A along a piecewise-linear path. Every observable row is
// an independent copy of the same d-dimensional system, and a common exact
// constant source carries the forcing. A must have nonnegative epsilon
// valuation. Each chart computes its Taylor coefficients once, by the
// finite-lag rational-circuit recurrence, and only then chooses its endpoint:
// as far as the retained tail stays below the target and the poles of the
// connection allow. Ball arithmetic encloses the retained polynomials; no
// omitted Taylor tail certificate is claimed.
namespace adjoint_detail {
// The last `probes` rows of `initial` are the uncertainty probes of
// AdjointOptions::Uncertainty::probe (zero forcing, zero initial value).
inline LaurentRows transport_rows(const ExactEpsilonMatrix& matrix,LaurentRows initial,
    const ExactEpsilonMatrix& forcing,const std::vector<Exact>& vertices,const AdjointOptions& options,unsigned probes) {
  using B=Jet::Ball;namespace rc=rational_circuit;
  const auto d=initial.columns(),r=initial.coefficients.size();const unsigned order=options.taylor_order;
  if(matrix.size()!=d || forcing.size()!=r || vertices.empty() || order<8 || order>1000 || !options.max_charts_per_leg ||
      options.max_conditioning_halvings>32 || options.truncation_bits>1000000 ||
      !options.max_circuit_terms || !options.max_circuit_operations)
    throw std::invalid_argument("adjoint transport dimensions/path/budget");
  for(const auto& row:matrix)if(row.size()!=d)throw std::invalid_argument("adjoint connection must be square");
  auto [xi,ei]=path_epsilon_variables(vertices[0]);int low=std::min(0,initial.low);
  for(const auto& row:forcing) {if(row.size()!=d)throw std::invalid_argument("adjoint forcing dimensions");for(const auto& c:row)if(!c.is_zero())low=std::min(low,static_cast<int>(*exact_epsilon_valuation(c,ei)));}
  if(low< -1000 || initial.high-low>1000)throw std::invalid_argument("adjoint epsilon depth budget");
  const auto epsilon_count=static_cast<std::size_t>(initial.high-low+1);
  if((r*d+1)*epsilon_count>20000000)throw std::invalid_argument("adjoint persistent coefficient storage budget");
  // max_taylor_cells bounds the Taylor arrays of one row batch, separately
  // from the persistent result. Never lower order or precision to fit: only
  // independent observable rows may be separated.
  const auto available_dimension=options.max_taylor_cells/epsilon_count/(order+1);
  if(available_dimension<=d)throw std::invalid_argument("adjoint Taylor workspace budget cannot fit one observable row");
  // The probes are a batch of their own, after the batches of the rows.
  const std::size_t transported=r-probes;
  if(!transported)throw std::invalid_argument("adjoint transport dimensions/path/budget");
  auto batch_rows=std::min(transported,(available_dimension-1)/d);
  if(options.max_rows_per_batch)batch_rows=std::min(batch_rows,options.max_rows_per_batch);
  // A probe is an estimate of an error: it is transported at half the Taylor
  // order (the same chart then leaves it a tail of half the bits) and at a
  // third of the working precision, which costs a fraction of a row.
  struct Batch {std::size_t first,rows;unsigned order;slong bits;bool probe;};
  const slong working_bits=B::precision();
  const unsigned probe_order=std::max(8u,(order+1)/2);
  const slong probe_bits=std::min<slong>(working_bits,std::max<slong>(64,(working_bits+2)/3));
  std::vector<Batch> batches;
  const auto arrange=[&] {
    batches.clear();
    for(std::size_t first=0;first<transported;first+=batch_rows)
      batches.push_back({first,std::min(batch_rows,transported-first),order,working_bits,false});
    if(probes)batches.push_back({transported,probes,probe_order,probe_bits,true});
  };
  arrange();
  // The probes are estimates in a lower precision: the accuracy of their
  // enclosures decides nothing. Without probes every component is judged.
  const std::size_t judged=probes?transported*d:static_cast<std::size_t>(-1);
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
  const slong bits=B::precision();auto* stats=options.conditioning_stats;
  // log2 of the admitted relative Taylor tail (TaylorTailEstimate) per chart.
  const double target=-(options.truncation_bits?static_cast<double>(options.truncation_bits):
      std::min(2.0*order,static_cast<double>(bits)-32));
  // The retained polynomial also acts on the input uncertainty, whose
  // directions converge only at the pole rate: this fraction of the pole
  // distance keeps their own tail below 2^-32 of the radius they carry.
  const double fraction=std::min(0.75,std::exp2(-32.0/order));
  // Chart centres are binary numbers no finer than 1/128 of the step that
  // reached them, so exact polynomial coefficients shifted there stay short.
  const auto shorten=[](double center,double endpoint) {
    if(!(endpoint<1))return 1.0;
    int exponent;std::frexp(endpoint-center,&exponent);
    const double grid=std::ldexp(1.0,exponent-7);
    return std::floor(endpoint/grid)*grid;
  };
  ExactEpsilonMatrix transpose(d,std::vector<Exact>(d,vertices[0].constant(0)));
  for(unsigned i=0;i<d;++i)for(unsigned j=0;j<d;++j)transpose[i][j]=-matrix[j][i];
  rc::PoleCache pole_cache; // the roots of the denominator factors, found once for all legs
  // The singular points that bound a chart, of a set of exact coefficients
  // in the variable xi: the roots of their distinct denominators.
  const auto singular_points=[&](const std::vector<const Exact*>& coefficients) {
    std::set<std::string> met,seen;std::vector<B> roots;
    for(const auto* coefficient:coefficients) {
      // Epsilon-zero denominator geometry includes poles that cancel only in
      // the leading numerator. Such poles still affect higher coefficients.
      auto p=coefficient->denominator();
      if(p.is_rational() || !met.insert(p.str()).second)continue;
      auto valuation=*exact_epsilon_valuation(p,ei);
      if(valuation)p=p/p.variable(ei).pow(valuation);
      p=p.substitute(exact_point(p,ei,p.constant(0)));
      const auto& names=p.variables();auto ii=std::find(names.begin(),names.end(),"I");
      if(ii!=names.end())p=polynomial_norm(p,ii-names.begin(),p.constant(-1));
      if(p.is_zero())throw std::domain_error("adjoint denominator vanishes on the selected complex sheet");
      if(p.is_rational())continue;p=p/p.constant(p.numerator_terms()[0].coefficient);
      if(seen.insert(p.str()).second){auto found=polynomial_roots(p,xi,bits);roots.insert(roots.end(),found.begin(),found.end());}
    }
    return roots;
  };
  // Every leg is an affine image of the same connection: its singular points
  // are those of the connection in x, found once, at the chart parameter
  // t=(x-from)/width of the leg. (The conjugate points, which a norm over I
  // adds for a leg off the real axis, are as far from a real t as these.)
  std::vector<B> singular_x;
  {
    std::vector<const Exact*> coefficients;
    for(const auto& row:transpose)for(const auto& entry:row)if(!entry.is_zero())coefficients.push_back(&entry);
    for(const auto& row:forcing)for(const auto& entry:row)if(!entry.is_zero())coefficients.push_back(&entry);
    singular_x=singular_points(coefficients);
  }
  for(unsigned leg=first_leg;leg+1<vertices.size();++leg) {
    const auto& from=vertices[leg];const auto& to=vertices[leg+1];if(from==to)continue;
    double center=leg==first_leg?first_parameter:0;
    unsigned charts=leg==first_leg?first_charts:0;
    if(center==1)continue;
    // The exact d x d block of this leg, shared by every observable row, and
    // the forcing entries, which address the common source of the full state.
    std::vector<RationalLineEntry> block,sources;
    const auto width=to-from;const auto point=exact_point(from,xi,from+width*from.variable(xi));
    const rc::AffinePullback pullback(pole_cache,from,width);
    for(unsigned j=0;j<d;++j)for(unsigned k=0;k<d;++k)if(!transpose[j][k].is_zero()) {
      if(*exact_epsilon_valuation(transpose[j][k],ei)<0)throw std::domain_error("ordinary connection requires an epsilon gauge before transport");
      block.push_back({j,k,0,width*transpose[j][k].substitute(point)});
    }
    for(unsigned i=0;i<r;++i)for(unsigned j=0;j<d;++j)if(!forcing[i][j].is_zero()) {
      auto normalized=low<0?forcing[i][j]*from.variable(ei).pow(-low):forcing[i][j]/from.variable(ei).pow(low);
      sources.push_back({static_cast<unsigned>(i*d+j),static_cast<unsigned>(r*d),0,width*normalized.substitute(point)});
    }
    std::vector<B> roots;
    if(const auto start=rc::gaussian(from),step=rc::gaussian(width);start && step) {
      // Chart lengths are doubles: sixty-four bits of a singular point decide them.
      const auto ball=[](const rc::Scalar& value){return B::from_strings(value.real.str(),value.imaginary.str());};
      const auto origin=ball(*start),scale=ball(*step);
      for(const auto& point:singular_x) {
        auto mapped=(point-origin)/scale;
        acb_set_round(mapped.raw(),mapped.raw(),64);
        roots.push_back(std::move(mapped));
      }
    } else {
      std::vector<const Exact*> coefficients;
      for(const auto* entries:{&block,&sources})for(const auto& e:*entries)coefficients.push_back(&e.coefficient);
      roots=singular_points(coefficients);
    }
    rc::Options limits;limits.max_cells=options.max_taylor_cells;
    limits.max_terms=options.max_circuit_terms;limits.max_operations=options.max_circuit_operations;
    limits.max_dimension=static_cast<unsigned>(r*d+1);
    limits.max_epsilon_degree=std::max(limits.max_epsilon_degree,static_cast<unsigned>(initial.high-low));
    std::vector<rc::Compiled> circuits;
    {
      AdjointPhaseTimer timer{stats?&stats->compilation_seconds:nullptr};
      // Storage, exact-term and operation budgets all grow with the rows of
      // a batch. Retry smaller independent-row batches transactionally if any
      // declared budget cannot be met.
      std::string overflow;
      for(;;) {
        circuits.clear();bool fits=true;
        for(std::size_t index=0;index<batches.size() && fits;++index) {
          const auto first=batches[index].first,rows=batches[index].rows,begin=first*d,count=rows*d;
          const auto order=batches[index].order;
          std::vector<RationalLineEntry> selected;
          for(std::size_t row=0;row<rows;++row)for(const auto& e:block)
            selected.push_back({static_cast<unsigned>(row*d+e.row),static_cast<unsigned>(row*d+e.column),0,e.coefficient});
          for(const auto& e:sources)if(e.row>=begin && e.row<begin+count)
            selected.push_back({static_cast<unsigned>(e.row-begin),static_cast<unsigned>(count),0,e.coefficient});
          try {
            if(selected.size()>limits.max_terms)throw std::length_error("rational circuit sparse entry budget exhausted");
            auto candidate=rc::compile(selected,count+1,xi,ei,initial.high-low,order,limits,&pullback);
            const auto& plan=candidate.data();
            (void)rc::operation_upper_bound(plan,order,static_cast<unsigned>(epsilon_count));
            // Upper bounds are small enough for uint64 by the admitted d,N,K and
            // compilation term caps; prepare() independently checks every sum.
            std::uint64_t rings=0,largest=0;
            for(const auto& stream:plan.streams)if(!plan.denominators[stream.denominator].one)
              rings+=std::min(order,stream.lag+1)*std::uint64_t(epsilon_count);
            for(const auto& poly:plan.polynomials)largest=std::max(largest,std::uint64_t(poly.size()));
            const auto shared=std::uint64_t(plan.coefficient_cells)+plan.denominators.size();
            const auto live=(std::uint64_t(order)+2)*plan.dimension*epsilon_count+rings+shared+3+
              2*rc::dot_capacity(plan,order,epsilon_count);
            if(std::max(live,shared+2*largest+3)>limits.max_cells)
              throw std::length_error("adjoint circuit cannot fit one observable row within the live-cell budget");
            circuits.push_back(std::move(candidate));
          } catch(const std::length_error& error) {fits=false;overflow=error.what();}
        }
        if(fits)break;
        if(batch_rows==1)throw std::length_error(overflow);
        batch_rows=(batch_rows+1)/2;arrange();
        if(stats)++stats->circuit_batch_reductions;
      }
    }
    // The homogeneous block alone, compiled when a chart first needs its map.
    std::optional<rc::Compiled> map_circuit;
    const bool map_fits=!probes && options.centered_input && d<=options.max_taylor_cells/(order+1)/epsilon_count/d;
    // Probe propagation is used from the first chart of a leg that would wrap
    // its input (the charts a chart map would otherwise serve); until then the
    // rows are plain interval evaluations and the probes stay exactly zero,
    // which costs nothing. `fold` moves every radius into the probes.
    bool folding=false;
    const auto fold=[&] {
      // Every row enters the chart as its midpoint. The largest discarded
      // radius of each (master, epsilon order) is added to both probes
      // along the sign they already have, so additions never cancel.
      const auto rows=r-probes;arf_t radius,largest;arf_init(radius);arf_init(largest);
      for(unsigned j=0;j<d;++j)for(std::size_t k=0;k<epsilon_count;++k) {
        arf_zero(largest);
        for(unsigned i=0;i<r;++i) {
          auto* value=state[i*d+j][k].raw();
          if(i<rows)for(auto* part:{acb_realref(value),acb_imagref(value)}) {
            arf_set_mag(radius,arb_radref(part));if(arf_cmp(radius,largest)>0)arf_set(largest,radius);
          }
          mag_zero(arb_radref(acb_realref(value)));mag_zero(arb_radref(acb_imagref(value)));
        }
        if(arf_is_zero(largest))continue;
        for(unsigned p=0;p<probes;++p) {
          auto* value=state[(rows+p)*d+j][k].raw();
          unsigned pattern=p;
          for(auto* part:{acb_realref(value),acb_imagref(value)}) {
            // The sign pattern a part starts with: all plus, or alternating
            // with the master. Two probes take one each, in both parts; a
            // single probe has the first in its real and the second in its
            // imaginary part, which a real flow transports separately.
            const int sign=arf_is_zero(arb_midref(part))?((pattern&1)&&(j&1)?-1:1):arf_sgn(arb_midref(part));
            if(probes==1)pattern=1;
            if(sign>0)arf_add(arb_midref(part),arb_midref(part),largest,bits,ARF_RND_UP);
            else arf_sub(arb_midref(part),arb_midref(part),largest,bits,ARF_RND_DOWN);
          }
        }
      }
      arf_clear(radius);arf_clear(largest);
    };
    while(center<1) {
      if(++charts>options.max_charts_per_leg)throw std::runtime_error("adjoint transport chart budget exhausted");
      // A probe is a midpoint: the enclosure an interval recurrence would
      // make of it says nothing about the error it estimates.
      for(std::size_t i=transported*d;i<r*d;++i)for(auto& value:state[i])acb_get_mid(value.raw(),value.raw());
      if(probes && (folding || beyond_reserve(state,judged))){fold();folding=true;}
      const double limit=shorten(center,clearance_endpoint(center,roots,fraction));
      if(!(limit>center) || limit>1)throw std::runtime_error("adjoint transport geometric chart made no progress");
      B c;acb_set_d(c.raw(),center);
      // Prepared coefficients and Taylor series depend on the centre only,
      // not on the endpoint: they survive every shortening of this chart.
      std::vector<std::optional<rc::Prepared>> prepared(circuits.size());
      // The Taylor series of a batch depends on the centre and on its input
      // rows, not on the endpoint: shortening the step reuses it.
      std::vector<std::optional<rc::Series>> expanded(circuits.size());
      std::optional<rc::Prepared> map_prepared;std::vector<std::optional<rc::Series>> map_columns(d);
      const auto prepare=[&](const rc::Compiled& circuit,unsigned taylor_order=0) {
        AdjointPhaseTimer timer{stats?&stats->preparation_seconds:nullptr};
        if(stats)++stats->circuit_preparations;
        return rc::prepare(circuit,c,taylor_order?taylor_order:order,epsilon_count);
      };
      const auto expand=[&](const rc::Prepared& plan,const Boundary& input,bool probe=false) {
        AdjointPhaseTimer timer{stats?&stats->recurrence_seconds:nullptr};
        // Midpoint transport: the chart's own rounding only (defect radii).
        // The acceptance below then concerns the transported midpoints, and
        // radii no longer grow with the absolute values of the connection.
        rc::Statistics work;auto series=rc::expand(plan,input,&work,probes && folding);
        if(stats){++stats->recurrences;(probe?stats->probe_terms:stats->circuit_terms)+=work.addmul_operations;}
        return series;
      };
      // Value of a series at the farthest endpoint not beyond `next` whose
      // retained tail stays below 2^bound. `next` only ever decreases.
      const auto fit=[&](const rc::Series& series,double& next,double bound) {
        for(unsigned attempt=0;;++attempt) {
          B end;acb_set_d(end.raw(),next);double tail=-std::numeric_limits<double>::infinity();
          auto values=rc::evaluate(series,end-c,&tail);
          if(!(tail>bound))return values;
          // The tail shrinks at least like step^(order-3).
          const double shorter=shorten(center,center+(next-center)*std::exp2(-(tail-bound+1)/(order-3)));
          if(attempt>=64 || !(shorter>center) || !(shorter<next))
            throw std::runtime_error("adjoint transport Taylor truncation step budget exhausted");
          next=shorter;
        }
      };
      // A Taylor chart is affine in each observable's initial vector and its
      // forcing source is exact. Split y=m+delta, evolve the exact midpoint m
      // including forcing, and apply the homogeneous chart map, computed from
      // exact unit vectors, to delta once. Every term remains a ball
      // enclosure of the same retained polynomial.
      const auto centered=[&](std::size_t index,const Boundary& input,double& next) {
        Boundary midpoint=input,delta=input;NativeTailMagnitude widest;
        for(std::size_t i=0;i<input.size();++i)for(unsigned k=0;k<epsilon_count;++k) {
          acb_get_mid(midpoint[i][k].raw(),input[i][k].raw());
          acb_sub(delta[i][k].raw(),input[i][k].raw(),midpoint[i][k].raw(),bits);
          widest=NativeTailMagnitude::maximum(widest,NativeTailMagnitude::upper_abs(delta[i][k]));
        }
        // The map only multiplies delta: its tail may be larger than the
        // target by 1/|delta|, but never worse than 2^-32 of what it carries.
        const double map_bound=std::min(-32.0,widest.is_zero()?-32.0:target-std::log2(widest.approximate_upper()));
        for(const auto& value:delta.back())if(!value.is_zero())throw std::logic_error("adjoint forcing source lost exactness");
        std::vector<unsigned> columns;
        for(unsigned j=0;j<d;++j) {
          bool uncertain=false;
          for(std::size_t first=0;first+1<input.size() && !uncertain;first+=d)
            for(const auto& value:delta[first+j])if(!value.is_zero()){uncertain=true;break;}
          if(uncertain)columns.push_back(j);
        }
        if(!map_circuit) {
          AdjointPhaseTimer timer{stats?&stats->compilation_seconds:nullptr};
          auto block_limits=limits;block_limits.max_dimension=static_cast<unsigned>(d);
          map_circuit=rc::compile(block,static_cast<unsigned>(d),xi,ei,initial.high-low,order,block_limits,&pullback);
        }
        if(!map_prepared)map_prepared=prepare(*map_circuit);
        for(auto j:columns)if(!map_columns[j]) {
          Boundary unit(d,std::vector<B>(epsilon_count,B(0)));unit[j][0]=B(1);
          map_columns[j]=expand(*map_prepared,unit);
        }
        const auto series=expand(*prepared[index],midpoint);
        Boundary values;std::vector<Boundary> maps(columns.size());
        // All polynomials of this chart are evaluated at one common endpoint.
        for(double agreed=-1;agreed!=next;) {
          agreed=next;values=fit(series,next,target);
          for(std::size_t column=0;column<columns.size();++column)maps[column]=fit(*map_columns[columns[column]],next,map_bound);
        }
        for(std::size_t first=0;first+1<input.size();first+=d)for(std::size_t column=0;column<columns.size();++column) {
          const auto& source=delta[first+columns[column]];
          for(unsigned i=0;i<d;++i)for(unsigned k=0;k<epsilon_count;++k)for(unsigned e=0;e<=k;++e)
            if(!maps[column][i][e].is_zero() && !source[k-e].is_zero())
              acb_addmul(values[first+i][k].raw(),maps[column][i][e].raw(),source[k-e].raw(),bits);
        }
        if(stats)++stats->centered_charts;
        return values;
      };
      const auto advance=[&](std::size_t index,const Boundary& input,double& next) {
        if(batches[index].probe) {
          // In their own precision and order; a probe never shortens a chart.
          kernel::PrecisionScope scope(batches[index].bits);
          if(!prepared[index])prepared[index]=prepare(circuits[index],batches[index].order);
          if(!expanded[index])expanded[index]=expand(*prepared[index],input,true);
          B end;acb_set_d(end.raw(),next);
          return rc::evaluate(*expanded[index],end-c);
        }
        if(!prepared[index])prepared[index]=prepare(circuits[index]);
        // Inherited uncertainty beyond the reserve would be wrapped by the
        // recurrence; a direct enclosure that lost the reserve was.
        if(map_fits && beyond_reserve(input))return centered(index,input,next);
        if(!expanded[index])expanded[index]=expand(*prepared[index],input);
        auto values=fit(*expanded[index],next,target);
        if(map_fits && beyond_reserve(values))return centered(index,input,next);
        return values;
      };
      double next=limit;
      for(unsigned halvings=0;;) {
        Boundary output(state.size());output.back()=state.back();
        // Row batches are independent but share the chart endpoint: a batch
        // evaluated before a later one shortened it is advanced again.
        std::vector<double> evaluated(circuits.size(),-1);
        for(bool complete=false;!complete;) {
          for(std::size_t index=0;index<batches.size();++index) {
            if(evaluated[index]==next)continue;
            const auto begin=batches[index].first*d,count=batches[index].rows*d;
            Boundary input(state.begin()+begin,state.begin()+begin+count);input.push_back(state.back());
            auto values=advance(index,input,next);
            for(std::size_t j=0;j<count;++j)output[begin+j]=std::move(values[j]);
            evaluated[index]=next;
          }
          complete=true;for(auto endpoint:evaluated)complete=complete && endpoint==next;
        }
        if(probes && !folding && beyond_reserve(output,judged)) {
          fold();folding=true;for(auto& series:expanded)series.reset(); // the inputs changed
          continue;
        }
        if(options.max_conditioning_halvings && adjoint_detail::needs_conditioning_subdivision(state,output,judged)) {
          if(halvings>=options.max_conditioning_halvings)
            throw ArithmeticConditioningFailure("adjoint arithmetic-conditioning step budget exhausted");
          const double smaller=center+(next-center)/2;
          if(!(smaller>center) || !(smaller<next))throw ArithmeticConditioningFailure("adjoint conditioning step reached floating-point geometry floor");
          next=smaller;++halvings;
          if(stats)++stats->conditioning_subdivisions;
          continue;
        }
        for(const auto& row:output)for(const auto& value:row)if(!value.is_finite())throw std::domain_error("adjoint transport produced a nonfinite coefficient");
        if(stats){++stats->charts;if(!halvings && next<limit)++stats->tail_limited_charts;}
        state=std::move(output);
        if(probes && options.tail_probe) {
          // The series a chart was evaluated from are those of its accepted
          // input; their last terms at the accepted endpoint go to the probes.
          B step;acb_set_d(step.raw(),next);acb_sub(step.raw(),step.raw(),c.raw(),bits);
          const auto rows=r-probes;B midpoint;
          for(std::size_t index=0;index<batches.size();++index) {
            if(batches[index].probe || !expanded[index])continue;
            const auto first=batches[index].first,count=batches[index].rows*d;
            const auto terms=rc::tail_estimate(*expanded[index],step);
            for(std::size_t j=0;j<count;++j) {
              const auto row=first+j/d,column=j%d;
              for(std::size_t k=0;k<epsilon_count;++k) {
                if(terms[j][k].is_zero() || !terms[j][k].is_finite())continue;
                acb_get_mid(midpoint.raw(),terms[j][k].raw());
                for(unsigned p=0;p<probes;++p) {
                  auto* probe=state[(rows+p)*d+column][k].raw();
                  if(p==1 && (row&1))acb_sub(probe,probe,midpoint.raw(),bits);
                  else acb_add(probe,probe,midpoint.raw(),bits);
                }
              }
            }
          }
        }
        break;
      }
      if(options.chart_observer || options.continuation_observer) {
        // One copy of the accepted rows serves both observers.
        AdjointContinuation accepted{leg,charts,next,LaurentRows{low,initial.high,std::vector(r,std::vector(d,std::vector<B>()))}};
        for(unsigned i=0;i<r;++i)for(unsigned j=0;j<d;++j)accepted.rows.coefficients[i][j]=state[i*d+j];
        if(options.continuation_observer)options.continuation_observer(accepted);
        if(options.chart_observer)options.chart_observer(leg,center,next,accepted.rows);
      }
      center=next;
    }
  }
  LaurentRows out{low,initial.high,std::vector(r,std::vector(d,std::vector<B>(initial.high-low+1)))};
  for(unsigned i=0;i<r;++i)for(unsigned j=0;j<d;++j)out.coefficients[i][j]=std::move(state[i*d+j]);
  return out;
}
} // namespace adjoint_detail
inline LaurentRows transport_adjoint_rows(const ExactEpsilonMatrix& matrix,LaurentRows initial,
    const ExactEpsilonMatrix& forcing,const std::vector<Exact>& vertices,const AdjointOptions& options={}) {
  using B=Jet::Ball;
  if(options.uncertainty==AdjointOptions::Uncertainty::map)
    return adjoint_detail::transport_rows(matrix,std::move(initial),forcing,vertices,options,0);
  const auto d=initial.columns(),r=initial.coefficients.size();
  // Two probe rows beside one transported row would triple the work of every
  // chart after folding starts. One complex row holds two real directions:
  // along a real path with a real connection its real and imaginary parts
  // are transported independently, exactly as the two probes are, and the
  // size below is taken over the parts separately. With a complex flow the
  // two patterns mix into one complex direction, which is still an estimate
  // of the same kind along a generic direction.
  const unsigned probes=r==1?1:2;
  if(vertices.empty())throw std::invalid_argument("adjoint transport dimensions/path/budget");
  initial.coefficients.insert(initial.coefficients.end(),probes,
    std::vector(d,std::vector<B>(initial.high-initial.low+1,B(0))));
  if(options.error && !options.error->coefficients.empty()) {
    const auto& given=*options.error;
    if(given.low!=initial.low || given.high!=initial.high || given.coefficients.size()>probes || given.columns()!=d)
      throw std::invalid_argument("adjoint transport error rows do not match the initial rows");
    for(std::size_t p=0;p<given.coefficients.size();++p)initial.coefficients[r+p]=given.coefficients[p];
  }
  auto wide=forcing;wide.insert(wide.end(),probes,std::vector<Exact>(d,vertices[0].constant(0)));
  auto out=adjoint_detail::transport_rows(matrix,std::move(initial),wide,vertices,options,probes);
  if(options.error) {
    // The caller carries the probes on; nothing of them enters the radii.
    *options.error=LaurentRows{out.low,out.high,{out.coefficients.begin()+r,out.coefficients.end()}};
    for(auto& row:options.error->coefficients)for(auto& column:row)for(auto& value:column)acb_get_mid(value.raw(),value.raw());
    out.coefficients.resize(r);
    return out;
  }
  // The probes hold the transported size of everything that was discarded.
  const slong bits=B::precision();
  arf_t size,largest;arf_init(size);arf_init(largest);
  for(std::size_t k=0;k<out.coefficients[0][0].size();++k) {
    arf_zero(largest);
    for(unsigned p=0;p<probes;++p)for(unsigned j=0;j<d;++j)
      for(const auto* part:{acb_realref(out.coefficients[r+p][j][k].raw()),acb_imagref(out.coefficients[r+p][j][k].raw())}) {
        arf_abs(size,arb_midref(part));if(arf_cmp(size,largest)>0)arf_set(largest,size);
      }
    if(arf_is_zero(largest))continue;
    arf_mul_2exp_si(largest,largest,3);
    // An exact zero is a structural zero of the transported operator: the
    // epsilon-demand analysis relies on it and no uncertainty reached it.
    for(unsigned i=0;i<r;++i)for(unsigned j=0;j<d;++j) {
      auto* value=out.coefficients[i][j][k].raw();
      if(acb_is_zero(value))continue;
      arb_add_error_arf(acb_realref(value),largest);arb_add_error_arf(acb_imagref(value),largest);
    }
  }
  arf_clear(size);arf_clear(largest);(void)bits;
  out.coefficients.resize(r);
  return out;
}
} // namespace diffexp
