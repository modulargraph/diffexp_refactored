#pragma once
// Cache-only bounded discovery of exact internal physical coordinates.
// scan_candidates exposes a demand-first frontier to a caller-supplied cost.
#include "diffexp/level_preparation.hpp"
#include "diffexp/recursion_physical_basis.hpp"
#include <chrono>
#include <functional>
#include <cmath>
#include <map>
#include <optional>
#include <set>

namespace diffexp::physical_basis_scan {
using Matrix = ExactEpsilonMatrix;
struct Options {
  double seconds = 1;
  std::size_t max_trials = 256, max_source_equations = 4096;
  std::size_t epsilon_variable = 1;
  bool use_source_identities = true;
  // Discovery bounds, separate from the caller's complete demand acceptance.
  std::size_t max_candidates = 8, max_candidate_rows = 256, max_depth = 2;
  std::vector<int> epsilon_powers{-1, 0, 1};
  bool sparse_combinations = true, allow_lower_sector_mixes = true;
  // Additional exact coordinate rows (in the ORIGINAL ordered basis).
  // These are algebraic coordinates, not claims about named IBP integrals.
  Matrix coordinate_rows;
};
struct Stats {
  std::size_t candidates = 0, numerator_candidates = 0, trials = 0;
  std::size_t accepted = 0, singular = 0, new_poles = 0, failures = 0;
  std::size_t original_characters = 0, selected_characters = 0;
  std::size_t recovered_coordinates = 0;
  double seconds = 0, source_recovery_seconds = 0;
  // budget_reached also covers the deterministic trial bound; time_expired is
  // the wall-clock part alone. Only the latter makes a search incomplete in a
  // machine-dependent way, so only it forbids persisting a negative decision.
  bool budget_reached = false, time_expired = false, exact_verified = false;
  std::string reason;
};
inline std::size_t expression_characters(const Matrix &a) {
  std::size_t count = 0;
  for (const auto &row : a)
    for (const auto &v : row)
      if (!v.is_zero())
        count += v.str().size();
  return count;
}
// Strip every original x-dependent denominator divisor, including repeated
// powers. Residual epsilon-only factors do not introduce a new x singularity.
inline std::set<std::string>
additional_parameter_poles(const Matrix &original,
                           const std::vector<const Matrix *> &matrices,
                           std::size_t xi) {
  std::vector<Exact> denominators;
  std::set<std::string> old_seen, checked, additional;
  for (const auto &row : original)
    for (const auto &v : row) {
      auto den = v.denominator();
      if (!den.derivative(xi).is_zero() && old_seen.insert(den.str()).second)
        denominators.push_back(std::move(den));
    }
  for (const auto *matrix : matrices)
    for (const auto &row : *matrix)
      for (const auto &v : row) {
        auto den = v.denominator();
        if (!checked.insert(den.str()).second)
          continue;
        for (const auto &old : denominators) {
          while (!den.derivative(xi).is_zero()) {
            auto gcd = den * old / den.polynomial_lcm(old);
            if (gcd.is_rational())
              break;
            auto smaller = den / gcd;
            if (smaller == den)
              break;
            den = std::move(smaller);
          }
        }
        if (!den.derivative(xi).is_zero())
          additional.insert(den.str());
      }
  return additional;
}
inline bool no_new_parameter_poles(const Matrix &original,
                                   const std::vector<const Matrix *> &matrices,
                                   std::size_t xi) {
  return additional_parameter_poles(original, matrices, xi).empty();
}
namespace detail {
using Sector = std::vector<bool>;
inline Sector sector(const ibp::Integral &a, std::size_t physical) {
  Sector s;
  for (std::size_t k = 0; k < physical; ++k)
    s.push_back(a.at(k) > 0);
  return s;
}
} // namespace detail

// Local metadata is diagnostic only: no epsilon profile here certifies the
// whole recursion's source demand, singular endpoint prescription or contour.
struct Metadata {
  std::string provenance;
  std::size_t depth = 0, connection_characters = 0, transform_characters = 0;
  std::optional<std::int64_t> connection_minimum, to_basis_minimum,
      from_basis_minimum, requested_minimum;
  bool exact_verified = false;
};
using DemandCost = std::vector<std::uint64_t>;
struct Candidate {
  recursion::PhysicalBasisTransform transform;
  Metadata metadata;
  std::optional<DemandCost> demand_cost;
};
using Evaluator = std::function<std::optional<DemandCost>(const Candidate &)>;
struct CandidateResult {
  // Pareto frontier if an evaluator is supplied; otherwise a bounded discovery
  // beam. The original basis is evaluated by the same callback as every change.
  std::vector<Candidate> candidates;
  // Callback nullopt means unresolved/not admissible, never an acceptance proof.
  std::vector<Candidate> unresolved;
  Stats stats;
  std::size_t evaluated = 0, dominated = 0, duplicates = 0, truncated = 0;
};
inline bool dominates(const DemandCost &a, const DemandCost &b) {
  if (a.size() != b.size() || a.empty())
    throw std::invalid_argument("basis demand cost axes must be fixed and nonempty");
  bool strict = false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (a[i] > b[i]) return false;
    strict |= a[i] < b[i];
  }
  return strict;
}
inline Metadata metadata(const recursion::PhysicalBasisTransform &basis,
                         const Matrix &requests, std::size_t ei,
                         std::string provenance, std::size_t depth) {
  Metadata m;
  m.provenance = std::move(provenance); m.depth = depth;
  m.connection_characters = expression_characters(basis.connection);
  m.transform_characters = expression_characters(basis.to_basis) +
                           expression_characters(basis.from_basis);
  auto minimum = [ei](const Matrix &matrix) {
    std::optional<std::int64_t> result;
    for (const auto &row : matrix) for (const auto &value : row)
      if (auto v = exact_epsilon_valuation(value, ei); v && (!result || *v < *result))
        result = v;
    return result;
  };
  m.connection_minimum = minimum(basis.connection);
  m.to_basis_minimum = minimum(basis.to_basis);
  m.from_basis_minimum = minimum(basis.from_basis);
  if (!requests.empty())
    m.requested_minimum = minimum(fuchsify::detail::multiply(requests, basis.from_basis));
  m.exact_verified = true;
  return m;
}

// Demand-first bounded discovery. Every distinct, pole-safe, exact candidate is
// offered to evaluate BEFORE cost pruning, including larger expressions and
// candidates with local epsilon poles. Only caller-supplied full-demand costs
// define Pareto dominance. The beam cap is an explicit incomplete search bound.
// No new IBP generation, numerical transport or epsilon specialization occurs.
inline CandidateResult scan_candidates(const level::Result &level,
    std::size_t physical_count, std::size_t parameter,
    const Options &options = {}, const Evaluator &evaluate = {}) {
  using namespace fuchsify::detail;
  CandidateResult result;
  const auto start = std::chrono::steady_clock::now();
  auto elapsed = [&] { return std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count(); };
  auto expired = [&] { return elapsed() >= options.seconds; };
  if (!std::isfinite(options.seconds) || options.seconds <= 0 || options.seconds > 60 ||
      !options.max_trials || !options.max_candidates || !options.max_candidate_rows ||
      !options.max_depth)
    throw std::invalid_argument("basis discovery requires finite positive bounds");
  const auto d = level.ordered_basis.size();
  if (!d || level.matrix.size() != d || level.matrix[0].size() != d ||
      level.requested_integrals.size() != level.target_rows.size())
    throw std::invalid_argument("basis discovery dimensions");
  const auto &z = level.matrix[0][0];
  if (parameter >= z.variable_count() || options.epsilon_variable >= z.variable_count())
    throw std::invalid_argument("basis discovery variables");
  for (const auto &a : level.ordered_basis)
    if (a.size() < physical_count) throw std::invalid_argument("basis discovery sectors");
  for (int p : options.epsilon_powers)
    if (p < -32 || p > 32) throw std::invalid_argument("basis discovery epsilon power bound");
  const auto unit = identity(d,z);
  recursion::PhysicalBasisTransform original{level.matrix,unit,unit,level.matrix};
  recursion::verify_physical_basis(original,level.matrix,parameter);
  result.stats.original_characters = expression_characters(level.matrix);
  std::set<std::string> seen;
  std::optional<std::size_t> cost_axes;
  auto submit = [&](recursion::PhysicalBasisTransform basis, std::string why,
                    std::size_t depth) {
    std::string key;
    for (const auto &row : basis.to_basis) for (const auto &v : row) {
      auto s = v.str(); key += std::to_string(s.size()) + ":" + s;
    }
    if (!seen.insert(key).second) { ++result.duplicates; return; }
    if (!no_new_parameter_poles(level.matrix,
        {&basis.to_basis,&basis.from_basis,&basis.connection},parameter)) {
      ++result.stats.new_poles; return;
    }
    recursion::verify_physical_basis(basis,level.matrix,parameter);
    if (!level.target_rows.empty() &&
        multiply(multiply(level.target_rows,basis.from_basis),basis.to_basis) != level.target_rows)
      throw std::logic_error("basis candidate requested-coordinate identity");
    Candidate candidate{std::move(basis),{},std::nullopt};
    candidate.metadata = metadata(candidate.transform,level.target_rows,
        options.epsilon_variable,std::move(why),depth);
    ++result.stats.accepted; result.stats.exact_verified = true;
    if (evaluate) {
      ++result.evaluated;
      candidate.demand_cost = evaluate(candidate);
      if (!candidate.demand_cost) {
        if (result.unresolved.size() < options.max_candidates)
          result.unresolved.push_back(std::move(candidate));
        else ++result.truncated;
        return;
      }
      if (candidate.demand_cost->empty() ||
          (cost_axes && *cost_axes != candidate.demand_cost->size()))
        throw std::invalid_argument("basis evaluator changed demand cost axes");
      cost_axes = candidate.demand_cost->size();
      for (const auto &old : result.candidates)
        if (dominates(*old.demand_cost,*candidate.demand_cost)) {
          ++result.dominated; return;
        }
      std::erase_if(result.candidates,[&](const Candidate &old) {
        return dominates(*candidate.demand_cost,*old.demand_cost);
      });
    }
    result.candidates.push_back(std::move(candidate));
    // Equal or incomparable costs retain deterministic discovery diversity.
    // Never use expression length to reject an unevaluated candidate.
    if (result.candidates.size() > options.max_candidates) {
      result.candidates.erase(result.candidates.begin()+1);
      ++result.truncated;
    }
  };
  submit(original,"original",0);
  // A gauge is offered separately: never silently normalize away a candidate's
  // epsilon powers before the consumer can assess its interface products.
  try {
    auto gauge = epsilon_diagonal_gauge(level.matrix,options.epsilon_variable);
    auto basis = original;
    for (std::size_t i=0;i<d;++i) {
      basis.to_basis[i][i] = epsilon_gauge_detail::multiply_power(z.constant(1),options.epsilon_variable,-gauge.shifts[i]);
      basis.from_basis[i][i] = epsilon_gauge_detail::multiply_power(z.constant(1),options.epsilon_variable,gauge.shifts[i]);
    }
    basis.connection = std::move(gauge.matrix);
    submit(std::move(basis),"epsilon diagonal gauge",0);
  } catch (const UnsupportedEpsilonGauge &) { ++result.stats.failures; }
  struct Row { std::vector<Exact> values; std::string why; };
  std::vector<Row> rows;
  auto add_row = [&](const std::vector<Exact> &row,std::string why) {
    if(row.size()!=d) throw std::invalid_argument("basis coordinate row width");
    for(const auto &v:row) if(v.variables()!=z.variables())
      throw std::invalid_argument("basis coordinate row field");
    if(rows.size()>=options.max_candidate_rows) return;
    if(std::none_of(rows.begin(),rows.end(),[&](const Row &r){return r.values==row;}))
      rows.push_back({row,std::move(why)});
  };
  for(const auto &row:options.coordinate_rows) add_row(row,"supplied exact coordinate");
  for(std::size_t i=0;i<level.target_rows.size();++i)
    add_row(level.target_rows[i],"cached requested integral");
  if(options.use_source_identities && !level.source_identities.empty() &&
      level.source_identities.size()<=options.max_source_equations && !expired()) {
    auto began=elapsed();
    ibp::ExactReducer reducer(z,level.source_identities.size());
    std::set<ibp::Integral> possible;
    for(const auto &relation:level.source_identities) {
      if(elapsed()>=options.seconds*.25) break;
      reducer.insert(relation);
      for(const auto &[a,c]:relation) possible.insert(a);
    }
    ibp::BasisReduction coordinates(reducer,level.ordered_basis,z);
    for(const auto &a:possible) {
      if(elapsed()>=options.seconds*.35 || rows.size()>=options.max_candidate_rows) break;
      try {
        auto row=coordinates.resolve({{a,z.constant(1)}});
        ibp::Relation residual{{a,z.constant(1)}};
        for(std::size_t k=0;k<d;++k) ibp::add(residual,level.ordered_basis[k],-row[k]);
        auto proof=reducer.reduce(residual);
        if(!proof.remainder.empty() || !reducer.verify(residual,proof)) continue;
        add_row(row,"cached IBP integral with exact witness");
        ++result.stats.recovered_coordinates;
        result.stats.numerator_candidates += std::any_of(a.begin(),a.end(),[](int p){return p<0;});
      } catch(const std::runtime_error &) { ++result.stats.failures; }
    }
    result.stats.source_recovery_seconds=elapsed()-began;
  }
  result.stats.candidates=rows.size();
  auto legitimate_mix = [&](std::size_t pos,std::size_t other) {
    auto a=detail::sector(level.ordered_basis[pos],physical_count);
    auto b=detail::sector(level.ordered_basis[other],physical_count);
    if(a==b) return true;
    if(!options.allow_lower_sector_mixes) return false;
    for(std::size_t k=0;k<a.size();++k) if(b[k]&&!a[k]) return false;
    return true;
  };
  auto trial = [&](const Candidate &base,std::size_t pos,
                   const std::vector<Exact> &q,const std::string &why) {
    if(expired() || result.stats.trials>=options.max_trials) return;
    ++result.stats.trials;
    if(q[pos].is_zero()) { ++result.stats.singular; return; }
    auto h=unit,hi=unit; h[pos]=q;
    for(std::size_t k=0;k<d;++k)
      hi[pos][k]=k==pos ? z.constant(1)/q[pos] : -q[k]/q[pos];
    auto rhs=multiply(h,base.transform.connection);
    for(std::size_t k=0;k<d;++k) rhs[pos][k]=rhs[pos][k]+q[k].derivative(parameter);
    recursion::PhysicalBasisTransform next{level.matrix,
      multiply(h,base.transform.to_basis),multiply(base.transform.from_basis,hi),multiply(rhs,hi)};
    submit(std::move(next),base.metadata.provenance+"; "+why,base.metadata.depth+1);
  };
  // Each depth is expanded from a fixed snapshot, so accepting a candidate
  // cannot invalidate an in-flight base. Later depths permit paired changes.
  std::vector<Candidate> bases{{original,metadata(original,level.target_rows,options.epsilon_variable,"original",0),{}}};
  for(std::size_t depth=0;depth<options.max_depth && !expired();++depth) {
    const auto depth_limit = std::min(options.max_trials,
        (depth+1)*(options.max_trials/options.max_depth) + options.max_trials%options.max_depth);
    for(const auto &base:bases) {
      auto available = [&] { return !expired() && result.stats.trials<depth_limit; };
      // Cover all component rescalings before allowing a large cached-row or
      // shear family at one position to exhaust the trial budget.
      for(int p:options.epsilon_powers) if(p)
        for(std::size_t pos=0;pos<d && available();++pos) {
          auto q=unit[pos];q[pos]=epsilon_gauge_detail::multiply_power(q[pos],options.epsilon_variable,p);
          trial(base,pos,q,"epsilon scale row "+std::to_string(pos)+" power "+std::to_string(p));
        }
      const auto rows_limit = result.stats.trials + (depth_limit-result.stats.trials+1)/2;
      for(const auto &row:rows) {
        if(!available() || result.stats.trials>=rows_limit) break;
        // The row in the base coordinates does not depend on the replaced position.
        const auto q=multiply(Matrix{row.values},base.transform.from_basis)[0];
        for(std::size_t pos=0;pos<d && available() && result.stats.trials<rows_limit;++pos)
          trial(base,pos,q,row.why);
        if(!available() || result.stats.trials>=rows_limit) break;
      }
      // Interleave destination rows at each distance to avoid first-row bias.
      if(options.sparse_combinations)
        for(std::size_t distance=1;distance<d && available();++distance)
          for(int p:options.epsilon_powers)
            for(int sign:{-1,1})
              for(std::size_t pos=0;pos<d && available();++pos) {
                const auto j=(pos+distance)%d;
                if(!legitimate_mix(pos,j)) continue;
                auto q=unit[pos];q[j]=epsilon_gauge_detail::multiply_power(z.constant(sign),options.epsilon_variable,p);
                trial(base,pos,q,"sparse row "+std::to_string(pos)+" plus "+std::to_string(sign)+" eps^"+std::to_string(p)+" row "+std::to_string(j));
              }
      if(!available()) break;
    }
    bases.clear();
    for(const auto &candidate:result.candidates)
      if(candidate.metadata.depth==depth+1) bases.push_back(candidate);
    if(bases.empty() || result.stats.trials>=options.max_trials) break;
  }
  result.stats.seconds=elapsed();
  result.stats.time_expired=expired();
  result.stats.budget_reached=result.stats.time_expired || result.stats.trials>=options.max_trials;
  result.stats.selected_characters=result.candidates.empty()?result.stats.original_characters:
      result.candidates.front().metadata.connection_characters;
  result.stats.reason="bounded exact candidate frontier; full demand and endpoint acceptance belongs to caller";
  return result;
}
} // namespace diffexp::physical_basis_scan
