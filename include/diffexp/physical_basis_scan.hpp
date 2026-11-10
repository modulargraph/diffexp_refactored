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
  bool budget_reached = false, time_expired = false;
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
// A coefficient introduces no new parameter singularity when, after stripping
// every divisor it shares with the original x-dependent denominators (repeated
// powers included), nothing x-dependent is left. Residual epsilon-only factors
// are harmless. The original denominators enter through one common multiple,
// and each distinct candidate denominator is decided once.
class PoleScreen {
  std::size_t xi_;
  std::optional<Exact> common_;
  std::map<std::string, bool> known_;

public:
  PoleScreen(const Matrix &original, std::size_t xi) : xi_(xi) {
    for (const auto &row : original)
      for (const auto &v : row) {
        auto den = v.denominator();
        if (den.derivative(xi).is_zero() || !known_.emplace(den.str(), true).second)
          continue;
        common_ = common_ ? common_->polynomial_lcm(den) : den;
      }
  }
  bool admissible(const Exact &value) {
    if (value.is_zero())
      return true;
    auto den = value.denominator();
    if (den.is_rational())
      return true;
    auto [found, inserted] = known_.try_emplace(den.str(), false);
    if (!inserted)
      return found->second;
    while (common_ && !den.derivative(xi_).is_zero()) {
      auto gcd = den * *common_ / den.polynomial_lcm(*common_);
      if (gcd.is_rational())
        break;
      auto smaller = den / gcd;
      if (smaller == den)
        break;
      den = std::move(smaller);
    }
    return found->second = den.derivative(xi_).is_zero();
  }
};
namespace detail {
using Sector = std::vector<bool>;
inline Sector sector(const ibp::Integral &a, std::size_t physical) {
  Sector s;
  for (std::size_t k = 0; k < physical; ++k)
    s.push_back(a.at(k) > 0);
  return s;
}
} // namespace detail

// A candidate is exact by construction (a verified base times one explicit
// row change); the consumer verifies the one it adopts. connection_characters
// measures the transported matrix, which is what a simpler basis saves.
struct Metadata {
  std::string provenance;
  std::size_t depth = 0, connection_characters = 0;
};
using DemandCost = std::vector<std::uint64_t>;
struct Candidate {
  recursion::PhysicalBasisTransform transform;
  Metadata metadata;
  std::optional<DemandCost> demand_cost;
  // Per-row identity of to_basis and size of the connection. A single row
  // change touches few rows, so its successor updates these, not all d^2.
  std::vector<std::string> row_keys;
  std::vector<std::size_t> row_characters;
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
// Bounded discovery of exact internal bases. Every distinct, pole-safe
// candidate is offered to evaluate BEFORE cost pruning, including larger
// expressions and candidates with local epsilon poles. Only caller-supplied
// costs define Pareto dominance. The beam cap is an explicit incomplete search
// bound. Candidates are exact by construction; they are not re-verified here.
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
  PoleScreen poles(level.matrix,parameter);
  std::set<std::string> seen;
  std::optional<std::size_t> cost_axes;
  const auto row_key = [](const std::vector<Exact> &row) {
    std::string key;
    for (const auto &v : row) { auto s = v.str(); key += std::to_string(s.size()) + ":" + s; }
    return key;
  };
  const auto row_size = [](const std::vector<Exact> &row) {
    std::size_t size = 0;
    for (const auto &v : row) if (!v.is_zero()) size += v.str().size();
    return size;
  };
  const auto row_admissible = [&](const std::vector<Exact> &row) {
    return std::all_of(row.begin(),row.end(),[&](const Exact &v){return poles.admissible(v);});
  };
  // The row masks name what differs from the candidate's base, whose keys,
  // sizes and pole screening it inherits; a null mask means every row.
  using Rows = std::vector<bool>;
  auto submit = [&](Candidate candidate, const Rows *to_rows, const Rows *from_rows,
                    const Rows *connection_rows) {
    const auto changed = [](const Rows *rows, std::size_t i) { return !rows || (*rows)[i]; };
    const auto &basis = candidate.transform;
    candidate.row_keys.resize(d); candidate.row_characters.resize(d);
    for (std::size_t i=0;i<d;++i) if (changed(to_rows,i)) candidate.row_keys[i] = row_key(basis.to_basis[i]);
    std::string key;
    for (const auto &part : candidate.row_keys) key += part;
    if (!seen.insert(std::move(key)).second) { ++result.duplicates; return; }
    for (std::size_t i=0;i<d;++i)
      if ((changed(to_rows,i) && !row_admissible(basis.to_basis[i])) ||
          (changed(from_rows,i) && !row_admissible(basis.from_basis[i])) ||
          (changed(connection_rows,i) && !row_admissible(basis.connection[i]))) {
        ++result.stats.new_poles; return;
      }
    candidate.metadata.connection_characters = 0;
    for (std::size_t i=0;i<d;++i) {
      if (changed(connection_rows,i)) candidate.row_characters[i] = row_size(basis.connection[i]);
      candidate.metadata.connection_characters += candidate.row_characters[i];
    }
    ++result.stats.accepted;
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
  Candidate root{original,{"original",0},std::nullopt,{},{}};
  for (std::size_t i=0;i<d;++i) {
    root.row_keys.push_back(row_key(unit[i]));
    root.row_characters.push_back(row_size(level.matrix[i]));
  }
  submit(root,nullptr,nullptr,nullptr);
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
    submit(Candidate{std::move(basis),{"epsilon diagonal gauge",0},std::nullopt,{},{}},nullptr,nullptr,nullptr);
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
  // J'=H J, where H is the identity except for row pos, which is q. Then
  // to_basis changes in that row only, and the right factor H^-1 (the identity
  // except for row pos) adds a multiple of column pos to the other columns:
  // every product below is a rank-one update, O(d^2) exact operations.
  auto trial = [&](const Candidate &base,std::size_t pos,
                   const std::vector<Exact> &q,const std::string &why) {
    if(expired() || result.stats.trials>=options.max_trials) return;
    ++result.stats.trials;
    if(q[pos].is_zero()) { ++result.stats.singular; return; }
    Candidate next{base.transform,{base.metadata.provenance+"; "+why,base.metadata.depth+1},
                   std::nullopt,base.row_keys,base.row_characters};
    auto &basis=next.transform;
    const auto row_times=[&](const Matrix &m) {
      std::vector<Exact> out(d,z.constant(0));
      for(std::size_t k=0;k<d;++k) if(!q[k].is_zero())
        for(std::size_t j=0;j<d;++j) if(!m[k][j].is_zero()) out[j]=out[j]+q[k]*m[k][j];
      return out;
    };
    std::vector<Exact> inverse_row(d,z.constant(0));
    for(std::size_t k=0;k<d;++k) inverse_row[k]=k==pos ? z.constant(1)/q[pos] : -q[k]/q[pos];
    const auto times_inverse=[&](Matrix &m,Rows &touched) {
      for(std::size_t i=0;i<d;++i) {
        const Exact pivot=m[i][pos];
        if(pivot.is_zero()) continue;
        touched[i]=true;
        for(std::size_t k=0;k<d;++k)
          if(k==pos) m[i][k]=pivot*inverse_row[k];
          else if(!inverse_row[k].is_zero()) m[i][k]=m[i][k]+pivot*inverse_row[k];
      }
    };
    Rows to_rows(d,false),from_rows(d,false),connection_rows(d,false);
    to_rows[pos]=connection_rows[pos]=true;
    auto connection_row=row_times(base.transform.connection);
    for(std::size_t k=0;k<d;++k) connection_row[k]=connection_row[k]+q[k].derivative(parameter);
    basis.to_basis[pos]=row_times(base.transform.to_basis);
    basis.connection[pos]=std::move(connection_row);
    times_inverse(basis.connection,connection_rows);
    times_inverse(basis.from_basis,from_rows);
    submit(std::move(next),&to_rows,&from_rows,&connection_rows);
  };
  // Each depth is expanded from a fixed snapshot, so accepting a candidate
  // cannot invalidate an in-flight base. Later depths permit paired changes.
  std::vector<Candidate> bases{root};
  for(std::size_t depth=0;depth<options.max_depth && !expired();++depth) {
    const auto depth_limit = std::min(options.max_trials,
        (depth+1)*(options.max_trials/options.max_depth) + options.max_trials%options.max_depth);
    for(const auto &base:bases) {
      auto available = [&] { return !expired() && result.stats.trials<depth_limit; };
      // Pure rescalings J_i=eps^p Y_i are not tried: structural child windows
      // are invariant under every diagonal epsilon gauge, and the connection
      // only trades epsilon powers between entries.
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
