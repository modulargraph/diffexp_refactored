#pragma once
#include "diffexp/adjoint_transport.hpp"
#include "diffexp/diagonal_log_gauge.hpp"
#include "diffexp/transport.hpp"

namespace diffexp::ft_spectral {
using B = Jet::Ball;
namespace sp = transport::spectral_detail;
struct Options {
  bool diagonal_gauge = true, endpoint_clustering = false, conservative = false;
  unsigned accuracy_goal = 40, max_nodes = 128, max_block_size = 4,
           max_block_nodes = 256;
  std::size_t max_cells = 2000000;
  double seconds_budget = 15;
  unsigned max_subdivisions = 0;
};
struct Diagnostics {
  std::vector<unsigned> block_sizes, nodes;
  unsigned normalized_diagonals = 0, clustered_legs = 0;
  unsigned legs = 0, factorizations = 0, absolute_stability_components = 0;
  double preparation_seconds = 0, numerical_seconds = 0;
  std::string reason, last_failure;
  unsigned subdivisions = 0, accepted_subsegments = 0;
};
namespace detail {
struct Edge {
  unsigned row, column, scalar;
  bool zero_order;
};
struct Source {
  unsigned observable, row, scalar;
};
struct Scalar {
  std::vector<Jet> top, bottom;
  unsigned numerator, denominator;
};
inline std::vector<Jet> compile_polynomial(const Exact &p, unsigned xi,
                                           unsigned ei, unsigned window,
                                           slong bits) {
  const auto terms = p.numerator_terms();
  unsigned xdegree = 0, edegree = 0;
  for (const auto &term : terms) {
    xdegree = std::max(xdegree, static_cast<unsigned>(term.powers[xi]));
    edegree = std::max(edegree, static_cast<unsigned>(std::min<std::uint64_t>(
                                    term.powers[ei], window - 1)));
  }
  std::vector<Jet> coefficients(edegree + 1, Jet(0, xdegree + 1, bits));
  for (const auto &term : terms) {
    if (term.powers[ei] >= window)
      continue;
    B coefficient = B::from_strings(term.coefficient.str());
    for (unsigned v = 0; v < term.powers.size(); ++v)
      if (v != xi && v != ei && term.powers[v]) {
        if (p.variables()[v] != "I")
          throw sp::Rejected("FT spectral unassigned coefficient symbol");
        for (unsigned power = 0; power < term.powers[v] % 4; ++power)
          acb_mul_onei(coefficient.raw(), coefficient.raw());
      }
    auto &polynomial = coefficients[term.powers[ei]];
    polynomial.set(term.powers[xi],
                   polynomial.at(term.powers[xi]) + coefficient);
  }
  return coefficients;
}
inline Jet sample_polynomial(const std::vector<Jet> &coefficients,
                             const B &point, unsigned window, slong bits) {
  Jet result(0, window, bits);
  for (unsigned e = 0; e < coefficients.size(); ++e)
    result.set(e, coefficients[e].evaluate_polynomial(point));
  return result;
}
struct Node {
  B point, coordinate, jacobian;
  std::vector<std::vector<B>> coefficients;
  std::vector<B> gauges, roots;
};
using Key = std::pair<unsigned, unsigned>;
inline std::vector<std::vector<unsigned>>
components(unsigned d, const std::vector<Edge> &edges) {
  std::vector<std::vector<unsigned>> dependencies(d), result;
  for (const auto &edge : edges)
    if (edge.zero_order)
      dependencies[edge.row].push_back(edge.column);
  std::vector<int> index(d, -1), low(d);
  std::vector<unsigned> stack;
  std::vector<bool> active(d);
  int next = 0;
  std::function<void(unsigned)> visit = [&](unsigned v) {
    index[v] = low[v] = next++;
    stack.push_back(v);
    active[v] = true;
    for (auto w : dependencies[v])
      if (index[w] < 0) {
        visit(w);
        low[v] = std::min(low[v], low[w]);
      } else if (active[w])
        low[v] = std::min(low[v], index[w]);
    if (low[v] == index[v]) {
      result.emplace_back();
      for (;;) {
        auto w = stack.back();
        stack.pop_back();
        active[w] = false;
        result.back().push_back(w);
        if (w == v)
          break;
      }
      std::sort(result.back().begin(), result.back().end());
    }
  };
  for (unsigned v = 0; v < d; ++v)
    if (index[v] < 0)
      visit(v);
  return result;
}
inline B scale(const B &value) { return B(1) + transport::magnitude(value); }
// The remainder of this namespace is the part of g'=forcing-g*A transport
// that does not depend on how one refinement level is solved; the collocation
// and coefficient backends share it.
inline void check_request(const ExactEpsilonMatrix &matrix,
                          const LaurentRows &initial,
                          const ExactEpsilonMatrix &forcing,
                          const std::vector<Exact> &vertices,
                          const Options &options, slong bits) {
  const auto d = initial.columns(), r = initial.coefficients.size();
  if (matrix.size() != d || forcing.size() != r || vertices.empty())
    throw sp::Rejected("FT spectral shape mismatch");
  if (!options.accuracy_goal || options.accuracy_goal>100000 || options.max_nodes < 16 ||
      options.max_nodes > 128 ||
      (!std::isfinite(options.seconds_budget) ||
       !(options.seconds_budget > 0)) ||
      (options.accuracy_goal + 2) * 3.321928094887362 > bits - 32)
    throw sp::Rejected("FT spectral accuracy/node/precision budget");
  if (!d || !r || initial.high < initial.low)
    throw sp::Rejected("FT spectral empty or invalid boundary");
  for (const auto &row : initial.coefficients) {
    if (row.size() != d)
      throw sp::Rejected("FT spectral ragged boundary");
    for (const auto &series : row)
      if (series.size() !=
          static_cast<std::size_t>(initial.high - initial.low + 1))
        throw sp::Rejected("FT spectral boundary epsilon shape");
  }
}
// Epsilon window [low,high] (w coefficients), the state carried from leg to
// leg, and the acceptance budgets.
struct Window {
  std::size_t xi = 0, ei = 0;
  int low = 0;
  unsigned w = 0;
  LaurentRows current;
  B tolerance, arithmetic_tolerance, floor;
};
inline Window open_window(const ExactEpsilonMatrix &matrix,
                          const LaurentRows &initial,
                          const ExactEpsilonMatrix &forcing,
                          const std::vector<Exact> &vertices,
                          const Options &options,
                          unsigned arithmetic_accuracy_goal) {
  const auto d = initial.columns(), r = initial.coefficients.size();
  auto [xi, ei] = path_epsilon_variables(vertices[0]);
  int low = std::min(0, initial.low);
  for (const auto &row : matrix)
    if (row.size() != d)
      throw sp::Rejected("FT spectral connection is not square");
  for (const auto &row : forcing) {
    if (row.size() != d)
      throw sp::Rejected("FT spectral forcing shape");
    for (const auto &q : row)
      if (!q.is_zero())
        low =
            std::min(low, static_cast<int>(*exact_epsilon_valuation(q, ei)));
  }
  if (low < -1000 || initial.high - low > 1000)
    throw sp::Rejected("FT spectral epsilon window budget");
  const unsigned w = initial.high - low + 1;
  if(d*r*w>options.max_cells)throw sp::Rejected("FT spectral boundary storage budget");
  Window window{xi, ei, low, w,
                LaurentRows{low, initial.high,
                            std::vector(r, std::vector(d, std::vector<B>(w, B(0))))},
                B(0), B(0), B(1)};
  for (unsigned a = 0; a < r; ++a)
    for (unsigned i = 0; i < d; ++i)
      for (int e = initial.low; e <= initial.high; ++e) {
        const auto &v = initial.coefficients[a][i][e - initial.low];
        if (!v.is_finite())
          throw sp::Rejected("FT spectral nonfinite boundary");
        window.current.coefficients[a][i][e - low] = v;
      }
  window.tolerance =
      B::from_strings("1e-" + std::to_string(options.accuracy_goal + 2));
  // A subdivided arm spends a stricter local discretization budget, but
  // inherited enclosures are tested against the unchanged whole-arm goal.
  // Never discard input radii or classify them as new discretization error.
  window.arithmetic_tolerance=arithmetic_accuracy_goal
      ? B::from_strings("1e-"+std::to_string(arithmetic_accuracy_goal+2))
      : window.tolerance;
  acb_mul_2exp_si(window.floor.raw(), window.floor.raw(), 16 - B::precision());
  return window;
}
// False when some inherited ball is already wider than the arm budget.
// Neither refinement nor subdivision of the path can narrow such a ball.
inline bool inherited_within(const LaurentRows &rows, const B &tolerance) {
  for (const auto &row : rows.coefficients)
    for (const auto &series : row)
      for (const auto &value : series)
        if (!sp::le(transport::arithmetic_error(value),
                    tolerance * sp::midpoint_scale(value)))
          return false;
  return true;
}
// Edges and sources of one straight leg in its unit parameter. intern(q,
// numerator, denominator) returns the backend's scalar id and whether q is
// new; removed_diagonal[i], when set, is the epsilon-zero diagonal already
// absorbed by a gauge. Every denominator is recorded at epsilon zero.
struct Coupling {
  std::vector<Edge> edges;
  std::vector<Source> sources;
  std::vector<Exact> denominators;
  unsigned minimum_degree = 0;
};
template <class Intern>
Coupling couple(const ExactEpsilonMatrix &matrix,
                const ExactEpsilonMatrix &forcing,
                const std::vector<Exact> &point, const Exact &width,
                const Window &window, unsigned max_nodes,
                const std::vector<const Exact *> &removed_diagonal,
                Intern intern) {
  const unsigned d = matrix.size(), r = forcing.size();
  const auto xi = window.xi, ei = window.ei;
  Coupling out;
  const auto scalar = [&](const Exact &q, unsigned numerator,
                          unsigned denominator) {
    for (const auto &terms : {q.numerator_terms(), q.denominator_terms()})
      for (const auto &term : terms) {
        if (term.powers[xi] > max_nodes)
          throw sp::Rejected("FT spectral polynomial degree budget");
        out.minimum_degree = std::max(out.minimum_degree,
                                      static_cast<unsigned>(term.powers[xi]));
      }
    const std::pair<unsigned, bool> interned =
        intern(q, numerator, denominator);
    if (interned.second)
      out.denominators.push_back(
          q.denominator().substitute(exact_point(q, ei, q.constant(0))));
    return interned.first;
  };
  for (unsigned i = 0; i < d; ++i)
    for (unsigned j = 0; j < d; ++j)
      if (!matrix[j][i].is_zero()) {
        out.denominators.push_back(
            matrix[j][i].substitute(point).denominator().substitute(
                exact_point(matrix[j][i], ei, matrix[j][i].constant(0))));
        auto physical = matrix[j][i];
        if (i == j && removed_diagonal[i])
          physical = physical - *removed_diagonal[i];
        if (physical.is_zero())
          continue;
        auto valuation = *exact_epsilon_valuation(physical, ei);
        if (valuation < 0)
          throw sp::Rejected(
              "FT spectral connection needs epsilon normalization");
        if (valuation >= window.w)
          continue;
        auto q = -width * physical.substitute(point);
        out.edges.push_back({i, j, scalar(q, j, i), valuation == 0});
      }
  for (unsigned a = 0; a < r; ++a)
    for (unsigned i = 0; i < d; ++i)
      if (!forcing[a][i].is_zero()) {
        auto q = forcing[a][i];
        q = window.low < 0 ? q * q.variable(ei).pow(-window.low)
                           : q / q.variable(ei).pow(window.low);
        if (*exact_epsilon_valuation(q, ei) < window.w)
          out.sources.push_back(
              {a, i, scalar(width * q.substitute(point), d, i)});
      }
  return out;
}
// Epsilon-zero strongly connected blocks in dependency order.
struct Blocks {
  std::vector<std::vector<unsigned>> members;
  std::vector<unsigned> block_id, position;
  unsigned largest = 0;
};
inline Blocks partition(unsigned d, const std::vector<Edge> &edges,
                        const Options &options, Diagnostics &diagnostics) {
  Blocks out{components(d, edges), std::vector<unsigned>(d),
             std::vector<unsigned>(d)};
  for (unsigned b = 0; b < out.members.size(); ++b) {
    const auto &block = out.members[b];
    diagnostics.block_sizes.push_back(block.size());
    if (block.size() > options.max_block_size)
      throw sp::Rejected("FT spectral epsilon-zero block exceeds budget");
    out.largest = std::max(out.largest, static_cast<unsigned>(block.size()));
    for (unsigned i = 0; i < block.size(); ++i) {
      out.block_id[block[i]] = b;
      out.position[block[i]] = i;
    }
  }
  for (const auto &edge : edges)
    if (edge.zero_order && out.block_id[edge.column] > out.block_id[edge.row])
      throw std::logic_error("FT spectral dependency order");
  return out;
}
// Known rational singularities guard against unresolved narrow features.
struct Geometry {
  std::vector<B> singularities;
  double left_distance = 1, right_distance = 1;
};
template <class Check>
Geometry path_geometry(const std::vector<Exact> &denominators, std::size_t xi,
                       slong bits, Check &check) {
  Geometry geometry;
  std::set<std::string> seen;
  for (auto q : denominators) {
    check();
    if (q.is_zero())
      throw sp::Rejected("FT spectral unresolved epsilon denominator");
    if (q.is_rational())
      continue;
    const auto &names = q.variables();
    auto im = std::find(names.begin(), names.end(), "I");
    if (im != names.end())
      q = polynomial_norm(q, im - names.begin(), q.constant(-1));
    q = q / q.constant(q.numerator_terms().front().coefficient);
    if (!seen.insert(q.str()).second)
      continue;
    auto roots = polynomial_roots(q, xi, bits);
    geometry.singularities.insert(geometry.singularities.end(), roots.begin(),
                                  roots.end());
  }
  B domain = B::from_strings("1/2");
  arb_add_error_2exp_si(acb_realref(domain.raw()), -1);
  for (const auto &pole : geometry.singularities)
    if (arb_contains_zero(acb_imagref(pole.raw())) &&
        arb_overlaps(acb_realref(pole.raw()), acb_realref(domain.raw())))
      throw sp::Rejected(
          "FT spectral path intersects singularity or branch point");
  for (const auto &pole : geometry.singularities) {
    const std::complex<double> z(
        arf_get_d(arb_midref(acb_realref(pole.raw())), ARF_RND_NEAR),
        arf_get_d(arb_midref(acb_imagref(pole.raw())), ARF_RND_NEAR));
    geometry.left_distance = std::min(geometry.left_distance, std::abs(z));
    geometry.right_distance =
        std::min(geometry.right_distance, std::abs(z - std::complex<double>(1)));
  }
  return geometry;
}
// Refinement state of one leg: the endpoint values of the last level and the
// last two differences between levels, in (observable, component, epsilon)
// order. A level is accepted only when every component's estimated tail and
// its arithmetic radius fit their budgets.
struct Refinement {
  std::vector<B> previous, previous_difference;
  unsigned previous_n = 0, previous_increment = 0;
  // Outcome of the last tested level, kept to classify an exhausted ladder.
  bool tested = false, discretization_met = false, arithmetic_excess = false;
  std::optional<std::vector<B>> converged(const std::vector<B> &values,
                                          const std::vector<B> &difference,
                                          unsigned n, const Window &window,
                                          unsigned &absolute) {
    std::vector<B> tails(values.size());
    bool all = true;
    tested = discretization_met = true;
    arithmetic_excess = false;
    absolute = 0;
    for (unsigned i = 0; i < values.size(); ++i) {
      auto roundoff = window.floor * (scale(values[i]) + scale(previous[i]));
      auto estimate =
          sp::estimate_tail(previous_difference[i], difference[i], roundoff,
                            previous_increment, n - previous_n);
      if (estimate.absolute)
        ++absolute;
      tails[i] = std::move(estimate.tail);
      const auto size = scale(values[i]);
      const auto radius = transport::arithmetic_error(values[i]);
      if (!sp::le(B(2)*tails[i],window.tolerance*size))
        all = discretization_met = false;
      if (!sp::le(radius+B(2)*tails[i],window.arithmetic_tolerance*size))
        all = false;
      if (!sp::le(radius, window.arithmetic_tolerance * size))
        arithmetic_excess = true;
    }
    if (!all)
      return std::nullopt;
    return tails;
  }
  void advance(std::vector<B> values, std::vector<B> difference, unsigned n) {
    if (!previous.empty())
      previous_difference = std::move(difference);
    previous = std::move(values);
    previous_increment = n - previous_n;
    previous_n = n;
  }
};
// Accepted values absorb their estimated tails and become the next state.
inline void commit(LaurentRows &current, std::vector<B> &values,
                   const std::vector<B> &tails) {
  std::size_t id = 0;
  for (auto &row : current.coefficients)
    for (auto &series : row)
      for (auto &coefficient : series) {
        arb_add_error(acb_realref(values[id].raw()),
                      acb_realref(tails[id].raw()));
        arb_add_error(acb_imagref(values[id].raw()),
                      acb_realref(tails[id].raw()));
        coefficient = std::move(values[id]);
        ++id;
      }
}
// The inverse of the scalar collocation operator depends only on the node
// count and the working precision, not on the leg or the system.
using PlainInverses =
    std::map<std::pair<unsigned, slong>, std::shared_ptr<sp::Matrix>>;
// Subdivision changes only the parametrization of each original straight leg.
// Failed attempts never update the accepted state. No pole/structural rejection
// is retried; every child goes through the complete original geometry checks.
// single(state, path, options, diagnostics, arithmetic_goal) transports one
// path; arithmetic_goal is zero for an unsubdivided arm. progress(leg,
// fraction) reports the accepted part of an original leg.
template <class BackendOptions, class Single, class Progress>
std::optional<LaurentRows> subdivide(const LaurentRows &initial,
                                     const std::vector<Exact> &vertices,
                                     const BackendOptions &options,
                                     Diagnostics &diagnostics, Single single,
                                     Progress progress) {
  if(!options.max_subdivisions)
    return single(initial,vertices,options,diagnostics,0u);
  diagnostics={};
  if(options.max_subdivisions>4096 || vertices.empty() ||
     !std::isfinite(options.seconds_budget) || options.seconds_budget<=0 ||
     !options.accuracy_goal || options.accuracy_goal>100000) {
    diagnostics.reason="FT spectral subdivision/options budget";return std::nullopt;
  }
  if(std::adjacent_find(vertices.begin(),vertices.end(),[](const Exact& a,const Exact& b){return a!=b;})==vertices.end())
    return single(initial,vertices,options,diagnostics,0u);
  const auto started=std::chrono::steady_clock::now();
  auto elapsed=[&]{return std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();};
  auto fail=[&](const std::string& reason)->std::optional<LaurentRows>{diagnostics.reason=reason;return std::nullopt;};
  BackendOptions local=options;local.max_subdivisions=0;
  // A decimal reserve distributes a conservative local accuracy target over
  // the maximum number of accepted segments; propagated balls remain intact.
  std::size_t segments=vertices.size()-1+options.max_subdivisions;
  unsigned reserve=1;for(auto n=segments;n>1;n=(n+9)/10)++reserve;
  local.accuracy_goal+=reserve;
  LaurentRows state=initial;
  try {
    for(std::size_t leg=0;leg+1<vertices.size();++leg) {
      if(vertices[leg]==vertices[leg+1])continue;
      struct Segment{Exact first,second;double from,to;};
      std::vector<Segment> pending{{vertices[leg],vertices[leg+1],0.,1.}};
      while(!pending.empty()) {
        local.seconds_budget=options.seconds_budget-elapsed();
        if(local.seconds_budget<=0)return fail("FT spectral total subdivision time budget exhausted");
        auto segment=std::move(pending.back());pending.pop_back();Diagnostics attempt;
        auto value=single(state,std::vector<Exact>{segment.first,segment.second},local,attempt,options.accuracy_goal);
        diagnostics.preparation_seconds+=attempt.preparation_seconds;
        diagnostics.numerical_seconds+=attempt.numerical_seconds;
        diagnostics.factorizations+=attempt.factorizations;
        if(elapsed()>options.seconds_budget)return fail("FT spectral total subdivision time budget exhausted");
        if(value) {
          state=std::move(*value);++diagnostics.accepted_subsegments;
          progress(static_cast<unsigned>(leg),segment.to);
          diagnostics.nodes.insert(diagnostics.nodes.end(),attempt.nodes.begin(),attempt.nodes.end());
          diagnostics.block_sizes.insert(diagnostics.block_sizes.end(),attempt.block_sizes.begin(),attempt.block_sizes.end());
          diagnostics.normalized_diagonals=std::max(diagnostics.normalized_diagonals,attempt.normalized_diagonals);
          diagnostics.clustered_legs+=attempt.clustered_legs;
          diagnostics.absolute_stability_components+=attempt.absolute_stability_components;
          continue;
        }
        diagnostics.last_failure=attempt.reason;
        const bool retry=attempt.reason=="FT spectral resolution or arithmetic accuracy exhausted" ||
                         attempt.reason=="FT spectral interval too close to singularities";
        if(!retry)return fail(attempt.reason);
        if(diagnostics.subdivisions>=options.max_subdivisions)
          return fail("FT spectral subdivision count exhausted: "+attempt.reason);
        const auto midpoint=(segment.first+segment.second)/segment.first.constant(2);
        ++diagnostics.subdivisions;
        const auto fraction=(segment.from+segment.to)/2;
        pending.push_back({midpoint,segment.second,fraction,segment.to});
        pending.push_back({segment.first,midpoint,segment.from,fraction});
      }
      ++diagnostics.legs;
    }
    return state;
  }catch(const std::exception& error){return fail(error.what());}
}
} // namespace detail
// Solve g'=forcing-g*A. Epsilon-zero strongly connected blocks are factored
// once per resolution; epsilon-positive terms and preceding blocks form known
// right-hand sides. The same operators act on every observable row.
inline std::optional<LaurentRows>
try_transport_single(const ExactEpsilonMatrix &matrix, const LaurentRows &initial,
              const ExactEpsilonMatrix &forcing,
              const std::vector<Exact> &vertices, const Options &options,
              Diagnostics &diagnostics, unsigned arithmetic_accuracy_goal = 0,
              detail::PlainInverses *shared_inverses = nullptr) {
  using namespace detail;
  diagnostics = {};
  auto started = std::chrono::steady_clock::now();
  double preparing = 0;
  auto elapsed = [&] {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         started)
        .count();
  };
  auto check = [&] {
    if (elapsed() > options.seconds_budget)
      throw sp::Rejected("FT spectral time budget exhausted");
  };
  try {
    const auto bits = B::precision();
    check_request(matrix, initial, forcing, vertices, options, bits);
    const auto d = initial.columns(), r = initial.coefficients.size();
    auto window = open_window(matrix, initial, forcing, vertices, options,
                              arithmetic_accuracy_goal);
    const auto xi = window.xi, ei = window.ei;
    const unsigned w = window.w;
    auto &current = window.current;
    std::vector<DiagonalLogGaugeResult> gauges(d);
    std::vector<Exact> diagonal_zero(d, vertices[0].constant(0));
    if (options.diagonal_gauge)
      for (unsigned i = 0; i < d; ++i) {
        diagonal_zero[i] = feynman::scalar_functional_detail::epsilon_series(
                               matrix[i][i], ei, 0)
                               .at(0);
        gauges[i] = diagonal_log_gauge(-diagonal_zero[i], xi);
        if (gauges[i].supported)
          for (const auto &factor : gauges[i].factors) {
            auto doubled = factor.exponent * Rational(2);
            if (doubled > Rational(64) || doubled < Rational(-64)) {
              gauges[i].supported = false;
              break;
            }
          }
        if (gauges[i].supported && !gauges[i].factors.empty())
          ++diagnostics.normalized_diagonals;
      }
    PlainInverses local_inverses;
    auto &plain_inverses = shared_inverses ? *shared_inverses : local_inverses;
    for (unsigned leg = 0; leg + 1 < vertices.size(); ++leg) {
      const auto prep = std::chrono::steady_clock::now();
      check();
      const auto &from = vertices[leg];
      const auto width = vertices[leg + 1] - from;
      if (width.is_zero())
        continue;
      auto point = exact_point(from, xi, from + width * from.variable(xi));
      const bool inherited_excess =
          !inherited_within(current, window.arithmetic_tolerance);
      std::vector<Scalar> scalars;
      std::map<std::string, unsigned> ids;
      std::vector<const Exact *> removed_diagonal(d, nullptr);
      for (unsigned i = 0; i < d; ++i)
        if (gauges[i].supported)
          removed_diagonal[i] = &diagonal_zero[i];
      auto coupling = couple(
          matrix, forcing, point, width, window, options.max_nodes,
          removed_diagonal,
          [&](const Exact &q, unsigned numerator, unsigned denominator) {
            auto text = q.str();
            auto key = text + "#" + std::to_string(numerator) + "/" +
                       std::to_string(denominator);
            auto [it, added] = ids.try_emplace(key, scalars.size());
            if (added)
              scalars.push_back(
                  {compile_polynomial(q.numerator(), xi, ei, w, bits),
                   compile_polynomial(q.denominator(), xi, ei, w, bits),
                   numerator, denominator});
            return std::pair<unsigned, bool>{it->second, added};
          });
      const auto &edges = coupling.edges;
      const auto &sources = coupling.sources;
      auto &denominators = coupling.denominators;
      const unsigned minimum_degree = coupling.minimum_degree;
      std::vector<Jet> gauge_polynomials;
      std::vector<B> initial_roots;
      std::vector<bool> root_needed;
      std::vector<std::vector<std::pair<unsigned, slong>>> gauge_factors(d);
      std::map<std::string, unsigned> factor_ids;
      for (unsigned i = 0; i < d; ++i)
        if (gauges[i].supported)
          for (const auto &factor : gauges[i].factors) {
            auto polynomial = factor.polynomial.substitute(point);
            denominators.push_back(polynomial);
            auto [it, added] = factor_ids.try_emplace(polynomial.str(),
                                                      gauge_polynomials.size());
            if (added) {
              unsigned degree = 0;
              for (const auto &term : polynomial.numerator_terms())
                degree =
                    std::max(degree, static_cast<unsigned>(term.powers[xi]));
              Jet x(0, degree + 2, bits);
              x.set(1, B(1));
              auto jet = evaluate(data::Reader(polynomial.str()).read(), x,
                                  {{"x", x}});
              B root;
              auto value = jet.at(0);
              acb_sqrt(root.raw(), value.raw(), bits);
              if (!root.is_finite() || root.contains_zero())
                throw sp::Rejected(
                    "FT spectral diagonal gauge singular at initial point");
              gauge_polynomials.push_back(std::move(jet));
              initial_roots.push_back(std::move(root));
              root_needed.push_back(false);
            }
            auto power = std::stol((factor.exponent * Rational(2)).str());
            gauge_factors[i].push_back({it->second, power});
            if (power % 2)
              root_needed[it->second] = true;
          }
      const auto partitioned = partition(d, edges, options, diagnostics);
      const auto &blocks = partitioned.members;
      const auto &block_id = partitioned.block_id;
      const auto &position = partitioned.position;
      std::vector<std::vector<unsigned>> incoming(d), inhomogeneous(d);
      for (unsigned e = 0; e < edges.size(); ++e)
        incoming[edges[e].row].push_back(e);
      for (unsigned s = 0; s < sources.size(); ++s)
        inhomogeneous[sources[s].row].push_back(s);
      const auto geometry = path_geometry(denominators, xi, bits, check);
      const double left_distance = geometry.left_distance,
                   right_distance = geometry.right_distance;
      const bool cluster = options.endpoint_clustering &&
                           std::min(left_distance, right_distance) < 0.1;
      const bool cluster_right = right_distance < left_distance;
      B strength(0), exponential_range(1);
      if (cluster) {
        ++diagnostics.clustered_legs;
        acb_set_d(strength.raw(),
                  std::min(8., std::log1p(1 / std::min(left_distance,
                                                       right_distance))));
        acb_expm1(exponential_range.raw(), strength.raw(), bits);
      }
      if (!cluster && sp::node_forecast(geometry.singularities, options.accuracy_goal) >
                          (options.conservative ? 2 * options.max_nodes : 512))
        throw sp::Rejected("FT spectral interval too close to singularities");
      std::map<Key, Node, sp::NodeLess> cache;
      auto node = [&](unsigned j, unsigned n) -> const Node & {
        unsigned divisor = std::gcd(j, n);
        Key key{j / divisor, n / divisor};
        if (auto it = cache.find(key); it != cache.end())
          return it->second;
        check();
        if ((cache.size() + 1) *
                (scalars.size() * w + d + gauge_polynomials.size()) >
            options.max_cells)
          throw sp::Rejected("FT spectral sample storage budget");
        B angle = B::from_strings(std::to_string(key.first) + "/" +
                                  std::to_string(key.second)),
          cosine;
        acb_cos_pi(cosine.raw(), angle.raw(), bits);
        Node out{(B(1) - cosine) / B(2), (B(1) - cosine) / B(2), B(1), {}};
        if (cluster) {
          auto argument = strength * (cluster_right ? B(1) - out.coordinate
                                                    : out.coordinate);
          B exponential;
          acb_exp(exponential.raw(), argument.raw(), bits);
          out.jacobian = strength * exponential / exponential_range;
          acb_expm1(exponential.raw(), argument.raw(), bits);
          out.point = exponential / exponential_range;
          if (cluster_right)
            out.point = B(1) - out.point;
          if (j == 0)
            out.point = B(0);
          if (j == n)
            out.point = B(1);
        }
        out.gauges.assign(d, B(1));
        if (j) {
          auto previous = std::prev(cache.lower_bound(key));
          for (unsigned i = 0; i < gauge_polynomials.size(); ++i)
            out.roots.push_back(
                root_needed[i]
                    ? sp::continue_root(gauge_polynomials[i],
                                        previous->second.point, out.point,
                                        previous->second.roots[i], check)
                    : B(1));
          for (unsigned i = 0; i < d; ++i)
            for (auto [factor, power] : gauge_factors[i]) {
              B ratio;
              slong exponent = power;
              if (power % 2)
                ratio = out.roots[factor] / initial_roots[factor];
              else {
                ratio =
                    gauge_polynomials[factor].evaluate_polynomial(out.point) /
                    gauge_polynomials[factor].at(0);
                exponent /= 2;
              }
              B value;
              acb_pow_si(value.raw(), ratio.raw(), exponent, bits);
              out.gauges[i] = out.gauges[i] * value;
            }
        } else
          out.roots = initial_roots;
        if (j) {
          for (const auto &scalar : scalars) {
            auto q = sample_polynomial(scalar.top, out.point, w, bits) /
                     sample_polynomial(scalar.bottom, out.point, w, bits);
            B ratio(1);
            if (scalar.numerator != scalar.denominator)
              ratio =
                  (scalar.numerator < d ? out.gauges[scalar.numerator] : B(1)) /
                  (scalar.denominator < d ? out.gauges[scalar.denominator]
                                          : B(1));
            out.coefficients.emplace_back();
            for (unsigned e = 0; e < w; ++e) {
              auto value = q.at(e) * ratio * out.jacobian;
              if (!value.is_finite())
                throw sp::Rejected("FT spectral singular sample");
              out.coefficients.back().push_back(std::move(value));
            }
          }
        }
        return cache.emplace(key, std::move(out)).first->second;
      };
      preparing +=
          std::chrono::duration<double>(std::chrono::steady_clock::now() - prep)
              .count();
      Refinement refinement;
      bool accepted = false;
      for (unsigned n : sp::ladder) {
        if (n > options.max_nodes)
          break;
        check();
        // The budget counts intervals: 4x64 and 2x128 are the intended
        // limits, and the closing node must not exclude them. A level beyond
        // the budget ends refinement like an exhausted ladder, so shorter
        // segments may still converge below it; a block too large for the
        // first level cannot be helped by subdivision.
        if (partitioned.largest * n > options.max_block_nodes) {
          if (!refinement.previous_n)
            throw sp::Rejected("FT spectral coupled-node budget");
          break;
        }
        const unsigned m = n + 1;
        if (static_cast<std::uint64_t>(m) * w * d * r > options.max_cells)
          throw sp::Rejected("FT spectral solution storage budget");
        diagnostics.nodes.push_back(n);
        std::vector<const Node *> nodes(m);
        for (unsigned j = 0; j < m; ++j)
          nodes[j] = &node(j, n);
        sp::Matrix derivative(m, m), state(w * m * d, r);
        std::vector<B> weights(m);
        for (unsigned j = 0; j < m; ++j) {
          weights[j] = B(j % 2 ? -1 : 1);
          if (j == 0 || j == n)
            weights[j] = weights[j] / B(2);
        }
        acb_one(acb_mat_entry(derivative.value, 0, 0));
        for (unsigned j = 1; j < m; ++j) {
          B diagonal(0);
          for (unsigned h = 0; h < m; ++h)
            if (j != h) {
              auto q = weights[h] / weights[j] /
                       (nodes[j]->coordinate - nodes[h]->coordinate);
              acb_set(acb_mat_entry(derivative.value, j, h), q.raw());
              diagonal -= q;
            }
          acb_set(acb_mat_entry(derivative.value, j, j), diagonal.raw());
        }
        std::vector<std::shared_ptr<sp::Matrix>> inverses;
        std::shared_ptr<sp::Matrix> plain = plain_inverses[{n, bits}];
        for (unsigned b = 0; b < blocks.size(); ++b) {
          check();
          unsigned size = blocks[b].size(), z = size * m;
          bool zero = true;
          for (auto i : blocks[b])
            for (auto e : incoming[i])
              if (edges[e].zero_order && block_id[edges[e].column] == b)
                zero = false;
          if (zero && plain) {
            inverses.push_back(plain);
            continue;
          }
          sp::Matrix op(z, z);
          auto inverse = std::make_shared<sp::Matrix>(z, z);
          for (unsigned i = 0; i < size; ++i) {
            acb_one(acb_mat_entry(op.value, i, i));
            for (unsigned j = 1; j < m; ++j) {
              for (unsigned h = 0; h < m; ++h)
                acb_set(acb_mat_entry(op.value, j * size + i, h * size + i),
                        acb_mat_entry(derivative.value, j, h));
              for (auto e : incoming[blocks[b][i]]) {
                const auto &edge = edges[e];
                if (edge.zero_order && block_id[edge.column] == b)
                  acb_sub(acb_mat_entry(op.value, j * size + i,
                                        j * size + position[edge.column]),
                          acb_mat_entry(op.value, j * size + i,
                                        j * size + position[edge.column]),
                          nodes[j]->coefficients[edge.scalar][0].raw(), bits);
              }
            }
          }
          if (!acb_mat_inv(inverse->value, op.value, bits))
            throw sp::Rejected("FT spectral block inverse unresolved");
          ++diagnostics.factorizations;
          check();
          inverses.push_back(inverse);
          if (zero) {
            plain = inverse;
            plain_inverses[{n, bits}] = inverse;
          }
        }
        for (unsigned e = 0; e < w; ++e)
          for (unsigned b = 0; b < blocks.size(); ++b) {
            check();
            const unsigned size = blocks[b].size(), z = size * m;
            sp::Matrix rhs(z, r), solution(z, r);
            for (unsigned i = 0; i < size; ++i) {
              unsigned row = blocks[b][i];
              for (unsigned a = 0; a < r; ++a)
                acb_set(acb_mat_entry(rhs.value, i, a),
                        current.coefficients[a][row][e].raw());
              for (unsigned j = 1; j < m; ++j) {
                for (auto s : inhomogeneous[row]) {
                  const auto &source = sources[s];
                  acb_add(
                      acb_mat_entry(rhs.value, j * size + i, source.observable),
                      acb_mat_entry(rhs.value, j * size + i, source.observable),
                      nodes[j]->coefficients[source.scalar][e].raw(), bits);
                }
                for (auto id : incoming[row]) {
                  const auto &edge = edges[id];
                  unsigned first = block_id[edge.column] == b ? 1 : 0;
                  for (unsigned q = first; q <= e; ++q) {
                    const auto &coefficient =
                        nodes[j]->coefficients[edge.scalar][q];
                    if (coefficient.is_zero())
                      continue;
                    for (unsigned a = 0; a < r; ++a)
                      acb_addmul(
                          acb_mat_entry(rhs.value, j * size + i, a),
                          coefficient.raw(),
                          acb_mat_entry(state.value,
                                        ((e - q) * m + j) * d + edge.column, a),
                          bits);
                  }
                }
              }
            }
            acb_mat_mul(solution.value, inverses[b]->value, rhs.value, bits);
            for (unsigned j = 0; j < m; ++j)
              for (unsigned i = 0; i < size; ++i)
                for (unsigned a = 0; a < r; ++a)
                  acb_set(acb_mat_entry(state.value,
                                        (e * m + j) * d + blocks[b][i], a),
                          acb_mat_entry(solution.value, j * size + i, a));
          }
        std::vector<B> values(r * d * w), difference(r * d * w);
        for (unsigned a = 0; a < r; ++a)
          for (unsigned i = 0; i < d; ++i)
            for (unsigned e = 0; e < w; ++e) {
              auto id = (a * d + i) * w + e;
              acb_mul(values[id].raw(),
                      acb_mat_entry(state.value, (e * m + n) * d + i, a),
                      nodes[n]->gauges[i].raw(), bits);
              if (!values[id].is_finite())
                throw sp::Rejected("FT spectral nonfinite endpoint");
              if (!refinement.previous.empty())
                difference[id] =
                    sp::difference(values[id], refinement.previous[id]);
            }
        if (!refinement.previous_difference.empty() && n > minimum_degree) {
          unsigned absolute = 0;
          if (auto tails = refinement.converged(values, difference, n, window,
                                                absolute)) {
            commit(current, values, *tails);
            accepted = true;
            diagnostics.absolute_stability_components += absolute;
            ++diagnostics.legs;
            break;
          }
        }
        refinement.advance(std::move(values), std::move(difference), n);
      }
      if (!accepted) {
        // Shorter segments reduce discretization error only. When the last
        // level resolved the solution and the radius inherited from the
        // boundary alone exceeds the arm budget, subdividing cannot help.
        if (inherited_excess && refinement.tested &&
            refinement.discretization_met && refinement.arithmetic_excess)
          throw sp::Rejected(
              "FT spectral inherited input uncertainty exceeds arm budget");
        throw sp::Rejected(
            "FT spectral resolution or arithmetic accuracy exhausted");
      }
    }
    diagnostics.preparation_seconds = preparing;
    diagnostics.numerical_seconds = elapsed() - preparing;
    return current;
  } catch (const std::exception &error) {
    diagnostics.reason = error.what();
    diagnostics.preparation_seconds = preparing;
    diagnostics.numerical_seconds = elapsed() - preparing;
    return std::nullopt;
  }
}
// Completed-arm transport with optional bisection of legs; see
// detail::subdivide for the retry contract.
inline std::optional<LaurentRows>
try_transport(const ExactEpsilonMatrix& matrix,const LaurentRows& initial,
              const ExactEpsilonMatrix& forcing,const std::vector<Exact>& vertices,
              const Options& options,Diagnostics& diagnostics) {
  detail::PlainInverses inverses;
  return detail::subdivide(initial,vertices,options,diagnostics,
    [&](const LaurentRows& state,const std::vector<Exact>& path,const Options& local,Diagnostics& attempt,unsigned goal) {
      return try_transport_single(matrix,state,forcing,path,local,attempt,goal,&inverses);
    },[](unsigned,double){});
}
} // namespace diffexp::ft_spectral
