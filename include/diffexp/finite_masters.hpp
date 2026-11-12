#pragma once
// Masters that are finite at epsilon = 0, where a level has such integrals.
//
// A level hands the level below one weight per master, and the requested
// integral is the sum of weight times master over them. If the masters have
// poles in epsilon that the requested integrals do not have, that sum
// cancels: in a sample every cancelled order costs log2(1/epsilon) bits, in
// an expanded run every one is an epsilon order that all levels below carry
// (docs/performance.md, "Cancellation between the levels"). The masters the
// reduction selects are the simplest integrals of each sector, and those are
// the most divergent ones: the four-loop banana in two dimensions is finite,
// its first level integrates the merged line squared, and the masters with
// the first power of that line have four poles.
//
// What is decided here from exact data: whether the parametric representation
// of an integral without numerators,
//
//   Gamma(a - L d/2) Int prod alpha_i^(a_i-1) U^(a-(L+1)d/2) F^(-(a-L d/2)),
//
// converges at epsilon = 0 under every scaling in which a subset of the
// parameters tends to zero (equivalently, the complement to infinity). These
// are the ultraviolet subdivergences and the infrared ones of sets of lines;
// a divergence that needs different powers for different lines (collinear
// regions) is not seen, so the test is necessary for finiteness in general
// and sufficient where no such region exists (Euclidean massive families).
// A master that fails is exchanged for an integral of its sector with more
// dots that passes, if the reduction of that integral has a coefficient at
// the master. Sectors without such an integral keep their masters.
#include "diffexp/moving_poles.hpp"
#include "diffexp/scaleless_sectors.hpp"
#include <chrono>

namespace diffexp::finite_masters {
// The convergence test of one family (a level's merged family).
class Convergence {
 public:
  Convergence(const ibp::PropagatorBasis &basis, std::size_t loops, Rational dimension_at_zero)
      : lines_(basis.physical_count), loops_(static_cast<long>(loops)), dimension_(std::move(dimension_at_zero)),
        symanzik_(scaleless::symanzik(basis)) {}
  std::size_t lines() const { return lines_; }
  // The smallest sum of powers with a finite overall Gamma function.
  long minimal_total() const {
    const Rational bound = Rational(loops_) * dimension_ / Rational(2);
    long total = 1;
    while (!(Rational(total) > bound))
      ++total;
    return total;
  }
  bool operator()(const ibp::Integral &integral) const {
    if (integral.size() < lines_)
      return false;
    std::size_t mask = 0;
    long total = 0;
    for (std::size_t i = 0; i < integral.size(); ++i) {
      if (integral[i] < 0 || (i >= lines_ && integral[i] != 0))
        return false; // numerators are not classified
      if (integral[i] > 0) {
        mask |= std::size_t(1) << i;
        total += integral[i];
      }
    }
    if (!mask)
      return false;
    const Rational q = Rational(total) - Rational(loops_) * dimension_ / Rational(2);
    const Rational p = q - dimension_ / Rational(2);
    if (!(q > Rational(0)))
      return false;
    // Monomials of U and F on the face where the absent lines' parameters vanish.
    const auto &face = faces(mask);
    if (face.first.empty() || face.second.empty())
      return false;
    for (std::size_t subset = (mask - 1) & mask; subset; subset = (subset - 1) & mask) {
      long powers = 0;
      for (std::size_t i = 0; i < lines_; ++i)
        if (subset >> i & 1)
          powers += integral[i];
      const auto lowest = [&](const std::vector<const scaleless::detail::Exponent *> &monomials) {
        long least = std::numeric_limits<long>::max();
        for (const auto *e : monomials) {
          long degree = 0;
          for (std::size_t i = 0; i < lines_; ++i)
            if (subset >> i & 1)
              degree += (*e)[i];
          least = std::min(least, degree);
        }
        return least;
      };
      const Rational exponent = Rational(powers) + p * Rational(lowest(face.first)) - q * Rational(lowest(face.second));
      if (!(exponent > Rational(0)))
        return false;
    }
    return true;
  }

 private:
  using Monomials = std::vector<const scaleless::detail::Exponent *>;
  const std::pair<Monomials, Monomials> &faces(std::size_t mask) const {
    auto [found, inserted] = faces_.try_emplace(mask);
    if (inserted) {
      const auto collect = [&](const scaleless::detail::Polynomial &polynomial, Monomials &out) {
        for (const auto &[e, c] : polynomial) {
          bool inside = true;
          for (std::size_t i = 0; i < lines_ && inside; ++i)
            inside = !e[i] || (mask >> i & 1);
          if (inside)
            out.push_back(&e);
        }
      };
      collect(symanzik_.U, found->second.first);
      collect(symanzik_.F, found->second.second);
    }
    return found->second;
  }
  std::size_t lines_;
  long loops_;
  Rational dimension_;
  scaleless::Symanzik symanzik_;
  mutable std::map<std::size_t, std::pair<Monomials, Monomials>> faces_;
};

inline std::vector<bool> sector(const ibp::Integral &a) {
  std::vector<bool> out;
  for (int power : a)
    out.push_back(power > 0);
  return out;
}
inline long total_power(const ibp::Integral &a) {
  long total = 0;
  for (int power : a)
    if (power > 0)
      total += power;
  return total;
}

// The convergent integrals without numerators of the sector of `pattern`
// whose powers add up to `total`, in the order of the reduction.
inline std::vector<ibp::Integral> convergent_integrals(const Convergence &convergent, const ibp::Integral &pattern, long total) {
  std::vector<std::size_t> present;
  for (std::size_t i = 0; i < convergent.lines() && i < pattern.size(); ++i)
    if (pattern[i] > 0)
      present.push_back(i);
  std::vector<ibp::Integral> out;
  if (present.empty() || total < static_cast<long>(present.size()))
    return out;
  ibp::Integral current(pattern.size(), 0);
  for (auto i : present)
    current[i] = 1;
  const std::function<void(std::size_t, long)> place = [&](std::size_t at, long left) {
    if (at + 1 == present.size()) {
      current[present[at]] = static_cast<int>(1 + left);
      if (convergent(current))
        out.push_back(current);
      current[present[at]] = 1;
      return;
    }
    for (long here = 0; here <= left; ++here) {
      current[present[at]] = static_cast<int>(1 + here);
      place(at + 1, left - here);
    }
    current[present[at]] = 1;
  };
  place(0, total - static_cast<long>(present.size()));
  std::sort(out.begin(), out.end(), ibp::IntegralOrder{});
  return out;
}

struct Step {
  std::size_t replaced;
  ibp::Integral from, to;
};
struct Report {
  std::size_t masters = 0, divergent_before = 0, divergent_after = 0, candidates = 0, reductions = 0;
  std::size_t characters_before = 0, characters_after = 0;
  std::vector<Step> steps;
  double seconds = 0, reduction_seconds = 0;
  std::string note;
};

// Exchange the masters that fail the test for integrals that pass it, among
// the integrals whose reduction the closure holds (its requested integrals
// and source identities). The closure of the result is the given one
// transformed by an exact, verified change of basis.
inline level::Result exchange(level::Result level, std::size_t xi, const Convergence &convergent, Report &report) {
  if (!level.success || level.matrix.empty())
    return level;
  const auto d = level.ordered_basis.size();
  const auto &z = level.matrix[0][0];
  moving_poles::Matrix unit(d, std::vector<Exact>(d, z.constant(0)));
  for (std::size_t i = 0; i < d; ++i)
    unit[i][i] = z.constant(1);
  recursion::PhysicalBasisTransform basis{level.matrix, unit, unit, level.matrix};
  std::vector<bool> open(d, false);
  std::size_t left = 0;
  for (std::size_t k = 0; k < d; ++k)
    if (!convergent(level.ordered_basis[k])) {
      open[k] = true;
      ++left;
    }
  if (!left)
    return level;
  auto rows = moving_poles::cached_integrals(level);
  std::erase_if(rows, [&](const moving_poles::Row &row) { return !convergent(row.integral); });
  std::sort(rows.begin(), rows.end(),
            [](const moving_poles::Row &a, const moving_poles::Row &b) { return ibp::IntegralOrder{}(a.integral, b.integral); });
  auto masters = level.ordered_basis;
  for (const auto &row : rows) {
    if (!left)
      break;
    const auto shape = sector(row.integral);
    bool wanted = false;
    for (std::size_t k = 0; k < d && !wanted; ++k)
      wanted = open[k] && sector(masters[k]) == shape;
    if (!wanted)
      continue;
    // The integral in the current masters.
    std::vector<Exact> q(d, z.constant(0));
    for (std::size_t k = 0; k < d; ++k)
      if (!row.coordinates[k].is_zero())
        for (std::size_t j = 0; j < d; ++j)
          if (!basis.from_basis[k][j].is_zero())
            q[j] = q[j] + row.coordinates[k] * basis.from_basis[k][j];
    for (std::size_t k = d; k-- > 0;)
      if (open[k] && !q[k].is_zero() && sector(masters[k]) == shape) {
        moving_poles::replace_row(basis, k, q, xi);
        report.steps.push_back({k, masters[k], row.integral});
        masters[k] = row.integral;
        open[k] = false;
        --left;
        break;
      }
  }
  if (report.steps.empty())
    return level;
  recursion::verify_physical_basis(basis, level.matrix, xi);
  if (!level.target_rows.empty())
    level.target_rows = fuchsify::detail::multiply(level.target_rows, basis.from_basis);
  level.matrix = basis.connection;
  level.ordered_basis = masters;
  level.differential_witnesses.clear();
  level.target_witnesses.clear();
  return level;
}

// The whole choice for one level. `reduce` returns the closure of the level
// for the given requests (the level's sources first); it is asked for the
// sources and the candidate integrals of one more total power per round,
// until no failing master is left, a round exchanges nothing, or a sector has
// no passing integral within `lines` powers of the first total that could
// have one. The returned closure has rows for the sources only.
inline level::Result choose(const ibp::PropagatorBasis &merged, std::size_t loops, const Rational &dimension_at_zero,
                            level::Result closure, const std::vector<ibp::Integral> &sources, std::size_t xi,
                            const std::function<level::Result(const std::vector<ibp::Integral> &)> &reduce, Report &report) {
  const auto start = std::chrono::steady_clock::now();
  const auto elapsed = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(); };
  report.masters = closure.ordered_basis.size();
  if (!closure.success || closure.matrix.empty() || merged.physical_count >= 8 * sizeof(std::size_t))
    return closure;
  const Convergence convergent(merged, loops, dimension_at_zero);
  const auto failing = [&](const level::Result &level) {
    std::size_t count = 0;
    for (const auto &master : level.ordered_basis)
      count += !convergent(master);
    return count;
  };
  report.divergent_before = report.divergent_after = failing(closure);
  report.characters_before = report.characters_after = moving_poles::characters(closure.matrix);
  if (!report.divergent_before)
    return closure;
  const auto trim = [&](level::Result level) {
    level.requested_integrals.resize(sources.size());
    level.target_rows.resize(sources.size());
    if (level.target_witnesses.size() > sources.size())
      level.target_witnesses.resize(sources.size());
    return level;
  };
  // Per sector with a failing master: the next total power to enumerate and
  // the last one worth trying.
  struct Search {
    ibp::Integral pattern;
    long next, last;
    bool found = false;
  };
  std::map<std::vector<bool>, Search> searches;
  for (const auto &master : closure.ordered_basis)
    if (!convergent(master)) {
      const auto shape = sector(master);
      if (searches.count(shape))
        continue;
      long lines = 0;
      for (bool present : shape)
        lines += present;
      const long first = std::max(lines, convergent.minimal_total());
      searches.emplace(shape, Search{master, first, first + lines});
    }
  const std::set<ibp::Integral> masters(closure.ordered_basis.begin(), closure.ordered_basis.end());
  std::vector<ibp::Integral> requests = sources;
  std::set<ibp::Integral> requested(sources.begin(), sources.end());
  level::Result best = closure;
  Report best_report = report;
  const auto attempt = [&](const level::Result &from) {
    Report scratch;
    auto now = exchange(from, xi, convergent, scratch);
    if (scratch.steps.size() <= best_report.steps.size())
      return;
    best_report.steps = std::move(scratch.steps);
    best_report.divergent_after = failing(now);
    best_report.characters_after = moving_poles::characters(now.matrix);
    best = std::move(now);
  };
  // First with the integrals the closure already holds.
  attempt(closure);
  for (;;) {
    // One more total power in every sector that still has a failing master.
    std::set<std::vector<bool>> wanted;
    for (const auto &master : best.ordered_basis)
      if (!convergent(master))
        wanted.insert(sector(master));
    std::vector<ibp::Integral> added;
    for (auto &[shape, search] : searches) {
      if (!wanted.count(shape))
        continue;
      while (search.next <= search.last) {
        const auto total = search.next++;
        const auto found = convergent_integrals(convergent, search.pattern, total);
        if (found.empty())
          continue;
        if (!search.found) {
          search.found = true;
          search.last = std::max(search.last, total + static_cast<long>(std::count(shape.begin(), shape.end(), true)));
        }
        bool fresh = false;
        for (const auto &integral : found)
          if (!masters.count(integral) && requested.insert(integral).second) {
            added.push_back(integral);
            fresh = true;
          }
        if (fresh)
          break;
      }
    }
    if (added.empty())
      break;
    requests.insert(requests.end(), added.begin(), added.end());
    report.candidates += added.size();
    const auto before = elapsed();
    level::Result wider;
    try {
      wider = reduce(requests);
    } catch (const std::exception &failure) {
      report.note = std::string("reduction of the candidate integrals failed: ") + failure.what();
      break;
    }
    report.reduction_seconds += elapsed() - before;
    ++report.reductions;
    if (!wider.success) {
      report.note = "reduction of the candidate integrals failed: " + wider.reason;
      break;
    }
    if (wider.ordered_basis != closure.ordered_basis) {
      report.note = "the reduction of the candidate integrals selected other masters";
      break;
    }
    attempt(wider);
  }
  report.steps = best_report.steps;
  report.divergent_after = best_report.divergent_after;
  report.characters_after = best_report.characters_after;
  report.seconds = elapsed();
  return trim(std::move(best));
}
} // namespace diffexp::finite_masters
