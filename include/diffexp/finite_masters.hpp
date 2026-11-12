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
#include "diffexp/epsilon_gauge.hpp"
#include "diffexp/level_cache.hpp"
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
  bool recorded = false;
  std::string note;
  // Where the choice is recorded (empty without a directory).
  std::filesystem::path record;
};

// Exchange the masters that fail the test for integrals that pass it, among
// the integrals whose reduction the closure holds (its requested integrals
// and source identities). The closure of the result is the given one
// transformed by an exact, verified change of basis, with rows for the first
// `sources` requests. With `script`, exactly those steps are taken (an
// exception if one does not apply). Without `finish_incomplete` a result
// that still has failing masters is not a closure: only its masters are set.
// The source identities of the result name
// the passing integrals in the new masters: they are what the internal basis
// scan of the numerical phase may mix in, and the exchanged masters are not
// offered to it again.
inline level::Result exchange(level::Result level, std::size_t xi, std::size_t ei, const Convergence &convergent,
                              std::vector<Step> &steps, std::size_t sources, const std::vector<Step> *script = nullptr,
                              bool guarded = true, bool finish_incomplete = true) {
  steps.clear();
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
  // The integral in the current masters.
  const auto current = [&](const std::vector<Exact> &coordinates) {
    std::vector<Exact> q(d, z.constant(0));
    for (std::size_t k = 0; k < d; ++k)
      if (!coordinates[k].is_zero())
        for (std::size_t j = 0; j < d; ++j)
          if (!basis.from_basis[k][j].is_zero())
            q[j] = q[j] + coordinates[k] * basis.from_basis[k][j];
    return q;
  };
  // Integrals that are finite at epsilon = 0 can still be linearly dependent
  // there. Their connection then has poles in epsilon that no rescaling of
  // single masters removes, and the expanded transport has no gauge for it:
  // a step that leads to such a connection is not taken.
  const auto sheared = [&](const moving_poles::Matrix &connection) {
    try {
      (void)epsilon_diagonal_gauge(connection, ei);
      return true;
    } catch (const UnsupportedEpsilonGauge &) {
      return false;
    }
  };
  const bool checked = !script && guarded, guard = checked && sheared(level.matrix);
  // Unchecked and all or none: the steps are found first, without the
  // connection, and taken only if they exchange every failing master.
  const bool light = !script && !guarded && !finish_incomplete;
  // Nor is a step taken that gives the connection a singular point that moves
  // with epsilon and that it did not have: finite integrals can be dependent
  // at single points, and a moving pole that meets an endpoint at epsilon = 0
  // leaves the level without an endpoint series (planar double box, lowest
  // level). A new singular point in x alone is admitted: it is an apparent
  // one, which the paths pass like any other (the finite masters of the
  // four-loop banana have one on each of the two upper levels).
  denominator_factors::Catalogue catalogue(xi, ei);
  std::set<std::size_t> singular;
  if (checked)
    for (const auto &matrix_row : level.matrix)
      for (const auto &value : matrix_row)
        if (!value.is_zero())
          for (const auto &[factor, power] : catalogue.split(value.denominator()))
            singular.insert(factor);
  const auto new_moving_pole = [&](const moving_poles::Matrix &connection) {
    for (const auto &matrix_row : connection)
      for (const auto &value : matrix_row)
        if (!value.is_zero())
          for (const auto &[factor, power] : catalogue.split(value.denominator()))
            if (catalogue.factors()[factor].kind == denominator_factors::Kind::mixed && !singular.count(factor))
              return true;
    return false;
  };
  const auto take = [&](std::size_t k, const std::vector<Exact> &q, const ibp::Integral &integral) {
    if (checked) {
      auto next = basis;
      moving_poles::replace_row(next, k, q, xi);
      if ((guard && !sheared(next.connection)) || new_moving_pole(next.connection))
        return false;
      basis = std::move(next);
    } else if (light) {
      // Which steps apply is decided by the old masters in the new ones alone.
      std::vector<Exact> inverse_row(d, z.constant(0));
      for (std::size_t j = 0; j < d; ++j)
        inverse_row[j] = j == k ? z.constant(1) / q[k] : z.constant(0) - q[j] / q[k];
      for (std::size_t i = 0; i < d; ++i) {
        const Exact pivot = basis.from_basis[i][k];
        if (pivot.is_zero())
          continue;
        for (std::size_t j = 0; j < d; ++j)
          if (j == k)
            basis.from_basis[i][j] = pivot * inverse_row[j];
          else if (!inverse_row[j].is_zero())
            basis.from_basis[i][j] = basis.from_basis[i][j] + pivot * inverse_row[j];
      }
    } else
      moving_poles::replace_row(basis, k, q, xi);
    steps.push_back({k, masters[k], integral});
    masters[k] = integral;
    open[k] = false;
    --left;
    return true;
  };
  if (script) {
    for (const auto &step : *script) {
      const auto row = std::find_if(rows.begin(), rows.end(), [&](const moving_poles::Row &r) { return r.integral == step.to; });
      if (row == rows.end() || step.replaced >= d || !open[step.replaced] || masters[step.replaced] != step.from)
        throw std::runtime_error("recorded master exchange does not apply");
      const auto q = current(row->coordinates);
      if (q[step.replaced].is_zero())
        throw std::runtime_error("recorded master exchange is singular");
      take(step.replaced, q, step.to);
    }
  } else
    for (const auto &row : rows) {
      if (!left)
        break;
      const auto shape = sector(row.integral);
      bool wanted = false;
      for (std::size_t k = 0; k < d && !wanted; ++k)
        wanted = open[k] && sector(masters[k]) == shape;
      if (!wanted)
        continue;
      const auto q = current(row.coordinates);
      for (std::size_t k = d; k-- > 0;)
        if (open[k] && !q[k].is_zero() && sector(masters[k]) == shape && take(k, q, row.integral))
          break;
    }
  if (steps.empty())
    return level;
  if (left && !finish_incomplete) {
    // Only the masters are reported: the caller wants all or none.
    level.ordered_basis = masters;
    return level;
  }
  if (light) {
    const auto planned = steps;
    return exchange(std::move(level), xi, ei, convergent, steps, sources, &planned);
  }
  recursion::verify_physical_basis(basis, level.matrix, xi);
  level.requested_integrals.resize(std::min(sources, level.requested_integrals.size()));
  level.target_rows.resize(std::min(sources, level.target_rows.size()));
  if (!level.target_rows.empty())
    level.target_rows = fuchsify::detail::multiply(level.target_rows, basis.from_basis);
  const std::set<ibp::Integral> chosen(masters.begin(), masters.end());
  std::vector<ibp::Relation> identities;
  for (const auto &row : rows) {
    if (chosen.count(row.integral))
      continue;
    const auto q = current(row.coordinates);
    ibp::Relation relation{{row.integral, z.constant(1)}};
    for (std::size_t j = 0; j < d; ++j)
      ibp::add(relation, masters[j], z.constant(0) - q[j]);
    identities.push_back(std::move(relation));
  }
  level.source_identities = std::move(identities);
  level.matrix = basis.connection;
  level.ordered_basis = masters;
  level.differential_witnesses.clear();
  level.target_witnesses.clear();
  return level;
}

// Whether `connection` has only the moving singular points of `original` and
// no epsilon poles beyond a diagonal shear if `original` has none: the
// conditions exchange() imposes per step, on a finished basis.
inline bool regular(const moving_poles::Matrix &original, const moving_poles::Matrix &connection, std::size_t xi, std::size_t ei) {
  const auto sheared = [&](const moving_poles::Matrix &matrix) {
    try {
      (void)epsilon_diagonal_gauge(matrix, ei);
      return true;
    } catch (const UnsupportedEpsilonGauge &) {
      return false;
    }
  };
  if (sheared(original) && !sheared(connection))
    return false;
  denominator_factors::Catalogue catalogue(xi, ei);
  std::set<std::size_t> singular;
  for (const auto &matrix_row : original)
    for (const auto &value : matrix_row)
      if (!value.is_zero())
        for (const auto &[factor, power] : catalogue.split(value.denominator()))
          singular.insert(factor);
  for (const auto &matrix_row : connection)
    for (const auto &value : matrix_row)
      if (!value.is_zero())
        for (const auto &[factor, power] : catalogue.split(value.denominator()))
          if (catalogue.factors()[factor].kind == denominator_factors::Kind::mixed && !singular.count(factor))
            return false;
  return true;
}

namespace detail {
inline void write_record(const std::filesystem::path &file, const boost::json::object &record) {
  std::filesystem::create_directories(file.parent_path());
  const auto temporary = file.parent_path() / (file.filename().string() + "." + std::to_string(::getpid()) + ".tmp");
  {
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    out << boost::json::serialize(record);
    if (!out)
      throw std::runtime_error("finite-master record not written");
  }
  std::filesystem::rename(temporary, file);
}
inline boost::json::array integral_json(const ibp::Integral &a) {
  boost::json::array out;
  for (int power : a)
    out.push_back(power);
  return out;
}
inline ibp::Integral integral_from(const boost::json::value &value) {
  ibp::Integral out;
  for (const auto &power : value.as_array())
    out.push_back(static_cast<int>(power.to_number<std::int64_t>()));
  return out;
}
} // namespace detail

// The whole choice for one level. `reduce` returns the closure of the level
// for the given requests (the level's sources first) within the given number
// of seconds; all such reductions of a level together are allowed `seconds`
// (a level that needs longer keeps its masters). It is asked for the
// sources and, per sector with failing masters, the passing integrals of the
// lowest total powers that are at least as many as those masters; then for
// the next two passing integrals per master still open, round by round,
// while its sector has passing integrals within reach (as many further
// powers as the sector has lines beyond the first total that has one). The returned
// closure has rows for the sources only. What was decided is recorded in
// `directory` (unless empty) with the exchanged closure, and a later
// preparation of the same level reads it there.
inline level::Result choose(const ibp::PropagatorBasis &merged, const ExactField &field, std::size_t loops,
                            const Rational &dimension_at_zero, level::Result closure, const std::vector<ibp::Integral> &sources, std::size_t xi, std::size_t ei,
                            const std::function<level::Result(const std::vector<ibp::Integral> &, double)> &reduce, Report &report,
                            const std::filesystem::path &directory = {}, double seconds = 120) {
  namespace json = boost::json;
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
  // Only where the integrals the level is asked for are finite themselves.
  // Masters that are more finite than the requests move the poles into the
  // weights instead, and a level with both kinds of masters carries more
  // epsilon orders than before (planar double box: 21.0 -> 33.9 million
  // products with its three lowest levels partly exchanged).
  for (const auto &source : sources)
    if (!convergent(source)) {
      report.note = "kept: the integrals requested from this level are not finite at epsilon = 0";
      report.seconds = elapsed();
      return closure;
    }
  const auto trim = [&](level::Result level) {
    level.requested_integrals.resize(std::min(sources.size(), level.requested_integrals.size()));
    level.target_rows.resize(std::min(sources.size(), level.target_rows.size()));
    if (level.target_witnesses.size() > sources.size())
      level.target_witnesses.resize(sources.size());
    return level;
  };
  const auto finish = [&](level::Result level, std::vector<Step> steps) {
    report.divergent_after = failing(level);
    report.characters_after = moving_poles::characters(level.matrix);
    report.steps = std::move(steps);
    report.seconds = elapsed();
    return trim(std::move(level));
  };
  // The record of an earlier preparation of this level.
  std::filesystem::path file;
  if (!directory.empty()) {
    std::string text = "finite-masters-v7|" + std::to_string(loops) + "|" + dimension_at_zero.str() + "|" + std::to_string(xi) + "," + std::to_string(ei) + "," + std::to_string(seconds) + "|";
    for (const auto &name : closure.matrix.at(0).at(0).variables())
      text += name + ",";
    for (const std::vector<ibp::Integral> *list : std::initializer_list<const std::vector<ibp::Integral> *>{&closure.ordered_basis, &sources}) {
      text += "|";
      for (const auto &integral : *list) {
        for (int power : integral)
          text += std::to_string(power) + ",";
        text += ";";
      }
    }
    for (const auto &row : closure.matrix)
      for (const auto &value : row)
        text += value.is_zero() ? std::string(";") : value.str() + ";";
    file = directory / (artifacts::detail::sha256(text) + ".json");
    report.record = file;
    try {
      std::ifstream in(file);
      if (in) {
        const std::string bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        const auto record = json::parse(bytes);
        const auto &object = record.as_object();
        std::vector<Step> script;
        for (const auto &item : object.at("steps").as_array()) {
          const auto &step = item.as_object();
          script.push_back({static_cast<std::size_t>(step.at("replaced").to_number<std::uint64_t>()),
                            detail::integral_from(step.at("from")), detail::integral_from(step.at("to"))});
        }
        report.recorded = true;
        report.candidates = static_cast<std::size_t>(object.at("candidates").to_number<std::uint64_t>());
        report.note = std::string(object.at("note").as_string());
        if (script.empty())
          return finish(std::move(closure), {});
        // The exchanged closure itself, as it was written.
        const auto &payload = object.at("closure");
        if (artifacts::detail::sha256(artifacts::detail::canonical(payload)) != std::string(object.at("closure_sha256").as_string()))
          throw std::runtime_error("finite-master record is damaged");
        auto level = cached_level::detail::decode(payload, merged, field);
        level.differential_witnesses.clear();
        level.target_witnesses.clear();
        if (level.requested_integrals != sources || level.ordered_basis.size() != closure.ordered_basis.size())
          throw std::runtime_error("finite-master record does not belong to this level");
        return finish(std::move(level), std::move(script));
      }
    } catch (const std::exception &) {
      // An unreadable or stale record: the search decides.
    }
    report = Report{};
    report.record = file;
    report.masters = closure.ordered_basis.size();
    report.divergent_before = report.divergent_after = failing(closure);
    report.characters_before = report.characters_after = moving_poles::characters(closure.matrix);
  }
  // Per sector with a failing master: the next total power to enumerate and
  // the last one worth trying.
  struct Search {
    ibp::Integral pattern;
    long next, last;
    bool found = false;
    // Passing integrals of the totals enumerated so far that were not requested yet.
    std::vector<ibp::Integral> pending;
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
      searches.emplace(shape, Search{master, first, first + lines, false, {}});
    }
  const std::set<ibp::Integral> masters(closure.ordered_basis.begin(), closure.ordered_basis.end());
  std::vector<ibp::Integral> requests = sources;
  std::set<ibp::Integral> requested(sources.begin(), sources.end());
  level::Result best = closure;
  std::vector<Step> best_steps;
  // The requests beyond the sources of the reduction the best choice came from.
  std::vector<ibp::Integral> best_requests;
  // All masters or none (see above). The exchange is first made without
  // looking at intermediate bases, which may have singular points that the
  // finished one has not, and the finished basis is checked; if that fails,
  // with the check at every step, which can pass over a bad integral.
  bool complete = false;
  std::size_t open_masters = report.divergent_before;
  std::map<std::vector<bool>, std::size_t> wanted;
  const auto attempt = [&](const level::Result &from) {
    for (const bool guarded : {false, true}) {
      std::vector<Step> steps;
      auto now = exchange(from, xi, ei, convergent, steps, sources.size(), nullptr, guarded, false);
      const auto left = failing(now);
      if (left && !guarded) {
        // More candidates are needed whatever the check says.
        wanted.clear();
        for (const auto &master : now.ordered_basis)
          if (!convergent(master))
            ++wanted[sector(master)];
        open_masters = left;
        return;
      }
      if (!left && (guarded || regular(from.matrix, now.matrix, xi, ei))) {
        complete = true;
        open_masters = 0;
        best_steps = std::move(steps);
        best = std::move(now);
        best_requests.assign(requests.begin() + static_cast<std::ptrdiff_t>(sources.size()), requests.end());
        return;
      }
      if (guarded) {
        wanted.clear();
        for (const auto &master : now.ordered_basis)
          if (!convergent(master))
            ++wanted[sector(master)];
        open_masters = left;
      }
    }
  };
  // First with the integrals the closure already holds.
  attempt(closure);
  for (bool first = true; !complete; first = false) {
    std::vector<ibp::Integral> added;
    for (auto &[shape, search] : searches) {
      const auto need = wanted.find(shape);
      if (need == wanted.end())
        continue;
      // The first reduction takes whole totals, as many as give one
      // candidate per failing master. Not all of them are independent, and
      // the reductions get dearer with every dot: a later one takes only
      // the next two candidates per master that is still open.
      const std::size_t count = first ? need->second : 2 * need->second;
      std::size_t fresh = 0;
      while (fresh < count) {
        if (search.pending.empty()) {
          if (search.next > search.last)
            break;
          const auto total = search.next++;
          for (auto &integral : convergent_integrals(convergent, search.pattern, total)) {
            if (!search.found) {
              search.found = true;
              search.last = std::max(search.last, total + static_cast<long>(std::count(shape.begin(), shape.end(), true)));
            }
            if (!masters.count(integral) && !requested.count(integral))
              search.pending.push_back(std::move(integral));
          }
          continue;
        }
        const std::size_t take = first ? search.pending.size() : std::min(search.pending.size(), count - fresh);
        for (std::size_t i = 0; i < take; ++i)
          if (requested.insert(search.pending[i]).second) {
            added.push_back(search.pending[i]);
            ++fresh;
          }
        search.pending.erase(search.pending.begin(), search.pending.begin() + static_cast<std::ptrdiff_t>(take));
      }
    }
    if (added.empty())
      break;
    const double remaining = seconds - report.reduction_seconds;
    if (!(remaining >= 1)) {
      report.note = "kept: the time allowed for reducing candidate integrals is spent";
      break;
    }
    requests.insert(requests.end(), added.begin(), added.end());
    report.candidates += added.size();
    const auto before = elapsed();
    level::Result wider;
    try {
      wider = reduce(requests, remaining);
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
  if (!complete) {
    if (report.note.empty())
      report.note = "kept: " + std::to_string(open_masters) + " masters have no finite replacement without a moving pole";
    else if (report.note.rfind("kept: ", 0) != 0)
      report.note = "kept: " + report.note;
    best = closure;
    best_steps.clear();
    best_requests.clear();
  }
  if (!file.empty())
    try {
      json::array steps, extra;
      for (const auto &step : best_steps)
        steps.push_back(json::object{{"replaced", step.replaced}, {"from", detail::integral_json(step.from)}, {"to", detail::integral_json(step.to)}});
      for (const auto &integral : best_requests)
        extra.push_back(detail::integral_json(integral));
      json::object record{{"schema", "DiffExp.FiniteMasters/v1"}, {"steps", std::move(steps)},
                          {"requests", std::move(extra)}, {"candidates", report.candidates}, {"note", report.note}};
      if (complete) {
        // The closure is stored without witnesses: it is the certified one
        // of the recorded requests transformed by the recorded steps.
        auto stored = trim(best);
        stored.differential_witnesses.assign(stored.ordered_basis.size(), {});
        stored.target_witnesses.assign(stored.requested_integrals.size(), {});
        auto payload = cached_level::detail::payload(stored);
        record["closure_sha256"] = artifacts::detail::sha256(artifacts::detail::canonical(payload));
        record["closure"] = std::move(payload);
      }
      detail::write_record(file, record);
    } catch (const std::exception &) {
      // The record is a convenience; the result does not depend on it.
    }
  return finish(std::move(best), std::move(best_steps));
}

// Withdraw a recorded choice: later preparations keep the level's masters.
// A level below an exchanged one that cannot follow (its requests are finite
// now, its masters are not) is worse off than before, so the levels above it
// are taken back (four-loop banana with the first level exchanged alone:
// 73 -> 125 million products).
inline void revoke(const std::filesystem::path &record, const std::string &note) {
  if (record.empty())
    return;
  try {
    detail::write_record(record, boost::json::object{{"schema", "DiffExp.FiniteMasters/v1"}, {"steps", boost::json::array{}},
                                                     {"requests", boost::json::array{}}, {"candidates", 0}, {"note", note}});
  } catch (const std::exception &) {
  }
}
} // namespace diffexp::finite_masters
