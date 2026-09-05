#pragma once

#include "diffexp/kernel/physical_ode.hpp"
#include "diffexp/kernel/recurrence.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace diffexp::kernel {

// An SCC block may need coefficients whose recurrence-validity edge lies
// below its first stored physical coefficient.  Such a rectangle is not an
// AssembledResult: it is a transaction-local candidate which may be retained
// only after the complete parent q/C equation and its canonical seed prove
// it.  The ordinary assembler remains the sole producer of non-deferred
// results.
template <typename Scalar> struct SCCAssemblyCandidate {
  AssembledResult<Scalar> coefficients;
  bool requires_parent_certificate = false;
  std::int32_t recurrence_top_valid = kCompleteInfinity;
};

template <typename Scalar>
std::int32_t scc_assembly_honest_complete_max(
    const PreparedRecurrenceOperator<Scalar> &prepared,
    const RecurrenceProblem<Scalar> &problem,
    const RecurrenceResult<Scalar> &recurrence) {
  if (!prepared.assembly_matrix.has_value())
    throw RecurrenceError("E5", "compiled assembly matrix is missing");
  const auto frame_top = local_detail::checked_i32(
      static_cast<std::int64_t>(prepared.frame_base) + prepared.frame_width - 1,
      "SCC assembly frame top");
  const auto logs = static_cast<std::size_t>(problem.log_max) + 2;
  const auto expected =
      (static_cast<std::size_t>(problem.nmax) + 1) * logs * prepared.dimension;
  if (recurrence.validity.size() != expected)
    throw RecurrenceError("E5", "malformed recurrence validity tensor");
  const auto validity_index = [&](std::uint32_t n, std::uint32_t log,
                                  std::uint32_t component) {
    return ((static_cast<std::size_t>(n) * logs + log) * prepared.dimension +
            component);
  };
  std::int32_t complete_max = kCompleteInfinity;
  const auto &matrix = *prepared.assembly_matrix;
  if (!prepared.epsilon_regular_principal &&
      matrix.valuations.size() !=
          static_cast<std::size_t>(prepared.dimension) * prepared.dimension)
    throw RecurrenceError("E5", "malformed assembly valuation tensor");
  for (std::uint32_t n = 0; n <= problem.nmax; ++n)
    for (std::uint32_t log = 0; log <= problem.log_max; ++log)
      for (std::uint32_t row = 0; row < prepared.dimension; ++row) {
        std::int32_t row_valid =
            prepared.epsilon_regular_principal
                ? recurrence.validity[validity_index(n, log, row)]
                : kCompleteInfinity;
        if (!prepared.epsilon_regular_principal)
          for (std::uint32_t column = 0; column < prepared.dimension;
               ++column) {
            const auto valuation =
                matrix.valuations[static_cast<std::size_t>(row) *
                                      prepared.dimension +
                                  column];
            if (valuation != kCompleteInfinity)
              row_valid = std::min(
                  row_valid,
                  detail::valid_shift(
                      recurrence.validity[validity_index(n, log, column)],
                      valuation, frame_top));
          }
        if (row_valid != kCompleteInfinity)
          complete_max = std::min(complete_max, row_valid);
      }
  return complete_max == kCompleteInfinity ? frame_top : complete_max;
}

template <typename Scalar>
SCCAssemblyCandidate<Scalar> assemble_scc_recurrence_candidate(
    const PreparedRecurrenceOperator<Scalar> &prepared,
    const RecurrenceProblem<Scalar> &problem,
    RecurrenceResult<Scalar> &recurrence) {
  const auto honest_complete_max =
      scc_assembly_honest_complete_max(prepared, problem, recurrence);
  try {
    auto assembled = assemble_recurrence(prepared, problem, recurrence);
    const bool deferred = assembled.min_power > honest_complete_max;
    return {std::move(assembled), deferred, recurrence.top_valid};
  } catch (const RecurrenceError &error) {
    if (error.id != "E6" ||
        std::string(error.what()) !=
            "compiled assembly exhausted completeness below stored content")
      throw;
  }

  // Reuse the canonical coefficient transformation without copying the
  // potentially large recurrence tensor.  Suppressing validity is confined
  // to this call and the original evidence is restored on every exit path.
  std::vector<std::int32_t> deferred_validity(recurrence.validity.size(),
                                              kCompleteInfinity);
  recurrence.validity.swap(deferred_validity);
  try {
    auto candidate = assemble_recurrence(prepared, problem, recurrence);
    recurrence.validity.swap(deferred_validity);
    return {std::move(candidate), true, recurrence.top_valid};
  } catch (...) {
    recurrence.validity.swap(deferred_validity);
    throw;
  }
}

template <typename Scalar>
void require_exact_domain_for_deferred_scc_candidate(
    const SCCAssemblyCandidate<Scalar> &candidate) {
  if (!candidate.requires_parent_certificate)
    return;
  if constexpr (std::is_same_v<Scalar, ComplexBall>) {
    throw RecurrenceError("E5", "native Acb deferred SCC completeness requires "
                                "the exact Rational shadow");
  } else if constexpr (!std::is_same_v<Scalar, Rational>) {
    throw RecurrenceError(
        "E5",
        "deferred SCC completeness supports only an exact Rational parent");
  }
}

struct SCCFormalResidualCertificate {
  EpsilonWindow epsilon;
  std::uint32_t taylor_complete_max = 0;
  std::size_t exact_tag_count = 0;
  std::size_t coefficient_rows = 0;
};

namespace scc_completeness_detail {

struct FormalTag {
  std::string a;
  std::string b;
  std::uint32_t log_power = 0;

  friend bool operator<(const FormalTag &left, const FormalTag &right) {
    return std::tie(left.a, left.b, left.log_power) <
           std::tie(right.a, right.b, right.log_power);
  }
};

struct ExactVector {
  EpsilonWindow epsilon;
  std::uint32_t dimension = 0;
  std::vector<Rational> coefficients;

  Rational &at(std::int32_t power, std::uint32_t component) {
    return coefficients[static_cast<std::size_t>(power - epsilon.min_power) *
                            dimension +
                        component];
  }
  const Rational &at(std::int32_t power, std::uint32_t component) const {
    return coefficients[static_cast<std::size_t>(power - epsilon.min_power) *
                            dimension +
                        component];
  }
};

using FormalTaylorSlab = std::vector<ExactVector>;
using FormalSlabs = std::map<FormalTag, FormalTaylorSlab>;

struct CanonicalAClass {
  std::string b;
  Rational base;
};

inline bool material_sector(const LocalSector<Rational> &sector) {
  return std::any_of(
      sector.coefficients.begin(), sector.coefficients.end(),
      [](const Rational &coefficient) { return !coefficient.is_zero(); });
}

inline bool exact_integer(const Rational &value) {
  return value.str().find('/') == std::string::npos;
}

inline std::uint32_t exact_nonnegative_integer_value(const Rational &value) {
  if (value.sign() < 0 || !exact_integer(value))
    throw std::logic_error(
        "formal tower alignment produced a noninteger negative shift");
  const auto text = value.str();
  const auto parsed = std::stoull(text);
  if (parsed > std::numeric_limits<std::uint32_t>::max())
    throw std::overflow_error(
        "formal tower alignment shift exceeds uint32 range");
  return static_cast<std::uint32_t>(parsed);
}

inline std::vector<CanonicalAClass>
canonical_a_classes(const LocalSolution<Rational> &solution) {
  std::vector<CanonicalAClass> classes;
  for (const auto &sector : solution.sectors) {
    // Exact-zero sectors are storage artifacts, not formal solutions.  In
    // particular, a zero high-log or CASE-P sector must not enlarge a tower
    // or create an all-future recurrence proof obligation.
    if (!material_sector(sector))
      continue;
    if (sector.a.domain != ExactDomain::Rational ||
        sector.b.domain != ExactDomain::Rational)
      throw std::invalid_argument(
          "SCC parent completeness requires exact Rational Frobenius tags");
    const Rational a(sector.a.canonical);
    auto found = std::find_if(classes.begin(), classes.end(),
                              [&](const auto &candidate) {
                                return candidate.b == sector.b.canonical &&
                                       exact_integer(a - candidate.base);
                              });
    if (found == classes.end()) {
      classes.push_back(CanonicalAClass{sector.b.canonical, a});
    } else if (a < found->base) {
      found->base = a;
    }
  }
  return classes;
}

inline const CanonicalAClass &
canonical_a_class(const std::vector<CanonicalAClass> &classes,
                  const LocalSector<Rational> &sector) {
  const Rational a(sector.a.canonical);
  const auto found =
      std::find_if(classes.begin(), classes.end(), [&](const auto &candidate) {
        return candidate.b == sector.b.canonical &&
               exact_integer(a - candidate.base);
      });
  if (found == classes.end())
    throw std::logic_error(
        "formal tower lost its canonical integer-shift class");
  return *found;
}

inline ExactVector zero_vector(EpsilonWindow epsilon, std::uint32_t dimension) {
  ExactVector result;
  result.epsilon = epsilon;
  result.dimension = dimension;
  result.coefficients.assign(epsilon.width() * dimension, Rational(0));
  return result;
}

inline FormalTaylorSlab &slab_for(FormalSlabs &slabs, const FormalTag &tag,
                                  EpsilonWindow epsilon,
                                  std::uint32_t dimension,
                                  std::uint32_t taylor_complete_max) {
  auto [found, inserted] = slabs.try_emplace(tag);
  if (inserted) {
    found->second.reserve(static_cast<std::size_t>(taylor_complete_max) + 1);
    for (std::uint32_t n = 0; n <= taylor_complete_max; ++n)
      found->second.push_back(zero_vector(epsilon, dimension));
  }
  return found->second;
}

inline FormalSlabs
collect_formal_slabs(const LocalSolution<Rational> &solution) {
  local_detail::validate_local_solution(solution, true);
  const auto classes = canonical_a_classes(solution);
  FormalSlabs result;
  for (const auto &sector : solution.sectors) {
    if (!material_sector(sector))
      continue;
    const auto &canonical = canonical_a_class(classes, sector);
    const auto shift = exact_nonnegative_integer_value(
        Rational(sector.a.canonical) - canonical.base);
    const FormalTag tag{canonical.base.str(), sector.b.canonical,
                        sector.log_power};
    auto &slab = slab_for(result, tag, solution.epsilon, solution.dimension,
                          solution.taylor_complete_max);
    for (std::int32_t power = solution.epsilon.min_power;
         power <= solution.epsilon.complete_max; ++power) {
      const auto epsilon_index =
          static_cast<std::size_t>(power - solution.epsilon.min_power);
      for (std::uint32_t n = 0; n <= solution.taylor_complete_max; ++n) {
        const auto canonical_n = static_cast<std::uint64_t>(n) + shift;
        if (canonical_n > solution.taylor_complete_max)
          break;
        for (std::uint32_t component = 0; component < solution.dimension;
             ++component)
          slab[static_cast<std::size_t>(canonical_n)].at(power, component) +=
              sector.coefficients[local_detail::sector_index(
                  solution, epsilon_index, n, component)];
      }
    }
  }
  return result;
}

inline FormalSlabs
collect_theta_slabs(const LocalSolution<Rational> &solution) {
  local_detail::validate_local_solution(solution, true);
  const auto classes = canonical_a_classes(solution);
  FormalSlabs result;
  for (const auto &sector : solution.sectors) {
    if (!material_sector(sector))
      continue;
    const auto &canonical = canonical_a_class(classes, sector);
    const Rational a(sector.a.canonical);
    const Rational b(sector.b.canonical);
    const auto shift = exact_nonnegative_integer_value(a - canonical.base);
    const FormalTag same_tag{canonical.base.str(), sector.b.canonical,
                             sector.log_power};
    auto &same = slab_for(result, same_tag, solution.epsilon,
                          solution.dimension, solution.taylor_complete_max);
    FormalTaylorSlab *lower = nullptr;
    if (sector.log_power > 0) {
      const FormalTag lower_tag{canonical.base.str(), sector.b.canonical,
                                sector.log_power - 1};
      lower = &slab_for(result, lower_tag, solution.epsilon, solution.dimension,
                        solution.taylor_complete_max);
    }
    for (std::int32_t power = solution.epsilon.min_power;
         power <= solution.epsilon.complete_max; ++power) {
      const auto epsilon_index =
          static_cast<std::size_t>(power - solution.epsilon.min_power);
      for (std::uint32_t n = 0; n <= solution.taylor_complete_max; ++n) {
        const auto canonical_n = static_cast<std::uint64_t>(n) + shift;
        if (canonical_n > solution.taylor_complete_max)
          break;
        const Rational a_plus_n = a + Rational(std::to_string(n));
        for (std::uint32_t component = 0; component < solution.dimension;
             ++component) {
          const auto value = sector.coefficients[local_detail::sector_index(
              solution, epsilon_index, n, component)];
          same[static_cast<std::size_t>(canonical_n)].at(power, component) +=
              a_plus_n * value;
          if (epsilon_index > 0) {
            const auto previous =
                sector.coefficients[local_detail::sector_index(
                    solution, epsilon_index - 1, n, component)];
            same[static_cast<std::size_t>(canonical_n)].at(power, component) +=
                b * previous;
            if (lower != nullptr)
              (*lower)[static_cast<std::size_t>(canonical_n)].at(
                  power, component) += previous;
          }
        }
      }
    }
  }
  return result;
}

inline std::string bounded_rational_witness(const Rational &value) {
  auto witness = value.str();
  constexpr std::size_t kWitnessLimit = 512;
  const auto witness_size = witness.size();
  if (witness_size > kWitnessLimit)
    witness = witness.substr(0, kWitnessLimit) +
              "...[decimal-bytes=" + std::to_string(witness_size) + "]";
  return witness;
}

inline void
accumulate_multiplier(const ExactEpsilonRational<Rational> &multiplier,
                      const ExactVector &source, std::uint32_t source_component,
                      ExactVector &target, std::uint32_t target_component,
                      const Rational &sign, const FormalTag &tag,
                      std::uint32_t taylor_index) {
  physical_ode_detail::validate_rational(multiplier,
                                         "SCC parent exact epsilon multiplier");
  if (multiplier.zero || source_component >= source.dimension ||
      target_component >= target.dimension)
    throw std::invalid_argument(
        "SCC parent exact formal multiplier has an invalid component");
  const auto product_complete = physical_ode_detail::shifted_power(
      source.epsilon.complete_max, multiplier.valuation,
      "SCC parent exact product complete maximum");
  if (target.epsilon.complete_max > product_complete)
    throw std::domain_error("SCC parent formal residual needs coefficients "
                            "above the retained epsilon reservoir for tag (a=" +
                            tag.a + ",b=" + tag.b +
                            ",log=" + std::to_string(tag.log_power) +
                            "), taylor=" + std::to_string(taylor_index));
  const auto width = source.epsilon.width();
  std::vector<Rational> quotient(width, Rational(0));
  for (std::size_t offset = 0; offset < width; ++offset) {
    Rational rhs(0);
    for (std::size_t degree = 0;
         degree < multiplier.numerator.size() && degree <= offset; ++degree) {
      const auto source_power = static_cast<std::int32_t>(
          static_cast<std::int64_t>(source.epsilon.min_power) +
          static_cast<std::int64_t>(offset - degree));
      rhs += multiplier.numerator[degree] *
             source.at(source_power, source_component);
    }
    for (std::size_t degree = 1;
         degree < multiplier.denominator.size() && degree <= offset; ++degree)
      rhs -= multiplier.denominator[degree] * quotient[offset - degree];
    quotient[offset] = rhs / multiplier.denominator.front();
    const auto target_power64 =
        static_cast<std::int64_t>(source.epsilon.min_power) +
        static_cast<std::int64_t>(offset) + multiplier.valuation;
    if (target_power64 < target.epsilon.min_power ||
        target_power64 > target.epsilon.complete_max)
      continue;
    target.at(static_cast<std::int32_t>(target_power64), target_component) +=
        sign * quotient[offset];
  }
}

} // namespace scc_completeness_detail

// Prove the stored Rational parent slab coefficientwise in the formal chart
// variable.  This is exact equality, not a ball-consistency test: every
// coefficient of q theta(F)-C F must be Rational zero.  Acb callers must use
// the paired Rational-shadow solve and specialize that exact result.
inline SCCFormalResidualCertificate certify_scc_parent_exact_formal_residual(
    const PreparedPhysicalClearedODE<Rational> &equation,
    const LocalSolution<Rational> &parent, EpsilonWindow claimed_epsilon,
    std::uint32_t claimed_taylor_complete_max,
    bool admit_maximal_exact_prefix = false) {
  physical_ode_detail::validate_ode(equation);
  local_detail::validate_local_solution(parent, true);
  if (parent.dimension != equation.dimension ||
      claimed_epsilon.min_power < parent.epsilon.min_power ||
      claimed_epsilon.complete_max > parent.epsilon.complete_max ||
      claimed_taylor_complete_max > parent.taylor_complete_max)
    throw std::invalid_argument("SCC parent formal residual claim lies outside "
                                "its retained local slab");

  const auto values = scc_completeness_detail::collect_formal_slabs(parent);
  const auto theta_values =
      scc_completeness_detail::collect_theta_slabs(parent);
  scc_completeness_detail::FormalSlabs residual;
  const Rational plus(1);
  const Rational minus(-1);

  for (const auto &[tag, slab] : theta_values) {
    auto &target = scc_completeness_detail::slab_for(
        residual, tag, claimed_epsilon, parent.dimension,
        claimed_taylor_complete_max);
    for (std::size_t lag = 0; lag < equation.q_lags.size(); ++lag) {
      const auto &q = equation.q_lags[lag];
      if (q.zero || lag > claimed_taylor_complete_max)
        continue;
      for (std::uint32_t source_n = 0;
           source_n + lag <= claimed_taylor_complete_max; ++source_n) {
        const auto output_n = static_cast<std::uint32_t>(source_n + lag);
        for (std::uint32_t component = 0; component < parent.dimension;
             ++component)
          scc_completeness_detail::accumulate_multiplier(
              q, slab[source_n], component, target[output_n], component, plus,
              tag, output_n);
      }
    }
  }

  for (const auto &[tag, slab] : values) {
    auto &target = scc_completeness_detail::slab_for(
        residual, tag, claimed_epsilon, parent.dimension,
        claimed_taylor_complete_max);
    for (std::size_t lag = 0; lag < equation.c_lags.size(); ++lag) {
      const auto &entries = equation.c_lags[lag];
      if (entries.empty() || lag > claimed_taylor_complete_max)
        continue;
      for (std::uint32_t source_n = 0;
           source_n + lag <= claimed_taylor_complete_max; ++source_n) {
        const auto output_n = static_cast<std::uint32_t>(source_n + lag);
        for (const auto &entry : entries)
          scc_completeness_detail::accumulate_multiplier(
              entry.value, slab[source_n], entry.column, target[output_n],
              entry.row, minus, tag, output_n);
      }
    }
  }

  if (residual.empty())
    throw std::domain_error(
        "SCC parent formal residual has no exact formal sectors");
  std::size_t rows = 0;
  std::optional<std::int32_t> first_nonzero_power;
  for (const auto &[tag, slab] : residual)
    for (std::uint32_t n = 0; n <= claimed_taylor_complete_max; ++n)
      for (std::int32_t power = claimed_epsilon.min_power;
           power <= claimed_epsilon.complete_max; ++power)
        for (std::uint32_t component = 0; component < parent.dimension;
             ++component) {
          ++rows;
          if (!slab[n].at(power, component).is_zero()) {
            if (admit_maximal_exact_prefix) {
              if (!first_nonzero_power.has_value() ||
                  power < *first_nonzero_power)
                first_nonzero_power = power;
              continue;
            }
            auto lhs = scc_completeness_detail::zero_vector(claimed_epsilon,
                                                            parent.dimension);
            auto rhs = scc_completeness_detail::zero_vector(claimed_epsilon,
                                                            parent.dimension);
            const auto theta_found = theta_values.find(tag);
            const auto value_found = values.find(tag);
            if (theta_found == theta_values.end() ||
                value_found == values.end())
              throw std::logic_error(
                  "SCC parent residual witness lost its exact formal tag");
            for (std::size_t lag = 0; lag < equation.q_lags.size() && lag <= n;
                 ++lag) {
              const auto &q = equation.q_lags[lag];
              if (q.zero)
                continue;
              const auto source_n = n - static_cast<std::uint32_t>(lag);
              for (std::uint32_t source_component = 0;
                   source_component < parent.dimension; ++source_component)
                scc_completeness_detail::accumulate_multiplier(
                    q, theta_found->second[source_n], source_component, lhs,
                    source_component, plus, tag, n);
            }
            for (std::size_t lag = 0; lag < equation.c_lags.size() && lag <= n;
                 ++lag) {
              const auto &entries = equation.c_lags[lag];
              if (entries.empty())
                continue;
              const auto source_n = n - static_cast<std::uint32_t>(lag);
              for (const auto &entry : entries)
                scc_completeness_detail::accumulate_multiplier(
                    entry.value, value_found->second[source_n], entry.column,
                    rhs, entry.row, plus, tag, n);
            }
            throw std::domain_error(
                "SCC parent physical q/C residual is not exact zero for tag "
                "(a=" +
                tag.a + ",b=" + tag.b +
                ",log=" + std::to_string(tag.log_power) + "), taylor=" +
                std::to_string(n) + ", epsilon=" + std::to_string(power) +
                ", component=" + std::to_string(component) + ", coefficient=" +
                scc_completeness_detail::bounded_rational_witness(
                    slab[n].at(power, component)) +
                ", q_theta=" +
                scc_completeness_detail::bounded_rational_witness(
                    lhs.at(power, component)) +
                ", c_value=" +
                scc_completeness_detail::bounded_rational_witness(
                    rhs.at(power, component)));
          }
        }
  if (first_nonzero_power.has_value()) {
    if (*first_nonzero_power <= claimed_epsilon.min_power) {
      // Re-run the one-row strict claim to retain the detailed exact
      // q-theta/C witness instead of returning an empty prefix.
      return certify_scc_parent_exact_formal_residual(
          equation, parent,
          {claimed_epsilon.min_power, claimed_epsilon.min_power},
          claimed_taylor_complete_max, false);
    }
    claimed_epsilon.complete_max =
        static_cast<std::int32_t>(*first_nonzero_power - 1);
    rows = residual.size() *
           (static_cast<std::size_t>(claimed_taylor_complete_max) + 1) *
           claimed_epsilon.width() * parent.dimension;
  }
  return {claimed_epsilon, claimed_taylor_complete_max, residual.size(), rows};
}

} // namespace diffexp::kernel
