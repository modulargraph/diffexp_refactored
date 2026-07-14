#pragma once

#include "diffexp2/recurrence.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace diffexp2 {

struct ExactAffineIndicialRoot {
  Rational a{0};
  Rational b{0};
};

inline bool operator==(const ExactAffineIndicialRoot& left,
                       const ExactAffineIndicialRoot& right) {
  return left.a == right.a && left.b == right.b;
}

/* A proof object for the exact prepared spectral frame.  `columns` keeps the
   prepared Jordan-chain order (eigenvector first), so the certificate binds
   both the affine root and the nilpotent-chain convention. */
struct ExactJordanBlockCertificate {
  std::uint32_t block_index = 0;
  std::vector<std::uint32_t> columns;
  ExactAffineIndicialRoot root;

  [[nodiscard]] std::uint32_t size() const {
    return static_cast<std::uint32_t>(columns.size());
  }
};

struct ExactJordanIndicialCertificate {
  std::uint32_t dimension = 0;
  std::vector<ExactJordanBlockCertificate> blocks;
  // One entry per prepared spectral column.  This makes later schedule and
  // seed validation independent of any assumed singleton/contiguous layout.
  std::vector<std::uint32_t> block_of_column;
  std::vector<std::uint32_t> position_in_block;
};

struct ExactJordanScheduleStepCertificate {
  std::uint32_t taylor_index = 0;
  std::uint32_t block_index = 0;
  StepCase kind = StepCase::Taylor;
  Rational d_a{0};
  Rational d_b{0};
  std::uint32_t jordan_size = 0;

  // For R, q is the exact Jordan-chain width coupled by the log ladder and
  // q-1 is its intrinsic homogeneous log degree.  For P, the inverse of
  // (d_b eps I-N_q) can contain eps^(-q), so q is the exact worst pole depth
  // at this Taylor layer.  These are structural facts, not numeric guesses.
  std::uint32_t resonant_jordan_chain_length = 0;
  std::uint32_t intrinsic_homogeneous_log_degree = 0;
  std::uint32_t pseudo_epsilon_pole_depth = 0;
};

struct ExactJordanScheduleCertificate {
  ExactAffineIndicialRoot target;
  std::vector<std::vector<ExactJordanScheduleStepCertificate>> steps;
  bool contains_true_resonance = false;
  bool contains_pseudo = false;
  std::uint32_t max_resonant_jordan_chain_length = 0;
  std::uint32_t max_pseudo_epsilon_pole_depth = 0;
};

namespace singular_indicial_detail {

using LaurentPolynomial = std::map<std::int32_t, Rational>;

[[noreturn]] inline void fail(const std::string& detail) {
  throw RecurrenceError("E5", detail);
}

inline void add_coefficient(LaurentPolynomial& polynomial,
                            std::int32_t power, const Rational& value) {
  if (value.is_zero()) return;
  auto [entry, inserted] = polynomial.try_emplace(power, Rational(0));
  entry->second += value;
  if (entry->second.is_zero()) polynomial.erase(entry);
}

inline std::int32_t checked_power_sum(std::int32_t left,
                                      std::int32_t right) {
  const auto result = static_cast<std::int64_t>(left) + right;
  if (result < std::numeric_limits<std::int32_t>::min() ||
      result > std::numeric_limits<std::int32_t>::max())
    fail("exact indicial Laurent power overflows int32");
  return static_cast<std::int32_t>(result);
}

inline LaurentPolynomial add(const LaurentPolynomial& left,
                             const LaurentPolynomial& right) {
  auto result = left;
  for (const auto& [power, value] : right)
    add_coefficient(result, power, value);
  return result;
}

inline LaurentPolynomial scale(const LaurentPolynomial& input,
                               const Rational& factor) {
  LaurentPolynomial result;
  if (factor.is_zero()) return result;
  for (const auto& [power, value] : input)
    add_coefficient(result, power, value * factor);
  return result;
}

inline LaurentPolynomial shift(const LaurentPolynomial& input,
                               std::int32_t amount) {
  LaurentPolynomial result;
  for (const auto& [power, value] : input)
    add_coefficient(result, checked_power_sum(power, amount), value);
  return result;
}

inline LaurentPolynomial multiply(const LaurentPolynomial& left,
                                  const LaurentPolynomial& right) {
  LaurentPolynomial result;
  for (const auto& [left_power, left_value] : left)
    for (const auto& [right_power, right_value] : right)
      add_coefficient(result, checked_power_sum(left_power, right_power),
                      left_value * right_value);
  return result;
}

inline Rational coefficient(const LaurentPolynomial& polynomial,
                            std::int32_t power) {
  const auto found = polynomial.find(power);
  return found == polynomial.end() ? Rational(0) : found->second;
}

struct LaurentFraction {
  LaurentPolynomial numerator;
  LaurentPolynomial denominator{{0, Rational(1)}};
};

inline void add_fraction(LaurentFraction& accumulator,
                         const LaurentPolynomial& numerator,
                         const LaurentPolynomial& denominator) {
  if (numerator.empty()) return;
  if (denominator.empty()) fail("exact indicial denominator is zero");
  if (accumulator.numerator.empty()) {
    accumulator.numerator = numerator;
    accumulator.denominator = denominator;
    return;
  }
  accumulator.numerator = add(
      multiply(accumulator.numerator, denominator),
      multiply(numerator, accumulator.denominator));
  accumulator.denominator =
      multiply(accumulator.denominator, denominator);
  if (accumulator.numerator.empty())
    accumulator.denominator = {{0, Rational(1)}};
}

inline LaurentPolynomial d0_polynomial(
    const PreparedRecurrenceOperator<Rational>& prepared) {
  if (prepared.d_lags.empty()) fail("exact indicial certificate has no d0 lag");
  LaurentPolynomial d0;
  for (const auto& term : prepared.d_lags.front())
    add_coefficient(d0, term.shift, term.value);
  if (d0.empty()) fail("exact indicial certificate has identically zero d0");

  // A retained scalar inverse is used directly by RecurrenceSolver.  Prove
  // it is the inverse of the same exact d0 instead of trusting a cache flag.
  if (prepared.d0_inverse_scalar.has_value()) {
    if (d0.size() != 1 || d0.begin()->first != 0 ||
        !(d0.begin()->second * *prepared.d0_inverse_scalar == Rational(1)))
      fail("retained scalar d0 inverse contradicts the exact prepared d0");
  }
  return d0;
}

inline void collect_matrix_shifts(
    const std::vector<MatrixShift<Rational>>& shifts,
    std::uint32_t dimension, std::vector<LaurentPolynomial>& entries,
    const char* label) {
  for (const auto& matrix : shifts) {
    for (const auto& entry : matrix.entries) {
      if (entry.row >= dimension || entry.col >= dimension)
        fail(std::string(label) + " entry is outside the prepared dimension");
      const auto index = static_cast<std::size_t>(entry.row) * dimension +
                         entry.col;
      add_coefficient(entries[index], matrix.shift, entry.value);
    }
  }
}

inline LaurentPolynomial denominator_polynomial(
    const std::vector<Rational>& coefficients) {
  if (coefficients.empty() || coefficients.front().is_zero())
    fail("exact indicial rational denominator has zero constant term");
  LaurentPolynomial denominator;
  for (std::size_t power = 0; power < coefficients.size(); ++power) {
    if (power > static_cast<std::size_t>(
                    std::numeric_limits<std::int32_t>::max()))
      fail("exact indicial rational denominator degree overflows int32");
    add_coefficient(denominator, static_cast<std::int32_t>(power),
                    coefficients[power]);
  }
  return denominator;
}

inline std::vector<LaurentFraction> nhat_zero_entries(
    const PreparedRecurrenceOperator<Rational>& prepared) {
  if (prepared.nhat_lags.empty())
    fail("exact indicial certificate has no Nhat_0 lag");
  if (prepared.epsilon_regular_principal &&
      !prepared.spectral_principal_lag.has_value())
    fail("epsilon-regular exact indicial certificate has no spectral Nhat_0 lag");
  const auto dimension = prepared.dimension;
  if (dimension == 0) fail("exact indicial certificate has zero dimension");
  if (dimension > std::numeric_limits<std::size_t>::max() / dimension)
    fail("exact indicial matrix size overflows size_t");
  const auto count = static_cast<std::size_t>(dimension) * dimension;
  const auto& lag = prepared.epsilon_regular_principal
      ? *prepared.spectral_principal_lag
      : prepared.nhat_lags.front();

  std::vector<LaurentPolynomial> polynomial_entries(count);
  collect_matrix_shifts(lag.polynomial, dimension, polynomial_entries,
                        "Nhat_0 polynomial");
  std::vector<LaurentFraction> result(count);
  for (std::size_t index = 0; index < count; ++index)
    result[index].numerator = std::move(polynomial_entries[index]);

  for (const auto& group : lag.rational) {
    if (group.denominator_index >= prepared.rational_denominators.size())
      fail("exact indicial rational denominator index is out of range");
    const auto denominator = denominator_polynomial(
        prepared.rational_denominators[group.denominator_index]);
    std::vector<LaurentPolynomial> numerator_entries(count);
    collect_matrix_shifts(group.numerator, dimension, numerator_entries,
                          "Nhat_0 rational numerator");
    for (std::size_t index = 0; index < count; ++index)
      add_fraction(result[index], numerator_entries[index], denominator);
  }
  return result;
}

inline ExactAffineIndicialRoot affine_quotient(
    const LaurentFraction& value, const LaurentPolynomial& d0,
    const std::string& label) {
  const auto base = multiply(value.denominator, d0);
  if (base.empty()) fail(label + " has a zero exact divisor");
  const auto minimum_power = base.begin()->first;
  const auto leading = base.begin()->second;
  const auto a = coefficient(value.numerator, minimum_power) / leading;
  const auto after_a = add(value.numerator, scale(base, -a));
  const auto next_power = checked_power_sum(minimum_power, 1);
  const auto b = coefficient(after_a, next_power) / leading;
  if (after_a != scale(shift(base, 1), b))
    fail(label + " is not exactly affine in epsilon after division by d0");
  return {a, b};
}

inline bool quotient_is_constant(const LaurentFraction& value,
                                 const LaurentPolynomial& d0,
                                 const Rational& expected) {
  return value.numerator ==
         scale(multiply(value.denominator, d0), expected);
}

inline StepCase classify_step(const Rational& d_a, const Rational& d_b) {
  if (!d_a.is_zero()) return StepCase::Taylor;
  if (!d_b.is_zero()) return StepCase::Pseudo;
  return StepCase::Resonant;
}

inline const char* step_name(StepCase kind) {
  switch (kind) {
    case StepCase::Taylor: return "T";
    case StepCase::Pseudo: return "P";
    case StepCase::Resonant: return "R";
  }
  return "?";
}

}  // namespace singular_indicial_detail

/* Derive and prove the complete exact Jordan form represented by a prepared
   Rational recurrence operator.  The proof checks

       Nhat_0(eps) / d0(eps) = direct_sum_i(
           (a_i + b_i eps) I + N_{q_i})

   with +1 on each declared superdiagonal and exact zero everywhere else.
   Rational groups are compared by exact cross multiplication; no sampled or
   midpoint coefficient can influence the result. */
inline ExactJordanIndicialCertificate certify_exact_affine_jordan_operator(
    const PreparedRecurrenceOperator<Rational>& prepared) {
  using namespace singular_indicial_detail;
  if (prepared.dimension == 0)
    fail("exact Jordan indicial certificate has zero dimension");
  if (prepared.blocks.empty())
    fail("exact Jordan indicial certificate has no declared blocks");

  ExactJordanIndicialCertificate certificate;
  certificate.dimension = prepared.dimension;
  certificate.block_of_column.assign(
      prepared.dimension, std::numeric_limits<std::uint32_t>::max());
  certificate.position_in_block.assign(
      prepared.dimension, std::numeric_limits<std::uint32_t>::max());
  certificate.blocks.reserve(prepared.blocks.size());

  for (std::size_t block_index = 0; block_index < prepared.blocks.size();
       ++block_index) {
    const auto& declared = prepared.blocks[block_index];
    if (declared.columns.empty())
      fail("exact Jordan indicial certificate contains an empty block");
    if (declared.columns.size() >
        static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
      fail("exact Jordan block size overflows uint32");
    if (block_index > std::numeric_limits<std::uint32_t>::max())
      fail("exact Jordan block index overflows uint32");
    for (std::size_t position = 0; position < declared.columns.size();
         ++position) {
      const auto column = declared.columns[position];
      if (column >= prepared.dimension)
        fail("exact Jordan block column is outside the prepared dimension");
      if (certificate.block_of_column[column] !=
          std::numeric_limits<std::uint32_t>::max())
        fail("exact Jordan blocks contain a duplicate column");
      if (position > std::numeric_limits<std::uint32_t>::max())
        fail("exact Jordan block position overflows uint32");
      certificate.block_of_column[column] =
          static_cast<std::uint32_t>(block_index);
      certificate.position_in_block[column] =
          static_cast<std::uint32_t>(position);
    }
  }
  if (std::any_of(certificate.block_of_column.begin(),
                  certificate.block_of_column.end(),
                  [](std::uint32_t value) {
                    return value == std::numeric_limits<std::uint32_t>::max();
                  }))
    fail("exact Jordan blocks do not partition the prepared dimension");

  const auto d0 = d0_polynomial(prepared);
  const auto entries = nhat_zero_entries(prepared);
  for (std::size_t block_index = 0; block_index < prepared.blocks.size();
       ++block_index) {
    const auto& declared = prepared.blocks[block_index];
    const auto first_column = declared.columns.front();
    const auto first_index = static_cast<std::size_t>(first_column) *
                             prepared.dimension + first_column;
    const auto root = affine_quotient(
        entries[first_index], d0,
        "exact Jordan block " + std::to_string(block_index) + " diagonal");
    for (const auto column : declared.columns) {
      const auto index = static_cast<std::size_t>(column) *
                         prepared.dimension + column;
      const auto candidate = affine_quotient(
          entries[index], d0,
          "exact Jordan block " + std::to_string(block_index) +
              " diagonal");
      if (!(candidate == root))
        fail("exact Jordan block diagonal entries have unequal affine roots");
    }
    certificate.blocks.push_back(
        {static_cast<std::uint32_t>(block_index), declared.columns, root});
  }

  for (std::uint32_t row = 0; row < prepared.dimension; ++row) {
    for (std::uint32_t column = 0; column < prepared.dimension; ++column) {
      if (row == column) continue;  // already certified as the block root
      const auto index = static_cast<std::size_t>(row) * prepared.dimension +
                         column;
      const auto same_block = certificate.block_of_column[row] ==
                              certificate.block_of_column[column];
      const auto unit_superdiagonal = same_block &&
          static_cast<std::uint64_t>(certificate.position_in_block[row]) + 1 ==
              certificate.position_in_block[column];
      if (unit_superdiagonal) {
        if (!quotient_is_constant(entries[index], d0, Rational(1)))
          fail("exact prepared Jordan superdiagonal is not the unit +1");
      } else if (!entries[index].numerator.empty()) {
        fail("Nhat_0/d0 is nonzero outside the declared block-Jordan structure");
      }
    }
  }
  return certificate;
}

/* Validate a submitted finite T/P/R schedule against the retained proof.
   Every step is reconstructed from exact roots as

       d_a = a_target + n - a_i,   d_b = b_target - b_i.

   P is returned as an exact classification (with its q-order Laurent risk),
   even if a caller's current execution capability elects to reject it. */
inline ExactJordanScheduleCertificate certify_exact_affine_jordan_schedule(
    const ExactJordanIndicialCertificate& indicial,
    const Rational& a_target, const Rational& b_target,
    const std::vector<std::vector<BlockStep<Rational>>>& schedule) {
  using namespace singular_indicial_detail;
  if (indicial.dimension == 0 || indicial.blocks.empty())
    fail("cannot validate a schedule against an empty indicial certificate");
  if (schedule.empty()) fail("exact Jordan recurrence schedule is empty");
  if (schedule.size() - 1 >
      static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
    fail("exact Jordan recurrence Taylor index overflows uint32");
  if (indicial.blocks.size() >
      static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
    fail("exact Jordan recurrence block index overflows uint32");

  ExactJordanScheduleCertificate certificate;
  certificate.target = {a_target, b_target};
  certificate.steps.resize(schedule.size());
  for (std::size_t n = 0; n < schedule.size(); ++n) {
    if (schedule[n].size() != indicial.blocks.size())
      fail("exact Jordan recurrence schedule row has the wrong block count");
    certificate.steps[n].reserve(indicial.blocks.size());
    for (std::size_t block_index = 0; block_index < indicial.blocks.size();
         ++block_index) {
      const auto& block = indicial.blocks[block_index];
      const auto& submitted = schedule[n][block_index];
      const auto expected_d_a =
          a_target + Rational(std::to_string(n)) - block.root.a;
      const auto expected_d_b = b_target - block.root.b;
      if (!(submitted.d_a == expected_d_a) ||
          !(submitted.d_b == expected_d_b))
        fail("submitted Jordan schedule offsets contradict the retained exact affine root at n=" +
             std::to_string(n) + ", block=" +
             std::to_string(block_index));
      const auto expected_kind = classify_step(expected_d_a, expected_d_b);
      if (submitted.kind != expected_kind)
        fail("submitted Jordan schedule case " +
             std::string(step_name(submitted.kind)) +
             " contradicts exact case " + step_name(expected_kind) +
             " at n=" + std::to_string(n) + ", block=" +
             std::to_string(block_index));

      ExactJordanScheduleStepCertificate step;
      step.taylor_index = static_cast<std::uint32_t>(n);
      step.block_index = static_cast<std::uint32_t>(block_index);
      step.kind = expected_kind;
      step.d_a = expected_d_a;
      step.d_b = expected_d_b;
      step.jordan_size = block.size();
      if (expected_kind == StepCase::Resonant) {
        step.resonant_jordan_chain_length = block.size();
        step.intrinsic_homogeneous_log_degree = block.size() - 1;
        certificate.contains_true_resonance = true;
        certificate.max_resonant_jordan_chain_length = std::max(
            certificate.max_resonant_jordan_chain_length, block.size());
      } else if (expected_kind == StepCase::Pseudo) {
        step.pseudo_epsilon_pole_depth = block.size();
        certificate.contains_pseudo = true;
        certificate.max_pseudo_epsilon_pole_depth = std::max(
            certificate.max_pseudo_epsilon_pole_depth, block.size());
      }
      certificate.steps[n].push_back(std::move(step));
    }
  }
  return certificate;
}

inline ExactJordanScheduleCertificate certify_exact_affine_jordan_schedule(
    const PreparedRecurrenceOperator<Rational>& prepared,
    const Rational& a_target, const Rational& b_target,
    const std::vector<std::vector<BlockStep<Rational>>>& schedule) {
  return certify_exact_affine_jordan_schedule(
      certify_exact_affine_jordan_operator(prepared), a_target, b_target,
      schedule);
}

}  // namespace diffexp2
