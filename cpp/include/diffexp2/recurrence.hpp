#pragma once

#include "diffexp2/scalar.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace diffexp2 {

constexpr std::int32_t kCompleteInfinity =
    std::numeric_limits<std::int32_t>::max();

enum class StepCase : std::uint8_t { Taylor, Pseudo, Resonant };

struct RecurrenceError : std::runtime_error {
  std::string id;
  std::int32_t frame_base = 0;
  std::int32_t shift = 0;

  RecurrenceError(std::string error_id, std::string detail,
                  std::int32_t fb = 0, std::int32_t eps_shift = 0)
      : std::runtime_error(std::move(detail)),
        id(std::move(error_id)),
        frame_base(fb),
        shift(eps_shift) {}
};

template <typename Scalar>
using Frame = std::vector<Scalar>;

template <typename Scalar>
using FrameBlock = std::vector<Frame<Scalar>>;  // [component][epsilon]

template <typename Scalar>
struct MatrixEntry {
  std::uint32_t row = 0;
  std::uint32_t col = 0;
  Scalar value;
};

template <typename Scalar>
struct MatrixShift {
  std::int32_t shift = 0;
  std::vector<MatrixEntry<Scalar>> entries;
};

template <typename Scalar>
struct ScalarShift {
  std::int32_t shift = 0;
  Scalar value;
};

template <typename Scalar>
struct RationalGroup {
  std::uint32_t denominator_index = 0;
  std::vector<MatrixShift<Scalar>> numerator;
};

template <typename Scalar>
struct PreparedLag {
  std::vector<MatrixShift<Scalar>> polynomial;
  std::vector<RationalGroup<Scalar>> rational;
  std::vector<std::int32_t> valuations;  // d*d, INF = structural zero
};

template <typename Scalar>
struct PreparedMatrix : PreparedLag<Scalar> {
  bool identity = false;
};

struct JordanBlock {
  std::vector<std::uint32_t> columns;
};

template <typename Scalar>
struct BlockStep {
  StepCase kind = StepCase::Taylor;
  Scalar d_a;
  Scalar d_b;
};

template <typename Scalar>
struct SourceData {
  // Flat tensors with epsilon fastest.
  std::vector<Scalar> frames;            // [n][log][component][eps]
  std::vector<std::int32_t> validity;    // [n][log][component]
  std::vector<std::uint8_t> present;     // [n][log]
};

template <typename Scalar>
struct RecurrenceProblem {
  std::uint32_t dimension = 0;
  std::uint32_t nmax = 0;
  std::uint32_t log_max = 0;
  std::int32_t frame_base = 0;
  std::uint32_t frame_width = 0;
  bool has_initial = true;
  bool adaptive_lower_frame_probe = false;

  Scalar a_target;
  Scalar b_target;
  // a_target+m, serialized after exact simplification so true zeros remain
  // exact. Entry zero corresponds to m=a_shift_min.
  std::int32_t a_shift_min = 0;
  std::vector<Scalar> a_shifts;

  std::vector<std::vector<ScalarShift<Scalar>>> d_lags;
  std::vector<PreparedLag<Scalar>> nhat_lags;  // index zero intentionally empty
  std::vector<std::vector<Scalar>> rational_denominators;
  std::optional<Scalar> d0_inverse_scalar;

  std::vector<JordanBlock> blocks;
  std::vector<std::vector<BlockStep<Scalar>>> schedule;  // [n][block]

  std::vector<Scalar> initial;              // [log][component][eps]
  std::vector<std::int32_t> initial_validity;  // [log][component]
  std::optional<SourceData<Scalar>> source;
  std::optional<PreparedMatrix<Scalar>> assembly_matrix;
  std::int32_t chop_digits = 0;
  bool return_u = true;
};

template <typename Scalar>
struct PseudoHit {
  std::uint32_t n = 0;
  std::vector<std::uint32_t> columns;
  Scalar delta_b;
  FrameBlock<Scalar> gamma_frames;
  std::vector<std::int32_t> gamma_validity;
};

template <typename Scalar>
struct RecurrenceResult {
  std::vector<Scalar> u;                  // [n][log][component][eps]
  std::vector<std::int32_t> validity;     // [n][log][component]
  std::vector<PseudoHit<Scalar>> hits;
  std::int32_t top_valid = kCompleteInfinity;
};

template <typename Scalar>
struct AssembledResult {
  std::int32_t min_power = 0;
  std::int32_t complete_max = 0;
  // [log][epsilon][n][component], matching SectorSeries Coeffs ordering.
  std::vector<Scalar> coefficients;
};

namespace detail {

inline std::int32_t valid_shift(std::int32_t value, std::int32_t shift,
                                std::int32_t frame_top) {
  if (value == kCompleteInfinity) return kCompleteInfinity;
  const auto shifted = static_cast<std::int64_t>(value) + shift;
  return static_cast<std::int32_t>(
      std::clamp<std::int64_t>(shifted,
          std::numeric_limits<std::int32_t>::min(), frame_top));
}

inline std::int32_t valid_min(std::int32_t a, std::int32_t b) {
  return std::min(a, b);
}

template <typename Scalar>
Frame<Scalar> zero_frame(std::size_t width) {
  return Frame<Scalar>(width, ScalarTraits<Scalar>::zero());
}

template <typename Scalar>
FrameBlock<Scalar> zero_block(std::size_t dimension, std::size_t width) {
  return FrameBlock<Scalar>(dimension, zero_frame<Scalar>(width));
}

template <typename Scalar>
std::int32_t frame_valuation(const Frame<Scalar>& frame,
                             std::int32_t frame_base) {
  for (std::size_t i = 0; i < frame.size(); ++i) {
    if (!ScalarTraits<Scalar>::is_zero(frame[i])) {
      return frame_base + static_cast<std::int32_t>(i);
    }
  }
  return kCompleteInfinity;
}

template <typename Scalar>
Frame<Scalar> shift_frame(const Frame<Scalar>& input, std::int32_t shift,
                          std::int32_t frame_base) {
  const auto width = input.size();
  auto out = zero_frame<Scalar>(width);
  if (shift == 0) return input;
  if (shift > 0) {
    if (static_cast<std::size_t>(shift) >= width) return out;
    for (std::size_t i = 0; i + static_cast<std::size_t>(shift) < width; ++i) {
      out[i + shift] = input[i];
    }
    return out;
  }
  const auto drop64 = -static_cast<std::int64_t>(shift);
  const auto drop = static_cast<std::size_t>(drop64);
  if (drop >= width) {
    for (const auto& value : input) {
      if (!ScalarTraits<Scalar>::is_zero(value)) {
        throw RecurrenceError("E4",
            "epsilon shift lies wholly below the work frame",
            frame_base, shift);
      }
    }
    return out;
  }
  for (std::size_t i = 0; i < drop; ++i) {
    if (!ScalarTraits<Scalar>::is_zero(input[i])) {
      throw RecurrenceError("E4",
          "epsilon shift would discard nonzero lower-frame content",
          frame_base, shift);
    }
  }
  for (std::size_t i = drop; i < width; ++i) out[i - drop] = input[i];
  return out;
}

template <typename Scalar>
FrameBlock<Scalar> shift_block(const FrameBlock<Scalar>& input,
                               std::int32_t shift,
                               std::int32_t frame_base) {
  FrameBlock<Scalar> out;
  out.reserve(input.size());
  for (const auto& frame : input) out.push_back(shift_frame(frame, shift, frame_base));
  return out;
}

template <typename Scalar>
void add_block_inplace(FrameBlock<Scalar>& target,
                       const FrameBlock<Scalar>& value) {
  for (std::size_t r = 0; r < target.size(); ++r) {
    for (std::size_t k = 0; k < target[r].size(); ++k) target[r][k] += value[r][k];
  }
}

template <typename Scalar>
void sub_block_inplace(FrameBlock<Scalar>& target,
                       const FrameBlock<Scalar>& value) {
  for (std::size_t r = 0; r < target.size(); ++r) {
    for (std::size_t k = 0; k < target[r].size(); ++k) target[r][k] -= value[r][k];
  }
}

template <typename Scalar>
FrameBlock<Scalar> scale_block(const FrameBlock<Scalar>& input,
                               const Scalar& scale) {
  auto out = input;
  for (auto& row : out) for (auto& value : row) value *= scale;
  return out;
}

template <typename Scalar>
FrameBlock<Scalar> apply_matrix(const MatrixShift<Scalar>& matrix,
                                const FrameBlock<Scalar>& input,
                                std::size_t dimension) {
  const auto width = input.front().size();
  auto out = zero_block<Scalar>(dimension, width);
  for (const auto& entry : matrix.entries) {
    for (std::size_t k = 0; k < width; ++k) {
      out[entry.row][k] += entry.value * input[entry.col][k];
    }
  }
  return out;
}

template <typename Scalar>
FrameBlock<Scalar> matrix_shift_product(const MatrixShift<Scalar>& matrix,
                                        const FrameBlock<Scalar>& input,
                                        const RecurrenceProblem<Scalar>& p) {
  if (p.adaptive_lower_frame_probe && matrix.shift < 0) {
    return shift_block(apply_matrix(matrix, input, p.dimension), matrix.shift,
                       p.frame_base);
  }
  return apply_matrix(matrix,
      shift_block(input, matrix.shift, p.frame_base), p.dimension);
}

template <typename Scalar>
FrameBlock<Scalar> divide_rational(const FrameBlock<Scalar>& rhs,
                                   const std::vector<Scalar>& q);

template <typename Scalar>
FrameBlock<Scalar> apply_prepared_matrix(
    const PreparedMatrix<Scalar>& matrix, const FrameBlock<Scalar>& input,
    const RecurrenceProblem<Scalar>& p) {
  if (matrix.identity) return input;
  std::vector<std::int32_t> input_valuations(p.dimension, kCompleteInfinity);
  for (std::uint32_t c = 0; c < p.dimension; ++c)
    input_valuations[c] = frame_valuation(input[c], p.frame_base);
  for (std::uint32_t r = 0; r < p.dimension; ++r) {
    for (std::uint32_t c = 0; c < p.dimension; ++c) {
      const auto mv = matrix.valuations[static_cast<std::size_t>(r) * p.dimension + c];
      const auto uv = input_valuations[c];
      if (mv != kCompleteInfinity && uv != kCompleteInfinity &&
          static_cast<std::int64_t>(mv) + uv < p.frame_base) {
        throw RecurrenceError("E4",
            "prepared matrix product would discard lower-epsilon content",
            p.frame_base, mv);
      }
    }
  }
  auto out = zero_block<Scalar>(p.dimension, p.frame_width);
  for (const auto& shifted_matrix : matrix.polynomial) {
    auto value = apply_matrix(shifted_matrix, input, p.dimension);
    value = shift_block(value, shifted_matrix.shift, p.frame_base);
    add_block_inplace(out, value);
  }
  for (const auto& group : matrix.rational) {
    auto numerator = zero_block<Scalar>(p.dimension, p.frame_width);
    for (const auto& shifted_matrix : group.numerator) {
      auto value = apply_matrix(shifted_matrix, input, p.dimension);
      value = shift_block(value, shifted_matrix.shift, p.frame_base);
      add_block_inplace(numerator, value);
    }
    add_block_inplace(out, divide_rational(numerator,
        p.rational_denominators.at(group.denominator_index)));
  }
  return out;
}

template <typename Scalar>
Frame<Scalar> divide_epsilon(const Frame<Scalar>& input,
                             std::int32_t frame_base) {
  if (!ScalarTraits<Scalar>::is_zero(input.front())) {
    throw RecurrenceError("E4",
        "division by epsilon would discard the lowest framed coefficient",
        frame_base, -1);
  }
  auto out = zero_frame<Scalar>(input.size());
  for (std::size_t i = 1; i < input.size(); ++i) out[i - 1] = input[i];
  return out;
}

template <typename Scalar>
Frame<Scalar> multiply_epsilon(const Frame<Scalar>& input) {
  auto out = zero_frame<Scalar>(input.size());
  for (std::size_t i = 1; i < input.size(); ++i) out[i] = input[i - 1];
  return out;
}

template <typename Scalar>
Frame<Scalar> convolve_frame(const Frame<Scalar>& a, const Frame<Scalar>& b,
                             std::int32_t frame_base) {
  const auto width = a.size();
  auto out = zero_frame<Scalar>(width);
  for (std::size_t ia = 0; ia < width; ++ia) {
    if (ScalarTraits<Scalar>::is_zero(a[ia])) continue;
    const auto pa = frame_base + static_cast<std::int32_t>(ia);
    for (std::size_t ib = 0; ib < width; ++ib) {
      if (ScalarTraits<Scalar>::is_zero(b[ib])) continue;
      const auto power = static_cast<std::int64_t>(pa) + frame_base +
                         static_cast<std::int64_t>(ib);
      if (power < frame_base) {
        throw RecurrenceError("E4",
            "framed convolution would discard nonzero lower-epsilon content",
            frame_base, power - frame_base);
      }
      const auto index = power - frame_base;
      if (index < static_cast<std::int64_t>(width)) {
        out[static_cast<std::size_t>(index)] += a[ia] * b[ib];
      }
    }
  }
  return out;
}

template <typename Scalar>
Frame<Scalar> invert_frame(const Frame<Scalar>& input,
                           std::int32_t frame_base) {
  const auto width = input.size();
  std::size_t first = width;
  for (std::size_t i = 0; i < width; ++i) {
    if (!ScalarTraits<Scalar>::is_zero(input[i])) { first = i; break; }
  }
  if (first == width) throw RecurrenceError("E5", "framed inverse of zero");
  const auto lead = frame_base + static_cast<std::int32_t>(first);
  const auto start64 = -static_cast<std::int64_t>(lead) - frame_base;
  if (start64 < 0 || start64 >= static_cast<std::int64_t>(width)) {
    throw RecurrenceError("E4", "framed inverse lies outside work frame",
                          frame_base, -lead);
  }
  const auto start = static_cast<std::size_t>(start64);
  auto out = zero_frame<Scalar>(width);
  const Scalar inv_lead = ScalarTraits<Scalar>::one() / input[first];
  out[start] = inv_lead;
  const auto max_relative = width - 1 - start;
  for (std::size_t k = 1; k <= max_relative; ++k) {
    Scalar sum = ScalarTraits<Scalar>::zero();
    for (std::size_t j = 1; j <= k; ++j) {
      if (first + j < width) sum += input[first + j] * out[start + k - j];
    }
    out[start + k] = -inv_lead * sum;
  }
  return out;
}

template <typename Scalar>
FrameBlock<Scalar> divide_rational(const FrameBlock<Scalar>& rhs,
                                   const std::vector<Scalar>& q) {
  if (q.empty() || ScalarTraits<Scalar>::is_zero(q.front())) {
    throw RecurrenceError("E5", "rational denominator has zero constant term");
  }
  const auto dimension = rhs.size();
  const auto width = rhs.front().size();
  auto out = zero_block<Scalar>(dimension, width);
  for (std::size_t r = 0; r < dimension; ++r) {
    out[r][0] = rhs[r][0] / q[0];
    for (std::size_t i = 1; i < width; ++i) {
      Scalar sum = ScalarTraits<Scalar>::zero();
      const auto last = std::min(i, q.size() - 1);
      for (std::size_t h = 1; h <= last; ++h) sum += q[h] * out[r][i - h];
      out[r][i] = (rhs[r][i] - sum) / q[0];
    }
  }
  return out;
}

template <typename Scalar>
Frame<Scalar> affine_divide(const Frame<Scalar>& rhs, const Scalar& d_a,
                            const Scalar& d_b, StepCase kind,
                            std::int32_t frame_base) {
  const auto width = rhs.size();
  auto out = zero_frame<Scalar>(width);
  if (kind == StepCase::Taylor) {
    out[0] = rhs[0] / d_a;
    for (std::size_t i = 1; i < width; ++i) {
      out[i] = (rhs[i] - d_b * out[i - 1]) / d_a;
    }
    return out;
  }
  if (kind == StepCase::Pseudo) {
    auto shifted = divide_epsilon(rhs, frame_base);
    for (std::size_t i = 0; i < width; ++i) out[i] = shifted[i] / d_b;
    return out;
  }
  throw RecurrenceError("E5", "resonant step passed to affine division");
}

template <typename Scalar>
FrameBlock<Scalar> apply_inv_d0(const FrameBlock<Scalar>& input,
                                const std::optional<Scalar>& inverse_scalar,
                                const Frame<Scalar>& inverse_frame,
                                std::int32_t frame_base) {
  auto out = input;
  if (inverse_scalar.has_value()) {
    for (auto& row : out) for (auto& value : row) value *= *inverse_scalar;
    return out;
  }
  for (auto& row : out) row = convolve_frame(row, inverse_frame, frame_base);
  return out;
}

template <typename Scalar>
FrameBlock<Scalar> solve_jordan(const FrameBlock<Scalar>& rhs,
                                const BlockStep<Scalar>& step,
                                const std::optional<Scalar>& inverse_scalar,
                                const Frame<Scalar>& inverse_frame,
                                std::int32_t frame_base) {
  const auto q = rhs.size();
  const auto width = rhs.front().size();
  if (step.kind == StepCase::Pseudo) {
    if (frame_base > -static_cast<std::int32_t>(q) ||
        frame_base + static_cast<std::int32_t>(width) - 1 < 1) {
      throw RecurrenceError("E4", "pseudo Jordan inverse exceeds work frame",
                            frame_base, -static_cast<std::int32_t>(q));
    }
    for (std::size_t j = 0; j < q; ++j) {
      const auto valuation = frame_valuation(rhs[j], frame_base);
      if (valuation != kCompleteInfinity &&
          static_cast<std::int64_t>(valuation) -
              static_cast<std::int64_t>(j + 1) < frame_base) {
        throw RecurrenceError("E4",
            "pseudo Jordan solve would discard lower-frame content",
            frame_base, -static_cast<std::int32_t>(j + 1));
      }
    }
  }
  auto z = zero_block<Scalar>(q, width);
  z[q - 1] = affine_divide(rhs[q - 1], step.d_a, step.d_b,
                           step.kind, frame_base);
  for (std::size_t rr = q - 1; rr-- > 0;) {
    auto combined = rhs[rr];
    for (std::size_t k = 0; k < width; ++k) combined[k] += z[rr + 1][k];
    z[rr] = affine_divide(combined, step.d_a, step.d_b,
                          step.kind, frame_base);
  }
  return apply_inv_d0(z, inverse_scalar, inverse_frame, frame_base);
}

}  // namespace detail

template <typename Scalar>
class RecurrenceSolver {
 public:
  explicit RecurrenceSolver(const RecurrenceProblem<Scalar>& problem)
      : p_(problem), d_(problem.dimension), nmax_(problem.nmax),
        logs_(0), width_(problem.frame_width),
        frame_top_(0), result_{} {
    validate_problem();
    logs_ = problem.log_max + 2;
    frame_top_ = static_cast<std::int32_t>(
        static_cast<std::int64_t>(problem.frame_base) +
        static_cast<std::int64_t>(problem.frame_width) - 1);
    const auto u_count = checked_tensor_size(
        {static_cast<std::size_t>(nmax_) + 1, logs_, d_, width_},
        "recurrence coefficient tensor");
    const auto validity_count = checked_tensor_size(
        {static_cast<std::size_t>(nmax_) + 1, logs_, d_},
        "recurrence validity tensor");
    result_.u.assign(u_count, ScalarTraits<Scalar>::zero());
    result_.validity.assign(validity_count, kCompleteInfinity);
    result_.top_valid = kCompleteInfinity;
    if (!p_.d0_inverse_scalar.has_value()) {
      Frame<Scalar> d0(width_, ScalarTraits<Scalar>::zero());
      for (const auto& item : p_.d_lags.front()) {
        const auto index = static_cast<std::int64_t>(item.shift) -
                           p_.frame_base;
        if (index < 0 || index >= static_cast<std::int64_t>(width_)) {
          throw RecurrenceError("E5", "d0 shift outside work frame");
        }
        d0[static_cast<std::size_t>(index)] = item.value;
      }
      inv_d0_frame_ = detail::invert_frame(d0, p_.frame_base);
    } else {
      inv_d0_frame_ = detail::zero_frame<Scalar>(width_);
    }
  }

  RecurrenceResult<Scalar> run() {
    seed_initial();
    const std::uint32_t n0 = p_.has_initial ? 1 : 0;
    for (std::uint32_t n = n0; n <= nmax_; ++n) step_n(n);
    result_.top_valid = kCompleteInfinity;
    for (const auto value : result_.validity) {
      if (value != kCompleteInfinity) result_.top_valid = std::min(result_.top_valid, value);
    }
    if (result_.top_valid == kCompleteInfinity) result_.top_valid = frame_top_;
    return std::move(result_);
  }

 private:
  std::size_t u_index(std::uint32_t n, std::uint32_t log,
                      std::uint32_t component, std::uint32_t eps) const {
    return ((((static_cast<std::size_t>(n) * logs_) + log) * d_ + component) *
            width_) + eps;
  }
  std::size_t v_index(std::uint32_t n, std::uint32_t log,
                      std::uint32_t component) const {
    return ((static_cast<std::size_t>(n) * logs_ + log) * d_ + component);
  }

  FrameBlock<Scalar> get_block(std::uint32_t n, std::uint32_t log) const {
    auto out = detail::zero_block<Scalar>(d_, width_);
    for (std::uint32_t r = 0; r < d_; ++r) {
      for (std::uint32_t k = 0; k < width_; ++k) out[r][k] = result_.u[u_index(n, log, r, k)];
    }
    return out;
  }
  void set_block(std::uint32_t n, std::uint32_t log,
                 const FrameBlock<Scalar>& block) {
    for (std::uint32_t r = 0; r < d_; ++r) {
      for (std::uint32_t k = 0; k < width_; ++k) result_.u[u_index(n, log, r, k)] = block[r][k];
    }
  }
  std::vector<std::int32_t> get_validity(std::uint32_t n,
                                         std::uint32_t log) const {
    std::vector<std::int32_t> out(d_);
    for (std::uint32_t r = 0; r < d_; ++r) out[r] = result_.validity[v_index(n, log, r)];
    return out;
  }
  void set_validity(std::uint32_t n, std::uint32_t log,
                    const std::vector<std::int32_t>& values) {
    for (std::uint32_t r = 0; r < d_; ++r) result_.validity[v_index(n, log, r)] = values[r];
  }

  const Scalar& a_shift(std::int32_t m) const {
    const auto index = m - p_.a_shift_min;
    if (index < 0 || index >= static_cast<std::int32_t>(p_.a_shifts.size())) {
      throw RecurrenceError("E5", "missing exact a-target shift schedule");
    }
    return p_.a_shifts[index];
  }

  void seed_initial() {
    if (!p_.has_initial) return;
    for (std::uint32_t log = 0; log <= p_.log_max; ++log) {
      for (std::uint32_t r = 0; r < d_; ++r) {
        for (std::uint32_t k = 0; k < width_; ++k) {
          const auto source = ((static_cast<std::size_t>(log) * d_ + r) * width_ + k);
          result_.u[u_index(0, log, r, k)] = p_.initial[source];
        }
        result_.validity[v_index(0, log, r)] =
            p_.initial_validity[static_cast<std::size_t>(log) * d_ + r];
      }
    }
  }

  FrameBlock<Scalar> apply_nhat_lag(const PreparedLag<Scalar>& lag,
                                    const FrameBlock<Scalar>& input) const {
    auto out = detail::zero_block<Scalar>(d_, width_);
    for (const auto& matrix : lag.polynomial) {
      detail::add_block_inplace(out, detail::matrix_shift_product(matrix, input, p_));
    }
    for (const auto& group : lag.rational) {
      auto numerator = detail::zero_block<Scalar>(d_, width_);
      for (const auto& matrix : group.numerator) {
        detail::add_block_inplace(numerator,
            detail::matrix_shift_product(matrix, input, p_));
      }
      detail::add_block_inplace(out,
          detail::divide_rational(numerator,
              p_.rational_denominators.at(group.denominator_index)));
    }
    return out;
  }

  void update_nhat_validity(std::vector<std::int32_t>& acc,
                            const PreparedLag<Scalar>& lag,
                            const std::vector<std::int32_t>& input) const {
    for (std::uint32_t r = 0; r < d_; ++r) {
      std::int32_t row_min = kCompleteInfinity;
      for (std::uint32_t c = 0; c < d_; ++c) {
        const auto valuation = lag.valuations[static_cast<std::size_t>(r) * d_ + c];
        if (valuation != kCompleteInfinity) {
          row_min = std::min(row_min,
              detail::valid_shift(input[c], valuation, frame_top_));
        }
      }
      acc[r] = std::min(acc[r], row_min);
    }
  }

  std::optional<std::pair<FrameBlock<Scalar>, std::vector<std::int32_t>>>
  source_at(std::uint32_t n, std::uint32_t log) const {
    if (!p_.source.has_value()) return std::nullopt;
    const auto& src = *p_.source;
    const auto pindex = static_cast<std::size_t>(n) * (p_.log_max + 1) + log;
    if (!src.present[pindex]) return std::nullopt;
    auto frames = detail::zero_block<Scalar>(d_, width_);
    std::vector<std::int32_t> validity(d_);
    for (std::uint32_t r = 0; r < d_; ++r) {
      validity[r] = src.validity[pindex * d_ + r];
      for (std::uint32_t k = 0; k < width_; ++k) {
        frames[r][k] = src.frames[(pindex * d_ + r) * width_ + k];
      }
    }
    return std::make_pair(std::move(frames), std::move(validity));
  }

  std::pair<FrameBlock<Scalar>, std::vector<std::int32_t>>
  build_rhs(std::uint32_t n, std::uint32_t log) const {
    auto acc = detail::zero_block<Scalar>(d_, width_);
    std::vector<std::int32_t> acc_valid(d_, kCompleteInfinity);

    // Polynomial Nhat terms and unfused rational groups. The JSON codec may
    // fuse equal denominators across lags in a later optimization; this form
    // is coefficient-identical and keeps every per-contribution underflow
    // witness before any cancellation.
    const auto max_nhat = std::min<std::uint32_t>(n, p_.nhat_lags.size() - 1);
    for (std::uint32_t j = 1; j <= max_nhat; ++j) {
      const auto input = get_block(n - j, log);
      detail::add_block_inplace(acc, apply_nhat_lag(p_.nhat_lags[j], input));
      update_nhat_validity(acc_valid, p_.nhat_lags[j], get_validity(n - j, log));
    }

    const auto max_d = std::min<std::uint32_t>(n, p_.d_lags.size() - 1);
    for (std::uint32_t j = 1; j <= max_d; ++j) {
      const auto input = get_block(n - j, log);
      const auto above = get_block(n - j, log + 1);
      auto term = detail::scale_block(input, a_shift(static_cast<std::int32_t>(n - j)));
      detail::add_block_inplace(term,
          detail::shift_block(detail::scale_block(input, p_.b_target), 1, p_.frame_base));
      detail::add_block_inplace(term,
          detail::shift_block(above, 1, p_.frame_base));

      const auto input_valid = get_validity(n - j, log);
      const auto above_valid = get_validity(n - j, log + 1);
      std::vector<std::int32_t> term_valid(d_, kCompleteInfinity);
      for (std::uint32_t r = 0; r < d_; ++r) {
        const auto va = ScalarTraits<Scalar>::is_zero(
                            a_shift(static_cast<std::int32_t>(n - j)))
                            ? kCompleteInfinity : input_valid[r];
        const auto vb = ScalarTraits<Scalar>::is_zero(p_.b_target)
                            ? kCompleteInfinity
                            : detail::valid_shift(input_valid[r], 1, frame_top_);
        const auto vl = detail::valid_shift(above_valid[r], 1, frame_top_);
        term_valid[r] = std::min({va, vb, vl});
      }

      for (const auto& scalar_shift : p_.d_lags[j]) {
        auto shifted = detail::shift_block(term, scalar_shift.shift, p_.frame_base);
        shifted = detail::scale_block(shifted, scalar_shift.value);
        detail::sub_block_inplace(acc, shifted);
        for (std::uint32_t r = 0; r < d_; ++r) {
          acc_valid[r] = std::min(acc_valid[r],
              detail::valid_shift(term_valid[r], scalar_shift.shift, frame_top_));
        }
      }
    }

    if (p_.source.has_value()) {
      for (std::uint32_t j = 0; j <= max_d; ++j) {
        if (n < j) continue;
        auto source = source_at(n - j, log);
        if (!source.has_value()) continue;
        for (const auto& scalar_shift : p_.d_lags[j]) {
          auto shifted = detail::shift_block(source->first,
                                              scalar_shift.shift, p_.frame_base);
          shifted = detail::scale_block(shifted, scalar_shift.value);
          detail::add_block_inplace(acc, shifted);
          for (std::uint32_t r = 0; r < d_; ++r) {
            acc_valid[r] = std::min(acc_valid[r],
                detail::valid_shift(source->second[r], scalar_shift.shift, frame_top_));
          }
        }
      }
    }
    return {std::move(acc), std::move(acc_valid)};
  }

  void solve_nonresonant(std::uint32_t n, std::uint32_t log,
                         const FrameBlock<Scalar>& r_block,
                         const std::vector<std::int32_t>& r_valid) {
    auto rhs = r_block;
    auto rhs_valid = r_valid;
    const auto above = get_block(n, log + 1);
    const auto above_valid = get_validity(n, log + 1);
    const auto eps_above = detail::shift_block(
        detail::shift_block(above, 1, p_.frame_base), 0, p_.frame_base);
    for (const auto& scalar_shift : p_.d_lags.front()) {
      auto shifted = detail::shift_block(eps_above, scalar_shift.shift, p_.frame_base);
      shifted = detail::scale_block(shifted, scalar_shift.value);
      detail::sub_block_inplace(rhs, shifted);
      for (std::uint32_t r = 0; r < d_; ++r) {
        rhs_valid[r] = std::min(rhs_valid[r],
            detail::valid_shift(above_valid[r], 1 + scalar_shift.shift, frame_top_));
      }
    }

    auto current = get_block(n, log);
    auto current_valid = get_validity(n, log);
    for (std::size_t bi = 0; bi < p_.blocks.size(); ++bi) {
      const auto& block = p_.blocks[bi];
      const auto& step = p_.schedule[n][bi];
      if (step.kind == StepCase::Resonant) continue;
      FrameBlock<Scalar> block_rhs;
      std::vector<std::int32_t> block_valid;
      for (const auto col : block.columns) {
        block_rhs.push_back(rhs[col]);
        block_valid.push_back(rhs_valid[col]);
      }
      const auto solved = detail::solve_jordan(block_rhs, step,
          p_.d0_inverse_scalar, inv_d0_frame_, p_.frame_base);
      for (std::size_t r = 0; r < block.columns.size(); ++r) {
        current[block.columns[r]] = solved[r];
        std::int32_t solved_valid = kCompleteInfinity;
        if (step.kind == StepCase::Taylor) {
          for (std::size_t m = r; m < block_valid.size(); ++m)
            solved_valid = std::min(solved_valid, block_valid[m]);
        } else {
          for (std::size_t m = r; m < block_valid.size(); ++m) {
            solved_valid = std::min(solved_valid,
                detail::valid_shift(block_valid[m],
                    -static_cast<std::int32_t>(m - r + 1), frame_top_));
          }
        }
        current_valid[block.columns[r]] = solved_valid;
      }
      if (step.kind == StepCase::Pseudo && log == 0) {
        PseudoHit<Scalar> hit;
        hit.n = n;
        hit.columns = block.columns;
        hit.delta_b = step.d_b;
        for (const auto col : block.columns) {
          hit.gamma_frames.push_back(current[col]);
          hit.gamma_validity.push_back(current_valid[col]);
        }
        result_.hits.push_back(std::move(hit));
      }
    }
    set_block(n, log, current);
    set_validity(n, log, current_valid);
  }

  void solve_resonant(std::uint32_t n,
                      const std::vector<FrameBlock<Scalar>>& r_blocks,
                      const std::vector<std::vector<std::int32_t>>& r_valid) {
    for (std::size_t bi = 0; bi < p_.blocks.size(); ++bi) {
      const auto& block = p_.blocks[bi];
      if (p_.schedule[n][bi].kind != StepCase::Resonant) continue;
      const auto q = block.columns.size();
      std::vector<FrameBlock<Scalar>> rt;
      std::vector<std::vector<std::int32_t>> rtv;
      for (std::uint32_t log = 0; log <= p_.log_max; ++log) {
        FrameBlock<Scalar> selected;
        std::vector<std::int32_t> selected_valid;
        for (const auto col : block.columns) {
          selected.push_back(r_blocks[log][col]);
          selected_valid.push_back(r_valid[log][col]);
        }
        rt.push_back(detail::apply_inv_d0(selected, p_.d0_inverse_scalar,
                                          inv_d0_frame_, p_.frame_base));
        rtv.push_back(std::move(selected_valid));
      }

      std::vector<std::vector<Frame<Scalar>>> assigned(
          p_.log_max + 2,
          std::vector<Frame<Scalar>>(q, detail::zero_frame<Scalar>(width_)));
      std::vector<std::vector<std::int32_t>> assigned_valid(
          p_.log_max + 2,
          std::vector<std::int32_t>(q, kCompleteInfinity));
      std::vector<std::vector<std::uint8_t>> present(
          p_.log_max + 2, std::vector<std::uint8_t>(q, 0));

      for (std::uint32_t log = 0; log < p_.log_max; ++log) {
        assigned[log + 1][q - 1] =
            detail::divide_epsilon(rt[log][q - 1], p_.frame_base);
        assigned_valid[log + 1][q - 1] =
            detail::valid_shift(rtv[log][q - 1], -1, frame_top_);
        present[log + 1][q - 1] = 1;
      }
      for (std::int32_t log = static_cast<std::int32_t>(p_.log_max);
           log >= 0; --log) {
        for (std::size_t r = 0; r + 1 < q; ++r) {
          if (present[log][r + 1]) continue;
          auto upper = assigned[log + 1][r];
          auto value = detail::multiply_epsilon(upper);
          for (std::size_t k = 0; k < width_; ++k) value[k] -= rt[log][r][k];
          assigned[log][r + 1] = std::move(value);
          assigned_valid[log][r + 1] = std::min(
              detail::valid_shift(assigned_valid[log + 1][r], 1, frame_top_),
              rtv[log][r]);
          present[log][r + 1] = 1;
        }
      }

      for (std::uint32_t log = 0; log <= p_.log_max; ++log) {
        auto current = get_block(n, log);
        auto current_valid = get_validity(n, log);
        for (std::size_t r = 0; r < q; ++r) {
          if (!present[log][r]) continue;
          current[block.columns[r]] = assigned[log][r];
          current_valid[block.columns[r]] = assigned_valid[log][r];
        }
        set_block(n, log, current);
        set_validity(n, log, current_valid);
      }
    }
  }

  void step_n(std::uint32_t n) {
    std::vector<FrameBlock<Scalar>> r_blocks;
    std::vector<std::vector<std::int32_t>> r_valid;
    r_blocks.reserve(p_.log_max + 1);
    r_valid.reserve(p_.log_max + 1);
    for (std::uint32_t log = 0; log <= p_.log_max; ++log) {
      auto rhs = build_rhs(n, log);
      r_blocks.push_back(std::move(rhs.first));
      r_valid.push_back(std::move(rhs.second));
    }
    for (std::int32_t log = static_cast<std::int32_t>(p_.log_max);
         log >= 0; --log) {
      solve_nonresonant(n, static_cast<std::uint32_t>(log),
                        r_blocks[log], r_valid[log]);
    }
    solve_resonant(n, r_blocks, r_valid);
  }

  void validate_problem() const {
    if (d_ == 0 || width_ == 0) throw RecurrenceError("E5", "zero recurrence dimension");
    if (p_.log_max > std::numeric_limits<std::uint32_t>::max() - 2)
      throw RecurrenceError("E5", "logarithmic depth exceeds uint32 range");
    const auto frame_top = static_cast<std::int64_t>(p_.frame_base) +
                           static_cast<std::int64_t>(width_) - 1;
    if (frame_top < std::numeric_limits<std::int32_t>::min() ||
        frame_top > std::numeric_limits<std::int32_t>::max())
      throw RecurrenceError("E5", "epsilon work frame exceeds int32 range");
    (void)checked_tensor_size(
        {static_cast<std::size_t>(nmax_) + 1,
         static_cast<std::size_t>(p_.log_max) + 2, d_, width_},
        "recurrence coefficient tensor");
    if (p_.d_lags.empty() || p_.nhat_lags.empty())
      throw RecurrenceError("E5", "missing prepared recurrence lags");
    if (p_.d_lags.front().empty())
      throw RecurrenceError("E5", "prepared d0 lag is empty");
    if (p_.a_shift_min > 0 ||
        static_cast<std::int64_t>(p_.a_shift_min) +
            static_cast<std::int64_t>(p_.a_shifts.size()) <=
            static_cast<std::int64_t>(nmax_))
      throw RecurrenceError("E5", "incomplete exact a-target shift schedule");
    if (p_.schedule.size() != nmax_ + 1)
      throw RecurrenceError("E5", "invalid block-step schedule height");
    for (const auto& row : p_.schedule) {
      if (row.size() != p_.blocks.size())
        throw RecurrenceError("E5", "invalid block-step schedule width");
    }
    if (p_.has_initial) {
      const auto frame_count = checked_tensor_size(
          {static_cast<std::size_t>(p_.log_max) + 1, d_},
          "initial validity tensor");
      const auto initial_count = checked_tensor_size(
          {frame_count, width_}, "initial recurrence tensor");
      if (p_.initial.size() != initial_count ||
          p_.initial_validity.size() != frame_count)
        throw RecurrenceError("E5", "malformed initial recurrence tensor");
    }
    if (p_.source.has_value()) {
      const auto points = checked_tensor_size(
          {static_cast<std::size_t>(nmax_) + 1,
           static_cast<std::size_t>(p_.log_max) + 1},
          "source point tensor");
      const auto source_validity_count = checked_tensor_size(
          {points, d_}, "source validity tensor");
      const auto source_frame_count = checked_tensor_size(
          {source_validity_count, width_}, "source recurrence tensor");
      if (p_.source->present.size() != points ||
          p_.source->validity.size() != source_validity_count ||
          p_.source->frames.size() != source_frame_count)
        throw RecurrenceError("E5", "malformed source recurrence tensor");
    }
    const auto matrix_size = checked_tensor_size({d_, d_},
                                                  "matrix valuation tensor");
    const auto validate_shift = [](std::int32_t shift) {
      if (shift == std::numeric_limits<std::int32_t>::min())
        throw RecurrenceError("E5", "epsilon shift is outside the supported range");
    };
    const auto validate_matrix_shifts = [&](const auto& owner) {
      for (const auto& matrix : owner.polynomial) validate_shift(matrix.shift);
      for (const auto& group : owner.rational)
        for (const auto& matrix : group.numerator) validate_shift(matrix.shift);
    };
    const auto validate_groups = [&](const auto& owner) {
      for (const auto& group : owner.rational) {
        if (group.denominator_index >= p_.rational_denominators.size())
          throw RecurrenceError("E5", "rational group denominator index is out of range");
      }
    };
    for (const auto& denominator : p_.rational_denominators) {
      if (denominator.empty() || ScalarTraits<Scalar>::is_zero(denominator.front()))
        throw RecurrenceError("E5", "rational denominator has zero constant term");
    }
    for (const auto& lag : p_.d_lags)
      for (const auto& item : lag) validate_shift(item.shift);
    for (const auto& lag : p_.nhat_lags) {
      if (lag.valuations.size() != matrix_size)
        throw RecurrenceError("E5", "malformed Nhat valuation tensor");
      validate_groups(lag);
      validate_matrix_shifts(lag);
    }
    if (p_.assembly_matrix.has_value()) {
      if (p_.assembly_matrix->valuations.size() != matrix_size)
        throw RecurrenceError("E5", "malformed assembly valuation tensor");
      validate_groups(*p_.assembly_matrix);
      validate_matrix_shifts(*p_.assembly_matrix);
    } else if (!p_.return_u) {
      throw RecurrenceError("E5", "request suppresses U without compiled assembly");
    }
    std::vector<std::uint8_t> covered(d_, 0);
    for (const auto& block : p_.blocks) {
      if (block.columns.empty()) throw RecurrenceError("E5", "empty Jordan block");
      for (const auto col : block.columns) {
        if (col >= d_)
          throw RecurrenceError("E5", "Jordan column outside system");
        if (covered[col])
          throw RecurrenceError("E5", "Jordan blocks contain a duplicate column");
        covered[col] = 1;
      }
    }
    if (std::any_of(covered.begin(), covered.end(),
                    [](std::uint8_t value) { return value == 0; }))
      throw RecurrenceError("E5", "Jordan blocks do not partition the system");
  }

  static std::size_t checked_tensor_size(
      std::initializer_list<std::size_t> dimensions, const char* label) {
    constexpr std::size_t kMaxTensorElements = 100'000'000;
    std::size_t count = 1;
    for (const auto dimension : dimensions) {
      if (dimension != 0 &&
          count > std::numeric_limits<std::size_t>::max() / dimension)
        throw RecurrenceError("E5", std::string(label) + " size overflows");
      count *= dimension;
      if (count > kMaxTensorElements)
        throw RecurrenceError("E5", std::string(label) + " is unreasonably large");
    }
    return count;
  }

  const RecurrenceProblem<Scalar>& p_;
  std::uint32_t d_, nmax_, logs_, width_;
  std::int32_t frame_top_;
  RecurrenceResult<Scalar> result_;
  Frame<Scalar> inv_d0_frame_;
};

template <typename Scalar>
AssembledResult<Scalar> assemble_recurrence(
    const RecurrenceProblem<Scalar>& p, const RecurrenceResult<Scalar>& result) {
  if (!p.assembly_matrix.has_value())
    throw RecurrenceError("E5", "compiled assembly matrix is missing");
  const auto& matrix = *p.assembly_matrix;
  const std::uint32_t logs = p.log_max + 2;
  const std::int32_t frame_top =
      p.frame_base + static_cast<std::int32_t>(p.frame_width) - 1;
  const auto u_index = [&](std::uint32_t n, std::uint32_t log,
                           std::uint32_t component, std::uint32_t eps) {
    return ((((static_cast<std::size_t>(n) * logs) + log) * p.dimension +
             component) * p.frame_width + eps);
  };
  const auto v_index = [&](std::uint32_t n, std::uint32_t log,
                           std::uint32_t component) {
    return ((static_cast<std::size_t>(n) * logs + log) * p.dimension + component);
  };

  // Store transformed frames only for the physical log levels 0..P.
  std::vector<Scalar> transformed(
      static_cast<std::size_t>(p.nmax + 1) * (p.log_max + 1) *
          p.dimension * p.frame_width,
      ScalarTraits<Scalar>::zero());
  std::vector<std::int32_t> transformed_validity(
      static_cast<std::size_t>(p.nmax + 1) * (p.log_max + 1) * p.dimension,
      kCompleteInfinity);
  const auto t_index = [&](std::uint32_t n, std::uint32_t log,
                           std::uint32_t component, std::uint32_t eps) {
    return ((((static_cast<std::size_t>(n) * (p.log_max + 1)) + log) *
             p.dimension + component) * p.frame_width + eps);
  };
  const auto tv_index = [&](std::uint32_t n, std::uint32_t log,
                            std::uint32_t component) {
    return ((static_cast<std::size_t>(n) * (p.log_max + 1) + log) *
            p.dimension + component);
  };

  std::int32_t complete_max = kCompleteInfinity;
  for (std::uint32_t n = 0; n <= p.nmax; ++n) {
    for (std::uint32_t log = 0; log <= p.log_max; ++log) {
      auto input = detail::zero_block<Scalar>(p.dimension, p.frame_width);
      for (std::uint32_t c = 0; c < p.dimension; ++c)
        for (std::uint32_t k = 0; k < p.frame_width; ++k)
          input[c][k] = result.u[u_index(n, log, c, k)];
      const auto output = detail::apply_prepared_matrix(matrix, input, p);
      for (std::uint32_t r = 0; r < p.dimension; ++r) {
        for (std::uint32_t k = 0; k < p.frame_width; ++k)
          transformed[t_index(n, log, r, k)] = output[r][k];
        std::int32_t row_valid = kCompleteInfinity;
        for (std::uint32_t c = 0; c < p.dimension; ++c) {
          const auto mv = matrix.valuations[static_cast<std::size_t>(r) * p.dimension + c];
          if (mv != kCompleteInfinity) {
            row_valid = std::min(row_valid, detail::valid_shift(
                result.validity[v_index(n, log, c)], mv, frame_top));
          }
        }
        transformed_validity[tv_index(n, log, r)] = row_valid;
        if (row_valid != kCompleteInfinity) complete_max = std::min(complete_max, row_valid);
      }
    }
  }
  if (complete_max == kCompleteInfinity) complete_max = frame_top;

  std::optional<std::uint32_t> first_nonzero;
  for (std::uint32_t k = 0; k < p.frame_width && !first_nonzero.has_value(); ++k) {
    for (std::uint32_t n = 0; n <= p.nmax && !first_nonzero.has_value(); ++n) {
      for (std::uint32_t log = 0; log <= p.log_max && !first_nonzero.has_value(); ++log) {
        for (std::uint32_t r = 0; r < p.dimension; ++r) {
          if (!ScalarTraits<Scalar>::certified_zero(
                  transformed[t_index(n, log, r, k)], p.chop_digits)) {
            first_nonzero = k;
            break;
          }
        }
      }
    }
  }
  const std::int32_t min_power = first_nonzero.has_value()
      ? p.frame_base + static_cast<std::int32_t>(*first_nonzero)
      : complete_max;
  if (first_nonzero.has_value() && complete_max < min_power) {
    throw RecurrenceError("E6",
        "compiled assembly exhausted completeness below stored content",
        p.frame_base, 0);
  }

  AssembledResult<Scalar> assembled;
  assembled.min_power = min_power;
  assembled.complete_max = complete_max;
  const auto count = static_cast<std::size_t>(complete_max - min_power + 1);
  assembled.coefficients.reserve(
      static_cast<std::size_t>(p.log_max + 1) * count *
      (p.nmax + 1) * p.dimension);
  for (std::uint32_t log = 0; log <= p.log_max; ++log) {
    for (std::int32_t power = min_power; power <= complete_max; ++power) {
      const auto k = power - p.frame_base;
      for (std::uint32_t n = 0; n <= p.nmax; ++n) {
        for (std::uint32_t r = 0; r < p.dimension; ++r) {
          Scalar value = ScalarTraits<Scalar>::zero();
          if (k >= 0 && k < static_cast<std::int32_t>(p.frame_width))
            value = transformed[t_index(n, log, r, static_cast<std::uint32_t>(k))];
          assembled.coefficients.push_back(
              ScalarTraits<Scalar>::canonicalized(value, p.chop_digits));
        }
      }
    }
  }
  return assembled;
}

}  // namespace diffexp2
