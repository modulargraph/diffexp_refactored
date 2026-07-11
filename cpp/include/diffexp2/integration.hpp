#pragma once

#include "diffexp2/local_solution.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace diffexp2 {

// A finite Laurent frame with the same honesty contract as EpsSeries.m:
// below min_power is certified structural zero, while anything above
// complete_max is unknown.  In particular add takes the intersection of the
// complete upper windows and multiplication uses the two-sided Cauchy bound.
template <typename Scalar>
class EpsilonFrame {
 public:
  EpsilonFrame(EpsilonWindow window, std::vector<Scalar> coefficients)
      : window_(window), coefficients_(std::move(coefficients)) {
    if (coefficients_.size() != window_.width())
      throw std::invalid_argument(
          "epsilon-frame coefficient count disagrees with its honest window");
  }

  EpsilonFrame(std::int32_t min_power, std::vector<Scalar> coefficients)
      : window_(window_from_width(min_power, coefficients.size())),
        coefficients_(std::move(coefficients)) {}

  static EpsilonFrame zero(std::int32_t complete_max) {
    return EpsilonFrame({complete_max, complete_max},
                        {ScalarTraits<Scalar>::zero()});
  }

  [[nodiscard]] const EpsilonWindow& window() const { return window_; }
  [[nodiscard]] std::int32_t min_power() const { return window_.min_power; }
  [[nodiscard]] std::int32_t complete_max() const {
    return window_.complete_max;
  }
  [[nodiscard]] const std::vector<Scalar>& coefficients() const {
    return coefficients_;
  }

  [[nodiscard]] Scalar coefficient(std::int32_t power) const {
    if (power < min_power()) return ScalarTraits<Scalar>::zero();
    if (power > complete_max())
      throw std::out_of_range(
          "epsilon coefficient requested above the complete window");
    return coefficients_.at(static_cast<std::size_t>(
        static_cast<std::int64_t>(power) - min_power()));
  }

  [[nodiscard]] EpsilonFrame shifted(std::int32_t amount) const {
    return EpsilonFrame(
        {checked_power(static_cast<std::int64_t>(min_power()) + amount,
                       "shifted epsilon minimum"),
         checked_power(static_cast<std::int64_t>(complete_max()) + amount,
                       "shifted epsilon complete maximum")},
        coefficients_);
  }

  [[nodiscard]] EpsilonFrame scaled(const Scalar& factor) const {
    auto result = coefficients_;
    for (auto& value : result) value *= factor;
    return EpsilonFrame(window_, std::move(result));
  }

  [[nodiscard]] EpsilonFrame truncated(std::int32_t new_complete_max) const {
    if (new_complete_max > complete_max() || new_complete_max < min_power())
      throw std::invalid_argument(
          "epsilon truncation must remain inside the complete window");
    auto result = coefficients_;
    result.resize(static_cast<std::size_t>(
        static_cast<std::int64_t>(new_complete_max) - min_power() + 1));
    return EpsilonFrame({min_power(), new_complete_max}, std::move(result));
  }

  friend EpsilonFrame operator+(const EpsilonFrame& left,
                                const EpsilonFrame& right) {
    const auto minimum = std::min(left.min_power(), right.min_power());
    const auto maximum =
        std::min(left.complete_max(), right.complete_max());
    std::vector<Scalar> result;
    result.reserve(EpsilonWindow{minimum, maximum}.width());
    for (std::int64_t k = minimum; k <= maximum; ++k) {
      const auto power = static_cast<std::int32_t>(k);
      result.push_back(left.coefficient(power) + right.coefficient(power));
    }
    return EpsilonFrame({minimum, maximum}, std::move(result));
  }

  friend EpsilonFrame operator-(const EpsilonFrame& value) {
    auto result = value.coefficients_;
    for (auto& coefficient : result) coefficient = -coefficient;
    return EpsilonFrame(value.window_, std::move(result));
  }

  friend EpsilonFrame operator-(const EpsilonFrame& left,
                                const EpsilonFrame& right) {
    return left + (-right);
  }

  friend EpsilonFrame operator*(const EpsilonFrame& left,
                                const EpsilonFrame& right) {
    const auto minimum = checked_power(
        static_cast<std::int64_t>(left.min_power()) + right.min_power(),
        "epsilon-product minimum");
    const auto maximum = checked_power(
        std::min(static_cast<std::int64_t>(left.complete_max()) +
                     right.min_power(),
                 static_cast<std::int64_t>(right.complete_max()) +
                     left.min_power()),
        "epsilon-product complete maximum");
    std::vector<Scalar> result(EpsilonWindow{minimum, maximum}.width(),
                               ScalarTraits<Scalar>::zero());
    for (std::int64_t k64 = minimum; k64 <= maximum; ++k64) {
      const auto first64 = std::max<std::int64_t>(
          left.min_power(), k64 - right.complete_max());
      const auto last64 = std::min<std::int64_t>(
          left.complete_max(), k64 - right.min_power());
      auto value = ScalarTraits<Scalar>::zero();
      for (std::int64_t i64 = first64; i64 <= last64; ++i64)
        value += left.coefficient(static_cast<std::int32_t>(i64)) *
                 right.coefficient(static_cast<std::int32_t>(k64 - i64));
      result[static_cast<std::size_t>(k64 - minimum)] = std::move(value);
    }
    return EpsilonFrame({minimum, maximum}, std::move(result));
  }

 private:
  static std::int32_t checked_power(std::int64_t value, const char* label) {
    if (value < std::numeric_limits<std::int32_t>::min() ||
        value > std::numeric_limits<std::int32_t>::max())
      throw std::overflow_error(std::string(label) + " exceeds int32 range");
    return static_cast<std::int32_t>(value);
  }

  static EpsilonWindow window_from_width(std::int32_t min_power,
                                         std::size_t width) {
    if (width == 0)
      throw std::invalid_argument("an epsilon frame cannot be empty");
    if (width - 1 >
        static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
      throw std::overflow_error("epsilon-frame width exceeds int64 range");
    return {min_power,
            checked_power(static_cast<std::int64_t>(min_power) +
                              static_cast<std::int64_t>(width - 1),
                          "epsilon-frame complete maximum")};
  }

  EpsilonWindow window_;
  std::vector<Scalar> coefficients_;
};

enum class NativeIntegrationErrorCode : std::uint8_t {
  DivergentEndpoint,
  MissingBranchPrescription,
  IncompleteTaylorWindow,
  UnsupportedExactTag,
  UnsupportedSymbolicRegulator,
  InvalidInterval,
  UncertifiedCancellation
};

class NativeIntegrationError : public std::runtime_error {
 public:
  NativeIntegrationError(NativeIntegrationErrorCode error_code,
                         std::string error_id, std::string detail)
      : std::runtime_error(std::move(detail)),
        code(error_code),
        id(std::move(error_id)) {}

  NativeIntegrationErrorCode code;
  std::string id;
  std::string absolute_power;
  std::uint32_t log_power = 0;
  std::int32_t epsilon_power = 0;
  std::uint32_t component = 0;
};

// The cell power and alpha0=m+1 are both retained as exact descriptors.  The
// latter avoids ever deciding a resonant denominator by inspecting a numeric
// ball.  The local-sector adapter below constructs both in the rational
// domain; future algebraic preparation may supply an equally exact pair.
struct SectorMonomialTag {
  ExactScalarDescriptor m;
  ExactScalarDescriptor alpha0;
  ExactScalarDescriptor b;
  std::uint32_t log_power = 0;

  static SectorMonomialTag rational(const std::string& m_value,
                                    const std::string& b_value,
                                    std::uint32_t p) {
    const Rational m_exact(m_value);
    return {ExactScalarDescriptor::rational(m_exact.str()),
            ExactScalarDescriptor::rational((m_exact + Rational(1)).str()),
            ExactScalarDescriptor::rational(b_value), p};
  }
};

template <typename Scalar>
SectorMonomialTag sector_monomial_tag(const LocalSector<Scalar>& sector,
                                      std::uint32_t taylor_order) {
  if (sector.a.domain != ExactDomain::Rational)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
        "native monomial preparation currently requires a rational local "
        "power; it will not classify an algebraic/symbolic power by sampling");
  const Rational m = Rational(sector.a.canonical) +
                     Rational(static_cast<long>(taylor_order));
  return {ExactScalarDescriptor::rational(m.str()),
          ExactScalarDescriptor::rational((m + Rational(1)).str()), sector.b,
          sector.log_power};
}

struct MonomialIntegrationOptions {
  std::int32_t complete_max = 0;
  // Required for a branch-sensitive negative arm.  This is exact topology,
  // never inferred from the sign of an Acb imaginary midpoint.
  std::optional<std::int32_t> imaginary_sign;
};

enum class EndpointCellClass : std::uint8_t {
  DropAnalyticRegulator,
  Vanish,
  Finite,
  Divergent
};

struct EndpointLimitOptions {
  // +1/-1 records the geometric arm.  The surviving value is independent of
  // it, but validating it keeps endpoint checkpoint state explicit.
  std::int32_t approach_direction = 1;
  std::optional<std::int32_t> imaginary_sign;
  std::int32_t cancellation_digits = 24;
  bool allow_certified_numeric_cancellation = true;
};

struct EndpointLimitResult {
  std::vector<EpsilonFrame<ComplexBall>> values;
  std::size_t dropped_regulated_sectors = 0;
  std::size_t cancelled_divergent_coefficients = 0;
};

namespace integration_detail {

inline NativeIntegrationError unsupported_symbolic_regulator() {
  return NativeIntegrationError(
      NativeIntegrationErrorCode::UnsupportedSymbolicRegulator, "E10",
      "symbolic regulator slopes are not numerically sampled by the native "
      "integrator; install an exact algebraic/rational specialization first");
}

inline void validate_descriptor_facts(const ExactScalarDescriptor& value,
                                      const char* label,
                                      bool reject_symbolic_regulator = false) {
  if (reject_symbolic_regulator &&
      value.domain == ExactDomain::SymbolicRational)
    throw unsupported_symbolic_regulator();
  if (value.is_zero == TruthValue::Unknown || value.sign == ExactSign::Unknown)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
        std::string(label) +
            " lacks the exact zero/sign facts required for integration");
  if (!value.specialization.has_value())
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
        std::string(label) +
            " has no explicit numeric specialization for Acb evaluation");
}

inline void validate_tag(const SectorMonomialTag& tag) {
  validate_descriptor_facts(tag.m, "monomial power");
  validate_descriptor_facts(tag.alpha0, "monomial antiderivative power");
  validate_descriptor_facts(tag.b, "regulator slope", true);
  if (tag.log_power >
      static_cast<std::uint32_t>(std::numeric_limits<std::int32_t>::max() / 2))
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
        "logarithmic depth exceeds the native finite-frame index range");
}

inline ComplexBall ball_pow_ui(const ComplexBall& base, std::uint32_t power) {
  ComplexBall out;
  acb_pow_ui(out.raw(), base.raw(), power, ComplexBall::precision());
  return out;
}

inline ComplexBall ball_pow(const ComplexBall& base,
                            const ComplexBall& exponent) {
  ComplexBall out;
  acb_pow(out.raw(), base.raw(), exponent.raw(), ComplexBall::precision());
  return out;
}

inline ComplexBall ball_log(const ComplexBall& value) {
  ComplexBall out;
  acb_log(out.raw(), value.raw(), ComplexBall::precision());
  return out;
}

inline ComplexBall ball_exp(const ComplexBall& value) {
  ComplexBall out;
  acb_exp(out.raw(), value.raw(), ComplexBall::precision());
  return out;
}

inline ComplexBall ball_from_fmpz(const fmpz_t value) {
  ComplexBall out;
  arb_set_fmpz(acb_realref(out.raw()), value);
  arb_zero(acb_imagref(out.raw()));
  return out;
}

inline ComplexBall factorial(std::uint32_t n) {
  fmpz_t value;
  fmpz_init(value);
  fmpz_fac_ui(value, n);
  auto out = ball_from_fmpz(value);
  fmpz_clear(value);
  return out;
}

inline ComplexBall binomial(std::uint32_t n, std::uint32_t k) {
  fmpz_t value;
  fmpz_init(value);
  fmpz_bin_uiui(value, n, k);
  auto out = ball_from_fmpz(value);
  fmpz_clear(value);
  return out;
}

inline ComplexBall signed_imaginary_pi(std::int32_t sigma) {
  if (sigma != 1 && sigma != -1)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
        "negative-arm imaginary sign must be exactly +1 or -1");
  ComplexBall out;
  acb_const_pi(out.raw(), ComplexBall::precision());
  acb_mul_onei(out.raw(), out.raw());
  if (sigma < 0) acb_neg(out.raw(), out.raw());
  return out;
}

inline bool descriptor_integer_odd(const ExactScalarDescriptor& value) {
  if (value.domain != ExactDomain::Rational ||
      value.is_integer != TruthValue::Yes)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
        "the real-log PV convention requires an exact integer power");
  const auto last = value.canonical.back();
  return ((last - '0') & 1) != 0;
}

inline void validate_point(const RealEvaluationPoint& point) {
  if (point.sign < -1 || point.sign > 1)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::InvalidInterval, "E9",
        "real interval point sign must be -1, 0 or +1");
  if (!arb_is_zero(acb_imagref(point.modulus.raw())) ||
      !arb_is_nonnegative(acb_realref(point.modulus.raw())) ||
      ((point.sign == 0) != acb_is_zero(point.modulus.raw())))
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::InvalidInterval, "E9",
        "real interval point has inconsistent exact sign and modulus");
}

inline void validate_interval(const RealEvaluationPoint& lower,
                              const RealEvaluationPoint& upper) {
  validate_point(lower);
  validate_point(upper);
  if (lower.sign == 0 && upper.sign == 1) return;
  if (lower.sign == -1 && upper.sign == 0) return;
  if (lower.sign == -1 && upper.sign == 1) return;
  if (lower.sign == 1 && upper.sign == 1 &&
      arb_lt(acb_realref(lower.modulus.raw()),
             acb_realref(upper.modulus.raw())))
    return;
  if (lower.sign == -1 && upper.sign == -1 &&
      arb_lt(acb_realref(upper.modulus.raw()),
             acb_realref(lower.modulus.raw())))
    return;
  throw NativeIntegrationError(
      NativeIntegrationErrorCode::InvalidInterval, "E9",
      "interval order is not certified by the exact signs and Acb moduli");
}

inline std::int32_t require_sigma(
    const MonomialIntegrationOptions& options) {
  if (!options.imaginary_sign.has_value()) {
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
        "branch-sensitive negative-arm integration requires an explicit "
        "imaginary sign");
  }
  const auto sigma = *options.imaginary_sign;
  if (sigma != 1 && sigma != -1)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
        "negative-arm imaginary sign must be exactly +1 or -1");
  return sigma;
}

inline EpsilonFrame<ComplexBall> log_difference_frame(
    std::uint32_t p, const ComplexBall& from_log, const ComplexBall& to_log,
    std::int32_t complete_max) {
  if (complete_max < static_cast<std::int32_t>(p))
    throw std::invalid_argument(
        "requested complete order lies below the log monomial minimum");
  std::vector<ComplexBall> coefficients(
      static_cast<std::size_t>(complete_max - static_cast<std::int32_t>(p) + 1),
      ComplexBall(0));
  coefficients.front() =
      (ball_pow_ui(to_log, p + 1) - ball_pow_ui(from_log, p + 1)) /
      (factorial(p) * ComplexBall(static_cast<long>(p + 1)));
  return EpsilonFrame<ComplexBall>(static_cast<std::int32_t>(p),
                                   std::move(coefficients));
}

// Combined m=-1 primitive difference.  The two 1/(b eps) endpoint rows are
// never materialized: this is regular and starts at eps^p.
inline EpsilonFrame<ComplexBall> paired_pole_frame(
    const SectorMonomialTag& tag, const ComplexBall& from_log,
    const ComplexBall& to_log, std::int32_t complete_max) {
  const auto p = tag.log_power;
  if (complete_max < static_cast<std::int32_t>(p))
    throw std::invalid_argument(
        "requested complete order lies below the paired-pole minimum");
  const auto& b = tag.b.numeric();
  const auto p_factorial = factorial(p);
  std::vector<ComplexBall> coefficients;
  coefficients.reserve(static_cast<std::size_t>(
      complete_max - static_cast<std::int32_t>(p) + 1));
  for (std::int64_t k = static_cast<std::int32_t>(p); k <= complete_max;
       ++k) {
    const auto q = static_cast<std::uint32_t>(k - p);
    const auto degree = p + q + 1;
    auto value = ball_pow_ui(b, q) *
                 (ball_pow_ui(to_log, degree) -
                  ball_pow_ui(from_log, degree));
    value = value /
            (p_factorial * factorial(q) *
             ComplexBall(static_cast<long>(degree)));
    coefficients.push_back(std::move(value));
  }
  return EpsilonFrame<ComplexBall>(static_cast<std::int32_t>(p),
                                   std::move(coefficients));
}

// F(t) at one nonzero boundary.  log_t already contains the prescribed
// negative-arm i*pi*sigma and phase0 is exp(i*pi*sigma*(m+1)).
inline EpsilonFrame<ComplexBall> antiderivative_at(
    const SectorMonomialTag& tag, const ComplexBall& modulus,
    const ComplexBall& log_t, const ComplexBall& phase0,
    std::int32_t complete_max) {
  const auto p = tag.log_power;
  const bool alpha_zero = tag.alpha0.is_zero == TruthValue::Yes;
  const bool b_zero = tag.b.is_zero == TruthValue::Yes;
  if (alpha_zero && b_zero)
    return log_difference_frame(p, ComplexBall(0), log_t, complete_max);

  const auto prefactor =
      phase0 * ball_pow(modulus, tag.alpha0.numeric());
  const auto& b = tag.b.numeric();

  if (alpha_zero) {
    if (complete_max < -1)
      throw std::invalid_argument(
          "requested complete order lies below the pole-cell minimum");
    std::vector<ComplexBall> coefficients;
    coefficients.reserve(static_cast<std::size_t>(complete_max + 2));
    for (std::int64_t k = -1; k <= complete_max; ++k) {
      ComplexBall value(0);
      for (std::uint32_t j = 0; j <= p; ++j) {
        const auto q_signed = static_cast<std::int64_t>(k) - p + j + 1;
        if (q_signed < 0) continue;
        const auto q = static_cast<std::uint32_t>(q_signed);
        auto term = ball_pow_ui(log_t, p - j) /
                    factorial(p - j);
        term *= ball_pow_ui(b, q) * ball_pow_ui(log_t, q) /
                factorial(q);
        term = term / ball_pow_ui(b, j + 1);
        if ((j & 1U) != 0) term = -term;
        value += term;
      }
      coefficients.push_back(prefactor * value);
    }
    return EpsilonFrame<ComplexBall>(-1, std::move(coefficients));
  }

  if (complete_max < static_cast<std::int32_t>(p))
    throw std::invalid_argument(
        "requested complete order lies below the regular monomial minimum");
  const auto& alpha0 = tag.alpha0.numeric();
  std::vector<ComplexBall> coefficients;
  coefficients.reserve(static_cast<std::size_t>(
      complete_max - static_cast<std::int32_t>(p) + 1));
  for (std::int64_t k = static_cast<std::int32_t>(p); k <= complete_max;
       ++k) {
    const auto r = static_cast<std::uint32_t>(k - p);
    if (b_zero && r > 0) {
      coefficients.emplace_back(0);
      continue;
    }
    ComplexBall outer_sum(0);
    for (std::uint32_t j = 0; j <= p; ++j) {
      ComplexBall inner_sum(0);
      for (std::uint32_t s = 0; s <= r; ++s) {
        auto inner = ball_pow_ui(log_t, r - s) /
                     factorial(r - s);
        inner *= binomial(j + s, s);
        inner = inner / ball_pow_ui(alpha0, j + 1 + s);
        if ((s & 1U) != 0) inner = -inner;
        inner_sum += inner;
      }
      auto term = ball_pow_ui(log_t, p - j) /
                  factorial(p - j);
      term *= inner_sum;
      if ((j & 1U) != 0) term = -term;
      outer_sum += term;
    }
    coefficients.push_back(prefactor * ball_pow_ui(b, r) * outer_sum);
  }
  return EpsilonFrame<ComplexBall>(static_cast<std::int32_t>(p),
                                   std::move(coefficients));
}

inline ComplexBall endpoint_log(const RealEvaluationPoint& point,
                                std::optional<std::int32_t> sigma) {
  auto result = ball_log(point.modulus);
  if (point.sign < 0) result += signed_imaginary_pi(*sigma);
  return result;
}

inline ComplexBall endpoint_phase(const SectorMonomialTag& tag,
                                  const RealEvaluationPoint& point,
                                  std::optional<std::int32_t> sigma) {
  if (point.sign >= 0) return ComplexBall(1);
  return ball_exp(signed_imaginary_pi(*sigma) * tag.alpha0.numeric());
}

inline bool branch_sensitive_on_negative_arm(const SectorMonomialTag& tag) {
  return tag.b.is_zero != TruthValue::Yes ||
         tag.m.is_integer != TruthValue::Yes || tag.log_power > 0;
}

inline EpsilonFrame<ComplexBall> endpoint_primitive(
    const SectorMonomialTag& tag, const RealEvaluationPoint& outer,
    const MonomialIntegrationOptions& options) {
  if (tag.b.is_zero == TruthValue::Yes &&
      tag.alpha0.sign != ExactSign::Positive)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::DivergentEndpoint, "E2",
        "a b=0 endpoint monomial with m+1<=0 reached the primitive before "
        "the assembled cancellation gate");

  std::optional<std::int32_t> sigma;
  if (outer.sign < 0 && branch_sensitive_on_negative_arm(tag))
    sigma = require_sigma(options);
  else if (outer.sign < 0)
    sigma = options.imaginary_sign.value_or(1);
  const auto log_t = endpoint_log(outer, sigma);
  const auto phase = endpoint_phase(tag, outer, sigma);
  auto result = antiderivative_at(tag, outer.modulus, log_t, phase,
                                  options.complete_max);
  return outer.sign > 0 ? result : -result;
}

template <typename Scalar>
bool material_sector(const LocalSector<Scalar>& sector) {
  return std::any_of(sector.coefficients.begin(), sector.coefficients.end(),
                     [](const Scalar& value) {
                       return !ScalarTraits<Scalar>::is_zero(value);
                     });
}

inline bool numeric_cancellation_certified(
    const ComplexBall& total, const std::vector<ComplexBall>& terms,
    const EndpointLimitOptions& options) {
  if (total.is_zero()) return true;
  if (!options.allow_certified_numeric_cancellation ||
      options.cancellation_digits < 0)
    return false;
  auto scale = Magnitude::one();
  for (const auto& term : terms)
    scale = Magnitude::maximum(scale, Magnitude::upper_abs(term));
  const auto tolerance =
      Magnitude::decimal("1e-" + std::to_string(options.cancellation_digits));
  return Magnitude::upper_abs(total) <= scale * tolerance;
}

template <typename Scalar>
ComplexBall coefficient_ball(const Scalar& value) {
  return local_detail::to_ball(value);
}

template <typename Scalar>
bool divergent_sum_passes(const std::vector<const Scalar*>& term_ptrs,
                          const EndpointLimitOptions& options,
                          bool* numeric_cancellation) {
  *numeric_cancellation = false;
  if constexpr (std::is_same_v<Scalar, Rational>) {
    Rational total(0);
    for (const auto* term : term_ptrs) total += *term;
    return total.is_zero();
  } else if constexpr (std::is_same_v<Scalar, ComplexBall>) {
    ComplexBall total(0);
    std::vector<ComplexBall> terms;
    terms.reserve(term_ptrs.size());
    for (const auto* term : term_ptrs) {
      total += *term;
      terms.push_back(*term);
    }
    const auto pass = numeric_cancellation_certified(total, terms, options);
    *numeric_cancellation = pass && !total.is_zero();
    return pass;
  } else {
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
        "endpoint limit evaluation does not sample symbolic coefficient "
        "fields; specialize the coefficient frame explicitly");
  }
}

}  // namespace integration_detail

inline EndpointCellClass classify_endpoint_limit_cell(
    const SectorMonomialTag& tag) {
  integration_detail::validate_tag(tag);
  if (tag.b.is_zero == TruthValue::No)
    return EndpointCellClass::DropAnalyticRegulator;
  if (tag.b.is_zero != TruthValue::Yes)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
        "endpoint classification requires an exact b==0 decision");
  if (tag.m.sign == ExactSign::Positive) return EndpointCellClass::Vanish;
  if (tag.m.sign == ExactSign::Zero && tag.log_power == 0)
    return EndpointCellClass::Finite;
  if (tag.m.sign == ExactSign::Negative ||
      (tag.m.sign == ExactSign::Zero && tag.log_power > 0))
    return EndpointCellClass::Divergent;
  throw NativeIntegrationError(
      NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
      "endpoint monomial power could not be classified exactly");
}

// Closed-form definite integral of one normalized sector monomial
//   t^(m+b eps) (eps Log t)^p / p!
// over a real interval.  Center endpoints use the dimensional-regulator
// boundary F(0):=0 for exact b!=0.  A crossing is always one paired owner.
inline EpsilonFrame<ComplexBall> integrate_sector_monomial(
    const SectorMonomialTag& tag, const RealEvaluationPoint& lower,
    const RealEvaluationPoint& upper,
    const MonomialIntegrationOptions& options) {
  using namespace integration_detail;
  validate_tag(tag);
  validate_interval(lower, upper);

  if (lower.sign == 0 || upper.sign == 0) {
    const auto& outer = lower.sign == 0 ? upper : lower;
    return endpoint_primitive(tag, outer, options);
  }

  const bool crossing = lower.sign < 0 && upper.sign > 0;
  if (crossing && tag.b.is_zero == TruthValue::Yes &&
      tag.m.is_integer == TruthValue::Yes) {
    const auto log_a = ball_log(lower.modulus);
    const auto log_b = ball_log(upper.modulus);
    if (tag.alpha0.is_zero == TruthValue::Yes)
      return log_difference_frame(tag.log_power, log_a, log_b,
                                  options.complete_max);

    auto positive = antiderivative_at(
        tag, upper.modulus, log_b, ComplexBall(1), options.complete_max);
    auto negative_modulus = antiderivative_at(
        tag, lower.modulus, log_a, ComplexBall(1), options.complete_max);
    // Real-log PV/Hadamard finite part.  The negative-arm substitution gives
    // (-1)^m; the formal center terms are dropped jointly, once.
    if (descriptor_integer_odd(tag.m)) negative_modulus = -negative_modulus;
    return positive + negative_modulus;
  }

  std::optional<std::int32_t> sigma;
  const bool has_negative_arm = lower.sign < 0 || upper.sign < 0;
  if (has_negative_arm && branch_sensitive_on_negative_arm(tag))
    sigma = require_sigma(options);
  else if (has_negative_arm)
    sigma = options.imaginary_sign.value_or(1);

  const auto lower_log = endpoint_log(lower, sigma);
  const auto upper_log = endpoint_log(upper, sigma);
  if (tag.alpha0.is_zero == TruthValue::Yes) {
    if (tag.b.is_zero == TruthValue::Yes)
      return log_difference_frame(tag.log_power, lower_log, upper_log,
                                  options.complete_max);
    // m=-1,b!=0: pair before frame arithmetic, eliminating the spurious
    // eps^-1 row on same-side and crossing intervals.
    return paired_pole_frame(tag, lower_log, upper_log,
                             options.complete_max);
  }

  auto from = antiderivative_at(
      tag, lower.modulus, lower_log, endpoint_phase(tag, lower, sigma),
      options.complete_max);
  auto to = antiderivative_at(
      tag, upper.modulus, upper_log, endpoint_phase(tag, upper, sigma),
      options.complete_max);
  return to - from;
}

template <typename Scalar>
EndpointLimitResult endpoint_sector_limit(
    const LocalSolution<Scalar>& solution,
    const EndpointLimitOptions& options = {}) {
  using namespace integration_detail;
  validate_local_solution(solution, false);
  if constexpr (std::is_same_v<Scalar, SymbolicRational>)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
        "endpoint evaluation requires an explicitly specialized coefficient "
        "frame; symbolic coefficients are not sampled");
  if (options.approach_direction != 1 && options.approach_direction != -1)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::InvalidInterval, "E9",
        "endpoint approach direction must be exactly +1 or -1");
  if (options.imaginary_sign.has_value() &&
      *options.imaginary_sign != 1 && *options.imaginary_sign != -1)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
        "endpoint imaginary sign must be exactly +1 or -1");

  bool branch_sensitive = false;
  for (const auto& sector : solution.sectors) {
    if (sector.b.domain == ExactDomain::SymbolicRational)
      throw unsupported_symbolic_regulator();
    if (sector.a.domain != ExactDomain::Rational)
      throw NativeIntegrationError(
          NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
          "endpoint classification currently requires rational local powers");
    if (!material_sector(sector)) continue;
    if (sector.b.is_zero != TruthValue::Yes || sector.log_power > 0 ||
        sector.a.is_integer != TruthValue::Yes)
      branch_sensitive = true;
  }
  const auto derived_sigma = derive_chart_imaginary_sign(solution);
  if (branch_sensitive && !derived_sigma.has_value() &&
      !options.imaginary_sign.has_value())
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
        "multivalued endpoint data has no explicit or chart-derived branch "
        "prescription");
  if (derived_sigma.has_value() && options.imaginary_sign.has_value() &&
      *derived_sigma != *options.imaginary_sign)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
        "explicit endpoint branch sign conflicts with the prepared chart");

  struct CellRef {
    const LocalSector<Scalar>* sector = nullptr;
    std::uint32_t taylor = 0;
    SectorMonomialTag tag;
  };
  std::vector<CellRef> finite;
  std::map<std::pair<std::string, std::uint32_t>, std::vector<CellRef>>
      divergent;
  EndpointLimitResult result;
  const auto taylor_width = solution.taylor_width();
  for (const auto& sector : solution.sectors) {
    if (sector.b.is_zero == TruthValue::No) {
      if (material_sector(sector)) ++result.dropped_regulated_sectors;
      continue;  // exact dimensional-regulator endpoint drop
    }
    if (sector.b.is_zero != TruthValue::Yes)
      throw NativeIntegrationError(
          NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
          "endpoint classification requires an exact regulator zero fact");

    const Rational a(sector.a.canonical);
    const Rational first_unseen =
        a + Rational(static_cast<long>(taylor_width));
    const auto unseen = ExactScalarDescriptor::rational(first_unseen.str());
    if (unseen.sign != ExactSign::Positive)
      throw NativeIntegrationError(
          NativeIntegrationErrorCode::IncompleteTaylorWindow, "E10",
          "endpoint limit is uncertified: the first unseen Taylor cell can "
          "still have nonpositive absolute power");

    for (std::uint32_t n = 0; n < taylor_width; ++n) {
      auto tag = sector_monomial_tag(sector, n);
      const auto classification = classify_endpoint_limit_cell(tag);
      CellRef cell{&sector, n, std::move(tag)};
      if (classification == EndpointCellClass::Finite)
        finite.push_back(std::move(cell));
      else if (classification == EndpointCellClass::Divergent) {
        const auto key =
            std::make_pair(cell.tag.m.canonical, cell.tag.log_power);
        divergent[key].push_back(std::move(cell));
      }
    }
  }

  const auto eps_width = solution.epsilon.width();
  for (const auto& [key, cells] : divergent) {
    for (std::size_t ei = 0; ei < eps_width; ++ei) {
      const auto epsilon_power = local_detail::checked_i32(
          static_cast<std::int64_t>(solution.epsilon.min_power) + ei,
          "endpoint epsilon power");
      for (std::uint32_t component = 0; component < solution.dimension;
           ++component) {
        std::vector<const Scalar*> terms;
        terms.reserve(cells.size());
        for (const auto& cell : cells)
          terms.push_back(&cell.sector->coefficients[local_detail::sector_index(
              solution, ei, cell.taylor, component)]);
        bool numeric_cancellation = false;
        if (!divergent_sum_passes(terms, options, &numeric_cancellation)) {
          NativeIntegrationError error(
              NativeIntegrationErrorCode::DivergentEndpoint, "E2",
              "divergent b=0 endpoint content does not cancel in its exact "
              "absolute-power/log-depth class");
          error.absolute_power = key.first;
          error.log_power = key.second;
          error.epsilon_power = epsilon_power;
          error.component = component;
          throw error;
        }
        if (numeric_cancellation)
          ++result.cancelled_divergent_coefficients;
      }
    }
  }

  if (finite.empty()) {
    result.values.reserve(solution.dimension);
    for (std::uint32_t component = 0; component < solution.dimension;
         ++component)
      result.values.push_back(
          EpsilonFrame<ComplexBall>::zero(solution.epsilon.complete_max));
    return result;
  }

  result.values.reserve(solution.dimension);
  for (std::uint32_t component = 0; component < solution.dimension;
       ++component) {
    std::vector<ComplexBall> coefficients(eps_width, ComplexBall(0));
    for (std::size_t ei = 0; ei < eps_width; ++ei)
      for (const auto& cell : finite)
        coefficients[ei] += coefficient_ball(
            cell.sector->coefficients[local_detail::sector_index(
                solution, ei, cell.taylor, component)]);
    result.values.emplace_back(solution.epsilon, std::move(coefficients));
  }
  return result;
}

}  // namespace diffexp2
