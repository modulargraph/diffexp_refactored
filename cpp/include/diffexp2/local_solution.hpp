#pragma once

#include "diffexp2/scalar.hpp"
#include "diffexp2/series_types.hpp"

#include <flint/mag.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace diffexp2 {

// Native representation of the SectorSeries LocalSolution contract.  The
// exact descriptor and its numeric specialization are deliberately separate:
// branch, window and endpoint decisions are made from exact facts, never from
// an Acb enclosure.  A symbolic/algebraic descriptor may therefore survive a
// checkpoint even when this process has not installed a numeric specialization
// for it yet.

enum class TruthValue : std::uint8_t { No, Yes, Unknown };
enum class ExactSign : std::int8_t { Negative = -1, Zero = 0, Positive = 1,
                                     Unknown = 2 };
enum class ExactDomain : std::uint8_t { Rational, SymbolicRational, Algebraic };

struct ExactScalarDescriptor {
  ExactDomain domain = ExactDomain::Rational;
  std::string canonical;
  std::vector<std::string> symbols;
  TruthValue is_zero = TruthValue::Unknown;
  TruthValue is_integer = TruthValue::Unknown;
  ExactSign sign = ExactSign::Unknown;
  std::optional<ComplexBall> specialization;

  static ExactScalarDescriptor rational(const std::string& value) {
    const Rational q(value);
    const auto canonical = q.str();
    ExactScalarDescriptor out;
    out.domain = ExactDomain::Rational;
    out.canonical = canonical;
    out.is_zero = q.is_zero() ? TruthValue::Yes : TruthValue::No;
    out.is_integer = canonical.find('/') == std::string::npos
                         ? TruthValue::Yes
                         : TruthValue::No;
    out.sign = q.is_zero() ? ExactSign::Zero
                           : (canonical.front() == '-' ? ExactSign::Negative
                                                       : ExactSign::Positive);
    out.specialization = ComplexBall::from_strings(canonical);
    return out;
  }

  static ExactScalarDescriptor symbolic(
      std::string expression, std::vector<std::string> variables,
      TruthValue zero, TruthValue integer, ExactSign exact_sign,
      std::optional<ComplexBall> numeric_specialization = std::nullopt) {
    ExactScalarDescriptor out;
    out.domain = ExactDomain::SymbolicRational;
    out.canonical = std::move(expression);
    out.symbols = std::move(variables);
    out.is_zero = zero;
    out.is_integer = integer;
    out.sign = exact_sign;
    out.specialization = std::move(numeric_specialization);
    return out;
  }

  static ExactScalarDescriptor algebraic(
      std::string expression, TruthValue zero, TruthValue integer,
      ExactSign exact_sign, ComplexBall numeric_specialization) {
    ExactScalarDescriptor out;
    out.domain = ExactDomain::Algebraic;
    out.canonical = std::move(expression);
    out.is_zero = zero;
    out.is_integer = integer;
    out.sign = exact_sign;
    out.specialization = std::move(numeric_specialization);
    return out;
  }

  [[nodiscard]] const ComplexBall& numeric() const {
    if (!specialization.has_value()) {
      throw std::domain_error(
          "exact tag has no numeric specialization: " + canonical);
    }
    return *specialization;
  }
};

class Magnitude {
 public:
  Magnitude() { mag_init(value_); }
  Magnitude(const Magnitude& other) : Magnitude() { mag_set(value_, other.value_); }
  Magnitude(Magnitude&& other) noexcept : Magnitude() { mag_swap(value_, other.value_); }
  ~Magnitude() { mag_clear(value_); }

  Magnitude& operator=(Magnitude other) noexcept {
    mag_swap(value_, other.value_);
    return *this;
  }

  static Magnitude zero() {
    Magnitude out;
    mag_zero(out.value_);
    return out;
  }
  static Magnitude one() {
    Magnitude out;
    mag_one(out.value_);
    return out;
  }
  static Magnitude from_ui(ulong value) {
    Magnitude out;
    mag_set_ui(out.value_, value);
    return out;
  }
  static Magnitude decimal(const std::string& value) {
    arb_t ball;
    arb_init(ball);
    if (arb_set_str(ball, value.c_str(), ComplexBall::precision()) != 0 ||
        !arb_is_nonnegative(ball)) {
      arb_clear(ball);
      throw std::invalid_argument("invalid nonnegative magnitude: " + value);
    }
    Magnitude out;
    arb_get_mag(out.value_, ball);
    arb_clear(ball);
    return out;
  }
  static Magnitude upper_abs(const ComplexBall& value) {
    Magnitude out;
    acb_get_mag(out.value_, value.raw());
    return out;
  }
  static Magnitude lower_abs(const ComplexBall& value) {
    Magnitude out;
    acb_get_mag_lower(out.value_, value.raw());
    return out;
  }

  [[nodiscard]] bool is_zero() const { return mag_is_zero(value_); }
  [[nodiscard]] bool is_finite() const { return !mag_is_inf(value_); }
  [[nodiscard]] double approximate_upper() const { return mag_get_d(value_); }

  // FLINT's dump/load representation preserves the exact mag exponent and
  // mantissa.  Checkpoints must use this instead of a diagnostic double so a
  // restored certificate is neither weakened nor accidentally strengthened.
  [[nodiscard]] std::string dump_exact() const {
    char* raw = mag_dump_str(value_);
    if (raw == nullptr) throw std::bad_alloc();
    std::string output(raw);
    flint_free(raw);
    return output;
  }

  static Magnitude from_exact_dump(const std::string& dump) {
    if (dump.empty())
      throw std::invalid_argument("empty exact magnitude checkpoint dump");
    Magnitude output;
    if (mag_load_str(output.value_, dump.c_str()) != 0)
      throw std::invalid_argument("invalid exact magnitude checkpoint dump");
    return output;
  }

  friend Magnitude operator+(const Magnitude& a, const Magnitude& b) {
    Magnitude out;
    mag_add(out.value_, a.value_, b.value_);
    return out;
  }
  friend Magnitude operator*(const Magnitude& a, const Magnitude& b) {
    Magnitude out;
    mag_mul(out.value_, a.value_, b.value_);
    return out;
  }
  friend Magnitude operator/(const Magnitude& a, const Magnitude& b) {
    if (b.is_zero()) throw std::domain_error("magnitude division by zero");
    Magnitude out;
    mag_div(out.value_, a.value_, b.value_);
    return out;
  }
  // The following helpers deliberately name the direction of their
  // enclosure.  In particular, positive_difference_lower is used when a
  // triangle inequality proves a denominator is separated from zero; feeding
  // that lower bound to reciprocal_upper then gives a safe upper bound.
  static Magnitude positive_difference_lower(const Magnitude& a,
                                             const Magnitude& b) {
    Magnitude out;
    mag_sub_lower(out.value_, a.value_, b.value_);
    return out;
  }
  [[nodiscard]] Magnitude reciprocal_upper() const {
    if (is_zero())
      throw std::domain_error("magnitude reciprocal of zero lower bound");
    Magnitude out;
    mag_inv(out.value_, value_);
    return out;
  }
  [[nodiscard]] Magnitude exponential_upper() const {
    Magnitude out;
    mag_exp(out.value_, value_);
    return out;
  }
  [[nodiscard]] Magnitude power_upper(ulong exponent) const {
    Magnitude out;
    mag_pow_ui(out.value_, value_, exponent);
    return out;
  }
  Magnitude& operator+=(const Magnitude& b) {
    mag_add(value_, value_, b.value_);
    return *this;
  }
  friend bool operator<=(const Magnitude& a, const Magnitude& b) {
    return mag_cmp(a.value_, b.value_) <= 0;
  }
  friend bool operator>(const Magnitude& a, const Magnitude& b) {
    return mag_cmp(a.value_, b.value_) > 0;
  }

  static Magnitude maximum(const Magnitude& a, const Magnitude& b) {
    return a > b ? a : b;
  }

 private:
  mag_t value_;
};

enum class ErrorGuarantee : std::uint8_t { None, Advisory, Certified };

struct ErrorEnvelope {
  EpsilonWindow frame;
  ErrorGuarantee guarantee = ErrorGuarantee::None;
  std::vector<Magnitude> absolute;
  std::string provenance;

  [[nodiscard]] bool empty() const { return absolute.empty(); }
};

struct Prescription {
  std::string factor_exact;
  std::int32_t sign = 0;
  std::uint32_t multiplicity = 0;
  std::int32_t leading_coefficient_sign = 0;
};

struct ChartGeometry {
  std::string center_exact;
  std::string scale_exact = "1";
  ComplexBall radius = ComplexBall(1);
  bool infinite_radius = false;
};

template <typename Scalar>
struct LocalSector {
  ExactScalarDescriptor a;
  ExactScalarDescriptor b;
  std::uint32_t log_power = 0;
  // Flat [epsilon][Taylor][component], component fastest.
  std::vector<Scalar> coefficients;
};

template <typename Scalar>
struct LocalSolution {
  ChartGeometry chart;
  EpsilonWindow epsilon;
  std::uint32_t taylor_complete_max = 0;
  std::uint32_t dimension = 0;
  std::vector<LocalSector<Scalar>> sectors;
  std::vector<Prescription> prescriptions;
  ErrorEnvelope error;
  std::string checkpoint_identity;

  [[nodiscard]] std::size_t taylor_width() const {
    return static_cast<std::size_t>(taylor_complete_max) + 1;
  }
  [[nodiscard]] std::size_t sector_size() const {
    return epsilon.width() * taylor_width() * dimension;
  }
};

struct RealEvaluationPoint {
  // Nonnegative |t|.  sign is exact topology and must not be inferred from
  // the ball: -1, 0 or +1.
  ComplexBall modulus;
  std::int32_t sign = 0;
  std::string exact_coordinate;

  static RealEvaluationPoint rational(const std::string& value) {
    const Rational q(value);
    const auto canonical = q.str();
    RealEvaluationPoint out{ComplexBall::from_strings(
                                canonical.front() == '-' ? canonical.substr(1)
                                                          : canonical),
                            q.is_zero() ? 0 : (canonical.front() == '-' ? -1 : 1),
                            canonical};
    return out;
  }
};

struct EpsilonVector {
  EpsilonWindow epsilon;
  std::uint32_t dimension = 0;
  // Flat [epsilon][component].
  std::vector<ComplexBall> coefficients;
  ErrorEnvelope error;

  [[nodiscard]] std::size_t index(std::int32_t power,
                                  std::uint32_t component) const {
    return (static_cast<std::size_t>(power - epsilon.min_power) * dimension) +
           component;
  }
  [[nodiscard]] const ComplexBall& at(std::int32_t power,
                                      std::uint32_t component) const {
    if (power < epsilon.min_power || power > epsilon.complete_max ||
        component >= dimension)
      throw std::out_of_range("epsilon vector coefficient outside window");
    return coefficients.at(index(power, component));
  }
  ComplexBall& at(std::int32_t power, std::uint32_t component) {
    if (power < epsilon.min_power || power > epsilon.complete_max ||
        component >= dimension)
      throw std::out_of_range("epsilon vector coefficient outside window");
    return coefficients.at(index(power, component));
  }
};

struct EpsilonMatrix {
  EpsilonWindow epsilon;
  std::uint32_t dimension = 0;
  // Flat [epsilon][row][column].
  std::vector<ComplexBall> coefficients;

  [[nodiscard]] std::size_t index(std::int32_t power, std::uint32_t row,
                                  std::uint32_t column) const {
    return ((static_cast<std::size_t>(power - epsilon.min_power) * dimension +
             row) * dimension) + column;
  }
  [[nodiscard]] const ComplexBall& at(std::int32_t power, std::uint32_t row,
                                      std::uint32_t column) const {
    if (power < epsilon.min_power || power > epsilon.complete_max ||
        row >= dimension || column >= dimension)
      throw std::out_of_range("epsilon matrix coefficient outside window");
    return coefficients.at(index(power, row, column));
  }
};

struct EvaluationOptions {
  std::optional<std::int32_t> imaginary_sign;
  std::uint32_t t_order_reduction = 0;
  bool compute_tail_estimate = true;
};

struct LocalEvaluation {
  EpsilonVector value;
  EpsilonVector theta_value;
  std::optional<std::int32_t> imaginary_sign;
  bool arithmetic_enclosed = true;
};

namespace local_detail {

inline std::int32_t checked_i32(std::int64_t value, const char* label) {
  if (value < std::numeric_limits<std::int32_t>::min() ||
      value > std::numeric_limits<std::int32_t>::max())
    throw std::overflow_error(std::string(label) + " exceeds int32 range");
  return static_cast<std::int32_t>(value);
}

inline bool exactly_real(const ComplexBall& value) {
  return arb_is_zero(acb_imagref(value.raw()));
}

inline ComplexBall cb_log(const ComplexBall& value) {
  ComplexBall out;
  acb_log(out.raw(), value.raw(), ComplexBall::precision());
  return out;
}

inline ComplexBall cb_exp(const ComplexBall& value) {
  ComplexBall out;
  acb_exp(out.raw(), value.raw(), ComplexBall::precision());
  return out;
}

inline ComplexBall cb_pow(const ComplexBall& base, const ComplexBall& exponent) {
  ComplexBall out;
  acb_pow(out.raw(), base.raw(), exponent.raw(), ComplexBall::precision());
  return out;
}

inline ComplexBall cb_pow_ui(const ComplexBall& base, std::uint32_t exponent) {
  ComplexBall out;
  acb_pow_ui(out.raw(), base.raw(), exponent, ComplexBall::precision());
  return out;
}

inline ComplexBall cb_div_ui(const ComplexBall& value, std::uint32_t divisor) {
  if (divisor == 0) throw std::domain_error("integer division by zero");
  ComplexBall out;
  acb_div_ui(out.raw(), value.raw(), divisor, ComplexBall::precision());
  return out;
}

inline ComplexBall imaginary_pi(std::int32_t sign) {
  ComplexBall out;
  acb_const_pi(out.raw(), ComplexBall::precision());
  acb_mul_onei(out.raw(), out.raw());
  if (sign < 0) acb_neg(out.raw(), out.raw());
  return out;
}

inline bool exact_zero(const ComplexBall& value) { return acb_is_zero(value.raw()); }

template <typename Scalar>
ComplexBall to_ball(const Scalar& value);

template <>
inline ComplexBall to_ball<Rational>(const Rational& value) {
  return ComplexBall::from_strings(value.str());
}

template <>
inline ComplexBall to_ball<ComplexBall>(const ComplexBall& value) {
  return value;
}

template <>
inline ComplexBall to_ball<SymbolicRational>(const SymbolicRational&) {
  throw std::domain_error(
      "symbolic coefficients require an explicit regulator specialization "
      "before native local evaluation");
}

template <typename Scalar>
void validate_local_solution(const LocalSolution<Scalar>& solution,
                             bool require_numeric_tags) {
  if (solution.dimension == 0) throw std::invalid_argument("zero local dimension");
  const auto expected = solution.sector_size();
  if (solution.sectors.empty()) throw std::invalid_argument("no local sectors");
  if (!solution.chart.infinite_radius) {
    if (!exactly_real(solution.chart.radius) ||
        !arb_is_positive(acb_realref(solution.chart.radius.raw())))
      throw std::invalid_argument("chart radius must be a provably positive real ball");
  }
  for (const auto& sector : solution.sectors) {
    if (sector.a.canonical.empty() || sector.b.canonical.empty())
      throw std::invalid_argument("sector tags must retain exact canonical forms");
    if (sector.coefficients.size() != expected)
      throw std::invalid_argument("sector coefficient tensor has wrong dimensions");
    if (require_numeric_tags &&
        (!sector.a.specialization.has_value() ||
         !sector.b.specialization.has_value()))
      throw std::domain_error("numeric local evaluation requires specialized tags");
    if (sector.a.specialization.has_value() &&
        !exactly_real(*sector.a.specialization))
      throw std::invalid_argument("sector a specialization must be real");
    if (sector.b.specialization.has_value() &&
        !exactly_real(*sector.b.specialization))
      throw std::invalid_argument("sector b specialization must be real");
  }
  for (const auto& prescription : solution.prescriptions) {
    if ((prescription.sign != 1 && prescription.sign != -1) ||
        (prescription.leading_coefficient_sign != 1 &&
         prescription.leading_coefficient_sign != -1) ||
        prescription.multiplicity == 0)
      throw std::invalid_argument("malformed analytic-continuation prescription");
  }
  if (!solution.error.empty() &&
      solution.error.absolute.size() != solution.error.frame.width())
    throw std::invalid_argument("local error envelope has wrong epsilon width");
}

template <typename Scalar>
bool material_sector(const LocalSector<Scalar>& sector) {
  return std::any_of(sector.coefficients.begin(), sector.coefficients.end(),
                     [](const Scalar& value) {
                       return !ScalarTraits<Scalar>::is_zero(value);
                     });
}

template <typename Scalar>
std::optional<std::int32_t> chart_imaginary_sign(
    const LocalSolution<Scalar>& solution) {
  std::optional<std::int32_t> derived;
  for (const auto& prescription : solution.prescriptions) {
    if ((prescription.multiplicity & 1U) == 0) continue;
    const auto candidate = prescription.sign *
                           prescription.leading_coefficient_sign;
    if (derived.has_value() && *derived != candidate)
      throw std::domain_error(
          "conflicting odd-multiplicity analytic-continuation prescriptions");
    derived = candidate;
  }
  return derived;
}

template <typename Scalar>
bool branch_sensitive(const LocalSolution<Scalar>& solution) {
  return std::any_of(solution.sectors.begin(), solution.sectors.end(),
      [](const LocalSector<Scalar>& sector) {
        if (!material_sector(sector)) return false;
        return sector.a.is_integer != TruthValue::Yes ||
               sector.b.is_zero != TruthValue::Yes || sector.log_power > 0;
      });
}

template <typename Scalar>
std::size_t sector_index(const LocalSolution<Scalar>& solution,
                         std::size_t epsilon_index, std::size_t taylor_index,
                         std::uint32_t component) {
  return ((epsilon_index * solution.taylor_width() + taylor_index) *
          solution.dimension) + component;
}

template <typename Scalar>
EpsilonVector evaluate_value(const LocalSolution<Scalar>& solution,
                             const RealEvaluationPoint& point,
                             const EvaluationOptions& options,
                             std::optional<std::int32_t>* used_sign) {
  local_detail::validate_local_solution(solution, true);
  if (point.sign < -1 || point.sign > 1)
    throw std::invalid_argument("evaluation point sign must be -1, 0 or +1");
  if (!exactly_real(point.modulus) ||
      !arb_is_nonnegative(acb_realref(point.modulus.raw())))
    throw std::invalid_argument("evaluation point modulus must be nonnegative real");
  if ((point.sign == 0) != acb_is_zero(point.modulus.raw()))
    throw std::invalid_argument("evaluation point sign/modulus mismatch");
  if (!solution.chart.infinite_radius &&
      !arb_lt(acb_realref(point.modulus.raw()),
              acb_realref(solution.chart.radius.raw())))
    throw std::domain_error("evaluation point is not provably inside chart radius");

  const bool needs_branch = branch_sensitive(solution);
  auto sigma = options.imaginary_sign.has_value()
                   ? options.imaginary_sign
                   : chart_imaginary_sign(solution);
  if (sigma.has_value() && *sigma != 1 && *sigma != -1)
    throw std::invalid_argument("imaginary sign must be +1 or -1");
  if (point.sign < 0 && needs_branch && !sigma.has_value())
    throw std::domain_error(
        "negative multivalued local evaluation has no derivable imaginary sign");
  if (used_sign != nullptr) *used_sign = point.sign < 0 ? sigma : std::nullopt;

  if (point.sign == 0) {
    for (const auto& sector : solution.sectors) {
      if (!material_sector(sector)) continue;
      if (sector.b.is_zero != TruthValue::Yes || sector.log_power != 0 ||
          sector.a.is_integer != TruthValue::Yes ||
          sector.a.sign == ExactSign::Negative ||
          sector.a.sign == ExactSign::Unknown)
        throw std::domain_error(
            "singular-center evaluation requires endpoint-limit semantics");
    }
  }

  const auto full_taylor_width = solution.taylor_width();
  if (options.t_order_reduction >= full_taylor_width)
    throw std::invalid_argument("t-order reduction removes every coefficient");
  const auto taylor_width = full_taylor_width - options.t_order_reduction;

  std::int64_t complete_max64 = std::numeric_limits<std::int64_t>::max();
  std::int64_t min_power64 = std::numeric_limits<std::int64_t>::max();
  for (const auto& sector : solution.sectors) {
    complete_max64 = std::min(complete_max64,
        static_cast<std::int64_t>(solution.epsilon.complete_max) +
            sector.log_power);
    std::int32_t first = solution.epsilon.complete_max;
    for (std::int32_t k = solution.epsilon.min_power;
         k <= solution.epsilon.complete_max; ++k) {
      const auto ei = static_cast<std::size_t>(k - solution.epsilon.min_power);
      bool row_nonzero = false;
      for (std::size_t n = 0; n < taylor_width && !row_nonzero; ++n)
        for (std::uint32_t c = 0; c < solution.dimension; ++c)
          if (!ScalarTraits<Scalar>::is_zero(
                  sector.coefficients[sector_index(solution, ei, n, c)])) {
            row_nonzero = true;
            break;
          }
      if (row_nonzero) {
        first = k;
        break;
      }
    }
    min_power64 = std::min(min_power64,
        static_cast<std::int64_t>(first) + sector.log_power);
  }
  const EpsilonWindow output_window{
      checked_i32(min_power64, "evaluated minimum epsilon power"),
      checked_i32(complete_max64, "evaluated complete epsilon maximum")};
  EpsilonVector output;
  output.epsilon = output_window;
  output.dimension = solution.dimension;
  output.coefficients.assign(output_window.width() * solution.dimension,
                             ComplexBall(0));

  ComplexBall logarithm(0);
  if (point.sign != 0) {
    logarithm = cb_log(point.modulus);
    if (point.sign < 0 && sigma.has_value())
      logarithm += imaginary_pi(*sigma);
  }
  ComplexBall signed_t = point.modulus;
  if (point.sign < 0) signed_t = -signed_t;
  std::vector<ComplexBall> t_powers(taylor_width, ComplexBall(1));
  for (std::size_t n = 1; n < taylor_width; ++n)
    t_powers[n] = t_powers[n - 1] * signed_t;

  std::vector<Magnitude> tail(output_window.width(), Magnitude::zero());
  ComplexBall geometric(0);
  if (options.compute_tail_estimate && point.sign != 0 &&
      !solution.chart.infinite_radius) {
    const auto q = point.modulus / solution.chart.radius;
    geometric = q / (ComplexBall(1) - q);
  }

  for (const auto& sector : solution.sectors) {
    const auto& a = sector.a.numeric();
    const auto& b = sector.b.numeric();
    ComplexBall t_to_a(0);
    if (point.sign == 0) {
      t_to_a = sector.a.is_zero == TruthValue::Yes ? ComplexBall(1)
                                                   : ComplexBall(0);
    } else {
      t_to_a = cb_pow(point.modulus, a);
      if (point.sign < 0) {
        // For wholly single-valued integer sectors no prescription is needed;
        // either rim gives the same integer phase.  +1 is therefore only a
        // computational representative here, never an inferred sheet.
        const auto phase = cb_exp(imaginary_pi(sigma.value_or(1)) * a);
        t_to_a *= phase;
      }
    }
    auto outer = t_to_a * cb_pow_ui(logarithm, sector.log_power);
    for (std::uint32_t d = 2; d <= sector.log_power; ++d)
      outer = cb_div_ui(outer, d);
    const auto b_log = b * logarithm;

    const auto eps_width = solution.epsilon.width();
    std::vector<ComplexBall> alpha(eps_width * solution.dimension,
                                   ComplexBall(0));
    std::vector<Magnitude> top(eps_width, Magnitude::zero());
    for (std::size_t ei = 0; ei < eps_width; ++ei) {
      for (std::uint32_t c = 0; c < solution.dimension; ++c) {
        ComplexBall sum(0);
        for (std::size_t n = 0; n < taylor_width; ++n)
          sum += to_ball(sector.coefficients[
                     sector_index(solution, ei, n, c)]) * t_powers[n];
        alpha[ei * solution.dimension + c] = std::move(sum);
        if (options.compute_tail_estimate) {
          top[ei] = Magnitude::maximum(top[ei], Magnitude::upper_abs(
              to_ball(sector.coefficients[sector_index(
                  solution, ei, taylor_width - 1, c)])));
        }
      }
    }

    for (std::int32_t K = output_window.min_power;
         K <= output_window.complete_max; ++K) {
      ComplexBall exponential_coefficient(1);
      const auto jmax64 = static_cast<std::int64_t>(K) -
                          sector.log_power - solution.epsilon.min_power;
      if (jmax64 < 0) continue;
      for (std::int64_t j = 0; j <= jmax64; ++j) {
        const auto input_power64 = static_cast<std::int64_t>(K) -
                                   sector.log_power - j;
        if (input_power64 >= solution.epsilon.min_power &&
            input_power64 <= solution.epsilon.complete_max) {
          const auto ei = static_cast<std::size_t>(
              input_power64 - solution.epsilon.min_power);
          for (std::uint32_t c = 0; c < solution.dimension; ++c)
            output.at(K, c) += outer * exponential_coefficient *
                               alpha[ei * solution.dimension + c];
          if (options.compute_tail_estimate && point.sign != 0 &&
              !solution.chart.infinite_radius) {
            const auto top_term = top[ei] *
                Magnitude::upper_abs(t_powers.back()) *
                Magnitude::upper_abs(geometric) *
                Magnitude::upper_abs(outer) *
                Magnitude::upper_abs(exponential_coefficient);
            tail[static_cast<std::size_t>(K - output_window.min_power)] +=
                top_term;
          }
        }
        exponential_coefficient = cb_div_ui(
            exponential_coefficient * b_log,
            static_cast<std::uint32_t>(j + 1));
      }
    }
  }

  if (options.compute_tail_estimate) {
    output.error.frame = output_window;
    output.error.guarantee = ErrorGuarantee::Advisory;
    output.error.absolute = std::move(tail);
    output.error.provenance =
        "geometric top-column model; Acb arithmetic is enclosed but the "
        "coefficient-ratio hypothesis is not yet certified";
  }
  return output;
}

template <typename Scalar>
LocalSolution<ComplexBall> theta_solution(const LocalSolution<Scalar>& input) {
  local_detail::validate_local_solution(input, true);
  LocalSolution<ComplexBall> output;
  output.chart = input.chart;
  output.epsilon = input.epsilon;
  output.taylor_complete_max = input.taylor_complete_max;
  output.dimension = input.dimension;
  output.prescriptions = input.prescriptions;
  output.checkpoint_identity = input.checkpoint_identity + ":theta";
  const auto eps_width = input.epsilon.width();
  const auto taylor_width = input.taylor_width();
  for (const auto& sector : input.sectors) {
    LocalSector<ComplexBall> same;
    same.a = sector.a;
    same.b = sector.b;
    same.log_power = sector.log_power;
    same.coefficients.assign(input.sector_size(), ComplexBall(0));
    const auto& a = sector.a.numeric();
    const auto& b = sector.b.numeric();
    for (std::size_t ei = 0; ei < eps_width; ++ei) {
      for (std::size_t n = 0; n < taylor_width; ++n) {
        const auto a_plus_n = a + ComplexBall(static_cast<long>(n));
        for (std::uint32_t c = 0; c < input.dimension; ++c) {
          const auto index = sector_index(input, ei, n, c);
          same.coefficients[index] =
              a_plus_n * to_ball(sector.coefficients[index]);
          if (ei > 0) {
            const auto previous = sector_index(input, ei - 1, n, c);
            same.coefficients[index] += b * to_ball(sector.coefficients[previous]);
          }
        }
      }
    }
    output.sectors.push_back(std::move(same));
    if (sector.log_power > 0) {
      LocalSector<ComplexBall> lower;
      lower.a = sector.a;
      lower.b = sector.b;
      lower.log_power = sector.log_power - 1;
      lower.coefficients.assign(input.sector_size(), ComplexBall(0));
      for (std::size_t ei = 1; ei < eps_width; ++ei)
        for (std::size_t n = 0; n < taylor_width; ++n)
          for (std::uint32_t c = 0; c < input.dimension; ++c)
            lower.coefficients[sector_index(input, ei, n, c)] =
                to_ball(sector.coefficients[
                    sector_index(input, ei - 1, n, c)]);
      output.sectors.push_back(std::move(lower));
    }
  }
  return output;
}

inline ComplexBall coefficient_or_zero(const EpsilonVector& value,
                                       std::int32_t power,
                                       std::uint32_t component) {
  if (power < value.epsilon.min_power ||
      power > value.epsilon.complete_max)
    return ComplexBall(0);
  return value.at(power, component);
}

inline EpsilonVector matrix_times_vector(const EpsilonMatrix& matrix,
                                         const EpsilonVector& vector) {
  if (matrix.dimension == 0 || matrix.dimension != vector.dimension)
    throw std::invalid_argument("matrix/vector dimensions disagree");
  if (matrix.coefficients.size() != matrix.epsilon.width() * matrix.dimension *
                                      matrix.dimension)
    throw std::invalid_argument("epsilon matrix tensor has wrong dimensions");
  const auto min64 = static_cast<std::int64_t>(matrix.epsilon.min_power) +
                     vector.epsilon.min_power;
  const auto max64 = std::min(
      static_cast<std::int64_t>(matrix.epsilon.complete_max) +
          vector.epsilon.min_power,
      static_cast<std::int64_t>(vector.epsilon.complete_max) +
          matrix.epsilon.min_power);
  EpsilonVector output;
  output.epsilon = {checked_i32(min64, "matrix product minimum"),
                    checked_i32(max64, "matrix product complete maximum")};
  output.dimension = vector.dimension;
  output.coefficients.assign(output.epsilon.width() * output.dimension,
                             ComplexBall(0));
  for (std::int32_t k = output.epsilon.min_power;
       k <= output.epsilon.complete_max; ++k) {
    for (std::int32_t ka = matrix.epsilon.min_power;
         ka <= matrix.epsilon.complete_max; ++ka) {
      const auto kf = k - ka;
      if (kf < vector.epsilon.min_power ||
          kf > vector.epsilon.complete_max)
        continue;
      for (std::uint32_t r = 0; r < matrix.dimension; ++r)
        for (std::uint32_t c = 0; c < matrix.dimension; ++c)
          output.at(k, r) += matrix.at(ka, r, c) * vector.at(kf, c);
    }
  }
  return output;
}

}  // namespace local_detail

template <typename Scalar>
void validate_local_solution(const LocalSolution<Scalar>& solution,
                             bool require_numeric_tags = false) {
  local_detail::validate_local_solution(solution, require_numeric_tags);
}

template <typename Scalar>
std::optional<std::int32_t> derive_chart_imaginary_sign(
    const LocalSolution<Scalar>& solution) {
  local_detail::validate_local_solution(solution, false);
  return local_detail::chart_imaginary_sign(solution);
}

template <typename Scalar>
LocalEvaluation evaluate_local_solution(
    const LocalSolution<Scalar>& solution, const RealEvaluationPoint& point,
    const EvaluationOptions& options = {}) {
  LocalEvaluation output;
  output.value = local_detail::evaluate_value(
      solution, point, options, &output.imaginary_sign);
  auto theta = local_detail::theta_solution(solution);
  output.theta_value = local_detail::evaluate_value(
      theta, point, options, nullptr);
  return output;
}

enum class ResidualScope : std::uint8_t { StoredTruncation, FullLocalSolution };
enum class ResidualVerdict : std::uint8_t { Pass, Fail, Inconclusive };

struct ResidualCertificate {
  EpsilonVector residual;
  std::vector<Magnitude> residual_upper;
  std::vector<Magnitude> scale_lower;
  std::vector<Magnitude> relative_upper;
  ResidualScope scope = ResidualScope::StoredTruncation;
  ResidualVerdict verdict = ResidualVerdict::Inconclusive;
  std::string detail;
};

inline ResidualCertificate certify_theta_residual(
    const LocalEvaluation& evaluation, const EpsilonMatrix& theta_operator,
    const std::optional<EpsilonVector>& source,
    const Magnitude& relative_tolerance,
    ResidualScope scope = ResidualScope::StoredTruncation) {
  if (evaluation.value.dimension != evaluation.theta_value.dimension ||
      theta_operator.dimension != evaluation.value.dimension)
    throw std::invalid_argument("residual dimensions disagree");
  if (source.has_value() && source->dimension != evaluation.value.dimension)
    throw std::invalid_argument("residual source dimension disagrees");
  auto rhs = local_detail::matrix_times_vector(theta_operator, evaluation.value);
  const auto min_power = std::min({evaluation.theta_value.epsilon.min_power,
                                   rhs.epsilon.min_power,
                                   source.has_value()
                                       ? source->epsilon.min_power
                                       : std::numeric_limits<std::int32_t>::max()});
  const auto complete_max = std::min({
      evaluation.theta_value.epsilon.complete_max, rhs.epsilon.complete_max,
      source.has_value() ? source->epsilon.complete_max
                         : std::numeric_limits<std::int32_t>::max()});
  if (complete_max < min_power)
    throw std::domain_error("residual has no common complete epsilon window");

  ResidualCertificate result;
  result.scope = scope;
  result.residual.epsilon = {min_power, complete_max};
  result.residual.dimension = evaluation.value.dimension;
  result.residual.coefficients.assign(
      result.residual.epsilon.width() * result.residual.dimension,
      ComplexBall(0));
  result.residual_upper.assign(result.residual.epsilon.width(),
                               Magnitude::zero());
  result.scale_lower.assign(result.residual.epsilon.width(), Magnitude::one());
  std::vector<Magnitude> scale_upper(result.residual.epsilon.width(),
                                     Magnitude::one());
  std::vector<Magnitude> residual_lower(result.residual.epsilon.width(),
                                        Magnitude::zero());
  result.relative_upper.assign(result.residual.epsilon.width(),
                               Magnitude::zero());

  bool all_pass = true;
  bool any_fail = false;
  for (std::int32_t k = min_power; k <= complete_max; ++k) {
    const auto row = static_cast<std::size_t>(k - min_power);
    for (std::uint32_t c = 0; c < result.residual.dimension; ++c) {
      const auto lhs = local_detail::coefficient_or_zero(
          evaluation.theta_value, k, c);
      const auto rhs_value = local_detail::coefficient_or_zero(rhs, k, c);
      const auto source_value = source.has_value()
          ? local_detail::coefficient_or_zero(*source, k, c)
          : ComplexBall(0);
      auto residual = lhs - rhs_value - source_value;
      result.residual.at(k, c) = residual;
      result.residual_upper[row] = Magnitude::maximum(
          result.residual_upper[row], Magnitude::upper_abs(residual));
      residual_lower[row] = Magnitude::maximum(
          residual_lower[row], Magnitude::lower_abs(residual));
      for (const auto& term : {lhs, rhs_value, source_value}) {
        result.scale_lower[row] = Magnitude::maximum(
            result.scale_lower[row], Magnitude::lower_abs(term));
        scale_upper[row] = Magnitude::maximum(
            scale_upper[row], Magnitude::upper_abs(term));
      }
    }
    result.relative_upper[row] =
        result.residual_upper[row] / result.scale_lower[row];
    const bool row_pass = result.residual_upper[row] <=
                          relative_tolerance * result.scale_lower[row];
    const bool row_fail = residual_lower[row] >
                          relative_tolerance * scale_upper[row];
    all_pass = all_pass && row_pass;
    any_fail = any_fail || row_fail;
  }

  if (scope == ResidualScope::FullLocalSolution) {
    // The current Wolfram geometric tail is a model, not a proof.  Until a
    // certified recurrence majorant is installed and propagated through A,
    // never promote a polynomial enclosure to a full-solution certificate.
    result.verdict = ResidualVerdict::Inconclusive;
    result.detail =
        "Acb encloses the stored-truncation residual, but full local-solution "
        "certification requires certified value/theta tail bounds";
  } else if (all_pass) {
    result.verdict = ResidualVerdict::Pass;
    result.detail = "Acb enclosure proves the stored-truncation residual bound";
  } else if (any_fail) {
    result.verdict = ResidualVerdict::Fail;
    result.detail = "Acb lower bound disproves the stored-truncation tolerance";
  } else {
    result.verdict = ResidualVerdict::Inconclusive;
    result.detail = "residual enclosure overlaps the requested tolerance";
  }
  return result;
}

}  // namespace diffexp2
