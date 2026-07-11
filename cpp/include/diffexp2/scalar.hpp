#pragma once

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz_mpoly_q.h>

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace diffexp2 {

class Rational {
 public:
  Rational() { fmpq_init(value_); }
  explicit Rational(long value) : Rational() { fmpq_set_si(value_, value, 1); }
  explicit Rational(const std::string& value) : Rational() {
    if (fmpq_set_str(value_, value.c_str(), 10) != 0) {
      throw std::invalid_argument("invalid rational: " + value);
    }
    fmpq_canonicalise(value_);
  }
  Rational(const Rational& other) : Rational() { fmpq_set(value_, other.value_); }
  Rational(Rational&& other) noexcept : Rational() { fmpq_swap(value_, other.value_); }
  ~Rational() { fmpq_clear(value_); }

  Rational& operator=(Rational other) noexcept {
    fmpq_swap(value_, other.value_);
    return *this;
  }

  [[nodiscard]] bool is_zero() const { return fmpq_is_zero(value_); }
  [[nodiscard]] std::string str() const {
    char* raw = fmpq_get_str(nullptr, 10, value_);
    if (raw == nullptr) throw std::bad_alloc();
    std::string result(raw);
    flint_free(raw);
    return result;
  }

  friend Rational operator+(const Rational& a, const Rational& b) {
    Rational out;
    fmpq_add(out.value_, a.value_, b.value_);
    return out;
  }
  friend Rational operator-(const Rational& a, const Rational& b) {
    Rational out;
    fmpq_sub(out.value_, a.value_, b.value_);
    return out;
  }
  friend Rational operator-(const Rational& a) {
    Rational out;
    fmpq_neg(out.value_, a.value_);
    return out;
  }
  friend Rational operator*(const Rational& a, const Rational& b) {
    Rational out;
    fmpq_mul(out.value_, a.value_, b.value_);
    return out;
  }
  friend Rational operator/(const Rational& a, const Rational& b) {
    if (b.is_zero()) throw std::domain_error("rational division by zero");
    Rational out;
    fmpq_div(out.value_, a.value_, b.value_);
    return out;
  }
  Rational& operator+=(const Rational& b) {
    fmpq_add(value_, value_, b.value_);
    return *this;
  }
  Rational& operator-=(const Rational& b) {
    fmpq_sub(value_, value_, b.value_);
    return *this;
  }
  Rational& operator*=(const Rational& b) {
    fmpq_mul(value_, value_, b.value_);
    return *this;
  }
  friend bool operator==(const Rational& a, const Rational& b) {
    return fmpq_equal(a.value_, b.value_);
  }

 private:
  fmpq_t value_;
};

class SymbolicRational {
 public:
  static void configure(const std::vector<std::string>& variables) {
    auto& holder = context_holder();
    if (holder.initialized && holder.names == variables) return;
    if (holder.initialized && holder.live_objects != 0) {
      throw std::logic_error(
          "cannot change symbolic coefficient field while values are alive");
    }
    if (holder.initialized) fmpz_mpoly_ctx_clear(holder.context);
    holder.names = variables;
    holder.c_names.clear();
    holder.c_names.reserve(holder.names.size());
    for (const auto& name : holder.names) holder.c_names.push_back(name.c_str());
    fmpz_mpoly_ctx_init(holder.context,
                        static_cast<slong>(holder.names.size()), ORD_DEGLEX);
    holder.initialized = true;
  }

  SymbolicRational() {
    ensure_context();
    fmpz_mpoly_q_init(value_, context_holder().context);
    ++context_holder().live_objects;
  }
  explicit SymbolicRational(long value) : SymbolicRational() {
    fmpz_mpoly_q_set_si(value_, value, context_holder().context);
  }
  explicit SymbolicRational(const std::string& value) : SymbolicRational() {
    auto& holder = context_holder();
    if (fmpz_mpoly_q_set_str_pretty(value_, value.c_str(),
          holder.c_names.data(), holder.context) != 0) {
      throw std::invalid_argument("invalid symbolic rational: " + value);
    }
    fmpz_mpoly_q_canonicalise(value_, holder.context);
  }
  SymbolicRational(const SymbolicRational& other) : SymbolicRational() {
    fmpz_mpoly_q_set(value_, other.value_, context_holder().context);
  }
  SymbolicRational(SymbolicRational&& other) noexcept : SymbolicRational() {
    fmpz_mpoly_q_swap(value_, other.value_, context_holder().context);
  }
  ~SymbolicRational() {
    fmpz_mpoly_q_clear(value_, context_holder().context);
    --context_holder().live_objects;
  }

  SymbolicRational& operator=(SymbolicRational other) noexcept {
    fmpz_mpoly_q_swap(value_, other.value_, context_holder().context);
    return *this;
  }

  [[nodiscard]] bool is_zero() const {
    return fmpz_mpoly_q_is_zero(value_, context_holder().context);
  }
  [[nodiscard]] std::string str() const {
    auto& holder = context_holder();
    char* raw = fmpz_mpoly_q_get_str_pretty(
        value_, holder.c_names.data(), holder.context);
    if (raw == nullptr) throw std::bad_alloc();
    std::string result(raw);
    flint_free(raw);
    return result;
  }

  friend SymbolicRational operator+(const SymbolicRational& a,
                                    const SymbolicRational& b) {
    SymbolicRational out;
    fmpz_mpoly_q_add(out.value_, a.value_, b.value_, context_holder().context);
    return out;
  }
  friend SymbolicRational operator-(const SymbolicRational& a,
                                    const SymbolicRational& b) {
    SymbolicRational out;
    fmpz_mpoly_q_sub(out.value_, a.value_, b.value_, context_holder().context);
    return out;
  }
  friend SymbolicRational operator-(const SymbolicRational& a) {
    SymbolicRational out;
    fmpz_mpoly_q_neg(out.value_, a.value_, context_holder().context);
    return out;
  }
  friend SymbolicRational operator*(const SymbolicRational& a,
                                    const SymbolicRational& b) {
    SymbolicRational out;
    fmpz_mpoly_q_mul(out.value_, a.value_, b.value_, context_holder().context);
    return out;
  }
  friend SymbolicRational operator/(const SymbolicRational& a,
                                    const SymbolicRational& b) {
    if (b.is_zero()) throw std::domain_error("symbolic rational division by zero");
    SymbolicRational out;
    fmpz_mpoly_q_div(out.value_, a.value_, b.value_, context_holder().context);
    return out;
  }
  SymbolicRational& operator+=(const SymbolicRational& b) {
    fmpz_mpoly_q_add(value_, value_, b.value_, context_holder().context);
    return *this;
  }
  SymbolicRational& operator-=(const SymbolicRational& b) {
    fmpz_mpoly_q_sub(value_, value_, b.value_, context_holder().context);
    return *this;
  }
  SymbolicRational& operator*=(const SymbolicRational& b) {
    fmpz_mpoly_q_mul(value_, value_, b.value_, context_holder().context);
    return *this;
  }
  friend bool operator==(const SymbolicRational& a,
                         const SymbolicRational& b) {
    return fmpz_mpoly_q_equal(a.value_, b.value_, context_holder().context);
  }

 private:
  struct ContextHolder {
    fmpz_mpoly_ctx_t context;
    bool initialized = false;
    std::size_t live_objects = 0;
    std::vector<std::string> names;
    std::vector<const char*> c_names;
    ~ContextHolder() {
      if (initialized) fmpz_mpoly_ctx_clear(context);
    }
  };
  static ContextHolder& context_holder() {
    static ContextHolder holder;
    return holder;
  }
  static void ensure_context() {
    if (!context_holder().initialized) configure({});
  }

  fmpz_mpoly_q_t value_;
};

class ComplexBall {
 public:
  ComplexBall() { acb_init(value_); }
  explicit ComplexBall(long value) : ComplexBall() { acb_set_si(value_, value); }
  ComplexBall(const ComplexBall& other) : ComplexBall() { acb_set(value_, other.value_); }
  ComplexBall(ComplexBall&& other) noexcept : ComplexBall() { acb_swap(value_, other.value_); }
  ~ComplexBall() { acb_clear(value_); }

  ComplexBall& operator=(ComplexBall other) noexcept {
    acb_swap(value_, other.value_);
    return *this;
  }

  static void set_precision(slong bits) {
    if (bits < 64) throw std::invalid_argument("Acb precision must be at least 64 bits");
    precision_bits_ = bits;
  }
  static slong precision() { return precision_bits_; }

  static ComplexBall from_strings(const std::string& real,
                                  const std::string& imag = "0") {
    ComplexBall out;
    set_arb_string(acb_realref(out.value_), real);
    set_arb_string(acb_imagref(out.value_), imag);
    return out;
  }

  [[nodiscard]] bool is_zero() const { return acb_is_zero(value_); }
  [[nodiscard]] bool contains_zero() const { return acb_contains_zero(value_); }
  [[nodiscard]] bool is_finite() const { return acb_is_finite(value_); }

  // This mirrors the Wolfram certified-frame chop: only an enclosure wholly
  // below the configured absolute floor is replaced by structural zero.
  bool chop_if_certified(int decimal_digits) {
    if (is_zero()) return true;
    arf_t upper;
    arf_init(upper);
    acb_get_abs_ubound_arf(upper, value_, precision_bits_);
    const auto bits = static_cast<slong>(
        std::ceil(static_cast<double>(decimal_digits) * std::log2(10.0)));
    const bool tiny = arf_cmp_2exp_si(upper, -bits) < 0;
    arf_clear(upper);
    if (tiny) acb_zero(value_);
    return tiny;
  }

  [[nodiscard]] std::string real_midpoint(int digits) const {
    return arb_midpoint_string(acb_realref(value_), digits);
  }
  [[nodiscard]] std::string imag_midpoint(int digits) const {
    return arb_midpoint_string(acb_imagref(value_), digits);
  }
  [[nodiscard]] std::string real_radius_exponent() const {
    return mag_radius_exponent(arb_radref(acb_realref(value_)));
  }
  [[nodiscard]] std::string imag_radius_exponent() const {
    return mag_radius_exponent(arb_radref(acb_imagref(value_)));
  }

  const acb_struct* raw() const { return value_; }
  acb_struct* raw() { return value_; }

  friend ComplexBall operator+(const ComplexBall& a, const ComplexBall& b) {
    ComplexBall out;
    acb_add(out.value_, a.value_, b.value_, precision_bits_);
    return out;
  }
  friend ComplexBall operator-(const ComplexBall& a, const ComplexBall& b) {
    ComplexBall out;
    acb_sub(out.value_, a.value_, b.value_, precision_bits_);
    return out;
  }
  friend ComplexBall operator-(const ComplexBall& a) {
    ComplexBall out;
    acb_neg(out.value_, a.value_);
    return out;
  }
  friend ComplexBall operator*(const ComplexBall& a, const ComplexBall& b) {
    ComplexBall out;
    acb_mul(out.value_, a.value_, b.value_, precision_bits_);
    return out;
  }
  friend ComplexBall operator/(const ComplexBall& a, const ComplexBall& b) {
    if (b.contains_zero()) throw std::domain_error("Acb division by enclosure containing zero");
    ComplexBall out;
    acb_div(out.value_, a.value_, b.value_, precision_bits_);
    return out;
  }
  ComplexBall& operator+=(const ComplexBall& b) {
    acb_add(value_, value_, b.value_, precision_bits_);
    return *this;
  }
  ComplexBall& operator-=(const ComplexBall& b) {
    acb_sub(value_, value_, b.value_, precision_bits_);
    return *this;
  }
  ComplexBall& operator*=(const ComplexBall& b) {
    acb_mul(value_, value_, b.value_, precision_bits_);
    return *this;
  }

 private:
  static void set_arb_string(arb_t out, const std::string& value) {
    if (!value.empty() && value.front() != '[' &&
        value.find('/') != std::string::npos) {
      fmpq_t q;
      fmpq_init(q);
      const int status = fmpq_set_str(q, value.c_str(), 10);
      if (status == 0) {
        fmpq_canonicalise(q);
        arb_set_fmpq(out, q, precision_bits_);
      }
      fmpq_clear(q);
      if (status != 0) throw std::invalid_argument("invalid rational scalar: " + value);
      return;
    }
    if (arb_set_str(out, value.c_str(), precision_bits_) != 0) {
      throw std::invalid_argument("invalid arbitrary-precision scalar: " + value);
    }
  }

  static std::string arb_midpoint_string(const arb_t value, int digits) {
    char* raw = arb_get_str(value, digits, ARB_STR_NO_RADIUS | ARB_STR_MORE);
    if (raw == nullptr) throw std::bad_alloc();
    std::string result(raw);
    flint_free(raw);
    return result;
  }

  // A finite nonzero mag is man * 2^(exp-MAG_BITS), with man < 2^MAG_BITS,
  // hence radius < 2^exp. Exporting only this exact integer exponent keeps
  // the bridge record compact at thousands of digits while giving Wolfram
  // a certified absolute-error bound. "zero" is the exact-radius sentinel.
  static std::string mag_radius_exponent(const mag_t value) {
    if (mag_is_zero(value)) return "zero";
    if (mag_is_inf(value))
      throw std::domain_error("non-finite Arb error radius");
    char* raw = fmpz_get_str(nullptr, 10, MAG_EXPREF(value));
    if (raw == nullptr) throw std::bad_alloc();
    std::string result(raw);
    flint_free(raw);
    return result;
  }

  inline static thread_local slong precision_bits_ = 256;
  acb_t value_;
};

template <typename Scalar>
struct ScalarTraits;

template <>
struct ScalarTraits<Rational> {
  static Rational zero() { return Rational(0); }
  static Rational one() { return Rational(1); }
  static Rational integer(long value) { return Rational(value); }
  static bool is_zero(const Rational& value) { return value.is_zero(); }
  static bool certified_zero(const Rational& value, int) {
    return value.is_zero();
  }
  static Rational canonicalized(const Rational& value, int) { return value; }
};

template <>
struct ScalarTraits<ComplexBall> {
  static ComplexBall zero() { return ComplexBall(0); }
  static ComplexBall one() { return ComplexBall(1); }
  static ComplexBall integer(long value) { return ComplexBall(value); }
  static bool is_zero(const ComplexBall& value) { return value.is_zero(); }
  static bool certified_zero(const ComplexBall& value, int chop_digits) {
    auto copy = value;
    return copy.chop_if_certified(chop_digits);
  }
  static ComplexBall canonicalized(const ComplexBall& value, int chop_digits) {
    auto copy = value;
    copy.chop_if_certified(chop_digits);
    return copy;
  }
};

template <>
struct ScalarTraits<SymbolicRational> {
  static SymbolicRational zero() { return SymbolicRational(0); }
  static SymbolicRational one() { return SymbolicRational(1); }
  static SymbolicRational integer(long value) { return SymbolicRational(value); }
  static bool is_zero(const SymbolicRational& value) { return value.is_zero(); }
  static bool certified_zero(const SymbolicRational& value, int) {
    return value.is_zero();
  }
  static SymbolicRational canonicalized(const SymbolicRational& value, int) {
    return value;
  }
};

}  // namespace diffexp2
