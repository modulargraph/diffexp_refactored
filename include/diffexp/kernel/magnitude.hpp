#pragma once

#include "diffexp/kernel/scalar.hpp"

#include <flint/arb.h>
#include <flint/mag.h>

#include <new>
#include <stdexcept>
#include <string>
#include <utility>

namespace diffexp::kernel {

// An upward-rounded nonnegative bound (FLINT mag_t) with value semantics.
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

  // Absorb this certified absolute bound into an Acb coefficient.  Keeping
  // the operation here avoids exposing FLINT's raw magnitude representation
  // to transport code and makes the one-way certificate-to-ball handoff
  // explicit.
  void add_error_to(ComplexBall& value) const {
    acb_add_error_mag(value.raw(), value_);
  }

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

}  // namespace diffexp::kernel
