#pragma once
#include "diffexp/kernel/scalar.hpp"
#include <cmath>
#include <stdexcept>

namespace diffexp {
// Arithmetic difference between retained coefficients, not a solver-tail or
// reference-error certificate. Reject indeterminate values before reduction:
// FLINT can return NaN, which std::max(previous, NaN) would silently ignore.
inline double finite_reference_error(const kernel::ComplexBall& value,
    const kernel::ComplexBall& reference, slong bits) {
  if(bits<2)throw std::invalid_argument("reference comparison precision must be at least two bits");
  if(!value.is_finite() || !reference.is_finite())
    throw std::domain_error("nonfinite coefficient in numerical reference comparison");
  kernel::ComplexBall difference;
  acb_sub(difference.raw(),value.raw(),reference.raw(),bits);
  if(!difference.is_finite())throw std::domain_error("nonfinite reference difference");
  arf_t upper;arf_init(upper);
  acb_get_abs_ubound_arf(upper,difference.raw(),bits);
  const double result=arf_get_d(upper,ARF_RND_CEIL);
  arf_clear(upper);
  if(!std::isfinite(result))throw std::domain_error("nonfinite reference error bound");
  return result;
}
} // namespace diffexp
