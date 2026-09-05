#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace diffexp::kernel {

// Honest finite Laurent frame shared by recurrence output, local evaluation,
// matching, endpoint limits and integration.  Coefficients below min_power
// are structural zero; coefficients above complete_max are unknown, not zero.
struct EpsilonWindow {
  std::int32_t min_power = 0;
  std::int32_t complete_max = 0;

  [[nodiscard]] std::size_t width() const {
    if (complete_max < min_power)
      throw std::invalid_argument("empty epsilon window");
    return static_cast<std::size_t>(
        static_cast<std::int64_t>(complete_max) - min_power + 1);
  }
};

}  // namespace diffexp::kernel
