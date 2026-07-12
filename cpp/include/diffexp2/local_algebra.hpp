#pragma once

#include "diffexp2/local_solution.hpp"
#include "diffexp2/matching.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace diffexp2 {

// Prepared exact rational multiplication for a finite local-sector slab.
// Wolfram still owns exact pole classification in the first SCC milestone:
// it expands c(t,eps) as
//
//   eps^epsilon_shift * t^-center_pole_order
//     * sum_j eps^j sum_n kernels[j][n] t^n.
//
// The native operation below owns only the finite triangular convolutions.
// `kernels` must contain at least the input epsilon width and every row must
// contain at least the input Taylor width.  This matches
// SectorSeries`MultiplyRational: the output has the same number of epsilon
// rows, with both ends shifted by the exact leading epsilon valuation.
template <typename Scalar>
struct PreparedRationalTaylorMultiplier {
  std::int32_t epsilon_shift = 0;
  std::uint32_t center_pole_order = 0;
  std::vector<std::vector<Scalar>> kernels;  // [epsilon][Taylor]
  std::string exact_identity;
  // Exact preparation provenance, never a conclusion drawn from a numeric
  // slab.  An Acb tensor whose entries currently enclose exact zero remains
  // active unless the rational multiplier was proved structurally zero
  // before specialization.
  bool proven_zero = false;

  [[nodiscard]] bool structurally_zero() const {
    return proven_zero;
  }
};

template <typename Scalar>
struct PreparedSparseLocalMultiplierMatrix {
  struct Entry {
    std::uint32_t row = 0;
    std::uint32_t column = 0;
    PreparedRationalTaylorMultiplier<Scalar> multiplier;
  };

  std::uint32_t rows = 0;
  std::uint32_t columns = 0;
  std::vector<Entry> entries;
  std::string exact_identity;
};

namespace local_algebra_detail {

inline std::int32_t checked_i32(std::int64_t value, const char* label) {
  if (value < std::numeric_limits<std::int32_t>::min() ||
      value > std::numeric_limits<std::int32_t>::max())
    throw std::overflow_error(std::string(label) + " exceeds int32 range");
  return static_cast<std::int32_t>(value);
}

inline bool same_descriptor(const ExactScalarDescriptor& left,
                            const ExactScalarDescriptor& right) {
  if (left.domain != right.domain || left.canonical != right.canonical ||
      left.symbols != right.symbols || left.is_zero != right.is_zero ||
      left.is_integer != right.is_integer || left.sign != right.sign ||
      left.specialization.has_value() != right.specialization.has_value())
    return false;
  return !left.specialization.has_value() ||
         acb_equal(left.specialization->raw(), right.specialization->raw());
}

inline bool same_prescriptions(const std::vector<Prescription>& left,
                               const std::vector<Prescription>& right) {
  if (left.size() != right.size()) return false;
  for (std::size_t i = 0; i < left.size(); ++i) {
    if (left[i].factor_exact != right[i].factor_exact ||
        left[i].sign != right[i].sign ||
        left[i].multiplicity != right[i].multiplicity ||
        left[i].leading_coefficient_sign !=
            right[i].leading_coefficient_sign)
      return false;
  }
  return true;
}

inline bool same_chart(const ChartGeometry& left,
                       const ChartGeometry& right) {
  if (left.center_exact != right.center_exact ||
      left.scale_exact != right.scale_exact ||
      left.infinite_radius != right.infinite_radius)
    return false;
  return left.infinite_radius || acb_equal(left.radius.raw(), right.radius.raw());
}

template <typename Scalar>
void require_same_local_space(const LocalSolution<Scalar>& left,
                              const LocalSolution<Scalar>& right) {
  if (!same_chart(left.chart, right.chart) ||
      !same_prescriptions(left.prescriptions, right.prescriptions) ||
      left.dimension != right.dimension)
    throw std::invalid_argument(
        "local algebra requires identical chart, prescription and component spaces");
}

inline std::size_t flat_index(std::size_t epsilon_index,
                              std::size_t taylor_index,
                              std::uint32_t component,
                              std::size_t taylor_width,
                              std::uint32_t dimension) {
  return ((epsilon_index * taylor_width + taylor_index) * dimension) +
         component;
}

inline ExactScalarDescriptor subtract_nonnegative_integer(
    const ExactScalarDescriptor& input, std::uint32_t amount) {
  if (amount == 0) return input;
  const auto amount_string = std::to_string(amount);
  if (input.domain == ExactDomain::Rational) {
    const auto value = Rational(input.canonical) - Rational(amount_string);
    return ExactScalarDescriptor::rational(value.str());
  }

  auto output = input;
  output.canonical = "(" + input.canonical + ")-" + amount_string;
  // Subtracting an integer preserves integer/noninteger status, but zero and
  // sign need an exact comparison which this descriptor layer must not guess.
  output.is_zero = TruthValue::Unknown;
  output.sign = ExactSign::Unknown;
  if (input.specialization.has_value())
    output.specialization = *input.specialization -
        ComplexBall::from_strings(amount_string);
  return output;
}

template <typename Scalar>
bool same_sector_tag(const LocalSector<Scalar>& left,
                     const LocalSector<Scalar>& right) {
  return left.log_power == right.log_power &&
         same_descriptor(left.a, right.a) &&
         same_descriptor(left.b, right.b);
}

template <typename Scalar>
bool sector_less(const LocalSector<Scalar>& left,
                 const LocalSector<Scalar>& right) {
  return std::tie(left.a.domain, left.a.canonical, left.b.domain,
                  left.b.canonical, left.log_power) <
         std::tie(right.a.domain, right.a.canonical, right.b.domain,
                  right.b.canonical, right.log_power);
}

template <typename Scalar>
LocalSolution<Scalar> with_selected_component(
    const LocalSolution<Scalar>& input, std::uint32_t component) {
  if (component >= input.dimension)
    throw std::out_of_range("selected local component is outside dimension");
  LocalSolution<Scalar> output;
  output.chart = input.chart;
  output.epsilon = input.epsilon;
  output.taylor_complete_max = input.taylor_complete_max;
  output.dimension = 1;
  output.prescriptions = input.prescriptions;
  output.checkpoint_identity = input.checkpoint_identity + ":component:" +
      std::to_string(component);
  output.sectors.reserve(input.sectors.size());
  for (const auto& sector : input.sectors) {
    LocalSector<Scalar> selected;
    selected.a = sector.a;
    selected.b = sector.b;
    selected.log_power = sector.log_power;
    selected.coefficients.assign(output.sector_size(),
                                 ScalarTraits<Scalar>::zero());
    for (std::size_t ei = 0; ei < input.epsilon.width(); ++ei)
      for (std::size_t n = 0; n < input.taylor_width(); ++n)
        selected.coefficients[flat_index(ei, n, 0, output.taylor_width(), 1)] =
            sector.coefficients[flat_index(
                ei, n, component, input.taylor_width(), input.dimension)];
    output.sectors.push_back(std::move(selected));
  }
  return output;
}

template <typename Scalar>
LocalSolution<Scalar> embedded_components(
    const LocalSolution<Scalar>& input,
    const std::vector<std::uint32_t>& components,
    std::uint32_t dimension) {
  if (input.dimension == 0 || input.dimension != components.size() ||
      dimension == 0)
    throw std::invalid_argument("invalid local component embedding");
  std::vector<std::uint8_t> seen(dimension, 0);
  for (const auto component : components) {
    if (component >= dimension || seen[component])
      throw std::invalid_argument("invalid local component embedding");
    seen[component] = 1;
  }
  LocalSolution<Scalar> output;
  output.chart = input.chart;
  output.epsilon = input.epsilon;
  output.taylor_complete_max = input.taylor_complete_max;
  output.dimension = dimension;
  output.prescriptions = input.prescriptions;
  output.checkpoint_identity = input.checkpoint_identity + ":embedded-block:" +
      std::to_string(input.dimension) + ":" + std::to_string(dimension);
  for (const auto component : components)
    output.checkpoint_identity += ":" + std::to_string(component);
  output.sectors.reserve(input.sectors.size());
  for (const auto& sector : input.sectors) {
    LocalSector<Scalar> embedded;
    embedded.a = sector.a;
    embedded.b = sector.b;
    embedded.log_power = sector.log_power;
    embedded.coefficients.assign(output.sector_size(),
                                 ScalarTraits<Scalar>::zero());
    for (std::size_t ei = 0; ei < input.epsilon.width(); ++ei)
      for (std::size_t n = 0; n < input.taylor_width(); ++n)
        for (std::uint32_t local = 0; local < input.dimension; ++local)
          embedded.coefficients[flat_index(
              ei, n, components[local], output.taylor_width(), dimension)] =
              sector.coefficients[flat_index(
                  ei, n, local, input.taylor_width(), input.dimension)];
    output.sectors.push_back(std::move(embedded));
  }
  return output;
}

template <typename Scalar>
LocalSolution<Scalar> embedded_component(const LocalSolution<Scalar>& input,
                                         std::uint32_t component,
                                         std::uint32_t dimension) {
  auto output = embedded_components(
      input, std::vector<std::uint32_t>{component}, dimension);
  output.checkpoint_identity = input.checkpoint_identity + ":embedded:" +
      std::to_string(component) + ":" + std::to_string(dimension);
  return output;
}

}  // namespace local_algebra_detail

// Merge only exactly identical tags.  Integer-spaced tower compaction is a
// representation optimization, not a semantic requirement, and is deferred
// until the native exact field can prove the integer difference and the
// finite Taylor slab can prove that shifting discards no nonzero tail.
template <typename Scalar>
LocalSolution<Scalar> canonicalize_identical_local_sectors(
    LocalSolution<Scalar> input) {
  validate_local_solution(input, false);
  std::vector<LocalSector<Scalar>> merged;
  for (auto& sector : input.sectors) {
    const auto found = std::find_if(merged.begin(), merged.end(),
        [&](const auto& candidate) {
          return local_algebra_detail::same_sector_tag(candidate, sector);
        });
    if (found == merged.end()) {
      merged.push_back(std::move(sector));
    } else {
      if (found->coefficients.size() != sector.coefficients.size())
        throw std::invalid_argument("identical local tags have unequal slabs");
      for (std::size_t i = 0; i < found->coefficients.size(); ++i)
        found->coefficients[i] += sector.coefficients[i];
    }
  }
  std::stable_sort(merged.begin(), merged.end(),
      local_algebra_detail::sector_less<Scalar>);
  input.sectors = std::move(merged);
  validate_local_solution(input, false);
  return input;
}

template <typename Scalar>
LocalSolution<Scalar> multiply_prepared_rational(
    const LocalSolution<Scalar>& input,
    const PreparedRationalTaylorMultiplier<Scalar>& multiplier,
    std::string checkpoint_identity = {}) {
  validate_local_solution(input, false);
  if (!input.error.empty())
    throw std::invalid_argument(
        "native rational multiplication needs explicit error-envelope propagation");
  const auto epsilon_width = input.epsilon.width();
  const auto taylor_width = input.taylor_width();
  if (multiplier.kernels.size() < epsilon_width)
    throw std::invalid_argument(
        "prepared rational multiplier has too few epsilon kernels");
  for (std::size_t j = 0; j < epsilon_width; ++j)
    if (multiplier.kernels[j].size() < taylor_width)
      throw std::invalid_argument(
          "prepared rational multiplier has too few Taylor coefficients");

  LocalSolution<Scalar> output;
  output.chart = input.chart;
  output.epsilon = {
      local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(input.epsilon.min_power) +
              multiplier.epsilon_shift,
          "rational-product epsilon minimum"),
      local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(input.epsilon.complete_max) +
              multiplier.epsilon_shift,
          "rational-product epsilon maximum")};
  output.taylor_complete_max = input.taylor_complete_max;
  output.dimension = input.dimension;
  output.prescriptions = input.prescriptions;
  output.checkpoint_identity = checkpoint_identity.empty()
      ? input.checkpoint_identity + ":mul:" + multiplier.exact_identity
      : std::move(checkpoint_identity);
  output.sectors.reserve(input.sectors.size());

  for (const auto& sector : input.sectors) {
    LocalSector<Scalar> product;
    product.a = local_algebra_detail::subtract_nonnegative_integer(
        sector.a, multiplier.center_pole_order);
    product.b = sector.b;
    product.log_power = sector.log_power;
    product.coefficients.assign(output.sector_size(),
                                ScalarTraits<Scalar>::zero());
    for (std::size_t out_ei = 0; out_ei < epsilon_width; ++out_ei) {
      for (std::size_t kernel_ei = 0; kernel_ei <= out_ei; ++kernel_ei) {
        const auto input_ei = out_ei - kernel_ei;
        const auto& kernel = multiplier.kernels[kernel_ei];
        for (std::size_t n = 0; n < taylor_width; ++n) {
          for (std::size_t m = 0; m <= n; ++m) {
            if (ScalarTraits<Scalar>::is_zero(kernel[m])) continue;
            for (std::uint32_t component = 0;
                 component < input.dimension; ++component) {
              product.coefficients[local_algebra_detail::flat_index(
                  out_ei, n, component, taylor_width, input.dimension)] +=
                  kernel[m] * sector.coefficients[
                      local_algebra_detail::flat_index(
                          input_ei, n - m, component, taylor_width,
                          input.dimension)];
            }
          }
        }
      }
    }
    output.sectors.push_back(std::move(product));
  }
  return canonicalize_identical_local_sectors(std::move(output));
}

// Restrict a finite local slab to a target recurrence frame.  Discarding
// known upper coefficients is safe because the target never consumes them.
// Discarding lower coefficients is only safe when every discarded exact
// coefficient is zero: nonzero lower content means the caller did not retain
// enough Laurent halo for the signed multiplier shift.
template <typename Scalar>
LocalSolution<Scalar> restrict_local_epsilon_frame_strict_lower(
    const LocalSolution<Scalar>& input, std::int32_t target_min,
    std::int32_t target_complete_max,
    std::string checkpoint_identity = {}) {
  validate_local_solution(input, false);
  if (!input.error.empty())
    throw std::invalid_argument(
        "native epsilon-frame restriction needs explicit error-envelope propagation");
  if (target_complete_max < target_min)
    throw std::invalid_argument("empty target epsilon frame");

  if (input.epsilon.min_power < target_min) {
    const auto discarded_max = std::min<std::int32_t>(
        input.epsilon.complete_max,
        local_algebra_detail::checked_i32(
            static_cast<std::int64_t>(target_min) - 1,
            "lower epsilon-frame boundary"));
    for (const auto& sector : input.sectors)
      for (std::int64_t power = input.epsilon.min_power;
           power <= discarded_max; ++power) {
        const auto epsilon_index = static_cast<std::size_t>(
            power - input.epsilon.min_power);
        for (std::size_t n = 0; n < input.taylor_width(); ++n)
          for (std::uint32_t component = 0;
               component < input.dimension; ++component)
            if (!ScalarTraits<Scalar>::is_zero(
                    sector.coefficients[local_algebra_detail::flat_index(
                        epsilon_index, n, component, input.taylor_width(),
                        input.dimension)]))
              throw std::invalid_argument(
                  "signed epsilon shift has nonzero content below the target frame; insufficient lower halo");
      }
  }

  if (input.epsilon.complete_max < target_min)
    throw std::invalid_argument(
        "signed epsilon shift leaves no complete coefficient in the target frame");

  const bool structurally_zero_on_target =
      input.epsilon.min_power > target_complete_max;
  const EpsilonWindow output_window = structurally_zero_on_target
      ? EpsilonWindow{target_min, target_complete_max}
      : EpsilonWindow{
            std::max(input.epsilon.min_power, target_min),
            std::min(input.epsilon.complete_max, target_complete_max)};

  LocalSolution<Scalar> output;
  output.chart = input.chart;
  output.epsilon = output_window;
  output.taylor_complete_max = input.taylor_complete_max;
  output.dimension = input.dimension;
  output.prescriptions = input.prescriptions;
  output.checkpoint_identity = checkpoint_identity.empty()
      ? input.checkpoint_identity + ":epsilon-frame"
      : std::move(checkpoint_identity);
  output.sectors.reserve(input.sectors.size());
  for (const auto& sector : input.sectors) {
    LocalSector<Scalar> restricted;
    restricted.a = sector.a;
    restricted.b = sector.b;
    restricted.log_power = sector.log_power;
    restricted.coefficients.assign(output.sector_size(),
                                    ScalarTraits<Scalar>::zero());
    if (!structurally_zero_on_target) {
      for (std::int64_t power = output.epsilon.min_power;
           power <= output.epsilon.complete_max; ++power) {
        const auto input_epsilon = static_cast<std::size_t>(
            power - input.epsilon.min_power);
        const auto output_epsilon = static_cast<std::size_t>(
            power - output.epsilon.min_power);
        for (std::size_t n = 0; n < output.taylor_width(); ++n)
          for (std::uint32_t component = 0;
               component < output.dimension; ++component)
            restricted.coefficients[local_algebra_detail::flat_index(
                output_epsilon, n, component, output.taylor_width(),
                output.dimension)] = sector.coefficients[
                    local_algebra_detail::flat_index(
                        input_epsilon, n, component, input.taylor_width(),
                        input.dimension)];
      }
    }
    output.sectors.push_back(std::move(restricted));
  }
  validate_local_solution(output, false);
  return output;
}

template <typename Scalar>
LocalSolution<Scalar> combine_local_solutions(
    const std::vector<LocalSolution<Scalar>>& inputs,
    std::string checkpoint_identity = {}) {
  if (inputs.empty())
    throw std::invalid_argument("cannot combine an empty local-solution list");
  for (const auto& input : inputs) {
    validate_local_solution(input, false);
    if (!input.error.empty())
      throw std::invalid_argument(
          "native local combination needs explicit error-envelope propagation");
    local_algebra_detail::require_same_local_space(inputs.front(), input);
  }

  const auto min_power = std::min_element(inputs.begin(), inputs.end(),
      [](const auto& left, const auto& right) {
        return left.epsilon.min_power < right.epsilon.min_power;
      })->epsilon.min_power;
  const auto complete_max = std::min_element(inputs.begin(), inputs.end(),
      [](const auto& left, const auto& right) {
        return left.epsilon.complete_max < right.epsilon.complete_max;
      })->epsilon.complete_max;
  if (complete_max < min_power)
    throw std::invalid_argument("combined local solutions have an empty epsilon window");
  const auto taylor_complete_max = std::min_element(
      inputs.begin(), inputs.end(), [](const auto& left, const auto& right) {
        return left.taylor_complete_max < right.taylor_complete_max;
      })->taylor_complete_max;

  LocalSolution<Scalar> output;
  output.chart = inputs.front().chart;
  output.epsilon = {min_power, complete_max};
  output.taylor_complete_max = taylor_complete_max;
  output.dimension = inputs.front().dimension;
  output.prescriptions = inputs.front().prescriptions;
  output.checkpoint_identity = checkpoint_identity.empty()
      ? inputs.front().checkpoint_identity + ":combined"
      : std::move(checkpoint_identity);

  const auto out_taylor_width = output.taylor_width();
  for (const auto& input : inputs) {
    for (const auto& sector : input.sectors) {
      LocalSector<Scalar> aligned;
      aligned.a = sector.a;
      aligned.b = sector.b;
      aligned.log_power = sector.log_power;
      aligned.coefficients.assign(output.sector_size(),
                                  ScalarTraits<Scalar>::zero());
      const auto copy_min = std::max(min_power, input.epsilon.min_power);
      const auto copy_max = std::min(complete_max, input.epsilon.complete_max);
      for (std::int32_t power = copy_min; power <= copy_max; ++power) {
        const auto input_ei = static_cast<std::size_t>(
            power - input.epsilon.min_power);
        const auto output_ei = static_cast<std::size_t>(power - min_power);
        for (std::size_t n = 0; n < out_taylor_width; ++n)
          for (std::uint32_t component = 0;
               component < output.dimension; ++component)
            aligned.coefficients[local_algebra_detail::flat_index(
                output_ei, n, component, out_taylor_width,
                output.dimension)] = sector.coefficients[
                    local_algebra_detail::flat_index(
                        input_ei, n, component, input.taylor_width(),
                        input.dimension)];
      }
      output.sectors.push_back(std::move(aligned));
    }
  }
  return canonicalize_identical_local_sectors(std::move(output));
}

// Materialize one retained matching state in its receiving chart without
// exporting either the basis slabs or the Laurent weights.  Each basis column
// is a full local vector S_j(t,eps), while weights[j] is a finite Laurent
// frame independent of t.  The complete upper edge is the first edge at which
// an unseen coefficient of either factor could contribute:
//
//   CompleteMax(S_j w_j) =
//     min(CompleteMax(S_j) + Min(w_j), Min(S_j) + CompleteMax(w_j)).
//
// Lower frame edges are exact structural bounds in both representations.
// No zero padding above a complete edge, midpoint zero test, or coefficient
// export is involved.
template <typename Scalar>
LocalSolution<Scalar> materialize_local_basis_weights(
    const std::vector<const LocalSolution<Scalar>*>& basis,
    const FiniteLaurentVector<Scalar>& weights,
    std::string checkpoint_identity) {
  if (basis.empty() || basis.size() != weights.size())
    throw std::invalid_argument(
        "local basis materialization requires one weight per nonempty basis column");
  if (checkpoint_identity.empty())
    throw std::invalid_argument(
        "local basis materialization checkpoint identity is empty");
  for (const auto* column : basis) {
    if (column == nullptr)
      throw std::invalid_argument(
          "local basis materialization received a null basis column");
    validate_local_solution(*column, false);
    if (!column->error.empty())
      throw std::invalid_argument(
          "local basis materialization cannot discard an error envelope");
    local_algebra_detail::require_same_local_space(*basis.front(), *column);
  }

  std::vector<EpsilonWindow> product_windows;
  product_windows.reserve(basis.size());
  for (std::size_t column = 0; column < basis.size(); ++column) {
    const auto minimum = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(basis[column]->epsilon.min_power) +
            weights[column].min_power(),
        "materialized local epsilon minimum");
    const auto basis_complete = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(basis[column]->epsilon.complete_max) +
            weights[column].min_power(),
        "materialized local basis-complete edge");
    const auto weight_complete = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(basis[column]->epsilon.min_power) +
            weights[column].complete_max(),
        "materialized local weight-complete edge");
    product_windows.push_back(
        {minimum, std::min(basis_complete, weight_complete)});
  }
  const auto minimum = std::min_element(
      product_windows.begin(), product_windows.end(),
      [](const auto& left, const auto& right) {
        return left.min_power < right.min_power;
      })->min_power;
  const auto complete_max = std::min_element(
      product_windows.begin(), product_windows.end(),
      [](const auto& left, const auto& right) {
        return left.complete_max < right.complete_max;
      })->complete_max;
  if (complete_max < minimum)
    throw std::invalid_argument(
        "materialized local basis has no common complete epsilon window");

  LocalSolution<Scalar> output;
  output.chart = basis.front()->chart;
  output.epsilon = {minimum, complete_max};
  output.taylor_complete_max = (*std::min_element(
      basis.begin(), basis.end(), [](const auto* left, const auto* right) {
        return left->taylor_complete_max < right->taylor_complete_max;
      }))->taylor_complete_max;
  output.dimension = basis.front()->dimension;
  output.prescriptions = basis.front()->prescriptions;
  output.checkpoint_identity = std::move(checkpoint_identity);

  const auto output_taylor_width = output.taylor_width();
  for (std::size_t column = 0; column < basis.size(); ++column) {
    const auto& source = *basis[column];
    const auto& weight = weights[column];
    for (const auto& sector : source.sectors) {
      LocalSector<Scalar> materialized;
      materialized.a = sector.a;
      materialized.b = sector.b;
      materialized.log_power = sector.log_power;
      materialized.coefficients.assign(
          output.sector_size(), ScalarTraits<Scalar>::zero());
      for (std::int64_t power = output.epsilon.min_power;
           power <= output.epsilon.complete_max; ++power) {
        const auto output_epsilon = static_cast<std::size_t>(
            power - output.epsilon.min_power);
        for (std::int64_t weight_power = weight.min_power();
             weight_power <= weight.complete_max(); ++weight_power) {
          const auto source_power = power - weight_power;
          if (source_power < source.epsilon.min_power ||
              source_power > source.epsilon.complete_max)
            continue;
          const auto source_epsilon = static_cast<std::size_t>(
              source_power - source.epsilon.min_power);
          const auto& scalar_weight = weight.coefficient(
              static_cast<std::int32_t>(weight_power));
          if (ScalarTraits<Scalar>::is_zero(scalar_weight)) continue;
          for (std::size_t n = 0; n < output_taylor_width; ++n)
            for (std::uint32_t component = 0;
                 component < output.dimension; ++component)
              materialized.coefficients[local_algebra_detail::flat_index(
                  output_epsilon, n, component, output_taylor_width,
                  output.dimension)] += scalar_weight * sector.coefficients[
                      local_algebra_detail::flat_index(
                          source_epsilon, n, component,
                          source.taylor_width(), source.dimension)];
        }
      }
      output.sectors.push_back(std::move(materialized));
    }
  }
  return canonicalize_identical_local_sectors(std::move(output));
}

template <typename Scalar>
std::optional<LocalSolution<Scalar>> apply_prepared_sparse_local_matrix(
    const PreparedSparseLocalMultiplierMatrix<Scalar>& matrix,
    const LocalSolution<Scalar>& input,
    std::string checkpoint_identity = {}) {
  validate_local_solution(input, false);
  if (!input.error.empty())
    throw std::invalid_argument(
        "native local matrix application needs explicit error-envelope propagation");
  if (matrix.rows == 0 || matrix.columns != input.dimension)
    throw std::invalid_argument("prepared local matrix dimensions disagree");
  std::vector<LocalSolution<Scalar>> terms;
  terms.reserve(matrix.entries.size());
  for (const auto& entry : matrix.entries) {
    if (entry.row >= matrix.rows || entry.column >= matrix.columns)
      throw std::invalid_argument("prepared local matrix entry is out of range");
    if (entry.multiplier.structurally_zero()) continue;
    auto selected = local_algebra_detail::with_selected_component(
        input, entry.column);
    auto product = multiply_prepared_rational(
        selected, entry.multiplier,
        input.checkpoint_identity + ":matrix-entry:" +
            std::to_string(entry.row) + ":" +
            std::to_string(entry.column));
    // Keep a zero finite slab from an exact nonzero matrix entry: in the
    // Wolfram algebra it still constrains the honest intersection window
    // when another entry of the same matrix product is nonzero.  Only the
    // fully combined matrix result may be discarded as a structural source.
    terms.push_back(local_algebra_detail::embedded_component(
        product, entry.row, matrix.rows));
  }
  if (terms.empty()) return std::nullopt;
  return combine_local_solutions(terms,
      checkpoint_identity.empty()
          ? input.checkpoint_identity + ":matrix:" + matrix.exact_identity
          : std::move(checkpoint_identity));
}

}  // namespace diffexp2
