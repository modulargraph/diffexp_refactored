#pragma once

#include "diffexp2/local_algebra.hpp"
#include "diffexp2/physical_ode.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace diffexp2 {

// Finite-epsilon Taylor problem for the backward observable equation
//
//             q(t,eps) theta(lambda) + C(t,eps)^T lambda = b(t,eps).
//
// The endpoint condition lambda(0)=0 is structural: coefficients start at
// t^1.  This is the composed alternative to separately truncating a terminal
// basis A, a line row J, and the ill-conditioned connection solve A^-1 u.
// Each q/C/forcing entry is an honest finite Laurent frame.  A resonant layer
// is rejected by the certified Laurent factorization; logarithmic completion
// must be supplied by the later singular extension, never guessed here.
struct BackwardAdjointTaylorProblem {
  std::uint32_t dimension = 0;
  std::uint32_t taylor_complete_max = 0;
  std::int32_t required_epsilon_complete_max = 0;
  std::vector<EpsilonFrame<ComplexBall>> q_lags;
  std::vector<FiniteLaurentMatrix<ComplexBall>> c_transpose_lags;
  // forcing[n-1] is the coefficient of t^n.  A short forcing is completed by
  // structural zeros, while coefficients beyond taylor_complete_max are not
  // consumed by this finite recurrence.
  std::vector<FiniteLaurentVector<ComplexBall>> forcing;
};

struct BackwardAdjointTaylorResult {
  std::uint32_t dimension = 0;
  std::uint32_t taylor_complete_max = 0;
  std::int32_t common_epsilon_complete_max = 0;
  // coefficients[n-1][component] is lambda_{component,n}, n >= 1.
  std::vector<FiniteLaurentVector<ComplexBall>> coefficients;

  // This is only a planning signal.  It is not a tail certificate and must
  // never be used to publish a full-local result.
  std::optional<double> last_coefficient_ratio_upper;
};

struct BackwardAdjointTailCertificate {
  Magnitude absolute_vector_tail_upper = Magnitude::zero();
  Magnitude coefficient_majorant_upper = Magnitude::zero();
  Magnitude recurrence_contraction_upper = Magnitude::zero();
  Magnitude evaluation_ratio_upper = Magnitude::zero();
  std::uint32_t certified_after_taylor_order = 0;
};

namespace adjoint_observable_detail {

inline EpsilonFrame<ComplexBall> structural_zero_frame(
    std::int32_t complete_max) {
  return EpsilonFrame<ComplexBall>::zero(complete_max);
}

inline FiniteLaurentVector<ComplexBall> structural_zero_vector(
    std::uint32_t dimension, std::int32_t complete_max) {
  FiniteLaurentVector<ComplexBall> result;
  result.reserve(dimension);
  for (std::uint32_t component = 0; component < dimension; ++component)
    result.push_back(structural_zero_frame(complete_max));
  return result;
}

inline void require_matrix_shape(
    const FiniteLaurentMatrix<ComplexBall>& matrix,
    std::uint32_t dimension, const char* label) {
  if (matrix.size() != dimension)
    throw std::invalid_argument(std::string(label) +
                                " has the wrong row count");
  for (const auto& row : matrix)
    if (row.size() != dimension)
      throw std::invalid_argument(std::string(label) +
                                  " has the wrong column count");
}

inline EpsilonFrame<ComplexBall> scaled_integer(
    const EpsilonFrame<ComplexBall>& value, std::int64_t factor) {
  return value.scaled(ComplexBall::from_strings(std::to_string(factor)));
}

inline FiniteLaurentVector<ComplexBall> apply_lag_operator(
    const EpsilonFrame<ComplexBall>& q,
    const FiniteLaurentMatrix<ComplexBall>* c_transpose,
    std::int64_t theta_factor,
    const FiniteLaurentVector<ComplexBall>& input,
    std::int32_t structural_complete_max) {
  const auto dimension = static_cast<std::uint32_t>(input.size());
  auto output = structural_zero_vector(dimension, structural_complete_max);
  for (std::uint32_t row = 0; row < dimension; ++row) {
    if (theta_factor != 0)
      output[row] = output[row] +
          scaled_integer(q, theta_factor) * input[row];
    if (c_transpose == nullptr) continue;
    for (std::uint32_t column = 0; column < dimension; ++column)
      output[row] = output[row] +
          (*c_transpose)[row][column] * input[column];
  }
  return output;
}

inline double vector_upper_norm(
    const FiniteLaurentVector<ComplexBall>& value) {
  double result = 0.0;
  for (const auto& component : value)
    for (std::int64_t raw_power = component.min_power();
         raw_power <= component.complete_max(); ++raw_power)
      result = std::max(
          result,
          Magnitude::upper_abs(component.coefficient(
              static_cast<std::int32_t>(raw_power))).approximate_upper());
  return result;
}

template <typename Scalar>
EpsilonFrame<ComplexBall> expand_epsilon_rational(
    const ExactEpsilonRational<Scalar>& value,
    std::int32_t complete_max, const std::string& context) {
  physical_ode_detail::validate_rational(
      value, "backward adjoint epsilon-rational coefficient");
  if (value.zero || value.valuation > complete_max)
    return structural_zero_frame(complete_max);
  const auto width = EpsilonWindow{
      value.valuation, complete_max}.width();
  std::vector<ComplexBall> coefficients(width, ComplexBall(0));
  const auto denominator0 = local_detail::to_ball(value.denominator.front());
  if (denominator0.contains_zero())
    throw std::domain_error(
        context + ": epsilon-rational denominator contains zero");
  for (std::size_t offset = 0; offset < width; ++offset) {
    ComplexBall coefficient(0);
    if (offset < value.numerator.size())
      coefficient = local_detail::to_ball(value.numerator[offset]);
    for (std::size_t degree = 1;
         degree < value.denominator.size() && degree <= offset; ++degree)
      coefficient -= local_detail::to_ball(value.denominator[degree]) *
                     coefficients[offset - degree];
    coefficients[offset] = coefficient / denominator0;
  }
  return EpsilonFrame<ComplexBall>(value.valuation,
                                   std::move(coefficients));
}

inline bool material_frame(const EpsilonFrame<ComplexBall>& value) {
  return std::any_of(
      value.coefficients().begin(), value.coefficients().end(),
      [](const auto& coefficient) { return !coefficient.is_zero(); });
}

inline Magnitude frame_l1_upper(
    const EpsilonFrame<ComplexBall>& value) {
  auto result = Magnitude::zero();
  for (const auto& coefficient : value.coefficients())
    result += Magnitude::upper_abs(coefficient);
  return result;
}

inline Magnitude vector_infinity_upper(
    const FiniteLaurentVector<ComplexBall>& value) {
  auto result = Magnitude::zero();
  for (const auto& component : value)
    for (const auto& coefficient : component.coefficients())
      result = Magnitude::maximum(
          result, Magnitude::upper_abs(coefficient));
  return result;
}

inline Magnitude matrix_infinity_upper(
    const FiniteLaurentMatrix<ComplexBall>& matrix) {
  auto result = Magnitude::zero();
  for (const auto& row : matrix) {
    auto row_sum = Magnitude::zero();
    for (const auto& entry : row) row_sum += frame_l1_upper(entry);
    result = Magnitude::maximum(result, row_sum);
  }
  return result;
}

inline Magnitude polynomial_disk_upper(
    const std::vector<ComplexBall>& coefficients,
    const Magnitude& radius_upper) {
  auto result = Magnitude::zero();
  auto power = Magnitude::one();
  for (const auto& coefficient : coefficients) {
    result += Magnitude::upper_abs(coefficient) * power;
    power = power * radius_upper;
  }
  return result;
}

inline Magnitude rational_disk_upper(
    const PreparedRationalAnalyticCoefficient<ComplexBall>& rational,
    const Magnitude& radius_upper, const std::string& context) {
  if (rational.numerator.empty() || rational.denominator.empty())
    throw std::invalid_argument(
        context + ": analytic rational coefficient is empty");
  const auto numerator =
      polynomial_disk_upper(rational.numerator, radius_upper);
  auto denominator_tail = Magnitude::zero();
  auto power = radius_upper;
  for (std::size_t degree = 1; degree < rational.denominator.size();
       ++degree) {
    denominator_tail +=
        Magnitude::upper_abs(rational.denominator[degree]) * power;
    power = power * radius_upper;
  }
  const auto denominator_lower = Magnitude::positive_difference_lower(
      Magnitude::lower_abs(rational.denominator.front()),
      denominator_tail);
  if (denominator_lower.is_zero())
    throw std::domain_error(
        context +
        ": witness disk does not prove the analytic denominator nonzero");
  return numerator / denominator_lower;
}

}  // namespace adjoint_observable_detail

template <typename Scalar>
BackwardAdjointTaylorProblem prepare_backward_adjoint_taylor_problem(
    const PreparedPhysicalClearedODE<Scalar>& ode,
    const PreparedSparseLocalMultiplierMatrix<ComplexBall>& row,
    std::uint32_t taylor_complete_max,
    std::int32_t coefficient_epsilon_complete_max,
    std::int32_t required_epsilon_complete_max,
    const ComplexBall& oriented_physical_jacobian,
    const std::string& context = "prepare backward Fuchsian adjoint") {
  using namespace adjoint_observable_detail;
  physical_ode_detail::validate_ode(ode);
  if (taylor_complete_max == 0 || row.rows != 1 ||
      row.columns != ode.dimension || row.entries.empty() ||
      coefficient_epsilon_complete_max < required_epsilon_complete_max)
    throw std::invalid_argument(
        context + ": inconsistent ODE, row, Taylor, or epsilon contract");

  BackwardAdjointTaylorProblem problem;
  problem.dimension = ode.dimension;
  problem.taylor_complete_max = taylor_complete_max;
  problem.required_epsilon_complete_max =
      required_epsilon_complete_max;
  problem.q_lags.reserve(
      std::min<std::size_t>(ode.q_lags.size(), taylor_complete_max + 1));
  for (std::size_t lag = 0;
       lag < ode.q_lags.size() && lag <= taylor_complete_max; ++lag)
    problem.q_lags.push_back(expand_epsilon_rational(
        ode.q_lags[lag], coefficient_epsilon_complete_max,
        context + ": q lag " + std::to_string(lag)));

  const auto c_lag_count = std::min<std::size_t>(
      ode.c_lags.size(), taylor_complete_max + 1);
  problem.c_transpose_lags.reserve(c_lag_count);
  for (std::size_t lag = 0; lag < c_lag_count; ++lag) {
    FiniteLaurentMatrix<ComplexBall> transposed(
        ode.dimension, FiniteLaurentVector<ComplexBall>());
    for (auto& target_row : transposed)
      for (std::uint32_t column = 0; column < ode.dimension; ++column)
        target_row.push_back(
            structural_zero_frame(coefficient_epsilon_complete_max));
    for (const auto& entry : ode.c_lags[lag])
      transposed[entry.column][entry.row] = expand_epsilon_rational(
          entry.value, coefficient_epsilon_complete_max,
          context + ": C lag " + std::to_string(lag));
    problem.c_transpose_lags.push_back(std::move(transposed));
  }

  problem.forcing.assign(
      taylor_complete_max,
      structural_zero_vector(ode.dimension,
                             coefficient_epsilon_complete_max));
  for (const auto& entry : row.entries) {
    if (entry.row != 0 || entry.column >= ode.dimension ||
        entry.multiplier.structurally_zero())
      throw std::invalid_argument(
          context + ": prepared row contains an invalid active entry");
    const auto& multiplier = entry.multiplier;
    if (multiplier.kernels.empty())
      throw std::invalid_argument(
          context + ": prepared row entry has no Taylor kernels");
    const auto taylor_width = multiplier.kernels.front().size();
    if (taylor_width == 0 ||
        std::any_of(multiplier.kernels.begin(), multiplier.kernels.end(),
                    [&](const auto& kernel) {
                      return kernel.size() != taylor_width;
                    }))
      throw std::invalid_argument(
          context + ": prepared row entry has a ragged Taylor kernel");

    for (std::size_t row_taylor = 0; row_taylor < taylor_width;
         ++row_taylor) {
      std::vector<ComplexBall> row_epsilon;
      row_epsilon.reserve(multiplier.kernels.size());
      for (const auto& epsilon_kernel : multiplier.kernels)
        row_epsilon.push_back(epsilon_kernel[row_taylor]);
      EpsilonFrame<ComplexBall> row_coefficient(
          multiplier.epsilon_shift, std::move(row_epsilon));
      if (!material_frame(row_coefficient)) continue;
      const auto row_power =
          static_cast<std::int64_t>(row_taylor) -
          multiplier.center_pole_order;
      for (std::size_t q_lag = 0; q_lag < problem.q_lags.size(); ++q_lag) {
        if (!material_frame(problem.q_lags[q_lag])) continue;
        const auto forcing_power = row_power +
            static_cast<std::int64_t>(q_lag) + 1;
        auto forcing_coefficient =
            problem.q_lags[q_lag] * row_coefficient;
        forcing_coefficient = forcing_coefficient.scaled(
            -oriented_physical_jacobian);
        if (!material_frame(forcing_coefficient)) continue;
        if (forcing_power <= 0)
          throw std::domain_error(
              context +
              ": center-pole forcing reaches t^0 or below; a Laurent/log "
              "backward Fuchsian completion is required");
        if (forcing_power > taylor_complete_max) continue;
        auto& target = problem.forcing[
            static_cast<std::size_t>(forcing_power - 1)][entry.column];
        target = target + forcing_coefficient;
      }
    }
  }
  return problem;
}

inline Magnitude backward_adjoint_forcing_cauchy_numerator_upper(
    const BackwardAdjointTaylorProblem& problem,
    const PreparedSparseLocalMultiplierMatrix<ComplexBall>& row,
    const ComplexBall& oriented_physical_jacobian,
    const ComplexBall& witness_radius,
    const std::string& context =
        "backward Fuchsian adjoint forcing Cauchy bound") {
  using namespace adjoint_observable_detail;
  const auto radius_lower = Magnitude::lower_abs(witness_radius);
  const auto radius_upper = Magnitude::upper_abs(witness_radius);
  if (radius_lower.is_zero() || !radius_upper.is_finite() ||
      row.rows != 1 || row.columns != problem.dimension)
    throw std::invalid_argument(
        context + ": invalid witness radius or row dimension");

  auto q_disk_upper = Magnitude::zero();
  auto t_power = Magnitude::one();
  for (const auto& q : problem.q_lags) {
    q_disk_upper += frame_l1_upper(q) * t_power;
    t_power = t_power * radius_upper;
  }

  auto row_disk_upper = Magnitude::zero();
  for (const auto& entry : row.entries) {
    const auto& multiplier = entry.multiplier;
    if (multiplier.center_pole_order != 0)
      throw std::domain_error(
          context +
          ": ordinary composed-tail theorem does not admit a center pole");
    if (!multiplier.analytic_coefficients.has_value() ||
        multiplier.analytic_coefficients->size() !=
            multiplier.kernels.size())
      throw std::domain_error(
          context +
          ": complete analytic row coefficients are unavailable");
    auto component_upper = Magnitude::zero();
    for (std::size_t epsilon = 0;
         epsilon < multiplier.analytic_coefficients->size(); ++epsilon)
      component_upper += rational_disk_upper(
          (*multiplier.analytic_coefficients)[epsilon], radius_upper,
          context + ": row component " +
              std::to_string(entry.column) + ": epsilon offset " +
              std::to_string(epsilon));
    row_disk_upper = Magnitude::maximum(
        row_disk_upper, component_upper);
  }
  if (row_disk_upper.is_zero()) return Magnitude::zero();
  return Magnitude::upper_abs(oriented_physical_jacobian) *
         radius_upper * q_disk_upper * row_disk_upper;
}

inline BackwardAdjointTailCertificate
certify_backward_adjoint_taylor_tail(
    const BackwardAdjointTaylorProblem& problem,
    const BackwardAdjointTaylorResult& solution,
    const ComplexBall& evaluation_point,
    const ComplexBall& witness_radius,
    const Magnitude& forcing_cauchy_numerator_upper,
    const std::string& context =
        "backward Fuchsian adjoint Taylor-tail certificate") {
  using namespace adjoint_observable_detail;
  if (solution.dimension != problem.dimension ||
      solution.taylor_complete_max != problem.taylor_complete_max ||
      solution.coefficients.size() != problem.taylor_complete_max ||
      problem.q_lags.empty() || problem.c_transpose_lags.empty())
    throw std::invalid_argument(
        context + ": problem and Taylor solution disagree");
  const auto evaluation_upper = Magnitude::upper_abs(evaluation_point);
  const auto radius_lower = Magnitude::lower_abs(witness_radius);
  const auto radius_upper = Magnitude::upper_abs(witness_radius);
  if (radius_lower.is_zero() || !radius_upper.is_finite())
    throw std::invalid_argument(context + ": invalid witness radius");
  const auto evaluation_ratio = evaluation_upper / radius_lower;
  if (Magnitude::one() <= evaluation_ratio)
    throw std::domain_error(
        context + ": evaluation point is not strictly inside witness disk");

  const auto& q0 = problem.q_lags.front();
  if (q0.min_power() != 0 || q0.coefficient(0).contains_zero())
    throw std::domain_error(
        context + ": q(0,epsilon) is not a certified formal unit");
  for (const auto& q : problem.q_lags)
    if (q.min_power() < 0)
      throw std::domain_error(
          context +
          ": negative-epsilon q coupling is not closed in the finite causal stack");
  for (const auto& matrix : problem.c_transpose_lags)
    for (const auto& row : matrix)
      for (const auto& entry : row)
        if (entry.min_power() < 0)
          throw std::domain_error(
              context +
              ": negative-epsilon C coupling is not closed in the finite causal stack");
  std::vector<ComplexBall> one_coefficients(q0.window().width(),
                                            ComplexBall(0));
  one_coefficients.front() = ComplexBall(1);
  const auto q0_inverse = finite_laurent_quotient(
      EpsilonFrame<ComplexBall>(q0.window(), std::move(one_coefficients)),
      q0, context + ": q0 inverse");
  const auto q_inverse_upper = frame_l1_upper(q0_inverse);
  const auto c0_upper =
      matrix_infinity_upper(problem.c_transpose_lags.front());
  const auto indicial_upper = q_inverse_upper * c0_upper;

  const auto next_order =
      static_cast<ulong>(problem.taylor_complete_max) + 1;
  const auto next_order_magnitude = Magnitude::from_ui(next_order);
  const auto inverse_denominator_lower =
      Magnitude::positive_difference_lower(
          next_order_magnitude, indicial_upper);
  if (inverse_denominator_lower.is_zero())
    throw std::domain_error(
        context +
        ": retained Taylor order does not dominate the indicial norm");

  auto q_tail_sum = Magnitude::zero();
  for (std::size_t lag = 1; lag < problem.q_lags.size(); ++lag)
    q_tail_sum += frame_l1_upper(problem.q_lags[lag]) *
                  radius_upper.power_upper(static_cast<ulong>(lag));
  auto c_tail_sum = Magnitude::zero();
  for (std::size_t lag = 1;
       lag < problem.c_transpose_lags.size(); ++lag)
    c_tail_sum += matrix_infinity_upper(
                      problem.c_transpose_lags[lag]) *
                  radius_upper.power_upper(static_cast<ulong>(lag));
  const auto recurrence_contraction =
      q_inverse_upper *
      (next_order_magnitude * q_tail_sum + c_tail_sum) /
      inverse_denominator_lower;
  if (Magnitude::one() <= recurrence_contraction)
    throw std::domain_error(
        context +
        ": recurrence majorant is not contractive on the witness disk");
  const auto contraction_gap_lower =
      Magnitude::positive_difference_lower(
          Magnitude::one(), recurrence_contraction);
  if (contraction_gap_lower.is_zero())
    throw std::domain_error(
        context + ": recurrence contraction has no positive gap");

  auto coefficient_majorant =
      q_inverse_upper * forcing_cauchy_numerator_upper /
      inverse_denominator_lower / contraction_gap_lower;
  const auto lag_memory = std::max(
      problem.q_lags.size(), problem.c_transpose_lags.size());
  const auto first_known = problem.taylor_complete_max + 1 > lag_memory
      ? problem.taylor_complete_max + 1 - lag_memory
      : 1;
  for (std::uint32_t order = first_known;
       order <= problem.taylor_complete_max; ++order) {
    const auto known = vector_infinity_upper(
        solution.coefficients.at(order - 1)) *
        radius_upper.power_upper(order);
    coefficient_majorant = Magnitude::maximum(
        coefficient_majorant, known);
  }

  const auto ratio_gap_lower = Magnitude::positive_difference_lower(
      Magnitude::one(), evaluation_ratio);
  if (ratio_gap_lower.is_zero())
    throw std::domain_error(
        context + ": evaluation/witness ratio has no positive gap");
  const auto tail = coefficient_majorant *
      evaluation_ratio.power_upper(next_order) / ratio_gap_lower;
  return {tail, coefficient_majorant, recurrence_contraction,
          evaluation_ratio, problem.taylor_complete_max};
}

inline BackwardAdjointTaylorResult solve_backward_adjoint_taylor(
    const BackwardAdjointTaylorProblem& problem,
    const std::string& context = "backward Fuchsian adjoint") {
  using namespace adjoint_observable_detail;
  if (problem.dimension == 0 || problem.taylor_complete_max == 0 ||
      problem.q_lags.empty() || problem.c_transpose_lags.empty())
    throw std::invalid_argument(
        context + ": incomplete dimension, Taylor order, q, or C payload");
  if (problem.forcing.size() > problem.taylor_complete_max)
    throw std::invalid_argument(
        context + ": forcing extends beyond the requested Taylor order");
  for (const auto& matrix : problem.c_transpose_lags)
    require_matrix_shape(matrix, problem.dimension,
                         "backward adjoint C lag");
  for (const auto& rhs : problem.forcing)
    if (rhs.size() != problem.dimension)
      throw std::invalid_argument(
          context + ": forcing vector has the wrong dimension");

  BackwardAdjointTaylorResult result;
  result.dimension = problem.dimension;
  result.taylor_complete_max = problem.taylor_complete_max;
  result.common_epsilon_complete_max =
      std::numeric_limits<std::int32_t>::max();
  result.coefficients.reserve(problem.taylor_complete_max);

  for (std::uint32_t n = 1; n <= problem.taylor_complete_max; ++n) {
    auto rhs = n <= problem.forcing.size()
        ? problem.forcing[n - 1]
        : structural_zero_vector(
              problem.dimension, problem.required_epsilon_complete_max);

    const auto maximum_lag = std::min<std::size_t>(
        n - 1,
        std::max(problem.q_lags.size(),
                 problem.c_transpose_lags.size()) - 1);
    for (std::size_t lag = 1; lag <= maximum_lag; ++lag) {
      const auto previous_n = n - static_cast<std::uint32_t>(lag);
      const auto& previous = result.coefficients.at(previous_n - 1);
      const auto& q = lag < problem.q_lags.size()
          ? problem.q_lags[lag]
          : structural_zero_frame(problem.required_epsilon_complete_max);
      const auto* c = lag < problem.c_transpose_lags.size()
          ? &problem.c_transpose_lags[lag]
          : nullptr;
      const auto contribution = apply_lag_operator(
          q, c, previous_n, previous,
          problem.required_epsilon_complete_max);
      for (std::uint32_t component = 0;
           component < problem.dimension; ++component)
        rhs[component] = rhs[component] - contribution[component];
    }

    FiniteLaurentMatrix<ComplexBall> layer(
        problem.dimension, FiniteLaurentVector<ComplexBall>());
    for (std::uint32_t row = 0; row < problem.dimension; ++row) {
      layer[row].reserve(problem.dimension);
      for (std::uint32_t column = 0; column < problem.dimension; ++column) {
        auto value = problem.c_transpose_lags.front()[row][column];
        if (row == column)
          value = value + scaled_integer(problem.q_lags.front(), n);
        layer[row].push_back(std::move(value));
      }
    }

    FiniteLaurentVector<ComplexBall> solved;
    try {
      auto factorization = factor_preconditioned_acb_finite_laurent_system(
          std::move(layer), context + ": Taylor layer " +
                                std::to_string(n));
      solved = solve_factorized_finite_laurent_system(
          factorization, std::move(rhs),
          context + ": Taylor coefficient " + std::to_string(n));
    } catch (const MatchingArithmeticError& error) {
      throw MatchingArithmeticError(
          error.code,
          context + ": Taylor layer " + std::to_string(n) +
              " is singular or epsilon-incomplete; a logarithmic/resonant "
              "Fuchsian completion is required; detail=" + error.what(),
          error.row, error.column, error.epsilon_power);
    }
    for (const auto& component : solved) {
      result.common_epsilon_complete_max = std::min(
          result.common_epsilon_complete_max, component.complete_max());
      if (component.complete_max() <
          problem.required_epsilon_complete_max)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            context + ": Taylor layer " + std::to_string(n) +
                " does not cover the required epsilon output",
            std::nullopt, std::nullopt,
            component.complete_max());
    }
    result.coefficients.push_back(std::move(solved));
  }

  if (result.coefficients.size() >= 2) {
    const auto previous = vector_upper_norm(
        result.coefficients[result.coefficients.size() - 2]);
    const auto last = vector_upper_norm(result.coefficients.back());
    if (previous > 0.0 && std::isfinite(previous) &&
        std::isfinite(last))
      result.last_coefficient_ratio_upper = last / previous;
  }
  return result;
}

inline FiniteLaurentVector<ComplexBall> evaluate_backward_adjoint_taylor(
    const BackwardAdjointTaylorResult& solution,
    const ComplexBall& point,
    const std::string& context = "backward Fuchsian adjoint evaluation") {
  using namespace adjoint_observable_detail;
  if (solution.dimension == 0 || solution.coefficients.empty() ||
      solution.coefficients.size() != solution.taylor_complete_max)
    throw std::invalid_argument(context + ": malformed Taylor solution");
  auto value = structural_zero_vector(
      solution.dimension, solution.common_epsilon_complete_max);
  for (std::size_t reverse = solution.coefficients.size(); reverse-- > 0;) {
    for (std::uint32_t component = 0;
         component < solution.dimension; ++component)
      value[component] = value[component].scaled(point) +
                         solution.coefficients[reverse][component];
  }
  for (auto& component : value) component = component.scaled(point);
  return value;
}

inline EpsilonFrame<ComplexBall> contract_backward_adjoint(
    const FiniteLaurentVector<ComplexBall>& adjoint,
    const FiniteLaurentVector<ComplexBall>& physical_value,
    const std::string& context = "backward Fuchsian adjoint contraction") {
  if (adjoint.empty() || adjoint.size() != physical_value.size())
    throw std::invalid_argument(context + ": vector dimensions disagree");
  std::optional<EpsilonFrame<ComplexBall>> result;
  for (std::size_t component = 0; component < adjoint.size(); ++component) {
    auto term = adjoint[component] * physical_value[component];
    result = result.has_value() ? *result + term : std::move(term);
  }
  if (!result.has_value())
    throw std::logic_error(context + ": empty contraction");
  return *result;
}

}  // namespace diffexp2
