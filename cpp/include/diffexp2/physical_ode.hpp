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

// Exact-in-epsilon coefficient retained from the Wolfram physical cleared
// equation.  Nonzero entries are normalized as
//
//                  eps^valuation P(eps) / Q(eps),
//
// with P(0) != 0 and Q(0) = 1.  Keeping this representation independent of
// the recurrence work frame is essential: a local may honestly expose a
// narrower epsilon window than the prepared operator, while causal division
// by the unit Q still determines every coefficient in that exposed window.
template <typename Scalar>
struct ExactEpsilonRational {
  bool zero = true;
  std::int32_t valuation = 0;
  std::vector<Scalar> numerator;
  std::vector<Scalar> denominator;
};

template <typename Scalar>
struct PhysicalODEMatrixEntry {
  std::uint32_t row = 0;
  std::uint32_t column = 0;
  ExactEpsilonRational<Scalar> value;
};

template <typename Scalar>
struct PreparedPhysicalClearedODE {
  std::uint32_t dimension = 0;
  std::vector<ExactEpsilonRational<Scalar>> q_lags;
  std::vector<std::vector<PhysicalODEMatrixEntry<Scalar>>> c_lags;
  std::string owner_signature_identity;
  std::string payload_identity;
  std::string exact_payload_record;
};

namespace physical_ode_detail {

template <typename Scalar>
void validate_rational(const ExactEpsilonRational<Scalar>& value,
                       const char* label) {
  if (value.zero) {
    if (!value.numerator.empty() || !value.denominator.empty())
      throw std::invalid_argument(std::string(label) +
                                  " zero entry carries polynomial data");
    return;
  }
  if (value.numerator.empty() || value.denominator.empty() ||
      ScalarTraits<Scalar>::is_zero(value.numerator.front()) ||
      ScalarTraits<Scalar>::is_zero(value.denominator.front()))
    throw std::invalid_argument(std::string(label) +
                                " has a non-unit causal representation");
  if constexpr (std::is_same_v<Scalar, Rational>) {
    if (!(value.denominator.front() == Rational(1)))
      throw std::invalid_argument(std::string(label) +
                                  " denominator constant is not exactly one");
  } else if constexpr (std::is_same_v<Scalar, ComplexBall>) {
    const ComplexBall one(1);
    if (value.numerator.front().contains_zero())
      throw std::invalid_argument(std::string(label) +
                                  " leading numerator does not prove its normalized valuation");
    if (value.denominator.front().contains_zero() ||
        !acb_equal(value.denominator.front().raw(), one.raw()))
      throw std::invalid_argument(std::string(label) +
                                  " denominator constant is not the exact singleton one");
  }
}

template <typename Scalar>
void validate_ode(const PreparedPhysicalClearedODE<Scalar>& ode) {
  if (ode.dimension == 0 || ode.q_lags.empty() || ode.c_lags.empty() ||
      ode.owner_signature_identity.empty() || ode.payload_identity.empty() ||
      ode.exact_payload_record.empty())
    throw std::invalid_argument(
        "physical cleared ODE lost its dimension, identity, or polynomial payload");
  bool nonzero_q = false;
  for (const auto& lag : ode.q_lags) {
    validate_rational(lag, "physical q coefficient");
    nonzero_q = nonzero_q || !lag.zero;
  }
  if (!nonzero_q)
    throw std::invalid_argument("physical cleared ODE has identically zero q");
  for (const auto& lag : ode.c_lags) {
    std::optional<std::pair<std::uint32_t, std::uint32_t>> previous;
    for (const auto& entry : lag) {
      if (entry.row >= ode.dimension || entry.column >= ode.dimension ||
          entry.value.zero)
        throw std::invalid_argument(
            "physical C lag has an invalid or structurally zero sparse entry");
      const auto position = std::pair{entry.row, entry.column};
      if (previous.has_value() && *previous >= position)
        throw std::invalid_argument(
            "physical C lag entries are not in canonical row-major order");
      previous = position;
      validate_rational(entry.value, "physical C coefficient");
    }
  }
}

inline std::int32_t shifted_power(std::int32_t power,
                                  std::int32_t shift,
                                  const char* label) {
  return local_detail::checked_i32(
      static_cast<std::int64_t>(power) + shift, label);
}

template <typename Scalar>
std::vector<ComplexBall> causal_product_component(
    const ExactEpsilonRational<Scalar>& value,
    const EpsilonVector& input, std::uint32_t component) {
  validate_rational(value, "physical epsilon-rational multiplier");
  if (value.zero || component >= input.dimension)
    throw std::invalid_argument(
        "causal physical multiplier received a zero entry or invalid component");
  const auto width = input.epsilon.width();
  std::vector<ComplexBall> output(width, ComplexBall(0));
  const auto q0 = local_detail::to_ball(value.denominator.front());
  if (q0.contains_zero())
    throw std::domain_error(
        "physical epsilon-rational denominator is not provably causal");
  for (std::size_t offset = 0; offset < width; ++offset) {
    const auto power = local_detail::checked_i32(
        static_cast<std::int64_t>(input.epsilon.min_power) + offset,
        "causal epsilon power");
    ComplexBall rhs(0);
    for (std::size_t degree = 0; degree < value.numerator.size(); ++degree) {
      const auto source64 = static_cast<std::int64_t>(power) -
                            static_cast<std::int64_t>(degree);
      if (source64 < input.epsilon.min_power ||
          source64 > input.epsilon.complete_max)
        continue;
      rhs += local_detail::to_ball(value.numerator[degree]) *
             input.at(static_cast<std::int32_t>(source64), component);
    }
    for (std::size_t degree = 1;
         degree < value.denominator.size() && degree <= offset; ++degree)
      rhs -= local_detail::to_ball(value.denominator[degree]) *
             output[offset - degree];
    output[offset] = rhs / q0;
  }
  return output;
}

template <typename Scalar>
struct WeightedTerm {
  const ExactEpsilonRational<Scalar>* value = nullptr;
  std::uint32_t row = 0;
  std::uint32_t column = 0;
  ComplexBall weight = ComplexBall(0);
};

template <typename Scalar>
EpsilonVector apply_terms(const std::vector<WeightedTerm<Scalar>>& raw_terms,
                          const EpsilonVector& input,
                          std::uint32_t output_dimension) {
  std::vector<const WeightedTerm<Scalar>*> terms;
  terms.reserve(raw_terms.size());
  for (const auto& term : raw_terms) {
    if (term.value == nullptr || term.value->zero || term.weight.is_zero())
      continue;
    if (term.row >= output_dimension || term.column >= input.dimension)
      throw std::invalid_argument("physical ODE term dimension is invalid");
    terms.push_back(&term);
  }
  if (terms.empty())
    throw std::domain_error(
        "physical ODE multiplier is structurally zero at the evaluation point");

  auto min_power = std::numeric_limits<std::int32_t>::max();
  auto complete_max = std::numeric_limits<std::int32_t>::max();
  for (const auto* term : terms) {
    min_power = std::min(min_power, shifted_power(
        input.epsilon.min_power, term->value->valuation,
        "physical product minimum"));
    complete_max = std::min(complete_max, shifted_power(
        input.epsilon.complete_max, term->value->valuation,
        "physical product complete maximum"));
  }
  if (complete_max < min_power)
    throw std::domain_error(
        "physical ODE product has no common complete epsilon window");

  EpsilonVector output;
  output.epsilon = {min_power, complete_max};
  output.dimension = output_dimension;
  output.coefficients.assign(output.epsilon.width() * output.dimension,
                             ComplexBall(0));
  for (const auto* term : terms) {
    auto causal = causal_product_component(
        *term->value, input, term->column);
    for (std::int32_t power = output.epsilon.min_power;
         power <= output.epsilon.complete_max; ++power) {
      const auto unshifted64 = static_cast<std::int64_t>(power) -
                               term->value->valuation;
      if (unshifted64 < input.epsilon.min_power ||
          unshifted64 > input.epsilon.complete_max)
        continue;
      const auto offset = static_cast<std::size_t>(
          unshifted64 - input.epsilon.min_power);
      output.at(power, term->row) += term->weight * causal[offset];
    }
  }
  return output;
}

template <typename Coefficient>
std::vector<Coefficient> convolve(
    const std::vector<Coefficient>& left,
    const std::vector<Coefficient>& right) {
  if (left.empty() || right.empty()) return {};
  std::vector<Coefficient> output(
      left.size() + right.size() - 1, ScalarTraits<Coefficient>::zero());
  for (std::size_t i = 0; i < left.size(); ++i)
    for (std::size_t j = 0; j < right.size(); ++j)
      output[i + j] += left[i] * right[j];
  return output;
}

// Prove that the clearing multiplier is not the zero rational function at
// this point.  A small q must not turn an arbitrary local into a vacuous pass.
// We form the common numerator exactly for Rational charts and with enclosure
// arithmetic for Acb charts; at least one Acb coefficient must exclude zero.
template <typename Scalar>
void require_nonvacuous_q(
    const PreparedPhysicalClearedODE<Scalar>& ode,
    const RealEvaluationPoint& point,
    const std::vector<ComplexBall>& t_powers) {
  if constexpr (std::is_same_v<Scalar, Rational>) {
    struct ExactActive {
      const ExactEpsilonRational<Rational>* value;
      Rational weight;
    };
    const Rational t(point.exact_coordinate);
    Rational t_power(1);
    std::vector<ExactActive> active;
    auto minimum_valuation = std::numeric_limits<std::int32_t>::max();
    for (const auto& value : ode.q_lags) {
      if (!value.zero && !t_power.is_zero()) {
        active.push_back({&value, t_power});
        minimum_valuation = std::min(minimum_valuation, value.valuation);
      }
      t_power *= t;
    }
    if (active.empty())
      throw std::domain_error(
          "physical clearing multiplier q vanishes at the evaluation point");
    std::map<std::int64_t, Rational> common_numerator;
    for (std::size_t index = 0; index < active.size(); ++index) {
      const auto& term = active[index];
      const auto shift =
          static_cast<std::int64_t>(term.value->valuation) -
          minimum_valuation;
      auto polynomial = term.value->numerator;
      for (std::size_t other = 0; other < active.size(); ++other)
        if (other != index)
          polynomial = convolve(
              polynomial, active[other].value->denominator);
      for (std::size_t degree = 0; degree < polynomial.size(); ++degree)
        common_numerator[shift + static_cast<std::int64_t>(degree)] +=
            term.weight * polynomial[degree];
    }
    if (std::any_of(common_numerator.begin(), common_numerator.end(),
                    [](const auto& item) {
                      return !item.second.is_zero();
                    }))
      return;
    throw std::domain_error(
        "physical clearing multiplier q is exactly zero at the evaluation point");
  }

  struct Active {
    const ExactEpsilonRational<Scalar>* value;
    ComplexBall weight;
  };
  std::vector<Active> active;
  std::int32_t minimum_valuation =
      std::numeric_limits<std::int32_t>::max();
  for (std::size_t lag = 0; lag < ode.q_lags.size(); ++lag) {
    const auto& value = ode.q_lags[lag];
    if (value.zero || t_powers.at(lag).is_zero()) continue;
    active.push_back({&value, t_powers[lag]});
    minimum_valuation = std::min(minimum_valuation, value.valuation);
  }
  if (active.empty())
    throw std::domain_error(
        "physical clearing multiplier q vanishes at the evaluation point");

  std::map<std::int64_t, ComplexBall> common_numerator;
  for (std::size_t index = 0; index < active.size(); ++index) {
    const auto& term = active[index];
    const auto shift =
        static_cast<std::int64_t>(term.value->valuation) -
        minimum_valuation;
    std::vector<ComplexBall> polynomial;
    for (const auto& coefficient : term.value->numerator)
      polynomial.push_back(local_detail::to_ball(coefficient));
    for (std::size_t other = 0; other < active.size(); ++other) {
      if (other == index) continue;
      std::vector<ComplexBall> denominator;
      denominator.reserve(active[other].value->denominator.size());
      for (const auto& coefficient : active[other].value->denominator)
        denominator.push_back(local_detail::to_ball(coefficient));
      polynomial = convolve(polynomial, denominator);
    }
    for (std::size_t degree = 0; degree < polynomial.size(); ++degree)
      common_numerator[shift + static_cast<std::int64_t>(degree)] +=
          term.weight * polynomial[degree];
  }
  if (std::any_of(common_numerator.begin(), common_numerator.end(),
                  [](const auto& item) {
                    return !Magnitude::lower_abs(item.second).is_zero();
                  }))
    return;
  if (std::all_of(common_numerator.begin(), common_numerator.end(),
                  [](const auto& item) { return item.second.is_zero(); }))
    throw std::domain_error(
        "physical clearing multiplier q is exactly zero at the evaluation point");
  throw std::domain_error(
      "physical clearing multiplier q is not provably nonzero at the evaluation point");
}

struct PreparedCausalEpsilonMultiplier {
  std::int32_t valuation = 0;
  std::vector<ComplexBall> numerator;
  std::vector<ComplexBall> denominator;
};

template <typename Scalar>
PreparedCausalEpsilonMultiplier prepare_causal_multiplier(
    const ExactEpsilonRational<Scalar>& value) {
  if (value.zero || value.valuation < 0)
    throw std::invalid_argument(
        "causal epsilon multiplier requires a nonzero nonnegative valuation");
  PreparedCausalEpsilonMultiplier prepared;
  prepared.valuation = value.valuation;
  prepared.numerator.reserve(value.numerator.size());
  prepared.denominator.reserve(value.denominator.size());
  for (const auto& coefficient : value.numerator)
    prepared.numerator.push_back(local_detail::to_ball(coefficient));
  for (const auto& coefficient : value.denominator)
    prepared.denominator.push_back(local_detail::to_ball(coefficient));
  if (prepared.numerator.empty() || prepared.denominator.empty() ||
      prepared.numerator.front().contains_zero() ||
      prepared.denominator.front().contains_zero())
    throw std::domain_error(
        "causal epsilon multiplier has no provably invertible leading data");
  return prepared;
}

inline EpsilonVector zero_epsilon_vector(EpsilonWindow window,
                                         std::uint32_t dimension) {
  EpsilonVector output;
  output.epsilon = window;
  output.dimension = dimension;
  const auto width = window.width();
  if (dimension != 0 &&
      width > std::numeric_limits<std::size_t>::max() / dimension)
    throw std::overflow_error("ordinary-center epsilon vector size overflow");
  output.coefficients.assign(width * dimension, ComplexBall(0));
  return output;
}

inline void accumulate_causal_multiplier(
    const PreparedCausalEpsilonMultiplier& multiplier,
    const EpsilonVector& source, std::uint32_t source_component,
    EpsilonVector& target, std::uint32_t target_component,
    const ComplexBall& weight) {
  if (source.epsilon.min_power != target.epsilon.min_power ||
      source.epsilon.complete_max != target.epsilon.complete_max ||
      source_component >= source.dimension ||
      target_component >= target.dimension || multiplier.valuation < 0 ||
      multiplier.numerator.empty() || multiplier.denominator.empty())
    throw std::invalid_argument(
        "ordinary-center causal product has inconsistent frames or dimensions");
  if (weight.is_zero()) return;
  const auto width = source.epsilon.width();
  std::vector<ComplexBall> quotient(width, ComplexBall(0));
  const auto denominator0 = multiplier.denominator.front();
  if (denominator0.contains_zero())
    throw std::domain_error(
        "ordinary-center causal denominator contains zero");
  for (std::size_t offset = 0; offset < width; ++offset) {
    ComplexBall rhs(0);
    for (std::size_t degree = 0;
         degree < multiplier.numerator.size() && degree <= offset;
         ++degree) {
      const auto source_power = local_detail::checked_i32(
          static_cast<std::int64_t>(source.epsilon.min_power) +
              static_cast<std::int64_t>(offset - degree),
          "ordinary-center causal source power");
      rhs += multiplier.numerator[degree] *
             source.at(source_power, source_component);
    }
    for (std::size_t degree = 1;
         degree < multiplier.denominator.size() && degree <= offset;
         ++degree)
      rhs -= multiplier.denominator[degree] * quotient[offset - degree];
    quotient[offset] = rhs / denominator0;

    const auto target_power64 =
        static_cast<std::int64_t>(source.epsilon.min_power) +
        static_cast<std::int64_t>(offset) + multiplier.valuation;
    if (target_power64 > target.epsilon.complete_max) continue;
    const auto target_power = local_detail::checked_i32(
        target_power64, "ordinary-center causal target power");
    target.at(target_power, target_component) += weight * quotient[offset];
  }
}

inline EpsilonVector solve_formal_unit_q0(
    const PreparedCausalEpsilonMultiplier& q0,
    const EpsilonVector& rhs, std::uint32_t taylor_index) {
  if (taylor_index == 0 || q0.valuation != 0 || q0.numerator.empty() ||
      q0.denominator.empty() || q0.numerator.front().contains_zero())
    throw std::invalid_argument(
        "ordinary-center q0 solve requires a positive Taylor index and a formal epsilon unit");
  auto output = zero_epsilon_vector(rhs.epsilon, rhs.dimension);
  const auto divisor = ComplexBall::from_strings(
      std::to_string(taylor_index));
  const auto width = rhs.epsilon.width();
  for (std::size_t offset = 0; offset < width; ++offset) {
    const auto power = local_detail::checked_i32(
        static_cast<std::int64_t>(rhs.epsilon.min_power) +
            static_cast<std::int64_t>(offset),
        "ordinary-center q0 solve power");
    for (std::uint32_t component = 0; component < rhs.dimension;
         ++component) {
      ComplexBall value(0);
      for (std::size_t degree = 0;
           degree < q0.denominator.size() && degree <= offset; ++degree) {
        const auto source_power = local_detail::checked_i32(
            static_cast<std::int64_t>(power) -
                static_cast<std::int64_t>(degree),
            "ordinary-center q0 numerator power");
        value += q0.denominator[degree] *
                 rhs.at(source_power, component);
      }
      value = value / divisor;
      for (std::size_t degree = 1;
           degree < q0.numerator.size() && degree <= offset; ++degree) {
        const auto previous_power = local_detail::checked_i32(
            static_cast<std::int64_t>(power) -
                static_cast<std::int64_t>(degree),
            "ordinary-center q0 previous power");
        value -= q0.numerator[degree] *
                 output.at(previous_power, component);
      }
      output.at(power, component) = value / q0.numerator.front();
    }
  }
  return output;
}

}  // namespace physical_ode_detail

struct OrdinaryCenterValueEvolution {
  bool eligible = false;
  std::string reason;
  EpsilonWindow epsilon;
  std::uint32_t dimension = 0;
  std::uint32_t taylor_complete_max = 0;
  // One full epsilon vector for each Taylor coefficient f_n, n=0..N.
  std::vector<EpsilonVector> taylor_coefficients;

  const EpsilonVector& at(std::uint32_t taylor_index) const {
    if (!eligible || taylor_index >= taylor_coefficients.size())
      throw std::out_of_range(
          "ordinary-center Taylor coefficient is unavailable");
    return taylor_coefficients[taylor_index];
  }
};

template <typename Scalar>
OrdinaryCenterValueEvolution evolve_ordinary_center_value(
    const PreparedPhysicalClearedODE<Scalar>& ode,
    const EpsilonVector& initial, std::uint32_t taylor_complete_max) {
  static_assert(std::is_same_v<Scalar, Rational> ||
                    std::is_same_v<Scalar, ComplexBall>,
                "ordinary-center value evolution supports Rational or Acb physical equations");
  physical_ode_detail::validate_ode(ode);
  const auto width = initial.epsilon.width();
  if (initial.dimension == 0 || initial.dimension != ode.dimension ||
      width > std::numeric_limits<std::size_t>::max() /
                  initial.dimension ||
      initial.coefficients.size() != width * initial.dimension ||
      !initial.error.empty())
    throw std::invalid_argument(
        "ordinary-center initial value has an invalid dimension, epsilon frame, or unabsorbed error envelope");

  OrdinaryCenterValueEvolution result;
  result.epsilon = initial.epsilon;
  result.dimension = initial.dimension;
  result.taylor_complete_max = taylor_complete_max;
  const auto ineligible = [&](std::string reason) {
    result.reason = std::move(reason);
    result.taylor_coefficients.clear();
    return result;
  };

  if (!ode.c_lags.front().empty())
    return ineligible(
        "physical equation is not ordinary at t=0 because C_0 is structurally nonzero");
  const auto& raw_q0 = ode.q_lags.front();
  if (raw_q0.zero || raw_q0.valuation != 0 || raw_q0.numerator.empty() ||
      raw_q0.denominator.empty() ||
      ScalarTraits<Scalar>::is_zero(raw_q0.numerator.front()) ||
      ScalarTraits<Scalar>::is_zero(raw_q0.denominator.front()))
    return ineligible(
        "physical equation q_0 is not a formal epsilon unit");
  for (std::size_t lag = 0; lag < ode.q_lags.size(); ++lag) {
    const auto& coefficient = ode.q_lags[lag];
    if (!coefficient.zero && coefficient.valuation < 0)
      return ineligible(
          "physical q lag " + std::to_string(lag) +
          " has a negative epsilon valuation and requires coefficients above the supplied window");
  }
  for (std::size_t lag = 0; lag < ode.c_lags.size(); ++lag)
    for (const auto& entry : ode.c_lags[lag])
      if (entry.value.valuation < 0)
        return ineligible(
            "physical C lag " + std::to_string(lag) + " entry (" +
            std::to_string(entry.row) + "," +
            std::to_string(entry.column) +
            ") has a negative epsilon valuation and requires coefficients above the supplied window");

  std::vector<std::optional<
      physical_ode_detail::PreparedCausalEpsilonMultiplier>> q_lags;
  q_lags.reserve(ode.q_lags.size());
  for (const auto& coefficient : ode.q_lags) {
    if (coefficient.zero)
      q_lags.push_back(std::nullopt);
    else
      q_lags.push_back(
          physical_ode_detail::prepare_causal_multiplier(coefficient));
  }
  struct PreparedCEntry {
    std::uint32_t row = 0;
    std::uint32_t column = 0;
    physical_ode_detail::PreparedCausalEpsilonMultiplier multiplier;
  };
  std::vector<std::vector<PreparedCEntry>> c_lags(ode.c_lags.size());
  for (std::size_t lag = 0; lag < ode.c_lags.size(); ++lag) {
    c_lags[lag].reserve(ode.c_lags[lag].size());
    for (const auto& entry : ode.c_lags[lag])
      c_lags[lag].push_back(PreparedCEntry{
          entry.row, entry.column,
          physical_ode_detail::prepare_causal_multiplier(entry.value)});
  }

  if (static_cast<std::uint64_t>(taylor_complete_max) + 1 >
      std::numeric_limits<std::size_t>::max())
    throw std::overflow_error(
        "ordinary-center Taylor coefficient count overflows size_t");
  result.taylor_coefficients.reserve(
      static_cast<std::size_t>(taylor_complete_max) + 1);
  result.taylor_coefficients.push_back(initial);
  const ComplexBall one(1);
  for (std::uint64_t raw_n = 1; raw_n <= taylor_complete_max; ++raw_n) {
    const auto n = static_cast<std::uint32_t>(raw_n);
    auto rhs = physical_ode_detail::zero_epsilon_vector(
        initial.epsilon, initial.dimension);
    const auto maximum_lag = std::min<std::size_t>(
        n, std::max(q_lags.size(), c_lags.size()) - 1);
    for (std::size_t lag = 1; lag <= maximum_lag; ++lag) {
      const auto& source = result.taylor_coefficients[n - lag];
      if (lag < c_lags.size())
        for (const auto& entry : c_lags[lag])
          physical_ode_detail::accumulate_causal_multiplier(
              entry.multiplier, source, entry.column, rhs, entry.row, one);
      if (lag < q_lags.size() && q_lags[lag].has_value() && n > lag) {
        const auto derivative_weight = -ComplexBall::from_strings(
            std::to_string(n - lag));
        for (std::uint32_t component = 0; component < initial.dimension;
             ++component)
          physical_ode_detail::accumulate_causal_multiplier(
              *q_lags[lag], source, component, rhs, component,
              derivative_weight);
      }
    }
    result.taylor_coefficients.push_back(
        physical_ode_detail::solve_formal_unit_q0(
            *q_lags.front(), rhs, n));
  }
  result.eligible = true;
  result.reason = "epsilon-causal ordinary-center physical evolution";
  return result;
}

inline ResidualCertificate certify_cleared_vector_residual(
    const EpsilonVector& lhs, const EpsilonVector& rhs,
    const Magnitude& relative_tolerance,
    ResidualScope scope = ResidualScope::StoredTruncation) {
  if (lhs.dimension == 0 || lhs.dimension != rhs.dimension)
    throw std::invalid_argument(
        "physical cleared residual vector dimensions disagree");
  const auto min_power = std::min(lhs.epsilon.min_power,
                                  rhs.epsilon.min_power);
  const auto complete_max = std::min(lhs.epsilon.complete_max,
                                     rhs.epsilon.complete_max);
  if (complete_max < min_power)
    throw std::domain_error(
        "physical cleared residual has no common complete epsilon window");

  ResidualCertificate result;
  result.scope = scope;
  result.residual.epsilon = {min_power, complete_max};
  result.residual.dimension = lhs.dimension;
  result.residual.coefficients.assign(
      result.residual.epsilon.width() * result.residual.dimension,
      ComplexBall(0));
  result.residual_upper.assign(result.residual.epsilon.width(),
                               Magnitude::zero());
  result.scale_lower.assign(result.residual.epsilon.width(),
                            Magnitude::one());
  std::vector<Magnitude> scale_upper(result.residual.epsilon.width(),
                                     Magnitude::one());
  std::vector<Magnitude> residual_lower(result.residual.epsilon.width(),
                                        Magnitude::zero());
  result.relative_upper.assign(result.residual.epsilon.width(),
                               Magnitude::zero());

  bool all_pass = true;
  bool any_fail = false;
  for (std::int32_t power = min_power; power <= complete_max; ++power) {
    const auto row = static_cast<std::size_t>(power - min_power);
    for (std::uint32_t component = 0; component < lhs.dimension;
         ++component) {
      const auto left = local_detail::coefficient_or_zero(
          lhs, power, component);
      const auto right = local_detail::coefficient_or_zero(
          rhs, power, component);
      const auto residual = left - right;
      result.residual.at(power, component) = residual;
      result.residual_upper[row] = Magnitude::maximum(
          result.residual_upper[row], Magnitude::upper_abs(residual));
      residual_lower[row] = Magnitude::maximum(
          residual_lower[row], Magnitude::lower_abs(residual));
      for (const auto& term : {left, right}) {
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
    result.verdict = ResidualVerdict::Inconclusive;
    result.detail =
        "Acb encloses the stored-truncation residual, but full local-solution "
        "certification requires certified value/theta tail bounds";
  } else if (all_pass) {
    result.verdict = ResidualVerdict::Pass;
    result.detail =
        "Acb enclosure proves the stored-truncation residual bound";
  } else if (any_fail) {
    result.verdict = ResidualVerdict::Fail;
    result.detail =
        "Acb lower bound disproves the stored-truncation tolerance";
  } else {
    result.verdict = ResidualVerdict::Inconclusive;
    result.detail =
        "residual enclosure overlaps the requested tolerance";
  }
  return result;
}

template <typename Scalar>
ResidualCertificate certify_physical_cleared_ode_residual(
    const PreparedPhysicalClearedODE<Scalar>& ode,
    const LocalEvaluation& evaluation,
    const RealEvaluationPoint& point,
    const Magnitude& relative_tolerance,
    ResidualScope scope = ResidualScope::StoredTruncation) {
  static_assert(std::is_same_v<Scalar, Rational> ||
                    std::is_same_v<Scalar, ComplexBall>,
                "physical residuals support Rational or Acb locals only");
  physical_ode_detail::validate_ode(ode);
  if (evaluation.value.dimension != ode.dimension ||
      evaluation.theta_value.dimension != ode.dimension)
    throw std::invalid_argument(
        "physical residual evaluation dimension differs from its owner");

  ComplexBall signed_t = point.modulus;
  if (point.sign < 0) signed_t = -signed_t;
  const auto t_width = std::max(ode.q_lags.size(), ode.c_lags.size());
  std::vector<ComplexBall> t_powers(t_width, ComplexBall(1));
  for (std::size_t lag = 1; lag < t_width; ++lag)
    t_powers[lag] = t_powers[lag - 1] * signed_t;
  physical_ode_detail::require_nonvacuous_q(ode, point, t_powers);

  std::vector<physical_ode_detail::WeightedTerm<Scalar>> q_terms;
  q_terms.reserve(ode.q_lags.size() * ode.dimension);
  for (std::size_t lag = 0; lag < ode.q_lags.size(); ++lag)
    for (std::uint32_t component = 0; component < ode.dimension; ++component)
      q_terms.push_back(
          {&ode.q_lags[lag], component, component, t_powers[lag]});
  auto lhs = physical_ode_detail::apply_terms(
      q_terms, evaluation.theta_value, ode.dimension);

  std::vector<physical_ode_detail::WeightedTerm<Scalar>> c_terms;
  for (std::size_t lag = 0; lag < ode.c_lags.size(); ++lag)
    for (const auto& entry : ode.c_lags[lag])
      c_terms.push_back(
          {&entry.value, entry.row, entry.column, t_powers[lag]});
  EpsilonVector rhs;
  if (c_terms.empty()) {
    rhs.epsilon = lhs.epsilon;
    rhs.dimension = ode.dimension;
    rhs.coefficients.assign(rhs.epsilon.width() * rhs.dimension,
                            ComplexBall(0));
  } else {
    rhs = physical_ode_detail::apply_terms(
        c_terms, evaluation.value, ode.dimension);
  }

  auto certificate = certify_cleared_vector_residual(
      lhs, rhs, relative_tolerance, scope);
  if (scope == ResidualScope::StoredTruncation)
    certificate.detail +=
        "; equation derived from retained physical owner payload "
        "q(t,eps) theta(f)=C(t,eps) f (homogeneous zero source)";
  else
    certificate.detail +=
        "; physical q/C ownership is proven, but the requested full-local "
        "residual still requires certified value/theta tail propagation";
  return certificate;
}

}  // namespace diffexp2
