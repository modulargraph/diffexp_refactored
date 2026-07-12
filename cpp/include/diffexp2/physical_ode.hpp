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

}  // namespace physical_ode_detail

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
