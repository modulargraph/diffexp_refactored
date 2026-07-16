#include "diffexp2/physical_ode.hpp"

#include <iostream>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using diffexp2::ComplexBall;
using diffexp2::EpsilonVector;
using diffexp2::EpsilonWindow;
using diffexp2::ExactEpsilonRational;
using diffexp2::PhysicalODEMatrixEntry;
using diffexp2::PreparedPhysicalClearedODE;
using diffexp2::Rational;

ExactEpsilonRational<Rational> rational(
    std::int32_t valuation, std::initializer_list<const char*> numerator,
    std::initializer_list<const char*> denominator = {"1"}) {
  ExactEpsilonRational<Rational> value;
  value.zero = false;
  value.valuation = valuation;
  for (const auto* coefficient : numerator)
    value.numerator.emplace_back(std::string(coefficient));
  for (const auto* coefficient : denominator)
    value.denominator.emplace_back(std::string(coefficient));
  return value;
}

template <typename Scalar>
PreparedPhysicalClearedODE<Scalar> equation_shell(std::uint32_t dimension) {
  PreparedPhysicalClearedODE<Scalar> equation;
  equation.dimension = dimension;
  equation.owner_signature_identity = "de2-equation-evolution-fixture";
  equation.payload_identity = "de2-physical-ode-evolution-fixture";
  equation.exact_payload_record = "physical-ode-evolution-fixture";
  return equation;
}

EpsilonVector vector(EpsilonWindow window,
                     std::initializer_list<const char*> coefficients) {
  EpsilonVector value;
  value.epsilon = window;
  value.dimension = 1;
  for (const auto* coefficient : coefficients)
    value.coefficients.push_back(ComplexBall::from_strings(coefficient));
  if (value.coefficients.size() != window.width())
    throw std::invalid_argument("test epsilon vector has the wrong width");
  return value;
}

void require_exact(const ComplexBall& actual, const char* expected,
                   const char* label) {
  if (!(actual - ComplexBall::from_strings(expected)).is_zero())
    throw std::runtime_error(std::string(label) + " differs from " +
                             expected);
}

PreparedPhysicalClearedODE<Rational> unit_division_equation() {
  auto equation = equation_shell<Rational>(1);
  equation.q_lags = {rational(0, {"1", "1"})};
  equation.c_lags.resize(2);
  equation.c_lags[1].push_back(
      PhysicalODEMatrixEntry<Rational>{0, 0, rational(0, {"1"})});
  return equation;
}

}  // namespace

int main() {
  try {
    ComplexBall::set_precision(256);

    // q0=1+eps, C=t gives
    //   f1=(1+eps)^-1 f0,
    //   f2=1/2 (1+eps)^-2 f0.
    // The finite lower edge is structural zero, so both solves are causal in
    // ascending epsilon and retain the complete supplied maximum.
    const auto divided = diffexp2::evolve_ordinary_center_value(
        unit_division_equation(), vector({-1, 3}, {"1", "0", "0", "0", "0"}),
        2);
    if (!divided.eligible || divided.dimension != 1 ||
        divided.epsilon.min_power != -1 ||
        divided.epsilon.complete_max != 3 ||
        divided.taylor_complete_max != 2 ||
        divided.taylor_coefficients.size() != 3)
      throw std::runtime_error(
          "ordinary-center unit division lost its output contract");
    const std::vector<const char*> expected_f1{"1", "-1", "1", "-1",
                                                "1"};
    const std::vector<const char*> expected_f2{"1/2", "-1", "3/2", "-2",
                                                "5/2"};
    for (std::size_t offset = 0; offset < expected_f1.size(); ++offset) {
      const auto power = static_cast<std::int32_t>(offset) - 1;
      require_exact(divided.at(1).at(power, 0), expected_f1[offset],
                    "unit-division f1");
      require_exact(divided.at(2).at(power, 0), expected_f2[offset],
                    "unit-division f2");
    }

    // A positive valuation is causal.  eps*f0 must leave the lower edge
    // exactly zero while preserving every coefficient through the supplied
    // complete maximum.
    auto positive = equation_shell<Rational>(1);
    positive.q_lags = {rational(0, {"1"})};
    positive.c_lags.resize(2);
    positive.c_lags[1].push_back(
        PhysicalODEMatrixEntry<Rational>{0, 0, rational(1, {"1"})});
    const auto shifted = diffexp2::evolve_ordinary_center_value(
        positive, vector({-1, 2}, {"2", "3", "5", "7"}), 1);
    if (!shifted.eligible || shifted.at(1).epsilon.min_power != -1 ||
        shifted.at(1).epsilon.complete_max != 2)
      throw std::runtime_error(
          "positive-valuation evolution changed its epsilon window");
    const std::vector<const char*> expected_shift{"0", "2", "3", "5"};
    for (std::size_t offset = 0; offset < expected_shift.size(); ++offset)
      require_exact(shifted.at(1).at(
                        static_cast<std::int32_t>(offset) - 1, 0),
                    expected_shift[offset], "positive-valuation lower edge");

    auto noncausal = equation_shell<Rational>(1);
    noncausal.q_lags = {rational(0, {"1"})};
    noncausal.c_lags.resize(2);
    noncausal.c_lags[1].push_back(
        PhysicalODEMatrixEntry<Rational>{0, 0, rational(-1, {"1"})});
    const auto rejected = diffexp2::evolve_ordinary_center_value(
        noncausal, vector({-1, 2}, {"2", "3", "5", "7"}), 1);
    if (rejected.eligible || !rejected.taylor_coefficients.empty() ||
        rejected.reason.find("negative epsilon valuation") ==
            std::string::npos)
      throw std::runtime_error(
          "negative-valuation physical evolution was not ineligible");

    auto nonordinary = unit_division_equation();
    nonordinary.c_lags[0].push_back(
        PhysicalODEMatrixEntry<Rational>{0, 0, rational(0, {"1"})});
    const auto nonordinary_rejected =
        diffexp2::evolve_ordinary_center_value(
            nonordinary, vector({0, 1}, {"1", "0"}), 1);
    if (nonordinary_rejected.eligible ||
        nonordinary_rejected.reason.find("not ordinary") ==
            std::string::npos)
      throw std::runtime_error(
          "center-nonvanishing C0 acquired an ordinary evolution");

    auto ambiguous = equation_shell<ComplexBall>(1);
    ExactEpsilonRational<ComplexBall> ambiguous_q0;
    ambiguous_q0.zero = false;
    ambiguous_q0.valuation = 0;
    ambiguous_q0.numerator = {
        ComplexBall::from_strings("[0 +/- 1]")};
    ambiguous_q0.denominator = {ComplexBall(1)};
    ambiguous.q_lags = {std::move(ambiguous_q0)};
    ambiguous.c_lags.resize(1);
    bool ambiguous_rejected = false;
    try {
      (void)diffexp2::evolve_ordinary_center_value(
          ambiguous, vector({0, 1}, {"1", "0"}), 1);
    } catch (const std::invalid_argument& error) {
      ambiguous_rejected =
          std::string(error.what()).find("leading numerator") !=
          std::string::npos;
    }
    if (!ambiguous_rejected)
      throw std::runtime_error(
          "Acb q0 constant containing zero passed physical validation");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
