#include "diffexp/kernel/singular_tail_majorant.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using diffexp::kernel::ComplexBall;
using diffexp::kernel::EpsilonWindow;
using diffexp::kernel::ExactEpsilonRational;
using diffexp::kernel::LocalSector;
using diffexp::kernel::LocalSolution;
using diffexp::kernel::Magnitude;
using diffexp::kernel::PhysicalODEMatrixEntry;
using diffexp::kernel::PreparedPhysicalClearedODE;
using diffexp::kernel::Rational;
using diffexp::kernel::TailMajorantStatus;

ExactEpsilonRational<Rational>
rational(std::int32_t valuation,
         std::initializer_list<const char *> numerator) {
  ExactEpsilonRational<Rational> value;
  value.zero = false;
  value.valuation = valuation;
  for (const auto *coefficient : numerator)
    value.numerator.emplace_back(std::string(coefficient));
  value.denominator.emplace_back("1");
  return value;
}

PreparedPhysicalClearedODE<Rational> exponential_equation(const char *a) {
  PreparedPhysicalClearedODE<Rational> equation;
  equation.dimension = 1;
  equation.q_lags = {rational(0, {"1"})};
  equation.c_lags.resize(2);
  equation.c_lags[0].push_back(
      PhysicalODEMatrixEntry<Rational>{0, 0, rational(0, {a})});
  equation.c_lags[1].push_back(
      PhysicalODEMatrixEntry<Rational>{0, 0, rational(0, {"1"})});
  equation.owner_signature_identity = "singular-tail-exponential-owner";
  equation.payload_identity = "singular-tail-exponential-payload";
  equation.exact_payload_record = "theta-f=(a+t)-f";
  return equation;
}

LocalSolution<Rational> exponential_local(const char *a, std::uint32_t order) {
  LocalSolution<Rational> solution;
  solution.chart.center_exact = "0";
  solution.chart.scale_exact = "1";
  solution.chart.radius_exact = "1";
  solution.chart.radius = ComplexBall(1);
  solution.epsilon = {0, 0};
  solution.taylor_complete_max = order;
  solution.dimension = 1;
  LocalSector<Rational> sector;
  sector.a = diffexp::kernel::ExactScalarDescriptor::rational(a);
  sector.b = diffexp::kernel::ExactScalarDescriptor::rational("0");
  sector.log_power = 0;
  Rational factorial(1);
  for (std::uint32_t n = 0; n <= order; ++n) {
    if (n > 0)
      factorial *= Rational(std::to_string(n));
    sector.coefficients.push_back(Rational(1) / factorial);
  }
  solution.sectors.push_back(std::move(sector));
  solution.checkpoint_identity = "singular-tail-exponential-local";
  return solution;
}

PreparedPhysicalClearedODE<Rational>
constant_equation(std::uint32_t dimension) {
  PreparedPhysicalClearedODE<Rational> equation;
  equation.dimension = dimension;
  equation.q_lags = {rational(0, {"1"})};
  equation.c_lags.resize(1);
  equation.owner_signature_identity = "singular-tail-constant-owner";
  equation.payload_identity = "singular-tail-constant-payload";
  equation.exact_payload_record = "constant-singular-tail-fixture";
  return equation;
}

LocalSolution<Rational> local_shell(std::uint32_t dimension,
                                    std::uint32_t order, std::string identity) {
  LocalSolution<Rational> solution;
  solution.chart.center_exact = "0";
  solution.chart.scale_exact = "1";
  solution.chart.radius_exact = "1";
  solution.chart.radius = ComplexBall(1);
  solution.epsilon = {0, 0};
  solution.taylor_complete_max = order;
  solution.dimension = dimension;
  solution.checkpoint_identity = std::move(identity);
  return solution;
}

LocalSector<Rational> zero_sector(const char *a, std::uint32_t log_power,
                                  std::uint32_t dimension,
                                  std::uint32_t order) {
  LocalSector<Rational> sector;
  sector.a = diffexp::kernel::ExactScalarDescriptor::rational(a);
  sector.b = diffexp::kernel::ExactScalarDescriptor::rational("0");
  sector.log_power = log_power;
  sector.coefficients.assign((static_cast<std::size_t>(order) + 1) * dimension,
                             Rational(0));
  return sector;
}

void require_contains(const ComplexBall &enclosure, const ComplexBall &exact,
                      const char *label) {
  if (!acb_contains(enclosure.raw(), exact.raw()))
    throw std::runtime_error(std::string(label) +
                             " does not contain the exact value");
}

} // namespace

int main() {
  try {
    ComplexBall::set_precision(256);

    const auto equation = exponential_equation("1");
    const auto local = exponential_local("1", 4);
    const auto model_result =
        diffexp::kernel::prepare_singular_rational_shadow_tail_model(
            equation, local, EpsilonWindow{0, 0}, 4, "1/2");
    if (model_result.status != TailMajorantStatus::Certified ||
        !model_result.model.has_value() ||
        model_result.model->towers.size() != 1 ||
        !(model_result.model->towers.front()
              .contraction_upper.front()
              .approximate_upper() < 1.0))
      throw std::runtime_error("scalar Frobenius recurrence did not acquire an "
                               "all-future certificate: " +
                               model_result.detail);

    const auto point = diffexp::kernel::RealEvaluationPoint::rational("1/4");
    const auto certified =
        diffexp::kernel::evaluate_singular_rational_shadow_with_certified_tail(
            *model_result.model, local, point);
    if (certified.tail.status != TailMajorantStatus::Certified ||
        certified.tail.value.absolute.size() != 1 ||
        certified.tail.value.absolute.front().is_zero())
      throw std::runtime_error("scalar Frobenius point tail was not certified");
    const auto inflated =
        diffexp::kernel::inflate_certified_singular_seed_value(certified);
    const auto quarter = ComplexBall::from_strings("1/4");
    const auto exact = quarter * diffexp::kernel::local_detail::cb_exp(quarter);
    require_contains(inflated.at(0, 0), exact,
                     "inflated scalar Frobenius seed");

    const auto prefix_certificate =
        diffexp::kernel::certify_singular_rational_shadow_prefix(
            equation, local, EpsilonWindow{0, 0}, 4);
    const auto bridge =
        diffexp::kernel::certify_singular_rational_shadow_ordinary_bridge(
            equation, local, EpsilonWindow{0, 0}, 4,
            diffexp::kernel::RealEvaluationPoint::rational("1/6"), "1/2", point, 16,
            "1/8", {}, &prefix_certificate);
    if (bridge.status != TailMajorantStatus::Certified)
      throw std::runtime_error(
          "scalar singular-to-ordinary bridge was not certified: " +
          bridge.detail);
    if (bridge.detail.find("ordinary_tail_method=all-future-q/C-recurrence") ==
        std::string::npos)
      throw std::runtime_error(
          "scalar singular-to-ordinary bridge did not use the recurrence "
          "tail: " +
          bridge.detail);
    require_contains(bridge.evaluation.value.at(0, 0), exact,
                     "bridged scalar Frobenius value");
    const auto outward_bridge =
        diffexp::kernel::certify_singular_rational_shadow_ordinary_bridge(
            equation, local, EpsilonWindow{0, 0}, 4,
            diffexp::kernel::RealEvaluationPoint::rational("1/6"), "1/2", point, 16,
            "31/192", {}, &prefix_certificate);
    if (outward_bridge.status != TailMajorantStatus::Certified)
      throw std::runtime_error(
          "outward singular-to-ordinary witness was not certified: " +
          outward_bridge.detail);
    if (!(outward_bridge.ordinary_tail.value.absolute.front()
              .approximate_upper() <
          bridge.ordinary_tail.value.absolute.front().approximate_upper()))
      throw std::runtime_error(
          "outward recurrence witness did not tighten the certified tail: "
          "midpoint=" +
          bridge.ordinary_tail.value.absolute.front().dump_exact() +
          "; outward=" +
          outward_bridge.ordinary_tail.value.absolute.front().dump_exact());
    require_contains(outward_bridge.evaluation.value.at(0, 0), exact,
                     "outward bridged scalar Frobenius value");

    // The singular tail and ordinary continuation need not use the same
    // frame.  h=exp(t) is Fuchsian with theta h=t h, while f=t h obeys the
    // physically irregular cleared equation
    //
    //                  t theta f = (t + t^2) f.
    //
    // Certify the tail in h, map the enclosed seed by f=t h, and continue
    // only then in the original physical equation.
    auto tail_equation = constant_equation(1);
    tail_equation.c_lags.resize(2);
    tail_equation.c_lags[1].push_back(
        PhysicalODEMatrixEntry<Rational>{0, 0, rational(0, {"1"})});
    const auto tail_local = exponential_local("0", 4);
    auto irregular_physical = constant_equation(1);
    irregular_physical.q_lags = {ExactEpsilonRational<Rational>{},
                                 rational(0, {"1"})};
    irregular_physical.c_lags.resize(3);
    irregular_physical.c_lags[1].push_back(
        PhysicalODEMatrixEntry<Rational>{0, 0, rational(0, {"1"})});
    irregular_physical.c_lags[2].push_back(
        PhysicalODEMatrixEntry<Rational>{0, 0, rational(0, {"1"})});
    const auto tail_prefix = diffexp::kernel::certify_singular_rational_shadow_prefix(
        tail_equation, tail_local, EpsilonWindow{0, 0}, 4);
    const auto framed_bridge = diffexp::kernel::
        certify_singular_rational_shadow_ordinary_bridge_with_seed_map(
            tail_equation, tail_local, irregular_physical,
            [](const diffexp::kernel::EpsilonVector &seed,
               const diffexp::kernel::RealEvaluationPoint &seed_point)
                -> std::optional<diffexp::kernel::EpsilonVector> {
              auto physical = seed;
              const auto factor =
                  ComplexBall::from_strings(seed_point.exact_coordinate);
              for (auto &coefficient : physical.coefficients)
                coefficient *= factor;
              return physical;
            },
            EpsilonWindow{0, 0}, 4,
            diffexp::kernel::RealEvaluationPoint::rational("1/5"), "1/2", point, 16,
            "3/50", {}, &tail_prefix);
    if (framed_bridge.status != TailMajorantStatus::Certified)
      throw std::runtime_error(
          "frame-changing singular-to-ordinary bridge was not certified: " +
          framed_bridge.detail);
    require_contains(framed_bridge.evaluation.value.at(0, 0), exact,
                     "frame-changing bridged physical value");

    // A private reservoir may prepend a structurally zero epsilon row.  Local
    // evaluation legitimately trims that row; the full physical tail theorem
    // must be sliced to the public evaluated frame before the bridge inflates
    // the result.
    const auto leading_zero_frame_bridge = diffexp::kernel::
        certify_singular_rational_shadow_ordinary_bridge_with_seed_map(
            equation, local, equation,
            [](const diffexp::kernel::EpsilonVector &seed,
               const diffexp::kernel::RealEvaluationPoint &)
                -> std::optional<diffexp::kernel::EpsilonVector> {
              diffexp::kernel::EpsilonVector physical;
              physical.epsilon = {-1, 0};
              physical.dimension = seed.dimension;
              physical.coefficients.assign(
                  physical.epsilon.width() * physical.dimension,
                  ComplexBall(0));
              physical.at(0, 0) = seed.at(0, 0);
              return physical;
            },
            EpsilonWindow{0, 0}, 4,
            diffexp::kernel::RealEvaluationPoint::rational("1/6"), "1/2", point, 16,
            "1/8", {}, &prefix_certificate);
    if (leading_zero_frame_bridge.status != TailMajorantStatus::Certified ||
        leading_zero_frame_bridge.evaluation.value.epsilon.min_power != 0 ||
        leading_zero_frame_bridge.evaluation.value.epsilon.complete_max != 0)
      throw std::runtime_error(
          "leading-zero private epsilon frame blocked the singular bridge: " +
          leading_zero_frame_bridge.detail);
    require_contains(leading_zero_frame_bridge.evaluation.value.at(0, 0),
                     exact, "leading-zero-frame bridged value");

    auto common_monomial_equation = equation;
    for (auto &value : common_monomial_equation.q_lags)
      if (!value.zero)
        value.valuation -= 2;
    for (auto &lag : common_monomial_equation.c_lags)
      for (auto &entry : lag)
        entry.value.valuation -= 2;
    const auto common_monomial_bridge =
        diffexp::kernel::certify_singular_rational_shadow_ordinary_bridge(
            common_monomial_equation, local, EpsilonWindow{0, 0}, 4,
            diffexp::kernel::RealEvaluationPoint::rational("1/6"), "1/2", point, 16,
            "1/8");
    if (common_monomial_bridge.status != TailMajorantStatus::Certified)
      throw std::runtime_error(
          "removable common epsilon monomial blocked the singular bridge: " +
          common_monomial_bridge.detail);
    require_contains(common_monomial_bridge.evaluation.value.at(0, 0), exact,
                     "common-monomial bridged Frobenius value");
    const auto negative_bridge =
        diffexp::kernel::certify_singular_rational_shadow_ordinary_bridge(
            equation, local, EpsilonWindow{0, 0}, 4,
            diffexp::kernel::RealEvaluationPoint::rational("-1/6"), "1/2",
            diffexp::kernel::RealEvaluationPoint::rational("-1/4"), 16, "1/8");
    if (negative_bridge.status != TailMajorantStatus::Certified)
      throw std::runtime_error(
          "negative-ray singular-to-ordinary bridge was not certified: " +
          negative_bridge.detail);
    const auto negative_quarter = ComplexBall::from_strings("-1/4");
    const auto negative_exact =
        negative_quarter * diffexp::kernel::local_detail::cb_exp(negative_quarter);
    require_contains(negative_bridge.evaluation.value.at(0, 0), negative_exact,
                     "negative-ray bridged scalar Frobenius value");

    auto high_epsilon_coupling = constant_equation(1);
    high_epsilon_coupling.c_lags.resize(2);
    high_epsilon_coupling.c_lags[1].push_back(
        PhysicalODEMatrixEntry<Rational>{0, 0, rational(1, {"8192"})});
    const auto high_weight_contraction = diffexp::kernel::
        singular_tail_majorant_detail::certify_interval_transfer_contraction(
            high_epsilon_coupling, EpsilonWindow{0, 1}, Rational(0),
            Rational(0), 0, 0, 5, Rational(1));
    if (high_weight_contraction.status != TailMajorantStatus::Certified ||
        high_weight_contraction.epsilon_weight_base <= 1024)
      throw std::runtime_error(
          "ordinary recurrence contraction did not search beyond the old "
          "epsilon-weight ceiling: " +
          high_weight_contraction.detail);
    auto radius_sensitive = constant_equation(1);
    radius_sensitive.c_lags.resize(2);
    radius_sensitive.c_lags[1].push_back(
        PhysicalODEMatrixEntry<Rational>{0, 0, rational(0, {"100"})});
    const auto loose_radius_contraction = diffexp::kernel::
        singular_tail_majorant_detail::certify_interval_transfer_contraction(
            radius_sensitive, EpsilonWindow{0, 0}, Rational(0), Rational(0),
            0, 0, 17, Rational("1/5"));
    const auto tight_radius_contraction = diffexp::kernel::
        singular_tail_majorant_detail::certify_interval_transfer_contraction(
            radius_sensitive, EpsilonWindow{0, 0}, Rational(0), Rational(0),
            0, 0, 17, Rational("1/10"));
    if (loose_radius_contraction.status != TailMajorantStatus::Certified ||
        loose_radius_contraction.contraction_upper.back()
                .approximate_upper() <=
            1.0 ||
        tight_radius_contraction.status != TailMajorantStatus::Certified ||
        tight_radius_contraction.contraction_upper.back()
                .approximate_upper() >=
            1.0)
      throw std::runtime_error(
          "ordinary recurrence contraction did not distinguish loose and "
          "tight witness radii");

    const auto too_large =
        diffexp::kernel::prepare_singular_rational_shadow_tail_model(
            equation, local, EpsilonWindow{0, 0}, 4, "6");
    if (too_large.status != TailMajorantStatus::Inconclusive ||
        too_large.model.has_value() ||
        too_large.detail.find("rho") == std::string::npos)
      throw std::runtime_error(
          "noncontracting Frobenius witness radius did not fail closed: "
          "status=" +
          std::to_string(static_cast<int>(too_large.status)) +
          "; detail=" + too_large.detail +
          (too_large.model.has_value()
               ? "; rho=" + std::to_string(too_large.model->towers.front()
                                               .contraction_upper.front()
                                               .approximate_upper())
               : ""));

    auto negative_valuation = exponential_equation("1");
    negative_valuation.c_lags[1][0].value.valuation = -1;
    auto structurally_zero_local = local;
    std::fill(structurally_zero_local.sectors.front().coefficients.begin(),
              structurally_zero_local.sectors.front().coefficients.end(),
              Rational(0));
    const auto negative = diffexp::kernel::prepare_singular_rational_shadow_tail_model(
        negative_valuation, structurally_zero_local, EpsilonWindow{0, 0}, 4,
        "1/2");
    if (negative.status == TailMajorantStatus::Certified ||
        negative.model.has_value())
      throw std::runtime_error(
          "noncausal epsilon stack was not rejected after common cancellation: "
          "status=" +
          std::to_string(static_cast<int>(negative.status)) +
          "; detail=" + negative.detail);

    // Integer-shift sectors are representations of one analytic tower, not
    // independent tails.  Here t from the a=0 sector cancels -t from the
    // a=1 sector.  The exact zero solution of theta f=0 must therefore
    // acquire one zero-tail tower.  Perturbing the cancellation must be
    // rejected by the exact q/C residual before any bound is issued.
    const auto constant = constant_equation(1);
    auto shifted = local_shell(1, 1, "singular-tail-integer-shift-local");
    auto a0 = zero_sector("0", 0, 1, 1);
    a0.coefficients[1] = Rational(1);
    auto a1 = zero_sector("1", 0, 1, 1);
    a1.coefficients[0] = Rational(-1);
    shifted.sectors = {a0, a1};
    const auto shifted_model =
        diffexp::kernel::prepare_singular_rational_shadow_tail_model(
            constant, shifted, EpsilonWindow{0, 0}, 1, "1/2");
    if (shifted_model.status != TailMajorantStatus::Certified ||
        !shifted_model.model.has_value() ||
        shifted_model.model->towers.size() != 1 ||
        !shifted_model.model->towers.front()
             .induction_m_upper.front()
             .is_zero())
      throw std::runtime_error(
          "integer-shift Frobenius sectors were not canonicalized together");
    auto perturbed = shifted;
    perturbed.sectors[1].coefficients[0] = Rational(-2);
    const auto rejected_perturbation =
        diffexp::kernel::prepare_singular_rational_shadow_tail_model(
            constant, perturbed, EpsilonWindow{0, 0}, 1, "1/2");
    if (rejected_perturbation.status == TailMajorantStatus::Certified ||
        rejected_perturbation.model.has_value())
      throw std::runtime_error(
          "perturbed integer-shift tower retained a false exact certificate");

    // The stored log^p/p! convention is tested with the Jordan solution
    // f=(log(t),1) of theta f = [[0,1],[0,0]] f.  Its finite prefix is the
    // whole solution, so the certified omitted tail is exactly zero on the
    // negative rim as well.
    auto jordan_equation = constant_equation(2);
    jordan_equation.c_lags[0].push_back(
        PhysicalODEMatrixEntry<Rational>{0, 1, rational(0, {"1"})});
    auto jordan = local_shell(2, 2, "singular-tail-jordan-log-local");
    jordan.epsilon = {0, 1};
    auto jordan_log0 = zero_sector("0", 0, 2, 2);
    jordan_log0.coefficients.assign(12, Rational(0));
    jordan_log0.coefficients[7] = Rational(1);
    auto jordan_log1 = zero_sector("0", 1, 2, 2);
    jordan_log1.coefficients.assign(12, Rational(0));
    jordan_log1.coefficients[0] = Rational(1);
    jordan.sectors = {jordan_log0, jordan_log1};
    jordan.prescriptions.push_back(diffexp::kernel::Prescription{"t", 1, 1, 1});
    const auto jordan_model =
        diffexp::kernel::prepare_singular_rational_shadow_tail_model(
            jordan_equation, jordan, EpsilonWindow{0, 1}, 2, "1/2");
    if (jordan_model.status != TailMajorantStatus::Certified ||
        !jordan_model.model.has_value() ||
        jordan_model.model->towers.size() != 1 ||
        jordan_model.model->towers.front().log_complete_max != 1)
      throw std::runtime_error(
          "Jordan/log Frobenius tower was not certified: " +
          jordan_model.detail);
    diffexp::kernel::EvaluationOptions lower_rim;
    lower_rim.imaginary_sign = -1;
    const auto jordan_seed =
        diffexp::kernel::evaluate_singular_rational_shadow_with_certified_tail(
            *jordan_model.model, jordan,
            diffexp::kernel::RealEvaluationPoint::rational("-1/4"), lower_rim);
    if (jordan_seed.tail.status != TailMajorantStatus::Certified ||
        !jordan_seed.tail.value.absolute.front().is_zero())
      throw std::runtime_error(
          "finite Jordan/log solution acquired a nonzero omitted tail");
    const auto jordan_value =
        diffexp::kernel::inflate_certified_singular_seed_value(jordan_seed);
    auto expected_log =
        diffexp::kernel::local_detail::cb_log(ComplexBall::from_strings("1/4"));
    expected_log += diffexp::kernel::local_detail::imaginary_pi(-1);
    require_contains(jordan_value.at(1, 0), expected_log,
                     "negative-rim Jordan logarithm");
    require_contains(jordan_value.at(1, 1), ComplexBall(1),
                     "negative-rim Jordan constant");

    // A large nilpotent coupling can make the old whole-matrix Neumann
    // condition ||G||/n0<1 fail even though every future block is
    // invertible.  The structured proof checks only the finite dangerous
    // base resolvents and handles the strict triangular coupling by forward
    // substitution.
    auto nonnormal = constant_equation(2);
    nonnormal.c_lags[0].push_back(
        PhysicalODEMatrixEntry<Rational>{1, 0, rational(0, {"100"})});
    const auto nonnormal_future = diffexp::kernel::singular_tail_majorant_detail::
        certify_structured_future_inverse(nonnormal, EpsilonWindow{0, 0},
                                          Rational(0), Rational(0), 0, 0, 5);
    if (nonnormal_future.status != TailMajorantStatus::Certified ||
        nonnormal_future.exact_resolvent_checks == 0 ||
        nonnormal_future.inverse_n_times_prefix_upper.size() != 1)
      throw std::runtime_error(
          "large strict-triangular future recurrence was not certified: " +
          nonnormal_future.detail);

    // Conversely, an actual future indicial resonance must be found by an
    // exact finite resolvent check, not hidden by the asymptotic argument.
    auto resonant = constant_equation(1);
    resonant.c_lags[0].push_back(
        PhysicalODEMatrixEntry<Rational>{0, 0, rational(0, {"10"})});
    const auto resonant_future = diffexp::kernel::singular_tail_majorant_detail::
        certify_structured_future_inverse(resonant, EpsilonWindow{0, 0},
                                          Rational(0), Rational(0), 0, 0, 5);
    if (resonant_future.status == TailMajorantStatus::Certified ||
        resonant_future.detail.find("resonance at Taylor index 10") ==
            std::string::npos)
      throw std::runtime_error(
          "true future Frobenius resonance was not rejected exactly");

    // b*epsilon and the logarithmic derivative ladder are nilpotent on a
    // finite causal prefix.  Arbitrarily large values can increase the
    // inverse bound, but they cannot create a false spectral obstruction.
    const auto nilpotent_future = diffexp::kernel::singular_tail_majorant_detail::
        certify_structured_future_inverse(constant, EpsilonWindow{0, 1},
                                          Rational(0), Rational(100), 3, 0, 5);
    if (nilpotent_future.status != TailMajorantStatus::Certified ||
        nilpotent_future.exact_resolvent_checks != 0 ||
        nilpotent_future.inverse_n_times_prefix_upper.size() != 2 ||
        !(nilpotent_future.inverse_n_times_prefix_upper.front() <=
          nilpotent_future.inverse_n_times_prefix_upper.back()))
      throw std::runtime_error(
          "causal epsilon/log nilpotents changed future invertibility: " +
          nilpotent_future.detail);

    // A causal epsilon ladder can have a very large ordinary infinity norm
    // even though its diagonal recurrence contracts.  A geometric epsilon
    // weight suppresses only the strictly triangular coefficients.  The
    // retained-prefix induction bound must be rescaled by the inverse weight
    // so the resulting ordinary coefficient enclosure remains rigorous.
    auto weighted_equation = constant_equation(1);
    weighted_equation.c_lags.resize(2);
    weighted_equation.c_lags[1].push_back(
        PhysicalODEMatrixEntry<Rational>{0, 0, rational(0, {"1", "100"})});
    auto weighted_local =
        local_shell(1, 4, "singular-tail-weighted-epsilon-local");
    weighted_local.epsilon = {0, 1};
    auto weighted_sector = zero_sector("0", 0, 1, 4);
    weighted_sector.coefficients.assign(10, Rational(0));
    Rational factorial(1);
    for (std::uint32_t n = 0; n <= 4; ++n) {
      if (n > 0)
        factorial *= Rational(std::to_string(n));
      weighted_sector.coefficients[n] = Rational(1) / factorial;
      if (n > 0)
        weighted_sector.coefficients[5 + n] =
            Rational(100) / (factorial / Rational(std::to_string(n)));
    }
    weighted_local.sectors = {weighted_sector};
    const auto unweighted_transfer = diffexp::kernel::singular_tail_majorant_detail::
        certify_interval_transfer_contraction_on_x_interval(
            weighted_equation, EpsilonWindow{0, 1}, Rational(0), Rational(0), 0,
            0, Rational(0), Rational("1/5"), Rational("1/2"), 1);
    if (!(unweighted_transfer.contraction_upper.back() > Magnitude::one()))
      throw std::runtime_error(
          "weighted epsilon fixture does not exercise noncontracting "
          "ordinary infinity norm");
    const auto degree_bounds = diffexp::kernel::singular_tail_majorant_detail::
        certify_interval_transfer_degree_bounds_on_x_interval(
            weighted_equation, EpsilonWindow{0, 1}, Rational(0), Rational(0), 0,
            0, Rational(0), Rational("1/5"), Rational("1/2"));
    const auto cached_weighted =
        diffexp::kernel::singular_tail_majorant_detail::
            weighted_interval_transfer_contraction(degree_bounds, 16);
    const auto direct_weighted = diffexp::kernel::singular_tail_majorant_detail::
        certify_interval_transfer_contraction_on_x_interval(
            weighted_equation, EpsilonWindow{0, 1}, Rational(0), Rational(0), 0,
            0, Rational(0), Rational("1/5"), Rational("1/2"), 16);
    if (cached_weighted.status != TailMajorantStatus::Certified ||
        cached_weighted.contraction_upper.size() !=
            direct_weighted.contraction_upper.size())
      throw std::runtime_error(
          "weight-independent interval transfer cache lost its prefix");
    for (std::size_t index = 0;
         index < cached_weighted.contraction_upper.size(); ++index)
      if (cached_weighted.contraction_upper[index].dump_exact() !=
          direct_weighted.contraction_upper[index].dump_exact())
        throw std::runtime_error(
            "cached interval transfer weighting changed its enclosure");
    const auto weighted_model =
        diffexp::kernel::prepare_singular_rational_shadow_tail_model(
            weighted_equation, weighted_local, EpsilonWindow{0, 1}, 4, "1/2");
    if (weighted_model.status != TailMajorantStatus::Certified ||
        !weighted_model.model.has_value() ||
        weighted_model.model->towers.front().epsilon_weight_base <= 1 ||
        weighted_model.model->towers.front().contraction_upper.back() >
            Magnitude::one())
      throw std::runtime_error(
          "causal epsilon weighting did not certify the triangular "
          "recurrence: " +
          weighted_model.detail);
    const auto weighted_prefix =
        diffexp::kernel::certify_singular_rational_shadow_prefix(
            weighted_equation, weighted_local, EpsilonWindow{0, 1}, 4);
    if (weighted_prefix.tower_transfers.empty() ||
        weighted_prefix.tower_transfers.front()
                .interval_transfer.status !=
            TailMajorantStatus::Certified)
      throw std::runtime_error(
          "exact Frobenius prefix did not retain its radius-independent "
          "transfer");
    const auto cached_weighted_model =
        diffexp::kernel::prepare_singular_rational_shadow_tail_model(
            weighted_equation, weighted_local, EpsilonWindow{0, 1}, 4, "1/2",
            &weighted_prefix);
    if (cached_weighted_model.status != TailMajorantStatus::Certified ||
        !cached_weighted_model.model.has_value() ||
        cached_weighted_model.model->towers.size() !=
            weighted_model.model->towers.size())
      throw std::runtime_error(
          "cached Frobenius prefix did not reproduce its tail model");
    for (std::size_t tower = 0;
         tower < weighted_model.model->towers.size(); ++tower) {
      const auto &direct = weighted_model.model->towers[tower];
      const auto &cached = cached_weighted_model.model->towers[tower];
      if (direct.epsilon_weight_base != cached.epsilon_weight_base ||
          direct.contraction_upper.size() !=
              cached.contraction_upper.size())
        throw std::runtime_error(
            "radius-independent Frobenius transfer cache changed its weight");
      for (std::size_t index = 0;
           index < direct.contraction_upper.size(); ++index)
        if (direct.contraction_upper[index].dump_exact() !=
            cached.contraction_upper[index].dump_exact())
          throw std::runtime_error(
              "radius-independent Frobenius transfer cache changed its "
              "contraction enclosure");
    }

    // Material CASE-P b towers stay independent, while exact-zero sectors
    // do not create towers or inflate the retained log ceiling.
    auto casep_equation = constant_equation(2);
    casep_equation.c_lags[0].push_back(
        PhysicalODEMatrixEntry<Rational>{0, 0, rational(1, {"2"})});
    casep_equation.c_lags[0].push_back(
        PhysicalODEMatrixEntry<Rational>{1, 1, rational(1, {"3"})});
    auto casep = local_shell(2, 0, "singular-tail-casep-material-local");
    casep.epsilon = {0, 1};
    auto b2 = zero_sector("0", 0, 2, 0);
    b2.b = diffexp::kernel::ExactScalarDescriptor::rational("2");
    b2.coefficients.assign(4, Rational(0));
    b2.coefficients[0] = Rational(1);
    auto b3 = zero_sector("0", 0, 2, 0);
    b3.b = diffexp::kernel::ExactScalarDescriptor::rational("3");
    b3.coefficients.assign(4, Rational(0));
    b3.coefficients[1] = Rational(1);
    auto zero_high_log = zero_sector("0", 7, 2, 0);
    zero_high_log.b = diffexp::kernel::ExactScalarDescriptor::rational("99");
    zero_high_log.coefficients.assign(4, Rational(0));
    casep.sectors = {b2, b3, zero_high_log};
    const auto casep_model =
        diffexp::kernel::prepare_singular_rational_shadow_tail_model(
            casep_equation, casep, EpsilonWindow{0, 1}, 0, "1/2");
    if (casep_model.status != TailMajorantStatus::Certified ||
        !casep_model.model.has_value() ||
        casep_model.model->towers.size() != 2 ||
        std::any_of(casep_model.model->towers.begin(),
                    casep_model.model->towers.end(), [](const auto &tower) {
                      return tower.log_complete_max != 0 ||
                             (tower.b_exact != "2" && tower.b_exact != "3");
                    }))
      throw std::runtime_error(
          "material/zero CASE-P tower pruning is incorrect: " +
          casep_model.detail);

    std::cout << "singular Frobenius tail-majorant tests passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
