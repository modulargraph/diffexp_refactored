#include "diffexp2/local_algebra.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {

using diffexp2::ExactScalarDescriptor;
using diffexp2::ComplexBall;
using diffexp2::EpsilonFrame;
using diffexp2::ExactLaurentMatrix;
using diffexp2::ExactLaurentPolynomial;
using diffexp2::FiniteLaurentVector;
using diffexp2::LocalSector;
using diffexp2::LocalSolution;
using diffexp2::PreparedRationalAnalyticCoefficient;
using diffexp2::PreparedRationalTaylorMultiplier;
using diffexp2::PreparedSparseLocalMultiplierMatrix;
using diffexp2::Rational;
using diffexp2::RealEvaluationPoint;

std::size_t index(std::size_t ei, std::size_t n, std::uint32_t c,
                  std::size_t tw, std::uint32_t d) {
  return ((ei * tw + n) * d) + c;
}

LocalSolution<Rational> sample() {
  LocalSolution<Rational> out;
  out.chart.center_exact = "0";
  out.chart.scale_exact = "1";
  out.chart.radius = diffexp2::ComplexBall::from_strings("2");
  out.epsilon = {-1, 1};
  out.taylor_complete_max = 2;
  out.dimension = 2;
  out.checkpoint_identity = "local-algebra-input";
  out.prescriptions.push_back({"x", 1, 1, 1});
  LocalSector<Rational> sector;
  sector.a = ExactScalarDescriptor::rational("1/2");
  sector.b = ExactScalarDescriptor::rational("0");
  sector.log_power = 0;
  sector.coefficients.assign(out.sector_size(), Rational(0));
  for (std::size_t ei = 0; ei < out.epsilon.width(); ++ei)
    for (std::size_t n = 0; n < out.taylor_width(); ++n)
      for (std::uint32_t c = 0; c < out.dimension; ++c)
        sector.coefficients[index(ei, n, c, out.taylor_width(), out.dimension)] =
            Rational(std::to_string(100 * ei + 10 * n + c + 1));
  out.sectors.push_back(std::move(sector));
  diffexp2::validate_local_solution(out, false);
  return out;
}

PreparedRationalTaylorMultiplier<Rational> multiplier(
    std::int32_t shift, std::uint32_t pole) {
  PreparedRationalTaylorMultiplier<Rational> out;
  out.epsilon_shift = shift;
  out.center_pole_order = pole;
  out.exact_identity = "(1+2t)+3eps";
  out.kernels = {
      {Rational(1), Rational(2), Rational(0)},
      {Rational(3), Rational(0), Rational(0)},
      {Rational(0), Rational(0), Rational(0)}};
  return out;
}

bool equal(const Rational& value, const std::string& expected) {
  return value == Rational(expected);
}

bool same_rational_local(const LocalSolution<Rational>& left,
                         const LocalSolution<Rational>& right) {
  if (left.chart.center_exact != right.chart.center_exact ||
      left.chart.scale_exact != right.chart.scale_exact ||
      left.epsilon.min_power != right.epsilon.min_power ||
      left.epsilon.complete_max != right.epsilon.complete_max ||
      left.taylor_complete_max != right.taylor_complete_max ||
      left.dimension != right.dimension)
    return false;
  std::vector<const LocalSector<Rational>*> left_material;
  std::vector<const LocalSector<Rational>*> right_material;
  for (const auto& sector : left.sectors)
    if (std::any_of(sector.coefficients.begin(), sector.coefficients.end(),
                    [](const auto& value) { return !value.is_zero(); }))
      left_material.push_back(&sector);
  for (const auto& sector : right.sectors)
    if (std::any_of(sector.coefficients.begin(), sector.coefficients.end(),
                    [](const auto& value) { return !value.is_zero(); }))
      right_material.push_back(&sector);
  if (left_material.size() != right_material.size()) return false;
  for (std::size_t sector = 0; sector < left_material.size(); ++sector) {
    const auto& a = *left_material[sector];
    const auto& b = *right_material[sector];
    if (a.a.canonical != b.a.canonical ||
        a.b.canonical != b.b.canonical ||
        a.log_power != b.log_power ||
        a.coefficients != b.coefficients)
      return false;
  }
  return true;
}

LocalSolution<ComplexBall> acb_sample() {
  const auto exact = sample();
  LocalSolution<ComplexBall> out;
  out.chart = exact.chart;
  out.epsilon = exact.epsilon;
  out.taylor_complete_max = exact.taylor_complete_max;
  out.dimension = exact.dimension;
  out.prescriptions = exact.prescriptions;
  out.checkpoint_identity = "local-algebra-acb-input";
  for (const auto& source : exact.sectors) {
    LocalSector<ComplexBall> sector;
    sector.a = source.a;
    sector.b = source.b;
    sector.log_power = source.log_power;
    for (const auto& coefficient : source.coefficients)
      sector.coefficients.push_back(
          ComplexBall::from_strings(coefficient.str()));
    out.sectors.push_back(std::move(sector));
  }
  diffexp2::validate_local_solution(out, false);
  return out;
}

PreparedRationalTaylorMultiplier<ComplexBall> acb_multiplier(
    std::int32_t shift, std::uint32_t pole) {
  const auto exact = multiplier(shift, pole);
  PreparedRationalTaylorMultiplier<ComplexBall> out;
  out.epsilon_shift = exact.epsilon_shift;
  out.center_pole_order = exact.center_pole_order;
  out.exact_identity = exact.exact_identity;
  for (const auto& exact_kernel : exact.kernels) {
    std::vector<ComplexBall> kernel;
    for (const auto& coefficient : exact_kernel)
      kernel.push_back(ComplexBall::from_strings(coefficient.str()));
    out.kernels.push_back(std::move(kernel));
  }
  return out;
}

bool overlapping_acb_local(const LocalSolution<ComplexBall>& left,
                           const LocalSolution<ComplexBall>& right) {
  if (left.epsilon.min_power != right.epsilon.min_power ||
      left.epsilon.complete_max != right.epsilon.complete_max ||
      left.taylor_complete_max != right.taylor_complete_max ||
      left.dimension != right.dimension)
    return false;
  std::vector<const LocalSector<ComplexBall>*> left_material;
  std::vector<const LocalSector<ComplexBall>*> right_material;
  for (const auto& sector : left.sectors)
    if (std::any_of(sector.coefficients.begin(), sector.coefficients.end(),
                    [](const auto& value) { return !value.is_zero(); }))
      left_material.push_back(&sector);
  for (const auto& sector : right.sectors)
    if (std::any_of(sector.coefficients.begin(), sector.coefficients.end(),
                    [](const auto& value) { return !value.is_zero(); }))
      right_material.push_back(&sector);
  if (left_material.size() != right_material.size()) return false;
  for (std::size_t sector = 0; sector < left_material.size(); ++sector) {
    const auto& a = *left_material[sector];
    const auto& b = *right_material[sector];
    if (a.a.canonical != b.a.canonical ||
        a.b.canonical != b.b.canonical ||
        a.log_power != b.log_power ||
        a.coefficients.size() != b.coefficients.size())
      return false;
    for (std::size_t coefficient = 0;
         coefficient < a.coefficients.size(); ++coefficient)
      if (!(a.coefficients[coefficient] - b.coefficients[coefficient])
               .contains_zero())
        return false;
  }
  return true;
}

PreparedRationalTaylorMultiplier<Rational> identity_multiplier() {
  auto out = multiplier(0, 0);
  for (auto& kernel : out.kernels)
    std::fill(kernel.begin(), kernel.end(), Rational(0));
  out.kernels.front().front() = Rational(1);
  out.exact_identity = "1";
  return out;
}

LocalSolution<Rational> component_sparse_sample(std::size_t sector_count) {
  auto out = sample();
  out.checkpoint_identity = "component-sparse-input";
  out.sectors.clear();
  out.sectors.reserve(sector_count);
  for (std::size_t sector_index = 0; sector_index < sector_count;
       ++sector_index) {
    LocalSector<Rational> sector;
    sector.a = ExactScalarDescriptor::rational(
        Rational(std::to_string(2 * sector_index + 1) + "/2").str());
    sector.b = ExactScalarDescriptor::rational("0");
    sector.log_power = 0;
    sector.coefficients.assign(out.sector_size(), Rational(0));
    // Component zero is material everywhere, but the projected component
    // one occurs in only the final sector.
    sector.coefficients[index(0, 0, 0, out.taylor_width(), out.dimension)] =
        Rational(1);
    if (sector_index + 1 == sector_count)
      sector.coefficients[index(
          0, 0, 1, out.taylor_width(), out.dimension)] = Rational(7);
    out.sectors.push_back(std::move(sector));
  }
  diffexp2::validate_local_solution(out, false);
  return out;
}

PreparedRationalTaylorMultiplier<Rational> casep_constant_multiplier() {
  PreparedRationalTaylorMultiplier<Rational> out;
  out.epsilon_shift = -2;
  out.center_pole_order = 0;
  out.exact_identity = "casep-taylor-constant";
  out.kernels.assign(3, std::vector<Rational>(3, Rational(0)));
  out.kernels[0][0] = Rational("2/3");
  out.kernels[2][0] = Rational("-5/7");
  return out;
}

}  // namespace

int main() {
  diffexp2::ComplexBall::set_precision(128);
  const auto input = sample();
  const auto matched_prefix =
      diffexp2::restrict_local_taylor_prefix(
          input, 1, "matched-prefix");
  const auto product = diffexp2::multiply_prepared_rational(
      input, multiplier(-1, 1));

  bool ok = matched_prefix.taylor_complete_max == 1 &&
            matched_prefix.epsilon.min_power ==
                input.epsilon.min_power &&
            matched_prefix.epsilon.complete_max ==
                input.epsilon.complete_max &&
            matched_prefix.dimension == input.dimension &&
            matched_prefix.prescriptions.size() ==
                input.prescriptions.size() &&
            matched_prefix.prescriptions.front().factor_exact ==
                input.prescriptions.front().factor_exact &&
            matched_prefix.checkpoint_identity == "matched-prefix" &&
            equal(matched_prefix.sectors.front().coefficients[
                index(2, 1, 1, 2, 2)], "212") &&
            product.epsilon.min_power == -2 &&
            product.epsilon.complete_max == 0 &&
            product.taylor_complete_max == 2 &&
            product.dimension == 2 && product.sectors.size() == 1 &&
            product.sectors.front().a.canonical == "-1/2" &&
            product.prescriptions.size() == 1;
  const auto& coefficients = product.sectors.front().coefficients;
  // out(eps index 0,t^1,c0) = in(0,t^1,c0) + 2 in(0,t^0,c0)
  ok = ok && equal(coefficients[index(0, 1, 0, 3, 2)], "13");
  // out(eps index 1,t^0,c0) = in(1,t^0,c0) + 3 in(0,t^0,c0)
  ok = ok && equal(coefficients[index(1, 0, 0, 3, 2)], "104");
  // out(eps index 2,t^2,c1) = in(2,t^2,c1) + 2 in(2,t^1,c1)
  //                              + 3 in(1,t^2,c1)
  ok = ok && equal(coefficients[index(2, 2, 1, 3, 2)], "1012");

  // CASE-P polar weights have no positive Taylor coefficients.  Exercise the
  // exact fast path with a negative epsilon shift, a sparse epsilon kernel,
  // multiple Taylor orders/components, and distinct Frobenius sectors.
  auto casep_input = sample();
  auto second_sector = casep_input.sectors.front();
  second_sector.a = ExactScalarDescriptor::rational("3/2");
  second_sector.b = ExactScalarDescriptor::rational("1/2");
  second_sector.log_power = 2;
  for (std::size_t coefficient = 0;
       coefficient < second_sector.coefficients.size(); ++coefficient)
    second_sector.coefficients[coefficient] = Rational(
        "-" + std::to_string(coefficient + 1) + "/11");
  casep_input.sectors.push_back(std::move(second_sector));
  const auto casep_multiplier = casep_constant_multiplier();
  const auto casep_product = diffexp2::multiply_prepared_rational(
      casep_input, casep_multiplier, "casep-constant-product");
  bool casep_constant_ok =
      casep_product.checkpoint_identity == "casep-constant-product" &&
      casep_product.epsilon.min_power == -3 &&
      casep_product.epsilon.complete_max == -1 &&
      casep_product.taylor_complete_max ==
          casep_input.taylor_complete_max &&
      casep_product.dimension == casep_input.dimension &&
      casep_product.prescriptions.size() == 1 &&
      casep_product.prescriptions.front().factor_exact == "x" &&
      casep_product.prescriptions.front().sign == 1 &&
      casep_product.prescriptions.front().multiplicity == 1 &&
      casep_product.prescriptions.front().leading_coefficient_sign == 1 &&
      casep_product.sectors.size() == casep_input.sectors.size();
  for (std::size_t sector = 0;
       sector < casep_input.sectors.size() && casep_constant_ok; ++sector) {
    const auto& source = casep_input.sectors[sector];
    const auto& actual = casep_product.sectors[sector];
    casep_constant_ok = actual.a.canonical == source.a.canonical &&
        actual.b.canonical == source.b.canonical &&
        actual.log_power == source.log_power;
    for (std::size_t out_ei = 0;
         out_ei < casep_input.epsilon.width() && casep_constant_ok; ++out_ei) {
      for (std::size_t n = 0;
           n < casep_input.taylor_width() && casep_constant_ok; ++n) {
        for (std::uint32_t component = 0;
             component < casep_input.dimension; ++component) {
          Rational expected(0);
          for (std::size_t kernel_ei = 0; kernel_ei <= out_ei;
               ++kernel_ei) {
            expected += casep_multiplier.kernels[kernel_ei][0] *
                source.coefficients[index(
                    out_ei - kernel_ei, n, component,
                    casep_input.taylor_width(), casep_input.dimension)];
          }
          if (actual.coefficients[index(
                  out_ei, n, component, casep_input.taylor_width(),
                  casep_input.dimension)] != expected) {
            casep_constant_ok = false;
            break;
          }
        }
      }
    }
  }
  ok = ok && casep_constant_ok;

  PreparedSparseLocalMultiplierMatrix<Rational> matrix;
  matrix.rows = 2;
  matrix.columns = 2;
  matrix.exact_identity = "one-coupling";
  matrix.entries.push_back({1, 0, multiplier(-1, 1)});
  const auto applied = diffexp2::apply_prepared_sparse_local_matrix(
      matrix, input);
  ok = ok && applied.has_value() && applied->dimension == 2 &&
       applied->epsilon.min_power == -2 &&
       applied->epsilon.complete_max == 0;
  if (applied.has_value()) {
    const auto& slab = applied->sectors.front().coefficients;
    ok = ok && equal(slab[index(1, 0, 0, 3, 2)], "0") &&
         equal(slab[index(1, 0, 1, 3, 2)], "104");
  }

  const auto unshifted = diffexp2::multiply_prepared_rational(
      input, multiplier(0, 0));
  const auto combined = diffexp2::combine_local_solutions(
      std::vector<LocalSolution<Rational>>{product, unshifted});
  ok = ok && combined.epsilon.min_power == -2 &&
       combined.epsilon.complete_max == 0;

  // A saturated solve retains a wider transformed-weight reservoir than its
  // published physical weights.  Materializing F*(T*y) after truncating the
  // physical weights can therefore lose a coefficient which a polar basis
  // needs at the public edge.  The exact-right association (F*T)*y keeps it.
  const auto exact_right_column = [](
      const std::string& checkpoint, const Rational& polar_value) {
    LocalSolution<Rational> local;
    local.chart.center_exact = "0";
    local.chart.scale_exact = "1";
    local.chart.radius = ComplexBall(2);
    local.epsilon = {-1, 1};
    local.taylor_complete_max = 0;
    local.dimension = 1;
    local.checkpoint_identity = checkpoint;
    LocalSector<Rational> sector;
    sector.a = ExactScalarDescriptor::rational("0");
    sector.b = ExactScalarDescriptor::rational("0");
    sector.log_power = 0;
    sector.coefficients.assign(local.sector_size(), Rational(0));
    sector.coefficients.front() = polar_value;
    local.sectors.push_back(std::move(sector));
    diffexp2::validate_local_solution(local, false);
    return local;
  };
  const auto exact_right_zero =
      exact_right_column("exact-right-zero", Rational(0));
  const auto exact_right_polar =
      exact_right_column("exact-right-polar", Rational(7));
  const std::vector<const LocalSolution<Rational>*> exact_right_basis{
      &exact_right_zero, &exact_right_polar};
  ExactLaurentMatrix<Rational> exact_right_transformation(
      2, std::vector<ExactLaurentPolynomial<Rational>>(2));
  exact_right_transformation[0][0].add_term(-1, Rational(1));
  exact_right_transformation[1][0].add_term(0, Rational(1));
  exact_right_transformation[1][1].add_term(0, Rational(1));
  const auto exact_right_transformed =
      diffexp2::right_transform_local_basis_exact(
          exact_right_basis, exact_right_transformation,
          "exact-right-transformed");
  std::vector<const LocalSolution<Rational>*> transformed_basis;
  for (const auto& column : exact_right_transformed)
    transformed_basis.push_back(&column);
  const FiniteLaurentVector<Rational> transformed_weights{
      EpsilonFrame<Rational>({1, 1}, {Rational(1)}),
      EpsilonFrame<Rational>({1, 1}, {Rational(0)})};
  const auto stable_exact_right =
      diffexp2::materialize_local_basis_weights(
          transformed_basis, transformed_weights,
          "exact-right-stable-output");
  const FiniteLaurentVector<Rational> truncated_physical_weights{
      EpsilonFrame<Rational>({0, 0}, {Rational(1)}),
      EpsilonFrame<Rational>({0, 0}, {Rational(0)})};
  const auto truncated_physical =
      diffexp2::materialize_local_basis_weights(
          exact_right_basis, truncated_physical_weights,
          "exact-right-truncated-physical-output");
  ok = ok &&
       stable_exact_right.epsilon.min_power <= 0 &&
       stable_exact_right.epsilon.complete_max >= 0 &&
       stable_exact_right.sectors.size() == 1 &&
       equal(stable_exact_right.sectors.front().coefficients[
                 static_cast<std::size_t>(
                     -stable_exact_right.epsilon.min_power)],
             "7") &&
       truncated_physical.sectors.size() == 1 &&
       std::all_of(
           truncated_physical.sectors.front().coefficients.begin(),
           truncated_physical.sectors.front().coefficients.end(),
           [](const auto& value) { return value.is_zero(); });

  // A finite denormalization frame R may agree with its exact inverse only
  // through a bounded epsilon prefix.  It cannot replace the exact physical
  // saturation T which was certified directly against F.  In this scalar
  // example R=(1+eps), R^-1=(1-eps+...), so the retained finite product is
  // complete only through eps^1.  Contracting it with an eps^-2 weight loses
  // the public eps^0 coefficient, while direct exact F*T retains it.
  LocalSolution<Rational> direct_t_source;
  direct_t_source.chart.center_exact = "0";
  direct_t_source.chart.scale_exact = "1";
  direct_t_source.chart.radius = ComplexBall(2);
  direct_t_source.epsilon = {0, 3};
  direct_t_source.taylor_complete_max = 0;
  direct_t_source.dimension = 1;
  direct_t_source.checkpoint_identity = "direct-physical-T-source";
  LocalSector<Rational> direct_t_sector;
  direct_t_sector.a = ExactScalarDescriptor::rational("0");
  direct_t_sector.b = ExactScalarDescriptor::rational("0");
  direct_t_sector.log_power = 0;
  direct_t_sector.coefficients = {
      Rational(1), Rational(2), Rational(3), Rational(4)};
  direct_t_source.sectors.push_back(std::move(direct_t_sector));
  diffexp2::validate_local_solution(direct_t_source, false);
  const std::vector<const LocalSolution<Rational>*> direct_t_basis{
      &direct_t_source};
  ExactLaurentMatrix<Rational> direct_t(
      1, std::vector<ExactLaurentPolynomial<Rational>>(1));
  direct_t[0][0].add_term(0, Rational(1));
  const auto direct_t_columns =
      diffexp2::right_transform_local_basis_exact(
          direct_t_basis, direct_t, "direct-physical-F*T");
  diffexp2::FiniteLaurentMatrix<Rational> bounded_r{
      {EpsilonFrame<Rational>(
          {0, 1}, {Rational(1), Rational(1)})}};
  ExactLaurentMatrix<Rational> exact_r_inverse(
      1, std::vector<ExactLaurentPolynomial<Rational>>(1));
  exact_r_inverse[0][0].add_term(0, Rational(1));
  exact_r_inverse[0][0].add_term(1, Rational(-1));
  const auto bounded_identity =
      diffexp2::right_multiply_finite_by_exact_laurent(
          bounded_r, exact_r_inverse);
  const auto bounded_t_column =
      diffexp2::materialize_local_basis_weights(
          direct_t_basis, {bounded_identity[0][0]},
          "bounded-R-times-R-inverse-times-T");
  const EpsilonFrame<Rational> polar_terminal_weight(
      {-2, 0}, {Rational(1), Rational(0), Rational(0)});
  const std::vector<const LocalSolution<Rational>*> direct_t_columns_view{
      &direct_t_columns.front()};
  const std::vector<const LocalSolution<Rational>*> bounded_t_columns_view{
      &bounded_t_column};
  const auto direct_t_result =
      diffexp2::materialize_local_basis_weights(
          direct_t_columns_view, {polar_terminal_weight},
          "direct-physical-F*T-result");
  const auto bounded_t_result =
      diffexp2::materialize_local_basis_weights(
          bounded_t_columns_view, {polar_terminal_weight},
          "bounded-frame-result");
  ok = ok &&
       direct_t_result.epsilon.complete_max >= 0 &&
       equal(direct_t_result.sectors.front().coefficients[
                 static_cast<std::size_t>(
                     -direct_t_result.epsilon.min_power)],
             "3") &&
       bounded_t_result.epsilon.complete_max == -1;

  // A regular physical value can require more than one private Laurent row
  // when distinct Frobenius exponents coalesce.  Here
  //
  //   eps^-2 (1 - 2 t^eps + t^(2 eps)) = log(t)^2 + O(eps).
  //
  // Exercise the exact-right and materialization path, not only the local
  // evaluator: dropping either pole row before evaluation destroys the
  // finite value and its theta derivative.
  const auto confluent_column = [](
      const std::string& checkpoint, const std::string& b,
      const std::string& polar_coefficient) {
    LocalSolution<ComplexBall> local;
    local.chart.center_exact = "0";
    local.chart.scale_exact = "1";
    local.chart.radius = ComplexBall(2);
    local.epsilon = {-2, 2};
    local.taylor_complete_max = 0;
    local.dimension = 1;
    local.checkpoint_identity = checkpoint;
    LocalSector<ComplexBall> sector;
    sector.a = ExactScalarDescriptor::rational("0");
    sector.b = ExactScalarDescriptor::rational(b);
    sector.log_power = 0;
    sector.coefficients.assign(local.sector_size(), ComplexBall(0));
    sector.coefficients.front() =
        ComplexBall::from_strings(polar_coefficient);
    local.sectors.push_back(std::move(sector));
    diffexp2::validate_local_solution(local, false);
    return local;
  };
  const auto confluent_zero =
      confluent_column("confluent-zero", "0", "1");
  const auto confluent_one =
      confluent_column("confluent-one", "1", "-2");
  const auto confluent_two =
      confluent_column("confluent-two", "2", "1");
  const std::vector<const LocalSolution<ComplexBall>*> confluent_basis{
      &confluent_zero, &confluent_one, &confluent_two};
  ExactLaurentMatrix<ComplexBall> confluent_identity(
      3, std::vector<ExactLaurentPolynomial<ComplexBall>>(3));
  for (std::size_t diagonal = 0; diagonal < 3; ++diagonal)
    confluent_identity[diagonal][diagonal].add_term(
        0, ComplexBall(1));
  const auto confluent_transformed =
      diffexp2::right_transform_local_basis_exact(
          confluent_basis, confluent_identity,
          "confluent-exact-right");
  std::vector<const LocalSolution<ComplexBall>*>
      confluent_transformed_basis;
  for (const auto& column : confluent_transformed)
    confluent_transformed_basis.push_back(&column);
  const auto complete_constant_weight = [] {
    return EpsilonFrame<ComplexBall>(
        0, {ComplexBall(1), ComplexBall(0), ComplexBall(0),
            ComplexBall(0), ComplexBall(0)});
  };
  const FiniteLaurentVector<ComplexBall> confluent_weights{
      complete_constant_weight(), complete_constant_weight(),
      complete_constant_weight()};
  const auto confluent_materialized =
      diffexp2::materialize_local_basis_weights(
          confluent_transformed_basis, confluent_weights,
          "confluent-materialized");
  const auto confluent_evaluation =
      diffexp2::evaluate_local_solution(
          confluent_materialized,
          RealEvaluationPoint::rational("1/2"));
  ComplexBall log_half;
  acb_log(log_half.raw(), ComplexBall::from_strings("1/2").raw(),
          ComplexBall::precision());
  bool confluent_lower_restriction_rejected = false;
  try {
    (void)diffexp2::restrict_local_epsilon_frame_strict_lower(
        confluent_materialized, 0,
        confluent_materialized.epsilon.complete_max,
        "invalid-confluent-lower-restriction");
  } catch (const std::invalid_argument&) {
    confluent_lower_restriction_rejected = true;
  }
  ok = ok &&
       confluent_evaluation.value.at(-2, 0).contains_zero() &&
       confluent_evaluation.value.at(-1, 0).contains_zero() &&
       (confluent_evaluation.value.at(0, 0) -
        log_half * log_half).contains_zero() &&
       (confluent_evaluation.theta_value.at(0, 0) -
        ComplexBall(2) * log_half).contains_zero() &&
       confluent_lower_restriction_rejected;

  PreparedSparseLocalMultiplierMatrix<Rational> empty;
  empty.rows = 2;
  empty.columns = 2;
  ok = ok && !diffexp2::apply_prepared_sparse_local_matrix(
      empty, input).has_value();

  PreparedSparseLocalMultiplierMatrix<Rational> scalar_row;
  scalar_row.rows = 1;
  scalar_row.columns = 2;
  scalar_row.exact_identity = "direct-scalar-row";
  scalar_row.entries.push_back({0, 0, multiplier(-1, 1)});
  scalar_row.entries.push_back({0, 1, multiplier(0, 0)});
  const auto scalar_full = diffexp2::apply_prepared_sparse_local_matrix(
      scalar_row, input, "scalar-full");
  const auto scalar_direct = diffexp2::apply_prepared_scalar_row_window(
      scalar_row, input, -1, "scalar-direct");
  ok = ok && scalar_full.has_value() && scalar_direct.has_value();
  if (scalar_full.has_value() && scalar_direct.has_value()) {
    const auto scalar_restricted =
        diffexp2::restrict_local_epsilon_frame_strict_lower(
            *scalar_full, scalar_full->epsilon.min_power, -1,
            "scalar-restricted");
    ok = ok && same_rational_local(
                   scalar_restricted, *scalar_direct) &&
         scalar_direct->epsilon.min_power == -2 &&
         scalar_direct->epsilon.complete_max == -1 &&
         scalar_direct->dimension == 1;
  }

  PreparedSparseLocalMultiplierMatrix<Rational> zero_guard_row;
  zero_guard_row.rows = 1;
  zero_guard_row.columns = 2;
  zero_guard_row.exact_identity = "direct-zero-frame-guard";
  zero_guard_row.entries.push_back({0, 0, multiplier(-1, 1)});
  auto zero_active = multiplier(-2, 0);
  for (auto& kernel : zero_active.kernels)
    std::fill(kernel.begin(), kernel.end(), Rational(0));
  zero_guard_row.entries.push_back({0, 1, std::move(zero_active)});
  const auto zero_guard_full =
      diffexp2::apply_prepared_sparse_local_matrix(zero_guard_row, input);
  const auto zero_guard_direct =
      diffexp2::apply_prepared_scalar_row_window(
          zero_guard_row, input, 0);
  ok = ok && zero_guard_full.has_value() &&
       zero_guard_direct.has_value() &&
       same_rational_local(*zero_guard_full, *zero_guard_direct) &&
       zero_guard_direct->epsilon.min_power == -3 &&
       zero_guard_direct->epsilon.complete_max == -1;

  PreparedSparseLocalMultiplierMatrix<Rational> all_zero_row;
  all_zero_row.rows = 1;
  all_zero_row.columns = 2;
  all_zero_row.exact_identity = "active-all-zero-row";
  auto all_zero_multiplier = multiplier(-2, 1);
  for (auto& kernel : all_zero_multiplier.kernels)
    std::fill(kernel.begin(), kernel.end(), Rational(0));
  all_zero_row.entries.push_back(
      {0, 1, std::move(all_zero_multiplier)});
  const auto all_zero_direct =
      diffexp2::apply_prepared_scalar_row_window(
          all_zero_row, input, 0, "active-all-zero-direct");
  ok = ok && all_zero_direct.has_value() &&
       all_zero_direct->epsilon.min_power == -3 &&
       all_zero_direct->epsilon.complete_max == -1 &&
       all_zero_direct->sectors.size() == 1 &&
       all_zero_direct->sectors.front().a.canonical == "-1/2" &&
       std::all_of(
           all_zero_direct->sectors.front().coefficients.begin(),
           all_zero_direct->sectors.front().coefficients.end(),
           [](const auto& value) { return value.is_zero(); });

  const auto sparse_input = component_sparse_sample(128);
  PreparedSparseLocalMultiplierMatrix<Rational> sparse_row;
  sparse_row.rows = 1;
  sparse_row.columns = 2;
  sparse_row.exact_identity = "component-sparse-row";
  sparse_row.entries.push_back({0, 1, identity_multiplier()});
  const auto sparse_full = diffexp2::apply_prepared_sparse_local_matrix(
      sparse_row, sparse_input, "component-sparse-full");
  const auto sparse_direct =
      diffexp2::apply_prepared_scalar_row_window(
          sparse_row, sparse_input, sparse_input.epsilon.complete_max,
          "component-sparse-direct");
  ok = ok && sparse_full.has_value() && sparse_direct.has_value() &&
       same_rational_local(*sparse_full, *sparse_direct) &&
       sparse_full->sectors.size() == 128 &&
       sparse_direct->sectors.size() == 1 &&
       sparse_direct->sectors.front().a.canonical == "255/2";

  PreparedSparseLocalMultiplierMatrix<ComplexBall> acb_row;
  acb_row.rows = 1;
  acb_row.columns = 2;
  acb_row.exact_identity = "direct-acb-row";
  acb_row.entries.push_back({0, 0, acb_multiplier(-1, 1)});
  const auto acb_input = acb_sample();
  const auto acb_full = diffexp2::apply_prepared_sparse_local_matrix(
      acb_row, acb_input, "acb-full");
  const auto acb_direct = diffexp2::apply_prepared_scalar_row_window(
      acb_row, acb_input, -1, "acb-direct");
  ok = ok && acb_full.has_value() && acb_direct.has_value();
  if (acb_full.has_value() && acb_direct.has_value()) {
    const auto acb_restricted =
        diffexp2::restrict_local_epsilon_frame_strict_lower(
            *acb_full, acb_full->epsilon.min_power, -1,
            "acb-restricted");
    ok = ok && overlapping_acb_local(acb_restricted, *acb_direct);
  }

  auto ambiguous_input = acb_sample();
  for (auto& sector : ambiguous_input.sectors)
    for (auto& coefficient : sector.coefficients)
      coefficient = ComplexBall(0);
  auto ambiguous = ComplexBall(0);
  arb_add_error_2exp_si(acb_realref(ambiguous.raw()), -80);
  ambiguous_input.sectors.front().coefficients[index(
      0, 0, 0, ambiguous_input.taylor_width(),
      ambiguous_input.dimension)] = ambiguous;
  PreparedSparseLocalMultiplierMatrix<ComplexBall> ambiguous_row;
  ambiguous_row.rows = 1;
  ambiguous_row.columns = 2;
  ambiguous_row.exact_identity = "ambiguous-acb-row";
  auto ambiguous_multiplier = acb_multiplier(0, 0);
  for (auto& kernel : ambiguous_multiplier.kernels)
    for (auto& coefficient : kernel)
      coefficient = ComplexBall(0);
  ambiguous_multiplier.kernels.front().front() = ComplexBall(1);
  ambiguous_row.entries.push_back(
      {0, 0, std::move(ambiguous_multiplier)});
  const auto ambiguous_direct =
      diffexp2::apply_prepared_scalar_row_window(
          ambiguous_row, ambiguous_input,
          ambiguous_input.epsilon.complete_max, "ambiguous-acb-direct");
  ok = ok && ambiguous_direct.has_value() &&
       ambiguous_direct->sectors.size() == 1 &&
       !ambiguous_direct->sectors.front().coefficients.front().is_zero() &&
       ambiguous_direct->sectors.front().coefficients.front().contains_zero();

  PreparedSparseLocalMultiplierMatrix<Rational> empty_scalar;
  empty_scalar.rows = 1;
  empty_scalar.columns = 2;
  ok = ok && !diffexp2::apply_prepared_scalar_row_window(
      empty_scalar, input, 0).has_value();

  PreparedRationalTaylorMultiplier<ComplexBall> exact_point_multiplier;
  exact_point_multiplier.epsilon_shift = -2;
  exact_point_multiplier.center_pole_order = 2;
  exact_point_multiplier.exact_identity =
      "eps^-2 t^-2 ((1+t)/(1-t) + eps (2-t)/(1+t))";
  // Deliberately unrelated finite kernels: point specialization must consume
  // the retained analytic rationals, never this center Taylor prefix.
  exact_point_multiplier.kernels = {
      {ComplexBall::from_strings("999")},
      {ComplexBall::from_strings("-999")}};
  exact_point_multiplier.analytic_coefficients =
      std::vector<PreparedRationalAnalyticCoefficient<ComplexBall>>{
          {{ComplexBall::from_strings("1"),
            ComplexBall::from_strings("1")},
           {ComplexBall::from_strings("1"),
            ComplexBall::from_strings("-1")}},
          {{ComplexBall::from_strings("2"),
            ComplexBall::from_strings("-1")},
           {ComplexBall::from_strings("1"),
            ComplexBall::from_strings("1")}}};
  const auto specialized_negative =
      diffexp2::local_algebra_detail::
          specialize_prepared_rational_multiplier_at_real_point(
              exact_point_multiplier,
              RealEvaluationPoint::rational("-1/2"));
  const auto specialized_positive =
      diffexp2::local_algebra_detail::
          specialize_prepared_rational_multiplier_at_real_point(
              exact_point_multiplier,
              RealEvaluationPoint::rational("1/2"));
  ok = ok && specialized_negative.has_value() &&
       specialized_negative->min_power() == -2 &&
       specialized_negative->complete_max() == -1 &&
       (specialized_negative->coefficient(-2) -
        ComplexBall::from_strings("4/3")).contains_zero() &&
       (specialized_negative->coefficient(-1) -
        ComplexBall::from_strings("20")).contains_zero() &&
       specialized_positive.has_value() &&
       (specialized_positive->coefficient(-2) -
        ComplexBall::from_strings("12")).contains_zero();

  auto missing_analytic = exact_point_multiplier;
  missing_analytic.analytic_coefficients.reset();
  const auto missing_analytic_diagnostic =
      diffexp2::local_algebra_detail::
          diagnose_prepared_rational_specialization_at_real_point(
              missing_analytic,
              RealEvaluationPoint::rational("-1/2"));
  ok = ok &&
       !diffexp2::local_algebra_detail::
            specialize_prepared_rational_multiplier_at_real_point(
                missing_analytic,
                RealEvaluationPoint::rational("-1/2")).has_value() &&
       missing_analytic_diagnostic.has_value() &&
       missing_analytic_diagnostic->find("unavailable") !=
           std::string::npos;
  auto singular_denominator = exact_point_multiplier;
  singular_denominator.analytic_coefficients->front().denominator = {
      ComplexBall::from_strings("1"),
      ComplexBall::from_strings("2")};
  const auto singular_denominator_diagnostic =
      diffexp2::local_algebra_detail::
          diagnose_prepared_rational_specialization_at_real_point(
              singular_denominator,
              RealEvaluationPoint::rational("-1/2"));
  ok = ok &&
       !diffexp2::local_algebra_detail::
            specialize_prepared_rational_multiplier_at_real_point(
                singular_denominator,
                RealEvaluationPoint::rational("-1/2")).has_value() &&
       !diffexp2::local_algebra_detail::
            specialize_prepared_rational_multiplier_at_real_point(
                exact_point_multiplier,
                RealEvaluationPoint::rational("0")).has_value() &&
       singular_denominator_diagnostic.has_value() &&
       singular_denominator_diagnostic->find("exactly zero") !=
           std::string::npos;

  const auto analytic_entry = [](
      std::vector<std::string> numerator) {
    PreparedRationalTaylorMultiplier<ComplexBall> entry;
    entry.epsilon_shift = 0;
    entry.center_pole_order = 0;
    entry.exact_identity = "point-matrix-entry";
    entry.kernels = {
        {ComplexBall::from_strings("777")},
        {ComplexBall::from_strings("0")},
        {ComplexBall::from_strings("0")}};
    std::vector<ComplexBall> parsed_numerator;
    for (const auto& coefficient : numerator)
      parsed_numerator.push_back(ComplexBall::from_strings(coefficient));
    entry.analytic_coefficients =
        std::vector<PreparedRationalAnalyticCoefficient<ComplexBall>>{
            {std::move(parsed_numerator),
             {ComplexBall::from_strings("1")}},
            {{ComplexBall::from_strings("0")},
             {ComplexBall::from_strings("1")}},
            {{ComplexBall::from_strings("0")},
             {ComplexBall::from_strings("1")}}};
    return entry;
  };
  PreparedSparseLocalMultiplierMatrix<ComplexBall> gauge_inverse;
  gauge_inverse.rows = gauge_inverse.columns = 2;
  gauge_inverse.entries.push_back({0, 0, analytic_entry({"1"})});
  gauge_inverse.entries.push_back({0, 1, analytic_entry({"0", "-1"})});
  gauge_inverse.entries.push_back({1, 1, analytic_entry({"1"})});
  PreparedSparseLocalMultiplierMatrix<ComplexBall> spectral_inverse;
  spectral_inverse.rows = spectral_inverse.columns = 2;
  spectral_inverse.entries.push_back({0, 1, analytic_entry({"1"})});
  spectral_inverse.entries.push_back({1, 0, analytic_entry({"1"})});
  PreparedSparseLocalMultiplierMatrix<ComplexBall> spectral =
      spectral_inverse;
  PreparedSparseLocalMultiplierMatrix<ComplexBall> gauge;
  gauge.rows = gauge.columns = 2;
  gauge.entries.push_back({0, 0, analytic_entry({"1"})});
  gauge.entries.push_back({0, 1, analytic_entry({"0", "1"})});
  gauge.entries.push_back({1, 1, analytic_entry({"1"})});
  FiniteLaurentVector<ComplexBall> physical_vector{
      EpsilonFrame<ComplexBall>(
          0, {ComplexBall::from_strings("2"), ComplexBall(0),
              ComplexBall(0)}),
      EpsilonFrame<ComplexBall>(
          0, {ComplexBall::from_strings("3"), ComplexBall(0),
              ComplexBall(0)})};
  const auto matrix_point = RealEvaluationPoint::rational("1/2");
  const auto reduced_vector =
      diffexp2::local_algebra_detail::
          apply_prepared_sparse_epsilon_matrix_at_real_point(
              gauge_inverse, physical_vector, matrix_point);
  const auto normal_vector = reduced_vector.has_value()
      ? diffexp2::local_algebra_detail::
            apply_prepared_sparse_epsilon_matrix_at_real_point(
                spectral_inverse, *reduced_vector, matrix_point)
      : std::nullopt;
  const auto assembled_vector = normal_vector.has_value()
      ? diffexp2::local_algebra_detail::
            apply_prepared_sparse_epsilon_matrix_at_real_point(
                spectral, *normal_vector, matrix_point)
      : std::nullopt;
  const auto restored_vector = assembled_vector.has_value()
      ? diffexp2::local_algebra_detail::
            apply_prepared_sparse_epsilon_matrix_at_real_point(
                gauge, *assembled_vector, matrix_point)
      : std::nullopt;
  ok = ok && normal_vector.has_value() &&
       normal_vector->front().complete_max() == 2 &&
       (normal_vector->at(0).coefficient(0) -
        ComplexBall::from_strings("3")).contains_zero() &&
       (normal_vector->at(1).coefficient(0) -
        ComplexBall::from_strings("1/2")).contains_zero() &&
       restored_vector.has_value() &&
       (restored_vector->at(0).coefficient(0) -
        ComplexBall::from_strings("2")).contains_zero() &&
       (restored_vector->at(1).coefficient(0) -
        ComplexBall::from_strings("3")).contains_zero();

  std::cout << (ok ? "PASS" : "FAIL")
            << ": native local rational/SCC and direct scalar-row algebra\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
