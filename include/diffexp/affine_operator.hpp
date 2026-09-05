#pragma once
#include "diffexp/affine_matching.hpp"

namespace diffexp::affine_operator {
using B = kernel::ComplexBall;
using Status = affine_matching::Status;
using Options = affine_matching::Options;
using Window = affine_matching::Window;
using Boundary = affine_matching::Boundary;
using LaurentMatrix = AffineFrobeniusSeries::LaurentMatrix;
// Lower bounds are structural, not estimates from small ball coefficients.
// They can be conservative when projection cancels epsilon poles.
struct Operator {
  Status status = Status::Unsupported;
  std::string reason;
  LaurentMatrix matrix;
  std::vector<int> row_lower_bounds;
  int required_inverse_high = 0;
  bool omitted_tail_certified = false;
  bool success() const { return status == Status::Success; }
};
struct Prepared {
  Status status = Status::Unsupported;
  std::string reason;
  LaurentMatrix inverse;
  std::vector<int> row_lower_bounds;
  long determinant_valuation = 0;
  unsigned normalization_steps = 0, row_pole_loss = 0;
  bool wronskian_verified = false;
  bool used_numeric_convolution = false;
  unsigned numeric_convolution_fallbacks = 0;
  bool omitted_tail_certified = false;
  bool success() const { return status == Status::Success; }
};
namespace detail {
using BallMatrix = std::vector<std::vector<B>>;
inline BallMatrix zeros(unsigned rows, unsigned cols) {
  return BallMatrix(rows, std::vector<B>(cols, B(0)));
}
inline BallMatrix inverse(const BallMatrix &matrix) {
  unsigned d = matrix.size();
  acb_mat_t a, b;
  acb_mat_init(a, d, d);
  acb_mat_init(b, d, d);
  for (unsigned i = 0; i < d; ++i)
    for (unsigned j = 0; j < d; ++j)
      acb_set(acb_mat_entry(a, i, j), matrix[i][j].raw());
  bool ok = acb_mat_inv(b, a, B::precision());
  auto out = zeros(d, d);
  if (ok)
    for (unsigned i = 0; i < d; ++i)
      for (unsigned j = 0; j < d; ++j)
        acb_set(out[i][j].raw(), acb_mat_entry(b, i, j));
  acb_mat_clear(a);
  acb_mat_clear(b);
  if (!ok)
    throw affine_matching::detail::PrecisionNeeded(
        "prepared affine unit-leading matrix is not verifiably invertible");
  return out;
}
// Arithmetic enclosure guard only; this does not estimate omitted x terms.
// Use an upper radius norm and a lower midpoint norm so the test is conservative.
inline bool convolution_needs_exact(const BallMatrix &matrix) {
  mag_t radius, scale, one, quality;
  mag_init(radius); mag_init(scale); mag_init(one); mag_init(quality);
  mag_one(one);
  bool fallback = false;
  for (const auto &row : matrix) {
    for (const auto &value : row) {
      if (!value.is_finite()) { fallback = true; break; }
      B midpoint;
      acb_get_mid(midpoint.raw(), value.raw());
      acb_get_mag_lower(scale, midpoint.raw());
      mag_max(scale, scale, one);
      mag_add(radius, arb_radref(acb_realref(value.raw())),
                      arb_radref(acb_imagref(value.raw())));
      mag_div(quality, radius, scale);
      if (mag_cmp_2exp_si(quality, -B::precision() / 2) > 0) {
        fallback = true; break;
      }
    }
    if (fallback) break;
  }
  mag_clear(radius); mag_clear(scale); mag_clear(one); mag_clear(quality);
  return fallback;
}
inline void product_add(BallMatrix &out, const BallMatrix &a,
                        const BallMatrix &b, bool subtract = false,
                        affine_matching::detail::Budget *budget = nullptr) {
  for (unsigned i = 0; i < a.size(); ++i)
    for (unsigned k = 0; k < b.size(); ++k)
      if (!a[i][k].is_zero())
        for (unsigned j = 0; j < b[k].size(); ++j)
          if (!b[k][j].is_zero()) {
            if (budget)
              budget->spend();
            auto term = a[i][k] * b[k][j];
            out[i][j] = subtract ? out[i][j] - term : out[i][j] + term;
          }
}
} // namespace detail
// F = H diag(x^a eps^shift), K = R H. The same exact row normalization
// and K recurrence solve every inverse column together; uncertain boundary
// data never enters preparation. No determinant subset expansion is used.
inline Prepared prepare(const AffineFrobeniusSeries &series,
                        const AffineFrobeniusSeries::Expansion &physical,
                        const B &point, int inverse_high,
                        const Options &options = {}) {
  namespace am = affine_matching::detail;
  Prepared out;
  try {
    unsigned d = physical.rows;
    if (!d || d != physical.columns || d != series.dimension() ||
        d > options.max_dimension || options.max_dimension > 256 ||
        point.contains_zero() || !options.max_epsilon_depth ||
        options.max_row_normalization_steps >= options.max_epsilon_depth)
      throw std::invalid_argument(
          "invalid prepared affine inverse dimensions, point or budgets");
    auto guard = series.wronskian_valuation(physical, point);
    if (!physical.coherent_x_frontier || !guard)
      throw std::domain_error("prepared affine inverse requires coherent "
                              "absolute x frontiers and Wronskian guard");
    if (!guard->second.is_finite() || guard->second.contains_zero())
      throw am::PrecisionNeeded(
          "prepared affine Wronskian prefactor contains zero");
    am::Budget budget{options};
    am::ExactPointEvaluator evaluate_at(point);
    auto shifts = series.valuation_metadata(physical).minimum_by_column;
    for (auto shift : shifts)
      if (shift == AffineFrobeniusSeries::Valuations::zero_valuation)
        throw std::domain_error(
            "prepared affine inverse has an exact zero column");
    ExactField field({"x", "ell"});
    Exact zero(field, 0);
    auto gauge =
        am::adaptive_row_gauge(series, physical, shifts, point, zero, budget);
    unsigned loss = gauge.inverse_epsilon_powers.size() - 1;
    out.row_pole_loss = loss;
    out.normalization_steps = gauge.steps;
    out.determinant_valuation = gauge.sheared_rows;
    for (auto shift : shifts)
      out.determinant_valuation += shift;
    if (out.determinant_valuation != guard->first) {
      out.status = Status::NeedMoreXOrder;
      out.reason = "prepared retained frame disagrees with exact Wronskian "
                   "epsilon valuation";
      return out;
    }
    out.wronskian_verified = true;
    int low = inverse_high;
    long max_shift = *std::max_element(shifts.begin(), shifts.end());
    for (auto shift : shifts) {
      int lower = am::integer(-shift - static_cast<long>(loss));
      out.row_lower_bounds.push_back(lower);
      low = std::min(low, lower);
    }
    long depth =
        std::max(0L, static_cast<long>(inverse_high) + max_shift + loss);
    if (depth + loss >= options.max_epsilon_depth ||
        static_cast<long>(inverse_high) - low + 1 > options.max_epsilon_depth)
      throw std::length_error(
          "prepared affine inverse exceeds epsilon-depth budget");
    auto original = am::exact_coefficients(series, physical, shifts,
                                           depth + loss, zero, budget);
    for (unsigned q = 1; q <= loss; ++q) {
      auto residual = fuchsify::detail::zeros(d, d, zero);
      for (unsigned p = q; p <= loss; ++p) {
        auto product = am::exact_product(gauge.inverse_epsilon_powers[p],
                                         original[p - q], budget);
        for (unsigned i = 0; i < d; ++i)
          for (unsigned j = 0; j < d; ++j)
            residual[i][j] = residual[i][j] + product[i][j];
      }
      if (!fuchsify::detail::zero(residual))
        throw std::logic_error("prepared affine row gauge has uncertified "
                               "negative epsilon coefficients");
    }
    std::vector<detail::BallMatrix> k(depth + 1, detail::zeros(d, d)),
        r(loss + 1, detail::zeros(d, d));
    for (unsigned p = 0; p <= loss; ++p)
      for (unsigned i = 0; i < d; ++i)
        for (unsigned j = 0; j < d; ++j)
          r[p][i][j] = evaluate_at(gauge.inverse_epsilon_powers[p][i][j]);
    std::vector<detail::BallMatrix> h;
    if (options.numeric_operator_convolution && depth) {
      h.assign(depth + loss + 1, detail::zeros(d, d));
      for (unsigned n = 1; n < h.size(); ++n)
        for (unsigned i = 0; i < d; ++i)
          for (unsigned j = 0; j < d; ++j)
            h[n][i][j] = evaluate_at(original[n][i][j]);
      out.used_numeric_convolution = true;
    }
    for (unsigned n = 0; n <= depth; ++n) {
      bool exact_required = !n || !options.numeric_operator_convolution;
      if (!exact_required) {
        // R and H are independently evaluated enclosures. This is safe for
        // nonnegative coefficients, whose exact zeros are not used as pivots
        // or valuation evidence. It can widen balls when factors cancel.
        for (unsigned p = 0; p <= loss; ++p)
          detail::product_add(k[n], r[p], h[n + p], false, &budget);
        exact_required = detail::convolution_needs_exact(k[n]);
        if (exact_required) ++out.numeric_convolution_fallbacks;
      }
      if (exact_required) {
        // K0 is always formed exactly to retain stable pivot/rank decisions.
        auto exact = fuchsify::detail::zeros(d, d, zero);
        for (unsigned p = 0; p <= loss; ++p) {
          auto product = am::exact_product(gauge.inverse_epsilon_powers[p],
                                           original[n + p], budget);
          for (unsigned i = 0; i < d; ++i)
            for (unsigned j = 0; j < d; ++j)
              exact[i][j] = exact[i][j] + product[i][j];
        }
        for (unsigned i = 0; i < d; ++i)
          for (unsigned j = 0; j < d; ++j)
            k[n][i][j] = evaluate_at(exact[i][j]);
      }
    }
    auto inverse = detail::inverse(k[0]);
    std::vector<detail::BallMatrix> w(depth + 1, detail::zeros(d, d));
    for (unsigned n = 0; n <= depth; ++n) {
      auto rhs = detail::zeros(d, d);
      if (n <= loss)
        rhs = r[loss - n];
      for (unsigned lag = 1; lag <= n; ++lag)
        detail::product_add(rhs, k[lag], w[n - lag], true);
      detail::product_add(w[n], inverse, rhs);
    }
    out.inverse = {
        low, inverse_high,
        std::vector(
            d, std::vector(
                   d, std::vector<B>(static_cast<long>(inverse_high) - low + 1,
                                     B(0))))};
    B logx;
    acb_log(logx.raw(), point.raw(), B::precision());
    for (unsigned i = 0; i < d; ++i) {
      B factor;
      auto exponent =
          -affine_frobenius_detail::ball(series.exponents()[i].power) * logx;
      acb_exp(factor.raw(), exponent.raw(), B::precision());
      for (int e = std::max(low, out.row_lower_bounds[i]); e <= inverse_high;
           ++e) {
        long n = static_cast<long>(e) + shifts[i] + loss;
        for (unsigned j = 0; j < d; ++j)
          out.inverse.coefficients[i][j][e - low] = factor * w.at(n)[i][j];
      }
    }
    out.status = Status::Success;
  } catch (const affine_matching::detail::PrecisionNeeded &e) {
    out.status = Status::NeedMorePrecision;
    out.reason = e.what();
  } catch (const std::exception &e) {
    out.status = Status::Unsupported;
    out.reason = e.what();
  }
  return out;
}
inline Prepared prepare(const AffineFrobeniusSeries &series, const B &point,
                        int inverse_high, const Options &options = {}) {
  return prepare(series, series.terms(), point, inverse_high, options);
}
// Compose before boundary multiplication. Request operator_high >= desired
// output high - boundary.low; negative boundary poles need these extra terms.
inline Operator compose(const Prepared &prepared,
                        const AffineFrobeniusSeries &series,
                        const AffineFrobeniusSeries::Expansion &functional,
                        const B &point, int operator_high,
                        const Options &options = {}) {
  namespace am = affine_matching::detail;
  Operator out;
  try {
    if (!prepared.success()) {
      out.status = prepared.status;
      out.reason = prepared.reason;
      return out;
    }
    unsigned d = prepared.inverse.coefficients.size();
    if (!functional.rows || functional.columns != d)
      throw std::invalid_argument("affine operator functional dimensions");
    const long infinity = AffineFrobeniusSeries::Valuations::zero_valuation;
    std::vector<std::vector<long>> valuations(functional.rows,
                                              std::vector<long>(d, infinity));
    for (const auto &t : functional.terms)
      if (!t.coefficient.is_zero()) {
        AffineFrobeniusSeries::Expansion singleton{
            functional.rows, functional.columns, {t}};
        valuations[t.row][t.column] = std::min(
            valuations[t.row][t.column],
            series.valuation_metadata(singleton).minimum_by_column[t.column]);
      }
    long minimum = infinity, required = std::numeric_limits<long>::min();
    int low = operator_high;
    for (unsigned i = 0; i < functional.rows; ++i) {
      long row = infinity;
      for (unsigned j = 0; j < d; ++j)
        if (valuations[i][j] != infinity) {
          minimum = std::min(minimum, valuations[i][j]);
          row = std::min(row, valuations[i][j] + prepared.row_lower_bounds[j]);
          required = std::max(required, static_cast<long>(operator_high) -
                                            valuations[i][j]);
        }
      int bound = row == infinity ? operator_high + 1 : am::integer(row);
      out.row_lower_bounds.push_back(bound);
      low = std::min(low, bound);
    }
    out.required_inverse_high = required == std::numeric_limits<long>::min()
                                    ? prepared.inverse.high
                                    : am::integer(required);
    if (required > prepared.inverse.high) {
      out.status = Status::NeedMoreSource;
      out.reason =
          "functional composition requires a higher prepared inverse window";
      return out;
    }
    if (static_cast<long>(operator_high) - low + 1 > options.max_epsilon_depth)
      throw std::length_error("affine operator epsilon window exceeds budget");
    am::Budget{options}.series_cells(
        functional.rows, d, static_cast<long>(operator_high) - low + 1);
    out.matrix = {low, operator_high,
                  std::vector(functional.rows,
                              std::vector(d, std::vector<B>(static_cast<long>(
                                                                operator_high) -
                                                                low + 1,
                                                            B(0))))};
    if (minimum != infinity) {
      int ptop =
          am::integer(static_cast<long>(operator_high) -
                      *std::min_element(prepared.row_lower_bounds.begin(),
                                        prepared.row_lower_bounds.end()));
      auto p = series.evaluate_laurent(functional, point, am::integer(minimum),
                                       std::max(am::integer(minimum), ptop));
      for (unsigned i = 0; i < functional.rows; ++i)
        for (unsigned j = 0; j < d; ++j)
          for (unsigned k = 0; k < d; ++k)
            for (int a = p.low; a <= p.high; ++a)
              if (!p.coefficients[i][k][a - p.low].is_zero())
                for (int b = std::max(prepared.inverse.low,
                                      prepared.row_lower_bounds[k]);
                     b <= prepared.inverse.high &&
                     static_cast<long>(a) + b <= operator_high;
                     ++b)
                  if (static_cast<long>(a) + b >= low)
                    out.matrix.coefficients[i][j][a + b - low] =
                        out.matrix.coefficients[i][j][a + b - low] +
                        p.coefficients[i][k][a - p.low] *
                            prepared.inverse
                                .coefficients[k][j][b - prepared.inverse.low];
    }
    out.status = Status::Success;
  } catch (const std::exception &e) {
    out.status = Status::Unsupported;
    out.reason = e.what();
  }
  return out;
}
inline affine_matching::Result apply_operator(const Operator &op,
                                              const Boundary &boundary,
                                              Window requested,
                                              const Options &options = {}) {
  affine_matching::Result out;
  try {
    if (!op.success()) {
      out.status = op.status;
      out.reason = op.reason;
      return out;
    }
    unsigned rows = op.matrix.coefficients.size(),
             cols = op.matrix.coefficients.at(0).size();
    affine_matching::detail::validate(boundary, cols);
    if (requested.low > requested.high ||
        static_cast<long>(requested.high) - requested.low + 1 >
            options.max_epsilon_depth)
      throw std::invalid_argument("affine operator output window");
    long required = static_cast<long>(requested.high) -
                    *std::min_element(op.row_lower_bounds.begin(),
                                      op.row_lower_bounds.end());
    out.required_boundary_high = affine_matching::detail::integer(required);
    out.required_source_high = affine_matching::detail::integer(
        static_cast<long>(requested.high) - boundary.low);
    if (out.required_source_high > op.matrix.high) {
      out.status = Status::NeedMoreSource;
      out.reason = "operator needs higher coefficients for the declared "
                   "boundary lower order";
      return out;
    }
    if (required > boundary.high) {
      out.status = Status::NeedMoreBoundary;
      out.reason =
          "operator epsilon poles require additional boundary coefficients";
      return out;
    }
    int low = affine_matching::detail::integer(
        std::min(static_cast<long>(requested.low),
                 static_cast<long>(boundary.low) + op.matrix.low));
    out.value = {
        low, requested.high,
        std::vector(
            rows,
            std::vector<B>(static_cast<long>(requested.high) - low + 1, B(0)))};
    for (unsigned i = 0; i < rows; ++i)
      for (unsigned j = 0; j < cols; ++j)
        for (int a = op.matrix.low; a <= op.matrix.high; ++a)
          for (int b = boundary.low;
               b <= boundary.high && static_cast<long>(a) + b <= requested.high;
               ++b)
            if (static_cast<long>(a) + b >= low)
              out.value.coefficients[i][a + b - low] =
                  out.value.coefficients[i][a + b - low] +
                  op.matrix.coefficients[i][j][a - op.matrix.low] *
                      boundary.coefficients[j][b - boundary.low];
    out.status = Status::Success;
  } catch (const std::exception &e) {
    out.status = Status::Unsupported;
    out.reason = e.what();
  }
  return out;
}
} // namespace diffexp::affine_operator
