#pragma once
#include "diffexp/affine_matching.hpp"
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/ulong_extras.h>
#include <memory>
#include <mutex>

namespace diffexp::affine_operator {
using B = kernel::ComplexBall;
using Status = affine_matching::Status;
using Options = affine_matching::Options;
using Window = affine_matching::Window;
using Boundary = affine_matching::Boundary;
using LaurentMatrix = AffineFrobeniusSeries::LaurentMatrix;
// Lower bounds are structural, not estimates from small ball coefficients.
// They can be conservative when projection cancels epsilon poles.
struct PrincipalPartCache {
  std::mutex mutex;
  affine_matching::detail::ExactSeries normalized;
};
struct PrincipalPartData {
  affine_matching::detail::ExactSeries original, row_gauge;
  std::vector<long> shifts;
  B point;
  std::optional<Rational> exact_x;
  std::shared_ptr<PrincipalPartCache> cache = std::make_shared<PrincipalPartCache>();
};
namespace detail {
using BallMatrix = std::vector<std::vector<B>>;
struct Recurrence;
} // namespace detail
struct Operator {
  Status status = Status::Unsupported;
  std::string reason;
  LaurentMatrix matrix;
  std::vector<int> row_lower_bounds;
  int required_inverse_high = 0;
  bool principal_part_attempted = false, principal_part_certified = false;
  unsigned principal_part_zero_coefficients = 0;
  // Leading coefficients a modular evaluation proved to vanish exactly
  // (impose_exact_zeros).
  unsigned exact_zero_coefficients = 0;
  std::size_t principal_part_exact_operations = 0;
  std::string principal_part_reason;
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
  bool specialized_exact_point = false;
  bool omitted_tail_certified = false;
  std::shared_ptr<const PrincipalPartData> principal_part_data;
  // State of the epsilon recurrence behind `inverse`; see extend().
  std::shared_ptr<detail::Recurrence> recurrence;
  bool success() const { return status == Status::Success; }
};
namespace detail {
inline BallMatrix zeros(unsigned rows, unsigned cols) {
  return BallMatrix(rows, std::vector<B>(cols, B(0)));
}
inline BallMatrix inverse(const BallMatrix &matrix) {
  auto out = affine_matching::detail::verified_inverse(matrix);
  if (!out)
    throw affine_matching::detail::PrecisionNeeded(
        "prepared affine unit-leading matrix is not verifiably invertible");
  return std::move(*out);
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
// O=P F^-1 = eps^v J K^-1 R, where J=P diag(x^-a eps^(-shift-v)).
// Only a finite principal part is needed. K0 is the exact identity produced by
// row normalization; require that identity rather than guessing a numeric rank.
// Any failed/over-budget certificate leaves both balls and bounds untouched.
inline void certify_principal_part(Operator &out, const Prepared &prepared,
    const AffineFrobeniusSeries &series,
    const AffineFrobeniusSeries::Expansion &functional, const B &point,
    const std::vector<std::vector<long>> &valuations, const Options &options) {
  namespace am = affine_matching::detail;
  if (!options.certify_operator_principal_part || out.matrix.low >= 0) return;
  out.principal_part_attempted = true;
  am::Budget budget{options};
  budget.options.max_exact_operations = std::min(options.max_exact_operations,
                                                options.principal_part_max_operations);
  budget.options.max_exact_monomials = std::min<std::size_t>(options.max_exact_monomials, 2048);
  budget.deadline = std::chrono::steady_clock::now() +
      std::chrono::milliseconds(options.principal_part_max_milliseconds);
  try {
    budget.spend(0);
    if (!prepared.principal_part_data)
      throw std::domain_error("inverse was prepared without exact principal-part data");
    const auto &data = *prepared.principal_part_data;
    if (!acb_equal(point.raw(), data.point.raw()))
      throw std::domain_error("principal-part preparation point mismatch");
    budget.exact_x = data.exact_x;
    unsigned d = data.shifts.size(), rows = functional.rows;
    unsigned loss = data.row_gauge.size() - 1;
    long v = AffineFrobeniusSeries::Valuations::zero_valuation;
    for (unsigned i = 0; i < rows; ++i)
      for (unsigned j = 0; j < d; ++j)
        if (valuations[i][j] != AffineFrobeniusSeries::Valuations::zero_valuation)
          v = std::min(v, valuations[i][j] - data.shifts[j]);
    long top = static_cast<long>(loss) - v - 1;
    if (top < 0 || top >= options.principal_part_max_depth ||
        static_cast<unsigned long>(top) + loss >= data.original.size())
      throw std::length_error("principal-part depth exceeds available exact coefficients or budget");
    budget.series_cells(2 * rows + d, d, top + 1);
    auto shifts = data.shifts;
    for (auto &shift : shifts) shift += v;
    const auto zero = data.original[0][0][0].constant(0);
    auto j = am::exact_coefficients(series, functional, shifts, top, zero, budget);
    auto guard = [&](const Exact &value) {
      budget.spend(0);
      if (value.numerator_terms().size() + value.denominator_terms().size() >
          budget.options.max_exact_monomials)
        throw std::length_error("principal-part exact coefficient size budget exhausted");
    };
    auto product_add = [&](am::ExactMatrix &out, const am::ExactMatrix &a,
                           const am::ExactMatrix &b, bool subtract = false) {
      for (unsigned i = 0; i < a.size(); ++i)
        for (unsigned k = 0; k < b.size(); ++k)
          if (!a[i][k].is_zero())
            for (unsigned l = 0; l < b[k].size(); ++l)
              if (!b[k][l].is_zero()) {
                budget.spend(); guard(a[i][k]); guard(b[k][l]);
                auto term = a[i][k] * b[k][l];
                out[i][l] = subtract ? out[i][l] - term : out[i][l] + term;
                guard(out[i][l]);
              }
    };
    // Reuse exact K across functional rows. Publish each coefficient only
    // after it is complete; a failed extension never exposes a partial sum.
    // Keep the lock while using references, so concurrent cache growth cannot
    // invalidate them. Endpoint certification itself remains bounded/serial.
    std::unique_lock cache_lock(data.cache->mutex);
    auto &k = data.cache->normalized;
    budget.spend(0);
    while (k.size() <= static_cast<unsigned>(top)) {
      unsigned n = k.size();
      auto coefficient = fuchsify::detail::zeros(d, d, zero);
      for (unsigned p = 0; p <= loss; ++p)
        product_add(coefficient, data.row_gauge[p], data.original[n + p]);
      k.push_back(std::move(coefficient));
    }
    for (unsigned i = 0; i < d; ++i)
      for (unsigned l = 0; l < d; ++l)
        if (!(k[0][i][l] - zero.constant(i == l ? 1 : 0)).is_zero())
          throw std::domain_error("principal-part unit-leading normalization is not exact identity");
    // Sparse right solves Y K=J. Only the functional rows are propagated.
    auto y = std::move(j);
    for (unsigned n = 0; n < y.size(); ++n)
      for (unsigned lag = 1; lag <= n; ++lag)
        product_add(y[n], y[n-lag], k[lag], true);
    am::ExactSeries principal(top + 1, fuchsify::detail::zeros(rows, d, zero));
    for (unsigned n = 0; n < y.size(); ++n)
      for (unsigned p = n + loss > static_cast<unsigned>(top) ? n + loss - top : 0;
           p <= loss; ++p)
        product_add(principal[n + loss - p], y[n], data.row_gauge[p]);
    budget.spend(0);
    // Commit only after the entire exact proof succeeds. No omitted x term is
    // asserted to vanish, and no unknown positive epsilon coefficient is added.
    for (unsigned i = 0; i < rows; ++i) {
      int bound = out.row_lower_bounds[i];
      for (int e = std::max(out.matrix.low, bound); e < 0 && e <= out.matrix.high; ++e) {
        long n = static_cast<long>(e) - v + loss;
        bool row_zero = true;
        for (unsigned l = 0; l < d; ++l) {
          bool is_zero = n < 0 || (n < static_cast<long>(principal.size()) && principal[n][i][l].is_zero());
          if (is_zero) {
            out.matrix.coefficients[i][l][e-out.matrix.low] = B(0);
            ++out.principal_part_zero_coefficients;
          } else row_zero = false;
        }
        if (e == bound && row_zero) ++bound;
      }
      out.row_lower_bounds[i] = bound;
    }
    int certified_low = std::min(out.matrix.high,
        *std::min_element(out.row_lower_bounds.begin(), out.row_lower_bounds.end()));
    if (certified_low > out.matrix.low) {
      for (auto &row : out.matrix.coefficients)
        for (auto &entry : row)
          entry.erase(entry.begin(), entry.begin() + (certified_low - out.matrix.low));
      out.matrix.low = certified_low;
    }
    out.principal_part_certified = true;
  } catch (const std::exception &e) { out.principal_part_reason = e.what(); }
  out.principal_part_exact_operations = budget.work;
}
// The recurrence K W = R is triangular in epsilon: W_n needs K_0..K_n and
// nothing above. A prepared inverse therefore keeps what it has solved, and a
// later, higher epsilon order costs only its own coefficients.
//
// Exact arithmetic is confined to what decides structure: the row gauge,
// which fixes ranks and epsilon valuations from the orders 0..loss. With
// numeric_operator_convolution every K_n, the leading one included, is a ball
// product of the evaluated gauge and the x^r log(x)^k coefficients evaluated
// at the point; nothing is assembled as an exact rational function unless an
// enclosure turns out too wide. K_0 is inverted by a certified ball inverse,
// so no pivot or rank is ever read off an interval.
struct Recurrence {
  affine_matching::detail::Matrix symbolic;  // exact x^r log^k coefficients, by epsilon order
  unsigned symbolic_orders = 0;
  affine_matching::detail::ExactSeries exact, row_gauge;  // contiguous exact prefix; R
  affine_matching::detail::ExactMatrix leading;  // exact K_0 left by the row gauge
  std::vector<long> shifts;
  unsigned loss = 0;
  // R_0 is a certified ball inverse of H_0, not an evaluated exact gauge.
  bool numeric_gauge = false;
  B point, logarithm;
  slong precision = 0;
  Exact zero;
  std::vector<BallMatrix> r, h, k, w;
  BallMatrix leading_inverse;
  std::map<long, B> x_powers;
  std::vector<B> log_powers;
  Recurrence(const B &at, const Exact &exact_zero)
      : point(at), precision(B::precision()), zero(exact_zero) {
    acb_log(logarithm.raw(), point.raw(), precision);
  }
  // One pass over all series terms per call: ask for everything foreseeable.
  void require_symbolic(const AffineFrobeniusSeries::Expansion &physical,
                        unsigned top, affine_matching::detail::Budget &budget) {
    if (top < symbolic_orders) return;
    auto more = affine_matching::detail::coefficients(physical, shifts, top, budget, symbolic_orders);
    if (!symbolic_orders) symbolic = std::move(more);
    else
      for (unsigned i = 0; i < symbolic.size(); ++i)
        for (unsigned j = 0; j < symbolic[i].size(); ++j) {
          symbolic[i][j].resize(top + 1);
          for (unsigned n = symbolic_orders; n <= top; ++n)
            symbolic[i][j][n] = std::move(more[i][j][n]);
        }
    symbolic_orders = top + 1;
  }
  const affine_matching::detail::ExactMatrix &exact_order(const AffineFrobeniusSeries &series,
      const AffineFrobeniusSeries::Expansion &physical, unsigned n,
      unsigned foreseen, affine_matching::detail::Budget &budget) {
    require_symbolic(physical, std::max(n, foreseen), budget);
    while (exact.size() <= n) {
      const unsigned order = exact.size();
      auto next = fuchsify::detail::zeros(physical.rows, physical.columns, zero);
      for (unsigned i = 0; i < physical.rows; ++i)
        for (unsigned j = 0; j < physical.columns; ++j)
          next[i][j] = affine_matching::detail::as_exact(symbolic[i][j][order],
                                    series.exponents()[j].power, zero, budget);
      exact.push_back(std::move(next));
    }
    return exact[n];
  }
  // The inverse of the order-zero matrix modulo a prime, at the exact point
  // x and with the logarithm, which is transcendental there, replaced by a
  // fixed residue. An exact zero of the inverse, or of a product of exact
  // data with it, vanishes for every such choice; an entry that vanishes for
  // two primes is an exact zero except with a probability of the order of
  // its degree over the prime (2^-50 and less).
  struct Modular {
    mp_limb_t prime = 0, base = 0, log_residue = 0;
    std::vector<mp_limb_t> inverse;  // d*d, by rows
  };
  // Of the two primes, once both were usable; see prepare().
  std::vector<Modular> modular;
  // Residue of a sum of q x^r log^k in column j, whose base exponent is
  // factored out; nothing if a term is not of that form modulo the prime.
  std::optional<mp_limb_t> residue_modulo(const affine_matching::detail::Polynomial &terms, unsigned j,
      const AffineFrobeniusSeries &series, const Modular &at) const {
    nmod_t field;
    nmod_init(&field, at.prime);
    mp_limb_t sum = 0;
    for (const auto &[key, q] : terms) {
      const auto exponent = key.first - series.exponents()[j].power;
      const mp_limb_t denominator = fmpz_fdiv_ui(fmpq_denref(q.raw()), at.prime);
      if (!affine_frobenius_detail::is_integer(exponent) || !denominator) return std::nullopt;
      const mp_limb_t coefficient = nmod_div(fmpz_fdiv_ui(fmpq_numref(q.raw()), at.prime), denominator, field);
      const long power = affine_frobenius_detail::integer_value(exponent);
      mp_limb_t term = nmod_pow_ui(at.base, static_cast<ulong>(power < 0 ? -power : power), field);
      if (power < 0) term = nmod_inv(term, field);
      term = nmod_mul(term, nmod_pow_ui(at.log_residue, key.second, field), field);
      sum = nmod_add(sum, nmod_mul(term, coefficient, field), field);
    }
    return sum;
  }
  // Empty if the prime is unusable: it divides a denominator, or the matrix
  // is singular modulo it.
  std::optional<Modular> inverse_modulo(const AffineFrobeniusSeries &series, const Rational &x,
      mp_limb_t prime, mp_limb_t log_residue) const {
    nmod_t field;
    nmod_init(&field, prime);
    const mp_limb_t denominator = fmpz_fdiv_ui(fmpq_denref(x.raw()), prime);
    if (!denominator) return std::nullopt;
    Modular out{prime, nmod_div(fmpz_fdiv_ui(fmpq_numref(x.raw()), prime), denominator, field), log_residue % prime, {}};
    if (!out.base) return std::nullopt;
    const unsigned d = symbolic.size();
    nmod_mat_t matrix, inverse;
    nmod_mat_init(matrix, d, d, prime);
    nmod_mat_init(inverse, d, d, prime);
    bool usable = true;
    for (unsigned i = 0; i < d && usable; ++i)
      for (unsigned j = 0; j < d && usable; ++j) {
        const auto value = residue_modulo(symbolic[i][j][0], j, series, out);
        if (!value) usable = false;
        else nmod_mat_set_entry(matrix, i, j, *value);
      }
    usable = usable && nmod_mat_inv(inverse, matrix);
    if (usable) {
      out.inverse.resize(std::size_t(d) * d);
      for (unsigned i = 0; i < d; ++i)
        for (unsigned j = 0; j < d; ++j)
          out.inverse[std::size_t(i) * d + j] = nmod_mat_entry(inverse, i, j);
    }
    nmod_mat_clear(matrix);
    nmod_mat_clear(inverse);
    if (!usable) return std::nullopt;
    return out;
  }
  BallMatrix numeric_order(const AffineFrobeniusSeries &series, unsigned n,
                           affine_matching::detail::Budget &budget) {
    BallMatrix out = zeros(symbolic.size(), symbolic.at(0).size());
    for (unsigned i = 0; i < out.size(); ++i)
      for (unsigned j = 0; j < out[i].size(); ++j)
        for (const auto &[key, q] : symbolic[i][j][n]) {
          budget.spend();
          const auto exponent = key.first - series.exponents()[j].power;
          if (!affine_frobenius_detail::is_integer(exponent))
            throw std::domain_error("physical affine column does not have integer x "
                                    "shifts after its base exponent is factored");
          const long power = affine_frobenius_detail::integer_value(exponent);
          auto x = x_powers.find(power);
          if (x == x_powers.end()) {
            B value;
            acb_pow_si(value.raw(), point.raw(), power, precision);
            x = x_powers.emplace(power, std::move(value)).first;
          }
          while (log_powers.size() <= key.second)
            log_powers.push_back(log_powers.empty() ? B(1) : log_powers.back() * logarithm);
          out[i][j] = out[i][j] + affine_frobenius_detail::ball(q) * x->second *
                                      log_powers[key.second];
        }
    return out;
  }
  // The gauge must leave no negative epsilon order in R H. It guarantees
  // this by construction; a violation is a defect, reported loudly: exactly
  // in exact mode, and as an enclosure excluding zero otherwise.
  void regular_gauge(const AffineFrobeniusSeries &series,
                     const AffineFrobeniusSeries::Expansion &physical,
                     bool numeric, affine_matching::detail::Budget &budget) {
    const unsigned d = physical.rows;
    for (unsigned q = 1; q <= loss; ++q) {
      bool regular = true;
      if (numeric) {
        auto residual = zeros(d, d);
        for (unsigned p = q; p <= loss; ++p)
          product_add(residual, r[p], h[p - q], false, &budget);
        for (const auto &row : residual)
          for (const auto &value : row)
            regular = regular && value.is_finite() && value.contains_zero();
      } else {
        exact_order(series, physical, loss, loss, budget);
        auto residual = fuchsify::detail::zeros(d, d, zero);
        for (unsigned p = q; p <= loss; ++p) {
          auto product = affine_matching::detail::exact_product(row_gauge[p], exact[p - q], budget);
          for (unsigned i = 0; i < d; ++i)
            for (unsigned j = 0; j < d; ++j)
              residual[i][j] = residual[i][j] + product[i][j];
        }
        regular = fuchsify::detail::zero(residual);
      }
      if (!regular)
        throw std::logic_error("prepared affine row gauge has uncertified "
                               "negative epsilon coefficients");
    }
  }
  // Solve the orders k.size()..depth. Each completed order is appended whole,
  // so an exhausted budget leaves a shorter but consistent recurrence.
  void advance(const AffineFrobeniusSeries &series,
               const AffineFrobeniusSeries::Expansion &physical, unsigned depth,
               affine_matching::detail::Budget &budget, Prepared &out) {
    const unsigned d = physical.rows;
    const bool numeric = budget.options.numeric_operator_convolution;
    affine_matching::detail::ExactPointEvaluator evaluate_at(point);
    require_symbolic(physical, depth + loss, budget);
    if (numeric)
      while (h.size() <= depth + loss) {
        auto next = numeric_order(series, h.size(), budget);
        h.push_back(std::move(next));
      }
    if (k.empty()) regular_gauge(series, physical, numeric, budget);
    while (k.size() <= depth) {
      const unsigned n = k.size();
      auto value = zeros(d, d);
      if (!n) {
        // The gauge ends on K_0 itself, an exact permutation matrix.
        for (unsigned i = 0; i < d; ++i)
          for (unsigned j = 0; j < d; ++j)
            value[i][j] = evaluate_at(leading[i][j]);
        leading_inverse = inverse(value);
        k.push_back(std::move(value));
        continue;
      }
      bool exact_required = !numeric;
      if (!exact_required) {
        // R and H are independently evaluated enclosures. Their exact zeros
        // are not used as pivots or valuation evidence. The product can widen
        // balls when factors cancel; then the order is formed exactly.
        for (unsigned p = 0; p <= loss; ++p)
          product_add(value, r[p], h[n + p], false, &budget);
        exact_required = convolution_needs_exact(value);
        if (exact_required) ++out.numeric_convolution_fallbacks;
      }
      if (exact_required) {
        if (numeric_gauge)
          throw affine_matching::detail::PrecisionNeeded(
              "numeric row gauge: an epsilon order of the prepared inverse is too wide");
        exact_order(series, physical, n + loss, depth + loss, budget);
        auto sum = fuchsify::detail::zeros(d, d, zero);
        for (unsigned p = 0; p <= loss; ++p) {
          auto product = affine_matching::detail::exact_product(row_gauge[p], exact[n + p], budget);
          for (unsigned i = 0; i < d; ++i)
            for (unsigned j = 0; j < d; ++j)
              sum[i][j] = sum[i][j] + product[i][j];
        }
        for (unsigned i = 0; i < d; ++i)
          for (unsigned j = 0; j < d; ++j)
            value[i][j] = evaluate_at(sum[i][j]);
      }
      k.push_back(std::move(value));
    }
    while (w.size() <= depth) {
      const unsigned n = w.size();
      auto rhs = n <= loss ? r[loss - n] : zeros(d, d);
      for (unsigned lag = 1; lag <= n; ++lag)
        product_add(rhs, k[lag], w[n - lag], true);
      auto value = zeros(d, d);
      product_add(value, leading_inverse, rhs);
      w.push_back(std::move(value));
    }
    out.used_numeric_convolution = numeric && depth;
  }
};
// Publish the window through inverse_high from the solved recurrence.
inline void publish(Prepared &out, const AffineFrobeniusSeries &series,
                    const AffineFrobeniusSeries::Expansion &physical,
                    int inverse_high, affine_matching::detail::Budget &budget) {
  auto &state = *out.recurrence;
  const auto &options = budget.options;
  const unsigned d = physical.rows, loss = state.loss;
  int low = inverse_high;
  const long max_shift = *std::max_element(state.shifts.begin(), state.shifts.end());
  for (auto lower : out.row_lower_bounds) low = std::min(low, lower);
  const long depth =
      std::max(0L, static_cast<long>(inverse_high) + max_shift + loss);
  if (depth + loss >= options.max_epsilon_depth ||
      static_cast<long>(inverse_high) - low + 1 > options.max_epsilon_depth)
    throw std::length_error(
        "prepared affine inverse exceeds epsilon-depth budget");
  state.advance(series, physical, depth, budget, out);
  out.inverse = {
      low, inverse_high,
      std::vector(
          d, std::vector(
                 d, std::vector<B>(static_cast<long>(inverse_high) - low + 1,
                                   B(0))))};
  for (unsigned i = 0; i < d; ++i) {
    B factor;
    auto exponent =
        -affine_frobenius_detail::ball(series.exponents()[i].power) * state.logarithm;
    acb_exp(factor.raw(), exponent.raw(), B::precision());
    for (int e = std::max(low, out.row_lower_bounds[i]); e <= inverse_high;
         ++e) {
      long n = static_cast<long>(e) + state.shifts[i] + loss;
      for (unsigned j = 0; j < d; ++j)
        out.inverse.coefficients[i][j][e - low] = factor * state.w.at(n)[i][j];
    }
  }
  if (options.certify_operator_principal_part && !state.numeric_gauge) {
    // The certificate reads exact coefficients through depth+loss.
    state.exact_order(series, physical, depth + loss, depth + loss, budget);
    out.principal_part_data = std::make_shared<PrincipalPartData>(
        PrincipalPartData{state.exact, state.row_gauge, state.shifts,
                          state.point, budget.exact_x});
  }
}
// A functional composed with an inverse that a ball inversion produced has
// balls of rounding size where the composition vanishes exactly, and rows
// built on it carry every such column through the transport. The lowest
// epsilon order of row i of the operator is
//   sum_k P_ik x^(-a_k) (H_0^-1)_kj
// with the coefficient P_ik of the functional at that order relative to the
// shift of column k: a product of exact data, evaluated modulo two primes
// like the zeros of the inverse itself, with the inverse modulo each prime
// that the preparation kept (Recurrence::modular). An entry
// whose ball contains zero and which vanishes for both is the exact zero. No
// size is compared with anything.
inline void impose_exact_zeros(Operator &out, const Prepared &prepared,
    const AffineFrobeniusSeries &series,
    const AffineFrobeniusSeries::Expansion &functional, const Options &options) {
  namespace am = affine_matching::detail;
  if (!prepared.recurrence || !prepared.recurrence->numeric_gauge || prepared.recurrence->loss ||
      prepared.recurrence->modular.size() != 2)
    return;
  const auto &state = *prepared.recurrence;
  try {
    am::Budget budget{options};
    const unsigned d = state.symbolic.size(), rows = functional.rows;
    if (functional.columns != d || out.row_lower_bounds.size() != rows) return;
    // One row at a time: its leading order differs from row to row.
    std::vector<AffineFrobeniusSeries::Expansion> split;
    if (rows > 1) {
      split.resize(rows);
      for (auto &row : split) { row.rows = 1; row.columns = d; }
      for (const auto &t : functional.terms)
        if (!t.coefficient.is_zero()) {
          split[t.row].terms.push_back(t);
          split[t.row].terms.back().row = 0;
        }
    }
    for (unsigned i = 0; i < rows; ++i) {
      const auto &row = rows > 1 ? split[i] : functional;
      const long bound = out.row_lower_bounds[i];
      if (row.terms.empty() || bound < out.matrix.low || bound > out.matrix.high) continue;
      std::vector<long> shifts(d);
      for (unsigned k = 0; k < d; ++k) shifts[k] = bound + state.shifts[k];
      const auto leading = am::coefficients(row, shifts, 0, budget);
      // The row times the inverse, modulo both primes.
      std::vector<char> vanishes(d, 1);
      bool usable = true;
      for (const auto &at : state.modular) {
        nmod_t field;
        nmod_init(&field, at.prime);
        std::vector<mp_limb_t> left(d);
        for (unsigned k = 0; k < d && usable; ++k) {
          const auto value = state.residue_modulo(leading[0][k].at(0), k, series, at);
          if (!value) usable = false;
          else left[k] = *value;
        }
        if (!usable) break;
        for (unsigned j = 0; j < d; ++j) {
          mp_limb_t sum = 0;
          for (unsigned k = 0; k < d; ++k)
            if (left[k]) sum = nmod_add(sum, nmod_mul(left[k], at.inverse[std::size_t(k) * d + j], field), field);
          if (sum) vanishes[j] = 0;
        }
      }
      if (!usable) continue;
      for (unsigned j = 0; j < d; ++j) {
        auto &value = out.matrix.coefficients[i][j][bound - out.matrix.low];
        if (vanishes[j] && !value.is_zero() && value.contains_zero()) {
          value = B(0);
          ++out.exact_zero_coefficients;
        }
      }
    }
  } catch (const std::exception &) {
    // An exhausted budget or an unusable term leaves the balls as they are.
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
    am::Budget budget{options};budget.specialize(point);
    out.specialized_exact_point=budget.exact_x.has_value();
    ExactField field({"x", "ell"});
    auto state = std::make_shared<detail::Recurrence>(point, Exact(field, 0));
    state->shifts = series.valuation_metadata(physical).minimum_by_column;
    const auto &shifts = state->shifts;
    for (auto shift : shifts)
      if (shift == AffineFrobeniusSeries::Valuations::zero_valuation)
        throw std::domain_error(
            "prepared affine inverse has an exact zero column");
    // Epsilon orders are converted to the exact field only when an exact zero
    // leading row requires another shear; the first pass over the series
    // already covers the window that a gauge without loss would need.
    const long max_shift = *std::max_element(shifts.begin(), shifts.end());
    const unsigned foreseen = static_cast<unsigned>(std::min<long>(
        options.max_epsilon_depth, std::max(0L, static_cast<long>(inverse_high) + max_shift)));
    std::optional<detail::BallMatrix> numeric_inverse;
    if (options.numeric_leading_gauge && options.numeric_operator_convolution) {
      state->require_symbolic(physical, foreseen, budget);
      numeric_inverse = am::verified_inverse(state->numeric_order(series, 0, budget));
      if (numeric_inverse && budget.exact_x) {
        // The exact gauge returns the exact zeros of the inverse, and the
        // rows built on it carry through the transport only the columns that
        // are theirs. A ball inverse returns small balls instead. An entry
        // is an exact zero again if its ball contains zero and it vanishes
        // in the inverse modulo two primes.
        auto first = state->inverse_modulo(series, *budget.exact_x, n_nextprime(UWORD(1) << 61, 1), UWORD(0x9E3779B97F4A7C15));
        auto second = state->inverse_modulo(series, *budget.exact_x, n_nextprime(UWORD(1) << 62, 1), UWORD(0xC2B2AE3D27D4EB4F));
        if (first && second) {
          for (unsigned i = 0; i < d; ++i)
            for (unsigned j = 0; j < d; ++j)
              if (!first->inverse[std::size_t(i) * d + j] && !second->inverse[std::size_t(i) * d + j] &&
                  (*numeric_inverse)[i][j].contains_zero())
                (*numeric_inverse)[i][j] = B(0);
          // Kept for the operators composed with this inverse.
          state->modular.push_back(std::move(*first));
          state->modular.push_back(std::move(*second));
        }
      }
    }
    auto gauge = numeric_inverse
        ? am::RowGauge{{fuchsify::detail::identity(d, state->zero)}, 0, 0, {},
                       fuchsify::detail::identity(d, state->zero)}
        : am::row_gauge_coefficients(
              d, state->zero, point, budget,
              [&](unsigned n) -> const am::ExactMatrix & {
                return state->exact_order(series, physical, n, foreseen, budget);
              });
    state->numeric_gauge = numeric_inverse.has_value();
    const unsigned loss = gauge.inverse_epsilon_powers.size() - 1;
    state->loss = loss;
    state->row_gauge = std::move(gauge.inverse_epsilon_powers);
    state->leading = std::move(gauge.leading);
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
    for (auto shift : shifts)
      out.row_lower_bounds.push_back(am::integer(-shift - static_cast<long>(loss)));
    am::ExactPointEvaluator evaluate_at(point);
    state->r.assign(loss + 1, detail::zeros(d, d));
    for (unsigned p = 0; p <= loss; ++p)
      for (unsigned i = 0; i < d; ++i)
        for (unsigned j = 0; j < d; ++j)
          state->r[p][i][j] = evaluate_at(state->row_gauge[p][i][j]);
    if (numeric_inverse) state->r[0] = std::move(*numeric_inverse);
    out.recurrence = std::move(state);
    detail::publish(out, series, physical, inverse_high, budget);
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
// Raise a prepared inverse to a higher epsilon order in place. Orders already
// solved are final, so this costs the new coefficients only. A lower or equal
// request leaves the (wider) window untouched.
inline void extend(Prepared &prepared, const AffineFrobeniusSeries &series,
                   const AffineFrobeniusSeries::Expansion &physical,
                   int inverse_high, const Options &options = {}) {
  if (!prepared.success() || inverse_high <= prepared.inverse.high) return;
  if (!prepared.recurrence || prepared.recurrence->precision != B::precision())
    throw std::logic_error("affine inverse extension requires the recurrence "
                           "state at its original working precision");
  try {
    affine_matching::detail::Budget budget{options};
    budget.specialize(prepared.recurrence->point);
    detail::publish(prepared, series, physical, inverse_high, budget);
  } catch (const affine_matching::detail::PrecisionNeeded &e) {
    prepared.status = Status::NeedMorePrecision;
    prepared.reason = e.what();
  } catch (const std::exception &e) {
    prepared.status = Status::Unsupported;
    prepared.reason = e.what();
  }
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
    for (const auto &t : functional.terms) {
      if (t.row >= functional.rows || t.column >= d)
        throw std::invalid_argument(
            "affine operator functional coordinate outside shape");
      if (!t.coefficient.is_zero())
        valuations[t.row][t.column] =
            std::min(valuations[t.row][t.column],
                     series.epsilon_valuation(t.coefficient));
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
    if (minimum != infinity && options.operator_exact_zeros)
      detail::impose_exact_zeros(out, prepared, series, functional, options);
    detail::certify_principal_part(out, prepared, series, functional, point,
                                   valuations, options);
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
