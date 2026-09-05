#include "diffexp/affine_operator.hpp"
#include <iostream>
using namespace diffexp;
using B = Jet::Ball;
namespace am = affine_matching;
namespace ao = affine_operator;
void require(bool v, const std::string &reason) {
  if (!v)
    throw std::runtime_error(reason);
}
void near(const B &a, const B &b) {
  require(acb_overlaps(a.raw(), b.raw()),
          "large affine analytic solution mismatch");
}
int main() {
  try {
    B::set_precision(256);
    ExactField field({"x", "eps"});
    Exact z(field, 0), eps(field, "eps");
    constexpr unsigned d = 109;
    auto connection = fuchsify::detail::zeros(d, d, z);
    AffineFrobeniusSeries::Options series_options;
    series_options.max_dimension = 256;
    auto series =
        AffineFrobeniusSeries::prepare(connection, 0, 1, 0, series_options);
    // A two-column epsilon collision inside109 sparse modes. The exact inverse
    // block is [[1,-1/eps],[0,1/eps]], with no determinant-subset expansion.
    auto transform = fuchsify::detail::identity(d, z);
    transform[0][1] = z.constant(1);
    transform[1][1] = eps;
    auto physical = series.project(transform);
    am::Options options;
    options.max_dimension = 256;
    options.max_determinant_dimension = 256;
    options.max_row_normalization_steps = 32;
    auto inverse =
        ao::prepare(series, physical, B::from_strings("0.125"), 2, options);
    require(inverse.success() && inverse.wronskian_verified &&
                inverse.normalization_steps == 1 &&
                inverse.determinant_valuation == 1,
            inverse.reason);
    near(inverse.inverse.coefficients[0][0][-inverse.inverse.low], B(1));
    near(inverse.inverse.coefficients[0][1][-1 - inverse.inverse.low], B(-1));
    near(inverse.inverse.coefficients[1][1][-1 - inverse.inverse.low], B(1));
    for (unsigned i = 2; i < d; ++i)
      near(inverse.inverse.coefficients[i][i][-inverse.inverse.low], B(1));
    am::Boundary boundary{0, 2, std::vector(d, std::vector<B>(3, B(0)))};
    for (unsigned i = 0; i < d; ++i)
      boundary.coefficients[i][0] = B(i + 1);
    auto matched = am::match(series, physical, B::from_strings("0.125"),
                             boundary, {-1, 0}, options);
    require(matched.success() && matched.used_row_normalization &&
                matched.wronskian_verified,
            matched.reason);
    near(matched.value.coefficients[0][-1 - matched.value.low], B(-2));
    near(matched.value.coefficients[0][-matched.value.low], B(1));
    near(matched.value.coefficients[1][-1 - matched.value.low], B(2));
    auto identity = ao::compose(inverse, series, physical,
                                B::from_strings("0.125"), 1, options);
    require(identity.success(), identity.reason);
    auto result = ao::apply_operator(identity, boundary, {0, 0}, options);
    require(result.success(), result.reason);
    for (unsigned i = 0; i < d; ++i)
      near(result.value.coefficients[i][-result.value.low], B(i + 1));
    // Distinct capacities: general row normalization admits109, the exponential
    // subset determinant stays bounded at12 even if callers raise its option.
    bool bounded = false;
    try {
      am::detail::Budget budget{options};
      am::detail::Matrix oversized(
          13, std::vector<am::detail::Series>(13, am::detail::Series(1)));
      (void)am::detail::determinant(oversized, 0, budget);
    } catch (const std::length_error &) {
      bounded = true;
    }
    require(bounded, "large determinant fallback did not reject");
    auto tiny = options;
    tiny.max_series_cells = 100;
    auto denied =
        ao::prepare(series, physical, B::from_strings("0.125"), 2, tiny);
    require(!denied.success() &&
                denied.reason.find("series-cell") != std::string::npos,
            "large coefficient allocation ignored finite budget");
    std::cout << "109-dimensional sparse epsilon resonance: inverse, matching, "
                 "projection and bounded allocation passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
