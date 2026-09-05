#pragma once
#include "diffexp/adjoint_transport.hpp"
#include <memory>

namespace diffexp::linear_boundary {
struct Options {
  // Intermediate master coordinates can be much wider than the independent
  // scalar leaf. In particular, Henn's 69/109-master maps must not inherit the
  // leaf's 64-component budget when composing through a single shared source.
  std::size_t max_rows = 5000, max_columns = 5000, max_coefficients = 20000000,
              max_operations = 50000000;
  int max_abs_order = 1000;
  std::size_t max_leaf_columns = 64;
};
// A transform and its shared, immutable source are kept together. Never
// replace this source with intermediate materialized componentwise balls.
struct Expression {
  LaurentRows transform;
  std::shared_ptr<const LaurentBoundary> leaf_source;
};
struct CompositionDemand : std::runtime_error {
  int required_outer_high, required_inner_high;
  CompositionDemand(int outer, int inner)
      : std::runtime_error("linear boundary composition requires additional "
                           "operator coefficients"),
        required_outer_high(outer), required_inner_high(inner) {}
};
struct OperatorDemand : std::runtime_error {
  int required_high;
  explicit OperatorDemand(int high)
      : std::runtime_error("linear boundary materialization requires "
                           "additional transform coefficients"),
        required_high(high) {}
};
namespace detail {
inline int order(long value, const Options &options) {
  if (value < -static_cast<long>(options.max_abs_order) ||
      value > options.max_abs_order)
    throw std::length_error("linear boundary epsilon-order budget exhausted");
  return static_cast<int>(value);
}
inline std::size_t product(std::initializer_list<std::size_t> factors,
                           std::size_t limit) {
  std::size_t out = 1;
  for (auto n : factors) {
    if (n && out > limit / n)
      throw std::length_error(
          "linear boundary finite resource budget exhausted");
    out *= n;
  }
  return out;
}
inline std::size_t validate(const LaurentRows &rows, const Options &options) {
  auto columns = rows.columns();
  if (!options.max_rows || !options.max_columns || !options.max_coefficients ||
      !options.max_operations || options.max_abs_order < 0 ||
      options.max_abs_order > 1000 ||
      rows.coefficients.size() > options.max_rows ||
      columns > options.max_columns)
    throw std::invalid_argument("linear boundary dimensions or budgets");
  order(rows.low, options);
  order(rows.high, options);
  product({rows.coefficients.size(), columns,
           static_cast<std::size_t>(rows.high - rows.low + 1)},
          options.max_coefficients);
  for (const auto &row : rows.coefficients)
    for (const auto &coefficients : row)
      for (const auto &value : coefficients)
        if (!value.is_finite())
          throw std::invalid_argument(
              "nonfinite linear boundary operator coefficient");
  return columns;
}
inline void source(const std::shared_ptr<const LaurentBoundary> &leaf,
                   const Options &options) {
  if (!leaf)
    throw std::invalid_argument(
        "linear boundary expression has no leaf source");
  if (!options.max_leaf_columns)
    throw std::invalid_argument("linear boundary leaf budget must be positive");
  auto high = leaf->high();
  order(leaf->low, options);
  order(high, options);
  if (leaf->values.size() > options.max_leaf_columns ||
      leaf->values.size() > options.max_columns)
    throw std::length_error("linear boundary leaf dimension budget exhausted");
  product({leaf->values.size(), static_cast<std::size_t>(high - leaf->low + 1)},
          options.max_coefficients);
  for (const auto &row : leaf->values)
    for (const auto &value : row)
      if (!value.is_finite())
        throw std::invalid_argument(
            "nonfinite linear boundary leaf coefficient");
}
} // namespace detail
inline Expression identity(std::shared_ptr<const LaurentBoundary> leaf,
                           int operator_high, const Options &options = {}) {
  detail::source(leaf, options);
  detail::order(operator_high, options);
  auto d = leaf->values.size();
  const int low = std::min(0, operator_high);
  detail::product({d, d, static_cast<std::size_t>(operator_high - low + 1)},
                  options.max_coefficients);
  LaurentRows rows{low, operator_high,
                   std::vector(d, std::vector(d, std::vector<Jet::Ball>(
                                                     operator_high - low + 1,
                                                     Jet::Ball(0))))};
  if (operator_high >= 0)
    for (unsigned i = 0; i < d; ++i)
      rows.coefficients[i][i][-low] = Jet::Ball(1);
  return {std::move(rows), std::move(leaf)};
}
// Orders below each declared low are structural zeros; orders above high
// are unknown. In particular, ball cancellations never tighten these bounds.
inline LaurentRows compose(const LaurentRows &outer, const LaurentRows &inner,
                           int high, const Options &options = {}) {
  using B = Jet::Ball;
  const auto middle = detail::validate(outer, options),
             columns = detail::validate(inner, options);
  if (middle != inner.coefficients.size())
    throw std::invalid_argument("linear boundary composition dimensions");
  detail::order(high, options);
  const int outer_high =
                detail::order(static_cast<long>(high) - inner.low, options),
            inner_high =
                detail::order(static_cast<long>(high) - outer.low, options);
  if (outer.high < outer_high || inner.high < inner_high)
    throw CompositionDemand(outer_high, inner_high);
  const int low =
      detail::order(std::min(static_cast<long>(high),
                             static_cast<long>(outer.low) + inner.low),
                    options);
  const auto rows = outer.coefficients.size();
  detail::product({rows, columns, static_cast<std::size_t>(high - low + 1)},
                  options.max_coefficients);
  LaurentRows out{
      low, high,
      std::vector(rows,
                  std::vector(columns, std::vector<B>(high - low + 1, B(0))))};
  std::size_t work = 0;
  for (unsigned i = 0; i < rows; ++i)
    for (unsigned k = 0; k < middle; ++k)
      for (int a = outer.low; a <= outer_high; ++a)
        if (!outer.coefficients[i][k][a - outer.low].is_zero())
          for (unsigned j = 0; j < columns; ++j)
            for (int b = inner.low;
                 b <= inner_high && static_cast<long>(a) + b <= high; ++b) {
              if (++work > options.max_operations)
                throw std::length_error(
                    "linear boundary composition operation budget exhausted");
              if (!inner.coefficients[k][j][b - inner.low].is_zero())
                out.coefficients[i][j][a + b - low] +=
                    outer.coefficients[i][k][a - outer.low] *
                    inner.coefficients[k][j][b - inner.low];
            }
  return out;
}
inline Expression compose(const LaurentRows &outer, const Expression &inner,
                          int high, const Options &options = {}) {
  detail::source(inner.leaf_source, options);
  if (detail::validate(inner.transform, options) !=
      inner.leaf_source->values.size())
    throw std::invalid_argument("linear boundary transform/source dimensions");
  return {compose(outer, inner.transform, high, options), inner.leaf_source};
}
inline LaurentBoundary materialize(const Expression &expression,
                                   int output_high,
                                   const Options &options = {}) {
  detail::source(expression.leaf_source, options);
  const auto columns = detail::validate(expression.transform, options);
  if (columns != expression.leaf_source->values.size())
    throw std::invalid_argument("linear boundary materialization dimensions");
  detail::order(output_high, options);
  const int transform_high = detail::order(static_cast<long>(output_high) -
                                               expression.leaf_source->low,
                                           options),
            source_high = detail::order(static_cast<long>(output_high) -
                                            expression.transform.low,
                                        options);
  if (expression.transform.high < transform_high)
    throw OperatorDemand(transform_high);
  if (expression.leaf_source->high() < source_high)
    throw BoundaryDemand(source_high, "linear boundary expression requires "
                                      "additional leaf-source coefficients");
  detail::product({expression.transform.coefficients.size(), columns,
                   static_cast<std::size_t>(expression.transform.high -
                                            expression.transform.low + 1),
                   static_cast<std::size_t>(expression.leaf_source->high() -
                                            expression.leaf_source->low + 1)},
                  options.max_operations);
  return apply_laurent_rows(expression.transform, *expression.leaf_source,
                            output_high);
}
} // namespace diffexp::linear_boundary
