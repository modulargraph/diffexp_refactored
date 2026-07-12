#include "diffexp2/json_codec.hpp"

#include "diffexp2/checkpoint.hpp"
#include "diffexp2/line_integration.hpp"
#include "diffexp2/local_algebra.hpp"
#include "diffexp2/local_solution.hpp"
#include "diffexp2/matching.hpp"
#include "diffexp2/path_planner.hpp"
#include "diffexp2/recurrence.hpp"
#include "diffexp2/singular_indicial.hpp"
#include "diffexp2/tail_majorant.hpp"

#include <boost/json.hpp>

#include <array>
#include <chrono>
#include <atomic>
#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>

namespace diffexp2 {
namespace json = boost::json;
namespace {

const json::object& as_object(const json::value& value, const char* label) {
  if (!value.is_object()) throw std::invalid_argument(std::string(label) + " must be an object");
  return value.as_object();
}

const json::array& as_array(const json::value& value, const char* label) {
  if (!value.is_array()) throw std::invalid_argument(std::string(label) + " must be an array");
  return value.as_array();
}

std::int64_t as_i64(const json::value& value, const char* label) {
  if (value.is_int64()) return value.as_int64();
  if (value.is_uint64() && value.as_uint64() <=
      static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
    return static_cast<std::int64_t>(value.as_uint64());
  throw std::invalid_argument(std::string(label) + " must be an integer");
}

std::uint32_t as_u32(const json::value& value, const char* label) {
  const auto number = as_i64(value, label);
  if (number < 0 || number > std::numeric_limits<std::uint32_t>::max())
    throw std::invalid_argument(std::string(label) + " is outside uint32 range");
  return static_cast<std::uint32_t>(number);
}

std::uint64_t as_u64(const json::value& value, const char* label) {
  if (value.is_uint64()) return value.as_uint64();
  if (value.is_int64() && value.as_int64() >= 0)
    return static_cast<std::uint64_t>(value.as_int64());
  throw std::invalid_argument(std::string(label) +
                              " must be a nonnegative integer");
}

double as_double(const json::value& value, const char* label) {
  if (value.is_double()) return value.as_double();
  if (value.is_int64()) return static_cast<double>(value.as_int64());
  if (value.is_uint64()) return static_cast<double>(value.as_uint64());
  throw std::invalid_argument(std::string(label) + " must be numeric");
}

std::int32_t as_i32(const json::value& value, const char* label) {
  const auto number = as_i64(value, label);
  if (number < std::numeric_limits<std::int32_t>::min() ||
      number > std::numeric_limits<std::int32_t>::max())
    throw std::invalid_argument(std::string(label) + " is outside int32 range");
  return static_cast<std::int32_t>(number);
}

std::int32_t parse_validity(const json::value& value) {
  return value.is_null() ? kCompleteInfinity : as_i32(value, "validity");
}

template <typename Scalar>
Scalar parse_scalar(const json::value& value);

template <>
Rational parse_scalar<Rational>(const json::value& value) {
  if (value.is_string()) return Rational(std::string(value.as_string()));
  if (value.is_int64()) return Rational(static_cast<long>(value.as_int64()));
  const auto& pair = as_array(value, "rational scalar");
  if (pair.size() != 2 || !pair[0].is_string() || !pair[1].is_string() ||
      std::string(pair[1].as_string()) != "0")
    throw std::invalid_argument("exact rational scalar must have zero imaginary part");
  return Rational(std::string(pair[0].as_string()));
}

template <>
ComplexBall parse_scalar<ComplexBall>(const json::value& value) {
  if (value.is_string()) return ComplexBall::from_strings(std::string(value.as_string()));
  if (value.is_int64()) return ComplexBall(static_cast<long>(value.as_int64()));
  const auto& pair = as_array(value, "Acb scalar");
  if (pair.size() != 2 || !pair[0].is_string() || !pair[1].is_string())
    throw std::invalid_argument("Acb scalar must be [real-string,imag-string]");
  return ComplexBall::from_strings(std::string(pair[0].as_string()),
                                   std::string(pair[1].as_string()));
}

template <>
SymbolicRational parse_scalar<SymbolicRational>(const json::value& value) {
  if (value.is_string())
    return SymbolicRational(std::string(value.as_string()));
  if (value.is_int64())
    return SymbolicRational(static_cast<long>(value.as_int64()));
  throw std::invalid_argument("symbolic scalar must be a rational-function string");
}

template <typename Scalar>
MatrixShift<Scalar> parse_matrix_shift(const json::value& value,
                                       std::uint32_t dimension) {
  const auto& object = as_object(value, "matrix shift");
  MatrixShift<Scalar> out;
  out.shift = as_i32(object.at("s"), "matrix shift exponent");
  for (const auto& raw_entry : as_array(object.at("e"), "matrix entries")) {
    const auto& entry = as_array(raw_entry, "matrix entry");
    if (entry.size() != 3) throw std::invalid_argument("matrix entry must be [row,col,value]");
    MatrixEntry<Scalar> item;
    item.row = as_u32(entry[0], "matrix row");
    item.col = as_u32(entry[1], "matrix column");
    if (item.row >= dimension || item.col >= dimension)
      throw std::invalid_argument("matrix entry outside recurrence dimension");
    item.value = parse_scalar<Scalar>(entry[2]);
    if (!ScalarTraits<Scalar>::is_zero(item.value)) out.entries.push_back(std::move(item));
  }
  return out;
}

template <typename Scalar>
RecurrenceProblem<Scalar> parse_problem(const json::object& root) {
  RecurrenceProblem<Scalar> p;
  p.dimension = as_u32(root.at("d"), "dimension");
  p.nmax = as_u32(root.at("nmax"), "nmax");
  p.log_max = as_u32(root.at("p"), "log maximum");
  p.frame_base = as_i32(root.at("fb"), "frame base");
  p.frame_width = as_u32(root.at("w"), "frame width");
  p.has_initial = root.if_contains("has_initial") == nullptr ||
                  root.at("has_initial").as_bool();
  p.adaptive_lower_frame_probe = root.if_contains("adaptive_probe") != nullptr &&
                                  root.at("adaptive_probe").as_bool();
  p.a_target = parse_scalar<Scalar>(root.at("a_target"));
  p.b_target = parse_scalar<Scalar>(root.at("b_target"));
  p.a_shift_min = as_i32(root.at("a_shift_min"), "a shift minimum");
  for (const auto& value : as_array(root.at("a_shifts"), "a shifts"))
    p.a_shifts.push_back(parse_scalar<Scalar>(value));

  for (const auto& raw_lag : as_array(root.at("d_lags"), "d lags")) {
    std::vector<ScalarShift<Scalar>> lag;
    for (const auto& raw_shift : as_array(raw_lag, "d lag")) {
      const auto& shift = as_object(raw_shift, "d scalar shift");
      lag.push_back({as_i32(shift.at("s"), "d shift"),
                     parse_scalar<Scalar>(shift.at("v"))});
    }
    p.d_lags.push_back(std::move(lag));
  }

  for (const auto& raw_den : as_array(root.at("denominators"), "denominators")) {
    std::vector<Scalar> den;
    for (const auto& value : as_array(raw_den, "denominator"))
      den.push_back(parse_scalar<Scalar>(value));
    p.rational_denominators.push_back(std::move(den));
  }

  for (const auto& raw_lag : as_array(root.at("nhat_lags"), "Nhat lags")) {
    const auto& lag_object = as_object(raw_lag, "Nhat lag");
    PreparedLag<Scalar> lag;
    for (const auto& raw_matrix : as_array(lag_object.at("poly"), "polynomial matrices"))
      lag.polynomial.push_back(parse_matrix_shift<Scalar>(raw_matrix, p.dimension));
    for (const auto& raw_group : as_array(lag_object.at("rat"), "rational groups")) {
      const auto& group_object = as_object(raw_group, "rational group");
      RationalGroup<Scalar> group;
      group.denominator_index = as_u32(group_object.at("q"), "denominator index");
      for (const auto& raw_matrix : as_array(group_object.at("num"), "rational numerator"))
        group.numerator.push_back(parse_matrix_shift<Scalar>(raw_matrix, p.dimension));
      lag.rational.push_back(std::move(group));
    }
    for (const auto& value : as_array(lag_object.at("val"), "Nhat valuations"))
      lag.valuations.push_back(parse_validity(value));
    p.nhat_lags.push_back(std::move(lag));
  }

  if (!root.at("d0_inverse").is_null())
    p.d0_inverse_scalar = parse_scalar<Scalar>(root.at("d0_inverse"));

  for (const auto& raw_block : as_array(root.at("blocks"), "Jordan blocks")) {
    JordanBlock block;
    for (const auto& col : as_array(raw_block, "Jordan block"))
      block.columns.push_back(as_u32(col, "Jordan column"));
    p.blocks.push_back(std::move(block));
  }

  for (const auto& raw_row : as_array(root.at("schedule"), "step schedule")) {
    std::vector<BlockStep<Scalar>> row;
    for (const auto& raw_step : as_array(raw_row, "schedule row")) {
      const auto& step_object = as_object(raw_step, "schedule step");
      const auto kind = std::string(step_object.at("case").as_string());
      StepCase step_case;
      if (kind == "T") step_case = StepCase::Taylor;
      else if (kind == "P") step_case = StepCase::Pseudo;
      else if (kind == "R") step_case = StepCase::Resonant;
      else throw std::invalid_argument("unknown recurrence step case: " + kind);
      row.push_back({step_case, parse_scalar<Scalar>(step_object.at("da")),
                     parse_scalar<Scalar>(step_object.at("db"))});
    }
    p.schedule.push_back(std::move(row));
  }

  for (const auto& value : as_array(root.at("initial"), "initial tensor"))
    p.initial.push_back(parse_scalar<Scalar>(value));
  for (const auto& value : as_array(root.at("initial_validity"), "initial validity"))
    p.initial_validity.push_back(parse_validity(value));

  if (const auto* raw_source = root.if_contains("source"); raw_source && !raw_source->is_null()) {
    const auto& source_object = as_object(*raw_source, "source");
    SourceData<Scalar> source;
    for (const auto& value : as_array(source_object.at("frames"), "source frames"))
      source.frames.push_back(parse_scalar<Scalar>(value));
    for (const auto& value : as_array(source_object.at("validity"), "source validity"))
      source.validity.push_back(parse_validity(value));
    for (const auto& value : as_array(source_object.at("present"), "source presence"))
      source.present.push_back(value.as_bool() ? 1 : 0);
    p.source = std::move(source);
  }
  if (const auto* raw_assembly = root.if_contains("assembly");
      raw_assembly && !raw_assembly->is_null()) {
    const auto& assembly_object = as_object(*raw_assembly, "assembly matrix");
    PreparedMatrix<Scalar> matrix;
    matrix.identity = assembly_object.if_contains("identity") != nullptr &&
                      assembly_object.at("identity").as_bool();
    for (const auto& raw_matrix : as_array(
             assembly_object.at("poly"), "assembly polynomial matrices"))
      matrix.polynomial.push_back(
          parse_matrix_shift<Scalar>(raw_matrix, p.dimension));
    for (const auto& raw_group : as_array(
             assembly_object.at("rat"), "assembly rational groups")) {
      const auto& group_object = as_object(raw_group, "assembly rational group");
      RationalGroup<Scalar> group;
      group.denominator_index = as_u32(group_object.at("q"), "assembly denominator index");
      for (const auto& raw_matrix : as_array(
               group_object.at("num"), "assembly rational numerator"))
        group.numerator.push_back(
            parse_matrix_shift<Scalar>(raw_matrix, p.dimension));
      matrix.rational.push_back(std::move(group));
    }
    for (const auto& value : as_array(
             assembly_object.at("val"), "assembly valuations"))
      matrix.valuations.push_back(parse_validity(value));
    p.assembly_matrix = std::move(matrix);
    p.chop_digits = as_i32(root.at("chop_digits"), "chop digits");
    p.return_u = root.if_contains("return_u") != nullptr &&
                 root.at("return_u").as_bool();
  }
  return p;
}

/* A prepared chart/window operator contains only the tensors that are
   invariant across homogeneous columns and inhomogeneous sectors.  The
   frame is part of the identity: every epsilon shift and valuation was
   prepared for exactly this [frame_base, frame_width] rectangle. */
template <typename Scalar>
PreparedRecurrenceOperator<Scalar> parse_prepared_operator(
    const json::object& root) {
  PreparedRecurrenceOperator<Scalar> prepared;
  prepared.dimension = as_u32(root.at("d"), "dimension");
  prepared.frame_base = as_i32(root.at("fb"), "frame base");
  prepared.frame_width = as_u32(root.at("w"), "frame width");
  for (const auto& raw_lag : as_array(root.at("d_lags"), "d lags")) {
    std::vector<ScalarShift<Scalar>> lag;
    for (const auto& raw_shift : as_array(raw_lag, "d lag")) {
      const auto& shift = as_object(raw_shift, "d scalar shift");
      lag.push_back({as_i32(shift.at("s"), "d shift"),
                     parse_scalar<Scalar>(shift.at("v"))});
    }
    prepared.d_lags.push_back(std::move(lag));
  }
  for (const auto& raw_den : as_array(
           root.at("denominators"), "denominators")) {
    std::vector<Scalar> denominator;
    for (const auto& value : as_array(raw_den, "denominator"))
      denominator.push_back(parse_scalar<Scalar>(value));
    prepared.rational_denominators.push_back(std::move(denominator));
  }
  for (const auto& raw_lag : as_array(
           root.at("nhat_lags"), "Nhat lags")) {
    const auto& lag_object = as_object(raw_lag, "Nhat lag");
    PreparedLag<Scalar> lag;
    for (const auto& raw_matrix : as_array(
             lag_object.at("poly"), "polynomial matrices"))
      lag.polynomial.push_back(
          parse_matrix_shift<Scalar>(raw_matrix, prepared.dimension));
    for (const auto& raw_group : as_array(
             lag_object.at("rat"), "rational groups")) {
      const auto& group_object = as_object(raw_group, "rational group");
      RationalGroup<Scalar> group;
      group.denominator_index = as_u32(
          group_object.at("q"), "denominator index");
      for (const auto& raw_matrix : as_array(
               group_object.at("num"), "rational numerator"))
        group.numerator.push_back(
            parse_matrix_shift<Scalar>(raw_matrix, prepared.dimension));
      lag.rational.push_back(std::move(group));
    }
    for (const auto& value : as_array(
             lag_object.at("val"), "Nhat valuations"))
      lag.valuations.push_back(parse_validity(value));
    prepared.nhat_lags.push_back(std::move(lag));
  }
  if (!root.at("d0_inverse").is_null())
    prepared.d0_inverse_scalar =
        parse_scalar<Scalar>(root.at("d0_inverse"));
  for (const auto& raw_block : as_array(root.at("blocks"), "Jordan blocks")) {
    JordanBlock block;
    for (const auto& col : as_array(raw_block, "Jordan block"))
      block.columns.push_back(as_u32(col, "Jordan column"));
    prepared.blocks.push_back(std::move(block));
  }
  if (const auto* raw_assembly = root.if_contains("assembly");
      raw_assembly && !raw_assembly->is_null()) {
    const auto& assembly_object = as_object(*raw_assembly, "assembly matrix");
    PreparedMatrix<Scalar> matrix;
    matrix.identity = assembly_object.if_contains("identity") != nullptr &&
                      assembly_object.at("identity").as_bool();
    for (const auto& raw_matrix : as_array(
             assembly_object.at("poly"), "assembly polynomial matrices"))
      matrix.polynomial.push_back(
          parse_matrix_shift<Scalar>(raw_matrix, prepared.dimension));
    for (const auto& raw_group : as_array(
             assembly_object.at("rat"), "assembly rational groups")) {
      const auto& group_object = as_object(
          raw_group, "assembly rational group");
      RationalGroup<Scalar> group;
      group.denominator_index = as_u32(
          group_object.at("q"), "assembly denominator index");
      for (const auto& raw_matrix : as_array(
               group_object.at("num"), "assembly rational numerator"))
        group.numerator.push_back(
            parse_matrix_shift<Scalar>(raw_matrix, prepared.dimension));
      matrix.rational.push_back(std::move(group));
    }
    for (const auto& value : as_array(
             assembly_object.at("val"), "assembly valuations"))
      matrix.valuations.push_back(parse_validity(value));
    prepared.assembly_matrix = std::move(matrix);
  }
  prepared.chop_digits = root.if_contains("chop_digits")
      ? as_i32(root.at("chop_digits"), "chop digits") : 0;
  // d0 is immutable for this chart/frame.  Retain its framed inverse now so
  // every later column/source run uses the same typed kernel without repeating
  // the quadratic series inversion.
  retain_framed_d0_inverse(prepared);
  return prepared;
}

template <typename Scalar>
void parse_run_state(const json::object& run,
                     const PreparedRecurrenceOperator<Scalar>& prepared,
                     RecurrenceProblem<Scalar>& problem) {
  // These fields are deliberately all required.  A persistent solve must
  // never inherit a previous column's seed, source, or resonance schedule.
  problem.dimension = prepared.dimension;
  problem.frame_base = prepared.frame_base;
  problem.frame_width = prepared.frame_width;
  problem.nmax = as_u32(run.at("nmax"), "nmax");
  problem.log_max = as_u32(run.at("p"), "log maximum");
  problem.has_initial = run.at("has_initial").as_bool();
  problem.adaptive_lower_frame_probe = run.at("adaptive_probe").as_bool();
  problem.a_target = parse_scalar<Scalar>(run.at("a_target"));
  problem.b_target = parse_scalar<Scalar>(run.at("b_target"));
  problem.a_shift_min = as_i32(run.at("a_shift_min"), "a shift minimum");
  for (const auto& value : as_array(run.at("a_shifts"), "a shifts"))
    problem.a_shifts.push_back(parse_scalar<Scalar>(value));

  for (const auto& raw_row : as_array(run.at("schedule"), "step schedule")) {
    std::vector<BlockStep<Scalar>> row;
    for (const auto& raw_step : as_array(raw_row, "schedule row")) {
      const auto& step = as_object(raw_step, "schedule step");
      const auto kind = std::string(step.at("case").as_string());
      const auto step_case = kind == "T" ? StepCase::Taylor
          : kind == "P" ? StepCase::Pseudo
          : kind == "R" ? StepCase::Resonant
          : throw std::invalid_argument("unknown recurrence step case: " + kind);
      row.push_back({step_case, parse_scalar<Scalar>(step.at("da")),
                     parse_scalar<Scalar>(step.at("db"))});
    }
    problem.schedule.push_back(std::move(row));
  }
  for (const auto& value : as_array(run.at("initial"), "initial tensor"))
    problem.initial.push_back(parse_scalar<Scalar>(value));
  for (const auto& value : as_array(
           run.at("initial_validity"), "initial validity"))
    problem.initial_validity.push_back(parse_validity(value));

  const auto& raw_source = run.at("source");
  if (!raw_source.is_null()) {
    const auto& source_object = as_object(raw_source, "source");
    SourceData<Scalar> source;
    for (const auto& value : as_array(
             source_object.at("frames"), "source frames"))
      source.frames.push_back(parse_scalar<Scalar>(value));
    for (const auto& value : as_array(
             source_object.at("validity"), "source validity"))
      source.validity.push_back(parse_validity(value));
    for (const auto& value : as_array(
             source_object.at("present"), "source presence"))
      source.present.push_back(value.as_bool() ? 1 : 0);
    problem.source = std::move(source);
  }
  problem.return_u = run.at("return_u").as_bool();
  if (!problem.return_u && !prepared.assembly_matrix.has_value())
    throw std::invalid_argument(
        "persistent run suppresses U without a prepared assembly matrix");
}

json::value encode_validity(std::int32_t value) {
  return value == kCompleteInfinity ? json::value(nullptr) : json::value(value);
}

json::value encode_scalar(const Rational& value, int) {
  return json::string(value.str());
}

json::value encode_scalar(const SymbolicRational& value, int) {
  return json::string(value.str());
}

json::value encode_scalar(const ComplexBall& value, int digits) {
  json::array result;
  result.emplace_back(value.real_midpoint(digits));
  result.emplace_back(value.imag_midpoint(digits));
  result.emplace_back(value.real_radius_exponent());
  result.emplace_back(value.imag_radius_exponent());
  return result;
}

template <typename Scalar>
json::object encode_result(RecurrenceResult<Scalar>&& result,
                           std::optional<AssembledResult<Scalar>>&& assembled,
                           bool return_u, int digits, double elapsed_ms) {
  json::object out;
  out["status"] = "ok";
  out["elapsed_ms"] = elapsed_ms;
  out["top_valid"] = encode_validity(result.top_valid);
  if (return_u) {
    json::array u;
    u.reserve(result.u.size());
    for (const auto& value : result.u) u.push_back(encode_scalar(value, digits));
    out["u"] = std::move(u);
    json::array validity;
    validity.reserve(result.validity.size());
    for (const auto value : result.validity) validity.push_back(encode_validity(value));
    out["validity"] = std::move(validity);
  }
  if (assembled.has_value()) {
    json::object encoded_assembly;
    encoded_assembly["min"] = assembled->min_power;
    encoded_assembly["max"] = assembled->complete_max;
    json::array coefficients;
    coefficients.reserve(assembled->coefficients.size());
    for (const auto& value : assembled->coefficients)
      coefficients.push_back(encode_scalar(value, digits));
    encoded_assembly["coefficients"] = std::move(coefficients);
    out["assembled"] = std::move(encoded_assembly);
  }
  json::array hits;
  for (const auto& hit : result.hits) {
    json::object encoded;
    encoded["n"] = hit.n;
    json::array cols;
    for (const auto col : hit.columns) cols.emplace_back(col);
    encoded["cols"] = std::move(cols);
    encoded["delta_b"] = encode_scalar(hit.delta_b, digits);
    json::array frames;
    for (const auto& row : hit.gamma_frames)
      for (const auto& value : row) frames.push_back(encode_scalar(value, digits));
    encoded["frames"] = std::move(frames);
    json::array hit_validity;
    for (const auto value : hit.gamma_validity)
      hit_validity.push_back(encode_validity(value));
    encoded["validity"] = std::move(hit_validity);
    hits.push_back(std::move(encoded));
  }
  out["hits"] = std::move(hits);
  return out;
}

template <typename Scalar>
json::object run_problem(RecurrenceProblem<Scalar>& problem, int digits) {
  const auto started = std::chrono::steady_clock::now();
  auto result = RecurrenceSolver<Scalar>(problem).run();
  std::optional<AssembledResult<Scalar>> assembled;
  if (problem.assembly_matrix.has_value())
    assembled = assemble_recurrence(problem, result);
  const auto ended = std::chrono::steady_clock::now();
  const double elapsed = std::chrono::duration<double, std::milli>(ended - started).count();
  return encode_result(std::move(result), std::move(assembled),
                       problem.return_u, digits, elapsed);
}

template <typename Scalar>
json::object run_prepared_problem(
    const PreparedRecurrenceOperator<Scalar>& prepared,
    RecurrenceProblem<Scalar>& problem, int digits) {
  const auto started = std::chrono::steady_clock::now();
  auto result = RecurrenceSolver<Scalar>(problem, prepared).run();
  std::optional<AssembledResult<Scalar>> assembled;
  if (prepared.assembly_matrix.has_value())
    assembled = assemble_recurrence(prepared, problem, result);
  const auto ended = std::chrono::steady_clock::now();
  const double elapsed =
      std::chrono::duration<double, std::milli>(ended - started).count();
  return encode_result(std::move(result), std::move(assembled),
                       problem.return_u, digits, elapsed);
}

template <typename Scalar>
json::object run_typed(const json::object& root, int digits) {
  auto problem = parse_problem<Scalar>(root);
  return run_problem(problem, digits);
}

struct SCCCertificate {
  std::uint32_t component_count = 0;
  std::uint32_t coupling_depth = 0;
  std::vector<std::vector<std::uint32_t>> components;
  std::vector<std::uint32_t> component_of;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> structural_edges;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> condensation_edges;
  std::vector<std::uint32_t> topological_order;
  std::string exact_record;
};

std::pair<std::uint32_t, std::uint32_t> parse_edge(
    const json::value& value, std::uint32_t bound, const char* label) {
  const auto& edge = as_array(value, label);
  if (edge.size() != 2)
    throw std::invalid_argument(std::string(label) + " must have two vertices");
  const auto source = as_u32(edge[0], "edge source");
  const auto target = as_u32(edge[1], "edge target");
  if (source >= bound || target >= bound)
    throw std::invalid_argument(std::string(label) + " vertex is outside range");
  return {source, target};
}

SCCCertificate validate_scc_certificate(const json::value& raw,
                                        std::uint32_t dimension) {
  const auto& object = as_object(raw, "SCC certificate");
  const auto& raw_components = as_array(
      object.at("components"), "SCC components");
  if (raw_components.empty())
    throw std::invalid_argument("SCC certificate has no components");

  std::vector<std::vector<std::uint32_t>> components;
  std::vector<std::uint32_t> component_of(dimension,
      std::numeric_limits<std::uint32_t>::max());
  for (std::uint32_t block = 0; block < raw_components.size(); ++block) {
    std::vector<std::uint32_t> component;
    for (const auto& raw_vertex : as_array(
             raw_components[block], "SCC component")) {
      const auto vertex = as_u32(raw_vertex, "SCC vertex");
      if (vertex >= dimension ||
          component_of[vertex] != std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument(
            "SCC components must partition the recurrence dimension");
      component_of[vertex] = block;
      component.push_back(vertex);
    }
    if (component.empty())
      throw std::invalid_argument("SCC certificate contains an empty component");
    std::sort(component.begin(), component.end());
    components.push_back(std::move(component));
  }
  if (std::any_of(component_of.begin(), component_of.end(), [](auto value) {
        return value == std::numeric_limits<std::uint32_t>::max();
      }))
    throw std::invalid_argument(
        "SCC components do not cover the recurrence dimension");

  std::vector<std::vector<std::uint32_t>> adjacency(dimension), reverse(dimension);
  std::set<std::pair<std::uint32_t, std::uint32_t>> structural_edges;
  for (const auto& raw_edge : as_array(
           object.at("structural_edges"), "exact structural edges")) {
    const auto edge = parse_edge(raw_edge, dimension, "exact structural edge");
    if (structural_edges.insert(edge).second) {
      adjacency[edge.first].push_back(edge.second);
      reverse[edge.second].push_back(edge.first);
    }
  }

  // Recompute SCCs from the exact structural nonzero graph.  The supplied
  // certificate is accepted only when it is the same partition, not merely
  // a topologically plausible grouping.
  std::vector<std::uint8_t> seen(dimension, 0);
  std::vector<std::uint32_t> order;
  std::function<void(std::uint32_t)> visit = [&](std::uint32_t vertex) {
    seen[vertex] = 1;
    for (const auto next : adjacency[vertex]) if (!seen[next]) visit(next);
    order.push_back(vertex);
  };
  for (std::uint32_t vertex = 0; vertex < dimension; ++vertex)
    if (!seen[vertex]) visit(vertex);
  std::fill(seen.begin(), seen.end(), 0);
  std::vector<std::vector<std::uint32_t>> derived_components;
  std::function<void(std::uint32_t, std::vector<std::uint32_t>&)> collect =
      [&](std::uint32_t vertex, std::vector<std::uint32_t>& component) {
        seen[vertex] = 1;
        component.push_back(vertex);
        for (const auto next : reverse[vertex])
          if (!seen[next]) collect(next, component);
      };
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    if (seen[*it]) continue;
    std::vector<std::uint32_t> component;
    collect(*it, component);
    std::sort(component.begin(), component.end());
    derived_components.push_back(std::move(component));
  }
  auto canonical_components = [](auto value) {
    std::sort(value.begin(), value.end());
    return value;
  };
  if (canonical_components(components) !=
      canonical_components(derived_components))
    throw std::invalid_argument(
        "SCC partition does not match the exact structural graph");

  const auto block_count = static_cast<std::uint32_t>(components.size());
  std::set<std::pair<std::uint32_t, std::uint32_t>> derived_condensation;
  for (const auto [source, target] : structural_edges) {
    const auto from = component_of[source];
    const auto to = component_of[target];
    if (from != to) derived_condensation.insert({from, to});
  }
  std::set<std::pair<std::uint32_t, std::uint32_t>> supplied_condensation;
  for (const auto& raw_edge : as_array(
           object.at("condensation_edges"), "condensation edges"))
    supplied_condensation.insert(parse_edge(
        raw_edge, block_count, "condensation edge"));
  if (derived_condensation != supplied_condensation)
    throw std::invalid_argument(
        "SCC condensation edges do not match the exact structural graph");

  std::vector<std::uint32_t> topological_order;
  std::vector<std::uint32_t> topological_position(
      block_count, std::numeric_limits<std::uint32_t>::max());
  for (const auto& raw_block : as_array(
           object.at("topological_order"), "SCC topological order")) {
    const auto block = as_u32(raw_block, "SCC topological block");
    if (block >= block_count ||
        topological_position[block] != std::numeric_limits<std::uint32_t>::max())
      throw std::invalid_argument("invalid SCC topological order");
    topological_position[block] =
        static_cast<std::uint32_t>(topological_order.size());
    topological_order.push_back(block);
  }
  if (topological_order.size() != block_count)
    throw std::invalid_argument("SCC topological order is incomplete");
  for (const auto [source, target] : derived_condensation)
    if (topological_position[source] >= topological_position[target])
      throw std::invalid_argument(
          "SCC topological order contradicts a condensation edge");

  std::vector<std::uint32_t> depth(block_count, 0);
  for (const auto block : topological_order) {
    for (const auto [source, target] : derived_condensation)
      if (source == block) depth[target] = std::max(depth[target], depth[source] + 1);
  }
  const auto coupling_depth = depth.empty()
      ? 0 : *std::max_element(depth.begin(), depth.end());
  if (as_u32(object.at("coupling_depth"), "SCC coupling depth") !=
      coupling_depth)
    throw std::invalid_argument(
        "SCC coupling depth does not match the exact condensation graph");

  SCCCertificate certificate;
  certificate.component_count = block_count;
  certificate.coupling_depth = coupling_depth;
  certificate.components = std::move(components);
  certificate.component_of = std::move(component_of);
  certificate.structural_edges.assign(
      structural_edges.begin(), structural_edges.end());
  certificate.condensation_edges.assign(
      derived_condensation.begin(), derived_condensation.end());
  certificate.topological_order = std::move(topological_order);
  certificate.exact_record = json::serialize(raw);
  return certificate;
}

json::array encode_indices(const std::vector<std::uint32_t>& values) {
  json::array encoded;
  encoded.reserve(values.size());
  for (const auto value : values) encoded.push_back(value);
  return encoded;
}

std::vector<std::string> parse_symbols(const json::object& object) {
  std::vector<std::string> symbols;
  if (const auto* raw_symbols = object.if_contains("symbols")) {
    for (const auto& value : as_array(*raw_symbols, "symbolic variables")) {
      if (!value.is_string())
        throw std::invalid_argument("symbolic variable names must be strings");
      symbols.emplace_back(value.as_string());
    }
  }
  return symbols;
}

std::string required_string(const json::object& object, const char* key) {
  const auto& value = object.at(key);
  if (!value.is_string() || value.as_string().empty())
    throw std::invalid_argument(std::string(key) + " must be a nonempty string");
  return std::string(value.as_string());
}

std::string canonical_chart_geometry_record(const json::value& raw) {
  const auto& geometry = as_object(raw, "exact chart geometry");
  const auto center = required_string(geometry, "center_exact");
  const auto scale = required_string(geometry, "scale_exact");
  if (scale == "0")
    throw std::invalid_argument(
        "exact chart geometry requires a nonzero scale");
  const auto infinite = geometry.at("infinite_radius").as_bool();

  json::array prescriptions;
  for (const auto& raw_prescription : as_array(
           geometry.at("prescriptions"), "exact chart prescriptions")) {
    const auto& prescription = as_object(
        raw_prescription, "exact chart prescription");
    const auto factor = required_string(prescription, "factor_exact");
    const auto sign = as_i32(prescription.at("sign"), "prescription sign");
    const auto multiplicity = as_u32(
        prescription.at("multiplicity"), "prescription multiplicity");
    const auto leading = as_i32(
        prescription.at("leading_coefficient_sign"),
        "prescription leading coefficient sign");
    if ((sign != -1 && sign != 1) || multiplicity == 0 ||
        (leading != -1 && leading != 1))
      throw std::invalid_argument(
          "malformed exact chart prescription");
    prescriptions.push_back(json::object{
        {"factor_exact", factor}, {"sign", sign},
        {"multiplicity", multiplicity},
        {"leading_coefficient_sign", leading}});
  }

  json::object canonical{{"center_exact", center},
                         {"scale_exact", scale},
                         {"infinite_radius", infinite}};
  if (infinite) {
    if (const auto* radius = geometry.if_contains("radius_exact");
        radius != nullptr && !radius->is_null() &&
        (!radius->is_string() || radius->as_string() != "Infinity"))
      throw std::invalid_argument(
          "infinite exact chart geometry radius identity must be Infinity");
    // The Wolfram producer historically emitted radius_exact even for an
    // infinite chart.  Canonicalize both its explicit form and an omitted
    // field to the same exact record so composite preparation remains
    // backward compatible without accepting a finite-radius mismatch.
    canonical["radius_exact"] = "Infinity";
  } else {
    canonical["radius_exact"] = required_string(geometry, "radius_exact");
  }
  canonical["prescriptions"] = std::move(prescriptions);
  return json::serialize(canonical);
}

void validate_first_slice_rational_geometry(const json::value& raw) {
  const auto& geometry = as_object(raw, "exact chart geometry");
  try {
    (void)Rational(required_string(geometry, "center_exact"));
  } catch (const std::invalid_argument&) {
    throw std::invalid_argument(
        "native SCC first slice requires an exact rational chart center");
  }
  Rational scale;
  try {
    scale = Rational(required_string(geometry, "scale_exact"));
  } catch (const std::invalid_argument&) {
    throw std::invalid_argument(
        "native SCC first slice requires an exact rational chart scale");
  }
  if (scale.is_zero())
    throw std::invalid_argument(
        "native SCC first slice requires a nonzero chart scale");
  if (geometry.at("infinite_radius").as_bool())
    throw std::invalid_argument(
        "native SCC first slice requires a positive finite rational radius");

  Rational radius;
  try {
    radius = Rational(required_string(geometry, "radius_exact"));
  } catch (const std::invalid_argument&) {
    throw std::invalid_argument(
        "native SCC first slice requires an exact rational finite radius");
  }
  const auto canonical_radius = radius.str();
  if (radius.is_zero() || canonical_radius.front() == '-')
    throw std::invalid_argument(
        "native SCC first slice requires a positive finite radius");
}

struct RetainedCompositeGeometry {
  ChartGeometry chart;
  std::string radius_exact;
  std::vector<Prescription> prescriptions;
};

RetainedCompositeGeometry parse_retained_composite_geometry(
    const json::value& raw) {
  const auto& geometry = as_object(raw, "exact chart geometry");
  RetainedCompositeGeometry retained;
  retained.chart.center_exact = required_string(geometry, "center_exact");
  retained.chart.scale_exact = required_string(geometry, "scale_exact");
  retained.chart.infinite_radius = geometry.at("infinite_radius").as_bool();
  retained.radius_exact = retained.chart.infinite_radius
      ? "Infinity"
      : Rational(required_string(geometry, "radius_exact")).str();
  if (!retained.chart.infinite_radius)
    retained.chart.radius = ComplexBall::from_strings(
        retained.radius_exact);
  for (const auto& raw_prescription : as_array(
           geometry.at("prescriptions"), "exact chart prescriptions")) {
    const auto& prescription = as_object(
        raw_prescription, "exact chart prescription");
    retained.prescriptions.push_back(Prescription{
        required_string(prescription, "factor_exact"),
        as_i32(prescription.at("sign"), "prescription sign"),
        as_u32(prescription.at("multiplicity"),
               "prescription multiplicity"),
        as_i32(prescription.at("leading_coefficient_sign"),
               "prescription leading coefficient sign")});
  }
  return retained;
}

std::string canonical_native_scc_capabilities(const json::value& raw) {
  const auto& capabilities = as_object(
      raw, "native SCC chart capabilities");
  return json::serialize(json::object{
      {"regular", capabilities.at("regular").as_bool()},
      {"identity_gauge", capabilities.at("identity_gauge").as_bool()},
      {"identity_v", capabilities.at("identity_v").as_bool()},
      {"no_pseudo", capabilities.at("no_pseudo").as_bool()}});
}

TruthValue parse_truth_value(const json::value& value, const char* label) {
  if (!value.is_string())
    throw std::invalid_argument(std::string(label) + " must be yes, no or unknown");
  const auto text = std::string(value.as_string());
  if (text == "yes") return TruthValue::Yes;
  if (text == "no") return TruthValue::No;
  if (text == "unknown") return TruthValue::Unknown;
  throw std::invalid_argument(std::string(label) + " must be yes, no or unknown");
}

ExactSign parse_exact_sign(const json::value& value, const char* label) {
  if (!value.is_string())
    throw std::invalid_argument(
        std::string(label) + " must be negative, zero, positive or unknown");
  const auto text = std::string(value.as_string());
  if (text == "negative") return ExactSign::Negative;
  if (text == "zero") return ExactSign::Zero;
  if (text == "positive") return ExactSign::Positive;
  if (text == "unknown") return ExactSign::Unknown;
  throw std::invalid_argument(
      std::string(label) + " must be negative, zero, positive or unknown");
}

void verify_optional_truth(const json::object& object, const char* key,
                           TruthValue expected) {
  if (const auto* raw = object.if_contains(key);
      raw != nullptr && parse_truth_value(*raw, key) != expected)
    throw std::invalid_argument(
        std::string("contradictory rational exact-tag fact: ") + key);
}

void verify_optional_sign(const json::object& object, const char* key,
                          ExactSign expected) {
  if (const auto* raw = object.if_contains(key);
      raw != nullptr && parse_exact_sign(*raw, key) != expected)
    throw std::invalid_argument(
        std::string("contradictory rational exact-tag fact: ") + key);
}

ExactScalarDescriptor parse_exact_descriptor(const json::value& raw,
                                              const char* label) {
  const auto& object = as_object(raw, label);
  const auto domain = required_string(object, "domain");
  const auto canonical = required_string(object, "canonical");
  if (domain == "rational") {
    auto descriptor = ExactScalarDescriptor::rational(canonical);
    verify_optional_truth(object, "is_zero", descriptor.is_zero);
    verify_optional_truth(object, "is_integer", descriptor.is_integer);
    verify_optional_sign(object, "sign", descriptor.sign);
    if (const auto* raw_specialization = object.if_contains("specialization")) {
      const auto specialization = parse_scalar<ComplexBall>(*raw_specialization);
      if (!acb_equal(specialization.raw(), descriptor.numeric().raw()))
        throw std::invalid_argument(
            std::string(label) + " specialization contradicts exact rational tag");
    }
    return descriptor;
  }

  const auto zero = parse_truth_value(object.at("is_zero"), "is_zero");
  const auto integer = parse_truth_value(object.at("is_integer"), "is_integer");
  const auto sign = parse_exact_sign(object.at("sign"), "sign");
  std::optional<ComplexBall> specialization;
  if (const auto* raw_specialization = object.if_contains("specialization"))
    specialization = parse_scalar<ComplexBall>(*raw_specialization);
  if (domain == "symbolic-rational") {
    std::vector<std::string> symbols;
    for (const auto& symbol : as_array(object.at("symbols"), "tag symbols")) {
      if (!symbol.is_string() || symbol.as_string().empty())
        throw std::invalid_argument("tag symbols must be nonempty strings");
      symbols.emplace_back(symbol.as_string());
    }
    return ExactScalarDescriptor::symbolic(
        canonical, std::move(symbols), zero, integer, sign,
        std::move(specialization));
  }
  if (domain == "algebraic") {
    if (!specialization.has_value())
      throw std::invalid_argument(
          std::string(label) + " algebraic tag requires a specialization");
    return ExactScalarDescriptor::algebraic(
        canonical, zero, integer, sign, std::move(*specialization));
  }
  throw std::invalid_argument(
      std::string(label) + " has unsupported exact domain: " + domain);
}

const char* encode_truth_value(TruthValue value) {
  if (value == TruthValue::Yes) return "yes";
  if (value == TruthValue::No) return "no";
  return "unknown";
}

const char* encode_exact_sign(ExactSign value) {
  if (value == ExactSign::Negative) return "negative";
  if (value == ExactSign::Zero) return "zero";
  if (value == ExactSign::Positive) return "positive";
  return "unknown";
}

json::object encode_exact_descriptor(const ExactScalarDescriptor& descriptor) {
  const char* domain = descriptor.domain == ExactDomain::Rational
      ? "rational"
      : descriptor.domain == ExactDomain::SymbolicRational
          ? "symbolic-rational" : "algebraic";
  json::array symbols;
  for (const auto& symbol : descriptor.symbols) symbols.emplace_back(symbol);
  return json::object{{"domain", domain}, {"canonical", descriptor.canonical},
                      {"symbols", std::move(symbols)},
                      {"is_zero", encode_truth_value(descriptor.is_zero)},
                      {"is_integer", encode_truth_value(descriptor.is_integer)},
                      {"sign", encode_exact_sign(descriptor.sign)},
                      {"has_specialization",
                       descriptor.specialization.has_value()}};
}

struct LocalMetadata {
  ChartGeometry chart;
  ExactScalarDescriptor a;
  ExactScalarDescriptor b;
  std::vector<Prescription> prescriptions;
  std::string checkpoint_identity;
};

LocalMetadata parse_local_metadata(const json::object& metadata) {
  LocalMetadata out;
  const auto& chart = as_object(metadata.at("chart"), "local chart metadata");
  out.chart.center_exact = required_string(chart, "center_exact");
  out.chart.scale_exact = required_string(chart, "scale_exact");
  out.chart.infinite_radius = chart.at("infinite_radius").as_bool();
  if (const auto* radius = chart.if_contains("radius"))
    out.chart.radius = parse_scalar<ComplexBall>(*radius);
  else if (!out.chart.infinite_radius)
    throw std::invalid_argument("finite local chart requires a radius");
  if (!out.chart.infinite_radius &&
      (!local_detail::exactly_real(out.chart.radius) ||
       !arb_is_positive(acb_realref(out.chart.radius.raw()))))
    throw std::invalid_argument(
        "finite local chart radius must be a provably positive real ball");

  const auto& tag = as_object(metadata.at("tag"), "local exact tag");
  out.a = parse_exact_descriptor(tag.at("a"), "local a tag");
  out.b = parse_exact_descriptor(tag.at("b"), "local b tag");
  for (const auto* descriptor : {&out.a, &out.b})
    if (descriptor->specialization.has_value() &&
        !local_detail::exactly_real(*descriptor->specialization))
      throw std::invalid_argument(
          "local exact-tag specialization must be real");
  for (const auto& raw_prescription : as_array(
           metadata.at("prescriptions"), "local prescriptions")) {
    const auto& prescription = as_object(
        raw_prescription, "local prescription");
    Prescription parsed{
        required_string(prescription, "factor_exact"),
        as_i32(prescription.at("sign"), "prescription sign"),
        as_u32(prescription.at("multiplicity"), "prescription multiplicity"),
        as_i32(prescription.at("leading_coefficient_sign"),
               "prescription leading coefficient sign")};
    if ((parsed.sign != -1 && parsed.sign != 1) ||
        (parsed.leading_coefficient_sign != -1 &&
         parsed.leading_coefficient_sign != 1) ||
        parsed.multiplicity == 0)
      throw std::invalid_argument(
          "malformed analytic-continuation prescription");
    out.prescriptions.push_back(std::move(parsed));
  }
  out.checkpoint_identity = required_string(metadata, "checkpoint_identity");
  return out;
}

template <typename Scalar>
void verify_tag_binding(const ExactScalarDescriptor& descriptor,
                        const Scalar& target, const char* label);

template <>
void verify_tag_binding<Rational>(const ExactScalarDescriptor& descriptor,
                                  const Rational& target,
                                  const char* label) {
  if (descriptor.domain != ExactDomain::Rational ||
      !(Rational(descriptor.canonical) == target))
    throw std::invalid_argument(
        std::string(label) + " exact tag does not equal recurrence target");
}

template <>
void verify_tag_binding<SymbolicRational>(
    const ExactScalarDescriptor& descriptor,
    const SymbolicRational& target, const char* label) {
  if (descriptor.domain == ExactDomain::Algebraic ||
      !(SymbolicRational(descriptor.canonical) == target))
    throw std::invalid_argument(
        std::string(label) + " exact tag does not equal recurrence target");
}

template <>
void verify_tag_binding<ComplexBall>(const ExactScalarDescriptor& descriptor,
                                     const ComplexBall& target,
                                     const char* label) {
  if (!descriptor.specialization.has_value())
    throw std::invalid_argument(
        std::string(label) + " exact tag needs a numeric specialization");
  if (!acb_equal(descriptor.specialization->raw(), target.raw()))
    throw std::invalid_argument(
        std::string(label) + " specialization does not equal recurrence target");
}

struct ChartStats {
  std::uint64_t runs = 0;
  std::uint64_t local_runs = 0;
  double prepare_parse_ms = 0.0;
  double run_parse_ms = 0.0;
  double kernel_ms = 0.0;
  double local_run_parse_ms = 0.0;
  double local_kernel_ms = 0.0;
};

class StoredLocalBase;
class StoredTilePlan;

class PreparedChartBase {
 public:
  PreparedChartBase(std::string handle, std::string key,
                    std::string exact_identity, std::string signature,
                    std::optional<std::string> geometry_record,
                    std::optional<std::string> principal_matrix_record,
                    std::optional<std::string> native_scc_capabilities,
                    SCCCertificate scc, double prepare_parse_ms)
      : handle_(std::move(handle)), key_(std::move(key)),
        exact_identity_(std::move(exact_identity)),
        signature_(std::move(signature)),
        geometry_record_(std::move(geometry_record)),
        principal_matrix_record_(std::move(principal_matrix_record)),
        native_scc_capabilities_(std::move(native_scc_capabilities)),
        scc_(std::move(scc)),
        prepare_parse_ms_(prepare_parse_ms) {}
  virtual ~PreparedChartBase() = default;

  virtual json::object solve(const json::object& run, int output_digits) = 0;
  virtual std::shared_ptr<StoredLocalBase> solve_local(
      const std::string& local_handle, const json::object& run,
      const json::object& metadata) = 0;
  virtual std::uint32_t dimension() const = 0;
  virtual std::int32_t frame_base() const = 0;
  virtual std::uint32_t frame_width() const = 0;
  virtual const char* d0_inverse_mode() const = 0;
  virtual ChartStats stats() const = 0;

  const std::string& handle() const { return handle_; }
  const std::string& key() const { return key_; }
  const std::string& exact_identity() const { return exact_identity_; }
  const std::string& signature() const { return signature_; }
  const std::optional<std::string>& geometry_record() const {
    return geometry_record_;
  }
  const std::optional<std::string>& principal_matrix_record() const {
    return principal_matrix_record_;
  }
  const std::optional<std::string>& native_scc_capabilities() const {
    return native_scc_capabilities_;
  }
  const SCCCertificate& scc() const { return scc_; }

 protected:
  std::string handle_;
  std::string key_;
  std::string exact_identity_;
  std::string signature_;
  std::optional<std::string> geometry_record_;
  std::optional<std::string> principal_matrix_record_;
  std::optional<std::string> native_scc_capabilities_;
  SCCCertificate scc_;
  double prepare_parse_ms_ = 0.0;
};

std::recursive_mutex& symbolic_run_mutex() {
  static std::recursive_mutex mutex;
  return mutex;
}

template <typename Scalar, typename Object, typename... Arguments>
std::shared_ptr<Object> make_retained_typed_shared(Arguments&&... arguments) {
  if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
    // SymbolicRational destruction touches the process-global FLINT context
    // and live-object count.  The deleter travels with the shared_ptr control
    // block through base conversions and dynamic casts, so releases, session
    // close, library reset, and delayed in-flight holders all serialize.  A
    // recursive mutex is required because destroying a composite can release
    // its last typed chart pointers under this same deleter.
    std::lock_guard<std::recursive_mutex> construction_lock(
        symbolic_run_mutex());
    return std::shared_ptr<Object>(
        new Object(std::forward<Arguments>(arguments)...),
        [](Object* object) {
          std::lock_guard<std::recursive_mutex> lock(symbolic_run_mutex());
          delete object;
        });
  } else {
    return std::make_shared<Object>(
        std::forward<Arguments>(arguments)...);
  }
}

struct AcbExecutionState {
  std::mutex mutex;
  std::condition_variable changed;
  std::optional<slong> precision_bits;
  std::size_t active = 0;
};

AcbExecutionState& acb_execution_state() {
  static AcbExecutionState state;
  return state;
}

class AcbPrecisionLease {
 public:
  explicit AcbPrecisionLease(slong precision_bits)
      : state_(acb_execution_state()) {
    std::unique_lock<std::mutex> lock(state_.mutex);
    state_.changed.wait(lock, [&] {
      return state_.active == 0 || state_.precision_bits == precision_bits;
    });
    if (state_.active == 0) state_.precision_bits = precision_bits;
    ++state_.active;
  }
  AcbPrecisionLease(const AcbPrecisionLease&) = delete;
  AcbPrecisionLease& operator=(const AcbPrecisionLease&) = delete;
  ~AcbPrecisionLease() {
    std::lock_guard<std::mutex> lock(state_.mutex);
    --state_.active;
    if (state_.active == 0) state_.changed.notify_all();
  }

 private:
  AcbExecutionState& state_;
};

json::object encode_epsilon_vector(const EpsilonVector& vector, int digits) {
  json::array coefficients;
  coefficients.reserve(vector.coefficients.size());
  for (const auto& coefficient : vector.coefficients)
    coefficients.push_back(encode_scalar(coefficient, digits));
  json::object encoded{{"min", vector.epsilon.min_power},
                       {"max", vector.epsilon.complete_max},
                       {"dimension", vector.dimension},
                       {"coefficients", std::move(coefficients)}};
  if (!vector.error.empty()) {
    json::array upper;
    upper.reserve(vector.error.absolute.size());
    for (const auto& bound : vector.error.absolute)
      upper.push_back(bound.approximate_upper());
    const char* guarantee = vector.error.guarantee == ErrorGuarantee::Certified
        ? "certified"
        : vector.error.guarantee == ErrorGuarantee::Advisory
            ? "advisory" : "none";
    encoded["error"] = json::object{
        {"min", vector.error.frame.min_power},
        {"max", vector.error.frame.complete_max},
        {"guarantee", guarantee},
        {"absolute_upper_approx", std::move(upper)},
        {"bound_encoding", "approximate-double"},
        {"provenance", vector.error.provenance}};
  }
  return encoded;
}

constexpr const char* kRegularTailMajorantCapability =
    "retained-regular-homogeneous-gronwall-cauchy-tail-v1";

const char* tail_majorant_status_name(TailMajorantStatus status) {
  switch (status) {
    case TailMajorantStatus::Certified:
      return "certified";
    case TailMajorantStatus::Inconclusive:
      return "inconclusive";
    case TailMajorantStatus::Unsupported:
      return "unsupported";
  }
  throw std::logic_error("unknown tail-majorant status");
}

RegularTaylorTailModelResult unavailable_tail_model(std::string detail) {
  return {TailMajorantStatus::Unsupported, std::nullopt,
          std::move(detail)};
}

struct TailModelCheckpointMarker {
  std::string saved_status;
  bool attached_before_save = false;
};

json::object encode_tail_model_status(
    const RegularTaylorTailModelResult& result) {
  json::object encoded{
      {"capability", kRegularTailMajorantCapability},
      {"status", tail_majorant_status_name(result.status)},
      {"attached", result.model.has_value()},
      {"detail", result.detail},
      {"checkpoint_serialized", result.model.has_value()}};
  if (result.model.has_value()) {
    encoded["operator_identity"] = result.model->operator_identity;
    encoded["local_checkpoint_identity"] =
        result.model->local_checkpoint_identity;
    encoded["epsilon"] = json::object{
        {"min", result.model->epsilon.min_power},
        {"max", result.model->epsilon.complete_max}};
    encoded["taylor_complete_max"] =
        result.model->taylor_complete_max;
    encoded["provenance"] = result.model->provenance;
  }
  return encoded;
}

std::optional<std::string> parse_certified_tail_witness(
    const json::object& request) {
  const auto* raw_options = request.if_contains("options");
  if (raw_options == nullptr) return std::nullopt;
  const auto& options = as_object(*raw_options, "local evaluation options");
  const auto* raw_radius =
      options.if_contains("certified_tail_radius_exact");
  if (raw_radius == nullptr || raw_radius->is_null()) return std::nullopt;
  if (!raw_radius->is_string() || raw_radius->as_string().empty())
    throw std::invalid_argument(
        "certified tail radius must be a nonempty exact rational string");
  const std::string radius(raw_radius->as_string());
  const Rational parsed(radius);
  if (parsed.sign() <= 0)
    throw std::invalid_argument("certified tail radius must be positive");
  return parsed.str();
}

json::object encode_point_tail_certificate(
    const RegularTaylorPointTailCertificate& certificate,
    const std::string& witness_radius) {
  json::object result{
      {"capability", kRegularTailMajorantCapability},
      {"requested", true},
      {"status", tail_majorant_status_name(certificate.status)},
      {"witness_radius_exact", witness_radius},
      {"detail", certificate.detail}};
  if (certificate.disk.status == TailMajorantStatus::Certified) {
    result["q_lower_approx"] =
        certificate.disk.q_lower.approximate_upper();
    result["ode_norm_upper_approx"] =
        certificate.disk.ode_norm_upper.approximate_upper();
    result["bound_encoding"] = "approximate-double-diagnostics";
  }
  return result;
}

EvaluationOptions parse_local_evaluation_options(
    const json::object& request, bool default_tail_estimate) {
  EvaluationOptions options;
  options.compute_tail_estimate = default_tail_estimate;
  if (const auto* raw_options = request.if_contains("options")) {
    const auto& object = as_object(*raw_options, "local evaluation options");
    if (const auto* raw_sign = object.if_contains("imaginary_sign");
        raw_sign != nullptr && !raw_sign->is_null()) {
      const auto sign = as_i32(*raw_sign, "imaginary sign");
      if (sign != -1 && sign != 1)
        throw std::invalid_argument("imaginary sign must be +1 or -1");
      options.imaginary_sign = sign;
    }
    if (const auto* reduction = object.if_contains("t_order_reduction"))
      options.t_order_reduction = as_u32(
          *reduction, "Taylor-order reduction");
    if (const auto* tail = object.if_contains("tail_estimate"))
      options.compute_tail_estimate = tail->as_bool();
  }
  return options;
}

RealEvaluationPoint parse_local_evaluation_point(const json::object& request) {
  const auto& point_object = as_object(
      request.at("point"), "local evaluation point");
  return RealEvaluationPoint::rational(
      required_string(point_object, "exact"));
}

EpsilonWindow parse_epsilon_window(const json::object& object,
                                   const char* label) {
  EpsilonWindow window{as_i32(object.at("min"), label),
                       as_i32(object.at("max"), label)};
  (void)window.width();
  return window;
}

std::size_t checked_flat_count(std::size_t left, std::size_t right,
                               const char* label) {
  if (right != 0 && left > std::numeric_limits<std::size_t>::max() / right)
    throw std::overflow_error(std::string(label) + " size overflow");
  return left * right;
}

EpsilonMatrix parse_epsilon_matrix(const json::value& raw,
                                   const char* label) {
  const auto& object = as_object(raw, label);
  EpsilonMatrix matrix;
  matrix.epsilon = parse_epsilon_window(object, label);
  matrix.dimension = as_u32(object.at("dimension"), label);
  if (matrix.dimension == 0)
    throw std::invalid_argument(std::string(label) +
                                " dimension must be positive");
  const auto matrix_size = checked_flat_count(
      matrix.dimension, matrix.dimension, label);
  const auto expected = checked_flat_count(
      matrix.epsilon.width(), matrix_size, label);
  const auto& coefficients = as_array(object.at("coefficients"), label);
  if (coefficients.size() != expected)
    throw std::invalid_argument(std::string(label) +
                                " coefficient tensor has the wrong size");
  matrix.coefficients.reserve(expected);
  for (const auto& coefficient : coefficients)
    matrix.coefficients.push_back(parse_scalar<ComplexBall>(coefficient));
  return matrix;
}

EpsilonVector parse_epsilon_vector(const json::value& raw,
                                   const char* label) {
  const auto& object = as_object(raw, label);
  EpsilonVector vector;
  vector.epsilon = parse_epsilon_window(object, label);
  vector.dimension = as_u32(object.at("dimension"), label);
  if (vector.dimension == 0)
    throw std::invalid_argument(std::string(label) +
                                " dimension must be positive");
  const auto expected = checked_flat_count(
      vector.epsilon.width(), vector.dimension, label);
  const auto& coefficients = as_array(object.at("coefficients"), label);
  if (coefficients.size() != expected)
    throw std::invalid_argument(std::string(label) +
                                " coefficient tensor has the wrong size");
  vector.coefficients.reserve(expected);
  for (const auto& coefficient : coefficients)
    vector.coefficients.push_back(parse_scalar<ComplexBall>(coefficient));
  return vector;
}

json::array encode_magnitude_diagnostics(
    const std::vector<Magnitude>& magnitudes) {
  json::array encoded;
  encoded.reserve(magnitudes.size());
  for (const auto& magnitude : magnitudes)
    encoded.push_back(magnitude.is_finite()
                          ? json::value(magnitude.approximate_upper())
                          : json::value(nullptr));
  return encoded;
}

struct StoredLocalStats {
  std::uint64_t evaluations = 0;
  std::uint64_t residual_certifications = 0;
  std::uint64_t endpoint_limits = 0;
  std::uint64_t line_integrations = 0;
  double evaluate_ms = 0.0;
  double residual_certify_ms = 0.0;
  double endpoint_limit_ms = 0.0;
  double line_integration_ms = 0.0;
  double create_parse_ms = 0.0;
  double create_kernel_ms = 0.0;
  std::size_t coefficient_count = 0;
  std::uint64_t tail_certificate_requests = 0;
  std::uint64_t tail_certificate_certified = 0;
  std::uint64_t tail_certificate_inconclusive = 0;
  std::uint64_t tail_certificate_unsupported = 0;
};

struct NativeLocalDiagnostics {
  std::int32_t top_valid = kCompleteInfinity;
  double parse_ms = 0.0;
  double kernel_ms = 0.0;
  std::uint64_t pseudo_hits = 0;
  std::uint64_t pseudo_compensations = 0;
  std::uint32_t max_pseudo_depth = 0;
  bool pseudo_value_certified = true;
};

void require_exact_keys(const json::object& object,
                        std::initializer_list<std::string_view> expected,
                        const char* label);

json::object checkpoint_ball_record(const ComplexBall& value) {
  const auto dump = checkpoint::dump_complex_ball_exact(value);
  return json::object{{"real", dump.real}, {"imaginary", dump.imaginary}};
}

ComplexBall parse_checkpoint_ball(const json::value& raw,
                                  const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(object, {"real", "imaginary"}, label);
  return checkpoint::load_complex_ball_exact(
      {required_string(object, "real"),
       required_string(object, "imaginary")});
}

Magnitude parse_checkpoint_magnitude(const json::value& raw,
                                     const char* label) {
  if (!raw.is_string())
    throw std::invalid_argument(std::string(label) +
                                " must be an exact dump string");
  auto result = Magnitude::from_exact_dump(std::string(raw.as_string()));
  if (!result.is_finite())
    throw std::invalid_argument(std::string(label) +
                                " must be finite");
  return result;
}

template <typename Scalar>
json::value checkpoint_scalar_record(const Scalar& value);

template <>
json::value checkpoint_scalar_record<Rational>(const Rational& value) {
  return json::string(value.str());
}

template <>
json::value checkpoint_scalar_record<ComplexBall>(const ComplexBall& value) {
  return checkpoint_ball_record(value);
}

template <typename Scalar>
Scalar parse_checkpoint_scalar(const json::value& raw, const char* label);

template <>
Rational parse_checkpoint_scalar<Rational>(const json::value& raw,
                                            const char* label) {
  if (!raw.is_string())
    throw std::invalid_argument(std::string(label) +
                                " must be an exact rational string");
  return Rational(std::string(raw.as_string()));
}

template <>
ComplexBall parse_checkpoint_scalar<ComplexBall>(const json::value& raw,
                                                  const char* label) {
  return parse_checkpoint_ball(raw, label);
}

json::object checkpoint_exact_descriptor_record(
    const ExactScalarDescriptor& descriptor) {
  auto output = encode_exact_descriptor(descriptor);
  output.erase("has_specialization");
  output["specialization"] = descriptor.specialization.has_value()
      ? json::value(checkpoint_ball_record(*descriptor.specialization))
      : json::value(nullptr);
  return output;
}

ExactScalarDescriptor parse_checkpoint_exact_descriptor(
    const json::value& raw, const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(object,
      {"domain", "canonical", "symbols", "is_zero", "is_integer",
       "sign", "specialization"}, label);
  const auto domain = required_string(object, "domain");
  const auto canonical = required_string(object, "canonical");
  if (canonical.empty())
    throw std::invalid_argument(std::string(label) +
                                " canonical form is empty");
  std::vector<std::string> symbols;
  for (const auto& raw_symbol : as_array(object.at("symbols"), label)) {
    if (!raw_symbol.is_string() || raw_symbol.as_string().empty())
      throw std::invalid_argument(std::string(label) +
                                  " symbols must be nonempty strings");
    symbols.emplace_back(raw_symbol.as_string());
  }
  const auto zero = parse_truth_value(object.at("is_zero"), "is_zero");
  const auto integer = parse_truth_value(object.at("is_integer"),
                                         "is_integer");
  const auto sign = parse_exact_sign(object.at("sign"), "sign");
  std::optional<ComplexBall> specialization;
  if (!object.at("specialization").is_null())
    specialization = parse_checkpoint_ball(object.at("specialization"),
                                           "exact tag specialization");

  if (domain == "rational") {
    if (!symbols.empty())
      throw std::invalid_argument(std::string(label) +
                                  " rational tag cannot name symbols");
    auto result = ExactScalarDescriptor::rational(canonical);
    if (result.is_zero != zero || result.is_integer != integer ||
        result.sign != sign)
      throw std::invalid_argument(std::string(label) +
                                  " rational facts contradict its value");
    fmpq_t exact;
    fmpq_init(exact);
    const auto parse_status = fmpq_set_str(exact, canonical.c_str(), 10);
    if (parse_status == 0) fmpq_canonicalise(exact);
    const bool consistent = parse_status == 0 && specialization.has_value() &&
        acb_contains_fmpq(specialization->raw(), exact);
    fmpq_clear(exact);
    if (!consistent)
      throw std::invalid_argument(std::string(label) +
                                  " rational specialization is missing or inconsistent");
    result.specialization = std::move(specialization);
    return result;
  }
  if (domain == "symbolic-rational")
    return ExactScalarDescriptor::symbolic(
        canonical, std::move(symbols), zero, integer, sign,
        std::move(specialization));
  if (domain == "algebraic") {
    if (!symbols.empty())
      throw std::invalid_argument(std::string(label) +
                                  " algebraic tag cannot name regulator symbols");
    if (!specialization.has_value())
      throw std::invalid_argument(std::string(label) +
                                  " algebraic tag lost its specialization");
    return ExactScalarDescriptor::algebraic(
        canonical, zero, integer, sign, std::move(*specialization));
  }
  throw std::invalid_argument(std::string(label) +
                              " has an unsupported exact domain");
}

const char* checkpoint_error_guarantee_name(ErrorGuarantee guarantee) {
  if (guarantee == ErrorGuarantee::Certified) return "certified";
  if (guarantee == ErrorGuarantee::Advisory) return "advisory";
  return "none";
}

json::object checkpoint_error_envelope_record(
    const ErrorEnvelope& envelope) {
  json::array absolute;
  absolute.reserve(envelope.absolute.size());
  for (const auto& magnitude : envelope.absolute)
    absolute.emplace_back(magnitude.dump_exact());
  return json::object{
      {"frame", json::object{{"min", envelope.frame.min_power},
                              {"max", envelope.frame.complete_max}}},
      {"guarantee", checkpoint_error_guarantee_name(envelope.guarantee)},
      {"absolute_exact", std::move(absolute)},
      {"provenance", envelope.provenance}};
}

ErrorEnvelope parse_checkpoint_error_envelope(const json::value& raw) {
  const auto& object = as_object(raw, "checkpoint error envelope");
  require_exact_keys(object,
      {"frame", "guarantee", "absolute_exact", "provenance"},
      "checkpoint error envelope");
  const auto& frame = as_object(object.at("frame"),
                                "checkpoint error frame");
  require_exact_keys(frame, {"min", "max"}, "checkpoint error frame");
  ErrorEnvelope result;
  result.frame = {as_i32(frame.at("min"), "checkpoint error minimum"),
                  as_i32(frame.at("max"), "checkpoint error maximum")};
  (void)result.frame.width();
  const auto guarantee = required_string(object, "guarantee");
  result.guarantee = guarantee == "none" ? ErrorGuarantee::None
      : guarantee == "advisory" ? ErrorGuarantee::Advisory
      : guarantee == "certified" ? ErrorGuarantee::Certified
      : throw std::invalid_argument(
            "checkpoint error guarantee is unsupported");
  for (const auto& magnitude : as_array(
           object.at("absolute_exact"), "checkpoint error magnitudes")) {
    if (!magnitude.is_string())
      throw std::invalid_argument(
          "checkpoint error magnitudes must be exact dump strings");
    result.absolute.push_back(Magnitude::from_exact_dump(
        std::string(magnitude.as_string())));
  }
  if (!object.at("provenance").is_string())
    throw std::invalid_argument(
        "checkpoint error provenance must be a string");
  result.provenance = std::string(object.at("provenance").as_string());
  if (!result.absolute.empty() &&
      result.absolute.size() != result.frame.width())
    throw std::invalid_argument(
        "checkpoint error envelope width is inconsistent");
  if (result.absolute.empty() && result.guarantee != ErrorGuarantee::None)
    throw std::invalid_argument(
        "empty checkpoint error envelope cannot claim a guarantee");
  return result;
}

template <typename Scalar>
json::object checkpoint_epsilon_frame_record(
    const EpsilonFrame<Scalar>& frame) {
  json::array coefficients;
  coefficients.reserve(frame.coefficients().size());
  for (const auto& coefficient : frame.coefficients())
    coefficients.push_back(checkpoint_scalar_record<Scalar>(coefficient));
  return json::object{{"min", frame.min_power()},
                      {"max", frame.complete_max()},
                      {"coefficients", std::move(coefficients)}};
}

template <typename Scalar>
EpsilonFrame<Scalar> parse_checkpoint_epsilon_frame(
    const json::value& raw, const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(object, {"min", "max", "coefficients"}, label);
  EpsilonWindow window{as_i32(object.at("min"), label),
                       as_i32(object.at("max"), label)};
  const auto& raw_coefficients = as_array(object.at("coefficients"), label);
  if (raw_coefficients.size() != window.width())
    throw std::invalid_argument(std::string(label) +
                                " coefficient count is inconsistent");
  std::vector<Scalar> coefficients;
  coefficients.reserve(raw_coefficients.size());
  for (const auto& coefficient : raw_coefficients)
    coefficients.push_back(parse_checkpoint_scalar<Scalar>(coefficient,
                                                             label));
  return EpsilonFrame<Scalar>(window, std::move(coefficients));
}

template <typename Scalar>
json::array checkpoint_frame_vector_record(
    const FiniteLaurentVector<Scalar>& frames) {
  json::array output;
  output.reserve(frames.size());
  for (const auto& frame : frames)
    output.push_back(checkpoint_epsilon_frame_record(frame));
  return output;
}

template <typename Scalar>
FiniteLaurentVector<Scalar> parse_checkpoint_frame_vector(
    const json::value& raw, std::size_t expected_size, const char* label) {
  const auto& values = as_array(raw, label);
  if (values.size() != expected_size)
    throw std::invalid_argument(std::string(label) +
                                " dimension is inconsistent");
  FiniteLaurentVector<Scalar> output;
  output.reserve(values.size());
  for (const auto& value : values)
    output.push_back(parse_checkpoint_epsilon_frame<Scalar>(value, label));
  return output;
}

json::array checkpoint_exact_laurent_matrix_record(
    const ExactLaurentMatrix<Rational>& matrix) {
  json::array rows;
  rows.reserve(matrix.size());
  for (const auto& row : matrix) {
    json::array encoded_row;
    encoded_row.reserve(row.size());
    for (const auto& polynomial : row) {
      json::array terms;
      terms.reserve(polynomial.terms().size());
      for (const auto& [power, coefficient] : polynomial.terms())
        terms.push_back(json::object{{"power", power},
                                     {"coefficient", coefficient.str()}});
      encoded_row.push_back(std::move(terms));
    }
    rows.push_back(std::move(encoded_row));
  }
  return rows;
}

ExactLaurentMatrix<Rational> parse_checkpoint_exact_laurent_matrix(
    const json::value& raw, std::uint32_t dimension, const char* label) {
  const auto& rows = as_array(raw, label);
  if (rows.size() != dimension)
    throw std::invalid_argument(std::string(label) +
                                " row count differs from its dimension");
  ExactLaurentMatrix<Rational> matrix;
  matrix.reserve(rows.size());
  for (const auto& raw_row : rows) {
    const auto& row = as_array(raw_row, label);
    if (row.size() != dimension)
      throw std::invalid_argument(std::string(label) +
                                  " is not a square dimension-by-dimension matrix");
    std::vector<ExactLaurentPolynomial<Rational>> parsed_row;
    parsed_row.reserve(row.size());
    for (const auto& raw_entry : row) {
      ExactLaurentPolynomial<Rational> polynomial;
      std::optional<std::int32_t> previous_power;
      for (const auto& raw_term : as_array(raw_entry, label)) {
        const auto& term = as_object(raw_term, label);
        require_exact_keys(term, {"power", "coefficient"}, label);
        const auto power = as_i32(term.at("power"), label);
        if (previous_power.has_value() && power <= *previous_power)
          throw std::invalid_argument(std::string(label) +
                                      " terms are not in strict power order");
        previous_power = power;
        auto coefficient = parse_checkpoint_scalar<Rational>(
            term.at("coefficient"), label);
        if (coefficient.is_zero())
          throw std::invalid_argument(std::string(label) +
                                      " explicitly stores a structural-zero term");
        polynomial.add_term(power, std::move(coefficient));
      }
      parsed_row.push_back(std::move(polynomial));
    }
    matrix.push_back(std::move(parsed_row));
  }
  return matrix;
}

json::object checkpoint_saturation_diagnostics_record(
    const EpsilonLatticeSaturationDiagnostics<Rational>& diagnostics) {
  json::array valuations;
  for (const auto value : diagnostics.initial_column_valuations)
    valuations.push_back(value);
  json::array shifts;
  for (const auto value : diagnostics.initial_column_shifts)
    shifts.push_back(value);
  json::array actions;
  for (const auto& action : diagnostics.actions) {
    json::array relation;
    for (const auto& value : action.null_relation)
      relation.emplace_back(value.str());
    actions.push_back(json::object{
        {"leading_rank_before", action.leading_rank_before},
        {"target_column", action.target_column},
        {"null_relation", std::move(relation)}});
  }
  return json::object{
      {"initial_column_valuations", std::move(valuations)},
      {"initial_column_shifts", std::move(shifts)},
      {"normalized_determinant",
       checkpoint_epsilon_frame_record(diagnostics.normalized_determinant)},
      {"normalized_determinant_valuation",
       diagnostics.normalized_determinant_valuation},
      {"initial_leading_rank", diagnostics.initial_leading_rank},
      {"final_leading_rank", diagnostics.final_leading_rank},
      {"actions", std::move(actions)}};
}

EpsilonLatticeSaturationDiagnostics<Rational>
parse_checkpoint_saturation_diagnostics(const json::value& raw,
                                        std::uint32_t dimension) {
  const auto& object = as_object(raw, "checkpoint saturation diagnostics");
  require_exact_keys(
      object,
      {"initial_column_valuations", "initial_column_shifts",
       "normalized_determinant", "normalized_determinant_valuation",
       "initial_leading_rank", "final_leading_rank", "actions"},
      "checkpoint saturation diagnostics");
  std::vector<std::int32_t> valuations;
  std::vector<std::int32_t> shifts;
  for (const auto& value : as_array(object.at("initial_column_valuations"),
                                     "checkpoint saturation valuations"))
    valuations.push_back(
        as_i32(value, "checkpoint saturation valuation"));
  for (const auto& value : as_array(object.at("initial_column_shifts"),
                                     "checkpoint saturation shifts"))
    shifts.push_back(
        as_i32(value, "checkpoint saturation shift"));
  if (valuations.size() != dimension || shifts.size() != dimension)
    throw std::invalid_argument(
        "checkpoint saturation column metadata differs from its dimension");
  auto normalized_determinant = parse_checkpoint_epsilon_frame<Rational>(
      object.at("normalized_determinant"),
      "checkpoint normalized determinant");
  const auto normalized_determinant_valuation = as_i32(
      object.at("normalized_determinant_valuation"),
      "checkpoint normalized determinant valuation");
  const auto leading = finite_laurent_leading_power(
      normalized_determinant,
      "checkpoint normalized determinant valuation");
  if (!leading.has_value() ||
      *leading != normalized_determinant_valuation)
    throw std::invalid_argument(
        "checkpoint normalized determinant valuation is inconsistent");
  const auto initial_leading_rank = static_cast<std::size_t>(as_u64(
      object.at("initial_leading_rank"),
      "checkpoint initial leading rank"));
  const auto final_leading_rank = static_cast<std::size_t>(as_u64(
      object.at("final_leading_rank"), "checkpoint final leading rank"));
  if (initial_leading_rank > dimension || final_leading_rank != dimension ||
      initial_leading_rank > final_leading_rank)
    throw std::invalid_argument(
        "checkpoint saturation leading ranks are inconsistent");
  std::vector<EpsilonLatticeSaturationAction<Rational>> actions;
  for (const auto& raw_action : as_array(object.at("actions"),
                                          "checkpoint saturation actions")) {
    const auto& action = as_object(raw_action,
                                   "checkpoint saturation action");
    require_exact_keys(action,
        {"leading_rank_before", "target_column", "null_relation"},
        "checkpoint saturation action");
    EpsilonLatticeSaturationAction<Rational> parsed;
    parsed.leading_rank_before = static_cast<std::size_t>(as_u64(
        action.at("leading_rank_before"),
        "checkpoint saturation action rank"));
    parsed.target_column = static_cast<std::size_t>(as_u64(
        action.at("target_column"),
        "checkpoint saturation target column"));
    for (const auto& value : as_array(action.at("null_relation"),
                                      "checkpoint saturation null relation"))
      parsed.null_relation.push_back(parse_checkpoint_scalar<Rational>(
          value, "checkpoint saturation null relation coefficient"));
    if (parsed.leading_rank_before >= dimension ||
        parsed.target_column >= dimension ||
        parsed.null_relation.size() != dimension)
      throw std::invalid_argument(
          "checkpoint saturation action is outside its exact dimension");
    actions.push_back(std::move(parsed));
  }
  if (normalized_determinant_valuation < 0 ||
      actions.size() !=
          static_cast<std::size_t>(normalized_determinant_valuation))
    throw std::invalid_argument(
        "checkpoint saturation action count does not reproduce its determinant valuation");
  return {std::move(valuations), std::move(shifts),
          std::move(normalized_determinant),
          normalized_determinant_valuation, initial_leading_rank,
          final_leading_rank, std::move(actions)};
}

json::object checkpoint_epsilon_vector_record(const EpsilonVector& value) {
  json::array coefficients;
  coefficients.reserve(value.coefficients.size());
  for (const auto& coefficient : value.coefficients)
    coefficients.push_back(checkpoint_ball_record(coefficient));
  return json::object{
      {"epsilon", json::object{{"min", value.epsilon.min_power},
                                {"max", value.epsilon.complete_max}}},
      {"dimension", value.dimension},
      {"coefficients", std::move(coefficients)},
      {"error", checkpoint_error_envelope_record(value.error)}};
}

EpsilonVector parse_checkpoint_epsilon_vector(const json::value& raw,
                                               const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(object, {"epsilon", "dimension", "coefficients", "error"},
                     label);
  const auto& epsilon = as_object(object.at("epsilon"), label);
  require_exact_keys(epsilon, {"min", "max"}, label);
  EpsilonVector result;
  result.epsilon = {as_i32(epsilon.at("min"), label),
                    as_i32(epsilon.at("max"), label)};
  const auto width = result.epsilon.width();
  result.dimension = as_u32(object.at("dimension"), label);
  if (result.dimension == 0)
    throw std::invalid_argument(std::string(label) +
                                " has zero dimension");
  if (width > std::numeric_limits<std::size_t>::max() /
                  static_cast<std::size_t>(result.dimension))
    throw std::invalid_argument(std::string(label) +
                                " tensor size overflows size_t");
  const auto expected =
      width * static_cast<std::size_t>(result.dimension);
  const auto& coefficients = as_array(object.at("coefficients"), label);
  if (coefficients.size() != expected)
    throw std::invalid_argument(std::string(label) +
                                " coefficient count differs from its tensor");
  result.coefficients.reserve(expected);
  for (const auto& coefficient : coefficients)
    result.coefficients.push_back(parse_checkpoint_ball(coefficient, label));
  result.error = parse_checkpoint_error_envelope(object.at("error"));
  return result;
}

template <typename Scalar>
json::object checkpoint_local_analytic_metadata_record(
    const LocalSolution<Scalar>& solution) {
  json::array sectors;
  sectors.reserve(solution.sectors.size());
  for (const auto& sector : solution.sectors)
    sectors.push_back(json::object{
        {"a", checkpoint_exact_descriptor_record(sector.a)},
        {"b", checkpoint_exact_descriptor_record(sector.b)},
        {"log_power", sector.log_power}});
  json::array prescriptions;
  prescriptions.reserve(solution.prescriptions.size());
  for (const auto& prescription : solution.prescriptions)
    prescriptions.push_back(json::object{
        {"factor_exact", prescription.factor_exact},
        {"sign", prescription.sign},
        {"multiplicity", prescription.multiplicity},
        {"leading_coefficient_sign",
         prescription.leading_coefficient_sign}});
  return json::object{
      {"schema", "diffexp2-exact-local-analytic-metadata-v2"},
      {"chart", json::object{
          {"center_exact", solution.chart.center_exact},
          {"scale_exact", solution.chart.scale_exact},
          {"radius_exact_ball", checkpoint_ball_record(solution.chart.radius)},
          {"infinite_radius", solution.chart.infinite_radius}}},
      {"sectors", std::move(sectors)},
      {"prescriptions", std::move(prescriptions)}};
}

template <typename Scalar>
json::object checkpoint_local_solution_record(
    const LocalSolution<Scalar>& solution) {
  json::array sectors;
  sectors.reserve(solution.sectors.size());
  for (const auto& sector : solution.sectors) {
    json::array coefficients;
    coefficients.reserve(sector.coefficients.size());
    for (const auto& coefficient : sector.coefficients)
      coefficients.push_back(checkpoint_scalar_record<Scalar>(coefficient));
    sectors.push_back(json::object{
        {"a", checkpoint_exact_descriptor_record(sector.a)},
        {"b", checkpoint_exact_descriptor_record(sector.b)},
        {"log_power", sector.log_power},
        {"coefficients", std::move(coefficients)}});
  }
  json::array prescriptions;
  prescriptions.reserve(solution.prescriptions.size());
  for (const auto& prescription : solution.prescriptions)
    prescriptions.push_back(json::object{
        {"factor_exact", prescription.factor_exact},
        {"sign", prescription.sign},
        {"multiplicity", prescription.multiplicity},
        {"leading_coefficient_sign",
         prescription.leading_coefficient_sign}});
  return json::object{
      {"chart", json::object{
          {"center_exact", solution.chart.center_exact},
          {"scale_exact", solution.chart.scale_exact},
          {"radius_exact_ball", checkpoint_ball_record(solution.chart.radius)},
          {"infinite_radius", solution.chart.infinite_radius}}},
      {"epsilon", json::object{{"min", solution.epsilon.min_power},
                                {"max", solution.epsilon.complete_max}}},
      {"taylor_complete_max", solution.taylor_complete_max},
      {"dimension", solution.dimension},
      {"sectors", std::move(sectors)},
      {"prescriptions", std::move(prescriptions)},
      {"error", checkpoint_error_envelope_record(solution.error)},
      {"checkpoint_identity", solution.checkpoint_identity}};
}

json::object checkpoint_regular_tail_model_record(
    const RegularTaylorTailModel& model) {
  json::array q_coefficients;
  q_coefficients.reserve(model.q_coefficients.size());
  for (const auto& coefficient : model.q_coefficients)
    q_coefficients.push_back(checkpoint_ball_record(coefficient));
  json::array n_coefficients;
  n_coefficients.reserve(model.n_coefficients.size());
  for (const auto& matrix : model.n_coefficients) {
    json::array encoded_matrix;
    encoded_matrix.reserve(matrix.size());
    for (const auto& coefficient : matrix)
      encoded_matrix.push_back(checkpoint_ball_record(coefficient));
    n_coefficients.push_back(std::move(encoded_matrix));
  }
  json::array n_row_sum_upper;
  n_row_sum_upper.reserve(model.n_row_sum_upper.size());
  for (const auto& magnitude : model.n_row_sum_upper)
    n_row_sum_upper.emplace_back(magnitude.dump_exact());
  json::array initial_row_upper;
  initial_row_upper.reserve(model.initial_row_upper.size());
  for (const auto& magnitude : model.initial_row_upper)
    initial_row_upper.emplace_back(magnitude.dump_exact());
  json::array prescriptions;
  prescriptions.reserve(model.prescriptions.size());
  for (const auto& prescription : model.prescriptions)
    prescriptions.push_back(json::object{
        {"factor_exact", prescription.factor_exact},
        {"sign", prescription.sign},
        {"multiplicity", prescription.multiplicity},
        {"leading_coefficient_sign",
         prescription.leading_coefficient_sign}});
  return json::object{
      {"schema", "diffexp2-regular-taylor-tail-model-v1"},
      {"epsilon", json::object{{"min", model.epsilon.min_power},
                                {"max", model.epsilon.complete_max}}},
      {"dimension", model.dimension},
      {"taylor_complete_max", model.taylor_complete_max},
      {"q_coefficients", std::move(q_coefficients)},
      {"n_coefficients", std::move(n_coefficients)},
      {"n_row_sum_upper_exact", std::move(n_row_sum_upper)},
      {"initial_row_upper_exact", std::move(initial_row_upper)},
      {"chart", json::object{
           {"center_exact", model.chart.center_exact},
           {"scale_exact", model.chart.scale_exact},
           {"radius_exact_ball", checkpoint_ball_record(model.chart.radius)},
           {"infinite_radius", model.chart.infinite_radius}}},
      {"prescriptions", std::move(prescriptions)},
      {"operator_identity", model.operator_identity},
      {"local_checkpoint_identity", model.local_checkpoint_identity},
      {"provenance", model.provenance}};
}

RegularTaylorTailModel parse_checkpoint_regular_tail_model(
    const json::value& raw) {
  const auto& object = as_object(raw, "checkpoint regular tail model");
  require_exact_keys(
      object,
      {"schema", "epsilon", "dimension", "taylor_complete_max",
       "q_coefficients", "n_coefficients", "n_row_sum_upper_exact",
       "initial_row_upper_exact", "chart", "prescriptions",
       "operator_identity", "local_checkpoint_identity", "provenance"},
      "checkpoint regular tail model");
  if (required_string(object, "schema") !=
      "diffexp2-regular-taylor-tail-model-v1")
    throw std::invalid_argument(
        "unsupported checkpoint regular tail-model schema");
  RegularTaylorTailModel model;
  const auto& epsilon = as_object(
      object.at("epsilon"), "checkpoint tail epsilon window");
  require_exact_keys(epsilon, {"min", "max"},
                     "checkpoint tail epsilon window");
  model.epsilon = {
      as_i32(epsilon.at("min"), "checkpoint tail epsilon minimum"),
      as_i32(epsilon.at("max"), "checkpoint tail epsilon maximum")};
  (void)model.epsilon.width();
  model.dimension = as_u32(
      object.at("dimension"), "checkpoint tail dimension");
  if (model.dimension == 0)
    throw std::invalid_argument("checkpoint tail dimension is zero");
  model.taylor_complete_max = as_u32(
      object.at("taylor_complete_max"),
      "checkpoint tail Taylor complete maximum");
  for (const auto& coefficient : as_array(
           object.at("q_coefficients"),
           "checkpoint tail q coefficients"))
    model.q_coefficients.push_back(parse_checkpoint_ball(
        coefficient, "checkpoint tail q coefficient"));
  for (const auto& raw_matrix : as_array(
           object.at("n_coefficients"),
           "checkpoint tail N coefficients")) {
    std::vector<ComplexBall> matrix;
    for (const auto& coefficient : as_array(
             raw_matrix, "checkpoint tail N matrix"))
      matrix.push_back(parse_checkpoint_ball(
          coefficient, "checkpoint tail N coefficient"));
    model.n_coefficients.push_back(std::move(matrix));
  }
  for (const auto& magnitude : as_array(
           object.at("n_row_sum_upper_exact"),
           "checkpoint tail N norms"))
    model.n_row_sum_upper.push_back(parse_checkpoint_magnitude(
        magnitude, "checkpoint tail N norm"));
  for (const auto& magnitude : as_array(
           object.at("initial_row_upper_exact"),
           "checkpoint tail initial magnitudes"))
    model.initial_row_upper.push_back(parse_checkpoint_magnitude(
        magnitude, "checkpoint tail initial magnitude"));
  const auto& chart = as_object(
      object.at("chart"), "checkpoint tail chart");
  require_exact_keys(
      chart,
      {"center_exact", "scale_exact", "radius_exact_ball",
       "infinite_radius"},
      "checkpoint tail chart");
  model.chart.center_exact = required_string(chart, "center_exact");
  model.chart.scale_exact = required_string(chart, "scale_exact");
  model.chart.radius = parse_checkpoint_ball(
      chart.at("radius_exact_ball"), "checkpoint tail chart radius");
  if (!chart.at("infinite_radius").is_bool())
    throw std::invalid_argument(
        "checkpoint tail infinite-radius flag must be Boolean");
  model.chart.infinite_radius = chart.at("infinite_radius").as_bool();
  for (const auto& raw_prescription : as_array(
           object.at("prescriptions"), "checkpoint tail prescriptions")) {
    const auto& prescription = as_object(
        raw_prescription, "checkpoint tail prescription");
    require_exact_keys(
        prescription,
        {"factor_exact", "sign", "multiplicity",
         "leading_coefficient_sign"},
        "checkpoint tail prescription");
    const auto sign = as_i32(
        prescription.at("sign"), "checkpoint tail prescription sign");
    const auto multiplicity = as_u32(
        prescription.at("multiplicity"),
        "checkpoint tail prescription multiplicity");
    const auto leading = as_i32(
        prescription.at("leading_coefficient_sign"),
        "checkpoint tail leading-coefficient sign");
    if ((sign != -1 && sign != 1) || multiplicity == 0 ||
        (leading != -1 && leading != 1))
      throw std::invalid_argument(
          "checkpoint tail prescription is malformed");
    model.prescriptions.push_back(Prescription{
        required_string(prescription, "factor_exact"), sign,
        multiplicity, leading});
  }
  model.operator_identity = required_string(object, "operator_identity");
  model.local_checkpoint_identity = required_string(
      object, "local_checkpoint_identity");
  model.provenance = required_string(object, "provenance");
  return model;
}

template <typename Scalar>
json::array checkpoint_pseudo_hits_record(
    const std::vector<PseudoHit<Scalar>>& hits) {
  json::array output;
  output.reserve(hits.size());
  for (const auto& hit : hits) {
    json::array columns;
    columns.reserve(hit.columns.size());
    for (const auto column : hit.columns) columns.emplace_back(column);
    json::array frames;
    frames.reserve(hit.gamma_frames.size());
    for (const auto& frame : hit.gamma_frames) {
      json::array coefficients;
      coefficients.reserve(frame.size());
      for (const auto& coefficient : frame)
        coefficients.push_back(checkpoint_scalar_record<Scalar>(coefficient));
      frames.push_back(std::move(coefficients));
    }
    json::array validity;
    validity.reserve(hit.gamma_validity.size());
    for (const auto value : hit.gamma_validity)
      validity.push_back(encode_validity(value));
    output.push_back(json::object{
        {"n", hit.n}, {"columns", std::move(columns)},
        {"delta_b", checkpoint_scalar_record<Scalar>(hit.delta_b)},
        {"gamma_frames", std::move(frames)},
        {"gamma_validity", std::move(validity)}});
  }
  return output;
}

template <typename Scalar>
struct NativeLocalRun {
  LocalSolution<Scalar> solution;
  std::vector<PseudoHit<Scalar>> pseudo_hits;
  NativeLocalDiagnostics diagnostics;
  RegularTaylorTailModelResult tail_model = unavailable_tail_model(
      "tail model was not prepared for this native local run");
};

json::value canonical_json_value(const json::value& value);
void require_exact_keys(const json::object& object,
                        std::initializer_list<std::string_view> required,
                        const char* label);
std::uint64_t scoped_handle_id(const std::string& handle,
                               std::string_view prefix,
                               const char* label);

struct SCCColumnProvenance {
  std::string scc_handle;
  std::string scc_exact_identity;
  std::uint32_t seed_block = 0;
  std::uint32_t basis_index = 0;
  std::string exact_column_identity;

  json::object encode() const {
    return json::object{{"scc", scc_handle},
                        {"scc_exact_identity", scc_exact_identity},
                        {"seed_block", seed_block},
                        {"basis_index", basis_index},
                        {"exact_column_identity", exact_column_identity}};
  }
};

class StoredLocalBase {
 public:
  StoredLocalBase(std::string handle, std::string source_chart,
                  std::string source_operator_identity,
                  double create_parse_ms, double create_kernel_ms,
                  std::optional<SCCColumnProvenance> column_provenance =
                      std::nullopt)
      : handle_(std::move(handle)), source_chart_(std::move(source_chart)),
        source_operator_identity_(std::move(source_operator_identity)),
        create_parse_ms_(create_parse_ms),
        create_kernel_ms_(create_kernel_ms),
        column_provenance_(std::move(column_provenance)) {}
  virtual ~StoredLocalBase() = default;

  virtual json::object evaluate(const json::object& request,
                                int output_digits) = 0;
  virtual json::object certify_residual(const json::object& request,
                                        int output_digits) = 0;
  virtual EndpointLimitResult endpoint_limit(
      const EndpointLimitOptions& options) = 0;
  virtual StoredLineIntegral integrate_planned_line(
      const ExactAffineChart& chart,
      const std::vector<Prescription>& prescriptions,
      const Rational& local_begin, const Rational& local_end,
      const EpsilonWindow& delivered_epsilon,
      std::optional<std::int32_t> exact_rim,
      bool certify_tail) = 0;
  virtual json::object endpoint_metadata() const = 0;
  virtual json::object exact_analytic_metadata() const = 0;
  virtual void require_exact_plan_binding(
      const ExactAffineChart& chart,
      const std::vector<Prescription>& prescriptions,
      const std::string& label) const = 0;
  virtual const std::string& checkpoint_identity() const = 0;
  virtual const char* scalar_domain() const = 0;
  virtual json::object summary() const = 0;
  virtual json::object stats_json() const = 0;
  virtual StoredLocalStats stats() const = 0;
  virtual json::object checkpoint_record() const = 0;
  virtual const std::optional<json::object>& retained_derivation() const = 0;
  virtual std::shared_ptr<void> retained_derivation_owner() const = 0;

  const std::string& handle() const { return handle_; }
  const std::string& source_chart() const { return source_chart_; }
  const std::string& source_operator_identity() const {
    return source_operator_identity_;
  }
  const std::optional<SCCColumnProvenance>& column_provenance() const {
    return column_provenance_;
  }

 protected:
  std::string handle_;
  std::string source_chart_;
  std::string source_operator_identity_;
  double create_parse_ms_ = 0.0;
  double create_kernel_ms_ = 0.0;
  std::optional<SCCColumnProvenance> column_provenance_;
};

template <typename Scalar>
class StoredLocal final : public StoredLocalBase {
 public:
  StoredLocal(std::string handle, std::string source_chart,
              std::string source_operator_identity,
              LocalSolution<Scalar>&& solution, slong precision_bits,
              std::vector<PseudoHit<Scalar>>&& pseudo_hits,
              NativeLocalDiagnostics diagnostics,
              std::optional<SCCColumnProvenance> column_provenance =
                  std::nullopt,
              std::optional<json::object> retained_derivation =
                  std::nullopt,
              std::shared_ptr<void> retained_owner = nullptr,
              RegularTaylorTailModelResult tail_model =
                  unavailable_tail_model(
                      "tail model is unavailable for this retained local"),
              std::optional<TailModelCheckpointMarker>
                  tail_checkpoint_marker = std::nullopt,
              bool serialize_tail_checkpoint_fields = true,
              bool serialize_derivation_checkpoint_fields = true)
      : StoredLocalBase(std::move(handle), std::move(source_chart),
                        std::move(source_operator_identity),
                        diagnostics.parse_ms, diagnostics.kernel_ms,
                        std::move(column_provenance)),
        solution_(std::move(solution)), precision_bits_(precision_bits),
        pseudo_hits_(std::move(pseudo_hits)), top_valid_(diagnostics.top_valid),
        retained_derivation_(std::move(retained_derivation)),
        retained_owner_(std::move(retained_owner)),
        tail_model_(std::move(tail_model)),
        tail_checkpoint_marker_(std::move(tail_checkpoint_marker)),
        serialize_tail_checkpoint_fields_(serialize_tail_checkpoint_fields),
        serialize_derivation_checkpoint_fields_(
            serialize_derivation_checkpoint_fields) {
    validate_local_solution(solution_, false);
  }

  json::object evaluate(const json::object& request,
                        int output_digits) override {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      throw std::domain_error(
          "local.evaluate rejects unresolved symbolic coefficients; solve a "
          "numerically specialized chart before native evaluation");
    } else {
      AcbPrecisionLease lease(precision_bits_);
      ComplexBall::set_precision(precision_bits_);
      const auto point = parse_local_evaluation_point(request);
      auto options = parse_local_evaluation_options(request, true);
      const auto tail_witness = parse_certified_tail_witness(request);
      if (tail_witness.has_value())
        options.compute_tail_estimate = false;
      const auto started = std::chrono::steady_clock::now();
      LocalEvaluation result;
      std::optional<RegularTaylorPointTailCertificate> tail_certificate;
      TailMajorantStatus requested_tail_status = tail_model_.status;
      std::string requested_tail_detail = tail_model_.detail;
      if (tail_witness.has_value() && tail_model_.model.has_value()) {
        auto certified = evaluate_local_solution_with_certified_tail(
            solution_, *tail_model_.model, point, *tail_witness, options);
        result = std::move(certified.evaluation);
        requested_tail_status = certified.tail.status;
        requested_tail_detail = certified.tail.detail;
        tail_certificate = std::move(certified.tail);
      } else {
        result = evaluate_local_solution(solution_, point, options);
      }
      const auto elapsed = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      evaluations_.fetch_add(1);
      if (tail_witness.has_value()) {
        tail_certificate_requests_.fetch_add(1);
        switch (requested_tail_status) {
          case TailMajorantStatus::Certified:
            tail_certificate_certified_.fetch_add(1);
            break;
          case TailMajorantStatus::Inconclusive:
            tail_certificate_inconclusive_.fetch_add(1);
            break;
          case TailMajorantStatus::Unsupported:
            tail_certificate_unsupported_.fetch_add(1);
            break;
        }
      }
      {
        std::lock_guard<std::mutex> lock(stats_mutex_);
        evaluate_ms_ += elapsed;
      }
      json::object response{
          {"point_exact", point.exact_coordinate},
          {"imaginary_sign", result.imaginary_sign.has_value()
               ? json::value(*result.imaginary_sign) : json::value(nullptr)},
          {"arithmetic_enclosed", result.arithmetic_enclosed},
          {"elapsed_ms", elapsed},
          {"value", encode_epsilon_vector(result.value, output_digits)},
          {"theta", encode_epsilon_vector(result.theta_value, output_digits)}};
      if (tail_witness.has_value()) {
        if (tail_certificate.has_value()) {
          response["tail_certificate"] = encode_point_tail_certificate(
              *tail_certificate, *tail_witness);
        } else {
          response["tail_certificate"] = json::object{
              {"capability", kRegularTailMajorantCapability},
              {"requested", true},
              {"status", tail_majorant_status_name(requested_tail_status)},
              {"model_status", tail_majorant_status_name(tail_model_.status)},
              {"witness_radius_exact", *tail_witness},
              {"detail", requested_tail_detail}};
        }
      }
      return response;
    }
  }

  json::object certify_residual(const json::object& request,
                                int output_digits) override {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      throw std::domain_error(
          "local.certify_residual rejects unresolved symbolic coefficients; "
          "solve a numerically specialized chart first");
    } else {
      AcbPrecisionLease lease(precision_bits_);
      ComplexBall::set_precision(precision_bits_);
      const auto checkpoint_identity = required_string(
          request, "checkpoint_identity");
      const auto operator_identity = required_string(
          request, "operator_identity");
      if (checkpoint_identity.empty() || operator_identity.empty())
        throw std::invalid_argument(
            "native residual identities must be nonempty");
      const auto point = parse_local_evaluation_point(request);
      const auto options = parse_local_evaluation_options(request, false);
      if (options.compute_tail_estimate)
        throw std::invalid_argument(
            "stored-truncation residual certification rejects advisory tail estimates");
      const auto theta_operator = parse_epsilon_matrix(
          request.at("theta_operator"), "theta operator");
      std::optional<EpsilonVector> source;
      std::optional<std::string> source_identity;
      if (const auto* raw_source = request.if_contains("source");
          raw_source != nullptr && !raw_source->is_null()) {
        source = parse_epsilon_vector(*raw_source, "residual source");
        source_identity = required_string(request, "source_identity");
        if (source_identity->empty())
          throw std::invalid_argument(
              "native residual source identity must be nonempty");
      } else if (const auto* raw_identity =
                     request.if_contains("source_identity");
                 raw_identity != nullptr && !raw_identity->is_null()) {
        throw std::invalid_argument(
            "source_identity requires an explicit residual source");
      }
      const auto tolerance_text = required_string(
          request, "relative_tolerance");
      const auto tolerance = Magnitude::decimal(tolerance_text);
      const auto scope_text = required_string(request, "scope");
      const auto scope = scope_text == "stored_truncation"
          ? ResidualScope::StoredTruncation
          : scope_text == "full_local_solution"
              ? ResidualScope::FullLocalSolution
              : throw std::invalid_argument(
                    "residual scope must be stored_truncation or full_local_solution");
      const auto include_residual =
          request.if_contains("include_residual") != nullptr &&
          request.at("include_residual").as_bool();

      const auto started = std::chrono::steady_clock::now();
      const auto evaluation = evaluate_local_solution(solution_, point, options);
      auto certificate = certify_theta_residual(
          evaluation, theta_operator, source, tolerance, scope);
      const auto elapsed = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      residual_certifications_.fetch_add(1);
      {
        std::lock_guard<std::mutex> lock(stats_mutex_);
        residual_certify_ms_ += elapsed;
      }

      const char* verdict = certificate.verdict == ResidualVerdict::Pass
          ? "pass"
          : certificate.verdict == ResidualVerdict::Fail
              ? "fail" : "inconclusive";
      json::object result{
          {"point_exact", point.exact_coordinate},
          {"imaginary_sign", evaluation.imaginary_sign.has_value()
               ? json::value(*evaluation.imaginary_sign)
               : json::value(nullptr)},
          {"arithmetic_enclosed", evaluation.arithmetic_enclosed},
          {"scope", scope_text}, {"verdict", verdict},
          {"detail", certificate.detail},
          {"relative_tolerance", tolerance_text},
          {"checkpoint_identity", checkpoint_identity},
          {"operator_identity", operator_identity},
          {"source_identity", source_identity.has_value()
               ? json::value(*source_identity) : json::value(nullptr)},
          {"epsilon_min", certificate.residual.epsilon.min_power},
          {"epsilon_max", certificate.residual.epsilon.complete_max},
          {"dimension", certificate.residual.dimension},
          {"residual_upper_approx",
           encode_magnitude_diagnostics(certificate.residual_upper)},
          {"scale_lower_approx",
           encode_magnitude_diagnostics(certificate.scale_lower)},
          {"relative_upper_approx",
           encode_magnitude_diagnostics(certificate.relative_upper)},
          {"bound_encoding", "approximate-double-diagnostics"},
          {"json_coefficients", include_residual
               ? certificate.residual.coefficients.size() : 0},
          {"elapsed_ms", elapsed}};
      if (include_residual)
        result["residual"] = encode_epsilon_vector(
            certificate.residual, output_digits);
      return result;
    }
  }

  EndpointLimitResult endpoint_limit(
      const EndpointLimitOptions& options) override {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      throw std::domain_error(
          "local.endpoint_limit rejects unresolved symbolic coefficients; "
          "solve an explicitly specialized chart before endpoint evaluation");
    } else {
      AcbPrecisionLease lease(precision_bits_);
      ComplexBall::set_precision(precision_bits_);
      const auto started = std::chrono::steady_clock::now();
      auto result = endpoint_sector_limit(solution_, options);
      const auto elapsed = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      endpoint_limits_.fetch_add(1);
      {
        std::lock_guard<std::mutex> lock(stats_mutex_);
        endpoint_limit_ms_ += elapsed;
      }
      return result;
    }
  }

  StoredLineIntegral integrate_planned_line(
      const ExactAffineChart& chart,
      const std::vector<Prescription>& prescriptions,
      const Rational& local_begin, const Rational& local_end,
      const EpsilonWindow& delivered_epsilon,
      std::optional<std::int32_t> exact_rim,
      bool certify_tail) override {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      throw std::domain_error(
          "native planned line integration rejects unresolved symbolic "
          "coefficients; specialize the retained local first");
    } else {
      AcbPrecisionLease lease(precision_bits_);
      ComplexBall::set_precision(precision_bits_);
      if (!(Rational(solution_.chart.center_exact) == chart.center) ||
          !(Rational(solution_.chart.scale_exact) == chart.scale) ||
          solution_.chart.infinite_radius ||
          !acb_equal(solution_.chart.radius.raw(),
                     ComplexBall::from_strings(chart.radius.str()).raw()))
        throw std::invalid_argument(
            "retained local chart geometry differs from its native tile plan");
      const auto same_prescription = [](const Prescription& left,
                                        const Prescription& right) {
        return left.factor_exact == right.factor_exact &&
               left.sign == right.sign &&
               left.multiplicity == right.multiplicity &&
               left.leading_coefficient_sign ==
                   right.leading_coefficient_sign;
      };
      if (solution_.prescriptions.size() != prescriptions.size() ||
          !std::equal(solution_.prescriptions.begin(),
                      solution_.prescriptions.end(), prescriptions.begin(),
                      same_prescription))
        throw std::invalid_argument(
            "retained local prescriptions differ from its native tile plan");

      StoredLineIntegrationOptions options;
      options.delivered_epsilon = delivered_epsilon;
      options.imaginary_sign = exact_rim;
      if (local_begin == local_end)
        throw std::invalid_argument(
            "native planned line tile has zero local length");
      const bool reverse_local_orientation = local_end < local_begin;
      const auto& primitive_begin =
          reverse_local_orientation ? local_end : local_begin;
      const auto& primitive_end =
          reverse_local_orientation ? local_begin : local_end;
      const auto started = std::chrono::steady_clock::now();
      const auto begin_point =
          RealEvaluationPoint::rational(primitive_begin.str());
      const auto end_point =
          RealEvaluationPoint::rational(primitive_end.str());
      StoredLineIntegral result;
      if (certify_tail && tail_model_.model.has_value()) {
        const auto begin_modulus = primitive_begin.sign() < 0
            ? -primitive_begin : primitive_begin;
        const auto end_modulus = primitive_end.sign() < 0
            ? -primitive_end : primitive_end;
        const auto outer = begin_modulus < end_modulus
            ? end_modulus : begin_modulus;
        if (!(outer < chart.radius))
          throw std::invalid_argument(
              "planned line endpoint is not strictly inside its exact chart radius");
        const auto witness =
            (outer + chart.radius) / Rational(2);
        auto certified = integrate_regular_local_line_with_certified_tail(
            solution_, *tail_model_.model, begin_point, end_point,
            options, witness.str());
        result = std::move(certified.integral);
      } else {
        result = integrate_stored_local_line(
            solution_, begin_point, end_point, options);
        if (certify_tail) {
          result.diagnostics.tail_certificate_requested = true;
          result.diagnostics.tail_certificate_status =
              tail_majorant_status_name(tail_model_.status);
          result.value.error.provenance =
              std::string(tail_majorant_status_name(tail_model_.status)) +
              ": " + tail_model_.detail +
              "; returned value remains stored Taylor truncation only";
          result.diagnostics.detail = result.value.error.provenance;
        }
      }
      // integrate_stored_local_line is deliberately a local-coordinate
      // primitive.  A retained physical tile has dx = scale dt, so apply the
      // exact affine Jacobian before publishing the physical line result.
      const auto oriented_jacobian = reverse_local_orientation
          ? -chart.scale : chart.scale;
      const auto jacobian =
          ComplexBall::from_strings(oriented_jacobian.str());
      for (auto& coefficient : result.value.coefficients)
        coefficient *= jacobian;
      if (!result.value.error.empty()) {
        const auto jacobian_upper = Magnitude::upper_abs(jacobian);
        for (auto& bound : result.value.error.absolute)
          bound = bound * jacobian_upper;
        result.value.error.provenance +=
            "; physical_jacobian_exact=" + oriented_jacobian.str();
      }
      const auto elapsed = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      line_integrations_.fetch_add(1);
      {
        std::lock_guard<std::mutex> lock(stats_mutex_);
        line_integration_ms_ += elapsed;
      }
      return result;
    }
  }

  json::object endpoint_metadata() const override { return metadata_json(); }

  json::object exact_analytic_metadata() const override {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      throw std::domain_error(
          "exact checkpoint analytic metadata does not serialize symbolic coefficients");
    } else {
      return checkpoint_local_analytic_metadata_record(solution_);
    }
  }

  void require_exact_plan_binding(
      const ExactAffineChart& chart,
      const std::vector<Prescription>& prescriptions,
      const std::string& label) const override {
    AcbPrecisionLease lease(precision_bits_);
    ComplexBall::set_precision(precision_bits_);
    if (!(Rational(solution_.chart.center_exact) == chart.center) ||
        !(Rational(solution_.chart.scale_exact) == chart.scale) ||
        solution_.chart.infinite_radius ||
        !acb_equal(solution_.chart.radius.raw(),
                   ComplexBall::from_strings(chart.radius.str()).raw()))
      throw std::invalid_argument(
          label +
          " retained local geometry differs from its exact tile-plan chart");
    const auto same_prescription = [](const Prescription& left,
                                      const Prescription& right) {
      return left.factor_exact == right.factor_exact &&
             left.sign == right.sign &&
             left.multiplicity == right.multiplicity &&
             left.leading_coefficient_sign ==
                 right.leading_coefficient_sign;
    };
    if (solution_.prescriptions.size() != prescriptions.size() ||
        !std::equal(solution_.prescriptions.begin(),
                    solution_.prescriptions.end(), prescriptions.begin(),
                    same_prescription))
      throw std::invalid_argument(
          label +
          " retained local prescriptions differ from its exact tile-plan chart");
  }

  const std::string& checkpoint_identity() const override {
    return solution_.checkpoint_identity;
  }

  const char* scalar_domain() const override {
    if constexpr (std::is_same_v<Scalar, Rational>) return "rational";
    if constexpr (std::is_same_v<Scalar, ComplexBall>) return "acb";
    return "symbolic";
  }

  json::object summary() const override {
    json::object result{
        {"local", handle_}, {"chart", source_chart_},
        {"source_operator_identity", source_operator_identity_},
        {"dimension", solution_.dimension},
        {"epsilon_min", solution_.epsilon.min_power},
        {"epsilon_max", solution_.epsilon.complete_max},
        {"taylor_complete_max", solution_.taylor_complete_max},
        {"sectors", solution_.sectors.size()},
        {"coefficient_count", coefficient_count()},
        {"pseudo_hit_count", pseudo_hits_.size()},
        {"top_valid", encode_validity(top_valid_)},
        {"checkpoint_identity", solution_.checkpoint_identity},
        {"tail_majorant", encode_tail_model_status(tail_model_)},
        {"metadata", metadata_json()},
        {"create_parse_ms", create_parse_ms_},
        {"create_kernel_ms", create_kernel_ms_}};
    if (column_provenance_.has_value())
      result["column_provenance"] = column_provenance_->encode();
    if (retained_derivation_.has_value()) {
      result["retained_derivation"] = *retained_derivation_;
      result["strong_derivation_ownership"] =
          retained_owner_ != nullptr;
    }
    return result;
  }

  json::object stats_json() const override {
    auto out = summary();
    const auto current = stats();
    out["evaluations"] = current.evaluations;
    out["residual_certifications"] = current.residual_certifications;
    out["endpoint_limits"] = current.endpoint_limits;
    out["line_integrations"] = current.line_integrations;
    out["evaluate_ms"] = current.evaluate_ms;
    out["residual_certify_ms"] = current.residual_certify_ms;
    out["endpoint_limit_ms"] = current.endpoint_limit_ms;
    out["line_integration_ms"] = current.line_integration_ms;
    out["tail_certificate_requests"] = tail_certificate_requests_.load();
    out["tail_certificate_certified"] =
        tail_certificate_certified_.load();
    out["tail_certificate_inconclusive"] =
        tail_certificate_inconclusive_.load();
    out["tail_certificate_unsupported"] =
        tail_certificate_unsupported_.load();
    return out;
  }

  StoredLocalStats stats() const override {
    std::lock_guard<std::mutex> lock(stats_mutex_);
    return {evaluations_.load(), residual_certifications_.load(),
            endpoint_limits_.load(), line_integrations_.load(), evaluate_ms_,
            residual_certify_ms_, endpoint_limit_ms_, line_integration_ms_,
            create_parse_ms_, create_kernel_ms_, coefficient_count(),
            tail_certificate_requests_.load(),
            tail_certificate_certified_.load(),
            tail_certificate_inconclusive_.load(),
            tail_certificate_unsupported_.load()};
  }

  json::object checkpoint_record() const override {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      throw std::domain_error(
          "native checkpoint does not serialize symbolic-coefficient local state");
    } else {
      AcbPrecisionLease lease(precision_bits_);
      ComplexBall::set_precision(precision_bits_);
      const auto current = stats();
      if (retained_derivation_.has_value()) {
        const auto schema = required_string(
            *retained_derivation_, "schema");
        if (schema !=
                "diffexp2-retained-plan-match-local-materialization-v1" &&
            schema !=
                "diffexp2-retained-rational-row-local-application-v1")
          throw std::domain_error(
              "native checkpoint does not serialize this retained local derivation kind");
        if (retained_owner_ == nullptr)
          throw std::logic_error(
              "derived local lost its strong derivation owner before checkpointing");
      } else if (retained_owner_ != nullptr) {
        throw std::logic_error(
            "primitive local unexpectedly retains a derivation owner before checkpointing");
      }
      json::value owner_lineage = nullptr;
      if (retained_derivation_.has_value()) {
        const auto& derivation = *retained_derivation_;
        if (required_string(derivation, "schema") ==
            "diffexp2-retained-plan-match-local-materialization-v1") {
          owner_lineage = json::object{
              {"match", derivation.at("source_match")},
              {"match_checkpoint_identity",
               derivation.at("source_match_checkpoint_identity")},
              {"match_provenance_identity",
               derivation.at("source_match_provenance_identity")},
              {"planned_hop_provenance_identity",
               derivation.at("planned_hop_provenance_identity")},
              {"derivation_provenance_identity",
               derivation.at("provenance_identity")}};
        } else {
          owner_lineage = rational_row_owner_lineage();
        }
      }
      json::object runtime{
          {"evaluations", current.evaluations},
          {"residual_certifications", current.residual_certifications},
          {"endpoint_limits", current.endpoint_limits},
          {"line_integrations", current.line_integrations},
          {"evaluate_ms", current.evaluate_ms},
          {"residual_certify_ms", current.residual_certify_ms},
          {"endpoint_limit_ms", current.endpoint_limit_ms},
          {"line_integration_ms", current.line_integration_ms},
          {"coefficient_count", current.coefficient_count}};
      if (serialize_tail_checkpoint_fields_) {
        runtime["tail_certificate_requests"] =
            current.tail_certificate_requests;
        runtime["tail_certificate_certified"] =
            current.tail_certificate_certified;
        runtime["tail_certificate_inconclusive"] =
            current.tail_certificate_inconclusive;
        runtime["tail_certificate_unsupported"] =
            current.tail_certificate_unsupported;
      }
      json::object record{
        {"schema", "diffexp2-retained-local-v2"},
        {"handle", handle_},
        {"source_chart", source_chart_},
        {"source_operator_identity", source_operator_identity_},
        {"scalar_domain", scalar_domain()},
        {"precision_bits", precision_bits_},
        {"solution", checkpoint_local_solution_record(solution_)},
        {"pseudo_hits", checkpoint_pseudo_hits_record(pseudo_hits_)},
        {"diagnostics",
         json::object{{"top_valid", encode_validity(top_valid_)},
                      {"create_parse_ms", create_parse_ms_},
                      {"create_kernel_ms", create_kernel_ms_}}},
        {"runtime_stats", std::move(runtime)},
        {"column_provenance", column_provenance_.has_value()
             ? json::value(column_provenance_->encode())
             : json::value(nullptr)}};
      if (serialize_derivation_checkpoint_fields_) {
        record["retained_derivation"] = retained_derivation_.has_value()
            ? json::value(*retained_derivation_) : json::value(nullptr);
        record["retained_owner_lineage"] = std::move(owner_lineage);
      }
      if (serialize_tail_checkpoint_fields_) {
        if (tail_model_.model.has_value()) {
          if (tail_model_.status != TailMajorantStatus::Certified ||
              tail_checkpoint_marker_.has_value())
            throw std::logic_error(
                "attached regular tail model has inconsistent checkpoint status");
          tail_majorant_detail::validate_restored_regular_taylor_tail_model(
              *tail_model_.model, solution_, source_operator_identity_);
          record["tail_model_restore"] = json::object{
              {"capability", kRegularTailMajorantCapability},
              {"serialized", true},
              {"status", "certified"},
              {"attached_before_save", true},
              {"model", checkpoint_regular_tail_model_record(
                   *tail_model_.model)}};
        } else {
          record["tail_model_restore"] = json::object{
              {"capability", kRegularTailMajorantCapability},
              {"serialized", false},
              {"status", tail_checkpoint_marker_.has_value()
                   ? tail_checkpoint_marker_->saved_status
                   : tail_majorant_status_name(tail_model_.status)},
              {"attached_before_save", tail_checkpoint_marker_.has_value()
                   ? tail_checkpoint_marker_->attached_before_save
                   : false}};
        }
      }
      return record;
    }
  }

  void restore_runtime_stats(const StoredLocalStats& state) {
    if (state.coefficient_count != coefficient_count())
      throw std::invalid_argument(
          "checkpoint local coefficient count does not match its tensor");
    evaluations_.store(state.evaluations);
    residual_certifications_.store(state.residual_certifications);
    endpoint_limits_.store(state.endpoint_limits);
    line_integrations_.store(state.line_integrations);
    tail_certificate_requests_.store(state.tail_certificate_requests);
    tail_certificate_certified_.store(state.tail_certificate_certified);
    tail_certificate_inconclusive_.store(
        state.tail_certificate_inconclusive);
    tail_certificate_unsupported_.store(state.tail_certificate_unsupported);
    std::lock_guard<std::mutex> lock(stats_mutex_);
    evaluate_ms_ = state.evaluate_ms;
    residual_certify_ms_ = state.residual_certify_ms;
    endpoint_limit_ms_ = state.endpoint_limit_ms;
    line_integration_ms_ = state.line_integration_ms;
  }

  const LocalSolution<Scalar>& solution() const { return solution_; }
  const RegularTaylorTailModelResult& tail_model() const {
    return tail_model_;
  }
  std::int32_t top_valid() const { return top_valid_; }
  const std::optional<json::object>& retained_derivation() const override {
    return retained_derivation_;
  }
  std::shared_ptr<void> retained_derivation_owner() const override {
    return retained_owner_;
  }
  const std::vector<PseudoHit<Scalar>>& pseudo_hits() const {
    return pseudo_hits_;
  }

 private:
  json::object metadata_json() const {
    json::array prescriptions;
    for (const auto& prescription : solution_.prescriptions) {
      prescriptions.push_back(json::object{
          {"factor_exact", prescription.factor_exact},
          {"sign", prescription.sign},
          {"multiplicity", prescription.multiplicity},
          {"leading_coefficient_sign",
           prescription.leading_coefficient_sign}});
    }
    json::array sectors;
    sectors.reserve(solution_.sectors.size());
    for (const auto& sector : solution_.sectors)
      sectors.push_back(json::object{
          {"a", encode_exact_descriptor(sector.a)},
          {"b", encode_exact_descriptor(sector.b)},
          {"log_power", sector.log_power}});
    const auto& sector = solution_.sectors.front();
    return json::object{
        {"chart", json::object{
            {"center_exact", solution_.chart.center_exact},
            {"scale_exact", solution_.chart.scale_exact},
            {"infinite_radius", solution_.chart.infinite_radius},
            {"radius_ball", encode_scalar(solution_.chart.radius, 30)}}},
        {"tag", json::object{{"a", encode_exact_descriptor(sector.a)},
                             {"b", encode_exact_descriptor(sector.b)}}},
        {"sectors", std::move(sectors)},
        {"prescriptions", std::move(prescriptions)}};
  }

  std::size_t coefficient_count() const {
    std::size_t count = 0;
    for (const auto& sector : solution_.sectors)
      count += sector.coefficients.size();
    return count;
  }

  json::object rational_row_owner_lineage() const {
    if (!retained_derivation_.has_value() || retained_owner_ == nullptr)
      throw std::logic_error(
          "rational-row local lost its derivation or source owner");
    const auto& derivation = *retained_derivation_;
    require_exact_keys(
        derivation,
        {"schema", "capability", "source", "row", "output",
         "analytic_prescriptions", "coefficient_transport",
         "provenance_identity"},
        "retained rational-row local derivation");
    if (required_string(derivation, "schema") !=
            "diffexp2-retained-rational-row-local-application-v1" ||
        required_string(derivation, "capability") !=
            "retained-native-rational-row-local-application-v1" ||
        required_string(derivation, "analytic_prescriptions") !=
            "preserved-exactly" ||
        required_string(derivation, "coefficient_transport") !=
            "native-retained-only")
      throw std::invalid_argument(
          "retained rational-row derivation changes its certified scope");

    auto erased_source =
        std::static_pointer_cast<StoredLocalBase>(retained_owner_);
    auto source = std::dynamic_pointer_cast<StoredLocal<Scalar>>(
        erased_source);
    if (!source || std::string(source->scalar_domain()) != scalar_domain())
      throw std::invalid_argument(
          "retained rational-row derivation source domain is inconsistent");
    const auto& source_record = as_object(
        derivation.at("source"), "retained rational-row source");
    require_exact_keys(
        source_record,
        {"local", "chart", "source_operator_identity",
         "checkpoint_identity", "dimension", "epsilon",
         "taylor_complete_max"},
        "retained rational-row source");
    const auto& source_epsilon = as_object(
        source_record.at("epsilon"), "retained rational-row source epsilon");
    require_exact_keys(source_epsilon, {"min", "max"},
                       "retained rational-row source epsilon");
    if (required_string(source_record, "local") != source->handle() ||
        required_string(source_record, "chart") != source->source_chart() ||
        required_string(source_record, "source_operator_identity") !=
            source->source_operator_identity() ||
        required_string(source_record, "checkpoint_identity") !=
            source->checkpoint_identity() ||
        as_u32(source_record.at("dimension"),
               "rational-row source dimension") !=
            source->solution().dimension ||
        as_i32(source_epsilon.at("min"),
               "rational-row source epsilon minimum") !=
            source->solution().epsilon.min_power ||
        as_i32(source_epsilon.at("max"),
               "rational-row source epsilon maximum") !=
            source->solution().epsilon.complete_max ||
        as_u32(source_record.at("taylor_complete_max"),
               "rational-row source Taylor maximum") !=
            source->solution().taylor_complete_max)
      throw std::invalid_argument(
          "retained rational-row source provenance disagrees with its strong owner");

    const auto& row = as_object(
        derivation.at("row"), "retained rational-row identity");
    require_exact_keys(
        row, {"exact_identity", "columns", "active_entries",
              "structurally_zero"},
        "retained rational-row identity");
    const auto row_identity = required_string(row, "exact_identity");
    if (as_u32(row.at("columns"), "retained rational-row columns") !=
        source->solution().dimension)
      throw std::invalid_argument(
          "retained rational-row column count differs from its source");
    if (!row.at("structurally_zero").is_bool())
      throw std::invalid_argument(
          "retained rational-row zero fact must be Boolean");
    const auto& active = as_array(
        row.at("active_entries"), "retained rational-row active entries");
    if (row.at("structurally_zero").as_bool() != active.empty())
      throw std::invalid_argument(
          "retained rational-row zero fact disagrees with its active entries");
    std::optional<std::uint32_t> previous_column;
    for (const auto& raw_entry : active) {
      const auto& entry = as_object(
          raw_entry, "retained rational-row active entry");
      require_exact_keys(
          entry, {"column", "epsilon_shift", "center_pole_order",
                  "exact_identity"},
          "retained rational-row active entry");
      const auto column = as_u32(
          entry.at("column"), "retained rational-row active column");
      if (column >= source->solution().dimension ||
          (previous_column.has_value() && *previous_column >= column))
        throw std::invalid_argument(
            "retained rational-row active columns are not canonical");
      previous_column = column;
      (void)as_i32(entry.at("epsilon_shift"),
                   "retained rational-row epsilon shift");
      (void)as_u32(entry.at("center_pole_order"),
                   "retained rational-row pole order");
      (void)required_string(entry, "exact_identity");
    }

    const auto& output = as_object(
        derivation.at("output"), "retained rational-row output");
    require_exact_keys(
        output, {"checkpoint_identity", "dimension", "epsilon",
                 "taylor_complete_max"},
        "retained rational-row output");
    const auto& output_epsilon = as_object(
        output.at("epsilon"), "retained rational-row output epsilon");
    require_exact_keys(output_epsilon, {"min", "max"},
                       "retained rational-row output epsilon");
    if (required_string(output, "checkpoint_identity") !=
            solution_.checkpoint_identity ||
        as_u32(output.at("dimension"),
               "rational-row output dimension") != solution_.dimension ||
        as_i32(output_epsilon.at("min"),
               "rational-row output epsilon minimum") !=
            solution_.epsilon.min_power ||
        as_i32(output_epsilon.at("max"),
               "rational-row output epsilon maximum") !=
            solution_.epsilon.complete_max ||
        as_u32(output.at("taylor_complete_max"),
               "rational-row output Taylor maximum") !=
            solution_.taylor_complete_max)
      throw std::invalid_argument(
          "retained rational-row output provenance disagrees with its tensor");
    if (source_chart_ != source->source_chart() ||
        !local_algebra_detail::same_chart(
            solution_.chart, source->solution().chart) ||
        !local_algebra_detail::same_prescriptions(
            solution_.prescriptions, source->solution().prescriptions))
      throw std::invalid_argument(
          "retained rational-row output left its source analytic chart");

    auto identity_input = derivation;
    const auto derivation_identity = required_string(
        derivation, "provenance_identity");
    identity_input.erase("provenance_identity");
    if (json::serialize(canonical_json_value(identity_input)) !=
        derivation_identity)
      throw std::invalid_argument(
          "retained rational-row derivation identity is inconsistent");
    const json::object operator_provenance{
        {"schema", "diffexp2-rational-row-derived-operator-v1"},
        {"source_operator_identity", source->source_operator_identity()},
        {"row_exact_identity", row_identity},
        {"provenance_identity", derivation_identity}};
    if (json::serialize(canonical_json_value(operator_provenance)) !=
        source_operator_identity_)
      throw std::invalid_argument(
          "retained rational-row derived operator identity is inconsistent");

    return json::object{
        {"source_local", source->handle()},
        {"source_chart", source->source_chart()},
        {"source_operator_identity", source->source_operator_identity()},
        {"source_checkpoint_identity", source->checkpoint_identity()},
        {"row_exact_identity", row_identity},
        {"derivation_provenance_identity", derivation_identity},
        {"derived_operator_identity", source_operator_identity_}};
  }

  LocalSolution<Scalar> solution_;
  slong precision_bits_ = 256;
  std::vector<PseudoHit<Scalar>> pseudo_hits_;
  std::int32_t top_valid_ = kCompleteInfinity;
  std::optional<json::object> retained_derivation_;
  std::shared_ptr<void> retained_owner_;
  RegularTaylorTailModelResult tail_model_ = unavailable_tail_model(
      "tail model is unavailable for this retained local");
  std::optional<TailModelCheckpointMarker> tail_checkpoint_marker_;
  bool serialize_tail_checkpoint_fields_ = true;
  bool serialize_derivation_checkpoint_fields_ = true;
  std::atomic<std::uint64_t> evaluations_{0};
  std::atomic<std::uint64_t> residual_certifications_{0};
  std::atomic<std::uint64_t> endpoint_limits_{0};
  std::atomic<std::uint64_t> line_integrations_{0};
  std::atomic<std::uint64_t> tail_certificate_requests_{0};
  std::atomic<std::uint64_t> tail_certificate_certified_{0};
  std::atomic<std::uint64_t> tail_certificate_inconclusive_{0};
  std::atomic<std::uint64_t> tail_certificate_unsupported_{0};
  mutable std::mutex stats_mutex_;
  double evaluate_ms_ = 0.0;
  double residual_certify_ms_ = 0.0;
  double endpoint_limit_ms_ = 0.0;
  double line_integration_ms_ = 0.0;
};

constexpr const char* kRetainedEndpointLimitCapability =
    "retained-native-endpoint-sector-limit-v1";
constexpr const char* kRetainedPlannedEndpointLimitCapability =
    "retained-native-plan-bound-endpoint-sector-limit-v1";
constexpr const char* kRetainedRationalRowCapability =
    "retained-native-rational-row-local-application-v1";

template <typename Scalar>
PreparedSparseLocalMultiplierMatrix<Scalar> parse_prepared_rational_row(
    const json::value& raw, const LocalSolution<Scalar>& source) {
  const auto& row = as_object(raw, "prepared rational local row");
  require_exact_keys(row,
      {"schema", "columns", "exact_identity", "entries"},
      "prepared rational local row");
  if (required_string(row, "schema") !=
      "diffexp2-prepared-rational-local-row-v1")
    throw std::invalid_argument(
        "unsupported prepared rational local-row schema");
  const auto columns = as_u32(
      row.at("columns"), "prepared rational local-row columns");
  if (columns == 0 || columns != source.dimension)
    throw std::invalid_argument(
        "prepared rational local-row dimension differs from its retained local");

  PreparedSparseLocalMultiplierMatrix<Scalar> matrix;
  matrix.rows = 1;
  matrix.columns = columns;
  matrix.exact_identity = required_string(row, "exact_identity");
  if (matrix.exact_identity.empty())
    throw std::invalid_argument(
        "prepared rational local-row identity must be nonempty");

  std::optional<std::uint32_t> previous_column;
  for (const auto& raw_entry : as_array(
           row.at("entries"), "prepared rational local-row entries")) {
    const auto& entry = as_object(
        raw_entry, "prepared rational local-row entry");
    require_exact_keys(entry, {"column", "multiplier"},
                       "prepared rational local-row entry");
    const auto column = as_u32(
        entry.at("column"), "prepared rational local-row column");
    if (column >= columns ||
        (previous_column.has_value() && *previous_column >= column))
      throw std::invalid_argument(
          "prepared rational local-row columns must be unique, in range, and strictly increasing");
    previous_column = column;

    const auto& raw_multiplier = as_object(
        entry.at("multiplier"), "prepared rational local-row multiplier");
    require_exact_keys(raw_multiplier,
        {"epsilon_shift", "center_pole_order", "kernels",
         "exact_identity", "proven_zero"},
        "prepared rational local-row multiplier");
    if (!raw_multiplier.at("proven_zero").is_bool())
      throw std::invalid_argument(
          "prepared rational local-row structural-zero fact must be boolean");
    if (raw_multiplier.at("proven_zero").as_bool())
      throw std::invalid_argument(
          "structurally zero rational-row entries must be omitted");

    PreparedRationalTaylorMultiplier<Scalar> multiplier;
    multiplier.epsilon_shift = as_i32(
        raw_multiplier.at("epsilon_shift"),
        "prepared rational local-row epsilon shift");
    multiplier.center_pole_order = as_u32(
        raw_multiplier.at("center_pole_order"),
        "prepared rational local-row center-pole order");
    multiplier.exact_identity = required_string(
        raw_multiplier, "exact_identity");
    if (multiplier.exact_identity.empty())
      throw std::invalid_argument(
          "prepared rational local-row multiplier identity must be nonempty");
    multiplier.proven_zero = false;

    const auto& raw_kernels = as_array(
        raw_multiplier.at("kernels"),
        "prepared rational local-row epsilon kernels");
    if (raw_kernels.size() < source.epsilon.width())
      throw std::invalid_argument(
          "prepared rational local-row multiplier does not cover the exact source epsilon width");
    multiplier.kernels.reserve(raw_kernels.size());
    for (const auto& raw_kernel : raw_kernels) {
      const auto& coefficients = as_array(
          raw_kernel, "prepared rational local-row Taylor kernel");
      if (coefficients.size() < source.taylor_width())
        throw std::invalid_argument(
          "prepared rational local-row multiplier does not cover the exact source Taylor width");
      std::vector<Scalar> kernel;
      kernel.reserve(coefficients.size());
      for (const auto& coefficient : coefficients)
        kernel.push_back(parse_scalar<Scalar>(coefficient));
      multiplier.kernels.push_back(std::move(kernel));
    }
    matrix.entries.push_back(
        typename PreparedSparseLocalMultiplierMatrix<Scalar>::Entry{
            0, column, std::move(multiplier)});
  }
  return matrix;
}

template <typename Scalar>
LocalSolution<Scalar> exact_zero_scalar_local_like(
    const LocalSolution<Scalar>& source,
    const std::string& checkpoint_identity) {
  auto result = local_algebra_detail::with_selected_component(source, 0);
  for (auto& sector : result.sectors)
    std::fill(sector.coefficients.begin(), sector.coefficients.end(),
              ScalarTraits<Scalar>::zero());
  result.checkpoint_identity = checkpoint_identity;
  validate_local_solution(result, false);
  return result;
}

std::int32_t shifted_local_validity(std::int32_t validity,
                                    std::int32_t shift) {
  if (validity == kCompleteInfinity) return kCompleteInfinity;
  return local_algebra_detail::checked_i32(
      static_cast<std::int64_t>(validity) + shift,
      "rational-row output validity");
}

template <typename Scalar>
std::shared_ptr<StoredLocalBase> build_rational_row_local(
    const std::string& local_handle, const json::object& request,
    slong precision_bits,
    const std::shared_ptr<StoredLocal<Scalar>>& source,
    const std::shared_ptr<StoredLocalBase>& erased_source) {
  if (source == nullptr || erased_source == nullptr ||
      source.get() != erased_source.get())
    throw std::logic_error(
        "retained rational-row application lost typed source ownership");
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  const auto source_checkpoint_identity = required_string(
      request, "source_checkpoint_identity");
  if (checkpoint_identity.empty() || source_checkpoint_identity.empty())
    throw std::invalid_argument(
        "retained rational-row checkpoint identities must be nonempty");
  if (source_checkpoint_identity != source->checkpoint_identity())
    throw std::invalid_argument(
        "rational-row source checkpoint identity differs from its retained local");
  if (!source->solution().error.empty())
    throw std::domain_error(
        "native rational-row application requires explicit source error-envelope propagation");

  std::unique_ptr<AcbPrecisionLease> acb_lease;
  if constexpr (std::is_same_v<Scalar, ComplexBall>) {
    acb_lease = std::make_unique<AcbPrecisionLease>(precision_bits);
    ComplexBall::set_precision(precision_bits);
  }
  const auto started = std::chrono::steady_clock::now();
  auto matrix = parse_prepared_rational_row<Scalar>(
      request.at("row"), source->solution());

  json::array entry_provenance;
  entry_provenance.reserve(matrix.entries.size());
  std::int32_t output_top_valid = matrix.entries.empty()
      ? source->top_valid() : kCompleteInfinity;
  for (const auto& entry : matrix.entries) {
    output_top_valid = std::min(
        output_top_valid,
        shifted_local_validity(source->top_valid(),
                               entry.multiplier.epsilon_shift));
    entry_provenance.push_back(json::object{
        {"column", entry.column},
        {"epsilon_shift", entry.multiplier.epsilon_shift},
        {"center_pole_order", entry.multiplier.center_pole_order},
        {"exact_identity", entry.multiplier.exact_identity}});
  }

  auto applied = apply_prepared_sparse_local_matrix(
      matrix, source->solution(), checkpoint_identity);
  auto solution = applied.has_value()
      ? std::move(*applied)
      : exact_zero_scalar_local_like(
            source->solution(), checkpoint_identity);
  output_top_valid = std::min(
      output_top_valid, solution.epsilon.complete_max);
  if (output_top_valid < solution.epsilon.min_power)
    throw std::domain_error(
        "rational-row application has no valid output epsilon coefficient");
  if (output_top_valid < solution.epsilon.complete_max)
    solution = restrict_local_epsilon_frame_strict_lower(
        solution, solution.epsilon.min_power, output_top_valid,
        checkpoint_identity);
  if (solution.dimension != 1)
    throw std::logic_error(
        "rational-row application did not produce a scalar local solution");

  json::object derivation{
      {"schema", "diffexp2-retained-rational-row-local-application-v1"},
      {"capability", kRetainedRationalRowCapability},
      {"source", json::object{
           {"local", source->handle()},
           {"chart", source->source_chart()},
           {"source_operator_identity",
            source->source_operator_identity()},
           {"checkpoint_identity", source_checkpoint_identity},
           {"dimension", source->solution().dimension},
           {"epsilon", json::object{
                {"min", source->solution().epsilon.min_power},
                {"max", source->solution().epsilon.complete_max}}},
           {"taylor_complete_max",
            source->solution().taylor_complete_max}}},
      {"row", json::object{
           {"exact_identity", matrix.exact_identity},
           {"columns", matrix.columns},
           {"active_entries", std::move(entry_provenance)},
           {"structurally_zero", matrix.entries.empty()}}},
      {"output", json::object{
           {"checkpoint_identity", checkpoint_identity},
           {"dimension", solution.dimension},
           {"epsilon", json::object{
                {"min", solution.epsilon.min_power},
                {"max", solution.epsilon.complete_max}}},
           {"taylor_complete_max", solution.taylor_complete_max}}},
      {"analytic_prescriptions", "preserved-exactly"},
      {"coefficient_transport", "native-retained-only"}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(derivation));
  derivation["provenance_identity"] = provenance_identity;

  json::object operator_provenance{
      {"schema", "diffexp2-rational-row-derived-operator-v1"},
      {"source_operator_identity", source->source_operator_identity()},
      {"row_exact_identity", matrix.exact_identity},
      {"provenance_identity", provenance_identity}};
  const auto derived_operator_identity = json::serialize(
      canonical_json_value(operator_provenance));
  const auto elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  NativeLocalDiagnostics diagnostics;
  diagnostics.top_valid = output_top_valid;
  diagnostics.kernel_ms = elapsed_ms;
  return make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
      local_handle, source->source_chart(), derived_operator_identity,
      std::move(solution), precision_bits,
      std::vector<PseudoHit<Scalar>>{}, diagnostics, std::nullopt,
      std::move(derivation),
      std::static_pointer_cast<void>(erased_source));
}

struct ParsedEndpointLimitPolicy {
  EndpointLimitOptions options;
  std::string cancellation_mode;
  std::optional<std::int32_t> requested_rim;
};

struct ParsedEndpointCancellationPolicy {
  bool allow_certified_numeric_cancellation = false;
  std::string cancellation_mode;
};

ParsedEndpointCancellationPolicy parse_endpoint_cancellation_policy(
    const json::object& request) {
  ParsedEndpointCancellationPolicy parsed;
  const auto& cancellation = as_object(
      request.at("cancellation"), "endpoint cancellation policy");
  parsed.cancellation_mode = required_string(cancellation, "mode");
  if (parsed.cancellation_mode == "exact-coefficient-field") {
    parsed.allow_certified_numeric_cancellation = false;
  } else if (parsed.cancellation_mode == "exact-or-acb-singleton") {
    parsed.allow_certified_numeric_cancellation = true;
  } else {
    throw std::invalid_argument(
        "endpoint cancellation mode must be exact-coefficient-field or "
        "exact-or-acb-singleton");
  }
  if (cancellation.size() != 1)
    throw std::invalid_argument(
        "endpoint cancellation policy accepts only its exact mode; "
        "tolerance-based cancellation is unsupported");
  return parsed;
}

ParsedEndpointLimitPolicy parse_endpoint_limit_policy(
    const json::object& request) {
  ParsedEndpointLimitPolicy parsed;
  parsed.options.approach_direction = as_i32(
      request.at("approach_direction"), "endpoint approach direction");
  if (parsed.options.approach_direction != 1 &&
      parsed.options.approach_direction != -1)
    throw std::invalid_argument(
        "endpoint approach direction must be exactly +1 or -1");

  if (const auto* raw_rim = request.if_contains("rim");
      raw_rim != nullptr && !raw_rim->is_null()) {
    const auto rim = as_i32(*raw_rim, "endpoint rim");
    if (rim != 1 && rim != -1)
      throw std::invalid_argument("endpoint rim must be exactly +1 or -1");
    parsed.requested_rim = rim;
    parsed.options.imaginary_sign = rim;
  }

  const auto cancellation = parse_endpoint_cancellation_policy(request);
  parsed.cancellation_mode = cancellation.cancellation_mode;
  parsed.options.allow_certified_numeric_cancellation =
      cancellation.allow_certified_numeric_cancellation;
  return parsed;
}

EpsilonWindow endpoint_value_window(const EndpointLimitResult& result) {
  if (result.values.empty())
    throw std::logic_error("retained endpoint result has no components");
  const auto window = result.values.front().window();
  for (const auto& component : result.values)
    if (component.window().min_power != window.min_power ||
        component.window().complete_max != window.complete_max)
      throw std::logic_error(
          "retained endpoint components have unequal epsilon windows");
  return window;
}

EpsilonVector endpoint_values_vector(const EndpointLimitResult& result) {
  const auto window = endpoint_value_window(result);
  EpsilonVector vector;
  vector.epsilon = window;
  vector.dimension = static_cast<std::uint32_t>(result.values.size());
  vector.coefficients.reserve(window.width() * vector.dimension);
  for (std::size_t ei = 0; ei < window.width(); ++ei)
    for (const auto& component : result.values)
      vector.coefficients.push_back(component.coefficients().at(ei));
  return vector;
}

class StoredEndpointResult {
 public:
  StoredEndpointResult(
      std::string handle, std::string checkpoint_identity,
      std::string provenance_identity, std::string source_local,
      std::string source_chart, std::string source_operator_identity,
      std::string source_checkpoint,
      std::string source_domain, std::int32_t approach_direction,
      std::optional<std::int32_t> requested_rim,
      std::string cancellation_mode, json::object analytic_metadata,
      EndpointLimitResult&& result, double elapsed_ms,
      std::optional<json::object> planned_source = std::nullopt,
      std::optional<std::int32_t> planned_effective_rim = std::nullopt,
      std::shared_ptr<StoredTilePlan> plan_owner = nullptr,
      std::shared_ptr<StoredLocalBase> local_owner = nullptr)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        source_local_(std::move(source_local)),
        source_chart_(std::move(source_chart)),
        source_operator_identity_(std::move(source_operator_identity)),
        source_checkpoint_(std::move(source_checkpoint)),
        source_domain_(std::move(source_domain)),
        approach_direction_(approach_direction),
        requested_rim_(requested_rim),
        cancellation_mode_(std::move(cancellation_mode)),
        analytic_metadata_(std::move(analytic_metadata)),
        result_(std::move(result)), elapsed_ms_(elapsed_ms),
        planned_source_(std::move(planned_source)),
        planned_effective_rim_(planned_effective_rim),
        plan_owner_(std::move(plan_owner)),
        local_owner_(std::move(local_owner)) {
    // Validate the retained public frame once, before publishing its handle.
    (void)endpoint_value_window(result_);
    if (planned_source_.has_value()
            ? (plan_owner_ == nullptr || local_owner_ == nullptr)
            : (plan_owner_ != nullptr || local_owner_ != nullptr))
      throw std::invalid_argument(
          "plan-bound endpoint result must retain both its plan and local owners");
  }

  const std::string& handle() const { return handle_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& provenance_identity() const {
    return provenance_identity_;
  }
  double elapsed_ms() const { return elapsed_ms_; }
  bool plan_bound() const { return planned_source_.has_value(); }
  const std::shared_ptr<StoredTilePlan>& plan_owner() const {
    return plan_owner_;
  }
  const std::shared_ptr<StoredLocalBase>& local_owner() const {
    return local_owner_;
  }

  json::object summary() const {
    const auto window = endpoint_value_window(result_);
    const auto cancellation_scope = source_domain_ == "rational"
        ? "exact-rational"
        : cancellation_mode_ == "exact-or-acb-singleton"
            ? "acb-exact-singleton-zero"
            : "exact-coefficient-field-only";
    json::object source{
        {"local", source_local_}, {"chart", source_chart_},
        {"source_operator_identity", source_operator_identity_},
        {"checkpoint_identity", source_checkpoint_},
        {"coefficient_domain", source_domain_}};
    if (planned_source_.has_value()) source = *planned_source_;
    json::object result{
        {"endpoint", handle_},
        {"capability", plan_bound()
             ? kRetainedPlannedEndpointLimitCapability
             : kRetainedEndpointLimitCapability},
        {"native_retained", true},
        {"retained_state", "specialized-acb-epsilon-vector"},
        {"json_coefficients", 0},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"execution_scope", plan_bound()
             ? "plan-bound-final-arm-endpoint"
             : "unplanned-low-level"},
        {"source", std::move(source)},
        {"dimension", result_.values.size()},
        {"epsilon_min", window.min_power},
        {"epsilon_max", window.complete_max},
        {"coefficient_field", "acb-specialized"},
        {"arithmetic_enclosed", true},
        {"approach_direction", approach_direction_},
        {"cancellation", json::object{
             {"mode", cancellation_mode_},
             {"effective_scope", cancellation_scope},
             {"numeric_singleton_cancellations",
              result_.cancelled_divergent_coefficients}}},
        {"analytic_regularization", json::object{
             {"regulator_slope_scope",
              "exact-zero-fact; certified-nonzero symbolic slopes allowed"},
             {"unregulated_power_scope", "exact-rational"},
             {"endpoint_rule", "drop-exact-nonzero-regulator-slope"},
              {"dropped_regulated_sectors",
               result_.dropped_regulated_sectors},
              {"metadata", analytic_metadata_}}},
        {"elapsed_ms", elapsed_ms_}};
    result["requested_rim"] = requested_rim_.has_value()
        ? json::value(*requested_rim_) : json::value(nullptr);
    if (plan_bound()) {
      result["derived_rim"] = planned_effective_rim_.has_value()
          ? json::value(*planned_effective_rim_) : json::value(nullptr);
      result["effective_rim"] = planned_effective_rim_.has_value()
          ? json::value(*planned_effective_rim_) : json::value(nullptr);
      result["rim_source"] =
          "final-chart-exact-odd-multiplicity-prescriptions";
    } else {
      result["effective_rim"] = result_.imaginary_sign;
      result["rim_source"] = "unplanned-caller-or-principal-default";
    }
    return result;
  }

  json::object stats_json() const {
    auto out = summary();
    std::lock_guard<std::mutex> lock(stats_mutex_);
    out["exports"] = exports_;
    out["export_ms"] = export_ms_;
    return out;
  }

  json::object checkpoint_record() const {
    json::array values;
    values.reserve(result_.values.size());
    for (const auto& value : result_.values)
      values.push_back(checkpoint_epsilon_frame_record(value));
    std::lock_guard<std::mutex> lock(stats_mutex_);
    json::object source{
        {"local", source_local_}, {"chart", source_chart_},
        {"source_operator_identity", source_operator_identity_},
        {"checkpoint_identity", source_checkpoint_},
        {"coefficient_domain", source_domain_}};
    if (planned_source_.has_value()) source = *planned_source_;
    json::object record{
        {"schema", plan_bound()
             ? "diffexp2-retained-plan-bound-endpoint-result-v1"
             : "diffexp2-retained-endpoint-result-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"source", std::move(source)},
        {"approach_direction", approach_direction_},
        {"cancellation_mode", cancellation_mode_},
        {"analytic_metadata", analytic_metadata_},
        {"result", json::object{
            {"values", std::move(values)},
            {"dropped_regulated_sectors",
             result_.dropped_regulated_sectors},
             {"cancelled_divergent_coefficients",
              result_.cancelled_divergent_coefficients},
            {"imaginary_sign", plan_bound()
                 ? (planned_effective_rim_.has_value()
                       ? json::value(*planned_effective_rim_)
                       : json::value(nullptr))
                 : json::value(result_.imaginary_sign)}}},
        {"elapsed_ms", elapsed_ms_},
        {"runtime_stats", json::object{{"exports", exports_},
                                        {"export_ms", export_ms_}}}};
    if (plan_bound()) {
      record["derived_rim"] = planned_effective_rim_.has_value()
          ? json::value(*planned_effective_rim_) : json::value(nullptr);
    } else {
      record["requested_rim"] = requested_rim_.has_value()
          ? json::value(*requested_rim_) : json::value(nullptr);
    }
    return record;
  }

  void restore_runtime_stats(std::uint64_t exports, double export_ms) {
    std::lock_guard<std::mutex> lock(stats_mutex_);
    exports_ = exports;
    export_ms_ = export_ms;
  }

  json::object export_values(const std::string& expected_checkpoint,
                             int output_digits) {
    if (expected_checkpoint.empty())
      throw std::invalid_argument(
          "endpoint export checkpoint identity must be nonempty");
    if (expected_checkpoint != checkpoint_identity_)
      throw std::invalid_argument(
          "endpoint export checkpoint identity does not match retained state");
    const auto started = std::chrono::steady_clock::now();
    auto encoded = encode_epsilon_vector(
        endpoint_values_vector(result_), output_digits);
    const auto elapsed = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    {
      std::lock_guard<std::mutex> lock(stats_mutex_);
      ++exports_;
      export_ms_ += elapsed;
    }
    return json::object{
        {"endpoint", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"compatibility_export", true},
        {"coefficient_field", "acb-specialized"},
        {"json_coefficients", encoded.at("coefficients").as_array().size()},
        {"value", std::move(encoded)},
        {"elapsed_ms", elapsed}};
  }

 private:
  std::string handle_;
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::string source_local_;
  std::string source_chart_;
  std::string source_operator_identity_;
  std::string source_checkpoint_;
  std::string source_domain_;
  std::int32_t approach_direction_ = 1;
  std::optional<std::int32_t> requested_rim_;
  std::string cancellation_mode_;
  json::object analytic_metadata_;
  EndpointLimitResult result_;
  double elapsed_ms_ = 0.0;
  std::optional<json::object> planned_source_;
  std::optional<std::int32_t> planned_effective_rim_;
  std::shared_ptr<StoredTilePlan> plan_owner_;
  std::shared_ptr<StoredLocalBase> local_owner_;
  mutable std::mutex stats_mutex_;
  std::uint64_t exports_ = 0;
  double export_ms_ = 0.0;
};

std::shared_ptr<StoredEndpointResult> build_endpoint_limit(
    const std::string& endpoint_handle, const json::object& request,
    const std::shared_ptr<StoredLocalBase>& local) {
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  const auto expected_source_checkpoint = required_string(
      request, "source_checkpoint_identity");
  if (checkpoint_identity.empty() || expected_source_checkpoint.empty())
    throw std::invalid_argument(
        "endpoint checkpoint identities must be nonempty");
  if (expected_source_checkpoint != local->checkpoint_identity())
    throw std::invalid_argument(
        "endpoint source checkpoint identity does not match retained local");
  const auto policy = parse_endpoint_limit_policy(request);
  auto analytic_metadata = local->exact_analytic_metadata();
  json::object provenance{
      {"schema", "diffexp2-retained-native-endpoint-sector-limit-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"source", json::object{
           {"local", local->handle()}, {"chart", local->source_chart()},
           {"source_operator_identity",
            local->source_operator_identity()},
           {"checkpoint_identity", expected_source_checkpoint},
           {"coefficient_domain", local->scalar_domain()}}},
      {"approach_direction", policy.options.approach_direction},
      {"rim", policy.requested_rim.has_value()
           ? json::value(*policy.requested_rim) : json::value(nullptr)},
      {"cancellation", json::object{{"mode", policy.cancellation_mode}}},
      {"analytic_metadata", analytic_metadata}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto started = std::chrono::steady_clock::now();
  auto result = local->endpoint_limit(policy.options);
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredEndpointResult>(
      endpoint_handle, checkpoint_identity, provenance_identity,
      local->handle(), local->source_chart(),
      local->source_operator_identity(), expected_source_checkpoint,
      local->scalar_domain(), policy.options.approach_direction,
      policy.requested_rim, policy.cancellation_mode,
      std::move(analytic_metadata), std::move(result), elapsed);
}

constexpr const char* kExactRegularLocalMatchCapability =
    "exact-rational-regular-local-match-v1";
constexpr const char* kRefinedAcbLocalMatchCapability =
    "exact-lattice-guided-acb-local-match-v1";
constexpr const char* kExactEvaluatedLatticeSchema =
    "diffexp2-exact-evaluated-epsilon-lattice-v1";
constexpr const char* kNativeUnitSaturationRequestSchema =
    "diffexp2-native-acb-unit-leading-saturation-request-v1";
constexpr const char* kNativeUnitSaturationProofSchema =
    "diffexp2-native-acb-unit-leading-saturation-proof-v1";
constexpr const char* kNativeSingularSCCSaturationRequestSchema =
    "diffexp2-native-acb-singular-scc-valuation-zero-saturation-request-v1";
constexpr const char* kNativeSingularSCCSaturationProofSchema =
    "diffexp2-native-acb-singular-scc-valuation-zero-saturation-proof-v1";
constexpr const char* kAcbSingularScalarSCCColumnCapability =
    "acb-regular-singular-scalar-block-dag-column-v1";
constexpr const char* kAcbSingularJordanSCCColumnCapability =
    "acb-regular-singular-jordan-block-dag-column-v1";

class StoredMatchBase {
 public:
  explicit StoredMatchBase(std::string handle) : handle_(std::move(handle)) {}
  virtual ~StoredMatchBase() = default;

  virtual json::object summary() const = 0;
  virtual json::object checkpoint_record() const = 0;
  const std::string& handle() const { return handle_; }

 protected:
  std::string handle_;
};

class StoredExactRegularMatch final : public StoredMatchBase {
 public:
  StoredExactRegularMatch(
      std::string handle, std::string checkpoint_identity,
      std::string provenance_identity,
      std::vector<json::object> basis_sources, json::object incoming_source,
      std::string basis_chart,
      std::string incoming_chart, std::string basis_point,
      std::string incoming_point, std::string physical_point,
      EpsilonWindow requested_window, std::int32_t required_complete_max,
      std::uint32_t dimension,
      ExactLaurentMatrix<Rational>&& transformation,
      FiniteLaurentVector<Rational>&& weights,
      EpsilonLatticeSaturationDiagnostics<Rational>&& diagnostics,
      EpsilonWindow residual_window, double elapsed_ms,
      std::vector<std::shared_ptr<StoredLocalBase>> basis_owners,
      std::shared_ptr<StoredLocalBase> incoming_owner)
      : StoredMatchBase(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        basis_sources_(std::move(basis_sources)),
        incoming_source_(std::move(incoming_source)),
        basis_chart_(std::move(basis_chart)),
        incoming_chart_(std::move(incoming_chart)),
        basis_point_(std::move(basis_point)),
        incoming_point_(std::move(incoming_point)),
        physical_point_(std::move(physical_point)),
        requested_window_(requested_window),
        required_complete_max_(required_complete_max),
        dimension_(dimension),
        transformation_(std::move(transformation)),
        weights_(std::move(weights)),
        diagnostics_(std::move(diagnostics)),
        residual_window_(residual_window),
        elapsed_ms_(elapsed_ms), basis_owners_(std::move(basis_owners)),
        incoming_owner_(std::move(incoming_owner)) {
    if (basis_sources_.size() != dimension_ ||
        basis_owners_.size() != dimension_ || !incoming_owner_)
      throw std::invalid_argument(
          "retained exact match source ownership differs from its dimension");
  }

  json::object summary() const override {
    json::array basis;
    basis.reserve(basis_sources_.size());
    for (const auto& source : basis_sources_)
      basis.push_back(json::object{
          {"column", source.at("column")},
          {"local", source.at("local")},
          {"checkpoint_identity", source.at("checkpoint_identity")}});

    json::array shifts;
    for (const auto shift : diagnostics_.initial_column_shifts)
      shifts.push_back(shift);
    json::array actions;
    for (const auto& action : diagnostics_.actions)
      actions.push_back(json::object{
          {"leading_rank_before", action.leading_rank_before},
          {"target_column", action.target_column},
          {"relation_support", std::count_if(
               action.null_relation.begin(), action.null_relation.end(),
               [](const Rational& value) { return !value.is_zero(); })}});

    json::array weight_windows;
    for (const auto& weight : weights_)
      weight_windows.push_back(json::object{{"min", weight.min_power()},
                                            {"max", weight.complete_max()}});

    std::size_t transformation_terms = 0;
    for (const auto& row : transformation_)
      for (const auto& entry : row)
        transformation_terms += entry.terms().size();

    return json::object{
        {"match", handle_},
        {"capability", kExactRegularLocalMatchCapability},
        {"retained_state",
         "exact-lattice-transformation-and-laurent-weights"},
        {"native_retained", true},
        {"json_coefficients", 0},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"dimension", dimension_},
        {"basis", std::move(basis)},
        {"incoming", incoming_source_.at("local")},
        {"incoming_checkpoint_identity",
         incoming_source_.at("checkpoint_identity")},
        {"basis_chart", basis_chart_},
        {"incoming_chart", incoming_chart_},
        {"basis_point_exact", basis_point_},
        {"incoming_point_exact", incoming_point_},
        {"physical_match_point_exact", physical_point_},
        {"epsilon", json::object{{"min", requested_window_.min_power},
                                  {"max", requested_window_.complete_max},
                                  {"required_complete_max",
                                   required_complete_max_}}},
        {"weight_windows", std::move(weight_windows)},
        {"transformation_terms", transformation_terms},
        {"initial_column_shifts", std::move(shifts)},
        {"normalized_determinant_valuation",
         diagnostics_.normalized_determinant_valuation},
        {"initial_leading_rank", diagnostics_.initial_leading_rank},
        {"final_leading_rank", diagnostics_.final_leading_rank},
        {"saturation_actions", std::move(actions)},
        {"residual", json::object{{"status", "exact-zero"},
                                   {"scope", "stored-taylor-truncation"},
                                   {"min", residual_window_.min_power},
                                   {"max", residual_window_.complete_max}}},
        {"elapsed_ms", elapsed_ms_}};
  }

  double elapsed_ms() const { return elapsed_ms_; }

  const std::vector<std::shared_ptr<StoredLocalBase>>& basis_owners() const {
    return basis_owners_;
  }
  const std::shared_ptr<StoredLocalBase>& incoming_owner() const {
    return incoming_owner_;
  }
  const FiniteLaurentVector<Rational>& weights() const { return weights_; }

  json::object checkpoint_record() const override {
    json::array basis;
    for (const auto& source : basis_sources_) basis.push_back(source);
    return json::object{
        {"schema", "diffexp2-retained-exact-rational-match-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"basis_sources", std::move(basis)},
        {"incoming_source", incoming_source_},
        {"basis_chart", basis_chart_},
        {"incoming_chart", incoming_chart_},
        {"basis_point_exact", basis_point_},
        {"incoming_point_exact", incoming_point_},
        {"physical_match_point_exact", physical_point_},
        {"epsilon", json::object{{"min", requested_window_.min_power},
                                  {"max", requested_window_.complete_max},
                                  {"required_complete_max",
                                   required_complete_max_}}},
        {"dimension", dimension_},
        {"transformation",
         checkpoint_exact_laurent_matrix_record(transformation_)},
        {"weights", checkpoint_frame_vector_record(weights_)},
        {"saturation",
         checkpoint_saturation_diagnostics_record(diagnostics_)},
        {"residual_window",
         json::object{{"min", residual_window_.min_power},
                      {"max", residual_window_.complete_max}}},
        {"elapsed_ms", elapsed_ms_}};
  }

 private:
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::vector<json::object> basis_sources_;
  json::object incoming_source_;
  std::string basis_chart_;
  std::string incoming_chart_;
  std::string basis_point_;
  std::string incoming_point_;
  std::string physical_point_;
  EpsilonWindow requested_window_;
  std::int32_t required_complete_max_ = 0;
  std::uint32_t dimension_ = 0;
  ExactLaurentMatrix<Rational> transformation_;
  FiniteLaurentVector<Rational> weights_;
  EpsilonLatticeSaturationDiagnostics<Rational> diagnostics_;
  EpsilonWindow residual_window_;
  double elapsed_ms_ = 0.0;
  std::vector<std::shared_ptr<StoredLocalBase>> basis_owners_;
  std::shared_ptr<StoredLocalBase> incoming_owner_;
};

void require_exact_regular_local(const LocalSolution<Rational>& solution,
                                 EpsilonWindow requested_window,
                                 const RealEvaluationPoint& point,
                                 const std::string& label) {
  validate_local_solution(solution, false);
  if (!solution.error.empty())
    throw std::invalid_argument(
        label +
        " carries an error envelope; exact regular local matching does not "
        "silently discard certificates");
  if (solution.sectors.size() != 1)
    throw std::invalid_argument(
        label + " must contain exactly one regular local sector");
  const auto& sector = solution.sectors.front();
  if (sector.a.domain != ExactDomain::Rational ||
      sector.b.domain != ExactDomain::Rational ||
      !(Rational(sector.a.canonical) == Rational(0)) ||
      !(Rational(sector.b.canonical) == Rational(0)) ||
      sector.log_power != 0)
    throw std::invalid_argument(
        label +
        " is not an exact-rational regular (a=0,b=0,log=0) local");
  if (requested_window.complete_max > solution.epsilon.complete_max)
    throw std::invalid_argument(
        label + " does not cover the requested complete epsilon upper edge");
  if (!solution.chart.infinite_radius &&
      !arb_lt(acb_realref(point.modulus.raw()),
              acb_realref(solution.chart.radius.raw())))
    throw std::invalid_argument(
        label + " match point is not provably inside its chart radius");
}

bool same_chart_geometry(const ChartGeometry& left,
                         const ChartGeometry& right) {
  return left.center_exact == right.center_exact &&
         left.scale_exact == right.scale_exact &&
         left.infinite_radius == right.infinite_radius &&
         (left.infinite_radius || acb_equal(left.radius.raw(), right.radius.raw()));
}

Rational physical_match_point(const ChartGeometry& chart,
                              const RealEvaluationPoint& local_point,
                              const std::string& label) {
  try {
    return Rational(chart.center_exact) +
           Rational(chart.scale_exact) *
               Rational(local_point.exact_coordinate);
  } catch (const std::invalid_argument&) {
    throw std::invalid_argument(
        label +
        " requires rational chart center/scale geometry; algebraic geometry "
        "is outside exact-regular-local-match-v1");
  }
}

FiniteLaurentVector<Rational> evaluate_exact_regular_local(
    const LocalSolution<Rational>& solution,
    const RealEvaluationPoint& point, EpsilonWindow window,
    const std::string& label) {
  const Rational t(point.exact_coordinate);
  std::vector<Rational> t_powers(solution.taylor_width(), Rational(1));
  for (std::size_t n = 1; n < t_powers.size(); ++n)
    t_powers[n] = t_powers[n - 1] * t;

  const auto& sector = solution.sectors.front();
  const auto coefficient_at = [&](std::int32_t power,
                                  std::uint32_t component) {
    if (power < solution.epsilon.min_power) return Rational(0);
    const auto epsilon_index = static_cast<std::size_t>(
        static_cast<std::int64_t>(power) - solution.epsilon.min_power);
    Rational coefficient(0);
    for (std::size_t n = 0; n < solution.taylor_width(); ++n)
      coefficient += sector.coefficients[local_detail::sector_index(
                         solution, epsilon_index, n, component)] *
                     t_powers[n];
    return coefficient;
  };
  for (std::int64_t raw_power = solution.epsilon.min_power;
       raw_power < window.min_power; ++raw_power) {
    const auto power = static_cast<std::int32_t>(raw_power);
    for (std::uint32_t component = 0; component < solution.dimension;
         ++component)
      if (!coefficient_at(power, component).is_zero())
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            label +
                " work minimum would discard a nonzero lower epsilon coefficient",
            component, std::nullopt, power);
  }

  FiniteLaurentVector<Rational> value;
  value.reserve(solution.dimension);
  for (std::uint32_t component = 0; component < solution.dimension;
       ++component) {
    std::vector<Rational> coefficients;
    coefficients.reserve(window.width());
    for (std::int64_t power = window.min_power;
         power <= window.complete_max; ++power) {
      coefficients.push_back(coefficient_at(
          static_cast<std::int32_t>(power), component));
    }
    value.emplace_back(window, std::move(coefficients));
  }
  return value;
}

std::shared_ptr<StoredExactRegularMatch> build_exact_regular_match(
    const std::string& match_handle, const json::object& request,
    const std::vector<std::string>& basis_handles,
    const std::vector<std::shared_ptr<StoredLocalBase>>& erased_basis,
    const std::string& incoming_handle,
    const std::shared_ptr<StoredLocalBase>& erased_incoming) {
  const auto started = std::chrono::steady_clock::now();
  const auto& raw_window = as_object(
      request.at("epsilon"), "exact regular match epsilon window");
  EpsilonWindow window{as_i32(raw_window.at("min"), "match epsilon minimum"),
                       as_i32(raw_window.at("max"), "match epsilon maximum")};
  (void)window.width();
  const auto required_complete_max = as_i32(
      raw_window.at("required_complete_max"),
      "required match residual complete maximum");
  if (required_complete_max < window.min_power ||
      required_complete_max > window.complete_max)
    throw std::invalid_argument(
        "required match residual maximum must lie inside the supplied work "
        "epsilon window");

  const auto basis_point = RealEvaluationPoint::rational(required_string(
      as_object(request.at("basis_point"), "basis match point"), "exact"));
  const auto incoming_point = RealEvaluationPoint::rational(required_string(
      as_object(request.at("incoming_point"), "incoming match point"),
      "exact"));
  const auto basis_chart = required_string(request, "basis_chart");
  const auto incoming_chart = required_string(request, "incoming_chart");
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  if (checkpoint_identity.empty())
    throw std::invalid_argument(
        "exact regular match checkpoint identity cannot be empty");

  std::vector<std::shared_ptr<StoredLocal<Rational>>> basis;
  basis.reserve(erased_basis.size());
  for (const auto& local : erased_basis) {
    auto typed = std::dynamic_pointer_cast<StoredLocal<Rational>>(local);
    if (!typed)
      throw std::invalid_argument(
          "exact regular local matching requires rational retained locals");
    basis.push_back(std::move(typed));
  }
  auto incoming =
      std::dynamic_pointer_cast<StoredLocal<Rational>>(erased_incoming);
  if (!incoming)
    throw std::invalid_argument(
        "exact regular local matching requires a rational incoming local");
  if (basis.empty())
    throw std::invalid_argument(
        "exact regular local matching requires a nonempty basis");
  const auto dimension = basis.front()->solution().dimension;
  if (basis.size() != dimension || incoming->solution().dimension != dimension)
    throw std::invalid_argument(
        "exact regular local matching requires d basis columns and a "
        "d-component incoming local");

  const auto& raw_basis_checkpoints = as_array(
      request.at("basis_checkpoint_identities"),
      "basis checkpoint identities");
  if (raw_basis_checkpoints.size() != basis.size())
    throw std::invalid_argument(
        "basis checkpoint identity count differs from the basis dimension");
  std::vector<std::string> basis_checkpoints;
  basis_checkpoints.reserve(basis.size());
  for (std::size_t column = 0; column < basis.size(); ++column) {
    if (!raw_basis_checkpoints[column].is_string())
      throw std::invalid_argument(
          "basis checkpoint identities must be strings");
    const std::string expected(raw_basis_checkpoints[column].as_string());
    if (basis[column]->solution().checkpoint_identity != expected)
      throw std::invalid_argument(
          "basis checkpoint provenance mismatch at column " +
          std::to_string(column));
    basis_checkpoints.push_back(expected);
    if (basis[column]->source_chart() != basis_chart)
      throw std::invalid_argument(
          "basis chart provenance mismatch at column " +
          std::to_string(column));
    if (!same_chart_geometry(basis.front()->solution().chart,
                             basis[column]->solution().chart))
      throw std::invalid_argument(
          "basis locals do not share identical retained chart geometry");
    require_exact_regular_local(
        basis[column]->solution(), window, basis_point,
        "basis local " + basis_handles[column]);
  }
  const auto expected_incoming_checkpoint = required_string(
      request, "incoming_checkpoint_identity");
  if (incoming->solution().checkpoint_identity !=
      expected_incoming_checkpoint)
    throw std::invalid_argument("incoming checkpoint provenance mismatch");
  if (incoming->source_chart() != incoming_chart)
    throw std::invalid_argument("incoming chart provenance mismatch");
  require_exact_regular_local(incoming->solution(), window, incoming_point,
                              "incoming local " + incoming_handle);
  const auto basis_physical_point = physical_match_point(
      basis.front()->solution().chart, basis_point, "basis match point");
  for (std::size_t column = 1; column < basis.size(); ++column)
    if (!(physical_match_point(
              basis[column]->solution().chart, basis_point,
              "basis match point at column " + std::to_string(column)) ==
          basis_physical_point))
      throw std::invalid_argument(
          "basis locals do not name one exact physical match point");
  const auto incoming_physical_point = physical_match_point(
      incoming->solution().chart, incoming_point, "incoming match point");
  if (!(basis_physical_point == incoming_physical_point))
    throw std::invalid_argument(
        "basis and incoming local coordinates do not name the same exact "
        "physical match point");

  FiniteLaurentMatrix<Rational> evaluated_basis(
      dimension, FiniteLaurentVector<Rational>());
  for (auto& row : evaluated_basis) row.reserve(dimension);
  for (std::size_t column = 0; column < basis.size(); ++column) {
    auto value = evaluate_exact_regular_local(
        basis[column]->solution(), basis_point, window,
        "basis local " + basis_handles[column]);
    for (std::uint32_t component = 0; component < dimension; ++component)
      evaluated_basis[component].push_back(std::move(value[component]));
  }
  auto incoming_value = evaluate_exact_regular_local(
      incoming->solution(), incoming_point, window,
      "incoming local " + incoming_handle);

  // Exact zero rows at the lower edge are certified structural zeros.  Trim
  // them once before both saturation and the final residual so honest
  // completeness is not lost merely because a caller supplied a wider
  // declared Laurent minimum than the evaluated entry actually needs.
  for (std::uint32_t row = 0; row < dimension; ++row) {
    incoming_value[row] = matching_detail::canonical_leading_frame(
        incoming_value[row], checkpoint_identity + ":incoming", row);
    for (std::uint32_t column = 0; column < dimension; ++column)
      evaluated_basis[row][column] =
          matching_detail::canonical_leading_frame(
              evaluated_basis[row][column],
              checkpoint_identity + ":basis", row, column);
  }

  auto saturated = saturate_finite_laurent_basis(
      evaluated_basis, checkpoint_identity + ":saturation");
  auto saturated_weights = solve_finite_laurent_system(
      saturated.basis_times_transformation, incoming_value,
      checkpoint_identity + ":solve");
  auto weights = apply_exact_laurent_matrix(
      saturated.transformation, saturated_weights);
  auto reconstructed = apply_finite_laurent_matrix(
      evaluated_basis, weights);

  std::int32_t residual_min = reconstructed.front().min_power();
  std::int32_t residual_max = reconstructed.front().complete_max();
  for (std::uint32_t component = 0; component < dimension; ++component) {
    auto residual = reconstructed[component] - incoming_value[component];
    if (const auto leading = finite_laurent_leading_power(
            residual, checkpoint_identity + ":residual");
        leading.has_value())
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::SaturationFailure,
          checkpoint_identity +
              ": exact regular match residual is nonzero in its complete "
              "window",
          component, std::nullopt, *leading);
    residual_min = std::min(residual_min, residual.min_power());
    residual_max = std::min(residual_max, residual.complete_max());
  }
  (void)EpsilonWindow{residual_min, residual_max}.width();
  if (residual_max < required_complete_max)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InsufficientCompleteWindow,
        checkpoint_identity +
            ": exact regular matching consumed the required complete "
            "epsilon window",
        std::nullopt, std::nullopt, residual_max);

  std::vector<json::object> basis_sources;
  basis_sources.reserve(basis.size());
  for (std::size_t column = 0; column < basis.size(); ++column) {
    json::object entry{{"column", column},
                       {"local", basis_handles[column]},
                       {"chart", basis_chart},
                       {"source_operator_identity",
                        basis[column]->source_operator_identity()},
                       {"checkpoint_identity", basis_checkpoints[column]}};
    entry["analytic_metadata"] =
        basis[column]->exact_analytic_metadata();
    if (basis[column]->column_provenance().has_value())
      entry["column_provenance"] =
          basis[column]->column_provenance()->encode();
    basis_sources.push_back(std::move(entry));
  }
  json::object incoming_source{
      {"local", incoming_handle}, {"chart", incoming_chart},
      {"source_operator_identity", incoming->source_operator_identity()},
      {"checkpoint_identity", expected_incoming_checkpoint},
      {"analytic_metadata", incoming->exact_analytic_metadata()}};
  if (incoming->column_provenance().has_value())
    incoming_source["column_provenance"] =
        incoming->column_provenance()->encode();
  json::array provenance_basis;
  for (const auto& source : basis_sources)
    provenance_basis.push_back(source);
  json::object provenance{
      {"schema", "diffexp2-native-exact-regular-local-match-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"basis", std::move(provenance_basis)},
      {"incoming", incoming_source},
      {"basis_point_exact", basis_point.exact_coordinate},
      {"incoming_point_exact", incoming_point.exact_coordinate},
      {"physical_match_point_exact", basis_physical_point.str()},
      {"epsilon", json::object{{"min", window.min_power},
                                {"max", window.complete_max},
                                {"required_complete_max",
                                 required_complete_max}}}};
  const auto provenance_identity =
      json::serialize(canonical_json_value(provenance));
  const auto elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredExactRegularMatch>(
      match_handle, checkpoint_identity, provenance_identity,
      std::move(basis_sources), std::move(incoming_source), basis_chart,
      incoming_chart,
      basis_point.exact_coordinate, incoming_point.exact_coordinate,
      basis_physical_point.str(), window, required_complete_max, dimension,
      std::move(saturated.transformation), std::move(weights),
      std::move(saturated.diagnostics),
      EpsilonWindow{residual_min, residual_max}, elapsed_ms, erased_basis,
      erased_incoming);
}

const char* acb_match_verdict_name(AcbMatchingResidualVerdict verdict) {
  if (verdict == AcbMatchingResidualVerdict::Pass) return "pass";
  if (verdict == AcbMatchingResidualVerdict::Fail) return "fail";
  return "inconclusive";
}

json::object encode_acb_match_residual_diagnostics(
    const AcbMatchingResidualDiagnostics& diagnostics) {
  std::size_t pass = 0, fail = 0, inconclusive = 0;
  for (const auto& coefficient : diagnostics.coefficients) {
    if (coefficient.verdict == AcbMatchingResidualVerdict::Pass)
      ++pass;
    else if (coefficient.verdict == AcbMatchingResidualVerdict::Fail)
      ++fail;
    else
      ++inconclusive;
  }
  return json::object{
      {"verdict", acb_match_verdict_name(diagnostics.verdict)},
      {"complete_window",
       json::object{{"min", diagnostics.complete_window.min_power},
                    {"max", diagnostics.complete_window.complete_max}}},
      {"required_complete_max", diagnostics.required_complete_max},
      {"complete_through_required",
       diagnostics.complete_through_required},
      {"coefficient_diagnostics", diagnostics.coefficients.size()},
      {"coefficient_verdicts",
       json::object{{"pass", pass}, {"fail", fail},
                    {"inconclusive", inconclusive}}},
      {"detail", diagnostics.detail}};
}

json::object checkpoint_acb_match_residual_record(
    const AcbMatchingResidualDiagnostics& diagnostics) {
  json::array coefficients;
  coefficients.reserve(diagnostics.coefficients.size());
  for (const auto& coefficient : diagnostics.coefficients)
    coefficients.push_back(json::object{
        {"row", coefficient.row},
        {"epsilon_power", coefficient.epsilon_power},
        {"residual_lower_exact", coefficient.residual_lower.dump_exact()},
        {"residual_upper_exact", coefficient.residual_upper.dump_exact()},
        {"scale_lower_exact", coefficient.scale_lower.dump_exact()},
        {"scale_upper_exact", coefficient.scale_upper.dump_exact()},
        {"verdict", acb_match_verdict_name(coefficient.verdict)}});
  return json::object{
      {"verdict", acb_match_verdict_name(diagnostics.verdict)},
      {"complete_window",
       json::object{{"min", diagnostics.complete_window.min_power},
                    {"max", diagnostics.complete_window.complete_max}}},
      {"required_complete_max", diagnostics.required_complete_max},
      {"complete_through_required",
       diagnostics.complete_through_required},
      {"coefficients", std::move(coefficients)},
      {"detail", diagnostics.detail}};
}

class StoredRefinedAcbMatch final : public StoredMatchBase {
 public:
  StoredRefinedAcbMatch(
      std::string handle, std::string checkpoint_identity,
      std::string provenance_identity, std::string exact_lattice_identity,
      std::string exact_lattice_provenance_identity,
      std::string exact_lattice_witness_record,
      std::string saturation_witness_schema,
      std::vector<json::object> basis_sources, json::object incoming_source,
      std::string basis_chart, std::string incoming_chart,
      std::string basis_point, std::string incoming_point,
      std::string physical_point, EpsilonWindow requested_window,
      std::int32_t required_complete_max, std::uint32_t dimension,
      std::string relative_tolerance, std::size_t max_refinement_steps,
      EpsilonLatticeSaturationResult<Rational>&& exact_saturation,
      RefinedAcbLaurentMatch&& refined, double elapsed_ms)
      : StoredMatchBase(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        exact_lattice_identity_(std::move(exact_lattice_identity)),
        exact_lattice_provenance_identity_(
            std::move(exact_lattice_provenance_identity)),
        exact_lattice_witness_record_(
            std::move(exact_lattice_witness_record)),
        saturation_witness_schema_(
            std::move(saturation_witness_schema)),
        basis_sources_(std::move(basis_sources)),
        incoming_source_(std::move(incoming_source)),
        basis_chart_(std::move(basis_chart)),
        incoming_chart_(std::move(incoming_chart)),
        basis_point_(std::move(basis_point)),
        incoming_point_(std::move(incoming_point)),
        physical_point_(std::move(physical_point)),
        requested_window_(requested_window),
        required_complete_max_(required_complete_max),
        dimension_(dimension),
        relative_tolerance_(std::move(relative_tolerance)),
        max_refinement_steps_(max_refinement_steps),
        exact_saturation_(std::move(exact_saturation)),
        refined_(std::move(refined)), elapsed_ms_(elapsed_ms) {}

  json::object summary() const override {
    json::array basis;
    basis.reserve(basis_sources_.size());
    for (const auto& source : basis_sources_) basis.push_back(source);

    json::array history;
    history.reserve(refined_.residual_history.size());
    for (std::size_t iteration = 0;
         iteration < refined_.residual_history.size(); ++iteration) {
      auto encoded = encode_acb_match_residual_diagnostics(
          refined_.residual_history[iteration]);
      encoded["iteration"] = iteration;
      history.push_back(std::move(encoded));
    }
    if (refined_.residual_history.empty())
      throw std::logic_error("retained Acb match has no residual history");

    json::array weight_windows;
    for (const auto& weight : refined_.weights)
      weight_windows.push_back(json::object{{"min", weight.min_power()},
                                            {"max", weight.complete_max()}});
    json::array transformed_weight_windows;
    for (const auto& weight : refined_.transformed_weights)
      transformed_weight_windows.push_back(
          json::object{{"min", weight.min_power()},
                       {"max", weight.complete_max()}});

    json::array shifts;
    for (const auto shift :
         exact_saturation_.diagnostics.initial_column_shifts)
      shifts.push_back(shift);
    json::array actions;
    for (const auto& action : exact_saturation_.diagnostics.actions)
      actions.push_back(json::object{
          {"leading_rank_before", action.leading_rank_before},
          {"target_column", action.target_column},
          {"relation_support", std::count_if(
               action.null_relation.begin(), action.null_relation.end(),
               [](const Rational& value) { return !value.is_zero(); })}});
    std::size_t transformation_terms = 0;
    std::optional<std::int32_t> transformation_minimum;
    for (const auto& row : exact_saturation_.transformation)
      for (const auto& entry : row) {
        transformation_terms += entry.terms().size();
        if (const auto minimum = entry.minimum_power(); minimum.has_value())
          transformation_minimum = !transformation_minimum.has_value()
              ? *minimum
              : std::min(*transformation_minimum, *minimum);
      }

    auto residual = encode_acb_match_residual_diagnostics(
        refined_.residual_history.back());
    residual["scope"] = "stored-taylor-truncation";
    residual["history"] = std::move(history);
    return json::object{
        {"match", handle_},
        {"capability", kRefinedAcbLocalMatchCapability},
        {"retained_state",
         "exact-lattice-transformation-acb-weights-and-residual"},
        {"native_retained", true},
        {"json_coefficients", 0},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"dimension", dimension_},
        {"basis", std::move(basis)},
        {"incoming", incoming_source_},
        {"basis_chart", basis_chart_},
        {"incoming_chart", incoming_chart_},
        {"basis_point_exact", basis_point_},
        {"incoming_point_exact", incoming_point_},
        {"physical_match_point_exact", physical_point_},
        {"epsilon",
         json::object{{"min", requested_window_.min_power},
                      {"max", requested_window_.complete_max},
                      {"required_complete_max", required_complete_max_}}},
        {"exact_lattice",
         json::object{
             {"schema", saturation_witness_schema_},
             {"identity", exact_lattice_identity_},
             {"provenance_identity",
              exact_lattice_provenance_identity_},
             {"canonical_witness_retained", true},
             {"canonical_witness_bytes",
              exact_lattice_witness_record_.size()},
             {"transformation_terms", transformation_terms},
             {"transformation_min_power",
              transformation_minimum.has_value()
                  ? json::value(*transformation_minimum)
                  : json::value(nullptr)},
             {"initial_column_shifts", std::move(shifts)},
             {"normalized_determinant_valuation",
              exact_saturation_.diagnostics
                  .normalized_determinant_valuation},
             {"initial_leading_rank",
              exact_saturation_.diagnostics.initial_leading_rank},
             {"final_leading_rank",
              exact_saturation_.diagnostics.final_leading_rank},
             {"saturation_actions", std::move(actions)}}},
        {"refinement",
         json::object{{"relative_tolerance", relative_tolerance_},
                      {"max_steps", max_refinement_steps_},
                      {"steps", refined_.refinement_steps},
                      {"factorizations", 1}}},
        {"weight_windows", std::move(weight_windows)},
        {"transformed_weight_windows",
         std::move(transformed_weight_windows)},
        {"residual", std::move(residual)},
        {"elapsed_ms", elapsed_ms_}};
  }

  const FiniteLaurentVector<ComplexBall>& weights() const {
    return refined_.weights;
  }

  bool certified_for_materialization() const {
    return !refined_.residual_history.empty() &&
        refined_.residual_history.back().verdict ==
            AcbMatchingResidualVerdict::Pass &&
        refined_.residual_history.back().complete_through_required;
  }

  double elapsed_ms() const { return elapsed_ms_; }

  json::object checkpoint_record() const override {
    if (refined_.residual_history.empty())
      throw std::logic_error(
          "cannot checkpoint an Acb match without residual history");
    json::array basis;
    basis.reserve(basis_sources_.size());
    for (const auto& source : basis_sources_) basis.push_back(source);
    json::array history;
    history.reserve(refined_.residual_history.size());
    for (const auto& residual : refined_.residual_history)
      history.push_back(checkpoint_acb_match_residual_record(residual));
    return json::object{
        {"schema", "diffexp2-retained-acb-match-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"exact_lattice_identity", exact_lattice_identity_},
        {"exact_lattice_provenance_identity",
         exact_lattice_provenance_identity_},
        {"exact_lattice_canonical_witness",
         exact_lattice_witness_record_},
        {"basis_sources", std::move(basis)},
        {"incoming_source", incoming_source_},
        {"basis_chart", basis_chart_},
        {"incoming_chart", incoming_chart_},
        {"basis_point_exact", basis_point_},
        {"incoming_point_exact", incoming_point_},
        {"physical_match_point_exact", physical_point_},
        {"epsilon",
         json::object{{"min", requested_window_.min_power},
                      {"max", requested_window_.complete_max},
                      {"required_complete_max", required_complete_max_}}},
        {"dimension", dimension_},
        {"relative_tolerance", relative_tolerance_},
        {"max_refinement_steps", max_refinement_steps_},
        {"refined",
         json::object{
             {"transformed_weights",
              checkpoint_frame_vector_record(refined_.transformed_weights)},
             {"weights", checkpoint_frame_vector_record(refined_.weights)},
             {"residual", checkpoint_frame_vector_record(refined_.residual)},
             {"residual_history", std::move(history)},
             {"refinement_steps", refined_.refinement_steps}}},
        {"elapsed_ms", elapsed_ms_}};
  }

 private:
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::string exact_lattice_identity_;
  std::string exact_lattice_provenance_identity_;
  std::string exact_lattice_witness_record_;
  std::string saturation_witness_schema_;
  std::vector<json::object> basis_sources_;
  json::object incoming_source_;
  std::string basis_chart_;
  std::string incoming_chart_;
  std::string basis_point_;
  std::string incoming_point_;
  std::string physical_point_;
  EpsilonWindow requested_window_;
  std::int32_t required_complete_max_ = 0;
  std::uint32_t dimension_ = 0;
  std::string relative_tolerance_;
  std::size_t max_refinement_steps_ = 0;
  EpsilonLatticeSaturationResult<Rational> exact_saturation_;
  RefinedAcbLaurentMatch refined_;
  double elapsed_ms_ = 0.0;
};

struct ParsedExactEvaluatedLattice {
  std::string witness_schema;
  std::string identity;
  std::string canonical_witness;
  EpsilonLatticeSaturationResult<Rational> saturation;
};

ParsedExactEvaluatedLattice parse_exact_evaluated_lattice(
    const json::value& raw, std::uint32_t dimension, EpsilonWindow window,
    const std::string& context) {
  const auto& object = as_object(raw, "exact evaluated epsilon lattice");
  if (object.size() != 3 || object.if_contains("schema") == nullptr ||
      object.if_contains("identity") == nullptr ||
      object.if_contains("evaluated_basis") == nullptr)
    throw std::invalid_argument(
        "exact_lattice accepts exactly schema, identity, and evaluated_basis");
  if (required_string(object, "schema") != kExactEvaluatedLatticeSchema)
    throw std::invalid_argument("unsupported exact evaluated lattice schema");
  const auto identity = required_string(object, "identity");
  if (identity.empty())
    throw std::invalid_argument("exact evaluated lattice identity is empty");
  const auto& rows = as_array(object.at("evaluated_basis"),
                              "exact evaluated lattice basis");
  if (rows.size() != dimension)
    throw std::invalid_argument(
        "exact evaluated lattice must have one row per local component");

  FiniteLaurentMatrix<Rational> basis;
  basis.reserve(dimension);
  json::array canonical_rows;
  canonical_rows.reserve(dimension);
  for (std::uint32_t row = 0; row < dimension; ++row) {
    const auto& columns = as_array(rows[row],
                                   "exact evaluated lattice basis row");
    if (columns.size() != dimension)
      throw std::invalid_argument(
          "exact evaluated lattice basis must be square");
    FiniteLaurentVector<Rational> parsed_row;
    parsed_row.reserve(dimension);
    json::array canonical_columns;
    canonical_columns.reserve(dimension);
    for (std::uint32_t column = 0; column < dimension; ++column) {
      const auto& frame = as_object(columns[column],
                                    "exact evaluated lattice frame");
      if (frame.size() != 3 || frame.if_contains("min") == nullptr ||
          frame.if_contains("max") == nullptr ||
          frame.if_contains("coefficients") == nullptr)
        throw std::invalid_argument(
            "exact evaluated lattice frame accepts exactly min, max, and coefficients");
      const EpsilonWindow frame_window{
          as_i32(frame.at("min"), "exact lattice epsilon minimum"),
          as_i32(frame.at("max"), "exact lattice epsilon maximum")};
      (void)frame_window.width();
      if (frame_window.min_power != window.min_power ||
          frame_window.complete_max != window.complete_max)
        throw std::invalid_argument(
            "every exact evaluated lattice frame must equal the matching work window");
      const auto& raw_coefficients = as_array(
          frame.at("coefficients"), "exact lattice coefficients");
      if (raw_coefficients.size() != frame_window.width())
        throw std::invalid_argument(
            "exact lattice coefficient count differs from its epsilon window");
      std::vector<Rational> coefficients;
      coefficients.reserve(raw_coefficients.size());
      json::array canonical_coefficients;
      canonical_coefficients.reserve(raw_coefficients.size());
      for (const auto& raw_coefficient : raw_coefficients) {
        auto coefficient = parse_scalar<Rational>(raw_coefficient);
        canonical_coefficients.emplace_back(coefficient.str());
        coefficients.push_back(std::move(coefficient));
      }
      parsed_row.emplace_back(frame_window, std::move(coefficients));
      canonical_columns.push_back(json::object{
          {"min", frame_window.min_power},
          {"max", frame_window.complete_max},
          {"coefficients", std::move(canonical_coefficients)}});
    }
    basis.push_back(std::move(parsed_row));
    canonical_rows.push_back(std::move(canonical_columns));
  }
  json::object canonical{
      {"schema", kExactEvaluatedLatticeSchema}, {"identity", identity},
      {"evaluated_basis", std::move(canonical_rows)}};
  auto saturation = saturate_finite_laurent_basis(
      basis, context + ":exact-lattice-saturation");
  return {kExactEvaluatedLatticeSchema, identity,
          json::serialize(canonical), std::move(saturation)};
}

EpsilonLatticeSaturationResult<Rational> unit_rational_saturation(
    std::uint32_t dimension, EpsilonWindow window,
    const std::string& context) {
  if (dimension == 0 || window.min_power > 0 || window.complete_max < 0)
    throw std::invalid_argument(
        context + ": the unit saturation requires epsilon^0 in a nonempty square window");
  FiniteLaurentMatrix<Rational> identity(
      dimension, FiniteLaurentVector<Rational>());
  for (std::uint32_t row = 0; row < dimension; ++row) {
    identity[row].reserve(dimension);
    for (std::uint32_t column = 0; column < dimension; ++column) {
      std::vector<Rational> coefficients(
          window.width(), Rational(0));
      if (row == column)
        coefficients[static_cast<std::size_t>(
            -static_cast<std::int64_t>(window.min_power))] =
            Rational(1);
      identity[row].emplace_back(window, std::move(coefficients));
    }
  }
  auto saturation = saturate_finite_laurent_basis(
      identity, context + ":unit-saturation");
  if (saturation.diagnostics.initial_leading_rank != dimension ||
      saturation.diagnostics.final_leading_rank != dimension ||
      saturation.diagnostics.normalized_determinant_valuation != 0 ||
      !saturation.diagnostics.actions.empty() ||
      std::any_of(
          saturation.diagnostics.initial_column_shifts.begin(),
          saturation.diagnostics.initial_column_shifts.end(),
          [](std::int32_t shift) { return shift != 0; }))
    throw std::logic_error(
        context + ": internally constructed unit saturation is not the identity transformation");
  return saturation;
}

json::object validate_native_unit_saturation_request(
    const json::value& raw, const std::string& context) {
  auto request = as_object(raw, "native Acb unit-saturation request");
  require_exact_keys(
      request,
      {"schema", "tile_plan", "tile_plan_checkpoint_identity",
       "tile_plan_provenance_identity", "arm", "match"},
      "native Acb unit-saturation request");
  if (required_string(request, "schema") !=
      kNativeUnitSaturationRequestSchema)
    throw std::invalid_argument(
        context + ": unsupported native unit-saturation request schema");
  if (required_string(request, "tile_plan").empty() ||
      required_string(request, "tile_plan_checkpoint_identity").empty() ||
      required_string(request, "tile_plan_provenance_identity").empty())
    throw std::invalid_argument(
        context + ": native unit-saturation request lost its plan binding");
  const auto arm = required_string(request, "arm");
  if (arm != "lower" && arm != "upper")
    throw std::invalid_argument(
        context + ": native unit-saturation request has an unknown arm");
  (void)as_u64(request.at("match"),
               "native unit-saturation match index");
  return request;
}

void require_ordinary_regular_basis_for_unit_saturation(
    const std::vector<std::shared_ptr<StoredLocalBase>>& basis) {
  if (basis.empty())
    throw std::invalid_argument(
        "native Acb unit-leading certification requires a nonempty basis");
  const auto dimension = as_u32(
      basis.front()->summary().at("dimension"),
      "native unit-leading basis dimension");
  if (basis.size() != dimension)
    throw std::invalid_argument(
        "native Acb unit-leading certification requires a square basis");
  for (const auto& column : basis) {
    const auto metadata = column->exact_analytic_metadata();
    const auto& sectors = as_array(
        metadata.at("sectors"), "native unit-leading sectors");
    if (sectors.size() != 1)
      throw std::domain_error(
          "native Acb unit-leading certification requires one ordinary sector per basis column");
    const auto& sector = as_object(
        sectors.front(), "native unit-leading sector");
    const auto a = parse_checkpoint_exact_descriptor(
        sector.at("a"), "native unit-leading a tag");
    const auto b = parse_checkpoint_exact_descriptor(
        sector.at("b"), "native unit-leading b tag");
    if (a.domain != ExactDomain::Rational ||
        b.domain != ExactDomain::Rational ||
        !(Rational(a.canonical) == Rational(0)) ||
        !(Rational(b.canonical) == Rational(0)) ||
        as_u32(sector.at("log_power"),
               "native unit-leading log power") != 0)
      throw std::domain_error(
          "native Acb unit-leading certification requires exact a=b=0, log_power=0 basis tags");
  }
}

ParsedExactEvaluatedLattice certify_native_unit_saturation(
    const json::value& raw_request,
    const FiniteLaurentMatrix<ComplexBall>& evaluated_basis,
    const std::vector<std::shared_ptr<StoredLocalBase>>& retained_basis,
    const std::vector<json::object>& basis_sources,
    const std::string& basis_point, const std::string& physical_point,
    EpsilonWindow window, const std::string& context) {
  require_ordinary_regular_basis_for_unit_saturation(retained_basis);
  const auto dimension = retained_basis.size();
  if (evaluated_basis.size() != dimension ||
      basis_sources.size() != dimension || window.min_power > 0 ||
      window.complete_max < 0)
    throw std::domain_error(
        context + ": native unit-leading certification requires a square actual basis complete through epsilon^0");
  for (std::size_t row = 0; row < dimension; ++row) {
    if (evaluated_basis[row].size() != dimension)
      throw std::domain_error(
          context + ": native unit-leading certification received a nonsquare actual basis");
    for (std::size_t column = 0; column < dimension; ++column) {
      const auto& frame = evaluated_basis[row][column];
      if (frame.complete_max() < 0)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            context + ": actual Acb basis is incomplete through epsilon^0",
            row, column, frame.complete_max());
      for (std::int64_t power = window.min_power; power < 0; ++power)
        if (!frame.coefficient(static_cast<std::int32_t>(power)).is_zero())
          throw MatchingArithmeticError(
              MatchingArithmeticErrorCode::InvalidSaturationLattice,
              context + ": actual Acb basis has a nonzero or zero-ambiguous negative epsilon coefficient",
              row, column, static_cast<std::int32_t>(power));
    }
  }
  const auto leading_rank =
      matching_detail::certify_full_rank_by_nonzero_pivots(
          matching_detail::epsilon_zero_matrix(
              evaluated_basis, context + ":actual-leading-frame"),
          context + ":actual-leading-rank");
  if (leading_rank != dimension)
    throw std::logic_error(
        context + ": full-rank proof returned the wrong dimension");

  auto native_request = validate_native_unit_saturation_request(
      raw_request, context);
  json::array proof_basis;
  proof_basis.reserve(basis_sources.size());
  for (const auto& source : basis_sources) proof_basis.push_back(source);
  json::object proof_without_identity{
      {"schema", kNativeUnitSaturationProofSchema},
      {"native_request", std::move(native_request)},
      {"coefficient_domain", "acb"},
      {"basis", std::move(proof_basis)},
      {"basis_point_exact", basis_point},
      {"physical_match_point_exact", physical_point},
      {"epsilon", json::object{{"min", window.min_power},
                                {"max", window.complete_max}}},
      {"negative_epsilon_coefficients", "exact-singleton-zero"},
      {"leading_power", 0},
      {"leading_rank", dimension},
      {"leading_rank_certificate",
       "full-pivot-acb-pivots-exclude-zero"},
      {"transformation", "identity"}};
  const auto identity = json::serialize(
      canonical_json_value(proof_without_identity));
  auto proof = proof_without_identity;
  proof["identity"] = identity;
  return {kNativeUnitSaturationProofSchema, identity,
          json::serialize(canonical_json_value(proof)),
          unit_rational_saturation(
              static_cast<std::uint32_t>(dimension), window, context)};
}

ParsedExactEvaluatedLattice parse_native_unit_saturation_proof(
    const json::value& raw, std::uint32_t dimension, EpsilonWindow window,
    const std::vector<json::object>& expected_basis_sources,
    const std::string& expected_basis_point,
    const std::string& expected_physical_point,
    const std::string& context) {
  const auto& proof = as_object(raw, "native Acb unit-saturation proof");
  require_exact_keys(
      proof,
      {"schema", "identity", "native_request", "coefficient_domain",
       "basis", "basis_point_exact", "physical_match_point_exact",
       "epsilon", "negative_epsilon_coefficients", "leading_power",
       "leading_rank", "leading_rank_certificate", "transformation"},
      "native Acb unit-saturation proof");
  if (required_string(proof, "schema") !=
          kNativeUnitSaturationProofSchema ||
      required_string(proof, "coefficient_domain") != "acb" ||
      required_string(proof, "negative_epsilon_coefficients") !=
          "exact-singleton-zero" ||
      as_i32(proof.at("leading_power"),
             "native unit-leading power") != 0 ||
      as_u32(proof.at("leading_rank"),
             "native unit-leading rank") != dimension ||
      required_string(proof, "leading_rank_certificate") !=
          "full-pivot-acb-pivots-exclude-zero" ||
      required_string(proof, "transformation") != "identity")
    throw std::invalid_argument(
        context + ": native unit-saturation proof facts are inconsistent");
  (void)validate_native_unit_saturation_request(
      proof.at("native_request"), context);
  const auto& epsilon = as_object(
      proof.at("epsilon"), "native unit-saturation proof epsilon");
  require_exact_keys(epsilon, {"min", "max"},
                     "native unit-saturation proof epsilon");
  if (as_i32(epsilon.at("min"), "native proof epsilon minimum") !=
          window.min_power ||
      as_i32(epsilon.at("max"), "native proof epsilon maximum") !=
          window.complete_max ||
      required_string(proof, "basis_point_exact") !=
          expected_basis_point ||
      required_string(proof, "physical_match_point_exact") !=
          expected_physical_point)
    throw std::invalid_argument(
        context + ": native unit-saturation proof changed its point or epsilon binding");
  json::array expected_basis;
  expected_basis.reserve(expected_basis_sources.size());
  for (const auto& source : expected_basis_sources)
    expected_basis.push_back(source);
  if (json::serialize(canonical_json_value(proof.at("basis"))) !=
      json::serialize(canonical_json_value(expected_basis)))
    throw std::invalid_argument(
        context + ": native unit-saturation proof changed its basis/checkpoint binding");
  auto identity_input = proof;
  const auto identity = required_string(proof, "identity");
  identity_input.erase("identity");
  if (identity.empty() ||
      json::serialize(canonical_json_value(identity_input)) != identity)
    throw std::invalid_argument(
        context + ": native unit-saturation proof identity is inconsistent");
  return {kNativeUnitSaturationProofSchema, identity,
          json::serialize(canonical_json_value(proof)),
          unit_rational_saturation(dimension, window, context)};
}

bool is_supported_acb_singular_scc_column_capability(
    const std::string& capability) {
  return capability == kAcbSingularScalarSCCColumnCapability ||
      capability == kAcbSingularJordanSCCColumnCapability;
}

json::object validate_native_singular_scc_saturation_request(
    const json::value& raw, const std::string& context,
    const std::optional<std::string>& expected_session_configuration =
        std::nullopt,
    const std::optional<json::object>& expected_request = std::nullopt) {
  auto request = as_object(
      raw, "native Acb singular-SCC valuation-zero request");
  require_exact_keys(
      request,
      {"schema", "session_configuration_identity", "tile_plan",
       "tile_plan_checkpoint_identity", "tile_plan_provenance_identity",
       "arm", "match", "match_checkpoint_identity", "receiving_scc",
       "receiving_scc_exact_identity", "receiving_execution_capability",
       "receiving_basis_point_exact", "physical_match_point_exact",
       "receiving_rim"},
      "native Acb singular-SCC valuation-zero request");
  if (required_string(request, "schema") !=
      kNativeSingularSCCSaturationRequestSchema)
    throw std::invalid_argument(
        context + ": unsupported native singular-SCC saturation request schema");
  for (const auto* key :
       {"session_configuration_identity", "tile_plan",
        "tile_plan_checkpoint_identity",
        "tile_plan_provenance_identity", "match_checkpoint_identity",
        "receiving_scc", "receiving_scc_exact_identity",
        "receiving_basis_point_exact", "physical_match_point_exact"})
    if (required_string(request, key).empty())
      throw std::invalid_argument(
          context + ": native singular-SCC saturation request lost its " +
          key + " binding");
  const auto arm = required_string(request, "arm");
  if (arm != "lower" && arm != "upper")
    throw std::invalid_argument(
        context + ": native singular-SCC saturation request has an unknown arm");
  (void)as_u64(request.at("match"),
               "native singular-SCC saturation match index");
  const auto capability = required_string(
      request, "receiving_execution_capability");
  if (!is_supported_acb_singular_scc_column_capability(capability))
    throw std::domain_error(
        context + ": receiving SCC lacks a supported affine-Jordan Acb column capability");
  if (!request.at("receiving_rim").is_null()) {
    const auto rim = as_i32(
        request.at("receiving_rim"), "native singular-SCC receiving rim");
    if (rim != -1 && rim != 1)
      throw std::invalid_argument(
          context + ": native singular-SCC receiving rim must be +1, -1, or null");
  }
  if (expected_session_configuration.has_value() &&
      required_string(request, "session_configuration_identity") !=
          *expected_session_configuration)
    throw std::invalid_argument(
        context + ": native singular-SCC saturation request changed its stable session-configuration binding");
  if (expected_request.has_value() &&
      json::serialize(canonical_json_value(request)) !=
          json::serialize(canonical_json_value(*expected_request)))
    throw std::invalid_argument(
        context + ": native singular-SCC saturation request changed its retained plan/SCC binding");
  return request;
}

void validate_singular_scc_basis_sources(
    const json::object& native_request,
    const std::vector<json::object>& basis_sources,
    std::uint32_t dimension, const std::string& context) {
  if (dimension == 0 || basis_sources.size() != dimension)
    throw std::domain_error(
        context + ": singular-SCC valuation-zero proof requires one complete square basis");
  const auto expected_scc = required_string(native_request, "receiving_scc");
  const auto expected_identity = required_string(
      native_request, "receiving_scc_exact_identity");
  const auto capability = required_string(
      native_request, "receiving_execution_capability");
  const bool scalar = capability == kAcbSingularScalarSCCColumnCapability;
  const Rational basis_point(required_string(
      native_request, "receiving_basis_point_exact"));
  if (basis_point.is_zero())
    throw std::invalid_argument(
        context + ": singular-SCC matching cannot certify a chart-center evaluation");
  const json::value expected_effective_rim = basis_point.sign() < 0
      ? native_request.at("receiving_rim") : json::value(nullptr);
  const auto* expected_column_schema = scalar
      ? "diffexp2-native-scc-acb-regular-singular-scalar-column-v1"
      : "diffexp2-native-scc-acb-regular-singular-jordan-column-v1";
  std::vector<std::uint8_t> seen(dimension, 0);
  for (std::size_t column = 0; column < basis_sources.size(); ++column) {
    const auto& source = basis_sources[column];
    if (as_u64(source.at("column"), "singular-SCC proof basis column") !=
        column)
      throw std::invalid_argument(
          context + ": singular-SCC basis sources are not in receiving column order");
    if (required_string(source, "chart") != expected_scc ||
        required_string(source, "source_operator_identity") !=
            expected_identity)
      throw std::invalid_argument(
          context + ": singular-SCC basis source is not owned by the receiving SCC");
    if (source.at("requested_imaginary_sign") !=
            native_request.at("receiving_rim") ||
        source.at("effective_imaginary_sign") != expected_effective_rim)
      throw std::invalid_argument(
          context + ": singular-SCC basis source changed its point-dependent requested/effective plan-selected rim");
    const auto* raw_provenance = source.if_contains("column_provenance");
    if (raw_provenance == nullptr)
      throw std::domain_error(
          context + ": singular-SCC basis column lacks certified SCC provenance");
    const auto& provenance = as_object(
        *raw_provenance, "singular-SCC basis column provenance");
    require_exact_keys(
        provenance,
        {"scc", "scc_exact_identity", "seed_block", "basis_index",
         "exact_column_identity"},
        "singular-SCC basis column provenance");
    const auto basis_index = as_u32(
        provenance.at("basis_index"), "singular-SCC provenance basis index");
    if (required_string(provenance, "scc") != expected_scc ||
        required_string(provenance, "scc_exact_identity") !=
            expected_identity ||
        basis_index >= dimension || seen[basis_index] != 0)
      throw std::domain_error(
          context + ": singular-SCC basis provenance is incomplete, duplicated, or belongs to another SCC");
    seen[basis_index] = 1;

    const auto exact_column_record = required_string(
        provenance, "exact_column_identity");
    if (exact_column_record.empty())
      throw std::invalid_argument(
          context + ": singular-SCC basis column lost its exact identity");
    const auto parsed_column = json::parse(exact_column_record);
    const auto& exact_column = as_object(
        parsed_column, "singular-SCC exact column identity");
    if (json::serialize(canonical_json_value(parsed_column)) !=
        exact_column_record)
      throw std::invalid_argument(
          context + ": singular-SCC exact column identity is not canonically encoded");
    if (scalar)
      require_exact_keys(
          exact_column,
          {"schema", "scc_exact_identity", "basis_index", "seed",
           "targets", "pseudo_compensation"},
          "singular-SCC scalar exact column identity");
    else
      require_exact_keys(
          exact_column,
          {"schema", "scc_exact_identity", "basis_index", "seed",
           "targets", "pseudo_compensation", "seed_local_component"},
          "singular-SCC Jordan exact column identity");
    if (required_string(exact_column, "schema") != expected_column_schema ||
        required_string(exact_column, "scc_exact_identity") !=
            expected_identity ||
        as_u32(exact_column.at("basis_index"),
               "singular-SCC exact column basis index") != basis_index ||
        required_string(exact_column, "pseudo_compensation") != "none")
      throw std::domain_error(
          context + ": singular-SCC basis column is not a supported no-pseudo affine-Jordan Acb column");
    const auto& seed = as_object(
        exact_column.at("seed"), "singular-SCC exact column seed");
    if (as_u32(seed.at("block"), "singular-SCC exact column seed block") !=
        as_u32(provenance.at("seed_block"),
               "singular-SCC provenance seed block"))
      throw std::invalid_argument(
          context + ": singular-SCC exact column changed its seed-block binding");
    (void)as_array(exact_column.at("targets"),
                   "singular-SCC exact column targets");
  }
  if (std::any_of(seen.begin(), seen.end(),
                  [](std::uint8_t value) { return value == 0; }))
    throw std::domain_error(
        context + ": singular-SCC basis provenance does not cover every canonical column");
}

ParsedExactEvaluatedLattice certify_native_singular_scc_saturation(
    const json::value& raw_request,
    const FiniteLaurentMatrix<ComplexBall>& evaluated_basis,
    const std::vector<std::shared_ptr<StoredLocalBase>>& retained_basis,
    const std::vector<json::object>& basis_sources,
    const std::string& basis_point, const std::string& physical_point,
    EpsilonWindow window,
    const std::string& expected_session_configuration,
    const json::object& expected_native_request,
    const std::string& expected_checkpoint_identity,
    const std::string& context) {
  auto native_request = validate_native_singular_scc_saturation_request(
      raw_request, context, expected_session_configuration,
      expected_native_request);
  if (required_string(native_request, "match_checkpoint_identity") !=
          expected_checkpoint_identity ||
      required_string(native_request, "receiving_basis_point_exact") !=
          basis_point ||
      required_string(native_request, "physical_match_point_exact") !=
          physical_point)
    throw std::invalid_argument(
        context + ": native singular-SCC proof request changed its checkpoint or point binding");
  const auto dimension = retained_basis.size();
  if (dimension == 0 || evaluated_basis.size() != dimension ||
      basis_sources.size() != dimension || window.min_power > 0 ||
      window.complete_max < 0)
    throw std::domain_error(
        context + ": singular-SCC valuation-zero certification requires a square actual basis complete through epsilon^0");
  validate_singular_scc_basis_sources(
      native_request, basis_sources, static_cast<std::uint32_t>(dimension),
      context);
  for (std::size_t column = 0; column < dimension; ++column) {
    const auto& provenance = retained_basis[column]->column_provenance();
    if (!provenance.has_value() ||
        json::serialize(canonical_json_value(provenance->encode())) !=
            json::serialize(canonical_json_value(
                basis_sources[column].at("column_provenance"))))
      throw std::invalid_argument(
          context + ": singular-SCC proof source disagrees with its retained column owner");
  }
  for (std::size_t row = 0; row < dimension; ++row) {
    if (evaluated_basis[row].size() != dimension)
      throw std::domain_error(
          context + ": singular-SCC valuation-zero certification received a nonsquare actual basis");
    for (std::size_t column = 0; column < dimension; ++column) {
      const auto& frame = evaluated_basis[row][column];
      if (frame.complete_max() < 0)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            context + ": actual singular-SCC Acb basis is incomplete through epsilon^0",
            row, column, frame.complete_max());
      for (std::int64_t power = window.min_power; power < 0; ++power)
        if (!frame.coefficient(static_cast<std::int32_t>(power)).is_zero())
          throw MatchingArithmeticError(
              MatchingArithmeticErrorCode::InvalidSaturationLattice,
              context + ": actual singular-SCC Acb basis has a nonzero or zero-ambiguous negative epsilon coefficient",
              row, column, static_cast<std::int32_t>(power));
    }
  }
  const auto leading_rank =
      matching_detail::certify_full_rank_by_nonzero_pivots(
          matching_detail::epsilon_zero_matrix(
              evaluated_basis, context + ":actual-leading-frame"),
          context + ":actual-leading-rank");
  if (leading_rank != dimension)
    throw std::logic_error(
        context + ": singular-SCC full-rank proof returned the wrong dimension");

  json::array proof_basis;
  proof_basis.reserve(basis_sources.size());
  for (const auto& source : basis_sources) proof_basis.push_back(source);
  json::object proof_without_identity{
      {"schema", kNativeSingularSCCSaturationProofSchema},
      {"native_request", std::move(native_request)},
      {"coefficient_domain", "acb"},
      {"basis", std::move(proof_basis)},
      {"basis_point_exact", basis_point},
      {"physical_match_point_exact", physical_point},
      {"epsilon", json::object{{"min", window.min_power},
                                {"max", window.complete_max}}},
      {"negative_epsilon_coefficients", "exact-singleton-zero"},
      {"leading_power", 0},
      {"leading_rank", dimension},
      {"leading_rank_certificate", "full-pivot-acb-pivots-exclude-zero"},
      {"column_provenance_certificate",
       "complete-one-receiving-scc-composite-affine-jordan-acb-no-pseudo"},
      {"determinant_valuation", 0},
      {"transformation", "identity"}};
  const auto identity = json::serialize(
      canonical_json_value(proof_without_identity));
  auto proof = proof_without_identity;
  proof["identity"] = identity;
  return {kNativeSingularSCCSaturationProofSchema, identity,
          json::serialize(canonical_json_value(proof)),
          unit_rational_saturation(
              static_cast<std::uint32_t>(dimension), window, context)};
}

ParsedExactEvaluatedLattice parse_native_singular_scc_saturation_proof(
    const json::value& raw, std::uint32_t dimension, EpsilonWindow window,
    const std::vector<json::object>& expected_basis_sources,
    const std::string& expected_basis_point,
    const std::string& expected_physical_point,
    const std::optional<std::string>& expected_session_configuration,
    const std::optional<json::object>& expected_native_request,
    const std::string& context) {
  const auto& proof = as_object(
      raw, "native Acb singular-SCC valuation-zero proof");
  require_exact_keys(
      proof,
      {"schema", "identity", "native_request", "coefficient_domain",
       "basis", "basis_point_exact", "physical_match_point_exact",
       "epsilon", "negative_epsilon_coefficients", "leading_power",
       "leading_rank", "leading_rank_certificate",
       "column_provenance_certificate", "determinant_valuation",
       "transformation"},
      "native Acb singular-SCC valuation-zero proof");
  if (required_string(proof, "schema") !=
          kNativeSingularSCCSaturationProofSchema ||
      required_string(proof, "coefficient_domain") != "acb" ||
      required_string(proof, "negative_epsilon_coefficients") !=
          "exact-singleton-zero" ||
      as_i32(proof.at("leading_power"),
             "singular-SCC proof leading power") != 0 ||
      as_u32(proof.at("leading_rank"),
             "singular-SCC proof leading rank") != dimension ||
      required_string(proof, "leading_rank_certificate") !=
          "full-pivot-acb-pivots-exclude-zero" ||
      required_string(proof, "column_provenance_certificate") !=
          "complete-one-receiving-scc-composite-affine-jordan-acb-no-pseudo" ||
      as_i32(proof.at("determinant_valuation"),
             "singular-SCC determinant valuation") != 0 ||
      required_string(proof, "transformation") != "identity")
    throw std::invalid_argument(
        context + ": native singular-SCC saturation proof facts are inconsistent");
  auto native_request = validate_native_singular_scc_saturation_request(
      proof.at("native_request"), context,
      expected_session_configuration,
      expected_native_request);
  const auto& epsilon = as_object(
      proof.at("epsilon"), "singular-SCC proof epsilon window");
  require_exact_keys(epsilon, {"min", "max"},
                     "singular-SCC proof epsilon window");
  if (as_i32(epsilon.at("min"), "singular-SCC proof epsilon minimum") !=
          window.min_power ||
      as_i32(epsilon.at("max"), "singular-SCC proof epsilon maximum") !=
          window.complete_max ||
      required_string(proof, "basis_point_exact") !=
          expected_basis_point ||
      required_string(proof, "physical_match_point_exact") !=
          expected_physical_point ||
      required_string(native_request, "receiving_basis_point_exact") !=
          expected_basis_point ||
      required_string(native_request, "physical_match_point_exact") !=
          expected_physical_point)
    throw std::invalid_argument(
        context + ": native singular-SCC proof changed its point or epsilon binding");
  json::array expected_basis;
  expected_basis.reserve(expected_basis_sources.size());
  for (const auto& source : expected_basis_sources)
    expected_basis.push_back(source);
  if (json::serialize(canonical_json_value(proof.at("basis"))) !=
      json::serialize(canonical_json_value(expected_basis)))
    throw std::invalid_argument(
        context + ": native singular-SCC proof changed its basis/checkpoint binding");
  validate_singular_scc_basis_sources(
      native_request, expected_basis_sources, dimension, context);
  auto identity_input = proof;
  const auto identity = required_string(proof, "identity");
  identity_input.erase("identity");
  if (identity.empty() ||
      json::serialize(canonical_json_value(identity_input)) != identity)
    throw std::invalid_argument(
        context + ": native singular-SCC saturation proof identity is inconsistent");
  return {kNativeSingularSCCSaturationProofSchema, identity,
          json::serialize(canonical_json_value(proof)),
          unit_rational_saturation(dimension, window, context)};
}

std::optional<std::int32_t> parse_optional_match_imaginary_sign(
    const json::object& request, const char* key) {
  const auto* raw = request.if_contains(key);
  if (raw == nullptr || raw->is_null()) return std::nullopt;
  const auto sign = as_i32(*raw, key);
  if (sign != -1 && sign != 1)
    throw std::invalid_argument(std::string(key) + " must be +1 or -1");
  return sign;
}

void require_acb_match_local(const LocalSolution<ComplexBall>& solution,
                             const RealEvaluationPoint& point,
                             const std::string& label) {
  validate_local_solution(solution, true);
  if (point.sign == 0)
    throw std::invalid_argument(
        label +
        " uses the chart center; refined Acb matching requires a nonzero interior point so singular powers are never assigned an artificial center value");
  if (!solution.error.empty())
    throw std::invalid_argument(
        label +
        " carries an error envelope not represented in Acb matching residuals");
  if (!solution.chart.infinite_radius &&
      !arb_lt(acb_realref(point.modulus.raw()),
              acb_realref(solution.chart.radius.raw())))
    throw std::invalid_argument(
        label + " match point is not provably inside its chart radius");
}

Rational acb_physical_match_point(const ChartGeometry& chart,
                                  const RealEvaluationPoint& local_point,
                                  const std::string& label) {
  try {
    return Rational(chart.center_exact) +
           Rational(chart.scale_exact) *
               Rational(local_point.exact_coordinate);
  } catch (const std::invalid_argument&) {
    throw std::invalid_argument(
        label +
        " requires exact rational chart center/scale geometry for a common physical point proof");
  }
}

FiniteLaurentVector<ComplexBall> acb_evaluation_frames(
    const EpsilonVector& value, EpsilonWindow window,
    const std::string& label) {
  if (value.epsilon.complete_max < window.complete_max)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InsufficientCompleteWindow,
        label + " does not cover the requested matching epsilon window",
        std::nullopt, std::nullopt, value.epsilon.complete_max);
  if (window.min_power > value.epsilon.min_power) {
    for (std::int64_t power = value.epsilon.min_power;
         power < window.min_power; ++power)
      for (std::uint32_t component = 0; component < value.dimension;
           ++component)
        if (!value.at(static_cast<std::int32_t>(power), component).is_zero())
          throw MatchingArithmeticError(
              MatchingArithmeticErrorCode::InsufficientCompleteWindow,
              label +
                  " work minimum would discard a nonzero or zero-ambiguous lower epsilon coefficient",
              component, std::nullopt,
              static_cast<std::int32_t>(power));
  }
  FiniteLaurentVector<ComplexBall> frames;
  frames.reserve(value.dimension);
  for (std::uint32_t component = 0; component < value.dimension;
       ++component) {
    std::vector<ComplexBall> coefficients;
    coefficients.reserve(window.width());
    for (std::int64_t power = window.min_power;
         power <= window.complete_max; ++power) {
      const auto epsilon_power = static_cast<std::int32_t>(power);
      coefficients.push_back(epsilon_power < value.epsilon.min_power
          ? ComplexBall(0)
          : value.at(epsilon_power, component));
    }
    frames.emplace_back(window, std::move(coefficients));
  }
  return frames;
}

json::value optional_match_sign_json(
    const std::optional<std::int32_t>& sign) {
  return sign.has_value() ? json::value(*sign) : json::value(nullptr);
}

std::shared_ptr<StoredRefinedAcbMatch> build_refined_acb_match(
    const std::string& match_handle, const json::object& request,
    const std::vector<std::string>& basis_handles,
    const std::vector<std::shared_ptr<StoredLocalBase>>& erased_basis,
    const std::string& incoming_handle,
    const std::shared_ptr<StoredLocalBase>& erased_incoming,
    slong precision_bits,
    const std::string& active_session_configuration_identity,
    const std::optional<json::object>& expected_singular_request =
        std::nullopt) {
  const auto started = std::chrono::steady_clock::now();
  if (request.if_contains("native_singular_scc_saturation") != nullptr &&
      !expected_singular_request.has_value())
    throw std::invalid_argument(
        "native singular-SCC Acb saturation is admitted only through a retained planned match");
  AcbPrecisionLease lease(precision_bits);
  ComplexBall::set_precision(precision_bits);

  const auto& raw_window = as_object(
      request.at("epsilon"), "Acb local match epsilon window");
  EpsilonWindow window{as_i32(raw_window.at("min"), "match epsilon minimum"),
                       as_i32(raw_window.at("max"), "match epsilon maximum")};
  (void)window.width();
  const auto required_complete_max = as_i32(
      raw_window.at("required_complete_max"),
      "required Acb match residual complete maximum");
  if (required_complete_max < window.min_power ||
      required_complete_max > window.complete_max)
    throw std::invalid_argument(
        "required Acb match residual maximum must lie inside the supplied work epsilon window");

  const auto& raw_refinement = as_object(
      request.at("refinement"), "Acb local match refinement policy");
  if (raw_refinement.size() != 2 ||
      raw_refinement.if_contains("relative_tolerance") == nullptr ||
      raw_refinement.if_contains("max_steps") == nullptr)
    throw std::invalid_argument(
        "Acb match refinement accepts exactly relative_tolerance and max_steps");
  const auto relative_tolerance = required_string(
      raw_refinement, "relative_tolerance");
  const auto max_refinement_steps = as_u32(
      raw_refinement.at("max_steps"), "Acb match refinement steps");
  if (max_refinement_steps > 32)
    throw std::invalid_argument(
        "Acb match refinement steps must lie in 0..32");
  AcbLaurentRefinementOptions refinement;
  refinement.relative_tolerance = Magnitude::decimal(relative_tolerance);
  refinement.required_complete_max = required_complete_max;
  refinement.max_refinement_steps = max_refinement_steps;

  const auto basis_point = RealEvaluationPoint::rational(required_string(
      as_object(request.at("basis_point"), "basis match point"), "exact"));
  const auto incoming_point = RealEvaluationPoint::rational(required_string(
      as_object(request.at("incoming_point"), "incoming match point"),
      "exact"));
  const auto requested_basis_sign = parse_optional_match_imaginary_sign(
      request, "basis_imaginary_sign");
  const auto requested_incoming_sign = parse_optional_match_imaginary_sign(
      request, "incoming_imaginary_sign");
  EvaluationOptions basis_options;
  basis_options.imaginary_sign = requested_basis_sign;
  basis_options.compute_tail_estimate = false;
  EvaluationOptions incoming_options;
  incoming_options.imaginary_sign = requested_incoming_sign;
  incoming_options.compute_tail_estimate = false;

  const auto basis_chart = required_string(request, "basis_chart");
  const auto incoming_chart = required_string(request, "incoming_chart");
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  if (checkpoint_identity.empty())
    throw std::invalid_argument(
        "refined Acb match checkpoint identity cannot be empty");

  std::vector<std::shared_ptr<StoredLocal<ComplexBall>>> basis;
  basis.reserve(erased_basis.size());
  for (const auto& local : erased_basis) {
    auto typed = std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(local);
    if (!typed)
      throw std::invalid_argument(
          "refined Acb matching requires Acb retained locals");
    basis.push_back(std::move(typed));
  }
  auto incoming =
      std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(erased_incoming);
  if (!incoming)
    throw std::invalid_argument(
        "refined Acb matching requires an Acb incoming local");
  if (basis.empty())
    throw std::invalid_argument(
        "refined Acb matching requires a nonempty basis");
  const auto dimension = basis.front()->solution().dimension;
  if (basis.size() != dimension || incoming->solution().dimension != dimension)
    throw std::invalid_argument(
        "refined Acb matching requires d basis columns and a d-component incoming local");

  const auto& raw_basis_checkpoints = as_array(
      request.at("basis_checkpoint_identities"),
      "Acb basis checkpoint identities");
  if (raw_basis_checkpoints.size() != basis.size())
    throw std::invalid_argument(
        "Acb basis checkpoint identity count differs from the basis dimension");
  std::vector<std::string> basis_checkpoints;
  basis_checkpoints.reserve(basis.size());
  std::vector<json::object> basis_sources;
  basis_sources.reserve(basis.size());
  for (std::size_t column = 0; column < basis.size(); ++column) {
    if (!raw_basis_checkpoints[column].is_string())
      throw std::invalid_argument(
          "Acb basis checkpoint identities must be strings");
    const std::string expected(raw_basis_checkpoints[column].as_string());
    if (expected.empty() ||
        basis[column]->solution().checkpoint_identity != expected)
      throw std::invalid_argument(
          "Acb basis checkpoint provenance mismatch at column " +
          std::to_string(column));
    if (basis[column]->source_chart() != basis_chart)
      throw std::invalid_argument(
          "Acb basis chart provenance mismatch at column " +
          std::to_string(column));
    if (!same_chart_geometry(basis.front()->solution().chart,
                             basis[column]->solution().chart))
      throw std::invalid_argument(
          "Acb basis locals do not share identical retained chart geometry");
    require_acb_match_local(
        basis[column]->solution(), basis_point,
        "Acb basis local " + basis_handles[column]);
    basis_checkpoints.push_back(expected);
  }
  const auto expected_incoming_checkpoint = required_string(
      request, "incoming_checkpoint_identity");
  if (expected_incoming_checkpoint.empty() ||
      incoming->solution().checkpoint_identity !=
          expected_incoming_checkpoint)
    throw std::invalid_argument(
        "Acb incoming checkpoint provenance mismatch");
  if (incoming->source_chart() != incoming_chart)
    throw std::invalid_argument("Acb incoming chart provenance mismatch");
  require_acb_match_local(incoming->solution(), incoming_point,
                          "Acb incoming local " + incoming_handle);

  const auto basis_physical_point = acb_physical_match_point(
      basis.front()->solution().chart, basis_point, "Acb basis match point");
  for (std::size_t column = 1; column < basis.size(); ++column)
    if (!(acb_physical_match_point(
              basis[column]->solution().chart, basis_point,
              "Acb basis match point at column " +
                  std::to_string(column)) == basis_physical_point))
      throw std::invalid_argument(
          "Acb basis locals do not name one exact physical match point");
  const auto incoming_physical_point = acb_physical_match_point(
      incoming->solution().chart, incoming_point,
      "Acb incoming match point");
  if (!(basis_physical_point == incoming_physical_point))
    throw std::invalid_argument(
        "Acb basis and incoming coordinates do not name the same exact physical match point");

  FiniteLaurentMatrix<ComplexBall> evaluated_basis(
      dimension, FiniteLaurentVector<ComplexBall>());
  for (auto& row : evaluated_basis) row.reserve(dimension);
  std::vector<std::optional<std::int32_t>> effective_basis_signs;
  effective_basis_signs.reserve(dimension);
  for (std::size_t column = 0; column < basis.size(); ++column) {
    const auto evaluation = evaluate_local_solution(
        basis[column]->solution(), basis_point, basis_options);
    effective_basis_signs.push_back(evaluation.imaginary_sign);
    auto frames = acb_evaluation_frames(
        evaluation.value, window,
        "Acb basis evaluation at column " + std::to_string(column));
    for (std::uint32_t component = 0; component < dimension; ++component)
      evaluated_basis[component].push_back(std::move(frames[component]));
  }
  const auto incoming_evaluation = evaluate_local_solution(
      incoming->solution(), incoming_point, incoming_options);
  auto incoming_value = acb_evaluation_frames(
      incoming_evaluation.value, window, "Acb incoming evaluation");

  for (std::size_t column = 0; column < basis.size(); ++column) {
    json::object source{
        {"column", column}, {"local", basis_handles[column]},
        {"chart", basis_chart},
        {"source_operator_identity",
         basis[column]->source_operator_identity()},
        {"checkpoint_identity", basis_checkpoints[column]},
        {"requested_imaginary_sign",
         optional_match_sign_json(requested_basis_sign)},
        {"effective_imaginary_sign",
         optional_match_sign_json(effective_basis_signs[column])},
        {"analytic_metadata", basis[column]->exact_analytic_metadata()}};
    if (basis[column]->column_provenance().has_value())
      source["column_provenance"] =
          basis[column]->column_provenance()->encode();
    basis_sources.push_back(std::move(source));
  }
  json::object incoming_source{
      {"local", incoming_handle}, {"chart", incoming_chart},
      {"source_operator_identity",
       incoming->source_operator_identity()},
      {"checkpoint_identity", expected_incoming_checkpoint},
      {"requested_imaginary_sign",
       optional_match_sign_json(requested_incoming_sign)},
      {"effective_imaginary_sign",
       optional_match_sign_json(incoming_evaluation.imaginary_sign)},
      {"analytic_metadata", incoming->exact_analytic_metadata()}};
  if (incoming->column_provenance().has_value())
    incoming_source["column_provenance"] =
        incoming->column_provenance()->encode();

  auto exact_lattice = [&]() -> ParsedExactEvaluatedLattice {
    const auto proof_request_count =
        (request.if_contains("exact_lattice") != nullptr ? 1U : 0U) +
        (request.if_contains("native_unit_saturation") != nullptr ? 1U : 0U) +
        (request.if_contains("native_singular_scc_saturation") != nullptr
             ? 1U : 0U);
    if (proof_request_count != 1)
      throw std::invalid_argument(
          "Acb matching requires exactly one exact lattice, ordinary native unit-leading request, or singular-SCC valuation-zero request");
    if (const auto* raw_exact = request.if_contains("exact_lattice")) {
      return parse_exact_evaluated_lattice(
          *raw_exact, dimension, window, checkpoint_identity);
    }
    if (const auto* raw_native =
            request.if_contains("native_unit_saturation"))
      return certify_native_unit_saturation(
          *raw_native, evaluated_basis, erased_basis, basis_sources,
          basis_point.exact_coordinate, basis_physical_point.str(), window,
          checkpoint_identity + ":native-unit-leading-proof");
    return certify_native_singular_scc_saturation(
        request.at("native_singular_scc_saturation"), evaluated_basis,
        erased_basis, basis_sources, basis_point.exact_coordinate,
        basis_physical_point.str(), window,
        active_session_configuration_identity,
        *expected_singular_request,
        checkpoint_identity,
        checkpoint_identity + ":native-singular-scc-valuation-zero-proof");
  }();
  json::array exact_binding_basis;
  for (const auto& source : basis_sources)
    exact_binding_basis.push_back(source);
  json::object exact_lattice_provenance{
      {"schema", "diffexp2-retained-exact-lattice-binding-v1"},
      {"witness_schema", exact_lattice.witness_schema},
      {"witness_identity", exact_lattice.identity},
      {"basis", std::move(exact_binding_basis)},
      {"basis_point_exact", basis_point.exact_coordinate},
      {"physical_match_point_exact", basis_physical_point.str()},
      {"epsilon", json::object{{"min", window.min_power},
                                {"max", window.complete_max}}}};
  const auto exact_lattice_provenance_identity = json::serialize(
      canonical_json_value(exact_lattice_provenance));

  auto refined = refine_acb_finite_laurent_match(
      evaluated_basis, incoming_value, exact_lattice.saturation, refinement,
      checkpoint_identity + ":refined-acb-match");

  json::array provenance_basis;
  for (const auto& source : basis_sources) provenance_basis.push_back(source);
  json::object provenance{
      {"schema", "diffexp2-native-refined-acb-local-match-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"basis", std::move(provenance_basis)},
      {"incoming", incoming_source},
      {"basis_point_exact", basis_point.exact_coordinate},
      {"incoming_point_exact", incoming_point.exact_coordinate},
      {"physical_match_point_exact", basis_physical_point.str()},
      {"epsilon", json::object{{"min", window.min_power},
                                {"max", window.complete_max},
                                {"required_complete_max",
                                 required_complete_max}}},
      {"exact_lattice_provenance_identity",
       exact_lattice_provenance_identity},
      {"refinement", json::object{{"relative_tolerance",
                                    relative_tolerance},
                                   {"max_steps",
                                    max_refinement_steps}}}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredRefinedAcbMatch>(
      match_handle, checkpoint_identity, provenance_identity,
      exact_lattice.identity, exact_lattice_provenance_identity,
      std::move(exact_lattice.canonical_witness),
      exact_lattice.witness_schema,
      std::move(basis_sources), std::move(incoming_source), basis_chart,
      incoming_chart, basis_point.exact_coordinate,
      incoming_point.exact_coordinate, basis_physical_point.str(), window,
      required_complete_max, dimension, relative_tolerance,
      max_refinement_steps, std::move(exact_lattice.saturation),
      std::move(refined), elapsed_ms);
}

double checkpoint_nonnegative_double(const json::value& raw,
                                     const char* label) {
  const auto value = as_double(raw, label);
  if (!std::isfinite(value) || value < 0.0)
    throw std::invalid_argument(std::string(label) +
                                " must be finite and nonnegative");
  return value;
}

SCCColumnProvenance parse_checkpoint_column_provenance(
    const json::value& raw) {
  const auto& object = as_object(raw, "checkpoint SCC-column provenance");
  require_exact_keys(object,
      {"scc", "scc_exact_identity", "seed_block", "basis_index",
       "exact_column_identity"}, "checkpoint SCC-column provenance");
  SCCColumnProvenance result{
      required_string(object, "scc"),
      required_string(object, "scc_exact_identity"),
      as_u32(object.at("seed_block"), "checkpoint seed block"),
      as_u32(object.at("basis_index"), "checkpoint basis index"),
      required_string(object, "exact_column_identity")};
  if (result.scc_handle.empty() || result.scc_exact_identity.empty() ||
      result.exact_column_identity.empty())
    throw std::invalid_argument(
        "checkpoint SCC-column provenance contains an empty identity");
  return result;
}

template <typename Scalar>
LocalSolution<Scalar> parse_checkpoint_local_solution(
    const json::value& raw) {
  const auto& object = as_object(raw, "checkpoint local solution");
  require_exact_keys(object,
      {"chart", "epsilon", "taylor_complete_max", "dimension", "sectors",
       "prescriptions", "error", "checkpoint_identity"},
      "checkpoint local solution");
  LocalSolution<Scalar> solution;
  const auto& chart = as_object(object.at("chart"),
                                "checkpoint local chart");
  require_exact_keys(chart,
      {"center_exact", "scale_exact", "radius_exact_ball",
       "infinite_radius"}, "checkpoint local chart");
  solution.chart.center_exact = required_string(chart, "center_exact");
  solution.chart.scale_exact = required_string(chart, "scale_exact");
  if (solution.chart.center_exact.empty() || solution.chart.scale_exact.empty())
    throw std::invalid_argument(
        "checkpoint local chart lost exact center or scale provenance");
  solution.chart.radius = parse_checkpoint_ball(
      chart.at("radius_exact_ball"), "checkpoint chart radius");
  if (!chart.at("infinite_radius").is_bool())
    throw std::invalid_argument(
        "checkpoint infinite-radius flag must be boolean");
  solution.chart.infinite_radius = chart.at("infinite_radius").as_bool();

  const auto& epsilon = as_object(object.at("epsilon"),
                                  "checkpoint local epsilon window");
  require_exact_keys(epsilon, {"min", "max"},
                     "checkpoint local epsilon window");
  solution.epsilon = {
      as_i32(epsilon.at("min"), "checkpoint local epsilon minimum"),
      as_i32(epsilon.at("max"), "checkpoint local epsilon maximum")};
  (void)solution.epsilon.width();
  solution.taylor_complete_max = as_u32(
      object.at("taylor_complete_max"), "checkpoint Taylor complete maximum");
  solution.dimension = as_u32(object.at("dimension"),
                              "checkpoint local dimension");
  if (solution.dimension == 0)
    throw std::invalid_argument("checkpoint local dimension is zero");
  const auto expected = checked_flat_count(
      checked_flat_count(solution.epsilon.width(), solution.taylor_width(),
                         "checkpoint local tensor"),
      solution.dimension, "checkpoint local tensor");

  for (const auto& raw_sector : as_array(object.at("sectors"),
                                          "checkpoint local sectors")) {
    const auto& sector_object = as_object(raw_sector,
                                          "checkpoint local sector");
    require_exact_keys(sector_object,
        {"a", "b", "log_power", "coefficients"},
        "checkpoint local sector");
    LocalSector<Scalar> sector;
    sector.a = parse_checkpoint_exact_descriptor(
        sector_object.at("a"), "checkpoint local a tag");
    sector.b = parse_checkpoint_exact_descriptor(
        sector_object.at("b"), "checkpoint local b tag");
    sector.log_power = as_u32(sector_object.at("log_power"),
                              "checkpoint local log power");
    const auto& coefficients = as_array(
        sector_object.at("coefficients"), "checkpoint sector coefficients");
    if (coefficients.size() != expected)
      throw std::invalid_argument(
          "checkpoint sector coefficient tensor has the wrong size");
    sector.coefficients.reserve(coefficients.size());
    for (const auto& coefficient : coefficients)
      sector.coefficients.push_back(parse_checkpoint_scalar<Scalar>(
          coefficient, "checkpoint local coefficient"));
    solution.sectors.push_back(std::move(sector));
  }
  if (solution.sectors.empty())
    throw std::invalid_argument("checkpoint local solution has no sectors");

  for (const auto& raw_prescription : as_array(
           object.at("prescriptions"), "checkpoint prescriptions")) {
    const auto& prescription = as_object(raw_prescription,
                                         "checkpoint prescription");
    require_exact_keys(prescription,
        {"factor_exact", "sign", "multiplicity",
         "leading_coefficient_sign"}, "checkpoint prescription");
    solution.prescriptions.push_back(Prescription{
        required_string(prescription, "factor_exact"),
        as_i32(prescription.at("sign"), "checkpoint prescription sign"),
        as_u32(prescription.at("multiplicity"),
               "checkpoint prescription multiplicity"),
        as_i32(prescription.at("leading_coefficient_sign"),
               "checkpoint leading-coefficient sign")});
  }
  solution.error = parse_checkpoint_error_envelope(object.at("error"));
  solution.checkpoint_identity = required_string(
      object, "checkpoint_identity");
  if (solution.checkpoint_identity.empty())
    throw std::invalid_argument(
        "checkpoint local solution identity is empty");
  validate_local_solution(solution, false);
  return solution;
}

template <typename Scalar>
std::vector<PseudoHit<Scalar>> parse_checkpoint_pseudo_hits(
    const json::value& raw, std::uint32_t dimension,
    std::uint32_t taylor_complete_max) {
  std::vector<PseudoHit<Scalar>> result;
  for (const auto& raw_hit : as_array(raw, "checkpoint pseudo hits")) {
    const auto& object = as_object(raw_hit, "checkpoint pseudo hit");
    require_exact_keys(object,
        {"n", "columns", "delta_b", "gamma_frames", "gamma_validity"},
        "checkpoint pseudo hit");
    PseudoHit<Scalar> hit;
    hit.n = as_u32(object.at("n"), "checkpoint pseudo-hit Taylor order");
    if (hit.n > taylor_complete_max)
      throw std::invalid_argument(
          "checkpoint pseudo hit lies above the retained Taylor window");
    std::set<std::uint32_t> unique_columns;
    for (const auto& raw_column : as_array(object.at("columns"),
                                            "checkpoint pseudo columns")) {
      const auto column = as_u32(raw_column, "checkpoint pseudo column");
      if (column >= dimension || !unique_columns.insert(column).second)
        throw std::invalid_argument(
            "checkpoint pseudo-hit columns are invalid or duplicated");
      hit.columns.push_back(column);
    }
    if (hit.columns.empty())
      throw std::invalid_argument(
          "checkpoint pseudo hit has an empty Jordan block");
    hit.delta_b = parse_checkpoint_scalar<Scalar>(
        object.at("delta_b"), "checkpoint pseudo delta-b");
    const auto& frames = as_array(object.at("gamma_frames"),
                                  "checkpoint pseudo gamma frames");
    const auto& validity = as_array(object.at("gamma_validity"),
                                    "checkpoint pseudo validity");
    if (frames.size() != hit.columns.size() ||
        validity.size() != hit.columns.size())
      throw std::invalid_argument(
          "checkpoint pseudo-hit block dimensions are inconsistent");
    std::optional<std::size_t> frame_width;
    for (const auto& raw_frame : frames) {
      const auto& coefficients = as_array(raw_frame,
                                          "checkpoint pseudo gamma frame");
      if (coefficients.empty() ||
          (frame_width.has_value() && *frame_width != coefficients.size()))
        throw std::invalid_argument(
            "checkpoint pseudo gamma frames have inconsistent widths");
      frame_width = coefficients.size();
      Frame<Scalar> frame;
      frame.reserve(coefficients.size());
      for (const auto& coefficient : coefficients)
        frame.push_back(parse_checkpoint_scalar<Scalar>(
            coefficient, "checkpoint pseudo gamma coefficient"));
      hit.gamma_frames.push_back(std::move(frame));
    }
    for (const auto& raw_validity : validity)
      hit.gamma_validity.push_back(parse_validity(raw_validity));
    result.push_back(std::move(hit));
  }
  return result;
}

template <typename Scalar>
std::shared_ptr<StoredLocalBase> restore_checkpoint_local_record(
    const json::value& raw, const std::string& expected_domain,
    slong expected_precision_bits,
    std::shared_ptr<void> retained_owner = nullptr) {
  AcbPrecisionLease lease(expected_precision_bits);
  ComplexBall::set_precision(expected_precision_bits);
  const auto& object = as_object(raw, "checkpoint retained local");
  const bool has_tail_restore =
      object.if_contains("tail_model_restore") != nullptr;
  const bool has_derivation_record =
      object.if_contains("retained_derivation") != nullptr;
  if (has_derivation_record !=
      (object.if_contains("retained_owner_lineage") != nullptr))
    throw std::invalid_argument(
        "checkpoint retained-local derivation fields are incomplete");
  if (has_tail_restore && has_derivation_record)
    require_exact_keys(object,
        {"schema", "handle", "source_chart", "source_operator_identity",
         "scalar_domain", "precision_bits", "solution", "pseudo_hits",
         "diagnostics", "runtime_stats", "column_provenance",
         "retained_derivation", "retained_owner_lineage",
         "tail_model_restore"},
        "checkpoint retained local");
  else if (has_tail_restore)
    require_exact_keys(object,
        {"schema", "handle", "source_chart", "source_operator_identity",
         "scalar_domain", "precision_bits", "solution", "pseudo_hits",
         "diagnostics", "runtime_stats", "column_provenance",
         "tail_model_restore"},
        "checkpoint retained local");
  else if (has_derivation_record)
    require_exact_keys(object,
        {"schema", "handle", "source_chart", "source_operator_identity",
         "scalar_domain", "precision_bits", "solution", "pseudo_hits",
         "diagnostics", "runtime_stats", "column_provenance",
         "retained_derivation", "retained_owner_lineage"},
        "checkpoint retained local");
  else
    require_exact_keys(object,
        {"schema", "handle", "source_chart", "source_operator_identity",
         "scalar_domain", "precision_bits", "solution", "pseudo_hits",
         "diagnostics", "runtime_stats", "column_provenance"},
        "checkpoint retained local");
  if (required_string(object, "schema") != "diffexp2-retained-local-v2")
    throw std::invalid_argument("unsupported retained-local checkpoint schema");
  const auto scalar_domain = required_string(object, "scalar_domain");
  const char* compile_domain = std::is_same_v<Scalar, Rational>
      ? "rational" : "acb";
  if (scalar_domain != compile_domain || scalar_domain != expected_domain)
    throw std::invalid_argument(
        "checkpoint local scalar domain is incompatible with its session");
  const auto precision_bits = as_i64(object.at("precision_bits"),
                                     "checkpoint local precision");
  if (precision_bits != expected_precision_bits)
    throw std::invalid_argument(
        "checkpoint local precision differs from its session precision");
  const auto handle = required_string(object, "handle");
  const auto source_chart = required_string(object, "source_chart");
  const auto source_operator_identity = required_string(
      object, "source_operator_identity");
  if (handle.empty() || source_chart.empty() ||
      source_operator_identity.empty())
    throw std::invalid_argument(
        "checkpoint local lost its handle or source-chart provenance");
  auto solution = parse_checkpoint_local_solution<Scalar>(
      object.at("solution"));
  auto pseudo_hits = parse_checkpoint_pseudo_hits<Scalar>(
      object.at("pseudo_hits"), solution.dimension,
      solution.taylor_complete_max);
  std::optional<json::object> retained_derivation;
  const bool has_derivation = has_derivation_record &&
      !object.at("retained_derivation").is_null();
  if (has_derivation_record && has_derivation !=
      !object.at("retained_owner_lineage").is_null())
    throw std::invalid_argument(
        "checkpoint materialized-local derivation and owner lineage disagree");
  if (has_derivation) {
    if (retained_owner == nullptr)
      throw std::invalid_argument(
          "checkpoint derived local lost its strong owner");
    auto derivation = as_object(
        object.at("retained_derivation"),
        "checkpoint retained-local derivation");
    const auto derivation_schema = required_string(derivation, "schema");
    if (derivation_schema ==
        "diffexp2-retained-rational-row-local-application-v1") {
      require_exact_keys(
          derivation,
          {"schema", "capability", "source", "row", "output",
           "analytic_prescriptions", "coefficient_transport",
           "provenance_identity"},
          "checkpoint rational-row local derivation");
      if (required_string(derivation, "capability") !=
              kRetainedRationalRowCapability ||
          required_string(derivation, "analytic_prescriptions") !=
              "preserved-exactly" ||
          required_string(derivation, "coefficient_transport") !=
              "native-retained-only")
        throw std::invalid_argument(
            "checkpoint rational-row derivation changes its retained scope");
      const auto& source = as_object(
          derivation.at("source"), "checkpoint rational-row source");
      require_exact_keys(
          source,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity", "dimension", "epsilon",
           "taylor_complete_max"},
          "checkpoint rational-row source");
      (void)scoped_handle_id(required_string(source, "local"), "l:",
                             "rational-row source local");
      const auto source_chart_identity = required_string(source, "chart");
      if ((!source_chart_identity.starts_with("c:") &&
           !source_chart_identity.starts_with("scc:")) ||
          source_chart_identity != source_chart)
        throw std::invalid_argument(
            "checkpoint rational-row source chart is inconsistent");
      const auto source_operator = required_string(
          source, "source_operator_identity");
      (void)required_string(source, "checkpoint_identity");
      const auto source_dimension = as_u32(
          source.at("dimension"), "checkpoint rational-row source dimension");
      if (source_dimension == 0)
        throw std::invalid_argument(
            "checkpoint rational-row source dimension is zero");
      const auto& source_epsilon = as_object(
          source.at("epsilon"), "checkpoint rational-row source epsilon");
      require_exact_keys(source_epsilon, {"min", "max"},
                         "checkpoint rational-row source epsilon");
      (void)EpsilonWindow{
          as_i32(source_epsilon.at("min"),
                 "checkpoint rational-row source epsilon minimum"),
          as_i32(source_epsilon.at("max"),
                 "checkpoint rational-row source epsilon maximum")}.width();
      (void)as_u32(source.at("taylor_complete_max"),
                   "checkpoint rational-row source Taylor maximum");

      const auto& row = as_object(
          derivation.at("row"), "checkpoint rational-row identity");
      require_exact_keys(
          row, {"exact_identity", "columns", "active_entries",
                "structurally_zero"},
          "checkpoint rational-row identity");
      const auto row_identity = required_string(row, "exact_identity");
      if (as_u32(row.at("columns"),
                 "checkpoint rational-row columns") != source_dimension ||
          !row.at("structurally_zero").is_bool())
        throw std::invalid_argument(
            "checkpoint rational-row dimension or zero fact is malformed");
      const auto& entries = as_array(
          row.at("active_entries"), "checkpoint rational-row entries");
      if (row.at("structurally_zero").as_bool() != entries.empty())
        throw std::invalid_argument(
            "checkpoint rational-row zero fact disagrees with its entries");
      std::optional<std::uint32_t> previous_column;
      for (const auto& raw_entry : entries) {
        const auto& entry = as_object(
            raw_entry, "checkpoint rational-row entry");
        require_exact_keys(
            entry, {"column", "epsilon_shift", "center_pole_order",
                    "exact_identity"},
            "checkpoint rational-row entry");
        const auto column = as_u32(
            entry.at("column"), "checkpoint rational-row entry column");
        if (column >= source_dimension ||
            (previous_column.has_value() && *previous_column >= column))
          throw std::invalid_argument(
              "checkpoint rational-row entry columns are not canonical");
        previous_column = column;
        (void)as_i32(entry.at("epsilon_shift"),
                     "checkpoint rational-row epsilon shift");
        (void)as_u32(entry.at("center_pole_order"),
                     "checkpoint rational-row pole order");
        (void)required_string(entry, "exact_identity");
      }

      const auto& output = as_object(
          derivation.at("output"), "checkpoint rational-row output");
      require_exact_keys(
          output, {"checkpoint_identity", "dimension", "epsilon",
                   "taylor_complete_max"},
          "checkpoint rational-row output");
      const auto& output_epsilon = as_object(
          output.at("epsilon"), "checkpoint rational-row output epsilon");
      require_exact_keys(output_epsilon, {"min", "max"},
                         "checkpoint rational-row output epsilon");
      if (required_string(output, "checkpoint_identity") !=
              solution.checkpoint_identity ||
          as_u32(output.at("dimension"),
                 "checkpoint rational-row output dimension") !=
              solution.dimension ||
          as_i32(output_epsilon.at("min"),
                 "checkpoint rational-row output epsilon minimum") !=
              solution.epsilon.min_power ||
          as_i32(output_epsilon.at("max"),
                 "checkpoint rational-row output epsilon maximum") !=
              solution.epsilon.complete_max ||
          as_u32(output.at("taylor_complete_max"),
                 "checkpoint rational-row output Taylor maximum") !=
              solution.taylor_complete_max)
        throw std::invalid_argument(
            "checkpoint rational-row output disagrees with its tensor");
      auto identity_input = derivation;
      const auto derivation_identity = required_string(
          derivation, "provenance_identity");
      identity_input.erase("provenance_identity");
      if (json::serialize(canonical_json_value(identity_input)) !=
          derivation_identity)
        throw std::invalid_argument(
            "checkpoint rational-row derivation identity is inconsistent");
      const json::object operator_provenance{
          {"schema", "diffexp2-rational-row-derived-operator-v1"},
          {"source_operator_identity", source_operator},
          {"row_exact_identity", row_identity},
          {"provenance_identity", derivation_identity}};
      if (json::serialize(canonical_json_value(operator_provenance)) !=
          source_operator_identity)
        throw std::invalid_argument(
            "checkpoint rational-row derived operator identity is inconsistent");
      const auto& lineage = as_object(
          object.at("retained_owner_lineage"),
          "checkpoint rational-row owner lineage");
      require_exact_keys(
          lineage,
          {"source_local", "source_chart", "source_operator_identity",
           "source_checkpoint_identity", "row_exact_identity",
           "derivation_provenance_identity", "derived_operator_identity"},
          "checkpoint rational-row owner lineage");
      if (lineage.at("source_local") != source.at("local") ||
          lineage.at("source_chart") != source.at("chart") ||
          lineage.at("source_operator_identity") !=
              source.at("source_operator_identity") ||
          lineage.at("source_checkpoint_identity") !=
              source.at("checkpoint_identity") ||
          lineage.at("row_exact_identity") != row.at("exact_identity") ||
          lineage.at("derivation_provenance_identity") !=
              derivation.at("provenance_identity") ||
          required_string(lineage, "derived_operator_identity") !=
              source_operator_identity)
        throw std::invalid_argument(
            "checkpoint rational-row owner lineage is inconsistent");
      retained_derivation = std::move(derivation);
    } else {
      require_exact_keys(
        derivation,
        {"schema", "capability", "source_match",
         "source_match_checkpoint_identity",
         "source_match_provenance_identity",
         "planned_hop_provenance_identity", "planned_hop",
         "weight_windows", "match_certified_complete_max", "output",
         "scope", "coefficient_transport", "whole_arm_complete",
         "provenance_identity"},
        "checkpoint materialized-local derivation");
      if (required_string(derivation, "schema") !=
            "diffexp2-retained-plan-match-local-materialization-v1" ||
        required_string(derivation, "capability") !=
            "retained-native-plan-match-local-materialization-v1" ||
        required_string(derivation, "scope") !=
            "single-match-receiving-local" ||
        required_string(derivation, "coefficient_transport") !=
            "native-retained-only" ||
        !derivation.at("whole_arm_complete").is_bool() ||
        derivation.at("whole_arm_complete").as_bool())
      throw std::invalid_argument(
          "checkpoint materialized-local derivation changes its certified scope");
    (void)scoped_handle_id(required_string(derivation, "source_match"),
                           "m:", "materialized-local source match");
    const auto& output = as_object(
        derivation.at("output"), "checkpoint materialized-local output");
    require_exact_keys(output,
        {"checkpoint_identity", "chart", "source_operator_identity",
         "epsilon", "taylor_complete_max", "dimension"},
        "checkpoint materialized-local output");
    const auto& output_epsilon = as_object(
        output.at("epsilon"), "checkpoint materialized-local epsilon");
    require_exact_keys(output_epsilon, {"min", "max"},
                       "checkpoint materialized-local epsilon");
    if (required_string(output, "checkpoint_identity") !=
            solution.checkpoint_identity ||
        required_string(output, "chart") != source_chart ||
        required_string(output, "source_operator_identity") !=
            source_operator_identity ||
        as_i32(output_epsilon.at("min"),
               "materialized-local epsilon minimum") !=
            solution.epsilon.min_power ||
        as_i32(output_epsilon.at("max"),
               "materialized-local epsilon maximum") !=
            solution.epsilon.complete_max ||
        as_u32(output.at("taylor_complete_max"),
               "materialized-local Taylor maximum") !=
            solution.taylor_complete_max ||
        as_u32(output.at("dimension"),
               "materialized-local dimension") != solution.dimension)
      throw std::invalid_argument(
          "checkpoint materialized-local output provenance disagrees with its tensor");
    for (const auto& raw_window : as_array(
             derivation.at("weight_windows"),
             "checkpoint materialization weight windows")) {
      const auto& window = as_object(
          raw_window, "checkpoint materialization weight window");
      require_exact_keys(window, {"min", "max"},
                         "checkpoint materialization weight window");
      (void)EpsilonWindow{
          as_i32(window.at("min"), "materialization weight minimum"),
          as_i32(window.at("max"), "materialization weight maximum")}.width();
    }
    auto identity_input = derivation;
    const auto derivation_identity = required_string(
        derivation, "provenance_identity");
    identity_input.erase("provenance_identity");
    if (json::serialize(canonical_json_value(identity_input)) !=
        derivation_identity)
      throw std::invalid_argument(
          "checkpoint materialized-local derivation identity is inconsistent");
    const auto& lineage = as_object(
        object.at("retained_owner_lineage"),
        "checkpoint materialized-local owner lineage");
    require_exact_keys(
        lineage,
        {"match", "match_checkpoint_identity", "match_provenance_identity",
         "planned_hop_provenance_identity",
         "derivation_provenance_identity"},
        "checkpoint materialized-local owner lineage");
    if (lineage.at("match") != derivation.at("source_match") ||
        lineage.at("match_checkpoint_identity") !=
            derivation.at("source_match_checkpoint_identity") ||
        lineage.at("match_provenance_identity") !=
            derivation.at("source_match_provenance_identity") ||
        lineage.at("planned_hop_provenance_identity") !=
            derivation.at("planned_hop_provenance_identity") ||
        lineage.at("derivation_provenance_identity") !=
            derivation.at("provenance_identity"))
      throw std::invalid_argument(
          "checkpoint materialized-local owner lineage is inconsistent");
    retained_derivation = std::move(derivation);
    }
  } else if (retained_owner != nullptr) {
    throw std::invalid_argument(
        "checkpoint primitive local unexpectedly acquired a derivation owner");
  }
  const auto& diagnostics = as_object(object.at("diagnostics"),
                                      "checkpoint local diagnostics");
  require_exact_keys(diagnostics,
      {"top_valid", "create_parse_ms", "create_kernel_ms"},
      "checkpoint local diagnostics");
  NativeLocalDiagnostics native{
      parse_validity(diagnostics.at("top_valid")),
      checkpoint_nonnegative_double(diagnostics.at("create_parse_ms"),
                                    "checkpoint local parse time"),
      checkpoint_nonnegative_double(diagnostics.at("create_kernel_ms"),
                                    "checkpoint local kernel time")};
  std::optional<SCCColumnProvenance> column_provenance;
  if (!object.at("column_provenance").is_null())
    column_provenance = parse_checkpoint_column_provenance(
        object.at("column_provenance"));
  std::string saved_tail_status = "unrecorded";
  bool saved_tail_attached = false;
  RegularTaylorTailModelResult restored_tail_model = unavailable_tail_model(
      "checkpoint has no serialized regular tail model");
  std::optional<TailModelCheckpointMarker> tail_checkpoint_marker;
  if (has_tail_restore) {
    const auto& tail = as_object(
        object.at("tail_model_restore"),
        "checkpoint tail-model restore marker");
    const auto* raw_serialized = tail.if_contains("serialized");
    const auto* raw_attached = tail.if_contains("attached_before_save");
    if (raw_serialized == nullptr || !raw_serialized->is_bool() ||
        raw_attached == nullptr || !raw_attached->is_bool() ||
        required_string(tail, "capability") !=
            kRegularTailMajorantCapability)
      throw std::invalid_argument(
          "checkpoint tail-model restore marker is incompatible");
    const auto serialized = raw_serialized->as_bool();
    saved_tail_attached = raw_attached->as_bool();
    saved_tail_status = required_string(tail, "status");
    if (saved_tail_status != "certified" &&
        saved_tail_status != "inconclusive" &&
        saved_tail_status != "unsupported")
      throw std::invalid_argument(
          "checkpoint tail-model restore marker has an unknown status");
    if (serialized) {
      require_exact_keys(
          tail,
          {"capability", "serialized", "status", "attached_before_save",
           "model"},
          "checkpoint serialized tail model");
      if (saved_tail_status != "certified" || !saved_tail_attached)
        throw std::invalid_argument(
            "serialized checkpoint tail model is not attached/certified");
      auto model = parse_checkpoint_regular_tail_model(tail.at("model"));
      tail_majorant_detail::validate_restored_regular_taylor_tail_model(
          model, solution, source_operator_identity);
      restored_tail_model = {
          TailMajorantStatus::Certified, std::move(model),
          "certified regular tail model restored from exact checkpoint state"};
    } else {
      require_exact_keys(
          tail,
          {"capability", "serialized", "status", "attached_before_save"},
          "checkpoint tail-model restore marker");
      restored_tail_model = unavailable_tail_model(
          "checkpoint does not serialize this regular tail-model state; "
          "saved model status was " + saved_tail_status +
          "; re-solve the retained local to reattach certification state");
      tail_checkpoint_marker = TailModelCheckpointMarker{
          saved_tail_status, saved_tail_attached};
    }
  }

  const auto& stats = as_object(object.at("runtime_stats"),
                                "checkpoint local runtime stats");
  if (has_tail_restore)
    require_exact_keys(stats,
        {"evaluations", "residual_certifications", "endpoint_limits",
         "line_integrations", "evaluate_ms", "residual_certify_ms",
         "endpoint_limit_ms", "line_integration_ms",
         "coefficient_count", "tail_certificate_requests",
         "tail_certificate_certified", "tail_certificate_inconclusive",
         "tail_certificate_unsupported"},
        "checkpoint local runtime stats");
  else
    require_exact_keys(stats,
        {"evaluations", "residual_certifications", "endpoint_limits",
         "line_integrations", "evaluate_ms", "residual_certify_ms",
         "endpoint_limit_ms", "line_integration_ms",
         "coefficient_count"}, "checkpoint local runtime stats");
  StoredLocalStats restored_stats;
  restored_stats.evaluations = as_u64(stats.at("evaluations"),
                                     "checkpoint local evaluations");
  restored_stats.residual_certifications = as_u64(
      stats.at("residual_certifications"),
      "checkpoint local residual certifications");
  restored_stats.endpoint_limits = as_u64(
      stats.at("endpoint_limits"), "checkpoint local endpoint limits");
  restored_stats.line_integrations = as_u64(
      stats.at("line_integrations"), "checkpoint local line integrations");
  restored_stats.evaluate_ms = checkpoint_nonnegative_double(
      stats.at("evaluate_ms"), "checkpoint local evaluation time");
  restored_stats.residual_certify_ms = checkpoint_nonnegative_double(
      stats.at("residual_certify_ms"),
      "checkpoint local residual-certification time");
  restored_stats.endpoint_limit_ms = checkpoint_nonnegative_double(
      stats.at("endpoint_limit_ms"), "checkpoint local endpoint time");
  restored_stats.line_integration_ms = checkpoint_nonnegative_double(
      stats.at("line_integration_ms"),
      "checkpoint local line-integration time");
  restored_stats.create_parse_ms = native.parse_ms;
  restored_stats.create_kernel_ms = native.kernel_ms;
  const auto coefficient_count = as_u64(
      stats.at("coefficient_count"), "checkpoint local coefficient count");
  if (coefficient_count > std::numeric_limits<std::size_t>::max())
    throw std::invalid_argument(
        "checkpoint local coefficient count exceeds size_t");
  restored_stats.coefficient_count =
      static_cast<std::size_t>(coefficient_count);
  if (has_tail_restore) {
    restored_stats.tail_certificate_requests = as_u64(
        stats.at("tail_certificate_requests"),
        "checkpoint tail certificate requests");
    restored_stats.tail_certificate_certified = as_u64(
        stats.at("tail_certificate_certified"),
        "checkpoint certified tail certificates");
    restored_stats.tail_certificate_inconclusive = as_u64(
        stats.at("tail_certificate_inconclusive"),
        "checkpoint inconclusive tail certificates");
    restored_stats.tail_certificate_unsupported = as_u64(
        stats.at("tail_certificate_unsupported"),
        "checkpoint unsupported tail certificates");
  }

  auto local = make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
      handle, source_chart, source_operator_identity, std::move(solution),
      expected_precision_bits,
      std::move(pseudo_hits), native, std::move(column_provenance),
      std::move(retained_derivation), std::move(retained_owner),
      std::move(restored_tail_model), std::move(tail_checkpoint_marker),
      has_tail_restore, has_derivation_record);
  local->restore_runtime_stats(restored_stats);
  return local;
}

AcbMatchingResidualVerdict parse_checkpoint_acb_match_verdict(
    const json::value& raw, const char* label) {
  if (!raw.is_string())
    throw std::invalid_argument(std::string(label) + " must be a string");
  const auto value = std::string(raw.as_string());
  if (value == "pass") return AcbMatchingResidualVerdict::Pass;
  if (value == "fail") return AcbMatchingResidualVerdict::Fail;
  if (value == "inconclusive")
    return AcbMatchingResidualVerdict::Inconclusive;
  throw std::invalid_argument(std::string(label) +
                              " has an unsupported verdict");
}

AcbMatchingResidualDiagnostics parse_checkpoint_acb_match_residual(
    const json::value& raw, std::uint32_t dimension,
    std::int32_t expected_required_complete_max) {
  const auto& object = as_object(raw, "checkpoint Acb match residual");
  require_exact_keys(object,
      {"verdict", "complete_window", "required_complete_max",
       "complete_through_required", "coefficients", "detail"},
      "checkpoint Acb match residual");
  AcbMatchingResidualDiagnostics result;
  result.verdict = parse_checkpoint_acb_match_verdict(
      object.at("verdict"), "checkpoint Acb residual verdict");
  const auto& window = as_object(object.at("complete_window"),
                                 "checkpoint Acb residual window");
  require_exact_keys(window, {"min", "max"},
                     "checkpoint Acb residual window");
  result.complete_window = {
      as_i32(window.at("min"), "checkpoint Acb residual minimum"),
      as_i32(window.at("max"), "checkpoint Acb residual maximum")};
  const auto width = result.complete_window.width();
  result.required_complete_max = as_i32(
      object.at("required_complete_max"),
      "checkpoint Acb required complete maximum");
  if (result.required_complete_max != expected_required_complete_max)
    throw std::invalid_argument(
        "checkpoint Acb residual requirement changed across history");
  if (!object.at("complete_through_required").is_bool())
    throw std::invalid_argument(
        "checkpoint Acb residual completeness flag must be boolean");
  result.complete_through_required =
      object.at("complete_through_required").as_bool();
  if (result.complete_through_required !=
      (result.complete_window.complete_max >= result.required_complete_max))
    throw std::invalid_argument(
        "checkpoint Acb residual completeness flag contradicts its window");
  const auto& coefficients = as_array(object.at("coefficients"),
                                      "checkpoint Acb residual coefficients");
  if (coefficients.size() != checked_flat_count(
          dimension, width, "checkpoint Acb residual diagnostics"))
    throw std::invalid_argument(
        "checkpoint Acb residual coefficient history is incomplete");
  bool any_fail = false;
  bool all_pass = true;
  result.coefficients.reserve(coefficients.size());
  for (std::size_t index = 0; index < coefficients.size(); ++index) {
    const auto& coefficient = as_object(
        coefficients[index], "checkpoint Acb residual coefficient");
    require_exact_keys(coefficient,
        {"row", "epsilon_power", "residual_lower_exact",
         "residual_upper_exact", "scale_lower_exact", "scale_upper_exact",
         "verdict"}, "checkpoint Acb residual coefficient");
    const auto row = as_u64(coefficient.at("row"),
                            "checkpoint Acb residual row");
    const auto power = as_i32(coefficient.at("epsilon_power"),
                              "checkpoint Acb residual power");
    const auto expected_row = index / width;
    const auto expected_power = static_cast<std::int32_t>(
        static_cast<std::int64_t>(result.complete_window.min_power) +
        static_cast<std::int64_t>(index % width));
    if (row != expected_row || power != expected_power)
      throw std::invalid_argument(
          "checkpoint Acb residual coefficients are not a complete row-major history");
    AcbMatchingCoefficientResidual parsed;
    parsed.row = static_cast<std::size_t>(row);
    parsed.epsilon_power = power;
    parsed.residual_lower = parse_checkpoint_magnitude(
        coefficient.at("residual_lower_exact"),
        "checkpoint residual lower bound");
    parsed.residual_upper = parse_checkpoint_magnitude(
        coefficient.at("residual_upper_exact"),
        "checkpoint residual upper bound");
    parsed.scale_lower = parse_checkpoint_magnitude(
        coefficient.at("scale_lower_exact"),
        "checkpoint residual scale lower bound");
    parsed.scale_upper = parse_checkpoint_magnitude(
        coefficient.at("scale_upper_exact"),
        "checkpoint residual scale upper bound");
    parsed.verdict = parse_checkpoint_acb_match_verdict(
        coefficient.at("verdict"),
        "checkpoint Acb coefficient verdict");
    if (power <= result.required_complete_max) {
      any_fail = any_fail ||
          parsed.verdict == AcbMatchingResidualVerdict::Fail;
      all_pass = all_pass &&
          parsed.verdict == AcbMatchingResidualVerdict::Pass;
    }
    result.coefficients.push_back(std::move(parsed));
  }
  const auto derived_verdict = any_fail
      ? AcbMatchingResidualVerdict::Fail
      : (all_pass && result.complete_through_required)
          ? AcbMatchingResidualVerdict::Pass
          : AcbMatchingResidualVerdict::Inconclusive;
  if (derived_verdict != result.verdict)
    throw std::invalid_argument(
        "checkpoint Acb residual aggregate verdict is inconsistent");
  result.detail = required_string(object, "detail");
  if (result.detail.empty())
    throw std::invalid_argument(
        "checkpoint Acb residual history lost its diagnostic detail");
  return result;
}

void validate_checkpoint_exact_analytic_metadata(const json::value& raw) {
  const auto& metadata = as_object(
      raw, "checkpoint exact local analytic metadata");
  require_exact_keys(metadata,
      {"schema", "chart", "sectors", "prescriptions"},
      "checkpoint exact local analytic metadata");
  if (required_string(metadata, "schema") !=
      "diffexp2-exact-local-analytic-metadata-v2")
    throw std::invalid_argument(
        "unsupported checkpoint local analytic metadata schema");
  const auto& chart = as_object(metadata.at("chart"),
                                "checkpoint exact analytic chart");
  require_exact_keys(chart,
      {"center_exact", "scale_exact", "radius_exact_ball",
       "infinite_radius"}, "checkpoint exact analytic chart");
  (void)required_string(chart, "center_exact");
  (void)required_string(chart, "scale_exact");
  (void)parse_checkpoint_ball(chart.at("radius_exact_ball"),
                              "checkpoint exact analytic radius");
  if (!chart.at("infinite_radius").is_bool())
    throw std::invalid_argument(
        "checkpoint exact analytic radius flag must be boolean");
  const auto& sectors = as_array(metadata.at("sectors"),
                                 "checkpoint exact analytic sectors");
  if (sectors.empty())
    throw std::invalid_argument(
        "checkpoint exact analytic metadata has no sectors");
  for (const auto& raw_sector : sectors) {
    const auto& sector = as_object(raw_sector,
                                   "checkpoint exact analytic sector");
    require_exact_keys(sector, {"a", "b", "log_power"},
                       "checkpoint exact analytic sector");
    (void)parse_checkpoint_exact_descriptor(
        sector.at("a"), "checkpoint exact analytic a tag");
    (void)parse_checkpoint_exact_descriptor(
        sector.at("b"), "checkpoint exact analytic b tag");
    (void)as_u32(sector.at("log_power"),
                 "checkpoint exact analytic log power");
  }
  for (const auto& raw_prescription : as_array(
           metadata.at("prescriptions"),
           "checkpoint exact analytic prescriptions")) {
    const auto& prescription = as_object(
        raw_prescription, "checkpoint exact analytic prescription");
    require_exact_keys(prescription,
        {"factor_exact", "sign", "multiplicity",
         "leading_coefficient_sign"},
        "checkpoint exact analytic prescription");
    (void)required_string(prescription, "factor_exact");
    const auto sign = as_i32(prescription.at("sign"),
                             "checkpoint analytic prescription sign");
    const auto leading = as_i32(
        prescription.at("leading_coefficient_sign"),
        "checkpoint analytic leading sign");
    const auto multiplicity = as_u32(
        prescription.at("multiplicity"),
        "checkpoint analytic prescription multiplicity");
    if ((sign != -1 && sign != 1) ||
        (leading != -1 && leading != 1) || multiplicity == 0)
      throw std::invalid_argument(
          "checkpoint exact analytic prescription is malformed");
  }
}

std::shared_ptr<StoredEndpointResult> restore_checkpoint_endpoint_record(
    const json::value& raw, const std::string& expected_domain) {
  const auto& object = as_object(raw, "checkpoint retained endpoint");
  require_exact_keys(object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "source", "approach_direction", "requested_rim",
       "cancellation_mode", "analytic_metadata", "result", "elapsed_ms",
       "runtime_stats"}, "checkpoint retained endpoint");
  if (required_string(object, "schema") !=
      "diffexp2-retained-endpoint-result-v2")
    throw std::invalid_argument(
        "unsupported retained endpoint checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  const auto& source = as_object(object.at("source"),
                                 "checkpoint endpoint source");
  require_exact_keys(source,
      {"local", "chart", "source_operator_identity",
       "checkpoint_identity", "coefficient_domain"},
      "checkpoint endpoint source");
  const auto source_local = required_string(source, "local");
  const auto source_chart = required_string(source, "chart");
  const auto source_operator_identity = required_string(
      source, "source_operator_identity");
  const auto source_checkpoint = required_string(
      source, "checkpoint_identity");
  const auto source_domain = required_string(source, "coefficient_domain");
  if (source_domain != expected_domain ||
      (source_domain != "rational" && source_domain != "acb"))
    throw std::invalid_argument(
        "checkpoint endpoint coefficient domain is incompatible with its session");
  const auto approach_direction = as_i32(
      object.at("approach_direction"),
      "checkpoint endpoint approach direction");
  if (approach_direction != -1 && approach_direction != 1)
    throw std::invalid_argument(
        "checkpoint endpoint approach direction must be +1 or -1");
  std::optional<std::int32_t> requested_rim;
  if (!object.at("requested_rim").is_null()) {
    requested_rim = as_i32(object.at("requested_rim"),
                           "checkpoint endpoint requested rim");
    if (*requested_rim != -1 && *requested_rim != 1)
      throw std::invalid_argument(
          "checkpoint endpoint rim must be +1 or -1");
  }
  const auto cancellation_mode = required_string(
      object, "cancellation_mode");
  if (cancellation_mode != "exact-coefficient-field" &&
      cancellation_mode != "exact-or-acb-singleton")
    throw std::invalid_argument(
        "checkpoint endpoint cancellation mode is unsupported");
  auto analytic_metadata = as_object(
      object.at("analytic_metadata"),
      "checkpoint endpoint analytic metadata");
  validate_checkpoint_exact_analytic_metadata(analytic_metadata);

  json::object provenance{
      {"schema", "diffexp2-retained-native-endpoint-sector-limit-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"source", source},
      {"approach_direction", approach_direction},
      {"rim", requested_rim.has_value()
           ? json::value(*requested_rim) : json::value(nullptr)},
      {"cancellation", json::object{{"mode", cancellation_mode}}},
      {"analytic_metadata", analytic_metadata}};
  if (json::serialize(canonical_json_value(provenance)) !=
      provenance_identity)
    throw std::invalid_argument(
        "checkpoint endpoint provenance identity is inconsistent");

  const auto& raw_result = as_object(object.at("result"),
                                     "checkpoint endpoint result");
  require_exact_keys(raw_result,
      {"values", "dropped_regulated_sectors",
       "cancelled_divergent_coefficients", "imaginary_sign"},
      "checkpoint endpoint result");
  EndpointLimitResult result;
  for (const auto& raw_value : as_array(raw_result.at("values"),
                                        "checkpoint endpoint values"))
    result.values.push_back(parse_checkpoint_epsilon_frame<ComplexBall>(
        raw_value, "checkpoint endpoint value"));
  if (result.values.empty())
    throw std::invalid_argument(
        "checkpoint endpoint result has no components");
  const auto checked_size = [](const json::value& value,
                               const char* label) {
    const auto parsed = as_u64(value, label);
    if (parsed > std::numeric_limits<std::size_t>::max())
      throw std::invalid_argument(std::string(label) +
                                  " exceeds size_t");
    return static_cast<std::size_t>(parsed);
  };
  result.dropped_regulated_sectors = checked_size(
      raw_result.at("dropped_regulated_sectors"),
      "checkpoint dropped regulated sectors");
  result.cancelled_divergent_coefficients = checked_size(
      raw_result.at("cancelled_divergent_coefficients"),
      "checkpoint cancelled divergent coefficients");
  result.imaginary_sign = as_i32(raw_result.at("imaginary_sign"),
                                 "checkpoint endpoint effective rim");
  if (result.imaginary_sign != -1 && result.imaginary_sign != 1)
    throw std::invalid_argument(
        "checkpoint endpoint effective rim must be +1 or -1");
  (void)endpoint_value_window(result);
  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint endpoint elapsed time");
  const auto& stats = as_object(object.at("runtime_stats"),
                                "checkpoint endpoint runtime stats");
  require_exact_keys(stats, {"exports", "export_ms"},
                     "checkpoint endpoint runtime stats");
  const auto exports = as_u64(stats.at("exports"),
                              "checkpoint endpoint exports");
  const auto export_ms = checkpoint_nonnegative_double(
      stats.at("export_ms"), "checkpoint endpoint export time");
  auto endpoint = std::make_shared<StoredEndpointResult>(
      handle, checkpoint_identity, provenance_identity, source_local,
      source_chart, source_operator_identity, source_checkpoint,
      source_domain, approach_direction, requested_rim, cancellation_mode,
      std::move(analytic_metadata), std::move(result), elapsed_ms);
  endpoint->restore_runtime_stats(exports, export_ms);
  return endpoint;
}

void validate_checkpoint_match_source(const json::object& source,
                                      bool basis,
                                      std::size_t expected_column = 0) {
  const bool has_column_provenance =
      source.if_contains("column_provenance") != nullptr;
  if (basis) {
    if (has_column_provenance)
      require_exact_keys(source,
          {"column", "local", "chart", "source_operator_identity",
           "checkpoint_identity",
           "requested_imaginary_sign", "effective_imaginary_sign",
           "analytic_metadata", "column_provenance"},
          "checkpoint Acb basis source");
    else
      require_exact_keys(source,
          {"column", "local", "chart", "source_operator_identity",
           "checkpoint_identity",
           "requested_imaginary_sign", "effective_imaginary_sign",
           "analytic_metadata"}, "checkpoint Acb basis source");
    if (as_u64(source.at("column"), "checkpoint Acb basis column") !=
        expected_column)
      throw std::invalid_argument(
          "checkpoint Acb basis sources are not in column order");
  } else {
    if (has_column_provenance)
      require_exact_keys(source,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity",
           "requested_imaginary_sign", "effective_imaginary_sign",
           "analytic_metadata", "column_provenance"},
          "checkpoint Acb incoming source");
    else
      require_exact_keys(source,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity",
           "requested_imaginary_sign", "effective_imaginary_sign",
           "analytic_metadata"}, "checkpoint Acb incoming source");
  }
  if (required_string(source, "local").empty() ||
      required_string(source, "chart").empty() ||
      required_string(source, "source_operator_identity").empty() ||
      required_string(source, "checkpoint_identity").empty())
    throw std::invalid_argument(
        "checkpoint Acb match source has empty provenance");
  for (const auto* key : {"requested_imaginary_sign",
                          "effective_imaginary_sign"}) {
    if (source.at(key).is_null()) continue;
    const auto sign = as_i32(source.at(key), key);
    if (sign != -1 && sign != 1)
      throw std::invalid_argument(
          "checkpoint Acb match source has an invalid branch sign");
  }
  validate_checkpoint_exact_analytic_metadata(
      source.at("analytic_metadata"));
  if (has_column_provenance)
    (void)parse_checkpoint_column_provenance(
        source.at("column_provenance"));
}

void validate_checkpoint_exact_match_source(
    const json::object& source, bool basis_source,
    std::size_t expected_column,
    const std::shared_ptr<StoredLocalBase>& owner) {
  const bool has_column_provenance =
      source.if_contains("column_provenance") != nullptr;
  if (basis_source) {
    if (has_column_provenance)
      require_exact_keys(
          source,
          {"column", "local", "chart", "source_operator_identity",
           "checkpoint_identity", "analytic_metadata",
           "column_provenance"},
          "checkpoint exact-match basis source");
    else
      require_exact_keys(
          source,
          {"column", "local", "chart", "source_operator_identity",
           "checkpoint_identity", "analytic_metadata"},
          "checkpoint exact-match basis source");
    if (as_u64(source.at("column"),
               "checkpoint exact-match basis column") != expected_column)
      throw std::invalid_argument(
          "checkpoint exact-match basis sources are not in column order");
  } else {
    if (has_column_provenance)
      require_exact_keys(
          source,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity", "analytic_metadata",
           "column_provenance"},
          "checkpoint exact-match incoming source");
    else
      require_exact_keys(
          source,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity", "analytic_metadata"},
          "checkpoint exact-match incoming source");
  }
  if (!owner || std::string(owner->scalar_domain()) != "rational")
    throw std::invalid_argument(
        "checkpoint exact-match source has no exact-rational owner");
  if (required_string(source, "local") != owner->handle() ||
      required_string(source, "chart") != owner->source_chart() ||
      required_string(source, "source_operator_identity") !=
          owner->source_operator_identity() ||
      required_string(source, "checkpoint_identity") !=
          owner->checkpoint_identity() ||
      source.at("analytic_metadata") != owner->exact_analytic_metadata())
    throw std::invalid_argument(
        "checkpoint exact-match source provenance disagrees with its retained local owner");
  validate_checkpoint_exact_analytic_metadata(
      source.at("analytic_metadata"));
  if (has_column_provenance)
    (void)parse_checkpoint_column_provenance(
        source.at("column_provenance"));
}

std::shared_ptr<StoredExactRegularMatch>
restore_checkpoint_exact_match_record(
    const json::value& raw,
    std::vector<std::shared_ptr<StoredLocalBase>> basis_owners,
    std::shared_ptr<StoredLocalBase> incoming_owner) {
  const auto& object = as_object(raw,
                                 "checkpoint retained exact-rational match");
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "basis_sources", "incoming_source", "basis_chart",
       "incoming_chart", "basis_point_exact", "incoming_point_exact",
       "physical_match_point_exact", "epsilon", "dimension",
       "transformation", "weights", "saturation", "residual_window",
       "elapsed_ms"},
      "checkpoint retained exact-rational match");
  if (required_string(object, "schema") !=
      "diffexp2-retained-exact-rational-match-v2")
    throw std::invalid_argument(
        "unsupported retained exact-rational-match checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  const auto basis_chart = required_string(object, "basis_chart");
  const auto incoming_chart = required_string(object, "incoming_chart");
  const auto basis_point = required_string(object, "basis_point_exact");
  const auto incoming_point = required_string(object,
                                               "incoming_point_exact");
  const auto physical_point = required_string(
      object, "physical_match_point_exact");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty() || basis_chart.empty() ||
      incoming_chart.empty() || basis_point.empty() ||
      incoming_point.empty() || physical_point.empty())
    throw std::invalid_argument(
        "checkpoint exact-rational match contains an empty identity");
  (void)Rational(basis_point);
  (void)Rational(incoming_point);
  (void)Rational(physical_point);
  const auto dimension = as_u32(object.at("dimension"),
                                "checkpoint exact-match dimension");
  if (dimension == 0 || basis_owners.size() != dimension ||
      !incoming_owner)
    throw std::invalid_argument(
        "checkpoint exact-match ownership differs from its dimension");
  std::vector<json::object> basis_sources;
  const auto& raw_basis = as_array(object.at("basis_sources"),
                                   "checkpoint exact-match basis sources");
  if (raw_basis.size() != dimension)
    throw std::invalid_argument(
        "checkpoint exact-match source count differs from its dimension");
  basis_sources.reserve(dimension);
  for (std::size_t column = 0; column < dimension; ++column) {
    auto source = as_object(raw_basis[column],
                            "checkpoint exact-match basis source");
    validate_checkpoint_exact_match_source(
        source, true, column, basis_owners[column]);
    if (required_string(source, "chart") != basis_chart)
      throw std::invalid_argument(
          "checkpoint exact-match basis chart provenance is inconsistent");
    basis_sources.push_back(std::move(source));
  }
  auto incoming_source = as_object(
      object.at("incoming_source"),
      "checkpoint exact-match incoming source");
  validate_checkpoint_exact_match_source(
      incoming_source, false, 0, incoming_owner);
  if (required_string(incoming_source, "chart") != incoming_chart)
    throw std::invalid_argument(
        "checkpoint exact-match incoming chart provenance is inconsistent");

  const auto& raw_epsilon = as_object(object.at("epsilon"),
                                      "checkpoint exact-match epsilon");
  require_exact_keys(raw_epsilon, {"min", "max", "required_complete_max"},
                     "checkpoint exact-match epsilon");
  EpsilonWindow window{
      as_i32(raw_epsilon.at("min"), "checkpoint exact-match epsilon minimum"),
      as_i32(raw_epsilon.at("max"), "checkpoint exact-match epsilon maximum")};
  (void)window.width();
  const auto required_complete_max = as_i32(
      raw_epsilon.at("required_complete_max"),
      "checkpoint exact-match required complete maximum");
  if (required_complete_max < window.min_power ||
      required_complete_max > window.complete_max)
    throw std::invalid_argument(
        "checkpoint exact-match requirement lies outside its work window");
  auto transformation = parse_checkpoint_exact_laurent_matrix(
      object.at("transformation"), dimension,
      "checkpoint exact-match transformation");
  auto weights = parse_checkpoint_frame_vector<Rational>(
      object.at("weights"), dimension, "checkpoint exact-match weights");
  auto diagnostics = parse_checkpoint_saturation_diagnostics(
      object.at("saturation"), dimension);
  const auto& raw_residual = as_object(
      object.at("residual_window"), "checkpoint exact-match residual window");
  require_exact_keys(raw_residual, {"min", "max"},
                     "checkpoint exact-match residual window");
  EpsilonWindow residual{
      as_i32(raw_residual.at("min"), "checkpoint residual minimum"),
      as_i32(raw_residual.at("max"), "checkpoint residual maximum")};
  (void)residual.width();
  if (residual.complete_max < required_complete_max)
    throw std::invalid_argument(
        "checkpoint exact-match residual lost its required complete window");

  json::array provenance_basis;
  for (const auto& source : basis_sources)
    provenance_basis.push_back(source);
  json::object provenance{
      {"schema", "diffexp2-native-exact-regular-local-match-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"basis", std::move(provenance_basis)},
      {"incoming", incoming_source},
      {"basis_point_exact", basis_point},
      {"incoming_point_exact", incoming_point},
      {"physical_match_point_exact", physical_point},
      {"epsilon", json::object{{"min", window.min_power},
                                {"max", window.complete_max},
                                {"required_complete_max",
                                 required_complete_max}}}};
  if (json::serialize(canonical_json_value(provenance)) !=
      provenance_identity)
    throw std::invalid_argument(
        "checkpoint exact-match provenance identity is inconsistent");
  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint exact-match elapsed time");
  return std::make_shared<StoredExactRegularMatch>(
      handle, checkpoint_identity, provenance_identity,
      std::move(basis_sources), std::move(incoming_source), basis_chart,
      incoming_chart, basis_point, incoming_point, physical_point, window,
      required_complete_max, dimension, std::move(transformation),
      std::move(weights), std::move(diagnostics), residual, elapsed_ms,
      std::move(basis_owners), std::move(incoming_owner));
}

std::shared_ptr<StoredRefinedAcbMatch> restore_checkpoint_acb_match_record(
    const json::value& raw,
    const std::optional<std::string>& expected_session_configuration =
        std::nullopt,
    const std::optional<json::object>& expected_singular_request =
        std::nullopt) {
  const auto& object = as_object(raw, "checkpoint retained Acb match");
  require_exact_keys(object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "exact_lattice_identity", "exact_lattice_provenance_identity",
       "exact_lattice_canonical_witness", "basis_sources",
       "incoming_source", "basis_chart", "incoming_chart",
       "basis_point_exact", "incoming_point_exact",
       "physical_match_point_exact", "epsilon", "dimension",
       "relative_tolerance", "max_refinement_steps", "refined",
       "elapsed_ms"}, "checkpoint retained Acb match");
  if (required_string(object, "schema") !=
      "diffexp2-retained-acb-match-v2")
    throw std::invalid_argument(
        "unsupported retained Acb-match checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  const auto exact_lattice_identity = required_string(
      object, "exact_lattice_identity");
  const auto exact_lattice_provenance_identity = required_string(
      object, "exact_lattice_provenance_identity");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty() || exact_lattice_identity.empty() ||
      exact_lattice_provenance_identity.empty())
    throw std::invalid_argument(
        "checkpoint retained Acb match contains an empty identity");
  const auto basis_chart = required_string(object, "basis_chart");
  const auto incoming_chart = required_string(object, "incoming_chart");
  const auto basis_point = required_string(object, "basis_point_exact");
  const auto incoming_point = required_string(object,
                                               "incoming_point_exact");
  const auto physical_point = required_string(
      object, "physical_match_point_exact");
  if (basis_chart.empty() || incoming_chart.empty() || basis_point.empty() ||
      incoming_point.empty() || physical_point.empty())
    throw std::invalid_argument(
        "checkpoint Acb match lost chart or point provenance");
  const auto& raw_epsilon = as_object(object.at("epsilon"),
                                      "checkpoint Acb match epsilon window");
  require_exact_keys(raw_epsilon,
      {"min", "max", "required_complete_max"},
      "checkpoint Acb match epsilon window");
  const EpsilonWindow window{
      as_i32(raw_epsilon.at("min"), "checkpoint match epsilon minimum"),
      as_i32(raw_epsilon.at("max"), "checkpoint match epsilon maximum")};
  (void)window.width();
  const auto required_complete_max = as_i32(
      raw_epsilon.at("required_complete_max"),
      "checkpoint match required complete maximum");
  if (required_complete_max < window.min_power ||
      required_complete_max > window.complete_max)
    throw std::invalid_argument(
        "checkpoint Acb match requirement lies outside its work window");
  const auto dimension = as_u32(object.at("dimension"),
                                "checkpoint Acb match dimension");
  if (dimension == 0)
    throw std::invalid_argument("checkpoint Acb match dimension is zero");
  const auto relative_tolerance = required_string(
      object, "relative_tolerance");
  auto parsed_tolerance = Magnitude::decimal(relative_tolerance);
  if (!parsed_tolerance.is_finite())
    throw std::invalid_argument(
        "checkpoint Acb match tolerance is not finite");
  const auto max_refinement_steps = as_u64(
      object.at("max_refinement_steps"),
      "checkpoint Acb maximum refinement steps");
  if (max_refinement_steps > 32)
    throw std::invalid_argument(
        "checkpoint Acb maximum refinement steps exceeds 32");

  std::vector<json::object> basis_sources;
  const auto& raw_basis_sources = as_array(
      object.at("basis_sources"), "checkpoint Acb basis sources");
  if (raw_basis_sources.size() != dimension)
    throw std::invalid_argument(
        "checkpoint Acb match basis source count differs from its dimension");
  basis_sources.reserve(raw_basis_sources.size());
  for (std::size_t column = 0; column < raw_basis_sources.size(); ++column) {
    auto source = as_object(raw_basis_sources[column],
                            "checkpoint Acb basis source");
    validate_checkpoint_match_source(source, true, column);
    basis_sources.push_back(std::move(source));
  }
  auto incoming_source = as_object(object.at("incoming_source"),
                                   "checkpoint Acb incoming source");
  validate_checkpoint_match_source(incoming_source, false);

  const auto exact_lattice_witness_record = required_string(
      object, "exact_lattice_canonical_witness");
  if (exact_lattice_witness_record.empty())
    throw std::invalid_argument(
        "checkpoint Acb match lost its exact lattice witness");
  const auto raw_saturation_witness = json::parse(
      exact_lattice_witness_record);
  const auto& saturation_witness = as_object(
      raw_saturation_witness, "checkpoint Acb saturation witness");
  const auto saturation_witness_schema = required_string(
      saturation_witness, "schema");
  auto exact_lattice = [&]() -> ParsedExactEvaluatedLattice {
    if (saturation_witness_schema == kExactEvaluatedLatticeSchema)
      return parse_exact_evaluated_lattice(
          raw_saturation_witness, dimension, window,
          checkpoint_identity + ":checkpoint-restore");
    if (saturation_witness_schema == kNativeUnitSaturationProofSchema)
      return parse_native_unit_saturation_proof(
          raw_saturation_witness, dimension, window, basis_sources,
          basis_point, physical_point,
          checkpoint_identity + ":checkpoint-restore");
    if (saturation_witness_schema ==
        kNativeSingularSCCSaturationProofSchema) {
      if (!expected_session_configuration.has_value() ||
          !expected_singular_request.has_value())
        throw std::invalid_argument(
            "checkpoint singular-SCC saturation witness lacks its retained planned-match binding");
      return parse_native_singular_scc_saturation_proof(
          raw_saturation_witness, dimension, window, basis_sources,
          basis_point, physical_point, expected_session_configuration,
          expected_singular_request,
          checkpoint_identity + ":checkpoint-restore");
    }
    throw std::invalid_argument(
        "checkpoint Acb match has an unsupported saturation witness schema");
  }();
  if (exact_lattice.identity != exact_lattice_identity ||
      exact_lattice.canonical_witness != exact_lattice_witness_record)
    throw std::invalid_argument(
        "checkpoint exact lattice identity or canonical witness is inconsistent");

  json::array exact_binding_basis;
  for (const auto& source : basis_sources)
    exact_binding_basis.push_back(source);
  json::object exact_lattice_provenance{
      {"schema", "diffexp2-retained-exact-lattice-binding-v1"},
      {"witness_schema", exact_lattice.witness_schema},
      {"witness_identity", exact_lattice_identity},
      {"basis", std::move(exact_binding_basis)},
      {"basis_point_exact", basis_point},
      {"physical_match_point_exact", physical_point},
      {"epsilon", json::object{{"min", window.min_power},
                                {"max", window.complete_max}}}};
  if (json::serialize(canonical_json_value(exact_lattice_provenance)) !=
      exact_lattice_provenance_identity)
    throw std::invalid_argument(
        "checkpoint exact lattice provenance identity is inconsistent");

  json::array provenance_basis;
  for (const auto& source : basis_sources)
    provenance_basis.push_back(source);
  json::object provenance{
      {"schema", "diffexp2-native-refined-acb-local-match-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"basis", std::move(provenance_basis)},
      {"incoming", incoming_source},
      {"basis_point_exact", basis_point},
      {"incoming_point_exact", incoming_point},
      {"physical_match_point_exact", physical_point},
      {"epsilon", json::object{{"min", window.min_power},
                                {"max", window.complete_max},
                                {"required_complete_max",
                                 required_complete_max}}},
      {"exact_lattice_provenance_identity",
       exact_lattice_provenance_identity},
      {"refinement", json::object{{"relative_tolerance",
                                    relative_tolerance},
                                   {"max_steps",
                                    max_refinement_steps}}}};
  if (json::serialize(canonical_json_value(provenance)) !=
      provenance_identity)
    throw std::invalid_argument(
        "checkpoint retained Acb match provenance identity is inconsistent");

  const auto& raw_refined = as_object(object.at("refined"),
                                      "checkpoint refined Acb state");
  require_exact_keys(raw_refined,
      {"transformed_weights", "weights", "residual", "residual_history",
       "refinement_steps"}, "checkpoint refined Acb state");
  RefinedAcbLaurentMatch refined;
  refined.transformed_weights = parse_checkpoint_frame_vector<ComplexBall>(
      raw_refined.at("transformed_weights"), dimension,
      "checkpoint transformed Acb weights");
  refined.weights = parse_checkpoint_frame_vector<ComplexBall>(
      raw_refined.at("weights"), dimension, "checkpoint Acb weights");
  refined.residual = parse_checkpoint_frame_vector<ComplexBall>(
      raw_refined.at("residual"), dimension, "checkpoint Acb residual");
  refined.refinement_steps = static_cast<std::size_t>(as_u64(
      raw_refined.at("refinement_steps"),
      "checkpoint Acb refinement steps"));
  if (refined.refinement_steps > max_refinement_steps)
    throw std::invalid_argument(
        "checkpoint Acb refinement count exceeds its policy");
  for (const auto& raw_history : as_array(
           raw_refined.at("residual_history"),
           "checkpoint Acb residual history"))
    refined.residual_history.push_back(parse_checkpoint_acb_match_residual(
        raw_history, dimension, required_complete_max));
  if (refined.residual_history.size() != refined.refinement_steps + 1)
    throw std::invalid_argument(
        "checkpoint Acb residual history does not cover every refinement step");
  auto residual_min = refined.residual.front().min_power();
  auto residual_max = refined.residual.front().complete_max();
  for (const auto& row : refined.residual) {
    residual_min = std::min(residual_min, row.min_power());
    residual_max = std::min(residual_max, row.complete_max());
  }
  const auto& final_diagnostics = refined.residual_history.back();
  if (final_diagnostics.complete_window.min_power != residual_min ||
      final_diagnostics.complete_window.complete_max != residual_max)
    throw std::invalid_argument(
        "checkpoint final Acb residual diagnostics do not match the retained residual frame");
  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint Acb match elapsed time");
  return std::make_shared<StoredRefinedAcbMatch>(
      handle, checkpoint_identity, provenance_identity,
      exact_lattice_identity, exact_lattice_provenance_identity,
      exact_lattice_witness_record, exact_lattice.witness_schema,
      std::move(basis_sources),
      std::move(incoming_source), basis_chart, incoming_chart, basis_point,
      incoming_point, physical_point, window, required_complete_max,
      dimension, relative_tolerance,
      static_cast<std::size_t>(max_refinement_steps),
      std::move(exact_lattice.saturation), std::move(refined), elapsed_ms);
}

template <typename Scalar>
LocalSolution<Scalar> make_local_solution(
    const RecurrenceProblem<Scalar>& problem,
    AssembledResult<Scalar>&& assembled, LocalMetadata&& metadata) {
  LocalSolution<Scalar> solution;
  solution.chart = std::move(metadata.chart);
  solution.epsilon = {assembled.min_power, assembled.complete_max};
  solution.taylor_complete_max = problem.nmax;
  solution.dimension = problem.dimension;
  solution.prescriptions = std::move(metadata.prescriptions);
  solution.checkpoint_identity = std::move(metadata.checkpoint_identity);

  const auto sector_size = solution.sector_size();
  const auto sector_count = static_cast<std::size_t>(problem.log_max) + 1;
  if (sector_size > std::numeric_limits<std::size_t>::max() / sector_count ||
      assembled.coefficients.size() != sector_size * sector_count)
    throw std::invalid_argument(
        "assembled recurrence tensor cannot form the declared local sectors");
  auto cursor = assembled.coefficients.begin();
  for (std::uint32_t log = 0; log <= problem.log_max; ++log) {
    LocalSector<Scalar> sector;
    sector.a = metadata.a;
    sector.b = metadata.b;
    sector.log_power = log;
    sector.coefficients.reserve(sector_size);
    auto end = cursor + static_cast<std::ptrdiff_t>(sector_size);
    sector.coefficients.insert(
        sector.coefficients.end(), std::make_move_iterator(cursor),
        std::make_move_iterator(end));
    cursor = end;
    solution.sectors.push_back(std::move(sector));
  }
  return solution;
}

template <typename Scalar>
class PreparedChart final : public PreparedChartBase {
 public:
  PreparedChart(std::string handle, std::string key,
                std::string exact_identity, std::string signature,
                std::optional<std::string> geometry_record,
                std::optional<std::string> principal_matrix_record,
                std::optional<std::string> native_scc_capabilities,
                SCCCertificate scc,
                PreparedRecurrenceOperator<Scalar>&& prepared,
                slong precision_bits, std::vector<std::string> symbols,
                double prepare_parse_ms)
      : PreparedChartBase(std::move(handle), std::move(key),
                          std::move(exact_identity), std::move(signature),
                          std::move(geometry_record),
                          std::move(principal_matrix_record),
                          std::move(native_scc_capabilities), std::move(scc),
                          prepare_parse_ms),
        prepared_(std::move(prepared)), precision_bits_(precision_bits),
        symbols_(std::move(symbols)) {
    if constexpr (std::is_same_v<Scalar, Rational>) {
      try {
        exact_jordan_indicial_ =
            certify_exact_affine_jordan_operator(prepared_);
      } catch (const RecurrenceError& error) {
        exact_jordan_indicial_error_ = error.what();
      }
    }
  }

  json::object solve(const json::object& run, int output_digits) override {
    const auto parse_started = std::chrono::steady_clock::now();
    std::unique_ptr<AcbPrecisionLease> acb_lease;
    if constexpr (std::is_same_v<Scalar, ComplexBall>) {
      acb_lease = std::make_unique<AcbPrecisionLease>(precision_bits_);
      ComplexBall::set_precision(precision_bits_);
    }
    std::unique_lock<std::recursive_mutex> symbolic_lock;
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      symbolic_lock =
          std::unique_lock<std::recursive_mutex>(symbolic_run_mutex());
      SymbolicRational::configure(symbols_);
    }
    RecurrenceProblem<Scalar> problem;
    parse_run_state(run, prepared_, problem);
    const auto parse_ended = std::chrono::steady_clock::now();
    const auto run_parse_ms = std::chrono::duration<double, std::milli>(
        parse_ended - parse_started).count();
    auto result = run_prepared_problem(prepared_, problem, output_digits);
    const auto kernel_ms = result.at("elapsed_ms").as_double();
    const auto run_index = runs_.fetch_add(1) + 1;
    {
      std::lock_guard<std::mutex> lock(stats_mutex_);
      run_parse_ms_ += run_parse_ms;
      kernel_ms_ += kernel_ms;
    }
    result["persistent"] = json::object{
        {"run", run_index}, {"prepare_parse_ms", prepare_parse_ms_},
        {"run_parse_ms", run_parse_ms}, {"static_tensor_copies", 0},
        {"scc_components", scc_.component_count},
        {"scc_coupling_depth", scc_.coupling_depth}};
    return result;
  }

  NativeLocalRun<Scalar> solve_native(
      const json::object& run, const json::object& metadata_object) {
    return solve_native_impl(run, metadata_object, std::nullopt, false);
  }

  NativeLocalRun<Scalar> solve_native_with_source(
      const json::object& run, const json::object& metadata_object,
      SourceData<Scalar>&& source) {
    if (!run.at("source").is_null())
      throw std::invalid_argument(
          "native SCC source injection rejects caller-supplied source data");
    return solve_native_impl(
        run, metadata_object, std::move(source), false);
  }

  void record_native_local_success(
      const NativeLocalDiagnostics& diagnostics) {
    local_runs_.fetch_add(1);
    std::lock_guard<std::mutex> lock(stats_mutex_);
    local_run_parse_ms_ += diagnostics.parse_ms;
    local_kernel_ms_ += diagnostics.kernel_ms;
  }

  std::uint32_t dimension() const override { return prepared_.dimension; }
  std::int32_t frame_base() const override { return prepared_.frame_base; }
  std::uint32_t frame_width() const override { return prepared_.frame_width; }
  const char* d0_inverse_mode() const override {
    return prepared_.d0_inverse_scalar.has_value()
        ? "retained-scalar" : "retained-frame";
  }
  slong precision_bits() const { return precision_bits_; }
  bool has_identity_assembly() const {
    return prepared_.assembly_matrix.has_value() &&
           prepared_.assembly_matrix->identity;
  }
  bool has_regular_singleton_partition() const {
    if (prepared_.blocks.size() != prepared_.dimension) return false;
    std::vector<std::uint8_t> seen(prepared_.dimension, 0);
    for (const auto& block : prepared_.blocks) {
      if (block.columns.size() != 1 ||
          block.columns.front() >= prepared_.dimension ||
          seen[block.columns.front()])
        return false;
      seen[block.columns.front()] = 1;
    }
    return true;
  }
  std::size_t jordan_block_count() const {
    return prepared_.blocks.size();
  }
  bool jordan_partition_matches(
      const ExactJordanIndicialCertificate& certificate) const {
    if (certificate.dimension != prepared_.dimension ||
        certificate.blocks.size() != prepared_.blocks.size())
      return false;
    for (std::size_t index = 0; index < prepared_.blocks.size(); ++index)
      if (certificate.blocks[index].block_index != index ||
          certificate.blocks[index].columns !=
              prepared_.blocks[index].columns)
        return false;
    return true;
  }
  const std::optional<ExactJordanIndicialCertificate>&
  exact_jordan_indicial() const {
    return exact_jordan_indicial_;
  }
  const std::optional<std::string>& exact_jordan_indicial_error() const {
    return exact_jordan_indicial_error_;
  }
  ChartStats stats() const override {
    std::lock_guard<std::mutex> lock(stats_mutex_);
    return {runs_.load(), local_runs_.load(), prepare_parse_ms_,
            run_parse_ms_, kernel_ms_, local_run_parse_ms_,
            local_kernel_ms_};
  }

 private:
  NativeLocalRun<Scalar> solve_native_impl(
      const json::object& run, const json::object& metadata_object,
      std::optional<SourceData<Scalar>> native_source,
      bool attach_tail_model) {
    if (precision_bits_ < 64)
      throw std::invalid_argument(
          "native local solutions require at least 64 bits of Acb precision");
    // LocalSolution always carries numeric chart geometry and exact-tag
    // specializations even when its coefficient field is exact.  Lease the
    // output precision before parsing any such ball.
    AcbPrecisionLease acb_lease(precision_bits_);
    ComplexBall::set_precision(precision_bits_);
    std::unique_lock<std::recursive_mutex> symbolic_lock;
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      symbolic_lock =
          std::unique_lock<std::recursive_mutex>(symbolic_run_mutex());
      SymbolicRational::configure(symbols_);
    }

    const auto parse_started = std::chrono::steady_clock::now();
    RecurrenceProblem<Scalar> problem;
    parse_run_state(run, prepared_, problem);
    if (native_source.has_value()) {
      if (problem.source.has_value())
        throw std::invalid_argument(
            "native source injection found an unexpected parsed source");
      problem.source = std::move(*native_source);
    }
    if (!prepared_.assembly_matrix.has_value())
      throw std::invalid_argument(
          "local.solve requires a retained chart with native assembly");
    auto metadata = parse_local_metadata(metadata_object);
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      for (const auto* descriptor : {&metadata.a, &metadata.b})
        if (descriptor->domain == ExactDomain::SymbolicRational &&
            descriptor->symbols != symbols_)
          throw std::invalid_argument(
              "local exact-tag regulator field differs from its chart session");
    }
    verify_tag_binding(metadata.a, problem.a_target, "local a");
    verify_tag_binding(metadata.b, problem.b_target, "local b");
    const auto parse_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - parse_started).count();

    const auto kernel_started = std::chrono::steady_clock::now();
    auto recurrence = RecurrenceSolver<Scalar>(problem, prepared_).run();
    auto assembled = assemble_recurrence(prepared_, problem, recurrence);
    auto solution = make_local_solution(
        problem, std::move(assembled), std::move(metadata));
    validate_local_solution(solution, false);
    auto tail_model = unavailable_tail_model(
        attach_tail_model
            ? "symbolic coefficient locals do not support numeric certified tail bounds"
            : "tail model was not requested for this internal native block run");
    if (attach_tail_model) {
      if constexpr (std::is_same_v<Scalar, Rational> ||
                    std::is_same_v<Scalar, ComplexBall>)
        tail_model = prepare_regular_homogeneous_tail_model(
            prepared_, problem, solution, exact_identity_);
    }
    auto pseudo_hits = std::move(recurrence.hits);
    const auto kernel_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - kernel_started).count();
    return {std::move(solution), std::move(pseudo_hits),
            {recurrence.top_valid, parse_ms, kernel_ms},
            std::move(tail_model)};
  }

  std::shared_ptr<StoredLocalBase> solve_local(
      const std::string& local_handle, const json::object& run,
      const json::object& metadata_object) override {
    auto native = solve_native_impl(
        run, metadata_object, std::nullopt, true);
    auto local = make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
        local_handle, handle_, exact_identity_, std::move(native.solution),
        precision_bits_,
        std::move(native.pseudo_hits), native.diagnostics, std::nullopt,
        std::nullopt, nullptr,
        std::move(native.tail_model));
    record_native_local_success(native.diagnostics);
    return local;
  }
  PreparedRecurrenceOperator<Scalar> prepared_;
  std::optional<ExactJordanIndicialCertificate> exact_jordan_indicial_;
  std::optional<std::string> exact_jordan_indicial_error_;
  slong precision_bits_ = 256;
  std::vector<std::string> symbols_;
  std::atomic<std::uint64_t> runs_{0};
  std::atomic<std::uint64_t> local_runs_{0};
  mutable std::mutex stats_mutex_;
  double run_parse_ms_ = 0.0;
  double kernel_ms_ = 0.0;
  double local_run_parse_ms_ = 0.0;
  double local_kernel_ms_ = 0.0;
};

struct CompositeColumnSolveResult {
  std::shared_ptr<StoredLocalBase> local;
  json::array block_diagnostics;
  double elapsed_ms = 0.0;
};

class CompositeSCCChartBase {
 public:
  CompositeSCCChartBase(std::string handle, std::string key,
                        std::string exact_identity, std::string signature)
      : handle_(std::move(handle)), key_(std::move(key)),
        exact_identity_(std::move(exact_identity)),
        signature_(std::move(signature)) {}
  virtual ~CompositeSCCChartBase() = default;

  virtual json::object stats_json() const = 0;
  virtual CompositeColumnSolveResult solve_column(
      const std::string& local_handle, const json::object& request) = 0;
  virtual const char* column_execution_capability() const = 0;
  virtual const std::string& geometry_record() const = 0;
  virtual std::vector<std::shared_ptr<PreparedChartBase>>
  dependency_charts() const = 0;
  const std::string& handle() const { return handle_; }
  const std::string& key() const { return key_; }
  const std::string& exact_identity() const { return exact_identity_; }
  const std::string& signature() const { return signature_; }

 protected:
  std::string handle_;
  std::string key_;
  std::string exact_identity_;
  std::string signature_;
};

struct CompositeWorkContract {
  std::int32_t work_min = 0;
  std::int32_t requested_min = 0;
  std::int32_t requested_max = 0;
  std::int32_t work_complete_max = 0;
  std::uint32_t public_t_order = 0;
  std::uint32_t work_t_order = 0;
  std::uint32_t wolfram_coupling_depth = 0;
};

template <typename Scalar>
struct CompositeSCCBlock {
  std::uint32_t block = 0;
  std::vector<std::uint32_t> vertices;
  std::string source_handle;
  std::string principal_identity;
  bool regular = true;
  bool no_pseudo = false;
  std::optional<ExactJordanIndicialCertificate> exact_jordan_indicial;
  std::shared_ptr<PreparedChart<Scalar>> chart;
};

ExactJordanIndicialCertificate parse_exact_jordan_indicial_record(
    const json::value& raw, std::uint32_t expected_dimension) {
  const auto& object = as_object(
      raw, "exact affine-Jordan indicial certificate");
  require_exact_keys(object, {"schema", "dimension", "blocks"},
                     "exact affine-Jordan indicial certificate");
  if (required_string(object, "schema") !=
      "diffexp2-exact-affine-jordan-indicial-v1")
    throw std::invalid_argument(
        "unsupported exact affine-Jordan indicial certificate schema");
  ExactJordanIndicialCertificate certificate;
  certificate.dimension = as_u32(
      object.at("dimension"), "exact affine-Jordan dimension");
  if (certificate.dimension == 0 ||
      certificate.dimension != expected_dimension)
    throw std::invalid_argument(
        "exact affine-Jordan certificate dimension differs from its block");
  certificate.block_of_column.assign(
      certificate.dimension, std::numeric_limits<std::uint32_t>::max());
  certificate.position_in_block.assign(
      certificate.dimension, std::numeric_limits<std::uint32_t>::max());
  const auto& blocks = as_array(
      object.at("blocks"), "exact affine-Jordan blocks");
  if (blocks.empty())
    throw std::invalid_argument(
        "exact affine-Jordan certificate has no Jordan blocks");
  certificate.blocks.reserve(blocks.size());
  for (std::size_t index = 0; index < blocks.size(); ++index) {
    const auto& block = as_object(
        blocks[index], "exact affine-Jordan block");
    require_exact_keys(block, {"block", "columns", "a", "b"},
                       "exact affine-Jordan block");
    const auto block_index = as_u32(
        block.at("block"), "exact affine-Jordan block index");
    if (block_index != index)
      throw std::invalid_argument(
          "exact affine-Jordan blocks are not in deterministic order");
    std::vector<std::uint32_t> columns;
    for (const auto& raw_column : as_array(
             block.at("columns"), "exact affine-Jordan columns")) {
      const auto column = as_u32(
          raw_column, "exact affine-Jordan column");
      if (column >= certificate.dimension)
        throw std::invalid_argument(
            "exact affine-Jordan column is outside its dimension");
      columns.push_back(column);
    }
    if (columns.empty())
      throw std::invalid_argument(
          "exact affine-Jordan certificate contains an empty block");
    for (std::size_t position = 0; position < columns.size(); ++position) {
      const auto column = columns[position];
      if (certificate.block_of_column[column] !=
          std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument(
            "exact affine-Jordan blocks contain a duplicate column");
      certificate.block_of_column[column] = block_index;
      certificate.position_in_block[column] =
          static_cast<std::uint32_t>(position);
    }
    certificate.blocks.push_back(ExactJordanBlockCertificate{
        block_index, std::move(columns),
        ExactAffineIndicialRoot{
            Rational(required_string(block, "a")),
            Rational(required_string(block, "b"))}});
  }
  if (std::any_of(
          certificate.block_of_column.begin(),
          certificate.block_of_column.end(), [](std::uint32_t value) {
            return value == std::numeric_limits<std::uint32_t>::max();
          }))
    throw std::invalid_argument(
        "exact affine-Jordan blocks do not partition their dimension");
  return certificate;
}

bool same_exact_jordan_indicial(
    const ExactJordanIndicialCertificate& left,
    const ExactJordanIndicialCertificate& right) {
  if (left.dimension != right.dimension ||
      left.block_of_column != right.block_of_column ||
      left.position_in_block != right.position_in_block ||
      left.blocks.size() != right.blocks.size())
    return false;
  for (std::size_t index = 0; index < left.blocks.size(); ++index) {
    const auto& a = left.blocks[index];
    const auto& b = right.blocks[index];
    if (a.block_index != b.block_index || a.columns != b.columns ||
        !(a.root == b.root))
      return false;
  }
  return true;
}

struct CompositeCouplingIdentity {
  std::uint32_t source_vertex = 0;
  std::uint32_t target_vertex = 0;
  std::string exact_original_entry;
  std::string exact_theta_entry;
  bool proven_zero = false;
};

template <typename Scalar>
struct CompositeSCCCoupling {
  std::uint32_t source_block = 0;
  std::uint32_t target_block = 0;
  std::string producer_identity;
  PreparedSparseLocalMultiplierMatrix<Scalar> matrix;
  std::vector<CompositeCouplingIdentity> identities;
};

template <typename Scalar>
LocalSolution<Scalar> cap_composite_public_local(
    const LocalSolution<Scalar>& input, std::int32_t complete_max,
    std::uint32_t taylor_complete_max, const ChartGeometry& parent_chart,
    const std::vector<Prescription>& parent_prescriptions,
    std::string checkpoint_identity) {
  validate_local_solution(input, false);
  if (!input.error.empty())
    throw std::invalid_argument(
        "native SCC delivery cannot discard an error envelope");
  if (input.epsilon.complete_max < complete_max ||
      input.epsilon.min_power > complete_max)
    throw std::invalid_argument(
        "native SCC work state cannot deliver the requested epsilon maximum");
  if (input.taylor_complete_max < taylor_complete_max)
    throw std::invalid_argument(
        "native SCC work state cannot deliver the requested Taylor order");

  LocalSolution<Scalar> output;
  output.chart = parent_chart;
  output.epsilon = {input.epsilon.min_power, complete_max};
  output.taylor_complete_max = taylor_complete_max;
  output.dimension = input.dimension;
  output.prescriptions = parent_prescriptions;
  output.checkpoint_identity = std::move(checkpoint_identity);
  output.sectors.reserve(input.sectors.size());
  for (const auto& sector : input.sectors) {
    LocalSector<Scalar> capped;
    capped.a = sector.a;
    capped.b = sector.b;
    capped.log_power = sector.log_power;
    capped.coefficients.assign(output.sector_size(),
                               ScalarTraits<Scalar>::zero());
    for (std::int64_t power = output.epsilon.min_power;
         power <= output.epsilon.complete_max; ++power) {
      const auto input_epsilon = static_cast<std::size_t>(
          power - input.epsilon.min_power);
      const auto output_epsilon = static_cast<std::size_t>(
          power - output.epsilon.min_power);
      for (std::size_t n = 0; n < output.taylor_width(); ++n)
        for (std::uint32_t component = 0; component < output.dimension;
             ++component)
          capped.coefficients[local_algebra_detail::flat_index(
              output_epsilon, n, component, output.taylor_width(),
              output.dimension)] = sector.coefficients[
                  local_algebra_detail::flat_index(
                      input_epsilon, n, component, input.taylor_width(),
                      input.dimension)];
    }
    output.sectors.push_back(std::move(capped));
  }
  validate_local_solution(output, false);
  return output;
}

template <typename Scalar>
SourceData<Scalar> local_solution_source_data(
    const LocalSolution<Scalar>& source, std::uint32_t nmax,
    std::uint32_t log_max, std::int32_t frame_base,
    std::uint32_t frame_width) {
  validate_local_solution(source, false);
  if (source.dimension == 0 || !source.error.empty())
    throw std::invalid_argument(
        "native SCC source injection requires an uncertified nonempty local");
  const auto frame_top_i64 = static_cast<std::int64_t>(frame_base) +
      frame_width - 1;
  if (frame_top_i64 > std::numeric_limits<std::int32_t>::max() ||
      source.epsilon.min_power < frame_base ||
      source.epsilon.complete_max > frame_top_i64)
    throw std::invalid_argument(
        "native SCC source window lies outside the target retained frame");
  if (source.taylor_complete_max < nmax)
    throw std::invalid_argument(
        "native SCC source has insufficient Taylor order for its target");
  const auto taylor_points = static_cast<std::size_t>(nmax) + 1;
  const auto log_points = static_cast<std::size_t>(log_max) + 1;
  if (frame_width == 0 ||
      taylor_points > std::numeric_limits<std::size_t>::max() / log_points)
    throw std::overflow_error("native SCC source point tensor size overflow");
  const auto points = taylor_points * log_points;
  if (points > std::numeric_limits<std::size_t>::max() / source.dimension)
    throw std::overflow_error("native SCC source validity size overflow");
  const auto component_points = points * source.dimension;
  if (component_points >
      std::numeric_limits<std::size_t>::max() / frame_width)
    throw std::overflow_error("native SCC source tensor size overflow");
  SourceData<Scalar> data;
  data.frames.assign(component_points * frame_width,
                     ScalarTraits<Scalar>::zero());
  data.validity.assign(component_points, kCompleteInfinity);
  data.present.assign(points, 0);
  for (const auto& sector : source.sectors) {
    if (sector.log_power > log_max)
      throw std::invalid_argument(
          "native SCC source log sector exceeds the target run depth");
    for (std::uint32_t n = 0; n <= nmax; ++n) {
      const auto point = static_cast<std::size_t>(n) * log_points +
          sector.log_power;
      if (data.present[point])
        throw std::invalid_argument(
            "native SCC source contains duplicate exact log sectors");
      data.present[point] = 1;
      // Coefficients above CompleteMax are unknown, not certified zeros.  The
      // dense work-frame storage remains zero there while this finite validity
      // bound propagates independently for every target-block component.
      for (std::uint32_t component = 0;
           component < source.dimension; ++component) {
        const auto component_point = point * source.dimension + component;
        data.validity[component_point] = source.epsilon.complete_max;
        for (std::int64_t power = source.epsilon.min_power;
             power <= source.epsilon.complete_max; ++power) {
          const auto input_epsilon = static_cast<std::size_t>(
              power - source.epsilon.min_power);
          const auto output_epsilon = static_cast<std::size_t>(
              power - frame_base);
          data.frames[component_point * frame_width + output_epsilon] =
              sector.coefficients[local_algebra_detail::flat_index(
                  input_epsilon, n, component, source.taylor_width(),
                  source.dimension)];
        }
      }
    }
  }
  return data;
}

bool exact_nonnegative_integer(const Rational& value, bool include_zero) {
  if (value.sign() < 0 || (!include_zero && value.is_zero())) return false;
  return value.str().find('/') == std::string::npos;
}

std::uint32_t exact_log_ceiling(
    const ExactJordanIndicialCertificate& indicial,
    const Rational& a, const Rational& b, std::uint32_t base,
    bool include_zero_offset) {
  std::uint64_t result = base;
  for (const auto& block : indicial.blocks) {
    if (!(block.root.b == b)) continue;
    const auto offset = block.root.a - a;
    if (exact_nonnegative_integer(offset, include_zero_offset))
      result += block.size();
  }
  if (result > std::numeric_limits<std::uint32_t>::max())
    throw RecurrenceError(
        "E5", "derived exact Jordan log ceiling exceeds uint32 range");
  return static_cast<std::uint32_t>(result);
}

json::object exact_derived_run(
    const json::object& prototype, const PreparedChart<Rational>& chart,
    const Rational& a, const Rational& b, std::uint32_t base_log,
    bool homogeneous, std::optional<std::uint32_t> seed_component) {
  const auto& indicial = chart.exact_jordan_indicial();
  if (!indicial.has_value())
    throw RecurrenceError(
        "E5", "cannot derive a pseudo-compensation run without an exact Jordan certificate");
  const auto nmax = as_u32(prototype.at("nmax"),
                           "derived pseudo-compensation Taylor order");
  const auto dimension = chart.dimension();
  const auto frame_base = chart.frame_base();
  const auto frame_width = chart.frame_width();
  const auto frame_top_i64 = static_cast<std::int64_t>(frame_base) +
                             frame_width - 1;
  if (frame_top_i64 < std::numeric_limits<std::int32_t>::min() ||
      frame_top_i64 > std::numeric_limits<std::int32_t>::max())
    throw RecurrenceError(
        "E5", "derived pseudo-compensation frame exceeds int32 range");
  const auto frame_top = static_cast<std::int32_t>(frame_top_i64);

  std::uint32_t position = 0;
  if (homogeneous) {
    if (!seed_component.has_value() || *seed_component >= dimension)
      throw RecurrenceError(
          "E5", "derived pseudo-compensation seed component is out of range");
    position = indicial->position_in_block[*seed_component];
    base_log = std::max(base_log, position);
  } else if (seed_component.has_value()) {
    throw RecurrenceError(
        "E5", "derived particular run unexpectedly carries a seed component");
  }
  const auto log_max = exact_log_ceiling(
      *indicial, a, b, base_log, !homogeneous);

  json::array a_shifts;
  json::array schedule;
  a_shifts.reserve(static_cast<std::size_t>(nmax) + 1);
  schedule.reserve(static_cast<std::size_t>(nmax) + 1);
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    const auto a_n = a + Rational(std::to_string(n));
    a_shifts.push_back(json::string(a_n.str()));
    json::array row;
    row.reserve(indicial->blocks.size());
    for (const auto& block : indicial->blocks) {
      const auto d_a = a_n - block.root.a;
      const auto d_b = b - block.root.b;
      const auto step = singular_indicial_detail::classify_step(d_a, d_b);
      row.push_back(json::object{
          {"case", singular_indicial_detail::step_name(step)},
          {"da", d_a.str()}, {"db", d_b.str()}});
    }
    schedule.push_back(std::move(row));
  }

  json::array initial;
  json::array validity;
  const auto points = static_cast<std::size_t>(log_max) + 1;
  if (points > std::numeric_limits<std::size_t>::max() / dimension ||
      points * dimension >
          std::numeric_limits<std::size_t>::max() / frame_width)
    throw RecurrenceError(
        "E5", "derived pseudo-compensation seed tensor size overflows");
  initial.reserve(points * dimension * frame_width);
  validity.reserve(points * dimension);
  for (std::uint32_t log = 0; log <= log_max; ++log) {
    std::optional<std::uint32_t> expected_component;
    std::optional<std::size_t> expected_epsilon;
    if (homogeneous && log <= position) {
      const auto block_index = indicial->block_of_column[*seed_component];
      const auto& block = indicial->blocks[block_index];
      expected_component = block.columns[position - log];
      const auto epsilon_i64 = -static_cast<std::int64_t>(log) - frame_base;
      if (epsilon_i64 < 0 || epsilon_i64 >= frame_width)
        throw RecurrenceError(
            "E4", "derived canonical Jordan seed exceeds the retained lower epsilon frame",
            frame_base, -static_cast<std::int32_t>(log));
      expected_epsilon = static_cast<std::size_t>(epsilon_i64);
    }
    for (std::uint32_t component = 0; component < dimension; ++component) {
      for (std::uint32_t epsilon = 0; epsilon < frame_width; ++epsilon) {
        const bool unit = expected_component.has_value() &&
                          component == *expected_component &&
                          epsilon == *expected_epsilon;
        initial.push_back(unit ? "1" : "0");
      }
      validity.push_back(homogeneous ? json::value(frame_top)
                                     : json::value(nullptr));
    }
  }

  auto run = prototype;
  run["p"] = log_max;
  run["has_initial"] = homogeneous;
  run["adaptive_probe"] = false;
  run["a_target"] = a.str();
  run["b_target"] = b.str();
  run["a_shift_min"] = 0;
  run["a_shifts"] = std::move(a_shifts);
  run["schedule"] = std::move(schedule);
  run["initial"] = std::move(initial);
  run["initial_validity"] = std::move(validity);
  run["source"] = nullptr;
  run["return_u"] = false;
  return run;
}

json::object exact_derived_metadata(
    const json::object& prototype, const Rational& a, const Rational& b,
    std::uint32_t log_max, const std::string& suffix) {
  auto metadata = prototype;
  auto& tag = metadata.at("tag").as_object();
  tag["a"] = json::object{{"domain", "rational"},
                           {"canonical", a.str()}};
  tag["b"] = json::object{{"domain", "rational"},
                           {"canonical", b.str()}};
  tag["p"] = json::object{{"domain", "integer"},
                           {"canonical", std::to_string(log_max)}};
  metadata["checkpoint_identity"] =
      required_string(prototype, "checkpoint_identity") + suffix;
  return metadata;
}

struct ExactFormalKey {
  std::string t_power;
  std::int32_t epsilon_power = 0;
  std::uint32_t log_power = 0;
  std::uint32_t component = 0;

  friend bool operator<(const ExactFormalKey& left,
                        const ExactFormalKey& right) {
    return std::tie(left.t_power, left.epsilon_power, left.log_power,
                    left.component) <
           std::tie(right.t_power, right.epsilon_power, right.log_power,
                    right.component);
  }
};

// Expand an exact local slab in the formal basis
//
//   t^(a+n) eps^K Log(t)^p
//
// using t^(b eps) = Sum_j (b eps Log(t))^j/j!.  This is the
// Rational-domain CASE-P certificate: neither an Acb midpoint nor a
// tolerance can decide whether a polar coefficient cancels.  The optional
// t ceiling removes only the unmatched high-Taylor tail introduced when a
// target root starts n>0 orders above the source; every overlapping stored
// coefficient remains an exact proof obligation.
void add_exact_formal_below(
    const LocalSolution<Rational>& solution, std::int32_t exclusive_top,
    const std::optional<Rational>& maximum_t_power,
    std::map<ExactFormalKey, Rational>& coefficients) {
  validate_local_solution(solution, false);
  for (const auto& sector : solution.sectors) {
    if (sector.a.domain != ExactDomain::Rational ||
        sector.b.domain != ExactDomain::Rational)
      throw RecurrenceError(
          "E5", "exact pseudo-compensation certificate requires rational sector tags");
    const Rational a(sector.a.canonical);
    const Rational b(sector.b.canonical);
    Rational log_normalization(1);
    for (std::uint32_t divisor = 2; divisor <= sector.log_power; ++divisor)
      log_normalization =
          log_normalization / Rational(std::to_string(divisor));
    for (std::size_t n = 0; n < solution.taylor_width(); ++n) {
      const auto total_t_power = a + Rational(std::to_string(n));
      if (maximum_t_power.has_value() &&
          total_t_power > *maximum_t_power)
        continue;
      for (std::int64_t epsilon = solution.epsilon.min_power;
           epsilon <= solution.epsilon.complete_max; ++epsilon) {
        const auto epsilon_index = static_cast<std::size_t>(
            epsilon - solution.epsilon.min_power);
        const auto base_power_i64 = epsilon + sector.log_power;
        if (base_power_i64 >= exclusive_top) continue;
        if (base_power_i64 < std::numeric_limits<std::int32_t>::min())
          throw RecurrenceError(
              "E5", "exact pseudo-compensation epsilon power underflows int32");
        for (std::uint32_t component = 0;
             component < solution.dimension; ++component) {
          const auto& value = sector.coefficients[
              local_algebra_detail::flat_index(
                  epsilon_index, n, component, solution.taylor_width(),
                  solution.dimension)];
          if (value.is_zero()) continue;
          Rational exponential_factor(1);
          for (std::uint64_t j = 0;
               base_power_i64 + static_cast<std::int64_t>(j) <
                   exclusive_top;
               ++j) {
            if (j > std::numeric_limits<std::uint32_t>::max() -
                        sector.log_power)
              throw RecurrenceError(
                  "E5", "exact pseudo-compensation log degree overflows uint32");
            ExactFormalKey key{
                total_t_power.str(),
                static_cast<std::int32_t>(base_power_i64 +
                                          static_cast<std::int64_t>(j)),
                sector.log_power + static_cast<std::uint32_t>(j),
                component};
            auto found = coefficients.try_emplace(key, Rational(0)).first;
            found->second += value * log_normalization * exponential_factor;
            if (found->second.is_zero()) coefficients.erase(found);
            if (b.is_zero()) break;
            if (j == std::numeric_limits<std::uint32_t>::max())
              throw RecurrenceError(
                  "E5", "exact pseudo-compensation exponential order overflows uint32");
            exponential_factor = exponential_factor * b /
                Rational(std::to_string(j + 1));
          }
        }
      }
    }
  }
}

std::int32_t exact_formal_value_floor(
    const LocalSolution<Rational>& solution) {
  std::map<ExactFormalKey, Rational> coefficients;
  add_exact_formal_below(solution, 0, std::nullopt, coefficients);
  if (coefficients.empty()) return 0;
  return std::min_element(
      coefficients.begin(), coefficients.end(),
      [](const auto& left, const auto& right) {
        return left.first.epsilon_power < right.first.epsilon_power;
      })->first.epsilon_power;
}

void certify_exact_pseudo_value_floor(
    const LocalSolution<Rational>& solution, std::int32_t allowed_floor,
    const Rational& source_a, std::uint32_t source_taylor_max) {
  // Homogeneous targets pass floor 0.  Particulars pass the exact formal
  // floor already present in their source tag, so CASE-P may preserve a
  // genuine dimensional pole but may never manufacture a deeper one.
  const auto checked_top = std::min<std::int32_t>(0, allowed_floor);
  std::map<ExactFormalKey, Rational> coefficients;
  add_exact_formal_below(
      solution, checked_top,
      source_a + Rational(std::to_string(source_taylor_max)), coefficients);
  if (coefficients.empty()) return;
  const auto& witness = coefficients.begin()->first;
  std::string sector_witnesses;
  for (const auto& sector : solution.sectors) {
    auto single = solution;
    single.sectors = {sector};
    std::map<ExactFormalKey, Rational> part;
    add_exact_formal_below(
        single, checked_top,
        source_a + Rational(std::to_string(source_taylor_max)), part);
    const auto found = part.find(witness);
    if (found != part.end()) {
      if (!sector_witnesses.empty()) sector_witnesses += ";";
      sector_witnesses += "a=" + sector.a.canonical + ",b=" +
          sector.b.canonical + ",p=" +
          std::to_string(sector.log_power) + ":" + found->second.str();
    }
  }
  throw RecurrenceError(
      "E5", "exact CASE-P compensation leaves a value pole below the input floor at eps^" +
                std::to_string(witness.epsilon_power) + ", t_power=" +
                witness.t_power + ", log_power=" +
                std::to_string(witness.log_power) + ", component=" +
                std::to_string(witness.component) + ", coefficient=" +
                coefficients.begin()->second.str() + ", sectors=[" +
                sector_witnesses + "]");
}

std::vector<LocalSolution<Rational>> split_exact_rational_tags(
    const LocalSolution<Rational>& source, const std::string& identity) {
  validate_local_solution(source, false);
  std::map<std::pair<std::string, std::string>,
           std::vector<LocalSector<Rational>>> grouped;
  for (const auto& sector : source.sectors) {
    if (sector.a.domain != ExactDomain::Rational ||
        sector.b.domain != ExactDomain::Rational)
      throw RecurrenceError(
          "E5", "native SCC pseudo propagation requires exact rational source tags");
    grouped[{Rational(sector.a.canonical).str(),
             Rational(sector.b.canonical).str()}].push_back(sector);
  }
  std::vector<LocalSolution<Rational>> result;
  result.reserve(grouped.size());
  for (auto& [tag, sectors] : grouped) {
    auto group = source;
    group.sectors = std::move(sectors);
    group.checkpoint_identity = identity + ":tag:" + tag.first + ":" +
                                tag.second;
    result.push_back(canonicalize_identical_local_sectors(std::move(group)));
  }
  return result;
}

class ExactPseudoCompensator {
 public:
  ExactPseudoCompensator(PreparedChart<Rational>& chart,
                         std::string identity)
      : chart_(chart), identity_(std::move(identity)) {}

  NativeLocalRun<Rational> solve(
      const json::object& run, const json::object& metadata,
      std::optional<SourceData<Rational>> source,
      std::int32_t allowed_value_floor) {
    NativeLocalDiagnostics diagnostics;
    auto result = solve_impl(run, metadata, std::move(source),
                             allowed_value_floor, diagnostics);
    result.diagnostics = diagnostics;
    return result;
  }

 private:
  static void accumulate(NativeLocalDiagnostics& total,
                         const NativeLocalDiagnostics& current) {
    total.top_valid = std::min(total.top_valid, current.top_valid);
    total.parse_ms += current.parse_ms;
    total.kernel_ms += current.kernel_ms;
    total.pseudo_hits += current.pseudo_hits;
    total.pseudo_compensations += current.pseudo_compensations;
    total.max_pseudo_depth = std::max(
        total.max_pseudo_depth, current.max_pseudo_depth);
    total.pseudo_value_certified =
        total.pseudo_value_certified && current.pseudo_value_certified;
  }

  std::optional<PreparedRationalTaylorMultiplier<Rational>> polar_weight(
      const PseudoHit<Rational>& hit, std::size_t row,
      const LocalSolution<Rational>& target) const {
    if (row >= hit.gamma_frames.size() || row >= hit.gamma_validity.size())
      throw RecurrenceError(
          "E5", "CASE-P hit has inconsistent gamma frame dimensions");
    if (hit.gamma_validity[row] != kCompleteInfinity &&
        hit.gamma_validity[row] < -1)
      throw RecurrenceError(
          "E4", "CASE-P polar frame is not complete through eps^-1",
          chart_.frame_base(), hit.gamma_validity[row]);
    const auto& gamma = hit.gamma_frames[row];
    if (gamma.size() != chart_.frame_width())
      throw RecurrenceError(
          "E5", "CASE-P gamma frame differs from its retained chart width");
    std::optional<std::int32_t> minimum;
    for (std::size_t index = 0; index < gamma.size(); ++index) {
      const auto power = chart_.frame_base() +
                         static_cast<std::int32_t>(index);
      if (power >= 0) break;
      if (!gamma[index].is_zero()) {
        minimum = power;
        break;
      }
    }
    if (!minimum.has_value()) return std::nullopt;
    PreparedRationalTaylorMultiplier<Rational> multiplier;
    multiplier.epsilon_shift = *minimum;
    multiplier.center_pole_order = 0;
    multiplier.exact_identity = identity_ + ":casep:" +
        std::to_string(hit.n) + ":" + std::to_string(row);
    multiplier.kernels.assign(
        target.epsilon.width(),
        std::vector<Rational>(target.taylor_width(), Rational(0)));
    for (std::size_t kernel = 0; kernel < multiplier.kernels.size();
         ++kernel) {
      const auto power = static_cast<std::int64_t>(*minimum) +
                         static_cast<std::int64_t>(kernel);
      if (power >= 0) break;
      const auto gamma_index = power - chart_.frame_base();
      if (gamma_index < 0 ||
          gamma_index >= static_cast<std::int64_t>(gamma.size()))
        throw RecurrenceError(
            "E4", "CASE-P polar weight lies outside its retained gamma frame",
            chart_.frame_base(), static_cast<std::int32_t>(power));
      multiplier.kernels[kernel][0] =
          -gamma[static_cast<std::size_t>(gamma_index)];
    }
    return multiplier;
  }

  const LocalSolution<Rational>& homogeneous_target(
      std::uint32_t component, const json::object& prototype_run,
      const json::object& prototype_metadata,
      NativeLocalDiagnostics& diagnostics) {
    if (const auto found = homogeneous_cache_.find(component);
        found != homogeneous_cache_.end())
      return found->second;
    if (!active_components_.insert(component).second)
      throw RecurrenceError(
          "E5", "exact CASE-P compensation dependency is cyclic; the retained family ordering is not well founded");
    const auto& indicial = chart_.exact_jordan_indicial();
    if (!indicial.has_value() || component >= indicial->dimension)
      throw RecurrenceError(
          "E5", "CASE-P target component has no retained exact Jordan root");
    const auto& block =
        indicial->blocks[indicial->block_of_column[component]];
    auto run = exact_derived_run(
        prototype_run, chart_, block.root.a, block.root.b,
        indicial->position_in_block[component], true, component);
    auto metadata = exact_derived_metadata(
        prototype_metadata, block.root.a, block.root.b,
        as_u32(run.at("p"), "derived homogeneous log maximum"),
        ":casep-target:" + std::to_string(component));
    auto solved = solve_impl(run, metadata, std::nullopt, 0, diagnostics);
    active_components_.erase(component);
    auto [stored, inserted] = homogeneous_cache_.emplace(
        component, std::move(solved.solution));
    if (!inserted)
      throw std::logic_error("CASE-P homogeneous target cache insertion failed");
    return stored->second;
  }

  NativeLocalRun<Rational> solve_impl(
      const json::object& run, const json::object& metadata,
      std::optional<SourceData<Rational>> source,
      std::int32_t allowed_value_floor,
      NativeLocalDiagnostics& diagnostics) {
    auto raw = source.has_value()
        ? chart_.solve_native_with_source(
              run, metadata, std::move(*source))
        : chart_.solve_native(run, metadata);
    accumulate(diagnostics, raw.diagnostics);
    if (raw.pseudo_hits.empty()) return raw;

    const auto source_a = parse_scalar<Rational>(run.at("a_target"));
    const auto nmax = as_u32(run.at("nmax"),
                             "CASE-P certificate Taylor order");
    diagnostics.pseudo_hits += raw.pseudo_hits.size();
    std::vector<LocalSolution<Rational>> terms;
    terms.push_back(std::move(raw.solution));
    for (const auto& hit : raw.pseudo_hits) {
      diagnostics.max_pseudo_depth = std::max<std::uint32_t>(
          diagnostics.max_pseudo_depth,
          static_cast<std::uint32_t>(hit.columns.size()));
      if (hit.columns.size() != hit.gamma_frames.size() ||
          hit.columns.size() != hit.gamma_validity.size())
        throw RecurrenceError(
            "E5", "CASE-P hit target and gamma dimensions disagree");
      for (std::size_t row = 0; row < hit.columns.size(); ++row) {
        const auto& target = homogeneous_target(
            hit.columns[row], run, metadata, diagnostics);
        auto multiplier = polar_weight(hit, row, target);
        if (!multiplier.has_value()) continue;
        auto product = multiply_prepared_rational(
            target, *multiplier,
            identity_ + ":casep-product:" + std::to_string(hit.n) + ":" +
                std::to_string(hit.columns[row]));
        terms.push_back(std::move(product));
        ++diagnostics.pseudo_compensations;
      }
    }
    auto compensated = terms.size() == 1
        ? std::move(terms.front())
        : combine_local_solutions(
              terms, identity_ + ":casep-compensated");
    certify_exact_pseudo_value_floor(
        compensated, allowed_value_floor, source_a, nmax);
    raw.solution = std::move(compensated);
    raw.pseudo_hits.clear();
    return raw;
  }

  PreparedChart<Rational>& chart_;
  std::string identity_;
  std::map<std::uint32_t, LocalSolution<Rational>> homogeneous_cache_;
  std::set<std::uint32_t> active_components_;
};

template <typename Scalar>
class CompositeSCCChart final : public CompositeSCCChartBase {
 public:
  CompositeSCCChart(std::string handle, std::string key,
                    std::string exact_identity, std::string signature,
                    std::uint32_t dimension, SCCCertificate graph,
                    std::string exact_system_record,
                    std::string exact_theta_record,
                    std::string geometry_record,
                    RetainedCompositeGeometry retained_geometry,
                    CompositeWorkContract work,
                    std::vector<CompositeSCCBlock<Scalar>> blocks,
                    std::vector<CompositeSCCCoupling<Scalar>> couplings)
      : CompositeSCCChartBase(std::move(handle), std::move(key),
                              std::move(exact_identity),
                              std::move(signature)),
        dimension_(dimension), graph_(std::move(graph)),
        exact_system_record_(std::move(exact_system_record)),
        exact_theta_record_(std::move(exact_theta_record)),
        geometry_record_(std::move(geometry_record)),
        retained_geometry_(std::move(retained_geometry)), work_(work),
        blocks_(std::move(blocks)), couplings_(std::move(couplings)) {}

  CompositeColumnSolveResult solve_column(
      const std::string& local_handle,
      const json::object& request) override {
    if constexpr (std::is_same_v<Scalar, ComplexBall>) {
      return solve_regular_acb_column(local_handle, request);
    } else if constexpr (!std::is_same_v<Scalar, Rational>) {
      throw std::invalid_argument(
          "native SCC column execution supports exact rational columns and regular Acb columns only");
    } else {
      const auto started = std::chrono::steady_clock::now();
      const bool regular_execution = regular_block_column_ready();
      const bool regular_singular_execution =
          regular_singular_jordan_column_ready();
      if (!regular_execution && !regular_singular_execution)
        throw std::invalid_argument(
            "retained SCC chart does not satisfy a native exact-rational SCC column capability");
      const auto checkpoint_identity = required_string(
          request, "checkpoint_identity");
      const auto& seed_request = as_object(
          request.at("seed"), "native SCC seed request");
      const auto seed_block = as_u32(
          seed_request.at("block"), "native SCC seed block");
      if (seed_block >= blocks_.size())
        throw std::invalid_argument("native SCC seed block is out of range");

      std::vector<std::uint8_t> reachable(blocks_.size(), 0);
      reachable[seed_block] = 1;
      for (const auto block : graph_.topological_order) {
        if (!reachable[block]) continue;
        for (const auto [source, target] : graph_.condensation_edges)
          if (source == block) reachable[target] = 1;
      }
      std::vector<std::uint32_t> expected_targets;
      for (const auto block : graph_.topological_order)
        if (block != seed_block && reachable[block])
          expected_targets.push_back(block);
      const auto& target_requests = as_array(
          request.at("targets"), "native SCC target requests");
      if (target_requests.size() != expected_targets.size())
        throw std::invalid_argument(
            "native SCC targets must cover every reachable descendant exactly");
      for (std::size_t index = 0; index < target_requests.size(); ++index) {
        const auto& target = as_object(
            target_requests[index], "native SCC target request");
        if (as_u32(target.at("block"), "native SCC target block") !=
            expected_targets[index])
          throw std::invalid_argument(
              "native SCC targets are not in deterministic topological order");
      }

      std::vector<std::optional<LocalSolution<Scalar>>> state(blocks_.size());
      json::array diagnostics;
      NativeLocalDiagnostics aggregate;
      aggregate.top_valid = kCompleteInfinity;
      std::vector<std::unique_ptr<ExactPseudoCompensator>>
          pseudo_compensators(blocks_.size());
      const auto pseudo_compensator = [&](std::uint32_t block)
          -> ExactPseudoCompensator& {
        if (!pseudo_compensators[block])
          pseudo_compensators[block] =
              std::make_unique<ExactPseudoCompensator>(
                  *blocks_[block].chart,
                  checkpoint_identity + ":block:" +
                      std::to_string(block));
        return *pseudo_compensators[block];
      };

      const auto& seed_run = checked_column_run(
          seed_request, seed_block, true, nullptr,
          regular_singular_execution);
      const auto seed_local_component = seed_component_from_run(
          seed_run, seed_block, regular_singular_execution);
      const auto basis_index =
          blocks_[seed_block].vertices[seed_local_component];
      auto seed_native = regular_singular_execution
          ? pseudo_compensator(seed_block).solve(
                seed_run,
                as_object(seed_request.at("metadata"),
                          "native SCC seed metadata"),
                std::nullopt, 0)
          : blocks_[seed_block].chart->solve_native(
                seed_run, as_object(seed_request.at("metadata"),
                                    "native SCC seed metadata"));
      validate_block_result(seed_native, seed_block, true,
                            regular_singular_execution);
      blocks_[seed_block].chart->record_native_local_success(
          seed_native.diagnostics);
      accumulate_diagnostics(aggregate, seed_native.diagnostics);
      diagnostics.push_back(block_diagnostic(
          seed_block, "seed", {}, nullptr, seed_native));
      state[seed_block] = std::move(seed_native.solution);

      for (std::size_t target_index = 0;
           target_index < target_requests.size(); ++target_index) {
        const auto target_block = expected_targets[target_index];
        const auto& target_request = as_object(
            target_requests[target_index], "native SCC target request");
        std::vector<LocalSolution<Scalar>> incoming;
        std::vector<std::uint32_t> predecessors;
        for (const auto& coupling : couplings_) {
          if (coupling.target_block != target_block ||
              !state[coupling.source_block].has_value())
            continue;
          validate_column_coupling(
              coupling, regular_singular_execution);
          auto contribution = apply_prepared_sparse_local_matrix(
              coupling.matrix, *state[coupling.source_block],
              checkpoint_identity + ":source:" +
                  std::to_string(coupling.source_block) + ":" +
                  std::to_string(target_block));
          if (!contribution.has_value())
            throw std::logic_error(
                "an exact nonzero SCC edge produced no structural source");
          predecessors.push_back(coupling.source_block);
          incoming.push_back(std::move(*contribution));
        }
        if (incoming.empty())
          throw std::invalid_argument(
              "reachable native SCC target has no available predecessor source");
        auto source = incoming.size() == 1
            ? std::move(incoming.front())
            : combine_local_solutions(
                  incoming, checkpoint_identity + ":combined-source:" +
                                std::to_string(target_block));
        // A signed-shift halo is a property of the complete target source.
        // Restricting predecessors separately would reject exact below-frame
        // pieces which cancel only after all incoming block edges are summed.
        source = restrict_local_epsilon_frame_strict_lower(
            source, work_.work_min, work_.work_complete_max,
            checkpoint_identity + ":source-frame:" +
                std::to_string(target_block));
        require_work_local(source, "combined coupling source");
        if (!regular_singular_execution) {
          const auto& target_run = checked_column_run(
              target_request, target_block, false, &source, false);
          require_source_tag_matches_run(source, target_run);
          auto source_data = local_solution_source_data(
              source, as_u32(target_run.at("nmax"), "target nmax"),
              as_u32(target_run.at("p"), "target log maximum"),
              blocks_[target_block].chart->frame_base(),
              blocks_[target_block].chart->frame_width());
          auto target_native =
              blocks_[target_block].chart->solve_native_with_source(
                  target_run,
                  as_object(target_request.at("metadata"),
                            "native SCC target metadata"),
                  std::move(source_data));
          validate_block_result(target_native, target_block, false, false);
          blocks_[target_block].chart->record_native_local_success(
              target_native.diagnostics);
          accumulate_diagnostics(aggregate, target_native.diagnostics);
          diagnostics.push_back(block_diagnostic(
              target_block, "particular", predecessors, &source,
              target_native));
          state[target_block] = std::move(target_native.solution);
          continue;
        }

        // CASE-P compensation introduces homogeneous target-root sectors.
        // They are linearly independent exact tags and must remain separate
        // through every later SCC edge.  Solve one exact particular per tag,
        // then recombine; collapsing them into a single recurrence would
        // silently apply the wrong affine schedule to all but one sector.
        auto groups = split_exact_rational_tags(
            source, checkpoint_identity + ":source-groups:" +
                        std::to_string(target_block));
        const auto& submitted_run = as_object(
            target_request.at("run"), "native SCC target run");
        const Rational submitted_a =
            parse_scalar<Rational>(submitted_run.at("a_target"));
        const Rational submitted_b =
            parse_scalar<Rational>(submitted_run.at("b_target"));
        bool submitted_used = false;
        std::vector<LocalSolution<Rational>> target_parts;
        target_parts.reserve(groups.size());
        for (std::size_t group_index = 0; group_index < groups.size();
             ++group_index) {
          auto& group = groups[group_index];
          if (group.sectors.empty())
            throw std::logic_error("exact SCC source tag group is empty");
          const Rational group_a(group.sectors.front().a.canonical);
          const Rational group_b(group.sectors.front().b.canonical);
          const bool use_submitted = group_a == submitted_a &&
                                     group_b == submitted_b;
          json::object derived_entry;
          const json::object* entry = &target_request;
          if (!use_submitted) {
            std::uint32_t source_log = 0;
            for (const auto& sector : group.sectors)
              source_log = std::max(source_log, sector.log_power);
            auto run = exact_derived_run(
                submitted_run, *blocks_[target_block].chart,
                group_a, group_b, source_log, false, std::nullopt);
            auto metadata = exact_derived_metadata(
                as_object(target_request.at("metadata"),
                          "native SCC target metadata"),
                group_a, group_b,
                as_u32(run.at("p"), "derived particular log maximum"),
                ":derived-tag:" + std::to_string(group_index));
            derived_entry = json::object{
                {"block", target_block}, {"run", std::move(run)},
                {"metadata", std::move(metadata)}};
            entry = &derived_entry;
          } else {
            if (submitted_used)
              throw std::logic_error(
                  "exact SCC source contains duplicate submitted tag groups");
            submitted_used = true;
          }

          const auto& target_run = checked_column_run(
              *entry, target_block, false, &group, true);
          require_source_tag_matches_run(group, target_run);
          const auto allowed_floor = exact_formal_value_floor(group);
          auto source_data = local_solution_source_data(
              group, as_u32(target_run.at("nmax"), "target nmax"),
              as_u32(target_run.at("p"), "target log maximum"),
              blocks_[target_block].chart->frame_base(),
              blocks_[target_block].chart->frame_width());
          auto target_native = pseudo_compensator(target_block).solve(
              target_run,
              as_object(entry->at("metadata"),
                        "native SCC target metadata"),
              std::move(source_data), allowed_floor);
          validate_block_result(target_native, target_block, false, true);
          blocks_[target_block].chart->record_native_local_success(
              target_native.diagnostics);
          accumulate_diagnostics(aggregate, target_native.diagnostics);
          auto diagnostic = block_diagnostic(
              target_block, "particular-tag", predecessors, &group,
              target_native);
          diagnostic["source_a"] = group_a.str();
          diagnostic["source_b"] = group_b.str();
          diagnostics.push_back(std::move(diagnostic));
          target_parts.push_back(std::move(target_native.solution));
        }
        if (!submitted_used)
          throw std::invalid_argument(
              "native SCC submitted target tag is absent from the exact propagated source groups");
        auto target_state = target_parts.size() == 1
            ? std::move(target_parts.front())
            : combine_local_solutions(
                  target_parts,
                  checkpoint_identity + ":target-tag-sum:" +
                      std::to_string(target_block));
        require_work_local(target_state, "combined target tag particulars");
        state[target_block] = std::move(target_state);
      }

      std::vector<LocalSolution<Scalar>> embedded;
      for (std::uint32_t block = 0; block < state.size(); ++block) {
        if (!state[block].has_value()) continue;
        embedded.push_back(local_algebra_detail::embedded_components(
            *state[block], blocks_[block].vertices, dimension_));
      }
      if (embedded.empty())
        throw std::logic_error("native SCC column produced no block state");
      auto parent = embedded.size() == 1
          ? std::move(embedded.front())
          : combine_local_solutions(
                embedded, checkpoint_identity + ":work-parent");
      require_work_local(parent, "combined parent work state");
      parent = cap_composite_public_local(
          parent, work_.requested_max, work_.public_t_order,
          retained_geometry_.chart, retained_geometry_.prescriptions,
          checkpoint_identity);
      validate_local_solution(parent, false);
      const auto scalar_execution = scalar_block_shape();
      json::object column_identity_record{
          {"schema", regular_singular_execution
               ? (scalar_execution
                     ? "diffexp2-native-scc-regular-singular-scalar-column-v1"
                     : "diffexp2-native-scc-regular-singular-jordan-column-v2")
               : (scalar_execution
                     ? "diffexp2-native-scc-column-v1"
                     : "diffexp2-native-scc-column-v2")},
          {"scc_exact_identity", exact_identity_},
          {"basis_index", basis_index},
          {"seed", seed_request},
          {"targets", target_requests}};
      if (!scalar_execution)
        column_identity_record["seed_local_component"] =
            seed_local_component;
      column_identity_record["pseudo_compensation"] =
          aggregate.pseudo_hits == 0
              ? "none"
              : "exact-rational-derived-jordan-targets-v1";
      SCCColumnProvenance column_provenance{
          handle_, exact_identity_, seed_block,
          basis_index,
          json::serialize(canonical_json_value(column_identity_record))};
      auto local = make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
          local_handle, handle_, exact_identity_, std::move(parent),
          blocks_[seed_block].chart->precision_bits(),
          std::vector<PseudoHit<Scalar>>{}, aggregate,
          std::move(column_provenance));
      const auto elapsed_ms = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      column_solves_.fetch_add(1);
      {
        std::lock_guard<std::mutex> lock(column_stats_mutex_);
        column_solve_ms_ += elapsed_ms;
      }
      return {std::move(local), std::move(diagnostics), elapsed_ms};
    }
  }

  double column_solve_ms() const {
    std::lock_guard<std::mutex> lock(column_stats_mutex_);
    return column_solve_ms_;
  }

  const char* column_execution_capability() const override {
    if (regular_singular_jordan_column_ready())
      return regular_singular_column_capability();
    if (!regular_block_column_ready())
      return "unsupported-native-scc-column";
    return regular_column_capability();
  }

  const std::string& geometry_record() const override {
    return geometry_record_;
  }

  std::vector<std::shared_ptr<PreparedChartBase>>
  dependency_charts() const override {
    std::vector<std::shared_ptr<PreparedChartBase>> result;
    result.reserve(blocks_.size());
    for (const auto& block : blocks_) result.push_back(block.chart);
    return result;
  }

  json::object stats_json() const override {
    std::size_t active_entries = 0, proven_zero_entries = 0;
    std::optional<std::int32_t> min_coupling_shift;
    std::optional<std::int32_t> max_coupling_shift;
    json::array block_handles;
    block_handles.reserve(blocks_.size());
    for (const auto& block : blocks_) {
      json::object block_record{
          {"block", block.block}, {"chart", block.source_handle},
          {"dimension", block.chart->dimension()},
          {"regular", block.regular},
          {"no_pseudo", block.no_pseudo},
          {"principal_identity", block.principal_identity}};
      if (const auto& indicial = block.exact_jordan_indicial;
          indicial.has_value()) {
        json::array indicial_blocks;
        std::uint32_t max_jordan_size = 0;
        for (const auto& spectral_block : indicial->blocks) {
          max_jordan_size = std::max(max_jordan_size,
                                     spectral_block.size());
          indicial_blocks.push_back(json::object{
              {"block", spectral_block.block_index},
              {"columns", encode_indices(spectral_block.columns)},
              {"jordan_size", spectral_block.size()},
              {"a", spectral_block.root.a.str()},
              {"b", spectral_block.root.b.str()}});
        }
        block_record["exact_affine_jordan_indicial"] = json::object{
            {"dimension", indicial->dimension},
            {"blocks", std::move(indicial_blocks)},
            {"max_jordan_size", max_jordan_size}};
        // Preserve the scalar-v1 diagnostic field while deriving it from the
        // same complete exact proof used by multidimensional admission.
        if (indicial->dimension == 1 && indicial->blocks.size() == 1)
          block_record["affine_indicial_root"] = json::object{
              {"a", indicial->blocks.front().root.a.str()},
              {"b", indicial->blocks.front().root.b.str()}};
      } else if (block.chart->exact_jordan_indicial_error().has_value()) {
        block_record["exact_affine_jordan_indicial_error"] =
            *block.chart->exact_jordan_indicial_error();
      }
      block_handles.push_back(std::move(block_record));
    }
    for (const auto& coupling : couplings_)
      for (std::size_t index = 0; index < coupling.identities.size(); ++index) {
        const auto& identity = coupling.identities[index];
        identity.proven_zero ? ++proven_zero_entries : ++active_entries;
        if (!identity.proven_zero) {
          const auto shift = coupling.matrix.entries[index]
                                 .multiplier.epsilon_shift;
          min_coupling_shift = min_coupling_shift.has_value()
              ? std::min(*min_coupling_shift, shift) : shift;
          max_coupling_shift = max_coupling_shift.has_value()
              ? std::max(*max_coupling_shift, shift) : shift;
        }
      }
    const auto regular_ready = regular_block_column_ready();
    const auto regular_singular_ready =
        regular_singular_jordan_column_ready();
    const auto execution_ready = regular_ready || regular_singular_ready;
    const auto scalar_shape = scalar_block_shape();
    const auto scalar_ready = scalar_column_ready();
    json::object result{
        {"scc", handle_}, {"key", key_}, {"identity", exact_identity_},
        {"dimension", dimension_}, {"blocks", blocks_.size()},
        {"coupling_groups", couplings_.size()},
        {"coupling_entries", active_entries + proven_zero_entries},
        {"active_coupling_entries", active_entries},
        {"proven_zero_coupling_entries", proven_zero_entries},
        {"frame_base", work_.work_min},
        {"frame_width", static_cast<std::uint32_t>(
             static_cast<std::int64_t>(work_.work_complete_max) -
             work_.work_min + 1)},
        {"requested_min", work_.requested_min},
        {"requested_max", work_.requested_max},
        {"work_complete_max", work_.work_complete_max},
        {"public_t_order", work_.public_t_order},
        {"work_t_order", work_.work_t_order},
        {"wolfram_coupling_depth", work_.wolfram_coupling_depth},
        {"native_coupling_depth", graph_.coupling_depth},
        {"min_coupling_shift", min_coupling_shift.has_value()
             ? json::value(*min_coupling_shift) : json::value(nullptr)},
        {"max_coupling_shift", max_coupling_shift.has_value()
             ? json::value(*max_coupling_shift) : json::value(nullptr)},
        {"execution_mode", "BlockSequentialStrict"},
        {"execution_implemented", execution_ready},
        {"execution_scope", regular_singular_ready
             ? regular_singular_column_capability()
             : (regular_ready
                   ? regular_column_capability()
                   : "unsupported")},
        {"general_scc_execution", false},
        {"scalar_block_dag_column_execution", scalar_ready},
        {"regular_singular_scalar_block_dag_column_execution",
         regular_singular_ready && scalar_shape},
        {"regular_singular_jordan_block_dag_column_execution",
         regular_singular_ready},
        {"scc_column_solves", column_solves_.load()},
        {"scc_column_solve_ms", column_solve_ms()},
        {"capability_evidence", json::object{
             {"identity_v", "native-retained-assembly"},
             {"regular", regular_ready
                  ? "collision-bound-producer-certificate"
                  : "not-required-by-selected-scope"},
             {"regular_or_regular_singular",
              "collision-bound-producer-certificate"},
             {"identity_gauge", "collision-bound-producer-certificate"},
             {"no_pseudo", regular_singular_ready
                  ? (std::is_same_v<Scalar, ComplexBall>
                        ? "producer-proven-and-exact-schedule-revalidated-no-case-p"
                        : "producer-provenance-only-execution-revalidated-by-exact-schedule-certificate")
                  : "collision-bound-producer-certificate"},
             {"jordan_indicial",
              regular_singular_ready
                  ? "retained-exact-rational-full-matrix-certificate"
                  : "not-required-by-selected-scope"},
             {"pseudo_schedule_execution",
              regular_singular_ready
                  ? (std::is_same_v<Scalar, ComplexBall>
                        ? "exact-rational-certificate-case-p-rejected-for-acb"
                        : "exact-rational-joint-compensation-and-formal-overlap-certificate")
                  : "not-required-by-selected-scope"},
             {"resonance_schedule",
              regular_singular_ready
                  ? "retained-affine-jordan-verified-exact-captured-run"
                  : "retained-affine-root-verified-exact-captured-run"}}},
        {"execution_must_revalidate_producer_capabilities", true},
        {"block_charts", std::move(block_handles)}};
    if (!scalar_shape)
      result["regular_block_dag_column_execution"] = regular_ready;
    return result;
  }

 private:
  CompositeColumnSolveResult solve_regular_acb_column(
      const std::string& local_handle, const json::object& request) {
    static_assert(std::is_same_v<Scalar, ComplexBall>);
    const bool regular_singular_execution =
        regular_singular_jordan_column_ready();
    if (!regular_block_column_ready() && !regular_singular_execution)
      throw std::invalid_argument(
          "retained Acb SCC chart does not satisfy a strict regular or exact affine-Jordan no-CASE-P block-DAG column capability");

    // Coupling products and block recombination happen between the individual
    // PreparedChart solves, so the complete composite operation must retain
    // the same Arb precision lease rather than relying on each nested solve's
    // shorter lease.
    AcbPrecisionLease acb_lease(blocks_.front().chart->precision_bits());
    ComplexBall::set_precision(blocks_.front().chart->precision_bits());
    const auto started = std::chrono::steady_clock::now();
    const auto checkpoint_identity = required_string(
        request, "checkpoint_identity");
    const auto& seed_request = as_object(
        request.at("seed"), "native Acb SCC seed request");
    const auto seed_block = as_u32(
        seed_request.at("block"), "native Acb SCC seed block");
    if (seed_block >= blocks_.size())
      throw std::invalid_argument("native Acb SCC seed block is out of range");

    std::vector<std::uint8_t> reachable(blocks_.size(), 0);
    reachable[seed_block] = 1;
    for (const auto block : graph_.topological_order) {
      if (!reachable[block]) continue;
      for (const auto [source, target] : graph_.condensation_edges)
        if (source == block) reachable[target] = 1;
    }
    std::vector<std::uint32_t> expected_targets;
    for (const auto block : graph_.topological_order)
      if (block != seed_block && reachable[block])
        expected_targets.push_back(block);
    const auto& target_requests = as_array(
        request.at("targets"), "native Acb SCC target requests");
    if (target_requests.size() != expected_targets.size())
      throw std::invalid_argument(
          "native Acb SCC targets must cover every reachable descendant exactly");
    for (std::size_t index = 0; index < target_requests.size(); ++index) {
      const auto& target = as_object(
          target_requests[index], "native Acb SCC target request");
      if (as_u32(target.at("block"), "native Acb SCC target block") !=
          expected_targets[index])
        throw std::invalid_argument(
            "native Acb SCC targets are not in deterministic topological order");
    }

    std::vector<std::optional<LocalSolution<Scalar>>> state(blocks_.size());
    json::array diagnostics;
    NativeLocalDiagnostics aggregate;
    aggregate.top_valid = kCompleteInfinity;

    const auto& seed_run = checked_column_run(
        seed_request, seed_block, true, nullptr,
        regular_singular_execution);
    const auto seed_local_component = seed_component_from_run(
        seed_run, seed_block, regular_singular_execution);
    const auto basis_index =
        blocks_[seed_block].vertices[seed_local_component];
    auto seed_native = blocks_[seed_block].chart->solve_native(
        seed_run, as_object(seed_request.at("metadata"),
                            "native Acb SCC seed metadata"));
    validate_block_result(seed_native, seed_block, true,
                          regular_singular_execution);
    blocks_[seed_block].chart->record_native_local_success(
        seed_native.diagnostics);
    accumulate_diagnostics(aggregate, seed_native.diagnostics);
    diagnostics.push_back(block_diagnostic(
        seed_block, "seed", {}, nullptr, seed_native));
    state[seed_block] = std::move(seed_native.solution);

    for (std::size_t target_index = 0;
         target_index < target_requests.size(); ++target_index) {
      const auto target_block = expected_targets[target_index];
      const auto& target_request = as_object(
          target_requests[target_index], "native Acb SCC target request");
      std::vector<LocalSolution<Scalar>> incoming;
      std::vector<std::uint32_t> predecessors;
      for (const auto& coupling : couplings_) {
        if (coupling.target_block != target_block ||
            !state[coupling.source_block].has_value())
          continue;
        validate_column_coupling(coupling, regular_singular_execution);
        auto contribution = apply_prepared_sparse_local_matrix(
            coupling.matrix, *state[coupling.source_block],
            checkpoint_identity + ":source:" +
                std::to_string(coupling.source_block) + ":" +
                std::to_string(target_block));
        if (!contribution.has_value())
          throw std::logic_error(
              "a certified nonzero Acb SCC edge produced no structural source");
        predecessors.push_back(coupling.source_block);
        incoming.push_back(std::move(*contribution));
      }
      if (incoming.empty())
        throw std::invalid_argument(
            "reachable native Acb SCC target has no available predecessor source");
      auto source = incoming.size() == 1
          ? std::move(incoming.front())
          : combine_local_solutions(
                incoming, checkpoint_identity + ":combined-source:" +
                              std::to_string(target_block));
      source = restrict_local_epsilon_frame_strict_lower(
          source, work_.work_min, work_.work_complete_max,
          checkpoint_identity + ":source-frame:" +
              std::to_string(target_block));
      require_work_local(source, "combined Acb coupling source");

      const auto& target_run = checked_column_run(
          target_request, target_block, false, &source,
          regular_singular_execution);
      require_source_tag_matches_run(source, target_run);
      auto source_data = local_solution_source_data(
          source, as_u32(target_run.at("nmax"), "target nmax"),
          as_u32(target_run.at("p"), "target log maximum"),
          blocks_[target_block].chart->frame_base(),
          blocks_[target_block].chart->frame_width());
      auto target_native =
          blocks_[target_block].chart->solve_native_with_source(
              target_run,
              as_object(target_request.at("metadata"),
                        "native Acb SCC target metadata"),
              std::move(source_data));
      validate_block_result(target_native, target_block, false,
                            regular_singular_execution);
      blocks_[target_block].chart->record_native_local_success(
          target_native.diagnostics);
      accumulate_diagnostics(aggregate, target_native.diagnostics);
      diagnostics.push_back(block_diagnostic(
          target_block, "particular", predecessors, &source,
          target_native));
      state[target_block] = std::move(target_native.solution);
    }

    std::vector<LocalSolution<Scalar>> embedded;
    for (std::uint32_t block = 0; block < state.size(); ++block) {
      if (!state[block].has_value()) continue;
      embedded.push_back(local_algebra_detail::embedded_components(
          *state[block], blocks_[block].vertices, dimension_));
    }
    if (embedded.empty())
      throw std::logic_error("native Acb SCC column produced no block state");
    auto parent = embedded.size() == 1
        ? std::move(embedded.front())
        : combine_local_solutions(
              embedded, checkpoint_identity + ":work-parent");
    require_work_local(parent, "combined Acb parent work state");
    parent = cap_composite_public_local(
        parent, work_.requested_max, work_.public_t_order,
        retained_geometry_.chart, retained_geometry_.prescriptions,
        checkpoint_identity);
    validate_local_solution(parent, false);

    const auto scalar_execution = scalar_block_shape();
    json::object column_identity_record{
        {"schema", regular_singular_execution
             ? (scalar_execution
                   ? "diffexp2-native-scc-acb-regular-singular-scalar-column-v1"
                   : "diffexp2-native-scc-acb-regular-singular-jordan-column-v1")
             : (scalar_execution
                   ? "diffexp2-native-scc-acb-regular-scalar-column-v1"
                   : "diffexp2-native-scc-acb-regular-column-v2")},
        {"scc_exact_identity", exact_identity_},
        {"basis_index", basis_index},
        {"seed", seed_request},
        {"targets", target_requests},
        {"pseudo_compensation", "none"}};
    if (!scalar_execution)
      column_identity_record["seed_local_component"] =
          seed_local_component;
    SCCColumnProvenance column_provenance{
        handle_, exact_identity_, seed_block, basis_index,
        json::serialize(canonical_json_value(column_identity_record))};
    auto local = make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
        local_handle, handle_, exact_identity_, std::move(parent),
        blocks_[seed_block].chart->precision_bits(),
        std::vector<PseudoHit<Scalar>>{}, aggregate,
        std::move(column_provenance));
    const auto elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    column_solves_.fetch_add(1);
    {
      std::lock_guard<std::mutex> lock(column_stats_mutex_);
      column_solve_ms_ += elapsed_ms;
    }
    return {std::move(local), std::move(diagnostics), elapsed_ms};
  }

  bool regular_block_column_ready() const {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      return false;
    } else {
      if (blocks_.size() < 2 || work_.work_min > 0 ||
          work_.requested_max < 0 || work_.work_complete_max < 0 ||
          std::any_of(blocks_.begin(), blocks_.end(), [](const auto& block) {
            return !block.regular || block.vertices.empty() ||
                   block.chart->dimension() != block.vertices.size() ||
                   !block.chart->has_identity_assembly() ||
                   !block.chart->has_regular_singleton_partition();
          }))
        return false;
      return std::all_of(
          couplings_.begin(), couplings_.end(), [&](const auto& coupling) {
            return regular_coupling_ready(coupling);
          });
    }
  }

  bool scalar_column_ready() const {
    return regular_block_column_ready() && scalar_block_shape();
  }

  const char* regular_column_capability() const {
    if constexpr (std::is_same_v<Scalar, ComplexBall>)
      return scalar_block_shape()
          ? "acb-regular-scalar-block-dag-column-v1"
          : "acb-regular-block-dag-column-v2";
    else
      return scalar_block_shape()
          ? "exact-rational-regular-scalar-block-dag-column-v1"
          : "exact-rational-regular-block-dag-column-v2";
  }

  bool regular_singular_jordan_column_ready() const {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      return false;
    } else {
      if (blocks_.empty() || work_.work_min > 0 ||
          work_.requested_max < 0 || work_.work_complete_max < 0 ||
          std::none_of(blocks_.begin(), blocks_.end(), [](const auto& block) {
            return !block.regular;
          }) ||
          std::any_of(blocks_.begin(), blocks_.end(), [](const auto& block) {
            return block.vertices.empty() ||
                   block.chart->dimension() != block.vertices.size() ||
                   !block.chart->has_identity_assembly() ||
                   !block.exact_jordan_indicial.has_value();
          }))
        return false;
      if constexpr (std::is_same_v<Scalar, ComplexBall>)
        if (std::any_of(blocks_.begin(), blocks_.end(),
                        [](const auto& block) {
                          return !block.no_pseudo;
                        }))
          return false;
      return std::all_of(
          couplings_.begin(), couplings_.end(), [&](const auto& coupling) {
            return sector_preserving_coupling_ready(coupling);
          });
    }
  }

  const char* regular_singular_column_capability() const {
    const auto scalar = scalar_block_shape();
    if constexpr (std::is_same_v<Scalar, ComplexBall>)
      return scalar
          ? "acb-regular-singular-scalar-block-dag-column-v1"
          : "acb-regular-singular-jordan-block-dag-column-v1";
    else
      return scalar
          ? "exact-rational-regular-singular-scalar-block-dag-column-v1"
          : "exact-rational-regular-singular-jordan-block-dag-column-v2";
  }

  bool scalar_block_shape() const {
    return std::all_of(
        blocks_.begin(), blocks_.end(), [](const auto& block) {
          return block.vertices.size() == 1;
        });
  }

  bool regular_coupling_ready(
      const CompositeSCCCoupling<Scalar>& coupling) const {
    if (!sector_preserving_coupling_ready(coupling)) return false;
    for (const auto& entry : coupling.matrix.entries) {
      if (entry.multiplier.proven_zero) continue;
      if (std::any_of(entry.multiplier.kernels.begin(),
                      entry.multiplier.kernels.end(),
                      [](const auto& kernel) {
                        return !ScalarTraits<Scalar>::is_zero(kernel.front());
                      }))
        return false;
    }
    return true;
  }

  bool sector_preserving_coupling_ready(
      const CompositeSCCCoupling<Scalar>& coupling) const {
    if (coupling.source_block >= blocks_.size() ||
        coupling.target_block >= blocks_.size() ||
        coupling.matrix.columns !=
            blocks_[coupling.source_block].vertices.size() ||
        coupling.matrix.rows !=
            blocks_[coupling.target_block].vertices.size())
      return false;
    bool active = false;
    for (const auto& entry : coupling.matrix.entries) {
      if (entry.row >= coupling.matrix.rows ||
          entry.column >= coupling.matrix.columns)
        return false;
      if (entry.multiplier.proven_zero) continue;
      active = true;
      if (entry.multiplier.center_pole_order != 0 ||
          entry.multiplier.kernels.empty() ||
          std::any_of(entry.multiplier.kernels.begin(),
                      entry.multiplier.kernels.end(),
                      [](const auto& kernel) {
                        return kernel.empty();
                      }))
        return false;
    }
    return active;
  }

  static bool scalar_identical(const Scalar& left, const Scalar& right) {
    if constexpr (std::is_same_v<Scalar, ComplexBall>)
      return acb_equal(left.raw(), right.raw());
    else
      return left == right;
  }

  static Scalar scalar_from_exact_rational(const Rational& value) {
    if constexpr (std::is_same_v<Scalar, Rational>)
      return value;
    else if constexpr (std::is_same_v<Scalar, ComplexBall>)
      return ComplexBall::from_strings(value.str());
    else
      return SymbolicRational(value.str());
  }

  const json::object& checked_column_run(
      const json::object& entry, std::uint32_t block_index, bool seed,
      const LocalSolution<Scalar>* source,
      bool regular_singular_execution) const {
    if (block_index >= blocks_.size())
      throw std::invalid_argument("native SCC run block is out of range");
    const auto block_dimension = blocks_[block_index].chart->dimension();
    const auto frame_width = blocks_[block_index].chart->frame_width();
    const ExactJordanIndicialCertificate* retained_indicial = nullptr;
    if (regular_singular_execution) {
      const auto& certificate =
          blocks_[block_index].exact_jordan_indicial;
      if (!certificate.has_value())
        throw std::invalid_argument(
            "native regular-singular SCC chart has no retained exact affine Jordan indicial certificate");
      retained_indicial = &*certificate;
    }
    const auto& run = as_object(entry.at("run"), "native SCC recurrence run");
    const auto& metadata_object = as_object(
        entry.at("metadata"), "native SCC local metadata");
    validate_metadata_geometry(metadata_object);
    auto metadata = parse_local_metadata(metadata_object);
    if (!run.at("source").is_null())
      throw std::invalid_argument(
          "native SCC column rejects caller-supplied recurrence source data");
    if (run.at("adaptive_probe").as_bool())
      throw std::invalid_argument(
          "native SCC column requires a fixed retained lower frame");
    if (run.at("return_u").as_bool())
      throw std::invalid_argument(
          "native SCC column requires retained assembly without U JSON");
    const auto nmax = as_u32(run.at("nmax"), "native SCC Taylor order");
    const auto log_max = as_u32(run.at("p"), "native SCC log maximum");
    if (nmax != work_.work_t_order)
      throw std::invalid_argument(
          "native SCC run Taylor order differs from its retained work contract");
    if (as_i32(run.at("a_shift_min"),
               "native SCC a-shift minimum") != 0)
      throw std::invalid_argument(
          "native SCC column requires a zero exact a-shift origin");
    const auto& a_shifts = as_array(
        run.at("a_shifts"), "native SCC exact a-shift schedule");
    if (a_shifts.size() != static_cast<std::size_t>(nmax) + 1)
      throw std::invalid_argument(
          "native SCC regular block run has an incomplete a-shift schedule");
    std::optional<Rational> exact_a_target;
    std::optional<Rational> exact_b_target;
    if (regular_singular_execution) {
      if (metadata.a.domain != ExactDomain::Rational ||
          metadata.b.domain != ExactDomain::Rational)
        throw std::invalid_argument(
            "native regular-singular SCC execution requires exact Rational task tags");
      exact_a_target = Rational(metadata.a.canonical);
      exact_b_target = Rational(metadata.b.canonical);
      const auto encoded_a = parse_scalar<Scalar>(run.at("a_target"));
      const auto encoded_b = parse_scalar<Scalar>(run.at("b_target"));
      if (!scalar_identical(
              encoded_a, scalar_from_exact_rational(*exact_a_target)) ||
          !scalar_identical(
              encoded_b, scalar_from_exact_rational(*exact_b_target)))
        throw std::invalid_argument(
            "native regular-singular SCC numeric targets differ from their exact task tags");
      const auto& raw_tag = as_object(
          metadata_object.at("tag"), "native SCC exact task tag");
      if (const auto* raw_p = raw_tag.if_contains("p")) {
        const auto& p = as_object(*raw_p, "native SCC exact log tag");
        if (required_string(p, "domain") != "integer" ||
            Rational(required_string(p, "canonical")) !=
                Rational(std::to_string(log_max)))
          throw std::invalid_argument(
              "native regular-singular SCC exact log tag differs from its recurrence run");
      } else if constexpr (std::is_same_v<Scalar, ComplexBall>) {
        throw std::invalid_argument(
            "native Acb regular-singular SCC metadata must retain the exact log tag");
      }
      for (std::size_t n = 0; n < a_shifts.size(); ++n) {
        const auto expected =
            *exact_a_target + Rational(std::to_string(n));
        if (!scalar_identical(parse_scalar<Scalar>(a_shifts[n]),
                              scalar_from_exact_rational(expected)))
          throw std::invalid_argument(
              "native SCC a-shift schedule must equal a_target plus the exact Taylor index");
      }
    } else {
      const auto a_target = parse_scalar<Scalar>(run.at("a_target"));
      const auto b_target = parse_scalar<Scalar>(run.at("b_target"));
      if (log_max != 0 || !ScalarTraits<Scalar>::is_zero(a_target) ||
          !ScalarTraits<Scalar>::is_zero(b_target))
        throw std::invalid_argument(
            "native SCC regular block runs require p=0 and exact a=b=0");
      for (std::size_t n = 0; n < a_shifts.size(); ++n) {
        const auto shift = parse_scalar<Scalar>(a_shifts[n]);
        const auto expected = ScalarTraits<Scalar>::integer(
            static_cast<long>(n));
        if (!scalar_identical(shift, expected))
          throw std::invalid_argument(
              "native SCC regular a-shift schedule must equal the exact Taylor index");
      }
    }
    if (seed) {
      if (!run.at("has_initial").as_bool())
        throw std::invalid_argument(
            "native SCC seed requires one initialized log-zero sector");
    } else {
      if (run.at("has_initial").as_bool())
        throw std::invalid_argument(
            "native SCC target rejects caller-supplied initial state");
      const auto& initial = as_array(
          run.at("initial"), "native SCC target initial tensor");
      const auto& validity = as_array(
          run.at("initial_validity"),
          "native SCC target initial validity");
      const auto log_points = static_cast<std::size_t>(log_max) + 1;
      if (log_points > std::numeric_limits<std::size_t>::max() /
                           block_dimension)
        throw std::overflow_error(
            "native SCC target initial validity size overflow");
      const auto expected_validity = log_points * block_dimension;
      if (expected_validity > std::numeric_limits<std::size_t>::max() /
                                  frame_width)
        throw std::overflow_error(
            "native SCC target initial tensor size overflow");
      const auto expected_initial_coefficients = expected_validity * frame_width;
      if (initial.size() != expected_initial_coefficients ||
          validity.size() != expected_validity ||
          std::any_of(initial.begin(), initial.end(), [](const auto& value) {
            return !ScalarTraits<Scalar>::is_zero(
                parse_scalar<Scalar>(value));
          }) ||
          std::any_of(validity.begin(), validity.end(), [](const auto& value) {
            return !value.is_null();
          }))
        throw std::invalid_argument(
            "native SCC target initial template must be explicit exact zero with unknown validity");
      if (source == nullptr)
        throw std::logic_error("native SCC target source was not constructed");
      if (source->dimension != block_dimension)
        throw std::invalid_argument(
            "native SCC coupling source dimension differs from its target block");
      for (const auto& sector : source->sectors)
        if (sector.log_power > log_max)
          throw std::invalid_argument(
              "native SCC target log template does not cover its exact source sectors");
    }
    const auto& schedule = as_array(
        run.at("schedule"), "native SCC recurrence schedule");
    if (schedule.size() != static_cast<std::size_t>(nmax) + 1)
      throw std::invalid_argument(
          "native SCC run has an inconsistent recurrence schedule height");
    if (regular_singular_execution) {
      std::vector<std::vector<BlockStep<Rational>>> exact_schedule;
      exact_schedule.reserve(schedule.size());
      for (std::size_t n = 0; n < schedule.size(); ++n) {
        const auto& raw_row = schedule[n];
        const auto& row = as_array(raw_row, "native SCC schedule row");
        if (row.size() != retained_indicial->blocks.size())
          throw std::invalid_argument(
              "native regular-singular SCC execution requires one step per retained exact Jordan block");
        std::vector<BlockStep<Rational>> exact_row;
        exact_row.reserve(row.size());
        for (std::size_t block = 0; block < row.size(); ++block) {
          const auto& raw_step = row[block];
          const auto& step = as_object(raw_step, "native SCC schedule step");
          const auto kind = required_string(step, "case");
          const auto step_case = kind == "T" ? StepCase::Taylor
              : kind == "P" ? StepCase::Pseudo
              : kind == "R" ? StepCase::Resonant
              : throw std::invalid_argument(
                    "native SCC schedule has an unknown recurrence case");
          if constexpr (std::is_same_v<Scalar, Rational>) {
            exact_row.push_back({step_case,
                parse_scalar<Rational>(step.at("da")),
                parse_scalar<Rational>(step.at("db"))});
          } else {
            const auto& exact_block = retained_indicial->blocks[block];
            const auto expected_da = *exact_a_target +
                Rational(std::to_string(n)) - exact_block.root.a;
            const auto expected_db =
                *exact_b_target - exact_block.root.b;
            if (!scalar_identical(
                    parse_scalar<Scalar>(step.at("da")),
                    scalar_from_exact_rational(expected_da)) ||
                !scalar_identical(
                    parse_scalar<Scalar>(step.at("db")),
                    scalar_from_exact_rational(expected_db)))
              throw std::invalid_argument(
                  "native Acb SCC schedule enclosure differs from its exact affine-Jordan offsets");
            exact_row.push_back(
                {step_case, expected_da, expected_db});
          }
        }
        exact_schedule.push_back(std::move(exact_row));
      }
      const auto schedule_certificate =
          certify_exact_affine_jordan_schedule(
          *retained_indicial, *exact_a_target, *exact_b_target,
          exact_schedule);
      if constexpr (std::is_same_v<Scalar, ComplexBall>)
        if (schedule_certificate.contains_pseudo)
          throw RecurrenceError(
              "E5",
              "native Acb regular-singular SCC execution rejects exact CASE-P collisions; use the exact Rational compensation path");
    } else {
      if (!blocks_[block_index].chart->has_regular_singleton_partition())
        throw std::invalid_argument(
            "native SCC regular block lost its retained singleton partition");
      for (std::size_t n = 0; n < schedule.size(); ++n) {
        const auto& row = as_array(schedule[n], "native SCC schedule row");
        if (row.size() != static_cast<std::size_t>(block_dimension))
          throw std::invalid_argument(
              "native SCC regular execution requires one step per retained Jordan singleton");
        const auto expected_kind = n == 0 ? "R" : "T";
        const auto expected_da = ScalarTraits<Scalar>::integer(
            static_cast<long>(n));
        for (const auto& raw_step : row) {
          const auto& step = as_object(raw_step, "native SCC schedule step");
          const auto kind = required_string(step, "case");
          const auto da = parse_scalar<Scalar>(step.at("da"));
          const auto db = parse_scalar<Scalar>(step.at("db"));
          if (kind != expected_kind ||
              !ScalarTraits<Scalar>::is_zero(db) ||
              !scalar_identical(da, expected_da))
            throw std::invalid_argument(
                "native SCC regular block schedule must be resonant at zero and Taylor by exact index for every Jordan singleton");
        }
      }
    }
    return run;
  }

  std::uint32_t seed_component_from_run(
      const json::object& run, std::uint32_t block_index,
      bool regular_singular_execution) const {
    const auto dimension = blocks_[block_index].chart->dimension();
    const auto frame_width = blocks_[block_index].chart->frame_width();
    const auto& initial = as_array(
        run.at("initial"), "native SCC seed initial tensor");
    const auto& validity = as_array(
        run.at("initial_validity"), "native SCC seed initial validity");
    const auto unit_index_i64 = -static_cast<std::int64_t>(work_.work_min);
    if (regular_singular_execution) {
      const auto& indicial =
          blocks_[block_index].exact_jordan_indicial;
      if (!indicial.has_value())
        throw std::logic_error(
            "native regular-singular seed lost its retained indicial certificate");
      const auto log_max = as_u32(
          run.at("p"), "native regular-singular seed log maximum");
      const auto log_points = static_cast<std::size_t>(log_max) + 1;
      if (unit_index_i64 < 0 || unit_index_i64 >= frame_width ||
          log_points > std::numeric_limits<std::size_t>::max() / dimension ||
          log_points * dimension >
              std::numeric_limits<std::size_t>::max() / frame_width ||
          initial.size() != log_points * dimension * frame_width ||
          validity.size() != log_points * dimension)
        throw std::invalid_argument(
            "native regular-singular SCC seed has a malformed finite Jordan/log tensor");
      for (const auto& value : validity)
        if (value.is_null() ||
            as_i32(value, "native regular-singular seed validity") !=
                work_.work_complete_max)
          throw std::invalid_argument(
              "native regular-singular SCC seed requires finite validity through the retained work maximum");

      const auto unit_index = static_cast<std::size_t>(unit_index_i64);
      const auto initial_index = [&](std::uint32_t log,
                                     std::uint32_t component,
                                     std::size_t epsilon) {
        return ((static_cast<std::size_t>(log) * dimension + component) *
                frame_width) + epsilon;
      };
      const auto one = ScalarTraits<Scalar>::one();
      std::optional<std::uint32_t> selected;
      for (std::uint32_t component = 0; component < dimension; ++component) {
        for (std::size_t epsilon = 0; epsilon < frame_width; ++epsilon) {
          const auto coefficient = parse_scalar<Scalar>(
              initial[initial_index(0, component, epsilon)]);
          if (epsilon == unit_index && scalar_identical(coefficient, one)) {
            if (selected.has_value())
              throw std::invalid_argument(
                  "native regular-singular SCC seed contains more than one log-zero eps^0 unit component");
            selected = component;
          } else if (!ScalarTraits<Scalar>::is_zero(coefficient)) {
            throw std::invalid_argument(
                "native regular-singular SCC seed log-zero frame is not one exact eps^0 unit component");
          }
        }
      }
      if (!selected.has_value())
        throw std::invalid_argument(
            "native regular-singular SCC seed contains no log-zero eps^0 unit component");

      const auto spectral_block_index = indicial->block_of_column[*selected];
      const auto position = indicial->position_in_block[*selected];
      const auto& spectral_block =
          indicial->blocks[spectral_block_index];
      const auto a_target = parse_scalar<Scalar>(run.at("a_target"));
      const auto b_target = parse_scalar<Scalar>(run.at("b_target"));
      if (!scalar_identical(
              a_target, scalar_from_exact_rational(spectral_block.root.a)) ||
          !scalar_identical(
              b_target, scalar_from_exact_rational(spectral_block.root.b)))
        throw std::invalid_argument(
            "native regular-singular SCC seed tag is not the exact affine root of its selected Jordan chain");
      if (log_max < position)
        throw std::invalid_argument(
            "native regular-singular SCC seed log ceiling truncates its exact Jordan chain");

      for (std::uint32_t log = 0; log <= log_max; ++log) {
        std::optional<std::uint32_t> expected_component;
        std::optional<std::size_t> expected_epsilon;
        if (log <= position) {
          expected_component = spectral_block.columns[position - log];
          const auto epsilon_i64 = -static_cast<std::int64_t>(log) -
                                   work_.work_min;
          if (epsilon_i64 < 0 || epsilon_i64 >= frame_width)
            throw RecurrenceError(
                "E4",
                "canonical Jordan seed exceeds the retained lower epsilon frame",
                work_.work_min, -static_cast<std::int32_t>(log));
          expected_epsilon = static_cast<std::size_t>(epsilon_i64);
        }
        for (std::uint32_t component = 0; component < dimension;
             ++component) {
          for (std::size_t epsilon = 0; epsilon < frame_width; ++epsilon) {
            const auto coefficient = parse_scalar<Scalar>(
                initial[initial_index(log, component, epsilon)]);
            const auto expected = expected_component.has_value() &&
                                  component == *expected_component &&
                                  epsilon == *expected_epsilon;
            if ((expected && !scalar_identical(coefficient, one)) ||
                (!expected && !ScalarTraits<Scalar>::is_zero(coefficient)))
              throw std::invalid_argument(
                  "native regular-singular SCC seed differs from the captured canonical Jordan/log normalization");
          }
        }
      }
      return *selected;
    }
    if (unit_index_i64 < 0 || unit_index_i64 >= frame_width ||
        dimension > std::numeric_limits<std::size_t>::max() / frame_width ||
        initial.size() != static_cast<std::size_t>(dimension) * frame_width ||
        validity.size() != dimension)
      throw std::invalid_argument(
          "native SCC regular block seed requires one honest finite eps^0 unit frame");
    for (const auto& value : validity)
      if (value.is_null() ||
          as_i32(value, "native SCC seed validity") !=
              work_.work_complete_max)
        throw std::invalid_argument(
            "native SCC regular block seed requires finite validity through the retained work maximum");

    const auto unit_index = static_cast<std::size_t>(unit_index_i64);
    const auto one = ScalarTraits<Scalar>::one();
    std::optional<std::uint32_t> selected;
    for (std::uint32_t component = 0; component < dimension; ++component) {
      for (std::size_t epsilon = 0; epsilon < frame_width; ++epsilon) {
        const auto coefficient = parse_scalar<Scalar>(
            initial[static_cast<std::size_t>(component) * frame_width +
                    epsilon]);
        if (epsilon == unit_index && scalar_identical(coefficient, one)) {
          if (selected.has_value())
            throw std::invalid_argument(
                "native SCC seed contains more than one eps^0 unit component");
          selected = component;
        } else if (!ScalarTraits<Scalar>::is_zero(coefficient)) {
          throw std::invalid_argument(
              "native SCC regular block seed must be exactly one eps^0 unit column");
        }
      }
    }
    if (!selected.has_value())
      throw std::invalid_argument(
          "native SCC regular block seed contains no eps^0 unit component");
    return *selected;
  }

  void validate_metadata_geometry(const json::object& metadata) const {
    const auto& chart = as_object(
        metadata.at("chart"), "native SCC local chart metadata");
    const auto infinite = chart.at("infinite_radius").as_bool();
    if (required_string(chart, "center_exact") !=
            retained_geometry_.chart.center_exact ||
        required_string(chart, "scale_exact") !=
            retained_geometry_.chart.scale_exact ||
        infinite != retained_geometry_.chart.infinite_radius)
      throw std::invalid_argument(
          "native SCC local metadata differs from retained parent chart geometry");
    if (!infinite) {
      const auto* radius = chart.if_contains("radius");
      if (radius == nullptr || !radius->is_string() ||
          Rational(std::string(radius->as_string())).str() !=
              retained_geometry_.radius_exact)
        throw std::invalid_argument(
            "native SCC local radius differs from retained exact parent radius");
    }

    const auto& prescriptions = as_array(
        metadata.at("prescriptions"),
        "native SCC local prescriptions");
    if (prescriptions.size() != retained_geometry_.prescriptions.size())
      throw std::invalid_argument(
          "native SCC local prescriptions differ from retained parent prescriptions");
    for (std::size_t index = 0; index < prescriptions.size(); ++index) {
      const auto& raw = as_object(
          prescriptions[index], "native SCC local prescription");
      const auto& retained = retained_geometry_.prescriptions[index];
      if (required_string(raw, "factor_exact") != retained.factor_exact ||
          as_i32(raw.at("sign"), "prescription sign") != retained.sign ||
          as_u32(raw.at("multiplicity"), "prescription multiplicity") !=
              retained.multiplicity ||
          as_i32(raw.at("leading_coefficient_sign"),
                 "prescription leading coefficient sign") !=
              retained.leading_coefficient_sign)
        throw std::invalid_argument(
            "native SCC local prescriptions differ from retained parent prescriptions");
    }
  }

  void require_work_local(const LocalSolution<Scalar>& solution,
                          const char* label) const {
    validate_local_solution(solution, false);
    if (!solution.error.empty() ||
        !local_algebra_detail::same_chart(
            solution.chart, retained_geometry_.chart) ||
        !local_algebra_detail::same_prescriptions(
            solution.prescriptions, retained_geometry_.prescriptions))
      throw std::invalid_argument(
          std::string(label) +
          " differs from retained composite geometry or carries an unsupported error envelope");
    if (solution.epsilon.min_power < work_.work_min ||
        solution.epsilon.complete_max > work_.work_complete_max ||
        solution.taylor_complete_max != work_.work_t_order)
      throw std::invalid_argument(
          std::string(label) + " lies outside the retained SCC work rectangle");
  }

  void validate_block_result(const NativeLocalRun<Scalar>& native,
                             std::uint32_t block_index, bool seed,
                             bool regular_singular_execution) const {
    if (!native.pseudo_hits.empty() && regular_singular_execution &&
        std::is_same_v<Scalar, ComplexBall>)
      throw RecurrenceError(
          "E5",
          "native Acb regular-singular SCC recurrence produced an unsupported pseudo hit after exact no-CASE-P certification");
    if (!native.pseudo_hits.empty())
      throw std::invalid_argument(
          "native SCC execution encountered pseudo hits after its exact schedule certificate");
    if (block_index >= blocks_.size() ||
        native.solution.dimension != blocks_[block_index].vertices.size())
      throw std::invalid_argument(
          "native SCC block solve returned the wrong retained dimension");
    require_work_local(native.solution, "native SCC block result");
    if (native.solution.sectors.empty())
      throw std::invalid_argument("native SCC block solve returned no sectors");
    if (seed && !regular_singular_execution &&
        (native.solution.sectors.size() != 1 ||
         native.solution.sectors.front().log_power != 0))
      throw std::invalid_argument(
          "native SCC seed must assemble exactly one log-zero sector");
    for (const auto& sector : native.solution.sectors)
      if (sector.a.domain != ExactDomain::Rational ||
          sector.b.domain != ExactDomain::Rational)
        throw std::invalid_argument(
            "native SCC column execution requires exact rational sector tags");
  }

  void validate_column_coupling(
      const CompositeSCCCoupling<Scalar>& coupling,
      bool regular_singular_execution) const {
    const bool ready = regular_singular_execution
        ? sector_preserving_coupling_ready(coupling)
        : regular_coupling_ready(coupling);
    if (!ready)
      throw std::invalid_argument(
          regular_singular_execution
              ? "native regular-singular SCC column requires a dimension-matched exact coupling matrix which preserves the sector tag"
              : "native SCC column requires a dimension-matched exact pole-free coupling matrix whose active entries vanish at chart center");
  }

  void require_source_tag_matches_run(
      const LocalSolution<Scalar>& source, const json::object& run) const {
    const auto a_target = parse_scalar<Scalar>(run.at("a_target"));
    const auto b_target = parse_scalar<Scalar>(run.at("b_target"));
    const auto log_max = as_u32(run.at("p"), "target log maximum");
    for (const auto& sector : source.sectors) {
      if (sector.a.domain != ExactDomain::Rational ||
          sector.b.domain != ExactDomain::Rational ||
          sector.log_power > log_max)
        throw std::invalid_argument(
            "native SCC coupling source tag/log sector differs from its target run");
      verify_tag_binding<Scalar>(sector.a, a_target,
                                 "native SCC coupling source a");
      verify_tag_binding<Scalar>(sector.b, b_target,
                                 "native SCC coupling source b");
    }
  }

  static void accumulate_diagnostics(
      NativeLocalDiagnostics& total,
      const NativeLocalDiagnostics& current) {
    total.top_valid = std::min(total.top_valid, current.top_valid);
    total.parse_ms += current.parse_ms;
    total.kernel_ms += current.kernel_ms;
  }

  static json::object block_diagnostic(
      std::uint32_t block, const char* role,
      const std::vector<std::uint32_t>& predecessors,
      const LocalSolution<Scalar>* source,
      const NativeLocalRun<Scalar>& result) {
    json::object diagnostic{
        {"block", block}, {"role", role},
        {"predecessors", encode_indices(predecessors)},
        {"result_epsilon_min", result.solution.epsilon.min_power},
        {"result_epsilon_max", result.solution.epsilon.complete_max},
        {"result_taylor_max", result.solution.taylor_complete_max},
        {"result_sectors", result.solution.sectors.size()},
        {"pseudo_hit_count", result.diagnostics.pseudo_hits},
        {"pseudo_compensation_count",
         result.diagnostics.pseudo_compensations},
        {"max_pseudo_depth", result.diagnostics.max_pseudo_depth},
        {"pseudo_value_certified",
         result.diagnostics.pseudo_value_certified},
        {"uncompensated_pseudo_hit_count", result.pseudo_hits.size()},
        {"top_valid", encode_validity(result.diagnostics.top_valid)},
        {"parse_ms", result.diagnostics.parse_ms},
        {"kernel_ms", result.diagnostics.kernel_ms}};
    if (source != nullptr) {
      diagnostic["source_epsilon_min"] = source->epsilon.min_power;
      diagnostic["source_epsilon_max"] = source->epsilon.complete_max;
      diagnostic["source_taylor_max"] = source->taylor_complete_max;
      diagnostic["source_sectors"] = source->sectors.size();
    }
    return diagnostic;
  }

  std::uint32_t dimension_ = 0;
  SCCCertificate graph_;
  std::string exact_system_record_;
  std::string exact_theta_record_;
  std::string geometry_record_;
  RetainedCompositeGeometry retained_geometry_;
  CompositeWorkContract work_;
  std::vector<CompositeSCCBlock<Scalar>> blocks_;
  std::vector<CompositeSCCCoupling<Scalar>> couplings_;
  std::atomic<std::uint64_t> column_solves_{0};
  mutable std::mutex column_stats_mutex_;
  double column_solve_ms_ = 0.0;
};

constexpr const char* kRetainedTilePlanCapability =
    "retained-exact-independent-arm-tile-plan-v1";
constexpr const char* kRetainedSingleArmTilePlanCapability =
    "retained-exact-single-arm-tile-plan-v1";
constexpr const char* kRetainedSingleArmTilePlanCheckpointSchema =
    "diffexp2-retained-single-arm-tile-plan-v1";
constexpr const char* kRetainedSingleArmTilePlanProvenanceSchema =
    "diffexp2-retained-exact-single-arm-tile-plan-v1";
constexpr const char* kRetainedPlannedMatchHopCapability =
    "retained-exact-plan-driven-local-match-hop-v1";
constexpr const char* kRetainedPlannedMatchMaterializationCapability =
    "retained-native-plan-match-local-materialization-v1";
constexpr const char* kRetainedStoredLineCapability =
    "retained-native-stored-truncation-physical-tile-integral-v1";
constexpr const char* kRetainedCertifiedLineCapability =
    "retained-native-certified-full-local-physical-tile-integral-v1";
constexpr const char* kRetainedParallelArmCapability =
    "retained-native-concurrent-two-arm-march-v1";
constexpr const char* kRetainedTransportArmStateCapability =
    "retained-native-transport-arm-state-v1";
constexpr const char* kRetainedLineAggregateCapability =
    "retained-native-line-aggregate-v1";

const char* line_integration_scope_name(LineIntegrationScope scope) {
  switch (scope) {
    case LineIntegrationScope::StoredTruncation:
      return "stored_truncation";
    case LineIntegrationScope::FullLocalWithCertifiedTail:
      return "full_local_with_certified_tail";
  }
  throw std::logic_error("unknown line integration scope");
}

const char* error_guarantee_name(ErrorGuarantee guarantee) {
  switch (guarantee) {
    case ErrorGuarantee::None:
      return "none";
    case ErrorGuarantee::Advisory:
      return "advisory";
    case ErrorGuarantee::Certified:
      return "certified";
  }
  throw std::logic_error("unknown error guarantee");
}

json::object encode_error_envelope_summary(const ErrorEnvelope& error) {
  json::object result{
      {"guarantee", error_guarantee_name(error.guarantee)},
      {"provenance", error.provenance}};
  if (!error.empty()) {
    json::array upper;
    upper.reserve(error.absolute.size());
    for (const auto& bound : error.absolute)
      upper.push_back(bound.approximate_upper());
    result["epsilon_min"] = error.frame.min_power;
    result["epsilon_max"] = error.frame.complete_max;
    result["absolute_upper_approx"] = std::move(upper);
    result["bound_encoding"] = "approximate-double-diagnostics";
  }
  return result;
}

struct RetainedPlanChartBinding {
  using Owner = std::variant<std::shared_ptr<PreparedChartBase>,
                             std::shared_ptr<CompositeSCCChartBase>>;

  std::string handle;
  std::string exact_identity;
  ExactAffineChart geometry;
  std::vector<Prescription> prescriptions;
  Owner owner;
};

struct RetainedArmPlan {
  ExactArmPlan exact;
  std::vector<RetainedPlanChartBinding> charts;
};

json::array encode_path_branch_sheets(
    const std::vector<ExactBranchSheet>& sheets) {
  json::array encoded;
  encoded.reserve(sheets.size());
  for (const auto& sheet : sheets)
    encoded.push_back(json::object{{"factor_exact", sheet.factor_exact},
                                   {"sign", sheet.imaginary_sign}});
  return encoded;
}

json::array encode_plan_prescriptions(
    const std::vector<Prescription>& prescriptions) {
  json::array encoded;
  encoded.reserve(prescriptions.size());
  for (const auto& prescription : prescriptions)
    encoded.push_back(json::object{
        {"factor_exact", prescription.factor_exact},
        {"sign", prescription.sign},
        {"multiplicity", prescription.multiplicity},
        {"leading_coefficient_sign",
         prescription.leading_coefficient_sign}});
  return encoded;
}

json::object encode_plan_chart(const RetainedPlanChartBinding& binding,
                               std::size_t index) {
  return json::object{
      {"index", index}, {"chart", binding.handle},
      {"identity", binding.exact_identity},
      {"center_exact", binding.geometry.center.str()},
      {"scale_exact", binding.geometry.scale.str()},
      {"radius_exact", binding.geometry.radius.str()},
      {"singular_center", binding.geometry.singular_center},
      {"prescriptions", encode_plan_prescriptions(binding.prescriptions)}};
}

const char* exact_match_kind_name(ExactMatchKind kind) {
  switch (kind) {
    case ExactMatchKind::SymmetricDivisionPoint:
      return "symmetric-division-point";
    case ExactMatchKind::BalancedSafeOverlap:
      return "balanced-safe-overlap";
    case ExactMatchKind::ForbiddenPointAvoidance:
      return "forbidden-point-avoidance";
  }
  throw std::logic_error("unknown exact match kind");
}

json::object encode_plan_match(const RetainedArmPlan& arm,
                               std::size_t index) {
  if (index >= arm.exact.matches.size())
    throw std::invalid_argument("native tile-plan match index is out of range");
  const auto& match = arm.exact.matches[index];
  return json::object{
      {"index", index},
      {"producing_chart_index", match.producing_chart},
      {"receiving_chart_index", match.receiving_chart},
      {"producing_chart", arm.charts.at(match.producing_chart).handle},
      {"receiving_chart", arm.charts.at(match.receiving_chart).handle},
      {"physical_exact", match.physical.str()},
      {"producing_local_exact", match.producing_local.str()},
      {"receiving_local_exact", match.receiving_local.str()},
      {"kind", exact_match_kind_name(match.kind)},
      {"branch_sheets", encode_path_branch_sheets(match.branch_sheets)}};
}

json::object encode_plan_tile(const RetainedArmPlan& arm,
                              std::size_t index) {
  if (index >= arm.exact.tiles.size())
    throw std::invalid_argument("native tile-plan tile index is out of range");
  const auto& tile = arm.exact.tiles[index];
  const auto& binding = arm.charts.at(tile.chart);
  return json::object{
      {"index", index}, {"chart_index", tile.chart},
      {"chart", binding.handle}, {"chart_identity", binding.exact_identity},
      {"physical_begin_exact", tile.physical_begin.str()},
      {"physical_end_exact", tile.physical_end.str()},
      {"local_begin_exact", tile.local_begin.str()},
      {"local_end_exact", tile.local_end.str()},
      {"jacobian_exact", binding.geometry.scale.str()},
      {"crosses_singular_center", tile.crosses_singular_center},
      {"branch_sheets", encode_path_branch_sheets(tile.branch_sheets)},
      {"prescriptions", encode_plan_prescriptions(binding.prescriptions)}};
}

json::object encode_plan_topology(const ExactPathTopology& topology) {
  json::array singular_points;
  for (const auto& point : topology.singular_points)
    singular_points.push_back(json::value(point.str()));
  json::array boundary_points;
  for (const auto& point : topology.boundary_points)
    boundary_points.push_back(json::value(point.str()));
  json::array projections;
  for (const auto& projection : topology.complex_projections)
    projections.push_back(json::object{
        {"source_identity", projection.source_identity},
        {"real_part_exact", projection.real_part.str()},
        {"imaginary_magnitude_exact", projection.imaginary_magnitude.str()},
        {"retain_minus_imaginary", projection.retain_minus_imaginary},
        {"retain_real_part", projection.retain_real_part},
        {"retain_plus_imaginary", projection.retain_plus_imaginary}});
  return json::object{
      {"singular_points", std::move(singular_points)},
      {"boundary_points", std::move(boundary_points)},
      {"complex_projections", std::move(projections)},
      {"branch_sheets", encode_path_branch_sheets(topology.branch_sheets)}};
}

json::object encode_retained_arm(const RetainedArmPlan& arm) {
  json::array charts;
  for (std::size_t index = 0; index < arm.charts.size(); ++index)
    charts.push_back(encode_plan_chart(arm.charts[index], index));
  json::array matches;
  for (std::size_t index = 0; index < arm.exact.matches.size(); ++index)
    matches.push_back(encode_plan_match(arm, index));
  json::array tiles;
  for (std::size_t index = 0; index < arm.exact.tiles.size(); ++index)
    tiles.push_back(encode_plan_tile(arm, index));
  return json::object{
      {"from_exact", arm.exact.from.str()},
      {"to_exact", arm.exact.to.str()},
      {"direction", arm.exact.direction},
      {"division_order", arm.exact.division_order},
      {"charts", std::move(charts)}, {"matches", std::move(matches)},
      {"tiles", std::move(tiles)},
      {"topology", encode_plan_topology(arm.exact.topology)}};
}

std::optional<std::int32_t> exact_plan_rim(
    const std::vector<Prescription>& prescriptions,
    const Rational& chart_scale) {
  if (chart_scale.is_zero())
    throw std::invalid_argument(
        "prepared tile chart has a zero exact scale");
  const auto scale_sign = chart_scale.sign();
  std::optional<std::int32_t> rim;
  for (const auto& prescription : prescriptions) {
    if ((prescription.multiplicity & 1U) == 0)
      throw std::invalid_argument(
          "prepared tile chart has an even-multiplicity tangential "
          "prescription; a one-sided real-axis rim is not defined");
    const auto candidate =
        prescription.sign * prescription.leading_coefficient_sign *
        scale_sign;
    if (rim.has_value() && *rim != candidate)
      throw std::invalid_argument(
          "prepared tile chart has conflicting exact odd-multiplicity prescriptions");
    rim = candidate;
  }
  return rim;
}

class StoredTilePlan {
 public:
  StoredTilePlan(std::string handle, std::string checkpoint_identity,
                 std::string provenance_identity, std::uint32_t division_order,
                 RetainedArmPlan lower, RetainedArmPlan upper,
                 double elapsed_ms)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        division_order_(division_order), lower_(std::move(lower)),
        upper_(std::move(upper)), elapsed_ms_(elapsed_ms) {
    validate_arm_set();
  }

  StoredTilePlan(std::string handle, std::string checkpoint_identity,
                 std::string provenance_identity, std::uint32_t division_order,
                 std::string arm_name, RetainedArmPlan arm,
                 double elapsed_ms)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        division_order_(division_order), elapsed_ms_(elapsed_ms) {
    if (arm_name == "lower")
      lower_ = std::move(arm);
    else if (arm_name == "upper")
      upper_ = std::move(arm);
    else
      throw std::invalid_argument(
          "single-arm tile plan name must be lower or upper");
    validate_arm_set();
  }

  const std::string& handle() const { return handle_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& provenance_identity() const {
    return provenance_identity_;
  }

  bool has_arm(const std::string& name) const {
    if (name == "lower") return lower_.has_value();
    if (name == "upper") return upper_.has_value();
    return false;
  }

  bool has_two_arms() const {
    return lower_.has_value() && upper_.has_value();
  }

  const RetainedArmPlan& arm(const std::string& name) const {
    if (name == "lower" && lower_.has_value()) return *lower_;
    if (name == "upper" && upper_.has_value()) return *upper_;
    if (name != "lower" && name != "upper")
      throw std::invalid_argument(
          "native tile-plan arm must be lower or upper");
    throw std::invalid_argument(
        "native tile plan does not retain the requested " + name + " arm");
  }

  json::object match_interval(const std::string& name,
                              std::size_t index) const {
    match_queries_.fetch_add(1);
    auto encoded = encode_plan_match(arm(name), index);
    encoded["arm"] = name;
    encoded["tile_plan"] = handle_;
    encoded["checkpoint_identity"] = checkpoint_identity_;
    return encoded;
  }

  json::object tile_interval(const std::string& name,
                             std::size_t index) const {
    tile_queries_.fetch_add(1);
    auto encoded = encode_plan_tile(arm(name), index);
    encoded["arm"] = name;
    encoded["tile_plan"] = handle_;
    encoded["checkpoint_identity"] = checkpoint_identity_;
    return encoded;
  }

  void note_integration() { integrations_.fetch_add(1); }

  void note_match_advance(const std::string& name) {
    (void)arm(name);
    if (name == "lower") {
      lower_match_advances_.fetch_add(1);
      return;
    }
    if (name == "upper") {
      upper_match_advances_.fetch_add(1);
      return;
    }
    throw std::invalid_argument("native tile-plan arm must be lower or upper");
  }

  std::vector<std::shared_ptr<PreparedChartBase>> dependency_charts() const {
    std::vector<std::shared_ptr<PreparedChartBase>> result;
    result.reserve((lower_.has_value() ? lower_->charts.size() : 0) +
                   (upper_.has_value() ? upper_->charts.size() : 0));
    const auto append = [&](const RetainedPlanChartBinding& binding) {
      std::visit(
          [&](const auto& owner) {
            using Owner = typename std::decay_t<decltype(owner)>::element_type;
            if constexpr (std::is_same_v<Owner, PreparedChartBase>) {
              result.push_back(owner);
            } else {
              auto dependencies = owner->dependency_charts();
              result.insert(result.end(), dependencies.begin(),
                            dependencies.end());
            }
          },
          binding.owner);
    };
    if (lower_.has_value())
      for (const auto& binding : lower_->charts) append(binding);
    if (upper_.has_value())
      for (const auto& binding : upper_->charts) append(binding);
    return result;
  }

  std::vector<std::shared_ptr<CompositeSCCChartBase>> dependency_sccs() const {
    std::vector<std::shared_ptr<CompositeSCCChartBase>> result;
    const auto append = [&](const RetainedPlanChartBinding& binding) {
      if (const auto* owner = std::get_if<
              std::shared_ptr<CompositeSCCChartBase>>(&binding.owner))
        result.push_back(*owner);
    };
    if (lower_.has_value())
      for (const auto& binding : lower_->charts) append(binding);
    if (upper_.has_value())
      for (const auto& binding : upper_->charts) append(binding);
    return result;
  }

  json::object summary(bool include_intervals = true) const {
    if (!has_two_arms()) return single_arm_summary(include_intervals);
    json::object result{
        {"tile_plan", handle_}, {"capability", kRetainedTilePlanCapability},
        {"native_retained", true}, {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"division_order", division_order_},
        {"independent_arms", true},
        {"concurrent_execution", "immutable-independent-arm-snapshots"},
        {"anchor", json::object{
             {"lower_chart", lower_->charts.front().handle},
             {"upper_chart", upper_->charts.front().handle},
             {"center_exact", lower_->exact.from.str()}}},
        {"lower_matches", lower_->exact.matches.size()},
        {"upper_matches", upper_->exact.matches.size()},
        {"lower_tiles", lower_->exact.tiles.size()},
        {"upper_tiles", upper_->exact.tiles.size()},
        {"match_interval_queries", match_queries_.load()},
        {"tile_interval_queries", tile_queries_.load()},
        {"lower_match_advances", lower_match_advances_.load()},
        {"upper_match_advances", upper_match_advances_.load()},
        {"integrations", integrations_.load()},
        {"elapsed_ms", elapsed_ms_}};
    if (include_intervals) {
      result["lower"] = encode_retained_arm(*lower_);
      result["upper"] = encode_retained_arm(*upper_);
    }
    return result;
  }

  json::object checkpoint_record() const {
    if (!has_two_arms()) {
      const auto name = single_arm_name();
      return json::object{
          {"schema", kRetainedSingleArmTilePlanCheckpointSchema},
          {"handle", handle_},
          {"checkpoint_identity", checkpoint_identity_},
          {"provenance_identity", provenance_identity_},
          {"division_order", division_order_},
          {"arm_name", name},
          {"arm", encode_retained_arm(arm(name))},
          {"elapsed_ms", elapsed_ms_},
          {"runtime_stats", runtime_stats_record()}};
    }
    return json::object{
        {"schema", "diffexp2-retained-tile-plan-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"division_order", division_order_},
        {"lower", encode_retained_arm(*lower_)},
        {"upper", encode_retained_arm(*upper_)},
        {"elapsed_ms", elapsed_ms_},
        {"runtime_stats", runtime_stats_record()}};
  }

  void restore_runtime_stats(std::uint64_t match_queries,
                             std::uint64_t tile_queries,
                             std::uint64_t lower_match_advances,
                             std::uint64_t upper_match_advances,
                             std::uint64_t integrations) {
    if ((!lower_.has_value() && lower_match_advances != 0) ||
        (!upper_.has_value() && upper_match_advances != 0))
      throw std::invalid_argument(
          "single-arm tile-plan checkpoint advances an absent arm");
    match_queries_.store(match_queries);
    tile_queries_.store(tile_queries);
    lower_match_advances_.store(lower_match_advances);
    upper_match_advances_.store(upper_match_advances);
    integrations_.store(integrations);
  }

 private:
  void validate_arm_set() const {
    if (!lower_.has_value() && !upper_.has_value())
      throw std::invalid_argument(
          "retained tile plan must own one or two arms");
    if (lower_.has_value()) validate_exact_arm_plan(lower_->exact);
    if (upper_.has_value()) validate_exact_arm_plan(upper_->exact);
    // Existing two-arm requests name the slots but historically only require
    // opposite directions; preserve that behavior exactly. A genuine
    // single-arm plan derives its retained name from the exact direction.
    if (!has_two_arms()) {
      const auto& retained = lower_.has_value() ? *lower_ : *upper_;
      const auto expected_direction = lower_.has_value() ? -1 : 1;
      if (retained.exact.direction != expected_direction)
        throw std::invalid_argument(
            "retained single tile-arm name differs from its exact direction");
    }
  }

  std::string single_arm_name() const {
    if (has_two_arms())
      throw std::logic_error(
          "two-arm tile plan has no single arm name");
    return lower_.has_value() ? "lower" : "upper";
  }

  json::object runtime_stats_record() const {
    return json::object{
        {"match_interval_queries", match_queries_.load()},
        {"tile_interval_queries", tile_queries_.load()},
        {"lower_match_advances", lower_match_advances_.load()},
        {"upper_match_advances", upper_match_advances_.load()},
        {"integrations", integrations_.load()}};
  }

  json::object single_arm_summary(bool include_intervals) const {
    const auto name = single_arm_name();
    const auto& retained = arm(name);
    json::object result{
        {"tile_plan", handle_},
        {"capability", kRetainedSingleArmTilePlanCapability},
        {"native_retained", true},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"division_order", division_order_},
        {"independent_arms", false},
        {"concurrent_execution", "single-immutable-arm-snapshot"},
        {"arm_name", name},
        {"anchor", json::object{
             {"chart", retained.charts.front().handle},
             {"center_exact", retained.exact.from.str()}}},
        {"matches", retained.exact.matches.size()},
        {"tiles", retained.exact.tiles.size()},
        {"match_interval_queries", match_queries_.load()},
        {"tile_interval_queries", tile_queries_.load()},
        {"lower_match_advances", lower_match_advances_.load()},
        {"upper_match_advances", upper_match_advances_.load()},
        {"integrations", integrations_.load()},
        {"elapsed_ms", elapsed_ms_}};
    if (include_intervals) {
      result["arm"] = encode_retained_arm(retained);
      result[name] = encode_retained_arm(retained);
    }
    return result;
  }

  std::string handle_;
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::uint32_t division_order_ = 3;
  std::optional<RetainedArmPlan> lower_;
  std::optional<RetainedArmPlan> upper_;
  double elapsed_ms_ = 0.0;
  mutable std::atomic<std::uint64_t> match_queries_{0};
  mutable std::atomic<std::uint64_t> tile_queries_{0};
  std::atomic<std::uint64_t> lower_match_advances_{0};
  std::atomic<std::uint64_t> upper_match_advances_{0};
  std::atomic<std::uint64_t> integrations_{0};
};

// A plan-driven match is one exact handoff, not a completed arm.  It owns the
// immutable plan snapshot and every local used to construct the retained
// matching weights.  Registry release of the public plan/local tokens can
// therefore never turn a published handoff into dangling provenance.
class StoredPlannedMatchHop final : public StoredMatchBase {
 public:
  StoredPlannedMatchHop(
      std::shared_ptr<StoredMatchBase> match,
      std::string checkpoint_identity, std::string provenance_identity,
      json::object handoff, double elapsed_ms,
      std::shared_ptr<StoredTilePlan> plan_owner,
      std::vector<std::shared_ptr<StoredLocalBase>> basis_owners,
      std::shared_ptr<StoredLocalBase> incoming_owner)
      : StoredMatchBase(match == nullptr ? std::string() : match->handle()),
        match_(std::move(match)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        handoff_(std::move(handoff)), elapsed_ms_(elapsed_ms),
        plan_owner_(std::move(plan_owner)),
        basis_owners_(std::move(basis_owners)),
        incoming_owner_(std::move(incoming_owner)) {
    if (match_ == nullptr || plan_owner_ == nullptr ||
        incoming_owner_ == nullptr || basis_owners_.empty())
      throw std::invalid_argument(
          "retained planned match hop requires all strong owners");
  }

  json::object summary() const override {
    auto result = match_->summary();
    if (required_string(result, "checkpoint_identity") !=
        checkpoint_identity_)
      throw std::logic_error(
          "retained planned match checkpoint identity changed");
    result["planned_hop_capability"] =
        kRetainedPlannedMatchHopCapability;
    result["plan_driven"] = true;
    result["planned_hop_provenance_identity"] = provenance_identity_;
    result["planned_hop"] = handoff_;
    result["strong_ownership"] = json::object{
        {"tile_plan", true}, {"basis_locals", basis_owners_.size()},
        {"incoming_local", true}};
    result["materializations"] = materializations_.load();
    result["elapsed_ms"] = elapsed_ms_;
    return result;
  }

  double elapsed_ms() const { return elapsed_ms_; }

  const std::shared_ptr<StoredMatchBase>& native_match() const {
    return match_;
  }
  const std::shared_ptr<StoredTilePlan>& plan_owner() const {
    return plan_owner_;
  }
  const std::vector<std::shared_ptr<StoredLocalBase>>& basis_owners() const {
    return basis_owners_;
  }
  const std::shared_ptr<StoredLocalBase>& incoming_owner() const {
    return incoming_owner_;
  }
  const json::object& handoff() const { return handoff_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& provenance_identity() const {
    return provenance_identity_;
  }

  void validate_materialized_derivation(
      const json::object& derivation, const char* scalar_domain) const {
    const auto native_summary = match_->summary();
    if (required_string(derivation, "source_match") != handle_ ||
        required_string(derivation, "source_match_checkpoint_identity") !=
            checkpoint_identity_ ||
        required_string(derivation, "source_match_provenance_identity") !=
            required_string(native_summary, "provenance_identity") ||
        required_string(derivation, "planned_hop_provenance_identity") !=
            provenance_identity_ ||
        derivation.at("planned_hop") != handoff_)
      throw std::invalid_argument(
          "checkpoint materialized-local lineage disagrees with its planned-hop owner");

    json::array expected_windows;
    std::int32_t expected_certified_max = 0;
    if (const auto exact =
            std::dynamic_pointer_cast<StoredExactRegularMatch>(match_)) {
      if (std::string(scalar_domain) != "rational")
        throw std::invalid_argument(
            "checkpoint materialized-local scalar domain differs from its exact match owner");
      for (const auto& weight : exact->weights())
        expected_windows.push_back(json::object{
            {"min", weight.min_power()}, {"max", weight.complete_max()}});
      expected_certified_max = as_i32(
          as_object(native_summary.at("residual"),
                    "exact retained match residual").at("max"),
          "exact retained match residual maximum");
    } else if (const auto acb =
                   std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_)) {
      if (std::string(scalar_domain) != "acb" ||
          !acb->certified_for_materialization())
        throw std::invalid_argument(
            "checkpoint materialized-local Acb owner lost its passing complete match certificate");
      for (const auto& weight : acb->weights())
        expected_windows.push_back(json::object{
            {"min", weight.min_power()}, {"max", weight.complete_max()}});
      expected_certified_max = as_i32(
          as_object(native_summary.at("epsilon"),
                    "Acb retained match epsilon").at(
                        "required_complete_max"),
          "Acb retained match required maximum");
    } else {
      throw std::invalid_argument(
          "checkpoint materialized-local owner embeds an unsupported native match");
    }
    if (derivation.at("weight_windows") != expected_windows ||
        as_i32(derivation.at("match_certified_complete_max"),
               "materialized-local certified maximum") !=
            expected_certified_max)
      throw std::invalid_argument(
          "checkpoint materialized-local derivation differs from its retained match weights/certificate");
  }

  std::shared_ptr<StoredLocalBase> materialize(
      const std::string& local_handle,
      const std::string& result_checkpoint_identity,
      slong precision_bits,
      const std::shared_ptr<StoredPlannedMatchHop>& self) {
    if (self.get() != this)
      throw std::logic_error(
          "retained plan-match materialization lost self ownership");
    std::shared_ptr<StoredLocalBase> result;
    if (const auto exact =
            std::dynamic_pointer_cast<StoredExactRegularMatch>(match_)) {
      result = materialize_typed<Rational>(
          local_handle, result_checkpoint_identity, precision_bits,
          exact->weights(), self);
    } else if (const auto acb =
                   std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_)) {
      if (!acb->certified_for_materialization())
        throw std::domain_error(
            "an Acb plan-match handoff must have a passing complete residual before materialization");
      result = materialize_typed<ComplexBall>(
          local_handle, result_checkpoint_identity, precision_bits,
          acb->weights(), self);
    } else {
      throw std::logic_error(
          "retained plan-match handoff has an unsupported matching state");
    }
    materializations_.fetch_add(1);
    return result;
  }

  json::object checkpoint_record() const override {
    return json::object{
        {"schema", "diffexp2-retained-planned-match-hop-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"handoff", handoff_},
        {"native_match", match_->checkpoint_record()},
        {"elapsed_ms", elapsed_ms_},
        {"runtime_stats",
         json::object{{"materializations", materializations_.load()}}}};
  }

  void restore_runtime_stats(std::uint64_t materializations) {
    materializations_.store(materializations);
  }

 private:
  template <typename Scalar>
  std::shared_ptr<StoredLocalBase> materialize_typed(
      const std::string& local_handle,
      const std::string& result_checkpoint_identity,
      slong precision_bits,
      const FiniteLaurentVector<Scalar>& weights,
      const std::shared_ptr<StoredPlannedMatchHop>& self) const {
    if (local_handle.empty() || result_checkpoint_identity.empty())
      throw std::invalid_argument(
          "plan-match local materialization identities must be nonempty");
    AcbPrecisionLease lease(precision_bits);
    ComplexBall::set_precision(precision_bits);
    const auto started = std::chrono::steady_clock::now();
    const auto native_match_summary = match_->summary();
    const auto& match_epsilon = as_object(
        native_match_summary.at("epsilon"),
        "retained match epsilon provenance");
    const auto match_work_max = as_i32(
        match_epsilon.at("max"), "retained match epsilon maximum");
    const auto match_certified_max = [&]() {
      if constexpr (std::is_same_v<Scalar, Rational>)
        return as_i32(
            as_object(native_match_summary.at("residual"),
                      "exact retained match residual").at("max"),
            "exact retained match residual maximum");
      else
        return as_i32(match_epsilon.at("required_complete_max"),
                      "Acb retained match required maximum");
    }();

    std::vector<std::shared_ptr<StoredLocal<Scalar>>> typed_basis;
    std::vector<const LocalSolution<Scalar>*> basis_solutions;
    typed_basis.reserve(basis_owners_.size());
    basis_solutions.reserve(basis_owners_.size());
    std::int32_t materialized_top = kCompleteInfinity;
    const auto receiving_chart = basis_owners_.front()->source_chart();
    const auto receiving_operator =
        basis_owners_.front()->source_operator_identity();
    for (std::size_t column = 0; column < basis_owners_.size(); ++column) {
      auto typed =
          std::dynamic_pointer_cast<StoredLocal<Scalar>>(basis_owners_[column]);
      if (!typed)
        throw std::logic_error(
            "retained plan-match basis coefficient domain changed");
      if (typed->source_chart() != receiving_chart ||
          typed->source_operator_identity() != receiving_operator)
        throw std::logic_error(
            "retained plan-match basis chart provenance changed");
      if (typed->top_valid() < match_work_max)
        throw std::domain_error(
            "retained plan-match basis validity does not cover its matching work window");
      if (column >= weights.size())
        throw std::logic_error(
            "retained plan-match weight count is smaller than its basis");
      const auto basis_valid = std::min(
          typed->top_valid(), typed->solution().epsilon.complete_max);
      const auto shifted_basis_valid = local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(basis_valid) +
              weights[column].min_power(),
          "materialized local top validity");
      const auto weight_valid = local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(typed->solution().epsilon.min_power) +
              weights[column].complete_max(),
          "materialized weight top validity");
      materialized_top = std::min(
          materialized_top, std::min(shifted_basis_valid, weight_valid));
      basis_solutions.push_back(&typed->solution());
      typed_basis.push_back(std::move(typed));
    }
    if (weights.size() != typed_basis.size())
      throw std::logic_error(
          "retained plan-match weight count differs from its basis");
    const auto typed_incoming =
        std::dynamic_pointer_cast<StoredLocal<Scalar>>(incoming_owner_);
    if (!typed_incoming || typed_incoming->top_valid() < match_work_max)
      throw std::domain_error(
          "retained plan-match incoming validity does not cover its matching work window");

    auto solution = materialize_local_basis_weights(
        basis_solutions, weights, result_checkpoint_identity);
    materialized_top = std::min(
        {materialized_top, solution.epsilon.complete_max,
         match_certified_max});
    if (materialized_top < solution.epsilon.min_power)
      throw std::domain_error(
          "plan-match materialization has no valid output epsilon coefficient");
    if (materialized_top < solution.epsilon.complete_max)
      solution = restrict_local_epsilon_frame_strict_lower(
          solution, solution.epsilon.min_power, materialized_top,
          result_checkpoint_identity);
    json::array weight_windows;
    weight_windows.reserve(weights.size());
    for (const auto& weight : weights)
      weight_windows.push_back(json::object{
          {"min", weight.min_power()},
          {"max", weight.complete_max()}});
    json::object derivation{
        {"schema", "diffexp2-retained-plan-match-local-materialization-v1"},
        {"capability", kRetainedPlannedMatchMaterializationCapability},
        {"source_match", handle_},
        {"source_match_checkpoint_identity", checkpoint_identity_},
        {"source_match_provenance_identity",
         required_string(native_match_summary, "provenance_identity")},
        {"planned_hop_provenance_identity", provenance_identity_},
        {"planned_hop", handoff_},
        {"weight_windows", std::move(weight_windows)},
        {"match_certified_complete_max", match_certified_max},
        {"output", json::object{
             {"checkpoint_identity", result_checkpoint_identity},
             {"chart", receiving_chart},
             {"source_operator_identity", receiving_operator},
             {"epsilon", json::object{
                  {"min", solution.epsilon.min_power},
                  {"max", solution.epsilon.complete_max}}},
             {"taylor_complete_max", solution.taylor_complete_max},
             {"dimension", solution.dimension}}},
        {"scope", "single-match-receiving-local"},
        {"coefficient_transport", "native-retained-only"},
        {"whole_arm_complete", false}};
    const auto derivation_identity = json::serialize(
        canonical_json_value(derivation));
    derivation["provenance_identity"] = derivation_identity;
    std::vector<const RegularTaylorTailModelResult*> basis_tail_models;
    basis_tail_models.reserve(typed_basis.size());
    for (const auto& column : typed_basis)
      basis_tail_models.push_back(&column->tail_model());
    auto tail_model = derive_materialized_regular_homogeneous_tail_model(
        basis_solutions, basis_tail_models, weights, solution,
        receiving_operator, checkpoint_identity_);
    const auto elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    NativeLocalDiagnostics diagnostics;
    diagnostics.top_valid = materialized_top;
    diagnostics.kernel_ms = elapsed_ms;
    return make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
        local_handle, receiving_chart, receiving_operator,
        std::move(solution), precision_bits,
        std::vector<PseudoHit<Scalar>>{}, diagnostics, std::nullopt,
        std::move(derivation), std::static_pointer_cast<void>(self),
        std::move(tail_model));
  }

  std::shared_ptr<StoredMatchBase> match_;
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  json::object handoff_;
  double elapsed_ms_ = 0.0;
  std::shared_ptr<StoredTilePlan> plan_owner_;
  std::vector<std::shared_ptr<StoredLocalBase>> basis_owners_;
  std::shared_ptr<StoredLocalBase> incoming_owner_;
  std::atomic<std::uint64_t> materializations_{0};
};

class StoredTransportArmState {
 public:
  StoredTransportArmState(
      std::string handle, std::string checkpoint_identity,
      std::string arm, std::shared_ptr<StoredTilePlan> plan_owner,
      std::shared_ptr<StoredLocalBase> anchor_owner,
      std::vector<std::vector<std::shared_ptr<StoredLocalBase>>> basis_owners,
      std::vector<std::shared_ptr<StoredPlannedMatchHop>> matches,
      std::vector<std::shared_ptr<StoredLocalBase>> tile_sources,
      EpsilonWindow work_epsilon,
      std::int32_t public_required_complete_max,
      std::int32_t match_required_complete_max,
      json::object refinement, double elapsed_ms)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        arm_(std::move(arm)), plan_owner_(std::move(plan_owner)),
        anchor_owner_(std::move(anchor_owner)),
        basis_owners_(std::move(basis_owners)),
        matches_(std::move(matches)),
        tile_sources_(std::move(tile_sources)),
        work_epsilon_(work_epsilon),
        public_required_complete_max_(public_required_complete_max),
        match_required_complete_max_(match_required_complete_max),
        refinement_(std::move(refinement)), elapsed_ms_(elapsed_ms) {
    validate();
    provenance_identity_ = json::serialize(
        canonical_json_value(provenance_record()));
  }

  const std::string& handle() const { return handle_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& provenance_identity() const {
    return provenance_identity_;
  }
  const std::string& arm_name() const { return arm_; }
  double elapsed_ms() const { return elapsed_ms_; }
  const std::shared_ptr<StoredTilePlan>& plan_owner() const {
    return plan_owner_;
  }
  const std::shared_ptr<StoredLocalBase>& anchor_owner() const {
    return anchor_owner_;
  }
  const std::vector<std::vector<std::shared_ptr<StoredLocalBase>>>&
  basis_owners() const {
    return basis_owners_;
  }
  const std::vector<std::shared_ptr<StoredPlannedMatchHop>>& matches() const {
    return matches_;
  }
  const std::vector<std::shared_ptr<StoredLocalBase>>& tile_sources() const {
    return tile_sources_;
  }
  const std::shared_ptr<StoredLocalBase>& final_local() const {
    return tile_sources_.back();
  }

  json::object summary() const {
    return json::object{
        {"transport_state", handle_},
        {"capability", kRetainedTransportArmStateCapability},
        {"native_retained", true},
        {"json_coefficients", 0},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"tile_plan", plan_owner_->handle()},
        {"tile_plan_checkpoint_identity",
         plan_owner_->checkpoint_identity()},
        {"arm", arm_},
        {"matches", matches_.size()},
        {"tiles", tile_sources_.size()},
        {"epsilon", epsilon_record()},
        {"refinement", refinement_},
        {"final_local", local_reference(final_local())},
        {"strong_ownership", json::object{
             {"tile_plan", true}, {"anchor", true},
             {"basis_locals", basis_owner_count()},
             {"matches", matches_.size()},
             {"tile_sources", tile_sources_.size()}}},
        {"elapsed_ms", elapsed_ms_}};
  }

  json::object stats_json() const {
    auto result = summary();
    result["stats_queries"] = stats_queries_.fetch_add(1) + 1;
    return result;
  }

  json::object checkpoint_record() const {
    return json::object{
        {"schema", "diffexp2-retained-transport-arm-state-v1"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"provenance", provenance_record()},
        {"elapsed_ms", elapsed_ms_},
        {"runtime_stats",
         json::object{{"stats_queries", stats_queries_.load()}}}};
  }

  void restore_runtime_stats(std::uint64_t stats_queries) {
    stats_queries_.store(stats_queries);
  }

 private:
  static json::object local_reference(
      const std::shared_ptr<StoredLocalBase>& local) {
    if (!local)
      throw std::logic_error(
          "transport-arm state contains a null local owner");
    return json::object{
        {"local", local->handle()}, {"chart", local->source_chart()},
        {"source_operator_identity", local->source_operator_identity()},
        {"checkpoint_identity", local->checkpoint_identity()},
        {"coefficient_domain", local->scalar_domain()}};
  }

  json::array basis_reference() const {
    json::array result;
    result.reserve(basis_owners_.size());
    for (const auto& basis : basis_owners_) {
      json::array columns;
      columns.reserve(basis.size());
      for (const auto& column : basis)
        columns.push_back(local_reference(column));
      result.push_back(std::move(columns));
    }
    return result;
  }

  json::array match_reference() const {
    json::array result;
    result.reserve(matches_.size());
    for (std::size_t index = 0; index < matches_.size(); ++index)
      result.push_back(json::object{
          {"index", index}, {"match", matches_[index]->handle()},
          {"checkpoint_identity", matches_[index]->checkpoint_identity()},
          {"provenance_identity", matches_[index]->provenance_identity()}});
    return result;
  }

  json::array tile_source_reference() const {
    json::array result;
    result.reserve(tile_sources_.size());
    for (std::size_t index = 0; index < tile_sources_.size(); ++index) {
      auto record = local_reference(tile_sources_[index]);
      record["tile"] = index;
      result.push_back(std::move(record));
    }
    return result;
  }

  json::object epsilon_record() const {
    return json::object{
        {"min", work_epsilon_.min_power},
        {"max", work_epsilon_.complete_max},
        {"required_complete_max", public_required_complete_max_},
        {"match_required_complete_max", match_required_complete_max_}};
  }

  json::object provenance_record() const {
    return json::object{
        {"schema", "diffexp2-retained-native-transport-arm-state-v1"},
        {"checkpoint_identity", checkpoint_identity_},
        {"tile_plan", json::object{
             {"handle", plan_owner_->handle()},
             {"checkpoint_identity", plan_owner_->checkpoint_identity()},
             {"provenance_identity", plan_owner_->provenance_identity()}}},
        {"arm", arm_},
        {"anchor", local_reference(anchor_owner_)},
        {"receiving_basis", basis_reference()},
        {"matches", match_reference()},
        {"tile_sources", tile_source_reference()},
        {"final_local", local_reference(final_local())},
        {"epsilon", epsilon_record()},
        {"refinement", refinement_}};
  }

  std::size_t basis_owner_count() const {
    std::size_t result = 0;
    for (const auto& basis : basis_owners_) {
      if (basis.size() > std::numeric_limits<std::size_t>::max() - result)
        throw std::overflow_error(
            "transport-arm basis owner count overflow");
      result += basis.size();
    }
    return result;
  }

  void validate() const {
    if (handle_.empty() || checkpoint_identity_.empty() ||
        !plan_owner_ || !anchor_owner_ || elapsed_ms_ < 0.0 ||
        !std::isfinite(elapsed_ms_))
      throw std::invalid_argument(
          "retained transport-arm state lost an identity or strong owner");
    const auto& retained = plan_owner_->arm(arm_);
    if (basis_owners_.size() != retained.exact.matches.size() ||
        matches_.size() != retained.exact.matches.size() ||
        tile_sources_.size() != retained.exact.tiles.size() ||
        tile_sources_.size() != matches_.size() + 1)
      throw std::invalid_argument(
          "retained transport-arm state does not reproduce its plan topology");
    if (tile_sources_.empty() || tile_sources_.front().get() !=
                                     anchor_owner_.get())
      throw std::invalid_argument(
          "retained transport-arm state lost its anchor tile source");
    (void)work_epsilon_.width();
    if (public_required_complete_max_ < work_epsilon_.min_power ||
        match_required_complete_max_ < public_required_complete_max_ ||
        match_required_complete_max_ > work_epsilon_.complete_max)
      throw std::invalid_argument(
          "retained transport-arm epsilon contract is inconsistent");
    require_exact_keys(refinement_, {"relative_tolerance", "max_steps"},
                       "retained transport-arm refinement policy");
    if (required_string(refinement_, "relative_tolerance").empty() ||
        as_u32(refinement_.at("max_steps"),
               "retained transport-arm refinement steps") > 32)
      throw std::invalid_argument(
          "retained transport-arm refinement policy is invalid");

    for (std::size_t tile = 0; tile < tile_sources_.size(); ++tile) {
      const auto& source = tile_sources_[tile];
      if (!source)
        throw std::invalid_argument(
            "retained transport-arm state contains a null tile source");
      const auto& exact_tile = retained.exact.tiles[tile];
      const auto& chart = retained.charts.at(exact_tile.chart);
      if (source->source_chart() != chart.handle)
        throw std::invalid_argument(
            "retained transport-arm tile source belongs to a different chart");
      source->require_exact_plan_binding(
          chart.geometry, chart.prescriptions,
          "retained transport-arm tile source");
    }
    for (std::size_t index = 0; index < matches_.size(); ++index) {
      const auto& match = matches_[index];
      const auto& basis = basis_owners_[index];
      if (!match || basis.empty() ||
          std::any_of(basis.begin(), basis.end(),
                      [](const auto& owner) { return owner == nullptr; }) ||
          match->plan_owner().get() != plan_owner_.get() ||
          match->incoming_owner().get() != tile_sources_[index].get() ||
          match->basis_owners().size() != basis.size())
        throw std::invalid_argument(
            "retained transport-arm match lost its exact owner set");
      for (std::size_t column = 0; column < basis.size(); ++column)
        if (match->basis_owners()[column].get() != basis[column].get())
          throw std::invalid_argument(
              "retained transport-arm basis differs from its match owner");
      const auto& handoff = match->handoff();
      if (required_string(handoff, "arm") != arm_ ||
          as_u64(handoff.at("match"),
                 "retained transport-arm match index") != index)
        throw std::invalid_argument(
            "retained transport-arm match provenance is out of order");
      const auto& next = tile_sources_[index + 1];
      if (!next->retained_derivation().has_value() ||
          next->retained_derivation_owner().get() != match.get() ||
          required_string(*next->retained_derivation(), "source_match") !=
              match->handle())
        throw std::invalid_argument(
            "retained transport-arm tile source is not materialized from its match");
    }
  }

  std::string handle_;
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::string arm_;
  std::shared_ptr<StoredTilePlan> plan_owner_;
  std::shared_ptr<StoredLocalBase> anchor_owner_;
  std::vector<std::vector<std::shared_ptr<StoredLocalBase>>> basis_owners_;
  std::vector<std::shared_ptr<StoredPlannedMatchHop>> matches_;
  std::vector<std::shared_ptr<StoredLocalBase>> tile_sources_;
  EpsilonWindow work_epsilon_;
  std::int32_t public_required_complete_max_ = 0;
  std::int32_t match_required_complete_max_ = 0;
  json::object refinement_;
  double elapsed_ms_ = 0.0;
  mutable std::atomic<std::uint64_t> stats_queries_{0};
};

class StoredLineResult {
 public:
  StoredLineResult(std::string handle, std::string checkpoint_identity,
                   std::string provenance_identity, std::string arm,
                   std::size_t tile_index, json::object interval,
                   std::string source_checkpoint,
                   StoredLineIntegral result, double elapsed_ms,
                   std::shared_ptr<StoredTilePlan> plan_owner,
                   std::shared_ptr<StoredLocalBase> local_owner)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        arm_(std::move(arm)), tile_index_(tile_index),
        interval_(std::move(interval)),
        source_checkpoint_(std::move(source_checkpoint)),
        result_(std::move(result)), elapsed_ms_(elapsed_ms),
        plan_owner_(std::move(plan_owner)),
        local_owners_{std::move(local_owner)} {}

  StoredLineResult(std::string handle, std::string checkpoint_identity,
                   std::string provenance_identity,
                   StoredLineIntegral result, double elapsed_ms,
                   std::shared_ptr<StoredTilePlan> plan_owner,
                   std::vector<std::shared_ptr<StoredLocalBase>> local_owners,
                   json::object aggregate_provenance)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        result_(std::move(result)), elapsed_ms_(elapsed_ms),
        plan_owner_(std::move(plan_owner)),
        local_owners_(std::move(local_owners)),
        aggregate_provenance_(std::move(aggregate_provenance)) {
    if (!plan_owner_ || local_owners_.empty() ||
        std::any_of(local_owners_.begin(), local_owners_.end(),
                    [](const auto& owner) { return owner == nullptr; }))
      throw std::invalid_argument(
          "retained line aggregate requires its plan and local owners");
    if (json::serialize(canonical_json_value(*aggregate_provenance_)) !=
        provenance_identity_)
      throw std::invalid_argument(
          "retained line aggregate provenance identity is inconsistent");
  }

  const std::string& handle() const { return handle_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& provenance_identity() const {
    return provenance_identity_;
  }
  double elapsed_ms() const { return elapsed_ms_; }
  const StoredLineIntegral& result() const { return result_; }
  bool is_aggregate() const { return aggregate_provenance_.has_value(); }

  const std::shared_ptr<StoredTilePlan>& plan_owner() const {
    return plan_owner_;
  }
  const std::shared_ptr<StoredLocalBase>& local_owner() const {
    return local_owners_.front();
  }
  const std::vector<std::shared_ptr<StoredLocalBase>>& local_owners() const {
    return local_owners_;
  }

  json::object summary() const {
    const auto& diagnostics = result_.diagnostics;
    json::object output{
        {"line", handle_},
        {"capability", is_aggregate()
             ? kRetainedLineAggregateCapability
             : result_.scope ==
                           LineIntegrationScope::FullLocalWithCertifiedTail
             ? kRetainedCertifiedLineCapability
             : kRetainedStoredLineCapability},
        {"native_retained", true}, {"json_coefficients", 0},
        {"scope", line_integration_scope_name(result_.scope)},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"dimension", result_.value.dimension},
        {"epsilon_min", result_.value.epsilon.min_power},
        {"epsilon_max", result_.value.epsilon.complete_max},
        {"effective_rim", result_.imaginary_sign.has_value()
             ? json::value(*result_.imaginary_sign) : json::value(nullptr)},
        {"error", encode_error_envelope_summary(result_.value.error)},
        {"diagnostics", json::object{
             {"input_monomial_cells", diagnostics.input_monomial_cells},
             {"grouped_monomials", diagnostics.grouped_monomials},
             {"zero_groups_skipped", diagnostics.zero_groups_skipped},
             {"cancelled_divergent_groups",
              diagnostics.cancelled_divergent_groups},
             {"primitive_evaluations", diagnostics.primitive_evaluations},
             {"primitive_component_applications",
              diagnostics.primitive_component_applications},
             {"primitive_component_reuses",
              diagnostics.primitive_component_reuses},
             {"has_center_endpoint", diagnostics.has_center_endpoint},
             {"tail_certificate_requested",
              diagnostics.tail_certificate_requested},
             {"tail_certificate_status",
              diagnostics.tail_certificate_status},
             {"tail_witness_radius_exact",
              diagnostics.tail_witness_radius_exact.empty()
                  ? json::value(nullptr)
                  : json::value(
                        diagnostics.tail_witness_radius_exact)},
             {"detail", diagnostics.detail}}},
        {"elapsed_ms", elapsed_ms_}};
    if (aggregate_provenance_.has_value()) {
      output["source"] = aggregate_provenance_->at("source");
      output["arm"] = aggregate_provenance_->at("arm");
      output["tile"] = nullptr;
      output["interval"] = aggregate_provenance_->at("interval");
      output["aggregate"] = aggregate_provenance_->at("aggregate");
    } else {
      const auto& owner = local_owners_.front();
      output["source"] = json::object{
          {"tile_plan", plan_owner_->handle()},
          {"tile_plan_checkpoint_identity",
           plan_owner_->checkpoint_identity()},
          {"local", owner->handle()},
          {"chart", owner->source_chart()},
          {"local_checkpoint_identity", source_checkpoint_}};
      output["arm"] = arm_;
      output["tile"] = tile_index_;
      output["interval"] = interval_;
    }
    return output;
  }

  json::object stats_json() const {
    auto result = summary();
    std::lock_guard<std::mutex> lock(stats_mutex_);
    result["exports"] = exports_;
    result["export_ms"] = export_ms_;
    return result;
  }

  json::object export_values(const std::string& expected_checkpoint,
                             int output_digits) {
    if (expected_checkpoint != checkpoint_identity_)
      throw std::invalid_argument(
          "line export checkpoint identity does not match retained result");
    const auto started = std::chrono::steady_clock::now();
    auto value = encode_epsilon_vector(result_.value, output_digits);
    const auto elapsed = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    {
      std::lock_guard<std::mutex> lock(stats_mutex_);
      ++exports_;
      export_ms_ += elapsed;
    }
    return json::object{
        {"line", handle_}, {"checkpoint_identity", checkpoint_identity_},
        {"compatibility_export", true},
        {"scope", line_integration_scope_name(result_.scope)},
        {"error_guarantee",
         error_guarantee_name(result_.value.error.guarantee)},
        {"json_coefficients", value.at("coefficients").as_array().size()},
        {"value", std::move(value)}, {"elapsed_ms", elapsed}};
  }

  json::object checkpoint_record() const {
    const auto& diagnostics = result_.diagnostics;
    std::lock_guard<std::mutex> lock(stats_mutex_);
    if (aggregate_provenance_.has_value()) {
      return json::object{
          {"schema", "diffexp2-retained-line-aggregate-v1"},
          {"handle", handle_},
          {"checkpoint_identity", checkpoint_identity_},
          {"provenance_identity", provenance_identity_},
          {"provenance", *aggregate_provenance_},
          {"result",
           json::object{
               {"value", checkpoint_epsilon_vector_record(result_.value)},
               {"scope", line_integration_scope_name(result_.scope)},
               {"imaginary_sign", result_.imaginary_sign.has_value()
                    ? json::value(*result_.imaginary_sign)
                    : json::value(nullptr)},
               {"diagnostics",
                json::object{
                    {"input_monomial_cells",
                     diagnostics.input_monomial_cells},
                    {"grouped_monomials", diagnostics.grouped_monomials},
                    {"zero_groups_skipped",
                     diagnostics.zero_groups_skipped},
                    {"cancelled_divergent_groups",
                     diagnostics.cancelled_divergent_groups},
                    {"primitive_evaluations",
                     diagnostics.primitive_evaluations},
                    {"primitive_component_applications",
                     diagnostics.primitive_component_applications},
                    {"primitive_component_reuses",
                     diagnostics.primitive_component_reuses},
                    {"has_center_endpoint",
                     diagnostics.has_center_endpoint},
                    {"tail_certificate_requested",
                     diagnostics.tail_certificate_requested},
                    {"tail_certificate_status",
                     diagnostics.tail_certificate_status},
                    {"tail_witness_radius_exact",
                     diagnostics.tail_witness_radius_exact},
                    {"detail", diagnostics.detail}}}}},
          {"elapsed_ms", elapsed_ms_},
          {"runtime_stats",
           json::object{{"exports", exports_}, {"export_ms", export_ms_}}}};
    }
    const auto& owner = local_owners_.front();
    return json::object{
        {"schema", "diffexp2-retained-line-result-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"arm", arm_},
        {"tile", tile_index_},
        {"interval", interval_},
        {"source",
         json::object{
             {"tile_plan", plan_owner_->handle()},
             {"tile_plan_checkpoint_identity",
              plan_owner_->checkpoint_identity()},
             {"local", owner->handle()},
             {"chart", owner->source_chart()},
             {"source_operator_identity",
              owner->source_operator_identity()},
             {"local_checkpoint_identity", source_checkpoint_},
             {"coefficient_domain", owner->scalar_domain()},
             {"analytic_metadata",
              owner->exact_analytic_metadata()}}},
        {"result",
         json::object{
             {"value", checkpoint_epsilon_vector_record(result_.value)},
             {"scope", line_integration_scope_name(result_.scope)},
             {"imaginary_sign", result_.imaginary_sign.has_value()
                  ? json::value(*result_.imaginary_sign)
                  : json::value(nullptr)},
             {"diagnostics",
              json::object{
                  {"input_monomial_cells",
                   diagnostics.input_monomial_cells},
                  {"grouped_monomials", diagnostics.grouped_monomials},
                  {"zero_groups_skipped", diagnostics.zero_groups_skipped},
                  {"cancelled_divergent_groups",
                   diagnostics.cancelled_divergent_groups},
                  {"primitive_evaluations",
                   diagnostics.primitive_evaluations},
                  {"primitive_component_applications",
                   diagnostics.primitive_component_applications},
                  {"primitive_component_reuses",
                   diagnostics.primitive_component_reuses},
                  {"has_center_endpoint",
                   diagnostics.has_center_endpoint},
                  {"tail_certificate_requested",
                   diagnostics.tail_certificate_requested},
                  {"tail_certificate_status",
                   diagnostics.tail_certificate_status},
                  {"tail_witness_radius_exact",
                   diagnostics.tail_witness_radius_exact.empty()
                       ? json::value(nullptr)
                       : json::value(
                             diagnostics.tail_witness_radius_exact)},
                  {"detail", diagnostics.detail}}}}},
        {"elapsed_ms", elapsed_ms_},
        {"runtime_stats",
         json::object{{"exports", exports_}, {"export_ms", export_ms_}}}};
  }

  void restore_runtime_stats(std::uint64_t exports, double export_ms) {
    std::lock_guard<std::mutex> lock(stats_mutex_);
    exports_ = exports;
    export_ms_ = export_ms;
  }

 private:
  std::string handle_;
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::string arm_;
  std::size_t tile_index_ = 0;
  json::object interval_;
  std::string source_checkpoint_;
  StoredLineIntegral result_;
  double elapsed_ms_ = 0.0;
  std::shared_ptr<StoredTilePlan> plan_owner_;
  std::vector<std::shared_ptr<StoredLocalBase>> local_owners_;
  std::optional<json::object> aggregate_provenance_;
  mutable std::mutex stats_mutex_;
  std::uint64_t exports_ = 0;
  double export_ms_ = 0.0;
};

struct SolverSession {
  std::string handle;
  std::string domain;
  slong precision_bits = 256;
  int output_digits = 50;
  std::vector<std::string> symbols;
  std::string analytic_identity;
  std::size_t chart_capacity = 256;
  std::size_t local_capacity = 1024;
  std::size_t scc_capacity = 128;
  std::size_t match_capacity = 1024;
  std::size_t endpoint_capacity = 1024;
  std::size_t tile_plan_capacity = 256;
  std::size_t transport_state_capacity = 256;
  std::size_t line_result_capacity = 2048;
  std::uint64_t next_chart = 1;
  std::uint64_t next_local = 1;
  std::uint64_t next_scc = 1;
  std::uint64_t next_match = 1;
  std::uint64_t next_endpoint = 1;
  std::uint64_t next_tile_plan = 1;
  std::uint64_t next_transport_state = 1;
  std::uint64_t next_line_result = 1;
  std::size_t pending_local_solves = 0;
  std::size_t pending_matches = 0;
  std::size_t pending_endpoint_limits = 0;
  std::size_t pending_tile_plans = 0;
  std::size_t pending_transport_states = 0;
  std::size_t pending_line_integrations = 0;
  std::uint64_t total_local_solves = 0;
  std::uint64_t total_scc_column_solves = 0;
  std::uint64_t total_local_matches = 0;
  std::uint64_t total_endpoint_limits = 0;
  std::uint64_t total_endpoint_exports = 0;
  std::uint64_t total_tile_plans = 0;
  std::uint64_t total_transport_arm_marches = 0;
  std::uint64_t total_line_integrations = 0;
  std::uint64_t total_line_exports = 0;
  double total_local_run_parse_ms = 0.0;
  double total_local_kernel_ms = 0.0;
  double total_local_match_ms = 0.0;
  std::uint64_t checkpoint_generation = 0;
  std::uint64_t checkpoint_restore_count = 0;
  std::string restored_from_checkpoint_identity;
  double total_endpoint_limit_ms = 0.0;
  double total_endpoint_export_ms = 0.0;
  double total_tile_plan_ms = 0.0;
  double total_transport_arm_ms = 0.0;
  double total_line_integration_ms = 0.0;
  double total_line_export_ms = 0.0;
  bool closed = false;
  mutable std::mutex mutex;
  std::unordered_map<std::string, std::shared_ptr<PreparedChartBase>> charts;
  std::unordered_map<std::string, std::string> handles_by_key;
  std::unordered_map<std::string, std::shared_ptr<StoredLocalBase>> locals;
  std::unordered_map<std::string, std::shared_ptr<StoredMatchBase>> matches;
  std::unordered_map<std::string, std::shared_ptr<StoredEndpointResult>>
      endpoints;
  std::unordered_map<std::string, std::shared_ptr<StoredTilePlan>> tile_plans;
  std::unordered_map<std::string, std::shared_ptr<StoredTransportArmState>>
      transport_states;
  std::unordered_map<std::string, std::shared_ptr<StoredLineResult>>
      line_results;
  std::unordered_map<std::string, std::shared_ptr<CompositeSCCChartBase>> sccs;
  std::unordered_map<std::string, std::string> scc_handles_by_key;
};

Rational parse_exact_path_rational(const json::value& raw,
                                   const char* label) {
  if (!raw.is_string() || raw.as_string().empty())
    throw std::invalid_argument(std::string(label) +
                                " must be a nonempty exact rational string");
  try {
    return Rational(std::string(raw.as_string()));
  } catch (const std::invalid_argument&) {
    throw std::invalid_argument(std::string(label) +
                                " is not an exact rational number");
  }
}

std::vector<Rational> parse_exact_path_points(const json::value& raw,
                                              const char* label) {
  std::vector<Rational> points;
  for (const auto& value : as_array(raw, label))
    points.push_back(parse_exact_path_rational(value, label));
  return points;
}

ExactPathTopology parse_exact_path_topology(const json::value& raw) {
  const auto& object = as_object(raw, "native exact path topology");
  require_exact_keys(object,
      {"singular_points", "boundary_points", "complex_projections",
       "branch_sheets"}, "native exact path topology");
  ExactPathTopology topology;
  topology.singular_points = parse_exact_path_points(
      object.at("singular_points"), "path singular point");
  topology.boundary_points = parse_exact_path_points(
      object.at("boundary_points"), "path boundary point");
  for (const auto& raw_projection : as_array(
           object.at("complex_projections"), "path complex projections")) {
    const auto& projection = as_object(
        raw_projection, "path complex projection");
    require_exact_keys(projection,
        {"source_identity", "real_part_exact",
         "imaginary_magnitude_exact", "retain_minus_imaginary",
         "retain_real_part", "retain_plus_imaginary"},
        "path complex projection");
    topology.complex_projections.push_back(ExactComplexProjection{
        required_string(projection, "source_identity"),
        parse_exact_path_rational(projection.at("real_part_exact"),
                                  "projection real part"),
        parse_exact_path_rational(
            projection.at("imaginary_magnitude_exact"),
            "projection imaginary magnitude"),
        projection.at("retain_minus_imaginary").as_bool(),
        projection.at("retain_real_part").as_bool(),
        projection.at("retain_plus_imaginary").as_bool()});
  }
  for (const auto& raw_sheet : as_array(
           object.at("branch_sheets"), "path branch sheets")) {
    const auto& sheet = as_object(raw_sheet, "path branch sheet");
    require_exact_keys(sheet, {"factor_exact", "sign"},
                       "path branch sheet");
    topology.branch_sheets.push_back(ExactBranchSheet{
        required_string(sheet, "factor_exact"),
        as_i32(sheet.at("sign"), "path branch sign")});
  }
  return topology;
}

std::string retained_plan_owner_handle(
    const RetainedPlanChartBinding::Owner& owner) {
  return std::visit(
      [](const auto& typed) {
        if (typed == nullptr)
          throw std::logic_error(
              "native tile planning received a null retained owner");
        return typed->handle();
      },
      owner);
}

std::string retained_plan_owner_identity(
    const RetainedPlanChartBinding::Owner& owner) {
  return std::visit(
      [](const auto& typed) {
        if (typed == nullptr)
          throw std::logic_error(
              "native tile planning received a null retained owner");
        return typed->exact_identity();
      },
      owner);
}

std::string retained_plan_owner_geometry_record(
    const RetainedPlanChartBinding::Owner& owner) {
  return std::visit(
      [](const auto& typed) -> std::string {
        using Owner = typename std::decay_t<decltype(typed)>::element_type;
        if (typed == nullptr)
          throw std::logic_error(
              "native tile planning received a null retained owner");
        if constexpr (std::is_same_v<Owner, PreparedChartBase>) {
          if (!typed->geometry_record().has_value())
            throw std::invalid_argument(
                "native tile planning requires retained exact chart geometry");
          return *typed->geometry_record();
        } else {
          return typed->geometry_record();
        }
      },
      owner);
}

RetainedPlanChartBinding bind_plan_chart(
    const RetainedPlanChartBinding::Owner& owner,
    const ExactPathTopology& topology) {
  const auto handle = retained_plan_owner_handle(owner);
  const auto exact_identity = retained_plan_owner_identity(owner);
  const auto geometry_value = json::parse(
      retained_plan_owner_geometry_record(owner));
  const auto& geometry = as_object(
      geometry_value, "retained native tile chart geometry");
  if (geometry.at("infinite_radius").as_bool())
    throw std::invalid_argument(
        "native exact tile planning currently requires finite chart radii");

  RetainedPlanChartBinding binding;
  binding.handle = handle;
  binding.exact_identity = exact_identity;
  binding.owner = owner;
  binding.geometry.identity = exact_identity;
  binding.geometry.center = parse_exact_path_rational(
      geometry.at("center_exact"), "tile chart center");
  binding.geometry.scale = parse_exact_path_rational(
      geometry.at("scale_exact"), "tile chart scale");
  binding.geometry.radius = parse_exact_path_rational(
      geometry.at("radius_exact"), "tile chart radius");
  binding.geometry.singular_center = std::any_of(
      topology.singular_points.begin(), topology.singular_points.end(),
      [&](const Rational& point) { return point == binding.geometry.center; });
  for (const auto& raw_prescription : as_array(
           geometry.at("prescriptions"), "tile chart prescriptions")) {
    const auto& prescription = as_object(
        raw_prescription, "tile chart prescription");
    binding.prescriptions.push_back(Prescription{
        required_string(prescription, "factor_exact"),
        as_i32(prescription.at("sign"), "tile prescription sign"),
        as_u32(prescription.at("multiplicity"),
               "tile prescription multiplicity"),
        as_i32(prescription.at("leading_coefficient_sign"),
               "tile prescription leading coefficient sign")});
  }
  // Every exact prepared prescription must be represented on the path.  The
  // topology may additionally retain factors which are inactive in this
  // chart, but it may never silently alter or drop an active sheet.
  for (const auto& prescription : binding.prescriptions) {
    const auto found = std::find_if(
        topology.branch_sheets.begin(), topology.branch_sheets.end(),
        [&](const ExactBranchSheet& sheet) {
          return sheet.factor_exact == prescription.factor_exact;
        });
    if (found == topology.branch_sheets.end() ||
        found->imaginary_sign != prescription.sign)
      throw std::invalid_argument(
          "native tile topology does not reproduce a prepared chart branch prescription");
  }
  (void)exact_plan_rim(binding.prescriptions, binding.geometry.scale);
  return binding;
}

std::vector<std::string> parse_plan_chart_handles(const json::object& arm) {
  std::vector<std::string> handles;
  for (const auto& raw : as_array(arm.at("charts"), "tile arm charts")) {
    if (!raw.is_string() || raw.as_string().empty())
      throw std::invalid_argument(
          "native tile arm chart handles must be nonempty strings");
    handles.emplace_back(raw.as_string());
  }
  if (handles.empty())
    throw std::invalid_argument("native tile arm requires at least one chart");
  return handles;
}

std::pair<ExactArmRequest, std::vector<RetainedPlanChartBinding>>
parse_retained_arm_request(
    const json::object& arm,
    const std::vector<RetainedPlanChartBinding::Owner>& charts) {
  require_exact_keys(arm, {"from_exact", "to_exact", "charts", "topology"},
                     "native tile arm");
  const auto handles = parse_plan_chart_handles(arm);
  if (handles.size() != charts.size())
    throw std::invalid_argument(
        "resolved native tile chart count differs from its request");
  ExactArmRequest request;
  request.from = parse_exact_path_rational(arm.at("from_exact"),
                                           "tile arm start");
  request.to = parse_exact_path_rational(arm.at("to_exact"),
                                         "tile arm end");
  request.topology = parse_exact_path_topology(arm.at("topology"));
  std::vector<RetainedPlanChartBinding> bindings;
  bindings.reserve(charts.size());
  request.charts.reserve(charts.size());
  for (std::size_t index = 0; index < charts.size(); ++index) {
    if (retained_plan_owner_handle(charts[index]) != handles[index])
      throw std::logic_error("resolved tile chart handle changed");
    auto binding = bind_plan_chart(charts[index], request.topology);
    request.charts.push_back(binding.geometry);
    bindings.push_back(std::move(binding));
  }
  return {std::move(request), std::move(bindings)};
}

std::shared_ptr<StoredTilePlan> build_tile_plan(
    const std::string& handle, const json::object& request,
    const std::vector<RetainedPlanChartBinding::Owner>& lower_charts,
    const std::vector<RetainedPlanChartBinding::Owner>& upper_charts) {
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  const auto division_order = as_u32(
      request.at("division_order"), "native tile division order");
  auto [lower_request, lower_bindings] = parse_retained_arm_request(
      as_object(request.at("lower"), "lower native tile arm"), lower_charts);
  auto [upper_request, upper_bindings] = parse_retained_arm_request(
      as_object(request.at("upper"), "upper native tile arm"), upper_charts);
  if (lower_bindings.front().handle != upper_bindings.front().handle)
    throw std::invalid_argument(
        "independent native tile arms must share one retained anchor chart");

  ExactPathPlanOptions options;
  options.division_order = division_order;
  const auto started = std::chrono::steady_clock::now();
  auto exact = plan_exact_independent_arms(
      lower_request, upper_request, options);
  RetainedArmPlan lower{std::move(exact.lower),
                        std::move(lower_bindings)};
  RetainedArmPlan upper{std::move(exact.upper),
                        std::move(upper_bindings)};
  json::object provenance{
      {"schema", "diffexp2-retained-exact-independent-arm-tile-plan-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"division_order", division_order},
      {"lower", encode_retained_arm(lower)},
      {"upper", encode_retained_arm(upper)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredTilePlan>(
      handle, checkpoint_identity, provenance_identity, division_order,
      std::move(lower), std::move(upper), elapsed);
}

std::shared_ptr<StoredTilePlan> build_single_arm_tile_plan(
    const std::string& handle, const json::object& request,
    const std::vector<RetainedPlanChartBinding::Owner>& charts) {
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  if (checkpoint_identity.empty())
    throw std::invalid_argument(
        "native single-arm tile-plan checkpoint identity cannot be empty");
  const auto division_order = as_u32(
      request.at("division_order"), "native single-arm tile division order");
  auto [arm_request, bindings] = parse_retained_arm_request(
      as_object(request.at("arm"), "native single tile arm"), charts);
  ExactPathPlanOptions options;
  options.division_order = division_order;
  const auto started = std::chrono::steady_clock::now();
  auto exact = plan_exact_arm(arm_request, options);
  const std::string arm_name = exact.direction < 0 ? "lower" : "upper";
  RetainedArmPlan retained{std::move(exact), std::move(bindings)};
  json::object provenance{
      {"schema", kRetainedSingleArmTilePlanProvenanceSchema},
      {"checkpoint_identity", checkpoint_identity},
      {"division_order", division_order},
      {"arm_name", arm_name},
      {"arm", encode_retained_arm(retained)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredTilePlan>(
      handle, checkpoint_identity, provenance_identity, division_order,
      arm_name, std::move(retained), elapsed);
}

json::value optional_plan_rim_json(
    const std::optional<std::int32_t>& rim) {
  return rim.has_value() ? json::value(*rim) : json::value(nullptr);
}

struct ResolvedPlannedEndpointBinding {
  json::object source;
  std::int32_t approach_direction = 0;
  std::optional<std::int32_t> rim;
};

ResolvedPlannedEndpointBinding resolve_planned_endpoint_binding(
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::string& arm_name,
    const std::shared_ptr<StoredLocalBase>& local) {
  if (!plan || !local)
    throw std::invalid_argument(
        "plan-bound endpoint evaluation requires retained plan and local owners");
  const auto& arm = plan->arm(arm_name);
  if (arm.exact.tiles.empty() || arm.charts.empty())
    throw std::invalid_argument(
        "plan-bound endpoint arm has no final tile/chart");
  const auto final_tile_index = arm.exact.tiles.size() - 1;
  const auto& final_tile = arm.exact.tiles.back();
  if (final_tile.chart >= arm.charts.size())
    throw std::logic_error(
        "plan-bound endpoint final tile has an invalid chart index");
  const auto& final_chart = arm.charts[final_tile.chart];
  if (!(final_tile.physical_end == arm.exact.to) ||
      !(final_chart.geometry.center == arm.exact.to) ||
      !final_tile.local_end.is_zero())
    throw std::invalid_argument(
        "plan-bound endpoint requires the final tile to end at the exact center of its retained final chart");
  if (arm.exact.direction != -1 && arm.exact.direction != 1)
    throw std::logic_error(
        "plan-bound endpoint arm has an invalid exact direction");
  if (final_chart.geometry.scale.is_zero())
    throw std::invalid_argument(
        "plan-bound endpoint final chart has a zero exact scale");
  if (local->source_chart() != final_chart.handle)
    throw std::invalid_argument(
        "plan-bound endpoint local does not name the retained final chart");
  local->require_exact_plan_binding(
      final_chart.geometry, final_chart.prescriptions,
      "plan-bound endpoint final local");

  ResolvedPlannedEndpointBinding resolved;
  resolved.approach_direction =
      -arm.exact.direction * final_chart.geometry.scale.sign();
  resolved.rim = exact_plan_rim(
      final_chart.prescriptions, final_chart.geometry.scale);
  resolved.source = json::object{
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
      {"tile_plan_provenance_identity", plan->provenance_identity()},
      {"arm", arm_name},
      {"endpoint_exact", arm.exact.to.str()},
      {"direction", arm.exact.direction},
      {"final_tile", final_tile_index},
      {"final_chart_index", final_tile.chart},
      {"final_chart", final_chart.handle},
      {"final_chart_identity", final_chart.exact_identity},
      {"local", local->handle()},
      {"chart", local->source_chart()},
      {"source_operator_identity", local->source_operator_identity()},
      {"checkpoint_identity", local->checkpoint_identity()},
      {"coefficient_domain", local->scalar_domain()},
      {"prescriptions", encode_plan_prescriptions(
           final_chart.prescriptions)}};
  return resolved;
}

json::object planned_endpoint_provenance(
    const std::string& checkpoint_identity,
    const ResolvedPlannedEndpointBinding& binding,
    const std::string& cancellation_mode,
    const json::object& analytic_metadata) {
  return json::object{
      {"schema",
       "diffexp2-retained-native-plan-bound-endpoint-sector-limit-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"source", binding.source},
      {"approach_direction", binding.approach_direction},
      {"rim", optional_plan_rim_json(binding.rim)},
      {"cancellation", json::object{{"mode", cancellation_mode}}},
      {"analytic_metadata", analytic_metadata}};
}

std::shared_ptr<StoredEndpointResult> build_planned_endpoint_limit(
    const std::string& endpoint_handle, const json::object& request,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& local) {
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  const auto expected_plan_checkpoint = required_string(
      request, "tile_plan_checkpoint_identity");
  const auto expected_source_checkpoint = required_string(
      request, "source_checkpoint_identity");
  if (checkpoint_identity.empty() || expected_plan_checkpoint.empty() ||
      expected_source_checkpoint.empty())
    throw std::invalid_argument(
        "plan-bound endpoint checkpoint identities must be nonempty");
  if (expected_plan_checkpoint != plan->checkpoint_identity())
    throw std::invalid_argument(
        "plan-bound endpoint tile-plan checkpoint identity is stale or mismatched");
  if (expected_source_checkpoint != local->checkpoint_identity())
    throw std::invalid_argument(
        "plan-bound endpoint source checkpoint identity is stale or mismatched");
  const auto arm_name = required_string(request, "arm");
  const auto binding = resolve_planned_endpoint_binding(
      plan, arm_name, local);
  const auto cancellation = parse_endpoint_cancellation_policy(request);
  auto analytic_metadata = local->exact_analytic_metadata();
  const auto provenance = planned_endpoint_provenance(
      checkpoint_identity, binding, cancellation.cancellation_mode,
      analytic_metadata);
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));

  EndpointLimitOptions options;
  options.approach_direction = binding.approach_direction;
  options.imaginary_sign = binding.rim;
  options.allow_certified_numeric_cancellation =
      cancellation.allow_certified_numeric_cancellation;
  const auto started = std::chrono::steady_clock::now();
  auto result = local->endpoint_limit(options);
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredEndpointResult>(
      endpoint_handle, checkpoint_identity, provenance_identity,
      local->handle(), local->source_chart(),
      local->source_operator_identity(), expected_source_checkpoint,
      local->scalar_domain(), binding.approach_direction, std::nullopt,
      cancellation.cancellation_mode, std::move(analytic_metadata),
      std::move(result), elapsed, binding.source, binding.rim,
      plan, local);
}

json::object planned_match_handoff_record(
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::string& arm_name, std::size_t match_index,
    const std::vector<std::string>& basis_handles,
    const std::vector<std::shared_ptr<StoredLocalBase>>& basis,
    const std::string& incoming_handle,
    const std::shared_ptr<StoredLocalBase>& incoming,
    const std::string& result_checkpoint,
    const json::object& native_summary) {
  if (!plan || !incoming || basis.empty() ||
      basis.size() != basis_handles.size())
    throw std::invalid_argument(
        "planned match handoff record requires its complete owner set");
  const auto& arm = plan->arm(arm_name);
  if (match_index >= arm.exact.matches.size())
    throw std::invalid_argument(
        "planned match handoff index is outside its retained arm");
  const auto& exact_match = arm.exact.matches[match_index];
  const auto& producing = arm.charts.at(exact_match.producing_chart);
  const auto& receiving = arm.charts.at(exact_match.receiving_chart);
  const auto producing_rim = exact_plan_rim(
      producing.prescriptions, producing.geometry.scale);
  const auto receiving_rim = exact_plan_rim(
      receiving.prescriptions, receiving.geometry.scale);
  json::array basis_sources;
  basis_sources.reserve(basis.size());
  for (std::size_t column = 0; column < basis.size(); ++column)
    basis_sources.push_back(json::object{
        {"column", column}, {"local", basis_handles[column]},
        {"checkpoint_identity", basis[column]->checkpoint_identity()},
        {"source_operator_identity",
         basis[column]->source_operator_identity()}});
  return json::object{
      {"schema", "diffexp2-retained-exact-plan-match-hop-v1"},
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
      {"tile_plan_provenance_identity", plan->provenance_identity()},
      {"arm", arm_name}, {"match", match_index},
      {"geometry", encode_plan_match(arm, match_index)},
      {"producing", json::object{
           {"tile", match_index},
           {"chart", producing.handle},
           {"chart_identity", producing.exact_identity},
           {"local_point_exact", exact_match.producing_local.str()},
           {"effective_rim", optional_plan_rim_json(producing_rim)},
           {"prescriptions",
            encode_plan_prescriptions(producing.prescriptions)},
           {"incoming", json::object{
                {"local", incoming_handle},
                {"checkpoint_identity", incoming->checkpoint_identity()},
                {"source_operator_identity",
                 incoming->source_operator_identity()}}}}},
      {"receiving", json::object{
           {"tile", match_index + 1},
           {"chart", receiving.handle},
           {"chart_identity", receiving.exact_identity},
           {"local_point_exact", exact_match.receiving_local.str()},
           {"effective_rim", optional_plan_rim_json(receiving_rim)},
           {"prescriptions",
            encode_plan_prescriptions(receiving.prescriptions)},
           {"basis", std::move(basis_sources)}}},
      {"result_checkpoint_identity", result_checkpoint},
      {"native_match_provenance_identity",
       required_string(native_summary, "provenance_identity")},
      {"advance", json::object{
           {"scope", "single-match-handoff"},
           {"state", "retained-receiving-basis-weights"},
           {"source_tile", match_index},
           {"receiving_tile", match_index + 1},
           {"whole_arm_complete", false}}}};
}

struct NativeAcbSaturationBinding {
  std::string request_key;
  json::object request;
};

NativeAcbSaturationBinding native_acb_saturation_binding(
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::string& session_configuration_identity,
    const std::string& arm,
    std::size_t match_index, const std::string& match_checkpoint_identity);

std::shared_ptr<StoredPlannedMatchHop> build_planned_match_hop(
    const std::string& match_handle, const json::object& request,
    const std::string& domain, slong precision_bits,
    const std::string& active_session_configuration_identity,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::vector<std::string>& basis_handles,
    const std::vector<std::shared_ptr<StoredLocalBase>>& basis,
    const std::string& incoming_handle,
    const std::shared_ptr<StoredLocalBase>& incoming) {
  const auto started = std::chrono::steady_clock::now();
  const auto arm_name = required_string(request, "arm");
  const auto match_index = static_cast<std::size_t>(
      as_u64(request.at("match"), "planned local match index"));
  const auto& arm = plan->arm(arm_name);
  if (match_index >= arm.exact.matches.size())
    throw std::invalid_argument(
        "planned local match index is out of range");
  const auto& exact_match = arm.exact.matches[match_index];
  const auto& producing = arm.charts.at(exact_match.producing_chart);
  const auto& receiving = arm.charts.at(exact_match.receiving_chart);
  if (basis.empty())
    throw std::invalid_argument("planned local match basis cannot be empty");

  // The plan, rather than the caller, binds every coordinate, chart,
  // prescription, rim and source checkpoint passed to the existing matching
  // kernels.  Locals must reproduce the prepared chart snapshot exactly.
  incoming->require_exact_plan_binding(
      producing.geometry, producing.prescriptions,
      "planned incoming " + incoming_handle);
  for (std::size_t column = 0; column < basis.size(); ++column)
    basis[column]->require_exact_plan_binding(
        receiving.geometry, receiving.prescriptions,
        "planned basis " + basis_handles[column]);

  json::array basis_checkpoints;
  basis_checkpoints.reserve(basis.size());
  for (const auto& local : basis)
    basis_checkpoints.emplace_back(local->checkpoint_identity());
  const auto result_checkpoint = required_string(
      request, "checkpoint_identity");
  if (result_checkpoint.empty())
    throw std::invalid_argument(
        "planned local match checkpoint identity cannot be empty");

  json::object kernel_request{
      {"basis", [&]() {
         json::array values;
         for (const auto& handle : basis_handles) values.emplace_back(handle);
         return values;
       }()},
      {"incoming", incoming_handle},
      {"basis_chart", receiving.handle},
      {"incoming_chart", producing.handle},
      {"basis_point", json::object{
           {"exact", exact_match.receiving_local.str()}}},
      {"incoming_point", json::object{
           {"exact", exact_match.producing_local.str()}}},
      {"epsilon", request.at("epsilon")},
      {"basis_checkpoint_identities", std::move(basis_checkpoints)},
      {"incoming_checkpoint_identity", incoming->checkpoint_identity()},
      {"checkpoint_identity", result_checkpoint}};

  const auto producing_rim = exact_plan_rim(
      producing.prescriptions, producing.geometry.scale);
  const auto receiving_rim = exact_plan_rim(
      receiving.prescriptions, receiving.geometry.scale);
  std::shared_ptr<StoredMatchBase> native_match;
  std::optional<json::object> expected_singular_request;
  if (domain == "rational") {
    native_match = build_exact_regular_match(
        match_handle, kernel_request, basis_handles, basis, incoming_handle,
        incoming);
  } else if (domain == "acb") {
    kernel_request["basis_imaginary_sign"] =
        optional_plan_rim_json(receiving_rim);
    kernel_request["incoming_imaginary_sign"] =
        optional_plan_rim_json(producing_rim);
    kernel_request["refinement"] = request.at("refinement");
    if (const auto* exact = request.if_contains("exact_lattice"))
      kernel_request["exact_lattice"] = *exact;
    else if (const auto* native =
                 request.if_contains("native_unit_saturation"))
      kernel_request["native_unit_saturation"] = *native;
    else if (const auto* singular =
                 request.if_contains("native_singular_scc_saturation")) {
      const auto expected = native_acb_saturation_binding(
          plan, active_session_configuration_identity, arm_name,
          match_index, result_checkpoint);
      if (expected.request_key != "native_singular_scc_saturation" ||
          json::serialize(canonical_json_value(*singular)) !=
              json::serialize(canonical_json_value(expected.request)))
        throw std::invalid_argument(
            "planned singular-SCC Acb saturation request does not match the retained receiving SCC");
      kernel_request["native_singular_scc_saturation"] = *singular;
      expected_singular_request = std::move(expected.request);
    } else
      throw std::invalid_argument(
          "planned Acb matching requires an exact lattice, ordinary native unit-leading request, or singular-SCC valuation-zero request");
    native_match = build_refined_acb_match(
        match_handle, kernel_request, basis_handles, basis, incoming_handle,
        incoming, precision_bits, active_session_configuration_identity,
        expected_singular_request);
  } else {
    throw std::invalid_argument(
        "plan-driven local matching requires rational or Acb coefficients");
  }

  const auto native_summary = native_match->summary();
  if (required_string(native_summary, "physical_match_point_exact") !=
      exact_match.physical.str())
    throw std::logic_error(
        "native local match physical point differs from its retained exact plan");

  auto handoff = planned_match_handoff_record(
      plan, arm_name, match_index, basis_handles, basis, incoming_handle,
      incoming, result_checkpoint, native_summary);
  const auto provenance_identity = json::serialize(
      canonical_json_value(handoff));
  const auto elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredPlannedMatchHop>(
      std::move(native_match), result_checkpoint, provenance_identity,
      std::move(handoff), elapsed_ms, plan, basis, incoming);
}

std::shared_ptr<StoredLineResult> build_planned_line_result(
    const std::string& handle, const json::object& request,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& local) {
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  const auto source_checkpoint = required_string(
      request, "source_checkpoint_identity");
  if (source_checkpoint != local->checkpoint_identity())
    throw std::invalid_argument(
        "planned line source checkpoint identity differs from its retained local");
  const auto expected_plan_checkpoint = required_string(
      request, "tile_plan_checkpoint_identity");
  if (expected_plan_checkpoint != plan->checkpoint_identity())
    throw std::invalid_argument(
        "planned line tile-plan checkpoint identity differs from retained state");
  const auto arm_name = required_string(request, "arm");
  const auto tile_index = static_cast<std::size_t>(
      as_u64(request.at("tile"), "planned line tile index"));
  const auto& arm = plan->arm(arm_name);
  if (tile_index >= arm.exact.tiles.size())
    throw std::invalid_argument("planned line tile index is out of range");
  const auto& tile = arm.exact.tiles[tile_index];
  const auto& binding = arm.charts.at(tile.chart);
  if (local->source_chart() != binding.handle)
    throw std::invalid_argument(
        "planned line local does not belong to the tile's retained chart");
  const auto& epsilon = as_object(request.at("epsilon"),
                                  "planned line epsilon window");
  const auto delivered = parse_epsilon_window(
      epsilon, "planned line epsilon window");
  const bool certify_tail =
      request.if_contains("certify_tail") != nullptr &&
      request.at("certify_tail").as_bool();
  auto interval = encode_plan_tile(arm, tile_index);
  const auto rim = exact_plan_rim(
      binding.prescriptions, binding.geometry.scale);
  const auto started = std::chrono::steady_clock::now();
  auto result = local->integrate_planned_line(
      binding.geometry, binding.prescriptions, tile.local_begin,
      tile.local_end, delivered, rim, certify_tail);
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  json::object provenance{
      {"schema",
       "diffexp2-retained-native-physical-tile-integral-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", expected_plan_checkpoint},
      {"arm", arm_name}, {"tile", tile_index},
      {"interval", interval},
      {"source", json::object{
           {"local", local->handle()}, {"chart", local->source_chart()},
           {"checkpoint_identity", source_checkpoint}}},
      {"epsilon", json::object{{"min", delivered.min_power},
                                {"max", delivered.complete_max}}},
      {"tail_certificate_requested", certify_tail},
      {"tail_certificate_status",
       result.diagnostics.tail_certificate_status},
      {"tail_witness_radius_exact",
       result.diagnostics.tail_witness_radius_exact.empty()
           ? json::value(nullptr)
           : json::value(result.diagnostics.tail_witness_radius_exact)},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)},
      {"error_provenance", result.value.error.provenance}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  return std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, arm_name, tile_index,
      std::move(interval), source_checkpoint, std::move(result), elapsed,
      plan, local);
}

std::size_t checked_diagnostic_sum(std::size_t left, std::size_t right,
                                   const char* label) {
  if (right > std::numeric_limits<std::size_t>::max() - left)
    throw std::overflow_error(std::string(label) + " overflow");
  return left + right;
}

StoredLineIntegral aggregate_retained_lines(
    const std::vector<std::shared_ptr<StoredLineResult>>& components,
    const std::vector<std::int32_t>& signs, const std::string& detail) {
  if (components.empty() || components.size() != signs.size())
    throw std::invalid_argument(
        "native line aggregation requires a nonempty signed component list");
  for (const auto sign : signs)
    if (sign != -1 && sign != 1)
      throw std::invalid_argument(
          "native line aggregation signs must be +1 or -1");

  auto epsilon_min = components.front()->result().value.epsilon.min_power;
  auto epsilon_max =
      components.front()->result().value.epsilon.complete_max;
  const auto dimension = components.front()->result().value.dimension;
  for (const auto& component : components) {
    const auto& value = component->result().value;
    if (value.dimension != dimension)
      throw std::invalid_argument(
          "native line aggregate component dimensions differ");
    // A finite Laurent frame's lower edge is an exact structural bound:
    // powers below it are zero, not unknown.  Aggregation therefore takes
    // the union lower edge while still intersecting complete upper edges.
    epsilon_min = std::min(epsilon_min, value.epsilon.min_power);
    epsilon_max = std::min(epsilon_max, value.epsilon.complete_max);
  }
  if (epsilon_min > epsilon_max)
    throw std::domain_error(
        "native line aggregate components have no common complete epsilon window");

  StoredLineIntegral result;
  result.value.epsilon = {epsilon_min, epsilon_max};
  result.value.dimension = dimension;
  result.value.coefficients.reserve(
      result.value.epsilon.width() * dimension);
  for (std::int64_t raw_power = epsilon_min; raw_power <= epsilon_max;
       ++raw_power) {
    const auto power = static_cast<std::int32_t>(raw_power);
    for (std::uint32_t component_index = 0;
         component_index < dimension; ++component_index) {
      ComplexBall sum(0);
      for (std::size_t index = 0; index < components.size(); ++index) {
        const auto& value = components[index]->result().value;
        if (power < value.epsilon.min_power) continue;
        const auto& coefficient = value.at(power, component_index);
        sum += signs[index] == 1 ? coefficient : -coefficient;
      }
      result.value.coefficients.push_back(std::move(sum));
    }
  }

  bool all_error_envelopes = true;
  bool all_certified_errors = true;
  bool any_advisory_error = false;
  json::array error_sources;
  for (std::size_t index = 0; index < components.size(); ++index) {
    const auto& error = components[index]->result().value.error;
    if (error.empty() || error.frame.complete_max < epsilon_max) {
      all_error_envelopes = false;
      break;
    }
    all_certified_errors &= error.guarantee == ErrorGuarantee::Certified;
    any_advisory_error |= error.guarantee == ErrorGuarantee::Advisory;
    error_sources.push_back(json::object{
        {"sign", signs[index]},
        {"guarantee", error_guarantee_name(error.guarantee)},
        {"provenance", error.provenance}});
  }
  if (all_error_envelopes) {
    result.value.error.frame = result.value.epsilon;
    result.value.error.guarantee = all_certified_errors
        ? ErrorGuarantee::Certified
        : any_advisory_error ? ErrorGuarantee::Advisory
                             : ErrorGuarantee::None;
    result.value.error.absolute.reserve(result.value.epsilon.width());
    for (std::int64_t raw_power = epsilon_min; raw_power <= epsilon_max;
         ++raw_power) {
      const auto power = static_cast<std::int32_t>(raw_power);
      auto error_sum = Magnitude::zero();
      for (const auto& component : components) {
        const auto& error = component->result().value.error;
        if (power < error.frame.min_power) continue;
        error_sum += error.absolute.at(static_cast<std::size_t>(
            power - error.frame.min_power));
      }
      result.value.error.absolute.push_back(std::move(error_sum));
    }
    result.value.error.provenance = json::serialize(json::object{
        {"schema", "diffexp2-native-line-error-sum-v1"},
        {"components", std::move(error_sources)}});
  }

  bool all_full_local = true;
  auto& diagnostics = result.diagnostics;
  diagnostics.detail = detail;
  diagnostics.tail_certificate_requested = true;
  diagnostics.tail_certificate_status = "aggregate-certified";
  for (const auto& component : components) {
    const auto& input = component->result();
    all_full_local &=
        input.scope == LineIntegrationScope::FullLocalWithCertifiedTail;
    const auto& source = input.diagnostics;
    diagnostics.input_monomial_cells = checked_diagnostic_sum(
        diagnostics.input_monomial_cells, source.input_monomial_cells,
        "aggregate input monomial count");
    diagnostics.grouped_monomials = checked_diagnostic_sum(
        diagnostics.grouped_monomials, source.grouped_monomials,
        "aggregate grouped monomial count");
    diagnostics.zero_groups_skipped = checked_diagnostic_sum(
        diagnostics.zero_groups_skipped, source.zero_groups_skipped,
        "aggregate skipped-zero count");
    diagnostics.cancelled_divergent_groups = checked_diagnostic_sum(
        diagnostics.cancelled_divergent_groups,
        source.cancelled_divergent_groups,
        "aggregate cancelled-divergence count");
    diagnostics.primitive_evaluations = checked_diagnostic_sum(
        diagnostics.primitive_evaluations, source.primitive_evaluations,
        "aggregate primitive evaluation count");
    diagnostics.primitive_component_applications = checked_diagnostic_sum(
        diagnostics.primitive_component_applications,
        source.primitive_component_applications,
        "aggregate primitive application count");
    diagnostics.primitive_component_reuses = checked_diagnostic_sum(
        diagnostics.primitive_component_reuses,
        source.primitive_component_reuses,
        "aggregate primitive reuse count");
    diagnostics.has_center_endpoint |= source.has_center_endpoint;
    diagnostics.tail_certificate_requested &=
        source.tail_certificate_requested;
  }
  if (all_full_local && result.value.error.guarantee ==
                            ErrorGuarantee::Certified) {
    result.scope = LineIntegrationScope::FullLocalWithCertifiedTail;
  } else {
    result.scope = LineIntegrationScope::StoredTruncation;
    diagnostics.tail_certificate_status = all_full_local
        ? "aggregate-missing-certified-error-envelope"
        : "aggregate-has-stored-truncation-component";
  }
  result.imaginary_sign = std::nullopt;
  return result;
}

std::vector<std::shared_ptr<StoredLocalBase>> unique_line_local_owners(
    const std::vector<std::shared_ptr<StoredLocalBase>>& owners) {
  std::set<std::string> seen;
  std::vector<std::shared_ptr<StoredLocalBase>> unique;
  unique.reserve(owners.size());
  for (const auto& owner : owners) {
    if (!owner)
      throw std::logic_error("native line aggregate lost a local owner");
    if (seen.insert(owner->handle()).second) unique.push_back(owner);
  }
  return unique;
}

struct RetainedLocalFrameContract {
  EpsilonWindow epsilon;
  std::int32_t top_valid = kCompleteInfinity;
  std::uint32_t dimension = 0;
  std::uint32_t taylor_complete_max = 0;
};

struct WholeArmEpsilonContract {
  EpsilonWindow work;
  std::int32_t public_required_complete_max = 0;
  std::int32_t match_required_complete_max = 0;
};

WholeArmEpsilonContract parse_whole_arm_epsilon_contract(
    const json::value& raw, const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(
      object,
      {"min", "max", "required_complete_max",
       "match_required_complete_max"},
      label);
  WholeArmEpsilonContract result{
      {as_i32(object.at("min"), "whole-arm epsilon minimum"),
       as_i32(object.at("max"), "whole-arm epsilon maximum")},
      as_i32(object.at("required_complete_max"),
             "whole-arm projected/public required epsilon maximum"),
      as_i32(object.at("match_required_complete_max"),
             "whole-arm source/match required epsilon maximum")};
  (void)result.work.width();
  if (result.public_required_complete_max < result.work.min_power ||
      result.match_required_complete_max <
          result.public_required_complete_max ||
      result.match_required_complete_max > result.work.complete_max)
    throw std::invalid_argument(
        "whole-arm public and match required epsilon maxima must lie in the work window with match_required_complete_max >= required_complete_max");
  return result;
}

RetainedLocalFrameContract retained_local_frame_contract(
    const std::shared_ptr<StoredLocalBase>& local) {
  if (!local)
    throw std::invalid_argument(
        "native arm frame intersection received a null local");
  const auto summary = local->summary();
  RetainedLocalFrameContract result{
      {as_i32(summary.at("epsilon_min"), "retained epsilon minimum"),
       as_i32(summary.at("epsilon_max"), "retained epsilon maximum")},
      parse_validity(summary.at("top_valid")),
      as_u32(summary.at("dimension"), "retained local dimension"),
      as_u32(summary.at("taylor_complete_max"),
             "retained local Taylor maximum")};
  (void)result.epsilon.width();
  return result;
}

EpsilonWindow live_match_epsilon_intersection(
    EpsilonWindow requested, std::int32_t required_complete_max,
    const std::shared_ptr<StoredLocalBase>& incoming,
    const std::vector<std::shared_ptr<StoredLocalBase>>& basis) {
  if (basis.empty())
    throw std::invalid_argument(
        "native whole-arm match basis cannot be empty");
  const auto incoming_frame = retained_local_frame_contract(incoming);
  auto union_minimum = incoming_frame.epsilon.min_power;
  auto complete_max = requested.complete_max;
  const auto dimension = incoming_frame.dimension;
  if (basis.size() != dimension)
    throw std::invalid_argument(
        "native whole-arm match requires one receiving column per component");
  const auto admit = [&](const RetainedLocalFrameContract& frame) {
    if (frame.dimension != dimension)
      throw std::invalid_argument(
          "native whole-arm matching local dimensions differ");
    // A finite Laurent lower edge is structural: coefficients below it are
    // exact zero.  Clip the caller's lower edge only to the union of actual
    // local frames, then let each matcher zero-pad locals which begin later.
    // Complete upper edges still intersect.
    union_minimum = std::min(union_minimum, frame.epsilon.min_power);
    complete_max = std::min(complete_max, frame.epsilon.complete_max);
    if (frame.top_valid != kCompleteInfinity)
      complete_max = std::min(complete_max, frame.top_valid);
  };
  admit(incoming_frame);
  for (const auto& column : basis)
    admit(retained_local_frame_contract(column));
  const auto minimum = std::max(requested.min_power, union_minimum);
  if (minimum > complete_max)
    throw std::domain_error(
        "native whole-arm match has no common complete epsilon window");
  if (required_complete_max < minimum ||
      required_complete_max > complete_max)
    throw std::domain_error(
        "native whole-arm live match intersection does not cover the globally required complete epsilon maximum");
  return {minimum, complete_max};
}

EpsilonWindow live_line_epsilon_intersection(
    EpsilonWindow requested, std::int32_t required_complete_max,
    const std::shared_ptr<StoredLocalBase>& local) {
  const auto frame = retained_local_frame_contract(local);
  auto minimum = std::max(requested.min_power, frame.epsilon.min_power);
  auto complete_max = std::min(requested.complete_max,
                               frame.epsilon.complete_max);
  if (frame.top_valid != kCompleteInfinity)
    complete_max = std::min(complete_max, frame.top_valid);
  if (minimum > complete_max || required_complete_max > complete_max)
    throw std::domain_error(
        "native whole-arm integrand row does not cover the globally required complete epsilon maximum");
  return {minimum, complete_max};
}

json::object native_unit_saturation_request(
    const std::shared_ptr<StoredTilePlan>& plan, const std::string& arm,
    std::size_t match_index) {
  if (!plan)
    throw std::invalid_argument(
        "native unit-saturation request requires its retained tile plan");
  return json::object{
      {"schema", kNativeUnitSaturationRequestSchema},
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
      {"tile_plan_provenance_identity", plan->provenance_identity()},
      {"arm", arm}, {"match", match_index}};
}

NativeAcbSaturationBinding native_acb_saturation_binding(
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::string& session_configuration_identity,
    const std::string& arm,
    std::size_t match_index, const std::string& match_checkpoint_identity) {
  if (!plan || session_configuration_identity.empty() ||
      match_checkpoint_identity.empty())
    throw std::invalid_argument(
        "native Acb saturation selection requires its session, retained plan, and match checkpoint");
  const auto& retained = plan->arm(arm);
  if (match_index >= retained.exact.matches.size())
    throw std::invalid_argument(
        "native Acb saturation selection match index is outside its retained arm");
  const auto& exact_match = retained.exact.matches[match_index];
  const auto& receiving = retained.charts.at(exact_match.receiving_chart);
  const auto* scc_owner = std::get_if<
      std::shared_ptr<CompositeSCCChartBase>>(&receiving.owner);
  if (scc_owner == nullptr || *scc_owner == nullptr ||
      !is_supported_acb_singular_scc_column_capability(
          (*scc_owner)->column_execution_capability()))
    return {"native_unit_saturation",
            native_unit_saturation_request(plan, arm, match_index)};
  const auto receiving_rim = exact_plan_rim(
      receiving.prescriptions, receiving.geometry.scale);
  return {
      "native_singular_scc_saturation",
      json::object{
          {"schema", kNativeSingularSCCSaturationRequestSchema},
          {"session_configuration_identity",
           session_configuration_identity},
          {"tile_plan", plan->handle()},
          {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
          {"tile_plan_provenance_identity", plan->provenance_identity()},
          {"arm", arm},
          {"match", match_index},
          {"match_checkpoint_identity", match_checkpoint_identity},
          {"receiving_scc", (*scc_owner)->handle()},
          {"receiving_scc_exact_identity", (*scc_owner)->exact_identity()},
          {"receiving_execution_capability",
           (*scc_owner)->column_execution_capability()},
          {"receiving_basis_point_exact",
           exact_match.receiving_local.str()},
          {"physical_match_point_exact", exact_match.physical.str()},
          {"receiving_rim", optional_plan_rim_json(receiving_rim)}}};
}

struct RetainedArmMarchInput {
  std::string name;
  std::vector<std::vector<std::string>> basis_handles;
  std::vector<std::vector<std::shared_ptr<StoredLocalBase>>> basis;
  std::vector<std::string> match_handles;
  std::vector<std::string> local_handles;
};

struct RetainedArmMarchResult {
  std::vector<std::shared_ptr<StoredPlannedMatchHop>> matches;
  std::vector<std::shared_ptr<StoredLocalBase>> tile_sources;
  double elapsed_ms = 0.0;

  const std::shared_ptr<StoredLocalBase>& final_local() const {
    if (tile_sources.empty())
      throw std::logic_error(
          "retained arm march has no tile source");
    return tile_sources.back();
  }
};

std::string arm_checkpoint_identity(const std::string& root,
                                    const std::string& arm,
                                    const char* kind,
                                    std::size_t one_based_index) {
  return root + ":" + arm + ":" + kind + ":" +
         std::to_string(one_based_index);
}

RetainedArmMarchResult march_retained_arm(
    const std::string& domain, slong precision_bits,
    const std::string& session_configuration_identity,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& anchor,
    const RetainedArmMarchInput& input, EpsilonWindow work_epsilon,
    std::int32_t match_required_complete_max,
    const json::object& refinement, const std::string& checkpoint_root) {
  if ((domain != "rational" && domain != "acb") || !plan || !anchor ||
      session_configuration_identity.empty() || checkpoint_root.empty())
    throw std::invalid_argument(
        "retained arm march requires a numeric domain, session identity, plan, anchor, and checkpoint root");
  const auto& retained = plan->arm(input.name);
  const auto match_count = retained.exact.matches.size();
  if (input.basis_handles.size() != match_count ||
      input.basis.size() != match_count ||
      input.match_handles.size() != match_count ||
      input.local_handles.size() != match_count)
    throw std::invalid_argument(
        "retained arm march input does not reproduce its plan match count");
  require_exact_keys(refinement, {"relative_tolerance", "max_steps"},
                     "retained arm march refinement policy");
  if (required_string(refinement, "relative_tolerance").empty() ||
      as_u32(refinement.at("max_steps"),
             "retained arm march refinement steps") > 32)
    throw std::invalid_argument(
        "retained arm march refinement policy is invalid");
  (void)work_epsilon.width();
  if (domain == "acb") ComplexBall::set_precision(precision_bits);

  const auto started = std::chrono::steady_clock::now();
  RetainedArmMarchResult output;
  output.matches.reserve(match_count);
  output.tile_sources.reserve(retained.exact.tiles.size());
  output.tile_sources.push_back(anchor);
  std::shared_ptr<StoredLocalBase> current = anchor;
  for (std::size_t match_index = 0; match_index < match_count;
       ++match_index) {
    const auto match_checkpoint = arm_checkpoint_identity(
        checkpoint_root, input.name, "match", match_index + 1);
    const auto match_epsilon = live_match_epsilon_intersection(
        work_epsilon, match_required_complete_max, current,
        input.basis[match_index]);
    json::object match_request{
        {"arm", input.name}, {"match", match_index},
        {"epsilon", json::object{
             {"min", match_epsilon.min_power},
             {"max", match_epsilon.complete_max},
             {"required_complete_max", match_required_complete_max}}},
        {"checkpoint_identity", match_checkpoint}};
    if (domain == "acb") {
      auto saturation = native_acb_saturation_binding(
          plan, session_configuration_identity, input.name, match_index,
          match_checkpoint);
      match_request[saturation.request_key] =
          std::move(saturation.request);
      match_request["refinement"] = refinement;
    }
    auto match = build_planned_match_hop(
        input.match_handles[match_index], match_request, domain,
        precision_bits, session_configuration_identity, plan,
        input.basis_handles[match_index], input.basis[match_index],
        current->handle(), current);
    auto next = match->materialize(
        input.local_handles[match_index],
        arm_checkpoint_identity(checkpoint_root, input.name, "local",
                                match_index + 1),
        precision_bits, match);
    output.matches.push_back(std::move(match));
    output.tile_sources.push_back(next);
    current = std::move(next);
  }
  if (output.tile_sources.size() != retained.exact.tiles.size())
    throw std::logic_error(
        "retained arm march did not produce one source local per tile");
  output.elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return output;
}

json::array line_aggregate_source_records(
    const std::vector<std::shared_ptr<StoredLocalBase>>& owners) {
  json::array records;
  records.reserve(owners.size());
  for (const auto& owner : owners)
    records.push_back(json::object{
        {"local", owner->handle()}, {"chart", owner->source_chart()},
        {"source_operator_identity", owner->source_operator_identity()},
        {"checkpoint_identity", owner->checkpoint_identity()},
        {"coefficient_domain", owner->scalar_domain()},
        {"analytic_metadata", owner->exact_analytic_metadata()},
        {"retained_derivation", owner->retained_derivation().has_value()
             ? json::value(*owner->retained_derivation())
             : json::value(nullptr)}});
  return records;
}

std::shared_ptr<StoredLineResult> build_retained_line_aggregate(
    const std::string& handle, const std::string& checkpoint_identity,
    const std::string& arm_name, json::object interval,
    json::object aggregate_record,
    const std::shared_ptr<StoredTilePlan>& plan,
    std::vector<std::shared_ptr<StoredLocalBase>> local_owners,
    const std::vector<std::shared_ptr<StoredLineResult>>& components,
    const std::vector<std::int32_t>& signs, double elapsed_ms) {
  if (checkpoint_identity.empty())
    throw std::invalid_argument(
        "native line aggregate checkpoint identity cannot be empty");
  local_owners = unique_line_local_owners(local_owners);
  auto result = aggregate_retained_lines(
      components, signs, "native retained line aggregate");
  json::array component_records;
  component_records.reserve(components.size());
  for (std::size_t index = 0; index < components.size(); ++index) {
    const auto summary = components[index]->summary();
    component_records.push_back(json::object{
        {"index", index}, {"sign", signs[index]},
        {"checkpoint_identity", components[index]->checkpoint_identity()},
        {"provenance_identity", components[index]->provenance_identity()},
        {"scope", summary.at("scope")},
        {"source", summary.at("source")},
        {"interval", summary.at("interval")},
        {"epsilon", json::object{
             {"min", components[index]->result().value.epsilon.min_power},
             {"max", components[index]->result().value.epsilon.complete_max}}},
        {"error", summary.at("error")}});
  }
  aggregate_record["component_count"] = components.size();
  aggregate_record["components"] = std::move(component_records);
  json::object provenance{
      {"schema", "diffexp2-retained-native-line-aggregate-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"arm", arm_name},
      {"interval", std::move(interval)},
      {"source", json::object{
           {"tile_plan", plan->handle()},
           {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
           {"locals", line_aggregate_source_records(local_owners)}}},
      {"aggregate", std::move(aggregate_record)},
      {"epsilon", json::object{{"min", result.value.epsilon.min_power},
                                {"max", result.value.epsilon.complete_max}}},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  return std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, std::move(result),
      elapsed_ms, plan, std::move(local_owners), std::move(provenance));
}

std::vector<std::uint32_t> parse_index_vector(
    const json::value& raw, std::uint32_t bound, const char* label) {
  std::vector<std::uint32_t> values;
  for (const auto& value : as_array(raw, label)) {
    const auto index = as_u32(value, label);
    if (index >= bound)
      throw std::invalid_argument(std::string(label) + " is outside range");
    values.push_back(index);
  }
  return values;
}

struct ExactParentMatrix {
  struct Cell {
    std::string exact;
    bool proven_zero = false;
  };
  std::uint32_t dimension = 0;
  std::vector<Cell> cells;  // row (target), then column (source)
  std::string canonical_record;

  const Cell& at(std::uint32_t row, std::uint32_t column) const {
    return cells.at(static_cast<std::size_t>(row) * dimension + column);
  }
};

ExactParentMatrix parse_exact_parent_matrix(
    const json::value& raw, std::uint32_t dimension, const char* label) {
  const auto& rows = as_array(raw, label);
  if (rows.size() != dimension)
    throw std::invalid_argument(
        std::string(label) + " must have parent-dimension rows");
  ExactParentMatrix result;
  result.dimension = dimension;
  result.cells.reserve(static_cast<std::size_t>(dimension) * dimension);
  json::array canonical_rows;
  canonical_rows.reserve(dimension);
  for (std::uint32_t row = 0; row < dimension; ++row) {
    const auto& columns = as_array(rows[row], label);
    if (columns.size() != dimension)
      throw std::invalid_argument(
          std::string(label) + " must be a square parent-dimension matrix");
    json::array canonical_columns;
    canonical_columns.reserve(dimension);
    for (std::uint32_t column = 0; column < dimension; ++column) {
      const auto& cell = as_object(columns[column], "exact parent matrix cell");
      const auto exact = required_string(cell, "exact");
      const auto proven_zero = cell.at("proven_zero").as_bool();
      result.cells.push_back({exact, proven_zero});
      canonical_columns.push_back(json::object{
          {"exact", exact}, {"proven_zero", proven_zero}});
    }
    canonical_rows.push_back(std::move(canonical_columns));
  }
  result.canonical_record = json::serialize(canonical_rows);
  return result;
}

std::string canonical_exact_submatrix(
    const ExactParentMatrix& matrix,
    const std::vector<std::uint32_t>& vertices) {
  json::array rows;
  rows.reserve(vertices.size());
  for (const auto target : vertices) {
    json::array columns;
    columns.reserve(vertices.size());
    for (const auto source : vertices) {
      const auto& cell = matrix.at(target, source);
      columns.push_back(json::object{
          {"exact", cell.exact}, {"proven_zero", cell.proven_zero}});
    }
    rows.push_back(std::move(columns));
  }
  return json::serialize(rows);
}

std::string derived_coupling_identity(
    std::uint32_t source_block, std::uint32_t target_block,
    const std::vector<std::uint32_t>& source_vertices,
    const std::vector<std::uint32_t>& target_vertices,
    const ExactParentMatrix& original, const ExactParentMatrix& theta) {
  auto rectangular_record = [&](const ExactParentMatrix& matrix) {
    json::array rows;
    rows.reserve(target_vertices.size());
    for (const auto target : target_vertices) {
      json::array columns;
      columns.reserve(source_vertices.size());
      for (const auto source : source_vertices) {
        const auto& cell = matrix.at(target, source);
        columns.push_back(json::object{
            {"exact", cell.exact}, {"proven_zero", cell.proven_zero}});
      }
      rows.push_back(std::move(columns));
    }
    return rows;
  };
  return json::serialize(json::object{
      {"schema", "diffexp2-native-scc-coupling-v1"},
      {"source_block", source_block}, {"target_block", target_block},
      {"source_vertices", encode_indices(source_vertices)},
      {"target_vertices", encode_indices(target_vertices)},
      {"exact_original", rectangular_record(original)},
      {"exact_theta", rectangular_record(theta)}});
}

template <typename Scalar>
std::shared_ptr<CompositeSCCChartBase> parse_composite_scc_chart(
    const std::shared_ptr<SolverSession>& session, const json::object& root,
    const std::string& handle, const std::string& key,
    const std::string& exact_identity, std::string signature,
    const std::vector<std::shared_ptr<PreparedChartBase>>& erased_charts) {
  // Composite geometry is retained as an exact-rational Acb ball even when
  // recurrence coefficients themselves are exact.  Parse it under the same
  // guarded precision later used by block LocalSolutions.
  AcbPrecisionLease acb_lease(session->precision_bits);
  ComplexBall::set_precision(session->precision_bits);
  std::unique_lock<std::recursive_mutex> symbolic_lock;
  if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
    symbolic_lock =
        std::unique_lock<std::recursive_mutex>(symbolic_run_mutex());
    SymbolicRational::configure(session->symbols);
  }

  const auto& parent = as_object(root.at("parent"), "SCC parent manifest");
  const auto dimension = as_u32(parent.at("dimension"), "parent dimension");
  if (dimension == 0)
    throw std::invalid_argument("SCC parent dimension must be positive");
  const auto exact_system = parse_exact_parent_matrix(
      parent.at("exact_system_record"), dimension,
      "exact parent system record");
  const auto exact_theta = parse_exact_parent_matrix(
      parent.at("exact_theta_record"), dimension,
      "exact parent theta record");
  validate_first_slice_rational_geometry(parent.at("chart"));
  const auto geometry_record = canonical_chart_geometry_record(
      parent.at("chart"));
  auto retained_geometry = parse_retained_composite_geometry(
      parent.at("chart"));
  auto graph = validate_scc_certificate(parent.at("scc"), dimension);
  const std::set<std::pair<std::uint32_t, std::uint32_t>> structural_edges(
      graph.structural_edges.begin(), graph.structural_edges.end());
  for (std::uint32_t target = 0; target < dimension; ++target) {
    for (std::uint32_t source = 0; source < dimension; ++source) {
      const auto original_zero = exact_system.at(target, source).proven_zero;
      const auto theta_zero = exact_theta.at(target, source).proven_zero;
      if (original_zero != theta_zero)
        throw std::invalid_argument(
            "parent original/theta structural-zero facts disagree");
      const auto graph_zero = !structural_edges.contains({source, target});
      if (original_zero != graph_zero)
        throw std::invalid_argument(
            "parent exact matrix zero facts do not reproduce the structural graph");
    }
  }

  const auto& execution = as_object(
      parent.at("execution"), "SCC execution contract");
  if (required_string(execution, "mode") != "BlockSequentialStrict")
    throw std::invalid_argument(
        "native SCC preparation only supports BlockSequentialStrict");
  CompositeWorkContract work;
  work.work_t_order = as_u32(
      execution.at("work_t_order"), "SCC work Taylor order");
  const auto& work_object = as_object(
      parent.at("work_contract"), "SCC work contract");
  work.work_min = as_i32(
      work_object.at("work_min"), "work epsilon minimum");
  work.requested_min = as_i32(
      work_object.at("requested_min"), "requested epsilon minimum");
  work.requested_max = as_i32(
      work_object.at("requested_max"), "requested epsilon maximum");
  work.work_complete_max = as_i32(
      work_object.at("work_complete_max"), "work epsilon maximum");
  work.public_t_order = as_u32(
      work_object.at("public_t_order"), "public Taylor order");
  work.wolfram_coupling_depth = as_u32(
      work_object.at("wolfram_coupling_depth"),
      "Wolfram coupling depth");
  if (work.work_min > work.requested_min ||
      work.requested_min > work.requested_max ||
      work.requested_max > work.work_complete_max)
    throw std::invalid_argument(
        "SCC epsilon work contract has inconsistent ordered bounds");
  if (work.wolfram_coupling_depth != graph.coupling_depth + 1)
    throw std::invalid_argument(
        "Wolfram coupling depth must equal native edge depth plus one");
  const auto expected_work_t_order =
      static_cast<std::uint64_t>(work.public_t_order) + 2 +
      2 * static_cast<std::uint64_t>(work.wolfram_coupling_depth);
  if (expected_work_t_order > std::numeric_limits<std::uint32_t>::max() ||
      work.work_t_order != expected_work_t_order)
    throw std::invalid_argument(
        "SCC work Taylor order does not match the exact depth budget");
  const auto frame_width_i64 =
      static_cast<std::int64_t>(work.work_complete_max) -
      work.work_min + 1;
  if (frame_width_i64 <= 0 ||
      frame_width_i64 > std::numeric_limits<std::uint32_t>::max())
    throw std::invalid_argument("SCC work frame width is invalid");
  const auto frame_width = static_cast<std::uint32_t>(frame_width_i64);

  const auto& raw_blocks = as_array(root.at("blocks"), "SCC blocks");
  if (raw_blocks.size() != graph.component_count ||
      erased_charts.size() != raw_blocks.size())
    throw std::invalid_argument(
        "SCC preparation requires exactly one chart per component");
  std::vector<CompositeSCCBlock<Scalar>> blocks(graph.component_count);
  std::vector<std::uint8_t> block_seen(graph.component_count, 0);
  std::set<std::string> chart_handles;
  for (std::size_t index = 0; index < raw_blocks.size(); ++index) {
    const auto& raw_block = as_object(raw_blocks[index], "SCC block");
    const auto block = as_u32(raw_block.at("block"), "SCC block index");
    if (block != index || block >= graph.component_count || block_seen[block])
      throw std::invalid_argument(
          "SCC blocks must be in deterministic component order");
    block_seen[block] = 1;
    const auto declared_capabilities =
        canonical_native_scc_capabilities(raw_block);
    const bool regular = raw_block.at("regular").as_bool();
    const bool no_pseudo = raw_block.at("no_pseudo").as_bool();
    if (!raw_block.at("identity_gauge").as_bool())
      throw std::invalid_argument(
          "native SCC preparation requires an exact identity gauge");
    if (!raw_block.at("identity_v").as_bool())
      throw std::invalid_argument(
          "native SCC preparation requires an exact identity spectral transform");
    // `no_pseudo` is retained as producer provenance, not an admission
    // decision.  Exact Rational execution reconstructs every T/P/R branch
    // from the retained affine Jordan certificate; a producer may therefore
    // truthfully advertise false here without disabling native execution.
    // Rational execution may compensate a revalidated CASE-P event.  Acb
    // execution below admits only the producer-proved no-collision class and
    // still reconstructs every submitted schedule from an exact certificate.

    auto vertices = parse_index_vector(
        raw_block.at("vertices"), dimension, "SCC block vertex");
    if (vertices != graph.components[block])
      throw std::invalid_argument(
          "SCC block vertices do not equal their parent component in exact order");

    // This is the only type-erasure boundary for a composite.  Every later
    // execution method owns typed pointers directly and performs no casts.
    auto chart = std::dynamic_pointer_cast<PreparedChart<Scalar>>(
        erased_charts[index]);
    if (!chart)
      throw std::invalid_argument(
          "SCC block chart scalar domain differs from its session");
    const auto source_handle = required_string(raw_block, "chart");
    if (chart->handle() != source_handle ||
        !chart_handles.insert(source_handle).second)
      throw std::invalid_argument(
          "SCC block chart handles must be exact and one-to-one");
    const auto principal_identity = required_string(
        raw_block, "principal_identity");
    if (chart->exact_identity() != principal_identity)
      throw std::invalid_argument(
          "SCC principal identity differs from its prepared chart identity");
    if (!chart->native_scc_capabilities().has_value())
      throw std::invalid_argument(
          "SCC block chart lacks analytic.native_scc_capabilities metadata");
    if (*chart->native_scc_capabilities() != declared_capabilities)
      throw std::invalid_argument(
          "SCC block capability claims differ from the retained chart metadata");
    if (!chart->has_identity_assembly())
      throw std::invalid_argument(
          "SCC identity_v capability contradicts the retained assembly operator");
    if (!chart->geometry_record().has_value())
      throw std::invalid_argument(
          "SCC block chart lacks analytic.geometry preparation metadata");
    if (*chart->geometry_record() != geometry_record)
      throw std::invalid_argument(
          "SCC block chart geometry differs from the parent manifest");
    if (!chart->principal_matrix_record().has_value())
      throw std::invalid_argument(
          "SCC block chart lacks analytic.principal_matrix metadata");
    if (*chart->principal_matrix_record() !=
        canonical_exact_submatrix(exact_system, vertices))
      throw std::invalid_argument(
          "SCC principal chart matrix differs from the indexed parent submatrix");
    if (chart->dimension() != vertices.size())
      throw std::invalid_argument(
          "SCC block chart dimension differs from its vertex count");
    if (chart->frame_base() != work.work_min ||
        chart->frame_width() != frame_width)
      throw std::invalid_argument(
          "SCC block chart frame differs from the exact work contract");

    std::optional<ExactJordanIndicialCertificate> exact_indicial;
    if (const auto* raw_indicial =
            raw_block.if_contains("exact_affine_jordan_indicial"))
      exact_indicial = parse_exact_jordan_indicial_record(
          *raw_indicial, chart->dimension());
    if (const auto& retained = chart->exact_jordan_indicial();
        retained.has_value()) {
      if (exact_indicial.has_value() &&
          !same_exact_jordan_indicial(*exact_indicial, *retained))
        throw std::invalid_argument(
            "SCC block exact indicial manifest differs from its retained Rational operator certificate");
      exact_indicial = *retained;
    }
    if (exact_indicial.has_value() &&
        !chart->jordan_partition_matches(*exact_indicial))
      throw std::invalid_argument(
          "SCC block exact indicial certificate differs from its prepared Jordan partition");

    const auto& local_graph = chart->scc();
    std::vector<std::uint32_t> local_vertices(vertices.size());
    for (std::uint32_t local = 0; local < vertices.size(); ++local)
      local_vertices[local] = local;
    if (local_graph.component_count != 1 ||
        local_graph.components.size() != 1 ||
        local_graph.components.front() != local_vertices ||
        local_graph.coupling_depth != 0)
      throw std::invalid_argument(
          "SCC principal chart must be exactly one retained component");
    std::vector<std::uint32_t> local_of(dimension,
        std::numeric_limits<std::uint32_t>::max());
    for (std::uint32_t local = 0; local < vertices.size(); ++local)
      local_of[vertices[local]] = local;
    std::set<std::pair<std::uint32_t, std::uint32_t>> expected_intra;
    for (const auto [source, target] : graph.structural_edges)
      if (graph.component_of[source] == block &&
          graph.component_of[target] == block)
        expected_intra.insert({local_of[source], local_of[target]});
    const std::set<std::pair<std::uint32_t, std::uint32_t>> supplied_intra(
        local_graph.structural_edges.begin(),
        local_graph.structural_edges.end());
    if (supplied_intra != expected_intra)
      throw std::invalid_argument(
          "SCC principal chart graph does not cover its exact intra-component edges");

    blocks[block] = CompositeSCCBlock<Scalar>{
        block, std::move(vertices), source_handle, principal_identity,
        regular, no_pseudo, std::move(exact_indicial), std::move(chart)};
  }

  std::set<std::pair<std::uint32_t, std::uint32_t>> expected_cross;
  for (const auto edge : graph.structural_edges)
    if (graph.component_of[edge.first] != graph.component_of[edge.second])
      expected_cross.insert(edge);
  const std::set<std::pair<std::uint32_t, std::uint32_t>>
      expected_condensation(graph.condensation_edges.begin(),
                            graph.condensation_edges.end());
  std::set<std::pair<std::uint32_t, std::uint32_t>> supplied_condensation;
  std::set<std::pair<std::uint32_t, std::uint32_t>> active_edges;
  std::set<std::pair<std::uint32_t, std::uint32_t>> prepared_pairs;
  std::vector<CompositeSCCCoupling<Scalar>> couplings;
  const auto& raw_couplings = as_array(
      root.at("couplings"), "SCC couplings");
  if (raw_couplings.size() != graph.condensation_edges.size())
    throw std::invalid_argument(
        "SCC coupling groups must equal the condensation edge count");
  couplings.reserve(raw_couplings.size());
  for (std::size_t coupling_index = 0;
       coupling_index < raw_couplings.size(); ++coupling_index) {
    const auto& raw_value = raw_couplings[coupling_index];
    const auto& raw = as_object(raw_value, "SCC coupling");
    const auto source_block = as_u32(
        raw.at("source_block"), "source SCC block");
    const auto target_block = as_u32(
        raw.at("target_block"), "target SCC block");
    if (source_block >= graph.component_count ||
        target_block >= graph.component_count || source_block == target_block)
      throw std::invalid_argument("invalid cross-SCC coupling block pair");
    const auto block_pair = std::make_pair(source_block, target_block);
    if (block_pair != graph.condensation_edges[coupling_index])
      throw std::invalid_argument(
          "SCC coupling groups must use deterministic condensation-edge order");
    if (!supplied_condensation.insert(block_pair).second ||
        !expected_condensation.contains(block_pair))
      throw std::invalid_argument(
          "SCC coupling groups must match condensation edges one-to-one");
    const auto source_vertices = parse_index_vector(
        raw.at("source_vertices"), dimension, "coupling source vertex");
    const auto target_vertices = parse_index_vector(
        raw.at("target_vertices"), dimension, "coupling target vertex");
    if (source_vertices != blocks[source_block].vertices ||
        target_vertices != blocks[target_block].vertices)
      throw std::invalid_argument(
          "SCC coupling vertex bases differ from their block manifests");
    const auto columns = as_u32(raw.at("columns"), "coupling columns");
    const auto rows = as_u32(raw.at("rows"), "coupling rows");
    if (columns != source_vertices.size() || rows != target_vertices.size())
      throw std::invalid_argument(
          "SCC coupling matrix dimensions disagree with its blocks");
    if (raw.if_contains("symbols") == nullptr ||
        required_string(raw, "domain") != session->domain ||
        parse_symbols(raw) != session->symbols)
      throw std::invalid_argument(
          "SCC coupling coefficient field differs from its session");

    CompositeSCCCoupling<Scalar> coupling;
    coupling.source_block = source_block;
    coupling.target_block = target_block;
    coupling.producer_identity = required_string(raw, "exact_identity");
    coupling.matrix.rows = rows;
    coupling.matrix.columns = columns;
    coupling.matrix.exact_identity = derived_coupling_identity(
        source_block, target_block, source_vertices, target_vertices,
        exact_system, exact_theta);
    std::optional<std::pair<std::uint32_t, std::uint32_t>> previous_local;
    for (const auto& raw_entry_value : as_array(
             raw.at("entries"), "SCC coupling entries")) {
      const auto& raw_entry = as_object(
          raw_entry_value, "SCC coupling entry");
      const auto row = as_u32(raw_entry.at("row"), "coupling row");
      const auto column = as_u32(
          raw_entry.at("column"), "coupling column");
      if (row >= rows || column >= columns)
        throw std::invalid_argument("SCC coupling entry is out of range");
      const auto local_pair = std::make_pair(row, column);
      if (previous_local.has_value() && *previous_local >= local_pair)
        throw std::invalid_argument(
            "SCC coupling entries must use deterministic row-column order");
      previous_local = local_pair;
      const auto source_vertex = as_u32(
          raw_entry.at("source_vertex"), "coupling source vertex");
      const auto target_vertex = as_u32(
          raw_entry.at("target_vertex"), "coupling target vertex");
      if (source_vertex != source_vertices[column] ||
          target_vertex != target_vertices[row] ||
          graph.component_of[source_vertex] != source_block ||
          graph.component_of[target_vertex] != target_block)
        throw std::invalid_argument(
            "SCC coupling local/global vertex identities disagree");
      const auto edge = std::make_pair(source_vertex, target_vertex);
      if (!prepared_pairs.insert(edge).second)
        throw std::invalid_argument(
            "SCC coupling contains a duplicate global matrix entry");
      const auto exact_original_entry = required_string(
          raw_entry, "exact_original_entry");
      const auto exact_theta_entry = required_string(
          raw_entry, "exact_theta_entry");
      const auto& parent_original = exact_system.at(
          target_vertex, source_vertex);
      const auto& parent_theta = exact_theta.at(
          target_vertex, source_vertex);
      if (exact_original_entry != parent_original.exact ||
          exact_theta_entry != parent_theta.exact)
        throw std::invalid_argument(
            "SCC coupling exact entries differ from their indexed parent cells");
      const auto& raw_multiplier = as_object(
          raw_entry.at("multiplier"), "prepared SCC multiplier");
      PreparedRationalTaylorMultiplier<Scalar> multiplier;
      multiplier.epsilon_shift = as_i32(
          raw_multiplier.at("epsilon_shift"), "multiplier epsilon shift");
      const auto shifted_work_min = static_cast<std::int64_t>(work.work_min) +
          multiplier.epsilon_shift;
      const auto shifted_work_max =
          static_cast<std::int64_t>(work.work_complete_max) +
          multiplier.epsilon_shift;
      if (shifted_work_min < std::numeric_limits<std::int32_t>::min() ||
          shifted_work_min > std::numeric_limits<std::int32_t>::max() ||
          shifted_work_max < std::numeric_limits<std::int32_t>::min() ||
          shifted_work_max > std::numeric_limits<std::int32_t>::max())
        throw std::invalid_argument(
            "SCC multiplier shift overflows the retained epsilon frame");
      multiplier.center_pole_order = as_u32(
          raw_multiplier.at("center_pole_order"),
          "multiplier center-pole order");
      multiplier.exact_identity = required_string(
          raw_multiplier, "exact_identity");
      if (multiplier.exact_identity != exact_theta_entry)
        throw std::invalid_argument(
            "prepared multiplier identity differs from the exact theta entry");
      multiplier.proven_zero = raw_multiplier.at("proven_zero").as_bool();
      if (multiplier.proven_zero != parent_original.proven_zero)
        throw std::invalid_argument(
            "prepared multiplier structural-zero fact differs from its exact parent cell");
      for (const auto& raw_kernel : as_array(
               raw_multiplier.at("kernels"), "multiplier epsilon kernels")) {
        std::vector<Scalar> kernel;
        for (const auto& coefficient : as_array(
                 raw_kernel, "multiplier Taylor kernel"))
          kernel.push_back(parse_scalar<Scalar>(coefficient));
        multiplier.kernels.push_back(std::move(kernel));
      }
      if (!multiplier.proven_zero) {
        if (multiplier.kernels.size() != frame_width ||
            std::any_of(multiplier.kernels.begin(), multiplier.kernels.end(),
                [&](const auto& kernel) {
                  return kernel.size() !=
                      static_cast<std::size_t>(work.work_t_order) + 1;
                }))
          throw std::invalid_argument(
              "active SCC multiplier kernels do not cover the exact work rectangle");
        if (!expected_cross.contains(edge))
          throw std::invalid_argument(
              "active SCC multiplier lacks an exact parent structural edge");
        active_edges.insert(edge);
      } else if (expected_cross.contains(edge)) {
        throw std::invalid_argument(
            "a proven-zero SCC multiplier contradicts a parent structural edge");
      }
      coupling.identities.push_back(CompositeCouplingIdentity{
          source_vertex, target_vertex, exact_original_entry,
          exact_theta_entry, multiplier.proven_zero});
      coupling.matrix.entries.push_back(
          typename PreparedSparseLocalMultiplierMatrix<Scalar>::Entry{
              row, column, std::move(multiplier)});
    }
    couplings.push_back(std::move(coupling));
  }
  if (supplied_condensation != expected_condensation)
    throw std::invalid_argument(
        "SCC coupling groups do not cover the condensation graph exactly");
  if (active_edges != expected_cross)
    throw std::invalid_argument(
        "active SCC coupling entries do not cover cross-component structural edges exactly");

  return make_retained_typed_shared<Scalar, CompositeSCCChart<Scalar>>(
      handle, key, exact_identity, std::move(signature), dimension,
      std::move(graph), exact_system.canonical_record,
      exact_theta.canonical_record,
      geometry_record, std::move(retained_geometry), work,
      std::move(blocks), std::move(couplings));
}

struct SessionRegistry {
  std::mutex mutex;
  std::uint64_t next_session = 1;
  std::unordered_map<std::string, std::shared_ptr<SolverSession>> sessions;
};

SessionRegistry& session_registry() {
  // LibraryLink's uninitialize hook clears this registry.  Deliberately leak
  // the trivial holder itself so static destruction can never run after the
  // FLINT context used by retained SymbolicRational values.
  static auto* registry = new SessionRegistry;
  return *registry;
}

std::shared_ptr<SolverSession> find_session(const std::string& handle) {
  auto& registry = session_registry();
  std::lock_guard<std::mutex> lock(registry.mutex);
  const auto found = registry.sessions.find(handle);
  if (found == registry.sessions.end())
    throw std::invalid_argument("unknown or closed persistent solver session");
  return found->second;
}

json::object solve_prepared_chart_safe(
    const std::shared_ptr<PreparedChartBase>& chart,
    const json::value& raw_run, int output_digits,
    const std::string& session_handle) {
  try {
    auto result = chart->solve(
        as_object(raw_run, "persistent recurrence run"), output_digits);
    result["session"] = session_handle;
    result["chart"] = chart->handle();
    return result;
  } catch (const RecurrenceError& error) {
    return json::object{{"status", "error"}, {"id", error.id},
                        {"detail", error.what()},
                        {"frame_base", error.frame_base},
                        {"shift", error.shift},
                        {"session", session_handle},
                        {"chart", chart->handle()}};
  } catch (const std::exception& error) {
    return json::object{{"status", "error"}, {"id", "CPP"},
                        {"detail", error.what()},
                        {"session", session_handle},
                        {"chart", chart->handle()}};
  }
}

constexpr std::uint32_t kMaxPersistentBatchThreads = 64;

json::value canonical_json_value(const json::value& value) {
  if (value.is_array()) {
    json::array output;
    output.reserve(value.as_array().size());
    for (const auto& item : value.as_array())
      output.push_back(canonical_json_value(item));
    return output;
  }
  if (value.is_object()) {
    std::vector<std::string> keys;
    keys.reserve(value.as_object().size());
    for (const auto& item : value.as_object())
      keys.emplace_back(item.key());
    std::sort(keys.begin(), keys.end());
    json::object output;
    for (const auto& key : keys)
      output[key] = canonical_json_value(value.as_object().at(key));
    return output;
  }
  return value;
}

std::string composite_scc_signature(const json::object& root) {
  json::object exact{{"identity", root.at("identity")},
                     {"parent", root.at("parent")},
                     {"blocks", root.at("blocks")},
                     {"couplings", root.at("couplings")}};
  return json::serialize(canonical_json_value(exact));
}

std::string static_problem_signature(const json::object& problem,
                                     const json::value& analytic,
                                     const SCCCertificate& scc,
                                     const std::string& identity) {
  json::object exact;
  for (const auto* key : {"domain", "symbols", "precision_bits", "d", "fb",
                          "w", "d_lags", "denominators", "nhat_lags",
                          "d0_inverse", "blocks", "assembly", "chop_digits"}) {
    if (const auto* value = problem.if_contains(key)) exact[key] = *value;
  }
  exact["identity"] = identity;
  exact["analytic"] = analytic;
  exact["scc"] = json::parse(scc.exact_record);
  return json::serialize(exact);
}

constexpr const char* kCheckpointFormat =
    "diffexp2-persistent-native-session";
constexpr std::uint32_t kCheckpointPayloadSchema = 3;

json::object run_session_command(const json::object& root);

void require_exact_keys(const json::object& object,
                        std::initializer_list<std::string_view> expected,
                        const char* label) {
  std::set<std::string> actual;
  for (const auto& entry : object) actual.emplace(entry.key());
  std::set<std::string> wanted;
  for (const auto key : expected) wanted.emplace(key);
  if (actual != wanted)
    throw std::invalid_argument(std::string(label) +
                                " has unknown or missing fields");
}

std::uint64_t scoped_handle_id(const std::string& handle,
                               std::string_view prefix,
                               const char* label) {
  if (!handle.starts_with(prefix) || handle.size() == prefix.size())
    throw std::invalid_argument(std::string(label) +
                                " has an invalid scoped handle");
  std::uint64_t value = 0;
  for (std::size_t index = prefix.size(); index < handle.size(); ++index) {
    const char digit = handle[index];
    if (digit < '0' || digit > '9')
      throw std::invalid_argument(std::string(label) +
                                  " has an invalid scoped handle");
    const auto unsigned_digit = static_cast<std::uint64_t>(digit - '0');
    if (value > (std::numeric_limits<std::uint64_t>::max() -
                 unsigned_digit) / 10)
      throw std::invalid_argument(std::string(label) +
                                  " scoped handle overflows uint64");
    value = value * 10 + unsigned_digit;
  }
  if (value == 0)
    throw std::invalid_argument(std::string(label) +
                                " scoped handle index must be positive");
  return value;
}

std::pair<ExactArmRequest, std::vector<RetainedPlanChartBinding>>
parse_checkpoint_retained_arm(
    const json::value& raw,
    const std::unordered_map<std::string,
                             std::shared_ptr<PreparedChartBase>>& charts,
    const std::unordered_map<std::string,
                             std::shared_ptr<CompositeSCCChartBase>>& sccs,
    const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(object,
      {"from_exact", "to_exact", "direction", "division_order", "charts",
       "matches", "tiles", "topology"}, label);
  ExactArmRequest request;
  request.from = parse_exact_path_rational(object.at("from_exact"), label);
  request.to = parse_exact_path_rational(object.at("to_exact"), label);
  request.topology = parse_exact_path_topology(object.at("topology"));
  std::vector<RetainedPlanChartBinding> bindings;
  const auto& raw_charts = as_array(object.at("charts"), label);
  if (raw_charts.empty())
    throw std::invalid_argument(std::string(label) +
                                " has no retained charts");
  bindings.reserve(raw_charts.size());
  request.charts.reserve(raw_charts.size());
  for (std::size_t index = 0; index < raw_charts.size(); ++index) {
    const auto& raw_chart = as_object(raw_charts[index], label);
    require_exact_keys(raw_chart,
        {"index", "chart", "identity", "center_exact", "scale_exact",
         "radius_exact", "singular_center", "prescriptions"}, label);
    if (as_u64(raw_chart.at("index"), label) != index)
      throw std::invalid_argument(std::string(label) +
                                  " chart bindings are not in index order");
    const auto handle = required_string(raw_chart, "chart");
    RetainedPlanChartBinding::Owner owner;
    if (handle.starts_with("c:")) {
      const auto found = charts.find(handle);
      if (found == charts.end())
        throw std::invalid_argument(
            std::string(label) +
            " references an absent prepared-chart owner");
      owner = found->second;
    } else if (handle.starts_with("scc:")) {
      const auto found = sccs.find(handle);
      if (found == sccs.end())
        throw std::invalid_argument(
            std::string(label) +
            " references an absent composite-SCC owner");
      owner = found->second;
    } else {
      throw std::invalid_argument(
          std::string(label) +
          " chart binding is neither a prepared-chart nor composite-SCC handle");
    }
    auto binding = bind_plan_chart(owner, request.topology);
    if (encode_plan_chart(binding, index) != raw_chart)
      throw std::invalid_argument(std::string(label) +
                                  " chart binding differs from its exact retained owner");
    request.charts.push_back(binding.geometry);
    bindings.push_back(std::move(binding));
  }
  // The plan is rederived below.  These fields are still parsed here to make
  // malformed scalar types fail before any planning work is attempted.
  (void)as_i32(object.at("direction"), label);
  (void)as_u32(object.at("division_order"), label);
  (void)as_array(object.at("matches"), label);
  (void)as_array(object.at("tiles"), label);
  return {std::move(request), std::move(bindings)};
}

std::shared_ptr<StoredTilePlan> restore_checkpoint_tile_plan_record(
    const json::value& raw,
    const std::unordered_map<std::string,
                             std::shared_ptr<PreparedChartBase>>& charts,
    const std::unordered_map<std::string,
                             std::shared_ptr<CompositeSCCChartBase>>& sccs) {
  const auto& object = as_object(raw, "checkpoint retained tile plan");
  if (required_string(object, "schema") ==
      kRetainedSingleArmTilePlanCheckpointSchema) {
    require_exact_keys(
        object,
        {"schema", "handle", "checkpoint_identity", "provenance_identity",
         "division_order", "arm_name", "arm", "elapsed_ms",
         "runtime_stats"},
        "checkpoint retained single-arm tile plan");
    const auto handle = required_string(object, "handle");
    const auto checkpoint_identity = required_string(
        object, "checkpoint_identity");
    const auto provenance_identity = required_string(
        object, "provenance_identity");
    const auto arm_name = required_string(object, "arm_name");
    if (handle.empty() || checkpoint_identity.empty() ||
        provenance_identity.empty() ||
        (arm_name != "lower" && arm_name != "upper"))
      throw std::invalid_argument(
          "checkpoint retained single-arm tile plan contains an invalid identity or arm name");
    const auto division_order = as_u32(
        object.at("division_order"),
        "checkpoint single-arm tile division order");
    auto [arm_request, bindings] = parse_checkpoint_retained_arm(
        object.at("arm"), charts, sccs,
        "checkpoint retained single tile arm");
    ExactPathPlanOptions options;
    options.division_order = division_order;
    auto exact = plan_exact_arm(arm_request, options);
    const std::string derived_name =
        exact.direction < 0 ? "lower" : "upper";
    if (derived_name != arm_name)
      throw std::invalid_argument(
          "checkpoint single-arm tile-plan name differs from its exact direction");
    RetainedArmPlan retained{std::move(exact), std::move(bindings)};
    if (encode_retained_arm(retained) != object.at("arm"))
      throw std::invalid_argument(
          "checkpoint single-arm tile intervals do not reproduce the exact planner result");
    json::object provenance{
        {"schema", kRetainedSingleArmTilePlanProvenanceSchema},
        {"checkpoint_identity", checkpoint_identity},
        {"division_order", division_order},
        {"arm_name", arm_name},
        {"arm", encode_retained_arm(retained)}};
    if (json::serialize(canonical_json_value(provenance)) !=
        provenance_identity)
      throw std::invalid_argument(
          "checkpoint single-arm tile-plan provenance identity is inconsistent");
    const auto elapsed_ms = checkpoint_nonnegative_double(
        object.at("elapsed_ms"),
        "checkpoint single-arm tile-plan elapsed time");
    const auto& stats = as_object(
        object.at("runtime_stats"),
        "checkpoint single-arm tile-plan runtime stats");
    require_exact_keys(
        stats,
        {"match_interval_queries", "tile_interval_queries",
         "lower_match_advances", "upper_match_advances", "integrations"},
        "checkpoint single-arm tile-plan runtime stats");
    auto plan = std::make_shared<StoredTilePlan>(
        handle, checkpoint_identity, provenance_identity, division_order,
        arm_name, std::move(retained), elapsed_ms);
    plan->restore_runtime_stats(
        as_u64(stats.at("match_interval_queries"),
               "checkpoint single-arm tile match queries"),
        as_u64(stats.at("tile_interval_queries"),
               "checkpoint single-arm tile interval queries"),
        as_u64(stats.at("lower_match_advances"),
               "checkpoint single-arm lower match advances"),
        as_u64(stats.at("upper_match_advances"),
               "checkpoint single-arm upper match advances"),
        as_u64(stats.at("integrations"),
               "checkpoint single-arm tile integrations"));
    return plan;
  }
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "division_order", "lower", "upper", "elapsed_ms",
       "runtime_stats"},
      "checkpoint retained tile plan");
  if (required_string(object, "schema") !=
      "diffexp2-retained-tile-plan-v2")
    throw std::invalid_argument(
        "unsupported retained tile-plan checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty())
    throw std::invalid_argument(
        "checkpoint retained tile plan contains an empty identity");
  const auto division_order = as_u32(
      object.at("division_order"), "checkpoint tile division order");
  auto [lower_request, lower_bindings] = parse_checkpoint_retained_arm(
      object.at("lower"), charts, sccs, "checkpoint lower tile arm");
  auto [upper_request, upper_bindings] = parse_checkpoint_retained_arm(
      object.at("upper"), charts, sccs, "checkpoint upper tile arm");
  if (lower_bindings.front().handle != upper_bindings.front().handle)
    throw std::invalid_argument(
        "checkpoint independent tile arms lost their shared anchor owner");
  ExactPathPlanOptions options;
  options.division_order = division_order;
  auto planned = plan_exact_independent_arms(
      lower_request, upper_request, options);
  RetainedArmPlan lower{std::move(planned.lower),
                        std::move(lower_bindings)};
  RetainedArmPlan upper{std::move(planned.upper),
                        std::move(upper_bindings)};
  if (encode_retained_arm(lower) != object.at("lower") ||
      encode_retained_arm(upper) != object.at("upper"))
    throw std::invalid_argument(
        "checkpoint tile intervals do not reproduce the exact planner result");
  json::object provenance{
      {"schema", "diffexp2-retained-exact-independent-arm-tile-plan-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"division_order", division_order},
      {"lower", encode_retained_arm(lower)},
      {"upper", encode_retained_arm(upper)}};
  if (json::serialize(canonical_json_value(provenance)) !=
      provenance_identity)
    throw std::invalid_argument(
        "checkpoint tile-plan provenance identity is inconsistent");
  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint tile-plan elapsed time");
  const auto& stats = as_object(object.at("runtime_stats"),
                                "checkpoint tile-plan runtime stats");
  require_exact_keys(stats,
      {"match_interval_queries", "tile_interval_queries",
       "lower_match_advances", "upper_match_advances", "integrations"},
      "checkpoint tile-plan runtime stats");
  auto plan = std::make_shared<StoredTilePlan>(
      handle, checkpoint_identity, provenance_identity, division_order,
      std::move(lower), std::move(upper), elapsed_ms);
  plan->restore_runtime_stats(
      as_u64(stats.at("match_interval_queries"),
             "checkpoint tile match queries"),
      as_u64(stats.at("tile_interval_queries"),
             "checkpoint tile interval queries"),
      as_u64(stats.at("lower_match_advances"),
             "checkpoint lower match advances"),
      as_u64(stats.at("upper_match_advances"),
             "checkpoint upper match advances"),
      as_u64(stats.at("integrations"),
             "checkpoint tile integrations"));
  return plan;
}

std::shared_ptr<StoredEndpointResult>
restore_checkpoint_planned_endpoint_record(
    const json::value& raw, const std::string& expected_domain,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& local) {
  const auto& object = as_object(
      raw, "checkpoint retained plan-bound endpoint");
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "source", "approach_direction", "derived_rim",
       "cancellation_mode", "analytic_metadata", "result", "elapsed_ms",
       "runtime_stats"},
      "checkpoint retained plan-bound endpoint");
  if (required_string(object, "schema") !=
      "diffexp2-retained-plan-bound-endpoint-result-v1")
    throw std::invalid_argument(
        "unsupported retained plan-bound endpoint checkpoint schema");
  if (!plan || !local)
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint lost a strongly owned plan or local");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty())
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint contains an empty identity");
  if (std::string(local->scalar_domain()) != expected_domain ||
      (expected_domain != "rational" && expected_domain != "acb"))
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint coefficient domain differs from its session/local");

  const auto& source = as_object(
      object.at("source"), "checkpoint plan-bound endpoint source");
  const auto arm_name = required_string(source, "arm");
  const auto binding = resolve_planned_endpoint_binding(
      plan, arm_name, local);
  if (source != binding.source)
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint source differs from its exact plan/local owners");
  const auto approach_direction = as_i32(
      object.at("approach_direction"),
      "checkpoint plan-bound endpoint approach direction");
  if (approach_direction != binding.approach_direction)
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint approach side differs from its exact arm/chart orientation");
  std::optional<std::int32_t> derived_rim;
  if (!object.at("derived_rim").is_null()) {
    derived_rim = as_i32(object.at("derived_rim"),
                         "checkpoint plan-bound endpoint rim");
    if (*derived_rim != -1 && *derived_rim != 1)
      throw std::invalid_argument(
          "checkpoint plan-bound endpoint rim must be +1 or -1");
  }
  if (derived_rim != binding.rim)
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint rim differs from its final-chart prescriptions");
  const auto cancellation_mode = required_string(
      object, "cancellation_mode");
  if (cancellation_mode != "exact-coefficient-field" &&
      cancellation_mode != "exact-or-acb-singleton")
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint cancellation mode is unsupported");
  auto analytic_metadata = as_object(
      object.at("analytic_metadata"),
      "checkpoint plan-bound endpoint analytic metadata");
  validate_checkpoint_exact_analytic_metadata(analytic_metadata);
  if (analytic_metadata != local->exact_analytic_metadata())
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint analytic metadata differs from its local owner");
  const auto provenance = planned_endpoint_provenance(
      checkpoint_identity, binding, cancellation_mode, analytic_metadata);
  if (json::serialize(canonical_json_value(provenance)) !=
      provenance_identity)
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint provenance identity is inconsistent");

  const auto& raw_result = as_object(
      object.at("result"), "checkpoint plan-bound endpoint result");
  require_exact_keys(
      raw_result,
      {"values", "dropped_regulated_sectors",
       "cancelled_divergent_coefficients", "imaginary_sign"},
      "checkpoint plan-bound endpoint result");
  EndpointLimitResult result;
  for (const auto& raw_value : as_array(
           raw_result.at("values"),
           "checkpoint plan-bound endpoint values"))
    result.values.push_back(parse_checkpoint_epsilon_frame<ComplexBall>(
        raw_value, "checkpoint plan-bound endpoint value"));
  if (result.values.empty())
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint result has no components");
  const auto checked_size = [](const json::value& value,
                               const char* label) {
    const auto parsed = as_u64(value, label);
    if (parsed > std::numeric_limits<std::size_t>::max())
      throw std::invalid_argument(std::string(label) + " exceeds size_t");
    return static_cast<std::size_t>(parsed);
  };
  result.dropped_regulated_sectors = checked_size(
      raw_result.at("dropped_regulated_sectors"),
      "checkpoint plan-bound dropped regulated sectors");
  result.cancelled_divergent_coefficients = checked_size(
      raw_result.at("cancelled_divergent_coefficients"),
      "checkpoint plan-bound cancelled divergent coefficients");
  std::optional<std::int32_t> result_rim;
  if (!raw_result.at("imaginary_sign").is_null()) {
    result_rim = as_i32(raw_result.at("imaginary_sign"),
                        "checkpoint plan-bound endpoint result rim");
    if (*result_rim != -1 && *result_rim != 1)
      throw std::invalid_argument(
          "checkpoint plan-bound endpoint result rim must be +1 or -1");
  }
  if (result_rim != binding.rim)
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint result rim differs from its derived branch");
  // The endpoint kernel's limit classification is branch independent.  It
  // currently carries an internal integer diagnostic, while the plan-bound
  // public/checkpoint contract deliberately keeps an unprescribed rim null.
  result.imaginary_sign = binding.rim.value_or(1);
  (void)endpoint_value_window(result);

  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"),
      "checkpoint plan-bound endpoint elapsed time");
  const auto& stats = as_object(
      object.at("runtime_stats"),
      "checkpoint plan-bound endpoint runtime stats");
  require_exact_keys(stats, {"exports", "export_ms"},
                     "checkpoint plan-bound endpoint runtime stats");
  const auto exports = as_u64(
      stats.at("exports"), "checkpoint plan-bound endpoint exports");
  const auto export_ms = checkpoint_nonnegative_double(
      stats.at("export_ms"),
      "checkpoint plan-bound endpoint export time");
  auto endpoint = std::make_shared<StoredEndpointResult>(
      handle, checkpoint_identity, provenance_identity,
      local->handle(), local->source_chart(),
      local->source_operator_identity(), local->checkpoint_identity(),
      local->scalar_domain(), approach_direction, std::nullopt,
      cancellation_mode, std::move(analytic_metadata), std::move(result),
      elapsed_ms, binding.source, binding.rim,
      plan, local);
  endpoint->restore_runtime_stats(exports, export_ms);
  return endpoint;
}

std::size_t checkpoint_size_t(const json::value& raw, const char* label);

std::shared_ptr<StoredPlannedMatchHop>
restore_checkpoint_planned_match_hop_record(
    const json::value& raw,
    const std::shared_ptr<StoredTilePlan>& plan,
    std::vector<std::shared_ptr<StoredLocalBase>> basis,
    std::shared_ptr<StoredLocalBase> incoming,
    const std::string& source_session_configuration_identity) {
  const auto& object = as_object(
      raw, "checkpoint retained planned match hop");
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "handoff", "native_match", "elapsed_ms", "runtime_stats"},
      "checkpoint retained planned match hop");
  if (required_string(object, "schema") !=
      "diffexp2-retained-planned-match-hop-v2")
    throw std::invalid_argument(
        "unsupported retained planned-match checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty() || !plan || !incoming || basis.empty())
    throw std::invalid_argument(
        "checkpoint planned match lost an identity or strong owner");
  const auto& handoff = as_object(
      object.at("handoff"), "checkpoint planned-match handoff");
  require_exact_keys(
      handoff,
      {"schema", "tile_plan", "tile_plan_checkpoint_identity",
       "tile_plan_provenance_identity", "arm", "match", "geometry",
       "producing", "receiving", "result_checkpoint_identity",
       "native_match_provenance_identity", "advance"},
      "checkpoint planned-match handoff");
  if (required_string(handoff, "schema") !=
          "diffexp2-retained-exact-plan-match-hop-v1" ||
      required_string(handoff, "tile_plan") != plan->handle() ||
      required_string(handoff, "tile_plan_checkpoint_identity") !=
          plan->checkpoint_identity() ||
      required_string(handoff, "tile_plan_provenance_identity") !=
          plan->provenance_identity() ||
      required_string(handoff, "result_checkpoint_identity") !=
          checkpoint_identity)
    throw std::invalid_argument(
        "checkpoint planned-match handoff lost its plan/checkpoint provenance");
  const auto arm_name = required_string(handoff, "arm");
  const auto match_index = checkpoint_size_t(
      handoff.at("match"), "checkpoint planned-match index");
  const auto& arm = plan->arm(arm_name);
  if (match_index >= arm.exact.matches.size())
    throw std::invalid_argument(
        "checkpoint planned-match index lies outside its retained plan");
  const auto& exact_match = arm.exact.matches[match_index];
  const auto& producing = arm.charts.at(exact_match.producing_chart);
  const auto& receiving = arm.charts.at(exact_match.receiving_chart);
  incoming->require_exact_plan_binding(
      producing.geometry, producing.prescriptions,
      "checkpoint planned-match incoming owner");
  for (const auto& local : basis)
    local->require_exact_plan_binding(
        receiving.geometry, receiving.prescriptions,
        "checkpoint planned-match basis owner");

  const auto& native_record = as_object(
      object.at("native_match"), "checkpoint embedded native match");
  const auto native_schema = required_string(native_record, "schema");
  std::shared_ptr<StoredMatchBase> native_match;
  if (native_schema == "diffexp2-retained-exact-rational-match-v2") {
    native_match = restore_checkpoint_exact_match_record(
        native_record, basis, incoming);
  } else if (native_schema == "diffexp2-retained-acb-match-v2") {
    const auto saturation = native_acb_saturation_binding(
        plan, source_session_configuration_identity, arm_name, match_index,
        checkpoint_identity);
    const std::optional<json::object> expected_singular_request =
        saturation.request_key == "native_singular_scc_saturation"
            ? std::optional<json::object>(saturation.request)
            : std::nullopt;
    native_match = restore_checkpoint_acb_match_record(
        native_record, source_session_configuration_identity,
        expected_singular_request);
    const auto cross_check = [](const json::object& source,
                                const std::shared_ptr<StoredLocalBase>& owner,
                                const char* label) {
      if (!owner || required_string(source, "local") != owner->handle() ||
          required_string(source, "chart") != owner->source_chart() ||
          required_string(source, "source_operator_identity") !=
              owner->source_operator_identity() ||
          required_string(source, "checkpoint_identity") !=
              owner->checkpoint_identity() ||
          source.at("analytic_metadata") != owner->exact_analytic_metadata())
        throw std::invalid_argument(
            std::string("checkpoint planned Acb match ") + label +
            " disagrees with its strong local owner");
    };
    const auto& sources = as_array(
        native_record.at("basis_sources"),
        "checkpoint planned Acb basis sources");
    if (sources.size() != basis.size())
      throw std::invalid_argument(
          "checkpoint planned Acb basis ownership differs from its dimension");
    for (std::size_t column = 0; column < basis.size(); ++column)
      cross_check(as_object(sources[column],
                            "checkpoint planned Acb basis source"),
                  basis[column], "basis source");
    cross_check(as_object(native_record.at("incoming_source"),
                          "checkpoint planned Acb incoming source"),
                incoming, "incoming source");
  } else {
    throw std::invalid_argument(
        "checkpoint planned hop embeds an unsupported native match kind");
  }
  if (native_match->handle() != handle ||
      native_match->checkpoint_record() != object.at("native_match"))
    throw std::invalid_argument(
        "checkpoint planned hop embedded match does not reproduce its exact payload");
  const auto native_summary = native_match->summary();
  if (required_string(native_summary, "checkpoint_identity") !=
          checkpoint_identity ||
      required_string(native_summary, "physical_match_point_exact") !=
          exact_match.physical.str())
    throw std::invalid_argument(
        "checkpoint planned hop embedded match changed its checkpoint/physical point");
  std::vector<std::string> basis_handles;
  basis_handles.reserve(basis.size());
  for (const auto& local : basis) basis_handles.push_back(local->handle());
  auto expected_handoff = planned_match_handoff_record(
      plan, arm_name, match_index, basis_handles, basis, incoming->handle(),
      incoming, checkpoint_identity, native_summary);
  if (handoff != expected_handoff ||
      json::serialize(canonical_json_value(expected_handoff)) !=
          provenance_identity)
    throw std::invalid_argument(
        "checkpoint planned-match handoff/provenance differs from its exact owners");
  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint planned-match elapsed time");
  const auto& stats = as_object(
      object.at("runtime_stats"),
      "checkpoint planned-match runtime stats");
  require_exact_keys(stats, {"materializations"},
                     "checkpoint planned-match runtime stats");
  auto hop = std::make_shared<StoredPlannedMatchHop>(
      std::move(native_match), checkpoint_identity, provenance_identity,
      expected_handoff, elapsed_ms, plan, std::move(basis), incoming);
  hop->restore_runtime_stats(as_u64(
      stats.at("materializations"),
      "checkpoint planned-match materializations"));
  return hop;
}

std::size_t checkpoint_size_t(const json::value& raw, const char* label) {
  const auto value = as_u64(raw, label);
  if (value > std::numeric_limits<std::size_t>::max())
    throw std::invalid_argument(std::string(label) + " exceeds size_t");
  return static_cast<std::size_t>(value);
}

std::shared_ptr<StoredTransportArmState>
restore_checkpoint_transport_arm_state_record(
    const json::value& raw,
    const std::unordered_map<std::string,
                             std::shared_ptr<StoredTilePlan>>& plans,
    const std::unordered_map<std::string,
                             std::shared_ptr<StoredLocalBase>>& locals,
    const std::unordered_map<std::string,
                             std::shared_ptr<StoredMatchBase>>& matches) {
  const auto& object = as_object(
      raw, "checkpoint retained transport-arm state");
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "provenance", "elapsed_ms", "runtime_stats"},
      "checkpoint retained transport-arm state");
  if (required_string(object, "schema") !=
      "diffexp2-retained-transport-arm-state-v1")
    throw std::invalid_argument(
        "unsupported retained transport-arm checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty())
    throw std::invalid_argument(
        "checkpoint transport-arm state contains an empty identity");
  const auto& provenance = as_object(
      object.at("provenance"),
      "checkpoint transport-arm state provenance");
  require_exact_keys(
      provenance,
      {"schema", "checkpoint_identity", "tile_plan", "arm", "anchor",
       "receiving_basis", "matches", "tile_sources", "final_local",
       "epsilon", "refinement"},
      "checkpoint transport-arm state provenance");
  if (required_string(provenance, "schema") !=
          "diffexp2-retained-native-transport-arm-state-v1" ||
      required_string(provenance, "checkpoint_identity") !=
          checkpoint_identity ||
      json::serialize(canonical_json_value(provenance)) !=
          provenance_identity)
    throw std::invalid_argument(
        "checkpoint transport-arm provenance identity is inconsistent");

  const auto& plan_reference = as_object(
      provenance.at("tile_plan"),
      "checkpoint transport-arm tile-plan reference");
  require_exact_keys(
      plan_reference,
      {"handle", "checkpoint_identity", "provenance_identity"},
      "checkpoint transport-arm tile-plan reference");
  const auto plan_found = plans.find(
      required_string(plan_reference, "handle"));
  if (plan_found == plans.end() || !plan_found->second ||
      required_string(plan_reference, "checkpoint_identity") !=
          plan_found->second->checkpoint_identity() ||
      required_string(plan_reference, "provenance_identity") !=
          plan_found->second->provenance_identity())
    throw std::invalid_argument(
        "checkpoint transport-arm state lost its retained tile plan");
  const auto plan = plan_found->second;
  const auto arm_name = required_string(provenance, "arm");
  (void)plan->arm(arm_name);

  const auto resolve_local_reference = [&](const json::value& raw_reference,
                                           const char* label,
                                           bool has_tile) {
    const auto& reference = as_object(raw_reference, label);
    if (has_tile)
      require_exact_keys(
          reference,
          {"tile", "local", "chart", "source_operator_identity",
           "checkpoint_identity", "coefficient_domain"},
          label);
    else
      require_exact_keys(
          reference,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity", "coefficient_domain"},
          label);
    const auto found = locals.find(required_string(reference, "local"));
    if (found == locals.end() || !found->second ||
        required_string(reference, "chart") !=
            found->second->source_chart() ||
        required_string(reference, "source_operator_identity") !=
            found->second->source_operator_identity() ||
        required_string(reference, "checkpoint_identity") !=
            found->second->checkpoint_identity() ||
        required_string(reference, "coefficient_domain") !=
            found->second->scalar_domain())
      throw std::invalid_argument(
          std::string(label) + " disagrees with its retained local owner");
    return found->second;
  };

  auto anchor = resolve_local_reference(
      provenance.at("anchor"),
      "checkpoint transport-arm anchor reference", false);
  std::vector<std::vector<std::shared_ptr<StoredLocalBase>>> basis;
  for (const auto& raw_basis : as_array(
           provenance.at("receiving_basis"),
           "checkpoint transport-arm receiving bases")) {
    std::vector<std::shared_ptr<StoredLocalBase>> columns;
    for (const auto& raw_column : as_array(
             raw_basis, "checkpoint transport-arm receiving basis"))
      columns.push_back(resolve_local_reference(
          raw_column, "checkpoint transport-arm basis reference", false));
    if (columns.empty())
      throw std::invalid_argument(
          "checkpoint transport-arm receiving basis cannot be empty");
    basis.push_back(std::move(columns));
  }

  std::vector<std::shared_ptr<StoredPlannedMatchHop>> planned_matches;
  const auto& raw_matches = as_array(
      provenance.at("matches"), "checkpoint transport-arm matches");
  planned_matches.reserve(raw_matches.size());
  for (std::size_t index = 0; index < raw_matches.size(); ++index) {
    const auto& reference = as_object(
        raw_matches[index], "checkpoint transport-arm match reference");
    require_exact_keys(
        reference,
        {"index", "match", "checkpoint_identity", "provenance_identity"},
        "checkpoint transport-arm match reference");
    if (checkpoint_size_t(reference.at("index"),
                          "checkpoint transport-arm match index") != index)
      throw std::invalid_argument(
          "checkpoint transport-arm matches are out of order");
    const auto found = matches.find(required_string(reference, "match"));
    if (found == matches.end())
      throw std::invalid_argument(
          "checkpoint transport-arm state lost a planned match");
    auto match = std::dynamic_pointer_cast<StoredPlannedMatchHop>(
        found->second);
    if (!match ||
        required_string(reference, "checkpoint_identity") !=
            match->checkpoint_identity() ||
        required_string(reference, "provenance_identity") !=
            match->provenance_identity())
      throw std::invalid_argument(
          "checkpoint transport-arm match reference is inconsistent");
    planned_matches.push_back(std::move(match));
  }

  std::vector<std::shared_ptr<StoredLocalBase>> tile_sources;
  const auto& raw_sources = as_array(
      provenance.at("tile_sources"),
      "checkpoint transport-arm tile sources");
  tile_sources.reserve(raw_sources.size());
  for (std::size_t tile = 0; tile < raw_sources.size(); ++tile) {
    const auto& reference = as_object(
        raw_sources[tile], "checkpoint transport-arm tile-source reference");
    if (checkpoint_size_t(reference.at("tile"),
                          "checkpoint transport-arm tile index") != tile)
      throw std::invalid_argument(
          "checkpoint transport-arm tile sources are out of order");
    tile_sources.push_back(resolve_local_reference(
        raw_sources[tile],
        "checkpoint transport-arm tile-source reference", true));
  }
  auto final_local = resolve_local_reference(
      provenance.at("final_local"),
      "checkpoint transport-arm final-local reference", false);
  if (tile_sources.empty() || tile_sources.back().get() != final_local.get())
    throw std::invalid_argument(
        "checkpoint transport-arm final local differs from its last tile source");

  const auto& epsilon = as_object(
      provenance.at("epsilon"), "checkpoint transport-arm epsilon contract");
  require_exact_keys(
      epsilon,
      {"min", "max", "required_complete_max",
       "match_required_complete_max"},
      "checkpoint transport-arm epsilon contract");
  EpsilonWindow work_epsilon{
      as_i32(epsilon.at("min"), "checkpoint transport-arm epsilon minimum"),
      as_i32(epsilon.at("max"), "checkpoint transport-arm epsilon maximum")};
  (void)work_epsilon.width();
  const auto public_required_complete_max = as_i32(
      epsilon.at("required_complete_max"),
      "checkpoint transport-arm public epsilon maximum");
  const auto match_required_complete_max = as_i32(
      epsilon.at("match_required_complete_max"),
      "checkpoint transport-arm match epsilon maximum");
  const auto refinement = as_object(
      provenance.at("refinement"),
      "checkpoint transport-arm refinement policy");
  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint transport-arm elapsed time");
  auto state = std::make_shared<StoredTransportArmState>(
      handle, checkpoint_identity, arm_name, plan, anchor,
      std::move(basis), std::move(planned_matches),
      std::move(tile_sources), work_epsilon,
      public_required_complete_max, match_required_complete_max,
      refinement, elapsed_ms);
  if (state->provenance_identity() != provenance_identity)
    throw std::invalid_argument(
        "restored transport-arm state changed its exact provenance");
  const auto& stats = as_object(
      object.at("runtime_stats"),
      "checkpoint transport-arm runtime stats");
  require_exact_keys(stats, {"stats_queries"},
                     "checkpoint transport-arm runtime stats");
  state->restore_runtime_stats(as_u64(
      stats.at("stats_queries"),
      "checkpoint transport-arm statistics queries"));
  if (state->checkpoint_record() != raw)
    throw std::invalid_argument(
        "restored transport-arm state does not reproduce its exact retained state");
  return state;
}

std::shared_ptr<StoredLineResult> restore_checkpoint_line_result_record(
    const json::value& raw,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& local) {
  const auto& object = as_object(raw, "checkpoint retained line result");
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "arm", "tile", "interval", "source", "result", "elapsed_ms",
       "runtime_stats"},
      "checkpoint retained line result");
  if (required_string(object, "schema") !=
      "diffexp2-retained-line-result-v2")
    throw std::invalid_argument(
        "unsupported retained line-result checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  const auto arm_name = required_string(object, "arm");
  const auto tile_index = checkpoint_size_t(
      object.at("tile"), "checkpoint line tile index");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty())
    throw std::invalid_argument(
        "checkpoint retained line result contains an empty identity");
  const auto& arm = plan->arm(arm_name);
  if (tile_index >= arm.exact.tiles.size())
    throw std::invalid_argument(
        "checkpoint line tile lies outside its retained plan owner");
  const auto expected_interval = encode_plan_tile(arm, tile_index);
  const auto interval = as_object(object.at("interval"),
                                  "checkpoint line interval");
  if (interval != expected_interval)
    throw std::invalid_argument(
        "checkpoint line interval differs from its exact tile-plan owner");
  const auto& source = as_object(object.at("source"),
                                 "checkpoint line source");
  require_exact_keys(
      source,
      {"tile_plan", "tile_plan_checkpoint_identity", "local", "chart",
       "source_operator_identity", "local_checkpoint_identity",
       "coefficient_domain", "analytic_metadata"},
      "checkpoint line source");
  if (!plan || !local ||
      required_string(source, "tile_plan") != plan->handle() ||
      required_string(source, "tile_plan_checkpoint_identity") !=
          plan->checkpoint_identity() ||
      required_string(source, "local") != local->handle() ||
      required_string(source, "chart") != local->source_chart() ||
      required_string(source, "source_operator_identity") !=
          local->source_operator_identity() ||
      required_string(source, "local_checkpoint_identity") !=
          local->checkpoint_identity() ||
      required_string(source, "coefficient_domain") !=
          local->scalar_domain() ||
      source.at("analytic_metadata") != local->exact_analytic_metadata())
    throw std::invalid_argument(
        "checkpoint line source provenance disagrees with its strong owners");
  if (std::string(local->scalar_domain()) == "symbolic")
    throw std::invalid_argument(
        "checkpoint line result cannot own a symbolic local");
  validate_checkpoint_exact_analytic_metadata(
      source.at("analytic_metadata"));
  const auto& tile = arm.exact.tiles[tile_index];
  const auto& binding = arm.charts.at(tile.chart);
  if (local->source_chart() != binding.handle)
    throw std::invalid_argument(
        "checkpoint line local does not own the tile's retained chart");

  const auto& raw_result = as_object(object.at("result"),
                                     "checkpoint line result state");
  require_exact_keys(raw_result, {"value", "scope", "imaginary_sign",
                                  "diagnostics"},
                     "checkpoint line result state");
  const auto scope = required_string(raw_result, "scope");
  StoredLineIntegral result;
  if (scope == "stored_truncation")
    result.scope = LineIntegrationScope::StoredTruncation;
  else if (scope == "full_local_with_certified_tail")
    result.scope = LineIntegrationScope::FullLocalWithCertifiedTail;
  else
    throw std::invalid_argument(
        "checkpoint line result has an unsupported integration scope");
  result.value = parse_checkpoint_epsilon_vector(
      raw_result.at("value"), "checkpoint line epsilon vector");
  if (!raw_result.at("imaginary_sign").is_null()) {
    result.imaginary_sign = as_i32(raw_result.at("imaginary_sign"),
                                   "checkpoint line rim");
    if (*result.imaginary_sign != -1 && *result.imaginary_sign != 1)
      throw std::invalid_argument(
          "checkpoint line rim must be +1 or -1");
  }
  const auto expected_rim = exact_plan_rim(
      binding.prescriptions, binding.geometry.scale);
  if (result.imaginary_sign != expected_rim)
    throw std::invalid_argument(
        "checkpoint line rim differs from its exact branch prescriptions");
  if (result.value.dimension !=
      as_u32(local->summary().at("dimension"),
             "checkpoint line source dimension"))
    throw std::invalid_argument(
        "checkpoint line result dimension differs from its local owner");
  const auto& diagnostics = as_object(raw_result.at("diagnostics"),
                                      "checkpoint line diagnostics");
  const bool has_tail_requested =
      diagnostics.if_contains("tail_certificate_requested") != nullptr;
  const bool has_tail_status =
      diagnostics.if_contains("tail_certificate_status") != nullptr;
  const bool has_tail_witness =
      diagnostics.if_contains("tail_witness_radius_exact") != nullptr;
  if (has_tail_requested != has_tail_status ||
      has_tail_requested != has_tail_witness)
    throw std::invalid_argument(
        "checkpoint line tail diagnostics are incomplete");
  if (has_tail_requested)
    require_exact_keys(
        diagnostics,
        {"input_monomial_cells", "grouped_monomials", "zero_groups_skipped",
         "cancelled_divergent_groups", "primitive_evaluations",
         "primitive_component_applications", "primitive_component_reuses",
         "has_center_endpoint", "tail_certificate_requested",
         "tail_certificate_status", "tail_witness_radius_exact", "detail"},
        "checkpoint line diagnostics");
  else
    require_exact_keys(
        diagnostics,
        {"input_monomial_cells", "grouped_monomials", "zero_groups_skipped",
         "cancelled_divergent_groups", "primitive_evaluations",
         "primitive_component_applications", "primitive_component_reuses",
         "has_center_endpoint", "detail"},
        "checkpoint line diagnostics");
  result.diagnostics.input_monomial_cells = checkpoint_size_t(
      diagnostics.at("input_monomial_cells"),
      "checkpoint line input monomials");
  result.diagnostics.grouped_monomials = checkpoint_size_t(
      diagnostics.at("grouped_monomials"),
      "checkpoint line grouped monomials");
  result.diagnostics.zero_groups_skipped = checkpoint_size_t(
      diagnostics.at("zero_groups_skipped"),
      "checkpoint line skipped zero groups");
  result.diagnostics.cancelled_divergent_groups = checkpoint_size_t(
      diagnostics.at("cancelled_divergent_groups"),
      "checkpoint line cancelled divergent groups");
  result.diagnostics.primitive_evaluations = checkpoint_size_t(
      diagnostics.at("primitive_evaluations"),
      "checkpoint line primitive evaluations");
  result.diagnostics.primitive_component_applications = checkpoint_size_t(
      diagnostics.at("primitive_component_applications"),
      "checkpoint line primitive component applications");
  result.diagnostics.primitive_component_reuses = checkpoint_size_t(
      diagnostics.at("primitive_component_reuses"),
      "checkpoint line primitive component reuses");
  if (!diagnostics.at("has_center_endpoint").is_bool())
    throw std::invalid_argument(
        "checkpoint line center-endpoint flag must be Boolean");
  result.diagnostics.has_center_endpoint =
      diagnostics.at("has_center_endpoint").as_bool();
  if (has_tail_requested) {
    if (!diagnostics.at("tail_certificate_requested").is_bool())
      throw std::invalid_argument(
          "checkpoint line tail-request flag must be Boolean");
    result.diagnostics.tail_certificate_requested =
        diagnostics.at("tail_certificate_requested").as_bool();
    result.diagnostics.tail_certificate_status = required_string(
        diagnostics, "tail_certificate_status");
    if (!diagnostics.at("tail_witness_radius_exact").is_null())
      result.diagnostics.tail_witness_radius_exact = required_string(
          diagnostics, "tail_witness_radius_exact");
  }
  result.diagnostics.detail = required_string(diagnostics, "detail");

  const auto& error = result.value.error;
  const auto& line_epsilon = result.value.epsilon;
  const auto& tail = result.diagnostics;
  if (!tail.tail_witness_radius_exact.empty()) {
    try {
      const Rational witness(tail.tail_witness_radius_exact);
      const auto begin_modulus = tile.local_begin.sign() < 0
          ? -tile.local_begin : tile.local_begin;
      const auto end_modulus = tile.local_end.sign() < 0
          ? -tile.local_end : tile.local_end;
      const auto outer = begin_modulus < end_modulus
          ? end_modulus : begin_modulus;
      if (!(outer < witness) || !(witness < binding.geometry.radius))
        throw std::invalid_argument("tail witness lies outside its annulus");
    } catch (const std::invalid_argument&) {
      throw std::invalid_argument(
          "checkpoint line tail witness radius must be an exact rational "
          "strictly outside its tile and inside its chart");
    }
  }
  if (result.scope == LineIntegrationScope::StoredTruncation) {
    if (!error.empty() || error.guarantee != ErrorGuarantee::None)
      throw std::invalid_argument(
          "checkpoint stored-truncation line cannot carry an error "
          "envelope or guarantee");
    if (tail.tail_certificate_requested) {
      if (tail.tail_certificate_status != "unsupported" &&
          tail.tail_certificate_status != "inconclusive")
        throw std::invalid_argument(
            "checkpoint stored-truncation line has an inconsistent "
            "tail-certificate status");
    } else if (tail.tail_certificate_status != "not-requested" ||
               !tail.tail_witness_radius_exact.empty()) {
      throw std::invalid_argument(
          "checkpoint stored-truncation line has unsolicited "
          "tail-certificate diagnostics");
    }
  } else {
    if (error.empty() || error.guarantee != ErrorGuarantee::Certified ||
        error.frame.min_power != line_epsilon.min_power ||
        error.frame.complete_max != line_epsilon.complete_max ||
        error.provenance.empty())
      throw std::invalid_argument(
          "checkpoint full-local line requires a frame-aligned certified "
          "error envelope");
    if (!tail.tail_certificate_requested ||
        tail.tail_certificate_status != "certified" ||
        tail.tail_witness_radius_exact.empty())
      throw std::invalid_argument(
          "checkpoint full-local line requires complete certified-tail "
          "diagnostics");
  }

  json::object legacy_provenance{
      {"schema",
       "diffexp2-retained-native-stored-truncation-physical-tile-integral-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
      {"arm", arm_name}, {"tile", tile_index},
      {"interval", interval},
      {"source", json::object{
           {"local", local->handle()}, {"chart", local->source_chart()},
           {"checkpoint_identity", local->checkpoint_identity()}}},
      {"epsilon", json::object{{"min", result.value.epsilon.min_power},
                                {"max", result.value.epsilon.complete_max}}},
      {"scope", "stored_truncation"},
      {"error_guarantee", "none"}};
  json::object provenance{
      {"schema",
       "diffexp2-retained-native-physical-tile-integral-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
      {"arm", arm_name}, {"tile", tile_index},
      {"interval", interval},
      {"source", json::object{
           {"local", local->handle()}, {"chart", local->source_chart()},
           {"checkpoint_identity", local->checkpoint_identity()}}},
      {"epsilon", json::object{{"min", result.value.epsilon.min_power},
                                {"max", result.value.epsilon.complete_max}}},
      {"tail_certificate_requested",
       result.diagnostics.tail_certificate_requested},
      {"tail_certificate_status",
       result.diagnostics.tail_certificate_status},
      {"tail_witness_radius_exact",
       result.diagnostics.tail_witness_radius_exact.empty()
           ? json::value(nullptr)
           : json::value(result.diagnostics.tail_witness_radius_exact)},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)},
      {"error_provenance", result.value.error.provenance}};
  const auto current_identity =
      json::serialize(canonical_json_value(provenance));
  const auto legacy_identity =
      json::serialize(canonical_json_value(legacy_provenance));
  if (provenance_identity != current_identity &&
      (has_tail_requested ||
       result.scope != LineIntegrationScope::StoredTruncation ||
       provenance_identity != legacy_identity))
    throw std::invalid_argument(
        "checkpoint line provenance identity is inconsistent");
  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint line elapsed time");
  const auto& stats = as_object(object.at("runtime_stats"),
                                "checkpoint line runtime stats");
  require_exact_keys(stats, {"exports", "export_ms"},
                     "checkpoint line runtime stats");
  auto stored = std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, arm_name, tile_index,
      interval, local->checkpoint_identity(), std::move(result), elapsed_ms,
      plan, local);
  stored->restore_runtime_stats(
      as_u64(stats.at("exports"), "checkpoint line exports"),
      checkpoint_nonnegative_double(stats.at("export_ms"),
                                    "checkpoint line export time"));
  return stored;
}

std::shared_ptr<StoredLineResult> restore_checkpoint_line_aggregate_record(
    const json::value& raw,
    const std::shared_ptr<StoredTilePlan>& plan,
    std::vector<std::shared_ptr<StoredLocalBase>> local_owners) {
  const auto& object = as_object(raw, "checkpoint retained line aggregate");
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "provenance", "result", "elapsed_ms", "runtime_stats"},
      "checkpoint retained line aggregate");
  if (required_string(object, "schema") !=
      "diffexp2-retained-line-aggregate-v1")
    throw std::invalid_argument(
        "unsupported retained line-aggregate checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty() || !plan || local_owners.empty())
    throw std::invalid_argument(
        "checkpoint line aggregate lost an identity or strong owner");
  local_owners = unique_line_local_owners(local_owners);

  const auto& provenance = as_object(
      object.at("provenance"), "checkpoint line aggregate provenance");
  require_exact_keys(
      provenance,
      {"schema", "checkpoint_identity", "arm", "interval", "source",
       "aggregate", "epsilon", "scope", "error_guarantee"},
      "checkpoint line aggregate provenance");
  if (required_string(provenance, "schema") !=
          "diffexp2-retained-native-line-aggregate-v1" ||
      required_string(provenance, "checkpoint_identity") !=
          checkpoint_identity ||
      json::serialize(canonical_json_value(provenance)) !=
          provenance_identity)
    throw std::invalid_argument(
        "checkpoint line aggregate provenance identity is inconsistent");
  const auto& source = as_object(
      provenance.at("source"), "checkpoint line aggregate source");
  require_exact_keys(source,
      {"tile_plan", "tile_plan_checkpoint_identity", "locals"},
      "checkpoint line aggregate source");
  if (required_string(source, "tile_plan") != plan->handle() ||
      required_string(source, "tile_plan_checkpoint_identity") !=
          plan->checkpoint_identity() ||
      source.at("locals") != line_aggregate_source_records(local_owners))
    throw std::invalid_argument(
        "checkpoint line aggregate source disagrees with its strong owners");
  for (const auto& owner : local_owners) {
    if (std::string(owner->scalar_domain()) == "symbolic")
      throw std::invalid_argument(
          "checkpoint line aggregate cannot own a symbolic local");
    validate_checkpoint_exact_analytic_metadata(
        owner->exact_analytic_metadata());
  }
  const auto& aggregate = as_object(
      provenance.at("aggregate"), "checkpoint line aggregate recipe");
  std::optional<WholeArmEpsilonContract> whole_arm_epsilon_contract;
  if (const auto* raw_contract = aggregate.if_contains("epsilon_contract"))
    whole_arm_epsilon_contract = parse_whole_arm_epsilon_contract(
        *raw_contract, "checkpoint whole-arm epsilon contract");
  const auto& components = as_array(
      aggregate.at("components"), "checkpoint line aggregate components");
  if (checkpoint_size_t(aggregate.at("component_count"),
                        "checkpoint line aggregate component count") !=
      components.size() ||
      components.empty())
    throw std::invalid_argument(
        "checkpoint line aggregate component manifest is inconsistent");
  for (std::size_t index = 0; index < components.size(); ++index) {
    const auto& component = as_object(
        components[index], "checkpoint line aggregate component");
    require_exact_keys(
        component,
        {"index", "sign", "checkpoint_identity", "provenance_identity",
         "scope", "source", "interval", "epsilon", "error"},
        "checkpoint line aggregate component");
    const auto sign = as_i32(component.at("sign"),
                             "checkpoint line aggregate sign");
    if (checkpoint_size_t(component.at("index"),
                          "checkpoint line aggregate component index") !=
            index ||
        (sign != -1 && sign != 1) ||
        required_string(component, "checkpoint_identity").empty() ||
        required_string(component, "provenance_identity").empty())
      throw std::invalid_argument(
          "checkpoint line aggregate component provenance is malformed");
  }

  const auto& raw_result = as_object(
      object.at("result"), "checkpoint line aggregate result");
  require_exact_keys(raw_result,
      {"value", "scope", "imaginary_sign", "diagnostics"},
      "checkpoint line aggregate result");
  StoredLineIntegral result;
  result.value = parse_checkpoint_epsilon_vector(
      raw_result.at("value"), "checkpoint line aggregate epsilon vector");
  const auto scope = required_string(raw_result, "scope");
  if (scope == "stored_truncation")
    result.scope = LineIntegrationScope::StoredTruncation;
  else if (scope == "full_local_with_certified_tail")
    result.scope = LineIntegrationScope::FullLocalWithCertifiedTail;
  else
    throw std::invalid_argument(
        "checkpoint line aggregate has an unsupported integration scope");
  if (!raw_result.at("imaginary_sign").is_null())
    throw std::invalid_argument(
        "checkpoint multi-chart line aggregate cannot carry one effective rim");
  result.imaginary_sign = std::nullopt;
  const auto& epsilon = as_object(
      provenance.at("epsilon"), "checkpoint line aggregate epsilon");
  require_exact_keys(epsilon, {"min", "max"},
                     "checkpoint line aggregate epsilon");
  if (result.value.epsilon.min_power !=
          as_i32(epsilon.at("min"), "checkpoint aggregate epsilon minimum") ||
      result.value.epsilon.complete_max !=
          as_i32(epsilon.at("max"), "checkpoint aggregate epsilon maximum") ||
      required_string(provenance, "scope") != scope ||
      required_string(provenance, "error_guarantee") !=
          error_guarantee_name(result.value.error.guarantee))
    throw std::invalid_argument(
        "checkpoint line aggregate result differs from its provenance");
  if (whole_arm_epsilon_contract.has_value() &&
      result.value.epsilon.complete_max <
          whole_arm_epsilon_contract->public_required_complete_max)
    throw std::invalid_argument(
        "checkpoint whole-arm aggregate no longer covers its public required epsilon maximum");
  const auto expected_dimension = as_u32(
      local_owners.front()->summary().at("dimension"),
      "checkpoint line aggregate source dimension");
  if (result.value.dimension != expected_dimension ||
      std::any_of(local_owners.begin(), local_owners.end(),
                  [&](const auto& owner) {
                    return as_u32(owner->summary().at("dimension"),
                                  "checkpoint aggregate owner dimension") !=
                           expected_dimension;
                  }))
    throw std::invalid_argument(
        "checkpoint line aggregate dimensions disagree with its owners");
  if (result.scope == LineIntegrationScope::FullLocalWithCertifiedTail &&
      result.value.error.guarantee != ErrorGuarantee::Certified)
    throw std::invalid_argument(
        "checkpoint certified line aggregate lost its certified error envelope");

  const auto& diagnostics = as_object(
      raw_result.at("diagnostics"),
      "checkpoint line aggregate diagnostics");
  require_exact_keys(
      diagnostics,
      {"input_monomial_cells", "grouped_monomials", "zero_groups_skipped",
       "cancelled_divergent_groups", "primitive_evaluations",
       "primitive_component_applications", "primitive_component_reuses",
       "has_center_endpoint", "tail_certificate_requested",
       "tail_certificate_status", "tail_witness_radius_exact", "detail"},
      "checkpoint line aggregate diagnostics");
  result.diagnostics.input_monomial_cells = checkpoint_size_t(
      diagnostics.at("input_monomial_cells"),
      "checkpoint aggregate input monomials");
  result.diagnostics.grouped_monomials = checkpoint_size_t(
      diagnostics.at("grouped_monomials"),
      "checkpoint aggregate grouped monomials");
  result.diagnostics.zero_groups_skipped = checkpoint_size_t(
      diagnostics.at("zero_groups_skipped"),
      "checkpoint aggregate skipped zero groups");
  result.diagnostics.cancelled_divergent_groups = checkpoint_size_t(
      diagnostics.at("cancelled_divergent_groups"),
      "checkpoint aggregate cancelled divergent groups");
  result.diagnostics.primitive_evaluations = checkpoint_size_t(
      diagnostics.at("primitive_evaluations"),
      "checkpoint aggregate primitive evaluations");
  result.diagnostics.primitive_component_applications = checkpoint_size_t(
      diagnostics.at("primitive_component_applications"),
      "checkpoint aggregate primitive applications");
  result.diagnostics.primitive_component_reuses = checkpoint_size_t(
      diagnostics.at("primitive_component_reuses"),
      "checkpoint aggregate primitive reuses");
  if (!diagnostics.at("has_center_endpoint").is_bool() ||
      !diagnostics.at("tail_certificate_requested").is_bool())
    throw std::invalid_argument(
        "checkpoint line aggregate diagnostic flags must be Boolean");
  result.diagnostics.has_center_endpoint =
      diagnostics.at("has_center_endpoint").as_bool();
  result.diagnostics.tail_certificate_requested =
      diagnostics.at("tail_certificate_requested").as_bool();
  result.diagnostics.tail_certificate_status = required_string(
      diagnostics, "tail_certificate_status");
  if (!diagnostics.at("tail_witness_radius_exact").is_string())
    throw std::invalid_argument(
        "checkpoint aggregate tail witness radius must be a string");
  result.diagnostics.tail_witness_radius_exact = std::string(
      diagnostics.at("tail_witness_radius_exact").as_string());
  result.diagnostics.detail = required_string(diagnostics, "detail");

  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint line aggregate elapsed time");
  const auto& stats = as_object(object.at("runtime_stats"),
                                "checkpoint line aggregate runtime stats");
  require_exact_keys(stats, {"exports", "export_ms"},
                     "checkpoint line aggregate runtime stats");
  auto stored = std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, std::move(result),
      elapsed_ms, plan, std::move(local_owners), provenance);
  stored->restore_runtime_stats(
      as_u64(stats.at("exports"), "checkpoint aggregate exports"),
      checkpoint_nonnegative_double(stats.at("export_ms"),
                                    "checkpoint aggregate export time"));
  return stored;
}

json::array encode_strings(const std::vector<std::string>& values) {
  json::array output;
  output.reserve(values.size());
  for (const auto& value : values) output.emplace_back(value);
  return output;
}

json::object checkpoint_configuration_record(const SolverSession& session) {
  return json::object{
      {"domain", session.domain},
      {"precision_bits", session.precision_bits},
      {"output_digits", session.output_digits},
      {"symbols", encode_strings(session.symbols)},
      {"analytic", json::parse(session.analytic_identity)},
      {"chart_capacity", session.chart_capacity},
      {"local_capacity", session.local_capacity},
      {"scc_capacity", session.scc_capacity},
      {"match_capacity", session.match_capacity},
      {"endpoint_capacity", session.endpoint_capacity},
      {"transport_state_capacity", session.transport_state_capacity}};
}

std::string checkpoint_configuration_identity(const SolverSession& session) {
  return json::serialize(canonical_json_value(
      checkpoint_configuration_record(session)));
}

json::object checkpoint_chart_item(
    const SolverSession& session,
    const std::shared_ptr<PreparedChartBase>& chart) {
  const auto signature_value = json::parse(chart->signature());
  const auto& signature = as_object(
      signature_value, "prepared chart exact signature");
  const auto& analytic = as_object(
      signature.at("analytic"), "prepared chart analytic signature");
  if (json::serialize(analytic.at("session")) != session.analytic_identity)
    throw std::logic_error(
        "prepared chart session analytic identity changed after retention");
  if (required_string(signature, "identity") != chart->exact_identity())
    throw std::logic_error(
        "prepared chart exact identity changed after retention");

  json::object problem = signature;
  problem.erase("identity");
  problem.erase("analytic");
  problem.erase("scc");
  json::object request{
      {"schema", 2}, {"op", "chart.prepare"},
      {"session", session.handle}, {"key", chart->key()},
      {"identity", chart->exact_identity()},
      {"analytic", analytic.at("chart")},
      {"scc", signature.at("scc")}, {"problem", std::move(problem)}};
  return json::object{{"handle", chart->handle()},
                      {"key", chart->key()},
                      {"identity", chart->exact_identity()},
                      {"signature", chart->signature()},
                      {"request", std::move(request)}};
}

json::object checkpoint_scc_item(
    const SolverSession& session,
    const std::shared_ptr<CompositeSCCChartBase>& composite) {
  const auto signature_value = json::parse(composite->signature());
  const auto& signature = as_object(
      signature_value, "retained SCC exact signature");
  if (required_string(signature, "identity") != composite->exact_identity())
    throw std::logic_error(
        "retained SCC exact identity changed after retention");
  json::object request = signature;
  request["schema"] = 2;
  request["op"] = "scc.prepare";
  request["session"] = session.handle;
  request["key"] = composite->key();
  return json::object{{"handle", composite->handle()},
                      {"key", composite->key()},
                      {"identity", composite->exact_identity()},
                      {"signature", composite->signature()},
                      {"request", std::move(request)}};
}

json::array checkpoint_identity_manifest(const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(value, "checkpoint retained item");
    manifest.push_back(json::object{{"handle", item.at("handle")},
                                    {"key", item.at("key")},
                                    {"identity", item.at("identity")}});
  }
  return manifest;
}

json::array checkpoint_local_identity_manifest(const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(value, "checkpoint retained local");
    const auto& solution = as_object(item.at("solution"),
                                     "checkpoint local solution");
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"source_chart", item.at("source_chart")},
        {"source_operator_identity", item.at("source_operator_identity")},
        {"scalar_domain", item.at("scalar_domain")},
        {"checkpoint_identity", solution.at("checkpoint_identity")}});
  }
  return manifest;
}

json::array checkpoint_acb_match_identity_manifest(const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(value, "checkpoint retained Acb match");
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"checkpoint_identity", item.at("checkpoint_identity")},
        {"provenance_identity", item.at("provenance_identity")},
        {"exact_lattice_identity", item.at("exact_lattice_identity")}});
  }
  return manifest;
}

json::array checkpoint_exact_match_identity_manifest(
    const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(
        value, "checkpoint retained exact-rational match");
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"checkpoint_identity", item.at("checkpoint_identity")},
        {"provenance_identity", item.at("provenance_identity")},
        {"basis_sources", item.at("basis_sources")},
        {"incoming_source", item.at("incoming_source")}});
  }
  return manifest;
}

json::array checkpoint_planned_match_identity_manifest(
    const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(
        value, "checkpoint retained planned match hop");
    const auto& handoff = as_object(
        item.at("handoff"), "checkpoint planned match handoff");
    const auto& native_match = as_object(
        item.at("native_match"), "checkpoint planned embedded match");
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"checkpoint_identity", item.at("checkpoint_identity")},
        {"provenance_identity", item.at("provenance_identity")},
        {"tile_plan", handoff.at("tile_plan")},
        {"tile_plan_checkpoint_identity",
         handoff.at("tile_plan_checkpoint_identity")},
        {"tile_plan_provenance_identity",
         handoff.at("tile_plan_provenance_identity")},
        {"native_match_schema", native_match.at("schema")},
        {"native_match_checkpoint_identity",
         native_match.at("checkpoint_identity")},
        {"native_match_provenance_identity",
         native_match.at("provenance_identity")}});
  }
  return manifest;
}

json::array checkpoint_endpoint_identity_manifest(const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(value, "checkpoint retained endpoint");
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"checkpoint_identity", item.at("checkpoint_identity")},
        {"provenance_identity", item.at("provenance_identity")},
        {"source", item.at("source")}});
  }
  return manifest;
}

json::array checkpoint_tile_identity_manifest(const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(value, "checkpoint retained tile plan");
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"checkpoint_identity", item.at("checkpoint_identity")},
        {"provenance_identity", item.at("provenance_identity")}});
  }
  return manifest;
}

json::array checkpoint_transport_state_identity_manifest(
    const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(
        value, "checkpoint retained transport-arm state");
    const auto& provenance = as_object(
        item.at("provenance"),
        "checkpoint transport-arm state provenance");
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"checkpoint_identity", item.at("checkpoint_identity")},
        {"provenance_identity", item.at("provenance_identity")},
        {"tile_plan", provenance.at("tile_plan")},
        {"arm", provenance.at("arm")},
        {"anchor", provenance.at("anchor")},
        {"final_local", provenance.at("final_local")}});
  }
  return manifest;
}

json::array checkpoint_line_identity_manifest(const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(value, "checkpoint retained line result");
    const auto& source = required_string(item, "schema") ==
            "diffexp2-retained-line-aggregate-v1"
        ? as_object(item.at("provenance"),
                    "checkpoint line aggregate provenance").at("source")
        : item.at("source");
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"checkpoint_identity", item.at("checkpoint_identity")},
        {"provenance_identity", item.at("provenance_identity")},
        {"source", source}});
  }
  return manifest;
}

struct SessionCheckpointSnapshot {
  json::object header;
  json::object payload;
  std::uint64_t generation = 0;
  std::size_t charts = 0;
  std::size_t sccs = 0;
  std::size_t locals = 0;
  std::size_t exact_matches = 0;
  std::size_t acb_matches = 0;
  std::size_t planned_matches = 0;
  std::size_t endpoints = 0;
  std::size_t tile_plans = 0;
  std::size_t transport_states = 0;
  std::size_t line_results = 0;
};

SessionCheckpointSnapshot make_checkpoint_snapshot(
    SolverSession& session, const std::string& checkpoint_identity) {
  if (session.closed)
    throw std::invalid_argument("cannot checkpoint a closed solver session");
  if (session.pending_local_solves != 0 || session.pending_matches != 0 ||
      session.pending_endpoint_limits != 0 ||
      session.pending_tile_plans != 0 ||
      session.pending_transport_states != 0 ||
      session.pending_line_integrations != 0)
    throw std::invalid_argument(
        "checkpoint requires a quiescent session with no pending local solve, match, endpoint limit, tile plan, transport arm, or line integration");

  // Serialize the strong-ownership closure, not only the public registries.
  // Retained lines, tile plans, matches, and materialized locals deliberately
  // survive release of their source objects.  Serialize that immutable
  // closure while recording public registry visibility separately.
  std::unordered_map<std::string, std::shared_ptr<PreparedChartBase>>
      chart_closure = session.charts;
  std::unordered_map<std::string, std::shared_ptr<CompositeSCCChartBase>>
      scc_closure = session.sccs;
  std::unordered_map<std::string, std::shared_ptr<StoredLocalBase>>
      local_closure;
  std::unordered_map<std::string, std::shared_ptr<StoredTilePlan>>
      tile_closure;
  std::unordered_map<std::string, std::shared_ptr<StoredMatchBase>>
      match_closure;
  std::unordered_map<std::string, std::shared_ptr<StoredTransportArmState>>
      transport_closure;
  std::function<void(const std::shared_ptr<PreparedChartBase>&)> add_chart;
  std::function<void(const std::shared_ptr<StoredTilePlan>&)> add_tile;
  std::function<void(const std::shared_ptr<StoredLocalBase>&)> add_local;
  std::function<void(const std::shared_ptr<StoredMatchBase>&)> add_match;
  std::function<void(const std::shared_ptr<StoredTransportArmState>&)>
      add_transport;
  std::function<void(const std::shared_ptr<CompositeSCCChartBase>&)> add_scc;
  add_chart = [&](const std::shared_ptr<PreparedChartBase>& chart) {
    if (!chart)
      throw std::logic_error("checkpoint ownership closure contains a null chart");
    const auto [found, inserted] = chart_closure.emplace(chart->handle(), chart);
    if (!inserted && found->second.get() != chart.get())
      throw std::logic_error(
          "checkpoint ownership closure contains distinct charts with one handle");
  };
  add_scc = [&](
      const std::shared_ptr<CompositeSCCChartBase>& composite) {
    if (!composite)
      throw std::logic_error(
          "checkpoint ownership closure contains a null SCC");
    const auto [found, inserted] =
        scc_closure.emplace(composite->handle(), composite);
      if (!inserted && found->second.get() != composite.get())
        throw std::logic_error(
            "checkpoint ownership closure contains distinct SCCs with one handle");
      if (inserted)
        for (const auto& chart : composite->dependency_charts())
          add_chart(chart);
  };
  add_tile = [&](const std::shared_ptr<StoredTilePlan>& plan) {
    if (!plan)
      throw std::logic_error("checkpoint ownership closure contains a null tile plan");
    const auto [found, inserted] = tile_closure.emplace(plan->handle(), plan);
    if (!inserted) {
      if (found->second.get() != plan.get())
        throw std::logic_error(
            "checkpoint ownership closure contains distinct tile plans with one handle");
      return;
    }
    for (const auto& composite : plan->dependency_sccs()) add_scc(composite);
    for (const auto& chart : plan->dependency_charts()) add_chart(chart);
  };
  add_match = [&](const std::shared_ptr<StoredMatchBase>& match) {
    if (!match)
      throw std::logic_error("checkpoint ownership closure contains a null match");
    const auto [found, inserted] = match_closure.emplace(
        match->handle(), match);
    if (!inserted) {
      if (found->second.get() != match.get())
        throw std::logic_error(
            "checkpoint ownership closure contains distinct matches with one handle");
      return;
    }
    if (const auto exact =
            std::dynamic_pointer_cast<StoredExactRegularMatch>(match)) {
      for (const auto& local : exact->basis_owners()) add_local(local);
      add_local(exact->incoming_owner());
      return;
    }
    if (const auto hop =
            std::dynamic_pointer_cast<StoredPlannedMatchHop>(match)) {
      add_tile(hop->plan_owner());
      for (const auto& local : hop->basis_owners()) add_local(local);
      add_local(hop->incoming_owner());
      return;
    }
    if (!std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match))
      throw std::logic_error(
          "checkpoint ownership closure contains an unknown retained match implementation");
  };
  add_local = [&](const std::shared_ptr<StoredLocalBase>& local) {
    if (!local)
      throw std::logic_error("checkpoint ownership closure contains a null local");
    const auto [found, inserted] = local_closure.emplace(local->handle(), local);
    if (!inserted) {
      if (found->second.get() != local.get())
        throw std::logic_error(
            "checkpoint ownership closure contains distinct locals with one handle");
      return;
    }
    if (local->retained_derivation().has_value()) {
      const auto schema = required_string(
          *local->retained_derivation(), "schema");
      if (schema !=
              "diffexp2-retained-plan-match-local-materialization-v1" &&
          schema !=
              "diffexp2-retained-rational-row-local-application-v1")
        throw std::domain_error(
            "native checkpoint does not serialize this retained local derivation kind");
      const auto opaque = local->retained_derivation_owner();
      if (!opaque)
        throw std::logic_error(
            "checkpoint derived local lost its strong owner");
      if (schema ==
          "diffexp2-retained-plan-match-local-materialization-v1") {
        auto hop = std::static_pointer_cast<StoredPlannedMatchHop>(opaque);
        if (required_string(*local->retained_derivation(), "source_match") !=
            hop->handle())
          throw std::logic_error(
              "checkpoint materialized local derivation names a different owner handle");
        add_match(std::move(hop));
      } else {
        auto source = std::static_pointer_cast<StoredLocalBase>(opaque);
        const auto& source_record = as_object(
            local->retained_derivation()->at("source"),
            "checkpoint rational-row source");
        if (!source || required_string(source_record, "local") !=
                           source->handle())
          throw std::logic_error(
              "checkpoint rational-row derivation names a different source local");
        add_local(std::move(source));
      }
    } else if (local->retained_derivation_owner() != nullptr) {
      throw std::logic_error(
          "checkpoint primitive local unexpectedly retains a derivation owner");
    }
  };
  add_transport = [&](const std::shared_ptr<StoredTransportArmState>& state) {
    if (!state)
      throw std::logic_error(
          "checkpoint ownership closure contains a null transport-arm state");
    const auto [found, inserted] =
        transport_closure.emplace(state->handle(), state);
    if (!inserted) {
      if (found->second.get() != state.get())
        throw std::logic_error(
            "checkpoint ownership closure contains distinct transport-arm states with one handle");
      return;
    }
    add_tile(state->plan_owner());
    add_local(state->anchor_owner());
    for (const auto& basis : state->basis_owners())
      for (const auto& local : basis) add_local(local);
    for (const auto& match : state->matches()) add_match(match);
    for (const auto& local : state->tile_sources()) add_local(local);
  };
  for (const auto& [ignored, chart] : session.charts) add_chart(chart);
  for (const auto& [ignored, plan] : session.tile_plans) add_tile(plan);
  for (const auto& [ignored, local] : session.locals) add_local(local);
  for (const auto& [ignored, match] : session.matches) add_match(match);
  for (const auto& [ignored, state] : session.transport_states)
    add_transport(state);
  for (const auto& [ignored, endpoint] : session.endpoints) {
    if (!endpoint->plan_bound()) continue;
    const auto& plan = endpoint->plan_owner();
    if (!plan || !endpoint->local_owner())
      throw std::logic_error(
          "checkpoint plan-bound endpoint lost a strong plan/local owner");
    add_tile(plan);
    add_local(endpoint->local_owner());
  }
  for (const auto& [ignored, line] : session.line_results) {
    add_tile(line->plan_owner());
    for (const auto& local : line->local_owners()) add_local(local);
  }
  for (const auto& [ignored, plan] : tile_closure) {
    for (const auto& composite : plan->dependency_sccs()) add_scc(composite);
    for (const auto& chart : plan->dependency_charts()) add_chart(chart);
  }
  for (const auto& [ignored, composite] : scc_closure)
    for (const auto& chart : composite->dependency_charts()) add_chart(chart);

  std::vector<std::shared_ptr<PreparedChartBase>> charts;
  charts.reserve(chart_closure.size());
  for (const auto& [ignored, chart] : chart_closure) charts.push_back(chart);
  std::sort(charts.begin(), charts.end(), [](const auto& left,
                                             const auto& right) {
    return scoped_handle_id(left->handle(), "c:", "chart") <
           scoped_handle_id(right->handle(), "c:", "chart");
  });
  std::vector<std::shared_ptr<CompositeSCCChartBase>> sccs;
  sccs.reserve(scc_closure.size());
  for (const auto& [ignored, composite] : scc_closure)
    sccs.push_back(composite);
  std::sort(sccs.begin(), sccs.end(), [](const auto& left,
                                         const auto& right) {
    return scoped_handle_id(left->handle(), "scc:", "SCC") <
           scoped_handle_id(right->handle(), "scc:", "SCC");
  });
  std::vector<std::shared_ptr<StoredLocalBase>> locals;
  locals.reserve(local_closure.size());
  for (const auto& [ignored, local] : local_closure)
    locals.push_back(local);
  std::sort(locals.begin(), locals.end(), [](const auto& left,
                                             const auto& right) {
    return scoped_handle_id(left->handle(), "l:", "local") <
           scoped_handle_id(right->handle(), "l:", "local");
  });
  std::vector<std::shared_ptr<StoredMatchBase>> matches;
  matches.reserve(match_closure.size());
  for (const auto& [ignored, match] : match_closure)
    matches.push_back(match);
  std::sort(matches.begin(), matches.end(), [](const auto& left,
                                               const auto& right) {
    return scoped_handle_id(left->handle(), "m:", "match") <
           scoped_handle_id(right->handle(), "m:", "match");
  });
  std::vector<std::shared_ptr<StoredEndpointResult>> endpoints;
  endpoints.reserve(session.endpoints.size());
  for (const auto& [ignored, endpoint] : session.endpoints)
    endpoints.push_back(endpoint);
  std::sort(endpoints.begin(), endpoints.end(), [](const auto& left,
                                                   const auto& right) {
    return scoped_handle_id(left->handle(), "e:", "endpoint") <
           scoped_handle_id(right->handle(), "e:", "endpoint");
  });
  std::vector<std::shared_ptr<StoredTilePlan>> tile_plans;
  tile_plans.reserve(tile_closure.size());
  for (const auto& [ignored, plan] : tile_closure) tile_plans.push_back(plan);
  std::sort(tile_plans.begin(), tile_plans.end(), [](const auto& left,
                                                     const auto& right) {
    return scoped_handle_id(left->handle(), "tile:", "tile plan") <
           scoped_handle_id(right->handle(), "tile:", "tile plan");
  });
  std::vector<std::shared_ptr<StoredTransportArmState>> transport_states;
  transport_states.reserve(transport_closure.size());
  for (const auto& [ignored, state] : transport_closure)
    transport_states.push_back(state);
  std::sort(transport_states.begin(), transport_states.end(),
            [](const auto& left, const auto& right) {
    return scoped_handle_id(left->handle(), "transport:",
                            "transport-arm state") <
           scoped_handle_id(right->handle(), "transport:",
                            "transport-arm state");
  });
  std::vector<std::shared_ptr<StoredLineResult>> line_results;
  line_results.reserve(session.line_results.size());
  for (const auto& [ignored, line] : session.line_results)
    line_results.push_back(line);
  std::sort(line_results.begin(), line_results.end(), [](const auto& left,
                                                         const auto& right) {
    return scoped_handle_id(left->handle(), "line:", "line result") <
           scoped_handle_id(right->handle(), "line:", "line result");
  });

  json::array chart_items;
  chart_items.reserve(charts.size());
  std::set<std::string> chart_handles;
  for (const auto& chart : charts) {
    chart_handles.insert(chart->handle());
    chart_items.push_back(checkpoint_chart_item(session, chart));
  }
  json::array scc_items;
  scc_items.reserve(sccs.size());
  for (const auto& composite : sccs) {
    auto item = checkpoint_scc_item(session, composite);
    const auto& request = as_object(item.at("request"), "SCC request");
    for (const auto& raw_block : as_array(request.at("blocks"), "SCC blocks")) {
      const auto& block = as_object(raw_block, "SCC block");
      if (!chart_handles.contains(required_string(block, "chart")))
        throw std::invalid_argument(
            "checkpoint cannot serialize an SCC whose retained diagonal chart was publicly released");
    }
    scc_items.push_back(std::move(item));
  }
  json::array local_items;
  local_items.reserve(locals.size());
  for (const auto& local : locals)
    local_items.push_back(local->checkpoint_record());
  json::array exact_match_items;
  json::array acb_match_items;
  json::array planned_match_items;
  for (const auto& match : matches) {
    auto record = match->checkpoint_record();
    if (std::dynamic_pointer_cast<StoredPlannedMatchHop>(match))
      planned_match_items.push_back(std::move(record));
    else if (std::dynamic_pointer_cast<StoredExactRegularMatch>(match))
      exact_match_items.push_back(std::move(record));
    else if (std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match))
      acb_match_items.push_back(std::move(record));
    else
      throw std::logic_error(
          "checkpoint encountered an unknown retained match implementation");
  }
  json::array endpoint_items;
  endpoint_items.reserve(endpoints.size());
  for (const auto& endpoint : endpoints)
    endpoint_items.push_back(endpoint->checkpoint_record());
  json::array tile_items;
  tile_items.reserve(tile_plans.size());
  for (const auto& plan : tile_plans)
    tile_items.push_back(plan->checkpoint_record());
  json::array transport_items;
  transport_items.reserve(transport_states.size());
  for (const auto& state : transport_states)
    transport_items.push_back(state->checkpoint_record());
  json::array line_items;
  line_items.reserve(line_results.size());
  for (const auto& line : line_results)
    line_items.push_back(line->checkpoint_record());

  json::array visible_charts;
  for (const auto& chart : charts)
    if (session.charts.contains(chart->handle()))
      visible_charts.emplace_back(chart->handle());
  json::array visible_sccs;
  for (const auto& composite : sccs)
    if (session.sccs.contains(composite->handle()))
      visible_sccs.emplace_back(composite->handle());
  json::array visible_locals;
  for (const auto& local : locals)
    if (session.locals.contains(local->handle()))
      visible_locals.emplace_back(local->handle());
  json::array visible_matches;
  for (const auto& match : matches)
    if (session.matches.contains(match->handle()))
      visible_matches.emplace_back(match->handle());
  json::array visible_tiles;
  for (const auto& plan : tile_plans)
    if (session.tile_plans.contains(plan->handle()))
      visible_tiles.emplace_back(plan->handle());
  json::array visible_transport_states;
  for (const auto& state : transport_states)
    if (session.transport_states.contains(state->handle()))
      visible_transport_states.emplace_back(state->handle());
  json::object registry_visibility{
      {"charts", std::move(visible_charts)},
      {"sccs", std::move(visible_sccs)},
      {"locals", std::move(visible_locals)},
      {"matches", std::move(visible_matches)},
      {"tile_plans", std::move(visible_tiles)},
      {"transport_states", std::move(visible_transport_states)}};

  if (session.checkpoint_generation ==
      std::numeric_limits<std::uint64_t>::max())
    throw std::overflow_error("checkpoint generation counter overflow");
  const auto generation = session.checkpoint_generation + 1;
  auto configuration = checkpoint_configuration_record(session);
  json::object counters{
      {"next_chart", session.next_chart},
      {"next_local", session.next_local},
      {"next_scc", session.next_scc},
      {"next_match", session.next_match},
      {"next_endpoint", session.next_endpoint},
      {"next_tile_plan", session.next_tile_plan},
      {"next_transport_state", session.next_transport_state},
      {"next_line_result", session.next_line_result},
      {"total_local_solves", session.total_local_solves},
      {"total_scc_column_solves", session.total_scc_column_solves},
      {"total_local_matches", session.total_local_matches},
      {"total_endpoint_limits", session.total_endpoint_limits},
      {"total_endpoint_exports", session.total_endpoint_exports},
      {"total_tile_plans", session.total_tile_plans},
      {"total_transport_arm_marches", session.total_transport_arm_marches},
      {"total_line_integrations", session.total_line_integrations},
      {"total_line_exports", session.total_line_exports},
      {"total_local_run_parse_ms", session.total_local_run_parse_ms},
      {"total_local_kernel_ms", session.total_local_kernel_ms},
      {"total_local_match_ms", session.total_local_match_ms},
      {"total_endpoint_limit_ms", session.total_endpoint_limit_ms},
      {"total_endpoint_export_ms", session.total_endpoint_export_ms},
      {"total_tile_plan_ms", session.total_tile_plan_ms},
      {"total_transport_arm_ms", session.total_transport_arm_ms},
      {"total_line_integration_ms", session.total_line_integration_ms},
      {"total_line_export_ms", session.total_line_export_ms},
      {"checkpoint_generation", generation},
      {"checkpoint_restore_count", session.checkpoint_restore_count}};
  json::object session_record{
      {"source_handle", session.handle},
      {"configuration", configuration},
      {"configuration_identity", checkpoint_configuration_identity(session)},
      {"registry_visibility", registry_visibility},
      {"counters", std::move(counters)}};
  json::object payload{
      {"schema", kCheckpointPayloadSchema},
      {"session", std::move(session_record)},
      {"prepared_charts", chart_items},
      {"prepared_scc", scc_items},
      {"retained_locals", local_items},
      {"retained_exact_matches", exact_match_items},
      {"retained_acb_matches", acb_match_items},
      {"retained_planned_match_hops", planned_match_items},
      {"retained_endpoints", endpoint_items},
      {"retained_tile_plans", tile_items},
      {"retained_transport_states", transport_items},
      {"retained_line_results", line_items}};
  json::array mandatory_sections{"session", "prepared_charts",
                                  "prepared_scc", "retained_locals",
                                  "retained_exact_matches",
                                  "retained_acb_matches",
                                  "retained_planned_match_hops",
                                  "retained_endpoints",
                                  "retained_tile_plans",
                                  "retained_transport_states",
                                  "retained_line_results"};
  json::array deferred_kinds{"symbolic-local"};
  json::object header{
      {"format", kCheckpointFormat},
      {"schema", kCheckpointPayloadSchema},
      {"build", checkpoint::kBuildIdentity},
      {"flint", flint_version},
      {"checkpoint_identity", checkpoint_identity},
      {"configuration_identity", checkpoint_configuration_identity(session)},
      {"analytic_identity", json::parse(session.analytic_identity)},
      {"mandatory_sections", std::move(mandatory_sections)},
      {"optional_sections", json::array{}},
      {"deferred_handle_kinds", std::move(deferred_kinds)},
      {"chart_identities", checkpoint_identity_manifest(chart_items)},
      {"scc_identities", checkpoint_identity_manifest(scc_items)},
      {"local_identities", checkpoint_local_identity_manifest(local_items)},
      {"exact_match_identities",
       checkpoint_exact_match_identity_manifest(exact_match_items)},
      {"acb_match_identities",
       checkpoint_acb_match_identity_manifest(acb_match_items)},
      {"planned_match_identities",
       checkpoint_planned_match_identity_manifest(planned_match_items)},
      {"endpoint_identities",
       checkpoint_endpoint_identity_manifest(endpoint_items)},
      {"tile_plan_identities",
       checkpoint_tile_identity_manifest(tile_items)},
      {"transport_state_identities",
       checkpoint_transport_state_identity_manifest(transport_items)},
      {"line_result_identities",
       checkpoint_line_identity_manifest(line_items)},
      {"generation", generation}};
  return {std::move(header), std::move(payload), generation,
          session.charts.size(), session.sccs.size(), session.locals.size(),
          exact_match_items.size(), acb_match_items.size(),
          planned_match_items.size(), endpoints.size(),
          session.tile_plans.size(), session.transport_states.size(),
          line_results.size()};
}

std::vector<std::string> checkpoint_string_array(const json::value& raw,
                                                 const char* label) {
  std::vector<std::string> result;
  for (const auto& value : as_array(raw, label)) {
    if (!value.is_string())
      throw std::invalid_argument(std::string(label) +
                                  " must contain only strings");
    result.emplace_back(value.as_string());
  }
  return result;
}

void validate_checkpoint_envelope(const json::object& header,
                                  const json::object& payload,
                                  const std::string& expected_identity) {
  require_exact_keys(header,
      {"format", "schema", "build", "flint", "checkpoint_identity",
       "configuration_identity", "analytic_identity", "mandatory_sections",
       "optional_sections", "deferred_handle_kinds", "chart_identities",
       "scc_identities", "local_identities", "exact_match_identities",
       "acb_match_identities", "planned_match_identities",
       "endpoint_identities",
       "tile_plan_identities", "transport_state_identities",
       "line_result_identities", "generation"},
      "checkpoint header");
  require_exact_keys(payload,
      {"schema", "session", "prepared_charts", "prepared_scc",
       "retained_locals", "retained_exact_matches",
       "retained_acb_matches", "retained_planned_match_hops",
       "retained_endpoints",
       "retained_tile_plans", "retained_transport_states",
       "retained_line_results"},
      "checkpoint payload");
  if (required_string(header, "format") != kCheckpointFormat ||
      as_u32(header.at("schema"), "checkpoint header schema") !=
          kCheckpointPayloadSchema ||
      as_u32(payload.at("schema"), "checkpoint payload schema") !=
          kCheckpointPayloadSchema)
    throw std::invalid_argument("unsupported native checkpoint schema");
  if (required_string(header, "build") != checkpoint::kBuildIdentity)
    throw std::invalid_argument(
        "native checkpoint solver build identity is incompatible");
  if (required_string(header, "flint") != flint_version)
    throw std::invalid_argument(
        "native checkpoint FLINT build identity is incompatible");
  if (required_string(header, "checkpoint_identity") != expected_identity)
    throw std::invalid_argument(
        "native checkpoint identity differs from the expected identity");
  auto mandatory = checkpoint_string_array(
      header.at("mandatory_sections"), "mandatory checkpoint sections");
  std::sort(mandatory.begin(), mandatory.end());
  const std::vector<std::string> expected_sections{
      "prepared_charts", "prepared_scc", "retained_acb_matches",
      "retained_endpoints", "retained_exact_matches",
      "retained_line_results", "retained_locals",
      "retained_planned_match_hops", "retained_tile_plans",
      "retained_transport_states",
      "session"};
  if (mandatory != expected_sections)
    throw std::invalid_argument(
        "native checkpoint contains unknown or missing mandatory sections");
  if (!as_array(header.at("optional_sections"),
                "optional checkpoint sections").empty())
    throw std::invalid_argument(
        "native checkpoint declares unsupported optional sections");
  auto deferred = checkpoint_string_array(
      header.at("deferred_handle_kinds"),
      "deferred checkpoint handle kinds");
  std::sort(deferred.begin(), deferred.end());
  const std::vector<std::string> expected_deferred{"symbolic-local"};
  if (deferred != expected_deferred)
    throw std::invalid_argument(
        "native checkpoint deferred-state contract is incompatible");

  const auto& session = as_object(payload.at("session"),
                                  "checkpoint session section");
  require_exact_keys(session,
      {"source_handle", "configuration", "configuration_identity",
       "registry_visibility", "counters"}, "checkpoint session section");
  const auto& visibility = as_object(
      session.at("registry_visibility"), "checkpoint registry visibility");
  require_exact_keys(
      visibility,
      {"charts", "sccs", "locals", "matches", "tile_plans",
       "transport_states"},
      "checkpoint registry visibility");
  (void)checkpoint_string_array(visibility.at("charts"),
                                "visible checkpoint charts");
  (void)checkpoint_string_array(visibility.at("sccs"),
                                "visible checkpoint SCCs");
  (void)checkpoint_string_array(visibility.at("locals"),
                                "visible checkpoint locals");
  (void)checkpoint_string_array(visibility.at("matches"),
                                "visible checkpoint matches");
  (void)checkpoint_string_array(visibility.at("tile_plans"),
                                "visible checkpoint tile plans");
  (void)checkpoint_string_array(visibility.at("transport_states"),
                                "visible checkpoint transport states");
  if (required_string(session, "configuration_identity") !=
      required_string(header, "configuration_identity"))
    throw std::invalid_argument(
        "checkpoint header and payload configuration identities disagree");
  const auto& configuration = as_object(
      session.at("configuration"), "checkpoint configuration");
  if (json::serialize(canonical_json_value(configuration)) !=
      required_string(header, "configuration_identity"))
    throw std::invalid_argument(
        "checkpoint configuration does not reproduce its exact identity");
  if (configuration.at("analytic") != header.at("analytic_identity"))
    throw std::invalid_argument(
        "checkpoint analytic-regularization identity is inconsistent");
  const auto& chart_items = as_array(payload.at("prepared_charts"),
                                     "checkpoint prepared charts");
  const auto& scc_items = as_array(payload.at("prepared_scc"),
                                   "checkpoint prepared SCC charts");
  const auto& local_items = as_array(payload.at("retained_locals"),
                                     "checkpoint retained locals");
  const auto& exact_match_items = as_array(
      payload.at("retained_exact_matches"),
      "checkpoint retained exact-rational matches");
  const auto& acb_match_items = as_array(
      payload.at("retained_acb_matches"),
      "checkpoint retained Acb matches");
  const auto& planned_match_items = as_array(
      payload.at("retained_planned_match_hops"),
      "checkpoint retained planned match hops");
  const auto& endpoint_items = as_array(
      payload.at("retained_endpoints"),
      "checkpoint retained endpoints");
  const auto& tile_items = as_array(
      payload.at("retained_tile_plans"),
      "checkpoint retained tile plans");
  const auto& transport_items = as_array(
      payload.at("retained_transport_states"),
      "checkpoint retained transport-arm states");
  const auto& line_items = as_array(
      payload.at("retained_line_results"),
      "checkpoint retained line results");
  if (checkpoint_identity_manifest(chart_items) !=
          as_array(header.at("chart_identities"),
                   "checkpoint chart identities") ||
      checkpoint_identity_manifest(scc_items) !=
          as_array(header.at("scc_identities"),
                   "checkpoint SCC identities") ||
      checkpoint_local_identity_manifest(local_items) !=
          as_array(header.at("local_identities"),
                   "checkpoint local identities") ||
      checkpoint_exact_match_identity_manifest(exact_match_items) !=
          as_array(header.at("exact_match_identities"),
                   "checkpoint exact-match identities") ||
      checkpoint_acb_match_identity_manifest(acb_match_items) !=
          as_array(header.at("acb_match_identities"),
                   "checkpoint Acb match identities") ||
      checkpoint_planned_match_identity_manifest(planned_match_items) !=
          as_array(header.at("planned_match_identities"),
                   "checkpoint planned-match identities") ||
      checkpoint_endpoint_identity_manifest(endpoint_items) !=
          as_array(header.at("endpoint_identities"),
                   "checkpoint endpoint identities") ||
      checkpoint_tile_identity_manifest(tile_items) !=
          as_array(header.at("tile_plan_identities"),
                   "checkpoint tile-plan identities") ||
      checkpoint_transport_state_identity_manifest(transport_items) !=
          as_array(header.at("transport_state_identities"),
                   "checkpoint transport-state identities") ||
      checkpoint_line_identity_manifest(line_items) !=
          as_array(header.at("line_result_identities"),
                   "checkpoint line-result identities"))
    throw std::invalid_argument(
        "checkpoint retained identity manifest is inconsistent");
  if (as_u64(header.at("generation"), "checkpoint generation") !=
      as_u64(as_object(session.at("counters"), "checkpoint counters")
                 .at("checkpoint_generation"),
             "checkpoint generation"))
    throw std::invalid_argument(
        "checkpoint generation differs between header and payload");
}

json::object restore_checkpoint(const std::string& path,
                                const std::string& expected_identity) {
  const auto container = checkpoint::read(path);
  const auto header_value = json::parse(container.header_json);
  const auto payload_value = json::parse(container.payload_json);
  const auto& header = as_object(header_value, "checkpoint JSON header");
  const auto& payload = as_object(payload_value, "checkpoint JSON payload");
  validate_checkpoint_envelope(header, payload, expected_identity);

  const auto& saved_session = as_object(payload.at("session"),
                                        "checkpoint session section");
  const auto& configuration = as_object(
      saved_session.at("configuration"), "checkpoint configuration");
  require_exact_keys(configuration,
      {"domain", "precision_bits", "output_digits", "symbols", "analytic",
       "chart_capacity", "local_capacity", "scc_capacity",
       "match_capacity", "endpoint_capacity", "transport_state_capacity"},
      "checkpoint configuration");
  const auto& raw_visibility = as_object(
      saved_session.at("registry_visibility"),
      "checkpoint registry visibility");
  const auto visibility_set = [&](const char* key, std::string_view prefix) {
    std::set<std::string> result;
    for (const auto& handle : checkpoint_string_array(
             raw_visibility.at(key), key)) {
      (void)scoped_handle_id(handle, prefix, key);
      if (!result.insert(handle).second)
        throw std::invalid_argument(
            std::string("checkpoint registry visibility duplicates a ") + key +
            " handle");
    }
    return result;
  };
  const auto visible_charts = visibility_set("charts", "c:");
  const auto visible_sccs = visibility_set("sccs", "scc:");
  const auto visible_locals = visibility_set("locals", "l:");
  const auto visible_matches = visibility_set("matches", "m:");
  const auto visible_tiles = visibility_set("tile_plans", "tile:");
  const auto visible_transport_states =
      visibility_set("transport_states", "transport:");
  const auto configured_chart_capacity = static_cast<std::size_t>(as_u64(
      configuration.at("chart_capacity"),
      "checkpoint configured chart capacity"));
  if (visible_charts.size() > configured_chart_capacity)
    throw std::invalid_argument(
        "checkpoint visible charts exceed the restored session capacity");
  const auto configured_scc_capacity = static_cast<std::size_t>(as_u64(
      configuration.at("scc_capacity"),
      "checkpoint configured SCC capacity"));
  if (visible_sccs.size() > configured_scc_capacity)
    throw std::invalid_argument(
        "checkpoint visible SCCs exceed the restored session capacity");
  if (visible_matches.size() > static_cast<std::size_t>(as_u64(
          configuration.at("match_capacity"),
          "checkpoint configured match capacity")))
    throw std::invalid_argument(
        "checkpoint visible matches exceed the restored session capacity");
  const auto configured_transport_state_capacity =
      static_cast<std::size_t>(as_u64(
          configuration.at("transport_state_capacity"),
          "checkpoint configured transport-state capacity"));
  if (visible_transport_states.size() > configured_transport_state_capacity)
    throw std::invalid_argument(
        "checkpoint visible transport states exceed the restored session capacity");
  json::object create{
      {"schema", 2}, {"op", "session.create"},
      {"domain", configuration.at("domain")},
      {"precision_bits", configuration.at("precision_bits")},
      {"output_digits", configuration.at("output_digits")},
      {"symbols", configuration.at("symbols")},
      {"analytic", configuration.at("analytic")},
      {"chart_capacity", configuration.at("chart_capacity")},
      {"local_capacity", configuration.at("local_capacity")},
      {"scc_capacity", configuration.at("scc_capacity")},
      {"match_capacity", configuration.at("match_capacity")},
      {"endpoint_capacity", configuration.at("endpoint_capacity")},
      {"transport_state_capacity",
       configuration.at("transport_state_capacity")}};
  const auto created = run_session_command(create);
  const auto restored_handle = required_string(created, "session");
  bool live = true;
  try {
    const auto restored = find_session(restored_handle);
    {
      // Dependency-only chart/SCC owners are replayed briefly so tile plans
      // can acquire their strong pointers, then removed from the public maps.
      // Their closure may legitimately exceed the public capacities.
      const auto chart_closure_size = as_array(
          payload.at("prepared_charts"),
          "checkpoint prepared chart closure").size();
      const auto scc_closure_size = as_array(
          payload.at("prepared_scc"),
          "checkpoint prepared SCC closure").size();
      std::lock_guard<std::mutex> lock(restored->mutex);
      restored->chart_capacity =
          std::max(configured_chart_capacity, chart_closure_size);
      restored->scc_capacity =
          std::max(configured_scc_capacity, scc_closure_size);
    }
    const auto source_handle = required_string(saved_session, "source_handle");
    const auto source_configuration_identity = required_string(
        saved_session, "configuration_identity");
    json::array restored_charts;
    std::set<std::string> all_chart_handles;
    std::uint64_t largest_chart = 0;
    for (const auto& raw_item : as_array(
             payload.at("prepared_charts"), "checkpoint prepared charts")) {
      const auto& item = as_object(raw_item, "checkpoint chart item");
      require_exact_keys(item,
          {"handle", "key", "identity", "signature", "request"},
          "checkpoint chart item");
      const auto old_handle = required_string(item, "handle");
      if (!all_chart_handles.insert(old_handle).second)
        throw std::invalid_argument(
            "checkpoint contains duplicate prepared chart handles");
      const auto handle_id = scoped_handle_id(old_handle, "c:", "chart");
      if (handle_id <= largest_chart)
        throw std::invalid_argument(
            "checkpoint chart handles are not in strict creation order");
      largest_chart = handle_id;
      auto request = as_object(item.at("request"),
                               "checkpoint chart request");
      if (required_string(request, "session") != source_handle ||
          required_string(request, "key") != required_string(item, "key") ||
          required_string(request, "identity") !=
              required_string(item, "identity"))
        throw std::invalid_argument(
            "checkpoint chart request provenance is inconsistent");
      request["session"] = restored_handle;
      {
        std::lock_guard<std::mutex> lock(restored->mutex);
        restored->next_chart = handle_id;
      }
      const auto result = run_session_command(request);
      if (required_string(result, "chart") != old_handle)
        throw std::logic_error(
            "checkpoint chart handle could not be restored exactly");
      {
        std::lock_guard<std::mutex> lock(restored->mutex);
        const auto found = restored->charts.find(old_handle);
        if (found == restored->charts.end() ||
            found->second->signature() != required_string(item, "signature"))
          throw std::invalid_argument(
              "restored chart does not reproduce its exact operator identity");
      }
      if (visible_charts.contains(old_handle))
        restored_charts.push_back(json::object{
            {"chart", old_handle}, {"key", item.at("key")},
            {"identity", item.at("identity")}});
    }
    if (!std::includes(all_chart_handles.begin(), all_chart_handles.end(),
                       visible_charts.begin(), visible_charts.end()))
      throw std::invalid_argument(
          "checkpoint chart visibility names an absent ownership object");

    json::array restored_sccs;
    std::set<std::string> all_scc_handles;
    std::uint64_t largest_scc = 0;
    for (const auto& raw_item : as_array(
             payload.at("prepared_scc"), "checkpoint prepared SCC charts")) {
      const auto& item = as_object(raw_item, "checkpoint SCC item");
      require_exact_keys(item,
          {"handle", "key", "identity", "signature", "request"},
          "checkpoint SCC item");
      const auto old_handle = required_string(item, "handle");
      if (!all_scc_handles.insert(old_handle).second)
        throw std::invalid_argument(
            "checkpoint contains duplicate retained SCC handles");
      const auto handle_id = scoped_handle_id(old_handle, "scc:", "SCC");
      if (handle_id <= largest_scc)
        throw std::invalid_argument(
            "checkpoint SCC handles are not in strict creation order");
      largest_scc = handle_id;
      auto request = as_object(item.at("request"),
                               "checkpoint SCC request");
      if (required_string(request, "session") != source_handle ||
          required_string(request, "key") != required_string(item, "key") ||
          required_string(request, "identity") !=
              required_string(item, "identity"))
        throw std::invalid_argument(
            "checkpoint SCC request provenance is inconsistent");
      request["session"] = restored_handle;
      {
        std::lock_guard<std::mutex> lock(restored->mutex);
        restored->next_scc = handle_id;
      }
      const auto result = run_session_command(request);
      if (required_string(result, "scc") != old_handle)
        throw std::logic_error(
            "checkpoint SCC handle could not be restored exactly");
      {
        std::lock_guard<std::mutex> lock(restored->mutex);
        const auto found = restored->sccs.find(old_handle);
        if (found == restored->sccs.end() ||
            found->second->signature() != required_string(item, "signature"))
          throw std::invalid_argument(
              "restored SCC does not reproduce its exact graph/operator identity");
      }
      if (visible_sccs.contains(old_handle))
        restored_sccs.push_back(json::object{
            {"scc", old_handle}, {"key", item.at("key")},
            {"identity", item.at("identity")}});
    }
    if (!std::includes(all_scc_handles.begin(), all_scc_handles.end(),
                       visible_sccs.begin(), visible_sccs.end()))
      throw std::invalid_argument(
          "checkpoint SCC visibility names an absent ownership object");

    std::unique_ptr<AcbPrecisionLease> checkpoint_acb_lease;
    if (restored->domain == "acb") {
      checkpoint_acb_lease =
          std::make_unique<AcbPrecisionLease>(restored->precision_bits);
      ComplexBall::set_precision(restored->precision_bits);
    }

    json::array restored_locals;
    json::array restored_exact_matches;
    json::array restored_acb_matches;
    json::array restored_planned_matches;
    json::array restored_transport_states;
    json::array restored_endpoints;
    json::array restored_tile_plans;
    json::array restored_line_results;
    std::uint64_t largest_local = 0, largest_match = 0;
    std::uint64_t largest_endpoint = 0, largest_tile_plan = 0;
    std::uint64_t largest_transport_state = 0;
    std::uint64_t largest_line_result = 0;
    {
      // Rebuild the ownership DAG under one publication lock. Primitive
      // locals and plans are roots; exact/planned matches and materialized
      // locals are admitted by a fixed-point walk over their strong-owner
      // lineage. No partially restored chain can become observable.
      std::lock_guard<std::mutex> restore_state_lock(restored->mutex);
      const auto& saved_locals = as_array(
          payload.at("retained_locals"), "checkpoint retained locals");
      const auto& saved_tiles = as_array(
          payload.at("retained_tile_plans"),
          "checkpoint retained tile plans");
      const auto& saved_exact = as_array(
          payload.at("retained_exact_matches"),
          "checkpoint retained exact-rational matches");
      const auto& saved_acb = as_array(
          payload.at("retained_acb_matches"),
          "checkpoint retained Acb matches");
      const auto& saved_planned = as_array(
          payload.at("retained_planned_match_hops"),
          "checkpoint retained planned match hops");
      const auto& saved_transport_states = as_array(
          payload.at("retained_transport_states"),
          "checkpoint retained transport-arm states");
      if (visible_locals.size() > restored->local_capacity)
        throw std::invalid_argument(
            "checkpoint visible locals exceed the restored session capacity");
      if (visible_tiles.size() > restored->tile_plan_capacity)
        throw std::invalid_argument(
            "checkpoint visible tile plans exceed the restored session capacity");
      if (visible_transport_states.size() >
          restored->transport_state_capacity)
        throw std::invalid_argument(
            "checkpoint visible transport states exceed the restored session capacity");

      std::set<std::string> all_tile_handles;
      for (const auto& raw_item : saved_tiles) {
        const auto& item = as_object(raw_item,
                                     "checkpoint retained tile plan");
        const auto handle = required_string(item, "handle");
        if (!all_tile_handles.insert(handle).second)
          throw std::invalid_argument(
              "checkpoint contains duplicate retained tile-plan handles");
        const auto id = scoped_handle_id(handle, "tile:", "tile plan");
        if (id <= largest_tile_plan)
          throw std::invalid_argument(
              "checkpoint tile-plan handles are not in strict creation order");
        largest_tile_plan = id;
        auto plan = restore_checkpoint_tile_plan_record(
            item, restored->charts, restored->sccs);
        if (plan->checkpoint_record() != raw_item ||
            !restored->tile_plans.emplace(handle, plan).second)
          throw std::invalid_argument(
              "restored tile plan does not reproduce its exact retained state");
        if (visible_tiles.contains(handle))
          restored_tile_plans.push_back(json::object{
              {"tile_plan", handle},
              {"checkpoint_identity", item.at("checkpoint_identity")},
              {"provenance_identity", item.at("provenance_identity")},
              {"generation", header.at("generation")}});
      }
      if (!std::includes(all_tile_handles.begin(), all_tile_handles.end(),
                         visible_tiles.begin(), visible_tiles.end()))
        throw std::invalid_argument(
            "checkpoint tile visibility names an absent ownership object");

      std::set<std::string> all_local_handles;
      std::vector<const json::value*> pending_locals;
      auto validate_local_source = [&](const std::shared_ptr<StoredLocalBase>& local) {
        const auto& source = local->source_chart();
        const bool rational_row_derived =
            local->retained_derivation().has_value() &&
            required_string(*local->retained_derivation(), "schema") ==
                "diffexp2-retained-rational-row-local-application-v1";
        if (source.starts_with("c:")) {
          (void)scoped_handle_id(source, "c:", "local source chart");
          const auto found = restored->charts.find(source);
          if (!rational_row_derived && found != restored->charts.end() &&
              found->second->exact_identity() !=
                  local->source_operator_identity())
            throw std::invalid_argument(
                "checkpoint local source identity disagrees with its restored chart");
        } else if (source.starts_with("scc:")) {
          (void)scoped_handle_id(source, "scc:", "local source SCC");
          const auto found = restored->sccs.find(source);
          if (!rational_row_derived && found != restored->sccs.end() &&
              found->second->exact_identity() !=
                  local->source_operator_identity())
            throw std::invalid_argument(
                "checkpoint local source identity disagrees with its restored SCC");
        } else {
          throw std::invalid_argument(
              "checkpoint local source is neither a chart nor an SCC handle");
        }
        if (local->column_provenance().has_value()) {
          const auto& column = *local->column_provenance();
          const auto found = restored->sccs.find(column.scc_handle);
          if (column.scc_handle != local->source_chart() ||
              column.scc_exact_identity != local->source_operator_identity() ||
              (found != restored->sccs.end() &&
               found->second->exact_identity() != column.scc_exact_identity))
            throw std::invalid_argument(
                "checkpoint local source and SCC-column provenance disagree");
        }
      };
      auto install_local = [&](const json::value& raw_item,
                               std::shared_ptr<void> owner) {
        const auto& item = as_object(raw_item, "checkpoint retained local");
        const auto handle = required_string(item, "handle");
        std::shared_ptr<StoredLocalBase> local;
        if (restored->domain == "rational")
          local = restore_checkpoint_local_record<Rational>(
              item, restored->domain, restored->precision_bits, owner);
        else if (restored->domain == "acb")
          local = restore_checkpoint_local_record<ComplexBall>(
              item, restored->domain, restored->precision_bits, owner);
        else
          throw std::invalid_argument(
              "native checkpoint cannot restore symbolic local state");
        if (local->checkpoint_record() != raw_item)
          throw std::invalid_argument(
              "restored local does not reproduce its exact retained state");
        validate_local_source(local);
        if (!restored->locals.emplace(handle, local).second)
          throw std::invalid_argument(
              "checkpoint contains duplicate retained local handles");
        if (visible_locals.contains(handle))
          restored_locals.push_back(json::object{
              {"local", handle}, {"chart", local->source_chart()},
              {"source_operator_identity", local->source_operator_identity()},
              {"checkpoint_identity", local->checkpoint_identity()},
              {"generation", header.at("generation")}});
        return local;
      };
      for (const auto& raw_item : saved_locals) {
        const auto& item = as_object(raw_item, "checkpoint retained local");
        const auto handle = required_string(item, "handle");
        if (!all_local_handles.insert(handle).second)
          throw std::invalid_argument(
              "checkpoint contains duplicate retained local handles");
        const auto id = scoped_handle_id(handle, "l:", "local");
        if (id <= largest_local)
          throw std::invalid_argument(
              "checkpoint local handles are not in strict creation order");
        largest_local = id;
        if (item.at("retained_derivation").is_null())
          (void)install_local(raw_item, nullptr);
        else
          pending_locals.push_back(&raw_item);
      }
      if (!std::includes(all_local_handles.begin(), all_local_handles.end(),
                         visible_locals.begin(), visible_locals.end()))
        throw std::invalid_argument(
            "checkpoint local visibility names an absent ownership object");

      std::set<std::string> all_match_handles;
      const auto register_match_handles = [&](const json::array& records,
                                               const char* label) {
        std::uint64_t previous = 0;
        for (const auto& raw_item : records) {
          const auto& item = as_object(raw_item, label);
          const auto handle = required_string(item, "handle");
          if (!all_match_handles.insert(handle).second)
            throw std::invalid_argument(
                "checkpoint contains duplicate retained match handles");
          const auto id = scoped_handle_id(handle, "m:", "match");
          if (id <= previous)
            throw std::invalid_argument(
                "checkpoint match handles are not in strict creation order");
          previous = id;
          largest_match = std::max(largest_match, id);
        }
      };
      register_match_handles(saved_exact, "checkpoint exact match");
      register_match_handles(saved_acb, "checkpoint Acb match");
      register_match_handles(saved_planned, "checkpoint planned match");
      if (!std::includes(all_match_handles.begin(), all_match_handles.end(),
                         visible_matches.begin(), visible_matches.end()))
        throw std::invalid_argument(
            "checkpoint match visibility names an absent ownership object");

      auto publish_match = [&](const json::value& raw_item,
                               const std::shared_ptr<StoredMatchBase>& match,
                               json::array& response) {
        const auto& item = as_object(raw_item, "checkpoint retained match");
        const auto handle = required_string(item, "handle");
        if (match->checkpoint_record() != raw_item ||
            !restored->matches.emplace(handle, match).second)
          throw std::invalid_argument(
              "restored match does not reproduce its exact retained state");
        if (visible_matches.contains(handle))
          response.push_back(json::object{
              {"match", handle},
              {"checkpoint_identity", item.at("checkpoint_identity")},
              {"provenance_identity", item.at("provenance_identity")},
              {"generation", header.at("generation")}});
      };

      if (!saved_acb.empty() && restored->domain != "acb")
        throw std::invalid_argument(
            "retained Acb match state requires an Acb checkpoint session");
      for (const auto& raw_item : saved_acb) {
        auto match = restore_checkpoint_acb_match_record(
            raw_item, source_configuration_identity);
        publish_match(raw_item, match, restored_acb_matches);
      }
      std::vector<const json::value*> pending_exact;
      std::vector<const json::value*> pending_planned;
      for (const auto& raw_item : saved_exact) pending_exact.push_back(&raw_item);
      for (const auto& raw_item : saved_planned) pending_planned.push_back(&raw_item);

      auto resolve_local = [&](const std::string& handle) {
        const auto found = restored->locals.find(handle);
        return found == restored->locals.end()
            ? std::shared_ptr<StoredLocalBase>() : found->second;
      };
      std::size_t remaining = pending_exact.size() + pending_planned.size() +
                              pending_locals.size();
      while (remaining != 0) {
        bool progress = false;
        for (auto*& raw_ptr : pending_exact) {
          if (raw_ptr == nullptr) continue;
          const auto& item = as_object(*raw_ptr, "checkpoint exact match");
          std::vector<std::shared_ptr<StoredLocalBase>> basis;
          bool ready = true;
          for (const auto& raw_source : as_array(
                   item.at("basis_sources"), "checkpoint exact basis")) {
            auto local = resolve_local(required_string(
                as_object(raw_source, "checkpoint exact basis source"),
                "local"));
            if (!local) { ready = false; break; }
            basis.push_back(std::move(local));
          }
          auto incoming = resolve_local(required_string(
              as_object(item.at("incoming_source"),
                        "checkpoint exact incoming source"), "local"));
          if (!ready || !incoming) continue;
          auto match = restore_checkpoint_exact_match_record(
              item, std::move(basis), incoming);
          publish_match(*raw_ptr, match, restored_exact_matches);
          raw_ptr = nullptr; --remaining; progress = true;
        }
        for (auto*& raw_ptr : pending_planned) {
          if (raw_ptr == nullptr) continue;
          const auto& item = as_object(*raw_ptr, "checkpoint planned match");
          const auto& handoff = as_object(item.at("handoff"),
                                          "checkpoint planned handoff");
          const auto plan_found = restored->tile_plans.find(
              required_string(handoff, "tile_plan"));
          if (plan_found == restored->tile_plans.end()) continue;
          const auto& producing = as_object(handoff.at("producing"),
                                             "checkpoint producing handoff");
          const auto& incoming_source = as_object(
              producing.at("incoming"), "checkpoint planned incoming");
          auto incoming = resolve_local(required_string(incoming_source, "local"));
          if (!incoming) continue;
          std::vector<std::shared_ptr<StoredLocalBase>> basis;
          bool ready = true;
          const auto& receiving = as_object(handoff.at("receiving"),
                                             "checkpoint receiving handoff");
          for (const auto& raw_source : as_array(
                   receiving.at("basis"), "checkpoint planned basis")) {
            auto local = resolve_local(required_string(
                as_object(raw_source, "checkpoint planned basis source"),
                "local"));
            if (!local) { ready = false; break; }
            basis.push_back(std::move(local));
          }
          if (!ready) continue;
          auto match = restore_checkpoint_planned_match_hop_record(
              item, plan_found->second, std::move(basis), incoming,
              source_configuration_identity);
          publish_match(*raw_ptr, match, restored_planned_matches);
          raw_ptr = nullptr; --remaining; progress = true;
        }
        for (auto*& raw_ptr : pending_locals) {
          if (raw_ptr == nullptr) continue;
          const auto& item = as_object(*raw_ptr,
                                       "checkpoint derived local");
          const auto& derivation = as_object(
              item.at("retained_derivation"),
              "checkpoint retained-local derivation");
          const auto derivation_schema = required_string(
              derivation, "schema");
          const auto& lineage = as_object(item.at("retained_owner_lineage"),
                                           "checkpoint local owner lineage");
          if (derivation_schema ==
              "diffexp2-retained-rational-row-local-application-v1") {
            const auto source_handle = required_string(
                lineage, "source_local");
            const auto source = restored->locals.find(source_handle);
            if (source == restored->locals.end()) continue;
            (void)install_local(
                *raw_ptr,
                std::static_pointer_cast<void>(source->second));
            raw_ptr = nullptr; --remaining; progress = true;
            continue;
          }
          if (derivation_schema !=
              "diffexp2-retained-plan-match-local-materialization-v1")
            throw std::invalid_argument(
                "checkpoint retained local has an unsupported derivation kind");
          const auto found = restored->matches.find(
              required_string(lineage, "match"));
          if (found == restored->matches.end()) continue;
          auto hop = std::dynamic_pointer_cast<StoredPlannedMatchHop>(
              found->second);
          if (!hop)
            throw std::invalid_argument(
                "checkpoint materialized local owner is not a planned match hop");
          auto local = install_local(
              *raw_ptr, std::static_pointer_cast<void>(hop));
          hop->validate_materialized_derivation(
              *local->retained_derivation(), local->scalar_domain());
          raw_ptr = nullptr; --remaining; progress = true;
        }
        if (!progress)
          throw std::invalid_argument(
              "checkpoint local/match ownership graph has a missing dependency or cycle");
      }

      std::set<std::string> all_transport_state_handles;
      for (const auto& raw_item : saved_transport_states) {
        const auto& item = as_object(
            raw_item, "checkpoint retained transport-arm state");
        const auto handle = required_string(item, "handle");
        if (!all_transport_state_handles.insert(handle).second)
          throw std::invalid_argument(
              "checkpoint contains duplicate retained transport-state handles");
        const auto id = scoped_handle_id(
            handle, "transport:", "transport-arm state");
        if (id <= largest_transport_state)
          throw std::invalid_argument(
              "checkpoint transport-state handles are not in strict creation order");
        largest_transport_state = id;
        auto state = restore_checkpoint_transport_arm_state_record(
            raw_item, restored->tile_plans, restored->locals,
            restored->matches);
        if (state->checkpoint_record() != raw_item ||
            !restored->transport_states.emplace(handle, state).second)
          throw std::invalid_argument(
              "restored transport-arm state does not reproduce its exact retained state");
        if (visible_transport_states.contains(handle))
          restored_transport_states.push_back(json::object{
              {"transport_state", handle},
              {"checkpoint_identity", item.at("checkpoint_identity")},
              {"provenance_identity", item.at("provenance_identity")},
              {"generation", header.at("generation")}});
      }
      if (!std::includes(all_transport_state_handles.begin(),
                         all_transport_state_handles.end(),
                         visible_transport_states.begin(),
                         visible_transport_states.end()))
        throw std::invalid_argument(
            "checkpoint transport-state visibility names an absent ownership object");

      const auto& saved_endpoints = as_array(
          payload.at("retained_endpoints"),
          "checkpoint retained endpoints");
      if (saved_endpoints.size() > restored->endpoint_capacity)
        throw std::invalid_argument(
            "checkpoint retained endpoints exceed the restored session capacity");
      for (const auto& raw_item : saved_endpoints) {
        const auto& item = as_object(raw_item, "checkpoint retained endpoint");
        const auto handle = required_string(item, "handle");
        const auto id = scoped_handle_id(handle, "e:", "endpoint");
        if (id <= largest_endpoint)
          throw std::invalid_argument(
              "checkpoint endpoint handles are not in strict creation order");
        largest_endpoint = id;
        const auto& source = as_object(item.at("source"),
                                       "checkpoint endpoint source");
        std::shared_ptr<StoredEndpointResult> endpoint;
        const auto endpoint_schema = required_string(item, "schema");
        if (endpoint_schema ==
            "diffexp2-retained-plan-bound-endpoint-result-v1") {
          const auto plan_found = restored->tile_plans.find(
              required_string(source, "tile_plan"));
          const auto local_found = restored->locals.find(
              required_string(source, "local"));
          if (plan_found == restored->tile_plans.end() ||
              local_found == restored->locals.end())
            throw std::invalid_argument(
                "checkpoint plan-bound endpoint lost a strongly owned plan or local");
          endpoint = restore_checkpoint_planned_endpoint_record(
              item, restored->domain, plan_found->second,
              local_found->second);
        } else {
          endpoint = restore_checkpoint_endpoint_record(
              item, restored->domain);
          const auto found = restored->locals.find(
              required_string(source, "local"));
          if (found != restored->locals.end() &&
              (found->second->source_chart() !=
                   required_string(source, "chart") ||
               found->second->source_operator_identity() !=
                   required_string(source, "source_operator_identity") ||
               found->second->checkpoint_identity() !=
                   required_string(source, "checkpoint_identity") ||
               found->second->exact_analytic_metadata() !=
                   item.at("analytic_metadata")))
            throw std::invalid_argument(
                "checkpoint endpoint provenance disagrees with its restored source local");
        }
        if (endpoint->checkpoint_record() != raw_item ||
            !restored->endpoints.emplace(handle, endpoint).second)
          throw std::invalid_argument(
              "restored endpoint does not reproduce its exact retained state");
        restored_endpoints.push_back(json::object{
            {"endpoint", handle},
            {"checkpoint_identity", item.at("checkpoint_identity")},
            {"provenance_identity", item.at("provenance_identity")},
            {"generation", header.at("generation")}});
      }

      const auto& saved_lines = as_array(
          payload.at("retained_line_results"),
          "checkpoint retained line results");
      if (saved_lines.size() > restored->line_result_capacity)
        throw std::invalid_argument(
            "checkpoint retained line results exceed the restored session capacity");
      for (const auto& raw_item : saved_lines) {
        const auto& item = as_object(raw_item, "checkpoint retained line result");
        const auto handle = required_string(item, "handle");
        const auto id = scoped_handle_id(handle, "line:", "line result");
        if (id <= largest_line_result)
          throw std::invalid_argument(
              "checkpoint line-result handles are not in strict creation order");
        largest_line_result = id;
        const auto schema = required_string(item, "schema");
        const json::object* source = nullptr;
        if (schema == "diffexp2-retained-line-aggregate-v1")
          source = &as_object(
              as_object(item.at("provenance"),
                        "checkpoint line aggregate provenance").at("source"),
              "checkpoint line aggregate source");
        else
          source = &as_object(item.at("source"), "checkpoint line source");
        const auto plan = restored->tile_plans.find(
            required_string(*source, "tile_plan"));
        if (plan == restored->tile_plans.end())
          throw std::invalid_argument(
              "checkpoint line result lost its strongly owned plan");
        std::shared_ptr<StoredLineResult> line;
        if (schema == "diffexp2-retained-line-aggregate-v1") {
          std::vector<std::shared_ptr<StoredLocalBase>> owners;
          for (const auto& raw_owner : as_array(
                   source->at("locals"),
                   "checkpoint line aggregate local owners")) {
            const auto& owner = as_object(
                raw_owner, "checkpoint line aggregate local owner");
            const auto local = restored->locals.find(
                required_string(owner, "local"));
            if (local == restored->locals.end())
              throw std::invalid_argument(
                  "checkpoint line aggregate lost a strongly owned local");
            owners.push_back(local->second);
          }
          line = restore_checkpoint_line_aggregate_record(
              item, plan->second, std::move(owners));
        } else {
          const auto local = restored->locals.find(
              required_string(*source, "local"));
          if (local == restored->locals.end())
            throw std::invalid_argument(
                "checkpoint line result lost its strongly owned local");
          line = restore_checkpoint_line_result_record(
              item, plan->second, local->second);
        }
        if (line->checkpoint_record() != raw_item ||
            !restored->line_results.emplace(handle, line).second)
          throw std::invalid_argument(
              "restored line result does not reproduce its exact retained state");
        restored_line_results.push_back(json::object{
            {"line", handle},
            {"checkpoint_identity", item.at("checkpoint_identity")},
            {"provenance_identity", item.at("provenance_identity")},
            {"generation", header.at("generation")}});
      }

      for (auto it = restored->transport_states.begin();
           it != restored->transport_states.end();)
        it = visible_transport_states.contains(it->first)
            ? std::next(it) : restored->transport_states.erase(it);
      for (auto it = restored->matches.begin(); it != restored->matches.end();)
        it = visible_matches.contains(it->first)
            ? std::next(it) : restored->matches.erase(it);
      for (auto it = restored->tile_plans.begin(); it != restored->tile_plans.end();)
        it = visible_tiles.contains(it->first)
            ? std::next(it) : restored->tile_plans.erase(it);
      for (auto it = restored->locals.begin(); it != restored->locals.end();)
        it = visible_locals.contains(it->first)
            ? std::next(it) : restored->locals.erase(it);
      for (auto iterator = restored->sccs.begin();
           iterator != restored->sccs.end();) {
        if (visible_sccs.contains(iterator->first)) {
          ++iterator;
        } else {
          restored->scc_handles_by_key.erase(iterator->second->key());
          iterator = restored->sccs.erase(iterator);
        }
      }
      for (auto it = restored->charts.begin(); it != restored->charts.end();) {
        if (visible_charts.contains(it->first)) {
          ++it;
        } else {
          restored->handles_by_key.erase(it->second->key());
          it = restored->charts.erase(it);
        }
      }
    }

    const auto& counters = as_object(saved_session.at("counters"),
                                     "checkpoint counters");
    require_exact_keys(counters,
        {"next_chart", "next_local", "next_scc", "next_match",
         "next_endpoint", "next_tile_plan", "next_transport_state",
         "next_line_result",
         "total_local_solves", "total_scc_column_solves",
         "total_local_matches", "total_endpoint_limits",
         "total_endpoint_exports", "total_tile_plans",
         "total_transport_arm_marches",
         "total_line_integrations", "total_line_exports",
         "total_local_run_parse_ms",
         "total_local_kernel_ms", "total_local_match_ms",
         "total_endpoint_limit_ms", "total_endpoint_export_ms",
         "total_tile_plan_ms", "total_transport_arm_ms",
         "total_line_integration_ms",
         "total_line_export_ms",
         "checkpoint_generation", "checkpoint_restore_count"},
        "checkpoint counters");
    const auto next_chart = as_u64(counters.at("next_chart"), "next chart");
    const auto next_scc = as_u64(counters.at("next_scc"), "next SCC");
    const auto next_local = as_u64(counters.at("next_local"), "next local");
    const auto next_match = as_u64(counters.at("next_match"), "next match");
    const auto next_endpoint = as_u64(
        counters.at("next_endpoint"), "next endpoint");
    const auto next_tile_plan = as_u64(
        counters.at("next_tile_plan"), "next tile plan");
    const auto next_transport_state = as_u64(
        counters.at("next_transport_state"), "next transport state");
    const auto next_line_result = as_u64(
        counters.at("next_line_result"), "next line result");
    const auto restore_count = as_u64(
        counters.at("checkpoint_restore_count"),
        "checkpoint restore count");
    if (next_chart <= largest_chart || next_scc <= largest_scc ||
        next_local <= largest_local || next_match <= largest_match ||
        next_endpoint <= largest_endpoint ||
        next_tile_plan <= largest_tile_plan ||
        next_transport_state <= largest_transport_state ||
        next_line_result <= largest_line_result)
      throw std::invalid_argument(
          "checkpoint next-handle counters do not follow retained handles");
    if (restore_count == std::numeric_limits<std::uint64_t>::max())
      throw std::overflow_error("checkpoint restore counter overflow");
    {
      std::lock_guard<std::mutex> lock(restored->mutex);
      restored->next_chart = next_chart;
      restored->next_local = next_local;
      restored->next_scc = next_scc;
      restored->next_match = next_match;
      restored->next_endpoint = next_endpoint;
      restored->next_tile_plan = next_tile_plan;
      restored->next_transport_state = next_transport_state;
      restored->next_line_result = next_line_result;
      restored->chart_capacity = configured_chart_capacity;
      restored->scc_capacity = configured_scc_capacity;
      restored->total_local_solves = as_u64(
          counters.at("total_local_solves"), "total local solves");
      restored->total_scc_column_solves = as_u64(
          counters.at("total_scc_column_solves"),
          "total SCC column solves");
      restored->total_local_matches = as_u64(
          counters.at("total_local_matches"), "total local matches");
      restored->total_endpoint_limits = as_u64(
          counters.at("total_endpoint_limits"), "total endpoint limits");
      restored->total_endpoint_exports = as_u64(
          counters.at("total_endpoint_exports"), "total endpoint exports");
      restored->total_tile_plans = as_u64(
          counters.at("total_tile_plans"), "total tile plans");
      restored->total_transport_arm_marches = as_u64(
          counters.at("total_transport_arm_marches"),
          "total transport-arm marches");
      restored->total_line_integrations = as_u64(
          counters.at("total_line_integrations"),
          "total line integrations");
      restored->total_line_exports = as_u64(
          counters.at("total_line_exports"), "total line exports");
      restored->total_local_run_parse_ms = as_double(
          counters.at("total_local_run_parse_ms"),
          "total local parse time");
      restored->total_local_kernel_ms = as_double(
          counters.at("total_local_kernel_ms"),
          "total local kernel time");
      restored->total_local_match_ms = as_double(
          counters.at("total_local_match_ms"),
          "total local match time");
      restored->total_endpoint_limit_ms = as_double(
          counters.at("total_endpoint_limit_ms"),
          "total endpoint limit time");
      restored->total_endpoint_export_ms = as_double(
          counters.at("total_endpoint_export_ms"),
          "total endpoint export time");
      restored->total_tile_plan_ms = checkpoint_nonnegative_double(
          counters.at("total_tile_plan_ms"), "total tile-plan time");
      restored->total_transport_arm_ms = checkpoint_nonnegative_double(
          counters.at("total_transport_arm_ms"),
          "total transport-arm time");
      restored->total_line_integration_ms = checkpoint_nonnegative_double(
          counters.at("total_line_integration_ms"),
          "total line-integration time");
      restored->total_line_export_ms = checkpoint_nonnegative_double(
          counters.at("total_line_export_ms"),
          "total line-export time");
      restored->checkpoint_generation = as_u64(
          counters.at("checkpoint_generation"), "checkpoint generation");
      restored->checkpoint_restore_count = restore_count + 1;
      restored->restored_from_checkpoint_identity = expected_identity;
    }
    live = false;
    return json::object{
        {"status", "ok"}, {"session", restored_handle},
        {"checkpoint_identity", expected_identity},
        {"generation", header.at("generation")},
        {"restore_count", restored->checkpoint_restore_count},
        {"configuration_identity", header.at("configuration_identity")},
        {"analytic_identity", header.at("analytic_identity")},
        {"charts", std::move(restored_charts)},
        {"sccs", std::move(restored_sccs)},
        {"locals", std::move(restored_locals)},
        {"exact_matches", std::move(restored_exact_matches)},
        {"acb_matches", std::move(restored_acb_matches)},
        {"planned_match_hops", std::move(restored_planned_matches)},
        {"endpoints", std::move(restored_endpoints)},
        {"tile_plans", std::move(restored_tile_plans)},
        {"transport_states", std::move(restored_transport_states)},
        {"line_results", std::move(restored_line_results)},
        {"deferred_handle_kinds", header.at("deferred_handle_kinds")},
        {"replayed_wolfram_preprocessing", false}};
  } catch (...) {
    if (live) {
      try {
        run_session_command(json::object{{"schema", 2},
                                         {"op", "session.close"},
                                         {"session", restored_handle}});
      } catch (...) {
      }
    }
    throw;
  }
}

template <typename Scalar>
std::shared_ptr<PreparedChartBase> parse_prepared_chart(
    const std::shared_ptr<SolverSession>& session, const json::object& root,
    const std::string& handle, const std::string& key,
    const std::string& exact_identity,
    std::optional<std::string> geometry_record,
    std::optional<std::string> principal_matrix_record,
    std::optional<std::string> native_scc_capabilities,
    SCCCertificate scc,
    std::string signature) {
  const auto started = std::chrono::steady_clock::now();
  const auto& problem = as_object(root.at("problem"), "prepared problem");
  std::unique_ptr<AcbPrecisionLease> acb_lease;
  if constexpr (std::is_same_v<Scalar, ComplexBall>) {
    acb_lease = std::make_unique<AcbPrecisionLease>(session->precision_bits);
    ComplexBall::set_precision(session->precision_bits);
  }
  std::unique_lock<std::recursive_mutex> symbolic_lock;
  if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
    symbolic_lock =
        std::unique_lock<std::recursive_mutex>(symbolic_run_mutex());
    SymbolicRational::configure(session->symbols);
  }
  auto prepared = parse_prepared_operator<Scalar>(problem);
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return make_retained_typed_shared<Scalar, PreparedChart<Scalar>>(
      handle, key, exact_identity, std::move(signature),
      std::move(geometry_record), std::move(principal_matrix_record),
      std::move(native_scc_capabilities), std::move(scc), std::move(prepared),
      session->precision_bits, session->symbols, elapsed);
}

json::object run_session_command(const json::object& root) {
  if (as_i64(root.at("schema"), "schema") != 2)
    throw std::invalid_argument("unsupported persistent solver schema");
  const auto operation = required_string(root, "op");

  if (operation == "checkpoint.restore") {
    require_exact_keys(root,
        {"schema", "op", "path", "expected_identity"},
        "checkpoint.restore request");
    return restore_checkpoint(required_string(root, "path"),
                              required_string(root, "expected_identity"));
  }

  if (operation == "session.create") {
    const auto domain = required_string(root, "domain");
    if (domain != "rational" && domain != "acb" && domain != "symbolic")
      throw std::invalid_argument("unsupported persistent scalar domain: " + domain);
    auto session = std::make_shared<SolverSession>();
    session->domain = domain;
    session->precision_bits = root.if_contains("precision_bits")
        ? static_cast<slong>(as_i64(root.at("precision_bits"), "precision bits"))
        : 256;
    if (domain == "acb" && session->precision_bits < 64)
      throw std::invalid_argument("Acb precision must be at least 64 bits");
    session->output_digits = root.if_contains("output_digits")
        ? static_cast<int>(as_i64(root.at("output_digits"), "output digits"))
        : 50;
    if (session->output_digits < 1)
      throw std::invalid_argument("output digits must be positive");
    session->symbols = parse_symbols(root);
    if (domain == "symbolic" && session->symbols.empty())
      throw std::invalid_argument(
          "symbolic persistent session requires declared regulator symbols");
    if (domain != "symbolic" && !session->symbols.empty())
      throw std::invalid_argument(
          "regulator symbols are only valid for the symbolic scalar domain");
    session->analytic_identity = root.if_contains("analytic")
        ? json::serialize(root.at("analytic")) : "null";
    if (root.if_contains("chart_capacity")) {
      const auto capacity = as_u32(root.at("chart_capacity"), "chart capacity");
      if (capacity == 0 || capacity > 4096)
        throw std::invalid_argument("chart capacity must lie in 1..4096");
      session->chart_capacity = capacity;
    }
    if (root.if_contains("local_capacity")) {
      const auto capacity = as_u32(root.at("local_capacity"), "local capacity");
      if (capacity == 0 || capacity > 16384)
        throw std::invalid_argument("local capacity must lie in 1..16384");
      session->local_capacity = capacity;
    }
    if (root.if_contains("scc_capacity")) {
      const auto capacity = as_u32(root.at("scc_capacity"), "SCC capacity");
      if (capacity == 0 || capacity > 4096)
        throw std::invalid_argument("SCC capacity must lie in 1..4096");
      session->scc_capacity = capacity;
    }
    if (root.if_contains("match_capacity")) {
      const auto capacity = as_u32(root.at("match_capacity"),
                                   "local match capacity");
      if (capacity == 0 || capacity > 16384)
        throw std::invalid_argument(
            "local match capacity must lie in 1..16384");
      session->match_capacity = capacity;
    }
    if (root.if_contains("endpoint_capacity")) {
      const auto capacity = as_u32(root.at("endpoint_capacity"),
                                   "endpoint result capacity");
      if (capacity == 0 || capacity > 16384)
        throw std::invalid_argument(
            "endpoint result capacity must lie in 1..16384");
      session->endpoint_capacity = capacity;
    }
    if (root.if_contains("transport_state_capacity")) {
      const auto capacity = as_u32(
          root.at("transport_state_capacity"),
          "transport-arm state capacity");
      if (capacity == 0 || capacity > 16384)
        throw std::invalid_argument(
            "transport-arm state capacity must lie in 1..16384");
      session->transport_state_capacity = capacity;
    }
    auto& registry = session_registry();
    {
      std::lock_guard<std::mutex> lock(registry.mutex);
      // Precision is thread-local, and every Acb solve/evaluation acquires a
      // lease and installs its session precision.  Live sessions may retain
      // different precisions; unequal active operations serialize as needed.
      session->handle = "s:" + std::to_string(registry.next_session++);
      registry.sessions.emplace(session->handle, session);
    }
    return json::object{{"status", "ok"}, {"session", session->handle},
                        {"domain", session->domain},
                        {"precision_bits", session->precision_bits},
                        {"chart_capacity", session->chart_capacity},
                        {"local_capacity", session->local_capacity},
                        {"scc_capacity", session->scc_capacity},
                        {"match_capacity", session->match_capacity},
                        {"endpoint_capacity", session->endpoint_capacity},
                        {"tile_plan_capacity", session->tile_plan_capacity},
                        {"transport_state_capacity",
                         session->transport_state_capacity},
                        {"line_result_capacity", session->line_result_capacity},
                        {"local_match_capability",
                         domain == "rational"
                             ? kExactRegularLocalMatchCapability
                             : "unsupported"},
                        {"acb_local_match_capability",
                         domain == "acb"
                             ? kRefinedAcbLocalMatchCapability
                             : "unsupported"},
                        {"endpoint_limit_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedEndpointLimitCapability},
                        {"planned_endpoint_limit_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedPlannedEndpointLimitCapability},
                        {"tile_plan_capability", kRetainedTilePlanCapability},
                        {"single_arm_tile_plan_capability",
                         kRetainedSingleArmTilePlanCapability},
                        {"planned_match_hop_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedPlannedMatchHopCapability},
                        {"planned_match_materialization_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedPlannedMatchMaterializationCapability},
                        {"rational_row_application_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedRationalRowCapability},
                        {"line_integration_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedStoredLineCapability},
                        {"parallel_arm_march_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedParallelArmCapability},
                        {"transport_arm_state_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedTransportArmStateCapability},
                        {"certified_tail_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRegularTailMajorantCapability},
                        {"certified_line_integration_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedCertifiedLineCapability}};
  }

  if (operation == "session.close") {
    const auto handle = required_string(root, "session");
    auto& registry = session_registry();
    std::shared_ptr<SolverSession> removed;
    {
      std::lock_guard<std::mutex> lock(registry.mutex);
      const auto found = registry.sessions.find(handle);
      if (found == registry.sessions.end())
        throw std::invalid_argument("unknown or already closed solver session");
      removed = std::move(found->second);
      registry.sessions.erase(found);
    }
    std::size_t charts = 0, locals = 0, matches = 0, endpoints = 0,
                tile_plans = 0, transport_states = 0,
                line_results = 0, sccs = 0;
    {
      std::lock_guard<std::mutex> lock(removed->mutex);
      removed->closed = true;
      // In-flight solve/match/transport/endpoint calls own their reservations
      // and decrement them on exactly one completion path. Do not reset
      // pending counters.
      charts = removed->charts.size();
      locals = removed->locals.size();
      matches = removed->matches.size();
      endpoints = removed->endpoints.size();
      tile_plans = removed->tile_plans.size();
      transport_states = removed->transport_states.size();
      line_results = removed->line_results.size();
      sccs = removed->sccs.size();
      removed->line_results.clear();
      removed->transport_states.clear();
      removed->tile_plans.clear();
      removed->endpoints.clear();
      removed->matches.clear();
      removed->locals.clear();
      removed->sccs.clear();
      removed->scc_handles_by_key.clear();
      removed->charts.clear();
      removed->handles_by_key.clear();
    }
    return json::object{{"status", "ok"}, {"closed", handle},
                        {"released_charts", charts},
                        {"released_locals", locals},
                        {"released_matches", matches},
                        {"released_endpoints", endpoints},
                        {"released_tile_plans", tile_plans},
                        {"released_transport_states", transport_states},
                        {"released_line_results", line_results},
                        {"released_scc_charts", sccs}};
  }

  const auto session = find_session(required_string(root, "session"));

  if (operation == "checkpoint.save") {
    require_exact_keys(root,
        {"schema", "op", "session", "path", "checkpoint_identity"},
        "checkpoint.save request");
    const auto path = required_string(root, "path");
    const auto checkpoint_identity = required_string(
        root, "checkpoint_identity");
    std::lock_guard<std::mutex> lock(session->mutex);
    auto snapshot = make_checkpoint_snapshot(*session, checkpoint_identity);
    checkpoint::write_atomic(path, json::serialize(snapshot.header),
                             json::serialize(snapshot.payload));
    session->checkpoint_generation = snapshot.generation;
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"path", path}, {"checkpoint_identity", checkpoint_identity},
        {"generation", snapshot.generation},
        {"charts", snapshot.charts}, {"sccs", snapshot.sccs},
        {"locals", snapshot.locals},
        {"exact_matches", snapshot.exact_matches},
        {"acb_matches", snapshot.acb_matches},
        {"planned_match_hops", snapshot.planned_matches},
        {"endpoints", snapshot.endpoints},
        {"tile_plans", snapshot.tile_plans},
        {"transport_states", snapshot.transport_states},
        {"line_results", snapshot.line_results},
        {"serialized_handle_kinds",
         json::array{"chart", "scc", "local", "exact-rational-match",
                     "acb-match", "planned-match-hop",
                     "materialized-local", "endpoint", "tile",
                     "transport-arm-state", "line"}},
        {"deferred_handle_kinds",
         json::array{"symbolic-local"}},
        {"atomic", true}};
  }

  if (operation == "tile.plan") {
    require_exact_keys(root,
        {"schema", "op", "session", "checkpoint_identity",
         "division_order", "lower", "upper"},
        "native tile.plan request");
    const auto& lower_request = as_object(
        root.at("lower"), "lower native tile arm");
    const auto& upper_request = as_object(
        root.at("upper"), "upper native tile arm");
    const auto lower_handles = parse_plan_chart_handles(lower_request);
    const auto upper_handles = parse_plan_chart_handles(upper_request);
    std::vector<RetainedPlanChartBinding::Owner> lower_charts;
    std::vector<RetainedPlanChartBinding::Owner> upper_charts;
    std::string plan_handle;
    {
      // Resolve every prepared-chart or composite-SCC owner and acquire strong
      // typed ownership in one admission section.  Public chart.release or
      // scc.release cannot invalidate either independently executable arm
      // after this point.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      if (session->tile_plans.size() + session->pending_tile_plans >=
          session->tile_plan_capacity)
        throw std::invalid_argument(
            "persistent native tile-plan capacity is exhausted");
      const auto resolve_owner = [&](const std::string& handle,
                                     const char* arm_name)
          -> RetainedPlanChartBinding::Owner {
        if (handle.starts_with("c:")) {
          const auto found = session->charts.find(handle);
          if (found != session->charts.end()) return found->second;
        } else if (handle.starts_with("scc:")) {
          const auto found = session->sccs.find(handle);
          if (found != session->sccs.end()) return found->second;
        }
        throw std::invalid_argument(
            std::string(
                "unknown or released prepared-chart/composite-SCC owner in ") +
            arm_name + " native tile arm: " + handle);
      };
      for (const auto& handle : lower_handles)
        lower_charts.push_back(resolve_owner(handle, "lower"));
      for (const auto& handle : upper_handles)
        upper_charts.push_back(resolve_owner(handle, "upper"));
      plan_handle = "tile:" +
          std::to_string(session->next_tile_plan++);
      ++session->pending_tile_plans;
    }
    std::shared_ptr<StoredTilePlan> plan;
    try {
      plan = build_tile_plan(plan_handle, root, lower_charts, upper_charts);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_tile_plans == 0)
        throw std::logic_error(
            "native tile-plan reservation accounting underflow");
      --session->pending_tile_plans;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_tile_plans == 0)
        throw std::logic_error(
            "native tile-plan reservation accounting underflow");
      --session->pending_tile_plans;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during native tile planning");
      session->tile_plans.emplace(plan_handle, plan);
      ++session->total_tile_plans;
      session->total_tile_plan_ms +=
          plan->summary(false).at("elapsed_ms").as_double();
    }
    auto result = plan->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "tile.plan_arm") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "checkpoint_identity",
         "division_order", "arm"},
        "native tile.plan_arm request");
    const auto& arm_request = as_object(
        root.at("arm"), "native single tile arm");
    const auto handles = parse_plan_chart_handles(arm_request);
    std::vector<RetainedPlanChartBinding::Owner> charts;
    std::string plan_handle;
    {
      // Acquire every chart/SCC owner and the publication reservation in one
      // admission section. The plan keeps these owners alive even if their
      // public handles are released while exact planning runs.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      if (session->tile_plans.size() + session->pending_tile_plans >=
          session->tile_plan_capacity)
        throw std::invalid_argument(
            "persistent native tile-plan capacity is exhausted");
      const auto resolve_owner = [&](const std::string& handle)
          -> RetainedPlanChartBinding::Owner {
        if (handle.starts_with("c:")) {
          const auto found = session->charts.find(handle);
          if (found != session->charts.end()) return found->second;
        } else if (handle.starts_with("scc:")) {
          const auto found = session->sccs.find(handle);
          if (found != session->sccs.end()) return found->second;
        }
        throw std::invalid_argument(
            "unknown or released prepared-chart/composite-SCC owner in native single tile arm: " +
            handle);
      };
      charts.reserve(handles.size());
      for (const auto& handle : handles)
        charts.push_back(resolve_owner(handle));
      plan_handle = "tile:" +
          std::to_string(session->next_tile_plan++);
      ++session->pending_tile_plans;
    }
    std::shared_ptr<StoredTilePlan> plan;
    try {
      plan = build_single_arm_tile_plan(
          plan_handle, root, charts);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_tile_plans == 0)
        throw std::logic_error(
            "native single-arm tile-plan reservation accounting underflow");
      --session->pending_tile_plans;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_tile_plans == 0)
        throw std::logic_error(
            "native single-arm tile-plan reservation accounting underflow");
      --session->pending_tile_plans;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during native single-arm tile planning");
      if (!session->tile_plans.emplace(plan_handle, plan).second)
        throw std::logic_error(
            "native single-arm tile-plan handle collided during publication");
      ++session->total_tile_plans;
      session->total_tile_plan_ms +=
          plan->summary(false).at("elapsed_ms").as_double();
    }
    auto result = plan->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "tile.stats" || operation == "tile.match_interval" ||
      operation == "tile.integration_interval") {
    const auto plan_handle = required_string(root, "tile_plan");
    std::shared_ptr<StoredTilePlan> plan;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->tile_plans.find(plan_handle);
      if (found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or released native tile plan");
      plan = found->second;
    }
    json::object result;
    if (operation == "tile.stats") {
      result = plan->summary();
    } else {
      const auto arm = required_string(root, "arm");
      const auto index = static_cast<std::size_t>(as_u64(
          root.at(operation == "tile.match_interval" ? "match" : "tile"),
          operation == "tile.match_interval" ? "tile-plan match index"
                                               : "tile-plan tile index"));
      result = operation == "tile.match_interval"
          ? plan->match_interval(arm, index)
          : plan->tile_interval(arm, index);
    }
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "tile.release") {
    const auto plan_handle = required_string(root, "tile_plan");
    std::shared_ptr<StoredTilePlan> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->tile_plans.find(plan_handle);
      if (found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or already released native tile plan");
      removed = std::move(found->second);
      session->tile_plans.erase(found);
    }
    return json::object{
        {"status", "ok"}, {"released", plan_handle},
        {"checkpoint_identity", removed->checkpoint_identity()}};
  }

  if (operation == "tile.match_advance") {
    if (session->domain == "rational")
      require_exact_keys(
          root,
          {"schema", "op", "session", "tile_plan", "arm", "match",
           "basis", "incoming", "epsilon", "checkpoint_identity"},
          "native tile.match_advance request");
    else if (session->domain == "acb")
      require_exact_keys(
          root,
          {"schema", "op", "session", "tile_plan", "arm", "match",
           "basis", "incoming", "epsilon", "refinement",
           "exact_lattice", "checkpoint_identity"},
          "native Acb tile.match_advance request");
    else
      throw std::invalid_argument(
          "native tile.match_advance requires rational or Acb coefficients");

    const auto plan_handle = required_string(root, "tile_plan");
    const auto incoming_handle = required_string(root, "incoming");
    const auto& raw_basis = as_array(
        root.at("basis"), "planned local match basis");
    if (raw_basis.empty())
      throw std::invalid_argument(
          "planned local match basis cannot be empty");
    std::vector<std::string> basis_handles;
    basis_handles.reserve(raw_basis.size());
    std::set<std::string> unique_handles;
    for (const auto& raw_handle : raw_basis) {
      if (!raw_handle.is_string() || raw_handle.as_string().empty())
        throw std::invalid_argument(
            "planned local match basis handles must be nonempty strings");
      std::string handle(raw_handle.as_string());
      if (!unique_handles.insert(handle).second)
        throw std::invalid_argument(
            "planned local match basis handles must be pairwise distinct");
      basis_handles.push_back(std::move(handle));
    }
    if (unique_handles.contains(incoming_handle))
      throw std::invalid_argument(
          "planned incoming local must be distinct from its basis");

    std::shared_ptr<StoredTilePlan> plan;
    std::vector<std::shared_ptr<StoredLocalBase>> basis;
    std::shared_ptr<StoredLocalBase> incoming;
    std::string match_handle;
    {
      // Admission is the only serialized section.  Each lower/upper hop owns
      // an immutable plan snapshot and strong local references while its
      // matching arithmetic runs independently outside the session lock.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto plan_found = session->tile_plans.find(plan_handle);
      if (plan_found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or released tile plan for planned local matching");
      if (session->matches.size() + session->pending_matches >=
          session->match_capacity)
        throw std::invalid_argument(
            "persistent local match capacity is exhausted");
      plan = plan_found->second;
      for (const auto& handle : basis_handles) {
        const auto found = session->locals.find(handle);
        if (found == session->locals.end())
          throw std::invalid_argument(
              "unknown or released native local in planned match basis: " +
              handle);
        basis.push_back(found->second);
      }
      const auto incoming_found = session->locals.find(incoming_handle);
      if (incoming_found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released incoming native local for planned match: " +
            incoming_handle);
      incoming = incoming_found->second;
      match_handle = "m:" + std::to_string(session->next_match++);
      ++session->pending_matches;
    }

    std::shared_ptr<StoredPlannedMatchHop> match;
    try {
      match = build_planned_match_hop(
          match_handle, root, session->domain, session->precision_bits,
          checkpoint_configuration_identity(*session), plan, basis_handles,
          basis, incoming_handle, incoming);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches == 0)
        throw std::logic_error(
            "native planned match reservation accounting underflow");
      --session->pending_matches;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches == 0)
        throw std::logic_error(
            "native planned match reservation accounting underflow");
      --session->pending_matches;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during plan-driven local matching");
      session->matches.emplace(match_handle, match);
      ++session->total_local_matches;
      session->total_local_match_ms += match->elapsed_ms();
      plan->note_match_advance(required_string(root, "arm"));
    }
    auto response = match->summary();
    response["status"] = "ok";
    response["session"] = session->handle;
    return response;
  }

  if (operation == "transport.run_arm") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "tile_plan", "anchor",
         "tile_plan_checkpoint_identity", "anchor_checkpoint_identity",
         "arm", "receiving_basis", "epsilon", "refinement",
         "checkpoint_policy"},
        "native transport.run_arm request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native transport-arm marching requires rational or Acb coefficients");
    const auto arm_name = required_string(root, "arm");
    if (arm_name != "lower" && arm_name != "upper")
      throw std::invalid_argument(
          "native transport-arm name must be lower or upper");
    const auto epsilon_contract = parse_whole_arm_epsilon_contract(
        root.at("epsilon"), "native transport-arm epsilon contract");
    const auto& refinement = as_object(
        root.at("refinement"), "native transport-arm refinement policy");
    require_exact_keys(refinement, {"relative_tolerance", "max_steps"},
                       "native transport-arm refinement policy");
    (void)required_string(refinement, "relative_tolerance");
    if (as_u32(refinement.at("max_steps"),
               "native transport-arm refinement steps") > 32)
      throw std::invalid_argument(
          "native transport-arm refinement steps must lie in 0..32");
    const auto& checkpoint_policy = as_object(
        root.at("checkpoint_policy"),
        "native transport-arm checkpoint policy");
    require_exact_keys(checkpoint_policy, {"schema", "root"},
                       "native transport-arm checkpoint policy");
    if (required_string(checkpoint_policy, "schema") !=
        "diffexp2-deterministic-arm-checkpoints-v1")
      throw std::invalid_argument(
          "unsupported native transport-arm checkpoint policy schema");
    const auto checkpoint_root = required_string(checkpoint_policy, "root");
    if (checkpoint_root.empty())
      throw std::invalid_argument(
          "native transport-arm checkpoint root cannot be empty");

    RetainedArmMarchInput input;
    input.name = arm_name;
    const auto& raw_basis_sets = as_array(
        root.at("receiving_basis"),
        "native transport-arm receiving basis sets");
    input.basis_handles.reserve(raw_basis_sets.size());
    for (const auto& raw_set : raw_basis_sets) {
      const auto& values = as_array(
          raw_set, "native transport-arm receiving basis set");
      if (values.empty())
        throw std::invalid_argument(
            "native transport-arm receiving basis sets cannot be empty");
      std::set<std::string> unique;
      std::vector<std::string> handles;
      handles.reserve(values.size());
      for (const auto& raw_handle : values) {
        if (!raw_handle.is_string() || raw_handle.as_string().empty())
          throw std::invalid_argument(
              "native transport-arm basis handles must be nonempty strings");
        std::string handle(raw_handle.as_string());
        if (!unique.insert(handle).second)
          throw std::invalid_argument(
              "native transport-arm basis handles must be pairwise distinct");
        handles.push_back(std::move(handle));
      }
      input.basis_handles.push_back(std::move(handles));
    }

    const auto plan_handle = required_string(root, "tile_plan");
    const auto anchor_handle = required_string(root, "anchor");
    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> anchor;
    std::string state_handle;
    const auto match_count = input.basis_handles.size();
    bool reservation_live = false;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto plan_found = session->tile_plans.find(plan_handle);
      if (plan_found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or released tile plan for transport-arm marching");
      plan = plan_found->second;
      if (required_string(root, "tile_plan_checkpoint_identity") !=
          plan->checkpoint_identity())
        throw std::invalid_argument(
            "transport-arm tile-plan checkpoint token is stale");
      const auto& retained = plan->arm(arm_name);
      if (retained.exact.matches.size() != match_count ||
          retained.exact.tiles.size() != match_count + 1)
        throw std::invalid_argument(
            "transport-arm basis count does not reproduce the retained plan topology");
      const auto anchor_found = session->locals.find(anchor_handle);
      if (anchor_found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released anchor local for transport-arm marching");
      anchor = anchor_found->second;
      if (required_string(root, "anchor_checkpoint_identity") !=
          anchor->checkpoint_identity())
        throw std::invalid_argument(
            "transport-arm anchor checkpoint token is stale");
      const auto& anchor_binding = retained.charts.front();
      if (anchor->source_chart() != anchor_binding.handle)
        throw std::invalid_argument(
            "transport-arm anchor belongs to a different retained chart");
      anchor->require_exact_plan_binding(
          anchor_binding.geometry, anchor_binding.prescriptions,
          "transport-arm anchor");

      input.basis.reserve(match_count);
      for (const auto& handles : input.basis_handles) {
        std::vector<std::shared_ptr<StoredLocalBase>> resolved;
        resolved.reserve(handles.size());
        for (const auto& handle : handles) {
          const auto found = session->locals.find(handle);
          if (found == session->locals.end())
            throw std::invalid_argument(
                "unknown or released native local in transport-arm receiving basis: " +
                handle);
          resolved.push_back(found->second);
        }
        input.basis.push_back(std::move(resolved));
      }

      if (match_count > session->match_capacity -
                            std::min(session->match_capacity,
                                     session->matches.size() +
                                         session->pending_matches))
        throw std::invalid_argument(
            "persistent local match capacity is exhausted by transport-arm marching");
      if (match_count > session->local_capacity -
                            std::min(session->local_capacity,
                                     session->locals.size() +
                                         session->pending_local_solves))
        throw std::invalid_argument(
            "persistent local capacity is exhausted by transport-arm marching");
      if (session->transport_states.size() +
              session->pending_transport_states >=
          session->transport_state_capacity)
        throw std::invalid_argument(
            "persistent transport-state capacity is exhausted");

      input.match_handles.reserve(match_count);
      input.local_handles.reserve(match_count);
      for (std::size_t index = 0; index < match_count; ++index) {
        input.match_handles.push_back(
            "m:" + std::to_string(session->next_match++));
        input.local_handles.push_back(
            "l:" + std::to_string(session->next_local++));
      }
      state_handle = "transport:" +
          std::to_string(session->next_transport_state++);
      session->pending_matches += match_count;
      session->pending_local_solves += match_count;
      ++session->pending_transport_states;
      reservation_live = true;
    }
    const auto release_reservation = [&]() {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (!reservation_live) return;
      if (session->pending_matches < match_count ||
          session->pending_local_solves < match_count ||
          session->pending_transport_states == 0)
        throw std::logic_error(
            "native transport-arm reservation accounting underflow");
      session->pending_matches -= match_count;
      session->pending_local_solves -= match_count;
      --session->pending_transport_states;
      reservation_live = false;
    };
    struct TransportArmReservationGuard final {
      decltype(release_reservation)& release;
      bool& live;
      ~TransportArmReservationGuard() noexcept {
        if (!live) return;
        try {
          release();
        } catch (...) {
          std::terminate();
        }
      }
    } reservation_guard{release_reservation, reservation_live};

    const auto operation_started = std::chrono::steady_clock::now();
    std::unique_ptr<AcbPrecisionLease> acb_lease;
    if (session->domain == "acb") {
      acb_lease = std::make_unique<AcbPrecisionLease>(
          session->precision_bits);
      ComplexBall::set_precision(session->precision_bits);
    }
    const auto session_configuration =
        checkpoint_configuration_identity(*session);
    auto marched = march_retained_arm(
        session->domain, session->precision_bits, session_configuration,
        plan, anchor, input, epsilon_contract.work,
        epsilon_contract.match_required_complete_max, refinement,
        checkpoint_root);
    auto state = std::make_shared<StoredTransportArmState>(
        state_handle, checkpoint_root + ":" + arm_name + ":state",
        arm_name, plan, anchor, input.basis, marched.matches,
        marched.tile_sources, epsilon_contract.work,
        epsilon_contract.public_required_complete_max,
        epsilon_contract.match_required_complete_max, refinement,
        marched.elapsed_ms);
    const auto final_local = state->final_local();

    try {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches < match_count ||
          session->pending_local_solves < match_count ||
          session->pending_transport_states == 0)
        throw std::logic_error(
            "native transport-arm reservation accounting underflow");
      session->pending_matches -= match_count;
      session->pending_local_solves -= match_count;
      --session->pending_transport_states;
      reservation_live = false;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during transport-arm marching");

      session->locals.reserve(session->locals.size() +
                              (match_count == 0 ? 0 : 1));
      session->transport_states.reserve(
          session->transport_states.size() + 1);
      bool inserted_local = false;
      try {
        const auto existing = session->locals.find(final_local->handle());
        if (existing == session->locals.end()) {
          if (!session->locals.emplace(final_local->handle(), final_local)
                   .second)
            throw std::logic_error(
                "transport-arm final-local handle collision at publication");
          inserted_local = true;
        } else if (existing->second.get() != final_local.get()) {
          throw std::logic_error(
              "transport-arm final-local handle names a different retained object");
        }
        if (!session->transport_states.emplace(state_handle, state).second)
          throw std::logic_error(
              "transport-arm state handle collision at publication");
      } catch (...) {
        session->transport_states.erase(state_handle);
        if (inserted_local) session->locals.erase(final_local->handle());
        throw;
      }
      for (const auto& match : state->matches()) {
        ++session->total_local_matches;
        session->total_local_match_ms += match->elapsed_ms();
        plan->note_match_advance(arm_name);
      }
      ++session->total_transport_arm_marches;
      session->total_transport_arm_ms += state->elapsed_ms();
    } catch (...) {
      if (reservation_live) release_reservation();
      throw;
    }

    auto response = state->summary();
    auto final_summary = final_local->summary();
    final_summary["session"] = session->handle;
    response["final_local"] = std::move(final_summary);
    response["status"] = "ok";
    response["session"] = session->handle;
    response["atomic_publication"] = true;
    response["checkpoint_policy"] = checkpoint_policy;
    response["operation_elapsed_ms"] = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - operation_started).count();
    return response;
  }

  if (operation == "transport.stats" ||
      operation == "transport.release") {
    require_exact_keys(root,
        {"schema", "op", "session", "transport_state"},
        operation == "transport.stats"
            ? "native transport.stats request"
            : "native transport.release request");
    const auto state_handle = required_string(root, "transport_state");
    if (operation == "transport.stats") {
      std::shared_ptr<StoredTransportArmState> state;
      {
        std::lock_guard<std::mutex> lock(session->mutex);
        const auto found = session->transport_states.find(state_handle);
        if (found == session->transport_states.end())
          throw std::invalid_argument(
              "unknown or released native transport-arm state");
        state = found->second;
      }
      auto result = state->stats_json();
      result["status"] = "ok";
      result["session"] = session->handle;
      return result;
    }
    std::shared_ptr<StoredTransportArmState> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->transport_states.find(state_handle);
      if (found == session->transport_states.end())
        throw std::invalid_argument(
            "unknown or already released native transport-arm state");
      removed = std::move(found->second);
      session->transport_states.erase(found);
    }
    return json::object{
        {"status", "ok"}, {"released", state_handle},
        {"checkpoint_identity", removed->checkpoint_identity()}};
  }

  if (operation == "integration.run_arms") {
    if (root.if_contains("certify_tail") != nullptr)
      require_exact_keys(
          root,
          {"schema", "op", "session", "tile_plan", "anchor",
           "tile_plan_checkpoint_identity", "anchor_checkpoint_identity",
           "epsilon", "refinement", "checkpoint_policy", "lower", "upper",
           "certify_tail"},
          "native integration.run_arms request");
    else
      require_exact_keys(
          root,
          {"schema", "op", "session", "tile_plan", "anchor",
           "tile_plan_checkpoint_identity", "anchor_checkpoint_identity",
           "epsilon", "refinement", "checkpoint_policy", "lower", "upper"},
          "native integration.run_arms request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native whole-arm marching requires rational or Acb coefficients");
    const bool certify_tail = root.if_contains("certify_tail") != nullptr
        ? root.at("certify_tail").as_bool()
        : false;

    const auto& raw_epsilon = as_object(
        root.at("epsilon"), "native whole-arm epsilon contract");
    const auto epsilon_contract = parse_whole_arm_epsilon_contract(
        raw_epsilon, "native whole-arm epsilon contract");
    const auto work_epsilon = epsilon_contract.work;
    const auto required_complete_max =
        epsilon_contract.public_required_complete_max;
    const auto match_required_complete_max =
        epsilon_contract.match_required_complete_max;

    const auto& refinement = as_object(
        root.at("refinement"), "native whole-arm refinement policy");
    require_exact_keys(refinement, {"relative_tolerance", "max_steps"},
                       "native whole-arm refinement policy");
    (void)required_string(refinement, "relative_tolerance");
    if (as_u32(refinement.at("max_steps"),
               "native whole-arm refinement steps") > 32)
      throw std::invalid_argument(
          "native whole-arm refinement steps must lie in 0..32");

    const auto& checkpoint_policy = as_object(
        root.at("checkpoint_policy"),
        "native whole-arm checkpoint policy");
    require_exact_keys(checkpoint_policy, {"schema", "root"},
                       "native whole-arm checkpoint policy");
    if (required_string(checkpoint_policy, "schema") !=
        "diffexp2-deterministic-arm-checkpoints-v1")
      throw std::invalid_argument(
          "unsupported native whole-arm checkpoint policy schema");
    const auto checkpoint_root = required_string(checkpoint_policy, "root");
    if (checkpoint_root.empty())
      throw std::invalid_argument(
          "native whole-arm checkpoint root cannot be empty");

    struct PendingArmMarch {
      std::string name;
      std::vector<std::vector<std::string>> basis_handles;
      std::vector<json::object> integrand_rows;
      std::vector<std::vector<std::shared_ptr<StoredLocalBase>>> basis;
      std::vector<std::string> match_handles;
      std::vector<std::string> local_handles;
      std::vector<std::string> row_local_handles;
      std::string aggregate_handle;
    };
    const auto parse_pending_arm = [&](const char* name) {
      const auto& raw_arm = as_object(root.at(name),
                                      "native whole-arm arm request");
      require_exact_keys(raw_arm, {"receiving_basis", "integrand_rows"},
                         "native whole-arm arm request");
      PendingArmMarch arm;
      arm.name = name;
      const auto& raw_basis_sets = as_array(
          raw_arm.at("receiving_basis"),
          "native whole-arm receiving basis sets");
      arm.basis_handles.reserve(raw_basis_sets.size());
      for (const auto& raw_set : raw_basis_sets) {
        const auto& values = as_array(
            raw_set, "native whole-arm receiving basis set");
        if (values.empty())
          throw std::invalid_argument(
              "native whole-arm receiving basis sets cannot be empty");
        std::set<std::string> unique;
        std::vector<std::string> handles;
        handles.reserve(values.size());
        for (const auto& raw_handle : values) {
          if (!raw_handle.is_string() || raw_handle.as_string().empty())
            throw std::invalid_argument(
                "native whole-arm basis handles must be nonempty strings");
          std::string handle(raw_handle.as_string());
          if (!unique.insert(handle).second)
            throw std::invalid_argument(
                "native whole-arm basis handles must be pairwise distinct");
          handles.push_back(std::move(handle));
        }
        arm.basis_handles.push_back(std::move(handles));
      }
      const auto& raw_rows = as_array(
          raw_arm.at("integrand_rows"),
          "native whole-arm integrand rows");
      arm.integrand_rows.reserve(raw_rows.size());
      for (const auto& raw_row : raw_rows)
        arm.integrand_rows.push_back(
            as_object(raw_row, "native whole-arm integrand row"));
      return arm;
    };
    std::array<PendingArmMarch, 2> arms{
        parse_pending_arm("lower"), parse_pending_arm("upper")};

    const auto plan_handle = required_string(root, "tile_plan");
    const auto anchor_handle = required_string(root, "anchor");
    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> anchor;
    std::string combined_handle;
    const std::size_t total_matches =
        arms[0].basis_handles.size() + arms[1].basis_handles.size();
    std::size_t total_tiles = 0;
    std::size_t retained_local_reservation = 0;
    constexpr std::size_t published_line_results = 3;
    bool reservation_live = false;
    {
      // Resolve every source token and reserve the complete publication set
      // before either worker exists.  Subsequent public releases cannot
      // invalidate the acquired strong owners.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto plan_found = session->tile_plans.find(plan_handle);
      if (plan_found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or released tile plan for whole-arm marching");
      plan = plan_found->second;
      if (required_string(root, "tile_plan_checkpoint_identity") !=
          plan->checkpoint_identity())
        throw std::invalid_argument(
            "whole-arm tile-plan checkpoint token is stale");
      const auto anchor_found = session->locals.find(anchor_handle);
      if (anchor_found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released anchor local for whole-arm marching");
      anchor = anchor_found->second;
      if (required_string(root, "anchor_checkpoint_identity") !=
          anchor->checkpoint_identity())
        throw std::invalid_argument(
            "whole-arm anchor checkpoint token is stale");

      for (auto& arm : arms) {
        const auto& retained = plan->arm(arm.name);
        if (retained.exact.matches.size() != arm.basis_handles.size() ||
            arm.integrand_rows.size() != retained.exact.tiles.size() ||
            retained.exact.tiles.size() != arm.basis_handles.size() + 1)
          throw std::invalid_argument(
              "whole-arm basis/row counts do not reproduce the retained plan topology for " +
              arm.name);
        total_tiles = checked_diagnostic_sum(
            total_tiles, retained.exact.tiles.size(),
            "whole-arm tile count");
        arm.basis.reserve(arm.basis_handles.size());
        for (const auto& handles : arm.basis_handles) {
          std::vector<std::shared_ptr<StoredLocalBase>> resolved;
          resolved.reserve(handles.size());
          for (const auto& handle : handles) {
            const auto found = session->locals.find(handle);
            if (found == session->locals.end())
              throw std::invalid_argument(
                  "unknown or released native local in whole-arm receiving basis: " +
                  handle);
            resolved.push_back(found->second);
          }
          arm.basis.push_back(std::move(resolved));
        }
      }

      if (total_matches > session->match_capacity -
                              std::min(session->match_capacity,
                                       session->matches.size() +
                                           session->pending_matches))
        throw std::invalid_argument(
            "persistent local match capacity is exhausted by whole-arm marching");
      retained_local_reservation = checked_diagnostic_sum(
          total_matches, total_tiles,
          "whole-arm retained local reservation");
      if (retained_local_reservation > session->local_capacity -
                              std::min(session->local_capacity,
                                       session->locals.size() +
                                           session->pending_local_solves))
        throw std::invalid_argument(
            "persistent local capacity is exhausted by whole-arm marching");
      if (published_line_results > session->line_result_capacity -
              std::min(session->line_result_capacity,
                       session->line_results.size() +
                           session->pending_line_integrations))
        throw std::invalid_argument(
            "persistent line-result capacity is exhausted by whole-arm marching");

      for (auto& arm : arms) {
        arm.match_handles.reserve(arm.basis_handles.size());
        arm.local_handles.reserve(arm.basis_handles.size());
        arm.row_local_handles.reserve(arm.integrand_rows.size());
        for (std::size_t index = 0; index < arm.basis_handles.size(); ++index) {
          arm.match_handles.push_back(
              "m:" + std::to_string(session->next_match++));
          arm.local_handles.push_back(
              "l:" + std::to_string(session->next_local++));
        }
        for (std::size_t index = 0; index < arm.integrand_rows.size(); ++index)
          arm.row_local_handles.push_back(
              "l:" + std::to_string(session->next_local++));
        arm.aggregate_handle =
            "line:" + std::to_string(session->next_line_result++);
      }
      combined_handle =
          "line:" + std::to_string(session->next_line_result++);
      session->pending_matches += total_matches;
      session->pending_local_solves += retained_local_reservation;
      session->pending_line_integrations += published_line_results;
      reservation_live = true;
    }
    const auto release_reservation = [&]() {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (!reservation_live) return;
      if (session->pending_matches < total_matches ||
          session->pending_local_solves < retained_local_reservation ||
          session->pending_line_integrations < published_line_results)
        throw std::logic_error(
            "native whole-arm reservation accounting underflow");
      session->pending_matches -= total_matches;
      session->pending_local_solves -= retained_local_reservation;
      session->pending_line_integrations -= published_line_results;
      reservation_live = false;
    };
    struct WholeArmReservationGuard final {
      decltype(release_reservation)& release;
      bool& live;

      ~WholeArmReservationGuard() noexcept {
        if (!live) return;
        try {
          release();
        } catch (...) {
          // Reservation underflow is an internal invariant failure.  It is
          // unsafe to continue unwinding with counters that may still admit
          // work beyond the configured capacities.
          std::terminate();
        }
      }
    } reservation_guard{release_reservation, reservation_live};

    struct CompletedArmMarch {
      RetainedArmMarchResult march;
      std::vector<std::shared_ptr<StoredLocalBase>> projected;
      std::vector<std::shared_ptr<StoredLineResult>> tile_lines;
      std::vector<std::shared_ptr<StoredLocalBase>> tile_sources;
      std::shared_ptr<StoredLineResult> aggregate;
      double elapsed_ms = 0.0;
    };
    std::array<CompletedArmMarch, 2> completed;
    std::array<std::exception_ptr, 2> failures;
    std::atomic<std::size_t> active_workers{0};
    std::atomic<std::size_t> max_active_workers{0};
    std::mutex start_mutex;
    std::condition_variable start_changed;
    std::size_t workers_ready = 0;
    bool workers_start = false;
    bool workers_cancel = false;
    const auto operation_started = std::chrono::steady_clock::now();
    std::unique_ptr<AcbPrecisionLease> whole_arm_acb_lease;
    if (session->domain == "acb") {
      whole_arm_acb_lease =
          std::make_unique<AcbPrecisionLease>(session->precision_bits);
      ComplexBall::set_precision(session->precision_bits);
    }
    const auto active_session_configuration_identity =
        checkpoint_configuration_identity(*session);

    const auto update_max_active = [&](std::size_t candidate) {
      auto observed = max_active_workers.load();
      while (observed < candidate &&
             !max_active_workers.compare_exchange_weak(observed, candidate)) {
      }
    };
    const auto run_arm = [&](std::size_t arm_index) {
      auto active = active_workers.fetch_add(1) + 1;
      update_max_active(active);
      {
        std::unique_lock<std::mutex> lock(start_mutex);
        ++workers_ready;
        start_changed.notify_all();
        start_changed.wait(lock, [&] { return workers_start; });
        if (workers_cancel) {
          active_workers.fetch_sub(1);
          return;
        }
      }
      try {
        if (session->domain == "acb")
          ComplexBall::set_precision(session->precision_bits);
        const auto started = std::chrono::steady_clock::now();
        auto& input = arms[arm_index];
        auto& output = completed[arm_index];
        const auto& retained = plan->arm(input.name);
        RetainedArmMarchInput march_input{
            input.name, input.basis_handles, input.basis,
            input.match_handles, input.local_handles};
        output.march = march_retained_arm(
            session->domain, session->precision_bits,
            active_session_configuration_identity, plan, anchor,
            march_input, work_epsilon, match_required_complete_max,
            refinement, checkpoint_root);
        output.projected.reserve(retained.exact.tiles.size());
        output.tile_lines.reserve(retained.exact.tiles.size());
        output.tile_sources.reserve(retained.exact.tiles.size());
        for (std::size_t tile = 0; tile < retained.exact.tiles.size(); ++tile) {
          const auto& current = output.march.tile_sources[tile];
          const auto row_identity = required_string(
              input.integrand_rows[tile], "exact_identity");
          json::object row_request{
              {"row", input.integrand_rows[tile]},
              {"source_checkpoint_identity", current->checkpoint_identity()},
              {"checkpoint_identity",
               arm_checkpoint_identity(checkpoint_root, input.name,
                                       "integrand", tile + 1) + ":" +
                   row_identity}};
          std::shared_ptr<StoredLocalBase> projected;
          if (session->domain == "rational") {
            const auto typed =
                std::dynamic_pointer_cast<StoredLocal<Rational>>(current);
            if (!typed)
              throw std::logic_error(
                  "whole-arm Rational integrand source changed coefficient domain");
            projected = build_rational_row_local<Rational>(
                input.row_local_handles[tile], row_request,
                session->precision_bits, typed, current);
          } else {
            const auto typed =
                std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(current);
            if (!typed)
              throw std::logic_error(
                  "whole-arm Acb integrand source changed coefficient domain");
            projected = build_rational_row_local<ComplexBall>(
                input.row_local_handles[tile], row_request,
                session->precision_bits, typed, current);
          }
          const auto line_epsilon = live_line_epsilon_intersection(
              work_epsilon, required_complete_max, projected);
          output.projected.push_back(projected);
          output.tile_sources.push_back(projected);
          json::object line_request{
              {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
              {"source_checkpoint_identity", projected->checkpoint_identity()},
              {"checkpoint_identity",
               arm_checkpoint_identity(checkpoint_root, input.name,
                                       "tile", tile + 1)},
              {"arm", input.name}, {"tile", tile},
              {"epsilon", json::object{{"min", line_epsilon.min_power},
                                        {"max", line_epsilon.complete_max}}}};
          if (certify_tail) line_request["certify_tail"] = true;
          output.tile_lines.push_back(build_planned_line_result(
              "private:" + arm_checkpoint_identity(
                  checkpoint_root, input.name, "tile", tile + 1),
              line_request, plan, projected));
        }
        const auto aggregate_started = std::chrono::steady_clock::now();
        std::vector<std::int32_t> signs(output.tile_lines.size(), 1);
        output.aggregate = build_retained_line_aggregate(
            input.aggregate_handle,
            checkpoint_root + ":" + input.name + ":aggregate",
            input.name,
            json::object{{"from_exact", retained.exact.from.str()},
                         {"to_exact", retained.exact.to.str()}},
            json::object{{"kind", "complete-retained-arm"},
                         {"combination", "sum-physical-tiles"},
                         {"epsilon_contract", raw_epsilon},
                         {"certify_tail_requested", certify_tail}},
            plan, output.tile_sources, output.tile_lines, signs,
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - aggregate_started).count());
        output.elapsed_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - started).count();
      } catch (...) {
        failures[arm_index] = std::current_exception();
      }
      active_workers.fetch_sub(1);
    };

    std::vector<std::jthread> workers;
    workers.reserve(2);
    try {
      workers.emplace_back([&] { run_arm(0); });
      workers.emplace_back([&] { run_arm(1); });
    } catch (...) {
      {
        std::lock_guard<std::mutex> lock(start_mutex);
        workers_cancel = true;
        workers_start = true;
      }
      start_changed.notify_all();
      for (auto& worker : workers)
        if (worker.joinable()) worker.join();
      release_reservation();
      throw;
    }
    {
      std::unique_lock<std::mutex> lock(start_mutex);
      start_changed.wait(lock, [&] { return workers_ready == 2; });
      workers_start = true;
    }
    start_changed.notify_all();
    for (auto& worker : workers) worker.join();

    if (failures[0] || failures[1]) {
      release_reservation();
      // No worker result has entered a session registry.  Prefer the lower
      // failure only for deterministic reporting when both arms fail.
      std::rethrow_exception(failures[0] ? failures[0] : failures[1]);
    }

    std::vector<std::shared_ptr<StoredLineResult>> arm_lines{
        completed[0].aggregate, completed[1].aggregate};
    std::vector<std::shared_ptr<StoredLocalBase>> combined_owners =
        completed[0].tile_sources;
    combined_owners.insert(combined_owners.end(),
                           completed[1].tile_sources.begin(),
                           completed[1].tile_sources.end());
    const auto combined_started = std::chrono::steady_clock::now();
    auto combined = build_retained_line_aggregate(
        combined_handle, checkpoint_root + ":combined", "combined",
        json::object{
            {"from_exact", plan->arm("lower").exact.to.str()},
            {"to_exact", plan->arm("upper").exact.to.str()}},
        json::object{
            {"kind", "complete-lower-to-upper-line"},
            {"combination", "negative-lower-anchor-arm-plus-upper-anchor-arm"},
            {"epsilon_contract", raw_epsilon},
            {"certify_tail_requested", certify_tail}},
        plan, std::move(combined_owners), arm_lines,
        std::vector<std::int32_t>{-1, 1},
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - combined_started).count());

    try {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches < total_matches ||
          session->pending_local_solves < retained_local_reservation ||
          session->pending_line_integrations < published_line_results)
        throw std::logic_error(
            "native whole-arm reservation accounting underflow");
      session->pending_matches -= total_matches;
      session->pending_local_solves -= retained_local_reservation;
      session->pending_line_integrations -= published_line_results;
      reservation_live = false;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during whole-arm marching");

      session->locals.reserve(session->locals.size() + 2);
      session->line_results.reserve(
          session->line_results.size() + published_line_results);
      std::vector<std::string> inserted_locals;
      std::vector<std::string> inserted_lines;
      try {
        for (std::size_t arm_index = 0; arm_index < 2; ++arm_index) {
          const auto& final_local =
              completed[arm_index].march.final_local();
          const auto existing = session->locals.find(final_local->handle());
          if (existing == session->locals.end()) {
            if (!session->locals.emplace(
                    final_local->handle(), final_local).second)
              throw std::logic_error(
                  "whole-arm final-local handle collision at publication");
            inserted_locals.push_back(final_local->handle());
          } else if (existing->second.get() != final_local.get()) {
            throw std::logic_error(
                "whole-arm final-local handle names a different retained object");
          }
          const auto& line = completed[arm_index].aggregate;
          if (!session->line_results.emplace(line->handle(), line).second)
            throw std::logic_error(
                "whole-arm aggregate handle collision at publication");
          inserted_lines.push_back(line->handle());
        }
        if (!session->line_results.emplace(combined->handle(), combined).second)
          throw std::logic_error(
              "whole-arm combined handle collision at publication");
        inserted_lines.push_back(combined->handle());
      } catch (...) {
        for (const auto& handle : inserted_lines)
          session->line_results.erase(handle);
        for (const auto& handle : inserted_locals)
          session->locals.erase(handle);
        throw;
      }

      for (std::size_t arm_index = 0; arm_index < 2; ++arm_index) {
        for (const auto& match : completed[arm_index].march.matches) {
          ++session->total_local_matches;
          session->total_local_match_ms += match->elapsed_ms();
          plan->note_match_advance(arms[arm_index].name);
        }
        for (const auto& line : completed[arm_index].tile_lines) {
          ++session->total_line_integrations;
          session->total_line_integration_ms += line->elapsed_ms();
          plan->note_integration();
        }
      }
    } catch (...) {
      if (reservation_live) release_reservation();
      throw;
    }

    json::object arm_response;
    for (std::size_t index = 0; index < 2; ++index) {
      auto final_local = completed[index].march.final_local()->summary();
      final_local["session"] = session->handle;
      auto line = completed[index].aggregate->summary();
      line["session"] = session->handle;
      arm_response[arms[index].name] = json::object{
          {"final_local", std::move(final_local)},
          {"line_result", std::move(line)},
          {"matches", completed[index].march.matches.size()},
          {"tiles", completed[index].tile_lines.size()},
          {"elapsed_ms", completed[index].elapsed_ms}};
    }
    auto combined_summary = combined->summary();
    combined_summary["session"] = session->handle;
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"capability", kRetainedParallelArmCapability},
        {"native_retained", true}, {"json_coefficients", 0},
        {"atomic_publication", true},
        {"workers", 2},
        {"max_parallel_arms", max_active_workers.load()},
        {"worker_overlap", max_active_workers.load() == 2},
        {"checkpoint_policy", checkpoint_policy},
        {"epsilon", raw_epsilon},
        {"arms", std::move(arm_response)},
        {"combined_line_result", std::move(combined_summary)},
        {"elapsed_ms", std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - operation_started).count()}};
  }

  if (operation == "tile.endpoint_limit") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "tile_plan",
         "tile_plan_checkpoint_identity", "arm", "local",
         "source_checkpoint_identity", "checkpoint_identity",
         "cancellation"},
        "native tile.endpoint_limit request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native plan-bound endpoint evaluation requires rational or Acb coefficients");
    const auto plan_handle = required_string(root, "tile_plan");
    const auto local_handle = required_string(root, "local");
    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> local;
    std::string endpoint_handle;
    {
      // Admission acquires strong ownership of both immutable dependencies.
      // Releasing either public token after this point cannot invalidate the
      // retained endpoint result or its checkpoint closure.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto plan_found = session->tile_plans.find(plan_handle);
      if (plan_found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or released tile plan for plan-bound endpoint limit");
      const auto local_found = session->locals.find(local_handle);
      if (local_found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released final local for plan-bound endpoint limit");
      if (session->endpoints.size() + session->pending_endpoint_limits >=
          session->endpoint_capacity)
        throw std::invalid_argument(
            "persistent endpoint result capacity is exhausted");
      plan = plan_found->second;
      local = local_found->second;
      endpoint_handle = "e:" +
          std::to_string(session->next_endpoint++);
      ++session->pending_endpoint_limits;
    }

    std::shared_ptr<StoredEndpointResult> endpoint;
    try {
      endpoint = build_planned_endpoint_limit(
          endpoint_handle, root, plan, local);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_endpoint_limits == 0)
        throw std::logic_error(
            "native plan-bound endpoint reservation accounting underflow");
      --session->pending_endpoint_limits;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_endpoint_limits == 0)
        throw std::logic_error(
            "native plan-bound endpoint reservation accounting underflow");
      --session->pending_endpoint_limits;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during plan-bound endpoint limit");
      session->endpoints.emplace(endpoint_handle, endpoint);
      ++session->total_endpoint_limits;
      session->total_endpoint_limit_ms += endpoint->elapsed_ms();
    }
    auto response = endpoint->summary();
    response["status"] = "ok";
    response["session"] = session->handle;
    return response;
  }

  if (operation == "integration.line") {
    if (root.if_contains("certify_tail") != nullptr)
      require_exact_keys(root,
          {"schema", "op", "session", "tile_plan", "local", "arm",
           "tile", "epsilon", "source_checkpoint_identity",
           "tile_plan_checkpoint_identity", "checkpoint_identity",
           "certify_tail"},
          "native integration.line request");
    else
      require_exact_keys(root,
          {"schema", "op", "session", "tile_plan", "local", "arm",
           "tile", "epsilon", "source_checkpoint_identity",
           "tile_plan_checkpoint_identity", "checkpoint_identity"},
          "native integration.line request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native planned line integration requires rational or Acb coefficients");
    const auto plan_handle = required_string(root, "tile_plan");
    const auto local_handle = required_string(root, "local");
    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> local;
    std::string line_handle;
    {
      // The admitted operation strongly owns both dependencies.  Lower and
      // upper calls take this lock only for admission/publication and execute
      // their immutable plan arms concurrently outside it.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto plan_found = session->tile_plans.find(plan_handle);
      if (plan_found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or released tile plan for line integration");
      const auto local_found = session->locals.find(local_handle);
      if (local_found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released local for line integration");
      if (session->line_results.size() +
              session->pending_line_integrations >=
          session->line_result_capacity)
        throw std::invalid_argument(
            "persistent native line-result capacity is exhausted");
      plan = plan_found->second;
      local = local_found->second;
      line_handle = "line:" +
          std::to_string(session->next_line_result++);
      ++session->pending_line_integrations;
    }
    std::shared_ptr<StoredLineResult> result;
    try {
      result = build_planned_line_result(
          line_handle, root, plan, local);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_line_integrations == 0)
        throw std::logic_error(
            "native line-integration reservation accounting underflow");
      --session->pending_line_integrations;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_line_integrations == 0)
        throw std::logic_error(
            "native line-integration reservation accounting underflow");
      --session->pending_line_integrations;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during native line integration");
      session->line_results.emplace(line_handle, result);
      ++session->total_line_integrations;
      session->total_line_integration_ms += result->elapsed_ms();
      plan->note_integration();
    }
    auto response = result->summary();
    response["status"] = "ok";
    response["session"] = session->handle;
    return response;
  }

  if (operation == "integration.stats" ||
      operation == "integration.export") {
    const auto line_handle = required_string(root, "line");
    std::shared_ptr<StoredLineResult> result;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->line_results.find(line_handle);
      if (found == session->line_results.end())
        throw std::invalid_argument(
            "unknown or released native line result");
      result = found->second;
    }
    json::object response;
    if (operation == "integration.stats") {
      response = result->stats_json();
    } else {
      const auto output_digits = root.if_contains("output_digits")
          ? static_cast<int>(as_i64(root.at("output_digits"),
                                    "line export output digits"))
          : session->output_digits;
      if (output_digits < 1)
        throw std::invalid_argument(
            "line export output digits must be positive");
      response = result->export_values(
          required_string(root, "checkpoint_identity"), output_digits);
      const auto elapsed = response.at("elapsed_ms").as_double();
      std::lock_guard<std::mutex> lock(session->mutex);
      ++session->total_line_exports;
      session->total_line_export_ms += elapsed;
    }
    response["status"] = "ok";
    response["session"] = session->handle;
    return response;
  }

  if (operation == "integration.release") {
    const auto line_handle = required_string(root, "line");
    std::shared_ptr<StoredLineResult> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->line_results.find(line_handle);
      if (found == session->line_results.end())
        throw std::invalid_argument(
            "unknown or already released native line result");
      removed = std::move(found->second);
      session->line_results.erase(found);
    }
    return json::object{
        {"status", "ok"}, {"released", line_handle},
        {"checkpoint_identity", removed->checkpoint_identity()}};
  }

  if (operation == "chart.prepare") {
    const auto key = required_string(root, "key");
    const auto identity = required_string(root, "identity");
    const auto& problem = as_object(root.at("problem"), "prepared problem");
    const auto problem_domain = required_string(problem, "domain");
    if (problem_domain != session->domain)
      throw std::invalid_argument(
          "prepared problem domain differs from its solver session");
    if (session->domain == "acb" &&
        as_i64(problem.at("precision_bits"), "precision bits") !=
            session->precision_bits)
      throw std::invalid_argument(
          "prepared problem precision differs from its solver session");
    if (session->domain == "symbolic" &&
        parse_symbols(problem) != session->symbols)
      throw std::invalid_argument(
          "prepared problem regulator field differs from its solver session");
    const auto analytic = root.if_contains("analytic")
        ? root.at("analytic") : json::value(nullptr);
    std::optional<std::string> geometry_record;
    std::optional<std::string> principal_matrix_record;
    std::optional<std::string> native_scc_capabilities;
    if (analytic.is_object()) {
      const auto& analytic_object = analytic.as_object();
      if (const auto* geometry = analytic_object.if_contains("geometry"))
        geometry_record = canonical_chart_geometry_record(*geometry);
      if (const auto* principal =
              analytic_object.if_contains("principal_matrix"))
        principal_matrix_record = parse_exact_parent_matrix(
            *principal, as_u32(problem.at("d"), "dimension"),
            "prepared chart principal matrix").canonical_record;
      if (const auto* capabilities =
              analytic_object.if_contains("native_scc_capabilities"))
        native_scc_capabilities =
            canonical_native_scc_capabilities(*capabilities);
    }
    json::object combined_analytic;
    combined_analytic["session"] = json::parse(session->analytic_identity);
    combined_analytic["chart"] = analytic;
    auto scc = validate_scc_certificate(
        root.at("scc"), as_u32(problem.at("d"), "dimension"));
    auto signature = static_problem_signature(
        problem, combined_analytic, scc, identity);

    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (const auto found = session->handles_by_key.find(key);
          found != session->handles_by_key.end()) {
        const auto chart = session->charts.at(found->second);
        if (chart->signature() != signature)
          throw std::invalid_argument(
              "persistent chart cache key collision with unequal exact identity");
        return json::object{{"status", "ok"}, {"chart", chart->handle()},
                            {"reused", true},
                            {"dimension", chart->dimension()},
                            {"frame_base", chart->frame_base()},
                            {"frame_width", chart->frame_width()},
                            {"d0_inverse_mode", chart->d0_inverse_mode()}};
      }
      if (session->charts.size() >= session->chart_capacity)
        throw std::invalid_argument("persistent chart capacity is exhausted");
    }
    std::string chart_handle;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      chart_handle = "c:" + std::to_string(session->next_chart++);
    }
    std::shared_ptr<PreparedChartBase> chart;
    if (session->domain == "rational")
      chart = parse_prepared_chart<Rational>(
          session, root, chart_handle, key, identity, geometry_record,
          principal_matrix_record, native_scc_capabilities,
          std::move(scc), std::move(signature));
    else if (session->domain == "acb")
      chart = parse_prepared_chart<ComplexBall>(
          session, root, chart_handle, key, identity, geometry_record,
          principal_matrix_record, native_scc_capabilities,
          std::move(scc), std::move(signature));
    else
      chart = parse_prepared_chart<SymbolicRational>(
          session, root, chart_handle, key, identity, geometry_record,
          principal_matrix_record, native_scc_capabilities,
          std::move(scc), std::move(signature));
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      // A concurrent duplicate prepare is harmless only when its complete
      // collision certificate is byte-identical.
      if (const auto found = session->handles_by_key.find(key);
          found != session->handles_by_key.end()) {
        const auto existing = session->charts.at(found->second);
        if (existing->signature() != chart->signature())
          throw std::invalid_argument(
              "concurrent chart cache key collision with unequal identity");
        chart = existing;
      } else {
        session->charts.emplace(chart->handle(), chart);
        session->handles_by_key.emplace(key, chart->handle());
      }
    }
    return json::object{{"status", "ok"}, {"chart", chart->handle()},
                        {"reused", chart->handle() != chart_handle},
                        {"dimension", chart->dimension()},
                        {"frame_base", chart->frame_base()},
                        {"frame_width", chart->frame_width()},
                        {"d0_inverse_mode", chart->d0_inverse_mode()},
                        {"scc_components", chart->scc().component_count},
                        {"scc_structural_edges",
                         chart->scc().structural_edges.size()},
                        {"scc_condensation_edges",
                         chart->scc().condensation_edges.size()},
                        {"scc_topological_order",
                         encode_indices(chart->scc().topological_order)},
                        {"scc_coupling_depth", chart->scc().coupling_depth}};
  }

  if (operation == "scc.prepare") {
    const auto key = required_string(root, "key");
    const auto identity = required_string(root, "identity");
    const auto signature = composite_scc_signature(root);
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (const auto found = session->scc_handles_by_key.find(key);
          found != session->scc_handles_by_key.end()) {
        const auto composite = session->sccs.at(found->second);
        if (composite->signature() != signature)
          throw std::invalid_argument(
              "persistent SCC cache key collision with unequal exact identity");
        auto result = composite->stats_json();
        result["status"] = "ok";
        result["session"] = session->handle;
        result["reused"] = true;
        return result;
      }
      if (session->sccs.size() >= session->scc_capacity)
        throw std::invalid_argument("persistent SCC capacity is exhausted");
    }

    const auto& raw_blocks = as_array(root.at("blocks"), "SCC blocks");
    std::vector<std::shared_ptr<PreparedChartBase>> erased_charts;
    erased_charts.reserve(raw_blocks.size());
    std::string scc_handle;
    {
      // Resolve and strongly retain the complete diagonal handle set before
      // leaving the session lock.  The composite remains valid if the public
      // chart handles are released after this preparation boundary.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      for (const auto& raw_block_value : raw_blocks) {
        const auto& raw_block = as_object(raw_block_value, "SCC block");
        const auto chart_handle = required_string(raw_block, "chart");
        const auto found = session->charts.find(chart_handle);
        if (found == session->charts.end())
          throw std::invalid_argument(
              "SCC preparation references an unknown or released chart");
        erased_charts.push_back(found->second);
      }
      scc_handle = "scc:" + std::to_string(session->next_scc++);
    }

    std::shared_ptr<CompositeSCCChartBase> composite;
    if (session->domain == "rational")
      composite = parse_composite_scc_chart<Rational>(
          session, root, scc_handle, key, identity, signature,
          erased_charts);
    else if (session->domain == "acb")
      composite = parse_composite_scc_chart<ComplexBall>(
          session, root, scc_handle, key, identity, signature,
          erased_charts);
    else
      composite = parse_composite_scc_chart<SymbolicRational>(
          session, root, scc_handle, key, identity, signature,
          erased_charts);

    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during SCC preparation");
      if (const auto found = session->scc_handles_by_key.find(key);
          found != session->scc_handles_by_key.end()) {
        const auto existing = session->sccs.at(found->second);
        if (existing->signature() != composite->signature())
          throw std::invalid_argument(
              "concurrent SCC cache key collision with unequal exact identity");
        composite = existing;
      } else {
        if (session->sccs.size() >= session->scc_capacity)
          throw std::invalid_argument(
              "persistent SCC capacity was exhausted during preparation");
        session->sccs.emplace(composite->handle(), composite);
        session->scc_handles_by_key.emplace(key, composite->handle());
      }
    }
    auto result = composite->stats_json();
    result["status"] = "ok";
    result["session"] = session->handle;
    result["reused"] = composite->handle() != scc_handle;
    return result;
  }

  if (operation == "chart.solve") {
    const auto chart_handle = required_string(root, "chart");
    std::shared_ptr<PreparedChartBase> chart;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->charts.find(chart_handle);
      if (found == session->charts.end())
        throw std::invalid_argument("unknown or released persistent chart");
      chart = found->second;
    }
    const auto output_digits = root.if_contains("output_digits")
        ? static_cast<int>(as_i64(root.at("output_digits"), "output digits"))
        : session->output_digits;
    if (output_digits < 1)
      throw std::invalid_argument("output digits must be positive");
    auto result = chart->solve(as_object(root.at("run"), "recurrence run"),
                               output_digits);
    result["session"] = session->handle;
    result["chart"] = chart->handle();
    return result;
  }

  if (operation == "chart.solve_batch") {
    const auto chart_handle = required_string(root, "chart");
    std::shared_ptr<PreparedChartBase> chart;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->charts.find(chart_handle);
      if (found == session->charts.end())
        throw std::invalid_argument("unknown or released persistent chart");
      chart = found->second;
    }
    const auto output_digits = root.if_contains("output_digits")
        ? static_cast<int>(as_i64(root.at("output_digits"), "output digits"))
        : session->output_digits;
    if (output_digits < 1)
      throw std::invalid_argument("output digits must be positive");
    const auto& runs = as_array(root.at("runs"), "persistent recurrence runs");
    const auto requested_threads = root.if_contains("threads")
        ? as_u32(root.at("threads"), "batch threads") : 1;
    if (requested_threads == 0)
      throw std::invalid_argument("batch threads must be positive");

    const bool symbolic_serialized = session->domain == "symbolic";
    const auto bounded_threads = std::min<std::size_t>(
        {static_cast<std::size_t>(requested_threads),
         static_cast<std::size_t>(kMaxPersistentBatchThreads), runs.size()});
    const auto worker_count = symbolic_serialized && bounded_threads != 0
        ? std::size_t{1} : bounded_threads;
    const auto started = std::chrono::steady_clock::now();
    std::vector<json::object> results(runs.size());
    std::atomic<std::size_t> next{0};
    auto worker = [&]() {
      while (true) {
        const auto index = next.fetch_add(1);
        if (index >= runs.size()) return;
        results[index] = solve_prepared_chart_safe(
            chart, runs[index], output_digits, session->handle);
      }
    };
    // jthread guarantees already-started workers are joined if a later
    // thread construction throws; destroying a joinable std::thread here
    // would otherwise terminate the host Wolfram kernel.
    std::vector<std::jthread> workers;
    workers.reserve(worker_count);
    for (std::size_t i = 0; i < worker_count; ++i)
      workers.emplace_back(worker);
    for (auto& thread : workers) thread.join();

    std::size_t succeeded = 0;
    json::array encoded;
    encoded.reserve(results.size());
    for (auto& result : results) {
      if (result.if_contains("status") != nullptr &&
          result.at("status") == "ok")
        ++succeeded;
      encoded.push_back(std::move(result));
    }
    const auto elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"chart", chart->handle()}, {"results", std::move(encoded)},
        {"attempted", runs.size()}, {"succeeded", succeeded},
        {"failed", runs.size() - succeeded},
        {"requested_threads", requested_threads},
        {"worker_threads", worker_count},
        {"thread_limit", kMaxPersistentBatchThreads},
        {"symbolic_serialized", symbolic_serialized},
        {"elapsed_ms", elapsed_ms}};
  }

  if (operation == "session.solve_many") {
    const auto& raw_jobs = as_array(
        root.at("jobs"), "persistent cross-chart recurrence jobs");
    const auto requested_threads = root.if_contains("threads")
        ? as_u32(root.at("threads"), "batch threads") : 1;
    if (requested_threads == 0)
      throw std::invalid_argument("batch threads must be positive");

    struct PendingJob {
      std::string chart_handle;
      const json::value* run = nullptr;
      int output_digits = 0;
    };
    std::vector<PendingJob> pending;
    pending.reserve(raw_jobs.size());
    for (std::size_t index = 0; index < raw_jobs.size(); ++index) {
      const auto& job = as_object(
          raw_jobs[index], "persistent cross-chart recurrence job");
      const auto* raw_run = job.if_contains("run");
      if (raw_run == nullptr)
        throw std::invalid_argument(
            "session.solve_many job " + std::to_string(index) +
            " is missing its complete run record");
      const auto output_digits = job.if_contains("output_digits")
          ? static_cast<int>(as_i64(
                job.at("output_digits"), "job output digits"))
          : session->output_digits;
      if (output_digits < 1)
        throw std::invalid_argument("job output digits must be positive");
      pending.push_back(PendingJob{
          required_string(job, "chart"), raw_run, output_digits});
    }

    struct ResolvedJob {
      std::shared_ptr<PreparedChartBase> chart;
      const json::value* run = nullptr;
      int output_digits = 0;
    };
    std::vector<ResolvedJob> jobs;
    jobs.reserve(pending.size());
    {
      // Resolve the complete handle set before any worker exists.  Apart
      // from making cross-session/released handles loud, the shared_ptrs
      // retain every selected chart for the whole batch even if a concurrent
      // chart.release follows this validation boundary.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      for (std::size_t index = 0; index < pending.size(); ++index) {
        const auto found = session->charts.find(pending[index].chart_handle);
        if (found == session->charts.end())
          throw std::invalid_argument(
              "unknown or released persistent chart in session.solve_many "
              "job " + std::to_string(index) + ": " +
              pending[index].chart_handle);
        jobs.push_back(ResolvedJob{
            found->second, pending[index].run, pending[index].output_digits});
      }
    }

    const bool symbolic_serialized = session->domain == "symbolic";
    const auto bounded_threads = std::min<std::size_t>(
        {static_cast<std::size_t>(requested_threads),
         static_cast<std::size_t>(kMaxPersistentBatchThreads), jobs.size()});
    const auto worker_count = symbolic_serialized && bounded_threads != 0
        ? std::size_t{1} : bounded_threads;
    const auto started = std::chrono::steady_clock::now();
    std::vector<json::object> results(jobs.size());
    std::atomic<std::size_t> next{0};
    auto worker = [&](std::stop_token stop) {
      while (!stop.stop_requested()) {
        const auto index = next.fetch_add(1);
        if (index >= jobs.size()) return;
        const auto& job = jobs[index];
        results[index] = solve_prepared_chart_safe(
            job.chart, *job.run, job.output_digits, session->handle);
      }
    };
    // If construction of a later worker fails, jthread destruction requests
    // stop and joins every worker that already started.  No joinable native
    // thread can escape into the host Wolfram kernel.
    std::vector<std::jthread> workers;
    workers.reserve(worker_count);
    for (std::size_t i = 0; i < worker_count; ++i)
      workers.emplace_back(worker);
    for (auto& thread : workers) thread.join();

    std::size_t succeeded = 0;
    json::array encoded;
    encoded.reserve(results.size());
    for (auto& result : results) {
      if (result.if_contains("status") != nullptr &&
          result.at("status") == "ok")
        ++succeeded;
      encoded.push_back(std::move(result));
    }
    const auto elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"results", std::move(encoded)}, {"attempted", jobs.size()},
        {"succeeded", succeeded}, {"failed", jobs.size() - succeeded},
        {"requested_threads", requested_threads},
        {"worker_threads", worker_count},
        {"thread_limit", kMaxPersistentBatchThreads},
        {"symbolic_serialized", symbolic_serialized},
        {"elapsed_ms", elapsed_ms}};
  }

  if (operation == "scc.solve_column") {
    const auto scc_handle = required_string(root, "scc");
    std::shared_ptr<CompositeSCCChartBase> composite;
    std::string local_handle;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->sccs.find(scc_handle);
      if (found == session->sccs.end())
        throw std::invalid_argument(
            "unknown or released persistent SCC chart");
      if (session->locals.size() + session->pending_local_solves >=
          session->local_capacity)
        throw std::invalid_argument("persistent local capacity is exhausted");
      composite = found->second;
      local_handle = "l:" + std::to_string(session->next_local++);
      ++session->pending_local_solves;
    }

    CompositeColumnSolveResult column;
    try {
      column = composite->solve_column(local_handle, root);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "native SCC local reservation accounting underflow");
      --session->pending_local_solves;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "native SCC local reservation accounting underflow");
      --session->pending_local_solves;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during SCC column solve");
      session->locals.emplace(local_handle, column.local);
      const auto local_stats = column.local->stats();
      ++session->total_local_solves;
      ++session->total_scc_column_solves;
      session->total_local_run_parse_ms += local_stats.create_parse_ms;
      session->total_local_kernel_ms += local_stats.create_kernel_ms;
    }
    auto response = column.local->summary();
    response["status"] = "ok";
    response["session"] = session->handle;
    response["scc"] = composite->handle();
    response["native_retained"] = true;
    response["json_coefficients"] = 0;
    response["execution_capability"] =
        composite->column_execution_capability();
    response["block_diagnostics"] = std::move(column.block_diagnostics);
    response["elapsed_ms"] = column.elapsed_ms;
    return response;
  }

  if (operation == "scc.solve_columns") {
    require_exact_keys(root,
        {"schema", "op", "session", "scc", "columns", "threads"},
        "native SCC column-batch request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native SCC column batches require rational or Acb coefficients");
    const auto requested_threads = as_u32(
        root.at("threads"), "native SCC column-batch threads");
    if (requested_threads == 0)
      throw std::invalid_argument(
          "native SCC column-batch threads must be positive");
    const auto& raw_columns = as_array(
        root.at("columns"), "native SCC column-batch columns");
    if (raw_columns.empty())
      throw std::invalid_argument(
          "native SCC column batch cannot be empty");
    for (std::size_t index = 0; index < raw_columns.size(); ++index)
      require_exact_keys(
          as_object(raw_columns[index], "native SCC batch column"),
          {"seed", "targets", "checkpoint_identity"},
          "native SCC batch column");

    const auto scc_handle = required_string(root, "scc");
    std::shared_ptr<CompositeSCCChartBase> composite;
    std::vector<std::string> local_handles;
    local_handles.reserve(raw_columns.size());
    {
      // Reserve the complete ordered result set before starting workers.
      // Public SCC release cannot invalidate the strongly owned composite,
      // and capacity cannot be consumed between individual columns.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto found = session->sccs.find(scc_handle);
      if (found == session->sccs.end())
        throw std::invalid_argument(
            "unknown or released persistent SCC chart");
      if (raw_columns.size() > session->local_capacity ||
          session->locals.size() + session->pending_local_solves >
              session->local_capacity - raw_columns.size())
        throw std::invalid_argument(
            "persistent local capacity cannot admit the complete SCC column batch");
      composite = found->second;
      for (std::size_t index = 0; index < raw_columns.size(); ++index)
        local_handles.push_back(
            "l:" + std::to_string(session->next_local++));
      session->pending_local_solves += raw_columns.size();
    }

    const auto worker_count = std::min<std::size_t>(
        {static_cast<std::size_t>(requested_threads),
         static_cast<std::size_t>(kMaxPersistentBatchThreads),
         raw_columns.size()});
    const auto started = std::chrono::steady_clock::now();
    std::vector<CompositeColumnSolveResult> columns(raw_columns.size());
    std::vector<std::exception_ptr> errors(raw_columns.size());
    std::atomic<std::size_t> next{0};
    auto worker = [&](std::stop_token stop) {
      while (!stop.stop_requested()) {
        const auto index = next.fetch_add(1);
        if (index >= raw_columns.size()) return;
        try {
          columns[index] = composite->solve_column(
              local_handles[index],
              as_object(raw_columns[index], "native SCC batch column"));
        } catch (...) {
          errors[index] = std::current_exception();
        }
      }
    };
    std::vector<std::jthread> workers;
    workers.reserve(worker_count);
    try {
      for (std::size_t index = 0; index < worker_count; ++index)
        workers.emplace_back(worker);
      for (auto& thread : workers) thread.join();
    } catch (...) {
      for (auto& thread : workers) thread.request_stop();
      for (auto& thread : workers)
        if (thread.joinable()) thread.join();
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves < raw_columns.size())
        throw std::logic_error(
            "native SCC column-batch reservation accounting underflow");
      session->pending_local_solves -= raw_columns.size();
      throw;
    }

    const auto failed = std::find_if(
        errors.begin(), errors.end(), [](const auto& error) {
          return error != nullptr;
        });
    if (failed != errors.end()) {
      const auto failure = *failed;
      {
        std::lock_guard<std::mutex> lock(session->mutex);
        if (session->pending_local_solves < raw_columns.size())
          throw std::logic_error(
              "native SCC column-batch reservation accounting underflow");
        session->pending_local_solves -= raw_columns.size();
      }
      // A batch is atomic at the retained-state boundary: successful worker
      // temporaries are discarded and the first ordered failure is loud.
      std::rethrow_exception(failure);
    }

    if (std::any_of(columns.begin(), columns.end(), [](const auto& column) {
          return column.local == nullptr;
        })) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves < raw_columns.size())
        throw std::logic_error(
            "native SCC column-batch reservation accounting underflow");
      session->pending_local_solves -= raw_columns.size();
      throw std::logic_error(
          "native SCC column batch completed without a retained local");
    }

    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves < raw_columns.size())
        throw std::logic_error(
            "native SCC column-batch reservation accounting underflow");
      session->pending_local_solves -= raw_columns.size();
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during SCC column batch");
      session->locals.reserve(session->locals.size() + columns.size());
      std::vector<std::string> inserted;
      inserted.reserve(columns.size());
      try {
        for (std::size_t index = 0; index < columns.size(); ++index) {
          if (!session->locals.emplace(
                  local_handles[index], columns[index].local).second)
            throw std::logic_error(
                "native SCC column batch produced a duplicate local handle");
          inserted.push_back(local_handles[index]);
        }
      } catch (...) {
        for (const auto& handle : inserted) session->locals.erase(handle);
        throw;
      }
      for (std::size_t index = 0; index < columns.size(); ++index) {
        const auto local_stats = columns[index].local->stats();
        ++session->total_local_solves;
        ++session->total_scc_column_solves;
        session->total_local_run_parse_ms += local_stats.create_parse_ms;
        session->total_local_kernel_ms += local_stats.create_kernel_ms;
      }
    }

    json::array responses;
    responses.reserve(columns.size());
    for (auto& column : columns) {
      auto response = column.local->summary();
      response["status"] = "ok";
      response["session"] = session->handle;
      response["scc"] = composite->handle();
      response["native_retained"] = true;
      response["json_coefficients"] = 0;
      response["execution_capability"] =
          composite->column_execution_capability();
      response["block_diagnostics"] =
          std::move(column.block_diagnostics);
      response["elapsed_ms"] = column.elapsed_ms;
      responses.push_back(std::move(response));
    }
    const auto elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"scc", composite->handle()},
        {"results", std::move(responses)},
        {"columns", raw_columns.size()},
        {"requested_threads", requested_threads},
        {"worker_threads", worker_count},
        {"thread_limit", kMaxPersistentBatchThreads},
        {"atomic_retention", true},
        {"json_coefficients", 0}, {"elapsed_ms", elapsed_ms}};
  }

  if (operation == "local.solve") {
    const auto chart_handle = required_string(root, "chart");
    std::shared_ptr<PreparedChartBase> chart;
    std::string local_handle;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->charts.find(chart_handle);
      if (found == session->charts.end())
        throw std::invalid_argument("unknown or released persistent chart");
      if (session->locals.size() + session->pending_local_solves >=
          session->local_capacity)
        throw std::invalid_argument("persistent local capacity is exhausted");
      chart = found->second;
      local_handle = "l:" + std::to_string(session->next_local++);
      ++session->pending_local_solves;
    }

    std::shared_ptr<StoredLocalBase> local;
    try {
      local = chart->solve_local(
          local_handle, as_object(root.at("run"), "recurrence run"),
          as_object(root.at("metadata"), "local metadata"));
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error("native local reservation accounting underflow");
      --session->pending_local_solves;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error("native local reservation accounting underflow");
      --session->pending_local_solves;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during local solve");
      session->locals.emplace(local_handle, local);
      const auto local_stats = local->stats();
      ++session->total_local_solves;
      session->total_local_run_parse_ms += local_stats.create_parse_ms;
      session->total_local_kernel_ms += local_stats.create_kernel_ms;
    }
    auto result = local->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    result["native_retained"] = true;
    result["json_coefficients"] = 0;
    return result;
  }

  if (operation == "local.apply_rational_row") {
    require_exact_keys(root,
        {"schema", "op", "session", "local", "row",
         "source_checkpoint_identity", "checkpoint_identity"},
        "native local.apply_rational_row request");
    if (session->domain == "symbolic")
      throw std::domain_error(
          "native rational-row application requires exact Rational or specialized Acb coefficients");
    const auto source_handle = required_string(root, "local");
    std::shared_ptr<StoredLocalBase> source;
    std::string local_handle;
    {
      // Materialization runs outside the registry lock, but owns the source
      // for its full lifetime.  The derived scalar local also retains that
      // source as provenance after publication.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto found = session->locals.find(source_handle);
      if (found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released retained local for rational-row application");
      if (session->locals.size() + session->pending_local_solves >=
          session->local_capacity)
        throw std::invalid_argument("persistent local capacity is exhausted");
      source = found->second;
      local_handle = "l:" + std::to_string(session->next_local++);
      ++session->pending_local_solves;
    }

    std::shared_ptr<StoredLocalBase> local;
    try {
      if (session->domain == "rational") {
        const auto typed =
            std::dynamic_pointer_cast<StoredLocal<Rational>>(source);
        if (!typed)
          throw std::logic_error(
              "rational-row source local differs from its Rational session");
        local = build_rational_row_local<Rational>(
            local_handle, root, session->precision_bits, typed, source);
      } else {
        const auto typed =
            std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(source);
        if (!typed)
          throw std::logic_error(
              "rational-row source local differs from its Acb session");
        local = build_rational_row_local<ComplexBall>(
            local_handle, root, session->precision_bits, typed, source);
      }
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "native rational-row reservation accounting underflow");
      --session->pending_local_solves;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "native rational-row reservation accounting underflow");
      --session->pending_local_solves;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during rational-row application");
      if (!session->locals.emplace(local_handle, local).second)
        throw std::logic_error(
            "native rational-row application produced a duplicate local handle");
    }
    auto result = local->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    result["source_local"] = source_handle;
    result["application_capability"] = kRetainedRationalRowCapability;
    result["native_retained"] = true;
    result["json_coefficients"] = 0;
    return result;
  }

  if (operation == "local.endpoint_limit") {
    if (root.if_contains("output_digits") != nullptr ||
        root.if_contains("include_coefficients") != nullptr)
      throw std::invalid_argument(
          "local.endpoint_limit never exports coefficients; use the explicit "
          "endpoint.export compatibility operation");
    const auto local_handle = required_string(root, "local");
    std::shared_ptr<StoredLocalBase> local;
    std::string endpoint_handle;
    {
      // Strong ownership is acquired under the session lock.  Public
      // local.release may remove the registry token after admission without
      // invalidating this endpoint computation.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto found = session->locals.find(local_handle);
      if (found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released native local for endpoint limit");
      if (session->endpoints.size() + session->pending_endpoint_limits >=
          session->endpoint_capacity)
        throw std::invalid_argument(
            "persistent endpoint result capacity is exhausted");
      local = found->second;
      endpoint_handle = "e:" + std::to_string(session->next_endpoint++);
      ++session->pending_endpoint_limits;
    }

    std::shared_ptr<StoredEndpointResult> endpoint;
    try {
      endpoint = build_endpoint_limit(endpoint_handle, root, local);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_endpoint_limits == 0)
        throw std::logic_error(
            "native endpoint reservation accounting underflow");
      --session->pending_endpoint_limits;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_endpoint_limits == 0)
        throw std::logic_error(
            "native endpoint reservation accounting underflow");
      --session->pending_endpoint_limits;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during endpoint limit");
      session->endpoints.emplace(endpoint_handle, endpoint);
      ++session->total_endpoint_limits;
      session->total_endpoint_limit_ms += endpoint->elapsed_ms();
    }
    auto result = endpoint->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "local.match") {
    if (session->domain != "rational")
      throw std::invalid_argument(
          "local.match currently supports only the exact rational regular "
          "local capability; Acb and symbolic saturation are not routed "
          "through it");
    const auto& raw_basis = as_array(
        root.at("basis"), "exact regular local match basis");
    if (raw_basis.empty())
      throw std::invalid_argument("local.match basis cannot be empty");
    std::vector<std::string> basis_handles;
    basis_handles.reserve(raw_basis.size());
    std::set<std::string> unique_handles;
    for (const auto& raw_handle : raw_basis) {
      if (!raw_handle.is_string())
        throw std::invalid_argument(
            "local.match basis handles must be strings");
      std::string handle(raw_handle.as_string());
      if (!unique_handles.insert(handle).second)
        throw std::invalid_argument(
            "local.match basis handles must be pairwise distinct");
      basis_handles.push_back(std::move(handle));
    }
    const auto incoming_handle = required_string(root, "incoming");
    if (unique_handles.contains(incoming_handle))
      throw std::invalid_argument(
          "local.match incoming handle must be distinct from its basis");

    std::vector<std::shared_ptr<StoredLocalBase>> basis;
    std::shared_ptr<StoredLocalBase> incoming;
    std::string match_handle;
    {
      // Resolve and strongly retain every local before releasing the session
      // lock.  A concurrent public local.release cannot invalidate an
      // already admitted native match operation.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      if (session->matches.size() + session->pending_matches >=
          session->match_capacity)
        throw std::invalid_argument(
            "persistent local match capacity is exhausted");
      for (const auto& handle : basis_handles) {
        const auto found = session->locals.find(handle);
        if (found == session->locals.end())
          throw std::invalid_argument(
              "unknown or released native local in match basis: " + handle);
        basis.push_back(found->second);
      }
      const auto found = session->locals.find(incoming_handle);
      if (found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released incoming native local: " +
            incoming_handle);
      incoming = found->second;
      match_handle = "m:" + std::to_string(session->next_match++);
      ++session->pending_matches;
    }

    std::shared_ptr<StoredExactRegularMatch> match;
    try {
      match = build_exact_regular_match(
          match_handle, root, basis_handles, basis, incoming_handle,
          incoming);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches == 0)
        throw std::logic_error(
            "native local match reservation accounting underflow");
      --session->pending_matches;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches == 0)
        throw std::logic_error(
            "native local match reservation accounting underflow");
      --session->pending_matches;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during local matching");
      session->matches.emplace(match_handle, match);
      ++session->total_local_matches;
      session->total_local_match_ms += match->elapsed_ms();
    }
    auto result = match->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "local.match_acb") {
    if (session->domain != "acb")
      throw std::invalid_argument(
          "local.match_acb requires an Acb persistent session; exact rational matching remains local.match v1");
    const auto& raw_basis = as_array(
        root.at("basis"), "refined Acb local match basis");
    if (raw_basis.empty())
      throw std::invalid_argument(
          "local.match_acb basis cannot be empty");
    std::vector<std::string> basis_handles;
    basis_handles.reserve(raw_basis.size());
    std::set<std::string> unique_handles;
    for (const auto& raw_handle : raw_basis) {
      if (!raw_handle.is_string())
        throw std::invalid_argument(
            "local.match_acb basis handles must be strings");
      std::string handle(raw_handle.as_string());
      if (!unique_handles.insert(handle).second)
        throw std::invalid_argument(
            "local.match_acb basis handles must be pairwise distinct");
      basis_handles.push_back(std::move(handle));
    }
    const auto incoming_handle = required_string(root, "incoming");
    if (unique_handles.contains(incoming_handle))
      throw std::invalid_argument(
          "local.match_acb incoming handle must be distinct from its basis");

    std::vector<std::shared_ptr<StoredLocalBase>> basis;
    std::shared_ptr<StoredLocalBase> incoming;
    std::string match_handle;
    {
      // Admission takes strong ownership of every source before releasing
      // the session lock.  Concurrent public releases cannot invalidate the
      // exact-point evaluation or the bounded refinement operation.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      if (session->matches.size() + session->pending_matches >=
          session->match_capacity)
        throw std::invalid_argument(
            "persistent local match capacity is exhausted");
      for (const auto& handle : basis_handles) {
        const auto found = session->locals.find(handle);
        if (found == session->locals.end())
          throw std::invalid_argument(
              "unknown or released native local in Acb match basis: " +
              handle);
        basis.push_back(found->second);
      }
      const auto found = session->locals.find(incoming_handle);
      if (found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released incoming native Acb local: " +
            incoming_handle);
      incoming = found->second;
      match_handle = "m:" + std::to_string(session->next_match++);
      ++session->pending_matches;
    }

    std::shared_ptr<StoredRefinedAcbMatch> match;
    try {
      match = build_refined_acb_match(
          match_handle, root, basis_handles, basis, incoming_handle,
          incoming, session->precision_bits,
          checkpoint_configuration_identity(*session));
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches == 0)
        throw std::logic_error(
            "native Acb local match reservation accounting underflow");
      --session->pending_matches;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches == 0)
        throw std::logic_error(
            "native Acb local match reservation accounting underflow");
      --session->pending_matches;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during Acb local matching");
      session->matches.emplace(match_handle, match);
      ++session->total_local_matches;
      session->total_local_match_ms += match->elapsed_ms();
    }
    auto result = match->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "match.materialize_local") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "match", "checkpoint_identity"},
        "native match.materialize_local request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native plan-match local materialization requires rational or Acb coefficients");
    const auto match_handle = required_string(root, "match");
    const auto checkpoint_identity = required_string(
        root, "checkpoint_identity");
    std::shared_ptr<StoredPlannedMatchHop> match;
    std::string local_handle;
    {
      // Admission strongly retains the complete handoff before releasing the
      // session lock.  The finite Laurent combination then runs natively and
      // independently of public match/plan/basis tokens.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto found = session->matches.find(match_handle);
      if (found == session->matches.end())
        throw std::invalid_argument(
            "unknown or released retained match for local materialization");
      match = std::dynamic_pointer_cast<StoredPlannedMatchHop>(found->second);
      if (!match)
        throw std::invalid_argument(
            "match.materialize_local requires a plan-driven match handoff");
      if (session->locals.size() + session->pending_local_solves >=
          session->local_capacity)
        throw std::invalid_argument(
            "persistent local capacity is exhausted");
      local_handle = "l:" + std::to_string(session->next_local++);
      ++session->pending_local_solves;
    }

    std::shared_ptr<StoredLocalBase> local;
    try {
      local = match->materialize(
          local_handle, checkpoint_identity, session->precision_bits, match);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "native match materialization reservation accounting underflow");
      --session->pending_local_solves;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "native match materialization reservation accounting underflow");
      --session->pending_local_solves;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during match materialization");
      session->locals.emplace(local_handle, local);
    }
    auto result = local->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    result["materialization_capability"] =
        kRetainedPlannedMatchMaterializationCapability;
    result["native_retained"] = true;
    result["json_coefficients"] = 0;
    return result;
  }

  if (operation == "local.evaluate") {
    const auto local_handle = required_string(root, "local");
    std::shared_ptr<StoredLocalBase> local;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->locals.find(local_handle);
      if (found == session->locals.end())
        throw std::invalid_argument("unknown or released native local solution");
      local = found->second;
    }
    const auto output_digits = root.if_contains("output_digits")
        ? static_cast<int>(as_i64(root.at("output_digits"), "output digits"))
        : session->output_digits;
    if (output_digits < 1)
      throw std::invalid_argument("output digits must be positive");
    auto result = local->evaluate(root, output_digits);
    result["status"] = "ok";
    result["session"] = session->handle;
    result["local"] = local->handle();
    result["chart"] = local->source_chart();
    return result;
  }

  if (operation == "local.certify_residual") {
    const auto local_handle = required_string(root, "local");
    std::shared_ptr<StoredLocalBase> local;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->locals.find(local_handle);
      if (found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released native local solution");
      local = found->second;
    }
    const auto output_digits = root.if_contains("output_digits")
        ? static_cast<int>(as_i64(root.at("output_digits"), "output digits"))
        : session->output_digits;
    if (output_digits < 1)
      throw std::invalid_argument("output digits must be positive");
    auto result = local->certify_residual(root, output_digits);
    result["status"] = "ok";
    result["session"] = session->handle;
    result["local"] = local->handle();
    result["chart"] = local->source_chart();
    return result;
  }

  if (operation == "endpoint.stats") {
    const auto endpoint_handle = required_string(root, "endpoint");
    std::shared_ptr<StoredEndpointResult> endpoint;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->endpoints.find(endpoint_handle);
      if (found == session->endpoints.end())
        throw std::invalid_argument(
            "unknown or released native endpoint result");
      endpoint = found->second;
    }
    auto result = endpoint->stats_json();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "endpoint.export") {
    const auto endpoint_handle = required_string(root, "endpoint");
    std::shared_ptr<StoredEndpointResult> endpoint;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->endpoints.find(endpoint_handle);
      if (found == session->endpoints.end())
        throw std::invalid_argument(
            "unknown or released native endpoint result");
      endpoint = found->second;
    }
    const auto output_digits = root.if_contains("output_digits")
        ? static_cast<int>(as_i64(root.at("output_digits"), "output digits"))
        : session->output_digits;
    if (output_digits < 1)
      throw std::invalid_argument("output digits must be positive");
    auto result = endpoint->export_values(
        required_string(root, "checkpoint_identity"), output_digits);
    const auto export_ms = result.at("elapsed_ms").as_double();
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      ++session->total_endpoint_exports;
      session->total_endpoint_export_ms += export_ms;
    }
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "endpoint.release") {
    const auto endpoint_handle = required_string(root, "endpoint");
    std::shared_ptr<StoredEndpointResult> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->endpoints.find(endpoint_handle);
      if (found == session->endpoints.end())
        throw std::invalid_argument(
            "unknown or already released native endpoint result");
      removed = std::move(found->second);
      session->endpoints.erase(found);
    }
    return json::object{{"status", "ok"}, {"released", endpoint_handle},
                        {"checkpoint_identity",
                         removed->checkpoint_identity()}};
  }

  if (operation == "match.stats") {
    const auto match_handle = required_string(root, "match");
    std::shared_ptr<StoredMatchBase> match;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->matches.find(match_handle);
      if (found == session->matches.end())
        throw std::invalid_argument(
            "unknown or released native local match");
      match = found->second;
    }
    auto result = match->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "match.release") {
    const auto match_handle = required_string(root, "match");
    std::shared_ptr<StoredMatchBase> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->matches.find(match_handle);
      if (found == session->matches.end())
        throw std::invalid_argument(
            "unknown or already released native local match");
      removed = std::move(found->second);
      session->matches.erase(found);
    }
    return json::object{{"status", "ok"}, {"released", match_handle}};
  }

  if (operation == "local.release") {
    const auto local_handle = required_string(root, "local");
    std::shared_ptr<StoredLocalBase> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->locals.find(local_handle);
      if (found == session->locals.end())
        throw std::invalid_argument(
            "unknown or already released native local solution");
      removed = std::move(found->second);
      session->locals.erase(found);
    }
    return json::object{{"status", "ok"}, {"released", local_handle},
                        {"chart", removed->source_chart()}};
  }

  if (operation == "local.stats") {
    const auto local_handle = required_string(root, "local");
    std::shared_ptr<StoredLocalBase> local;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->locals.find(local_handle);
      if (found == session->locals.end())
        throw std::invalid_argument("unknown or released native local solution");
      local = found->second;
    }
    auto result = local->stats_json();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "scc.stats") {
    const auto scc_handle = required_string(root, "scc");
    std::shared_ptr<CompositeSCCChartBase> composite;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->sccs.find(scc_handle);
      if (found == session->sccs.end())
        throw std::invalid_argument(
            "unknown or released persistent SCC chart");
      composite = found->second;
    }
    auto result = composite->stats_json();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "scc.release") {
    const auto scc_handle = required_string(root, "scc");
    std::shared_ptr<CompositeSCCChartBase> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->sccs.find(scc_handle);
      if (found == session->sccs.end())
        throw std::invalid_argument(
            "unknown or already released persistent SCC chart");
      removed = std::move(found->second);
      session->scc_handles_by_key.erase(removed->key());
      session->sccs.erase(found);
    }
    return json::object{{"status", "ok"}, {"released", scc_handle}};
  }

  if (operation == "chart.release") {
    const auto chart_handle = required_string(root, "chart");
    std::shared_ptr<PreparedChartBase> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->charts.find(chart_handle);
      if (found == session->charts.end())
        throw std::invalid_argument("unknown or already released chart");
      removed = found->second;
      session->handles_by_key.erase(removed->key());
      session->charts.erase(found);
    }
    return json::object{{"status", "ok"}, {"released", chart_handle}};
  }

  if (operation == "session.stats") {
    std::vector<std::shared_ptr<PreparedChartBase>> charts;
    std::vector<std::shared_ptr<StoredLocalBase>> locals;
    std::vector<std::shared_ptr<StoredMatchBase>> matches;
    std::vector<std::shared_ptr<StoredEndpointResult>> endpoints;
    std::vector<std::shared_ptr<StoredTilePlan>> tile_plans;
    std::vector<std::shared_ptr<StoredTransportArmState>> transport_states;
    std::vector<std::shared_ptr<StoredLineResult>> line_results;
    std::vector<std::shared_ptr<CompositeSCCChartBase>> sccs;
    std::size_t pending_local_solves = 0;
    std::size_t pending_matches = 0;
    std::size_t pending_endpoint_limits = 0;
    std::size_t pending_tile_plans = 0;
    std::size_t pending_transport_states = 0;
    std::size_t pending_line_integrations = 0;
    std::uint64_t total_local_solves = 0;
    std::uint64_t total_scc_column_solves = 0;
    std::uint64_t total_local_matches = 0;
    std::uint64_t checkpoint_generation = 0;
    std::uint64_t checkpoint_restore_count = 0;
    std::string restored_from_checkpoint_identity;
    std::uint64_t total_endpoint_limits = 0;
    std::uint64_t total_endpoint_exports = 0;
    std::uint64_t total_tile_plans = 0;
    std::uint64_t total_transport_arm_marches = 0;
    std::uint64_t total_line_integrations = 0;
    std::uint64_t total_line_exports = 0;
    double total_local_run_parse_ms = 0.0, total_local_kernel_ms = 0.0;
    double total_local_match_ms = 0.0;
    double total_endpoint_limit_ms = 0.0;
    double total_endpoint_export_ms = 0.0;
    double total_tile_plan_ms = 0.0;
    double total_transport_arm_ms = 0.0;
    double total_line_integration_ms = 0.0;
    double total_line_export_ms = 0.0;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      for (const auto& [ignored, chart] : session->charts)
        charts.push_back(chart);
      for (const auto& [ignored, local] : session->locals)
        locals.push_back(local);
      for (const auto& [ignored, match] : session->matches)
        matches.push_back(match);
      for (const auto& [ignored, endpoint] : session->endpoints)
        endpoints.push_back(endpoint);
      for (const auto& [ignored, plan] : session->tile_plans)
        tile_plans.push_back(plan);
      for (const auto& [ignored, state] : session->transport_states)
        transport_states.push_back(state);
      for (const auto& [ignored, result] : session->line_results)
        line_results.push_back(result);
      for (const auto& [ignored, composite] : session->sccs)
        sccs.push_back(composite);
      pending_local_solves = session->pending_local_solves;
      pending_matches = session->pending_matches;
      pending_endpoint_limits = session->pending_endpoint_limits;
      pending_tile_plans = session->pending_tile_plans;
      pending_transport_states = session->pending_transport_states;
      pending_line_integrations = session->pending_line_integrations;
      total_local_solves = session->total_local_solves;
      total_scc_column_solves = session->total_scc_column_solves;
      total_local_matches = session->total_local_matches;
      total_endpoint_limits = session->total_endpoint_limits;
      total_endpoint_exports = session->total_endpoint_exports;
      total_tile_plans = session->total_tile_plans;
      total_transport_arm_marches = session->total_transport_arm_marches;
      total_line_integrations = session->total_line_integrations;
      total_line_exports = session->total_line_exports;
      total_local_run_parse_ms = session->total_local_run_parse_ms;
      total_local_kernel_ms = session->total_local_kernel_ms;
      total_local_match_ms = session->total_local_match_ms;
      checkpoint_generation = session->checkpoint_generation;
      checkpoint_restore_count = session->checkpoint_restore_count;
      restored_from_checkpoint_identity =
          session->restored_from_checkpoint_identity;
      total_endpoint_limit_ms = session->total_endpoint_limit_ms;
      total_endpoint_export_ms = session->total_endpoint_export_ms;
      total_tile_plan_ms = session->total_tile_plan_ms;
      total_transport_arm_ms = session->total_transport_arm_ms;
      total_line_integration_ms = session->total_line_integration_ms;
      total_line_export_ms = session->total_line_export_ms;
    }
    std::uint64_t runs = 0;
    double prepare_parse_ms = 0.0, run_parse_ms = 0.0, kernel_ms = 0.0;
    json::array chart_stats;
    for (const auto& chart : charts) {
      const auto stats = chart->stats();
      runs += stats.runs;
      prepare_parse_ms += stats.prepare_parse_ms;
      run_parse_ms += stats.run_parse_ms;
      kernel_ms += stats.kernel_ms;
      chart_stats.push_back(json::object{
          {"chart", chart->handle()}, {"key", chart->key()},
          {"dimension", chart->dimension()},
          {"frame_base", chart->frame_base()},
          {"frame_width", chart->frame_width()}, {"runs", stats.runs},
          {"d0_inverse_mode", chart->d0_inverse_mode()},
          {"local_solves", stats.local_runs},
          {"scc_components", chart->scc().component_count},
          {"scc_structural_edges", chart->scc().structural_edges.size()},
          {"scc_condensation_edges", chart->scc().condensation_edges.size()},
          {"scc_topological_order",
           encode_indices(chart->scc().topological_order)},
          {"scc_coupling_depth", chart->scc().coupling_depth},
          {"prepare_parse_ms", stats.prepare_parse_ms},
          {"run_parse_ms", stats.run_parse_ms},
          {"kernel_ms", stats.kernel_ms},
          {"local_run_parse_ms", stats.local_run_parse_ms},
          {"local_kernel_ms", stats.local_kernel_ms}});
    }
    std::uint64_t local_evaluations = 0;
    std::uint64_t local_residual_certifications = 0;
    std::uint64_t local_endpoint_limits = 0;
    std::uint64_t local_line_integrations = 0;
    std::uint64_t local_tail_certificate_requests = 0;
    std::uint64_t local_tail_certificate_certified = 0;
    std::uint64_t local_tail_certificate_inconclusive = 0;
    std::uint64_t local_tail_certificate_unsupported = 0;
    std::size_t local_coefficients = 0;
    double local_evaluate_ms = 0.0, local_residual_certify_ms = 0.0;
    double local_endpoint_limit_ms = 0.0;
    double local_line_integration_ms = 0.0;
    json::array local_stats;
    for (const auto& local : locals) {
      const auto stats = local->stats();
      local_evaluations += stats.evaluations;
      local_residual_certifications += stats.residual_certifications;
      local_endpoint_limits += stats.endpoint_limits;
      local_line_integrations += stats.line_integrations;
      local_tail_certificate_requests += stats.tail_certificate_requests;
      local_tail_certificate_certified += stats.tail_certificate_certified;
      local_tail_certificate_inconclusive +=
          stats.tail_certificate_inconclusive;
      local_tail_certificate_unsupported +=
          stats.tail_certificate_unsupported;
      local_coefficients += stats.coefficient_count;
      local_evaluate_ms += stats.evaluate_ms;
      local_residual_certify_ms += stats.residual_certify_ms;
      local_endpoint_limit_ms += stats.endpoint_limit_ms;
      local_line_integration_ms += stats.line_integration_ms;
      auto encoded = local->stats_json();
      local_stats.push_back(std::move(encoded));
    }
    json::array scc_stats;
    for (const auto& composite : sccs)
      scc_stats.push_back(composite->stats_json());
    json::array match_stats;
    for (const auto& match : matches)
      match_stats.push_back(match->summary());
    json::array endpoint_stats;
    for (const auto& endpoint : endpoints)
      endpoint_stats.push_back(endpoint->stats_json());
    json::array tile_plan_stats;
    for (const auto& plan : tile_plans)
      tile_plan_stats.push_back(plan->summary(false));
    json::array transport_state_stats;
    for (const auto& state : transport_states)
      transport_state_stats.push_back(state->stats_json());
    json::array line_result_stats;
    for (const auto& result : line_results)
      line_result_stats.push_back(result->stats_json());
    return json::object{{"status", "ok"}, {"session", session->handle},
                        {"domain", session->domain},
                        {"precision_bits", session->precision_bits},
                        {"charts", charts.size()}, {"runs", runs},
                        {"locals", locals.size()},
                        {"matches", matches.size()},
                        {"endpoints", endpoints.size()},
                        {"tile_plans", tile_plans.size()},
                        {"transport_states", transport_states.size()},
                        {"line_results", line_results.size()},
                        {"scc_charts", sccs.size()},
                        {"pending_local_solves", pending_local_solves},
                        {"pending_matches", pending_matches},
                        {"pending_endpoint_limits", pending_endpoint_limits},
                        {"pending_tile_plans", pending_tile_plans},
                        {"pending_transport_states",
                         pending_transport_states},
                        {"pending_line_integrations",
                         pending_line_integrations},
                        {"local_solves", total_local_solves},
                        {"local_matches", total_local_matches},
                        {"endpoint_limits", total_endpoint_limits},
                        {"endpoint_exports", total_endpoint_exports},
                        {"tile_plans_created", total_tile_plans},
                        {"transport_arm_marches",
                         total_transport_arm_marches},
                        {"line_integrations", total_line_integrations},
                        {"line_exports", total_line_exports},
                        {"local_match_capability",
                         session->domain == "rational"
                             ? kExactRegularLocalMatchCapability
                             : "unsupported"},
                        {"acb_local_match_capability",
                         session->domain == "acb"
                             ? kRefinedAcbLocalMatchCapability
                             : "unsupported"},
                        {"endpoint_limit_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedEndpointLimitCapability},
                        {"planned_endpoint_limit_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedPlannedEndpointLimitCapability},
                        {"tile_plan_capability", kRetainedTilePlanCapability},
                        {"single_arm_tile_plan_capability",
                         kRetainedSingleArmTilePlanCapability},
                        {"planned_match_hop_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedPlannedMatchHopCapability},
                        {"planned_match_materialization_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedPlannedMatchMaterializationCapability},
                        {"rational_row_application_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedRationalRowCapability},
                        {"line_integration_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedStoredLineCapability},
                        {"transport_arm_state_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedTransportArmStateCapability},
                        {"certified_tail_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRegularTailMajorantCapability},
                        {"certified_line_integration_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedCertifiedLineCapability},
                        {"scc_column_solves", total_scc_column_solves},
                        {"local_evaluations", local_evaluations},
                        {"local_residual_certifications",
                         local_residual_certifications},
                        {"local_endpoint_limits", local_endpoint_limits},
                        {"local_line_integrations",
                         local_line_integrations},
                        {"local_tail_certificate_requests",
                         local_tail_certificate_requests},
                        {"local_tail_certificate_certified",
                         local_tail_certificate_certified},
                        {"local_tail_certificate_inconclusive",
                         local_tail_certificate_inconclusive},
                        {"local_tail_certificate_unsupported",
                         local_tail_certificate_unsupported},
                        {"local_coefficient_count", local_coefficients},
                        {"static_tensor_copies", 0},
                        {"prepare_parse_ms", prepare_parse_ms},
                        {"run_parse_ms", run_parse_ms},
                        {"kernel_ms", kernel_ms},
                        {"local_run_parse_ms", total_local_run_parse_ms},
                        {"local_kernel_ms", total_local_kernel_ms},
                        {"local_evaluate_ms", local_evaluate_ms},
                        {"local_residual_certify_ms",
                         local_residual_certify_ms},
                        {"local_endpoint_limit_ms", local_endpoint_limit_ms},
                        {"local_line_integration_ms",
                         local_line_integration_ms},
                        {"local_match_ms", total_local_match_ms},
                        {"checkpoint_generation", checkpoint_generation},
                        {"checkpoint_restore_count",
                         checkpoint_restore_count},
                        {"restored_from_checkpoint_identity",
                         restored_from_checkpoint_identity.empty()
                             ? json::value(nullptr)
                             : json::value(
                                   restored_from_checkpoint_identity)},
                        {"endpoint_limit_ms", total_endpoint_limit_ms},
                        {"endpoint_export_ms", total_endpoint_export_ms},
                        {"tile_plan_ms", total_tile_plan_ms},
                        {"transport_arm_ms", total_transport_arm_ms},
                        {"line_integration_ms", total_line_integration_ms},
                        {"line_export_ms", total_line_export_ms},
                        {"chart_stats", std::move(chart_stats)},
                        {"local_stats", std::move(local_stats)},
                        {"match_stats", std::move(match_stats)},
                        {"endpoint_stats", std::move(endpoint_stats)},
                        {"tile_plan_stats", std::move(tile_plan_stats)},
                        {"transport_state_stats",
                         std::move(transport_state_stats)},
                        {"line_result_stats", std::move(line_result_stats)},
                        {"scc_stats", std::move(scc_stats)}};
  }

  throw std::invalid_argument("unknown persistent solver operation: " + operation);
}

json::object run_one(const json::object& root) {
  if (as_i64(root.at("schema"), "schema") != 1)
    throw std::invalid_argument("unsupported recurrence schema");
  const auto domain = std::string(root.at("domain").as_string());
  const int digits = root.if_contains("output_digits")
                         ? static_cast<int>(as_i64(root.at("output_digits"), "output digits"))
                         : 50;
  if (domain == "rational") return run_typed<Rational>(root, digits);
  if (domain == "acb") {
    const auto precision_bits = static_cast<slong>(
        as_i64(root.at("precision_bits"), "precision bits"));
    AcbPrecisionLease lease(precision_bits);
    ComplexBall::set_precision(precision_bits);
    return run_typed<ComplexBall>(root, digits);
  }
  if (domain == "symbolic") {
    std::lock_guard<std::recursive_mutex> lock(symbolic_run_mutex());
    std::vector<std::string> variables;
    for (const auto& value : as_array(root.at("symbols"), "symbolic variables")) {
      if (!value.is_string())
        throw std::invalid_argument("symbolic variable names must be strings");
      variables.emplace_back(value.as_string());
    }
    SymbolicRational::configure(variables);
    return run_typed<SymbolicRational>(root, digits);
  }
  throw std::invalid_argument("unsupported scalar domain: " + domain);
}

json::object run_one_safe(const json::object& root) {
  try {
    return run_one(root);
  } catch (const RecurrenceError& error) {
    return json::object{{"status", "error"}, {"id", error.id},
      {"detail", error.what()}, {"frame_base", error.frame_base},
      {"shift", error.shift}};
  } catch (const std::exception& error) {
    return json::object{{"status", "error"}, {"id", "CPP"},
                        {"detail", error.what()}};
  }
}

}  // namespace

std::string run_recurrence_json(std::string_view input) {
  try {
    const auto parsed = json::parse(input);
    const auto& root = as_object(parsed, "root");
    if (const auto* schema = root.if_contains("schema");
        schema != nullptr && as_i64(*schema, "schema") == 2)
      return json::serialize(run_session_command(root));
    if (const auto* raw_batch = root.if_contains("batch"); raw_batch != nullptr) {
      const auto& requests = as_array(*raw_batch, "batch requests");
      const auto requested_threads = root.if_contains("threads")
          ? as_u32(root.at("threads"), "batch threads") : 1;
      auto thread_count = std::max<std::uint32_t>(1,
          std::min<std::uint32_t>(requested_threads,
                                 static_cast<std::uint32_t>(requests.size())));
      // Acb precision is process-global in the light RAII wrapper. All batch
      // requests must therefore agree before workers start.
      std::optional<std::int64_t> acb_precision;
      std::optional<std::vector<std::string>> symbolic_variables;
      for (const auto& raw_request : requests) {
        const auto& request = as_object(raw_request, "batch request");
        const auto domain = std::string(request.at("domain").as_string());
        if (domain == "acb") {
          const auto bits = as_i64(request.at("precision_bits"), "precision bits");
          if (acb_precision.has_value() && *acb_precision != bits)
            throw std::invalid_argument("all Acb batch requests must use one precision");
          acb_precision = bits;
        } else if (domain == "symbolic") {
          std::vector<std::string> variables;
          for (const auto& value : as_array(
                   request.at("symbols"), "symbolic variables"))
            variables.emplace_back(value.as_string());
          if (symbolic_variables.has_value() &&
              *symbolic_variables != variables)
            throw std::invalid_argument(
                "all symbolic batch requests must use one coefficient field");
          symbolic_variables = std::move(variables);
        }
      }
      if (acb_precision.has_value())
        ComplexBall::set_precision(static_cast<slong>(*acb_precision));
      if (symbolic_variables.has_value()) {
        SymbolicRational::configure(*symbolic_variables);
        // FLINT multivariate contexts are shared read-only, but keep this
        // first symbolic milestone single-threaded until its allocator-level
        // parallel behavior is separately certified.
        thread_count = 1;
      }
      std::vector<json::object> results(requests.size());
      std::atomic<std::size_t> next{0};
      auto worker = [&]() {
        while (true) {
          const auto index = next.fetch_add(1);
          if (index >= requests.size()) return;
          results[index] = run_one_safe(
              as_object(requests[index], "batch request"));
        }
      };
      std::vector<std::jthread> workers;
      workers.reserve(thread_count);
      for (std::uint32_t i = 0; i < thread_count; ++i)
        workers.emplace_back(worker);
      for (auto& thread : workers) thread.join();
      json::array encoded;
      encoded.reserve(results.size());
      for (auto& result : results) encoded.push_back(std::move(result));
      return json::serialize(json::object{{"status", "ok"},
                                          {"results", std::move(encoded)}});
    }
    return json::serialize(run_one_safe(root));
  } catch (const NativeIntegrationError& error) {
    json::object result{{"status", "error"}, {"id", error.id},
                        {"detail", error.what()}};
    if (!error.absolute_power.empty()) {
      result["absolute_power"] = error.absolute_power;
      result["log_power"] = error.log_power;
      result["epsilon_power"] = error.epsilon_power;
      result["component"] = error.component;
    }
    return json::serialize(result);
  } catch (const RecurrenceError& error) {
    return json::serialize(json::object{{"status", "error"}, {"id", error.id},
      {"detail", error.what()}, {"frame_base", error.frame_base}, {"shift", error.shift}});
  } catch (const std::exception& error) {
    return json::serialize(json::object{{"status", "error"}, {"id", "CPP"},
                                        {"detail", error.what()}});
  }
}

void reset_solver_sessions() {
  auto& registry = session_registry();
  std::unordered_map<std::string, std::shared_ptr<SolverSession>> sessions;
  {
    std::lock_guard<std::mutex> lock(registry.mutex);
    sessions.swap(registry.sessions);
  }
  // Destroy retained FLINT/Arb objects while the library and coefficient
  // contexts are still alive.  In-flight calls keep shared ownership until
  // they leave the backend.
  sessions.clear();
}

std::string backend_info_json() {
  return json::serialize(json::object{{"schema", 1},
                                      {"schemas", json::array{1, 2}},
                                      {"persistent_sessions", true},
                                      {"persistent_local_solutions", true},
                                      {"persistent_local_residual_certification",
                                       true},
                                      {"persistent_scc_prepare", true},
                                      {"persistent_scc_execute", false},
                                      {"persistent_scc_scalar_block_dag_column",
                                       true},
                                      {"persistent_scc_regular_block_dag_column",
                                       true},
                                      {"persistent_exact_regular_local_match",
                                       true},
                                      {"persistent_endpoint_limits", true},
                                      {"persistent_endpoint_limit_capability",
                                       kRetainedEndpointLimitCapability},
                                      {"persistent_plan_bound_endpoint_limits",
                                       true},
                                      {"persistent_plan_bound_endpoint_limit_capability",
                                       kRetainedPlannedEndpointLimitCapability},
                                      {"persistent_symbolic_endpoint_limits",
                                       false},
                                      {"persistent_exact_tile_plans", true},
                                      {"persistent_exact_tile_plan_capability",
                                       kRetainedTilePlanCapability},
                                      {"persistent_exact_single_arm_tile_plans",
                                       true},
                                      {"persistent_exact_single_arm_tile_plan_capability",
                                       kRetainedSingleArmTilePlanCapability},
                                      {"persistent_plan_driven_match_hop",
                                       true},
                                      {"persistent_plan_driven_match_hop_capability",
                                       kRetainedPlannedMatchHopCapability},
                                      {"persistent_plan_match_local_materialization",
                                       true},
                                      {"persistent_plan_match_local_materialization_capability",
                                       kRetainedPlannedMatchMaterializationCapability},
                                      {"persistent_stored_line_integration",
                                       true},
                                      {"persistent_stored_line_integration_capability",
                                       kRetainedStoredLineCapability},
                                      {"persistent_parallel_arm_march", true},
                                      {"persistent_parallel_arm_march_capability",
                                       kRetainedParallelArmCapability},
                                      {"persistent_transport_arm_state", true},
                                      {"persistent_transport_arm_state_capability",
                                       kRetainedTransportArmStateCapability},
                                      {"persistent_certified_tail_majorant",
                                       true},
                                      {"persistent_certified_tail_majorant_capability",
                                       kRegularTailMajorantCapability},
                                      {"persistent_certified_line_integration_capability",
                                       kRetainedCertifiedLineCapability},
                                      {"persistent_symbolic_line_integration",
                                       false},
                                      {"persistent_acb_local_match", true},
                                      {"persistent_acb_local_match_capability",
                                       kRefinedAcbLocalMatchCapability},
                                      {"persistent_checkpoint", true},
                                      {"persistent_checkpoint_schema",
                                       kCheckpointPayloadSchema},
                                      {"persistent_checkpoint_handle_scope",
                                       "complete-retained-native-ownership-closure"},
                                      {"persistent_scc_regular_singular_scalar_block_dag_column",
                                       true},
                                      {"persistent_scc_regular_singular_jordan_block_dag_column",
                                       true},
                                      {"persistent_scc_pseudo_compensation",
                                       false},
                                      {"backend", "DiffExp2 C++"},
                                      {"flint", flint_version},
                                      {"librarylink", true}});
}

}  // namespace diffexp2
