#include "diffexp2/json_codec.hpp"

#include "diffexp2/checkpoint.hpp"
#include "diffexp2/local_algebra.hpp"
#include "diffexp2/local_solution.hpp"
#include "diffexp2/matching.hpp"
#include "diffexp2/recurrence.hpp"

#include <boost/json.hpp>

#include <chrono>
#include <atomic>
#include <algorithm>
#include <condition_variable>
#include <cstdint>
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
  double evaluate_ms = 0.0;
  double residual_certify_ms = 0.0;
  double endpoint_limit_ms = 0.0;
  double create_parse_ms = 0.0;
  double create_kernel_ms = 0.0;
  std::size_t coefficient_count = 0;
};

struct NativeLocalDiagnostics {
  std::int32_t top_valid = kCompleteInfinity;
  double parse_ms = 0.0;
  double kernel_ms = 0.0;
};

template <typename Scalar>
struct NativeLocalRun {
  LocalSolution<Scalar> solution;
  std::vector<PseudoHit<Scalar>> pseudo_hits;
  NativeLocalDiagnostics diagnostics;
};

json::value canonical_json_value(const json::value& value);

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
                  double create_parse_ms, double create_kernel_ms,
                  std::optional<SCCColumnProvenance> column_provenance =
                      std::nullopt)
      : handle_(std::move(handle)), source_chart_(std::move(source_chart)),
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
  virtual json::object endpoint_metadata() const = 0;
  virtual const std::string& checkpoint_identity() const = 0;
  virtual const char* scalar_domain() const = 0;
  virtual json::object summary() const = 0;
  virtual json::object stats_json() const = 0;
  virtual StoredLocalStats stats() const = 0;

  const std::string& handle() const { return handle_; }
  const std::string& source_chart() const { return source_chart_; }
  const std::optional<SCCColumnProvenance>& column_provenance() const {
    return column_provenance_;
  }

 protected:
  std::string handle_;
  std::string source_chart_;
  double create_parse_ms_ = 0.0;
  double create_kernel_ms_ = 0.0;
  std::optional<SCCColumnProvenance> column_provenance_;
};

template <typename Scalar>
class StoredLocal final : public StoredLocalBase {
 public:
  StoredLocal(std::string handle, std::string source_chart,
              LocalSolution<Scalar>&& solution, slong precision_bits,
              std::vector<PseudoHit<Scalar>>&& pseudo_hits,
              NativeLocalDiagnostics diagnostics,
              std::optional<SCCColumnProvenance> column_provenance =
                  std::nullopt)
      : StoredLocalBase(std::move(handle), std::move(source_chart),
                        diagnostics.parse_ms, diagnostics.kernel_ms,
                        std::move(column_provenance)),
        solution_(std::move(solution)), precision_bits_(precision_bits),
        pseudo_hits_(std::move(pseudo_hits)), top_valid_(diagnostics.top_valid) {
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
      const auto options = parse_local_evaluation_options(request, true);
      const auto started = std::chrono::steady_clock::now();
      auto result = evaluate_local_solution(solution_, point, options);
      const auto elapsed = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      evaluations_.fetch_add(1);
      {
        std::lock_guard<std::mutex> lock(stats_mutex_);
        evaluate_ms_ += elapsed;
      }
      return json::object{
          {"point_exact", point.exact_coordinate},
          {"imaginary_sign", result.imaginary_sign.has_value()
               ? json::value(*result.imaginary_sign) : json::value(nullptr)},
          {"arithmetic_enclosed", result.arithmetic_enclosed},
          {"elapsed_ms", elapsed},
          {"value", encode_epsilon_vector(result.value, output_digits)},
          {"theta", encode_epsilon_vector(result.theta_value, output_digits)}};
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

  json::object endpoint_metadata() const override { return metadata_json(); }

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
        {"dimension", solution_.dimension},
        {"epsilon_min", solution_.epsilon.min_power},
        {"epsilon_max", solution_.epsilon.complete_max},
        {"taylor_complete_max", solution_.taylor_complete_max},
        {"sectors", solution_.sectors.size()},
        {"coefficient_count", coefficient_count()},
        {"pseudo_hit_count", pseudo_hits_.size()},
        {"top_valid", encode_validity(top_valid_)},
        {"checkpoint_identity", solution_.checkpoint_identity},
        {"metadata", metadata_json()},
        {"create_parse_ms", create_parse_ms_},
        {"create_kernel_ms", create_kernel_ms_}};
    if (column_provenance_.has_value())
      result["column_provenance"] = column_provenance_->encode();
    return result;
  }

  json::object stats_json() const override {
    auto out = summary();
    const auto current = stats();
    out["evaluations"] = current.evaluations;
    out["residual_certifications"] = current.residual_certifications;
    out["endpoint_limits"] = current.endpoint_limits;
    out["evaluate_ms"] = current.evaluate_ms;
    out["residual_certify_ms"] = current.residual_certify_ms;
    out["endpoint_limit_ms"] = current.endpoint_limit_ms;
    return out;
  }

  StoredLocalStats stats() const override {
    std::lock_guard<std::mutex> lock(stats_mutex_);
    return {evaluations_.load(), residual_certifications_.load(),
            endpoint_limits_.load(), evaluate_ms_, residual_certify_ms_,
            endpoint_limit_ms_, create_parse_ms_, create_kernel_ms_,
            coefficient_count()};
  }

  const LocalSolution<Scalar>& solution() const { return solution_; }
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

  LocalSolution<Scalar> solution_;
  slong precision_bits_ = 256;
  std::vector<PseudoHit<Scalar>> pseudo_hits_;
  std::int32_t top_valid_ = kCompleteInfinity;
  std::atomic<std::uint64_t> evaluations_{0};
  std::atomic<std::uint64_t> residual_certifications_{0};
  std::atomic<std::uint64_t> endpoint_limits_{0};
  mutable std::mutex stats_mutex_;
  double evaluate_ms_ = 0.0;
  double residual_certify_ms_ = 0.0;
  double endpoint_limit_ms_ = 0.0;
};

constexpr const char* kRetainedEndpointLimitCapability =
    "retained-native-endpoint-sector-limit-v1";

struct ParsedEndpointLimitPolicy {
  EndpointLimitOptions options;
  std::string cancellation_mode;
  std::optional<std::int32_t> requested_rim;
};

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

  const auto& cancellation = as_object(
      request.at("cancellation"), "endpoint cancellation policy");
  parsed.cancellation_mode = required_string(cancellation, "mode");
  if (parsed.cancellation_mode == "exact-coefficient-field") {
    parsed.options.allow_certified_numeric_cancellation = false;
  } else if (parsed.cancellation_mode == "exact-or-acb-singleton") {
    parsed.options.allow_certified_numeric_cancellation = true;
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
      std::string source_chart, std::string source_checkpoint,
      std::string source_domain, std::int32_t approach_direction,
      std::optional<std::int32_t> requested_rim,
      std::string cancellation_mode, json::object analytic_metadata,
      EndpointLimitResult&& result, double elapsed_ms)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        source_local_(std::move(source_local)),
        source_chart_(std::move(source_chart)),
        source_checkpoint_(std::move(source_checkpoint)),
        source_domain_(std::move(source_domain)),
        approach_direction_(approach_direction),
        requested_rim_(requested_rim),
        cancellation_mode_(std::move(cancellation_mode)),
        analytic_metadata_(std::move(analytic_metadata)),
        result_(std::move(result)), elapsed_ms_(elapsed_ms) {
    // Validate the retained public frame once, before publishing its handle.
    (void)endpoint_value_window(result_);
  }

  const std::string& handle() const { return handle_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  double elapsed_ms() const { return elapsed_ms_; }

  json::object summary() const {
    const auto window = endpoint_value_window(result_);
    const auto cancellation_scope = source_domain_ == "rational"
        ? "exact-rational"
        : cancellation_mode_ == "exact-or-acb-singleton"
            ? "acb-exact-singleton-zero"
            : "exact-coefficient-field-only";
    return json::object{
        {"endpoint", handle_},
        {"capability", kRetainedEndpointLimitCapability},
        {"native_retained", true},
        {"retained_state", "specialized-acb-epsilon-vector"},
        {"json_coefficients", 0},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"source", json::object{
             {"local", source_local_}, {"chart", source_chart_},
             {"checkpoint_identity", source_checkpoint_},
             {"coefficient_domain", source_domain_}}},
        {"dimension", result_.values.size()},
        {"epsilon_min", window.min_power},
        {"epsilon_max", window.complete_max},
        {"coefficient_field", "acb-specialized"},
        {"arithmetic_enclosed", true},
        {"approach_direction", approach_direction_},
        {"requested_rim", requested_rim_.has_value()
             ? json::value(*requested_rim_) : json::value(nullptr)},
        {"effective_rim", result_.imaginary_sign},
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
  }

  json::object stats_json() const {
    auto out = summary();
    std::lock_guard<std::mutex> lock(stats_mutex_);
    out["exports"] = exports_;
    out["export_ms"] = export_ms_;
    return out;
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
  std::string source_checkpoint_;
  std::string source_domain_;
  std::int32_t approach_direction_ = 1;
  std::optional<std::int32_t> requested_rim_;
  std::string cancellation_mode_;
  json::object analytic_metadata_;
  EndpointLimitResult result_;
  double elapsed_ms_ = 0.0;
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
  auto analytic_metadata = local->endpoint_metadata();
  json::object provenance{
      {"schema", "diffexp2-retained-native-endpoint-sector-limit-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"source", json::object{
           {"local", local->handle()}, {"chart", local->source_chart()},
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
      local->handle(), local->source_chart(), expected_source_checkpoint,
      local->scalar_domain(), policy.options.approach_direction,
      policy.requested_rim, policy.cancellation_mode,
      std::move(analytic_metadata), std::move(result), elapsed);
}

constexpr const char* kExactRegularLocalMatchCapability =
    "exact-rational-regular-local-match-v1";

class StoredMatchBase {
 public:
  explicit StoredMatchBase(std::string handle) : handle_(std::move(handle)) {}
  virtual ~StoredMatchBase() = default;

  virtual json::object summary() const = 0;
  const std::string& handle() const { return handle_; }

 protected:
  std::string handle_;
};

class StoredExactRegularMatch final : public StoredMatchBase {
 public:
  StoredExactRegularMatch(
      std::string handle, std::string checkpoint_identity,
      std::string provenance_identity, std::vector<std::string> basis_handles,
      std::vector<std::string> basis_checkpoints, std::string incoming_handle,
      std::string incoming_checkpoint, std::string basis_chart,
      std::string incoming_chart, std::string basis_point,
      std::string incoming_point, std::string physical_point,
      EpsilonWindow requested_window, std::int32_t required_complete_max,
      std::uint32_t dimension,
      ExactLaurentMatrix<Rational>&& transformation,
      FiniteLaurentVector<Rational>&& weights,
      EpsilonLatticeSaturationDiagnostics<Rational>&& diagnostics,
      EpsilonWindow residual_window, double elapsed_ms)
      : StoredMatchBase(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        basis_handles_(std::move(basis_handles)),
        basis_checkpoints_(std::move(basis_checkpoints)),
        incoming_handle_(std::move(incoming_handle)),
        incoming_checkpoint_(std::move(incoming_checkpoint)),
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
        elapsed_ms_(elapsed_ms) {}

  json::object summary() const override {
    json::array basis;
    basis.reserve(basis_handles_.size());
    for (std::size_t column = 0; column < basis_handles_.size(); ++column)
      basis.push_back(json::object{{"column", column},
                                   {"local", basis_handles_[column]},
                                   {"checkpoint_identity",
                                    basis_checkpoints_[column]}});

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
        {"incoming", incoming_handle_},
        {"incoming_checkpoint_identity", incoming_checkpoint_},
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

 private:
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::vector<std::string> basis_handles_;
  std::vector<std::string> basis_checkpoints_;
  std::string incoming_handle_;
  std::string incoming_checkpoint_;
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
  if (requested_window.min_power < solution.epsilon.min_power ||
      requested_window.complete_max > solution.epsilon.complete_max)
    throw std::invalid_argument(
        label + " does not cover the requested complete epsilon window");
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
    const RealEvaluationPoint& point, EpsilonWindow window) {
  const Rational t(point.exact_coordinate);
  std::vector<Rational> t_powers(solution.taylor_width(), Rational(1));
  for (std::size_t n = 1; n < t_powers.size(); ++n)
    t_powers[n] = t_powers[n - 1] * t;

  const auto& sector = solution.sectors.front();
  FiniteLaurentVector<Rational> value;
  value.reserve(solution.dimension);
  for (std::uint32_t component = 0; component < solution.dimension;
       ++component) {
    std::vector<Rational> coefficients;
    coefficients.reserve(window.width());
    for (std::int64_t power = window.min_power;
         power <= window.complete_max; ++power) {
      const auto epsilon_index = static_cast<std::size_t>(
          power - solution.epsilon.min_power);
      Rational coefficient(0);
      for (std::size_t n = 0; n < solution.taylor_width(); ++n)
        coefficient += sector.coefficients[local_detail::sector_index(
                           solution, epsilon_index, n, component)] *
                       t_powers[n];
      coefficients.push_back(std::move(coefficient));
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
  for (const auto& column : basis) {
    auto value = evaluate_exact_regular_local(
        column->solution(), basis_point, window);
    for (std::uint32_t component = 0; component < dimension; ++component)
      evaluated_basis[component].push_back(std::move(value[component]));
  }
  auto incoming_value = evaluate_exact_regular_local(
      incoming->solution(), incoming_point, window);

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

  json::array provenance_basis;
  for (std::size_t column = 0; column < basis.size(); ++column) {
    json::object entry{{"column", column},
                       {"local", basis_handles[column]},
                       {"chart", basis_chart},
                       {"checkpoint_identity", basis_checkpoints[column]}};
    if (basis[column]->column_provenance().has_value())
      entry["column_provenance"] =
          basis[column]->column_provenance()->encode();
    provenance_basis.push_back(std::move(entry));
  }
  json::object provenance{
      {"schema", "diffexp2-native-exact-regular-local-match-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"basis", std::move(provenance_basis)},
      {"incoming", json::object{
           {"local", incoming_handle}, {"chart", incoming_chart},
           {"checkpoint_identity", expected_incoming_checkpoint}}},
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
      match_handle, checkpoint_identity, provenance_identity, basis_handles,
      std::move(basis_checkpoints), incoming_handle,
      expected_incoming_checkpoint, basis_chart, incoming_chart,
      basis_point.exact_coordinate, incoming_point.exact_coordinate,
      basis_physical_point.str(), window, required_complete_max, dimension,
      std::move(saturated.transformation), std::move(weights),
      std::move(saturated.diagnostics),
      EpsilonWindow{residual_min, residual_max}, elapsed_ms);
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
        symbols_(std::move(symbols)) {}

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
    return solve_native_impl(run, metadata_object, std::nullopt);
  }

  NativeLocalRun<Scalar> solve_native_with_source(
      const json::object& run, const json::object& metadata_object,
      SourceData<Scalar>&& source) {
    if (!run.at("source").is_null())
      throw std::invalid_argument(
          "native SCC source injection rejects caller-supplied source data");
    return solve_native_impl(
        run, metadata_object, std::move(source));
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
  std::optional<std::pair<Rational, Rational>>
  exact_scalar_affine_indicial_root() const {
    if constexpr (!std::is_same_v<Scalar, Rational>) {
      return std::nullopt;
    } else {
      if (prepared_.dimension != 1 || prepared_.nhat_lags.empty() ||
          !prepared_.d0_inverse_scalar.has_value() ||
          !prepared_.nhat_lags.front().rational.empty())
        return std::nullopt;
      std::map<std::int32_t, Rational> coefficients;
      for (const auto& matrix : prepared_.nhat_lags.front().polynomial) {
        for (const auto& entry : matrix.entries) {
          if (entry.row != 0 || entry.col != 0) return std::nullopt;
          auto found = coefficients.try_emplace(
              matrix.shift, Rational(0)).first;
          found->second += entry.value;
        }
      }
      for (auto iterator = coefficients.begin();
           iterator != coefficients.end();) {
        iterator->second *= *prepared_.d0_inverse_scalar;
        if (iterator->second.is_zero())
          iterator = coefficients.erase(iterator);
        else
          ++iterator;
      }
      if (std::any_of(coefficients.begin(), coefficients.end(),
                      [](const auto& item) {
                        return item.first != 0 && item.first != 1;
                      }))
        return std::nullopt;
      const auto coefficient = [&](std::int32_t shift) {
        const auto found = coefficients.find(shift);
        return found == coefficients.end() ? Rational(0) : found->second;
      };
      return std::make_pair(coefficient(0), coefficient(1));
    }
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
      std::optional<SourceData<Scalar>> native_source) {
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
    auto pseudo_hits = std::move(recurrence.hits);
    const auto kernel_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - kernel_started).count();
    return {std::move(solution), std::move(pseudo_hits),
            {recurrence.top_valid, parse_ms, kernel_ms}};
  }

  std::shared_ptr<StoredLocalBase> solve_local(
      const std::string& local_handle, const json::object& run,
      const json::object& metadata_object) override {
    auto native = solve_native(run, metadata_object);
    auto local = make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
        local_handle, handle_, std::move(native.solution), precision_bits_,
        std::move(native.pseudo_hits), native.diagnostics);
    record_native_local_success(native.diagnostics);
    return local;
  }
  PreparedRecurrenceOperator<Scalar> prepared_;
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
  std::shared_ptr<PreparedChart<Scalar>> chart;
};

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
    if constexpr (!std::is_same_v<Scalar, Rational>) {
      throw std::invalid_argument(
          "native SCC column execution currently supports only exact rational coefficients");
    } else {
      const auto started = std::chrono::steady_clock::now();
      const bool regular_execution = regular_block_column_ready();
      const bool regular_singular_execution =
          regular_singular_scalar_column_ready();
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

      const auto& seed_run = checked_column_run(
          seed_request, seed_block, true, nullptr,
          regular_singular_execution);
      const auto seed_local_component = seed_component_from_run(
          seed_run, seed_block, regular_singular_execution);
      const auto basis_index =
          blocks_[seed_block].vertices[seed_local_component];
      auto seed_native = blocks_[seed_block].chart->solve_native(
          seed_run, as_object(seed_request.at("metadata"),
                              "native SCC seed metadata"));
      validate_block_result(seed_native, seed_block, true);
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
        const auto& target_run = checked_column_run(
            target_request, target_block, false, &source,
            regular_singular_execution);
        require_source_tag_matches_run(source, target_run);
        auto source_data = local_solution_source_data(
            source, as_u32(target_run.at("nmax"), "target nmax"),
            as_u32(target_run.at("p"), "target log maximum"),
            blocks_[target_block].chart->frame_base(),
            blocks_[target_block].chart->frame_width());
        auto target_native = blocks_[target_block].chart->solve_native_with_source(
            target_run,
            as_object(target_request.at("metadata"),
                      "native SCC target metadata"),
            std::move(source_data));
        validate_block_result(target_native, target_block, false);
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
               ? "diffexp2-native-scc-regular-singular-scalar-column-v1"
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
      SCCColumnProvenance column_provenance{
          handle_, exact_identity_, seed_block,
          basis_index,
          json::serialize(canonical_json_value(column_identity_record))};
      auto local = make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
          local_handle, handle_, std::move(parent),
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
    if (regular_singular_scalar_column_ready())
      return "exact-rational-regular-singular-scalar-block-dag-column-v1";
    if (!regular_block_column_ready())
      return "unsupported-native-scc-column";
    return scalar_block_shape()
        ? "exact-rational-regular-scalar-block-dag-column-v1"
        : "exact-rational-regular-block-dag-column-v2";
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
          {"principal_identity", block.principal_identity}};
      if (const auto root =
              block.chart->exact_scalar_affine_indicial_root();
          root.has_value())
        block_record["affine_indicial_root"] = json::object{
            {"a", root->first.str()}, {"b", root->second.str()}};
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
        regular_singular_scalar_column_ready();
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
             ? "exact-rational-regular-singular-scalar-block-dag-column-v1"
             : (regular_ready
                   ? (scalar_shape
                         ? "exact-rational-regular-scalar-block-dag-column-v1"
                         : "exact-rational-regular-block-dag-column-v2")
                   : "unsupported")},
        {"general_scc_execution", false},
        {"scalar_block_dag_column_execution", scalar_ready},
        {"regular_singular_scalar_block_dag_column_execution",
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
             {"no_pseudo", "collision-bound-producer-certificate"},
             {"resonance_schedule",
              "retained-affine-root-verified-exact-captured-run"}}},
        {"execution_must_revalidate_producer_capabilities", true},
        {"block_charts", std::move(block_handles)}};
    if (!scalar_shape)
      result["regular_block_dag_column_execution"] = regular_ready;
    return result;
  }

 private:
  bool regular_block_column_ready() const {
    if constexpr (!std::is_same_v<Scalar, Rational>) {
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

  bool regular_singular_scalar_column_ready() const {
    if constexpr (!std::is_same_v<Scalar, Rational>) {
      return false;
    } else {
      if (blocks_.size() < 2 || work_.work_min > 0 ||
          work_.requested_max < 0 || work_.work_complete_max < 0 ||
          !scalar_block_shape() ||
          std::none_of(blocks_.begin(), blocks_.end(), [](const auto& block) {
            return !block.regular;
          }) ||
          std::any_of(blocks_.begin(), blocks_.end(), [](const auto& block) {
            return block.chart->dimension() != 1 ||
                   !block.chart->has_identity_assembly() ||
                   !block.chart->has_regular_singleton_partition() ||
                   block.chart->jordan_block_count() != 1 ||
                   !block.chart->exact_scalar_affine_indicial_root().has_value();
          }))
        return false;
      return std::all_of(
          couplings_.begin(), couplings_.end(), [&](const auto& coupling) {
            return sector_preserving_coupling_ready(coupling);
          });
    }
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

  const json::object& checked_column_run(
      const json::object& entry, std::uint32_t block_index, bool seed,
      const LocalSolution<Scalar>* source,
      bool regular_singular_execution) const {
    if (block_index >= blocks_.size())
      throw std::invalid_argument("native SCC run block is out of range");
    const auto block_dimension = blocks_[block_index].chart->dimension();
    const auto frame_width = blocks_[block_index].chart->frame_width();
    const auto retained_root = regular_singular_execution
        ? blocks_[block_index].chart->exact_scalar_affine_indicial_root()
        : std::nullopt;
    if (regular_singular_execution && !retained_root.has_value())
      throw std::invalid_argument(
          "native regular-singular SCC chart has no retained exact affine scalar indicial root");
    const auto& run = as_object(entry.at("run"), "native SCC recurrence run");
    validate_metadata_geometry(as_object(
        entry.at("metadata"), "native SCC local metadata"));
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
    const auto a_target = parse_scalar<Rational>(run.at("a_target"));
    const auto b_target = parse_scalar<Rational>(run.at("b_target"));
    if (as_i32(run.at("a_shift_min"),
               "native SCC a-shift minimum") != 0)
      throw std::invalid_argument(
          "native SCC column requires a zero exact a-shift origin");
    if (!regular_singular_execution &&
        (log_max != 0 || !a_target.is_zero() || !b_target.is_zero()))
      throw std::invalid_argument(
          "native SCC regular block runs require p=0 and a=b=0");
    if (regular_singular_execution &&
        (block_dimension != 1 ||
         blocks_[block_index].chart->jordan_block_count() != 1 ||
         (seed && log_max != 0)))
      throw std::invalid_argument(
          "native regular-singular SCC execution requires scalar blocks and a canonical log-zero homogeneous seed");
    const auto& a_shifts = as_array(
        run.at("a_shifts"), "native SCC exact a-shift schedule");
    if (a_shifts.size() != static_cast<std::size_t>(nmax) + 1)
      throw std::invalid_argument(
          "native SCC regular block run has an incomplete a-shift schedule");
    for (std::size_t n = 0; n < a_shifts.size(); ++n)
      if (!(parse_scalar<Rational>(a_shifts[n]) ==
            a_target + Rational(std::to_string(n))))
        throw std::invalid_argument(
            "native SCC a-shift schedule must equal a_target plus the exact Taylor index");
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
            return !parse_scalar<Rational>(value).is_zero();
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
    for (std::size_t n = 0; n < schedule.size(); ++n) {
      const auto& raw_row = schedule[n];
      const auto& row = as_array(raw_row, "native SCC schedule row");
      if (!blocks_[block_index].chart->has_regular_singleton_partition() ||
          row.size() != block_dimension)
        throw std::invalid_argument(
            "native SCC regular block execution requires one step per retained Jordan singleton");
      for (const auto& raw_step : row) {
        const auto& step = as_object(raw_step, "native SCC schedule step");
        const auto kind = required_string(step, "case");
        const auto da = parse_scalar<Rational>(step.at("da"));
        const auto db = parse_scalar<Rational>(step.at("db"));
        if (regular_singular_execution) {
          const auto expected_da =
              a_target + Rational(std::to_string(n)) - retained_root->first;
          const auto expected_db = b_target - retained_root->second;
          const auto expected_kind = !da.is_zero() ? "T"
              : !db.is_zero() ? "P" : "R";
          if (!(da == expected_da) || !(db == expected_db))
            throw std::invalid_argument(
                "native regular-singular SCC schedule offsets differ from the retained exact affine indicial root");
          if (kind != expected_kind)
            throw std::invalid_argument(
                "native regular-singular SCC schedule case contradicts its exact da/db offsets");
          if (kind == "P")
            throw std::invalid_argument(
                "native regular-singular SCC scalar execution does not yet implement pseudo-resonant compensation");
          if (seed && n == 0 && kind != "R")
            throw std::invalid_argument(
                "native regular-singular SCC seed must begin at its exact indicial root");
        } else {
          const auto expected_kind = n == 0 ? "R" : "T";
          if (kind != expected_kind || !db.is_zero() ||
              !(da == Rational(std::to_string(n))))
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
      if (dimension != 1 || unit_index_i64 < 0 ||
          unit_index_i64 >= frame_width ||
          initial.size() != frame_width || validity.size() != 1 ||
          validity.front().is_null() ||
          as_i32(validity.front(), "native SCC seed validity") !=
              work_.work_complete_max)
        throw std::invalid_argument(
            "native regular-singular SCC seed requires one honest scalar eps^0 unit frame");
      const auto unit_index = static_cast<std::size_t>(unit_index_i64);
      for (std::size_t epsilon = 0; epsilon < frame_width; ++epsilon) {
        const auto coefficient = parse_scalar<Rational>(initial[epsilon]);
        if ((epsilon == unit_index && !(coefficient == Rational(1))) ||
            (epsilon != unit_index && !coefficient.is_zero()))
          throw std::invalid_argument(
              "native regular-singular SCC seed is not the captured canonical eps^0 scalar normalization");
      }
      return 0;
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
    std::optional<std::uint32_t> selected;
    for (std::uint32_t component = 0; component < dimension; ++component) {
      for (std::size_t epsilon = 0; epsilon < frame_width; ++epsilon) {
        const auto coefficient = parse_scalar<Rational>(
            initial[static_cast<std::size_t>(component) * frame_width +
                    epsilon]);
        if (epsilon == unit_index && coefficient == Rational(1)) {
          if (selected.has_value())
            throw std::invalid_argument(
                "native SCC seed contains more than one eps^0 unit component");
          selected = component;
        } else if (!coefficient.is_zero()) {
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
                             std::uint32_t block_index, bool seed) const {
    if (!native.pseudo_hits.empty())
      throw std::invalid_argument(
          "native SCC regular block execution encountered unsupported pseudo hits");
    if (block_index >= blocks_.size() ||
        native.solution.dimension != blocks_[block_index].vertices.size())
      throw std::invalid_argument(
          "native SCC block solve returned the wrong retained dimension");
    require_work_local(native.solution, "native SCC block result");
    if (native.solution.sectors.empty())
      throw std::invalid_argument("native SCC block solve returned no sectors");
    if (seed && (native.solution.sectors.size() != 1 ||
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
    const auto a_target = parse_scalar<Rational>(run.at("a_target"));
    const auto b_target = parse_scalar<Rational>(run.at("b_target"));
    const auto log_max = as_u32(run.at("p"), "target log maximum");
    for (const auto& sector : source.sectors) {
      if (sector.a.domain != ExactDomain::Rational ||
          sector.b.domain != ExactDomain::Rational ||
          !(Rational(sector.a.canonical) == a_target) ||
          !(Rational(sector.b.canonical) == b_target) ||
          sector.log_power > log_max)
        throw std::invalid_argument(
            "native SCC coupling source tag/log sector differs from its target run");
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
        {"pseudo_hit_count", result.pseudo_hits.size()},
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
  std::uint64_t next_chart = 1;
  std::uint64_t next_local = 1;
  std::uint64_t next_scc = 1;
  std::uint64_t next_match = 1;
  std::uint64_t next_endpoint = 1;
  std::size_t pending_local_solves = 0;
  std::size_t pending_matches = 0;
  std::size_t pending_endpoint_limits = 0;
  std::uint64_t total_local_solves = 0;
  std::uint64_t total_scc_column_solves = 0;
  std::uint64_t total_local_matches = 0;
  std::uint64_t total_endpoint_limits = 0;
  std::uint64_t total_endpoint_exports = 0;
  double total_local_run_parse_ms = 0.0;
  double total_local_kernel_ms = 0.0;
  double total_local_match_ms = 0.0;
  std::uint64_t checkpoint_generation = 0;
  std::uint64_t checkpoint_restore_count = 0;
  std::string restored_from_checkpoint_identity;
  double total_endpoint_limit_ms = 0.0;
  double total_endpoint_export_ms = 0.0;
  bool closed = false;
  mutable std::mutex mutex;
  std::unordered_map<std::string, std::shared_ptr<PreparedChartBase>> charts;
  std::unordered_map<std::string, std::string> handles_by_key;
  std::unordered_map<std::string, std::shared_ptr<StoredLocalBase>> locals;
  std::unordered_map<std::string, std::shared_ptr<StoredMatchBase>> matches;
  std::unordered_map<std::string, std::shared_ptr<StoredEndpointResult>>
      endpoints;
  std::unordered_map<std::string, std::shared_ptr<CompositeSCCChartBase>> sccs;
  std::unordered_map<std::string, std::string> scc_handles_by_key;
};

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
    if (!raw_block.at("identity_gauge").as_bool())
      throw std::invalid_argument(
          "native SCC preparation requires an exact identity gauge");
    if (!raw_block.at("identity_v").as_bool())
      throw std::invalid_argument(
          "native SCC preparation requires an exact identity spectral transform");
    if (!raw_block.at("no_pseudo").as_bool())
      throw std::invalid_argument(
          "native SCC preparation does not yet execute pseudo compensation");

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
        regular, std::move(chart)};
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
constexpr std::uint32_t kCheckpointPayloadSchema = 1;

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
      {"endpoint_capacity", session.endpoint_capacity}};
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

struct SessionCheckpointSnapshot {
  json::object header;
  json::object payload;
  std::uint64_t generation = 0;
  std::size_t charts = 0;
  std::size_t sccs = 0;
};

SessionCheckpointSnapshot make_checkpoint_snapshot(
    SolverSession& session, const std::string& checkpoint_identity) {
  if (session.closed)
    throw std::invalid_argument("cannot checkpoint a closed solver session");
  if (session.pending_local_solves != 0 || session.pending_matches != 0 ||
      session.pending_endpoint_limits != 0)
    throw std::invalid_argument(
        "checkpoint requires a quiescent session with no pending local solve, match, or endpoint limit");
  if (!session.locals.empty() || !session.matches.empty() ||
      !session.endpoints.empty())
    throw std::invalid_argument(
        "checkpoint schema v1 does not serialize retained local, match, or endpoint handles; release them before saving");

  std::vector<std::shared_ptr<PreparedChartBase>> charts;
  charts.reserve(session.charts.size());
  for (const auto& [ignored, chart] : session.charts) charts.push_back(chart);
  std::sort(charts.begin(), charts.end(), [](const auto& left,
                                             const auto& right) {
    return scoped_handle_id(left->handle(), "c:", "chart") <
           scoped_handle_id(right->handle(), "c:", "chart");
  });
  std::vector<std::shared_ptr<CompositeSCCChartBase>> sccs;
  sccs.reserve(session.sccs.size());
  for (const auto& [ignored, composite] : session.sccs)
    sccs.push_back(composite);
  std::sort(sccs.begin(), sccs.end(), [](const auto& left,
                                         const auto& right) {
    return scoped_handle_id(left->handle(), "scc:", "SCC") <
           scoped_handle_id(right->handle(), "scc:", "SCC");
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
      {"total_local_solves", session.total_local_solves},
      {"total_scc_column_solves", session.total_scc_column_solves},
      {"total_local_matches", session.total_local_matches},
      {"total_endpoint_limits", session.total_endpoint_limits},
      {"total_endpoint_exports", session.total_endpoint_exports},
      {"total_local_run_parse_ms", session.total_local_run_parse_ms},
      {"total_local_kernel_ms", session.total_local_kernel_ms},
      {"total_local_match_ms", session.total_local_match_ms},
      {"total_endpoint_limit_ms", session.total_endpoint_limit_ms},
      {"total_endpoint_export_ms", session.total_endpoint_export_ms},
      {"checkpoint_generation", generation},
      {"checkpoint_restore_count", session.checkpoint_restore_count}};
  json::object session_record{
      {"source_handle", session.handle},
      {"configuration", configuration},
      {"configuration_identity", checkpoint_configuration_identity(session)},
      {"counters", std::move(counters)}};
  json::object payload{
      {"schema", kCheckpointPayloadSchema},
      {"session", std::move(session_record)},
      {"prepared_charts", chart_items},
      {"prepared_scc", scc_items}};
  json::array mandatory_sections{"session", "prepared_charts",
                                  "prepared_scc"};
  json::array deferred_kinds{"local", "match", "endpoint", "line", "tile"};
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
      {"generation", generation}};
  return {std::move(header), std::move(payload), generation,
          charts.size(), sccs.size()};
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
       "scc_identities", "generation"}, "checkpoint header");
  require_exact_keys(payload,
      {"schema", "session", "prepared_charts", "prepared_scc"},
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
      "prepared_charts", "prepared_scc", "session"};
  if (mandatory != expected_sections)
    throw std::invalid_argument(
        "native checkpoint contains unknown or missing mandatory sections");
  if (!as_array(header.at("optional_sections"),
                "optional checkpoint sections").empty())
    throw std::invalid_argument(
        "native checkpoint declares unsupported optional sections");

  const auto& session = as_object(payload.at("session"),
                                  "checkpoint session section");
  require_exact_keys(session,
      {"source_handle", "configuration", "configuration_identity",
       "counters"}, "checkpoint session section");
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
  if (checkpoint_identity_manifest(chart_items) !=
          as_array(header.at("chart_identities"),
                   "checkpoint chart identities") ||
      checkpoint_identity_manifest(scc_items) !=
          as_array(header.at("scc_identities"),
                   "checkpoint SCC identities"))
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
       "match_capacity", "endpoint_capacity"},
      "checkpoint configuration");
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
      {"endpoint_capacity", configuration.at("endpoint_capacity")}};
  const auto created = run_session_command(create);
  const auto restored_handle = required_string(created, "session");
  bool live = true;
  try {
    const auto restored = find_session(restored_handle);
    const auto source_handle = required_string(saved_session, "source_handle");
    json::array restored_charts;
    std::uint64_t largest_chart = 0;
    for (const auto& raw_item : as_array(
             payload.at("prepared_charts"), "checkpoint prepared charts")) {
      const auto& item = as_object(raw_item, "checkpoint chart item");
      require_exact_keys(item,
          {"handle", "key", "identity", "signature", "request"},
          "checkpoint chart item");
      const auto old_handle = required_string(item, "handle");
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
      restored_charts.push_back(json::object{
          {"chart", old_handle}, {"key", item.at("key")},
          {"identity", item.at("identity")}});
    }

    json::array restored_sccs;
    std::uint64_t largest_scc = 0;
    for (const auto& raw_item : as_array(
             payload.at("prepared_scc"), "checkpoint prepared SCC charts")) {
      const auto& item = as_object(raw_item, "checkpoint SCC item");
      require_exact_keys(item,
          {"handle", "key", "identity", "signature", "request"},
          "checkpoint SCC item");
      const auto old_handle = required_string(item, "handle");
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
      restored_sccs.push_back(json::object{
          {"scc", old_handle}, {"key", item.at("key")},
          {"identity", item.at("identity")}});
    }

    const auto& counters = as_object(saved_session.at("counters"),
                                     "checkpoint counters");
    require_exact_keys(counters,
        {"next_chart", "next_local", "next_scc", "next_match",
         "next_endpoint",
         "total_local_solves", "total_scc_column_solves",
         "total_local_matches", "total_endpoint_limits",
         "total_endpoint_exports", "total_local_run_parse_ms",
         "total_local_kernel_ms", "total_local_match_ms",
         "total_endpoint_limit_ms", "total_endpoint_export_ms",
         "checkpoint_generation", "checkpoint_restore_count"},
        "checkpoint counters");
    const auto next_chart = as_u64(counters.at("next_chart"), "next chart");
    const auto next_scc = as_u64(counters.at("next_scc"), "next SCC");
    const auto next_local = as_u64(counters.at("next_local"), "next local");
    const auto next_match = as_u64(counters.at("next_match"), "next match");
    const auto next_endpoint = as_u64(
        counters.at("next_endpoint"), "next endpoint");
    const auto restore_count = as_u64(
        counters.at("checkpoint_restore_count"),
        "checkpoint restore count");
    if (next_chart <= largest_chart || next_scc <= largest_scc ||
        next_local == 0 || next_match == 0 || next_endpoint == 0)
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
                        {"local_match_capability",
                         domain == "rational"
                             ? kExactRegularLocalMatchCapability
                             : "unsupported"},
                        {"endpoint_limit_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedEndpointLimitCapability}};
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
    std::size_t charts = 0, locals = 0, matches = 0, endpoints = 0, sccs = 0;
    {
      std::lock_guard<std::mutex> lock(removed->mutex);
      removed->closed = true;
      // In-flight solve/match/endpoint calls own their reservations and
      // decrement them on exactly one completion path.  Do not reset pending
      // counters.
      charts = removed->charts.size();
      locals = removed->locals.size();
      matches = removed->matches.size();
      endpoints = removed->endpoints.size();
      sccs = removed->sccs.size();
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
        {"serialized_handle_kinds", json::array{"chart", "scc"}},
        {"deferred_handle_kinds",
         json::array{"local", "match", "endpoint", "line", "tile"}},
        {"atomic", true}};
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
    std::vector<std::shared_ptr<CompositeSCCChartBase>> sccs;
    std::size_t pending_local_solves = 0;
    std::size_t pending_matches = 0;
    std::size_t pending_endpoint_limits = 0;
    std::uint64_t total_local_solves = 0;
    std::uint64_t total_scc_column_solves = 0;
    std::uint64_t total_local_matches = 0;
    std::uint64_t checkpoint_generation = 0;
    std::uint64_t checkpoint_restore_count = 0;
    std::string restored_from_checkpoint_identity;
    std::uint64_t total_endpoint_limits = 0;
    std::uint64_t total_endpoint_exports = 0;
    double total_local_run_parse_ms = 0.0, total_local_kernel_ms = 0.0;
    double total_local_match_ms = 0.0;
    double total_endpoint_limit_ms = 0.0;
    double total_endpoint_export_ms = 0.0;
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
      for (const auto& [ignored, composite] : session->sccs)
        sccs.push_back(composite);
      pending_local_solves = session->pending_local_solves;
      pending_matches = session->pending_matches;
      pending_endpoint_limits = session->pending_endpoint_limits;
      total_local_solves = session->total_local_solves;
      total_scc_column_solves = session->total_scc_column_solves;
      total_local_matches = session->total_local_matches;
      total_endpoint_limits = session->total_endpoint_limits;
      total_endpoint_exports = session->total_endpoint_exports;
      total_local_run_parse_ms = session->total_local_run_parse_ms;
      total_local_kernel_ms = session->total_local_kernel_ms;
      total_local_match_ms = session->total_local_match_ms;
      checkpoint_generation = session->checkpoint_generation;
      checkpoint_restore_count = session->checkpoint_restore_count;
      restored_from_checkpoint_identity =
          session->restored_from_checkpoint_identity;
      total_endpoint_limit_ms = session->total_endpoint_limit_ms;
      total_endpoint_export_ms = session->total_endpoint_export_ms;
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
    std::size_t local_coefficients = 0;
    double local_evaluate_ms = 0.0, local_residual_certify_ms = 0.0;
    double local_endpoint_limit_ms = 0.0;
    json::array local_stats;
    for (const auto& local : locals) {
      const auto stats = local->stats();
      local_evaluations += stats.evaluations;
      local_residual_certifications += stats.residual_certifications;
      local_endpoint_limits += stats.endpoint_limits;
      local_coefficients += stats.coefficient_count;
      local_evaluate_ms += stats.evaluate_ms;
      local_residual_certify_ms += stats.residual_certify_ms;
      local_endpoint_limit_ms += stats.endpoint_limit_ms;
      local_stats.push_back(local->stats_json());
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
    return json::object{{"status", "ok"}, {"session", session->handle},
                        {"domain", session->domain},
                        {"precision_bits", session->precision_bits},
                        {"charts", charts.size()}, {"runs", runs},
                        {"locals", locals.size()},
                        {"matches", matches.size()},
                        {"endpoints", endpoints.size()},
                        {"scc_charts", sccs.size()},
                        {"pending_local_solves", pending_local_solves},
                        {"pending_matches", pending_matches},
                        {"pending_endpoint_limits", pending_endpoint_limits},
                        {"local_solves", total_local_solves},
                        {"local_matches", total_local_matches},
                        {"endpoint_limits", total_endpoint_limits},
                        {"endpoint_exports", total_endpoint_exports},
                        {"local_match_capability",
                         session->domain == "rational"
                             ? kExactRegularLocalMatchCapability
                             : "unsupported"},
                        {"endpoint_limit_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedEndpointLimitCapability},
                        {"scc_column_solves", total_scc_column_solves},
                        {"local_evaluations", local_evaluations},
                        {"local_residual_certifications",
                         local_residual_certifications},
                        {"local_endpoint_limits", local_endpoint_limits},
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
                        {"chart_stats", std::move(chart_stats)},
                        {"local_stats", std::move(local_stats)},
                        {"match_stats", std::move(match_stats)},
                        {"endpoint_stats", std::move(endpoint_stats)},
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
                                      {"persistent_symbolic_endpoint_limits",
                                       false},
                                      {"persistent_acb_local_match", false},
                                      {"persistent_checkpoint", true},
                                      {"persistent_checkpoint_schema", 1},
                                      {"persistent_checkpoint_handle_scope",
                                       "prepared-chart-and-scc"},
                                      {"persistent_scc_regular_singular_scalar_block_dag_column",
                                       true},
                                      {"backend", "DiffExp2 C++"},
                                      {"flint", flint_version},
                                      {"librarylink", true}});
}

}  // namespace diffexp2
