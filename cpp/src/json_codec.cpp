#include "diffexp2/json_codec.hpp"

#include "diffexp2/local_algebra.hpp"
#include "diffexp2/local_solution.hpp"
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

struct StoredLocalStats {
  std::uint64_t evaluations = 0;
  double evaluate_ms = 0.0;
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
      const auto& point_object = as_object(
          request.at("point"), "local evaluation point");
      const auto point = RealEvaluationPoint::rational(
          required_string(point_object, "exact"));
      EvaluationOptions options;
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
    out["evaluate_ms"] = current.evaluate_ms;
    return out;
  }

  StoredLocalStats stats() const override {
    std::lock_guard<std::mutex> lock(stats_mutex_);
    return {evaluations_.load(), evaluate_ms_, create_parse_ms_,
            create_kernel_ms_, coefficient_count()};
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
    const auto& sector = solution_.sectors.front();
    return json::object{
        {"chart", json::object{
            {"center_exact", solution_.chart.center_exact},
            {"scale_exact", solution_.chart.scale_exact},
            {"infinite_radius", solution_.chart.infinite_radius},
            {"radius_ball", encode_scalar(solution_.chart.radius, 30)}}},
        {"tag", json::object{{"a", encode_exact_descriptor(sector.a)},
                             {"b", encode_exact_descriptor(sector.b)}}},
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
  mutable std::mutex stats_mutex_;
  double evaluate_ms_ = 0.0;
};

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
  slong precision_bits() const { return precision_bits_; }
  bool has_identity_assembly() const {
    return prepared_.assembly_matrix.has_value() &&
           prepared_.assembly_matrix->identity;
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
  if (source.dimension != 1 || !source.error.empty())
    throw std::invalid_argument(
        "native SCC source injection requires one uncertified scalar local");
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
  if (points > std::numeric_limits<std::size_t>::max() / frame_width)
    throw std::overflow_error("native SCC source tensor size overflow");
  SourceData<Scalar> data;
  data.frames.assign(points * frame_width, ScalarTraits<Scalar>::zero());
  data.validity.assign(points, kCompleteInfinity);
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
      // bound propagates the honest source window through the recurrence.
      data.validity[point] = source.epsilon.complete_max;
      for (std::int64_t power = source.epsilon.min_power;
           power <= source.epsilon.complete_max; ++power) {
        const auto input_epsilon = static_cast<std::size_t>(
            power - source.epsilon.min_power);
        const auto output_epsilon = static_cast<std::size_t>(
            power - frame_base);
        data.frames[point * frame_width + output_epsilon] =
            sector.coefficients[local_algebra_detail::flat_index(
                input_epsilon, n, 0, source.taylor_width(), 1)];
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
      if (!scalar_column_ready())
        throw std::invalid_argument(
            "retained SCC chart does not satisfy the exact-rational regular scalar-DAG column capability");
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
          seed_request, true, nullptr);
      auto seed_native = blocks_[seed_block].chart->solve_native(
          seed_run, as_object(seed_request.at("metadata"),
                              "native SCC seed metadata"));
      validate_block_result(seed_native, true);
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
          validate_scalar_coupling(coupling);
          auto contribution = apply_prepared_sparse_local_matrix(
              coupling.matrix, *state[coupling.source_block],
              checkpoint_identity + ":source:" +
                  std::to_string(coupling.source_block) + ":" +
                  std::to_string(target_block));
          if (!contribution.has_value())
            throw std::logic_error(
                "an exact nonzero scalar SCC edge produced no structural source");
          *contribution = restrict_local_epsilon_frame_strict_lower(
              *contribution, work_.work_min, work_.work_complete_max,
              checkpoint_identity + ":source-frame:" +
                  std::to_string(coupling.source_block) + ":" +
                  std::to_string(target_block));
          require_work_local(*contribution, "coupling source");
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
        require_work_local(source, "combined coupling source");
        const auto& target_run = checked_column_run(
            target_request, false, &source);
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
        validate_block_result(target_native, false);
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
        embedded.push_back(local_algebra_detail::embedded_component(
            *state[block], blocks_[block].vertices.front(), dimension_));
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
      json::object column_identity_record{
          {"schema", "diffexp2-native-scc-column-v1"},
          {"scc_exact_identity", exact_identity_},
          {"basis_index", blocks_[seed_block].vertices.front()},
          {"seed", seed_request},
          {"targets", target_requests}};
      SCCColumnProvenance column_provenance{
          handle_, exact_identity_, seed_block,
          blocks_[seed_block].vertices.front(),
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

  json::object stats_json() const override {
    std::size_t active_entries = 0, proven_zero_entries = 0;
    std::optional<std::int32_t> min_coupling_shift;
    std::optional<std::int32_t> max_coupling_shift;
    json::array block_handles;
    block_handles.reserve(blocks_.size());
    for (const auto& block : blocks_) {
      block_handles.push_back(json::object{
          {"block", block.block}, {"chart", block.source_handle},
          {"dimension", block.chart->dimension()},
          {"principal_identity", block.principal_identity}});
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
    return json::object{
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
        {"execution_implemented", scalar_column_ready()},
        {"execution_scope",
         "exact-rational-regular-scalar-block-dag-column-v1"},
        {"general_scc_execution", false},
        {"scalar_block_dag_column_execution", scalar_column_ready()},
        {"scc_column_solves", column_solves_.load()},
        {"scc_column_solve_ms", column_solve_ms()},
        {"capability_evidence", json::object{
             {"identity_v", "native-retained-assembly"},
             {"regular", "collision-bound-producer-certificate"},
             {"identity_gauge", "collision-bound-producer-certificate"},
             {"no_pseudo", "collision-bound-producer-certificate"}}},
        {"execution_must_revalidate_producer_capabilities", true},
        {"block_charts", std::move(block_handles)}};
  }

 private:
  bool scalar_column_ready() const {
    if constexpr (!std::is_same_v<Scalar, Rational>) {
      return false;
    } else {
      if (blocks_.size() < 2 || work_.work_min > 0 ||
          work_.requested_max < 0 || work_.work_complete_max < 0 ||
          std::any_of(blocks_.begin(), blocks_.end(), [](const auto& block) {
            return block.vertices.size() != 1 ||
                   block.chart->dimension() != 1 ||
                   !block.chart->has_identity_assembly();
          }))
        return false;
      return std::all_of(
          couplings_.begin(), couplings_.end(), [](const auto& coupling) {
            return coupling.matrix.rows == 1 &&
                   coupling.matrix.columns == 1 &&
                   coupling.matrix.entries.size() == 1 &&
                   coupling.matrix.entries.front().row == 0 &&
                   coupling.matrix.entries.front().column == 0 &&
                   !coupling.matrix.entries.front().multiplier.proven_zero &&
                   coupling.matrix.entries.front()
                           .multiplier.center_pole_order == 0 &&
                   !coupling.matrix.entries.front()
                        .multiplier.kernels.empty() &&
                   std::all_of(
                       coupling.matrix.entries.front()
                           .multiplier.kernels.begin(),
                       coupling.matrix.entries.front()
                           .multiplier.kernels.end(),
                       [](const auto& kernel) {
                         return !kernel.empty() &&
                                ScalarTraits<Scalar>::is_zero(kernel.front());
                       });
          });
    }
  }

  const json::object& checked_column_run(
      const json::object& entry, bool seed,
      const LocalSolution<Scalar>* source) const {
    const auto& run = as_object(entry.at("run"), "native SCC recurrence run");
    validate_metadata_geometry(as_object(
        entry.at("metadata"), "native SCC local metadata"));
    if (!run.at("source").is_null())
      throw std::invalid_argument(
          "native SCC column rejects caller-supplied recurrence source data");
    if (run.at("adaptive_probe").as_bool())
      throw std::invalid_argument(
          "native SCC regular scalar column requires a fixed retained lower frame");
    if (run.at("return_u").as_bool())
      throw std::invalid_argument(
          "native SCC column requires retained assembly without U JSON");
    const auto nmax = as_u32(run.at("nmax"), "native SCC Taylor order");
    const auto log_max = as_u32(run.at("p"), "native SCC log maximum");
    if (nmax != work_.work_t_order)
      throw std::invalid_argument(
          "native SCC run Taylor order differs from its retained work contract");
    if (log_max != 0 ||
        !parse_scalar<Rational>(run.at("a_target")).is_zero() ||
        !parse_scalar<Rational>(run.at("b_target")).is_zero() ||
        as_i32(run.at("a_shift_min"), "native SCC a-shift minimum") != 0)
      throw std::invalid_argument(
          "native SCC regular scalar runs require p=0, a=b=0 and zero a-shift origin");
    const auto& a_shifts = as_array(
        run.at("a_shifts"), "native SCC exact a-shift schedule");
    if (a_shifts.size() != static_cast<std::size_t>(nmax) + 1)
      throw std::invalid_argument(
          "native SCC regular scalar run has an incomplete a-shift schedule");
    for (std::size_t n = 0; n < a_shifts.size(); ++n)
      if (!(parse_scalar<Rational>(a_shifts[n]) ==
            Rational(std::to_string(n))))
        throw std::invalid_argument(
            "native SCC regular scalar a-shift schedule must equal the Taylor index");
    if (seed) {
      if (!run.at("has_initial").as_bool())
        throw std::invalid_argument(
            "native SCC seed requires one initialized log-zero sector");
      const auto& initial = as_array(
          run.at("initial"), "native SCC seed initial tensor");
      const auto& validity = as_array(
          run.at("initial_validity"),
          "native SCC seed initial validity");
      const auto frame_width = blocks_.front().chart->frame_width();
      const auto unit_index_i64 = -static_cast<std::int64_t>(work_.work_min);
      if (unit_index_i64 < 0 || unit_index_i64 >= frame_width ||
          initial.size() != frame_width || validity.size() != 1 ||
          validity.front().is_null() ||
          as_i32(validity.front(), "native SCC seed validity") !=
              work_.work_complete_max)
        throw std::invalid_argument(
            "native SCC regular scalar seed requires one honest finite eps^0 unit frame");
      const auto unit_index = static_cast<std::size_t>(unit_index_i64);
      for (std::size_t index = 0; index < initial.size(); ++index) {
        const auto coefficient = parse_scalar<Rational>(initial[index]);
        if ((index == unit_index && coefficient != Rational(1)) ||
            (index != unit_index && !coefficient.is_zero()))
          throw std::invalid_argument(
              "native SCC regular scalar seed must be the exact eps^0 unit column");
      }
    } else {
      if (run.at("has_initial").as_bool())
        throw std::invalid_argument(
            "native SCC target rejects caller-supplied initial state");
      const auto& initial = as_array(
          run.at("initial"), "native SCC target initial tensor");
      const auto& validity = as_array(
          run.at("initial_validity"),
          "native SCC target initial validity");
      const auto expected_initial =
          static_cast<std::size_t>(log_max) + 1;
      if (expected_initial > std::numeric_limits<std::size_t>::max() /
                                 blocks_.front().chart->frame_width())
        throw std::overflow_error(
            "native SCC target initial tensor size overflow");
      const auto expected_initial_coefficients =
          expected_initial * blocks_.front().chart->frame_width();
      if (initial.size() != expected_initial_coefficients ||
          validity.size() != static_cast<std::size_t>(log_max) + 1 ||
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
      if (row.size() != 1)
        throw std::invalid_argument(
            "native SCC scalar execution requires one Jordan block step");
      const auto& step = as_object(row.front(), "native SCC schedule step");
      const auto kind = required_string(step, "case");
      const auto da = parse_scalar<Rational>(step.at("da"));
      const auto db = parse_scalar<Rational>(step.at("db"));
      const auto expected_kind = n == 0 ? "R" : "T";
      if (kind != expected_kind || !db.is_zero() ||
          !(da == Rational(std::to_string(n))))
        throw std::invalid_argument(
            "native SCC regular scalar schedule must be resonant at zero and Taylor by exact index thereafter");
    }
    return run;
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
                             bool seed) const {
    if (!native.pseudo_hits.empty())
      throw std::invalid_argument(
          "native SCC scalar execution encountered unsupported pseudo hits");
    if (native.solution.dimension != 1)
      throw std::invalid_argument(
          "native SCC block solve did not return one scalar component");
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
            "native SCC first execution slice requires exact rational sector tags");
  }

  void validate_scalar_coupling(
      const CompositeSCCCoupling<Scalar>& coupling) const {
    if (coupling.matrix.rows != 1 || coupling.matrix.columns != 1 ||
        coupling.matrix.entries.size() != 1 ||
        coupling.matrix.entries.front().row != 0 ||
        coupling.matrix.entries.front().column != 0 ||
        coupling.matrix.entries.front().multiplier.proven_zero ||
        coupling.matrix.entries.front().multiplier.center_pole_order != 0 ||
        coupling.matrix.entries.front().multiplier.kernels.empty() ||
        std::any_of(
            coupling.matrix.entries.front().multiplier.kernels.begin(),
            coupling.matrix.entries.front().multiplier.kernels.end(),
            [](const auto& kernel) {
              return kernel.empty() ||
                     !ScalarTraits<Scalar>::is_zero(kernel.front());
            }))
      throw std::invalid_argument(
          "native SCC column requires one exact nonzero pole-free scalar coupling vanishing at chart center per edge");
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
  std::uint64_t next_chart = 1;
  std::uint64_t next_local = 1;
  std::uint64_t next_scc = 1;
  std::size_t pending_local_solves = 0;
  std::uint64_t total_local_solves = 0;
  std::uint64_t total_scc_column_solves = 0;
  double total_local_run_parse_ms = 0.0;
  double total_local_kernel_ms = 0.0;
  bool closed = false;
  mutable std::mutex mutex;
  std::unordered_map<std::string, std::shared_ptr<PreparedChartBase>> charts;
  std::unordered_map<std::string, std::string> handles_by_key;
  std::unordered_map<std::string, std::shared_ptr<StoredLocalBase>> locals;
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
    if (!raw_block.at("regular").as_bool())
      throw std::invalid_argument(
          "native SCC preparation does not yet support singular blocks");
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
        std::move(chart)};
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
                        {"scc_capacity", session->scc_capacity}};
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
    std::size_t charts = 0, locals = 0, sccs = 0;
    {
      std::lock_guard<std::mutex> lock(removed->mutex);
      removed->closed = true;
      // In-flight local.solve calls own these reservations and decrement them
      // on exactly one completion path.  Do not reset pending_local_solves.
      charts = removed->charts.size();
      locals = removed->locals.size();
      sccs = removed->sccs.size();
      removed->locals.clear();
      removed->sccs.clear();
      removed->scc_handles_by_key.clear();
      removed->charts.clear();
      removed->handles_by_key.clear();
    }
    return json::object{{"status", "ok"}, {"closed", handle},
                        {"released_charts", charts},
                        {"released_locals", locals},
                        {"released_scc_charts", sccs}};
  }

  const auto session = find_session(required_string(root, "session"));

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
                            {"frame_width", chart->frame_width()}};
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
        "exact-rational-regular-scalar-block-dag-column-v1";
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
    std::vector<std::shared_ptr<CompositeSCCChartBase>> sccs;
    std::size_t pending_local_solves = 0;
    std::uint64_t total_local_solves = 0;
    std::uint64_t total_scc_column_solves = 0;
    double total_local_run_parse_ms = 0.0, total_local_kernel_ms = 0.0;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      for (const auto& [ignored, chart] : session->charts)
        charts.push_back(chart);
      for (const auto& [ignored, local] : session->locals)
        locals.push_back(local);
      for (const auto& [ignored, composite] : session->sccs)
        sccs.push_back(composite);
      pending_local_solves = session->pending_local_solves;
      total_local_solves = session->total_local_solves;
      total_scc_column_solves = session->total_scc_column_solves;
      total_local_run_parse_ms = session->total_local_run_parse_ms;
      total_local_kernel_ms = session->total_local_kernel_ms;
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
    std::size_t local_coefficients = 0;
    double local_evaluate_ms = 0.0;
    json::array local_stats;
    for (const auto& local : locals) {
      const auto stats = local->stats();
      local_evaluations += stats.evaluations;
      local_coefficients += stats.coefficient_count;
      local_evaluate_ms += stats.evaluate_ms;
      local_stats.push_back(local->stats_json());
    }
    json::array scc_stats;
    for (const auto& composite : sccs)
      scc_stats.push_back(composite->stats_json());
    return json::object{{"status", "ok"}, {"session", session->handle},
                        {"domain", session->domain},
                        {"precision_bits", session->precision_bits},
                        {"charts", charts.size()}, {"runs", runs},
                        {"locals", locals.size()},
                        {"scc_charts", sccs.size()},
                        {"pending_local_solves", pending_local_solves},
                        {"local_solves", total_local_solves},
                        {"scc_column_solves", total_scc_column_solves},
                        {"local_evaluations", local_evaluations},
                        {"local_coefficient_count", local_coefficients},
                        {"static_tensor_copies", 0},
                        {"prepare_parse_ms", prepare_parse_ms},
                        {"run_parse_ms", run_parse_ms},
                        {"kernel_ms", kernel_ms},
                        {"local_run_parse_ms", total_local_run_parse_ms},
                        {"local_kernel_ms", total_local_kernel_ms},
                        {"local_evaluate_ms", local_evaluate_ms},
                        {"chart_stats", std::move(chart_stats)},
                        {"local_stats", std::move(local_stats)},
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
                                      {"persistent_scc_prepare", true},
                                      {"persistent_scc_execute", false},
                                      {"persistent_scc_scalar_block_dag_column",
                                       true},
                                      {"backend", "DiffExp2 C++"},
                                      {"flint", flint_version},
                                      {"librarylink", true}});
}

}  // namespace diffexp2
