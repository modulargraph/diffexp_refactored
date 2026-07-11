#include "diffexp2/json_codec.hpp"

#include "diffexp2/recurrence.hpp"

#include <boost/json.hpp>

#include <chrono>
#include <atomic>
#include <algorithm>
#include <condition_variable>
#include <cstdint>
#include <functional>
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
  if (!problem.return_u && !problem.assembly_matrix.has_value())
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

struct ChartStats {
  std::uint64_t runs = 0;
  double prepare_parse_ms = 0.0;
  double run_parse_ms = 0.0;
  double kernel_ms = 0.0;
};

class PreparedChartBase {
 public:
  PreparedChartBase(std::string handle, std::string key,
                    std::string signature, SCCCertificate scc,
                    double prepare_parse_ms)
      : handle_(std::move(handle)), key_(std::move(key)),
        signature_(std::move(signature)), scc_(std::move(scc)),
        prepare_parse_ms_(prepare_parse_ms) {}
  virtual ~PreparedChartBase() = default;

  virtual json::object solve(const json::object& run, int output_digits) = 0;
  virtual std::uint32_t dimension() const = 0;
  virtual std::int32_t frame_base() const = 0;
  virtual std::uint32_t frame_width() const = 0;
  virtual ChartStats stats() const = 0;

  const std::string& handle() const { return handle_; }
  const std::string& key() const { return key_; }
  const std::string& signature() const { return signature_; }
  const SCCCertificate& scc() const { return scc_; }

 protected:
  std::string handle_;
  std::string key_;
  std::string signature_;
  SCCCertificate scc_;
  double prepare_parse_ms_ = 0.0;
};

std::mutex& symbolic_run_mutex() {
  static std::mutex mutex;
  return mutex;
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

template <typename Scalar>
class PreparedChart final : public PreparedChartBase {
 public:
  PreparedChart(std::string handle, std::string key, std::string signature,
                SCCCertificate scc,
                PreparedRecurrenceOperator<Scalar>&& prepared,
                slong precision_bits, std::vector<std::string> symbols,
                double prepare_parse_ms)
      : PreparedChartBase(std::move(handle), std::move(key),
                          std::move(signature), std::move(scc),
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
    std::unique_lock<std::mutex> symbolic_lock;
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      symbolic_lock = std::unique_lock<std::mutex>(symbolic_run_mutex());
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

  std::uint32_t dimension() const override { return prepared_.dimension; }
  std::int32_t frame_base() const override { return prepared_.frame_base; }
  std::uint32_t frame_width() const override { return prepared_.frame_width; }
  ChartStats stats() const override {
    std::lock_guard<std::mutex> lock(stats_mutex_);
    return {runs_.load(), prepare_parse_ms_, run_parse_ms_, kernel_ms_};
  }

 private:
  PreparedRecurrenceOperator<Scalar> prepared_;
  slong precision_bits_ = 256;
  std::vector<std::string> symbols_;
  std::atomic<std::uint64_t> runs_{0};
  mutable std::mutex stats_mutex_;
  double run_parse_ms_ = 0.0;
  double kernel_ms_ = 0.0;
};

struct SolverSession {
  std::string handle;
  std::string domain;
  slong precision_bits = 256;
  int output_digits = 50;
  std::vector<std::string> symbols;
  std::string analytic_identity;
  std::size_t chart_capacity = 256;
  std::uint64_t next_chart = 1;
  mutable std::mutex mutex;
  std::unordered_map<std::string, std::shared_ptr<PreparedChartBase>> charts;
  std::unordered_map<std::string, std::string> handles_by_key;
};

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
    SCCCertificate scc, std::string signature) {
  const auto started = std::chrono::steady_clock::now();
  const auto& problem = as_object(root.at("problem"), "prepared problem");
  std::unique_ptr<AcbPrecisionLease> acb_lease;
  if constexpr (std::is_same_v<Scalar, ComplexBall>) {
    acb_lease = std::make_unique<AcbPrecisionLease>(session->precision_bits);
    ComplexBall::set_precision(session->precision_bits);
  }
  std::unique_lock<std::mutex> symbolic_lock;
  if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
    symbolic_lock = std::unique_lock<std::mutex>(symbolic_run_mutex());
    SymbolicRational::configure(session->symbols);
  }
  auto prepared = parse_prepared_operator<Scalar>(problem);
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<PreparedChart<Scalar>>(
      handle, key, std::move(signature), std::move(scc),
      std::move(prepared), session->precision_bits, session->symbols, elapsed);
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

    auto& registry = session_registry();
    {
      std::lock_guard<std::mutex> lock(registry.mutex);
      // ComplexBall's configured precision is thread-local in this backend,
      // so equal-precision sessions can execute concurrently.  Retaining one
      // precision across all live sessions also protects any future FLINT
      // helper that acquires process-global precision state.
      if (domain == "acb") {
        for (const auto& [ignored, existing] : registry.sessions) {
          if (existing->domain == "acb" &&
              existing->precision_bits != session->precision_bits)
            throw std::invalid_argument(
                "simultaneous Acb sessions must use one working precision");
        }
      }
      session->handle = "s:" + std::to_string(registry.next_session++);
      registry.sessions.emplace(session->handle, session);
    }
    return json::object{{"status", "ok"}, {"session", session->handle},
                        {"domain", session->domain},
                        {"precision_bits", session->precision_bits},
                        {"chart_capacity", session->chart_capacity}};
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
    std::size_t charts = 0;
    {
      std::lock_guard<std::mutex> lock(removed->mutex);
      charts = removed->charts.size();
      removed->charts.clear();
      removed->handles_by_key.clear();
    }
    return json::object{{"status", "ok"}, {"closed", handle},
                        {"released_charts", charts}};
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
          session, root, chart_handle, key,
          std::move(scc), std::move(signature));
    else if (session->domain == "acb")
      chart = parse_prepared_chart<ComplexBall>(
          session, root, chart_handle, key,
          std::move(scc), std::move(signature));
    else
      chart = parse_prepared_chart<SymbolicRational>(
          session, root, chart_handle, key,
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
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      for (const auto& [ignored, chart] : session->charts)
        charts.push_back(chart);
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
          {"scc_components", chart->scc().component_count},
          {"scc_structural_edges", chart->scc().structural_edges.size()},
          {"scc_condensation_edges", chart->scc().condensation_edges.size()},
          {"scc_topological_order",
           encode_indices(chart->scc().topological_order)},
          {"scc_coupling_depth", chart->scc().coupling_depth},
          {"prepare_parse_ms", stats.prepare_parse_ms},
          {"run_parse_ms", stats.run_parse_ms},
          {"kernel_ms", stats.kernel_ms}});
    }
    return json::object{{"status", "ok"}, {"session", session->handle},
                        {"domain", session->domain},
                        {"precision_bits", session->precision_bits},
                        {"charts", charts.size()}, {"runs", runs},
                        {"static_tensor_copies", 0},
                        {"prepare_parse_ms", prepare_parse_ms},
                        {"run_parse_ms", run_parse_ms},
                        {"kernel_ms", kernel_ms},
                        {"chart_stats", std::move(chart_stats)}};
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
    std::lock_guard<std::mutex> lock(symbolic_run_mutex());
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
      std::vector<std::thread> workers;
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
                                      {"backend", "DiffExp2 C++"},
                                      {"flint", flint_version},
                                      {"librarylink", true}});
}

}  // namespace diffexp2
