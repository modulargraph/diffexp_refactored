#include "json_codec_recurrence.hpp"

#include "json_codec_support.hpp"

#include <chrono>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace diffexp2::json_codec_detail {
namespace json = boost::json;

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
    throw std::invalid_argument(
        "exact rational scalar must have zero imaginary part");
  return Rational(std::string(pair[0].as_string()));
}

template <>
ComplexBall parse_scalar<ComplexBall>(const json::value& value) {
  if (value.is_string())
    return ComplexBall::from_strings(std::string(value.as_string()));
  if (value.is_int64())
    return ComplexBall(static_cast<long>(value.as_int64()));
  const auto& pair = as_array(value, "Acb scalar");
  if (pair.size() != 2 || !pair[0].is_string() || !pair[1].is_string())
    throw std::invalid_argument(
        "Acb scalar must be [real-string,imag-string]");
  return ComplexBall::from_strings(std::string(pair[0].as_string()),
                                   std::string(pair[1].as_string()));
}

template <>
SymbolicRational parse_scalar<SymbolicRational>(const json::value& value) {
  if (value.is_string())
    return SymbolicRational(std::string(value.as_string()));
  if (value.is_int64())
    return SymbolicRational(static_cast<long>(value.as_int64()));
  throw std::invalid_argument(
      "symbolic scalar must be a rational-function string");
}

namespace {

template <typename Scalar>
MatrixShift<Scalar> parse_matrix_shift(const json::value& value,
                                       std::uint32_t dimension) {
  const auto& object = as_object(value, "matrix shift");
  MatrixShift<Scalar> out;
  out.shift = as_i32(object.at("s"), "matrix shift exponent");
  for (const auto& raw_entry : as_array(object.at("e"), "matrix entries")) {
    const auto& entry = as_array(raw_entry, "matrix entry");
    if (entry.size() != 3)
      throw std::invalid_argument("matrix entry must be [row,col,value]");
    MatrixEntry<Scalar> item;
    item.row = as_u32(entry[0], "matrix row");
    item.col = as_u32(entry[1], "matrix column");
    if (item.row >= dimension || item.col >= dimension)
      throw std::invalid_argument(
          "matrix entry outside recurrence dimension");
    item.value = parse_scalar<Scalar>(entry[2]);
    if (!ScalarTraits<Scalar>::is_zero(item.value))
      out.entries.push_back(std::move(item));
  }
  return out;
}

template <typename Scalar>
RecurrenceProblem<Scalar> parse_problem(const json::object& root) {
  RecurrenceProblem<Scalar> problem;
  problem.dimension = as_u32(root.at("d"), "dimension");
  problem.nmax = as_u32(root.at("nmax"), "nmax");
  problem.log_max = as_u32(root.at("p"), "log maximum");
  problem.frame_base = as_i32(root.at("fb"), "frame base");
  problem.frame_width = as_u32(root.at("w"), "frame width");
  problem.has_initial = root.if_contains("has_initial") == nullptr ||
                        root.at("has_initial").as_bool();
  problem.adaptive_lower_frame_probe =
      root.if_contains("adaptive_probe") != nullptr &&
      root.at("adaptive_probe").as_bool();
  problem.a_target = parse_scalar<Scalar>(root.at("a_target"));
  problem.b_target = parse_scalar<Scalar>(root.at("b_target"));
  problem.a_shift_min = as_i32(root.at("a_shift_min"), "a shift minimum");
  for (const auto& value : as_array(root.at("a_shifts"), "a shifts"))
    problem.a_shifts.push_back(parse_scalar<Scalar>(value));

  for (const auto& raw_lag : as_array(root.at("d_lags"), "d lags")) {
    std::vector<ScalarShift<Scalar>> lag;
    for (const auto& raw_shift : as_array(raw_lag, "d lag")) {
      const auto& shift = as_object(raw_shift, "d scalar shift");
      lag.push_back({as_i32(shift.at("s"), "d shift"),
                     parse_scalar<Scalar>(shift.at("v"))});
    }
    problem.d_lags.push_back(std::move(lag));
  }

  for (const auto& raw_denominator :
       as_array(root.at("denominators"), "denominators")) {
    std::vector<Scalar> denominator;
    for (const auto& value : as_array(raw_denominator, "denominator"))
      denominator.push_back(parse_scalar<Scalar>(value));
    problem.rational_denominators.push_back(std::move(denominator));
  }

  for (const auto& raw_lag : as_array(root.at("nhat_lags"), "Nhat lags")) {
    const auto& lag_object = as_object(raw_lag, "Nhat lag");
    PreparedLag<Scalar> lag;
    for (const auto& raw_matrix :
         as_array(lag_object.at("poly"), "polynomial matrices"))
      lag.polynomial.push_back(
          parse_matrix_shift<Scalar>(raw_matrix, problem.dimension));
    for (const auto& raw_group :
         as_array(lag_object.at("rat"), "rational groups")) {
      const auto& group_object = as_object(raw_group, "rational group");
      RationalGroup<Scalar> group;
      group.denominator_index =
          as_u32(group_object.at("q"), "denominator index");
      for (const auto& raw_matrix :
           as_array(group_object.at("num"), "rational numerator"))
        group.numerator.push_back(
            parse_matrix_shift<Scalar>(raw_matrix, problem.dimension));
      lag.rational.push_back(std::move(group));
    }
    for (const auto& value :
         as_array(lag_object.at("val"), "Nhat valuations"))
      lag.valuations.push_back(parse_validity(value));
    problem.nhat_lags.push_back(std::move(lag));
  }

  if (!root.at("d0_inverse").is_null())
    problem.d0_inverse_scalar =
        parse_scalar<Scalar>(root.at("d0_inverse"));

  for (const auto& raw_block : as_array(root.at("blocks"), "Jordan blocks")) {
    JordanBlock block;
    for (const auto& column : as_array(raw_block, "Jordan block"))
      block.columns.push_back(as_u32(column, "Jordan column"));
    problem.blocks.push_back(std::move(block));
  }

  for (const auto& raw_row : as_array(root.at("schedule"), "step schedule")) {
    std::vector<BlockStep<Scalar>> row;
    for (const auto& raw_step : as_array(raw_row, "schedule row")) {
      const auto& step = as_object(raw_step, "schedule step");
      const auto kind = std::string(step.at("case").as_string());
      StepCase step_case;
      if (kind == "T")
        step_case = StepCase::Taylor;
      else if (kind == "P")
        step_case = StepCase::Pseudo;
      else if (kind == "R")
        step_case = StepCase::Resonant;
      else
        throw std::invalid_argument("unknown recurrence step case: " + kind);
      row.push_back({step_case, parse_scalar<Scalar>(step.at("da")),
                     parse_scalar<Scalar>(step.at("db"))});
    }
    problem.schedule.push_back(std::move(row));
  }

  for (const auto& value : as_array(root.at("initial"), "initial tensor"))
    problem.initial.push_back(parse_scalar<Scalar>(value));
  for (const auto& value :
       as_array(root.at("initial_validity"), "initial validity"))
    problem.initial_validity.push_back(parse_validity(value));

  if (const auto* raw_source = root.if_contains("source");
      raw_source != nullptr && !raw_source->is_null()) {
    const auto& source_object = as_object(*raw_source, "source");
    SourceData<Scalar> source;
    for (const auto& value :
         as_array(source_object.at("frames"), "source frames"))
      source.frames.push_back(parse_scalar<Scalar>(value));
    for (const auto& value :
         as_array(source_object.at("validity"), "source validity"))
      source.validity.push_back(parse_validity(value));
    for (const auto& value :
         as_array(source_object.at("present"), "source presence"))
      source.present.push_back(value.as_bool() ? 1 : 0);
    problem.source = std::move(source);
  }

  if (const auto* raw_assembly = root.if_contains("assembly");
      raw_assembly != nullptr && !raw_assembly->is_null()) {
    const auto& assembly = as_object(*raw_assembly, "assembly matrix");
    PreparedMatrix<Scalar> matrix;
    matrix.identity = assembly.if_contains("identity") != nullptr &&
                      assembly.at("identity").as_bool();
    if (const auto* identity = assembly.if_contains("exact_identity"))
      matrix.exact_identity = std::string(identity->as_string());
    for (const auto& raw_matrix :
         as_array(assembly.at("poly"), "assembly polynomial matrices"))
      matrix.polynomial.push_back(
          parse_matrix_shift<Scalar>(raw_matrix, problem.dimension));
    for (const auto& raw_group :
         as_array(assembly.at("rat"), "assembly rational groups")) {
      const auto& group_object =
          as_object(raw_group, "assembly rational group");
      RationalGroup<Scalar> group;
      group.denominator_index =
          as_u32(group_object.at("q"), "assembly denominator index");
      for (const auto& raw_matrix :
           as_array(group_object.at("num"), "assembly rational numerator"))
        group.numerator.push_back(
            parse_matrix_shift<Scalar>(raw_matrix, problem.dimension));
      matrix.rational.push_back(std::move(group));
    }
    for (const auto& value :
         as_array(assembly.at("val"), "assembly valuations"))
      matrix.valuations.push_back(parse_validity(value));
    problem.assembly_matrix = std::move(matrix);
    problem.chop_digits = as_i32(root.at("chop_digits"), "chop digits");
    problem.return_u = root.if_contains("return_u") != nullptr &&
                       root.at("return_u").as_bool();
  }
  return problem;
}

template <typename Scalar>
json::object encode_result(RecurrenceResult<Scalar>&& result,
                           std::optional<AssembledResult<Scalar>>&& assembled,
                           bool return_u, int digits, double elapsed_ms) {
  json::object output;
  output["status"] = "ok";
  output["elapsed_ms"] = elapsed_ms;
  output["top_valid"] = encode_validity(result.top_valid);
  if (return_u) {
    json::array coefficients;
    coefficients.reserve(result.u.size());
    for (const auto& value : result.u)
      coefficients.push_back(encode_scalar(value, digits));
    output["u"] = std::move(coefficients);
    json::array validity;
    validity.reserve(result.validity.size());
    for (const auto value : result.validity)
      validity.push_back(encode_validity(value));
    output["validity"] = std::move(validity);
  }
  if (assembled.has_value()) {
    json::object encoded;
    encoded["min"] = assembled->min_power;
    encoded["max"] = assembled->complete_max;
    json::array coefficients;
    coefficients.reserve(assembled->coefficients.size());
    for (const auto& value : assembled->coefficients)
      coefficients.push_back(encode_scalar(value, digits));
    encoded["coefficients"] = std::move(coefficients);
    output["assembled"] = std::move(encoded);
  }
  json::array hits;
  for (const auto& hit : result.hits) {
    json::object encoded;
    encoded["n"] = hit.n;
    json::array columns;
    for (const auto column : hit.columns) columns.emplace_back(column);
    encoded["cols"] = std::move(columns);
    encoded["delta_b"] = encode_scalar(hit.delta_b, digits);
    json::array frames;
    for (const auto& row : hit.gamma_frames)
      for (const auto& value : row)
        frames.push_back(encode_scalar(value, digits));
    encoded["frames"] = std::move(frames);
    json::array validity;
    for (const auto value : hit.gamma_validity)
      validity.push_back(encode_validity(value));
    encoded["validity"] = std::move(validity);
    hits.push_back(std::move(encoded));
  }
  output["hits"] = std::move(hits);
  return output;
}

template <typename Scalar>
json::object run_problem(RecurrenceProblem<Scalar>& problem, int digits) {
  const auto started = std::chrono::steady_clock::now();
  auto result = RecurrenceSolver<Scalar>(problem).run();
  std::optional<AssembledResult<Scalar>> assembled;
  if (problem.assembly_matrix.has_value())
    assembled = assemble_recurrence(problem, result);
  const auto ended = std::chrono::steady_clock::now();
  const double elapsed =
      std::chrono::duration<double, std::milli>(ended - started).count();
  return encode_result(std::move(result), std::move(assembled),
                       problem.return_u, digits, elapsed);
}

}  // namespace

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
  for (const auto& raw_denominator :
       as_array(root.at("denominators"), "denominators")) {
    std::vector<Scalar> denominator;
    for (const auto& value : as_array(raw_denominator, "denominator"))
      denominator.push_back(parse_scalar<Scalar>(value));
    prepared.rational_denominators.push_back(std::move(denominator));
  }
  for (const auto& raw_lag : as_array(root.at("nhat_lags"), "Nhat lags")) {
    const auto& lag_object = as_object(raw_lag, "Nhat lag");
    PreparedLag<Scalar> lag;
    for (const auto& raw_matrix :
         as_array(lag_object.at("poly"), "polynomial matrices"))
      lag.polynomial.push_back(
          parse_matrix_shift<Scalar>(raw_matrix, prepared.dimension));
    for (const auto& raw_group :
         as_array(lag_object.at("rat"), "rational groups")) {
      const auto& group_object = as_object(raw_group, "rational group");
      RationalGroup<Scalar> group;
      group.denominator_index =
          as_u32(group_object.at("q"), "denominator index");
      for (const auto& raw_matrix :
           as_array(group_object.at("num"), "rational numerator"))
        group.numerator.push_back(
            parse_matrix_shift<Scalar>(raw_matrix, prepared.dimension));
      lag.rational.push_back(std::move(group));
    }
    for (const auto& value :
         as_array(lag_object.at("val"), "Nhat valuations"))
      lag.valuations.push_back(parse_validity(value));
    prepared.nhat_lags.push_back(std::move(lag));
  }
  if (!root.at("d0_inverse").is_null())
    prepared.d0_inverse_scalar =
        parse_scalar<Scalar>(root.at("d0_inverse"));
  for (const auto& raw_block : as_array(root.at("blocks"), "Jordan blocks")) {
    JordanBlock block;
    for (const auto& column : as_array(raw_block, "Jordan block"))
      block.columns.push_back(as_u32(column, "Jordan column"));
    prepared.blocks.push_back(std::move(block));
  }
  if (const auto* raw_assembly = root.if_contains("assembly");
      raw_assembly != nullptr && !raw_assembly->is_null()) {
    const auto& assembly = as_object(*raw_assembly, "assembly matrix");
    PreparedMatrix<Scalar> matrix;
    matrix.identity = assembly.if_contains("identity") != nullptr &&
                      assembly.at("identity").as_bool();
    if (const auto* identity = assembly.if_contains("exact_identity"))
      matrix.exact_identity = std::string(identity->as_string());
    for (const auto& raw_matrix :
         as_array(assembly.at("poly"), "assembly polynomial matrices"))
      matrix.polynomial.push_back(
          parse_matrix_shift<Scalar>(raw_matrix, prepared.dimension));
    for (const auto& raw_group :
         as_array(assembly.at("rat"), "assembly rational groups")) {
      const auto& group_object =
          as_object(raw_group, "assembly rational group");
      RationalGroup<Scalar> group;
      group.denominator_index =
          as_u32(group_object.at("q"), "assembly denominator index");
      for (const auto& raw_matrix :
           as_array(group_object.at("num"), "assembly rational numerator"))
        group.numerator.push_back(
            parse_matrix_shift<Scalar>(raw_matrix, prepared.dimension));
      matrix.rational.push_back(std::move(group));
    }
    for (const auto& value :
         as_array(assembly.at("val"), "assembly valuations"))
      matrix.valuations.push_back(parse_validity(value));
    prepared.assembly_matrix = std::move(matrix);
  }
  prepared.chop_digits = root.if_contains("chop_digits")
      ? as_i32(root.at("chop_digits"), "chop digits")
      : 0;
  retain_framed_d0_inverse(prepared);
  return prepared;
}

template <typename Scalar>
void parse_run_state(const json::object& run,
                     const PreparedRecurrenceOperator<Scalar>& prepared,
                     RecurrenceProblem<Scalar>& problem) {
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
          : kind == "P"             ? StepCase::Pseudo
          : kind == "R"             ? StepCase::Resonant
                                      : throw std::invalid_argument(
                                            "unknown recurrence step case: " +
                                            kind);
      row.push_back({step_case, parse_scalar<Scalar>(step.at("da")),
                     parse_scalar<Scalar>(step.at("db"))});
    }
    problem.schedule.push_back(std::move(row));
  }
  for (const auto& value : as_array(run.at("initial"), "initial tensor"))
    problem.initial.push_back(parse_scalar<Scalar>(value));
  for (const auto& value :
       as_array(run.at("initial_validity"), "initial validity"))
    problem.initial_validity.push_back(parse_validity(value));

  const auto& raw_source = run.at("source");
  if (!raw_source.is_null()) {
    const auto& source_object = as_object(raw_source, "source");
    SourceData<Scalar> source;
    for (const auto& value :
         as_array(source_object.at("frames"), "source frames"))
      source.frames.push_back(parse_scalar<Scalar>(value));
    for (const auto& value :
         as_array(source_object.at("validity"), "source validity"))
      source.validity.push_back(parse_validity(value));
    for (const auto& value :
         as_array(source_object.at("present"), "source presence"))
      source.present.push_back(value.as_bool() ? 1 : 0);
    problem.source = std::move(source);
  }
  problem.return_u = run.at("return_u").as_bool();
  if (!problem.return_u && !prepared.assembly_matrix.has_value())
    throw std::invalid_argument(
        "persistent run suppresses U without a prepared assembly matrix");
}

json::value encode_validity(std::int32_t value) {
  return value == kCompleteInfinity ? json::value(nullptr)
                                    : json::value(value);
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

template PreparedRecurrenceOperator<Rational>
parse_prepared_operator<Rational>(const json::object&);
template PreparedRecurrenceOperator<ComplexBall>
parse_prepared_operator<ComplexBall>(const json::object&);
template PreparedRecurrenceOperator<SymbolicRational>
parse_prepared_operator<SymbolicRational>(const json::object&);

template void parse_run_state<Rational>(
    const json::object&, const PreparedRecurrenceOperator<Rational>&,
    RecurrenceProblem<Rational>&);
template void parse_run_state<ComplexBall>(
    const json::object&, const PreparedRecurrenceOperator<ComplexBall>&,
    RecurrenceProblem<ComplexBall>&);
template void parse_run_state<SymbolicRational>(
    const json::object&, const PreparedRecurrenceOperator<SymbolicRational>&,
    RecurrenceProblem<SymbolicRational>&);

template json::object run_prepared_problem<Rational>(
    const PreparedRecurrenceOperator<Rational>&, RecurrenceProblem<Rational>&,
    int);
template json::object run_prepared_problem<ComplexBall>(
    const PreparedRecurrenceOperator<ComplexBall>&,
    RecurrenceProblem<ComplexBall>&, int);
template json::object run_prepared_problem<SymbolicRational>(
    const PreparedRecurrenceOperator<SymbolicRational>&,
    RecurrenceProblem<SymbolicRational>&, int);

template json::object run_typed<Rational>(const json::object&, int);
template json::object run_typed<ComplexBall>(const json::object&, int);
template json::object run_typed<SymbolicRational>(const json::object&, int);

}  // namespace diffexp2::json_codec_detail
