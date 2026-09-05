#pragma once

#include "diffexp/kernel/recurrence.hpp"

#include <boost/json.hpp>

#include <cstdint>

namespace diffexp::kernel::json_codec_detail {

std::int32_t parse_validity(const boost::json::value& value);

template <typename Scalar>
Scalar parse_scalar(const boost::json::value& value);

template <>
Rational parse_scalar<Rational>(const boost::json::value& value);
template <>
ComplexBall parse_scalar<ComplexBall>(const boost::json::value& value);
template <>
SymbolicRational parse_scalar<SymbolicRational>(
    const boost::json::value& value);

template <typename Scalar>
PreparedRecurrenceOperator<Scalar> parse_prepared_operator(
    const boost::json::object& root);

template <typename Scalar>
void parse_run_state(const boost::json::object& run,
                     const PreparedRecurrenceOperator<Scalar>& prepared,
                     RecurrenceProblem<Scalar>& problem);

boost::json::value encode_validity(std::int32_t value);
boost::json::value encode_scalar(const Rational& value, int digits);
boost::json::value encode_scalar(const SymbolicRational& value, int digits);
boost::json::value encode_scalar(const ComplexBall& value, int digits);
template <typename Scalar>
boost::json::object run_prepared_problem(
    const PreparedRecurrenceOperator<Scalar>& prepared,
    RecurrenceProblem<Scalar>& problem, int digits);

template <typename Scalar>
boost::json::object run_typed(const boost::json::object& root, int digits);

}  // namespace diffexp::kernel::json_codec_detail
