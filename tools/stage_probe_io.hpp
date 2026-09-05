#pragma once
#include "diffexp/recursion_pipeline.hpp"
#include "diffexp/numerical_rows_io.hpp"
#include <bit>
#include <iomanip>
#include <limits>
#include <sstream>
namespace diffexp::stage_probe {
namespace json=boost::json;using B=Jet::Ball;
inline std::string approximate(double value){std::ostringstream out;out<<std::scientific<<std::setprecision(9)<<value;return out.str();}
using numerical_rows_io::exact_ball_part;
using numerical_rows_io::exact_matrix;
using numerical_rows_io::exact_rows;
using numerical_rows_io::read_rows;
inline json::object stage_payload(const recursion::PreparedAdjointStage& stage,const Exact& dimension,const std::vector<ibp::Integral>& ordered_masters) {
  json::array names,lower_path,upper_path,shifts,masters;
  for(const auto& name:dimension.variables())names.emplace_back(name);
  for(const auto& value:stage.lower_path)lower_path.emplace_back(value.str());
  for(const auto& value:stage.upper_path)upper_path.emplace_back(value.str());
  for(auto shift:stage.epsilon_gauge_shifts)shifts.emplace_back(shift);
  for(const auto& integral:ordered_masters){json::array powers;for(auto power:integral)powers.emplace_back(power);masters.push_back(std::move(powers));}
  return {{"schema","DiffExp3.PreparedAdjointStage/v1"},{"field_symbols",names},{"dimension",dimension.str()},
    {"ordered_master_basis",masters},{"connection",exact_matrix(stage.connection)},
    {"lower_forcing",exact_matrix(stage.lower_forcing)},{"upper_forcing",exact_matrix(stage.upper_forcing)},
    {"lower_endpoint",exact_rows(stage.lower_endpoint)},{"upper_endpoint",exact_rows(stage.upper_endpoint)},
    {"lower_path",lower_path},{"upper_path",upper_path},{"epsilon_gauge_shifts",shifts},
    {"ball_encoding","lossless arb_dump_str real/imaginary pairs"},{"certificate","uncertified retained numerical checkpoint"}};
}
inline json::object stage_payload(const recursion::PreparedAdjointStage& stage,const recursion::Graph& graph,std::size_t depth) {
  return stage_payload(stage,graph.dimension,graph.nodes.at(depth).closure.ordered_basis);
}
inline json::object quality(const LaurentRows& rows) {
  const auto columns=rows.columns();NativeTailMagnitude maximum;
  double radius_log2=-std::numeric_limits<double>::infinity();bool finite=true;
  for(const auto& row:rows.coefficients)for(const auto& series:row)for(const auto& value:series) {
    finite&=value.is_finite();maximum=NativeTailMagnitude::maximum(maximum,adjoint_detail::enclosure_quality(value));
    radius_log2=std::max({radius_log2,mag_get_d_log2_approx(arb_radref(acb_realref(value.raw()))),mag_get_d_log2_approx(arb_radref(acb_imagref(value.raw())))});
  }
  return {{"rows",rows.coefficients.size()},{"columns",columns},{"epsilon_low",rows.low},{"epsilon_high",rows.high},
    {"finite",finite},{"maximum_normalized_arithmetic_radius_approx",approximate(maximum.approximate_upper())},
    {"maximum_component_radius_log2_approx",approximate(radius_log2)},{"omitted_tails_included",false}};
}
inline json::object exact_expression(const linear_boundary::Expression& expression) {
  linear_boundary::Options limits;linear_boundary::detail::source(expression.leaf_source,limits);
  if(linear_boundary::detail::validate(expression.transform,limits)!=expression.leaf_source->values.size())
    throw std::invalid_argument("checkpoint expression transform/source dimensions");
  const auto& leaf=*expression.leaf_source;LaurentRows source{leaf.low,leaf.high(),{}};
  for(const auto& row:leaf.values)source.coefficients.push_back({row});
  return {{"schema","DiffExp3.LinearExpressionCheckpoint/v1"},{"transform",exact_rows(expression.transform)},
    {"leaf_source",exact_rows(source)},{"producer_leaf_tail_certified",leaf.taylor_tail_certified},
    {"certificate","uncertified numerical expression checkpoint"}};
}
inline linear_boundary::Expression read_expression(const json::value& payload,const std::string& checksum) {
  if(artifacts::detail::sha256(artifacts::detail::canonical(payload))!=checksum)
    throw std::invalid_argument("linear expression checksum mismatch");
  const auto& o=payload.as_object();
  artifacts::detail::keys(o,{"schema","transform","leaf_source","producer_leaf_tail_certified","certificate"});
  if(o.at("schema")!="DiffExp3.LinearExpressionCheckpoint/v1" ||
      o.at("certificate")!="uncertified numerical expression checkpoint" || !o.at("producer_leaf_tail_certified").is_bool())
    throw std::invalid_argument("linear expression schema/certificate");
  auto transform=read_rows(o.at("transform")),source=read_rows(o.at("leaf_source"));
  if(source.columns()!=1 || transform.columns()!=source.coefficients.size())
    throw std::invalid_argument("checkpoint expression transform/source dimensions");
  // A producer's certification flag is provenance, not a proof we verified.
  LaurentBoundary leaf{source.low,{},false};
  for(auto& row:source.coefficients)leaf.values.push_back(std::move(row.front()));
  linear_boundary::Options limits;linear_boundary::detail::validate(transform,limits);
  auto owned=std::make_shared<const LaurentBoundary>(std::move(leaf));
  linear_boundary::detail::source(owned,limits);return {std::move(transform),std::move(owned)};
}
inline ExactEpsilonMatrix read_matrix(const json::value& value,const ExactField& field,std::size_t rows,std::size_t columns) {
  if(value.as_array().size()!=rows)throw std::invalid_argument("checkpoint exact matrix row mismatch");
  ExactEpsilonMatrix result;
  for(const auto& row:value.as_array()) {
    if(row.as_array().size()!=columns)throw std::invalid_argument("checkpoint exact matrix column mismatch");
    result.emplace_back();for(const auto& entry:row.as_array())result.back().emplace_back(field,artifacts::detail::string(entry));
  }
  return result;
}
inline recursion::PreparedAdjointStage read_stage(const json::object& report) {
  const auto& payload=report.at("prepared_stage");
  if(artifacts::detail::sha256(artifacts::detail::canonical(payload))!=artifacts::detail::string(report.at("prepared_stage_sha256")))
    throw std::invalid_argument("prepared stage checksum mismatch");
  const auto& o=payload.as_object();
  artifacts::detail::keys(o,{"schema","field_symbols","dimension","ordered_master_basis","connection","lower_forcing","upper_forcing","lower_endpoint","upper_endpoint","lower_path","upper_path","epsilon_gauge_shifts","ball_encoding","certificate"});
  if(o.at("schema")!="DiffExp3.PreparedAdjointStage/v1" || o.at("ball_encoding")!="lossless arb_dump_str real/imaginary pairs" ||
      o.at("certificate")!="uncertified retained numerical checkpoint")throw std::invalid_argument("prepared stage schema or unsupported certificate");
  std::vector<std::string> names;for(const auto& name:o.at("field_symbols").as_array())names.push_back(artifacts::detail::string(name));
  if(names.empty() || names.size()>16)throw std::invalid_argument("checkpoint exact field budget");
  ExactField field(names);Exact dimension(field,artifacts::detail::string(o.at("dimension")));(void)path_epsilon_variables(dimension);
  auto lower=read_rows(o.at("lower_endpoint")),upper=read_rows(o.at("upper_endpoint"));
  const auto d=lower.columns(),r=lower.coefficients.size();
  if(d>256 || upper.columns()!=d || upper.coefficients.size()!=r)throw std::invalid_argument("checkpoint endpoint dimensions");
  auto connection=read_matrix(o.at("connection"),field,d,d),lower_forcing=read_matrix(o.at("lower_forcing"),field,r,d),upper_forcing=read_matrix(o.at("upper_forcing"),field,r,d);
  const auto path=[&](const json::value& value) {
    if(value.as_array().size()<2 || value.as_array().size()>32)throw std::invalid_argument("checkpoint path budget");
    std::vector<Exact> result;for(const auto& vertex:value.as_array())result.emplace_back(field,artifacts::detail::string(vertex));return result;
  };
  std::vector<std::int64_t> shifts;for(const auto& value:o.at("epsilon_gauge_shifts").as_array())shifts.push_back(artifacts::detail::integer(value));
  if(shifts.size()!=d || o.at("ordered_master_basis").as_array().size()!=d)throw std::invalid_argument("checkpoint ordered basis/gauge dimensions");
  return {std::move(connection),std::move(lower_forcing),std::move(upper_forcing),std::move(lower),std::move(upper),path(o.at("lower_path")),path(o.at("upper_path")),std::move(shifts)};
}
inline std::vector<Exact> remaining_path(const std::vector<Exact>& vertices,unsigned leg,double parameter) {
  if(leg+1>=vertices.size() || !std::isfinite(parameter) || parameter<=0 || parameter>1)
    throw std::invalid_argument("checkpoint chart position");
  if(parameter==1)return {vertices.begin()+leg+1,vertices.end()};
  const auto bits=std::bit_cast<std::uint64_t>(parameter);
  const auto exponent_bits=(bits>>52)&2047,mantissa_bits=bits&((std::uint64_t(1)<<52)-1);
  const auto mantissa=exponent_bits?mantissa_bits|(std::uint64_t(1)<<52):mantissa_bits;
  const int exponent=exponent_bits?static_cast<int>(exponent_bits)-1023-52:-1074;
  const auto& from=vertices[leg];
  const auto scale=from.constant(2).pow(static_cast<unsigned>(std::abs(exponent)));
  const auto numerator=from.constant(Rational(std::to_string(mantissa)));
  const auto exact_parameter=exponent<0?numerator/scale:numerator*scale;
  std::vector<Exact> result{from+(vertices[leg+1]-from)*exact_parameter};
  result.insert(result.end(),vertices.begin()+leg+1,vertices.end());return result;
}
} // namespace diffexp::stage_probe
