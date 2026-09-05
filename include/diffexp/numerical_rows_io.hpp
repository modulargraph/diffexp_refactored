#pragma once
#include "diffexp/adjoint_transport.hpp"
#include "diffexp/artifact_store.hpp"

// Lossless numerical data only. This codec confers no tail certificate.
namespace diffexp::numerical_rows_io {
namespace json=boost::json;using B=Jet::Ball;
inline std::string exact_ball_part(const arb_t part){char* text=arb_dump_str(part);std::string result(text);flint_free(text);return result;}
inline json::array exact_matrix(const ExactEpsilonMatrix& matrix) {
  json::array result;for(const auto& row:matrix){json::array values;for(const auto& value:row)values.emplace_back(value.str());result.push_back(std::move(values));}return result;
}
inline json::object exact_rows(const LaurentRows& rows) {
  (void)rows.columns();json::array result;
  for(const auto& row:rows.coefficients){json::array columns;for(const auto& series:row){json::array coefficients;
    for(const auto& value:series)coefficients.push_back(json::array{exact_ball_part(acb_realref(value.raw())),exact_ball_part(acb_imagref(value.raw()))});
    columns.push_back(std::move(coefficients));}result.push_back(std::move(columns));}
  return {{"low",rows.low},{"high",rows.high},{"coefficients",result}};
}
inline LaurentRows read_rows(const json::value& value) {
  const auto& o=value.as_object();artifacts::detail::keys(o,{"low","high","coefficients"});
  auto low=artifacts::detail::integer(o.at("low")),high=artifacts::detail::integer(o.at("high"));
  const auto& rows=o.at("coefficients").as_array();
  if(low< -1000 || high>1000 || low>high || high-low>1000 || rows.empty() || rows.size()>5000)
    throw std::invalid_argument("checkpoint Laurent row budget");
  const auto columns=rows.front().as_array().size(),width=static_cast<std::size_t>(high-low+1);
  if(!columns || columns>5000 || rows.size()>20000000/columns/width)throw std::invalid_argument("checkpoint coefficient storage budget");
  LaurentRows result{static_cast<int>(low),static_cast<int>(high),{}};
  for(const auto& row:rows) {
    if(row.as_array().size()!=columns)throw std::invalid_argument("checkpoint row shape mismatch");
    result.coefficients.emplace_back();
    for(const auto& series:row.as_array()) {
      if(series.as_array().size()!=width)throw std::invalid_argument("checkpoint epsilon window mismatch");
      result.coefficients.back().emplace_back();
      for(const auto& encoded:series.as_array()) {
        const auto& pair=encoded.as_array();if(pair.size()!=2)throw std::invalid_argument("checkpoint complex ball shape");
        B value;
        for(unsigned component=0;component<2;++component) {
          const auto text=artifacts::detail::string(pair[component]);
          if(text.size()>512000)throw std::length_error("checkpoint ball text budget");
          auto part=component?acb_imagref(value.raw()):acb_realref(value.raw());
          if(arb_load_str(part,text.c_str()) || !arb_is_finite(part) || exact_ball_part(part)!=text)
            throw std::invalid_argument("checkpoint ball is malformed, nonfinite or noncanonical");
        }
        result.coefficients.back().back().push_back(std::move(value));
      }
    }
  }
  result.columns();return result;
}
} // namespace diffexp::numerical_rows_io
