#pragma once
#include "diffexp/system.hpp"
#include <boost/json.hpp>
#include <fstream>
#include <sstream>
#include <iostream>

namespace diffexp {
inline boost::json::value series_request(const boost::json::value& input) {
  const auto& request=input.as_object();
  ExactField field({"x","eps"});
  RationalSystem::Matrix matrix;
  for(const auto& row:request.at("matrix").as_array()) {
    std::vector<Exact> values;
    for(const auto& entry:row.as_array())values.emplace_back(field,std::string(entry.as_string()));
    matrix.push_back(std::move(values));
  }
  const auto bounded=[&](const char* name,long minimum,long maximum) {
    const auto v=request.at(name).as_int64();
    if(v<minimum || v>maximum)throw std::invalid_argument(std::string(name)+" outside supported resource range");
    return static_cast<int>(v);
  };
  const int low=bounded("epsilon_low",-10000,10000),high=bounded("epsilon_high",-10000,10000);
  const unsigned order=bounded("taylor_order",0,1000000);
  if(high<low || (std::uint64_t(order)+1)*matrix.size()*std::uint64_t(high-low+1)>10000000)
    throw std::invalid_argument("series request exceeds coefficient budget");
  std::vector<std::vector<Rational>> boundary;
  for(const auto& row:request.at("boundary").as_array()) {
    std::vector<Rational> values;
    for(const auto& value:row.as_array())values.push_back(Exact(field,std::string(value.as_string())).rational());
    boundary.push_back(std::move(values));
  }
  const auto center=Exact(field,std::string(request.at("center").as_string())).rational();
  auto local=regular_series(RationalSystem(std::move(matrix),0,1),center,order,low,high,boundary);
  boost::json::array coefficients;
  for(unsigned n=0;n<=order;++n) {
    boost::json::array rows;
    for(unsigned i=0;i<local.dimension;++i) {
      boost::json::array powers;
      for(int e=low;e<=high;++e)powers.emplace_back(local.at(n,i,e).str());
      rows.push_back(std::move(powers));
    }
    coefficients.push_back(std::move(rows));
  }
  return boost::json::object{{"schema","DiffExp.TaylorSeries/v1"},{"center",center.str()},
    {"epsilon_low",low},{"epsilon_high",high},{"taylor_order",order},
    {"dimension",local.dimension},{"coefficient_order","taylor,component,epsilon"},
    {"coefficients",std::move(coefficients)},{"infinite_tail_certified",false}};
}
// Decimal interval strings retain arbitrary precision without JSON doubles.
inline boost::json::object ball_json(const kernel::ComplexBall& value, long digits) {
  auto part=[&](const arb_t x) {
    char* raw=arb_get_str(x,digits,0);if(!raw)throw std::bad_alloc();
    std::string result(raw);flint_free(raw);return result;
  };
  return {{"real",part(acb_realref(value.raw()))},{"imaginary",part(acb_imagref(value.raw()))}};
}
inline int run_series_file(const std::string& path) {
  std::ifstream file;
  if(path!="-") {file.open(path);if(!file)throw std::runtime_error("cannot open request: "+path);}
  auto& input=path=="-"?std::cin:file;
  constexpr std::size_t max_bytes=16*1024*1024;
  std::string text;char block[4096];
  while(input.read(block,sizeof(block)) || input.gcount()) {
    text.append(block,static_cast<std::size_t>(input.gcount()));
    if(text.size()>max_bytes)throw std::invalid_argument("series request exceeds 16 MiB input budget");
  }
  if(input.bad())throw std::runtime_error("cannot read series request");
  std::cout<<boost::json::serialize(series_request(boost::json::parse(text)))<<'\n';return 0;
}
} // namespace diffexp
