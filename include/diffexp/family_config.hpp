#pragma once
#include "diffexp/recursion_pipeline.hpp"
#include <boost/json.hpp>
#include <fstream>

namespace diffexp::family_config {
namespace json=boost::json;
struct Configuration {
  feynman::ExampleFamily family;
  std::vector<ibp::Integral> integrals;
  recursion::Options preparation;
  recursion::NumericalOptions numerical;
  unsigned epsilon_order=0,ibp_dots=1,ibp_numerators=2;
  std::string ibp_provider="auto";
};
inline bool same_geometry(const feynman::ExampleFamily& a,const feynman::ExampleFamily& b) {
  if(a.physical_count!=b.physical_count || a.dimension_at_epsilon_zero!=b.dimension_at_epsilon_zero ||
     a.momenta.loops!=b.momenta.loops || a.momenta.external_gram!=b.momenta.external_gram ||
     a.momenta.lines.size()!=b.momenta.lines.size())return false;
  for(unsigned i=0;i<a.momenta.lines.size();++i) {
    const auto& x=a.momenta.lines[i];const auto& y=b.momenta.lines[i];
    if(x.loop_coefficients!=y.loop_coefficients || x.external_coefficients!=y.external_coefficients || x.mass_squared!=y.mass_squared)return false;
  }
  return true;
}
inline Rational rational(const json::value& v) {
  if(v.is_string())return Rational(std::string(v.as_string()));
  if(v.is_int64())return Rational(std::to_string(v.as_int64()));
  if(v.is_uint64())return Rational(std::to_string(v.as_uint64()));
  throw std::invalid_argument("family coefficients must be exact integers or rational strings, e.g. 3/2");
}
inline long integer(const json::value& v,long low,long high,const char* name) {
  if(v.is_uint64() && v.as_uint64()<=static_cast<std::uint64_t>(high) && (low<=0 || v.as_uint64()>=static_cast<std::uint64_t>(low)))return static_cast<long>(v.as_uint64());
  if(!v.is_int64() || v.as_int64()<low || v.as_int64()>high)
    throw std::invalid_argument(std::string(name)+" outside its supported integer range");
  return v.as_int64();
}
inline std::vector<Rational> rationals(const json::value& v) {
  std::vector<Rational> out;for(const auto& x:v.as_array())out.push_back(rational(x));return out;
}
inline void known_keys(const json::object& o,std::initializer_list<std::string_view> keys,const char* where) {
  for(const auto& kv:o)if(std::find(keys.begin(),keys.end(),std::string_view(kv.key()))==keys.end())
    throw std::invalid_argument(std::string("unknown ")+where+" field: "+std::string(kv.key()));
}
inline Configuration parse(const json::value& value) {
  const auto& o=value.as_object();
  known_keys(o,{"schema","name","loops","external_gram","propagators","physical_count","dimension_at_epsilon_zero",
    "integrals","epsilon_order","preparation","numerical","causal","description"},"family");
  if(auto p=o.if_contains("schema");p && *p!="DiffExp.FeynmanFamily/v1")throw std::invalid_argument("unsupported family schema");
  Configuration c;auto& f=c.family;
  f.name=o.if_contains("name")?std::string(o.at("name").as_string()):"configured-family";
  if(f.name.empty() || f.name.size()>256)throw std::invalid_argument("family name must have 1..256 characters");
  f.momenta.loops=integer(o.at("loops"),1,16,"loop count");
  for(const auto& row:o.at("external_gram").as_array())f.momenta.external_gram.push_back(rationals(row));
  if(f.momenta.external_gram.size()>16)throw std::invalid_argument("external Gram dimension exceeds 16");
  for(const auto& p:o.at("propagators").as_array()) {
    const auto& line=p.as_object();known_keys(line,{"loop_coefficients","external_coefficients","mass_squared"},"propagator");
    f.momenta.lines.push_back({rationals(line.at("loop_coefficients")),rationals(line.at("external_coefficients")),rational(line.at("mass_squared"))});
  }
  if(f.momenta.lines.size()>256)throw std::invalid_argument("propagator count exceeds 256");
  f.physical_count=o.if_contains("physical_count")?integer(o.at("physical_count"),2,f.momenta.lines.size(),"physical count"):f.momenta.lines.size();
  f.dimension_at_epsilon_zero=integer(o.at("dimension_at_epsilon_zero"),1,32,"base dimension");f.momenta.validate();
  if(f.physical_count<2)throw std::invalid_argument("at least two physical propagators are required");
  const auto slots=f.momenta.loops*(f.momenta.loops+1)/2+f.momenta.loops*f.momenta.external_gram.size();
  if(auto p=o.if_contains("integrals"))for(const auto& row:p->as_array()) {
    ibp::Integral indices;for(const auto& v:row.as_array())indices.push_back(integer(v,-100,100,"integral power"));
    if(indices.size()<f.physical_count || indices.size()>slots)throw std::invalid_argument("integral indices must cover physical slots and fit the scalar-product basis");
    indices.resize(slots,0);c.integrals.push_back(std::move(indices));
  }
  if(c.integrals.size()>5000)throw std::invalid_argument("integral request count exceeds 5000");
  if(auto p=o.if_contains("epsilon_order"))c.epsilon_order=integer(*p,0,100,"epsilon order");
  if(auto p=o.if_contains("preparation")) {
    const auto& a=p->as_object();known_keys(a,{"anchors","merges","total_seconds","level_seconds","fire_seconds","max_sources_per_level","ibp_provider","ibp_dots","ibp_numerators"},"preparation");
    if(auto v=a.if_contains("ibp_provider")){c.ibp_provider=std::string(v->as_string());if(c.ibp_provider!="auto"&&c.ibp_provider!="ibp-solver"&&c.ibp_provider!="fire"&&c.ibp_provider!="fire-modular")throw std::invalid_argument("unknown IBP provider");}
    if(auto v=a.if_contains("ibp_dots"))c.ibp_dots=integer(*v,0,8,"IBP dots");
    if(auto v=a.if_contains("ibp_numerators"))c.ibp_numerators=integer(*v,0,8,"IBP numerators");
    if(auto v=a.if_contains("anchors"))c.preparation.anchors=rationals(*v);
    if(auto v=a.if_contains("merges"))for(const auto& pair:v->as_array()) {
      if(pair.as_array().size()!=2)throw std::invalid_argument("each merge has two zero-based positions");
      c.preparation.merges.emplace_back(integer(pair.as_array()[0],0,255,"merge index"),integer(pair.as_array()[1],0,255,"merge index"));
    }
    if(auto v=a.if_contains("total_seconds"))c.preparation.total_timeout_seconds=integer(*v,1,86400,"total seconds");
    if(auto v=a.if_contains("level_seconds"))c.preparation.reduction.total_timeout_seconds=integer(*v,1,86400,"level seconds");
    if(auto v=a.if_contains("fire_seconds"))c.preparation.reduction.provider.timeout_seconds=integer(*v,1,86400,"FIRE seconds");
    if(auto v=a.if_contains("max_sources_per_level"))c.preparation.max_sources_per_level=integer(*v,1,100000,"source budget");
  }
  if(auto p=o.if_contains("numerical")) {
    const auto& a=p->as_object();known_keys(a,{"endpoint_order","ordinary_order","working_bits","leaf_digits","method","transport","transport_digits","contour_height","overlap"},"numerical");
    if(auto v=a.if_contains("endpoint_order"))c.numerical.endpoint_order=integer(*v,1,1000,"endpoint order");
    if(auto v=a.if_contains("ordinary_order"))c.numerical.ordinary_order=integer(*v,8,1000,"ordinary order");
    if(auto v=a.if_contains("working_bits"))c.numerical.working_bits=integer(*v,64,1000000,"working bits");
    if(auto v=a.if_contains("leaf_digits"))c.numerical.leaf_digits=integer(*v,1,100000,"leaf digits");
    if(auto v=a.if_contains("contour_height"))c.numerical.contour_height=rational(*v);
    if(auto v=a.if_contains("overlap"))c.numerical.overlap=rational(*v);
    if(auto v=a.if_contains("transport")) {
      const auto method=v->as_string();if(method!="auto" && method!="spectral" && method!="taylor")throw std::invalid_argument("unknown FT ordinary transport");
      c.numerical.ordinary_method=method=="auto"?recursion::OrdinaryMethod::automatic:method=="spectral"?recursion::OrdinaryMethod::spectral:recursion::OrdinaryMethod::taylor;
    }
    if(auto v=a.if_contains("transport_digits"))c.numerical.spectral.accuracy_goal=integer(*v,1,100000,"FT transport digits");
    if(auto v=a.if_contains("method")) {
      const auto method=v->as_string();
      if(method!="adjoint" && method!="factored" && method!="auto" && method!="values")throw std::invalid_argument("unknown numerical method");
      c.numerical.observable_adjoint=method!="values";
      c.numerical.linear_method=method=="factored"?recursion::LinearMethod::factored:method=="auto"?recursion::LinearMethod::automatic:recursion::LinearMethod::adjoint;
    }
  }
  if(auto p=o.if_contains("causal")) {
    const auto& a=p->as_object();known_keys(a,{"f_rim","level_signs","provenance"},"causal");
    causal::Prescription prescription;prescription.f_rim=integer(a.at("f_rim"),-1,1,"F rim");
    for(const auto& sign:a.at("level_signs").as_array())prescription.levels.push_back({static_cast<int>(integer(sign,-1,1,"contour sign"))});
    prescription.provenance=std::string(a.at("provenance").as_string());prescription.validate(f.physical_count-1);
    c.numerical.causal_prescription=std::move(prescription);
  }
  return c;
}
inline json::value describe(const feynman::ExampleFamily& f) {
  json::array gram,lines;
  for(const auto& row:f.momenta.external_gram) {json::array out;for(const auto& q:row)out.emplace_back(q.str());gram.push_back(std::move(out));}
  for(const auto& line:f.momenta.lines) {
    json::array loops,external;for(const auto& q:line.loop_coefficients)loops.emplace_back(q.str());for(const auto& q:line.external_coefficients)external.emplace_back(q.str());
    lines.push_back(json::object{{"loop_coefficients",std::move(loops)},{"external_coefficients",std::move(external)},{"mass_squared",line.mass_squared.str()}});
  }
  json::object out{{"schema","DiffExp.FeynmanFamily/v1"},{"name",f.name},{"loops",f.momenta.loops},
    {"external_gram",std::move(gram)},{"propagators",std::move(lines)},{"physical_count",f.physical_count},
    {"dimension_at_epsilon_zero",f.dimension_at_epsilon_zero},{"epsilon_order",f.name=="bubble"?4:f.name=="pentagon_massive"?2:0},
    {"numerical",json::object{{"endpoint_order",32},{"ordinary_order",80},{"working_bits",384},{"leaf_digits",28},{"method","adjoint"},{"transport","auto"}}}};
  if(f.name=="henn_double_pentagon_x0") {
    json::array anchors,signs;for(const auto& q:causal::henn_anchors())anchors.emplace_back(q.str());
    for(unsigned i=0;i<7;++i)signs.push_back(i==0?1:-1);
    out["preparation"]=json::object{{"anchors",std::move(anchors)}};
    out["causal"]=json::object{{"f_rim",-1},{"level_signs",std::move(signs)},
      {"provenance","Published original ordered ladder; experimental full reconstruction is frozen."}};
  }
  return out;
}
inline Configuration load(const std::string& path) {
  if(path!="-" && !std::filesystem::exists(path)) {
    if(path.find('/')!=std::string::npos || path.ends_with(".json"))throw std::runtime_error("cannot open family configuration: "+path);
    // Kept for command-line convenience; the generic parser never dispatches by name.
    return parse(describe(feynman::example_family(path)));
  }
  std::ifstream file;if(path!="-"){file.open(path);if(!file)throw std::runtime_error("cannot open family configuration");}
  auto& input=path=="-"?std::cin:file;std::string text;char block[4096];
  while(input.read(block,sizeof(block))||input.gcount()) {
    text.append(block,static_cast<std::size_t>(input.gcount()));
    if(text.size()>16*1024*1024)throw std::invalid_argument("family configuration exceeds 16 MiB");
  }
  return parse(json::parse(text));
}
} // namespace diffexp::family_config
