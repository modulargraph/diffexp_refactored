#pragma once
#include "diffexp/artifact_store.hpp"
#include "diffexp/level_preparation.hpp"

namespace diffexp::cached_level {
namespace json=boost::json;
struct Conventions {
  json::object normalization{{"measure","native scalar integral"},{"propagators","exact affine polynomials as supplied"}};
  json::object branch{{"status","formal rational IBP and differential identities"}};
  json::object boundary{{"status","not-applicable to exact level closure"}};
};
struct Result {
  level::Result result;
  bool cache_hit=false;
  std::string semantic_id,content_id;
};
namespace detail {
inline json::array integral_json(const ibp::Integral& integral) {
  json::array result;for(int n:integral)result.push_back(n);return result;
}
inline json::array integrals_json(const std::vector<ibp::Integral>& integrals) {
  json::array result;for(const auto& a:integrals)result.push_back(integral_json(a));return result;
}
inline void validate_integral(const ibp::Integral& integral,const ibp::PropagatorBasis& basis) {
  if(integral.size()!=basis.denominators.size())throw std::invalid_argument("cached level integral arity mismatch");
  for(std::size_t i=0;i<integral.size();++i)
    if(integral[i] < -1000 || integral[i]>1000 || (i>=basis.physical_count && integral[i]>0))
      throw std::invalid_argument("cached level integral index violates denominator/numerator contract");
}
inline ibp::Integral integral_from(const json::value& value,const ibp::PropagatorBasis& basis) {
  ibp::Integral result;for(const auto& n:value.as_array())result.push_back(artifacts::detail::integer(n));
  validate_integral(result,basis);return result;
}
inline json::array exact_row(const std::vector<Exact>& row) {
  json::array result;for(const auto& c:row)result.emplace_back(c.str());return result;
}
inline json::array matrix_json(const std::vector<std::vector<Exact>>& matrix) {
  json::array result;for(const auto& row:matrix)result.push_back(exact_row(row));return result;
}
inline std::vector<std::vector<Exact>> matrix_from(const json::value& value,const ExactField& field,
    std::size_t rows,std::size_t columns) {
  if(value.as_array().size()!=rows)throw std::invalid_argument("cached level matrix row count mismatch");
  std::vector<std::vector<Exact>> result;
  for(const auto& row:value.as_array()) {
    if(row.as_array().size()!=columns)throw std::invalid_argument("cached level matrix column count mismatch");
    result.emplace_back();for(const auto& c:row.as_array())result.back().emplace_back(field,artifacts::detail::string(c));
  }
  return result;
}
inline json::array relations_json(const std::vector<ibp::Relation>& relations) {
  json::array result;
  for(const auto& row:relations) {
    json::array terms;for(const auto& [integral,c]:row)
      terms.push_back(json::object{{"integral",integral_json(integral)},{"coefficient",c.str()}});
    result.push_back(std::move(terms));
  }
  return result;
}
inline json::array witnesses_json(const std::vector<ibp::ExactReducer::Witness>& witnesses) {
  json::array result;for(const auto& witness:witnesses) {
    json::array terms;for(const auto& [id,c]:witness)
      terms.push_back(json::object{{"source",static_cast<std::int64_t>(id)},{"coefficient",c.str()}});
    result.push_back(std::move(terms));
  }
  return result;
}
inline json::object payload(const level::Result& result) {
  return {{"schema","DiffExp3.ExactLevelClosure/v1"},{"ordered_basis",integrals_json(result.ordered_basis)},
    {"requested_integrals",integrals_json(result.requested_integrals)},
    {"matrix",matrix_json(result.matrix)},{"target_rows",matrix_json(result.target_rows)},
    {"source_identities",relations_json(result.source_identities)},
    {"differential_witnesses",witnesses_json(result.differential_witnesses)},
    {"target_witnesses",witnesses_json(result.target_witnesses)},
    {"passes",static_cast<std::int64_t>(result.passes)},{"demands",static_cast<std::int64_t>(result.demands)}};
}
inline level::Result decode(const json::value& payload,const ibp::PropagatorBasis& basis,const ExactField& field) {
  const auto& o=payload.as_object();artifacts::detail::keys(o,{"schema","ordered_basis","requested_integrals","matrix","target_rows","source_identities","differential_witnesses","target_witnesses","passes","demands"});
  if(artifacts::detail::string(o.at("schema"))!="DiffExp3.ExactLevelClosure/v1")throw std::invalid_argument("unsupported exact level closure schema");
  level::Result result;
  const auto& masters=o.at("ordered_basis").as_array();const auto& requested=o.at("requested_integrals").as_array();
  const auto& sources=o.at("source_identities").as_array();
  if(masters.size()>1024 || requested.size()>10000 || sources.size()>200000)
    throw std::invalid_argument("cached level closure exceeds verification size limit");
  for(const auto& a:masters)result.ordered_basis.push_back(integral_from(a,basis));
  for(const auto& a:requested)result.requested_integrals.push_back(integral_from(a,basis));
  result.matrix=matrix_from(o.at("matrix"),field,masters.size(),masters.size());
  result.target_rows=matrix_from(o.at("target_rows"),field,requested.size(),masters.size());
  for(const auto& source:sources) {
    ibp::Relation row;
    if(source.as_array().empty() || source.as_array().size()>100000)throw std::invalid_argument("invalid cached source identity size");
    for(const auto& item:source.as_array()) {
      const auto& term=item.as_object();artifacts::detail::keys(term,{"integral","coefficient"});
      auto integral=integral_from(term.at("integral"),basis);Exact c(field,artifacts::detail::string(term.at("coefficient")));
      if(c.is_zero() || !row.emplace(std::move(integral),std::move(c)).second)
        throw std::invalid_argument("cached source identity contains duplicate or zero terms");
    }
    result.source_identities.push_back(std::move(row));
  }
  auto witnesses=[&](const json::value& value,std::size_t rows) {
    if(value.as_array().size()!=rows)throw std::invalid_argument("cached witness row count mismatch");
    std::vector<ibp::ExactReducer::Witness> result;
    for(const auto& row:value.as_array()) {
      result.emplace_back();for(const auto& item:row.as_array()) {
        const auto& term=item.as_object();artifacts::detail::keys(term,{"source","coefficient"});
        const auto id=artifacts::detail::integer(term.at("source"));Exact c(field,artifacts::detail::string(term.at("coefficient")));
        if(id<0 || static_cast<std::size_t>(id)>=sources.size() || c.is_zero() ||
            !result.back().emplace(id,std::move(c)).second)throw std::invalid_argument("invalid cached witness source or coefficient");
      }
    }
    return result;
  };
  result.differential_witnesses=witnesses(o.at("differential_witnesses"),masters.size());
  result.target_witnesses=witnesses(o.at("target_witnesses"),requested.size());
  auto passes=artifacts::detail::integer(o.at("passes")),demands=artifacts::detail::integer(o.at("demands"));
  if(passes<0 || demands<0)throw std::invalid_argument("invalid cached preparation counters");
  result.passes=passes;result.demands=demands;result.equations=sources.size();result.success=true;return result;
}
inline void verify(const level::Result& result,const ibp::PropagatorBasis& basis,const Exact& dimension,
    std::size_t parameter,const std::vector<ibp::Integral>& requested) {
  if(!result.success || result.requested_integrals!=requested || result.matrix.size()!=result.ordered_basis.size() ||
      result.target_rows.size()!=requested.size() || result.differential_witnesses.size()!=result.matrix.size() ||
      result.target_witnesses.size()!=requested.size())throw std::runtime_error("cached closure shape/request mismatch");
  std::set<ibp::Integral> unique;
  for(const auto& a:result.ordered_basis) {
    validate_integral(a,basis);if(!unique.insert(a).second)throw std::runtime_error("duplicate cached master");
  }
  ibp::Generator generator(basis,dimension);
  ibp::ExactReducer reducer(dimension,std::max<std::size_t>(1,result.source_identities.size()));
  for(const auto& row:result.source_identities)reducer.insert(row);
  // Independently reconstruct the selected span, not just the claimed rows.
  if(!result.ordered_basis.empty()) {
    ibp::BasisReduction coordinates(reducer,result.ordered_basis,dimension);
    (void)coordinates;
  }
  const auto one=dimension.constant(1);
  auto check=[&](ibp::Relation input,const std::vector<Exact>& coordinates,const ibp::ExactReducer::Witness& witness) {
    if(coordinates.size()!=result.ordered_basis.size())throw std::runtime_error("cached coordinate width mismatch");
    for(std::size_t j=0;j<coordinates.size();++j)ibp::add(input,result.ordered_basis[j],-coordinates[j]);
    for(const auto& [id,c]:witness) {
      if(id>=result.source_identities.size())throw std::runtime_error("cached witness source outside payload");
      ibp::add_scaled(input,result.source_identities[id],-c);
    }
    if(!input.empty())throw std::runtime_error("cached exact closure witness reconstruction failed");
  };
  for(std::size_t i=0;i<result.ordered_basis.size();++i)
    check(generator.derivative(result.ordered_basis[i],parameter),result.matrix[i],result.differential_witnesses[i]);
  for(std::size_t i=0;i<requested.size();++i)check({{requested[i],one}},result.target_rows[i],result.target_witnesses[i]);
}
} // namespace detail

inline artifacts::Identity identity(const ibp::PropagatorBasis& basis,const Exact& dimension,
    const ExactField& field,std::size_t parameter,const std::vector<ibp::Integral>& requested,
    const level::Options& options,const Conventions& conventions={}) {
  if(parameter>=field.variables().size() || requested.empty() || requested.size()>10000 ||
      basis.denominators.empty() || basis.denominators.size()>256 || !basis.physical_count ||
      basis.physical_count>basis.denominators.size())throw std::invalid_argument("invalid cached level input shape");
  auto exact=[&](const Exact& c) {
    if(c.variables()!=field.variables())throw std::invalid_argument("cached level exact field context mismatch");return c.str();
  };
  (void)exact(dimension);for(const auto& a:requested)detail::validate_integral(a,basis);
  for(const auto& sector:options.provider.zero_sectors)detail::validate_integral(sector,basis);
  json::array gram,pairs,denominators,rewrites,symbols;
  for(const auto& row:basis.space.external_gram) {
    json::array r;for(const auto& c:row)r.emplace_back(exact(c));gram.push_back(std::move(r));
  }
  for(const auto& [a,b]:basis.space.pairs)pairs.push_back(json::array{a,b});
  auto affine=[&](const ibp::Affine& d) {
    if(d.linear.size()!=basis.space.size())throw std::invalid_argument("cached level affine shape mismatch");
    json::array linear;for(const auto& c:d.linear)linear.emplace_back(exact(c));
    return json::object{{"constant",exact(d.constant)},{"linear",linear}};
  };
  for(const auto& d:basis.denominators)denominators.push_back(affine(d));
  for(const auto& d:basis.scalar_products_in_denominators)rewrites.push_back(affine(d));
  for(const auto& name:field.variables())symbols.emplace_back(name);
  // The public basis owns a derived inverse map. Bind it as well, so mutable
  // caller structures cannot reuse a cache entry after changing that map.
  artifacts::Identity id;id.kind="exact_equation";id.algorithm_version="native-exact-level-closure-v1";
  id.family={{"loops",basis.space.loops},{"external_gram",gram},{"scalar_product_pairs",pairs},
    {"physical_count",static_cast<std::int64_t>(basis.physical_count)},{"denominators",denominators},
    {"scalar_products_in_denominators",rewrites}};
  id.ordered_basis={json::object{{"input_denominators",denominators}},json::object{{"requested_integrals",detail::integrals_json(requested)}}};
  id.normalization=conventions.normalization;id.branch=conventions.branch;id.boundary=conventions.boundary;
  id.geometry={{"parameter",field.variables()[parameter]},{"scalar_product_pairs",pairs},{"external_gram",gram}};
  id.scientific_inputs={{"dimension",exact(dimension)},{"ordered_field_symbols",symbols},
    {"proven_zero_sectors",detail::integrals_json(options.provider.zero_sectors)}};
  id.json_value();return id;
}

inline Result prepare(artifacts::Store& store,const ibp::PropagatorBasis& basis,const Exact& dimension,
    const ExactField& field,std::size_t parameter,std::vector<ibp::Integral> requested,
    const level::Options& options,level::Provider provider={},const Conventions& conventions={}) {
  Result output;auto key=identity(basis,dimension,field,parameter,requested,options,conventions);output.semantic_id=key.key();
  // Exact rational closure has no epsilon/Taylor/precision truncation. The
  // fixed store resource record is intentional; runtime refinement is not an
  // input to this exact operation, and even a missing FIRE executable can hit.
  const Demand exact_resource{0,0,0,64,0};
  constexpr const char* verifier="native-exact-level-witness-v1";
  if(auto cached=store.lookup(key,exact_resource,{"exact",verifier,"exact_closure"})) {
    output.result=detail::decode(cached->payload,basis,field);
    detail::verify(output.result,basis,dimension,parameter,requested);
    output.cache_hit=true;output.content_id=cached->content_id;return output;
  }
  output.result=level::prepare(basis,dimension,field,parameter,requested,options,std::move(provider));
  if(!output.result.success)return output;
  detail::verify(output.result,basis,dimension,parameter,requested);
  // Producer and cache-hit paths use the same exact residual verification.
  artifacts::Certificate certificate{"exact",verifier,"exact_closure",
    {{"statement","differential and requested rows close on the ordered output basis"},
     {"assumptions","conditional on the imported source IBP identities and supplied zero sectors"},
     {"check","regenerated derivatives, exact source-witness reconstruction and basis independence"}}};
  auto stored=store.put(key,exact_resource,detail::payload(output.result),certificate);
  output.content_id=stored.content_id;return output;
}
} // namespace diffexp::cached_level
