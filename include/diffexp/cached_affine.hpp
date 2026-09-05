#pragma once
#include "diffexp/affine_frobenius.hpp"
#include "diffexp/artifact_store.hpp"
#include "diffexp/geometry.hpp"

namespace diffexp::cached_affine {
namespace json=boost::json;
using Series=AffineFrobeniusSeries;
using Matrix=Series::Matrix;
// Verification resources are receiving limits, not scientific identity axes.
struct VerificationLimits {
  std::size_t max_polynomial_terms=1000000, max_shared_terms=4000000, max_total_column_terms=4000000,
              max_term_products=2000000000;
  // Cumulative conservative monomial-product estimate, not elapsed CPU work.
  std::function<void(unsigned,unsigned,std::size_t)> column_progress;
};
struct Result {
  Series series;
  bool cache_hit;
  std::string semantic_id,content_id;
  unsigned columns_reused=0,columns_prepared=0;
};
namespace detail {
inline constexpr const char* algorithm="native-affine-frobenius-canonical-frontier-v1";
inline constexpr const char* verifier="native-affine-polynomial-residual-normalization-v1";
inline constexpr const char* scope="exact_retained_frobenius_coefficients";
inline void require(bool condition,const char* why) {
  if(!condition)throw std::runtime_error(why);
}
inline void input(const Matrix& a,std::size_t xi,std::size_t ei,unsigned n,const Series::Options& options) {
  require(!a.empty() && a.size()<=256 && a.size()<=options.max_dimension &&
      n<=1000 && n<=options.max_x_order && options.max_terms && options.max_epsilon_depth,
      "cached affine input exceeds receiving dimensions/order/term budgets");
  require(!a[0].empty(),"cached affine empty source row");
  const auto& sample=a[0][0];
  require(xi!=ei && xi<sample.variable_count() && ei<sample.variable_count(),"cached affine variable indices");
  for(const auto& row:a) {
    require(row.size()==a.size(),"cached affine source matrix must be square");
    for(const auto& c:row) {
      (void)(sample.constant(0)+c); // Also rejects distinct exact field instances.
      for(const auto& terms:{c.numerator_terms(),c.denominator_terms()})
        for(const auto& term:terms)for(unsigned j=0;j<term.powers.size();++j)
          require(!term.powers[j] || j==xi || j==ei,"cached affine source has an unsubstituted parameter");
    }
  }
}
inline json::array matrix_json(const Matrix& matrix) {
  json::array out;for(const auto& row:matrix){json::array r;for(const auto& c:row)r.emplace_back(c.str());out.push_back(std::move(r));}return out;
}
inline json::object term_json(const Series::Term& t) {
  return {{"row",static_cast<std::int64_t>(t.row)},{"column",static_cast<std::int64_t>(t.column)},
    {"log_degree",static_cast<std::int64_t>(t.log_degree)},{"power",t.power.str()},
    {"slope",t.slope.str()},{"coefficient",t.coefficient.str()}};
}
inline json::object payload(const Series& series,bool include_terms=true) {
  json::array exponents,successors,frontiers,terms;
  for(const auto& e:series.exponents())exponents.push_back(json::array{e.power.str(),e.slope.str()});
  for(auto next:series.jordan_successors())successors.push_back(next);
  for(const auto& f:series.absolute_x_frontiers())frontiers.emplace_back(f.str());
  if(include_terms)for(const auto& t:series.terms().terms)terms.push_back(term_json(t));
  return {{"schema","DiffExp3.ExactAffineFrobenius/v1"},{"dimension",static_cast<std::int64_t>(series.dimension())},{"x_order",static_cast<std::int64_t>(series.x_order())},
    {"frame",matrix_json(series.residue_frame())},{"exponents",exponents},{"successors",successors},{"frontiers",frontiers},
    {"terms",terms},{"omitted_tail_certified",series.terms().omitted_tail_certified},
    {"coherent_x_frontier",series.terms().coherent_x_frontier},
    {"wronskian_prefactor",series.terms().wronskian_prefactor?json::value(series.terms().wronskian_prefactor->str()):json::value(nullptr)}};
}
struct State {
  Matrix frame;
  std::vector<Series::Exponent> exponents;
  std::vector<int> successors;
  std::vector<Rational> frontiers;
  Series::Expansion expansion;
};
inline State decode(const json::value& value,const Exact& sample,unsigned d,unsigned n,const Series::Options& options) {
  const auto& o=value.as_object();artifacts::detail::keys(o,{"schema","dimension","x_order","frame","exponents","successors","frontiers","terms","omitted_tail_certified","coherent_x_frontier","wronskian_prefactor"});
  require(artifacts::detail::string(o.at("schema"))=="DiffExp3.ExactAffineFrobenius/v1" &&
      artifacts::detail::integer(o.at("dimension"))==d && artifacts::detail::integer(o.at("x_order"))==n,
      "cached affine payload dimension/order/schema mismatch");
  require(!o.at("omitted_tail_certified").as_bool() && o.at("coherent_x_frontier").as_bool(),"cached affine invalid frontier or omitted-tail claim");
  require(o.at("frame").as_array().size()==d && o.at("exponents").as_array().size()==d &&
      o.at("successors").as_array().size()==d && o.at("frontiers").as_array().size()==d &&
      o.at("terms").as_array().size()<=options.max_terms,"cached affine payload exceeds shape/term budget");
  std::map<std::string,Exact> variables;for(unsigned j=0;j<sample.variable_count();++j)variables.emplace(sample.variables()[j],sample.variable(j));
  auto exact=[&](const json::value& v){return diffexp::evaluate_exact(data::Reader(artifacts::detail::string(v)).read(),sample,variables);};
  State out{{},{},{},{},{d,d,{},false,std::nullopt,true}};
  for(const auto& row:o.at("frame").as_array()) {
    require(row.as_array().size()==d,"cached affine residue frame width");
    out.frame.emplace_back();for(const auto& c:row.as_array())out.frame.back().push_back(exact(c));
  }
  for(unsigned j=0;j<d;++j) {
    const auto& e=o.at("exponents").as_array()[j].as_array();require(e.size()==2,"cached affine exponent shape");
    out.exponents.push_back({Rational(artifacts::detail::string(e[0])),Rational(artifacts::detail::string(e[1]))});
    auto successor=artifacts::detail::integer(o.at("successors").as_array()[j]);
    require(successor>=-1 && successor<d,"cached affine successor index");out.successors.push_back(successor);
    out.frontiers.emplace_back(artifacts::detail::string(o.at("frontiers").as_array()[j]));
  }
  for(const auto& entry:o.at("terms").as_array()) {
    const auto& t=entry.as_object();artifacts::detail::keys(t,{"row","column","log_degree","power","slope","coefficient"});
    auto row=artifacts::detail::integer(t.at("row")),column=artifacts::detail::integer(t.at("column")),log=artifacts::detail::integer(t.at("log_degree"));
    require(row>=0 && row<d && column>=0 && column<d && log>=0 && static_cast<std::size_t>(log)<=static_cast<std::size_t>(d)*(n+1),"cached affine term index/log budget");
    out.expansion.terms.push_back({static_cast<unsigned>(row),static_cast<unsigned>(column),static_cast<unsigned>(log),
        Rational(artifacts::detail::string(t.at("power"))),Rational(artifacts::detail::string(t.at("slope"))),exact(t.at("coefficient"))});
  }
  if(!o.at("wronskian_prefactor").is_null())out.expansion.wronskian_prefactor=exact(o.at("wronskian_prefactor"));
  return out;
}
struct PolynomialBudget {
  VerificationLimits limits;
  std::size_t products=0;
  std::size_t check(const Exact& p) const {
    require(p.denominator().is_rational(),"cached affine residual expected a polynomial");
    auto size=p.numerator_terms().size();require(size<=limits.max_polynomial_terms,"cached affine polynomial term budget exhausted");return size;
  }
  Exact multiply(const Exact& a,const Exact& b) {
    auto na=check(a),nb=check(b);
    require(!nb || na<=(limits.max_term_products-products)/nb,"cached affine residual work budget exhausted");
    products+=na*nb;auto out=a*b;check(out);return out;
  }
};
inline long integer_offset(const Rational& offset,unsigned n) {
  const auto text=offset.str();require(text.find('/')==std::string::npos && offset.sign()>=0 && offset<=Rational(n),"cached affine power outside retained integer frontier");
  return std::stol(text);
}
// Independent verification uses polynomial residuals, not the preparation
// recurrence. A common epsilon denominator removes rational-expression gcds
// from the large residual products. No residual beyond the frontier is claimed.
inline void verify(const Matrix& a,std::size_t xi,std::size_t ei,unsigned n,const Series::Options& options,
    const State& state,const VerificationLimits& limits={}) {
  input(a,xi,ei,n,options);
  require(limits.max_polynomial_terms && limits.max_shared_terms && limits.max_total_column_terms && limits.max_term_products,"cached affine zero verification budget");
  PolynomialBudget budget{limits};const unsigned d=a.size();const auto z=a[0][0].constant(0),x=z.variable(xi);
  auto connection=fuchsify::detail::zeros(d,d,z),residue=connection;
  for(unsigned i=0;i<d;++i)for(unsigned j=0;j<d;++j){connection[i][j]=x*a[i][j];residue[i][j]=affine_frobenius_detail::taylor(connection[i][j],xi,0)[0];}
  const auto canonical=affine_frobenius_detail::frame(residue,xi,ei);
  require(state.frame==canonical.transform && state.successors==canonical.successor && state.exponents.size()==d && state.frontiers.size()==d,"cached affine canonical Jordan frame mismatch");
  std::vector<unsigned> cutoffs;
  for(unsigned c=0;c<d;++c) {
    require(state.exponents[c].power==canonical.exponents[c].power && state.exponents[c].slope==canonical.exponents[c].slope,"cached affine canonical spectrum mismatch");
    auto minimum=canonical.exponents[c].power;
    for(unsigned j=0;j<d;++j)if((canonical.exponents[j].power-canonical.exponents[c].power).str().find('/')==std::string::npos)minimum=std::min(minimum,canonical.exponents[j].power);
    require(state.frontiers[c]==minimum+Rational(n),"cached affine absolute frontier mismatch");
    cutoffs.push_back(integer_offset(state.frontiers[c]-canonical.exponents[c].power,n));
  }
  require(state.expansion.rows==d && state.expansion.columns==d && state.expansion.coherent_x_frontier && !state.expansion.omitted_tail_certified && state.expansion.terms.size()<=options.max_terms,"cached affine expansion shape/certificate budget");
  using Coordinate=std::tuple<unsigned,unsigned,unsigned,unsigned>; // column,k,row,log
  std::map<Coordinate,const Series::Term*> coordinates;
  std::vector<std::vector<const Series::Term*>> columns(d);
  for(const auto& term:state.expansion.terms) {
    require(term.row<d && term.column<d && term.log_degree<=static_cast<std::size_t>(d)*(n+1) && !term.coefficient.is_zero(),"cached affine invalid nonzero term coordinate");
    require(term.slope==canonical.exponents[term.column].slope,"cached affine term slope mismatch");
    auto k=integer_offset(term.power-canonical.exponents[term.column].power,cutoffs[term.column]);
    for(const auto& terms:{term.coefficient.numerator_terms(),term.coefficient.denominator_terms()})
      for(const auto& t:terms)for(unsigned j=0;j<t.powers.size();++j)require(!t.powers[j] || j==ei,"cached affine coefficient depends on a non-epsilon variable");
    require(coordinates.emplace(Coordinate{term.column,k,term.row,term.log_degree},&term).second,"cached affine duplicate term coordinate");
    columns[term.column].push_back(&term);
  }
  // Polynomial ODE data are shared across all retained columns and log powers.
  std::vector<Exact> q(d,z.constant(1));auto p=connection;
  std::size_t shared_terms=0;
  auto count_shared=[&](const Exact& value){auto count=budget.check(value);require(count<=limits.max_shared_terms-shared_terms,"cached affine shared polynomial storage budget exhausted");shared_terms+=count;};
  for(unsigned i=0;i<d;++i) {
    for(unsigned j=0;j<d;++j){q[i]=q[i].polynomial_lcm(connection[i][j].denominator());budget.check(q[i]);}
    count_shared(q[i]);
    for(unsigned j=0;j<d;++j){p[i][j]=q[i]*connection[i][j];count_shared(p[i][j]);}
  }
  const auto inverse=fuchsify::detail::inverse(canonical.transform);
  std::vector<Exact> powers(n+1,z.constant(1));for(unsigned k=1;k<=n;++k)powers[k]=powers[k-1]*x;
  for(unsigned c=0;c<d;++c) {
    auto common=z.constant(1);unsigned log_high=0;
    // Denominators repeat across physical rows and log powers. Intern their
    // LCM contributions and exact quotients once per column.
    std::map<std::string,Exact> denominators;
    std::vector<std::pair<const Series::Term*,std::string>> indexed;
    for(const auto* t:columns[c]) {
      auto denominator=t->coefficient.denominator();auto name=denominator.str();
      denominators.try_emplace(name,std::move(denominator));indexed.emplace_back(t,std::move(name));
      log_high=std::max(log_high,t->log_degree);
    }
    for(const auto& [name,denominator]:denominators){common=common.polynomial_lcm(denominator);budget.check(common);}
    for(auto& [name,denominator]:denominators)denominator=common/denominator;
    require(d<=limits.max_total_column_terms/(static_cast<std::size_t>(log_high)+1),"cached affine column polynomial slot budget exhausted");
    std::vector<std::vector<Exact>> numerator(d,std::vector<Exact>(log_high+1,z));
    for(const auto& [t,name]:indexed) {
      const auto k=integer_offset(t->power-canonical.exponents[c].power,cutoffs[c]);
      auto regular=budget.multiply(t->coefficient.numerator(),denominators.at(name));
      auto term=budget.multiply(regular,powers[k]);
      numerator[t->row][t->log_degree]=numerator[t->row][t->log_degree]+term;
      budget.check(numerator[t->row][t->log_degree]);
    }
    std::size_t total=0;for(const auto& row:numerator)for(const auto& v:row){auto terms=budget.check(v);require(terms<=limits.max_total_column_terms-total,"cached affine total column polynomial budget exhausted");total+=terms;}
    const auto lambda=z.constant(canonical.exponents[c].power)+z.constant(canonical.exponents[c].slope)*z.variable(ei);
    for(unsigned i=0;i<d;++i)for(unsigned l=0;l<=log_high;++l) {
      auto derivative=budget.multiply(x,numerator[i][l].derivative(xi))+budget.multiply(lambda,numerator[i][l]);
      if(l<log_high)derivative=derivative+numerator[i][l+1]*z.constant(l+1);
      auto residual=budget.multiply(q[i],derivative);
      for(unsigned j=0;j<d;++j)if(!p[i][j].is_zero() && !numerator[j][l].is_zero()) {
        residual=residual-budget.multiply(p[i][j],numerator[j][l]);budget.check(residual);
      }
      for(const auto& term:residual.numerator_terms())require(term.powers[xi]>cutoffs[c],"cached affine defining polynomial residual is nonzero within frontier");
    }
    // A polynomial ODE alone permits arbitrary resonant homogeneous constants.
    // These checks bind precisely the zero-integration-constant convention used
    // by prepare, including each seeded diagonal constant at k=0.
    for(unsigned i=0;i<d;++i)if(canonical.exponents[i].slope==canonical.exponents[c].slope) {
      const auto offset=canonical.exponents[i].power-canonical.exponents[c].power;
      if(offset.sign()<0 || offset>Rational(cutoffs[c]) || offset.str().find('/')!=std::string::npos)continue;
      const auto k=integer_offset(offset,cutoffs[c]);auto value=z;
      for(unsigned j=0;j<d;++j)if(auto found=coordinates.find(Coordinate{c,k,j,0});found!=coordinates.end())value=value+inverse[i][j]*found->second->coefficient;
      require(value==z.constant(k==0 && i==c?1:0),"cached affine resonant integration constant mismatch");
    }
    if(limits.column_progress)limits.column_progress(c+1,d,budget.products);
  }
  auto regular_trace=z;
  for(unsigned i=0;i<d;++i)regular_trace=regular_trace+a[i][i]-(z.constant(canonical.exponents[i].power)+z.constant(canonical.exponents[i].slope)*z.variable(ei))/x;
  std::vector<Exact> origin;for(unsigned j=0;j<z.variable_count();++j)origin.push_back(j==xi || j==ei?z:z.variable(j));
  std::optional<Exact> wronskian;
  if(!regular_trace.denominator().substitute(origin).is_zero())wronskian=affine_frobenius_detail::determinant(canonical.transform);
  require(state.expansion.wronskian_prefactor==wronskian,"cached affine Wronskian metadata mismatch");
}
} // namespace detail

// Narrow internal access to private series state. The public cache path always
// verifies restored/assembled columns before returning a usable series.
struct StateAccess {
  static detail::State release(Series&& series) {
    return {std::move(series.frame_),std::move(series.eigen_),std::move(series.successor_),
      std::move(series.frontiers_),std::move(series.expansion_)};
  }
  static Series build_unverified_columns(const Matrix& a,std::size_t xi,std::size_t ei,unsigned n,
      const Series::Options& options,
      std::function<std::optional<std::vector<Series::Term>>(unsigned)> load,
      std::function<void(unsigned,std::span<const Series::Term>)> save) {
    Series::ColumnCallbacks callbacks{std::move(load),std::move(save)};
    return Series::prepare_impl(a,xi,ei,n,options,&callbacks);
  }
  static Series restore(const Matrix& a,std::size_t xi,std::size_t ei,unsigned n,const Series::Options& options,
      detail::State state,const VerificationLimits& limits={}) {
    detail::verify(a,xi,ei,n,options,state,limits);
    Series out(a[0][0].constant(0));out.d_=a.size();out.n_=n;out.xi_=xi;out.ei_=ei;out.options_=options;
    out.frame_=std::move(state.frame);out.eigen_=std::move(state.exponents);out.successor_=std::move(state.successors);
    out.frontiers_=std::move(state.frontiers);out.expansion_=std::move(state.expansion);return out;
  }
};
inline artifacts::Identity identity(const Matrix& a,std::size_t xi,std::size_t ei,unsigned n,const Series::Options& options={}) {
  detail::input(a,xi,ei,n,options);json::array symbols,basis;
  for(const auto& symbol:a[0][0].variables())symbols.emplace_back(symbol);
  for(unsigned i=0;i<a.size();++i)basis.emplace_back(i);
  artifacts::Identity id;id.kind="exact_equation";id.algorithm_version=detail::algorithm;
  id.family={{"source_matrix",detail::matrix_json(a)}};id.ordered_basis=std::move(basis);
  id.normalization={{"frame","canonical exact affine Jordan chains"},{"resonance","zero integration constants except seeded k=0 diagonal"},{"frontier","common absolute cutoff for integer-related powers"}};
  id.branch={{"status","formal powers and logarithms; no evaluated branch"}};
  id.geometry={{"expansion_origin","0"},{"x_variable_index",static_cast<std::int64_t>(xi)},{"epsilon_variable_index",static_cast<std::int64_t>(ei)}};
  id.boundary={{"status","not-applicable; exact fundamental series"}};
  id.scientific_inputs={{"ordered_field_symbols",symbols},{"retained_x_order",n}};
  id.json_value();return id;
}
inline Result prepare_legacy(const Matrix& a,std::size_t xi,std::size_t ei,unsigned n,const Series::Options& options,
    artifacts::Store& store,const VerificationLimits& limits={}) {
  auto key=identity(a,xi,ei,n,options);const Demand resources{0,0,n,64,0};
  const artifacts::Certificate exact{"exact",detail::verifier,detail::scope,
    {{"check","canonical frame/frontier, polynomial ODE residual, resonant constants and Wronskian"},
     {"scope","exact retained coefficients only; no omitted Taylor or epsilon tails"}}};
  auto hit=store.lookup(key,resources,{"exact",detail::verifier,detail::scope});
  if(!hit)hit=store.lookup(key,resources);
  if(hit) {
    auto state=detail::decode(hit->payload,a[0][0],a.size(),n,options);
    auto restored=StateAccess::restore(a,xi,ei,n,options,std::move(state),limits);
    if(!artifacts::CertificateRequirement{"exact",detail::verifier,detail::scope}.accepts(hit->certificate))
      hit=store.put(key,resources,hit->payload,exact);
    return {std::move(restored),true,key.key(),hit->content_id};
  }
  auto fresh=Series::prepare(a,xi,ei,n,options);auto payload=detail::payload(fresh);
  // Persist preparation before verification, so an interrupted verifier cannot
  // lose the expensive exact solve. The checkpoint makes no exact certificate
  // claim, and a restart must still execute the complete mathematical verifier.
  const artifacts::Certificate checkpoint{"uncertified_numerical","","unverified_exact_preparation_checkpoint",
      {{"status","exact-data checkpoint pending independent residual verification; no error or tail claim"}}};
  (void)store.put(key,resources,payload,checkpoint);
  auto state=detail::decode(payload,a[0][0],a.size(),n,options);
  detail::verify(a,xi,ei,n,options,state,limits);
  auto stored=store.put(key,resources,payload,exact);
  return {std::move(fresh),false,key.key(),stored.content_id};
}
// Each completed column is durable independently. A small manifest binds the
// ordered concrete column records; only the complete independently verified
// series receives an exact certificate. Column records alone are checkpoints.
inline artifacts::Identity column_identity(const artifacts::Identity& base,unsigned column) {
  auto id=base;id.algorithm_version+="-column-v1";id.scientific_inputs["column_index"]=column;return id;
}
inline artifacts::Identity manifest_identity(const artifacts::Identity& base,unsigned dimension) {
  auto id=base;id.algorithm_version+="-column-manifest-v1";
  for(unsigned c=0;c<dimension;++c)id.parents.push_back({"column"+std::to_string(c),column_identity(base,c).key()});
  return id;
}
namespace detail {
inline json::object column_payload(unsigned column,std::span<const Series::Term> terms) {
  json::array entries;for(const auto& term:terms)entries.push_back(term_json(term));
  return {{"schema","DiffExp3.ExactAffineFrobeniusColumn/v1"},{"column",column},{"terms",std::move(entries)}};
}
inline std::vector<Series::Term> decode_column(const json::value& value,const Exact& sample,unsigned d,
    unsigned n,unsigned column,std::size_t term_limit) {
  const auto& o=value.as_object();artifacts::detail::keys(o,{"schema","column","terms"});
  require(artifacts::detail::string(o.at("schema"))=="DiffExp3.ExactAffineFrobeniusColumn/v1" &&
    artifacts::detail::integer(o.at("column"))==column,"cached affine column schema/index mismatch");
  const auto& terms=o.at("terms").as_array();require(terms.size()<=term_limit,"cached affine column term budget exhausted");
  std::map<std::string,Exact> variables;for(unsigned j=0;j<sample.variable_count();++j)variables.emplace(sample.variables()[j],sample.variable(j));
  std::vector<Series::Term> result;result.reserve(terms.size());
  for(const auto& entry:terms) {
    const auto& t=entry.as_object();artifacts::detail::keys(t,{"row","column","log_degree","power","slope","coefficient"});
    const auto row=artifacts::detail::integer(t.at("row")),c=artifacts::detail::integer(t.at("column")),l=artifacts::detail::integer(t.at("log_degree"));
    require(row>=0 && row<d && c==column && l>=0 && static_cast<std::size_t>(l)<=static_cast<std::size_t>(d)*(n+1),"cached affine column term index/log budget");
    auto coefficient=evaluate_exact(data::Reader(artifacts::detail::string(t.at("coefficient"))).read(),sample,variables);
    require(!coefficient.is_zero(),"cached affine column zero term");
    result.push_back({static_cast<unsigned>(row),column,static_cast<unsigned>(l),
      Rational(artifacts::detail::string(t.at("power"))),Rational(artifacts::detail::string(t.at("slope"))),std::move(coefficient)});
  }
  return result;
}
inline json::object manifest_payload(const Series& series) {
  auto out=payload(series,false);out.erase("terms");out["schema"]="DiffExp3.ExactAffineFrobeniusColumnManifest/v1";return out;
}
inline State decode_manifest(const json::value& value,const Exact& sample,unsigned d,unsigned n,const Series::Options& options) {
  auto metadata=value.as_object();require(!metadata.contains("terms") &&
    artifacts::detail::string(metadata.at("schema"))=="DiffExp3.ExactAffineFrobeniusColumnManifest/v1","cached affine manifest schema mismatch");
  metadata["schema"]="DiffExp3.ExactAffineFrobenius/v1";metadata["terms"]=json::array{};
  return decode(metadata,sample,d,n,options);
}
} // namespace detail
inline Result prepare(const Matrix& a,std::size_t xi,std::size_t ei,unsigned n,const Series::Options& options,
    artifacts::Store& store,const VerificationLimits& limits={}) {
  const auto base=identity(a,xi,ei,n,options);const Demand resources{0,0,n,64,0};
  // Preserve previously completed v1 series and their independent verification.
  if(store.lookup(base,resources))return prepare_legacy(a,xi,ei,n,options,store,limits);
  const unsigned d=a.size();const auto key=manifest_identity(base,d);
  const artifacts::Certificate exact{"exact",detail::verifier,detail::scope,
    {{"check","canonical frame/frontier, polynomial ODE residual, resonant constants and Wronskian"},
     {"scope","exact retained coefficients only; no omitted Taylor or epsilon tails"}}};
  const artifacts::Certificate checkpoint{"uncertified_numerical","","unverified_exact_preparation_checkpoint",
    {{"status","exact-data checkpoint pending independent residual verification; no error or tail claim"}}};
  const artifacts::CertificateRequirement verified{"exact",detail::verifier,detail::scope};
  auto hit=store.lookup(key,resources,verified);if(!hit)hit=store.lookup(key,resources);
  if(hit) {
    auto state=detail::decode_manifest(hit->payload,a[0][0],d,n,options);
    for(unsigned c=0;c<d;++c) {
      auto record=store.read(column_identity(base,c),hit->parent_content_ids.at(c));
      detail::require(record.guarantee.dominates(resources),"cached affine column resources mismatch");
      auto terms=detail::decode_column(record.payload,a[0][0],d,n,c,options.max_terms-state.expansion.terms.size());
      for(auto& term:terms)state.expansion.terms.push_back(std::move(term));
    }
    auto series=StateAccess::restore(a,xi,ei,n,options,std::move(state),limits);
    if(!verified.accepts(hit->certificate))hit=store.put(key,resources,hit->payload,exact,hit->parent_content_ids);
    return {std::move(series),true,key.key(),hit->content_id,d,0};
  }
  std::vector<std::string> contents(d);unsigned reused=0,prepared=0;std::size_t loaded_terms=0;
  auto load=[&](unsigned c)->std::optional<std::vector<Series::Term>> {
    auto record=store.lookup(column_identity(base,c),resources);if(!record)return std::nullopt;
    auto terms=detail::decode_column(record->payload,a[0][0],d,n,c,options.max_terms-loaded_terms);
    loaded_terms+=terms.size();contents[c]=record->content_id;++reused;return terms;
  };
  auto save=[&](unsigned c,std::span<const Series::Term> terms) {
    auto record=store.put(column_identity(base,c),resources,detail::column_payload(c,terms),checkpoint);
    contents[c]=record.content_id;loaded_terms+=terms.size();++prepared;
  };
  auto fresh=StateAccess::build_unverified_columns(a,xi,ei,n,options,load,save);
  auto payload=detail::manifest_payload(fresh);
  (void)store.put(key,resources,payload,checkpoint,contents);
  auto restored=StateAccess::restore(a,xi,ei,n,options,StateAccess::release(std::move(fresh)),limits);
  auto record=store.put(key,resources,payload,exact,contents);
  return {std::move(restored),false,key.key(),record.content_id,reused,prepared};
}
} // namespace diffexp::cached_affine
