#pragma once
#include "diffexp/recursion_graph.hpp"
#include "diffexp/recursion_leaf.hpp"
#include "diffexp/laurent_transport.hpp"
#include "diffexp/affine_matching.hpp"
#include "diffexp/causal.hpp"
#include "diffexp/adjoint_transport.hpp"
#include "diffexp/affine_operator.hpp"
#include "diffexp/linear_boundary.hpp"
#include "diffexp/geometry.hpp"
#include <sstream>
#include "diffexp/factored_transport.hpp"
#include "diffexp/adjoint_checkpoint.hpp"
#include "diffexp/direct_adjoint_endpoint.hpp"
#include "diffexp/recursion_physical_basis.hpp"
#include "diffexp/physical_basis_scan.hpp"
#include "diffexp/level_cache.hpp"
#include "diffexp/ft_functional_reduction.hpp"
#include "diffexp/exact_jet_composition.hpp"
#include "diffexp/endpoint_boundary.hpp"
#include <map>

namespace diffexp::recursion {
enum class LinearMethod { adjoint, factored, automatic };
// How the chain of local operators is evaluated. `rows`: bottom-up, every
// level transports one row per integral its parent needs and composes the
// result with its child's expression. `weighted`: top-down (reverse-mode
// accumulation), every level transports one row per top-level target, the
// combination of its rows with the weights its parent hands down, and passes
// the weights on its masters to its child; the leaf values are contracted
// last. Applies to the adjoint method with few top-level targets; any other
// request is evaluated by rows.
enum class Accumulation { rows, weighted };
struct BasisSelection;
struct NumericalOptions {
  bool component_epsilon_demands=true;
  bool exact_functional_reduction=true;
  bool automatic_basis_scan=true;
  double basis_scan_seconds=30;
  std::function<void(const BasisSelection&)> basis_selection_observer;
  // Exact, explicitly supplied internal bases, keyed by zero-based FT depth.
  // Every transform is checked against the graph before endpoint preparation.
  std::map<std::size_t,PhysicalBasisTransform> physical_bases;
  // A supplied basis is final by default. With this flag the automatic scan
  // starts from it instead of from the level's masters (a basis supplied to
  // remove moving poles is not yet the smallest one); its choice is then
  // not persisted, because the cache identity names the level alone.
  bool scan_physical_bases=false;
  slong working_bits=384;
  // The exact endpoint series is summed at the overlap point, at most
  // 1/overlap_ratio of the way to the nearest other singular point, where the
  // transport takes over: its tail is overlap_ratio^-endpoint_order. Each
  // endpoint has its own overlap point (EndpointGeometry). A wider
  // ratio buys a shorter series (whose exact cost grows faster than its
  // order) for a few more charts next to the endpoint. Zero selects the pair
  // for a 2^-128 tail: 64 and 22, or 16 and 32 for a recursion specialised
  // to a numeric dimension, whose exact work is cheap. With only the order
  // given the ratio follows from it, between 16 and 64.
  unsigned endpoint_order=0,overlap_ratio=0;
  unsigned ordinary_order=80,max_epsilon=100,max_refinements=8,max_overlap_halvings=128;
  // Most orders of the leaf's Laurent series at epsilon=0 summed beyond the
  // requested top when the recursion is expanded around a shifted epsilon.
  unsigned shifted_leaf_orders=60;
  Rational overlap{Rational("1/16")},contour_height{Rational("1/10")};
  unsigned leaf_digits=28;
  // Consistency of retained physical endpoint constraints, not a proof about
  // omitted epsilon/Taylor coefficients or exact vanishing of a period.
  Rational endpoint_constraint_tolerance{Rational("1/100000000000000000000")};
  unsigned max_endpoint_constraints=256;
  // Guarded experimental route: build endpoint integration rows directly.
  // Unsupported fixed sectors retain the complete Frobenius construction.
  bool direct_endpoint_functionals=false;
  LinearMethod linear_method=LinearMethod::adjoint;
  // The library default keeps one row per requested integral; `diffexp ft`
  // selects weighted. More top-level targets than weighted_max_targets are
  // evaluated by rows, which then costs no more.
  Accumulation accumulation=Accumulation::rows;
  unsigned weighted_max_targets=4;
  // The recursion propagates inherited uncertainty by probes (see
  // AdjointOptions::Uncertainty); the transport's own default is the map.
  AdjointOptions adjoint=[]{AdjointOptions value;value.uncertainty=AdjointOptions::Uncertainty::probe;return value;}();
  std::optional<causal::Prescription> causal_prescription;
  std::function<void(std::size_t,const std::string&,int)> progress;
  // Accepted fraction of the current ordinary arm's original path. This is
  // geometric progress, not an estimate of remaining computing time.
  std::function<void(double)> ordinary_fraction_progress;
  std::function<void(const ExactEpsilonMatrix&,const LaurentRows&,
      const ExactEpsilonMatrix&,const std::vector<Exact>&)> ordinary_problem_observer;
  // Optional read-only diagnostics on completed coefficient maps. This lets
  // callers locate arithmetic loss without materializing the shared source.
  std::function<void(std::size_t,const std::string&,const LaurentRows&)> operator_observer;
  // Persists the exact automatic internal-basis choice of every level, so a
  // resumed run transports in the same coordinates as its checkpoints. Empty
  // disables it; explicit physical_bases always take priority.
  std::filesystem::path basis_cache_directory;
  std::filesystem::path ordinary_cache_directory;
  std::size_t ordinary_cache_max_bytes=64*1024*1024;
  std::size_t basis_cache_max_bytes=64*1024*1024;
  affine_matching::Options matching=[] {affine_matching::Options value;value.max_dimension=256;value.max_epsilon_depth=512;value.numeric_operator_convolution=true;return value;}();
  AffineFrobeniusSeries::Options endpoint=[] {AffineFrobeniusSeries::Options value;value.max_dimension=256;return value;}();
  fuchsify::Options fuchsification=[] {fuchsify::Options value;value.max_dimension=256;return value;}();
  gaussian::Budget gaussian;
};
struct EndpointGeometry {
  // Each endpoint has its own overlap, the distance from it of the point
  // where its exact series is summed: it follows from the distance of that
  // endpoint to the nearest other singular point (lower_/upper_nearest_pole,
  // a lower bound; zero if its coefficients have no such point). A singular
  // point next to one end does not shorten the series' reach at the other.
  Rational lower_overlap,upper_overlap;
  // Of both endpoints, each in the coordinate that has its endpoint at zero.
  std::vector<kernel::ComplexBall> nonzero_poles;
  // The smaller of the two distances.
  NativeTailMagnitude nearest_pole_lower;
  NativeTailMagnitude lower_nearest_pole,upper_nearest_pole;
  const Rational& overlap(bool upper)const {return upper?upper_overlap:lower_overlap;}
};
struct NumericalStatistics {
  std::size_t component_transports=0,component_slots=0,component_rectangular_slots=0;
  std::size_t functional_reduction_attempts=0,functional_reduction_certificates=0;
  std::size_t basis_scans=0,basis_selections=0,basis_cache_hits=0,basis_cache_writes=0;
  adjoint_checkpoint::Statistics ordinary_checkpoints;
  std::size_t exact_plans=0,leaf_evaluations=0,numeric_evaluations=0,refinements=0,cache_hits=0;
  std::size_t endpoint_series_built=0;
  std::size_t direct_endpoints=0;
  std::size_t adjoint_selections=0,factored_selections=0;
  // Levels transported as weighted rows, and evaluations that asked for
  // weighted accumulation but were evaluated by rows.
  std::size_t weighted_levels=0,weighted_fallbacks=0;
  std::size_t demand_preflights=0,local_operator_reuses=0;
  std::size_t conditioning_method_fallbacks=0;
  std::size_t endpoint_constraint_rows=0,endpoint_constraint_coefficients=0;
  double maximum_endpoint_constraint_residual=0;
  double ordinary_seconds=0;
  // Wall time of the phases that precede transport. basis_selection_seconds
  // covers discovery, screening and installation, not the original plan.
  double endpoint_series_seconds=0,endpoint_operator_seconds=0,basis_selection_seconds=0,leaf_seconds=0;
  // The rest of an endpoint's exact preparation: the gauge to Fuchsian form,
  // the frames (the series in the masters' coordinates), and the projection
  // and regulated integration of every observable row.
  double endpoint_gauge_seconds=0,endpoint_frame_seconds=0,endpoint_observable_seconds=0;
  // Local-operator coefficients below a child component's structural window:
  // rounding noise replaced by the exact zero, and (never expected) values
  // excluding zero, which are kept together with their wider child demand.
  std::size_t structural_zero_coefficients=0,structural_zero_conflicts=0;
};
// A boundary-independent local problem. Endpoint rows act on the
// epsilon-gauged master vector. Each arm satisfies g'=forcing-g*connection
// along its own path; the upper forcing sign is already included.
struct PreparedAdjointStage {
  ExactEpsilonMatrix connection,lower_forcing,upper_forcing;
  LaurentRows lower_endpoint,upper_endpoint;
  std::vector<Exact> lower_path,upper_path;
  std::vector<std::int64_t> epsilon_gauge_shifts;
};
struct PreparedFactoredConsumers {
  ExactEpsilonMatrix connection,beta;
  LaurentRows right;
};
struct EpsilonDemandProfile {
  std::size_t depth,masters;
  int lower_endpoint_low,upper_endpoint_low,beta_low,anchor_inverse_pole_loss;
  int adjoint_child_loss,factored_child_loss;
  epsilon_demand::Highs child_component_loss;
};
struct BasisSelection {
  std::size_t depth=0,masters=0,original_characters=0,candidate_characters=0;
  unsigned candidates_checked=0;
  unsigned candidates_screened=0,full_candidates_checked=0;
  int original_child_loss=0,candidate_child_loss=0;
  bool accepted=false,cache_hit=false;
  double scan_seconds=0,total_selection_seconds=0;
  std::string reason;
};

// Graph-independent option contract of the recursive evaluator. A caller that
// prepares the exact graph first (minutes of IBP reduction) checks this BEFORE
// that work, so a mistyped numerical option is not reported after it.
inline void validate_numerical_options(const NumericalOptions& options) {
  if(options.adjoint.continuation)
    throw std::invalid_argument("a continuation belongs to one adjoint arm; use ordinary_cache_directory for recursive continuation");
  const auto fail=[](const std::string& what){throw std::invalid_argument("recursive numerical options or finite budgets: "+what);};
  if(options.working_bits<64 || options.working_bits>1000000)fail("working precision must be 64..1000000 bits");
  if(options.overlap_ratio && (options.overlap_ratio<4 || options.overlap_ratio>65536))fail("overlap ratio must lie in 4..65536");
  if(options.ordinary_order<8 || options.ordinary_order>1000)fail("ordinary order must be 8..1000");
  if(!options.max_epsilon || options.max_epsilon>100)fail("epsilon budget must be 1..100");
  if(!options.max_refinements || !options.max_overlap_halvings)fail("refinement and overlap-halving budgets must be positive");
  if(options.overlap<=Rational(0) || options.overlap>=Rational("1/2"))fail("overlap must lie strictly between 0 and 1/2");
  if(options.contour_height<=Rational(0))fail("contour height must be positive");
  if(options.endpoint_constraint_tolerance<=Rational(0))fail("endpoint constraint tolerance must be positive");
  if(!options.max_endpoint_constraints || options.max_endpoint_constraints>5000)fail("endpoint constraint budget must be 1..5000");
  if(!std::isfinite(options.basis_scan_seconds) || options.basis_scan_seconds<=0 || options.basis_scan_seconds>60)
    fail("basis scan seconds must lie in (0,60]");
  if(options.leaf_digits>100000)fail("leaf digits must not exceed 100000");
}

// Only exact coordinates are durable. Endpoint acceptance, epsilon demands,
// precision and approximation scores are deliberately absent from the payload.
namespace basis_persistence {
namespace json=boost::json;
// v5: candidates are ranked by the size of their exact connection; the adopted
// one is installed by exact pullback of the original endpoint functionals and
// must not widen the child's epsilon window. Undecided discoveries (wall-clock
// cut, error) are not stored as "original retained".
inline constexpr const char* revision="physical-basis-connection-size-v5";
inline constexpr const char* verifier="physical-basis-two-inverses-and-gauge-v1";
inline constexpr const char* scope="exact_physical_basis_change";
inline artifacts::Identity identity(const Graph& graph,std::size_t depth,
    std::size_t xi,std::size_t ei,const causal::Prescription& prescription) {
  const auto& node=graph.nodes.at(depth);
  if(node.scalar_leaf || node.closure.matrix.empty() || ei>=graph.field.variables().size())
    throw std::invalid_argument("physical basis cache requires a nonleaf closure and epsilon variable");
  // Reuse the exact propagator/domain identity, including derived coordinates.
  auto id=cached_level::identity(node.merged,graph.dimension,graph.field,xi,
                                 node.closure.ordered_basis,level::Options{});
  id.algorithm_version=revision;
  id.ordered_basis=cached_level::detail::integrals_json(node.closure.ordered_basis);
  const auto hash=[](const json::value& value) {
    return artifacts::detail::sha256(artifacts::detail::canonical(value));
  };
  id.scientific_inputs["original_connection_sha256"]=hash(cached_level::detail::matrix_json(node.closure.matrix));
  id.scientific_inputs["observable_rows_sha256"]=hash(cached_level::detail::matrix_json(node.observable_rows));
  id.scientific_inputs["epsilon_variable"]=graph.field.variables()[ei];
  id.scientific_inputs["selection_policy"]="smallest exact connection that does not widen the child epsilon window; exact pullback of the original endpoint functionals";
  // Scan runtime is intentionally NOT an identity axis: restart must reuse the
  // completed choice despite a different machine or wall-clock allowance.
  id.scientific_inputs["scan_revision"]=revision;
  json::array operations;
  for(const auto& op:node.operations)
    operations.push_back(json::object{{"operation",static_cast<int>(op.operation)},
      {"left_power",op.left_power},{"right_power",op.right_power},
      {"normalization",op.normalization.str()}});
  id.scientific_inputs["operations"]=std::move(operations);
  id.geometry["depth"]=depth;id.geometry["left"]=node.left;id.geometry["right"]=node.right;
  id.geometry["anchor"]=node.anchor.str();
  id.branch={{"convention",prescription.convention},{"f_rim",prescription.f_rim},
             {"x_detour_sign",prescription.levels.at(depth).x_detour_sign}};
  // Derived orientations are a function of the level's propagators; what was
  // derived at its pinch candidates is part of the identity all the same.
  if(prescription.derive_orientations) {
    std::string derived;
    try {derived=causal::Landau(node.merged,graph.field.variables()[xi]).signature(prescription.f_rim);}
    catch(const std::exception& error) {derived=std::string(causal::derivation_revision)+";unavailable: "+error.what();}
    id.branch["derived_orientations"]=derived;
  }
  if(!prescription.points.empty()) {
    json::array stated;
    for(const auto& point:prescription.points)if(point.level==depth+1)
      stated.push_back(json::object{{"x",point.x.str()},{"side",point.side}});
    id.branch["point_sides"]=std::move(stated);
  }
  id.boundary={{"status","basis selection independent of numerical boundary and requested epsilon high"}};
  id.json_value();return id;
}
struct Choice {PhysicalBasisTransform basis;bool selected=false;};
inline json::object payload(const Choice& choice) {
  return {{"schema","DiffExp21.ExactPhysicalBasis/v1"},{"selected",choice.selected},
    {"original_connection",cached_level::detail::matrix_json(choice.basis.original_connection)},
    {"to_basis",cached_level::detail::matrix_json(choice.basis.to_basis)},
    {"from_basis",cached_level::detail::matrix_json(choice.basis.from_basis)},
    {"connection",cached_level::detail::matrix_json(choice.basis.connection)}};
}
inline void verify(const Choice& choice,const ExactEpsilonMatrix& original,std::size_t xi) {
  verify_physical_basis(choice.basis,original,xi);
  if(!choice.selected) {
    const auto unit=fuchsify::detail::identity(original.size(),original[0][0]);
    if(choice.basis.to_basis!=unit || choice.basis.from_basis!=unit || choice.basis.connection!=original)
      throw std::invalid_argument("cached original-basis decision has a nonidentity transformation");
  }
}
inline Choice decode(const json::value& value,const ExactEpsilonMatrix& original,std::size_t xi) {
  if(original.empty() || original.size()>256 || original[0].empty())
    throw std::invalid_argument("physical basis cache original dimension");
  const auto& object=value.as_object();
  artifacts::detail::keys(object,{"schema","selected","original_connection","to_basis","from_basis","connection"});
  if(artifacts::detail::string(object.at("schema"))!="DiffExp21.ExactPhysicalBasis/v1")
    throw std::invalid_argument("physical basis cache schema");
  auto matrix=[&](const char* name) {
    const auto& encoded=object.at(name).as_array();
    if(encoded.size()!=original.size())throw std::invalid_argument("physical basis cache matrix height");
    ExactEpsilonMatrix result;
    for(const auto& row:encoded) {
      if(row.as_array().size()!=original.size())throw std::invalid_argument("physical basis cache matrix width");
      result.emplace_back();
      for(const auto& item:row.as_array())result.back().push_back(original[0][0].parse(artifacts::detail::string(item)));
    }
    return result;
  };
  Choice choice{{matrix("original_connection"),matrix("to_basis"),matrix("from_basis"),matrix("connection")},
                object.at("selected").as_bool()};
  verify(choice,original,xi);return choice;
}
inline std::optional<Choice> load(const artifacts::Store& store,const artifacts::Identity& id,
    const ExactEpsilonMatrix& original,std::size_t xi) {
  if(auto record=store.lookup(id,Demand{0,0,0,64,0},{"exact",verifier,scope}))
    return decode(record->payload,original,xi);
  return std::nullopt;
}
inline void save(artifacts::Store& store,const artifacts::Identity& id,
    const Choice& choice,const ExactEpsilonMatrix& original,std::size_t xi) {
  verify(choice,original,xi);
  // Immutable content-addressed publication includes fsync and atomic link.
  // Never replace a completed choice because a later scan happened to differ.
  if(load(store,id,original,xi))return;
  store.put(id,Demand{0,0,0,64,0},payload(choice),
    artifacts::Certificate{"exact",verifier,scope,
      {{"check","both inverse products and exact differential gauge identity"},
       {"scope","exact coordinates only; endpoint and epsilon demand gates must be recomputed"}}});
}
} // namespace basis_persistence

// Executes an immutable exact request graph. Only child numerical epsilon
// demands grow during refinement; graph closure and endpoint plans are reused.
// Nonleaf results are retained local-series approximations and deliberately do
// not claim omitted-tail certificates. The F rim and every level's contour
// orientation are explicit physical input; unsupported endpoint sectors and
// failed matching stop clearly.
class Evaluator {
  using Matrix=ExactEpsilonMatrix;
  using Expansion=AffineFrobeniusSeries::Expansion;
  using B=kernel::ComplexBall;
  struct PrecisionScope {
    slong previous;
    explicit PrecisionScope(slong bits):previous(B::precision()){B::set_precision(bits);}
    ~PrecisionScope(){B::set_precision(previous);}
  };
  struct Endpoint {
    std::optional<AffineFrobeniusSeries> series;
    Expansion matching_frame;
    std::vector<std::optional<Expansion>> functionals;
    std::vector<Exact> clearance_coefficients;
    std::vector<std::pair<Expansion,std::string>> pending_constraints;
    std::optional<Matrix> direct_operator;
    bool direct_operator_exact=false;
    // The numerical inverse frame at the cut. It is prepared once and only
    // extended to higher epsilon orders afterwards.
    std::optional<affine_operator::Prepared> inverse;
    std::optional<Rational> inverse_point;
    // Of the endpoint that carries the level's boundary: the fundamental
    // solution in the masters' own coordinates, and the order of the pole at
    // the endpoint of the gauge that leads back from them.
    std::optional<Expansion> physical_frame;
    unsigned physical_inverse_pole_order=0;
  };
  struct SingularPoint {
    Rational low,high;
    int side=0;
    // Why no side could be derived; such a point may not lie on a leg.
    std::string undetermined;
  };
  struct ComplexPoint {Rational low,high,height;};
  struct Plan {
    EpsilonGaugeResult gauge;
    Matrix diagonal,inverse_diagonal;
    Endpoint lower,upper;
    Matrix beta_rows;
    std::vector<std::size_t> beta_indices;
    std::int64_t largest_shift=0;
    EndpointGeometry geometry;
    std::optional<LaurentRows> lower_operator,upper_operator,local_operator;
    std::optional<int> adjoint_child_loss,factored_child_loss;
    bool factored_conditioning_failed=false;
    // Rational enclosures of the real singular points of the transported
    // connection and its forcing, ascending, each with the side on which the
    // contour passes it (zero: free); see path().
    std::optional<std::vector<SingularPoint>> real_singularities;
    // The complex singular points next to the segment: enclosures of the
    // real part and a lower bound of |Im|. Kept when points are treated one
    // by one, to keep them out of every detour.
    std::vector<ComplexPoint> complex_singularities;
    std::vector<std::string> constraint_labels;
    bool has_physical_basis=false;
    epsilon_demand::Highs component_child_loss;
    std::optional<Matrix> rewritten_lower,rewritten_upper;
    // An internal basis change Y=D_base U_base=D U leaves every endpoint
    // functional unchanged: E=E_base V(cut) with the exact V=D_base^-1 D.
    // Such a plan owns the plan whose endpoint frames it reuses; its own
    // lower/upper members only carry the functional count.
    std::unique_ptr<Plan> pullback_base;
    Matrix pullback_map;
    // Exact valuations read by every demand propagation of this plan.
    std::optional<epsilon_demand::Valuations> connection_valuations,anchor_inverse_valuations;
    // A level with its boundary at an endpoint: the masters at the overlap
    // point as an exact matrix applied to their values at the endpoint, which
    // the child provides. It is stated in the masters' coordinates and so is
    // the same for every internal basis of the level.
    std::optional<Matrix> boundary_map;
  };
  struct DirectDemand {int low=0;epsilon_demand::Highs needs;};
  struct Seconds {
    double& total;std::chrono::steady_clock::time_point start=std::chrono::steady_clock::now();
    explicit Seconds(double& destination):total(destination){}
    Seconds(const Seconds&)=delete;Seconds& operator=(const Seconds&)=delete;
    ~Seconds(){total+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
  };
 public:
  explicit Evaluator(const Graph& graph,NumericalOptions options={})
      :graph_(graph),options_(std::move(options)),plans_(graph.nodes.size()),cache_(graph.nodes.size()),expressions_(graph.nodes.size()),direct_demands_(graph.nodes.size()) {
    validate_numerical_options(options_);
    if(graph_.nodes.empty())throw std::invalid_argument("recursive numerical options or finite budgets: empty recursion graph");
    // The scalar leaf is asked for leaf_digits digits, and nothing built on it
    // is more accurate: a Taylor tail ten digits below that per chart keeps
    // the transport from limiting the result without computing charts for
    // digits the leaf does not have.
    if(!options_.adjoint.truncation_bits)options_.adjoint.truncation_bits=static_cast<unsigned>(std::min<double>(
      std::min<double>(2.0*options_.ordinary_order,static_cast<double>(options_.working_bits)-32),
      std::ceil((options_.leaf_digits+10.0)*std::log2(10.0))));
    // The endpoint constraints are checked on transported values: with a
    // Taylor tail admitted above 2^-87 the check cannot hold to 1e-20, and
    // its tolerance follows the tail at half its bits: a level of 112 rows
    // on 47 masters was seen to lose 25 bits against the tail.
    if(options_.adjoint.truncation_bits<87) {
      const unsigned bits=std::max(8u,options_.adjoint.truncation_bits/2);
      std::string power="1";for(unsigned i=0;i<bits;++i){int carry=0;for(auto digit=power.rbegin();digit!=power.rend();++digit){const int v=(*digit-'0')*2+carry;*digit=static_cast<char>('0'+v%10);carry=v/10;}if(carry)power.insert(power.begin(),'1');}
      const Rational loose("1/"+power);
      if(options_.endpoint_constraint_tolerance<loose)options_.endpoint_constraint_tolerance=loose;
    }
    if(graph_.specialised)options_.matching.numeric_leading_gauge=true;
    if(!options_.overlap_ratio) {
      options_.overlap_ratio=graph_.specialised?16:64;
      if(options_.endpoint_order)
        options_.overlap_ratio=options_.endpoint_order>=32?16:options_.endpoint_order>=26?32:64;
    }
    if(!options_.endpoint_order)options_.endpoint_order=static_cast<unsigned>(std::ceil(128.0/std::log2(static_cast<double>(options_.overlap_ratio))));
    auto [xi,ei]=path_epsilon_variables(graph_.dimension);xi_=xi;ei_=ei;
    for(const auto& [depth,basis]:options_.physical_bases)
      if(depth>=graph_.nodes.size() || graph_.nodes[depth].scalar_leaf)
        throw std::invalid_argument("general physical basis requires a nonleaf depth");
    auto d=graph_.dimension.substitute(exact_point(graph_.dimension,ei_,graph_.dimension.constant(0))).rational()
      +Rational(2)*graph_.epsilon_shift;
    if(d.str().find('/')!=std::string::npos)throw std::invalid_argument("recursive numerical base dimension must be integral");
    d0_=std::stoi(d.str());
    if(!(graph_.dimension==graph_.dimension.constant(d0_)-graph_.dimension.constant(2)*
        (graph_.dimension.variable(ei_)+graph_.dimension.constant(graph_.epsilon_shift))))
      throw std::invalid_argument("recursive numerical dimension must be d0-2eps");
    prescription_=options_.causal_prescription?*options_.causal_prescription:example_prescription();
    if(prescription_.derive_orientations)prescription_.levels.assign(graph_.nodes.size(),causal::LevelPlan{1});
    prescription_.validate(graph_.nodes.size());
    landau_.resize(graph_.nodes.size());
    if(options_.automatic_basis_scan && !options_.basis_cache_directory.empty())
      basis_store_=std::make_unique<artifacts::Store>(options_.basis_cache_directory,options_.basis_cache_max_bytes);
  }
  const NumericalStatistics& statistics()const{return statistics_;}
  // The next evaluation runs at `bits` of working precision. Every number
  // computed so far is dropped; the exact part of the endpoint plans (series,
  // frames, functionals, overlaps, internal bases, contour sides) is kept, so
  // the evaluation repeats the endpoint operators, the transports and the
  // leaf, not the exact preparation.
  void raise_precision(slong bits) {
    if(bits<=options_.working_bits)return;
    if(bits>1000000)throw std::invalid_argument("recursive numerical options or finite budgets: working precision must be 64..1000000 bits");
    options_.working_bits=bits;
    const std::function<void(Plan&)> reset=[&](Plan& prepared) {
      prepared.lower_operator.reset();prepared.upper_operator.reset();prepared.local_operator.reset();
      for(auto* endpoint:{&prepared.lower,&prepared.upper}){endpoint->inverse.reset();endpoint->inverse_point.reset();}
      if(prepared.pullback_base)reset(*prepared.pullback_base);
    };
    for(auto& prepared:plans_)if(prepared)reset(*prepared);
    for(auto& value:cache_)value.reset();
    for(auto& expression:expressions_)expression.reset();
  }
  slong working_bits()const{return options_.working_bits;}
  const std::vector<BasisSelection>& basis_selections()const{return basis_selections_;}
  // Snapshot a completed transform together with its immutable leaf source.
  // Copying the transform prevents a later refinement (or caller edit) from
  // changing the saved expression. No evaluation or materialization occurs.
  linear_boundary::Expression linear_expression(std::size_t depth=0)const {
    if(depth>=expressions_.size() || !expressions_[depth])
      throw std::logic_error("recursive linear expression is not available");
    return *expressions_[depth];
  }
  const causal::Prescription& causal_prescription()const{return prescription_;}
  // The real singular points strictly inside (0,1) of the levels evaluated so
  // far, with the side each was given (zero: free).
  struct CausalPoint {std::size_t depth;std::string x;int side;};
  std::vector<CausalPoint> causal_points()const {
    std::vector<CausalPoint> out;
    for(std::size_t depth=0;depth<plans_.size();++depth) {
      if(!plans_[depth] || !plans_[depth]->real_singularities)continue;
      for(const auto& point:*plans_[depth]->real_singularities) {
        if(!(Rational(0)<point.low) || !(point.high<Rational(1)))continue;
        B middle;arb_set_fmpq(acb_realref(middle.raw()),((point.low+point.high)/Rational(2)).raw(),128);
        out.push_back({depth,middle.real_midpoint(20),point.side});
      }
    }
    return out;
  }
  const EndpointGeometry& endpoint_geometry(std::size_t depth)const {
    if(depth>=plans_.size() || !plans_[depth])throw std::out_of_range("recursive endpoint geometry is not prepared");
    return plans_[depth]->geometry;
  }
  // Read-only consumer fixture: no child evaluation or ordinary transport.
  // The exact connection binds the returned numerical endpoint to its basis.
  PreparedFactoredConsumers prepare_factored_consumers(std::size_t depth,unsigned high) {
    if(depth>=graph_.nodes.size() || graph_.nodes[depth].scalar_leaf || high>options_.max_epsilon)
      throw std::invalid_argument("factored consumer fixture depth or order");
    require_anchor(depth,"factored consumer fixture");
    PrecisionScope precision(options_.working_bits);auto& p=plan(depth,0);
    Matrix beta(p.lower.functionals.size(),std::vector<Exact>(p.gauge.matrix.size(),graph_.dimension.constant(0)));
    for(std::size_t i=0;i<p.beta_indices.size();++i)beta[p.beta_indices[i]]=p.beta_rows[i];
    return {p.gauge.matrix,std::move(beta),cached_endpoint_operator(p,true,high)};
  }
  // Prepare one local endpoint/transport problem without evaluating a child.
  // high is the endpoint/operator epsilon top before the final inverse-D
  // gauge shift, not the requested physical integral's epsilon order.
  PreparedAdjointStage prepare_adjoint_stage(std::size_t depth,unsigned high) {
    if(depth>=graph_.nodes.size() || high>options_.max_epsilon ||
        graph_.nodes[depth].scalar_leaf || graph_.nodes[depth].closure.ordered_basis.empty())
      throw std::invalid_argument("local adjoint preparation depth or epsilon budget");
    require_anchor(depth,"standalone adjoint-stage export");
    PrecisionScope precision(options_.working_bits);
    if(!plans_[depth])report(depth,"endpoint preparation",high);
    // The standalone v1 stage contract exports the original ordered masters.
    // Automatic coordinate changes belong to the complete recursive evaluator.
    struct RestoreScan {bool& flag;bool original;~RestoreScan(){flag=original;}} restore_scan{options_.automatic_basis_scan,options_.automatic_basis_scan};
    options_.automatic_basis_scan=false;
    auto& prepared=plan(depth,high);
    if(prepared.has_physical_basis)throw std::domain_error(
      "legacy standalone adjoint-stage export cannot encode a general physical basis; use the recursive evaluator");
    if(!prepared.constraint_labels.empty())throw std::domain_error(
      "standalone adjoint-stage export cannot discharge physical endpoint constraints; use the full recursive evaluator");
    return adjoint_stage(depth,prepared,high);
  }
  // Exact structural demand inspection without evaluating the child or
  // transporting an ordinary segment. Endpoint artifacts can be reused.
  EpsilonDemandProfile epsilon_demand_profile(std::size_t depth) {
    if(depth>=graph_.nodes.size() || graph_.nodes[depth].scalar_leaf || graph_.nodes[depth].closure.ordered_basis.empty())
      throw std::invalid_argument("epsilon demand profile requires a nonleaf level");
    PrecisionScope precision(options_.working_bits);
    auto& prepared=plan(depth,0);preflight_linear_demands(depth,prepared);
    const int l=cached_endpoint_operator(prepared,false,0).low;
    const int r=cached_endpoint_operator(prepared,true,0).low;
    int q=0;for(const auto& row:prepared.beta_rows)for(const auto& value:row)
      if(!value.is_zero())q=std::min(q,checked(*exact_epsilon_valuation(value,ei_)));
    const int loss=checked(std::max<std::int64_t>(0,prepared.largest_shift));
    return {depth,prepared.gauge.matrix.size(),l,r,q,loss,
      *prepared.adjoint_child_loss,*prepared.factored_child_loss,prepared.component_child_loss};
  }
  LaurentBoundary evaluate(unsigned desired_top=0){return evaluate(0,desired_top);}
  LaurentBoundary evaluate(std::size_t depth,unsigned desired_top) {
    auto result=evaluate_impl(depth,desired_top);
    report(depth,"completed",desired_top);
    return result;
  }
  LaurentBoundary evaluate_impl(std::size_t depth,unsigned desired_top) {
    if(depth>=graph_.nodes.size() || desired_top>options_.max_epsilon)
      throw std::invalid_argument("recursive numerical depth or epsilon budget");
    if(cache_[depth] && cache_[depth]->high()>=static_cast<int>(desired_top)) {
      ++statistics_.cache_hits;report(depth,"cache hit",desired_top);return *cache_[depth];
    }
    PrecisionScope precision(options_.working_bits);
    if(depth==0 && options_.accumulation==Accumulation::weighted)
      if(auto weighted=evaluate_weighted(desired_top))return *weighted;
    const auto& node=graph_.nodes[depth];++statistics_.numeric_evaluations;
    if(node.scalar_leaf) {
      report(depth,"certified scalar leaf",desired_top);
      feynman::CertifiedDeepestBetaOptions leaf;leaf.working_bits=options_.working_bits;
      leaf.f_rim=prescription_.f_rim;
      leaf.requested_digits=options_.leaf_digits;leaf.taylor_order=options_.ordinary_order;
      auto value=[&]{Seconds timer{statistics_.leaf_seconds};return shifted_leaf(node,desired_top,leaf);}();
      ++statistics_.leaf_evaluations;
      auto result=LaurentBoundary{value.epsilon_low,std::move(value.values),value.taylor_tail_certified};
      expressions_[depth]=linear_boundary::identity(std::make_shared<const LaurentBoundary>(result),
        checked(static_cast<long>(desired_top)-result.low));
      cache_[depth]=std::move(result);return *cache_[depth];
    }
    if(node.closure.ordered_basis.empty()) {
      auto source=LaurentBoundary{0,Boundary(1,std::vector<B>(desired_top+1,B(0))),true};source.values[0][0]=B(1);
      expressions_[depth]=linear_boundary::Expression{
        LaurentRows{0,static_cast<int>(desired_top),std::vector(node.requested.size(),std::vector(1,std::vector<B>(desired_top+1,B(0))))},
        std::make_shared<const LaurentBoundary>(std::move(source))};
      cache_[depth]=LaurentBoundary{0,Boundary(node.requested.size(),std::vector<B>(desired_top+1,B(0))),true};return *cache_[depth];
    }
    if(depth+1>=graph_.nodes.size() || graph_.nodes[depth+1].requested!=node.closure.ordered_basis)
      throw std::logic_error("recursive child does not own the parent's ordered master demand");
    if(std::all_of(node.operations.begin(),node.operations.end(),[](const auto& operation){return operation.operation==feynman::Operation::Direct;})) {
      // No endpoint, no transport: the local operator is the exact rows at the
      // anchor, whose pole depth is the child's epsilon loss.
      int demand=checked(static_cast<long>(desired_top)-direct_demand(depth).low);
      for(unsigned attempt=0;attempt<options_.max_refinements;++attempt) {
        if(demand>static_cast<int>(options_.max_epsilon))throw std::runtime_error("recursive direct child epsilon demand exceeds finite budget");
        auto child=evaluate_child(depth+1,static_cast<unsigned>(demand));
        try {
          auto local=exact_laurent_rows(node.observable_rows,graph_.dimension.constant(node.anchor),checked(static_cast<long>(desired_top)-child_window_low(depth+1,child)));
          return compose_child_expression(depth,std::move(local),desired_top);
        }catch(const BoundaryDemand& request) {
          if(request.required_high<=demand)throw std::runtime_error("recursive direct demand made no progress");
          demand=request.required_high;++statistics_.refinements;report(depth,"child refinement",demand);
        }
      }
      throw std::runtime_error("recursive direct epsilon refinement budget exhausted");
    }
    if(!plans_[depth])report(depth,"endpoint preparation",desired_top);
    auto& prepared=plan(depth,desired_top);
    // The boundary map of an endpoint level is composed with transported rows.
    auto method=node.boundary==BoundarySite::anchor?options_.linear_method:LinearMethod::adjoint;
    if(method==LinearMethod::automatic) {
      const auto width=shared_source_width(depth+1);
      std::size_t connections=0,observables=0;
      for(const auto& row:prepared.gauge.matrix)for(const auto& c:row)connections+=!c.is_zero();
      for(const auto& row:prepared.beta_rows)for(const auto& c:row)observables+=!c.is_zero();
      const long double forward=static_cast<long double>(width)*(connections+observables);
      const long double adjoint=static_cast<long double>(prepared.lower.functionals.size())*connections+observables;
      const auto caps=factored_transport::Options{};
      const bool supported=width<=caps.max_leaf_width &&
          prepared.gauge.matrix.size()+prepared.lower.functionals.size()<=caps.max_augmented_dimension;
      method=!prepared.factored_conditioning_failed && supported && forward<adjoint?LinearMethod::factored:LinearMethod::adjoint;
    }
    if(method==LinearMethod::factored) {
      ++statistics_.factored_selections;
      try {return evaluate_factored(depth,desired_top,prepared);}
      catch(const ArithmeticConditioningFailure& failure) {
        if(options_.linear_method!=LinearMethod::automatic ||
            (failure.recursion_depth && *failure.recursion_depth!=depth))throw;
        // One method change per exact plan, retaining endpoints, child
        // expressions and their immutable sources. Future order extensions
        // do not repeat a route that already failed at these settings.
        prepared.factored_conditioning_failed=true;++statistics_.conditioning_method_fallbacks;
        report(depth,std::string("factored conditioning fallback to adjoint: ")+failure.what(),desired_top);
      }
    }
    ++statistics_.adjoint_selections;return evaluate_adjoint(depth,desired_top,prepared);
  }
 private:
  LaurentBoundary evaluate_child(std::size_t depth,unsigned high) {
    try {return evaluate(depth,high);}
    catch(ArithmeticConditioningFailure& failure) {
      // A failed descendant cannot be repaired by changing its parent's
      // transport direction. Attribute it before an ancestor considers retry.
      if(!failure.recursion_depth)failure.recursion_depth=depth;
      throw;
    }
  }
  const Graph& graph_;NumericalOptions options_;std::size_t xi_,ei_;int d0_;
  causal::Prescription prescription_;
  mutable std::vector<std::optional<causal::Landau>> landau_;
  NumericalStatistics statistics_;
  std::unique_ptr<artifacts::Store> basis_store_;
  std::vector<std::unique_ptr<Plan>> plans_;
  std::vector<BasisSelection> basis_selections_;
  std::vector<std::optional<LaurentBoundary>> cache_;
  std::vector<std::optional<linear_boundary::Expression>> expressions_;
  std::vector<std::optional<DirectDemand>> direct_demands_;
  void require_anchor(std::size_t depth,const char* what)const {
    if(graph_.nodes[depth].boundary!=BoundarySite::anchor)
      throw std::domain_error(std::string(what)+" does not cover a level whose boundary is at an endpoint; use the recursive evaluator");
  }
  // Where the two transported arms meet and the child's values enter: the
  // anchor, or the overlap point next to the endpoint that carries the
  // boundary, where its exact series is summed.
  Rational matching_point(std::size_t depth,const Plan& prepared)const {
    const auto& node=graph_.nodes[depth];
    return node.boundary==BoundarySite::lower?prepared.geometry.lower_overlap:
      node.boundary==BoundarySite::upper?Rational(1)-prepared.geometry.upper_overlap:node.anchor;
  }
  // Lowest epsilon order of the boundary map, at most zero.
  int boundary_map_low(const Plan& prepared)const {
    int low=0;
    if(prepared.boundary_map)for(const auto& row:*prepared.boundary_map)for(const auto& value:row)
      if(!value.is_zero())low=std::min(low,checked(*exact_epsilon_valuation(value,ei_)));
    return low;
  }
  LaurentRows ordinary_transport(const Matrix& matrix,LaurentRows initial,const Matrix& forcing,
      const std::vector<Exact>& vertices,const AdjointOptions& settings) {
    if(options_.ordinary_problem_observer)options_.ordinary_problem_observer(matrix,initial,forcing,vertices);
    // Time spent in a transport that throws (a conditioning failure followed
    // by a method fallback, a rejected budget) is still transport time.
    Seconds timer{statistics_.ordinary_seconds};
    if(options_.ordinary_fraction_progress)options_.ordinary_fraction_progress(0);
    auto observed=settings;
    if(options_.ordinary_fraction_progress) {
      const auto caller=observed.chart_observer;
      observed.chart_observer=[&,caller](unsigned leg,double from,double to,const LaurentRows& rows) {
        if(caller)caller(leg,from,to,rows);
        options_.ordinary_fraction_progress((leg+to)/static_cast<double>(std::max<std::size_t>(1,vertices.size()-1)));
      };
    }
    auto result=adjoint_checkpoint::transport(matrix,std::move(initial),forcing,vertices,observed,
      {options_.ordinary_cache_directory,options_.ordinary_cache_max_bytes,&statistics_.ordinary_checkpoints});
    if(options_.ordinary_fraction_progress)options_.ordinary_fraction_progress(1);
    return result;
  }
  // The scalar leaf is a closed form whose Laurent series at epsilon_total=0
  // the leaf integrator provides to any order. Around a shifted point its
  // Taylor coefficients follow by re-expansion: (shift+eps)^k contributes
  // binomial(k,m)*shift^(k-m) to eps^m, for negative k as well. The series is
  // cut where the retained terms have decayed below the working precision or
  // at the epsilon budget.
  LeafBoundary shifted_leaf(const Node& node,unsigned desired_top,const feynman::CertifiedDeepestBetaOptions& leaf) {
    if(graph_.epsilon_shift.is_zero())return evaluate_leaf(node,graph_.dimension,d0_,desired_top,leaf,options_.gaussian);
    const auto epsilon=graph_.dimension.variable(ei_),shift=graph_.dimension.constant(graph_.epsilon_shift);
    const auto back=exact_point(graph_.dimension,ei_,epsilon-shift);
    Node unshifted=node;
    for(auto& operation:unshifted.operations)for(auto& [integral,value]:operation.source_row)if(!value.is_zero())value=value.substitute(back);
    const auto dimension=graph_.dimension.constant(d0_)-graph_.dimension.constant(2)*epsilon;
    // The leaf series converges like (shift*2(L+1))^k for L loops (its nearest
    // gamma-function poles); sum it until that is below the leaf tolerance.
    double magnitude=0;{B value;arb_set_fmpq(acb_realref(value.raw()),graph_.epsilon_shift.raw(),64);
      magnitude=std::abs(arf_get_d(arb_midref(acb_realref(value.raw())),ARF_RND_NEAR));}
    const double convergence=magnitude*2.0*(node.merged.space.loops+1);
    unsigned orders=options_.shifted_leaf_orders;
    if(convergence>0 && convergence<0.5)
      orders=std::min<unsigned>(orders,8+static_cast<unsigned>(std::ceil((leaf.requested_digits+6)*std::log(10.)/-std::log(convergence))));
    const unsigned top=std::min<unsigned>(options_.max_epsilon,desired_top+orders);
    // The leaf integrator applies its absolute tolerance to every order of the
    // wide series, although order k enters weighted by shift^k: a deep family
    // may need fewer leaf digits here than in the unshifted recursion.
    const auto wide=[&] {
      try{return evaluate_leaf(unshifted,dimension,d0_,top,leaf,options_.gaussian);}
      catch(const std::runtime_error& error) {
        throw std::runtime_error(std::string(error.what())+" (the shifted recursion sums the leaf series through order "
          +std::to_string(top)+": lower the leaf digits or raise the Taylor order)");
      }
    }();
    const slong bits=B::precision();
    B at;{arb_set_fmpq(acb_realref(at.raw()),graph_.epsilon_shift.raw(),bits);}
    LeafBoundary result{0,Boundary(wide.values.size(),std::vector<B>(desired_top+1,B(0))),false};
    for(std::size_t row=0;row<wide.values.size();++row)for(unsigned m=0;m<=desired_top;++m) {
      auto& sum=result.values[row][m];
      for(std::size_t index=0;index<wide.values[row].size();++index) {
        const long k=wide.epsilon_low+static_cast<long>(index);
        if(wide.values[row][index].is_zero())continue;
        // binomial(k,m)=k(k-1)...(k-m+1)/m!, exact for every integer k.
        fmpq_t binomial;fmpq_init(binomial);fmpq_one(binomial);
        for(unsigned j=0;j<m;++j){fmpq_t factor;fmpq_init(factor);fmpq_set_si(factor,k-static_cast<long>(j),static_cast<ulong>(j)+1);fmpq_mul(binomial,binomial,factor);fmpq_clear(factor);}
        B term,power;
        acb_set_fmpq(term.raw(),binomial,bits);fmpq_clear(binomial);
        if(term.is_zero())continue;
        acb_pow_si(power.raw(),at.raw(),k-static_cast<long>(m),bits);
        acb_mul(term.raw(),term.raw(),power.raw(),bits);
        acb_addmul(sum.raw(),term.raw(),wide.values[row][index].raw(),bits);
      }
    }
    return result;
  }
  // Without a supplied prescription the contour is derived. A family whose
  // physical U and F have positive monomials has no cut on the open simplex:
  // every real singular point is free and takes the default side. Any other
  // kinematics orient each real singular point from the level's own F.
  causal::Prescription example_prescription()const {
    // A familiar name cannot authorize a prescription for altered kinematics.
    // Compare the actual incoming coordinates with the definition before
    // using its positivity proof or the supplied Henn contour data.
    std::optional<feynman::ExampleFamily> example;
    if(graph_.definition)example=*graph_.definition;
    else try{example=feynman::example_family(graph_.family_name);}catch(const std::exception&){}
    if(example) {
      const ibp::PropagatorBasis expected(ibp::quadratic_family(example->momenta,graph_.dimension,example->physical_count));
      const auto& actual=graph_.nodes.front().incoming;
      bool same=actual.physical_count==expected.physical_count && actual.space.loops==expected.space.loops &&
        actual.space.external_gram==expected.space.external_gram && actual.denominators.size()==expected.denominators.size() &&
        d0_==example->dimension_at_epsilon_zero;
      if(same)for(std::size_t i=0;i<actual.denominators.size();++i)
        same=same&&pullback::equal(actual.denominators[i],expected.denominators[i]);
      if(same && !graph_.definition && graph_.family_name=="henn_double_pentagon_x0") {
        std::vector<Rational> anchors;bool standard=true;
        for(const auto& node:graph_.nodes) {anchors.push_back(node.anchor);standard=standard&&node.left==0&&node.right==1;}
        return causal::current_example(graph_.family_name,graph_.nodes.size(),anchors,standard);
      }
      if(same && causal::euclidean(*example))return causal::euclidean_family(*example,graph_.nodes.size());
    }
    return causal::derived(graph_.nodes.size());
  }
  // The real singular points of a level's transported system, with the side
  // of each.
  const std::vector<SingularPoint>& singular_points(std::size_t depth,Plan& prepared)const {
    if(prepared.real_singularities)return *prepared.real_singularities;
    std::vector<B> complex;auto* others=prescription_.per_point()?&complex:nullptr;
    auto roots=real_singular_points(prepared.gauge.matrix,xi_,ei_,others);
    for(auto& root:real_singular_points(prepared.beta_rows,xi_,ei_,others))roots.push_back(std::move(root));
    prepared.complex_singularities.clear();
    for(const auto& root:complex) {
      struct Real {arb_t value;Real(){arb_init(value);}~Real(){arb_clear(value);}} magnitude;
      arb_abs(magnitude.value,acb_imagref(root.raw()));
      ComplexPoint point{rational_bound(acb_realref(root.raw()),false),rational_bound(acb_realref(root.raw()),true),rational_bound(magnitude.value,false)};
      if(point.high<Rational(0) || Rational(1)<point.low)continue;
      prepared.complex_singularities.push_back(std::move(point));
    }
    std::vector<std::pair<Rational,Rational>> intervals;
    for(const auto& root:roots)intervals.emplace_back(rational_bound(acb_realref(root.raw()),false),rational_bound(acb_realref(root.raw()),true));
    std::sort(intervals.begin(),intervals.end(),[](const auto& a,const auto& b){return a.first<b.first;});
    // Enclosures of one root found through two denominators overlap.
    std::vector<SingularPoint> merged;
    for(auto& interval:intervals) {
      if(!merged.empty() && !(merged.back().high<interval.first)) {
        if(merged.back().high<interval.second)merged.back().high=interval.second;
      } else merged.push_back({std::move(interval.first),std::move(interval.second),0,{}});
    }
    if(prescription_.derive_orientations) {
      auto& landau=landau_.at(depth);std::string unavailable;
      if(!landau)try {landau.emplace(graph_.nodes[depth].merged,graph_.field.variables()[xi_]);}
        catch(const std::exception& error) {unavailable=std::string("its Symanzik polynomials could not be formed: ")+error.what();}
      for(auto& point:merged) {
        if(!(Rational(0)<point.low) || !(point.high<Rational(1)))continue;
        if(!landau){point.undetermined=unavailable;continue;}
        const auto verdict=landau->orient(point.low,point.high,prescription_.f_rim);
        if(verdict.kind==causal::Orientation::Kind::Oriented)point.side=verdict.side;
        else if(verdict.kind==causal::Orientation::Kind::Undetermined)point.undetermined=verdict.reason.empty()?"undecided":verdict.reason;
      }
    } else if(!prescription_.points.empty()) {
      // A stated level side binds every point not stated on its own.
      for(auto& point:merged)point.side=prescription_.levels.at(depth).x_detour_sign;
    }
    const Rational tolerance("1/1000000");
    for(const auto& stated:prescription_.points) {
      if(stated.level!=depth+1)continue;
      bool found=false;
      for(auto& point:merged)if(!(stated.x<point.low-tolerance) && !(point.high+tolerance<stated.x)) {
        point.side=stated.side;point.undetermined.clear();found=true;
      }
      if(!found) {
        std::string known;
        for(const auto& point:merged)if(Rational(0)<point.low && point.high<Rational(1)) {
          B middle;arb_set_fmpq(acb_realref(middle.raw()),((point.low+point.high)/Rational(2)).raw(),128);
          known+=(known.empty()?"":", ")+middle.real_midpoint(12);
        }
        throw std::invalid_argument("causal.points: level "+std::to_string(depth+1)+" has no real singular point at x = "+stated.x.str()+
          " (its points in (0,1): "+(known.empty()?"none":known)+")");
      }
    }
    prepared.real_singularities=std::move(merged);
    return *prepared.real_singularities;
  }
  [[noreturn]] void undetermined(std::size_t depth,const SingularPoint& point)const {
    B middle;arb_set_fmpq(acb_realref(middle.raw()),((point.low+point.high)/Rational(2)).raw(),128);
    throw std::domain_error("the causal side of the real singular point x = "+middle.real_midpoint(20)+" of level "+
      std::to_string(depth+1)+" could not be derived from the level's Symanzik F ("+point.undetermined+
      "); supply an explicit `causal` object that states its side, \"causal\": {\"points\": [{\"level\": "+std::to_string(depth+1)+
      ", \"x\": \""+middle.real_midpoint(12)+"\", \"side\": 1}]}, with side 1 (above the real axis) or -1");
  }
  void report(std::size_t depth,const std::string& phase,int high)const {if(options_.progress)options_.progress(depth,phase,high);}
  static int checked(long value) {
    if(value< -1000 || value>1000)throw std::overflow_error("recursive epsilon index budget");return static_cast<int>(value);
  }
  static Expansion constant_expansion(const Matrix& matrix) {
    Expansion out{static_cast<unsigned>(matrix.size()),static_cast<unsigned>(matrix[0].size()),{}};
    for(unsigned i=0;i<matrix.size();++i)for(unsigned j=0;j<matrix[i].size();++j)
      if(!matrix[i][j].is_zero())out.terms.push_back({i,j,0,Rational(0),Rational(0),matrix[i][j]});
    return out;
  }
  Matrix reverse(Matrix matrix,bool connection=false)const {
    const auto& sample=graph_.dimension;auto point=exact_point(sample,xi_,sample.constant(1)-sample.variable(xi_));
    for(auto& row:matrix)for(auto& c:row){c=c.substitute(point);if(connection)c=-c;}return matrix;
  }
  std::optional<Endpoint> direct_endpoint(const Matrix& matrix,const Matrix& diagonal,
      const fuchsify::Result& fuchs,const Node& node,bool upper)const {
    if(!options_.direct_endpoint_functionals || options_.endpoint_order>256)
      return std::nullopt;
    // This stronger guard covers both integral rows and endpoint-limit rows:
    // the latter have no fixed Laurent/log sector and their formal DR constant
    // is exactly zero. It does not set the physical nonanalytic modes to zero.
    if(!direct_adjoint_endpoint::only_epsilon_dependent_exponents(fuchs.matrix,xi_,ei_))
      return std::nullopt;
    const auto d=matrix.size();const auto x=graph_.dimension.variable(xi_);
    auto weights=fuchsify::detail::zeros(node.operations.size(),d,graph_.dimension);
    Endpoint result;
    result.functionals.resize(node.operations.size());
    const auto collect=[&](const Matrix& coefficients){for(const auto& row:coefficients)
      for(const auto& c:row)if(!c.is_zero())result.clearance_coefficients.push_back(c);};
    const auto physical_gauge=fuchsify::detail::multiply(diagonal,fuchs.transform);
    collect(matrix);collect(fuchs.matrix);collect(fuchs.transform);collect(fuchs.inverse_transform);collect(physical_gauge);
    for(unsigned i=0;i<node.operations.size();++i) {
      const auto& operation=node.operations[i];
      if(operation.operation==feynman::Operation::Direct ||
          (upper && operation.operation==feynman::Operation::LowerLimit) ||
          (!upper && operation.operation==feynman::Operation::UpperLimit))continue;
      Matrix row{node.observable_rows.at(i)};
      if(operation.operation==feynman::Operation::BetaIntegral)for(auto& c:row[0])
        c=c*x.constant(operation.normalization)*x.pow(operation.left_power)*(x.constant(1)-x).pow(operation.right_power);
      if(upper)row=reverse(std::move(row));
      collect(row);
      const auto transformed=fuchsify::detail::multiply(row,physical_gauge);collect(transformed);
      if(operation.operation==feynman::Operation::BetaIntegral)weights[i]=transformed[0];
    }
    auto primitive=direct_adjoint_endpoint::prepare(fuchs.matrix,weights,xi_,ei_,options_.endpoint_order);
    if(!primitive.success())return std::nullopt;
    // q acts on the Fuchsian state; transport acts on the epsilon-gauged state.
    const auto q=direct_adjoint_endpoint::polynomial(primitive,xi_);
    if(options_.exact_functional_reduction) {
      try {result.direct_operator_exact=direct_adjoint_endpoint::exact_residual_zero(q,fuchs.matrix,weights,xi_,ei_);}
      catch(const std::length_error&) {result.direct_operator_exact=false;}
    }
    result.direct_operator=fuchsify::detail::multiply(q,fuchs.inverse_transform);
    collect(*result.direct_operator);
    return result;
  }
  Endpoint endpoint(const Matrix& matrix,const Matrix& diagonal,const Matrix& inverse_diagonal,const Node& node,bool upper,std::size_t depth,int high) {
    auto fuchs=[&]{Seconds timer{statistics_.endpoint_gauge_seconds};return fuchsify::prepare(matrix,xi_,options_.fuchsification);}();
    if(!fuchs.success)throw std::domain_error("recursive endpoint fuchsification unsupported: "+fuchs.reason);
    const std::string side=upper?"upper":"lower";
    // The endpoint that carries the boundary needs its fundamental series.
    const bool boundary=node.boundary==(upper?BoundarySite::upper:BoundarySite::lower);
    if(auto result=boundary?std::optional<Endpoint>{}:direct_endpoint(matrix,diagonal,fuchs,node,upper)) {
      ++statistics_.direct_endpoints;
      report(depth,side+" direct endpoint functionals; no fundamental series",high);
      return std::move(*result);
    }
    report(depth,side+" exact endpoint series",high);
    auto endpoint_options=options_.endpoint;
    endpoint_options.regulate_noninteger_powers=graph_.specialised;
    const auto caller_progress=endpoint_options.column_progress;
    endpoint_options.column_progress=[&,caller_progress](unsigned completed,unsigned total){
      if(caller_progress)caller_progress(completed,total);
      report(depth,side+" exact endpoint columns "+std::to_string(completed)+"/"+std::to_string(total),high);
    };
    auto series=[&] {
      Seconds timer{statistics_.endpoint_series_seconds};++statistics_.endpoint_series_built;
      return AffineFrobeniusSeries::prepare(fuchs.matrix,xi_,ei_,options_.endpoint_order,endpoint_options);
    }();
    std::optional<Seconds> phase;phase.emplace(statistics_.endpoint_frame_seconds);
    auto physical_gauge=fuchsify::detail::multiply(diagonal,fuchs.transform);
    // P*(D*T*F)^-1*D = P*(T*F)^-1: D stays in the physical observable P and
    // is cancelled from the receiving frame before epsilon normalization and
    // numerical inversion.
    auto matching_frame=series.project(fuchs.transform);
    if(!matching_frame.wronskian_prefactor){
      matching_frame.wronskian_prefactor=series.gauged_wronskian_prefactor(fuchs.transform,fuchs.matrix,matrix);
      if(matching_frame.wronskian_prefactor)report(depth,side+" Wronskian certified from exact trace and rational factors",high);
    }
    Endpoint result{std::move(series),std::move(matching_frame),std::vector<std::optional<Expansion>>(node.operations.size()),{}};
    if(boundary) {
      result.physical_frame=result.series->project(physical_gauge);
      result.physical_inverse_pole_order=fuchsify::detail::pole_order(fuchs.inverse_transform,xi_)+fuchsify::detail::pole_order(inverse_diagonal,xi_);
    }
    const auto collect=[&](const Matrix& coefficients){for(const auto& row:coefficients)for(const auto& c:row)if(!c.is_zero())result.clearance_coefficients.push_back(c);};
    collect(matrix);collect(fuchs.matrix);collect(fuchs.transform);collect(fuchs.inverse_transform);collect(physical_gauge);
    const auto x=graph_.dimension.variable(xi_);
    phase.emplace(statistics_.endpoint_observable_seconds);
    for(std::size_t i=0;i<node.operations.size();++i) {
      const auto& op=node.operations[i];
      if(op.operation==feynman::Operation::Direct ||
          (upper && op.operation==feynman::Operation::LowerLimit) ||
          (!upper && op.operation==feynman::Operation::UpperLimit))continue;
      Matrix row{node.observable_rows.at(i)};
      if(op.operation==feynman::Operation::BetaIntegral)for(auto& c:row[0])
        c=c*x.constant(op.normalization)*x.pow(op.left_power)*(x.constant(1)-x).pow(op.right_power);
      if(upper)row=reverse(std::move(row));
      collect(row);
      auto transformed=fuchsify::detail::multiply(row,physical_gauge);collect(transformed);
      try {
      const bool integral=op.operation==feynman::Operation::BetaIntegral;
      auto domain=integral?result.series->dr_domain(result.series->project(transformed),true):
        result.series->project_endpoint_domain(transformed);
      for(const auto& constraint:domain.zero_constraints) {
        if(result.pending_constraints.size()>=options_.max_endpoint_constraints)
          throw std::length_error("recursive endpoint constraint budget exhausted");
        result.pending_constraints.push_back({constant_expansion(Matrix{constraint.coefficients}),
          "level "+std::to_string(depth+1)+" "+side+" operation "+std::to_string(i)+
          " fixed power "+constraint.power.str()+" log degree "+std::to_string(constraint.log_degree)});
      }
      if(integral)result.functionals[i]=result.series->dr_integral_from_zero(domain.admissible);
      else result.functionals[i]=constant_expansion(result.series->dr_endpoint_constant(domain.admissible));
      }catch(const std::length_error& error) {
        throw std::length_error(side+" endpoint operation "+std::to_string(i)+": "+error.what());
      }catch(const std::domain_error& error) {
        // E.g. the retained x-order guard: name the endpoint and operation it concerns.
        throw std::domain_error(side+" endpoint operation "+std::to_string(i)+": "+error.what());
      }
      if((i+1)%32==0 || i+1==node.operations.size())
        report(depth,side+" endpoint observables "+std::to_string(i+1)+"/"+std::to_string(node.operations.size()),high);
    }
    return result;
  }
  // The SAME endpoint functionals in other coordinates. Y=D_base U_base=D U
  // gives U_base=V U with the exact V=D_base^-1 D, hence E=E_base V at the
  // ordinary cut and beta=beta_base V on the contour. Nothing about the
  // endpoint frames, their domain constraints or the overlap changes, so no
  // endpoint series is rebuilt. Returns null if the coordinates coincide.
  std::unique_ptr<Plan> pullback_basis_plan(Plan& base,const PhysicalBasisTransform& basis) {
    auto candidate=std::make_unique<Plan>();
    candidate->gauge=epsilon_diagonal_gauge(basis.connection,ei_);
    const auto d=candidate->gauge.matrix.size();
    candidate->diagonal=basis.from_basis;
    candidate->inverse_diagonal=basis.to_basis;
    for(std::size_t i=0;i<d;++i)for(std::size_t j=0;j<d;++j) {
      candidate->diagonal[i][j]=epsilon_gauge_detail::multiply_power(
          candidate->diagonal[i][j],ei_,candidate->gauge.shifts[j]);
      candidate->inverse_diagonal[i][j]=epsilon_gauge_detail::multiply_power(
          candidate->inverse_diagonal[i][j],ei_,-candidate->gauge.shifts[i]);
    }
    if(candidate->inverse_diagonal==base.inverse_diagonal)return {};
    candidate->pullback_map=fuchsify::detail::multiply(base.inverse_diagonal,candidate->diagonal);
    candidate->geometry=base.geometry;
    candidate->lower.functionals.resize(base.lower.functionals.size());
    candidate->upper.functionals.resize(base.upper.functionals.size());
    candidate->constraint_labels=base.constraint_labels;
    candidate->beta_indices=base.beta_indices;
    if(!base.beta_rows.empty())candidate->beta_rows=fuchsify::detail::multiply(base.beta_rows,candidate->pullback_map);
    candidate->has_physical_basis=true;
    candidate->boundary_map=base.boundary_map;
    for(const auto& row:candidate->inverse_diagonal)for(const auto& entry:row)if(!entry.is_zero())
      candidate->largest_shift=std::max(candidate->largest_shift,-*exact_epsilon_valuation(entry,ei_));
    candidate->lower_operator=pulled_back_operator(base,candidate->pullback_map,false,0);
    candidate->upper_operator=pulled_back_operator(base,candidate->pullback_map,true,0);
    return candidate;
  }
  // E_base V at the cut, through epsilon order high. V's exact pole depth at
  // the cut decides how far the base operator must be known, and the base
  // operator's low how far V must be expanded.
  LaurentRows pulled_back_operator(Plan& base,const Matrix& map,bool upper,int high) {
    const auto h=graph_.dimension.constant(base.geometry.overlap(upper));
    const auto at=upper?graph_.dimension.constant(1)-h:h;
    const auto& known=upper?base.upper_operator:base.lower_operator;
    // An operator's low is structural, so a retained one predicts the V window.
    auto inner=exact_laurent_rows(map,at,std::max(0,checked(static_cast<long>(high)-(known?known->low:0))));
    const auto& outer=cached_endpoint_operator(base,upper,std::max(0,checked(static_cast<long>(high)-inner.low)));
    const auto inner_high=std::max(0,checked(static_cast<long>(high)-outer.low));
    if(inner.high<inner_high)inner=exact_laurent_rows(map,at,inner_high);
    return linear_boundary::compose(outer,inner,high);
  }
  Plan& plan(std::size_t depth,unsigned desired_top) {
    if(plans_[depth])return *plans_[depth];const auto& node=graph_.nodes[depth];
    const auto explicit_basis=options_.physical_bases.find(depth);
    const bool supplied=explicit_basis!=options_.physical_bases.end();
    auto selected=build_plan(depth,desired_top,supplied?&explicit_basis->second:nullptr);
    if(options_.automatic_basis_scan && (!supplied || options_.scan_physical_bases)) {
      // Starts after the original plan's own preflight, which every execution
      // needs. Operator work done for candidates is selection time.
      std::optional<Seconds> timer;double operators_before=0;
      BasisSelection selection;selection.depth=depth;selection.masters=selected->gauge.matrix.size();
      const auto characters=[](const Matrix& matrix) {std::size_t size=0;for(const auto& row:matrix)for(const auto& value:row)if(!value.is_zero())size+=value.str().size();return size;};
      selection.original_characters=characters(selected->gauge.matrix);
      std::optional<artifacts::Identity> basis_identity;
      std::optional<basis_persistence::Choice> cached_basis;
      std::optional<PhysicalBasisTransform> accepted_basis;
      if(basis_store_ && !supplied) {
        basis_identity=basis_persistence::identity(graph_,depth,xi_,ei_,prescription_);
        // Corruption or a failed exact identity is fatal, not a fresh scan that
        // silently changes coordinates underneath an existing continuation.
        cached_basis=basis_persistence::load(*basis_store_,*basis_identity,node.closure.matrix,xi_);
        if(cached_basis){selection.cache_hit=true;++statistics_.basis_cache_hits;}
      }
      const auto began=std::chrono::steady_clock::now();
      // Only a decision that another run would reach again is durable: an
      // installed choice, or a discovery that ran to its deterministic end.
      // A wall-clock cut or an error says nothing about the coordinates.
      bool durable=false;
      try {
        preflight_linear_demands(depth,*selected);
        timer.emplace(statistics_.basis_selection_seconds);operators_before=statistics_.endpoint_operator_seconds;
        selection.original_child_loss=*selected->factored_child_loss;
        const auto cost=[&](const Plan& value) {
          // The public child coordinates are pinned, so this child lookahead
          // propagates into every remaining level. Preserve the per-component
          // ledger as a secondary score instead of adding unrelated pole maxima.
          std::uint64_t demand=0;
          for(const auto& h:value.component_child_loss)if(h)demand+=std::max(0,*h);
          return physical_basis_scan::DemandCost{
            static_cast<std::uint64_t>(*value.factored_child_loss),demand,
            static_cast<std::uint64_t>(characters(value.gauge.matrix))};
        };
        // Keep this original plan fixed during discovery. Cached endpoint
        // jets may grow, but neither its coordinates nor the ranking baseline
        // changes until the complete cheap frontier has been collected.
        const auto original_cost=cost(*selected);
        auto best_cost=original_cost;
        const auto describe_cost=[](const physical_basis_scan::DemandCost& value) {
          return "child="+std::to_string(value.at(0))+", component-sum="+
            std::to_string(value.at(1))+", characters="+std::to_string(value.at(2));
        };
        // Installing a basis costs one exact pullback of the retained endpoint
        // functionals. Its demand is the plan's own, not a prediction.
        const auto install=[&](const PhysicalBasisTransform& basis,bool persisted=false) -> bool {
          ++selection.full_candidates_checked;
          if(persisted)++selection.candidates_checked;
          try {
            // Discovery builds candidates by exact updates without checking
            // them; the one that is adopted is verified in full here.
            verify_physical_basis(basis,node.closure.matrix,xi_);
            auto candidate=pullback_basis_plan(*selected,basis);
            if(!candidate) {
              // The persisted coordinates are those of the original gauge.
              if(persisted){selection.accepted=true;selection.candidate_characters=selection.original_characters;
                selection.candidate_child_loss=selection.original_child_loss;accepted_basis=basis;}
              return persisted;
            }
            preflight_linear_demands(depth,*candidate,false);
            const auto candidate_cost=cost(*candidate);
            if(!persisted && !(candidate_cost<best_cost)) {
              report(depth,"basis candidate rejected: "+describe_cost(candidate_cost)+" does not improve "+describe_cost(best_cost),desired_top);return false;
            }
            selection.candidate_characters=candidate_cost.at(2);
            selection.candidate_child_loss=*candidate->factored_child_loss;
            accepted_basis=basis;best_cost=candidate_cost;
            candidate->pullback_base=std::move(selected);selected=std::move(candidate);selection.accepted=true;
            report(depth,std::string(persisted?"persisted basis installed: ":"basis candidate installed: ")+
              describe_cost(best_cost)+", screened="+std::to_string(selection.candidates_screened)+
              "; endpoint functionals pulled back exactly, no endpoint series rebuilt",desired_top);
            return true;
          } catch(const std::exception& error) {
            if(persisted)throw;
            selection.reason=std::string("candidate unresolved: ")+error.what();
            report(depth,selection.reason,desired_top);return false;
          }
        };
        if(cached_basis) {
          // A persisted selected coordinate system is a continuation contract,
          // not a candidate to rank again under today's approximation scores.
          report(depth,"rechecking persisted exact master basis",desired_top);
          if(cached_basis->selected)(void)install(cached_basis->basis,true);
          if(!selection.accepted)selection.reason="persisted original basis retained";
        } else {
          report(depth,"scanning simpler master bases",desired_top);
          ++statistics_.basis_scans;
          physical_basis_scan::Options scan_options;
          scan_options.seconds=options_.basis_scan_seconds;scan_options.epsilon_variable=ei_;
          scan_options.max_candidates=8;
          // Discovery ranks candidates by the size of their exact connection
          // alone: that is what transport pays for, and it needs no endpoint
          // data. Whether a candidate keeps the child's epsilon window is
          // decided once, on the plan pulled back for the best of them.
          // From a supplied basis the scan sees the level in those coordinates
          // (requested integrals included, cached identities not), and each of
          // its candidates is composed with the basis before installation.
          level::Result in_basis;
          if(supplied) {
            const auto& basis=explicit_basis->second;
            in_basis.ordered_basis=node.closure.ordered_basis;in_basis.requested_integrals=node.closure.requested_integrals;
            in_basis.matrix=basis.connection;
            if(!node.closure.target_rows.empty())in_basis.target_rows=fuchsify::detail::multiply(node.closure.target_rows,basis.from_basis);
          }
          const auto& scanned=supplied?in_basis:node.closure;
          const auto composed=[&](const PhysicalBasisTransform& change) {
            if(!supplied)return change;
            const auto& basis=explicit_basis->second;
            return PhysicalBasisTransform{basis.original_connection,fuchsify::detail::multiply(change.to_basis,basis.to_basis),
              fuchsify::detail::multiply(basis.from_basis,change.from_basis),change.connection};
          };
          const std::uint64_t original_size=physical_basis_scan::expression_characters(scanned.matrix);
          auto frontier=physical_basis_scan::scan_candidates(scanned,node.merged.physical_count,xi_,scan_options,
            [&](const physical_basis_scan::Candidate& candidate) -> std::optional<physical_basis_scan::DemandCost> {
              if(candidate.metadata.provenance!="original"){++selection.candidates_screened;++selection.candidates_checked;}
              return physical_basis_scan::DemandCost{candidate.metadata.connection_characters};
            });
          selection.scan_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
          report(depth,"basis discovery complete: screened="+std::to_string(selection.candidates_screened)+
            ", frontier="+std::to_string(frontier.candidates.size())+", original "+describe_cost(original_cost),desired_top);
          for(const auto& candidate:frontier.candidates) {
            if(candidate.demand_cost->at(0)>=original_size || selection.full_candidates_checked==8)break;
            if(install(composed(candidate.transform)))break;
          }
          durable=selection.accepted || !frontier.stats.time_expired;
          if(!durable)report(depth,"basis discovery cut by its time budget; the undecided original basis is not persisted",desired_top);
        }
        if(selection.accepted) {
          ++statistics_.basis_selections;
          selection.reason="smaller transported connection at no wider child epsilon window; endpoint functionals pulled back exactly";
        }
        else if(selection.reason.empty())selection.reason="no smaller connection within the bounded search";
      }catch(const std::exception& error) {
        if(cached_basis && cached_basis->selected)throw;
        if(!cached_basis && !selection.scan_seconds)selection.scan_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
        selection.reason=std::string("keeping original basis: ")+error.what();
        durable=selection.accepted;
      }
      if(timer)statistics_.endpoint_operator_seconds=operators_before;
      if(basis_store_ && !supplied && !cached_basis && durable) {
        auto unit=fuchsify::detail::identity(node.closure.matrix.size(),graph_.dimension);
        basis_persistence::Choice choice{accepted_basis?*accepted_basis:
          PhysicalBasisTransform{node.closure.matrix,unit,unit,node.closure.matrix},selection.accepted};
        basis_persistence::save(*basis_store_,*basis_identity,choice,node.closure.matrix,xi_);
        ++statistics_.basis_cache_writes;
      }
      selection.total_selection_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
      report(depth,std::string(selection.accepted?"simpler basis selected: ":"original basis retained: ")+selection.reason,desired_top);
      basis_selections_.push_back(std::move(selection));
      if(options_.basis_selection_observer)options_.basis_selection_observer(basis_selections_.back());
    }
    plans_[depth]=std::move(selected);++statistics_.exact_plans;return *plans_[depth];
  }
  std::unique_ptr<Plan> build_plan(std::size_t depth,unsigned desired_top,const PhysicalBasisTransform* chosen) {
    const auto& node=graph_.nodes[depth];
    if(node.observable_rows.size()!=node.operations.size())throw std::invalid_argument("recursive observable plan shape");
    if(chosen)verify_physical_basis(*chosen,node.closure.matrix,xi_);
    auto gauge=epsilon_diagonal_gauge(chosen?chosen->connection:node.closure.matrix,ei_);const auto d=gauge.matrix.size();
    auto diagonal=fuchsify::detail::identity(d,graph_.dimension),inverse=diagonal;
    for(std::size_t i=0;i<d;++i) {
      diagonal[i][i]=graph_.dimension.variable(ei_).pow(gauge.shifts[i]);
      inverse[i][i]=graph_.dimension.constant(1)/diagonal[i][i];
    }
    if(chosen) {
      exact_jet_composition::Limits fusion_limits;fusion_limits.max_term_work=20000000;
      diagonal=exact_jet_composition::fuse(chosen->from_basis,diagonal,ei_,fusion_limits).product;
      inverse=exact_jet_composition::fuse(inverse,chosen->to_basis,ei_,fusion_limits).product;
    }
    auto lower=endpoint(gauge.matrix,diagonal,inverse,node,false,depth,desired_top),
      upper=endpoint(reverse(gauge.matrix,true),reverse(diagonal),reverse(inverse),node,true,depth,desired_top);
    const auto constraints=lower.pending_constraints.size()+upper.pending_constraints.size();
    if(constraints>options_.max_endpoint_constraints || node.operations.size()+constraints>5000)
      throw std::length_error("recursive combined endpoint constraint budget exhausted");
    std::vector<std::string> labels;
    lower.functionals.resize(node.operations.size()+constraints);upper.functionals.resize(node.operations.size()+constraints);
    for(auto* endpoint:{&lower,&upper})for(auto& [functional,label]:endpoint->pending_constraints) {
      endpoint->functionals[node.operations.size()+labels.size()]=std::move(functional);labels.push_back(std::move(label));
    }
    lower.pending_constraints.clear();upper.pending_constraints.clear();
    Matrix beta;std::vector<std::size_t> indices;auto x=graph_.dimension.variable(xi_);
    for(std::size_t i=0;i<node.operations.size();++i)if(node.operations[i].operation==feynman::Operation::BetaIntegral) {
      auto row=node.observable_rows[i];const auto& op=node.operations[i];
      for(auto& c:row)c=c*x.constant(op.normalization)*x.pow(op.left_power)*(x.constant(1)-x).pow(op.right_power);
      beta.push_back(std::move(row));indices.push_back(i);
    }
    if(!beta.empty())beta=fuchsify::detail::multiply(beta,diagonal);
    auto geometry=endpoint_clearance(lower,upper);
    lower.clearance_coefficients.clear();upper.clearance_coefficients.clear();
    std::int64_t largest=0;
    for(const auto& row:inverse)for(const auto& entry:row)if(!entry.is_zero())
      largest=std::max(largest,-*exact_epsilon_valuation(entry,ei_));
    auto result=std::make_unique<Plan>(Plan{std::move(gauge),std::move(diagonal),std::move(inverse),std::move(lower),std::move(upper),std::move(beta),std::move(indices),largest,std::move(geometry)});
    result->constraint_labels=std::move(labels);
    result->has_physical_basis=chosen!=nullptr;
    if(node.boundary!=BoundarySite::anchor) {
      const auto& carrier=node.boundary==BoundarySite::upper?result->upper:result->lower;
      const auto determination=endpoint_boundary::determine(*carrier.series,*carrier.physical_frame,ei_,carrier.physical_inverse_pole_order);
      // The exact preparation chose this endpoint on the same test.
      if(!determination.determined)throw std::domain_error("level "+std::to_string(depth+1)+
        ": the endpoint values do not fix the level's solution: "+determination.reason);
      // At a numeric dimension an exponent that depends on epsilon is told
      // from an integer by its value: a dimension at which one of them is an
      // integer would be taken for an analytic solution.
      if(node.analytic_solutions && determination.analytic_dimension!=node.analytic_solutions)
        throw std::domain_error("level "+std::to_string(depth+1)+": "+std::to_string(determination.analytic_dimension)+
          " local solutions are analytic at the endpoint in this dimension, "+std::to_string(node.analytic_solutions)+" in general");
      result->boundary_map=endpoint_boundary::taylor_map(*carrier.physical_frame,determination,
        result->geometry.overlap(node.boundary==BoundarySite::upper));
      report(depth,std::string(node.boundary==BoundarySite::upper?"upper":"lower")+" endpoint boundary: "+
        std::to_string(determination.analytic_dimension)+" analytic local solutions of "+std::to_string(d)+
        " fixed by the masters at the endpoint; epsilon order of the map from "+std::to_string(boundary_map_low(*result)),desired_top);
    }
    // Exact endpoint functionals already encode the prescribed singular limit.
    // At ordinary cuts their rational maps can participate in a jointly
    // certified rewrite. Numerical Frobenius operators never enter this path.
    if(options_.exact_functional_reduction && result->lower.direct_operator_exact && result->upper.direct_operator_exact) {
      namespace fr=ft_functional_reduction;
      Matrix weight(result->lower.functionals.size(),std::vector<Exact>(d,graph_.dimension.constant(0)));
      for(std::size_t i=0;i<result->beta_indices.size();++i)weight[result->beta_indices[i]]=result->beta_rows[i];
      const auto h=graph_.dimension.constant(result->geometry.lower_overlap);
      const auto upper_h=graph_.dimension.constant(result->geometry.upper_overlap);
      fr::Functional functional{result->gauge.matrix,
        fr::detail::at(*result->lower.direct_operator,xi_,h),weight,
        fr::detail::at(*result->upper.direct_operator,xi_,upper_h),
        h,graph_.dimension.constant(1)-upper_h,xi_,ei_};
      if(!fr::epsilon_regular(functional)) {
        ++statistics_.functional_reduction_attempts;
        fr::SearchOptions limits;limits.seconds=1;limits.pole_order=1;limits.spatial_degree=2;
        for(const auto* matrix:{&functional.lower,&functional.weight,&functional.upper})
          for(const auto& row:*matrix)for(const auto& value:row)if(auto v=exact_epsilon_valuation(value,ei_);v && *v<0)
            limits.pole_order=std::max(limits.pole_order,static_cast<unsigned>(std::min<std::int64_t>(32,-*v)));
        auto search=fr::search_counterterm(functional,graph_.dimension.constant(1),limits);
        if(search.certificate) {
          fr::verify_counterterm(functional,*search.certificate);
          result->rewritten_lower=search.certificate->rewritten.lower;
          result->rewritten_upper=search.certificate->rewritten.upper;
          result->beta_rows.clear();result->beta_indices.clear();
          for(std::size_t i=0;i<weight.size();++i) {
            const auto& row=search.certificate->rewritten.weight[i];
            if(std::any_of(row.begin(),row.end(),[](const Exact& x){return !x.is_zero();})) {
              result->beta_indices.push_back(i);result->beta_rows.push_back(row);
            }
          }
          ++statistics_.functional_reduction_certificates;
          report(depth,"exact endpoint and integral epsilon poles cancelled together",desired_top);
        }
      }
    }
    if(constraints)report(depth,"physical endpoint zero constraints "+std::to_string(constraints),desired_top);
    report(depth,"endpoint overlap="+result->geometry.lower_overlap.str()+" (lower), "+result->geometry.upper_overlap.str()+" (upper)",desired_top);
    return result;
  }
  EndpointGeometry endpoint_clearance(const Endpoint& lower,const Endpoint& upper)const {
    using M=NativeTailMagnitude;
    EndpointGeometry result{Rational("1/2"),Rational("1/2"),{},M::zero(),M::zero(),M::zero()};
    bool found[2]={false,false};
    const auto& sample=graph_.dimension;auto epsilon=sample.variable(ei_),x=sample.variable(xi_);
    const auto& names=sample.variables();auto ii=std::find(names.begin(),names.end(),"I");
    // A polynomial may occur at both endpoints, in their two coordinates.
    std::set<std::pair<bool,std::string>> seen;
    for(const auto* endpoint:{&lower,&upper})for(const auto& coefficient:endpoint->clearance_coefficients) {
      const bool above=endpoint==&upper;
      // Keep the denominator before specializing: a numerator cancellation at
      // epsilon=0 must not hide poles present in higher epsilon coefficients.
      auto p=coefficient.denominator();
      auto ep=*exact_epsilon_valuation(p,ei_);
      if(ep)p=p/epsilon.pow(ep);
      p=p.substitute(exact_point(sample,ei_,sample.constant(0)));
      if(p.is_zero())throw std::domain_error("endpoint denominator has no epsilon-regular specialization");
      if(ii!=names.end())p=polynomial_norm(p,ii-names.begin(),sample.constant(-1));
      if(p.is_zero())throw std::domain_error("endpoint denominator vanishes on the specified complex sheet");
      auto origin=fuchsify::detail::valuation(p,xi_);if(origin)p=p/x.pow(origin);
      if(p.is_rational())continue;
      p=p/p.constant(p.numerator_terms()[0].coefficient);
      if(!seen.insert({above,p.str()}).second)continue;
      auto roots=polynomial_roots(p,xi_,128);
      for(auto& root:roots) {
        auto bound=M::lower_abs(root);
        if(bound.is_zero() || !bound.is_finite())throw std::domain_error("nonzero endpoint pole cannot be separated from zero at the root-isolation precision");
        if(result.nonzero_poles.empty() || result.nearest_pole_lower>bound)result.nearest_pole_lower=bound;
        auto& nearest=above?result.upper_nearest_pole:result.lower_nearest_pole;
        if(!found[above] || nearest>bound)nearest=bound;
        found[above]=true;
        result.nonzero_poles.push_back(std::move(root));
      }
    }
    // This proves a geometric convergence margin for the formal epsilon
    // coefficients. It does not certify the omitted affine Taylor tail.
    // Each end separately, from the distance to its own nearest singular point.
    const unsigned long ratio=options_.overlap_ratio;
    for(const bool above:{false,true}) {
      auto& overlap=above?result.upper_overlap:result.lower_overlap;
      const auto& nearest=above?result.upper_nearest_pole:result.lower_nearest_pole;
      bool settled=false;
      for(unsigned halvings=0;halvings<options_.max_overlap_halvings && !settled;++halvings) {
        const auto size=M::upper_abs(B::from_strings(overlap.str()));
        if(overlap<=options_.overlap && overlap<=Rational(1)/Rational(static_cast<long>(ratio)) &&
            (!found[above] || size*M::from_ui(ratio)<=nearest))settled=true;
        else overlap=overlap/Rational(2);
      }
      if(!settled)throw std::runtime_error("recursive geometric endpoint-overlap budget exhausted");
    }
    return result;
  }
  void check_endpoint_constraint(const LaurentBoundary& value,unsigned row,const std::string& label) {
    const auto tolerance=NativeTailMagnitude::lower_abs(B::from_strings(options_.endpoint_constraint_tolerance.str()));
    for(unsigned k=0;k<value.values.at(row).size();++k) {
      const auto& coefficient=value.values[row][k];
      auto bound=NativeTailMagnitude::upper_abs(coefficient);
      if(!coefficient.is_finite() || !bound.is_finite() || bound>tolerance)
        throw std::domain_error("physical endpoint constraint failed: "+label+" at epsilon "+
          std::to_string(value.low+static_cast<int>(k))+"; midpoint=("+
          coefficient.real_midpoint(16)+","+coefficient.imag_midpoint(16)+"); radius2exp=("+
          coefficient.real_radius_exponent()+","+coefficient.imag_radius_exponent()+
          "); contains_zero="+(coefficient.contains_zero()?"true":"false")+
          "; absolute tolerance="+options_.endpoint_constraint_tolerance.str());
      statistics_.maximum_endpoint_constraint_residual=std::max(statistics_.maximum_endpoint_constraint_residual,bound.approximate_upper());
      ++statistics_.endpoint_constraint_coefficients;
    }
    ++statistics_.endpoint_constraint_rows;
  }
  void discharge_endpoint_constraints(std::size_t depth,const Plan& prepared,
      linear_boundary::Expression& expression,LaurentBoundary& result) {
    const auto rows=graph_.nodes[depth].operations.size();
    if(result.values.size()!=rows+prepared.constraint_labels.size() ||
        expression.transform.coefficients.size()!=result.values.size())
      throw std::logic_error("recursive constrained expression shape mismatch");
    for(unsigned i=0;i<prepared.constraint_labels.size();++i)
      check_endpoint_constraint(result,rows+i,prepared.constraint_labels[i]);
    // Auxiliary domain rows are never exposed as requested integrals or passed
    // to a parent. The shared-source contraction happens before this split.
    expression.transform.coefficients.resize(rows);result.values.resize(rows);
  }
  // Composition measures the outer window against the child's stored map, a
  // materialized value against its own trimmed low. A child map produced by
  // another method (factored/component packing) can retain exactly-zero
  // leading orders, so its declared low may lie below the value low.
  int child_window_low(std::size_t child_depth,const LaurentBoundary& child) const {
    if(child_depth>=expressions_.size() || !expressions_[child_depth] || !expressions_[child_depth]->leaf_source)return child.low;
    const auto& expression=*expressions_[child_depth];
    return std::min(child.low,checked(static_cast<long>(expression.transform.low)+expression.leaf_source->low));
  }
  LaurentBoundary compose_child_expression(std::size_t depth,const LaurentRows& local,unsigned desired_top,const Plan* prepared=nullptr) {
    if(depth+1>=expressions_.size() || !expressions_[depth+1] || !cache_[depth+1])
      throw std::logic_error("recursive adjoint child has no synchronized linear expression");
    const auto& child=*expressions_[depth+1];
    const int child_required=checked(static_cast<long>(desired_top)-local.low);
    if(cache_[depth+1]->high()<child_required)throw BoundaryDemand(child_required,"recursive expression requires additional child epsilon coefficients");
    const int high=checked(static_cast<long>(desired_top)-child.leaf_source->low);
    linear_boundary::Expression expression;
    try {expression=linear_boundary::compose(local,child,high);}
    catch(const linear_boundary::CompositionDemand& demand) {
      if(demand.required_outer_high>local.high)throw std::logic_error("recursive local operator lookahead invariant failed");
      throw BoundaryDemand(checked(static_cast<long>(demand.required_inner_high)+child.leaf_source->low),"recursive expression requires a higher child operator window");
    }
    // Materialization drops whole-map exact-zero leading orders from the
    // value. Drop them from the stored map as well, so a parent that sizes its
    // operator from this value's low agrees with the map it composes.
    linear_boundary::tighten(expression.transform);
    LaurentBoundary result;
    try {result=linear_boundary::materialize(expression,desired_top);}
    catch(const BoundaryDemand& demand) {
      throw BoundaryDemand(checked(static_cast<long>(demand.required_high)+child.transform.low),"recursive expression requires a higher shared leaf window");
    }
    result.taylor_tail_certified=false;
    if(prepared)discharge_endpoint_constraints(depth,*prepared,expression,result);
    observe_operator(depth,"composed shared-source map",expression.transform);
    expressions_[depth]=std::move(expression);cache_[depth]=std::move(result);return *cache_[depth];
  }
  LaurentRows endpoint_operator(Endpoint& endpoint,const Rational& exact_point_value,int high) {
    Seconds timer(statistics_.endpoint_operator_seconds);
    if(endpoint.direct_operator) {
      auto direct=exact_laurent_rows(*endpoint.direct_operator,graph_.dimension.constant(exact_point_value),high);
      // The direct primitive has one row per operation. Zero-constraint rows
      // contributed by the OTHER endpoint's series have no functional here:
      // they are exact zero rows, as on the series route below.
      if(direct.coefficients.size()<endpoint.functionals.size())
        direct.coefficients.resize(endpoint.functionals.size(),
          std::vector(direct.columns(),std::vector<B>(direct.high-direct.low+1,B(0))));
      return direct;
    }
    const auto point=B::from_strings(exact_point_value.str());
    const auto rows=endpoint.functionals.size();const auto d=endpoint.matching_frame.columns;
    if(std::none_of(endpoint.functionals.begin(),endpoint.functionals.end(),[](const auto& value){return value.has_value();}))
      return {0,high,std::vector(rows,std::vector(d,std::vector<B>(high+1,B(0))))};
    const int inverse_high=constant_demand(endpoint,high);
    if(endpoint.inverse && endpoint.inverse->success() && endpoint.inverse_point==exact_point_value)
      affine_operator::extend(*endpoint.inverse,*endpoint.series,endpoint.matching_frame,inverse_high,options_.matching);
    else {
      endpoint.inverse=affine_operator::prepare(*endpoint.series,endpoint.matching_frame,point,inverse_high,options_.matching);
      endpoint.inverse_point=exact_point_value;
    }
    const auto& inverse=*endpoint.inverse;
    if(!inverse.success())throw std::domain_error("recursive endpoint operator inverse unsupported: "+inverse.reason);
    std::vector<std::optional<affine_operator::Operator>> operators(rows);int low=0;
    for(unsigned i=0;i<rows;++i)if(endpoint.functionals[i]) {
      operators[i]=affine_operator::compose(inverse,*endpoint.series,*endpoint.functionals[i],point,high,options_.matching);
      if(!operators[i]->success())throw std::domain_error("recursive endpoint operator composition unsupported: "+operators[i]->reason);
      low=std::min(low,operators[i]->matrix.low);
    }
    LaurentRows result{low,high,std::vector(rows,std::vector(d,std::vector<B>(high-low+1,B(0))))};
    for(unsigned i=0;i<rows;++i)if(operators[i]) {
      const auto& matrix=operators[i]->matrix;
      for(unsigned j=0;j<d;++j)for(int k=matrix.low;k<=high;++k)result.coefficients[i][j][k-low]=matrix.coefficients[0][j][k-matrix.low];
    }
    return result;
  }
  // This is the actual leaf-source width prescribed by the immutable graph:
  // scalar leaves return one row per request, and exact-zero nodes use one
  // constant source. It avoids materializing a low-order child just to select
  // a method; an already prepared expression remains authoritative.
  std::size_t shared_source_width(std::size_t depth)const {
    for(;depth<graph_.nodes.size();++depth) {
      if(expressions_[depth])return expressions_[depth]->transform.columns();
      const auto& node=graph_.nodes[depth];
      if(node.scalar_leaf)return node.requested.size();
      if(node.closure.ordered_basis.empty())return 1;
    }
    throw std::logic_error("recursive linear source has no terminal node");
  }
  void observe_operator(std::size_t depth,const std::string& phase,const LaurentRows& rows)const {
    if(options_.operator_observer)options_.operator_observer(depth,phase,rows);
  }
  const LaurentRows& cached_endpoint_operator(Plan& prepared,bool upper,int high) {
    auto& saved=upper?prepared.upper_operator:prepared.lower_operator;
    if(!saved || saved->high<high) {
      const auto& rewritten=upper?prepared.rewritten_upper:prepared.rewritten_lower;
      if(prepared.pullback_base)saved=pulled_back_operator(*prepared.pullback_base,prepared.pullback_map,upper,high);
      else if(!prepared.pullback_map.empty())throw std::logic_error("a screening plan cannot extend its endpoint operators");
      else if(rewritten)saved=exact_laurent_rows(*rewritten,graph_.dimension.constant(prepared.geometry.overlap(upper)),high);
      else saved=endpoint_operator(upper?prepared.upper:prepared.lower,prepared.geometry.overlap(upper),high);
    }
    return *saved;
  }
  // Exact rows applied to the child at the anchor: their pole depth and the
  // child-component highs they need for an output at order zero. They depend
  // on the graph node only, never on the internal basis of a plan.
  const DirectDemand& direct_demand(std::size_t depth,const Plan* prepared=nullptr) {
    auto& saved=direct_demands_[depth];if(saved)return *saved;
    const auto& node=graph_.nodes[depth];DirectDemand result;
    if(node.boundary!=BoundarySite::anchor && !prepared)throw std::logic_error("direct rows of an endpoint level need its plan");
    const Rational point=prepared?matching_point(depth,*prepared):node.anchor;
    result.needs.assign(node.closure.ordered_basis.size(),std::nullopt);
    for(std::size_t i=0;i<node.operations.size();++i)if(node.operations[i].operation==feynman::Operation::Direct)
      for(std::size_t j=0;j<node.observable_rows[i].size();++j) {
        const auto& entry=node.observable_rows[i][j];if(entry.is_zero())continue;
        const auto at=entry.substitute(exact_point(entry,xi_,entry.constant(point)));
        if(at.is_zero())continue;
        const int valuation=checked(*exact_epsilon_valuation(at,ei_));
        result.low=std::min(result.low,valuation);
        epsilon_demand::raise(result.needs.at(j),epsilon_demand::checked(-static_cast<long>(valuation)));
      }
    saved=std::move(result);return *saved;
  }
  // For every output at epsilon order zero: the order through which each child
  // component is consumed, from the endpoint rows retained through
  // operator_high. Exact valuations of the connection, the forcing rows and
  // the anchor gauge propagate the windows; a zero-containing ball of nonzero
  // radius counts as nonzero. Equivalently, column j of the local operator
  // vanishes identically below order -needs[j].
  epsilon_demand::Highs component_child_needs(std::size_t depth,Plan& prepared,int operator_high) {
    namespace ed=epsilon_demand;
    const auto d=prepared.gauge.matrix.size();
    const ed::Highs outputs(prepared.lower.functionals.size(),0);
    if(!prepared.connection_valuations)prepared.connection_valuations=ed::valuations(prepared.gauge.matrix,ei_);
    if(!prepared.anchor_inverse_valuations) {
      auto inverse=prepared.inverse_diagonal;const auto& node=graph_.nodes[depth];
      for(auto& row:inverse)for(auto& value:row)if(!value.is_zero())value=value.substitute(exact_point(value,xi_,value.constant(node.anchor)));
      prepared.anchor_inverse_valuations=ed::valuations(inverse,ei_);
    }
    auto needs=ed::merge(ed::pullback(cached_endpoint_operator(prepared,false,operator_high),outputs),
                         ed::pullback(cached_endpoint_operator(prepared,true,operator_high),outputs));
    if(needs.size()!=d)throw std::logic_error("endpoint operator width differs from the connection");
    for(const auto& row:prepared.beta_rows)for(std::size_t j=0;j<d;++j)if(!row[j].is_zero())
      ed::raise(needs[j],ed::checked(-static_cast<long>(checked(*exact_epsilon_valuation(row[j],ei_)))));
    // Both arms and the forcing reach the anchor through the same connection.
    const auto initial=ed::close(*prepared.connection_valuations,std::move(needs));
    return ed::merge(ed::pullback(*prepared.anchor_inverse_valuations,initial),direct_demand(depth).needs);
  }
  void preflight_linear_demands(std::size_t depth,Plan& prepared,bool announce=true) {
    if(prepared.adjoint_child_loss)return;
    if(announce)report(depth,"structural child-window preflight",0);
    const auto l=cached_endpoint_operator(prepared,false,0).low;
    const auto r=cached_endpoint_operator(prepared,true,0).low;
    int q=0;
    for(const auto& row:prepared.beta_rows)for(const auto& entry:row)if(!entry.is_zero())
      q=std::min(q,checked(*exact_epsilon_valuation(entry,ei_)));
    const int direct=direct_demand(depth,&prepared).low;
    const long loss=std::max<std::int64_t>(0,prepared.largest_shift);
    // Epsilon-regular homogeneous transport preserves lower bounds, while
    // forcing can lower them to q; the inverse D gauge follows. This is the
    // bound of one common window for all child components.
    const int uniform=checked(std::max(loss+std::max({-static_cast<long>(l),-static_cast<long>(r),-static_cast<long>(q)}),-static_cast<long>(direct)));
    int tight=uniform;epsilon_demand::Highs component;
    if(graph_.nodes[depth].boundary!=BoundarySite::anchor) {
      // The child enters through the boundary map: one common window, deeper
      // by the map's own pole order.
      tight=checked(static_cast<long>(uniform)-boundary_map_low(prepared));
      prepared.component_child_loss.clear();
      prepared.factored_child_loss=tight;prepared.adjoint_child_loss=tight;
      ++statistics_.demand_preflights;return;
    }
    if(options_.component_epsilon_demands) {
      component=component_child_needs(depth,prepared,0);
      tight=std::min(tight,std::max(0,epsilon_demand::maximum(component)));
    }
    // Publish only a complete profile: a failed propagation must fail again on
    // the next call instead of leaving the common-window half behind.
    prepared.component_child_loss=std::move(component);
    prepared.factored_child_loss=tight;prepared.adjoint_child_loss=uniform;
    ++statistics_.demand_preflights;
  }
  // Replace rounding noise below each child component's structural window by
  // the exact zero it represents, then drop the orders that became empty. A
  // column shift or an exact gauge product declares one common low for all
  // columns; without this the child would be asked for orders that no
  // coefficient consumes. A value excluding zero would contradict the
  // structural bound: it is kept, and the wider child demand with it.
  void impose_structural_columns(std::size_t depth,Plan& prepared,LaurentRows& local,int operator_high) {
    const auto needs=component_child_needs(depth,prepared,operator_high);
    for(std::size_t j=0;j<needs.size();++j)if(needs[j])
      for(int k=local.low;k<=local.high && k< -*needs[j];++k)for(auto& row:local.coefficients) {
        auto& value=row[j][k-local.low];if(value.is_zero())continue;
        if(value.contains_zero()){value=B(0);++statistics_.structural_zero_coefficients;}
        else ++statistics_.structural_zero_conflicts;
      }
    linear_boundary::tighten(local);
  }
  PreparedAdjointStage adjoint_stage(std::size_t depth,Plan& prepared,int high) {
    const auto& node=graph_.nodes[depth];const auto sample=graph_.dimension;
    report(depth,"lower observable operator",high);auto left=cached_endpoint_operator(prepared,false,high);
    observe_operator(depth,"lower endpoint operator",left);
    report(depth,"upper observable operator",high);auto right=cached_endpoint_operator(prepared,true,high);
    observe_operator(depth,"upper endpoint operator",right);
    Matrix forcing(prepared.lower.functionals.size(),std::vector<Exact>(prepared.gauge.matrix.size(),sample.constant(0)));
    for(unsigned i=0;i<prepared.beta_indices.size();++i)forcing[prepared.beta_indices[i]]=prepared.beta_rows[i];
    auto upper_forcing=forcing;for(auto& row:upper_forcing)for(auto& entry:row)entry=-entry;
    (void)node;
    const auto anchor=sample.constant(matching_point(depth,prepared)),h=sample.constant(prepared.geometry.lower_overlap),
      end=sample.constant(Rational(1)-prepared.geometry.upper_overlap);
    return {prepared.gauge.matrix,std::move(forcing),std::move(upper_forcing),
      std::move(left),std::move(right),path(h,anchor,depth,prepared),path(end,anchor,depth,prepared),prepared.gauge.shifts};
  }
  LaurentBoundary evaluate_adjoint(std::size_t depth,unsigned desired_top,Plan& prepared) {
    const auto& node=graph_.nodes[depth];auto sample=graph_.dimension;
    const auto anchor=sample.constant(matching_point(depth,prepared));
    const bool at_endpoint=node.boundary!=BoundarySite::anchor;
    preflight_linear_demands(depth,prepared);
    // With component demands the local operator is tightened to its structural
    // columns below, so the child is needed through the tight bound only. The
    // refinement loop still decides: a lower retained low raises the demand.
    int demand=checked(static_cast<long>(desired_top)+*prepared.factored_child_loss);
    auto& combined=prepared.local_operator;
    for(unsigned refinement=0;refinement<options_.max_refinements;++refinement) {
      if(demand>static_cast<int>(options_.max_epsilon))throw std::runtime_error("recursive adjoint child epsilon demand exceeds finite budget");
      report(depth,"child boundary",demand);auto child=evaluate_child(depth+1,demand);
      const int needed_operator=checked(static_cast<long>(desired_top)-child_window_low(depth+1,child));
      if(!combined || combined->high<needed_operator) {
        // Rows on the masters at the matching point; at an endpoint level the
        // boundary map follows, whose poles ask for more of them.
        const int local_high=checked(static_cast<long>(needed_operator)-boundary_map_low(prepared));
        const int high=checked(local_high+prepared.largest_shift);
        auto stage=adjoint_stage(depth,prepared,high);
        auto settings=options_.adjoint;settings.taylor_order=options_.ordinary_order;
        // What a leg cost and why: its size, path and chart counts.
        const auto summarized=[&](const char* side,LaurentRows input,const Matrix& forcing,const std::vector<Exact>& vertices) {
          const auto* stats=settings.conditioning_stats;
          const AdjointConditioningStats before=stats?*stats:AdjointConditioningStats{};
          // The arm of the endpoint that carries the boundary starts at the
          // matching point.
          if(vertices.size()<2)return input;
          const auto rows=input.coefficients.size(),columns=input.coefficients.empty()?0:input.coefficients[0].size();
          const int low=input.low,top=input.high;
          const auto began=std::chrono::steady_clock::now();
          auto result=ordinary_transport(stage.connection,std::move(input),forcing,vertices,settings);
          std::ostringstream text;
          text<<side<<" adjoint transport done: "<<rows<<" rows x "<<columns<<" masters, epsilon "<<low<<".."<<top
              <<", "<<vertices.size()-1<<" path legs";
          if(stats)text<<", "<<stats->charts-before.charts<<" charts ("<<stats->centered_charts-before.centered_charts
              <<" centred, "<<stats->conditioning_subdivisions-before.conditioning_subdivisions<<" halved)";
          text<<", "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<" s";
          if(stats)text<<" (compile "<<stats->compilation_seconds-before.compilation_seconds<<", prepare "
              <<stats->preparation_seconds-before.preparation_seconds<<", recurrence "<<stats->recurrence_seconds-before.recurrence_seconds
              <<"; "<<stats->circuit_terms-before.circuit_terms<<" products)";
          report(depth,text.str(),high);
          return result;
        };
        report(depth,"lower adjoint transport",high);auto left=summarized("lower",std::move(stage.lower_endpoint),stage.lower_forcing,stage.lower_path);
        observe_operator(depth,"lower transported operator",left);
        report(depth,"upper adjoint transport",high);auto right=summarized("upper",std::move(stage.upper_endpoint),stage.upper_forcing,stage.upper_path);
        observe_operator(depth,"upper transported operator",right);
        auto transported=add_laurent_rows(left,right);
        if(prepared.has_physical_basis) {
          const auto inverse=exact_laurent_rows(prepared.inverse_diagonal,anchor,
            std::max(0,checked(static_cast<long>(local_high)-transported.low)));
          combined=linear_boundary::compose(transported,inverse,local_high);
        } else combined=shift_laurent_columns(transported,prepared.gauge.shifts,true);
        Matrix direct(prepared.lower.functionals.size(),std::vector<Exact>(prepared.gauge.matrix.size(),sample.constant(0)));
        for(unsigned i=0;i<node.operations.size();++i)if(node.operations[i].operation==feynman::Operation::Direct)direct[i]=node.observable_rows[i];
        combined=add_laurent_rows(*combined,exact_laurent_rows(direct,anchor,combined->high));
        if(at_endpoint) {
          const auto map=exact_laurent_rows(*prepared.boundary_map,anchor,
            std::max(0,checked(static_cast<long>(needed_operator)-combined->low)));
          combined=linear_boundary::compose(*combined,map,needed_operator);
          linear_boundary::tighten(*combined);
        } else
        if(options_.component_epsilon_demands)impose_structural_columns(depth,prepared,*combined,high);
        observe_operator(depth,"combined local operator",*combined);
      } else ++statistics_.local_operator_reuses;
      try {
        report(depth,"shared boundary application",desired_top);
        return compose_child_expression(depth,*combined,desired_top,&prepared);
      }catch(const BoundaryDemand& request) {
        if(request.required_high<=demand)throw std::runtime_error("recursive adjoint epsilon refinement made no progress");
        demand=request.required_high;++statistics_.refinements;report(depth,"child refinement",demand);
      }
    }
    throw std::runtime_error("recursive adjoint child epsilon refinement budget exhausted");
  }
  LaurentBoundary evaluate_factored(std::size_t depth,unsigned desired_top,Plan& prepared) {
    const auto& node=graph_.nodes[depth];const auto sample=graph_.dimension;
    const auto anchor=sample.constant(node.anchor),h=sample.constant(prepared.geometry.lower_overlap),end=sample.constant(Rational(1)-prepared.geometry.upper_overlap);
    Matrix beta(prepared.lower.functionals.size(),std::vector<Exact>(prepared.gauge.matrix.size(),sample.constant(0)));
    for(unsigned i=0;i<prepared.beta_indices.size();++i)beta[prepared.beta_indices[i]]=prepared.beta_rows[i];
    int beta_low=0;
    for(const auto& row:beta)for(const auto& entry:row)if(!entry.is_zero())
      beta_low=std::min(beta_low,checked(*exact_epsilon_valuation(entry,ei_)));
    preflight_linear_demands(depth,prepared);
    int demand=checked(static_cast<long>(desired_top)+*prepared.factored_child_loss);
    for(unsigned refinement=0;refinement<options_.max_refinements;++refinement) {
      if(demand>static_cast<int>(options_.max_epsilon))throw std::runtime_error("recursive factored child epsilon demand exceeds finite budget");
      report(depth,"child boundary",demand);(void)evaluate_child(depth+1,demand);
      if(!expressions_[depth+1])throw std::logic_error("recursive factored child has no shared expression");
      const auto& child=*expressions_[depth+1];const int leaf_low=child.leaf_source->low;
      const int high=checked(static_cast<long>(desired_top)-leaf_low);
      const int gauge_low=checked(static_cast<long>(child.transform.low)-std::max<std::int64_t>(0,prepared.largest_shift));
      const int endpoint_high=std::max(0,checked(static_cast<long>(high)-gauge_low));
      report(depth,"lower observable operator",endpoint_high);
      const auto& left=cached_endpoint_operator(prepared,false,endpoint_high);
      report(depth,"upper observable operator",endpoint_high);
      const auto& right=cached_endpoint_operator(prepared,true,endpoint_high);
      observe_operator(depth,"lower endpoint operator",left);
      observe_operator(depth,"upper endpoint operator",right);
      try {
        if(options_.component_epsilon_demands) {
          namespace ed=epsilon_demand;
          const ed::Highs outputs(beta.size(),high);
          const auto left_needs=ed::pullback(left,outputs),right_needs=ed::pullback(right,outputs);
          const auto middle_needs=ed::close(prepared.gauge.matrix,
            ed::merge(right_needs,ed::pullback(beta,outputs,ei_)),ei_);
          const auto initial_needs=ed::close(prepared.gauge.matrix,ed::merge(left_needs,middle_needs),ei_);
          auto inverse=exact_laurent_rows(prepared.inverse_diagonal,anchor,
            std::max(0,checked(static_cast<long>(ed::maximum(initial_needs))-child.transform.low)));
          ed::Expression gauged;
          try {gauged=ed::compose(inverse,ed::rectangular(child),initial_needs);}
          catch(const linear_boundary::CompositionDemand& request) {
            throw BoundaryDemand(checked(static_cast<long>(request.required_inner_high)+leaf_low),"component anchor map requires additional child coefficients");
          }
          factored_transport::Options settings;settings.transport=options_.adjoint;settings.transport.taylor_order=options_.ordinary_order;
          settings.transport_dispatch=[this](const Matrix& a,LaurentRows input,const Matrix& b,const std::vector<Exact>& vertices,const AdjointOptions& options) {
            return ordinary_transport(a,std::move(input),b,vertices,options);
          };
          report(depth,"lower factored transport",ed::maximum(initial_needs));
          auto lower=factored_transport::evolve_components(prepared.gauge.matrix,gauged,{},path(anchor,h,depth,prepared),initial_needs,{},settings);
          observe_operator(depth,"lower factored map",lower.physical.value.transform);
          report(depth,"middle factored transport",ed::maximum(middle_needs));
          auto middle=factored_transport::evolve_components(prepared.gauge.matrix,lower.physical,beta,path(h,end,depth,prepared),right_needs,outputs,settings);
          observe_operator(depth,"upper factored map",middle.physical.value.transform);
          observe_operator(depth,"middle integrated map",middle.integrated.value.transform);
          statistics_.component_transports+=2;
          statistics_.component_slots+=lower.retained_slots+middle.retained_slots;
          statistics_.component_rectangular_slots+=lower.rectangular_slots+middle.rectangular_slots;
          auto lower_integral=ed::compose(left,lower.physical,high);
          auto upper_integral=ed::compose(right,middle.physical,high);
          auto combined=add_laurent_rows(add_laurent_rows(lower_integral.transform,upper_integral.transform),middle.integrated.value.transform);
          Matrix direct(prepared.lower.functionals.size(),std::vector<Exact>(prepared.gauge.matrix.size(),sample.constant(0)));
          for(std::size_t i=0;i<node.operations.size();++i)if(node.operations[i].operation==feynman::Operation::Direct)direct[i]=node.observable_rows[i];
          auto direct_rows=exact_laurent_rows(direct,anchor,std::max(0,checked(static_cast<long>(high)-child.transform.low)));
          auto direct_expression=ed::compose(direct_rows,ed::rectangular(child),high);
          combined=add_laurent_rows(combined,direct_expression.transform);
          linear_boundary::Expression expression{std::move(combined),child.leaf_source};
          observe_operator(depth,"composed shared-source map",expression.transform);
          LaurentBoundary value;
          try {value=linear_boundary::materialize(expression,desired_top);}
          catch(const BoundaryDemand& request) {throw BoundaryDemand(checked(static_cast<long>(request.required_high)+child.transform.low),"component result needs additional shared leaf coefficients");}
          value.taylor_tail_certified=false;discharge_endpoint_constraints(depth,prepared,expression,value);
          expressions_[depth]=std::move(expression);cache_[depth]=std::move(value);return *cache_[depth];
        }
        const int middle_high=std::max(high,checked(static_cast<long>(high)-right.low));
        const int initial_high=std::max({checked(static_cast<long>(high)-left.low),middle_high,checked(static_cast<long>(high)-beta_low)});
        auto inverse_gauge=exact_laurent_rows(prepared.inverse_diagonal,anchor,
            std::max(0,checked(static_cast<long>(initial_high)-child.transform.low)));
        linear_boundary::Expression gauged;
        try {gauged=linear_boundary::compose(inverse_gauge,child,initial_high);}
        catch(const linear_boundary::CompositionDemand& request) {
          if(request.required_outer_high>inverse_gauge.high)throw std::logic_error("factored anchor gauge lookahead invariant failed");
          throw BoundaryDemand(checked(static_cast<long>(request.required_inner_high)+leaf_low),"factored anchor map requires additional child coefficients");
        }
        auto settings=factored_transport::Options{};settings.transport=options_.adjoint;settings.transport.taylor_order=options_.ordinary_order;
        settings.transport_dispatch=[this](const Matrix& a,LaurentRows input,const Matrix& b,
            const std::vector<Exact>& vertices,const AdjointOptions& options) {
          return ordinary_transport(a,std::move(input),b,vertices,options);
        };
        report(depth,"lower factored transport",initial_high);
        // A single exact zero accumulator allows reuse of the same homogeneous
        // forward-map backend on the anchor-to-overlap segment.
        auto lower=factored_transport::evolve(prepared.gauge.matrix,gauged,
            Matrix(1,std::vector<Exact>(prepared.gauge.matrix.size(),sample.constant(0))),
            path(anchor,h,depth,prepared),initial_high,settings).physical;
        observe_operator(depth,"lower factored map",lower.transform);
        report(depth,"middle factored transport",middle_high);
        auto middle=factored_transport::evolve(prepared.gauge.matrix,lower,beta,path(h,end,depth,prepared),middle_high,high,settings);
        observe_operator(depth,"upper factored map",middle.physical.transform);
        observe_operator(depth,"middle integrated map",middle.integrated.transform);
        auto lower_integral=linear_boundary::compose(left,lower,high);
        auto upper_integral=linear_boundary::compose(right,middle.physical,high);
        auto combined=add_laurent_rows(add_laurent_rows(lower_integral.transform,upper_integral.transform),middle.integrated.transform);
        Matrix direct(prepared.lower.functionals.size(),std::vector<Exact>(prepared.gauge.matrix.size(),sample.constant(0)));
        for(unsigned i=0;i<node.operations.size();++i)if(node.operations[i].operation==feynman::Operation::Direct)direct[i]=node.observable_rows[i];
        auto direct_rows=exact_laurent_rows(direct,anchor,std::max(0,checked(static_cast<long>(high)-child.transform.low)));
        linear_boundary::Expression direct_expression;
        try {direct_expression=linear_boundary::compose(direct_rows,child,high);}
        catch(const linear_boundary::CompositionDemand& request) {
          if(request.required_outer_high>direct_rows.high)throw std::logic_error("factored direct operator lookahead invariant failed");
          throw BoundaryDemand(checked(static_cast<long>(request.required_inner_high)+leaf_low),"factored direct map requires additional child coefficients");
        }
        combined=add_laurent_rows(combined,direct_expression.transform);
        observe_operator(depth,"composed shared-source map",combined);
        linear_boundary::Expression expression{std::move(combined),child.leaf_source};
        LaurentBoundary value;
        try {value=linear_boundary::materialize(expression,desired_top);}
        catch(const BoundaryDemand& request) {
          throw BoundaryDemand(checked(static_cast<long>(request.required_high)+child.transform.low),"factored result requires additional shared leaf coefficients");
        }
        value.taylor_tail_certified=false;discharge_endpoint_constraints(depth,prepared,expression,value);
        expressions_[depth]=std::move(expression);cache_[depth]=std::move(value);return *cache_[depth];
      }catch(const factored_transport::MapDemand& request) {
        const int next=checked(static_cast<long>(request.required_high)+std::max<std::int64_t>(0,prepared.largest_shift)+leaf_low);
        if(next<=demand)throw std::runtime_error("recursive factored map refinement made no progress");
        demand=next;++statistics_.refinements;report(depth,"child refinement",demand);
      }catch(const BoundaryDemand& request) {
        if(request.required_high<=demand)throw std::runtime_error("recursive factored epsilon refinement made no progress");
        demand=request.required_high;++statistics_.refinements;report(depth,"child refinement",demand);
      }
    }
    throw std::runtime_error("recursive factored child epsilon refinement budget exhausted");
  }
  // ---- Weighted accumulation: one row per top-level target at every level ----
  //
  // The requested integrals are w0 L0 L1 ... L(n-1) v with the local operators
  // Lk (rows: the integrals level k is asked for, columns: its masters at the
  // matching point, which are the integrals its child is asked for), the
  // values v of the terminal node and w0 the identity on the top-level
  // targets. By rows the product is formed from the right, every Lk in full.
  // Here it is formed from the left: level k receives weights w(k) on its
  // requested integrals and returns w(k+1)=w(k) Lk, for which it transports
  // the single combined row w(k) E (E: its endpoint rows) with the combined
  // forcing w(k) F instead of all rows of E.
  //
  // The weights of a level whose input is not exact are short of an order.
  struct WeightedDemand {std::size_t depth;int deficit;};
  static bool direct_only(const Node& node) {
    return std::all_of(node.operations.begin(),node.operations.end(),
      [](const auto& operation){return operation.operation==feynman::Operation::Direct;});
  }
  // The weights through `high`, with zeros above their known window (the
  // caller accounts for what that leaves valid), followed by one row and one
  // column per endpoint constraint of the level with unit weight: a
  // constraint is a further target that starts at this level.
  static LaurentRows extended_weights(const LaurentRows& weights,std::size_t constraints,int high) {
    const auto rows=weights.coefficients.size(),columns=weights.columns();
    const int low=constraints?std::min(weights.low,0):weights.low;
    high=std::max(high,low);
    LaurentRows out{low,high,std::vector(rows+constraints,std::vector(columns+constraints,std::vector<B>(high-low+1,B(0))))};
    for(std::size_t a=0;a<rows;++a)for(std::size_t i=0;i<columns;++i)
      for(int k=weights.low;k<=std::min(weights.high,high);++k)
        out.coefficients[a][i][k-low]=weights.coefficients[a][i][k-weights.low];
    if(low<=0 && high>=0)for(std::size_t c=0;c<constraints;++c)out.coefficients[rows+c][columns+c][-low]=B(1);
    return out;
  }
  // One level. `weights` (rows: the top-level targets, then the endpoint
  // constraints of the levels above; columns: the level's requested integrals)
  // are known through weights.high, to every order if `exact`. Returns the
  // weights on the level's masters at its matching point through `target`,
  // with the level's own constraint rows appended, or reports by how many
  // orders its input is short.
  LaurentRows weighted_level(std::size_t depth,const LaurentRows& weights,bool exact,int target,
      std::vector<std::string>& labels,unsigned desired_top) {
    const auto& node=graph_.nodes[depth];const auto sample=graph_.dimension;
    const auto requested=node.operations.size();
    if(weights.columns()!=requested)throw std::logic_error("weighted accumulation: the weights do not match the level's requested integrals");
    ++statistics_.numeric_evaluations;
    const long known=weights.high;
    if(direct_only(node)) {
      // No endpoint, no transport: the exact rows at the anchor.
      const auto local=exact_laurent_rows(node.observable_rows,sample.constant(node.anchor),
        std::max(0,checked(static_cast<long>(target)-weights.low)));
      if(!exact && known+local.low<target)throw WeightedDemand{depth,checked(static_cast<long>(target)-known-local.low)};
      auto out=linear_boundary::compose(extended_weights(weights,0,checked(static_cast<long>(target)-local.low)),local,target);
      linear_boundary::tighten(out);
      observe_operator(depth,"weights on the child's integrals",out);
      return out;
    }
    auto& prepared=plan(depth,desired_top);
    preflight_linear_demands(depth,prepared);
    ++statistics_.adjoint_selections;++statistics_.weighted_levels;
    const bool at_endpoint=node.boundary!=BoundarySite::anchor;
    const auto anchor=sample.constant(matching_point(depth,prepared));
    const auto constraints=prepared.constraint_labels.size();
    const auto d=prepared.gauge.matrix.size();
    if(prepared.lower.functionals.size()!=requested+constraints)
      throw std::logic_error("weighted accumulation: endpoint rows differ from the level's operations and constraints");
    // As by rows: rows on the masters at the matching point through
    // local_high (the boundary map of an endpoint level follows and its poles
    // ask for more of them), transported through `high` before the gauge.
    const int map_low=at_endpoint?boundary_map_low(prepared):0;
    const int local_high=checked(static_cast<long>(target)-map_low);
    const int high=checked(local_high+prepared.largest_shift);
    const int weights_low=constraints?std::min(weights.low,0):weights.low;
    const int operator_high=std::max(0,checked(static_cast<long>(high)-weights_low));
    auto stage=adjoint_stage(depth,prepared,operator_high);
    // Every row with a forcing becomes one constant component of the
    // transported state, which holds its weight.
    std::vector<std::size_t> sources;int forcing_low=0;
    for(std::size_t i=0;i<stage.lower_forcing.size();++i) {
      bool forced=false;
      for(const auto& entry:stage.lower_forcing[i])if(!entry.is_zero()) {
        forced=true;forcing_low=std::min(forcing_low,checked(*exact_epsilon_valuation(entry,ei_)));
      }
      if(forced)sources.push_back(i);
    }
    Matrix direct(requested+constraints,std::vector<Exact>(d,sample.constant(0)));
    std::optional<LaurentRows> direct_rows;
    {
      bool any=false;
      for(std::size_t i=0;i<requested;++i)if(node.operations[i].operation==feynman::Operation::Direct){direct[i]=node.observable_rows[i];any=true;}
      if(any)direct_rows=exact_laurent_rows(direct,anchor,std::max(0,checked(static_cast<long>(local_high)-weights_low)));
    }
    // Lowest order of anything a weight is multiplied by before the gauge.
    const int reach=std::min({stage.lower_endpoint.low,stage.upper_endpoint.low,forcing_low});
    long padded=static_cast<long>(high)-reach;
    if(direct_rows)padded=std::max(padded,static_cast<long>(local_high)-direct_rows->low);
    const auto extended=extended_weights(weights,constraints,checked(std::max(0L,padded)));
    const auto rows=extended.coefficients.size(),s=sources.size();
    const Matrix no_forcing(rows,std::vector<Exact>(d+s,sample.constant(0)));
    auto settings=options_.adjoint;settings.taylor_order=options_.ordinary_order;
    // g'=w F-g A for the combined row g=w E is the homogeneous equation of the
    // row (g,u), u=epsilon^q w constant, with the connection [[A,0],[-F/epsilon^q,0]]:
    // q is the lowest epsilon order of the forcing, so that the block keeps a
    // nonnegative epsilon valuation. Uncertainty of the weights is then
    // carried like that of any transported row, and a checkpoint is keyed by
    // the weights through its initial rows.
    const auto arm=[&](const char* side,const LaurentRows& endpoint,const Matrix& forcing,const std::vector<Exact>& vertices) {
      auto initial=linear_boundary::compose(extended,endpoint,high);
      // The arm of the endpoint that carries the boundary starts at the
      // matching point.
      if(vertices.size()<2)return initial;
      const int low=s?std::min(initial.low,checked(static_cast<long>(extended.low)+forcing_low)):initial.low;
      LaurentRows wide{low,high,std::vector(rows,std::vector(d+s,std::vector<B>(high-low+1,B(0))))};
      for(std::size_t a=0;a<rows;++a) {
        for(std::size_t j=0;j<d;++j)for(int k=initial.low;k<=high;++k)
          wide.coefficients[a][j][k-low]=initial.coefficients[a][j][k-initial.low];
        for(std::size_t c=0;c<s;++c)for(int k=low;k<=high;++k) {
          const long order=static_cast<long>(k)-forcing_low;
          if(order>=extended.low && order<=extended.high)
            wide.coefficients[a][d+c][k-low]=extended.coefficients[a][sources[c]][order-extended.low];
        }
      }
      Matrix connection(d+s,std::vector<Exact>(d+s,sample.constant(0)));
      for(std::size_t j=0;j<d;++j)for(std::size_t k=0;k<d;++k)connection[j][k]=stage.connection[j][k];
      const auto scale=sample.variable(ei_).pow(static_cast<unsigned long>(-forcing_low));
      for(std::size_t c=0;c<s;++c)for(std::size_t j=0;j<d;++j) {
        const auto& entry=forcing[sources[c]][j];
        if(!entry.is_zero())connection[d+c][j]=-(entry*scale);
      }
      const auto* stats=settings.conditioning_stats;
      const AdjointConditioningStats before=stats?*stats:AdjointConditioningStats{};
      const auto began=std::chrono::steady_clock::now();
      auto result=ordinary_transport(connection,std::move(wide),no_forcing,vertices,settings);
      std::ostringstream text;
      text<<side<<" adjoint transport done: "<<rows<<" rows x "<<d<<" masters, weighted ("<<requested<<" requested integrals, "
          <<s<<" forcing sources), epsilon "<<low<<".."<<high<<", "<<vertices.size()-1<<" path legs";
      if(stats)text<<", "<<stats->charts-before.charts<<" charts ("<<stats->centered_charts-before.centered_charts
          <<" centred, "<<stats->conditioning_subdivisions-before.conditioning_subdivisions<<" halved)";
      text<<", "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<" s";
      if(stats)text<<" (compile "<<stats->compilation_seconds-before.compilation_seconds<<", prepare "
          <<stats->preparation_seconds-before.preparation_seconds<<", recurrence "<<stats->recurrence_seconds-before.recurrence_seconds
          <<"; "<<stats->circuit_terms-before.circuit_terms<<" products)";
      report(depth,text.str(),high);
      // The source components did their work.
      for(auto& row:result.coefficients)row.resize(d);
      return result;
    };
    report(depth,"lower adjoint transport",high);
    const auto left=arm("lower",stage.lower_endpoint,stage.lower_forcing,stage.lower_path);
    observe_operator(depth,"lower transported weighted rows",left);
    report(depth,"upper adjoint transport",high);
    const auto right=arm("upper",stage.upper_endpoint,stage.upper_forcing,stage.upper_path);
    observe_operator(depth,"upper transported weighted rows",right);
    const auto transported=add_laurent_rows(left,right);
    // `valid`: the order through which the result follows from the declared
    // windows alone, were the weights unknown above `known`.
    long valid=known+reach;
    LaurentRows combined;
    if(prepared.has_physical_basis) {
      const auto inverse=exact_laurent_rows(prepared.inverse_diagonal,anchor,
        std::max(0,checked(static_cast<long>(local_high)-transported.low)));
      combined=linear_boundary::compose(transported,inverse,local_high);
      valid+=inverse.low;
    } else {
      combined=shift_laurent_columns(transported,prepared.gauge.shifts,true);
      valid+=static_cast<long>(combined.high)-transported.high;
    }
    if(direct_rows) {
      combined=add_laurent_rows(combined,linear_boundary::compose(extended,*direct_rows,combined.high));
      valid=std::min(valid,known+direct_rows->low);
    }
    if(at_endpoint) {
      // The child's values enter through the exact boundary map: the weights
      // on the masters at the overlap point are multiplied by it.
      const auto map=exact_laurent_rows(*prepared.boundary_map,anchor,
        std::max(0,checked(static_cast<long>(target)-combined.low)));
      combined=linear_boundary::compose(combined,map,target);
      valid+=map.low;
    } else if(options_.component_epsilon_demands) {
      // Column j of the local operator vanishes identically below -needs[j]
      // (component_child_needs), so column j of the result vanishes below
      // weights_low-needs[j]: rounding noise there is the exact zero. For the
      // same reason weights unknown above `known` cannot reach column j at
      // or below known-needs[j], whatever the declared windows say.
      const auto needs=component_child_needs(depth,prepared,operator_high);
      for(std::size_t j=0;j<needs.size();++j)if(needs[j])
        for(int k=combined.low;k<=combined.high && k<static_cast<long>(weights_low)-*needs[j];++k)for(auto& row:combined.coefficients) {
          auto& value=row[j][k-combined.low];if(value.is_zero())continue;
          if(value.contains_zero()){value=B(0);++statistics_.structural_zero_coefficients;}
          else ++statistics_.structural_zero_conflicts;
        }
      valid=std::max(valid,known-epsilon_demand::maximum(needs,-1000));
    }
    if(!exact && valid<target)throw WeightedDemand{depth,checked(static_cast<long>(target)-valid)};
    if(combined.high<target)throw std::logic_error("weighted accumulation: the local window is below its target");
    linear_boundary::tighten(combined);
    labels.insert(labels.end(),prepared.constraint_labels.begin(),prepared.constraint_labels.end());
    observe_operator(depth,"weights on the child's integrals",combined);
    return combined;
  }
  // The requested integrals of the top level through desired_top, or nothing
  // if the request is not covered: the caller then evaluates by rows.
  // Covered: the adjoint method, a top level that is transported or direct,
  // at most weighted_max_targets requested integrals, and a chain in which
  // every level's child provides exactly its masters down to a scalar leaf or
  // an exact-zero node. Levels with an endpoint boundary, internal bases and
  // endpoint constraints are covered. Only the top level's expression and
  // value are published; a later request for an inner level is evaluated by
  // rows.
  std::optional<LaurentBoundary> evaluate_weighted(unsigned desired_top) {
    const auto& top=graph_.nodes.front();
    if(top.scalar_leaf || top.closure.ordered_basis.empty())return std::nullopt;
    const auto decline=[&](const std::string& reason) -> std::optional<LaurentBoundary> {
      ++statistics_.weighted_fallbacks;
      report(0,"weighted accumulation not applicable ("+reason+"): one row per requested integral",static_cast<int>(desired_top));
      return std::nullopt;
    };
    if(options_.linear_method!=LinearMethod::adjoint)return decline("the linear method is not adjoint");
    const auto targets=top.requested.size();
    if(!targets || targets>options_.weighted_max_targets)
      return decline(std::to_string(targets)+" top-level targets, at most "+std::to_string(options_.weighted_max_targets)+" are weighted");
    std::size_t terminal=0;
    for(;;++terminal) {
      const auto& node=graph_.nodes[terminal];
      if(node.requested.size()!=node.operations.size())return decline("a level's operations differ from its requested integrals");
      if(node.scalar_leaf || node.closure.ordered_basis.empty())break;
      if(terminal+1>=graph_.nodes.size() || graph_.nodes[terminal+1].requested!=node.closure.ordered_basis)
        throw std::logic_error("recursive child does not own the parent's ordered master demand");
    }
    // The structural losses of all levels, as by rows: the same exact plans
    // and preflights, and the terminal node asked for the same first order.
    std::vector<int> loss(terminal,0);long first=desired_top;
    for(std::size_t depth=0;depth<terminal;++depth) {
      const auto& node=graph_.nodes[depth];
      if(direct_only(node))loss[depth]=checked(-static_cast<long>(direct_demand(depth).low));
      else {
        if(!plans_[depth])report(depth,"endpoint preparation",static_cast<int>(desired_top));
        auto& prepared=plan(depth,desired_top);
        preflight_linear_demands(depth,prepared);
        loss[depth]=*prepared.factored_child_loss;
      }
      first+=loss[depth];
      if(first>static_cast<long>(options_.max_epsilon))throw std::runtime_error("recursive adjoint child epsilon demand exceeds finite budget");
    }
    report(0,"weighted accumulation: "+std::to_string(targets)+(targets==1?" row":" rows")+" per level through "+
      std::to_string(terminal)+" levels",static_cast<int>(desired_top));
    report(terminal,"child boundary",static_cast<int>(first));
    (void)evaluate_child(terminal,static_cast<unsigned>(std::max(0L,first)));
    // extra[k]: orders level k lost beyond its structural loss.
    std::vector<int> extra(terminal,0);
    for(unsigned attempt=0;attempt<options_.max_refinements;++attempt) {
      // The weights entering level k are needed through target[k]: those on
      // the terminal values through desired_top minus the lowest order of the
      // values, and every level above by its loss more.
      const int lambda=child_window_low(terminal,*cache_[terminal]);
      std::vector<int> target(terminal+1,0);
      target[terminal]=checked(static_cast<long>(desired_top)-lambda);
      for(std::size_t depth=terminal;depth-->0;)
        target[depth]=checked(static_cast<long>(target[depth+1])+loss[depth]+extra[depth]);
      const int unit_high=std::max(0,target[0]);
      LaurentRows weights{0,unit_high,std::vector(targets,std::vector(targets,std::vector<B>(unit_high+1,B(0))))};
      for(std::size_t i=0;i<targets;++i)weights.coefficients[i][i][0]=B(1);
      std::vector<std::string> labels;
      try {
        for(std::size_t depth=0;depth<terminal;++depth)
          weights=weighted_level(depth,weights,depth==0,target[depth+1],labels,desired_top);
      } catch(const WeightedDemand& request) {
        extra[request.depth]+=request.deficit;++statistics_.refinements;
        report(request.depth,"weighted window refinement",request.deficit);
        continue;
      }
      // The terminal values through the order the weights' poles ask for.
      int request=std::max(0,checked(static_cast<long>(desired_top)-weights.low));
      bool restart=false;
      for(unsigned inner=0;inner<options_.max_refinements;++inner) {
        if(request>static_cast<int>(options_.max_epsilon))throw std::runtime_error("recursive adjoint child epsilon demand exceeds finite budget");
        if(cache_[terminal]->high()<request) {
          ++statistics_.refinements;report(terminal-1,"child refinement",request);
          (void)evaluate_child(terminal,static_cast<unsigned>(request));
        }
        // Deeper poles of the terminal values than the windows were sized for.
        if(child_window_low(terminal,*cache_[terminal])<lambda){restart=true;break;}
        if(!expressions_[terminal])throw std::logic_error("weighted accumulation: the terminal node has no linear expression");
        const auto& child=*expressions_[terminal];
        int next=request;
        try {
          report(0,"shared boundary application",static_cast<int>(desired_top));
          auto expression=linear_boundary::compose(weights,child,checked(static_cast<long>(desired_top)-child.leaf_source->low));
          linear_boundary::tighten(expression.transform);
          auto result=linear_boundary::materialize(expression,static_cast<int>(desired_top));
          result.taylor_tail_certified=false;
          {
            const auto size=[](const B& value){mag_t m;mag_init(m);acb_get_mag(m,value.raw());const double out=mag_is_zero(m)?-std::numeric_limits<double>::infinity():mag_get_d_log2_approx(m);mag_clear(m);return out;};
            const auto& source=*expression.leaf_source;const auto& map=expression.transform;
            for(const auto& row:map.coefficients) {
              double largest=-std::numeric_limits<double>::infinity();
              for(std::size_t column=0;column<row.size();++column)for(std::size_t a=0;a<row[column].size();++a) {
                const double weight=size(row[column][a]);
                if(!std::isfinite(weight))continue;
                for(std::size_t b=0;b<source.values[column].size();++b) {
                  if(static_cast<long>(map.low)+static_cast<long>(a)+source.low+static_cast<long>(b)>static_cast<long>(desired_top))break;
                  largest=std::max(largest,weight+size(source.values[column][b]));
                }
              }
              result.term_log2.push_back(largest);
            }
          }
          if(result.values.size()!=targets+labels.size() || expression.transform.coefficients.size()!=result.values.size())
            throw std::logic_error("weighted constrained expression shape mismatch");
          // A constraint row is the constraint functional of its level applied
          // to that level's solution, contracted with the terminal values like
          // a target; it must vanish. Checked through desired_top.
          for(std::size_t i=0;i<labels.size();++i)check_endpoint_constraint(result,static_cast<unsigned>(targets+i),labels[i]);
          expression.transform.coefficients.resize(targets);result.values.resize(targets);result.term_log2.resize(targets);
          observe_operator(0,"composed shared-source map",expression.transform);
          expressions_[0]=std::move(expression);cache_[0]=std::move(result);
          return *cache_[0];
        } catch(const linear_boundary::CompositionDemand& demand) {
          if(demand.required_outer_high>weights.high)throw std::logic_error("weighted accumulation: the weights are short of the terminal window");
          next=checked(static_cast<long>(demand.required_inner_high)+child.leaf_source->low);
        } catch(const BoundaryDemand& demand) {
          next=checked(static_cast<long>(demand.required_high)+child.transform.low);
        }
        if(next<=request)throw std::runtime_error("recursive weighted epsilon refinement made no progress");
        request=next;
      }
      if(!restart)throw std::runtime_error("recursive weighted terminal epsilon refinement budget exhausted");
      ++statistics_.refinements;report(0,"weighted window refinement",lambda);
    }
    throw std::runtime_error("recursive weighted epsilon refinement budget exhausted");
  }
  int constant_demand(const Endpoint& endpoint,unsigned high)const {
    long result=high;for(const auto& functional:endpoint.functionals)if(functional)
      for(auto demand:endpoint.series->valuation_metadata(*functional).required_source_top(high))result=std::max(result,demand);
    return checked(result);
  }
  // A leg follows the real axis and leaves it only to pass a real singular
  // point of the level, over a triangle on the causal side whose half-width
  // is half the free distance around that point (at most the contour
  // height): passing a pole costs about a dozen charts whatever the size of
  // the detour, and a wide one saves the charts of the approach. The integral is defined on the real segment, so this is the
  // path that crosses no other singularity; and wherever centre and rows are
  // real the Taylor recurrence multiplies real numbers. Vertices are binary
  // fractions of the triangle size, so shifted coefficients stay short.
  std::vector<Exact> path(const Exact& from,const Exact& to,std::size_t depth,Plan& prepared)const {
    if(from==to)return {from};
    const auto& names=graph_.dimension.variables();auto ii=std::find(names.begin(),names.end(),"I");
    if(ii==names.end())throw std::invalid_argument("recursive contour requires exact I variable");
    const int free_side=prescription_.levels.at(depth).x_detour_sign;
    const auto unit=graph_.dimension.variable(ii-names.begin());
    const auto imaginary=[&](int side){return unit*graph_.dimension.constant(Rational(side));};
    // At a numeric dimension the moving poles of the master basis are real
    // points in nearly coincident pairs all along the segment, and the
    // coefficients are long whatever the centre: there the whole leg runs at
    // the contour height, away from them.
    const bool guard=prescription_.per_point();
    if(graph_.specialised && !guard) {
      const auto lift=imaginary(free_side)*graph_.dimension.constant(options_.contour_height);
      return {from,from+lift,to+lift,to};
    }
    const Rational start=from.rational(),finish=to.rational();
    const bool ascending=start<finish;const Rational& low=ascending?start:finish;const Rational& high=ascending?finish:start;
    std::vector<SingularPoint> inside;
    for(const auto& point:singular_points(depth,prepared)) {
      if(!(low<point.high) || !(point.low<high))continue;
      if(!(low<point.low) || !(point.high<high))
        throw std::domain_error("a transport leg ends on a real singular point of its level");
      if(!point.undetermined.empty())undetermined(depth,point);
      inside.push_back(point);
    }
    // No detour may hold a complex singular point of the level between
    // itself and the real axis: the path is then homotopic, in the domain
    // of the system, to the real segment indented at its real singular
    // points. Tested on the box [left,right] x (0,top].
    const auto blocked=[&](const Rational& left,const Rational& right,const Rational& top) {
      if(!guard)return false;
      for(const auto& point:prepared.complex_singularities)
        if(!(point.high<left) && !(right<point.low) && !(top<point.height))return true;
      return false;
    };
    if(graph_.specialised) {
      // The lifted leg stays on one side while the oriented points ask for
      // it; between two that ask for opposite sides it crosses the real
      // axis, vertically, in the middle of the widest gap between the
      // singular points that separate them.
      Rational lifted=options_.contour_height;
      for(unsigned attempt=0;blocked(low,high,lifted+lifted);++attempt) {
        if(attempt==64)throw std::runtime_error("a complex singular point of level "+std::to_string(depth+1)+" lies too close to the real segment for a lifted path");
        lifted=lifted/Rational(2);
      }
      const auto lift=graph_.dimension.constant(lifted);
      int side=0;for(const auto& point:inside)if(point.side){side=point.side;break;}
      if(!side)side=free_side;
      std::vector<Exact> vertices{graph_.dimension.constant(low),graph_.dimension.constant(low)+imaginary(side)*lift};
      for(std::size_t first=0;first<inside.size();++first) {
        if(!inside[first].side || inside[first].side==side)continue;
        std::size_t previous=first;while(previous-->0)if(inside[previous].side)break;
        Rational left=inside[previous].high,right=inside[previous+1].low;
        for(std::size_t k=previous+1;k<first;++k)
          if(right-left<inside[k+1].low-inside[k].high){left=inside[k].high;right=inside[k+1].low;}
        // A short binary fraction in the middle half of the gap.
        const Rational quarter=(right-left)/Rational(4);
        Rational grid(1),crossing(0);
        for(unsigned attempt=0;attempt<4096;++attempt,grid=grid/Rational(2)) {
          const auto ratio=(left+quarter)/grid;fmpz_t whole;fmpz_init(whole);
          fmpz_cdiv_q(whole,fmpq_numref(ratio.raw()),fmpq_denref(ratio.raw()));
          crossing=Rational::from_integer(whole)*grid;fmpz_clear(whole);
          if(!(right-quarter<crossing))break;
        }
        if(!(left<crossing) || !(crossing<right))throw std::runtime_error("no regular point between singular points of opposite causal sides");
        vertices.push_back(graph_.dimension.constant(crossing)+imaginary(side)*lift);
        side=inside[first].side;
        vertices.push_back(graph_.dimension.constant(crossing)+imaginary(side)*lift);
      }
      vertices.push_back(graph_.dimension.constant(high)+imaginary(side)*lift);
      vertices.push_back(graph_.dimension.constant(high));
      if(!ascending)std::reverse(vertices.begin(),vertices.end());
      return vertices;
    }
    std::vector<Exact> vertices{graph_.dimension.constant(low)};
    const auto floor_to=[](const Rational& value,const Rational& grid,bool up) {
      const auto ratio=value/grid;fmpz_t whole;fmpz_init(whole);
      if(up)fmpz_cdiv_q(whole,fmpq_numref(ratio.raw()),fmpq_denref(ratio.raw()));
      else fmpz_fdiv_q(whole,fmpq_numref(ratio.raw()),fmpq_denref(ratio.raw()));
      char* text=fmpz_get_str(nullptr,10,whole);Rational result=Rational{std::string(text)}*grid;
      flint_free(text);fmpz_clear(whole);return result;
    };
    // Singular points closer together than twice the contour height share
    // one detour: passing a pole costs about a dozen charts however small
    // the detour, so a row of poles is cheaper to pass at a distance.
    // Only points that agree on the side share a detour (a free point
    // agrees with either); between two groups of opposite sides the path
    // is back on the real axis, at regular points.
    const Rational height=options_.contour_height,reach=height*Rational(2);
    Rational passed=low;bool squeezed=false;
    for(std::size_t first=0;first<inside.size();) {
      std::size_t last=first;int side=inside[first].side;
      while(last+1<inside.size() && inside[last+1].low-inside[last].high<reach &&
          (!side || !inside[last+1].side || inside[last+1].side==side)) {
        ++last;if(!side)side=inside[last].side;
      }
      Rational size(1),grid(1),left(0),right(0);
      for(;;) {
        if(first==last)side=inside[first].side;
        if(!side)side=free_side;
        const Rational before=inside[first].low-(first?inside[first-1].high:low);
        const Rational after=(last+1<inside.size()?inside[last+1].low:high)-inside[last].high;
        Rational free=(before<after?before:after)/Rational(2);
        if(height<free)free=height;
        size=Rational(1);while(free<size)size=size/Rational(2);
        grid=size/Rational(4);
        left=floor_to(inside[first].low-size,grid,false);right=floor_to(inside[last].high+size,grid,true);
        const auto shrink=[&] {
          size=size/Rational(2);grid=size/Rational(4);
          left=floor_to(inside[first].low-size,grid,false);right=floor_to(inside[last].high+size,grid,true);
        };
        // Next to a separate detour within reach the two must not overlap.
        const bool tight=squeezed || (last+1<inside.size() && inside[last+1].low-inside[last].high<reach);
        if(tight) {
          const Rational limit=last+1<inside.size()?(inside[last].high+inside[last+1].low)/Rational(2):high;
          for(unsigned attempt=0;attempt<4096 && (!(passed<left) || !(right<limit));++attempt)shrink();
          if(!(passed<left) || !(right<limit))throw std::runtime_error("no room for neighbouring detours of opposite causal sides");
        }
        if(guard && first!=last) {
          // A group whose rectangle would hold a complex singular point is
          // given up: its first point is passed alone.
          Rational lift(1);while(height<lift)lift=lift/Rational(2);
          if(blocked(left,right,lift)){last=first;side=0;continue;}
        }
        if(guard && first==last) {
          unsigned attempt=0;
          for(;attempt<200 && blocked(left,right,(right-left)/Rational(2));++attempt)shrink();
          if(attempt==200)throw std::runtime_error("a complex singular point of level "+std::to_string(depth+1)+" lies too close to a real one for a detour");
        }
        break;
      }
      squeezed=last+1<inside.size() && inside[last+1].low-inside[last].high<reach;
      passed=right;
      vertices.push_back(graph_.dimension.constant(left));
      if(first==last) {
        vertices.push_back(graph_.dimension.constant((left+right)/Rational(2))+imaginary(side)*graph_.dimension.constant((right-left)/Rational(2)));
      } else {
        // A rectangle over the group, as high as the contour height allows.
        Rational lift(1);while(height<lift)lift=lift/Rational(2);
        vertices.push_back(graph_.dimension.constant(left)+imaginary(side)*graph_.dimension.constant(lift));
        vertices.push_back(graph_.dimension.constant(right)+imaginary(side)*graph_.dimension.constant(lift));
      }
      vertices.push_back(graph_.dimension.constant(right));
      first=last+1;
    }
    vertices.push_back(graph_.dimension.constant(high));
    if(!ascending)std::reverse(vertices.begin(),vertices.end());
    return vertices;
  }
};
} // namespace diffexp::recursion
