#pragma once
#include "diffexp/recursion_graph.hpp"
#include "diffexp/recursion_leaf.hpp"
#include "diffexp/laurent_transport.hpp"
#include "diffexp/affine_matching.hpp"
#include "diffexp/causal.hpp"
#include "diffexp/adjoint_transport.hpp"
#include "diffexp/affine_operator.hpp"
#include "diffexp/linear_boundary.hpp"
#include "diffexp/cached_affine.hpp"
#include "diffexp/factored_transport.hpp"
#include "diffexp/adjoint_checkpoint.hpp"
#include "diffexp/ft_spectral_checkpoint.hpp"
#include "diffexp/ft_ultraspherical_checkpoint.hpp"
#include "diffexp/direct_adjoint_endpoint.hpp"
#include "diffexp/recursion_physical_basis.hpp"
#include "diffexp/physical_basis_scan.hpp"
#include "diffexp/level_cache.hpp"
#include "diffexp/ft_functional_reduction.hpp"
#include "diffexp/exact_jet_composition.hpp"
#include <map>

namespace diffexp::recursion {
enum class OrdinaryMethod { automatic, spectral, taylor, ultraspherical };
enum class LinearMethod { adjoint, factored, automatic };
struct BasisSelection;
struct NumericalOptions {
  bool component_epsilon_demands=true;
  bool exact_functional_reduction=true;
  bool automatic_basis_scan=true;
  double basis_scan_seconds=30;
  unsigned principal_part_endpoint_milliseconds=5000;
  std::function<void(const BasisSelection&)> basis_selection_observer;
  // Exact, explicitly supplied internal bases, keyed by zero-based FT depth.
  // Every transform is checked against the graph before endpoint preparation.
  std::map<std::size_t,PhysicalBasisTransform> physical_bases;
  slong working_bits=384;
  unsigned endpoint_order=32,ordinary_order=80,max_epsilon=100,max_refinements=8,max_overlap_halvings=128;
  Rational overlap{Rational("1/16")},contour_height{Rational("1/10")};
  unsigned leaf_digits=28;
  // Consistency of retained physical endpoint constraints, not a proof about
  // omitted epsilon/Taylor coefficients or exact vanishing of a period.
  Rational endpoint_constraint_tolerance{Rational("1/100000000000000000000")};
  unsigned max_endpoint_constraints=256;
  bool observable_adjoint=true;
  // Guarded experimental route: build endpoint integration rows directly.
  // Unsupported fixed sectors retain the complete Frobenius construction.
  bool direct_endpoint_functionals=false;
  LinearMethod linear_method=LinearMethod::adjoint;
  AdjointOptions adjoint;
  OrdinaryMethod ordinary_method=OrdinaryMethod::automatic;
  // Zero derives the local spectral goal from leaf_digits + 12.
  ft_spectral::Options spectral=[] {ft_spectral::Options value;value.accuracy_goal=0;return value;}();
  std::optional<causal::Prescription> causal_prescription;
  std::function<void(std::size_t,const std::string&,int)> progress;
  std::function<void(const std::string&)> ordinary_progress;
  // Accepted fraction of the current ordinary arm's original path. This is
  // geometric progress, not an estimate of remaining computing time.
  std::function<void(double)> ordinary_fraction_progress;
  std::function<void(const ExactEpsilonMatrix&,const LaurentRows&,
      const ExactEpsilonMatrix&,const std::vector<Exact>&)> ordinary_problem_observer;
  // Optional read-only diagnostics on completed coefficient maps. This lets
  // callers locate arithmetic loss without materializing the shared source.
  std::function<void(std::size_t,const std::string&,const LaurentRows&)> operator_observer;
  std::filesystem::path endpoint_cache_directory;
  // Empty derives endpoint_cache_directory/physical-bases-v1. This persists
  // exact automatic choices only; explicit physical_bases always take priority.
  std::filesystem::path basis_cache_directory;
  std::filesystem::path ordinary_cache_directory;
  std::size_t ordinary_cache_max_bytes=64*1024*1024;
  std::size_t endpoint_cache_max_bytes=64*1024*1024;
  cached_affine::VerificationLimits endpoint_cache_verification;
  affine_matching::Options matching=[] {affine_matching::Options value;value.max_dimension=256;value.max_epsilon_depth=512;value.numeric_operator_convolution=true;return value;}();
  AffineFrobeniusSeries::Options endpoint=[] {AffineFrobeniusSeries::Options value;value.max_dimension=256;return value;}();
  fuchsify::Options fuchsification=[] {fuchsify::Options value;value.max_dimension=256;return value;}();
  gaussian::Budget gaussian;
};
struct EndpointGeometry {
  Rational overlap;
  std::vector<kernel::ComplexBall> nonzero_poles;
  NativeTailMagnitude nearest_pole_lower;
};
struct NumericalStatistics {
  std::size_t component_transports=0,component_slots=0,component_rectangular_slots=0;
  std::size_t functional_reduction_attempts=0,functional_reduction_certificates=0;
  std::size_t basis_scans=0,basis_selections=0,basis_cache_hits=0,basis_cache_writes=0;
  std::size_t principal_part_attempts=0,principal_part_certificates=0,principal_part_reused=0,principal_part_zero_coefficients=0;
  adjoint_checkpoint::Statistics ordinary_checkpoints;
  std::size_t exact_plans=0,leaf_evaluations=0,numeric_evaluations=0,refinements=0,cache_hits=0;
  std::size_t endpoint_series_built=0,endpoint_series_reused=0;
  std::size_t direct_endpoints=0;
  std::size_t adjoint_selections=0,factored_selections=0;
  std::size_t demand_preflights=0,local_operator_reuses=0;
  std::size_t conditioning_method_fallbacks=0;
  std::size_t endpoint_constraint_rows=0,endpoint_constraint_coefficients=0;
  double maximum_endpoint_constraint_residual=0;
  double ordinary_seconds=0,spectral_seconds=0;
  std::size_t spectral_attempts=0,spectral_accepted=0,spectral_legs=0,spectral_reused=0;
  std::size_t spectral_subdivisions=0,spectral_accepted_subsegments=0;
  std::vector<std::string> spectral_rejections;
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
  int adjoint_child_loss,factored_child_loss,common_high_child_loss;
  epsilon_demand::Highs child_component_loss;
};
struct EndpointEpsilonBasis {
  PhysicalBasisTransform basis;
  std::vector<int> additional_shifts;
  int original_child_loss,estimated_child_loss;
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

// Only exact coordinates are durable. Endpoint acceptance, epsilon demands,
// precision and approximation scores are deliberately absent from the payload.
namespace basis_persistence {
namespace json=boost::json;
inline constexpr const char* revision="physical-basis-demand-selection-v3";
inline constexpr const char* verifier="physical-basis-two-inverses-and-gauge-v1";
inline constexpr const char* scope="exact_physical_basis_change";
inline std::filesystem::path directory(const NumericalOptions& options) {
  if(!options.basis_cache_directory.empty())return options.basis_cache_directory;
  return options.endpoint_cache_directory.empty()?std::filesystem::path{}:
    options.endpoint_cache_directory/"physical-bases-v1";
}
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
  id.scientific_inputs["selection_policy"]="component demand first; expression cost second; exact endpoint and contour recheck";
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
    // Certificates belong to this immutable retained frame and functional set.
    // Preserve them across higher epsilon requests; a later resource fallback
    // must not undo an already proved principal-part cancellation.
    std::optional<Rational> certified_operator_point;
    std::vector<std::optional<int>> certified_operator_lows;
  };
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
    std::vector<std::string> constraint_labels;
    bool has_physical_basis=false;
    epsilon_demand::Highs component_child_loss;
    std::optional<Matrix> rewritten_lower,rewritten_upper;
  };
 public:
  explicit Evaluator(const Graph& graph,NumericalOptions options={})
      :graph_(graph),options_(std::move(options)),plans_(graph.nodes.size()),cache_(graph.nodes.size()),expressions_(graph.nodes.size()) {
    if(options_.adjoint.continuation)
      throw std::invalid_argument("a continuation belongs to one adjoint arm; use ordinary_cache_directory for recursive continuation");
    if(graph_.nodes.empty() || options_.working_bits<64 || options_.working_bits>1000000 ||
        options_.endpoint_order<1 || options_.ordinary_order<8 || options_.ordinary_order>1000 ||
        !options_.max_epsilon || options_.max_epsilon>100 || !options_.max_refinements || !options_.max_overlap_halvings ||
        options_.overlap<=Rational(0) || options_.overlap>=Rational("1/2") ||
        options_.contour_height<=Rational(0) || options_.endpoint_constraint_tolerance<=Rational(0) ||
        !options_.max_endpoint_constraints || options_.max_endpoint_constraints>5000 ||
        !std::isfinite(options_.basis_scan_seconds) || options_.basis_scan_seconds<=0 || options_.basis_scan_seconds>60)
      throw std::invalid_argument("recursive numerical options or finite budgets");
    auto [xi,ei]=path_epsilon_variables(graph_.dimension);xi_=xi;ei_=ei;
    if(!options_.physical_bases.empty() && !options_.observable_adjoint)
      throw std::invalid_argument("general physical bases require the adjoint or factored recursive route");
    for(const auto& [depth,basis]:options_.physical_bases)
      if(depth>=graph_.nodes.size() || graph_.nodes[depth].scalar_leaf)
        throw std::invalid_argument("general physical basis requires a nonleaf depth");
    auto d=graph_.dimension.substitute(exact_point(graph_.dimension,ei_,graph_.dimension.constant(0))).rational();
    if(d.str().find('/')!=std::string::npos)throw std::invalid_argument("recursive numerical base dimension must be integral");
    d0_=std::stoi(d.str());
    if(!(graph_.dimension==graph_.dimension.constant(d0_)-graph_.dimension.constant(2)*graph_.dimension.variable(ei_)))
      throw std::invalid_argument("recursive numerical dimension must be d0-2eps");
    prescription_=options_.causal_prescription?*options_.causal_prescription:example_prescription();
    prescription_.validate(graph_.nodes.size());
    if(!options_.endpoint_cache_directory.empty())
      endpoint_store_=std::make_unique<artifacts::Store>(options_.endpoint_cache_directory,options_.endpoint_cache_max_bytes);
    const auto basis_directory=basis_persistence::directory(options_);
    if(options_.automatic_basis_scan && options_.observable_adjoint && !basis_directory.empty())
      basis_store_=std::make_unique<artifacts::Store>(basis_directory,options_.endpoint_cache_max_bytes);
  }
  const NumericalStatistics& statistics()const{return statistics_;}
  const std::vector<BasisSelection>& basis_selections()const{return basis_selections_;}
  // Snapshot a completed transform together with its immutable leaf source.
  // Copying the transform prevents a later refinement (or caller edit) from
  // changing the saved expression. No evaluation or materialization occurs.
  linear_boundary::Expression linear_expression(std::size_t depth=0)const {
    if(!options_.observable_adjoint || depth>=expressions_.size() || !expressions_[depth])
      throw std::logic_error("recursive linear expression is not available");
    return *expressions_[depth];
  }
  const causal::Prescription& causal_prescription()const{return prescription_;}
  const EndpointGeometry& endpoint_geometry(std::size_t depth)const {
    if(depth>=plans_.size() || !plans_[depth])throw std::out_of_range("recursive endpoint geometry is not prepared");
    return plans_[depth]->geometry;
  }
  // Read-only consumer fixture: no child evaluation or ordinary transport.
  // The exact connection binds the returned numerical endpoint to its basis.
  PreparedFactoredConsumers prepare_factored_consumers(std::size_t depth,unsigned high) {
    if(depth>=graph_.nodes.size() || graph_.nodes[depth].scalar_leaf || high>options_.max_epsilon)
      throw std::invalid_argument("factored consumer fixture depth or order");
    PrecisionScope precision(options_.working_bits);auto& p=plan(depth,0);
    Matrix beta(p.lower.functionals.size(),std::vector<Exact>(p.gauge.matrix.size(),graph_.dimension.constant(0)));
    for(std::size_t i=0;i<p.beta_indices.size();++i)beta[p.beta_indices[i]]=p.beta_rows[i];
    return {p.gauge.matrix,std::move(beta),cached_endpoint_operator(p,true,high)};
  }
  // Prepare one local endpoint/transport problem without evaluating a child.
  // high is the endpoint/operator epsilon top before the final inverse-D
  // gauge shift, not the requested physical integral's epsilon order.
  PreparedAdjointStage prepare_adjoint_stage(std::size_t depth,unsigned high) {
    if(!options_.observable_adjoint || depth>=graph_.nodes.size() || high>options_.max_epsilon ||
        graph_.nodes[depth].scalar_leaf || graph_.nodes[depth].closure.ordered_basis.empty())
      throw std::invalid_argument("local adjoint preparation depth, method or epsilon budget");
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
    if(!options_.observable_adjoint || depth>=graph_.nodes.size() || graph_.nodes[depth].scalar_leaf ||
        graph_.nodes[depth].closure.ordered_basis.empty())
      throw std::invalid_argument("epsilon demand profile requires a nonleaf adjoint or factored level");
    PrecisionScope precision(options_.working_bits);
    auto& prepared=plan(depth,0);preflight_linear_demands(depth,prepared);
    const int l=cached_endpoint_operator(prepared,false,0).low;
    const int r=cached_endpoint_operator(prepared,true,0).low;
    int q=0;for(const auto& row:prepared.beta_rows)for(const auto& value:row)
      if(!value.is_zero())q=std::min(q,checked(*exact_epsilon_valuation(value,ei_)));
    const int loss=checked(std::max<std::int64_t>(0,prepared.largest_shift));
    return {depth,prepared.gauge.matrix.size(),l,r,q,loss,
      *prepared.adjoint_child_loss,*prepared.factored_child_loss,
      std::max(*prepared.factored_child_loss,checked(loss+std::max(-static_cast<long>(l),-static_cast<long>(r)-q))),
      prepared.component_child_loss};
  }
  // Optimize integer epsilon powers using the actual endpoint consumers.
  // A retained nonzero ball is treated conservatively as nonzero; no small
  // numerical coefficient is discarded. Positive orders are bounded by zero.
  // The proposed gauge is exact, but the demand estimate must be rechecked
  // by rebuilding endpoints before adopting it for a physical calculation.
  EndpointEpsilonBasis optimize_endpoint_epsilon_basis(std::size_t depth) {
    const auto profile=epsilon_demand_profile(depth);
    PrecisionScope precision(options_.working_bits);auto& p=*plans_[depth];
    const auto d=p.gauge.matrix.size();const auto sample=graph_.dimension;
    std::vector<int> consumer(d,0),input(d,0);
    for(bool upper:{false,true}) {
      const auto& rows=cached_endpoint_operator(p,upper,0);
      for(const auto& row:rows.coefficients)for(unsigned j=0;j<d;++j)
        for(int k=rows.low;k<=rows.high;++k)if(!row[j][k-rows.low].is_zero())
          consumer[j]=std::min(consumer[j],k);
    }
    for(const auto& row:p.beta_rows)for(unsigned j=0;j<d;++j)if(!row[j].is_zero())
      consumer[j]=std::min(consumer[j],checked(*exact_epsilon_valuation(row[j],ei_)));
    for(unsigned i=0;i<d;++i)for(const auto& entry:p.inverse_diagonal[i])if(!entry.is_zero()) {
      const auto at=entry.substitute(exact_point(entry,xi_,sample.constant(graph_.nodes[depth].anchor)));
      if(!at.is_zero())input[i]=std::min(input[i],checked(*exact_epsilon_valuation(at,ei_)));
    }
    struct Edge {unsigned from,to;int weight;};std::vector<Edge> connection;
    for(unsigned i=0;i<d;++i)for(unsigned j=0;j<d;++j)if(!p.gauge.matrix[i][j].is_zero())
      connection.push_back({j,i,checked(*exact_epsilon_valuation(p.gauge.matrix[i][j],ei_))});
    int direct_loss=0;const auto& node=graph_.nodes[depth];
    for(unsigned i=0;i<node.operations.size();++i)if(node.operations[i].operation==feynman::Operation::Direct)
      for(const auto& entry:node.observable_rows[i])if(!entry.is_zero()) {
        const auto at=entry.substitute(exact_point(entry,xi_,sample.constant(node.anchor)));
        if(!at.is_zero())direct_loss=std::max(direct_loss,checked(-*exact_epsilon_valuation(at,ei_)));
      }
    std::vector<int> shifts(d,0);int estimate=profile.factored_child_loss;bool found=false;
    for(int total=0;total<=profile.factored_child_loss && !found;++total)
      for(int a=0;a<=total && !found;++a) {
        auto edges=connection;const auto anchor=static_cast<unsigned>(d);
        for(unsigned i=0;i<d;++i) {
          edges.push_back({anchor,i,a+input[i]});
          edges.push_back({i,anchor,total-a+consumer[i]});
        }
        std::vector<int> distance(d+1,0);bool feasible=true;
        for(unsigned pass=0;pass<=d;++pass) {
          bool changed=false;
          for(const auto& e:edges)if(distance[e.to]>distance[e.from]+e.weight) {
            distance[e.to]=distance[e.from]+e.weight;changed=true;
          }
          if(!changed)break;if(pass==d)feasible=false;
        }
        if(feasible) {
          for(unsigned i=0;i<d;++i)shifts[i]=distance[i]-distance[anchor];
          estimate=std::max(total,direct_loss);found=true;
        }
      }
    if(!found)throw std::logic_error("endpoint epsilon power constraints lost the original feasible gauge");
    auto to=p.inverse_diagonal,from=p.diagonal,connection_new=p.gauge.matrix;
    for(unsigned i=0;i<d;++i)for(unsigned j=0;j<d;++j) {
      to[i][j]=epsilon_gauge_detail::multiply_power(to[i][j],ei_,-shifts[i]);
      from[i][j]=epsilon_gauge_detail::multiply_power(from[i][j],ei_,shifts[j]);
      connection_new[i][j]=epsilon_gauge_detail::multiply_power(connection_new[i][j],ei_,shifts[j]-shifts[i]);
    }
    PhysicalBasisTransform basis{graph_.nodes[depth].closure.matrix,std::move(to),std::move(from),std::move(connection_new)};
    verify_physical_basis(basis,graph_.nodes[depth].closure.matrix,xi_);
    return {std::move(basis),std::move(shifts),profile.factored_child_loss,estimate};
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
    const auto& node=graph_.nodes[depth];++statistics_.numeric_evaluations;
    if(node.scalar_leaf) {
      report(depth,"certified scalar leaf",desired_top);
      feynman::CertifiedDeepestBetaOptions leaf;leaf.working_bits=options_.working_bits;
      leaf.f_rim=prescription_.f_rim;
      leaf.requested_digits=options_.leaf_digits;leaf.taylor_order=options_.ordinary_order;
      auto value=evaluate_leaf(node,graph_.dimension,d0_,desired_top,leaf,options_.gaussian);
      ++statistics_.leaf_evaluations;
      auto result=LaurentBoundary{value.epsilon_low,std::move(value.values),value.taylor_tail_certified};
      if(options_.observable_adjoint) {
        auto source=std::make_shared<const LaurentBoundary>(result);
        expressions_[depth]=linear_boundary::identity(source,checked(static_cast<long>(desired_top)-result.low));
      }
      cache_[depth]=std::move(result);return *cache_[depth];
    }
    if(node.closure.ordered_basis.empty()) {
      if(options_.observable_adjoint) {
        auto source=LaurentBoundary{0,Boundary(1,std::vector<B>(desired_top+1,B(0))),true};source.values[0][0]=B(1);
        expressions_[depth]=linear_boundary::Expression{
          LaurentRows{0,static_cast<int>(desired_top),std::vector(node.requested.size(),std::vector(1,std::vector<B>(desired_top+1,B(0))))},
          std::make_shared<const LaurentBoundary>(std::move(source))};
      }
      cache_[depth]=LaurentBoundary{0,Boundary(node.requested.size(),std::vector<B>(desired_top+1,B(0))),true};return *cache_[depth];
    }
    if(depth+1>=graph_.nodes.size() || graph_.nodes[depth+1].requested!=node.closure.ordered_basis)
      throw std::logic_error("recursive child does not own the parent's ordered master demand");
    if(std::all_of(node.operations.begin(),node.operations.end(),[](const auto& operation){return operation.operation==feynman::Operation::Direct;})) {
      int demand=desired_top;
      // Direct rows already expose their exact pole demand; discover it before
      // recursively transporting a child at an insufficient materialized top.
      for(const auto& row:node.observable_rows)for(const auto& coefficient:row)if(!coefficient.is_zero()) {
        const auto at=coefficient.substitute(exact_point(coefficient,xi_,coefficient.constant(node.anchor)));
        if(!at.is_zero())demand=std::max(demand,checked(static_cast<long>(desired_top)-*exact_epsilon_valuation(at,ei_)));
      }
      for(unsigned attempt=0;attempt<options_.max_refinements;++attempt) {
        if(demand>static_cast<int>(options_.max_epsilon))throw std::runtime_error("recursive direct child epsilon demand exceeds finite budget");
        auto child=evaluate_child(depth+1,static_cast<unsigned>(demand));
        try {
          if(options_.observable_adjoint) {
            auto local=exact_laurent_rows(node.observable_rows,graph_.dimension.constant(node.anchor),checked(static_cast<long>(desired_top)-child.low));
            return compose_child_expression(depth,std::move(local),desired_top);
          }
          auto result=apply_rational_rows(node.observable_rows,graph_.dimension.constant(node.anchor),child,desired_top);
          result.taylor_tail_certified=false;cache_[depth]=std::move(result);return *cache_[depth];
        }catch(const BoundaryDemand& request) {
          if(request.required_high<=demand)throw std::runtime_error("recursive direct demand made no progress");
          demand=request.required_high;++statistics_.refinements;report(depth,"child refinement",demand);
        }
      }
      throw std::runtime_error("recursive direct epsilon refinement budget exhausted");
    }
    if(!plans_[depth])report(depth,"endpoint preparation",desired_top);
    auto& prepared=plan(depth,desired_top);
    if(options_.observable_adjoint) {
      auto method=options_.linear_method;
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
    const auto x=graph_.dimension.variable(xi_);
    auto anchor=x.constant(node.anchor),h=x.constant(prepared.geometry.overlap),end=x.constant(Rational(1)-prepared.geometry.overlap);
    const auto h_ball=B::from_strings(prepared.geometry.overlap.str());
    int lower_top=constant_demand(prepared.lower,desired_top),upper_top=constant_demand(prepared.upper,desired_top);
    int transport_top=std::max({static_cast<int>(desired_top),lower_top,upper_top});
    for(const auto& row:prepared.beta_rows)for(const auto& c:row)if(!c.is_zero())
      transport_top=std::max(transport_top,checked(static_cast<long>(desired_top)-*exact_epsilon_valuation(c,ei_)));
    for(unsigned refinement=0;refinement<options_.max_refinements;++refinement) {
      const int child_top=checked(transport_top+prepared.largest_shift);
      if(child_top>static_cast<int>(options_.max_epsilon))throw std::runtime_error("recursive child epsilon demand exceeds finite budget");
      report(depth,"child boundary",child_top);
      auto child=evaluate_child(depth+1,static_cast<unsigned>(child_top));
      try {
        std::vector<std::optional<LaurentBoundary>> outputs(node.operations.size());
        for(std::size_t i=0;i<node.operations.size();++i)if(node.operations[i].operation==feynman::Operation::Direct)
          outputs[i]=apply_rational_rows(Matrix{node.observable_rows[i]},anchor,child,desired_top);
        auto gauged=apply_rational_rows(prepared.inverse_diagonal,anchor,child,transport_top);
        report(depth,"lower transport",transport_top);
        auto lower_state=transport_laurent(prepared.gauge.matrix,std::move(gauged),path(anchor,h,depth),options_.ordinary_order);
        auto lower_physical=apply_rational_rows(prepared.diagonal,h,lower_state,transport_top);
        report(depth,"lower matching",lower_top);
        auto lower_constants=match(prepared.lower,lower_physical,h_ball,lower_top);
        std::optional<RegularIntegrals> middle;
        LaurentBoundary upper_state;
        report(depth,"middle transport",transport_top);
        if(!prepared.beta_rows.empty()) {
          middle=integrate_regular_rows(prepared.gauge.matrix,prepared.beta_rows,lower_state,path(h,end,depth),desired_top,options_.ordinary_order);
          upper_state=middle->endpoint;
        }else upper_state=transport_laurent(prepared.gauge.matrix,std::move(lower_state),path(h,end,depth),options_.ordinary_order);
        auto upper_physical=apply_rational_rows(prepared.diagonal,end,upper_state,transport_top);
        report(depth,"upper matching",upper_top);
        auto upper_constants=match(prepared.upper,upper_physical,h_ball,upper_top);
        for(std::size_t i=node.operations.size();i<prepared.lower.functionals.size();++i) {
          if(prepared.lower.functionals[i])check_endpoint_constraint(
            apply(prepared.lower,*prepared.lower.functionals[i],lower_constants,h_ball,desired_top),0,
            prepared.constraint_labels[i-node.operations.size()]);
          if(prepared.upper.functionals[i])check_endpoint_constraint(
            apply(prepared.upper,*prepared.upper.functionals[i],upper_constants,h_ball,desired_top),0,
            prepared.constraint_labels[i-node.operations.size()]);
        }
        for(std::size_t i=0;i<node.operations.size();++i) {
          const auto operation=node.operations[i].operation;
          if(operation==feynman::Operation::BetaIntegral) {
            auto lower=apply(prepared.lower,*prepared.lower.functionals[i],lower_constants,h_ball,desired_top);
            auto upper=apply(prepared.upper,*prepared.upper.functionals[i],upper_constants,h_ball,desired_top);
            auto position=std::find(prepared.beta_indices.begin(),prepared.beta_indices.end(),i)-prepared.beta_indices.begin();
            LaurentBoundary mid{middle->integrals.low,Boundary{middle->integrals.values.at(position)},false};
            outputs[i]=sum({lower,mid,upper},desired_top);
          }else if(operation==feynman::Operation::LowerLimit)
            outputs[i]=apply(prepared.lower,*prepared.lower.functionals[i],lower_constants,h_ball,desired_top);
          else if(operation==feynman::Operation::UpperLimit)
            outputs[i]=apply(prepared.upper,*prepared.upper.functionals[i],upper_constants,h_ball,desired_top);
        }
        int low=0;for(const auto& output:outputs) {if(!output)throw std::logic_error("recursive operation was not executed");low=std::min(low,output->low);}
        LaurentBoundary result{low,Boundary(outputs.size(),std::vector<B>(static_cast<int>(desired_top)-low+1,B(0))),false};
        for(std::size_t i=0;i<outputs.size();++i)for(int k=outputs[i]->low;k<=static_cast<int>(desired_top);++k)
          result.values[i][k-low]=outputs[i]->values[0][k-outputs[i]->low];
        cache_[depth]=std::move(result);return *cache_[depth];
      }catch(const BoundaryDemand& demand) {
        const auto next=std::max(transport_top+1,demand.required_high);
        if(next<=transport_top)throw std::runtime_error("recursive epsilon refinement made no progress");
        transport_top=checked(next);++statistics_.refinements;report(depth,"child refinement",transport_top);
      }
    }
    throw std::runtime_error("recursive child epsilon refinement budget exhausted");
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
  NumericalStatistics statistics_;
  std::unique_ptr<artifacts::Store> endpoint_store_,basis_store_;
  std::vector<std::unique_ptr<Plan>> plans_;
  std::vector<BasisSelection> basis_selections_;
  std::vector<std::optional<LaurentBoundary>> cache_;
  std::vector<std::optional<linear_boundary::Expression>> expressions_;
  LaurentRows ordinary_transport(const Matrix& matrix,LaurentRows initial,const Matrix& forcing,
      const std::vector<Exact>& vertices,const AdjointOptions& settings) {
    if(options_.ordinary_problem_observer)options_.ordinary_problem_observer(matrix,initial,forcing,vertices);
    const auto start=std::chrono::steady_clock::now();
    if(options_.ordinary_fraction_progress)options_.ordinary_fraction_progress(0);
    const auto finish=[&]{statistics_.ordinary_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
      if(options_.ordinary_fraction_progress)options_.ordinary_fraction_progress(1);};
    if(options_.ordinary_method!=OrdinaryMethod::taylor && !settings.continuation && !settings.continuation_observer && !settings.chart_observer) {
      // Existing checkpoints remain authoritative and are verified by their original adapter.
      bool saved=false;
      if(!options_.ordinary_cache_directory.empty())saved=std::filesystem::exists(options_.ordinary_cache_directory/(adjoint_checkpoint::identity(matrix,initial,forcing,vertices,settings)+".json"));
      if(!saved){++statistics_.spectral_attempts;ft_spectral::Diagnostics diagnostics;
        auto spectral=options_.spectral;
        if(!spectral.accuracy_goal)spectral.accuracy_goal=options_.leaf_digits+12;
        if(options_.ordinary_method==OrdinaryMethod::automatic){spectral.conservative=true;spectral.max_block_size=1;spectral.max_nodes=spectral.accuracy_goal>50?128:64;spectral.endpoint_clustering=false;spectral.seconds_budget=std::min(spectral.seconds_budget,2.0);}
        std::optional<LaurentRows> result;
        if(options_.ordinary_method==OrdinaryMethod::ultraspherical) {
          ft_ultraspherical::Options coefficient;
          static_cast<ft_spectral::Options&>(coefficient)=spectral;
          coefficient.diagonal_gauge=false;coefficient.endpoint_clustering=false;
          if(options_.ordinary_fraction_progress)coefficient.progress=[&](unsigned leg,double fraction) {
            options_.ordinary_fraction_progress((leg+fraction)/static_cast<double>(std::max<std::size_t>(1,vertices.size()-1)));
          };
          result=ft_ultraspherical_checkpoint::try_transport(matrix,initial,forcing,vertices,coefficient,diagnostics,settings,
            {options_.ordinary_cache_directory,options_.ordinary_cache_max_bytes,&statistics_.spectral_reused});
        } else result=ft_spectral_checkpoint::try_transport(matrix,initial,forcing,vertices,spectral,diagnostics,settings,
            {options_.ordinary_cache_directory,options_.ordinary_cache_max_bytes,&statistics_.spectral_reused});
        statistics_.spectral_seconds+=diagnostics.preparation_seconds+diagnostics.numerical_seconds;
        if(result){++statistics_.spectral_accepted;statistics_.spectral_legs+=diagnostics.legs;
          statistics_.spectral_subdivisions+=diagnostics.subdivisions;
          statistics_.spectral_accepted_subsegments+=diagnostics.accepted_subsegments;
          if(options_.ordinary_progress)options_.ordinary_progress("spectral transport accepted: "+std::to_string(diagnostics.legs)+" contour legs, "+std::to_string(diagnostics.subdivisions)+" subdivisions");
          finish();return std::move(*result);}
        statistics_.spectral_rejections.push_back(diagnostics.reason);
        if(options_.ordinary_progress){
          const auto largest=diagnostics.block_sizes.empty()?0:*std::max_element(diagnostics.block_sizes.begin(),diagnostics.block_sizes.end());
          const auto nodes=diagnostics.nodes.empty()?0:*std::max_element(diagnostics.nodes.begin(),diagnostics.nodes.end());
          options_.ordinary_progress("spectral transport fallback: "+diagnostics.reason+"; largest block "+std::to_string(largest)+", maximum nodes "+std::to_string(nodes));
        }
      }
    }
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
    finish();return result;
  }
  causal::Prescription example_prescription()const {
    // A familiar name cannot authorize a prescription for altered kinematics.
    // Compare the actual incoming coordinates with the native fixture before
    // using its positivity proof or its supplied Henn contour data.
    auto example=graph_.definition?*graph_.definition:feynman::example_family(graph_.family_name);
    const ibp::PropagatorBasis expected(ibp::quadratic_family(example.momenta,graph_.dimension,example.physical_count));
    const auto& actual=graph_.nodes.front().incoming;
    bool same=actual.physical_count==expected.physical_count && actual.space.loops==expected.space.loops &&
      actual.space.external_gram==expected.space.external_gram && actual.denominators.size()==expected.denominators.size() &&
      d0_==example.dimension_at_epsilon_zero;
    if(same)for(std::size_t i=0;i<actual.denominators.size();++i)
      same=same&&pullback::equal(actual.denominators[i],expected.denominators[i]);
    if(!same)throw std::invalid_argument("altered example kinematics require an explicit causal prescription");
    std::vector<Rational> anchors;bool standard=true;
    for(const auto& node:graph_.nodes) {anchors.push_back(node.anchor);standard=standard&&node.left==0&&node.right==1;}
    if(!graph_.definition && graph_.family_name=="henn_double_pentagon_x0")
      return causal::current_example(graph_.family_name,graph_.nodes.size(),anchors,standard);
    return causal::euclidean_family(example,graph_.nodes.size());
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
    if(!options_.observable_adjoint || !options_.direct_endpoint_functionals || options_.endpoint_order>256)
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
  Endpoint endpoint(const Matrix& matrix,const Matrix& diagonal,const Node& node,bool upper,std::size_t depth,int high) {
    auto fuchs=fuchsify::prepare(matrix,xi_,options_.fuchsification);
    if(!fuchs.success)throw std::domain_error("recursive endpoint fuchsification unsupported: "+fuchs.reason);
    const std::string side=upper?"upper":"lower";
    if(auto result=direct_endpoint(matrix,diagonal,fuchs,node,upper)) {
      ++statistics_.direct_endpoints;
      report(depth,side+" direct endpoint functionals; no fundamental series",high);
      return std::move(*result);
    }
    report(depth,side+" exact endpoint series",high);
    auto endpoint_options=options_.endpoint;
    const auto caller_progress=endpoint_options.column_progress;
    endpoint_options.column_progress=[&,caller_progress](unsigned completed,unsigned total){
      if(caller_progress)caller_progress(completed,total);
      report(depth,side+" exact endpoint columns "+std::to_string(completed)+"/"+std::to_string(total),high);
    };
    auto series=[&] {
      if(endpoint_store_) {
        auto verification=options_.endpoint_cache_verification;
        const auto caller_verification=verification.column_progress;
        verification.column_progress=[&,caller_verification](unsigned completed,unsigned total,std::size_t products){
          if(caller_verification)caller_verification(completed,total,products);
          report(depth,side+" verified endpoint columns "+std::to_string(completed)+"/"+std::to_string(total)+
            "; polynomial product estimate "+std::to_string(products),high);
        };
        auto cached=cached_affine::prepare(fuchs.matrix,xi_,ei_,options_.endpoint_order,endpoint_options,
            *endpoint_store_,verification);
        if(cached.cache_hit)++statistics_.endpoint_series_reused;else ++statistics_.endpoint_series_built;
        report(depth,side+(cached.cache_hit?" endpoint series verified from cache":" endpoint series verified and saved"),high);
        return std::move(cached.series);
      }
      ++statistics_.endpoint_series_built;
      return AffineFrobeniusSeries::prepare(fuchs.matrix,xi_,ei_,options_.endpoint_order,endpoint_options);
    }();
    auto physical_gauge=fuchsify::detail::multiply(diagonal,fuchs.transform);
    // In the adjoint route P*(D*T*F)^-1*D = P*(T*F)^-1.
    // Keep D in the physical observable P, while cancelling it from the
    // receiving frame before epsilon normalization and numerical inversion.
    const auto& receiving_gauge=options_.observable_adjoint?fuchs.transform:physical_gauge;
    auto matching_frame=series.project(receiving_gauge);
    if(!matching_frame.wronskian_prefactor){
      matching_frame.wronskian_prefactor=series.gauged_wronskian_prefactor(receiving_gauge,fuchs.matrix,matrix);
      if(matching_frame.wronskian_prefactor)report(depth,side+" Wronskian certified from exact trace and rational factors",high);
    }
    Endpoint result{std::move(series),std::move(matching_frame),std::vector<std::optional<Expansion>>(node.operations.size()),{}};
    const auto collect=[&](const Matrix& coefficients){for(const auto& row:coefficients)for(const auto& c:row)if(!c.is_zero())result.clearance_coefficients.push_back(c);};
    collect(matrix);collect(fuchs.matrix);collect(fuchs.transform);collect(fuchs.inverse_transform);collect(physical_gauge);
    const auto x=graph_.dimension.variable(xi_);
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
      }
      if((i+1)%32==0 || i+1==node.operations.size())
        report(depth,side+" endpoint observables "+std::to_string(i+1)+"/"+std::to_string(node.operations.size()),high);
    }
    return result;
  }
  // Screening only: pull the ORIGINAL plan's ordinary-cut functionals into
  // proposed coordinates. This does not construct or certify new endpoint
  // frames, and must never be installed as an execution plan.
  std::unique_ptr<Plan> screen_basis_plan(Plan& base,const PhysicalBasisTransform& basis) {
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
    // Y = D_base U_base = D_new U_new, hence U_base=V U_new.
    const auto v=fuchsify::detail::multiply(base.inverse_diagonal,candidate->diagonal);
    candidate->geometry=base.geometry;
    candidate->lower.functionals.resize(base.lower.functionals.size());
    candidate->upper.functionals.resize(base.upper.functionals.size());
    candidate->constraint_labels=base.constraint_labels;
    candidate->beta_indices=base.beta_indices;
    if(!base.beta_rows.empty())candidate->beta_rows=fuchsify::detail::multiply(base.beta_rows,v);
    candidate->has_physical_basis=true;
    for(const auto& row:candidate->inverse_diagonal)for(const auto& entry:row)if(!entry.is_zero())
      candidate->largest_shift=std::max(candidate->largest_shift,-*exact_epsilon_valuation(entry,ei_));
    auto endpoint_screen=[&](bool upper) {
      const auto h=graph_.dimension.constant(base.geometry.overlap);
      const auto at=upper?graph_.dimension.constant(1)-h:h;
      // Discover V's exact pole depth at the actual cut before requesting E.
      auto inner=exact_laurent_rows(v,at,0);
      const auto& outer=cached_endpoint_operator(base,upper,std::max(0,checked(-static_cast<long>(inner.low))));
      const auto inner_high=std::max(0,checked(-static_cast<long>(outer.low)));
      if(inner.high<inner_high)inner=exact_laurent_rows(v,at,inner_high);
      return linear_boundary::compose(outer,inner,0);
    };
    candidate->lower_operator=endpoint_screen(false);
    candidate->upper_operator=endpoint_screen(true);
    return candidate;
  }
  Plan& plan(std::size_t depth,unsigned desired_top) {
    if(plans_[depth])return *plans_[depth];const auto& node=graph_.nodes[depth];
    const auto explicit_basis=options_.physical_bases.find(depth);
    auto selected=build_plan(depth,desired_top,explicit_basis==options_.physical_bases.end()?nullptr:&explicit_basis->second);
    if(options_.automatic_basis_scan && options_.observable_adjoint && explicit_basis==options_.physical_bases.end()) {
      BasisSelection selection;selection.depth=depth;selection.masters=selected->gauge.matrix.size();
      const auto characters=[](const Matrix& matrix) {std::size_t size=0;for(const auto& row:matrix)for(const auto& value:row)if(!value.is_zero())size+=value.str().size();return size;};
      selection.original_characters=characters(selected->gauge.matrix);
      std::optional<artifacts::Identity> basis_identity;
      std::optional<basis_persistence::Choice> cached_basis;
      std::optional<PhysicalBasisTransform> accepted_basis;
      if(basis_store_) {
        basis_identity=basis_persistence::identity(graph_,depth,xi_,ei_,prescription_);
        // Corruption or a failed exact identity is fatal, not a fresh scan that
        // silently changes coordinates underneath an existing continuation.
        cached_basis=basis_persistence::load(*basis_store_,*basis_identity,node.closure.matrix,xi_);
        if(cached_basis){selection.cache_hit=true;++statistics_.basis_cache_hits;}
      }
      report(depth,cached_basis?"rechecking persisted exact master basis":"scanning simpler master bases",desired_top);
      const auto began=std::chrono::steady_clock::now();
      try {
        preflight_linear_demands(depth,*selected);
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
        const auto accept=[&](const PhysicalBasisTransform& basis,bool persisted=false) -> bool {
          ++selection.full_candidates_checked;
          if(persisted)++selection.candidates_checked;
          try {
            auto candidate=build_plan(depth,desired_top,&basis);
            preflight_linear_demands(depth,*candidate);
            const auto candidate_cost=cost(*candidate);
            report(depth,"basis full validation "+std::to_string(selection.full_candidates_checked)+
              ": "+describe_cost(candidate_cost)+", overlap="+candidate->geometry.overlap.str()+
              "; current "+describe_cost(best_cost)+", overlap="+selected->geometry.overlap.str(),desired_top);
            if(!persisted && candidate->geometry.overlap<selected->geometry.overlap) {
              report(depth,"basis candidate rejected: smaller certified overlap",desired_top);return false;
            }
            if(!persisted && !(candidate_cost<best_cost)) {
              report(depth,"basis candidate rejected: full cost did not improve",desired_top);return false;
            }
            accepted_basis=basis;best_cost=candidate_cost;
            selected=std::move(candidate);selection.accepted=true;
            selection.candidate_characters=characters(selected->gauge.matrix);
            selection.candidate_child_loss=*selected->factored_child_loss;
            selection.reason="lower component epsilon demand or simpler equal-demand matrix; exact endpoint and contour checks passed";
            report(depth,std::string(persisted?"persisted basis validated: ":"basis candidate certified and accepted: ")+
              describe_cost(best_cost)+", screened="+std::to_string(selection.candidates_screened)+
              ", full="+std::to_string(selection.full_candidates_checked),desired_top);
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
          if(cached_basis->selected)(void)accept(cached_basis->basis,true);
          if(!selection.accepted)selection.reason="persisted original basis retained";
        } else {
          ++statistics_.basis_scans;
          physical_basis_scan::Options scan_options;
          scan_options.seconds=options_.basis_scan_seconds;scan_options.epsilon_variable=ei_;
          scan_options.max_candidates=8;
          auto frontier=physical_basis_scan::scan_candidates(node.closure,node.merged.physical_count,xi_,scan_options,
            [&](const physical_basis_scan::Candidate& candidate) -> std::optional<physical_basis_scan::DemandCost> {
              if(candidate.metadata.provenance=="original")return original_cost;
              ++selection.candidates_screened;++selection.candidates_checked;
              try {
                auto screen=screen_basis_plan(*selected,candidate.transform);
                if(!screen)return original_cost;
                preflight_linear_demands(depth,*screen);
                return cost(*screen);
              } catch(const std::exception& error) {
                selection.reason=std::string("candidate screen unresolved: ")+error.what();
                return std::nullopt;
              }
            });
          selection.scan_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
          std::stable_sort(frontier.candidates.begin(),frontier.candidates.end(),
            [](const auto& a,const auto& b){return *a.demand_cost<*b.demand_cost;});
          report(depth,"basis discovery complete: screened="+std::to_string(selection.candidates_screened)+
            ", frontier="+std::to_string(frontier.candidates.size())+", original "+describe_cost(original_cost),desired_top);
          const auto validation_began=std::chrono::steady_clock::now();
          const double validation_seconds=std::min(30.0,std::max(1.0,options_.basis_scan_seconds));
          std::size_t rebuilt=0;
          for(const auto& candidate:frontier.candidates) {
            if(rebuilt==8)break;
            if(!candidate.demand_cost || !(*candidate.demand_cost<best_cost))continue;
            // Cooperative bound: do not start another expensive exact rebuild
            // after the budget. An in-flight exact operation is not interrupted.
            if(std::chrono::duration<double>(std::chrono::steady_clock::now()-validation_began).count()>=validation_seconds) {
              report(depth,"basis full-validation budget reached; retaining last certified coordinates",desired_top);break;
            }
            ++rebuilt;
            report(depth,"basis full candidate "+std::to_string(rebuilt)+": predicted "+
              describe_cost(*candidate.demand_cost)+"; "+candidate.metadata.provenance,desired_top);
            // The frontier is ordered by predicted cost. Once a fully checked
            // improvement is found, stop this bounded search and persist it.
            if(accept(candidate.transform))break;
          }
        }
        if(selection.accepted) {
          ++statistics_.basis_selections;
          selection.reason="lower component epsilon lookahead or simpler equal-demand matrix; exact endpoint and contour checks passed";
        }
        else if(selection.reason.empty())selection.reason="no certified demand improvement within the bounded search";
      }catch(const std::exception& error) {
        if(cached_basis && cached_basis->selected)throw;
        if(!cached_basis && !selection.scan_seconds)selection.scan_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
        selection.reason=std::string("keeping original basis: ")+error.what();
      }
      if(basis_store_ && !cached_basis) {
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
    auto lower=endpoint(gauge.matrix,diagonal,node,false,depth,desired_top),upper=endpoint(reverse(gauge.matrix,true),reverse(diagonal),node,true,depth,desired_top);
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
    // Exact endpoint functionals already encode the prescribed singular limit.
    // At ordinary cuts their rational maps can participate in a jointly
    // certified rewrite. Numerical Frobenius operators never enter this path.
    if(options_.exact_functional_reduction && result->lower.direct_operator_exact && result->upper.direct_operator_exact) {
      namespace fr=ft_functional_reduction;
      Matrix weight(result->lower.functionals.size(),std::vector<Exact>(d,graph_.dimension.constant(0)));
      for(std::size_t i=0;i<result->beta_indices.size();++i)weight[result->beta_indices[i]]=result->beta_rows[i];
      const auto h=graph_.dimension.constant(result->geometry.overlap);
      fr::Functional functional{result->gauge.matrix,
        fr::detail::at(*result->lower.direct_operator,xi_,h),weight,
        fr::detail::at(*result->upper.direct_operator,xi_,h),
        h,graph_.dimension.constant(1)-h,xi_,ei_};
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
    report(depth,"endpoint overlap="+result->geometry.overlap.str(),desired_top);
    return result;
  }
  EndpointGeometry endpoint_clearance(const Endpoint& lower,const Endpoint& upper)const {
    using M=NativeTailMagnitude;
    EndpointGeometry result{Rational("1/2"),{},M::zero()};
    const auto& sample=graph_.dimension;auto epsilon=sample.variable(ei_),x=sample.variable(xi_);
    const auto& names=sample.variables();auto ii=std::find(names.begin(),names.end(),"I");
    std::set<std::string> seen;
    for(const auto* endpoint:{&lower,&upper})for(const auto& coefficient:endpoint->clearance_coefficients) {
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
      if(!seen.insert(p.str()).second)continue;
      auto roots=polynomial_roots(p,xi_,128);
      for(auto& root:roots) {
        auto bound=M::lower_abs(root);
        if(bound.is_zero() || !bound.is_finite())throw std::domain_error("nonzero endpoint pole cannot be separated from zero at the root-isolation precision");
        if(result.nonzero_poles.empty() || result.nearest_pole_lower>bound)result.nearest_pole_lower=bound;
        result.nonzero_poles.push_back(std::move(root));
      }
    }
    // This proves a geometric convergence margin for the formal epsilon
    // coefficients. It does not certify the omitted affine Taylor tail.
    for(unsigned halvings=0;halvings<options_.max_overlap_halvings;++halvings) {
      const auto upper=M::upper_abs(B::from_strings(result.overlap.str()));
      if(result.overlap<=options_.overlap && (result.nonzero_poles.empty() || upper*M::from_ui(16)<=result.nearest_pole_lower))return result;
      result.overlap=result.overlap/Rational(2);
    }
    throw std::runtime_error("recursive geometric endpoint-overlap budget exhausted");
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
    if(endpoint.direct_operator)
      return exact_laurent_rows(*endpoint.direct_operator,graph_.dimension.constant(exact_point_value),high);
    const auto point=B::from_strings(exact_point_value.str());
    const auto rows=endpoint.functionals.size();const auto d=endpoint.matching_frame.columns;
    if(std::none_of(endpoint.functionals.begin(),endpoint.functionals.end(),[](const auto& value){return value.has_value();}))
      return {0,high,std::vector(rows,std::vector(d,std::vector<B>(high+1,B(0))))};
    const int inverse_high=constant_demand(endpoint,high);
    auto inverse=affine_operator::prepare(*endpoint.series,endpoint.matching_frame,point,inverse_high,options_.matching);
    if(!inverse.success())throw std::domain_error("recursive endpoint operator inverse unsupported: "+inverse.reason);
    if(!endpoint.certified_operator_point || *endpoint.certified_operator_point!=exact_point_value) {
      endpoint.certified_operator_point=exact_point_value;
      endpoint.certified_operator_lows.assign(rows,std::nullopt);
    }
    const auto certificate_deadline=std::chrono::steady_clock::now()+std::chrono::milliseconds(options_.principal_part_endpoint_milliseconds);
    std::vector<std::optional<affine_operator::Operator>> operators(rows);int low=0;
    for(unsigned i=0;i<rows;++i)if(endpoint.functionals[i]) {
      auto matching=options_.matching;
      const auto remaining=std::chrono::duration_cast<std::chrono::milliseconds>(certificate_deadline-std::chrono::steady_clock::now()).count();
      if(remaining<=0 || endpoint.certified_operator_lows[i])matching.certify_operator_principal_part=false;
      else matching.principal_part_max_milliseconds=std::min<unsigned>(matching.principal_part_max_milliseconds,static_cast<unsigned>(remaining));
      operators[i]=affine_operator::compose(inverse,*endpoint.series,*endpoint.functionals[i],point,high,matching);
      if(!operators[i]->success())throw std::domain_error("recursive endpoint operator composition unsupported: "+operators[i]->reason);
      auto& op=*operators[i];statistics_.principal_part_attempts+=op.principal_part_attempted;
      statistics_.principal_part_certificates+=op.principal_part_certified;
      statistics_.principal_part_zero_coefficients+=op.principal_part_zero_coefficients;
      if(op.principal_part_certified)endpoint.certified_operator_lows[i]=op.row_lower_bounds.at(0);
      else if(endpoint.certified_operator_lows[i]) {
        const int certified=std::min(high,*endpoint.certified_operator_lows[i]);
        if(certified>op.matrix.low) {
          for(auto& row:op.matrix.coefficients)for(auto& column:row)
            column.erase(column.begin(),column.begin()+(certified-op.matrix.low));
          op.matrix.low=certified;
        }
        op.row_lower_bounds[0]=std::max(op.row_lower_bounds[0],*endpoint.certified_operator_lows[i]);
        ++statistics_.principal_part_reused;
      }
      low=std::min(low,operators[i]->matrix.low);
    }
    LaurentRows result{low,high,std::vector(rows,std::vector(d,std::vector<B>(high-low+1,B(0))))};
    for(unsigned i=0;i<rows;++i)if(operators[i]) {
      const auto& matrix=operators[i]->matrix;
      for(unsigned j=0;j<d;++j)for(int k=matrix.low;k<=high;++k)result.coefficients[i][j][k-low]=matrix.coefficients[0][j][k-matrix.low];
    }
    return result;
  }
  static LaurentRows shift_operator_columns(const LaurentRows& rows,const std::vector<std::int64_t>& shifts,bool inverse=false) {
    return shift_laurent_columns(rows,shifts,inverse);
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
      if(rewritten)saved=exact_laurent_rows(*rewritten,graph_.dimension.constant(prepared.geometry.overlap),high);
      else saved=endpoint_operator(upper?prepared.upper:prepared.lower,prepared.geometry.overlap,high);
    }
    return *saved;
  }
  void preflight_linear_demands(std::size_t depth,Plan& prepared) {
    if(prepared.adjoint_child_loss)return;
    report(depth,"structural child-window preflight",0);
    const auto l=cached_endpoint_operator(prepared,false,0).low;
    const auto r=cached_endpoint_operator(prepared,true,0).low;
    int q=0,direct=0;
    for(const auto& row:prepared.beta_rows)for(const auto& entry:row)if(!entry.is_zero())
      q=std::min(q,checked(*exact_epsilon_valuation(entry,ei_)));
    const auto& node=graph_.nodes[depth];
    for(unsigned i=0;i<node.operations.size();++i)if(node.operations[i].operation==feynman::Operation::Direct)
      for(const auto& entry:node.observable_rows[i])if(!entry.is_zero()) {
        const auto at=entry.substitute(exact_point(entry,xi_,entry.constant(node.anchor)));
        if(!at.is_zero())direct=std::min(direct,checked(*exact_epsilon_valuation(at,ei_)));
      }
    const long loss=std::max<std::int64_t>(0,prepared.largest_shift);
    // Epsilon-regular homogeneous transport preserves lower bounds, while
    // forcing can lower them to q. No numerical cancellation is used here.
    prepared.adjoint_child_loss=checked(-std::min(static_cast<long>(std::min({l,r,q}))-loss,static_cast<long>(direct)));
    // Endpoint composition and middle integration independently request the
    // child map: max(H-l,H-r,H-q), followed by the inverse D gauge.
    prepared.factored_child_loss=checked(std::max(loss+std::max({-static_cast<long>(l),-static_cast<long>(r),-static_cast<long>(q)}),-static_cast<long>(direct)));
    if(options_.component_epsilon_demands) {
      namespace ed=epsilon_demand;
      const auto d=prepared.gauge.matrix.size();
      const ed::Highs outputs(prepared.lower.functionals.size(),0);
      auto left=ed::pullback(cached_endpoint_operator(prepared,false,0),outputs);
      auto right=ed::pullback(cached_endpoint_operator(prepared,true,0),outputs);
      Matrix beta(outputs.size(),std::vector<Exact>(d,graph_.dimension.constant(0)));
      for(std::size_t i=0;i<prepared.beta_indices.size();++i)beta[prepared.beta_indices[i]]=prepared.beta_rows[i];
      auto middle=ed::close(prepared.gauge.matrix,ed::merge(right,ed::pullback(beta,outputs,ei_)),ei_);
      auto initial=ed::close(prepared.gauge.matrix,ed::merge(left,middle),ei_);
      auto inverse=prepared.inverse_diagonal;
      for(auto& row:inverse)for(auto& value:row)value=value.substitute(exact_point(value,xi_,value.constant(node.anchor)));
      prepared.component_child_loss=ed::pullback(inverse,initial,ei_);
      Matrix direct_map(outputs.size(),std::vector<Exact>(d,graph_.dimension.constant(0)));
      for(std::size_t i=0;i<node.operations.size();++i)if(node.operations[i].operation==feynman::Operation::Direct)
        for(std::size_t j=0;j<d;++j)direct_map[i][j]=node.observable_rows[i][j].substitute(exact_point(graph_.dimension,xi_,graph_.dimension.constant(node.anchor)));
      prepared.component_child_loss=ed::merge(prepared.component_child_loss,ed::pullback(direct_map,outputs,ei_));
      prepared.factored_child_loss=std::min(*prepared.factored_child_loss,std::max(0,ed::maximum(prepared.component_child_loss)));
    }
    ++statistics_.demand_preflights;
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
    const auto anchor=sample.constant(node.anchor),h=sample.constant(prepared.geometry.overlap),
      end=sample.constant(Rational(1)-prepared.geometry.overlap);
    return {prepared.gauge.matrix,std::move(forcing),std::move(upper_forcing),
      std::move(left),std::move(right),path(h,anchor,depth),path(end,anchor,depth),prepared.gauge.shifts};
  }
  LaurentBoundary evaluate_adjoint(std::size_t depth,unsigned desired_top,Plan& prepared) {
    const auto& node=graph_.nodes[depth];auto sample=graph_.dimension;
    const auto anchor=sample.constant(node.anchor),h=sample.constant(prepared.geometry.overlap),end=sample.constant(Rational(1)-prepared.geometry.overlap);
    const auto h_ball=B::from_strings(prepared.geometry.overlap.str());
    preflight_linear_demands(depth,prepared);
    int demand=checked(static_cast<long>(desired_top)+*prepared.adjoint_child_loss);
    auto& combined=prepared.local_operator;
    for(unsigned refinement=0;refinement<options_.max_refinements;++refinement) {
      if(demand>static_cast<int>(options_.max_epsilon))throw std::runtime_error("recursive adjoint child epsilon demand exceeds finite budget");
      report(depth,"child boundary",demand);auto child=evaluate_child(depth+1,demand);
      const int needed_operator=checked(static_cast<long>(desired_top)-child.low);
      if(!combined || combined->high<needed_operator) {
        const int high=checked(needed_operator+prepared.largest_shift);
        auto stage=adjoint_stage(depth,prepared,high);
        auto settings=options_.adjoint;settings.taylor_order=options_.ordinary_order;
        report(depth,"lower adjoint transport",high);auto left=ordinary_transport(stage.connection,std::move(stage.lower_endpoint),stage.lower_forcing,stage.lower_path,settings);
        observe_operator(depth,"lower transported operator",left);
        report(depth,"upper adjoint transport",high);auto right=ordinary_transport(stage.connection,std::move(stage.upper_endpoint),stage.upper_forcing,stage.upper_path,settings);
        observe_operator(depth,"upper transported operator",right);
        auto transported=add_laurent_rows(left,right);
        if(prepared.has_physical_basis) {
          const auto inverse=exact_laurent_rows(prepared.inverse_diagonal,anchor,
            std::max(0,checked(static_cast<long>(needed_operator)-transported.low)));
          combined=linear_boundary::compose(transported,inverse,needed_operator);
        } else combined=shift_operator_columns(transported,prepared.gauge.shifts,true);
        Matrix direct(prepared.lower.functionals.size(),std::vector<Exact>(prepared.gauge.matrix.size(),sample.constant(0)));
        for(unsigned i=0;i<node.operations.size();++i)if(node.operations[i].operation==feynman::Operation::Direct)direct[i]=node.observable_rows[i];
        combined=add_laurent_rows(*combined,exact_laurent_rows(direct,anchor,combined->high));
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
    const auto anchor=sample.constant(node.anchor),h=sample.constant(prepared.geometry.overlap),end=sample.constant(Rational(1)-prepared.geometry.overlap);
    const auto h_ball=B::from_strings(prepared.geometry.overlap.str());
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
          auto lower=factored_transport::evolve_components(prepared.gauge.matrix,gauged,{},path(anchor,h,depth),initial_needs,{},settings);
          observe_operator(depth,"lower factored map",lower.physical.value.transform);
          report(depth,"middle factored transport",ed::maximum(middle_needs));
          auto middle=factored_transport::evolve_components(prepared.gauge.matrix,lower.physical,beta,path(h,end,depth),right_needs,outputs,settings);
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
            path(anchor,h,depth),initial_high,settings).physical;
        observe_operator(depth,"lower factored map",lower.transform);
        report(depth,"middle factored transport",middle_high);
        auto middle=factored_transport::evolve(prepared.gauge.matrix,lower,beta,path(h,end,depth),middle_high,high,settings);
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
  int constant_demand(const Endpoint& endpoint,unsigned high)const {
    long result=high;for(const auto& functional:endpoint.functionals)if(functional)
      for(auto demand:endpoint.series->valuation_metadata(*functional).required_source_top(high))result=std::max(result,demand);
    return checked(result);
  }
  affine_matching::Boundary match(const Endpoint& endpoint,const LaurentBoundary& boundary,const B& point,int high)const {
    auto result=affine_matching::match(*endpoint.series,endpoint.matching_frame,point,
      {boundary.low,boundary.high(),boundary.values},{0,high},options_.matching);
    if(result.status==affine_matching::Status::NeedMoreBoundary)throw BoundaryDemand(result.required_boundary_high,result.reason);
    if(!result.success())throw std::domain_error("recursive endpoint matching unsupported: "+result.reason);
    return std::move(result.value);
  }
  LaurentBoundary apply(const Endpoint& endpoint,const Expansion& functional,const affine_matching::Boundary& constants,
      const B& point,unsigned high)const {
    auto result=affine_matching::apply(*endpoint.series,functional,point,constants,{0,static_cast<int>(high)},options_.matching);
    if(!result.success())throw std::domain_error("recursive endpoint functional unsupported: "+result.reason);
    return {result.value.low,std::move(result.value.coefficients),false};
  }
  std::vector<Exact> path(const Exact& from,const Exact& to,std::size_t depth)const {
    if(from==to)return {from};
    const auto& names=graph_.dimension.variables();auto ii=std::find(names.begin(),names.end(),"I");
    if(ii==names.end())throw std::invalid_argument("recursive contour requires exact I variable");
    auto lift=graph_.dimension.variable(ii-names.begin())*graph_.dimension.constant(options_.contour_height)*
      graph_.dimension.constant(prescription_.levels.at(depth).x_detour_sign);
    return {from,from+lift,to+lift,to};
  }
  static LaurentBoundary sum(const std::vector<LaurentBoundary>& values,unsigned high) {
    int low=0;for(const auto& value:values)low=std::min(low,value.low);
    LaurentBoundary out{low,Boundary(1,std::vector<B>(static_cast<int>(high)-low+1,B(0))),false};
    for(const auto& value:values)for(int k=value.low;k<=static_cast<int>(high);++k)out.values[0][k-low]+=value.values[0][k-value.low];return out;
  }
};
} // namespace diffexp::recursion
