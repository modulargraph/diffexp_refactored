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

namespace diffexp::recursion {
enum class LinearMethod { adjoint, factored, automatic };
struct NumericalOptions {
  slong working_bits=384;
  unsigned endpoint_order=32,ordinary_order=80,max_epsilon=100,max_refinements=8,max_overlap_halvings=128;
  Rational overlap{Rational("1/16")},contour_height{Rational("1/10")};
  unsigned leaf_digits=28;
  // Consistency of retained physical endpoint constraints, not a proof about
  // omitted epsilon/Taylor coefficients or exact vanishing of a period.
  Rational endpoint_constraint_tolerance{Rational("1/100000000000000000000")};
  unsigned max_endpoint_constraints=256;
  bool observable_adjoint=true;
  LinearMethod linear_method=LinearMethod::adjoint;
  AdjointOptions adjoint;
  std::optional<causal::Prescription> causal_prescription;
  std::function<void(std::size_t,const std::string&,int)> progress;
  // Optional read-only diagnostics on completed coefficient maps. This lets
  // callers locate arithmetic loss without materializing the shared source.
  std::function<void(std::size_t,const std::string&,const LaurentRows&)> operator_observer;
  std::filesystem::path endpoint_cache_directory;
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
  adjoint_checkpoint::Statistics ordinary_checkpoints;
  std::size_t exact_plans=0,leaf_evaluations=0,numeric_evaluations=0,refinements=0,cache_hits=0;
  std::size_t endpoint_series_built=0,endpoint_series_reused=0;
  std::size_t adjoint_selections=0,factored_selections=0;
  std::size_t demand_preflights=0,local_operator_reuses=0;
  std::size_t conditioning_method_fallbacks=0;
  std::size_t endpoint_constraint_rows=0,endpoint_constraint_coefficients=0;
  double maximum_endpoint_constraint_residual=0;
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
    AffineFrobeniusSeries series;
    Expansion matching_frame;
    std::vector<std::optional<Expansion>> functionals;
    std::vector<Exact> clearance_coefficients;
    std::vector<std::pair<Expansion,std::string>> pending_constraints;
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
        !options_.max_endpoint_constraints || options_.max_endpoint_constraints>5000)
      throw std::invalid_argument("recursive numerical options or finite budgets");
    auto [xi,ei]=path_epsilon_variables(graph_.dimension);xi_=xi;ei_=ei;
    auto d=graph_.dimension.substitute(exact_point(graph_.dimension,ei_,graph_.dimension.constant(0))).rational();
    if(d.str().find('/')!=std::string::npos)throw std::invalid_argument("recursive numerical base dimension must be integral");
    d0_=std::stoi(d.str());
    if(!(graph_.dimension==graph_.dimension.constant(d0_)-graph_.dimension.constant(2)*graph_.dimension.variable(ei_)))
      throw std::invalid_argument("recursive numerical dimension must be d0-2eps");
    prescription_=options_.causal_prescription?*options_.causal_prescription:example_prescription();
    prescription_.validate(graph_.nodes.size());
    if(!options_.endpoint_cache_directory.empty())
      endpoint_store_=std::make_unique<artifacts::Store>(options_.endpoint_cache_directory,options_.endpoint_cache_max_bytes);
  }
  const NumericalStatistics& statistics()const{return statistics_;}
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
  // Prepare one local endpoint/transport problem without evaluating a child.
  // high is the endpoint/operator epsilon top before the final inverse-D
  // gauge shift, not the requested physical integral's epsilon order.
  PreparedAdjointStage prepare_adjoint_stage(std::size_t depth,unsigned high) {
    if(!options_.observable_adjoint || depth>=graph_.nodes.size() || high>options_.max_epsilon ||
        graph_.nodes[depth].scalar_leaf || graph_.nodes[depth].closure.ordered_basis.empty())
      throw std::invalid_argument("local adjoint preparation depth, method or epsilon budget");
    PrecisionScope precision(options_.working_bits);
    if(!plans_[depth])report(depth,"endpoint preparation",high);
    auto& prepared=plan(depth,high);
    if(!prepared.constraint_labels.empty())throw std::domain_error(
      "standalone adjoint-stage export cannot discharge physical endpoint constraints; use the full recursive evaluator");
    return adjoint_stage(depth,prepared,high);
  }
  LaurentBoundary evaluate(unsigned desired_top=0){return evaluate(0,desired_top);}
  LaurentBoundary evaluate(std::size_t depth,unsigned desired_top) {
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
  std::unique_ptr<artifacts::Store> endpoint_store_;
  std::vector<std::unique_ptr<Plan>> plans_;
  std::vector<std::optional<LaurentBoundary>> cache_;
  std::vector<std::optional<linear_boundary::Expression>> expressions_;
  LaurentRows ordinary_transport(const Matrix& matrix,LaurentRows initial,const Matrix& forcing,
      const std::vector<Exact>& vertices,const AdjointOptions& settings) {
    return adjoint_checkpoint::transport(matrix,std::move(initial),forcing,vertices,settings,
      {options_.ordinary_cache_directory,options_.ordinary_cache_max_bytes,&statistics_.ordinary_checkpoints});
  }
  causal::Prescription example_prescription()const {
    // A familiar name cannot authorize a prescription for altered kinematics.
    // Compare the actual incoming coordinates with the native fixture before
    // using its positivity proof or its supplied Henn contour data.
    auto example=feynman::example_family(graph_.family_name);
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
    return causal::current_example(graph_.family_name,graph_.nodes.size(),anchors,standard);
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
  Endpoint endpoint(const Matrix& matrix,const Matrix& diagonal,const Node& node,bool upper,std::size_t depth,int high) {
    auto fuchs=fuchsify::prepare(matrix,xi_,options_.fuchsification);
    if(!fuchs.success)throw std::domain_error("recursive endpoint fuchsification unsupported: "+fuchs.reason);
    const std::string side=upper?"upper":"lower";
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
    auto matching_frame=series.project(options_.observable_adjoint?fuchs.transform:physical_gauge);
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
      auto domain=integral?result.series.dr_domain(result.series.project(transformed),true):
        result.series.project_endpoint_domain(transformed);
      for(const auto& constraint:domain.zero_constraints) {
        if(result.pending_constraints.size()>=options_.max_endpoint_constraints)
          throw std::length_error("recursive endpoint constraint budget exhausted");
        result.pending_constraints.push_back({constant_expansion(Matrix{constraint.coefficients}),
          "level "+std::to_string(depth+1)+" "+side+" operation "+std::to_string(i)+
          " fixed power "+constraint.power.str()+" log degree "+std::to_string(constraint.log_degree)});
      }
      if(integral)result.functionals[i]=result.series.dr_integral_from_zero(domain.admissible);
      else result.functionals[i]=constant_expansion(result.series.dr_endpoint_constant(domain.admissible));
      }catch(const std::length_error& error) {
        throw std::length_error(side+" endpoint operation "+std::to_string(i)+": "+error.what());
      }
      if((i+1)%32==0 || i+1==node.operations.size())
        report(depth,side+" endpoint observables "+std::to_string(i+1)+"/"+std::to_string(node.operations.size()),high);
    }
    return result;
  }
  Plan& plan(std::size_t depth,unsigned desired_top) {
    if(plans_[depth])return *plans_[depth];const auto& node=graph_.nodes[depth];
    if(node.observable_rows.size()!=node.operations.size())throw std::invalid_argument("recursive observable plan shape");
    auto gauge=epsilon_diagonal_gauge(node.closure.matrix,ei_);const auto d=gauge.matrix.size();
    auto diagonal=fuchsify::detail::identity(d,graph_.dimension),inverse=diagonal;
    for(std::size_t i=0;i<d;++i) {
      diagonal[i][i]=graph_.dimension.variable(ei_).pow(gauge.shifts[i]);
      inverse[i][i]=graph_.dimension.constant(1)/diagonal[i][i];
    }
    auto lower=endpoint(gauge.matrix,diagonal,node,false,depth,desired_top),upper=endpoint(reverse(gauge.matrix,true),diagonal,node,true,depth,desired_top);
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
    if(!beta.empty())beta=gauge_observable_columns(gauge,beta);
    auto geometry=endpoint_clearance(lower,upper);
    lower.clearance_coefficients.clear();upper.clearance_coefficients.clear();
    auto largest=*std::max_element(gauge.shifts.begin(),gauge.shifts.end());
    plans_[depth]=std::make_unique<Plan>(Plan{std::move(gauge),std::move(diagonal),std::move(inverse),std::move(lower),std::move(upper),std::move(beta),std::move(indices),largest,std::move(geometry)});
    plans_[depth]->constraint_labels=std::move(labels);
    if(constraints)report(depth,"physical endpoint zero constraints "+std::to_string(constraints),desired_top);
    ++statistics_.exact_plans;
    report(depth,"endpoint overlap="+plans_[depth]->geometry.overlap.str(),desired_top);
    return *plans_[depth];
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
          std::to_string(value.low+static_cast<int>(k)));
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
  LaurentRows endpoint_operator(const Endpoint& endpoint,const B& point,int high)const {
    const auto rows=endpoint.functionals.size();const auto d=endpoint.matching_frame.columns;
    if(std::none_of(endpoint.functionals.begin(),endpoint.functionals.end(),[](const auto& value){return value.has_value();}))
      return {0,high,std::vector(rows,std::vector(d,std::vector<B>(high+1,B(0))))};
    const int inverse_high=constant_demand(endpoint,high);
    auto inverse=affine_operator::prepare(endpoint.series,endpoint.matching_frame,point,inverse_high,options_.matching);
    if(!inverse.success())throw std::domain_error("recursive endpoint operator inverse unsupported: "+inverse.reason);
    std::vector<std::optional<affine_operator::Operator>> operators(rows);int low=0;
    for(unsigned i=0;i<rows;++i)if(endpoint.functionals[i]) {
      operators[i]=affine_operator::compose(inverse,endpoint.series,*endpoint.functionals[i],point,high,options_.matching);
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
    if(!saved || saved->high<high)
      saved=endpoint_operator(upper?prepared.upper:prepared.lower,
          B::from_strings(prepared.geometry.overlap.str()),high);
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
    // The forward API carries physical and integrated maps at a common high:
    // upper map H-r, lower map max(H-l,H-r-q), then the inverse D gauge.
    prepared.factored_child_loss=checked(std::max(loss+std::max(-static_cast<long>(l),-static_cast<long>(r)-q),-static_cast<long>(direct)));
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
        combined=shift_operator_columns(add_laurent_rows(left,right),prepared.gauge.shifts,true);
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
        const int middle_high=std::max(high,checked(static_cast<long>(high)-right.low));
        const int initial_high=std::max(checked(static_cast<long>(high)-left.low),checked(static_cast<long>(middle_high)-beta_low));
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
        auto middle=factored_transport::evolve(prepared.gauge.matrix,lower,beta,path(h,end,depth),middle_high,settings);
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
      for(auto demand:endpoint.series.valuation_metadata(*functional).required_source_top(high))result=std::max(result,demand);
    return checked(result);
  }
  affine_matching::Boundary match(const Endpoint& endpoint,const LaurentBoundary& boundary,const B& point,int high)const {
    auto result=affine_matching::match(endpoint.series,endpoint.matching_frame,point,
      {boundary.low,boundary.high(),boundary.values},{0,high},options_.matching);
    if(result.status==affine_matching::Status::NeedMoreBoundary)throw BoundaryDemand(result.required_boundary_high,result.reason);
    if(!result.success())throw std::domain_error("recursive endpoint matching unsupported: "+result.reason);
    return std::move(result.value);
  }
  LaurentBoundary apply(const Endpoint& endpoint,const Expansion& functional,const affine_matching::Boundary& constants,
      const B& point,unsigned high)const {
    auto result=affine_matching::apply(endpoint.series,functional,point,constants,{0,static_cast<int>(high)},options_.matching);
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
