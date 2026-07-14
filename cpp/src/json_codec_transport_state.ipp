const char* line_integration_scope_name(LineIntegrationScope scope) {
  switch (scope) {
    case LineIntegrationScope::StoredTruncation:
      return "stored_truncation";
    case LineIntegrationScope::FullLocalWithCertifiedTail:
      return "full_local_with_certified_tail";
  }
  throw std::logic_error("unknown line integration scope");
}

const char* error_guarantee_name(ErrorGuarantee guarantee) {
  switch (guarantee) {
    case ErrorGuarantee::None:
      return "none";
    case ErrorGuarantee::Advisory:
      return "advisory";
    case ErrorGuarantee::Certified:
      return "certified";
  }
  throw std::logic_error("unknown error guarantee");
}

json::object encode_error_envelope_summary(const ErrorEnvelope& error) {
  json::object result{
      {"guarantee", error_guarantee_name(error.guarantee)},
      {"provenance", error.provenance}};
  if (!error.empty()) {
    json::array upper;
    upper.reserve(error.absolute.size());
    for (const auto& bound : error.absolute)
      upper.push_back(bound.approximate_upper());
    result["epsilon_min"] = error.frame.min_power;
    result["epsilon_max"] = error.frame.complete_max;
    result["absolute_upper_approx"] = std::move(upper);
    result["bound_encoding"] = "approximate-double-diagnostics";
  }
  return result;
}

struct RetainedPlanChartBinding {
  using Owner = std::variant<std::shared_ptr<PreparedChartBase>,
                             std::shared_ptr<CompositeSCCChartBase>>;

  std::string handle;
  std::string exact_identity;
  ExactAffineChart geometry;
  std::vector<Prescription> prescriptions;
  Owner owner;
};

struct RetainedArmPlan {
  ExactArmPlan exact;
  std::vector<RetainedPlanChartBinding> charts;
};

json::array encode_path_branch_sheets(
    const std::vector<ExactBranchSheet>& sheets) {
  json::array encoded;
  encoded.reserve(sheets.size());
  for (const auto& sheet : sheets)
    encoded.push_back(json::object{{"factor_exact", sheet.factor_exact},
                                   {"sign", sheet.imaginary_sign}});
  return encoded;
}

json::array encode_plan_prescriptions(
    const std::vector<Prescription>& prescriptions) {
  json::array encoded;
  encoded.reserve(prescriptions.size());
  for (const auto& prescription : prescriptions)
    encoded.push_back(json::object{
        {"factor_exact", prescription.factor_exact},
        {"sign", prescription.sign},
        {"multiplicity", prescription.multiplicity},
        {"leading_coefficient_sign",
         prescription.leading_coefficient_sign}});
  return encoded;
}

json::object encode_plan_chart(const RetainedPlanChartBinding& binding,
                               std::size_t index) {
  return json::object{
      {"index", index}, {"chart", binding.handle},
      {"identity", binding.exact_identity},
      {"center_exact", binding.geometry.center.str()},
      {"scale_exact", binding.geometry.scale.str()},
      {"radius_exact", binding.geometry.radius.str()},
      {"singular_center", binding.geometry.singular_center},
      {"prescriptions", encode_plan_prescriptions(binding.prescriptions)}};
}

const char* exact_match_kind_name(ExactMatchKind kind) {
  switch (kind) {
    case ExactMatchKind::SymmetricDivisionPoint:
      return "symmetric-division-point";
    case ExactMatchKind::BalancedSafeOverlap:
      return "balanced-safe-overlap";
    case ExactMatchKind::SingularBalancedApproach:
      return "singular-balanced-approach";
    case ExactMatchKind::ForbiddenPointAvoidance:
      return "forbidden-point-avoidance";
  }
  throw std::logic_error("unknown exact match kind");
}

json::object encode_plan_match(const RetainedArmPlan& arm,
                               std::size_t index) {
  if (index >= arm.exact.matches.size())
    throw std::invalid_argument("native tile-plan match index is out of range");
  const auto& match = arm.exact.matches[index];
  return json::object{
      {"index", index},
      {"producing_chart_index", match.producing_chart},
      {"receiving_chart_index", match.receiving_chart},
      {"producing_chart", arm.charts.at(match.producing_chart).handle},
      {"receiving_chart", arm.charts.at(match.receiving_chart).handle},
      {"physical_exact", match.physical.str()},
      {"producing_local_exact", match.producing_local.str()},
      {"receiving_local_exact", match.receiving_local.str()},
      {"kind", exact_match_kind_name(match.kind)},
      {"branch_sheets", encode_path_branch_sheets(match.branch_sheets)}};
}

json::object encode_plan_tile(const RetainedArmPlan& arm,
                              std::size_t index) {
  if (index >= arm.exact.tiles.size())
    throw std::invalid_argument("native tile-plan tile index is out of range");
  const auto& tile = arm.exact.tiles[index];
  const auto& binding = arm.charts.at(tile.chart);
  return json::object{
      {"index", index}, {"chart_index", tile.chart},
      {"chart", binding.handle}, {"chart_identity", binding.exact_identity},
      {"physical_begin_exact", tile.physical_begin.str()},
      {"physical_end_exact", tile.physical_end.str()},
      {"local_begin_exact", tile.local_begin.str()},
      {"local_end_exact", tile.local_end.str()},
      {"jacobian_exact", binding.geometry.scale.str()},
      {"crosses_singular_center", tile.crosses_singular_center},
      {"branch_sheets", encode_path_branch_sheets(tile.branch_sheets)},
      {"prescriptions", encode_plan_prescriptions(binding.prescriptions)}};
}

json::object encode_plan_topology(const ExactPathTopology& topology) {
  json::array singular_points;
  for (const auto& point : topology.singular_points)
    singular_points.push_back(json::value(point.str()));
  json::array boundary_points;
  for (const auto& point : topology.boundary_points)
    boundary_points.push_back(json::value(point.str()));
  json::array projections;
  for (const auto& projection : topology.complex_projections)
    projections.push_back(json::object{
        {"source_identity", projection.source_identity},
        {"real_part_exact", projection.real_part.str()},
        {"imaginary_magnitude_exact", projection.imaginary_magnitude.str()},
        {"retain_minus_imaginary", projection.retain_minus_imaginary},
        {"retain_real_part", projection.retain_real_part},
        {"retain_plus_imaginary", projection.retain_plus_imaginary}});
  return json::object{
      {"singular_points", std::move(singular_points)},
      {"boundary_points", std::move(boundary_points)},
      {"complex_projections", std::move(projections)},
      {"branch_sheets", encode_path_branch_sheets(topology.branch_sheets)}};
}

json::object encode_retained_arm(const RetainedArmPlan& arm) {
  json::array charts;
  for (std::size_t index = 0; index < arm.charts.size(); ++index)
    charts.push_back(encode_plan_chart(arm.charts[index], index));
  json::array matches;
  for (std::size_t index = 0; index < arm.exact.matches.size(); ++index)
    matches.push_back(encode_plan_match(arm, index));
  json::array tiles;
  for (std::size_t index = 0; index < arm.exact.tiles.size(); ++index)
    tiles.push_back(encode_plan_tile(arm, index));
  return json::object{
      {"from_exact", arm.exact.from.str()},
      {"to_exact", arm.exact.to.str()},
      {"direction", arm.exact.direction},
      {"division_order", arm.exact.division_order},
      {"charts", std::move(charts)}, {"matches", std::move(matches)},
      {"tiles", std::move(tiles)},
      {"topology", encode_plan_topology(arm.exact.topology)}};
}

std::optional<std::int32_t> exact_plan_rim(
    const std::vector<Prescription>& prescriptions,
    const Rational& chart_scale) {
  if (chart_scale.is_zero())
    throw std::invalid_argument(
        "prepared tile chart has a zero exact scale");
  const auto scale_sign = chart_scale.sign();
  std::optional<std::int32_t> rim;
  for (const auto& prescription : prescriptions) {
    if ((prescription.multiplicity & 1U) == 0)
      throw std::invalid_argument(
          "prepared tile chart has an even-multiplicity tangential "
          "prescription; a one-sided real-axis rim is not defined");
    const auto candidate =
        prescription.sign * prescription.leading_coefficient_sign *
        scale_sign;
    if (rim.has_value() && *rim != candidate)
      throw std::invalid_argument(
          "prepared tile chart has conflicting exact odd-multiplicity prescriptions");
    rim = candidate;
  }
  return rim;
}

class StoredTilePlan {
 public:
  StoredTilePlan(std::string handle, std::string checkpoint_identity,
                 std::string provenance_identity, std::uint32_t division_order,
                 RetainedArmPlan lower, RetainedArmPlan upper,
                 double elapsed_ms)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        division_order_(division_order), lower_(std::move(lower)),
        upper_(std::move(upper)), elapsed_ms_(elapsed_ms) {
    validate_arm_set();
  }

  StoredTilePlan(std::string handle, std::string checkpoint_identity,
                 std::string provenance_identity, std::uint32_t division_order,
                 std::string arm_name, RetainedArmPlan arm,
                 double elapsed_ms)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        division_order_(division_order), elapsed_ms_(elapsed_ms) {
    if (arm_name == "lower")
      lower_ = std::move(arm);
    else if (arm_name == "upper")
      upper_ = std::move(arm);
    else
      throw std::invalid_argument(
          "single-arm tile plan name must be lower or upper");
    validate_arm_set();
  }

  const std::string& handle() const { return handle_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& provenance_identity() const {
    return provenance_identity_;
  }

  bool has_arm(const std::string& name) const {
    if (name == "lower") return lower_.has_value();
    if (name == "upper") return upper_.has_value();
    return false;
  }

  bool has_two_arms() const {
    return lower_.has_value() && upper_.has_value();
  }

  const RetainedArmPlan& arm(const std::string& name) const {
    if (name == "lower" && lower_.has_value()) return *lower_;
    if (name == "upper" && upper_.has_value()) return *upper_;
    if (name != "lower" && name != "upper")
      throw std::invalid_argument(
          "native tile-plan arm must be lower or upper");
    throw std::invalid_argument(
        "native tile plan does not retain the requested " + name + " arm");
  }

  json::object match_interval(const std::string& name,
                              std::size_t index) const {
    match_queries_.fetch_add(1);
    auto encoded = encode_plan_match(arm(name), index);
    encoded["arm"] = name;
    encoded["tile_plan"] = handle_;
    encoded["checkpoint_identity"] = checkpoint_identity_;
    return encoded;
  }

  json::object tile_interval(const std::string& name,
                             std::size_t index) const {
    tile_queries_.fetch_add(1);
    auto encoded = encode_plan_tile(arm(name), index);
    encoded["arm"] = name;
    encoded["tile_plan"] = handle_;
    encoded["checkpoint_identity"] = checkpoint_identity_;
    return encoded;
  }

  void note_integration() { integrations_.fetch_add(1); }

  void note_match_advance(const std::string& name) {
    (void)arm(name);
    if (name == "lower") {
      lower_match_advances_.fetch_add(1);
      return;
    }
    if (name == "upper") {
      upper_match_advances_.fetch_add(1);
      return;
    }
    throw std::invalid_argument("native tile-plan arm must be lower or upper");
  }

  void require_two_arm_match_advance_capacity(
      std::size_t lower, std::size_t upper) const {
    if (!has_two_arms())
      throw std::invalid_argument(
          "parallel match accounting requires a retained two-arm plan");
    const auto fits = [](std::uint64_t current, std::size_t increment) {
      return increment <= std::numeric_limits<std::uint64_t>::max() - current;
    };
    if (!fits(lower_match_advances_.load(), lower) ||
        !fits(upper_match_advances_.load(), upper))
      throw std::overflow_error(
          "retained tile-plan match-advance counter overflow");
  }

  void note_two_arm_match_advances(
      std::size_t lower, std::size_t upper) noexcept {
    lower_match_advances_.fetch_add(static_cast<std::uint64_t>(lower));
    upper_match_advances_.fetch_add(static_cast<std::uint64_t>(upper));
  }

  std::vector<std::shared_ptr<PreparedChartBase>> dependency_charts() const {
    std::vector<std::shared_ptr<PreparedChartBase>> result;
    result.reserve((lower_.has_value() ? lower_->charts.size() : 0) +
                   (upper_.has_value() ? upper_->charts.size() : 0));
    const auto append = [&](const RetainedPlanChartBinding& binding) {
      std::visit(
          [&](const auto& owner) {
            using Owner = typename std::decay_t<decltype(owner)>::element_type;
            if constexpr (std::is_same_v<Owner, PreparedChartBase>) {
              result.push_back(owner);
            } else {
              auto dependencies = owner->dependency_charts();
              result.insert(result.end(), dependencies.begin(),
                            dependencies.end());
            }
          },
          binding.owner);
    };
    if (lower_.has_value())
      for (const auto& binding : lower_->charts) append(binding);
    if (upper_.has_value())
      for (const auto& binding : upper_->charts) append(binding);
    return result;
  }

  std::vector<std::shared_ptr<CompositeSCCChartBase>> dependency_sccs() const {
    std::vector<std::shared_ptr<CompositeSCCChartBase>> result;
    const auto append = [&](const RetainedPlanChartBinding& binding) {
      if (const auto* owner = std::get_if<
              std::shared_ptr<CompositeSCCChartBase>>(&binding.owner))
        result.push_back(*owner);
    };
    if (lower_.has_value())
      for (const auto& binding : lower_->charts) append(binding);
    if (upper_.has_value())
      for (const auto& binding : upper_->charts) append(binding);
    return result;
  }

  json::object summary(bool include_intervals = true) const {
    if (!has_two_arms()) return single_arm_summary(include_intervals);
    json::object result{
        {"tile_plan", handle_}, {"capability", kRetainedTilePlanCapability},
        {"native_retained", true}, {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"division_order", division_order_},
        {"independent_arms", true},
        {"concurrent_execution", "immutable-independent-arm-snapshots"},
        {"anchor", json::object{
             {"lower_chart", lower_->charts.front().handle},
             {"upper_chart", upper_->charts.front().handle},
             {"center_exact", lower_->exact.from.str()}}},
        {"lower_matches", lower_->exact.matches.size()},
        {"upper_matches", upper_->exact.matches.size()},
        {"lower_tiles", lower_->exact.tiles.size()},
        {"upper_tiles", upper_->exact.tiles.size()},
        {"match_interval_queries", match_queries_.load()},
        {"tile_interval_queries", tile_queries_.load()},
        {"lower_match_advances", lower_match_advances_.load()},
        {"upper_match_advances", upper_match_advances_.load()},
        {"integrations", integrations_.load()},
        {"elapsed_ms", elapsed_ms_}};
    if (include_intervals) {
      result["lower"] = encode_retained_arm(*lower_);
      result["upper"] = encode_retained_arm(*upper_);
    }
    return result;
  }

  json::object checkpoint_record() const {
    if (!has_two_arms()) {
      const auto name = single_arm_name();
      return json::object{
          {"schema", kRetainedSingleArmTilePlanCheckpointSchema},
          {"handle", handle_},
          {"checkpoint_identity", checkpoint_identity_},
          {"provenance_identity", provenance_identity_},
          {"division_order", division_order_},
          {"arm_name", name},
          {"arm", encode_retained_arm(arm(name))},
          {"elapsed_ms", elapsed_ms_},
          {"runtime_stats", runtime_stats_record()}};
    }
    return json::object{
        {"schema", "diffexp2-retained-tile-plan-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"division_order", division_order_},
        {"lower", encode_retained_arm(*lower_)},
        {"upper", encode_retained_arm(*upper_)},
        {"elapsed_ms", elapsed_ms_},
        {"runtime_stats", runtime_stats_record()}};
  }

  void restore_runtime_stats(std::uint64_t match_queries,
                             std::uint64_t tile_queries,
                             std::uint64_t lower_match_advances,
                             std::uint64_t upper_match_advances,
                             std::uint64_t integrations) {
    if ((!lower_.has_value() && lower_match_advances != 0) ||
        (!upper_.has_value() && upper_match_advances != 0))
      throw std::invalid_argument(
          "single-arm tile-plan checkpoint advances an absent arm");
    match_queries_.store(match_queries);
    tile_queries_.store(tile_queries);
    lower_match_advances_.store(lower_match_advances);
    upper_match_advances_.store(upper_match_advances);
    integrations_.store(integrations);
  }

 private:
  void validate_arm_set() const {
    if (!lower_.has_value() && !upper_.has_value())
      throw std::invalid_argument(
          "retained tile plan must own one or two arms");
    if (lower_.has_value()) validate_exact_arm_plan(lower_->exact);
    if (upper_.has_value()) validate_exact_arm_plan(upper_->exact);
    // Existing two-arm requests name the slots but historically only require
    // opposite directions; preserve that behavior exactly. A genuine
    // single-arm plan derives its retained name from the exact direction.
    if (!has_two_arms()) {
      const auto& retained = lower_.has_value() ? *lower_ : *upper_;
      const auto expected_direction = lower_.has_value() ? -1 : 1;
      if (retained.exact.direction != expected_direction)
        throw std::invalid_argument(
            "retained single tile-arm name differs from its exact direction");
    }
  }

  std::string single_arm_name() const {
    if (has_two_arms())
      throw std::logic_error(
          "two-arm tile plan has no single arm name");
    return lower_.has_value() ? "lower" : "upper";
  }

  json::object runtime_stats_record() const {
    return json::object{
        {"match_interval_queries", match_queries_.load()},
        {"tile_interval_queries", tile_queries_.load()},
        {"lower_match_advances", lower_match_advances_.load()},
        {"upper_match_advances", upper_match_advances_.load()},
        {"integrations", integrations_.load()}};
  }

  json::object single_arm_summary(bool include_intervals) const {
    const auto name = single_arm_name();
    const auto& retained = arm(name);
    json::object result{
        {"tile_plan", handle_},
        {"capability", kRetainedSingleArmTilePlanCapability},
        {"native_retained", true},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"division_order", division_order_},
        {"independent_arms", false},
        {"concurrent_execution", "single-immutable-arm-snapshot"},
        {"arm_name", name},
        {"anchor", json::object{
             {"chart", retained.charts.front().handle},
             {"center_exact", retained.exact.from.str()}}},
        {"matches", retained.exact.matches.size()},
        {"tiles", retained.exact.tiles.size()},
        {"match_interval_queries", match_queries_.load()},
        {"tile_interval_queries", tile_queries_.load()},
        {"lower_match_advances", lower_match_advances_.load()},
        {"upper_match_advances", upper_match_advances_.load()},
        {"integrations", integrations_.load()},
        {"elapsed_ms", elapsed_ms_}};
    if (include_intervals) {
      result["arm"] = encode_retained_arm(retained);
      result[name] = encode_retained_arm(retained);
    }
    return result;
  }

  std::string handle_;
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::uint32_t division_order_ = 3;
  std::optional<RetainedArmPlan> lower_;
  std::optional<RetainedArmPlan> upper_;
  double elapsed_ms_ = 0.0;
  mutable std::atomic<std::uint64_t> match_queries_{0};
  mutable std::atomic<std::uint64_t> tile_queries_{0};
  std::atomic<std::uint64_t> lower_match_advances_{0};
  std::atomic<std::uint64_t> upper_match_advances_{0};
  std::atomic<std::uint64_t> integrations_{0};
};

// A plan-driven match is one exact handoff, not a completed arm.  It owns the
// immutable plan snapshot and every local used to construct the retained
// matching weights.  Registry release of the public plan/local tokens can
// therefore never turn a published handoff into dangling provenance.
class StoredPlannedMatchHop final : public StoredMatchBase {
 public:
  StoredPlannedMatchHop(
      std::shared_ptr<StoredMatchBase> match,
      std::string checkpoint_identity, std::string provenance_identity,
      json::object handoff, double elapsed_ms,
      std::shared_ptr<StoredTilePlan> plan_owner,
      std::vector<std::shared_ptr<StoredLocalBase>> basis_owners,
      std::shared_ptr<StoredLocalBase> incoming_owner)
      : StoredMatchBase(match == nullptr ? std::string() : match->handle()),
        match_(std::move(match)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        handoff_(std::move(handoff)), elapsed_ms_(elapsed_ms),
        plan_owner_(std::move(plan_owner)),
        basis_owners_(std::move(basis_owners)),
        incoming_owner_(std::move(incoming_owner)) {
    if (match_ == nullptr || plan_owner_ == nullptr ||
        incoming_owner_ == nullptr || basis_owners_.empty())
      throw std::invalid_argument(
          "retained planned match hop requires all strong owners");
  }

  json::object summary() const override {
    auto result = match_->summary();
    if (required_string(result, "checkpoint_identity") !=
        checkpoint_identity_)
      throw std::logic_error(
          "retained planned match checkpoint identity changed");
    result["planned_hop_capability"] =
        kRetainedPlannedMatchHopCapability;
    result["plan_driven"] = true;
    result["planned_hop_provenance_identity"] = provenance_identity_;
    result["planned_hop"] = handoff_;
    result["strong_ownership"] = json::object{
        {"tile_plan", true}, {"basis_locals", basis_owners_.size()},
        {"incoming_local", true}};
    result["materializations"] = materializations_.load();
    result["elapsed_ms"] = elapsed_ms_;
    return result;
  }

  double elapsed_ms() const { return elapsed_ms_; }

  const std::shared_ptr<StoredMatchBase>& native_match() const {
    return match_;
  }
  const std::shared_ptr<StoredTilePlan>& plan_owner() const {
    return plan_owner_;
  }
  const std::vector<std::shared_ptr<StoredLocalBase>>& basis_owners() const {
    return basis_owners_;
  }
  const std::shared_ptr<StoredLocalBase>& incoming_owner() const {
    return incoming_owner_;
  }
  const json::object& handoff() const { return handoff_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& provenance_identity() const {
    return provenance_identity_;
  }

  std::optional<json::object> incomplete_acb_summary() const {
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    if (!acb || acb->certified_for_materialization()) return std::nullopt;
    return acb->summary();
  }

  std::shared_ptr<PhysicalEquationOwnerBase>
  inheritable_basis_equation_owner() const {
    if (basis_owners_.empty()) return nullptr;
    auto result = basis_owners_.front()->retained_equation_owner();
    for (const auto& basis : basis_owners_) {
      const auto candidate = basis->retained_equation_owner();
      if (candidate.get() != result.get()) return nullptr;
    }
    return result;
  }

  void validate_materialized_equation_owner(
      const std::shared_ptr<StoredLocalBase>& local) const {
    if (!local || local->retained_equation_owner().get() !=
                      inheritable_basis_equation_owner().get())
      throw std::invalid_argument(
          "checkpoint materialized-local equation owner differs from its planned-hop basis");
  }

  void validate_materialized_derivation(
      const json::object& derivation, const char* scalar_domain) const {
    const auto native_summary = match_->summary();
    if (required_string(derivation, "source_match") != handle_ ||
        required_string(derivation, "source_match_checkpoint_identity") !=
            checkpoint_identity_ ||
        required_string(derivation, "source_match_provenance_identity") !=
            required_string(native_summary, "provenance_identity") ||
        required_string(derivation, "planned_hop_provenance_identity") !=
            provenance_identity_ ||
        derivation.at("planned_hop") != handoff_)
      throw std::invalid_argument(
          "checkpoint materialized-local lineage disagrees with its planned-hop owner");

    json::array expected_windows;
    std::int32_t expected_certified_max = 0;
    if (const auto exact =
            std::dynamic_pointer_cast<StoredExactRegularMatch>(match_)) {
      if (std::string(scalar_domain) != "rational")
        throw std::invalid_argument(
            "checkpoint materialized-local scalar domain differs from its exact match owner");
      for (const auto& weight : exact->weights())
        expected_windows.push_back(json::object{
            {"min", weight.min_power()}, {"max", weight.complete_max()}});
      expected_certified_max = as_i32(
          as_object(native_summary.at("residual"),
                    "exact retained match residual").at("max"),
          "exact retained match residual maximum");
    } else if (const auto acb =
                   std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_)) {
      if (std::string(scalar_domain) != "acb" ||
          !acb->certified_for_materialization())
        throw std::invalid_argument(
            "checkpoint materialized-local Acb owner lost its passing complete match certificate");
      for (const auto& weight : acb->weights())
        expected_windows.push_back(json::object{
            {"min", weight.min_power()}, {"max", weight.complete_max()}});
      expected_certified_max = acb->certified_complete_max();
    } else {
      throw std::invalid_argument(
          "checkpoint materialized-local owner embeds an unsupported native match");
    }
    if (derivation.at("weight_windows") != expected_windows ||
        as_i32(derivation.at("match_certified_complete_max"),
               "materialized-local certified maximum") !=
            expected_certified_max)
      throw std::invalid_argument(
          "checkpoint materialized-local derivation differs from its retained match weights/certificate");
  }

  std::shared_ptr<StoredLocalBase> materialize(
      const std::string& local_handle,
      const std::string& result_checkpoint_identity,
      slong precision_bits,
      const std::shared_ptr<StoredPlannedMatchHop>& self) {
    if (self.get() != this)
      throw std::logic_error(
          "retained plan-match materialization lost self ownership");
    std::shared_ptr<StoredLocalBase> result;
    if (const auto exact =
            std::dynamic_pointer_cast<StoredExactRegularMatch>(match_)) {
      result = materialize_typed<Rational>(
          local_handle, result_checkpoint_identity, precision_bits,
          exact->weights(), self);
    } else if (const auto acb =
                   std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_)) {
      if (!acb->certified_for_materialization())
        {
          const auto summary = acb->summary();
          const auto& lattice = as_object(
              summary.at("exact_lattice"),
              "failed Acb match exact lattice summary");
          const json::object compact_lattice{
              {"normalized_determinant_valuation",
               lattice.at("normalized_determinant_valuation")},
              {"transformation_min_power",
               lattice.at("transformation_min_power")},
              {"transformation_terms",
               lattice.at("transformation_terms")},
              {"initial_column_shifts",
               lattice.at("initial_column_shifts")},
              {"initial_leading_rank",
               lattice.at("initial_leading_rank")},
              {"final_leading_rank",
               lattice.at("final_leading_rank")}};
          throw std::domain_error(
              "an Acb plan-match handoff must have a passing complete residual before materialization; residual=" +
              json::serialize(summary.at("residual")) +
              "; refinement=" +
              json::serialize(summary.at("refinement")) +
              "; epsilon=" +
              json::serialize(summary.at("epsilon")) +
              "; weight_windows=" +
              json::serialize(summary.at("weight_windows")) +
              "; exact_lattice=" +
              json::serialize(compact_lattice));
        }
      result = materialize_typed<ComplexBall>(
          local_handle, result_checkpoint_identity, precision_bits,
          acb->weights(), self);
    } else {
      throw std::logic_error(
          "retained plan-match handoff has an unsupported matching state");
    }
    materializations_.fetch_add(1);
    return result;
  }

  json::object checkpoint_record() const override {
    return json::object{
        {"schema", "diffexp2-retained-planned-match-hop-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"handoff", handoff_},
        {"native_match", match_->checkpoint_record()},
        {"elapsed_ms", elapsed_ms_},
        {"runtime_stats",
         json::object{{"materializations", materializations_.load()}}}};
  }

  void restore_runtime_stats(std::uint64_t materializations) {
    materializations_.store(materializations);
  }

 private:
  template <typename Scalar>
  std::shared_ptr<StoredLocalBase> materialize_typed(
      const std::string& local_handle,
      const std::string& result_checkpoint_identity,
      slong precision_bits,
      const FiniteLaurentVector<Scalar>& weights,
      const std::shared_ptr<StoredPlannedMatchHop>& self) const {
    if (local_handle.empty() || result_checkpoint_identity.empty())
      throw std::invalid_argument(
          "plan-match local materialization identities must be nonempty");
    AcbPrecisionLease lease(precision_bits);
    ComplexBall::set_precision(precision_bits);
    const auto started = std::chrono::steady_clock::now();
    const auto native_match_summary = match_->summary();
    const auto& match_epsilon = as_object(
        native_match_summary.at("epsilon"),
        "retained match epsilon provenance");
    const auto match_work_max = as_i32(
        match_epsilon.at("max"), "retained match epsilon maximum");
    const auto match_certified_max = [&]() {
      if constexpr (std::is_same_v<Scalar, Rational>)
        return as_i32(
            as_object(native_match_summary.at("residual"),
                      "exact retained match residual").at("max"),
            "exact retained match residual maximum");
      else {
        const auto acb =
            std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
        if (!acb)
          throw std::logic_error(
              "Acb materialization lost its retained match owner");
        return acb->certified_complete_max();
      }
    }();

    std::vector<std::shared_ptr<StoredLocal<Scalar>>> typed_basis;
    std::vector<const LocalSolution<Scalar>*> basis_solutions;
    typed_basis.reserve(basis_owners_.size());
    basis_solutions.reserve(basis_owners_.size());
    std::int32_t materialized_top = kCompleteInfinity;
    const auto receiving_chart = basis_owners_.front()->source_chart();
    const auto receiving_operator =
        basis_owners_.front()->source_operator_identity();
    auto equation_owner = inheritable_basis_equation_owner();
    for (std::size_t column = 0; column < basis_owners_.size(); ++column) {
      auto typed =
          std::dynamic_pointer_cast<StoredLocal<Scalar>>(basis_owners_[column]);
      if (!typed)
        throw std::logic_error(
            "retained plan-match basis coefficient domain changed");
      if (typed->source_chart() != receiving_chart ||
          typed->source_operator_identity() != receiving_operator)
        throw std::logic_error(
            "retained plan-match basis chart provenance changed");
      if (typed->top_valid() < match_work_max)
        throw std::domain_error(
            "retained plan-match basis validity does not cover its matching work window");
      if (column >= weights.size())
        throw std::logic_error(
            "retained plan-match weight count is smaller than its basis");
      const auto basis_valid = std::min(
          typed->top_valid(), typed->solution().epsilon.complete_max);
      const auto shifted_basis_valid = local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(basis_valid) +
              weights[column].min_power(),
          "materialized local top validity");
      const auto weight_valid = local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(typed->solution().epsilon.min_power) +
              weights[column].complete_max(),
          "materialized weight top validity");
      materialized_top = std::min(
          materialized_top, std::min(shifted_basis_valid, weight_valid));
      basis_solutions.push_back(&typed->solution());
      typed_basis.push_back(std::move(typed));
    }
    if (weights.size() != typed_basis.size())
      throw std::logic_error(
          "retained plan-match weight count differs from its basis");
    const auto typed_incoming =
        std::dynamic_pointer_cast<StoredLocal<Scalar>>(incoming_owner_);
    if (!typed_incoming || typed_incoming->top_valid() < match_work_max)
      throw std::domain_error(
          "retained plan-match incoming validity does not cover its matching work window");

    auto solution = materialize_local_basis_weights(
        basis_solutions, weights, result_checkpoint_identity);
    materialized_top = std::min(
        {materialized_top, solution.epsilon.complete_max,
         match_certified_max});
    if (materialized_top < solution.epsilon.min_power)
      throw std::domain_error(
          "plan-match materialization has no valid output epsilon coefficient");
    if (materialized_top < solution.epsilon.complete_max)
      solution = restrict_local_epsilon_frame_strict_lower(
          solution, solution.epsilon.min_power, materialized_top,
          result_checkpoint_identity);
    json::array weight_windows;
    weight_windows.reserve(weights.size());
    for (const auto& weight : weights)
      weight_windows.push_back(json::object{
          {"min", weight.min_power()},
          {"max", weight.complete_max()}});
    json::object output_record{
        {"checkpoint_identity", result_checkpoint_identity},
        {"chart", receiving_chart},
        {"source_operator_identity", receiving_operator},
        {"epsilon", json::object{
             {"min", solution.epsilon.min_power},
             {"max", solution.epsilon.complete_max}}},
        {"taylor_complete_max", solution.taylor_complete_max},
        {"dimension", solution.dimension}};
    const bool compact_derivation = required_string(handoff_, "schema") ==
        "diffexp2-retained-exact-plan-match-hop-v2";
    json::object derivation;
    if (compact_derivation) {
      if (!plan_owner_ || equation_owner == nullptr ||
          required_string(handoff_, "tile_plan") != plan_owner_->handle() ||
          required_string(handoff_, "tile_plan_checkpoint_identity") !=
              plan_owner_->checkpoint_identity() ||
          required_string(handoff_, "result_checkpoint_identity") !=
              checkpoint_identity_)
        throw std::logic_error(
            "compact plan-match derivation lost its live plan, match, or equation owner");
      const auto& producing = as_object(
          handoff_.at("producing"), "compact plan-match producing record");
      const auto& receiving = as_object(
          handoff_.at("receiving"), "compact plan-match receiving record");
      derivation = json::object{
          {"schema", "diffexp2-retained-plan-match-local-materialization-v2"},
          {"capability", kRetainedPlannedMatchMaterializationCapability},
          {"source_match", handle_},
          {"source_match_checkpoint_identity", checkpoint_identity_},
          {"tile_plan", handoff_.at("tile_plan")},
          {"tile_plan_checkpoint_identity",
           handoff_.at("tile_plan_checkpoint_identity")},
          {"arm", handoff_.at("arm")},
          {"match", handoff_.at("match")},
          {"incoming", producing.at("incoming")},
          {"basis", receiving.at("basis")},
          {"weight_windows", std::move(weight_windows)},
          {"match_certified_complete_max", match_certified_max},
          {"output", std::move(output_record)},
          {"equation_owner_signature_identity",
           equation_owner->owner_signature_identity()},
          {"equation_payload_identity",
           equation_owner->physical_payload_identity()},
          {"scope", "single-match-receiving-local"},
          {"coefficient_transport", "native-retained-only"},
          {"whole_arm_complete", false}};
    } else {
      derivation = json::object{
          {"schema", "diffexp2-retained-plan-match-local-materialization-v1"},
          {"capability", kRetainedPlannedMatchMaterializationCapability},
          {"source_match", handle_},
          {"source_match_checkpoint_identity", checkpoint_identity_},
          {"source_match_provenance_identity",
           required_string(native_match_summary, "provenance_identity")},
          {"planned_hop_provenance_identity", provenance_identity_},
          {"planned_hop", handoff_},
          {"weight_windows", std::move(weight_windows)},
          {"match_certified_complete_max", match_certified_max},
          {"output", std::move(output_record)},
          {"scope", "single-match-receiving-local"},
          {"coefficient_transport", "native-retained-only"},
          {"whole_arm_complete", false}};
    }
    const auto derivation_identity = json::serialize(
        canonical_json_value(derivation));
    derivation["provenance_identity"] = derivation_identity;
    std::vector<const RegularTaylorTailModelResult*> basis_tail_models;
    basis_tail_models.reserve(typed_basis.size());
    for (const auto& column : typed_basis)
      basis_tail_models.push_back(&column->tail_model());
    auto tail_model = derive_materialized_regular_homogeneous_tail_model(
        basis_solutions, basis_tail_models, weights, solution,
        receiving_operator, checkpoint_identity_);
    std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>>
        physical_equation;
    if (equation_owner != nullptr) {
      if (equation_owner->equation_owner_handle() != receiving_chart ||
          equation_owner->equation_operator_identity() !=
              receiving_operator ||
          std::string(equation_owner->equation_scalar_domain()) !=
              (std::is_same_v<Scalar, Rational> ? "rational" : "acb"))
        throw std::logic_error(
            "retained plan-match basis equation owner disagrees with its receiving chart");
      const auto erased = equation_owner->physical_ode_erased();
      if (!erased)
        throw std::logic_error(
            "retained plan-match basis equation owner lost its physical q/C payload");
      physical_equation =
          std::static_pointer_cast<const PreparedPhysicalClearedODE<Scalar>>(
              erased);
      physical_ode_detail::validate_ode(*physical_equation);
      if (physical_equation->exact_payload_record !=
              equation_owner->physical_payload_record() ||
          physical_equation->payload_identity !=
              equation_owner->physical_payload_identity() ||
          physical_equation->owner_signature_identity !=
              equation_owner->owner_signature_identity())
        throw std::logic_error(
            "retained plan-match physical q/C payload differs from its shared basis owner");
    }
    const auto elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    NativeLocalDiagnostics diagnostics;
    diagnostics.top_valid = materialized_top;
    diagnostics.kernel_ms = elapsed_ms;
    return make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
        local_handle, receiving_chart, receiving_operator,
        std::move(solution), precision_bits,
        std::vector<PseudoHit<Scalar>>{}, diagnostics, std::nullopt,
        std::move(derivation), std::static_pointer_cast<void>(self),
        std::move(tail_model), std::nullopt, true, true,
        std::move(equation_owner), std::move(physical_equation));
  }

  std::shared_ptr<StoredMatchBase> match_;
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  json::object handoff_;
  double elapsed_ms_ = 0.0;
  std::shared_ptr<StoredTilePlan> plan_owner_;
  std::vector<std::shared_ptr<StoredLocalBase>> basis_owners_;
  std::shared_ptr<StoredLocalBase> incoming_owner_;
  std::atomic<std::uint64_t> materializations_{0};
};

class StoredTransportArmState {
 public:
  StoredTransportArmState(
      std::string handle, std::string checkpoint_identity,
      std::string arm, std::shared_ptr<StoredTilePlan> plan_owner,
      std::shared_ptr<StoredLocalBase> anchor_owner,
      std::vector<std::vector<std::shared_ptr<StoredLocalBase>>> basis_owners,
      std::vector<std::shared_ptr<StoredPlannedMatchHop>> matches,
      std::vector<std::shared_ptr<StoredLocalBase>> tile_sources,
      EpsilonWindow work_epsilon,
      std::int32_t public_required_complete_max,
      std::int32_t match_required_complete_max,
      json::object refinement, double elapsed_ms,
      bool compact_provenance = true)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        arm_(std::move(arm)), plan_owner_(std::move(plan_owner)),
        anchor_owner_(std::move(anchor_owner)),
        basis_owners_(std::move(basis_owners)),
        matches_(std::move(matches)),
        tile_sources_(std::move(tile_sources)),
        work_epsilon_(work_epsilon),
        public_required_complete_max_(public_required_complete_max),
        match_required_complete_max_(match_required_complete_max),
        refinement_(std::move(refinement)), elapsed_ms_(elapsed_ms),
        compact_provenance_(compact_provenance) {
    validate();
    provenance_identity_ = json::serialize(
        canonical_json_value(provenance_record()));
  }

  // Ownership-taking streamed publication already has the complete ordered
  // chain of materialized tile locals.  Each non-anchor local carries a
  // sealed v2 plan-match derivation that was validated before its (large)
  // basis and native match owners were released.  Rebind that internal
  // certificate directly instead of round-tripping basis/match provenance
  // through Wolfram and copying it into the transport state.
  StoredTransportArmState(
      std::string handle, std::string checkpoint_identity,
      std::string arm, std::shared_ptr<StoredTilePlan> plan_owner,
      std::shared_ptr<StoredLocalBase> anchor_owner,
      std::vector<std::shared_ptr<StoredLocalBase>> tile_sources,
      EpsilonWindow work_epsilon,
      std::int32_t public_required_complete_max,
      std::int32_t match_required_complete_max,
      json::object refinement, double elapsed_ms)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        arm_(std::move(arm)), plan_owner_(std::move(plan_owner)),
        anchor_owner_(std::move(anchor_owner)),
        tile_sources_(std::move(tile_sources)),
        work_epsilon_(work_epsilon),
        public_required_complete_max_(public_required_complete_max),
        match_required_complete_max_(match_required_complete_max),
        refinement_(std::move(refinement)), elapsed_ms_(elapsed_ms),
        compact_provenance_(true), consumed_compact_(true),
        consumed_certificate_only_(true) {
    validate();
    provenance_identity_ = json::serialize(
        canonical_json_value(provenance_record()));
  }

  StoredTransportArmState(
      std::string handle, std::string checkpoint_identity,
      std::string arm, std::shared_ptr<StoredTilePlan> plan_owner,
      std::shared_ptr<StoredLocalBase> anchor_owner,
      json::array basis_references, json::array match_references,
      std::vector<std::shared_ptr<StoredLocalBase>> tile_sources,
      EpsilonWindow work_epsilon,
      std::int32_t public_required_complete_max,
      std::int32_t match_required_complete_max,
      json::object refinement, double elapsed_ms)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        arm_(std::move(arm)), plan_owner_(std::move(plan_owner)),
        anchor_owner_(std::move(anchor_owner)),
        tile_sources_(std::move(tile_sources)),
        work_epsilon_(work_epsilon),
        public_required_complete_max_(public_required_complete_max),
        match_required_complete_max_(match_required_complete_max),
        refinement_(std::move(refinement)), elapsed_ms_(elapsed_ms),
        compact_provenance_(true), consumed_compact_(true),
        cached_basis_references_(std::move(basis_references)),
        cached_match_references_(std::move(match_references)) {
    validate();
    provenance_identity_ = json::serialize(
        canonical_json_value(provenance_record()));
  }

  const std::string& handle() const { return handle_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& provenance_identity() const {
    return provenance_identity_;
  }
  const std::string& arm_name() const { return arm_; }
  double elapsed_ms() const { return elapsed_ms_; }
  const std::shared_ptr<StoredTilePlan>& plan_owner() const {
    return plan_owner_;
  }
  const std::shared_ptr<StoredLocalBase>& anchor_owner() const {
    return anchor_owner_;
  }
  const std::vector<std::vector<std::shared_ptr<StoredLocalBase>>>&
  basis_owners() const {
    return basis_owners_;
  }
  const std::vector<std::shared_ptr<StoredPlannedMatchHop>>& matches() const {
    return matches_;
  }
  const std::vector<std::shared_ptr<StoredLocalBase>>& tile_sources() const {
    return tile_sources_;
  }
  const std::shared_ptr<StoredLocalBase>& final_local() const {
    return tile_sources_.back();
  }
  std::int32_t public_required_complete_max() const {
    return public_required_complete_max_;
  }

  void require_contraction_counter_capacity(
      std::size_t observables) const {
    const auto operation_count = contraction_operations_.load();
    const auto observable_count = contracted_observables_.load();
    if (operation_count == std::numeric_limits<std::uint64_t>::max() ||
        observables > std::numeric_limits<std::uint64_t>::max() -
                          observable_count)
      throw std::overflow_error(
          "retained transport-arm contraction counter overflow");
  }

  void note_contraction_success(std::size_t observables) noexcept {
    contraction_operations_.fetch_add(1);
    contracted_observables_.fetch_add(
        static_cast<std::uint64_t>(observables));
  }

  void require_endpoint_batch_counter_capacity(std::size_t rows) const {
    const auto operation_count = endpoint_batch_operations_.load();
    const auto row_count = endpoint_rows_.load();
    if (operation_count == std::numeric_limits<std::uint64_t>::max() ||
        rows > std::numeric_limits<std::uint64_t>::max() - row_count)
      throw std::overflow_error(
          "retained transport endpoint-batch counter overflow");
  }

  void note_endpoint_batch_success(std::size_t rows) noexcept {
    endpoint_batch_operations_.fetch_add(1);
    endpoint_rows_.fetch_add(static_cast<std::uint64_t>(rows));
  }

  json::object summary() const {
    return json::object{
        {"transport_state", handle_},
        {"capability", kRetainedTransportArmStateCapability},
        {"native_retained", true},
        {"json_coefficients", 0},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"tile_plan", plan_owner_->handle()},
        {"tile_plan_checkpoint_identity",
         plan_owner_->checkpoint_identity()},
        {"arm", arm_},
        {"matches", transition_count()},
        {"value_hops", value_hop_count()},
        {"basis_matches", basis_match_count()},
        {"tiles", tile_sources_.size()},
        {"contraction_operations", contraction_operations_.load()},
        {"contracted_observables", contracted_observables_.load()},
        {"endpoint_batch_operations", endpoint_batch_operations_.load()},
        {"endpoint_rows", endpoint_rows_.load()},
        {"epsilon", epsilon_record()},
        {"refinement", refinement_},
        {"final_local", local_reference(final_local())},
        {"strong_ownership", json::object{
             {"tile_plan", true}, {"anchor", true},
             {"basis_locals", consumed_compact_ ? 0 : basis_owner_count()},
             {"matches", consumed_compact_ ? 0 : matches_.size()},
             {"tile_sources", tile_sources_.size()}}},
        {"elapsed_ms", elapsed_ms_}};
  }

  json::object stats_json() const {
    auto result = summary();
    result["stats_queries"] = stats_queries_.fetch_add(1) + 1;
    return result;
  }

  json::object checkpoint_record() const {
    return json::object{
        {"schema", consumed_certificate_only_
             ? "diffexp2-retained-transport-arm-state-v5"
             : consumed_compact_
             ? "diffexp2-retained-transport-arm-state-v4"
             : compact_provenance_
             ? "diffexp2-retained-transport-arm-state-v3"
             : "diffexp2-retained-transport-arm-state-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"provenance", provenance_record()},
        {"elapsed_ms", elapsed_ms_},
        {"runtime_stats",
         json::object{{"stats_queries", stats_queries_.load()},
                      {"contraction_operations",
                       contraction_operations_.load()},
                      {"contracted_observables",
                       contracted_observables_.load()},
                      {"endpoint_batch_operations",
                       endpoint_batch_operations_.load()},
                      {"endpoint_rows", endpoint_rows_.load()}}}};
  }

  void restore_runtime_stats(std::uint64_t stats_queries,
                             std::uint64_t contraction_operations,
                             std::uint64_t contracted_observables,
                             std::uint64_t endpoint_batch_operations,
                             std::uint64_t endpoint_rows) {
    stats_queries_.store(stats_queries);
    contraction_operations_.store(contraction_operations);
    contracted_observables_.store(contracted_observables);
    endpoint_batch_operations_.store(endpoint_batch_operations);
    endpoint_rows_.store(endpoint_rows);
  }

 private:
  static json::object local_reference(
      const std::shared_ptr<StoredLocalBase>& local) {
    if (!local)
      throw std::logic_error(
          "transport-arm state contains a null local owner");
    return json::object{
        {"local", local->handle()}, {"chart", local->source_chart()},
        {"source_operator_identity", local->source_operator_identity()},
        {"checkpoint_identity", local->checkpoint_identity()},
        {"coefficient_domain", local->scalar_domain()}};
  }

  json::array basis_reference() const {
    if (consumed_compact_) return cached_basis_references_;
    json::array result;
    result.reserve(basis_owners_.size());
    for (const auto& basis : basis_owners_) {
      json::array columns;
      columns.reserve(basis.size());
      for (const auto& column : basis)
        columns.push_back(local_reference(column));
      result.push_back(std::move(columns));
    }
    return result;
  }

  json::array match_reference() const {
    if (consumed_compact_) return cached_match_references_;
    json::array result;
    result.reserve(matches_.size());
    for (std::size_t index = 0; index < matches_.size(); ++index) {
      auto reference = json::object{
          {"index", index}, {"match", matches_[index]->handle()},
          {"checkpoint_identity", matches_[index]->checkpoint_identity()}};
      if (!compact_provenance_)
        reference["provenance_identity"] =
            matches_[index]->provenance_identity();
      result.push_back(std::move(reference));
    }
    return result;
  }

  json::array tile_source_reference() const {
    json::array result;
    result.reserve(tile_sources_.size());
    for (std::size_t index = 0; index < tile_sources_.size(); ++index) {
      auto record = local_reference(tile_sources_[index]);
      record["tile"] = index;
      result.push_back(std::move(record));
    }
    return result;
  }

  static json::object compact_tile_checkpoint_reference(
      const std::shared_ptr<StoredLocalBase>& local, std::size_t tile) {
    if (!local)
      throw std::logic_error(
          "certificate-only transport state contains a null tile local");
    return json::object{
        {"tile", tile}, {"local", local->handle()},
        {"chart", local->source_chart()},
        {"checkpoint_identity", local->checkpoint_identity()},
        {"coefficient_domain", local->scalar_domain()}};
  }

  json::array consumed_certificate_chain() const {
    json::array chain;
    chain.reserve(tile_sources_.size());
    for (std::size_t tile = 0; tile < tile_sources_.size(); ++tile) {
      auto record = compact_tile_checkpoint_reference(
          tile_sources_[tile], tile);
      if (tile == 0) {
        record["derivation"] = nullptr;
      } else {
        const auto& derivation = *tile_sources_[tile]->retained_derivation();
        const auto schema = required_string(derivation, "schema");
        if (schema == "diffexp2-retained-plan-value-handoff-v1") {
          record["derivation"] = json::object{
              {"schema",
               "diffexp2-consumed-plan-value-handoff-certificate-v1"},
              {"match", tile - 1},
              {"handoff_provenance_identity",
               derivation.at("provenance_identity")},
              {"incoming_checkpoint_identity",
               as_object(derivation.at("incoming"),
                         "certificate-only value incoming derivation")
                   .at("checkpoint_identity")},
              {"output_checkpoint_identity",
               as_object(derivation.at("output"),
                         "certificate-only value output derivation")
                   .at("checkpoint_identity")}};
        } else {
          record["derivation"] = json::object{
              {"schema", "diffexp2-consumed-plan-match-certificate-v1"},
              {"match", tile - 1},
              {"source_match_checkpoint_identity",
               derivation.at("source_match_checkpoint_identity")},
              {"incoming_checkpoint_identity",
               as_object(derivation.at("incoming"),
                         "certificate-only incoming derivation")
                   .at("checkpoint_identity")},
              {"output_checkpoint_identity",
               as_object(derivation.at("output"),
                         "certificate-only output derivation")
                   .at("checkpoint_identity")}};
        }
      }
      chain.push_back(std::move(record));
    }
    return chain;
  }

  json::object epsilon_record() const {
    return json::object{
        {"min", work_epsilon_.min_power},
        {"max", work_epsilon_.complete_max},
        {"required_complete_max", public_required_complete_max_},
        {"match_required_complete_max", match_required_complete_max_}};
  }

  json::object provenance_record() const {
    auto plan_reference = json::object{
        {"handle", plan_owner_->handle()},
        {"checkpoint_identity", plan_owner_->checkpoint_identity()}};
    if (consumed_certificate_only_)
      return json::object{
          {"schema", "diffexp2-retained-native-transport-arm-state-v4"},
          {"checkpoint_identity", checkpoint_identity_},
          {"tile_plan", std::move(plan_reference)},
          {"arm", arm_},
          {"tile_checkpoint_chain", consumed_certificate_chain()},
          {"epsilon", epsilon_record()},
          {"refinement", refinement_}};
    if (!compact_provenance_)
      plan_reference["provenance_identity"] =
          plan_owner_->provenance_identity();
    return json::object{
        {"schema", compact_provenance_
             ? (consumed_compact_
                    ? "diffexp2-retained-native-transport-arm-state-v3"
                    : "diffexp2-retained-native-transport-arm-state-v2")
             : "diffexp2-retained-native-transport-arm-state-v1"},
        {"checkpoint_identity", checkpoint_identity_},
        {"tile_plan", std::move(plan_reference)},
        {"arm", arm_},
        {"anchor", local_reference(anchor_owner_)},
        {"receiving_basis", basis_reference()},
        {"matches", match_reference()},
        {"tile_sources", tile_source_reference()},
        {"final_local", local_reference(final_local())},
        {"epsilon", epsilon_record()},
        {"refinement", refinement_}};
  }

  std::size_t basis_owner_count() const {
    std::size_t result = 0;
    for (const auto& basis : basis_owners_) {
      if (basis.size() > std::numeric_limits<std::size_t>::max() - result)
        throw std::overflow_error(
            "transport-arm basis owner count overflow");
      result += basis.size();
    }
    return result;
  }

  std::size_t transition_count() const {
    return consumed_certificate_only_
        ? tile_sources_.size() - 1
        : consumed_compact_ ? cached_match_references_.size()
                            : matches_.size();
  }

  std::size_t value_hop_count() const {
    if (!consumed_certificate_only_) return 0;
    std::size_t count = 0;
    for (std::size_t tile = 1; tile < tile_sources_.size(); ++tile) {
      const auto& derivation = tile_sources_[tile]->retained_derivation();
      if (derivation.has_value() &&
          required_string(*derivation, "schema") ==
              "diffexp2-retained-plan-value-handoff-v1")
        ++count;
    }
    return count;
  }

  std::size_t basis_match_count() const {
    return transition_count() - value_hop_count();
  }

  void validate() const {
    if (handle_.empty() || checkpoint_identity_.empty() ||
        !plan_owner_ || !anchor_owner_ || elapsed_ms_ < 0.0 ||
        !std::isfinite(elapsed_ms_))
      throw std::invalid_argument(
          "retained transport-arm state lost an identity or strong owner");
    const auto& retained = plan_owner_->arm(arm_);
    if ((!consumed_compact_ &&
         (basis_owners_.size() != retained.exact.matches.size() ||
          matches_.size() != retained.exact.matches.size())) ||
        (consumed_compact_ && !consumed_certificate_only_ &&
         (cached_basis_references_.size() != retained.exact.matches.size() ||
          cached_match_references_.size() != retained.exact.matches.size())) ||
        tile_sources_.size() != retained.exact.tiles.size() ||
        tile_sources_.size() != retained.exact.matches.size() + 1)
      throw std::invalid_argument(
          "retained transport-arm state does not reproduce its plan topology");
    if (tile_sources_.empty() || tile_sources_.front().get() !=
                                     anchor_owner_.get())
      throw std::invalid_argument(
          "retained transport-arm state lost its anchor tile source");
    (void)work_epsilon_.width();
    if (public_required_complete_max_ < work_epsilon_.min_power ||
        match_required_complete_max_ < public_required_complete_max_ ||
        match_required_complete_max_ > work_epsilon_.complete_max)
      throw std::invalid_argument(
          "retained transport-arm epsilon contract is inconsistent");
    require_exact_keys(refinement_, {"relative_tolerance", "max_steps"},
                       "retained transport-arm refinement policy");
    if (required_string(refinement_, "relative_tolerance").empty() ||
        as_u32(refinement_.at("max_steps"),
               "retained transport-arm refinement steps") > 32)
      throw std::invalid_argument(
          "retained transport-arm refinement policy is invalid");

    for (std::size_t tile = 0; tile < tile_sources_.size(); ++tile) {
      const auto& source = tile_sources_[tile];
      if (!source)
        throw std::invalid_argument(
            "retained transport-arm state contains a null tile source");
      const auto& exact_tile = retained.exact.tiles[tile];
      const auto& chart = retained.charts.at(exact_tile.chart);
      if (source->source_chart() != chart.handle)
        throw std::invalid_argument(
            "retained transport-arm tile source belongs to a different chart");
      source->require_exact_plan_binding(
          chart.geometry, chart.prescriptions,
          "retained transport-arm tile source");
    }
    if (consumed_certificate_only_) {
      validate_consumed_certificate_chain(retained);
      return;
    }
    if (consumed_compact_) {
      validate_consumed_references(retained);
      return;
    }
    for (std::size_t index = 0; index < matches_.size(); ++index) {
      const auto& match = matches_[index];
      const auto& basis = basis_owners_[index];
      if (!match || basis.empty() ||
          std::any_of(basis.begin(), basis.end(),
                      [](const auto& owner) { return owner == nullptr; }) ||
          match->plan_owner().get() != plan_owner_.get() ||
          match->incoming_owner().get() != tile_sources_[index].get() ||
          match->basis_owners().size() != basis.size())
        throw std::invalid_argument(
            "retained transport-arm match lost its exact owner set");
      for (std::size_t column = 0; column < basis.size(); ++column)
        if (match->basis_owners()[column].get() != basis[column].get())
          throw std::invalid_argument(
              "retained transport-arm basis differs from its match owner");
      const auto& handoff = match->handoff();
      if (required_string(handoff, "arm") != arm_ ||
          as_u64(handoff.at("match"),
                 "retained transport-arm match index") != index)
        throw std::invalid_argument(
            "retained transport-arm match provenance is out of order");
      const auto& next = tile_sources_[index + 1];
      if (!next->retained_derivation().has_value() ||
          next->retained_derivation_owner().get() != match.get() ||
          required_string(*next->retained_derivation(), "source_match") !=
              match->handle())
        throw std::invalid_argument(
            "retained transport-arm tile source is not materialized from its match");
    }
  }

  void validate_consumed_certificate_chain(
      const RetainedArmPlan& retained) const {
    if (tile_sources_.size() != retained.exact.matches.size() + 1)
      throw std::invalid_argument(
          "certificate-only transport state has the wrong tile chain length");
    for (std::size_t tile = 1; tile < tile_sources_.size(); ++tile) {
      const auto& incoming = tile_sources_[tile - 1];
      const auto& output = tile_sources_[tile];
      if (!output || !output->has_sealed_plan_match_lineage() ||
          output->retained_derivation_owner() != nullptr ||
          !output->retained_derivation().has_value())
        throw std::invalid_argument(
            "certificate-only transport tile has no sealed plan-match derivation");
      const auto& derivation = *output->retained_derivation();
      const auto schema = required_string(derivation, "schema");
      const bool value_handoff = schema ==
          "diffexp2-retained-plan-value-handoff-v1";
      if ((!value_handoff && schema !=
              "diffexp2-retained-plan-match-local-materialization-v2") ||
          required_string(derivation, "tile_plan") !=
              plan_owner_->handle() ||
          required_string(derivation, "tile_plan_checkpoint_identity") !=
              plan_owner_->checkpoint_identity() ||
          required_string(derivation, "arm") != arm_ ||
          as_u64(derivation.at("match"),
                 "certificate-only plan-match index") != tile - 1)
        throw std::invalid_argument(
            "certificate-only transport derivation differs from its retained plan");
      const auto& source = as_object(
          derivation.at("incoming"),
          "certificate-only transport incoming derivation");
      const auto& target = as_object(
          derivation.at("output"),
          "certificate-only transport output derivation");
      if (required_string(source, "local") != incoming->handle() ||
          required_string(source, "checkpoint_identity") !=
              incoming->checkpoint_identity() ||
          required_string(source, "source_operator_identity") !=
              incoming->source_operator_identity() ||
          required_string(target, "checkpoint_identity") !=
              output->checkpoint_identity() ||
          required_string(target, "chart") != output->source_chart() ||
          required_string(target, "source_operator_identity") !=
              output->source_operator_identity())
        throw std::invalid_argument(
            "certificate-only transport derivation breaks its ordered tile chain");
      const auto equation_owner = output->retained_equation_owner();
      if (!equation_owner ||
          required_string(derivation,
                          "equation_owner_signature_identity") !=
              equation_owner->owner_signature_identity() ||
          required_string(derivation, "equation_payload_identity") !=
              equation_owner->physical_payload_identity())
        throw std::invalid_argument(
            "certificate-only transport derivation lost its match/equation certificate");
      if (!value_handoff &&
          required_string(derivation,
                          "source_match_checkpoint_identity").empty())
        throw std::invalid_argument(
            "certificate-only plan match lost its source-match checkpoint");
    }
  }

  void validate_consumed_references(const RetainedArmPlan& retained) const {
    for (std::size_t index = 0; index < cached_basis_references_.size();
         ++index) {
      const auto& raw_basis = as_array(
          cached_basis_references_[index],
          "consumed transport-arm cached basis");
      if (raw_basis.empty())
        throw std::invalid_argument(
            "consumed transport-arm cached basis cannot be empty");
      for (const auto& raw_column : raw_basis) {
        const auto& column = as_object(
            raw_column, "consumed transport-arm cached basis column");
        require_exact_keys(
            column,
            {"local", "chart", "source_operator_identity",
             "checkpoint_identity", "coefficient_domain"},
            "consumed transport-arm cached basis column");
        (void)scoped_handle_id(required_string(column, "local"), "l:",
                               "consumed basis local");
        if (required_string(column, "chart") !=
                retained.charts.at(
                    retained.exact.matches[index].receiving_chart)
                    .handle ||
            required_string(column, "source_operator_identity").empty() ||
            required_string(column, "checkpoint_identity").empty() ||
            required_string(column, "coefficient_domain") !=
                tile_sources_[index + 1]->scalar_domain())
          throw std::invalid_argument(
              "consumed transport-arm cached basis identity is inconsistent");
      }

      const auto& match = as_object(
          cached_match_references_[index],
          "consumed transport-arm cached match");
      require_exact_keys(
          match,
          {"index", "checkpoint_identity", "provenance_identity",
           "planned_hop", "sealed_local_lineage"},
          "consumed transport-arm cached match");
      if (as_u64(match.at("index"),
                 "consumed transport-arm match index") != index)
        throw std::invalid_argument(
            "consumed transport-arm cached matches are out of order");
      const auto& handoff = as_object(
          match.at("planned_hop"),
          "consumed transport-arm planned hop");
      const auto handoff_schema = required_string(handoff, "schema");
      const bool compact_handoff = handoff_schema ==
          "diffexp2-retained-exact-plan-match-hop-v2";
      if (compact_handoff)
        require_exact_keys(
            handoff,
            {"schema", "tile_plan", "tile_plan_checkpoint_identity",
             "arm", "match", "geometry", "producing", "receiving",
             "result_checkpoint_identity", "advance"},
            "consumed compact transport-arm planned hop");
      else
        require_exact_keys(
            handoff,
            {"schema", "tile_plan", "tile_plan_checkpoint_identity",
             "tile_plan_provenance_identity", "arm", "match", "geometry",
             "producing", "receiving", "result_checkpoint_identity",
             "native_match_provenance_identity", "advance"},
            "consumed transport-arm planned hop");
      if ((!compact_handoff && handoff_schema !=
              "diffexp2-retained-exact-plan-match-hop-v1") ||
          required_string(handoff, "tile_plan") != plan_owner_->handle() ||
          required_string(handoff, "tile_plan_checkpoint_identity") !=
              plan_owner_->checkpoint_identity() ||
          (!compact_handoff &&
           required_string(handoff, "tile_plan_provenance_identity") !=
              plan_owner_->provenance_identity()) ||
          required_string(handoff, "arm") != arm_ ||
          as_u64(handoff.at("match"),
                 "consumed planned-hop match index") != index ||
          handoff.at("geometry") != encode_plan_match(retained, index) ||
          required_string(handoff, "result_checkpoint_identity") !=
              required_string(match, "checkpoint_identity") ||
          (compact_handoff &&
           json::serialize(canonical_json_value(handoff)) !=
              required_string(match, "provenance_identity")) ||
          (!compact_handoff &&
           required_string(handoff, "native_match_provenance_identity") !=
              required_string(match, "provenance_identity")))
        throw std::invalid_argument(
            "consumed transport-arm planned hop differs from its plan");
      const auto& producing = as_object(
          handoff.at("producing"),
          "consumed transport-arm producing handoff");
      const auto& incoming = as_object(
          producing.at("incoming"),
          "consumed transport-arm incoming handoff");
      if (required_string(incoming, "local") !=
              tile_sources_[index]->handle() ||
          required_string(incoming, "checkpoint_identity") !=
              tile_sources_[index]->checkpoint_identity() ||
          required_string(incoming, "source_operator_identity") !=
              tile_sources_[index]->source_operator_identity())
        throw std::invalid_argument(
            "consumed transport-arm incoming lineage is inconsistent");
      const auto& handoff_basis = as_array(
          as_object(handoff.at("receiving"),
                    "consumed transport-arm receiving handoff")
              .at("basis"),
          "consumed transport-arm receiving handoff basis");
      if (handoff_basis.size() != raw_basis.size())
        throw std::invalid_argument(
            "consumed transport-arm cached basis dimension changed");
      for (std::size_t column = 0; column < raw_basis.size(); ++column) {
        const auto& cached = as_object(
            raw_basis[column], "consumed transport-arm cached basis column");
        const auto& handed = as_object(
            handoff_basis[column],
            "consumed transport-arm handed-off basis column");
        if (as_u64(handed.at("column"),
                   "consumed transport-arm basis column") != column ||
            handed.at("local") != cached.at("local") ||
            handed.at("checkpoint_identity") !=
                cached.at("checkpoint_identity") ||
            handed.at("source_operator_identity") !=
                cached.at("source_operator_identity"))
          throw std::invalid_argument(
              "consumed transport-arm cached basis differs from its planned hop");
      }

      const auto& local = tile_sources_[index + 1];
      if (!local || !local->has_sealed_plan_match_lineage() ||
          local->retained_derivation_owner() != nullptr ||
          !local->retained_derivation().has_value())
        throw std::invalid_argument(
            "consumed transport-arm materialized local lineage is not sealed");
      const auto& derivation = *local->retained_derivation();
      const auto& sealed = as_object(
          match.at("sealed_local_lineage"),
          "consumed transport-arm sealed local lineage");
      const auto equation_owner = local->retained_equation_owner();
      const auto sealed_schema = required_string(sealed, "schema");
      if (sealed_schema ==
          "diffexp2-sealed-plan-match-local-lineage-v2") {
        require_exact_keys(
            sealed,
            {"schema", "local", "local_checkpoint_identity",
             "source_operator_identity", "match",
             "match_checkpoint_identity", "tile_plan",
             "tile_plan_checkpoint_identity", "arm", "match_index",
             "incoming_checkpoint_identity",
             "equation_owner_signature_identity",
             "equation_payload_identity"},
            "consumed compact transport-arm sealed local lineage");
        if (required_string(derivation, "schema") !=
                "diffexp2-retained-plan-match-local-materialization-v2" ||
            required_string(sealed, "local") != local->handle() ||
            required_string(sealed, "local_checkpoint_identity") !=
                local->checkpoint_identity() ||
            required_string(sealed, "source_operator_identity") !=
                local->source_operator_identity() ||
            sealed.at("match") != derivation.at("source_match") ||
            sealed.at("match_checkpoint_identity") !=
                match.at("checkpoint_identity") ||
            sealed.at("match_checkpoint_identity") !=
                derivation.at("source_match_checkpoint_identity") ||
            sealed.at("tile_plan") != derivation.at("tile_plan") ||
            required_string(sealed, "tile_plan") !=
                plan_owner_->handle() ||
            sealed.at("tile_plan_checkpoint_identity") !=
                derivation.at("tile_plan_checkpoint_identity") ||
            required_string(sealed, "tile_plan_checkpoint_identity") !=
                plan_owner_->checkpoint_identity() ||
            sealed.at("arm") != derivation.at("arm") ||
            required_string(sealed, "arm") != arm_ ||
            sealed.at("match_index") != derivation.at("match") ||
            as_u64(sealed.at("match_index"),
                   "consumed compact match index") != index ||
            sealed.at("incoming_checkpoint_identity") !=
                as_object(derivation.at("incoming"),
                          "consumed compact derivation incoming")
                    .at("checkpoint_identity") ||
            required_string(sealed, "incoming_checkpoint_identity") !=
                tile_sources_[index]->checkpoint_identity() ||
            derivation.at("basis") != handoff_basis ||
            !equation_owner ||
            required_string(sealed,
                "equation_owner_signature_identity") !=
                equation_owner->owner_signature_identity() ||
            required_string(sealed, "equation_payload_identity") !=
                equation_owner->physical_payload_identity())
          throw std::invalid_argument(
              "consumed compact transport-arm sealed lineage is inconsistent");
      } else {
        require_exact_keys(
            sealed,
            {"schema", "local", "local_checkpoint_identity",
             "source_operator_identity", "match",
             "match_checkpoint_identity", "match_provenance_identity",
             "planned_hop_provenance_identity",
             "derivation_provenance_identity",
             "equation_owner_signature_identity",
             "equation_payload_identity"},
            "consumed transport-arm sealed local lineage");
        if (sealed_schema !=
                "diffexp2-sealed-plan-match-local-lineage-v1" ||
          required_string(sealed, "local") != local->handle() ||
          required_string(sealed, "local_checkpoint_identity") !=
              local->checkpoint_identity() ||
          required_string(sealed, "source_operator_identity") !=
              local->source_operator_identity() ||
          sealed.at("match") != derivation.at("source_match") ||
          sealed.at("match_checkpoint_identity") !=
              match.at("checkpoint_identity") ||
          sealed.at("match_checkpoint_identity") !=
              derivation.at("source_match_checkpoint_identity") ||
          sealed.at("match_provenance_identity") !=
              match.at("provenance_identity") ||
          sealed.at("match_provenance_identity") !=
              derivation.at("source_match_provenance_identity") ||
          sealed.at("planned_hop_provenance_identity") !=
              derivation.at("planned_hop_provenance_identity") ||
          required_string(
              derivation, "planned_hop_provenance_identity") !=
              json::serialize(canonical_json_value(handoff)) ||
          sealed.at("derivation_provenance_identity") !=
              derivation.at("provenance_identity") ||
          derivation.at("planned_hop") != match.at("planned_hop") ||
          !equation_owner ||
          required_string(sealed, "equation_owner_signature_identity") !=
              equation_owner->owner_signature_identity() ||
          required_string(sealed, "equation_payload_identity") !=
              equation_owner->physical_payload_identity())
          throw std::invalid_argument(
              "consumed transport-arm sealed lineage identity is inconsistent");
      }
    }
  }

  std::string handle_;
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::string arm_;
  std::shared_ptr<StoredTilePlan> plan_owner_;
  std::shared_ptr<StoredLocalBase> anchor_owner_;
  std::vector<std::vector<std::shared_ptr<StoredLocalBase>>> basis_owners_;
  std::vector<std::shared_ptr<StoredPlannedMatchHop>> matches_;
  std::vector<std::shared_ptr<StoredLocalBase>> tile_sources_;
  EpsilonWindow work_epsilon_;
  std::int32_t public_required_complete_max_ = 0;
  std::int32_t match_required_complete_max_ = 0;
  json::object refinement_;
  double elapsed_ms_ = 0.0;
  bool compact_provenance_ = true;
  bool consumed_compact_ = false;
  bool consumed_certificate_only_ = false;
  json::array cached_basis_references_;
  json::array cached_match_references_;
  mutable std::atomic<std::uint64_t> stats_queries_{0};
  std::atomic<std::uint64_t> contraction_operations_{0};
  std::atomic<std::uint64_t> contracted_observables_{0};
  std::atomic<std::uint64_t> endpoint_batch_operations_{0};
  std::atomic<std::uint64_t> endpoint_rows_{0};
};

void require_transport_pair_compatibility(
    const std::shared_ptr<StoredTransportArmState>& lower,
    const std::shared_ptr<StoredTransportArmState>& upper,
    const std::string& expected_domain) {
  if (!lower || !upper || lower.get() == upper.get())
    throw std::invalid_argument(
        "retained transport pair requires two distinct state owners");
  if (lower->arm_name() != "lower" || upper->arm_name() != "upper")
    throw std::invalid_argument(
        "retained transport pair requires lower and upper state directions");
  if (expected_domain != "rational" && expected_domain != "acb")
    throw std::invalid_argument(
        "retained transport pair requires one numeric session domain");
  const auto& lower_anchor = lower->anchor_owner();
  const auto& upper_anchor = upper->anchor_owner();
  if (!lower_anchor || !upper_anchor ||
      lower_anchor.get() != upper_anchor.get() ||
      lower_anchor->checkpoint_identity() !=
          upper_anchor->checkpoint_identity() ||
      lower_anchor->source_operator_identity() !=
          upper_anchor->source_operator_identity() ||
      std::string(lower_anchor->scalar_domain()) != expected_domain ||
      std::string(upper_anchor->scalar_domain()) != expected_domain ||
      lower_anchor->exact_analytic_metadata() !=
          upper_anchor->exact_analytic_metadata())
    throw std::invalid_argument(
        "retained transport pair does not share one exact anchor owner, checkpoint, operator, domain, and analytic metadata");
  const auto& lower_plan = lower->plan_owner()->arm("lower").exact;
  const auto& upper_plan = upper->plan_owner()->arm("upper").exact;
  if (lower_plan.direction != -1 || upper_plan.direction != 1 ||
      lower_plan.from.str() != upper_plan.from.str())
    throw std::invalid_argument(
        "retained transport pair plans do not share one compatible exact physical anchor");
}

struct ResolvedTransportEndpointBinding {
  json::object source;
  std::string arm;
  Rational local_end;
  std::int32_t approach_direction = 0;
  std::optional<std::int32_t> rim;
  bool centered = false;
};

ResolvedTransportEndpointBinding resolve_transport_endpoint_binding(
    const std::shared_ptr<StoredTransportArmState>& state) {
  if (!state)
    throw std::invalid_argument(
        "transport endpoint binding requires a retained arm state");
  const auto& plan = state->plan_owner();
  const auto& arm_name = state->arm_name();
  const auto& arm = plan->arm(arm_name);
  if (arm.exact.tiles.empty() || arm.charts.empty())
    throw std::invalid_argument(
        "transport endpoint arm has no final tile/chart");
  const auto final_tile_index = arm.exact.tiles.size() - 1;
  const auto& final_tile = arm.exact.tiles.back();
  if (final_tile.chart >= arm.charts.size())
    throw std::logic_error(
        "transport endpoint final tile has an invalid chart index");
  const auto& final_chart = arm.charts[final_tile.chart];
  if (!(final_tile.physical_end == arm.exact.to) ||
      final_chart.geometry.scale.is_zero() ||
      arm.exact.direction != (arm_name == "lower" ? -1 : 1))
    throw std::invalid_argument(
        "transport endpoint plan has an inconsistent final physical binding");
  const auto mapped_endpoint = final_chart.geometry.center +
      final_chart.geometry.scale * final_tile.local_end;
  if (!(mapped_endpoint == arm.exact.to))
    throw std::invalid_argument(
        "transport endpoint final local coordinate does not map to its exact physical endpoint");
  const auto& local = state->final_local();
  if (!local || local->source_chart() != final_chart.handle)
    throw std::invalid_argument(
        "transport endpoint final local does not name its final retained chart");
  local->require_exact_plan_binding(
      final_chart.geometry, final_chart.prescriptions,
      "transport endpoint final local");

  ResolvedTransportEndpointBinding binding;
  binding.arm = arm_name;
  binding.local_end = final_tile.local_end;
  binding.centered = final_tile.local_end.is_zero();
  binding.approach_direction =
      -arm.exact.direction * final_chart.geometry.scale.sign();
  binding.rim = exact_plan_rim(
      final_chart.prescriptions, final_chart.geometry.scale);
  binding.source = json::object{
      {"transport_state", json::object{
           {"handle", state->handle()},
           {"checkpoint_identity", state->checkpoint_identity()},
           {"provenance_identity", state->provenance_identity()}}},
      {"tile_plan", json::object{
           {"handle", plan->handle()},
           {"checkpoint_identity", plan->checkpoint_identity()},
           {"provenance_identity", plan->provenance_identity()}}},
      {"arm", arm_name},
      {"endpoint_exact", arm.exact.to.str()},
      {"direction", arm.exact.direction},
      {"final_tile", final_tile_index},
      {"final_chart_index", final_tile.chart},
      {"final_chart", final_chart.handle},
      {"final_chart_identity", final_chart.exact_identity},
      {"local_endpoint_exact", final_tile.local_end.str()},
      {"centered", binding.centered},
      {"approach_direction", binding.approach_direction},
      {"derived_rim", binding.rim.has_value()
           ? json::value(*binding.rim) : json::value(nullptr)},
      {"final_local", json::object{
           {"handle", local->handle()}, {"chart", local->source_chart()},
           {"source_operator_identity", local->source_operator_identity()},
           {"checkpoint_identity", local->checkpoint_identity()},
           {"coefficient_domain", local->scalar_domain()}}},
      {"prescriptions", encode_plan_prescriptions(
           final_chart.prescriptions)}};
  return binding;
}

void validate_prepared_rational_row_structure(
    const json::value& raw, std::uint32_t expected_columns,
    const char* label) {
  const auto& row = as_object(raw, label);
  require_exact_keys(row,
                     {"schema", "columns", "exact_identity", "entries"},
                     label);
  if (required_string(row, "schema") !=
          "diffexp2-prepared-rational-local-row-v1" ||
      as_u32(row.at("columns"), label) != expected_columns ||
      required_string(row, "exact_identity").empty())
    throw std::invalid_argument(
        std::string(label) + " has an incompatible schema, dimension, or identity");
  std::optional<std::uint32_t> previous;
  for (const auto& raw_entry : as_array(row.at("entries"), label)) {
    const auto& entry = as_object(raw_entry, label);
    require_exact_keys(entry, {"column", "multiplier"}, label);
    const auto column = as_u32(entry.at("column"), label);
    if (column >= expected_columns ||
        (previous.has_value() && *previous >= column))
      throw std::invalid_argument(
          std::string(label) + " columns are not strictly ordered in range");
    previous = column;
    const auto& multiplier = as_object(entry.at("multiplier"), label);
    const bool has_kernels =
        multiplier.if_contains("kernels") != nullptr;
    const bool has_analytic_coefficients =
        multiplier.if_contains("analytic_coefficients") != nullptr;
    if (has_kernels && has_analytic_coefficients)
      require_exact_keys(
          multiplier,
          {"epsilon_shift", "center_pole_order", "kernels",
           "exact_identity", "proven_zero", "analytic_coefficients"},
          label);
    else if (has_kernels)
      require_exact_keys(
          multiplier,
          {"epsilon_shift", "center_pole_order", "kernels",
           "exact_identity", "proven_zero"}, label);
    else if (has_analytic_coefficients)
      require_exact_keys(
          multiplier,
          {"epsilon_shift", "center_pole_order", "exact_identity",
           "proven_zero", "analytic_coefficients"}, label);
    else
      throw std::invalid_argument(
          std::string(label) +
          " multiplier needs Taylor kernels or an analytic rational source");
    (void)as_i32(multiplier.at("epsilon_shift"), label);
    (void)as_u32(multiplier.at("center_pole_order"), label);
    if (required_string(multiplier, "exact_identity").empty() ||
        !multiplier.at("proven_zero").is_bool() ||
        multiplier.at("proven_zero").as_bool())
      throw std::invalid_argument(
          std::string(label) + " contains an invalid active multiplier");
    std::size_t kernel_count = 0;
    if (has_kernels) {
      const auto& kernels = as_array(multiplier.at("kernels"), label);
      if (kernels.empty())
        throw std::invalid_argument(
            std::string(label) + " contains no epsilon kernels");
      kernel_count = kernels.size();
      for (const auto& raw_kernel : kernels)
        if (as_array(raw_kernel, label).empty())
          throw std::invalid_argument(
              std::string(label) + " contains an empty Taylor kernel");
    }
    if (has_analytic_coefficients) {
      const auto& coefficients = as_array(
          multiplier.at("analytic_coefficients"), label);
      if (coefficients.empty() ||
          (has_kernels && coefficients.size() != kernel_count))
        throw std::invalid_argument(
            std::string(label) +
            " analytic coefficient count differs from its epsilon kernels");
      for (const auto& raw_coefficient : coefficients) {
        const auto& coefficient = as_object(raw_coefficient, label);
        require_exact_keys(coefficient, {"numerator", "denominator"}, label);
        if (as_array(coefficient.at("numerator"), label).empty() ||
            as_array(coefficient.at("denominator"), label).empty())
          throw std::invalid_argument(
              std::string(label) +
              " contains an empty analytic numerator or denominator");
      }
    }
  }
}

class StoredLineResult {
 public:
  StoredLineResult(std::string handle, std::string checkpoint_identity,
                   std::string provenance_identity, std::string arm,
                   std::size_t tile_index, json::object interval,
                   std::string source_checkpoint,
                   StoredLineIntegral result, double elapsed_ms,
                   std::shared_ptr<StoredTilePlan> plan_owner,
                   std::shared_ptr<StoredLocalBase> local_owner)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        arm_(std::move(arm)), tile_index_(tile_index),
        interval_(std::move(interval)),
        source_checkpoint_(std::move(source_checkpoint)),
        result_(std::move(result)), elapsed_ms_(elapsed_ms),
        plan_owners_{std::move(plan_owner)},
        local_owners_{std::move(local_owner)} {}

  StoredLineResult(std::string handle, std::string checkpoint_identity,
                   std::string provenance_identity,
                   StoredLineIntegral result, double elapsed_ms,
                   std::shared_ptr<StoredTilePlan> plan_owner,
                   std::vector<std::shared_ptr<StoredLocalBase>> local_owners,
                   json::object aggregate_provenance,
                   std::shared_ptr<StoredTransportArmState> transport_owner =
                       nullptr)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        result_(std::move(result)), elapsed_ms_(elapsed_ms),
        plan_owners_{std::move(plan_owner)},
        local_owners_(std::move(local_owners)),
        aggregate_provenance_(std::move(aggregate_provenance)) {
    if (transport_owner) transport_owners_.push_back(
        std::move(transport_owner));
    validate_aggregate_ownership();
  }

  StoredLineResult(
      std::string handle, std::string checkpoint_identity,
      std::string provenance_identity, StoredLineIntegral result,
      double elapsed_ms,
      std::vector<std::shared_ptr<StoredTilePlan>> plan_owners,
      std::vector<std::shared_ptr<StoredLocalBase>> local_owners,
      json::object aggregate_provenance,
      std::vector<std::shared_ptr<StoredTransportArmState>>
          transport_owners)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        result_(std::move(result)), elapsed_ms_(elapsed_ms),
        plan_owners_(std::move(plan_owners)),
        local_owners_(std::move(local_owners)),
        aggregate_provenance_(std::move(aggregate_provenance)),
        transport_owners_(std::move(transport_owners)) {
    validate_aggregate_ownership();
  }

  const std::string& handle() const { return handle_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& provenance_identity() const {
    return provenance_identity_;
  }
  double elapsed_ms() const { return elapsed_ms_; }
  const StoredLineIntegral& result() const { return result_; }
  bool is_aggregate() const { return aggregate_provenance_.has_value(); }

  const std::shared_ptr<StoredTilePlan>& plan_owner() const {
    if (plan_owners_.size() != 1)
      throw std::logic_error(
          "paired compact line has no singular retained plan owner");
    return plan_owners_.front();
  }
  const std::vector<std::shared_ptr<StoredTilePlan>>& plan_owners() const {
    return plan_owners_;
  }
  const std::shared_ptr<StoredLocalBase>& local_owner() const {
    if (local_owners_.empty())
      throw std::logic_error(
          "compact transport observable line has no retained local owner");
    return local_owners_.front();
  }
  const std::vector<std::shared_ptr<StoredLocalBase>>& local_owners() const {
    return local_owners_;
  }
  const std::shared_ptr<StoredTransportArmState>& transport_owner() const {
    if (transport_owners_.size() > 1)
      throw std::logic_error(
          "paired compact line has no singular transport-state owner");
    static const std::shared_ptr<StoredTransportArmState> empty;
    return transport_owners_.empty() ? empty : transport_owners_.front();
  }
  const std::vector<std::shared_ptr<StoredTransportArmState>>&
  transport_owners() const {
    return transport_owners_;
  }

  json::object summary() const {
    const auto& diagnostics = result_.diagnostics;
    json::object output{
        {"line", handle_},
        {"capability", transport_owners_.size() == 2
             ? kRetainedTransportPairContractionCapability
             : transport_owners_.size() == 1
             ? kRetainedTransportArmContractionCapability
             : is_aggregate()
             ? kRetainedLineAggregateCapability
             : result_.scope ==
                           LineIntegrationScope::FullLocalWithCertifiedTail
             ? kRetainedCertifiedLineCapability
             : kRetainedStoredLineCapability},
        {"native_retained", true}, {"json_coefficients", 0},
        {"scope", line_integration_scope_name(result_.scope)},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"dimension", result_.value.dimension},
        {"epsilon_min", result_.value.epsilon.min_power},
        {"epsilon_max", result_.value.epsilon.complete_max},
        {"effective_rim", result_.imaginary_sign.has_value()
             ? json::value(*result_.imaginary_sign) : json::value(nullptr)},
        {"error", encode_error_envelope_summary(result_.value.error)},
        {"diagnostics", json::object{
             {"input_monomial_cells", diagnostics.input_monomial_cells},
             {"grouped_monomials", diagnostics.grouped_monomials},
             {"zero_groups_skipped", diagnostics.zero_groups_skipped},
             {"cancelled_divergent_groups",
              diagnostics.cancelled_divergent_groups},
             {"bounded_cancelled_divergent_coefficients",
              diagnostics.bounded_cancelled_divergent_coefficients},
             {"divergent_cancellation_mode",
              diagnostics.divergent_cancellation_mode},
             {"divergent_relative_tolerance",
              diagnostics.divergent_relative_tolerance.empty()
                  ? json::value(nullptr)
                  : json::value(
                        diagnostics.divergent_relative_tolerance)},
             {"divergent_cancellation_provenance",
              diagnostics.divergent_cancellation_provenance.empty()
                  ? json::value(nullptr)
                  : json::value(
                        diagnostics.divergent_cancellation_provenance)},
             {"primitive_evaluations", diagnostics.primitive_evaluations},
             {"primitive_component_applications",
              diagnostics.primitive_component_applications},
             {"primitive_component_reuses",
              diagnostics.primitive_component_reuses},
             {"has_center_endpoint", diagnostics.has_center_endpoint},
             {"tail_certificate_requested",
              diagnostics.tail_certificate_requested},
             {"tail_certificate_status",
              diagnostics.tail_certificate_status},
             {"tail_witness_radius_exact",
              diagnostics.tail_witness_radius_exact.empty()
                  ? json::value(nullptr)
                  : json::value(
                        diagnostics.tail_witness_radius_exact)},
             {"detail", diagnostics.detail}}},
        {"elapsed_ms", elapsed_ms_}};
    if (aggregate_provenance_.has_value()) {
      output["source"] = aggregate_provenance_->at("source");
      output["arm"] = aggregate_provenance_->at("arm");
      output["tile"] = nullptr;
      output["interval"] = aggregate_provenance_->at("interval");
      output["aggregate"] = aggregate_provenance_->at("aggregate");
    } else {
      const auto& owner = local_owners_.front();
      output["source"] = json::object{
          {"tile_plan", plan_owners_.front()->handle()},
          {"tile_plan_checkpoint_identity",
           plan_owners_.front()->checkpoint_identity()},
          {"local", owner->handle()},
          {"chart", owner->source_chart()},
          {"local_checkpoint_identity", source_checkpoint_}};
      output["arm"] = arm_;
      output["tile"] = tile_index_;
      output["interval"] = interval_;
    }
    return output;
  }

  json::object stats_json() const {
    auto result = summary();
    std::lock_guard<std::mutex> lock(stats_mutex_);
    result["exports"] = exports_;
    result["export_ms"] = export_ms_;
    return result;
  }

  json::object export_values(const std::string& expected_checkpoint,
                             int output_digits) {
    if (expected_checkpoint != checkpoint_identity_)
      throw std::invalid_argument(
          "line export checkpoint identity does not match retained result");
    const auto started = std::chrono::steady_clock::now();
    auto value = encode_epsilon_vector(result_.value, output_digits);
    const auto elapsed = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    {
      std::lock_guard<std::mutex> lock(stats_mutex_);
      ++exports_;
      export_ms_ += elapsed;
    }
    return json::object{
        {"line", handle_}, {"checkpoint_identity", checkpoint_identity_},
        {"compatibility_export", true},
        {"scope", line_integration_scope_name(result_.scope)},
        {"error_guarantee",
         error_guarantee_name(result_.value.error.guarantee)},
        {"json_coefficients", value.at("coefficients").as_array().size()},
        {"value", std::move(value)}, {"elapsed_ms", elapsed}};
  }

  json::object checkpoint_record() const {
    const auto& diagnostics = result_.diagnostics;
    std::lock_guard<std::mutex> lock(stats_mutex_);
    if (aggregate_provenance_.has_value()) {
      const auto aggregate_schema = required_string(
          *aggregate_provenance_, "schema");
      const bool compact_v2 =
          aggregate_schema ==
              "diffexp2-retained-native-transport-observable-line-v2" ||
          aggregate_schema ==
              "diffexp2-retained-native-transport-pair-observable-line-v2";
      return json::object{
          {"schema", transport_owners_.size() == 2
               ? compact_v2
                   ? "diffexp2-retained-transport-pair-observable-line-v2"
                   : "diffexp2-retained-transport-pair-observable-line-v1"
               : transport_owners_.size() == 1
               ? compact_v2
                   ? "diffexp2-retained-transport-observable-line-v2"
                   : "diffexp2-retained-transport-observable-line-v1"
               : "diffexp2-retained-line-aggregate-v1"},
          {"handle", handle_},
          {"checkpoint_identity", checkpoint_identity_},
          {"provenance_identity", provenance_identity_},
          {"provenance", *aggregate_provenance_},
          {"result",
           json::object{
               {"value", checkpoint_epsilon_vector_record(result_.value)},
               {"scope", line_integration_scope_name(result_.scope)},
               {"imaginary_sign", result_.imaginary_sign.has_value()
                    ? json::value(*result_.imaginary_sign)
                    : json::value(nullptr)},
               {"diagnostics",
                json::object{
                    {"input_monomial_cells",
                     diagnostics.input_monomial_cells},
                    {"grouped_monomials", diagnostics.grouped_monomials},
                    {"zero_groups_skipped",
                     diagnostics.zero_groups_skipped},
                    {"cancelled_divergent_groups",
                     diagnostics.cancelled_divergent_groups},
                    {"bounded_cancelled_divergent_coefficients",
                     diagnostics.bounded_cancelled_divergent_coefficients},
                    {"divergent_cancellation_mode",
                     diagnostics.divergent_cancellation_mode},
                    {"divergent_relative_tolerance",
                     diagnostics.divergent_relative_tolerance},
                    {"divergent_cancellation_provenance",
                     diagnostics.divergent_cancellation_provenance},
                    {"primitive_evaluations",
                     diagnostics.primitive_evaluations},
                    {"primitive_component_applications",
                     diagnostics.primitive_component_applications},
                    {"primitive_component_reuses",
                     diagnostics.primitive_component_reuses},
                    {"has_center_endpoint",
                     diagnostics.has_center_endpoint},
                    {"tail_certificate_requested",
                     diagnostics.tail_certificate_requested},
                    {"tail_certificate_status",
                     diagnostics.tail_certificate_status},
                    {"tail_witness_radius_exact",
                     diagnostics.tail_witness_radius_exact},
                    {"detail", diagnostics.detail}}}}},
          {"elapsed_ms", elapsed_ms_},
          {"runtime_stats",
           json::object{{"exports", exports_}, {"export_ms", export_ms_}}}};
    }
    const auto& owner = local_owners_.front();
    return json::object{
        {"schema", "diffexp2-retained-line-result-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"arm", arm_},
        {"tile", tile_index_},
        {"interval", interval_},
        {"source",
         json::object{
             {"tile_plan", plan_owners_.front()->handle()},
             {"tile_plan_checkpoint_identity",
              plan_owners_.front()->checkpoint_identity()},
             {"local", owner->handle()},
             {"chart", owner->source_chart()},
             {"source_operator_identity",
              owner->source_operator_identity()},
             {"local_checkpoint_identity", source_checkpoint_},
             {"coefficient_domain", owner->scalar_domain()},
             {"analytic_metadata",
              owner->exact_analytic_metadata()}}},
        {"result",
         json::object{
             {"value", checkpoint_epsilon_vector_record(result_.value)},
             {"scope", line_integration_scope_name(result_.scope)},
             {"imaginary_sign", result_.imaginary_sign.has_value()
                  ? json::value(*result_.imaginary_sign)
                  : json::value(nullptr)},
             {"diagnostics",
              json::object{
                  {"input_monomial_cells",
                   diagnostics.input_monomial_cells},
                  {"grouped_monomials", diagnostics.grouped_monomials},
                  {"zero_groups_skipped", diagnostics.zero_groups_skipped},
                  {"cancelled_divergent_groups",
                   diagnostics.cancelled_divergent_groups},
                  {"bounded_cancelled_divergent_coefficients",
                   diagnostics.bounded_cancelled_divergent_coefficients},
                  {"divergent_cancellation_mode",
                   diagnostics.divergent_cancellation_mode},
                  {"divergent_relative_tolerance",
                   diagnostics.divergent_relative_tolerance},
                  {"divergent_cancellation_provenance",
                   diagnostics.divergent_cancellation_provenance},
                  {"primitive_evaluations",
                   diagnostics.primitive_evaluations},
                  {"primitive_component_applications",
                   diagnostics.primitive_component_applications},
                  {"primitive_component_reuses",
                   diagnostics.primitive_component_reuses},
                  {"has_center_endpoint",
                   diagnostics.has_center_endpoint},
                  {"tail_certificate_requested",
                   diagnostics.tail_certificate_requested},
                  {"tail_certificate_status",
                   diagnostics.tail_certificate_status},
                  {"tail_witness_radius_exact",
                   diagnostics.tail_witness_radius_exact.empty()
                       ? json::value(nullptr)
                       : json::value(
                             diagnostics.tail_witness_radius_exact)},
                  {"detail", diagnostics.detail}}}}},
        {"elapsed_ms", elapsed_ms_},
        {"runtime_stats",
         json::object{{"exports", exports_}, {"export_ms", export_ms_}}}};
  }

  void restore_runtime_stats(std::uint64_t exports, double export_ms) {
    std::lock_guard<std::mutex> lock(stats_mutex_);
    exports_ = exports;
    export_ms_ = export_ms;
  }

 private:
  static void validate_transport_reference(
      const json::value& raw,
      const std::shared_ptr<StoredTransportArmState>& owner,
      const char* label, bool compact_v2 = false) {
    const auto& reference = as_object(raw, label);
    if (compact_v2)
      require_exact_keys(reference, {"handle", "checkpoint_identity"},
                         label);
    else
      require_exact_keys(reference,
                         {"handle", "checkpoint_identity",
                          "provenance_identity"}, label);
    if (!owner || required_string(reference, "handle") != owner->handle() ||
        required_string(reference, "checkpoint_identity") !=
            owner->checkpoint_identity() ||
        (!compact_v2 &&
         required_string(reference, "provenance_identity") !=
             owner->provenance_identity()))
      throw std::invalid_argument(
          std::string(label) + " is stale");
  }

  void validate_aggregate_ownership() const {
    if (!aggregate_provenance_.has_value() || plan_owners_.empty() ||
        std::any_of(plan_owners_.begin(), plan_owners_.end(),
                    [](const auto& owner) { return owner == nullptr; }) ||
        std::any_of(local_owners_.begin(), local_owners_.end(),
                    [](const auto& owner) { return owner == nullptr; }) ||
        std::any_of(transport_owners_.begin(), transport_owners_.end(),
                    [](const auto& owner) { return owner == nullptr; }))
      throw std::invalid_argument(
          "retained line aggregate requires its exact strong owners");
    if (json::serialize(canonical_json_value(*aggregate_provenance_)) !=
        provenance_identity_)
      throw std::invalid_argument(
          "retained line aggregate provenance identity is inconsistent");
    const auto schema = required_string(
        *aggregate_provenance_, "schema");
    const auto& source = as_object(
        aggregate_provenance_->at("source"),
        "retained line aggregate source");
    const auto& aggregate = as_object(
        aggregate_provenance_->at("aggregate"),
        "retained line aggregate recipe");
    if (schema ==
            "diffexp2-retained-native-transport-observable-line-v1" ||
        schema ==
            "diffexp2-retained-native-transport-observable-line-v2") {
      const bool compact_v2 = schema.ends_with("-v2");
      if (plan_owners_.size() != 1 || transport_owners_.size() != 1 ||
          !local_owners_.empty())
        throw std::invalid_argument(
            "compact transport line requires one plan and state owner");
      validate_transport_reference(
          aggregate.at("transport_state"), transport_owners_.front(),
          "retained line aggregate transport reference", compact_v2);
      return;
    }
    if (schema ==
            "diffexp2-retained-native-transport-pair-observable-line-v1" ||
        schema ==
            "diffexp2-retained-native-transport-pair-observable-line-v2") {
      const bool compact_v2 = schema.ends_with("-v2");
      if (plan_owners_.size() != 2 || transport_owners_.size() != 2 ||
          !local_owners_.empty())
        throw std::invalid_argument(
            "compact transport-pair line requires two ordered plan and state owners");
      require_transport_pair_compatibility(
          transport_owners_[0], transport_owners_[1],
          transport_owners_[0]->anchor_owner()->scalar_domain());
      if (transport_owners_[0]->plan_owner().get() !=
              plan_owners_[0].get() ||
          transport_owners_[1]->plan_owner().get() !=
              plan_owners_[1].get())
        throw std::invalid_argument(
            "compact transport-pair plans differ from their ordered state owners");
      const auto& lower_source = as_object(
          source.at("lower"), "compact transport-pair lower source");
      const auto& upper_source = as_object(
          source.at("upper"), "compact transport-pair upper source");
      validate_transport_reference(
          lower_source.at("transport_state"), transport_owners_[0],
          "compact transport-pair lower state reference", compact_v2);
      validate_transport_reference(
          upper_source.at("transport_state"), transport_owners_[1],
          "compact transport-pair upper state reference", compact_v2);
      const auto& lower_recipe = as_object(
          aggregate.at("lower"), "compact transport-pair lower recipe");
      const auto& upper_recipe = as_object(
          aggregate.at("upper"), "compact transport-pair upper recipe");
      validate_transport_reference(
          lower_recipe.at("transport_state"), transport_owners_[0],
          "compact transport-pair lower recipe state", compact_v2);
      validate_transport_reference(
          upper_recipe.at("transport_state"), transport_owners_[1],
          "compact transport-pair upper recipe state", compact_v2);
      return;
    }
    if (schema != "diffexp2-retained-native-line-aggregate-v1" ||
        plan_owners_.size() != 1 || !transport_owners_.empty() ||
        local_owners_.empty())
      throw std::invalid_argument(
          "retained line aggregate ownership differs from its provenance schema");
  }

  std::string handle_;
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::string arm_;
  std::size_t tile_index_ = 0;
  json::object interval_;
  std::string source_checkpoint_;
  StoredLineIntegral result_;
  double elapsed_ms_ = 0.0;
  std::vector<std::shared_ptr<StoredTilePlan>> plan_owners_;
  std::vector<std::shared_ptr<StoredLocalBase>> local_owners_;
  std::optional<json::object> aggregate_provenance_;
  std::vector<std::shared_ptr<StoredTransportArmState>> transport_owners_;
  mutable std::mutex stats_mutex_;
  std::uint64_t exports_ = 0;
  double export_ms_ = 0.0;
};
