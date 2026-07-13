Rational parse_exact_path_rational(const json::value& raw,
                                   const char* label) {
  if (!raw.is_string() || raw.as_string().empty())
    throw std::invalid_argument(std::string(label) +
                                " must be a nonempty exact rational string");
  try {
    return Rational(std::string(raw.as_string()));
  } catch (const std::invalid_argument&) {
    throw std::invalid_argument(std::string(label) +
                                " is not an exact rational number");
  }
}

std::vector<Rational> parse_exact_path_points(const json::value& raw,
                                              const char* label) {
  std::vector<Rational> points;
  for (const auto& value : as_array(raw, label))
    points.push_back(parse_exact_path_rational(value, label));
  return points;
}

ExactPathTopology parse_exact_path_topology(const json::value& raw) {
  const auto& object = as_object(raw, "native exact path topology");
  require_exact_keys(object,
      {"singular_points", "boundary_points", "complex_projections",
       "branch_sheets"}, "native exact path topology");
  ExactPathTopology topology;
  topology.singular_points = parse_exact_path_points(
      object.at("singular_points"), "path singular point");
  topology.boundary_points = parse_exact_path_points(
      object.at("boundary_points"), "path boundary point");
  for (const auto& raw_projection : as_array(
           object.at("complex_projections"), "path complex projections")) {
    const auto& projection = as_object(
        raw_projection, "path complex projection");
    require_exact_keys(projection,
        {"source_identity", "real_part_exact",
         "imaginary_magnitude_exact", "retain_minus_imaginary",
         "retain_real_part", "retain_plus_imaginary"},
        "path complex projection");
    topology.complex_projections.push_back(ExactComplexProjection{
        required_string(projection, "source_identity"),
        parse_exact_path_rational(projection.at("real_part_exact"),
                                  "projection real part"),
        parse_exact_path_rational(
            projection.at("imaginary_magnitude_exact"),
            "projection imaginary magnitude"),
        projection.at("retain_minus_imaginary").as_bool(),
        projection.at("retain_real_part").as_bool(),
        projection.at("retain_plus_imaginary").as_bool()});
  }
  for (const auto& raw_sheet : as_array(
           object.at("branch_sheets"), "path branch sheets")) {
    const auto& sheet = as_object(raw_sheet, "path branch sheet");
    require_exact_keys(sheet, {"factor_exact", "sign"},
                       "path branch sheet");
    topology.branch_sheets.push_back(ExactBranchSheet{
        required_string(sheet, "factor_exact"),
        as_i32(sheet.at("sign"), "path branch sign")});
  }
  return topology;
}

std::string retained_plan_owner_handle(
    const RetainedPlanChartBinding::Owner& owner) {
  return std::visit(
      [](const auto& typed) {
        if (typed == nullptr)
          throw std::logic_error(
              "native tile planning received a null retained owner");
        return typed->handle();
      },
      owner);
}

std::string retained_plan_owner_identity(
    const RetainedPlanChartBinding::Owner& owner) {
  return std::visit(
      [](const auto& typed) {
        if (typed == nullptr)
          throw std::logic_error(
              "native tile planning received a null retained owner");
        return typed->exact_identity();
      },
      owner);
}

std::string retained_plan_owner_geometry_record(
    const RetainedPlanChartBinding::Owner& owner) {
  return std::visit(
      [](const auto& typed) -> std::string {
        using Owner = typename std::decay_t<decltype(typed)>::element_type;
        if (typed == nullptr)
          throw std::logic_error(
              "native tile planning received a null retained owner");
        if constexpr (std::is_same_v<Owner, PreparedChartBase>) {
          if (!typed->geometry_record().has_value())
            throw std::invalid_argument(
                "native tile planning requires retained exact chart geometry");
          return *typed->geometry_record();
        } else {
          return typed->geometry_record();
        }
      },
      owner);
}

RetainedPlanChartBinding bind_plan_chart(
    const RetainedPlanChartBinding::Owner& owner,
    const ExactPathTopology& topology) {
  const auto handle = retained_plan_owner_handle(owner);
  const auto exact_identity = retained_plan_owner_identity(owner);
  const auto geometry_value = json::parse(
      retained_plan_owner_geometry_record(owner));
  const auto& geometry = as_object(
      geometry_value, "retained native tile chart geometry");
  if (geometry.at("infinite_radius").as_bool())
    throw std::invalid_argument(
        "native exact tile planning currently requires finite chart radii");

  RetainedPlanChartBinding binding;
  binding.handle = handle;
  binding.exact_identity = exact_identity;
  binding.owner = owner;
  binding.geometry.identity = exact_identity;
  binding.geometry.center = parse_exact_path_rational(
      geometry.at("center_exact"), "tile chart center");
  binding.geometry.scale = parse_exact_path_rational(
      geometry.at("scale_exact"), "tile chart scale");
  binding.geometry.radius = parse_exact_path_rational(
      geometry.at("radius_exact"), "tile chart radius");
  binding.geometry.singular_center = std::any_of(
      topology.singular_points.begin(), topology.singular_points.end(),
      [&](const Rational& point) { return point == binding.geometry.center; });
  for (const auto& raw_prescription : as_array(
           geometry.at("prescriptions"), "tile chart prescriptions")) {
    const auto& prescription = as_object(
        raw_prescription, "tile chart prescription");
    binding.prescriptions.push_back(Prescription{
        required_string(prescription, "factor_exact"),
        as_i32(prescription.at("sign"), "tile prescription sign"),
        as_u32(prescription.at("multiplicity"),
               "tile prescription multiplicity"),
        as_i32(prescription.at("leading_coefficient_sign"),
               "tile prescription leading coefficient sign")});
  }
  // Every exact prepared prescription must be represented on the path.  The
  // topology may additionally retain factors which are inactive in this
  // chart, but it may never silently alter or drop an active sheet.
  for (const auto& prescription : binding.prescriptions) {
    const auto found = std::find_if(
        topology.branch_sheets.begin(), topology.branch_sheets.end(),
        [&](const ExactBranchSheet& sheet) {
          return sheet.factor_exact == prescription.factor_exact;
        });
    if (found == topology.branch_sheets.end() ||
        found->imaginary_sign != prescription.sign)
      throw std::invalid_argument(
          "native tile topology does not reproduce a prepared chart branch prescription");
  }
  (void)exact_plan_rim(binding.prescriptions, binding.geometry.scale);
  return binding;
}

std::vector<std::string> parse_plan_chart_handles(const json::object& arm) {
  std::vector<std::string> handles;
  for (const auto& raw : as_array(arm.at("charts"), "tile arm charts")) {
    if (!raw.is_string() || raw.as_string().empty())
      throw std::invalid_argument(
          "native tile arm chart handles must be nonempty strings");
    handles.emplace_back(raw.as_string());
  }
  if (handles.empty())
    throw std::invalid_argument("native tile arm requires at least one chart");
  return handles;
}

std::pair<ExactArmRequest, std::vector<RetainedPlanChartBinding>>
parse_retained_arm_request(
    const json::object& arm,
    const std::vector<RetainedPlanChartBinding::Owner>& charts) {
  require_exact_keys(arm, {"from_exact", "to_exact", "charts", "topology"},
                     "native tile arm");
  const auto handles = parse_plan_chart_handles(arm);
  if (handles.size() != charts.size())
    throw std::invalid_argument(
        "resolved native tile chart count differs from its request");
  ExactArmRequest request;
  request.from = parse_exact_path_rational(arm.at("from_exact"),
                                           "tile arm start");
  request.to = parse_exact_path_rational(arm.at("to_exact"),
                                         "tile arm end");
  request.topology = parse_exact_path_topology(arm.at("topology"));
  std::vector<RetainedPlanChartBinding> bindings;
  bindings.reserve(charts.size());
  request.charts.reserve(charts.size());
  for (std::size_t index = 0; index < charts.size(); ++index) {
    if (retained_plan_owner_handle(charts[index]) != handles[index])
      throw std::logic_error("resolved tile chart handle changed");
    auto binding = bind_plan_chart(charts[index], request.topology);
    request.charts.push_back(binding.geometry);
    bindings.push_back(std::move(binding));
  }
  return {std::move(request), std::move(bindings)};
}

std::shared_ptr<StoredTilePlan> build_tile_plan(
    const std::string& handle, const json::object& request,
    const std::vector<RetainedPlanChartBinding::Owner>& lower_charts,
    const std::vector<RetainedPlanChartBinding::Owner>& upper_charts) {
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  const auto division_order = as_u32(
      request.at("division_order"), "native tile division order");
  auto [lower_request, lower_bindings] = parse_retained_arm_request(
      as_object(request.at("lower"), "lower native tile arm"), lower_charts);
  auto [upper_request, upper_bindings] = parse_retained_arm_request(
      as_object(request.at("upper"), "upper native tile arm"), upper_charts);
  if (lower_bindings.front().handle != upper_bindings.front().handle)
    throw std::invalid_argument(
        "independent native tile arms must share one retained anchor chart");

  ExactPathPlanOptions options;
  options.division_order = division_order;
  const auto started = std::chrono::steady_clock::now();
  auto exact = plan_exact_independent_arms(
      lower_request, upper_request, options);
  RetainedArmPlan lower{std::move(exact.lower),
                        std::move(lower_bindings)};
  RetainedArmPlan upper{std::move(exact.upper),
                        std::move(upper_bindings)};
  json::object provenance{
      {"schema", "diffexp2-retained-exact-independent-arm-tile-plan-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"division_order", division_order},
      {"lower", encode_retained_arm(lower)},
      {"upper", encode_retained_arm(upper)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredTilePlan>(
      handle, checkpoint_identity, provenance_identity, division_order,
      std::move(lower), std::move(upper), elapsed);
}

std::shared_ptr<StoredTilePlan> build_single_arm_tile_plan(
    const std::string& handle, const json::object& request,
    const std::vector<RetainedPlanChartBinding::Owner>& charts) {
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  if (checkpoint_identity.empty())
    throw std::invalid_argument(
        "native single-arm tile-plan checkpoint identity cannot be empty");
  const auto division_order = as_u32(
      request.at("division_order"), "native single-arm tile division order");
  auto [arm_request, bindings] = parse_retained_arm_request(
      as_object(request.at("arm"), "native single tile arm"), charts);
  ExactPathPlanOptions options;
  options.division_order = division_order;
  const auto started = std::chrono::steady_clock::now();
  auto exact = plan_exact_arm(arm_request, options);
  const std::string arm_name = exact.direction < 0 ? "lower" : "upper";
  RetainedArmPlan retained{std::move(exact), std::move(bindings)};
  json::object provenance{
      {"schema", kRetainedSingleArmTilePlanProvenanceSchema},
      {"checkpoint_identity", checkpoint_identity},
      {"division_order", division_order},
      {"arm_name", arm_name},
      {"arm", encode_retained_arm(retained)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredTilePlan>(
      handle, checkpoint_identity, provenance_identity, division_order,
      arm_name, std::move(retained), elapsed);
}

json::value optional_plan_rim_json(
    const std::optional<std::int32_t>& rim) {
  return rim.has_value() ? json::value(*rim) : json::value(nullptr);
}

struct ResolvedPlannedEndpointBinding {
  json::object source;
  std::int32_t approach_direction = 0;
  std::optional<std::int32_t> rim;
};

ResolvedPlannedEndpointBinding resolve_planned_endpoint_binding(
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::string& arm_name,
    const std::shared_ptr<StoredLocalBase>& local) {
  if (!plan || !local)
    throw std::invalid_argument(
        "plan-bound endpoint evaluation requires retained plan and local owners");
  const auto& arm = plan->arm(arm_name);
  if (arm.exact.tiles.empty() || arm.charts.empty())
    throw std::invalid_argument(
        "plan-bound endpoint arm has no final tile/chart");
  const auto final_tile_index = arm.exact.tiles.size() - 1;
  const auto& final_tile = arm.exact.tiles.back();
  if (final_tile.chart >= arm.charts.size())
    throw std::logic_error(
        "plan-bound endpoint final tile has an invalid chart index");
  const auto& final_chart = arm.charts[final_tile.chart];
  if (!(final_tile.physical_end == arm.exact.to) ||
      !(final_chart.geometry.center == arm.exact.to) ||
      !final_tile.local_end.is_zero())
    throw std::invalid_argument(
        "plan-bound endpoint requires the final tile to end at the exact center of its retained final chart");
  if (arm.exact.direction != -1 && arm.exact.direction != 1)
    throw std::logic_error(
        "plan-bound endpoint arm has an invalid exact direction");
  if (final_chart.geometry.scale.is_zero())
    throw std::invalid_argument(
        "plan-bound endpoint final chart has a zero exact scale");
  if (local->source_chart() != final_chart.handle)
    throw std::invalid_argument(
        "plan-bound endpoint local does not name the retained final chart");
  local->require_exact_plan_binding(
      final_chart.geometry, final_chart.prescriptions,
      "plan-bound endpoint final local");

  ResolvedPlannedEndpointBinding resolved;
  resolved.approach_direction =
      -arm.exact.direction * final_chart.geometry.scale.sign();
  resolved.rim = exact_plan_rim(
      final_chart.prescriptions, final_chart.geometry.scale);
  resolved.source = json::object{
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
      {"tile_plan_provenance_identity", plan->provenance_identity()},
      {"arm", arm_name},
      {"endpoint_exact", arm.exact.to.str()},
      {"direction", arm.exact.direction},
      {"final_tile", final_tile_index},
      {"final_chart_index", final_tile.chart},
      {"final_chart", final_chart.handle},
      {"final_chart_identity", final_chart.exact_identity},
      {"local", local->handle()},
      {"chart", local->source_chart()},
      {"source_operator_identity", local->source_operator_identity()},
      {"checkpoint_identity", local->checkpoint_identity()},
      {"coefficient_domain", local->scalar_domain()},
      {"prescriptions", encode_plan_prescriptions(
           final_chart.prescriptions)}};
  return resolved;
}

json::object planned_endpoint_provenance(
    const std::string& checkpoint_identity,
    const ResolvedPlannedEndpointBinding& binding,
    const std::string& cancellation_mode,
    const json::object& analytic_metadata) {
  return json::object{
      {"schema",
       "diffexp2-retained-native-plan-bound-endpoint-sector-limit-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"source", binding.source},
      {"approach_direction", binding.approach_direction},
      {"rim", optional_plan_rim_json(binding.rim)},
      {"cancellation", json::object{{"mode", cancellation_mode}}},
      {"analytic_metadata", analytic_metadata}};
}

std::shared_ptr<StoredEndpointResult> build_planned_endpoint_limit(
    const std::string& endpoint_handle, const json::object& request,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& local) {
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  const auto expected_plan_checkpoint = required_string(
      request, "tile_plan_checkpoint_identity");
  const auto expected_source_checkpoint = required_string(
      request, "source_checkpoint_identity");
  if (checkpoint_identity.empty() || expected_plan_checkpoint.empty() ||
      expected_source_checkpoint.empty())
    throw std::invalid_argument(
        "plan-bound endpoint checkpoint identities must be nonempty");
  if (expected_plan_checkpoint != plan->checkpoint_identity())
    throw std::invalid_argument(
        "plan-bound endpoint tile-plan checkpoint identity is stale or mismatched");
  if (expected_source_checkpoint != local->checkpoint_identity())
    throw std::invalid_argument(
        "plan-bound endpoint source checkpoint identity is stale or mismatched");
  const auto arm_name = required_string(request, "arm");
  const auto binding = resolve_planned_endpoint_binding(
      plan, arm_name, local);
  const auto cancellation = parse_endpoint_cancellation_policy(request);
  auto analytic_metadata = local->exact_analytic_metadata();
  const auto provenance = planned_endpoint_provenance(
      checkpoint_identity, binding, cancellation.cancellation_mode,
      analytic_metadata);
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));

  EndpointLimitOptions options;
  options.approach_direction = binding.approach_direction;
  options.imaginary_sign = binding.rim;
  options.allow_certified_numeric_cancellation =
      cancellation.allow_certified_numeric_cancellation;
  const auto started = std::chrono::steady_clock::now();
  auto result = local->endpoint_limit(options);
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredEndpointResult>(
      endpoint_handle, checkpoint_identity, provenance_identity,
      local->handle(), local->source_chart(),
      local->source_operator_identity(), expected_source_checkpoint,
      local->scalar_domain(), binding.approach_direction, std::nullopt,
      cancellation.cancellation_mode, std::move(analytic_metadata),
      std::move(result), elapsed, binding.source, binding.rim,
      plan, local);
}

json::object planned_match_handoff_record(
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::string& arm_name, std::size_t match_index,
    const std::vector<std::string>& basis_handles,
    const std::vector<std::shared_ptr<StoredLocalBase>>& basis,
    const std::string& incoming_handle,
    const std::shared_ptr<StoredLocalBase>& incoming,
    const std::string& result_checkpoint,
    const json::object& native_summary,
    bool compact_plan_reference = false) {
  if (!plan || !incoming || basis.empty() ||
      basis.size() != basis_handles.size())
    throw std::invalid_argument(
        "planned match handoff record requires its complete owner set");
  const auto& arm = plan->arm(arm_name);
  if (match_index >= arm.exact.matches.size())
    throw std::invalid_argument(
        "planned match handoff index is outside its retained arm");
  const auto& exact_match = arm.exact.matches[match_index];
  const auto& producing = arm.charts.at(exact_match.producing_chart);
  const auto& receiving = arm.charts.at(exact_match.receiving_chart);
  const auto producing_rim = exact_plan_rim(
      producing.prescriptions, producing.geometry.scale);
  const auto receiving_rim = exact_plan_rim(
      receiving.prescriptions, receiving.geometry.scale);
  json::array basis_sources;
  basis_sources.reserve(basis.size());
  for (std::size_t column = 0; column < basis.size(); ++column)
    basis_sources.push_back(json::object{
        {"column", column}, {"local", basis_handles[column]},
        {"checkpoint_identity", basis[column]->checkpoint_identity()},
        {"source_operator_identity",
         basis[column]->source_operator_identity()}});
  if (compact_plan_reference)
    return json::object{
        {"schema", "diffexp2-retained-exact-plan-match-hop-v2"},
        {"tile_plan", plan->handle()},
        {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
        {"arm", arm_name}, {"match", match_index},
        {"geometry", encode_plan_match(arm, match_index)},
        {"producing", json::object{
             {"tile", match_index},
             {"chart", producing.handle},
             {"chart_identity", producing.exact_identity},
             {"local_point_exact", exact_match.producing_local.str()},
             {"effective_rim", optional_plan_rim_json(producing_rim)},
             {"prescriptions",
              encode_plan_prescriptions(producing.prescriptions)},
             {"incoming", json::object{
                  {"local", incoming_handle},
                  {"checkpoint_identity", incoming->checkpoint_identity()},
                  {"source_operator_identity",
                   incoming->source_operator_identity()}}}}},
        {"receiving", json::object{
             {"tile", match_index + 1},
             {"chart", receiving.handle},
             {"chart_identity", receiving.exact_identity},
             {"local_point_exact", exact_match.receiving_local.str()},
             {"effective_rim", optional_plan_rim_json(receiving_rim)},
             {"prescriptions",
              encode_plan_prescriptions(receiving.prescriptions)},
             {"basis", std::move(basis_sources)}}},
        {"result_checkpoint_identity", result_checkpoint},
        {"advance", json::object{
             {"scope", "single-match-handoff"},
             {"state", "retained-receiving-basis-weights"},
             {"source_tile", match_index},
             {"receiving_tile", match_index + 1},
             {"whole_arm_complete", false}}}};
  return json::object{
      {"schema", "diffexp2-retained-exact-plan-match-hop-v1"},
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
      {"tile_plan_provenance_identity", plan->provenance_identity()},
      {"arm", arm_name}, {"match", match_index},
      {"geometry", encode_plan_match(arm, match_index)},
      {"producing", json::object{
           {"tile", match_index},
           {"chart", producing.handle},
           {"chart_identity", producing.exact_identity},
           {"local_point_exact", exact_match.producing_local.str()},
           {"effective_rim", optional_plan_rim_json(producing_rim)},
           {"prescriptions",
            encode_plan_prescriptions(producing.prescriptions)},
           {"incoming", json::object{
                {"local", incoming_handle},
                {"checkpoint_identity", incoming->checkpoint_identity()},
                {"source_operator_identity",
                 incoming->source_operator_identity()}}}}},
      {"receiving", json::object{
           {"tile", match_index + 1},
           {"chart", receiving.handle},
           {"chart_identity", receiving.exact_identity},
           {"local_point_exact", exact_match.receiving_local.str()},
           {"effective_rim", optional_plan_rim_json(receiving_rim)},
           {"prescriptions",
            encode_plan_prescriptions(receiving.prescriptions)},
           {"basis", std::move(basis_sources)}}},
      {"result_checkpoint_identity", result_checkpoint},
      {"native_match_provenance_identity",
       required_string(native_summary, "provenance_identity")},
      {"advance", json::object{
           {"scope", "single-match-handoff"},
           {"state", "retained-receiving-basis-weights"},
           {"source_tile", match_index},
           {"receiving_tile", match_index + 1},
           {"whole_arm_complete", false}}}};
}

struct NativeAcbSaturationBinding {
  std::string request_key;
  json::object request;
};

NativeAcbSaturationBinding native_acb_saturation_binding(
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::string& session_configuration_identity,
    const std::string& arm,
    std::size_t match_index, const std::string& match_checkpoint_identity,
    bool compact_plan_reference = false);

std::shared_ptr<StoredPlannedMatchHop> build_planned_match_hop(
    const std::string& match_handle, const json::object& request,
    const std::string& domain, slong precision_bits,
    const std::string& active_session_configuration_identity,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::vector<std::string>& basis_handles,
    const std::vector<std::shared_ptr<StoredLocalBase>>& basis,
    const std::string& incoming_handle,
    const std::shared_ptr<StoredLocalBase>& incoming,
    bool compact_plan_reference = false) {
  const auto started = std::chrono::steady_clock::now();
  const auto arm_name = required_string(request, "arm");
  const auto match_index = static_cast<std::size_t>(
      as_u64(request.at("match"), "planned local match index"));
  const auto& arm = plan->arm(arm_name);
  if (match_index >= arm.exact.matches.size())
    throw std::invalid_argument(
        "planned local match index is out of range");
  const auto& exact_match = arm.exact.matches[match_index];
  const auto& producing = arm.charts.at(exact_match.producing_chart);
  const auto& receiving = arm.charts.at(exact_match.receiving_chart);
  if (basis.empty())
    throw std::invalid_argument("planned local match basis cannot be empty");

  // The plan, rather than the caller, binds every coordinate, chart,
  // prescription, rim and source checkpoint passed to the existing matching
  // kernels.  Locals must reproduce the prepared chart snapshot exactly.
  incoming->require_exact_plan_binding(
      producing.geometry, producing.prescriptions,
      "planned incoming " + incoming_handle);
  for (std::size_t column = 0; column < basis.size(); ++column)
    basis[column]->require_exact_plan_binding(
        receiving.geometry, receiving.prescriptions,
        "planned basis " + basis_handles[column]);

  json::array basis_checkpoints;
  basis_checkpoints.reserve(basis.size());
  for (const auto& local : basis)
    basis_checkpoints.emplace_back(local->checkpoint_identity());
  const auto result_checkpoint = required_string(
      request, "checkpoint_identity");
  if (result_checkpoint.empty())
    throw std::invalid_argument(
        "planned local match checkpoint identity cannot be empty");

  json::object kernel_request{
      {"basis", [&]() {
         json::array values;
         for (const auto& handle : basis_handles) values.emplace_back(handle);
         return values;
       }()},
      {"incoming", incoming_handle},
      {"basis_chart", receiving.handle},
      {"incoming_chart", producing.handle},
      {"basis_point", json::object{
           {"exact", exact_match.receiving_local.str()}}},
      {"incoming_point", json::object{
           {"exact", exact_match.producing_local.str()}}},
      {"epsilon", request.at("epsilon")},
      {"basis_checkpoint_identities", std::move(basis_checkpoints)},
      {"incoming_checkpoint_identity", incoming->checkpoint_identity()},
      {"checkpoint_identity", result_checkpoint}};

  const auto producing_rim = exact_plan_rim(
      producing.prescriptions, producing.geometry.scale);
  const auto receiving_rim = exact_plan_rim(
      receiving.prescriptions, receiving.geometry.scale);
  std::shared_ptr<StoredMatchBase> native_match;
  std::optional<json::object> expected_singular_request;
  if (domain == "rational") {
    native_match = build_exact_regular_match(
        match_handle, kernel_request, basis_handles, basis, incoming_handle,
        incoming);
  } else if (domain == "acb") {
    kernel_request["basis_imaginary_sign"] =
        optional_plan_rim_json(receiving_rim);
    kernel_request["incoming_imaginary_sign"] =
        optional_plan_rim_json(producing_rim);
    kernel_request["refinement"] = request.at("refinement");
    if (const auto* exact = request.if_contains("exact_lattice"))
      kernel_request["exact_lattice"] = *exact;
    else if (const auto* native =
                 request.if_contains("native_unit_saturation"))
      kernel_request["native_unit_saturation"] = *native;
    else if (const auto* singular =
                 request.if_contains("native_singular_scc_saturation")) {
      const auto expected = native_acb_saturation_binding(
          plan, active_session_configuration_identity, arm_name,
          match_index, result_checkpoint, compact_plan_reference);
      if (expected.request_key != "native_singular_scc_saturation" ||
          json::serialize(canonical_json_value(*singular)) !=
              json::serialize(canonical_json_value(expected.request)))
        throw std::invalid_argument(
            "planned singular-SCC Acb saturation request does not match the retained receiving SCC");
      kernel_request["native_singular_scc_saturation"] = *singular;
      expected_singular_request = std::move(expected.request);
    } else
      throw std::invalid_argument(
          "planned Acb matching requires an exact lattice, ordinary native unit-leading request, or singular-SCC valuation-zero request");
    native_match = build_refined_acb_match(
        match_handle, kernel_request, basis_handles, basis, incoming_handle,
        incoming, precision_bits, active_session_configuration_identity,
        expected_singular_request);
  } else {
    throw std::invalid_argument(
        "plan-driven local matching requires rational or Acb coefficients");
  }

  const auto native_summary = native_match->summary();
  if (required_string(native_summary, "physical_match_point_exact") !=
      exact_match.physical.str())
    throw std::logic_error(
        "native local match physical point differs from its retained exact plan");

  auto handoff = planned_match_handoff_record(
      plan, arm_name, match_index, basis_handles, basis, incoming_handle,
      incoming, result_checkpoint, native_summary,
      compact_plan_reference);
  const auto provenance_identity = json::serialize(
      canonical_json_value(handoff));
  const auto elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredPlannedMatchHop>(
      std::move(native_match), result_checkpoint, provenance_identity,
      std::move(handoff), elapsed_ms, plan, basis, incoming);
}

using BoundedDivergentCancellation =
    StoredLineIntegrationOptions::BoundedDivergentCancellation;

BoundedDivergentCancellation parse_bounded_divergent_cancellation(
    const json::value& raw, const char* label) {
  const auto& policy = as_object(raw, label);
  require_exact_keys(policy,
                     {"mode", "relative_tolerance", "provenance"},
                     label);
  if (required_string(policy, "mode") != "bounded-relative-acb")
    throw std::invalid_argument(
        std::string(label) +
        " mode must be exactly bounded-relative-acb");
  const auto tolerance_text = required_string(
      policy, "relative_tolerance");
  const auto provenance = required_string(policy, "provenance");
  if (tolerance_text.empty() || provenance.empty())
    throw std::invalid_argument(
        std::string(label) +
        " requires a nonempty tolerance and producer provenance");
  auto tolerance = Magnitude::decimal(tolerance_text);
  if (!tolerance.is_finite() || tolerance.is_zero() ||
      Magnitude::one() <= tolerance)
    throw std::invalid_argument(
        std::string(label) +
        " tolerance must be finite and strictly between zero and one");
  return {std::move(tolerance), tolerance_text, provenance};
}

json::object encode_bounded_divergent_cancellation(
    const BoundedDivergentCancellation& policy) {
  return json::object{
      {"mode", "bounded-relative-acb"},
      {"relative_tolerance", policy.relative_tolerance_text},
      {"provenance", policy.provenance}};
}

std::shared_ptr<StoredLineResult> build_planned_line_result(
    const std::string& handle, const json::object& request,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& local) {
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  const auto source_checkpoint = required_string(
      request, "source_checkpoint_identity");
  if (source_checkpoint != local->checkpoint_identity())
    throw std::invalid_argument(
        "planned line source checkpoint identity differs from its retained local");
  const auto expected_plan_checkpoint = required_string(
      request, "tile_plan_checkpoint_identity");
  if (expected_plan_checkpoint != plan->checkpoint_identity())
    throw std::invalid_argument(
        "planned line tile-plan checkpoint identity differs from retained state");
  const auto arm_name = required_string(request, "arm");
  const auto tile_index = static_cast<std::size_t>(
      as_u64(request.at("tile"), "planned line tile index"));
  const auto& arm = plan->arm(arm_name);
  if (tile_index >= arm.exact.tiles.size())
    throw std::invalid_argument("planned line tile index is out of range");
  const auto& tile = arm.exact.tiles[tile_index];
  const auto& binding = arm.charts.at(tile.chart);
  if (local->source_chart() != binding.handle)
    throw std::invalid_argument(
        "planned line local does not belong to the tile's retained chart");
  const auto& epsilon = as_object(request.at("epsilon"),
                                  "planned line epsilon window");
  const auto delivered = parse_epsilon_window(
      epsilon, "planned line epsilon window");
  const bool certify_tail =
      request.if_contains("certify_tail") != nullptr &&
      request.at("certify_tail").as_bool();
  std::optional<BoundedDivergentCancellation> divergent_cancellation;
  if (const auto* policy = request.if_contains("divergent_cancellation"))
    divergent_cancellation = parse_bounded_divergent_cancellation(
        *policy, "planned line divergent-cancellation policy");
  auto interval = encode_plan_tile(arm, tile_index);
  const auto rim = exact_plan_rim(
      binding.prescriptions, binding.geometry.scale);
  const auto started = std::chrono::steady_clock::now();
  StoredLineIntegral result;
  try {
    result = local->integrate_planned_line(
        binding.geometry, binding.prescriptions, tile.local_begin,
        tile.local_end, delivered, rim, certify_tail,
        divergent_cancellation);
  } catch (const NativeIntegrationError& error) {
    std::ostringstream detail;
    detail << error.what() << "; arm=" << arm_name
           << "; tile=" << tile_index
           << "; physical_interval=[" << tile.physical_begin.str()
           << "," << tile.physical_end.str() << "]"
           << "; local_interval=[" << tile.local_begin.str()
           << "," << tile.local_end.str() << "]"
           << "; chart=" << binding.handle
           << "; chart_center=" << binding.geometry.center.str()
           << "; chart_scale=" << binding.geometry.scale.str()
           << "; rim="
           << (rim.has_value() ? std::to_string(*rim) : "none")
           << "; source_local=" << local->handle();
    NativeIntegrationError contextual(error.code, error.id, detail.str());
    contextual.absolute_power = error.absolute_power;
    contextual.log_power = error.log_power;
    contextual.epsilon_power = error.epsilon_power;
    contextual.component = error.component;
    throw contextual;
  }
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  json::object provenance{
      {"schema",
       "diffexp2-retained-native-physical-tile-integral-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", expected_plan_checkpoint},
      {"arm", arm_name}, {"tile", tile_index},
      {"interval", interval},
      {"source", json::object{
           {"local", local->handle()}, {"chart", local->source_chart()},
           {"checkpoint_identity", source_checkpoint}}},
      {"epsilon", json::object{{"min", delivered.min_power},
                                {"max", delivered.complete_max}}},
      {"tail_certificate_requested", certify_tail},
      {"tail_certificate_status",
       result.diagnostics.tail_certificate_status},
      {"tail_witness_radius_exact",
       result.diagnostics.tail_witness_radius_exact.empty()
           ? json::value(nullptr)
           : json::value(result.diagnostics.tail_witness_radius_exact)},
      {"divergent_cancellation",
       divergent_cancellation.has_value()
           ? json::value(encode_bounded_divergent_cancellation(
                 *divergent_cancellation))
           : json::value(json::object{{"mode", "exact-singleton"}})},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)},
      {"error_provenance", result.value.error.provenance}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  return std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, arm_name, tile_index,
      std::move(interval), source_checkpoint, std::move(result), elapsed,
      plan, local);
}

std::size_t checked_diagnostic_sum(std::size_t left, std::size_t right,
                                   const char* label) {
  if (right > std::numeric_limits<std::size_t>::max() - left)
    throw std::overflow_error(std::string(label) + " overflow");
  return left + right;
}

StoredLineIntegral aggregate_retained_lines(
    const std::vector<std::shared_ptr<StoredLineResult>>& components,
    const std::vector<std::int32_t>& signs, const std::string& detail) {
  if (components.empty() || components.size() != signs.size())
    throw std::invalid_argument(
        "native line aggregation requires a nonempty signed component list");
  for (const auto sign : signs)
    if (sign != -1 && sign != 1)
      throw std::invalid_argument(
          "native line aggregation signs must be +1 or -1");

  auto epsilon_min = components.front()->result().value.epsilon.min_power;
  auto epsilon_max =
      components.front()->result().value.epsilon.complete_max;
  const auto dimension = components.front()->result().value.dimension;
  for (const auto& component : components) {
    const auto& value = component->result().value;
    if (value.dimension != dimension)
      throw std::invalid_argument(
          "native line aggregate component dimensions differ");
    // A finite Laurent frame's lower edge is an exact structural bound:
    // powers below it are zero, not unknown.  Aggregation therefore takes
    // the union lower edge while still intersecting complete upper edges.
    epsilon_min = std::min(epsilon_min, value.epsilon.min_power);
    epsilon_max = std::min(epsilon_max, value.epsilon.complete_max);
  }
  if (epsilon_min > epsilon_max)
    throw std::domain_error(
        "native line aggregate components have no common complete epsilon window");

  StoredLineIntegral result;
  result.value.epsilon = {epsilon_min, epsilon_max};
  result.value.dimension = dimension;
  result.value.coefficients.reserve(
      result.value.epsilon.width() * dimension);
  for (std::int64_t raw_power = epsilon_min; raw_power <= epsilon_max;
       ++raw_power) {
    const auto power = static_cast<std::int32_t>(raw_power);
    for (std::uint32_t component_index = 0;
         component_index < dimension; ++component_index) {
      ComplexBall sum(0);
      for (std::size_t index = 0; index < components.size(); ++index) {
        const auto& value = components[index]->result().value;
        if (power < value.epsilon.min_power) continue;
        const auto& coefficient = value.at(power, component_index);
        sum += signs[index] == 1 ? coefficient : -coefficient;
      }
      result.value.coefficients.push_back(std::move(sum));
    }
  }

  bool all_error_envelopes = true;
  bool all_certified_errors = true;
  bool any_advisory_error = false;
  json::array error_sources;
  for (std::size_t index = 0; index < components.size(); ++index) {
    const auto& error = components[index]->result().value.error;
    if (error.empty() || error.frame.complete_max < epsilon_max) {
      all_error_envelopes = false;
      break;
    }
    all_certified_errors &= error.guarantee == ErrorGuarantee::Certified;
    any_advisory_error |= error.guarantee == ErrorGuarantee::Advisory;
    error_sources.push_back(json::object{
        {"sign", signs[index]},
        {"guarantee", error_guarantee_name(error.guarantee)},
        {"provenance", error.provenance}});
  }
  if (all_error_envelopes) {
    result.value.error.frame = result.value.epsilon;
    result.value.error.guarantee = all_certified_errors
        ? ErrorGuarantee::Certified
        : any_advisory_error ? ErrorGuarantee::Advisory
                             : ErrorGuarantee::None;
    result.value.error.absolute.reserve(result.value.epsilon.width());
    for (std::int64_t raw_power = epsilon_min; raw_power <= epsilon_max;
         ++raw_power) {
      const auto power = static_cast<std::int32_t>(raw_power);
      auto error_sum = Magnitude::zero();
      for (const auto& component : components) {
        const auto& error = component->result().value.error;
        if (power < error.frame.min_power) continue;
        error_sum += error.absolute.at(static_cast<std::size_t>(
            power - error.frame.min_power));
      }
      result.value.error.absolute.push_back(std::move(error_sum));
    }
    result.value.error.provenance = json::serialize(json::object{
        {"schema", "diffexp2-native-line-error-sum-v1"},
        {"components", std::move(error_sources)}});
  }

  bool all_full_local = true;
  auto& diagnostics = result.diagnostics;
  diagnostics.detail = detail;
  diagnostics.tail_certificate_requested = true;
  diagnostics.tail_certificate_status = "aggregate-certified";
  for (const auto& component : components) {
    const auto& input = component->result();
    all_full_local &=
        input.scope == LineIntegrationScope::FullLocalWithCertifiedTail;
    const auto& source = input.diagnostics;
    diagnostics.input_monomial_cells = checked_diagnostic_sum(
        diagnostics.input_monomial_cells, source.input_monomial_cells,
        "aggregate input monomial count");
    diagnostics.grouped_monomials = checked_diagnostic_sum(
        diagnostics.grouped_monomials, source.grouped_monomials,
        "aggregate grouped monomial count");
    diagnostics.zero_groups_skipped = checked_diagnostic_sum(
        diagnostics.zero_groups_skipped, source.zero_groups_skipped,
        "aggregate skipped-zero count");
    diagnostics.cancelled_divergent_groups = checked_diagnostic_sum(
        diagnostics.cancelled_divergent_groups,
        source.cancelled_divergent_groups,
        "aggregate cancelled-divergence count");
    diagnostics.bounded_cancelled_divergent_coefficients =
        checked_diagnostic_sum(
            diagnostics.bounded_cancelled_divergent_coefficients,
            source.bounded_cancelled_divergent_coefficients,
            "aggregate bounded-cancelled-divergence count");
    if (source.divergent_cancellation_mode == "bounded-relative-acb") {
      if (diagnostics.divergent_cancellation_mode == "exact-singleton") {
        diagnostics.divergent_cancellation_mode =
            source.divergent_cancellation_mode;
        diagnostics.divergent_relative_tolerance =
            source.divergent_relative_tolerance;
        diagnostics.divergent_cancellation_provenance =
            source.divergent_cancellation_provenance;
      } else if (diagnostics.divergent_relative_tolerance !=
                     source.divergent_relative_tolerance ||
                 diagnostics.divergent_cancellation_provenance !=
                     source.divergent_cancellation_provenance) {
        throw std::invalid_argument(
            "native line aggregate mixes divergent-cancellation policies");
      }
    }
    diagnostics.primitive_evaluations = checked_diagnostic_sum(
        diagnostics.primitive_evaluations, source.primitive_evaluations,
        "aggregate primitive evaluation count");
    diagnostics.primitive_component_applications = checked_diagnostic_sum(
        diagnostics.primitive_component_applications,
        source.primitive_component_applications,
        "aggregate primitive application count");
    diagnostics.primitive_component_reuses = checked_diagnostic_sum(
        diagnostics.primitive_component_reuses,
        source.primitive_component_reuses,
        "aggregate primitive reuse count");
    diagnostics.has_center_endpoint |= source.has_center_endpoint;
    diagnostics.tail_certificate_requested &=
        source.tail_certificate_requested;
  }
  if (all_full_local && result.value.error.guarantee ==
                            ErrorGuarantee::Certified) {
    result.scope = LineIntegrationScope::FullLocalWithCertifiedTail;
  } else {
    result.scope = LineIntegrationScope::StoredTruncation;
    diagnostics.tail_certificate_status = all_full_local
        ? "aggregate-missing-certified-error-envelope"
        : "aggregate-has-stored-truncation-component";
  }
  result.imaginary_sign = std::nullopt;
  return result;
}

// Incremental equivalent of aggregate_retained_lines for the all-positive
// physical-tile sum used by transport observables.  It retains only the
// accumulated epsilon vector and the small per-tile error envelopes; the
// coefficient-bearing projected locals and tile line objects can therefore
// die immediately after each tile.
class StreamingStoredLineAccumulator {
 public:
  void add(const StoredLineIntegral& input) {
    if (!initialized_) {
      result_ = input;
      result_.imaginary_sign = std::nullopt;
      result_.diagnostics.tail_witness_radius_exact.clear();
      initialized_ = true;
    } else {
      if (input.value.dimension != result_.value.dimension)
        throw std::invalid_argument(
            "streaming native line aggregate component dimensions differ");
      const auto minimum = std::min(
          result_.value.epsilon.min_power, input.value.epsilon.min_power);
      const auto maximum = std::min(
          result_.value.epsilon.complete_max,
          input.value.epsilon.complete_max);
      if (minimum > maximum)
        throw std::domain_error(
            "streaming native line aggregate has no common complete epsilon window");
      std::vector<ComplexBall> coefficients;
      coefficients.reserve(EpsilonWindow{minimum, maximum}.width() *
                           result_.value.dimension);
      for (std::int64_t raw_power = minimum; raw_power <= maximum;
           ++raw_power) {
        const auto power = static_cast<std::int32_t>(raw_power);
        for (std::uint32_t component = 0;
             component < result_.value.dimension; ++component) {
          ComplexBall value(0);
          if (power >= result_.value.epsilon.min_power)
            value += result_.value.at(power, component);
          if (power >= input.value.epsilon.min_power)
            value += input.value.at(power, component);
          coefficients.push_back(std::move(value));
        }
      }
      result_.value.epsilon = {minimum, maximum};
      result_.value.coefficients = std::move(coefficients);
      accumulate_diagnostics(input.diagnostics);
    }
    errors_.push_back(input.value.error);
    all_full_local_ &=
        input.scope == LineIntegrationScope::FullLocalWithCertifiedTail;
  }

  StoredLineIntegral finish(const std::string& detail) {
    if (!initialized_)
      throw std::invalid_argument(
          "streaming native line aggregate has no tile values");
    finalize_errors();
    auto& diagnostics = result_.diagnostics;
    diagnostics.detail = detail;
    diagnostics.tail_certificate_status = "aggregate-certified";
    if (all_full_local_ && result_.value.error.guarantee ==
                               ErrorGuarantee::Certified) {
      result_.scope = LineIntegrationScope::FullLocalWithCertifiedTail;
    } else {
      result_.scope = LineIntegrationScope::StoredTruncation;
      diagnostics.tail_certificate_status = all_full_local_
          ? "aggregate-missing-certified-error-envelope"
          : "aggregate-has-stored-truncation-component";
    }
    result_.imaginary_sign = std::nullopt;
    return std::move(result_);
  }

 private:
  void accumulate_diagnostics(
      const StoredLineIntegrationDiagnostics& source) {
    auto& diagnostics = result_.diagnostics;
    diagnostics.input_monomial_cells = checked_diagnostic_sum(
        diagnostics.input_monomial_cells, source.input_monomial_cells,
        "streaming aggregate input monomial count");
    diagnostics.grouped_monomials = checked_diagnostic_sum(
        diagnostics.grouped_monomials, source.grouped_monomials,
        "streaming aggregate grouped monomial count");
    diagnostics.zero_groups_skipped = checked_diagnostic_sum(
        diagnostics.zero_groups_skipped, source.zero_groups_skipped,
        "streaming aggregate skipped-zero count");
    diagnostics.cancelled_divergent_groups = checked_diagnostic_sum(
        diagnostics.cancelled_divergent_groups,
        source.cancelled_divergent_groups,
        "streaming aggregate cancelled-divergence count");
    diagnostics.bounded_cancelled_divergent_coefficients =
        checked_diagnostic_sum(
            diagnostics.bounded_cancelled_divergent_coefficients,
            source.bounded_cancelled_divergent_coefficients,
            "streaming aggregate bounded-cancelled-divergence count");
    if (source.divergent_cancellation_mode == "bounded-relative-acb") {
      if (diagnostics.divergent_cancellation_mode == "exact-singleton") {
        diagnostics.divergent_cancellation_mode =
            source.divergent_cancellation_mode;
        diagnostics.divergent_relative_tolerance =
            source.divergent_relative_tolerance;
        diagnostics.divergent_cancellation_provenance =
            source.divergent_cancellation_provenance;
      } else if (diagnostics.divergent_relative_tolerance !=
                     source.divergent_relative_tolerance ||
                 diagnostics.divergent_cancellation_provenance !=
                     source.divergent_cancellation_provenance) {
        throw std::invalid_argument(
            "streaming native line aggregate mixes divergent-cancellation policies");
      }
    }
    diagnostics.primitive_evaluations = checked_diagnostic_sum(
        diagnostics.primitive_evaluations, source.primitive_evaluations,
        "streaming aggregate primitive evaluation count");
    diagnostics.primitive_component_applications = checked_diagnostic_sum(
        diagnostics.primitive_component_applications,
        source.primitive_component_applications,
        "streaming aggregate primitive application count");
    diagnostics.primitive_component_reuses = checked_diagnostic_sum(
        diagnostics.primitive_component_reuses,
        source.primitive_component_reuses,
        "streaming aggregate primitive reuse count");
    diagnostics.has_center_endpoint |= source.has_center_endpoint;
    diagnostics.tail_certificate_requested &=
        source.tail_certificate_requested;
  }

  void finalize_errors() {
    const auto epsilon = result_.value.epsilon;
    bool all_error_envelopes = true;
    bool all_certified_errors = true;
    bool any_advisory_error = false;
    json::array error_sources;
    error_sources.reserve(errors_.size());
    for (const auto& error : errors_) {
      if (error.empty() || error.frame.complete_max < epsilon.complete_max) {
        all_error_envelopes = false;
        break;
      }
      all_certified_errors &=
          error.guarantee == ErrorGuarantee::Certified;
      any_advisory_error |=
          error.guarantee == ErrorGuarantee::Advisory;
      error_sources.push_back(json::object{
          {"sign", 1},
          {"guarantee", error_guarantee_name(error.guarantee)},
          {"provenance", error.provenance}});
    }
    result_.value.error = ErrorEnvelope{};
    if (!all_error_envelopes) return;
    auto& output = result_.value.error;
    output.frame = epsilon;
    output.guarantee = all_certified_errors
        ? ErrorGuarantee::Certified
        : any_advisory_error ? ErrorGuarantee::Advisory
                             : ErrorGuarantee::None;
    output.absolute.reserve(epsilon.width());
    for (std::int64_t raw_power = epsilon.min_power;
         raw_power <= epsilon.complete_max; ++raw_power) {
      const auto power = static_cast<std::int32_t>(raw_power);
      auto sum = Magnitude::zero();
      for (const auto& error : errors_) {
        if (power < error.frame.min_power) continue;
        sum += error.absolute.at(static_cast<std::size_t>(
            power - error.frame.min_power));
      }
      output.absolute.push_back(std::move(sum));
    }
    output.provenance = json::serialize(json::object{
        {"schema", "diffexp2-native-line-error-sum-v1"},
        {"components", std::move(error_sources)}});
  }

  bool initialized_ = false;
  bool all_full_local_ = true;
  StoredLineIntegral result_;
  std::vector<ErrorEnvelope> errors_;
};

std::vector<std::shared_ptr<StoredLocalBase>> unique_line_local_owners(
    const std::vector<std::shared_ptr<StoredLocalBase>>& owners) {
  std::set<std::string> seen;
  std::vector<std::shared_ptr<StoredLocalBase>> unique;
  unique.reserve(owners.size());
  for (const auto& owner : owners) {
    if (!owner)
      throw std::logic_error("native line aggregate lost a local owner");
    if (seen.insert(owner->handle()).second) unique.push_back(owner);
  }
  return unique;
}

struct RetainedLocalFrameContract {
  EpsilonWindow epsilon;
  std::int32_t top_valid = kCompleteInfinity;
  std::uint32_t dimension = 0;
  std::uint32_t taylor_complete_max = 0;
};

struct WholeArmEpsilonContract {
  EpsilonWindow work;
  std::int32_t public_required_complete_max = 0;
  std::int32_t match_required_complete_max = 0;
};

struct ObservableEpsilonContract {
  EpsilonWindow requested;
  std::int32_t required_complete_max = 0;
};

ObservableEpsilonContract parse_observable_epsilon_contract(
    const json::value& raw, const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(object, {"min", "max", "required_complete_max"},
                     label);
  ObservableEpsilonContract result{
      {as_i32(object.at("min"), "observable epsilon minimum"),
       as_i32(object.at("max"), "observable epsilon maximum")},
      as_i32(object.at("required_complete_max"),
             "observable required epsilon maximum")};
  (void)result.requested.width();
  if (result.required_complete_max < result.requested.min_power ||
      result.required_complete_max > result.requested.complete_max)
    throw std::invalid_argument(
        "observable required epsilon maximum must lie in its exact output window");
  return result;
}

EndpointLimitResult restrict_endpoint_result_epsilon(
    EndpointLimitResult result,
    const ObservableEpsilonContract& contract,
    const char* label) {
  const auto available = endpoint_value_window(result);
  const auto complete_max = std::min(
      contract.requested.complete_max, available.complete_max);
  if (contract.requested.min_power > complete_max ||
      contract.required_complete_max > complete_max)
    throw std::domain_error(
        std::string(label) +
        " does not cover its required exact epsilon output contract");
  std::vector<EpsilonFrame<ComplexBall>> restricted;
  restricted.reserve(result.values.size());
  for (const auto& component : result.values) {
    std::vector<ComplexBall> coefficients;
    coefficients.reserve(static_cast<std::size_t>(
        static_cast<std::int64_t>(complete_max) -
        contract.requested.min_power + 1));
    for (std::int64_t raw_power = contract.requested.min_power;
         raw_power <= complete_max; ++raw_power)
      coefficients.push_back(component.coefficient(
          static_cast<std::int32_t>(raw_power)));
    restricted.emplace_back(
        EpsilonWindow{contract.requested.min_power, complete_max},
        std::move(coefficients));
  }
  result.values = std::move(restricted);
  return result;
}

EndpointLimitResult endpoint_result_from_retained_evaluation(
    const LocalEvaluation& evaluation,
    std::optional<std::int32_t> derived_rim) {
  if (evaluation.value.dimension != 1)
    throw std::logic_error(
        "transport endpoint row evaluation did not remain scalar");
  if (!evaluation.value.error.empty())
    throw std::domain_error(
        "transport endpoint regular-point evaluation produced an error envelope whose endpoint propagation is not implemented");
  std::vector<ComplexBall> coefficients;
  coefficients.reserve(evaluation.value.epsilon.width());
  for (std::int64_t raw_power = evaluation.value.epsilon.min_power;
       raw_power <= evaluation.value.epsilon.complete_max; ++raw_power)
    coefficients.push_back(evaluation.value.at(
        static_cast<std::int32_t>(raw_power), 0));
  EndpointLimitResult result;
  result.values.emplace_back(evaluation.value.epsilon,
                             std::move(coefficients));
  result.imaginary_sign = derived_rim.value_or(1);
  return result;
}

json::object transport_endpoint_provenance(
    const std::string& checkpoint_identity, const json::object& source,
    const json::object& analytic_metadata) {
  return json::object{
      {"schema",
       "diffexp2-retained-native-transport-endpoint-result-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"source", source},
      {"cancellation", json::object{
           {"mode", "exact-or-acb-singleton"}}},
      {"analytic_metadata", analytic_metadata}};
}

json::object transport_centered_projection_record(
    PreparedScalarRowEndpointMode mode,
    std::uint32_t source_taylor_complete_max,
    std::uint32_t projected_taylor_complete_max,
    std::int32_t projected_top_valid,
    const std::string& fallback_reason) {
  json::object record{
      {"schema",
       "diffexp2-centered-scalar-row-endpoint-projection-v1"},
      {"mode", prepared_scalar_row_endpoint_mode_name(mode)},
      {"source_taylor_complete_max", source_taylor_complete_max},
      {"projected_taylor_complete_max", projected_taylor_complete_max},
      {"projected_top_valid", encode_validity(projected_top_valid)}};
  record["fallback_reason"] = fallback_reason.empty()
      ? json::value(nullptr) : json::value(fallback_reason);
  return record;
}

template <typename Scalar>
json::object transport_centered_projection_record(
    const PreparedScalarRowEndpointResult<Scalar>& result) {
  return transport_centered_projection_record(
      result.mode, result.source_taylor_complete_max,
      result.projected_taylor_complete_max,
      result.projected_top_valid, result.fallback_reason);
}

json::object transport_centered_projection_record(
    const PreparedScalarRowEndpointPlan& plan) {
  return transport_centered_projection_record(
      plan.mode, plan.source_taylor_complete_max,
      plan.projected_taylor_complete_max,
      plan.projected_top_valid, plan.fallback_reason);
}

struct CenteredTransportEndpointRow {
  EndpointLimitResult endpoint;
  json::object analytic_metadata;
  json::object projection_record;
};

template <typename Scalar>
CenteredTransportEndpointRow build_centered_transport_endpoint_row(
    const json::object& row,
    const std::shared_ptr<StoredLocal<Scalar>>& source,
    const EndpointLimitOptions& options,
    const std::string& projected_checkpoint_identity) {
  if (!source)
    throw std::logic_error(
        "centered transport endpoint lost its typed final local");
  if (!source->solution().error.empty())
    throw std::domain_error(
        "native rational-row application requires explicit source error-envelope propagation");
  const auto matrix = parse_prepared_rational_row<Scalar>(
      row, source->solution());
  auto computation = centered_prepared_scalar_row_endpoint_limit(
      matrix, source->solution(), source->top_valid(), options,
      projected_checkpoint_identity);
  auto analytic_metadata = checkpoint_local_analytic_metadata_record(
      computation.projected);
  auto projection_record = transport_centered_projection_record(
      computation);
  return {std::move(computation.endpoint),
          std::move(analytic_metadata), std::move(projection_record)};
}

std::shared_ptr<StoredEndpointResult> build_transport_endpoint_row(
    const std::string& endpoint_handle,
    const std::string& checkpoint_identity,
    const std::string& observable_identity,
    const json::object& row,
    const ObservableEpsilonContract& epsilon_contract,
    const json::object& epsilon_record,
    const std::string& projected_handle,
    const std::string& projected_checkpoint_identity,
    const std::string& domain, slong precision_bits,
    const std::shared_ptr<StoredTransportArmState>& state,
    const ResolvedTransportEndpointBinding& binding) {
  if (endpoint_handle.empty() || checkpoint_identity.empty() ||
      observable_identity.empty() || projected_handle.empty() ||
      projected_checkpoint_identity.empty() || !state)
    throw std::invalid_argument(
        "transport endpoint row lost an identity or retained state owner");
  const auto& source = state->final_local();
  const auto source_summary = source->summary();
  const auto started = std::chrono::steady_clock::now();
  validate_prepared_rational_row_structure(
      row, as_u32(source_summary.at("dimension"),
                  "transport endpoint source dimension"),
      "transport endpoint prepared rational row");
  EndpointLimitResult result;
  json::object analytic_metadata;
  std::optional<json::object> projection_record;
  if (binding.centered) {
    EndpointLimitOptions options;
    options.approach_direction = binding.approach_direction;
    options.imaginary_sign = binding.rim;
    options.allow_certified_numeric_cancellation = true;
    if (domain == "rational") {
      const auto typed =
          std::dynamic_pointer_cast<StoredLocal<Rational>>(source);
      if (!typed)
        throw std::logic_error(
            "transport endpoint Rational source changed coefficient domain");
      auto centered = build_centered_transport_endpoint_row(
          row, typed, options, projected_checkpoint_identity);
      result = std::move(centered.endpoint);
      analytic_metadata = std::move(centered.analytic_metadata);
      projection_record = std::move(centered.projection_record);
    } else if (domain == "acb") {
      const auto typed =
          std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(source);
      if (!typed)
        throw std::logic_error(
            "transport endpoint Acb source changed coefficient domain");
      auto centered = build_centered_transport_endpoint_row(
          row, typed, options, projected_checkpoint_identity);
      result = std::move(centered.endpoint);
      analytic_metadata = std::move(centered.analytic_metadata);
      projection_record = std::move(centered.projection_record);
    } else {
      throw std::invalid_argument(
          "transport endpoint row requires one numeric coefficient domain");
    }
  } else {
    json::object row_request{
        {"row", row},
        {"source_checkpoint_identity", source->checkpoint_identity()},
        {"checkpoint_identity", projected_checkpoint_identity}};
    std::shared_ptr<StoredLocalBase> projected;
    if (domain == "rational") {
      const auto typed =
          std::dynamic_pointer_cast<StoredLocal<Rational>>(source);
      if (!typed)
        throw std::logic_error(
            "transport endpoint Rational source changed coefficient domain");
      projected = build_rational_row_local<Rational>(
          projected_handle, row_request, precision_bits, typed, source,
          false);
    } else if (domain == "acb") {
      const auto typed =
          std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(source);
      if (!typed)
        throw std::logic_error(
            "transport endpoint Acb source changed coefficient domain");
      projected = build_rational_row_local<ComplexBall>(
          projected_handle, row_request, precision_bits, typed, source,
          false);
    } else {
      throw std::invalid_argument(
          "transport endpoint row requires one numeric coefficient domain");
    }
    const auto point = RealEvaluationPoint::rational(
        binding.local_end.str());
    const auto evaluation = projected->evaluate_retained_point(
        point, binding.rim);
    result = endpoint_result_from_retained_evaluation(
        evaluation, binding.rim);
    analytic_metadata = projected->exact_analytic_metadata();
  }
  result = restrict_endpoint_result_epsilon(
      std::move(result), epsilon_contract,
      "transport endpoint row result");
  if (result.values.size() != 1)
    throw std::logic_error(
        "transport endpoint row result did not remain scalar");
  auto endpoint_source = binding.source;
  endpoint_source["observable"] = json::object{
      {"identity", observable_identity},
      {"checkpoint_identity", checkpoint_identity}};
  json::object row_record{
      {"exact_identity", required_string(row, "exact_identity")},
      {"prepared_row", row}};
  if (projection_record.has_value())
    row_record["projection"] = std::move(*projection_record);
  endpoint_source["row"] = std::move(row_record);
  endpoint_source["output_epsilon_contract"] = epsilon_record;
  const auto provenance = transport_endpoint_provenance(
      checkpoint_identity, endpoint_source, analytic_metadata);
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredEndpointResult>(
      endpoint_handle, checkpoint_identity, provenance_identity,
      source->handle(), source->source_chart(),
      source->source_operator_identity(), source->checkpoint_identity(),
      source->scalar_domain(), binding.approach_direction, std::nullopt,
      "exact-or-acb-singleton", std::move(analytic_metadata),
      std::move(result), elapsed_ms, std::move(endpoint_source),
      binding.rim, nullptr, nullptr, state);
}

WholeArmEpsilonContract parse_whole_arm_epsilon_contract(
    const json::value& raw, const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(
      object,
      {"min", "max", "required_complete_max",
       "match_required_complete_max"},
      label);
  WholeArmEpsilonContract result{
      {as_i32(object.at("min"), "whole-arm epsilon minimum"),
       as_i32(object.at("max"), "whole-arm epsilon maximum")},
      as_i32(object.at("required_complete_max"),
             "whole-arm projected/public required epsilon maximum"),
      as_i32(object.at("match_required_complete_max"),
             "whole-arm source/match required epsilon maximum")};
  (void)result.work.width();
  if (result.public_required_complete_max < result.work.min_power ||
      result.match_required_complete_max <
          result.public_required_complete_max ||
      result.match_required_complete_max > result.work.complete_max)
    throw std::invalid_argument(
        "whole-arm public and match required epsilon maxima must lie in the work window with match_required_complete_max >= required_complete_max");
  return result;
}

RetainedLocalFrameContract retained_local_frame_contract(
    const std::shared_ptr<StoredLocalBase>& local) {
  if (!local)
    throw std::invalid_argument(
        "native arm frame intersection received a null local");
  const auto summary = local->summary();
  RetainedLocalFrameContract result{
      {as_i32(summary.at("epsilon_min"), "retained epsilon minimum"),
       as_i32(summary.at("epsilon_max"), "retained epsilon maximum")},
      parse_validity(summary.at("top_valid")),
      as_u32(summary.at("dimension"), "retained local dimension"),
      as_u32(summary.at("taylor_complete_max"),
             "retained local Taylor maximum")};
  (void)result.epsilon.width();
  return result;
}

EpsilonWindow live_match_epsilon_intersection(
    EpsilonWindow requested, std::int32_t required_complete_max,
    const std::shared_ptr<StoredLocalBase>& incoming,
    const std::vector<std::shared_ptr<StoredLocalBase>>& basis) {
  if (basis.empty())
    throw std::invalid_argument(
        "native whole-arm match basis cannot be empty");
  const auto incoming_frame = retained_local_frame_contract(incoming);
  auto union_minimum = incoming_frame.epsilon.min_power;
  auto complete_max = requested.complete_max;
  const auto dimension = incoming_frame.dimension;
  if (basis.size() != dimension)
    throw std::invalid_argument(
        "native whole-arm match requires one receiving column per component");
  const auto admit = [&](const RetainedLocalFrameContract& frame) {
    if (frame.dimension != dimension)
      throw std::invalid_argument(
          "native whole-arm matching local dimensions differ");
    // A finite Laurent lower edge is structural: coefficients below it are
    // exact zero.  Clip the caller's lower edge only to the union of actual
    // local frames, then let each matcher zero-pad locals which begin later.
    // Complete upper edges still intersect.
    union_minimum = std::min(union_minimum, frame.epsilon.min_power);
    complete_max = std::min(complete_max, frame.epsilon.complete_max);
    if (frame.top_valid != kCompleteInfinity)
      complete_max = std::min(complete_max, frame.top_valid);
  };
  admit(incoming_frame);
  for (const auto& column : basis)
    admit(retained_local_frame_contract(column));
  const auto minimum = std::max(requested.min_power, union_minimum);
  if (minimum > complete_max)
    throw std::domain_error(
        "native whole-arm match has no common complete epsilon window");
  if (required_complete_max < minimum ||
      required_complete_max > complete_max)
    throw std::domain_error(
        "native whole-arm live match intersection does not cover the globally required complete epsilon maximum");
  return {minimum, complete_max};
}

EpsilonWindow live_line_epsilon_intersection(
    EpsilonWindow requested, std::int32_t required_complete_max,
    const std::shared_ptr<StoredLocalBase>& local,
    std::int32_t primitive_halo) {
  if (primitive_halo < 0)
    throw std::invalid_argument(
        "native whole-arm primitive halo must be nonnegative");
  const auto frame = retained_local_frame_contract(local);
  // A regulated centre-endpoint primitive can create Laurent coefficients
  // below the integrand frame (at most `primitive_halo` rows).  Preserve
  // those rows here just as we already reserve the corresponding upper
  // halo; otherwise the orchestration layer silently turns genuine poles
  // into structural zeros before the line integrator can emit them.
  auto minimum = std::max(
      requested.min_power,
      local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(frame.epsilon.min_power) -
              primitive_halo,
          "native whole-arm line deliverable minimum"));
  auto complete_max = std::min(
      requested.complete_max,
      local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(frame.epsilon.complete_max) -
              primitive_halo,
          "native whole-arm line deliverable maximum"));
  if (frame.top_valid != kCompleteInfinity)
    complete_max = std::min(
        complete_max,
        local_algebra_detail::checked_i32(
            static_cast<std::int64_t>(frame.top_valid) - primitive_halo,
            "native whole-arm line valid deliverable maximum"));
  if (minimum > complete_max || required_complete_max > complete_max)
    throw std::domain_error(
        "native whole-arm integrand row does not cover the globally required complete epsilon maximum");
  return {minimum, complete_max};
}

json::object native_unit_saturation_request(
    const std::shared_ptr<StoredTilePlan>& plan, const std::string& arm,
    std::size_t match_index, bool compact_plan_reference = false) {
  if (!plan)
    throw std::invalid_argument(
        "native unit-saturation request requires its retained tile plan");
  if (compact_plan_reference)
    return json::object{
        {"schema", kNativeUnitSaturationCompactRequestSchema},
        {"tile_plan", plan->handle()},
        {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
        {"arm", arm}, {"match", match_index}};
  return json::object{
        {"schema", kNativeUnitSaturationRequestSchema},
        {"tile_plan", plan->handle()},
        {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
        {"tile_plan_provenance_identity", plan->provenance_identity()},
        {"arm", arm}, {"match", match_index}};
}

NativeAcbSaturationBinding native_acb_saturation_binding(
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::string& session_configuration_identity,
    const std::string& arm,
    std::size_t match_index, const std::string& match_checkpoint_identity,
    bool compact_plan_reference) {
  if (!plan || session_configuration_identity.empty() ||
      match_checkpoint_identity.empty())
    throw std::invalid_argument(
        "native Acb saturation selection requires its session, retained plan, and match checkpoint");
  const auto& retained = plan->arm(arm);
  if (match_index >= retained.exact.matches.size())
    throw std::invalid_argument(
        "native Acb saturation selection match index is outside its retained arm");
  const auto& exact_match = retained.exact.matches[match_index];
  const auto& receiving = retained.charts.at(exact_match.receiving_chart);
  const auto* scc_owner = std::get_if<
      std::shared_ptr<CompositeSCCChartBase>>(&receiving.owner);
  if (scc_owner == nullptr || *scc_owner == nullptr ||
      !is_supported_acb_singular_scc_column_capability(
          (*scc_owner)->column_execution_capability()))
    return {"native_unit_saturation",
            native_unit_saturation_request(
                plan, arm, match_index, compact_plan_reference)};
  const auto receiving_rim = exact_plan_rim(
      receiving.prescriptions, receiving.geometry.scale);
  json::object request{
          {"schema", compact_plan_reference
              ? kNativeSingularSCCSaturationCompactRequestSchema
              : kNativeSingularSCCSaturationRequestSchema},
          {"session_configuration_identity",
           session_configuration_identity},
          {"tile_plan", plan->handle()},
          {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
          {"arm", arm},
          {"match", match_index},
          {"match_checkpoint_identity", match_checkpoint_identity},
          {"receiving_scc", (*scc_owner)->handle()},
          {"receiving_scc_exact_identity", (*scc_owner)->exact_identity()},
          {"receiving_execution_capability",
           (*scc_owner)->column_execution_capability()},
          {"receiving_basis_point_exact",
           exact_match.receiving_local.str()},
          {"physical_match_point_exact", exact_match.physical.str()},
          {"receiving_rim", optional_plan_rim_json(receiving_rim)}};
  if (!compact_plan_reference)
    request["tile_plan_provenance_identity"] =
        plan->provenance_identity();
  return {"native_singular_scc_saturation", std::move(request)};
}

struct RetainedArmMarchInput {
  std::string name;
  std::vector<std::vector<std::string>> basis_handles;
  std::vector<std::vector<std::shared_ptr<StoredLocalBase>>> basis;
  std::vector<std::string> match_handles;
  std::vector<std::string> local_handles;
};

struct RetainedArmMarchResult {
  std::vector<std::shared_ptr<StoredPlannedMatchHop>> matches;
  std::vector<std::shared_ptr<StoredLocalBase>> tile_sources;
  double elapsed_ms = 0.0;

  const std::shared_ptr<StoredLocalBase>& final_local() const {
    if (tile_sources.empty())
      throw std::logic_error(
          "retained arm march has no tile source");
    return tile_sources.back();
  }
};

struct ConsumingArmMarchResult {
  json::array basis_references;
  json::array match_references;
  json::array release_diagnostics;
  std::vector<std::shared_ptr<StoredLocalBase>> tile_sources;
  double elapsed_ms = 0.0;
};

json::object compact_transport_local_reference(
    const std::shared_ptr<StoredLocalBase>& local) {
  if (!local)
    throw std::invalid_argument(
        "compact transport provenance received a null local");
  return json::object{
      {"local", local->handle()},
      {"chart", local->source_chart()},
      {"source_operator_identity", local->source_operator_identity()},
      {"checkpoint_identity", local->checkpoint_identity()},
      {"coefficient_domain", local->scalar_domain()}};
}

std::string arm_checkpoint_identity(const std::string& root,
                                    const std::string& arm,
                                    const char* kind,
                                    std::size_t one_based_index) {
  return root + ":" + arm + ":" + kind + ":" +
         std::to_string(one_based_index);
}

RetainedArmMarchResult march_retained_arm(
    const std::string& domain, slong precision_bits,
    const std::string& session_configuration_identity,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& anchor,
    const RetainedArmMarchInput& input, EpsilonWindow work_epsilon,
    std::int32_t match_required_complete_max,
    const json::object& refinement, const std::string& checkpoint_root) {
  if ((domain != "rational" && domain != "acb") || !plan || !anchor ||
      session_configuration_identity.empty() || checkpoint_root.empty())
    throw std::invalid_argument(
        "retained arm march requires a numeric domain, session identity, plan, anchor, and checkpoint root");
  const auto& retained = plan->arm(input.name);
  const auto match_count = retained.exact.matches.size();
  if (input.basis_handles.size() != match_count ||
      input.basis.size() != match_count ||
      input.match_handles.size() != match_count ||
      input.local_handles.size() != match_count)
    throw std::invalid_argument(
        "retained arm march input does not reproduce its plan match count");
  require_exact_keys(refinement, {"relative_tolerance", "max_steps"},
                     "retained arm march refinement policy");
  if (required_string(refinement, "relative_tolerance").empty() ||
      as_u32(refinement.at("max_steps"),
             "retained arm march refinement steps") > 32)
    throw std::invalid_argument(
        "retained arm march refinement policy is invalid");
  (void)work_epsilon.width();
  if (domain == "acb") ComplexBall::set_precision(precision_bits);

  const auto started = std::chrono::steady_clock::now();
  RetainedArmMarchResult output;
  output.matches.reserve(match_count);
  output.tile_sources.reserve(retained.exact.tiles.size());
  output.tile_sources.push_back(anchor);
  std::shared_ptr<StoredLocalBase> current = anchor;
  for (std::size_t match_index = 0; match_index < match_count;
       ++match_index) {
    const auto match_checkpoint = arm_checkpoint_identity(
        checkpoint_root, input.name, "match", match_index + 1);
    const auto match_epsilon = live_match_epsilon_intersection(
        work_epsilon, match_required_complete_max, current,
        input.basis[match_index]);
    json::object match_request{
        {"arm", input.name}, {"match", match_index},
        {"epsilon", json::object{
             {"min", match_epsilon.min_power},
             {"max", match_epsilon.complete_max},
             {"required_complete_max", match_required_complete_max}}},
        {"checkpoint_identity", match_checkpoint}};
    if (domain == "acb") {
      auto saturation = native_acb_saturation_binding(
          plan, session_configuration_identity, input.name, match_index,
          match_checkpoint);
      match_request[saturation.request_key] =
          std::move(saturation.request);
      match_request["refinement"] = refinement;
    }
    auto match = build_planned_match_hop(
        input.match_handles[match_index], match_request, domain,
        precision_bits, session_configuration_identity, plan,
        input.basis_handles[match_index], input.basis[match_index],
        current->handle(), current);
    std::shared_ptr<StoredLocalBase> next;
    try {
      next = match->materialize(
          input.local_handles[match_index],
          arm_checkpoint_identity(checkpoint_root, input.name, "local",
                                  match_index + 1),
          precision_bits, match);
    } catch (const std::exception& error) {
      const auto incoming_frame = retained_local_frame_contract(current);
      std::string frames =
          "incoming=" +
          std::to_string(incoming_frame.epsilon.complete_max) +
          "/valid=" + std::to_string(incoming_frame.top_valid) +
          ",basis=[";
      for (std::size_t column = 0;
           column < input.basis[match_index].size(); ++column) {
        if (column != 0) frames += ",";
        const auto frame = retained_local_frame_contract(
            input.basis[match_index][column]);
        frames += std::to_string(frame.epsilon.complete_max) +
                  "/valid=" + std::to_string(frame.top_valid);
      }
      frames += "]";
      throw std::domain_error(
          std::string(error.what()) + "; match_window=[" +
          std::to_string(match_epsilon.min_power) + "," +
          std::to_string(match_epsilon.complete_max) +
          "]; match_required=" +
          std::to_string(match_required_complete_max) + "; " + frames);
    }
    output.matches.push_back(std::move(match));
    output.tile_sources.push_back(next);
    current = std::move(next);
  }
  if (output.tile_sources.size() != retained.exact.tiles.size())
    throw std::logic_error(
        "retained arm march did not produce one source local per tile");
  output.elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return output;
}

ConsumingArmMarchResult march_retained_arm_consuming(
    const std::shared_ptr<SolverSession>& session,
    const std::string& session_configuration_identity,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& anchor,
    RetainedArmMarchInput& input, EpsilonWindow work_epsilon,
    std::int32_t match_required_complete_max,
    const json::object& refinement, const std::string& checkpoint_root,
    std::unordered_map<std::string, std::size_t>& remaining_basis_uses) {
  if (!session || session->domain == "symbolic" || !plan || !anchor)
    throw std::invalid_argument(
        "consuming arm march requires one live numeric session, plan, and anchor");
  const auto& retained = plan->arm(input.name);
  const auto match_count = retained.exact.matches.size();
  if (input.basis_handles.size() != match_count ||
      input.match_handles.size() != match_count ||
      input.local_handles.size() != match_count)
    throw std::invalid_argument(
        "consuming arm march input does not reproduce its exact plan");
  require_exact_keys(refinement, {"relative_tolerance", "max_steps"},
                     "consuming arm refinement policy");
  (void)work_epsilon.width();
  if (session->domain == "acb")
    ComplexBall::set_precision(session->precision_bits);

  const auto started = std::chrono::steady_clock::now();
  ConsumingArmMarchResult output;
  output.basis_references.reserve(match_count);
  output.match_references.reserve(match_count);
  output.tile_sources.reserve(retained.exact.tiles.size());
  output.tile_sources.push_back(anchor);
  std::shared_ptr<StoredLocalBase> current = anchor;

  for (std::size_t match_index = 0; match_index < match_count;
       ++match_index) {
    std::vector<std::shared_ptr<StoredLocalBase>> basis;
    json::array basis_reference;
    const auto& handles = input.basis_handles[match_index];
    std::size_t retained_local_count = 0;
    std::uint64_t retained_coefficient_count = 0;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during consuming arm march");
      basis.reserve(handles.size());
      basis_reference.reserve(handles.size());
      for (const auto& handle : handles) {
        const auto found = session->locals.find(handle);
        if (found == session->locals.end() || !found->second)
          throw std::invalid_argument(
              "unknown or already consumed receiving-basis local: " +
              handle);
        basis.push_back(found->second);
        basis_reference.push_back(
            compact_transport_local_reference(found->second));
      }
    }
    if (basis.empty())
      throw std::invalid_argument(
          "consuming arm receiving basis cannot be empty");

    const auto match_checkpoint = arm_checkpoint_identity(
        checkpoint_root, input.name, "match", match_index + 1);
    const auto match_epsilon = live_match_epsilon_intersection(
        work_epsilon, match_required_complete_max, current, basis);
    json::object match_request{
        {"arm", input.name}, {"match", match_index},
        {"epsilon", json::object{
             {"min", match_epsilon.min_power},
             {"max", match_epsilon.complete_max},
             {"required_complete_max", match_required_complete_max}}},
        {"checkpoint_identity", match_checkpoint}};
    if (session->domain == "acb") {
      auto saturation = native_acb_saturation_binding(
          plan, session_configuration_identity, input.name, match_index,
          match_checkpoint, true);
      match_request[saturation.request_key] = std::move(saturation.request);
      match_request["refinement"] = refinement;
    }
    auto match = build_planned_match_hop(
        input.match_handles[match_index], match_request, session->domain,
        session->precision_bits, session_configuration_identity, plan,
        handles, basis, current->handle(), current, true);
    auto next = match->materialize(
        input.local_handles[match_index],
        arm_checkpoint_identity(checkpoint_root, input.name, "local",
                                match_index + 1),
        session->precision_bits, match);

    // At this point the full match, materialized local, tail model and
    // physical residual owner have all been validated.  Seal an immutable
    // lineage record before dropping the coefficient-bearing match owner.
    auto sealed_lineage = next->seal_plan_match_lineage();
    output.basis_references.push_back(std::move(basis_reference));
    output.match_references.push_back(json::object{
        {"index", match_index},
        {"checkpoint_identity", match->checkpoint_identity()},
        {"provenance_identity", match->provenance_identity()},
        {"planned_hop", match->handoff()},
        {"sealed_local_lineage", sealed_lineage}});

    // Consumption is committed per completed exact hop.  Shared basis
    // handles are erased only at their last use across both arms.
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      for (std::size_t column = 0; column < handles.size(); ++column) {
        const auto& handle = handles[column];
        const auto use = remaining_basis_uses.find(handle);
        if (use == remaining_basis_uses.end() || use->second == 0)
          throw std::logic_error(
              "consuming basis last-use accounting underflow");
        --use->second;
        if (use->second == 0) {
          const auto found = session->locals.find(handle);
          if (found == session->locals.end() ||
              found->second.get() != basis[column].get())
            throw std::logic_error(
                "consuming basis registry owner changed before last use");
          session->locals.erase(found);
        }
      }
      retained_local_count = session->locals.size();
      for (const auto& [handle, local] : session->locals) {
        (void)handle;
        const auto count = static_cast<std::uint64_t>(
            local->stats().coefficient_count);
        if (count > std::numeric_limits<std::uint64_t>::max() -
                        retained_coefficient_count)
          throw std::overflow_error(
              "consuming transport coefficient diagnostic overflow");
        retained_coefficient_count += count;
      }
    }
    basis.clear();
    match.reset();
    output.release_diagnostics.push_back(json::object{
        {"arm", input.name}, {"match", match_index},
        {"consumed_handles", handles.size()},
        {"session_locals_after_hop", retained_local_count},
        {"session_local_coefficient_count_after_hop",
         retained_coefficient_count}});
    output.tile_sources.push_back(next);
    current = std::move(next);
  }
  if (output.tile_sources.size() != retained.exact.tiles.size())
    throw std::logic_error(
        "consuming arm march did not produce one local per exact tile");
  output.elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return output;
}

json::array line_aggregate_source_records(
    const std::vector<std::shared_ptr<StoredLocalBase>>& owners) {
  json::array records;
  records.reserve(owners.size());
  for (const auto& owner : owners)
    records.push_back(json::object{
        {"local", owner->handle()}, {"chart", owner->source_chart()},
        {"source_operator_identity", owner->source_operator_identity()},
        {"checkpoint_identity", owner->checkpoint_identity()},
        {"coefficient_domain", owner->scalar_domain()},
        {"analytic_metadata", owner->exact_analytic_metadata()},
        {"retained_derivation", owner->retained_derivation().has_value()
             ? json::value(*owner->retained_derivation())
             : json::value(nullptr)}});
  return records;
}

std::shared_ptr<StoredLineResult> build_retained_line_aggregate(
    const std::string& handle, const std::string& checkpoint_identity,
    const std::string& arm_name, json::object interval,
    json::object aggregate_record,
    const std::shared_ptr<StoredTilePlan>& plan,
    std::vector<std::shared_ptr<StoredLocalBase>> local_owners,
    const std::vector<std::shared_ptr<StoredLineResult>>& components,
    const std::vector<std::int32_t>& signs, double elapsed_ms) {
  if (checkpoint_identity.empty())
    throw std::invalid_argument(
        "native line aggregate checkpoint identity cannot be empty");
  local_owners = unique_line_local_owners(local_owners);
  auto result = aggregate_retained_lines(
      components, signs, "native retained line aggregate");
  json::array component_records;
  component_records.reserve(components.size());
  for (std::size_t index = 0; index < components.size(); ++index) {
    const auto summary = components[index]->summary();
    component_records.push_back(json::object{
        {"index", index}, {"sign", signs[index]},
        {"checkpoint_identity", components[index]->checkpoint_identity()},
        {"provenance_identity", components[index]->provenance_identity()},
        {"scope", summary.at("scope")},
        {"source", summary.at("source")},
        {"interval", summary.at("interval")},
        {"epsilon", json::object{
             {"min", components[index]->result().value.epsilon.min_power},
             {"max", components[index]->result().value.epsilon.complete_max}}},
        {"error", summary.at("error")}});
  }
  aggregate_record["component_count"] = components.size();
  aggregate_record["components"] = std::move(component_records);
  json::object provenance{
      {"schema", "diffexp2-retained-native-line-aggregate-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"arm", arm_name},
      {"interval", std::move(interval)},
      {"source", json::object{
           {"tile_plan", plan->handle()},
           {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
           {"locals", line_aggregate_source_records(local_owners)}}},
      {"aggregate", std::move(aggregate_record)},
      {"epsilon", json::object{{"min", result.value.epsilon.min_power},
                                {"max", result.value.epsilon.complete_max}}},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  return std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, std::move(result),
      elapsed_ms, plan, std::move(local_owners), std::move(provenance));
}

std::shared_ptr<StoredLineResult> build_compact_transport_observable_line(
    const std::string& handle, const std::string& checkpoint_identity,
    const std::string& arm_name, json::object interval,
    json::object aggregate_record,
    const std::shared_ptr<StoredTransportArmState>& transport_state,
    const std::vector<std::shared_ptr<StoredLineResult>>& tile_lines,
    double elapsed_ms) {
  if (!transport_state || checkpoint_identity.empty() || tile_lines.empty())
    throw std::invalid_argument(
        "compact transport observable requires its state, checkpoint, and tile values");
  auto result = aggregate_retained_lines(
      tile_lines, std::vector<std::int32_t>(tile_lines.size(), 1),
      "compact native transport observable aggregate");
  aggregate_record["tile_count"] = tile_lines.size();
  json::object provenance{
      {"schema", "diffexp2-retained-native-transport-observable-line-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"arm", arm_name},
      {"interval", std::move(interval)},
      {"source", json::object{
           {"tile_plan", transport_state->plan_owner()->handle()},
           {"tile_plan_checkpoint_identity",
            transport_state->plan_owner()->checkpoint_identity()},
           {"transport_state", transport_state->handle()},
           {"transport_state_checkpoint_identity",
            transport_state->checkpoint_identity()}}},
      {"aggregate", std::move(aggregate_record)},
      {"epsilon", json::object{{"min", result.value.epsilon.min_power},
                                {"max", result.value.epsilon.complete_max}}},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  return std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, std::move(result),
      elapsed_ms, transport_state->plan_owner(),
      std::vector<std::shared_ptr<StoredLocalBase>>{}, std::move(provenance),
      transport_state);
}

std::shared_ptr<StoredLineResult>
build_compact_transport_observable_line_from_result(
    const std::string& handle, const std::string& checkpoint_identity,
    const std::string& arm_name, json::object interval,
    json::object aggregate_record,
    const std::shared_ptr<StoredTransportArmState>& transport_state,
    StoredLineIntegral result, std::size_t tile_count,
    double elapsed_ms) {
  if (!transport_state || checkpoint_identity.empty() || tile_count == 0)
    throw std::invalid_argument(
        "streaming compact transport observable requires its state, checkpoint, and tile count");
  aggregate_record["tile_count"] = tile_count;
  json::object provenance{
      {"schema", "diffexp2-retained-native-transport-observable-line-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"arm", arm_name},
      {"interval", std::move(interval)},
      {"source", json::object{
           {"tile_plan", transport_state->plan_owner()->handle()},
           {"tile_plan_checkpoint_identity",
            transport_state->plan_owner()->checkpoint_identity()},
           {"transport_state", transport_state->handle()},
           {"transport_state_checkpoint_identity",
            transport_state->checkpoint_identity()}}},
      {"aggregate", std::move(aggregate_record)},
      {"epsilon", json::object{{"min", result.value.epsilon.min_power},
                                {"max", result.value.epsilon.complete_max}}},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  return std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, std::move(result),
      elapsed_ms, transport_state->plan_owner(),
      std::vector<std::shared_ptr<StoredLocalBase>>{}, std::move(provenance),
      transport_state);
}

json::object compact_transport_state_reference(
    const std::shared_ptr<StoredTransportArmState>& state) {
  if (!state)
    throw std::invalid_argument(
        "compact transport reference requires a retained state");
  return json::object{
      {"handle", state->handle()},
      {"checkpoint_identity", state->checkpoint_identity()}};
}

json::object compact_transport_pair_source_side(
    const std::shared_ptr<StoredTransportArmState>& state) {
  return json::object{
      {"tile_plan", state->plan_owner()->handle()},
      {"tile_plan_checkpoint_identity",
       state->plan_owner()->checkpoint_identity()},
      {"transport_state", compact_transport_state_reference(state)}};
}

std::shared_ptr<StoredLineResult>
build_compact_transport_pair_observable_line(
    const std::string& handle, const std::string& checkpoint_identity,
    json::object aggregate_record,
    const std::shared_ptr<StoredTransportArmState>& lower_state,
    const std::shared_ptr<StoredTransportArmState>& upper_state,
    const std::shared_ptr<StoredLineResult>& lower_line,
    const std::shared_ptr<StoredLineResult>& upper_line) {
  const auto started = std::chrono::steady_clock::now();
  if (checkpoint_identity.empty() || !lower_line || !upper_line)
    throw std::invalid_argument(
        "compact transport-pair observable requires its checkpoint and both arm values");
  require_transport_pair_compatibility(
      lower_state, upper_state,
      lower_state ? lower_state->anchor_owner()->scalar_domain() : "");
  auto result = aggregate_retained_lines(
      {lower_line, upper_line}, {-1, 1},
      "compact native transport-pair observable aggregate");
  const auto& lower_exact =
      lower_state->plan_owner()->arm("lower").exact;
  const auto& upper_exact =
      upper_state->plan_owner()->arm("upper").exact;
  json::object provenance{
      {"schema",
       "diffexp2-retained-native-transport-pair-observable-line-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"arm", "combined"},
      {"interval", json::object{
           {"from_exact", lower_exact.to.str()},
           {"anchor_exact", lower_exact.from.str()},
           {"to_exact", upper_exact.to.str()}}},
      {"source", json::object{
           {"lower", compact_transport_pair_source_side(lower_state)},
           {"upper", compact_transport_pair_source_side(upper_state)},
           {"common_anchor", json::object{
                {"local", lower_state->anchor_owner()->handle()},
                {"checkpoint_identity",
                 lower_state->anchor_owner()->checkpoint_identity()},
                {"source_operator_identity",
                 lower_state->anchor_owner()->source_operator_identity()},
                {"physical_exact", lower_exact.from.str()}}}}},
      {"aggregate", std::move(aggregate_record)},
      {"epsilon", json::object{
           {"min", result.value.epsilon.min_power},
           {"max", result.value.epsilon.complete_max}}},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, std::move(result),
      elapsed_ms,
      std::vector<std::shared_ptr<StoredTilePlan>>{
          lower_state->plan_owner(), upper_state->plan_owner()},
      std::vector<std::shared_ptr<StoredLocalBase>>{},
      std::move(provenance),
      std::vector<std::shared_ptr<StoredTransportArmState>>{
          lower_state, upper_state});
}

enum class TransportTailPolicy { Stored, Attempt, Require };

TransportTailPolicy parse_transport_tail_policy(const json::value& raw,
                                                 const char* label) {
  if (!raw.is_string())
    throw std::invalid_argument(std::string(label) + " must be a string");
  const std::string value(raw.as_string());
  if (value == "stored") return TransportTailPolicy::Stored;
  if (value == "attempt") return TransportTailPolicy::Attempt;
  if (value == "require") return TransportTailPolicy::Require;
  throw std::invalid_argument(
      std::string(label) + " must be stored, attempt, or require");
}

const char* transport_tail_policy_name(TransportTailPolicy policy) {
  switch (policy) {
    case TransportTailPolicy::Stored:
      return "stored";
    case TransportTailPolicy::Attempt:
      return "attempt";
    case TransportTailPolicy::Require:
      return "require";
  }
  throw std::logic_error("unknown transport tail policy");
}

const char* transport_projection_mode_name(TransportTailPolicy policy) {
  return policy == TransportTailPolicy::Stored
      ? "fused-stored-hash-monomial-stream-v2"
      : "materialized-row-local-v1";
}

void require_transport_projection_mode(
    const json::object& aggregate, TransportTailPolicy policy,
    const char* label) {
  if (required_string(aggregate, "projection_mode") !=
      transport_projection_mode_name(policy))
    throw std::invalid_argument(
        std::string(label) +
        " does not match its stored-tail projection algorithm");
}

json::array compact_prepared_row_entry_facts(
    const json::object& prepared_row) {
  json::array facts;
  const auto& entries = as_array(
      prepared_row.at("entries"), "compact prepared-row entries");
  facts.reserve(entries.size());
  for (const auto& raw_entry : entries) {
    const auto& entry = as_object(
        raw_entry, "compact prepared-row entry");
    const auto& multiplier = as_object(
        entry.at("multiplier"), "compact prepared-row multiplier");
    facts.push_back(json::object{
        {"column", as_u32(entry.at("column"),
                          "compact prepared-row column")},
        {"epsilon_shift",
         as_i32(multiplier.at("epsilon_shift"),
                "compact prepared-row epsilon shift")},
        {"center_pole_order",
         as_u32(multiplier.at("center_pole_order"),
                "compact prepared-row center-pole order")},
        {"exact_identity",
         required_string(multiplier, "exact_identity")}});
  }
  return facts;
}

struct TransportObservableContractionInput {
  std::string identity;
  std::string checkpoint_identity;
  std::string checkpoint_root;
  std::vector<json::object> rows;
  // Pair contraction borrows prepared rows directly from the immutable
  // synchronous request tree.  This avoids copying the very large exact
  // Taylor-kernel JSON through parse and per-arm staging.  Other seams keep
  // their existing owned rows; exactly one representation may be active.
  std::vector<const json::object*> borrowed_rows;
  ObservableEpsilonContract epsilon;
  json::object epsilon_record;
  TransportTailPolicy tail_policy = TransportTailPolicy::Stored;
  std::optional<BoundedDivergentCancellation> divergent_cancellation;
  std::vector<std::string> projected_local_handles;
  std::string aggregate_handle;
  json::object aggregate_record;

  std::size_t row_count() const {
    if (!rows.empty() && !borrowed_rows.empty())
      throw std::logic_error(
          "transport observable cannot mix owned and borrowed rows");
    return borrowed_rows.empty() ? rows.size() : borrowed_rows.size();
  }

  const json::object& row(std::size_t index) const {
    if (!rows.empty() && !borrowed_rows.empty())
      throw std::logic_error(
          "transport observable cannot mix owned and borrowed rows");
    if (borrowed_rows.empty()) return rows.at(index);
    const auto* borrowed = borrowed_rows.at(index);
    if (borrowed == nullptr)
      throw std::logic_error("transport observable borrowed a null row");
    return *borrowed;
  }
};

struct TransportObservableContractionResult {
  std::string identity;
  // Populated only by the legacy whole-arm seam, whose aggregate owns the
  // projected locals explicitly.  Transport-state contraction streams and
  // leaves both vectors empty.
  std::vector<std::shared_ptr<StoredLocalBase>> projected;
  std::vector<std::shared_ptr<StoredLineResult>> tile_lines;
  std::shared_ptr<StoredLineResult> aggregate;
  std::size_t tile_integrations = 0;
  double tile_integration_ms = 0.0;
  double elapsed_ms = 0.0;
};

template <typename Scalar>
StoredLineIntegral integrate_transport_stored_row_tile(
    slong precision_bits,
    const std::shared_ptr<StoredLocal<Scalar>>& source,
    const json::object& prepared_row,
    const RetainedArmPlan& arm, const std::string& arm_name,
    std::size_t tile_index,
    const ObservableEpsilonContract& epsilon_contract,
    const std::optional<BoundedDivergentCancellation>&
        divergent_cancellation) {
  if (!source || tile_index >= arm.exact.tiles.size())
    throw std::invalid_argument(
        "fused transport row integration lost its source or tile");
  if constexpr (std::is_same_v<Scalar, ComplexBall>)
    ComplexBall::set_precision(precision_bits);
  const auto& tile = arm.exact.tiles[tile_index];
  const auto& binding = arm.charts.at(tile.chart);
  source->require_exact_plan_binding(
      binding.geometry, binding.prescriptions,
      "fused transport row source");
  auto matrix = parse_prepared_rational_row<Scalar>(
      prepared_row, source->solution());
  const auto projection_cap = local_algebra_detail::checked_i32(
      static_cast<std::int64_t>(
          epsilon_contract.requested.complete_max) + 1,
      "fused transport primitive halo");

  auto projected_min = source->solution().epsilon.min_power;
  auto projected_complete = source->solution().epsilon.complete_max;
  auto projected_top_valid = source->top_valid();
  bool active = false;
  for (const auto& entry : matrix.entries) {
    if (entry.multiplier.structurally_zero()) continue;
    const auto term_min = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(source->solution().epsilon.min_power) +
            entry.multiplier.epsilon_shift,
        "fused transport projected minimum");
    const auto term_complete = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(
            source->solution().epsilon.complete_max) +
            entry.multiplier.epsilon_shift,
        "fused transport projected complete maximum");
    const auto term_valid = shifted_local_validity(
        source->top_valid(), entry.multiplier.epsilon_shift);
    if (!active) {
      projected_min = term_min;
      projected_complete = term_complete;
      projected_top_valid = term_valid;
      active = true;
    } else {
      projected_min = std::min(projected_min, term_min);
      projected_complete = std::min(projected_complete, term_complete);
      projected_top_valid = std::min(projected_top_valid, term_valid);
    }
  }
  projected_complete = std::min(projected_complete, projection_cap);
  if (projected_top_valid != kCompleteInfinity)
    projected_complete = std::min(projected_complete, projected_top_valid);
  if (projected_complete < projected_min)
    throw std::domain_error(
        "fused transport row has no valid projected epsilon coefficient");

  // Exact integration at a regulated centre endpoint may deepen the
  // Laurent frame by one power of epsilon.  The low-level integrator
  // zero-pads this extra row for tiles which do not generate it, so keeping
  // it here is both safe and necessary for singular endpoint primitives.
  const auto line_min = std::max(
      epsilon_contract.requested.min_power,
      local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(projected_min) - 1,
          "fused transport line primitive minimum"));
  const auto line_complete = std::min(
      epsilon_contract.requested.complete_max,
      local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(projected_complete) - 1,
          "fused transport line primitive maximum"));
  if (line_min > line_complete ||
      epsilon_contract.required_complete_max > line_complete)
    throw std::domain_error(
        "fused transport row does not cover the globally required epsilon maximum");

  const auto rim = exact_plan_rim(
      binding.prescriptions, binding.geometry.scale);
  const bool reverse_local_orientation =
      tile.local_end < tile.local_begin;
  const auto& primitive_begin = reverse_local_orientation
      ? tile.local_end : tile.local_begin;
  const auto& primitive_end = reverse_local_orientation
      ? tile.local_begin : tile.local_end;
  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {line_min, line_complete};
  options.imaginary_sign = rim;
  options.divergent_cancellation = divergent_cancellation;

  StoredLineIntegral result;
  try {
    result = integrate_prepared_scalar_row_stored(
        matrix, source->solution(), projected_complete,
        RealEvaluationPoint::rational(primitive_begin.str()),
        RealEvaluationPoint::rational(primitive_end.str()), options);
  } catch (const NativeIntegrationError& error) {
    std::ostringstream detail;
    detail << error.what() << "; arm=" << arm_name
           << "; tile=" << tile_index
           << "; physical_interval=[" << tile.physical_begin.str()
           << "," << tile.physical_end.str() << "]"
           << "; local_interval=[" << tile.local_begin.str()
           << "," << tile.local_end.str() << "]"
           << "; chart=" << binding.handle
           << "; chart_center=" << binding.geometry.center.str()
           << "; chart_scale=" << binding.geometry.scale.str()
           << "; rim="
           << (rim.has_value() ? std::to_string(*rim) : "none")
           << "; source_local=" << source->handle()
           << "; fused_row_projection=true";
    NativeIntegrationError contextual(error.code, error.id, detail.str());
    contextual.absolute_power = error.absolute_power;
    contextual.log_power = error.log_power;
    contextual.epsilon_power = error.epsilon_power;
    contextual.component = error.component;
    throw contextual;
  }

  const auto oriented_jacobian = reverse_local_orientation
      ? -binding.geometry.scale : binding.geometry.scale;
  const auto jacobian =
      ComplexBall::from_strings(oriented_jacobian.str());
  for (auto& coefficient : result.value.coefficients)
    coefficient *= jacobian;
  if (!result.value.error.empty()) {
    const auto jacobian_upper = Magnitude::upper_abs(jacobian);
    for (auto& bound : result.value.error.absolute)
      bound = bound * jacobian_upper;
    result.value.error.provenance +=
        "; physical_jacobian_exact=" + oriented_jacobian.str();
  }
  return result;
}

class TransportPairObservableStream final {
 public:
  enum class Status { Active, Poisoned, Finishing, Finished, Aborted };

  struct FinishResult {
    std::shared_ptr<StoredLineResult> line;
    std::array<std::size_t, 2> tiles{0, 0};
    std::array<double, 2> arm_integration_ms{0.0, 0.0};
    double elapsed_ms = 0.0;
  };

  TransportPairObservableStream(
      std::string handle, std::string stream_checkpoint_identity,
      std::string line_handle, std::string checkpoint_root,
      std::string domain, slong precision_bits,
      std::array<std::shared_ptr<StoredTransportArmState>, 2> states,
      std::string identity, std::string checkpoint_identity,
      ObservableEpsilonContract epsilon, json::object epsilon_record,
      TransportTailPolicy tail_policy,
      std::optional<BoundedDivergentCancellation> divergent_cancellation)
      : handle_(std::move(handle)),
        stream_checkpoint_identity_(
            std::move(stream_checkpoint_identity)),
        line_handle_(std::move(line_handle)),
        checkpoint_root_(std::move(checkpoint_root)),
        domain_(std::move(domain)), precision_bits_(precision_bits),
        states_(std::move(states)), identity_(std::move(identity)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        epsilon_(epsilon), epsilon_record_(std::move(epsilon_record)),
        tail_policy_(tail_policy),
        divergent_cancellation_(std::move(divergent_cancellation)) {
    if (handle_.empty() || stream_checkpoint_identity_.empty() ||
        line_handle_.empty() || checkpoint_root_.empty() ||
        identity_.empty() || checkpoint_identity_.empty())
      throw std::invalid_argument(
          "transport-pair stream identities cannot be empty");
    if (domain_ != "rational" && domain_ != "acb")
      throw std::invalid_argument(
          "transport-pair stream requires a numeric session domain");
    if (tail_policy_ != TransportTailPolicy::Stored)
      throw std::invalid_argument(
          "transport-pair tile streaming currently requires stored tail policy");
    if (divergent_cancellation_.has_value() && domain_ != "acb")
      throw std::invalid_argument(
          "bounded divergent cancellation is restricted to Acb transport-pair streams");
    require_transport_pair_compatibility(states_[0], states_[1], domain_);
    for (std::size_t side = 0; side < 2; ++side) {
      expected_tiles_[side] = states_[side]->tile_sources().size();
      if (expected_tiles_[side] == 0)
        throw std::invalid_argument(
            "transport-pair stream state has no retained tiles");
      if (epsilon_.required_complete_max >
          states_[side]->public_required_complete_max())
        throw std::invalid_argument(
            "transport-pair stream epsilon target exceeds a state public target");
      states_[side]->require_contraction_counter_capacity(1);
      row_records_[side].reserve(expected_tiles_[side]);
    }
  }

  const std::string& handle() const { return handle_; }
  const std::string& stream_checkpoint_identity() const {
    return stream_checkpoint_identity_;
  }
  const std::string& identity() const { return identity_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& line_handle() const { return line_handle_; }
  const std::string& checkpoint_root() const { return checkpoint_root_; }
  const auto& states() const { return states_; }
  const auto& expected_tiles() const { return expected_tiles_; }
  const auto& next_tiles() const { return next_tiles_; }

  json::object add_tile(std::size_t side, std::size_t tile_index,
                        const json::object& prepared_row) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (status_ == Status::Poisoned)
      throw std::invalid_argument(
          "transport-pair stream is poisoned and must be aborted");
    if (status_ != Status::Active)
      throw std::invalid_argument(
          "transport-pair stream no longer accepts tiles");
    const std::size_t expected_side =
        next_tiles_[0] < expected_tiles_[0] ? 0 : 1;
    if (side > 1 || side != expected_side)
      throw std::invalid_argument(
          "transport-pair stream requires all lower tiles before upper tiles");
    if (next_tiles_[side] >= expected_tiles_[side] ||
        tile_index != next_tiles_[side])
      throw std::invalid_argument(
          "transport-pair stream tile index is missing, duplicated, or out of order");

    try {
      const auto started = std::chrono::steady_clock::now();
      StoredLineIntegral tile;
      std::unique_ptr<AcbPrecisionLease> acb_lease;
      if (domain_ == "acb") {
        acb_lease = std::make_unique<AcbPrecisionLease>(precision_bits_);
        ComplexBall::set_precision(precision_bits_);
        const auto source =
            std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(
                states_[side]->tile_sources().at(tile_index));
        if (!source)
          throw std::logic_error(
              "transport-pair stream Acb source changed coefficient domain");
        tile = integrate_transport_stored_row_tile<ComplexBall>(
            precision_bits_, source, prepared_row,
            states_[side]->plan_owner()->arm(
                side == 0 ? "lower" : "upper"),
            side == 0 ? "lower" : "upper", tile_index, epsilon_,
            divergent_cancellation_);
      } else {
        const auto source =
            std::dynamic_pointer_cast<StoredLocal<Rational>>(
                states_[side]->tile_sources().at(tile_index));
        if (!source)
          throw std::logic_error(
              "transport-pair stream Rational source changed coefficient domain");
        tile = integrate_transport_stored_row_tile<Rational>(
            precision_bits_, source, prepared_row,
            states_[side]->plan_owner()->arm(
                side == 0 ? "lower" : "upper"),
            side == 0 ? "lower" : "upper", tile_index, epsilon_,
            divergent_cancellation_);
      }
      accumulators_[side].add(tile);
      const auto row_identity = required_string(
          prepared_row, "exact_identity");
      row_records_[side].push_back(json::object{
          {"tile", tile_index}, {"exact_identity", row_identity},
          {"entries", compact_prepared_row_entry_facts(prepared_row)}});
      ++next_tiles_[side];
      const auto elapsed_ms = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      arm_integration_ms_[side] += elapsed_ms;
      return json::object{
          {"side", side == 0 ? "lower" : "upper"},
          {"tile", tile_index}, {"row_exact_identity", row_identity},
          {"lower_complete", next_tiles_[0]},
          {"upper_complete", next_tiles_[1]},
          {"next_side", next_tiles_[0] < expected_tiles_[0]
               ? json::value("lower")
               : next_tiles_[1] < expected_tiles_[1]
                   ? json::value("upper") : json::value(nullptr)},
          {"next_tile", next_tiles_[0] < expected_tiles_[0]
               ? json::value(next_tiles_[0])
               : next_tiles_[1] < expected_tiles_[1]
                   ? json::value(next_tiles_[1]) : json::value(nullptr)},
          {"elapsed_ms", elapsed_ms}};
    } catch (...) {
      status_ = Status::Poisoned;
      throw;
    }
  }

  FinishResult finish() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (status_ == Status::Poisoned)
      throw std::invalid_argument(
          "transport-pair stream is poisoned and cannot finish");
    if (status_ != Status::Active)
      throw std::invalid_argument(
          "transport-pair stream is not active");
    if (next_tiles_ != expected_tiles_)
      throw std::invalid_argument(
          "transport-pair stream cannot finish before every ordered tile was added exactly once");
    status_ = Status::Finishing;
    try {
    const auto finish_started = std::chrono::steady_clock::now();
    const auto observable_root = checkpoint_root_ + ":observable:1";
    std::array<std::shared_ptr<StoredLineResult>, 2> arm_lines;
    for (std::size_t side = 0; side < 2; ++side) {
      const std::string arm_name = side == 0 ? "lower" : "upper";
      const auto arm_identity = identity_ + ":" + arm_name;
      const auto arm_checkpoint =
          observable_root + ":" + arm_name + ":scratch";
      json::object arm_record{
          {"kind", "transport-state-observable-arm"},
          {"combination", "sum-physical-tiles"},
          {"request_index", 0},
          {"observable_identity", arm_identity},
          {"observable_checkpoint_identity", arm_checkpoint},
          {"transport_state",
           compact_transport_state_reference(states_[side])},
          {"output_epsilon_contract", epsilon_record_},
          {"tail_policy", transport_tail_policy_name(tail_policy_)},
          {"projection_mode",
           transport_projection_mode_name(tail_policy_)},
          {"divergent_cancellation",
           divergent_cancellation_.has_value()
               ? json::value(encode_bounded_divergent_cancellation(
                     *divergent_cancellation_))
               : json::value(json::object{{"mode", "exact-singleton"}})},
          {"rows", row_records_[side]},
          {"tile_count", expected_tiles_[side]}};
      const auto& exact = states_[side]->plan_owner()->arm(arm_name).exact;
      arm_lines[side] = build_compact_transport_observable_line_from_result(
          "private:" + observable_root + ":" + arm_name + ":aggregate",
          arm_checkpoint, arm_name,
          json::object{{"from_exact", exact.from.str()},
                       {"to_exact", exact.to.str()}},
          std::move(arm_record), states_[side],
          accumulators_[side].finish(
              "compact streamed native transport observable aggregate"),
          expected_tiles_[side], arm_integration_ms_[side]);
    }

    json::object pair_record{
        {"kind", "transport-state-observable-pair"},
        {"combination", "negative-lower-plus-upper"},
        {"no_remarching", true}, {"no_rematching", true},
        {"concurrent_arms", true}, {"request_index", 0},
        {"observable_identity", identity_},
        {"observable_checkpoint_identity", checkpoint_identity_},
        {"output_epsilon_contract", epsilon_record_},
        {"tail_policy", transport_tail_policy_name(tail_policy_)},
        {"projection_mode", transport_projection_mode_name(tail_policy_)},
        {"divergent_cancellation",
         divergent_cancellation_.has_value()
             ? json::value(encode_bounded_divergent_cancellation(
                   *divergent_cancellation_))
             : json::value(json::object{{"mode", "exact-singleton"}})},
        {"lower", json::object{
             {"transport_state",
              compact_transport_state_reference(states_[0])},
             {"rows", row_records_[0]},
             {"tile_count", expected_tiles_[0]}}},
        {"upper", json::object{
             {"transport_state",
              compact_transport_state_reference(states_[1])},
             {"rows", row_records_[1]},
             {"tile_count", expected_tiles_[1]}}}};
    auto line = build_compact_transport_pair_observable_line(
        line_handle_, checkpoint_identity_, std::move(pair_record),
        states_[0], states_[1], arm_lines[0], arm_lines[1]);
    if (line->result().value.epsilon.complete_max <
        epsilon_.required_complete_max)
      throw std::domain_error(
          "streamed transport-pair aggregate does not cover its required epsilon maximum");
    FinishResult output;
    output.line = std::move(line);
    output.tiles = expected_tiles_;
    output.arm_integration_ms = arm_integration_ms_;
    output.elapsed_ms = arm_integration_ms_[0] + arm_integration_ms_[1] +
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - finish_started).count();
    status_ = Status::Finished;
    return output;
    } catch (...) {
      status_ = Status::Poisoned;
      throw;
    }
  }

  void abort() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (status_ == Status::Finished)
      throw std::invalid_argument(
          "finished transport-pair stream cannot be aborted");
    status_ = Status::Aborted;
  }

 private:
  std::string handle_;
  std::string stream_checkpoint_identity_;
  std::string line_handle_;
  std::string checkpoint_root_;
  std::string domain_;
  slong precision_bits_ = 256;
  std::array<std::shared_ptr<StoredTransportArmState>, 2> states_;
  std::string identity_;
  std::string checkpoint_identity_;
  ObservableEpsilonContract epsilon_;
  json::object epsilon_record_;
  TransportTailPolicy tail_policy_ = TransportTailPolicy::Stored;
  std::optional<BoundedDivergentCancellation> divergent_cancellation_;
  std::array<std::size_t, 2> expected_tiles_{0, 0};
  std::array<std::size_t, 2> next_tiles_{0, 0};
  std::array<StreamingStoredLineAccumulator, 2> accumulators_;
  std::array<json::array, 2> row_records_;
  std::array<double, 2> arm_integration_ms_{0.0, 0.0};
  Status status_ = Status::Active;
  std::mutex mutex_;
};

std::vector<TransportObservableContractionResult> contract_transport_arm(
    const std::string& domain, slong precision_bits,
    const std::shared_ptr<StoredTilePlan>& plan, const std::string& arm_name,
    const std::vector<std::shared_ptr<StoredLocalBase>>& tile_sources,
    const std::vector<TransportObservableContractionInput>& observables,
    const std::shared_ptr<StoredTransportArmState>& transport_owner =
        nullptr) {
  if ((domain != "rational" && domain != "acb") || !plan)
    throw std::invalid_argument(
        "native transport contraction requires a numeric domain and retained plan");
  const auto& retained = plan->arm(arm_name);
  if (tile_sources.size() != retained.exact.tiles.size() ||
      tile_sources.empty())
    throw std::invalid_argument(
        "native transport contraction sources do not reproduce the retained arm tiles");
  if (domain == "acb") ComplexBall::set_precision(precision_bits);

  std::vector<TransportObservableContractionResult> results;
  results.reserve(observables.size());
  for (const auto& observable : observables) {
    if (observable.identity.empty() || observable.checkpoint_identity.empty() ||
        observable.checkpoint_root.empty() ||
        observable.aggregate_handle.empty() ||
        observable.row_count() != tile_sources.size() ||
        observable.projected_local_handles.size() != tile_sources.size())
      throw std::invalid_argument(
          "native observable contraction does not reproduce its retained tile topology or identities");
    const auto started = std::chrono::steady_clock::now();
    TransportObservableContractionResult output;
    output.identity = observable.identity;
    StreamingStoredLineAccumulator accumulator;
    for (std::size_t tile = 0; tile < tile_sources.size(); ++tile) {
      const auto& source = tile_sources[tile];
      const auto& prepared_row = observable.row(tile);
      const auto row_identity = required_string(
          prepared_row, "exact_identity");
      (void)row_identity;
      if (transport_owner &&
          observable.tail_policy == TransportTailPolicy::Stored) {
        const auto fused_started = std::chrono::steady_clock::now();
        StoredLineIntegral fused;
        if (domain == "rational") {
          const auto typed =
              std::dynamic_pointer_cast<StoredLocal<Rational>>(source);
          if (!typed)
            throw std::logic_error(
                "fused transport Rational source changed coefficient domain");
          fused = integrate_transport_stored_row_tile<Rational>(
              precision_bits, typed, prepared_row, retained, arm_name, tile,
              observable.epsilon, observable.divergent_cancellation);
        } else {
          const auto typed =
              std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(source);
          if (!typed)
            throw std::logic_error(
                "fused transport Acb source changed coefficient domain");
          fused = integrate_transport_stored_row_tile<ComplexBall>(
              precision_bits, typed, prepared_row, retained, arm_name, tile,
              observable.epsilon, observable.divergent_cancellation);
        }
        output.tile_integration_ms +=
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - fused_started)
                .count();
        accumulator.add(fused);
        ++output.tile_integrations;
        continue;
      }
      json::object row_request{
          {"source_checkpoint_identity", source->checkpoint_identity()},
          {"checkpoint_identity",
           arm_checkpoint_identity(observable.checkpoint_root, arm_name,
                                   "integrand", tile + 1) +
               ":" + row_identity}};
      const auto projection_complete_max =
          local_algebra_detail::checked_i32(
              static_cast<std::int64_t>(
                  observable.epsilon.requested.complete_max) + 1,
              "transport scalar-row primitive halo");
      std::shared_ptr<StoredLocalBase> projected;
      if (domain == "rational") {
        const auto typed =
            std::dynamic_pointer_cast<StoredLocal<Rational>>(source);
        if (!typed)
          throw std::logic_error(
              "transport contraction Rational source changed coefficient domain");
        projected = build_rational_row_local<Rational>(
            observable.projected_local_handles[tile], row_request,
            precision_bits, typed, source,
            observable.tail_policy != TransportTailPolicy::Stored,
            projection_complete_max, &prepared_row);
      } else {
        const auto typed =
            std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(source);
        if (!typed)
          throw std::logic_error(
              "transport contraction Acb source changed coefficient domain");
        projected = build_rational_row_local<ComplexBall>(
            observable.projected_local_handles[tile], row_request,
            precision_bits, typed, source,
            observable.tail_policy != TransportTailPolicy::Stored,
            projection_complete_max, &prepared_row);
      }
      const auto line_epsilon = live_line_epsilon_intersection(
          observable.epsilon.requested,
          observable.epsilon.required_complete_max, projected,
          1 /* exact native primitive bound: min epsilon power >= -1 */);
      json::object line_request{
          {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
          {"source_checkpoint_identity", projected->checkpoint_identity()},
          {"checkpoint_identity",
           arm_checkpoint_identity(observable.checkpoint_root, arm_name,
                                   "tile", tile + 1)},
          {"arm", arm_name},
          {"tile", tile},
          {"epsilon", json::object{{"min", line_epsilon.min_power},
                                    {"max", line_epsilon.complete_max}}}};
      if (observable.tail_policy != TransportTailPolicy::Stored)
        line_request["certify_tail"] = true;
      if (observable.divergent_cancellation.has_value())
        line_request["divergent_cancellation"] =
            encode_bounded_divergent_cancellation(
                *observable.divergent_cancellation);
      auto tile_line = build_planned_line_result(
          "private:" + arm_checkpoint_identity(
                           observable.checkpoint_root, arm_name, "tile",
                           tile + 1),
          line_request, plan, projected);
      output.tile_integration_ms += tile_line->elapsed_ms();
      if (transport_owner) {
        accumulator.add(tile_line->result());
      } else {
        output.projected.push_back(projected);
        output.tile_lines.push_back(tile_line);
      }
      ++output.tile_integrations;
      // The line owns the projected local and the projected local owns its
      // source.  Drop both before advancing so contraction peak memory is
      // independent of the number of tiles.
      if (transport_owner) {
        tile_line.reset();
        projected.reset();
      }
    }
    const auto aggregate_started = std::chrono::steady_clock::now();
    const auto interval = json::object{
        {"from_exact", retained.exact.from.str()},
        {"to_exact", retained.exact.to.str()}};
    if (transport_owner) {
      output.aggregate = build_compact_transport_observable_line_from_result(
          observable.aggregate_handle, observable.checkpoint_identity,
          arm_name, interval, observable.aggregate_record, transport_owner,
          accumulator.finish("compact native transport observable aggregate"),
          output.tile_integrations,
          std::chrono::duration<double, std::milli>(
              std::chrono::steady_clock::now() - aggregate_started)
              .count());
    } else {
      output.aggregate = build_retained_line_aggregate(
          observable.aggregate_handle, observable.checkpoint_identity,
          arm_name, interval, observable.aggregate_record, plan,
          output.projected, output.tile_lines,
          std::vector<std::int32_t>(output.tile_lines.size(), 1),
          std::chrono::duration<double, std::milli>(
              std::chrono::steady_clock::now() - aggregate_started)
              .count());
    }
    if (output.aggregate->result().value.epsilon.complete_max <
        observable.epsilon.required_complete_max)
      throw std::domain_error(
          "native observable aggregate does not cover its required epsilon maximum");
    if (observable.tail_policy == TransportTailPolicy::Require &&
        output.aggregate->result().scope !=
            LineIntegrationScope::FullLocalWithCertifiedTail)
      throw std::domain_error(
          "native observable contraction requires a certified full-local tail but certification was unavailable");
    output.elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    results.push_back(std::move(output));
  }
  return results;
}
