constexpr std::uint32_t kMaxPersistentBatchThreads = 64;

std::string composite_scc_signature(const json::object& root) {
  json::object exact{{"identity", root.at("identity")},
                     {"parent", root.at("parent")},
                     {"blocks", root.at("blocks")},
                     {"couplings", root.at("couplings")}};
  if (const auto* shadow = root.if_contains("rational_shadow_identity"))
    exact["rational_shadow_identity"] = *shadow;
  if (const auto* physical = root.if_contains("physical_ode"))
    exact["physical_ode"] = *physical;
  return json::serialize(canonical_json_value(exact));
}

std::string static_problem_signature(const json::object& problem,
                                     const json::value& analytic,
                                     const SCCCertificate& scc,
                                     const std::string& identity) {
  json::object exact;
  for (const auto* key : {"domain", "symbols", "precision_bits", "d", "fb",
                          "w", "d_lags", "denominators", "nhat_lags",
                          "d0_inverse", "blocks", "epsilon_regular_principal",
                          "spectral_principal", "spectral_source",
                          "assembly", "physical_ode",
                          "chop_digits"}) {
    if (const auto* value = problem.if_contains(key)) exact[key] = *value;
  }
  exact["identity"] = identity;
  exact["analytic"] = analytic;
  exact["scc"] = json::parse(scc.exact_record);
  return json::serialize(exact);
}

json::object run_session_command(const json::object& root);

std::uint64_t scoped_handle_id(const std::string& handle,
                               std::string_view prefix,
                               const char* label) {
  if (!handle.starts_with(prefix) || handle.size() == prefix.size())
    throw std::invalid_argument(std::string(label) +
                                " has an invalid scoped handle");
  std::uint64_t value = 0;
  for (std::size_t index = prefix.size(); index < handle.size(); ++index) {
    const char digit = handle[index];
    if (digit < '0' || digit > '9')
      throw std::invalid_argument(std::string(label) +
                                  " has an invalid scoped handle");
    const auto unsigned_digit = static_cast<std::uint64_t>(digit - '0');
    if (value > (std::numeric_limits<std::uint64_t>::max() -
                 unsigned_digit) / 10)
      throw std::invalid_argument(std::string(label) +
                                  " scoped handle overflows uint64");
    value = value * 10 + unsigned_digit;
  }
  if (value == 0)
    throw std::invalid_argument(std::string(label) +
                                " scoped handle index must be positive");
  return value;
}

std::pair<ExactArmRequest, std::vector<RetainedPlanChartBinding>>
parse_checkpoint_retained_arm(
    const json::value& raw,
    const std::unordered_map<std::string,
                             std::shared_ptr<PreparedChartBase>>& charts,
    const std::unordered_map<std::string,
                             std::shared_ptr<CompositeSCCChartBase>>& sccs,
    const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(object,
      {"from_exact", "to_exact", "direction", "division_order", "charts",
       "matches", "tiles", "planning_charts", "certified_geometry",
       "topology"}, label);
  ExactArmRequest request;
  request.from = parse_exact_path_rational(object.at("from_exact"), label);
  request.to = parse_exact_path_rational(object.at("to_exact"), label);
  request.topology = parse_exact_path_topology(object.at("topology"));
  std::vector<RetainedPlanChartBinding> bindings;
  const auto& raw_charts = as_array(object.at("charts"), label);
  if (raw_charts.empty())
    throw std::invalid_argument(std::string(label) +
                                " has no retained charts");
  bindings.reserve(raw_charts.size());
  request.charts.reserve(raw_charts.size());
  const auto& raw_planning = as_array(
      object.at("planning_charts"), "checkpoint tile planning charts");
  const auto& raw_certified = as_object(
      object.at("certified_geometry"),
      "checkpoint certified algebraic arm geometry");
  const auto& certified_charts = as_array(
      raw_certified.at("charts"),
      "checkpoint certified algebraic arm charts");
  if (raw_planning.size() != raw_charts.size() ||
      certified_charts.size() != raw_charts.size())
    throw std::invalid_argument(std::string(label) +
                                " planning/certified chart count differs from its retained owners");
  for (std::size_t index = 0; index < raw_charts.size(); ++index) {
    const auto& raw_chart = as_object(raw_charts[index], label);
    require_exact_keys(raw_chart,
        {"index", "chart", "identity", "center_exact", "scale_exact",
         "radius_exact", "singular_center", "prescriptions"}, label);
    if (as_u64(raw_chart.at("index"), label) != index)
      throw std::invalid_argument(std::string(label) +
                                  " chart bindings are not in index order");
    const auto handle = required_string(raw_chart, "chart");
    RetainedPlanChartBinding::Owner owner;
    if (handle.starts_with("c:")) {
      const auto found = charts.find(handle);
      if (found == charts.end())
        throw std::invalid_argument(
            std::string(label) +
            " references an absent prepared-chart owner");
      owner = found->second;
    } else if (handle.starts_with("scc:")) {
      const auto found = sccs.find(handle);
      if (found == sccs.end())
        throw std::invalid_argument(
            std::string(label) +
            " references an absent composite-SCC owner");
      owner = found->second;
    } else {
      throw std::invalid_argument(
          std::string(label) +
          " chart binding is neither a prepared-chart nor composite-SCC handle");
    }
    const auto& planning = as_object(
        raw_planning[index], "checkpoint tile planning chart");
    const auto& certified = as_object(
        certified_charts[index], "checkpoint certified algebraic chart");
    if (as_u64(certified.at("index"),
               "checkpoint certified chart index") != index)
      throw std::invalid_argument(std::string(label) +
                                  " certified charts are not in owner order");
    auto binding = bind_plan_chart(owner, request.topology,
                                   planning, certified);
    if (encode_plan_chart(binding, index) != raw_chart)
      throw std::invalid_argument(std::string(label) +
                                  " chart binding differs from its exact retained owner");
    request.charts.push_back(binding.geometry);
    bindings.push_back(std::move(binding));
  }
  // The plan is rederived below.  These fields are still parsed here to make
  // malformed scalar types fail before any planning work is attempted.
  (void)as_i32(object.at("direction"), label);
  (void)as_u32(object.at("division_order"), label);
  (void)as_array(object.at("matches"), label);
  (void)as_array(object.at("tiles"), label);
  return {std::move(request), std::move(bindings)};
}

std::shared_ptr<StoredTilePlan> restore_checkpoint_tile_plan_record(
    const json::value& raw,
    const std::unordered_map<std::string,
                             std::shared_ptr<PreparedChartBase>>& charts,
    const std::unordered_map<std::string,
                             std::shared_ptr<CompositeSCCChartBase>>& sccs) {
  const auto& object = as_object(raw, "checkpoint retained tile plan");
  if (required_string(object, "schema") ==
      kRetainedSingleArmTilePlanCheckpointSchema) {
    require_exact_keys(
        object,
        {"schema", "handle", "checkpoint_identity", "provenance_identity",
         "division_order", "arm_name", "arm", "elapsed_ms",
         "runtime_stats"},
        "checkpoint retained single-arm tile plan");
    const auto handle = required_string(object, "handle");
    const auto checkpoint_identity = required_string(
        object, "checkpoint_identity");
    const auto provenance_identity = required_string(
        object, "provenance_identity");
    const auto arm_name = required_string(object, "arm_name");
    if (handle.empty() || checkpoint_identity.empty() ||
        provenance_identity.empty() ||
        (arm_name != "lower" && arm_name != "upper"))
      throw std::invalid_argument(
          "checkpoint retained single-arm tile plan contains an invalid identity or arm name");
    const auto division_order = as_u32(
        object.at("division_order"),
        "checkpoint single-arm tile division order");
    auto [arm_request, bindings] = parse_checkpoint_retained_arm(
        object.at("arm"), charts, sccs,
        "checkpoint retained single tile arm");
    ExactPathPlanOptions options;
    options.division_order = division_order;
    auto exact = plan_exact_arm(arm_request, options);
    const std::string derived_name =
        exact.direction < 0 ? "lower" : "upper";
    if (derived_name != arm_name)
      throw std::invalid_argument(
          "checkpoint single-arm tile-plan name differs from its exact direction");
    RetainedArmPlan retained{std::move(exact), std::move(bindings)};
    attach_certified_arm_geometry(
        as_object(object.at("arm"),
                  "checkpoint retained single tile arm"), retained);
    if (encode_retained_arm(retained) != object.at("arm"))
      throw std::invalid_argument(
          "checkpoint single-arm tile intervals do not reproduce the exact planner result");
    json::object provenance{
        {"schema", kRetainedSingleArmTilePlanProvenanceSchema},
        {"checkpoint_identity", checkpoint_identity},
        {"division_order", division_order},
        {"arm_name", arm_name},
        {"arm", encode_retained_arm(retained)}};
    if (json::serialize(canonical_json_value(provenance)) !=
        provenance_identity)
      throw std::invalid_argument(
          "checkpoint single-arm tile-plan provenance identity is inconsistent");
    const auto elapsed_ms = checkpoint_nonnegative_double(
        object.at("elapsed_ms"),
        "checkpoint single-arm tile-plan elapsed time");
    const auto& stats = as_object(
        object.at("runtime_stats"),
        "checkpoint single-arm tile-plan runtime stats");
    require_exact_keys(
        stats,
        {"match_interval_queries", "tile_interval_queries",
         "lower_match_advances", "upper_match_advances", "integrations"},
        "checkpoint single-arm tile-plan runtime stats");
    auto plan = std::make_shared<StoredTilePlan>(
        handle, checkpoint_identity, provenance_identity, division_order,
        arm_name, std::move(retained), elapsed_ms);
    plan->restore_runtime_stats(
        as_u64(stats.at("match_interval_queries"),
               "checkpoint single-arm tile match queries"),
        as_u64(stats.at("tile_interval_queries"),
               "checkpoint single-arm tile interval queries"),
        as_u64(stats.at("lower_match_advances"),
               "checkpoint single-arm lower match advances"),
        as_u64(stats.at("upper_match_advances"),
               "checkpoint single-arm upper match advances"),
        as_u64(stats.at("integrations"),
               "checkpoint single-arm tile integrations"));
    return plan;
  }
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "division_order", "lower", "upper", "elapsed_ms",
       "runtime_stats"},
      "checkpoint retained tile plan");
  if (required_string(object, "schema") !=
      "diffexp2-retained-tile-plan-v2")
    throw std::invalid_argument(
        "unsupported retained tile-plan checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty())
    throw std::invalid_argument(
        "checkpoint retained tile plan contains an empty identity");
  const auto division_order = as_u32(
      object.at("division_order"), "checkpoint tile division order");
  auto [lower_request, lower_bindings] = parse_checkpoint_retained_arm(
      object.at("lower"), charts, sccs, "checkpoint lower tile arm");
  auto [upper_request, upper_bindings] = parse_checkpoint_retained_arm(
      object.at("upper"), charts, sccs, "checkpoint upper tile arm");
  if (lower_bindings.front().handle != upper_bindings.front().handle)
    throw std::invalid_argument(
        "checkpoint independent tile arms lost their shared anchor owner");
  ExactPathPlanOptions options;
  options.division_order = division_order;
  auto planned = plan_exact_independent_arms(
      lower_request, upper_request, options);
  RetainedArmPlan lower{std::move(planned.lower),
                        std::move(lower_bindings)};
  RetainedArmPlan upper{std::move(planned.upper),
                        std::move(upper_bindings)};
  attach_certified_arm_geometry(
      as_object(object.at("lower"), "checkpoint lower tile arm"), lower);
  attach_certified_arm_geometry(
      as_object(object.at("upper"), "checkpoint upper tile arm"), upper);
  if (encode_retained_arm(lower) != object.at("lower") ||
      encode_retained_arm(upper) != object.at("upper"))
    throw std::invalid_argument(
        "checkpoint tile intervals do not reproduce the exact planner result");
  json::object provenance{
      {"schema", "diffexp2-retained-exact-independent-arm-tile-plan-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"division_order", division_order},
      {"lower", encode_retained_arm(lower)},
      {"upper", encode_retained_arm(upper)}};
  if (json::serialize(canonical_json_value(provenance)) !=
      provenance_identity)
    throw std::invalid_argument(
        "checkpoint tile-plan provenance identity is inconsistent");
  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint tile-plan elapsed time");
  const auto& stats = as_object(object.at("runtime_stats"),
                                "checkpoint tile-plan runtime stats");
  require_exact_keys(stats,
      {"match_interval_queries", "tile_interval_queries",
       "lower_match_advances", "upper_match_advances", "integrations"},
      "checkpoint tile-plan runtime stats");
  auto plan = std::make_shared<StoredTilePlan>(
      handle, checkpoint_identity, provenance_identity, division_order,
      std::move(lower), std::move(upper), elapsed_ms);
  plan->restore_runtime_stats(
      as_u64(stats.at("match_interval_queries"),
             "checkpoint tile match queries"),
      as_u64(stats.at("tile_interval_queries"),
             "checkpoint tile interval queries"),
      as_u64(stats.at("lower_match_advances"),
             "checkpoint lower match advances"),
      as_u64(stats.at("upper_match_advances"),
             "checkpoint upper match advances"),
      as_u64(stats.at("integrations"),
             "checkpoint tile integrations"));
  return plan;
}

std::shared_ptr<StoredEndpointResult>
restore_checkpoint_planned_endpoint_record(
    const json::value& raw, const std::string& expected_domain,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& local) {
  const auto& object = as_object(
      raw, "checkpoint retained plan-bound endpoint");
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "source", "approach_direction", "derived_rim",
       "cancellation_mode", "analytic_metadata", "result", "elapsed_ms",
       "runtime_stats"},
      "checkpoint retained plan-bound endpoint");
  if (required_string(object, "schema") !=
      "diffexp2-retained-plan-bound-endpoint-result-v1")
    throw std::invalid_argument(
        "unsupported retained plan-bound endpoint checkpoint schema");
  if (!plan || !local)
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint lost a strongly owned plan or local");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty())
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint contains an empty identity");
  if (std::string(local->scalar_domain()) != expected_domain ||
      (expected_domain != "rational" && expected_domain != "acb"))
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint coefficient domain differs from its session/local");

  const auto& source = as_object(
      object.at("source"), "checkpoint plan-bound endpoint source");
  const auto arm_name = required_string(source, "arm");
  const auto binding = resolve_planned_endpoint_binding(
      plan, arm_name, local);
  if (source != binding.source)
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint source differs from its exact plan/local owners");
  const auto approach_direction = as_i32(
      object.at("approach_direction"),
      "checkpoint plan-bound endpoint approach direction");
  if (approach_direction != binding.approach_direction)
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint approach side differs from its exact arm/chart orientation");
  std::optional<std::int32_t> derived_rim;
  if (!object.at("derived_rim").is_null()) {
    derived_rim = as_i32(object.at("derived_rim"),
                         "checkpoint plan-bound endpoint rim");
    if (*derived_rim != -1 && *derived_rim != 1)
      throw std::invalid_argument(
          "checkpoint plan-bound endpoint rim must be +1 or -1");
  }
  if (derived_rim != binding.rim)
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint rim differs from its final-chart prescriptions");
  const auto cancellation_mode = required_string(
      object, "cancellation_mode");
  if (cancellation_mode != "exact-coefficient-field" &&
      cancellation_mode != "exact-or-acb-singleton")
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint cancellation mode is unsupported");
  auto analytic_metadata = as_object(
      object.at("analytic_metadata"),
      "checkpoint plan-bound endpoint analytic metadata");
  validate_checkpoint_exact_analytic_metadata(analytic_metadata);
  if (analytic_metadata != local->exact_analytic_metadata())
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint analytic metadata differs from its local owner");
  const auto provenance = planned_endpoint_provenance(
      checkpoint_identity, binding, cancellation_mode, analytic_metadata);
  if (json::serialize(canonical_json_value(provenance)) !=
      provenance_identity)
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint provenance identity is inconsistent");

  const auto& raw_result = as_object(
      object.at("result"), "checkpoint plan-bound endpoint result");
  require_exact_keys(
      raw_result,
      {"values", "dropped_regulated_sectors",
       "cancelled_divergent_coefficients", "imaginary_sign"},
      "checkpoint plan-bound endpoint result");
  EndpointLimitResult result;
  for (const auto& raw_value : as_array(
           raw_result.at("values"),
           "checkpoint plan-bound endpoint values"))
    result.values.push_back(parse_checkpoint_epsilon_frame<ComplexBall>(
        raw_value, "checkpoint plan-bound endpoint value"));
  if (result.values.empty())
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint result has no components");
  const auto checked_size = [](const json::value& value,
                               const char* label) {
    const auto parsed = as_u64(value, label);
    if (parsed > std::numeric_limits<std::size_t>::max())
      throw std::invalid_argument(std::string(label) + " exceeds size_t");
    return static_cast<std::size_t>(parsed);
  };
  result.dropped_regulated_sectors = checked_size(
      raw_result.at("dropped_regulated_sectors"),
      "checkpoint plan-bound dropped regulated sectors");
  result.cancelled_divergent_coefficients = checked_size(
      raw_result.at("cancelled_divergent_coefficients"),
      "checkpoint plan-bound cancelled divergent coefficients");
  std::optional<std::int32_t> result_rim;
  if (!raw_result.at("imaginary_sign").is_null()) {
    result_rim = as_i32(raw_result.at("imaginary_sign"),
                        "checkpoint plan-bound endpoint result rim");
    if (*result_rim != -1 && *result_rim != 1)
      throw std::invalid_argument(
          "checkpoint plan-bound endpoint result rim must be +1 or -1");
  }
  if (result_rim != binding.rim)
    throw std::invalid_argument(
        "checkpoint plan-bound endpoint result rim differs from its derived branch");
  // The endpoint kernel's limit classification is branch independent.  It
  // currently carries an internal integer diagnostic, while the plan-bound
  // public/checkpoint contract deliberately keeps an unprescribed rim null.
  result.imaginary_sign = binding.rim.value_or(1);
  (void)endpoint_value_window(result);

  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"),
      "checkpoint plan-bound endpoint elapsed time");
  const auto& stats = as_object(
      object.at("runtime_stats"),
      "checkpoint plan-bound endpoint runtime stats");
  require_exact_keys(stats, {"exports", "export_ms"},
                     "checkpoint plan-bound endpoint runtime stats");
  const auto exports = as_u64(
      stats.at("exports"), "checkpoint plan-bound endpoint exports");
  const auto export_ms = checkpoint_nonnegative_double(
      stats.at("export_ms"),
      "checkpoint plan-bound endpoint export time");
  auto endpoint = std::make_shared<StoredEndpointResult>(
      handle, checkpoint_identity, provenance_identity,
      local->handle(), local->source_chart(),
      local->source_operator_identity(), local->checkpoint_identity(),
      local->scalar_domain(), approach_direction, std::nullopt,
      cancellation_mode, std::move(analytic_metadata), std::move(result),
      elapsed_ms, binding.source, binding.rim,
      plan, local);
  endpoint->restore_runtime_stats(exports, export_ms);
  return endpoint;
}

std::shared_ptr<StoredEndpointResult>
restore_checkpoint_transport_endpoint_record(
    const json::value& raw, const std::string& expected_domain,
    const std::shared_ptr<StoredTransportArmState>& state) {
  const auto& object = as_object(
      raw, "checkpoint retained transport endpoint");
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "source", "approach_direction", "derived_rim",
       "cancellation_mode", "analytic_metadata", "result", "elapsed_ms",
       "runtime_stats"},
      "checkpoint retained transport endpoint");
  if (required_string(object, "schema") !=
      "diffexp2-retained-transport-endpoint-result-v1")
    throw std::invalid_argument(
        "unsupported retained transport endpoint checkpoint schema");
  if (!state || (expected_domain != "rational" && expected_domain != "acb") ||
      std::string(state->final_local()->scalar_domain()) != expected_domain)
    throw std::invalid_argument(
        "checkpoint transport endpoint lost its numeric state owner");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty())
    throw std::invalid_argument(
        "checkpoint transport endpoint contains an empty identity");

  const auto binding = resolve_transport_endpoint_binding(state);
  const auto& source = as_object(
      object.at("source"), "checkpoint transport endpoint source");
  auto expected_source = binding.source;
  const auto& observable = as_object(
      source.at("observable"), "checkpoint transport endpoint observable");
  require_exact_keys(observable, {"identity", "checkpoint_identity"},
                     "checkpoint transport endpoint observable");
  if (required_string(observable, "identity").empty() ||
      required_string(observable, "checkpoint_identity") !=
          checkpoint_identity)
    throw std::invalid_argument(
        "checkpoint transport endpoint observable identity is inconsistent");
  const auto& row_record = as_object(
      source.at("row"), "checkpoint transport endpoint row record");
  require_exact_keys(row_record, {"exact_identity", "prepared_row"},
                     "checkpoint transport endpoint row record");
  const auto& prepared_row = as_object(
      row_record.at("prepared_row"),
      "checkpoint transport endpoint prepared row");
  if (required_string(row_record, "exact_identity") !=
      required_string(prepared_row, "exact_identity"))
    throw std::invalid_argument(
        "checkpoint transport endpoint row identity is stale");
  const auto final_summary = state->final_local()->summary();
  validate_prepared_rational_row_structure(
      prepared_row,
      as_u32(final_summary.at("dimension"),
             "checkpoint transport endpoint source dimension"),
      "checkpoint transport endpoint prepared row");
  if (expected_domain == "rational") {
    const auto typed =
        std::dynamic_pointer_cast<StoredLocal<Rational>>(
            state->final_local());
    if (!typed)
      throw std::invalid_argument(
          "checkpoint transport endpoint lost its Rational final local");
    (void)parse_prepared_rational_row<Rational>(
        prepared_row, typed->solution());
  } else {
    const auto typed =
        std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(
            state->final_local());
    if (!typed)
      throw std::invalid_argument(
          "checkpoint transport endpoint lost its Acb final local");
    (void)parse_prepared_rational_row<ComplexBall>(
        prepared_row, typed->solution());
  }
  const auto& epsilon_record = as_object(
      source.at("output_epsilon_contract"),
      "checkpoint transport endpoint epsilon contract");
  const auto epsilon_contract = parse_observable_epsilon_contract(
      epsilon_record, "checkpoint transport endpoint epsilon contract");
  expected_source["observable"] = observable;
  expected_source["row"] = row_record;
  expected_source["output_epsilon_contract"] = epsilon_record;
  if (source != expected_source)
    throw std::invalid_argument(
        "checkpoint transport endpoint source differs from its exact state/plan/row binding");

  const auto approach_direction = as_i32(
      object.at("approach_direction"),
      "checkpoint transport endpoint approach direction");
  if (approach_direction != binding.approach_direction)
    throw std::invalid_argument(
        "checkpoint transport endpoint approach differs from its state plan");
  std::optional<std::int32_t> derived_rim;
  if (!object.at("derived_rim").is_null()) {
    derived_rim = as_i32(object.at("derived_rim"),
                         "checkpoint transport endpoint rim");
    if (*derived_rim != -1 && *derived_rim != 1)
      throw std::invalid_argument(
          "checkpoint transport endpoint rim must be +1 or -1");
  }
  if (derived_rim != binding.rim ||
      required_string(object, "cancellation_mode") !=
          "exact-or-acb-singleton")
    throw std::invalid_argument(
        "checkpoint transport endpoint rim or cancellation mode changed");
  auto analytic_metadata = as_object(
      object.at("analytic_metadata"),
      "checkpoint transport endpoint analytic metadata");
  validate_checkpoint_exact_analytic_metadata(analytic_metadata);
  if (analytic_metadata !=
      state->final_local()->exact_analytic_metadata())
    throw std::invalid_argument(
        "checkpoint transport endpoint analytic metadata differs from its retained final local");
  const auto provenance = transport_endpoint_provenance(
      checkpoint_identity, source, analytic_metadata);
  if (json::serialize(canonical_json_value(provenance)) !=
      provenance_identity)
    throw std::invalid_argument(
        "checkpoint transport endpoint provenance identity is inconsistent");

  const auto& raw_result = as_object(
      object.at("result"), "checkpoint transport endpoint result");
  require_exact_keys(
      raw_result,
      {"values", "dropped_regulated_sectors",
       "cancelled_divergent_coefficients", "imaginary_sign"},
      "checkpoint transport endpoint result");
  EndpointLimitResult result;
  for (const auto& raw_value : as_array(
           raw_result.at("values"),
           "checkpoint transport endpoint values"))
    result.values.push_back(parse_checkpoint_epsilon_frame<ComplexBall>(
        raw_value, "checkpoint transport endpoint value"));
  if (result.values.size() != 1)
    throw std::invalid_argument(
        "checkpoint transport endpoint result is not scalar");
  const auto checked_size = [](const json::value& value,
                               const char* label) {
    const auto parsed = as_u64(value, label);
    if (parsed > std::numeric_limits<std::size_t>::max())
      throw std::invalid_argument(std::string(label) + " exceeds size_t");
    return static_cast<std::size_t>(parsed);
  };
  result.dropped_regulated_sectors = checked_size(
      raw_result.at("dropped_regulated_sectors"),
      "checkpoint transport endpoint dropped sectors");
  result.cancelled_divergent_coefficients = checked_size(
      raw_result.at("cancelled_divergent_coefficients"),
      "checkpoint transport endpoint cancellations");
  std::optional<std::int32_t> result_rim;
  if (!raw_result.at("imaginary_sign").is_null())
    result_rim = as_i32(raw_result.at("imaginary_sign"),
                        "checkpoint transport endpoint result rim");
  if (result_rim != binding.rim)
    throw std::invalid_argument(
        "checkpoint transport endpoint result rim changed");
  result.imaginary_sign = binding.rim.value_or(1);
  const auto result_window = endpoint_value_window(result);
  if (result_window.min_power != epsilon_contract.requested.min_power ||
      result_window.complete_max >
          epsilon_contract.requested.complete_max ||
      result_window.complete_max <
          epsilon_contract.required_complete_max)
    throw std::invalid_argument(
        "checkpoint transport endpoint result violates its exact epsilon contract");

  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"),
      "checkpoint transport endpoint elapsed time");
  const auto& stats = as_object(
      object.at("runtime_stats"),
      "checkpoint transport endpoint runtime stats");
  require_exact_keys(stats, {"exports", "export_ms"},
                     "checkpoint transport endpoint runtime stats");
  auto endpoint = std::make_shared<StoredEndpointResult>(
      handle, checkpoint_identity, provenance_identity,
      state->final_local()->handle(), state->final_local()->source_chart(),
      state->final_local()->source_operator_identity(),
      state->final_local()->checkpoint_identity(),
      state->final_local()->scalar_domain(), approach_direction,
      std::nullopt, "exact-or-acb-singleton",
      std::move(analytic_metadata), std::move(result), elapsed_ms,
      source, binding.rim, nullptr, nullptr, state);
  endpoint->restore_runtime_stats(
      as_u64(stats.at("exports"),
             "checkpoint transport endpoint exports"),
      checkpoint_nonnegative_double(
          stats.at("export_ms"),
          "checkpoint transport endpoint export time"));
  return endpoint;
}

std::size_t checkpoint_size_t(const json::value& raw, const char* label);

std::shared_ptr<StoredPlannedMatchHop>
restore_checkpoint_planned_match_hop_record(
    const json::value& raw,
    const std::shared_ptr<StoredTilePlan>& plan,
    std::vector<std::shared_ptr<StoredLocalBase>> basis,
    std::shared_ptr<StoredLocalBase> incoming,
    const std::string& source_session_configuration_identity) {
  const auto& object = as_object(
      raw, "checkpoint retained planned match hop");
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "handoff", "native_match", "elapsed_ms", "runtime_stats"},
      "checkpoint retained planned match hop");
  if (required_string(object, "schema") !=
      "diffexp2-retained-planned-match-hop-v2")
    throw std::invalid_argument(
        "unsupported retained planned-match checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty() || !plan || !incoming || basis.empty())
    throw std::invalid_argument(
        "checkpoint planned match lost an identity or strong owner");
  const auto& handoff = as_object(
      object.at("handoff"), "checkpoint planned-match handoff");
  require_exact_keys(
      handoff,
      {"schema", "tile_plan", "tile_plan_checkpoint_identity",
       "tile_plan_provenance_identity", "arm", "match", "geometry",
       "producing", "receiving", "result_checkpoint_identity",
       "native_match_provenance_identity", "advance"},
      "checkpoint planned-match handoff");
  if (required_string(handoff, "schema") !=
          "diffexp2-retained-exact-plan-match-hop-v1" ||
      required_string(handoff, "tile_plan") != plan->handle() ||
      required_string(handoff, "tile_plan_checkpoint_identity") !=
          plan->checkpoint_identity() ||
      required_string(handoff, "tile_plan_provenance_identity") !=
          plan->provenance_identity() ||
      required_string(handoff, "result_checkpoint_identity") !=
          checkpoint_identity)
    throw std::invalid_argument(
        "checkpoint planned-match handoff lost its plan/checkpoint provenance");
  const auto arm_name = required_string(handoff, "arm");
  const auto match_index = checkpoint_size_t(
      handoff.at("match"), "checkpoint planned-match index");
  const auto& arm = plan->arm(arm_name);
  if (match_index >= arm.exact.matches.size())
    throw std::invalid_argument(
        "checkpoint planned-match index lies outside its retained plan");
  const auto& exact_match = arm.exact.matches[match_index];
  const auto& producing = arm.charts.at(exact_match.producing_chart);
  const auto& receiving = arm.charts.at(exact_match.receiving_chart);
  incoming->require_exact_plan_binding(
      producing.local_geometry, producing.prescriptions,
      "checkpoint planned-match incoming owner");
  for (const auto& local : basis)
    local->require_exact_plan_binding(
        receiving.local_geometry, receiving.prescriptions,
        "checkpoint planned-match basis owner");

  const auto& native_record = as_object(
      object.at("native_match"), "checkpoint embedded native match");
  const auto native_schema = required_string(native_record, "schema");
  std::shared_ptr<StoredMatchBase> native_match;
  if (native_schema == "diffexp2-retained-exact-rational-match-v2") {
    native_match = restore_checkpoint_exact_match_record(
        native_record, basis, incoming);
  } else if (native_schema == "diffexp2-retained-acb-match-v2") {
    const auto saturation = native_acb_saturation_binding(
        plan, source_session_configuration_identity, arm_name, match_index,
        checkpoint_identity);
    const std::optional<json::object> expected_singular_request =
        saturation.request_key == "native_singular_scc_saturation"
            ? std::optional<json::object>(saturation.request)
            : std::nullopt;
    native_match = restore_checkpoint_acb_match_record(
        native_record, source_session_configuration_identity,
        expected_singular_request);
    const auto cross_check = [](const json::object& source,
                                const std::shared_ptr<StoredLocalBase>& owner,
                                const char* label) {
      if (!owner || required_string(source, "local") != owner->handle() ||
          required_string(source, "chart") != owner->source_chart() ||
          required_string(source, "source_operator_identity") !=
              owner->source_operator_identity() ||
          required_string(source, "checkpoint_identity") !=
              owner->checkpoint_identity() ||
          source.at("analytic_metadata") != owner->exact_analytic_metadata())
        throw std::invalid_argument(
            std::string("checkpoint planned Acb match ") + label +
            " disagrees with its strong local owner");
    };
    const auto& sources = as_array(
        native_record.at("basis_sources"),
        "checkpoint planned Acb basis sources");
    if (sources.size() != basis.size())
      throw std::invalid_argument(
          "checkpoint planned Acb basis ownership differs from its dimension");
    for (std::size_t column = 0; column < basis.size(); ++column)
      cross_check(as_object(sources[column],
                            "checkpoint planned Acb basis source"),
                  basis[column], "basis source");
    cross_check(as_object(native_record.at("incoming_source"),
                          "checkpoint planned Acb incoming source"),
                incoming, "incoming source");
  } else {
    throw std::invalid_argument(
        "checkpoint planned hop embeds an unsupported native match kind");
  }
  if (native_match->handle() != handle ||
      native_match->checkpoint_record() != object.at("native_match"))
    throw std::invalid_argument(
        "checkpoint planned hop embedded match does not reproduce its exact payload");
  const auto native_summary = native_match->summary();
  const auto& certified_match = arm.certified_matches.at(match_index);
  if (required_string(native_summary, "checkpoint_identity") !=
          checkpoint_identity ||
      required_string(native_summary, "physical_match_point_exact") !=
          certified_match.receiving.physical_exact)
    throw std::invalid_argument(
        "checkpoint planned hop embedded match changed its checkpoint/physical point");
  std::vector<std::string> basis_handles;
  basis_handles.reserve(basis.size());
  for (const auto& local : basis) basis_handles.push_back(local->handle());
  auto expected_handoff = planned_match_handoff_record(
      plan, arm_name, match_index, basis_handles, basis, incoming->handle(),
      incoming, checkpoint_identity, native_summary);
  if (handoff != expected_handoff ||
      json::serialize(canonical_json_value(expected_handoff)) !=
          provenance_identity)
    throw std::invalid_argument(
        "checkpoint planned-match handoff/provenance differs from its exact owners");
  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint planned-match elapsed time");
  const auto& stats = as_object(
      object.at("runtime_stats"),
      "checkpoint planned-match runtime stats");
  require_exact_keys(stats, {"materializations"},
                     "checkpoint planned-match runtime stats");
  auto hop = std::make_shared<StoredPlannedMatchHop>(
      std::move(native_match), checkpoint_identity, provenance_identity,
      expected_handoff, elapsed_ms, plan, std::move(basis), incoming);
  hop->restore_runtime_stats(as_u64(
      stats.at("materializations"),
      "checkpoint planned-match materializations"));
  return hop;
}

std::size_t checkpoint_size_t(const json::value& raw, const char* label) {
  const auto value = as_u64(raw, label);
  if (value > std::numeric_limits<std::size_t>::max())
    throw std::invalid_argument(std::string(label) + " exceeds size_t");
  return static_cast<std::size_t>(value);
}

std::shared_ptr<StoredTransportArmState>
restore_checkpoint_transport_arm_state_record(
    const json::value& raw,
    const std::unordered_map<std::string,
                             std::shared_ptr<StoredTilePlan>>& plans,
    const std::unordered_map<std::string,
                             std::shared_ptr<StoredLocalBase>>& locals,
    const std::unordered_map<std::string,
                             std::shared_ptr<StoredMatchBase>>& matches) {
  const auto& object = as_object(
      raw, "checkpoint retained transport-arm state");
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "provenance", "elapsed_ms", "runtime_stats"},
      "checkpoint retained transport-arm state");
  const auto checkpoint_schema = required_string(object, "schema");
  const bool consumed_certificate_only =
      checkpoint_schema == "diffexp2-retained-transport-arm-state-v5";
  const bool consumed_compact =
      checkpoint_schema == "diffexp2-retained-transport-arm-state-v4";
  const bool compact_provenance =
      consumed_certificate_only || consumed_compact ||
      checkpoint_schema == "diffexp2-retained-transport-arm-state-v3";
  if (!compact_provenance && checkpoint_schema !=
      "diffexp2-retained-transport-arm-state-v2")
    throw std::invalid_argument(
        "unsupported retained transport-arm checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty())
    throw std::invalid_argument(
        "checkpoint transport-arm state contains an empty identity");
  const auto& provenance = as_object(
      object.at("provenance"),
      "checkpoint transport-arm state provenance");
  if (consumed_certificate_only)
    require_exact_keys(
        provenance,
        {"schema", "checkpoint_identity", "tile_plan", "arm",
         "tile_checkpoint_chain", "epsilon", "refinement"},
        "checkpoint certificate-only transport-arm state provenance");
  else
    require_exact_keys(
        provenance,
        {"schema", "checkpoint_identity", "tile_plan", "arm", "anchor",
         "receiving_basis", "matches", "tile_sources", "final_local",
         "epsilon", "refinement"},
        "checkpoint transport-arm state provenance");
  if (required_string(provenance, "schema") !=
          (compact_provenance
               ? (consumed_certificate_only
                      ? "diffexp2-retained-native-transport-arm-state-v4"
                      : consumed_compact
                      ? "diffexp2-retained-native-transport-arm-state-v3"
                      : "diffexp2-retained-native-transport-arm-state-v2")
               : "diffexp2-retained-native-transport-arm-state-v1") ||
      required_string(provenance, "checkpoint_identity") !=
          checkpoint_identity ||
      json::serialize(canonical_json_value(provenance)) !=
          provenance_identity)
    throw std::invalid_argument(
        "checkpoint transport-arm provenance identity is inconsistent");

  const auto& plan_reference = as_object(
      provenance.at("tile_plan"),
      "checkpoint transport-arm tile-plan reference");
  if (compact_provenance)
    require_exact_keys(
        plan_reference, {"handle", "checkpoint_identity"},
        "checkpoint transport-arm tile-plan reference");
  else
    require_exact_keys(
        plan_reference,
        {"handle", "checkpoint_identity", "provenance_identity"},
        "checkpoint transport-arm tile-plan reference");
  const auto plan_found = plans.find(
      required_string(plan_reference, "handle"));
  if (plan_found == plans.end() || !plan_found->second ||
      required_string(plan_reference, "checkpoint_identity") !=
          plan_found->second->checkpoint_identity() ||
      (!compact_provenance &&
       required_string(plan_reference, "provenance_identity") !=
           plan_found->second->provenance_identity()))
    throw std::invalid_argument(
        "checkpoint transport-arm state lost its retained tile plan");
  const auto plan = plan_found->second;
  const auto arm_name = required_string(provenance, "arm");
  (void)plan->arm(arm_name);

  const auto resolve_local_reference = [&](const json::value& raw_reference,
                                           const char* label,
                                           bool has_tile) {
    const auto& reference = as_object(raw_reference, label);
    if (has_tile)
      require_exact_keys(
          reference,
          {"tile", "local", "chart", "source_operator_identity",
           "checkpoint_identity", "coefficient_domain"},
          label);
    else
      require_exact_keys(
          reference,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity", "coefficient_domain"},
          label);
    const auto found = locals.find(required_string(reference, "local"));
    if (found == locals.end() || !found->second ||
        required_string(reference, "chart") !=
            found->second->source_chart() ||
        required_string(reference, "source_operator_identity") !=
            found->second->source_operator_identity() ||
        required_string(reference, "checkpoint_identity") !=
            found->second->checkpoint_identity() ||
        required_string(reference, "coefficient_domain") !=
            found->second->scalar_domain())
      throw std::invalid_argument(
          std::string(label) + " disagrees with its retained local owner");
    return found->second;
  };

  if (consumed_certificate_only) {
    std::vector<std::shared_ptr<StoredLocalBase>> tile_sources;
    const auto& chain = as_array(
        provenance.at("tile_checkpoint_chain"),
        "checkpoint certificate-only tile chain");
    const auto& retained = plan->arm(arm_name);
    if (chain.size() != retained.exact.tiles.size() || chain.empty())
      throw std::invalid_argument(
          "checkpoint certificate-only tile chain has the wrong topology");
    tile_sources.reserve(chain.size());
    for (std::size_t tile = 0; tile < chain.size(); ++tile) {
      const auto& reference = as_object(
          chain[tile], "checkpoint certificate-only tile reference");
      require_exact_keys(
          reference,
          {"tile", "local", "chart", "checkpoint_identity",
           "coefficient_domain", "derivation"},
          "checkpoint certificate-only tile reference");
      if (checkpoint_size_t(reference.at("tile"),
                            "checkpoint certificate-only tile index") !=
          tile)
        throw std::invalid_argument(
            "checkpoint certificate-only tiles are out of order");
      const auto found = locals.find(required_string(reference, "local"));
      if (found == locals.end() || !found->second ||
          required_string(reference, "chart") !=
              found->second->source_chart() ||
          required_string(reference, "checkpoint_identity") !=
              found->second->checkpoint_identity() ||
          required_string(reference, "coefficient_domain") !=
              found->second->scalar_domain())
        throw std::invalid_argument(
            "checkpoint certificate-only tile reference disagrees with its local");
      if (tile == 0) {
        if (!reference.at("derivation").is_null())
          throw std::invalid_argument(
              "checkpoint certificate-only anchor unexpectedly has a derivation");
      } else {
        const auto& certificate = as_object(
            reference.at("derivation"),
            "checkpoint certificate-only hop derivation");
        const auto certificate_schema = required_string(
            certificate, "schema");
        if (certificate_schema ==
            "diffexp2-consumed-plan-value-handoff-certificate-v1") {
          require_exact_keys(
              certificate,
              {"schema", "match", "handoff_provenance_identity",
               "incoming_checkpoint_identity", "output_checkpoint_identity"},
              "checkpoint certificate-only value-hop derivation");
          const auto& local_derivation = found->second->retained_derivation();
          if (!local_derivation.has_value() ||
              !is_retained_plan_value_handoff_schema(
                  required_string(*local_derivation, "schema")) ||
              required_string(certificate,
                              "handoff_provenance_identity") !=
                  required_string(*local_derivation,
                                  "provenance_identity"))
            throw std::invalid_argument(
                "checkpoint certificate-only value hop lost its sealed handoff identity");
        } else {
          require_exact_keys(
              certificate,
              {"schema", "match", "source_match_checkpoint_identity",
               "incoming_checkpoint_identity", "output_checkpoint_identity"},
              "checkpoint certificate-only match-hop derivation");
          if (certificate_schema !=
                  "diffexp2-consumed-plan-match-certificate-v1" ||
              required_string(certificate,
                              "source_match_checkpoint_identity").empty())
            throw std::invalid_argument(
                "checkpoint certificate-only match hop lost its match identity");
        }
        if (checkpoint_size_t(certificate.at("match"),
                              "checkpoint certificate-only match index") !=
                tile - 1 ||
            required_string(certificate,
                            "incoming_checkpoint_identity") !=
                tile_sources.back()->checkpoint_identity() ||
            required_string(certificate,
                            "output_checkpoint_identity") !=
                found->second->checkpoint_identity())
          throw std::invalid_argument(
              "checkpoint certificate-only hop breaks its checkpoint chain");
      }
      tile_sources.push_back(found->second);
    }
    auto anchor = tile_sources.front();
    const auto& epsilon = as_object(
        provenance.at("epsilon"),
        "checkpoint certificate-only epsilon contract");
    require_exact_keys(
        epsilon,
        {"min", "max", "required_complete_max",
         "match_required_complete_max"},
        "checkpoint certificate-only epsilon contract");
    EpsilonWindow work_epsilon{
        as_i32(epsilon.at("min"),
               "checkpoint certificate-only epsilon minimum"),
        as_i32(epsilon.at("max"),
               "checkpoint certificate-only epsilon maximum")};
    (void)work_epsilon.width();
    const auto public_required_complete_max = as_i32(
        epsilon.at("required_complete_max"),
        "checkpoint certificate-only public epsilon maximum");
    const auto match_required_complete_max = as_i32(
        epsilon.at("match_required_complete_max"),
        "checkpoint certificate-only match epsilon maximum");
    const auto refinement = as_object(
        provenance.at("refinement"),
        "checkpoint certificate-only refinement policy");
    const auto elapsed_ms = checkpoint_nonnegative_double(
        object.at("elapsed_ms"),
        "checkpoint certificate-only transport elapsed time");
    auto state = std::make_shared<StoredTransportArmState>(
        handle, checkpoint_identity, arm_name, plan, anchor,
        std::move(tile_sources), work_epsilon,
        public_required_complete_max, match_required_complete_max,
        refinement, elapsed_ms);
    if (state->provenance_identity() != provenance_identity)
      throw std::invalid_argument(
          "restored certificate-only transport state changed provenance");
    const auto& stats = as_object(
        object.at("runtime_stats"),
        "checkpoint certificate-only transport runtime stats");
    require_exact_keys(stats,
                       {"stats_queries", "contraction_operations",
                        "contracted_observables",
                        "endpoint_batch_operations", "endpoint_rows"},
                       "checkpoint certificate-only runtime stats");
    state->restore_runtime_stats(
        as_u64(stats.at("stats_queries"),
               "checkpoint certificate-only stats queries"),
        as_u64(stats.at("contraction_operations"),
               "checkpoint certificate-only contractions"),
        as_u64(stats.at("contracted_observables"),
               "checkpoint certificate-only observables"),
        as_u64(stats.at("endpoint_batch_operations"),
               "checkpoint certificate-only endpoint operations"),
        as_u64(stats.at("endpoint_rows"),
               "checkpoint certificate-only endpoint rows"));
    if (state->checkpoint_record() != raw)
      throw std::invalid_argument(
          "restored certificate-only transport state changed its checkpoint record");
    return state;
  }

  auto anchor = resolve_local_reference(
      provenance.at("anchor"),
      "checkpoint transport-arm anchor reference", false);
  std::vector<std::vector<std::shared_ptr<StoredLocalBase>>> basis;
  json::array cached_basis_references;
  for (const auto& raw_basis : as_array(
           provenance.at("receiving_basis"),
           "checkpoint transport-arm receiving bases")) {
    const auto& raw_columns = as_array(
        raw_basis, "checkpoint transport-arm receiving basis");
    if (raw_columns.empty())
      throw std::invalid_argument(
          "checkpoint transport-arm receiving basis cannot be empty");
    if (consumed_compact) {
      for (const auto& raw_column : raw_columns) {
        const auto& reference = as_object(
            raw_column, "checkpoint consumed transport-arm basis reference");
        const bool bounded_operator =
            reference.if_contains("source_operator_reference") != nullptr;
        if (bounded_operator)
          require_exact_keys(
              reference,
              {"local", "chart", "source_operator_reference",
               "checkpoint_identity", "coefficient_domain"},
              "checkpoint consumed bounded transport-arm basis reference");
        else
          require_exact_keys(
              reference,
              {"local", "chart", "source_operator_identity",
               "checkpoint_identity", "coefficient_domain"},
              "checkpoint consumed legacy transport-arm basis reference");
        (void)scoped_handle_id(required_string(reference, "local"), "l:",
                               "checkpoint consumed basis local");
        (void)required_string(reference, "chart");
        if (bounded_operator)
          (void)as_object(
              reference.at("source_operator_reference"),
              "checkpoint consumed basis source-operator reference");
        else
          (void)required_string(reference, "source_operator_identity");
        (void)required_string(reference, "checkpoint_identity");
        (void)required_string(reference, "coefficient_domain");
      }
      cached_basis_references.push_back(raw_basis);
    } else {
      std::vector<std::shared_ptr<StoredLocalBase>> columns;
      for (const auto& raw_column : raw_columns)
        columns.push_back(resolve_local_reference(
            raw_column, "checkpoint transport-arm basis reference", false));
      basis.push_back(std::move(columns));
    }
  }

  std::vector<std::shared_ptr<StoredPlannedMatchHop>> planned_matches;
  json::array cached_match_references;
  const auto& raw_matches = as_array(
      provenance.at("matches"), "checkpoint transport-arm matches");
  planned_matches.reserve(raw_matches.size());
  for (std::size_t index = 0; index < raw_matches.size(); ++index) {
    const auto& reference = as_object(
        raw_matches[index], "checkpoint transport-arm match reference");
    if (consumed_compact)
      require_exact_keys(
          reference,
          {"index", "checkpoint_identity", "provenance_identity",
           "planned_hop", "sealed_local_lineage"},
          "checkpoint consumed transport-arm match reference");
    else if (compact_provenance)
      require_exact_keys(
          reference, {"index", "match", "checkpoint_identity"},
          "checkpoint transport-arm match reference");
    else
      require_exact_keys(
          reference,
          {"index", "match", "checkpoint_identity", "provenance_identity"},
          "checkpoint transport-arm match reference");
    if (checkpoint_size_t(reference.at("index"),
                          "checkpoint transport-arm match index") != index)
      throw std::invalid_argument(
          "checkpoint transport-arm matches are out of order");
    if (consumed_compact) {
      (void)required_string(reference, "checkpoint_identity");
      (void)required_string(reference, "provenance_identity");
      (void)as_object(reference.at("planned_hop"),
                      "checkpoint consumed transport-arm planned hop");
      (void)as_object(reference.at("sealed_local_lineage"),
                      "checkpoint consumed transport-arm sealed lineage");
      cached_match_references.push_back(raw_matches[index]);
      continue;
    }
    const auto found = matches.find(required_string(reference, "match"));
    if (found == matches.end())
      throw std::invalid_argument(
          "checkpoint transport-arm state lost a planned match");
    auto match = std::dynamic_pointer_cast<StoredPlannedMatchHop>(
        found->second);
    if (!match ||
        required_string(reference, "checkpoint_identity") !=
            match->checkpoint_identity() ||
        (!compact_provenance &&
         required_string(reference, "provenance_identity") !=
             match->provenance_identity()))
      throw std::invalid_argument(
          "checkpoint transport-arm match reference is inconsistent");
    planned_matches.push_back(std::move(match));
  }

  std::vector<std::shared_ptr<StoredLocalBase>> tile_sources;
  const auto& raw_sources = as_array(
      provenance.at("tile_sources"),
      "checkpoint transport-arm tile sources");
  tile_sources.reserve(raw_sources.size());
  for (std::size_t tile = 0; tile < raw_sources.size(); ++tile) {
    const auto& reference = as_object(
        raw_sources[tile], "checkpoint transport-arm tile-source reference");
    if (checkpoint_size_t(reference.at("tile"),
                          "checkpoint transport-arm tile index") != tile)
      throw std::invalid_argument(
          "checkpoint transport-arm tile sources are out of order");
    tile_sources.push_back(resolve_local_reference(
        raw_sources[tile],
        "checkpoint transport-arm tile-source reference", true));
  }
  auto final_local = resolve_local_reference(
      provenance.at("final_local"),
      "checkpoint transport-arm final-local reference", false);
  if (tile_sources.empty() || tile_sources.back().get() != final_local.get())
    throw std::invalid_argument(
        "checkpoint transport-arm final local differs from its last tile source");

  const auto& epsilon = as_object(
      provenance.at("epsilon"), "checkpoint transport-arm epsilon contract");
  require_exact_keys(
      epsilon,
      {"min", "max", "required_complete_max",
       "match_required_complete_max"},
      "checkpoint transport-arm epsilon contract");
  EpsilonWindow work_epsilon{
      as_i32(epsilon.at("min"), "checkpoint transport-arm epsilon minimum"),
      as_i32(epsilon.at("max"), "checkpoint transport-arm epsilon maximum")};
  (void)work_epsilon.width();
  const auto public_required_complete_max = as_i32(
      epsilon.at("required_complete_max"),
      "checkpoint transport-arm public epsilon maximum");
  const auto match_required_complete_max = as_i32(
      epsilon.at("match_required_complete_max"),
      "checkpoint transport-arm match epsilon maximum");
  const auto refinement = as_object(
      provenance.at("refinement"),
      "checkpoint transport-arm refinement policy");
  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint transport-arm elapsed time");
  std::shared_ptr<StoredTransportArmState> state;
  if (consumed_compact)
    state = std::make_shared<StoredTransportArmState>(
        handle, checkpoint_identity, arm_name, plan, anchor,
        std::move(cached_basis_references),
        std::move(cached_match_references), std::move(tile_sources),
        work_epsilon, public_required_complete_max,
        match_required_complete_max, refinement, elapsed_ms);
  else
    state = std::make_shared<StoredTransportArmState>(
        handle, checkpoint_identity, arm_name, plan, anchor,
        std::move(basis), std::move(planned_matches),
        std::move(tile_sources), work_epsilon,
        public_required_complete_max, match_required_complete_max,
        refinement, elapsed_ms, compact_provenance);
  if (state->provenance_identity() != provenance_identity)
    throw std::invalid_argument(
        "restored transport-arm state changed its exact provenance");
  const auto& stats = as_object(
      object.at("runtime_stats"),
      "checkpoint transport-arm runtime stats");
  require_exact_keys(stats,
                     {"stats_queries", "contraction_operations",
                      "contracted_observables",
                      "endpoint_batch_operations", "endpoint_rows"},
                     "checkpoint transport-arm runtime stats");
  state->restore_runtime_stats(
      as_u64(stats.at("stats_queries"),
             "checkpoint transport-arm statistics queries"),
      as_u64(stats.at("contraction_operations"),
             "checkpoint transport-arm contraction operations"),
      as_u64(stats.at("contracted_observables"),
             "checkpoint transport-arm contracted observables"),
      as_u64(stats.at("endpoint_batch_operations"),
             "checkpoint transport-arm endpoint-batch operations"),
      as_u64(stats.at("endpoint_rows"),
             "checkpoint transport-arm endpoint rows"));
  if (state->checkpoint_record() != raw)
    throw std::invalid_argument(
        "restored transport-arm state does not reproduce its exact retained state");
  return state;
}

std::shared_ptr<StoredLineResult> restore_checkpoint_line_result_record(
    const json::value& raw,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& local) {
  const auto& object = as_object(raw, "checkpoint retained line result");
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "arm", "tile", "interval", "source", "result", "elapsed_ms",
       "runtime_stats"},
      "checkpoint retained line result");
  if (required_string(object, "schema") !=
      "diffexp2-retained-line-result-v2")
    throw std::invalid_argument(
        "unsupported retained line-result checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  const auto arm_name = required_string(object, "arm");
  const auto tile_index = checkpoint_size_t(
      object.at("tile"), "checkpoint line tile index");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty())
    throw std::invalid_argument(
        "checkpoint retained line result contains an empty identity");
  const auto& arm = plan->arm(arm_name);
  if (tile_index >= arm.exact.tiles.size())
    throw std::invalid_argument(
        "checkpoint line tile lies outside its retained plan owner");
  const auto expected_interval = encode_plan_tile(arm, tile_index);
  const auto interval = as_object(object.at("interval"),
                                  "checkpoint line interval");
  if (interval != expected_interval)
    throw std::invalid_argument(
        "checkpoint line interval differs from its exact tile-plan owner");
  const auto& source = as_object(object.at("source"),
                                 "checkpoint line source");
  require_exact_keys(
      source,
      {"tile_plan", "tile_plan_checkpoint_identity", "local", "chart",
       "source_operator_identity", "local_checkpoint_identity",
       "coefficient_domain", "analytic_metadata"},
      "checkpoint line source");
  if (!plan || !local ||
      required_string(source, "tile_plan") != plan->handle() ||
      required_string(source, "tile_plan_checkpoint_identity") !=
          plan->checkpoint_identity() ||
      required_string(source, "local") != local->handle() ||
      required_string(source, "chart") != local->source_chart() ||
      required_string(source, "source_operator_identity") !=
          local->source_operator_identity() ||
      required_string(source, "local_checkpoint_identity") !=
          local->checkpoint_identity() ||
      required_string(source, "coefficient_domain") !=
          local->scalar_domain() ||
      source.at("analytic_metadata") != local->exact_analytic_metadata())
    throw std::invalid_argument(
        "checkpoint line source provenance disagrees with its strong owners");
  if (std::string(local->scalar_domain()) == "symbolic")
    throw std::invalid_argument(
        "checkpoint line result cannot own a symbolic local");
  validate_checkpoint_exact_analytic_metadata(
      source.at("analytic_metadata"));
  const auto& tile = arm.exact.tiles[tile_index];
  const auto& binding = arm.charts.at(tile.chart);
  if (local->source_chart() != binding.handle)
    throw std::invalid_argument(
        "checkpoint line local does not own the tile's retained chart");

  const auto& raw_result = as_object(object.at("result"),
                                     "checkpoint line result state");
  require_exact_keys(raw_result, {"value", "scope", "imaginary_sign",
                                  "diagnostics"},
                     "checkpoint line result state");
  const auto scope = required_string(raw_result, "scope");
  StoredLineIntegral result;
  if (scope == "stored_truncation")
    result.scope = LineIntegrationScope::StoredTruncation;
  else if (scope == "full_local_with_certified_tail")
    result.scope = LineIntegrationScope::FullLocalWithCertifiedTail;
  else
    throw std::invalid_argument(
        "checkpoint line result has an unsupported integration scope");
  result.value = parse_checkpoint_epsilon_vector(
      raw_result.at("value"), "checkpoint line epsilon vector");
  if (!raw_result.at("imaginary_sign").is_null()) {
    result.imaginary_sign = as_i32(raw_result.at("imaginary_sign"),
                                   "checkpoint line rim");
    if (*result.imaginary_sign != -1 && *result.imaginary_sign != 1)
      throw std::invalid_argument(
          "checkpoint line rim must be +1 or -1");
  }
  const auto expected_rim = exact_plan_rim(
      binding.prescriptions, binding.geometry.scale);
  if (result.imaginary_sign != expected_rim)
    throw std::invalid_argument(
        "checkpoint line rim differs from its exact branch prescriptions");
  if (result.value.dimension !=
      as_u32(local->summary().at("dimension"),
             "checkpoint line source dimension"))
    throw std::invalid_argument(
        "checkpoint line result dimension differs from its local owner");
  const auto& diagnostics = as_object(raw_result.at("diagnostics"),
                                      "checkpoint line diagnostics");
  const bool has_tail_requested =
      diagnostics.if_contains("tail_certificate_requested") != nullptr;
  const bool has_tail_status =
      diagnostics.if_contains("tail_certificate_status") != nullptr;
  const bool has_tail_witness =
      diagnostics.if_contains("tail_witness_radius_exact") != nullptr;
  const std::array<const char*, 4> cancellation_diagnostic_keys{
      "bounded_cancelled_divergent_coefficients",
      "divergent_cancellation_mode", "divergent_relative_tolerance",
      "divergent_cancellation_provenance"};
  const auto cancellation_diagnostic_count = std::count_if(
      cancellation_diagnostic_keys.begin(),
      cancellation_diagnostic_keys.end(), [&](const char* key) {
        return diagnostics.if_contains(key) != nullptr;
      });
  const bool has_cancellation_diagnostics =
      cancellation_diagnostic_count == cancellation_diagnostic_keys.size();
  if (cancellation_diagnostic_count != 0 &&
      !has_cancellation_diagnostics)
    throw std::invalid_argument(
        "checkpoint line divergent-cancellation diagnostics are incomplete");
  if (has_tail_requested != has_tail_status ||
      has_tail_requested != has_tail_witness)
    throw std::invalid_argument(
        "checkpoint line tail diagnostics are incomplete");
  if (has_tail_requested)
    if (has_cancellation_diagnostics)
      require_exact_keys(
          diagnostics,
          {"input_monomial_cells", "grouped_monomials", "zero_groups_skipped",
           "cancelled_divergent_groups",
           "bounded_cancelled_divergent_coefficients",
           "divergent_cancellation_mode", "divergent_relative_tolerance",
           "divergent_cancellation_provenance", "primitive_evaluations",
           "primitive_component_applications", "primitive_component_reuses",
           "has_center_endpoint", "tail_certificate_requested",
           "tail_certificate_status", "tail_witness_radius_exact", "detail"},
          "checkpoint line diagnostics");
    else
      require_exact_keys(
          diagnostics,
          {"input_monomial_cells", "grouped_monomials", "zero_groups_skipped",
           "cancelled_divergent_groups", "primitive_evaluations",
           "primitive_component_applications", "primitive_component_reuses",
           "has_center_endpoint", "tail_certificate_requested",
           "tail_certificate_status", "tail_witness_radius_exact", "detail"},
          "checkpoint legacy line diagnostics");
  else
    if (has_cancellation_diagnostics)
      require_exact_keys(
          diagnostics,
          {"input_monomial_cells", "grouped_monomials", "zero_groups_skipped",
           "cancelled_divergent_groups",
           "bounded_cancelled_divergent_coefficients",
           "divergent_cancellation_mode", "divergent_relative_tolerance",
           "divergent_cancellation_provenance", "primitive_evaluations",
           "primitive_component_applications", "primitive_component_reuses",
           "has_center_endpoint", "detail"},
          "checkpoint line diagnostics");
    else
      require_exact_keys(
          diagnostics,
          {"input_monomial_cells", "grouped_monomials", "zero_groups_skipped",
           "cancelled_divergent_groups", "primitive_evaluations",
           "primitive_component_applications", "primitive_component_reuses",
           "has_center_endpoint", "detail"},
          "checkpoint legacy line diagnostics");
  result.diagnostics.input_monomial_cells = checkpoint_size_t(
      diagnostics.at("input_monomial_cells"),
      "checkpoint line input monomials");
  result.diagnostics.grouped_monomials = checkpoint_size_t(
      diagnostics.at("grouped_monomials"),
      "checkpoint line grouped monomials");
  result.diagnostics.zero_groups_skipped = checkpoint_size_t(
      diagnostics.at("zero_groups_skipped"),
      "checkpoint line skipped zero groups");
  result.diagnostics.cancelled_divergent_groups = checkpoint_size_t(
      diagnostics.at("cancelled_divergent_groups"),
      "checkpoint line cancelled divergent groups");
  if (has_cancellation_diagnostics) {
    result.diagnostics.bounded_cancelled_divergent_coefficients =
        checkpoint_size_t(
            diagnostics.at("bounded_cancelled_divergent_coefficients"),
            "checkpoint line bounded-cancelled divergent coefficients");
    result.diagnostics.divergent_cancellation_mode = required_string(
        diagnostics, "divergent_cancellation_mode");
    if (!diagnostics.at("divergent_relative_tolerance").is_string() ||
        !diagnostics.at("divergent_cancellation_provenance").is_string())
      throw std::invalid_argument(
          "checkpoint line divergent-cancellation strings are malformed");
    result.diagnostics.divergent_relative_tolerance = std::string(
        diagnostics.at("divergent_relative_tolerance").as_string());
    result.diagnostics.divergent_cancellation_provenance = std::string(
        diagnostics.at("divergent_cancellation_provenance").as_string());
    if (result.diagnostics.divergent_cancellation_mode != "exact-singleton" &&
        result.diagnostics.divergent_cancellation_mode !=
            "bounded-relative-acb")
      throw std::invalid_argument(
          "checkpoint line has an unsupported divergent-cancellation mode");
    const bool bounded_cancellation =
        result.diagnostics.divergent_cancellation_mode ==
        "bounded-relative-acb";
    const bool has_cancellation_tolerance =
        !result.diagnostics.divergent_relative_tolerance.empty();
    const bool has_cancellation_provenance =
        !result.diagnostics.divergent_cancellation_provenance.empty();
    if ((bounded_cancellation &&
         (!has_cancellation_tolerance || !has_cancellation_provenance)) ||
        (!bounded_cancellation &&
         (has_cancellation_tolerance || has_cancellation_provenance)))
      throw std::invalid_argument(
          "checkpoint line divergent-cancellation policy is incomplete");
    if (!bounded_cancellation &&
        result.diagnostics.bounded_cancelled_divergent_coefficients != 0)
      throw std::invalid_argument(
          "checkpoint exact-singleton line cannot report bounded-cancelled "
          "divergent coefficients");
  }
  result.diagnostics.primitive_evaluations = checkpoint_size_t(
      diagnostics.at("primitive_evaluations"),
      "checkpoint line primitive evaluations");
  result.diagnostics.primitive_component_applications = checkpoint_size_t(
      diagnostics.at("primitive_component_applications"),
      "checkpoint line primitive component applications");
  result.diagnostics.primitive_component_reuses = checkpoint_size_t(
      diagnostics.at("primitive_component_reuses"),
      "checkpoint line primitive component reuses");
  if (!diagnostics.at("has_center_endpoint").is_bool())
    throw std::invalid_argument(
        "checkpoint line center-endpoint flag must be Boolean");
  result.diagnostics.has_center_endpoint =
      diagnostics.at("has_center_endpoint").as_bool();
  if (has_tail_requested) {
    if (!diagnostics.at("tail_certificate_requested").is_bool())
      throw std::invalid_argument(
          "checkpoint line tail-request flag must be Boolean");
    result.diagnostics.tail_certificate_requested =
        diagnostics.at("tail_certificate_requested").as_bool();
    result.diagnostics.tail_certificate_status = required_string(
        diagnostics, "tail_certificate_status");
    if (!diagnostics.at("tail_witness_radius_exact").is_null())
      result.diagnostics.tail_witness_radius_exact = required_string(
          diagnostics, "tail_witness_radius_exact");
  }
  result.diagnostics.detail = required_string(diagnostics, "detail");

  const auto& error = result.value.error;
  const auto& line_epsilon = result.value.epsilon;
  const auto& tail = result.diagnostics;
  if (!tail.tail_witness_radius_exact.empty()) {
    try {
      const Rational witness(tail.tail_witness_radius_exact);
      const auto begin_modulus = tile.local_begin.sign() < 0
          ? -tile.local_begin : tile.local_begin;
      const auto end_modulus = tile.local_end.sign() < 0
          ? -tile.local_end : tile.local_end;
      const auto outer = begin_modulus < end_modulus
          ? end_modulus : begin_modulus;
      if (!(outer < witness) || !(witness < binding.geometry.radius))
        throw std::invalid_argument("tail witness lies outside its annulus");
    } catch (const std::invalid_argument&) {
      throw std::invalid_argument(
          "checkpoint line tail witness radius must be an exact rational "
          "strictly outside its tile and inside its chart");
    }
  }
  if (result.scope == LineIntegrationScope::StoredTruncation) {
    if (!error.empty() || error.guarantee != ErrorGuarantee::None)
      throw std::invalid_argument(
          "checkpoint stored-truncation line cannot carry an error "
          "envelope or guarantee");
    if (tail.tail_certificate_requested) {
      if (tail.tail_certificate_status != "unsupported" &&
          tail.tail_certificate_status != "inconclusive")
        throw std::invalid_argument(
            "checkpoint stored-truncation line has an inconsistent "
            "tail-certificate status");
    } else if (tail.tail_certificate_status != "not-requested" ||
               !tail.tail_witness_radius_exact.empty()) {
      throw std::invalid_argument(
          "checkpoint stored-truncation line has unsolicited "
          "tail-certificate diagnostics");
    }
  } else {
    if (error.empty() || error.guarantee != ErrorGuarantee::Certified ||
        error.frame.min_power != line_epsilon.min_power ||
        error.frame.complete_max != line_epsilon.complete_max ||
        error.provenance.empty())
      throw std::invalid_argument(
          "checkpoint full-local line requires a frame-aligned certified "
          "error envelope");
    if (!tail.tail_certificate_requested ||
        tail.tail_certificate_status != "certified" ||
        tail.tail_witness_radius_exact.empty())
      throw std::invalid_argument(
          "checkpoint full-local line requires complete certified-tail "
          "diagnostics");
  }

  json::object legacy_provenance{
      {"schema",
       "diffexp2-retained-native-stored-truncation-physical-tile-integral-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
      {"arm", arm_name}, {"tile", tile_index},
      {"interval", interval},
      {"source", json::object{
           {"local", local->handle()}, {"chart", local->source_chart()},
           {"checkpoint_identity", local->checkpoint_identity()}}},
      {"epsilon", json::object{{"min", result.value.epsilon.min_power},
                                {"max", result.value.epsilon.complete_max}}},
      {"scope", "stored_truncation"},
      {"error_guarantee", "none"}};
  json::object provenance{
      {"schema",
       "diffexp2-retained-native-physical-tile-integral-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
      {"arm", arm_name}, {"tile", tile_index},
      {"interval", interval},
      {"source", json::object{
           {"local", local->handle()}, {"chart", local->source_chart()},
           {"checkpoint_identity", local->checkpoint_identity()}}},
      {"epsilon", json::object{{"min", result.value.epsilon.min_power},
                                {"max", result.value.epsilon.complete_max}}},
      {"tail_certificate_requested",
       result.diagnostics.tail_certificate_requested},
      {"tail_certificate_status",
       result.diagnostics.tail_certificate_status},
      {"tail_witness_radius_exact",
       result.diagnostics.tail_witness_radius_exact.empty()
           ? json::value(nullptr)
           : json::value(result.diagnostics.tail_witness_radius_exact)},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)},
      {"error_provenance", result.value.error.provenance}};
  if (has_cancellation_diagnostics) {
    if (result.diagnostics.divergent_cancellation_mode ==
        "exact-singleton") {
      provenance["divergent_cancellation"] =
          json::object{{"mode", "exact-singleton"}};
    } else {
      json::object policy{
          {"mode", "bounded-relative-acb"},
          {"relative_tolerance",
           result.diagnostics.divergent_relative_tolerance},
          {"provenance",
           result.diagnostics.divergent_cancellation_provenance}};
      provenance["divergent_cancellation"] =
          encode_bounded_divergent_cancellation(
              parse_bounded_divergent_cancellation(
                  policy,
                  "checkpoint line divergent-cancellation policy"));
    }
  }
  const auto current_identity =
      json::serialize(canonical_json_value(provenance));
  const auto legacy_identity =
      json::serialize(canonical_json_value(legacy_provenance));
  if (provenance_identity != current_identity &&
      (has_tail_requested ||
       result.scope != LineIntegrationScope::StoredTruncation ||
       provenance_identity != legacy_identity))
    throw std::invalid_argument(
        "checkpoint line provenance identity is inconsistent");
  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint line elapsed time");
  const auto& stats = as_object(object.at("runtime_stats"),
                                "checkpoint line runtime stats");
  require_exact_keys(stats, {"exports", "export_ms"},
                     "checkpoint line runtime stats");
  auto stored = std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, arm_name, tile_index,
      interval, local->checkpoint_identity(), std::move(result), elapsed_ms,
      plan, local);
  stored->restore_runtime_stats(
      as_u64(stats.at("exports"), "checkpoint line exports"),
      checkpoint_nonnegative_double(stats.at("export_ms"),
                                    "checkpoint line export time"));
  return stored;
}

std::shared_ptr<StoredLineResult> restore_checkpoint_line_aggregate_record(
    const json::value& raw,
    std::vector<std::shared_ptr<StoredTilePlan>> plans,
    std::vector<std::shared_ptr<StoredLocalBase>> local_owners,
    std::vector<std::shared_ptr<StoredTransportArmState>>
        transport_owners) {
  const auto& object = as_object(raw, "checkpoint retained line aggregate");
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "provenance", "result", "elapsed_ms", "runtime_stats"},
      "checkpoint retained line aggregate");
  const auto schema = required_string(object, "schema");
  const bool compact_single_transport =
      schema == "diffexp2-retained-transport-observable-line-v1" ||
      schema == "diffexp2-retained-transport-observable-line-v2";
  const bool compact_pair_transport =
      schema ==
          "diffexp2-retained-transport-pair-observable-line-v1" ||
      schema ==
          "diffexp2-retained-transport-pair-observable-line-v2";
  const bool compact_transport =
      compact_single_transport || compact_pair_transport;
  const bool compact_v2 =
      schema == "diffexp2-retained-transport-observable-line-v2" ||
      schema == "diffexp2-retained-transport-pair-observable-line-v2";
  if (!compact_transport && schema != "diffexp2-retained-line-aggregate-v1")
    throw std::invalid_argument(
        "unsupported retained line-aggregate checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty() || plans.empty() ||
      std::any_of(plans.begin(), plans.end(),
                  [](const auto& owner) { return owner == nullptr; }) ||
      (local_owners.empty() && !compact_transport))
    throw std::invalid_argument(
        "checkpoint line aggregate lost an identity or strong owner");
  local_owners = unique_line_local_owners(local_owners);

  const auto& provenance = as_object(
      object.at("provenance"), "checkpoint line aggregate provenance");
  require_exact_keys(
      provenance,
      {"schema", "checkpoint_identity", "arm", "interval", "source",
       "aggregate", "epsilon", "scope", "error_guarantee"},
      "checkpoint line aggregate provenance");
  if (required_string(provenance, "schema") !=
          (compact_pair_transport
               ? compact_v2
                   ? "diffexp2-retained-native-transport-pair-observable-line-v2"
                   : "diffexp2-retained-native-transport-pair-observable-line-v1"
               : compact_single_transport
               ? compact_v2
                   ? "diffexp2-retained-native-transport-observable-line-v2"
                   : "diffexp2-retained-native-transport-observable-line-v1"
               : "diffexp2-retained-native-line-aggregate-v1") ||
      required_string(provenance, "checkpoint_identity") !=
          checkpoint_identity ||
      json::serialize(canonical_json_value(provenance)) !=
          provenance_identity)
    throw std::invalid_argument(
        "checkpoint line aggregate provenance identity is inconsistent");
  const auto& source = as_object(
      provenance.at("source"), "checkpoint line aggregate source");
  if (compact_single_transport)
    if (compact_v2)
      require_exact_keys(
          source,
          {"tile_plan", "tile_plan_checkpoint_identity", "transport_state",
           "transport_state_checkpoint_identity"},
          "checkpoint compact transport line source");
    else
      require_exact_keys(
          source,
          {"tile_plan", "tile_plan_checkpoint_identity", "transport_state",
           "transport_state_checkpoint_identity",
           "transport_state_provenance_identity"},
          "checkpoint compact transport line source");
  else if (compact_pair_transport)
    require_exact_keys(
        source, {"lower", "upper", "common_anchor"},
        "checkpoint compact transport-pair line source");
  else
    require_exact_keys(source,
        {"tile_plan", "tile_plan_checkpoint_identity", "locals"},
        "checkpoint line aggregate source");
  if (!compact_pair_transport &&
      (plans.size() != 1 ||
       required_string(source, "tile_plan") != plans.front()->handle() ||
       required_string(source, "tile_plan_checkpoint_identity") !=
           plans.front()->checkpoint_identity() ||
       (!compact_single_transport &&
        source.at("locals") != line_aggregate_source_records(local_owners))))
    throw std::invalid_argument(
        "checkpoint line aggregate source disagrees with its strong owners");
  if (compact_single_transport &&
      (transport_owners.size() != 1 ||
       required_string(source, "transport_state") !=
           transport_owners.front()->handle() ||
       required_string(source, "transport_state_checkpoint_identity") !=
           transport_owners.front()->checkpoint_identity() ||
       (!compact_v2 &&
        required_string(source, "transport_state_provenance_identity") !=
            transport_owners.front()->provenance_identity()) ||
       transport_owners.front()->plan_owner().get() != plans.front().get() ||
       transport_owners.front()->arm_name() !=
           required_string(provenance, "arm")))
    throw std::invalid_argument(
        "checkpoint compact transport line lost its exact retained state owner");
  if (compact_pair_transport) {
    if (plans.size() != 2 || transport_owners.size() != 2)
      throw std::invalid_argument(
          "checkpoint compact transport-pair line lost its ordered owners");
    require_transport_pair_compatibility(
        transport_owners[0], transport_owners[1],
        transport_owners[0]->anchor_owner()->scalar_domain());
    const auto validate_side = [&](const json::value& raw_side,
                                   const std::shared_ptr<StoredTilePlan>& plan,
                                   const std::shared_ptr<StoredTransportArmState>& state,
                                   const char* label) {
      const auto& side = as_object(raw_side, label);
      if (compact_v2)
        require_exact_keys(
            side,
            {"tile_plan", "tile_plan_checkpoint_identity",
             "transport_state"}, label);
      else
        require_exact_keys(
            side,
            {"tile_plan", "tile_plan_checkpoint_identity",
             "tile_plan_provenance_identity", "transport_state"},
            label);
      if (required_string(side, "tile_plan") != plan->handle() ||
          required_string(side, "tile_plan_checkpoint_identity") !=
              plan->checkpoint_identity() ||
          (!compact_v2 &&
           required_string(side, "tile_plan_provenance_identity") !=
               plan->provenance_identity()) ||
          state->plan_owner().get() != plan.get())
        throw std::invalid_argument(
            std::string(label) + " differs from its retained plan owner");
      const auto& state_reference = as_object(
          side.at("transport_state"), label);
      if (compact_v2)
        require_exact_keys(
            state_reference, {"handle", "checkpoint_identity"}, label);
      else
        require_exact_keys(
            state_reference,
            {"handle", "checkpoint_identity", "provenance_identity"},
            label);
      if (required_string(state_reference, "handle") != state->handle() ||
          required_string(state_reference, "checkpoint_identity") !=
              state->checkpoint_identity() ||
          (!compact_v2 &&
           required_string(state_reference, "provenance_identity") !=
               state->provenance_identity()))
        throw std::invalid_argument(
            std::string(label) + " differs from its retained state owner");
    };
    validate_side(source.at("lower"), plans[0], transport_owners[0],
                  "checkpoint transport-pair lower source");
    validate_side(source.at("upper"), plans[1], transport_owners[1],
                  "checkpoint transport-pair upper source");
    const auto& anchor = as_object(
        source.at("common_anchor"),
        "checkpoint transport-pair common anchor");
    require_exact_keys(
        anchor,
        {"local", "checkpoint_identity", "source_operator_identity",
         "physical_exact"},
        "checkpoint transport-pair common anchor");
    const auto& owner = transport_owners[0]->anchor_owner();
    if (required_string(anchor, "local") != owner->handle() ||
        required_string(anchor, "checkpoint_identity") !=
            owner->checkpoint_identity() ||
        required_string(anchor, "source_operator_identity") !=
            owner->source_operator_identity() ||
        required_string(anchor, "physical_exact") !=
            plans[0]->arm("lower").exact.from.str())
      throw std::invalid_argument(
          "checkpoint transport-pair common anchor is stale");
  } else if (compact_transport != (transport_owners.size() == 1)) {
    throw std::invalid_argument(
        "checkpoint compact transport ownership differs from its schema");
  }
  for (const auto& owner : local_owners) {
    if (std::string(owner->scalar_domain()) == "symbolic")
      throw std::invalid_argument(
          "checkpoint line aggregate cannot own a symbolic local");
    validate_checkpoint_exact_analytic_metadata(
        owner->exact_analytic_metadata());
  }
  const auto& aggregate = as_object(
      provenance.at("aggregate"), "checkpoint line aggregate recipe");
  std::optional<WholeArmEpsilonContract> whole_arm_epsilon_contract;
  std::optional<ObservableEpsilonContract> observable_epsilon_contract;
  std::optional<TransportTailPolicy> transport_tail_policy;
  bool has_aggregate_divergent_cancellation = false;
  std::string aggregate_divergent_cancellation_mode = "exact-singleton";
  std::optional<BoundedDivergentCancellation>
      aggregate_bounded_divergent_cancellation;
  if (!compact_transport) {
    const auto* raw_contract = aggregate.if_contains("epsilon_contract");
    if (raw_contract != nullptr)
      whole_arm_epsilon_contract = parse_whole_arm_epsilon_contract(
          *raw_contract, "checkpoint whole-arm epsilon contract");
  }
  const auto validate_observable_rows = [&](
      const json::value& raw_rows, const json::value& raw_tile_count,
      const std::shared_ptr<StoredTransportArmState>& owner,
      const char* label) {
    const auto& rows = as_array(raw_rows, label);
    if (!owner || rows.size() != owner->tile_sources().size() ||
        checkpoint_size_t(raw_tile_count, label) != rows.size())
      throw std::invalid_argument(
          std::string(label) +
          " does not reproduce its retained transport topology");
    for (std::size_t index = 0; index < rows.size(); ++index) {
      const auto& row = as_object(rows[index], label);
      const bool legacy_full_row =
          row.if_contains("prepared_row") != nullptr;
      if (!legacy_full_row && !compact_v2)
        throw std::invalid_argument(
            std::string(label) +
            " compact row identities require a SHA-bound v2 provenance");
      if (legacy_full_row)
        require_exact_keys(
            row, {"tile", "exact_identity", "prepared_row"}, label);
      else
        require_exact_keys(
            row, {"tile", "exact_identity", "entries"}, label);
      const auto exact_identity = required_string(row, "exact_identity");
      if (checkpoint_size_t(row.at("tile"), label) != index ||
          exact_identity.empty())
        throw std::invalid_argument(
            std::string(label) + " is out of order or stale");
      if (legacy_full_row) {
        const auto& prepared = as_object(row.at("prepared_row"), label);
        if (exact_identity != required_string(prepared, "exact_identity"))
          throw std::invalid_argument(
              std::string(label) + " legacy prepared row is stale");
      } else {
        std::optional<std::uint32_t> previous_column;
        for (const auto& raw_entry : as_array(row.at("entries"), label)) {
          const auto& entry = as_object(raw_entry, label);
          require_exact_keys(
              entry,
              {"column", "epsilon_shift", "center_pole_order",
               "exact_identity"},
              label);
          const auto column = as_u32(entry.at("column"), label);
          (void)as_i32(entry.at("epsilon_shift"), label);
          (void)as_u32(entry.at("center_pole_order"), label);
          if ((previous_column.has_value() && *previous_column >= column) ||
              required_string(entry, "exact_identity").empty())
            throw std::invalid_argument(
                std::string(label) +
                " has stale or unordered compact entry facts");
          previous_column = column;
        }
      }
    }
  };
  if (compact_single_transport) {
    if (compact_v2)
      require_exact_keys(
          aggregate,
          {"kind", "combination", "request_index", "observable_identity",
           "observable_checkpoint_identity", "transport_state",
           "output_epsilon_contract", "tail_policy", "projection_mode",
           "rows", "tile_count"},
          "checkpoint compact transport observable recipe");
    else
      require_exact_keys(
          aggregate,
          {"kind", "combination", "request_index", "observable_identity",
           "observable_checkpoint_identity", "transport_state",
           "output_epsilon_contract", "tail_policy", "rows", "tile_count"},
          "checkpoint legacy compact transport observable recipe");
    if (required_string(aggregate, "kind") !=
            "transport-state-observable-arm" ||
        required_string(aggregate, "combination") !=
            "sum-physical-tiles" ||
        required_string(aggregate, "observable_identity").empty() ||
        required_string(aggregate, "observable_checkpoint_identity") !=
            checkpoint_identity)
      throw std::invalid_argument(
          "checkpoint compact transport observable identity is inconsistent");
    (void)checkpoint_size_t(aggregate.at("request_index"),
                            "checkpoint observable request index");
    const auto& state_reference = as_object(
        aggregate.at("transport_state"),
        "checkpoint observable transport state reference");
    if (compact_v2)
      require_exact_keys(
          state_reference, {"handle", "checkpoint_identity"},
          "checkpoint observable transport state reference");
    else
      require_exact_keys(state_reference,
                         {"handle", "checkpoint_identity",
                          "provenance_identity"},
                         "checkpoint observable transport state reference");
    const auto& transport_owner = transport_owners.front();
    if (!transport_owner ||
        required_string(state_reference, "handle") !=
            transport_owner->handle() ||
        required_string(state_reference, "checkpoint_identity") !=
            transport_owner->checkpoint_identity() ||
        (!compact_v2 &&
         required_string(state_reference, "provenance_identity") !=
             transport_owner->provenance_identity()))
      throw std::invalid_argument(
          "checkpoint observable recipe names a different transport state");
    observable_epsilon_contract = parse_observable_epsilon_contract(
        aggregate.at("output_epsilon_contract"),
        "checkpoint observable output epsilon contract");
    transport_tail_policy = parse_transport_tail_policy(
        aggregate.at("tail_policy"),
        "checkpoint observable tail policy");
    if (compact_v2)
      require_transport_projection_mode(
          aggregate, *transport_tail_policy,
          "checkpoint compact transport observable recipe");
    validate_observable_rows(
        aggregate.at("rows"), aggregate.at("tile_count"), transport_owner,
        "checkpoint observable prepared rows");
  } else if (compact_pair_transport) {
    has_aggregate_divergent_cancellation =
        aggregate.if_contains("divergent_cancellation") != nullptr;
    if (has_aggregate_divergent_cancellation) {
      if (compact_v2)
        require_exact_keys(
            aggregate,
            {"kind", "combination", "request_index", "observable_identity",
             "observable_checkpoint_identity", "output_epsilon_contract",
             "tail_policy", "projection_mode", "divergent_cancellation",
             "no_remarching", "no_rematching", "concurrent_arms", "lower",
             "upper"},
            "checkpoint compact transport-pair observable recipe");
      else
        require_exact_keys(
            aggregate,
            {"kind", "combination", "request_index", "observable_identity",
             "observable_checkpoint_identity", "output_epsilon_contract",
             "tail_policy", "divergent_cancellation", "no_remarching",
             "no_rematching", "concurrent_arms", "lower", "upper"},
            "checkpoint legacy compact transport-pair observable recipe");
    } else {
      if (compact_v2)
        require_exact_keys(
            aggregate,
            {"kind", "combination", "request_index", "observable_identity",
             "observable_checkpoint_identity", "output_epsilon_contract",
             "tail_policy", "projection_mode", "no_remarching",
             "no_rematching", "concurrent_arms", "lower", "upper"},
            "checkpoint compact transport-pair observable recipe");
      else
        require_exact_keys(
            aggregate,
            {"kind", "combination", "request_index", "observable_identity",
             "observable_checkpoint_identity", "output_epsilon_contract",
             "tail_policy", "no_remarching", "no_rematching",
             "concurrent_arms", "lower", "upper"},
            "checkpoint legacy compact transport-pair observable recipe");
    }
    if (required_string(aggregate, "kind") !=
            "transport-state-observable-pair" ||
        required_string(aggregate, "combination") !=
            "negative-lower-plus-upper" ||
        required_string(aggregate, "observable_identity").empty() ||
        required_string(aggregate, "observable_checkpoint_identity") !=
            checkpoint_identity ||
        required_string(provenance, "arm") != "combined" ||
        !aggregate.at("no_remarching").is_bool() ||
        !aggregate.at("no_rematching").is_bool() ||
        !aggregate.at("concurrent_arms").is_bool() ||
        !aggregate.at("no_remarching").as_bool() ||
        !aggregate.at("no_rematching").as_bool() ||
        !aggregate.at("concurrent_arms").as_bool())
      throw std::invalid_argument(
          "checkpoint compact transport-pair observable identity or fixed combination is inconsistent");
    (void)checkpoint_size_t(
        aggregate.at("request_index"),
        "checkpoint transport-pair observable request index");
    observable_epsilon_contract = parse_observable_epsilon_contract(
        aggregate.at("output_epsilon_contract"),
        "checkpoint transport-pair output epsilon contract");
    transport_tail_policy = parse_transport_tail_policy(
        aggregate.at("tail_policy"),
        "checkpoint transport-pair tail policy");
    if (compact_v2)
      require_transport_projection_mode(
          aggregate, *transport_tail_policy,
          "checkpoint compact transport-pair observable recipe");
    if (has_aggregate_divergent_cancellation) {
      const auto& cancellation = as_object(
          aggregate.at("divergent_cancellation"),
          "checkpoint transport-pair divergent-cancellation policy");
      aggregate_divergent_cancellation_mode = required_string(
          cancellation, "mode");
      if (aggregate_divergent_cancellation_mode == "exact-singleton") {
        require_exact_keys(
            cancellation, {"mode"},
            "checkpoint exact-singleton divergent-cancellation policy");
      } else if (aggregate_divergent_cancellation_mode ==
                 "bounded-relative-acb") {
        aggregate_bounded_divergent_cancellation =
            parse_bounded_divergent_cancellation(
                cancellation,
                "checkpoint bounded divergent-cancellation policy");
      } else {
        throw std::invalid_argument(
            "checkpoint transport-pair has an unsupported divergent-cancellation mode");
      }
    }
    const auto validate_side_recipe = [&](
        const json::value& raw_side,
        const std::shared_ptr<StoredTransportArmState>& owner,
        const char* label) {
      const auto& side = as_object(raw_side, label);
      require_exact_keys(
          side, {"transport_state", "rows", "tile_count"}, label);
      const auto& reference = as_object(
          side.at("transport_state"), label);
      if (compact_v2)
        require_exact_keys(reference, {"handle", "checkpoint_identity"},
                           label);
      else
        require_exact_keys(reference,
                           {"handle", "checkpoint_identity",
                            "provenance_identity"}, label);
      if (required_string(reference, "handle") != owner->handle() ||
          required_string(reference, "checkpoint_identity") !=
              owner->checkpoint_identity() ||
          (!compact_v2 &&
           required_string(reference, "provenance_identity") !=
               owner->provenance_identity()))
        throw std::invalid_argument(
            std::string(label) + " names a different transport state");
      validate_observable_rows(side.at("rows"), side.at("tile_count"),
                               owner, label);
    };
    validate_side_recipe(aggregate.at("lower"), transport_owners[0],
                         "checkpoint transport-pair lower recipe");
    validate_side_recipe(aggregate.at("upper"), transport_owners[1],
                         "checkpoint transport-pair upper recipe");
  } else {
    const auto& components = as_array(
        aggregate.at("components"), "checkpoint line aggregate components");
    if (checkpoint_size_t(aggregate.at("component_count"),
                          "checkpoint line aggregate component count") !=
            components.size() ||
        components.empty())
      throw std::invalid_argument(
          "checkpoint line aggregate component manifest is inconsistent");
    for (std::size_t index = 0; index < components.size(); ++index) {
      const auto& component = as_object(
          components[index], "checkpoint line aggregate component");
      require_exact_keys(
          component,
          {"index", "sign", "checkpoint_identity", "provenance_identity",
           "scope", "source", "interval", "epsilon", "error"},
          "checkpoint line aggregate component");
      const auto sign = as_i32(component.at("sign"),
                               "checkpoint line aggregate sign");
      if (checkpoint_size_t(component.at("index"),
                            "checkpoint line aggregate component index") !=
              index ||
          (sign != -1 && sign != 1) ||
          required_string(component, "checkpoint_identity").empty() ||
          required_string(component, "provenance_identity").empty())
        throw std::invalid_argument(
            "checkpoint line aggregate component provenance is malformed");
    }
  }

  const auto& raw_result = as_object(
      object.at("result"), "checkpoint line aggregate result");
  require_exact_keys(raw_result,
      {"value", "scope", "imaginary_sign", "diagnostics"},
      "checkpoint line aggregate result");
  StoredLineIntegral result;
  result.value = parse_checkpoint_epsilon_vector(
      raw_result.at("value"), "checkpoint line aggregate epsilon vector");
  const auto scope = required_string(raw_result, "scope");
  if (scope == "stored_truncation")
    result.scope = LineIntegrationScope::StoredTruncation;
  else if (scope == "full_local_with_certified_tail")
    result.scope = LineIntegrationScope::FullLocalWithCertifiedTail;
  else
    throw std::invalid_argument(
        "checkpoint line aggregate has an unsupported integration scope");
  if (!raw_result.at("imaginary_sign").is_null())
    throw std::invalid_argument(
        "checkpoint multi-chart line aggregate cannot carry one effective rim");
  result.imaginary_sign = std::nullopt;
  const auto& epsilon = as_object(
      provenance.at("epsilon"), "checkpoint line aggregate epsilon");
  require_exact_keys(epsilon, {"min", "max"},
                     "checkpoint line aggregate epsilon");
  if (result.value.epsilon.min_power !=
          as_i32(epsilon.at("min"), "checkpoint aggregate epsilon minimum") ||
      result.value.epsilon.complete_max !=
          as_i32(epsilon.at("max"), "checkpoint aggregate epsilon maximum") ||
      required_string(provenance, "scope") != scope ||
      required_string(provenance, "error_guarantee") !=
          error_guarantee_name(result.value.error.guarantee))
    throw std::invalid_argument(
        "checkpoint line aggregate result differs from its provenance");
  if (whole_arm_epsilon_contract.has_value() &&
      result.value.epsilon.complete_max <
          whole_arm_epsilon_contract->public_required_complete_max)
    throw std::invalid_argument(
        "checkpoint whole-arm aggregate no longer covers its public required epsilon maximum");
  if (observable_epsilon_contract.has_value() &&
      result.value.epsilon.complete_max <
          observable_epsilon_contract->required_complete_max)
    throw std::invalid_argument(
        "checkpoint observable aggregate no longer covers its required epsilon maximum");
  if (transport_tail_policy == TransportTailPolicy::Require &&
      result.scope != LineIntegrationScope::FullLocalWithCertifiedTail)
    throw std::invalid_argument(
        "checkpoint required-tail observable is not certified full-local");
  const auto expected_dimension = compact_transport
      ? 1u
      : as_u32(local_owners.front()->summary().at("dimension"),
               "checkpoint line aggregate source dimension");
  if (result.value.dimension != expected_dimension ||
      (!compact_transport &&
       std::any_of(local_owners.begin(), local_owners.end(),
                   [&](const auto& owner) {
                     return as_u32(
                                owner->summary().at("dimension"),
                                "checkpoint aggregate owner dimension") !=
                            expected_dimension;
                   })))
    throw std::invalid_argument(
        "checkpoint line aggregate dimensions disagree with its owners");
  if (result.scope == LineIntegrationScope::FullLocalWithCertifiedTail &&
      result.value.error.guarantee != ErrorGuarantee::Certified)
    throw std::invalid_argument(
        "checkpoint certified line aggregate lost its certified error envelope");

  const auto& diagnostics = as_object(
      raw_result.at("diagnostics"),
      "checkpoint line aggregate diagnostics");
  const std::array<const char*, 4> cancellation_diagnostic_keys{
      "bounded_cancelled_divergent_coefficients",
      "divergent_cancellation_mode", "divergent_relative_tolerance",
      "divergent_cancellation_provenance"};
  const auto cancellation_diagnostic_count = std::count_if(
      cancellation_diagnostic_keys.begin(),
      cancellation_diagnostic_keys.end(), [&](const char* key) {
        return diagnostics.if_contains(key) != nullptr;
      });
  const bool has_cancellation_diagnostics =
      cancellation_diagnostic_count == cancellation_diagnostic_keys.size();
  if (cancellation_diagnostic_count != 0 &&
      !has_cancellation_diagnostics)
    throw std::invalid_argument(
        "checkpoint aggregate divergent-cancellation diagnostics are incomplete");
  if (has_cancellation_diagnostics)
    require_exact_keys(
        diagnostics,
        {"input_monomial_cells", "grouped_monomials", "zero_groups_skipped",
         "cancelled_divergent_groups",
         "bounded_cancelled_divergent_coefficients",
         "divergent_cancellation_mode", "divergent_relative_tolerance",
         "divergent_cancellation_provenance", "primitive_evaluations",
         "primitive_component_applications", "primitive_component_reuses",
         "has_center_endpoint", "tail_certificate_requested",
         "tail_certificate_status", "tail_witness_radius_exact", "detail"},
        "checkpoint line aggregate diagnostics");
  else
    require_exact_keys(
        diagnostics,
        {"input_monomial_cells", "grouped_monomials", "zero_groups_skipped",
         "cancelled_divergent_groups", "primitive_evaluations",
         "primitive_component_applications", "primitive_component_reuses",
         "has_center_endpoint", "tail_certificate_requested",
         "tail_certificate_status", "tail_witness_radius_exact", "detail"},
        "checkpoint legacy line aggregate diagnostics");
  result.diagnostics.input_monomial_cells = checkpoint_size_t(
      diagnostics.at("input_monomial_cells"),
      "checkpoint aggregate input monomials");
  result.diagnostics.grouped_monomials = checkpoint_size_t(
      diagnostics.at("grouped_monomials"),
      "checkpoint aggregate grouped monomials");
  result.diagnostics.zero_groups_skipped = checkpoint_size_t(
      diagnostics.at("zero_groups_skipped"),
      "checkpoint aggregate skipped zero groups");
  result.diagnostics.cancelled_divergent_groups = checkpoint_size_t(
      diagnostics.at("cancelled_divergent_groups"),
      "checkpoint aggregate cancelled divergent groups");
  if (has_cancellation_diagnostics) {
    result.diagnostics.bounded_cancelled_divergent_coefficients =
        checkpoint_size_t(
            diagnostics.at("bounded_cancelled_divergent_coefficients"),
            "checkpoint aggregate bounded-cancelled divergent coefficients");
    result.diagnostics.divergent_cancellation_mode = required_string(
        diagnostics, "divergent_cancellation_mode");
    if (!diagnostics.at("divergent_relative_tolerance").is_string() ||
        !diagnostics.at("divergent_cancellation_provenance").is_string())
      throw std::invalid_argument(
          "checkpoint aggregate divergent-cancellation strings are malformed");
    result.diagnostics.divergent_relative_tolerance = std::string(
        diagnostics.at("divergent_relative_tolerance").as_string());
    result.diagnostics.divergent_cancellation_provenance = std::string(
        diagnostics.at("divergent_cancellation_provenance").as_string());
    if (result.diagnostics.divergent_cancellation_mode != "exact-singleton" &&
        result.diagnostics.divergent_cancellation_mode !=
            "bounded-relative-acb")
      throw std::invalid_argument(
          "checkpoint aggregate has an unsupported divergent-cancellation mode");
    const bool bounded_cancellation =
        result.diagnostics.divergent_cancellation_mode ==
        "bounded-relative-acb";
    const bool has_cancellation_tolerance =
        !result.diagnostics.divergent_relative_tolerance.empty();
    const bool has_cancellation_provenance =
        !result.diagnostics.divergent_cancellation_provenance.empty();
    if ((bounded_cancellation &&
         (!has_cancellation_tolerance || !has_cancellation_provenance)) ||
        (!bounded_cancellation &&
         (has_cancellation_tolerance || has_cancellation_provenance)))
      throw std::invalid_argument(
          "checkpoint aggregate divergent-cancellation policy is incomplete");
    if (!bounded_cancellation &&
        result.diagnostics.bounded_cancelled_divergent_coefficients != 0)
      throw std::invalid_argument(
          "checkpoint exact-singleton aggregate cannot report "
          "bounded-cancelled divergent coefficients");
  }
  result.diagnostics.primitive_evaluations = checkpoint_size_t(
      diagnostics.at("primitive_evaluations"),
      "checkpoint aggregate primitive evaluations");
  result.diagnostics.primitive_component_applications = checkpoint_size_t(
      diagnostics.at("primitive_component_applications"),
      "checkpoint aggregate primitive applications");
  result.diagnostics.primitive_component_reuses = checkpoint_size_t(
      diagnostics.at("primitive_component_reuses"),
      "checkpoint aggregate primitive reuses");
  if (!diagnostics.at("has_center_endpoint").is_bool() ||
      !diagnostics.at("tail_certificate_requested").is_bool())
    throw std::invalid_argument(
        "checkpoint line aggregate diagnostic flags must be Boolean");
  result.diagnostics.has_center_endpoint =
      diagnostics.at("has_center_endpoint").as_bool();
  result.diagnostics.tail_certificate_requested =
      diagnostics.at("tail_certificate_requested").as_bool();
  result.diagnostics.tail_certificate_status = required_string(
      diagnostics, "tail_certificate_status");
  if (!diagnostics.at("tail_witness_radius_exact").is_string())
    throw std::invalid_argument(
        "checkpoint aggregate tail witness radius must be a string");
  result.diagnostics.tail_witness_radius_exact = std::string(
      diagnostics.at("tail_witness_radius_exact").as_string());
  result.diagnostics.detail = required_string(diagnostics, "detail");

  if (compact_pair_transport) {
    if (!has_aggregate_divergent_cancellation) {
      if (result.diagnostics.divergent_cancellation_mode !=
          "exact-singleton")
        throw std::invalid_argument(
            "legacy checkpoint pair lacks the bounded-cancellation provenance named by its diagnostics");
    } else if (result.diagnostics.divergent_cancellation_mode !=
                   aggregate_divergent_cancellation_mode ||
               (aggregate_bounded_divergent_cancellation.has_value() &&
                (result.diagnostics.divergent_relative_tolerance !=
                     aggregate_bounded_divergent_cancellation
                         ->relative_tolerance_text ||
                 result.diagnostics.divergent_cancellation_provenance !=
                     aggregate_bounded_divergent_cancellation
                         ->provenance))) {
      throw std::invalid_argument(
          "checkpoint pair divergent-cancellation diagnostics differ from aggregate provenance");
    }
  }

  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint line aggregate elapsed time");
  const auto& stats = as_object(object.at("runtime_stats"),
                                "checkpoint line aggregate runtime stats");
  require_exact_keys(stats, {"exports", "export_ms"},
                     "checkpoint line aggregate runtime stats");
  auto stored = std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, std::move(result),
      elapsed_ms, std::move(plans), std::move(local_owners), provenance,
      std::move(transport_owners));
  stored->restore_runtime_stats(
      as_u64(stats.at("exports"), "checkpoint aggregate exports"),
      checkpoint_nonnegative_double(stats.at("export_ms"),
                                    "checkpoint aggregate export time"));
  return stored;
}

json::object normalized_checkpoint_line_record_for_roundtrip(
    const json::value& raw) {
  auto normalized = as_object(
      raw, "checkpoint retained line result");
  auto& diagnostics = normalized.at("result").as_object()
      .at("diagnostics").as_object();
  const std::array<const char*, 4> cancellation_diagnostic_keys{
      "bounded_cancelled_divergent_coefficients",
      "divergent_cancellation_mode", "divergent_relative_tolerance",
      "divergent_cancellation_provenance"};
  const auto present = std::count_if(
      cancellation_diagnostic_keys.begin(),
      cancellation_diagnostic_keys.end(), [&](const char* key) {
        return diagnostics.if_contains(key) != nullptr;
      });
  // Checkpoints written before divergent-cancellation diagnostics were
  // retained encode the same exact-singleton policy by omitting all four
  // fields.  The restore routines validate that the group is either wholly
  // present or wholly absent; normalize only the legacy all-absent form so
  // the exact reconstructed-state comparison remains meaningful.
  if (present == 0) {
    diagnostics["bounded_cancelled_divergent_coefficients"] = 0;
    diagnostics["divergent_cancellation_mode"] = "exact-singleton";
    diagnostics["divergent_relative_tolerance"] = "";
    diagnostics["divergent_cancellation_provenance"] = "";
  }
  const std::array<const char*, 3> tail_diagnostic_keys{
      "tail_certificate_requested", "tail_certificate_status",
      "tail_witness_radius_exact"};
  const auto tail_present = std::count_if(
      tail_diagnostic_keys.begin(), tail_diagnostic_keys.end(),
      [&](const char* key) {
        return diagnostics.if_contains(key) != nullptr;
      });
  // The older stored-truncation record predates retained tail diagnostics.
  // Its decoded state is the current explicit not-requested state.
  if (tail_present == 0) {
    diagnostics["tail_certificate_requested"] = false;
    diagnostics["tail_certificate_status"] = "not-requested";
    diagnostics["tail_witness_radius_exact"] = nullptr;
  }
  return normalized;
}

json::array encode_strings(const std::vector<std::string>& values) {
  json::array output;
  output.reserve(values.size());
  for (const auto& value : values) output.emplace_back(value);
  return output;
}

json::object checkpoint_configuration_record(const SolverSession& session) {
  return json::object{
      {"domain", session.domain},
      {"precision_bits", session.precision_bits},
      {"output_digits", session.output_digits},
      {"symbols", encode_strings(session.symbols)},
      {"analytic", json::parse(session.analytic_identity)},
      {"chart_capacity", session.chart_capacity},
      {"local_capacity", session.local_capacity},
      {"scc_capacity", session.scc_capacity},
      {"match_capacity", session.match_capacity},
      {"endpoint_capacity", session.endpoint_capacity},
      {"transport_state_capacity", session.transport_state_capacity}};
}

std::string checkpoint_configuration_identity(const SolverSession& session) {
  return json::serialize(canonical_json_value(
      checkpoint_configuration_record(session)));
}

json::object checkpoint_chart_item(
    const SolverSession& session,
    const std::shared_ptr<PreparedChartBase>& chart) {
  const auto signature_value = json::parse(chart->signature());
  const auto& signature = as_object(
      signature_value, "prepared chart exact signature");
  const auto& analytic = as_object(
      signature.at("analytic"), "prepared chart analytic signature");
  if (json::serialize(analytic.at("session")) != session.analytic_identity)
    throw std::logic_error(
        "prepared chart session analytic identity changed after retention");
  if (required_string(signature, "identity") != chart->exact_identity())
    throw std::logic_error(
        "prepared chart exact identity changed after retention");
  if (const auto* physical = signature.if_contains("physical_ode")) {
    if (chart->physical_payload_record().empty() ||
        json::serialize(canonical_json_value(*physical)) !=
            chart->physical_payload_record())
      throw std::logic_error(
          "prepared chart physical q/C payload changed after retention");
  } else if (!chart->physical_payload_record().empty()) {
    throw std::logic_error(
        "prepared chart signature lost its retained physical q/C payload");
  }

  json::object problem = signature;
  problem.erase("identity");
  problem.erase("analytic");
  problem.erase("scc");
  json::object request{
      {"schema", 2}, {"op", "chart.prepare"},
      {"session", session.handle}, {"key", chart->key()},
      {"identity", chart->exact_identity()},
      {"analytic", analytic.at("chart")},
      {"scc", signature.at("scc")}, {"problem", std::move(problem)}};
  return json::object{{"handle", chart->handle()},
                      {"key", chart->key()},
                      {"identity", chart->exact_identity()},
                      {"signature", chart->signature()},
                      {"request", std::move(request)}};
}

json::object checkpoint_scc_item(
    const SolverSession& session,
    const std::shared_ptr<CompositeSCCChartBase>& composite) {
  const auto signature_value = json::parse(composite->signature());
  const auto& signature = as_object(
      signature_value, "retained SCC exact signature");
  if (required_string(signature, "identity") != composite->exact_identity())
    throw std::logic_error(
        "retained SCC exact identity changed after retention");
  if (const auto* physical = signature.if_contains("physical_ode")) {
    if (composite->physical_payload_record().empty() ||
        json::serialize(canonical_json_value(*physical)) !=
            composite->physical_payload_record() ||
        composite->owner_signature_identity() !=
            composite->exact_identity())
      throw std::logic_error(
          "retained SCC full-parent physical q/C payload changed after retention");
  } else if (!composite->physical_payload_record().empty()) {
    throw std::logic_error(
        "retained SCC signature lost its full-parent physical q/C payload");
  }
  json::object request = signature;
  request["schema"] = 2;
  request["op"] = "scc.prepare";
  request["session"] = session.handle;
  request["key"] = composite->key();
  return json::object{{"handle", composite->handle()},
                      {"key", composite->key()},
                      {"identity", composite->exact_identity()},
                      {"signature", composite->signature()},
                      {"request", std::move(request)}};
}

json::array checkpoint_identity_manifest(const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(value, "checkpoint retained item");
    manifest.push_back(json::object{{"handle", item.at("handle")},
                                    {"key", item.at("key")},
                                    {"identity", item.at("identity")}});
  }
  return manifest;
}

json::array checkpoint_local_identity_manifest(const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(value, "checkpoint retained local");
    const auto& solution = as_object(item.at("solution"),
                                     "checkpoint local solution");
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"source_chart", item.at("source_chart")},
        {"source_operator_identity", item.at("source_operator_identity")},
        {"scalar_domain", item.at("scalar_domain")},
        {"checkpoint_identity", solution.at("checkpoint_identity")}});
  }
  return manifest;
}

json::array checkpoint_acb_match_identity_manifest(const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(value, "checkpoint retained Acb match");
    json::object identity{
        {"handle", item.at("handle")},
        {"checkpoint_identity", item.at("checkpoint_identity")},
        {"provenance_identity", item.at("provenance_identity")},
        {"exact_lattice_identity", item.at("exact_lattice_identity")},
        {"matching_frame_identity",
         item.at("matching_frame_identity")}};
    if (const auto* residual_frame =
            item.if_contains("residual_frame_identity"))
      identity["residual_frame_identity"] = *residual_frame;
    manifest.push_back(std::move(identity));
  }
  return manifest;
}

json::array checkpoint_exact_match_identity_manifest(
    const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(
        value, "checkpoint retained exact-rational match");
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"checkpoint_identity", item.at("checkpoint_identity")},
        {"provenance_identity", item.at("provenance_identity")},
        {"basis_sources", item.at("basis_sources")},
        {"incoming_source", item.at("incoming_source")}});
  }
  return manifest;
}

json::array checkpoint_planned_match_identity_manifest(
    const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(
        value, "checkpoint retained planned match hop");
    const auto& handoff = as_object(
        item.at("handoff"), "checkpoint planned match handoff");
    const auto& native_match = as_object(
        item.at("native_match"), "checkpoint planned embedded match");
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"checkpoint_identity", item.at("checkpoint_identity")},
        {"provenance_identity", item.at("provenance_identity")},
        {"tile_plan", handoff.at("tile_plan")},
        {"tile_plan_checkpoint_identity",
         handoff.at("tile_plan_checkpoint_identity")},
        {"tile_plan_provenance_identity",
         handoff.at("tile_plan_provenance_identity")},
        {"native_match_schema", native_match.at("schema")},
        {"native_match_checkpoint_identity",
         native_match.at("checkpoint_identity")},
        {"native_match_provenance_identity",
         native_match.at("provenance_identity")}});
  }
  return manifest;
}

json::array checkpoint_endpoint_identity_manifest(const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(value, "checkpoint retained endpoint");
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"checkpoint_identity", item.at("checkpoint_identity")},
        {"provenance_identity", item.at("provenance_identity")},
        {"source", item.at("source")}});
  }
  return manifest;
}

json::array checkpoint_tile_identity_manifest(const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(value, "checkpoint retained tile plan");
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"checkpoint_identity", item.at("checkpoint_identity")},
        {"provenance_identity", item.at("provenance_identity")}});
  }
  return manifest;
}

json::array checkpoint_transport_state_identity_manifest(
    const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(
        value, "checkpoint retained transport-arm state");
    const auto& provenance = as_object(
        item.at("provenance"),
        "checkpoint transport-arm state provenance");
    json::value anchor;
    json::value final_local;
    if (required_string(provenance, "schema") ==
        "diffexp2-retained-native-transport-arm-state-v4") {
      const auto& chain = as_array(
          provenance.at("tile_checkpoint_chain"),
          "checkpoint certificate-only transport tile chain");
      if (chain.empty())
        throw std::logic_error(
            "checkpoint certificate-only transport state has no tile chain");
      anchor = chain.front();
      final_local = chain.back();
    } else {
      anchor = provenance.at("anchor");
      final_local = provenance.at("final_local");
    }
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"checkpoint_identity", item.at("checkpoint_identity")},
        {"provenance_identity", item.at("provenance_identity")},
        {"tile_plan", provenance.at("tile_plan")},
        {"arm", provenance.at("arm")},
        {"anchor", std::move(anchor)},
        {"final_local", std::move(final_local)}});
  }
  return manifest;
}

json::array checkpoint_line_identity_manifest(const json::array& items) {
  json::array manifest;
  manifest.reserve(items.size());
  for (const auto& value : items) {
    const auto& item = as_object(value, "checkpoint retained line result");
    const auto schema = required_string(item, "schema");
    const auto& source = (schema == "diffexp2-retained-line-aggregate-v1" ||
                          schema ==
                              "diffexp2-retained-transport-observable-line-v1" ||
                          schema ==
                              "diffexp2-retained-transport-observable-line-v2" ||
                          schema ==
                              "diffexp2-retained-transport-pair-observable-line-v1" ||
                          schema ==
                              "diffexp2-retained-transport-pair-observable-line-v2")
        ? as_object(item.at("provenance"),
                    "checkpoint line aggregate provenance").at("source")
        : item.at("source");
    manifest.push_back(json::object{
        {"handle", item.at("handle")},
        {"checkpoint_identity", item.at("checkpoint_identity")},
        {"provenance_identity", item.at("provenance_identity")},
        {"source", source}});
  }
  return manifest;
}
