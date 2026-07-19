template <typename Scalar>
std::shared_ptr<PreparedChartBase> parse_prepared_chart(
    const std::shared_ptr<SolverSession>& session, const json::object& root,
    const std::string& handle, const std::string& key,
    const std::string& exact_identity,
    std::optional<std::string> geometry_record,
    std::optional<std::string> principal_matrix_record,
    std::optional<std::string> native_scc_capabilities,
    std::optional<std::string> regular_value_relative_accuracy_max_exact,
    SCCCertificate scc,
    std::string signature) {
  const auto started = std::chrono::steady_clock::now();
  const auto& problem = as_object(root.at("problem"), "prepared problem");
  std::unique_ptr<AcbPrecisionLease> acb_lease;
  if constexpr (std::is_same_v<Scalar, ComplexBall>) {
    acb_lease = std::make_unique<AcbPrecisionLease>(session->precision_bits);
    ComplexBall::set_precision(session->precision_bits);
  }
  std::unique_lock<std::recursive_mutex> symbolic_lock;
  if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
    symbolic_lock =
        std::unique_lock<std::recursive_mutex>(symbolic_run_mutex());
    SymbolicRational::configure(session->symbols);
  }
  auto prepared = parse_prepared_operator<Scalar>(problem);
  std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>> physical_equation;
  if (const auto* raw_physical = problem.if_contains("physical_ode"))
    physical_equation = parse_prepared_physical_ode<Scalar>(
        *raw_physical, prepared.dimension);
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return make_retained_typed_shared<Scalar, PreparedChart<Scalar>>(
      handle, key, exact_identity, std::move(signature),
      std::move(geometry_record), std::move(principal_matrix_record),
      std::move(native_scc_capabilities),
      std::move(regular_value_relative_accuracy_max_exact), std::move(scc),
      std::move(prepared),
      std::move(physical_equation), session->precision_bits,
      session->symbols, elapsed);
}

json::object run_session_command(const json::object& root) {
  if (as_i64(root.at("schema"), "schema") != 2)
    throw std::invalid_argument("unsupported persistent solver schema");
  const auto operation = required_string(root, "op");

  if (operation == "checkpoint.restore") {
    require_exact_keys(root,
        {"schema", "op", "path", "expected_identity"},
        "checkpoint.restore request");
    return restore_checkpoint(required_string(root, "path"),
                              required_string(root, "expected_identity"));
  }

  if (operation == "session.create") {
    const auto domain = required_string(root, "domain");
    if (domain != "rational" && domain != "acb" && domain != "symbolic")
      throw std::invalid_argument("unsupported persistent scalar domain: " + domain);
    auto session = std::make_shared<SolverSession>();
    session->domain = domain;
    session->precision_bits = root.if_contains("precision_bits")
        ? static_cast<slong>(as_i64(root.at("precision_bits"), "precision bits"))
        : 256;
    if (domain == "acb" && session->precision_bits < 64)
      throw std::invalid_argument("Acb precision must be at least 64 bits");
    session->output_digits = root.if_contains("output_digits")
        ? static_cast<int>(as_i64(root.at("output_digits"), "output digits"))
        : 50;
    if (session->output_digits < 1)
      throw std::invalid_argument("output digits must be positive");
    session->symbols = parse_symbols(root);
    if (domain == "symbolic" && session->symbols.empty())
      throw std::invalid_argument(
          "symbolic persistent session requires declared regulator symbols");
    if (domain != "symbolic" && !session->symbols.empty())
      throw std::invalid_argument(
          "regulator symbols are only valid for the symbolic scalar domain");
    session->analytic_identity = root.if_contains("analytic")
        ? json::serialize(root.at("analytic")) : "null";
    if (root.if_contains("chart_capacity")) {
      const auto capacity = as_u32(root.at("chart_capacity"), "chart capacity");
      if (capacity == 0 || capacity > 4096)
        throw std::invalid_argument("chart capacity must lie in 1..4096");
      session->chart_capacity = capacity;
    }
    if (root.if_contains("local_capacity")) {
      const auto capacity = as_u32(root.at("local_capacity"), "local capacity");
      if (capacity == 0 || capacity > 16384)
        throw std::invalid_argument("local capacity must lie in 1..16384");
      session->local_capacity = capacity;
    }
    if (root.if_contains("scc_capacity")) {
      const auto capacity = as_u32(root.at("scc_capacity"), "SCC capacity");
      if (capacity == 0 || capacity > 4096)
        throw std::invalid_argument("SCC capacity must lie in 1..4096");
      session->scc_capacity = capacity;
    }
    if (root.if_contains("match_capacity")) {
      const auto capacity = as_u32(root.at("match_capacity"),
                                   "local match capacity");
      if (capacity == 0 || capacity > 16384)
        throw std::invalid_argument(
            "local match capacity must lie in 1..16384");
      session->match_capacity = capacity;
    }
    if (root.if_contains("endpoint_capacity")) {
      const auto capacity = as_u32(root.at("endpoint_capacity"),
                                   "endpoint result capacity");
      if (capacity == 0 || capacity > 16384)
        throw std::invalid_argument(
            "endpoint result capacity must lie in 1..16384");
      session->endpoint_capacity = capacity;
    }
    if (root.if_contains("transport_state_capacity")) {
      const auto capacity = as_u32(
          root.at("transport_state_capacity"),
          "transport-arm state capacity");
      if (capacity == 0 || capacity > 16384)
        throw std::invalid_argument(
            "transport-arm state capacity must lie in 1..16384");
      session->transport_state_capacity = capacity;
    }
    auto& registry = session_registry();
    {
      std::lock_guard<std::mutex> lock(registry.mutex);
      // Precision is thread-local, and every Acb solve/evaluation acquires a
      // lease and installs its session precision.  Live sessions may retain
      // different precisions; unequal active operations serialize as needed.
      session->handle = "s:" + std::to_string(registry.next_session++);
      registry.sessions.emplace(session->handle, session);
    }
    return json::object{{"status", "ok"}, {"session", session->handle},
                        {"domain", session->domain},
                        {"precision_bits", session->precision_bits},
                        {"chart_capacity", session->chart_capacity},
                        {"local_capacity", session->local_capacity},
                        {"scc_capacity", session->scc_capacity},
                        {"match_capacity", session->match_capacity},
                        {"endpoint_capacity", session->endpoint_capacity},
                        {"tile_plan_capacity", session->tile_plan_capacity},
                        {"transport_state_capacity",
                         session->transport_state_capacity},
                        {"line_result_capacity", session->line_result_capacity},
                        {"regular_equation_owner_capability",
                         kFrameIndependentRegularEquationOwnerCapability},
                        {"local_match_capability",
                         domain == "rational"
                             ? kExactRegularLocalMatchCapability
                             : "unsupported"},
                        {"acb_local_match_capability",
                         domain == "acb"
                             ? kRefinedAcbLocalMatchCapability
                             : "unsupported"},
                        {"endpoint_limit_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedEndpointLimitCapability},
                        {"planned_endpoint_limit_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedPlannedEndpointLimitCapability},
                        {"tile_plan_capability", kRetainedTilePlanCapability},
                        {"single_arm_tile_plan_capability",
                         kRetainedSingleArmTilePlanCapability},
                        {"planned_match_hop_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedPlannedMatchHopCapability},
                        {"planned_match_materialization_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedPlannedMatchMaterializationCapability},
                        {"rational_row_application_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedRationalRowCapability},
                        {"line_integration_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedStoredLineCapability},
                        {"parallel_arm_march_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedParallelArmCapability},
                        {"parallel_transport_arm_state_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedParallelTransportArmStateCapability},
                        {"transport_arm_state_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedTransportArmStateCapability},
                        {"transport_arm_contraction_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedTransportArmContractionCapability},
                        {"transport_pair_contraction_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedTransportPairContractionCapability},
                        {"transport_pair_stream_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedTransportPairStreamCapability},
                        {"transport_endpoint_batch_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedTransportEndpointBatchCapability},
                        {"certified_tail_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRegularTailMajorantCapability},
                        {"certified_line_integration_capability",
                         domain == "symbolic"
                             ? "unsupported"
                             : kRetainedCertifiedLineCapability}};
  }

  if (operation == "session.close") {
    const auto handle = required_string(root, "session");
    auto& registry = session_registry();
    std::shared_ptr<SolverSession> removed;
    {
      std::lock_guard<std::mutex> lock(registry.mutex);
      const auto found = registry.sessions.find(handle);
      if (found == registry.sessions.end())
        throw std::invalid_argument("unknown or already closed solver session");
      removed = std::move(found->second);
      registry.sessions.erase(found);
    }
    std::size_t charts = 0, regular_equation_owners = 0,
                locals = 0, matches = 0, endpoints = 0,
                tile_plans = 0, transport_states = 0,
                line_results = 0, transport_pair_streams = 0, sccs = 0;
    std::vector<std::shared_ptr<TransportPairObservableStream>> streams;
    {
      std::lock_guard<std::mutex> lock(removed->mutex);
      removed->closed = true;
      // In-flight solve/match/transport/endpoint calls own their reservations
      // and decrement them on exactly one completion path. Do not reset
      // pending counters.
      charts = removed->charts.size();
      regular_equation_owners =
          removed->regular_equation_owners.size();
      locals = removed->locals.size();
      matches = removed->matches.size();
      endpoints = removed->endpoints.size();
      tile_plans = removed->tile_plans.size();
      transport_states = removed->transport_states.size();
      line_results = removed->line_results.size();
      transport_pair_streams = removed->transport_pair_streams.size();
      sccs = removed->sccs.size();
      streams.reserve(transport_pair_streams);
      for (auto& [ignored, stream] : removed->transport_pair_streams)
        streams.push_back(std::move(stream));
      if (removed->pending_line_integrations < transport_pair_streams)
        throw std::logic_error(
            "session close found inconsistent transport-pair stream reservations");
      removed->pending_line_integrations -= transport_pair_streams;
      removed->transport_pair_streams.clear();
      removed->line_results.clear();
      removed->transport_states.clear();
      removed->tile_plans.clear();
      removed->endpoints.clear();
      removed->matches.clear();
      removed->locals.clear();
      removed->sccs.clear();
      removed->scc_handles_by_key.clear();
      removed->charts.clear();
      removed->handles_by_key.clear();
      removed->regular_equation_owners.clear();
      removed->regular_equation_owner_handles_by_key.clear();
    }
    for (const auto& stream : streams)
      stream->abort();
    return json::object{{"status", "ok"}, {"closed", handle},
                        {"released_charts", charts},
                        {"released_regular_equation_owners",
                         regular_equation_owners},
                        {"released_locals", locals},
                        {"released_matches", matches},
                        {"released_endpoints", endpoints},
                        {"released_tile_plans", tile_plans},
                        {"released_transport_states", transport_states},
                        {"released_line_results", line_results},
                        {"released_transport_pair_streams",
                         transport_pair_streams},
                        {"released_scc_charts", sccs}};
  }

  const auto session = find_session(required_string(root, "session"));

  if (operation == "local.specialize_rational_shadow") {
    require_exact_keys(root,
        {"schema", "op", "session", "source_session", "source_local",
         "source_checkpoint_identity", "target_scc",
         "rational_shadow_identity", "checkpoint_identity"},
        "local.specialize_rational_shadow request");
    if (session->domain != "acb")
      throw std::domain_error(
          "Rational-shadow specialization requires a target Acb session");
    const auto source_session = find_session(
        required_string(root, "source_session"));
    if (source_session.get() == session.get() ||
        source_session->domain != "rational")
      throw std::domain_error(
          "Rational-shadow specialization requires a distinct Rational source session");

    std::shared_ptr<StoredLocal<Rational>> source;
    {
      std::lock_guard<std::mutex> lock(source_session->mutex);
      if (source_session->closed)
        throw std::invalid_argument(
            "Rational-shadow source session is closed");
      const auto found = source_session->locals.find(
          required_string(root, "source_local"));
      if (found == source_session->locals.end())
        throw std::invalid_argument(
            "unknown or released Rational-shadow source local");
      source = std::dynamic_pointer_cast<StoredLocal<Rational>>(found->second);
      if (!source)
        throw std::logic_error(
            "Rational-shadow source local differs from its Rational session");
    }
    if (source->checkpoint_identity() !=
            required_string(root, "source_checkpoint_identity") ||
        !source->pseudo_hits().empty() ||
        !source->column_provenance().has_value() ||
        source->retained_derivation().has_value() ||
        !source->residual_binding().has_value())
      throw std::invalid_argument(
          "Rational-shadow source is not a completed, floor-certified homogeneous SCC column");
    auto source_owner = std::dynamic_pointer_cast<CompositeSCCChartBase>(
        source->retained_equation_owner());
    if (!source_owner ||
        std::string(source_owner->equation_scalar_domain()) != "rational")
      throw std::invalid_argument(
          "Rational-shadow source lacks its exact composite SCC owner");

    std::shared_ptr<CompositeSCCChartBase> target;
    std::string local_handle;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("target Acb session is closed");
      const auto found = session->sccs.find(
          required_string(root, "target_scc"));
      if (found == session->sccs.end())
        throw std::invalid_argument(
            "unknown or released target Acb SCC");
      if (session->locals.size() + session->pending_local_solves >=
          session->local_capacity)
        throw std::invalid_argument("persistent local capacity is exhausted");
      target = found->second;
      local_handle = "l:" + std::to_string(session->next_local++);
      ++session->pending_local_solves;
    }

    std::shared_ptr<StoredLocalBase> imported;
    try {
      const auto shadow_identity = required_string(
          root, "rational_shadow_identity");
      if (shadow_identity.empty() ||
          source_owner->rational_shadow_identity() != shadow_identity ||
          target->rational_shadow_identity() != shadow_identity ||
          source_owner->geometry_record() != target->geometry_record() ||
          source_owner->dimension() != target->dimension())
        throw std::invalid_argument(
            "Rational shadow and target Acb SCC lack an identical exact domain-independent owner binding");
      const auto& source_work = source_owner->work_contract();
      const auto& target_work = target->work_contract();
      if (source_work.work_min != target_work.work_min ||
          source_work.requested_min != target_work.requested_min ||
          source_work.requested_max != target_work.requested_max ||
          source_work.work_complete_max != target_work.work_complete_max ||
          source_work.public_t_order != target_work.public_t_order ||
          source_work.work_t_order != target_work.work_t_order ||
          source->solution().epsilon.complete_max <
              target_work.requested_max ||
          source->solution().epsilon.complete_max >
              target_work.work_complete_max ||
          (source->top_valid() != kCompleteInfinity &&
           source->top_valid() <
               source->solution().epsilon.complete_max) ||
          source->solution().epsilon.min_power <
              target_work.work_min ||
          source->solution().taylor_complete_max !=
              target_work.public_t_order)
        throw std::invalid_argument(
            "Rational shadow and target Acb SCC work/public windows differ");

      const auto source_physical =
          std::static_pointer_cast<const PreparedPhysicalClearedODE<Rational>>(
              source_owner->physical_ode_erased());
      const auto target_physical =
          std::static_pointer_cast<const PreparedPhysicalClearedODE<ComplexBall>>(
              target->physical_ode_erased());
      if (!source_physical || !target_physical ||
          !acb_physical_ode_contains_rational_shadow(
              *target_physical, *source_physical))
        throw std::invalid_argument(
            "target Acb physical equation does not enclose the exact Rational shadow equation");

      AcbPrecisionLease lease(session->precision_bits);
      ComplexBall::set_precision(session->precision_bits);
      const auto target_geometry = parse_retained_composite_geometry(
          json::parse(target->geometry_record()));
      const auto checkpoint_identity = required_string(
          root, "checkpoint_identity");
      auto solution = specialize_rational_local_to_acb(
          source->solution(), target_geometry, checkpoint_identity);
      const auto& source_column = *source->column_provenance();
      json::object column_record{
          {"schema", "diffexp2-native-scc-acb-rational-shadow-column-v1"},
          {"scc_exact_identity", target->exact_identity()},
          {"basis_index", source_column.basis_index},
          {"seed_block", source_column.seed_block},
          {"rational_shadow_identity", shadow_identity},
          {"rational_source_column_identity",
           source_column.exact_column_identity},
          {"pseudo_compensation",
           "exact-rational-shadow-case-p-floor-certified-v1"}};
      SCCColumnProvenance target_column{
          target->handle(), target->exact_identity(),
          source_column.seed_block, source_column.basis_index,
          json::serialize(canonical_json_value(column_record))};
      NativeLocalDiagnostics diagnostics;
      diagnostics.top_valid = source->top_valid();
      auto shadow_witness = std::make_shared<RationalShadowColumnWitness>(
          RationalShadowColumnWitness{
              std::make_shared<LocalSolution<Rational>>(source->solution()),
              shadow_identity, source_column.exact_column_identity,
              target_column.exact_column_identity});
      imported = make_retained_typed_shared<ComplexBall,
          StoredLocal<ComplexBall>>(
          local_handle, target->handle(), target->exact_identity(),
          std::move(solution), session->precision_bits,
          std::vector<PseudoHit<ComplexBall>>{}, diagnostics,
          std::move(target_column), std::nullopt, nullptr,
          unavailable_tail_model(
              "Rational-shadow specialization does not import a numeric tail model"),
          std::nullopt, true, true, target, target_physical, "",
          std::move(shadow_witness));
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "Rational-shadow specialization reservation underflow");
      --session->pending_local_solves;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "Rational-shadow specialization reservation underflow");
      --session->pending_local_solves;
      if (session->closed)
        throw std::invalid_argument(
            "target Acb session closed during Rational-shadow specialization");
      if (!session->locals.emplace(local_handle, imported).second)
        throw std::logic_error(
            "Rational-shadow specialization produced a duplicate local handle");
      ++session->total_local_solves;
    }
    auto response = imported->summary();
    response["status"] = "ok";
    response["session"] = session->handle;
    response["scc"] = target->handle();
    response["source_session"] = source_session->handle;
    response["source_local"] = source->handle();
    response["rational_shadow_identity"] =
        target->rational_shadow_identity();
    response["execution_capability"] =
        target->column_execution_capability();
    response["specialization_capability"] =
        "exact-rational-shadow-to-acb-local-v1";
    response["native_retained"] = true;
    response["json_coefficients"] = 0;
    return response;
  }

  if (operation == "checkpoint.save") {
    require_exact_keys(root,
        {"schema", "op", "session", "path", "checkpoint_identity"},
        "checkpoint.save request");
    const auto path = required_string(root, "path");
    const auto checkpoint_identity = required_string(
        root, "checkpoint_identity");
    std::lock_guard<std::mutex> lock(session->mutex);
    auto snapshot = make_checkpoint_snapshot(*session, checkpoint_identity);
    checkpoint::write_atomic(path, snapshot.header, snapshot.payload);
    session->checkpoint_generation = snapshot.generation;
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"path", path}, {"checkpoint_identity", checkpoint_identity},
        {"generation", snapshot.generation},
        {"charts", snapshot.charts},
        {"regular_equation_owners",
         snapshot.regular_equation_owners},
        {"sccs", snapshot.sccs},
        {"locals", snapshot.locals},
        {"exact_matches", snapshot.exact_matches},
        {"acb_matches", snapshot.acb_matches},
        {"planned_match_hops", snapshot.planned_matches},
        {"endpoints", snapshot.endpoints},
        {"tile_plans", snapshot.tile_plans},
        {"transport_states", snapshot.transport_states},
        {"line_results", snapshot.line_results},
        {"serialized_handle_kinds",
         json::array{"chart", "regular-equation-owner", "scc", "local",
                     "exact-rational-match",
                     "acb-match", "planned-match-hop",
                     "materialized-local", "endpoint", "tile",
                     "transport-arm-state", "line"}},
        {"deferred_handle_kinds",
         json::array{"symbolic-local"}},
        {"atomic", true}};
  }

  if (operation == "tile.plan") {
    require_exact_keys(root,
        {"schema", "op", "session", "checkpoint_identity",
         "division_order", "lower", "upper"},
        "native tile.plan request");
    const auto& lower_request = as_object(
        root.at("lower"), "lower native tile arm");
    const auto& upper_request = as_object(
        root.at("upper"), "upper native tile arm");
    const auto lower_handles = parse_plan_chart_handles(lower_request);
    const auto upper_handles = parse_plan_chart_handles(upper_request);
    std::vector<RetainedPlanChartBinding::Owner> lower_charts;
    std::vector<RetainedPlanChartBinding::Owner> upper_charts;
    std::string plan_handle;
    {
      // Resolve every prepared-chart, composite-SCC, or frame-independent
      // regular equation owner and acquire strong
      // typed ownership in one admission section.  Public chart.release or
      // scc.release cannot invalidate either independently executable arm
      // after this point.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      if (session->tile_plans.size() + session->pending_tile_plans >=
          session->tile_plan_capacity)
        throw std::invalid_argument(
            "persistent native tile-plan capacity is exhausted");
      const auto resolve_owner = [&](const std::string& handle,
                                     const char* arm_name)
          -> RetainedPlanChartBinding::Owner {
        if (handle.starts_with("c:")) {
          const auto found = session->charts.find(handle);
          if (found != session->charts.end()) return found->second;
        } else if (handle.starts_with("scc:")) {
          const auto found = session->sccs.find(handle);
          if (found != session->sccs.end()) return found->second;
        } else if (handle.starts_with("eq:")) {
          const auto found = session->regular_equation_owners.find(handle);
          if (found != session->regular_equation_owners.end())
            return found->second;
        }
        throw std::invalid_argument(
            std::string(
                "unknown or released prepared-chart/composite-SCC/regular-equation owner in ") +
            arm_name + " native tile arm: " + handle);
      };
      for (const auto& handle : lower_handles)
        lower_charts.push_back(resolve_owner(handle, "lower"));
      for (const auto& handle : upper_handles)
        upper_charts.push_back(resolve_owner(handle, "upper"));
      plan_handle = "tile:" +
          std::to_string(session->next_tile_plan++);
      ++session->pending_tile_plans;
    }
    std::shared_ptr<StoredTilePlan> plan;
    try {
      plan = build_tile_plan(plan_handle, root, lower_charts, upper_charts);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_tile_plans == 0)
        throw std::logic_error(
            "native tile-plan reservation accounting underflow");
      --session->pending_tile_plans;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_tile_plans == 0)
        throw std::logic_error(
            "native tile-plan reservation accounting underflow");
      --session->pending_tile_plans;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during native tile planning");
      session->tile_plans.emplace(plan_handle, plan);
      ++session->total_tile_plans;
      session->total_tile_plan_ms +=
          plan->summary(false).at("elapsed_ms").as_double();
    }
    auto result = plan->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "tile.plan_arm") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "checkpoint_identity",
         "division_order", "arm"},
        "native tile.plan_arm request");
    const auto& arm_request = as_object(
        root.at("arm"), "native single tile arm");
    const auto handles = parse_plan_chart_handles(arm_request);
    std::vector<RetainedPlanChartBinding::Owner> charts;
    std::string plan_handle;
    {
      // Acquire every chart/SCC owner and the publication reservation in one
      // admission section. The plan keeps these owners alive even if their
      // public handles are released while exact planning runs.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      if (session->tile_plans.size() + session->pending_tile_plans >=
          session->tile_plan_capacity)
        throw std::invalid_argument(
            "persistent native tile-plan capacity is exhausted");
      const auto resolve_owner = [&](const std::string& handle)
          -> RetainedPlanChartBinding::Owner {
        if (handle.starts_with("c:")) {
          const auto found = session->charts.find(handle);
          if (found != session->charts.end()) return found->second;
        } else if (handle.starts_with("scc:")) {
          const auto found = session->sccs.find(handle);
          if (found != session->sccs.end()) return found->second;
        } else if (handle.starts_with("eq:")) {
          const auto found = session->regular_equation_owners.find(handle);
          if (found != session->regular_equation_owners.end())
            return found->second;
        }
        throw std::invalid_argument(
            "unknown or released prepared-chart/composite-SCC/regular-equation owner in native single tile arm: " +
            handle);
      };
      charts.reserve(handles.size());
      for (const auto& handle : handles)
        charts.push_back(resolve_owner(handle));
      plan_handle = "tile:" +
          std::to_string(session->next_tile_plan++);
      ++session->pending_tile_plans;
    }
    std::shared_ptr<StoredTilePlan> plan;
    try {
      plan = build_single_arm_tile_plan(
          plan_handle, root, charts);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_tile_plans == 0)
        throw std::logic_error(
            "native single-arm tile-plan reservation accounting underflow");
      --session->pending_tile_plans;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_tile_plans == 0)
        throw std::logic_error(
            "native single-arm tile-plan reservation accounting underflow");
      --session->pending_tile_plans;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during native single-arm tile planning");
      if (!session->tile_plans.emplace(plan_handle, plan).second)
        throw std::logic_error(
            "native single-arm tile-plan handle collided during publication");
      ++session->total_tile_plans;
      session->total_tile_plan_ms +=
          plan->summary(false).at("elapsed_ms").as_double();
    }
    auto result = plan->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "tile.stats" || operation == "tile.match_interval" ||
      operation == "tile.integration_interval") {
    const auto plan_handle = required_string(root, "tile_plan");
    std::shared_ptr<StoredTilePlan> plan;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->tile_plans.find(plan_handle);
      if (found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or released native tile plan");
      plan = found->second;
    }
    json::object result;
    if (operation == "tile.stats") {
      result = plan->summary();
    } else {
      const auto arm = required_string(root, "arm");
      const auto index = static_cast<std::size_t>(as_u64(
          root.at(operation == "tile.match_interval" ? "match" : "tile"),
          operation == "tile.match_interval" ? "tile-plan match index"
                                               : "tile-plan tile index"));
      result = operation == "tile.match_interval"
          ? plan->match_interval(arm, index)
          : plan->tile_interval(arm, index);
    }
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "tile.release") {
    const auto plan_handle = required_string(root, "tile_plan");
    std::shared_ptr<StoredTilePlan> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->tile_plans.find(plan_handle);
      if (found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or already released native tile plan");
      removed = std::move(found->second);
      session->tile_plans.erase(found);
    }
    return json::object{
        {"status", "ok"}, {"released", plan_handle},
        {"checkpoint_identity", removed->checkpoint_identity()}};
  }

  if (operation == "tile.match_advance") {
    if (session->domain == "rational")
      require_exact_keys(
          root,
          {"schema", "op", "session", "tile_plan", "arm", "match",
           "basis", "incoming", "epsilon", "checkpoint_identity"},
          "native tile.match_advance request");
    else if (session->domain == "acb")
      require_exact_keys(
          root,
          {"schema", "op", "session", "tile_plan", "arm", "match",
           "basis", "incoming", "epsilon", "refinement",
           "exact_lattice", "checkpoint_identity"},
          "native Acb tile.match_advance request");
    else
      throw std::invalid_argument(
          "native tile.match_advance requires rational or Acb coefficients");

    const auto plan_handle = required_string(root, "tile_plan");
    const auto incoming_handle = required_string(root, "incoming");
    const auto& raw_basis = as_array(
        root.at("basis"), "planned local match basis");
    if (raw_basis.empty())
      throw std::invalid_argument(
          "planned local match basis cannot be empty");
    std::vector<std::string> basis_handles;
    basis_handles.reserve(raw_basis.size());
    std::set<std::string> unique_handles;
    for (const auto& raw_handle : raw_basis) {
      if (!raw_handle.is_string() || raw_handle.as_string().empty())
        throw std::invalid_argument(
            "planned local match basis handles must be nonempty strings");
      std::string handle(raw_handle.as_string());
      if (!unique_handles.insert(handle).second)
        throw std::invalid_argument(
            "planned local match basis handles must be pairwise distinct");
      basis_handles.push_back(std::move(handle));
    }
    if (unique_handles.contains(incoming_handle))
      throw std::invalid_argument(
          "planned incoming local must be distinct from its basis");

    std::shared_ptr<StoredTilePlan> plan;
    std::vector<std::shared_ptr<StoredLocalBase>> basis;
    std::shared_ptr<StoredLocalBase> incoming;
    std::string match_handle;
    {
      // Admission is the only serialized section.  Each lower/upper hop owns
      // an immutable plan snapshot and strong local references while its
      // matching arithmetic runs independently outside the session lock.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto plan_found = session->tile_plans.find(plan_handle);
      if (plan_found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or released tile plan for planned local matching");
      if (session->matches.size() + session->pending_matches >=
          session->match_capacity)
        throw std::invalid_argument(
            "persistent local match capacity is exhausted");
      plan = plan_found->second;
      for (const auto& handle : basis_handles) {
        const auto found = session->locals.find(handle);
        if (found == session->locals.end())
          throw std::invalid_argument(
              "unknown or released native local in planned match basis: " +
              handle);
        basis.push_back(found->second);
      }
      const auto incoming_found = session->locals.find(incoming_handle);
      if (incoming_found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released incoming native local for planned match: " +
            incoming_handle);
      incoming = incoming_found->second;
      match_handle = "m:" + std::to_string(session->next_match++);
      ++session->pending_matches;
    }

    std::shared_ptr<StoredPlannedMatchHop> match;
    try {
      match = build_planned_match_hop(
          match_handle, root, session->domain, session->precision_bits,
          checkpoint_configuration_identity(*session), plan, basis_handles,
          basis, incoming_handle, incoming);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches == 0)
        throw std::logic_error(
            "native planned match reservation accounting underflow");
      --session->pending_matches;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches == 0)
        throw std::logic_error(
            "native planned match reservation accounting underflow");
      --session->pending_matches;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during plan-driven local matching");
      session->matches.emplace(match_handle, match);
      ++session->total_local_matches;
      session->total_local_match_ms += match->elapsed_ms();
      plan->note_match_advance(required_string(root, "arm"));
    }
    auto response = match->summary();
    response["status"] = "ok";
    response["session"] = session->handle;
    return response;
  }
