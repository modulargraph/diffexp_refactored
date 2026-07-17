struct SessionCheckpointSnapshot {
  json::object header;
  json::object payload;
  std::uint64_t generation = 0;
  std::size_t charts = 0;
  std::size_t sccs = 0;
  std::size_t locals = 0;
  std::size_t exact_matches = 0;
  std::size_t acb_matches = 0;
  std::size_t planned_matches = 0;
  std::size_t endpoints = 0;
  std::size_t tile_plans = 0;
  std::size_t transport_states = 0;
  std::size_t line_results = 0;
};

SessionCheckpointSnapshot make_checkpoint_snapshot(
    SolverSession& session, const std::string& checkpoint_identity) {
  if (session.closed)
    throw std::invalid_argument("cannot checkpoint a closed solver session");
  if (!session.regular_equation_owners.empty())
    throw std::invalid_argument(
        "checkpoint does not yet support frame-independent regular physical equation owners; release every eq: handle before saving");
  if (session.pending_local_solves != 0 || session.pending_matches != 0 ||
      session.pending_endpoint_limits != 0 ||
      session.pending_tile_plans != 0 ||
      session.pending_transport_states != 0 ||
      session.pending_line_integrations != 0)
    throw std::invalid_argument(
        "checkpoint requires a quiescent session with no pending local solve, match, endpoint limit, tile plan, transport arm, or line integration");

  // Serialize the strong-ownership closure, not only the public registries.
  // Retained lines, tile plans, matches, and materialized locals deliberately
  // survive release of their source objects.  Serialize that immutable
  // closure while recording public registry visibility separately.
  std::unordered_map<std::string, std::shared_ptr<PreparedChartBase>>
      chart_closure = session.charts;
  std::unordered_map<std::string, std::shared_ptr<CompositeSCCChartBase>>
      scc_closure = session.sccs;
  std::unordered_map<std::string, std::shared_ptr<StoredLocalBase>>
      local_closure;
  std::unordered_map<std::string, std::shared_ptr<StoredTilePlan>>
      tile_closure;
  std::unordered_map<std::string, std::shared_ptr<StoredMatchBase>>
      match_closure;
  std::unordered_map<std::string, std::shared_ptr<StoredTransportArmState>>
      transport_closure;
  std::function<void(const std::shared_ptr<PreparedChartBase>&)> add_chart;
  std::function<void(const std::shared_ptr<StoredTilePlan>&)> add_tile;
  std::function<void(const std::shared_ptr<StoredLocalBase>&)> add_local;
  std::function<void(const std::shared_ptr<StoredMatchBase>&)> add_match;
  std::function<void(const std::shared_ptr<StoredTransportArmState>&)>
      add_transport;
  std::function<void(const std::shared_ptr<CompositeSCCChartBase>&)> add_scc;
  add_chart = [&](const std::shared_ptr<PreparedChartBase>& chart) {
    if (!chart)
      throw std::logic_error("checkpoint ownership closure contains a null chart");
    const auto [found, inserted] = chart_closure.emplace(chart->handle(), chart);
    if (!inserted && found->second.get() != chart.get())
      throw std::logic_error(
          "checkpoint ownership closure contains distinct charts with one handle");
  };
  add_scc = [&](
      const std::shared_ptr<CompositeSCCChartBase>& composite) {
    if (!composite)
      throw std::logic_error(
          "checkpoint ownership closure contains a null SCC");
    const auto [found, inserted] =
        scc_closure.emplace(composite->handle(), composite);
      if (!inserted && found->second.get() != composite.get())
        throw std::logic_error(
            "checkpoint ownership closure contains distinct SCCs with one handle");
      if (inserted)
        for (const auto& chart : composite->dependency_charts())
          add_chart(chart);
  };
  add_tile = [&](const std::shared_ptr<StoredTilePlan>& plan) {
    if (!plan)
      throw std::logic_error("checkpoint ownership closure contains a null tile plan");
    const auto [found, inserted] = tile_closure.emplace(plan->handle(), plan);
    if (!inserted) {
      if (found->second.get() != plan.get())
        throw std::logic_error(
            "checkpoint ownership closure contains distinct tile plans with one handle");
      return;
    }
    if (!plan->dependency_regular_equation_owners().empty())
      throw std::domain_error(
          "checkpoint does not yet support a tile plan retaining frame-independent regular physical equation owners");
    for (const auto& composite : plan->dependency_sccs()) add_scc(composite);
    for (const auto& chart : plan->dependency_charts()) add_chart(chart);
  };
  add_match = [&](const std::shared_ptr<StoredMatchBase>& match) {
    if (!match)
      throw std::logic_error("checkpoint ownership closure contains a null match");
    const auto [found, inserted] = match_closure.emplace(
        match->handle(), match);
    if (!inserted) {
      if (found->second.get() != match.get())
        throw std::logic_error(
            "checkpoint ownership closure contains distinct matches with one handle");
      return;
    }
    if (const auto exact =
            std::dynamic_pointer_cast<StoredExactRegularMatch>(match)) {
      for (const auto& local : exact->basis_owners()) add_local(local);
      add_local(exact->incoming_owner());
      return;
    }
    if (const auto hop =
            std::dynamic_pointer_cast<StoredPlannedMatchHop>(match)) {
      add_tile(hop->plan_owner());
      for (const auto& local : hop->basis_owners()) add_local(local);
      add_local(hop->incoming_owner());
      return;
    }
    if (!std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match))
      throw std::logic_error(
          "checkpoint ownership closure contains an unknown retained match implementation");
  };
  add_local = [&](const std::shared_ptr<StoredLocalBase>& local) {
    if (!local)
      throw std::logic_error("checkpoint ownership closure contains a null local");
    const auto [found, inserted] = local_closure.emplace(local->handle(), local);
    if (!inserted) {
      if (found->second.get() != local.get())
        throw std::logic_error(
            "checkpoint ownership closure contains distinct locals with one handle");
      return;
    }
    if (const auto equation_owner = local->retained_equation_owner()) {
      const auto kind = std::string(equation_owner->equation_owner_kind());
      if (kind == "prepared-chart") {
        auto chart = std::dynamic_pointer_cast<PreparedChartBase>(
            equation_owner);
        if (!chart)
          throw std::logic_error(
              "checkpoint local equation owner kind is not a prepared chart");
        add_chart(std::move(chart));
      } else if (kind == "composite-scc") {
        auto composite = std::dynamic_pointer_cast<CompositeSCCChartBase>(
            equation_owner);
        if (!composite)
          throw std::logic_error(
              "checkpoint local equation owner kind is not a CompositeSCC");
        add_scc(std::move(composite));
      } else {
        throw std::domain_error(
            "native checkpoint does not yet serialize this physical equation owner kind");
      }
    }
    if (local->retained_derivation().has_value()) {
      const auto schema = required_string(
          *local->retained_derivation(), "schema");
      if (schema !=
              "diffexp2-retained-plan-match-local-materialization-v1" &&
          schema !=
              "diffexp2-retained-plan-match-local-materialization-v2" &&
          !is_retained_plan_value_handoff_schema(schema) &&
          schema !=
              "diffexp2-retained-rational-row-local-application-v1")
        throw std::domain_error(
            "native checkpoint does not serialize this retained local derivation kind");
      const auto opaque = local->retained_derivation_owner();
      if (!opaque) {
        if ((schema ==
                 "diffexp2-retained-plan-match-local-materialization-v1" ||
             schema ==
                 "diffexp2-retained-plan-match-local-materialization-v2" ||
             is_retained_plan_value_handoff_schema(schema)) &&
            local->has_sealed_plan_match_lineage())
          return;
        throw std::logic_error(
            "checkpoint derived local lost its strong owner");
      }
      if (schema ==
          "diffexp2-retained-plan-match-local-materialization-v1") {
        auto hop = std::static_pointer_cast<StoredPlannedMatchHop>(opaque);
        if (required_string(*local->retained_derivation(), "source_match") !=
            hop->handle())
          throw std::logic_error(
              "checkpoint materialized local derivation names a different owner handle");
        add_match(std::move(hop));
      } else {
        auto source = std::static_pointer_cast<StoredLocalBase>(opaque);
        const auto& source_record = as_object(
            local->retained_derivation()->at("source"),
            "checkpoint rational-row source");
        if (!source || required_string(source_record, "local") !=
                           source->handle())
          throw std::logic_error(
              "checkpoint rational-row derivation names a different source local");
        add_local(std::move(source));
      }
    } else if (local->retained_derivation_owner() != nullptr) {
      throw std::logic_error(
          "checkpoint primitive local unexpectedly retains a derivation owner");
    }
  };
  add_transport = [&](const std::shared_ptr<StoredTransportArmState>& state) {
    if (!state)
      throw std::logic_error(
          "checkpoint ownership closure contains a null transport-arm state");
    const auto [found, inserted] =
        transport_closure.emplace(state->handle(), state);
    if (!inserted) {
      if (found->second.get() != state.get())
        throw std::logic_error(
            "checkpoint ownership closure contains distinct transport-arm states with one handle");
      return;
    }
    add_tile(state->plan_owner());
    add_local(state->anchor_owner());
    for (const auto& basis : state->basis_owners())
      for (const auto& local : basis) add_local(local);
    for (const auto& match : state->matches()) add_match(match);
    for (const auto& local : state->tile_sources()) add_local(local);
  };
  for (const auto& [ignored, chart] : session.charts) add_chart(chart);
  for (const auto& [ignored, plan] : session.tile_plans) add_tile(plan);
  for (const auto& [ignored, local] : session.locals) add_local(local);
  for (const auto& [ignored, match] : session.matches) add_match(match);
  for (const auto& [ignored, state] : session.transport_states)
    add_transport(state);
  for (const auto& [ignored, endpoint] : session.endpoints) {
    if (endpoint->transport_bound()) {
      add_transport(endpoint->transport_owner());
      continue;
    }
    if (!endpoint->plan_bound()) continue;
    const auto& plan = endpoint->plan_owner();
    if (!plan || !endpoint->local_owner())
      throw std::logic_error(
          "checkpoint plan-bound endpoint lost a strong plan/local owner");
    add_tile(plan);
    add_local(endpoint->local_owner());
  }
  for (const auto& [ignored, line] : session.line_results) {
    for (const auto& plan : line->plan_owners()) add_tile(plan);
    for (const auto& local : line->local_owners()) add_local(local);
    for (const auto& state : line->transport_owners())
      add_transport(state);
  }
  for (const auto& [ignored, plan] : tile_closure) {
    for (const auto& composite : plan->dependency_sccs()) add_scc(composite);
    for (const auto& chart : plan->dependency_charts()) add_chart(chart);
  }
  for (const auto& [ignored, composite] : scc_closure)
    for (const auto& chart : composite->dependency_charts()) add_chart(chart);

  std::vector<std::shared_ptr<PreparedChartBase>> charts;
  charts.reserve(chart_closure.size());
  for (const auto& [ignored, chart] : chart_closure) charts.push_back(chart);
  std::sort(charts.begin(), charts.end(), [](const auto& left,
                                             const auto& right) {
    return scoped_handle_id(left->handle(), "c:", "chart") <
           scoped_handle_id(right->handle(), "c:", "chart");
  });
  std::vector<std::shared_ptr<CompositeSCCChartBase>> sccs;
  sccs.reserve(scc_closure.size());
  for (const auto& [ignored, composite] : scc_closure)
    sccs.push_back(composite);
  std::sort(sccs.begin(), sccs.end(), [](const auto& left,
                                         const auto& right) {
    return scoped_handle_id(left->handle(), "scc:", "SCC") <
           scoped_handle_id(right->handle(), "scc:", "SCC");
  });
  std::vector<std::shared_ptr<StoredLocalBase>> locals;
  locals.reserve(local_closure.size());
  for (const auto& [ignored, local] : local_closure)
    locals.push_back(local);
  std::sort(locals.begin(), locals.end(), [](const auto& left,
                                             const auto& right) {
    return scoped_handle_id(left->handle(), "l:", "local") <
           scoped_handle_id(right->handle(), "l:", "local");
  });
  std::vector<std::shared_ptr<StoredMatchBase>> matches;
  matches.reserve(match_closure.size());
  for (const auto& [ignored, match] : match_closure)
    matches.push_back(match);
  std::sort(matches.begin(), matches.end(), [](const auto& left,
                                               const auto& right) {
    return scoped_handle_id(left->handle(), "m:", "match") <
           scoped_handle_id(right->handle(), "m:", "match");
  });
  std::vector<std::shared_ptr<StoredEndpointResult>> endpoints;
  endpoints.reserve(session.endpoints.size());
  for (const auto& [ignored, endpoint] : session.endpoints)
    endpoints.push_back(endpoint);
  std::sort(endpoints.begin(), endpoints.end(), [](const auto& left,
                                                   const auto& right) {
    return scoped_handle_id(left->handle(), "e:", "endpoint") <
           scoped_handle_id(right->handle(), "e:", "endpoint");
  });
  std::vector<std::shared_ptr<StoredTilePlan>> tile_plans;
  tile_plans.reserve(tile_closure.size());
  for (const auto& [ignored, plan] : tile_closure) tile_plans.push_back(plan);
  std::sort(tile_plans.begin(), tile_plans.end(), [](const auto& left,
                                                     const auto& right) {
    return scoped_handle_id(left->handle(), "tile:", "tile plan") <
           scoped_handle_id(right->handle(), "tile:", "tile plan");
  });
  std::vector<std::shared_ptr<StoredTransportArmState>> transport_states;
  transport_states.reserve(transport_closure.size());
  for (const auto& [ignored, state] : transport_closure)
    transport_states.push_back(state);
  std::sort(transport_states.begin(), transport_states.end(),
            [](const auto& left, const auto& right) {
    return scoped_handle_id(left->handle(), "transport:",
                            "transport-arm state") <
           scoped_handle_id(right->handle(), "transport:",
                            "transport-arm state");
  });
  std::vector<std::shared_ptr<StoredLineResult>> line_results;
  line_results.reserve(session.line_results.size());
  for (const auto& [ignored, line] : session.line_results)
    line_results.push_back(line);
  std::sort(line_results.begin(), line_results.end(), [](const auto& left,
                                                         const auto& right) {
    return scoped_handle_id(left->handle(), "line:", "line result") <
           scoped_handle_id(right->handle(), "line:", "line result");
  });

  json::array chart_items;
  chart_items.reserve(charts.size());
  std::set<std::string> chart_handles;
  for (const auto& chart : charts) {
    chart_handles.insert(chart->handle());
    chart_items.push_back(checkpoint_chart_item(session, chart));
  }
  json::array scc_items;
  scc_items.reserve(sccs.size());
  for (const auto& composite : sccs) {
    auto item = checkpoint_scc_item(session, composite);
    const auto& request = as_object(item.at("request"), "SCC request");
    for (const auto& raw_block : as_array(request.at("blocks"), "SCC blocks")) {
      const auto& block = as_object(raw_block, "SCC block");
      if (!chart_handles.contains(required_string(block, "chart")))
        throw std::invalid_argument(
            "checkpoint cannot serialize an SCC whose retained diagonal chart was publicly released");
    }
    scc_items.push_back(std::move(item));
  }
  json::array local_items;
  local_items.reserve(locals.size());
  for (const auto& local : locals)
    local_items.push_back(local->checkpoint_record());
  json::array exact_match_items;
  json::array acb_match_items;
  json::array planned_match_items;
  for (const auto& match : matches) {
    auto record = match->checkpoint_record();
    if (std::dynamic_pointer_cast<StoredPlannedMatchHop>(match))
      planned_match_items.push_back(std::move(record));
    else if (std::dynamic_pointer_cast<StoredExactRegularMatch>(match))
      exact_match_items.push_back(std::move(record));
    else if (std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match))
      acb_match_items.push_back(std::move(record));
    else
      throw std::logic_error(
          "checkpoint encountered an unknown retained match implementation");
  }
  json::array endpoint_items;
  endpoint_items.reserve(endpoints.size());
  for (const auto& endpoint : endpoints)
    endpoint_items.push_back(endpoint->checkpoint_record());
  json::array tile_items;
  tile_items.reserve(tile_plans.size());
  for (const auto& plan : tile_plans)
    tile_items.push_back(plan->checkpoint_record());
  json::array transport_items;
  transport_items.reserve(transport_states.size());
  for (const auto& state : transport_states)
    transport_items.push_back(state->checkpoint_record());
  json::array line_items;
  line_items.reserve(line_results.size());
  for (const auto& line : line_results)
    line_items.push_back(line->checkpoint_record());

  json::array visible_charts;
  for (const auto& chart : charts)
    if (session.charts.contains(chart->handle()))
      visible_charts.emplace_back(chart->handle());
  json::array visible_sccs;
  for (const auto& composite : sccs)
    if (session.sccs.contains(composite->handle()))
      visible_sccs.emplace_back(composite->handle());
  json::array visible_locals;
  for (const auto& local : locals)
    if (session.locals.contains(local->handle()))
      visible_locals.emplace_back(local->handle());
  json::array visible_matches;
  for (const auto& match : matches)
    if (session.matches.contains(match->handle()))
      visible_matches.emplace_back(match->handle());
  json::array visible_tiles;
  for (const auto& plan : tile_plans)
    if (session.tile_plans.contains(plan->handle()))
      visible_tiles.emplace_back(plan->handle());
  json::array visible_transport_states;
  for (const auto& state : transport_states)
    if (session.transport_states.contains(state->handle()))
      visible_transport_states.emplace_back(state->handle());
  json::object registry_visibility{
      {"charts", std::move(visible_charts)},
      {"sccs", std::move(visible_sccs)},
      {"locals", std::move(visible_locals)},
      {"matches", std::move(visible_matches)},
      {"tile_plans", std::move(visible_tiles)},
      {"transport_states", std::move(visible_transport_states)}};

  if (session.checkpoint_generation ==
      std::numeric_limits<std::uint64_t>::max())
    throw std::overflow_error("checkpoint generation counter overflow");
  const auto generation = session.checkpoint_generation + 1;
  auto configuration = checkpoint_configuration_record(session);
  json::object counters{
      {"next_chart", session.next_chart},
      {"next_local", session.next_local},
      {"next_scc", session.next_scc},
      {"next_match", session.next_match},
      {"next_endpoint", session.next_endpoint},
      {"next_tile_plan", session.next_tile_plan},
      {"next_transport_state", session.next_transport_state},
      {"next_line_result", session.next_line_result},
      {"total_local_solves", session.total_local_solves},
      {"total_scc_column_solves", session.total_scc_column_solves},
      {"total_local_matches", session.total_local_matches},
      {"total_endpoint_limits", session.total_endpoint_limits},
      {"total_endpoint_exports", session.total_endpoint_exports},
      {"total_tile_plans", session.total_tile_plans},
      {"total_transport_arm_marches", session.total_transport_arm_marches},
      {"total_transport_contractions",
       session.total_transport_contractions},
      {"total_transport_observables", session.total_transport_observables},
      {"total_transport_pair_contractions",
       session.total_transport_pair_contractions},
      {"total_transport_pair_observables",
       session.total_transport_pair_observables},
      {"total_transport_endpoint_batches",
       session.total_transport_endpoint_batches},
      {"total_transport_endpoint_rows",
       session.total_transport_endpoint_rows},
      {"total_line_integrations", session.total_line_integrations},
      {"total_line_exports", session.total_line_exports},
      {"total_local_run_parse_ms", session.total_local_run_parse_ms},
      {"total_local_kernel_ms", session.total_local_kernel_ms},
      {"total_local_match_ms", session.total_local_match_ms},
      {"total_endpoint_limit_ms", session.total_endpoint_limit_ms},
      {"total_endpoint_export_ms", session.total_endpoint_export_ms},
      {"total_tile_plan_ms", session.total_tile_plan_ms},
      {"total_transport_arm_ms", session.total_transport_arm_ms},
      {"total_transport_contraction_ms",
       session.total_transport_contraction_ms},
      {"total_transport_pair_contraction_ms",
       session.total_transport_pair_contraction_ms},
      {"total_transport_endpoint_batch_ms",
       session.total_transport_endpoint_batch_ms},
      {"total_line_integration_ms", session.total_line_integration_ms},
      {"total_line_export_ms", session.total_line_export_ms},
      {"checkpoint_generation", generation},
      {"checkpoint_restore_count", session.checkpoint_restore_count}};
  json::object session_record{
      {"source_handle", session.handle},
      {"configuration", configuration},
      {"configuration_identity", checkpoint_configuration_identity(session)},
      {"registry_visibility", registry_visibility},
      {"counters", std::move(counters)}};
  json::object payload{
      {"schema", kCheckpointPayloadSchema},
      {"session", std::move(session_record)},
      {"prepared_charts", chart_items},
      {"prepared_scc", scc_items},
      {"retained_locals", local_items},
      {"retained_exact_matches", exact_match_items},
      {"retained_acb_matches", acb_match_items},
      {"retained_planned_match_hops", planned_match_items},
      {"retained_endpoints", endpoint_items},
      {"retained_tile_plans", tile_items},
      {"retained_transport_states", transport_items},
      {"retained_line_results", line_items}};
  json::array mandatory_sections{"session", "prepared_charts",
                                  "prepared_scc", "retained_locals",
                                  "retained_exact_matches",
                                  "retained_acb_matches",
                                  "retained_planned_match_hops",
                                  "retained_endpoints",
                                  "retained_tile_plans",
                                  "retained_transport_states",
                                  "retained_line_results"};
  json::array deferred_kinds{"symbolic-local"};
  json::object header{
      {"format", kCheckpointFormat},
      {"schema", kCheckpointPayloadSchema},
      {"build", checkpoint::kBuildIdentity},
      {"flint", flint_version},
      {"checkpoint_identity", checkpoint_identity},
      {"configuration_identity", checkpoint_configuration_identity(session)},
      {"analytic_identity", json::parse(session.analytic_identity)},
      {"mandatory_sections", std::move(mandatory_sections)},
      {"optional_sections", json::array{}},
      {"deferred_handle_kinds", std::move(deferred_kinds)},
      {"chart_identities", checkpoint_identity_manifest(chart_items)},
      {"scc_identities", checkpoint_identity_manifest(scc_items)},
      {"local_identities", checkpoint_local_identity_manifest(local_items)},
      {"exact_match_identities",
       checkpoint_exact_match_identity_manifest(exact_match_items)},
      {"acb_match_identities",
       checkpoint_acb_match_identity_manifest(acb_match_items)},
      {"planned_match_identities",
       checkpoint_planned_match_identity_manifest(planned_match_items)},
      {"endpoint_identities",
       checkpoint_endpoint_identity_manifest(endpoint_items)},
      {"tile_plan_identities",
       checkpoint_tile_identity_manifest(tile_items)},
      {"transport_state_identities",
       checkpoint_transport_state_identity_manifest(transport_items)},
      {"line_result_identities",
       checkpoint_line_identity_manifest(line_items)},
      {"generation", generation}};
  return {std::move(header), std::move(payload), generation,
          session.charts.size(), session.sccs.size(), session.locals.size(),
          exact_match_items.size(), acb_match_items.size(),
          planned_match_items.size(), endpoints.size(),
          session.tile_plans.size(), session.transport_states.size(),
          line_results.size()};
}

std::vector<std::string> checkpoint_string_array(const json::value& raw,
                                                 const char* label) {
  std::vector<std::string> result;
  for (const auto& value : as_array(raw, label)) {
    if (!value.is_string())
      throw std::invalid_argument(std::string(label) +
                                  " must contain only strings");
    result.emplace_back(value.as_string());
  }
  return result;
}

void validate_checkpoint_envelope(const json::object& header,
                                  const json::object& payload,
                                  const std::string& expected_identity) {
  require_exact_keys(header,
      {"format", "schema", "build", "flint", "checkpoint_identity",
       "configuration_identity", "analytic_identity", "mandatory_sections",
       "optional_sections", "deferred_handle_kinds", "chart_identities",
       "scc_identities", "local_identities", "exact_match_identities",
       "acb_match_identities", "planned_match_identities",
       "endpoint_identities",
       "tile_plan_identities", "transport_state_identities",
       "line_result_identities", "generation"},
      "checkpoint header");
  require_exact_keys(payload,
      {"schema", "session", "prepared_charts", "prepared_scc",
       "retained_locals", "retained_exact_matches",
       "retained_acb_matches", "retained_planned_match_hops",
       "retained_endpoints",
       "retained_tile_plans", "retained_transport_states",
       "retained_line_results"},
      "checkpoint payload");
  if (required_string(header, "format") != kCheckpointFormat ||
      as_u32(header.at("schema"), "checkpoint header schema") !=
          kCheckpointPayloadSchema ||
      as_u32(payload.at("schema"), "checkpoint payload schema") !=
          kCheckpointPayloadSchema)
    throw std::invalid_argument("unsupported native checkpoint schema");
  if (required_string(header, "build") != checkpoint::kBuildIdentity)
    throw std::invalid_argument(
        "native checkpoint solver build identity is incompatible");
  if (required_string(header, "flint") != flint_version)
    throw std::invalid_argument(
        "native checkpoint FLINT build identity is incompatible");
  if (required_string(header, "checkpoint_identity") != expected_identity)
    throw std::invalid_argument(
        "native checkpoint identity differs from the expected identity");
  auto mandatory = checkpoint_string_array(
      header.at("mandatory_sections"), "mandatory checkpoint sections");
  std::sort(mandatory.begin(), mandatory.end());
  const std::vector<std::string> expected_sections{
      "prepared_charts", "prepared_scc", "retained_acb_matches",
      "retained_endpoints", "retained_exact_matches",
      "retained_line_results", "retained_locals",
      "retained_planned_match_hops", "retained_tile_plans",
      "retained_transport_states",
      "session"};
  if (mandatory != expected_sections)
    throw std::invalid_argument(
        "native checkpoint contains unknown or missing mandatory sections");
  if (!as_array(header.at("optional_sections"),
                "optional checkpoint sections").empty())
    throw std::invalid_argument(
        "native checkpoint declares unsupported optional sections");
  auto deferred = checkpoint_string_array(
      header.at("deferred_handle_kinds"),
      "deferred checkpoint handle kinds");
  std::sort(deferred.begin(), deferred.end());
  const std::vector<std::string> expected_deferred{"symbolic-local"};
  if (deferred != expected_deferred)
    throw std::invalid_argument(
        "native checkpoint deferred-state contract is incompatible");

  const auto& session = as_object(payload.at("session"),
                                  "checkpoint session section");
  require_exact_keys(session,
      {"source_handle", "configuration", "configuration_identity",
       "registry_visibility", "counters"}, "checkpoint session section");
  const auto& visibility = as_object(
      session.at("registry_visibility"), "checkpoint registry visibility");
  require_exact_keys(
      visibility,
      {"charts", "sccs", "locals", "matches", "tile_plans",
       "transport_states"},
      "checkpoint registry visibility");
  (void)checkpoint_string_array(visibility.at("charts"),
                                "visible checkpoint charts");
  (void)checkpoint_string_array(visibility.at("sccs"),
                                "visible checkpoint SCCs");
  (void)checkpoint_string_array(visibility.at("locals"),
                                "visible checkpoint locals");
  (void)checkpoint_string_array(visibility.at("matches"),
                                "visible checkpoint matches");
  (void)checkpoint_string_array(visibility.at("tile_plans"),
                                "visible checkpoint tile plans");
  (void)checkpoint_string_array(visibility.at("transport_states"),
                                "visible checkpoint transport states");
  if (required_string(session, "configuration_identity") !=
      required_string(header, "configuration_identity"))
    throw std::invalid_argument(
        "checkpoint header and payload configuration identities disagree");
  const auto& configuration = as_object(
      session.at("configuration"), "checkpoint configuration");
  if (json::serialize(canonical_json_value(configuration)) !=
      required_string(header, "configuration_identity"))
    throw std::invalid_argument(
        "checkpoint configuration does not reproduce its exact identity");
  if (configuration.at("analytic") != header.at("analytic_identity"))
    throw std::invalid_argument(
        "checkpoint analytic-regularization identity is inconsistent");
  const auto& chart_items = as_array(payload.at("prepared_charts"),
                                     "checkpoint prepared charts");
  const auto& scc_items = as_array(payload.at("prepared_scc"),
                                   "checkpoint prepared SCC charts");
  const auto& local_items = as_array(payload.at("retained_locals"),
                                     "checkpoint retained locals");
  const auto& exact_match_items = as_array(
      payload.at("retained_exact_matches"),
      "checkpoint retained exact-rational matches");
  const auto& acb_match_items = as_array(
      payload.at("retained_acb_matches"),
      "checkpoint retained Acb matches");
  const auto& planned_match_items = as_array(
      payload.at("retained_planned_match_hops"),
      "checkpoint retained planned match hops");
  const auto& endpoint_items = as_array(
      payload.at("retained_endpoints"),
      "checkpoint retained endpoints");
  const auto& tile_items = as_array(
      payload.at("retained_tile_plans"),
      "checkpoint retained tile plans");
  const auto& transport_items = as_array(
      payload.at("retained_transport_states"),
      "checkpoint retained transport-arm states");
  const auto& line_items = as_array(
      payload.at("retained_line_results"),
      "checkpoint retained line results");
  if (checkpoint_identity_manifest(chart_items) !=
          as_array(header.at("chart_identities"),
                   "checkpoint chart identities") ||
      checkpoint_identity_manifest(scc_items) !=
          as_array(header.at("scc_identities"),
                   "checkpoint SCC identities") ||
      checkpoint_local_identity_manifest(local_items) !=
          as_array(header.at("local_identities"),
                   "checkpoint local identities") ||
      checkpoint_exact_match_identity_manifest(exact_match_items) !=
          as_array(header.at("exact_match_identities"),
                   "checkpoint exact-match identities") ||
      checkpoint_acb_match_identity_manifest(acb_match_items) !=
          as_array(header.at("acb_match_identities"),
                   "checkpoint Acb match identities") ||
      checkpoint_planned_match_identity_manifest(planned_match_items) !=
          as_array(header.at("planned_match_identities"),
                   "checkpoint planned-match identities") ||
      checkpoint_endpoint_identity_manifest(endpoint_items) !=
          as_array(header.at("endpoint_identities"),
                   "checkpoint endpoint identities") ||
      checkpoint_tile_identity_manifest(tile_items) !=
          as_array(header.at("tile_plan_identities"),
                   "checkpoint tile-plan identities") ||
      checkpoint_transport_state_identity_manifest(transport_items) !=
          as_array(header.at("transport_state_identities"),
                   "checkpoint transport-state identities") ||
      checkpoint_line_identity_manifest(line_items) !=
          as_array(header.at("line_result_identities"),
                   "checkpoint line-result identities"))
    throw std::invalid_argument(
        "checkpoint retained identity manifest is inconsistent");
  if (as_u64(header.at("generation"), "checkpoint generation") !=
      as_u64(as_object(session.at("counters"), "checkpoint counters")
                 .at("checkpoint_generation"),
             "checkpoint generation"))
    throw std::invalid_argument(
        "checkpoint generation differs between header and payload");
}
