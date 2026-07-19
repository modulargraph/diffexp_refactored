json::object restore_checkpoint(const std::string& path,
                                const std::string& expected_identity) {
  const auto container = checkpoint::read_json(path);
  const auto& header_value = container.header;
  const auto& payload_value = container.payload;
  const auto& header = as_object(header_value, "checkpoint JSON header");
  const auto& payload = as_object(payload_value, "checkpoint JSON payload");
  validate_checkpoint_envelope(header, payload, expected_identity);

  const auto& saved_session = as_object(payload.at("session"),
                                        "checkpoint session section");
  const auto& configuration = as_object(
      saved_session.at("configuration"), "checkpoint configuration");
  require_exact_keys(configuration,
      {"domain", "precision_bits", "output_digits", "symbols", "analytic",
       "chart_capacity", "local_capacity", "scc_capacity",
       "match_capacity", "endpoint_capacity", "transport_state_capacity"},
      "checkpoint configuration");
  const auto& raw_visibility = as_object(
      saved_session.at("registry_visibility"),
      "checkpoint registry visibility");
  const auto visibility_set = [&](const char* key, std::string_view prefix) {
    std::set<std::string> result;
    for (const auto& handle : checkpoint_string_array(
             raw_visibility.at(key), key)) {
      (void)scoped_handle_id(handle, prefix, key);
      if (!result.insert(handle).second)
        throw std::invalid_argument(
            std::string("checkpoint registry visibility duplicates a ") + key +
            " handle");
    }
    return result;
  };
  const auto visible_charts = visibility_set("charts", "c:");
  const auto visible_regular_equation_owners =
      visibility_set("regular_equation_owners", "eq:");
  const auto visible_sccs = visibility_set("sccs", "scc:");
  const auto visible_locals = visibility_set("locals", "l:");
  const auto visible_matches = visibility_set("matches", "m:");
  const auto visible_tiles = visibility_set("tile_plans", "tile:");
  const auto visible_transport_states =
      visibility_set("transport_states", "transport:");
  const auto configured_chart_capacity = static_cast<std::size_t>(as_u64(
      configuration.at("chart_capacity"),
      "checkpoint configured chart capacity"));
  if (visible_charts.size() > configured_chart_capacity)
    throw std::invalid_argument(
        "checkpoint visible charts exceed the restored session capacity");
  if (visible_regular_equation_owners.size() >
      configured_chart_capacity)
    throw std::invalid_argument(
        "checkpoint visible regular equation owners exceed the restored session chart capacity");
  const auto configured_scc_capacity = static_cast<std::size_t>(as_u64(
      configuration.at("scc_capacity"),
      "checkpoint configured SCC capacity"));
  if (visible_sccs.size() > configured_scc_capacity)
    throw std::invalid_argument(
        "checkpoint visible SCCs exceed the restored session capacity");
  if (visible_matches.size() > static_cast<std::size_t>(as_u64(
          configuration.at("match_capacity"),
          "checkpoint configured match capacity")))
    throw std::invalid_argument(
        "checkpoint visible matches exceed the restored session capacity");
  const auto configured_transport_state_capacity =
      static_cast<std::size_t>(as_u64(
          configuration.at("transport_state_capacity"),
          "checkpoint configured transport-state capacity"));
  if (visible_transport_states.size() > configured_transport_state_capacity)
    throw std::invalid_argument(
        "checkpoint visible transport states exceed the restored session capacity");
  json::object create{
      {"schema", 2}, {"op", "session.create"},
      {"domain", configuration.at("domain")},
      {"precision_bits", configuration.at("precision_bits")},
      {"output_digits", configuration.at("output_digits")},
      {"symbols", configuration.at("symbols")},
      {"analytic", configuration.at("analytic")},
      {"chart_capacity", configuration.at("chart_capacity")},
      {"local_capacity", configuration.at("local_capacity")},
      {"scc_capacity", configuration.at("scc_capacity")},
      {"match_capacity", configuration.at("match_capacity")},
      {"endpoint_capacity", configuration.at("endpoint_capacity")},
      {"transport_state_capacity",
       configuration.at("transport_state_capacity")}};
  const auto created = run_session_command(create);
  const auto restored_handle = required_string(created, "session");
  bool live = true;
  try {
    const auto restored = find_session(restored_handle);
    {
      // Dependency-only chart/SCC owners are replayed briefly so tile plans
      // can acquire their strong pointers, then removed from the public maps.
      // Their closure may legitimately exceed the public capacities.
      const auto chart_closure_size = as_array(
          payload.at("prepared_charts"),
          "checkpoint prepared chart closure").size();
      const auto regular_equation_owner_closure_size = as_array(
          payload.at("regular_equation_owners"),
          "checkpoint regular equation-owner closure").size();
      const auto scc_closure_size = as_array(
          payload.at("prepared_scc"),
          "checkpoint prepared SCC closure").size();
      std::lock_guard<std::mutex> lock(restored->mutex);
      restored->chart_capacity =
          std::max({configured_chart_capacity, chart_closure_size,
                    regular_equation_owner_closure_size});
      restored->scc_capacity =
          std::max(configured_scc_capacity, scc_closure_size);
    }
    const auto source_handle = required_string(saved_session, "source_handle");
    const auto source_configuration_identity = required_string(
        saved_session, "configuration_identity");
    json::array restored_charts;
    std::set<std::string> all_chart_handles;
    std::uint64_t largest_chart = 0;
    for (const auto& raw_item : as_array(
             payload.at("prepared_charts"), "checkpoint prepared charts")) {
      const auto& item = as_object(raw_item, "checkpoint chart item");
      require_exact_keys(item,
          {"handle", "key", "identity", "signature", "request"},
          "checkpoint chart item");
      const auto old_handle = required_string(item, "handle");
      if (!all_chart_handles.insert(old_handle).second)
        throw std::invalid_argument(
            "checkpoint contains duplicate prepared chart handles");
      const auto handle_id = scoped_handle_id(old_handle, "c:", "chart");
      if (handle_id <= largest_chart)
        throw std::invalid_argument(
            "checkpoint chart handles are not in strict creation order");
      largest_chart = handle_id;
      auto request = as_object(item.at("request"),
                               "checkpoint chart request");
      if (required_string(request, "session") != source_handle ||
          required_string(request, "key") != required_string(item, "key") ||
          required_string(request, "identity") !=
              required_string(item, "identity"))
        throw std::invalid_argument(
            "checkpoint chart request provenance is inconsistent");
      request["session"] = restored_handle;
      {
        std::lock_guard<std::mutex> lock(restored->mutex);
        restored->next_chart = handle_id;
      }
      const auto result = run_session_command(request);
      if (required_string(result, "chart") != old_handle)
        throw std::logic_error(
            "checkpoint chart handle could not be restored exactly");
      {
        std::lock_guard<std::mutex> lock(restored->mutex);
        const auto found = restored->charts.find(old_handle);
        if (found == restored->charts.end() ||
            found->second->signature() != required_string(item, "signature"))
          throw std::invalid_argument(
              "restored chart does not reproduce its exact operator identity");
      }
      if (visible_charts.contains(old_handle))
        restored_charts.push_back(json::object{
            {"chart", old_handle}, {"key", item.at("key")},
            {"identity", item.at("identity")}});
    }
    if (!std::includes(all_chart_handles.begin(), all_chart_handles.end(),
                       visible_charts.begin(), visible_charts.end()))
      throw std::invalid_argument(
          "checkpoint chart visibility names an absent ownership object");

    json::array restored_sccs;
    std::set<std::string> all_scc_handles;
    std::uint64_t largest_scc = 0;
    for (const auto& raw_item : as_array(
             payload.at("prepared_scc"), "checkpoint prepared SCC charts")) {
      const auto& item = as_object(raw_item, "checkpoint SCC item");
      require_exact_keys(item,
          {"handle", "key", "identity", "signature", "request"},
          "checkpoint SCC item");
      const auto old_handle = required_string(item, "handle");
      if (!all_scc_handles.insert(old_handle).second)
        throw std::invalid_argument(
            "checkpoint contains duplicate retained SCC handles");
      const auto handle_id = scoped_handle_id(old_handle, "scc:", "SCC");
      if (handle_id <= largest_scc)
        throw std::invalid_argument(
            "checkpoint SCC handles are not in strict creation order");
      largest_scc = handle_id;
      auto request = as_object(item.at("request"),
                               "checkpoint SCC request");
      if (required_string(request, "session") != source_handle ||
          required_string(request, "key") != required_string(item, "key") ||
          required_string(request, "identity") !=
              required_string(item, "identity"))
        throw std::invalid_argument(
            "checkpoint SCC request provenance is inconsistent");
      request["session"] = restored_handle;
      {
        std::lock_guard<std::mutex> lock(restored->mutex);
        restored->next_scc = handle_id;
      }
      const auto result = run_session_command(request);
      if (required_string(result, "scc") != old_handle)
        throw std::logic_error(
            "checkpoint SCC handle could not be restored exactly");
      {
        std::lock_guard<std::mutex> lock(restored->mutex);
        const auto found = restored->sccs.find(old_handle);
        if (found == restored->sccs.end() ||
            found->second->signature() != required_string(item, "signature"))
          throw std::invalid_argument(
              "restored SCC does not reproduce its exact graph/operator identity");
      }
      if (visible_sccs.contains(old_handle))
        restored_sccs.push_back(json::object{
            {"scc", old_handle}, {"key", item.at("key")},
            {"identity", item.at("identity")}});
    }
    if (!std::includes(all_scc_handles.begin(), all_scc_handles.end(),
                       visible_sccs.begin(), visible_sccs.end()))
      throw std::invalid_argument(
          "checkpoint SCC visibility names an absent ownership object");

    json::array restored_regular_equation_owners;
    std::set<std::string> all_regular_equation_owner_handles;
    std::uint64_t largest_regular_equation_owner = 0;
    for (const auto& raw_item : as_array(
             payload.at("regular_equation_owners"),
             "checkpoint regular equation owners")) {
      const auto& item = as_object(
          raw_item, "checkpoint regular equation-owner item");
      require_exact_keys(
          item, {"handle", "key", "identity", "signature", "request"},
          "checkpoint regular equation-owner item");
      const auto old_handle = required_string(item, "handle");
      if (!all_regular_equation_owner_handles.insert(old_handle).second)
        throw std::invalid_argument(
            "checkpoint contains duplicate regular equation-owner handles");
      const auto handle_id = scoped_handle_id(
          old_handle, "eq:", "regular equation owner");
      if (handle_id <= largest_regular_equation_owner)
        throw std::invalid_argument(
            "checkpoint regular equation-owner handles are not in strict creation order");
      largest_regular_equation_owner = handle_id;
      auto request = as_object(
          item.at("request"), "checkpoint regular equation-owner request");
      if (required_string(request, "session") != source_handle ||
          required_string(request, "key") != required_string(item, "key") ||
          required_string(request, "identity") !=
              required_string(item, "identity"))
        throw std::invalid_argument(
            "checkpoint regular equation-owner request provenance is inconsistent");
      request["session"] = restored_handle;
      {
        std::lock_guard<std::mutex> lock(restored->mutex);
        restored->next_regular_equation_owner = handle_id;
      }
      const auto result = run_session_command(request);
      if (required_string(result, "equation_owner") != old_handle)
        throw std::logic_error(
            "checkpoint regular equation-owner handle could not be restored exactly");
      {
        std::lock_guard<std::mutex> lock(restored->mutex);
        const auto found =
            restored->regular_equation_owners.find(old_handle);
        if (found == restored->regular_equation_owners.end() ||
            found->second->signature() !=
                required_string(item, "signature"))
          throw std::invalid_argument(
              "restored regular equation owner does not reproduce its exact physical identity");
      }
      if (visible_regular_equation_owners.contains(old_handle))
        restored_regular_equation_owners.push_back(json::object{
            {"equation_owner", old_handle},
            {"key", item.at("key")},
            {"identity", item.at("identity")}});
    }
    if (!std::includes(
            all_regular_equation_owner_handles.begin(),
            all_regular_equation_owner_handles.end(),
            visible_regular_equation_owners.begin(),
            visible_regular_equation_owners.end()))
      throw std::invalid_argument(
          "checkpoint regular equation-owner visibility names an absent ownership object");

    std::unique_ptr<AcbPrecisionLease> checkpoint_acb_lease;
    if (restored->domain == "acb") {
      checkpoint_acb_lease =
          std::make_unique<AcbPrecisionLease>(restored->precision_bits);
      ComplexBall::set_precision(restored->precision_bits);
    }

    json::array restored_locals;
    json::array restored_exact_matches;
    json::array restored_acb_matches;
    json::array restored_planned_matches;
    json::array restored_transport_states;
    json::array restored_endpoints;
    json::array restored_tile_plans;
    json::array restored_line_results;
    std::uint64_t largest_local = 0, largest_match = 0;
    std::uint64_t largest_endpoint = 0, largest_tile_plan = 0;
    std::uint64_t largest_transport_state = 0;
    std::uint64_t largest_line_result = 0;
    {
      // Rebuild the ownership DAG under one publication lock. Primitive
      // locals and plans are roots; exact/planned matches and materialized
      // locals are admitted by a fixed-point walk over their strong-owner
      // lineage. No partially restored chain can become observable.
      std::lock_guard<std::mutex> restore_state_lock(restored->mutex);
      const auto& saved_locals = as_array(
          payload.at("retained_locals"), "checkpoint retained locals");
      const auto& saved_tiles = as_array(
          payload.at("retained_tile_plans"),
          "checkpoint retained tile plans");
      const auto& saved_exact = as_array(
          payload.at("retained_exact_matches"),
          "checkpoint retained exact-rational matches");
      const auto& saved_acb = as_array(
          payload.at("retained_acb_matches"),
          "checkpoint retained Acb matches");
      const auto& saved_planned = as_array(
          payload.at("retained_planned_match_hops"),
          "checkpoint retained planned match hops");
      const auto& saved_transport_states = as_array(
          payload.at("retained_transport_states"),
          "checkpoint retained transport-arm states");
      if (visible_locals.size() > restored->local_capacity)
        throw std::invalid_argument(
            "checkpoint visible locals exceed the restored session capacity");
      if (visible_tiles.size() > restored->tile_plan_capacity)
        throw std::invalid_argument(
            "checkpoint visible tile plans exceed the restored session capacity");
      if (visible_transport_states.size() >
          restored->transport_state_capacity)
        throw std::invalid_argument(
            "checkpoint visible transport states exceed the restored session capacity");

      std::set<std::string> all_tile_handles;
      for (const auto& raw_item : saved_tiles) {
        const auto& item = as_object(raw_item,
                                     "checkpoint retained tile plan");
        const auto handle = required_string(item, "handle");
        if (!all_tile_handles.insert(handle).second)
          throw std::invalid_argument(
              "checkpoint contains duplicate retained tile-plan handles");
        const auto id = scoped_handle_id(handle, "tile:", "tile plan");
        if (id <= largest_tile_plan)
          throw std::invalid_argument(
              "checkpoint tile-plan handles are not in strict creation order");
        largest_tile_plan = id;
        auto plan = restore_checkpoint_tile_plan_record(
            item, restored->charts,
            restored->regular_equation_owners, restored->sccs);
        if (plan->checkpoint_record() != raw_item ||
            !restored->tile_plans.emplace(handle, plan).second)
          throw std::invalid_argument(
              "restored tile plan does not reproduce its exact retained state");
        if (visible_tiles.contains(handle))
          restored_tile_plans.push_back(json::object{
              {"tile_plan", handle},
              {"checkpoint_identity", item.at("checkpoint_identity")},
              {"provenance_identity", item.at("provenance_identity")},
              {"generation", header.at("generation")}});
      }
      if (!std::includes(all_tile_handles.begin(), all_tile_handles.end(),
                         visible_tiles.begin(), visible_tiles.end()))
        throw std::invalid_argument(
            "checkpoint tile visibility names an absent ownership object");

      std::set<std::string> all_local_handles;
      std::vector<const json::value*> pending_locals;
      auto validate_local_source = [&](const std::shared_ptr<StoredLocalBase>& local) {
        const auto& source = local->source_chart();
        const bool rational_row_derived =
            local->retained_derivation().has_value() &&
            required_string(*local->retained_derivation(), "schema") ==
                "diffexp2-retained-rational-row-local-application-v1";
        if (source.starts_with("c:")) {
          (void)scoped_handle_id(source, "c:", "local source chart");
          const auto found = restored->charts.find(source);
          if (!rational_row_derived && found != restored->charts.end() &&
              found->second->exact_identity() !=
                  local->source_operator_identity())
            throw std::invalid_argument(
                "checkpoint local source identity disagrees with its restored chart");
        } else if (source.starts_with("scc:")) {
          (void)scoped_handle_id(source, "scc:", "local source SCC");
          const auto found = restored->sccs.find(source);
          if (!rational_row_derived && found != restored->sccs.end() &&
              found->second->exact_identity() !=
                  local->source_operator_identity())
            throw std::invalid_argument(
                "checkpoint local source identity disagrees with its restored SCC");
        } else if (source.starts_with("eq:")) {
          (void)scoped_handle_id(
              source, "eq:", "local source regular equation owner");
          const auto found =
              restored->regular_equation_owners.find(source);
          if (!rational_row_derived &&
              found != restored->regular_equation_owners.end() &&
              found->second->exact_identity() !=
                  local->source_operator_identity())
            throw std::invalid_argument(
                "checkpoint local source identity disagrees with its restored regular equation owner");
        } else {
          throw std::invalid_argument(
              "checkpoint local source is neither a chart, regular equation owner, nor an SCC handle");
        }
        if (local->column_provenance().has_value()) {
          const auto& column = *local->column_provenance();
          const auto found = restored->sccs.find(column.scc_handle);
          if (column.scc_handle != local->source_chart() ||
              column.scc_exact_identity != local->source_operator_identity() ||
              (found != restored->sccs.end() &&
               found->second->exact_identity() != column.scc_exact_identity))
            throw std::invalid_argument(
                "checkpoint local source and SCC-column provenance disagree");
        }
      };
      auto install_local = [&](const json::value& raw_item,
                               std::shared_ptr<void> owner) {
        const auto& item = as_object(raw_item, "checkpoint retained local");
        const auto handle = required_string(item, "handle");
        std::shared_ptr<PhysicalEquationOwnerBase> equation_owner;
        if (!item.at("equation_owner_restore").is_null()) {
          const auto& owner_record = as_object(
              item.at("equation_owner_restore"),
              "checkpoint local physical equation owner");
          const auto owner_kind = required_string(owner_record, "owner_kind");
          const auto owner_handle = required_string(owner_record, "owner_handle");
          if (owner_kind == "prepared-chart") {
            const auto found = restored->charts.find(owner_handle);
            if (found == restored->charts.end())
              throw std::invalid_argument(
                  "checkpoint local physical equation owner chart is missing");
            equation_owner = found->second;
          } else if (owner_kind == "composite-scc") {
            const auto found = restored->sccs.find(owner_handle);
            if (found == restored->sccs.end())
              throw std::invalid_argument(
                  "checkpoint local physical equation owner SCC is missing");
            equation_owner = found->second;
          } else if (owner_kind ==
                     "regular-physical-equation-v1") {
            const auto found =
                restored->regular_equation_owners.find(owner_handle);
            if (found ==
                restored->regular_equation_owners.end())
              throw std::invalid_argument(
                  "checkpoint local regular physical equation owner is missing");
            equation_owner = found->second;
          } else {
            throw std::invalid_argument(
                "checkpoint local has an unsupported physical equation owner kind");
          }
        }
        std::optional<CheckpointValueHandoffPlanBinding>
            value_handoff_plan;
        if (!item.at("retained_derivation").is_null()) {
          const auto& derivation = as_object(
              item.at("retained_derivation"),
              "checkpoint retained-local derivation");
          if (is_retained_plan_value_handoff_schema(
                  required_string(derivation, "schema"))) {
            const auto plan_handle = required_string(
                derivation, "tile_plan");
            const auto plan_found = restored->tile_plans.find(plan_handle);
            if (plan_found == restored->tile_plans.end() ||
                required_string(derivation,
                                "tile_plan_checkpoint_identity") !=
                    plan_found->second->checkpoint_identity() ||
                required_string(derivation,
                                "tile_plan_provenance_identity") !=
                    plan_found->second->provenance_identity())
              throw std::invalid_argument(
                  "checkpoint value handoff differs from its restored tile plan");
            const auto arm_name = required_string(derivation, "arm");
            const auto match_index = checkpoint_size_t(
                derivation.at("match"),
                "checkpoint value handoff match index");
            const auto& arm = plan_found->second->arm(arm_name);
            if (match_index >= arm.exact.matches.size())
              throw std::invalid_argument(
                  "checkpoint value handoff match is outside its restored arm");
            const auto& match = arm.exact.matches[match_index];
            const auto& producing = arm.charts.at(match.producing_chart);
            const auto& receiving = arm.charts.at(match.receiving_chart);
            std::shared_ptr<StoredLocalBase> restored_incoming;
            if (required_string(derivation, "schema") ==
                "diffexp2-retained-plan-value-handoff-v2") {
              const auto& incoming_record = as_object(
                  derivation.at("incoming"),
                  "checkpoint value handoff incoming");
              const auto incoming_handle = required_string(
                  incoming_record, "local");
              const auto incoming_found = restored->locals.find(
                  incoming_handle);
              if (incoming_found == restored->locals.end() ||
                  required_string(incoming_record,
                                  "checkpoint_identity") !=
                      incoming_found->second->checkpoint_identity())
                throw std::invalid_argument(
                    "checkpoint certified value handoff incoming local was not restored before its output");
              restored_incoming = incoming_found->second;
            }
            const auto owner_reference = [](
                const RetainedPlanChartBinding& binding) {
              return std::visit(
                  [&](const auto& owner) {
                    if (!owner)
                      throw std::invalid_argument(
                          "checkpoint value handoff plan chart has a null physical equation owner");
                    return json::object{
                        {"kind", owner->equation_owner_kind()},
                        {"handle", owner->equation_owner_handle()},
                        {"operator_identity",
                         owner->equation_operator_identity()},
                        {"plan_exact_identity",
                         binding.exact_identity},
                        {"owner_signature_identity",
                         owner->owner_signature_identity()},
                        {"physical_payload_identity",
                         owner->physical_payload_identity()}};
                  },
                  binding.owner);
            };
            value_handoff_plan = CheckpointValueHandoffPlanBinding{
                encode_plan_chart(producing, match.producing_chart),
                owner_reference(producing),
                encode_plan_chart(receiving, match.receiving_chart),
                owner_reference(receiving), std::move(restored_incoming),
                exact_plan_rim(producing.prescriptions,
                               producing.geometry.scale)};
          }
        }
        std::shared_ptr<StoredLocalBase> local;
        if (restored->domain == "rational")
          local = restore_checkpoint_local_record<Rational>(
              item, restored->domain, restored->precision_bits, owner,
              equation_owner, value_handoff_plan.has_value()
                  ? &*value_handoff_plan : nullptr);
        else if (restored->domain == "acb")
          local = restore_checkpoint_local_record<ComplexBall>(
              item, restored->domain, restored->precision_bits, owner,
              equation_owner, value_handoff_plan.has_value()
                  ? &*value_handoff_plan : nullptr);
        else
          throw std::invalid_argument(
              "native checkpoint cannot restore symbolic local state");
        if (local->checkpoint_record() != raw_item)
          throw std::invalid_argument(
              "restored local does not reproduce its exact retained state");
        validate_local_source(local);
        if (!restored->locals.emplace(handle, local).second)
          throw std::invalid_argument(
              "checkpoint contains duplicate retained local handles");
        if (visible_locals.contains(handle))
          restored_locals.push_back(json::object{
              {"local", handle}, {"chart", local->source_chart()},
              {"source_operator_identity", local->source_operator_identity()},
              {"checkpoint_identity", local->checkpoint_identity()},
              {"generation", header.at("generation")}});
        return local;
      };
      for (const auto& raw_item : saved_locals) {
        const auto& item = as_object(raw_item, "checkpoint retained local");
        const auto handle = required_string(item, "handle");
        if (!all_local_handles.insert(handle).second)
          throw std::invalid_argument(
              "checkpoint contains duplicate retained local handles");
        const auto id = scoped_handle_id(handle, "l:", "local");
        if (id <= largest_local)
          throw std::invalid_argument(
              "checkpoint local handles are not in strict creation order");
        largest_local = id;
        if (item.at("retained_derivation").is_null() ||
            required_string(item, "schema") ==
                "diffexp2-retained-local-v5")
          (void)install_local(raw_item, nullptr);
        else
          pending_locals.push_back(&raw_item);
      }
      if (!std::includes(all_local_handles.begin(), all_local_handles.end(),
                         visible_locals.begin(), visible_locals.end()))
        throw std::invalid_argument(
            "checkpoint local visibility names an absent ownership object");

      std::set<std::string> all_match_handles;
      const auto register_match_handles = [&](const json::array& records,
                                               const char* label) {
        std::uint64_t previous = 0;
        for (const auto& raw_item : records) {
          const auto& item = as_object(raw_item, label);
          const auto handle = required_string(item, "handle");
          if (!all_match_handles.insert(handle).second)
            throw std::invalid_argument(
                "checkpoint contains duplicate retained match handles");
          const auto id = scoped_handle_id(handle, "m:", "match");
          if (id <= previous)
            throw std::invalid_argument(
                "checkpoint match handles are not in strict creation order");
          previous = id;
          largest_match = std::max(largest_match, id);
        }
      };
      register_match_handles(saved_exact, "checkpoint exact match");
      register_match_handles(saved_acb, "checkpoint Acb match");
      register_match_handles(saved_planned, "checkpoint planned match");
      if (!std::includes(all_match_handles.begin(), all_match_handles.end(),
                         visible_matches.begin(), visible_matches.end()))
        throw std::invalid_argument(
            "checkpoint match visibility names an absent ownership object");

      auto publish_match = [&](const json::value& raw_item,
                               const std::shared_ptr<StoredMatchBase>& match,
                               json::array& response) {
        const auto& item = as_object(raw_item, "checkpoint retained match");
        const auto handle = required_string(item, "handle");
        if (match->checkpoint_record() != raw_item ||
            !restored->matches.emplace(handle, match).second)
          throw std::invalid_argument(
              "restored match does not reproduce its exact retained state");
        if (visible_matches.contains(handle))
          response.push_back(json::object{
              {"match", handle},
              {"checkpoint_identity", item.at("checkpoint_identity")},
              {"provenance_identity", item.at("provenance_identity")},
              {"generation", header.at("generation")}});
      };

      if (!saved_acb.empty() && restored->domain != "acb")
        throw std::invalid_argument(
            "retained Acb match state requires an Acb checkpoint session");
      for (const auto& raw_item : saved_acb) {
        auto match = restore_checkpoint_acb_match_record(
            raw_item, source_configuration_identity);
        publish_match(raw_item, match, restored_acb_matches);
      }
      std::vector<const json::value*> pending_exact;
      std::vector<const json::value*> pending_planned;
      for (const auto& raw_item : saved_exact) pending_exact.push_back(&raw_item);
      for (const auto& raw_item : saved_planned) pending_planned.push_back(&raw_item);

      auto resolve_local = [&](const std::string& handle) {
        const auto found = restored->locals.find(handle);
        return found == restored->locals.end()
            ? std::shared_ptr<StoredLocalBase>() : found->second;
      };
      std::size_t remaining = pending_exact.size() + pending_planned.size() +
                              pending_locals.size();
      while (remaining != 0) {
        bool progress = false;
        for (auto*& raw_ptr : pending_exact) {
          if (raw_ptr == nullptr) continue;
          const auto& item = as_object(*raw_ptr, "checkpoint exact match");
          std::vector<std::shared_ptr<StoredLocalBase>> basis;
          bool ready = true;
          for (const auto& raw_source : as_array(
                   item.at("basis_sources"), "checkpoint exact basis")) {
            auto local = resolve_local(required_string(
                as_object(raw_source, "checkpoint exact basis source"),
                "local"));
            if (!local) { ready = false; break; }
            basis.push_back(std::move(local));
          }
          auto incoming = resolve_local(required_string(
              as_object(item.at("incoming_source"),
                        "checkpoint exact incoming source"), "local"));
          if (!ready || !incoming) continue;
          auto match = restore_checkpoint_exact_match_record(
              item, std::move(basis), incoming);
          publish_match(*raw_ptr, match, restored_exact_matches);
          raw_ptr = nullptr; --remaining; progress = true;
        }
        for (auto*& raw_ptr : pending_planned) {
          if (raw_ptr == nullptr) continue;
          const auto& item = as_object(*raw_ptr, "checkpoint planned match");
          const auto& handoff = as_object(item.at("handoff"),
                                          "checkpoint planned handoff");
          const auto plan_found = restored->tile_plans.find(
              required_string(handoff, "tile_plan"));
          if (plan_found == restored->tile_plans.end()) continue;
          const auto& producing = as_object(handoff.at("producing"),
                                             "checkpoint producing handoff");
          const auto& incoming_source = as_object(
              producing.at("incoming"), "checkpoint planned incoming");
          auto incoming = resolve_local(required_string(incoming_source, "local"));
          if (!incoming) continue;
          std::vector<std::shared_ptr<StoredLocalBase>> basis;
          bool ready = true;
          const auto& receiving = as_object(handoff.at("receiving"),
                                             "checkpoint receiving handoff");
          for (const auto& raw_source : as_array(
                   receiving.at("basis"), "checkpoint planned basis")) {
            auto local = resolve_local(required_string(
                as_object(raw_source, "checkpoint planned basis source"),
                "local"));
            if (!local) { ready = false; break; }
            basis.push_back(std::move(local));
          }
          if (!ready) continue;
          auto match = restore_checkpoint_planned_match_hop_record(
              item, plan_found->second, std::move(basis), incoming,
              source_configuration_identity);
          publish_match(*raw_ptr, match, restored_planned_matches);
          raw_ptr = nullptr; --remaining; progress = true;
        }
        for (auto*& raw_ptr : pending_locals) {
          if (raw_ptr == nullptr) continue;
          const auto& item = as_object(*raw_ptr,
                                       "checkpoint derived local");
          const auto& derivation = as_object(
              item.at("retained_derivation"),
              "checkpoint retained-local derivation");
          const auto derivation_schema = required_string(
              derivation, "schema");
          const auto& lineage = as_object(item.at("retained_owner_lineage"),
                                           "checkpoint local owner lineage");
          if (derivation_schema ==
              "diffexp2-retained-rational-row-local-application-v1") {
            const auto source_handle = required_string(
                lineage, "source_local");
            const auto source = restored->locals.find(source_handle);
            if (source == restored->locals.end()) continue;
            (void)install_local(
                *raw_ptr,
                std::static_pointer_cast<void>(source->second));
            raw_ptr = nullptr; --remaining; progress = true;
            continue;
          }
          if (derivation_schema !=
                  "diffexp2-retained-plan-match-local-materialization-v1" &&
              derivation_schema !=
                  "diffexp2-retained-plan-match-local-materialization-v2")
            throw std::invalid_argument(
                "checkpoint retained local has an unsupported derivation kind");
          const auto found = restored->matches.find(
              required_string(lineage, "match"));
          if (found == restored->matches.end()) continue;
          auto hop = std::dynamic_pointer_cast<StoredPlannedMatchHop>(
              found->second);
          if (!hop)
            throw std::invalid_argument(
                "checkpoint materialized local owner is not a planned match hop");
          auto local = install_local(
              *raw_ptr, std::static_pointer_cast<void>(hop));
          hop->validate_materialized_derivation(
              *local->retained_derivation(), local->scalar_domain());
          hop->validate_materialized_equation_owner(local);
          raw_ptr = nullptr; --remaining; progress = true;
        }
        if (!progress)
          throw std::invalid_argument(
              "checkpoint local/match ownership graph has a missing dependency or cycle");
      }

      std::set<std::string> all_transport_state_handles;
      for (const auto& raw_item : saved_transport_states) {
        const auto& item = as_object(
            raw_item, "checkpoint retained transport-arm state");
        const auto handle = required_string(item, "handle");
        if (!all_transport_state_handles.insert(handle).second)
          throw std::invalid_argument(
              "checkpoint contains duplicate retained transport-state handles");
        const auto id = scoped_handle_id(
            handle, "transport:", "transport-arm state");
        if (id <= largest_transport_state)
          throw std::invalid_argument(
              "checkpoint transport-state handles are not in strict creation order");
        largest_transport_state = id;
        auto state = restore_checkpoint_transport_arm_state_record(
            raw_item, restored->tile_plans, restored->locals,
            restored->matches);
        if (state->checkpoint_record() != raw_item ||
            !restored->transport_states.emplace(handle, state).second)
          throw std::invalid_argument(
              "restored transport-arm state does not reproduce its exact retained state");
        if (visible_transport_states.contains(handle))
          restored_transport_states.push_back(json::object{
              {"transport_state", handle},
              {"checkpoint_identity", item.at("checkpoint_identity")},
              {"provenance_identity", item.at("provenance_identity")},
              {"generation", header.at("generation")}});
      }
      if (!std::includes(all_transport_state_handles.begin(),
                         all_transport_state_handles.end(),
                         visible_transport_states.begin(),
                         visible_transport_states.end()))
        throw std::invalid_argument(
            "checkpoint transport-state visibility names an absent ownership object");

      const auto& saved_endpoints = as_array(
          payload.at("retained_endpoints"),
          "checkpoint retained endpoints");
      if (saved_endpoints.size() > restored->endpoint_capacity)
        throw std::invalid_argument(
            "checkpoint retained endpoints exceed the restored session capacity");
      for (const auto& raw_item : saved_endpoints) {
        const auto& item = as_object(raw_item, "checkpoint retained endpoint");
        const auto handle = required_string(item, "handle");
        const auto id = scoped_handle_id(handle, "e:", "endpoint");
        if (id <= largest_endpoint)
          throw std::invalid_argument(
              "checkpoint endpoint handles are not in strict creation order");
        largest_endpoint = id;
        const auto& source = as_object(item.at("source"),
                                       "checkpoint endpoint source");
        std::shared_ptr<StoredEndpointResult> endpoint;
        const auto endpoint_schema = required_string(item, "schema");
        if (endpoint_schema ==
            "diffexp2-retained-transport-endpoint-result-v1") {
          const auto& state_reference = as_object(
              source.at("transport_state"),
              "checkpoint transport endpoint state reference");
          const auto state_found = restored->transport_states.find(
              required_string(state_reference, "handle"));
          if (state_found == restored->transport_states.end())
            throw std::invalid_argument(
                "checkpoint transport endpoint lost its strongly owned state");
          endpoint = restore_checkpoint_transport_endpoint_record(
              item, restored->domain, state_found->second);
        } else if (endpoint_schema ==
            "diffexp2-retained-plan-bound-endpoint-result-v1") {
          const auto plan_found = restored->tile_plans.find(
              required_string(source, "tile_plan"));
          const auto local_found = restored->locals.find(
              required_string(source, "local"));
          if (plan_found == restored->tile_plans.end() ||
              local_found == restored->locals.end())
            throw std::invalid_argument(
                "checkpoint plan-bound endpoint lost a strongly owned plan or local");
          endpoint = restore_checkpoint_planned_endpoint_record(
              item, restored->domain, plan_found->second,
              local_found->second);
        } else {
          endpoint = restore_checkpoint_endpoint_record(
              item, restored->domain);
          const auto found = restored->locals.find(
              required_string(source, "local"));
          if (found != restored->locals.end() &&
              (found->second->source_chart() !=
                   required_string(source, "chart") ||
               found->second->source_operator_identity() !=
                   required_string(source, "source_operator_identity") ||
               found->second->checkpoint_identity() !=
                   required_string(source, "checkpoint_identity") ||
               found->second->exact_analytic_metadata() !=
                   item.at("analytic_metadata")))
            throw std::invalid_argument(
                "checkpoint endpoint provenance disagrees with its restored source local");
        }
        if (endpoint->checkpoint_record() != raw_item ||
            !restored->endpoints.emplace(handle, endpoint).second)
          throw std::invalid_argument(
              "restored endpoint does not reproduce its exact retained state");
        restored_endpoints.push_back(json::object{
            {"endpoint", handle},
            {"checkpoint_identity", item.at("checkpoint_identity")},
            {"provenance_identity", item.at("provenance_identity")},
            {"generation", header.at("generation")}});
      }

      const auto& saved_lines = as_array(
          payload.at("retained_line_results"),
          "checkpoint retained line results");
      if (saved_lines.size() > restored->line_result_capacity)
        throw std::invalid_argument(
            "checkpoint retained line results exceed the restored session capacity");
      for (const auto& raw_item : saved_lines) {
        const auto& item = as_object(raw_item, "checkpoint retained line result");
        const auto handle = required_string(item, "handle");
        const auto id = scoped_handle_id(handle, "line:", "line result");
        if (id <= largest_line_result)
          throw std::invalid_argument(
              "checkpoint line-result handles are not in strict creation order");
        largest_line_result = id;
        const auto schema = required_string(item, "schema");
        const json::object* source = nullptr;
        const bool compact_pair_line =
            schema ==
                "diffexp2-retained-transport-pair-observable-line-v1" ||
            schema ==
                "diffexp2-retained-transport-pair-observable-line-v2";
        const bool aggregate_line =
            schema == "diffexp2-retained-line-aggregate-v1" ||
            schema == "diffexp2-retained-transport-observable-line-v1" ||
            schema == "diffexp2-retained-transport-observable-line-v2" ||
            compact_pair_line;
        if (aggregate_line)
          source = &as_object(
              as_object(item.at("provenance"),
                        "checkpoint line aggregate provenance").at("source"),
              "checkpoint line aggregate source");
        else
          source = &as_object(item.at("source"), "checkpoint line source");
        std::shared_ptr<StoredLineResult> line;
        if (aggregate_line) {
          std::vector<std::shared_ptr<StoredTilePlan>> plans;
          std::vector<std::shared_ptr<StoredLocalBase>> owners;
          std::vector<std::shared_ptr<StoredTransportArmState>>
              transport_owners;
          if (compact_pair_line) {
            for (const auto* name : {"lower", "upper"}) {
              const auto& side = as_object(
                  source->at(name),
                  "checkpoint transport-pair line source side");
              const auto plan = restored->tile_plans.find(
                  required_string(side, "tile_plan"));
              const auto& state_reference = as_object(
                  side.at("transport_state"),
                  "checkpoint transport-pair state reference");
              const auto state = restored->transport_states.find(
                  required_string(state_reference, "handle"));
              if (plan == restored->tile_plans.end() ||
                  state == restored->transport_states.end())
                throw std::invalid_argument(
                    "checkpoint compact transport-pair line lost an ordered plan or state owner");
              plans.push_back(plan->second);
              transport_owners.push_back(state->second);
            }
          } else {
            const auto plan = restored->tile_plans.find(
                required_string(*source, "tile_plan"));
            if (plan == restored->tile_plans.end())
              throw std::invalid_argument(
                  "checkpoint line result lost its strongly owned plan");
            plans.push_back(plan->second);
          }
          if (schema == "diffexp2-retained-line-aggregate-v1") {
            for (const auto& raw_owner : as_array(
                     source->at("locals"),
                     "checkpoint line aggregate local owners")) {
              const auto& owner = as_object(
                  raw_owner, "checkpoint line aggregate local owner");
              const auto local = restored->locals.find(
                  required_string(owner, "local"));
              if (local == restored->locals.end())
                throw std::invalid_argument(
                    "checkpoint line aggregate lost a strongly owned local");
              owners.push_back(local->second);
            }
          } else if (!compact_pair_line) {
            const auto state = restored->transport_states.find(
                required_string(*source, "transport_state"));
            if (state == restored->transport_states.end())
              throw std::invalid_argument(
                  "checkpoint compact transport line lost its strongly owned state");
            transport_owners.push_back(state->second);
          }
          line = restore_checkpoint_line_aggregate_record(
              item, std::move(plans), std::move(owners),
              std::move(transport_owners));
        } else {
          const auto plan = restored->tile_plans.find(
              required_string(*source, "tile_plan"));
          if (plan == restored->tile_plans.end())
            throw std::invalid_argument(
                "checkpoint line result lost its strongly owned plan");
          const auto local = restored->locals.find(
              required_string(*source, "local"));
          if (local == restored->locals.end())
            throw std::invalid_argument(
                "checkpoint line result lost its strongly owned local");
          line = restore_checkpoint_line_result_record(
              item, plan->second, local->second);
        }
        if (line->checkpoint_record() !=
                normalized_checkpoint_line_record_for_roundtrip(raw_item) ||
            !restored->line_results.emplace(handle, line).second)
          throw std::invalid_argument(
              "restored line result does not reproduce its exact retained state");
        restored_line_results.push_back(json::object{
            {"line", handle},
            {"checkpoint_identity", item.at("checkpoint_identity")},
            {"provenance_identity", item.at("provenance_identity")},
            {"generation", header.at("generation")}});
      }

      for (auto it = restored->transport_states.begin();
           it != restored->transport_states.end();)
        it = visible_transport_states.contains(it->first)
            ? std::next(it) : restored->transport_states.erase(it);
      for (auto it = restored->matches.begin(); it != restored->matches.end();)
        it = visible_matches.contains(it->first)
            ? std::next(it) : restored->matches.erase(it);
      for (auto it = restored->tile_plans.begin(); it != restored->tile_plans.end();)
        it = visible_tiles.contains(it->first)
            ? std::next(it) : restored->tile_plans.erase(it);
      for (auto it = restored->locals.begin(); it != restored->locals.end();)
        it = visible_locals.contains(it->first)
            ? std::next(it) : restored->locals.erase(it);
      for (auto iterator =
               restored->regular_equation_owners.begin();
           iterator != restored->regular_equation_owners.end();) {
        if (visible_regular_equation_owners.contains(iterator->first)) {
          ++iterator;
        } else {
          restored->regular_equation_owner_handles_by_key.erase(
              iterator->second->key());
          iterator =
              restored->regular_equation_owners.erase(iterator);
        }
      }
      for (auto iterator = restored->sccs.begin();
           iterator != restored->sccs.end();) {
        if (visible_sccs.contains(iterator->first)) {
          ++iterator;
        } else {
          restored->scc_handles_by_key.erase(iterator->second->key());
          iterator = restored->sccs.erase(iterator);
        }
      }
      for (auto it = restored->charts.begin(); it != restored->charts.end();) {
        if (visible_charts.contains(it->first)) {
          ++it;
        } else {
          restored->handles_by_key.erase(it->second->key());
          it = restored->charts.erase(it);
        }
      }
    }

    const auto& counters = as_object(saved_session.at("counters"),
                                     "checkpoint counters");
    require_exact_keys(counters,
        {"next_chart", "next_regular_equation_owner", "next_local",
         "next_scc", "next_match",
         "next_endpoint", "next_tile_plan", "next_transport_state",
         "next_line_result",
         "total_local_solves", "total_scc_column_solves",
         "total_local_matches", "total_endpoint_limits",
         "total_endpoint_exports", "total_tile_plans",
         "total_transport_arm_marches",
         "total_transport_contractions", "total_transport_observables",
         "total_transport_pair_contractions",
         "total_transport_pair_observables",
         "total_transport_endpoint_batches",
         "total_transport_endpoint_rows",
         "total_line_integrations", "total_line_exports",
         "total_local_run_parse_ms",
         "total_local_kernel_ms", "total_local_match_ms",
         "total_endpoint_limit_ms", "total_endpoint_export_ms",
         "total_tile_plan_ms", "total_transport_arm_ms",
         "total_transport_contraction_ms",
         "total_transport_pair_contraction_ms",
         "total_transport_endpoint_batch_ms",
         "total_line_integration_ms",
         "total_line_export_ms",
         "checkpoint_generation", "checkpoint_restore_count"},
        "checkpoint counters");
    const auto next_chart = as_u64(counters.at("next_chart"), "next chart");
    const auto next_regular_equation_owner = as_u64(
        counters.at("next_regular_equation_owner"),
        "next regular equation owner");
    const auto next_scc = as_u64(counters.at("next_scc"), "next SCC");
    const auto next_local = as_u64(counters.at("next_local"), "next local");
    const auto next_match = as_u64(counters.at("next_match"), "next match");
    const auto next_endpoint = as_u64(
        counters.at("next_endpoint"), "next endpoint");
    const auto next_tile_plan = as_u64(
        counters.at("next_tile_plan"), "next tile plan");
    const auto next_transport_state = as_u64(
        counters.at("next_transport_state"), "next transport state");
    const auto next_line_result = as_u64(
        counters.at("next_line_result"), "next line result");
    const auto restore_count = as_u64(
        counters.at("checkpoint_restore_count"),
        "checkpoint restore count");
    if (next_chart <= largest_chart ||
        next_regular_equation_owner <=
            largest_regular_equation_owner ||
        next_scc <= largest_scc ||
        next_local <= largest_local || next_match <= largest_match ||
        next_endpoint <= largest_endpoint ||
        next_tile_plan <= largest_tile_plan ||
        next_transport_state <= largest_transport_state ||
        next_line_result <= largest_line_result)
      throw std::invalid_argument(
          "checkpoint next-handle counters do not follow retained handles");
    if (restore_count == std::numeric_limits<std::uint64_t>::max())
      throw std::overflow_error("checkpoint restore counter overflow");
    {
      std::lock_guard<std::mutex> lock(restored->mutex);
      restored->next_chart = next_chart;
      restored->next_regular_equation_owner =
          next_regular_equation_owner;
      restored->next_local = next_local;
      restored->next_scc = next_scc;
      restored->next_match = next_match;
      restored->next_endpoint = next_endpoint;
      restored->next_tile_plan = next_tile_plan;
      restored->next_transport_state = next_transport_state;
      restored->next_line_result = next_line_result;
      restored->chart_capacity = configured_chart_capacity;
      restored->scc_capacity = configured_scc_capacity;
      restored->total_local_solves = as_u64(
          counters.at("total_local_solves"), "total local solves");
      restored->total_scc_column_solves = as_u64(
          counters.at("total_scc_column_solves"),
          "total SCC column solves");
      restored->total_local_matches = as_u64(
          counters.at("total_local_matches"), "total local matches");
      restored->total_endpoint_limits = as_u64(
          counters.at("total_endpoint_limits"), "total endpoint limits");
      restored->total_endpoint_exports = as_u64(
          counters.at("total_endpoint_exports"), "total endpoint exports");
      restored->total_tile_plans = as_u64(
          counters.at("total_tile_plans"), "total tile plans");
      restored->total_transport_arm_marches = as_u64(
          counters.at("total_transport_arm_marches"),
          "total transport-arm marches");
      restored->total_transport_contractions = as_u64(
          counters.at("total_transport_contractions"),
          "total transport contractions");
      restored->total_transport_observables = as_u64(
          counters.at("total_transport_observables"),
          "total transport observables");
      restored->total_transport_pair_contractions = as_u64(
          counters.at("total_transport_pair_contractions"),
          "total transport-pair contractions");
      restored->total_transport_pair_observables = as_u64(
          counters.at("total_transport_pair_observables"),
          "total transport-pair observables");
      restored->total_transport_endpoint_batches = as_u64(
          counters.at("total_transport_endpoint_batches"),
          "total transport endpoint batches");
      restored->total_transport_endpoint_rows = as_u64(
          counters.at("total_transport_endpoint_rows"),
          "total transport endpoint rows");
      restored->total_line_integrations = as_u64(
          counters.at("total_line_integrations"),
          "total line integrations");
      restored->total_line_exports = as_u64(
          counters.at("total_line_exports"), "total line exports");
      restored->total_local_run_parse_ms = as_double(
          counters.at("total_local_run_parse_ms"),
          "total local parse time");
      restored->total_local_kernel_ms = as_double(
          counters.at("total_local_kernel_ms"),
          "total local kernel time");
      restored->total_local_match_ms = as_double(
          counters.at("total_local_match_ms"),
          "total local match time");
      restored->total_endpoint_limit_ms = as_double(
          counters.at("total_endpoint_limit_ms"),
          "total endpoint limit time");
      restored->total_endpoint_export_ms = as_double(
          counters.at("total_endpoint_export_ms"),
          "total endpoint export time");
      restored->total_tile_plan_ms = checkpoint_nonnegative_double(
          counters.at("total_tile_plan_ms"), "total tile-plan time");
      restored->total_transport_arm_ms = checkpoint_nonnegative_double(
          counters.at("total_transport_arm_ms"),
          "total transport-arm time");
      restored->total_transport_contraction_ms =
          checkpoint_nonnegative_double(
              counters.at("total_transport_contraction_ms"),
              "total transport-contraction time");
      restored->total_transport_pair_contraction_ms =
          checkpoint_nonnegative_double(
              counters.at("total_transport_pair_contraction_ms"),
              "total transport-pair contraction time");
      restored->total_transport_endpoint_batch_ms =
          checkpoint_nonnegative_double(
              counters.at("total_transport_endpoint_batch_ms"),
              "total transport endpoint-batch time");
      restored->total_line_integration_ms = checkpoint_nonnegative_double(
          counters.at("total_line_integration_ms"),
          "total line-integration time");
      restored->total_line_export_ms = checkpoint_nonnegative_double(
          counters.at("total_line_export_ms"),
          "total line-export time");
      restored->checkpoint_generation = as_u64(
          counters.at("checkpoint_generation"), "checkpoint generation");
      restored->checkpoint_restore_count = restore_count + 1;
      restored->restored_from_checkpoint_identity = expected_identity;
    }
    live = false;
    return json::object{
        {"status", "ok"}, {"session", restored_handle},
        {"checkpoint_identity", expected_identity},
        {"generation", header.at("generation")},
        {"restore_count", restored->checkpoint_restore_count},
        {"configuration_identity", header.at("configuration_identity")},
        {"analytic_identity", header.at("analytic_identity")},
        {"charts", std::move(restored_charts)},
        {"regular_equation_owners",
         std::move(restored_regular_equation_owners)},
        {"sccs", std::move(restored_sccs)},
        {"locals", std::move(restored_locals)},
        {"exact_matches", std::move(restored_exact_matches)},
        {"acb_matches", std::move(restored_acb_matches)},
        {"planned_match_hops", std::move(restored_planned_matches)},
        {"endpoints", std::move(restored_endpoints)},
        {"tile_plans", std::move(restored_tile_plans)},
        {"transport_states", std::move(restored_transport_states)},
        {"line_results", std::move(restored_line_results)},
        {"deferred_handle_kinds", header.at("deferred_handle_kinds")},
        {"replayed_wolfram_preprocessing", false}};
  } catch (...) {
    if (live) {
      try {
        run_session_command(json::object{{"schema", 2},
                                         {"op", "session.close"},
                                         {"session", restored_handle}});
      } catch (...) {
      }
    }
    throw;
  }
}
