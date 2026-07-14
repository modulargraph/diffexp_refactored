  if (operation == "transport.consume_hop") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "tile_plan",
         "tile_plan_checkpoint_identity", "arm", "match",
         "receiving_basis", "incoming", "incoming_checkpoint_identity",
         "epsilon", "refinement", "checkpoint_policy"},
        "native transport.consume_hop request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "transport.consume_hop requires Rational or Acb coefficients");
    const auto arm_name = required_string(root, "arm");
    const auto match_index = checkpoint_size_t(
        root.at("match"), "consuming transport hop match index");
    const auto& raw_basis = as_array(
        root.at("receiving_basis"),
        "consuming transport hop receiving basis");
    if (raw_basis.empty())
      throw std::invalid_argument(
          "consuming transport hop receiving basis cannot be empty");
    std::vector<std::string> basis_handles;
    std::set<std::string> unique_basis;
    basis_handles.reserve(raw_basis.size());
    for (const auto& raw_handle : raw_basis) {
      if (!raw_handle.is_string() || raw_handle.as_string().empty())
        throw std::invalid_argument(
            "consuming transport hop basis handle must be nonempty");
      std::string handle(raw_handle.as_string());
      if (!unique_basis.insert(handle).second)
        throw std::invalid_argument(
            "consuming transport hop basis handles must be pairwise distinct");
      basis_handles.push_back(std::move(handle));
    }
    const auto incoming_handle = required_string(root, "incoming");
    if (unique_basis.contains(incoming_handle))
      throw std::invalid_argument(
          "consuming transport hop incoming local cannot be in its receiving basis");
    const auto& epsilon = as_object(
        root.at("epsilon"), "consuming transport hop epsilon contract");
    require_exact_keys(epsilon, {"min", "max", "required_complete_max"},
                       "consuming transport hop epsilon contract");
    EpsilonWindow requested_epsilon{
        as_i32(epsilon.at("min"),
               "consuming transport hop epsilon minimum"),
        as_i32(epsilon.at("max"),
               "consuming transport hop epsilon maximum")};
    (void)requested_epsilon.width();
    const auto required_complete_max = as_i32(
        epsilon.at("required_complete_max"),
        "consuming transport hop required epsilon maximum");
    if (required_complete_max < requested_epsilon.min_power ||
        required_complete_max > requested_epsilon.complete_max)
      throw std::invalid_argument(
          "consuming transport hop epsilon contract is inconsistent");
    const auto& refinement = as_object(
        root.at("refinement"),
        "consuming transport hop refinement policy");
    require_exact_keys(refinement, {"relative_tolerance", "max_steps"},
                       "consuming transport hop refinement policy");
    if (required_string(refinement, "relative_tolerance").empty() ||
        as_u32(refinement.at("max_steps"),
               "consuming transport hop refinement steps") > 32)
      throw std::invalid_argument(
          "consuming transport hop refinement policy is invalid");
    const auto& checkpoint_policy = as_object(
        root.at("checkpoint_policy"),
        "consuming transport hop checkpoint policy");
    require_exact_keys(checkpoint_policy, {"schema", "root"},
                       "consuming transport hop checkpoint policy");
    if (required_string(checkpoint_policy, "schema") !=
        "diffexp2-deterministic-arm-checkpoints-v1")
      throw std::invalid_argument(
          "unsupported consuming transport hop checkpoint policy");
    const auto checkpoint_root =
        required_string(checkpoint_policy, "root");

    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> incoming;
    std::vector<std::shared_ptr<StoredLocalBase>> basis;
    json::array basis_reference;
    std::string match_handle;
    std::string local_handle;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto plan_found = session->tile_plans.find(
          required_string(root, "tile_plan"));
      if (plan_found == session->tile_plans.end() ||
          required_string(root, "tile_plan_checkpoint_identity") !=
              plan_found->second->checkpoint_identity())
        throw std::invalid_argument(
            "consuming transport hop tile-plan binding is stale");
      plan = plan_found->second;
      const auto& retained = plan->arm(arm_name);
      if (match_index >= retained.exact.matches.size())
        throw std::invalid_argument(
            "consuming transport hop match lies outside its exact arm");
      const auto incoming_found = session->locals.find(incoming_handle);
      if (incoming_found == session->locals.end() ||
          required_string(root, "incoming_checkpoint_identity") !=
              incoming_found->second->checkpoint_identity())
        throw std::invalid_argument(
            "consuming transport hop incoming-local binding is stale");
      incoming = incoming_found->second;
      basis.reserve(basis_handles.size());
      basis_reference.reserve(basis_handles.size());
      for (const auto& handle : basis_handles) {
        const auto found = session->locals.find(handle);
        if (found == session->locals.end())
          throw std::invalid_argument(
              "unknown or already consumed receiving-basis local: " +
              handle);
        basis.push_back(found->second);
        basis_reference.push_back(
            compact_transport_local_reference(found->second));
      }
      const auto retained_after_commit =
          session->locals.size() - basis.size() +
          session->pending_local_solves + 1;
      if (retained_after_commit > session->local_capacity)
        throw std::invalid_argument(
            "persistent local capacity is exhausted by consuming transport hop");
      match_handle = "m:" + std::to_string(session->next_match++);
      local_handle = "l:" + std::to_string(session->next_local++);
      ++session->pending_local_solves;
    }

    const auto started = std::chrono::steady_clock::now();
    std::shared_ptr<StoredPlannedMatchHop> match;
    std::shared_ptr<StoredLocalBase> next;
    json::object sealed_lineage;
    try {
      const auto live_epsilon = live_match_epsilon_intersection(
          requested_epsilon, required_complete_max, incoming, basis);
      const auto match_checkpoint = arm_checkpoint_identity(
          checkpoint_root, arm_name, "match", match_index + 1);
      json::object match_request{
          {"arm", arm_name}, {"match", match_index},
          {"epsilon", json::object{
               {"min", live_epsilon.min_power},
               {"max", live_epsilon.complete_max},
               {"required_complete_max", required_complete_max}}},
          {"checkpoint_identity", match_checkpoint}};
      if (session->domain == "acb") {
        AcbPrecisionLease lease(session->precision_bits);
        ComplexBall::set_precision(session->precision_bits);
        auto saturation = native_acb_saturation_binding(
            plan, checkpoint_configuration_identity(*session), arm_name,
            match_index, match_checkpoint, true);
        match_request[saturation.request_key] = std::move(saturation.request);
        match_request["refinement"] = refinement;
        match = build_planned_match_hop(
            match_handle, match_request, session->domain,
            session->precision_bits,
            checkpoint_configuration_identity(*session), plan,
            basis_handles, basis, incoming_handle, incoming, true);
        next = match->materialize(
            local_handle,
            arm_checkpoint_identity(checkpoint_root, arm_name, "local",
                                    match_index + 1),
            session->precision_bits, match);
      } else {
        match = build_planned_match_hop(
            match_handle, match_request, session->domain,
            session->precision_bits,
            checkpoint_configuration_identity(*session), plan,
            basis_handles, basis, incoming_handle, incoming, true);
        next = match->materialize(
            local_handle,
            arm_checkpoint_identity(checkpoint_root, arm_name, "local",
                                    match_index + 1),
            session->precision_bits, match);
      }
      sealed_lineage = next->seal_plan_match_lineage();
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "consuming transport hop reservation accounting underflow");
      --session->pending_local_solves;
      throw;
    }

    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "consuming transport hop reservation accounting underflow");
      --session->pending_local_solves;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during consuming transport hop");
      for (std::size_t column = 0; column < basis_handles.size(); ++column) {
        const auto found = session->locals.find(basis_handles[column]);
        if (found == session->locals.end() ||
            found->second.get() != basis[column].get())
          throw std::invalid_argument(
              "consuming transport hop basis owner changed before publication");
      }
      for (const auto& handle : basis_handles)
        session->locals.erase(handle);
      if (!session->locals.emplace(local_handle, next).second)
        throw std::logic_error(
            "consuming transport hop local handle collision");
      ++session->total_local_matches;
      session->total_local_match_ms += match->elapsed_ms();
      plan->note_match_advance(arm_name);
    }

    json::object match_reference{
        {"index", match_index},
        {"checkpoint_identity", match->checkpoint_identity()},
        {"provenance_identity", match->provenance_identity()},
        {"planned_hop", match->handoff()},
        {"sealed_local_lineage", sealed_lineage}};
    auto next_summary = compact_transport_local_reference(next);
    next_summary["release_via"] = "local";
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"capability", "consuming-transport-hop-v1"},
        {"native_retained", true}, {"json_coefficients", 0},
        {"arm", arm_name}, {"match", match_index},
        {"next_local", std::move(next_summary)},
        {"basis_reference", std::move(basis_reference)},
        {"match_reference", std::move(match_reference)},
        {"consumed_basis_handles", raw_basis},
        {"elapsed_ms", std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - started).count()}};
  }

  if (operation == "transport.publish_consumed_states") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "tile_plan",
         "tile_plan_checkpoint_identity", "anchor",
         "anchor_checkpoint_identity", "epsilon", "refinement",
         "checkpoint_policy", "lower", "upper"},
        "native transport.publish_consumed_states request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "transport.publish_consumed_states requires Rational or Acb coefficients");
    const auto epsilon_contract = parse_whole_arm_epsilon_contract(
        root.at("epsilon"), "published consumed-state epsilon contract");
    const auto& refinement = as_object(
        root.at("refinement"),
        "published consumed-state refinement policy");
    require_exact_keys(refinement, {"relative_tolerance", "max_steps"},
                       "published consumed-state refinement policy");
    if (required_string(refinement, "relative_tolerance").empty() ||
        as_u32(refinement.at("max_steps"),
               "published consumed-state refinement steps") > 32)
      throw std::invalid_argument(
          "published consumed-state refinement policy is invalid");
    const auto& checkpoint_policy = as_object(
        root.at("checkpoint_policy"),
        "published consumed-state checkpoint policy");
    require_exact_keys(checkpoint_policy, {"schema", "root"},
                       "published consumed-state checkpoint policy");
    if (required_string(checkpoint_policy, "schema") !=
        "diffexp2-deterministic-arm-checkpoints-v1")
      throw std::invalid_argument(
          "unsupported published consumed-state checkpoint policy");
    const auto checkpoint_root =
        required_string(checkpoint_policy, "root");

    struct ConsumedStateInput {
      std::string arm;
      std::vector<std::string> tile_handles;
      std::vector<std::shared_ptr<StoredLocalBase>> tile_sources;
    };
    const auto parse_side = [&](const char* arm) {
      const auto& side = as_object(
          root.at(arm), "published consumed transport-arm input");
      require_exact_keys(
          side, {"tile_sources"},
          "published consumed transport-arm input");
      ConsumedStateInput input;
      input.arm = arm;
      const auto& raw_tiles = as_array(
          side.at("tile_sources"),
          "published consumed transport-arm tile sources");
      if (raw_tiles.empty())
        throw std::invalid_argument(
            "published consumed transport arm has no tile sources");
      for (const auto& raw_handle : raw_tiles) {
        if (!raw_handle.is_string() || raw_handle.as_string().empty())
          throw std::invalid_argument(
              "published consumed transport tile-source handle must be nonempty");
        input.tile_handles.emplace_back(raw_handle.as_string());
      }
      return input;
    };
    std::array<ConsumedStateInput, 2> inputs{
        parse_side("lower"), parse_side("upper")};
    const auto plan_handle = required_string(root, "tile_plan");
    const auto anchor_handle = required_string(root, "anchor");
    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> anchor;
    std::array<std::string, 2> state_handles;
    std::set<std::string> unique_consumed_handles;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto plan_found = session->tile_plans.find(plan_handle);
      if (plan_found == session->tile_plans.end() ||
          required_string(root, "tile_plan_checkpoint_identity") !=
              plan_found->second->checkpoint_identity())
        throw std::invalid_argument(
            "published consumed-state tile-plan binding is stale");
      plan = plan_found->second;
      const auto anchor_found = session->locals.find(anchor_handle);
      if (anchor_found == session->locals.end() ||
          required_string(root, "anchor_checkpoint_identity") !=
              anchor_found->second->checkpoint_identity())
        throw std::invalid_argument(
            "published consumed-state anchor binding is stale");
      anchor = anchor_found->second;
      for (auto& input : inputs) {
        const auto& retained = plan->arm(input.arm);
        if (input.tile_handles.size() != retained.exact.tiles.size() ||
            input.tile_handles.front() != anchor_handle)
          throw std::invalid_argument(
              "published consumed-state input does not reproduce its exact arm topology");
        input.tile_sources.reserve(input.tile_handles.size());
        for (std::size_t tile = 0; tile < input.tile_handles.size(); ++tile) {
          const auto found = session->locals.find(input.tile_handles[tile]);
          if (found == session->locals.end())
            throw std::invalid_argument(
                "unknown or already consumed transport tile-source local: " +
                input.tile_handles[tile]);
          if (tile == 0 && found->second.get() != anchor.get())
            throw std::invalid_argument(
                "published consumed-state arm changed its common anchor owner");
          if (tile != 0 &&
              !unique_consumed_handles.insert(input.tile_handles[tile]).second)
            throw std::invalid_argument(
                "published consumed states require distinct non-anchor tile locals");
          input.tile_sources.push_back(found->second);
        }
      }
      if (session->transport_states.size() +
              session->pending_transport_states + 2 >
          session->transport_state_capacity)
        throw std::invalid_argument(
            "persistent transport-state capacity is exhausted");
      for (auto& handle : state_handles)
        handle = "transport:" +
            std::to_string(session->next_transport_state++);
      session->pending_transport_states += 2;
    }

    const auto started = std::chrono::steady_clock::now();
    std::array<std::shared_ptr<StoredTransportArmState>, 2> states;
    try {
      for (std::size_t index = 0; index < inputs.size(); ++index) {
        auto& input = inputs[index];
        states[index] = std::make_shared<StoredTransportArmState>(
            state_handles[index],
            checkpoint_root + ":" + input.arm + ":state", input.arm,
            plan, anchor, std::move(input.tile_sources), epsilon_contract.work,
            epsilon_contract.public_required_complete_max,
            epsilon_contract.match_required_complete_max, refinement, 0.0);
      }
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_transport_states < 2)
        throw std::logic_error(
            "published consumed-state reservation accounting underflow");
      session->pending_transport_states -= 2;
      throw;
    }

    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_transport_states < 2)
        throw std::logic_error(
            "published consumed-state reservation accounting underflow");
      session->pending_transport_states -= 2;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during consumed-state publication");
      for (const auto& handle : unique_consumed_handles)
        if (session->locals.find(handle) == session->locals.end())
          throw std::invalid_argument(
              "published consumed-state tile local changed before publication");
      std::size_t inserted = 0;
      try {
        for (; inserted < states.size(); ++inserted)
          if (!session->transport_states.emplace(
                  state_handles[inserted], states[inserted]).second)
            throw std::logic_error(
                "published consumed-state handle collision");
      } catch (...) {
        for (std::size_t index = 0; index < inserted; ++index)
          session->transport_states.erase(state_handles[index]);
        throw;
      }
      for (const auto& handle : unique_consumed_handles)
        session->locals.erase(handle);
      session->total_transport_arm_marches += 2;
    }

    json::object response_states;
    for (std::size_t index = 0; index < states.size(); ++index) {
      auto summary = states[index]->summary();
      summary["session"] = session->handle;
      response_states[inputs[index].arm] = std::move(summary);
    }
    json::array consumed_handles;
    for (const auto& handle : unique_consumed_handles)
      consumed_handles.push_back(json::value(handle));
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"capability", "published-consumed-transport-arm-states-v1"},
        {"native_retained", true}, {"json_coefficients", 0},
        {"atomic_publication", true},
        {"consumed_tile_local_handles", std::move(consumed_handles)},
        {"states", std::move(response_states)},
        {"elapsed_ms", std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - started).count()}};
  }

  if (operation == "transport.run_arms_consuming") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "tile_plan", "anchor",
         "tile_plan_checkpoint_identity", "anchor_checkpoint_identity",
         "epsilon", "refinement", "checkpoint_policy", "lower", "upper"},
        "native consuming transport request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "consuming transport requires Rational or Acb coefficients");
    const auto epsilon_contract = parse_whole_arm_epsilon_contract(
        root.at("epsilon"), "consuming transport epsilon contract");
    const auto& refinement = as_object(
        root.at("refinement"), "consuming transport refinement policy");
    require_exact_keys(refinement, {"relative_tolerance", "max_steps"},
                       "consuming transport refinement policy");
    (void)required_string(refinement, "relative_tolerance");
    if (as_u32(refinement.at("max_steps"),
               "consuming transport refinement steps") > 32)
      throw std::invalid_argument(
          "consuming transport refinement steps must lie in 0..32");
    const auto& checkpoint_policy = as_object(
        root.at("checkpoint_policy"),
        "consuming transport checkpoint policy");
    require_exact_keys(checkpoint_policy, {"schema", "root"},
                       "consuming transport checkpoint policy");
    if (required_string(checkpoint_policy, "schema") !=
        "diffexp2-deterministic-arm-checkpoints-v1")
      throw std::invalid_argument(
          "unsupported consuming transport checkpoint policy");
    const auto checkpoint_root =
        required_string(checkpoint_policy, "root");

    const auto parse_consuming_arm = [&](const char* name) {
      const auto& raw_arm = as_object(
          root.at(name), "consuming transport arm request");
      require_exact_keys(raw_arm, {"receiving_basis"},
                         "consuming transport arm request");
      RetainedArmMarchInput input;
      input.name = name;
      for (const auto& raw_set : as_array(
               raw_arm.at("receiving_basis"),
               "consuming transport receiving bases")) {
        const auto& values = as_array(
            raw_set, "consuming transport receiving basis");
        if (values.empty())
          throw std::invalid_argument(
              "consuming transport basis cannot be empty");
        std::set<std::string> unique;
        std::vector<std::string> handles;
        for (const auto& raw_handle : values) {
          if (!raw_handle.is_string() || raw_handle.as_string().empty())
            throw std::invalid_argument(
                "consuming transport basis handle must be nonempty");
          std::string handle(raw_handle.as_string());
          if (!unique.insert(handle).second)
            throw std::invalid_argument(
                "one consuming transport basis contains a duplicate handle");
          handles.push_back(std::move(handle));
        }
        input.basis_handles.push_back(std::move(handles));
      }
      return input;
    };
    std::array<RetainedArmMarchInput, 2> inputs{
        parse_consuming_arm("lower"),
        parse_consuming_arm("upper")};
    const auto plan_handle = required_string(root, "tile_plan");
    const auto anchor_handle = required_string(root, "anchor");
    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> anchor;
    std::array<std::string, 2> state_handles;
    std::unordered_map<std::string, std::size_t> remaining_basis_uses;
    const auto total_matches = checked_diagnostic_sum(
        inputs[0].basis_handles.size(), inputs[1].basis_handles.size(),
        "consuming transport match count");
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto plan_found = session->tile_plans.find(plan_handle);
      if (plan_found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown tile plan for consuming transport");
      plan = plan_found->second;
      if (required_string(root, "tile_plan_checkpoint_identity") !=
          plan->checkpoint_identity())
        throw std::invalid_argument(
            "consuming transport tile-plan checkpoint is stale");
      const auto anchor_found = session->locals.find(anchor_handle);
      if (anchor_found == session->locals.end())
        throw std::invalid_argument(
            "unknown anchor for consuming transport");
      anchor = anchor_found->second;
      if (required_string(root, "anchor_checkpoint_identity") !=
          anchor->checkpoint_identity())
        throw std::invalid_argument(
            "consuming transport anchor checkpoint is stale");
      for (auto& input : inputs) {
        const auto& retained = plan->arm(input.name);
        if (retained.exact.matches.size() != input.basis_handles.size() ||
            retained.exact.tiles.size() != input.basis_handles.size() + 1)
          throw std::invalid_argument(
              "consuming transport bases do not reproduce plan topology");
        for (const auto& handles : input.basis_handles)
          for (const auto& handle : handles) {
            if (session->locals.find(handle) == session->locals.end())
              throw std::invalid_argument(
                  "unknown receiving basis for consuming transport: " +
                  handle);
            ++remaining_basis_uses[handle];
          }
        input.match_handles.reserve(input.basis_handles.size());
        input.local_handles.reserve(input.basis_handles.size());
        for (std::size_t i = 0; i < input.basis_handles.size(); ++i) {
          input.match_handles.push_back(
              "m:" + std::to_string(session->next_match++));
          input.local_handles.push_back(
              "l:" + std::to_string(session->next_local++));
        }
      }
      for (auto& handle : state_handles)
        handle = "transport:" +
            std::to_string(session->next_transport_state++);
      if (2 > session->transport_state_capacity -
                  std::min(session->transport_state_capacity,
                           session->transport_states.size()))
        throw std::invalid_argument(
            "transport-state capacity exhausted by consuming march");
    }

    const auto started = std::chrono::steady_clock::now();
    std::array<ConsumingArmMarchResult, 2> marched;
    std::array<std::shared_ptr<StoredTransportArmState>, 2> states;
    std::unique_ptr<AcbPrecisionLease> acb_lease;
    if (session->domain == "acb") {
      acb_lease = std::make_unique<AcbPrecisionLease>(
          session->precision_bits);
      ComplexBall::set_precision(session->precision_bits);
    }
    const auto session_configuration =
        checkpoint_configuration_identity(*session);
    // Sequential by construction: a completed hop can reclaim its receiving
    // basis before the next hop allocates a materialized local.
    for (std::size_t arm_index = 0; arm_index < inputs.size(); ++arm_index) {
      marched[arm_index] = march_retained_arm_consuming(
          session, session_configuration, plan, anchor, inputs[arm_index],
          epsilon_contract.work,
          epsilon_contract.match_required_complete_max, refinement,
          checkpoint_root, remaining_basis_uses);
      states[arm_index] = std::make_shared<StoredTransportArmState>(
          state_handles[arm_index],
          checkpoint_root + ":" + inputs[arm_index].name + ":state",
          inputs[arm_index].name, plan, anchor,
          std::move(marched[arm_index].basis_references),
          std::move(marched[arm_index].match_references),
          std::move(marched[arm_index].tile_sources),
          epsilon_contract.work,
          epsilon_contract.public_required_complete_max,
          epsilon_contract.match_required_complete_max, refinement,
          marched[arm_index].elapsed_ms);
    }
    if (std::any_of(remaining_basis_uses.begin(),
                    remaining_basis_uses.end(),
                    [](const auto& item) { return item.second != 0; }))
      throw std::logic_error(
          "consuming transport left nonzero basis last-use counts");

    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument(
            "persistent session closed before consuming-state publication");
      for (std::size_t i = 0; i < states.size(); ++i)
        if (!session->transport_states.emplace(
                state_handles[i], states[i]).second)
          throw std::logic_error(
              "consuming transport-state handle collision");
      session->total_local_matches +=
          static_cast<std::uint64_t>(total_matches);
      session->total_transport_arm_marches += 2;
      plan->note_two_arm_match_advances(
          inputs[0].basis_handles.size(),
          inputs[1].basis_handles.size());
    }
    json::object response_states;
    for (std::size_t i = 0; i < states.size(); ++i) {
      auto summary = states[i]->summary();
      summary["session"] = session->handle;
      response_states[inputs[i].name] = std::move(summary);
    }
    json::array consumption;
    for (auto& result : marched)
      for (auto& diagnostic : result.release_diagnostics)
        consumption.push_back(std::move(diagnostic));
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"capability", "consuming-sequential-transport-arm-state-v1"},
        {"native_retained", true}, {"json_coefficients", 0},
        {"atomic_publication", true},
        {"consuming_basis_handles", true},
        {"workers", 1}, {"max_parallel_arms", 1},
        {"consumption", std::move(consumption)},
        {"marches", 2}, {"states", std::move(response_states)},
        {"elapsed_ms", std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - started).count()}};
  }

  if (operation == "transport.run_arms") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "tile_plan", "anchor",
         "tile_plan_checkpoint_identity", "anchor_checkpoint_identity",
         "epsilon", "refinement", "checkpoint_policy", "lower", "upper"},
        "native transport.run_arms request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native two-arm transport marching requires rational or Acb coefficients");
    const auto epsilon_contract = parse_whole_arm_epsilon_contract(
        root.at("epsilon"), "native two-arm transport epsilon contract");
    const auto& refinement = as_object(
        root.at("refinement"),
        "native two-arm transport refinement policy");
    require_exact_keys(refinement, {"relative_tolerance", "max_steps"},
                       "native two-arm transport refinement policy");
    (void)required_string(refinement, "relative_tolerance");
    if (as_u32(refinement.at("max_steps"),
               "native two-arm transport refinement steps") > 32)
      throw std::invalid_argument(
          "native two-arm transport refinement steps must lie in 0..32");
    const auto& checkpoint_policy = as_object(
        root.at("checkpoint_policy"),
        "native two-arm transport checkpoint policy");
    require_exact_keys(checkpoint_policy, {"schema", "root"},
                       "native two-arm transport checkpoint policy");
    if (required_string(checkpoint_policy, "schema") !=
        "diffexp2-deterministic-arm-checkpoints-v1")
      throw std::invalid_argument(
          "unsupported native two-arm transport checkpoint policy schema");
    const auto checkpoint_root = required_string(checkpoint_policy, "root");
    if (checkpoint_root.empty())
      throw std::invalid_argument(
          "native two-arm transport checkpoint root cannot be empty");

    const auto parse_arm = [&](const char* name) {
      const auto& raw_arm = as_object(
          root.at(name), "native two-arm transport arm request");
      require_exact_keys(raw_arm, {"receiving_basis"},
                         "native two-arm transport arm request");
      RetainedArmMarchInput input;
      input.name = name;
      const auto& raw_basis_sets = as_array(
          raw_arm.at("receiving_basis"),
          "native two-arm transport receiving basis sets");
      input.basis_handles.reserve(raw_basis_sets.size());
      for (const auto& raw_set : raw_basis_sets) {
        const auto& values = as_array(
            raw_set, "native two-arm transport receiving basis set");
        if (values.empty())
          throw std::invalid_argument(
              "native two-arm transport receiving basis sets cannot be empty");
        std::set<std::string> unique;
        std::vector<std::string> handles;
        handles.reserve(values.size());
        for (const auto& raw_handle : values) {
          if (!raw_handle.is_string() || raw_handle.as_string().empty())
            throw std::invalid_argument(
                "native two-arm transport basis handles must be nonempty strings");
          std::string handle(raw_handle.as_string());
          if (!unique.insert(handle).second)
            throw std::invalid_argument(
                "native two-arm transport basis handles must be pairwise distinct");
          handles.push_back(std::move(handle));
        }
        input.basis_handles.push_back(std::move(handles));
      }
      return input;
    };
    std::array<RetainedArmMarchInput, 2> inputs{
        parse_arm("lower"), parse_arm("upper")};

    const auto plan_handle = required_string(root, "tile_plan");
    const auto anchor_handle = required_string(root, "anchor");
    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> anchor;
    std::array<std::string, 2> state_handles;
    const auto total_matches = checked_diagnostic_sum(
        inputs[0].basis_handles.size(), inputs[1].basis_handles.size(),
        "two-arm transport match count");
    constexpr std::size_t published_transport_states = 2;
    bool reservation_live = false;
    {
      // Resolve both complete owner sets and reserve both marches before a
      // worker starts. Public release calls cannot invalidate these owners.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto plan_found = session->tile_plans.find(plan_handle);
      if (plan_found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or released tile plan for two-arm transport marching");
      plan = plan_found->second;
      if (required_string(root, "tile_plan_checkpoint_identity") !=
          plan->checkpoint_identity())
        throw std::invalid_argument(
            "two-arm transport tile-plan checkpoint token is stale");
      if (!plan->has_two_arms())
        throw std::invalid_argument(
            "two-arm transport marching requires a retained two-arm plan");
      const auto anchor_found = session->locals.find(anchor_handle);
      if (anchor_found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released anchor local for two-arm transport marching");
      anchor = anchor_found->second;
      if (required_string(root, "anchor_checkpoint_identity") !=
          anchor->checkpoint_identity())
        throw std::invalid_argument(
            "two-arm transport anchor checkpoint token is stale");

      for (auto& input : inputs) {
        const auto& retained = plan->arm(input.name);
        if (retained.exact.matches.size() != input.basis_handles.size() ||
            retained.exact.tiles.size() != input.basis_handles.size() + 1)
          throw std::invalid_argument(
              "two-arm transport basis count does not reproduce the retained plan topology for " +
              input.name);
        const auto& anchor_binding = retained.charts.front();
        if (anchor->source_chart() != anchor_binding.handle)
          throw std::invalid_argument(
              "two-arm transport anchor belongs to a different retained chart for " +
              input.name);
        anchor->require_exact_plan_binding(
            anchor_binding.geometry, anchor_binding.prescriptions,
            "two-arm transport anchor");

        input.basis.reserve(input.basis_handles.size());
        for (const auto& handles : input.basis_handles) {
          std::vector<std::shared_ptr<StoredLocalBase>> resolved;
          resolved.reserve(handles.size());
          for (const auto& handle : handles) {
            const auto found = session->locals.find(handle);
            if (found == session->locals.end())
              throw std::invalid_argument(
                  "unknown or released native local in two-arm transport receiving basis: " +
                  handle);
            resolved.push_back(found->second);
          }
          input.basis.push_back(std::move(resolved));
        }
      }

      if (total_matches > session->match_capacity -
                              std::min(session->match_capacity,
                                       session->matches.size() +
                                           session->pending_matches))
        throw std::invalid_argument(
            "persistent local match capacity is exhausted by two-arm transport marching");
      if (total_matches > session->local_capacity -
                              std::min(session->local_capacity,
                                       session->locals.size() +
                                           session->pending_local_solves))
        throw std::invalid_argument(
            "persistent local capacity is exhausted by two-arm transport marching");
      if (published_transport_states >
          session->transport_state_capacity -
              std::min(session->transport_state_capacity,
                       session->transport_states.size() +
                           session->pending_transport_states))
        throw std::invalid_argument(
            "persistent transport-state capacity is exhausted by two-arm transport marching");

      for (std::size_t arm_index = 0; arm_index < inputs.size(); ++arm_index) {
        auto& input = inputs[arm_index];
        input.match_handles.reserve(input.basis_handles.size());
        input.local_handles.reserve(input.basis_handles.size());
        for (std::size_t index = 0; index < input.basis_handles.size();
             ++index) {
          input.match_handles.push_back(
              "m:" + std::to_string(session->next_match++));
          input.local_handles.push_back(
              "l:" + std::to_string(session->next_local++));
        }
        state_handles[arm_index] = "transport:" +
            std::to_string(session->next_transport_state++);
      }
      session->pending_matches += total_matches;
      session->pending_local_solves += total_matches;
      session->pending_transport_states += published_transport_states;
      reservation_live = true;
    }
    const auto release_reservation = [&]() {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (!reservation_live) return;
      if (session->pending_matches < total_matches ||
          session->pending_local_solves < total_matches ||
          session->pending_transport_states < published_transport_states)
        throw std::logic_error(
            "native two-arm transport reservation accounting underflow");
      session->pending_matches -= total_matches;
      session->pending_local_solves -= total_matches;
      session->pending_transport_states -= published_transport_states;
      reservation_live = false;
    };
    struct TwoArmTransportReservationGuard final {
      decltype(release_reservation)& release;
      bool& live;
      ~TwoArmTransportReservationGuard() noexcept {
        if (!live) return;
        try {
          release();
        } catch (...) {
          std::terminate();
        }
      }
    } reservation_guard{release_reservation, reservation_live};

    struct CompletedTransportArm {
      RetainedArmMarchResult march;
      std::shared_ptr<StoredTransportArmState> state;
    };
    std::array<CompletedTransportArm, 2> completed;
    std::array<std::exception_ptr, 2> failures;
    std::atomic<std::size_t> active_workers{0};
    std::atomic<std::size_t> max_active_workers{0};
    std::mutex start_mutex;
    std::condition_variable start_changed;
    std::size_t workers_ready = 0;
    bool workers_start = false;
    bool workers_cancel = false;
    const auto operation_started = std::chrono::steady_clock::now();
    std::unique_ptr<AcbPrecisionLease> acb_lease;
    if (session->domain == "acb") {
      acb_lease = std::make_unique<AcbPrecisionLease>(
          session->precision_bits);
      ComplexBall::set_precision(session->precision_bits);
    }
    const auto session_configuration =
        checkpoint_configuration_identity(*session);
    const auto update_max_active = [&](std::size_t candidate) {
      auto observed = max_active_workers.load();
      while (observed < candidate &&
             !max_active_workers.compare_exchange_weak(observed, candidate)) {
      }
    };
    const auto run_arm = [&](std::size_t arm_index) {
      const auto active = active_workers.fetch_add(1) + 1;
      update_max_active(active);
      {
        std::unique_lock<std::mutex> lock(start_mutex);
        ++workers_ready;
        start_changed.notify_all();
        start_changed.wait(lock, [&] { return workers_start; });
        if (workers_cancel) {
          active_workers.fetch_sub(1);
          return;
        }
      }
      try {
        if (session->domain == "acb")
          ComplexBall::set_precision(session->precision_bits);
        auto& input = inputs[arm_index];
        auto& output = completed[arm_index];
        output.march = march_retained_arm(
            session->domain, session->precision_bits,
            session_configuration, plan, anchor, input,
            epsilon_contract.work,
            epsilon_contract.match_required_complete_max, refinement,
            checkpoint_root);
        output.state = std::make_shared<StoredTransportArmState>(
            state_handles[arm_index],
            checkpoint_root + ":" + input.name + ":state", input.name,
            plan, anchor, input.basis, output.march.matches,
            output.march.tile_sources, epsilon_contract.work,
            epsilon_contract.public_required_complete_max,
            epsilon_contract.match_required_complete_max, refinement,
            output.march.elapsed_ms);
      } catch (...) {
        failures[arm_index] = std::current_exception();
      }
      active_workers.fetch_sub(1);
    };

    std::vector<std::jthread> workers;
    workers.reserve(2);
    try {
      workers.emplace_back([&] { run_arm(0); });
      workers.emplace_back([&] { run_arm(1); });
    } catch (...) {
      {
        std::lock_guard<std::mutex> lock(start_mutex);
        workers_cancel = true;
        workers_start = true;
      }
      start_changed.notify_all();
      for (auto& worker : workers)
        if (worker.joinable()) worker.join();
      release_reservation();
      throw;
    }
    {
      std::unique_lock<std::mutex> lock(start_mutex);
      start_changed.wait(lock, [&] { return workers_ready == 2; });
      workers_start = true;
    }
    start_changed.notify_all();
    for (auto& worker : workers) worker.join();

    if (failures[0] || failures[1]) {
      release_reservation();
      // Deterministically prefer the lower failure when both workers fail.
      std::rethrow_exception(failures[0] ? failures[0] : failures[1]);
    }

    double completed_match_ms = 0.0;
    for (const auto& output : completed)
      for (const auto& match : output.state->matches())
        completed_match_ms += match->elapsed_ms();
    const double completed_arm_ms =
        completed[0].state->elapsed_ms() +
        completed[1].state->elapsed_ms();

    try {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches < total_matches ||
          session->pending_local_solves < total_matches ||
          session->pending_transport_states < published_transport_states)
        throw std::logic_error(
            "native two-arm transport reservation accounting underflow");
      plan->require_two_arm_match_advance_capacity(
          completed[0].state->matches().size(),
          completed[1].state->matches().size());
      if (total_matches >
              std::numeric_limits<std::uint64_t>::max() -
                  session->total_local_matches ||
          session->total_transport_arm_marches >
              std::numeric_limits<std::uint64_t>::max() - 2)
        throw std::overflow_error(
            "native two-arm transport session counter overflow");
      session->pending_matches -= total_matches;
      session->pending_local_solves -= total_matches;
      session->pending_transport_states -= published_transport_states;
      reservation_live = false;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during two-arm transport marching");

      session->transport_states.reserve(
          session->transport_states.size() + published_transport_states);
      std::vector<std::string> inserted;
      inserted.reserve(published_transport_states);
      try {
        for (std::size_t arm_index = 0; arm_index < completed.size();
             ++arm_index) {
          if (!session->transport_states.emplace(
                  state_handles[arm_index], completed[arm_index].state).second)
            throw std::logic_error(
                "two-arm transport state handle collision at publication");
          inserted.push_back(state_handles[arm_index]);
        }
      } catch (...) {
        for (const auto& handle : inserted)
          session->transport_states.erase(handle);
        throw;
      }
      session->total_local_matches +=
          static_cast<std::uint64_t>(total_matches);
      session->total_local_match_ms += completed_match_ms;
      plan->note_two_arm_match_advances(
          completed[0].state->matches().size(),
          completed[1].state->matches().size());
      session->total_transport_arm_marches += 2;
      session->total_transport_arm_ms += completed_arm_ms;
    } catch (...) {
      if (reservation_live) release_reservation();
      throw;
    }

    json::object state_response;
    for (std::size_t arm_index = 0; arm_index < completed.size();
         ++arm_index) {
      auto summary = completed[arm_index].state->summary();
      summary["session"] = session->handle;
      auto& final_local = summary.at("final_local").as_object();
      final_local["dependency_only"] = true;
      final_local["public_token"] = false;
      final_local["release_via"] = "transport_state";
      state_response[inputs[arm_index].name] = std::move(summary);
    }
    return json::object{
        {"status", "ok"},
        {"session", session->handle},
        {"capability", kRetainedParallelTransportArmStateCapability},
        {"native_retained", true},
        {"json_coefficients", 0},
        {"atomic_publication", true},
        {"public_result_tokens", "transport_states_only"},
        {"dependency_only_final_locals", true},
        {"workers", 2},
        {"max_parallel_arms", max_active_workers.load()},
        {"worker_overlap", max_active_workers.load() == 2},
        {"checkpoint_policy", checkpoint_policy},
        {"epsilon", root.at("epsilon")},
        {"matches", json::object{
             {"lower", completed[0].march.matches.size()},
             {"upper", completed[1].march.matches.size()},
             {"total", total_matches}}},
        {"marches", 2},
        {"plan_stats", plan->summary(false)},
        {"states", std::move(state_response)},
        {"elapsed_ms", std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - operation_started).count()}};
  }

  if (operation == "transport.run_arm") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "tile_plan", "anchor",
         "tile_plan_checkpoint_identity", "anchor_checkpoint_identity",
         "arm", "receiving_basis", "epsilon", "refinement",
         "checkpoint_policy"},
        "native transport.run_arm request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native transport-arm marching requires rational or Acb coefficients");
    const auto arm_name = required_string(root, "arm");
    if (arm_name != "lower" && arm_name != "upper")
      throw std::invalid_argument(
          "native transport-arm name must be lower or upper");
    const auto epsilon_contract = parse_whole_arm_epsilon_contract(
        root.at("epsilon"), "native transport-arm epsilon contract");
    const auto& refinement = as_object(
        root.at("refinement"), "native transport-arm refinement policy");
    require_exact_keys(refinement, {"relative_tolerance", "max_steps"},
                       "native transport-arm refinement policy");
    (void)required_string(refinement, "relative_tolerance");
    if (as_u32(refinement.at("max_steps"),
               "native transport-arm refinement steps") > 32)
      throw std::invalid_argument(
          "native transport-arm refinement steps must lie in 0..32");
    const auto& checkpoint_policy = as_object(
        root.at("checkpoint_policy"),
        "native transport-arm checkpoint policy");
    require_exact_keys(checkpoint_policy, {"schema", "root"},
                       "native transport-arm checkpoint policy");
    if (required_string(checkpoint_policy, "schema") !=
        "diffexp2-deterministic-arm-checkpoints-v1")
      throw std::invalid_argument(
          "unsupported native transport-arm checkpoint policy schema");
    const auto checkpoint_root = required_string(checkpoint_policy, "root");
    if (checkpoint_root.empty())
      throw std::invalid_argument(
          "native transport-arm checkpoint root cannot be empty");

    RetainedArmMarchInput input;
    input.name = arm_name;
    const auto& raw_basis_sets = as_array(
        root.at("receiving_basis"),
        "native transport-arm receiving basis sets");
    input.basis_handles.reserve(raw_basis_sets.size());
    for (const auto& raw_set : raw_basis_sets) {
      const auto& values = as_array(
          raw_set, "native transport-arm receiving basis set");
      if (values.empty())
        throw std::invalid_argument(
            "native transport-arm receiving basis sets cannot be empty");
      std::set<std::string> unique;
      std::vector<std::string> handles;
      handles.reserve(values.size());
      for (const auto& raw_handle : values) {
        if (!raw_handle.is_string() || raw_handle.as_string().empty())
          throw std::invalid_argument(
              "native transport-arm basis handles must be nonempty strings");
        std::string handle(raw_handle.as_string());
        if (!unique.insert(handle).second)
          throw std::invalid_argument(
              "native transport-arm basis handles must be pairwise distinct");
        handles.push_back(std::move(handle));
      }
      input.basis_handles.push_back(std::move(handles));
    }

    const auto plan_handle = required_string(root, "tile_plan");
    const auto anchor_handle = required_string(root, "anchor");
    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> anchor;
    std::string state_handle;
    const auto match_count = input.basis_handles.size();
    bool reservation_live = false;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto plan_found = session->tile_plans.find(plan_handle);
      if (plan_found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or released tile plan for transport-arm marching");
      plan = plan_found->second;
      if (required_string(root, "tile_plan_checkpoint_identity") !=
          plan->checkpoint_identity())
        throw std::invalid_argument(
            "transport-arm tile-plan checkpoint token is stale");
      const auto& retained = plan->arm(arm_name);
      if (retained.exact.matches.size() != match_count ||
          retained.exact.tiles.size() != match_count + 1)
        throw std::invalid_argument(
            "transport-arm basis count does not reproduce the retained plan topology");
      const auto anchor_found = session->locals.find(anchor_handle);
      if (anchor_found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released anchor local for transport-arm marching");
      anchor = anchor_found->second;
      if (required_string(root, "anchor_checkpoint_identity") !=
          anchor->checkpoint_identity())
        throw std::invalid_argument(
            "transport-arm anchor checkpoint token is stale");
      const auto& anchor_binding = retained.charts.front();
      if (anchor->source_chart() != anchor_binding.handle)
        throw std::invalid_argument(
            "transport-arm anchor belongs to a different retained chart");
      anchor->require_exact_plan_binding(
          anchor_binding.geometry, anchor_binding.prescriptions,
          "transport-arm anchor");

      input.basis.reserve(match_count);
      for (const auto& handles : input.basis_handles) {
        std::vector<std::shared_ptr<StoredLocalBase>> resolved;
        resolved.reserve(handles.size());
        for (const auto& handle : handles) {
          const auto found = session->locals.find(handle);
          if (found == session->locals.end())
            throw std::invalid_argument(
                "unknown or released native local in transport-arm receiving basis: " +
                handle);
          resolved.push_back(found->second);
        }
        input.basis.push_back(std::move(resolved));
      }

      if (match_count > session->match_capacity -
                            std::min(session->match_capacity,
                                     session->matches.size() +
                                         session->pending_matches))
        throw std::invalid_argument(
            "persistent local match capacity is exhausted by transport-arm marching");
      if (match_count > session->local_capacity -
                            std::min(session->local_capacity,
                                     session->locals.size() +
                                         session->pending_local_solves))
        throw std::invalid_argument(
            "persistent local capacity is exhausted by transport-arm marching");
      if (session->transport_states.size() +
              session->pending_transport_states >=
          session->transport_state_capacity)
        throw std::invalid_argument(
            "persistent transport-state capacity is exhausted");

      input.match_handles.reserve(match_count);
      input.local_handles.reserve(match_count);
      for (std::size_t index = 0; index < match_count; ++index) {
        input.match_handles.push_back(
            "m:" + std::to_string(session->next_match++));
        input.local_handles.push_back(
            "l:" + std::to_string(session->next_local++));
      }
      state_handle = "transport:" +
          std::to_string(session->next_transport_state++);
      session->pending_matches += match_count;
      session->pending_local_solves += match_count;
      ++session->pending_transport_states;
      reservation_live = true;
    }
    const auto release_reservation = [&]() {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (!reservation_live) return;
      if (session->pending_matches < match_count ||
          session->pending_local_solves < match_count ||
          session->pending_transport_states == 0)
        throw std::logic_error(
            "native transport-arm reservation accounting underflow");
      session->pending_matches -= match_count;
      session->pending_local_solves -= match_count;
      --session->pending_transport_states;
      reservation_live = false;
    };
    struct TransportArmReservationGuard final {
      decltype(release_reservation)& release;
      bool& live;
      ~TransportArmReservationGuard() noexcept {
        if (!live) return;
        try {
          release();
        } catch (...) {
          std::terminate();
        }
      }
    } reservation_guard{release_reservation, reservation_live};

    const auto operation_started = std::chrono::steady_clock::now();
    std::unique_ptr<AcbPrecisionLease> acb_lease;
    if (session->domain == "acb") {
      acb_lease = std::make_unique<AcbPrecisionLease>(
          session->precision_bits);
      ComplexBall::set_precision(session->precision_bits);
    }
    const auto session_configuration =
        checkpoint_configuration_identity(*session);
    auto marched = march_retained_arm(
        session->domain, session->precision_bits, session_configuration,
        plan, anchor, input, epsilon_contract.work,
        epsilon_contract.match_required_complete_max, refinement,
        checkpoint_root);
    auto state = std::make_shared<StoredTransportArmState>(
        state_handle, checkpoint_root + ":" + arm_name + ":state",
        arm_name, plan, anchor, input.basis, marched.matches,
        marched.tile_sources, epsilon_contract.work,
        epsilon_contract.public_required_complete_max,
        epsilon_contract.match_required_complete_max, refinement,
        marched.elapsed_ms);
    const auto final_local = state->final_local();

    try {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches < match_count ||
          session->pending_local_solves < match_count ||
          session->pending_transport_states == 0)
        throw std::logic_error(
            "native transport-arm reservation accounting underflow");
      session->pending_matches -= match_count;
      session->pending_local_solves -= match_count;
      --session->pending_transport_states;
      reservation_live = false;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during transport-arm marching");

      session->locals.reserve(session->locals.size() +
                              (match_count == 0 ? 0 : 1));
      session->transport_states.reserve(
          session->transport_states.size() + 1);
      bool inserted_local = false;
      try {
        const auto existing = session->locals.find(final_local->handle());
        if (existing == session->locals.end()) {
          if (!session->locals.emplace(final_local->handle(), final_local)
                   .second)
            throw std::logic_error(
                "transport-arm final-local handle collision at publication");
          inserted_local = true;
        } else if (existing->second.get() != final_local.get()) {
          throw std::logic_error(
              "transport-arm final-local handle names a different retained object");
        }
        if (!session->transport_states.emplace(state_handle, state).second)
          throw std::logic_error(
              "transport-arm state handle collision at publication");
      } catch (...) {
        session->transport_states.erase(state_handle);
        if (inserted_local) session->locals.erase(final_local->handle());
        throw;
      }
      for (const auto& match : state->matches()) {
        ++session->total_local_matches;
        session->total_local_match_ms += match->elapsed_ms();
        plan->note_match_advance(arm_name);
      }
      ++session->total_transport_arm_marches;
      session->total_transport_arm_ms += state->elapsed_ms();
    } catch (...) {
      if (reservation_live) release_reservation();
      throw;
    }

    auto response = state->summary();
    auto final_summary = final_local->summary();
    final_summary["session"] = session->handle;
    response["final_local"] = std::move(final_summary);
    response["status"] = "ok";
    response["session"] = session->handle;
    response["atomic_publication"] = true;
    response["checkpoint_policy"] = checkpoint_policy;
    response["operation_elapsed_ms"] = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - operation_started).count();
    return response;
  }

  if (operation == "transport.stats" ||
      operation == "transport.release") {
    require_exact_keys(root,
        {"schema", "op", "session", "transport_state"},
        operation == "transport.stats"
            ? "native transport.stats request"
            : "native transport.release request");
    const auto state_handle = required_string(root, "transport_state");
    if (operation == "transport.stats") {
      std::shared_ptr<StoredTransportArmState> state;
      {
        std::lock_guard<std::mutex> lock(session->mutex);
        const auto found = session->transport_states.find(state_handle);
        if (found == session->transport_states.end())
          throw std::invalid_argument(
              "unknown or released native transport-arm state");
        state = found->second;
      }
      auto result = state->stats_json();
      result["status"] = "ok";
      result["session"] = session->handle;
      return result;
    }
    std::shared_ptr<StoredTransportArmState> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->transport_states.find(state_handle);
      if (found == session->transport_states.end())
        throw std::invalid_argument(
            "unknown or already released native transport-arm state");
      removed = std::move(found->second);
      session->transport_states.erase(found);
    }
    return json::object{
        {"status", "ok"}, {"released", state_handle},
        {"checkpoint_identity", removed->checkpoint_identity()}};
  }

  if (operation == "transport.contract") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "transport_state",
         "transport_state_checkpoint_identity",
         "transport_state_provenance_identity", "checkpoint_policy",
         "observables"},
        "native transport.contract request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native transport contraction requires rational or Acb coefficients");
    const auto state_handle = required_string(root, "transport_state");
    const auto& checkpoint_policy = as_object(
        root.at("checkpoint_policy"),
        "native transport contraction checkpoint policy");
    require_exact_keys(checkpoint_policy, {"schema", "root"},
                       "native transport contraction checkpoint policy");
    if (required_string(checkpoint_policy, "schema") !=
        "diffexp2-deterministic-transport-contraction-checkpoints-v1")
      throw std::invalid_argument(
          "unsupported native transport contraction checkpoint policy schema");
    const auto checkpoint_root = required_string(checkpoint_policy, "root");
    if (checkpoint_root.empty())
      throw std::invalid_argument(
          "native transport contraction checkpoint root cannot be empty");

    struct PendingObservable {
      std::string identity;
      std::string checkpoint_identity;
      std::vector<json::object> rows;
      ObservableEpsilonContract epsilon;
      json::object epsilon_record;
      TransportTailPolicy tail_policy = TransportTailPolicy::Stored;
    };
    std::vector<PendingObservable> pending_observables;
    const auto& raw_observables = as_array(
        root.at("observables"), "native transport observables");
    pending_observables.reserve(raw_observables.size());
    std::set<std::string> observable_identities;
    std::set<std::string> observable_checkpoints;
    for (const auto& raw_observable : raw_observables) {
      const auto& observable = as_object(
          raw_observable, "native transport observable");
      require_exact_keys(
          observable,
          {"identity", "checkpoint_identity", "integrand_rows", "epsilon",
           "tail_policy"},
          "native transport observable");
      PendingObservable parsed;
      parsed.identity = required_string(observable, "identity");
      parsed.checkpoint_identity = required_string(
          observable, "checkpoint_identity");
      if (parsed.identity.empty() || parsed.checkpoint_identity.empty() ||
          !observable_identities.insert(parsed.identity).second ||
          !observable_checkpoints.insert(parsed.checkpoint_identity).second)
        throw std::invalid_argument(
            "native transport observable identities and checkpoints must be nonempty and pairwise unique");
      const auto& raw_rows = as_array(
          observable.at("integrand_rows"),
          "native transport observable rows");
      parsed.rows.reserve(raw_rows.size());
      for (const auto& raw_row : raw_rows)
        parsed.rows.push_back(as_object(
            raw_row, "native transport observable rational row"));
      parsed.epsilon = parse_observable_epsilon_contract(
          observable.at("epsilon"),
          "native transport observable output epsilon contract");
      parsed.epsilon_record = as_object(
          observable.at("epsilon"),
          "native transport observable output epsilon contract");
      parsed.tail_policy = parse_transport_tail_policy(
          observable.at("tail_policy"),
          "native transport observable tail policy");
      pending_observables.push_back(std::move(parsed));
    }

    std::shared_ptr<StoredTransportArmState> state;
    std::vector<TransportObservableContractionInput> contraction_inputs;
    bool reservation_live = false;
    const auto observable_count = pending_observables.size();
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto found = session->transport_states.find(state_handle);
      if (found == session->transport_states.end())
        throw std::invalid_argument(
            "unknown or released native transport-arm state for contraction");
      state = found->second;
      if (required_string(root, "transport_state_checkpoint_identity") !=
              state->checkpoint_identity() ||
          required_string(root, "transport_state_provenance_identity") !=
              state->provenance_identity())
        throw std::invalid_argument(
            "native transport contraction state binding is stale");
      state->require_contraction_counter_capacity(observable_count);
      const auto tile_count = state->tile_sources().size();
      for (const auto& observable : pending_observables) {
        if (observable.rows.size() != tile_count)
          throw std::invalid_argument(
              "native transport observable must provide one prepared rational row per retained tile");
        if (observable.epsilon.required_complete_max >
            state->public_required_complete_max())
          throw std::invalid_argument(
              "native transport observable required epsilon maximum exceeds the state public target");
      }
      if (observable_count > session->line_result_capacity -
                                 std::min(
                                     session->line_result_capacity,
                                     session->line_results.size() +
                                         session->pending_line_integrations))
        throw std::invalid_argument(
            "persistent line-result capacity is exhausted by transport contraction");

      contraction_inputs.reserve(observable_count);
      for (std::size_t index = 0; index < observable_count; ++index) {
        const auto& observable = pending_observables[index];
        TransportObservableContractionInput input;
        input.identity = observable.identity;
        input.checkpoint_identity = observable.checkpoint_identity;
        input.checkpoint_root = checkpoint_root + ":observable:" +
                                std::to_string(index + 1);
        input.rows = observable.rows;
        input.epsilon = observable.epsilon;
        input.epsilon_record = observable.epsilon_record;
        input.tail_policy = observable.tail_policy;
        input.aggregate_handle =
            "line:" + std::to_string(session->next_line_result++);
        input.projected_local_handles.reserve(tile_count);
        json::array row_records;
        row_records.reserve(tile_count);
        for (std::size_t tile = 0; tile < tile_count; ++tile) {
          input.projected_local_handles.push_back(
              "private:" + arm_checkpoint_identity(
                               input.checkpoint_root, state->arm_name(),
                               "projected", tile + 1));
          const auto row_identity = required_string(
              observable.rows[tile], "exact_identity");
          row_records.push_back(json::object{
              {"tile", tile}, {"exact_identity", row_identity},
              {"prepared_row", observable.rows[tile]}});
        }
        input.aggregate_record = json::object{
            {"kind", "transport-state-observable-arm"},
            {"combination", "sum-physical-tiles"},
            {"request_index", index},
            {"observable_identity", observable.identity},
            {"observable_checkpoint_identity",
             observable.checkpoint_identity},
            {"transport_state",
             compact_transport_state_reference(state)},
            {"output_epsilon_contract", observable.epsilon_record},
            {"tail_policy",
             transport_tail_policy_name(observable.tail_policy)},
            {"projection_mode",
             transport_projection_mode_name(observable.tail_policy)},
            {"rows", std::move(row_records)}};
        contraction_inputs.push_back(std::move(input));
      }
      session->pending_line_integrations += observable_count;
      reservation_live = true;
    }

    const auto release_reservation = [&]() {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (!reservation_live) return;
      if (session->pending_line_integrations < observable_count)
        throw std::logic_error(
            "native transport contraction reservation accounting underflow");
      session->pending_line_integrations -= observable_count;
      reservation_live = false;
    };
    struct TransportContractionReservationGuard final {
      decltype(release_reservation)& release;
      bool& live;
      ~TransportContractionReservationGuard() noexcept {
        if (!live) return;
        try {
          release();
        } catch (...) {
          std::terminate();
        }
      }
    } reservation_guard{release_reservation, reservation_live};

    const auto operation_started = std::chrono::steady_clock::now();
    std::unique_ptr<AcbPrecisionLease> acb_lease;
    if (session->domain == "acb") {
      acb_lease = std::make_unique<AcbPrecisionLease>(
          session->precision_bits);
      ComplexBall::set_precision(session->precision_bits);
    }
    auto contracted = contract_transport_arm(
        session->domain, session->precision_bits, state->plan_owner(),
        state->arm_name(), state->tile_sources(), contraction_inputs, state);
    const auto contraction_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - operation_started).count();

    try {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_line_integrations < observable_count)
        throw std::logic_error(
            "native transport contraction reservation accounting underflow");
      session->pending_line_integrations -= observable_count;
      reservation_live = false;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during transport contraction");
      state->require_contraction_counter_capacity(observable_count);
      session->line_results.reserve(
          session->line_results.size() + observable_count);
      std::vector<std::string> inserted;
      inserted.reserve(observable_count);
      try {
        for (const auto& result : contracted) {
          if (!session->line_results.emplace(
                  result.aggregate->handle(), result.aggregate).second)
            throw std::logic_error(
                "transport contraction line handle collided during publication");
          inserted.push_back(result.aggregate->handle());
        }
      } catch (...) {
        for (const auto& handle : inserted)
          session->line_results.erase(handle);
        throw;
      }
      state->note_contraction_success(observable_count);
      ++session->total_transport_contractions;
      session->total_transport_observables += observable_count;
      session->total_transport_contraction_ms += contraction_ms;
      for (const auto& result : contracted) {
        session->total_line_integrations += result.tile_integrations;
        session->total_line_integration_ms += result.tile_integration_ms;
        for (std::size_t tile = 0; tile < result.tile_integrations; ++tile)
          state->plan_owner()->note_integration();
      }
    } catch (...) {
      if (reservation_live) release_reservation();
      throw;
    }

    json::array output_lines;
    output_lines.reserve(contracted.size());
    for (std::size_t index = 0; index < contracted.size(); ++index) {
      const auto& result = contracted[index];
      output_lines.push_back(json::object{
          {"request_index", index},
          {"observable_identity", result.identity},
          {"session", session->handle},
          {"line", result.aggregate->handle()},
          {"checkpoint_identity", result.aggregate->checkpoint_identity()},
          {"provenance_identity", result.aggregate->provenance_identity()},
          {"scope", line_integration_scope_name(
                        result.aggregate->result().scope)},
          {"epsilon", json::object{
               {"min", result.aggregate->result().value.epsilon.min_power},
               {"max",
                result.aggregate->result().value.epsilon.complete_max}}},
          {"error", encode_error_envelope_summary(
                        result.aggregate->result().value.error)},
          {"tiles", result.tile_integrations},
          {"elapsed_ms", result.elapsed_ms}});
    }
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"capability", kRetainedTransportArmContractionCapability},
        {"native_retained", true}, {"json_coefficients", 0},
        {"transport_state", state->handle()},
        {"transport_state_checkpoint_identity",
         state->checkpoint_identity()},
        {"transport_state_provenance_identity",
         state->provenance_identity()},
        {"arm", state->arm_name()},
        {"observables", observable_count},
        {"lines", std::move(output_lines)},
        {"no_rematching", true}, {"atomic_publication", true},
        {"compact_outputs", true},
        {"streaming_tile_contraction", true},
        {"checkpoint_policy", checkpoint_policy},
        {"elapsed_ms", contraction_ms}};
  }

  if (operation == "transport.contract_pair_stream_begin") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "lower", "upper",
         "checkpoint_policy", "observable"},
        "native transport.contract_pair_stream_begin request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native transport-pair streaming requires rational or Acb coefficients");

    struct StateReference {
      std::string handle;
      std::string checkpoint_identity;
      std::string provenance_identity;
    };
    const auto parse_state_reference = [&](const json::value& raw,
                                           const char* label) {
      const auto& reference = as_object(raw, label);
      require_exact_keys(reference,
                         {"transport_state", "checkpoint_identity",
                          "provenance_identity"}, label);
      StateReference parsed{
          required_string(reference, "transport_state"),
          required_string(reference, "checkpoint_identity"),
          required_string(reference, "provenance_identity")};
      if (parsed.handle.empty() || parsed.checkpoint_identity.empty() ||
          parsed.provenance_identity.empty())
        throw std::invalid_argument(
            std::string(label) + " contains an empty exact state token");
      return parsed;
    };
    const std::array<StateReference, 2> state_references{
        parse_state_reference(root.at("lower"),
                              "streamed transport-pair lower state"),
        parse_state_reference(root.at("upper"),
                              "streamed transport-pair upper state")};

    const auto& checkpoint_policy = as_object(
        root.at("checkpoint_policy"),
        "streamed transport-pair checkpoint policy");
    require_exact_keys(checkpoint_policy, {"schema", "root"},
                       "streamed transport-pair checkpoint policy");
    if (required_string(checkpoint_policy, "schema") !=
        "diffexp2-deterministic-transport-pair-contraction-checkpoints-v1")
      throw std::invalid_argument(
          "unsupported streamed transport-pair checkpoint policy schema");
    const auto checkpoint_root = required_string(checkpoint_policy, "root");
    if (checkpoint_root.empty())
      throw std::invalid_argument(
          "streamed transport-pair checkpoint root cannot be empty");

    const auto& observable = as_object(
        root.at("observable"), "streamed transport-pair observable");
    if (observable.if_contains("divergent_cancellation") != nullptr)
      require_exact_keys(
          observable,
          {"identity", "checkpoint_identity", "epsilon", "tail_policy",
           "divergent_cancellation"},
          "streamed transport-pair observable");
    else
      require_exact_keys(
          observable,
          {"identity", "checkpoint_identity", "epsilon", "tail_policy"},
          "streamed transport-pair observable");
    const auto identity = required_string(observable, "identity");
    const auto observable_checkpoint = required_string(
        observable, "checkpoint_identity");
    if (identity.empty() || observable_checkpoint.empty())
      throw std::invalid_argument(
          "streamed transport-pair observable identities cannot be empty");
    const auto epsilon = parse_observable_epsilon_contract(
        observable.at("epsilon"),
        "streamed transport-pair output epsilon contract");
    const auto epsilon_record = as_object(
        observable.at("epsilon"),
        "streamed transport-pair output epsilon contract");
    const auto tail_policy = parse_transport_tail_policy(
        observable.at("tail_policy"),
        "streamed transport-pair tail policy");
    if (tail_policy != TransportTailPolicy::Stored)
      throw std::invalid_argument(
          "transport-pair tile streaming currently supports only stored tails");
    std::optional<BoundedDivergentCancellation> divergent_cancellation;
    if (const auto* policy =
            observable.if_contains("divergent_cancellation")) {
      if (session->domain != "acb")
        throw std::invalid_argument(
            "bounded divergent cancellation is restricted to Acb transport-pair streams");
      divergent_cancellation = parse_bounded_divergent_cancellation(
          *policy, "streamed transport-pair divergent-cancellation policy");
    }

    std::array<std::shared_ptr<StoredTransportArmState>, 2> states;
    std::shared_ptr<TransportPairObservableStream> stream;
    std::string stream_handle;
    std::string stream_checkpoint;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      if (!session->transport_pair_streams.empty())
        throw std::invalid_argument(
            "persistent solver session already owns an active transport-pair stream");
      for (std::size_t side = 0; side < 2; ++side) {
        const auto found = session->transport_states.find(
            state_references[side].handle);
        if (found == session->transport_states.end())
          throw std::invalid_argument(
              "unknown or released transport-arm state for streamed pair contraction");
        states[side] = found->second;
        if (state_references[side].checkpoint_identity !=
                states[side]->checkpoint_identity() ||
            state_references[side].provenance_identity !=
                states[side]->provenance_identity())
          throw std::invalid_argument(
              "streamed transport-pair state binding is stale");
      }
      require_transport_pair_compatibility(
          states[0], states[1], session->domain);
      states[0]->require_contraction_counter_capacity(1);
      states[1]->require_contraction_counter_capacity(1);
      if (session->total_transport_contractions >
              std::numeric_limits<std::uint64_t>::max() - 2 ||
          session->total_transport_observables >
              std::numeric_limits<std::uint64_t>::max() - 2 ||
          session->total_transport_pair_contractions ==
              std::numeric_limits<std::uint64_t>::max() ||
          session->total_transport_pair_observables ==
              std::numeric_limits<std::uint64_t>::max())
        throw std::overflow_error(
            "streamed transport-pair session counter overflow");
      if (session->line_results.size() +
              session->pending_line_integrations >=
          session->line_result_capacity)
        throw std::invalid_argument(
            "persistent line-result capacity is exhausted by streamed transport-pair contraction");
      if (session->next_transport_pair_stream ==
              std::numeric_limits<std::uint64_t>::max() ||
          session->next_line_result ==
              std::numeric_limits<std::uint64_t>::max())
        throw std::overflow_error(
            "streamed transport-pair handle counter overflow");
      stream_handle = "pair-stream:" + session->handle + ":" +
          std::to_string(session->next_transport_pair_stream++);
      stream_checkpoint = checkpoint_root + ":stream";
      const auto line_handle =
          "line:" + std::to_string(session->next_line_result++);
      stream = std::make_shared<TransportPairObservableStream>(
          stream_handle, stream_checkpoint, line_handle, checkpoint_root,
          session->domain, session->precision_bits, states, identity,
          observable_checkpoint, epsilon, epsilon_record, tail_policy,
          divergent_cancellation);
      if (!session->transport_pair_streams.emplace(
              stream_handle, stream).second)
        throw std::logic_error(
            "streamed transport-pair handle collided at publication");
      ++session->pending_line_integrations;
    }
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"capability", kRetainedTransportPairStreamCapability},
        {"stream", stream_handle},
        {"stream_checkpoint_identity", stream_checkpoint},
        {"observable_identity", identity},
        {"observable_checkpoint_identity", observable_checkpoint},
        {"lower", json::object{
             {"transport_state", states[0]->handle()},
             {"checkpoint_identity", states[0]->checkpoint_identity()},
             {"provenance_identity", states[0]->provenance_identity()},
             {"tiles", stream->expected_tiles()[0]}}},
        {"upper", json::object{
             {"transport_state", states[1]->handle()},
             {"checkpoint_identity", states[1]->checkpoint_identity()},
             {"provenance_identity", states[1]->provenance_identity()},
             {"tiles", stream->expected_tiles()[1]}}},
        {"next_side", "lower"}, {"next_tile", 0},
        {"tail_policy", "stored"},
        {"atomic_publication", true},
        {"checkpoint_policy", checkpoint_policy}};
  }

  if (operation == "transport.contract_pair_stream_add_tile") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "stream",
         "stream_checkpoint_identity", "side", "tile", "row"},
        "native transport.contract_pair_stream_add_tile request");
    const auto stream_handle = required_string(root, "stream");
    const auto stream_checkpoint = required_string(
        root, "stream_checkpoint_identity");
    const auto side_name = required_string(root, "side");
    const std::size_t side = side_name == "lower" ? 0 :
        side_name == "upper" ? 1 : 2;
    if (side > 1)
      throw std::invalid_argument(
          "streamed transport-pair tile side must be lower or upper");
    const auto tile = static_cast<std::size_t>(
        as_u32(root.at("tile"), "streamed transport-pair tile index"));
    const auto& row = as_object(
        root.at("row"), "streamed transport-pair prepared row");
    std::shared_ptr<TransportPairObservableStream> stream;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto found = session->transport_pair_streams.find(stream_handle);
      if (found == session->transport_pair_streams.end())
        throw std::invalid_argument(
            "unknown, aborted, or finished transport-pair stream");
      stream = found->second;
      if (stream_checkpoint != stream->stream_checkpoint_identity())
        throw std::invalid_argument(
            "streamed transport-pair checkpoint binding is stale");
    }
    auto progress = stream->add_tile(side, tile, row);
    progress["status"] = "ok";
    progress["session"] = session->handle;
    progress["capability"] = kRetainedTransportPairStreamCapability;
    progress["stream"] = stream_handle;
    progress["stream_checkpoint_identity"] = stream_checkpoint;
    progress["atomic_publication"] = true;
    return progress;
  }

  if (operation == "transport.contract_pair_stream_abort") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "stream",
         "stream_checkpoint_identity"},
        "native transport.contract_pair_stream_abort request");
    const auto stream_handle = required_string(root, "stream");
    const auto stream_checkpoint = required_string(
        root, "stream_checkpoint_identity");
    const auto prefix = "pair-stream:" + session->handle + ":";
    if (stream_handle.rfind(prefix, 0) != 0)
      throw std::invalid_argument(
          "transport-pair stream belongs to a different session");
    std::shared_ptr<TransportPairObservableStream> stream;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->transport_pair_streams.find(stream_handle);
      if (found == session->transport_pair_streams.end())
        return json::object{
            {"status", "ok"}, {"session", session->handle},
            {"capability", kRetainedTransportPairStreamCapability},
            {"stream", stream_handle}, {"aborted", false},
            {"already_absent", true}};
      stream = found->second;
      if (stream_checkpoint != stream->stream_checkpoint_identity())
        throw std::invalid_argument(
            "streamed transport-pair checkpoint binding is stale");
      session->transport_pair_streams.erase(found);
    }
    stream->abort();
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_line_integrations == 0)
        throw std::logic_error(
            "streamed transport-pair abort reservation underflow");
      --session->pending_line_integrations;
    }
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"capability", kRetainedTransportPairStreamCapability},
        {"stream", stream_handle},
        {"stream_checkpoint_identity", stream_checkpoint},
        {"aborted", true}, {"atomic_publication", true}};
  }

  if (operation == "transport.contract_pair_stream_finish") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "stream",
         "stream_checkpoint_identity"},
        "native transport.contract_pair_stream_finish request");
    const auto stream_handle = required_string(root, "stream");
    const auto stream_checkpoint = required_string(
        root, "stream_checkpoint_identity");
    std::shared_ptr<TransportPairObservableStream> stream;
    bool reservation_live = false;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto found = session->transport_pair_streams.find(stream_handle);
      if (found == session->transport_pair_streams.end())
        throw std::invalid_argument(
            "unknown, aborted, or finished transport-pair stream");
      stream = found->second;
      if (stream_checkpoint != stream->stream_checkpoint_identity())
        throw std::invalid_argument(
            "streamed transport-pair checkpoint binding is stale");
      session->transport_pair_streams.erase(found);
      reservation_live = true;
    }
    const auto release_reservation = [&]() {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (!reservation_live) return;
      if (session->pending_line_integrations == 0)
        throw std::logic_error(
            "streamed transport-pair finish reservation underflow");
      --session->pending_line_integrations;
      reservation_live = false;
    };
    struct StreamFinishReservationGuard final {
      decltype(release_reservation)& release;
      bool& live;
      ~StreamFinishReservationGuard() noexcept {
        if (!live) return;
        try {
          release();
        } catch (...) {
          std::terminate();
        }
      }
    } reservation_guard{release_reservation, reservation_live};

    std::unique_ptr<AcbPrecisionLease> acb_lease;
    if (session->domain == "acb") {
      acb_lease = std::make_unique<AcbPrecisionLease>(
          session->precision_bits);
      ComplexBall::set_precision(session->precision_bits);
    }
    auto finished = stream->finish();
    const auto& states = stream->states();
    const auto tile_integrations = static_cast<std::uint64_t>(
        finished.tiles[0]) + static_cast<std::uint64_t>(finished.tiles[1]);
    try {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_line_integrations == 0)
        throw std::logic_error(
            "streamed transport-pair finish reservation underflow");
      --session->pending_line_integrations;
      reservation_live = false;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during streamed transport-pair contraction");
      states[0]->require_contraction_counter_capacity(1);
      states[1]->require_contraction_counter_capacity(1);
      if (session->total_transport_contractions >
              std::numeric_limits<std::uint64_t>::max() - 2 ||
          session->total_transport_observables >
              std::numeric_limits<std::uint64_t>::max() - 2 ||
          session->total_transport_pair_contractions ==
              std::numeric_limits<std::uint64_t>::max() ||
          session->total_transport_pair_observables ==
              std::numeric_limits<std::uint64_t>::max() ||
          tile_integrations >
              std::numeric_limits<std::uint64_t>::max() -
                  session->total_line_integrations)
        throw std::overflow_error(
            "streamed transport-pair session counter overflow at publication");
      if (!session->line_results.emplace(
              finished.line->handle(), finished.line).second)
        throw std::logic_error(
            "streamed transport-pair line handle collided during publication");
      states[0]->note_contraction_success(1);
      states[1]->note_contraction_success(1);
      session->total_transport_contractions += 2;
      session->total_transport_observables += 2;
      ++session->total_transport_pair_contractions;
      ++session->total_transport_pair_observables;
      session->total_transport_contraction_ms +=
          finished.arm_integration_ms[0] +
          finished.arm_integration_ms[1];
      session->total_transport_pair_contraction_ms += finished.elapsed_ms;
      session->total_line_integrations += tile_integrations;
      session->total_line_integration_ms +=
          finished.arm_integration_ms[0] +
          finished.arm_integration_ms[1];
      for (std::size_t side = 0; side < 2; ++side)
        for (std::size_t tile = 0; tile < finished.tiles[side]; ++tile)
          states[side]->plan_owner()->note_integration();
    } catch (...) {
      if (reservation_live) release_reservation();
      throw;
    }

    const auto& line = finished.line;
    json::array output_lines;
    output_lines.push_back(json::object{
        {"request_index", 0},
        {"observable_identity", stream->identity()},
        {"session", session->handle}, {"line", line->handle()},
        {"checkpoint_identity", line->checkpoint_identity()},
        {"provenance_identity", line->provenance_identity()},
        {"combination", "negative-lower-plus-upper"},
        {"scope", line_integration_scope_name(line->result().scope)},
        {"epsilon", json::object{
             {"min", line->result().value.epsilon.min_power},
             {"max", line->result().value.epsilon.complete_max}}},
        {"error", encode_error_envelope_summary(
                      line->result().value.error)},
        {"lower_tiles", finished.tiles[0]},
        {"upper_tiles", finished.tiles[1]},
        {"elapsed_ms", finished.elapsed_ms}});
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"capability", kRetainedTransportPairStreamCapability},
        {"native_retained", true}, {"json_coefficients", 0},
        {"stream", stream_handle},
        {"stream_checkpoint_identity", stream_checkpoint},
        {"lower", json::object{
             {"transport_state", states[0]->handle()},
             {"checkpoint_identity", states[0]->checkpoint_identity()},
             {"provenance_identity", states[0]->provenance_identity()},
             {"tiles", finished.tiles[0]},
             {"elapsed_ms", finished.arm_integration_ms[0]}}},
        {"upper", json::object{
             {"transport_state", states[1]->handle()},
             {"checkpoint_identity", states[1]->checkpoint_identity()},
             {"provenance_identity", states[1]->provenance_identity()},
             {"tiles", finished.tiles[1]},
             {"elapsed_ms", finished.arm_integration_ms[1]}}},
        {"combination", "negative-lower-plus-upper"},
        {"observables", 1}, {"lines", std::move(output_lines)},
        {"tile_integrations", tile_integrations},
        {"no_remarching", true}, {"no_rematching", true},
        {"concurrent_arms", false}, {"max_parallel_arms", 1},
        {"streaming_tile_contraction", true},
        {"persistent_tile_stream", true},
        {"atomic_publication", true}, {"compact_outputs", true},
        {"checkpoint_policy", json::object{
             {"schema",
              "diffexp2-deterministic-transport-pair-contraction-checkpoints-v1"},
             {"root", stream->checkpoint_root()}}},
        {"elapsed_ms", finished.elapsed_ms}};
  }

  if (operation == "transport.contract_pair") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "lower", "upper",
         "checkpoint_policy", "observables"},
        "native transport.contract_pair request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native transport-pair contraction requires rational or Acb coefficients");

    struct StateReference {
      std::string handle;
      std::string checkpoint_identity;
      std::string provenance_identity;
    };
    const auto parse_state_reference = [&](const json::value& raw,
                                           const char* label) {
      const auto& reference = as_object(raw, label);
      require_exact_keys(reference,
                         {"transport_state", "checkpoint_identity",
                          "provenance_identity"}, label);
      StateReference parsed{
          required_string(reference, "transport_state"),
          required_string(reference, "checkpoint_identity"),
          required_string(reference, "provenance_identity")};
      if (parsed.handle.empty() || parsed.checkpoint_identity.empty() ||
          parsed.provenance_identity.empty())
        throw std::invalid_argument(
            std::string(label) + " contains an empty exact state token");
      return parsed;
    };
    const std::array<StateReference, 2> state_references{
        parse_state_reference(root.at("lower"),
                              "native transport-pair lower state"),
        parse_state_reference(root.at("upper"),
                              "native transport-pair upper state")};

    const auto& checkpoint_policy = as_object(
        root.at("checkpoint_policy"),
        "native transport-pair contraction checkpoint policy");
    require_exact_keys(checkpoint_policy, {"schema", "root"},
                       "native transport-pair contraction checkpoint policy");
    if (required_string(checkpoint_policy, "schema") !=
        "diffexp2-deterministic-transport-pair-contraction-checkpoints-v1")
      throw std::invalid_argument(
          "unsupported native transport-pair contraction checkpoint policy schema");
    const auto checkpoint_root = required_string(checkpoint_policy, "root");
    if (checkpoint_root.empty())
      throw std::invalid_argument(
          "native transport-pair contraction checkpoint root cannot be empty");

    struct PendingPairObservable {
      std::string identity;
      std::string checkpoint_identity;
      // The request root outlives this synchronous operation.  Borrow its
      // immutable prepared rows instead of deep-copying every exact kernel
      // through pending and per-arm records.
      std::array<std::vector<const json::object*>, 2> rows;
      ObservableEpsilonContract epsilon;
      json::object epsilon_record;
      TransportTailPolicy tail_policy = TransportTailPolicy::Stored;
      std::optional<BoundedDivergentCancellation> divergent_cancellation;
    };
    std::vector<PendingPairObservable> pending_observables;
    const auto& raw_observables = as_array(
        root.at("observables"), "native transport-pair observables");
    pending_observables.reserve(raw_observables.size());
    std::set<std::string> observable_identities;
    std::set<std::string> observable_checkpoints;
    for (const auto& raw_observable : raw_observables) {
      const auto& observable = as_object(
          raw_observable, "native transport-pair observable");
      if (observable.if_contains("divergent_cancellation") != nullptr)
        require_exact_keys(
            observable,
            {"identity", "checkpoint_identity", "lower_integrand_rows",
             "upper_integrand_rows", "epsilon", "tail_policy",
             "divergent_cancellation"},
            "native transport-pair observable");
      else
        require_exact_keys(
            observable,
            {"identity", "checkpoint_identity", "lower_integrand_rows",
             "upper_integrand_rows", "epsilon", "tail_policy"},
            "native transport-pair observable");
      PendingPairObservable parsed;
      parsed.identity = required_string(observable, "identity");
      parsed.checkpoint_identity = required_string(
          observable, "checkpoint_identity");
      if (parsed.identity.empty() || parsed.checkpoint_identity.empty() ||
          !observable_identities.insert(parsed.identity).second ||
          !observable_checkpoints.insert(parsed.checkpoint_identity).second)
        throw std::invalid_argument(
            "native transport-pair observable identities and checkpoints must be nonempty and pairwise unique");
      for (std::size_t side = 0; side < 2; ++side) {
        const auto* key = side == 0 ? "lower_integrand_rows"
                                    : "upper_integrand_rows";
        for (const auto& raw_row : as_array(observable.at(key), key))
          parsed.rows[side].push_back(&as_object(raw_row, key));
      }
      parsed.epsilon = parse_observable_epsilon_contract(
          observable.at("epsilon"),
          "native transport-pair output epsilon contract");
      parsed.epsilon_record = as_object(
          observable.at("epsilon"),
          "native transport-pair output epsilon contract");
      parsed.tail_policy = parse_transport_tail_policy(
          observable.at("tail_policy"),
          "native transport-pair tail policy");
      if (const auto* policy =
              observable.if_contains("divergent_cancellation")) {
        if (session->domain != "acb")
          throw std::invalid_argument(
              "bounded divergent cancellation is restricted to Acb transport-pair observables");
        parsed.divergent_cancellation =
            parse_bounded_divergent_cancellation(
                *policy,
                "native transport-pair divergent-cancellation policy");
      }
      pending_observables.push_back(std::move(parsed));
    }

    const auto observable_count = pending_observables.size();
    std::array<std::shared_ptr<StoredTransportArmState>, 2> states;
    std::array<std::vector<TransportObservableContractionInput>, 2>
        arm_inputs;
    std::vector<json::object> pair_aggregate_records;
    pair_aggregate_records.reserve(observable_count);
    bool reservation_live = false;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      for (std::size_t side = 0; side < 2; ++side) {
        const auto found = session->transport_states.find(
            state_references[side].handle);
        if (found == session->transport_states.end())
          throw std::invalid_argument(
              "unknown or released native transport-arm state for paired contraction");
        states[side] = found->second;
        if (state_references[side].checkpoint_identity !=
                states[side]->checkpoint_identity() ||
            state_references[side].provenance_identity !=
                states[side]->provenance_identity())
          throw std::invalid_argument(
              "native transport-pair state binding is stale");
      }
      require_transport_pair_compatibility(
          states[0], states[1], session->domain);
      states[0]->require_contraction_counter_capacity(observable_count);
      states[1]->require_contraction_counter_capacity(observable_count);
      if (session->total_transport_contractions >
              std::numeric_limits<std::uint64_t>::max() - 2 ||
          observable_count >
              std::numeric_limits<std::uint64_t>::max() / 2 ||
          2 * static_cast<std::uint64_t>(observable_count) >
              std::numeric_limits<std::uint64_t>::max() -
                  session->total_transport_observables ||
          session->total_transport_pair_contractions ==
              std::numeric_limits<std::uint64_t>::max() ||
          observable_count >
              std::numeric_limits<std::uint64_t>::max() -
                  session->total_transport_pair_observables)
        throw std::overflow_error(
            "native transport-pair session counter overflow");
      for (const auto& observable : pending_observables) {
        for (std::size_t side = 0; side < 2; ++side) {
          if (observable.rows[side].size() !=
              states[side]->tile_sources().size())
            throw std::invalid_argument(
                side == 0
                    ? "native transport-pair observable must provide one lower row per lower retained tile"
                    : "native transport-pair observable must provide one upper row per upper retained tile");
          if (observable.epsilon.required_complete_max >
              states[side]->public_required_complete_max())
            throw std::invalid_argument(
                "native transport-pair observable required epsilon maximum exceeds an arm-state public target");
        }
      }
      if (observable_count > session->line_result_capacity -
                                 std::min(
                                     session->line_result_capacity,
                                     session->line_results.size() +
                                         session->pending_line_integrations))
        throw std::invalid_argument(
            "persistent line-result capacity is exhausted by transport-pair contraction");

      arm_inputs[0].reserve(observable_count);
      arm_inputs[1].reserve(observable_count);
      for (std::size_t index = 0; index < observable_count; ++index) {
        const auto& observable = pending_observables[index];
        const auto observable_root = checkpoint_root + ":observable:" +
                                     std::to_string(index + 1);
        std::array<json::array, 2> pair_row_records;
        for (std::size_t side = 0; side < 2; ++side) {
          const auto arm_name = side == 0 ? "lower" : "upper";
          TransportObservableContractionInput input;
          input.identity = observable.identity + ":" + arm_name;
          input.checkpoint_identity =
              observable_root + ":" + arm_name + ":scratch";
          input.checkpoint_root = observable_root + ":" + arm_name;
          input.borrowed_rows = observable.rows[side];
          input.epsilon = observable.epsilon;
          input.epsilon_record = observable.epsilon_record;
          input.tail_policy = observable.tail_policy;
          input.divergent_cancellation =
              observable.divergent_cancellation;
          input.aggregate_handle =
              "private:" + observable_root + ":" + arm_name + ":aggregate";
          input.projected_local_handles.reserve(input.row_count());
          pair_row_records[side].reserve(input.row_count());
          for (std::size_t tile = 0; tile < input.row_count(); ++tile) {
            input.projected_local_handles.push_back(
                "private:" + arm_checkpoint_identity(
                                 input.checkpoint_root, arm_name,
                                 "projected", tile + 1));
            const auto row_identity = required_string(
                input.row(tile), "exact_identity");
            pair_row_records[side].push_back(json::object{
                {"tile", tile}, {"exact_identity", row_identity},
                {"entries",
                 compact_prepared_row_entry_facts(input.row(tile))}});
          }
          input.aggregate_record = json::object{
              {"kind", "transport-state-observable-arm"},
              {"combination", "sum-physical-tiles"},
              {"request_index", index},
              {"observable_identity", input.identity},
              {"observable_checkpoint_identity",
               input.checkpoint_identity},
              {"transport_state",
               compact_transport_state_reference(states[side])},
              {"output_epsilon_contract", observable.epsilon_record},
            {"tail_policy",
               transport_tail_policy_name(observable.tail_policy)},
            {"projection_mode",
             transport_projection_mode_name(observable.tail_policy)},
            {"divergent_cancellation",
             observable.divergent_cancellation.has_value()
                 ? json::value(encode_bounded_divergent_cancellation(
                       *observable.divergent_cancellation))
                 : json::value(
                       json::object{{"mode", "exact-singleton"}})},
              {"rows", pair_row_records[side]},
              {"tile_count", input.row_count()}};
          arm_inputs[side].push_back(std::move(input));
        }
        const auto line_handle =
            "line:" + std::to_string(session->next_line_result++);
        pair_aggregate_records.push_back(json::object{
            {"kind", "transport-state-observable-pair"},
            {"combination", "negative-lower-plus-upper"},
            {"no_remarching", true}, {"no_rematching", true},
            {"concurrent_arms", true},
            {"request_index", index},
            {"observable_identity", observable.identity},
            {"observable_checkpoint_identity",
             observable.checkpoint_identity},
            {"output_epsilon_contract", observable.epsilon_record},
            {"tail_policy",
             transport_tail_policy_name(observable.tail_policy)},
            {"projection_mode",
             transport_projection_mode_name(observable.tail_policy)},
            {"divergent_cancellation",
             observable.divergent_cancellation.has_value()
                 ? json::value(encode_bounded_divergent_cancellation(
                       *observable.divergent_cancellation))
                 : json::value(
                       json::object{{"mode", "exact-singleton"}})},
            {"lower", json::object{
                 {"transport_state",
                  compact_transport_state_reference(states[0])},
                 {"rows", pair_row_records[0]},
                 {"tile_count", observable.rows[0].size()}}},
            {"upper", json::object{
                 {"transport_state",
                  compact_transport_state_reference(states[1])},
                 {"rows", pair_row_records[1]},
                 {"tile_count", observable.rows[1].size()}}},
            {"line_handle", line_handle}});
      }
      session->pending_line_integrations += observable_count;
      reservation_live = true;
    }

    const auto release_reservation = [&]() {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (!reservation_live) return;
      if (session->pending_line_integrations < observable_count)
        throw std::logic_error(
            "native transport-pair reservation accounting underflow");
      session->pending_line_integrations -= observable_count;
      reservation_live = false;
    };
    struct TransportPairReservationGuard final {
      decltype(release_reservation)& release;
      bool& live;
      ~TransportPairReservationGuard() noexcept {
        if (!live) return;
        try {
          release();
        } catch (...) {
          std::terminate();
        }
      }
    } reservation_guard{release_reservation, reservation_live};

    const auto operation_started = std::chrono::steady_clock::now();
    std::unique_ptr<AcbPrecisionLease> acb_lease;
    if (session->domain == "acb") {
      acb_lease = std::make_unique<AcbPrecisionLease>(
          session->precision_bits);
      ComplexBall::set_precision(session->precision_bits);
    }
    std::array<std::vector<TransportObservableContractionResult>, 2>
        arm_results;
    std::array<double, 2> arm_elapsed_ms{0.0, 0.0};
    if (observable_count != 0) {
      // Pair contraction is deliberately sequential.  Each arm now streams
      // one projected tile at a time; running both arms together would still
      // double the largest projection/integration scratch allocation.
      try {
        for (std::size_t side = 0; side < 2; ++side) {
          if (session->domain == "acb")
            ComplexBall::set_precision(session->precision_bits);
          const auto started = std::chrono::steady_clock::now();
          arm_results[side] = contract_transport_arm(
              session->domain, session->precision_bits,
              states[side]->plan_owner(),
              side == 0 ? "lower" : "upper",
              states[side]->tile_sources(), arm_inputs[side], states[side]);
          arm_elapsed_ms[side] =
              std::chrono::duration<double, std::milli>(
                  std::chrono::steady_clock::now() - started)
                  .count();
        }
      } catch (...) {
        release_reservation();
        throw;
      }
    }

    std::vector<std::shared_ptr<StoredLineResult>> combined_lines;
    combined_lines.reserve(observable_count);
    for (std::size_t index = 0; index < observable_count; ++index) {
      if (arm_results[0].size() != observable_count ||
          arm_results[1].size() != observable_count)
        throw std::logic_error(
            "native transport-pair contraction returned the wrong arm result count");
      auto aggregate_record = std::move(pair_aggregate_records[index]);
      const auto line_handle = required_string(
          aggregate_record, "line_handle");
      aggregate_record.erase("line_handle");
      auto combined = build_compact_transport_pair_observable_line(
          line_handle, pending_observables[index].checkpoint_identity,
          std::move(aggregate_record), states[0], states[1],
          arm_results[0][index].aggregate,
          arm_results[1][index].aggregate);
      if (combined->result().value.epsilon.complete_max <
          pending_observables[index].epsilon.required_complete_max)
        throw std::domain_error(
            "native transport-pair aggregate does not cover its required epsilon maximum");
      if (pending_observables[index].tail_policy ==
              TransportTailPolicy::Require &&
          combined->result().scope !=
              LineIntegrationScope::FullLocalWithCertifiedTail)
        throw std::domain_error(
            "native transport-pair contraction requires certified full-local tails on both arms");
      combined_lines.push_back(std::move(combined));
    }
    const auto pair_wall_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - operation_started).count();

    std::uint64_t tile_integrations = 0;
    double tile_integration_ms = 0.0;
    for (std::size_t side = 0; side < 2; ++side)
      for (const auto& result : arm_results[side]) {
        tile_integrations += result.tile_integrations;
        tile_integration_ms += result.tile_integration_ms;
      }
    try {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_line_integrations < observable_count)
        throw std::logic_error(
            "native transport-pair reservation accounting underflow");
      session->pending_line_integrations -= observable_count;
      reservation_live = false;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during transport-pair contraction");
      states[0]->require_contraction_counter_capacity(observable_count);
      states[1]->require_contraction_counter_capacity(observable_count);
      if (session->total_transport_contractions >
              std::numeric_limits<std::uint64_t>::max() - 2 ||
          observable_count >
              std::numeric_limits<std::uint64_t>::max() / 2 ||
          2 * static_cast<std::uint64_t>(observable_count) >
              std::numeric_limits<std::uint64_t>::max() -
                  session->total_transport_observables ||
          session->total_transport_pair_contractions ==
              std::numeric_limits<std::uint64_t>::max() ||
          observable_count >
              std::numeric_limits<std::uint64_t>::max() -
                  session->total_transport_pair_observables)
        throw std::overflow_error(
            "native transport-pair session counter overflow at publication");
      session->line_results.reserve(
          session->line_results.size() + observable_count);
      std::vector<std::string> inserted;
      inserted.reserve(observable_count);
      try {
        for (const auto& line : combined_lines) {
          if (!session->line_results.emplace(line->handle(), line).second)
            throw std::logic_error(
                "transport-pair line handle collided during publication");
          inserted.push_back(line->handle());
        }
      } catch (...) {
        for (const auto& handle : inserted)
          session->line_results.erase(handle);
        throw;
      }
      states[0]->note_contraction_success(observable_count);
      states[1]->note_contraction_success(observable_count);
      session->total_transport_contractions += 2;
      session->total_transport_observables +=
          2 * static_cast<std::uint64_t>(observable_count);
      ++session->total_transport_pair_contractions;
      session->total_transport_pair_observables += observable_count;
      session->total_transport_contraction_ms +=
          arm_elapsed_ms[0] + arm_elapsed_ms[1];
      session->total_transport_pair_contraction_ms += pair_wall_ms;
      session->total_line_integrations += tile_integrations;
      session->total_line_integration_ms += tile_integration_ms;
      for (std::size_t side = 0; side < 2; ++side)
        for (const auto& result : arm_results[side])
          for (std::size_t tile = 0; tile < result.tile_integrations; ++tile)
            states[side]->plan_owner()->note_integration();
    } catch (...) {
      if (reservation_live) release_reservation();
      throw;
    }

    json::array output_lines;
    output_lines.reserve(combined_lines.size());
    for (std::size_t index = 0; index < combined_lines.size(); ++index) {
      const auto& line = combined_lines[index];
      output_lines.push_back(json::object{
          {"request_index", index},
          {"observable_identity", pending_observables[index].identity},
          {"session", session->handle},
          {"line", line->handle()},
          {"checkpoint_identity", line->checkpoint_identity()},
          {"provenance_identity", line->provenance_identity()},
          {"combination", "negative-lower-plus-upper"},
          {"scope", line_integration_scope_name(line->result().scope)},
          {"epsilon", json::object{
               {"min", line->result().value.epsilon.min_power},
               {"max", line->result().value.epsilon.complete_max}}},
          {"error", encode_error_envelope_summary(
                        line->result().value.error)},
          {"lower_tiles", arm_results[0][index].tile_integrations},
          {"upper_tiles", arm_results[1][index].tile_integrations},
          {"elapsed_ms", arm_results[0][index].elapsed_ms +
                             arm_results[1][index].elapsed_ms +
                             line->elapsed_ms()}});
    }
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"capability", kRetainedTransportPairContractionCapability},
        {"native_retained", true}, {"json_coefficients", 0},
        {"lower", json::object{
             {"transport_state", states[0]->handle()},
             {"checkpoint_identity", states[0]->checkpoint_identity()},
             {"provenance_identity", states[0]->provenance_identity()},
             {"tiles", states[0]->tile_sources().size()},
             {"elapsed_ms", arm_elapsed_ms[0]}}},
        {"upper", json::object{
             {"transport_state", states[1]->handle()},
             {"checkpoint_identity", states[1]->checkpoint_identity()},
             {"provenance_identity", states[1]->provenance_identity()},
             {"tiles", states[1]->tile_sources().size()},
             {"elapsed_ms", arm_elapsed_ms[1]}}},
        {"combination", "negative-lower-plus-upper"},
        {"observables", observable_count},
        {"lines", std::move(output_lines)},
        {"tile_integrations", tile_integrations},
        {"no_remarching", true}, {"no_rematching", true},
        {"concurrent_arms", false},
        {"max_parallel_arms", observable_count == 0 ? 0 : 1},
        {"streaming_tile_contraction", true},
        {"atomic_publication", true}, {"compact_outputs", true},
        {"checkpoint_policy", checkpoint_policy},
        {"elapsed_ms", pair_wall_ms}};
  }

  if (operation == "transport.endpoint_batch") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "transport_state",
         "transport_state_checkpoint_identity",
         "transport_state_provenance_identity", "checkpoint_policy",
         "observables"},
        "native transport.endpoint_batch request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native transport endpoint batch requires rational or Acb coefficients");
    const auto state_handle = required_string(root, "transport_state");
    const auto& checkpoint_policy = as_object(
        root.at("checkpoint_policy"),
        "native transport endpoint-batch checkpoint policy");
    require_exact_keys(checkpoint_policy, {"schema", "root"},
                       "native transport endpoint-batch checkpoint policy");
    if (required_string(checkpoint_policy, "schema") !=
        "diffexp2-deterministic-transport-endpoint-checkpoints-v1")
      throw std::invalid_argument(
          "unsupported native transport endpoint-batch checkpoint policy schema");
    const auto checkpoint_root = required_string(checkpoint_policy, "root");
    if (checkpoint_root.empty())
      throw std::invalid_argument(
          "native transport endpoint-batch checkpoint root cannot be empty");

    struct PendingEndpointObservable {
      std::string identity;
      std::string checkpoint_identity;
      json::object row;
      ObservableEpsilonContract epsilon;
      json::object epsilon_record;
      std::string endpoint_handle;
      std::string projected_handle;
      std::string projected_checkpoint_identity;
    };
    std::vector<PendingEndpointObservable> pending_observables;
    const auto& raw_observables = as_array(
        root.at("observables"),
        "native transport endpoint observables");
    pending_observables.reserve(raw_observables.size());
    std::set<std::string> identities;
    std::set<std::string> checkpoints;
    for (const auto& raw_observable : raw_observables) {
      const auto& observable = as_object(
          raw_observable, "native transport endpoint observable");
      require_exact_keys(
          observable,
          {"identity", "checkpoint_identity", "integrand_row", "epsilon"},
          "native transport endpoint observable");
      PendingEndpointObservable parsed;
      parsed.identity = required_string(observable, "identity");
      parsed.checkpoint_identity = required_string(
          observable, "checkpoint_identity");
      if (parsed.identity.empty() || parsed.checkpoint_identity.empty() ||
          !identities.insert(parsed.identity).second ||
          !checkpoints.insert(parsed.checkpoint_identity).second)
        throw std::invalid_argument(
            "native transport endpoint observable identities and checkpoints must be nonempty and pairwise unique");
      parsed.row = as_object(
          observable.at("integrand_row"),
          "native transport endpoint prepared row");
      parsed.epsilon = parse_observable_epsilon_contract(
          observable.at("epsilon"),
          "native transport endpoint epsilon contract");
      parsed.epsilon_record = as_object(
          observable.at("epsilon"),
          "native transport endpoint epsilon contract");
      pending_observables.push_back(std::move(parsed));
    }

    const auto observable_count = pending_observables.size();
    std::shared_ptr<StoredTransportArmState> state;
    ResolvedTransportEndpointBinding binding;
    bool reservation_live = false;
    const auto require_session_counter_capacity = [&]() {
      if (session->total_transport_endpoint_batches ==
              std::numeric_limits<std::uint64_t>::max() ||
          observable_count >
              std::numeric_limits<std::uint64_t>::max() -
                  session->total_transport_endpoint_rows ||
          observable_count >
              std::numeric_limits<std::uint64_t>::max() -
                  session->total_endpoint_limits)
        throw std::overflow_error(
            "native transport endpoint-batch session counter overflow");
    };
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto found = session->transport_states.find(state_handle);
      if (found == session->transport_states.end())
        throw std::invalid_argument(
            "unknown or released native transport-arm state for endpoint batch");
      state = found->second;
      if (required_string(root, "transport_state_checkpoint_identity") !=
              state->checkpoint_identity() ||
          required_string(root, "transport_state_provenance_identity") !=
              state->provenance_identity())
        throw std::invalid_argument(
            "native transport endpoint-batch state binding is stale");
      if (std::string(state->final_local()->scalar_domain()) !=
          session->domain)
        throw std::invalid_argument(
            "native transport endpoint-batch state domain differs from its session");
      binding = resolve_transport_endpoint_binding(state);
      state->require_endpoint_batch_counter_capacity(observable_count);
      require_session_counter_capacity();
      const auto source_dimension = as_u32(
          state->final_local()->summary().at("dimension"),
          "transport endpoint source dimension");
      for (const auto& observable : pending_observables) {
        validate_prepared_rational_row_structure(
            observable.row, source_dimension,
            "native transport endpoint prepared row");
        if (observable.epsilon.required_complete_max >
            state->public_required_complete_max())
          throw std::invalid_argument(
              "native transport endpoint required epsilon maximum exceeds the state public target");
      }
      if (observable_count > session->endpoint_capacity -
                                 std::min(
                                     session->endpoint_capacity,
                                     session->endpoints.size() +
                                         session->pending_endpoint_limits))
        throw std::invalid_argument(
            "persistent endpoint result capacity is exhausted by transport endpoint batch");
      for (std::size_t index = 0; index < observable_count; ++index) {
        auto& observable = pending_observables[index];
        observable.endpoint_handle =
            "e:" + std::to_string(session->next_endpoint++);
        const auto scratch_root = checkpoint_root + ":observable:" +
                                  std::to_string(index + 1);
        observable.projected_handle =
            "private:" + scratch_root + ":projected";
        observable.projected_checkpoint_identity =
            scratch_root + ":projected";
      }
      session->pending_endpoint_limits += observable_count;
      reservation_live = true;
    }

    const auto release_reservation = [&]() {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (!reservation_live) return;
      if (session->pending_endpoint_limits < observable_count)
        throw std::logic_error(
            "native transport endpoint-batch reservation accounting underflow");
      session->pending_endpoint_limits -= observable_count;
      reservation_live = false;
    };
    struct TransportEndpointReservationGuard final {
      decltype(release_reservation)& release;
      bool& live;
      ~TransportEndpointReservationGuard() noexcept {
        if (!live) return;
        try {
          release();
        } catch (...) {
          std::terminate();
        }
      }
    } reservation_guard{release_reservation, reservation_live};

    const auto operation_started = std::chrono::steady_clock::now();
    std::unique_ptr<AcbPrecisionLease> acb_lease;
    if (session->domain == "acb") {
      acb_lease = std::make_unique<AcbPrecisionLease>(
          session->precision_bits);
      ComplexBall::set_precision(session->precision_bits);
    }
    std::vector<std::shared_ptr<StoredEndpointResult>> endpoints;
    endpoints.reserve(observable_count);
    for (const auto& observable : pending_observables)
      endpoints.push_back(build_transport_endpoint_row(
          observable.endpoint_handle, observable.checkpoint_identity,
          observable.identity, observable.row, observable.epsilon,
          observable.epsilon_record, observable.projected_handle,
          observable.projected_checkpoint_identity, session->domain,
          session->precision_bits, state, binding));
    const auto operation_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - operation_started).count();
    double endpoint_ms = 0.0;
    for (const auto& endpoint : endpoints)
      endpoint_ms += endpoint->elapsed_ms();

    try {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_endpoint_limits < observable_count)
        throw std::logic_error(
            "native transport endpoint-batch reservation accounting underflow");
      session->pending_endpoint_limits -= observable_count;
      reservation_live = false;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during transport endpoint batch");
      state->require_endpoint_batch_counter_capacity(observable_count);
      require_session_counter_capacity();
      if (observable_count > session->endpoint_capacity -
                                 std::min(
                                     session->endpoint_capacity,
                                     session->endpoints.size() +
                                         session->pending_endpoint_limits))
        throw std::invalid_argument(
            "persistent endpoint result capacity changed during transport endpoint batch");
      session->endpoints.reserve(
          session->endpoints.size() + observable_count);
      std::vector<std::string> inserted;
      inserted.reserve(observable_count);
      try {
        for (const auto& endpoint : endpoints) {
          if (!session->endpoints.emplace(
                  endpoint->handle(), endpoint).second)
            throw std::logic_error(
                "transport endpoint handle collided during publication");
          inserted.push_back(endpoint->handle());
        }
      } catch (...) {
        for (const auto& handle : inserted)
          session->endpoints.erase(handle);
        throw;
      }
      state->note_endpoint_batch_success(observable_count);
      ++session->total_transport_endpoint_batches;
      session->total_transport_endpoint_rows += observable_count;
      session->total_transport_endpoint_batch_ms += operation_ms;
      session->total_endpoint_limits += observable_count;
      session->total_endpoint_limit_ms += endpoint_ms;
    } catch (...) {
      if (reservation_live) release_reservation();
      throw;
    }

    json::array output_endpoints;
    output_endpoints.reserve(endpoints.size());
    for (std::size_t index = 0; index < endpoints.size(); ++index) {
      const auto& endpoint = endpoints[index];
      auto summary = endpoint->summary();
      output_endpoints.push_back(json::object{
          {"request_index", index},
          {"observable_identity", pending_observables[index].identity},
          {"session", session->handle},
          {"endpoint", endpoint->handle()},
          {"checkpoint_identity", endpoint->checkpoint_identity()},
          {"provenance_identity", endpoint->provenance_identity()},
          {"epsilon", json::object{
               {"min", summary.at("epsilon_min")},
               {"max", summary.at("epsilon_max")}}},
          {"centered", binding.centered},
          {"elapsed_ms", endpoint->elapsed_ms()}});
    }
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"capability", kRetainedTransportEndpointBatchCapability},
        {"native_retained", true}, {"json_coefficients", 0},
        {"transport_state", state->handle()},
        {"transport_state_checkpoint_identity",
         state->checkpoint_identity()},
        {"transport_state_provenance_identity",
         state->provenance_identity()},
        {"arm", binding.arm},
        {"endpoint_exact", binding.source.at("endpoint_exact")},
        {"local_endpoint_exact", binding.local_end.str()},
        {"centered", binding.centered},
        {"approach_direction", binding.approach_direction},
        {"derived_rim", binding.rim.has_value()
             ? json::value(*binding.rim) : json::value(nullptr)},
        {"observables", observable_count},
        {"endpoints", std::move(output_endpoints)},
        {"no_projected_local_publication", true},
        {"atomic_publication", true}, {"compact_outputs", true},
        {"checkpoint_policy", checkpoint_policy},
        {"elapsed_ms", operation_ms}};
  }

