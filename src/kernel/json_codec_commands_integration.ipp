  if (operation == "integration.run_arms") {
    if (root.if_contains("certify_tail") != nullptr)
      require_exact_keys(
          root,
          {"schema", "op", "session", "tile_plan", "anchor",
           "tile_plan_checkpoint_identity", "anchor_checkpoint_identity",
           "epsilon", "refinement", "checkpoint_policy", "lower", "upper",
           "certify_tail"},
          "native integration.run_arms request");
    else
      require_exact_keys(
          root,
          {"schema", "op", "session", "tile_plan", "anchor",
           "tile_plan_checkpoint_identity", "anchor_checkpoint_identity",
           "epsilon", "refinement", "checkpoint_policy", "lower", "upper"},
          "native integration.run_arms request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native whole-arm marching requires rational or Acb coefficients");
    const bool certify_tail = root.if_contains("certify_tail") != nullptr
        ? root.at("certify_tail").as_bool()
        : false;

    const auto& raw_epsilon = as_object(
        root.at("epsilon"), "native whole-arm epsilon contract");
    const auto epsilon_contract = parse_whole_arm_epsilon_contract(
        raw_epsilon, "native whole-arm epsilon contract");
    const auto work_epsilon = epsilon_contract.work;
    const auto required_complete_max =
        epsilon_contract.public_required_complete_max;
    const auto match_required_complete_max =
        epsilon_contract.match_required_complete_max;

    const auto& refinement = as_object(
        root.at("refinement"), "native whole-arm refinement policy");
    require_exact_keys(refinement, {"relative_tolerance", "max_steps"},
                       "native whole-arm refinement policy");
    (void)required_string(refinement, "relative_tolerance");
    if (as_u32(refinement.at("max_steps"),
               "native whole-arm refinement steps") > 32)
      throw std::invalid_argument(
          "native whole-arm refinement steps must lie in 0..32");

    const auto& checkpoint_policy = as_object(
        root.at("checkpoint_policy"),
        "native whole-arm checkpoint policy");
    require_exact_keys(checkpoint_policy, {"schema", "root"},
                       "native whole-arm checkpoint policy");
    if (required_string(checkpoint_policy, "schema") !=
        "diffexp3-deterministic-arm-checkpoints-v1")
      throw std::invalid_argument(
          "unsupported native whole-arm checkpoint policy schema");
    const auto checkpoint_root = required_string(checkpoint_policy, "root");
    if (checkpoint_root.empty())
      throw std::invalid_argument(
          "native whole-arm checkpoint root cannot be empty");

    struct PendingArmMarch {
      std::string name;
      std::vector<std::vector<std::string>> basis_handles;
      std::vector<json::object> integrand_rows;
      std::vector<std::vector<std::shared_ptr<StoredLocalBase>>> basis;
      std::vector<std::string> match_handles;
      std::vector<std::string> local_handles;
      std::vector<std::string> row_local_handles;
      std::string aggregate_handle;
    };
    const auto parse_pending_arm = [&](const char* name) {
      const auto& raw_arm = as_object(root.at(name),
                                      "native whole-arm arm request");
      require_exact_keys(raw_arm, {"receiving_basis", "integrand_rows"},
                         "native whole-arm arm request");
      PendingArmMarch arm;
      arm.name = name;
      const auto& raw_basis_sets = as_array(
          raw_arm.at("receiving_basis"),
          "native whole-arm receiving basis sets");
      arm.basis_handles.reserve(raw_basis_sets.size());
      for (const auto& raw_set : raw_basis_sets) {
        const auto& values = as_array(
            raw_set, "native whole-arm receiving basis set");
        if (values.empty())
          throw std::invalid_argument(
              "native whole-arm receiving basis sets cannot be empty");
        std::set<std::string> unique;
        std::vector<std::string> handles;
        handles.reserve(values.size());
        for (const auto& raw_handle : values) {
          if (!raw_handle.is_string() || raw_handle.as_string().empty())
            throw std::invalid_argument(
                "native whole-arm basis handles must be nonempty strings");
          std::string handle(raw_handle.as_string());
          if (!unique.insert(handle).second)
            throw std::invalid_argument(
                "native whole-arm basis handles must be pairwise distinct");
          handles.push_back(std::move(handle));
        }
        arm.basis_handles.push_back(std::move(handles));
      }
      const auto& raw_rows = as_array(
          raw_arm.at("integrand_rows"),
          "native whole-arm integrand rows");
      arm.integrand_rows.reserve(raw_rows.size());
      for (const auto& raw_row : raw_rows)
        arm.integrand_rows.push_back(
            as_object(raw_row, "native whole-arm integrand row"));
      return arm;
    };
    std::array<PendingArmMarch, 2> arms{
        parse_pending_arm("lower"), parse_pending_arm("upper")};

    const auto plan_handle = required_string(root, "tile_plan");
    const auto anchor_handle = required_string(root, "anchor");
    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> anchor;
    std::string combined_handle;
    const std::size_t total_matches =
        arms[0].basis_handles.size() + arms[1].basis_handles.size();
    std::size_t total_tiles = 0;
    std::size_t retained_local_reservation = 0;
    constexpr std::size_t published_line_results = 3;
    bool reservation_live = false;
    {
      // Resolve every source token and reserve the complete publication set
      // before either worker exists.  Subsequent public releases cannot
      // invalidate the acquired strong owners.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto plan_found = session->tile_plans.find(plan_handle);
      if (plan_found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or released tile plan for whole-arm marching");
      plan = plan_found->second;
      if (required_string(root, "tile_plan_checkpoint_identity") !=
          plan->checkpoint_identity())
        throw std::invalid_argument(
            "whole-arm tile-plan checkpoint token is stale");
      const auto anchor_found = session->locals.find(anchor_handle);
      if (anchor_found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released anchor local for whole-arm marching");
      anchor = anchor_found->second;
      if (required_string(root, "anchor_checkpoint_identity") !=
          anchor->checkpoint_identity())
        throw std::invalid_argument(
            "whole-arm anchor checkpoint token is stale");

      for (auto& arm : arms) {
        const auto& retained = plan->arm(arm.name);
        if (retained.exact.matches.size() != arm.basis_handles.size() ||
            arm.integrand_rows.size() != retained.exact.tiles.size() ||
            retained.exact.tiles.size() != arm.basis_handles.size() + 1)
          throw std::invalid_argument(
              "whole-arm basis/row counts do not reproduce the retained plan topology for " +
              arm.name);
        total_tiles = checked_diagnostic_sum(
            total_tiles, retained.exact.tiles.size(),
            "whole-arm tile count");
        arm.basis.reserve(arm.basis_handles.size());
        for (const auto& handles : arm.basis_handles) {
          std::vector<std::shared_ptr<StoredLocalBase>> resolved;
          resolved.reserve(handles.size());
          for (const auto& handle : handles) {
            const auto found = session->locals.find(handle);
            if (found == session->locals.end())
              throw std::invalid_argument(
                  "unknown or released native local in whole-arm receiving basis: " +
                  handle);
            resolved.push_back(found->second);
          }
          arm.basis.push_back(std::move(resolved));
        }
      }

      if (total_matches > session->match_capacity -
                              std::min(session->match_capacity,
                                       session->matches.size() +
                                           session->pending_matches))
        throw std::invalid_argument(
            "persistent local match capacity is exhausted by whole-arm marching");
      retained_local_reservation = checked_diagnostic_sum(
          total_matches, total_tiles,
          "whole-arm retained local reservation");
      if (retained_local_reservation > session->local_capacity -
                              std::min(session->local_capacity,
                                       session->locals.size() +
                                           session->pending_local_solves))
        throw std::invalid_argument(
            "persistent local capacity is exhausted by whole-arm marching");
      if (published_line_results > session->line_result_capacity -
              std::min(session->line_result_capacity,
                       session->line_results.size() +
                           session->pending_line_integrations))
        throw std::invalid_argument(
            "persistent line-result capacity is exhausted by whole-arm marching");

      for (auto& arm : arms) {
        arm.match_handles.reserve(arm.basis_handles.size());
        arm.local_handles.reserve(arm.basis_handles.size());
        arm.row_local_handles.reserve(arm.integrand_rows.size());
        for (std::size_t index = 0; index < arm.basis_handles.size(); ++index) {
          arm.match_handles.push_back(
              "m:" + std::to_string(session->next_match++));
          arm.local_handles.push_back(
              "l:" + std::to_string(session->next_local++));
        }
        for (std::size_t index = 0; index < arm.integrand_rows.size(); ++index)
          arm.row_local_handles.push_back(
              "l:" + std::to_string(session->next_local++));
        arm.aggregate_handle =
            "line:" + std::to_string(session->next_line_result++);
      }
      combined_handle =
          "line:" + std::to_string(session->next_line_result++);
      session->pending_matches += total_matches;
      session->pending_local_solves += retained_local_reservation;
      session->pending_line_integrations += published_line_results;
      reservation_live = true;
    }
    const auto release_reservation = [&]() {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (!reservation_live) return;
      if (session->pending_matches < total_matches ||
          session->pending_local_solves < retained_local_reservation ||
          session->pending_line_integrations < published_line_results)
        throw std::logic_error(
            "native whole-arm reservation accounting underflow");
      session->pending_matches -= total_matches;
      session->pending_local_solves -= retained_local_reservation;
      session->pending_line_integrations -= published_line_results;
      reservation_live = false;
    };
    struct WholeArmReservationGuard final {
      decltype(release_reservation)& release;
      bool& live;

      ~WholeArmReservationGuard() noexcept {
        if (!live) return;
        try {
          release();
        } catch (...) {
          // Reservation underflow is an internal invariant failure.  It is
          // unsafe to continue unwinding with counters that may still admit
          // work beyond the configured capacities.
          std::terminate();
        }
      }
    } reservation_guard{release_reservation, reservation_live};

    struct CompletedArmMarch {
      RetainedArmMarchResult march;
      std::vector<std::shared_ptr<StoredLineResult>> tile_lines;
      std::vector<std::shared_ptr<StoredLocalBase>> tile_sources;
      std::shared_ptr<StoredLineResult> aggregate;
      double elapsed_ms = 0.0;
    };
    std::array<CompletedArmMarch, 2> completed;
    std::array<std::exception_ptr, 2> failures;
    std::atomic<std::size_t> active_workers{0};
    std::atomic<std::size_t> max_active_workers{0};
    std::mutex start_mutex;
    std::condition_variable start_changed;
    std::size_t workers_ready = 0;
    bool workers_start = false;
    bool workers_cancel = false;
    const auto operation_started = std::chrono::steady_clock::now();
    std::unique_ptr<AcbPrecisionLease> whole_arm_acb_lease;
    if (session->domain == "acb") {
      whole_arm_acb_lease =
          std::make_unique<AcbPrecisionLease>(session->precision_bits);
      ComplexBall::set_precision(session->precision_bits);
    }
    const auto active_session_configuration_identity =
        checkpoint_configuration_identity(*session);

    const auto update_max_active = [&](std::size_t candidate) {
      auto observed = max_active_workers.load();
      while (observed < candidate &&
             !max_active_workers.compare_exchange_weak(observed, candidate)) {
      }
    };
    const auto run_arm = [&](std::size_t arm_index) {
      auto active = active_workers.fetch_add(1) + 1;
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
        const auto started = std::chrono::steady_clock::now();
        auto& input = arms[arm_index];
        auto& output = completed[arm_index];
        RetainedArmMarchInput march_input{
            input.name, input.basis_handles, input.basis,
            input.match_handles, input.local_handles};
        output.march = march_retained_arm(
            session->domain, session->precision_bits,
            active_session_configuration_identity, plan, anchor,
            march_input, work_epsilon, match_required_complete_max,
            refinement, checkpoint_root, true);
        TransportObservableContractionInput observable;
        observable.identity = checkpoint_root + ":" + input.name;
        observable.checkpoint_identity =
            checkpoint_root + ":" + input.name + ":aggregate";
        observable.checkpoint_root = checkpoint_root;
        observable.rows = input.integrand_rows;
        observable.epsilon = {work_epsilon, required_complete_max};
        observable.epsilon_record = json::object{
            {"min", work_epsilon.min_power},
            {"max", work_epsilon.complete_max},
            {"required_complete_max", required_complete_max}};
        observable.tail_policy = certify_tail
            ? TransportTailPolicy::Attempt
            : TransportTailPolicy::Stored;
        observable.projected_local_handles = input.row_local_handles;
        observable.aggregate_handle = input.aggregate_handle;
        observable.aggregate_record = json::object{
            {"kind", "complete-retained-arm"},
            {"combination", "sum-physical-tiles"},
            {"epsilon_contract", raw_epsilon},
            {"certify_tail_requested", certify_tail}};
        auto contractions = contract_transport_arm(
            session->domain, session->precision_bits, plan, input.name,
            output.march.tile_sources,
            std::vector<TransportObservableContractionInput>{observable});
        if (contractions.size() != 1)
          throw std::logic_error(
              "whole-arm contraction kernel returned the wrong result count");
        auto contracted = std::move(contractions.front());
        output.tile_sources = std::move(contracted.projected);
        output.tile_lines = std::move(contracted.tile_lines);
        output.aggregate = std::move(contracted.aggregate);
        output.elapsed_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - started).count();
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
      // No worker result has entered a session registry.  Prefer the lower
      // failure only for deterministic reporting when both arms fail.
      std::rethrow_exception(failures[0] ? failures[0] : failures[1]);
    }

    std::vector<std::shared_ptr<StoredLineResult>> arm_lines{
        completed[0].aggregate, completed[1].aggregate};
    std::vector<std::shared_ptr<StoredLocalBase>> combined_owners =
        completed[0].tile_sources;
    combined_owners.insert(combined_owners.end(),
                           completed[1].tile_sources.begin(),
                           completed[1].tile_sources.end());
    const auto combined_started = std::chrono::steady_clock::now();
    auto combined = build_retained_line_aggregate(
        combined_handle, checkpoint_root + ":combined", "combined",
        json::object{
            {"from_exact", plan->arm("lower").exact.to.str()},
            {"to_exact", plan->arm("upper").exact.to.str()}},
        json::object{
            {"kind", "complete-lower-to-upper-line"},
            {"combination", "negative-lower-anchor-arm-plus-upper-anchor-arm"},
            {"epsilon_contract", raw_epsilon},
            {"certify_tail_requested", certify_tail}},
        plan, std::move(combined_owners), arm_lines,
        std::vector<std::int32_t>{-1, 1},
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - combined_started).count());

    try {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches < total_matches ||
          session->pending_local_solves < retained_local_reservation ||
          session->pending_line_integrations < published_line_results)
        throw std::logic_error(
            "native whole-arm reservation accounting underflow");
      session->pending_matches -= total_matches;
      session->pending_local_solves -= retained_local_reservation;
      session->pending_line_integrations -= published_line_results;
      reservation_live = false;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during whole-arm marching");

      session->locals.reserve(session->locals.size() + 2);
      session->line_results.reserve(
          session->line_results.size() + published_line_results);
      std::vector<std::string> inserted_locals;
      std::vector<std::string> inserted_lines;
      try {
        for (std::size_t arm_index = 0; arm_index < 2; ++arm_index) {
          const auto& final_local =
              completed[arm_index].march.final_local();
          const auto existing = session->locals.find(final_local->handle());
          if (existing == session->locals.end()) {
            if (!session->locals.emplace(
                    final_local->handle(), final_local).second)
              throw std::logic_error(
                  "whole-arm final-local handle collision at publication");
            inserted_locals.push_back(final_local->handle());
          } else if (existing->second.get() != final_local.get()) {
            throw std::logic_error(
                "whole-arm final-local handle names a different retained object");
          }
          const auto& line = completed[arm_index].aggregate;
          if (!session->line_results.emplace(line->handle(), line).second)
            throw std::logic_error(
                "whole-arm aggregate handle collision at publication");
          inserted_lines.push_back(line->handle());
        }
        if (!session->line_results.emplace(combined->handle(), combined).second)
          throw std::logic_error(
              "whole-arm combined handle collision at publication");
        inserted_lines.push_back(combined->handle());
      } catch (...) {
        for (const auto& handle : inserted_lines)
          session->line_results.erase(handle);
        for (const auto& handle : inserted_locals)
          session->locals.erase(handle);
        throw;
      }

      for (std::size_t arm_index = 0; arm_index < 2; ++arm_index) {
        for (const auto& match : completed[arm_index].march.matches) {
          ++session->total_local_matches;
          session->total_local_match_ms += match->elapsed_ms();
          plan->note_match_advance(arms[arm_index].name);
        }
        for (const auto& line : completed[arm_index].tile_lines) {
          ++session->total_line_integrations;
          session->total_line_integration_ms += line->elapsed_ms();
          plan->note_integration();
        }
      }
    } catch (...) {
      if (reservation_live) release_reservation();
      throw;
    }

    json::object arm_response;
    for (std::size_t index = 0; index < 2; ++index) {
      auto final_local = completed[index].march.final_local()->summary();
      final_local["session"] = session->handle;
      auto line = completed[index].aggregate->summary();
      line["session"] = session->handle;
      arm_response[arms[index].name] = json::object{
          {"final_local", std::move(final_local)},
          {"line_result", std::move(line)},
          {"matches", completed[index].march.matches.size()},
          {"tiles", completed[index].tile_lines.size()},
          {"elapsed_ms", completed[index].elapsed_ms}};
    }
    auto combined_summary = combined->summary();
    combined_summary["session"] = session->handle;
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"capability", kRetainedParallelArmCapability},
        {"native_retained", true}, {"json_coefficients", 0},
        {"atomic_publication", true},
        {"workers", 2},
        {"max_parallel_arms", max_active_workers.load()},
        {"worker_overlap", max_active_workers.load() == 2},
        {"checkpoint_policy", checkpoint_policy},
        {"epsilon", raw_epsilon},
        {"arms", std::move(arm_response)},
        {"combined_line_result", std::move(combined_summary)},
        {"elapsed_ms", std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - operation_started).count()}};
  }

  if (operation == "tile.endpoint_limit") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "tile_plan",
         "tile_plan_checkpoint_identity", "arm", "local",
         "source_checkpoint_identity", "checkpoint_identity",
         "cancellation"},
        "native tile.endpoint_limit request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native plan-bound endpoint evaluation requires rational or Acb coefficients");
    const auto plan_handle = required_string(root, "tile_plan");
    const auto local_handle = required_string(root, "local");
    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> local;
    std::string endpoint_handle;
    {
      // Admission acquires strong ownership of both immutable dependencies.
      // Releasing either public token after this point cannot invalidate the
      // retained endpoint result or its checkpoint closure.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto plan_found = session->tile_plans.find(plan_handle);
      if (plan_found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or released tile plan for plan-bound endpoint limit");
      const auto local_found = session->locals.find(local_handle);
      if (local_found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released final local for plan-bound endpoint limit");
      if (session->endpoints.size() + session->pending_endpoint_limits >=
          session->endpoint_capacity)
        throw std::invalid_argument(
            "persistent endpoint result capacity is exhausted");
      plan = plan_found->second;
      local = local_found->second;
      endpoint_handle = "e:" +
          std::to_string(session->next_endpoint++);
      ++session->pending_endpoint_limits;
    }

    std::shared_ptr<StoredEndpointResult> endpoint;
    try {
      endpoint = build_planned_endpoint_limit(
          endpoint_handle, root, plan, local);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_endpoint_limits == 0)
        throw std::logic_error(
            "native plan-bound endpoint reservation accounting underflow");
      --session->pending_endpoint_limits;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_endpoint_limits == 0)
        throw std::logic_error(
            "native plan-bound endpoint reservation accounting underflow");
      --session->pending_endpoint_limits;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during plan-bound endpoint limit");
      session->endpoints.emplace(endpoint_handle, endpoint);
      ++session->total_endpoint_limits;
      session->total_endpoint_limit_ms += endpoint->elapsed_ms();
    }
    auto response = endpoint->summary();
    response["status"] = "ok";
    response["session"] = session->handle;
    return response;
  }

  if (operation == "integration.line") {
    if (root.if_contains("certify_tail") != nullptr)
      require_exact_keys(root,
          {"schema", "op", "session", "tile_plan", "local", "arm",
           "tile", "epsilon", "source_checkpoint_identity",
           "tile_plan_checkpoint_identity", "checkpoint_identity",
           "certify_tail"},
          "native integration.line request");
    else
      require_exact_keys(root,
          {"schema", "op", "session", "tile_plan", "local", "arm",
           "tile", "epsilon", "source_checkpoint_identity",
           "tile_plan_checkpoint_identity", "checkpoint_identity"},
          "native integration.line request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native planned line integration requires rational or Acb coefficients");
    const auto plan_handle = required_string(root, "tile_plan");
    const auto local_handle = required_string(root, "local");
    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> local;
    std::string line_handle;
    {
      // The admitted operation strongly owns both dependencies.  Lower and
      // upper calls take this lock only for admission/publication and execute
      // their immutable plan arms concurrently outside it.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto plan_found = session->tile_plans.find(plan_handle);
      if (plan_found == session->tile_plans.end())
        throw std::invalid_argument(
            "unknown or released tile plan for line integration");
      const auto local_found = session->locals.find(local_handle);
      if (local_found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released local for line integration");
      if (session->line_results.size() +
              session->pending_line_integrations >=
          session->line_result_capacity)
        throw std::invalid_argument(
            "persistent native line-result capacity is exhausted");
      plan = plan_found->second;
      local = local_found->second;
      line_handle = "line:" +
          std::to_string(session->next_line_result++);
      ++session->pending_line_integrations;
    }
    std::shared_ptr<StoredLineResult> result;
    try {
      result = build_planned_line_result(
          line_handle, root, plan, local);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_line_integrations == 0)
        throw std::logic_error(
            "native line-integration reservation accounting underflow");
      --session->pending_line_integrations;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_line_integrations == 0)
        throw std::logic_error(
            "native line-integration reservation accounting underflow");
      --session->pending_line_integrations;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during native line integration");
      session->line_results.emplace(line_handle, result);
      ++session->total_line_integrations;
      session->total_line_integration_ms += result->elapsed_ms();
      plan->note_integration();
    }
    auto response = result->summary();
    response["status"] = "ok";
    response["session"] = session->handle;
    return response;
  }

  if (operation == "integration.stats" ||
      operation == "integration.export") {
    const auto line_handle = required_string(root, "line");
    std::shared_ptr<StoredLineResult> result;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->line_results.find(line_handle);
      if (found == session->line_results.end())
        throw std::invalid_argument(
            "unknown or released native line result");
      result = found->second;
    }
    json::object response;
    if (operation == "integration.stats") {
      response = result->stats_json();
    } else {
      const auto output_digits = root.if_contains("output_digits")
          ? static_cast<int>(as_i64(root.at("output_digits"),
                                    "line export output digits"))
          : session->output_digits;
      if (output_digits < 1)
        throw std::invalid_argument(
            "line export output digits must be positive");
      response = result->export_values(
          required_string(root, "checkpoint_identity"), output_digits);
      const auto elapsed = response.at("elapsed_ms").as_double();
      std::lock_guard<std::mutex> lock(session->mutex);
      ++session->total_line_exports;
      session->total_line_export_ms += elapsed;
    }
    response["status"] = "ok";
    response["session"] = session->handle;
    return response;
  }

  if (operation == "integration.release") {
    const auto line_handle = required_string(root, "line");
    std::shared_ptr<StoredLineResult> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->line_results.find(line_handle);
      if (found == session->line_results.end())
        throw std::invalid_argument(
            "unknown or already released native line result");
      removed = std::move(found->second);
      session->line_results.erase(found);
    }
    return json::object{
        {"status", "ok"}, {"released", line_handle},
        {"checkpoint_identity", removed->checkpoint_identity()}};
  }
