  if (operation == "chart.prepare") {
    const auto key = required_string(root, "key");
    const auto identity = required_string(root, "identity");
    const auto& problem = as_object(root.at("problem"), "prepared problem");
    const auto problem_domain = required_string(problem, "domain");
    if (problem_domain != session->domain)
      throw std::invalid_argument(
          "prepared problem domain differs from its solver session");
    if (session->domain == "acb" &&
        as_i64(problem.at("precision_bits"), "precision bits") !=
            session->precision_bits)
      throw std::invalid_argument(
          "prepared problem precision differs from its solver session");
    if (session->domain == "symbolic" &&
        parse_symbols(problem) != session->symbols)
      throw std::invalid_argument(
          "prepared problem regulator field differs from its solver session");
    const auto analytic = root.if_contains("analytic")
        ? root.at("analytic") : json::value(nullptr);
    std::optional<std::string> geometry_record;
    std::optional<std::string> principal_matrix_record;
    std::optional<std::string> native_scc_capabilities;
    std::optional<std::string> regular_value_relative_accuracy_max_exact;
    if (analytic.is_object()) {
      const auto& analytic_object = analytic.as_object();
      if (const auto* geometry = analytic_object.if_contains("geometry"))
        geometry_record = canonical_chart_geometry_record(*geometry);
      if (const auto* principal =
              analytic_object.if_contains("principal_matrix"))
        principal_matrix_record = parse_exact_parent_matrix(
            *principal, as_u32(problem.at("d"), "dimension"),
            "prepared chart principal matrix").canonical_record;
      if (const auto* capabilities =
              analytic_object.if_contains("native_scc_capabilities"))
        native_scc_capabilities =
            canonical_native_scc_capabilities(*capabilities);
      if (const auto* raw_accuracy = analytic_object.if_contains(
              "regular_value_relative_accuracy_max_exact")) {
        if (!raw_accuracy->is_string())
          throw std::invalid_argument(
              "prepared chart regular-value relative-accuracy contract must be an exact rational string");
        const auto accuracy_text = std::string(raw_accuracy->as_string());
        const Rational accuracy(accuracy_text);
        if (accuracy.sign() <= 0 || !(accuracy < Rational(1)) ||
            accuracy.str() != accuracy_text)
          throw std::invalid_argument(
              "prepared chart regular-value relative-accuracy contract must be a canonical exact rational strictly between zero and one");
        regular_value_relative_accuracy_max_exact = accuracy_text;
      }
    }
    json::object combined_analytic;
    combined_analytic["session"] = json::parse(session->analytic_identity);
    combined_analytic["chart"] = analytic;
    auto scc = validate_scc_certificate(
        root.at("scc"), as_u32(problem.at("d"), "dimension"));
    auto signature = static_problem_signature(
        problem, combined_analytic, scc, identity);

    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (const auto found = session->handles_by_key.find(key);
          found != session->handles_by_key.end()) {
        const auto chart = session->charts.at(found->second);
        if (chart->signature() != signature)
          throw std::invalid_argument(
              "persistent chart cache key collision with unequal exact identity");
        return json::object{{"status", "ok"}, {"chart", chart->handle()},
                            {"reused", true},
                            {"dimension", chart->dimension()},
                            {"frame_base", chart->frame_base()},
                            {"frame_width", chart->frame_width()},
                            {"d0_inverse_mode", chart->d0_inverse_mode()}};
      }
      if (session->charts.size() >= session->chart_capacity)
        throw std::invalid_argument("persistent chart capacity is exhausted");
    }
    std::string chart_handle;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      chart_handle = "c:" + std::to_string(session->next_chart++);
    }
    std::shared_ptr<PreparedChartBase> chart;
    if (session->domain == "rational")
      chart = parse_prepared_chart<Rational>(
          session, root, chart_handle, key, identity, geometry_record,
          principal_matrix_record, native_scc_capabilities,
          regular_value_relative_accuracy_max_exact,
          std::move(scc), std::move(signature));
    else if (session->domain == "acb")
      chart = parse_prepared_chart<ComplexBall>(
          session, root, chart_handle, key, identity, geometry_record,
          principal_matrix_record, native_scc_capabilities,
          regular_value_relative_accuracy_max_exact,
          std::move(scc), std::move(signature));
    else
      chart = parse_prepared_chart<SymbolicRational>(
          session, root, chart_handle, key, identity, geometry_record,
          principal_matrix_record, native_scc_capabilities,
          regular_value_relative_accuracy_max_exact,
          std::move(scc), std::move(signature));
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      // A concurrent duplicate prepare is harmless only when its complete
      // collision certificate is byte-identical.
      if (const auto found = session->handles_by_key.find(key);
          found != session->handles_by_key.end()) {
        const auto existing = session->charts.at(found->second);
        if (existing->signature() != chart->signature())
          throw std::invalid_argument(
              "concurrent chart cache key collision with unequal identity");
        chart = existing;
      } else {
        session->charts.emplace(chart->handle(), chart);
        session->handles_by_key.emplace(key, chart->handle());
      }
    }
    return json::object{{"status", "ok"}, {"chart", chart->handle()},
                        {"reused", chart->handle() != chart_handle},
                        {"dimension", chart->dimension()},
                        {"frame_base", chart->frame_base()},
                        {"frame_width", chart->frame_width()},
                        {"d0_inverse_mode", chart->d0_inverse_mode()},
                        {"scc_components", chart->scc().component_count},
                        {"scc_structural_edges",
                         chart->scc().structural_edges.size()},
                        {"scc_condensation_edges",
                         chart->scc().condensation_edges.size()},
                        {"scc_topological_order",
                         encode_indices(chart->scc().topological_order)},
                        {"scc_coupling_depth", chart->scc().coupling_depth}};
  }

  if (operation == "scc.prepare") {
    const auto key = required_string(root, "key");
    const auto identity = required_string(root, "identity");
    const auto signature = composite_scc_signature(root);
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (const auto found = session->scc_handles_by_key.find(key);
          found != session->scc_handles_by_key.end()) {
        const auto composite = session->sccs.at(found->second);
        if (composite->signature() != signature)
          throw std::invalid_argument(
              "persistent SCC cache key collision with unequal exact identity");
        auto result = composite->stats_json();
        result["status"] = "ok";
        result["session"] = session->handle;
        result["reused"] = true;
        return result;
      }
      if (session->sccs.size() >= session->scc_capacity)
        throw std::invalid_argument("persistent SCC capacity is exhausted");
    }

    const auto& raw_blocks = as_array(root.at("blocks"), "SCC blocks");
    std::vector<std::shared_ptr<PreparedChartBase>> erased_charts;
    erased_charts.reserve(raw_blocks.size());
    std::string scc_handle;
    {
      // Resolve and strongly retain the complete diagonal handle set before
      // leaving the session lock.  The composite remains valid if the public
      // chart handles are released after this preparation boundary.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      for (const auto& raw_block_value : raw_blocks) {
        const auto& raw_block = as_object(raw_block_value, "SCC block");
        const auto chart_handle = required_string(raw_block, "chart");
        const auto found = session->charts.find(chart_handle);
        if (found == session->charts.end())
          throw std::invalid_argument(
              "SCC preparation references an unknown or released chart");
        erased_charts.push_back(found->second);
      }
      scc_handle = "scc:" + std::to_string(session->next_scc++);
    }

    std::shared_ptr<CompositeSCCChartBase> composite;
    if (session->domain == "rational")
      composite = parse_composite_scc_chart<Rational>(
          session, root, scc_handle, key, identity, signature,
          erased_charts);
    else if (session->domain == "acb")
      composite = parse_composite_scc_chart<ComplexBall>(
          session, root, scc_handle, key, identity, signature,
          erased_charts);
    else
      composite = parse_composite_scc_chart<SymbolicRational>(
          session, root, scc_handle, key, identity, signature,
          erased_charts);

    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during SCC preparation");
      if (const auto found = session->scc_handles_by_key.find(key);
          found != session->scc_handles_by_key.end()) {
        const auto existing = session->sccs.at(found->second);
        if (existing->signature() != composite->signature())
          throw std::invalid_argument(
              "concurrent SCC cache key collision with unequal exact identity");
        composite = existing;
      } else {
        if (session->sccs.size() >= session->scc_capacity)
          throw std::invalid_argument(
              "persistent SCC capacity was exhausted during preparation");
        session->sccs.emplace(composite->handle(), composite);
        session->scc_handles_by_key.emplace(key, composite->handle());
      }
    }
    auto result = composite->stats_json();
    result["status"] = "ok";
    result["session"] = session->handle;
    result["reused"] = composite->handle() != scc_handle;
    return result;
  }

  if (operation == "chart.solve") {
    const auto chart_handle = required_string(root, "chart");
    std::shared_ptr<PreparedChartBase> chart;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->charts.find(chart_handle);
      if (found == session->charts.end())
        throw std::invalid_argument("unknown or released persistent chart");
      chart = found->second;
    }
    const auto output_digits = root.if_contains("output_digits")
        ? static_cast<int>(as_i64(root.at("output_digits"), "output digits"))
        : session->output_digits;
    if (output_digits < 1)
      throw std::invalid_argument("output digits must be positive");
    auto result = chart->solve(as_object(root.at("run"), "recurrence run"),
                               output_digits);
    result["session"] = session->handle;
    result["chart"] = chart->handle();
    return result;
  }

  if (operation == "chart.solve_batch") {
    const auto chart_handle = required_string(root, "chart");
    std::shared_ptr<PreparedChartBase> chart;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->charts.find(chart_handle);
      if (found == session->charts.end())
        throw std::invalid_argument("unknown or released persistent chart");
      chart = found->second;
    }
    const auto output_digits = root.if_contains("output_digits")
        ? static_cast<int>(as_i64(root.at("output_digits"), "output digits"))
        : session->output_digits;
    if (output_digits < 1)
      throw std::invalid_argument("output digits must be positive");
    const auto& runs = as_array(root.at("runs"), "persistent recurrence runs");
    const auto requested_threads = root.if_contains("threads")
        ? as_u32(root.at("threads"), "batch threads") : 1;
    if (requested_threads == 0)
      throw std::invalid_argument("batch threads must be positive");

    const bool symbolic_serialized = session->domain == "symbolic";
    const auto bounded_threads = std::min<std::size_t>(
        {static_cast<std::size_t>(requested_threads),
         static_cast<std::size_t>(kMaxPersistentBatchThreads), runs.size()});
    const auto worker_count = symbolic_serialized && bounded_threads != 0
        ? std::size_t{1} : bounded_threads;
    const auto started = std::chrono::steady_clock::now();
    std::vector<json::object> results(runs.size());
    std::atomic<std::size_t> next{0};
    auto worker = [&]() {
      while (true) {
        const auto index = next.fetch_add(1);
        if (index >= runs.size()) return;
        results[index] = solve_prepared_chart_safe(
            chart, runs[index], output_digits, session->handle);
      }
    };
    // jthread guarantees already-started workers are joined if a later
    // thread construction throws; destroying a joinable std::thread here
    // would otherwise terminate the host Wolfram kernel.
    std::vector<std::jthread> workers;
    workers.reserve(worker_count);
    for (std::size_t i = 0; i < worker_count; ++i)
      workers.emplace_back(worker);
    for (auto& thread : workers) thread.join();

    std::size_t succeeded = 0;
    json::array encoded;
    encoded.reserve(results.size());
    for (auto& result : results) {
      if (result.if_contains("status") != nullptr &&
          result.at("status") == "ok")
        ++succeeded;
      encoded.push_back(std::move(result));
    }
    const auto elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"chart", chart->handle()}, {"results", std::move(encoded)},
        {"attempted", runs.size()}, {"succeeded", succeeded},
        {"failed", runs.size() - succeeded},
        {"requested_threads", requested_threads},
        {"worker_threads", worker_count},
        {"thread_limit", kMaxPersistentBatchThreads},
        {"symbolic_serialized", symbolic_serialized},
        {"elapsed_ms", elapsed_ms}};
  }

  if (operation == "session.solve_many") {
    const auto& raw_jobs = as_array(
        root.at("jobs"), "persistent cross-chart recurrence jobs");
    const auto requested_threads = root.if_contains("threads")
        ? as_u32(root.at("threads"), "batch threads") : 1;
    if (requested_threads == 0)
      throw std::invalid_argument("batch threads must be positive");

    struct PendingJob {
      std::string chart_handle;
      const json::value* run = nullptr;
      int output_digits = 0;
    };
    std::vector<PendingJob> pending;
    pending.reserve(raw_jobs.size());
    for (std::size_t index = 0; index < raw_jobs.size(); ++index) {
      const auto& job = as_object(
          raw_jobs[index], "persistent cross-chart recurrence job");
      const auto* raw_run = job.if_contains("run");
      if (raw_run == nullptr)
        throw std::invalid_argument(
            "session.solve_many job " + std::to_string(index) +
            " is missing its complete run record");
      const auto output_digits = job.if_contains("output_digits")
          ? static_cast<int>(as_i64(
                job.at("output_digits"), "job output digits"))
          : session->output_digits;
      if (output_digits < 1)
        throw std::invalid_argument("job output digits must be positive");
      pending.push_back(PendingJob{
          required_string(job, "chart"), raw_run, output_digits});
    }

    struct ResolvedJob {
      std::shared_ptr<PreparedChartBase> chart;
      const json::value* run = nullptr;
      int output_digits = 0;
    };
    std::vector<ResolvedJob> jobs;
    jobs.reserve(pending.size());
    {
      // Resolve the complete handle set before any worker exists.  Apart
      // from making cross-session/released handles loud, the shared_ptrs
      // retain every selected chart for the whole batch even if a concurrent
      // chart.release follows this validation boundary.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      for (std::size_t index = 0; index < pending.size(); ++index) {
        const auto found = session->charts.find(pending[index].chart_handle);
        if (found == session->charts.end())
          throw std::invalid_argument(
              "unknown or released persistent chart in session.solve_many "
              "job " + std::to_string(index) + ": " +
              pending[index].chart_handle);
        jobs.push_back(ResolvedJob{
            found->second, pending[index].run, pending[index].output_digits});
      }
    }

    const bool symbolic_serialized = session->domain == "symbolic";
    const auto bounded_threads = std::min<std::size_t>(
        {static_cast<std::size_t>(requested_threads),
         static_cast<std::size_t>(kMaxPersistentBatchThreads), jobs.size()});
    const auto worker_count = symbolic_serialized && bounded_threads != 0
        ? std::size_t{1} : bounded_threads;
    const auto started = std::chrono::steady_clock::now();
    std::vector<json::object> results(jobs.size());
    std::atomic<std::size_t> next{0};
    auto worker = [&](std::stop_token stop) {
      while (!stop.stop_requested()) {
        const auto index = next.fetch_add(1);
        if (index >= jobs.size()) return;
        const auto& job = jobs[index];
        results[index] = solve_prepared_chart_safe(
            job.chart, *job.run, job.output_digits, session->handle);
      }
    };
    // If construction of a later worker fails, jthread destruction requests
    // stop and joins every worker that already started.  No joinable native
    // thread can escape into the host Wolfram kernel.
    std::vector<std::jthread> workers;
    workers.reserve(worker_count);
    for (std::size_t i = 0; i < worker_count; ++i)
      workers.emplace_back(worker);
    for (auto& thread : workers) thread.join();

    std::size_t succeeded = 0;
    json::array encoded;
    encoded.reserve(results.size());
    for (auto& result : results) {
      if (result.if_contains("status") != nullptr &&
          result.at("status") == "ok")
        ++succeeded;
      encoded.push_back(std::move(result));
    }
    const auto elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"results", std::move(encoded)}, {"attempted", jobs.size()},
        {"succeeded", succeeded}, {"failed", jobs.size() - succeeded},
        {"requested_threads", requested_threads},
        {"worker_threads", worker_count},
        {"thread_limit", kMaxPersistentBatchThreads},
        {"symbolic_serialized", symbolic_serialized},
        {"elapsed_ms", elapsed_ms}};
  }

  if (operation == "scc.solve_column") {
    const auto scc_handle = required_string(root, "scc");
    std::shared_ptr<CompositeSCCChartBase> composite;
    std::string local_handle;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->sccs.find(scc_handle);
      if (found == session->sccs.end())
        throw std::invalid_argument(
            "unknown or released persistent SCC chart");
      if (session->locals.size() + session->pending_local_solves >=
          session->local_capacity)
        throw std::invalid_argument("persistent local capacity is exhausted");
      composite = found->second;
      local_handle = "l:" + std::to_string(session->next_local++);
      ++session->pending_local_solves;
    }

    CompositeColumnSolveResult column;
    try {
      column = composite->solve_column(local_handle, root, composite);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "native SCC local reservation accounting underflow");
      --session->pending_local_solves;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "native SCC local reservation accounting underflow");
      --session->pending_local_solves;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during SCC column solve");
      session->locals.emplace(local_handle, column.local);
      const auto local_stats = column.local->stats();
      ++session->total_local_solves;
      ++session->total_scc_column_solves;
      session->total_local_run_parse_ms += local_stats.create_parse_ms;
      session->total_local_kernel_ms += local_stats.create_kernel_ms;
    }
    auto response = column.local->summary();
    response["status"] = "ok";
    response["session"] = session->handle;
    response["scc"] = composite->handle();
    response["native_retained"] = true;
    response["json_coefficients"] = 0;
    response["execution_capability"] =
        composite->column_execution_capability();
    response["block_diagnostics"] = std::move(column.block_diagnostics);
    response["elapsed_ms"] = column.elapsed_ms;
    return response;
  }

  if (operation == "scc.solve_columns") {
    require_exact_keys(root,
        {"schema", "op", "session", "scc", "columns", "threads"},
        "native SCC column-batch request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native SCC column batches require rational or Acb coefficients");
    const auto requested_threads = as_u32(
        root.at("threads"), "native SCC column-batch threads");
    if (requested_threads == 0)
      throw std::invalid_argument(
          "native SCC column-batch threads must be positive");
    const auto& raw_columns = as_array(
        root.at("columns"), "native SCC column-batch columns");
    if (raw_columns.empty())
      throw std::invalid_argument(
          "native SCC column batch cannot be empty");
    for (std::size_t index = 0; index < raw_columns.size(); ++index)
      require_exact_keys(
          as_object(raw_columns[index], "native SCC batch column"),
          {"seed", "targets", "checkpoint_identity"},
          "native SCC batch column");

    const auto scc_handle = required_string(root, "scc");
    std::shared_ptr<CompositeSCCChartBase> composite;
    std::vector<std::string> local_handles;
    local_handles.reserve(raw_columns.size());
    {
      // Reserve the complete ordered result set before starting workers.
      // Public SCC release cannot invalidate the strongly owned composite,
      // and capacity cannot be consumed between individual columns.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto found = session->sccs.find(scc_handle);
      if (found == session->sccs.end())
        throw std::invalid_argument(
            "unknown or released persistent SCC chart");
      if (raw_columns.size() > session->local_capacity ||
          session->locals.size() + session->pending_local_solves >
              session->local_capacity - raw_columns.size())
        throw std::invalid_argument(
            "persistent local capacity cannot admit the complete SCC column batch");
      composite = found->second;
      for (std::size_t index = 0; index < raw_columns.size(); ++index)
        local_handles.push_back(
            "l:" + std::to_string(session->next_local++));
      session->pending_local_solves += raw_columns.size();
    }

    const auto worker_count = std::min<std::size_t>(
        {static_cast<std::size_t>(requested_threads),
         static_cast<std::size_t>(kMaxPersistentBatchThreads),
         raw_columns.size()});
    const auto started = std::chrono::steady_clock::now();
    std::vector<CompositeColumnSolveResult> columns(raw_columns.size());
    std::vector<std::exception_ptr> errors(raw_columns.size());
    std::atomic<std::size_t> next{0};
    auto worker = [&](std::stop_token stop) {
      while (!stop.stop_requested()) {
        const auto index = next.fetch_add(1);
        if (index >= raw_columns.size()) return;
        try {
          columns[index] = composite->solve_column(
              local_handles[index],
              as_object(raw_columns[index], "native SCC batch column"),
              composite);
        } catch (...) {
          errors[index] = std::current_exception();
        }
      }
    };
    std::vector<std::jthread> workers;
    workers.reserve(worker_count);
    try {
      for (std::size_t index = 0; index < worker_count; ++index)
        workers.emplace_back(worker);
      for (auto& thread : workers) thread.join();
    } catch (...) {
      for (auto& thread : workers) thread.request_stop();
      for (auto& thread : workers)
        if (thread.joinable()) thread.join();
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves < raw_columns.size())
        throw std::logic_error(
            "native SCC column-batch reservation accounting underflow");
      session->pending_local_solves -= raw_columns.size();
      throw;
    }

    const auto failed = std::find_if(
        errors.begin(), errors.end(), [](const auto& error) {
          return error != nullptr;
        });
    if (failed != errors.end()) {
      const auto failure = *failed;
      {
        std::lock_guard<std::mutex> lock(session->mutex);
        if (session->pending_local_solves < raw_columns.size())
          throw std::logic_error(
              "native SCC column-batch reservation accounting underflow");
        session->pending_local_solves -= raw_columns.size();
      }
      // A batch is atomic at the retained-state boundary: successful worker
      // temporaries are discarded and the first ordered failure is loud.
      std::rethrow_exception(failure);
    }

    if (std::any_of(columns.begin(), columns.end(), [](const auto& column) {
          return column.local == nullptr;
        })) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves < raw_columns.size())
        throw std::logic_error(
            "native SCC column-batch reservation accounting underflow");
      session->pending_local_solves -= raw_columns.size();
      throw std::logic_error(
          "native SCC column batch completed without a retained local");
    }

    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves < raw_columns.size())
        throw std::logic_error(
            "native SCC column-batch reservation accounting underflow");
      session->pending_local_solves -= raw_columns.size();
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during SCC column batch");
      session->locals.reserve(session->locals.size() + columns.size());
      std::vector<std::string> inserted;
      inserted.reserve(columns.size());
      try {
        for (std::size_t index = 0; index < columns.size(); ++index) {
          if (!session->locals.emplace(
                  local_handles[index], columns[index].local).second)
            throw std::logic_error(
                "native SCC column batch produced a duplicate local handle");
          inserted.push_back(local_handles[index]);
        }
      } catch (...) {
        for (const auto& handle : inserted) session->locals.erase(handle);
        throw;
      }
      for (std::size_t index = 0; index < columns.size(); ++index) {
        const auto local_stats = columns[index].local->stats();
        ++session->total_local_solves;
        ++session->total_scc_column_solves;
        session->total_local_run_parse_ms += local_stats.create_parse_ms;
        session->total_local_kernel_ms += local_stats.create_kernel_ms;
      }
    }

    json::array responses;
    responses.reserve(columns.size());
    for (auto& column : columns) {
      auto response = column.local->summary();
      response["status"] = "ok";
      response["session"] = session->handle;
      response["scc"] = composite->handle();
      response["native_retained"] = true;
      response["json_coefficients"] = 0;
      response["execution_capability"] =
          composite->column_execution_capability();
      response["block_diagnostics"] =
          std::move(column.block_diagnostics);
      response["elapsed_ms"] = column.elapsed_ms;
      responses.push_back(std::move(response));
    }
    const auto elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"scc", composite->handle()},
        {"results", std::move(responses)},
        {"columns", raw_columns.size()},
        {"requested_threads", requested_threads},
        {"worker_threads", worker_count},
        {"thread_limit", kMaxPersistentBatchThreads},
        {"atomic_retention", true},
        {"json_coefficients", 0}, {"elapsed_ms", elapsed_ms}};
  }

  if (operation == "local.solve") {
    const auto chart_handle = required_string(root, "chart");
    std::shared_ptr<PreparedChartBase> chart;
    std::string local_handle;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->charts.find(chart_handle);
      if (found == session->charts.end())
        throw std::invalid_argument("unknown or released persistent chart");
      if (session->locals.size() + session->pending_local_solves >=
          session->local_capacity)
        throw std::invalid_argument("persistent local capacity is exhausted");
      chart = found->second;
      local_handle = "l:" + std::to_string(session->next_local++);
      ++session->pending_local_solves;
    }

    std::shared_ptr<StoredLocalBase> local;
    try {
      local = chart->solve_local(
          local_handle, as_object(root.at("run"), "recurrence run"),
          as_object(root.at("metadata"), "local metadata"), chart);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error("native local reservation accounting underflow");
      --session->pending_local_solves;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error("native local reservation accounting underflow");
      --session->pending_local_solves;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during local solve");
      session->locals.emplace(local_handle, local);
      const auto local_stats = local->stats();
      ++session->total_local_solves;
      session->total_local_run_parse_ms += local_stats.create_parse_ms;
      session->total_local_kernel_ms += local_stats.create_kernel_ms;
    }
    auto result = local->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    result["native_retained"] = true;
    result["json_coefficients"] = 0;
    return result;
  }

  if (operation == "local.apply_rational_row") {
    require_exact_keys(root,
        {"schema", "op", "session", "local", "row",
         "source_checkpoint_identity", "checkpoint_identity"},
        "native local.apply_rational_row request");
    if (session->domain == "symbolic")
      throw std::domain_error(
          "native rational-row application requires exact Rational or specialized Acb coefficients");
    const auto source_handle = required_string(root, "local");
    std::shared_ptr<StoredLocalBase> source;
    std::string local_handle;
    {
      // Materialization runs outside the registry lock, but owns the source
      // for its full lifetime.  The derived scalar local also retains that
      // source as provenance after publication.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto found = session->locals.find(source_handle);
      if (found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released retained local for rational-row application");
      if (session->locals.size() + session->pending_local_solves >=
          session->local_capacity)
        throw std::invalid_argument("persistent local capacity is exhausted");
      source = found->second;
      local_handle = "l:" + std::to_string(session->next_local++);
      ++session->pending_local_solves;
    }

    std::shared_ptr<StoredLocalBase> local;
    try {
      if (session->domain == "rational") {
        const auto typed =
            std::dynamic_pointer_cast<StoredLocal<Rational>>(source);
        if (!typed)
          throw std::logic_error(
              "rational-row source local differs from its Rational session");
        local = build_rational_row_local<Rational>(
            local_handle, root, session->precision_bits, typed, source);
      } else {
        const auto typed =
            std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(source);
        if (!typed)
          throw std::logic_error(
              "rational-row source local differs from its Acb session");
        local = build_rational_row_local<ComplexBall>(
            local_handle, root, session->precision_bits, typed, source);
      }
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "native rational-row reservation accounting underflow");
      --session->pending_local_solves;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "native rational-row reservation accounting underflow");
      --session->pending_local_solves;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during rational-row application");
      if (!session->locals.emplace(local_handle, local).second)
        throw std::logic_error(
            "native rational-row application produced a duplicate local handle");
    }
    auto result = local->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    result["source_local"] = source_handle;
    result["application_capability"] = kRetainedRationalRowCapability;
    result["native_retained"] = true;
    result["json_coefficients"] = 0;
    return result;
  }

  if (operation == "local.endpoint_limit") {
    if (root.if_contains("output_digits") != nullptr ||
        root.if_contains("include_coefficients") != nullptr)
      throw std::invalid_argument(
          "local.endpoint_limit never exports coefficients; use the explicit "
          "endpoint.export compatibility operation");
    const auto local_handle = required_string(root, "local");
    std::shared_ptr<StoredLocalBase> local;
    std::string endpoint_handle;
    {
      // Strong ownership is acquired under the session lock.  Public
      // local.release may remove the registry token after admission without
      // invalidating this endpoint computation.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto found = session->locals.find(local_handle);
      if (found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released native local for endpoint limit");
      if (session->endpoints.size() + session->pending_endpoint_limits >=
          session->endpoint_capacity)
        throw std::invalid_argument(
            "persistent endpoint result capacity is exhausted");
      local = found->second;
      endpoint_handle = "e:" + std::to_string(session->next_endpoint++);
      ++session->pending_endpoint_limits;
    }

    std::shared_ptr<StoredEndpointResult> endpoint;
    try {
      endpoint = build_endpoint_limit(endpoint_handle, root, local);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_endpoint_limits == 0)
        throw std::logic_error(
            "native endpoint reservation accounting underflow");
      --session->pending_endpoint_limits;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_endpoint_limits == 0)
        throw std::logic_error(
            "native endpoint reservation accounting underflow");
      --session->pending_endpoint_limits;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during endpoint limit");
      session->endpoints.emplace(endpoint_handle, endpoint);
      ++session->total_endpoint_limits;
      session->total_endpoint_limit_ms += endpoint->elapsed_ms();
    }
    auto result = endpoint->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "local.match") {
    if (session->domain != "rational")
      throw std::invalid_argument(
          "local.match currently supports only the exact rational regular "
          "local capability; Acb and symbolic saturation are not routed "
          "through it");
    const auto& raw_basis = as_array(
        root.at("basis"), "exact regular local match basis");
    if (raw_basis.empty())
      throw std::invalid_argument("local.match basis cannot be empty");
    std::vector<std::string> basis_handles;
    basis_handles.reserve(raw_basis.size());
    std::set<std::string> unique_handles;
    for (const auto& raw_handle : raw_basis) {
      if (!raw_handle.is_string())
        throw std::invalid_argument(
            "local.match basis handles must be strings");
      std::string handle(raw_handle.as_string());
      if (!unique_handles.insert(handle).second)
        throw std::invalid_argument(
            "local.match basis handles must be pairwise distinct");
      basis_handles.push_back(std::move(handle));
    }
    const auto incoming_handle = required_string(root, "incoming");
    if (unique_handles.contains(incoming_handle))
      throw std::invalid_argument(
          "local.match incoming handle must be distinct from its basis");

    std::vector<std::shared_ptr<StoredLocalBase>> basis;
    std::shared_ptr<StoredLocalBase> incoming;
    std::string match_handle;
    {
      // Resolve and strongly retain every local before releasing the session
      // lock.  A concurrent public local.release cannot invalidate an
      // already admitted native match operation.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      if (session->matches.size() + session->pending_matches >=
          session->match_capacity)
        throw std::invalid_argument(
            "persistent local match capacity is exhausted");
      for (const auto& handle : basis_handles) {
        const auto found = session->locals.find(handle);
        if (found == session->locals.end())
          throw std::invalid_argument(
              "unknown or released native local in match basis: " + handle);
        basis.push_back(found->second);
      }
      const auto found = session->locals.find(incoming_handle);
      if (found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released incoming native local: " +
            incoming_handle);
      incoming = found->second;
      match_handle = "m:" + std::to_string(session->next_match++);
      ++session->pending_matches;
    }

    std::shared_ptr<StoredExactRegularMatch> match;
    try {
      match = build_exact_regular_match(
          match_handle, root, basis_handles, basis, incoming_handle,
          incoming);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches == 0)
        throw std::logic_error(
            "native local match reservation accounting underflow");
      --session->pending_matches;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches == 0)
        throw std::logic_error(
            "native local match reservation accounting underflow");
      --session->pending_matches;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during local matching");
      session->matches.emplace(match_handle, match);
      ++session->total_local_matches;
      session->total_local_match_ms += match->elapsed_ms();
    }
    auto result = match->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "local.match_acb") {
    if (session->domain != "acb")
      throw std::invalid_argument(
          "local.match_acb requires an Acb persistent session; exact rational matching remains local.match v1");
    const auto& raw_basis = as_array(
        root.at("basis"), "refined Acb local match basis");
    if (raw_basis.empty())
      throw std::invalid_argument(
          "local.match_acb basis cannot be empty");
    std::vector<std::string> basis_handles;
    basis_handles.reserve(raw_basis.size());
    std::set<std::string> unique_handles;
    for (const auto& raw_handle : raw_basis) {
      if (!raw_handle.is_string())
        throw std::invalid_argument(
            "local.match_acb basis handles must be strings");
      std::string handle(raw_handle.as_string());
      if (!unique_handles.insert(handle).second)
        throw std::invalid_argument(
            "local.match_acb basis handles must be pairwise distinct");
      basis_handles.push_back(std::move(handle));
    }
    const auto incoming_handle = required_string(root, "incoming");
    if (unique_handles.contains(incoming_handle))
      throw std::invalid_argument(
          "local.match_acb incoming handle must be distinct from its basis");

    std::vector<std::shared_ptr<StoredLocalBase>> basis;
    std::shared_ptr<StoredLocalBase> incoming;
    std::string match_handle;
    {
      // Admission takes strong ownership of every source before releasing
      // the session lock.  Concurrent public releases cannot invalidate the
      // exact-point evaluation or the bounded refinement operation.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      if (session->matches.size() + session->pending_matches >=
          session->match_capacity)
        throw std::invalid_argument(
            "persistent local match capacity is exhausted");
      for (const auto& handle : basis_handles) {
        const auto found = session->locals.find(handle);
        if (found == session->locals.end())
          throw std::invalid_argument(
              "unknown or released native local in Acb match basis: " +
              handle);
        basis.push_back(found->second);
      }
      const auto found = session->locals.find(incoming_handle);
      if (found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released incoming native Acb local: " +
            incoming_handle);
      incoming = found->second;
      match_handle = "m:" + std::to_string(session->next_match++);
      ++session->pending_matches;
    }

    std::shared_ptr<StoredRefinedAcbMatch> match;
    try {
      match = build_refined_acb_match(
          match_handle, root, basis_handles, basis, incoming_handle,
          incoming, session->precision_bits,
          checkpoint_configuration_identity(*session));
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches == 0)
        throw std::logic_error(
            "native Acb local match reservation accounting underflow");
      --session->pending_matches;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_matches == 0)
        throw std::logic_error(
            "native Acb local match reservation accounting underflow");
      --session->pending_matches;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during Acb local matching");
      session->matches.emplace(match_handle, match);
      ++session->total_local_matches;
      session->total_local_match_ms += match->elapsed_ms();
    }
    auto result = match->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "match.materialize_local") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "match", "checkpoint_identity"},
        "native match.materialize_local request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native plan-match local materialization requires rational or Acb coefficients");
    const auto match_handle = required_string(root, "match");
    const auto checkpoint_identity = required_string(
        root, "checkpoint_identity");
    std::shared_ptr<StoredPlannedMatchHop> match;
    std::string local_handle;
    {
      // Admission strongly retains the complete handoff before releasing the
      // session lock.  The finite Laurent combination then runs natively and
      // independently of public match/plan/basis tokens.
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      const auto found = session->matches.find(match_handle);
      if (found == session->matches.end())
        throw std::invalid_argument(
            "unknown or released retained match for local materialization");
      match = std::dynamic_pointer_cast<StoredPlannedMatchHop>(found->second);
      if (!match)
        throw std::invalid_argument(
            "match.materialize_local requires a plan-driven match handoff");
      if (session->locals.size() + session->pending_local_solves >=
          session->local_capacity)
        throw std::invalid_argument(
            "persistent local capacity is exhausted");
      local_handle = "l:" + std::to_string(session->next_local++);
      ++session->pending_local_solves;
    }

    std::shared_ptr<StoredLocalBase> local;
    try {
      local = match->materialize(
          local_handle, checkpoint_identity, session->precision_bits, match);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "native match materialization reservation accounting underflow");
      --session->pending_local_solves;
      throw;
    }
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "native match materialization reservation accounting underflow");
      --session->pending_local_solves;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during match materialization");
      session->locals.emplace(local_handle, local);
    }
    auto result = local->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    result["materialization_capability"] =
        kRetainedPlannedMatchMaterializationCapability;
    result["native_retained"] = true;
    result["json_coefficients"] = 0;
    return result;
  }

  if (operation == "local.evaluate") {
    const auto local_handle = required_string(root, "local");
    std::shared_ptr<StoredLocalBase> local;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->locals.find(local_handle);
      if (found == session->locals.end())
        throw std::invalid_argument("unknown or released native local solution");
      local = found->second;
    }
    const auto output_digits = root.if_contains("output_digits")
        ? static_cast<int>(as_i64(root.at("output_digits"), "output digits"))
        : session->output_digits;
    if (output_digits < 1)
      throw std::invalid_argument("output digits must be positive");
    auto result = local->evaluate(root, output_digits);
    result["status"] = "ok";
    result["session"] = session->handle;
    result["local"] = local->handle();
    result["chart"] = local->source_chart();
    return result;
  }

  if (operation == "local.certify_residual") {
    const auto local_handle = required_string(root, "local");
    std::shared_ptr<StoredLocalBase> local;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->locals.find(local_handle);
      if (found == session->locals.end())
        throw std::invalid_argument(
            "unknown or released native local solution");
      local = found->second;
    }
    const auto output_digits = root.if_contains("output_digits")
        ? static_cast<int>(as_i64(root.at("output_digits"), "output digits"))
        : session->output_digits;
    if (output_digits < 1)
      throw std::invalid_argument("output digits must be positive");
    auto result = local->certify_residual(root, output_digits);
    result["status"] = "ok";
    result["session"] = session->handle;
    result["local"] = local->handle();
    result["chart"] = local->source_chart();
    return result;
  }

  if (operation == "endpoint.stats") {
    const auto endpoint_handle = required_string(root, "endpoint");
    std::shared_ptr<StoredEndpointResult> endpoint;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->endpoints.find(endpoint_handle);
      if (found == session->endpoints.end())
        throw std::invalid_argument(
            "unknown or released native endpoint result");
      endpoint = found->second;
    }
    auto result = endpoint->stats_json();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "endpoint.export") {
    const auto endpoint_handle = required_string(root, "endpoint");
    std::shared_ptr<StoredEndpointResult> endpoint;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->endpoints.find(endpoint_handle);
      if (found == session->endpoints.end())
        throw std::invalid_argument(
            "unknown or released native endpoint result");
      endpoint = found->second;
    }
    const auto output_digits = root.if_contains("output_digits")
        ? static_cast<int>(as_i64(root.at("output_digits"), "output digits"))
        : session->output_digits;
    if (output_digits < 1)
      throw std::invalid_argument("output digits must be positive");
    auto result = endpoint->export_values(
        required_string(root, "checkpoint_identity"), output_digits);
    const auto export_ms = result.at("elapsed_ms").as_double();
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      ++session->total_endpoint_exports;
      session->total_endpoint_export_ms += export_ms;
    }
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "endpoint.release") {
    const auto endpoint_handle = required_string(root, "endpoint");
    std::shared_ptr<StoredEndpointResult> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->endpoints.find(endpoint_handle);
      if (found == session->endpoints.end())
        throw std::invalid_argument(
            "unknown or already released native endpoint result");
      removed = std::move(found->second);
      session->endpoints.erase(found);
    }
    return json::object{{"status", "ok"}, {"released", endpoint_handle},
                        {"checkpoint_identity",
                         removed->checkpoint_identity()}};
  }

  if (operation == "match.stats") {
    const auto match_handle = required_string(root, "match");
    std::shared_ptr<StoredMatchBase> match;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->matches.find(match_handle);
      if (found == session->matches.end())
        throw std::invalid_argument(
            "unknown or released native local match");
      match = found->second;
    }
    auto result = match->summary();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "match.release") {
    const auto match_handle = required_string(root, "match");
    std::shared_ptr<StoredMatchBase> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->matches.find(match_handle);
      if (found == session->matches.end())
        throw std::invalid_argument(
            "unknown or already released native local match");
      removed = std::move(found->second);
      session->matches.erase(found);
    }
    return json::object{{"status", "ok"}, {"released", match_handle}};
  }

  if (operation == "local.release") {
    const auto local_handle = required_string(root, "local");
    std::shared_ptr<StoredLocalBase> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->locals.find(local_handle);
      if (found == session->locals.end())
        throw std::invalid_argument(
            "unknown or already released native local solution");
      removed = std::move(found->second);
      session->locals.erase(found);
    }
    return json::object{{"status", "ok"}, {"released", local_handle},
                        {"chart", removed->source_chart()}};
  }

  if (operation == "local.stats") {
    const auto local_handle = required_string(root, "local");
    std::shared_ptr<StoredLocalBase> local;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->locals.find(local_handle);
      if (found == session->locals.end())
        throw std::invalid_argument("unknown or released native local solution");
      local = found->second;
    }
    auto result = local->stats_json();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "scc.stats") {
    const auto scc_handle = required_string(root, "scc");
    std::shared_ptr<CompositeSCCChartBase> composite;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->sccs.find(scc_handle);
      if (found == session->sccs.end())
        throw std::invalid_argument(
            "unknown or released persistent SCC chart");
      composite = found->second;
    }
    auto result = composite->stats_json();
    result["status"] = "ok";
    result["session"] = session->handle;
    return result;
  }

  if (operation == "scc.release") {
    const auto scc_handle = required_string(root, "scc");
    std::shared_ptr<CompositeSCCChartBase> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->sccs.find(scc_handle);
      if (found == session->sccs.end())
        throw std::invalid_argument(
            "unknown or already released persistent SCC chart");
      removed = std::move(found->second);
      session->scc_handles_by_key.erase(removed->key());
      session->sccs.erase(found);
    }
    return json::object{{"status", "ok"}, {"released", scc_handle}};
  }

  if (operation == "chart.release") {
    const auto chart_handle = required_string(root, "chart");
    std::shared_ptr<PreparedChartBase> removed;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->charts.find(chart_handle);
      if (found == session->charts.end())
        throw std::invalid_argument("unknown or already released chart");
      removed = found->second;
      session->handles_by_key.erase(removed->key());
      session->charts.erase(found);
    }
    return json::object{{"status", "ok"}, {"released", chart_handle}};
  }

  if (operation == "session.stats") {
    std::vector<std::shared_ptr<PreparedChartBase>> charts;
    std::vector<std::shared_ptr<StoredLocalBase>> locals;
    std::vector<std::shared_ptr<StoredMatchBase>> matches;
    std::vector<std::shared_ptr<StoredEndpointResult>> endpoints;
    std::vector<std::shared_ptr<StoredTilePlan>> tile_plans;
    std::vector<std::shared_ptr<StoredTransportArmState>> transport_states;
    std::vector<std::shared_ptr<StoredLineResult>> line_results;
    std::vector<std::shared_ptr<CompositeSCCChartBase>> sccs;
    std::size_t pending_local_solves = 0;
    std::size_t pending_matches = 0;
    std::size_t pending_endpoint_limits = 0;
    std::size_t pending_tile_plans = 0;
    std::size_t pending_transport_states = 0;
    std::size_t pending_line_integrations = 0;
    std::size_t transport_pair_streams = 0;
    std::uint64_t total_local_solves = 0;
    std::uint64_t total_scc_column_solves = 0;
    std::uint64_t total_local_matches = 0;
    std::uint64_t checkpoint_generation = 0;
    std::uint64_t checkpoint_restore_count = 0;
    std::string restored_from_checkpoint_identity;
    std::uint64_t total_endpoint_limits = 0;
    std::uint64_t total_endpoint_exports = 0;
    std::uint64_t total_tile_plans = 0;
    std::uint64_t total_transport_arm_marches = 0;
    std::uint64_t total_transport_contractions = 0;
    std::uint64_t total_transport_observables = 0;
    std::uint64_t total_transport_pair_contractions = 0;
    std::uint64_t total_transport_pair_observables = 0;
    std::uint64_t total_transport_endpoint_batches = 0;
    std::uint64_t total_transport_endpoint_rows = 0;
    std::uint64_t total_line_integrations = 0;
    std::uint64_t total_line_exports = 0;
    double total_local_run_parse_ms = 0.0, total_local_kernel_ms = 0.0;
    double total_local_match_ms = 0.0;
    double total_endpoint_limit_ms = 0.0;
    double total_endpoint_export_ms = 0.0;
    double total_tile_plan_ms = 0.0;
    double total_transport_arm_ms = 0.0;
    double total_transport_contraction_ms = 0.0;
    double total_transport_pair_contraction_ms = 0.0;
    double total_transport_endpoint_batch_ms = 0.0;
    double total_line_integration_ms = 0.0;
    double total_line_export_ms = 0.0;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      for (const auto& [ignored, chart] : session->charts)
        charts.push_back(chart);
      for (const auto& [ignored, local] : session->locals)
        locals.push_back(local);
      for (const auto& [ignored, match] : session->matches)
        matches.push_back(match);
      for (const auto& [ignored, endpoint] : session->endpoints)
        endpoints.push_back(endpoint);
      for (const auto& [ignored, plan] : session->tile_plans)
        tile_plans.push_back(plan);
      for (const auto& [ignored, state] : session->transport_states)
        transport_states.push_back(state);
      for (const auto& [ignored, result] : session->line_results)
        line_results.push_back(result);
      for (const auto& [ignored, composite] : session->sccs)
        sccs.push_back(composite);
      pending_local_solves = session->pending_local_solves;
      pending_matches = session->pending_matches;
      pending_endpoint_limits = session->pending_endpoint_limits;
      pending_tile_plans = session->pending_tile_plans;
      pending_transport_states = session->pending_transport_states;
      pending_line_integrations = session->pending_line_integrations;
      transport_pair_streams = session->transport_pair_streams.size();
      total_local_solves = session->total_local_solves;
      total_scc_column_solves = session->total_scc_column_solves;
      total_local_matches = session->total_local_matches;
      total_endpoint_limits = session->total_endpoint_limits;
      total_endpoint_exports = session->total_endpoint_exports;
      total_tile_plans = session->total_tile_plans;
      total_transport_arm_marches = session->total_transport_arm_marches;
      total_transport_contractions = session->total_transport_contractions;
      total_transport_observables = session->total_transport_observables;
      total_transport_pair_contractions =
          session->total_transport_pair_contractions;
      total_transport_pair_observables =
          session->total_transport_pair_observables;
      total_transport_endpoint_batches =
          session->total_transport_endpoint_batches;
      total_transport_endpoint_rows =
          session->total_transport_endpoint_rows;
      total_line_integrations = session->total_line_integrations;
      total_line_exports = session->total_line_exports;
      total_local_run_parse_ms = session->total_local_run_parse_ms;
      total_local_kernel_ms = session->total_local_kernel_ms;
      total_local_match_ms = session->total_local_match_ms;
      checkpoint_generation = session->checkpoint_generation;
      checkpoint_restore_count = session->checkpoint_restore_count;
      restored_from_checkpoint_identity =
          session->restored_from_checkpoint_identity;
      total_endpoint_limit_ms = session->total_endpoint_limit_ms;
      total_endpoint_export_ms = session->total_endpoint_export_ms;
      total_tile_plan_ms = session->total_tile_plan_ms;
      total_transport_arm_ms = session->total_transport_arm_ms;
      total_transport_contraction_ms =
          session->total_transport_contraction_ms;
      total_transport_pair_contraction_ms =
          session->total_transport_pair_contraction_ms;
      total_transport_endpoint_batch_ms =
          session->total_transport_endpoint_batch_ms;
      total_line_integration_ms = session->total_line_integration_ms;
      total_line_export_ms = session->total_line_export_ms;
    }
    std::uint64_t runs = 0;
    double prepare_parse_ms = 0.0, run_parse_ms = 0.0, kernel_ms = 0.0;
    json::array chart_stats;
    for (const auto& chart : charts) {
      const auto stats = chart->stats();
      runs += stats.runs;
      prepare_parse_ms += stats.prepare_parse_ms;
      run_parse_ms += stats.run_parse_ms;
      kernel_ms += stats.kernel_ms;
      chart_stats.push_back(json::object{
          {"chart", chart->handle()}, {"key", chart->key()},
          {"dimension", chart->dimension()},
          {"frame_base", chart->frame_base()},
          {"frame_width", chart->frame_width()}, {"runs", stats.runs},
          {"d0_inverse_mode", chart->d0_inverse_mode()},
          {"local_solves", stats.local_runs},
          {"scc_components", chart->scc().component_count},
          {"scc_structural_edges", chart->scc().structural_edges.size()},
          {"scc_condensation_edges", chart->scc().condensation_edges.size()},
          {"scc_topological_order",
           encode_indices(chart->scc().topological_order)},
          {"scc_coupling_depth", chart->scc().coupling_depth},
          {"prepare_parse_ms", stats.prepare_parse_ms},
          {"run_parse_ms", stats.run_parse_ms},
          {"kernel_ms", stats.kernel_ms},
          {"local_run_parse_ms", stats.local_run_parse_ms},
          {"local_kernel_ms", stats.local_kernel_ms}});
    }
    std::uint64_t local_evaluations = 0;
    std::uint64_t local_residual_certifications = 0;
    std::uint64_t local_endpoint_limits = 0;
    std::uint64_t local_line_integrations = 0;
    std::uint64_t local_tail_certificate_requests = 0;
    std::uint64_t local_tail_certificate_certified = 0;
    std::uint64_t local_tail_certificate_inconclusive = 0;
    std::uint64_t local_tail_certificate_unsupported = 0;
    std::size_t local_coefficients = 0;
    double local_evaluate_ms = 0.0, local_residual_certify_ms = 0.0;
    double local_endpoint_limit_ms = 0.0;
    double local_line_integration_ms = 0.0;
    json::array local_stats;
    for (const auto& local : locals) {
      const auto stats = local->stats();
      local_evaluations += stats.evaluations;
      local_residual_certifications += stats.residual_certifications;
      local_endpoint_limits += stats.endpoint_limits;
      local_line_integrations += stats.line_integrations;
      local_tail_certificate_requests += stats.tail_certificate_requests;
      local_tail_certificate_certified += stats.tail_certificate_certified;
      local_tail_certificate_inconclusive +=
          stats.tail_certificate_inconclusive;
      local_tail_certificate_unsupported +=
          stats.tail_certificate_unsupported;
      local_coefficients += stats.coefficient_count;
      local_evaluate_ms += stats.evaluate_ms;
      local_residual_certify_ms += stats.residual_certify_ms;
      local_endpoint_limit_ms += stats.endpoint_limit_ms;
      local_line_integration_ms += stats.line_integration_ms;
      auto encoded = local->stats_json();
      local_stats.push_back(std::move(encoded));
    }
    json::array scc_stats;
    for (const auto& composite : sccs)
      scc_stats.push_back(composite->stats_json());
    json::array match_stats;
    for (const auto& match : matches)
      match_stats.push_back(match->summary());
    json::array endpoint_stats;
    for (const auto& endpoint : endpoints)
      endpoint_stats.push_back(endpoint->stats_json());
    json::array tile_plan_stats;
    for (const auto& plan : tile_plans)
      tile_plan_stats.push_back(plan->summary(false));
    json::array transport_state_stats;
    for (const auto& state : transport_states)
      transport_state_stats.push_back(state->stats_json());
    json::array line_result_stats;
    for (const auto& result : line_results)
      line_result_stats.push_back(result->stats_json());
    return json::object{{"status", "ok"}, {"session", session->handle},
                        {"domain", session->domain},
                        {"precision_bits", session->precision_bits},
                        {"charts", charts.size()}, {"runs", runs},
                        {"locals", locals.size()},
                        {"matches", matches.size()},
                        {"endpoints", endpoints.size()},
                        {"tile_plans", tile_plans.size()},
                        {"transport_states", transport_states.size()},
                        {"line_results", line_results.size()},
                        {"transport_pair_streams", transport_pair_streams},
                        {"scc_charts", sccs.size()},
                        {"pending_local_solves", pending_local_solves},
                        {"pending_matches", pending_matches},
                        {"pending_endpoint_limits", pending_endpoint_limits},
                        {"pending_tile_plans", pending_tile_plans},
                        {"pending_transport_states",
                         pending_transport_states},
                        {"pending_line_integrations",
                         pending_line_integrations},
                        {"local_solves", total_local_solves},
                        {"local_matches", total_local_matches},
                        {"endpoint_limits", total_endpoint_limits},
                        {"endpoint_exports", total_endpoint_exports},
                        {"tile_plans_created", total_tile_plans},
                        {"transport_arm_marches",
                         total_transport_arm_marches},
                        {"transport_contractions",
                         total_transport_contractions},
                        {"transport_observables",
                         total_transport_observables},
                        {"transport_pair_contractions",
                         total_transport_pair_contractions},
                        {"transport_pair_observables",
                         total_transport_pair_observables},
                        {"transport_endpoint_batches",
                         total_transport_endpoint_batches},
                        {"transport_endpoint_rows",
                         total_transport_endpoint_rows},
                        {"line_integrations", total_line_integrations},
                        {"line_exports", total_line_exports},
                        {"local_match_capability",
                         session->domain == "rational"
                             ? kExactRegularLocalMatchCapability
                             : "unsupported"},
                        {"acb_local_match_capability",
                         session->domain == "acb"
                             ? kRefinedAcbLocalMatchCapability
                             : "unsupported"},
                        {"endpoint_limit_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedEndpointLimitCapability},
                        {"planned_endpoint_limit_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedPlannedEndpointLimitCapability},
                        {"tile_plan_capability", kRetainedTilePlanCapability},
                        {"single_arm_tile_plan_capability",
                         kRetainedSingleArmTilePlanCapability},
                        {"planned_match_hop_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedPlannedMatchHopCapability},
                        {"planned_match_materialization_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedPlannedMatchMaterializationCapability},
                        {"rational_row_application_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedRationalRowCapability},
                        {"line_integration_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedStoredLineCapability},
                        {"parallel_transport_arm_state_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedParallelTransportArmStateCapability},
                        {"transport_arm_state_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedTransportArmStateCapability},
                        {"transport_arm_contraction_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedTransportArmContractionCapability},
                        {"transport_pair_contraction_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedTransportPairContractionCapability},
                        {"transport_pair_stream_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedTransportPairStreamCapability},
                        {"transport_endpoint_batch_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedTransportEndpointBatchCapability},
                        {"certified_tail_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRegularTailMajorantCapability},
                        {"certified_line_integration_capability",
                         session->domain == "symbolic"
                             ? "unsupported"
                             : kRetainedCertifiedLineCapability},
                        {"scc_column_solves", total_scc_column_solves},
                        {"local_evaluations", local_evaluations},
                        {"local_residual_certifications",
                         local_residual_certifications},
                        {"local_endpoint_limits", local_endpoint_limits},
                        {"local_line_integrations",
                         local_line_integrations},
                        {"local_tail_certificate_requests",
                         local_tail_certificate_requests},
                        {"local_tail_certificate_certified",
                         local_tail_certificate_certified},
                        {"local_tail_certificate_inconclusive",
                         local_tail_certificate_inconclusive},
                        {"local_tail_certificate_unsupported",
                         local_tail_certificate_unsupported},
                        {"local_coefficient_count", local_coefficients},
                        {"static_tensor_copies", 0},
                        {"prepare_parse_ms", prepare_parse_ms},
                        {"run_parse_ms", run_parse_ms},
                        {"kernel_ms", kernel_ms},
                        {"local_run_parse_ms", total_local_run_parse_ms},
                        {"local_kernel_ms", total_local_kernel_ms},
                        {"local_evaluate_ms", local_evaluate_ms},
                        {"local_residual_certify_ms",
                         local_residual_certify_ms},
                        {"local_endpoint_limit_ms", local_endpoint_limit_ms},
                        {"local_line_integration_ms",
                         local_line_integration_ms},
                        {"local_match_ms", total_local_match_ms},
                        {"checkpoint_generation", checkpoint_generation},
                        {"checkpoint_restore_count",
                         checkpoint_restore_count},
                        {"restored_from_checkpoint_identity",
                         restored_from_checkpoint_identity.empty()
                             ? json::value(nullptr)
                             : json::value(
                                   restored_from_checkpoint_identity)},
                        {"endpoint_limit_ms", total_endpoint_limit_ms},
                        {"endpoint_export_ms", total_endpoint_export_ms},
                        {"tile_plan_ms", total_tile_plan_ms},
                        {"transport_arm_ms", total_transport_arm_ms},
                        {"transport_contraction_ms",
                         total_transport_contraction_ms},
                        {"transport_pair_contraction_ms",
                         total_transport_pair_contraction_ms},
                        {"transport_endpoint_batch_ms",
                         total_transport_endpoint_batch_ms},
                        {"line_integration_ms", total_line_integration_ms},
                        {"line_export_ms", total_line_export_ms},
                        {"chart_stats", std::move(chart_stats)},
                        {"local_stats", std::move(local_stats)},
                        {"match_stats", std::move(match_stats)},
                        {"endpoint_stats", std::move(endpoint_stats)},
                        {"tile_plan_stats", std::move(tile_plan_stats)},
                        {"transport_state_stats",
                         std::move(transport_state_stats)},
                        {"line_result_stats", std::move(line_result_stats)},
                        {"scc_stats", std::move(scc_stats)}};
  }

  throw std::invalid_argument("unknown persistent solver operation: " + operation);
}
