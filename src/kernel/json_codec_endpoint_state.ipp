template <typename Scalar>
PreparedSparseLocalMultiplierMatrix<Scalar> parse_prepared_rational_row(
    const json::object& row, const LocalSolution<Scalar>& source,
    std::optional<std::int32_t> projected_complete_cap = std::nullopt) {
  require_exact_keys(row,
      {"schema", "columns", "exact_identity", "entries"},
      "prepared rational local row");
  if (required_string(row, "schema") !=
      "diffexp3-prepared-rational-local-row-v1")
    throw std::invalid_argument(
        "unsupported prepared rational local-row schema");
  const auto columns = as_u32(
      row.at("columns"), "prepared rational local-row columns");
  if (columns == 0 || columns != source.dimension)
    throw std::invalid_argument(
        "prepared rational local-row dimension differs from its retained local");

  PreparedSparseLocalMultiplierMatrix<Scalar> matrix;
  matrix.rows = 1;
  matrix.columns = columns;
  matrix.exact_identity = required_string(row, "exact_identity");
  if (matrix.exact_identity.empty())
    throw std::invalid_argument(
        "prepared rational local-row identity must be nonempty");

  std::optional<std::uint32_t> previous_column;
  for (const auto& raw_entry : as_array(
           row.at("entries"), "prepared rational local-row entries")) {
    const auto& entry = as_object(
        raw_entry, "prepared rational local-row entry");
    require_exact_keys(entry, {"column", "multiplier"},
                       "prepared rational local-row entry");
    const auto column = as_u32(
        entry.at("column"), "prepared rational local-row column");
    if (column >= columns ||
        (previous_column.has_value() && *previous_column >= column))
      throw std::invalid_argument(
          "prepared rational local-row columns must be unique, in range, and strictly increasing");
    previous_column = column;

    const auto& raw_multiplier = as_object(
        entry.at("multiplier"), "prepared rational local-row multiplier");
    const auto available_width = [&]() {
      if (const auto* analytic =
              raw_multiplier.if_contains("analytic_coefficients"))
        return as_array(
            *analytic,
            "prepared rational local-row analytic coefficients").size();
      return as_array(
          raw_multiplier.at("kernels"),
          "prepared rational local-row epsilon kernels").size();
    }();
    // With no requested cap, consume the full honest intersection of the
    // retained local and prepared multiplier prefixes.  A multiplier is not
    // required to be as wide as its source.  Explicit caps remain strict
    // below so a caller cannot silently receive less than it requested.
    auto epsilon_width =
        std::min(source.epsilon.width(), available_width);
    if (projected_complete_cap.has_value()) {
      const auto shift = as_i32(
          raw_multiplier.at("epsilon_shift"),
          "prepared rational local-row multiplier epsilon shift");
      const auto term_min = local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(source.epsilon.min_power) + shift,
          "prepared rational local-row projected minimum");
      const auto relative_needed =
          static_cast<std::int64_t>(*projected_complete_cap) - term_min;
      const auto required_width = relative_needed < 0
          ? std::size_t{1}
          : static_cast<std::size_t>(relative_needed) + 1;
      epsilon_width = std::min(source.epsilon.width(), required_width);
      if (available_width < epsilon_width)
        throw std::invalid_argument(
            "prepared rational local-row multiplier consumer prefix is too short"
            "; source_epsilon=[" +
            std::to_string(source.epsilon.min_power) + "," +
            std::to_string(source.epsilon.complete_max) + "]" +
            "; epsilon_shift=" + std::to_string(shift) +
            "; actual_width=" + std::to_string(available_width) +
            "; required_width=" + std::to_string(epsilon_width) +
            "; projected_complete_cap=" +
            std::to_string(*projected_complete_cap));
    }
    auto multiplier = parse_prepared_rational_taylor_multiplier<Scalar>(
        raw_multiplier, epsilon_width, source.taylor_width(), false,
        "prepared rational local-row multiplier");
    if (multiplier.proven_zero)
      throw std::invalid_argument(
          "structurally zero rational-row entries must be omitted");
    matrix.entries.push_back(
        typename PreparedSparseLocalMultiplierMatrix<Scalar>::Entry{
            0, column, std::move(multiplier)});
  }
  return matrix;
}

template <typename Scalar>
PreparedSparseLocalMultiplierMatrix<Scalar> parse_prepared_rational_row(
    const json::value& raw, const LocalSolution<Scalar>& source,
    std::optional<std::int32_t> projected_complete_cap = std::nullopt) {
  return parse_prepared_rational_row<Scalar>(
      as_object(raw, "prepared rational local row"), source,
      projected_complete_cap);
}

template <typename Scalar>
LocalSolution<Scalar> exact_zero_scalar_local_like(
    const LocalSolution<Scalar>& source,
    const std::string& checkpoint_identity) {
  auto result = local_algebra_detail::with_selected_component(source, 0);
  for (auto& sector : result.sectors)
    std::fill(sector.coefficients.begin(), sector.coefficients.end(),
              ScalarTraits<Scalar>::zero());
  result.checkpoint_identity = checkpoint_identity;
  validate_local_solution(result, false);
  return result;
}

std::int32_t shifted_local_validity(std::int32_t validity,
                                    std::int32_t shift) {
  if (validity == kCompleteInfinity) return kCompleteInfinity;
  return local_algebra_detail::checked_i32(
      static_cast<std::int64_t>(validity) + shift,
      "rational-row output validity");
}

template <typename Scalar>
std::shared_ptr<StoredLocalBase> build_rational_row_local(
    const std::string& local_handle, const json::object& request,
    slong precision_bits,
    const std::shared_ptr<StoredLocal<Scalar>>& source,
    const std::shared_ptr<StoredLocalBase>& erased_source,
    bool prepare_line_tail = true,
    std::optional<std::int32_t> projection_complete_max = std::nullopt,
    const json::object* prepared_row_override = nullptr) {
  if (source == nullptr || erased_source == nullptr ||
      source.get() != erased_source.get())
    throw std::logic_error(
        "retained rational-row application lost typed source ownership");
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  const auto source_checkpoint_identity = required_string(
      request, "source_checkpoint_identity");
  if (checkpoint_identity.empty() || source_checkpoint_identity.empty())
    throw std::invalid_argument(
        "retained rational-row checkpoint identities must be nonempty");
  if (source_checkpoint_identity != source->checkpoint_identity())
    throw std::invalid_argument(
        "rational-row source checkpoint identity differs from its retained local");
  if (!source->solution().error.empty())
    throw std::domain_error(
        "native rational-row application requires explicit source error-envelope propagation");

  std::unique_ptr<AcbPrecisionLease> acb_lease;
  if constexpr (std::is_same_v<Scalar, ComplexBall>) {
    acb_lease = std::make_unique<AcbPrecisionLease>(precision_bits);
    ComplexBall::set_precision(precision_bits);
  }
  const auto started = std::chrono::steady_clock::now();
  auto matrix = prepared_row_override != nullptr
      ? parse_prepared_rational_row<Scalar>(
            *prepared_row_override, source->solution(),
            projection_complete_max)
      : parse_prepared_rational_row<Scalar>(
            request.at("row"), source->solution(),
            projection_complete_max);

  json::array entry_provenance;
  entry_provenance.reserve(matrix.entries.size());
  std::int32_t output_top_valid = matrix.entries.empty()
      ? source->top_valid() : kCompleteInfinity;
  for (const auto& entry : matrix.entries) {
    output_top_valid = std::min(
        output_top_valid,
        shifted_local_validity(source->top_valid(),
                               entry.multiplier.epsilon_shift));
    entry_provenance.push_back(json::object{
        {"column", entry.column},
        {"epsilon_shift", entry.multiplier.epsilon_shift},
        {"center_pole_order", entry.multiplier.center_pole_order},
        {"exact_identity", entry.multiplier.exact_identity}});
  }

  auto applied = projection_complete_max.has_value()
      ? apply_prepared_scalar_row_window(
            matrix, source->solution(), *projection_complete_max,
            checkpoint_identity)
      : apply_prepared_sparse_local_matrix(
            matrix, source->solution(), checkpoint_identity);
  auto solution = applied.has_value()
      ? std::move(*applied)
      : exact_zero_scalar_local_like(
            source->solution(), checkpoint_identity);
  if (projection_complete_max.has_value() &&
      solution.epsilon.complete_max > *projection_complete_max)
    solution = restrict_local_epsilon_frame_strict_lower(
        solution, solution.epsilon.min_power, *projection_complete_max,
        checkpoint_identity);
  output_top_valid = std::min(
      output_top_valid, solution.epsilon.complete_max);
  if (output_top_valid < solution.epsilon.min_power)
    throw std::domain_error(
        "rational-row application has no valid output epsilon coefficient");
  if (output_top_valid < solution.epsilon.complete_max)
    solution = restrict_local_epsilon_frame_strict_lower(
        solution, solution.epsilon.min_power, output_top_valid,
        checkpoint_identity);
  if (solution.dimension != 1)
    throw std::logic_error(
        "rational-row application did not produce a scalar local solution");

  json::object derivation{
      {"schema", "diffexp3-retained-rational-row-local-application-v1"},
      {"capability", kRetainedRationalRowCapability},
      {"source", json::object{
           {"local", source->handle()},
           {"chart", source->source_chart()},
           {"source_operator_identity",
            source->source_operator_identity()},
           {"checkpoint_identity", source_checkpoint_identity},
           {"dimension", source->solution().dimension},
           {"epsilon", json::object{
                {"min", source->solution().epsilon.min_power},
                {"max", source->solution().epsilon.complete_max}}},
           {"taylor_complete_max",
            source->solution().taylor_complete_max}}},
      {"row", json::object{
           {"exact_identity", matrix.exact_identity},
           {"columns", matrix.columns},
           {"active_entries", std::move(entry_provenance)},
           {"structurally_zero", matrix.entries.empty()}}},
      {"output", json::object{
           {"checkpoint_identity", checkpoint_identity},
           {"dimension", solution.dimension},
           {"epsilon", json::object{
                {"min", solution.epsilon.min_power},
                {"max", solution.epsilon.complete_max}}},
           {"taylor_complete_max", solution.taylor_complete_max}}},
      {"analytic_prescriptions", "preserved-exactly"},
      {"coefficient_transport", "native-retained-only"}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(derivation));
  derivation["provenance_identity"] = provenance_identity;

  json::object operator_provenance{
      {"schema", "diffexp3-rational-row-derived-operator-v1"},
      {"source_operator_identity", source->source_operator_identity()},
      {"row_exact_identity", matrix.exact_identity},
      {"provenance_identity", provenance_identity}};
  const auto derived_operator_identity = json::serialize(
      canonical_json_value(operator_provenance));
  std::optional<RationalRowLineTailModelResult<Scalar>> rational_row_tail;
  if (prepare_line_tail)
    rational_row_tail = derive_rational_row_line_tail_model(
        matrix, source->solution(), source->tail_model(), solution);
  const auto elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  NativeLocalDiagnostics diagnostics;
  diagnostics.top_valid = output_top_valid;
  diagnostics.kernel_ms = elapsed_ms;
  auto retained = make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
      local_handle, source->source_chart(), derived_operator_identity,
      std::move(solution), precision_bits,
      std::vector<PseudoHit<Scalar>>{}, diagnostics, std::nullopt,
      std::move(derivation),
      std::static_pointer_cast<void>(erased_source));
  if (rational_row_tail.has_value())
    retained->attach_rational_row_line_tail_model(
        std::move(*rational_row_tail));
  return retained;
}

struct ParsedEndpointLimitPolicy {
  EndpointLimitOptions options;
  std::string cancellation_mode;
  std::optional<std::int32_t> requested_rim;
};

struct ParsedEndpointCancellationPolicy {
  bool allow_certified_numeric_cancellation = false;
  std::string cancellation_mode;
};

ParsedEndpointCancellationPolicy parse_endpoint_cancellation_policy(
    const json::object& request) {
  ParsedEndpointCancellationPolicy parsed;
  const auto& cancellation = as_object(
      request.at("cancellation"), "endpoint cancellation policy");
  parsed.cancellation_mode = required_string(cancellation, "mode");
  if (parsed.cancellation_mode == "exact-coefficient-field") {
    parsed.allow_certified_numeric_cancellation = false;
  } else if (parsed.cancellation_mode == "exact-or-acb-singleton") {
    parsed.allow_certified_numeric_cancellation = true;
  } else {
    throw std::invalid_argument(
        "endpoint cancellation mode must be exact-coefficient-field or "
        "exact-or-acb-singleton");
  }
  if (cancellation.size() != 1)
    throw std::invalid_argument(
        "endpoint cancellation policy accepts only its exact mode; "
        "tolerance-based cancellation is unsupported");
  return parsed;
}

ParsedEndpointLimitPolicy parse_endpoint_limit_policy(
    const json::object& request) {
  ParsedEndpointLimitPolicy parsed;
  parsed.options.approach_direction = as_i32(
      request.at("approach_direction"), "endpoint approach direction");
  if (parsed.options.approach_direction != 1 &&
      parsed.options.approach_direction != -1)
    throw std::invalid_argument(
        "endpoint approach direction must be exactly +1 or -1");

  if (const auto* raw_rim = request.if_contains("rim");
      raw_rim != nullptr && !raw_rim->is_null()) {
    const auto rim = as_i32(*raw_rim, "endpoint rim");
    if (rim != 1 && rim != -1)
      throw std::invalid_argument("endpoint rim must be exactly +1 or -1");
    parsed.requested_rim = rim;
    parsed.options.imaginary_sign = rim;
  }

  const auto cancellation = parse_endpoint_cancellation_policy(request);
  parsed.cancellation_mode = cancellation.cancellation_mode;
  parsed.options.allow_certified_numeric_cancellation =
      cancellation.allow_certified_numeric_cancellation;
  return parsed;
}

EpsilonWindow endpoint_value_window(const EndpointLimitResult& result) {
  if (result.values.empty())
    throw std::logic_error("retained endpoint result has no components");
  const auto window = result.values.front().window();
  for (const auto& component : result.values)
    if (component.window().min_power != window.min_power ||
        component.window().complete_max != window.complete_max)
      throw std::logic_error(
          "retained endpoint components have unequal epsilon windows");
  return window;
}

EpsilonVector endpoint_values_vector(const EndpointLimitResult& result) {
  const auto window = endpoint_value_window(result);
  EpsilonVector vector;
  vector.epsilon = window;
  vector.dimension = static_cast<std::uint32_t>(result.values.size());
  vector.coefficients.reserve(window.width() * vector.dimension);
  for (std::size_t ei = 0; ei < window.width(); ++ei)
    for (const auto& component : result.values)
      vector.coefficients.push_back(component.coefficients().at(ei));
  return vector;
}

class StoredEndpointResult {
 public:
  StoredEndpointResult(
      std::string handle, std::string checkpoint_identity,
      std::string provenance_identity, std::string source_local,
      std::string source_chart, std::string source_operator_identity,
      std::string source_checkpoint,
      std::string source_domain, std::int32_t approach_direction,
      std::optional<std::int32_t> requested_rim,
      std::string cancellation_mode, json::object analytic_metadata,
      EndpointLimitResult&& result, double elapsed_ms,
      std::optional<json::object> planned_source = std::nullopt,
      std::optional<std::int32_t> planned_effective_rim = std::nullopt,
      std::shared_ptr<StoredTilePlan> plan_owner = nullptr,
      std::shared_ptr<StoredLocalBase> local_owner = nullptr,
      std::shared_ptr<StoredTransportArmState> transport_owner = nullptr,
      std::optional<json::object> prepared_row_checkpoint = std::nullopt)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        source_local_(std::move(source_local)),
        source_chart_(std::move(source_chart)),
        source_operator_identity_(std::move(source_operator_identity)),
        source_checkpoint_(std::move(source_checkpoint)),
        source_domain_(std::move(source_domain)),
        approach_direction_(approach_direction),
        requested_rim_(requested_rim),
        cancellation_mode_(std::move(cancellation_mode)),
        analytic_metadata_(std::move(analytic_metadata)),
        result_(std::move(result)), elapsed_ms_(elapsed_ms),
        planned_source_(std::move(planned_source)),
        planned_effective_rim_(planned_effective_rim),
        plan_owner_(std::move(plan_owner)),
        local_owner_(std::move(local_owner)),
        transport_owner_(std::move(transport_owner)),
        prepared_row_checkpoint_(std::move(prepared_row_checkpoint)) {
    // Validate the retained public frame once, before publishing its handle.
    (void)endpoint_value_window(result_);
    const bool transport_bound = transport_owner_ != nullptr;
    if (transport_bound
            ? (!planned_source_.has_value() || plan_owner_ != nullptr ||
               local_owner_ != nullptr)
            : planned_source_.has_value()
            ? (plan_owner_ == nullptr || local_owner_ == nullptr)
            : (plan_owner_ != nullptr || local_owner_ != nullptr))
      throw std::invalid_argument(
          "retained endpoint ownership mode is inconsistent");
    if (prepared_row_checkpoint_.has_value() &&
        !(transport_bound || rolling_transport_bound()))
      throw std::invalid_argument(
          "transport endpoint prepared-row checkpoint ownership is inconsistent");
  }

  const std::string& handle() const { return handle_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& provenance_identity() const {
    return provenance_identity_;
  }
  double elapsed_ms() const { return elapsed_ms_; }
  bool plan_bound() const {
    return planned_source_.has_value() && !transport_bound();
  }
  bool transport_bound() const { return transport_owner_ != nullptr; }
  bool direct_plan_bound() const { return plan_bound(); }
  bool rolling_transport_bound() const {
    if (!plan_bound()) return false;
    const auto* mode = planned_source_->if_contains("binding_mode");
    return mode != nullptr && mode->is_string() &&
        mode->as_string() == "rolling-terminal-local";
  }
  bool derived_plan_bound() const { return planned_source_.has_value(); }
  const std::shared_ptr<StoredTilePlan>& plan_owner() const {
    return plan_owner_;
  }
  const std::shared_ptr<StoredLocalBase>& local_owner() const {
    return local_owner_;
  }
  const std::shared_ptr<StoredTransportArmState>& transport_owner() const {
    return transport_owner_;
  }

  json::object summary() const {
    const auto window = endpoint_value_window(result_);
    const auto cancellation_scope = source_domain_ == "rational"
        ? "exact-rational"
        : cancellation_mode_ == "exact-or-acb-singleton"
            ? "acb-exact-singleton-zero"
            : "exact-coefficient-field-only";
    json::object source{
        {"local", source_local_}, {"chart", source_chart_},
        {"source_operator_identity", source_operator_identity_},
        {"checkpoint_identity", source_checkpoint_},
        {"coefficient_domain", source_domain_}};
    if (planned_source_.has_value()) source = *planned_source_;
    json::object result{
        {"endpoint", handle_},
        {"capability", transport_bound()
             ? kRetainedTransportEndpointBatchCapability
             : rolling_transport_bound()
             ? "rolling-transport-endpoint-batch-v1"
             : direct_plan_bound()
             ? kRetainedPlannedEndpointLimitCapability
             : kRetainedEndpointLimitCapability},
        {"native_retained", true},
        {"retained_state", "specialized-acb-epsilon-vector"},
        {"json_coefficients", 0},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"execution_scope", transport_bound()
             ? "transport-state-final-vector-row-endpoint"
             : direct_plan_bound()
             ? "plan-bound-final-arm-endpoint"
             : "unplanned-low-level"},
        {"source", std::move(source)},
        {"dimension", result_.values.size()},
        {"epsilon_min", window.min_power},
        {"epsilon_max", window.complete_max},
        {"coefficient_field", "acb-specialized"},
        {"arithmetic_enclosed", true},
        {"approach_direction", approach_direction_},
        {"cancellation", json::object{
             {"mode", cancellation_mode_},
             {"effective_scope", cancellation_scope},
             {"numeric_singleton_cancellations",
              result_.cancelled_divergent_coefficients}}},
        {"analytic_regularization", json::object{
             {"regulator_slope_scope",
              "exact-zero-fact; certified-nonzero symbolic slopes allowed"},
             {"unregulated_power_scope", "exact-rational"},
             {"endpoint_rule", "drop-exact-nonzero-regulator-slope"},
              {"dropped_regulated_sectors",
               result_.dropped_regulated_sectors},
              {"metadata", analytic_metadata_}}},
        {"elapsed_ms", elapsed_ms_}};
    result["requested_rim"] = requested_rim_.has_value()
        ? json::value(*requested_rim_) : json::value(nullptr);
    if (derived_plan_bound()) {
      result["derived_rim"] = planned_effective_rim_.has_value()
          ? json::value(*planned_effective_rim_) : json::value(nullptr);
      result["effective_rim"] = planned_effective_rim_.has_value()
          ? json::value(*planned_effective_rim_) : json::value(nullptr);
      result["rim_source"] =
          "final-chart-exact-odd-multiplicity-prescriptions";
    } else {
      result["effective_rim"] = result_.imaginary_sign;
      result["rim_source"] = "unplanned-caller-or-principal-default";
    }
    return result;
  }

  json::object stats_json() const {
    auto out = summary();
    std::lock_guard<std::mutex> lock(stats_mutex_);
    out["exports"] = exports_;
    out["export_ms"] = export_ms_;
    return out;
  }

  json::object checkpoint_record() const {
    json::array values;
    values.reserve(result_.values.size());
    for (const auto& value : result_.values)
      values.push_back(checkpoint_epsilon_frame_record(value));
    std::lock_guard<std::mutex> lock(stats_mutex_);
    json::object source{
        {"local", source_local_}, {"chart", source_chart_},
        {"source_operator_identity", source_operator_identity_},
        {"checkpoint_identity", source_checkpoint_},
        {"coefficient_domain", source_domain_}};
    if (planned_source_.has_value()) source = *planned_source_;
    json::object record{
        {"schema", transport_bound()
             ? prepared_row_checkpoint_.has_value()
                 ? "diffexp3-retained-transport-endpoint-result-v2"
                 : "diffexp3-retained-transport-endpoint-result-v1"
             : rolling_transport_bound()
             ? prepared_row_checkpoint_.has_value()
                 ? "diffexp3-retained-rolling-transport-endpoint-result-v2"
                 : "diffexp3-retained-rolling-transport-endpoint-result-v1"
             : direct_plan_bound()
             ? "diffexp3-retained-plan-bound-endpoint-result-v1"
             : "diffexp3-retained-endpoint-result-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"source", std::move(source)},
        {"approach_direction", approach_direction_},
        {"cancellation_mode", cancellation_mode_},
        {"analytic_metadata", analytic_metadata_},
        {"result", json::object{
            {"values", std::move(values)},
            {"dropped_regulated_sectors",
             result_.dropped_regulated_sectors},
             {"cancelled_divergent_coefficients",
              result_.cancelled_divergent_coefficients},
            {"imaginary_sign", derived_plan_bound()
                 ? (planned_effective_rim_.has_value()
                       ? json::value(*planned_effective_rim_)
                       : json::value(nullptr))
                 : json::value(result_.imaginary_sign)}}},
        {"elapsed_ms", elapsed_ms_},
        {"runtime_stats", json::object{{"exports", exports_},
                                        {"export_ms", export_ms_}}}};
    if (derived_plan_bound()) {
      record["derived_rim"] = planned_effective_rim_.has_value()
          ? json::value(*planned_effective_rim_) : json::value(nullptr);
    } else {
      record["requested_rim"] = requested_rim_.has_value()
          ? json::value(*requested_rim_) : json::value(nullptr);
    }
    if (prepared_row_checkpoint_.has_value())
      record["prepared_row"] = *prepared_row_checkpoint_;
    return record;
  }

  void restore_runtime_stats(std::uint64_t exports, double export_ms) {
    std::lock_guard<std::mutex> lock(stats_mutex_);
    exports_ = exports;
    export_ms_ = export_ms;
  }

  json::object export_values(const std::string& expected_checkpoint,
                             int output_digits) {
    if (expected_checkpoint.empty())
      throw std::invalid_argument(
          "endpoint export checkpoint identity must be nonempty");
    if (expected_checkpoint != checkpoint_identity_)
      throw std::invalid_argument(
          "endpoint export checkpoint identity does not match retained state");
    const auto started = std::chrono::steady_clock::now();
    auto encoded = encode_epsilon_vector(
        endpoint_values_vector(result_), output_digits);
    const auto elapsed = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    {
      std::lock_guard<std::mutex> lock(stats_mutex_);
      ++exports_;
      export_ms_ += elapsed;
    }
    return json::object{
        {"endpoint", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"compatibility_export", true},
        {"coefficient_field", "acb-specialized"},
        {"json_coefficients", encoded.at("coefficients").as_array().size()},
        {"value", std::move(encoded)},
        {"elapsed_ms", elapsed}};
  }

 private:
  std::string handle_;
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::string source_local_;
  std::string source_chart_;
  std::string source_operator_identity_;
  std::string source_checkpoint_;
  std::string source_domain_;
  std::int32_t approach_direction_ = 1;
  std::optional<std::int32_t> requested_rim_;
  std::string cancellation_mode_;
  json::object analytic_metadata_;
  EndpointLimitResult result_;
  double elapsed_ms_ = 0.0;
  std::optional<json::object> planned_source_;
  std::optional<std::int32_t> planned_effective_rim_;
  std::shared_ptr<StoredTilePlan> plan_owner_;
  std::shared_ptr<StoredLocalBase> local_owner_;
  std::shared_ptr<StoredTransportArmState> transport_owner_;
  std::optional<json::object> prepared_row_checkpoint_;
  mutable std::mutex stats_mutex_;
  std::uint64_t exports_ = 0;
  double export_ms_ = 0.0;
};

std::shared_ptr<StoredEndpointResult> build_endpoint_limit(
    const std::string& endpoint_handle, const json::object& request,
    const std::shared_ptr<StoredLocalBase>& local) {
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  const auto expected_source_checkpoint = required_string(
      request, "source_checkpoint_identity");
  if (checkpoint_identity.empty() || expected_source_checkpoint.empty())
    throw std::invalid_argument(
        "endpoint checkpoint identities must be nonempty");
  if (expected_source_checkpoint != local->checkpoint_identity())
    throw std::invalid_argument(
        "endpoint source checkpoint identity does not match retained local");
  const auto policy = parse_endpoint_limit_policy(request);
  auto analytic_metadata = local->exact_analytic_metadata();
  json::object provenance{
      {"schema", "diffexp3-retained-native-endpoint-sector-limit-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"source", json::object{
           {"local", local->handle()}, {"chart", local->source_chart()},
           {"source_operator_identity",
            local->source_operator_identity()},
           {"checkpoint_identity", expected_source_checkpoint},
           {"coefficient_domain", local->scalar_domain()}}},
      {"approach_direction", policy.options.approach_direction},
      {"rim", policy.requested_rim.has_value()
           ? json::value(*policy.requested_rim) : json::value(nullptr)},
      {"cancellation", json::object{{"mode", policy.cancellation_mode}}},
      {"analytic_metadata", analytic_metadata}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto started = std::chrono::steady_clock::now();
  auto result = local->endpoint_limit(policy.options);
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredEndpointResult>(
      endpoint_handle, checkpoint_identity, provenance_identity,
      local->handle(), local->source_chart(),
      local->source_operator_identity(), expected_source_checkpoint,
      local->scalar_domain(), policy.options.approach_direction,
      policy.requested_rim, policy.cancellation_mode,
      std::move(analytic_metadata), std::move(result), elapsed);
}
