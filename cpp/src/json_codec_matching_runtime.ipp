json::object compact_matching_identity_reference(
    std::string_view identity) {
  if (identity.empty())
    throw std::invalid_argument(
        "compact matching identity reference cannot bind an empty identity");
  return json::object{
      {"algorithm", "fnv1a64-v1"},
      {"fingerprint", public_provenance_fingerprint(identity)},
      {"identity_bytes", identity.size()}};
}

bool compact_matching_identity_reference_matches(
    const json::value& raw, std::string_view expected_identity,
    const char* label) {
  const auto& reference = as_object(raw, label);
  require_exact_keys(
      reference, {"algorithm", "fingerprint", "identity_bytes"}, label);
  return reference ==
      compact_matching_identity_reference(expected_identity);
}

std::string_view matching_required_string_view(
    const json::object& object, const char* key) {
  const auto& value = object.at(key);
  if (!value.is_string() || value.as_string().empty())
    throw std::invalid_argument(
        std::string("compact matching ") + key +
        " must be a nonempty string");
  const auto& text = value.as_string();
  return std::string_view(text.data(), text.size());
}

// Exact source records are retained in StoredLocal and in checkpoint object
// records.  Matching identities need a bounded reference to those live
// authorities, not another copy of every SCC request and column schedule.
// The checkpoint/local handle, plan binding, and direct retained-owner checks
// remain authoritative; the FNV fields below are diagnostics, not a
// replacement proof.
json::object compact_matching_source_reference(
    const json::object& source) {
  json::object operator_reference;
  if (const auto* raw_reference =
          source.if_contains("source_operator_reference")) {
    operator_reference = as_object(
        *raw_reference, "compact matching source-operator reference");
    require_exact_keys(
        operator_reference,
        {"algorithm", "fingerprint", "identity_bytes"},
        "compact matching source-operator reference");
  } else {
    operator_reference = compact_matching_identity_reference(
        matching_required_string_view(
            source, "source_operator_identity"));
  }
  json::object compact{
      {"schema", "diffexp2-retained-match-source-reference-v1"},
      {"authority", "retained-native-local-binding-validated"},
      {"local", source.at("local")},
      {"chart", source.at("chart")},
      {"checkpoint_identity", source.at("checkpoint_identity")},
      {"source_operator_reference",
       std::move(operator_reference)}};
  if (const auto* column = source.if_contains("column"))
    compact["column"] = *column;
  if (const auto* sign =
          source.if_contains("requested_imaginary_sign"))
    compact["requested_imaginary_sign"] = *sign;
  if (const auto* sign =
          source.if_contains("effective_imaginary_sign"))
    compact["effective_imaginary_sign"] = *sign;
  if (const auto* width = source.if_contains("matching_taylor_width"))
    compact["matching_taylor_width"] = *width;
  if (const auto* raw_metadata = source.if_contains("analytic_metadata")) {
    const auto identity = json::serialize(
        canonical_json_value(*raw_metadata));
    compact["analytic_metadata_reference"] =
        compact_matching_identity_reference(identity);
  }
  if (const auto* raw_provenance =
          source.if_contains("column_provenance")) {
    const auto& provenance = as_object(
        *raw_provenance, "compact matching SCC-column provenance");
    require_exact_keys(
        provenance,
        {"scc", "scc_exact_identity", "seed_block", "basis_index",
         "exact_column_identity"},
        "compact matching SCC-column provenance");
    compact["column_provenance"] = json::object{
        {"schema", "diffexp2-retained-scc-column-reference-v1"},
        {"authority", "retained-native-exact-column-owner"},
        {"scc", provenance.at("scc")},
        {"seed_block", provenance.at("seed_block")},
        {"basis_index", provenance.at("basis_index")},
        {"identity_diagnostics",
         json::object{
             {"algorithm", "fnv1a64-v1"},
             {"scc_exact_identity_fingerprint",
              public_provenance_fingerprint(
                  matching_required_string_view(
                      provenance, "scc_exact_identity"))},
             {"scc_exact_identity_bytes",
              matching_required_string_view(
                  provenance, "scc_exact_identity").size()},
             {"exact_column_identity_fingerprint",
              public_provenance_fingerprint(
                  matching_required_string_view(
                      provenance, "exact_column_identity"))},
             {"exact_column_identity_bytes",
              matching_required_string_view(
                  provenance, "exact_column_identity").size()}}}};
  }
  return compact;
}

json::array compact_matching_source_references(
    const std::vector<json::object>& sources) {
  json::array compact;
  compact.reserve(sources.size());
  for (const auto& source : sources)
    compact.push_back(compact_matching_source_reference(source));
  return compact;
}

json::object checkpoint_matching_source_record(
    const json::object& source) {
  auto record = source;
  // This is a cached derivative of source_operator_identity, used only to
  // avoid rehashing multi-megabyte live identities.  Checkpoints and legacy
  // full-source proof records keep their established authoritative schema.
  record.erase("source_operator_reference");
  return record;
}

class StoredMatchBase {
 public:
  explicit StoredMatchBase(std::string handle) : handle_(std::move(handle)) {}
  virtual ~StoredMatchBase() = default;

  virtual json::object summary() const = 0;
  virtual json::object checkpoint_record() const = 0;
  const std::string& handle() const { return handle_; }

 protected:
  std::string handle_;
};

class StoredExactRegularMatch final : public StoredMatchBase {
 public:
  StoredExactRegularMatch(
      std::string handle, std::string checkpoint_identity,
      std::string provenance_identity,
      std::vector<json::object> basis_sources, json::object incoming_source,
      std::string basis_chart,
      std::string incoming_chart, std::string basis_point,
      std::string incoming_point, std::string physical_point,
      EpsilonWindow requested_window, std::int32_t required_complete_max,
      std::uint32_t dimension,
      ExactLaurentMatrix<Rational>&& transformation,
      FiniteLaurentVector<Rational>&& weights,
      EpsilonLatticeSaturationDiagnostics<Rational>&& diagnostics,
      EpsilonWindow residual_window, double elapsed_ms,
      std::vector<std::shared_ptr<StoredLocalBase>> basis_owners,
      std::shared_ptr<StoredLocalBase> incoming_owner)
      : StoredMatchBase(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        basis_sources_(std::move(basis_sources)),
        incoming_source_(std::move(incoming_source)),
        basis_chart_(std::move(basis_chart)),
        incoming_chart_(std::move(incoming_chart)),
        basis_point_(std::move(basis_point)),
        incoming_point_(std::move(incoming_point)),
        physical_point_(std::move(physical_point)),
        requested_window_(requested_window),
        required_complete_max_(required_complete_max),
        dimension_(dimension),
        transformation_(std::move(transformation)),
        weights_(std::move(weights)),
        diagnostics_(std::move(diagnostics)),
        residual_window_(residual_window),
        elapsed_ms_(elapsed_ms), basis_owners_(std::move(basis_owners)),
        incoming_owner_(std::move(incoming_owner)) {
    if (basis_sources_.size() != dimension_ ||
        basis_owners_.size() != dimension_ || !incoming_owner_)
      throw std::invalid_argument(
          "retained exact match source ownership differs from its dimension");
  }

  json::object summary() const override {
    json::array basis;
    basis.reserve(basis_sources_.size());
    for (const auto& source : basis_sources_)
      basis.push_back(json::object{
          {"column", source.at("column")},
          {"local", source.at("local")},
          {"checkpoint_identity", source.at("checkpoint_identity")}});

    json::array shifts;
    for (const auto shift : diagnostics_.initial_column_shifts)
      shifts.push_back(shift);
    json::array actions;
    for (const auto& action : diagnostics_.actions)
      actions.push_back(json::object{
          {"leading_rank_before", action.leading_rank_before},
          {"target_column", action.target_column},
          {"relation_support", std::count_if(
               action.null_relation.begin(), action.null_relation.end(),
               [](const Rational& value) { return !value.is_zero(); })}});

    json::array weight_windows;
    for (const auto& weight : weights_)
      weight_windows.push_back(json::object{{"min", weight.min_power()},
                                            {"max", weight.complete_max()}});

    std::size_t transformation_terms = 0;
    for (const auto& row : transformation_)
      for (const auto& entry : row)
        transformation_terms += entry.terms().size();

    return json::object{
        {"match", handle_},
        {"capability", kExactRegularLocalMatchCapability},
        {"retained_state",
         "exact-lattice-transformation-and-laurent-weights"},
        {"native_retained", true},
        {"json_coefficients", 0},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"dimension", dimension_},
        {"basis", std::move(basis)},
        {"incoming", incoming_source_.at("local")},
        {"incoming_checkpoint_identity",
         incoming_source_.at("checkpoint_identity")},
        {"basis_chart", basis_chart_},
        {"incoming_chart", incoming_chart_},
        {"basis_point_exact", basis_point_},
        {"incoming_point_exact", incoming_point_},
        {"physical_match_point_exact", physical_point_},
        {"epsilon", json::object{{"min", requested_window_.min_power},
                                  {"max", requested_window_.complete_max},
                                  {"required_complete_max",
                                   required_complete_max_}}},
        {"weight_windows", std::move(weight_windows)},
        {"transformation_terms", transformation_terms},
        {"initial_column_shifts", std::move(shifts)},
        {"normalized_determinant_valuation",
         diagnostics_.normalized_determinant_valuation},
        {"initial_leading_rank", diagnostics_.initial_leading_rank},
        {"final_leading_rank", diagnostics_.final_leading_rank},
        {"saturation_actions", std::move(actions)},
        {"residual", json::object{{"status", "exact-zero"},
                                   {"scope", "stored-taylor-truncation"},
                                   {"min", residual_window_.min_power},
                                   {"max", residual_window_.complete_max}}},
        {"elapsed_ms", elapsed_ms_}};
  }

  double elapsed_ms() const { return elapsed_ms_; }

  const std::vector<std::shared_ptr<StoredLocalBase>>& basis_owners() const {
    return basis_owners_;
  }
  const std::shared_ptr<StoredLocalBase>& incoming_owner() const {
    return incoming_owner_;
  }
  const FiniteLaurentVector<Rational>& weights() const { return weights_; }

  json::object checkpoint_record() const override {
    json::array basis;
    for (const auto& source : basis_sources_) {
      auto checkpoint_source = source;
      checkpoint_source.erase("source_operator_reference");
      basis.push_back(std::move(checkpoint_source));
    }
    auto checkpoint_incoming = incoming_source_;
    checkpoint_incoming.erase("source_operator_reference");
    return json::object{
        {"schema", "diffexp2-retained-exact-rational-match-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"basis_sources", std::move(basis)},
        {"incoming_source", std::move(checkpoint_incoming)},
        {"basis_chart", basis_chart_},
        {"incoming_chart", incoming_chart_},
        {"basis_point_exact", basis_point_},
        {"incoming_point_exact", incoming_point_},
        {"physical_match_point_exact", physical_point_},
        {"epsilon", json::object{{"min", requested_window_.min_power},
                                  {"max", requested_window_.complete_max},
                                  {"required_complete_max",
                                   required_complete_max_}}},
        {"dimension", dimension_},
        {"transformation",
         checkpoint_exact_laurent_matrix_record(transformation_)},
        {"weights", checkpoint_frame_vector_record(weights_)},
        {"saturation",
         checkpoint_saturation_diagnostics_record(diagnostics_)},
        {"residual_window",
         json::object{{"min", residual_window_.min_power},
                      {"max", residual_window_.complete_max}}},
        {"elapsed_ms", elapsed_ms_}};
  }

 private:
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::vector<json::object> basis_sources_;
  json::object incoming_source_;
  std::string basis_chart_;
  std::string incoming_chart_;
  std::string basis_point_;
  std::string incoming_point_;
  std::string physical_point_;
  EpsilonWindow requested_window_;
  std::int32_t required_complete_max_ = 0;
  std::uint32_t dimension_ = 0;
  ExactLaurentMatrix<Rational> transformation_;
  FiniteLaurentVector<Rational> weights_;
  EpsilonLatticeSaturationDiagnostics<Rational> diagnostics_;
  EpsilonWindow residual_window_;
  double elapsed_ms_ = 0.0;
  std::vector<std::shared_ptr<StoredLocalBase>> basis_owners_;
  std::shared_ptr<StoredLocalBase> incoming_owner_;
};

void require_exact_regular_local(const LocalSolution<Rational>& solution,
                                 EpsilonWindow requested_window,
                                 const RealEvaluationPoint& point,
                                 const std::string& label) {
  validate_local_solution(solution, false);
  if (!solution.error.empty())
    throw std::invalid_argument(
        label +
        " carries an error envelope; exact regular local matching does not "
        "silently discard certificates");
  if (solution.sectors.size() != 1)
    throw std::invalid_argument(
        label + " must contain exactly one regular local sector");
  const auto& sector = solution.sectors.front();
  if (sector.a.domain != ExactDomain::Rational ||
      sector.b.domain != ExactDomain::Rational ||
      !(Rational(sector.a.canonical) == Rational(0)) ||
      !(Rational(sector.b.canonical) == Rational(0)) ||
      sector.log_power != 0)
    throw std::invalid_argument(
        label +
        " is not an exact-rational regular (a=0,b=0,log=0) local");
  if (requested_window.complete_max > solution.epsilon.complete_max)
    throw std::invalid_argument(
        label + " does not cover the requested complete epsilon upper edge");
  if (!solution.chart.infinite_radius &&
      !arb_lt(acb_realref(point.modulus.raw()),
              acb_realref(solution.chart.radius.raw())))
    throw std::invalid_argument(
        label + " match point is not provably inside its chart radius");
}

bool same_chart_geometry(const ChartGeometry& left,
                         const ChartGeometry& right) {
  const bool exact_radius_available =
      !left.radius_exact.empty() && !right.radius_exact.empty();
  return left.center_exact == right.center_exact &&
         left.scale_exact == right.scale_exact &&
         left.infinite_radius == right.infinite_radius &&
         (left.infinite_radius ||
          (exact_radius_available
               ? left.radius_exact == right.radius_exact
               : acb_equal(left.radius.raw(), right.radius.raw())));
}

Rational physical_match_point(const ChartGeometry& chart,
                              const RealEvaluationPoint& local_point,
                              const std::string& label) {
  try {
    return Rational(chart.center_exact) +
           Rational(chart.scale_exact) *
               Rational(local_point.exact_coordinate);
  } catch (const std::invalid_argument&) {
    throw std::invalid_argument(
        label +
        " requires rational chart center/scale geometry; algebraic geometry "
        "is outside exact-regular-local-match-v1");
  }
}

FiniteLaurentVector<Rational> evaluate_exact_regular_local(
    const LocalSolution<Rational>& solution,
    const RealEvaluationPoint& point, EpsilonWindow window,
    const std::string& label) {
  const Rational t(point.exact_coordinate);
  std::vector<Rational> t_powers(solution.taylor_width(), Rational(1));
  for (std::size_t n = 1; n < t_powers.size(); ++n)
    t_powers[n] = t_powers[n - 1] * t;

  const auto& sector = solution.sectors.front();
  const auto coefficient_at = [&](std::int32_t power,
                                  std::uint32_t component) {
    if (power < solution.epsilon.min_power) return Rational(0);
    const auto epsilon_index = static_cast<std::size_t>(
        static_cast<std::int64_t>(power) - solution.epsilon.min_power);
    Rational coefficient(0);
    for (std::size_t n = 0; n < solution.taylor_width(); ++n)
      coefficient += sector.coefficients[local_detail::sector_index(
                         solution, epsilon_index, n, component)] *
                     t_powers[n];
    return coefficient;
  };
  for (std::int64_t raw_power = solution.epsilon.min_power;
       raw_power < window.min_power; ++raw_power) {
    const auto power = static_cast<std::int32_t>(raw_power);
    for (std::uint32_t component = 0; component < solution.dimension;
         ++component)
      if (!coefficient_at(power, component).is_zero())
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            label +
                " work minimum would discard a nonzero lower epsilon coefficient",
            component, std::nullopt, power);
  }

  FiniteLaurentVector<Rational> value;
  value.reserve(solution.dimension);
  for (std::uint32_t component = 0; component < solution.dimension;
       ++component) {
    std::vector<Rational> coefficients;
    coefficients.reserve(window.width());
    for (std::int64_t power = window.min_power;
         power <= window.complete_max; ++power) {
      coefficients.push_back(coefficient_at(
          static_cast<std::int32_t>(power), component));
    }
    value.emplace_back(window, std::move(coefficients));
  }
  return value;
}

std::shared_ptr<StoredExactRegularMatch> build_exact_regular_match(
    const std::string& match_handle, const json::object& request,
    const std::vector<std::string>& basis_handles,
    const std::vector<std::shared_ptr<StoredLocalBase>>& erased_basis,
    const std::string& incoming_handle,
    const std::shared_ptr<StoredLocalBase>& erased_incoming) {
  const auto started = std::chrono::steady_clock::now();
  const auto& raw_window = as_object(
      request.at("epsilon"), "exact regular match epsilon window");
  EpsilonWindow window{as_i32(raw_window.at("min"), "match epsilon minimum"),
                       as_i32(raw_window.at("max"), "match epsilon maximum")};
  (void)window.width();
  const auto required_complete_max = as_i32(
      raw_window.at("required_complete_max"),
      "required match residual complete maximum");
  if (required_complete_max < window.min_power ||
      required_complete_max > window.complete_max)
    throw std::invalid_argument(
        "required match residual maximum must lie inside the supplied work "
        "epsilon window");

  const auto basis_point = RealEvaluationPoint::rational(required_string(
      as_object(request.at("basis_point"), "basis match point"), "exact"));
  const auto incoming_point = RealEvaluationPoint::rational(required_string(
      as_object(request.at("incoming_point"), "incoming match point"),
      "exact"));
  const auto basis_chart = required_string(request, "basis_chart");
  const auto incoming_chart = required_string(request, "incoming_chart");
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  if (checkpoint_identity.empty())
    throw std::invalid_argument(
        "exact regular match checkpoint identity cannot be empty");

  std::vector<std::shared_ptr<StoredLocal<Rational>>> basis;
  basis.reserve(erased_basis.size());
  for (const auto& local : erased_basis) {
    auto typed = std::dynamic_pointer_cast<StoredLocal<Rational>>(local);
    if (!typed)
      throw std::invalid_argument(
          "exact regular local matching requires rational retained locals");
    basis.push_back(std::move(typed));
  }
  auto incoming =
      std::dynamic_pointer_cast<StoredLocal<Rational>>(erased_incoming);
  if (!incoming)
    throw std::invalid_argument(
        "exact regular local matching requires a rational incoming local");
  if (basis.empty())
    throw std::invalid_argument(
        "exact regular local matching requires a nonempty basis");
  const auto dimension = basis.front()->solution().dimension;
  if (basis.size() != dimension || incoming->solution().dimension != dimension)
    throw std::invalid_argument(
        "exact regular local matching requires d basis columns and a "
        "d-component incoming local");

  const auto& raw_basis_checkpoints = as_array(
      request.at("basis_checkpoint_identities"),
      "basis checkpoint identities");
  if (raw_basis_checkpoints.size() != basis.size())
    throw std::invalid_argument(
        "basis checkpoint identity count differs from the basis dimension");
  std::vector<std::string> basis_checkpoints;
  basis_checkpoints.reserve(basis.size());
  for (std::size_t column = 0; column < basis.size(); ++column) {
    if (!raw_basis_checkpoints[column].is_string())
      throw std::invalid_argument(
          "basis checkpoint identities must be strings");
    const std::string expected(raw_basis_checkpoints[column].as_string());
    if (basis[column]->solution().checkpoint_identity != expected)
      throw std::invalid_argument(
          "basis checkpoint provenance mismatch at column " +
          std::to_string(column));
    basis_checkpoints.push_back(expected);
    if (basis[column]->source_chart() != basis_chart)
      throw std::invalid_argument(
          "basis chart provenance mismatch at column " +
          std::to_string(column));
    if (!same_chart_geometry(basis.front()->solution().chart,
                             basis[column]->solution().chart))
      throw std::invalid_argument(
          "basis locals do not share identical retained chart geometry");
    require_exact_regular_local(
        basis[column]->solution(), window, basis_point,
        "basis local " + basis_handles[column]);
  }
  const auto expected_incoming_checkpoint = required_string(
      request, "incoming_checkpoint_identity");
  if (incoming->solution().checkpoint_identity !=
      expected_incoming_checkpoint)
    throw std::invalid_argument("incoming checkpoint provenance mismatch");
  if (incoming->source_chart() != incoming_chart)
    throw std::invalid_argument("incoming chart provenance mismatch");
  require_exact_regular_local(incoming->solution(), window, incoming_point,
                              "incoming local " + incoming_handle);
  const auto basis_physical_point = physical_match_point(
      basis.front()->solution().chart, basis_point, "basis match point");
  for (std::size_t column = 1; column < basis.size(); ++column)
    if (!(physical_match_point(
              basis[column]->solution().chart, basis_point,
              "basis match point at column " + std::to_string(column)) ==
          basis_physical_point))
      throw std::invalid_argument(
          "basis locals do not name one exact physical match point");
  const auto incoming_physical_point = physical_match_point(
      incoming->solution().chart, incoming_point, "incoming match point");
  if (!(basis_physical_point == incoming_physical_point))
    throw std::invalid_argument(
        "basis and incoming local coordinates do not name the same exact "
        "physical match point");

  FiniteLaurentMatrix<Rational> evaluated_basis(
      dimension, FiniteLaurentVector<Rational>());
  for (auto& row : evaluated_basis) row.reserve(dimension);
  for (std::size_t column = 0; column < basis.size(); ++column) {
    auto value = evaluate_exact_regular_local(
        basis[column]->solution(), basis_point, window,
        "basis local " + basis_handles[column]);
    for (std::uint32_t component = 0; component < dimension; ++component)
      evaluated_basis[component].push_back(std::move(value[component]));
  }
  auto incoming_value = evaluate_exact_regular_local(
      incoming->solution(), incoming_point, window,
      "incoming local " + incoming_handle);

  // Exact zero rows at the lower edge are certified structural zeros.  Trim
  // them once before both saturation and the final residual so honest
  // completeness is not lost merely because a caller supplied a wider
  // declared Laurent minimum than the evaluated entry actually needs.
  for (std::uint32_t row = 0; row < dimension; ++row) {
    incoming_value[row] = matching_detail::canonical_leading_frame(
        incoming_value[row], checkpoint_identity + ":incoming", row);
    for (std::uint32_t column = 0; column < dimension; ++column)
      evaluated_basis[row][column] =
          matching_detail::canonical_leading_frame(
              evaluated_basis[row][column],
              checkpoint_identity + ":basis", row, column);
  }

  auto saturated = saturate_finite_laurent_basis(
      evaluated_basis, checkpoint_identity + ":saturation");
  auto saturated_weights = solve_finite_laurent_system(
      saturated.basis_times_transformation, incoming_value,
      checkpoint_identity + ":solve");
  auto weights = apply_exact_laurent_matrix(
      saturated.transformation, saturated_weights);
  auto reconstructed = apply_finite_laurent_matrix(
      evaluated_basis, weights);

  std::int32_t residual_min = reconstructed.front().min_power();
  std::int32_t residual_max = reconstructed.front().complete_max();
  for (std::uint32_t component = 0; component < dimension; ++component) {
    auto residual = reconstructed[component] - incoming_value[component];
    if (const auto leading = finite_laurent_leading_power(
            residual, checkpoint_identity + ":residual");
        leading.has_value())
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::SaturationFailure,
          checkpoint_identity +
              ": exact regular match residual is nonzero in its complete "
              "window",
          component, std::nullopt, *leading);
    residual_min = std::min(residual_min, residual.min_power());
    residual_max = std::min(residual_max, residual.complete_max());
  }
  (void)EpsilonWindow{residual_min, residual_max}.width();
  if (residual_max < required_complete_max)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InsufficientCompleteWindow,
        checkpoint_identity +
            ": exact regular matching consumed the required complete "
            "epsilon window",
        std::nullopt, std::nullopt, residual_max);

  std::vector<json::object> basis_sources;
  basis_sources.reserve(basis.size());
  for (std::size_t column = 0; column < basis.size(); ++column) {
    json::object entry{{"column", column},
                       {"local", basis_handles[column]},
                       {"chart", basis_chart},
                       {"source_operator_identity",
                        basis[column]->source_operator_identity()},
                       {"source_operator_reference",
                        basis[column]->source_operator_reference()},
                       {"checkpoint_identity", basis_checkpoints[column]}};
    entry["analytic_metadata"] =
        basis[column]->exact_analytic_metadata();
    if (basis[column]->column_provenance().has_value())
      entry["column_provenance"] =
          basis[column]->column_provenance()->encode();
    basis_sources.push_back(std::move(entry));
  }
  json::object incoming_source{
      {"local", incoming_handle}, {"chart", incoming_chart},
      {"source_operator_identity", incoming->source_operator_identity()},
      {"source_operator_reference",
       incoming->source_operator_reference()},
      {"checkpoint_identity", expected_incoming_checkpoint},
      {"analytic_metadata", incoming->exact_analytic_metadata()}};
  if (incoming->column_provenance().has_value())
    incoming_source["column_provenance"] =
        incoming->column_provenance()->encode();
  auto provenance_basis =
      compact_matching_source_references(basis_sources);
  json::object provenance{
      {"schema", "diffexp2-native-exact-regular-local-match-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"basis", std::move(provenance_basis)},
      {"incoming",
       compact_matching_source_reference(incoming_source)},
      {"basis_point_exact", basis_point.exact_coordinate},
      {"incoming_point_exact", incoming_point.exact_coordinate},
      {"physical_match_point_exact", basis_physical_point.str()},
      {"epsilon", json::object{{"min", window.min_power},
                                {"max", window.complete_max},
                                {"required_complete_max",
                                 required_complete_max}}}};
  const auto provenance_identity =
      json::serialize(canonical_json_value(provenance));
  const auto elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredExactRegularMatch>(
      match_handle, checkpoint_identity, provenance_identity,
      std::move(basis_sources), std::move(incoming_source), basis_chart,
      incoming_chart,
      basis_point.exact_coordinate, incoming_point.exact_coordinate,
      basis_physical_point.str(), window, required_complete_max, dimension,
      std::move(saturated.transformation), std::move(weights),
      std::move(saturated.diagnostics),
      EpsilonWindow{residual_min, residual_max}, elapsed_ms, erased_basis,
      erased_incoming);
}

const char* acb_match_verdict_name(AcbMatchingResidualVerdict verdict) {
  if (verdict == AcbMatchingResidualVerdict::Pass) return "pass";
  if (verdict == AcbMatchingResidualVerdict::Fail) return "fail";
  return "inconclusive";
}

json::object encode_acb_match_residual_diagnostics(
    const AcbMatchingResidualDiagnostics& diagnostics) {
  std::size_t pass = 0, fail = 0, inconclusive = 0;
  std::size_t required_pass = 0, required_fail = 0,
              required_inconclusive = 0;
  std::size_t required_zero_overlaps = 0;
  double maximum_required_residual_to_scale_upper = 0.0;
  json::array inconclusive_examples;
  json::array required_inconclusive_examples;
  for (const auto& coefficient : diagnostics.coefficients) {
    const bool required =
        coefficient.epsilon_power <= diagnostics.required_complete_max;
    if (required) {
      if (coefficient.residual_lower.is_zero())
        ++required_zero_overlaps;
      maximum_required_residual_to_scale_upper = std::max(
          maximum_required_residual_to_scale_upper,
          (coefficient.residual_upper / coefficient.scale_lower)
              .approximate_upper());
    }
    if (coefficient.verdict == AcbMatchingResidualVerdict::Pass) {
      ++pass;
      if (required) ++required_pass;
    } else if (coefficient.verdict == AcbMatchingResidualVerdict::Fail) {
      ++fail;
      if (required) ++required_fail;
    } else {
      ++inconclusive;
      const auto example = [&]() {
        return json::object{
            {"row", coefficient.row},
            {"epsilon_power", coefficient.epsilon_power},
            {"residual_upper_exact",
             coefficient.residual_upper.dump_exact()},
            {"scale_lower_exact", coefficient.scale_lower.dump_exact()}};
      };
      if (inconclusive_examples.size() < 6)
        inconclusive_examples.push_back(example());
      if (required) {
        ++required_inconclusive;
        if (required_inconclusive_examples.size() < 6)
          required_inconclusive_examples.push_back(example());
      }
    }
  }
  return json::object{
      {"verdict", acb_match_verdict_name(diagnostics.verdict)},
      {"complete_window",
       json::object{{"min", diagnostics.complete_window.min_power},
                    {"max", diagnostics.complete_window.complete_max}}},
      {"required_complete_max", diagnostics.required_complete_max},
      {"complete_through_required",
       diagnostics.complete_through_required},
      {"coefficient_diagnostics", diagnostics.coefficients.size()},
      {"coefficient_verdicts",
       json::object{{"pass", pass}, {"fail", fail},
                    {"inconclusive", inconclusive}}},
      {"required_coefficient_verdicts",
       json::object{{"pass", required_pass},
                    {"fail", required_fail},
                    {"inconclusive", required_inconclusive}}},
      {"required_zero_overlapping_coefficients",
       required_zero_overlaps},
      {"maximum_required_residual_to_scale_upper_approx",
       maximum_required_residual_to_scale_upper},
      {"inconclusive_examples", std::move(inconclusive_examples)},
      {"required_inconclusive_examples",
       std::move(required_inconclusive_examples)},
      {"detail", diagnostics.detail}};
}

json::object diagnose_acb_laurent_matrix_overlap(
    const FiniteLaurentMatrix<ComplexBall>& left,
    const FiniteLaurentMatrix<ComplexBall>& right,
    const std::string& label) {
  if (left.size() != right.size())
    throw std::invalid_argument(
        label + ": compared Laurent matrices have different row counts");
  std::size_t compared = 0;
  std::size_t overlapping = 0;
  std::size_t disjoint = 0;
  double maximum_relative_discrepancy_upper = 0.0;
  json::array disjoint_examples;
  for (std::size_t row = 0; row < left.size(); ++row) {
    if (left[row].size() != right[row].size())
      throw std::invalid_argument(
          label + ": compared Laurent matrices have different column counts");
    for (std::size_t column = 0; column < left[row].size(); ++column) {
      const auto minimum = std::max(
          left[row][column].min_power(),
          right[row][column].min_power());
      const auto complete_max = std::min(
          left[row][column].complete_max(),
          right[row][column].complete_max());
      for (std::int64_t raw_power = minimum;
           raw_power <= complete_max; ++raw_power) {
        const auto power = matching_detail::checked_power(
            raw_power, "diagnostic Laurent-matrix overlap power");
        const auto& left_value =
            left[row][column].coefficient(power);
        const auto& right_value =
            right[row][column].coefficient(power);
        const auto discrepancy = left_value - right_value;
        const bool overlaps = discrepancy.contains_zero();
        ++compared;
        if (overlaps) {
          ++overlapping;
        } else {
          ++disjoint;
          if (disjoint_examples.size() < 8)
            disjoint_examples.push_back(json::object{
                {"row", row},
                {"column", column},
                {"epsilon_power", power},
                {"discrepancy_absolute_upper_approx",
                 Magnitude::upper_abs(discrepancy)
                     .approximate_upper()},
                {"left_absolute_upper_approx",
                 Magnitude::upper_abs(left_value)
                     .approximate_upper()},
                {"right_absolute_upper_approx",
                 Magnitude::upper_abs(right_value)
                     .approximate_upper()},
                {"discrepancy_real_radius_exponent",
                 discrepancy.real_radius_exponent()},
                {"discrepancy_imag_radius_exponent",
                 discrepancy.imag_radius_exponent()}});
        }
        const auto scale = std::max(
            Magnitude::upper_abs(left_value).approximate_upper(),
            Magnitude::upper_abs(right_value).approximate_upper());
        const auto discrepancy_upper =
            Magnitude::upper_abs(discrepancy).approximate_upper();
        const auto relative =
            scale > 0.0 ? discrepancy_upper / scale
                        : (discrepancy_upper == 0.0 ? 0.0
                                                    : std::numeric_limits<double>::infinity());
        maximum_relative_discrepancy_upper = std::max(
            maximum_relative_discrepancy_upper, relative);
      }
    }
  }
  return json::object{
      {"schema", "diffexp2-acb-laurent-matrix-overlap-diagnostic-v1"},
      {"label", label},
      {"compared_coefficients", compared},
      {"overlapping_coefficients", overlapping},
      {"disjoint_coefficients", disjoint},
      {"all_coefficients_overlap", disjoint == 0},
      {"maximum_relative_discrepancy_upper_approx",
       maximum_relative_discrepancy_upper},
      {"disjoint_examples", std::move(disjoint_examples)}};
}

json::object checkpoint_acb_match_residual_record(
    const AcbMatchingResidualDiagnostics& diagnostics) {
  json::array coefficients;
  coefficients.reserve(diagnostics.coefficients.size());
  for (const auto& coefficient : diagnostics.coefficients)
    coefficients.push_back(json::object{
        {"row", coefficient.row},
        {"epsilon_power", coefficient.epsilon_power},
        {"residual_lower_exact", coefficient.residual_lower.dump_exact()},
        {"residual_upper_exact", coefficient.residual_upper.dump_exact()},
        {"scale_lower_exact", coefficient.scale_lower.dump_exact()},
        {"scale_upper_exact", coefficient.scale_upper.dump_exact()},
        {"verdict", acb_match_verdict_name(coefficient.verdict)}});
  return json::object{
      {"verdict", acb_match_verdict_name(diagnostics.verdict)},
      {"complete_window",
       json::object{{"min", diagnostics.complete_window.min_power},
                    {"max", diagnostics.complete_window.complete_max}}},
      {"required_complete_max", diagnostics.required_complete_max},
      {"complete_through_required",
       diagnostics.complete_through_required},
      {"coefficients", std::move(coefficients)},
      {"detail", diagnostics.detail}};
}

json::array checkpoint_acb_laurent_matrix_record(
    const ExactLaurentMatrix<ComplexBall>& matrix);

LocalSolution<ComplexBall> specialize_rational_local_solution_to_acb(
    const LocalSolution<Rational>& exact,
    const std::string& checkpoint_identity) {
  LocalSolution<ComplexBall> result;
  result.chart = exact.chart;
  result.epsilon = exact.epsilon;
  result.taylor_complete_max = exact.taylor_complete_max;
  result.dimension = exact.dimension;
  result.prescriptions = exact.prescriptions;
  result.error = exact.error;
  result.checkpoint_identity = checkpoint_identity;
  result.sectors.reserve(exact.sectors.size());
  for (const auto& exact_sector : exact.sectors) {
    LocalSector<ComplexBall> sector;
    sector.a = exact_sector.a;
    sector.b = exact_sector.b;
    sector.log_power = exact_sector.log_power;
    sector.coefficients.reserve(exact_sector.coefficients.size());
    for (const auto& coefficient : exact_sector.coefficients)
      sector.coefficients.push_back(
          ComplexBall::from_strings(coefficient.str()));
    result.sectors.push_back(std::move(sector));
  }
  return result;
}

class StoredRefinedAcbMatch final : public StoredMatchBase {
 public:
  StoredRefinedAcbMatch(
      std::string handle, std::string checkpoint_identity,
      std::string provenance_identity, std::string exact_lattice_identity,
      std::string exact_lattice_provenance_identity,
      std::string exact_lattice_witness_record,
      std::string saturation_witness_schema,
      std::vector<json::object> basis_sources, json::object incoming_source,
      std::string basis_chart, std::string incoming_chart,
      std::string basis_point, std::string incoming_point,
      std::string physical_point, std::string matching_frame_identity,
      std::string residual_frame_identity,
      std::string residual_certificate_identity,
      EpsilonWindow requested_window,
      std::int32_t required_complete_max, std::uint32_t dimension,
      std::string relative_tolerance, std::size_t max_refinement_steps,
      EpsilonLatticeSaturationResult<Rational>&& exact_saturation,
      std::optional<ExactLaurentMatrix<Rational>>
          exact_right_materialization_transformation,
      std::optional<ExactLaurentMatrix<ComplexBall>>
          acb_right_materialization_transformation,
      std::optional<ExactLaurentMatrix<ComplexBall>>
          acb_right_materialization_preconditioner,
      std::optional<FiniteLaurentMatrix<ComplexBall>>
          terminal_normal_frame_right_transformation,
      std::optional<std::vector<LocalSolution<ComplexBall>>>
          exact_shadow_factorized_basis,
      std::size_t exact_shadow_extra_precision_bits,
      std::shared_ptr<const PhysicalEquationOwnerBase>
          exact_shadow_equation_owner,
      RefinedAcbLaurentMatch&& refined, double elapsed_ms)
      : StoredMatchBase(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        exact_lattice_identity_(std::move(exact_lattice_identity)),
        exact_lattice_provenance_identity_(
            std::move(exact_lattice_provenance_identity)),
        exact_lattice_witness_record_(
            std::move(exact_lattice_witness_record)),
        saturation_witness_schema_(
            std::move(saturation_witness_schema)),
        basis_sources_(std::move(basis_sources)),
        incoming_source_(std::move(incoming_source)),
        basis_chart_(std::move(basis_chart)),
        incoming_chart_(std::move(incoming_chart)),
        basis_point_(std::move(basis_point)),
        incoming_point_(std::move(incoming_point)),
        physical_point_(std::move(physical_point)),
        matching_frame_identity_(std::move(matching_frame_identity)),
        residual_frame_identity_(std::move(residual_frame_identity)),
        residual_certificate_identity_(
            std::move(residual_certificate_identity)),
        requested_window_(requested_window),
        required_complete_max_(required_complete_max),
        dimension_(dimension),
        relative_tolerance_(std::move(relative_tolerance)),
        max_refinement_steps_(max_refinement_steps),
        exact_saturation_(std::move(exact_saturation)),
        exact_right_materialization_transformation_(
            std::move(exact_right_materialization_transformation)),
        acb_right_materialization_preconditioner_(
            std::move(acb_right_materialization_preconditioner)),
        terminal_normal_frame_right_transformation_(
            std::move(terminal_normal_frame_right_transformation)),
        exact_shadow_factorized_basis_(
            std::move(exact_shadow_factorized_basis)),
        exact_shadow_extra_precision_bits_(
            exact_shadow_extra_precision_bits),
        exact_shadow_equation_owner_(
            std::move(exact_shadow_equation_owner)),
        refined_(std::move(refined)), elapsed_ms_(elapsed_ms) {
    if (exact_right_materialization_transformation_.has_value()) {
      if (acb_right_materialization_transformation.has_value())
        throw std::logic_error(
            "retained Acb match received both exact and Acb right materialization transformations");
      acb_materialization_right_transformation_ =
          specialize_exact_rational_laurent_matrix_to_acb(
              *exact_right_materialization_transformation_);
    } else {
      acb_materialization_right_transformation_ =
          std::move(acb_right_materialization_transformation);
    }
    const bool normalized_matching_frame =
        matching_frame_identity_ != "physical-parent-frame";
    if (normalized_matching_frame &&
        !terminal_normal_frame_right_transformation_.has_value() &&
        acb_materialization_right_transformation_.has_value())
      throw std::logic_error(
          "nonterminal normalized Acb match exported an internal right "
          "coordinate transformation for physical materialization");
    if (!normalized_matching_frame &&
        terminal_normal_frame_right_transformation_.has_value())
      throw std::logic_error(
          "physical Acb match retained a terminal normal-frame right "
          "transformation");
    if (exact_shadow_factorized_basis_.has_value()) {
      if (!terminal_normal_frame_right_transformation_.has_value() ||
          exact_shadow_factorized_basis_->size() != dimension_ ||
          exact_shadow_extra_precision_bits_ == 0)
        throw std::logic_error(
            "retained exact-shadow factorized basis is not bound to one "
            "terminal normal-frame match");
      for (const auto& column : *exact_shadow_factorized_basis_)
        if (column.dimension != dimension_ ||
            column.sectors.empty())
          throw std::logic_error(
              "retained exact-shadow factorized column has invalid shape");
    } else if (exact_shadow_extra_precision_bits_ != 0) {
      throw std::logic_error(
          "retained match has exact-shadow precision without its basis");
    }
  }

  json::object summary() const override {
    auto basis = compact_matching_source_references(
        basis_sources_);

    json::array history;
    history.reserve(refined_.residual_history.size());
    for (std::size_t iteration = 0;
         iteration < refined_.residual_history.size(); ++iteration) {
      auto encoded = encode_acb_match_residual_diagnostics(
          refined_.residual_history[iteration]);
      encoded["iteration"] = iteration;
      const auto final = iteration + 1 == refined_.residual_history.size();
      encoded["frame_identity"] = final
          ? residual_frame_identity_ : matching_frame_identity_;
      encoded["publication_authority"] = final;
      history.push_back(std::move(encoded));
    }
    if (refined_.residual_history.empty())
      throw std::logic_error("retained Acb match has no residual history");

    json::array weight_windows;
    for (const auto& weight : refined_.weights)
      weight_windows.push_back(json::object{{"min", weight.min_power()},
                                            {"max", weight.complete_max()}});
    json::array transformed_weight_windows;
    for (const auto& weight : refined_.transformed_weights)
      transformed_weight_windows.push_back(
          json::object{{"min", weight.min_power()},
                       {"max", weight.complete_max()}});

    json::array shifts;
    for (const auto shift :
         exact_saturation_.diagnostics.initial_column_shifts)
      shifts.push_back(shift);
    json::array actions;
    for (const auto& action : exact_saturation_.diagnostics.actions)
      actions.push_back(json::object{
          {"leading_rank_before", action.leading_rank_before},
          {"target_column", action.target_column},
          {"relation_support", std::count_if(
               action.null_relation.begin(), action.null_relation.end(),
               [](const Rational& value) { return !value.is_zero(); })}});
    std::size_t transformation_terms = 0;
    std::optional<std::int32_t> transformation_minimum;
    for (const auto& row : exact_saturation_.transformation)
      for (const auto& entry : row) {
        transformation_terms += entry.terms().size();
        if (const auto minimum = entry.minimum_power(); minimum.has_value())
          transformation_minimum = !transformation_minimum.has_value()
              ? *minimum
              : std::min(*transformation_minimum, *minimum);
      }

    auto residual = encode_acb_match_residual_diagnostics(
        refined_.residual_history.back());
    residual["scope"] = "stored-taylor-truncation";
    residual["frame_identity"] = residual_frame_identity_;
    residual["history"] = std::move(history);
    json::object record{
        {"match", handle_},
        {"capability", kRefinedAcbLocalMatchCapability},
        {"retained_state",
         "exact-lattice-transformation-acb-weights-and-residual"},
        {"native_retained", true},
        {"json_coefficients", 0},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"dimension", dimension_},
        {"basis", std::move(basis)},
        {"incoming",
         compact_matching_source_reference(incoming_source_)},
        {"basis_chart", basis_chart_},
        {"incoming_chart", incoming_chart_},
        {"basis_point_exact", basis_point_},
        {"incoming_point_exact", incoming_point_},
        {"physical_match_point_exact", physical_point_},
        {"matching_frame_identity", matching_frame_identity_},
        {"residual_frame_identity", residual_frame_identity_},
        {"epsilon",
         json::object{{"min", requested_window_.min_power},
                      {"max", requested_window_.complete_max},
                      {"required_complete_max", required_complete_max_}}},
        {"exact_lattice",
         json::object{
             {"schema", saturation_witness_schema_},
             {"identity", exact_lattice_identity_},
             {"provenance_identity",
              exact_lattice_provenance_identity_},
             {"canonical_witness_retained", true},
             {"canonical_witness_bytes",
              exact_lattice_witness_record_.size()},
             {"transformation_terms", transformation_terms},
             {"transformation_min_power",
              transformation_minimum.has_value()
                  ? json::value(*transformation_minimum)
                  : json::value(nullptr)},
             {"initial_column_shifts", std::move(shifts)},
             {"normalized_determinant_valuation",
              exact_saturation_.diagnostics
                  .normalized_determinant_valuation},
             {"initial_leading_rank",
              exact_saturation_.diagnostics.initial_leading_rank},
             {"final_leading_rank",
              exact_saturation_.diagnostics.final_leading_rank},
             {"saturation_actions", std::move(actions)}}},
        {"refinement",
         json::object{{"relative_tolerance", relative_tolerance_},
                      {"max_steps", max_refinement_steps_},
                      {"steps", refined_.refinement_steps},
                      {"factorization_preconditioner",
                       refined_.factorization_preconditioner},
                      {"factorized_authoritative_rhs",
                       refined_.factorized_authoritative_rhs},
                      {"factorized_authoritative_rhs_columns",
                       refined_
                           .factorized_authoritative_rhs_columns},
                      {"factorizations", 1}}},
        {"weight_windows", std::move(weight_windows)},
        {"transformed_weight_windows",
         std::move(transformed_weight_windows)},
        {"normal_frame_attempt", normal_frame_attempt_},
        {"materialization_association",
         exact_shadow_factorized_basis_.has_value()
             ? "terminal-exact-shadow-(F*T)*P*y"
         : terminal_normal_frame_right_transformation_.has_value()
             ? "terminal-normal-frame-F*(V*T)*y"
         : acb_materialization_right_transformation_.has_value()
             ? exact_right_materialization_transformation_.has_value()
             ? acb_right_materialization_preconditioner_.has_value()
                 ? "conditioned-exact-right-((F*T)*P)*y"
                 : "exact-right-(F*T)*y"
             : "certified-acb-right-(F*T_acb)*y"
             : "physical-F*w"},
        {"residual", std::move(residual)},
        {"elapsed_ms", elapsed_ms_}};
    if (!residual_certificate_identity_.empty())
      record["residual_certificate_identity"] =
          residual_certificate_identity_;
    return record;
  }

  json::object compact_terminal_diagnostic_summary() const {
    if (refined_.residual_history.empty())
      throw std::logic_error(
          "retained Acb terminal match has no residual history");
    const auto taylor_complete_max = matching_taylor_complete_max();
    return json::object{
        {"schema", "diffexp2-compact-terminal-match-diagnostic-v1"},
        {"basis_point_exact", basis_point_},
        {"incoming_point_exact", incoming_point_},
        {"physical_match_point_exact", physical_point_},
        {"matching_taylor_width",
         taylor_complete_max.has_value()
             ? json::value(
                   static_cast<std::uint64_t>(*taylor_complete_max) + 1)
             : json::value(nullptr)},
        {"matching_frame_reference",
         json::object{
             {"algorithm", "fnv1a64-v1"},
             {"fingerprint",
              public_provenance_fingerprint(
                  matching_frame_identity_)},
             {"identity_bytes", matching_frame_identity_.size()}}},
        {"residual_frame_reference",
         json::object{
             {"algorithm", "fnv1a64-v1"},
             {"fingerprint",
              public_provenance_fingerprint(
                  residual_frame_identity_)},
             {"identity_bytes", residual_frame_identity_.size()}}},
        {"epsilon",
         json::object{{"min", requested_window_.min_power},
                      {"max", requested_window_.complete_max},
                      {"required_complete_max", required_complete_max_}}},
        {"exact_lattice",
         json::object{
             {"schema", saturation_witness_schema_},
             {"normalized_determinant_valuation",
              exact_saturation_.diagnostics
                  .normalized_determinant_valuation},
             {"initial_leading_rank",
              exact_saturation_.diagnostics.initial_leading_rank},
             {"final_leading_rank",
              exact_saturation_.diagnostics.final_leading_rank}}},
        {"refinement",
         json::object{{"relative_tolerance", relative_tolerance_},
                      {"max_steps", max_refinement_steps_},
                      {"steps", refined_.refinement_steps}}},
        {"normal_frame_attempt", normal_frame_attempt_},
        {"materialization_association",
         exact_shadow_factorized_basis_.has_value()
             ? "terminal-exact-shadow-(F*T)*P*y"
         : terminal_normal_frame_right_transformation_.has_value()
             ? "terminal-normal-frame-F*(V*T)*y"
         : acb_materialization_right_transformation_.has_value()
             ? exact_right_materialization_transformation_.has_value()
             ? acb_right_materialization_preconditioner_.has_value()
                 ? "conditioned-exact-right-((F*T)*P)*y"
                 : "exact-right-(F*T)*y"
             : "certified-acb-right-(F*T_acb)*y"
             : "physical-F*w"},
        {"residual", encode_acb_match_residual_diagnostics(
             refined_.residual_history.back())},
        {"elapsed_ms", elapsed_ms_}};
  }

  const FiniteLaurentVector<ComplexBall>& weights() const {
    return refined_.weights;
  }

  const FiniteLaurentVector<ComplexBall>& transformed_weights() const {
    return refined_.transformed_weights;
  }

  json::object retained_singular_saturation_request() const {
    const auto proof = json::parse(exact_lattice_witness_record_);
    const auto& object = as_object(
        proof, "retained Acb singular saturation proof");
    const auto* request = object.if_contains("native_request");
    if (request == nullptr)
      throw std::invalid_argument(
          "retained Acb match has no singular-SCC saturation request");
    return as_object(
        *request, "retained Acb singular saturation request");
  }

  const std::optional<ExactLaurentMatrix<Rational>>&
  exact_right_materialization_transformation() const {
    return exact_right_materialization_transformation_;
  }

  const std::optional<ExactLaurentMatrix<ComplexBall>>&
  acb_materialization_right_transformation() const {
    return acb_materialization_right_transformation_;
  }

  const std::optional<ExactLaurentMatrix<ComplexBall>>&
  acb_right_materialization_preconditioner() const {
    return acb_right_materialization_preconditioner_;
  }

  const std::optional<FiniteLaurentMatrix<ComplexBall>>&
  terminal_normal_frame_right_transformation() const {
    return terminal_normal_frame_right_transformation_;
  }

  const std::optional<std::vector<LocalSolution<ComplexBall>>>&
  exact_shadow_factorized_basis() const {
    return exact_shadow_factorized_basis_;
  }

  std::size_t exact_shadow_extra_precision_bits() const {
    return exact_shadow_extra_precision_bits_;
  }

  const std::shared_ptr<const PhysicalEquationOwnerBase>&
  exact_shadow_equation_owner() const {
    return exact_shadow_equation_owner_;
  }

  const ExactLaurentMatrix<Rational>&
  exact_saturation_transformation() const {
    return exact_saturation_.transformation;
  }

  const std::string& basis_point_exact() const {
    return basis_point_;
  }

  const std::string& incoming_point_exact() const {
    return incoming_point_;
  }

  const EpsilonWindow& requested_window() const {
    return requested_window_;
  }

  const std::string& relative_tolerance_text() const {
    return relative_tolerance_;
  }

  std::optional<std::int32_t> effective_basis_imaginary_sign() const {
    if (basis_sources_.empty())
      throw std::logic_error(
          "retained Acb match has no basis-source branch record");
    const auto& raw =
        basis_sources_.front().at("effective_imaginary_sign");
    if (raw.is_null()) return std::nullopt;
    const auto sign = as_i32(
        raw, "retained Acb basis effective imaginary sign");
    if (sign != -1 && sign != 1)
      throw std::logic_error(
          "retained Acb basis effective imaginary sign is invalid");
    return sign;
  }

  std::optional<std::int32_t> effective_incoming_imaginary_sign() const {
    const auto& raw =
        incoming_source_.at("effective_imaginary_sign");
    if (raw.is_null()) return std::nullopt;
    const auto sign = as_i32(
        raw, "retained Acb incoming effective imaginary sign");
    if (sign != -1 && sign != 1)
      throw std::logic_error(
          "retained Acb incoming effective imaginary sign is invalid");
    return sign;
  }

  std::optional<std::uint32_t> matching_taylor_complete_max() const {
    std::optional<std::uint64_t> common_width;
    for (const auto& source : basis_sources_) {
      const auto* raw = source.if_contains("matching_taylor_width");
      if (raw == nullptr) return std::nullopt;
      const auto width =
          as_u64(*raw, "retained Acb matching Taylor width");
      if (width == 0 || (common_width.has_value() &&
                         *common_width != width))
        throw std::logic_error(
            "retained Acb basis sources lost their common Taylor prefix");
      common_width = width;
    }
    if (!common_width.has_value() ||
        *common_width - 1 >
            std::numeric_limits<std::uint32_t>::max())
      throw std::logic_error(
          "retained Acb matching Taylor prefix is invalid");
    return static_cast<std::uint32_t>(*common_width - 1);
  }

  std::optional<std::uint32_t>
  incoming_matching_taylor_complete_max() const {
    const auto* raw =
        incoming_source_.if_contains("matching_taylor_width");
    if (raw == nullptr) return std::nullopt;
    const auto width =
        as_u64(*raw, "retained Acb incoming matching Taylor width");
    if (width == 0 ||
        width - 1 > std::numeric_limits<std::uint32_t>::max())
      throw std::logic_error(
          "retained Acb incoming matching Taylor prefix is invalid");
    return static_cast<std::uint32_t>(width - 1);
  }

  bool certified_for_materialization() const {
    return !refined_.residual_history.empty() &&
        refined_.residual_history.back().verdict ==
            AcbMatchingResidualVerdict::Pass &&
        refined_.residual_history.back().complete_through_required;
  }

  AcbMatchingResidualVerdict final_residual_verdict() const {
    if (refined_.residual_history.empty())
      throw std::logic_error("retained Acb match has no residual history");
    return refined_.residual_history.back().verdict;
  }

  bool residual_complete_through_required() const {
    return !refined_.residual_history.empty() &&
        refined_.residual_history.back().complete_through_required;
  }

  std::int64_t contiguous_residual_pass_prefix() const {
    if (refined_.residual_history.empty())
      throw std::logic_error("retained Acb match has no residual history");
    const auto& diagnostics = refined_.residual_history.back();
    const auto upper = std::min(diagnostics.complete_window.complete_max,
                                diagnostics.required_complete_max);
    std::int64_t prefix =
        static_cast<std::int64_t>(diagnostics.complete_window.min_power) - 1;
    for (std::int64_t power = diagnostics.complete_window.min_power;
         power <= upper; ++power) {
      std::size_t passed = 0;
      for (const auto& coefficient : diagnostics.coefficients)
        if (coefficient.epsilon_power == power &&
            coefficient.verdict == AcbMatchingResidualVerdict::Pass)
          ++passed;
      if (passed != dimension_) break;
      prefix = power;
    }
    return prefix;
  }

  bool has_better_inconclusive_certificate_than(
      const StoredRefinedAcbMatch& other) const {
    const auto statistics = [](const StoredRefinedAcbMatch& match) {
      const auto& diagnostics = match.refined_.residual_history.back();
      std::size_t passed = 0;
      std::size_t failed = 0;
      for (const auto& coefficient : diagnostics.coefficients) {
        if (coefficient.epsilon_power >
            diagnostics.required_complete_max)
          continue;
        passed += coefficient.verdict ==
            AcbMatchingResidualVerdict::Pass;
        failed += coefficient.verdict ==
            AcbMatchingResidualVerdict::Fail;
      }
      const auto* raw_width =
          match.incoming_source_.if_contains("matching_taylor_width");
      const auto width = raw_width == nullptr
          ? std::size_t{0}
          : static_cast<std::size_t>(as_u64(
                *raw_width, "retained Acb matching Taylor width"));
      return std::tuple{
          match.contiguous_residual_pass_prefix(), failed, passed,
          diagnostics.complete_window.complete_max, width};
    };
    const auto [prefix, failed, passed, complete_max, width] =
        statistics(*this);
    const auto [other_prefix, other_failed, other_passed,
                other_complete_max, other_width] = statistics(other);
    if (prefix != other_prefix) return prefix > other_prefix;
    if (failed != other_failed) return failed < other_failed;
    if (passed != other_passed) return passed > other_passed;
    if (complete_max != other_complete_max)
      return complete_max > other_complete_max;
    return width > other_width;
  }

  void replace_elapsed_ms(double elapsed_ms) { elapsed_ms_ = elapsed_ms; }

  void replace_normal_frame_attempt(json::object report) {
    if (report.if_contains("schema") == nullptr ||
        report.if_contains("status") == nullptr)
      throw std::logic_error(
          "Acb normal-frame attempt report lacks its schema or status");
    normal_frame_attempt_ = std::move(report);
  }

  std::int32_t certified_complete_max() const {
    if (!certified_for_materialization())
      throw std::domain_error(
          "Acb match has no passing required residual certificate");
    const auto& diagnostics = refined_.residual_history.back();
    auto certified = diagnostics.complete_window.min_power - 1;
    for (std::int64_t power = diagnostics.complete_window.min_power;
         power <= diagnostics.complete_window.complete_max; ++power) {
      std::size_t passed = 0;
      for (const auto& coefficient : diagnostics.coefficients)
        if (coefficient.epsilon_power == power &&
            coefficient.verdict == AcbMatchingResidualVerdict::Pass)
          ++passed;
      if (passed != dimension_) break;
      certified = static_cast<std::int32_t>(power);
    }
    if (certified < diagnostics.required_complete_max)
      throw std::logic_error(
          "passing Acb residual did not retain its required certified prefix");
    return certified;
  }

  double elapsed_ms() const { return elapsed_ms_; }

  json::object checkpoint_record() const override {
    if (refined_.residual_history.empty())
      throw std::logic_error(
          "cannot checkpoint an Acb match without residual history");
    json::array basis;
    basis.reserve(basis_sources_.size());
    for (const auto& source : basis_sources_) {
      auto checkpoint_source = source;
      checkpoint_source.erase("source_operator_reference");
      basis.push_back(std::move(checkpoint_source));
    }
    auto checkpoint_incoming = incoming_source_;
    checkpoint_incoming.erase("source_operator_reference");
    json::value checkpoint_exact_shadow_basis = nullptr;
    if (exact_shadow_factorized_basis_.has_value()) {
      json::array columns;
      columns.reserve(exact_shadow_factorized_basis_->size());
      for (const auto& column : *exact_shadow_factorized_basis_)
        columns.push_back(checkpoint_local_solution_record(column));
      checkpoint_exact_shadow_basis = std::move(columns);
    }
    json::array history;
    history.reserve(refined_.residual_history.size());
    for (const auto& residual : refined_.residual_history)
      history.push_back(checkpoint_acb_match_residual_record(residual));
    json::object record{
        {"schema", "diffexp2-retained-acb-match-v6"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"exact_lattice_identity", exact_lattice_identity_},
        {"exact_lattice_provenance_identity",
         exact_lattice_provenance_identity_},
        {"exact_lattice_canonical_witness",
         exact_lattice_witness_record_},
        {"basis_sources", std::move(basis)},
        {"incoming_source", std::move(checkpoint_incoming)},
        {"basis_chart", basis_chart_},
        {"incoming_chart", incoming_chart_},
        {"basis_point_exact", basis_point_},
        {"incoming_point_exact", incoming_point_},
        {"physical_match_point_exact", physical_point_},
        {"matching_frame_identity", matching_frame_identity_},
        {"residual_frame_identity", residual_frame_identity_},
        {"epsilon",
         json::object{{"min", requested_window_.min_power},
                      {"max", requested_window_.complete_max},
                      {"required_complete_max", required_complete_max_}}},
        {"dimension", dimension_},
        {"relative_tolerance", relative_tolerance_},
        {"max_refinement_steps", max_refinement_steps_},
        {"exact_right_materialization_transformation",
         exact_right_materialization_transformation_.has_value()
             ? json::value(checkpoint_exact_laurent_matrix_record(
                   *exact_right_materialization_transformation_))
             : json::value(nullptr)},
        {"acb_right_materialization_preconditioner",
         acb_right_materialization_preconditioner_.has_value()
             ? json::value(checkpoint_acb_laurent_matrix_record(
                   *acb_right_materialization_preconditioner_))
             : json::value(nullptr)},
        {"terminal_normal_frame_right_transformation",
         terminal_normal_frame_right_transformation_.has_value()
             ? json::value(checkpoint_frame_matrix_record(
                   *terminal_normal_frame_right_transformation_))
             : json::value(nullptr)},
        {"terminal_normal_frame_exact_right_transformation",
         terminal_normal_frame_right_transformation_.has_value()
             ? json::value(checkpoint_exact_laurent_matrix_record(
                   exact_saturation_.transformation))
             : json::value(nullptr)},
        {"exact_shadow_factorized_basis",
         std::move(checkpoint_exact_shadow_basis)},
        {"exact_shadow_extra_precision_bits",
         exact_shadow_extra_precision_bits_},
        {"refined",
         json::object{
             {"transformed_weights",
              checkpoint_frame_vector_record(refined_.transformed_weights)},
             {"weights", checkpoint_frame_vector_record(refined_.weights)},
             {"residual", checkpoint_frame_vector_record(refined_.residual)},
             {"residual_history", std::move(history)},
             {"refinement_steps", refined_.refinement_steps},
             {"factorization_preconditioner",
              refined_.factorization_preconditioner},
             {"factorized_authoritative_rhs",
              refined_.factorized_authoritative_rhs},
             {"factorized_authoritative_rhs_columns",
              refined_.factorized_authoritative_rhs_columns}}},
        {"elapsed_ms", elapsed_ms_}};
    if (!residual_certificate_identity_.empty())
      record["residual_certificate_identity"] =
          residual_certificate_identity_;
    return record;
  }

 private:
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::string exact_lattice_identity_;
  std::string exact_lattice_provenance_identity_;
  std::string exact_lattice_witness_record_;
  std::string saturation_witness_schema_;
  std::vector<json::object> basis_sources_;
  json::object incoming_source_;
  std::string basis_chart_;
  std::string incoming_chart_;
  std::string basis_point_;
  std::string incoming_point_;
  std::string physical_point_;
  std::string matching_frame_identity_;
  std::string residual_frame_identity_;
  std::string residual_certificate_identity_;
  EpsilonWindow requested_window_;
  std::int32_t required_complete_max_ = 0;
  std::uint32_t dimension_ = 0;
  std::string relative_tolerance_;
  std::size_t max_refinement_steps_ = 0;
  EpsilonLatticeSaturationResult<Rational> exact_saturation_;
  std::optional<ExactLaurentMatrix<Rational>>
      exact_right_materialization_transformation_;
  std::optional<ExactLaurentMatrix<ComplexBall>>
      acb_materialization_right_transformation_;
  std::optional<ExactLaurentMatrix<ComplexBall>>
      acb_right_materialization_preconditioner_;
  std::optional<FiniteLaurentMatrix<ComplexBall>>
      terminal_normal_frame_right_transformation_;
  std::optional<std::vector<LocalSolution<ComplexBall>>>
      exact_shadow_factorized_basis_;
  std::size_t exact_shadow_extra_precision_bits_ = 0;
  std::shared_ptr<const PhysicalEquationOwnerBase>
      exact_shadow_equation_owner_;
  RefinedAcbLaurentMatch refined_;
  json::object normal_frame_attempt_{
      {"schema", "diffexp2-acb-normal-frame-attempt-v1"},
      {"status", "runtime-diagnostic-unavailable"}};
  double elapsed_ms_ = 0.0;
};

struct ParsedExactEvaluatedLattice {
  std::string witness_schema;
  std::string identity;
  std::string canonical_witness;
  EpsilonLatticeSaturationResult<Rational> saturation;
  bool exact_formal_negative_coefficients_are_zero = false;
  std::optional<ExactLaurentMatrix<ComplexBall>> acb_transformation;
  // T is certified directly against the retained physical Rational-shadow
  // columns F.  Preserve the exact local association F*T before any
  // coefficient is specialized.  The SCC right frame R is only a bounded
  // numerical proposal preconditioner: composing its finite truncation with
  // R^-1*T must never replace the exact physical transformation T.
  std::optional<std::vector<LocalSolution<Rational>>>
      exact_shadow_fused_local_basis;
  std::shared_ptr<const PhysicalEquationOwnerBase>
      exact_shadow_equation_owner;
};

json::array checkpoint_acb_laurent_matrix_record(
    const ExactLaurentMatrix<ComplexBall>& matrix) {
  json::array rows;
  rows.reserve(matrix.size());
  for (const auto& row : matrix) {
    json::array encoded_row;
    encoded_row.reserve(row.size());
    for (const auto& polynomial : row) {
      json::array terms;
      terms.reserve(polynomial.terms().size());
      for (const auto& [power, coefficient] : polynomial.terms())
        terms.push_back(json::object{
            {"power", power},
            {"coefficient", checkpoint_ball_record(coefficient)}});
      encoded_row.push_back(std::move(terms));
    }
    rows.push_back(std::move(encoded_row));
  }
  return rows;
}

ExactLaurentMatrix<ComplexBall> parse_checkpoint_acb_laurent_matrix(
    const json::value& raw, std::uint32_t dimension, const char* label) {
  const auto& rows = as_array(raw, label);
  if (rows.size() != dimension)
    throw std::invalid_argument(
        std::string(label) + " row count differs from its dimension");
  ExactLaurentMatrix<ComplexBall> matrix;
  matrix.reserve(rows.size());
  for (const auto& raw_row : rows) {
    const auto& row = as_array(raw_row, label);
    if (row.size() != dimension)
      throw std::invalid_argument(
          std::string(label) + " is not square");
    std::vector<ExactLaurentPolynomial<ComplexBall>> parsed_row;
    parsed_row.reserve(row.size());
    for (const auto& raw_entry : row) {
      ExactLaurentPolynomial<ComplexBall> polynomial;
      std::optional<std::int32_t> previous_power;
      for (const auto& raw_term : as_array(raw_entry, label)) {
        const auto& term = as_object(raw_term, label);
        require_exact_keys(term, {"power", "coefficient"}, label);
        const auto power = as_i32(term.at("power"), label);
        if (previous_power.has_value() && power <= *previous_power)
          throw std::invalid_argument(
              std::string(label) + " powers are not strictly increasing");
        polynomial.add_term(
            power, parse_checkpoint_ball(term.at("coefficient"), label));
        previous_power = power;
      }
      parsed_row.push_back(std::move(polynomial));
    }
    matrix.push_back(std::move(parsed_row));
  }
  return matrix;
}

ParsedExactEvaluatedLattice parse_exact_evaluated_lattice(
    const json::value& raw, std::uint32_t dimension, EpsilonWindow window,
    const std::string& context) {
  const auto& object = as_object(raw, "exact evaluated epsilon lattice");
  if (object.size() != 3 || object.if_contains("schema") == nullptr ||
      object.if_contains("identity") == nullptr ||
      object.if_contains("evaluated_basis") == nullptr)
    throw std::invalid_argument(
        "exact_lattice accepts exactly schema, identity, and evaluated_basis");
  if (required_string(object, "schema") != kExactEvaluatedLatticeSchema)
    throw std::invalid_argument("unsupported exact evaluated lattice schema");
  const auto identity = required_string(object, "identity");
  if (identity.empty())
    throw std::invalid_argument("exact evaluated lattice identity is empty");
  const auto& rows = as_array(object.at("evaluated_basis"),
                              "exact evaluated lattice basis");
  if (rows.size() != dimension)
    throw std::invalid_argument(
        "exact evaluated lattice must have one row per local component");

  FiniteLaurentMatrix<Rational> basis;
  basis.reserve(dimension);
  json::array canonical_rows;
  canonical_rows.reserve(dimension);
  for (std::uint32_t row = 0; row < dimension; ++row) {
    const auto& columns = as_array(rows[row],
                                   "exact evaluated lattice basis row");
    if (columns.size() != dimension)
      throw std::invalid_argument(
          "exact evaluated lattice basis must be square");
    FiniteLaurentVector<Rational> parsed_row;
    parsed_row.reserve(dimension);
    json::array canonical_columns;
    canonical_columns.reserve(dimension);
    for (std::uint32_t column = 0; column < dimension; ++column) {
      const auto& frame = as_object(columns[column],
                                    "exact evaluated lattice frame");
      if (frame.size() != 3 || frame.if_contains("min") == nullptr ||
          frame.if_contains("max") == nullptr ||
          frame.if_contains("coefficients") == nullptr)
        throw std::invalid_argument(
            "exact evaluated lattice frame accepts exactly min, max, and coefficients");
      const EpsilonWindow frame_window{
          as_i32(frame.at("min"), "exact lattice epsilon minimum"),
          as_i32(frame.at("max"), "exact lattice epsilon maximum")};
      (void)frame_window.width();
      if (frame_window.min_power != window.min_power ||
          frame_window.complete_max != window.complete_max)
        throw std::invalid_argument(
            "every exact evaluated lattice frame must equal the matching work window");
      const auto& raw_coefficients = as_array(
          frame.at("coefficients"), "exact lattice coefficients");
      if (raw_coefficients.size() != frame_window.width())
        throw std::invalid_argument(
            "exact lattice coefficient count differs from its epsilon window");
      std::vector<Rational> coefficients;
      coefficients.reserve(raw_coefficients.size());
      json::array canonical_coefficients;
      canonical_coefficients.reserve(raw_coefficients.size());
      for (const auto& raw_coefficient : raw_coefficients) {
        auto coefficient = parse_scalar<Rational>(raw_coefficient);
        canonical_coefficients.emplace_back(coefficient.str());
        coefficients.push_back(std::move(coefficient));
      }
      parsed_row.emplace_back(frame_window, std::move(coefficients));
      canonical_columns.push_back(json::object{
          {"min", frame_window.min_power},
          {"max", frame_window.complete_max},
          {"coefficients", std::move(canonical_coefficients)}});
    }
    basis.push_back(std::move(parsed_row));
    canonical_rows.push_back(std::move(canonical_columns));
  }
  json::object canonical{
      {"schema", kExactEvaluatedLatticeSchema}, {"identity", identity},
      {"evaluated_basis", std::move(canonical_rows)}};
  auto saturation = saturate_finite_laurent_basis(
      basis, context + ":exact-lattice-saturation");
  return {kExactEvaluatedLatticeSchema, identity,
          json::serialize(canonical), std::move(saturation)};
}

EpsilonLatticeSaturationResult<Rational> unit_rational_saturation(
    std::uint32_t dimension, EpsilonWindow window,
    const std::string& context) {
  if (dimension == 0 || window.min_power > 0 || window.complete_max < 0)
    throw std::invalid_argument(
        context + ": the unit saturation requires epsilon^0 in a nonempty square window");
  FiniteLaurentMatrix<Rational> identity(
      dimension, FiniteLaurentVector<Rational>());
  for (std::uint32_t row = 0; row < dimension; ++row) {
    identity[row].reserve(dimension);
    for (std::uint32_t column = 0; column < dimension; ++column) {
      std::vector<Rational> coefficients(
          window.width(), Rational(0));
      if (row == column)
        coefficients[static_cast<std::size_t>(
            -static_cast<std::int64_t>(window.min_power))] =
            Rational(1);
      identity[row].emplace_back(window, std::move(coefficients));
    }
  }
  auto saturation = saturate_finite_laurent_basis(
      identity, context + ":unit-saturation");
  if (saturation.diagnostics.initial_leading_rank != dimension ||
      saturation.diagnostics.final_leading_rank != dimension ||
      saturation.diagnostics.normalized_determinant_valuation != 0 ||
      !saturation.diagnostics.actions.empty() ||
      std::any_of(
          saturation.diagnostics.initial_column_shifts.begin(),
          saturation.diagnostics.initial_column_shifts.end(),
          [](std::int32_t shift) { return shift != 0; }))
    throw std::logic_error(
        context + ": internally constructed unit saturation is not the identity transformation");
  return saturation;
}

json::object validate_native_unit_saturation_request(
    const json::value& raw, const std::string& context) {
  auto request = as_object(raw, "native Acb unit-saturation request");
  const auto schema = required_string(request, "schema");
  const bool compact = schema == kNativeUnitSaturationCompactRequestSchema;
  if (compact)
    require_exact_keys(
        request,
        {"schema", "tile_plan", "tile_plan_checkpoint_identity",
         "arm", "match"},
        "compact native Acb unit-saturation request");
  else
    require_exact_keys(
        request,
        {"schema", "tile_plan", "tile_plan_checkpoint_identity",
         "tile_plan_provenance_identity", "arm", "match"},
        "native Acb unit-saturation request");
  if (!compact && schema != kNativeUnitSaturationRequestSchema)
    throw std::invalid_argument(
        context + ": unsupported native unit-saturation request schema");
  if (required_string(request, "tile_plan").empty() ||
      required_string(request, "tile_plan_checkpoint_identity").empty() ||
      (!compact &&
       required_string(request, "tile_plan_provenance_identity").empty()))
    throw std::invalid_argument(
        context + ": native unit-saturation request lost its plan binding");
  const auto arm = required_string(request, "arm");
  if (arm != "lower" && arm != "upper")
    throw std::invalid_argument(
        context + ": native unit-saturation request has an unknown arm");
  (void)as_u64(request.at("match"),
               "native unit-saturation match index");
  return request;
}

void require_ordinary_regular_basis_for_unit_saturation(
    const std::vector<std::shared_ptr<StoredLocalBase>>& basis) {
  if (basis.empty())
    throw std::invalid_argument(
        "native Acb unit-leading certification requires a nonempty basis");
  const auto dimension = as_u32(
      basis.front()->summary().at("dimension"),
      "native unit-leading basis dimension");
  if (basis.size() != dimension)
    throw std::invalid_argument(
        "native Acb unit-leading certification requires a square basis");
  for (const auto& column : basis) {
    const auto metadata = column->exact_analytic_metadata();
    const auto& sectors = as_array(
        metadata.at("sectors"), "native unit-leading sectors");
    if (sectors.size() != 1)
      throw std::domain_error(
          "native Acb unit-leading certification requires one ordinary sector per basis column");
    const auto& sector = as_object(
        sectors.front(), "native unit-leading sector");
    const auto a = parse_checkpoint_exact_descriptor(
        sector.at("a"), "native unit-leading a tag");
    const auto b = parse_checkpoint_exact_descriptor(
        sector.at("b"), "native unit-leading b tag");
    if (a.domain != ExactDomain::Rational ||
        b.domain != ExactDomain::Rational ||
        !(Rational(a.canonical) == Rational(0)) ||
        !(Rational(b.canonical) == Rational(0)) ||
        as_u32(sector.at("log_power"),
               "native unit-leading log power") != 0)
      throw std::domain_error(
          "native Acb unit-leading certification requires exact a=b=0, log_power=0 basis tags");
  }
}

ParsedExactEvaluatedLattice certify_native_unit_saturation(
    const json::value& raw_request,
    const FiniteLaurentMatrix<ComplexBall>& evaluated_basis,
    const std::vector<std::shared_ptr<StoredLocalBase>>& retained_basis,
    const std::vector<json::object>& basis_sources,
    const std::string& basis_point, const std::string& physical_point,
    EpsilonWindow window, const std::string& context) {
  require_ordinary_regular_basis_for_unit_saturation(retained_basis);
  const auto dimension = retained_basis.size();
  if (evaluated_basis.size() != dimension ||
      basis_sources.size() != dimension || window.min_power > 0 ||
      window.complete_max < 0)
    throw std::domain_error(
        context + ": native unit-leading certification requires a square actual basis complete through epsilon^0");
  for (std::size_t row = 0; row < dimension; ++row) {
    if (evaluated_basis[row].size() != dimension)
      throw std::domain_error(
          context + ": native unit-leading certification received a nonsquare actual basis");
    for (std::size_t column = 0; column < dimension; ++column) {
      const auto& frame = evaluated_basis[row][column];
      if (frame.complete_max() < 0)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            context + ": actual Acb basis is incomplete through epsilon^0",
            row, column, frame.complete_max());
      for (std::int64_t power = window.min_power; power < 0; ++power)
        if (!frame.coefficient(static_cast<std::int32_t>(power)).is_zero())
          throw MatchingArithmeticError(
              MatchingArithmeticErrorCode::InvalidSaturationLattice,
              context + ": actual Acb basis has a nonzero or zero-ambiguous negative epsilon coefficient",
              row, column, static_cast<std::int32_t>(power));
    }
  }
  const auto leading_rank =
      matching_detail::certify_full_rank_by_nonzero_pivots(
          matching_detail::epsilon_zero_matrix(
              evaluated_basis, context + ":actual-leading-frame"),
          context + ":actual-leading-rank");
  if (leading_rank != dimension)
    throw std::logic_error(
        context + ": full-rank proof returned the wrong dimension");

  auto native_request = validate_native_unit_saturation_request(
      raw_request, context);
  json::array proof_basis;
  proof_basis.reserve(basis_sources.size());
  for (const auto& source : basis_sources)
    proof_basis.push_back(checkpoint_matching_source_record(source));
  json::object proof_without_identity{
      {"schema", kNativeUnitSaturationProofSchema},
      {"native_request", std::move(native_request)},
      {"coefficient_domain", "acb"},
      {"basis", std::move(proof_basis)},
      {"basis_point_exact", basis_point},
      {"physical_match_point_exact", physical_point},
      {"epsilon", json::object{{"min", window.min_power},
                                {"max", window.complete_max}}},
      {"negative_epsilon_coefficients", "exact-singleton-zero"},
      {"leading_power", 0},
      {"leading_rank", dimension},
      {"leading_rank_certificate",
       "full-pivot-acb-pivots-exclude-zero"},
      {"transformation", "identity"}};
  const auto identity = json::serialize(
      canonical_json_value(proof_without_identity));
  auto proof = proof_without_identity;
  proof["identity"] = identity;
  return {kNativeUnitSaturationProofSchema, identity,
          json::serialize(canonical_json_value(proof)),
          unit_rational_saturation(
              static_cast<std::uint32_t>(dimension), window, context)};
}

ParsedExactEvaluatedLattice parse_native_unit_saturation_proof(
    const json::value& raw, std::uint32_t dimension, EpsilonWindow window,
    const std::vector<json::object>& expected_basis_sources,
    const std::string& expected_basis_point,
    const std::string& expected_physical_point,
    const std::string& context) {
  const auto& proof = as_object(raw, "native Acb unit-saturation proof");
  require_exact_keys(
      proof,
      {"schema", "identity", "native_request", "coefficient_domain",
       "basis", "basis_point_exact", "physical_match_point_exact",
       "epsilon", "negative_epsilon_coefficients", "leading_power",
       "leading_rank", "leading_rank_certificate", "transformation"},
      "native Acb unit-saturation proof");
  if (required_string(proof, "schema") !=
          kNativeUnitSaturationProofSchema ||
      required_string(proof, "coefficient_domain") != "acb" ||
      required_string(proof, "negative_epsilon_coefficients") !=
          "exact-singleton-zero" ||
      as_i32(proof.at("leading_power"),
             "native unit-leading power") != 0 ||
      as_u32(proof.at("leading_rank"),
             "native unit-leading rank") != dimension ||
      required_string(proof, "leading_rank_certificate") !=
          "full-pivot-acb-pivots-exclude-zero" ||
      required_string(proof, "transformation") != "identity")
    throw std::invalid_argument(
        context + ": native unit-saturation proof facts are inconsistent");
  (void)validate_native_unit_saturation_request(
      proof.at("native_request"), context);
  const auto& epsilon = as_object(
      proof.at("epsilon"), "native unit-saturation proof epsilon");
  require_exact_keys(epsilon, {"min", "max"},
                     "native unit-saturation proof epsilon");
  if (as_i32(epsilon.at("min"), "native proof epsilon minimum") !=
          window.min_power ||
      as_i32(epsilon.at("max"), "native proof epsilon maximum") !=
          window.complete_max ||
      required_string(proof, "basis_point_exact") !=
          expected_basis_point ||
      required_string(proof, "physical_match_point_exact") !=
          expected_physical_point)
    throw std::invalid_argument(
        context + ": native unit-saturation proof changed its point or epsilon binding");
  json::array expected_basis;
  expected_basis.reserve(expected_basis_sources.size());
  for (const auto& source : expected_basis_sources)
    expected_basis.push_back(checkpoint_matching_source_record(source));
  if (json::serialize(canonical_json_value(proof.at("basis"))) !=
      json::serialize(canonical_json_value(expected_basis)))
    throw std::invalid_argument(
        context + ": native unit-saturation proof changed its basis/checkpoint binding");
  auto identity_input = proof;
  const auto identity = required_string(proof, "identity");
  identity_input.erase("identity");
  if (identity.empty() ||
      json::serialize(canonical_json_value(identity_input)) != identity)
    throw std::invalid_argument(
        context + ": native unit-saturation proof identity is inconsistent");
  return {kNativeUnitSaturationProofSchema, identity,
          json::serialize(canonical_json_value(proof)),
          unit_rational_saturation(dimension, window, context)};
}

bool is_supported_acb_singular_scc_column_capability(
    const std::string& capability) {
  return capability == kAcbSingularScalarSCCColumnCapability ||
      capability == kAcbSingularJordanSCCColumnCapability;
}

json::object compact_native_singular_scc_saturation_request(
    const json::object& request) {
  const auto schema = required_string(request, "schema");
  if (schema == kNativeSingularSCCSaturationBoundedRequestSchema)
    return request;
  if (schema != kNativeSingularSCCSaturationRequestSchema &&
      schema != kNativeSingularSCCSaturationCompactRequestSchema)
    throw std::invalid_argument(
        "cannot compact an unsupported native singular-SCC request");
  const auto identity = matching_required_string_view(
      request, "receiving_scc_exact_identity");
  auto compact = request;
  compact["schema"] =
      kNativeSingularSCCSaturationBoundedRequestSchema;
  compact.erase("tile_plan_provenance_identity");
  compact.erase("receiving_scc_exact_identity");
  auto identity_reference =
      compact_matching_identity_reference(identity);
  identity_reference["schema"] =
      "diffexp2-retained-scc-identity-reference-v1";
  identity_reference["authority"] =
      "retained-native-scc-owner-validated";
  compact["receiving_scc_exact_identity_reference"] =
      std::move(identity_reference);
  return compact;
}

json::object validate_native_singular_scc_saturation_request(
    const json::value& raw, const std::string& context,
    const std::optional<std::string>& expected_session_configuration =
        std::nullopt,
    const std::optional<json::object>& expected_request = std::nullopt) {
  auto request = as_object(
      raw, "native Acb singular-SCC valuation-zero request");
  const auto schema = required_string(request, "schema");
  const bool compact =
      schema == kNativeSingularSCCSaturationCompactRequestSchema;
  const bool bounded =
      schema == kNativeSingularSCCSaturationBoundedRequestSchema;
  if (bounded)
    require_exact_keys(
        request,
        {"schema", "session_configuration_identity", "tile_plan",
         "tile_plan_checkpoint_identity", "arm", "match",
         "match_checkpoint_identity", "receiving_scc",
         "receiving_scc_exact_identity_reference",
         "receiving_execution_capability",
         "receiving_basis_point_exact", "receiving_basis_point_sign",
         "physical_match_point_exact", "receiving_rim"},
        "bounded native Acb singular-SCC valuation-zero request");
  else if (compact)
    require_exact_keys(
        request,
        {"schema", "session_configuration_identity", "tile_plan",
         "tile_plan_checkpoint_identity", "arm", "match",
         "match_checkpoint_identity", "receiving_scc",
         "receiving_scc_exact_identity", "receiving_execution_capability",
         "receiving_basis_point_exact", "receiving_basis_point_sign",
         "physical_match_point_exact", "receiving_rim"},
        "compact native Acb singular-SCC valuation-zero request");
  else
    require_exact_keys(
        request,
        {"schema", "session_configuration_identity", "tile_plan",
         "tile_plan_checkpoint_identity", "tile_plan_provenance_identity",
         "arm", "match", "match_checkpoint_identity", "receiving_scc",
         "receiving_scc_exact_identity", "receiving_execution_capability",
         "receiving_basis_point_exact", "receiving_basis_point_sign",
         "physical_match_point_exact", "receiving_rim"},
        "native Acb singular-SCC valuation-zero request");
  if (!compact && !bounded &&
      schema != kNativeSingularSCCSaturationRequestSchema)
    throw std::invalid_argument(
        context + ": unsupported native singular-SCC saturation request schema");
  for (const auto* key :
       {"session_configuration_identity", "tile_plan",
        "tile_plan_checkpoint_identity", "match_checkpoint_identity",
        "receiving_scc",
        "receiving_basis_point_exact", "physical_match_point_exact"})
    if (required_string(request, key).empty())
      throw std::invalid_argument(
          context + ": native singular-SCC saturation request lost its " +
          key + " binding");
  if (bounded) {
    const auto& reference = as_object(
        request.at("receiving_scc_exact_identity_reference"),
        "bounded native singular-SCC exact identity reference");
    require_exact_keys(
        reference,
        {"algorithm", "fingerprint", "identity_bytes", "schema",
         "authority"},
        "bounded native singular-SCC exact identity reference");
    if (required_string(reference, "schema") !=
            "diffexp2-retained-scc-identity-reference-v1" ||
        required_string(reference, "authority") !=
            "retained-native-scc-owner-validated" ||
        required_string(reference, "algorithm") != "fnv1a64-v1" ||
        required_string(reference, "fingerprint").empty() ||
        as_u64(reference.at("identity_bytes"),
               "bounded native singular-SCC identity bytes") == 0)
      throw std::invalid_argument(
          context +
          ": bounded native singular-SCC identity reference is invalid");
  } else if (required_string(
                 request, "receiving_scc_exact_identity").empty()) {
    throw std::invalid_argument(
        context +
        ": native singular-SCC saturation request lost its receiving_scc_exact_identity binding");
  }
  if (!compact && !bounded &&
      required_string(request, "tile_plan_provenance_identity").empty())
    throw std::invalid_argument(
        context + ": native singular-SCC saturation request lost its tile-plan provenance binding");
  const auto arm = required_string(request, "arm");
  if (arm != "lower" && arm != "upper")
    throw std::invalid_argument(
        context + ": native singular-SCC saturation request has an unknown arm");
  (void)as_u64(request.at("match"),
               "native singular-SCC saturation match index");
  const auto basis_point_sign = as_i32(
      request.at("receiving_basis_point_sign"),
      "native singular-SCC receiving basis-point sign");
  if (basis_point_sign != -1 && basis_point_sign != 1)
    throw std::invalid_argument(
        context +
        ": native singular-SCC receiving basis-point sign must be +1 or -1");
  const auto capability = required_string(
      request, "receiving_execution_capability");
  if (!is_supported_acb_singular_scc_column_capability(capability))
    throw std::domain_error(
        context + ": receiving SCC lacks a supported affine-Jordan Acb column capability");
  if (!request.at("receiving_rim").is_null()) {
    const auto rim = as_i32(
        request.at("receiving_rim"), "native singular-SCC receiving rim");
    if (rim != -1 && rim != 1)
      throw std::invalid_argument(
          context + ": native singular-SCC receiving rim must be +1, -1, or null");
  }
  if (expected_session_configuration.has_value() &&
      required_string(request, "session_configuration_identity") !=
          *expected_session_configuration)
    throw std::invalid_argument(
        context + ": native singular-SCC saturation request changed its stable session-configuration binding");
  if (expected_request.has_value()) {
    const auto expected = bounded
        ? compact_native_singular_scc_saturation_request(
              *expected_request)
        : *expected_request;
    if (json::serialize(canonical_json_value(request)) !=
        json::serialize(canonical_json_value(expected)))
      throw std::invalid_argument(
          context + ": native singular-SCC saturation request changed its retained plan/SCC binding");
  } else if (bounded) {
    throw std::invalid_argument(
        context +
        ": bounded native singular-SCC request requires its retained owner binding");
  }
  return request;
}

void validate_singular_scc_basis_sources(
    const json::object& native_request,
    const std::vector<json::object>& basis_sources,
    std::uint32_t dimension, const std::string& context) {
  if (dimension == 0 || basis_sources.size() != dimension)
    throw std::domain_error(
        context + ": singular-SCC valuation-zero proof requires one complete square basis");
  const auto expected_scc = required_string(native_request, "receiving_scc");
  const auto expected_identity = required_string(
      native_request, "receiving_scc_exact_identity");
  const auto capability = required_string(
      native_request, "receiving_execution_capability");
  const bool scalar = capability == kAcbSingularScalarSCCColumnCapability;
  const auto basis_point_sign = as_i32(
      native_request.at("receiving_basis_point_sign"),
      "native singular-SCC receiving basis-point sign");
  const json::value expected_effective_rim = basis_point_sign < 0
      ? native_request.at("receiving_rim") : json::value(nullptr);
  const auto* expected_column_schema = scalar
      ? "diffexp2-native-scc-acb-regular-singular-scalar-column-v1"
      : "diffexp2-native-scc-acb-regular-singular-jordan-column-v1";
  std::vector<std::uint8_t> seen(dimension, 0);
  for (std::size_t column = 0; column < basis_sources.size(); ++column) {
    const auto& source = basis_sources[column];
    if (as_u64(source.at("column"), "singular-SCC proof basis column") !=
        column)
      throw std::invalid_argument(
          context + ": singular-SCC basis sources are not in receiving column order");
    if (required_string(source, "chart") != expected_scc ||
        required_string(source, "source_operator_identity") !=
            expected_identity)
      throw std::invalid_argument(
          context + ": singular-SCC basis source is not owned by the receiving SCC");
    if (source.at("requested_imaginary_sign") !=
            native_request.at("receiving_rim") ||
        source.at("effective_imaginary_sign") != expected_effective_rim)
      throw std::invalid_argument(
          context + ": singular-SCC basis source changed its point-dependent requested/effective plan-selected rim");
    const auto* raw_provenance = source.if_contains("column_provenance");
    if (raw_provenance == nullptr)
      throw std::domain_error(
          context + ": singular-SCC basis column lacks certified SCC provenance");
    const auto& provenance = as_object(
        *raw_provenance, "singular-SCC basis column provenance");
    require_exact_keys(
        provenance,
        {"scc", "scc_exact_identity", "seed_block", "basis_index",
         "exact_column_identity"},
        "singular-SCC basis column provenance");
    const auto basis_index = as_u32(
        provenance.at("basis_index"), "singular-SCC provenance basis index");
    if (required_string(provenance, "scc") != expected_scc ||
        required_string(provenance, "scc_exact_identity") !=
            expected_identity ||
        basis_index >= dimension || seen[basis_index] != 0)
      throw std::domain_error(
          context + ": singular-SCC basis provenance is incomplete, duplicated, or belongs to another SCC");
    seen[basis_index] = 1;

    const auto exact_column_record = required_string(
        provenance, "exact_column_identity");
    if (exact_column_record.empty())
      throw std::invalid_argument(
          context + ": singular-SCC basis column lost its exact identity");
    const auto parsed_column = json::parse(exact_column_record);
    const auto& exact_column = as_object(
        parsed_column, "singular-SCC exact column identity");
    if (json::serialize(canonical_json_value(parsed_column)) !=
        exact_column_record)
      throw std::invalid_argument(
          context + ": singular-SCC exact column identity is not canonically encoded");
    const auto column_schema = required_string(exact_column, "schema");
    if (column_schema ==
        "diffexp2-native-scc-acb-rational-shadow-column-v1") {
      require_exact_keys(
          exact_column,
          {"schema", "scc_exact_identity", "basis_index", "seed_block",
           "rational_shadow_identity", "rational_source_column_identity",
           "pseudo_compensation"},
          "singular-SCC Rational-shadow exact column identity");
      if (required_string(exact_column, "scc_exact_identity") !=
              expected_identity ||
          as_u32(exact_column.at("basis_index"),
                 "singular-SCC shadow basis index") != basis_index ||
          as_u32(exact_column.at("seed_block"),
                 "singular-SCC shadow seed block") !=
              as_u32(provenance.at("seed_block"),
                     "singular-SCC provenance seed block") ||
          required_string(exact_column,
                          "rational_shadow_identity").empty() ||
          required_string(exact_column,
                          "rational_source_column_identity").empty() ||
          required_string(exact_column, "pseudo_compensation") !=
              "exact-rational-shadow-case-p-floor-certified-v1")
        throw std::domain_error(
            context + ": singular-SCC Rational-shadow column lost its exact CASE-P owner/floor certificate");
      const auto source_column_record = required_string(
          exact_column, "rational_source_column_identity");
      const auto source_column = json::parse(source_column_record);
      if (json::serialize(canonical_json_value(source_column)) !=
          source_column_record)
        throw std::invalid_argument(
            context + ": Rational-shadow source column identity is not canonical JSON");
    } else {
      if (scalar)
        require_exact_keys(
            exact_column,
            {"schema", "scc_exact_identity", "basis_index", "seed",
             "targets", "pseudo_compensation"},
            "singular-SCC scalar exact column identity");
      else
        require_exact_keys(
            exact_column,
            {"schema", "scc_exact_identity", "basis_index", "seed",
             "targets", "pseudo_compensation", "seed_local_component"},
            "singular-SCC Jordan exact column identity");
      if (column_schema != expected_column_schema ||
          required_string(exact_column, "scc_exact_identity") !=
              expected_identity ||
          as_u32(exact_column.at("basis_index"),
                 "singular-SCC exact column basis index") != basis_index ||
          (required_string(exact_column, "pseudo_compensation") != "none" &&
           required_string(exact_column, "pseudo_compensation") !=
               "exact-schedule-acb-ball-floor-certified-v1"))
        throw std::domain_error(
            context + ": singular-SCC basis column lacks a supported exact-schedule Acb compensation certificate");
      const auto& seed = as_object(
          exact_column.at("seed"), "singular-SCC exact column seed");
      if (as_u32(seed.at("block"),
                 "singular-SCC exact column seed block") !=
          as_u32(provenance.at("seed_block"),
                 "singular-SCC provenance seed block"))
        throw std::invalid_argument(
            context + ": singular-SCC exact column changed its seed-block binding");
      (void)as_array(exact_column.at("targets"),
                     "singular-SCC exact column targets");
    }
  }
  if (std::any_of(seen.begin(), seen.end(),
                  [](std::uint8_t value) { return value == 0; }))
    throw std::domain_error(
        context + ": singular-SCC basis provenance does not cover every canonical column");
}

std::shared_ptr<PhysicalEquationOwnerBase>
validate_retained_singular_scc_basis_sources(
    const json::object& native_request,
    const std::vector<std::shared_ptr<StoredLocalBase>>& retained_basis,
    const std::vector<json::object>& basis_sources,
    std::uint32_t dimension, const std::string& context) {
  if (dimension == 0 || retained_basis.size() != dimension ||
      basis_sources.size() != dimension)
    throw std::domain_error(
        context +
        ": retained singular-SCC proof requires one complete square basis");
  const auto expected_scc = required_string(native_request, "receiving_scc");
  const bool bounded_request =
      required_string(native_request, "schema") ==
      kNativeSingularSCCSaturationBoundedRequestSchema;
  const auto expected_identity = bounded_request
      ? retained_basis.front()->source_operator_identity()
      : required_string(native_request, "receiving_scc_exact_identity");
  if (expected_identity.empty())
    throw std::invalid_argument(
        context +
        ": retained singular-SCC basis lost its exact operator identity");
  if (bounded_request) {
    const auto& reference = as_object(
        native_request.at("receiving_scc_exact_identity_reference"),
        "bounded retained singular-SCC exact identity reference");
    if (required_string(reference, "algorithm") != "fnv1a64-v1" ||
        required_string(reference, "fingerprint") !=
            public_provenance_fingerprint(expected_identity) ||
        as_u64(reference.at("identity_bytes"),
               "bounded retained singular-SCC identity bytes") !=
            expected_identity.size())
      throw std::invalid_argument(
          context +
          ": bounded singular-SCC request does not match its retained exact operator owner");
  }
  const auto expected_capability = required_string(
      native_request, "receiving_execution_capability");
  const auto basis_point_sign = as_i32(
      native_request.at("receiving_basis_point_sign"),
      "retained singular-SCC receiving basis-point sign");
  const json::value expected_effective_rim = basis_point_sign < 0
      ? native_request.at("receiving_rim") : json::value(nullptr);
  std::vector<std::uint8_t> seen(dimension, 0);
  std::shared_ptr<PhysicalEquationOwnerBase> common_scc_owner;
  bool ownerless_basis = false;
  for (std::size_t column = 0; column < basis_sources.size(); ++column) {
    const auto& source = basis_sources[column];
    const auto& local = retained_basis[column];
    if (!local ||
        as_u64(source.at("column"),
               "retained singular-SCC proof basis column") != column ||
        required_string(source, "local") != local->handle() ||
        required_string(source, "chart") != local->source_chart() ||
        required_string(source, "source_operator_identity") !=
            local->source_operator_identity() ||
        required_string(source, "checkpoint_identity") !=
            local->checkpoint_identity())
      throw std::invalid_argument(
          context +
          ": retained singular-SCC basis source changed its strong local binding");
    const auto scc_owner = local->retained_equation_owner();
    const auto owner_dimension = scc_owner
        ? scc_owner->matching_scc_dimension() : std::nullopt;
    const auto* owner_capability = scc_owner
        ? scc_owner->matching_scc_column_execution_capability() : nullptr;
    if (!scc_owner) {
      if (common_scc_owner)
        throw std::invalid_argument(
            context +
            ": retained singular-SCC basis has inconsistent equation ownership");
      ownerless_basis = true;
    } else {
      if (ownerless_basis)
        throw std::invalid_argument(
            context +
            ": retained singular-SCC basis has inconsistent equation ownership");
      if (std::string(scc_owner->equation_owner_kind()) != "composite-scc")
        throw std::invalid_argument(
            context +
            ": retained singular-SCC basis equation owner is not a CompositeSCC");
      if (common_scc_owner &&
          common_scc_owner.get() != scc_owner.get())
        throw std::invalid_argument(
            context +
            ": retained singular-SCC basis columns do not share one live equation owner");
      if (scc_owner->equation_owner_handle() != expected_scc)
        throw std::invalid_argument(
            context +
            ": retained singular-SCC basis live owner changed its handle");
      if (scc_owner->equation_operator_identity() != expected_identity)
        throw std::invalid_argument(
            context +
            ": retained singular-SCC basis live owner changed its exact identity");
      if (!owner_dimension.has_value() || *owner_dimension != dimension)
        throw std::invalid_argument(
            context +
            ": retained singular-SCC basis live owner changed its dimension");
      if (owner_capability == nullptr ||
          std::string(owner_capability) != expected_capability)
        throw std::invalid_argument(
            context +
            ": retained singular-SCC basis live owner changed its execution capability");
      if (!common_scc_owner) common_scc_owner = scc_owner;
    }
    if (local->source_chart() != expected_scc ||
        local->source_operator_identity() != expected_identity ||
        source.at("requested_imaginary_sign") !=
            native_request.at("receiving_rim") ||
        source.at("effective_imaginary_sign") != expected_effective_rim ||
        source.at("analytic_metadata") != local->exact_analytic_metadata())
      throw std::invalid_argument(
          context +
          ": retained singular-SCC basis source changed its SCC, rim, or analytic binding");
    const auto& provenance = local->column_provenance();
    const auto* raw_provenance = source.if_contains("column_provenance");
    if (!provenance.has_value() || raw_provenance == nullptr ||
        *raw_provenance != provenance->encode() ||
        provenance->scc_handle != expected_scc ||
        provenance->scc_exact_identity != expected_identity ||
        provenance->basis_index >= dimension ||
        seen[provenance->basis_index] != 0 ||
        provenance->exact_column_identity.empty())
      throw std::domain_error(
          context +
          ": retained singular-SCC column owner is incomplete, duplicated, or changed");
    seen[provenance->basis_index] = 1;
  }
  if (std::any_of(seen.begin(), seen.end(),
                  [](std::uint8_t value) { return value == 0; }))
    throw std::domain_error(
        context +
        ": retained singular-SCC basis does not cover every canonical column");
  return common_scc_owner;
}

ParsedExactEvaluatedLattice certify_native_singular_scc_saturation(
    const json::value& raw_request,
    const FiniteLaurentMatrix<ComplexBall>& evaluated_basis,
    const std::vector<std::shared_ptr<StoredLocalBase>>& retained_basis,
    const std::vector<json::object>& basis_sources,
    const std::string& basis_point, const std::string& physical_point,
    EpsilonWindow window,
    const std::string& expected_session_configuration,
    const json::object& expected_native_request,
    const std::string& expected_checkpoint_identity,
    const std::string& context,
    bool prefer_retained_rational_shadow = true,
    const std::optional<ExactLaurentMatrix<ComplexBall>>&
        normalized_right_inverse = std::nullopt) {
  auto native_request = validate_native_singular_scc_saturation_request(
      raw_request, context, expected_session_configuration,
      expected_native_request);
  const auto native_request_schema =
      required_string(native_request, "schema");
  const bool compact_request =
      native_request_schema ==
          kNativeSingularSCCSaturationCompactRequestSchema ||
      native_request_schema ==
          kNativeSingularSCCSaturationBoundedRequestSchema;
  if (required_string(native_request, "match_checkpoint_identity") !=
          expected_checkpoint_identity ||
      required_string(native_request, "receiving_basis_point_exact") !=
          basis_point ||
      required_string(native_request, "physical_match_point_exact") !=
          physical_point)
    throw std::invalid_argument(
        context + ": native singular-SCC proof request changed its checkpoint or point binding");
  const auto dimension = retained_basis.size();
  if (dimension == 0 || evaluated_basis.size() != dimension ||
      basis_sources.size() != dimension || window.min_power > 0 ||
      window.complete_max < 0)
    throw std::domain_error(
        context + ": singular-SCC valuation-zero certification requires a square actual basis complete through epsilon^0");
  std::shared_ptr<PhysicalEquationOwnerBase> common_scc_owner;
  if (compact_request)
    common_scc_owner = validate_retained_singular_scc_basis_sources(
        native_request, retained_basis, basis_sources,
        static_cast<std::uint32_t>(dimension), context);
  else {
    validate_singular_scc_basis_sources(
        native_request, basis_sources,
        static_cast<std::uint32_t>(dimension), context);
    for (std::size_t column = 0; column < dimension; ++column) {
      const auto& provenance = retained_basis[column]->column_provenance();
      if (!provenance.has_value() ||
          json::serialize(canonical_json_value(provenance->encode())) !=
              json::serialize(canonical_json_value(
                  basis_sources[column].at("column_provenance"))))
        throw std::invalid_argument(
            context +
            ": singular-SCC proof source disagrees with its retained column owner");
    }
  }

  std::vector<std::shared_ptr<const RationalShadowColumnWitness>>
      shadow_witnesses;
  shadow_witnesses.reserve(dimension);
  for (const auto& column : retained_basis)
    shadow_witnesses.push_back(column->rational_shadow_witness());
  const bool shadow_basis = std::all_of(
      shadow_witnesses.begin(), shadow_witnesses.end(),
      [](const auto& witness) { return witness != nullptr; });
  std::size_t shadow_provenance_columns = 0;
  for (const auto& column : retained_basis) {
    const auto& provenance = column->column_provenance();
    if (!provenance.has_value()) continue;
    const auto parsed_column =
        json::parse(provenance->exact_column_identity);
    const auto& exact_column = as_object(
        parsed_column,
        "retained singular-SCC exact column identity");
    shadow_provenance_columns +=
        required_string(exact_column, "schema") ==
        "diffexp2-native-scc-acb-rational-shadow-column-v1";
  }
  if (prefer_retained_rational_shadow &&
      shadow_provenance_columns != 0 &&
      shadow_provenance_columns != dimension)
    throw std::invalid_argument(
        context +
        ": singular-SCC basis mixes Rational-shadow and native Acb column provenance");
  // A post-match checkpoint intentionally retains compact column provenance
  // but not another copy of every Rational local.  In that case recompute an
  // Acb proposal from the restored columns below.  Fresh transport still
  // takes the live exact-shadow route.
  if (shadow_basis && prefer_retained_rational_shadow) {
    for (std::size_t column = 0; column < dimension; ++column) {
      const auto& witness = *shadow_witnesses[column];
      if (!witness.solution || witness.rational_shadow_identity.empty() ||
          witness.source_column_identity.empty() ||
          witness.target_column_identity.empty() ||
          !witness.exact_equation_owner)
        throw std::invalid_argument(
            context + ": Rational-shadow saturation witness is incomplete");
      if (column != 0 &&
          witness.exact_equation_owner.get() !=
              shadow_witnesses.front()->exact_equation_owner.get())
        throw std::invalid_argument(
            context +
            ": Rational-shadow columns disagree on their exact equation owner");
      const auto& exact_column_identity =
          retained_basis[column]->column_provenance()
              ->exact_column_identity;
      if (compact_request) {
        const auto* owner_shadow_identity = common_scc_owner
            ? common_scc_owner->matching_scc_rational_shadow_identity()
            : nullptr;
        if (owner_shadow_identity == nullptr ||
            witness.rational_shadow_identity !=
                *owner_shadow_identity ||
            exact_column_identity != witness.target_column_identity)
          throw std::invalid_argument(
              context +
              ": retained Rational-shadow witness changed its target-column binding");
      } else {
        const auto parsed_column = json::parse(exact_column_identity);
        const auto& exact_column = as_object(
            parsed_column,
            "Rational-shadow target column identity");
        if (required_string(exact_column, "rational_shadow_identity") !=
                witness.rational_shadow_identity ||
            required_string(exact_column,
                            "rational_source_column_identity") !=
                witness.source_column_identity)
          throw std::invalid_argument(
              context +
              ": Rational-shadow saturation witness changed its target provenance binding");
      }
      validate_local_solution(*witness.solution, false);
    }
    auto saturation = rational_shadow_formal_saturation(
        shadow_witnesses, window,
        context + ":Rational-shadow-formal-saturation");
    auto proof_basis =
        compact_matching_source_references(basis_sources);
    json::array encoded_valuations;
    for (const auto valuation :
         saturation.diagnostics.initial_column_valuations)
      encoded_valuations.push_back(valuation);
    json::object proof_without_identity{
        {"schema",
         kNativeSingularSCCSaturationCompactProofSchema},
        {"native_request",
         compact_native_singular_scc_saturation_request(
             native_request)},
        {"coefficient_domain", "acb"},
        {"basis", std::move(proof_basis)},
        {"basis_point_exact", basis_point},
        {"physical_match_point_exact", physical_point},
        {"epsilon", json::object{{"min", window.min_power},
                                  {"max", window.complete_max}}},
        {"negative_epsilon_coefficients",
         "exact-rational-shadow-formal-valuation"},
        {"leading_power", 0},
        {"leading_rank", dimension},
        {"leading_rank_certificate",
         "deferred-acb-factorization-after-exact-column-shifts"},
        {"column_provenance_certificate",
         "complete-target-bound-rational-shadow-case-p-floor-certified"},
        {"determinant_valuation",
         saturation.diagnostics.normalized_determinant_valuation},
        {"transformation", "exact-rational-shadow-column-monomials"},
        {"column_valuations", std::move(encoded_valuations)}};
    const auto identity = json::serialize(
        canonical_json_value(proof_without_identity));
    auto proof = proof_without_identity;
    proof["identity"] = identity;
    std::optional<ExactLaurentMatrix<ComplexBall>>
        normalized_transformation;
    std::optional<std::vector<LocalSolution<Rational>>>
        exact_fused_local_basis;
    if (normalized_right_inverse.has_value()) {
      normalized_transformation = multiply_exact_laurent_matrices(
          *normalized_right_inverse,
          specialize_exact_rational_laurent_matrix_to_acb(
              saturation.transformation));
      std::vector<const LocalSolution<Rational>*> exact_basis;
      exact_basis.reserve(dimension);
      for (const auto& witness : shadow_witnesses)
        exact_basis.push_back(witness->solution.get());
      exact_fused_local_basis =
          right_transform_local_basis_exact(
              exact_basis, saturation.transformation,
              context + ":exact-physical-F*T");
    }
    return {kNativeSingularSCCSaturationCompactProofSchema, identity,
            json::serialize(canonical_json_value(proof)),
            std::move(saturation),
            !normalized_right_inverse.has_value(),
            std::move(normalized_transformation),
            std::move(exact_fused_local_basis),
            shadow_witnesses.front()->exact_equation_owner};
  }
  if (prefer_retained_rational_shadow &&
      std::any_of(shadow_witnesses.begin(), shadow_witnesses.end(),
                  [](const auto& witness) { return witness != nullptr; }))
    throw std::invalid_argument(
        context + ": singular-SCC basis mixes Rational-shadow and native Acb columns");
  for (std::size_t row = 0; row < dimension; ++row)
    if (evaluated_basis[row].size() != dimension)
      throw std::domain_error(
          context + ": singular-SCC Laurent certification received a nonsquare actual basis");
  for (std::size_t row = 0; row < dimension; ++row)
    for (std::size_t column = 0; column < dimension; ++column)
      if (evaluated_basis[row][column].complete_max() < 0)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            context + ": actual singular-SCC Acb basis is incomplete through epsilon^0",
            row, column, evaluated_basis[row][column].complete_max());

  // Canonicalize only enclosures whose *whole ball* lies below the structural
  // floor, then determine the column Laurent valuations from that completed
  // evaluated basis.  This is not a midpoint zero decision: every chopped
  // ball is rigorously negligible, and the final matching residual below is
  // still evaluated against `evaluated_basis`, not this canonical copy.
  //
  // Normalizing the column valuations before factorization is essential for
  // honest finite-window arithmetic.  Direct Gaussian elimination on the
  // raw polar basis repeatedly combines negative minima and can discard
  // several perfectly known upper coefficients even when the normalized
  // determinant valuation is zero.
  // This copy is only a proposal for an invertible Levelt coordinate
  // transformation; the authoritative residual is evaluated against the
  // untouched Acb basis below.  Keep the conservative historical floor for
  // the numeric fallback.  Large singular systems must use their retained
  // exact Rational-shadow actions: lowering this floor merely postpones
  // catastrophic interval-nullspace loss (banana4 advances from action 18
  // to action 24 and then widens to radius 2^196).
  constexpr int structural_chop_digits = 100;
  auto structural_basis = evaluated_basis;
  for (auto& row : structural_basis)
    for (auto& frame : row) {
      auto coefficients = frame.coefficients();
      for (auto& coefficient : coefficients)
        coefficient = ScalarTraits<ComplexBall>::canonicalized(
            coefficient, structural_chop_digits);
      frame = EpsilonFrame<ComplexBall>(
          frame.window(), std::move(coefficients));
    }
  std::vector<std::int32_t> column_valuations;
  std::optional<EpsilonLatticeSaturationResult<Rational>>
      rational_saturation;
  std::optional<ExactLaurentMatrix<ComplexBall>> acb_transformation;
  try {
    column_valuations.resize(dimension);
    for (std::size_t column = 0; column < dimension; ++column) {
      std::optional<std::int32_t> valuation;
      std::optional<std::int32_t> ambiguous_floor;
      for (std::size_t row = 0; row < dimension; ++row) {
        const auto leading = matching_detail::certified_laurent_leading_power(
            structural_basis[row][column]);
        if (leading.power.has_value() &&
            (!valuation.has_value() || *leading.power < *valuation))
          valuation = leading.power;
        if (leading.first_ambiguous_power.has_value() &&
            (!ambiguous_floor.has_value() ||
             *leading.first_ambiguous_power < *ambiguous_floor))
          ambiguous_floor = leading.first_ambiguous_power;
      }
      if (!valuation.has_value())
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::SingularOrIncompleteSystem,
            context +
                ": no certified nonzero coefficient determines an Acb column valuation",
            std::nullopt, column);
      // An unresolved ball below this certified nonzero candidate does not
      // invalidate the monomial change of coordinates: every integer column
      // shift is exactly invertible.  It only prevents calling the candidate
      // an exact valuation.  The transformed Laurent factorization and final
      // residual remain the proof authorities.
      (void)ambiguous_floor;
      column_valuations[column] = *valuation;
    }

    auto candidate_acb_saturation = saturate_finite_laurent_basis(
        structural_basis, context + ":certified-Acb-Levelt-saturation",
        structural_chop_digits);
    if (candidate_acb_saturation.diagnostics.final_leading_rank !=
            dimension)
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::InvalidSaturationLattice,
          context +
              ": certified Acb Levelt saturation did not produce a valuation-zero full-rank lattice; determinant_valuation=" +
              std::to_string(candidate_acb_saturation.diagnostics
                                 .normalized_determinant_valuation) +
              "; initial_rank=" +
              std::to_string(candidate_acb_saturation.diagnostics
                                 .initial_leading_rank) +
              "; final_rank=" +
              std::to_string(candidate_acb_saturation.diagnostics
                                 .final_leading_rank) +
              "; actions=" +
              std::to_string(
                  candidate_acb_saturation.diagnostics.actions.size()));

    FiniteLaurentMatrix<Rational> monomial_basis(
        dimension, FiniteLaurentVector<Rational>());
    for (std::size_t row = 0; row < dimension; ++row) {
      monomial_basis[row].reserve(dimension);
      for (std::size_t column = 0; column < dimension; ++column) {
        const auto valuation = column_valuations[column];
        if (valuation < window.min_power || valuation > window.complete_max)
          throw MatchingArithmeticError(
              MatchingArithmeticErrorCode::InsufficientCompleteWindow,
              context +
                  ": certified Acb column valuation lies outside the matching window",
              row, column, valuation);
        std::vector<Rational> coefficients(window.width(), Rational(0));
        if (row == column)
          coefficients[static_cast<std::size_t>(
              static_cast<std::int64_t>(valuation) - window.min_power)] =
              Rational(1);
        monomial_basis[row].emplace_back(window, std::move(coefficients));
      }
    }
    auto candidate_rational_saturation = saturate_finite_laurent_basis(
        monomial_basis, context + ":certified-Acb-monomial-provenance");
    if (candidate_rational_saturation.diagnostics
                .normalized_determinant_valuation != 0 ||
        candidate_rational_saturation.diagnostics.initial_leading_rank !=
            dimension ||
        candidate_rational_saturation.diagnostics.final_leading_rank !=
            dimension ||
        !candidate_rational_saturation.diagnostics.actions.empty())
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::InvalidSaturationLattice,
          context +
              ": certified column valuations do not define a monomial valuation-zero lattice");

    rational_saturation = std::move(candidate_rational_saturation);
    acb_transformation = std::move(
        candidate_acb_saturation.transformation);
  } catch (const MatchingArithmeticError& error) {
    // Some otherwise well-conditioned small SCCs still have correlated
    // leading balls whose rank is inconclusive after the strict structural
    // chop.  Their direct Laurent factorization remains a fully certified
    // fallback; only the final residual, never this attempted normalization,
    // authorizes materialization.
    if (std::getenv("DE2_DIAGNOSTIC_TERMINAL_STATE") != nullptr)
      std::fprintf(
          stderr,
          "[diffexp2 terminal diagnostic] Acb Levelt saturation rejected: %s\n",
          error.what());
  }
  const bool monomial_saturation = rational_saturation.has_value() &&
                                   acb_transformation.has_value();
  if (!monomial_saturation) {
    (void)factor_preconditioned_acb_finite_laurent_system(
        evaluated_basis, context + ":certified-direct-Acb-Laurent-pivots");
    auto proof_basis =
        compact_matching_source_references(basis_sources);
    json::object proof_without_identity{
        {"schema",
         kNativeSingularSCCSaturationCompactProofSchema},
        {"native_request",
         compact_native_singular_scc_saturation_request(
             native_request)},
        {"coefficient_domain", "acb"},
        {"basis", std::move(proof_basis)},
        {"basis_point_exact", basis_point},
        {"physical_match_point_exact", physical_point},
        {"epsilon", json::object{{"min", window.min_power},
                                  {"max", window.complete_max}}},
        {"negative_epsilon_coefficients",
         "retained-certified-Acb-Laurent-frames"},
        {"leading_power", 0},
        {"leading_rank", dimension},
        {"leading_rank_certificate",
         "full-pivot-acb-Laurent-pivots-exclude-zero"},
        {"column_provenance_certificate",
         "complete-exact-schedule-receiving-scc-acb"},
        {"determinant_valuation", 0},
        {"transformation", "direct-acb-Laurent-pivoting"}};
    const auto identity = json::serialize(
        canonical_json_value(proof_without_identity));
    auto proof = proof_without_identity;
    proof["identity"] = identity;
    return {kNativeSingularSCCSaturationCompactProofSchema, identity,
            json::serialize(canonical_json_value(proof)),
            unit_rational_saturation(
                static_cast<std::uint32_t>(dimension), window, context)};
  }

  auto proof_basis =
      compact_matching_source_references(basis_sources);
  json::object proof_without_identity{
      {"schema",
       kNativeSingularSCCSaturationCompactProofSchema},
      {"native_request",
       compact_native_singular_scc_saturation_request(
           native_request)},
      {"coefficient_domain", "acb"},
      {"basis", std::move(proof_basis)},
      {"basis_point_exact", basis_point},
      {"physical_match_point_exact", physical_point},
      {"epsilon", json::object{{"min", window.min_power},
                                {"max", window.complete_max}}},
      {"negative_epsilon_coefficients",
       "acb-certified-finite-laurent-lattice"},
      {"leading_power", 0},
      {"leading_rank", dimension},
      {"leading_rank_certificate",
       "full-acb-Levelt-saturation-pivots-exclude-zero"},
      {"column_provenance_certificate",
       "complete-exact-schedule-receiving-scc-acb"},
      {"determinant_valuation", 0},
      {"transformation", "certified-acb-Levelt-lattice-saturation"},
      {"column_valuations", json::array{}},
      {"structural_chop_digits", structural_chop_digits},
      {"acb_transformation", checkpoint_acb_laurent_matrix_record(
           *acb_transformation)}};
  auto& encoded_valuations =
      proof_without_identity.at("column_valuations").as_array();
  for (const auto valuation : column_valuations)
    encoded_valuations.push_back(valuation);
  const auto identity = json::serialize(
      canonical_json_value(proof_without_identity));
  auto proof = proof_without_identity;
  proof["identity"] = identity;
  return {kNativeSingularSCCSaturationCompactProofSchema, identity,
          json::serialize(canonical_json_value(proof)),
          std::move(*rational_saturation), false,
          std::move(*acb_transformation)};
}

ParsedExactEvaluatedLattice parse_native_singular_scc_saturation_proof(
    const json::value& raw, std::uint32_t dimension, EpsilonWindow window,
    const std::vector<json::object>& expected_basis_sources,
    const std::string& expected_basis_point,
    const std::string& expected_physical_point,
    const std::optional<std::string>& expected_session_configuration,
    const std::optional<json::object>& expected_native_request,
    const std::string& context) {
  const auto& proof = as_object(
      raw, "native Acb singular-SCC valuation-zero proof");
  const auto proof_schema = required_string(proof, "schema");
  const bool compact_proof =
      proof_schema == kNativeSingularSCCSaturationCompactProofSchema;
  if (!compact_proof &&
      proof_schema != kNativeSingularSCCSaturationProofSchema)
    throw std::invalid_argument(
        context +
        ": native singular-SCC saturation proof schema is unsupported");
  const bool shadow_proof =
      required_string(proof, "transformation") ==
      "exact-rational-shadow-column-monomials";
  const bool acb_monomial_proof =
      required_string(proof, "transformation") ==
      "certified-acb-column-monomials";
  const bool acb_levelt_proof =
      required_string(proof, "transformation") ==
      "certified-acb-Levelt-lattice-saturation";
  const bool direct_acb_proof =
      required_string(proof, "transformation") ==
      "direct-acb-Laurent-pivoting";
  const bool monomial_proof =
      shadow_proof || acb_monomial_proof || acb_levelt_proof;
  if (shadow_proof)
    require_exact_keys(
        proof,
        {"schema", "identity", "native_request", "coefficient_domain",
         "basis", "basis_point_exact", "physical_match_point_exact",
         "epsilon", "negative_epsilon_coefficients", "leading_power",
         "leading_rank", "leading_rank_certificate",
         "column_provenance_certificate", "determinant_valuation",
         "transformation", "column_valuations"},
        "native Acb monomial singular-SCC proof");
  else if (acb_monomial_proof || acb_levelt_proof)
    require_exact_keys(
        proof,
        {"schema", "identity", "native_request", "coefficient_domain",
         "basis", "basis_point_exact", "physical_match_point_exact",
         "epsilon", "negative_epsilon_coefficients", "leading_power",
         "leading_rank", "leading_rank_certificate",
         "column_provenance_certificate", "determinant_valuation",
         "transformation", "column_valuations", "structural_chop_digits",
         "acb_transformation"},
        "native Acb saturated singular-SCC proof");
  else
    require_exact_keys(
        proof,
        {"schema", "identity", "native_request", "coefficient_domain",
         "basis", "basis_point_exact", "physical_match_point_exact",
         "epsilon", "negative_epsilon_coefficients", "leading_power",
         "leading_rank", "leading_rank_certificate",
         "column_provenance_certificate", "determinant_valuation",
         "transformation"},
        "native Acb singular-SCC valuation-zero proof");
  const auto coefficient_domain =
      required_string(proof, "coefficient_domain");
  const auto leading_power = as_i32(
      proof.at("leading_power"), "singular-SCC proof leading power");
  const auto leading_rank = as_u32(
      proof.at("leading_rank"), "singular-SCC proof leading rank");
  const auto determinant_valuation = as_i32(
      proof.at("determinant_valuation"),
      "singular-SCC determinant valuation");
  if (coefficient_domain != "acb" || leading_power != 0 ||
      leading_rank != dimension ||
      (!shadow_proof && determinant_valuation != 0))
    throw std::invalid_argument(
        context +
        ": native singular-SCC saturation proof facts are inconsistent "
        "(coefficient_domain=" + coefficient_domain +
        ", leading_power=" + std::to_string(leading_power) +
        ", leading_rank=" + std::to_string(leading_rank) +
        ", expected_dimension=" + std::to_string(dimension) +
        ", determinant_valuation=" +
        std::to_string(determinant_valuation) + ")");
  if (shadow_proof) {
    if (required_string(proof, "negative_epsilon_coefficients") !=
            "exact-rational-shadow-formal-valuation" ||
        required_string(proof, "leading_rank_certificate") !=
            "deferred-acb-factorization-after-exact-column-shifts" ||
        required_string(proof, "column_provenance_certificate") !=
            "complete-target-bound-rational-shadow-case-p-floor-certified")
      throw std::invalid_argument(
          context + ": Rational-shadow saturation proof facts are inconsistent");
  } else if (acb_monomial_proof || acb_levelt_proof) {
    const auto expected_negative = acb_levelt_proof
        ? "acb-certified-finite-laurent-lattice"
        : "acb-certified-column-valuation";
    const auto expected_rank = acb_levelt_proof
        ? "full-acb-Levelt-saturation-pivots-exclude-zero"
        : "full-pivot-acb-pivots-exclude-zero";
    if (required_string(proof, "negative_epsilon_coefficients") !=
            expected_negative ||
        required_string(proof, "leading_rank_certificate") !=
            expected_rank ||
        required_string(proof, "column_provenance_certificate") !=
            "complete-exact-schedule-receiving-scc-acb" ||
        as_i32(proof.at("structural_chop_digits"),
               "singular-SCC structural chop digits") < 32)
      throw std::invalid_argument(
          context + ": native singular-SCC Acb monomial proof facts are inconsistent");
  } else if (direct_acb_proof) {
    if (required_string(proof, "negative_epsilon_coefficients") !=
            "retained-certified-Acb-Laurent-frames" ||
        required_string(proof, "leading_rank_certificate") !=
            "full-pivot-acb-Laurent-pivots-exclude-zero" ||
        required_string(proof, "column_provenance_certificate") !=
            "complete-exact-schedule-receiving-scc-acb")
      throw std::invalid_argument(
          context + ": native singular-SCC direct Acb Laurent proof facts are inconsistent");
  } else if (required_string(proof, "negative_epsilon_coefficients") !=
                 "exact-singleton-zero" ||
             required_string(proof, "leading_rank_certificate") !=
                 "full-pivot-acb-pivots-exclude-zero" ||
             required_string(proof, "column_provenance_certificate") !=
                 "complete-one-receiving-scc-composite-affine-jordan-acb-no-pseudo" ||
             required_string(proof, "transformation") != "identity")
    throw std::invalid_argument(
        context + ": native singular-SCC identity proof facts are inconsistent");
  auto native_request = validate_native_singular_scc_saturation_request(
      proof.at("native_request"), context,
      expected_session_configuration,
      expected_native_request);
  const auto& epsilon = as_object(
      proof.at("epsilon"), "singular-SCC proof epsilon window");
  require_exact_keys(epsilon, {"min", "max"},
                     "singular-SCC proof epsilon window");
  if (as_i32(epsilon.at("min"), "singular-SCC proof epsilon minimum") !=
          window.min_power ||
      as_i32(epsilon.at("max"), "singular-SCC proof epsilon maximum") !=
          window.complete_max ||
      required_string(proof, "basis_point_exact") !=
          expected_basis_point ||
      required_string(proof, "physical_match_point_exact") !=
          expected_physical_point ||
      required_string(native_request, "receiving_basis_point_exact") !=
          expected_basis_point ||
      required_string(native_request, "physical_match_point_exact") !=
          expected_physical_point)
    throw std::invalid_argument(
        context + ": native singular-SCC proof changed its point or epsilon binding");
  json::array expected_basis;
  if (compact_proof) {
    expected_basis =
        compact_matching_source_references(expected_basis_sources);
  } else {
    expected_basis.reserve(expected_basis_sources.size());
    for (const auto& source : expected_basis_sources)
      expected_basis.push_back(checkpoint_matching_source_record(source));
  }
  if (json::serialize(canonical_json_value(proof.at("basis"))) !=
      json::serialize(canonical_json_value(expected_basis)))
    throw std::invalid_argument(
        context + ": native singular-SCC proof changed its basis/checkpoint binding");
  const bool bounded_request =
      required_string(native_request, "schema") ==
      kNativeSingularSCCSaturationBoundedRequestSchema;
  if (bounded_request && !expected_native_request.has_value())
    throw std::invalid_argument(
        context +
        ": bounded singular-SCC proof lacks its retained exact request owner");
  validate_singular_scc_basis_sources(
      bounded_request ? *expected_native_request : native_request,
      expected_basis_sources, dimension, context);
  auto identity_input = proof;
  const auto identity = required_string(proof, "identity");
  identity_input.erase("identity");
  if (identity.empty() ||
      json::serialize(canonical_json_value(identity_input)) != identity)
    throw std::invalid_argument(
        context + ": native singular-SCC saturation proof identity is inconsistent");
  if (!monomial_proof)
    return {proof_schema, identity,
            json::serialize(canonical_json_value(proof)),
            unit_rational_saturation(dimension, window, context)};
  const auto& raw_valuations = as_array(
      proof.at("column_valuations"),
      "singular-SCC monomial saturation column valuations");
  if (raw_valuations.size() != dimension)
    throw std::invalid_argument(
        context + ": monomial saturation valuation count differs from its basis dimension");
  FiniteLaurentMatrix<Rational> monomial_basis(
      dimension, FiniteLaurentVector<Rational>());
  for (std::uint32_t row = 0; row < dimension; ++row) {
    monomial_basis[row].reserve(dimension);
    for (std::uint32_t column = 0; column < dimension; ++column) {
      const auto valuation = as_i32(
          raw_valuations[column],
          "singular-SCC monomial saturation column valuation");
      if (valuation < window.min_power || valuation > window.complete_max)
        throw std::invalid_argument(
            context + ": restored monomial valuation lies outside its epsilon window");
      std::vector<Rational> coefficients(window.width(), Rational(0));
      if (row == column)
        coefficients[static_cast<std::size_t>(
            static_cast<std::int64_t>(valuation) - window.min_power)] =
            Rational(1);
      monomial_basis[row].emplace_back(window, std::move(coefficients));
    }
  }
  auto saturation = saturate_finite_laurent_basis(
      monomial_basis, context + ":restored-singular-SCC-valuations");
  std::optional<ExactLaurentMatrix<ComplexBall>> acb_transformation;
  if (acb_monomial_proof || acb_levelt_proof)
    acb_transformation = parse_checkpoint_acb_laurent_matrix(
        proof.at("acb_transformation"), dimension,
        "singular-SCC Acb saturation transformation");
  return {proof_schema, identity,
          json::serialize(canonical_json_value(proof)),
          std::move(saturation), shadow_proof,
          std::move(acb_transformation)};
}

std::optional<std::int32_t> parse_optional_match_imaginary_sign(
    const json::object& request, const char* key) {
  const auto* raw = request.if_contains(key);
  if (raw == nullptr || raw->is_null()) return std::nullopt;
  const auto sign = as_i32(*raw, key);
  if (sign != -1 && sign != 1)
    throw std::invalid_argument(std::string(key) + " must be +1 or -1");
  return sign;
}

void require_acb_match_local(const LocalSolution<ComplexBall>& solution,
                             const RealEvaluationPoint& point,
                             const std::string& label) {
  validate_local_solution(solution, true);
  if (!solution.error.empty())
    throw std::invalid_argument(
        label +
        " carries an error envelope not represented in Acb matching residuals");
  if (!solution.chart.infinite_radius &&
      !arb_lt(acb_realref(point.modulus.raw()),
              acb_realref(solution.chart.radius.raw())))
    throw std::invalid_argument(
        label + " match point is not provably inside its chart radius");
}

RealEvaluationPoint parse_acb_match_point(const json::value& raw,
                                          const char* label) {
  const auto& point = as_object(raw, label);
  if (point.if_contains("value") == nullptr &&
      point.if_contains("sign") == nullptr)
    return RealEvaluationPoint::rational(required_string(point, "exact"));
  require_exact_keys(point, {"exact", "value", "sign"}, label);
  return RealEvaluationPoint::certified(
      required_string(point, "exact"),
      parse_scalar<ComplexBall>(point.at("value")),
      as_i32(point.at("sign"), label));
}

Rational acb_physical_match_point(const ChartGeometry& chart,
                                  const RealEvaluationPoint& local_point,
                                  const std::string& label) {
  try {
    return Rational(chart.center_exact) +
           Rational(chart.scale_exact) *
               Rational(local_point.exact_coordinate);
  } catch (const std::invalid_argument&) {
    throw std::invalid_argument(
        label +
        " requires exact rational chart center/scale geometry for a common physical point proof");
  }
}

FiniteLaurentVector<ComplexBall> acb_evaluation_frames(
    const EpsilonVector& value, EpsilonWindow window,
    const std::string& label) {
  if (value.epsilon.complete_max < window.complete_max)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InsufficientCompleteWindow,
        label + " does not cover the requested matching epsilon window",
        std::nullopt, std::nullopt, value.epsilon.complete_max);
  if (window.min_power > value.epsilon.min_power) {
    for (std::int64_t power = value.epsilon.min_power;
         power < window.min_power; ++power)
      for (std::uint32_t component = 0; component < value.dimension;
           ++component)
        if (!value.at(static_cast<std::int32_t>(power), component).is_zero())
          throw MatchingArithmeticError(
              MatchingArithmeticErrorCode::InsufficientCompleteWindow,
              label +
                  " work minimum would discard a nonzero or zero-ambiguous lower epsilon coefficient",
              component, std::nullopt,
              static_cast<std::int32_t>(power));
  }
  FiniteLaurentVector<ComplexBall> frames;
  frames.reserve(value.dimension);
  for (std::uint32_t component = 0; component < value.dimension;
       ++component) {
    std::vector<ComplexBall> coefficients;
    coefficients.reserve(window.width());
    for (std::int64_t power = window.min_power;
         power <= window.complete_max; ++power) {
      const auto epsilon_power = static_cast<std::int32_t>(power);
      coefficients.push_back(epsilon_power < value.epsilon.min_power
          ? ComplexBall(0)
          : value.at(epsilon_power, component));
    }
    frames.emplace_back(window, std::move(coefficients));
  }
  return frames;
}

json::value optional_match_sign_json(
    const std::optional<std::int32_t>& sign) {
  return sign.has_value() ? json::value(*sign) : json::value(nullptr);
}

std::shared_ptr<StoredRefinedAcbMatch> build_refined_acb_match_once(
    const std::string& match_handle, const json::object& request,
    const std::vector<std::string>& basis_handles,
    const std::vector<std::shared_ptr<StoredLocalBase>>& erased_basis,
    const std::string& incoming_handle,
    const std::shared_ptr<StoredLocalBase>& erased_incoming,
    slong precision_bits,
    const std::string& active_session_configuration_identity,
    const std::optional<json::object>& expected_singular_request,
    bool require_normalized_singular_frame,
    const std::optional<std::size_t>& common_taylor_width,
    bool allow_exact_recurrence_diagnostic) {
  const auto started = std::chrono::steady_clock::now();
  if (common_taylor_width.has_value() && *common_taylor_width == 0)
    throw std::invalid_argument(
        "Acb matching Taylor prefix must retain at least one coefficient");
  if (request.if_contains("native_singular_scc_saturation") != nullptr &&
      !expected_singular_request.has_value())
    throw std::invalid_argument(
        "native singular-SCC Acb saturation is admitted only through a retained planned match");
  AcbPrecisionLease lease(precision_bits);
  ComplexBall::set_precision(precision_bits);

  const auto& raw_window = as_object(
      request.at("epsilon"), "Acb local match epsilon window");
  EpsilonWindow window{as_i32(raw_window.at("min"), "match epsilon minimum"),
                       as_i32(raw_window.at("max"), "match epsilon maximum")};
  (void)window.width();
  const auto required_complete_max = as_i32(
      raw_window.at("required_complete_max"),
      "required Acb match residual complete maximum");
  if (required_complete_max < window.min_power ||
      required_complete_max > window.complete_max)
    throw std::invalid_argument(
        "required Acb match residual maximum must lie inside the supplied work epsilon window");

  const auto& raw_refinement = as_object(
      request.at("refinement"), "Acb local match refinement policy");
  if (raw_refinement.size() != 2 ||
      raw_refinement.if_contains("relative_tolerance") == nullptr ||
      raw_refinement.if_contains("max_steps") == nullptr)
    throw std::invalid_argument(
        "Acb match refinement accepts exactly relative_tolerance and max_steps");
  const auto relative_tolerance = required_string(
      raw_refinement, "relative_tolerance");
  const auto max_refinement_steps = as_u32(
      raw_refinement.at("max_steps"), "Acb match refinement steps");
  if (max_refinement_steps > 32)
    throw std::invalid_argument(
        "Acb match refinement steps must lie in 0..32");
  AcbLaurentRefinementOptions refinement;
  refinement.relative_tolerance = Magnitude::decimal(relative_tolerance);
  refinement.required_complete_max = required_complete_max;
  refinement.max_refinement_steps = max_refinement_steps;

  const auto basis_point = parse_acb_match_point(
      request.at("basis_point"), "basis match point");
  const auto incoming_point = parse_acb_match_point(
      request.at("incoming_point"), "incoming match point");
  const auto requested_basis_sign = parse_optional_match_imaginary_sign(
      request, "basis_imaginary_sign");
  const auto requested_incoming_sign = parse_optional_match_imaginary_sign(
      request, "incoming_imaginary_sign");
  EvaluationOptions basis_options;
  basis_options.imaginary_sign = requested_basis_sign;
  basis_options.compute_tail_estimate = false;
  EvaluationOptions incoming_options;
  incoming_options.imaginary_sign = requested_incoming_sign;
  incoming_options.compute_tail_estimate = false;

  const auto basis_chart = required_string(request, "basis_chart");
  const auto incoming_chart = required_string(request, "incoming_chart");
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  if (checkpoint_identity.empty())
    throw std::invalid_argument(
        "refined Acb match checkpoint identity cannot be empty");

  std::vector<std::shared_ptr<StoredLocal<ComplexBall>>> basis;
  basis.reserve(erased_basis.size());
  for (const auto& local : erased_basis) {
    auto typed = std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(local);
    if (!typed)
      throw std::invalid_argument(
          "refined Acb matching requires Acb retained locals");
    basis.push_back(std::move(typed));
  }
  auto incoming =
      std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(erased_incoming);
  if (!incoming)
    throw std::invalid_argument(
        "refined Acb matching requires an Acb incoming local");
  if (basis.empty())
    throw std::invalid_argument(
        "refined Acb matching requires a nonempty basis");
  const auto dimension = basis.front()->solution().dimension;
  if (basis.size() != dimension || incoming->solution().dimension != dimension)
    throw std::invalid_argument(
        "refined Acb matching requires d basis columns and a d-component incoming local");

  const auto& raw_basis_checkpoints = as_array(
      request.at("basis_checkpoint_identities"),
      "Acb basis checkpoint identities");
  if (raw_basis_checkpoints.size() != basis.size())
    throw std::invalid_argument(
        "Acb basis checkpoint identity count differs from the basis dimension");
  std::vector<std::string> basis_checkpoints;
  basis_checkpoints.reserve(basis.size());
  std::vector<json::object> basis_sources;
  basis_sources.reserve(basis.size());
  for (std::size_t column = 0; column < basis.size(); ++column) {
    if (!raw_basis_checkpoints[column].is_string())
      throw std::invalid_argument(
          "Acb basis checkpoint identities must be strings");
    const std::string expected(raw_basis_checkpoints[column].as_string());
    if (expected.empty() ||
        basis[column]->solution().checkpoint_identity != expected)
      throw std::invalid_argument(
          "Acb basis checkpoint provenance mismatch at column " +
          std::to_string(column));
    if (basis[column]->source_chart() != basis_chart)
      throw std::invalid_argument(
          "Acb basis chart provenance mismatch at column " +
          std::to_string(column));
    if (!same_chart_geometry(basis.front()->solution().chart,
                             basis[column]->solution().chart))
      throw std::invalid_argument(
          "Acb basis locals do not share identical retained chart geometry");
    require_acb_match_local(
        basis[column]->solution(), basis_point,
        "Acb basis local " + basis_handles[column]);
    basis_checkpoints.push_back(expected);
  }
  const auto expected_incoming_checkpoint = required_string(
      request, "incoming_checkpoint_identity");
  if (expected_incoming_checkpoint.empty() ||
      incoming->solution().checkpoint_identity !=
          expected_incoming_checkpoint)
    throw std::invalid_argument(
        "Acb incoming checkpoint provenance mismatch");
  if (incoming->source_chart() != incoming_chart)
    throw std::invalid_argument("Acb incoming chart provenance mismatch");
  require_acb_match_local(incoming->solution(), incoming_point,
                          "Acb incoming local " + incoming_handle);

  std::string basis_physical_point;
  if (const auto* certified = request.if_contains(
          "certified_physical_match_point_exact")) {
    if (!certified->is_string() || certified->as_string().empty() ||
        !basis_point.certified_algebraic ||
        !incoming_point.certified_algebraic)
      throw std::invalid_argument(
          "certified physical match identity requires two certified algebraic local points");
    basis_physical_point = std::string(certified->as_string());
  } else {
    const auto exact_basis_physical = acb_physical_match_point(
        basis.front()->solution().chart, basis_point,
        "Acb basis match point");
    for (std::size_t column = 1; column < basis.size(); ++column)
      if (!(acb_physical_match_point(
                basis[column]->solution().chart, basis_point,
                "Acb basis match point at column " +
                    std::to_string(column)) == exact_basis_physical))
        throw std::invalid_argument(
            "Acb basis locals do not name one exact physical match point");
    const auto incoming_physical_point = acb_physical_match_point(
        incoming->solution().chart, incoming_point,
        "Acb incoming match point");
    if (!(exact_basis_physical == incoming_physical_point))
      throw std::invalid_argument(
          "Acb basis and incoming coordinates do not name the same exact physical match point");
    basis_physical_point = exact_basis_physical.str();
  }

  FiniteLaurentMatrix<ComplexBall> evaluated_basis(
      dimension, FiniteLaurentVector<ComplexBall>());
  for (auto& row : evaluated_basis) row.reserve(dimension);
  std::vector<std::optional<std::int32_t>> effective_basis_signs;
  effective_basis_signs.reserve(dimension);
  std::vector<std::size_t> basis_taylor_widths;
  basis_taylor_widths.reserve(dimension);
  std::vector<LocalEvaluation> basis_evaluations;
  basis_evaluations.reserve(dimension);
  for (std::size_t column = 0; column < basis.size(); ++column) {
    auto column_options = basis_options;
    const auto full_width = basis[column]->solution().taylor_width();
    const auto retained_width = common_taylor_width.value_or(full_width);
    if (retained_width > full_width)
      throw std::invalid_argument(
          "Acb matching common Taylor prefix exceeds a basis local");
    column_options.t_order_reduction = static_cast<std::uint32_t>(
        full_width - retained_width);
    basis_evaluations.push_back(evaluate_local_solution(
        basis[column]->solution(), basis_point, column_options));
    basis_taylor_widths.push_back(retained_width);
    effective_basis_signs.push_back(
        basis_evaluations.back().imaginary_sign);
  }
  const auto incoming_full_width = incoming->solution().taylor_width();
  const auto incoming_taylor_width =
      common_taylor_width.value_or(incoming_full_width);
  if (incoming_taylor_width > incoming_full_width)
    throw std::invalid_argument(
        "Acb matching common Taylor prefix exceeds the incoming local");
  incoming_options.t_order_reduction = static_cast<std::uint32_t>(
      incoming_full_width - incoming_taylor_width);
  const auto incoming_evaluation = evaluate_local_solution(
      incoming->solution(), incoming_point, incoming_options);
  auto evaluation_window = window;
  if (expected_singular_request.has_value()) {
    auto common_basis_max = basis_evaluations.front().value.epsilon.complete_max;
    for (const auto& evaluation : basis_evaluations)
      common_basis_max = std::min(
          common_basis_max, evaluation.value.epsilon.complete_max);
    common_basis_max = std::min(
        common_basis_max, incoming_evaluation.value.epsilon.complete_max);
    const auto desired_basis_max = matching_detail::checked_power(
        static_cast<std::int64_t>(window.complete_max) +
            2 * static_cast<std::int64_t>(dimension),
        "singular matching basis halo maximum");
    evaluation_window.complete_max = std::max(
        window.complete_max, std::min(common_basis_max, desired_basis_max));
  }
  for (std::size_t column = 0; column < basis.size(); ++column) {
    auto frames = acb_evaluation_frames(
        basis_evaluations[column].value, evaluation_window,
        "Acb basis evaluation at column " + std::to_string(column));
    for (std::uint32_t component = 0; component < dimension; ++component)
      evaluated_basis[component].push_back(std::move(frames[component]));
  }
  auto incoming_value = acb_evaluation_frames(
      incoming_evaluation.value, evaluation_window,
      "Acb incoming evaluation");

  // In an epsilon-regular singular SCC, the physical fundamental columns can
  // be strongly confluent even though the reduced system is Fuchsian.  Match
  // after applying the receiving SCC's exact V^-1 normal frame to both sides.
  // This is an invertible left transformation, so it cannot change the
  // weights; it only removes the artificial finite-precision conditioning.
  auto matching_basis = evaluated_basis;
  auto matching_incoming = incoming_value;
  auto matching_window = evaluation_window;
  std::string matching_frame_identity = "physical-parent-frame";
  bool normalized_matching_frame = false;
  std::optional<ExactLaurentMatrix<ComplexBall>>
      normalized_matching_right_inverse;
  std::optional<FiniteLaurentMatrix<ComplexBall>>
      left_normalized_matching_basis;
  std::optional<FiniteLaurentMatrix<ComplexBall>>
      matching_right_normalization_matrix;
  json::object normal_frame_attempt{
      {"schema", "diffexp2-acb-normal-frame-attempt-v1"},
      {"applicable", expected_singular_request.has_value()},
      {"status", expected_singular_request.has_value()
           ? "receiving-owner-unavailable" : "not-applicable"},
      {"physical_window",
       json::object{{"min", evaluation_window.min_power},
                    {"max", evaluation_window.complete_max}}},
      {"required_complete_max", required_complete_max}};
  if (expected_singular_request.has_value()) {
    const auto receiving_owner = basis.front()->retained_equation_owner();
    if (receiving_owner) {
      normal_frame_attempt["status"] = "normalizing-basis";
      for (std::size_t column = 1; column < basis.size(); ++column)
        if (basis[column]->retained_equation_owner().get() !=
            receiving_owner.get())
          throw std::invalid_argument(
              "singular matching basis columns do not share one receiving SCC owner");

      FiniteLaurentMatrix<ComplexBall> candidate(
          dimension, FiniteLaurentVector<ComplexBall>());
      for (auto& row : candidate) row.reserve(dimension);
      std::optional<std::string> candidate_identity;
      bool available = true;
      for (std::size_t column = 0; column < dimension; ++column) {
        FiniteLaurentVector<ComplexBall> physical_column;
        physical_column.reserve(dimension);
        for (std::size_t row = 0; row < dimension; ++row)
          physical_column.push_back(evaluated_basis[row][column]);
        auto transformed = receiving_owner->normalize_acb_matching_vector(
            physical_column, basis_point);
        if (!transformed.has_value()) {
          available = false;
          normal_frame_attempt["status"] =
              "basis-normalizer-unavailable";
          normal_frame_attempt["basis_column"] = column;
          normal_frame_attempt["owner_diagnostic"] =
              receiving_owner
                  ->acb_matching_normal_frame_diagnostic(
                      basis_point);
          break;
        }
        if (candidate_identity.has_value() &&
            *candidate_identity != transformed->second)
          throw std::logic_error(
              "receiving SCC matching normal-frame identity changed between columns");
        candidate_identity = transformed->second;
        for (std::size_t row = 0; row < dimension; ++row)
          candidate[row].push_back(std::move(transformed->first[row]));
      }
      if (!available && require_normalized_singular_frame)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InvalidSaturationLattice,
            checkpoint_identity +
                ": terminal singular match requires an available Fuchsian "
                "normal frame; owner_diagnostic=" +
                json::serialize(
                    normal_frame_attempt.at("owner_diagnostic")));
      if (available) {
        auto transformed_incoming =
            receiving_owner->normalize_acb_matching_vector(
                incoming_value, basis_point);
        if (!transformed_incoming.has_value() ||
            !candidate_identity.has_value() ||
            transformed_incoming->second != *candidate_identity)
          throw std::logic_error(
              "incoming value lost the receiving SCC matching normal frame");
        auto two_sided_basis =
            receiving_owner->right_normalize_acb_matching_basis(candidate);
        if (!two_sided_basis.has_value())
          throw std::logic_error(
              "receiving SCC supplied V^-1 but not the matching right V frame");
        normalized_matching_right_inverse =
            receiving_owner
                ->right_normalization_acb_matching_transformation();
        if (!normalized_matching_right_inverse.has_value())
          throw std::logic_error(
              "receiving SCC supplied its numeric right V frame but not the retained exact-support V^-1 transformation");
        auto retained_right = receiving_owner
            ->right_denormalization_acb_matching_matrix();
        if (!retained_right.has_value())
          throw std::logic_error(
              "receiving SCC supplied its right V frame without retaining its finite matching matrix");
        left_normalized_matching_basis = candidate;
        matching_right_normalization_matrix =
            std::move(*retained_right);
        matching_basis = std::move(*two_sided_basis);
        matching_incoming = std::move(transformed_incoming->first);
        matching_frame_identity = std::move(*candidate_identity);
        normalized_matching_frame = true;
        normal_frame_attempt["candidate_frame"] =
            compact_matching_identity_reference(
                matching_frame_identity);

        std::optional<std::int32_t> common_min;
        std::optional<std::int32_t> common_max;
        const auto include_frame = [&](const auto& frame) {
          common_min = !common_min.has_value()
              ? frame.min_power()
              : std::min(*common_min, frame.min_power());
          common_max = !common_max.has_value()
              ? frame.complete_max()
              : std::min(*common_max, frame.complete_max());
        };
        for (const auto& row : matching_basis)
          for (const auto& frame : row) include_frame(frame);
        for (const auto& frame : matching_incoming) include_frame(frame);
        normal_frame_attempt["candidate_window"] =
            common_min.has_value() && common_max.has_value()
            ? json::value(json::object{{"min", *common_min},
                                      {"max", *common_max}})
            : json::value(nullptr);
        if (!common_min.has_value() || !common_max.has_value() ||
            *common_min > *common_max ||
            *common_max < required_complete_max) {
          if (require_normalized_singular_frame) {
            if (!common_min.has_value() || !common_max.has_value() ||
                *common_min > *common_max)
              throw MatchingArithmeticError(
                  MatchingArithmeticErrorCode::InvalidSaturationLattice,
                  checkpoint_identity +
                      ": terminal singular Fuchsian frame has no common "
                      "finite epsilon window");
            throw MatchingArithmeticError(
                MatchingArithmeticErrorCode::InsufficientCompleteWindow,
                checkpoint_identity +
                    ": terminal singular Fuchsian frame does not reach the required epsilon prefix",
                std::nullopt, std::nullopt, *common_max);
          }
          // A low public-order request need not recursively manufacture the
          // private Laurent halo consumed by V^-1/V.  In that case retain the
          // physical frame: its residual certificate is still authoritative
          // through the requested prefix, and no surplus-order condition is
          // imposed merely to enable the conditioning optimization.
          matching_basis = evaluated_basis;
          matching_incoming = incoming_value;
          matching_window = evaluation_window;
          matching_frame_identity = "physical-parent-frame";
          normalized_matching_frame = false;
          normalized_matching_right_inverse.reset();
          normal_frame_attempt["status"] =
              "transformed-window-insufficient";
        } else {
          matching_window = {*common_min, *common_max};
          normal_frame_attempt["status"] =
              "normalized-refinement-started";
        }
      }
    } else if (require_normalized_singular_frame) {
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::InvalidSaturationLattice,
          checkpoint_identity +
              ": terminal singular match has no retained receiving equation "
              "owner for its required Fuchsian normal frame");
    }
  }
  if (const auto* diagnostics =
          std::getenv("DE2_DIAGNOSTIC_TERMINAL_STATE");
      diagnostics != nullptr && std::string(diagnostics) == "1" &&
      expected_singular_request.has_value()) {
    const auto frame_scales = [](const auto& frames) {
      std::map<std::int32_t, ComplexBall> largest;
      const auto include = [&](const EpsilonFrame<ComplexBall>& frame) {
        for (std::int64_t raw_power = frame.min_power();
             raw_power <= frame.complete_max(); ++raw_power) {
          const auto power = matching_detail::checked_power(
              raw_power, "matching-frame scale diagnostic power");
          const auto& coefficient = frame.coefficient(power);
          const auto magnitude = Magnitude::upper_abs(coefficient);
          if (magnitude.is_zero()) continue;
          const auto found = largest.find(power);
          if (found == largest.end() ||
              !(magnitude <=
                Magnitude::upper_abs(found->second)))
            largest.insert_or_assign(power, coefficient);
        }
      };
      for (const auto& row : frames)
        if constexpr (std::is_same_v<
                          std::decay_t<decltype(row)>,
                          EpsilonFrame<ComplexBall>>)
          include(row);
        else
          for (const auto& frame : row) include(frame);
      json::array encoded;
      encoded.reserve(largest.size());
      for (const auto& [power, coefficient] : largest)
        encoded.push_back(json::object{
            {"power", power},
            {"absolute_upper_exact",
             Magnitude::upper_abs(coefficient).dump_exact()},
            {"real_radius_exponent",
             coefficient.real_radius_exponent()},
            {"imag_radius_exponent",
             coefficient.imag_radius_exponent()}});
      return encoded;
    };
    normal_frame_attempt["physical_basis_scales"] =
        frame_scales(evaluated_basis);
    normal_frame_attempt["physical_incoming_scales"] =
        frame_scales(incoming_value);
    normal_frame_attempt["matching_basis_scales"] =
        frame_scales(matching_basis);
    normal_frame_attempt["matching_incoming_scales"] =
        frame_scales(matching_incoming);
  }
  refinement.required_min_power = matching_window.min_power;
  std::string residual_certificate_identity;
  std::optional<FiniteLaurentMatrix<ComplexBall>>
      selected_terminal_normal_frame_right_transformation;

  const auto reset_to_physical_matching_frame = [&]() {
    matching_basis = evaluated_basis;
    matching_incoming = incoming_value;
    matching_window = evaluation_window;
    matching_frame_identity = "physical-parent-frame";
    normalized_matching_frame = false;
    refinement.required_min_power = matching_window.min_power;
    residual_certificate_identity.clear();
    selected_terminal_normal_frame_right_transformation.reset();
    normalized_matching_right_inverse.reset();
    left_normalized_matching_basis.reset();
    matching_right_normalization_matrix.reset();
  };
  if (require_normalized_singular_frame &&
      expected_singular_request.has_value() &&
      !normalized_matching_frame)
    throw std::logic_error(
        checkpoint_identity +
        ": strict terminal singular match escaped normal-frame selection");

  for (std::size_t column = 0; column < basis.size(); ++column) {
    json::object source{
        {"column", column}, {"local", basis_handles[column]},
        {"chart", basis_chart},
        {"source_operator_identity",
         basis[column]->source_operator_identity()},
        {"source_operator_reference",
         basis[column]->source_operator_reference()},
        {"checkpoint_identity", basis_checkpoints[column]},
        {"requested_imaginary_sign",
         optional_match_sign_json(requested_basis_sign)},
        {"effective_imaginary_sign",
         optional_match_sign_json(effective_basis_signs[column])},
        {"matching_taylor_width", basis_taylor_widths[column]},
        {"analytic_metadata", basis[column]->exact_analytic_metadata()}};
    if (basis[column]->column_provenance().has_value())
      source["column_provenance"] =
          basis[column]->column_provenance()->encode();
    basis_sources.push_back(std::move(source));
  }
  json::object incoming_source{
      {"local", incoming_handle}, {"chart", incoming_chart},
      {"source_operator_identity",
       incoming->source_operator_identity()},
      {"source_operator_reference",
       incoming->source_operator_reference()},
      {"checkpoint_identity", expected_incoming_checkpoint},
      {"requested_imaginary_sign",
       optional_match_sign_json(requested_incoming_sign)},
      {"effective_imaginary_sign",
       optional_match_sign_json(incoming_evaluation.imaginary_sign)},
      {"matching_taylor_width", incoming_taylor_width},
      {"analytic_metadata", incoming->exact_analytic_metadata()}};
  if (incoming->column_provenance().has_value())
    incoming_source["column_provenance"] =
        incoming->column_provenance()->encode();

  const auto make_exact_lattice = [&]() -> ParsedExactEvaluatedLattice {
    const auto proof_request_count =
        (request.if_contains("exact_lattice") != nullptr ? 1U : 0U) +
        (request.if_contains("native_unit_saturation") != nullptr ? 1U : 0U) +
        (request.if_contains("native_singular_scc_saturation") != nullptr
             ? 1U : 0U);
    if (proof_request_count != 1)
      throw std::invalid_argument(
          "Acb matching requires exactly one exact lattice, ordinary native unit-leading request, or singular-SCC valuation-zero request");
    if (const auto* raw_exact = request.if_contains("exact_lattice")) {
      return parse_exact_evaluated_lattice(
          *raw_exact, dimension, matching_window, checkpoint_identity);
    }
    if (const auto* raw_native =
            request.if_contains("native_unit_saturation"))
      return certify_native_unit_saturation(
          *raw_native, matching_basis, erased_basis, basis_sources,
          basis_point.exact_coordinate, basis_physical_point,
          matching_window,
          checkpoint_identity + ":native-unit-leading-proof");
    return certify_native_singular_scc_saturation(
        request.at("native_singular_scc_saturation"), matching_basis,
        erased_basis, basis_sources, basis_point.exact_coordinate,
        basis_physical_point, matching_window,
        active_session_configuration_identity,
        *expected_singular_request,
        checkpoint_identity,
        checkpoint_identity + ":native-singular-scc-valuation-zero-proof",
        true, normalized_matching_frame
                  ? normalized_matching_right_inverse
                  : std::nullopt);
  };
  auto exact_lattice = make_exact_lattice();
  if (const auto* diagnostic =
          std::getenv("DE2_DIAGNOSTIC_RATIONAL_SHADOW_QC");
      diagnostic != nullptr && std::string(diagnostic) == "1" &&
      allow_exact_recurrence_diagnostic &&
      require_normalized_singular_frame &&
      exact_lattice.exact_shadow_fused_local_basis.has_value()) {
    const auto& owner = exact_lattice.exact_shadow_equation_owner;
    if (!owner ||
        std::string(owner->equation_scalar_domain()) != "rational")
      throw std::logic_error(
          "exact Rational-shadow q/C diagnostic lost its Rational equation owner");
    const auto equation =
        std::static_pointer_cast<
            const PreparedPhysicalClearedODE<Rational>>(
            owner->physical_ode_erased());
    if (!equation)
      throw std::logic_error(
          "exact Rational-shadow q/C diagnostic lost its physical equation");
    auto diagnostic_taylor_width =
        basis_taylor_widths.front();
    if (const auto* raw_width = std::getenv(
            "DE2_DIAGNOSTIC_RATIONAL_SHADOW_QC_TAYLOR_WIDTH");
        raw_width != nullptr && std::string(raw_width).size() != 0) {
      std::size_t parsed = 0;
      const auto requested = std::stoull(raw_width, &parsed);
      if (parsed != std::string(raw_width).size() ||
          requested == 0 ||
          requested > diagnostic_taylor_width)
        throw std::invalid_argument(
            "DE2_DIAGNOSTIC_RATIONAL_SHADOW_QC_TAYLOR_WIDTH must lie inside the retained matching Taylor width");
      diagnostic_taylor_width =
          static_cast<std::size_t>(requested);
    }
    const EpsilonWindow requested_diagnostic_window{
        evaluation_window.min_power, required_complete_max};
    json::array columns;
    columns.reserve(
        exact_lattice.exact_shadow_fused_local_basis->size());
    for (std::size_t column = 0;
         column <
         exact_lattice.exact_shadow_fused_local_basis->size();
         ++column) {
      const auto& exact_local =
          (*exact_lattice.exact_shadow_fused_local_basis)[column];
      if (exact_local.epsilon.complete_max < required_complete_max)
        throw std::logic_error(
            "exact Rational-shadow q/C diagnostic column " +
            std::to_string(column) +
            " does not reach the required epsilon maximum: retained=[" +
            std::to_string(exact_local.epsilon.min_power) + "," +
            std::to_string(exact_local.epsilon.complete_max) +
            "], required_complete_max=" +
            std::to_string(required_complete_max));
      if (exact_local.taylor_complete_max + 1 <
          diagnostic_taylor_width)
        throw std::logic_error(
            "exact Rational-shadow q/C diagnostic column " +
            std::to_string(column) +
            " does not reach the requested Taylor width: retained_width=" +
            std::to_string(
                static_cast<std::uint64_t>(
                    exact_local.taylor_complete_max) + 1) +
            ", requested_width=" +
            std::to_string(diagnostic_taylor_width));
      // The saturated exact column can have a higher structural epsilon
      // minimum than the unsaturated Acb evaluation frame.  Rows below that
      // exact minimum are identically absent, not an uncomputed prefix.  Audit
      // the nonzero retained slab while continuing to require the complete
      // public upper edge.
      const EpsilonWindow diagnostic_window{
          std::max(requested_diagnostic_window.min_power,
                   exact_local.epsilon.min_power),
          requested_diagnostic_window.complete_max};
      const auto certificate =
          certify_scc_parent_exact_formal_residual(
              *equation, exact_local, diagnostic_window,
              static_cast<std::uint32_t>(
                  diagnostic_taylor_width - 1));
      columns.push_back(json::object{
          {"column", column},
          {"retained_epsilon",
           json::object{
               {"min", exact_local.epsilon.min_power},
               {"complete_max",
                exact_local.epsilon.complete_max}}},
          {"retained_taylor_complete_max",
           exact_local.taylor_complete_max},
          {"epsilon",
           json::object{
               {"min", certificate.epsilon.min_power},
               {"complete_max",
                certificate.epsilon.complete_max}}},
          {"taylor_complete_max",
           certificate.taylor_complete_max},
          {"exact_tag_count", certificate.exact_tag_count},
          {"coefficient_rows", certificate.coefficient_rows}});
    }
    std::fprintf(
        stderr,
        "[diffexp2 exact Rational-shadow q/C diagnostic] %s\n",
        json::serialize(json::object{
            {"schema",
             "diffexp2-exact-rational-shadow-qc-diagnostic-v1"},
            {"physical_point_exact", basis_physical_point},
            {"equation_payload_reference",
             json::object{
                 {"algorithm", "fnv1a64-v1"},
                 {"fingerprint",
                  public_provenance_fingerprint(
                      equation->payload_identity)},
                 {"identity_bytes",
                  equation->payload_identity.size()}}},
            {"columns", std::move(columns)}}).c_str());
  }
  constexpr std::size_t kExactShadowExtraPrecisionBits = 2048;
  std::optional<matching_detail::ScopedAcbPrecision>
      exact_shadow_precision;
  std::optional<FiniteLaurentMatrix<ComplexBall>>
      exact_shadow_physical_transformed_basis;
  std::optional<FiniteLaurentMatrix<ComplexBall>>
      exact_shadow_normalized_transformed_basis;
  std::optional<std::vector<LocalSolution<ComplexBall>>>
      exact_shadow_factorized_basis;
  if (normalized_matching_frame &&
      exact_lattice.exact_shadow_fused_local_basis.has_value()) {
    // The exact transformed local columns are the authority which removes
    // the terminal right-frame cancellation.  Specialize them, evaluate the
    // overlap matrices, solve, and certify under one precision lease.  A
    // later precision increase cannot recover coefficients already rounded
    // during Rational-to-Acb specialization.
    exact_shadow_precision.emplace(
        kExactShadowExtraPrecisionBits);
    const auto retained_taylor_width = basis_taylor_widths.front();
    if (std::any_of(
            basis_taylor_widths.begin(), basis_taylor_widths.end(),
            [&](const auto width) {
              return width != retained_taylor_width;
            }))
      throw std::logic_error(
          "exact Rational-shadow transformed basis has nonuniform retained Taylor widths");
    exact_shadow_physical_transformed_basis.emplace(
        dimension, FiniteLaurentVector<ComplexBall>());
    for (auto& row : *exact_shadow_physical_transformed_basis)
      row.reserve(dimension);
    exact_shadow_factorized_basis.emplace();
    exact_shadow_factorized_basis->reserve(dimension);
    for (std::size_t column = 0; column < dimension; ++column) {
      const auto& exact_local =
          (*exact_lattice.exact_shadow_fused_local_basis)[column];
      if (retained_taylor_width > exact_local.taylor_width())
        throw std::logic_error(
            "exact Rational-shadow transformed local is narrower than its Acb target");
      auto exact_options = basis_options;
      exact_options.t_order_reduction =
          static_cast<std::uint32_t>(
              exact_local.taylor_width() - retained_taylor_width);
      exact_shadow_factorized_basis->push_back(
          specialize_rational_local_solution_to_acb(
              exact_local,
              checkpoint_identity +
                  ":exact-shadow-factorized-column:" +
                  std::to_string(column)));
      const auto evaluated = evaluate_local_solution(
          exact_shadow_factorized_basis->back(),
          basis_point, exact_options);
      const EpsilonWindow exact_evaluation_window{
          evaluated.value.epsilon.min_power,
          evaluation_window.complete_max};
      auto frames = acb_evaluation_frames(
          evaluated.value, exact_evaluation_window,
          "exact Rational-shadow transformed basis at column " +
              std::to_string(column));
      for (std::size_t row = 0; row < dimension; ++row)
        (*exact_shadow_physical_transformed_basis)[row].push_back(
            std::move(frames[row]));
    }
    if (const auto* diagnostics =
            std::getenv("DE2_DIAGNOSTIC_TERMINAL_STATE");
        diagnostics != nullptr && std::string(diagnostics) == "1") {
      const auto specialized_saturation =
          specialize_exact_rational_laurent_matrix_to_acb(
              exact_lattice.saturation.transformation);
      const auto numeric_physical_transformed_basis =
          right_multiply_finite_by_exact_laurent(
              evaluated_basis, specialized_saturation);
      normal_frame_attempt[
          "exact_shadow_vs_numeric_physical_transformed_basis"] =
          diagnose_acb_laurent_matrix_overlap(
              *exact_shadow_physical_transformed_basis,
              numeric_physical_transformed_basis,
              "exact Rational-shadow F*T versus retained Acb F*T");
    }

    const auto receiving_owner =
        basis.front()->retained_equation_owner();
    if (!receiving_owner)
      throw std::logic_error(
          "exact Rational-shadow transformed basis lost its receiving SCC owner");
    exact_shadow_normalized_transformed_basis.emplace(
        dimension, FiniteLaurentVector<ComplexBall>());
    for (auto& row : *exact_shadow_normalized_transformed_basis)
      row.reserve(dimension);
    for (std::size_t column = 0; column < dimension; ++column) {
      FiniteLaurentVector<ComplexBall> physical_column;
      physical_column.reserve(dimension);
      for (std::size_t row = 0; row < dimension; ++row)
        physical_column.push_back(
            (*exact_shadow_physical_transformed_basis)[row][column]);
      auto normalized =
          receiving_owner->normalize_acb_matching_vector(
              physical_column, basis_point);
      if (!normalized.has_value() ||
          normalized->second != matching_frame_identity)
        throw std::logic_error(
            "exact Rational-shadow transformed basis lost its selected left normal frame");
      for (std::size_t row = 0; row < dimension; ++row)
        (*exact_shadow_normalized_transformed_basis)[row].push_back(
            std::move(normalized->first[row]));
    }
    // The high-precision Acb columns now retain the exact-shadow association
    // needed by terminal materialization.  Release only the Rational copies.
    exact_lattice.exact_shadow_fused_local_basis.reset();
  }
  if (expected_singular_request.has_value()) {
    const auto proof = json::parse(exact_lattice.canonical_witness);
    const auto& proof_object = as_object(
        proof, "selected singular-SCC lattice proof");
    normal_frame_attempt["lattice_transformation"] =
        required_string(proof_object, "transformation");
    normal_frame_attempt["lattice_coordinate_composition"] =
        normalized_matching_frame &&
            exact_lattice.acb_transformation.has_value()
        ? "right-normalization-inverse-times-retained-lattice"
        : "direct-selected-lattice";
  }
  std::optional<ExactLaurentMatrix<ComplexBall>>
      selected_acb_right_materialization_preconditioner;
  const bool diagnostic_normalized_match_authority = [] {
    const auto* raw = std::getenv(
        "DE2_DIAGNOSTIC_NORMALIZED_MATCH_AUTHORITY");
    return raw != nullptr && std::string(raw) == "1";
  }();
  const bool diagnostic_transformed_weight_floor = [] {
    const auto* raw = std::getenv(
        "DE2_DIAGNOSTIC_TRANSFORMED_WEIGHT_FLOOR");
    return raw != nullptr && std::string(raw) == "1";
  }();
  const auto transformed_weight_structural_min_power =
      diagnostic_transformed_weight_floor &&
              normalized_matching_frame &&
              exact_shadow_physical_transformed_basis.has_value()
          ? std::optional<std::int32_t>(
                incoming->solution().epsilon.min_power)
          : std::nullopt;
  if (transformed_weight_structural_min_power.has_value()) {
    normal_frame_attempt["transformed_weight_structural_min_power"] =
        *transformed_weight_structural_min_power;
    normal_frame_attempt["transformed_weight_floor_authority"] =
        "exact-valuation-zero-saturated-basis";
  }
  const auto fused_base_transformed_basis = [&](
      const ExactLaurentMatrix<ComplexBall>& base_transformation)
      -> std::optional<FiniteLaurentMatrix<ComplexBall>> {
    if (!normalized_matching_frame ||
        !left_normalized_matching_basis.has_value() ||
        !matching_right_normalization_matrix.has_value())
      return std::nullopt;
    if (exact_shadow_normalized_transformed_basis.has_value()) {
      normal_frame_attempt["right_coordinate_association"] =
          "exact-rational-physical-F*T-high-precision-midpoint-proposal-physical-residual";
      normal_frame_attempt["midpoint_proposal_extra_precision_bits"] =
          2048;
      return *exact_shadow_normalized_transformed_basis;
    }
    normal_frame_attempt["right_coordinate_association"] =
        "finite-acb-frame-fallback";
    auto right_composition =
        right_multiply_finite_by_exact_laurent(
            *matching_right_normalization_matrix,
            base_transformation);
    return right_multiply_finite_laurent_matrices(
        *left_normalized_matching_basis, right_composition);
  };
  const auto exact_physical_residual_correction_map = [&]()
      -> std::optional<std::function<FiniteLaurentVector<ComplexBall>(
          const FiniteLaurentVector<ComplexBall>&)>> {
    if (!normalized_matching_frame)
      return std::function<FiniteLaurentVector<ComplexBall>(
          const FiniteLaurentVector<ComplexBall>&)>{
          [](const FiniteLaurentVector<ComplexBall>& residual) {
            return residual;
          }};
    const auto receiving_owner =
        basis.front()->retained_equation_owner();
    if (!receiving_owner)
      throw std::logic_error(
          "exact Rational-shadow physical residual lost its receiving SCC owner");
    return [receiving_owner, basis_point, matching_frame_identity,
            checkpoint_identity](
               const FiniteLaurentVector<ComplexBall>& physical_residual) {
      auto normalized =
          receiving_owner->normalize_acb_matching_vector(
              physical_residual, basis_point);
      if (!normalized.has_value() ||
          normalized->second != matching_frame_identity)
        throw std::logic_error(
            checkpoint_identity +
            ": physical residual correction lost its selected left normal frame");
      return std::move(normalized->first);
    };
  };
  const auto finite_physical_transformed_basis = [&](
      const ExactLaurentMatrix<ComplexBall>& transformation)
      -> std::optional<FiniteLaurentMatrix<ComplexBall>> {
    if (!normalized_matching_frame)
      return right_multiply_finite_by_exact_laurent(
          evaluated_basis, transformation);
    if (!matching_right_normalization_matrix.has_value())
      return std::nullopt;
    auto physical_right_composition =
        right_multiply_finite_by_exact_laurent(
            *matching_right_normalization_matrix,
            transformation);
    return right_multiply_finite_laurent_matrices(
        evaluated_basis, physical_right_composition);
  };
  const auto run_refinement = [&](
      bool force_midpoint_proposal = false,
      bool factorized_authoritative_rhs_proposal = false) {
    selected_acb_right_materialization_preconditioner.reset();
    const auto refinement_context =
        checkpoint_identity + ":refined-acb-match[formal-negative-zero=" +
        (exact_lattice.exact_formal_negative_coefficients_are_zero
             ? std::string("true")
             : std::string("false")) +
        ",acb-transformation=" +
        (exact_lattice.acb_transformation.has_value()
             ? std::string("true")
             : std::string("false")) +
        "]";
    // Exact epsilon-lattice saturation T removes formal degeneracy, but its
    // epsilon^0 value matrix A=(F*T)|_0 can still be catastrophically
    // ill-conditioned.  A left preconditioner stabilizes the coefficient
    // solve only; returning its weights to F*T recreates the same destructive
    // cancellation when the receiving LocalSolution is materialized.  For a
    // singular receiver, normalize the retained right frame at the actual
    // match point with the fixed midpoint inverse P=A_mid^-1.  The physical
    // residual below remains authoritative, while downstream materialization
    // can preserve the numerically stable association ((F*T)*P)*y.
    std::optional<ExactLaurentMatrix<ComplexBall>>
        specialized_exact_transformation;
    const auto& base_transformation =
        exact_lattice.acb_transformation.has_value()
        ? *exact_lattice.acb_transformation
        : specialized_exact_transformation.emplace(
              specialize_exact_rational_laurent_matrix_to_acb(
                  exact_lattice.saturation.transformation));
    if (expected_singular_request.has_value()) {
      try {
        normal_frame_attempt["right_preconditioner_status"] =
            "building-leading-midpoint-inverse";
        auto retained_base_basis =
            fused_base_transformed_basis(base_transformation);
        const auto saturated_basis =
            retained_base_basis.has_value()
            ? std::move(*retained_base_basis)
            : right_multiply_finite_by_exact_laurent(
                  matching_basis, base_transformation);
        matching_detail::DenseScalarMatrix<ComplexBall> leading(
            dimension, std::vector<ComplexBall>(dimension));
        bool complete_at_zero = true;
        for (std::size_t row = 0; row < dimension; ++row)
          for (std::size_t column = 0; column < dimension; ++column) {
            const auto& frame = saturated_basis[row][column];
            if (frame.min_power() > 0 || frame.complete_max() < 0) {
              complete_at_zero = false;
              continue;
            }
            leading[row][column] = frame.coefficient(0);
          }
        auto dense_preconditioner = complete_at_zero
            ? matching_detail::acb_midpoint_inverse_preconditioner(
                  leading)
            : std::nullopt;
        if (!complete_at_zero)
          normal_frame_attempt["right_preconditioner_status"] =
              "epsilon-zero-leading-matrix-incomplete";
        else if (!dense_preconditioner.has_value())
          normal_frame_attempt["right_preconditioner_status"] =
              "midpoint-inverse-unavailable";
        if (dense_preconditioner.has_value()) {
          const auto normalized_leading =
              matching_detail::multiply_dense_matrices(
                  leading, *dense_preconditioner,
                  refinement_context +
                      ":constant-right-leading-verification");
          if (!matching_detail::certified_full_rank_plan(
                   normalized_leading, 0).has_value())
            throw MatchingArithmeticError(
                MatchingArithmeticErrorCode::AmbiguousZero,
                refinement_context +
                    ": constant right frame did not leave a certified full-rank epsilon-zero value matrix");
          ExactLaurentMatrix<ComplexBall> right_preconditioner(
              dimension,
              std::vector<ExactLaurentPolynomial<ComplexBall>>(
                  dimension));
          for (std::size_t row = 0; row < dimension; ++row)
            for (std::size_t column = 0; column < dimension; ++column)
              if (!(*dense_preconditioner)[row][column].is_zero())
                right_preconditioner[row][column] =
                    ExactLaurentPolynomial<ComplexBall>::monomial(
                        0, (*dense_preconditioner)[row][column]);
          auto conditioned_transformation =
              multiply_exact_laurent_matrices(
                  base_transformation, right_preconditioner);
          auto conditioned_transformed_basis =
              fused_base_transformed_basis(base_transformation);
          if (conditioned_transformed_basis.has_value())
            conditioned_transformed_basis =
                right_multiply_finite_by_exact_laurent(
                    *conditioned_transformed_basis,
                    right_preconditioner);
          auto conditioned_physical_residual_basis =
              normalized_matching_frame
              ? exact_shadow_physical_transformed_basis
              : std::nullopt;
          if (conditioned_physical_residual_basis.has_value())
            conditioned_physical_residual_basis =
                right_multiply_finite_by_exact_laurent(
                    *conditioned_physical_residual_basis,
                    right_preconditioner);
          else if (factorized_authoritative_rhs_proposal)
            conditioned_physical_residual_basis =
                finite_physical_transformed_basis(
                    conditioned_transformation);
          auto physical_residual_correction_map =
              exact_physical_residual_correction_map();
          auto conditioned = refine_acb_finite_laurent_match(
              matching_basis, matching_incoming,
              conditioned_transformation, refinement,
              refinement_context + ":constant-right-conditioned",
              exact_lattice
                  .exact_formal_negative_coefficients_are_zero,
              true,
              conditioned_transformed_basis.has_value()
                  ? &*conditioned_transformed_basis : nullptr,
              conditioned_physical_residual_basis.has_value() &&
                      !diagnostic_normalized_match_authority
                  ? &*conditioned_physical_residual_basis : nullptr,
              conditioned_physical_residual_basis.has_value() &&
                      !diagnostic_normalized_match_authority
                  ? &incoming_value : nullptr,
              physical_residual_correction_map.has_value() &&
                      !diagnostic_normalized_match_authority
                  ? &*physical_residual_correction_map : nullptr,
              force_midpoint_proposal ||
                  conditioned_physical_residual_basis.has_value(),
              std::size_t{0},
              transformed_weight_structural_min_power,
              factorized_authoritative_rhs_proposal);
          normal_frame_attempt["right_preconditioner_residual"] =
              conditioned.residual_history.empty()
              ? json::value(nullptr)
              : json::value(
                    encode_acb_match_residual_diagnostics(
                        conditioned.residual_history.back()));
          if (!conditioned.residual_history.empty() &&
              conditioned.residual_history.back().verdict ==
                  AcbMatchingResidualVerdict::Pass &&
              conditioned.residual_history.back()
                  .complete_through_required) {
            selected_acb_right_materialization_preconditioner =
                std::move(right_preconditioner);
            normal_frame_attempt["right_preconditioner_status"] =
                "certified";
            return conditioned;
          }
          normal_frame_attempt["right_preconditioner_status"] =
              "residual-not-certified";
        }
      } catch (const MatchingArithmeticError& error) {
        normal_frame_attempt["right_preconditioner_status"] =
            "arithmetic-rejected";
        normal_frame_attempt["right_preconditioner_detail"] =
            error.what();
        // A midpoint-derived P is only a conditioning proposal.  The exact
        // lattice and the unconditioned physical residual remain the
        // authoritative fallback.
      }
    }
    if (exact_lattice.acb_transformation.has_value()) {
      auto stable_transformed_basis =
          fused_base_transformed_basis(
              *exact_lattice.acb_transformation);
      auto physical_residual_correction_map =
          exact_physical_residual_correction_map();
      auto factorized_physical_basis =
          factorized_authoritative_rhs_proposal &&
                  !exact_shadow_physical_transformed_basis.has_value()
          ? finite_physical_transformed_basis(
                *exact_lattice.acb_transformation)
          : std::nullopt;
      const auto* authoritative_physical_basis =
          normalized_matching_frame &&
                  !diagnostic_normalized_match_authority &&
                  exact_shadow_physical_transformed_basis.has_value()
          ? &*exact_shadow_physical_transformed_basis
          : factorized_physical_basis.has_value()
          ? &*factorized_physical_basis
          : nullptr;
      return refine_acb_finite_laurent_match(
          matching_basis, matching_incoming,
          *exact_lattice.acb_transformation, refinement,
          refinement_context, false, false,
          stable_transformed_basis.has_value()
              ? &*stable_transformed_basis : nullptr,
          authoritative_physical_basis,
          authoritative_physical_basis != nullptr
              ? &incoming_value : nullptr,
          authoritative_physical_basis != nullptr &&
                  physical_residual_correction_map.has_value()
              ? &*physical_residual_correction_map : nullptr,
          force_midpoint_proposal ||
              authoritative_physical_basis != nullptr,
          std::size_t{0},
          transformed_weight_structural_min_power,
          factorized_authoritative_rhs_proposal);
    }
    auto specialized_transformation =
        specialize_exact_rational_laurent_matrix_to_acb(
            exact_lattice.saturation.transformation);
    auto stable_transformed_basis =
        fused_base_transformed_basis(specialized_transformation);
    auto physical_residual_correction_map =
        exact_physical_residual_correction_map();
    auto factorized_physical_basis =
        factorized_authoritative_rhs_proposal &&
                !exact_shadow_physical_transformed_basis.has_value()
        ? finite_physical_transformed_basis(
              specialized_transformation)
        : std::nullopt;
    const auto* authoritative_physical_basis =
        normalized_matching_frame &&
                !diagnostic_normalized_match_authority &&
                exact_shadow_physical_transformed_basis.has_value()
        ? &*exact_shadow_physical_transformed_basis
        : factorized_physical_basis.has_value()
        ? &*factorized_physical_basis
        : nullptr;
    return refine_acb_finite_laurent_match(
        matching_basis, matching_incoming,
        specialized_transformation, refinement,
        refinement_context,
        exact_lattice.exact_formal_negative_coefficients_are_zero,
        true,
        stable_transformed_basis.has_value()
            ? &*stable_transformed_basis : nullptr,
        authoritative_physical_basis,
        authoritative_physical_basis != nullptr
            ? &incoming_value : nullptr,
        authoritative_physical_basis != nullptr &&
                physical_residual_correction_map.has_value()
            ? &*physical_residual_correction_map : nullptr,
        force_midpoint_proposal ||
            authoritative_physical_basis != nullptr,
        std::size_t{0},
        transformed_weight_structural_min_power,
        factorized_authoritative_rhs_proposal);
  };
  RefinedAcbLaurentMatch refined;
  try {
    refined = run_refinement();
  } catch (const MatchingArithmeticError& error) {
    if (!normalized_matching_frame ||
        error.code != MatchingArithmeticErrorCode::InsufficientCompleteWindow)
      throw;
    // A terminal singular tile is an unbounded consumer of the physical
    // matching coordinates: a tiny coefficient error in a non-integrable
    // mode can remain invisible at the handoff and dominate the endpoint
    // integral.  Falling back to the epsilon-singular physical frame in that
    // case turns a private-reservoir deficit into a confidently wrong
    // result.  Propagate the honest complete-window error so the existing
    // ladder retry adds the missing epsilon orders and rebuilds this match in
    // its Fuchsian normal frame.
    if (require_normalized_singular_frame) throw;
    // The Fuchsian normal frame is a numerical conditioning optimization,
    // not part of the mathematical matching contract.  Its finite Laurent
    // transforms can consume all private overlap late in a long arm even
    // though the original physical frames still cover the requested result.
    // Fall back transactionally and rebuild the exact lattice in that frame;
    // never turn a conditioning aid into a spurious hard failure.
    normal_frame_attempt["status"] =
        "normalized-refinement-insufficient-window";
    normal_frame_attempt["arithmetic_detail"] = error.what();
    normal_frame_attempt["arithmetic_epsilon_power"] =
        error.epsilon_power.has_value()
        ? json::value(*error.epsilon_power)
        : json::value(nullptr);
    reset_to_physical_matching_frame();
    exact_lattice = make_exact_lattice();
    refined = run_refinement();
  }

  // Reaching factorization is not enough to admit the optional SCC normal
  // frame.  A finite V^-1/F/V replay may return an honest residual which no
  // longer reaches the requested prefix, even though the untransformed
  // physical columns still do.  Retrying with a wider owner can then change
  // the prepared V/V^-1 rectangle and make the reported complete edge move
  // backwards.  That violates the reservoir monotonicity required by the
  // caller: for one physical operator/chart, widening the requested top may
  // append coefficients but must not invalidate the old prefix.
  //
  // Treat the normal frame strictly as a conditioning optimization.  If its
  // first residual does not cover the required prefix, fail closed to the
  // physical-parent frame before publishing retry metadata.  This is a real
  // recomputation with the physical exact lattice, never a diagnostic clamp.
  if (normalized_matching_frame &&
      (refined.residual_history.empty() ||
       !refined.residual_history.back().complete_through_required ||
       refined.residual_history.back().verdict !=
           AcbMatchingResidualVerdict::Pass)) {
    if (require_normalized_singular_frame &&
        !refined.residual_history.empty() &&
        !refined.residual_history.back().complete_through_required)
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::InsufficientCompleteWindow,
          checkpoint_identity +
              ": terminal singular Fuchsian residual does not reach the required epsilon prefix",
          std::nullopt, std::nullopt,
          refined.residual_history.back().complete_window.complete_max);
    normal_frame_attempt["status"] =
        "normalized-residual-not-certified";
    if (!refined.residual_history.empty())
      normal_frame_attempt["normalized_residual"] =
          encode_acb_match_residual_diagnostics(
              refined.residual_history.back());
    if (!require_normalized_singular_frame) {
      reset_to_physical_matching_frame();
      exact_lattice = make_exact_lattice();
      refined = run_refinement();
    }
  }

  if (normalized_matching_frame) {
    const auto receiving_owner = basis.front()->retained_equation_owner();
    auto terminal_right_transformation = receiving_owner
        ? receiving_owner
              ->right_denormalization_acb_matching_matrix()
        : std::nullopt;
    if (!terminal_right_transformation.has_value())
      throw std::logic_error(
          "receiving SCC could not retain its terminal right normal-frame transformation");
    auto physical_weights = receiving_owner
        ? receiving_owner->denormalize_acb_matching_weights(refined.weights)
        : std::nullopt;
    if (!physical_weights.has_value())
      throw std::logic_error(
          "receiving SCC could not return Fuchsian matching weights to the physical basis");
    auto physical_options = refinement;
    physical_options.required_min_power = evaluation_window.min_power;
    bool physical_prefix_preserved = false;
    std::optional<matching_detail::AcbResidualEvaluation>
        physical_certificate;
    try {
      if (exact_shadow_physical_transformed_basis.has_value()) {
        auto physical_transformed_basis =
            *exact_shadow_physical_transformed_basis;
        if (selected_acb_right_materialization_preconditioner
                .has_value())
          physical_transformed_basis =
              right_multiply_finite_by_exact_laurent(
                  physical_transformed_basis,
                  *selected_acb_right_materialization_preconditioner);
        auto certified =
            matching_detail::evaluate_acb_matching_residual(
                physical_transformed_basis,
                refined.transformed_weights, incoming_value,
                physical_options,
                checkpoint_identity +
                    ":physical-via-exact-rational-local-basis");
        physical_prefix_preserved =
            certified.diagnostics.complete_through_required &&
            certified.diagnostics.verdict ==
                AcbMatchingResidualVerdict::Pass;
        if (physical_prefix_preserved)
          residual_certificate_identity =
              exact_lattice.identity +
              ":exact-rational-local-basis-physical-residual";
        if (!physical_prefix_preserved) {
          normal_frame_attempt["status"] =
              "exact-local-physical-residual-not-certified";
          normal_frame_attempt["physical_residual"] =
              encode_acb_match_residual_diagnostics(
                  certified.diagnostics);
        }
        // Retain every complete physical-authority verdict.  In particular,
        // a complete Inconclusive candidate must survive as an uncertified
        // StoredRefinedAcbMatch so the caller can retry shorter common Taylor
        // prefixes.  Only Pass receives a publication identity.
        physical_certificate = std::move(certified);
      } else if (auto pushed_residual =
                     receiving_owner
                         ->pushforward_acb_matching_residual(
                             refined.residual,
                             matching_frame_identity, basis_point);
                 pushed_residual.has_value()) {
        if (pushed_residual->normal_frame_identity !=
                matching_frame_identity ||
            pushed_residual->certificate_identity.empty())
          throw std::logic_error(
              "receiving SCC residual pushforward lost its exact frame binding");
        auto certified =
            matching_detail::certify_precomputed_acb_matching_residual(
                evaluated_basis, *physical_weights, incoming_value,
                std::move(pushed_residual->residual), physical_options,
                checkpoint_identity +
                    ":physical-prefix-via-exact-normal-frame");
        physical_prefix_preserved =
            certified.diagnostics.complete_through_required &&
            certified.diagnostics.verdict ==
              AcbMatchingResidualVerdict::Pass;
        if (physical_prefix_preserved)
          residual_certificate_identity =
              std::move(pushed_residual->certificate_identity);
        if (!physical_prefix_preserved) {
          normal_frame_attempt["status"] =
              "pushed-physical-residual-not-certified";
          normal_frame_attempt["physical_residual"] =
              encode_acb_match_residual_diagnostics(
                  certified.diagnostics);
        }
        physical_certificate = std::move(certified);
      } else {
        normal_frame_attempt["status"] =
            "residual-pushforward-unavailable";
      }
    } catch (const MatchingArithmeticError& error) {
      if (error.code !=
          MatchingArithmeticErrorCode::InsufficientCompleteWindow)
        throw;
      normal_frame_attempt["status"] =
          "pushed-physical-residual-insufficient-window";
      normal_frame_attempt["arithmetic_detail"] = error.what();
      normal_frame_attempt["arithmetic_epsilon_power"] =
          error.epsilon_power.has_value()
          ? json::value(*error.epsilon_power)
          : json::value(nullptr);
    }
    if (physical_certificate.has_value() &&
        physical_certificate->diagnostics.complete_through_required &&
        physical_certificate->diagnostics.verdict ==
            AcbMatchingResidualVerdict::Inconclusive) {
      try {
        const auto midpoint_basis =
            matching_detail::acb_midpoint_matrix(evaluated_basis);
        const auto midpoint_weights =
            matching_detail::acb_midpoint_vector(*physical_weights);
        const auto midpoint_incoming =
            matching_detail::acb_midpoint_vector(incoming_value);
        const auto midpoint_probe =
            matching_detail::evaluate_acb_matching_residual(
                midpoint_basis, midpoint_weights, midpoint_incoming,
                physical_options,
                checkpoint_identity +
                    ": physical residual zero-radius midpoint probe");
        normal_frame_attempt["physical_midpoint_residual"] =
            encode_acb_match_residual_diagnostics(
                midpoint_probe.diagnostics);
        normal_frame_attempt["physical_clearance_source"] =
            midpoint_probe.diagnostics.complete_through_required &&
                    midpoint_probe.diagnostics.verdict ==
                        AcbMatchingResidualVerdict::Pass
                ? "propagated-enclosure"
                : "candidate-midpoint-defect";
        if (normal_frame_attempt.at("physical_clearance_source") ==
            "propagated-enclosure") {
          json::object source_probes;
          const auto record_source =
              [&](const char* source,
                  const FiniteLaurentMatrix<ComplexBall>& probe_basis,
                  const FiniteLaurentVector<ComplexBall>& probe_weights,
                  const FiniteLaurentVector<ComplexBall>& probe_incoming) {
                const auto probe =
                    matching_detail::evaluate_acb_matching_residual(
                        probe_basis, probe_weights, probe_incoming,
                        physical_options,
                        checkpoint_identity +
                            ": physical residual " + source +
                            " enclosure probe");
                source_probes[source] =
                    encode_acb_match_residual_diagnostics(
                        probe.diagnostics);
                return probe.diagnostics;
              };
          const auto basis_probe =
              record_source("basis", evaluated_basis, midpoint_weights,
                            midpoint_incoming);
          const auto weights_probe =
              record_source("weights", midpoint_basis, *physical_weights,
                            midpoint_incoming);
          const auto incoming_probe =
              record_source("incoming", midpoint_basis, midpoint_weights,
                            incoming_value);
          normal_frame_attempt[
              "physical_clearance_source_probes"] =
              std::move(source_probes);
          if (acb_factorized_rhs_proposal_applicable(
                  physical_certificate->diagnostics,
                  midpoint_probe.diagnostics, basis_probe,
                  weights_probe, incoming_probe)) {
            auto factorized_candidate =
                run_refinement(true, true);
            normal_frame_attempt[
                "factorized_authoritative_rhs_columns"] =
                factorized_candidate
                    .factorized_authoritative_rhs_columns;
            if (!factorized_candidate.residual_history.empty() &&
                acb_matching_residual_certified_pass(
                    factorized_candidate.residual_history.back())) {
              refined = std::move(factorized_candidate);
              physical_weights =
                  receiving_owner
                      ->denormalize_acb_matching_weights(
                          refined.weights);
              if (!physical_weights.has_value())
                throw std::logic_error(
                    "factorized authoritative-rhs proposal could not return matching weights to the physical basis");
              physical_certificate =
                  matching_detail::AcbResidualEvaluation{
                      refined.residual,
                      refined.residual_history.back()};
              physical_prefix_preserved = true;
              residual_certificate_identity =
                  exact_lattice.identity +
                  ":factorized-authoritative-rhs-physical-residual";
              normal_frame_attempt[
                  "factorized_authoritative_rhs_proposal"] =
                  "certified";
              normal_frame_attempt["status"] =
                  "certified-factorized-authoritative-rhs";
              normal_frame_attempt["physical_residual"] =
                  encode_acb_match_residual_diagnostics(
                      refined.residual_history.back());
            } else {
              normal_frame_attempt[
                  "factorized_authoritative_rhs_proposal"] =
                  "full-ball-operator-residual-not-certified";
              if (!factorized_candidate.residual_history.empty())
                normal_frame_attempt[
                    "factorized_authoritative_rhs_residual"] =
                    encode_acb_match_residual_diagnostics(
                        factorized_candidate
                            .residual_history.back());
            }
          }
        }
      } catch (const std::exception& error) {
        normal_frame_attempt["physical_clearance_source"] =
            "midpoint-probe-unavailable";
        normal_frame_attempt["physical_midpoint_detail"] =
            error.what();
      }
    }
    if (diagnostic_normalized_match_authority) {
      // Diagnostic experiment: keep the candidate obtained in the retained
      // Fuchsian frame and treat the physical stored-Taylor residual as an
      // observed truncation defect, not as a correction source.  This path
      // is deliberately opt-in until the residual is bounded by a certified
      // Taylor-tail envelope.
      if (require_normalized_singular_frame)
        selected_terminal_normal_frame_right_transformation =
            std::move(*terminal_right_transformation);
      refined.weights = std::move(*physical_weights);
      normal_frame_attempt["status"] =
          "diagnostic-normalized-authority";
      normal_frame_attempt["physical_residual_authority"] = false;
    } else if (!physical_prefix_preserved) {
      if (require_normalized_singular_frame) {
        if (!physical_certificate.has_value())
          throw MatchingArithmeticError(
              MatchingArithmeticErrorCode::InvalidSaturationLattice,
              checkpoint_identity +
                  ": terminal singular Fuchsian pushforward did not retain a physical residual certificate; normal_frame_attempt=" +
                  json::serialize(normal_frame_attempt));
        if (!physical_certificate->diagnostics.complete_through_required)
          throw MatchingArithmeticError(
              MatchingArithmeticErrorCode::InsufficientCompleteWindow,
              checkpoint_identity +
                  ": terminal singular Fuchsian pushforward does not reach the required physical epsilon prefix",
              std::nullopt, std::nullopt,
              physical_certificate->diagnostics.complete_window.complete_max);
        // Complete Inconclusive and complete Fail are valid diagnostic
        // outcomes, not malformed lattices.  Preserve the physical authority
        // and terminal frame but leave the match uncertified.  The outer
        // prefix search retries only Inconclusive; Fail remains final, and
        // every materialization consumer independently requires Pass.
        selected_terminal_normal_frame_right_transformation =
            std::move(*terminal_right_transformation);
        refined.weights = std::move(*physical_weights);
        refined.residual = std::move(physical_certificate->residual);
        refined.residual_history.back() =
            std::move(physical_certificate->diagnostics);
      } else {
        reset_to_physical_matching_frame();
        exact_lattice = make_exact_lattice();
        refined = run_refinement();
      }
    } else {
      if (require_normalized_singular_frame)
        selected_terminal_normal_frame_right_transformation =
            std::move(*terminal_right_transformation);
      refined.weights = std::move(*physical_weights);
      refined.residual = std::move(physical_certificate->residual);
      refined.residual_history.back() =
          std::move(physical_certificate->diagnostics);
      normal_frame_attempt["status"] = "certified";
      normal_frame_attempt["physical_residual"] =
          encode_acb_match_residual_diagnostics(
              refined.residual_history.back());
    }
  }
  if (!normalized_matching_frame &&
      !refined.residual_history.empty() &&
      refined.residual_history.back().complete_through_required &&
      refined.residual_history.back().verdict ==
          AcbMatchingResidualVerdict::Inconclusive) {
    // Regular and already-physical singular receivers do not enter the
    // Fuchsian pushforward block above, but their complete residual can fail
    // for exactly the same two distinct reasons: a defective midpoint
    // candidate, or an otherwise valid candidate whose input balls are too
    // wide.  Preserve the same source attribution so the ladder never sends
    // an incoming-boundary problem to a Taylor retry.
    try {
      auto physical_options = refinement;
      physical_options.required_min_power =
          evaluation_window.min_power;
      const auto midpoint_basis =
          matching_detail::acb_midpoint_matrix(evaluated_basis);
      const auto midpoint_weights =
          matching_detail::acb_midpoint_vector(refined.weights);
      const auto midpoint_incoming =
          matching_detail::acb_midpoint_vector(incoming_value);
      const auto midpoint_probe =
          matching_detail::evaluate_acb_matching_residual(
              midpoint_basis, midpoint_weights, midpoint_incoming,
              physical_options,
              checkpoint_identity +
                  ": physical-frame residual zero-radius midpoint probe");
      normal_frame_attempt["physical_midpoint_residual"] =
          encode_acb_match_residual_diagnostics(
              midpoint_probe.diagnostics);
      normal_frame_attempt["physical_clearance_source"] =
          midpoint_probe.diagnostics.complete_through_required &&
                  midpoint_probe.diagnostics.verdict ==
                      AcbMatchingResidualVerdict::Pass
              ? "propagated-enclosure"
              : "candidate-midpoint-defect";
      if (normal_frame_attempt.at("physical_clearance_source") ==
          "propagated-enclosure") {
        json::object source_probes;
        const auto evaluate_source =
            [&](const char* source,
                const FiniteLaurentMatrix<ComplexBall>& probe_basis,
                const FiniteLaurentVector<ComplexBall>& probe_weights,
                const FiniteLaurentVector<ComplexBall>& probe_incoming) {
              const auto probe =
                  matching_detail::evaluate_acb_matching_residual(
                      probe_basis, probe_weights, probe_incoming,
                      physical_options,
                      checkpoint_identity +
                          ": physical-frame residual " + source +
                          " enclosure probe");
              source_probes[source] =
                  encode_acb_match_residual_diagnostics(
                      probe.diagnostics);
              return probe.diagnostics;
            };
        const auto basis_probe =
            evaluate_source("basis", evaluated_basis, midpoint_weights,
                            midpoint_incoming);
        const auto weights_probe =
            evaluate_source("weights", midpoint_basis, refined.weights,
                            midpoint_incoming);
        const auto incoming_probe =
            evaluate_source("incoming", midpoint_basis, midpoint_weights,
                            incoming_value);
        normal_frame_attempt["physical_clearance_source_probes"] =
            std::move(source_probes);
        const auto authoritative_probe =
            refined.residual_history.back();
        if (acb_midpoint_weight_proposal_applicable(
                authoritative_probe,
                midpoint_probe.diagnostics, basis_probe, weights_probe,
                incoming_probe)) {
          // The solve has either magnified harmless input radii into the
          // candidate weights or copied an otherwise acceptable incoming
          // radius into them.  Re-solve the same finite Laurent problem using
          // zero-radius midpoint data only to propose coordinates.  The
          // untouched full-ball basis and incoming boundary still form the
          // authoritative forward residual inside run_refinement(), so this
          // cannot publish a midpoint-only decision.
          auto midpoint_candidate = run_refinement(true);
          if (!midpoint_candidate.residual_history.empty() &&
              acb_matching_residual_certified_pass(
                  midpoint_candidate.residual_history.back())) {
            refined = std::move(midpoint_candidate);
            normal_frame_attempt[
                "regular_midpoint_weight_proposal"] = "certified";
            normal_frame_attempt[
                "regular_midpoint_weight_proposal_residual"] =
                encode_acb_match_residual_diagnostics(
                    refined.residual_history.back());
          } else {
            normal_frame_attempt[
                "regular_midpoint_weight_proposal"] =
                "full-ball-residual-not-certified";
            if (!midpoint_candidate.residual_history.empty())
              normal_frame_attempt[
                  "regular_midpoint_weight_proposal_residual"] =
                  encode_acb_match_residual_diagnostics(
                      midpoint_candidate.residual_history.back());
          }
        }
        if (!acb_matching_residual_certified_pass(
                refined.residual_history.back()) &&
            acb_factorized_rhs_proposal_applicable(
                authoritative_probe, midpoint_probe.diagnostics,
                basis_probe, weights_probe, incoming_probe)) {
          auto factorized_candidate =
              run_refinement(true, true);
          normal_frame_attempt[
              "factorized_authoritative_rhs_columns"] =
              factorized_candidate
                  .factorized_authoritative_rhs_columns;
          if (!factorized_candidate.residual_history.empty() &&
              acb_matching_residual_certified_pass(
                  factorized_candidate.residual_history.back())) {
            refined = std::move(factorized_candidate);
            normal_frame_attempt[
                "factorized_authoritative_rhs_proposal"] =
                "certified";
            normal_frame_attempt["status"] =
                "certified-factorized-authoritative-rhs";
            normal_frame_attempt[
                "factorized_authoritative_rhs_residual"] =
                encode_acb_match_residual_diagnostics(
                    refined.residual_history.back());
          } else {
            normal_frame_attempt[
                "factorized_authoritative_rhs_proposal"] =
                "full-ball-operator-residual-not-certified";
            if (!factorized_candidate.residual_history.empty())
              normal_frame_attempt[
                  "factorized_authoritative_rhs_residual"] =
                  encode_acb_match_residual_diagnostics(
                      factorized_candidate
                          .residual_history.back());
          }
        }
      }
    } catch (const std::exception& error) {
      normal_frame_attempt["physical_clearance_source"] =
          "midpoint-probe-unavailable";
      normal_frame_attempt["physical_midpoint_detail"] =
          error.what();
    }
  }
  const std::string residual_frame_identity = normalized_matching_frame
      ? "physical-parent-frame"
      : matching_frame_identity;
  // A Rational lattice transformation used directly in the physical parent
  // frame is an exact t-independent identity for the complete receiving
  // local basis.  Retain it separately from the matching proof so downstream
  // materialization can preserve the stable (F*T)*y association.
  //
  // A two-sided SCC normal frame solves against L*F*R and returns certified
  // physical weights w=R*T*y.  R is retained as a finite frame with an
  // honest completeness edge, not as an exact Laurent polynomial.  It is
  // therefore unsound to publish a finite R*T product as an exact downstream
  // right transformation: doing so would turn every coefficient above that
  // edge into structural zero.  Keep the normal frame strictly internal to
  // matching and materialize the certified physical association F*w.
  //
  // Exact-right terminal factorization remains available in the physical
  // frame, where the Rational lattice transformation really is exact and
  // t-independent.
  std::optional<ExactLaurentMatrix<Rational>>
      exact_right_materialization_transformation;
  if (!normalized_matching_frame &&
      !exact_lattice.acb_transformation.has_value())
    exact_right_materialization_transformation =
        exact_lattice.saturation.transformation;
  std::optional<ExactLaurentMatrix<ComplexBall>>
      acb_right_materialization_transformation;
  if (normalized_matching_frame) {
    normal_frame_attempt["terminal_materialization"] =
        require_normalized_singular_frame
            ? exact_lattice.acb_transformation.has_value()
                ? "finite-normal-frame-acb-Levelt-right"
                : "finite-normal-frame-exact-formal-right"
            : "physical-F*w";
    // R^-1*T is a coordinate transformation for the normalized matching
    // basis B=L*F*R: B*(R^-1*T)=L*F*T.  It is not a right transformation of
    // the physical receiving basis F.  Only a terminal factorized consumer
    // retains both R and R^-1*T and can replay that association.  An
    // ordinary streamed hop must materialize the certified physical weights
    // w=R*T*y as F*w.
    if (require_normalized_singular_frame)
      acb_right_materialization_transformation =
          std::move(exact_lattice.acb_transformation);
  } else {
    acb_right_materialization_transformation =
        std::move(exact_lattice.acb_transformation);
  }

  auto exact_binding_basis =
      compact_matching_source_references(basis_sources);
  json::object exact_lattice_provenance{
      {"schema", "diffexp2-retained-exact-lattice-binding-v1"},
      {"witness_schema", exact_lattice.witness_schema},
      {"witness_identity", exact_lattice.identity},
      {"basis", std::move(exact_binding_basis)},
      {"basis_point_exact", basis_point.exact_coordinate},
      {"physical_match_point_exact", basis_physical_point},
      {"matching_frame_identity", matching_frame_identity},
      {"epsilon", json::object{{"min", matching_window.min_power},
                                {"max", matching_window.complete_max}}}};
  const auto exact_lattice_provenance_identity = json::serialize(
      canonical_json_value(exact_lattice_provenance));

  auto provenance_basis =
      compact_matching_source_references(basis_sources);
  json::object provenance{
      {"schema", "diffexp2-native-refined-acb-local-match-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"basis", std::move(provenance_basis)},
      {"incoming",
       compact_matching_source_reference(incoming_source)},
      {"basis_point_exact", basis_point.exact_coordinate},
      {"incoming_point_exact", incoming_point.exact_coordinate},
      {"physical_match_point_exact", basis_physical_point},
      {"matching_frame_identity", matching_frame_identity},
      {"residual_frame_identity", residual_frame_identity},
      {"epsilon", json::object{{"min", matching_window.min_power},
                                {"max", matching_window.complete_max},
                                {"required_complete_max",
                                 required_complete_max}}},
      {"exact_lattice_provenance_identity",
       exact_lattice_provenance_identity},
      {"refinement", json::object{{"relative_tolerance",
                                    relative_tolerance},
                                   {"max_steps",
                                    max_refinement_steps}}}};
  if (!residual_certificate_identity.empty())
    provenance["residual_certificate_identity"] =
        residual_certificate_identity;
  if (selected_acb_right_materialization_preconditioner.has_value())
    provenance["right_materialization_preconditioner"] =
        checkpoint_acb_laurent_matrix_record(
            *selected_acb_right_materialization_preconditioner);
  if (selected_terminal_normal_frame_right_transformation.has_value())
    provenance["terminal_normal_frame_right_transformation"] =
        checkpoint_frame_matrix_record(
            *selected_terminal_normal_frame_right_transformation);
  if (selected_terminal_normal_frame_right_transformation.has_value())
    provenance["terminal_normal_frame_exact_right_transformation"] =
        checkpoint_exact_laurent_matrix_record(
            exact_lattice.saturation.transformation);
  if (exact_shadow_factorized_basis.has_value() &&
      selected_terminal_normal_frame_right_transformation.has_value())
    provenance["terminal_exact_shadow_factorized_basis"] =
        json::object{
            {"columns", exact_shadow_factorized_basis->size()},
            {"extra_precision_bits",
             kExactShadowExtraPrecisionBits}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  if (std::getenv("DE2_DIAGNOSTIC_TERMINAL_STATE") != nullptr) {
    const auto* association =
        normal_frame_attempt.if_contains(
            "right_coordinate_association");
    if (association != nullptr) {
      std::map<std::int32_t, Magnitude> transformed_weight_scales;
      for (const auto& weight : refined.transformed_weights)
        for (std::int64_t power = weight.min_power();
             power <= std::min<std::int64_t>(
                          weight.complete_max(), 0);
             ++power) {
          const auto checked = matching_detail::checked_power(
              power, "terminal transformed-weight diagnostic power");
          const auto magnitude = Magnitude::upper_abs(
              weight.coefficient(checked));
          const auto found = transformed_weight_scales.find(checked);
          if (found == transformed_weight_scales.end() ||
              !(magnitude <= found->second))
            transformed_weight_scales.insert_or_assign(
                checked, magnitude);
        }
      json::array encoded_weight_scales;
      for (const auto& [power, magnitude] :
           transformed_weight_scales)
        encoded_weight_scales.push_back(json::object{
            {"power", power},
            {"absolute_upper_exact", magnitude.dump_exact()},
            {"absolute_upper_approx",
             magnitude.approximate_upper()}});
      std::fprintf(
          stderr,
          "[diffexp2 terminal diagnostic] matching right coordinate association=%s incoming_structural_min=%d transformed_weight_scales=%s elapsed_ms=%.3f\n",
          json::serialize(*association).c_str(),
          incoming->solution().epsilon.min_power,
          json::serialize(encoded_weight_scales).c_str(), elapsed_ms);
    }
  }
  if (!selected_terminal_normal_frame_right_transformation.has_value())
    exact_shadow_factorized_basis.reset();
  const auto retained_exact_shadow_extra_precision_bits =
      exact_shadow_factorized_basis.has_value()
      ? kExactShadowExtraPrecisionBits : std::size_t{0};
  auto result = std::make_shared<StoredRefinedAcbMatch>(
      match_handle, checkpoint_identity, provenance_identity,
      exact_lattice.identity, exact_lattice_provenance_identity,
      std::move(exact_lattice.canonical_witness),
      exact_lattice.witness_schema,
      std::move(basis_sources), std::move(incoming_source), basis_chart,
      incoming_chart, basis_point.exact_coordinate,
      incoming_point.exact_coordinate, basis_physical_point,
      matching_frame_identity, residual_frame_identity,
      residual_certificate_identity,
      matching_window,
      required_complete_max, dimension, relative_tolerance,
      max_refinement_steps, std::move(exact_lattice.saturation),
      std::move(exact_right_materialization_transformation),
      std::move(acb_right_materialization_transformation),
      std::move(selected_acb_right_materialization_preconditioner),
      std::move(selected_terminal_normal_frame_right_transformation),
      std::move(exact_shadow_factorized_basis),
      retained_exact_shadow_extra_precision_bits,
      std::move(exact_lattice.exact_shadow_equation_owner),
      std::move(refined), elapsed_ms);
  result->replace_normal_frame_attempt(
      std::move(normal_frame_attempt));
  return result;
}

struct TerminalMatchDiagnosticScanSpec {
  std::string label;
  std::optional<Rational> receiving_local;
  std::optional<std::size_t> taylor_width;
};

std::vector<TerminalMatchDiagnosticScanSpec>
parse_terminal_match_diagnostic_scan_specs(const char* raw) {
  std::vector<TerminalMatchDiagnosticScanSpec> specs;
  if (raw == nullptr || std::string(raw).empty()) return specs;
  const std::string input(raw);
  std::size_t begin = 0;
  while (begin <= input.size()) {
    const auto end = input.find(',', begin);
    auto token = input.substr(
        begin, end == std::string::npos ? std::string::npos
                                        : end - begin);
    const auto first = token.find_first_not_of(" \t");
    const auto last = token.find_last_not_of(" \t");
    if (first == std::string::npos)
      throw std::invalid_argument(
          "DE2_DIAGNOSTIC_TERMINAL_MATCH_SCAN contains an empty item");
    token = token.substr(first, last - first + 1);
    const auto separator = token.find('@');
    if (separator != std::string::npos &&
        token.find('@', separator + 1) != std::string::npos)
      throw std::invalid_argument(
          "terminal match scan items accept at most one @TaylorWidth suffix");
    const auto point_text = token.substr(0, separator);
    const auto width_text = separator == std::string::npos
        ? std::string() : token.substr(separator + 1);
    TerminalMatchDiagnosticScanSpec spec;
    spec.label = token;
    if (point_text != "same") {
      if (point_text.empty())
        throw std::invalid_argument(
            "terminal match scan requires an exact receiving local point or 'same'");
      spec.receiving_local = Rational(point_text);
    }
    if (separator != std::string::npos) {
      if (width_text.empty())
        throw std::invalid_argument(
            "terminal match scan @ suffix requires a Taylor width");
      std::size_t parsed = 0;
      const auto width = std::stoull(width_text, &parsed);
      if (parsed != width_text.size() || width == 0 ||
          width > std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument(
            "terminal match scan Taylor width must be a positive integer");
      spec.taylor_width = static_cast<std::size_t>(width);
    }
    specs.push_back(std::move(spec));
    if (specs.size() > 16)
      throw std::invalid_argument(
          "terminal match diagnostic scan accepts at most 16 items");
    if (end == std::string::npos) break;
    begin = end + 1;
  }
  return specs;
}

void run_failed_terminal_match_diagnostic_scans(
    const json::object& request,
    const std::vector<std::string>& basis_handles,
    const std::vector<std::shared_ptr<StoredLocalBase>>& erased_basis,
    const std::string& incoming_handle,
    const std::shared_ptr<StoredLocalBase>& erased_incoming,
    slong precision_bits,
    const std::string& active_session_configuration_identity,
    const std::optional<json::object>& expected_singular_request,
    bool require_normalized_singular_frame) {
  const auto specs = parse_terminal_match_diagnostic_scan_specs(
      std::getenv("DE2_DIAGNOSTIC_TERMINAL_MATCH_SCAN"));
  if (specs.empty()) return;
  if (!expected_singular_request.has_value() ||
      request.if_contains("native_singular_scc_saturation") == nullptr)
    throw std::invalid_argument(
        "terminal match diagnostic scan requires a retained singular-SCC request");
  if (erased_basis.empty())
    throw std::invalid_argument(
        "terminal match diagnostic scan requires a receiving basis");
  const auto basis =
      std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(
          erased_basis.front());
  const auto incoming =
      std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(
          erased_incoming);
  if (!basis || !incoming)
    throw std::invalid_argument(
        "terminal match diagnostic scan requires Acb locals");

  const auto original_checkpoint =
      required_string(request, "checkpoint_identity");
  json::array results;
  results.reserve(specs.size());
  for (std::size_t index = 0; index < specs.size(); ++index) {
    const auto& spec = specs[index];
    json::object item{
        {"scan", spec.label},
        {"requested_taylor_width",
         spec.taylor_width.has_value()
             ? json::value(
                   static_cast<std::uint64_t>(*spec.taylor_width))
             : json::value(nullptr)}};
    try {
      auto scan_request = request;
      auto scan_expected = *expected_singular_request;
      const auto scan_checkpoint =
          original_checkpoint + ":failed-terminal-scan:" +
          std::to_string(index);
      scan_request["checkpoint_identity"] = scan_checkpoint;
      scan_expected["match_checkpoint_identity"] =
          scan_checkpoint;

      auto scan_singular = as_object(
          scan_request.at("native_singular_scc_saturation"),
          "terminal match diagnostic singular request");
      scan_singular["match_checkpoint_identity"] =
          scan_checkpoint;
      if (spec.receiving_local.has_value()) {
        const auto& receiving_local = *spec.receiving_local;
        const auto receiving_center =
            Rational(basis->solution().chart.center_exact);
        const auto receiving_scale =
            Rational(basis->solution().chart.scale_exact);
        const auto producing_center =
            Rational(incoming->solution().chart.center_exact);
        const auto producing_scale =
            Rational(incoming->solution().chart.scale_exact);
        const auto physical =
            receiving_center + receiving_scale * receiving_local;
        const auto incoming_local =
            (physical - producing_center) / producing_scale;
        scan_request["basis_point"] =
            json::object{{"exact", receiving_local.str()}};
        scan_request["incoming_point"] =
            json::object{{"exact", incoming_local.str()}};
        // These diagnostic points are exact rationals constructed locally,
        // not retained CertifiedPlanPoint balls.  Let the ordinary exact
        // rational geometry check below recompute the common physical point
        // instead of falsely claiming a paired algebraic-point certificate.
        scan_request.erase(
            "certified_physical_match_point_exact");
        scan_expected["receiving_basis_point_exact"] =
            receiving_local.str();
        scan_expected["receiving_basis_point_sign"] =
            receiving_local.sign();
        scan_expected["physical_match_point_exact"] =
            physical.str();
        scan_singular["receiving_basis_point_exact"] =
            receiving_local.str();
        scan_singular["receiving_basis_point_sign"] =
            receiving_local.sign();
        scan_singular["physical_match_point_exact"] =
            physical.str();
        item["receiving_local"] = receiving_local.str();
        item["incoming_local"] = incoming_local.str();
        item["physical"] = physical.str();
      } else {
        item["receiving_local"] =
            as_object(scan_request.at("basis_point"),
                      "terminal match scan basis point")
                .at("exact");
        item["incoming_local"] =
            as_object(scan_request.at("incoming_point"),
                      "terminal match scan incoming point")
                .at("exact");
        item["physical"] =
            scan_request.at(
                "certified_physical_match_point_exact");
      }
      scan_request["native_singular_scc_saturation"] =
          scan_singular;
      const auto rebuilt = build_refined_acb_match_once(
          original_checkpoint + ":failed-terminal-scan-match:" +
              std::to_string(index),
          scan_request, basis_handles, erased_basis,
          incoming_handle, erased_incoming, precision_bits,
          active_session_configuration_identity, scan_expected,
          require_normalized_singular_frame,
          spec.taylor_width, false);
      item["status"] = "ok";
      item["match"] =
          rebuilt->compact_terminal_diagnostic_summary();
    } catch (const std::exception& error) {
      item["status"] = "error";
      item["detail"] = error.what();
    }
    results.push_back(std::move(item));
  }
  std::fprintf(
      stderr,
      "[diffexp2 failed terminal match scan] %s\n",
      json::serialize(json::object{
          {"schema",
           "diffexp2-failed-terminal-match-scan-v1"},
          {"checkpoint_identity", original_checkpoint},
          {"scan", std::move(results)}}).c_str());
}

std::shared_ptr<StoredRefinedAcbMatch> build_refined_acb_match(
    const std::string& match_handle, const json::object& request,
    const std::vector<std::string>& basis_handles,
    const std::vector<std::shared_ptr<StoredLocalBase>>& erased_basis,
    const std::string& incoming_handle,
    const std::shared_ptr<StoredLocalBase>& erased_incoming,
    slong precision_bits,
    const std::string& active_session_configuration_identity,
    const std::optional<json::object>& expected_singular_request =
        std::nullopt,
    bool require_normalized_singular_frame = false) {
  const auto started = std::chrono::steady_clock::now();
  std::size_t common_width = std::numeric_limits<std::size_t>::max();
  for (const auto& erased : erased_basis) {
    const auto local =
        std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(erased);
    if (!local)
      throw std::invalid_argument(
          "refined Acb matching requires Acb retained locals");
    common_width = std::min(common_width,
                            local->solution().taylor_width());
  }
  const auto incoming =
      std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(erased_incoming);
  if (!incoming)
    throw std::invalid_argument(
        "refined Acb matching requires an Acb incoming local");
  common_width = std::min(common_width,
                          incoming->solution().taylor_width());

  std::shared_ptr<StoredRefinedAcbMatch> full;
  std::exception_ptr full_arithmetic_failure;
  try {
    full = build_refined_acb_match_once(
        match_handle, request, basis_handles, erased_basis, incoming_handle,
        erased_incoming, precision_bits,
        active_session_configuration_identity, expected_singular_request,
        require_normalized_singular_frame, std::nullopt, true);
  } catch (const MatchingArithmeticError& error) {
    if (!acb_taylor_prefix_retryable(
            error.code, require_normalized_singular_frame))
      throw;
    full_arithmetic_failure = std::current_exception();
  }
  if (full && full->certified_for_materialization()) {
    full->replace_elapsed_ms(std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count());
    return full;
  }
  if (full &&
      full->final_residual_verdict() ==
          AcbMatchingResidualVerdict::Fail &&
      full->residual_complete_through_required()) {
    run_failed_terminal_match_diagnostic_scans(
        request, basis_handles, erased_basis, incoming_handle,
        erased_incoming, precision_bits,
        active_session_configuration_identity,
        expected_singular_request,
        require_normalized_singular_frame);
  }

  // A terminal singular match ultimately defines an endpoint functional.
  // A shorter finite Taylor sum is not interchangeable with that functional
  // unless the omitted singular and incoming tails have both been bounded.
  // build_refined_acb_match_once deliberately disables tail estimates, so a
  // prefix-only pass here would merely replace an inconclusive full sum by an
  // accidental equality of two truncations.  Keep ordinary regular matching
  // behavior unchanged while the general prefix path is migrated to a
  // certified-tail contract, but fail closed for the terminal normal frame.
  if (require_normalized_singular_frame) {
    if (full) {
      full->replace_elapsed_ms(std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count());
      return full;
    }
    std::rethrow_exception(full_arithmetic_failure);
  }

  // A stored Taylor order is a maximum work reservoir, not an obligation to
  // sum a numerically poisoned suffix.  A complete Inconclusive residual may
  // retry, as may a proof/factorization ambiguity caused by a pivot enclosure
  // overlapping zero.  A genuine residual Fail, an incomplete epsilon frame,
  // or a structural/invalid-lattice error remains authoritative.  Test every
  // shorter common prefix from largest to smallest, so increasing the stored
  // ExpansionOrder cannot make a prefix accepted at a lower order
  // unavailable.  The selected width is bound into every source record and
  // therefore into the exact-lattice and match provenance.
  if (full &&
      !acb_taylor_prefix_retryable(
          full->final_residual_verdict(),
          full->residual_complete_through_required(),
          require_normalized_singular_frame)) {
    full->replace_elapsed_ms(std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count());
    return full;
  }
  auto best = full;
  if (common_width <= 1) {
    if (best) {
      best->replace_elapsed_ms(std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count());
      return best;
    }
    std::rethrow_exception(full_arithmetic_failure);
  }
  for (auto width = common_width - 1; width > 0; --width) {
    std::shared_ptr<StoredRefinedAcbMatch> candidate;
    try {
      candidate = build_refined_acb_match_once(
          match_handle, request, basis_handles, erased_basis, incoming_handle,
          erased_incoming, precision_bits,
          active_session_configuration_identity, expected_singular_request,
          require_normalized_singular_frame, width, false);
    } catch (const MatchingArithmeticError& error) {
      if (!acb_taylor_prefix_retryable(error.code, false)) throw;
      // A narrower speculative prefix can lose the full-rank or complete
      // Laurent frame that the original attempt possessed.  It is not a new
      // failure of the requested match; keep searching and retain the best
      // completed diagnostic if no shorter prefix certifies.
      continue;
    }
    if (candidate->certified_for_materialization()) {
      candidate->replace_elapsed_ms(std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count());
      return candidate;
    }
    if (candidate->final_residual_verdict() ==
            AcbMatchingResidualVerdict::Inconclusive &&
        (!best ||
         candidate->has_better_inconclusive_certificate_than(*best)))
      best = std::move(candidate);
  }
  if (!best) std::rethrow_exception(full_arithmetic_failure);
  best->replace_elapsed_ms(std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count());
  return best;
}

double checkpoint_nonnegative_double(const json::value& raw,
                                     const char* label) {
  const auto value = as_double(raw, label);
  if (!std::isfinite(value) || value < 0.0)
    throw std::invalid_argument(std::string(label) +
                                " must be finite and nonnegative");
  return value;
}

SCCColumnProvenance parse_checkpoint_column_provenance(
    const json::value& raw) {
  const auto& object = as_object(raw, "checkpoint SCC-column provenance");
  require_exact_keys(object,
      {"scc", "scc_exact_identity", "seed_block", "basis_index",
       "exact_column_identity"}, "checkpoint SCC-column provenance");
  SCCColumnProvenance result{
      required_string(object, "scc"),
      required_string(object, "scc_exact_identity"),
      as_u32(object.at("seed_block"), "checkpoint seed block"),
      as_u32(object.at("basis_index"), "checkpoint basis index"),
      required_string(object, "exact_column_identity")};
  if (result.scc_handle.empty() || result.scc_exact_identity.empty() ||
      result.exact_column_identity.empty())
    throw std::invalid_argument(
        "checkpoint SCC-column provenance contains an empty identity");
  return result;
}
