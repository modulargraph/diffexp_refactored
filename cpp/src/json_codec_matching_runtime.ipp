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
    for (const auto& source : basis_sources_) basis.push_back(source);
    return json::object{
        {"schema", "diffexp2-retained-exact-rational-match-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"basis_sources", std::move(basis)},
        {"incoming_source", incoming_source_},
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
  return left.center_exact == right.center_exact &&
         left.scale_exact == right.scale_exact &&
         left.infinite_radius == right.infinite_radius &&
         (left.infinite_radius || acb_equal(left.radius.raw(), right.radius.raw()));
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
      {"checkpoint_identity", expected_incoming_checkpoint},
      {"analytic_metadata", incoming->exact_analytic_metadata()}};
  if (incoming->column_provenance().has_value())
    incoming_source["column_provenance"] =
        incoming->column_provenance()->encode();
  json::array provenance_basis;
  for (const auto& source : basis_sources)
    provenance_basis.push_back(source);
  json::object provenance{
      {"schema", "diffexp2-native-exact-regular-local-match-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"basis", std::move(provenance_basis)},
      {"incoming", incoming_source},
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
  json::array inconclusive_examples;
  for (const auto& coefficient : diagnostics.coefficients) {
    if (coefficient.verdict == AcbMatchingResidualVerdict::Pass)
      ++pass;
    else if (coefficient.verdict == AcbMatchingResidualVerdict::Fail)
      ++fail;
    else {
      ++inconclusive;
      if (inconclusive_examples.size() < 6)
        inconclusive_examples.push_back(json::object{
            {"row", coefficient.row},
            {"epsilon_power", coefficient.epsilon_power},
            {"residual_upper_exact",
             coefficient.residual_upper.dump_exact()},
            {"scale_lower_exact", coefficient.scale_lower.dump_exact()}});
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
      {"inconclusive_examples", std::move(inconclusive_examples)},
      {"detail", diagnostics.detail}};
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
      EpsilonWindow requested_window,
      std::int32_t required_complete_max, std::uint32_t dimension,
      std::string relative_tolerance, std::size_t max_refinement_steps,
      EpsilonLatticeSaturationResult<Rational>&& exact_saturation,
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
        requested_window_(requested_window),
        required_complete_max_(required_complete_max),
        dimension_(dimension),
        relative_tolerance_(std::move(relative_tolerance)),
        max_refinement_steps_(max_refinement_steps),
        exact_saturation_(std::move(exact_saturation)),
        refined_(std::move(refined)), elapsed_ms_(elapsed_ms) {}

  json::object summary() const override {
    json::array basis;
    basis.reserve(basis_sources_.size());
    for (const auto& source : basis_sources_) basis.push_back(source);

    json::array history;
    history.reserve(refined_.residual_history.size());
    for (std::size_t iteration = 0;
         iteration < refined_.residual_history.size(); ++iteration) {
      auto encoded = encode_acb_match_residual_diagnostics(
          refined_.residual_history[iteration]);
      encoded["iteration"] = iteration;
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
    return json::object{
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
        {"incoming", incoming_source_},
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
                      {"factorizations", 1}}},
        {"weight_windows", std::move(weight_windows)},
        {"transformed_weight_windows",
         std::move(transformed_weight_windows)},
        {"residual", std::move(residual)},
        {"elapsed_ms", elapsed_ms_}};
  }

  const FiniteLaurentVector<ComplexBall>& weights() const {
    return refined_.weights;
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
    for (const auto& source : basis_sources_) basis.push_back(source);
    json::array history;
    history.reserve(refined_.residual_history.size());
    for (const auto& residual : refined_.residual_history)
      history.push_back(checkpoint_acb_match_residual_record(residual));
    return json::object{
        {"schema", "diffexp2-retained-acb-match-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"exact_lattice_identity", exact_lattice_identity_},
        {"exact_lattice_provenance_identity",
         exact_lattice_provenance_identity_},
        {"exact_lattice_canonical_witness",
         exact_lattice_witness_record_},
        {"basis_sources", std::move(basis)},
        {"incoming_source", incoming_source_},
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
        {"refined",
         json::object{
             {"transformed_weights",
              checkpoint_frame_vector_record(refined_.transformed_weights)},
             {"weights", checkpoint_frame_vector_record(refined_.weights)},
             {"residual", checkpoint_frame_vector_record(refined_.residual)},
             {"residual_history", std::move(history)},
             {"refinement_steps", refined_.refinement_steps},
             {"factorization_preconditioner",
              refined_.factorization_preconditioner}}},
        {"elapsed_ms", elapsed_ms_}};
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
  EpsilonWindow requested_window_;
  std::int32_t required_complete_max_ = 0;
  std::uint32_t dimension_ = 0;
  std::string relative_tolerance_;
  std::size_t max_refinement_steps_ = 0;
  EpsilonLatticeSaturationResult<Rational> exact_saturation_;
  RefinedAcbLaurentMatch refined_;
  double elapsed_ms_ = 0.0;
};

struct ParsedExactEvaluatedLattice {
  std::string witness_schema;
  std::string identity;
  std::string canonical_witness;
  EpsilonLatticeSaturationResult<Rational> saturation;
  bool exact_formal_negative_coefficients_are_zero = false;
  std::optional<ExactLaurentMatrix<ComplexBall>> acb_transformation;
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
  for (const auto& source : basis_sources) proof_basis.push_back(source);
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
    expected_basis.push_back(source);
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
  if (compact)
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
  if (!compact && schema != kNativeSingularSCCSaturationRequestSchema)
    throw std::invalid_argument(
        context + ": unsupported native singular-SCC saturation request schema");
  for (const auto* key :
       {"session_configuration_identity", "tile_plan",
        "tile_plan_checkpoint_identity", "match_checkpoint_identity",
        "receiving_scc", "receiving_scc_exact_identity",
        "receiving_basis_point_exact", "physical_match_point_exact"})
    if (required_string(request, key).empty())
      throw std::invalid_argument(
          context + ": native singular-SCC saturation request lost its " +
          key + " binding");
  if (!compact &&
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
  if (expected_request.has_value() &&
      json::serialize(canonical_json_value(request)) !=
          json::serialize(canonical_json_value(*expected_request)))
    throw std::invalid_argument(
        context + ": native singular-SCC saturation request changed its retained plan/SCC binding");
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
    bool prefer_retained_rational_shadow = true) {
  auto native_request = validate_native_singular_scc_saturation_request(
      raw_request, context, expected_session_configuration,
      expected_native_request);
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
  validate_singular_scc_basis_sources(
      native_request, basis_sources, static_cast<std::uint32_t>(dimension),
      context);
  for (std::size_t column = 0; column < dimension; ++column) {
    const auto& provenance = retained_basis[column]->column_provenance();
    if (!provenance.has_value() ||
        json::serialize(canonical_json_value(provenance->encode())) !=
            json::serialize(canonical_json_value(
                basis_sources[column].at("column_provenance"))))
      throw std::invalid_argument(
          context + ": singular-SCC proof source disagrees with its retained column owner");
  }

  std::vector<std::shared_ptr<const RationalShadowColumnWitness>>
      shadow_witnesses;
  shadow_witnesses.reserve(dimension);
  for (const auto& column : retained_basis)
    shadow_witnesses.push_back(column->rational_shadow_witness());
  const bool shadow_basis = std::all_of(
      shadow_witnesses.begin(), shadow_witnesses.end(),
      [](const auto& witness) { return witness != nullptr; });
  if (shadow_basis && prefer_retained_rational_shadow) {
    for (std::size_t column = 0; column < dimension; ++column) {
      const auto& witness = *shadow_witnesses[column];
      if (!witness.solution || witness.rational_shadow_identity.empty() ||
          witness.source_column_identity.empty())
        throw std::invalid_argument(
            context + ": Rational-shadow saturation witness is incomplete");
      const auto parsed_column = json::parse(
          retained_basis[column]->column_provenance()
              ->exact_column_identity);
      const auto& exact_column = as_object(
          parsed_column,
          "Rational-shadow target column identity");
      if (required_string(exact_column, "rational_shadow_identity") !=
              witness.rational_shadow_identity ||
          required_string(exact_column,
                          "rational_source_column_identity") !=
              witness.source_column_identity)
        throw std::invalid_argument(
            context + ": Rational-shadow saturation witness changed its target provenance binding");
      validate_local_solution(*witness.solution, false);
    }
    auto saturation = rational_shadow_formal_saturation(
        shadow_witnesses, window,
        context + ":Rational-shadow-formal-saturation");
    json::array proof_basis;
    for (const auto& source : basis_sources) proof_basis.push_back(source);
    json::array encoded_valuations;
    for (const auto valuation :
         saturation.diagnostics.initial_column_valuations)
      encoded_valuations.push_back(valuation);
    json::object proof_without_identity{
        {"schema", kNativeSingularSCCSaturationProofSchema},
        {"native_request", std::move(native_request)},
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
    return {kNativeSingularSCCSaturationProofSchema, identity,
            json::serialize(canonical_json_value(proof)),
            std::move(saturation), true};
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
        structural_basis, context + ":certified-Acb-Levelt-saturation");
    if (candidate_acb_saturation.diagnostics
                .normalized_determinant_valuation != 0 ||
        candidate_acb_saturation.diagnostics.final_leading_rank !=
            dimension)
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::InvalidSaturationLattice,
          context +
              ": certified Acb Levelt saturation did not produce a valuation-zero full-rank lattice");

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
  } catch (const MatchingArithmeticError&) {
    // Some otherwise well-conditioned small SCCs still have correlated
    // leading balls whose rank is inconclusive after the strict structural
    // chop.  Their direct Laurent factorization remains a fully certified
    // fallback; only the final residual, never this attempted normalization,
    // authorizes materialization.
  }
  const bool monomial_saturation = rational_saturation.has_value() &&
                                   acb_transformation.has_value();
  if (!monomial_saturation) {
    (void)factor_preconditioned_acb_finite_laurent_system(
        evaluated_basis, context + ":certified-direct-Acb-Laurent-pivots");
    json::array proof_basis;
    proof_basis.reserve(basis_sources.size());
    for (const auto& source : basis_sources) proof_basis.push_back(source);
    json::object proof_without_identity{
        {"schema", kNativeSingularSCCSaturationProofSchema},
        {"native_request", std::move(native_request)},
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
    return {kNativeSingularSCCSaturationProofSchema, identity,
            json::serialize(canonical_json_value(proof)),
            unit_rational_saturation(
                static_cast<std::uint32_t>(dimension), window, context)};
  }

  json::array proof_basis;
  proof_basis.reserve(basis_sources.size());
  for (const auto& source : basis_sources) proof_basis.push_back(source);
  json::object proof_without_identity{
      {"schema", kNativeSingularSCCSaturationProofSchema},
      {"native_request", std::move(native_request)},
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
  return {kNativeSingularSCCSaturationProofSchema, identity,
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
  if (required_string(proof, "schema") !=
          kNativeSingularSCCSaturationProofSchema ||
      required_string(proof, "coefficient_domain") != "acb" ||
      as_i32(proof.at("leading_power"),
             "singular-SCC proof leading power") != 0 ||
      as_u32(proof.at("leading_rank"),
             "singular-SCC proof leading rank") != dimension ||
      as_i32(proof.at("determinant_valuation"),
             "singular-SCC determinant valuation") != 0)
    throw std::invalid_argument(
        context + ": native singular-SCC saturation proof facts are inconsistent");
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
  expected_basis.reserve(expected_basis_sources.size());
  for (const auto& source : expected_basis_sources)
    expected_basis.push_back(source);
  if (json::serialize(canonical_json_value(proof.at("basis"))) !=
      json::serialize(canonical_json_value(expected_basis)))
    throw std::invalid_argument(
        context + ": native singular-SCC proof changed its basis/checkpoint binding");
  validate_singular_scc_basis_sources(
      native_request, expected_basis_sources, dimension, context);
  auto identity_input = proof;
  const auto identity = required_string(proof, "identity");
  identity_input.erase("identity");
  if (identity.empty() ||
      json::serialize(canonical_json_value(identity_input)) != identity)
    throw std::invalid_argument(
        context + ": native singular-SCC saturation proof identity is inconsistent");
  if (!monomial_proof)
    return {kNativeSingularSCCSaturationProofSchema, identity,
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
  return {kNativeSingularSCCSaturationProofSchema, identity,
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
    const std::optional<std::size_t>& common_taylor_width) {
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
  if (expected_singular_request.has_value()) {
    const auto receiving_owner = basis.front()->retained_equation_owner();
    if (receiving_owner) {
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
            physical_column);
        if (!transformed.has_value()) {
          available = false;
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
      if (available) {
        auto transformed_incoming =
            receiving_owner->normalize_acb_matching_vector(incoming_value);
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
        matching_basis = std::move(*two_sided_basis);
        matching_incoming = std::move(transformed_incoming->first);
        matching_frame_identity = std::move(*candidate_identity);
        normalized_matching_frame = true;

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
        if (!common_min.has_value() || !common_max.has_value() ||
            *common_min > *common_max ||
            *common_max < required_complete_max) {
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
        } else {
          matching_window = {*common_min, *common_max};
        }
      }
    }
  }
  refinement.required_min_power = matching_window.min_power;

  const auto reset_to_physical_matching_frame = [&]() {
    matching_basis = evaluated_basis;
    matching_incoming = incoming_value;
    matching_window = evaluation_window;
    matching_frame_identity = "physical-parent-frame";
    normalized_matching_frame = false;
    refinement.required_min_power = matching_window.min_power;
  };

  for (std::size_t column = 0; column < basis.size(); ++column) {
    json::object source{
        {"column", column}, {"local", basis_handles[column]},
        {"chart", basis_chart},
        {"source_operator_identity",
         basis[column]->source_operator_identity()},
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
        !normalized_matching_frame);
  };
  auto exact_lattice = make_exact_lattice();
  const auto run_refinement = [&]() {
    return exact_lattice.acb_transformation.has_value()
        ? refine_acb_finite_laurent_match(
              matching_basis, matching_incoming,
              *exact_lattice.acb_transformation, refinement,
              checkpoint_identity + ":refined-acb-match", false)
        : refine_acb_finite_laurent_match(
              matching_basis, matching_incoming, exact_lattice.saturation,
              refinement, checkpoint_identity + ":refined-acb-match",
              exact_lattice.exact_formal_negative_coefficients_are_zero);
  };
  RefinedAcbLaurentMatch refined;
  try {
    refined = run_refinement();
  } catch (const MatchingArithmeticError& error) {
    if (!normalized_matching_frame ||
        error.code != MatchingArithmeticErrorCode::InsufficientCompleteWindow)
      throw;
    // The Fuchsian normal frame is a numerical conditioning optimization,
    // not part of the mathematical matching contract.  Its finite Laurent
    // transforms can consume all private overlap late in a long arm even
    // though the original physical frames still cover the requested result.
    // Fall back transactionally and rebuild the exact lattice in that frame;
    // never turn a conditioning aid into a spurious hard failure.
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
    reset_to_physical_matching_frame();
    exact_lattice = make_exact_lattice();
    refined = run_refinement();
  }

  if (normalized_matching_frame) {
    const auto receiving_owner = basis.front()->retained_equation_owner();
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
      auto physical_residual =
          matching_detail::evaluate_acb_matching_residual(
              evaluated_basis, *physical_weights, incoming_value,
              physical_options,
              checkpoint_identity + ":physical-prefix-check");
      physical_prefix_preserved =
          physical_residual.diagnostics.complete_through_required &&
          physical_residual.diagnostics.verdict ==
              AcbMatchingResidualVerdict::Pass;
      if (physical_prefix_preserved)
        physical_certificate = std::move(physical_residual);
    } catch (const MatchingArithmeticError& error) {
      if (error.code !=
          MatchingArithmeticErrorCode::InsufficientCompleteWindow)
        throw;
    }
    if (!physical_prefix_preserved) {
      reset_to_physical_matching_frame();
      exact_lattice = make_exact_lattice();
      refined = run_refinement();
    } else {
      refined.weights = std::move(*physical_weights);
      refined.residual = std::move(physical_certificate->residual);
      refined.residual_history.back() =
          std::move(physical_certificate->diagnostics);
    }
  }
  const std::string residual_frame_identity = normalized_matching_frame
      ? "physical-parent-frame"
      : matching_frame_identity;

  json::array exact_binding_basis;
  for (const auto& source : basis_sources)
    exact_binding_basis.push_back(source);
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

  json::array provenance_basis;
  for (const auto& source : basis_sources) provenance_basis.push_back(source);
  json::object provenance{
      {"schema", "diffexp2-native-refined-acb-local-match-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"basis", std::move(provenance_basis)},
      {"incoming", incoming_source},
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
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredRefinedAcbMatch>(
      match_handle, checkpoint_identity, provenance_identity,
      exact_lattice.identity, exact_lattice_provenance_identity,
      std::move(exact_lattice.canonical_witness),
      exact_lattice.witness_schema,
      std::move(basis_sources), std::move(incoming_source), basis_chart,
      incoming_chart, basis_point.exact_coordinate,
      incoming_point.exact_coordinate, basis_physical_point,
      matching_frame_identity, residual_frame_identity,
      matching_window,
      required_complete_max, dimension, relative_tolerance,
      max_refinement_steps, std::move(exact_lattice.saturation),
      std::move(refined), elapsed_ms);
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
        std::nullopt) {
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

  const auto retryable_prefix_arithmetic = [](const auto code) {
    return code == MatchingArithmeticErrorCode::AmbiguousZero ||
        code == MatchingArithmeticErrorCode::ZeroDivisor ||
        code == MatchingArithmeticErrorCode::SingularOrIncompleteSystem ||
        code == MatchingArithmeticErrorCode::SaturationFailure ||
        code == MatchingArithmeticErrorCode::SearchBudgetExhausted;
  };
  std::shared_ptr<StoredRefinedAcbMatch> full;
  std::exception_ptr full_arithmetic_failure;
  try {
    full = build_refined_acb_match_once(
        match_handle, request, basis_handles, erased_basis, incoming_handle,
        erased_incoming, precision_bits,
        active_session_configuration_identity, expected_singular_request,
        std::nullopt);
  } catch (const MatchingArithmeticError& error) {
    if (!retryable_prefix_arithmetic(error.code)) throw;
    full_arithmetic_failure = std::current_exception();
  }
  if (full && full->certified_for_materialization()) {
    full->replace_elapsed_ms(std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count());
    return full;
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
      (full->final_residual_verdict() !=
           AcbMatchingResidualVerdict::Inconclusive ||
       !full->residual_complete_through_required())) {
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
          width);
    } catch (const MatchingArithmeticError& error) {
      if (!retryable_prefix_arithmetic(error.code)) throw;
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
