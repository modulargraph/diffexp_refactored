const char* line_integration_scope_name(LineIntegrationScope scope) {
  switch (scope) {
    case LineIntegrationScope::StoredTruncation:
      return "stored_truncation";
    case LineIntegrationScope::FullLocalWithCertifiedTail:
      return "full_local_with_certified_tail";
  }
  throw std::logic_error("unknown line integration scope");
}

const char* error_guarantee_name(ErrorGuarantee guarantee) {
  switch (guarantee) {
    case ErrorGuarantee::None:
      return "none";
    case ErrorGuarantee::Advisory:
      return "advisory";
    case ErrorGuarantee::Certified:
      return "certified";
  }
  throw std::logic_error("unknown error guarantee");
}

json::object encode_error_envelope_summary(const ErrorEnvelope& error) {
  json::object result{
      {"guarantee", error_guarantee_name(error.guarantee)},
      {"provenance", error.provenance}};
  if (!error.empty()) {
    json::array upper;
    upper.reserve(error.absolute.size());
    for (const auto& bound : error.absolute)
      upper.push_back(bound.approximate_upper());
    result["epsilon_min"] = error.frame.min_power;
    result["epsilon_max"] = error.frame.complete_max;
    result["absolute_upper_approx"] = std::move(upper);
    result["bound_encoding"] = "approximate-double-diagnostics";
  }
  return result;
}

void diagnose_plain_regular_taylor_evaluation(
    const LocalSolution<ComplexBall>& solution,
    const RealEvaluationPoint& point,
    const LocalEvaluation& evaluation,
    const std::string& context) {
  const auto* requested =
      std::getenv("DIFFEXP_DIAGNOSTIC_TERMINAL_INCOMING_TAYLOR");
  if (requested == nullptr) return;
  const std::string request(requested);
  const auto separator = request.find(':');
  if (separator == std::string::npos ||
      request.find(':', separator + 1) != std::string::npos)
    throw std::invalid_argument(
        "DIFFEXP_DIAGNOSTIC_TERMINAL_INCOMING_TAYLOR must be "
        "<epsilon-power>:<component>");
  const auto epsilon_power =
      static_cast<std::int32_t>(std::stoll(request.substr(0, separator)));
  const auto component_wide =
      std::stoull(request.substr(separator + 1));
  if (component_wide >= solution.dimension)
    throw std::invalid_argument(
        "DIFFEXP_DIAGNOSTIC_TERMINAL_INCOMING_TAYLOR component is outside "
        "the incoming local");
  const auto component = static_cast<std::uint32_t>(component_wide);
  if (epsilon_power < solution.epsilon.min_power ||
      epsilon_power > solution.epsilon.complete_max ||
      epsilon_power < evaluation.value.epsilon.min_power ||
      epsilon_power > evaluation.value.epsilon.complete_max) {
    std::cerr
        << "terminal-incoming-taylor-summary context=" << context
        << " status=outside-retained-window"
        << " epsilon_power=" << epsilon_power
        << " component=" << component << '\n';
    return;
  }

  for (const auto& sector : solution.sectors)
    if (sector.a.is_zero != TruthValue::Yes ||
        sector.b.is_zero != TruthValue::Yes ||
        sector.log_power != 0) {
      std::cerr
          << "terminal-incoming-taylor-summary context=" << context
          << " status=unsupported-nonordinary-sector"
          << " epsilon_power=" << epsilon_power
          << " component=" << component << '\n';
      return;
    }

  ComplexBall signed_t = point.modulus;
  if (point.sign < 0) signed_t = -signed_t;
  std::vector<ComplexBall> t_powers(solution.taylor_width(),
                                    ComplexBall(1));
  for (std::size_t taylor = 1;
       taylor < solution.taylor_width(); ++taylor)
    t_powers[taylor] = t_powers[taylor - 1] * signed_t;

  const auto epsilon_index = static_cast<std::size_t>(
      epsilon_power - solution.epsilon.min_power);
  ComplexBall recomposed(0);
  auto summed_uncertainty = Magnitude::zero();
  auto dominant_uncertainty = Magnitude::zero();
  std::optional<std::size_t> dominant_sector;
  std::optional<std::size_t> dominant_taylor;
  std::optional<ComplexBall> dominant_coefficient;
  std::optional<ComplexBall> dominant_term;
  for (std::size_t sector_index = 0;
       sector_index < solution.sectors.size(); ++sector_index) {
    const auto& sector = solution.sectors[sector_index];
    for (std::size_t taylor = 0;
         taylor < solution.taylor_width(); ++taylor) {
      const auto& coefficient = sector.coefficients[
          local_detail::sector_index(
              solution, epsilon_index, taylor, component)];
      const auto term = coefficient * t_powers[taylor];
      recomposed += term;
      const auto uncertainty = Magnitude::upper_abs(
          term - matching_detail::acb_midpoint_value(term));
      summed_uncertainty += uncertainty;
      if (!dominant_taylor.has_value() ||
          !(uncertainty <= dominant_uncertainty)) {
        dominant_uncertainty = uncertainty;
        dominant_sector = sector_index;
        dominant_taylor = taylor;
        dominant_coefficient = coefficient;
        dominant_term = term;
      }
      std::cerr
          << "terminal-incoming-taylor-term context=" << context
          << " epsilon_power=" << epsilon_power
          << " component=" << component
          << " sector=" << sector_index
          << " taylor=" << taylor
          << " coefficient_midpoint=("
          << coefficient.real_midpoint(12) << ","
          << coefficient.imag_midpoint(12) << ")"
          << " coefficient_radius2exp=("
          << coefficient.real_radius_exponent() << ","
          << coefficient.imag_radius_exponent() << ")"
          << " term_midpoint=("
          << term.real_midpoint(12) << ","
          << term.imag_midpoint(12) << ")"
          << " term_radius2exp=("
          << term.real_radius_exponent() << ","
          << term.imag_radius_exponent() << ")"
          << " uncertainty_upper="
          << uncertainty.approximate_upper()
          << '\n';
    }
  }
  if (!dominant_sector.has_value() || !dominant_taylor.has_value() ||
      !dominant_coefficient.has_value() || !dominant_term.has_value())
    throw std::logic_error(
        "terminal incoming Taylor diagnostic found no retained term");
  const auto& authoritative =
      evaluation.value.at(epsilon_power, component);
  std::cerr
      << "terminal-incoming-taylor-summary context=" << context
      << " status=ordinary"
      << " epsilon_power=" << epsilon_power
      << " component=" << component
      << " point_exact=" << point.exact_coordinate
      << " retained_taylor_order="
      << solution.taylor_complete_max
      << " dominant_sector=" << *dominant_sector
      << " dominant_taylor=" << *dominant_taylor
      << " dominant_coefficient_radius2exp=("
      << dominant_coefficient->real_radius_exponent() << ","
      << dominant_coefficient->imag_radius_exponent() << ")"
      << " dominant_term_radius2exp=("
      << dominant_term->real_radius_exponent() << ","
      << dominant_term->imag_radius_exponent() << ")"
      << " recomposed_radius2exp=("
      << recomposed.real_radius_exponent() << ","
      << recomposed.imag_radius_exponent() << ")"
      << " authoritative_radius2exp=("
      << authoritative.real_radius_exponent() << ","
      << authoritative.imag_radius_exponent() << ")"
      << " dominant_uncertainty_upper="
      << dominant_uncertainty.approximate_upper()
      << " summed_uncertainty_upper="
      << summed_uncertainty.approximate_upper()
      << '\n';
}

EpsilonVector inflate_certified_physical_evaluation(
    const CertifiedLocalEvaluation& certified,
    const EpsilonWindow& required_frame,
    const std::string& context) {
  if (certified.tail.status != TailMajorantStatus::Certified ||
      certified.tail.value.guarantee != ErrorGuarantee::Certified ||
      !tail_majorant_detail::same_epsilon_window(
          certified.tail.value.frame, required_frame) ||
      certified.tail.value.absolute.size() != required_frame.width())
    throw std::logic_error(
        context + ": certified physical evaluation has the wrong tail frame");
  auto value = certified.evaluation.value;
  value.error = ErrorEnvelope{};
  if (value.epsilon.complete_max != required_frame.complete_max ||
      value.epsilon.min_power < required_frame.min_power)
    throw std::logic_error(
        context +
        ": certified physical evaluation changed its retained epsilon frame");
  if (required_frame.min_power < value.epsilon.min_power) {
    auto widened = physical_ode_detail::zero_epsilon_vector(
        required_frame, value.dimension);
    for (std::int64_t raw_power = value.epsilon.min_power;
         raw_power <= value.epsilon.complete_max; ++raw_power) {
      const auto power = static_cast<std::int32_t>(raw_power);
      for (std::uint32_t component = 0;
           component < value.dimension; ++component)
        widened.at(power, component) = value.at(power, component);
    }
    value = std::move(widened);
  }
  for (std::int64_t raw_power = required_frame.min_power;
       raw_power <= required_frame.complete_max; ++raw_power) {
    const auto power = static_cast<std::int32_t>(raw_power);
    const auto row = static_cast<std::size_t>(
        raw_power - required_frame.min_power);
    for (std::uint32_t component = 0;
         component < value.dimension; ++component)
      certified.tail.value.absolute[row].add_error_to(
          value.at(power, component));
  }
  return value;
}

std::optional<LocalEvaluation>
diagnose_two_step_terminal_physical_recenter(
    const LocalSolution<ComplexBall>& source,
    const std::shared_ptr<
        const PreparedPhysicalClearedODE<ComplexBall>>& equation,
    const RealEvaluationPoint& target,
    const EvaluationOptions& options,
    const std::string& context) {
  const auto* mode =
      std::getenv("DIFFEXP_DIAGNOSTIC_TERMINAL_PHYSICAL_RECENTER");
  if (std::getenv("DIFFEXP_DIAGNOSTIC_TERMINAL_INCOMING_TAYLOR") != nullptr)
    std::cerr
        << "terminal-physical-recenter-env context=" << context
        << " value=" << (mode == nullptr ? "unset" : mode) << '\n';
  if (mode == nullptr) return std::nullopt;
  if (std::string(mode) != "report")
    throw std::invalid_argument(
        "DIFFEXP_DIAGNOSTIC_TERMINAL_PHYSICAL_RECENTER must be report");
  const auto fail = [&](const std::string& detail)
      -> std::optional<LocalEvaluation> {
    std::cerr
        << "terminal-physical-recenter context=" << context
        << " status=inconclusive detail=" << detail << '\n';
    return std::nullopt;
  };
  if (!equation)
    return fail("no-retained-physical-equation");
  if (!target.certified_algebraic && target.exact_coordinate.empty())
    return fail("target-has-no-exact-coordinate");
  const Rational target_exact(target.exact_coordinate);
  if (target_exact.is_zero())
    return fail("target-is-already-the-source-center");
  if (source.chart.infinite_radius || source.chart.radius_exact.empty())
    return fail("source-chart-has-no-finite-exact-radius");
  const Rational source_radius(source.chart.radius_exact);
  const Rational intermediate = target_exact / Rational(2);
  const auto intermediate_modulus =
      exact_path_detail::abs(intermediate);
  if (!(intermediate_modulus < source_radius))
    return fail("intermediate-point-is-outside-source-chart");

  auto source_prepared =
      prepare_physical_regular_homogeneous_tail_model(
          *equation, source);
  if (source_prepared.status != TailMajorantStatus::Certified ||
      !source_prepared.model.has_value())
    return fail(
        "source-tail-model-" +
        std::string(tail_majorant_status_name(source_prepared.status)) +
        ":" + source_prepared.detail);
  const auto find_witness = [](
      const PhysicalRegularTaylorTailModel& model,
      const RealEvaluationPoint& point,
      const Rational& point_modulus,
      const Rational& radius,
      const EvaluationOptions& evaluation_options)
      -> std::optional<CertifiedLocalEvaluation> {
    Rational denominator(1);
    for (std::uint32_t exponent = 1; exponent <= 16; ++exponent) {
      denominator *= Rational(2);
      const auto witness =
          point_modulus + (radius - point_modulus) / denominator;
      auto evaluated =
          evaluate_physical_local_solution_with_certified_tail(
              model, point, witness.str(), evaluation_options);
      if (evaluated.tail.status == TailMajorantStatus::Certified)
        return evaluated;
    }
    return std::nullopt;
  };
  const auto intermediate_point =
      RealEvaluationPoint::rational(intermediate.str());
  auto source_certified = find_witness(
      *source_prepared.model, intermediate_point,
      intermediate_modulus, source_radius, options);
  if (!source_certified.has_value())
    return fail("source-to-intermediate-has-no-certified-dyadic-witness");
  auto intermediate_value = inflate_certified_physical_evaluation(
      *source_certified, source_prepared.model->epsilon,
      context + ": source-to-intermediate");

  auto recentered =
      recenter_physical_cleared_ode(*equation, intermediate);
  if (!recentered.eligible || !recentered.equation.has_value())
    return fail("recentered-equation-ineligible:" + recentered.reason);
  const auto recentered_radius =
      source_radius - intermediate_modulus;
  const auto return_exact = target_exact - intermediate;
  const auto return_modulus =
      exact_path_detail::abs(return_exact);
  if (!(return_modulus < recentered_radius))
    return fail("target-is-outside-conservative-recentered-chart");
  const auto evolved = evolve_ordinary_center_value<ComplexBall>(
      *recentered.equation, intermediate_value,
      source.taylor_complete_max);
  if (!evolved.eligible)
    return fail("recentered-evolution-ineligible:" + evolved.reason);
  ChartGeometry recentered_chart;
  recentered_chart.center_exact = "0";
  recentered_chart.scale_exact = "1";
  recentered_chart.radius_exact = recentered_radius.str();
  recentered_chart.radius =
      ComplexBall::from_strings(recentered_radius.str());
  recentered_chart.infinite_radius = false;
  auto recentered_solution = ordinary_evolution_local_solution(
      evolved, std::move(recentered_chart), source.prescriptions,
      source.checkpoint_identity +
          ":diagnostic-terminal-physical-half-recenter");
  auto recentered_prepared =
      prepare_physical_regular_homogeneous_tail_model(
          *recentered.equation, recentered_solution);
  if (recentered_prepared.status != TailMajorantStatus::Certified ||
      !recentered_prepared.model.has_value())
    return fail(
        "recentered-tail-model-" +
        std::string(tail_majorant_status_name(
            recentered_prepared.status)) +
        ":" + recentered_prepared.detail);
  EvaluationOptions recentered_options = options;
  recentered_options.imaginary_sign.reset();
  const auto return_point =
      RealEvaluationPoint::rational(return_exact.str());
  auto return_certified = find_witness(
      *recentered_prepared.model, return_point,
      return_modulus, recentered_radius, recentered_options);
  if (!return_certified.has_value())
    return fail("recentered-return-has-no-certified-dyadic-witness");
  auto final_value = inflate_certified_physical_evaluation(
      *return_certified, recentered_prepared.model->epsilon,
      context + ": recentered-return");
  auto result = return_certified->evaluation;
  result.value = std::move(final_value);
  result.value.error = ErrorEnvelope{};

  double source_tail_max = 0.0;
  for (const auto& bound :
       source_certified->tail.value.absolute)
    source_tail_max =
        std::max(source_tail_max, bound.approximate_upper());
  double return_tail_max = 0.0;
  for (const auto& bound :
       return_certified->tail.value.absolute)
    return_tail_max =
        std::max(return_tail_max, bound.approximate_upper());
  std::cerr
      << "terminal-physical-recenter context=" << context
      << " status=certified"
      << " target_exact=" << target_exact.str()
      << " intermediate_exact=" << intermediate.str()
      << " recentered_return_exact=" << return_exact.str()
      << " source_taylor_order="
      << source_prepared.model->taylor_complete_max
      << " recentered_taylor_order="
      << recentered_prepared.model->taylor_complete_max
      << " source_tail_upper_max=" << source_tail_max
      << " return_tail_upper_max=" << return_tail_max
      << '\n';
  return result;
}

std::optional<LocalEvaluation>
diagnose_factorized_terminal_physical_evaluation(
    const LocalSolution<ComplexBall>& source,
    const std::shared_ptr<
        const PreparedPhysicalClearedODE<ComplexBall>>& equation,
    const RealEvaluationPoint& target,
    const EvaluationOptions& options,
    const std::string& context) {
  const auto* mode =
      std::getenv("DIFFEXP_DIAGNOSTIC_TERMINAL_FACTORIZED_EVALUATION");
  if (mode == nullptr) return std::nullopt;
  if (std::string(mode) != "report")
    throw std::invalid_argument(
        "DIFFEXP_DIAGNOSTIC_TERMINAL_FACTORIZED_EVALUATION must be report");
  const auto fail = [&](const std::string& detail)
      -> std::optional<LocalEvaluation> {
    std::cerr
        << "terminal-factorized-physical-evaluation context=" << context
        << " status=inconclusive detail=" << detail << '\n';
    return std::nullopt;
  };
  if (!equation)
    return fail("no-retained-physical-equation");
  const auto started = std::chrono::steady_clock::now();
  auto factorized = evaluate_ordinary_center_value_factorized(
      *equation, source, target, options);
  if (!factorized.eligible)
    return fail(factorized.reason);
  const auto elapsed_ms =
      std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
  std::cerr
      << "terminal-factorized-physical-evaluation context=" << context
      << " status=complete"
      << " operator_columns=" << factorized.operator_columns
      << " retained_taylor_order="
      << source.taylor_complete_max
      << " elapsed_ms=" << elapsed_ms
      << '\n';
  return factorized.evaluation;
}

LocalEvaluation terminal_incoming_evaluation_with_factorized_fallback(
    const LocalSolution<ComplexBall>& source,
    const std::shared_ptr<
        const PreparedPhysicalClearedODE<ComplexBall>>& equation,
    const RealEvaluationPoint& point,
    const EvaluationOptions& options,
    LocalEvaluation direct,
    const std::string& context) {
  return certified_ordinary_center_evaluation_with_factorized_fallback(
      source, equation, point, options, std::move(direct), context);
}

struct RetainedPlanChartBinding {
  using Owner = std::variant<std::shared_ptr<PreparedChartBase>,
                             std::shared_ptr<CompositeSCCChartBase>,
                             std::shared_ptr<RegularPhysicalEquationOwnerBase>>;

  std::string handle;
  std::string exact_identity;
  // Rational surrogate used only by the exact combinatorial tile planner.
  ExactAffineChart geometry;
  // Prepared equation geometry.  It may be genuinely algebraic; its exact
  // strings are identities and the balls are rigorous specializations.
  ChartGeometry local_geometry;
  ComplexBall center_numeric;
  ComplexBall scale_numeric;
  std::int32_t scale_sign = 0;
  std::string planning_certificate_identity;
  std::vector<Prescription> prescriptions;
  Owner owner;
};

struct CertifiedPlanPoint {
  RealEvaluationPoint local;
  std::string physical_exact;
  ComplexBall physical_numeric;
  // Preserve the producer's rigorous two-string Acb request encoding.  The
  // parsed ball is used for validation/integration; this original encoding
  // is replayed when a retained planned match invokes the matching kernel.
  json::value local_numeric_encoding;
};

struct CertifiedPlanMatch {
  CertifiedPlanPoint producing;
  CertifiedPlanPoint receiving;
};

struct CertifiedPlanTile {
  std::string physical_begin_exact;
  std::string physical_end_exact;
  ComplexBall physical_begin_numeric;
  ComplexBall physical_end_numeric;
  RealEvaluationPoint local_begin;
  RealEvaluationPoint local_end;
};

struct RetainedArmPlan {
  ExactArmPlan exact;
  std::vector<RetainedPlanChartBinding> charts;
  std::string certified_geometry_identity;
  std::vector<CertifiedPlanMatch> certified_matches;
  std::vector<CertifiedPlanTile> certified_tiles;
  // These canonical request records are part of the retained/checkpointed
  // plan identity.  They let restore rebind the rational planning surrogate
  // and the Wolfram-certified algebraic geometry to the exact same owners.
  json::array planning_charts;
  json::object certified_geometry;
};

json::array encode_path_branch_sheets(
    const std::vector<ExactBranchSheet>& sheets) {
  json::array encoded;
  encoded.reserve(sheets.size());
  for (const auto& sheet : sheets)
    encoded.push_back(json::object{{"factor_exact", sheet.factor_exact},
                                   {"sign", sheet.imaginary_sign}});
  return encoded;
}

json::array encode_plan_prescriptions(
    const std::vector<Prescription>& prescriptions) {
  json::array encoded;
  encoded.reserve(prescriptions.size());
  for (const auto& prescription : prescriptions)
    encoded.push_back(json::object{
        {"factor_exact", prescription.factor_exact},
        {"sign", prescription.sign},
        {"multiplicity", prescription.multiplicity},
        {"leading_coefficient_sign",
         prescription.leading_coefficient_sign}});
  return encoded;
}

json::object encode_plan_chart(const RetainedPlanChartBinding& binding,
                               std::size_t index) {
  return json::object{
      {"index", index}, {"chart", binding.handle},
      {"identity", binding.exact_identity},
      {"center_exact", binding.geometry.center.str()},
      {"scale_exact", binding.geometry.scale.str()},
      {"radius_exact", binding.geometry.radius.str()},
      {"singular_center", binding.geometry.singular_center},
      {"prescriptions", encode_plan_prescriptions(binding.prescriptions)}};
}

const char* exact_match_kind_name(ExactMatchKind kind) {
  switch (kind) {
    case ExactMatchKind::SymmetricDivisionPoint:
      return "symmetric-division-point";
    case ExactMatchKind::BalancedSafeOverlap:
      return "balanced-safe-overlap";
    case ExactMatchKind::SingularBalancedApproach:
      return "singular-balanced-approach";
    case ExactMatchKind::ForbiddenPointAvoidance:
      return "forbidden-point-avoidance";
  }
  throw std::logic_error("unknown exact match kind");
}

json::object encode_plan_match(const RetainedArmPlan& arm,
                               std::size_t index) {
  if (index >= arm.exact.matches.size())
    throw std::invalid_argument("native tile-plan match index is out of range");
  const auto& match = arm.exact.matches[index];
  return json::object{
      {"index", index},
      {"producing_chart_index", match.producing_chart},
      {"receiving_chart_index", match.receiving_chart},
      {"producing_chart", arm.charts.at(match.producing_chart).handle},
      {"receiving_chart", arm.charts.at(match.receiving_chart).handle},
      {"physical_exact", match.physical.str()},
      {"producing_local_exact", match.producing_local.str()},
      {"receiving_local_exact", match.receiving_local.str()},
      {"kind", exact_match_kind_name(match.kind)},
      {"branch_sheets", encode_path_branch_sheets(match.branch_sheets)}};
}

json::object encode_plan_tile(const RetainedArmPlan& arm,
                              std::size_t index) {
  if (index >= arm.exact.tiles.size())
    throw std::invalid_argument("native tile-plan tile index is out of range");
  const auto& tile = arm.exact.tiles[index];
  const auto& binding = arm.charts.at(tile.chart);
  return json::object{
      {"index", index}, {"chart_index", tile.chart},
      {"chart", binding.handle}, {"chart_identity", binding.exact_identity},
      {"physical_begin_exact", tile.physical_begin.str()},
      {"physical_end_exact", tile.physical_end.str()},
      {"local_begin_exact", tile.local_begin.str()},
      {"local_end_exact", tile.local_end.str()},
      {"jacobian_exact", binding.geometry.scale.str()},
      {"crosses_singular_center", tile.crosses_singular_center},
      {"branch_sheets", encode_path_branch_sheets(tile.branch_sheets)},
      {"prescriptions", encode_plan_prescriptions(binding.prescriptions)}};
}

json::object encode_plan_topology(const ExactPathTopology& topology) {
  json::array singular_points;
  for (const auto& point : topology.singular_points)
    singular_points.push_back(json::value(point.str()));
  json::array boundary_points;
  for (const auto& point : topology.boundary_points)
    boundary_points.push_back(json::value(point.str()));
  json::array projections;
  for (const auto& projection : topology.complex_projections)
    projections.push_back(json::object{
        {"source_identity", projection.source_identity},
        {"real_part_exact", projection.real_part.str()},
        {"imaginary_magnitude_exact", projection.imaginary_magnitude.str()},
        {"retain_minus_imaginary", projection.retain_minus_imaginary},
        {"retain_real_part", projection.retain_real_part},
        {"retain_plus_imaginary", projection.retain_plus_imaginary}});
  return json::object{
      {"singular_points", std::move(singular_points)},
      {"boundary_points", std::move(boundary_points)},
      {"complex_projections", std::move(projections)},
      {"branch_sheets", encode_path_branch_sheets(topology.branch_sheets)}};
}

json::object encode_retained_arm(const RetainedArmPlan& arm) {
  json::array charts;
  for (std::size_t index = 0; index < arm.charts.size(); ++index)
    charts.push_back(encode_plan_chart(arm.charts[index], index));
  json::array matches;
  for (std::size_t index = 0; index < arm.exact.matches.size(); ++index)
    matches.push_back(encode_plan_match(arm, index));
  json::array tiles;
  for (std::size_t index = 0; index < arm.exact.tiles.size(); ++index)
    tiles.push_back(encode_plan_tile(arm, index));
  return json::object{
      {"from_exact", arm.exact.from.str()},
      {"to_exact", arm.exact.to.str()},
      {"direction", arm.exact.direction},
      {"division_order", arm.exact.division_order},
      {"charts", std::move(charts)}, {"matches", std::move(matches)},
      {"tiles", std::move(tiles)},
      {"planning_charts", arm.planning_charts},
      {"certified_geometry", arm.certified_geometry},
      {"topology", encode_plan_topology(arm.exact.topology)}};
}

// The hot tile.plan response is an opaque retained-handle DTO.  Publishing
// exact algebraic interval strings here made response size depend on the
// complexity of Root objects (and duplicated data that the server already
// owns).  Downstream transport commands consume the retained plan handle;
// exact intervals remain available through the diagnostic tile.stats and
// pointwise tile.*_interval routes.  This summary is therefore bounded by
// scalar counters and server-generated handles, not mathematical payload
// size.
json::object encode_published_arm_summary(const RetainedArmPlan& arm,
                                          const std::string& name) {
  return json::object{
      {"arm_name", name},
      {"direction", arm.exact.direction},
      {"chart_count", arm.charts.size()},
      {"match_count", arm.exact.matches.size()},
      {"tile_count", arm.exact.tiles.size()},
      {"exact_intervals_published", false},
      {"detail", "opaque-retained-plan-arm"}};
}

std::optional<std::int32_t> exact_plan_rim(
    const std::vector<Prescription>& prescriptions,
    const Rational& chart_scale) {
  if (chart_scale.is_zero())
    throw std::invalid_argument(
        "prepared tile chart has a zero exact scale");
  const auto scale_sign = chart_scale.sign();
  std::optional<std::int32_t> rim;
  for (const auto& prescription : prescriptions) {
    if ((prescription.multiplicity & 1U) == 0)
      throw std::invalid_argument(
          "prepared tile chart has an even-multiplicity tangential "
          "prescription; a one-sided real-axis rim is not defined");
    const auto candidate =
        prescription.sign * prescription.leading_coefficient_sign *
        scale_sign;
    if (rim.has_value() && *rim != candidate)
      throw std::invalid_argument(
          "prepared tile chart has conflicting exact odd-multiplicity prescriptions");
    rim = candidate;
  }
  return rim;
}

class StoredTilePlan {
 public:
  StoredTilePlan(std::string handle, std::string checkpoint_identity,
                 std::string provenance_identity, std::uint32_t division_order,
                 RetainedArmPlan lower, RetainedArmPlan upper,
                 double elapsed_ms)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        division_order_(division_order), lower_(std::move(lower)),
        upper_(std::move(upper)), elapsed_ms_(elapsed_ms) {
    validate_arm_set();
  }

  StoredTilePlan(std::string handle, std::string checkpoint_identity,
                 std::string provenance_identity, std::uint32_t division_order,
                 std::string arm_name, RetainedArmPlan arm,
                 double elapsed_ms)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        division_order_(division_order), elapsed_ms_(elapsed_ms) {
    if (arm_name == "lower")
      lower_ = std::move(arm);
    else if (arm_name == "upper")
      upper_ = std::move(arm);
    else
      throw std::invalid_argument(
          "single-arm tile plan name must be lower or upper");
    validate_arm_set();
  }

  const std::string& handle() const { return handle_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& provenance_identity() const {
    return provenance_identity_;
  }

  bool has_arm(const std::string& name) const {
    if (name == "lower") return lower_.has_value();
    if (name == "upper") return upper_.has_value();
    return false;
  }

  bool has_two_arms() const {
    return lower_.has_value() && upper_.has_value();
  }

  const RetainedArmPlan& arm(const std::string& name) const {
    if (name == "lower" && lower_.has_value()) return *lower_;
    if (name == "upper" && upper_.has_value()) return *upper_;
    if (name != "lower" && name != "upper")
      throw std::invalid_argument(
          "native tile-plan arm must be lower or upper");
    throw std::invalid_argument(
        "native tile plan does not retain the requested " + name + " arm");
  }

  json::object match_interval(const std::string& name,
                              std::size_t index) const {
    match_queries_.fetch_add(1);
    auto encoded = encode_plan_match(arm(name), index);
    encoded["arm"] = name;
    encoded["tile_plan"] = handle_;
    encoded["checkpoint_identity"] = checkpoint_identity_;
    return encoded;
  }

  json::object tile_interval(const std::string& name,
                             std::size_t index) const {
    tile_queries_.fetch_add(1);
    auto encoded = encode_plan_tile(arm(name), index);
    encoded["arm"] = name;
    encoded["tile_plan"] = handle_;
    encoded["checkpoint_identity"] = checkpoint_identity_;
    return encoded;
  }

  void note_integration() { integrations_.fetch_add(1); }

  void note_match_advance(const std::string& name) {
    (void)arm(name);
    if (name == "lower") {
      lower_match_advances_.fetch_add(1);
      return;
    }
    if (name == "upper") {
      upper_match_advances_.fetch_add(1);
      return;
    }
    throw std::invalid_argument("native tile-plan arm must be lower or upper");
  }

  void require_two_arm_match_advance_capacity(
      std::size_t lower, std::size_t upper) const {
    if (!has_two_arms())
      throw std::invalid_argument(
          "parallel match accounting requires a retained two-arm plan");
    const auto fits = [](std::uint64_t current, std::size_t increment) {
      return increment <= std::numeric_limits<std::uint64_t>::max() - current;
    };
    if (!fits(lower_match_advances_.load(), lower) ||
        !fits(upper_match_advances_.load(), upper))
      throw std::overflow_error(
          "retained tile-plan match-advance counter overflow");
  }

  void note_two_arm_match_advances(
      std::size_t lower, std::size_t upper) noexcept {
    lower_match_advances_.fetch_add(static_cast<std::uint64_t>(lower));
    upper_match_advances_.fetch_add(static_cast<std::uint64_t>(upper));
  }

  std::vector<std::shared_ptr<PreparedChartBase>> dependency_charts() const {
    std::vector<std::shared_ptr<PreparedChartBase>> result;
    result.reserve((lower_.has_value() ? lower_->charts.size() : 0) +
                   (upper_.has_value() ? upper_->charts.size() : 0));
    const auto append = [&](const RetainedPlanChartBinding& binding) {
      std::visit(
          [&](const auto& owner) {
            using Owner = typename std::decay_t<decltype(owner)>::element_type;
            if constexpr (std::is_same_v<Owner, PreparedChartBase>) {
              result.push_back(owner);
            } else if constexpr (
                std::is_same_v<Owner, CompositeSCCChartBase>) {
              auto dependencies = owner->dependency_charts();
              result.insert(result.end(), dependencies.begin(),
                            dependencies.end());
            }
          },
          binding.owner);
    };
    if (lower_.has_value())
      for (const auto& binding : lower_->charts) append(binding);
    if (upper_.has_value())
      for (const auto& binding : upper_->charts) append(binding);
    return result;
  }

  std::vector<std::shared_ptr<CompositeSCCChartBase>> dependency_sccs() const {
    std::vector<std::shared_ptr<CompositeSCCChartBase>> result;
    const auto append = [&](const RetainedPlanChartBinding& binding) {
      if (const auto* owner = std::get_if<
              std::shared_ptr<CompositeSCCChartBase>>(&binding.owner))
        result.push_back(*owner);
    };
    if (lower_.has_value())
      for (const auto& binding : lower_->charts) append(binding);
    if (upper_.has_value())
      for (const auto& binding : upper_->charts) append(binding);
    return result;
  }

  std::vector<std::shared_ptr<RegularPhysicalEquationOwnerBase>>
  dependency_regular_equation_owners() const {
    std::vector<std::shared_ptr<RegularPhysicalEquationOwnerBase>> result;
    const auto append = [&](const RetainedPlanChartBinding& binding) {
      if (const auto* owner = std::get_if<
              std::shared_ptr<RegularPhysicalEquationOwnerBase>>(
              &binding.owner))
        result.push_back(*owner);
    };
    if (lower_.has_value())
      for (const auto& binding : lower_->charts) append(binding);
    if (upper_.has_value())
      for (const auto& binding : upper_->charts) append(binding);
    return result;
  }

  json::object summary(bool include_intervals = true) const {
    if (!has_two_arms()) return single_arm_summary(include_intervals);
    json::object result{
        {"tile_plan", handle_}, {"capability", kRetainedTilePlanCapability},
        {"native_retained", true}, {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"division_order", division_order_},
        {"independent_arms", true},
        {"concurrent_execution", "immutable-independent-arm-snapshots"},
        {"anchor", json::object{
             {"lower_chart", lower_->charts.front().handle},
             {"upper_chart", upper_->charts.front().handle},
             {"center_exact", lower_->exact.from.str()}}},
        {"lower_matches", lower_->exact.matches.size()},
        {"upper_matches", upper_->exact.matches.size()},
        {"lower_tiles", lower_->exact.tiles.size()},
        {"upper_tiles", upper_->exact.tiles.size()},
        {"match_interval_queries", match_queries_.load()},
        {"tile_interval_queries", tile_queries_.load()},
        {"lower_match_advances", lower_match_advances_.load()},
        {"upper_match_advances", upper_match_advances_.load()},
        {"integrations", integrations_.load()},
        {"elapsed_ms", elapsed_ms_}};
    if (include_intervals) {
      result["lower"] = encode_retained_arm(*lower_);
      result["upper"] = encode_retained_arm(*upper_);
    }
    return result;
  }

  json::object publication_summary() const {
    if (!has_two_arms()) {
      const auto name = single_arm_name();
      const auto& retained = arm(name);
      json::object result{
          {"tile_plan", handle_},
          {"capability", kRetainedSingleArmTilePlanCapability},
          {"native_retained", true},
          {"checkpoint_identity", checkpoint_identity_},
          {"provenance_reference",
           compact_matching_identity_reference(provenance_identity_)},
          {"division_order", division_order_},
          {"independent_arms", false},
          {"concurrent_execution", "single-immutable-arm-snapshot"},
          {"arm_name", name},
          {"anchor", json::object{{"chart", retained.charts.front().handle}}},
          {"matches", retained.exact.matches.size()},
          {"tiles", retained.exact.tiles.size()},
          {"match_interval_queries", match_queries_.load()},
          {"tile_interval_queries", tile_queries_.load()},
          {"lower_match_advances", lower_match_advances_.load()},
          {"upper_match_advances", upper_match_advances_.load()},
          {"integrations", integrations_.load()},
          {"elapsed_ms", elapsed_ms_}};
      auto published = encode_published_arm_summary(retained, name);
      result["arm"] = published;
      result[name] = std::move(published);
      result["diagnostic_detail"] = "bounded-opaque-plan-publication";
      return result;
    }
    json::object result{
        {"tile_plan", handle_},
        {"capability", kRetainedTilePlanCapability},
        {"native_retained", true},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_reference",
         compact_matching_identity_reference(provenance_identity_)},
        {"division_order", division_order_},
        {"independent_arms", true},
        {"concurrent_execution", "immutable-independent-arm-snapshots"},
        {"anchor", json::object{
             {"lower_chart", lower_->charts.front().handle},
             {"upper_chart", upper_->charts.front().handle}}},
        {"lower_matches", lower_->exact.matches.size()},
        {"upper_matches", upper_->exact.matches.size()},
        {"lower_tiles", lower_->exact.tiles.size()},
        {"upper_tiles", upper_->exact.tiles.size()},
        {"match_interval_queries", match_queries_.load()},
        {"tile_interval_queries", tile_queries_.load()},
        {"lower_match_advances", lower_match_advances_.load()},
        {"upper_match_advances", upper_match_advances_.load()},
        {"integrations", integrations_.load()},
        {"elapsed_ms", elapsed_ms_}};
    result["lower"] = encode_published_arm_summary(*lower_, "lower");
    result["upper"] = encode_published_arm_summary(*upper_, "upper");
    result["diagnostic_detail"] = "bounded-opaque-plan-publication";
    return result;
  }

  json::object checkpoint_record() const {
    if (!has_two_arms()) {
      const auto name = single_arm_name();
      return json::object{
          {"schema", kRetainedSingleArmTilePlanCheckpointSchema},
          {"handle", handle_},
          {"checkpoint_identity", checkpoint_identity_},
          {"provenance_identity", provenance_identity_},
          {"division_order", division_order_},
          {"arm_name", name},
          {"arm", encode_retained_arm(arm(name))},
          {"elapsed_ms", elapsed_ms_},
          {"runtime_stats", runtime_stats_record()}};
    }
    return json::object{
        {"schema", "diffexp3-retained-tile-plan-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"division_order", division_order_},
        {"lower", encode_retained_arm(*lower_)},
        {"upper", encode_retained_arm(*upper_)},
        {"elapsed_ms", elapsed_ms_},
        {"runtime_stats", runtime_stats_record()}};
  }

  void restore_runtime_stats(std::uint64_t match_queries,
                             std::uint64_t tile_queries,
                             std::uint64_t lower_match_advances,
                             std::uint64_t upper_match_advances,
                             std::uint64_t integrations) {
    if ((!lower_.has_value() && lower_match_advances != 0) ||
        (!upper_.has_value() && upper_match_advances != 0))
      throw std::invalid_argument(
          "single-arm tile-plan checkpoint advances an absent arm");
    match_queries_.store(match_queries);
    tile_queries_.store(tile_queries);
    lower_match_advances_.store(lower_match_advances);
    upper_match_advances_.store(upper_match_advances);
    integrations_.store(integrations);
  }

 private:
  void validate_arm_set() const {
    if (!lower_.has_value() && !upper_.has_value())
      throw std::invalid_argument(
          "retained tile plan must own one or two arms");
    if (lower_.has_value()) validate_exact_arm_plan(lower_->exact);
    if (upper_.has_value()) validate_exact_arm_plan(upper_->exact);
    // Existing two-arm requests name the slots but historically only require
    // opposite directions; preserve that behavior exactly. A genuine
    // single-arm plan derives its retained name from the exact direction.
    if (!has_two_arms()) {
      const auto& retained = lower_.has_value() ? *lower_ : *upper_;
      const auto expected_direction = lower_.has_value() ? -1 : 1;
      if (retained.exact.direction != expected_direction)
        throw std::invalid_argument(
            "retained single tile-arm name differs from its exact direction");
    }
  }

  std::string single_arm_name() const {
    if (has_two_arms())
      throw std::logic_error(
          "two-arm tile plan has no single arm name");
    return lower_.has_value() ? "lower" : "upper";
  }

  json::object runtime_stats_record() const {
    return json::object{
        {"match_interval_queries", match_queries_.load()},
        {"tile_interval_queries", tile_queries_.load()},
        {"lower_match_advances", lower_match_advances_.load()},
        {"upper_match_advances", upper_match_advances_.load()},
        {"integrations", integrations_.load()}};
  }

  json::object single_arm_summary(bool include_intervals) const {
    const auto name = single_arm_name();
    const auto& retained = arm(name);
    json::object result{
        {"tile_plan", handle_},
        {"capability", kRetainedSingleArmTilePlanCapability},
        {"native_retained", true},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"division_order", division_order_},
        {"independent_arms", false},
        {"concurrent_execution", "single-immutable-arm-snapshot"},
        {"arm_name", name},
        {"anchor", json::object{
             {"chart", retained.charts.front().handle},
             {"center_exact", retained.exact.from.str()}}},
        {"matches", retained.exact.matches.size()},
        {"tiles", retained.exact.tiles.size()},
        {"match_interval_queries", match_queries_.load()},
        {"tile_interval_queries", tile_queries_.load()},
        {"lower_match_advances", lower_match_advances_.load()},
        {"upper_match_advances", upper_match_advances_.load()},
        {"integrations", integrations_.load()},
        {"elapsed_ms", elapsed_ms_}};
    if (include_intervals) {
      result["arm"] = encode_retained_arm(retained);
      result[name] = encode_retained_arm(retained);
    }
    return result;
  }

  std::string handle_;
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::uint32_t division_order_ = 3;
  std::optional<RetainedArmPlan> lower_;
  std::optional<RetainedArmPlan> upper_;
  double elapsed_ms_ = 0.0;
  mutable std::atomic<std::uint64_t> match_queries_{0};
  mutable std::atomic<std::uint64_t> tile_queries_{0};
  std::atomic<std::uint64_t> lower_match_advances_{0};
  std::atomic<std::uint64_t> upper_match_advances_{0};
  std::atomic<std::uint64_t> integrations_{0};
};

// A plan-driven match is one exact handoff, not a completed arm.  It owns the
// immutable plan snapshot and every local used to construct the retained
// matching weights.  Registry release of the public plan/local tokens can
// therefore never turn a published handoff into dangling provenance.
class StoredPlannedMatchHop final : public StoredMatchBase {
 public:
  StoredPlannedMatchHop(
      std::shared_ptr<StoredMatchBase> match,
      std::string checkpoint_identity, std::string provenance_identity,
      json::object handoff, double elapsed_ms,
      std::shared_ptr<StoredTilePlan> plan_owner,
      std::vector<std::shared_ptr<StoredLocalBase>> basis_owners,
      std::shared_ptr<StoredLocalBase> incoming_owner)
      : StoredMatchBase(match == nullptr ? std::string() : match->handle()),
        match_(std::move(match)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        handoff_(std::move(handoff)), elapsed_ms_(elapsed_ms),
        plan_owner_(std::move(plan_owner)),
        basis_owners_(std::move(basis_owners)),
        incoming_owner_(std::move(incoming_owner)) {
    if (match_ == nullptr || plan_owner_ == nullptr ||
        incoming_owner_ == nullptr || basis_owners_.empty())
      throw std::invalid_argument(
          "retained planned match hop requires all strong owners");
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    if (acb &&
        acb->exact_shadow_factorized_basis().has_value()) {
      const auto& factorized =
          *acb->exact_shadow_factorized_basis();
      if (factorized.size() != basis_owners_.size())
        throw std::logic_error(
            "retained planned match exact-shadow basis disagrees with its "
            "receiving basis dimension");
      for (std::size_t column = 0;
           column < basis_owners_.size(); ++column) {
        const auto receiving =
            std::dynamic_pointer_cast<
                StoredLocal<ComplexBall>>(
                basis_owners_[column]);
        if (!receiving)
          throw std::logic_error(
              "retained planned match exact-shadow basis lost its Acb "
              "receiving local");
        local_algebra_detail::
            require_factorized_receiving_local_compatibility(
                factorized[column], receiving->solution(),
                "retained planned match exact-shadow basis column " +
                    std::to_string(column));
      }
    }
  }

  json::object summary() const override {
    auto result = match_->summary();
    if (required_string(result, "checkpoint_identity") !=
        checkpoint_identity_)
      throw std::logic_error(
          "retained planned match checkpoint identity changed");
    result["planned_hop_capability"] =
        kRetainedPlannedMatchHopCapability;
    result["plan_driven"] = true;
    result["planned_hop_provenance_identity"] = provenance_identity_;
    result["planned_hop"] = handoff_;
    result["strong_ownership"] = json::object{
        {"tile_plan", true}, {"basis_locals", basis_owners_.size()},
        {"incoming_local", true}};
    result["materializations"] = materializations_.load();
    result["elapsed_ms"] = elapsed_ms_;
    return result;
  }

  double elapsed_ms() const { return elapsed_ms_; }

  const std::shared_ptr<StoredMatchBase>& native_match() const {
    return match_;
  }
  bool has_acb_match() const {
    return std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_) !=
        nullptr;
  }
  json::object compact_terminal_match_diagnostic() const {
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    if (!acb)
      throw std::invalid_argument(
          "compact terminal match diagnostic requires one retained Acb match");
    return acb->compact_terminal_diagnostic_summary();
  }
  const std::string& terminal_acb_relative_tolerance_text() const {
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    if (!acb)
      throw std::invalid_argument(
          "terminal matching tolerance requires one retained Acb match");
    return acb->relative_tolerance_text();
  }
  const std::shared_ptr<StoredTilePlan>& plan_owner() const {
    return plan_owner_;
  }
  const std::vector<std::shared_ptr<StoredLocalBase>>& basis_owners() const {
    return basis_owners_;
  }
  const std::shared_ptr<StoredLocalBase>& incoming_owner() const {
    return incoming_owner_;
  }
  const json::object& handoff() const { return handoff_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }

  json::object terminal_mode_diagnostic_summary(
      int output_digits = 18) const {
    if (!has_terminal_acb_factorization())
      throw std::invalid_argument(
          "terminal mode diagnostics require one certified Acb exact-right match");
    const auto basis = terminal_acb_basis_owners();
    const auto& weights = terminal_acb_physical_weights();
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    if (!acb)
      throw std::logic_error(
          "terminal mode diagnostics lost their retained Acb match");
    const auto& transformed_weights = acb->transformed_weights();
    if (basis.size() != weights.size() ||
        basis.size() != transformed_weights.size())
      throw std::logic_error(
          "terminal mode diagnostics found different basis and weight dimensions");
    const auto compact_match = compact_terminal_match_diagnostic();
    json::array replay_basis;
    replay_basis.reserve(basis.size());
    for (const auto& owner : basis)
      replay_basis.push_back(json::object{
          {"local", owner->handle()},
          {"checkpoint_identity", owner->checkpoint_identity()}});
    const auto arm_name = required_string(handoff_, "arm");
    const auto match_index = static_cast<std::size_t>(
        as_u64(handoff_.at("match"),
               "terminal mode diagnostic replay match index"));
    const auto& replay_epsilon =
        as_object(compact_match.at("epsilon"),
                  "terminal mode diagnostic replay epsilon");
    const auto& replay_refinement =
        as_object(compact_match.at("refinement"),
                  "terminal mode diagnostic replay refinement");
    json::object replay_record{
        {"schema", "diffexp3-terminal-match-replay-v1"},
        {"plan", json::object{
             {"tile_plan", plan_owner_->handle()},
             {"checkpoint_identity",
              plan_owner_->checkpoint_identity()}}},
        {"arm", arm_name},
        {"index", match_index + 1},
        {"basis", std::move(replay_basis)},
        {"incoming", json::object{
             {"local", incoming_owner_->handle()},
             {"checkpoint_identity",
              incoming_owner_->checkpoint_identity()}}},
        {"epsilon", json::object{
             {"min", replay_epsilon.at("min")},
             {"max", replay_epsilon.at("max")},
             {"required_complete_max",
              replay_epsilon.at("required_complete_max")}}},
        {"refinement", json::object{
             {"relative_tolerance",
              replay_refinement.at("relative_tolerance")},
             {"max_steps", replay_refinement.at("max_steps")}}}};
    json::value owner_normal_frame_diagnostic = nullptr;
    try {
      const auto& retained = plan_owner_->arm(arm_name);
      if (match_index < retained.certified_matches.size()) {
        const auto owner = basis.front()->retained_equation_owner();
        if (owner)
          owner_normal_frame_diagnostic =
              owner->acb_matching_normal_frame_diagnostic(
                  retained.certified_matches[match_index]
                      .receiving.local);
      }
    } catch (const std::exception& error) {
      owner_normal_frame_diagnostic = json::object{
          {"schema", "diffexp3-acb-normal-frame-owner-diagnostic-v1"},
          {"status", "diagnostic-error"},
          {"detail", error.what()}};
    }
    const auto encode_weight =
        [output_digits](const EpsilonFrame<ComplexBall>& weight) {
          json::array coefficients;
          coefficients.reserve(weight.coefficients().size());
          for (std::int64_t raw_power = weight.min_power();
               raw_power <= weight.complete_max(); ++raw_power) {
            const auto power = static_cast<std::int32_t>(raw_power);
            coefficients.push_back(json::object{
                {"power", power},
                {"value", encode_scalar(
                    weight.coefficient(power), output_digits)},
                {"absolute_upper",
                 Magnitude::upper_abs(weight.coefficient(power))
                     .approximate_upper()}});
          }
          return json::object{
              {"min", weight.min_power()},
              {"max", weight.complete_max()},
              {"coefficients", std::move(coefficients)}};
        };

    json::array columns;
    columns.reserve(basis.size());
    for (std::size_t column = 0; column < basis.size(); ++column) {
      const auto& solution = basis[column]->solution();
      json::array sectors;
      sectors.reserve(solution.sectors.size());
      for (const auto& sector : solution.sectors)
        sectors.push_back(json::object{
            {"a", encode_exact_descriptor(sector.a)},
            {"b", encode_exact_descriptor(sector.b)},
            {"log_power", sector.log_power},
            {"material", local_detail::material_sector(sector)}});

      const auto& weight = weights[column];
      columns.push_back(json::object{
          {"column", column},
          {"local", basis[column]->handle()},
          {"epsilon",
           json::object{{"min", solution.epsilon.min_power},
                        {"max", solution.epsilon.complete_max}}},
          {"taylor_complete_max", solution.taylor_complete_max},
          {"sectors", std::move(sectors)},
          {"physical_weight", encode_weight(weight)},
          {"transformed_weight",
           encode_weight(transformed_weights[column])}});
    }

    // Matching residuals are measured at an interior handoff.  They do not
    // reveal whether the selected linear combination leaves a tiny
    // nonvanishing Frobenius mode at a singular endpoint.  Materialize the
    // physical combination only for this opt-in diagnostic and report every
    // coefficient whose absolute endpoint power is nonpositive.  Exact mag
    // dumps keep subnormal scales such as 10^-900 observable without
    // underflowing through a diagnostic double.
    std::vector<const LocalSolution<ComplexBall>*> physical_sources;
    physical_sources.reserve(basis.size());
    for (const auto& owner : basis)
      physical_sources.push_back(&owner->solution());
    const auto physical_solution = materialize_local_basis_weights(
        physical_sources, weights,
        checkpoint_identity_ + ":terminal-mode-diagnostic");
    json::array endpoint_modes;
    for (std::size_t sector_index = 0;
         sector_index < physical_solution.sectors.size();
         ++sector_index) {
      const auto& sector = physical_solution.sectors[sector_index];
      if (sector.a.domain != ExactDomain::Rational) continue;
      const Rational a(sector.a.canonical);
      for (std::uint32_t taylor = 0;
           taylor < physical_solution.taylor_width(); ++taylor) {
        const auto absolute_power = ExactScalarDescriptor::rational(
            (a + Rational(static_cast<long>(taylor))).str());
        if (absolute_power.sign == ExactSign::Positive) continue;
        json::array coefficient_scales;
        Magnitude overall = Magnitude::zero();
        std::optional<std::uint32_t> largest_component;
        std::optional<std::int32_t> largest_epsilon_power;
        for (std::size_t epsilon_index = 0;
             epsilon_index < physical_solution.epsilon.width();
             ++epsilon_index) {
          Magnitude epsilon_maximum = Magnitude::zero();
          std::optional<std::uint32_t> epsilon_component;
          for (std::uint32_t component = 0;
               component < physical_solution.dimension; ++component) {
            const auto magnitude = Magnitude::upper_abs(
                sector.coefficients[local_detail::sector_index(
                    physical_solution, epsilon_index, taylor,
                    component)]);
            if (!(magnitude <= epsilon_maximum)) {
              epsilon_maximum = magnitude;
              epsilon_component = component;
            }
            if (!(magnitude <= overall)) {
              overall = magnitude;
              largest_component = component;
              largest_epsilon_power =
                  local_detail::checked_i32(
                      static_cast<std::int64_t>(
                          physical_solution.epsilon.min_power) +
                          epsilon_index,
                      "terminal mode diagnostic epsilon power");
            }
          }
          if (!epsilon_maximum.is_zero())
            coefficient_scales.push_back(json::object{
                {"power",
                 local_detail::checked_i32(
                     static_cast<std::int64_t>(
                         physical_solution.epsilon.min_power) +
                         epsilon_index,
                     "terminal mode diagnostic coefficient power")},
                {"component", *epsilon_component},
                {"absolute_upper_exact",
                 epsilon_maximum.dump_exact()},
                {"absolute_upper_approx",
                 epsilon_maximum.approximate_upper()}});
        }
        if (overall.is_zero()) continue;
        const auto mode_class =
            sector.b.is_zero == TruthValue::No
            ? "regulated-nonpositive"
            : absolute_power.sign == ExactSign::Negative
            ? "unregulated-power-divergent"
            : sector.log_power == 0
            ? "unregulated-finite"
            : "unregulated-log-divergent";
        endpoint_modes.push_back(json::object{
            {"sector", sector_index},
            {"a", encode_exact_descriptor(sector.a)},
            {"b", encode_exact_descriptor(sector.b)},
            {"log_power", sector.log_power},
            {"taylor", taylor},
            {"absolute_power", absolute_power.canonical},
            {"class", mode_class},
            {"largest_component", *largest_component},
            {"largest_epsilon_power", *largest_epsilon_power},
            {"absolute_upper_exact", overall.dump_exact()},
            {"absolute_upper_approx", overall.approximate_upper()},
            {"coefficient_scales", std::move(coefficient_scales)}});
      }
    }
    return json::object{
        {"schema", "diffexp3-terminal-mode-diagnostic-v1"},
        {"match", compact_match},
        {"replay", std::move(replay_record)},
        {"owner_normal_frame_diagnostic",
         std::move(owner_normal_frame_diagnostic)},
        {"columns", std::move(columns)},
        {"physical_endpoint_modes", std::move(endpoint_modes)}};
  }
  const std::string& provenance_identity() const {
    return provenance_identity_;
  }

  bool has_terminal_acb_factorization() const {
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    return acb != nullptr &&
        acb->certified_for_materialization() &&
        (acb->acb_materialization_right_transformation().has_value() ||
         acb->terminal_normal_frame_right_transformation().has_value());
  }

  std::vector<std::shared_ptr<StoredLocal<ComplexBall>>>
  terminal_acb_basis_owners() const {
    if (!has_terminal_acb_factorization())
      throw std::invalid_argument(
          "terminal factorized consumption requires one certified Acb exact-right match");
    std::vector<std::shared_ptr<StoredLocal<ComplexBall>>> result;
    result.reserve(basis_owners_.size());
    for (const auto& owner : basis_owners_) {
      auto typed =
          std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(owner);
      if (!typed)
        throw std::logic_error(
            "terminal factorized Acb basis changed coefficient domain");
      result.push_back(std::move(typed));
    }
    return result;
  }

  std::shared_ptr<const PreparedPhysicalClearedODE<ComplexBall>>
  terminal_acb_receiving_physical_equation() const {
    auto owners = terminal_acb_basis_owners();
    if (owners.empty())
      throw std::logic_error(
          "terminal composed observable has no receiving physical basis owner");
    const auto equation = owners.front()->physical_equation();
    if (!equation)
      throw std::domain_error(
          "terminal composed observable has no retained receiving physical equation");
    for (std::size_t column = 1; column < owners.size(); ++column) {
      const auto candidate = owners[column]->physical_equation();
      if (!candidate ||
          candidate->payload_identity != equation->payload_identity ||
          candidate->owner_signature_identity !=
              equation->owner_signature_identity)
        throw std::logic_error(
            "terminal composed observable basis columns disagree on their physical equation");
    }
    return equation;
  }

  std::shared_ptr<const PreparedPhysicalClearedODE<Rational>>
  terminal_rational_shadow_physical_equation() const {
    const auto owners = terminal_acb_basis_owners();
    if (owners.empty())
      throw std::logic_error(
          "terminal composed observable has no Rational-shadow basis owner");
    std::vector<std::shared_ptr<const RationalShadowColumnWitness>> witnesses;
    witnesses.reserve(owners.size());
    bool any_witness = false;
    bool all_witnesses = true;
    for (const auto& owner : owners) {
      auto witness = owner->rational_shadow_witness();
      any_witness = any_witness || static_cast<bool>(witness);
      all_witnesses = all_witnesses && static_cast<bool>(witness);
      witnesses.push_back(std::move(witness));
    }
    if (!any_witness) {
      const auto acb =
          std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
      auto exact_owner = acb ? acb->exact_shadow_equation_owner() : nullptr;
      if (!exact_owner) {
        const auto receiving_owner =
            std::dynamic_pointer_cast<CompositeSCCChartBase>(
                owners.front()->retained_equation_owner());
        if (receiving_owner) {
          auto equation =
              receiving_owner->rational_shadow_physical_equation();
          if (equation) {
            if (equation->dimension != owners.size())
              throw std::logic_error(
                  "terminal composed observable manifest Rational-shadow equation has the wrong dimension");
            return equation;
          }
        }
      }
      if (!exact_owner) return nullptr;
      if (std::string(exact_owner->equation_scalar_domain()) != "rational")
        throw std::logic_error(
            "terminal composed observable retained shadow owner is not Rational");
      auto erased = exact_owner->physical_ode_erased();
      if (!erased)
        throw std::domain_error(
            "terminal composed observable retained shadow owner has no physical equation");
      auto equation = std::static_pointer_cast<
          const PreparedPhysicalClearedODE<Rational>>(std::move(erased));
      if (!equation || equation->dimension != owners.size())
        throw std::logic_error(
            "terminal composed observable retained shadow equation has the wrong dimension");
      return equation;
    }
    if (!all_witnesses)
      throw std::logic_error(
          "terminal composed observable has only a partial Rational-shadow basis");
    std::shared_ptr<const PhysicalEquationOwnerBase> exact_owner;
    std::string shadow_identity;
    for (std::size_t column = 0; column < owners.size(); ++column) {
      const auto& witness = witnesses[column];
      if (!witness || !witness->exact_equation_owner ||
          witness->rational_shadow_identity.empty())
        throw std::domain_error(
            "terminal composed observable has no retained exact Rational-shadow equation");
      if (column == 0) {
        exact_owner = witness->exact_equation_owner;
        shadow_identity = witness->rational_shadow_identity;
      } else if (witness->exact_equation_owner.get() != exact_owner.get() ||
                 witness->rational_shadow_identity != shadow_identity) {
        throw std::logic_error(
            "terminal composed observable columns disagree on their exact Rational-shadow owner");
      }
    }
    if (!exact_owner ||
        std::string(exact_owner->equation_scalar_domain()) != "rational")
      throw std::logic_error(
          "terminal composed observable exact shadow is not Rational");
    const auto* retained_shadow =
        exact_owner->matching_scc_rational_shadow_identity();
    if (retained_shadow == nullptr || *retained_shadow != shadow_identity)
      throw std::logic_error(
          "terminal composed observable exact owner has a different Rational-shadow identity");
    auto erased = exact_owner->physical_ode_erased();
    if (!erased)
      throw std::domain_error(
          "terminal composed observable exact owner has no physical equation");
    auto equation = std::static_pointer_cast<
        const PreparedPhysicalClearedODE<Rational>>(std::move(erased));
    if (!equation || equation->dimension != owners.size())
      throw std::logic_error(
          "terminal composed observable exact physical equation has the wrong dimension");
    return equation;
  }

  std::shared_ptr<const
      adjoint_observable_detail::NormalizedBackwardAdjointExactODE>
  terminal_normalized_backward_adjoint_exact_equation(
      const std::shared_ptr<const PreparedPhysicalClearedODE<Rational>>&
          exact_equation,
      std::uint32_t taylor_complete_max,
      std::int32_t epsilon_complete_max,
      const std::string& context) const {
    if (!exact_equation)
      throw std::invalid_argument(
          context +
          ": normalized terminal adjoint cache needs an exact equation");
    const auto key =
        std::make_pair(taylor_complete_max, epsilon_complete_max);
    const auto contract = exact_equation->payload_identity + ":" +
        exact_equation->owner_signature_identity +
        ":terminal-normalized-adjoint-v1";
    return terminal_normalized_adjoint_cache_.get_or_build(
        key, contract, [&] {
          return adjoint_observable_detail::
              normalize_backward_adjoint_exact_ode_by_q(
                  *exact_equation, taylor_complete_max,
                  epsilon_complete_max, context);
        }).value;
  }

  adjoint_observable_detail::BackwardAdjointRealRayOperatorCache&
  terminal_backward_adjoint_real_ray_operator_cache() const {
    return terminal_real_ray_operator_cache_;
  }

  std::uint32_t terminal_composed_taylor_order_floor() const {
    return terminal_composed_taylor_order_floor_.load(
        std::memory_order_relaxed);
  }

  void learn_terminal_composed_taylor_order_floor(
      std::uint32_t order) const {
    auto current = terminal_composed_taylor_order_floor_.load(
        std::memory_order_relaxed);
    while (current < order &&
           !terminal_composed_taylor_order_floor_.compare_exchange_weak(
               current, order, std::memory_order_relaxed,
               std::memory_order_relaxed)) {}
  }

  FiniteLaurentVector<ComplexBall>
  terminal_acb_incoming_physical_value() const {
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    auto incoming =
        std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(incoming_owner_);
    if (!acb || !incoming || !has_terminal_acb_factorization())
      throw std::invalid_argument(
          "terminal composed observable requires one retained Acb match and incoming local");
    EvaluationOptions options;
    options.imaginary_sign = acb->effective_incoming_imaginary_sign();
    options.compute_tail_estimate = false;
    const auto evaluation = evaluate_local_solution(
        incoming->solution(),
        RealEvaluationPoint::rational(acb->incoming_point_exact()), options);
    if (!evaluation.value.error.empty() ||
        evaluation.value.dimension != incoming->solution().dimension)
      throw std::domain_error(
          "terminal composed observable incoming evaluation has an unsupported error envelope or dimension");
    FiniteLaurentVector<ComplexBall> value;
    value.reserve(evaluation.value.dimension);
    for (std::uint32_t component = 0;
         component < evaluation.value.dimension; ++component) {
      std::vector<ComplexBall> coefficients;
      coefficients.reserve(evaluation.value.epsilon.width());
      for (std::int64_t raw_power = evaluation.value.epsilon.min_power;
           raw_power <= evaluation.value.epsilon.complete_max; ++raw_power)
        coefficients.push_back(evaluation.value.at(
            static_cast<std::int32_t>(raw_power), component));
      value.emplace_back(evaluation.value.epsilon, std::move(coefficients));
    }
    return value;
  }

  const FiniteLaurentVector<ComplexBall>&
  terminal_acb_physical_weights() const {
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    if (!acb || !has_terminal_acb_factorization())
      throw std::invalid_argument(
          "terminal physical weights require one certified Acb exact-right match");
    return acb->weights();
  }

  const FiniteLaurentVector<ComplexBall>&
  terminal_acb_transformed_weights() const {
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    if (!acb || !has_terminal_acb_factorization())
      throw std::invalid_argument(
          "terminal transformed weights require one certified Acb exact-right match");
    return acb->transformed_weights();
  }

  std::int32_t terminal_acb_physical_weight_min_power() const {
    const auto& weights = terminal_acb_physical_weights();
    if (weights.empty())
      throw std::logic_error(
          "terminal physical contraction has no match weights");
    auto result = weights.front().min_power();
    for (const auto& weight : weights)
      result = std::min(result, weight.min_power());
    return result;
  }

  std::size_t terminal_acb_extra_precision_bits() const {
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    if (!acb || !has_terminal_acb_factorization())
      throw std::invalid_argument(
          "terminal extra precision requires one certified Acb factorization");
    return acb->exact_shadow_extra_precision_bits();
  }

  std::shared_ptr<const std::vector<LocalSolution<ComplexBall>>>
  terminal_acb_factorized_basis(
      const std::string& checkpoint_identity_root) const {
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    if (!acb || !has_terminal_acb_factorization() ||
        checkpoint_identity_root.empty())
      throw std::invalid_argument(
          "terminal factorized basis requires one certified Acb exact-right match and checkpoint root");
    const auto precision_bits = ComplexBall::precision();
    std::lock_guard<std::mutex> cache_lock(
        terminal_factorized_basis_cache_mutex_);
    if (terminal_factorized_basis_cache_ &&
        terminal_factorized_basis_cache_precision_bits_ == precision_bits)
      return terminal_factorized_basis_cache_;
    auto owners = terminal_acb_basis_owners();
    std::vector<LocalSolution<ComplexBall>> valid_physical_basis;
    valid_physical_basis.reserve(owners.size());
    std::vector<const LocalSolution<ComplexBall>*> physical_sources;
    physical_sources.reserve(owners.size());
    for (std::size_t column = 0; column < owners.size(); ++column) {
      auto source = owners[column]->solution();
      const auto complete_max = std::min(
          owners[column]->top_valid(),
          source.epsilon.complete_max);
      if (complete_max < source.epsilon.min_power)
        throw std::domain_error(
            "terminal physical basis has no valid epsilon coefficient");
      if (complete_max < source.epsilon.complete_max)
        source = restrict_local_epsilon_frame_strict_lower(
            source, source.epsilon.min_power, complete_max,
            checkpoint_identity_root + ":physical-valid:column:" +
                std::to_string(column));
      valid_physical_basis.push_back(std::move(source));
      physical_sources.push_back(&valid_physical_basis.back());
    }

    std::vector<LocalSolution<ComplexBall>> factorized;
    if (const auto& exact_shadow =
            acb->exact_shadow_factorized_basis();
        exact_shadow.has_value()) {
      factorized = *exact_shadow;
    } else if (const auto& normal_right =
            acb->terminal_normal_frame_right_transformation();
        normal_right.has_value()) {
      if (normal_right->size() != physical_sources.size() ||
          std::any_of(normal_right->begin(), normal_right->end(),
              [&](const auto& row) {
                return row.size() != physical_sources.size();
              }))
        throw std::logic_error(
            "terminal normal-frame right transformation is not square");
      std::vector<LocalSolution<ComplexBall>> normal_basis;
      normal_basis.reserve(physical_sources.size());
      for (std::size_t output = 0;
           output < physical_sources.size(); ++output) {
        FiniteLaurentVector<ComplexBall> weights;
        weights.reserve(physical_sources.size());
        for (std::size_t input = 0;
             input < physical_sources.size(); ++input)
          weights.push_back((*normal_right)[input][output]);
        normal_basis.push_back(materialize_local_basis_weights(
            physical_sources, weights,
            checkpoint_identity_root + ":finite-normal-right:column:" +
                std::to_string(output)));
      }
      std::vector<const LocalSolution<ComplexBall>*> normal_sources;
      normal_sources.reserve(normal_basis.size());
      for (const auto& column : normal_basis)
        normal_sources.push_back(&column);
      if (const auto& acb_right =
              acb->acb_materialization_right_transformation();
          acb_right.has_value())
        factorized = right_transform_local_basis_exact<ComplexBall>(
            normal_sources, *acb_right,
            checkpoint_identity_root +
                ":normal-frame-acb-Levelt-right");
      else {
        const auto exact_right =
            specialize_exact_rational_laurent_matrix_to_acb(
                acb->exact_saturation_transformation());
        factorized = right_transform_local_basis_exact<ComplexBall>(
            normal_sources, exact_right,
            checkpoint_identity_root +
                ":normal-frame-exact-right");
      }
    } else {
      factorized =
          right_transform_local_basis_exact<ComplexBall>(
              physical_sources,
              *acb->acb_materialization_right_transformation(),
              checkpoint_identity_root + ":exact-right");
    }
    if (const auto& preconditioner =
            acb->acb_right_materialization_preconditioner();
        preconditioner.has_value()) {
      std::vector<const LocalSolution<ComplexBall>*> factorized_sources;
      factorized_sources.reserve(factorized.size());
      for (const auto& column : factorized)
        factorized_sources.push_back(&column);
      factorized =
          right_transform_local_basis_exact<ComplexBall>(
              factorized_sources, *preconditioner,
              checkpoint_identity_root + ":conditioned-exact-right");
    }
    terminal_factorized_basis_cache_ =
        std::make_shared<const std::vector<LocalSolution<ComplexBall>>>(
            std::move(factorized));
    terminal_factorized_basis_cache_precision_bits_ = precision_bits;
    return terminal_factorized_basis_cache_;
  }

  FiniteLaurentVector<ComplexBall>
  contract_terminal_acb_physical_functionals(
      const FiniteLaurentMatrix<ComplexBall>& physical_rows,
      const std::string& context) const {
    if (physical_rows.empty())
      throw std::invalid_argument(
          context + ": terminal physical functional batch is empty");
    const auto& weights = terminal_acb_physical_weights();
    for (const auto& row : physical_rows)
      if (row.size() != weights.size())
        throw std::invalid_argument(
            context +
            ": terminal physical functional differs from the match dimension");
    return apply_finite_laurent_matrix(physical_rows, weights);
  }

  FiniteLaurentVector<ComplexBall>
  contract_terminal_acb_factorized_functionals(
      const FiniteLaurentMatrix<ComplexBall>& factorized_rows,
      const std::string& context) const {
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    if (!acb || !has_terminal_acb_factorization())
      throw std::invalid_argument(
          context +
          ": terminal factorized contraction requires one certified Acb exact-right match");
    if (factorized_rows.empty())
      throw std::invalid_argument(
          context + ": terminal factorized functional batch is empty");
    for (const auto& row : factorized_rows)
      if (row.size() != basis_owners_.size())
        throw std::invalid_argument(
            context +
            ": terminal factorized functional differs from the transformed basis dimension");
    return apply_finite_laurent_matrix(
        factorized_rows, acb->transformed_weights());
  }

  FiniteLaurentVector<ComplexBall>
  adjoint_contract_terminal_acb_functionals(
      const FiniteLaurentMatrix<ComplexBall>& functional_rows,
      std::int32_t required_output_complete_max,
      const std::string& context,
      bool factorized_coordinates = false,
      std::optional<Magnitude> publication_relative_tolerance =
          std::nullopt) const {
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    if (!acb || !has_terminal_acb_factorization())
      throw std::invalid_argument(
          context +
          ": terminal adjoint contraction requires one certified Acb exact-right match");
    if (!factorized_coordinates &&
        acb->terminal_normal_frame_right_transformation().has_value())
      throw std::invalid_argument(
          context +
          ": adjoint diagnostic is not implemented for a finite normal-frame terminal transformation");
    if (functional_rows.empty())
      throw std::invalid_argument(
          context + ": terminal adjoint functional batch is empty");
    auto physical_basis_owners = factorized_coordinates
        ? std::vector<std::shared_ptr<StoredLocal<ComplexBall>>>{}
        : terminal_acb_basis_owners();
    auto factorized_basis = factorized_coordinates
        ? terminal_acb_factorized_basis(checkpoint_identity_)
        : std::shared_ptr<const std::vector<LocalSolution<ComplexBall>>>{};
    std::vector<const LocalSolution<ComplexBall>*> basis;
    if (factorized_coordinates) {
      basis.reserve(factorized_basis->size());
      for (const auto& column : *factorized_basis)
        basis.push_back(&column);
    } else {
      basis.reserve(physical_basis_owners.size());
      for (const auto& owner : physical_basis_owners)
        basis.push_back(&owner->solution());
    }
    const auto dimension = basis.size();
    for (const auto& row : functional_rows)
      if (row.size() != dimension)
        throw std::invalid_argument(
          context +
              ": terminal adjoint functional differs from its selected basis dimension");
    auto incoming =
        std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(
            incoming_owner_);
    if (!incoming ||
        incoming->solution().dimension != dimension)
      throw std::logic_error(
          context +
          ": terminal adjoint contraction lost its incoming physical local");

    const auto basis_point =
        RealEvaluationPoint::rational(acb->basis_point_exact());
    const auto incoming_point =
        RealEvaluationPoint::rational(acb->incoming_point_exact());
    EvaluationOptions basis_options;
    basis_options.imaginary_sign =
        acb->effective_basis_imaginary_sign();
    basis_options.compute_tail_estimate = false;
    EvaluationOptions incoming_options;
    incoming_options.imaginary_sign =
        acb->effective_incoming_imaginary_sign();
    incoming_options.compute_tail_estimate = false;
    std::vector<LocalEvaluation> basis_evaluations;
    basis_evaluations.reserve(dimension);
    auto evaluation_window = acb->requested_window();
    // The published match window is only the residual-certification
    // contract.  A terminal adjoint can start several epsilon powers below
    // that window, and its coefficient recurrence then needs the positive
    // tail of the retained local evaluations.  Re-evaluate the same certified
    // Taylor prefix at its full epsilon width instead of silently clipping it
    // back to the public match maximum.
    evaluation_window.complete_max =
        std::numeric_limits<std::int32_t>::max();
    for (std::size_t column = 0; column < dimension; ++column) {
      auto options = basis_options;
      // A shorter Taylor prefix may be selected to keep the pointwise
      // connection-weight solve well conditioned.  It is not a new local
      // solution and cannot truncate the terminal observable.  The adjoint
      // evaluates the complete retained basis and incoming local here, just
      // as project_terminal_acb_basis_row does before line integration.
      basis_evaluations.push_back(evaluate_local_solution(
          *basis[column], basis_point, options));
      evaluation_window.min_power = std::min(
          evaluation_window.min_power,
          basis_evaluations.back().value.epsilon.min_power);
      evaluation_window.complete_max = std::min(
          evaluation_window.complete_max,
          basis_evaluations.back().value.epsilon.complete_max);
    }
    auto incoming_evaluation = evaluate_local_solution(
        incoming->solution(), incoming_point,
        incoming_options);
    diagnose_plain_regular_taylor_evaluation(
        incoming->solution(), incoming_point,
        incoming_evaluation, context);
    const auto recentered_diagnostic =
        diagnose_two_step_terminal_physical_recenter(
            incoming->solution(), incoming->physical_equation(),
            incoming_point, incoming_options, context);
    if (recentered_diagnostic.has_value()) {
      const auto* requested =
          std::getenv("DIFFEXP_DIAGNOSTIC_TERMINAL_INCOMING_TAYLOR");
      if (requested != nullptr) {
        const std::string request(requested);
        const auto separator = request.find(':');
        if (separator == std::string::npos)
          throw std::invalid_argument(
              "DIFFEXP_DIAGNOSTIC_TERMINAL_INCOMING_TAYLOR must be "
              "<epsilon-power>:<component>");
        const auto epsilon_power = static_cast<std::int32_t>(
            std::stoll(request.substr(0, separator)));
        const auto component = static_cast<std::uint32_t>(
            std::stoull(request.substr(separator + 1)));
        if (epsilon_power >=
                recentered_diagnostic->value.epsilon.min_power &&
            epsilon_power <=
                recentered_diagnostic->value.epsilon.complete_max &&
            component < recentered_diagnostic->value.dimension) {
          const auto& direct =
              incoming_evaluation.value.at(epsilon_power, component);
          const auto& recentered =
              recentered_diagnostic->value.at(epsilon_power, component);
          const auto discrepancy = recentered - direct;
          std::cerr
              << "terminal-physical-recenter-comparison context="
              << context
              << " epsilon_power=" << epsilon_power
              << " component=" << component
              << " direct_midpoint=("
              << direct.real_midpoint(12) << ","
              << direct.imag_midpoint(12) << ")"
              << " direct_radius2exp=("
              << direct.real_radius_exponent() << ","
              << direct.imag_radius_exponent() << ")"
              << " recentered_midpoint=("
              << recentered.real_midpoint(12) << ","
              << recentered.imag_midpoint(12) << ")"
              << " recentered_radius2exp=("
              << recentered.real_radius_exponent() << ","
              << recentered.imag_radius_exponent() << ")"
              << " overlaps=" << discrepancy.contains_zero()
              << " discrepancy_radius2exp=("
              << discrepancy.real_radius_exponent() << ","
              << discrepancy.imag_radius_exponent() << ")"
              << '\n';
        }
      }
    }
    const auto factorized_evaluation_diagnostic =
        diagnose_factorized_terminal_physical_evaluation(
            incoming->solution(), incoming->physical_equation(),
            incoming_point, incoming_options, context);
    if (factorized_evaluation_diagnostic.has_value()) {
      const auto* requested =
          std::getenv("DIFFEXP_DIAGNOSTIC_TERMINAL_INCOMING_TAYLOR");
      if (requested != nullptr) {
        const std::string request(requested);
        const auto separator = request.find(':');
        if (separator == std::string::npos)
          throw std::invalid_argument(
              "DIFFEXP_DIAGNOSTIC_TERMINAL_INCOMING_TAYLOR must be "
              "<epsilon-power>:<component>");
        const auto epsilon_power = static_cast<std::int32_t>(
            std::stoll(request.substr(0, separator)));
        const auto component = static_cast<std::uint32_t>(
            std::stoull(request.substr(separator + 1)));
        if (epsilon_power >= factorized_evaluation_diagnostic
                                 ->value.epsilon.min_power &&
            epsilon_power <= factorized_evaluation_diagnostic
                                 ->value.epsilon.complete_max &&
            component <
                factorized_evaluation_diagnostic->value.dimension) {
          const auto& direct =
              incoming_evaluation.value.at(epsilon_power, component);
          const auto& factorized =
              factorized_evaluation_diagnostic->value.at(
                  epsilon_power, component);
          const auto discrepancy = factorized - direct;
          std::cerr
              << "terminal-factorized-physical-comparison context="
              << context
              << " epsilon_power=" << epsilon_power
              << " component=" << component
              << " direct_midpoint=("
              << direct.real_midpoint(12) << ","
              << direct.imag_midpoint(12) << ")"
              << " direct_radius2exp=("
              << direct.real_radius_exponent() << ","
              << direct.imag_radius_exponent() << ")"
              << " factorized_midpoint=("
              << factorized.real_midpoint(12) << ","
              << factorized.imag_midpoint(12) << ")"
              << " factorized_radius2exp=("
              << factorized.real_radius_exponent() << ","
              << factorized.imag_radius_exponent() << ")"
              << " overlaps=" << discrepancy.contains_zero()
              << " discrepancy_midpoint=("
              << discrepancy.real_midpoint(12) << ","
              << discrepancy.imag_midpoint(12) << ")"
              << " discrepancy_radius2exp=("
              << discrepancy.real_radius_exponent() << ","
              << discrepancy.imag_radius_exponent() << ")"
              << '\n';
        }
      }
    }
    auto factorized_incoming_transfer =
        certified_ordinary_center_factorization(
            incoming->solution(), incoming->physical_equation(),
            incoming_point, incoming_options,
            incoming_evaluation, context);
    if (factorized_incoming_transfer.has_value())
      incoming_evaluation =
          factorized_incoming_transfer->evaluation;
    std::optional<ErrorEnvelope> certified_incoming_tail;
    const auto* terminal_tail_mode =
        std::getenv("DIFFEXP_DIAGNOSTIC_TERMINAL_PHYSICAL_TAIL");
    if (terminal_tail_mode != nullptr) {
      const std::string mode(terminal_tail_mode);
      if (mode != "report" && mode != "inflate")
        throw std::invalid_argument(
            "DIFFEXP_DIAGNOSTIC_TERMINAL_PHYSICAL_TAIL must be report or inflate");
      const auto& equation = incoming->physical_equation();
      if (!equation) {
        std::cerr
            << "terminal-physical-tail context=" << context
            << " status=unsupported detail=no-retained-physical-equation\n";
      } else {
        const auto prepared =
            prepare_physical_regular_homogeneous_tail_model(
                *equation, incoming->solution());
        if (prepared.status != TailMajorantStatus::Certified ||
            !prepared.model.has_value()) {
          std::cerr
              << "terminal-physical-tail context=" << context
              << " status=" << tail_majorant_status_name(prepared.status)
              << " detail=" << prepared.detail << '\n';
        } else {
          const auto arm_name = required_string(handoff_, "arm");
          const auto match_index = static_cast<std::size_t>(
              as_u64(handoff_.at("match"),
                     "terminal physical-tail match index"));
          const auto& arm = plan_owner_->arm(arm_name);
          if (match_index >= arm.exact.matches.size())
            throw std::logic_error(
                context +
                ": terminal physical-tail match is outside its retained arm");
          const auto& exact_match = arm.exact.matches[match_index];
          if (exact_match.producing_local.str() !=
              acb->incoming_point_exact())
            throw std::logic_error(
                context +
                ": terminal physical-tail point differs from retained match geometry");
          const auto& producing =
              arm.charts.at(exact_match.producing_chart);
          const auto point_modulus =
              exact_path_detail::abs(exact_match.producing_local);
          const auto witness_gap =
              producing.geometry.radius - point_modulus;
          constexpr std::uint32_t kWitnessSearchCap = 16;
          std::optional<Rational> witness_radius;
          Rational dyadic_denominator(1);
          for (std::uint32_t exponent = 1;
               exponent <= kWitnessSearchCap; ++exponent) {
            dyadic_denominator *= Rational(2);
            const auto candidate =
                point_modulus + witness_gap / dyadic_denominator;
            const auto candidate_certificate =
                certify_physical_regular_taylor_point_tail(
                    *prepared.model, incoming_point, candidate.str());
            if (candidate_certificate.status ==
                TailMajorantStatus::Certified) {
              witness_radius = candidate;
              break;
            }
          }
          if (!witness_radius.has_value()) {
            std::cerr
                << "terminal-physical-tail context=" << context
                << " status=inconclusive detail=no-certified-dyadic-witness\n";
          } else {
            auto certified =
                evaluate_physical_local_solution_with_certified_tail(
                    *prepared.model, incoming_point,
                    witness_radius->str(), incoming_options);
            if (certified.tail.status !=
                TailMajorantStatus::Certified) {
              std::cerr
                  << "terminal-physical-tail context=" << context
                  << " status="
                  << tail_majorant_status_name(certified.tail.status)
                  << " detail=" << certified.tail.detail << '\n';
            } else {
              if (mode == "inflate")
                incoming_evaluation = std::move(certified.evaluation);
              certified_incoming_tail = certified.tail.value;
              double maximum_tail_upper = 0.0;
              for (const auto& bound :
                   certified_incoming_tail->absolute)
                maximum_tail_upper = std::max(
                    maximum_tail_upper, bound.approximate_upper());
              std::cerr
                  << "terminal-physical-tail context=" << context
                  << " status=certified mode=" << mode
                  << " witness_radius_exact="
                  << witness_radius->str()
                  << " retained_taylor_order="
                  << prepared.model->taylor_complete_max
                  << " tail_upper_max=" << maximum_tail_upper
                  << '\n';
            }
          }
        }
      }
    }
    evaluation_window.min_power = std::min(
        evaluation_window.min_power,
        incoming_evaluation.value.epsilon.min_power);
    evaluation_window.complete_max = std::min(
        evaluation_window.complete_max,
        incoming_evaluation.value.epsilon.complete_max);
    if (evaluation_window.complete_max <
        acb->certified_complete_max())
      throw MatchingArithmeticError(
          MatchingArithmeticErrorCode::InsufficientCompleteWindow,
          context +
              ": terminal adjoint evaluation does not cover the certified match prefix",
          std::nullopt, std::nullopt,
          evaluation_window.complete_max);

    FiniteLaurentMatrix<ComplexBall> evaluated_basis(
        dimension, FiniteLaurentVector<ComplexBall>());
    for (auto& row : evaluated_basis) row.reserve(dimension);
    for (std::size_t column = 0; column < dimension; ++column) {
      auto frames = acb_evaluation_frames(
          basis_evaluations[column].value,
          evaluation_window,
          context + ": basis evaluation");
      for (std::size_t component = 0;
           component < dimension; ++component)
        evaluated_basis[component].push_back(
            std::move(frames[component]));
    }
    auto incoming_value = acb_evaluation_frames(
        incoming_evaluation.value, evaluation_window,
        context + ": incoming evaluation");
    if (certified_incoming_tail.has_value() &&
        std::string(terminal_tail_mode) == "inflate") {
      if (!tail_majorant_detail::same_epsilon_window(
              certified_incoming_tail->frame,
              incoming_evaluation.value.epsilon) ||
          certified_incoming_tail->absolute.size() !=
              incoming_evaluation.value.epsilon.width())
        throw std::logic_error(
            context +
            ": certified terminal incoming tail changed its epsilon frame");
      for (auto& frame : incoming_value) {
        auto coefficients = frame.coefficients();
        for (std::int64_t raw_power = frame.min_power();
             raw_power <= frame.complete_max(); ++raw_power) {
          if (raw_power < certified_incoming_tail->frame.min_power ||
              raw_power > certified_incoming_tail->frame.complete_max)
            continue;
          const auto power = local_algebra_detail::checked_i32(
              raw_power, "terminal certified-tail epsilon power");
          const auto row = static_cast<std::size_t>(
              raw_power - certified_incoming_tail->frame.min_power);
          certified_incoming_tail->absolute[row].add_error_to(
              coefficients.at(static_cast<std::size_t>(
                  raw_power - frame.min_power())));
        }
        frame = EpsilonFrame<ComplexBall>(
            frame.window(), std::move(coefficients));
      }
    }
    // `evaluation_window` is shared with the transformed basis so that A's
    // rows align, but the physical incoming vector u need not begin at A's
    // lower edge.  Remove only exact zero padding.  An ambiguous Acb ball is
    // retained as possible support; it is never guessed to be zero.
    for (auto& frame : incoming_value)
      frame = matching_detail::
          trim_leading_exact_zeros_preserve_ambiguous(frame);

    auto conditioned_basis = evaluated_basis;
    auto conditioned_rows = functional_rows;
    if (!factorized_coordinates) {
      conditioned_basis =
          right_multiply_finite_by_exact_laurent(
              evaluated_basis,
              *acb->acb_materialization_right_transformation());
      conditioned_rows =
          right_multiply_finite_by_exact_laurent(
              functional_rows,
              *acb->acb_materialization_right_transformation());
      if (const auto& preconditioner =
              acb->acb_right_materialization_preconditioner();
          preconditioner.has_value()) {
        conditioned_basis =
            right_multiply_finite_by_exact_laurent(
                conditioned_basis, *preconditioner);
        conditioned_rows =
            right_multiply_finite_by_exact_laurent(
                conditioned_rows, *preconditioner);
      }
    }

    auto conditioned_basis_complete_max =
        conditioned_basis.front().front().complete_max();
    for (const auto& row : conditioned_basis)
      for (const auto& frame : row)
        conditioned_basis_complete_max = std::min(
            conditioned_basis_complete_max, frame.complete_max());
    auto conditioned_functional_min =
        conditioned_rows.front().front().min_power();
    for (const auto& row : conditioned_rows)
      for (const auto& frame : row)
        conditioned_functional_min = std::min(
            conditioned_functional_min, frame.min_power());
    auto incoming_min = incoming_value.front().min_power();
    auto incoming_complete_max =
        incoming_value.front().complete_max();
    for (const auto& frame : incoming_value) {
      incoming_min = std::min(incoming_min, frame.min_power());
      incoming_complete_max = std::min(
          incoming_complete_max, frame.complete_max());
    }
    // The legacy physical-coordinate diagnostic reconstructs the exact-right
    // frame numerically and may leave tiny enclosures around zeros proved by
    // the retained lattice.  Its historical structural chop stays confined to
    // that proposal.  The production factorized adjoint consumes the retained
    // G=(F*T)*P columns directly and must preserve every ball: a locally tiny
    // coefficient can be endpoint-significant.
    if (!factorized_coordinates) {
      constexpr int structural_chop_digits = 100;
      for (auto& row : conditioned_basis)
        for (auto& frame : row) {
          auto coefficients = frame.coefficients();
          for (auto& coefficient : coefficients)
            coefficient =
                ScalarTraits<ComplexBall>::canonicalized(
                    coefficient, structural_chop_digits);
          frame = EpsilonFrame<ComplexBall>(
              frame.window(), std::move(coefficients));
        }
    }

    FiniteLaurentMatrix<ComplexBall> transposed_basis(
        dimension, FiniteLaurentVector<ComplexBall>());
    for (auto& row : transposed_basis)
      row.reserve(dimension);
    for (std::size_t row = 0; row < dimension; ++row)
      for (std::size_t column = 0;
           column < dimension; ++column)
        transposed_basis[column].push_back(
            conditioned_basis[row][column]);

    // The exact saturation witness proves that the transformed basis has no
    // negative epsilon powers.  Numerical evaluation can retain wide balls
    // around those algebraic zeros, especially in the final CASE-P row.  The
    // coefficient recurrence is therefore factored on this projected
    // matrix, and its forward certificate must use that same exact-proven
    // matrix.  Certifying against `transposed_basis` instead would reject the
    // solve for numerical remnants which the lattice proof has already shown
    // to be identically zero.
    auto certified_nonnegative_transposed_basis =
        exact_nonnegative_finite_laurent_matrix(
            transposed_basis,
            context + ": certified adjoint nonnegative basis");
    std::optional<
        ExactNonnegativePowerSeriesFactorization<ComplexBall>>
        power_series_factorization;
    std::optional<FiniteLaurentFactorization<ComplexBall>>
        laurent_factorization;
    std::string power_series_fallback_reason;
    try {
      power_series_factorization =
          factor_exact_nonnegative_power_series_system(
              certified_nonnegative_transposed_basis,
              context +
                  ": certified adjoint power-series factorization");
    } catch (const MatchingArithmeticError& error) {
      switch (error.code) {
        case MatchingArithmeticErrorCode::AmbiguousZero:
        case MatchingArithmeticErrorCode::ZeroDivisor:
        case MatchingArithmeticErrorCode::SingularOrIncompleteSystem:
        case MatchingArithmeticErrorCode::UnresolvedDeterminantTail:
          break;
        default:
          throw;
      }
      power_series_fallback_reason = error.what();
      laurent_factorization =
          factor_preconditioned_acb_finite_laurent_system(
              std::move(certified_nonnegative_transposed_basis),
              context + ": adjoint Laurent factorization");
    }

    const bool diagnostic_terminal_state =
        environment_flag_is_one("DIFFEXP_DIAGNOSTIC_TERMINAL_STATE");
    if (diagnostic_terminal_state) {
      std::cerr
          << "terminal-adjoint-window context=" << context
          << " factorization="
          << (power_series_factorization.has_value()
                  ? "nonnegative-power-series"
                  : "finite-laurent")
          << " basis_complete=" << conditioned_basis_complete_max
          << " functional_min=" << conditioned_functional_min
          << " incoming=[" << incoming_min << ","
          << incoming_complete_max << "]"
          << " required_output=" << required_output_complete_max
          << " power_series_fallback_reason="
          << (power_series_fallback_reason.empty()
                  ? "none"
                  : power_series_fallback_reason)
          << '\n';
    }

    FiniteLaurentVector<ComplexBall> result;
    result.reserve(conditioned_rows.size());
    for (std::size_t functional = 0;
         functional < conditioned_rows.size(); ++functional) {
      auto adjoint = power_series_factorization.has_value()
          ? solve_factorized_exact_nonnegative_power_series_system(
                *power_series_factorization,
                conditioned_rows[functional],
                context + ": adjoint solve")
          : solve_factorized_finite_laurent_system(
                *laurent_factorization,
                conditioned_rows[functional],
                context + ": adjoint solve");
      AcbLaurentRefinementOptions adjoint_certificate_options;
      adjoint_certificate_options.relative_tolerance =
          Magnitude::decimal(acb->relative_tolerance_text());
      adjoint_certificate_options.required_min_power =
          conditioned_functional_min;
      adjoint_certificate_options.required_complete_max =
          local_algebra_detail::checked_i32(
              static_cast<std::int64_t>(
                  required_output_complete_max) - incoming_min,
              "terminal adjoint certificate maximum");
      adjoint_certificate_options.max_refinement_steps = 0;
      const auto adjoint_certificate =
          matching_detail::evaluate_acb_matching_residual(
              certified_nonnegative_transposed_basis, adjoint,
              conditioned_rows[functional],
              adjoint_certificate_options,
              context + ": adjoint forward certificate");
      if (!adjoint_certificate.diagnostics.complete_through_required)
      {
        const auto additional = std::max<std::int32_t>(
            0,
            adjoint_certificate.diagnostics.required_complete_max -
                adjoint_certificate.diagnostics.complete_window.complete_max);
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            context +
                ": terminal adjoint forward residual does not cover the epsilon prefix needed by the output; common_complete_max=" +
                std::to_string(
                    adjoint_certificate.diagnostics.complete_window.complete_max) +
                "; required_complete_max=" +
                std::to_string(
                    adjoint_certificate.diagnostics.required_complete_max) +
                "; required_additional_epsilon_orders=" +
                std::to_string(additional),
            static_cast<std::size_t>(additional), functional,
            adjoint_certificate.diagnostics.complete_window.complete_max);
      }
      if (adjoint_certificate.diagnostics.verdict !=
          AcbMatchingResidualVerdict::Pass) {
        std::string coefficient_examples;
        std::size_t examples = 0;
        for (const auto& coefficient :
             adjoint_certificate.diagnostics.coefficients) {
          if (coefficient.epsilon_power >
                  adjoint_certificate.diagnostics.required_complete_max ||
              coefficient.verdict == AcbMatchingResidualVerdict::Pass)
            continue;
          if (!coefficient_examples.empty()) coefficient_examples += ",";
          coefficient_examples +=
              "{row=" + std::to_string(coefficient.row) +
              ";power=" + std::to_string(coefficient.epsilon_power) +
              ";verdict=" +
              std::string(
                  coefficient.verdict == AcbMatchingResidualVerdict::Fail
                      ? "fail" : "inconclusive") +
              ";residual_lower=" +
              coefficient.residual_lower.dump_exact() +
              ";residual_upper=" +
              coefficient.residual_upper.dump_exact() +
              ";scale_upper=" + coefficient.scale_upper.dump_exact() +
              "}";
          if (++examples == 4) break;
        }
        throw std::domain_error(
            context +
            ": terminal adjoint forward residual is not certified; verdict=" +
            std::string(
                adjoint_certificate.diagnostics.verdict ==
                        AcbMatchingResidualVerdict::Fail
                    ? "fail" : "inconclusive") +
            "; complete_window=[" +
            std::to_string(
                adjoint_certificate.diagnostics.complete_window.min_power) +
            "," +
            std::to_string(
                adjoint_certificate.diagnostics.complete_window.complete_max) +
            "]; required_complete_max=" +
            std::to_string(
                adjoint_certificate.diagnostics.required_complete_max) +
            "; coefficient_examples=[" + coefficient_examples + "]");
      }
      if (diagnostic_terminal_state) {
        double maximum_residual_upper = 0.0;
        double maximum_scale_upper = 0.0;
        for (const auto& coefficient :
             adjoint_certificate.diagnostics.coefficients) {
          if (coefficient.epsilon_power >
              adjoint_certificate.diagnostics.required_complete_max)
            continue;
          maximum_residual_upper = std::max(
              maximum_residual_upper,
              coefficient.residual_upper.approximate_upper());
          maximum_scale_upper = std::max(
              maximum_scale_upper,
              coefficient.scale_upper.approximate_upper());
        }
        std::cerr
            << "terminal-adjoint-certificate context=" << context
            << " functional=" << functional
            << " residual_upper_max=" << maximum_residual_upper
            << " scale_upper_max=" << maximum_scale_upper
            << " certified_complete_max="
            << adjoint_certificate.diagnostics.required_complete_max
            << '\n';
      }
      if (certified_incoming_tail.has_value()) {
        const auto output_min = local_algebra_detail::checked_i32(
            static_cast<std::int64_t>(conditioned_functional_min) +
                certified_incoming_tail->frame.min_power,
            "terminal certified-tail output minimum");
        const auto output_max = std::min(
            required_output_complete_max,
            local_algebra_detail::checked_i32(
                static_cast<std::int64_t>(
                    adjoint.front().complete_max()) +
                    certified_incoming_tail->frame.complete_max,
                "terminal certified-tail output maximum"));
        for (std::int64_t raw_output = output_min;
             raw_output <= output_max; ++raw_output) {
          const auto output_power = local_algebra_detail::checked_i32(
              raw_output, "terminal certified-tail output power");
          auto tail_output_upper = Magnitude::zero();
          for (std::size_t component = 0;
               component < dimension; ++component) {
            const auto first_adjoint = std::max<std::int64_t>(
                adjoint[component].min_power(),
                raw_output -
                    certified_incoming_tail->frame.complete_max);
            const auto last_adjoint = std::min<std::int64_t>(
                adjoint[component].complete_max(),
                raw_output -
                    certified_incoming_tail->frame.min_power);
            for (std::int64_t raw_adjoint = first_adjoint;
                 raw_adjoint <= last_adjoint; ++raw_adjoint) {
              const auto adjoint_power =
                  local_algebra_detail::checked_i32(
                      raw_adjoint,
                      "terminal certified-tail adjoint power");
              const auto tail_power = raw_output - raw_adjoint;
              const auto tail_row = static_cast<std::size_t>(
                  tail_power -
                      certified_incoming_tail->frame.min_power);
              tail_output_upper += Magnitude::upper_abs(
                  adjoint[component].coefficient(adjoint_power)) *
                  certified_incoming_tail->absolute.at(tail_row);
            }
          }
          std::cerr
              << "terminal-physical-tail-amplification context="
              << context
              << " functional=" << functional
              << " epsilon_power=" << output_power
              << " absolute_upper="
              << tail_output_upper.approximate_upper()
              << " absolute_upper_exact="
              << tail_output_upper.dump_exact()
              << '\n';
        }
      }
      std::optional<EpsilonFrame<ComplexBall>> value;
      for (std::size_t component = 0;
           component < dimension; ++component) {
        auto term = adjoint[component] *
            incoming_value[component];
        if (diagnostic_terminal_state) {
          double adjoint_upper_max = 0.0;
          std::int32_t adjoint_upper_power = adjoint[component].min_power();
          for (std::int64_t raw_power = adjoint[component].min_power();
               raw_power <= adjoint[component].complete_max(); ++raw_power) {
            const auto power = local_algebra_detail::checked_i32(
                raw_power, "terminal adjoint sensitivity power");
            const auto upper = Magnitude::upper_abs(
                adjoint[component].coefficient(power)).approximate_upper();
            if (upper > adjoint_upper_max) {
              adjoint_upper_max = upper;
              adjoint_upper_power = power;
            }
          }
          double contribution_upper_max = 0.0;
          std::int32_t contribution_upper_power = term.min_power();
          const auto contribution_complete = std::min(
              term.complete_max(), required_output_complete_max);
          for (std::int64_t raw_power = term.min_power();
               raw_power <= contribution_complete; ++raw_power) {
            const auto power = local_algebra_detail::checked_i32(
                raw_power, "terminal adjoint contribution power");
            const auto upper = Magnitude::upper_abs(
                term.coefficient(power)).approximate_upper();
            if (upper > contribution_upper_max) {
              contribution_upper_max = upper;
              contribution_upper_power = power;
            }
          }
          std::cerr
              << "terminal-adjoint-sensitivity context=" << context
              << " functional=" << functional
              << " component=" << component
              << " adjoint_upper_max=" << adjoint_upper_max
              << " adjoint_upper_power=" << adjoint_upper_power
              << " contribution_upper_max=" << contribution_upper_max
              << " contribution_upper_power="
              << contribution_upper_power
              << '\n';

          // Attribute the visible interval width of each scalar convolution
          // coefficient to the individual epsilon-product which contributes
          // the largest uncertainty.  This is deliberately diagnostic-only:
          // the untouched ball convolution below remains authoritative.
          for (std::int64_t raw_output = term.min_power();
               raw_output <= contribution_complete; ++raw_output) {
            const auto first_adjoint = std::max<std::int64_t>(
                adjoint[component].min_power(),
                raw_output - incoming_value[component].complete_max());
            const auto last_adjoint = std::min<std::int64_t>(
                adjoint[component].complete_max(),
                raw_output - incoming_value[component].min_power());
            if (first_adjoint > last_adjoint) continue;

            auto dominant_uncertainty = Magnitude::zero();
            auto summed_uncertainty = Magnitude::zero();
            std::optional<std::int32_t> dominant_adjoint_power;
            std::optional<std::int32_t> dominant_incoming_power;
            std::optional<ComplexBall> dominant_product;
            for (std::int64_t raw_adjoint = first_adjoint;
                 raw_adjoint <= last_adjoint; ++raw_adjoint) {
              const auto adjoint_power =
                  local_algebra_detail::checked_i32(
                      raw_adjoint,
                      "terminal adjoint radius-source adjoint power");
              const auto incoming_power =
                  local_algebra_detail::checked_i32(
                      raw_output - raw_adjoint,
                      "terminal adjoint radius-source incoming power");
              const auto product =
                  adjoint[component].coefficient(adjoint_power) *
                  incoming_value[component].coefficient(incoming_power);
              const auto product_midpoint =
                  matching_detail::acb_midpoint_value(product);
              const auto uncertainty = Magnitude::upper_abs(
                  product - product_midpoint);
              summed_uncertainty += uncertainty;
              if (!dominant_adjoint_power.has_value() ||
                  !(uncertainty <= dominant_uncertainty)) {
                dominant_uncertainty = uncertainty;
                dominant_adjoint_power = adjoint_power;
                dominant_incoming_power = incoming_power;
                dominant_product = product;
              }
            }
            if (!dominant_adjoint_power.has_value() ||
                !dominant_incoming_power.has_value() ||
                !dominant_product.has_value())
              continue;

            const auto output_power =
                local_algebra_detail::checked_i32(
                    raw_output,
                    "terminal adjoint radius-source output power");
            const auto& dominant_adjoint =
                adjoint[component].coefficient(
                    *dominant_adjoint_power);
            const auto& dominant_incoming =
                incoming_value[component].coefficient(
                    *dominant_incoming_power);
            const auto& component_output =
                term.coefficient(output_power);
            std::cerr
                << "terminal-adjoint-radius-source context=" << context
                << " functional=" << functional
                << " component=" << component
                << " output_power=" << output_power
                << " adjoint_power=" << *dominant_adjoint_power
                << " incoming_power=" << *dominant_incoming_power
                << " adjoint_midpoint=("
                << dominant_adjoint.real_midpoint(12) << ","
                << dominant_adjoint.imag_midpoint(12) << ")"
                << " adjoint_radius2exp=("
                << dominant_adjoint.real_radius_exponent() << ","
                << dominant_adjoint.imag_radius_exponent() << ")"
                << " incoming_midpoint=("
                << dominant_incoming.real_midpoint(12) << ","
                << dominant_incoming.imag_midpoint(12) << ")"
                << " incoming_radius2exp=("
                << dominant_incoming.real_radius_exponent() << ","
                << dominant_incoming.imag_radius_exponent() << ")"
                << " product_midpoint=("
                << dominant_product->real_midpoint(12) << ","
                << dominant_product->imag_midpoint(12) << ")"
                << " product_radius2exp=("
                << dominant_product->real_radius_exponent() << ","
                << dominant_product->imag_radius_exponent() << ")"
                << " component_output_radius2exp=("
                << component_output.real_radius_exponent() << ","
                << component_output.imag_radius_exponent() << ")"
                << " dominant_uncertainty_upper="
                << dominant_uncertainty.approximate_upper()
                << " dominant_uncertainty_upper_exact="
                << dominant_uncertainty.dump_exact()
                << " summed_uncertainty_upper="
                << summed_uncertainty.approximate_upper()
                << " summed_uncertainty_upper_exact="
                << summed_uncertainty.dump_exact()
                << '\n';
          }
        }
        value = value.has_value()
            ? *value + term
            : std::move(term);
      }
      if (!value.has_value())
        throw std::logic_error(
            context +
            ": terminal adjoint contraction produced no scalar value");
      if (factorized_incoming_transfer.has_value() &&
          !certified_incoming_tail.has_value()) {
        auto factorized_value =
            contract_factorized_ordinary_center_adjoint(
                *factorized_incoming_transfer, adjoint,
                required_output_complete_max, context);
        const auto overlap_min = std::max(
            value->min_power(), factorized_value.min_power());
        const auto overlap_max = std::min({
            value->complete_max(),
            factorized_value.complete_max(),
            required_output_complete_max});
        if (overlap_min > overlap_max ||
            factorized_value.complete_max() <
                required_output_complete_max)
          throw std::logic_error(
              context +
              ": factorized incoming adjoint contraction changed its required epsilon coverage; direct=[" +
              std::to_string(value->min_power()) + "," +
              std::to_string(value->complete_max()) +
              "]; factorized=[" +
              std::to_string(factorized_value.min_power()) + "," +
              std::to_string(factorized_value.complete_max()) +
              "]; required_complete_max=" +
              std::to_string(required_output_complete_max));
        for (std::int64_t raw_power = overlap_min;
             raw_power <= overlap_max; ++raw_power) {
          const auto power =
              local_algebra_detail::checked_i32(
                  raw_power,
                  "factorized incoming adjoint overlap power");
          if (!acb_overlaps(
                  value->coefficient(power).raw(),
                  factorized_value.coefficient(power).raw()))
            throw std::domain_error(
                context +
                ": factorized incoming adjoint contraction is disjoint from the retained direct contraction; epsilon_power=" +
                std::to_string(power));
        }
        if (std::getenv(
                "DIFFEXP_DIAGNOSTIC_FACTORIZED_ORDINARY_EVALUATION") !=
                nullptr ||
            diagnostic_terminal_state)
          std::cerr
              << "factorized-ordinary-adjoint-authority context="
              << context
              << " functional=" << functional
              << " operator_columns="
              << factorized_incoming_transfer->operator_columns
              << " policy=certified-q/C-prefix-replay-and-full-overlap-v1"
              << '\n';
        value = std::move(factorized_value);
      }
      if (value->complete_max() < required_output_complete_max) {
        const auto additional =
            required_output_complete_max - value->complete_max();
        std::ostringstream detail;
        detail
            << context
            << ": terminal adjoint result is incomplete after the actual "
               "factorization and convolution; required_output_complete_max="
            << required_output_complete_max
            << "; actual_output=[" << value->min_power() << ","
            << value->complete_max() << "]"
            << "; factorization="
            << (power_series_factorization.has_value()
                    ? "nonnegative-power-series"
                    : "finite-laurent")
            << "; conditioned_functional_min="
            << conditioned_functional_min
            << "; incoming=[" << incoming_min << ","
            << incoming_complete_max << "]"
            << "; conditioned_basis_complete_max="
            << conditioned_basis_complete_max
            << "; required_additional_epsilon_orders="
            << additional;
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            detail.str(), functional, std::nullopt,
            value->complete_max());
      }
      if (publication_relative_tolerance.has_value())
        for (std::int64_t raw_power = value->min_power();
             raw_power <= std::min<std::int64_t>(
                              value->complete_max(),
                              required_output_complete_max);
             ++raw_power) {
          const auto power = local_algebra_detail::checked_i32(
              raw_power, "terminal adjoint output accuracy power");
          const auto coefficient = value->coefficient(power);
          const auto publication =
              certify_acb_publication_accuracy(
                  coefficient, *publication_relative_tolerance);
          if (!publication.acceptable)
            throw MatchingArithmeticError(
                MatchingArithmeticErrorCode::
                    TerminalOutputInconclusive,
                context +
                    ": terminal adjoint output ball is too wide to publish; functional=" +
                    std::to_string(functional) +
                    "; epsilon_power=" + std::to_string(power) +
                    "; midpoint=(" + coefficient.real_midpoint(16) + "," +
                    coefficient.imag_midpoint(16) + ")" +
                    "; radius2exp=(" +
                    coefficient.real_radius_exponent() + "," +
                    coefficient.imag_radius_exponent() + ")",
                functional, std::nullopt, power,
                publication.required_additional_digits);
        }
      if (diagnostic_terminal_state) {
        auto adjoint_min = adjoint.front().min_power();
        auto adjoint_complete = adjoint.front().complete_max();
        for (const auto& frame : adjoint) {
          adjoint_min = std::min(adjoint_min, frame.min_power());
          adjoint_complete =
              std::min(adjoint_complete, frame.complete_max());
        }
        std::cerr
            << "terminal-adjoint-result context=" << context
            << " functional=" << functional
            << " functional_window=["
            << conditioned_rows[functional].front().min_power()
            << ","
            << conditioned_rows[functional].front().complete_max()
            << "] adjoint=[" << adjoint_min << ","
            << adjoint_complete << "] output=["
            << value->min_power() << ","
            << value->complete_max() << "]\n";
      }
      result.push_back(std::move(*value));
    }
    return result;
  }

  std::int32_t terminal_acb_factorized_input_min_power() const {
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    if (!acb || !has_terminal_acb_factorization())
      throw std::invalid_argument(
          "terminal factorized input valuation requires one certified Acb exact-right match");
    std::vector<std::int64_t> input_minima(
        basis_owners_.size(), 0);
    const auto propagate_exact_right =
        [&](const auto& matrix, const char* context) {
          if (matrix.size() != input_minima.size() ||
              matrix.empty() || matrix.front().empty())
            throw std::logic_error(
                std::string(context) +
                " differs from the terminal physical basis dimension");
          const auto output_columns = matrix.front().size();
          std::vector<std::optional<std::int64_t>> output_minima(
              output_columns);
          for (std::size_t input = 0; input < matrix.size();
               ++input) {
            if (matrix[input].size() != output_columns)
              throw std::logic_error(
                  std::string(context) + " is not rectangular");
            for (std::size_t output = 0;
                 output < output_columns; ++output) {
              const auto factor_minimum =
                  matrix[input][output].minimum_power();
              if (!factor_minimum.has_value()) continue;
              const auto candidate =
                  input_minima[input] + *factor_minimum;
              if (!output_minima[output].has_value() ||
                  candidate < *output_minima[output])
                output_minima[output] = candidate;
            }
          }
          std::vector<std::int64_t> propagated;
          propagated.reserve(output_columns);
          for (const auto minimum : output_minima) {
            if (!minimum.has_value())
              throw std::logic_error(
                  std::string(context) +
                  " has a structurally zero output column");
            propagated.push_back(*minimum);
          }
          input_minima = std::move(propagated);
        };
    if (const auto& normal_right =
            acb->terminal_normal_frame_right_transformation();
        normal_right.has_value()) {
      const auto propagate_finite_right =
          [&](const FiniteLaurentMatrix<ComplexBall>& matrix,
              const char* context) {
            if (matrix.size() != input_minima.size() ||
                matrix.empty() || matrix.front().empty())
              throw std::logic_error(
                  std::string(context) +
                  " differs from the terminal physical basis dimension");
            const auto output_columns = matrix.front().size();
            std::vector<std::optional<std::int64_t>> output_minima(
                output_columns);
            for (std::size_t input = 0; input < matrix.size(); ++input) {
              if (matrix[input].size() != output_columns)
                throw std::logic_error(
                    std::string(context) + " is not rectangular");
              for (std::size_t output = 0;
                   output < output_columns; ++output) {
                const auto& frame = matrix[input][output];
                std::optional<std::int32_t> factor_minimum;
                for (std::int64_t raw_power = frame.min_power();
                     raw_power <= frame.complete_max(); ++raw_power) {
                  const auto power =
                      static_cast<std::int32_t>(raw_power);
                  if (!frame.coefficient(power).is_zero()) {
                    factor_minimum = power;
                    break;
                  }
                }
                if (!factor_minimum.has_value()) continue;
                const auto candidate =
                    input_minima[input] + *factor_minimum;
                if (!output_minima[output].has_value() ||
                    candidate < *output_minima[output])
                  output_minima[output] = candidate;
              }
            }
            std::vector<std::int64_t> propagated;
            propagated.reserve(output_columns);
            for (const auto minimum : output_minima) {
              if (!minimum.has_value())
                throw std::logic_error(
                    std::string(context) +
                    " has a structurally zero output column");
              propagated.push_back(*minimum);
            }
            input_minima = std::move(propagated);
          };
      propagate_finite_right(
          *normal_right,
          "terminal finite normal-frame right transformation");
      if (const auto& acb_right =
              acb->acb_materialization_right_transformation();
          acb_right.has_value())
        propagate_exact_right(
            *acb_right,
            "terminal normal-frame Acb exact-support transformation");
      else {
        const auto exact_right =
            specialize_exact_rational_laurent_matrix_to_acb(
                acb->exact_saturation_transformation());
        propagate_exact_right(
            exact_right,
            "terminal normal-frame exact-right transformation");
      }
    } else {
      propagate_exact_right(
          *acb->acb_materialization_right_transformation(),
          "terminal exact-right transformation");
    }
    if (const auto& preconditioner =
            acb->acb_right_materialization_preconditioner();
        preconditioner.has_value())
      propagate_exact_right(
          *preconditioner,
          "terminal exact-right Acb preconditioner");

    std::optional<std::int64_t> result;
    for (const auto minimum : input_minima)
      if (!result.has_value() || minimum < *result)
        result = minimum;
    if (!result.has_value())
      throw std::logic_error(
          "terminal factorized input valuation has no structural path");
    return local_algebra_detail::checked_i32(
        *result, "terminal factorized input epsilon shift");
  }

  std::int32_t terminal_acb_transformed_weight_min_power() const {
    const auto& transformed_weights = terminal_acb_transformed_weights();
    if (transformed_weights.empty())
      throw std::logic_error(
          "terminal factorized contraction has no transformed weights");
    auto result = transformed_weights.front().min_power();
    for (const auto& weight : transformed_weights)
      result = std::min(result, weight.min_power());
    return result;
  }

  FiniteLaurentVector<ComplexBall>
  contract_terminal_acb_functionals(
      const FiniteLaurentMatrix<ComplexBall>& physical_rows,
      const std::string& context) const {
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    if (!acb || !acb->certified_for_materialization())
      throw std::invalid_argument(
          context +
          ": terminal functional contraction requires one certified Acb match");
    const auto& right =
        acb->acb_materialization_right_transformation();
    if (!right.has_value())
      throw std::invalid_argument(
          context +
          ": terminal functional contraction has no certified right factorization");
    if (physical_rows.empty())
      throw std::invalid_argument(
          context + ": terminal functional batch is empty");
    for (const auto& row : physical_rows)
      if (row.size() != basis_owners_.size())
        throw std::invalid_argument(
            context +
            ": terminal functional row differs from the physical basis dimension");

    auto transformed =
        right_multiply_finite_by_exact_laurent(
            physical_rows, *right);
    if (const auto& preconditioner =
            acb->acb_right_materialization_preconditioner();
        preconditioner.has_value())
      transformed = right_multiply_finite_by_exact_laurent(
          transformed, *preconditioner);
    return apply_finite_laurent_matrix(
        transformed, acb->transformed_weights());
  }

  std::optional<json::object> incomplete_acb_summary() const {
    const auto acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
    if (!acb || acb->certified_for_materialization()) return std::nullopt;
    return acb->summary();
  }

  std::shared_ptr<PhysicalEquationOwnerBase>
  inheritable_basis_equation_owner() const {
    if (basis_owners_.empty()) return nullptr;
    auto result = basis_owners_.front()->retained_equation_owner();
    for (const auto& basis : basis_owners_) {
      const auto candidate = basis->retained_equation_owner();
      if (candidate.get() != result.get()) return nullptr;
    }
    return result;
  }

  std::shared_ptr<PhysicalEquationOwnerBase>
  planned_receiving_equation_owner() const {
    if (!plan_owner_) return nullptr;
    const auto arm_name = required_string(handoff_, "arm");
    const auto match_index = static_cast<std::size_t>(
        as_u64(handoff_.at("match"),
               "planned-hop receiving equation-owner match index"));
    const auto& arm = plan_owner_->arm(arm_name);
    if (match_index >= arm.exact.matches.size())
      throw std::logic_error(
          "planned-hop receiving equation owner lies outside its retained arm");
    const auto& binding = arm.charts.at(
        arm.exact.matches[match_index].receiving_chart);
    return std::visit(
        [](const auto& owner)
            -> std::shared_ptr<PhysicalEquationOwnerBase> {
          return std::static_pointer_cast<PhysicalEquationOwnerBase>(
              owner);
        },
        binding.owner);
  }

  static std::string planned_equation_owner_signature_binding(
      const PhysicalEquationOwnerBase& owner) {
    if (!owner.owner_signature_identity().empty())
      return owner.owner_signature_identity();
    const auto& identity = owner.equation_operator_identity();
    return "owner-operator-reference-v1:" +
        public_provenance_fingerprint(identity) + ":" +
        std::to_string(identity.size());
  }

  static std::string planned_equation_payload_binding(
      const PhysicalEquationOwnerBase& owner) {
    return owner.physical_payload_identity().empty()
        ? "no-physical-qc-payload-v1"
        : owner.physical_payload_identity();
  }

  void validate_materialized_equation_owner(
      const std::shared_ptr<StoredLocalBase>& local) const {
    if (!local || local->retained_equation_owner().get() !=
                      inheritable_basis_equation_owner().get())
      throw std::invalid_argument(
          "checkpoint materialized-local equation owner differs from its planned-hop basis");
  }

  void validate_materialized_derivation(
      const json::object& derivation, const char* scalar_domain) const {
    const auto native_summary = match_->summary();
    const auto derivation_schema = required_string(derivation, "schema");
    if (derivation_schema ==
        "diffexp3-retained-plan-match-local-materialization-v2") {
      const auto equation_owner = inheritable_basis_equation_owner();
      const auto planned_equation_owner =
          planned_receiving_equation_owner();
      const auto& producing = as_object(
          handoff_.at("producing"),
          "checkpoint compact planned-hop producing record");
      const auto& receiving = as_object(
          handoff_.at("receiving"),
          "checkpoint compact planned-hop receiving record");
      if (!plan_owner_ || !planned_equation_owner ||
          (equation_owner &&
           equation_owner.get() != planned_equation_owner.get()) ||
          required_string(derivation, "source_match") != handle_ ||
          required_string(
              derivation, "source_match_checkpoint_identity") !=
              checkpoint_identity_ ||
          required_string(derivation, "tile_plan") !=
              plan_owner_->handle() ||
          required_string(
              derivation, "tile_plan_checkpoint_identity") !=
              plan_owner_->checkpoint_identity() ||
          derivation.at("arm") != handoff_.at("arm") ||
          derivation.at("match") != handoff_.at("match") ||
          derivation.at("incoming") != producing.at("incoming") ||
          derivation.at("basis") != receiving.at("basis") ||
          required_string(
              derivation, "equation_owner_signature_identity") !=
              planned_equation_owner_signature_binding(
                  *planned_equation_owner) ||
          required_string(derivation, "equation_payload_identity") !=
              planned_equation_payload_binding(
                  *planned_equation_owner))
        throw std::invalid_argument(
            "checkpoint compact materialized-local lineage disagrees with its live plan, match, or equation owner");
    } else if (derivation_schema ==
               "diffexp3-retained-plan-match-local-materialization-v1") {
      if (required_string(derivation, "source_match") != handle_ ||
          required_string(
              derivation, "source_match_checkpoint_identity") !=
              checkpoint_identity_ ||
          required_string(
              derivation, "source_match_provenance_identity") !=
              required_string(native_summary, "provenance_identity") ||
          required_string(
              derivation, "planned_hop_provenance_identity") !=
              provenance_identity_ ||
          derivation.at("planned_hop") != handoff_)
        throw std::invalid_argument(
            "checkpoint materialized-local lineage disagrees with its planned-hop owner");
    } else {
      throw std::invalid_argument(
          "checkpoint materialized-local lineage has an unsupported schema");
    }

    json::array expected_windows;
    std::int32_t expected_certified_max = 0;
    if (const auto exact =
            std::dynamic_pointer_cast<StoredExactRegularMatch>(match_)) {
      if (std::string(scalar_domain) != "rational")
        throw std::invalid_argument(
            "checkpoint materialized-local scalar domain differs from its exact match owner");
      for (const auto& weight : exact->weights())
        expected_windows.push_back(json::object{
            {"min", weight.min_power()}, {"max", weight.complete_max()}});
      expected_certified_max = as_i32(
          as_object(native_summary.at("residual"),
                    "exact retained match residual").at("max"),
          "exact retained match residual maximum");
    } else if (const auto acb =
                   std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_)) {
      if (std::string(scalar_domain) != "acb" ||
          !acb->certified_for_materialization())
        throw std::invalid_argument(
            "checkpoint materialized-local Acb owner lost its passing complete match certificate");
      for (const auto& weight : acb->weights())
        expected_windows.push_back(json::object{
            {"min", weight.min_power()}, {"max", weight.complete_max()}});
      expected_certified_max = acb->certified_complete_max();
    } else {
      throw std::invalid_argument(
          "checkpoint materialized-local owner embeds an unsupported native match");
    }
    if (derivation.at("weight_windows") != expected_windows ||
        as_i32(derivation.at("match_certified_complete_max"),
               "materialized-local certified maximum") !=
            expected_certified_max)
      throw std::invalid_argument(
          "checkpoint materialized-local derivation differs from its retained match weights/certificate");
  }

  std::shared_ptr<StoredLocalBase> materialize(
      const std::string& local_handle,
      const std::string& result_checkpoint_identity,
      slong precision_bits,
      const std::shared_ptr<StoredPlannedMatchHop>& self,
      bool allow_terminal_factorized_proxy = false) {
    if (self.get() != this)
      throw std::logic_error(
          "retained plan-match materialization lost self ownership");
    std::shared_ptr<StoredLocalBase> result;
    if (const auto exact =
            std::dynamic_pointer_cast<StoredExactRegularMatch>(match_)) {
      result = materialize_typed<Rational>(
          local_handle, result_checkpoint_identity, precision_bits,
          exact->weights(), self, allow_terminal_factorized_proxy);
    } else if (const auto acb =
                   std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_)) {
      if (!acb->certified_for_materialization())
        {
          const auto summary = acb->summary();
          const auto& lattice = as_object(
              summary.at("exact_lattice"),
              "failed Acb match exact lattice summary");
          const json::object compact_lattice{
              {"normalized_determinant_valuation",
               lattice.at("normalized_determinant_valuation")},
              {"transformation_min_power",
               lattice.at("transformation_min_power")},
              {"transformation_terms",
               lattice.at("transformation_terms")},
              {"initial_column_shifts",
               lattice.at("initial_column_shifts")},
              {"initial_leading_rank",
               lattice.at("initial_leading_rank")},
              {"final_leading_rank",
               lattice.at("final_leading_rank")}};
          throw std::domain_error(
              "an Acb plan-match handoff must have a passing complete residual before materialization; residual=" +
              json::serialize(summary.at("residual")) +
              "; refinement=" +
              json::serialize(summary.at("refinement")) +
              "; epsilon=" +
              json::serialize(summary.at("epsilon")) +
              "; weight_windows=" +
              json::serialize(summary.at("weight_windows")) +
              "; exact_lattice=" +
              json::serialize(compact_lattice));
        }
      result = materialize_typed<ComplexBall>(
          local_handle, result_checkpoint_identity, precision_bits,
          acb->weights(), self, allow_terminal_factorized_proxy);
    } else {
      throw std::logic_error(
          "retained plan-match handoff has an unsupported matching state");
    }
    materializations_.fetch_add(1);
    return result;
  }

  json::object checkpoint_record() const override {
    return json::object{
        {"schema", "diffexp3-retained-planned-match-hop-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"handoff", handoff_},
        {"native_match", match_->checkpoint_record()},
        {"elapsed_ms", elapsed_ms_},
        {"runtime_stats",
         json::object{{"materializations", materializations_.load()}}}};
  }

  void restore_runtime_stats(std::uint64_t materializations) {
    materializations_.store(materializations);
  }

 private:
  template <typename Scalar>
  std::shared_ptr<StoredLocalBase> materialize_typed(
      const std::string& local_handle,
      const std::string& result_checkpoint_identity,
      slong precision_bits,
      const FiniteLaurentVector<Scalar>& weights,
      const std::shared_ptr<StoredPlannedMatchHop>& self,
      bool allow_terminal_factorized_proxy) const {
    if (local_handle.empty() || result_checkpoint_identity.empty())
      throw std::invalid_argument(
          "plan-match local materialization identities must be nonempty");
    AcbPrecisionLease lease(precision_bits);
    ComplexBall::set_precision(precision_bits);
    std::optional<matching_detail::ScopedAcbPrecision>
        exact_shadow_precision;
    if constexpr (std::is_same_v<Scalar, ComplexBall>) {
      const auto acb =
          std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
      if (acb && acb->exact_shadow_factorized_basis().has_value())
        exact_shadow_precision.emplace(
            acb->exact_shadow_extra_precision_bits());
    }
    const auto started = std::chrono::steady_clock::now();
    const auto native_match_summary = match_->summary();
    const auto& match_epsilon = as_object(
        native_match_summary.at("epsilon"),
        "retained match epsilon provenance");
    const auto match_work_max = as_i32(
        match_epsilon.at("max"), "retained match epsilon maximum");
    const auto match_certified_max = [&]() {
      if constexpr (std::is_same_v<Scalar, Rational>)
        return as_i32(
            as_object(native_match_summary.at("residual"),
                      "exact retained match residual").at("max"),
            "exact retained match residual maximum");
      else {
        const auto acb =
            std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
        if (!acb)
          throw std::logic_error(
              "Acb materialization lost its retained match owner");
        return acb->certified_complete_max();
      }
    }();

    std::vector<std::shared_ptr<StoredLocal<Scalar>>> typed_basis;
    std::vector<const LocalSolution<Scalar>*> basis_solutions;
    typed_basis.reserve(basis_owners_.size());
    basis_solutions.reserve(basis_owners_.size());
    std::int32_t materialized_top = kCompleteInfinity;
    const auto receiving_chart = basis_owners_.front()->source_chart();
    const auto receiving_operator =
        basis_owners_.front()->source_operator_identity();
    auto equation_owner = inheritable_basis_equation_owner();
    const auto planned_equation_owner =
        planned_receiving_equation_owner();
    for (std::size_t column = 0; column < basis_owners_.size(); ++column) {
      auto typed =
          std::dynamic_pointer_cast<StoredLocal<Scalar>>(basis_owners_[column]);
      if (!typed)
        throw std::logic_error(
            "retained plan-match basis coefficient domain changed");
      if (typed->source_chart() != receiving_chart ||
          typed->source_operator_identity() != receiving_operator)
        throw std::logic_error(
            "retained plan-match basis chart provenance changed");
      if (typed->top_valid() < match_work_max)
        throw std::domain_error(
            "retained plan-match basis validity does not cover its matching work window");
      if (column >= weights.size())
        throw std::logic_error(
            "retained plan-match weight count is smaller than its basis");
      const auto basis_valid = std::min(
          typed->top_valid(), typed->solution().epsilon.complete_max);
      const auto shifted_basis_valid = local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(basis_valid) +
              weights[column].min_power(),
          "materialized local top validity");
      const auto weight_valid = local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(typed->solution().epsilon.min_power) +
              weights[column].complete_max(),
          "materialized weight top validity");
      materialized_top = std::min(
          materialized_top, std::min(shifted_basis_valid, weight_valid));
      basis_solutions.push_back(&typed->solution());
      typed_basis.push_back(std::move(typed));
    }
    if (weights.size() != typed_basis.size())
      throw std::logic_error(
          "retained plan-match weight count differs from its basis");
    const auto typed_incoming =
        std::dynamic_pointer_cast<StoredLocal<Scalar>>(incoming_owner_);
    if (!typed_incoming || typed_incoming->top_valid() < match_work_max)
      throw std::domain_error(
          "retained plan-match incoming validity does not cover its matching work window");

    std::optional<std::vector<LocalSolution<ComplexBall>>>
        exact_right_basis;
    std::vector<const LocalSolution<Scalar>*> materialization_basis =
        basis_solutions;
    const FiniteLaurentVector<Scalar>* materialization_weights = &weights;
    bool exact_right_physical_association = false;
    bool exact_right_used_rational_shadow = false;
    std::optional<std::uint32_t> matching_taylor_complete_max;
    std::optional<std::uint32_t>
        incoming_matching_taylor_complete_max;
    std::optional<std::vector<LocalSolution<ComplexBall>>>
        matching_prefix_basis;
    if constexpr (std::is_same_v<Scalar, ComplexBall>) {
      const auto acb =
          std::dynamic_pointer_cast<StoredRefinedAcbMatch>(match_);
      if (!acb)
        throw std::logic_error(
            "Acb materialization lost its refined match owner");
      matching_taylor_complete_max =
          acb->matching_taylor_complete_max();
      incoming_matching_taylor_complete_max =
          acb->incoming_matching_taylor_complete_max();
      if (matching_taylor_complete_max.has_value()) {
        matching_prefix_basis.emplace();
        matching_prefix_basis->reserve(basis_solutions.size());
        for (std::size_t column = 0;
             column < basis_solutions.size(); ++column)
          matching_prefix_basis->push_back(
              restrict_local_taylor_prefix(
                  *basis_solutions[column],
                  *matching_taylor_complete_max,
                  result_checkpoint_identity +
                      ":matching-taylor-prefix:column:" +
                      std::to_string(column)));
        materialization_basis.clear();
        materialization_basis.reserve(
            matching_prefix_basis->size());
        for (const auto& column : *matching_prefix_basis)
          materialization_basis.push_back(&column);
      }
      if (const auto& exact_shadow =
              acb->exact_shadow_factorized_basis();
          exact_shadow.has_value()) {
        // The final normal-frame match was certified against the exact
        // physical Rational-shadow columns F*T. Reusing F and the separate
        // bounded V/V^-1 maps
        // would recreate the catastrophic cancellation which the exact
        // coefficientwise association removed.
        exact_right_basis = *exact_shadow;
        if (const auto& preconditioner =
                acb->acb_right_materialization_preconditioner();
            preconditioner.has_value()) {
          std::vector<const LocalSolution<ComplexBall>*>
              exact_shadow_sources;
          exact_shadow_sources.reserve(exact_right_basis->size());
          for (const auto& column : *exact_right_basis)
            exact_shadow_sources.push_back(&column);
          exact_right_basis =
              right_transform_local_basis_exact<ComplexBall>(
                  exact_shadow_sources, *preconditioner,
                  result_checkpoint_identity +
                      ":conditioned-exact-shadow-basis");
        }
        materialization_basis.clear();
        materialization_basis.reserve(exact_right_basis->size());
        for (const auto& column : *exact_right_basis)
          materialization_basis.push_back(&column);
        materialization_weights = &acb->transformed_weights();
        exact_right_physical_association = true;
        exact_right_used_rational_shadow = true;
      } else if (const auto& right =
              acb->acb_materialization_right_transformation();
          right.has_value() &&
          !acb->terminal_normal_frame_right_transformation().has_value()) {
        // The certified match was solved against these retained Acb columns,
        // so they remain the sole materialization authority.  Preserve the
        // certified right association whether T came from an exact Rational
        // shadow or from the retained singular-SCC Acb lattice witness.
        // A terminal normal-frame match is different: its Acb right acts
        // after the retained finite R frame.  Applying it directly to the
        // physical basis would form F*(R^-1*T), not the certified F*w.
        // The ordinary LocalSolution proxy therefore uses physical weights;
        // terminal contraction replays R and R^-1*T sequentially.
        // Private pole rows must remain attached to epsilon-dependent
        // Frobenius sectors: they can feed finite and positive orders when
        // the local is evaluated, even when the final observable is regular.
        exact_right_basis =
            right_transform_local_basis_exact<ComplexBall>(
                materialization_basis, *right,
                result_checkpoint_identity +
                    ":exact-right-basis");
        if (const auto& preconditioner =
                acb->acb_right_materialization_preconditioner();
            preconditioner.has_value()) {
          std::vector<const LocalSolution<ComplexBall>*>
              exact_right_basis_sources;
          exact_right_basis_sources.reserve(exact_right_basis->size());
          for (const auto& column : *exact_right_basis)
            exact_right_basis_sources.push_back(&column);
          exact_right_basis =
              right_transform_local_basis_exact<ComplexBall>(
                  exact_right_basis_sources, *preconditioner,
                  result_checkpoint_identity +
                      ":conditioned-exact-right-basis");
        }
        materialization_basis.clear();
        materialization_basis.reserve(exact_right_basis->size());
        for (const auto& column : *exact_right_basis)
          materialization_basis.push_back(&column);
        materialization_weights = &acb->transformed_weights();
        exact_right_physical_association = true;
        exact_right_used_rational_shadow =
            acb->exact_right_materialization_transformation().has_value();
      }
    }

    // Re-evaluate the honest epsilon edge in the association that will
    // actually be materialized.  In particular, transformed weights are not
    // interchangeable with the physical F,w windows: exact saturation and
    // the constant right frame can shift both Laurent edges.  Retain the
    // earlier physical cap as an additional conservative bound so
    // zero-extended proposal coefficients can never become publishable.
    std::int32_t associated_materialized_top = kCompleteInfinity;
    if (materialization_basis.size() != materialization_weights->size())
      throw std::logic_error(
          "selected materialization basis and weight counts differ");
    for (std::size_t column = 0;
         column < materialization_basis.size(); ++column) {
      const auto& local = *materialization_basis[column];
      const auto& weight = (*materialization_weights)[column];
      const auto shifted_basis_valid =
          local_algebra_detail::checked_i32(
              static_cast<std::int64_t>(
                  local.epsilon.complete_max) +
                  weight.min_power(),
              "associated materialized local top validity");
      const auto weight_valid =
          local_algebra_detail::checked_i32(
              static_cast<std::int64_t>(local.epsilon.min_power) +
                  weight.complete_max(),
              "associated materialized weight top validity");
      associated_materialized_top = std::min(
          associated_materialized_top,
          std::min(shifted_basis_valid, weight_valid));
    }
    materialized_top = std::min(
        materialized_top, associated_materialized_top);

    auto solution = materialize_local_basis_weights(
        materialization_basis, *materialization_weights,
        result_checkpoint_identity);
    materialized_top = std::min(materialized_top,
                                solution.epsilon.complete_max);
    if constexpr (std::is_same_v<Scalar, Rational>) {
      // For exact arithmetic the residual edge is also the algebraic finite
      // Laurent validity edge.
      materialized_top = std::min(materialized_top, match_certified_max);
    } else {
      // An Acb residual pass is an accuracy certificate, not a proof that
      // coefficients above that prefix are invalid.  The factorized ball
      // solve still encloses every coefficient in its honest finite work
      // frame.  Retain that private reservoir so later epsilon losses can
      // consume it; its radii remain visible and the eventual public
      // residual/output checks still fail closed if it is too imprecise.
    }
    if (materialized_top < solution.epsilon.min_power)
      throw std::domain_error(
          "plan-match materialization has no valid output epsilon coefficient");
    if (materialized_top < solution.epsilon.complete_max)
      solution = restrict_local_epsilon_frame_strict_lower(
          solution, solution.epsilon.min_power, materialized_top,
          result_checkpoint_identity);
    if constexpr (std::is_same_v<Scalar, ComplexBall>) {
      const auto arm_name = required_string(handoff_, "arm");
      const auto match_index = static_cast<std::size_t>(
          as_u64(handoff_.at("match"),
                 "materialized continuity match index"));
      const auto& retained_arm = plan_owner_->arm(arm_name);
      if (match_index >= retained_arm.certified_matches.size())
        throw std::logic_error(
            "materialized continuity match lies outside its certified plan");
      const auto& certified_match =
          retained_arm.certified_matches[match_index];
      const auto& exact_match =
          retained_arm.exact.matches.at(match_index);
      const auto producing_rim = exact_plan_rim(
          retained_arm.charts.at(exact_match.producing_chart)
              .prescriptions,
          retained_arm.charts.at(exact_match.producing_chart)
              .geometry.scale);
      const auto receiving_rim = exact_plan_rim(
          retained_arm.charts.at(exact_match.receiving_chart)
              .prescriptions,
          retained_arm.charts.at(exact_match.receiving_chart)
              .geometry.scale);
      auto continuity_incoming = typed_incoming->solution();
      if (incoming_matching_taylor_complete_max.has_value())
        continuity_incoming = restrict_local_taylor_prefix(
            continuity_incoming,
            *incoming_matching_taylor_complete_max,
            result_checkpoint_identity +
                ":continuity-incoming-prefix");
      EvaluationOptions receiving_options;
      receiving_options.imaginary_sign = receiving_rim;
      receiving_options.compute_tail_estimate = false;
      EvaluationOptions incoming_options;
      incoming_options.imaginary_sign = producing_rim;
      incoming_options.compute_tail_estimate = false;
      const auto receiving_evaluation = evaluate_local_solution(
          solution, certified_match.receiving.local,
          receiving_options);
      const auto incoming_evaluation = evaluate_local_solution(
          continuity_incoming, certified_match.producing.local,
          incoming_options);
      const auto& summary_epsilon = as_object(
          native_match_summary.at("epsilon"),
          "materialized continuity match epsilon");
      const auto required_continuity_max = as_i32(
          summary_epsilon.at("required_complete_max"),
          "materialized continuity required epsilon maximum");
      const auto continuity_min = std::max({
          as_i32(summary_epsilon.at("min"),
                 "materialized continuity epsilon minimum"),
          receiving_evaluation.value.epsilon.min_power,
          incoming_evaluation.value.epsilon.min_power});
      const auto continuity_max = std::min({
          match_certified_max,
          receiving_evaluation.value.epsilon.complete_max,
          incoming_evaluation.value.epsilon.complete_max});
      if (continuity_max < required_continuity_max)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            "materialized receiving local does not cover its certified "
            "handoff continuity window: required_max=" +
                std::to_string(required_continuity_max) +
                "; overlap=[" + std::to_string(continuity_min) + "," +
                std::to_string(continuity_max) + "]" +
                "; receiving=[" +
                std::to_string(
                    receiving_evaluation.value.epsilon.min_power) +
                "," +
                std::to_string(
                    receiving_evaluation.value.epsilon.complete_max) +
                "]; incoming=[" +
                std::to_string(
                    incoming_evaluation.value.epsilon.min_power) +
                "," +
                std::to_string(
                    incoming_evaluation.value.epsilon.complete_max) +
                "]; match_certified_max=" +
                std::to_string(match_certified_max),
            std::nullopt, std::nullopt, continuity_max);
      if (continuity_min > required_continuity_max)
        throw std::domain_error(
            "materialized receiving local does not cover its certified "
            "handoff continuity window: required_max=" +
            std::to_string(required_continuity_max) +
            "; overlap=[" + std::to_string(continuity_min) + "," +
            std::to_string(continuity_max) + "]" +
            "; receiving=[" +
            std::to_string(receiving_evaluation.value.epsilon.min_power) +
            "," +
            std::to_string(receiving_evaluation.value.epsilon.complete_max) +
            "]; incoming=[" +
            std::to_string(incoming_evaluation.value.epsilon.min_power) +
            "," +
            std::to_string(incoming_evaluation.value.epsilon.complete_max) +
            "]; match_certified_max=" +
            std::to_string(match_certified_max));
      const auto& refinement_record = as_object(
          native_match_summary.at("refinement"),
          "materialized continuity refinement");
      const auto tolerance = Magnitude::decimal(
          required_string(refinement_record, "relative_tolerance"));
      // The match certifies the requested handoff prefix, while higher
      // coefficients can be a deliberately speculative reservoir retained
      // for downstream epsilon losses.  Requiring those private coefficients
      // to satisfy the public residual tolerance turns useful wide balls into
      // a false continuity failure.  Coverage above proves that the requested
      // edge is present; validate exactly through that certified edge.
      bool continuity_authoritative = true;
      for (std::int64_t raw_power = continuity_min;
           raw_power <= required_continuity_max &&
           continuity_authoritative; ++raw_power) {
        const auto power = static_cast<std::int32_t>(raw_power);
        for (std::uint32_t component = 0;
             component < solution.dimension; ++component) {
          const auto& expected =
              incoming_evaluation.value.at(power, component);
          const auto continuity =
              certify_acb_handoff_continuity(
                  receiving_evaluation.value.at(power, component),
                  expected, tolerance);
          if (!continuity.acceptable) {
            if (allow_terminal_factorized_proxy &&
                exact_right_physical_association) {
              // A terminal consumer keeps the certified factorization
              // ((L*F)*T)*P*y intact.  The coefficientwise LocalSolution is
              // retained only as a geometry/provenance proxy: interval
              // wrapping in that deliberately unused association must not
              // reject an otherwise certified terminal handoff.
              continuity_authoritative = false;
              break;
            }
            const auto& normal_frame = as_object(
                native_match_summary.at("normal_frame_attempt"),
                "materialized continuity normal-frame attempt");
            const auto compact_ball = [](const ComplexBall& value) {
              return json::object{
                  {"real_midpoint", value.real_midpoint(20)},
                  {"imag_midpoint", value.imag_midpoint(20)},
                  {"real_radius_exponent",
                   value.real_radius_exponent()},
                  {"imag_radius_exponent",
                   value.imag_radius_exponent()}};
            };
            json::object reconstruction_routes;
            try {
              const auto acb_match =
                  std::dynamic_pointer_cast<StoredRefinedAcbMatch>(
                      match_);
              if (!acb_match)
                throw std::logic_error(
                    "continuity route audit lost its Acb match");
              std::vector<const LocalSolution<ComplexBall>*>
                  physical_matching_basis;
              if (matching_prefix_basis.has_value()) {
                physical_matching_basis.reserve(
                    matching_prefix_basis->size());
                for (const auto& column : *matching_prefix_basis)
                  physical_matching_basis.push_back(&column);
              } else {
                physical_matching_basis = basis_solutions;
              }
              const auto audit_route = [&](
                  const LocalSolution<ComplexBall>& candidate) {
                const auto evaluated = evaluate_local_solution(
                    candidate, certified_match.receiving.local,
                    receiving_options);
                const auto& value =
                    evaluated.value.at(power, component);
                const auto route_residual = value - expected;
                return json::object{
                    {"value", compact_ball(value)},
                    {"residual_upper",
                     Magnitude::upper_abs(route_residual)
                         .approximate_upper()},
                    {"epsilon",
                     json::object{
                         {"min", candidate.epsilon.min_power},
                         {"max",
                          candidate.epsilon.complete_max}}}};
              };

              const auto physical_candidate =
                  materialize_local_basis_weights(
                      physical_matching_basis,
                      acb_match->weights(),
                      result_checkpoint_identity +
                          ":continuity-audit:physical");
              reconstruction_routes["physical-F*w"] =
                  audit_route(physical_candidate);

              const auto& exact_right =
                  acb_match
                      ->exact_right_materialization_transformation();
              if (exact_right.has_value()) {
                const auto numeric_right =
                    specialize_exact_rational_laurent_matrix_to_acb(
                        *exact_right);
                auto numeric_right_basis =
                    right_transform_local_basis_exact<ComplexBall>(
                        physical_matching_basis, numeric_right,
                        result_checkpoint_identity +
                            ":continuity-audit:numeric-F*T");
                if (const auto& preconditioner =
                        acb_match
                            ->acb_right_materialization_preconditioner();
                    preconditioner.has_value()) {
                  std::vector<
                      const LocalSolution<ComplexBall>*>
                      numeric_right_sources;
                  numeric_right_sources.reserve(
                      numeric_right_basis.size());
                  for (const auto& column :
                       numeric_right_basis)
                    numeric_right_sources.push_back(&column);
                  auto numeric_sequential_basis =
                      right_transform_local_basis_exact<ComplexBall>(
                          numeric_right_sources,
                          *preconditioner,
                          result_checkpoint_identity +
                              ":continuity-audit:numeric-(F*T)*P");
                  std::vector<
                      const LocalSolution<ComplexBall>*>
                      numeric_sequential_sources;
                  numeric_sequential_sources.reserve(
                      numeric_sequential_basis.size());
                  for (const auto& column :
                       numeric_sequential_basis)
                    numeric_sequential_sources.push_back(
                        &column);
                  const auto numeric_sequential_candidate =
                      materialize_local_basis_weights(
                          numeric_sequential_sources,
                          acb_match->transformed_weights(),
                          result_checkpoint_identity +
                              ":continuity-audit:"
                              "numeric-((F*T)*P)*y");
                  reconstruction_routes[
                      "numeric-((F*T)*P)*y"] =
                      audit_route(
                          numeric_sequential_candidate);

                  const auto combined =
                      multiply_exact_laurent_matrices(
                          numeric_right, *preconditioner);
                  auto numeric_combined_basis =
                      right_transform_local_basis_exact<
                          ComplexBall>(
                          physical_matching_basis, combined,
                          result_checkpoint_identity +
                              ":continuity-audit:"
                              "numeric-F*(T*P)");
                  std::vector<
                      const LocalSolution<ComplexBall>*>
                      numeric_combined_sources;
                  numeric_combined_sources.reserve(
                      numeric_combined_basis.size());
                  for (const auto& column :
                       numeric_combined_basis)
                    numeric_combined_sources.push_back(
                        &column);
                  const auto numeric_combined_candidate =
                      materialize_local_basis_weights(
                          numeric_combined_sources,
                          acb_match->transformed_weights(),
                          result_checkpoint_identity +
                              ":continuity-audit:"
                              "numeric-(F*(T*P))*y");
                  reconstruction_routes[
                      "numeric-(F*(T*P))*y"] =
                      audit_route(numeric_combined_candidate);
                }
              }
            } catch (const std::exception& audit_error) {
              reconstruction_routes["audit_error"] =
                  audit_error.what();
            }
            json::object diagnostics{
                {"association", required_string(
                     native_match_summary,
                     "materialization_association")},
                {"normal_frame_status",
                 required_string(normal_frame, "status")},
                {"rational_shadow_basis",
                 exact_right_used_rational_shadow},
                {"epsilon", power},
                {"component", component},
                {"residual_upper",
                 continuity.residual_upper.approximate_upper()},
                {"allowed_upper",
                 continuity.allowed_upper.approximate_upper()},
                {"scale_lower",
                 continuity.scale_lower.approximate_upper()},
                {"overlaps", continuity.overlaps},
                {"relative_tolerance",
                 required_string(refinement_record,
                                 "relative_tolerance")},
                {"incoming", compact_ball(expected)},
                {"receiving", compact_ball(
                     receiving_evaluation.value.at(power,
                                                   component))},
                {"receiving_epsilon",
                 json::object{
                     {"min",
                      receiving_evaluation.value.epsilon.min_power},
                     {"max",
                      receiving_evaluation.value.epsilon.complete_max}}},
                {"incoming_epsilon",
                 json::object{
                     {"min",
                      incoming_evaluation.value.epsilon.min_power},
                     {"max",
                      incoming_evaluation.value.epsilon.complete_max}}},
                {"matching_taylor_complete_max",
                 matching_taylor_complete_max.has_value()
                     ? json::value(*matching_taylor_complete_max)
                     : json::value(nullptr)},
                {"incoming_matching_taylor_complete_max",
                 incoming_matching_taylor_complete_max.has_value()
                     ? json::value(
                           *incoming_matching_taylor_complete_max)
                     : json::value(nullptr)},
                {"reconstruction_routes",
                 std::move(reconstruction_routes)}};
            throw MatchingArithmeticError(
                MatchingArithmeticErrorCode::
                    MaterializedContinuityInconclusive,
                "materialized receiving local fails its certified handoff continuity at epsilon=" +
                std::to_string(power) + ", component=" +
                std::to_string(component) + "; diagnostics=" +
                json::serialize(diagnostics),
                component, std::nullopt, power);
          }
        }
      }
    }
    json::array weight_windows;
    weight_windows.reserve(weights.size());
    for (const auto& weight : weights)
      weight_windows.push_back(json::object{
          {"min", weight.min_power()},
          {"max", weight.complete_max()}});
    json::object output_record{
        {"checkpoint_identity", result_checkpoint_identity},
        {"chart", receiving_chart},
        {"source_operator_identity", receiving_operator},
        {"epsilon", json::object{
             {"min", solution.epsilon.min_power},
             {"max", solution.epsilon.complete_max}}},
        {"taylor_complete_max", solution.taylor_complete_max},
        {"dimension", solution.dimension}};
    const auto handoff_schema = required_string(handoff_, "schema");
    const bool compact_derivation =
        handoff_schema ==
            "diffexp3-retained-exact-plan-match-hop-v2" ||
        handoff_schema ==
            "diffexp3-retained-exact-plan-match-hop-v3";
    json::object derivation;
    if (compact_derivation) {
      if (!plan_owner_ || !planned_equation_owner ||
          (equation_owner &&
           equation_owner.get() != planned_equation_owner.get()) ||
          required_string(handoff_, "tile_plan") != plan_owner_->handle() ||
          required_string(handoff_, "tile_plan_checkpoint_identity") !=
              plan_owner_->checkpoint_identity() ||
          required_string(handoff_, "result_checkpoint_identity") !=
              checkpoint_identity_)
        throw std::logic_error(
            "compact plan-match derivation lost its live plan, match, or equation owner");
      const auto& producing = as_object(
          handoff_.at("producing"), "compact plan-match producing record");
      const auto& receiving = as_object(
          handoff_.at("receiving"), "compact plan-match receiving record");
      derivation = json::object{
          {"schema", "diffexp3-retained-plan-match-local-materialization-v2"},
          {"capability", kRetainedPlannedMatchMaterializationCapability},
          {"source_match", handle_},
          {"source_match_checkpoint_identity", checkpoint_identity_},
          {"tile_plan", handoff_.at("tile_plan")},
          {"tile_plan_checkpoint_identity",
           handoff_.at("tile_plan_checkpoint_identity")},
          {"arm", handoff_.at("arm")},
          {"match", handoff_.at("match")},
          {"incoming", producing.at("incoming")},
          {"basis", receiving.at("basis")},
          {"weight_windows", std::move(weight_windows)},
          {"match_certified_complete_max", match_certified_max},
          {"output", std::move(output_record)},
          {"equation_owner_signature_identity",
           planned_equation_owner_signature_binding(
               *planned_equation_owner)},
          {"equation_payload_identity",
           planned_equation_payload_binding(
               *planned_equation_owner)},
          {"scope", "single-match-receiving-local"},
          {"coefficient_transport", "native-retained-only"},
          {"whole_arm_complete", false}};
    } else {
      derivation = json::object{
          {"schema", "diffexp3-retained-plan-match-local-materialization-v1"},
          {"capability", kRetainedPlannedMatchMaterializationCapability},
          {"source_match", handle_},
          {"source_match_checkpoint_identity", checkpoint_identity_},
          {"source_match_provenance_identity",
           required_string(native_match_summary, "provenance_identity")},
          {"planned_hop_provenance_identity", provenance_identity_},
          {"planned_hop", handoff_},
          {"weight_windows", std::move(weight_windows)},
          {"match_certified_complete_max", match_certified_max},
          {"output", std::move(output_record)},
          {"scope", "single-match-receiving-local"},
          {"coefficient_transport", "native-retained-only"},
          {"whole_arm_complete", false}};
    }
    const auto derivation_identity = json::serialize(
        canonical_json_value(derivation));
    derivation["provenance_identity"] = derivation_identity;
    std::vector<const RegularTaylorTailModelResult*> basis_tail_models;
    basis_tail_models.reserve(typed_basis.size());
    for (const auto& column : typed_basis)
      basis_tail_models.push_back(&column->tail_model());
    std::vector<RegularTaylorTailModelResult>
        transformed_basis_tail_models;
    bool transformed_tail_association = false;
    if constexpr (std::is_same_v<Scalar, ComplexBall>) {
      if (exact_right_physical_association) {
        const RegularTaylorTailModel* reference_model = nullptr;
        bool can_rebind = true;
        for (const auto* candidate : basis_tail_models) {
          if (candidate == nullptr ||
              candidate->status != TailMajorantStatus::Certified ||
              !candidate->model.has_value()) {
            can_rebind = false;
            break;
          }
          if (reference_model == nullptr)
            reference_model = &*candidate->model;
          else if (!tail_majorant_detail::same_operator_payload(
                       *reference_model, *candidate->model)) {
            can_rebind = false;
            break;
          }
        }
        if (can_rebind && reference_model != nullptr) {
          transformed_basis_tail_models.reserve(
              materialization_basis.size());
          for (const auto* transformed : materialization_basis) {
            RegularTaylorTailModel model = *reference_model;
            model.epsilon = transformed->epsilon;
            model.taylor_complete_max =
                transformed->taylor_complete_max;
            model.chart = transformed->chart;
            model.prescriptions = transformed->prescriptions;
            model.local_checkpoint_identity =
                transformed->checkpoint_identity;
            model.initial_row_upper.assign(
                transformed->epsilon.width(), Magnitude::zero());
            if (transformed->sectors.size() != 1 ||
                transformed->sectors.front().log_power != 0 ||
                transformed->sectors.front().a.is_zero !=
                    TruthValue::Yes ||
                transformed->sectors.front().b.is_zero !=
                    TruthValue::Yes) {
              can_rebind = false;
              break;
            }
            const auto& sector = transformed->sectors.front();
            for (std::size_t epsilon_index = 0;
                 epsilon_index < transformed->epsilon.width();
                 ++epsilon_index)
              for (std::uint32_t component = 0;
                   component < transformed->dimension; ++component)
                model.initial_row_upper[epsilon_index] =
                    Magnitude::maximum(
                        model.initial_row_upper[epsilon_index],
                        Magnitude::upper_abs(
                            sector.coefficients[
                                local_detail::sector_index(
                                    *transformed, epsilon_index, 0,
                                    component)]));
            model.provenance =
                "regular homogeneous tail model rebound to the exact-right "
                "constant match frame; the receiving basis transformation "
                "is t-independent and every source column owns the same "
                "epsilon-decoupled q/N payload; " +
                reference_model->provenance;
            transformed_basis_tail_models.push_back(
                {TailMajorantStatus::Certified, std::move(model),
                 "certified regular tail model rebound through the retained exact-right match frame"});
          }
          if (can_rebind) {
            basis_tail_models.clear();
            basis_tail_models.reserve(
                transformed_basis_tail_models.size());
            for (const auto& rebound :
                 transformed_basis_tail_models)
              basis_tail_models.push_back(&rebound);
            basis_solutions = materialization_basis;
            transformed_tail_association = true;
          }
        }
      }
    }
    auto tail_model = derive_materialized_regular_homogeneous_tail_model(
        basis_solutions, basis_tail_models,
        transformed_tail_association
            ? *materialization_weights : weights,
        solution,
        receiving_operator, checkpoint_identity_);
    std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>>
        physical_equation;
    if (equation_owner != nullptr) {
      if (equation_owner->equation_owner_handle() != receiving_chart ||
          equation_owner->equation_operator_identity() !=
              receiving_operator ||
          std::string(equation_owner->equation_scalar_domain()) !=
              (std::is_same_v<Scalar, Rational> ? "rational" : "acb"))
        throw std::logic_error(
            "retained plan-match basis equation owner disagrees with its receiving chart");
      const auto erased = equation_owner->physical_ode_erased();
      if (!erased)
        throw std::logic_error(
            "retained plan-match basis equation owner lost its physical q/C payload");
      physical_equation =
          std::static_pointer_cast<const PreparedPhysicalClearedODE<Scalar>>(
              erased);
      physical_ode_detail::validate_ode(*physical_equation);
      if (physical_equation->exact_payload_record !=
              equation_owner->physical_payload_record() ||
          physical_equation->payload_identity !=
              equation_owner->physical_payload_identity() ||
          physical_equation->owner_signature_identity !=
              equation_owner->owner_signature_identity())
        throw std::logic_error(
            "retained plan-match physical q/C payload differs from its shared basis owner");
    }
    const auto elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    NativeLocalDiagnostics diagnostics;
    diagnostics.top_valid = materialized_top;
    diagnostics.kernel_ms = elapsed_ms;
    return make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
        local_handle, receiving_chart, receiving_operator,
        std::move(solution), precision_bits,
        std::vector<PseudoHit<Scalar>>{}, diagnostics, std::nullopt,
        std::move(derivation), std::static_pointer_cast<void>(self),
        std::move(tail_model), std::nullopt, true, true,
        std::move(equation_owner), std::move(physical_equation));
  }

  std::shared_ptr<StoredMatchBase> match_;
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  json::object handoff_;
  double elapsed_ms_ = 0.0;
  std::shared_ptr<StoredTilePlan> plan_owner_;
  std::vector<std::shared_ptr<StoredLocalBase>> basis_owners_;
  std::shared_ptr<StoredLocalBase> incoming_owner_;
  std::atomic<std::uint64_t> materializations_{0};
  mutable std::mutex terminal_factorized_basis_cache_mutex_;
  mutable slong terminal_factorized_basis_cache_precision_bits_ = 0;
  mutable std::shared_ptr<
      const std::vector<LocalSolution<ComplexBall>>>
      terminal_factorized_basis_cache_;
  mutable detail::ImmutableRecursiveCache<
      std::pair<std::uint32_t, std::int32_t>,
      adjoint_observable_detail::NormalizedBackwardAdjointExactODE>
      terminal_normalized_adjoint_cache_;
  mutable adjoint_observable_detail::BackwardAdjointRealRayOperatorCache
      terminal_real_ray_operator_cache_;
  mutable std::atomic<std::uint32_t>
      terminal_composed_taylor_order_floor_{0};
};

json::object transport_local_consumer_epsilon_contract(
    const std::shared_ptr<StoredLocalBase>& source,
    std::shared_ptr<StoredPlannedMatchHop> terminal_match = nullptr) {
  if (!source)
    throw std::invalid_argument(
        "transport consumer epsilon contract requires one retained local");
  const auto source_summary = source->summary();
  const auto source_window = json::object{
      {"min", source_summary.at("epsilon_min")},
      {"max", source_summary.at("epsilon_max")}};
  if (!terminal_match) {
    if (const auto erased = source->terminal_factorized_owner();
        erased != nullptr)
      terminal_match =
          std::static_pointer_cast<StoredPlannedMatchHop>(erased);
  }

  json::array consumer_sources;
  std::int32_t output_min_shift = 0;
  auto basis_min_power = as_i32(
      source_summary.at("epsilon_min"),
      "transport tile source epsilon minimum");
  std::int32_t projection_weight_min_power = 0;
  if (terminal_match && terminal_match->has_terminal_acb_factorization()) {
    const auto physical_basis =
        terminal_match->terminal_acb_basis_owners();
    consumer_sources.reserve(physical_basis.size());
    std::optional<std::int32_t> physical_min_power;
    for (const auto& physical_source : physical_basis) {
      const auto physical_summary = physical_source->summary();
      const auto minimum = as_i32(
          physical_summary.at("epsilon_min"),
          "terminal physical source epsilon minimum");
      if (!physical_min_power.has_value() ||
          minimum < *physical_min_power)
        physical_min_power = minimum;
      consumer_sources.push_back(json::object{
          {"min", physical_summary.at("epsilon_min")},
          {"max", physical_summary.at("epsilon_max")}});
    }
    if (!physical_min_power.has_value())
      throw std::logic_error(
          "terminal factorized consumer has no physical source minimum");
    output_min_shift =
        terminal_match->terminal_acb_factorized_input_min_power();
    basis_min_power = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(*physical_min_power) +
            output_min_shift,
        "terminal factorized basis conservative epsilon minimum");
    projection_weight_min_power =
        terminal_match->terminal_acb_transformed_weight_min_power();
  } else {
    consumer_sources.push_back(source_window);
  }
  return json::object{
      {"sources", std::move(consumer_sources)},
      {"output_min_shift", output_min_shift},
      {"basis_min_power", basis_min_power},
      {"projection_weight_min_power", projection_weight_min_power}};
}

class StoredTransportArmState {
 public:
  StoredTransportArmState(
      std::string handle, std::string checkpoint_identity,
      std::string arm, std::shared_ptr<StoredTilePlan> plan_owner,
      std::shared_ptr<StoredLocalBase> anchor_owner,
      std::vector<std::vector<std::shared_ptr<StoredLocalBase>>> basis_owners,
      std::vector<std::shared_ptr<StoredPlannedMatchHop>> matches,
      std::vector<std::shared_ptr<StoredLocalBase>> tile_sources,
      EpsilonWindow work_epsilon,
      std::int32_t public_required_complete_max,
      std::int32_t match_required_complete_max,
      json::object refinement, double elapsed_ms,
      bool compact_provenance = true)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        arm_(std::move(arm)), plan_owner_(std::move(plan_owner)),
        anchor_owner_(std::move(anchor_owner)),
        basis_owners_(std::move(basis_owners)),
        matches_(std::move(matches)),
        tile_sources_(std::move(tile_sources)),
        work_epsilon_(work_epsilon),
        public_required_complete_max_(public_required_complete_max),
        match_required_complete_max_(match_required_complete_max),
        refinement_(std::move(refinement)), elapsed_ms_(elapsed_ms),
        compact_provenance_(compact_provenance) {
    validate();
    provenance_identity_ = json::serialize(
        canonical_json_value(provenance_record()));
  }

  // Ownership-taking streamed publication already has the complete ordered
  // chain of materialized tile locals.  Each non-anchor local carries a
  // sealed v2 plan-match derivation that was validated before its (large)
  // basis and native match owners were released.  Rebind that internal
  // certificate directly instead of round-tripping basis/match provenance
  // through Wolfram and copying it into the transport state.
  StoredTransportArmState(
      std::string handle, std::string checkpoint_identity,
      std::string arm, std::shared_ptr<StoredTilePlan> plan_owner,
      std::shared_ptr<StoredLocalBase> anchor_owner,
      std::vector<std::shared_ptr<StoredLocalBase>> tile_sources,
      EpsilonWindow work_epsilon,
      std::int32_t public_required_complete_max,
      std::int32_t match_required_complete_max,
      json::object refinement, double elapsed_ms,
      std::shared_ptr<StoredPlannedMatchHop>
          terminal_factorized_match = nullptr)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        arm_(std::move(arm)), plan_owner_(std::move(plan_owner)),
        anchor_owner_(std::move(anchor_owner)),
        tile_sources_(std::move(tile_sources)),
        terminal_factorized_match_(
            std::move(terminal_factorized_match)),
        work_epsilon_(work_epsilon),
        public_required_complete_max_(public_required_complete_max),
        match_required_complete_max_(match_required_complete_max),
        refinement_(std::move(refinement)), elapsed_ms_(elapsed_ms),
        compact_provenance_(true), consumed_compact_(true),
        consumed_certificate_only_(true) {
    validate();
    provenance_identity_ = json::serialize(
        canonical_json_value(provenance_record()));
  }

  StoredTransportArmState(
      std::string handle, std::string checkpoint_identity,
      std::string arm, std::shared_ptr<StoredTilePlan> plan_owner,
      std::shared_ptr<StoredLocalBase> anchor_owner,
      json::array basis_references, json::array match_references,
      std::vector<std::shared_ptr<StoredLocalBase>> tile_sources,
      EpsilonWindow work_epsilon,
      std::int32_t public_required_complete_max,
      std::int32_t match_required_complete_max,
      json::object refinement, double elapsed_ms)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        arm_(std::move(arm)), plan_owner_(std::move(plan_owner)),
        anchor_owner_(std::move(anchor_owner)),
        tile_sources_(std::move(tile_sources)),
        work_epsilon_(work_epsilon),
        public_required_complete_max_(public_required_complete_max),
        match_required_complete_max_(match_required_complete_max),
        refinement_(std::move(refinement)), elapsed_ms_(elapsed_ms),
        compact_provenance_(true), consumed_compact_(true),
        cached_basis_references_(std::move(basis_references)),
        cached_match_references_(std::move(match_references)) {
    validate();
    provenance_identity_ = json::serialize(
        canonical_json_value(provenance_record()));
  }

  const std::string& handle() const { return handle_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& provenance_identity() const {
    return provenance_identity_;
  }
  const std::string& arm_name() const { return arm_; }
  double elapsed_ms() const { return elapsed_ms_; }
  const std::shared_ptr<StoredTilePlan>& plan_owner() const {
    return plan_owner_;
  }
  const std::shared_ptr<StoredLocalBase>& anchor_owner() const {
    return anchor_owner_;
  }
  const std::vector<std::vector<std::shared_ptr<StoredLocalBase>>>&
  basis_owners() const {
    return basis_owners_;
  }
  const std::vector<std::shared_ptr<StoredPlannedMatchHop>>& matches() const {
    return matches_;
  }
  const std::vector<std::shared_ptr<StoredLocalBase>>& tile_sources() const {
    return tile_sources_;
  }
  const std::shared_ptr<StoredLocalBase>& final_local() const {
    return tile_sources_.back();
  }
  const std::shared_ptr<StoredPlannedMatchHop>&
  terminal_factorized_match() const {
    return terminal_factorized_match_;
  }
  std::int32_t public_required_complete_max() const {
    return public_required_complete_max_;
  }

  void require_contraction_counter_capacity(
      std::size_t observables) const {
    const auto operation_count = contraction_operations_.load();
    const auto observable_count = contracted_observables_.load();
    if (operation_count == std::numeric_limits<std::uint64_t>::max() ||
        observables > std::numeric_limits<std::uint64_t>::max() -
                          observable_count)
      throw std::overflow_error(
          "retained transport-arm contraction counter overflow");
  }

  void note_contraction_success(std::size_t observables) noexcept {
    contraction_operations_.fetch_add(1);
    contracted_observables_.fetch_add(
        static_cast<std::uint64_t>(observables));
  }

  void require_endpoint_batch_counter_capacity(std::size_t rows) const {
    const auto operation_count = endpoint_batch_operations_.load();
    const auto row_count = endpoint_rows_.load();
    if (operation_count == std::numeric_limits<std::uint64_t>::max() ||
        rows > std::numeric_limits<std::uint64_t>::max() - row_count)
      throw std::overflow_error(
          "retained transport endpoint-batch counter overflow");
  }

  void note_endpoint_batch_success(std::size_t rows) noexcept {
    endpoint_batch_operations_.fetch_add(1);
    endpoint_rows_.fetch_add(static_cast<std::uint64_t>(rows));
  }

  json::object summary() const {
    json::array tile_source_epsilon;
    json::array tile_consumer_epsilon;
    tile_source_epsilon.reserve(tile_sources_.size());
    tile_consumer_epsilon.reserve(tile_sources_.size());
    for (std::size_t tile = 0; tile < tile_sources_.size(); ++tile) {
      const auto& source = tile_sources_[tile];
      if (!source)
        throw std::logic_error(
            "retained transport state owns a null tile source");
      const auto source_summary = source->summary();
      const auto source_window = json::object{
          {"min", source_summary.at("epsilon_min")},
          {"max", source_summary.at("epsilon_max")}};
      tile_source_epsilon.push_back(source_window);
      tile_consumer_epsilon.push_back(
          transport_local_consumer_epsilon_contract(
              source,
              tile + 1 == tile_sources_.size()
                  ? terminal_factorized_match_
                  : nullptr));
    }
    return json::object{
        {"transport_state", handle_},
        {"capability", kRetainedTransportArmStateCapability},
        {"native_retained", true},
        {"json_coefficients", 0},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"tile_plan", plan_owner_->handle()},
        {"tile_plan_checkpoint_identity",
         plan_owner_->checkpoint_identity()},
        {"arm", arm_},
        {"matches", transition_count()},
        {"value_hops", value_hop_count()},
        {"basis_matches", basis_match_count()},
        {"tiles", tile_sources_.size()},
        {"tile_source_epsilon", std::move(tile_source_epsilon)},
        {"tile_consumer_epsilon", std::move(tile_consumer_epsilon)},
        {"contraction_operations", contraction_operations_.load()},
        {"contracted_observables", contracted_observables_.load()},
        {"endpoint_batch_operations", endpoint_batch_operations_.load()},
        {"endpoint_rows", endpoint_rows_.load()},
        {"terminal_factorized_match",
         terminal_factorized_match_ != nullptr},
        {"epsilon", epsilon_record()},
        {"refinement", refinement_},
        {"final_local", local_reference(final_local())},
        {"strong_ownership", json::object{
             {"tile_plan", true}, {"anchor", true},
             {"basis_locals", consumed_compact_ ? 0 : basis_owner_count()},
             {"matches", consumed_compact_ ? 0 : matches_.size()},
             {"terminal_factorized_match",
              terminal_factorized_match_ != nullptr},
             {"tile_sources", tile_sources_.size()}}},
        {"elapsed_ms", elapsed_ms_}};
  }

  json::object stats_json() const {
    auto result = summary();
    const auto* terminal_diagnostics =
        std::getenv("DIFFEXP_DIAGNOSTIC_TERMINAL_STATE");
    if (terminal_factorized_match_ != nullptr &&
        terminal_diagnostics != nullptr &&
        std::string(terminal_diagnostics) == "1")
      result["terminal_diagnostic"] =
          terminal_factorized_match_
              ->terminal_mode_diagnostic_summary();
    else if (!matches_.empty() && terminal_diagnostics != nullptr &&
             std::string(terminal_diagnostics) == "1" &&
             matches_.back()->has_acb_match())
      result["final_match_diagnostic"] =
          matches_.back()->compact_terminal_match_diagnostic();
    result["stats_queries"] = stats_queries_.fetch_add(1) + 1;
    return result;
  }

  json::object checkpoint_record() const {
    return json::object{
        {"schema", consumed_certificate_only_
             ? terminal_factorized_match_ != nullptr
                 ? "diffexp3-retained-transport-arm-state-v6"
                 : "diffexp3-retained-transport-arm-state-v5"
             : consumed_compact_
             ? "diffexp3-retained-transport-arm-state-v4"
             : compact_provenance_
             ? "diffexp3-retained-transport-arm-state-v3"
             : "diffexp3-retained-transport-arm-state-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"provenance", provenance_record()},
        {"elapsed_ms", elapsed_ms_},
        {"runtime_stats",
         json::object{{"stats_queries", stats_queries_.load()},
                      {"contraction_operations",
                       contraction_operations_.load()},
                      {"contracted_observables",
                       contracted_observables_.load()},
                      {"endpoint_batch_operations",
                       endpoint_batch_operations_.load()},
                      {"endpoint_rows", endpoint_rows_.load()}}}};
  }

  void restore_runtime_stats(std::uint64_t stats_queries,
                             std::uint64_t contraction_operations,
                             std::uint64_t contracted_observables,
                             std::uint64_t endpoint_batch_operations,
                             std::uint64_t endpoint_rows) {
    stats_queries_.store(stats_queries);
    contraction_operations_.store(contraction_operations);
    contracted_observables_.store(contracted_observables);
    endpoint_batch_operations_.store(endpoint_batch_operations);
    endpoint_rows_.store(endpoint_rows);
  }

 private:
  static json::object local_reference(
      const std::shared_ptr<StoredLocalBase>& local) {
    if (!local)
      throw std::logic_error(
          "transport-arm state contains a null local owner");
    return json::object{
        {"local", local->handle()}, {"chart", local->source_chart()},
        {"source_operator_identity", local->source_operator_identity()},
        {"checkpoint_identity", local->checkpoint_identity()},
        {"coefficient_domain", local->scalar_domain()}};
  }

  json::array basis_reference() const {
    if (consumed_compact_) return cached_basis_references_;
    json::array result;
    result.reserve(basis_owners_.size());
    for (const auto& basis : basis_owners_) {
      json::array columns;
      columns.reserve(basis.size());
      for (const auto& column : basis)
        columns.push_back(local_reference(column));
      result.push_back(std::move(columns));
    }
    return result;
  }

  json::array match_reference() const {
    if (consumed_compact_) return cached_match_references_;
    json::array result;
    result.reserve(matches_.size());
    for (std::size_t index = 0; index < matches_.size(); ++index) {
      auto reference = json::object{
          {"index", index}, {"match", matches_[index]->handle()},
          {"checkpoint_identity", matches_[index]->checkpoint_identity()}};
      if (!compact_provenance_)
        reference["provenance_identity"] =
            matches_[index]->provenance_identity();
      result.push_back(std::move(reference));
    }
    return result;
  }

  json::array tile_source_reference() const {
    json::array result;
    result.reserve(tile_sources_.size());
    for (std::size_t index = 0; index < tile_sources_.size(); ++index) {
      auto record = local_reference(tile_sources_[index]);
      record["tile"] = index;
      result.push_back(std::move(record));
    }
    return result;
  }

  static json::object compact_tile_checkpoint_reference(
      const std::shared_ptr<StoredLocalBase>& local, std::size_t tile) {
    if (!local)
      throw std::logic_error(
          "certificate-only transport state contains a null tile local");
    return json::object{
        {"tile", tile}, {"local", local->handle()},
        {"chart", local->source_chart()},
        {"checkpoint_identity", local->checkpoint_identity()},
        {"coefficient_domain", local->scalar_domain()}};
  }

  json::array consumed_certificate_chain() const {
    json::array chain;
    chain.reserve(tile_sources_.size());
    for (std::size_t tile = 0; tile < tile_sources_.size(); ++tile) {
      auto record = compact_tile_checkpoint_reference(
          tile_sources_[tile], tile);
      if (tile == 0) {
        record["derivation"] = nullptr;
      } else {
        const auto& derivation = *tile_sources_[tile]->retained_derivation();
        const auto schema = required_string(derivation, "schema");
        if (is_retained_plan_value_handoff_schema(schema)) {
          record["derivation"] = json::object{
              {"schema",
               "diffexp3-consumed-plan-value-handoff-certificate-v1"},
              {"match", tile - 1},
              {"handoff_provenance_identity",
               derivation.at("provenance_identity")},
              {"incoming_checkpoint_identity",
               as_object(derivation.at("incoming"),
                         "certificate-only value incoming derivation")
                   .at("checkpoint_identity")},
              {"output_checkpoint_identity",
               as_object(derivation.at("output"),
                         "certificate-only value output derivation")
                   .at("checkpoint_identity")}};
        } else {
          record["derivation"] = json::object{
              {"schema", "diffexp3-consumed-plan-match-certificate-v1"},
              {"match", tile - 1},
              {"source_match_checkpoint_identity",
               derivation.at("source_match_checkpoint_identity")},
              {"incoming_checkpoint_identity",
               as_object(derivation.at("incoming"),
                         "certificate-only incoming derivation")
                   .at("checkpoint_identity")},
              {"output_checkpoint_identity",
               as_object(derivation.at("output"),
                         "certificate-only output derivation")
                   .at("checkpoint_identity")}};
        }
      }
      chain.push_back(std::move(record));
    }
    return chain;
  }

  json::object epsilon_record() const {
    return json::object{
        {"min", work_epsilon_.min_power},
        {"max", work_epsilon_.complete_max},
        {"required_complete_max", public_required_complete_max_},
        {"match_required_complete_max", match_required_complete_max_}};
  }

  json::object provenance_record() const {
    auto plan_reference = json::object{
        {"handle", plan_owner_->handle()},
        {"checkpoint_identity", plan_owner_->checkpoint_identity()}};
    if (consumed_certificate_only_) {
      json::object result{
          {"schema", terminal_factorized_match_ != nullptr
              ? "diffexp3-retained-native-transport-arm-state-v5"
              : "diffexp3-retained-native-transport-arm-state-v4"},
          {"checkpoint_identity", checkpoint_identity_},
          {"tile_plan", std::move(plan_reference)},
          {"arm", arm_},
          {"tile_checkpoint_chain", consumed_certificate_chain()},
          {"epsilon", epsilon_record()},
          {"refinement", refinement_}};
      if (terminal_factorized_match_ != nullptr)
        result["terminal_match"] = json::object{
            {"match", terminal_factorized_match_->handle()},
            {"checkpoint_identity",
             terminal_factorized_match_->checkpoint_identity()},
            {"provenance_identity",
             terminal_factorized_match_->provenance_identity()}};
      return result;
    }
    if (!compact_provenance_)
      plan_reference["provenance_identity"] =
          plan_owner_->provenance_identity();
    return json::object{
        {"schema", compact_provenance_
             ? (consumed_compact_
                    ? "diffexp3-retained-native-transport-arm-state-v3"
                    : "diffexp3-retained-native-transport-arm-state-v2")
             : "diffexp3-retained-native-transport-arm-state-v1"},
        {"checkpoint_identity", checkpoint_identity_},
        {"tile_plan", std::move(plan_reference)},
        {"arm", arm_},
        {"anchor", local_reference(anchor_owner_)},
        {"receiving_basis", basis_reference()},
        {"matches", match_reference()},
        {"tile_sources", tile_source_reference()},
        {"final_local", local_reference(final_local())},
        {"epsilon", epsilon_record()},
        {"refinement", refinement_}};
  }

  std::size_t basis_owner_count() const {
    std::size_t result = 0;
    for (const auto& basis : basis_owners_) {
      if (basis.size() > std::numeric_limits<std::size_t>::max() - result)
        throw std::overflow_error(
            "transport-arm basis owner count overflow");
      result += basis.size();
    }
    return result;
  }

  std::size_t transition_count() const {
    return consumed_certificate_only_
        ? tile_sources_.size() - 1
        : consumed_compact_ ? cached_match_references_.size()
                            : matches_.size();
  }

  std::size_t value_hop_count() const {
    if (!consumed_certificate_only_) return 0;
    std::size_t count = 0;
    for (std::size_t tile = 1; tile < tile_sources_.size(); ++tile) {
      const auto& derivation = tile_sources_[tile]->retained_derivation();
      if (derivation.has_value() &&
          is_retained_plan_value_handoff_schema(
              required_string(*derivation, "schema")))
        ++count;
    }
    return count;
  }

  std::size_t basis_match_count() const {
    return transition_count() - value_hop_count();
  }

  void validate() const {
    if (handle_.empty() || checkpoint_identity_.empty() ||
        !plan_owner_ || !anchor_owner_ || elapsed_ms_ < 0.0 ||
        !std::isfinite(elapsed_ms_))
      throw std::invalid_argument(
          "retained transport-arm state lost an identity or strong owner");
    const auto& retained = plan_owner_->arm(arm_);
    if ((!consumed_compact_ &&
         (basis_owners_.size() != retained.exact.matches.size() ||
          matches_.size() != retained.exact.matches.size())) ||
        (consumed_compact_ && !consumed_certificate_only_ &&
         (cached_basis_references_.size() != retained.exact.matches.size() ||
          cached_match_references_.size() != retained.exact.matches.size())) ||
        tile_sources_.size() != retained.exact.tiles.size() ||
        tile_sources_.size() != retained.exact.matches.size() + 1)
      throw std::invalid_argument(
          "retained transport-arm state does not reproduce its plan topology");
    if (tile_sources_.empty() || tile_sources_.front().get() !=
                                     anchor_owner_.get())
      throw std::invalid_argument(
          "retained transport-arm state lost its anchor tile source");
    (void)work_epsilon_.width();
    if (public_required_complete_max_ < work_epsilon_.min_power ||
        match_required_complete_max_ < public_required_complete_max_ ||
        match_required_complete_max_ > work_epsilon_.complete_max)
      throw std::invalid_argument(
          "retained transport-arm epsilon contract is inconsistent");
    require_exact_keys(refinement_, {"relative_tolerance", "max_steps"},
                       "retained transport-arm refinement policy");
    if (required_string(refinement_, "relative_tolerance").empty() ||
        as_u32(refinement_.at("max_steps"),
               "retained transport-arm refinement steps") > 32)
      throw std::invalid_argument(
          "retained transport-arm refinement policy is invalid");

    for (std::size_t tile = 0; tile < tile_sources_.size(); ++tile) {
      const auto& source = tile_sources_[tile];
      if (!source)
        throw std::invalid_argument(
            "retained transport-arm state contains a null tile source");
      const auto source_summary = source->summary();
      const auto source_epsilon_max = as_i32(
          source_summary.at("epsilon_max"),
          "retained transport-arm tile-source epsilon maximum");
      const auto source_top_valid = parse_validity(
          source_summary.at("top_valid"));
      const auto source_complete_max =
          source_top_valid == kCompleteInfinity
              ? source_epsilon_max
              : std::min(source_epsilon_max, source_top_valid);
      if (source_complete_max < public_required_complete_max_)
        throw std::invalid_argument(
            "retained transport-arm tile source does not cover its public "
            "required complete epsilon maximum");
      const auto& exact_tile = retained.exact.tiles[tile];
      const auto& chart = retained.charts.at(exact_tile.chart);
      if (source->source_chart() != chart.handle)
        throw std::invalid_argument(
            "retained transport-arm tile source belongs to a different chart");
      source->require_exact_plan_binding(
          chart.local_geometry, chart.prescriptions,
          "retained transport-arm tile source");
    }
    if (consumed_certificate_only_) {
      validate_consumed_certificate_chain(retained);
      validate_terminal_factorized_match(retained);
      return;
    }
    if (consumed_compact_) {
      validate_consumed_references(retained);
      return;
    }
    for (std::size_t index = 0; index < matches_.size(); ++index) {
      const auto& match = matches_[index];
      const auto& basis = basis_owners_[index];
      if (!match || basis.empty() ||
          std::any_of(basis.begin(), basis.end(),
                      [](const auto& owner) { return owner == nullptr; }) ||
          match->plan_owner().get() != plan_owner_.get() ||
          match->incoming_owner().get() != tile_sources_[index].get() ||
          match->basis_owners().size() != basis.size())
        throw std::invalid_argument(
            "retained transport-arm match lost its exact owner set");
      for (std::size_t column = 0; column < basis.size(); ++column)
        if (match->basis_owners()[column].get() != basis[column].get())
          throw std::invalid_argument(
              "retained transport-arm basis differs from its match owner");
      const auto& handoff = match->handoff();
      if (required_string(handoff, "arm") != arm_ ||
          as_u64(handoff.at("match"),
                 "retained transport-arm match index") != index)
        throw std::invalid_argument(
            "retained transport-arm match provenance is out of order");
      const auto& next = tile_sources_[index + 1];
      if (!next->retained_derivation().has_value() ||
          next->retained_derivation_owner().get() != match.get() ||
          required_string(*next->retained_derivation(), "source_match") !=
              match->handle())
        throw std::invalid_argument(
            "retained transport-arm tile source is not materialized from its match");
    }
  }

  void validate_consumed_certificate_chain(
      const RetainedArmPlan& retained) const {
    if (tile_sources_.size() != retained.exact.matches.size() + 1)
      throw std::invalid_argument(
          "certificate-only transport state has the wrong tile chain length");
    for (std::size_t tile = 1; tile < tile_sources_.size(); ++tile) {
      const auto& incoming = tile_sources_[tile - 1];
      const auto& output = tile_sources_[tile];
      if (!output || !output->has_sealed_plan_match_lineage() ||
          output->retained_derivation_owner() != nullptr ||
          !output->retained_derivation().has_value())
        throw std::invalid_argument(
            "certificate-only transport tile has no sealed plan-match derivation");
      const auto& derivation = *output->retained_derivation();
      const auto schema = required_string(derivation, "schema");
      const bool value_handoff =
          is_retained_plan_value_handoff_schema(schema);
      if ((!value_handoff && schema !=
              "diffexp3-retained-plan-match-local-materialization-v2") ||
          required_string(derivation, "tile_plan") !=
              plan_owner_->handle() ||
          required_string(derivation, "tile_plan_checkpoint_identity") !=
              plan_owner_->checkpoint_identity() ||
          required_string(derivation, "arm") != arm_ ||
          as_u64(derivation.at("match"),
                 "certificate-only plan-match index") != tile - 1)
        throw std::invalid_argument(
            "certificate-only transport derivation differs from its retained plan");
      const auto& source = as_object(
          derivation.at("incoming"),
          "certificate-only transport incoming derivation");
      const auto& target = as_object(
          derivation.at("output"),
          "certificate-only transport output derivation");
      const bool bounded_source =
          source.if_contains("source_operator_reference") != nullptr;
      if (required_string(source, "local") != incoming->handle() ||
          required_string(source, "checkpoint_identity") !=
              incoming->checkpoint_identity() ||
          (bounded_source
               ? !compact_matching_identity_reference_matches(
                     source.at("source_operator_reference"),
                     incoming->source_operator_identity(),
                     "certificate-only incoming source-operator reference")
               : required_string(source, "source_operator_identity") !=
                     incoming->source_operator_identity()) ||
          required_string(target, "checkpoint_identity") !=
              output->checkpoint_identity() ||
          required_string(target, "chart") != output->source_chart() ||
          required_string(target, "source_operator_identity") !=
              output->source_operator_identity())
        throw std::invalid_argument(
            "certificate-only transport derivation breaks its ordered tile chain");
      const auto equation_owner = output->retained_equation_owner();
      if (!equation_owner ||
          required_string(derivation,
                          "equation_owner_signature_identity") !=
              equation_owner->owner_signature_identity() ||
          required_string(derivation, "equation_payload_identity") !=
              equation_owner->physical_payload_identity())
        throw std::invalid_argument(
            "certificate-only transport derivation lost its match/equation certificate");
      if (!value_handoff &&
          required_string(derivation,
                          "source_match_checkpoint_identity").empty())
        throw std::invalid_argument(
            "certificate-only plan match lost its source-match checkpoint");
    }
  }

  void validate_terminal_factorized_match(
      const RetainedArmPlan& retained) const {
    if (!terminal_factorized_match_) return;
    if (retained.exact.matches.empty() || tile_sources_.size() < 2)
      throw std::invalid_argument(
          "terminal factorized match has no terminal arm transition");
    const auto terminal_index = retained.exact.matches.size() - 1;
    const auto& match = terminal_factorized_match_;
    const auto& handoff = match->handoff();
    if (match->plan_owner().get() != plan_owner_.get() ||
        match->incoming_owner().get() !=
            tile_sources_[terminal_index].get() ||
        required_string(handoff, "arm") != arm_ ||
        as_u64(handoff.at("match"),
               "terminal factorized match index") != terminal_index)
      throw std::invalid_argument(
          "terminal factorized match differs from its arm, plan, or incoming tile");
    const auto& output = tile_sources_.back();
    if (!output->has_sealed_plan_match_lineage() ||
        !output->retained_derivation().has_value() ||
        required_string(*output->retained_derivation(),
                        "source_match") != match->handle() ||
        required_string(*output->retained_derivation(),
                        "source_match_checkpoint_identity") !=
            match->checkpoint_identity())
      throw std::invalid_argument(
          "terminal factorized match differs from its sealed output lineage");
    const auto attached = output->terminal_factorized_owner();
    if (attached != nullptr &&
        (attached.get() != match.get() ||
         output->terminal_factorized_owner_checkpoint_identity() !=
             match->checkpoint_identity()))
      throw std::invalid_argument(
          "terminal factorized output-local attachment changed before state publication");
  }

  void validate_consumed_references(const RetainedArmPlan& retained) const {
    for (std::size_t index = 0; index < cached_basis_references_.size();
         ++index) {
      const auto& raw_basis = as_array(
          cached_basis_references_[index],
          "consumed transport-arm cached basis");
      if (raw_basis.empty())
        throw std::invalid_argument(
            "consumed transport-arm cached basis cannot be empty");
      for (const auto& raw_column : raw_basis) {
        const auto& column = as_object(
            raw_column, "consumed transport-arm cached basis column");
        const bool bounded_operator =
            column.if_contains("source_operator_reference") != nullptr;
        if (bounded_operator)
          require_exact_keys(
              column,
              {"local", "chart", "source_operator_reference",
               "checkpoint_identity", "coefficient_domain"},
              "consumed bounded transport-arm cached basis column");
        else
          require_exact_keys(
              column,
              {"local", "chart", "source_operator_identity",
               "checkpoint_identity", "coefficient_domain"},
              "consumed legacy transport-arm cached basis column");
        (void)scoped_handle_id(required_string(column, "local"), "l:",
                               "consumed basis local");
        if (required_string(column, "chart") !=
                retained.charts.at(
                    retained.exact.matches[index].receiving_chart)
                    .handle ||
            (bounded_operator
                 ? !compact_matching_identity_reference_matches(
                       column.at("source_operator_reference"),
                       tile_sources_[index + 1]
                           ->source_operator_identity(),
                       "consumed basis source-operator reference")
                 : required_string(
                       column, "source_operator_identity") !=
                       tile_sources_[index + 1]
                           ->source_operator_identity()) ||
            required_string(column, "checkpoint_identity").empty() ||
            required_string(column, "coefficient_domain") !=
                tile_sources_[index + 1]->scalar_domain())
          throw std::invalid_argument(
              "consumed transport-arm cached basis identity is inconsistent");
      }

      const auto& match = as_object(
          cached_match_references_[index],
          "consumed transport-arm cached match");
      require_exact_keys(
          match,
          {"index", "checkpoint_identity", "provenance_identity",
           "planned_hop", "sealed_local_lineage"},
          "consumed transport-arm cached match");
      if (as_u64(match.at("index"),
                 "consumed transport-arm match index") != index)
        throw std::invalid_argument(
            "consumed transport-arm cached matches are out of order");
      const auto& handoff = as_object(
          match.at("planned_hop"),
          "consumed transport-arm planned hop");
      const auto handoff_schema = required_string(handoff, "schema");
      const bool compact_handoff =
          handoff_schema ==
              "diffexp3-retained-exact-plan-match-hop-v2" ||
          handoff_schema ==
              "diffexp3-retained-exact-plan-match-hop-v3";
      const bool bounded_handoff = handoff_schema ==
          "diffexp3-retained-exact-plan-match-hop-v3";
      if (compact_handoff)
        require_exact_keys(
            handoff,
            {"schema", "tile_plan", "tile_plan_checkpoint_identity",
             "arm", "match", "geometry", "producing", "receiving",
             "result_checkpoint_identity", "advance"},
            "consumed compact transport-arm planned hop");
      else
        require_exact_keys(
            handoff,
            {"schema", "tile_plan", "tile_plan_checkpoint_identity",
             "tile_plan_provenance_identity", "arm", "match", "geometry",
             "producing", "receiving", "result_checkpoint_identity",
             "native_match_provenance_identity", "advance"},
            "consumed transport-arm planned hop");
      if ((!compact_handoff && handoff_schema !=
              "diffexp3-retained-exact-plan-match-hop-v1") ||
          required_string(handoff, "tile_plan") != plan_owner_->handle() ||
          required_string(handoff, "tile_plan_checkpoint_identity") !=
              plan_owner_->checkpoint_identity() ||
          (!compact_handoff &&
           required_string(handoff, "tile_plan_provenance_identity") !=
              plan_owner_->provenance_identity()) ||
          required_string(handoff, "arm") != arm_ ||
          as_u64(handoff.at("match"),
                 "consumed planned-hop match index") != index ||
          handoff.at("geometry") != encode_plan_match(retained, index) ||
          required_string(handoff, "result_checkpoint_identity") !=
              required_string(match, "checkpoint_identity") ||
          (compact_handoff &&
           json::serialize(canonical_json_value(handoff)) !=
              required_string(match, "provenance_identity")) ||
          (!compact_handoff &&
           required_string(handoff, "native_match_provenance_identity") !=
              required_string(match, "provenance_identity")))
        throw std::invalid_argument(
            "consumed transport-arm planned hop differs from its plan");
      const auto& producing = as_object(
          handoff.at("producing"),
          "consumed transport-arm producing handoff");
      const auto& incoming = as_object(
          producing.at("incoming"),
          "consumed transport-arm incoming handoff");
      if (required_string(incoming, "local") !=
              tile_sources_[index]->handle() ||
          required_string(incoming, "checkpoint_identity") !=
              tile_sources_[index]->checkpoint_identity() ||
          (bounded_handoff
               ? !compact_matching_identity_reference_matches(
                     incoming.at("source_operator_reference"),
                     tile_sources_[index]->source_operator_identity(),
                     "consumed incoming source-operator reference")
               : required_string(incoming, "source_operator_identity") !=
                     tile_sources_[index]->source_operator_identity()))
        throw std::invalid_argument(
            "consumed transport-arm incoming lineage is inconsistent");
      const auto& producing_chart =
          retained.charts.at(retained.exact.matches[index].producing_chart);
      const auto& receiving_chart =
          retained.charts.at(retained.exact.matches[index].receiving_chart);
      if (bounded_handoff &&
          (!compact_matching_identity_reference_matches(
               producing.at("chart_identity_reference"),
               producing_chart.exact_identity,
               "consumed producing chart-identity reference") ||
           !compact_matching_identity_reference_matches(
               as_object(handoff.at("receiving"),
                         "consumed transport-arm receiving handoff")
                   .at("chart_identity_reference"),
               receiving_chart.exact_identity,
               "consumed receiving chart-identity reference")))
        throw std::invalid_argument(
            "consumed transport-arm compact chart identity is inconsistent");
      const auto& handoff_basis = as_array(
          as_object(handoff.at("receiving"),
                    "consumed transport-arm receiving handoff")
              .at("basis"),
          "consumed transport-arm receiving handoff basis");
      if (handoff_basis.size() != raw_basis.size())
        throw std::invalid_argument(
            "consumed transport-arm cached basis dimension changed");
      for (std::size_t column = 0; column < raw_basis.size(); ++column) {
        const auto& cached = as_object(
            raw_basis[column], "consumed transport-arm cached basis column");
        const auto& handed = as_object(
            handoff_basis[column],
            "consumed transport-arm handed-off basis column");
        if (as_u64(handed.at("column"),
                   "consumed transport-arm basis column") != column ||
            handed.at("local") != cached.at("local") ||
            handed.at("checkpoint_identity") !=
                cached.at("checkpoint_identity") ||
            (bounded_handoff
                 ? !compact_matching_identity_reference_matches(
                       handed.at("source_operator_reference"),
                       tile_sources_[index + 1]
                           ->source_operator_identity(),
                       "consumed handed basis source-operator reference")
                 : handed.at("source_operator_identity") !=
                       cached.at("source_operator_identity")))
          throw std::invalid_argument(
              "consumed transport-arm cached basis differs from its planned hop");
      }

      const auto& local = tile_sources_[index + 1];
      if (!local || !local->has_sealed_plan_match_lineage() ||
          local->retained_derivation_owner() != nullptr ||
          !local->retained_derivation().has_value())
        throw std::invalid_argument(
            "consumed transport-arm materialized local lineage is not sealed");
      const auto& derivation = *local->retained_derivation();
      const auto& sealed = as_object(
          match.at("sealed_local_lineage"),
          "consumed transport-arm sealed local lineage");
      const auto equation_owner = local->retained_equation_owner();
      const auto sealed_schema = required_string(sealed, "schema");
      if (sealed_schema ==
          "diffexp3-sealed-plan-match-local-lineage-v2") {
        require_exact_keys(
            sealed,
            {"schema", "local", "local_checkpoint_identity",
             "source_operator_identity", "match",
             "match_checkpoint_identity", "tile_plan",
             "tile_plan_checkpoint_identity", "arm", "match_index",
             "incoming_checkpoint_identity",
             "equation_owner_signature_identity",
             "equation_payload_identity"},
            "consumed compact transport-arm sealed local lineage");
        if (required_string(derivation, "schema") !=
                "diffexp3-retained-plan-match-local-materialization-v2" ||
            required_string(sealed, "local") != local->handle() ||
            required_string(sealed, "local_checkpoint_identity") !=
                local->checkpoint_identity() ||
            required_string(sealed, "source_operator_identity") !=
                local->source_operator_identity() ||
            sealed.at("match") != derivation.at("source_match") ||
            sealed.at("match_checkpoint_identity") !=
                match.at("checkpoint_identity") ||
            sealed.at("match_checkpoint_identity") !=
                derivation.at("source_match_checkpoint_identity") ||
            sealed.at("tile_plan") != derivation.at("tile_plan") ||
            required_string(sealed, "tile_plan") !=
                plan_owner_->handle() ||
            sealed.at("tile_plan_checkpoint_identity") !=
                derivation.at("tile_plan_checkpoint_identity") ||
            required_string(sealed, "tile_plan_checkpoint_identity") !=
                plan_owner_->checkpoint_identity() ||
            sealed.at("arm") != derivation.at("arm") ||
            required_string(sealed, "arm") != arm_ ||
            sealed.at("match_index") != derivation.at("match") ||
            as_u64(sealed.at("match_index"),
                   "consumed compact match index") != index ||
            sealed.at("incoming_checkpoint_identity") !=
                as_object(derivation.at("incoming"),
                          "consumed compact derivation incoming")
                    .at("checkpoint_identity") ||
            required_string(sealed, "incoming_checkpoint_identity") !=
                tile_sources_[index]->checkpoint_identity() ||
            derivation.at("basis") != handoff_basis ||
            !equation_owner ||
            required_string(sealed,
                "equation_owner_signature_identity") !=
                equation_owner->owner_signature_identity() ||
            required_string(sealed, "equation_payload_identity") !=
                equation_owner->physical_payload_identity())
          throw std::invalid_argument(
              "consumed compact transport-arm sealed lineage is inconsistent");
      } else {
        require_exact_keys(
            sealed,
            {"schema", "local", "local_checkpoint_identity",
             "source_operator_identity", "match",
             "match_checkpoint_identity", "match_provenance_identity",
             "planned_hop_provenance_identity",
             "derivation_provenance_identity",
             "equation_owner_signature_identity",
             "equation_payload_identity"},
            "consumed transport-arm sealed local lineage");
        if (sealed_schema !=
                "diffexp3-sealed-plan-match-local-lineage-v1" ||
          required_string(sealed, "local") != local->handle() ||
          required_string(sealed, "local_checkpoint_identity") !=
              local->checkpoint_identity() ||
          required_string(sealed, "source_operator_identity") !=
              local->source_operator_identity() ||
          sealed.at("match") != derivation.at("source_match") ||
          sealed.at("match_checkpoint_identity") !=
              match.at("checkpoint_identity") ||
          sealed.at("match_checkpoint_identity") !=
              derivation.at("source_match_checkpoint_identity") ||
          sealed.at("match_provenance_identity") !=
              match.at("provenance_identity") ||
          sealed.at("match_provenance_identity") !=
              derivation.at("source_match_provenance_identity") ||
          sealed.at("planned_hop_provenance_identity") !=
              derivation.at("planned_hop_provenance_identity") ||
          required_string(
              derivation, "planned_hop_provenance_identity") !=
              json::serialize(canonical_json_value(handoff)) ||
          sealed.at("derivation_provenance_identity") !=
              derivation.at("provenance_identity") ||
          derivation.at("planned_hop") != match.at("planned_hop") ||
          !equation_owner ||
          required_string(sealed, "equation_owner_signature_identity") !=
              equation_owner->owner_signature_identity() ||
          required_string(sealed, "equation_payload_identity") !=
              equation_owner->physical_payload_identity())
          throw std::invalid_argument(
              "consumed transport-arm sealed lineage identity is inconsistent");
      }
    }
  }

  std::string handle_;
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::string arm_;
  std::shared_ptr<StoredTilePlan> plan_owner_;
  std::shared_ptr<StoredLocalBase> anchor_owner_;
  std::vector<std::vector<std::shared_ptr<StoredLocalBase>>> basis_owners_;
  std::vector<std::shared_ptr<StoredPlannedMatchHop>> matches_;
  std::vector<std::shared_ptr<StoredLocalBase>> tile_sources_;
  std::shared_ptr<StoredPlannedMatchHop>
      terminal_factorized_match_;
  EpsilonWindow work_epsilon_;
  std::int32_t public_required_complete_max_ = 0;
  std::int32_t match_required_complete_max_ = 0;
  json::object refinement_;
  double elapsed_ms_ = 0.0;
  bool compact_provenance_ = true;
  bool consumed_compact_ = false;
  bool consumed_certificate_only_ = false;
  json::array cached_basis_references_;
  json::array cached_match_references_;
  mutable std::atomic<std::uint64_t> stats_queries_{0};
  std::atomic<std::uint64_t> contraction_operations_{0};
  std::atomic<std::uint64_t> contracted_observables_{0};
  std::atomic<std::uint64_t> endpoint_batch_operations_{0};
  std::atomic<std::uint64_t> endpoint_rows_{0};
};

void require_transport_pair_compatibility(
    const std::shared_ptr<StoredTransportArmState>& lower,
    const std::shared_ptr<StoredTransportArmState>& upper,
    const std::string& expected_domain) {
  if (!lower || !upper || lower.get() == upper.get())
    throw std::invalid_argument(
        "retained transport pair requires two distinct state owners");
  if (lower->arm_name() != "lower" || upper->arm_name() != "upper")
    throw std::invalid_argument(
        "retained transport pair requires lower and upper state directions");
  if (expected_domain != "rational" && expected_domain != "acb")
    throw std::invalid_argument(
        "retained transport pair requires one numeric session domain");
  const auto& lower_anchor = lower->anchor_owner();
  const auto& upper_anchor = upper->anchor_owner();
  if (!lower_anchor || !upper_anchor ||
      lower_anchor.get() != upper_anchor.get() ||
      lower_anchor->checkpoint_identity() !=
          upper_anchor->checkpoint_identity() ||
      lower_anchor->source_operator_identity() !=
          upper_anchor->source_operator_identity() ||
      std::string(lower_anchor->scalar_domain()) != expected_domain ||
      std::string(upper_anchor->scalar_domain()) != expected_domain ||
      lower_anchor->exact_analytic_metadata() !=
          upper_anchor->exact_analytic_metadata())
    throw std::invalid_argument(
        "retained transport pair does not share one exact anchor owner, checkpoint, operator, domain, and analytic metadata");
  const auto& lower_plan = lower->plan_owner()->arm("lower").exact;
  const auto& upper_plan = upper->plan_owner()->arm("upper").exact;
  if (lower_plan.direction != -1 || upper_plan.direction != 1 ||
      lower_plan.from.str() != upper_plan.from.str())
    throw std::invalid_argument(
        "retained transport pair plans do not share one compatible exact physical anchor");
}

struct ResolvedTransportEndpointBinding {
  json::object source;
  std::string arm;
  Rational local_end;
  std::int32_t approach_direction = 0;
  std::optional<std::int32_t> rim;
  bool centered = false;
};

ResolvedTransportEndpointBinding resolve_transport_endpoint_binding(
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::string& arm_name,
    const std::shared_ptr<StoredLocalBase>& local) {
  if (!plan || !local)
    throw std::invalid_argument(
        "transport endpoint binding requires a retained plan and terminal local");
  const auto& arm = plan->arm(arm_name);
  if (arm.exact.tiles.empty() || arm.charts.empty())
    throw std::invalid_argument(
        "transport endpoint arm has no final tile/chart");
  const auto final_tile_index = arm.exact.tiles.size() - 1;
  const auto& final_tile = arm.exact.tiles.back();
  if (final_tile.chart >= arm.charts.size())
    throw std::logic_error(
        "transport endpoint final tile has an invalid chart index");
  const auto& final_chart = arm.charts[final_tile.chart];
  if (!(final_tile.physical_end == arm.exact.to) ||
      final_chart.geometry.scale.is_zero() ||
      arm.exact.direction != (arm_name == "lower" ? -1 : 1))
    throw std::invalid_argument(
        "transport endpoint plan has an inconsistent final physical binding");
  const auto mapped_endpoint = final_chart.geometry.center +
      final_chart.geometry.scale * final_tile.local_end;
  if (!(mapped_endpoint == arm.exact.to))
    throw std::invalid_argument(
        "transport endpoint final local coordinate does not map to its exact physical endpoint");
  if (!local || local->source_chart() != final_chart.handle)
    throw std::invalid_argument(
        "transport endpoint final local does not name its final retained chart");
  local->require_exact_plan_binding(
      final_chart.local_geometry, final_chart.prescriptions,
      "transport endpoint final local");

  ResolvedTransportEndpointBinding binding;
  binding.arm = arm_name;
  binding.local_end = final_tile.local_end;
  binding.centered = final_tile.local_end.is_zero();
  binding.approach_direction =
      -arm.exact.direction * final_chart.geometry.scale.sign();
  binding.rim = exact_plan_rim(
      final_chart.prescriptions, final_chart.geometry.scale);
  binding.source = json::object{
      {"tile_plan", json::object{
           {"handle", plan->handle()},
           {"checkpoint_identity", plan->checkpoint_identity()},
           {"provenance_identity", plan->provenance_identity()}}},
      {"arm", arm_name},
      {"endpoint_exact", arm.exact.to.str()},
      {"direction", arm.exact.direction},
      {"final_tile", final_tile_index},
      {"final_chart_index", final_tile.chart},
      {"final_chart", final_chart.handle},
      {"final_chart_identity", final_chart.exact_identity},
      {"local_endpoint_exact", final_tile.local_end.str()},
      {"centered", binding.centered},
      {"approach_direction", binding.approach_direction},
      {"derived_rim", binding.rim.has_value()
           ? json::value(*binding.rim) : json::value(nullptr)},
      {"final_local", json::object{
           {"handle", local->handle()}, {"chart", local->source_chart()},
           {"source_operator_identity", local->source_operator_identity()},
           {"checkpoint_identity", local->checkpoint_identity()},
           {"coefficient_domain", local->scalar_domain()}}},
      {"prescriptions", encode_plan_prescriptions(
           final_chart.prescriptions)}};
  return binding;
}

ResolvedTransportEndpointBinding resolve_transport_endpoint_binding(
    const std::shared_ptr<StoredTransportArmState>& state) {
  if (!state)
    throw std::invalid_argument(
        "transport endpoint binding requires a retained arm state");
  auto binding = resolve_transport_endpoint_binding(
      state->plan_owner(), state->arm_name(), state->final_local());
  binding.source["transport_state"] = json::object{
      {"handle", state->handle()},
      {"checkpoint_identity", state->checkpoint_identity()},
      {"provenance_identity", state->provenance_identity()}};
  return binding;
}

json::object rolling_transport_endpoint_source(
    const ResolvedTransportEndpointBinding& binding,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& local) {
  if (!plan || !local)
    throw std::invalid_argument(
        "rolling transport endpoint source requires retained plan/local owners");
  const auto& source = binding.source;
  return json::object{
      {"binding_mode", "rolling-terminal-local"},
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
      {"tile_plan_provenance_identity", plan->provenance_identity()},
      {"arm", binding.arm},
      {"endpoint_exact", source.at("endpoint_exact")},
      {"direction", source.at("direction")},
      {"final_tile", source.at("final_tile")},
      {"final_chart_index", source.at("final_chart_index")},
      {"final_chart", source.at("final_chart")},
      {"final_chart_identity", source.at("final_chart_identity")},
      {"local_endpoint_exact", binding.local_end.str()},
      {"centered", binding.centered},
      {"approach_direction", binding.approach_direction},
      {"derived_rim", binding.rim.has_value()
           ? json::value(*binding.rim) : json::value(nullptr)},
      {"local", local->handle()},
      {"chart", local->source_chart()},
      {"source_operator_identity", local->source_operator_identity()},
      {"checkpoint_identity", local->checkpoint_identity()},
      {"coefficient_domain", local->scalar_domain()},
      {"prescriptions", source.at("prescriptions")}};
}

// Transport endpoint results strongly own the exact plan/local/state records.
// Their hot provenance therefore binds those owners by deterministic content
// references instead of recursively embedding multi-megabyte identities.
// The full prepared row is stored separately in the v2 native checkpoint and
// validated against the row references during restore.
json::object bounded_transport_endpoint_source(
    const ResolvedTransportEndpointBinding& binding,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& local,
    const std::shared_ptr<StoredTransportArmState>& state,
    bool rolling) {
  json::object result = rolling
      ? rolling_transport_endpoint_source(binding, plan, local)
      : binding.source;
  const auto replace_identity = [](json::object& object,
                                   const char* identity_key,
                                   const char* reference_key) {
    const auto identity = required_string(object, identity_key);
    object.erase(identity_key);
    object[reference_key] = compact_matching_identity_reference(identity);
  };
  if (rolling) {
    replace_identity(result, "tile_plan_provenance_identity",
                     "tile_plan_provenance_reference");
    replace_identity(result, "source_operator_identity",
                     "source_operator_reference");
    return result;
  }
  auto& plan_reference = result.at("tile_plan").as_object();
  replace_identity(plan_reference, "provenance_identity",
                   "provenance_reference");
  auto& local_reference = result.at("final_local").as_object();
  replace_identity(local_reference, "source_operator_identity",
                   "source_operator_reference");
  auto& state_reference = result.at("transport_state").as_object();
  replace_identity(state_reference, "provenance_identity",
                   "provenance_reference");
  if (!state)
    throw std::invalid_argument(
        "bounded state-backed transport endpoint lost its state owner");
  return result;
}

void validate_prepared_rational_row_structure(
    const json::value& raw, std::uint32_t expected_columns,
    const char* label) {
  const auto& row = as_object(raw, label);
  require_exact_keys(row,
                     {"schema", "columns", "exact_identity", "entries"},
                     label);
  if (required_string(row, "schema") !=
          "diffexp3-prepared-rational-local-row-v1" ||
      as_u32(row.at("columns"), label) != expected_columns ||
      required_string(row, "exact_identity").empty())
    throw std::invalid_argument(
        std::string(label) + " has an incompatible schema, dimension, or identity");
  std::optional<std::uint32_t> previous;
  for (const auto& raw_entry : as_array(row.at("entries"), label)) {
    const auto& entry = as_object(raw_entry, label);
    require_exact_keys(entry, {"column", "multiplier"}, label);
    const auto column = as_u32(entry.at("column"), label);
    if (column >= expected_columns ||
        (previous.has_value() && *previous >= column))
      throw std::invalid_argument(
          std::string(label) + " columns are not strictly ordered in range");
    previous = column;
    const auto& multiplier = as_object(entry.at("multiplier"), label);
    const bool has_kernels =
        multiplier.if_contains("kernels") != nullptr;
    const bool has_analytic_coefficients =
        multiplier.if_contains("analytic_coefficients") != nullptr;
    if (has_kernels && has_analytic_coefficients)
      require_exact_keys(
          multiplier,
          {"epsilon_shift", "center_pole_order", "kernels",
           "exact_identity", "proven_zero", "analytic_coefficients"},
          label);
    else if (has_kernels)
      require_exact_keys(
          multiplier,
          {"epsilon_shift", "center_pole_order", "kernels",
           "exact_identity", "proven_zero"}, label);
    else if (has_analytic_coefficients)
      require_exact_keys(
          multiplier,
          {"epsilon_shift", "center_pole_order", "exact_identity",
           "proven_zero", "analytic_coefficients"}, label);
    else
      throw std::invalid_argument(
          std::string(label) +
          " multiplier needs Taylor kernels or an analytic rational source");
    (void)as_i32(multiplier.at("epsilon_shift"), label);
    (void)as_u32(multiplier.at("center_pole_order"), label);
    if (required_string(multiplier, "exact_identity").empty() ||
        !multiplier.at("proven_zero").is_bool() ||
        multiplier.at("proven_zero").as_bool())
      throw std::invalid_argument(
          std::string(label) + " contains an invalid active multiplier");
    std::size_t kernel_count = 0;
    if (has_kernels) {
      const auto& kernels = as_array(multiplier.at("kernels"), label);
      if (kernels.empty())
        throw std::invalid_argument(
            std::string(label) + " contains no epsilon kernels");
      kernel_count = kernels.size();
      for (const auto& raw_kernel : kernels)
        if (as_array(raw_kernel, label).empty())
          throw std::invalid_argument(
              std::string(label) + " contains an empty Taylor kernel");
    }
    if (has_analytic_coefficients) {
      const auto& coefficients = as_array(
          multiplier.at("analytic_coefficients"), label);
      if (coefficients.empty() ||
          (has_kernels && coefficients.size() != kernel_count))
        throw std::invalid_argument(
            std::string(label) +
            " analytic coefficient count differs from its epsilon kernels");
      for (const auto& raw_coefficient : coefficients) {
        const auto& coefficient = as_object(raw_coefficient, label);
        require_exact_keys(coefficient, {"numerator", "denominator"}, label);
        if (as_array(coefficient.at("numerator"), label).empty() ||
            as_array(coefficient.at("denominator"), label).empty())
          throw std::invalid_argument(
              std::string(label) +
              " contains an empty analytic numerator or denominator");
      }
    }
  }
}

json::array bounded_line_aggregate_source_records(
    const std::vector<std::shared_ptr<StoredLocalBase>>& owners);

class StoredLineResult {
 public:
  StoredLineResult(std::string handle, std::string checkpoint_identity,
                   std::string provenance_identity, std::string arm,
                   std::size_t tile_index, json::object interval,
                   std::string source_checkpoint,
                   StoredLineIntegral result, double elapsed_ms,
                   std::shared_ptr<StoredTilePlan> plan_owner,
                   std::shared_ptr<StoredLocalBase> local_owner)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        arm_(std::move(arm)), tile_index_(tile_index),
        interval_(std::move(interval)),
        source_checkpoint_(std::move(source_checkpoint)),
        result_(std::move(result)), elapsed_ms_(elapsed_ms),
        plan_owners_{std::move(plan_owner)},
        local_owners_{std::move(local_owner)} {}

  StoredLineResult(std::string handle, std::string checkpoint_identity,
                   std::string provenance_identity,
                   StoredLineIntegral result, double elapsed_ms,
                   std::shared_ptr<StoredTilePlan> plan_owner,
                   std::vector<std::shared_ptr<StoredLocalBase>> local_owners,
                   json::object aggregate_provenance,
                   std::shared_ptr<StoredTransportArmState> transport_owner =
                       nullptr)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        result_(std::move(result)), elapsed_ms_(elapsed_ms),
        plan_owners_{std::move(plan_owner)},
        local_owners_(std::move(local_owners)),
        aggregate_provenance_(std::move(aggregate_provenance)) {
    if (transport_owner) transport_owners_.push_back(
        std::move(transport_owner));
    validate_aggregate_ownership();
  }

  StoredLineResult(
      std::string handle, std::string checkpoint_identity,
      std::string provenance_identity, StoredLineIntegral result,
      double elapsed_ms,
      std::vector<std::shared_ptr<StoredTilePlan>> plan_owners,
      std::vector<std::shared_ptr<StoredLocalBase>> local_owners,
      json::object aggregate_provenance,
      std::vector<std::shared_ptr<StoredTransportArmState>>
          transport_owners)
      : handle_(std::move(handle)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        provenance_identity_(std::move(provenance_identity)),
        result_(std::move(result)), elapsed_ms_(elapsed_ms),
        plan_owners_(std::move(plan_owners)),
        local_owners_(std::move(local_owners)),
        aggregate_provenance_(std::move(aggregate_provenance)),
        transport_owners_(std::move(transport_owners)) {
    validate_aggregate_ownership();
  }

  const std::string& handle() const { return handle_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& provenance_identity() const {
    return provenance_identity_;
  }
  double elapsed_ms() const { return elapsed_ms_; }
  const StoredLineIntegral& result() const { return result_; }
  bool is_aggregate() const { return aggregate_provenance_.has_value(); }

  const std::shared_ptr<StoredTilePlan>& plan_owner() const {
    if (plan_owners_.size() != 1)
      throw std::logic_error(
          "paired compact line has no singular retained plan owner");
    return plan_owners_.front();
  }
  const std::vector<std::shared_ptr<StoredTilePlan>>& plan_owners() const {
    return plan_owners_;
  }
  const std::shared_ptr<StoredLocalBase>& local_owner() const {
    if (local_owners_.empty())
      throw std::logic_error(
          "compact transport observable line has no retained local owner");
    return local_owners_.front();
  }
  const std::vector<std::shared_ptr<StoredLocalBase>>& local_owners() const {
    return local_owners_;
  }
  const std::shared_ptr<StoredTransportArmState>& transport_owner() const {
    if (transport_owners_.size() > 1)
      throw std::logic_error(
          "paired compact line has no singular transport-state owner");
    static const std::shared_ptr<StoredTransportArmState> empty;
    return transport_owners_.empty() ? empty : transport_owners_.front();
  }
  const std::vector<std::shared_ptr<StoredTransportArmState>>&
  transport_owners() const {
    return transport_owners_;
  }

  json::object summary() const {
    const auto& diagnostics = result_.diagnostics;
    json::object output{
        {"line", handle_},
        {"capability", transport_owners_.size() == 2
             ? kRetainedTransportPairContractionCapability
             : transport_owners_.size() == 1
             ? kRetainedTransportArmContractionCapability
             : is_aggregate()
             ? kRetainedLineAggregateCapability
             : result_.scope ==
                           LineIntegrationScope::FullLocalWithCertifiedTail
             ? kRetainedCertifiedLineCapability
             : kRetainedStoredLineCapability},
        {"native_retained", true}, {"json_coefficients", 0},
        {"scope", line_integration_scope_name(result_.scope)},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"dimension", result_.value.dimension},
        {"epsilon_min", result_.value.epsilon.min_power},
        {"epsilon_max", result_.value.epsilon.complete_max},
        {"effective_rim", result_.imaginary_sign.has_value()
             ? json::value(*result_.imaginary_sign) : json::value(nullptr)},
        {"error", encode_error_envelope_summary(result_.value.error)},
        {"diagnostics", json::object{
             {"input_monomial_cells", diagnostics.input_monomial_cells},
             {"grouped_monomials", diagnostics.grouped_monomials},
             {"zero_groups_skipped", diagnostics.zero_groups_skipped},
             {"cancelled_divergent_groups",
              diagnostics.cancelled_divergent_groups},
             {"bounded_cancelled_divergent_coefficients",
              diagnostics.bounded_cancelled_divergent_coefficients},
             {"divergent_cancellation_mode",
              diagnostics.divergent_cancellation_mode},
             {"divergent_relative_tolerance",
              diagnostics.divergent_relative_tolerance.empty()
                  ? json::value(nullptr)
                  : json::value(
                        diagnostics.divergent_relative_tolerance)},
             {"divergent_cancellation_provenance",
              diagnostics.divergent_cancellation_provenance.empty()
                  ? json::value(nullptr)
                  : json::value(
                        diagnostics.divergent_cancellation_provenance)},
             {"primitive_evaluations", diagnostics.primitive_evaluations},
             {"primitive_component_applications",
              diagnostics.primitive_component_applications},
             {"primitive_component_reuses",
              diagnostics.primitive_component_reuses},
             {"has_center_endpoint", diagnostics.has_center_endpoint},
             {"tail_certificate_requested",
              diagnostics.tail_certificate_requested},
             {"tail_certificate_status",
              diagnostics.tail_certificate_status},
             {"tail_witness_radius_exact",
              diagnostics.tail_witness_radius_exact.empty()
                  ? json::value(nullptr)
                  : json::value(
                        diagnostics.tail_witness_radius_exact)},
             {"detail", diagnostics.detail}}},
        {"elapsed_ms", elapsed_ms_}};
    if (aggregate_provenance_.has_value()) {
      output["source"] = aggregate_provenance_->at("source");
      output["arm"] = aggregate_provenance_->at("arm");
      output["tile"] = nullptr;
      output["interval"] = aggregate_provenance_->at("interval");
      output["aggregate"] = aggregate_provenance_->at("aggregate");
    } else {
      const auto& owner = local_owners_.front();
      output["source"] = json::object{
          {"tile_plan", plan_owners_.front()->handle()},
          {"tile_plan_checkpoint_identity",
           plan_owners_.front()->checkpoint_identity()},
          {"local", owner->handle()},
          {"chart", owner->source_chart()},
          {"local_checkpoint_identity", source_checkpoint_}};
      output["arm"] = arm_;
      output["tile"] = tile_index_;
      output["interval"] = interval_;
    }
    return output;
  }

  json::object stats_json() const {
    auto result = summary();
    std::lock_guard<std::mutex> lock(stats_mutex_);
    result["exports"] = exports_;
    result["export_ms"] = export_ms_;
    return result;
  }

  json::object export_values(const std::string& expected_checkpoint,
                             int output_digits) {
    if (expected_checkpoint != checkpoint_identity_)
      throw std::invalid_argument(
          "line export checkpoint identity does not match retained result");
    const auto started = std::chrono::steady_clock::now();
    auto value = encode_epsilon_vector(result_.value, output_digits);
    const auto elapsed = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    {
      std::lock_guard<std::mutex> lock(stats_mutex_);
      ++exports_;
      export_ms_ += elapsed;
    }
    return json::object{
        {"line", handle_}, {"checkpoint_identity", checkpoint_identity_},
        {"compatibility_export", true},
        {"scope", line_integration_scope_name(result_.scope)},
        {"error_guarantee",
         error_guarantee_name(result_.value.error.guarantee)},
        {"json_coefficients", value.at("coefficients").as_array().size()},
        {"value", std::move(value)}, {"elapsed_ms", elapsed}};
  }

  json::object checkpoint_record() const {
    const auto& diagnostics = result_.diagnostics;
    std::lock_guard<std::mutex> lock(stats_mutex_);
    if (aggregate_provenance_.has_value()) {
      const auto aggregate_schema = required_string(
          *aggregate_provenance_, "schema");
      const bool compact_v2 =
          aggregate_schema ==
              "diffexp3-retained-native-transport-observable-line-v2" ||
          aggregate_schema ==
              "diffexp3-retained-native-transport-pair-observable-line-v2";
      const bool bounded_aggregate_v2 = aggregate_schema ==
          "diffexp3-retained-native-line-aggregate-v2";
      return json::object{
          {"schema", transport_owners_.size() == 2
               ? compact_v2
                   ? "diffexp3-retained-transport-pair-observable-line-v2"
                   : "diffexp3-retained-transport-pair-observable-line-v1"
               : transport_owners_.size() == 1
               ? compact_v2
                   ? "diffexp3-retained-transport-observable-line-v2"
                   : "diffexp3-retained-transport-observable-line-v1"
               : bounded_aggregate_v2
                   ? "diffexp3-retained-line-aggregate-v2"
                   : "diffexp3-retained-line-aggregate-v1"},
          {"handle", handle_},
          {"checkpoint_identity", checkpoint_identity_},
          {"provenance_identity", provenance_identity_},
          {"provenance", *aggregate_provenance_},
          {"result",
           json::object{
               {"value", checkpoint_epsilon_vector_record(result_.value)},
               {"scope", line_integration_scope_name(result_.scope)},
               {"imaginary_sign", result_.imaginary_sign.has_value()
                    ? json::value(*result_.imaginary_sign)
                    : json::value(nullptr)},
               {"diagnostics",
                json::object{
                    {"input_monomial_cells",
                     diagnostics.input_monomial_cells},
                    {"grouped_monomials", diagnostics.grouped_monomials},
                    {"zero_groups_skipped",
                     diagnostics.zero_groups_skipped},
                    {"cancelled_divergent_groups",
                     diagnostics.cancelled_divergent_groups},
                    {"bounded_cancelled_divergent_coefficients",
                     diagnostics.bounded_cancelled_divergent_coefficients},
                    {"divergent_cancellation_mode",
                     diagnostics.divergent_cancellation_mode},
                    {"divergent_relative_tolerance",
                     diagnostics.divergent_relative_tolerance},
                    {"divergent_cancellation_provenance",
                     diagnostics.divergent_cancellation_provenance},
                    {"primitive_evaluations",
                     diagnostics.primitive_evaluations},
                    {"primitive_component_applications",
                     diagnostics.primitive_component_applications},
                    {"primitive_component_reuses",
                     diagnostics.primitive_component_reuses},
                    {"has_center_endpoint",
                     diagnostics.has_center_endpoint},
                    {"tail_certificate_requested",
                     diagnostics.tail_certificate_requested},
                    {"tail_certificate_status",
                     diagnostics.tail_certificate_status},
                    {"tail_witness_radius_exact",
                     diagnostics.tail_witness_radius_exact},
                    {"detail", diagnostics.detail}}}}},
          {"elapsed_ms", elapsed_ms_},
          {"runtime_stats",
           json::object{{"exports", exports_}, {"export_ms", export_ms_}}}};
    }
    const auto& owner = local_owners_.front();
    return json::object{
        {"schema", "diffexp3-retained-line-result-v2"},
        {"handle", handle_},
        {"checkpoint_identity", checkpoint_identity_},
        {"provenance_identity", provenance_identity_},
        {"arm", arm_},
        {"tile", tile_index_},
        {"interval", interval_},
        {"source",
         json::object{
             {"tile_plan", plan_owners_.front()->handle()},
             {"tile_plan_checkpoint_identity",
              plan_owners_.front()->checkpoint_identity()},
             {"local", owner->handle()},
             {"chart", owner->source_chart()},
             {"source_operator_identity",
              owner->source_operator_identity()},
             {"local_checkpoint_identity", source_checkpoint_},
             {"coefficient_domain", owner->scalar_domain()},
             {"analytic_metadata",
              owner->exact_analytic_metadata()}}},
        {"result",
         json::object{
             {"value", checkpoint_epsilon_vector_record(result_.value)},
             {"scope", line_integration_scope_name(result_.scope)},
             {"imaginary_sign", result_.imaginary_sign.has_value()
                  ? json::value(*result_.imaginary_sign)
                  : json::value(nullptr)},
             {"diagnostics",
              json::object{
                  {"input_monomial_cells",
                   diagnostics.input_monomial_cells},
                  {"grouped_monomials", diagnostics.grouped_monomials},
                  {"zero_groups_skipped", diagnostics.zero_groups_skipped},
                  {"cancelled_divergent_groups",
                   diagnostics.cancelled_divergent_groups},
                  {"bounded_cancelled_divergent_coefficients",
                   diagnostics.bounded_cancelled_divergent_coefficients},
                  {"divergent_cancellation_mode",
                   diagnostics.divergent_cancellation_mode},
                  {"divergent_relative_tolerance",
                   diagnostics.divergent_relative_tolerance},
                  {"divergent_cancellation_provenance",
                   diagnostics.divergent_cancellation_provenance},
                  {"primitive_evaluations",
                   diagnostics.primitive_evaluations},
                  {"primitive_component_applications",
                   diagnostics.primitive_component_applications},
                  {"primitive_component_reuses",
                   diagnostics.primitive_component_reuses},
                  {"has_center_endpoint",
                   diagnostics.has_center_endpoint},
                  {"tail_certificate_requested",
                   diagnostics.tail_certificate_requested},
                  {"tail_certificate_status",
                   diagnostics.tail_certificate_status},
                  {"tail_witness_radius_exact",
                   diagnostics.tail_witness_radius_exact.empty()
                       ? json::value(nullptr)
                       : json::value(
                             diagnostics.tail_witness_radius_exact)},
                  {"detail", diagnostics.detail}}}}},
        {"elapsed_ms", elapsed_ms_},
        {"runtime_stats",
         json::object{{"exports", exports_}, {"export_ms", export_ms_}}}};
  }

  void restore_runtime_stats(std::uint64_t exports, double export_ms) {
    std::lock_guard<std::mutex> lock(stats_mutex_);
    exports_ = exports;
    export_ms_ = export_ms;
  }

 private:
  static void validate_transport_reference(
      const json::value& raw,
      const std::shared_ptr<StoredTransportArmState>& owner,
      const char* label, bool compact_v2 = false) {
    const auto& reference = as_object(raw, label);
    if (compact_v2)
      require_exact_keys(reference, {"handle", "checkpoint_identity"},
                         label);
    else
      require_exact_keys(reference,
                         {"handle", "checkpoint_identity",
                          "provenance_identity"}, label);
    if (!owner || required_string(reference, "handle") != owner->handle() ||
        required_string(reference, "checkpoint_identity") !=
            owner->checkpoint_identity() ||
        (!compact_v2 &&
         required_string(reference, "provenance_identity") !=
             owner->provenance_identity()))
      throw std::invalid_argument(
          std::string(label) + " is stale");
  }

  void validate_aggregate_ownership() const {
    if (!aggregate_provenance_.has_value() || plan_owners_.empty() ||
        std::any_of(plan_owners_.begin(), plan_owners_.end(),
                    [](const auto& owner) { return owner == nullptr; }) ||
        std::any_of(local_owners_.begin(), local_owners_.end(),
                    [](const auto& owner) { return owner == nullptr; }) ||
        std::any_of(transport_owners_.begin(), transport_owners_.end(),
                    [](const auto& owner) { return owner == nullptr; }))
      throw std::invalid_argument(
          "retained line aggregate requires its exact strong owners");
    if (json::serialize(canonical_json_value(*aggregate_provenance_)) !=
        provenance_identity_)
      throw std::invalid_argument(
          "retained line aggregate provenance identity is inconsistent");
    const auto schema = required_string(
        *aggregate_provenance_, "schema");
    const auto& source = as_object(
        aggregate_provenance_->at("source"),
        "retained line aggregate source");
    const auto& aggregate = as_object(
        aggregate_provenance_->at("aggregate"),
        "retained line aggregate recipe");
    if (schema ==
            "diffexp3-retained-native-transport-observable-line-v1" ||
        schema ==
            "diffexp3-retained-native-transport-observable-line-v2") {
      const bool compact_v2 = schema.ends_with("-v2");
      if (plan_owners_.size() != 1 || transport_owners_.size() != 1 ||
          !local_owners_.empty())
        throw std::invalid_argument(
            "compact transport line requires one plan and state owner");
      validate_transport_reference(
          aggregate.at("transport_state"), transport_owners_.front(),
          "retained line aggregate transport reference", compact_v2);
      return;
    }
    if (schema ==
            "diffexp3-retained-native-transport-pair-observable-line-v1" ||
        schema ==
            "diffexp3-retained-native-transport-pair-observable-line-v2") {
      const bool compact_v2 = schema.ends_with("-v2");
      if (plan_owners_.size() != 2 || transport_owners_.size() != 2 ||
          !local_owners_.empty())
        throw std::invalid_argument(
            "compact transport-pair line requires two ordered plan and state owners");
      require_transport_pair_compatibility(
          transport_owners_[0], transport_owners_[1],
          transport_owners_[0]->anchor_owner()->scalar_domain());
      if (transport_owners_[0]->plan_owner().get() !=
              plan_owners_[0].get() ||
          transport_owners_[1]->plan_owner().get() !=
              plan_owners_[1].get())
        throw std::invalid_argument(
            "compact transport-pair plans differ from their ordered state owners");
      const auto& lower_source = as_object(
          source.at("lower"), "compact transport-pair lower source");
      const auto& upper_source = as_object(
          source.at("upper"), "compact transport-pair upper source");
      validate_transport_reference(
          lower_source.at("transport_state"), transport_owners_[0],
          "compact transport-pair lower state reference", compact_v2);
      validate_transport_reference(
          upper_source.at("transport_state"), transport_owners_[1],
          "compact transport-pair upper state reference", compact_v2);
      const auto& lower_recipe = as_object(
          aggregate.at("lower"), "compact transport-pair lower recipe");
      const auto& upper_recipe = as_object(
          aggregate.at("upper"), "compact transport-pair upper recipe");
      validate_transport_reference(
          lower_recipe.at("transport_state"), transport_owners_[0],
          "compact transport-pair lower recipe state", compact_v2);
      validate_transport_reference(
          upper_recipe.at("transport_state"), transport_owners_[1],
          "compact transport-pair upper recipe state", compact_v2);
      return;
    }
    if ((schema != "diffexp3-retained-native-line-aggregate-v1" &&
         schema != "diffexp3-retained-native-line-aggregate-v2") ||
        plan_owners_.size() != 1 || !transport_owners_.empty() ||
        local_owners_.empty())
      throw std::invalid_argument(
          "retained line aggregate ownership differs from its provenance schema");
    if (schema == "diffexp3-retained-native-line-aggregate-v2") {
      const auto& locals = as_object(
          aggregate_provenance_->at("source"),
          "bounded retained line aggregate source").at("locals");
      if (locals != bounded_line_aggregate_source_records(local_owners_))
        throw std::invalid_argument(
            "bounded retained line aggregate source differs from its strong owners");
    }
  }

  std::string handle_;
  std::string checkpoint_identity_;
  std::string provenance_identity_;
  std::string arm_;
  std::size_t tile_index_ = 0;
  json::object interval_;
  std::string source_checkpoint_;
  StoredLineIntegral result_;
  double elapsed_ms_ = 0.0;
  std::vector<std::shared_ptr<StoredTilePlan>> plan_owners_;
  std::vector<std::shared_ptr<StoredLocalBase>> local_owners_;
  std::optional<json::object> aggregate_provenance_;
  std::vector<std::shared_ptr<StoredTransportArmState>> transport_owners_;
  mutable std::mutex stats_mutex_;
  std::uint64_t exports_ = 0;
  double export_ms_ = 0.0;
};
