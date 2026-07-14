template <typename Scalar>
LocalSolution<Scalar> parse_checkpoint_local_solution(
    const json::value& raw) {
  const auto& object = as_object(raw, "checkpoint local solution");
  require_exact_keys(object,
      {"chart", "epsilon", "taylor_complete_max", "dimension", "sectors",
       "prescriptions", "error", "checkpoint_identity"},
      "checkpoint local solution");
  LocalSolution<Scalar> solution;
  const auto& chart = as_object(object.at("chart"),
                                "checkpoint local chart");
  require_exact_keys(chart,
      {"center_exact", "scale_exact", "radius_exact_ball",
       "infinite_radius"}, "checkpoint local chart");
  solution.chart.center_exact = required_string(chart, "center_exact");
  solution.chart.scale_exact = required_string(chart, "scale_exact");
  if (solution.chart.center_exact.empty() || solution.chart.scale_exact.empty())
    throw std::invalid_argument(
        "checkpoint local chart lost exact center or scale provenance");
  solution.chart.radius = parse_checkpoint_ball(
      chart.at("radius_exact_ball"), "checkpoint chart radius");
  if (!chart.at("infinite_radius").is_bool())
    throw std::invalid_argument(
        "checkpoint infinite-radius flag must be boolean");
  solution.chart.infinite_radius = chart.at("infinite_radius").as_bool();

  const auto& epsilon = as_object(object.at("epsilon"),
                                  "checkpoint local epsilon window");
  require_exact_keys(epsilon, {"min", "max"},
                     "checkpoint local epsilon window");
  solution.epsilon = {
      as_i32(epsilon.at("min"), "checkpoint local epsilon minimum"),
      as_i32(epsilon.at("max"), "checkpoint local epsilon maximum")};
  (void)solution.epsilon.width();
  solution.taylor_complete_max = as_u32(
      object.at("taylor_complete_max"), "checkpoint Taylor complete maximum");
  solution.dimension = as_u32(object.at("dimension"),
                              "checkpoint local dimension");
  if (solution.dimension == 0)
    throw std::invalid_argument("checkpoint local dimension is zero");
  const auto expected = checked_flat_count(
      checked_flat_count(solution.epsilon.width(), solution.taylor_width(),
                         "checkpoint local tensor"),
      solution.dimension, "checkpoint local tensor");

  for (const auto& raw_sector : as_array(object.at("sectors"),
                                          "checkpoint local sectors")) {
    const auto& sector_object = as_object(raw_sector,
                                          "checkpoint local sector");
    require_exact_keys(sector_object,
        {"a", "b", "log_power", "coefficients"},
        "checkpoint local sector");
    LocalSector<Scalar> sector;
    sector.a = parse_checkpoint_exact_descriptor(
        sector_object.at("a"), "checkpoint local a tag");
    sector.b = parse_checkpoint_exact_descriptor(
        sector_object.at("b"), "checkpoint local b tag");
    sector.log_power = as_u32(sector_object.at("log_power"),
                              "checkpoint local log power");
    const auto& coefficients = as_array(
        sector_object.at("coefficients"), "checkpoint sector coefficients");
    if (coefficients.size() != expected)
      throw std::invalid_argument(
          "checkpoint sector coefficient tensor has the wrong size");
    sector.coefficients.reserve(coefficients.size());
    for (const auto& coefficient : coefficients)
      sector.coefficients.push_back(parse_checkpoint_scalar<Scalar>(
          coefficient, "checkpoint local coefficient"));
    solution.sectors.push_back(std::move(sector));
  }
  if (solution.sectors.empty())
    throw std::invalid_argument("checkpoint local solution has no sectors");

  for (const auto& raw_prescription : as_array(
           object.at("prescriptions"), "checkpoint prescriptions")) {
    const auto& prescription = as_object(raw_prescription,
                                         "checkpoint prescription");
    require_exact_keys(prescription,
        {"factor_exact", "sign", "multiplicity",
         "leading_coefficient_sign"}, "checkpoint prescription");
    solution.prescriptions.push_back(Prescription{
        required_string(prescription, "factor_exact"),
        as_i32(prescription.at("sign"), "checkpoint prescription sign"),
        as_u32(prescription.at("multiplicity"),
               "checkpoint prescription multiplicity"),
        as_i32(prescription.at("leading_coefficient_sign"),
               "checkpoint leading-coefficient sign")});
  }
  solution.error = parse_checkpoint_error_envelope(object.at("error"));
  solution.checkpoint_identity = required_string(
      object, "checkpoint_identity");
  if (solution.checkpoint_identity.empty())
    throw std::invalid_argument(
        "checkpoint local solution identity is empty");
  validate_local_solution(solution, false);
  return solution;
}

template <typename Scalar>
std::vector<PseudoHit<Scalar>> parse_checkpoint_pseudo_hits(
    const json::value& raw, std::uint32_t dimension,
    std::uint32_t taylor_complete_max) {
  std::vector<PseudoHit<Scalar>> result;
  for (const auto& raw_hit : as_array(raw, "checkpoint pseudo hits")) {
    const auto& object = as_object(raw_hit, "checkpoint pseudo hit");
    require_exact_keys(object,
        {"n", "columns", "delta_b", "gamma_frames", "gamma_validity"},
        "checkpoint pseudo hit");
    PseudoHit<Scalar> hit;
    hit.n = as_u32(object.at("n"), "checkpoint pseudo-hit Taylor order");
    if (hit.n > taylor_complete_max)
      throw std::invalid_argument(
          "checkpoint pseudo hit lies above the retained Taylor window");
    std::set<std::uint32_t> unique_columns;
    for (const auto& raw_column : as_array(object.at("columns"),
                                            "checkpoint pseudo columns")) {
      const auto column = as_u32(raw_column, "checkpoint pseudo column");
      if (column >= dimension || !unique_columns.insert(column).second)
        throw std::invalid_argument(
            "checkpoint pseudo-hit columns are invalid or duplicated");
      hit.columns.push_back(column);
    }
    if (hit.columns.empty())
      throw std::invalid_argument(
          "checkpoint pseudo hit has an empty Jordan block");
    hit.delta_b = parse_checkpoint_scalar<Scalar>(
        object.at("delta_b"), "checkpoint pseudo delta-b");
    const auto& frames = as_array(object.at("gamma_frames"),
                                  "checkpoint pseudo gamma frames");
    const auto& validity = as_array(object.at("gamma_validity"),
                                    "checkpoint pseudo validity");
    if (frames.size() != hit.columns.size() ||
        validity.size() != hit.columns.size())
      throw std::invalid_argument(
          "checkpoint pseudo-hit block dimensions are inconsistent");
    std::optional<std::size_t> frame_width;
    for (const auto& raw_frame : frames) {
      const auto& coefficients = as_array(raw_frame,
                                          "checkpoint pseudo gamma frame");
      if (coefficients.empty() ||
          (frame_width.has_value() && *frame_width != coefficients.size()))
        throw std::invalid_argument(
            "checkpoint pseudo gamma frames have inconsistent widths");
      frame_width = coefficients.size();
      Frame<Scalar> frame;
      frame.reserve(coefficients.size());
      for (const auto& coefficient : coefficients)
        frame.push_back(parse_checkpoint_scalar<Scalar>(
            coefficient, "checkpoint pseudo gamma coefficient"));
      hit.gamma_frames.push_back(std::move(frame));
    }
    for (const auto& raw_validity : validity)
      hit.gamma_validity.push_back(parse_validity(raw_validity));
    result.push_back(std::move(hit));
  }
  return result;
}

template <typename Scalar>
std::shared_ptr<StoredLocalBase> restore_checkpoint_local_record(
    const json::value& raw, const std::string& expected_domain,
    slong expected_precision_bits,
    std::shared_ptr<void> retained_owner = nullptr,
    std::shared_ptr<PhysicalEquationOwnerBase> equation_owner = nullptr) {
  AcbPrecisionLease lease(expected_precision_bits);
  ComplexBall::set_precision(expected_precision_bits);
  const auto& object = as_object(raw, "checkpoint retained local");
  const bool has_tail_restore =
      object.if_contains("tail_model_restore") != nullptr;
  const bool has_derivation_record =
      object.if_contains("retained_derivation") != nullptr;
  const bool has_residual_restore =
      object.if_contains("residual_operator_restore") != nullptr;
  if (has_derivation_record !=
      (object.if_contains("retained_owner_lineage") != nullptr))
    throw std::invalid_argument(
        "checkpoint retained-local derivation fields are incomplete");
  const auto record_schema = required_string(object, "schema");
  const bool sealed_plan_match_lineage =
      record_schema == "diffexp2-retained-local-v5";
  if (!sealed_plan_match_lineage &&
      record_schema != "diffexp2-retained-local-v4")
    throw std::invalid_argument("unsupported retained-local checkpoint schema");
  if (!has_residual_restore)
    throw std::invalid_argument(
        "retained-local v4 checkpoint lacks its residual ownership record");
  std::set<std::string> actual_keys;
  for (const auto& entry : object) actual_keys.emplace(entry.key());
  std::set<std::string> expected_keys{
      "schema", "handle", "source_chart", "source_operator_identity",
      "scalar_domain", "precision_bits", "solution", "pseudo_hits",
      "diagnostics", "runtime_stats", "column_provenance",
      "equation_owner_restore"};
  if (sealed_plan_match_lineage)
    expected_keys.emplace("derivation_owner_restore");
  if (has_tail_restore) expected_keys.emplace("tail_model_restore");
  if (has_derivation_record) {
    expected_keys.emplace("retained_derivation");
    expected_keys.emplace("retained_owner_lineage");
  }
  if (has_residual_restore)
    expected_keys.emplace("residual_operator_restore");
  if (actual_keys != expected_keys)
    throw std::invalid_argument(
        "checkpoint retained local has unknown or missing fields");
  if (sealed_plan_match_lineage &&
      required_string(object, "derivation_owner_restore") !=
          "sealed-plan-match-lineage")
    throw std::invalid_argument(
        "checkpoint sealed local has an unsupported derivation-owner restore policy");
  const auto scalar_domain = required_string(object, "scalar_domain");
  const char* compile_domain = std::is_same_v<Scalar, Rational>
      ? "rational" : "acb";
  if (scalar_domain != compile_domain || scalar_domain != expected_domain)
    throw std::invalid_argument(
        "checkpoint local scalar domain is incompatible with its session");
  const auto precision_bits = as_i64(object.at("precision_bits"),
                                     "checkpoint local precision");
  if (precision_bits != expected_precision_bits)
    throw std::invalid_argument(
        "checkpoint local precision differs from its session precision");
  const auto handle = required_string(object, "handle");
  const auto source_chart = required_string(object, "source_chart");
  const auto source_operator_identity = required_string(
      object, "source_operator_identity");
  if (handle.empty() || source_chart.empty() ||
      source_operator_identity.empty())
    throw std::invalid_argument(
        "checkpoint local lost its handle or source-chart provenance");
  std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>>
      physical_equation;
  const bool saved_equation_owner =
      !object.at("equation_owner_restore").is_null();
  if (saved_equation_owner) {
    if (!equation_owner)
      throw std::invalid_argument(
          "checkpoint primitive local lost its physical equation owner");
    const auto& owner_record = as_object(
        object.at("equation_owner_restore"),
        "checkpoint local physical equation owner");
    require_exact_keys(
        owner_record,
        {"owner_kind", "owner_handle", "operator_identity",
         "owner_signature_identity", "physical_payload_identity"},
        "checkpoint local physical equation owner");
    if (required_string(owner_record, "owner_kind") !=
            equation_owner->equation_owner_kind() ||
        required_string(owner_record, "owner_handle") != source_chart ||
        equation_owner->equation_owner_handle() != source_chart ||
        required_string(owner_record, "operator_identity") !=
            source_operator_identity ||
        equation_owner->equation_operator_identity() !=
            source_operator_identity ||
        required_string(owner_record, "owner_signature_identity") !=
            equation_owner->owner_signature_identity() ||
        required_string(owner_record, "physical_payload_identity") !=
            equation_owner->physical_payload_identity() ||
        std::string(equation_owner->equation_scalar_domain()) !=
            expected_domain)
      throw std::invalid_argument(
          "checkpoint physical equation owner differs from its retained local provenance");
    const auto erased = equation_owner->physical_ode_erased();
    if (!erased)
      throw std::invalid_argument(
          "checkpoint physical equation owner has no typed q/C payload");
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
      throw std::invalid_argument(
          "checkpoint physical q/C payload did not round-trip byte-exactly through its owner");
  } else if (equation_owner) {
    throw std::invalid_argument(
        "checkpoint local unexpectedly acquired a physical equation owner");
  }
  auto solution = parse_checkpoint_local_solution<Scalar>(
      object.at("solution"));
  auto pseudo_hits = parse_checkpoint_pseudo_hits<Scalar>(
      object.at("pseudo_hits"), solution.dimension,
      solution.taylor_complete_max);
  std::optional<json::object> retained_derivation;
  const bool has_derivation = has_derivation_record &&
      !object.at("retained_derivation").is_null();
  if (has_derivation_record && has_derivation !=
      !object.at("retained_owner_lineage").is_null())
    throw std::invalid_argument(
        "checkpoint materialized-local derivation and owner lineage disagree");
  if (has_derivation) {
    if (retained_owner == nullptr && !sealed_plan_match_lineage)
      throw std::invalid_argument(
          "checkpoint derived local lost its strong owner");
    auto derivation = as_object(
        object.at("retained_derivation"),
        "checkpoint retained-local derivation");
    const auto derivation_schema = required_string(derivation, "schema");
    if (sealed_plan_match_lineage && derivation_schema !=
            "diffexp2-retained-plan-match-local-materialization-v1" &&
        derivation_schema !=
            "diffexp2-retained-plan-match-local-materialization-v2")
      throw std::invalid_argument(
          "checkpoint sealed lineage is not a plan-match materialization");
    if (derivation_schema ==
        "diffexp2-retained-rational-row-local-application-v1") {
      require_exact_keys(
          derivation,
          {"schema", "capability", "source", "row", "output",
           "analytic_prescriptions", "coefficient_transport",
           "provenance_identity"},
          "checkpoint rational-row local derivation");
      if (required_string(derivation, "capability") !=
              kRetainedRationalRowCapability ||
          required_string(derivation, "analytic_prescriptions") !=
              "preserved-exactly" ||
          required_string(derivation, "coefficient_transport") !=
              "native-retained-only")
        throw std::invalid_argument(
            "checkpoint rational-row derivation changes its retained scope");
      const auto& source = as_object(
          derivation.at("source"), "checkpoint rational-row source");
      require_exact_keys(
          source,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity", "dimension", "epsilon",
           "taylor_complete_max"},
          "checkpoint rational-row source");
      (void)scoped_handle_id(required_string(source, "local"), "l:",
                             "rational-row source local");
      const auto source_chart_identity = required_string(source, "chart");
      if ((!source_chart_identity.starts_with("c:") &&
           !source_chart_identity.starts_with("scc:")) ||
          source_chart_identity != source_chart)
        throw std::invalid_argument(
            "checkpoint rational-row source chart is inconsistent");
      const auto source_operator = required_string(
          source, "source_operator_identity");
      (void)required_string(source, "checkpoint_identity");
      const auto source_dimension = as_u32(
          source.at("dimension"), "checkpoint rational-row source dimension");
      if (source_dimension == 0)
        throw std::invalid_argument(
            "checkpoint rational-row source dimension is zero");
      const auto& source_epsilon = as_object(
          source.at("epsilon"), "checkpoint rational-row source epsilon");
      require_exact_keys(source_epsilon, {"min", "max"},
                         "checkpoint rational-row source epsilon");
      (void)EpsilonWindow{
          as_i32(source_epsilon.at("min"),
                 "checkpoint rational-row source epsilon minimum"),
          as_i32(source_epsilon.at("max"),
                 "checkpoint rational-row source epsilon maximum")}.width();
      (void)as_u32(source.at("taylor_complete_max"),
                   "checkpoint rational-row source Taylor maximum");

      const auto& row = as_object(
          derivation.at("row"), "checkpoint rational-row identity");
      require_exact_keys(
          row, {"exact_identity", "columns", "active_entries",
                "structurally_zero"},
          "checkpoint rational-row identity");
      const auto row_identity = required_string(row, "exact_identity");
      if (as_u32(row.at("columns"),
                 "checkpoint rational-row columns") != source_dimension ||
          !row.at("structurally_zero").is_bool())
        throw std::invalid_argument(
            "checkpoint rational-row dimension or zero fact is malformed");
      const auto& entries = as_array(
          row.at("active_entries"), "checkpoint rational-row entries");
      if (row.at("structurally_zero").as_bool() != entries.empty())
        throw std::invalid_argument(
            "checkpoint rational-row zero fact disagrees with its entries");
      std::optional<std::uint32_t> previous_column;
      for (const auto& raw_entry : entries) {
        const auto& entry = as_object(
            raw_entry, "checkpoint rational-row entry");
        require_exact_keys(
            entry, {"column", "epsilon_shift", "center_pole_order",
                    "exact_identity"},
            "checkpoint rational-row entry");
        const auto column = as_u32(
            entry.at("column"), "checkpoint rational-row entry column");
        if (column >= source_dimension ||
            (previous_column.has_value() && *previous_column >= column))
          throw std::invalid_argument(
              "checkpoint rational-row entry columns are not canonical");
        previous_column = column;
        (void)as_i32(entry.at("epsilon_shift"),
                     "checkpoint rational-row epsilon shift");
        (void)as_u32(entry.at("center_pole_order"),
                     "checkpoint rational-row pole order");
        (void)required_string(entry, "exact_identity");
      }

      const auto& output = as_object(
          derivation.at("output"), "checkpoint rational-row output");
      require_exact_keys(
          output, {"checkpoint_identity", "dimension", "epsilon",
                   "taylor_complete_max"},
          "checkpoint rational-row output");
      const auto& output_epsilon = as_object(
          output.at("epsilon"), "checkpoint rational-row output epsilon");
      require_exact_keys(output_epsilon, {"min", "max"},
                         "checkpoint rational-row output epsilon");
      if (required_string(output, "checkpoint_identity") !=
              solution.checkpoint_identity ||
          as_u32(output.at("dimension"),
                 "checkpoint rational-row output dimension") !=
              solution.dimension ||
          as_i32(output_epsilon.at("min"),
                 "checkpoint rational-row output epsilon minimum") !=
              solution.epsilon.min_power ||
          as_i32(output_epsilon.at("max"),
                 "checkpoint rational-row output epsilon maximum") !=
              solution.epsilon.complete_max ||
          as_u32(output.at("taylor_complete_max"),
                 "checkpoint rational-row output Taylor maximum") !=
              solution.taylor_complete_max)
        throw std::invalid_argument(
            "checkpoint rational-row output disagrees with its tensor");
      auto identity_input = derivation;
      const auto derivation_identity = required_string(
          derivation, "provenance_identity");
      identity_input.erase("provenance_identity");
      if (json::serialize(canonical_json_value(identity_input)) !=
          derivation_identity)
        throw std::invalid_argument(
            "checkpoint rational-row derivation identity is inconsistent");
      const json::object operator_provenance{
          {"schema", "diffexp2-rational-row-derived-operator-v1"},
          {"source_operator_identity", source_operator},
          {"row_exact_identity", row_identity},
          {"provenance_identity", derivation_identity}};
      if (json::serialize(canonical_json_value(operator_provenance)) !=
          source_operator_identity)
        throw std::invalid_argument(
            "checkpoint rational-row derived operator identity is inconsistent");
      const auto& lineage = as_object(
          object.at("retained_owner_lineage"),
          "checkpoint rational-row owner lineage");
      require_exact_keys(
          lineage,
          {"source_local", "source_chart", "source_operator_identity",
           "source_checkpoint_identity", "row_exact_identity",
           "derivation_provenance_identity", "derived_operator_identity"},
          "checkpoint rational-row owner lineage");
      if (lineage.at("source_local") != source.at("local") ||
          lineage.at("source_chart") != source.at("chart") ||
          lineage.at("source_operator_identity") !=
              source.at("source_operator_identity") ||
          lineage.at("source_checkpoint_identity") !=
              source.at("checkpoint_identity") ||
          lineage.at("row_exact_identity") != row.at("exact_identity") ||
          lineage.at("derivation_provenance_identity") !=
              derivation.at("provenance_identity") ||
          required_string(lineage, "derived_operator_identity") !=
              source_operator_identity)
        throw std::invalid_argument(
            "checkpoint rational-row owner lineage is inconsistent");
      retained_derivation = std::move(derivation);
    } else if (derivation_schema ==
        "diffexp2-retained-plan-match-local-materialization-v2") {
      require_exact_keys(
          derivation,
          {"schema", "capability", "source_match",
           "source_match_checkpoint_identity", "tile_plan",
           "tile_plan_checkpoint_identity", "arm", "match",
           "incoming", "basis", "weight_windows",
           "match_certified_complete_max", "output",
           "equation_owner_signature_identity",
           "equation_payload_identity", "scope",
           "coefficient_transport", "whole_arm_complete",
           "provenance_identity"},
          "checkpoint compact materialized-local derivation");
      if (required_string(derivation, "capability") !=
              "retained-native-plan-match-local-materialization-v1" ||
          required_string(derivation, "scope") !=
              "single-match-receiving-local" ||
          required_string(derivation, "coefficient_transport") !=
              "native-retained-only" ||
          !derivation.at("whole_arm_complete").is_bool() ||
          derivation.at("whole_arm_complete").as_bool() ||
          !equation_owner ||
          required_string(derivation,
              "equation_owner_signature_identity") !=
              equation_owner->owner_signature_identity() ||
          required_string(derivation, "equation_payload_identity") !=
              equation_owner->physical_payload_identity())
        throw std::invalid_argument(
            "checkpoint compact materialized-local derivation changes its certified scope or equation owner");
      (void)scoped_handle_id(required_string(derivation, "source_match"),
                             "m:", "compact materialized-local source match");
      (void)required_string(derivation,
                            "source_match_checkpoint_identity");
      (void)required_string(derivation, "tile_plan");
      (void)required_string(derivation,
                            "tile_plan_checkpoint_identity");
      const auto arm = required_string(derivation, "arm");
      if (arm != "lower" && arm != "upper")
        throw std::invalid_argument(
            "checkpoint compact materialized-local arm is invalid");
      (void)as_u64(derivation.at("match"),
                   "checkpoint compact materialized-local match index");
      const auto& incoming = as_object(
          derivation.at("incoming"),
          "checkpoint compact materialized-local incoming");
      require_exact_keys(
          incoming,
          {"local", "checkpoint_identity", "source_operator_identity"},
          "checkpoint compact materialized-local incoming");
      (void)scoped_handle_id(required_string(incoming, "local"), "l:",
                             "compact materialized-local incoming");
      (void)required_string(incoming, "checkpoint_identity");
      (void)required_string(incoming, "source_operator_identity");
      const auto& basis = as_array(
          derivation.at("basis"),
          "checkpoint compact materialized-local basis");
      if (basis.size() != solution.dimension || basis.empty())
        throw std::invalid_argument(
            "checkpoint compact materialized-local basis dimension changed");
      for (std::size_t column = 0; column < basis.size(); ++column) {
        const auto& source = as_object(
            basis[column], "checkpoint compact basis source");
        require_exact_keys(
            source,
            {"column", "local", "checkpoint_identity",
             "source_operator_identity"},
            "checkpoint compact basis source");
        if (as_u64(source.at("column"),
                   "checkpoint compact basis column") != column ||
            required_string(source, "source_operator_identity") !=
                source_operator_identity)
          throw std::invalid_argument(
              "checkpoint compact materialized-local basis order or operator changed");
        (void)scoped_handle_id(required_string(source, "local"), "l:",
                               "checkpoint compact basis local");
        (void)required_string(source, "checkpoint_identity");
      }
      const auto& windows = as_array(
          derivation.at("weight_windows"),
          "checkpoint compact materialization weight windows");
      if (windows.size() != basis.size())
        throw std::invalid_argument(
            "checkpoint compact materialization weight count changed");
      for (const auto& raw_window : windows) {
        const auto& window = as_object(
            raw_window, "checkpoint compact materialization window");
        require_exact_keys(window, {"min", "max"},
                           "checkpoint compact materialization window");
        (void)EpsilonWindow{
            as_i32(window.at("min"), "compact weight minimum"),
            as_i32(window.at("max"), "compact weight maximum")}.width();
      }
      (void)as_i32(derivation.at("match_certified_complete_max"),
                   "checkpoint compact certified maximum");
      const auto& output = as_object(
          derivation.at("output"),
          "checkpoint compact materialized-local output");
      require_exact_keys(
          output,
          {"checkpoint_identity", "chart", "source_operator_identity",
           "epsilon", "taylor_complete_max", "dimension"},
          "checkpoint compact materialized-local output");
      const auto& output_epsilon = as_object(
          output.at("epsilon"),
          "checkpoint compact materialized-local epsilon");
      require_exact_keys(output_epsilon, {"min", "max"},
                         "checkpoint compact materialized-local epsilon");
      if (required_string(output, "checkpoint_identity") !=
              solution.checkpoint_identity ||
          required_string(output, "chart") != source_chart ||
          required_string(output, "source_operator_identity") !=
              source_operator_identity ||
          as_i32(output_epsilon.at("min"), "compact output minimum") !=
              solution.epsilon.min_power ||
          as_i32(output_epsilon.at("max"), "compact output maximum") !=
              solution.epsilon.complete_max ||
          as_u32(output.at("taylor_complete_max"),
                 "compact output Taylor maximum") !=
              solution.taylor_complete_max ||
          as_u32(output.at("dimension"), "compact output dimension") !=
              solution.dimension)
        throw std::invalid_argument(
            "checkpoint compact materialized-local output disagrees with its tensor");
      auto identity_input = derivation;
      const auto derivation_identity = required_string(
          derivation, "provenance_identity");
      identity_input.erase("provenance_identity");
      if (json::serialize(canonical_json_value(identity_input)) !=
          derivation_identity)
        throw std::invalid_argument(
            "checkpoint compact materialized-local derivation identity is inconsistent");
      const auto& lineage = as_object(
          object.at("retained_owner_lineage"),
          "checkpoint compact materialized-local owner lineage");
      require_exact_keys(
          lineage,
          {"match", "match_checkpoint_identity", "tile_plan",
           "tile_plan_checkpoint_identity", "arm", "match_index"},
          "checkpoint compact materialized-local owner lineage");
      if (lineage.at("match") != derivation.at("source_match") ||
          lineage.at("match_checkpoint_identity") !=
              derivation.at("source_match_checkpoint_identity") ||
          lineage.at("tile_plan") != derivation.at("tile_plan") ||
          lineage.at("tile_plan_checkpoint_identity") !=
              derivation.at("tile_plan_checkpoint_identity") ||
          lineage.at("arm") != derivation.at("arm") ||
          lineage.at("match_index") != derivation.at("match"))
        throw std::invalid_argument(
            "checkpoint compact materialized-local owner lineage is inconsistent");
      retained_derivation = std::move(derivation);
    } else {
      require_exact_keys(
        derivation,
        {"schema", "capability", "source_match",
         "source_match_checkpoint_identity",
         "source_match_provenance_identity",
         "planned_hop_provenance_identity", "planned_hop",
         "weight_windows", "match_certified_complete_max", "output",
         "scope", "coefficient_transport", "whole_arm_complete",
         "provenance_identity"},
        "checkpoint materialized-local derivation");
      if (required_string(derivation, "schema") !=
            "diffexp2-retained-plan-match-local-materialization-v1" ||
        required_string(derivation, "capability") !=
            "retained-native-plan-match-local-materialization-v1" ||
        required_string(derivation, "scope") !=
            "single-match-receiving-local" ||
        required_string(derivation, "coefficient_transport") !=
            "native-retained-only" ||
        !derivation.at("whole_arm_complete").is_bool() ||
        derivation.at("whole_arm_complete").as_bool())
      throw std::invalid_argument(
          "checkpoint materialized-local derivation changes its certified scope");
    (void)scoped_handle_id(required_string(derivation, "source_match"),
                           "m:", "materialized-local source match");
    const auto& output = as_object(
        derivation.at("output"), "checkpoint materialized-local output");
    require_exact_keys(output,
        {"checkpoint_identity", "chart", "source_operator_identity",
         "epsilon", "taylor_complete_max", "dimension"},
        "checkpoint materialized-local output");
    const auto& output_epsilon = as_object(
        output.at("epsilon"), "checkpoint materialized-local epsilon");
    require_exact_keys(output_epsilon, {"min", "max"},
                       "checkpoint materialized-local epsilon");
    if (required_string(output, "checkpoint_identity") !=
            solution.checkpoint_identity ||
        required_string(output, "chart") != source_chart ||
        required_string(output, "source_operator_identity") !=
            source_operator_identity ||
        as_i32(output_epsilon.at("min"),
               "materialized-local epsilon minimum") !=
            solution.epsilon.min_power ||
        as_i32(output_epsilon.at("max"),
               "materialized-local epsilon maximum") !=
            solution.epsilon.complete_max ||
        as_u32(output.at("taylor_complete_max"),
               "materialized-local Taylor maximum") !=
            solution.taylor_complete_max ||
        as_u32(output.at("dimension"),
               "materialized-local dimension") != solution.dimension)
      throw std::invalid_argument(
          "checkpoint materialized-local output provenance disagrees with its tensor");
    for (const auto& raw_window : as_array(
             derivation.at("weight_windows"),
             "checkpoint materialization weight windows")) {
      const auto& window = as_object(
          raw_window, "checkpoint materialization weight window");
      require_exact_keys(window, {"min", "max"},
                         "checkpoint materialization weight window");
      (void)EpsilonWindow{
          as_i32(window.at("min"), "materialization weight minimum"),
          as_i32(window.at("max"), "materialization weight maximum")}.width();
    }
    auto identity_input = derivation;
    const auto derivation_identity = required_string(
        derivation, "provenance_identity");
    identity_input.erase("provenance_identity");
    if (json::serialize(canonical_json_value(identity_input)) !=
        derivation_identity)
      throw std::invalid_argument(
          "checkpoint materialized-local derivation identity is inconsistent");
    const auto& lineage = as_object(
        object.at("retained_owner_lineage"),
        "checkpoint materialized-local owner lineage");
    require_exact_keys(
        lineage,
        {"match", "match_checkpoint_identity", "match_provenance_identity",
         "planned_hop_provenance_identity",
         "derivation_provenance_identity"},
        "checkpoint materialized-local owner lineage");
    if (lineage.at("match") != derivation.at("source_match") ||
        lineage.at("match_checkpoint_identity") !=
            derivation.at("source_match_checkpoint_identity") ||
        lineage.at("match_provenance_identity") !=
            derivation.at("source_match_provenance_identity") ||
        lineage.at("planned_hop_provenance_identity") !=
            derivation.at("planned_hop_provenance_identity") ||
        lineage.at("derivation_provenance_identity") !=
            derivation.at("provenance_identity"))
      throw std::invalid_argument(
          "checkpoint materialized-local owner lineage is inconsistent");
    retained_derivation = std::move(derivation);
    }
  } else if (retained_owner != nullptr) {
    throw std::invalid_argument(
        "checkpoint primitive local unexpectedly acquired a derivation owner");
  }
  const auto& diagnostics = as_object(object.at("diagnostics"),
                                      "checkpoint local diagnostics");
  require_exact_keys(diagnostics,
      {"top_valid", "create_parse_ms", "create_kernel_ms"},
      "checkpoint local diagnostics");
  NativeLocalDiagnostics native{
      parse_validity(diagnostics.at("top_valid")),
      checkpoint_nonnegative_double(diagnostics.at("create_parse_ms"),
                                    "checkpoint local parse time"),
      checkpoint_nonnegative_double(diagnostics.at("create_kernel_ms"),
                                    "checkpoint local kernel time")};
  std::optional<SCCColumnProvenance> column_provenance;
  if (!object.at("column_provenance").is_null())
    column_provenance = parse_checkpoint_column_provenance(
        object.at("column_provenance"));
  std::string saved_tail_status = "unrecorded";
  bool saved_tail_attached = false;
  RegularTaylorTailModelResult restored_tail_model = unavailable_tail_model(
      "checkpoint has no serialized regular tail model");
  std::optional<TailModelCheckpointMarker> tail_checkpoint_marker;
  if (has_tail_restore) {
    const auto& tail = as_object(
        object.at("tail_model_restore"),
        "checkpoint tail-model restore marker");
    const auto* raw_serialized = tail.if_contains("serialized");
    const auto* raw_attached = tail.if_contains("attached_before_save");
    if (raw_serialized == nullptr || !raw_serialized->is_bool() ||
        raw_attached == nullptr || !raw_attached->is_bool() ||
        required_string(tail, "capability") !=
            kRegularTailMajorantCapability)
      throw std::invalid_argument(
          "checkpoint tail-model restore marker is incompatible");
    const auto serialized = raw_serialized->as_bool();
    saved_tail_attached = raw_attached->as_bool();
    saved_tail_status = required_string(tail, "status");
    if (saved_tail_status != "certified" &&
        saved_tail_status != "inconclusive" &&
        saved_tail_status != "unsupported")
      throw std::invalid_argument(
          "checkpoint tail-model restore marker has an unknown status");
    if (serialized) {
      require_exact_keys(
          tail,
          {"capability", "serialized", "status", "attached_before_save",
           "model"},
          "checkpoint serialized tail model");
      if (saved_tail_status != "certified" || !saved_tail_attached)
        throw std::invalid_argument(
            "serialized checkpoint tail model is not attached/certified");
      auto model = parse_checkpoint_regular_tail_model(tail.at("model"));
      tail_majorant_detail::validate_restored_regular_taylor_tail_model(
          model, solution, source_operator_identity);
      restored_tail_model = {
          TailMajorantStatus::Certified, std::move(model),
          "certified regular tail model restored from exact checkpoint state"};
    } else {
      require_exact_keys(
          tail,
          {"capability", "serialized", "status", "attached_before_save"},
          "checkpoint tail-model restore marker");
      restored_tail_model = unavailable_tail_model(
          "checkpoint does not serialize this regular tail-model state; "
          "saved model status was " + saved_tail_status +
          "; re-solve the retained local to reattach certification state");
      tail_checkpoint_marker = TailModelCheckpointMarker{
          saved_tail_status, saved_tail_attached};
    }
  }

  bool saved_residual_available = false;
  std::optional<json::object> saved_residual_binding;
  std::optional<std::string> saved_residual_reason;
  if (has_residual_restore) {
    const auto& residual = as_object(
        object.at("residual_operator_restore"),
        "checkpoint residual-operator restore record");
    require_exact_keys(
        residual,
        {"capability", "kind", "serialized", "status",
         "operator_payload_owner", "binding", "reason"},
        "checkpoint residual-operator restore record");
    if (required_string(residual, "capability") !=
            kOwnerBoundPhysicalResidualCapability ||
        !residual.at("serialized").is_bool())
      throw std::invalid_argument(
          "checkpoint residual-operator capability or serialization flag is malformed");
    const auto serialized = residual.at("serialized").as_bool();
    const auto status = required_string(residual, "status");
    const auto kind = required_string(residual, "kind");
    if (serialized) {
      if (status != "available" ||
          kind != "physical-homogeneous-cleared-ode" ||
          required_string(residual, "operator_payload_owner") !=
              "retained-immutable-physical-equation-owner" ||
          residual.at("binding").is_null() ||
          !residual.at("reason").is_null() ||
          !saved_equation_owner || equation_owner == nullptr ||
          physical_equation == nullptr)
        throw std::invalid_argument(
            "serialized checkpoint residual binding lost its physical q/C payload or ownership kind");
      const auto& binding = as_object(
          residual.at("binding"), "checkpoint residual binding");
      require_exact_keys(
          binding,
          {"kind", "capability", "operator_identity", "source_identity",
           "local_checkpoint_identity", "analytic_metadata",
           "owner_signature_identity", "physical_payload_identity",
           "provenance_identity", "equation", "basis", "source_kind",
           "operator_payload_owner",
           "supported_scopes"},
          "checkpoint residual binding");
      if (required_string(binding, "kind") != kind ||
          required_string(binding, "capability") !=
              kOwnerBoundPhysicalResidualCapability ||
          required_string(binding, "equation") !=
              "q(t,eps) theta(f) = C(t,eps) f" ||
          required_string(binding, "basis") !=
              "physical-original-master" ||
          required_string(binding, "source_kind") != "homogeneous-zero" ||
          required_string(binding, "operator_payload_owner") !=
              "retained-immutable-physical-equation-owner" ||
          required_string(binding, "owner_signature_identity") !=
              equation_owner->owner_signature_identity() ||
          required_string(binding, "physical_payload_identity") !=
              equation_owner->physical_payload_identity())
        throw std::invalid_argument(
            "checkpoint residual binding changes its certified equation scope");
      saved_residual_available = true;
      saved_residual_binding = binding;
    } else {
      if (status != "unsupported" || kind != "none" ||
          !residual.at("operator_payload_owner").is_null() ||
          !residual.at("binding").is_null() ||
          !residual.at("reason").is_string())
        throw std::invalid_argument(
            "unsupported checkpoint residual marker is inconsistent");
      saved_residual_reason = std::string(residual.at("reason").as_string());
      if (saved_residual_reason->empty())
        throw std::invalid_argument(
            "unsupported checkpoint residual marker has an empty reason");
    }
  }

  const auto& stats = as_object(object.at("runtime_stats"),
                                "checkpoint local runtime stats");
  if (has_tail_restore)
    require_exact_keys(stats,
        {"evaluations", "residual_certifications", "endpoint_limits",
         "line_integrations", "evaluate_ms", "residual_certify_ms",
         "endpoint_limit_ms", "line_integration_ms",
         "coefficient_count", "tail_certificate_requests",
         "tail_certificate_certified", "tail_certificate_inconclusive",
         "tail_certificate_unsupported"},
        "checkpoint local runtime stats");
  else
    require_exact_keys(stats,
        {"evaluations", "residual_certifications", "endpoint_limits",
         "line_integrations", "evaluate_ms", "residual_certify_ms",
         "endpoint_limit_ms", "line_integration_ms",
         "coefficient_count"}, "checkpoint local runtime stats");
  StoredLocalStats restored_stats;
  restored_stats.evaluations = as_u64(stats.at("evaluations"),
                                     "checkpoint local evaluations");
  restored_stats.residual_certifications = as_u64(
      stats.at("residual_certifications"),
      "checkpoint local residual certifications");
  restored_stats.endpoint_limits = as_u64(
      stats.at("endpoint_limits"), "checkpoint local endpoint limits");
  restored_stats.line_integrations = as_u64(
      stats.at("line_integrations"), "checkpoint local line integrations");
  restored_stats.evaluate_ms = checkpoint_nonnegative_double(
      stats.at("evaluate_ms"), "checkpoint local evaluation time");
  restored_stats.residual_certify_ms = checkpoint_nonnegative_double(
      stats.at("residual_certify_ms"),
      "checkpoint local residual-certification time");
  restored_stats.endpoint_limit_ms = checkpoint_nonnegative_double(
      stats.at("endpoint_limit_ms"), "checkpoint local endpoint time");
  restored_stats.line_integration_ms = checkpoint_nonnegative_double(
      stats.at("line_integration_ms"),
      "checkpoint local line-integration time");
  restored_stats.create_parse_ms = native.parse_ms;
  restored_stats.create_kernel_ms = native.kernel_ms;
  const auto coefficient_count = as_u64(
      stats.at("coefficient_count"), "checkpoint local coefficient count");
  if (coefficient_count > std::numeric_limits<std::size_t>::max())
    throw std::invalid_argument(
        "checkpoint local coefficient count exceeds size_t");
  restored_stats.coefficient_count =
      static_cast<std::size_t>(coefficient_count);
  if (has_tail_restore) {
    restored_stats.tail_certificate_requests = as_u64(
        stats.at("tail_certificate_requests"),
        "checkpoint tail certificate requests");
    restored_stats.tail_certificate_certified = as_u64(
        stats.at("tail_certificate_certified"),
        "checkpoint certified tail certificates");
    restored_stats.tail_certificate_inconclusive = as_u64(
        stats.at("tail_certificate_inconclusive"),
        "checkpoint inconclusive tail certificates");
    restored_stats.tail_certificate_unsupported = as_u64(
        stats.at("tail_certificate_unsupported"),
        "checkpoint unsupported tail certificates");
  }

  auto local = make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
      handle, source_chart, source_operator_identity, std::move(solution),
      expected_precision_bits,
      std::move(pseudo_hits), native, std::move(column_provenance),
      std::move(retained_derivation), std::move(retained_owner),
      std::move(restored_tail_model), std::move(tail_checkpoint_marker),
      has_tail_restore, has_derivation_record, std::move(equation_owner),
      std::move(physical_equation), std::string(), nullptr,
      sealed_plan_match_lineage);
  if (saved_residual_available) {
    if (!local->residual_binding().has_value() ||
        local->residual_binding()->encode() != *saved_residual_binding)
      throw std::invalid_argument(
          "checkpoint residual binding does not reproduce from its restored physical q/C owner payload");
  } else if (has_residual_restore) {
    if (local->residual_binding().has_value() ||
        !saved_residual_reason.has_value())
      throw std::invalid_argument(
          "checkpoint unsupported residual state unexpectedly reproduced an owner binding");
    local->restore_residual_binding_detail(*saved_residual_reason);
  }
  local->restore_runtime_stats(restored_stats);
  return local;
}

AcbMatchingResidualVerdict parse_checkpoint_acb_match_verdict(
    const json::value& raw, const char* label) {
  if (!raw.is_string())
    throw std::invalid_argument(std::string(label) + " must be a string");
  const auto value = std::string(raw.as_string());
  if (value == "pass") return AcbMatchingResidualVerdict::Pass;
  if (value == "fail") return AcbMatchingResidualVerdict::Fail;
  if (value == "inconclusive")
    return AcbMatchingResidualVerdict::Inconclusive;
  throw std::invalid_argument(std::string(label) +
                              " has an unsupported verdict");
}

AcbMatchingResidualDiagnostics parse_checkpoint_acb_match_residual(
    const json::value& raw, std::uint32_t dimension,
    std::int32_t expected_required_complete_max) {
  const auto& object = as_object(raw, "checkpoint Acb match residual");
  require_exact_keys(object,
      {"verdict", "complete_window", "required_complete_max",
       "complete_through_required", "coefficients", "detail"},
      "checkpoint Acb match residual");
  AcbMatchingResidualDiagnostics result;
  result.verdict = parse_checkpoint_acb_match_verdict(
      object.at("verdict"), "checkpoint Acb residual verdict");
  const auto& window = as_object(object.at("complete_window"),
                                 "checkpoint Acb residual window");
  require_exact_keys(window, {"min", "max"},
                     "checkpoint Acb residual window");
  result.complete_window = {
      as_i32(window.at("min"), "checkpoint Acb residual minimum"),
      as_i32(window.at("max"), "checkpoint Acb residual maximum")};
  const auto width = result.complete_window.width();
  result.required_complete_max = as_i32(
      object.at("required_complete_max"),
      "checkpoint Acb required complete maximum");
  if (result.required_complete_max != expected_required_complete_max)
    throw std::invalid_argument(
        "checkpoint Acb residual requirement changed across history");
  if (!object.at("complete_through_required").is_bool())
    throw std::invalid_argument(
        "checkpoint Acb residual completeness flag must be boolean");
  result.complete_through_required =
      object.at("complete_through_required").as_bool();
  if (result.complete_through_required !=
      (result.complete_window.complete_max >= result.required_complete_max))
    throw std::invalid_argument(
        "checkpoint Acb residual completeness flag contradicts its window");
  const auto& coefficients = as_array(object.at("coefficients"),
                                      "checkpoint Acb residual coefficients");
  if (coefficients.size() != checked_flat_count(
          dimension, width, "checkpoint Acb residual diagnostics"))
    throw std::invalid_argument(
        "checkpoint Acb residual coefficient history is incomplete");
  bool any_fail = false;
  bool all_pass = true;
  result.coefficients.reserve(coefficients.size());
  for (std::size_t index = 0; index < coefficients.size(); ++index) {
    const auto& coefficient = as_object(
        coefficients[index], "checkpoint Acb residual coefficient");
    require_exact_keys(coefficient,
        {"row", "epsilon_power", "residual_lower_exact",
         "residual_upper_exact", "scale_lower_exact", "scale_upper_exact",
         "verdict"}, "checkpoint Acb residual coefficient");
    const auto row = as_u64(coefficient.at("row"),
                            "checkpoint Acb residual row");
    const auto power = as_i32(coefficient.at("epsilon_power"),
                              "checkpoint Acb residual power");
    const auto expected_row = index / width;
    const auto expected_power = static_cast<std::int32_t>(
        static_cast<std::int64_t>(result.complete_window.min_power) +
        static_cast<std::int64_t>(index % width));
    if (row != expected_row || power != expected_power)
      throw std::invalid_argument(
          "checkpoint Acb residual coefficients are not a complete row-major history");
    AcbMatchingCoefficientResidual parsed;
    parsed.row = static_cast<std::size_t>(row);
    parsed.epsilon_power = power;
    parsed.residual_lower = parse_checkpoint_magnitude(
        coefficient.at("residual_lower_exact"),
        "checkpoint residual lower bound");
    parsed.residual_upper = parse_checkpoint_magnitude(
        coefficient.at("residual_upper_exact"),
        "checkpoint residual upper bound");
    parsed.scale_lower = parse_checkpoint_magnitude(
        coefficient.at("scale_lower_exact"),
        "checkpoint residual scale lower bound");
    parsed.scale_upper = parse_checkpoint_magnitude(
        coefficient.at("scale_upper_exact"),
        "checkpoint residual scale upper bound");
    parsed.verdict = parse_checkpoint_acb_match_verdict(
        coefficient.at("verdict"),
        "checkpoint Acb coefficient verdict");
    if (power <= result.required_complete_max) {
      any_fail = any_fail ||
          parsed.verdict == AcbMatchingResidualVerdict::Fail;
      all_pass = all_pass &&
          parsed.verdict == AcbMatchingResidualVerdict::Pass;
    }
    result.coefficients.push_back(std::move(parsed));
  }
  const auto derived_verdict = any_fail
      ? AcbMatchingResidualVerdict::Fail
      : (all_pass && result.complete_through_required)
          ? AcbMatchingResidualVerdict::Pass
          : AcbMatchingResidualVerdict::Inconclusive;
  if (derived_verdict != result.verdict)
    throw std::invalid_argument(
        "checkpoint Acb residual aggregate verdict is inconsistent");
  result.detail = required_string(object, "detail");
  if (result.detail.empty())
    throw std::invalid_argument(
        "checkpoint Acb residual history lost its diagnostic detail");
  return result;
}

void validate_checkpoint_exact_analytic_metadata(const json::value& raw) {
  const auto& metadata = as_object(
      raw, "checkpoint exact local analytic metadata");
  require_exact_keys(metadata,
      {"schema", "chart", "sectors", "prescriptions"},
      "checkpoint exact local analytic metadata");
  if (required_string(metadata, "schema") !=
      "diffexp2-exact-local-analytic-metadata-v2")
    throw std::invalid_argument(
        "unsupported checkpoint local analytic metadata schema");
  const auto& chart = as_object(metadata.at("chart"),
                                "checkpoint exact analytic chart");
  require_exact_keys(chart,
      {"center_exact", "scale_exact", "radius_exact_ball",
       "infinite_radius"}, "checkpoint exact analytic chart");
  (void)required_string(chart, "center_exact");
  (void)required_string(chart, "scale_exact");
  (void)parse_checkpoint_ball(chart.at("radius_exact_ball"),
                              "checkpoint exact analytic radius");
  if (!chart.at("infinite_radius").is_bool())
    throw std::invalid_argument(
        "checkpoint exact analytic radius flag must be boolean");
  const auto& sectors = as_array(metadata.at("sectors"),
                                 "checkpoint exact analytic sectors");
  if (sectors.empty())
    throw std::invalid_argument(
        "checkpoint exact analytic metadata has no sectors");
  for (const auto& raw_sector : sectors) {
    const auto& sector = as_object(raw_sector,
                                   "checkpoint exact analytic sector");
    require_exact_keys(sector, {"a", "b", "log_power"},
                       "checkpoint exact analytic sector");
    (void)parse_checkpoint_exact_descriptor(
        sector.at("a"), "checkpoint exact analytic a tag");
    (void)parse_checkpoint_exact_descriptor(
        sector.at("b"), "checkpoint exact analytic b tag");
    (void)as_u32(sector.at("log_power"),
                 "checkpoint exact analytic log power");
  }
  for (const auto& raw_prescription : as_array(
           metadata.at("prescriptions"),
           "checkpoint exact analytic prescriptions")) {
    const auto& prescription = as_object(
        raw_prescription, "checkpoint exact analytic prescription");
    require_exact_keys(prescription,
        {"factor_exact", "sign", "multiplicity",
         "leading_coefficient_sign"},
        "checkpoint exact analytic prescription");
    (void)required_string(prescription, "factor_exact");
    const auto sign = as_i32(prescription.at("sign"),
                             "checkpoint analytic prescription sign");
    const auto leading = as_i32(
        prescription.at("leading_coefficient_sign"),
        "checkpoint analytic leading sign");
    const auto multiplicity = as_u32(
        prescription.at("multiplicity"),
        "checkpoint analytic prescription multiplicity");
    if ((sign != -1 && sign != 1) ||
        (leading != -1 && leading != 1) || multiplicity == 0)
      throw std::invalid_argument(
          "checkpoint exact analytic prescription is malformed");
  }
}

std::shared_ptr<StoredEndpointResult> restore_checkpoint_endpoint_record(
    const json::value& raw, const std::string& expected_domain) {
  const auto& object = as_object(raw, "checkpoint retained endpoint");
  require_exact_keys(object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "source", "approach_direction", "requested_rim",
       "cancellation_mode", "analytic_metadata", "result", "elapsed_ms",
       "runtime_stats"}, "checkpoint retained endpoint");
  if (required_string(object, "schema") !=
      "diffexp2-retained-endpoint-result-v2")
    throw std::invalid_argument(
        "unsupported retained endpoint checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  const auto& source = as_object(object.at("source"),
                                 "checkpoint endpoint source");
  require_exact_keys(source,
      {"local", "chart", "source_operator_identity",
       "checkpoint_identity", "coefficient_domain"},
      "checkpoint endpoint source");
  const auto source_local = required_string(source, "local");
  const auto source_chart = required_string(source, "chart");
  const auto source_operator_identity = required_string(
      source, "source_operator_identity");
  const auto source_checkpoint = required_string(
      source, "checkpoint_identity");
  const auto source_domain = required_string(source, "coefficient_domain");
  if (source_domain != expected_domain ||
      (source_domain != "rational" && source_domain != "acb"))
    throw std::invalid_argument(
        "checkpoint endpoint coefficient domain is incompatible with its session");
  const auto approach_direction = as_i32(
      object.at("approach_direction"),
      "checkpoint endpoint approach direction");
  if (approach_direction != -1 && approach_direction != 1)
    throw std::invalid_argument(
        "checkpoint endpoint approach direction must be +1 or -1");
  std::optional<std::int32_t> requested_rim;
  if (!object.at("requested_rim").is_null()) {
    requested_rim = as_i32(object.at("requested_rim"),
                           "checkpoint endpoint requested rim");
    if (*requested_rim != -1 && *requested_rim != 1)
      throw std::invalid_argument(
          "checkpoint endpoint rim must be +1 or -1");
  }
  const auto cancellation_mode = required_string(
      object, "cancellation_mode");
  if (cancellation_mode != "exact-coefficient-field" &&
      cancellation_mode != "exact-or-acb-singleton")
    throw std::invalid_argument(
        "checkpoint endpoint cancellation mode is unsupported");
  auto analytic_metadata = as_object(
      object.at("analytic_metadata"),
      "checkpoint endpoint analytic metadata");
  validate_checkpoint_exact_analytic_metadata(analytic_metadata);

  json::object provenance{
      {"schema", "diffexp2-retained-native-endpoint-sector-limit-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"source", source},
      {"approach_direction", approach_direction},
      {"rim", requested_rim.has_value()
           ? json::value(*requested_rim) : json::value(nullptr)},
      {"cancellation", json::object{{"mode", cancellation_mode}}},
      {"analytic_metadata", analytic_metadata}};
  if (json::serialize(canonical_json_value(provenance)) !=
      provenance_identity)
    throw std::invalid_argument(
        "checkpoint endpoint provenance identity is inconsistent");

  const auto& raw_result = as_object(object.at("result"),
                                     "checkpoint endpoint result");
  require_exact_keys(raw_result,
      {"values", "dropped_regulated_sectors",
       "cancelled_divergent_coefficients", "imaginary_sign"},
      "checkpoint endpoint result");
  EndpointLimitResult result;
  for (const auto& raw_value : as_array(raw_result.at("values"),
                                        "checkpoint endpoint values"))
    result.values.push_back(parse_checkpoint_epsilon_frame<ComplexBall>(
        raw_value, "checkpoint endpoint value"));
  if (result.values.empty())
    throw std::invalid_argument(
        "checkpoint endpoint result has no components");
  const auto checked_size = [](const json::value& value,
                               const char* label) {
    const auto parsed = as_u64(value, label);
    if (parsed > std::numeric_limits<std::size_t>::max())
      throw std::invalid_argument(std::string(label) +
                                  " exceeds size_t");
    return static_cast<std::size_t>(parsed);
  };
  result.dropped_regulated_sectors = checked_size(
      raw_result.at("dropped_regulated_sectors"),
      "checkpoint dropped regulated sectors");
  result.cancelled_divergent_coefficients = checked_size(
      raw_result.at("cancelled_divergent_coefficients"),
      "checkpoint cancelled divergent coefficients");
  result.imaginary_sign = as_i32(raw_result.at("imaginary_sign"),
                                 "checkpoint endpoint effective rim");
  if (result.imaginary_sign != -1 && result.imaginary_sign != 1)
    throw std::invalid_argument(
        "checkpoint endpoint effective rim must be +1 or -1");
  (void)endpoint_value_window(result);
  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint endpoint elapsed time");
  const auto& stats = as_object(object.at("runtime_stats"),
                                "checkpoint endpoint runtime stats");
  require_exact_keys(stats, {"exports", "export_ms"},
                     "checkpoint endpoint runtime stats");
  const auto exports = as_u64(stats.at("exports"),
                              "checkpoint endpoint exports");
  const auto export_ms = checkpoint_nonnegative_double(
      stats.at("export_ms"), "checkpoint endpoint export time");
  auto endpoint = std::make_shared<StoredEndpointResult>(
      handle, checkpoint_identity, provenance_identity, source_local,
      source_chart, source_operator_identity, source_checkpoint,
      source_domain, approach_direction, requested_rim, cancellation_mode,
      std::move(analytic_metadata), std::move(result), elapsed_ms);
  endpoint->restore_runtime_stats(exports, export_ms);
  return endpoint;
}

void validate_checkpoint_match_source(const json::object& source,
                                      bool basis,
                                      std::size_t expected_column = 0) {
  const bool has_column_provenance =
      source.if_contains("column_provenance") != nullptr;
  if (basis) {
    if (has_column_provenance)
      require_exact_keys(source,
          {"column", "local", "chart", "source_operator_identity",
           "checkpoint_identity",
           "requested_imaginary_sign", "effective_imaginary_sign",
           "analytic_metadata", "column_provenance"},
          "checkpoint Acb basis source");
    else
      require_exact_keys(source,
          {"column", "local", "chart", "source_operator_identity",
           "checkpoint_identity",
           "requested_imaginary_sign", "effective_imaginary_sign",
           "analytic_metadata"}, "checkpoint Acb basis source");
    if (as_u64(source.at("column"), "checkpoint Acb basis column") !=
        expected_column)
      throw std::invalid_argument(
          "checkpoint Acb basis sources are not in column order");
  } else {
    if (has_column_provenance)
      require_exact_keys(source,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity",
           "requested_imaginary_sign", "effective_imaginary_sign",
           "analytic_metadata", "column_provenance"},
          "checkpoint Acb incoming source");
    else
      require_exact_keys(source,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity",
           "requested_imaginary_sign", "effective_imaginary_sign",
           "analytic_metadata"}, "checkpoint Acb incoming source");
  }
  if (required_string(source, "local").empty() ||
      required_string(source, "chart").empty() ||
      required_string(source, "source_operator_identity").empty() ||
      required_string(source, "checkpoint_identity").empty())
    throw std::invalid_argument(
        "checkpoint Acb match source has empty provenance");
  for (const auto* key : {"requested_imaginary_sign",
                          "effective_imaginary_sign"}) {
    if (source.at(key).is_null()) continue;
    const auto sign = as_i32(source.at(key), key);
    if (sign != -1 && sign != 1)
      throw std::invalid_argument(
          "checkpoint Acb match source has an invalid branch sign");
  }
  validate_checkpoint_exact_analytic_metadata(
      source.at("analytic_metadata"));
  if (has_column_provenance)
    (void)parse_checkpoint_column_provenance(
        source.at("column_provenance"));
}

void validate_checkpoint_exact_match_source(
    const json::object& source, bool basis_source,
    std::size_t expected_column,
    const std::shared_ptr<StoredLocalBase>& owner) {
  const bool has_column_provenance =
      source.if_contains("column_provenance") != nullptr;
  if (basis_source) {
    if (has_column_provenance)
      require_exact_keys(
          source,
          {"column", "local", "chart", "source_operator_identity",
           "checkpoint_identity", "analytic_metadata",
           "column_provenance"},
          "checkpoint exact-match basis source");
    else
      require_exact_keys(
          source,
          {"column", "local", "chart", "source_operator_identity",
           "checkpoint_identity", "analytic_metadata"},
          "checkpoint exact-match basis source");
    if (as_u64(source.at("column"),
               "checkpoint exact-match basis column") != expected_column)
      throw std::invalid_argument(
          "checkpoint exact-match basis sources are not in column order");
  } else {
    if (has_column_provenance)
      require_exact_keys(
          source,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity", "analytic_metadata",
           "column_provenance"},
          "checkpoint exact-match incoming source");
    else
      require_exact_keys(
          source,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity", "analytic_metadata"},
          "checkpoint exact-match incoming source");
  }
  if (!owner || std::string(owner->scalar_domain()) != "rational")
    throw std::invalid_argument(
        "checkpoint exact-match source has no exact-rational owner");
  if (required_string(source, "local") != owner->handle() ||
      required_string(source, "chart") != owner->source_chart() ||
      required_string(source, "source_operator_identity") !=
          owner->source_operator_identity() ||
      required_string(source, "checkpoint_identity") !=
          owner->checkpoint_identity() ||
      source.at("analytic_metadata") != owner->exact_analytic_metadata())
    throw std::invalid_argument(
        "checkpoint exact-match source provenance disagrees with its retained local owner");
  validate_checkpoint_exact_analytic_metadata(
      source.at("analytic_metadata"));
  if (has_column_provenance)
    (void)parse_checkpoint_column_provenance(
        source.at("column_provenance"));
}

std::shared_ptr<StoredExactRegularMatch>
restore_checkpoint_exact_match_record(
    const json::value& raw,
    std::vector<std::shared_ptr<StoredLocalBase>> basis_owners,
    std::shared_ptr<StoredLocalBase> incoming_owner) {
  const auto& object = as_object(raw,
                                 "checkpoint retained exact-rational match");
  require_exact_keys(
      object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "basis_sources", "incoming_source", "basis_chart",
       "incoming_chart", "basis_point_exact", "incoming_point_exact",
       "physical_match_point_exact", "epsilon", "dimension",
       "transformation", "weights", "saturation", "residual_window",
       "elapsed_ms"},
      "checkpoint retained exact-rational match");
  if (required_string(object, "schema") !=
      "diffexp2-retained-exact-rational-match-v2")
    throw std::invalid_argument(
        "unsupported retained exact-rational-match checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  const auto basis_chart = required_string(object, "basis_chart");
  const auto incoming_chart = required_string(object, "incoming_chart");
  const auto basis_point = required_string(object, "basis_point_exact");
  const auto incoming_point = required_string(object,
                                               "incoming_point_exact");
  const auto physical_point = required_string(
      object, "physical_match_point_exact");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty() || basis_chart.empty() ||
      incoming_chart.empty() || basis_point.empty() ||
      incoming_point.empty() || physical_point.empty())
    throw std::invalid_argument(
        "checkpoint exact-rational match contains an empty identity");
  (void)Rational(basis_point);
  (void)Rational(incoming_point);
  (void)Rational(physical_point);
  const auto dimension = as_u32(object.at("dimension"),
                                "checkpoint exact-match dimension");
  if (dimension == 0 || basis_owners.size() != dimension ||
      !incoming_owner)
    throw std::invalid_argument(
        "checkpoint exact-match ownership differs from its dimension");
  std::vector<json::object> basis_sources;
  const auto& raw_basis = as_array(object.at("basis_sources"),
                                   "checkpoint exact-match basis sources");
  if (raw_basis.size() != dimension)
    throw std::invalid_argument(
        "checkpoint exact-match source count differs from its dimension");
  basis_sources.reserve(dimension);
  for (std::size_t column = 0; column < dimension; ++column) {
    auto source = as_object(raw_basis[column],
                            "checkpoint exact-match basis source");
    validate_checkpoint_exact_match_source(
        source, true, column, basis_owners[column]);
    if (required_string(source, "chart") != basis_chart)
      throw std::invalid_argument(
          "checkpoint exact-match basis chart provenance is inconsistent");
    basis_sources.push_back(std::move(source));
  }
  auto incoming_source = as_object(
      object.at("incoming_source"),
      "checkpoint exact-match incoming source");
  validate_checkpoint_exact_match_source(
      incoming_source, false, 0, incoming_owner);
  if (required_string(incoming_source, "chart") != incoming_chart)
    throw std::invalid_argument(
        "checkpoint exact-match incoming chart provenance is inconsistent");

  const auto& raw_epsilon = as_object(object.at("epsilon"),
                                      "checkpoint exact-match epsilon");
  require_exact_keys(raw_epsilon, {"min", "max", "required_complete_max"},
                     "checkpoint exact-match epsilon");
  EpsilonWindow window{
      as_i32(raw_epsilon.at("min"), "checkpoint exact-match epsilon minimum"),
      as_i32(raw_epsilon.at("max"), "checkpoint exact-match epsilon maximum")};
  (void)window.width();
  const auto required_complete_max = as_i32(
      raw_epsilon.at("required_complete_max"),
      "checkpoint exact-match required complete maximum");
  if (required_complete_max < window.min_power ||
      required_complete_max > window.complete_max)
    throw std::invalid_argument(
        "checkpoint exact-match requirement lies outside its work window");
  auto transformation = parse_checkpoint_exact_laurent_matrix(
      object.at("transformation"), dimension,
      "checkpoint exact-match transformation");
  auto weights = parse_checkpoint_frame_vector<Rational>(
      object.at("weights"), dimension, "checkpoint exact-match weights");
  auto diagnostics = parse_checkpoint_saturation_diagnostics(
      object.at("saturation"), dimension);
  const auto& raw_residual = as_object(
      object.at("residual_window"), "checkpoint exact-match residual window");
  require_exact_keys(raw_residual, {"min", "max"},
                     "checkpoint exact-match residual window");
  EpsilonWindow residual{
      as_i32(raw_residual.at("min"), "checkpoint residual minimum"),
      as_i32(raw_residual.at("max"), "checkpoint residual maximum")};
  (void)residual.width();
  if (residual.complete_max < required_complete_max)
    throw std::invalid_argument(
        "checkpoint exact-match residual lost its required complete window");

  json::array provenance_basis;
  for (const auto& source : basis_sources)
    provenance_basis.push_back(source);
  json::object provenance{
      {"schema", "diffexp2-native-exact-regular-local-match-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"basis", std::move(provenance_basis)},
      {"incoming", incoming_source},
      {"basis_point_exact", basis_point},
      {"incoming_point_exact", incoming_point},
      {"physical_match_point_exact", physical_point},
      {"epsilon", json::object{{"min", window.min_power},
                                {"max", window.complete_max},
                                {"required_complete_max",
                                 required_complete_max}}}};
  if (json::serialize(canonical_json_value(provenance)) !=
      provenance_identity)
    throw std::invalid_argument(
        "checkpoint exact-match provenance identity is inconsistent");
  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint exact-match elapsed time");
  return std::make_shared<StoredExactRegularMatch>(
      handle, checkpoint_identity, provenance_identity,
      std::move(basis_sources), std::move(incoming_source), basis_chart,
      incoming_chart, basis_point, incoming_point, physical_point, window,
      required_complete_max, dimension, std::move(transformation),
      std::move(weights), std::move(diagnostics), residual, elapsed_ms,
      std::move(basis_owners), std::move(incoming_owner));
}

std::shared_ptr<StoredRefinedAcbMatch> restore_checkpoint_acb_match_record(
    const json::value& raw,
    const std::optional<std::string>& expected_session_configuration =
        std::nullopt,
    const std::optional<json::object>& expected_singular_request =
        std::nullopt) {
  const auto& object = as_object(raw, "checkpoint retained Acb match");
  require_exact_keys(object,
      {"schema", "handle", "checkpoint_identity", "provenance_identity",
       "exact_lattice_identity", "exact_lattice_provenance_identity",
       "exact_lattice_canonical_witness", "basis_sources",
       "incoming_source", "basis_chart", "incoming_chart",
       "basis_point_exact", "incoming_point_exact",
       "physical_match_point_exact", "matching_frame_identity",
       "epsilon", "dimension",
       "relative_tolerance", "max_refinement_steps", "refined",
       "elapsed_ms"}, "checkpoint retained Acb match");
  if (required_string(object, "schema") !=
      "diffexp2-retained-acb-match-v2")
    throw std::invalid_argument(
        "unsupported retained Acb-match checkpoint schema");
  const auto handle = required_string(object, "handle");
  const auto checkpoint_identity = required_string(
      object, "checkpoint_identity");
  const auto provenance_identity = required_string(
      object, "provenance_identity");
  const auto exact_lattice_identity = required_string(
      object, "exact_lattice_identity");
  const auto exact_lattice_provenance_identity = required_string(
      object, "exact_lattice_provenance_identity");
  if (handle.empty() || checkpoint_identity.empty() ||
      provenance_identity.empty() || exact_lattice_identity.empty() ||
      exact_lattice_provenance_identity.empty())
    throw std::invalid_argument(
        "checkpoint retained Acb match contains an empty identity");
  const auto basis_chart = required_string(object, "basis_chart");
  const auto incoming_chart = required_string(object, "incoming_chart");
  const auto basis_point = required_string(object, "basis_point_exact");
  const auto incoming_point = required_string(object,
                                               "incoming_point_exact");
  const auto physical_point = required_string(
      object, "physical_match_point_exact");
  const auto matching_frame_identity = required_string(
      object, "matching_frame_identity");
  if (basis_chart.empty() || incoming_chart.empty() || basis_point.empty() ||
      incoming_point.empty() || physical_point.empty() ||
      matching_frame_identity.empty())
    throw std::invalid_argument(
        "checkpoint Acb match lost chart or point provenance");
  const auto& raw_epsilon = as_object(object.at("epsilon"),
                                      "checkpoint Acb match epsilon window");
  require_exact_keys(raw_epsilon,
      {"min", "max", "required_complete_max"},
      "checkpoint Acb match epsilon window");
  const EpsilonWindow window{
      as_i32(raw_epsilon.at("min"), "checkpoint match epsilon minimum"),
      as_i32(raw_epsilon.at("max"), "checkpoint match epsilon maximum")};
  (void)window.width();
  const auto required_complete_max = as_i32(
      raw_epsilon.at("required_complete_max"),
      "checkpoint match required complete maximum");
  if (required_complete_max < window.min_power ||
      required_complete_max > window.complete_max)
    throw std::invalid_argument(
        "checkpoint Acb match requirement lies outside its work window");
  const auto dimension = as_u32(object.at("dimension"),
                                "checkpoint Acb match dimension");
  if (dimension == 0)
    throw std::invalid_argument("checkpoint Acb match dimension is zero");
  const auto relative_tolerance = required_string(
      object, "relative_tolerance");
  auto parsed_tolerance = Magnitude::decimal(relative_tolerance);
  if (!parsed_tolerance.is_finite())
    throw std::invalid_argument(
        "checkpoint Acb match tolerance is not finite");
  const auto max_refinement_steps = as_u64(
      object.at("max_refinement_steps"),
      "checkpoint Acb maximum refinement steps");
  if (max_refinement_steps > 32)
    throw std::invalid_argument(
        "checkpoint Acb maximum refinement steps exceeds 32");

  std::vector<json::object> basis_sources;
  const auto& raw_basis_sources = as_array(
      object.at("basis_sources"), "checkpoint Acb basis sources");
  if (raw_basis_sources.size() != dimension)
    throw std::invalid_argument(
        "checkpoint Acb match basis source count differs from its dimension");
  basis_sources.reserve(raw_basis_sources.size());
  for (std::size_t column = 0; column < raw_basis_sources.size(); ++column) {
    auto source = as_object(raw_basis_sources[column],
                            "checkpoint Acb basis source");
    validate_checkpoint_match_source(source, true, column);
    basis_sources.push_back(std::move(source));
  }
  auto incoming_source = as_object(object.at("incoming_source"),
                                   "checkpoint Acb incoming source");
  validate_checkpoint_match_source(incoming_source, false);

  const auto exact_lattice_witness_record = required_string(
      object, "exact_lattice_canonical_witness");
  if (exact_lattice_witness_record.empty())
    throw std::invalid_argument(
        "checkpoint Acb match lost its exact lattice witness");
  const auto raw_saturation_witness = json::parse(
      exact_lattice_witness_record);
  const auto& saturation_witness = as_object(
      raw_saturation_witness, "checkpoint Acb saturation witness");
  const auto saturation_witness_schema = required_string(
      saturation_witness, "schema");
  auto exact_lattice = [&]() -> ParsedExactEvaluatedLattice {
    if (saturation_witness_schema == kExactEvaluatedLatticeSchema)
      return parse_exact_evaluated_lattice(
          raw_saturation_witness, dimension, window,
          checkpoint_identity + ":checkpoint-restore");
    if (saturation_witness_schema == kNativeUnitSaturationProofSchema)
      return parse_native_unit_saturation_proof(
          raw_saturation_witness, dimension, window, basis_sources,
          basis_point, physical_point,
          checkpoint_identity + ":checkpoint-restore");
    if (saturation_witness_schema ==
        kNativeSingularSCCSaturationProofSchema) {
      if (!expected_session_configuration.has_value() ||
          !expected_singular_request.has_value())
        throw std::invalid_argument(
            "checkpoint singular-SCC saturation witness lacks its retained planned-match binding");
      return parse_native_singular_scc_saturation_proof(
          raw_saturation_witness, dimension, window, basis_sources,
          basis_point, physical_point, expected_session_configuration,
          expected_singular_request,
          checkpoint_identity + ":checkpoint-restore");
    }
    throw std::invalid_argument(
        "checkpoint Acb match has an unsupported saturation witness schema");
  }();
  if (exact_lattice.identity != exact_lattice_identity ||
      exact_lattice.canonical_witness != exact_lattice_witness_record)
    throw std::invalid_argument(
        "checkpoint exact lattice identity or canonical witness is inconsistent");

  json::array exact_binding_basis;
  for (const auto& source : basis_sources)
    exact_binding_basis.push_back(source);
  json::object exact_lattice_provenance{
      {"schema", "diffexp2-retained-exact-lattice-binding-v1"},
      {"witness_schema", exact_lattice.witness_schema},
      {"witness_identity", exact_lattice_identity},
      {"basis", std::move(exact_binding_basis)},
      {"basis_point_exact", basis_point},
      {"physical_match_point_exact", physical_point},
      {"matching_frame_identity", matching_frame_identity},
      {"epsilon", json::object{{"min", window.min_power},
                                {"max", window.complete_max}}}};
  if (json::serialize(canonical_json_value(exact_lattice_provenance)) !=
      exact_lattice_provenance_identity)
    throw std::invalid_argument(
        "checkpoint exact lattice provenance identity is inconsistent");

  json::array provenance_basis;
  for (const auto& source : basis_sources)
    provenance_basis.push_back(source);
  json::object provenance{
      {"schema", "diffexp2-native-refined-acb-local-match-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"basis", std::move(provenance_basis)},
      {"incoming", incoming_source},
      {"basis_point_exact", basis_point},
      {"incoming_point_exact", incoming_point},
      {"physical_match_point_exact", physical_point},
      {"matching_frame_identity", matching_frame_identity},
      {"epsilon", json::object{{"min", window.min_power},
                                {"max", window.complete_max},
                                {"required_complete_max",
                                 required_complete_max}}},
      {"exact_lattice_provenance_identity",
       exact_lattice_provenance_identity},
      {"refinement", json::object{{"relative_tolerance",
                                    relative_tolerance},
                                   {"max_steps",
                                    max_refinement_steps}}}};
  if (json::serialize(canonical_json_value(provenance)) !=
      provenance_identity)
    throw std::invalid_argument(
        "checkpoint retained Acb match provenance identity is inconsistent");

  const auto& raw_refined = as_object(object.at("refined"),
                                      "checkpoint refined Acb state");
  require_exact_keys(raw_refined,
      {"transformed_weights", "weights", "residual", "residual_history",
       "refinement_steps", "factorization_preconditioner"},
      "checkpoint refined Acb state");
  RefinedAcbLaurentMatch refined;
  refined.transformed_weights = parse_checkpoint_frame_vector<ComplexBall>(
      raw_refined.at("transformed_weights"), dimension,
      "checkpoint transformed Acb weights");
  refined.weights = parse_checkpoint_frame_vector<ComplexBall>(
      raw_refined.at("weights"), dimension, "checkpoint Acb weights");
  refined.residual = parse_checkpoint_frame_vector<ComplexBall>(
      raw_refined.at("residual"), dimension, "checkpoint Acb residual");
  refined.refinement_steps = static_cast<std::size_t>(as_u64(
      raw_refined.at("refinement_steps"),
      "checkpoint Acb refinement steps"));
  refined.factorization_preconditioner = required_string(
      raw_refined, "factorization_preconditioner");
  if (refined.refinement_steps > max_refinement_steps)
    throw std::invalid_argument(
        "checkpoint Acb refinement count exceeds its policy");
  for (const auto& raw_history : as_array(
           raw_refined.at("residual_history"),
           "checkpoint Acb residual history"))
    refined.residual_history.push_back(parse_checkpoint_acb_match_residual(
        raw_history, dimension, required_complete_max));
  if (refined.residual_history.size() != refined.refinement_steps + 1)
    throw std::invalid_argument(
        "checkpoint Acb residual history does not cover every refinement step");
  auto residual_min = refined.residual.front().min_power();
  auto residual_max = refined.residual.front().complete_max();
  for (const auto& row : refined.residual) {
    residual_min = std::min(residual_min, row.min_power());
    residual_max = std::min(residual_max, row.complete_max());
  }
  const auto& final_diagnostics = refined.residual_history.back();
  if (final_diagnostics.complete_window.min_power < residual_min ||
      final_diagnostics.complete_window.complete_max > residual_max ||
      final_diagnostics.complete_window.min_power > required_complete_max ||
      final_diagnostics.complete_window.complete_max < required_complete_max)
    throw std::invalid_argument(
        "checkpoint final Acb residual diagnostics do not match the retained residual frame");
  const auto elapsed_ms = checkpoint_nonnegative_double(
      object.at("elapsed_ms"), "checkpoint Acb match elapsed time");
  return std::make_shared<StoredRefinedAcbMatch>(
      handle, checkpoint_identity, provenance_identity,
      exact_lattice_identity, exact_lattice_provenance_identity,
      exact_lattice_witness_record, exact_lattice.witness_schema,
      std::move(basis_sources),
      std::move(incoming_source), basis_chart, incoming_chart, basis_point,
      incoming_point, physical_point, matching_frame_identity, window,
      required_complete_max,
      dimension, relative_tolerance,
      static_cast<std::size_t>(max_refinement_steps),
      std::move(exact_lattice.saturation), std::move(refined), elapsed_ms);
}
