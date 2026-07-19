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

struct CheckpointValueHandoffPlanBinding {
  json::object producing_chart;
  json::object producing_owner;
  json::object receiving_chart;
  json::object receiving_owner;
  std::shared_ptr<StoredLocalBase> incoming;
  std::optional<std::int32_t> producing_rim;
};

template <typename Scalar>
std::shared_ptr<StoredLocalBase> restore_checkpoint_local_record(
    const json::value& raw, const std::string& expected_domain,
    slong expected_precision_bits,
    std::shared_ptr<void> retained_owner = nullptr,
    std::shared_ptr<PhysicalEquationOwnerBase> equation_owner = nullptr,
    const CheckpointValueHandoffPlanBinding* value_handoff_plan = nullptr) {
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
  const auto sealed_owner_restore = sealed_plan_match_lineage
      ? required_string(object, "derivation_owner_restore")
      : std::string();
  if (sealed_plan_match_lineage &&
      sealed_owner_restore != "sealed-plan-match-lineage" &&
      sealed_owner_restore != "sealed-plan-value-handoff-lineage")
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
            "diffexp2-retained-plan-match-local-materialization-v2" &&
        !is_retained_plan_value_handoff_schema(derivation_schema))
      throw std::invalid_argument(
          "checkpoint sealed lineage is not a supported transport materialization");
    if (sealed_plan_match_lineage &&
        (is_retained_plan_value_handoff_schema(derivation_schema) !=
         (sealed_owner_restore ==
              "sealed-plan-value-handoff-lineage")))
      throw std::invalid_argument(
          "checkpoint sealed transport derivation and restore policy disagree");
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
    } else if (is_retained_plan_value_handoff_schema(
                   derivation_schema)) {
      const bool version_two_handoff = derivation_schema ==
          "diffexp2-retained-plan-value-handoff-v2";
      require_exact_keys(
          derivation,
          {"schema", "capability", "tile_plan",
           "tile_plan_checkpoint_identity", "tile_plan_provenance_identity",
           "arm", "match", "producing", "receiving",
           "receiver_center_physical_exact", "producing_local_exact",
           "prototype_identity", "tail_contract", "accuracy_contract",
           "epsilon", "incoming", "scope", "coefficient_transport",
           "whole_arm_complete", "evaluated_epsilon", "output",
           "equation_owner_signature_identity",
           "equation_payload_identity", "provenance_identity"},
          "checkpoint plan-value handoff derivation");
      if (required_string(derivation, "capability") !=
              (version_two_handoff
                   ? "retained-native-regular-value-handoff-v2"
                   : "retained-native-regular-value-handoff-v1") ||
          required_string(derivation, "scope") !=
              "single-regular-to-regular-transport-hop" ||
          required_string(derivation, "coefficient_transport") !=
              "native-retained-only" ||
          !derivation.at("whole_arm_complete").is_bool() ||
          derivation.at("whole_arm_complete").as_bool() ||
          !equation_owner ||
          (std::string(equation_owner->equation_owner_kind()) !=
               "prepared-chart" &&
           std::string(equation_owner->equation_owner_kind()) !=
               "regular-physical-equation-v1") ||
          required_string(derivation,
                          "equation_owner_signature_identity") !=
              equation_owner->owner_signature_identity() ||
          required_string(derivation, "equation_payload_identity") !=
              equation_owner->physical_payload_identity())
        throw std::invalid_argument(
            "checkpoint plan-value handoff changed its sealed scope or equation owner");
      const auto arm = required_string(derivation, "arm");
      if (arm != "lower" && arm != "upper")
        throw std::invalid_argument(
            "checkpoint plan-value handoff arm is invalid");
      (void)as_u64(derivation.at("match"),
                   "checkpoint plan-value handoff match index");

      const auto& producing = as_object(
          derivation.at("producing"),
          "checkpoint value producing binding");
      const auto& receiving = as_object(
          derivation.at("receiving"),
          "checkpoint value receiving binding");
      require_exact_keys(producing, {"chart", "owner"},
                         "checkpoint value producing binding");
      require_exact_keys(receiving, {"chart", "owner"},
                         "checkpoint value receiving binding");
      if (value_handoff_plan == nullptr ||
          producing.at("chart") != value_handoff_plan->producing_chart ||
          producing.at("owner") != value_handoff_plan->producing_owner ||
          receiving.at("chart") != value_handoff_plan->receiving_chart ||
          receiving.at("owner") != value_handoff_plan->receiving_owner)
        throw std::invalid_argument(
            "checkpoint value chart/owner bindings differ from their restored tile plan");
      const auto& producing_chart = as_object(
          producing.at("chart"), "checkpoint value producing chart");
      const auto& receiving_chart = as_object(
          receiving.at("chart"), "checkpoint value receiving chart");
      const auto& receiving_owner = as_object(
          receiving.at("owner"), "checkpoint value receiving owner");
      if (required_string(receiving_chart, "chart") != source_chart ||
          required_string(receiving_chart, "identity") !=
              source_operator_identity ||
          required_string(receiving_owner, "kind") !=
              equation_owner->equation_owner_kind() ||
          required_string(receiving_owner, "handle") != source_chart ||
          required_string(receiving_owner, "operator_identity") !=
              source_operator_identity ||
          required_string(receiving_owner, "owner_signature_identity") !=
              equation_owner->owner_signature_identity() ||
          required_string(receiving_owner, "physical_payload_identity") !=
              equation_owner->physical_payload_identity())
        throw std::invalid_argument(
            "checkpoint value receiving owner differs from its restored physical chart");
      const Rational producing_center(
          required_string(producing_chart, "center_exact"));
      const Rational producing_scale(
          required_string(producing_chart, "scale_exact"));
      const Rational producing_radius(
          required_string(producing_chart, "radius_exact"));
      const Rational receiving_center(
          required_string(receiving_chart, "center_exact"));
      if (producing_scale.is_zero() || producing_radius.sign() <= 0)
        throw std::invalid_argument(
            "checkpoint value producing chart has invalid geometry");
      const auto producing_local =
          (receiving_center - producing_center) / producing_scale;
      const auto center_ratio =
          exact_path_detail::abs(producing_local) / producing_radius;
      if (Rational(required_string(
              derivation, "receiver_center_physical_exact")) !=
              receiving_center ||
          Rational(required_string(derivation, "producing_local_exact")) !=
              producing_local)
        throw std::invalid_argument(
            "checkpoint plan-value handoff changed its exact center geometry");

      const auto prototype_identity = required_string(
          derivation, "prototype_identity");
      const auto prototype_value = json::parse(prototype_identity);
      const auto& prototype = as_object(
          prototype_value, "checkpoint value-solver prototype identity");
      const auto prototype_schema =
          required_string(prototype, "schema");
      const bool physical_value_prototype =
          prototype_schema ==
          "diffexp2-native-ordinary-physical-value-solver-v1";
      if (physical_value_prototype)
        require_exact_keys(
            prototype,
            {"schema", "taylor_complete_max", "metadata",
             "relative_accuracy_max_exact"},
            "checkpoint physical value-solver prototype identity");
      else
        require_exact_keys(
            prototype,
            {"schema", "run", "metadata", "tail_proxy_max_exact",
             "relative_accuracy_max_exact"},
            "checkpoint value-solver prototype identity");
      if ((!physical_value_prototype &&
           prototype_schema !=
               "diffexp2-native-regular-value-solver-prototype-v1") ||
          json::serialize(canonical_json_value(prototype_value)) !=
              prototype_identity)
        throw std::invalid_argument(
            "checkpoint value-solver prototype identity is not canonical");
      const auto prototype_taylor_complete_max =
          physical_value_prototype
          ? as_u32(prototype.at("taylor_complete_max"),
                   "checkpoint physical value-solver Taylor maximum")
          : as_u32(
                as_object(prototype.at("run"),
                          "checkpoint value-solver run")
                    .at("nmax"),
                "checkpoint value-solver Taylor maximum");
      const auto& tail = as_object(
          derivation.at("tail_contract"),
          "checkpoint value tail contract");
      const auto tail_mode = tail.if_contains("mode") != nullptr
          ? required_string(tail, "mode") : std::string();
      const bool physical_evolution_handoff =
          version_two_handoff &&
          tail_mode == "physical-evolution-tail-unavailable-v1";
      const bool certified_tail_handoff =
          version_two_handoff &&
          tail_mode == "certified-regular-taylor-point-tail-acb-v1";
      if (version_two_handoff && !physical_evolution_handoff &&
          !certified_tail_handoff)
        throw std::invalid_argument(
            "checkpoint version-two value handoff has an unsupported tail mode");
      if (physical_value_prototype != physical_evolution_handoff)
        throw std::invalid_argument(
            "checkpoint value-solver prototype and retained tail mode disagree");
      if (certified_tail_handoff)
        require_exact_keys(
            tail,
            {"mode", "producer_taylor_complete_max",
             "receiver_taylor_complete_max", "center_ratio_exact",
             "producing_point_exact", "producing_chart_radius_exact",
             "witness_radius_exact", "witness_dyadic_inward_exponent",
             "source_model", "certificate", "inflation"},
            "checkpoint certified value tail contract");
      else if (physical_evolution_handoff)
        require_exact_keys(
            tail,
            {"mode", "source_model", "witness_radius_exact",
             "witness_dyadic_inward_exponent", "output_status",
             "output_reason"},
            "checkpoint physical-evolution value tail contract");
      else
        require_exact_keys(
            tail,
            {"producer_taylor_complete_max", "receiver_taylor_complete_max",
             "center_ratio_exact", "tail_proxy_exact",
             "tail_proxy_max_exact"},
            "checkpoint value tail contract");
      const auto producer_order = physical_evolution_handoff
          ? as_u32(as_object(
                       derivation.at("incoming"),
                       "checkpoint physical-evolution incoming")
                       .at("taylor_complete_max"),
                   "checkpoint value producer Taylor maximum")
          : as_u32(tail.at("producer_taylor_complete_max"),
                   "checkpoint value producer Taylor maximum");
      const auto receiver_order = physical_evolution_handoff
          ? prototype_taylor_complete_max
          : as_u32(tail.at("receiver_taylor_complete_max"),
                   "checkpoint value receiver Taylor maximum");
      if (prototype_taylor_complete_max != receiver_order ||
          receiver_order != solution.taylor_complete_max)
        throw std::invalid_argument(
            "checkpoint value receiver order differs from its prototype");
      const auto prepared_owner =
          std::dynamic_pointer_cast<PreparedChartBase>(equation_owner);
      const auto regular_equation_owner =
          std::dynamic_pointer_cast<RegularPhysicalEquationOwnerBase>(
              equation_owner);
      if ((!physical_evolution_handoff &&
           (!prepared_owner ||
            Rational(required_string(tail, "center_ratio_exact")) !=
                center_ratio ||
            required_string(prototype, "tail_proxy_max_exact") !=
                prepared_owner->regular_value_tail_proxy_max_exact())) ||
          (physical_evolution_handoff && !regular_equation_owner))
        throw std::invalid_argument(
            "checkpoint plan-value handoff tail contract is inconsistent");
      std::optional<EpsilonVector> certified_inflated_value;
      if (!version_two_handoff) {
        Rational tail_proxy(1);
        for (std::uint64_t power = 0;
             power < static_cast<std::uint64_t>(producer_order) + 1;
             ++power)
          tail_proxy *= center_ratio;
        const Rational tail_proxy_max(
            required_string(tail, "tail_proxy_max_exact"));
        if (Rational(required_string(tail, "tail_proxy_exact")) !=
                tail_proxy ||
            required_string(tail, "tail_proxy_max_exact") !=
                required_string(prototype, "tail_proxy_max_exact") ||
            required_string(tail, "tail_proxy_max_exact") !=
                prepared_owner->regular_value_tail_proxy_max_exact() ||
            !(tail_proxy < tail_proxy_max))
          throw std::invalid_argument(
              "checkpoint plan-value handoff proxy tail contract is inconsistent");
      } else if (certified_tail_handoff) {
        if (expected_domain != "acb" ||
            required_string(tail, "mode") !=
                "certified-regular-taylor-point-tail-acb-v1" ||
            Rational(required_string(tail, "producing_point_exact")) !=
                producing_local ||
            Rational(required_string(
                tail, "producing_chart_radius_exact")) !=
                producing_radius)
          throw std::invalid_argument(
              "checkpoint certified value handoff changed its Acb geometry");
        const Rational witness(
            required_string(tail, "witness_radius_exact"));
        const auto witness_exponent = as_u32(
            tail.at("witness_dyadic_inward_exponent"),
            "checkpoint certified value witness exponent");
        Rational witness_denominator(1);
        for (std::uint32_t exponent = 0; exponent < witness_exponent;
             ++exponent)
          witness_denominator *= Rational(2);
        const auto point_modulus = exact_path_detail::abs(producing_local);
        const auto expected_witness = point_modulus +
            (producing_radius - point_modulus) / witness_denominator;
        if (witness_exponent == 0 || witness_exponent > 16 ||
            witness != expected_witness ||
            !(point_modulus < witness) || !(witness < producing_radius))
          throw std::invalid_argument(
              "checkpoint certified value witness is not strictly interior");
        auto model = parse_checkpoint_regular_tail_model(
            tail.at("source_model"));
        json::array model_prescriptions;
        for (const auto& prescription : model.prescriptions)
          model_prescriptions.push_back(json::object{
              {"factor_exact", prescription.factor_exact},
              {"sign", prescription.sign},
              {"multiplicity", prescription.multiplicity},
              {"leading_coefficient_sign",
               prescription.leading_coefficient_sign}});
        if (model.q_coefficients.empty() ||
            model.n_coefficients.empty() ||
            model.n_row_sum_upper.size() != model.n_coefficients.size() ||
            model.initial_row_upper.size() != model.epsilon.width() ||
            model.dimension != solution.dimension ||
            model.taylor_complete_max != producer_order ||
            model.operator_identity !=
                required_string(producing_chart, "identity") ||
            model.local_checkpoint_identity != required_string(
                as_object(derivation.at("incoming"),
                          "checkpoint value incoming"),
                "checkpoint_identity") ||
            Rational(model.chart.center_exact) != producing_center ||
            Rational(model.chart.scale_exact) != producing_scale ||
            model.chart.infinite_radius ||
            !acb_equal(model.chart.radius.raw(),
                       ComplexBall::from_strings(
                           producing_radius.str()).raw()) ||
            model_prescriptions !=
                producing_chart.at("prescriptions") ||
            std::any_of(model.q_coefficients.begin(),
                        model.q_coefficients.end(),
                        [](const ComplexBall& value) {
                          return !value.is_finite();
                        }) ||
            model.q_coefficients.front().contains_zero() ||
            std::any_of(model.initial_row_upper.begin(),
                        model.initial_row_upper.end(),
                        [](const Magnitude& value) {
                          return !value.is_finite();
                        }))
          throw std::invalid_argument(
              "checkpoint certified source model changed its retained binding");
        const auto matrix_size =
            static_cast<std::size_t>(model.dimension) * model.dimension;
        for (std::size_t lag = 0; lag < model.n_coefficients.size(); ++lag) {
          if (model.n_coefficients[lag].size() != matrix_size ||
              std::any_of(model.n_coefficients[lag].begin(),
                          model.n_coefficients[lag].end(),
                          [](const ComplexBall& value) {
                            return !value.is_finite();
                          }) ||
              (lag == 0 && std::any_of(
                  model.n_coefficients[lag].begin(),
                  model.n_coefficients[lag].end(),
                  [](const ComplexBall& value) { return !value.is_zero(); })) ||
              tail_majorant_detail::matrix_infinity_norm_upper(
                  model.n_coefficients[lag], model.dimension).dump_exact() !=
                  model.n_row_sum_upper[lag].dump_exact())
            throw std::invalid_argument(
                "checkpoint certified source model has an invalid N payload");
        }
        Rational replay_denominator(1);
        for (std::uint32_t exponent = 1;
             exponent <= witness_exponent; ++exponent) {
          replay_denominator *= Rational(2);
          const auto candidate = point_modulus +
              (producing_radius - point_modulus) / replay_denominator;
          const auto candidate_certificate =
              certify_regular_taylor_point_tail(
                  model,
                  RealEvaluationPoint::rational(producing_local.str()),
                  candidate.str(), producer_order);
          if ((exponent < witness_exponent &&
               (candidate_certificate.status !=
                    TailMajorantStatus::Inconclusive ||
                candidate_certificate.disk.status !=
                    TailMajorantStatus::Inconclusive)) ||
              (exponent == witness_exponent &&
               candidate_certificate.status !=
                   TailMajorantStatus::Certified))
            throw std::invalid_argument(
                "checkpoint certified value witness is not the first certifying dyadic candidate");
        }
        const auto replayed = certify_regular_taylor_point_tail(
            model, RealEvaluationPoint::rational(producing_local.str()),
            witness.str(), producer_order);
        const auto& certificate = as_object(
            tail.at("certificate"),
            "checkpoint certified value-tail certificate");
        require_exact_keys(
            certificate, {"status", "value", "theta", "disk", "detail"},
            "checkpoint certified value-tail certificate");
        const auto saved_value = parse_checkpoint_error_envelope(
            certificate.at("value"));
        const auto saved_theta = parse_checkpoint_error_envelope(
            certificate.at("theta"));
        const auto& disk = as_object(
            certificate.at("disk"),
            "checkpoint certified value-tail disk");
        require_exact_keys(
            disk,
            {"witness_radius_exact", "q_lower_exact",
             "ode_norm_upper_exact", "cauchy_circle_upper_exact",
             "detail"},
            "checkpoint certified value-tail disk");
        json::array replayed_circle;
        for (const auto& value : replayed.disk.cauchy_circle_upper)
          replayed_circle.emplace_back(value.dump_exact());
        if (replayed.status != TailMajorantStatus::Certified ||
            required_string(certificate, "status") != "certified" ||
            checkpoint_error_envelope_record(replayed.value) !=
                certificate.at("value") ||
            checkpoint_error_envelope_record(replayed.theta) !=
                certificate.at("theta") ||
            required_string(certificate, "detail") != replayed.detail ||
            required_string(disk, "witness_radius_exact") !=
                replayed.disk.witness_radius_exact ||
            required_string(disk, "q_lower_exact") !=
                replayed.disk.q_lower.dump_exact() ||
            required_string(disk, "ode_norm_upper_exact") !=
                replayed.disk.ode_norm_upper.dump_exact() ||
            disk.at("cauchy_circle_upper_exact") != replayed_circle ||
            required_string(disk, "detail") != replayed.disk.detail ||
            saved_value.guarantee != ErrorGuarantee::Certified ||
            saved_theta.guarantee != ErrorGuarantee::Certified ||
            !tail_majorant_detail::same_epsilon_window(
                saved_value.frame, model.epsilon) ||
            !tail_majorant_detail::same_epsilon_window(
                saved_theta.frame, model.epsilon))
          throw std::invalid_argument(
              "checkpoint certified value-tail theorem does not replay exactly");
        const auto& inflation = as_object(
            tail.at("inflation"),
            "checkpoint certified value inflation");
        require_exact_keys(
            inflation, {"gate", "retained_value", "inflated_value"},
            "checkpoint certified value inflation");
        if (required_string(inflation, "gate") !=
            "each-component-acb-add-error-mag-by-epsilon-row-v1")
          throw std::invalid_argument(
              "checkpoint certified value inflation gate changed");
        auto retained_value = parse_checkpoint_epsilon_vector(
            inflation.at("retained_value"),
            "checkpoint retained value before tail inflation");
        auto inflated_value = parse_checkpoint_epsilon_vector(
            inflation.at("inflated_value"),
            "checkpoint retained value after tail inflation");
        if (!retained_value.error.empty() || !inflated_value.error.empty() ||
            retained_value.dimension != model.dimension ||
            inflated_value.dimension != model.dimension ||
            !tail_majorant_detail::same_epsilon_window(
                retained_value.epsilon, model.epsilon) ||
            !tail_majorant_detail::same_epsilon_window(
                inflated_value.epsilon, model.epsilon))
          throw std::invalid_argument(
              "checkpoint certified value inflation frame changed");
        auto restored_incoming = std::dynamic_pointer_cast<
            StoredLocal<ComplexBall>>(value_handoff_plan->incoming);
        if (!restored_incoming ||
            restored_incoming->checkpoint_identity() !=
                model.local_checkpoint_identity ||
            restored_incoming->source_operator_identity() !=
                model.operator_identity ||
            restored_incoming->tail_model().status !=
                TailMajorantStatus::Certified ||
            !restored_incoming->tail_model().model.has_value() ||
            checkpoint_regular_tail_model_record(
                *restored_incoming->tail_model().model) !=
                tail.at("source_model"))
          throw std::invalid_argument(
              "checkpoint certified value handoff lost its restored incoming Acb theorem owner");
        EvaluationOptions replay_options;
        replay_options.imaginary_sign = value_handoff_plan->producing_rim;
        replay_options.compute_tail_estimate = false;
        const auto restored_evaluation =
            evaluate_local_solution_with_certified_tail(
                restored_incoming->solution(),
                *restored_incoming->tail_model().model,
                RealEvaluationPoint::rational(producing_local.str()),
                witness.str(), replay_options);
        const auto expected_rim = producing_local.sign() < 0
            ? value_handoff_plan->producing_rim : std::nullopt;
        if (restored_evaluation.tail.status !=
                TailMajorantStatus::Certified ||
            restored_evaluation.evaluation.imaginary_sign != expected_rim ||
            restored_evaluation.evaluation.value.dimension !=
                retained_value.dimension ||
            !tail_majorant_detail::same_epsilon_window(
                restored_evaluation.evaluation.value.epsilon,
                retained_value.epsilon) ||
            restored_evaluation.evaluation.value.coefficients.size() !=
                retained_value.coefficients.size())
          throw std::invalid_argument(
              "checkpoint certified retained center evaluation changed its frame or branch");
        for (std::size_t coefficient = 0;
             coefficient < retained_value.coefficients.size();
             ++coefficient)
          if (!acb_equal(
                  restored_evaluation.evaluation.value
                      .coefficients[coefficient].raw(),
                  retained_value.coefficients[coefficient].raw()))
            throw std::invalid_argument(
                "checkpoint certified retained center value differs from its restored incoming local");
        for (std::int64_t raw_power = model.epsilon.min_power;
             raw_power <= model.epsilon.complete_max; ++raw_power) {
          const auto power = static_cast<std::int32_t>(raw_power);
          const auto row = static_cast<std::size_t>(
              raw_power - model.epsilon.min_power);
          for (std::uint32_t component = 0; component < model.dimension;
               ++component) {
            auto expected = retained_value.at(power, component);
            replayed.value.absolute[row].add_error_to(expected);
            if (!acb_equal(expected.raw(),
                           inflated_value.at(power, component).raw()))
              throw std::invalid_argument(
                  "checkpoint certified value inflation does not replay exactly");
          }
        }
        certified_inflated_value = std::move(inflated_value);
      } else {
        if (expected_domain != "acb" ||
            required_string(tail, "output_status") != "unsupported" ||
            required_string(tail, "output_reason") !=
                "causal physical evolution has no framed recurrence tail theorem; the next hop must use exact framed fallback")
          throw std::invalid_argument(
              "checkpoint physical-evolution value handoff changed its output-tail contract");
        const Rational witness(
            required_string(tail, "witness_radius_exact"));
        const auto witness_exponent = as_u32(
            tail.at("witness_dyadic_inward_exponent"),
            "checkpoint physical-evolution witness exponent");
        Rational witness_denominator(1);
        for (std::uint32_t exponent = 0; exponent < witness_exponent;
             ++exponent)
          witness_denominator *= Rational(2);
        const auto point_modulus =
            exact_path_detail::abs(producing_local);
        const auto expected_witness = point_modulus +
            (producing_radius - point_modulus) / witness_denominator;
        if (witness_exponent == 0 || witness_exponent > 16 ||
            witness != expected_witness ||
            !(point_modulus < witness) ||
            !(witness < producing_radius))
          throw std::invalid_argument(
              "checkpoint physical-evolution witness is not strictly interior");
        auto model = parse_checkpoint_regular_tail_model(
            tail.at("source_model"));
        auto restored_incoming = std::dynamic_pointer_cast<
            StoredLocal<ComplexBall>>(value_handoff_plan->incoming);
        if (!restored_incoming ||
            restored_incoming->checkpoint_identity() !=
                model.local_checkpoint_identity ||
            restored_incoming->source_operator_identity() !=
                model.operator_identity ||
            restored_incoming->tail_model().status !=
                TailMajorantStatus::Certified ||
            !restored_incoming->tail_model().model.has_value() ||
            checkpoint_regular_tail_model_record(
                *restored_incoming->tail_model().model) !=
                tail.at("source_model") ||
            model.dimension != solution.dimension ||
            model.taylor_complete_max != producer_order ||
            model.operator_identity !=
                required_string(producing_chart, "identity") ||
            model.local_checkpoint_identity != required_string(
                as_object(derivation.at("incoming"),
                          "checkpoint physical-evolution incoming"),
                "checkpoint_identity"))
          throw std::invalid_argument(
              "checkpoint physical-evolution handoff lost its restored incoming Acb theorem owner");
        Rational replay_denominator(1);
        for (std::uint32_t exponent = 1;
             exponent <= witness_exponent; ++exponent) {
          replay_denominator *= Rational(2);
          const auto candidate = point_modulus +
              (producing_radius - point_modulus) /
                  replay_denominator;
          const auto candidate_certificate =
              certify_regular_taylor_point_tail(
                  model,
                  RealEvaluationPoint::rational(
                      producing_local.str()),
                  candidate.str(), producer_order);
          if ((exponent < witness_exponent &&
               (candidate_certificate.status !=
                    TailMajorantStatus::Inconclusive ||
                candidate_certificate.disk.status !=
                    TailMajorantStatus::Inconclusive)) ||
              (exponent == witness_exponent &&
               candidate_certificate.status !=
                   TailMajorantStatus::Certified))
            throw std::invalid_argument(
                "checkpoint physical-evolution witness is not the first certifying dyadic candidate");
        }
        EvaluationOptions replay_options;
        replay_options.imaginary_sign =
            value_handoff_plan->producing_rim;
        replay_options.compute_tail_estimate = false;
        auto restored_evaluation =
            evaluate_local_solution_with_certified_tail(
                restored_incoming->solution(),
                *restored_incoming->tail_model().model,
                RealEvaluationPoint::rational(
                    producing_local.str()),
                witness.str(), replay_options);
        const auto expected_rim = producing_local.sign() < 0
            ? value_handoff_plan->producing_rim : std::nullopt;
        if (restored_evaluation.tail.status !=
                TailMajorantStatus::Certified ||
            restored_evaluation.evaluation.imaginary_sign !=
                expected_rim ||
            restored_evaluation.evaluation.value.dimension !=
                model.dimension ||
            !tail_majorant_detail::same_epsilon_window(
                restored_evaluation.evaluation.value.epsilon,
                model.epsilon) ||
            restored_evaluation.tail.value.absolute.size() !=
                model.epsilon.width())
          throw std::invalid_argument(
              "checkpoint physical-evolution center evaluation changed its certified frame or branch");
        auto inflated =
            std::move(restored_evaluation.evaluation.value);
        inflated.error = ErrorEnvelope{};
        for (std::int64_t raw_power = model.epsilon.min_power;
             raw_power <= model.epsilon.complete_max; ++raw_power) {
          const auto power =
              static_cast<std::int32_t>(raw_power);
          const auto row = static_cast<std::size_t>(
              raw_power - model.epsilon.min_power);
          for (std::uint32_t component = 0;
               component < model.dimension; ++component)
            restored_evaluation.tail.value.absolute[row]
                .add_error_to(inflated.at(power, component));
        }
        if (solution.epsilon.complete_max !=
                inflated.epsilon.complete_max ||
            solution.epsilon.min_power >
                inflated.epsilon.min_power)
          throw std::invalid_argument(
              "checkpoint physical-evolution output cannot contain its restored center value");
        if (solution.epsilon.min_power <
            inflated.epsilon.min_power) {
          auto widened =
              physical_ode_detail::zero_epsilon_vector(
                  EpsilonWindow{
                      solution.epsilon.min_power,
                      inflated.epsilon.complete_max},
                  inflated.dimension);
          for (std::int64_t raw_power =
                   inflated.epsilon.min_power;
               raw_power <= inflated.epsilon.complete_max;
               ++raw_power) {
            const auto power =
                static_cast<std::int32_t>(raw_power);
            for (std::uint32_t component = 0;
                 component < inflated.dimension; ++component)
              widened.at(power, component) =
                  inflated.at(power, component);
          }
          inflated = std::move(widened);
        }
        certified_inflated_value = std::move(inflated);
      }
      const auto& accuracy = as_object(
          derivation.at("accuracy_contract"),
          "checkpoint value accuracy contract");
      require_exact_keys(
          accuracy,
          {"relative_error_max_exact", "gate", "acb_preflight_required"},
          "checkpoint value accuracy contract");
      const Rational relative_error_max(
          required_string(accuracy, "relative_error_max_exact"));
      std::optional<std::string> owner_relative_accuracy_max;
      if (prepared_owner)
        owner_relative_accuracy_max =
            prepared_owner
                ->regular_value_relative_accuracy_max_exact();
      else if (regular_equation_owner)
        owner_relative_accuracy_max =
            regular_equation_owner
                ->regular_value_relative_accuracy_max_exact();
      if (relative_error_max.sign() <= 0 ||
          !(relative_error_max < Rational(1)) ||
          required_string(accuracy, "relative_error_max_exact") !=
              relative_error_max.str() ||
          required_string(accuracy, "relative_error_max_exact") !=
              required_string(prototype, "relative_accuracy_max_exact") ||
          !owner_relative_accuracy_max.has_value() ||
          required_string(accuracy, "relative_error_max_exact") !=
              *owner_relative_accuracy_max ||
          required_string(accuracy, "gate") !=
              "component-radii-lte-threshold-times-max-one-upper-magnitude-v1" ||
          !accuracy.at("acb_preflight_required").is_bool() ||
          accuracy.at("acb_preflight_required").as_bool() !=
              (expected_domain == "acb") ||
          (certified_inflated_value.has_value() &&
           std::any_of(
               certified_inflated_value->coefficients.begin(),
               certified_inflated_value->coefficients.end(),
               [&](const ComplexBall& value) {
                 return !value_handoff_accurate(value,
                                                relative_error_max);
               })))
        throw std::invalid_argument(
            "checkpoint plan-value handoff accuracy contract is inconsistent");

      const auto parse_epsilon = [](const json::value& raw,
                                    const char* label) {
        const auto& epsilon = as_object(raw, label);
        require_exact_keys(epsilon, {"min", "max", "required_complete_max"},
                           label);
        EpsilonWindow window{as_i32(epsilon.at("min"), label),
                             as_i32(epsilon.at("max"), label)};
        (void)window.width();
        const auto required = as_i32(epsilon.at("required_complete_max"),
                                     label);
        if (required < window.min_power || required > window.complete_max)
          throw std::invalid_argument(
              std::string(label) + " has an inconsistent complete edge");
        return std::pair{window, required};
      };
      const auto [requested, requested_required] = parse_epsilon(
          derivation.at("epsilon"), "checkpoint value requested epsilon");
      const auto [evaluated, evaluated_required] = parse_epsilon(
          derivation.at("evaluated_epsilon"),
          "checkpoint value evaluated epsilon");
      if (evaluated_required != requested_required ||
          evaluated.min_power > requested.min_power ||
          evaluated.complete_max < requested.complete_max)
        throw std::invalid_argument(
            "checkpoint plan-value handoff lost evaluated epsilon completeness");
      const auto& incoming = as_object(
          derivation.at("incoming"), "checkpoint value incoming");
      const auto& output = as_object(
          derivation.at("output"), "checkpoint value output");
      require_exact_keys(
          incoming,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity", "epsilon", "taylor_complete_max",
           "top_valid", "analytic_metadata"},
          "checkpoint value incoming");
      require_exact_keys(
          output,
          {"checkpoint_identity", "chart", "source_operator_identity",
           "epsilon", "taylor_complete_max", "top_valid", "dimension"},
          "checkpoint value output");
      const auto& incoming_epsilon = as_object(
          incoming.at("epsilon"), "checkpoint value incoming epsilon");
      require_exact_keys(incoming_epsilon, {"min", "max"},
                         "checkpoint value incoming epsilon");
      const auto& output_epsilon = as_object(
          output.at("epsilon"), "checkpoint value output epsilon");
      if (required_string(incoming, "chart") !=
              required_string(producing_chart, "chart") ||
          required_string(incoming, "source_operator_identity") !=
              required_string(producing_chart, "identity") ||
          required_string(output, "checkpoint_identity") !=
              solution.checkpoint_identity ||
          required_string(output, "chart") != source_chart ||
          required_string(output, "source_operator_identity") !=
              source_operator_identity ||
          as_i32(output_epsilon.at("min"),
                 "checkpoint value output epsilon minimum") !=
              solution.epsilon.min_power ||
          as_i32(output_epsilon.at("max"),
                 "checkpoint value output epsilon maximum") !=
              solution.epsilon.complete_max ||
          as_u32(output.at("dimension"),
                 "checkpoint value output dimension") !=
              solution.dimension ||
          as_u32(output.at("taylor_complete_max"),
                 "checkpoint value output Taylor maximum") !=
              solution.taylor_complete_max ||
          parse_validity(output.at("top_valid")) < requested_required)
        throw std::invalid_argument(
            "checkpoint plan-value output disagrees with its retained tensor");
      if (certified_inflated_value.has_value()) {
        const auto saved_incoming_min = as_i32(
            incoming_epsilon.at("min"),
            "checkpoint value incoming epsilon minimum");
        if ((physical_evolution_handoff
                 ? saved_incoming_min <
                       certified_inflated_value->epsilon.min_power
                 : saved_incoming_min !=
                       certified_inflated_value->epsilon.min_power) ||
            as_i32(incoming_epsilon.at("max"),
                   "checkpoint value incoming epsilon maximum") !=
                certified_inflated_value->epsilon.complete_max ||
            certified_inflated_value->dimension != solution.dimension)
          throw std::invalid_argument(
              "checkpoint certified value source frame changed");
        if constexpr (!std::is_same_v<Scalar, ComplexBall>) {
          throw std::invalid_argument(
              "checkpoint certified value handoff is not an Acb local");
        } else {
          if (solution.sectors.size() != 1 ||
              solution.sectors.front().log_power != 0)
            throw std::invalid_argument(
                "checkpoint certified value output is not ordinary");
          for (std::int64_t raw_power = solution.epsilon.min_power;
               raw_power <= solution.epsilon.complete_max; ++raw_power) {
            const auto power = static_cast<std::int32_t>(raw_power);
            for (std::uint32_t component = 0;
                 component < solution.dimension; ++component) {
              const auto epsilon_index = static_cast<std::size_t>(
                  raw_power - solution.epsilon.min_power);
              const auto& output_constant =
                  solution.sectors.front().coefficients[
                      local_detail::sector_index(
                          solution, epsilon_index, 0, component)];
              if (!acb_equal(
                      output_constant.raw(),
                      certified_inflated_value->at(
                          power, component).raw()))
                throw std::invalid_argument(
                    "checkpoint certified inflated value differs from the receiver initial tensor");
            }
          }
        }
      }
      const auto& lineage = as_object(
          object.at("retained_owner_lineage"),
          "checkpoint plan-value owner lineage");
      require_exact_keys(
          lineage,
          {"tile_plan", "tile_plan_checkpoint_identity", "arm",
           "match_index", "incoming_checkpoint_identity",
           "handoff_provenance_identity"},
          "checkpoint plan-value owner lineage");
      if (lineage.at("tile_plan") != derivation.at("tile_plan") ||
          lineage.at("tile_plan_checkpoint_identity") !=
              derivation.at("tile_plan_checkpoint_identity") ||
          lineage.at("arm") != derivation.at("arm") ||
          lineage.at("match_index") != derivation.at("match") ||
          lineage.at("incoming_checkpoint_identity") !=
              incoming.at("checkpoint_identity") ||
          lineage.at("handoff_provenance_identity") !=
              derivation.at("provenance_identity"))
        throw std::invalid_argument(
            "checkpoint plan-value owner lineage is inconsistent");
      auto identity_input = derivation;
      const auto identity = required_string(
          derivation, "provenance_identity");
      identity_input.erase("provenance_identity");
      if (json::serialize(canonical_json_value(identity_input)) != identity)
        throw std::invalid_argument(
            "checkpoint plan-value derivation identity is inconsistent");
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
          (equation_owner
               ? required_string(
                     derivation,
                     "equation_owner_signature_identity") !=
                         equation_owner->owner_signature_identity() ||
                     required_string(
                         derivation, "equation_payload_identity") !=
                         equation_owner->physical_payload_identity()
               : !required_string(
                      derivation,
                      "equation_owner_signature_identity").starts_with(
                          "owner-operator-reference-v1:") ||
                     required_string(
                         derivation, "equation_payload_identity") !=
                         "no-physical-qc-payload-v1"))
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
      const bool bounded_incoming =
          incoming.if_contains("source_operator_reference") != nullptr;
      if (bounded_incoming)
        require_exact_keys(
            incoming,
            {"local", "checkpoint_identity",
             "source_operator_reference"},
            "checkpoint bounded materialized-local incoming");
      else
        require_exact_keys(
            incoming,
            {"local", "checkpoint_identity",
             "source_operator_identity"},
            "checkpoint legacy compact materialized-local incoming");
      (void)scoped_handle_id(required_string(incoming, "local"), "l:",
                             "compact materialized-local incoming");
      (void)required_string(incoming, "checkpoint_identity");
      if (bounded_incoming)
        (void)as_object(
            incoming.at("source_operator_reference"),
            "checkpoint incoming source-operator reference");
      else
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
        const bool bounded_source =
            source.if_contains("source_operator_reference") != nullptr;
        if (bounded_source)
          require_exact_keys(
              source,
              {"column", "local", "checkpoint_identity",
               "source_operator_reference"},
              "checkpoint bounded compact basis source");
        else
          require_exact_keys(
              source,
              {"column", "local", "checkpoint_identity",
               "source_operator_identity"},
              "checkpoint legacy compact basis source");
        if (as_u64(source.at("column"),
                   "checkpoint compact basis column") != column ||
            (bounded_source
                 ? !compact_matching_identity_reference_matches(
                       source.at("source_operator_reference"),
                       source_operator_identity,
                       "checkpoint compact basis source-operator reference")
                 : required_string(
                       source, "source_operator_identity") !=
                       source_operator_identity))
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
  const bool has_taylor_width =
      source.if_contains("matching_taylor_width") != nullptr;
  if (basis) {
    if (has_column_provenance && has_taylor_width)
      require_exact_keys(source,
          {"column", "local", "chart", "source_operator_identity",
           "checkpoint_identity",
           "requested_imaginary_sign", "effective_imaginary_sign",
           "matching_taylor_width",
           "analytic_metadata", "column_provenance"},
          "checkpoint Acb basis source");
    else if (has_column_provenance)
      require_exact_keys(source,
          {"column", "local", "chart", "source_operator_identity",
           "checkpoint_identity",
           "requested_imaginary_sign", "effective_imaginary_sign",
           "analytic_metadata", "column_provenance"},
          "checkpoint Acb basis source");
    else if (has_taylor_width)
      require_exact_keys(source,
          {"column", "local", "chart", "source_operator_identity",
           "checkpoint_identity",
           "requested_imaginary_sign", "effective_imaginary_sign",
           "matching_taylor_width",
           "analytic_metadata"}, "checkpoint Acb basis source");
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
    if (has_column_provenance && has_taylor_width)
      require_exact_keys(source,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity",
           "requested_imaginary_sign", "effective_imaginary_sign",
           "matching_taylor_width",
           "analytic_metadata", "column_provenance"},
          "checkpoint Acb incoming source");
    else if (has_column_provenance)
      require_exact_keys(source,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity",
           "requested_imaginary_sign", "effective_imaginary_sign",
           "analytic_metadata", "column_provenance"},
          "checkpoint Acb incoming source");
    else if (has_taylor_width)
      require_exact_keys(source,
          {"local", "chart", "source_operator_identity",
           "checkpoint_identity",
           "requested_imaginary_sign", "effective_imaginary_sign",
           "matching_taylor_width",
           "analytic_metadata"}, "checkpoint Acb incoming source");
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
  if (has_taylor_width &&
      as_u64(source.at("matching_taylor_width"),
             "checkpoint Acb matching Taylor width") == 0)
    throw std::invalid_argument(
        "checkpoint Acb match source has a zero Taylor-prefix width");
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

  const auto exact_match_provenance = [&](bool compact) {
    json::array provenance_basis;
    if (compact) {
      provenance_basis =
          compact_matching_source_references(basis_sources);
    } else {
      for (const auto& source : basis_sources)
        provenance_basis.push_back(source);
    }
    return json::object{
        {"schema", "diffexp2-native-exact-regular-local-match-v1"},
        {"checkpoint_identity", checkpoint_identity},
        {"basis", std::move(provenance_basis)},
        {"incoming",
         compact
             ? json::value(
                   compact_matching_source_reference(incoming_source))
             : json::value(incoming_source)},
        {"basis_point_exact", basis_point},
        {"incoming_point_exact", incoming_point},
        {"physical_match_point_exact", physical_point},
        {"epsilon", json::object{{"min", window.min_power},
                                  {"max", window.complete_max},
                                  {"required_complete_max",
                                   required_complete_max}}}};
  };
  const auto old_provenance = exact_match_provenance(false);
  const auto compact_provenance = exact_match_provenance(true);
  if (json::serialize(canonical_json_value(old_provenance)) !=
          provenance_identity &&
      json::serialize(canonical_json_value(compact_provenance)) !=
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
  const auto checkpoint_schema = required_string(object, "schema");
  const bool has_terminal_normal_frame_materialization =
      checkpoint_schema == "diffexp2-retained-acb-match-v5";
  const bool has_acb_right_materialization =
      has_terminal_normal_frame_materialization ||
      checkpoint_schema == "diffexp2-retained-acb-match-v4";
  const bool has_exact_right_materialization =
      has_acb_right_materialization ||
      checkpoint_schema == "diffexp2-retained-acb-match-v3";
  if (!has_exact_right_materialization &&
      checkpoint_schema != "diffexp2-retained-acb-match-v2")
    throw std::invalid_argument(
        "unsupported retained Acb-match checkpoint schema");
  auto shape = object;
  if (has_exact_right_materialization) {
    if (shape.if_contains(
            "exact_right_materialization_transformation") == nullptr)
      throw std::invalid_argument(
          "checkpoint Acb match v3 lost its exact-right materialization field");
    shape.erase("exact_right_materialization_transformation");
  }
  if (has_acb_right_materialization) {
    if (shape.if_contains(
            "acb_right_materialization_preconditioner") == nullptr)
      throw std::invalid_argument(
          "checkpoint Acb match v4 lost its right materialization preconditioner field");
    shape.erase("acb_right_materialization_preconditioner");
  }
  if (has_terminal_normal_frame_materialization) {
    if (shape.if_contains(
            "terminal_normal_frame_right_transformation") == nullptr ||
        shape.if_contains(
            "terminal_normal_frame_exact_right_transformation") == nullptr)
      throw std::invalid_argument(
          "checkpoint Acb match v5 lost a terminal normal-frame transformation field");
    shape.erase("terminal_normal_frame_right_transformation");
    shape.erase("terminal_normal_frame_exact_right_transformation");
  }
  const bool has_residual_frame_identity =
      object.if_contains("residual_frame_identity") != nullptr;
  const bool has_residual_certificate_identity =
      object.if_contains("residual_certificate_identity") != nullptr;
  if (has_residual_frame_identity && has_residual_certificate_identity)
    require_exact_keys(shape,
        {"schema", "handle", "checkpoint_identity", "provenance_identity",
         "exact_lattice_identity", "exact_lattice_provenance_identity",
         "exact_lattice_canonical_witness", "basis_sources",
         "incoming_source", "basis_chart", "incoming_chart",
         "basis_point_exact", "incoming_point_exact",
         "physical_match_point_exact", "matching_frame_identity",
         "residual_frame_identity", "residual_certificate_identity",
         "epsilon", "dimension", "relative_tolerance",
         "max_refinement_steps", "refined", "elapsed_ms"},
        "checkpoint retained Acb match");
  else if (has_residual_frame_identity)
    require_exact_keys(shape,
        {"schema", "handle", "checkpoint_identity", "provenance_identity",
         "exact_lattice_identity", "exact_lattice_provenance_identity",
         "exact_lattice_canonical_witness", "basis_sources",
         "incoming_source", "basis_chart", "incoming_chart",
         "basis_point_exact", "incoming_point_exact",
         "physical_match_point_exact", "matching_frame_identity",
         "residual_frame_identity", "epsilon", "dimension",
         "relative_tolerance", "max_refinement_steps", "refined",
         "elapsed_ms"}, "checkpoint retained Acb match");
  else {
    if (has_residual_certificate_identity)
      throw std::invalid_argument(
          "checkpoint Acb residual certificate lacks its frame identity");
    require_exact_keys(shape,
        {"schema", "handle", "checkpoint_identity", "provenance_identity",
         "exact_lattice_identity", "exact_lattice_provenance_identity",
         "exact_lattice_canonical_witness", "basis_sources",
         "incoming_source", "basis_chart", "incoming_chart",
         "basis_point_exact", "incoming_point_exact",
         "physical_match_point_exact", "matching_frame_identity",
         "epsilon", "dimension", "relative_tolerance",
         "max_refinement_steps", "refined", "elapsed_ms"},
        "checkpoint retained Acb match");
  }
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
  const auto residual_frame_identity = has_residual_frame_identity
      ? required_string(object, "residual_frame_identity")
      : matching_frame_identity;
  const auto residual_certificate_identity =
      has_residual_certificate_identity
      ? required_string(object, "residual_certificate_identity")
      : std::string();
  if (basis_chart.empty() || incoming_chart.empty() || basis_point.empty() ||
      incoming_point.empty() || physical_point.empty() ||
      matching_frame_identity.empty() || residual_frame_identity.empty())
    throw std::invalid_argument(
        "checkpoint Acb match lost chart or point provenance");
  if (has_residual_certificate_identity &&
      matching_frame_identity != residual_frame_identity &&
      residual_certificate_identity.empty())
    throw std::invalid_argument(
        "checkpoint Acb transformed residual lost its pushforward certificate");
  if (!residual_certificate_identity.empty()) {
    const auto raw_certificate =
        json::parse(residual_certificate_identity);
    const auto& certificate = as_object(
        raw_certificate,
        "checkpoint Acb residual pushforward certificate");
    const auto certificate_schema =
        required_string(certificate, "schema");
    const bool point_specialized_certificate =
        certificate_schema ==
            "diffexp2-acb-matching-residual-pushforward-v2";
    if (point_specialized_certificate)
      require_exact_keys(certificate,
          {"schema", "equation_operator_identity",
           "normal_frame_identity", "physical_frame_identity",
           "receiving_local_point", "blocks"},
          "checkpoint Acb residual pushforward certificate");
    else
      require_exact_keys(certificate,
          {"schema", "equation_operator_identity",
           "normal_frame_identity", "physical_frame_identity", "blocks"},
          "checkpoint Acb residual pushforward certificate");
    if ((certificate_schema !=
             "diffexp2-acb-matching-residual-pushforward-v1" &&
         !point_specialized_certificate) ||
        required_string(certificate, "normal_frame_identity") !=
            matching_frame_identity ||
        required_string(certificate, "physical_frame_identity") !=
            residual_frame_identity ||
        required_string(certificate, "equation_operator_identity").empty())
      throw std::invalid_argument(
          "checkpoint Acb residual pushforward certificate is inconsistent");
    if (json::serialize(canonical_json_value(certificate)) !=
        residual_certificate_identity)
      throw std::invalid_argument(
          "checkpoint Acb residual pushforward certificate is not canonical");
    const auto& certificate_blocks = as_array(
        certificate.at("blocks"),
        "checkpoint Acb residual pushforward blocks");
    if (certificate_blocks.empty())
      throw std::invalid_argument(
          "checkpoint Acb residual pushforward certificate has no exact blocks");
    std::set<std::uint32_t> block_indices;
    for (const auto& raw_block : certificate_blocks) {
      const auto& block = as_object(
          raw_block, "checkpoint Acb residual pushforward block");
      if (point_specialized_certificate)
        require_exact_keys(block,
            {"block", "vertices", "principal_identity",
             "source_transform_identity", "v_exact_identity",
             "vinv_exact_identity", "det_exact_identity",
             "to_reduced_transform_identity",
             "to_physical_transform_identity",
             "gauge_exact_identity", "gauge_inverse_exact_identity",
             "gauge_det_exact_identity", "assembly_exact_identity"},
            "checkpoint Acb residual pushforward block");
      else
        require_exact_keys(block,
            {"block", "principal_identity",
             "source_transform_identity", "v_exact_identity",
             "vinv_exact_identity", "det_exact_identity",
             "assembly_exact_identity"},
            "checkpoint Acb residual pushforward block");
      const auto index = as_u32(
          block.at("block"),
          "checkpoint Acb residual pushforward block index");
      if (!block_indices.insert(index).second ||
          required_string(block, "principal_identity").empty())
        throw std::invalid_argument(
            "checkpoint Acb residual pushforward block binding is invalid");
    }
    if (point_specialized_certificate) {
      const auto& point = as_object(
          certificate.at("receiving_local_point"),
          "checkpoint Acb residual pushforward point");
      require_exact_keys(
          point, {"exact", "sign"},
          "checkpoint Acb residual pushforward point");
      const auto point_sign = as_i32(
          point.at("sign"),
          "checkpoint Acb residual pushforward point sign");
      if (required_string(point, "exact") != basis_point ||
          point_sign < -1 || point_sign > 1)
        throw std::invalid_argument(
            "checkpoint Acb residual pushforward point binding is inconsistent");

      const auto raw_normal_frame =
          json::parse(matching_frame_identity);
      const auto& normal_frame = as_object(
          raw_normal_frame,
          "checkpoint Acb SCC matching normal frame");
      require_exact_keys(normal_frame,
          {"schema", "equation_operator_identity",
           "receiving_local_point", "blocks"},
          "checkpoint Acb SCC matching normal frame");
      if (required_string(normal_frame, "schema") !=
              "diffexp2-acb-scc-matching-normal-frame-v2" ||
          required_string(normal_frame, "equation_operator_identity") !=
              required_string(certificate,
                              "equation_operator_identity") ||
          normal_frame.at("receiving_local_point") !=
              certificate.at("receiving_local_point") ||
          normal_frame.at("blocks") != certificate.at("blocks"))
        throw std::invalid_argument(
            "checkpoint Acb residual pushforward differs from its point-specialized normal frame");
    }
  }
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
  const bool has_taylor_binding =
      incoming_source.if_contains("matching_taylor_width") != nullptr;
  for (const auto& source : basis_sources)
    if ((source.if_contains("matching_taylor_width") != nullptr) !=
        has_taylor_binding)
      throw std::invalid_argument(
          "checkpoint Acb match sources disagree on Taylor-prefix provenance");

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
            kNativeSingularSCCSaturationProofSchema ||
        saturation_witness_schema ==
            kNativeSingularSCCSaturationCompactProofSchema) {
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
  std::optional<ExactLaurentMatrix<Rational>>
      exact_right_materialization_transformation;
  if (has_exact_right_materialization) {
    const auto& raw_transformation =
        object.at("exact_right_materialization_transformation");
    if (!raw_transformation.is_null())
      exact_right_materialization_transformation =
          parse_checkpoint_exact_laurent_matrix(
              raw_transformation, dimension,
              "checkpoint Acb exact-right materialization transformation");
  }
  std::optional<ExactLaurentMatrix<Rational>>
      terminal_normal_frame_exact_right_transformation;
  if (has_terminal_normal_frame_materialization) {
    const auto& raw_terminal_exact =
        object.at("terminal_normal_frame_exact_right_transformation");
    if (!raw_terminal_exact.is_null())
      terminal_normal_frame_exact_right_transformation =
          parse_checkpoint_exact_laurent_matrix(
              raw_terminal_exact, dimension,
              "checkpoint Acb terminal normal-frame exact-right transformation");
  }
  // A compact Rational-shadow proof records its exact formal determinant
  // valuation and column valuations, but the latter do not reproduce the
  // full rectangular formal-module saturation transformation.  The v4 match
  // record stores that authoritative exact-right transformation separately.
  // Rebind the reconstructed proof shell to it before the ordinary identity
  // check below; otherwise every nontrivial CASE-P saturation is writable but
  // cannot be restored.
  if (exact_lattice.exact_formal_negative_coefficients_are_zero &&
      exact_right_materialization_transformation.has_value())
    exact_lattice.saturation.transformation =
        *exact_right_materialization_transformation;
  if (terminal_normal_frame_exact_right_transformation.has_value())
    exact_lattice.saturation.transformation =
        *terminal_normal_frame_exact_right_transformation;
  std::optional<ExactLaurentMatrix<ComplexBall>>
      acb_right_materialization_preconditioner;
  if (has_acb_right_materialization) {
    const auto& raw_preconditioner =
        object.at("acb_right_materialization_preconditioner");
    if (!raw_preconditioner.is_null())
      acb_right_materialization_preconditioner =
          parse_checkpoint_acb_laurent_matrix(
              raw_preconditioner, dimension,
              "checkpoint Acb right materialization preconditioner");
  }
  std::optional<FiniteLaurentMatrix<ComplexBall>>
      terminal_normal_frame_right_transformation;
  if (has_terminal_normal_frame_materialization) {
    const auto& raw_terminal =
        object.at("terminal_normal_frame_right_transformation");
    if (!raw_terminal.is_null())
      terminal_normal_frame_right_transformation =
          parse_checkpoint_frame_matrix<ComplexBall>(
              raw_terminal, dimension, dimension,
              "checkpoint Acb terminal normal-frame right transformation");
  }
  // A normalized terminal match may use either the exact formal saturation
  // or the certified Acb Levelt transformation reconstructed from its
  // canonical witness.  Keep the latter when present: the endpoint
  // association is then V*T_acb*P, not V*T_exact*P.
  if (exact_right_materialization_transformation.has_value() &&
      checkpoint_exact_laurent_matrix_record(
          *exact_right_materialization_transformation) !=
          checkpoint_exact_laurent_matrix_record(
              exact_lattice.saturation.transformation))
    throw std::invalid_argument(
        "checkpoint Acb exact-right materialization transformation differs from its restored exact lattice");
  if (acb_right_materialization_preconditioner.has_value()) {
    const bool physical_exact_right =
        exact_right_materialization_transformation.has_value() &&
        expected_singular_request.has_value() &&
        matching_frame_identity == "physical-parent-frame" &&
        residual_frame_identity == "physical-parent-frame";
    const bool terminal_normal_right =
        terminal_normal_frame_right_transformation.has_value() &&
        expected_singular_request.has_value() &&
        matching_frame_identity != "physical-parent-frame" &&
        residual_frame_identity == "physical-parent-frame";
    if (!physical_exact_right && !terminal_normal_right)
      throw std::invalid_argument(
          "checkpoint Acb right materialization preconditioner is not bound to an eligible singular physical or terminal normal frame");
    for (const auto& row :
         *acb_right_materialization_preconditioner)
      for (const auto& entry : row)
        for (const auto& [power, coefficient] : entry.terms()) {
          (void)coefficient;
          if (power != 0)
            throw std::invalid_argument(
                "checkpoint Acb right materialization preconditioner is not epsilon-constant");
        }
  }
  if (terminal_normal_frame_right_transformation.has_value() &&
      (!expected_singular_request.has_value() ||
       matching_frame_identity == "physical-parent-frame" ||
       residual_frame_identity != "physical-parent-frame" ||
       exact_right_materialization_transformation.has_value()))
    throw std::invalid_argument(
        "checkpoint Acb terminal normal-frame transformation is not bound to a certified singular SCC normal frame");
  if (terminal_normal_frame_right_transformation.has_value() !=
          terminal_normal_frame_exact_right_transformation.has_value())
    throw std::invalid_argument(
        "checkpoint Acb terminal normal-frame finite and exact transformations must be retained together");

  const auto lattice_provenance = [&](bool compact) {
    json::array exact_binding_basis;
    if (compact) {
      exact_binding_basis =
          compact_matching_source_references(basis_sources);
    } else {
      for (const auto& source : basis_sources)
        exact_binding_basis.push_back(source);
    }
    return json::object{
        {"schema", "diffexp2-retained-exact-lattice-binding-v1"},
        {"witness_schema", exact_lattice.witness_schema},
        {"witness_identity", exact_lattice_identity},
        {"basis", std::move(exact_binding_basis)},
        {"basis_point_exact", basis_point},
        {"physical_match_point_exact", physical_point},
        {"matching_frame_identity", matching_frame_identity},
        {"epsilon", json::object{{"min", window.min_power},
                                  {"max", window.complete_max}}}};
  };
  const auto old_lattice_provenance = lattice_provenance(false);
  const auto compact_lattice_provenance = lattice_provenance(true);
  if (json::serialize(canonical_json_value(old_lattice_provenance)) !=
          exact_lattice_provenance_identity &&
      json::serialize(canonical_json_value(compact_lattice_provenance)) !=
          exact_lattice_provenance_identity)
    throw std::invalid_argument(
        "checkpoint exact lattice provenance identity is inconsistent");

  const auto acb_match_provenance = [&](bool compact) {
    json::array provenance_basis;
    if (compact) {
      provenance_basis =
          compact_matching_source_references(basis_sources);
    } else {
      for (const auto& source : basis_sources)
        provenance_basis.push_back(source);
    }
    json::object provenance{
        {"schema", "diffexp2-native-refined-acb-local-match-v1"},
        {"checkpoint_identity", checkpoint_identity},
        {"basis", std::move(provenance_basis)},
        {"incoming",
         compact
             ? json::value(
                   compact_matching_source_reference(incoming_source))
             : json::value(incoming_source)},
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
    if (has_residual_frame_identity)
      provenance["residual_frame_identity"] = residual_frame_identity;
    if (!residual_certificate_identity.empty())
      provenance["residual_certificate_identity"] =
          residual_certificate_identity;
    if (acb_right_materialization_preconditioner.has_value())
      provenance["right_materialization_preconditioner"] =
          checkpoint_acb_laurent_matrix_record(
              *acb_right_materialization_preconditioner);
    if (terminal_normal_frame_right_transformation.has_value())
      provenance["terminal_normal_frame_right_transformation"] =
          checkpoint_frame_matrix_record(
              *terminal_normal_frame_right_transformation);
    if (terminal_normal_frame_exact_right_transformation.has_value())
      provenance["terminal_normal_frame_exact_right_transformation"] =
          checkpoint_exact_laurent_matrix_record(
              *terminal_normal_frame_exact_right_transformation);
    return provenance;
  };
  const auto old_provenance = acb_match_provenance(false);
  const auto compact_provenance = acb_match_provenance(true);
  if (json::serialize(canonical_json_value(old_provenance)) !=
          provenance_identity &&
      json::serialize(canonical_json_value(compact_provenance)) !=
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
      incoming_point, physical_point, matching_frame_identity,
      residual_frame_identity, residual_certificate_identity, window,
      required_complete_max,
      dimension, relative_tolerance,
      static_cast<std::size_t>(max_refinement_steps),
      std::move(exact_lattice.saturation),
      std::move(exact_right_materialization_transformation),
      std::move(exact_lattice.acb_transformation),
      std::move(acb_right_materialization_preconditioner),
      std::move(terminal_normal_frame_right_transformation),
      std::move(refined), elapsed_ms);
}
