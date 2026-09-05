Rational parse_exact_path_rational(const json::value& raw,
                                   const char* label) {
  if (!raw.is_string() || raw.as_string().empty())
    throw std::invalid_argument(std::string(label) +
                                " must be a nonempty exact rational string");
  try {
    return Rational(std::string(raw.as_string()));
  } catch (const std::invalid_argument&) {
    throw std::invalid_argument(std::string(label) +
                                " is not an exact rational number");
  }
}

std::vector<Rational> parse_exact_path_points(const json::value& raw,
                                              const char* label) {
  std::vector<Rational> points;
  for (const auto& value : as_array(raw, label))
    points.push_back(parse_exact_path_rational(value, label));
  return points;
}

ExactPathTopology parse_exact_path_topology(const json::value& raw) {
  const auto& object = as_object(raw, "native exact path topology");
  require_exact_keys(object,
      {"singular_points", "boundary_points", "complex_projections",
       "branch_sheets"}, "native exact path topology");
  ExactPathTopology topology;
  topology.singular_points = parse_exact_path_points(
      object.at("singular_points"), "path singular point");
  topology.boundary_points = parse_exact_path_points(
      object.at("boundary_points"), "path boundary point");
  for (const auto& raw_projection : as_array(
           object.at("complex_projections"), "path complex projections")) {
    const auto& projection = as_object(
        raw_projection, "path complex projection");
    require_exact_keys(projection,
        {"source_identity", "real_part_exact",
         "imaginary_magnitude_exact", "retain_minus_imaginary",
         "retain_real_part", "retain_plus_imaginary"},
        "path complex projection");
    topology.complex_projections.push_back(ExactComplexProjection{
        required_string(projection, "source_identity"),
        parse_exact_path_rational(projection.at("real_part_exact"),
                                  "projection real part"),
        parse_exact_path_rational(
            projection.at("imaginary_magnitude_exact"),
            "projection imaginary magnitude"),
        projection.at("retain_minus_imaginary").as_bool(),
        projection.at("retain_real_part").as_bool(),
        projection.at("retain_plus_imaginary").as_bool()});
  }
  for (const auto& raw_sheet : as_array(
           object.at("branch_sheets"), "path branch sheets")) {
    const auto& sheet = as_object(raw_sheet, "path branch sheet");
    require_exact_keys(sheet, {"factor_exact", "sign"},
                       "path branch sheet");
    topology.branch_sheets.push_back(ExactBranchSheet{
        required_string(sheet, "factor_exact"),
        as_i32(sheet.at("sign"), "path branch sign")});
  }
  return topology;
}

std::string retained_plan_owner_handle(
    const RetainedPlanChartBinding::Owner& owner) {
  return std::visit(
      [](const auto& typed) {
        if (typed == nullptr)
          throw std::logic_error(
              "native tile planning received a null retained owner");
        return typed->handle();
      },
      owner);
}

std::string retained_plan_owner_identity(
    const RetainedPlanChartBinding::Owner& owner) {
  return std::visit(
      [](const auto& typed) {
        if (typed == nullptr)
          throw std::logic_error(
              "native tile planning received a null retained owner");
        return typed->exact_identity();
      },
      owner);
}

std::string retained_plan_owner_geometry_record(
    const RetainedPlanChartBinding::Owner& owner) {
  return std::visit(
      [](const auto& typed) -> std::string {
        using Owner = typename std::decay_t<decltype(typed)>::element_type;
        if (typed == nullptr)
          throw std::logic_error(
              "native tile planning received a null retained owner");
        if constexpr (std::is_same_v<Owner, PreparedChartBase>) {
          if (!typed->geometry_record().has_value())
            throw std::invalid_argument(
                "native tile planning requires retained exact chart geometry");
          return *typed->geometry_record();
        } else {
          return typed->geometry_record();
        }
      },
      owner);
}

RetainedPlanChartBinding bind_plan_chart(
    const RetainedPlanChartBinding::Owner& owner,
    const ExactPathTopology& topology, const json::object& planning,
    const json::object& certified) {
  const auto handle = retained_plan_owner_handle(owner);
  const auto exact_identity = retained_plan_owner_identity(owner);
  const auto geometry_value = json::parse(
      retained_plan_owner_geometry_record(owner));
  const auto& geometry = as_object(
      geometry_value, "retained native tile chart geometry");
  if (geometry.at("infinite_radius").as_bool())
    throw std::invalid_argument(
        "native exact tile planning currently requires finite chart radii");

  RetainedPlanChartBinding binding;
  binding.handle = handle;
  binding.exact_identity = exact_identity;
  binding.owner = owner;
  binding.geometry.identity = exact_identity;
  require_exact_keys(planning,
      {"center_exact", "scale_exact", "radius_exact",
       "certificate_identity"}, "tile planning chart");
  binding.geometry.center = parse_exact_path_rational(
      planning.at("center_exact"), "tile planning chart center");
  binding.geometry.scale = parse_exact_path_rational(
      planning.at("scale_exact"), "tile planning chart scale");
  binding.geometry.radius = parse_exact_path_rational(
      planning.at("radius_exact"), "tile planning chart radius");
  binding.planning_certificate_identity = required_string(
      planning, "certificate_identity");
  if (binding.planning_certificate_identity.empty())
    throw std::invalid_argument(
        "tile planning chart requires a nonempty exact certificate identity");
  require_exact_keys(certified,
      {"index", "center", "scale", "radius"},
      "certified algebraic tile chart");
  const auto& center = as_object(certified.at("center"),
                                 "certified chart center");
  const auto& scale = as_object(certified.at("scale"),
                                "certified chart scale");
  const auto& radius = as_object(certified.at("radius"),
                                 "certified chart radius");
  for (const auto* point : {&center, &scale, &radius})
    require_exact_keys(*point, {"exact", "value", "sign"},
                       "certified chart scalar");
  binding.local_geometry.center_exact = required_string(center, "exact");
  binding.local_geometry.scale_exact = required_string(scale, "exact");
  binding.local_geometry.radius_exact = required_string(radius, "exact");
  binding.local_geometry.infinite_radius = false;
  binding.center_numeric = parse_scalar<ComplexBall>(center.at("value"));
  binding.scale_numeric = parse_scalar<ComplexBall>(scale.at("value"));
  binding.local_geometry.radius = parse_scalar<ComplexBall>(
      radius.at("value"));
  binding.scale_sign = as_i32(scale.at("sign"), "certified scale sign");
  if (binding.scale_sign != -1 && binding.scale_sign != 1)
    throw std::invalid_argument(
        "certified chart scale requires exact sign +/-1");
  if (required_string(geometry, "center_exact") !=
          binding.local_geometry.center_exact ||
      required_string(geometry, "scale_exact") !=
          binding.local_geometry.scale_exact ||
      required_string(geometry, "radius_exact") !=
          required_string(radius, "exact"))
    throw std::invalid_argument(
        "certified algebraic chart identity differs from its retained equation owner");
  if (!local_detail::exactly_real(binding.center_numeric) ||
      !local_detail::exactly_real(binding.scale_numeric) ||
      !local_detail::exactly_real(binding.local_geometry.radius) ||
      binding.scale_numeric.contains_zero() ||
      !arb_is_positive(acb_realref(binding.local_geometry.radius.raw())))
    throw std::invalid_argument(
        "certified algebraic chart specializations are not provably real/nondegenerate");
  const auto numeric_scale_sign =
      arb_is_positive(acb_realref(binding.scale_numeric.raw())) ? 1
      : arb_is_negative(acb_realref(binding.scale_numeric.raw())) ? -1 : 0;
  if (numeric_scale_sign != binding.scale_sign ||
      binding.geometry.scale.sign() != binding.scale_sign)
    throw std::invalid_argument(
        "certified chart scale sign contradicts its specialization or planning surrogate");
  std::optional<Rational> exact_scale;
  try {
    exact_scale = Rational(binding.local_geometry.scale_exact);
  } catch (const std::invalid_argument&) {
    // A genuinely algebraic exact identity is intentionally opaque.
  }
  if (exact_scale.has_value() && exact_scale->sign() != binding.scale_sign)
    throw std::invalid_argument(
        "certified rational chart scale sign contradicts its exact identity");
  binding.geometry.singular_center = std::any_of(
      topology.singular_points.begin(), topology.singular_points.end(),
      [&](const Rational& point) { return point == binding.geometry.center; });
  for (const auto& raw_prescription : as_array(
           geometry.at("prescriptions"), "tile chart prescriptions")) {
    const auto& prescription = as_object(
        raw_prescription, "tile chart prescription");
    binding.prescriptions.push_back(Prescription{
        required_string(prescription, "factor_exact"),
        as_i32(prescription.at("sign"), "tile prescription sign"),
        as_u32(prescription.at("multiplicity"),
               "tile prescription multiplicity"),
        as_i32(prescription.at("leading_coefficient_sign"),
               "tile prescription leading coefficient sign")});
  }
  // Every exact prepared prescription must be represented on the path.  The
  // topology may additionally retain factors which are inactive in this
  // chart, but it may never silently alter or drop an active sheet.
  for (const auto& prescription : binding.prescriptions) {
    const auto found = std::find_if(
        topology.branch_sheets.begin(), topology.branch_sheets.end(),
        [&](const ExactBranchSheet& sheet) {
          return sheet.factor_exact == prescription.factor_exact;
        });
    if (found == topology.branch_sheets.end() ||
        found->imaginary_sign != prescription.sign)
      throw std::invalid_argument(
          "native tile topology does not reproduce a prepared chart branch prescription");
  }
  (void)exact_plan_rim(binding.prescriptions, binding.geometry.scale);
  return binding;
}

json::object certified_rational_scalar(const Rational& value) {
  return json::object{
      {"exact", value.str()},
      {"value", json::array{value.str(), "0"}},
      {"sign", value.sign()}};
}

// Compatibility path for the existing exact-rational protocol.  It is
// deliberately implemented by constructing the same certified records as
// the algebraic protocol, so all downstream integration and checkpoint code
// has one representation and one validation path.
RetainedPlanChartBinding bind_plan_chart(
    const RetainedPlanChartBinding::Owner& owner,
    const ExactPathTopology& topology) {
  const auto geometry_value = json::parse(
      retained_plan_owner_geometry_record(owner));
  const auto& geometry = as_object(
      geometry_value, "retained native tile chart geometry");
  if (geometry.at("infinite_radius").as_bool())
    throw std::invalid_argument(
        "native exact tile planning currently requires finite chart radii");
  const Rational center(required_string(geometry, "center_exact"));
  const Rational scale(required_string(geometry, "scale_exact"));
  const Rational radius(required_string(geometry, "radius_exact"));
  json::object planning{
      {"center_exact", center.str()}, {"scale_exact", scale.str()},
      {"radius_exact", radius.str()},
      {"certificate_identity",
       "legacy-rational:" + retained_plan_owner_identity(owner) + ":" +
           center.str() + ":" + scale.str() + ":" + radius.str()}};
  json::object certified{
      {"index", 0}, {"center", certified_rational_scalar(center)},
      {"scale", certified_rational_scalar(scale)},
      {"radius", certified_rational_scalar(radius)}};
  return bind_plan_chart(owner, topology, planning, certified);
}

std::vector<std::string> parse_plan_chart_handles(const json::object& arm) {
  std::vector<std::string> handles;
  for (const auto& raw : as_array(arm.at("charts"), "tile arm charts")) {
    if (!raw.is_string() || raw.as_string().empty())
      throw std::invalid_argument(
          "native tile arm chart handles must be nonempty strings");
    handles.emplace_back(raw.as_string());
  }
  if (handles.empty())
    throw std::invalid_argument("native tile arm requires at least one chart");
  return handles;
}

std::pair<ExactArmRequest, std::vector<RetainedPlanChartBinding>>
parse_retained_arm_request(
    const json::object& arm,
    const std::vector<RetainedPlanChartBinding::Owner>& charts) {
  const bool has_planning = arm.if_contains("planning_charts") != nullptr;
  const bool has_certified = arm.if_contains("certified_geometry") != nullptr;
  if (has_planning != has_certified)
    throw std::invalid_argument(
        "native tile arm must provide planning and certified geometry together");
  if (has_planning)
    require_exact_keys(arm,
        {"from_exact", "to_exact", "charts", "planning_charts",
         "certified_geometry", "topology"}, "native tile arm");
  else
    require_exact_keys(arm,
        {"from_exact", "to_exact", "charts", "topology"},
        "native tile arm");
  const auto handles = parse_plan_chart_handles(arm);
  if (handles.size() != charts.size())
    throw std::invalid_argument(
        "resolved native tile chart count differs from its request");
  ExactArmRequest request;
  request.from = parse_exact_path_rational(arm.at("from_exact"),
                                           "tile arm start");
  request.to = parse_exact_path_rational(arm.at("to_exact"),
                                         "tile arm end");
  request.topology = parse_exact_path_topology(arm.at("topology"));
  std::vector<RetainedPlanChartBinding> bindings;
  bindings.reserve(charts.size());
  request.charts.reserve(charts.size());
  for (std::size_t index = 0; index < charts.size(); ++index) {
    if (retained_plan_owner_handle(charts[index]) != handles[index])
      throw std::logic_error("resolved tile chart handle changed");
    RetainedPlanChartBinding binding;
    if (has_planning) {
      const auto& raw_planning = as_array(
          arm.at("planning_charts"), "tile planning charts");
      const auto& raw_certified = as_object(
          arm.at("certified_geometry"), "certified algebraic arm geometry");
      const auto& certified_charts = as_array(
          raw_certified.at("charts"), "certified algebraic arm charts");
      if (raw_planning.size() != charts.size() ||
          certified_charts.size() != charts.size())
        throw std::invalid_argument(
            "planning/certified chart count differs from retained owners");
      const auto& planning = as_object(raw_planning[index],
                                       "tile planning chart");
      const auto& certified = as_object(certified_charts[index],
                                        "certified algebraic chart");
      if (as_u64(certified.at("index"), "certified chart index") != index)
        throw std::invalid_argument(
            "certified algebraic charts are not in exact owner order");
      binding = bind_plan_chart(charts[index], request.topology,
                                planning, certified);
    } else {
      binding = bind_plan_chart(charts[index], request.topology);
    }
    request.charts.push_back(binding.geometry);
    bindings.push_back(std::move(binding));
  }
  return {std::move(request), std::move(bindings)};
}

RealEvaluationPoint parse_certified_local_point(const json::value& raw,
                                                const char* label) {
  const auto& point = as_object(raw, label);
  require_exact_keys(point, {"exact", "value", "sign"}, label);
  const auto exact = required_string(point, "exact");
  const auto numeric = parse_scalar<ComplexBall>(point.at("value"));
  const auto sign = as_i32(point.at("sign"), label);
  std::optional<Rational> rational;
  try {
    rational = Rational(exact);
  } catch (const std::invalid_argument&) {
    // A genuine algebraic exact identity is intentionally opaque to C++.
  }
  if (rational.has_value()) {
    auto result = RealEvaluationPoint::rational(rational->str());
    const auto signed_result = result.sign < 0
        ? -result.modulus : result.modulus;
    if (result.sign != sign ||
        !acb_overlaps(signed_result.raw(), numeric.raw()))
      throw std::invalid_argument(
          std::string(label) +
          " rational identity contradicts its certified sign/specialization");
    return result;
  }
  return RealEvaluationPoint::certified(exact, numeric, sign);
}

ComplexBall signed_certified_value(const RealEvaluationPoint& point) {
  return point.sign < 0 ? -point.modulus : point.modulus;
}

void validate_certified_local_binding(
    const RetainedPlanChartBinding& chart,
    const RealEvaluationPoint& local, const ComplexBall& physical,
    const char* label) {
  if (!arb_lt(acb_realref(local.modulus.raw()),
              acb_realref(chart.local_geometry.radius.raw())))
    throw std::invalid_argument(std::string(label) +
        " is not provably inside its retained algebraic chart radius");
  const auto mapped = chart.center_numeric +
      chart.scale_numeric * signed_certified_value(local);
  if (!acb_overlaps(mapped.raw(), physical.raw()))
    throw std::invalid_argument(std::string(label) +
        " specialization does not map to its certified physical point");
}

bool certified_tile_zero_length(const RetainedArmPlan& arm,
                                std::size_t tile_index) {
  const auto& tile = arm.certified_tiles.at(tile_index);
  const bool physical_zero =
      tile.physical_begin_exact == tile.physical_end_exact;
  const bool local_zero =
      tile.local_begin.exact_coordinate == tile.local_end.exact_coordinate;
  if (physical_zero != local_zero)
    throw std::invalid_argument(
        "certified transport tile has inconsistent physical/local zero-length identities");
  if (physical_zero) return true;

  // Every nonzero certified tile must preserve the arm orientation.  Exact
  // rational physical endpoints are authoritative; genuinely algebraic
  // endpoints require strict separation of their rigorous real balls.
  std::optional<Rational> begin_exact;
  std::optional<Rational> end_exact;
  try {
    begin_exact = Rational(tile.physical_begin_exact);
    end_exact = Rational(tile.physical_end_exact);
  } catch (const std::invalid_argument&) {
    begin_exact.reset();
    end_exact.reset();
  }
  std::int32_t direction = 0;
  if (begin_exact.has_value() && end_exact.has_value()) {
    direction = *begin_exact < *end_exact
        ? 1 : *end_exact < *begin_exact ? -1 : 0;
  } else {
    if (!local_detail::exactly_real(tile.physical_begin_numeric) ||
        !local_detail::exactly_real(tile.physical_end_numeric))
      throw std::invalid_argument(
          "certified transport tile physical endpoints are not rigorously real");
    direction = arb_lt(acb_realref(tile.physical_begin_numeric.raw()),
                       acb_realref(tile.physical_end_numeric.raw()))
        ? 1 : arb_lt(acb_realref(tile.physical_end_numeric.raw()),
                     acb_realref(tile.physical_begin_numeric.raw()))
        ? -1 : 0;
  }
  if (direction == 0 || direction != arm.exact.direction)
    throw std::invalid_argument(
        "certified transport tile direction contradicts its retained arm");
  return false;
}

void attach_certified_arm_geometry(const json::object& raw_arm,
                                   RetainedArmPlan& arm) {
  const auto& certified = as_object(raw_arm.at("certified_geometry"),
                                    "certified algebraic arm geometry");
  require_exact_keys(certified,
      {"schema", "exact_identity", "charts", "matches", "tiles"},
      "certified algebraic arm geometry");
  if (required_string(certified, "schema") !=
          "diffexp3-wolfram-certified-algebraic-arm-v1")
    throw std::invalid_argument(
        "unsupported certified algebraic arm geometry schema");
  arm.certified_geometry_identity = required_string(
      certified, "exact_identity");
  if (arm.certified_geometry_identity.empty())
    throw std::invalid_argument(
        "certified algebraic arm geometry identity cannot be empty");
  arm.planning_charts = as_array(
      raw_arm.at("planning_charts"), "tile planning charts");
  arm.certified_geometry = certified;
  const auto& raw_matches = as_array(certified.at("matches"),
                                     "certified algebraic matches");
  const auto& raw_tiles = as_array(certified.at("tiles"),
                                   "certified algebraic tiles");
  if (raw_matches.size() != arm.exact.matches.size() ||
      raw_tiles.size() != arm.exact.tiles.size())
    throw std::invalid_argument(
        "certified algebraic geometry does not reproduce the exact plan topology");
  arm.certified_matches.reserve(raw_matches.size());
  for (std::size_t index = 0; index < raw_matches.size(); ++index) {
    const auto& raw = as_object(raw_matches[index],
                                "certified algebraic match");
    require_exact_keys(raw,
        {"index", "physical", "producing_local", "receiving_local"},
        "certified algebraic match");
    if (as_u64(raw.at("index"), "certified match index") != index)
      throw std::invalid_argument(
          "certified algebraic matches are not in exact plan order");
    const auto& physical_record = as_object(raw.at("physical"),
                                            "certified match physical point");
    const auto physical_exact = required_string(physical_record, "exact");
    const auto physical_numeric = parse_scalar<ComplexBall>(
        physical_record.at("value"));
    const auto producing_local = parse_certified_local_point(
        raw.at("producing_local"), "certified producing local point");
    const auto receiving_local = parse_certified_local_point(
        raw.at("receiving_local"), "certified receiving local point");
    const auto& planned = arm.exact.matches[index];
    const auto& producing = arm.charts.at(planned.producing_chart);
    const auto& receiving = arm.charts.at(planned.receiving_chart);
    validate_certified_local_binding(producing, producing_local,
                                     physical_numeric,
                                     "certified producing match point");
    validate_certified_local_binding(receiving, receiving_local,
                                     physical_numeric,
                                     "certified receiving match point");
    arm.certified_matches.push_back(CertifiedPlanMatch{
        {producing_local, physical_exact, physical_numeric,
         raw.at("producing_local").as_object().at("value")},
        {receiving_local, physical_exact, physical_numeric,
         raw.at("receiving_local").as_object().at("value")}});
  }
  arm.certified_tiles.reserve(raw_tiles.size());
  for (std::size_t index = 0; index < raw_tiles.size(); ++index) {
    const auto& raw = as_object(raw_tiles[index],
                                "certified algebraic tile");
    require_exact_keys(raw,
        {"index", "chart_index", "physical_begin", "physical_end",
         "local_begin", "local_end"}, "certified algebraic tile");
    if (as_u64(raw.at("index"), "certified tile index") != index ||
        as_u64(raw.at("chart_index"), "certified tile chart index") !=
            arm.exact.tiles[index].chart)
      throw std::invalid_argument(
          "certified algebraic tile order differs from the exact plan topology");
    const auto& begin = as_object(raw.at("physical_begin"),
                                  "certified tile physical begin");
    const auto& end = as_object(raw.at("physical_end"),
                                "certified tile physical end");
    CertifiedPlanTile tile;
    tile.physical_begin_exact = required_string(begin, "exact");
    tile.physical_end_exact = required_string(end, "exact");
    tile.physical_begin_numeric = parse_scalar<ComplexBall>(
        begin.at("value"));
    tile.physical_end_numeric = parse_scalar<ComplexBall>(end.at("value"));
    tile.local_begin = parse_certified_local_point(
        raw.at("local_begin"), "certified tile local begin");
    tile.local_end = parse_certified_local_point(
        raw.at("local_end"), "certified tile local end");
    const auto& chart = arm.charts.at(arm.exact.tiles[index].chart);
    validate_certified_local_binding(chart, tile.local_begin,
                                     tile.physical_begin_numeric,
                                     "certified tile begin");
    validate_certified_local_binding(chart, tile.local_end,
                                     tile.physical_end_numeric,
                                     "certified tile end");
    arm.certified_tiles.push_back(std::move(tile));
  }
}

void attach_rational_arm_geometry(RetainedArmPlan& arm) {
  json::array planning;
  json::array charts;
  planning.reserve(arm.charts.size());
  charts.reserve(arm.charts.size());
  for (std::size_t index = 0; index < arm.charts.size(); ++index) {
    const auto& binding = arm.charts[index];
    planning.push_back(json::object{
        {"center_exact", binding.geometry.center.str()},
        {"scale_exact", binding.geometry.scale.str()},
        {"radius_exact", binding.geometry.radius.str()},
        {"certificate_identity", binding.planning_certificate_identity}});
    charts.push_back(json::object{
        {"index", index},
        {"center", certified_rational_scalar(binding.geometry.center)},
        {"scale", certified_rational_scalar(binding.geometry.scale)},
        {"radius", certified_rational_scalar(binding.geometry.radius)}});
  }
  json::array matches;
  matches.reserve(arm.exact.matches.size());
  for (std::size_t index = 0; index < arm.exact.matches.size(); ++index) {
    const auto& match = arm.exact.matches[index];
    matches.push_back(json::object{
        {"index", index},
        {"physical", certified_rational_scalar(match.physical)},
        {"producing_local",
         certified_rational_scalar(match.producing_local)},
        {"receiving_local",
         certified_rational_scalar(match.receiving_local)}});
  }
  json::array tiles;
  tiles.reserve(arm.exact.tiles.size());
  for (std::size_t index = 0; index < arm.exact.tiles.size(); ++index) {
    const auto& tile = arm.exact.tiles[index];
    tiles.push_back(json::object{
        {"index", index}, {"chart_index", tile.chart},
        {"physical_begin",
         certified_rational_scalar(tile.physical_begin)},
        {"physical_end", certified_rational_scalar(tile.physical_end)},
        {"local_begin", certified_rational_scalar(tile.local_begin)},
        {"local_end", certified_rational_scalar(tile.local_end)}});
  }
  json::object certified{
      {"schema", "diffexp3-wolfram-certified-algebraic-arm-v1"},
      {"exact_identity",
       "legacy-rational-arm:" + arm.exact.from.str() + ":" +
           arm.exact.to.str()},
      {"charts", std::move(charts)}, {"matches", std::move(matches)},
      {"tiles", std::move(tiles)}};
  json::object raw_arm{{"planning_charts", std::move(planning)},
                       {"certified_geometry", std::move(certified)}};
  attach_certified_arm_geometry(raw_arm, arm);
}

std::shared_ptr<StoredTilePlan> build_tile_plan(
    const std::string& handle, const json::object& request,
    const std::vector<RetainedPlanChartBinding::Owner>& lower_charts,
    const std::vector<RetainedPlanChartBinding::Owner>& upper_charts) {
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  const auto division_order = as_u32(
      request.at("division_order"), "native tile division order");
  auto [lower_request, lower_bindings] = parse_retained_arm_request(
      as_object(request.at("lower"), "lower native tile arm"), lower_charts);
  auto [upper_request, upper_bindings] = parse_retained_arm_request(
      as_object(request.at("upper"), "upper native tile arm"), upper_charts);
  if (lower_bindings.front().handle != upper_bindings.front().handle)
    throw std::invalid_argument(
        "independent native tile arms must share one retained anchor chart");

  ExactPathPlanOptions options;
  options.division_order = division_order;
  const auto started = std::chrono::steady_clock::now();
  auto exact = plan_exact_independent_arms(
      lower_request, upper_request, options);
  RetainedArmPlan lower{std::move(exact.lower),
                        std::move(lower_bindings)};
  RetainedArmPlan upper{std::move(exact.upper),
                        std::move(upper_bindings)};
  const auto& lower_raw = as_object(request.at("lower"),
                                    "lower native tile arm");
  const auto& upper_raw = as_object(request.at("upper"),
                                    "upper native tile arm");
  if (lower_raw.if_contains("certified_geometry"))
    attach_certified_arm_geometry(lower_raw, lower);
  else
    attach_rational_arm_geometry(lower);
  if (upper_raw.if_contains("certified_geometry"))
    attach_certified_arm_geometry(upper_raw, upper);
  else
    attach_rational_arm_geometry(upper);
  json::object provenance{
      {"schema", "diffexp3-retained-exact-independent-arm-tile-plan-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"division_order", division_order},
      {"lower", encode_retained_arm(lower)},
      {"upper", encode_retained_arm(upper)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredTilePlan>(
      handle, checkpoint_identity, provenance_identity, division_order,
      std::move(lower), std::move(upper), elapsed);
}

std::shared_ptr<StoredTilePlan> build_single_arm_tile_plan(
    const std::string& handle, const json::object& request,
    const std::vector<RetainedPlanChartBinding::Owner>& charts) {
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  if (checkpoint_identity.empty())
    throw std::invalid_argument(
        "native single-arm tile-plan checkpoint identity cannot be empty");
  const auto division_order = as_u32(
      request.at("division_order"), "native single-arm tile division order");
  auto [arm_request, bindings] = parse_retained_arm_request(
      as_object(request.at("arm"), "native single tile arm"), charts);
  ExactPathPlanOptions options;
  options.division_order = division_order;
  const auto started = std::chrono::steady_clock::now();
  auto exact = plan_exact_arm(arm_request, options);
  const std::string arm_name = exact.direction < 0 ? "lower" : "upper";
  RetainedArmPlan retained{std::move(exact), std::move(bindings)};
  const auto& raw_arm = as_object(request.at("arm"),
                                  "single native tile arm");
  if (raw_arm.if_contains("certified_geometry"))
    attach_certified_arm_geometry(raw_arm, retained);
  else
    attach_rational_arm_geometry(retained);
  json::object provenance{
      {"schema", kRetainedSingleArmTilePlanProvenanceSchema},
      {"checkpoint_identity", checkpoint_identity},
      {"division_order", division_order},
      {"arm_name", arm_name},
      {"arm", encode_retained_arm(retained)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredTilePlan>(
      handle, checkpoint_identity, provenance_identity, division_order,
      arm_name, std::move(retained), elapsed);
}

json::value optional_plan_rim_json(
    const std::optional<std::int32_t>& rim) {
  return rim.has_value() ? json::value(*rim) : json::value(nullptr);
}

struct ResolvedPlannedEndpointBinding {
  json::object source;
  std::int32_t approach_direction = 0;
  std::optional<std::int32_t> rim;
};

ResolvedPlannedEndpointBinding resolve_planned_endpoint_binding(
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::string& arm_name,
    const std::shared_ptr<StoredLocalBase>& local) {
  if (!plan || !local)
    throw std::invalid_argument(
        "plan-bound endpoint evaluation requires retained plan and local owners");
  const auto& arm = plan->arm(arm_name);
  if (arm.exact.tiles.empty() || arm.charts.empty())
    throw std::invalid_argument(
        "plan-bound endpoint arm has no final tile/chart");
  const auto final_tile_index = arm.exact.tiles.size() - 1;
  const auto& final_tile = arm.exact.tiles.back();
  if (final_tile.chart >= arm.charts.size())
    throw std::logic_error(
        "plan-bound endpoint final tile has an invalid chart index");
  const auto& final_chart = arm.charts[final_tile.chart];
  if (!(final_tile.physical_end == arm.exact.to) ||
      !(final_chart.geometry.center == arm.exact.to) ||
      !final_tile.local_end.is_zero())
    throw std::invalid_argument(
        "plan-bound endpoint requires the final tile to end at the exact center of its retained final chart");
  if (arm.exact.direction != -1 && arm.exact.direction != 1)
    throw std::logic_error(
        "plan-bound endpoint arm has an invalid exact direction");
  if (final_chart.geometry.scale.is_zero())
    throw std::invalid_argument(
        "plan-bound endpoint final chart has a zero exact scale");
  if (local->source_chart() != final_chart.handle)
    throw std::invalid_argument(
        "plan-bound endpoint local does not name the retained final chart");
  local->require_exact_plan_binding(
      final_chart.local_geometry, final_chart.prescriptions,
      "plan-bound endpoint final local");

  ResolvedPlannedEndpointBinding resolved;
  resolved.approach_direction =
      -arm.exact.direction * final_chart.scale_sign;
  resolved.rim = exact_plan_rim(
      final_chart.prescriptions, final_chart.geometry.scale);
  resolved.source = json::object{
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
      {"tile_plan_provenance_identity", plan->provenance_identity()},
      {"arm", arm_name},
      {"endpoint_exact", arm.exact.to.str()},
      {"direction", arm.exact.direction},
      {"final_tile", final_tile_index},
      {"final_chart_index", final_tile.chart},
      {"final_chart", final_chart.handle},
      {"final_chart_identity", final_chart.exact_identity},
      {"local", local->handle()},
      {"chart", local->source_chart()},
      {"source_operator_identity", local->source_operator_identity()},
      {"checkpoint_identity", local->checkpoint_identity()},
      {"coefficient_domain", local->scalar_domain()},
      {"prescriptions", encode_plan_prescriptions(
           final_chart.prescriptions)}};
  return resolved;
}

json::object planned_endpoint_provenance(
    const std::string& checkpoint_identity,
    const ResolvedPlannedEndpointBinding& binding,
    const std::string& cancellation_mode,
    const json::object& analytic_metadata) {
  return json::object{
      {"schema",
       "diffexp3-retained-native-plan-bound-endpoint-sector-limit-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"source", binding.source},
      {"approach_direction", binding.approach_direction},
      {"rim", optional_plan_rim_json(binding.rim)},
      {"cancellation", json::object{{"mode", cancellation_mode}}},
      {"analytic_metadata", analytic_metadata}};
}

std::shared_ptr<StoredEndpointResult> build_planned_endpoint_limit(
    const std::string& endpoint_handle, const json::object& request,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& local) {
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  const auto expected_plan_checkpoint = required_string(
      request, "tile_plan_checkpoint_identity");
  const auto expected_source_checkpoint = required_string(
      request, "source_checkpoint_identity");
  if (checkpoint_identity.empty() || expected_plan_checkpoint.empty() ||
      expected_source_checkpoint.empty())
    throw std::invalid_argument(
        "plan-bound endpoint checkpoint identities must be nonempty");
  if (expected_plan_checkpoint != plan->checkpoint_identity())
    throw std::invalid_argument(
        "plan-bound endpoint tile-plan checkpoint identity is stale or mismatched");
  if (expected_source_checkpoint != local->checkpoint_identity())
    throw std::invalid_argument(
        "plan-bound endpoint source checkpoint identity is stale or mismatched");
  const auto arm_name = required_string(request, "arm");
  const auto binding = resolve_planned_endpoint_binding(
      plan, arm_name, local);
  const auto cancellation = parse_endpoint_cancellation_policy(request);
  auto analytic_metadata = local->exact_analytic_metadata();
  const auto provenance = planned_endpoint_provenance(
      checkpoint_identity, binding, cancellation.cancellation_mode,
      analytic_metadata);
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));

  EndpointLimitOptions options;
  options.approach_direction = binding.approach_direction;
  options.imaginary_sign = binding.rim;
  options.allow_certified_numeric_cancellation =
      cancellation.allow_certified_numeric_cancellation;
  const auto started = std::chrono::steady_clock::now();
  auto result = local->endpoint_limit(options);
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredEndpointResult>(
      endpoint_handle, checkpoint_identity, provenance_identity,
      local->handle(), local->source_chart(),
      local->source_operator_identity(), expected_source_checkpoint,
      local->scalar_domain(), binding.approach_direction, std::nullopt,
      cancellation.cancellation_mode, std::move(analytic_metadata),
      std::move(result), elapsed, binding.source, binding.rim,
      plan, local);
}

json::object planned_match_handoff_record(
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::string& arm_name, std::size_t match_index,
    const std::vector<std::string>& basis_handles,
    const std::vector<std::shared_ptr<StoredLocalBase>>& basis,
    const std::string& incoming_handle,
    const std::shared_ptr<StoredLocalBase>& incoming,
    const std::string& result_checkpoint,
    const json::object& native_summary,
    bool compact_plan_reference = false) {
  if (!plan || !incoming || basis.empty() ||
      basis.size() != basis_handles.size())
    throw std::invalid_argument(
        "planned match handoff record requires its complete owner set");
  const auto& arm = plan->arm(arm_name);
  if (match_index >= arm.exact.matches.size())
    throw std::invalid_argument(
        "planned match handoff index is outside its retained arm");
  const auto& exact_match = arm.exact.matches[match_index];
  const auto& producing = arm.charts.at(exact_match.producing_chart);
  const auto& receiving = arm.charts.at(exact_match.receiving_chart);
  const auto producing_rim = exact_plan_rim(
      producing.prescriptions, producing.geometry.scale);
  const auto receiving_rim = exact_plan_rim(
      receiving.prescriptions, receiving.geometry.scale);
  json::array basis_sources;
  basis_sources.reserve(basis.size());
  for (std::size_t column = 0; column < basis.size(); ++column) {
    json::object source{
        {"column", column}, {"local", basis_handles[column]},
        {"checkpoint_identity", basis[column]->checkpoint_identity()}};
    if (compact_plan_reference)
      source["source_operator_reference"] =
          basis[column]->source_operator_reference();
    else
      source["source_operator_identity"] =
          basis[column]->source_operator_identity();
    basis_sources.push_back(std::move(source));
  }
  if (compact_plan_reference)
    return json::object{
        {"schema", "diffexp3-retained-exact-plan-match-hop-v3"},
        {"tile_plan", plan->handle()},
        {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
        {"arm", arm_name}, {"match", match_index},
        {"geometry", encode_plan_match(arm, match_index)},
        {"producing", json::object{
             {"tile", match_index},
             {"chart", producing.handle},
             {"chart_identity_reference",
              compact_matching_identity_reference(
                  producing.exact_identity)},
             {"local_point_exact", exact_match.producing_local.str()},
             {"effective_rim", optional_plan_rim_json(producing_rim)},
             {"prescriptions",
              encode_plan_prescriptions(producing.prescriptions)},
             {"incoming", json::object{
                  {"local", incoming_handle},
                  {"checkpoint_identity", incoming->checkpoint_identity()},
                  {"source_operator_reference",
                   incoming->source_operator_reference()}}}}},
        {"receiving", json::object{
             {"tile", match_index + 1},
             {"chart", receiving.handle},
             {"chart_identity_reference",
              compact_matching_identity_reference(
                  receiving.exact_identity)},
             {"local_point_exact", exact_match.receiving_local.str()},
             {"effective_rim", optional_plan_rim_json(receiving_rim)},
             {"prescriptions",
              encode_plan_prescriptions(receiving.prescriptions)},
             {"basis", std::move(basis_sources)}}},
        {"result_checkpoint_identity", result_checkpoint},
        {"advance", json::object{
             {"scope", "single-match-handoff"},
             {"state", "retained-receiving-basis-weights"},
             {"source_tile", match_index},
             {"receiving_tile", match_index + 1},
             {"whole_arm_complete", false}}}};
  return json::object{
      {"schema", "diffexp3-retained-exact-plan-match-hop-v1"},
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
      {"tile_plan_provenance_identity", plan->provenance_identity()},
      {"arm", arm_name}, {"match", match_index},
      {"geometry", encode_plan_match(arm, match_index)},
      {"producing", json::object{
           {"tile", match_index},
           {"chart", producing.handle},
           {"chart_identity", producing.exact_identity},
           {"local_point_exact", exact_match.producing_local.str()},
           {"effective_rim", optional_plan_rim_json(producing_rim)},
           {"prescriptions",
            encode_plan_prescriptions(producing.prescriptions)},
           {"incoming", json::object{
                {"local", incoming_handle},
                {"checkpoint_identity", incoming->checkpoint_identity()},
                {"source_operator_identity",
                 incoming->source_operator_identity()}}}}},
      {"receiving", json::object{
           {"tile", match_index + 1},
           {"chart", receiving.handle},
           {"chart_identity", receiving.exact_identity},
           {"local_point_exact", exact_match.receiving_local.str()},
           {"effective_rim", optional_plan_rim_json(receiving_rim)},
           {"prescriptions",
            encode_plan_prescriptions(receiving.prescriptions)},
           {"basis", std::move(basis_sources)}}},
      {"result_checkpoint_identity", result_checkpoint},
      {"native_match_provenance_identity",
       required_string(native_summary, "provenance_identity")},
      {"advance", json::object{
           {"scope", "single-match-handoff"},
           {"state", "retained-receiving-basis-weights"},
           {"source_tile", match_index},
           {"receiving_tile", match_index + 1},
           {"whole_arm_complete", false}}}};
}

struct NativeAcbSaturationBinding {
  std::string request_key;
  json::object request;
};

NativeAcbSaturationBinding native_acb_saturation_binding(
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::string& session_configuration_identity,
    const std::string& arm,
    std::size_t match_index, const std::string& match_checkpoint_identity,
    bool compact_plan_reference = false);

std::shared_ptr<StoredPlannedMatchHop> build_planned_match_hop(
    const std::string& match_handle, const json::object& request,
    const std::string& domain, slong precision_bits,
    const std::string& active_session_configuration_identity,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::vector<std::string>& basis_handles,
    const std::vector<std::shared_ptr<StoredLocalBase>>& basis,
    const std::string& incoming_handle,
    const std::shared_ptr<StoredLocalBase>& incoming,
    bool compact_plan_reference = false,
    bool require_normalized_singular_frame = false,
    const std::optional<json::object>&
        retained_singular_request = std::nullopt) {
  const auto started = std::chrono::steady_clock::now();
  const auto arm_name = required_string(request, "arm");
  const auto match_index = static_cast<std::size_t>(
      as_u64(request.at("match"), "planned local match index"));
  const auto& arm = plan->arm(arm_name);
  if (match_index >= arm.exact.matches.size())
    throw std::invalid_argument(
        "planned local match index is out of range");
  const auto& exact_match = arm.exact.matches[match_index];
  const auto& producing = arm.charts.at(exact_match.producing_chart);
  const auto& receiving = arm.charts.at(exact_match.receiving_chart);
  if (basis.empty())
    throw std::invalid_argument("planned local match basis cannot be empty");

  // The plan, rather than the caller, binds every coordinate, chart,
  // prescription, rim and source checkpoint passed to the existing matching
  // kernels.  Locals must reproduce the prepared chart snapshot exactly.
  incoming->require_exact_plan_binding(
      producing.local_geometry, producing.prescriptions,
      "planned incoming " + incoming_handle);
  for (std::size_t column = 0; column < basis.size(); ++column)
    basis[column]->require_exact_plan_binding(
      receiving.local_geometry, receiving.prescriptions,
        "planned basis " + basis_handles[column]);

  json::array basis_checkpoints;
  basis_checkpoints.reserve(basis.size());
  for (const auto& local : basis)
    basis_checkpoints.emplace_back(local->checkpoint_identity());
  const auto result_checkpoint = required_string(
      request, "checkpoint_identity");
  if (result_checkpoint.empty())
    throw std::invalid_argument(
        "planned local match checkpoint identity cannot be empty");

  const auto& certified_match = arm.certified_matches.at(match_index);
  const auto encode_point = [&](const CertifiedPlanPoint& point) {
    return json::object{
        {"exact", point.local.exact_coordinate},
        {"value", point.local_numeric_encoding},
        {"sign", point.local.sign}};
  };
  json::object kernel_request{
      {"basis", [&]() {
         json::array values;
         for (const auto& handle : basis_handles) values.emplace_back(handle);
         return values;
       }()},
      {"incoming", incoming_handle},
      {"basis_chart", receiving.handle},
      {"incoming_chart", producing.handle},
      {"basis_point", encode_point(certified_match.receiving)},
      {"incoming_point", encode_point(certified_match.producing)},
      {"certified_physical_match_point_exact",
       certified_match.receiving.physical_exact},
      {"epsilon", request.at("epsilon")},
      {"basis_checkpoint_identities", std::move(basis_checkpoints)},
      {"incoming_checkpoint_identity", incoming->checkpoint_identity()},
      {"checkpoint_identity", result_checkpoint}};

  const auto producing_rim = exact_plan_rim(
      producing.prescriptions, producing.geometry.scale);
  const auto receiving_rim = exact_plan_rim(
      receiving.prescriptions, receiving.geometry.scale);
  std::shared_ptr<StoredMatchBase> native_match;
  std::optional<json::object> expected_singular_request;
  if (domain == "rational") {
    native_match = build_exact_regular_match(
        match_handle, kernel_request, basis_handles, basis, incoming_handle,
        incoming);
  } else if (domain == "acb") {
    kernel_request["basis_imaginary_sign"] =
        optional_plan_rim_json(receiving_rim);
    kernel_request["incoming_imaginary_sign"] =
        optional_plan_rim_json(producing_rim);
    kernel_request["refinement"] = request.at("refinement");
    if (const auto* exact = request.if_contains("exact_lattice"))
      kernel_request["exact_lattice"] = *exact;
    else if (const auto* native =
                 request.if_contains("native_unit_saturation"))
      kernel_request["native_unit_saturation"] = *native;
    else if (const auto* singular =
                 request.if_contains("native_singular_scc_saturation")) {
      const auto expected = retained_singular_request.has_value()
          ? NativeAcbSaturationBinding{
                "native_singular_scc_saturation",
                *retained_singular_request}
          : native_acb_saturation_binding(
                plan, active_session_configuration_identity, arm_name,
                match_index, result_checkpoint,
                compact_plan_reference);
      if (expected.request_key !=
              "native_singular_scc_saturation" ||
          json::serialize(canonical_json_value(*singular)) !=
              json::serialize(
                  canonical_json_value(expected.request)))
        throw std::invalid_argument(
            "planned singular-SCC Acb saturation request does not match the retained receiving SCC");
      kernel_request["native_singular_scc_saturation"] = *singular;
      expected_singular_request = std::move(expected.request);
    } else
      throw std::invalid_argument(
          "planned Acb matching requires an exact lattice, ordinary native unit-leading request, or singular-SCC valuation-zero request");
    native_match = build_refined_acb_match(
        match_handle, kernel_request, basis_handles, basis, incoming_handle,
        incoming, precision_bits, active_session_configuration_identity,
        expected_singular_request, require_normalized_singular_frame);
  } else {
    throw std::invalid_argument(
        "plan-driven local matching requires rational or Acb coefficients");
  }

  const auto native_summary = native_match->summary();
  if (required_string(native_summary, "physical_match_point_exact") !=
      certified_match.receiving.physical_exact)
    throw std::logic_error(
        "native local match physical point differs from its retained certified geometry");

  auto handoff = planned_match_handoff_record(
      plan, arm_name, match_index, basis_handles, basis, incoming_handle,
      incoming, result_checkpoint, native_summary,
      compact_plan_reference);
  const auto provenance_identity = json::serialize(
      canonical_json_value(handoff));
  const auto elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredPlannedMatchHop>(
      std::move(native_match), result_checkpoint, provenance_identity,
      std::move(handoff), elapsed_ms, plan, basis, incoming);
}

using BoundedDivergentCancellation =
    StoredLineIntegrationOptions::BoundedDivergentCancellation;

BoundedDivergentCancellation parse_bounded_divergent_cancellation(
    const json::value& raw, const char* label) {
  const auto& policy = as_object(raw, label);
  require_exact_keys(policy,
                     {"mode", "relative_tolerance", "provenance"},
                     label);
  if (required_string(policy, "mode") != "bounded-relative-acb")
    throw std::invalid_argument(
        std::string(label) +
        " mode must be exactly bounded-relative-acb");
  const auto tolerance_text = required_string(
      policy, "relative_tolerance");
  const auto provenance = required_string(policy, "provenance");
  if (tolerance_text.empty() || provenance.empty())
    throw std::invalid_argument(
        std::string(label) +
        " requires a nonempty tolerance and producer provenance");
  auto tolerance = Magnitude::decimal(tolerance_text);
  if (!tolerance.is_finite() || tolerance.is_zero() ||
      Magnitude::one() <= tolerance)
    throw std::invalid_argument(
        std::string(label) +
        " tolerance must be finite and strictly between zero and one");
  return {std::move(tolerance), tolerance_text, provenance};
}

json::object encode_bounded_divergent_cancellation(
    const BoundedDivergentCancellation& policy) {
  return json::object{
      {"mode", "bounded-relative-acb"},
      {"relative_tolerance", policy.relative_tolerance_text},
      {"provenance", policy.provenance}};
}

std::shared_ptr<StoredLineResult> build_planned_line_result(
    const std::string& handle, const json::object& request,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& local) {
  const auto checkpoint_identity = required_string(
      request, "checkpoint_identity");
  const auto source_checkpoint = required_string(
      request, "source_checkpoint_identity");
  if (source_checkpoint != local->checkpoint_identity())
    throw std::invalid_argument(
        "planned line source checkpoint identity differs from its retained local");
  const auto expected_plan_checkpoint = required_string(
      request, "tile_plan_checkpoint_identity");
  if (expected_plan_checkpoint != plan->checkpoint_identity())
    throw std::invalid_argument(
        "planned line tile-plan checkpoint identity differs from retained state");
  const auto arm_name = required_string(request, "arm");
  const auto tile_index = static_cast<std::size_t>(
      as_u64(request.at("tile"), "planned line tile index"));
  const auto& arm = plan->arm(arm_name);
  if (tile_index >= arm.exact.tiles.size())
    throw std::invalid_argument("planned line tile index is out of range");
  const auto& tile = arm.exact.tiles[tile_index];
  const auto& certified_tile = arm.certified_tiles.at(tile_index);
  const auto& binding = arm.charts.at(tile.chart);
  if (local->source_chart() != binding.handle)
    throw std::invalid_argument(
        "planned line local does not belong to the tile's retained chart");
  const auto& epsilon = as_object(request.at("epsilon"),
                                  "planned line epsilon window");
  const auto delivered = parse_epsilon_window(
      epsilon, "planned line epsilon window");
  const bool certify_tail =
      request.if_contains("certify_tail") != nullptr &&
      request.at("certify_tail").as_bool();
  std::optional<BoundedDivergentCancellation> divergent_cancellation;
  if (const auto* policy = request.if_contains("divergent_cancellation"))
    divergent_cancellation = parse_bounded_divergent_cancellation(
        *policy, "planned line divergent-cancellation policy");
  auto interval = encode_plan_tile(arm, tile_index);
  const auto rim = exact_plan_rim(
      binding.prescriptions, binding.geometry.scale);
  const bool certified_zero =
      certified_tile_zero_length(arm, tile_index);
  const auto started = std::chrono::steady_clock::now();
  StoredLineIntegral result;
  try {
    result = local->integrate_planned_line(
        binding.geometry, binding.local_geometry, binding.scale_numeric,
        binding.local_geometry.scale_exact, binding.scale_sign,
        binding.prescriptions,
        certified_tile.local_begin, certified_tile.local_end,
        tile.local_end < tile.local_begin, certified_zero,
        delivered, rim, certify_tail,
        divergent_cancellation);
  } catch (const NativeIntegrationError& error) {
    std::ostringstream detail;
    detail << error.what() << "; arm=" << arm_name
           << "; tile=" << tile_index
           << "; physical_interval=[" << tile.physical_begin.str()
           << "," << tile.physical_end.str() << "]"
           << "; local_interval=[" << tile.local_begin.str()
           << "," << tile.local_end.str() << "]"
           << "; chart=" << binding.handle
           << "; chart_center=" << binding.geometry.center.str()
           << "; chart_scale=" << binding.geometry.scale.str()
           << "; rim="
           << (rim.has_value() ? std::to_string(*rim) : "none")
           << "; source_local=" << local->handle();
    NativeIntegrationError contextual(error.code, error.id, detail.str());
    contextual.absolute_power = error.absolute_power;
    contextual.log_power = error.log_power;
    contextual.epsilon_power = error.epsilon_power;
    contextual.component = error.component;
    throw contextual;
  }
  const auto elapsed = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  json::object provenance{
      {"schema",
       "diffexp3-retained-native-physical-tile-integral-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"tile_plan", plan->handle()},
      {"tile_plan_checkpoint_identity", expected_plan_checkpoint},
      {"arm", arm_name}, {"tile", tile_index},
      {"interval", interval},
      {"source", json::object{
           {"local", local->handle()}, {"chart", local->source_chart()},
           {"checkpoint_identity", source_checkpoint}}},
      {"epsilon", json::object{{"min", delivered.min_power},
                                {"max", delivered.complete_max}}},
      {"tail_certificate_requested", certify_tail},
      {"tail_certificate_status",
       result.diagnostics.tail_certificate_status},
      {"tail_witness_radius_exact",
       result.diagnostics.tail_witness_radius_exact.empty()
           ? json::value(nullptr)
           : json::value(result.diagnostics.tail_witness_radius_exact)},
      {"divergent_cancellation",
       divergent_cancellation.has_value()
           ? json::value(encode_bounded_divergent_cancellation(
                 *divergent_cancellation))
           : json::value(json::object{{"mode", "exact-singleton"}})},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)},
      {"error_provenance", result.value.error.provenance}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  return std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, arm_name, tile_index,
      std::move(interval), source_checkpoint, std::move(result), elapsed,
      plan, local);
}

std::size_t checked_diagnostic_sum(std::size_t left, std::size_t right,
                                   const char* label) {
  if (right > std::numeric_limits<std::size_t>::max() - left)
    throw std::overflow_error(std::string(label) + " overflow");
  return left + right;
}

StoredLineIntegral aggregate_retained_lines(
    const std::vector<std::shared_ptr<StoredLineResult>>& components,
    const std::vector<std::int32_t>& signs, const std::string& detail) {
  if (components.empty() || components.size() != signs.size())
    throw std::invalid_argument(
        "native line aggregation requires a nonempty signed component list");
  for (const auto sign : signs)
    if (sign != -1 && sign != 1)
      throw std::invalid_argument(
          "native line aggregation signs must be +1 or -1");

  auto epsilon_min = components.front()->result().value.epsilon.min_power;
  auto epsilon_max =
      components.front()->result().value.epsilon.complete_max;
  const auto dimension = components.front()->result().value.dimension;
  for (const auto& component : components) {
    const auto& value = component->result().value;
    if (value.dimension != dimension)
      throw std::invalid_argument(
          "native line aggregate component dimensions differ");
    // A finite Laurent frame's lower edge is an exact structural bound:
    // powers below it are zero, not unknown.  Aggregation therefore takes
    // the union lower edge while still intersecting complete upper edges.
    epsilon_min = std::min(epsilon_min, value.epsilon.min_power);
    epsilon_max = std::min(epsilon_max, value.epsilon.complete_max);
  }
  if (epsilon_min > epsilon_max)
    throw std::domain_error(
        "native line aggregate components have no common complete epsilon window");

  StoredLineIntegral result;
  result.value.epsilon = {epsilon_min, epsilon_max};
  result.value.dimension = dimension;
  result.value.coefficients.reserve(
      result.value.epsilon.width() * dimension);
  for (std::int64_t raw_power = epsilon_min; raw_power <= epsilon_max;
       ++raw_power) {
    const auto power = static_cast<std::int32_t>(raw_power);
    for (std::uint32_t component_index = 0;
         component_index < dimension; ++component_index) {
      ComplexBall sum(0);
      for (std::size_t index = 0; index < components.size(); ++index) {
        const auto& value = components[index]->result().value;
        if (power < value.epsilon.min_power) continue;
        const auto& coefficient = value.at(power, component_index);
        sum += signs[index] == 1 ? coefficient : -coefficient;
      }
      result.value.coefficients.push_back(std::move(sum));
    }
  }

  bool all_error_envelopes = true;
  bool all_certified_errors = true;
  bool any_advisory_error = false;
  json::array error_sources;
  for (std::size_t index = 0; index < components.size(); ++index) {
    const auto& error = components[index]->result().value.error;
    if (error.empty() || error.frame.complete_max < epsilon_max) {
      all_error_envelopes = false;
      break;
    }
    all_certified_errors &= error.guarantee == ErrorGuarantee::Certified;
    any_advisory_error |= error.guarantee == ErrorGuarantee::Advisory;
    error_sources.push_back(json::object{
        {"sign", signs[index]},
        {"guarantee", error_guarantee_name(error.guarantee)},
        {"provenance", error.provenance}});
  }
  if (all_error_envelopes) {
    result.value.error.frame = result.value.epsilon;
    result.value.error.guarantee = all_certified_errors
        ? ErrorGuarantee::Certified
        : any_advisory_error ? ErrorGuarantee::Advisory
                             : ErrorGuarantee::None;
    result.value.error.absolute.reserve(result.value.epsilon.width());
    for (std::int64_t raw_power = epsilon_min; raw_power <= epsilon_max;
         ++raw_power) {
      const auto power = static_cast<std::int32_t>(raw_power);
      auto error_sum = Magnitude::zero();
      for (const auto& component : components) {
        const auto& error = component->result().value.error;
        if (power < error.frame.min_power) continue;
        error_sum += error.absolute.at(static_cast<std::size_t>(
            power - error.frame.min_power));
      }
      result.value.error.absolute.push_back(std::move(error_sum));
    }
    result.value.error.provenance = json::serialize(json::object{
        {"schema", "diffexp3-native-line-error-sum-v1"},
        {"components", std::move(error_sources)}});
  }

  bool all_full_local = true;
  auto& diagnostics = result.diagnostics;
  diagnostics.detail = detail;
  diagnostics.tail_certificate_requested = true;
  diagnostics.tail_certificate_status = "aggregate-certified";
  for (const auto& component : components) {
    const auto& input = component->result();
    all_full_local &=
        input.scope == LineIntegrationScope::FullLocalWithCertifiedTail;
    const auto& source = input.diagnostics;
    diagnostics.input_monomial_cells = checked_diagnostic_sum(
        diagnostics.input_monomial_cells, source.input_monomial_cells,
        "aggregate input monomial count");
    diagnostics.grouped_monomials = checked_diagnostic_sum(
        diagnostics.grouped_monomials, source.grouped_monomials,
        "aggregate grouped monomial count");
    diagnostics.zero_groups_skipped = checked_diagnostic_sum(
        diagnostics.zero_groups_skipped, source.zero_groups_skipped,
        "aggregate skipped-zero count");
    diagnostics.cancelled_divergent_groups = checked_diagnostic_sum(
        diagnostics.cancelled_divergent_groups,
        source.cancelled_divergent_groups,
        "aggregate cancelled-divergence count");
    diagnostics.bounded_cancelled_divergent_coefficients =
        checked_diagnostic_sum(
            diagnostics.bounded_cancelled_divergent_coefficients,
            source.bounded_cancelled_divergent_coefficients,
            "aggregate bounded-cancelled-divergence count");
    if (source.divergent_cancellation_mode == "bounded-relative-acb") {
      if (diagnostics.divergent_cancellation_mode == "exact-singleton") {
        diagnostics.divergent_cancellation_mode =
            source.divergent_cancellation_mode;
        diagnostics.divergent_relative_tolerance =
            source.divergent_relative_tolerance;
        diagnostics.divergent_cancellation_provenance =
            source.divergent_cancellation_provenance;
      } else if (diagnostics.divergent_relative_tolerance !=
                     source.divergent_relative_tolerance ||
                 diagnostics.divergent_cancellation_provenance !=
                     source.divergent_cancellation_provenance) {
        throw std::invalid_argument(
            "native line aggregate mixes divergent-cancellation policies");
      }
    }
    diagnostics.primitive_evaluations = checked_diagnostic_sum(
        diagnostics.primitive_evaluations, source.primitive_evaluations,
        "aggregate primitive evaluation count");
    diagnostics.primitive_component_applications = checked_diagnostic_sum(
        diagnostics.primitive_component_applications,
        source.primitive_component_applications,
        "aggregate primitive application count");
    diagnostics.primitive_component_reuses = checked_diagnostic_sum(
        diagnostics.primitive_component_reuses,
        source.primitive_component_reuses,
        "aggregate primitive reuse count");
    diagnostics.has_center_endpoint |= source.has_center_endpoint;
    diagnostics.tail_certificate_requested &=
        source.tail_certificate_requested;
  }
  if (all_full_local && result.value.error.guarantee ==
                            ErrorGuarantee::Certified) {
    result.scope = LineIntegrationScope::FullLocalWithCertifiedTail;
  } else {
    result.scope = LineIntegrationScope::StoredTruncation;
    diagnostics.tail_certificate_status = all_full_local
        ? "aggregate-missing-certified-error-envelope"
        : "aggregate-has-stored-truncation-component";
  }
  result.imaginary_sign = std::nullopt;
  return result;
}

// Incremental equivalent of aggregate_retained_lines for the all-positive
// physical-tile sum used by transport observables.  It retains only the
// accumulated epsilon vector and the small per-tile error envelopes; the
// coefficient-bearing projected locals and tile line objects can therefore
// die immediately after each tile.
class StreamingStoredLineAccumulator {
 public:
  void add(const StoredLineIntegral& input) {
    if (input.diagnostics.detail.find(
            "terminal factorized production contraction") !=
            std::string::npos ||
        input.diagnostics.detail.find(
            "terminal direct-physical diagnostic contraction") !=
            std::string::npos ||
        input.diagnostics.detail.find(
            "terminal adjoint diagnostic contraction") !=
            std::string::npos)
      terminal_contraction_detail_ = input.diagnostics.detail;
    if (!initialized_) {
      result_ = input;
      result_.imaginary_sign = std::nullopt;
      result_.diagnostics.tail_witness_radius_exact.clear();
      initialized_ = true;
    } else {
      if (input.value.dimension != result_.value.dimension)
        throw std::invalid_argument(
            "streaming native line aggregate component dimensions differ");
      const auto minimum = std::min(
          result_.value.epsilon.min_power, input.value.epsilon.min_power);
      const auto maximum = std::min(
          result_.value.epsilon.complete_max,
          input.value.epsilon.complete_max);
      if (minimum > maximum)
        throw std::domain_error(
            "streaming native line aggregate has no common complete epsilon window");
      std::vector<ComplexBall> coefficients;
      coefficients.reserve(EpsilonWindow{minimum, maximum}.width() *
                           result_.value.dimension);
      for (std::int64_t raw_power = minimum; raw_power <= maximum;
           ++raw_power) {
        const auto power = static_cast<std::int32_t>(raw_power);
        for (std::uint32_t component = 0;
             component < result_.value.dimension; ++component) {
          ComplexBall value(0);
          if (power >= result_.value.epsilon.min_power)
            value += result_.value.at(power, component);
          if (power >= input.value.epsilon.min_power)
            value += input.value.at(power, component);
          coefficients.push_back(std::move(value));
        }
      }
      result_.value.epsilon = {minimum, maximum};
      result_.value.coefficients = std::move(coefficients);
      accumulate_diagnostics(input.diagnostics);
    }
    errors_.push_back(input.value.error);
    all_full_local_ &=
        input.scope == LineIntegrationScope::FullLocalWithCertifiedTail;
  }

  StoredLineIntegral finish(const std::string& detail) {
    if (!initialized_)
      throw std::invalid_argument(
          "streaming native line aggregate has no tile values");
    finalize_errors();
    auto& diagnostics = result_.diagnostics;
    diagnostics.detail = detail;
    if (!terminal_contraction_detail_.empty())
      diagnostics.detail += "; " + terminal_contraction_detail_;
    diagnostics.tail_certificate_status = "aggregate-certified";
    if (all_full_local_ && result_.value.error.guarantee ==
                               ErrorGuarantee::Certified) {
      result_.scope = LineIntegrationScope::FullLocalWithCertifiedTail;
    } else {
      result_.scope = LineIntegrationScope::StoredTruncation;
      diagnostics.tail_certificate_status = all_full_local_
          ? "aggregate-missing-certified-error-envelope"
          : "aggregate-has-stored-truncation-component";
    }
    result_.imaginary_sign = std::nullopt;
    return std::move(result_);
  }

 private:
  std::string terminal_contraction_detail_;

  void accumulate_diagnostics(
      const StoredLineIntegrationDiagnostics& source) {
    auto& diagnostics = result_.diagnostics;
    diagnostics.input_monomial_cells = checked_diagnostic_sum(
        diagnostics.input_monomial_cells, source.input_monomial_cells,
        "streaming aggregate input monomial count");
    diagnostics.grouped_monomials = checked_diagnostic_sum(
        diagnostics.grouped_monomials, source.grouped_monomials,
        "streaming aggregate grouped monomial count");
    diagnostics.zero_groups_skipped = checked_diagnostic_sum(
        diagnostics.zero_groups_skipped, source.zero_groups_skipped,
        "streaming aggregate skipped-zero count");
    diagnostics.cancelled_divergent_groups = checked_diagnostic_sum(
        diagnostics.cancelled_divergent_groups,
        source.cancelled_divergent_groups,
        "streaming aggregate cancelled-divergence count");
    diagnostics.bounded_cancelled_divergent_coefficients =
        checked_diagnostic_sum(
            diagnostics.bounded_cancelled_divergent_coefficients,
            source.bounded_cancelled_divergent_coefficients,
            "streaming aggregate bounded-cancelled-divergence count");
    if (source.divergent_cancellation_mode == "bounded-relative-acb") {
      if (diagnostics.divergent_cancellation_mode == "exact-singleton") {
        diagnostics.divergent_cancellation_mode =
            source.divergent_cancellation_mode;
        diagnostics.divergent_relative_tolerance =
            source.divergent_relative_tolerance;
        diagnostics.divergent_cancellation_provenance =
            source.divergent_cancellation_provenance;
      } else if (diagnostics.divergent_relative_tolerance !=
                     source.divergent_relative_tolerance ||
                 diagnostics.divergent_cancellation_provenance !=
                     source.divergent_cancellation_provenance) {
        throw std::invalid_argument(
            "streaming native line aggregate mixes divergent-cancellation policies");
      }
    }
    diagnostics.primitive_evaluations = checked_diagnostic_sum(
        diagnostics.primitive_evaluations, source.primitive_evaluations,
        "streaming aggregate primitive evaluation count");
    diagnostics.primitive_component_applications = checked_diagnostic_sum(
        diagnostics.primitive_component_applications,
        source.primitive_component_applications,
        "streaming aggregate primitive application count");
    diagnostics.primitive_component_reuses = checked_diagnostic_sum(
        diagnostics.primitive_component_reuses,
        source.primitive_component_reuses,
        "streaming aggregate primitive reuse count");
    diagnostics.has_center_endpoint |= source.has_center_endpoint;
    diagnostics.tail_certificate_requested &=
        source.tail_certificate_requested;
  }

  void finalize_errors() {
    const auto epsilon = result_.value.epsilon;
    bool all_error_envelopes = true;
    bool all_certified_errors = true;
    bool any_advisory_error = false;
    json::array error_sources;
    error_sources.reserve(errors_.size());
    for (const auto& error : errors_) {
      if (error.empty() || error.frame.complete_max < epsilon.complete_max) {
        all_error_envelopes = false;
        break;
      }
      all_certified_errors &=
          error.guarantee == ErrorGuarantee::Certified;
      any_advisory_error |=
          error.guarantee == ErrorGuarantee::Advisory;
      error_sources.push_back(json::object{
          {"sign", 1},
          {"guarantee", error_guarantee_name(error.guarantee)},
          {"provenance", error.provenance}});
    }
    result_.value.error = ErrorEnvelope{};
    if (!all_error_envelopes) return;
    auto& output = result_.value.error;
    output.frame = epsilon;
    output.guarantee = all_certified_errors
        ? ErrorGuarantee::Certified
        : any_advisory_error ? ErrorGuarantee::Advisory
                             : ErrorGuarantee::None;
    output.absolute.reserve(epsilon.width());
    for (std::int64_t raw_power = epsilon.min_power;
         raw_power <= epsilon.complete_max; ++raw_power) {
      const auto power = static_cast<std::int32_t>(raw_power);
      auto sum = Magnitude::zero();
      for (const auto& error : errors_) {
        if (power < error.frame.min_power) continue;
        sum += error.absolute.at(static_cast<std::size_t>(
            power - error.frame.min_power));
      }
      output.absolute.push_back(std::move(sum));
    }
    output.provenance = json::serialize(json::object{
        {"schema", "diffexp3-native-line-error-sum-v1"},
        {"components", std::move(error_sources)}});
  }

  bool initialized_ = false;
  bool all_full_local_ = true;
  StoredLineIntegral result_;
  std::vector<ErrorEnvelope> errors_;
};

std::vector<std::shared_ptr<StoredLocalBase>> unique_line_local_owners(
    const std::vector<std::shared_ptr<StoredLocalBase>>& owners) {
  std::set<std::string> seen;
  std::vector<std::shared_ptr<StoredLocalBase>> unique;
  unique.reserve(owners.size());
  for (const auto& owner : owners) {
    if (!owner)
      throw std::logic_error("native line aggregate lost a local owner");
    if (seen.insert(owner->handle()).second) unique.push_back(owner);
  }
  return unique;
}

struct RetainedLocalFrameContract {
  EpsilonWindow epsilon;
  std::int32_t top_valid = kCompleteInfinity;
  std::uint32_t dimension = 0;
  std::uint32_t taylor_complete_max = 0;
};

struct WholeArmEpsilonContract {
  EpsilonWindow work;
  std::int32_t public_required_complete_max = 0;
  std::int32_t match_required_complete_max = 0;
};

struct ObservableEpsilonContract {
  EpsilonWindow requested;
  std::int32_t required_complete_max = 0;
};

ObservableEpsilonContract parse_observable_epsilon_contract(
    const json::value& raw, const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(object, {"min", "max", "required_complete_max"},
                     label);
  ObservableEpsilonContract result{
      {as_i32(object.at("min"), "observable epsilon minimum"),
       as_i32(object.at("max"), "observable epsilon maximum")},
      as_i32(object.at("required_complete_max"),
             "observable required epsilon maximum")};
  (void)result.requested.width();
  if (result.required_complete_max < result.requested.min_power ||
      result.required_complete_max > result.requested.complete_max)
    throw std::invalid_argument(
        "observable required epsilon maximum must lie in its exact output window");
  return result;
}

EndpointLimitResult restrict_endpoint_result_epsilon(
    EndpointLimitResult result,
    const ObservableEpsilonContract& contract,
    const char* label) {
  const auto available = endpoint_value_window(result);
  const auto complete_max = std::min(
      contract.requested.complete_max, available.complete_max);
  if (contract.requested.min_power > complete_max ||
      contract.required_complete_max > complete_max)
    throw std::domain_error(
        std::string(label) +
        " does not cover its required exact epsilon output contract");
  std::vector<EpsilonFrame<ComplexBall>> restricted;
  restricted.reserve(result.values.size());
  for (const auto& component : result.values) {
    std::vector<ComplexBall> coefficients;
    coefficients.reserve(static_cast<std::size_t>(
        static_cast<std::int64_t>(complete_max) -
        contract.requested.min_power + 1));
    for (std::int64_t raw_power = contract.requested.min_power;
         raw_power <= complete_max; ++raw_power)
      coefficients.push_back(component.coefficient(
          static_cast<std::int32_t>(raw_power)));
    restricted.emplace_back(
        EpsilonWindow{contract.requested.min_power, complete_max},
        std::move(coefficients));
  }
  result.values = std::move(restricted);
  return result;
}

EndpointLimitResult endpoint_result_from_retained_evaluation(
    const LocalEvaluation& evaluation,
    std::optional<std::int32_t> derived_rim) {
  if (evaluation.value.dimension != 1)
    throw std::logic_error(
        "transport endpoint row evaluation did not remain scalar");
  if (!evaluation.value.error.empty())
    throw std::domain_error(
        "transport endpoint regular-point evaluation produced an error envelope whose endpoint propagation is not implemented");
  std::vector<ComplexBall> coefficients;
  coefficients.reserve(evaluation.value.epsilon.width());
  for (std::int64_t raw_power = evaluation.value.epsilon.min_power;
       raw_power <= evaluation.value.epsilon.complete_max; ++raw_power)
    coefficients.push_back(evaluation.value.at(
        static_cast<std::int32_t>(raw_power), 0));
  EndpointLimitResult result;
  result.values.emplace_back(evaluation.value.epsilon,
                             std::move(coefficients));
  result.imaginary_sign = derived_rim.value_or(1);
  return result;
}

json::object transport_endpoint_provenance(
    const std::string& checkpoint_identity, const json::object& source,
    const json::object& analytic_metadata) {
  return json::object{
      {"schema",
       "diffexp3-retained-native-transport-endpoint-result-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"source", source},
      {"cancellation", json::object{
           {"mode", "exact-or-acb-singleton"}}},
      {"analytic_metadata", analytic_metadata}};
}

std::vector<LocalSolution<ComplexBall>>
project_terminal_acb_basis_row(
    const StoredPlannedMatchHop& match,
    const json::object& row,
    std::optional<std::int32_t> projected_complete_cap,
    const std::string& context) {
  // Preserve the same stable association used by the certified match:
  //
  //                         G = (F T) P.
  //
  // Projecting/integrating the confluent physical columns F independently
  // and multiplying the resulting row by T P afterwards is algebraically
  // equivalent, but it loses the exact coefficientwise cancellations which
  // make G epsilon-regular.  Build G first, then apply the observable to its
  // columns and finally contract with the already retained transformed
  // weights.
  auto basis = match.terminal_acb_factorized_basis(
      context + ":factorized-basis");
  std::vector<LocalSolution<ComplexBall>> projected;
  projected.reserve(basis->size());
  for (std::size_t column = 0; column < basis->size(); ++column) {
    // A shorter Taylor prefix may have been selected to make the pointwise
    // matching solve numerically stable.  That prefix records how the
    // constant connection weights were obtained; it is not a new local
    // solution and cannot truncate a downstream line or endpoint functional.
    // Apply the certified weights to the complete retained physical basis.
    auto source = (*basis)[column];
    auto matrix = parse_prepared_rational_row<ComplexBall>(
        row, source, projected_complete_cap);
    auto output = apply_prepared_scalar_row_window(
        matrix, source,
        projected_complete_cap.value_or(
            std::numeric_limits<std::int32_t>::max()),
        context + ":physical-column:" +
            std::to_string(column));
    auto scalar = output.has_value()
        ? std::move(*output)
        : exact_zero_scalar_local_like(
              source, context + ":physical-zero-column:" +
                  std::to_string(column));
    if (scalar.dimension != 1)
      throw std::logic_error(
          context +
          ": terminal factorized row projection did not remain scalar");
    projected.push_back(std::move(scalar));
  }
  return projected;
}

std::vector<LocalSolution<ComplexBall>>
project_terminal_acb_physical_basis_row(
    const StoredPlannedMatchHop& match,
    const json::object& row,
    std::int32_t projected_complete_cap,
    const std::string& context) {
  auto owners = match.terminal_acb_basis_owners();
  std::vector<LocalSolution<ComplexBall>> projected;
  projected.reserve(owners.size());
  for (std::size_t column = 0; column < owners.size(); ++column) {
    auto source = owners[column]->solution();
    const auto complete_max = std::min(
        owners[column]->top_valid(),
        source.epsilon.complete_max);
    if (complete_max < source.epsilon.min_power)
      throw std::domain_error(
          context +
          ": terminal physical basis has no valid epsilon coefficient");
    if (complete_max < source.epsilon.complete_max)
      source = restrict_local_epsilon_frame_strict_lower(
          source, source.epsilon.min_power, complete_max,
          context + ":physical-valid:column:" +
              std::to_string(column));
    auto matrix = parse_prepared_rational_row<ComplexBall>(
        row, source, projected_complete_cap);
    auto output = apply_prepared_scalar_row_window(
        matrix, source, projected_complete_cap,
        context + ":physical-column:" +
            std::to_string(column));
    auto scalar = output.has_value()
        ? std::move(*output)
        : exact_zero_scalar_local_like(
              source, context + ":physical-zero-column:" +
                  std::to_string(column));
    if (scalar.dimension != 1)
      throw std::logic_error(
          context +
          ": terminal physical row projection did not remain scalar");
    projected.push_back(std::move(scalar));
  }
  return projected;
}

EpsilonFrame<ComplexBall> scalar_epsilon_frame(
    const EpsilonVector& value, const std::string& context) {
  if (value.dimension != 1)
    throw std::logic_error(
        context + ": expected one scalar epsilon vector");
  std::vector<ComplexBall> coefficients;
  coefficients.reserve(value.epsilon.width());
  for (std::int64_t raw_power = value.epsilon.min_power;
       raw_power <= value.epsilon.complete_max; ++raw_power)
    coefficients.push_back(
        value.at(static_cast<std::int32_t>(raw_power), 0));
  return EpsilonFrame<ComplexBall>(
      value.epsilon, std::move(coefficients));
}

std::int32_t terminal_factorized_physical_complete_cap(
    const StoredPlannedMatchHop& match,
    std::int32_t requested_output_complete_max,
    std::int32_t functional_epsilon_loss,
    const std::string& context,
    bool contract_physical_weights = false) {
  if (functional_epsilon_loss < 0)
    throw std::invalid_argument(
        context +
        ": terminal functional epsilon loss cannot be negative");
  const auto& weights = match.terminal_acb_physical_weights();
  if (weights.empty())
    throw std::logic_error(
        context + ": terminal match has no physical weights");
  // If the final contraction is sum_j L(G_j) y_j, coefficient q needs
  // L(G_j) through q-min(y_j) for every column.  The leading order of the
  // product G_j*y_j is not a substitute: a +N basis column multiplied by a
  // -N weight can have product valuation zero while still requiring N extra
  // functional coefficients.  Using the product valuation here silently
  // truncated precisely the columns responsible for terminal cancellation.
  const auto minimum_weight_power = contract_physical_weights
      ? match.terminal_acb_physical_weight_min_power()
      : match.terminal_acb_transformed_weight_min_power();
  const auto label =
      context + ": terminal physical input complete maximum";
  return local_algebra_detail::checked_i32(
      static_cast<std::int64_t>(
          requested_output_complete_max) -
          minimum_weight_power +
          functional_epsilon_loss,
      label.c_str());
}

EndpointLimitResult terminal_factorized_endpoint_limit(
    const StoredPlannedMatchHop& match,
    const json::object& row,
    const ObservableEpsilonContract& epsilon_contract,
    const Magnitude& publication_relative_tolerance,
    const ResolvedTransportEndpointBinding& binding,
    const std::string& context) {
  matching_detail::ScopedAcbPrecision exact_shadow_precision(
      match.terminal_acb_extra_precision_bits());
  const auto* route_raw =
      std::getenv("DIFFEXP_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE");
  const std::string route =
      route_raw == nullptr ? std::string() : std::string(route_raw);
  if (!route.empty() && route != "physical" && route != "factorized" &&
      route != "adjoint" && route != "compare")
    throw std::invalid_argument(
        "DIFFEXP_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE must be "
        "factorized, physical, adjoint, or compare");
  const auto* composed_mode_raw =
      std::getenv("DIFFEXP_DIAGNOSTIC_TERMINAL_COMPOSED_ADJOINT");
  const bool authoritative_center_physical =
      route.empty() && binding.centered &&
      composed_mode_raw != nullptr &&
      std::string(composed_mode_raw) == "authoritative";
  const bool direct_physical =
      route == "physical" || authoritative_center_physical;
  if (authoritative_center_physical)
    std::cerr
        << "terminal-endpoint-selected status="
           "classified-direct-physical-center-fallback"
        << " detail=exact-center-endpoint-has-no-ordinary-composed-"
           "adjoint-pairing\n";
  // The adjoint factorization can be a genuinely Laurent system.  Its
  // inverse valuation and elimination-tail loss are properties of the live
  // endpoint matrix, not a universal integer guard derivable from the
  // requested output window.  This is a final consumer, so retain every
  // honestly available prepared-row coefficient and let the adjoint's exact
  // finite-window checks decide whether the upstream reservoir is sufficient.
  const std::optional<std::int32_t> projection_cap = direct_physical
      ? std::optional<std::int32_t>(terminal_factorized_physical_complete_cap(
            match, epsilon_contract.requested.complete_max, 0,
            context + ": direct physical endpoint window", true))
      : std::nullopt;
  auto projected = direct_physical
      ? project_terminal_acb_physical_basis_row(
            match, row, *projection_cap,
            context + ":physical-row")
      : project_terminal_acb_basis_row(
            match, row, projection_cap,
            context + ":row");
  FiniteLaurentMatrix<ComplexBall> physical_rows;
  EndpointLimitResult result;
  result.imaginary_sign = binding.rim.value_or(1);

  if (!binding.centered) {
    FiniteLaurentVector<ComplexBall> point_row;
    point_row.reserve(projected.size());
    const auto point =
        RealEvaluationPoint::rational(binding.local_end.str());
    for (std::size_t column = 0; column < projected.size(); ++column) {
      EvaluationOptions options;
      options.imaginary_sign = binding.rim;
      options.compute_tail_estimate = false;
      const auto evaluation =
          evaluate_local_solution(projected[column], point, options);
      if (!evaluation.value.error.empty())
        throw std::domain_error(
            context +
            ": terminal factorized point evaluation produced an unsupported error envelope");
      point_row.push_back(scalar_epsilon_frame(
          evaluation.value,
          context + ":point-column:" +
              std::to_string(column)));
    }
    physical_rows.push_back(std::move(point_row));
    auto contracted = direct_physical
        ? match.contract_terminal_acb_physical_functionals(
              physical_rows,
              context + ":direct-physical-point-contraction")
        : match.adjoint_contract_terminal_acb_functionals(
              physical_rows, epsilon_contract.required_complete_max,
              context + ":factorized-adjoint-point-contraction", true,
              publication_relative_tolerance);
    if (contracted.size() != 1)
      throw std::logic_error(
          context +
          ": terminal point contraction changed its scalar row count");
    result.values.push_back(std::move(contracted.front()));
    return result;
  }

  using DivergentKey = std::pair<std::string, std::uint32_t>;
  std::vector<EpsilonFrame<ComplexBall>> finite_columns;
  std::vector<std::map<DivergentKey, EpsilonFrame<ComplexBall>>>
      divergent_columns(projected.size());
  std::set<DivergentKey> divergent_keys;
  finite_columns.reserve(projected.size());
  for (std::size_t column = 0; column < projected.size(); ++column) {
    const auto& local = projected[column];
    std::vector<ComplexBall> finite(
        local.epsilon.width(), ComplexBall(0));
    std::map<DivergentKey, std::vector<ComplexBall>>
        divergent_coefficients;
    for (const auto& sector : local.sectors) {
      if (sector.b.is_zero == TruthValue::No) {
        if (integration_detail::material_sector(sector))
          ++result.dropped_regulated_sectors;
        continue;
      }
      if (sector.b.domain == ExactDomain::SymbolicRational)
        throw integration_detail::unsupported_symbolic_regulator();
      if (sector.b.is_zero != TruthValue::Yes)
        throw NativeIntegrationError(
            NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
            "terminal endpoint classification requires an exact regulator zero fact");
      if (sector.a.domain != ExactDomain::Rational)
        throw NativeIntegrationError(
            NativeIntegrationErrorCode::UnsupportedExactTag, "E10",
            "terminal endpoint classification requires rational local powers");
      const Rational a(sector.a.canonical);
      const auto first_unseen = ExactScalarDescriptor::rational(
          (a + Rational(
                   static_cast<long>(local.taylor_width())))
              .str());
      if (first_unseen.sign != ExactSign::Positive)
        throw NativeIntegrationError(
            NativeIntegrationErrorCode::IncompleteTaylorWindow, "E10",
            "terminal endpoint limit is uncertified: the first unseen Taylor cell can still have nonpositive absolute power");
      for (std::uint32_t taylor = 0;
           taylor < local.taylor_width(); ++taylor) {
        const auto tag = sector_monomial_tag(sector, taylor);
        const auto classification =
            classify_endpoint_limit_cell(tag);
        if (classification != EndpointCellClass::Finite &&
            classification != EndpointCellClass::Divergent)
          continue;
        std::vector<ComplexBall>* target = nullptr;
        DivergentKey key;
        if (classification == EndpointCellClass::Finite) {
          target = &finite;
        } else {
          key = {tag.m.canonical, tag.log_power};
          divergent_keys.insert(key);
          auto [found, inserted] =
              divergent_coefficients.try_emplace(
                  key, local.epsilon.width(), ComplexBall(0));
          (void)inserted;
          target = &found->second;
        }
        for (std::size_t epsilon = 0;
             epsilon < local.epsilon.width(); ++epsilon)
          (*target)[epsilon] +=
              sector.coefficients[local_detail::sector_index(
                  local, epsilon, taylor, 0)];
      }
    }
    finite_columns.emplace_back(
        local.epsilon, std::move(finite));
    for (auto& [key, coefficients] : divergent_coefficients)
      divergent_columns[column].emplace(
          key, EpsilonFrame<ComplexBall>(
                   local.epsilon, std::move(coefficients)));
  }

  physical_rows.push_back(std::move(finite_columns));
  for (const auto& key : divergent_keys) {
    FiniteLaurentVector<ComplexBall> row_frames;
    row_frames.reserve(projected.size());
    for (std::size_t column = 0; column < projected.size(); ++column) {
      const auto found = divergent_columns[column].find(key);
      if (found != divergent_columns[column].end()) {
        row_frames.push_back(found->second);
      } else {
        row_frames.emplace_back(
            projected[column].epsilon,
            std::vector<ComplexBall>(
                projected[column].epsilon.width(), ComplexBall(0)));
      }
    }
    physical_rows.push_back(std::move(row_frames));
  }
  auto contracted = direct_physical
      ? match.contract_terminal_acb_physical_functionals(
            physical_rows,
            context + ":direct-physical-center-contraction")
      : match.adjoint_contract_terminal_acb_functionals(
            physical_rows, epsilon_contract.required_complete_max,
            context + ":factorized-adjoint-center-contraction", true,
            publication_relative_tolerance);
  if (contracted.size() != 1 + divergent_keys.size())
    throw std::logic_error(
        context +
        ": terminal endpoint contraction changed its tagged row count");
  result.values.push_back(std::move(contracted.front()));
  std::size_t row_index = 1;
  for (const auto& key : divergent_keys) {
    const auto& frame = contracted[row_index++];
    for (std::int64_t raw_power = frame.min_power();
         raw_power <= frame.complete_max(); ++raw_power) {
      const auto power = static_cast<std::int32_t>(raw_power);
      const auto coefficient = frame.coefficient(power);
      if (coefficient.is_zero()) {
        ++result.cancelled_divergent_coefficients;
        continue;
      }
      const bool uncertified = coefficient.contains_zero();
      NativeIntegrationError error(
          uncertified
              ? NativeIntegrationErrorCode::UncertifiedCancellation
              : NativeIntegrationErrorCode::DivergentEndpoint,
          uncertified ? "E10" : "E2",
          uncertified
              ? "factorized terminal endpoint divergence contains zero but is not the exact singleton zero"
              : "factorized terminal endpoint divergent coefficient is nonzero");
      error.absolute_power = key.first;
      error.log_power = key.second;
      error.epsilon_power = power;
      error.component = 0;
      throw error;
    }
  }
  return result;
}

json::object compact_prepared_row_provenance_record(
    std::size_t tile, const json::object& prepared_row);

std::shared_ptr<StoredEndpointResult>
build_transport_endpoint_row_from_terminal_local(
    const std::string& endpoint_handle,
    const std::string& checkpoint_identity,
    const std::string& observable_identity,
    const json::object& row,
    const ObservableEpsilonContract& epsilon_contract,
    const json::object& epsilon_record,
    const Magnitude& publication_relative_tolerance,
    const std::string& publication_relative_tolerance_text,
    const std::string& projected_handle,
    const std::string& projected_checkpoint_identity,
    const std::string& domain, slong precision_bits,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& source,
    const std::shared_ptr<StoredPlannedMatchHop>& terminal_match,
    const std::shared_ptr<StoredTransportArmState>& transport_owner,
    const ResolvedTransportEndpointBinding& binding) {
  if (endpoint_handle.empty() || checkpoint_identity.empty() ||
      observable_identity.empty() || projected_handle.empty() ||
      projected_checkpoint_identity.empty() || !plan || !source)
    throw std::invalid_argument(
        "transport endpoint row lost an identity, plan, or terminal local owner");
  const auto source_summary = source->summary();
  const auto started = std::chrono::steady_clock::now();
  validate_prepared_rational_row_structure(
      row, as_u32(source_summary.at("dimension"),
                  "transport endpoint source dimension"),
      "transport endpoint prepared rational row");
  json::object row_request{
      {"row", row},
      {"source_checkpoint_identity", source->checkpoint_identity()},
      {"checkpoint_identity", projected_checkpoint_identity}};
  std::shared_ptr<StoredLocalBase> projected;
  const bool terminal_factorized =
      domain == "acb" &&
      terminal_match != nullptr &&
      terminal_match->has_terminal_acb_factorization();
  if (!terminal_factorized) {
    if (domain == "rational") {
      const auto typed =
          std::dynamic_pointer_cast<StoredLocal<Rational>>(source);
      if (!typed)
        throw std::logic_error(
            "transport endpoint Rational source changed coefficient domain");
      projected = build_rational_row_local<Rational>(
          projected_handle, row_request, precision_bits, typed, source,
          false, epsilon_contract.requested.complete_max);
    } else if (domain == "acb") {
      const auto typed =
          std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(source);
      if (!typed)
        throw std::logic_error(
            "transport endpoint Acb source changed coefficient domain");
      projected = build_rational_row_local<ComplexBall>(
          projected_handle, row_request, precision_bits, typed, source,
          false, epsilon_contract.requested.complete_max);
    } else {
      throw std::invalid_argument(
          "transport endpoint row requires one numeric coefficient domain");
    }
  }

  EndpointLimitResult result;
  if (terminal_factorized) {
    result = terminal_factorized_endpoint_limit(
        *terminal_match, row,
        epsilon_contract, publication_relative_tolerance, binding,
        checkpoint_identity + ":terminal-factorized");
  } else if (binding.centered) {
    EndpointLimitOptions options;
    options.approach_direction = binding.approach_direction;
    options.imaginary_sign = binding.rim;
    options.allow_certified_numeric_cancellation = true;
    result = projected->endpoint_limit(options);
  } else {
    const auto point = RealEvaluationPoint::rational(
        binding.local_end.str());
    const auto evaluation = projected->evaluate_retained_point(
        point, binding.rim);
    result = endpoint_result_from_retained_evaluation(
        evaluation, binding.rim);
  }
  result = restrict_endpoint_result_epsilon(
      std::move(result), epsilon_contract,
      "transport endpoint row result");
  if (result.values.size() != 1)
    throw std::logic_error(
        "transport endpoint row result did not remain scalar");
  const auto analytic_metadata = terminal_factorized
      ? source->exact_analytic_metadata()
      : projected->exact_analytic_metadata();
  auto endpoint_source = bounded_transport_endpoint_source(
      binding, plan, source, transport_owner, transport_owner == nullptr);
  endpoint_source["observable"] = json::object{
      {"identity", observable_identity},
      {"checkpoint_identity", checkpoint_identity}};
  endpoint_source["row"] =
      compact_prepared_row_provenance_record(0, row);
  endpoint_source["output_epsilon_contract"] = epsilon_record;
  endpoint_source["publication_relative_tolerance"] =
      publication_relative_tolerance_text;
  const auto provenance = transport_endpoint_provenance(
      checkpoint_identity, endpoint_source, analytic_metadata);
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredEndpointResult>(
      endpoint_handle, checkpoint_identity, provenance_identity,
      source->handle(), source->source_chart(),
      source->source_operator_identity(), source->checkpoint_identity(),
      source->scalar_domain(), binding.approach_direction, std::nullopt,
      "exact-or-acb-singleton", std::move(analytic_metadata),
      std::move(result), elapsed_ms, std::move(endpoint_source),
      binding.rim,
      transport_owner ? nullptr : plan,
      transport_owner ? nullptr : source,
      transport_owner, row);
}

std::shared_ptr<StoredEndpointResult> build_transport_endpoint_row(
    const std::string& endpoint_handle,
    const std::string& checkpoint_identity,
    const std::string& observable_identity,
    const json::object& row,
    const ObservableEpsilonContract& epsilon_contract,
    const json::object& epsilon_record,
    const Magnitude& publication_relative_tolerance,
    const std::string& publication_relative_tolerance_text,
    const std::string& projected_handle,
    const std::string& projected_checkpoint_identity,
    const std::string& domain, slong precision_bits,
    const std::shared_ptr<StoredTransportArmState>& state,
    const ResolvedTransportEndpointBinding& binding) {
  if (!state)
    throw std::invalid_argument(
        "transport endpoint row requires a retained state");
  return build_transport_endpoint_row_from_terminal_local(
      endpoint_handle, checkpoint_identity, observable_identity, row,
      epsilon_contract, epsilon_record, publication_relative_tolerance,
      publication_relative_tolerance_text, projected_handle,
      projected_checkpoint_identity, domain, precision_bits,
      state->plan_owner(), state->final_local(),
      state->terminal_factorized_match(), state, binding);
}

WholeArmEpsilonContract parse_whole_arm_epsilon_contract(
    const json::value& raw, const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(
      object,
      {"min", "max", "required_complete_max",
       "match_required_complete_max"},
      label);
  WholeArmEpsilonContract result{
      {as_i32(object.at("min"), "whole-arm epsilon minimum"),
       as_i32(object.at("max"), "whole-arm epsilon maximum")},
      as_i32(object.at("required_complete_max"),
             "whole-arm projected/public required epsilon maximum"),
      as_i32(object.at("match_required_complete_max"),
             "whole-arm source/match required epsilon maximum")};
  (void)result.work.width();
  if (result.public_required_complete_max < result.work.min_power ||
      result.match_required_complete_max <
          result.public_required_complete_max ||
      result.match_required_complete_max > result.work.complete_max)
    throw std::invalid_argument(
        "whole-arm public and match required epsilon maxima must lie in the work window with match_required_complete_max >= required_complete_max");
  return result;
}

RetainedLocalFrameContract retained_local_frame_contract(
    const std::shared_ptr<StoredLocalBase>& local) {
  if (!local)
    throw std::invalid_argument(
        "native arm frame intersection received a null local");
  const auto summary = local->summary();
  RetainedLocalFrameContract result{
      {as_i32(summary.at("epsilon_min"), "retained epsilon minimum"),
       as_i32(summary.at("epsilon_max"), "retained epsilon maximum")},
      parse_validity(summary.at("top_valid")),
      as_u32(summary.at("dimension"), "retained local dimension"),
      as_u32(summary.at("taylor_complete_max"),
             "retained local Taylor maximum")};
  (void)result.epsilon.width();
  return result;
}

std::int32_t retained_local_complete_max(
    const std::shared_ptr<StoredLocalBase>& local) {
  const auto frame = retained_local_frame_contract(local);
  return frame.top_valid == kCompleteInfinity
      ? frame.epsilon.complete_max
      : std::min(frame.epsilon.complete_max, frame.top_valid);
}

void require_retained_local_complete_max(
    const std::shared_ptr<StoredLocalBase>& local,
    std::int32_t required_complete_max, const char* label) {
  const auto complete_max = retained_local_complete_max(local);
  if (complete_max < required_complete_max)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InsufficientCompleteWindow,
        std::string(label) +
            " materialized physical source does not cover the globally "
            "required complete epsilon maximum",
        std::nullopt, std::nullopt, complete_max);
}

EpsilonWindow live_match_epsilon_intersection(
    EpsilonWindow requested, std::int32_t required_complete_max,
    const std::shared_ptr<StoredLocalBase>& incoming,
    const std::vector<std::shared_ptr<StoredLocalBase>>& basis) {
  if (basis.empty())
    throw std::invalid_argument(
        "native whole-arm match basis cannot be empty");
  const auto incoming_frame = retained_local_frame_contract(incoming);
  auto union_minimum = incoming_frame.epsilon.min_power;
  auto complete_max = requested.complete_max;
  const auto dimension = incoming_frame.dimension;
  if (basis.size() != dimension)
    throw std::invalid_argument(
        "native whole-arm match requires one receiving column per component");
  const auto admit = [&](const RetainedLocalFrameContract& frame) {
    if (frame.dimension != dimension)
      throw std::invalid_argument(
          "native whole-arm matching local dimensions differ");
    // A finite Laurent lower edge is structural: coefficients below it are
    // exact zero.  Clip the caller's lower edge only to the union of actual
    // local frames, then let each matcher zero-pad locals which begin later.
    // Complete upper edges still intersect.
    union_minimum = std::min(union_minimum, frame.epsilon.min_power);
    complete_max = std::min(complete_max, frame.epsilon.complete_max);
    if (frame.top_valid != kCompleteInfinity)
      complete_max = std::min(complete_max, frame.top_valid);
  };
  admit(incoming_frame);
  for (const auto& column : basis)
    admit(retained_local_frame_contract(column));
  const auto minimum = std::max(requested.min_power, union_minimum);
  if (minimum > complete_max)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InsufficientCompleteWindow,
        "native whole-arm match has no common complete epsilon window",
        std::nullopt, std::nullopt, complete_max);
  if (required_complete_max < minimum ||
      required_complete_max > complete_max)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InsufficientCompleteWindow,
        "native whole-arm live match intersection does not cover the globally required complete epsilon maximum",
        std::nullopt, std::nullopt, complete_max);
  return {minimum, complete_max};
}

EpsilonWindow live_line_epsilon_intersection(
    EpsilonWindow requested, std::int32_t required_complete_max,
    const std::shared_ptr<StoredLocalBase>& local,
    std::int32_t primitive_halo) {
  if (primitive_halo < 0)
    throw std::invalid_argument(
        "native whole-arm primitive halo must be nonnegative");
  const auto frame = retained_local_frame_contract(local);
  // A regulated centre-endpoint primitive can create Laurent coefficients
  // below the integrand frame (at most `primitive_halo` rows).  Preserve
  // those rows here.  Do not pessimistically subtract the same bound from
  // the upper edge: the retained line integrator inspects every exact
  // monomial tag and rejects an insufficient coefficient halo precisely
  // when its primitive really begins at eps^-1.  Applying the worst-case
  // loss here as well discards a valid order from ordinary primitives and
  // double-charges the caller's integrand halo.
  auto minimum = std::max(
      requested.min_power,
      local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(frame.epsilon.min_power) -
              primitive_halo,
          "native whole-arm line deliverable minimum"));
  auto complete_max =
      std::min(requested.complete_max, frame.epsilon.complete_max);
  if (frame.top_valid != kCompleteInfinity)
    complete_max = std::min(complete_max, frame.top_valid);
  if (minimum > complete_max || required_complete_max > complete_max)
    throw std::domain_error(
        "native whole-arm integrand row does not cover the globally required complete epsilon maximum");
  return {minimum, complete_max};
}

json::object native_unit_saturation_request(
    const std::shared_ptr<StoredTilePlan>& plan, const std::string& arm,
    std::size_t match_index, bool compact_plan_reference = false) {
  if (!plan)
    throw std::invalid_argument(
        "native unit-saturation request requires its retained tile plan");
  if (compact_plan_reference)
    return json::object{
        {"schema", kNativeUnitSaturationCompactRequestSchema},
        {"tile_plan", plan->handle()},
        {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
        {"arm", arm}, {"match", match_index}};
  return json::object{
        {"schema", kNativeUnitSaturationRequestSchema},
        {"tile_plan", plan->handle()},
        {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
        {"tile_plan_provenance_identity", plan->provenance_identity()},
        {"arm", arm}, {"match", match_index}};
}

NativeAcbSaturationBinding native_acb_saturation_binding(
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::string& session_configuration_identity,
    const std::string& arm,
    std::size_t match_index, const std::string& match_checkpoint_identity,
    bool compact_plan_reference) {
  if (!plan || session_configuration_identity.empty() ||
      match_checkpoint_identity.empty())
    throw std::invalid_argument(
        "native Acb saturation selection requires its session, retained plan, and match checkpoint");
  const auto& retained = plan->arm(arm);
  if (match_index >= retained.exact.matches.size())
    throw std::invalid_argument(
        "native Acb saturation selection match index is outside its retained arm");
  const auto& exact_match = retained.exact.matches[match_index];
  const auto& certified_match = retained.certified_matches.at(match_index);
  const auto& receiving = retained.charts.at(exact_match.receiving_chart);
  const auto* scc_owner = std::get_if<
      std::shared_ptr<CompositeSCCChartBase>>(&receiving.owner);
  if (scc_owner == nullptr || *scc_owner == nullptr ||
      !is_supported_acb_singular_scc_column_capability(
          (*scc_owner)->column_execution_capability()))
    return {"native_unit_saturation",
            native_unit_saturation_request(
                plan, arm, match_index, compact_plan_reference)};
  const auto receiving_rim = exact_plan_rim(
      receiving.prescriptions, receiving.geometry.scale);
  json::object request{
          {"schema", compact_plan_reference
              ? kNativeSingularSCCSaturationCompactRequestSchema
              : kNativeSingularSCCSaturationRequestSchema},
          {"session_configuration_identity",
           session_configuration_identity},
          {"tile_plan", plan->handle()},
          {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
          {"arm", arm},
          {"match", match_index},
          {"match_checkpoint_identity", match_checkpoint_identity},
          {"receiving_scc", (*scc_owner)->handle()},
          {"receiving_scc_exact_identity", (*scc_owner)->exact_identity()},
          {"receiving_execution_capability",
           (*scc_owner)->column_execution_capability()},
          {"receiving_basis_point_exact",
           certified_match.receiving.local.exact_coordinate},
          {"receiving_basis_point_sign",
           certified_match.receiving.local.sign},
          {"physical_match_point_exact",
           certified_match.receiving.physical_exact},
          {"receiving_rim", optional_plan_rim_json(receiving_rim)}};
  if (!compact_plan_reference)
    request["tile_plan_provenance_identity"] =
        plan->provenance_identity();
  return {"native_singular_scc_saturation", std::move(request)};
}

struct RetainedArmMarchInput {
  std::string name;
  std::vector<std::vector<std::string>> basis_handles;
  std::vector<std::vector<std::shared_ptr<StoredLocalBase>>> basis;
  std::vector<std::string> match_handles;
  std::vector<std::string> local_handles;
};

struct RetainedArmMarchResult {
  std::vector<std::shared_ptr<StoredPlannedMatchHop>> matches;
  std::vector<std::shared_ptr<StoredLocalBase>> tile_sources;
  double elapsed_ms = 0.0;

  const std::shared_ptr<StoredLocalBase>& final_local() const {
    if (tile_sources.empty())
      throw std::logic_error(
          "retained arm march has no tile source");
    return tile_sources.back();
  }
};

struct ConsumingArmMarchResult {
  json::array basis_references;
  json::array match_references;
  json::array release_diagnostics;
  std::vector<std::shared_ptr<StoredLocalBase>> tile_sources;
  double elapsed_ms = 0.0;
};

json::object compact_transport_local_reference(
    const std::shared_ptr<StoredLocalBase>& local) {
  if (!local)
    throw std::invalid_argument(
        "compact transport provenance received a null local");
  return json::object{
      {"local", local->handle()},
      {"chart", local->source_chart()},
      {"source_operator_reference",
       local->source_operator_reference()},
      {"checkpoint_identity", local->checkpoint_identity()},
      {"coefficient_domain", local->scalar_domain()}};
}

std::string arm_checkpoint_identity(const std::string& root,
                                    const std::string& arm,
                                    const char* kind,
                                    std::size_t one_based_index) {
  return root + ":" + arm + ":" + kind + ":" +
         std::to_string(one_based_index);
}

RetainedArmMarchResult march_retained_arm(
    const std::string& domain, slong precision_bits,
    const std::string& session_configuration_identity,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& anchor,
    const RetainedArmMarchInput& input, EpsilonWindow work_epsilon,
    std::int32_t match_required_complete_max,
    const json::object& refinement, const std::string& checkpoint_root,
    bool compact_plan_reference = false) {
  if ((domain != "rational" && domain != "acb") || !plan || !anchor ||
      session_configuration_identity.empty() || checkpoint_root.empty())
    throw std::invalid_argument(
        "retained arm march requires a numeric domain, session identity, plan, anchor, and checkpoint root");
  const auto& retained = plan->arm(input.name);
  const auto match_count = retained.exact.matches.size();
  if (input.basis_handles.size() != match_count ||
      input.basis.size() != match_count ||
      input.match_handles.size() != match_count ||
      input.local_handles.size() != match_count)
    throw std::invalid_argument(
        "retained arm march input does not reproduce its plan match count");
  require_exact_keys(refinement, {"relative_tolerance", "max_steps"},
                     "retained arm march refinement policy");
  if (required_string(refinement, "relative_tolerance").empty() ||
      as_u32(refinement.at("max_steps"),
             "retained arm march refinement steps") > 32)
    throw std::invalid_argument(
        "retained arm march refinement policy is invalid");
  (void)work_epsilon.width();
  if (domain == "acb") ComplexBall::set_precision(precision_bits);

  const auto started = std::chrono::steady_clock::now();
  RetainedArmMarchResult output;
  output.matches.reserve(match_count);
  output.tile_sources.reserve(retained.exact.tiles.size());
  output.tile_sources.push_back(anchor);
  std::shared_ptr<StoredLocalBase> current = anchor;
  for (std::size_t match_index = 0; match_index < match_count;
       ++match_index) {
    const bool terminal_basis_match =
        match_index + 1 == match_count;
    const auto match_checkpoint = arm_checkpoint_identity(
        checkpoint_root, input.name, "match", match_index + 1);
    const auto match_epsilon = live_match_epsilon_intersection(
        work_epsilon, match_required_complete_max, current,
        input.basis[match_index]);
    json::object match_request{
        {"arm", input.name}, {"match", match_index},
        {"epsilon", json::object{
             {"min", match_epsilon.min_power},
             {"max", match_epsilon.complete_max},
             {"required_complete_max", match_required_complete_max}}},
        {"checkpoint_identity", match_checkpoint}};
    if (domain == "acb") {
      auto saturation = native_acb_saturation_binding(
          plan, session_configuration_identity, input.name, match_index,
          match_checkpoint, compact_plan_reference);
      match_request[saturation.request_key] =
          std::move(saturation.request);
      match_request["refinement"] = refinement;
    }
    auto match = build_planned_match_hop(
        input.match_handles[match_index], match_request, domain,
        precision_bits, session_configuration_identity, plan,
        input.basis_handles[match_index], input.basis[match_index],
        current->handle(), current, compact_plan_reference,
        terminal_basis_match);
    std::shared_ptr<StoredLocalBase> next;
    try {
      next = match->materialize(
          input.local_handles[match_index],
          arm_checkpoint_identity(checkpoint_root, input.name, "local",
                                  match_index + 1),
          precision_bits, match, terminal_basis_match);
      require_retained_local_complete_max(
          next, match_required_complete_max, "native retained arm match");
    } catch (const std::exception& error) {
      const auto incoming_frame = retained_local_frame_contract(current);
      std::string frames =
          "incoming=" +
          std::to_string(incoming_frame.epsilon.complete_max) +
          "/valid=" + std::to_string(incoming_frame.top_valid) +
          ",basis=[";
      for (std::size_t column = 0;
           column < input.basis[match_index].size(); ++column) {
        if (column != 0) frames += ",";
        const auto frame = retained_local_frame_contract(
            input.basis[match_index][column]);
        frames += std::to_string(frame.epsilon.complete_max) +
                  "/valid=" + std::to_string(frame.top_valid);
      }
      frames += "]";
      throw std::domain_error(
          std::string(error.what()) + "; match_window=[" +
          std::to_string(match_epsilon.min_power) + "," +
          std::to_string(match_epsilon.complete_max) +
          "]; match_required=" +
          std::to_string(match_required_complete_max) + "; " + frames);
    }
    output.matches.push_back(std::move(match));
    output.tile_sources.push_back(next);
    current = std::move(next);
  }
  if (output.tile_sources.size() != retained.exact.tiles.size())
    throw std::logic_error(
        "retained arm march did not produce one source local per tile");
  output.elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return output;
}

ConsumingArmMarchResult march_retained_arm_consuming(
    const std::shared_ptr<SolverSession>& session,
    const std::string& session_configuration_identity,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& anchor,
    RetainedArmMarchInput& input, EpsilonWindow work_epsilon,
    std::int32_t match_required_complete_max,
    const json::object& refinement, const std::string& checkpoint_root,
    std::unordered_map<std::string, std::size_t>& remaining_basis_uses) {
  if (!session || session->domain == "symbolic" || !plan || !anchor)
    throw std::invalid_argument(
        "consuming arm march requires one live numeric session, plan, and anchor");
  const auto& retained = plan->arm(input.name);
  const auto match_count = retained.exact.matches.size();
  if (input.basis_handles.size() != match_count ||
      input.match_handles.size() != match_count ||
      input.local_handles.size() != match_count)
    throw std::invalid_argument(
        "consuming arm march input does not reproduce its exact plan");
  require_exact_keys(refinement, {"relative_tolerance", "max_steps"},
                     "consuming arm refinement policy");
  (void)work_epsilon.width();
  if (session->domain == "acb")
    ComplexBall::set_precision(session->precision_bits);

  const auto started = std::chrono::steady_clock::now();
  ConsumingArmMarchResult output;
  output.basis_references.reserve(match_count);
  output.match_references.reserve(match_count);
  output.tile_sources.reserve(retained.exact.tiles.size());
  output.tile_sources.push_back(anchor);
  std::shared_ptr<StoredLocalBase> current = anchor;

  for (std::size_t match_index = 0; match_index < match_count;
       ++match_index) {
    std::vector<std::shared_ptr<StoredLocalBase>> basis;
    json::array basis_reference;
    const auto& handles = input.basis_handles[match_index];
    std::size_t retained_local_count = 0;
    std::uint64_t retained_coefficient_count = 0;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during consuming arm march");
      basis.reserve(handles.size());
      basis_reference.reserve(handles.size());
      for (const auto& handle : handles) {
        const auto found = session->locals.find(handle);
        if (found == session->locals.end() || !found->second)
          throw std::invalid_argument(
              "unknown or already consumed receiving-basis local: " +
              handle);
        basis.push_back(found->second);
        basis_reference.push_back(
            compact_transport_local_reference(found->second));
      }
    }
    if (basis.empty())
      throw std::invalid_argument(
          "consuming arm receiving basis cannot be empty");

    const auto match_checkpoint = arm_checkpoint_identity(
        checkpoint_root, input.name, "match", match_index + 1);
    const auto match_epsilon = live_match_epsilon_intersection(
        work_epsilon, match_required_complete_max, current, basis);
    json::object match_request{
        {"arm", input.name}, {"match", match_index},
        {"epsilon", json::object{
             {"min", match_epsilon.min_power},
             {"max", match_epsilon.complete_max},
             {"required_complete_max", match_required_complete_max}}},
        {"checkpoint_identity", match_checkpoint}};
    if (session->domain == "acb") {
      auto saturation = native_acb_saturation_binding(
          plan, session_configuration_identity, input.name, match_index,
          match_checkpoint, true);
      match_request[saturation.request_key] = std::move(saturation.request);
      match_request["refinement"] = refinement;
    }
    auto match = build_planned_match_hop(
        input.match_handles[match_index], match_request, session->domain,
        session->precision_bits, session_configuration_identity, plan,
        handles, basis, current->handle(), current, true);
    auto next = match->materialize(
        input.local_handles[match_index],
        arm_checkpoint_identity(checkpoint_root, input.name, "local",
                                match_index + 1),
        session->precision_bits, match);
    require_retained_local_complete_max(
        next, match_required_complete_max,
        "native consuming retained arm match");

    // At this point the full match, materialized local, tail model and
    // physical residual owner have all been validated.  Seal an immutable
    // lineage record before dropping the coefficient-bearing match owner.
    auto sealed_lineage = next->seal_plan_match_lineage();
    output.basis_references.push_back(std::move(basis_reference));
    output.match_references.push_back(json::object{
        {"index", match_index},
        {"checkpoint_identity", match->checkpoint_identity()},
        {"provenance_identity", match->provenance_identity()},
        {"planned_hop", match->handoff()},
        {"sealed_local_lineage", sealed_lineage}});

    // Consumption is committed per completed exact hop.  Shared basis
    // handles are erased only at their last use across both arms.
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      for (std::size_t column = 0; column < handles.size(); ++column) {
        const auto& handle = handles[column];
        const auto use = remaining_basis_uses.find(handle);
        if (use == remaining_basis_uses.end() || use->second == 0)
          throw std::logic_error(
              "consuming basis last-use accounting underflow");
        --use->second;
        if (use->second == 0) {
          const auto found = session->locals.find(handle);
          if (found == session->locals.end() ||
              found->second.get() != basis[column].get())
            throw std::logic_error(
                "consuming basis registry owner changed before last use");
          session->locals.erase(found);
        }
      }
      retained_local_count = session->locals.size();
      for (const auto& [handle, local] : session->locals) {
        (void)handle;
        const auto count = static_cast<std::uint64_t>(
            local->stats().coefficient_count);
        if (count > std::numeric_limits<std::uint64_t>::max() -
                        retained_coefficient_count)
          throw std::overflow_error(
              "consuming transport coefficient diagnostic overflow");
        retained_coefficient_count += count;
      }
    }
    basis.clear();
    match.reset();
    output.release_diagnostics.push_back(json::object{
        {"arm", input.name}, {"match", match_index},
        {"consumed_handles", handles.size()},
        {"session_locals_after_hop", retained_local_count},
        {"session_local_coefficient_count_after_hop",
         retained_coefficient_count}});
    output.tile_sources.push_back(next);
    current = std::move(next);
  }
  if (output.tile_sources.size() != retained.exact.tiles.size())
    throw std::logic_error(
        "consuming arm march did not produce one local per exact tile");
  output.elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return output;
}

json::array line_aggregate_source_records(
    const std::vector<std::shared_ptr<StoredLocalBase>>& owners) {
  json::array records;
  records.reserve(owners.size());
  for (const auto& owner : owners)
    records.push_back(json::object{
        {"local", owner->handle()}, {"chart", owner->source_chart()},
        {"source_operator_identity", owner->source_operator_identity()},
        {"checkpoint_identity", owner->checkpoint_identity()},
        {"coefficient_domain", owner->scalar_domain()},
        {"analytic_metadata", owner->exact_analytic_metadata()},
        {"retained_derivation", owner->retained_derivation().has_value()
             ? json::value(*owner->retained_derivation())
             : json::value(nullptr)}});
  return records;
}

json::array bounded_line_aggregate_source_records(
    const std::vector<std::shared_ptr<StoredLocalBase>>& owners) {
  json::array records;
  records.reserve(owners.size());
  for (const auto& owner : owners) {
    const auto analytic_identity = json::serialize(
        canonical_json_value(owner->exact_analytic_metadata()));
    json::value derivation_reference = nullptr;
    if (owner->retained_derivation().has_value()) {
      const auto derivation_identity = json::serialize(
          canonical_json_value(*owner->retained_derivation()));
      derivation_reference =
          compact_matching_identity_reference(derivation_identity);
    }
    records.push_back(json::object{
        {"local", owner->handle()}, {"chart", owner->source_chart()},
        {"source_operator_reference",
         compact_matching_identity_reference(
             owner->source_operator_identity())},
        {"checkpoint_identity", owner->checkpoint_identity()},
        {"coefficient_domain", owner->scalar_domain()},
        {"analytic_metadata_reference",
         compact_matching_identity_reference(analytic_identity)},
        {"retained_derivation_reference",
         std::move(derivation_reference)}});
  }
  return records;
}

json::object bounded_line_component_record(
    std::size_t index, std::int32_t sign,
    const std::shared_ptr<StoredLineResult>& component) {
  return json::object{
      {"index", index}, {"sign", sign},
      {"checkpoint_identity", component->checkpoint_identity()},
      {"provenance_reference",
       compact_matching_identity_reference(
           component->provenance_identity())}};
}

std::shared_ptr<StoredLineResult> build_retained_line_aggregate(
    const std::string& handle, const std::string& checkpoint_identity,
    const std::string& arm_name, json::object interval,
    json::object aggregate_record,
    const std::shared_ptr<StoredTilePlan>& plan,
    std::vector<std::shared_ptr<StoredLocalBase>> local_owners,
    const std::vector<std::shared_ptr<StoredLineResult>>& components,
    const std::vector<std::int32_t>& signs, double elapsed_ms,
    bool bounded_provenance = false) {
  if (checkpoint_identity.empty())
    throw std::invalid_argument(
        "native line aggregate checkpoint identity cannot be empty");
  local_owners = unique_line_local_owners(local_owners);
  auto result = aggregate_retained_lines(
      components, signs, "native retained line aggregate");
  json::array component_records;
  component_records.reserve(components.size());
  for (std::size_t index = 0; index < components.size(); ++index) {
    if (bounded_provenance) {
      component_records.push_back(bounded_line_component_record(
          index, signs[index], components[index]));
      continue;
    }
    const auto summary = components[index]->summary();
    component_records.push_back(json::object{
        {"index", index}, {"sign", signs[index]},
        {"checkpoint_identity", components[index]->checkpoint_identity()},
        {"provenance_identity", components[index]->provenance_identity()},
        {"scope", summary.at("scope")},
        {"source", summary.at("source")},
        {"interval", summary.at("interval")},
        {"epsilon", json::object{
             {"min", components[index]->result().value.epsilon.min_power},
             {"max", components[index]->result().value.epsilon.complete_max}}},
        {"error", summary.at("error")}});
  }
  aggregate_record["component_count"] = components.size();
  aggregate_record["components"] = std::move(component_records);
  json::object provenance{
      {"schema", bounded_provenance
           ? "diffexp3-retained-native-line-aggregate-v2"
           : "diffexp3-retained-native-line-aggregate-v1"},
      {"checkpoint_identity", checkpoint_identity},
      {"arm", arm_name},
      {"interval", std::move(interval)},
      {"source", json::object{
           {"tile_plan", plan->handle()},
           {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
           {"locals", bounded_provenance
                ? bounded_line_aggregate_source_records(local_owners)
                : line_aggregate_source_records(local_owners)}}},
      {"aggregate", std::move(aggregate_record)},
      {"epsilon", json::object{{"min", result.value.epsilon.min_power},
                                {"max", result.value.epsilon.complete_max}}},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  return std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, std::move(result),
      elapsed_ms, plan, std::move(local_owners), std::move(provenance));
}

std::shared_ptr<StoredLineResult> build_compact_transport_observable_line(
    const std::string& handle, const std::string& checkpoint_identity,
    const std::string& arm_name, json::object interval,
    json::object aggregate_record,
    const std::shared_ptr<StoredTransportArmState>& transport_state,
    const std::vector<std::shared_ptr<StoredLineResult>>& tile_lines,
    double elapsed_ms) {
  if (!transport_state || checkpoint_identity.empty() || tile_lines.empty())
    throw std::invalid_argument(
        "compact transport observable requires its state, checkpoint, and tile values");
  auto result = aggregate_retained_lines(
      tile_lines, std::vector<std::int32_t>(tile_lines.size(), 1),
      "compact native transport observable aggregate");
  aggregate_record["tile_count"] = tile_lines.size();
  json::object provenance{
      {"schema", "diffexp3-retained-native-transport-observable-line-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"arm", arm_name},
      {"interval", std::move(interval)},
      {"source", json::object{
           {"tile_plan", transport_state->plan_owner()->handle()},
           {"tile_plan_checkpoint_identity",
            transport_state->plan_owner()->checkpoint_identity()},
           {"transport_state", transport_state->handle()},
           {"transport_state_checkpoint_identity",
            transport_state->checkpoint_identity()}}},
      {"aggregate", std::move(aggregate_record)},
      {"epsilon", json::object{{"min", result.value.epsilon.min_power},
                                {"max", result.value.epsilon.complete_max}}},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  return std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, std::move(result),
      elapsed_ms, transport_state->plan_owner(),
      std::vector<std::shared_ptr<StoredLocalBase>>{}, std::move(provenance),
      transport_state);
}

std::shared_ptr<StoredLineResult>
build_compact_transport_observable_line_from_result(
    const std::string& handle, const std::string& checkpoint_identity,
    const std::string& arm_name, json::object interval,
    json::object aggregate_record,
    const std::shared_ptr<StoredTransportArmState>& transport_state,
    StoredLineIntegral result, std::size_t tile_count,
    double elapsed_ms) {
  if (!transport_state || checkpoint_identity.empty() || tile_count == 0)
    throw std::invalid_argument(
        "streaming compact transport observable requires its state, checkpoint, and tile count");
  aggregate_record["tile_count"] = tile_count;
  json::object provenance{
      {"schema", "diffexp3-retained-native-transport-observable-line-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"arm", arm_name},
      {"interval", std::move(interval)},
      {"source", json::object{
           {"tile_plan", transport_state->plan_owner()->handle()},
           {"tile_plan_checkpoint_identity",
            transport_state->plan_owner()->checkpoint_identity()},
           {"transport_state", transport_state->handle()},
           {"transport_state_checkpoint_identity",
            transport_state->checkpoint_identity()}}},
      {"aggregate", std::move(aggregate_record)},
      {"epsilon", json::object{{"min", result.value.epsilon.min_power},
                                {"max", result.value.epsilon.complete_max}}},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  return std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, std::move(result),
      elapsed_ms, transport_state->plan_owner(),
      std::vector<std::shared_ptr<StoredLocalBase>>{}, std::move(provenance),
      transport_state);
}

json::object compact_transport_state_reference(
    const std::shared_ptr<StoredTransportArmState>& state) {
  if (!state)
    throw std::invalid_argument(
        "compact transport reference requires a retained state");
  return json::object{
      {"handle", state->handle()},
      {"checkpoint_identity", state->checkpoint_identity()}};
}

json::object compact_transport_pair_source_side(
    const std::shared_ptr<StoredTransportArmState>& state) {
  return json::object{
      {"tile_plan", state->plan_owner()->handle()},
      {"tile_plan_checkpoint_identity",
       state->plan_owner()->checkpoint_identity()},
      {"transport_state", compact_transport_state_reference(state)}};
}

std::shared_ptr<StoredLineResult>
build_compact_transport_pair_observable_line(
    const std::string& handle, const std::string& checkpoint_identity,
    json::object aggregate_record,
    const std::shared_ptr<StoredTransportArmState>& lower_state,
    const std::shared_ptr<StoredTransportArmState>& upper_state,
    const std::shared_ptr<StoredLineResult>& lower_line,
    const std::shared_ptr<StoredLineResult>& upper_line) {
  const auto started = std::chrono::steady_clock::now();
  if (checkpoint_identity.empty() || !lower_line || !upper_line)
    throw std::invalid_argument(
        "compact transport-pair observable requires its checkpoint and both arm values");
  require_transport_pair_compatibility(
      lower_state, upper_state,
      lower_state ? lower_state->anchor_owner()->scalar_domain() : "");
  auto result = aggregate_retained_lines(
      {lower_line, upper_line}, {-1, 1},
      "compact native transport-pair observable aggregate");
  const auto& lower_exact =
      lower_state->plan_owner()->arm("lower").exact;
  const auto& upper_exact =
      upper_state->plan_owner()->arm("upper").exact;
  json::object provenance{
      {"schema",
       "diffexp3-retained-native-transport-pair-observable-line-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"arm", "combined"},
      {"interval", json::object{
           {"from_exact", lower_exact.to.str()},
           {"anchor_exact", lower_exact.from.str()},
           {"to_exact", upper_exact.to.str()}}},
      {"source", json::object{
           {"lower", compact_transport_pair_source_side(lower_state)},
           {"upper", compact_transport_pair_source_side(upper_state)},
           {"common_anchor", json::object{
                {"local", lower_state->anchor_owner()->handle()},
                {"checkpoint_identity",
                 lower_state->anchor_owner()->checkpoint_identity()},
                {"source_operator_identity",
                 lower_state->anchor_owner()->source_operator_identity()},
                {"physical_exact", lower_exact.from.str()}}}}},
      {"aggregate", std::move(aggregate_record)},
      {"epsilon", json::object{
           {"min", result.value.epsilon.min_power},
           {"max", result.value.epsilon.complete_max}}},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  const auto elapsed_ms = std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - started).count();
  return std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, std::move(result),
      elapsed_ms,
      std::vector<std::shared_ptr<StoredTilePlan>>{
          lower_state->plan_owner(), upper_state->plan_owner()},
      std::vector<std::shared_ptr<StoredLocalBase>>{},
      std::move(provenance),
      std::vector<std::shared_ptr<StoredTransportArmState>>{
          lower_state, upper_state});
}

json::object encode_transport_line_value_diagnostics(
    const StoredLineIntegral& line) {
  const auto& result = line.value;
  if (result.dimension == 0)
    throw std::logic_error(
        "transport line-value diagnostics lost their scalar dimension");
  json::array entries;
  entries.reserve(result.epsilon.width() * result.dimension);
  for (std::int64_t raw_power = result.epsilon.min_power;
       raw_power <= result.epsilon.complete_max; ++raw_power) {
    const auto power = static_cast<std::int32_t>(raw_power);
    for (std::uint32_t component = 0; component < result.dimension;
         ++component) {
      const auto& value = result.at(power, component);
      entries.push_back(json::object{
          {"power", power},
          {"component", component},
          {"midpoint", json::array{
               value.real_midpoint(18), value.imag_midpoint(18)}},
          {"radius2exp", json::array{
               value.real_radius_exponent(),
               value.imag_radius_exponent()}},
          {"abs_upper_approx",
           Magnitude::upper_abs(value).approximate_upper()},
          {"contains_zero", value.contains_zero()}});
    }
  }
  return json::object{
      {"epsilon", json::object{
           {"min", result.epsilon.min_power},
           {"max", result.epsilon.complete_max}}},
      {"dimension", result.dimension},
      {"entries", std::move(entries)}};
}

// A rolling stream needs enough per-tile information to localize a bad
// contribution, but retaining one JSON object for every epsilon coefficient
// of every observable defeats the coefficient-streaming memory contract.
// The completed arm/pair still publishes full coefficient-wise conditioning;
// historical tiles retain only this bounded diagnostic summary.
json::object encode_compact_transport_line_value_diagnostics(
    const StoredLineIntegral& line) {
  const auto& result = line.value;
  if (result.dimension == 0)
    throw std::logic_error(
        "compact transport line-value diagnostics lost their scalar dimension");
  double max_abs_upper_approx = 0.0;
  std::uint64_t contains_zero = 0;
  std::uint64_t coefficients = 0;
  for (std::int64_t raw_power = result.epsilon.min_power;
       raw_power <= result.epsilon.complete_max; ++raw_power) {
    const auto power = static_cast<std::int32_t>(raw_power);
    for (std::uint32_t component = 0; component < result.dimension;
         ++component) {
      const auto& value = result.at(power, component);
      max_abs_upper_approx = std::max(
          max_abs_upper_approx,
          Magnitude::upper_abs(value).approximate_upper());
      if (value.contains_zero()) ++contains_zero;
      ++coefficients;
    }
  }
  return json::object{
      {"schema", "diffexp3-compact-rolling-tile-diagnostic-v1"},
      {"epsilon", json::object{
           {"min", result.epsilon.min_power},
           {"max", result.epsilon.complete_max}}},
      {"dimension", result.dimension},
      {"coefficients", coefficients},
      {"contains_zero", contains_zero},
      {"max_abs_upper_approx", max_abs_upper_approx}};
}

json::object encode_transport_pair_conditioning_diagnostics(
    const StoredLineIntegral& lower, const StoredLineIntegral& upper,
    const StoredLineIntegral& combined,
    json::array lower_tiles = {}, json::array upper_tiles = {}) {
  const auto& result = combined.value;
  if (lower.value.dimension != result.dimension ||
      upper.value.dimension != result.dimension ||
      result.dimension == 0)
    throw std::logic_error(
        "transport-pair conditioning diagnostics lost their scalar dimensions");

  json::array entries;
  entries.reserve(result.epsilon.width() * result.dimension);
  for (std::int64_t raw_power = result.epsilon.min_power;
       raw_power <= result.epsilon.complete_max; ++raw_power) {
    const auto power = static_cast<std::int32_t>(raw_power);
    for (std::uint32_t component = 0; component < result.dimension;
         ++component) {
      const auto lower_value =
          power < lower.value.epsilon.min_power
          ? ComplexBall(0)
          : lower.value.at(power, component);
      const auto upper_value =
          power < upper.value.epsilon.min_power
          ? ComplexBall(0)
          : upper.value.at(power, component);
      const auto& combined_value = result.at(power, component);
      const auto lower_upper =
          Magnitude::upper_abs(lower_value).approximate_upper();
      const auto upper_upper =
          Magnitude::upper_abs(upper_value).approximate_upper();
      const auto combined_upper =
          Magnitude::upper_abs(combined_value).approximate_upper();
      const auto combined_lower =
          Magnitude::lower_abs(combined_value).approximate_upper();
      const auto cancellation_scale = lower_upper + upper_upper;
      const auto condition = combined_lower > 0.0
          ? cancellation_scale / combined_lower
          : std::numeric_limits<double>::infinity();
      entries.push_back(json::object{
          {"power", power},
          {"component", component},
          {"lower_midpoint", json::array{
               lower_value.real_midpoint(18),
               lower_value.imag_midpoint(18)}},
          {"upper_midpoint", json::array{
               upper_value.real_midpoint(18),
               upper_value.imag_midpoint(18)}},
          {"combined_midpoint", json::array{
               combined_value.real_midpoint(18),
               combined_value.imag_midpoint(18)}},
          {"lower_radius2exp", json::array{
               lower_value.real_radius_exponent(),
               lower_value.imag_radius_exponent()}},
          {"upper_radius2exp", json::array{
               upper_value.real_radius_exponent(),
               upper_value.imag_radius_exponent()}},
          {"combined_radius2exp", json::array{
               combined_value.real_radius_exponent(),
               combined_value.imag_radius_exponent()}},
          {"lower_abs_upper_approx", lower_upper},
          {"upper_abs_upper_approx", upper_upper},
          {"combined_abs_upper_approx", combined_upper},
          {"combined_contains_zero", combined_value.contains_zero()},
          {"cancellation_condition_upper_approx",
           std::isfinite(condition) ? json::value(condition)
                                    : json::value(nullptr)}});
    }
  }
  return json::object{
      {"schema", "diffexp3-transport-pair-conditioning-diagnostics-v1"},
      {"epsilon", json::object{
           {"min", result.epsilon.min_power},
           {"max", result.epsilon.complete_max}}},
      {"dimension", result.dimension},
      {"entries", std::move(entries)},
      {"lower_tiles", std::move(lower_tiles)},
      {"upper_tiles", std::move(upper_tiles)}};
}

enum class TransportTailPolicy { Stored, Attempt, Require };

TransportTailPolicy parse_transport_tail_policy(const json::value& raw,
                                                 const char* label) {
  if (!raw.is_string())
    throw std::invalid_argument(std::string(label) + " must be a string");
  const std::string value(raw.as_string());
  if (value == "stored") return TransportTailPolicy::Stored;
  if (value == "attempt") return TransportTailPolicy::Attempt;
  if (value == "require") return TransportTailPolicy::Require;
  throw std::invalid_argument(
      std::string(label) + " must be stored, attempt, or require");
}

const char* transport_tail_policy_name(TransportTailPolicy policy) {
  switch (policy) {
    case TransportTailPolicy::Stored:
      return "stored";
    case TransportTailPolicy::Attempt:
      return "attempt";
    case TransportTailPolicy::Require:
      return "require";
  }
  throw std::logic_error("unknown transport tail policy");
}

const char* transport_projection_mode_name(TransportTailPolicy policy) {
  return policy == TransportTailPolicy::Stored
      ? "fused-stored-hash-monomial-stream-v2"
      : "materialized-row-local-v1";
}

void require_transport_projection_mode(
    const json::object& aggregate, TransportTailPolicy policy,
    const char* label) {
  if (required_string(aggregate, "projection_mode") !=
      transport_projection_mode_name(policy))
    throw std::invalid_argument(
        std::string(label) +
        " does not match its stored-tail projection algorithm");
}

json::array compact_prepared_row_entry_facts(
    const json::object& prepared_row) {
  json::array facts;
  const auto& entries = as_array(
      prepared_row.at("entries"), "compact prepared-row entries");
  facts.reserve(entries.size());
  for (const auto& raw_entry : entries) {
    const auto& entry = as_object(
        raw_entry, "compact prepared-row entry");
    const auto& multiplier = as_object(
        entry.at("multiplier"), "compact prepared-row multiplier");
    facts.push_back(json::object{
        {"column", as_u32(entry.at("column"),
                          "compact prepared-row column")},
        {"epsilon_shift",
         as_i32(multiplier.at("epsilon_shift"),
                "compact prepared-row epsilon shift")},
        {"center_pole_order",
         as_u32(multiplier.at("center_pole_order"),
                "compact prepared-row center-pole order")},
        {"exact_identity",
         required_string(multiplier, "exact_identity")}});
  }
  return facts;
}

json::object compact_prepared_row_provenance_record(
    std::size_t tile, const json::object& prepared_row) {
  const auto row_identity = required_string(prepared_row, "exact_identity");
  const auto facts = compact_prepared_row_entry_facts(prepared_row);
  const auto facts_identity = json::serialize(canonical_json_value(facts));
  return json::object{
      {"tile", tile},
      {"exact_identity_reference",
       compact_matching_identity_reference(row_identity)},
      {"entry_facts_reference",
       compact_matching_identity_reference(facts_identity)},
      {"entry_count", facts.size()}};
}

struct TransportObservableContractionInput {
  std::string identity;
  std::string checkpoint_identity;
  std::string checkpoint_root;
  std::vector<json::object> rows;
  // Pair contraction borrows prepared rows directly from the immutable
  // synchronous request tree.  This avoids copying the very large exact
  // Taylor-kernel JSON through parse and per-arm staging.  Other seams keep
  // their existing owned rows; exactly one representation may be active.
  std::vector<const json::object*> borrowed_rows;
  ObservableEpsilonContract epsilon;
  json::object epsilon_record;
  TransportTailPolicy tail_policy = TransportTailPolicy::Stored;
  std::optional<BoundedDivergentCancellation> divergent_cancellation;
  std::optional<Magnitude> publication_relative_tolerance;
  // A direct stored-row integral is already a rigorous enclosure. Paired
  // contraction can first test that cheap enclosure against its final
  // publication contract and enable the much more expensive center-value
  // reassociation only for an observable whose combined line is too wide.
  bool factorize_ordinary_stored_rows = true;
  std::vector<std::string> projected_local_handles;
  std::string aggregate_handle;
  json::object aggregate_record;

  std::size_t row_count() const {
    if (!rows.empty() && !borrowed_rows.empty())
      throw std::logic_error(
          "transport observable cannot mix owned and borrowed rows");
    return borrowed_rows.empty() ? rows.size() : borrowed_rows.size();
  }

  const json::object& row(std::size_t index) const {
    if (!rows.empty() && !borrowed_rows.empty())
      throw std::logic_error(
          "transport observable cannot mix owned and borrowed rows");
    if (borrowed_rows.empty()) return rows.at(index);
    const auto* borrowed = borrowed_rows.at(index);
    if (borrowed == nullptr)
      throw std::logic_error("transport observable borrowed a null row");
    return *borrowed;
  }
};

struct TransportObservableContractionResult {
  std::string identity;
  // Populated only by the legacy whole-arm seam, whose aggregate owns the
  // projected locals explicitly.  Transport-state contraction streams and
  // leaves both vectors empty.
  std::vector<std::shared_ptr<StoredLocalBase>> projected;
  std::vector<std::shared_ptr<StoredLineResult>> tile_lines;
  std::shared_ptr<StoredLineResult> aggregate;
  json::array tile_values;
  std::size_t tile_integrations = 0;
  double tile_integration_ms = 0.0;
  double elapsed_ms = 0.0;
};

bool transport_contraction_timing_enabled() {
  return std::getenv("DIFFEXP_DIAGNOSTIC_CONTRACTION_TIMING") != nullptr;
}

void emit_transport_contraction_timing(std::string message) {
  if (!transport_contraction_timing_enabled()) return;
  static std::mutex mutex;
  std::lock_guard<std::mutex> lock(mutex);
  std::cerr << message << '\n';
}

double elapsed_milliseconds(
    const std::chrono::steady_clock::time_point& started) {
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - started)
      .count();
}

template <typename Scalar>
StoredLineIntegral integrate_transport_stored_row_tile(
    slong precision_bits,
    const std::shared_ptr<StoredLocal<Scalar>>& source,
    const json::object& prepared_row,
    const RetainedArmPlan& arm, const std::string& arm_name,
    std::size_t tile_index,
    const ObservableEpsilonContract& epsilon_contract,
    const std::optional<BoundedDivergentCancellation>&
        divergent_cancellation,
    bool factorize_ordinary_row = true) {
  if (!source || tile_index >= arm.exact.tiles.size())
    throw std::invalid_argument(
        "fused transport row integration lost its source or tile");
  if constexpr (std::is_same_v<Scalar, ComplexBall>)
    ComplexBall::set_precision(precision_bits);
  const auto& tile = arm.exact.tiles[tile_index];
  const auto& certified_tile = arm.certified_tiles.at(tile_index);
  const auto& binding = arm.charts.at(tile.chart);
  source->require_exact_plan_binding(
      binding.local_geometry, binding.prescriptions,
      "fused transport row source");
  const auto projection_cap = local_algebra_detail::checked_i32(
      static_cast<std::int64_t>(
          epsilon_contract.requested.complete_max) + 1,
      "fused transport primitive halo");
  // Matching may retain a private epsilon reservoir above the public
  // observable contract.  Decode only the per-entry multiplier prefix needed
  // through the public output plus the exact one-row primitive halo.  The
  // fused convolution below computes and checks the honest product edge from
  // the full source and each multiplier's own retained width.
  auto matrix = parse_prepared_rational_row<Scalar>(
      prepared_row, source->solution(), projection_cap);

  auto projected_min = source->solution().epsilon.min_power;
  auto projected_complete = source->solution().epsilon.complete_max;
  auto projected_top_valid = source->top_valid();
  bool active = false;
  for (const auto& entry : matrix.entries) {
    if (entry.multiplier.structurally_zero()) continue;
    const auto term_min = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(source->solution().epsilon.min_power) +
            entry.multiplier.epsilon_shift,
        "fused transport projected minimum");
    const auto term_complete = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(
            source->solution().epsilon.complete_max) +
            entry.multiplier.epsilon_shift,
        "fused transport projected complete maximum");
    const auto term_valid = shifted_local_validity(
        source->top_valid(), entry.multiplier.epsilon_shift);
    if (!active) {
      projected_min = term_min;
      projected_complete = term_complete;
      projected_top_valid = term_valid;
      active = true;
    } else {
      projected_min = std::min(projected_min, term_min);
      projected_complete = std::min(projected_complete, term_complete);
      projected_top_valid = std::min(projected_top_valid, term_valid);
    }
  }
  projected_complete = std::min(projected_complete, projection_cap);
  if (projected_top_valid != kCompleteInfinity)
    projected_complete = std::min(projected_complete, projected_top_valid);
  if (projected_complete < projected_min)
    throw std::domain_error(
        "fused transport row has no valid projected epsilon coefficient");

  // Exact integration at a regulated centre endpoint may deepen the
  // Laurent frame by one power of epsilon.  Keep the possible lower row,
  // but let the low-level exact-monomial preflight decide whether the upper
  // coefficient halo is actually consumed.  Blindly subtracting one here
  // loses a valid order for every ordinary primitive.
  const auto line_min = std::max(
      epsilon_contract.requested.min_power,
      local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(projected_min) - 1,
          "fused transport line primitive minimum"));
  const auto line_complete = std::min(
      epsilon_contract.requested.complete_max, projected_complete);
  if (line_min > line_complete ||
      epsilon_contract.required_complete_max > line_complete)
    throw std::domain_error(
        "fused transport row does not cover the globally required epsilon maximum");

  const auto rim = exact_plan_rim(
      binding.prescriptions, binding.geometry.scale);
  const bool reverse_local_orientation =
      tile.local_end < tile.local_begin;
  const auto& primitive_begin = reverse_local_orientation
      ? certified_tile.local_end : certified_tile.local_begin;
  const auto& primitive_end = reverse_local_orientation
      ? certified_tile.local_begin : certified_tile.local_end;
  StoredLineIntegrationOptions options;
  options.delivered_epsilon = {line_min, line_complete};
  options.required_complete_max =
      epsilon_contract.required_complete_max;
  options.imaginary_sign = rim;
  options.certified_chart_scale_sign = binding.scale_sign;
  options.divergent_cancellation = divergent_cancellation;

  // Scale bridging can retain an original exact handoff while the rational
  // topology surrogate moves it inward, leaving a zero integration tile
  // which still performs a basis transfer.
  const bool certified_zero =
      certified_tile_zero_length(arm, tile_index);

  // Even a zero primitive must validate the retained branch prescription;
  // otherwise a malformed chart could hide behind a skipped integration.
  std::optional<std::int32_t> chart_sign;
  try {
    chart_sign = derive_chart_imaginary_sign(
        source->solution(), binding.scale_sign);
  } catch (const std::domain_error& error) {
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
        std::string("invalid prepared chart branch prescription: ") +
            error.what());
  }
  if (chart_sign.has_value() && rim.has_value() &&
      *chart_sign != *rim)
    throw NativeIntegrationError(
        NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
        "explicit fused-line branch sign conflicts with the prepared chart");

  if (certified_zero) {
    return certified_zero_physical_line(
        {line_min, line_complete}, 1,
        rim.has_value() ? rim : chart_sign,
        certified_tile.local_begin.sign == 0, false);
  }

  StoredLineIntegral result;
  const auto direct_started = std::chrono::steady_clock::now();
  double direct_ms = 0.0;
  double replay_ms = 0.0;
  double factorized_ms = 0.0;
  bool factorized_eligible = false;
  std::string factorized_reason =
      factorize_ordinary_row ? "not-attempted" : "disabled";
  try {
    result = integrate_prepared_scalar_row_stored(
        matrix, source->solution(), projected_complete,
        primitive_begin, primitive_end, options);
    direct_ms = elapsed_milliseconds(direct_started);
    if constexpr (std::is_same_v<Scalar, ComplexBall>) {
      // A retained ordinary local is a linear map of its center value.  The
      // historical fused integral expands the same uncertain center ball
      // into every Taylor coefficient and then treats those correlated
      // copies independently.  Reassociate that finite map only when the
      // retained q/C equation independently reconstructs the complete
      // stored prefix.  Unsupported/large cases keep the direct result.
      constexpr std::size_t kMaximumOperatorColumns = 256;
      const auto& equation = source->physical_equation();
      const bool column_count_fits =
          source->solution().epsilon.width() <=
          std::numeric_limits<std::size_t>::max() /
              source->solution().dimension;
      const auto column_count = column_count_fits
          ? source->solution().epsilon.width() *
                source->solution().dimension
          : std::numeric_limits<std::size_t>::max();
      if (factorize_ordinary_row && !equation)
        factorized_reason = "missing-physical-equation";
      else if (factorize_ordinary_row && !column_count_fits)
        factorized_reason = "operator-column-count-overflow";
      else if (factorize_ordinary_row &&
               column_count > kMaximumOperatorColumns)
        factorized_reason = "operator-column-cap";
      else if (factorize_ordinary_row) {
        const auto replay_started = std::chrono::steady_clock::now();
        const auto replay =
            prepare_physical_regular_homogeneous_tail_model(
                *equation, source->solution());
        replay_ms = elapsed_milliseconds(replay_started);
        if (replay.status == TailMajorantStatus::Certified &&
            replay.model.has_value()) {
          const auto factorized_started =
              std::chrono::steady_clock::now();
          auto factorized =
              integrate_ordinary_center_stored_row_factorized(
                  *equation, replay.model->reconstructed, matrix,
                  projected_complete, primitive_begin, primitive_end,
                  options, result, kMaximumOperatorColumns);
          factorized_ms = elapsed_milliseconds(factorized_started);
          factorized_eligible = factorized.eligible;
          factorized_reason = factorized.eligible
              ? "eligible" : factorized.reason;
          if (factorized.eligible) {
            if (factorized.integral.value.dimension !=
                    result.value.dimension ||
                factorized.integral.value.epsilon.min_power !=
                    result.value.epsilon.min_power ||
                factorized.integral.value.epsilon.complete_max !=
                    result.value.epsilon.complete_max)
              throw std::logic_error(
                  "factorized ordinary transport row changed its retained "
                  "epsilon coverage");
            for (std::int64_t raw_power =
                     result.value.epsilon.min_power;
                 raw_power <=
                     result.value.epsilon.complete_max;
                 ++raw_power) {
              const auto power =
                  static_cast<std::int32_t>(raw_power);
              for (std::uint32_t component = 0;
                   component < result.value.dimension;
                   ++component)
                if (!acb_overlaps(
                        result.value.at(power, component).raw(),
                        factorized.integral.value
                            .at(power, component).raw()))
                  throw std::domain_error(
                      "factorized ordinary transport row is disjoint from "
                      "the retained direct integral; arm=" +
                      arm_name + "; tile=" +
                      std::to_string(tile_index) +
                      "; epsilon_power=" +
                      std::to_string(power) +
                      "; component=" +
                      std::to_string(component));
            }
            if (std::getenv(
                    "DIFFEXP_DIAGNOSTIC_FACTORIZED_ORDINARY_ROW") != nullptr ||
                environment_flag_is_one(
                    "DIFFEXP_DIAGNOSTIC_TERMINAL_STATE"))
              std::cerr
                  << "factorized-ordinary-stored-row-authority"
                  << " arm=" << arm_name
                  << " tile=" << tile_index
                  << " source_local=" << source->handle()
                  << " operator_columns="
                  << factorized.operator_columns
                  << " policy=certified-q/C-prefix-replay-and-full-overlap-v1"
                  << '\n';
            result = std::move(factorized.integral);
          }
        } else {
          factorized_reason =
              std::string("physical-prefix-replay-") +
              tail_majorant_status_name(replay.status);
        }
      }
      if (factorize_ordinary_row && !factorized_eligible &&
          (std::getenv(
               "DIFFEXP_DIAGNOSTIC_FACTORIZED_ORDINARY_ROW") != nullptr ||
           environment_flag_is_one("DIFFEXP_DIAGNOSTIC_TERMINAL_STATE")))
        std::cerr
            << "factorized-ordinary-stored-row-ineligible"
            << " arm=" << arm_name
            << " tile=" << tile_index
            << " source_local=" << source->handle()
            << " reason=" << factorized_reason
            << '\n';
    }
  } catch (const NativeIntegrationError& error) {
    std::ostringstream detail;
    detail << error.what() << "; arm=" << arm_name
           << "; tile=" << tile_index
           << "; physical_interval=[" << tile.physical_begin.str()
           << "," << tile.physical_end.str() << "]"
           << "; local_interval=[" << tile.local_begin.str()
           << "," << tile.local_end.str() << "]"
           << "; certified_primitive_interval=["
           << primitive_begin.exact_coordinate << ","
           << primitive_end.exact_coordinate << "]"
           << "; certified_primitive_signs=["
           << primitive_begin.sign << "," << primitive_end.sign << "]"
           << "; chart=" << binding.handle
           << "; chart_center=" << binding.geometry.center.str()
           << "; chart_scale=" << binding.geometry.scale.str()
           << "; rim="
           << (rim.has_value() ? std::to_string(*rim) : "none")
           << "; source_local=" << source->handle()
           << "; fused_row_projection=true";
    NativeIntegrationError contextual(error.code, error.id, detail.str());
    contextual.absolute_power = error.absolute_power;
    contextual.log_power = error.log_power;
    contextual.epsilon_power = error.epsilon_power;
    contextual.component = error.component;
    throw contextual;
  }
  if (transport_contraction_timing_enabled()) {
    std::ostringstream timing;
    timing << "transport-contraction-timing"
           << " phase=stored-row"
           << " arm=" << arm_name
           << " tile=" << tile_index
           << " source_local=" << source->handle()
           << " direct_ms=" << direct_ms
           << " replay_ms=" << replay_ms
           << " factorized_ms=" << factorized_ms
           << " factorized_eligible="
           << (factorized_eligible ? "true" : "false")
           << " factorized_reason=" << factorized_reason;
    emit_transport_contraction_timing(timing.str());
  }

  const auto jacobian = reverse_local_orientation
      ? -binding.scale_numeric : binding.scale_numeric;
  const auto oriented_jacobian_exact = reverse_local_orientation
      ? "-(" + binding.local_geometry.scale_exact + ")"
      : binding.local_geometry.scale_exact;
  for (auto& coefficient : result.value.coefficients)
    coefficient *= jacobian;
  if (!result.value.error.empty()) {
    const auto jacobian_upper = Magnitude::upper_abs(jacobian);
    for (auto& bound : result.value.error.absolute)
      bound = bound * jacobian_upper;
    result.value.error.provenance +=
        "; physical_jacobian_exact=" + oriented_jacobian_exact;
  }
  return result;
}

Magnitude terminal_factorized_divergent_scale(
    const FiniteLaurentVector<ComplexBall>& weights,
    const std::vector<
        const StoredLineIntegral::DeferredDivergentGroup*>&
        selected_groups,
    std::int32_t output_power) {
  if (selected_groups.size() != weights.size())
    throw std::logic_error(
        "terminal divergent scale differs from the selected match dimension");
  Magnitude result = Magnitude::zero();
  for (std::size_t column = 0; column < weights.size(); ++column) {
    const auto* group = selected_groups[column];
    if (group == nullptr) continue;
    const auto& weight = weights[column];
    for (std::int64_t raw_weight_power = weight.min_power();
         raw_weight_power <= weight.complete_max();
         ++raw_weight_power) {
      const auto weight_power =
          static_cast<std::int32_t>(raw_weight_power);
      const auto coefficient_power =
          static_cast<std::int64_t>(output_power) -
          raw_weight_power;
      if (coefficient_power < group->coefficients.min_power() ||
          coefficient_power >
              group->coefficients.complete_max())
        continue;
      const auto coefficient_index = static_cast<std::size_t>(
          coefficient_power -
          group->coefficients.min_power());
      result = result +
          group->contribution_scale_uppers.at(
              coefficient_index) *
              Magnitude::upper_abs(
                  weight.coefficient(weight_power));
    }
  }
  return result;
}

void accumulate_terminal_line_diagnostics(
    StoredLineIntegrationDiagnostics& target,
    const StoredLineIntegrationDiagnostics& source) {
  target.input_monomial_cells = checked_diagnostic_sum(
      target.input_monomial_cells, source.input_monomial_cells,
      "terminal line input-monomial diagnostics");
  target.grouped_monomials = checked_diagnostic_sum(
      target.grouped_monomials, source.grouped_monomials,
      "terminal line grouped-monomial diagnostics");
  target.zero_groups_skipped = checked_diagnostic_sum(
      target.zero_groups_skipped, source.zero_groups_skipped,
      "terminal line zero-group diagnostics");
  target.cancelled_divergent_groups = checked_diagnostic_sum(
      target.cancelled_divergent_groups,
      source.cancelled_divergent_groups,
      "terminal line divergent-group diagnostics");
  target.bounded_cancelled_divergent_coefficients =
      checked_diagnostic_sum(
          target.bounded_cancelled_divergent_coefficients,
          source.bounded_cancelled_divergent_coefficients,
          "terminal line bounded-divergence diagnostics");
  target.primitive_evaluations = checked_diagnostic_sum(
      target.primitive_evaluations, source.primitive_evaluations,
      "terminal line primitive diagnostics");
  target.primitive_component_applications =
      checked_diagnostic_sum(
          target.primitive_component_applications,
          source.primitive_component_applications,
          "terminal line primitive-application diagnostics");
  target.primitive_component_reuses = checked_diagnostic_sum(
      target.primitive_component_reuses,
      source.primitive_component_reuses,
      "terminal line primitive-reuse diagnostics");
  target.has_center_endpoint =
      target.has_center_endpoint || source.has_center_endpoint;
}

struct TerminalComposedAdjointDiagnostic {
  EpsilonFrame<ComplexBall> value;
  std::vector<Magnitude> coefficient_tail_uppers;
  std::vector<Magnitude> coefficient_scale_lowers;
  bool can_extend_taylor_exactly = false;
  std::uint32_t taylor_complete_max = 0;
  std::int32_t adjoint_epsilon_input_complete_max = 0;
  std::int32_t adjoint_epsilon_complete_max = 0;
  std::int32_t coefficientwise_first_input_complete_max = 0;
  std::int32_t coefficientwise_last_input_complete_max = 0;
  std::int32_t coefficientwise_min_input_complete_max = 0;
  std::int32_t coefficientwise_max_input_complete_max = 0;
  std::optional<double> last_coefficient_ratio_upper;
  Magnitude certified_output_tail_upper = Magnitude::zero();
  Magnitude recurrence_contraction_upper = Magnitude::zero();
};

bool terminal_composed_adjoint_meets_relative_accuracy(
    const TerminalComposedAdjointDiagnostic& diagnostic,
    const Magnitude& relative_tolerance) {
  if (diagnostic.coefficient_tail_uppers.size() !=
      diagnostic.coefficient_scale_lowers.size())
    throw std::logic_error(
        "terminal composed adjoint tail/scale diagnostics differ in size");
  for (std::size_t index = 0;
       index < diagnostic.coefficient_tail_uppers.size(); ++index)
    if (diagnostic.coefficient_tail_uppers[index] >
        diagnostic.coefficient_scale_lowers[index] * relative_tolerance)
      return false;
  return true;
}

Magnitude terminal_composed_adjoint_selection_tolerance(
    const Magnitude& retained_match_tolerance) {
  // A terminal value is input to the next ladder level, not the final public
  // answer.  Selecting a tail which consumes essentially the entire match
  // tolerance leaves no budget for the following basis reconstruction and
  // turns a certified value into an inconclusive later handoff.  The value
  // can cross several remaining ladder levels, so a merely local two-digit
  // margin is insufficient (banana4 consumed it by level 2).  Reserve eight
  // decimal digits at every published boundary; the adaptive exact-shadow
  // route normally gains tens of digits with one additional Taylor step.
  return retained_match_tolerance * Magnitude::decimal("1e-8");
}

Magnitude terminal_composed_adjoint_relative_tail_upper(
    const TerminalComposedAdjointDiagnostic& diagnostic) {
  if (diagnostic.coefficient_tail_uppers.size() !=
      diagnostic.coefficient_scale_lowers.size())
    throw std::logic_error(
        "terminal composed adjoint tail/scale diagnostics differ in size");
  auto result = Magnitude::zero();
  for (std::size_t index = 0;
       index < diagnostic.coefficient_tail_uppers.size(); ++index)
    result = Magnitude::maximum(
        result, diagnostic.coefficient_tail_uppers[index] /
                    diagnostic.coefficient_scale_lowers[index]);
  return result;
}

std::string terminal_composed_adjoint_indicial_summary(
    const StoredPlannedMatchHop& match) {
  const auto owners = match.terminal_acb_basis_owners();
  if (owners.empty()) return "unavailable:no-basis-owner";
  const auto equation_owner = owners.front()->retained_equation_owner();
  const auto composite =
      std::dynamic_pointer_cast<CompositeSCCChartBase>(equation_owner);
  if (!composite) return "unavailable:not-composite-scc";
  const auto stats = composite->stats_json();
  const auto found = stats.find("block_charts");
  if (found == stats.end() || !found->value().is_array())
    return "unavailable:no-block-charts";
  json::array blocks;
  for (const auto& raw : found->value().as_array()) {
    const auto& block = as_object(raw, "terminal adjoint SCC block");
    json::object summary;
    if (const auto index = block.find("block"); index != block.end())
      summary["block"] = index->value();
    if (const auto vertices = block.find("vertices");
        vertices != block.end())
      summary["vertices"] = vertices->value();
    if (const auto indicial = block.find("exact_affine_jordan_indicial");
        indicial != block.end())
      summary["exact_affine_jordan_indicial"] = indicial->value();
    else if (const auto root = block.find("affine_indicial_root");
             root != block.end())
      summary["affine_indicial_root"] = root->value();
    else
      summary["exact_affine_jordan_indicial"] = nullptr;
    blocks.push_back(std::move(summary));
  }
  return json::serialize(canonical_json_value(blocks));
}

TerminalComposedAdjointDiagnostic
compute_terminal_composed_adjoint_diagnostic(
    const StoredPlannedMatchHop& match,
    const json::object& prepared_row,
    const RetainedArmPlan& arm, const std::string& arm_name,
    std::size_t tile_index,
    const ObservableEpsilonContract& epsilon_contract,
    std::optional<std::uint32_t> requested_taylor_complete_max =
        std::nullopt) {
  if (tile_index >= arm.exact.tiles.size())
    throw std::invalid_argument(
        "terminal composed adjoint tile index is out of range");
  const auto& tile = arm.exact.tiles[tile_index];
  const auto& certified_tile = arm.certified_tiles.at(tile_index);
  const auto& binding = arm.charts.at(tile.chart);
  if (certified_tile.local_end.sign != 0)
    throw std::domain_error(
        "terminal composed adjoint currently requires a center endpoint; "
        "arm=" + arm_name + ", tile=" + std::to_string(tile_index) +
        ", local_interval=[" + tile.local_begin.str() + "," +
        tile.local_end.str() + "]");
  auto owners = match.terminal_acb_basis_owners();
  if (owners.empty())
    throw std::logic_error(
        "terminal composed adjoint has no physical receiving owner");
  const auto& prototype = owners.front()->solution();
  auto row = parse_prepared_rational_row<ComplexBall>(
      prepared_row, prototype, std::nullopt);
  auto incoming = match.terminal_acb_incoming_physical_value();
  if (incoming.size() != prototype.dimension)
    throw std::logic_error(
        "terminal composed adjoint incoming value changed dimension");
  auto incoming_min = incoming.front().min_power();
  for (const auto& component : incoming)
    incoming_min = std::min(incoming_min, component.min_power());
  const auto required_adjoint_complete =
      local_algebra_detail::checked_i32(
          static_cast<std::int64_t>(
              epsilon_contract.required_complete_max) - incoming_min,
          "terminal composed adjoint epsilon maximum");

  auto row_complete = std::numeric_limits<std::int32_t>::max();
  for (const auto& entry : row.entries) {
    if (entry.multiplier.kernels.empty()) continue;
    row_complete = std::min(
        row_complete,
        local_algebra_detail::checked_i32(
            static_cast<std::int64_t>(entry.multiplier.epsilon_shift) +
                static_cast<std::int64_t>(
                    entry.multiplier.kernels.size()) - 1,
            "terminal composed row epsilon maximum"));
  }
  if (row_complete == std::numeric_limits<std::int32_t>::max() ||
      row_complete < required_adjoint_complete)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InsufficientCompleteWindow,
        "terminal composed adjoint row does not cover the required epsilon "
        "window", std::nullopt, std::nullopt, row_complete);

  // This recurrence is anchored by lambda(0)=0 and publishes
  // lambda(m)^T f(m) for the actual tile m -> 0.  Its forcing is therefore
  // -beta*t*q*r.  Do not reuse the sign with which direct integration sorts
  // its local endpoints: signed evaluation at m already carries that
  // orientation, independently of whether m is positive or negative.
  const auto oriented_jacobian = binding.scale_numeric;
  const auto retained_taylor_complete_max = static_cast<std::uint32_t>(
      prototype.taylor_width() - 1);
  auto taylor_complete_max = retained_taylor_complete_max;
  if (requested_taylor_complete_max.has_value()) {
    if (*requested_taylor_complete_max < retained_taylor_complete_max ||
        *requested_taylor_complete_max > 2000)
      throw std::invalid_argument(
          "requested terminal composed Taylor order must lie between the retained order and 2000");
    taylor_complete_max = *requested_taylor_complete_max;
  } else if (const auto* raw_order = std::getenv(
                 "DIFFEXP_DIAGNOSTIC_TERMINAL_COMPOSED_TAYLOR_ORDER")) {
    char* end = nullptr;
    const auto parsed = std::strtoul(raw_order, &end, 10);
    if (raw_order[0] == '\0' || end == raw_order || *end != '\0' ||
        parsed < retained_taylor_complete_max || parsed > 2000)
      throw std::invalid_argument(
          "DIFFEXP_DIAGNOSTIC_TERMINAL_COMPOSED_TAYLOR_ORDER must be an integer between the retained order and 2000");
    taylor_complete_max = static_cast<std::uint32_t>(parsed);
  }
  const auto physical_equation =
      match.terminal_acb_receiving_physical_equation();
  const auto exact_physical_equation =
      match.terminal_rational_shadow_physical_equation();
  if (!exact_physical_equation &&
      std::dynamic_pointer_cast<CompositeSCCChartBase>(
          owners.front()->retained_equation_owner()))
    throw std::domain_error(
        "terminal composed adjoint composite owner has no bound exact Rational-shadow physical equation");
  std::optional<PreparedSparseLocalMultiplierMatrix<Rational>> exact_row;
  if (exact_physical_equation) {
    PreparedSparseLocalMultiplierMatrix<Rational> parsed;
    parsed.rows = 1;
    parsed.columns = row.columns;
    parsed.exact_identity = row.exact_identity;
    const auto& raw_entries = as_array(
        prepared_row.at("entries"),
        "terminal composed exact rational-row entries");
    if (raw_entries.size() != row.entries.size())
      throw std::logic_error(
          "terminal composed exact/numeric rational rows differ in size");
    parsed.entries.reserve(raw_entries.size());
    for (std::size_t index = 0; index < raw_entries.size(); ++index) {
      const auto& raw_entry = as_object(
          raw_entries[index], "terminal composed exact rational-row entry");
      const auto column = as_u32(
          raw_entry.at("column"),
          "terminal composed exact rational-row column");
      if (column != row.entries[index].column)
        throw std::logic_error(
            "terminal composed exact/numeric rational-row columns differ");
      const auto& raw_multiplier = as_object(
          raw_entry.at("multiplier"),
          "terminal composed exact rational-row multiplier");
      auto multiplier =
          parse_prepared_rational_taylor_multiplier<Rational>(
              raw_multiplier,
              row.entries[index].multiplier.kernels.size(),
              prototype.taylor_width(), false,
              "terminal composed exact rational-row multiplier");
      parsed.entries.push_back({0, column, std::move(multiplier)});
    }
    exact_row = std::move(parsed);
  }
  if (taylor_complete_max > retained_taylor_complete_max) {
    if (!exact_row.has_value())
      throw std::domain_error(
          "terminal composed adjoint cannot extend its Taylor order without an exact row");
    row = adjoint_observable_detail::
        specialize_exact_backward_adjoint_row_taylor(
            *exact_row, taylor_complete_max,
            "terminal composed extended exact row:" + arm_name + ":" +
                std::to_string(tile_index));
  }
  std::shared_ptr<const
      adjoint_observable_detail::NormalizedBackwardAdjointExactODE>
      normalized_exact_equation;
  if (exact_physical_equation)
    normalized_exact_equation =
        match.terminal_normalized_backward_adjoint_exact_equation(
            exact_physical_equation, taylor_complete_max, row_complete,
            "terminal composed normalized adjoint:" + arm_name + ":" +
                std::to_string(tile_index));
  const auto* recurrence_exact_equation = normalized_exact_equation
      ? &normalized_exact_equation->truncated_ode
      : exact_physical_equation.get();
  const auto solve_context =
      "terminal composed adjoint:" + arm_name + ":" +
      std::to_string(tile_index);
  auto reservoir =
      solve_backward_adjoint_taylor_with_epsilon_reservoir(
          [&](std::int32_t input_complete_max) {
            return normalized_exact_equation
                ? prepare_backward_adjoint_taylor_problem(
                      normalized_exact_equation->truncated_ode, row,
                      taylor_complete_max, input_complete_max,
                      required_adjoint_complete, oriented_jacobian,
                      solve_context)
                : prepare_backward_adjoint_taylor_problem(
                      *physical_equation, row, taylor_complete_max,
                      input_complete_max, required_adjoint_complete,
                      oriented_jacobian, solve_context);
          },
          required_adjoint_complete, row_complete,
          recurrence_exact_equation, solve_context);
  const auto& solved = reservoir.result;
  const auto signed_point =
      local_algebra_detail::signed_real_evaluation_ball(
          certified_tile.local_begin);
  if (!signed_point.has_value())
    throw std::domain_error(
        "terminal composed adjoint match point has no signed real ball");
  auto logarithm = local_detail::cb_log(certified_tile.local_begin.modulus);
  if (certified_tile.local_begin.sign < 0) {
    const auto rim = exact_plan_rim(
        binding.prescriptions, binding.geometry.scale);
    if (!rim.has_value())
      throw std::domain_error(
          "terminal composed adjoint negative match point has no exact rim prescription");
    logarithm += local_detail::imaginary_pi(*rim);
  }
  const auto adjoint = evaluate_backward_adjoint_taylor(
      solved, *signed_point, logarithm,
      "terminal composed adjoint evaluation:" + arm_name + ":" +
          std::to_string(tile_index));
  auto point_modulus = tile.local_begin.sign() < 0
      ? -tile.local_begin : tile.local_begin;
  if (!(point_modulus < binding.geometry.radius))
    throw std::domain_error(
        "terminal composed adjoint match point is not inside the chart radius");
  std::cerr
      << "terminal-composed-adjoint-geometry arm=" << arm_name
      << " tile=" << tile_index
      << " planning_local_begin=" << tile.local_begin.str()
      << " certified_local_midpoint="
      << signed_point->real_midpoint(32)
      << " certified_local_radius_2exp="
      << signed_point->real_radius_exponent()
      << " exact_point_modulus=" << point_modulus.str()
      << " exact_chart_radius=" << binding.geometry.radius.str()
      << '\n';
  auto contracted = contract_backward_adjoint(
      adjoint, incoming,
      "terminal composed adjoint contraction:" + arm_name + ":" +
          std::to_string(tile_index));
  const auto certified_output_min = std::max(
      contracted.min_power(), epsilon_contract.requested.min_power);
  const auto certified_output_max = std::min(
      contracted.complete_max(), epsilon_contract.required_complete_max);
  if (certified_output_min > certified_output_max)
    throw std::domain_error(
        "terminal composed adjoint has no requested output epsilon window");
  std::vector<Magnitude> output_tails(
      EpsilonWindow{certified_output_min, certified_output_max}.width(),
      Magnitude::zero());
  auto recurrence_contraction = Magnitude::zero();
  auto coefficientwise_first_input_complete_max =
      std::numeric_limits<std::int32_t>::max();
  auto coefficientwise_last_input_complete_max =
      std::numeric_limits<std::int32_t>::min();
  auto coefficientwise_min_input_complete_max =
      std::numeric_limits<std::int32_t>::max();
  auto coefficientwise_max_input_complete_max =
      std::numeric_limits<std::int32_t>::min();
  for (std::int64_t raw_output = certified_output_min;
       raw_output <= certified_output_max; ++raw_output) {
    const auto output_power = static_cast<std::int32_t>(raw_output);
    const auto output_required_adjoint_complete =
        local_algebra_detail::checked_i32(
            raw_output - static_cast<std::int64_t>(incoming_min),
            "terminal composed coefficientwise adjoint epsilon maximum");
    const auto output_context =
        solve_context + ": output epsilon " +
        std::to_string(output_power);
    // The maximal private-reservoir solve is a formal Laurent solution.
    // Its lower epsilon prefixes are identical to independently re-solving
    // each smaller request; higher input coefficients cannot change a
    // coefficient that is already inside the honest complete edge.  Reuse
    // that one solution and expose only the prefix cap needed by this output.
    std::optional<BackwardAdjointReservoirSolve> output_reservoir;
    const BackwardAdjointTaylorProblem* output_problem =
        &reservoir.problem;
    const BackwardAdjointTaylorResult* output_solution = &solved;
    FiniteLaurentVector<ComplexBall> output_adjoint_storage;
    const FiniteLaurentVector<ComplexBall>* output_adjoint = &adjoint;
    std::int32_t output_input_complete_max;
    if (normalized_exact_equation && exact_row) {
      output_input_complete_max =
          backward_adjoint_prefix_input_complete_max(
              reservoir, output_required_adjoint_complete,
              output_context + ": maximal-reservoir slice");
    } else {
      output_reservoir =
          solve_backward_adjoint_taylor_with_epsilon_reservoir(
              [&](std::int32_t input_complete_max) {
                return prepare_backward_adjoint_taylor_problem(
                    *physical_equation, row, taylor_complete_max,
                    input_complete_max,
                    output_required_adjoint_complete,
                    oriented_jacobian, output_context);
              },
              output_required_adjoint_complete, row_complete,
              recurrence_exact_equation, output_context);
      output_input_complete_max =
          output_reservoir->input_epsilon_complete_max;
      output_problem = &output_reservoir->problem;
      output_solution = &output_reservoir->result;
      output_adjoint_storage = evaluate_backward_adjoint_taylor(
          *output_solution, *signed_point, logarithm,
          output_context + ": evaluation");
      output_adjoint = &output_adjoint_storage;
    }
    if (raw_output == certified_output_min)
      coefficientwise_first_input_complete_max =
          output_input_complete_max;
    coefficientwise_last_input_complete_max =
        output_input_complete_max;
    coefficientwise_min_input_complete_max = std::min(
        coefficientwise_min_input_complete_max,
        output_input_complete_max);
    coefficientwise_max_input_complete_max = std::max(
        coefficientwise_max_input_complete_max,
        output_input_complete_max);
    Magnitude adjoint_tail;
    if (normalized_exact_equation && exact_row) {
      const auto real_tail = certify_backward_adjoint_real_ray_tail(
          *normalized_exact_equation, *exact_row,
          solved, oriented_jacobian,
          tile.local_begin,
          output_input_complete_max, 128,
          "terminal composed real-ray adjoint tail:" + arm_name + ":" +
              std::to_string(tile_index) + ": output epsilon " +
              std::to_string(output_power),
          &match.terminal_backward_adjoint_real_ray_operator_cache());
      adjoint_tail = real_tail.absolute_vector_tail_upper;
      const auto stability_ratio = real_tail.operator_norm_upper /
          Magnitude::from_ui(
              static_cast<ulong>(taylor_complete_max) + 1);
      recurrence_contraction = Magnitude::maximum(
          recurrence_contraction, stability_ratio);
    } else {
      const auto disk_tail =
          certify_backward_adjoint_taylor_tail_adaptive_witness(
              *output_problem, *output_solution, row,
              oriented_jacobian, *signed_point, point_modulus,
              binding.geometry.radius, 128,
              "terminal composed adjoint tail:" + arm_name + ":" +
                  std::to_string(tile_index) + ": output epsilon " +
                  std::to_string(output_power),
              output_input_complete_max,
              exact_row ? &*exact_row : nullptr,
              recurrence_exact_equation);
      adjoint_tail = disk_tail.tail.absolute_vector_tail_upper;
      recurrence_contraction = Magnitude::maximum(
          recurrence_contraction,
          disk_tail.tail.recurrence_contraction_upper);
    }
    const auto contracted_tail =
        backward_adjoint_contracted_tail_by_output(
            adjoint_tail,
            *output_adjoint, incoming,
            {output_power, output_power},
            "terminal composed adjoint coefficientwise tail:" + arm_name +
                ":" + std::to_string(tile_index) + ": output epsilon " +
                std::to_string(output_power));
    output_tails[static_cast<std::size_t>(
        raw_output - certified_output_min)] = contracted_tail.front();
  }
  auto output_tail = Magnitude::zero();
  std::vector<Magnitude> output_scales;
  output_scales.reserve(output_tails.size());
  std::vector<ComplexBall> widened_coefficients;
  widened_coefficients.reserve(output_tails.size());
  for (std::size_t offset = 0; offset < output_tails.size();
       ++offset) {
    const auto power = local_algebra_detail::checked_i32(
        static_cast<std::int64_t>(certified_output_min) +
            static_cast<std::int64_t>(offset),
        "terminal composed requested output power");
    auto coefficient = contracted.coefficient(power);
    output_scales.push_back(Magnitude::maximum(
        Magnitude::one(), Magnitude::lower_abs(coefficient)));
    output_tails[offset].add_error_to(coefficient);
    widened_coefficients.push_back(std::move(coefficient));
    output_tail = Magnitude::maximum(output_tail, output_tails[offset]);
  }
  contracted = EpsilonFrame<ComplexBall>(certified_output_min,
                                         std::move(widened_coefficients));
  return {std::move(contracted), std::move(output_tails),
          std::move(output_scales), exact_physical_equation != nullptr,
          taylor_complete_max,
          reservoir.input_epsilon_complete_max,
          solved.common_epsilon_complete_max,
          coefficientwise_first_input_complete_max,
          coefficientwise_last_input_complete_max,
          coefficientwise_min_input_complete_max,
          coefficientwise_max_input_complete_max,
          solved.last_coefficient_ratio_upper, output_tail,
          recurrence_contraction};
}

StoredLineIntegral integrate_transport_terminal_factorized_acb_row_tile(
    slong precision_bits,
    const std::shared_ptr<StoredPlannedMatchHop>& terminal_match,
    std::size_t expected_tile_count,
    const json::object& prepared_row,
    const RetainedArmPlan& arm, const std::string& arm_name,
    std::size_t tile_index,
    const ObservableEpsilonContract& epsilon_contract,
    const std::optional<BoundedDivergentCancellation>&
        divergent_cancellation) {
  if (!terminal_match ||
      tile_index + 1 != expected_tile_count ||
      tile_index >= arm.exact.tiles.size())
    throw std::invalid_argument(
        "terminal factorized line integration lost its final tile or match owner");
  AcbPrecisionLease lease(precision_bits);
  ComplexBall::set_precision(precision_bits);
  const auto& match = *terminal_match;
  if (!match.has_terminal_acb_factorization())
    throw std::invalid_argument(
        "terminal factorized line integration has no certified Acb factorization");
  matching_detail::ScopedAcbPrecision exact_shadow_precision(
      match.terminal_acb_extra_precision_bits());
  const auto& tile = arm.exact.tiles[tile_index];
  const auto& certified_tile = arm.certified_tiles.at(tile_index);
  const auto& binding = arm.charts.at(tile.chart);
  const auto terminal_started = std::chrono::steady_clock::now();

  const auto* diagnostic_route_raw =
      std::getenv("DIFFEXP_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE");
  const std::string diagnostic_route =
      diagnostic_route_raw == nullptr
      ? std::string()
      : std::string(diagnostic_route_raw);
  if (!diagnostic_route.empty() &&
      diagnostic_route != "factorized" &&
      diagnostic_route != "physical" &&
      diagnostic_route != "adjoint" &&
      diagnostic_route != "compare")
    throw std::invalid_argument(
        "DIFFEXP_DIAGNOSTIC_TERMINAL_CONTRACTION_ROUTE must be "
        "factorized, physical, adjoint, or compare");
  bool direct_physical =
      diagnostic_route == "physical";
  const bool direct_factorized =
      diagnostic_route == "factorized";
  const bool compare_factorized =
      diagnostic_route == "compare";
  std::string contraction_provenance =
      diagnostic_route == "physical"
      ? "terminal direct-physical diagnostic contraction"
      : diagnostic_route == "factorized"
      ? "terminal direct-factorized diagnostic contraction"
      : diagnostic_route == "compare"
      ? "terminal factorized-adjoint production contraction with "
        "direct-factorized comparison"
      : "terminal factorized-adjoint production contraction";
  std::optional<TerminalComposedAdjointDiagnostic> composed_diagnostic;
  bool composed_authoritative = false;
  const auto composed_started = std::chrono::steady_clock::now();
  if (const auto* composed_mode =
          std::getenv("DIFFEXP_DIAGNOSTIC_TERMINAL_COMPOSED_ADJOINT")) {
    const std::string mode(composed_mode);
    if (mode != "report" && mode != "require" &&
        mode != "authoritative")
      throw std::invalid_argument(
          "DIFFEXP_DIAGNOSTIC_TERMINAL_COMPOSED_ADJOINT must be report, require, or authoritative");
    if (certified_tile.local_end.sign != 0) {
      std::cerr
          << "terminal-composed-adjoint arm=" << arm_name
          << " tile=" << tile_index
          << " status=not-applicable detail=tile-does-not-end-at-chart-center"
          << " local_interval=[" << tile.local_begin.str() << ","
          << tile.local_end.str() << "]\n";
    } else try {
      std::optional<std::uint32_t> requested_taylor_complete_max;
      const auto learned_taylor_order =
          match.terminal_composed_taylor_order_floor();
      if (learned_taylor_order > 0)
        requested_taylor_complete_max = learned_taylor_order;
      const auto relative_tolerance = Magnitude::decimal(
          match.terminal_acb_relative_tolerance_text());
      const auto selection_tolerance =
          terminal_composed_adjoint_selection_tolerance(
              relative_tolerance);
      for (;;) {
        composed_diagnostic = compute_terminal_composed_adjoint_diagnostic(
            match, prepared_row, arm, arm_name, tile_index,
            epsilon_contract, requested_taylor_complete_max);
        if (mode != "authoritative")
          break;
        if (terminal_composed_adjoint_meets_relative_accuracy(
                *composed_diagnostic, selection_tolerance)) {
          match.learn_terminal_composed_taylor_order_floor(
              composed_diagnostic->taylor_complete_max);
          break;
        }
        if (!composed_diagnostic->can_extend_taylor_exactly) {
          std::cerr
              << "terminal-composed-adjoint-selected arm=" << arm_name
              << " tile=" << tile_index
              << " status=authoritative-fallback-insufficient-accuracy"
              << " detail=no-exact-equation-for-taylor-extension"
              << " relative_tail_upper="
              << terminal_composed_adjoint_relative_tail_upper(
                     *composed_diagnostic).approximate_upper()
              << " relative_tolerance="
              << match.terminal_acb_relative_tolerance_text()
              << " selection_tolerance_factor=1e-8\n";
          break;
        }
        const auto completed = composed_diagnostic->taylor_complete_max;
        if (completed >= 2000)
          throw std::domain_error(
              "terminal composed adjoint did not meet the retained match accuracy by Taylor order 2000");
        const auto next = static_cast<std::uint32_t>(std::min<std::uint64_t>(
            2000, static_cast<std::uint64_t>(completed) +
                      std::max<std::uint64_t>(
                          32,
                          static_cast<std::uint64_t>(completed) / 2)));
        std::cerr
            << "terminal-composed-adjoint-retry arm=" << arm_name
            << " tile=" << tile_index
            << " reason=insufficient-authoritative-accuracy"
            << " completed_taylor=" << completed
            << " next_taylor=" << next
            << " relative_tail_upper="
            << terminal_composed_adjoint_relative_tail_upper(
                   *composed_diagnostic).approximate_upper()
            << " relative_tolerance="
            << match.terminal_acb_relative_tolerance_text()
            << " selection_tolerance_factor=1e-8\n";
        requested_taylor_complete_max = next;
      }
      composed_authoritative =
          mode == "authoritative" &&
          terminal_composed_adjoint_meets_relative_accuracy(
              *composed_diagnostic, selection_tolerance);
      std::cerr
          << "terminal-composed-adjoint arm=" << arm_name
          << " tile=" << tile_index
          << " status=certified-composed-tail"
          << " taylor_complete_max="
          << composed_diagnostic->taylor_complete_max
          << " adjoint_epsilon_input_complete_max="
          << composed_diagnostic->adjoint_epsilon_input_complete_max
          << " adjoint_epsilon_complete_max="
          << composed_diagnostic->adjoint_epsilon_complete_max
          << " coefficientwise_first_input_complete_max="
          << composed_diagnostic->coefficientwise_first_input_complete_max
          << " coefficientwise_last_input_complete_max="
          << composed_diagnostic->coefficientwise_last_input_complete_max
          << " coefficientwise_input_complete_range=["
          << composed_diagnostic->coefficientwise_min_input_complete_max
          << ","
          << composed_diagnostic->coefficientwise_max_input_complete_max
          << "]"
          << " last_coefficient_ratio_upper=";
      if (composed_diagnostic->last_coefficient_ratio_upper.has_value())
        std::cerr << *composed_diagnostic->last_coefficient_ratio_upper;
      else
        std::cerr << "unavailable";
      std::cerr
          << " certified_output_tail_upper="
          << composed_diagnostic->certified_output_tail_upper.approximate_upper()
          << " relative_output_tail_upper="
          << terminal_composed_adjoint_relative_tail_upper(
                 *composed_diagnostic).approximate_upper()
          << " recurrence_contraction_upper="
          << composed_diagnostic->recurrence_contraction_upper.approximate_upper()
          << '\n';
      if (composed_authoritative) {
        contraction_provenance =
            "terminal certified composed backward-adjoint production contraction";
        std::cerr
            << "terminal-composed-adjoint-selected arm=" << arm_name
            << " tile=" << tile_index
            << " status=authoritative-certified-value\n";
      }
    } catch (const BackwardAdjointCenterAnchoringError& error) {
      std::cerr
          << "terminal-composed-adjoint arm=" << arm_name
          << " tile=" << tile_index
          << " status=not-applicable"
          << " detail=row-requires-laurent-log-center-adjoint"
          << " forcing_power=" << error.forcing_power << '\n';
      if (mode == "authoritative" &&
          diagnostic_route.empty()) {
        direct_physical = true;
        contraction_provenance =
            "terminal certified direct-physical fallback for "
            "Laurent/log center-adjoint row";
        std::cerr
            << "terminal-composed-adjoint-selected arm=" << arm_name
            << " tile=" << tile_index
            << " status=classified-direct-physical-fallback"
            << " detail=row-requires-laurent-log-center-adjoint\n";
      }
    } catch (const BackwardAdjointCenterUnitError& error) {
      std::cerr
          << "terminal-composed-adjoint arm=" << arm_name
          << " tile=" << tile_index
          << " status=not-applicable"
          << " detail=q-epsilon0-is-not-center-unit"
          << " q_t_valuation=";
      if (error.t_valuation.has_value())
        std::cerr << *error.t_valuation;
      else
        std::cerr << "absent";
      std::cerr << '\n';
      if (mode == "authoritative" &&
          diagnostic_route.empty()) {
        direct_physical = true;
        contraction_provenance =
            "terminal certified direct-physical fallback for "
            "nonunit center-adjoint row";
        std::cerr
            << "terminal-composed-adjoint-selected arm=" << arm_name
            << " tile=" << tile_index
            << " status=classified-direct-physical-fallback"
            << " detail=q-epsilon0-is-not-center-unit\n";
      }
    } catch (const std::exception& error) {
      if (mode != "report") throw;
      std::cerr
          << "terminal-composed-adjoint arm=" << arm_name
          << " tile=" << tile_index
          << " status=unsupported detail=" << error.what();
      try {
        std::cerr << " indicial="
                  << terminal_composed_adjoint_indicial_summary(match);
      } catch (const std::exception& diagnostic_error) {
        std::cerr << " indicial=unavailable:" << diagnostic_error.what();
      }
      std::cerr << '\n';
    }
  }
  const auto composed_ms = elapsed_milliseconds(composed_started);
  const std::optional<std::int32_t> projection_cap =
      direct_physical || direct_factorized
      ? std::optional<std::int32_t>(terminal_factorized_physical_complete_cap(
            match, epsilon_contract.requested.complete_max, 1,
            "terminal factorized primitive window",
            direct_physical))
      : std::nullopt;
  const auto projection_started = std::chrono::steady_clock::now();
  auto projected = direct_physical
      ? project_terminal_acb_physical_basis_row(
            match, prepared_row, *projection_cap,
            "terminal-physical-line:" + arm_name + ":" +
                std::to_string(tile_index))
      : project_terminal_acb_basis_row(
            match, prepared_row, projection_cap,
            "terminal-factorized-line:" + arm_name + ":" +
                std::to_string(tile_index));
  const auto projection_ms = elapsed_milliseconds(projection_started);

  const auto rim = exact_plan_rim(
      binding.prescriptions, binding.geometry.scale);
  const bool reverse_local_orientation =
      tile.local_end < tile.local_begin;
  // `physical_rows` below are integrated in monotonically increasing local
  // t.  The published line is instead oriented in physical x, so this same
  // Jacobian belongs both in the final result and in any comparison with the
  // composed adjoint (whose forcing already contains dx = beta dt).
  const auto output_jacobian = reverse_local_orientation
      ? -binding.scale_numeric : binding.scale_numeric;
  const auto& primitive_begin = reverse_local_orientation
      ? certified_tile.local_end : certified_tile.local_begin;
  const auto& primitive_end = reverse_local_orientation
      ? certified_tile.local_begin : certified_tile.local_end;
  if (certified_tile_zero_length(arm, tile_index)) {
    auto zero = certified_zero_physical_line(
        epsilon_contract.requested, 1, rim,
        certified_tile.local_begin.sign == 0, false);
    zero.diagnostics.detail += "; selected " + contraction_provenance;
    return zero;
  }

  std::vector<StoredLineIntegral> physical;
  physical.reserve(projected.size());
  const auto physical_integration_started =
      std::chrono::steady_clock::now();
  for (std::size_t column = 0; column < projected.size(); ++column) {
    const auto& source = projected[column];
    const auto delivered_min =
        local_algebra_detail::checked_i32(
            static_cast<std::int64_t>(
                source.epsilon.min_power) - 1,
            "terminal physical primitive minimum");
    const auto delivered_max =
        local_algebra_detail::checked_i32(
            static_cast<std::int64_t>(
                source.epsilon.complete_max) - 1,
            "terminal physical primitive maximum");
    if (delivered_max < delivered_min)
      throw std::domain_error(
          "terminal physical column has no primitive epsilon window");
    StoredLineIntegrationOptions options;
    options.delivered_epsilon =
        {delivered_min, delivered_max};
    options.required_complete_max = delivered_min;
    options.imaginary_sign = rim;
    options.certified_chart_scale_sign = binding.scale_sign;
    options.divergent_cancellation = divergent_cancellation;
    options.defer_divergent_groups = true;
    try {
      physical.push_back(integrate_stored_local_line(
          source, primitive_begin, primitive_end, options));
    } catch (const NativeIntegrationError& error) {
      std::ostringstream detail;
      detail << error.what() << "; arm=" << arm_name
             << "; tile=" << tile_index
             << "; physical_basis_column=" << column
             << "; terminal_factorized=true";
      NativeIntegrationError contextual(
          error.code, error.id, detail.str());
      contextual.absolute_power = error.absolute_power;
      contextual.log_power = error.log_power;
      contextual.epsilon_power = error.epsilon_power;
      contextual.component = error.component;
      throw contextual;
    }
  }
  const auto physical_integration_ms =
      elapsed_milliseconds(physical_integration_started);

  using DeferredKey = line_integration_detail::MonomialKey;
  std::map<DeferredKey, SectorMonomialTag> deferred_tags;
  std::vector<std::map<
      DeferredKey,
      const StoredLineIntegral::DeferredDivergentGroup*>>
      deferred_by_column(physical.size());
  FiniteLaurentMatrix<ComplexBall> physical_rows(1);
  physical_rows.front().reserve(physical.size());
  for (std::size_t column = 0; column < physical.size(); ++column) {
    physical_rows.front().push_back(scalar_epsilon_frame(
        physical[column].value,
        "terminal physical line column"));
    for (const auto& group :
         physical[column].deferred_divergent_groups) {
      const auto key =
          line_integration_detail::monomial_key(group.tag);
      deferred_tags.try_emplace(key, group.tag);
      deferred_by_column[column].emplace(key, &group);
    }
  }
  for (const auto& [key, tag] : deferred_tags) {
    (void)tag;
    FiniteLaurentVector<ComplexBall> row;
    row.reserve(physical.size());
    for (std::size_t column = 0; column < physical.size();
         ++column) {
      const auto found = deferred_by_column[column].find(key);
      if (found != deferred_by_column[column].end()) {
        row.push_back(found->second->coefficients);
      } else {
        row.emplace_back(
            projected[column].epsilon,
            std::vector<ComplexBall>(
                projected[column].epsilon.width(),
                ComplexBall(0)));
      }
    }
    physical_rows.push_back(std::move(row));
  }

  // Production projects G=(F*T)*P, evaluates the full endpoint map through
  // A^T z = J^T, and contracts z^T u.  This certifies the observable itself
  // instead of assuming an interior backward residual remains harmless at a
  // singular endpoint.  Direct physical/factorized weight contractions are
  // explicit diagnostics only.
  FiniteLaurentVector<ComplexBall> contracted;
  const auto functional_contraction_started =
      std::chrono::steady_clock::now();
  if (composed_authoritative) {
    if (!composed_diagnostic.has_value())
      throw std::logic_error(
          "authoritative terminal composed adjoint lost its certified value");
    // The composed equation already includes dx=beta*dt and therefore
    // returns the physically oriented finite line.  `contracted` below is
    // stored before the common output-Jacobian application, so undo exactly
    // that final scalar here.  Divergent endpoint tags are separate rows;
    // retain their established factorized-adjoint cancellation certificate.
    contracted.push_back(composed_diagnostic->value.scaled(
        ComplexBall(1) / output_jacobian));
    if (physical_rows.size() > 1) {
      FiniteLaurentMatrix<ComplexBall> deferred_rows(
          physical_rows.begin() + 1, physical_rows.end());
      auto deferred = match.adjoint_contract_terminal_acb_functionals(
          deferred_rows,
          epsilon_contract.required_complete_max,
          "terminal-factorized-adjoint-deferred-contraction:" + arm_name +
              ":" + std::to_string(tile_index),
          true);
      contracted.insert(
          contracted.end(),
          std::make_move_iterator(deferred.begin()),
          std::make_move_iterator(deferred.end()));
    }
  } else if (direct_physical) {
    contracted = match.contract_terminal_acb_physical_functionals(
        physical_rows,
        "terminal-direct-physical-line-contraction:" + arm_name + ":" +
            std::to_string(tile_index));
  } else if (direct_factorized) {
    contracted = match.contract_terminal_acb_factorized_functionals(
        physical_rows,
        "terminal-direct-factorized-line-contraction:" + arm_name + ":" +
            std::to_string(tile_index));
  } else {
    contracted = match.adjoint_contract_terminal_acb_functionals(
        physical_rows,
        epsilon_contract.required_complete_max,
        "terminal-factorized-adjoint-line-contraction:" + arm_name + ":" +
            std::to_string(tile_index),
        true);
  }
  const auto functional_contraction_ms =
      elapsed_milliseconds(functional_contraction_started);
  if (composed_diagnostic.has_value() && !composed_authoritative) {
    const auto& composed = composed_diagnostic->value;
    const auto legacy = contracted.front().scaled(output_jacobian);
    const auto common_min = std::max({
        epsilon_contract.requested.min_power,
        composed.min_power(), legacy.min_power()});
    const auto common_max = std::min({
        epsilon_contract.requested.complete_max,
        composed.complete_max(), legacy.complete_max()});
    for (std::int64_t raw_power = common_min;
         raw_power <= common_max; ++raw_power) {
      const auto power = static_cast<std::int32_t>(raw_power);
      const auto discrepancy =
          composed.coefficient(power) - legacy.coefficient(power);
      std::cerr
          << "terminal-composed-adjoint-compare arm=" << arm_name
          << " tile=" << tile_index
          << " epsilon_power=" << power
          << " composed_midpoint=("
          << composed.coefficient(power).real_midpoint(16) << ","
          << composed.coefficient(power).imag_midpoint(16) << ")"
          << " legacy_midpoint=("
          << legacy.coefficient(power).real_midpoint(16) << ","
          << legacy.coefficient(power).imag_midpoint(16) << ")"
          << " discrepancy_upper="
          << Magnitude::upper_abs(discrepancy).approximate_upper()
          << '\n';
    }
  }
  if (compare_factorized) {
    const auto direct =
        match.contract_terminal_acb_factorized_functionals(
            physical_rows,
            "terminal-compare-direct-factorized-line-contraction:" +
                arm_name + ":" + std::to_string(tile_index));
    if (direct.size() != contracted.size())
      throw std::logic_error(
          "terminal direct/adjoint comparison changed its tagged row count");
    for (std::size_t functional = 0;
         functional < contracted.size(); ++functional) {
      const auto common_min = std::max({
          epsilon_contract.requested.min_power,
          direct[functional].min_power(),
          contracted[functional].min_power()});
      const auto common_max = std::min({
          epsilon_contract.requested.complete_max,
          direct[functional].complete_max(),
          contracted[functional].complete_max()});
      if (common_min > common_max)
        throw std::domain_error(
            "terminal direct/adjoint comparison has no common requested "
            "epsilon window");
      for (std::int64_t raw_power = common_min;
           raw_power <= common_max; ++raw_power) {
        const auto power = static_cast<std::int32_t>(raw_power);
        const auto direct_value = direct[functional].coefficient(power);
        const auto adjoint_value =
            contracted[functional].coefficient(power);
        const auto discrepancy = adjoint_value - direct_value;
        std::cerr
            << "terminal-contraction-compare arm=" << arm_name
            << " tile=" << tile_index
            << " functional=" << functional
            << " epsilon_power=" << power
            << " direct_midpoint=("
            << direct_value.real_midpoint(16) << ","
            << direct_value.imag_midpoint(16) << ")"
            << " adjoint_midpoint=("
            << adjoint_value.real_midpoint(16) << ","
            << adjoint_value.imag_midpoint(16) << ")"
            << " direct_radius2exp=("
            << direct_value.real_radius_exponent() << ","
            << direct_value.imag_radius_exponent() << ")"
            << " adjoint_radius2exp=("
            << adjoint_value.real_radius_exponent() << ","
            << adjoint_value.imag_radius_exponent() << ")"
            << " discrepancy_radius2exp=("
            << discrepancy.real_radius_exponent() << ","
            << discrepancy.imag_radius_exponent() << ")"
            << " discrepancy_upper="
            << Magnitude::upper_abs(discrepancy).approximate_upper()
            << " direct_upper="
            << Magnitude::upper_abs(direct_value).approximate_upper()
            << " adjoint_upper="
            << Magnitude::upper_abs(adjoint_value).approximate_upper()
            << '\n';
      }
    }
  }
  if (contracted.size() != 1 + deferred_tags.size())
    throw std::logic_error(
        "terminal factorized line contraction changed its tagged row count");

  StoredLineIntegral result;
  for (const auto& column : physical)
    accumulate_terminal_line_diagnostics(
        result.diagnostics, column.diagnostics);
  result.diagnostics.divergent_cancellation_mode =
      divergent_cancellation.has_value()
      ? "bounded-relative-acb" : "exact-singleton";
  if (divergent_cancellation.has_value()) {
    result.diagnostics.divergent_relative_tolerance =
        divergent_cancellation->relative_tolerance_text;
    result.diagnostics.divergent_cancellation_provenance =
        divergent_cancellation->provenance;
  }

  std::size_t contracted_index = 1;
  // The deferred groups have the same column coordinates as `projected`.
  // These deferred groups are J's per-G-column contributions.  Even though
  // production evaluates the equivalent adjoint expression z^T u, its
  // cancellation scale must therefore use the transformed G coordinates y.
  // A direct-physical route projects F and uses physical weights w.
  const auto& selected_divergent_weights = direct_physical
      ? match.terminal_acb_physical_weights()
      : match.terminal_acb_transformed_weights();
  for (const auto& [key, tag] : deferred_tags) {
    const auto& frame = contracted[contracted_index++];
    std::vector<
        const StoredLineIntegral::DeferredDivergentGroup*>
        physical_groups(physical.size(), nullptr);
    for (std::size_t column = 0; column < physical.size();
         ++column) {
      const auto found = deferred_by_column[column].find(key);
      if (found != deferred_by_column[column].end())
        physical_groups[column] = found->second;
    }
    const bool had_material = std::any_of(
        physical_groups.begin(), physical_groups.end(),
        [](const auto* group) {
          return group != nullptr && group->had_material_input;
        });
    for (std::int64_t raw_power = frame.min_power();
         raw_power <= frame.complete_max(); ++raw_power) {
      const auto power = static_cast<std::int32_t>(raw_power);
      const auto coefficient = frame.coefficient(power);
      if (coefficient.is_zero()) continue;
      bool bounded = false;
      Magnitude scale = Magnitude::one();
      Magnitude bound = Magnitude::zero();
      const auto coefficient_upper =
          Magnitude::upper_abs(coefficient);
      if (divergent_cancellation.has_value()) {
        scale = Magnitude::maximum(
            Magnitude::one(),
            terminal_factorized_divergent_scale(
                selected_divergent_weights, physical_groups, power));
        bound = scale *
            divergent_cancellation->relative_tolerance;
        bounded = coefficient_upper <= bound;
      }
      if (bounded) {
        ++result.diagnostics
              .bounded_cancelled_divergent_coefficients;
        continue;
      }
      const bool uncertified = coefficient.contains_zero();
      std::ostringstream detail;
      detail << (uncertified
          ? "factorized terminal divergent coefficient contains zero but is not certified by the active cancellation policy"
          : "factorized terminal divergent coefficient is nonzero");
      if (divergent_cancellation.has_value())
        detail << "; coefficient_upper="
               << coefficient_upper.dump_exact()
               << "; contribution_scale_upper="
               << scale.dump_exact()
               << "; relative_tolerance="
               << divergent_cancellation
                      ->relative_tolerance_text
               << "; tolerance_bound_upper="
               << bound.dump_exact();
      NativeIntegrationError error(
          uncertified
              ? NativeIntegrationErrorCode::UncertifiedCancellation
              : NativeIntegrationErrorCode::DivergentEndpoint,
          uncertified ? "E10" : "E2", detail.str());
      error.absolute_power = tag.m.canonical;
      error.log_power = tag.log_power;
      error.epsilon_power = power;
      error.component = 0;
      throw error;
    }
    if (had_material)
      ++result.diagnostics.cancelled_divergent_groups;
  }

  const auto& finite = contracted.front();
  const auto complete_max = std::min(
      epsilon_contract.requested.complete_max,
      finite.complete_max());
  const auto min_power = std::max(
      epsilon_contract.requested.min_power,
      finite.min_power());
  if (min_power > complete_max ||
      epsilon_contract.required_complete_max > complete_max) {
    std::ostringstream detail;
    detail
        << "terminal factorized line does not cover its required output "
           "epsilon window; requested=["
        << epsilon_contract.requested.min_power << ","
        << epsilon_contract.requested.complete_max
        << "]; required_complete_max="
        << epsilon_contract.required_complete_max
        << "; projection_cap=";
    if (projection_cap.has_value())
      detail << *projection_cap;
    else
      detail << "all-available";
    detail
        << "; contracted=[" << finite.min_power() << ","
        << finite.complete_max() << "]; projected=[";
    for (std::size_t column = 0; column < projected.size();
         ++column) {
      if (column != 0) detail << ",";
      detail << projected[column].epsilon.min_power << ":"
             << projected[column].epsilon.complete_max;
    }
    detail << "]; physical_integrals=[";
    for (std::size_t column = 0; column < physical.size();
         ++column) {
      if (column != 0) detail << ",";
      detail << physical[column].value.epsilon.min_power << ":"
             << physical[column].value.epsilon.complete_max;
    }
    detail << "]; physical_weights=[";
    const auto& weights = match.terminal_acb_physical_weights();
    for (std::size_t column = 0; column < weights.size();
         ++column) {
      if (column != 0) detail << ",";
      detail << weights[column].min_power() << ":"
             << weights[column].complete_max();
    }
    detail << "]";
    throw std::domain_error(detail.str());
  }
  result.value.epsilon = {min_power, complete_max};
  result.value.dimension = 1;
  result.value.coefficients.reserve(
      result.value.epsilon.width());
  for (std::int64_t raw_power = min_power;
       raw_power <= complete_max; ++raw_power)
    result.value.coefficients.push_back(
        finite.coefficient(
            static_cast<std::int32_t>(raw_power)) *
        output_jacobian);
  result.value.error.guarantee = ErrorGuarantee::None;
  result.value.error.provenance =
      "stored Taylor truncation only; " + contraction_provenance +
      "; no unseen-tail majorant or full-local certificate";
  result.scope = LineIntegrationScope::StoredTruncation;
  result.imaginary_sign = rim;
  result.diagnostics.detail =
      result.value.error.provenance;
  if (transport_contraction_timing_enabled()) {
    std::ostringstream timing;
    timing << "transport-contraction-timing"
           << " phase=terminal-row"
           << " arm=" << arm_name
           << " tile=" << tile_index
           << " composed_ms=" << composed_ms
           << " projection_ms=" << projection_ms
           << " physical_integration_ms="
           << physical_integration_ms
           << " functional_contraction_ms="
           << functional_contraction_ms
           << " projected_columns=" << projected.size()
           << " composed_authoritative="
           << (composed_authoritative ? "true" : "false")
           << " total_ms="
           << elapsed_milliseconds(terminal_started);
    emit_transport_contraction_timing(timing.str());
  }
  return result;
}

StoredLineIntegral integrate_transport_terminal_factorized_acb_row_tile(
    slong precision_bits,
    const std::shared_ptr<StoredTransportArmState>& state,
    const json::object& prepared_row,
    const RetainedArmPlan& arm, const std::string& arm_name,
    std::size_t tile_index,
    const ObservableEpsilonContract& epsilon_contract,
    const std::optional<BoundedDivergentCancellation>&
        divergent_cancellation) {
  if (!state)
    throw std::invalid_argument(
        "terminal factorized line integration lost its transport state");
  return integrate_transport_terminal_factorized_acb_row_tile(
      precision_bits, state->terminal_factorized_match(),
      state->tile_sources().size(), prepared_row, arm, arm_name, tile_index,
      epsilon_contract, divergent_cancellation);
}

void require_transport_line_publication_accuracy(
    const StoredLineIntegral& line,
    const ObservableEpsilonContract& epsilon,
    const Magnitude& publication_relative_tolerance,
    const std::string& context,
    std::optional<std::size_t> request_index = std::nullopt) {
  const auto first = std::max(
      epsilon.requested.min_power, line.value.epsilon.min_power);
  const auto last = std::min(
      epsilon.required_complete_max,
      line.value.epsilon.complete_max);
  if (last < first)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::InsufficientCompleteWindow,
        context +
            ": paired line has no complete publication coefficient",
        std::nullopt, request_index, line.value.epsilon.complete_max);
  for (std::size_t component = 0;
       component < line.value.dimension; ++component)
    for (std::int64_t raw_power = first;
         raw_power <= last; ++raw_power) {
      const auto power = local_algebra_detail::checked_i32(
          raw_power, "paired line publication power");
      const auto coefficient = line.value.at(power, component);
      const auto publication = certify_acb_publication_accuracy(
          coefficient, publication_relative_tolerance);
      if (!publication.acceptable)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::TerminalOutputInconclusive,
            context +
                ": final paired line ball is too wide to publish; functional=" +
                std::to_string(component) +
                "; epsilon_power=" + std::to_string(power) +
                "; midpoint=(" + coefficient.real_midpoint(16) + "," +
                coefficient.imag_midpoint(16) + ")" +
                "; radius2exp=(" +
                coefficient.real_radius_exponent() + "," +
                coefficient.imag_radius_exponent() + ")",
            component, request_index, power,
            publication.required_additional_digits);
    }
}

std::shared_ptr<StoredLineResult>
build_rolling_transport_observable_line_from_result(
    const std::string& handle, const std::string& checkpoint_identity,
    const std::string& arm_name, json::object interval,
    json::object aggregate_record,
    const std::shared_ptr<StoredTilePlan>& plan,
    const std::shared_ptr<StoredLocalBase>& terminal,
    StoredLineIntegral result, double elapsed_ms) {
  if (handle.empty() || checkpoint_identity.empty() || !plan || !terminal)
    throw std::invalid_argument(
        "rolling transport observable line lost an identity or compact owner");
  json::object provenance{
      {"schema", "diffexp3-retained-native-line-aggregate-v2"},
      {"checkpoint_identity", checkpoint_identity},
      {"arm", arm_name},
      {"interval", std::move(interval)},
      {"source", json::object{
           {"tile_plan", plan->handle()},
           {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
           {"locals", bounded_line_aggregate_source_records({terminal})}}},
      {"aggregate", std::move(aggregate_record)},
      {"epsilon", json::object{{"min", result.value.epsilon.min_power},
                                {"max", result.value.epsilon.complete_max}}},
      {"scope", line_integration_scope_name(result.scope)},
      {"error_guarantee",
       error_guarantee_name(result.value.error.guarantee)}};
  const auto provenance_identity = json::serialize(
      canonical_json_value(provenance));
  return std::make_shared<StoredLineResult>(
      handle, checkpoint_identity, provenance_identity, std::move(result),
      elapsed_ms, plan,
      std::vector<std::shared_ptr<StoredLocalBase>>{terminal},
      std::move(provenance));
}

class TransportPairObservableStream final {
 public:
  enum class Status { Active, Poisoned, Finishing, Finished, Aborted };

  struct FinishResult {
    std::shared_ptr<StoredLineResult> line;
    std::array<std::size_t, 2> tiles{0, 0};
    std::array<double, 2> arm_integration_ms{0.0, 0.0};
    json::object conditioning;
    double elapsed_ms = 0.0;
  };

  TransportPairObservableStream(
      std::string handle, std::string stream_checkpoint_identity,
      std::string line_handle, std::string checkpoint_root,
      std::string domain, slong precision_bits,
      std::array<std::shared_ptr<StoredTransportArmState>, 2> states,
      std::string identity, std::string checkpoint_identity,
      ObservableEpsilonContract epsilon, json::object epsilon_record,
      TransportTailPolicy tail_policy,
      std::optional<BoundedDivergentCancellation> divergent_cancellation,
      std::optional<Magnitude> publication_relative_tolerance)
      : handle_(std::move(handle)),
        stream_checkpoint_identity_(
            std::move(stream_checkpoint_identity)),
        line_handle_(std::move(line_handle)),
        checkpoint_root_(std::move(checkpoint_root)),
        domain_(std::move(domain)), precision_bits_(precision_bits),
        states_(std::move(states)), identity_(std::move(identity)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        epsilon_(epsilon), epsilon_record_(std::move(epsilon_record)),
        tail_policy_(tail_policy),
        divergent_cancellation_(std::move(divergent_cancellation)),
        publication_relative_tolerance_(
            std::move(publication_relative_tolerance)) {
    if (handle_.empty() || stream_checkpoint_identity_.empty() ||
        line_handle_.empty() || checkpoint_root_.empty() ||
        identity_.empty() || checkpoint_identity_.empty())
      throw std::invalid_argument(
          "transport-pair stream identities cannot be empty");
    if (domain_ != "rational" && domain_ != "acb")
      throw std::invalid_argument(
          "transport-pair stream requires a numeric session domain");
    if (tail_policy_ != TransportTailPolicy::Stored)
      throw std::invalid_argument(
          "transport-pair tile streaming currently requires stored tail policy");
    if (divergent_cancellation_.has_value() && domain_ != "acb")
      throw std::invalid_argument(
          "bounded divergent cancellation is restricted to Acb transport-pair streams");
    require_transport_pair_compatibility(states_[0], states_[1], domain_);
    for (std::size_t side = 0; side < 2; ++side) {
      expected_tiles_[side] = states_[side]->tile_sources().size();
      if (expected_tiles_[side] == 0)
        throw std::invalid_argument(
            "transport-pair stream state has no retained tiles");
      if (epsilon_.required_complete_max >
          states_[side]->public_required_complete_max())
        throw std::invalid_argument(
            "transport-pair stream epsilon target exceeds a state public target");
      states_[side]->require_contraction_counter_capacity(1);
      row_records_[side].reserve(expected_tiles_[side]);
    }
  }

  // Rolling mode contracts a tile while it is still the current march
  // local.  It therefore owns the immutable plan and common anchor, but not
  // a coefficient-bearing transport state containing every historical
  // source.  Only the two terminal locals are retained when the stream is
  // finished so the published line keeps a compact exact owner closure.
  TransportPairObservableStream(
      std::string handle, std::string stream_checkpoint_identity,
      std::string line_handle, std::string checkpoint_root,
      std::string domain, slong precision_bits,
      std::shared_ptr<StoredTilePlan> plan,
      std::shared_ptr<StoredLocalBase> anchor,
      std::string identity, std::string checkpoint_identity,
      ObservableEpsilonContract epsilon, json::object epsilon_record,
      TransportTailPolicy tail_policy,
      std::optional<BoundedDivergentCancellation> divergent_cancellation,
      std::optional<Magnitude> publication_relative_tolerance)
      : handle_(std::move(handle)),
        stream_checkpoint_identity_(
            std::move(stream_checkpoint_identity)),
        line_handle_(std::move(line_handle)),
        checkpoint_root_(std::move(checkpoint_root)),
        domain_(std::move(domain)), precision_bits_(precision_bits),
        rolling_plan_(std::move(plan)), rolling_anchor_(std::move(anchor)),
        identity_(std::move(identity)),
        checkpoint_identity_(std::move(checkpoint_identity)),
        epsilon_(epsilon), epsilon_record_(std::move(epsilon_record)),
        tail_policy_(tail_policy),
        divergent_cancellation_(std::move(divergent_cancellation)),
        publication_relative_tolerance_(
            std::move(publication_relative_tolerance)),
        rolling_(true) {
    if (handle_.empty() || stream_checkpoint_identity_.empty() ||
        line_handle_.empty() || checkpoint_root_.empty() ||
        identity_.empty() || checkpoint_identity_.empty() ||
        !rolling_plan_ || !rolling_anchor_)
      throw std::invalid_argument(
          "rolling transport-pair stream identities or owners cannot be empty");
    if (domain_ != "rational" && domain_ != "acb")
      throw std::invalid_argument(
          "rolling transport-pair stream requires a numeric session domain");
    if (std::string(rolling_anchor_->scalar_domain()) != domain_)
      throw std::invalid_argument(
          "rolling transport-pair anchor changed coefficient domain");
    if (tail_policy_ != TransportTailPolicy::Stored)
      throw std::invalid_argument(
          "rolling transport-pair streaming currently requires stored tails");
    if (divergent_cancellation_.has_value() && domain_ != "acb")
      throw std::invalid_argument(
          "bounded divergent cancellation is restricted to Acb rolling streams");
    const auto& lower = rolling_plan_->arm("lower").exact;
    const auto& upper = rolling_plan_->arm("upper").exact;
    if (lower.direction != -1 || upper.direction != 1 ||
        lower.from.str() != upper.from.str())
      throw std::invalid_argument(
          "rolling transport-pair plan has incompatible exact arms");
    expected_tiles_ = {lower.tiles.size(), upper.tiles.size()};
    for (std::size_t side = 0; side < 2; ++side) {
      if (expected_tiles_[side] == 0)
        throw std::invalid_argument(
            "rolling transport-pair arm has no exact tiles");
      row_records_[side].reserve(expected_tiles_[side]);
    }
  }

  const std::string& handle() const { return handle_; }
  const std::string& stream_checkpoint_identity() const {
    return stream_checkpoint_identity_;
  }
  const std::string& identity() const { return identity_; }
  const std::string& checkpoint_identity() const {
    return checkpoint_identity_;
  }
  const std::string& line_handle() const { return line_handle_; }
  const std::string& checkpoint_root() const { return checkpoint_root_; }
  const auto& states() const { return states_; }
  bool rolling() const { return rolling_; }
  const auto& rolling_plan() const { return rolling_plan_; }
  const auto& rolling_terminals() const { return rolling_terminals_; }
  const auto& expected_tiles() const { return expected_tiles_; }
  const auto& next_tiles() const { return next_tiles_; }
  const auto& epsilon() const { return epsilon_; }
  const auto& publication_relative_tolerance() const {
    return publication_relative_tolerance_;
  }

  json::object add_tile(std::size_t side, std::size_t tile_index,
                        const json::object& prepared_row) {
    if (rolling_)
      throw std::invalid_argument(
          "rolling transport-pair stream requires an explicit current local");
    const std::string arm_name = side == 0 ? "lower" : "upper";
    if (side > 1)
      throw std::invalid_argument(
          "transport-pair stream tile side is invalid");
    const auto source = states_[side]->tile_sources().at(tile_index);
    const auto terminal_match =
        tile_index + 1 == expected_tiles_[side]
        ? states_[side]->terminal_factorized_match()
        : std::shared_ptr<StoredPlannedMatchHop>{};
    return add_tile_from_source(side, tile_index, prepared_row, source,
                                terminal_match, arm_name,
                                states_[side]->plan_owner());
  }

  json::object add_rolling_tile(
      std::size_t side, std::size_t tile_index,
      const std::shared_ptr<StoredLocalBase>& source,
      const json::object& prepared_row) {
    if (!rolling_)
      throw std::invalid_argument(
          "state-backed transport-pair stream does not accept rolling locals");
    if (side > 1 || !source ||
        std::string(source->scalar_domain()) != domain_)
      throw std::invalid_argument(
          "rolling transport-pair tile lost its numeric source local");
    const std::string arm_name = side == 0 ? "lower" : "upper";
    std::shared_ptr<StoredPlannedMatchHop> terminal_match;
    if (tile_index + 1 == expected_tiles_[side]) {
      if (const auto erased = source->terminal_factorized_owner();
          erased != nullptr)
        terminal_match =
            std::static_pointer_cast<StoredPlannedMatchHop>(erased);
    }
    auto response = add_tile_from_source(
        side, tile_index, prepared_row, source, terminal_match, arm_name,
        rolling_plan_);
    if (tile_index + 1 == expected_tiles_[side])
      rolling_terminals_[side] = source;
    return response;
  }

 private:
  json::object add_tile_from_source(
      std::size_t side, std::size_t tile_index,
      const json::object& prepared_row,
      const std::shared_ptr<StoredLocalBase>& source,
      const std::shared_ptr<StoredPlannedMatchHop>& terminal_match,
      const std::string& arm_name,
      const std::shared_ptr<StoredTilePlan>& plan) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (status_ == Status::Poisoned)
      throw std::invalid_argument(
          "transport-pair stream is poisoned and must be aborted");
    if (status_ != Status::Active)
      throw std::invalid_argument(
          "transport-pair stream no longer accepts tiles");
    const std::size_t expected_side =
        next_tiles_[0] < expected_tiles_[0] ? 0 : 1;
    if (side > 1 || side != expected_side)
      throw std::invalid_argument(
          "transport-pair stream requires all lower tiles before upper tiles");
    if (next_tiles_[side] >= expected_tiles_[side] ||
        tile_index != next_tiles_[side])
      throw std::invalid_argument(
          "transport-pair stream tile index is missing, duplicated, or out of order");

    try {
      const auto started = std::chrono::steady_clock::now();
      StoredLineIntegral tile;
      std::unique_ptr<AcbPrecisionLease> acb_lease;
      if (!plan || !source)
        throw std::invalid_argument(
            "transport-pair stream tile lost its plan or source owner");
      const auto& arm = plan->arm(arm_name);
      if (domain_ == "acb") {
        acb_lease = std::make_unique<AcbPrecisionLease>(precision_bits_);
        ComplexBall::set_precision(precision_bits_);
        if (tile_index + 1 == expected_tiles_[side] &&
            terminal_match != nullptr) {
          tile =
              integrate_transport_terminal_factorized_acb_row_tile(
                  precision_bits_, terminal_match, expected_tiles_[side],
                  prepared_row,
                  arm, arm_name, tile_index, epsilon_,
                  divergent_cancellation_);
        } else {
          const auto typed_source =
              std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(
                  source);
          if (!typed_source)
            throw std::logic_error(
                "transport-pair stream Acb source changed coefficient domain");
          tile = integrate_transport_stored_row_tile<ComplexBall>(
              precision_bits_, typed_source, prepared_row, arm,
              arm_name, tile_index, epsilon_,
              divergent_cancellation_);
        }
      } else {
        const auto typed_source =
            std::dynamic_pointer_cast<StoredLocal<Rational>>(
                source);
        if (!typed_source)
          throw std::logic_error(
              "transport-pair stream Rational source changed coefficient domain");
        tile = integrate_transport_stored_row_tile<Rational>(
            precision_bits_, typed_source, prepared_row,
            arm, arm_name, tile_index, epsilon_,
            divergent_cancellation_);
      }
      accumulators_[side].add(tile);
      json::object tile_value{
          {"tile", tile_index},
          {"value", rolling_
               ? encode_compact_transport_line_value_diagnostics(tile)
               : encode_transport_line_value_diagnostics(tile)}};
      const auto row_identity = required_string(
          prepared_row, "exact_identity");
      auto row_record = compact_prepared_row_provenance_record(
          tile_index, prepared_row);
      retained_tile_diagnostic_bytes_[side] +=
          json::serialize(tile_value).size();
      retained_row_record_bytes_[side] +=
          json::serialize(row_record).size();
      tile_values_[side].push_back(std::move(tile_value));
      row_records_[side].push_back(std::move(row_record));
      ++next_tiles_[side];
      const auto elapsed_ms = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      arm_integration_ms_[side] += elapsed_ms;
      return json::object{
          {"side", side == 0 ? "lower" : "upper"},
          {"tile", tile_index}, {"row_exact_identity", row_identity},
          {"lower_complete", next_tiles_[0]},
          {"upper_complete", next_tiles_[1]},
          {"next_side", next_tiles_[0] < expected_tiles_[0]
               ? json::value("lower")
               : next_tiles_[1] < expected_tiles_[1]
                   ? json::value("upper") : json::value(nullptr)},
          {"next_tile", next_tiles_[0] < expected_tiles_[0]
               ? json::value(next_tiles_[0])
               : next_tiles_[1] < expected_tiles_[1]
                   ? json::value(next_tiles_[1]) : json::value(nullptr)},
          {"retained_row_record_bytes",
           retained_row_record_bytes_[0] +
               retained_row_record_bytes_[1]},
          {"retained_tile_diagnostic_bytes",
           retained_tile_diagnostic_bytes_[0] +
               retained_tile_diagnostic_bytes_[1]},
          {"elapsed_ms", elapsed_ms}};
    } catch (...) {
      status_ = Status::Poisoned;
      throw;
    }
  }

 public:

  FinishResult finish() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (status_ == Status::Poisoned)
      throw std::invalid_argument(
          "transport-pair stream is poisoned and cannot finish");
    if (status_ != Status::Active)
      throw std::invalid_argument(
          "transport-pair stream is not active");
    if (next_tiles_ != expected_tiles_)
      throw std::invalid_argument(
          "transport-pair stream cannot finish before every ordered tile was added exactly once");
    status_ = Status::Finishing;
    try {
    const auto finish_started = std::chrono::steady_clock::now();
    const auto observable_root = checkpoint_root_ + ":observable:1";
    std::array<std::shared_ptr<StoredLineResult>, 2> arm_lines;
    if (rolling_) {
      if (!rolling_plan_ || !rolling_terminals_[0] ||
          !rolling_terminals_[1])
        throw std::logic_error(
            "rolling transport-pair stream lost a terminal compact owner");
      for (std::size_t side = 0; side < 2; ++side) {
        const std::string arm_name = side == 0 ? "lower" : "upper";
        const auto& exact = rolling_plan_->arm(arm_name).exact;
        const auto arm_checkpoint = observable_root + ":" + arm_name +
                                    ":rolling-scratch";
        json::object arm_record{
            {"kind", "rolling-transport-observable-arm"},
            {"combination", "sum-physical-tiles"},
            {"request_index", 0},
            {"observable_identity", identity_ + ":" + arm_name},
            {"observable_checkpoint_identity", arm_checkpoint},
            {"output_epsilon_contract", epsilon_record_},
            {"tail_policy", transport_tail_policy_name(tail_policy_)},
            {"projection_mode",
             transport_projection_mode_name(tail_policy_)},
            {"divergent_cancellation",
             divergent_cancellation_.has_value()
                 ? json::value(encode_bounded_divergent_cancellation(
                       *divergent_cancellation_))
                 : json::value(
                       json::object{{"mode", "exact-singleton"}})},
            {"rows", row_records_[side]},
            {"tile_count", expected_tiles_[side]},
            {"coefficient_retention", "terminal-local-only"}};
        arm_lines[side] =
            build_rolling_transport_observable_line_from_result(
                "private:" + arm_checkpoint + ":aggregate",
                arm_checkpoint, arm_name,
                json::object{{"from_exact", exact.from.str()},
                             {"to_exact", exact.to.str()}},
                std::move(arm_record), rolling_plan_,
                rolling_terminals_[side],
                accumulators_[side].finish(
                    "rolling native transport observable aggregate"),
                arm_integration_ms_[side]);
      }
      const auto& lower_exact = rolling_plan_->arm("lower").exact;
      const auto& upper_exact = rolling_plan_->arm("upper").exact;
      json::object pair_record{
          {"kind", "rolling-transport-observable-pair"},
          {"combination", "negative-lower-plus-upper"},
          {"no_remarching", true}, {"no_rematching", true},
          {"request_index", 0},
          {"observable_identity", identity_},
          {"observable_checkpoint_identity", checkpoint_identity_},
          {"output_epsilon_contract", epsilon_record_},
          {"tail_policy", transport_tail_policy_name(tail_policy_)},
          {"lower_rows", row_records_[0]},
          {"upper_rows", row_records_[1]},
          {"lower_tiles", expected_tiles_[0]},
          {"upper_tiles", expected_tiles_[1]},
          {"coefficient_retention", "two-terminal-locals-only"}};
      auto line = build_retained_line_aggregate(
          line_handle_, checkpoint_identity_, "combined",
          json::object{{"from_exact", lower_exact.to.str()},
                       {"anchor_exact", lower_exact.from.str()},
                       {"to_exact", upper_exact.to.str()}},
          std::move(pair_record), rolling_plan_,
          {rolling_terminals_[0], rolling_terminals_[1]},
          {arm_lines[0], arm_lines[1]}, {-1, 1},
          elapsed_milliseconds(finish_started), true);
      if (line->result().value.epsilon.complete_max <
          epsilon_.required_complete_max)
        throw std::domain_error(
            "rolling transport-pair aggregate does not cover its required epsilon maximum");
      FinishResult output;
      output.line = std::move(line);
      output.tiles = expected_tiles_;
      output.arm_integration_ms = arm_integration_ms_;
      output.conditioning = encode_transport_pair_conditioning_diagnostics(
          arm_lines[0]->result(), arm_lines[1]->result(),
          output.line->result(), std::move(tile_values_[0]),
          std::move(tile_values_[1]));
      output.elapsed_ms = arm_integration_ms_[0] +
          arm_integration_ms_[1] + elapsed_milliseconds(finish_started);
      status_ = Status::Finished;
      return output;
    }
    for (std::size_t side = 0; side < 2; ++side) {
      const std::string arm_name = side == 0 ? "lower" : "upper";
      const auto arm_identity = identity_ + ":" + arm_name;
      const auto arm_checkpoint =
          observable_root + ":" + arm_name + ":scratch";
      json::object arm_record{
          {"kind", "transport-state-observable-arm"},
          {"combination", "sum-physical-tiles"},
          {"request_index", 0},
          {"observable_identity", arm_identity},
          {"observable_checkpoint_identity", arm_checkpoint},
          {"transport_state",
           compact_transport_state_reference(states_[side])},
          {"output_epsilon_contract", epsilon_record_},
          {"tail_policy", transport_tail_policy_name(tail_policy_)},
          {"projection_mode",
           transport_projection_mode_name(tail_policy_)},
          {"divergent_cancellation",
           divergent_cancellation_.has_value()
               ? json::value(encode_bounded_divergent_cancellation(
                     *divergent_cancellation_))
               : json::value(json::object{{"mode", "exact-singleton"}})},
          {"rows", row_records_[side]},
          {"tile_count", expected_tiles_[side]}};
      const auto& exact = states_[side]->plan_owner()->arm(arm_name).exact;
      arm_lines[side] = build_compact_transport_observable_line_from_result(
          "private:" + observable_root + ":" + arm_name + ":aggregate",
          arm_checkpoint, arm_name,
          json::object{{"from_exact", exact.from.str()},
                       {"to_exact", exact.to.str()}},
          std::move(arm_record), states_[side],
          accumulators_[side].finish(
              "compact streamed native transport observable aggregate"),
          expected_tiles_[side], arm_integration_ms_[side]);
    }

    json::object pair_record{
        {"kind", "transport-state-observable-pair"},
        {"combination", "negative-lower-plus-upper"},
        {"no_remarching", true}, {"no_rematching", true},
        {"concurrent_arms", true}, {"request_index", 0},
        {"observable_identity", identity_},
        {"observable_checkpoint_identity", checkpoint_identity_},
        {"output_epsilon_contract", epsilon_record_},
        {"tail_policy", transport_tail_policy_name(tail_policy_)},
        {"projection_mode", transport_projection_mode_name(tail_policy_)},
        {"divergent_cancellation",
         divergent_cancellation_.has_value()
             ? json::value(encode_bounded_divergent_cancellation(
                   *divergent_cancellation_))
             : json::value(json::object{{"mode", "exact-singleton"}})},
        {"lower", json::object{
             {"transport_state",
              compact_transport_state_reference(states_[0])},
             {"rows", row_records_[0]},
             {"tile_count", expected_tiles_[0]}}},
        {"upper", json::object{
             {"transport_state",
              compact_transport_state_reference(states_[1])},
             {"rows", row_records_[1]},
             {"tile_count", expected_tiles_[1]}}}};
    auto line = build_compact_transport_pair_observable_line(
        line_handle_, checkpoint_identity_, std::move(pair_record),
        states_[0], states_[1], arm_lines[0], arm_lines[1]);
    if (line->result().value.epsilon.complete_max <
        epsilon_.required_complete_max)
      throw std::domain_error(
          "streamed transport-pair aggregate does not cover its required epsilon maximum");
    FinishResult output;
    output.line = std::move(line);
    output.tiles = expected_tiles_;
    output.arm_integration_ms = arm_integration_ms_;
    output.conditioning = encode_transport_pair_conditioning_diagnostics(
        arm_lines[0]->result(), arm_lines[1]->result(),
        output.line->result(), std::move(tile_values_[0]),
        std::move(tile_values_[1]));
    output.elapsed_ms = arm_integration_ms_[0] + arm_integration_ms_[1] +
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - finish_started).count();
    status_ = Status::Finished;
    return output;
    } catch (...) {
      status_ = Status::Poisoned;
      throw;
    }
  }

  void abort() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (status_ == Status::Finished)
      throw std::invalid_argument(
          "finished transport-pair stream cannot be aborted");
    status_ = Status::Aborted;
  }

 private:
  std::string handle_;
  std::string stream_checkpoint_identity_;
  std::string line_handle_;
  std::string checkpoint_root_;
  std::string domain_;
  slong precision_bits_ = 256;
  std::array<std::shared_ptr<StoredTransportArmState>, 2> states_;
  std::shared_ptr<StoredTilePlan> rolling_plan_;
  std::shared_ptr<StoredLocalBase> rolling_anchor_;
  std::array<std::shared_ptr<StoredLocalBase>, 2> rolling_terminals_;
  std::string identity_;
  std::string checkpoint_identity_;
  ObservableEpsilonContract epsilon_;
  json::object epsilon_record_;
  TransportTailPolicy tail_policy_ = TransportTailPolicy::Stored;
  std::optional<BoundedDivergentCancellation> divergent_cancellation_;
  std::optional<Magnitude> publication_relative_tolerance_;
  std::array<std::size_t, 2> expected_tiles_{0, 0};
  std::array<std::size_t, 2> next_tiles_{0, 0};
  std::array<StreamingStoredLineAccumulator, 2> accumulators_;
  std::array<json::array, 2> tile_values_;
  std::array<json::array, 2> row_records_;
  std::array<std::size_t, 2> retained_tile_diagnostic_bytes_{0, 0};
  std::array<std::size_t, 2> retained_row_record_bytes_{0, 0};
  std::array<double, 2> arm_integration_ms_{0.0, 0.0};
  bool rolling_ = false;
  Status status_ = Status::Active;
  std::mutex mutex_;
};

std::vector<TransportObservableContractionResult> contract_transport_arm(
    const std::string& domain, slong precision_bits,
    const std::shared_ptr<StoredTilePlan>& plan, const std::string& arm_name,
    const std::vector<std::shared_ptr<StoredLocalBase>>& tile_sources,
    const std::vector<TransportObservableContractionInput>& observables,
    const std::shared_ptr<StoredTransportArmState>& transport_owner =
        nullptr) {
  if ((domain != "rational" && domain != "acb") || !plan)
    throw std::invalid_argument(
        "native transport contraction requires a numeric domain and retained plan");
  const auto& retained = plan->arm(arm_name);
  if (tile_sources.size() != retained.exact.tiles.size() ||
      tile_sources.empty())
    throw std::invalid_argument(
        "native transport contraction sources do not reproduce the retained arm tiles");
  if (domain == "acb") ComplexBall::set_precision(precision_bits);

  std::vector<TransportObservableContractionResult> results;
  results.reserve(observables.size());
  for (std::size_t observable_index = 0;
       observable_index < observables.size(); ++observable_index) {
    const auto& observable = observables[observable_index];
    if (observable.identity.empty() || observable.checkpoint_identity.empty() ||
        observable.checkpoint_root.empty() ||
        observable.aggregate_handle.empty() ||
        observable.row_count() != tile_sources.size() ||
        observable.projected_local_handles.size() != tile_sources.size())
      throw std::invalid_argument(
          "native observable contraction does not reproduce its retained tile topology or identities");
    const auto started = std::chrono::steady_clock::now();
    TransportObservableContractionResult output;
    output.identity = observable.identity;
    StreamingStoredLineAccumulator accumulator;
    for (std::size_t tile = 0; tile < tile_sources.size(); ++tile) {
      const auto& source = tile_sources[tile];
      const auto& prepared_row = observable.row(tile);
      const auto row_identity = required_string(
          prepared_row, "exact_identity");
      (void)row_identity;
      const auto tile_started = std::chrono::steady_clock::now();
      if (transport_contraction_timing_enabled()) {
        std::ostringstream timing;
        timing << "transport-contraction-timing"
               << " phase=tile-start"
               << " arm=" << arm_name
               << " observable_index=" << observable_index
               << " observable_identity=" << observable.identity
               << " tile=" << tile
               << " terminal="
               << ((domain == "acb" &&
                    tile + 1 == tile_sources.size() &&
                    transport_owner &&
                    transport_owner->terminal_factorized_match() != nullptr)
                       ? "true"
                       : "false");
        emit_transport_contraction_timing(timing.str());
      }
      if (transport_owner &&
          observable.tail_policy == TransportTailPolicy::Stored) {
        const auto fused_started = std::chrono::steady_clock::now();
        StoredLineIntegral fused;
        if (domain == "rational") {
          const auto typed =
              std::dynamic_pointer_cast<StoredLocal<Rational>>(source);
          if (!typed)
            throw std::logic_error(
                "fused transport Rational source changed coefficient domain");
          fused = integrate_transport_stored_row_tile<Rational>(
              precision_bits, typed, prepared_row, retained, arm_name, tile,
              observable.epsilon, observable.divergent_cancellation);
        } else {
          if (tile + 1 == tile_sources.size() &&
              transport_owner->terminal_factorized_match() !=
                  nullptr) {
            try {
              fused =
                  integrate_transport_terminal_factorized_acb_row_tile(
                      precision_bits, transport_owner, prepared_row,
                      retained, arm_name, tile, observable.epsilon,
                      observable.divergent_cancellation);
            } catch (const MatchingArithmeticError& error) {
              if (error.code != MatchingArithmeticErrorCode::
                                    TerminalOutputInconclusive)
                throw;
              throw MatchingArithmeticError(
                  error.code, error.what(), error.row, observable_index,
                  error.epsilon_power);
            }
          } else {
            const auto typed =
                std::dynamic_pointer_cast<
                    StoredLocal<ComplexBall>>(source);
            if (!typed)
              throw std::logic_error(
                  "fused transport Acb source changed coefficient domain");
            fused =
                integrate_transport_stored_row_tile<ComplexBall>(
                    precision_bits, typed, prepared_row, retained,
                    arm_name, tile, observable.epsilon,
                    observable.divergent_cancellation,
                    observable.factorize_ordinary_stored_rows);
          }
        }
        output.tile_integration_ms +=
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - fused_started)
                .count();
        output.tile_values.push_back(json::object{
            {"tile", tile},
            {"value", encode_transport_line_value_diagnostics(fused)}});
        accumulator.add(fused);
        ++output.tile_integrations;
        if (transport_contraction_timing_enabled()) {
          std::ostringstream timing;
          timing << "transport-contraction-timing"
                 << " phase=tile-done"
                 << " arm=" << arm_name
                 << " observable_index=" << observable_index
                 << " observable_identity=" << observable.identity
                 << " tile=" << tile
                 << " elapsed_ms="
                 << elapsed_milliseconds(tile_started);
          emit_transport_contraction_timing(timing.str());
        }
        continue;
      }
      json::object row_request{
          {"source_checkpoint_identity", source->checkpoint_identity()},
          {"checkpoint_identity",
           arm_checkpoint_identity(observable.checkpoint_root, arm_name,
                                   "integrand", tile + 1) +
               ":" + row_identity}};
      const auto projection_complete_max =
          local_algebra_detail::checked_i32(
              static_cast<std::int64_t>(
                  observable.epsilon.requested.complete_max) + 1,
              "transport scalar-row primitive halo");
      std::shared_ptr<StoredLocalBase> projected;
      if (domain == "rational") {
        const auto typed =
            std::dynamic_pointer_cast<StoredLocal<Rational>>(source);
        if (!typed)
          throw std::logic_error(
              "transport contraction Rational source changed coefficient domain");
        projected = build_rational_row_local<Rational>(
            observable.projected_local_handles[tile], row_request,
            precision_bits, typed, source,
            observable.tail_policy != TransportTailPolicy::Stored,
            projection_complete_max, &prepared_row);
      } else {
        const auto typed =
            std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(source);
        if (!typed)
          throw std::logic_error(
              "transport contraction Acb source changed coefficient domain");
        projected = build_rational_row_local<ComplexBall>(
            observable.projected_local_handles[tile], row_request,
            precision_bits, typed, source,
            observable.tail_policy != TransportTailPolicy::Stored,
            projection_complete_max, &prepared_row);
      }
      const auto line_epsilon = live_line_epsilon_intersection(
          observable.epsilon.requested,
          observable.epsilon.required_complete_max, projected,
          1 /* exact native primitive bound: min epsilon power >= -1 */);
      json::object line_request{
          {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
          {"source_checkpoint_identity", projected->checkpoint_identity()},
          {"checkpoint_identity",
           arm_checkpoint_identity(observable.checkpoint_root, arm_name,
                                   "tile", tile + 1)},
          {"arm", arm_name},
          {"tile", tile},
          {"epsilon", json::object{{"min", line_epsilon.min_power},
                                    {"max", line_epsilon.complete_max}}}};
      if (observable.tail_policy != TransportTailPolicy::Stored)
        line_request["certify_tail"] = true;
      if (observable.divergent_cancellation.has_value())
        line_request["divergent_cancellation"] =
            encode_bounded_divergent_cancellation(
                *observable.divergent_cancellation);
      auto tile_line = build_planned_line_result(
          "private:" + arm_checkpoint_identity(
                           observable.checkpoint_root, arm_name, "tile",
                           tile + 1),
          line_request, plan, projected);
      output.tile_integration_ms += tile_line->elapsed_ms();
      output.tile_values.push_back(json::object{
          {"tile", tile},
          {"value",
           encode_transport_line_value_diagnostics(tile_line->result())}});
      if (transport_owner) {
        accumulator.add(tile_line->result());
      } else {
        output.projected.push_back(projected);
        output.tile_lines.push_back(tile_line);
      }
      ++output.tile_integrations;
      if (transport_contraction_timing_enabled()) {
        std::ostringstream timing;
        timing << "transport-contraction-timing"
               << " phase=tile-done"
               << " arm=" << arm_name
               << " observable_index=" << observable_index
               << " observable_identity=" << observable.identity
               << " tile=" << tile
               << " elapsed_ms="
               << elapsed_milliseconds(tile_started);
        emit_transport_contraction_timing(timing.str());
      }
      // The line owns the projected local and the projected local owns its
      // source.  Drop both before advancing so contraction peak memory is
      // independent of the number of tiles.
      if (transport_owner) {
        tile_line.reset();
        projected.reset();
      }
    }
    const auto aggregate_started = std::chrono::steady_clock::now();
    const auto interval = json::object{
        {"from_exact", retained.exact.from.str()},
        {"to_exact", retained.exact.to.str()}};
    if (transport_owner) {
      output.aggregate = build_compact_transport_observable_line_from_result(
          observable.aggregate_handle, observable.checkpoint_identity,
          arm_name, interval, observable.aggregate_record, transport_owner,
          accumulator.finish("compact native transport observable aggregate"),
          output.tile_integrations,
          std::chrono::duration<double, std::milli>(
              std::chrono::steady_clock::now() - aggregate_started)
              .count());
    } else {
      output.aggregate = build_retained_line_aggregate(
          observable.aggregate_handle, observable.checkpoint_identity,
          arm_name, interval, observable.aggregate_record, plan,
          output.projected, output.tile_lines,
          std::vector<std::int32_t>(output.tile_lines.size(), 1),
          std::chrono::duration<double, std::milli>(
              std::chrono::steady_clock::now() - aggregate_started)
              .count());
    }
    if (output.aggregate->result().value.epsilon.complete_max <
        observable.epsilon.required_complete_max)
      throw std::domain_error(
          "native observable aggregate does not cover its required epsilon maximum");
    if (observable.tail_policy == TransportTailPolicy::Require &&
        output.aggregate->result().scope !=
            LineIntegrationScope::FullLocalWithCertifiedTail)
      throw std::domain_error(
          "native observable contraction requires a certified full-local tail but certification was unavailable");
    output.elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    results.push_back(std::move(output));
  }
  return results;
}
