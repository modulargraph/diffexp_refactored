  if (operation == "transport.consume_physical_value_hop") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "tile_plan",
         "tile_plan_checkpoint_identity", "arm", "match", "value_solver",
         "incoming", "incoming_checkpoint_identity", "epsilon",
         "checkpoint_policy"},
        "native transport.consume_physical_value_hop request");
    if (session->domain != "acb")
      throw std::invalid_argument(
          "transport.consume_physical_value_hop requires Acb coefficients");
    const auto arm_name = required_string(root, "arm");
    const auto match_index = checkpoint_size_t(
        root.at("match"), "physical value transport hop match index");
    const auto incoming_handle = required_string(root, "incoming");
    const auto& value_solver = as_object(
        root.at("value_solver"), "physical regular value-solver prototype");
    require_exact_keys(value_solver,
        {"schema", "taylor_complete_max", "metadata",
         "relative_accuracy_max_exact"},
        "physical regular value-solver prototype");
    if (required_string(value_solver, "schema") !=
            "diffexp2-native-ordinary-physical-value-solver-v1")
      throw std::invalid_argument(
          "unsupported physical regular value-solver prototype schema");
    const auto receiver_taylor_complete_max = as_u32(
        value_solver.at("taylor_complete_max"),
        "physical value-solver Taylor maximum");
    const auto& metadata_prototype = as_object(
        value_solver.at("metadata"),
        "physical regular value-solver metadata prototype");
    const auto relative_accuracy_text = required_string(
        value_solver, "relative_accuracy_max_exact");
    const Rational relative_accuracy_max(relative_accuracy_text);
    if (relative_accuracy_max.sign() <= 0 ||
        !(relative_accuracy_max < Rational(1)) ||
        relative_accuracy_max.str() != relative_accuracy_text)
      throw std::invalid_argument(
          "physical value-solver relative-accuracy threshold must be a canonical exact rational strictly between zero and one");

    const auto& epsilon = as_object(
        root.at("epsilon"), "physical value transport hop epsilon contract");
    require_exact_keys(epsilon, {"min", "max", "required_complete_max"},
                       "physical value transport hop epsilon contract");
    EpsilonWindow requested_epsilon{
        as_i32(epsilon.at("min"), "physical value hop epsilon minimum"),
        as_i32(epsilon.at("max"), "physical value hop epsilon maximum")};
    (void)requested_epsilon.width();
    const auto required_complete_max = as_i32(
        epsilon.at("required_complete_max"),
        "physical value hop required epsilon maximum");
    if (required_complete_max < requested_epsilon.min_power ||
        required_complete_max > requested_epsilon.complete_max)
      throw std::invalid_argument(
          "physical value transport hop epsilon contract is inconsistent");
    const auto& checkpoint_policy = as_object(
        root.at("checkpoint_policy"),
        "physical value transport hop checkpoint policy");
    require_exact_keys(checkpoint_policy, {"schema", "root"},
                       "physical value transport hop checkpoint policy");
    if (required_string(checkpoint_policy, "schema") !=
            "diffexp2-deterministic-arm-checkpoints-v1")
      throw std::invalid_argument(
          "unsupported physical value transport hop checkpoint policy");
    const auto checkpoint_root = required_string(checkpoint_policy, "root");

    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> incoming;
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
            "physical value transport hop tile-plan binding is stale");
      plan = plan_found->second;
      const auto& retained = plan->arm(arm_name);
      if (match_index >= retained.exact.matches.size())
        throw std::invalid_argument(
            "physical value transport hop match lies outside its exact arm");
      const auto incoming_found = session->locals.find(incoming_handle);
      if (incoming_found == session->locals.end() ||
          required_string(root, "incoming_checkpoint_identity") !=
              incoming_found->second->checkpoint_identity())
        throw std::invalid_argument(
            "physical value transport hop incoming-local binding is stale");
      incoming = incoming_found->second;
      ++session->total_transport_physical_value_hop_attempts;
    }

    const auto ineligible = [&](std::string reason,
                                std::string detail = {}) {
      {
        std::lock_guard<std::mutex> lock(session->mutex);
        ++session->total_transport_physical_value_hop_ineligible;
      }
      json::object response{
          {"status", "ok"}, {"session", session->handle},
          {"capability", "ordinary-physical-value-transport-hop-v1"},
          {"execution_mode", "framed-fallback-required"},
          {"used", false}, {"reason", std::move(reason)},
          {"arm", arm_name}, {"match", match_index},
          {"value_hops", 0}, {"basis_matches", 0}};
      if (!detail.empty())
        response["detail"] = std::move(detail);
      return response;
    };
    const auto& retained = plan->arm(arm_name);
    const auto& exact_match = retained.exact.matches[match_index];
    const auto& producing = retained.charts.at(exact_match.producing_chart);
    const auto& receiving = retained.charts.at(exact_match.receiving_chart);
    if (producing.geometry.singular_center ||
        receiving.geometry.singular_center)
      return ineligible("singular-chart-crossing");
    const auto producing_equation_owner = std::visit(
        [](const auto& owner) ->
            std::shared_ptr<PhysicalEquationOwnerBase> {
          using Owner = typename std::decay_t<decltype(owner)>::element_type;
          if constexpr (std::is_same_v<Owner, PreparedChartBase> ||
                        std::is_same_v<Owner,
                                       RegularPhysicalEquationOwnerBase>)
            return owner;
          return nullptr;
        }, producing.owner);
    const auto* receiving_owner = std::get_if<
        std::shared_ptr<RegularPhysicalEquationOwnerBase>>(&receiving.owner);
    if (!producing_equation_owner)
      return ineligible(
          "producing-owner-is-not-a-regular-physical-equation");
    if (!receiving_owner || !*receiving_owner)
      return ineligible(
          "receiving-owner-is-not-a-frame-independent-regular-equation");
    if ((*receiving_owner)->regular_value_relative_accuracy_max_exact() !=
        relative_accuracy_text)
      throw std::invalid_argument(
          "physical value-solver relative-accuracy threshold differs from its equation owner");

    std::optional<Rational> producing_center;
    std::optional<Rational> producing_scale;
    std::optional<Rational> receiving_center;
    try {
      producing_center = Rational(producing.local_geometry.center_exact);
      producing_scale = Rational(producing.local_geometry.scale_exact);
      receiving_center = Rational(receiving.local_geometry.center_exact);
    } catch (const std::invalid_argument&) {
      return ineligible(
          "certified-algebraic-chart-requires-basis-match");
    }
    if (!(*producing_center == producing.geometry.center) ||
        !(*producing_scale == producing.geometry.scale) ||
        !(*receiving_center == receiving.geometry.center) ||
        !acb_equal(producing.local_geometry.radius.raw(),
                   ComplexBall::from_strings(
                       producing.geometry.radius.str()).raw()) ||
        !acb_equal(receiving.local_geometry.radius.raw(),
                   ComplexBall::from_strings(
                       receiving.geometry.radius.str()).raw()))
      return ineligible(
          "certified-algebraic-chart-requires-basis-match");
    if (producing_scale->is_zero())
      throw std::invalid_argument(
          "physical value transport hop producing chart has a zero exact scale");
    const auto producing_local =
        (*receiving_center - *producing_center) / *producing_scale;
    const auto center_ratio = exact_path_detail::abs(producing_local) /
                              producing.geometry.radius;
    if (!(center_ratio < Rational(1)))
      return ineligible(
          "receiver-center-lies-outside-producing-certified-disk");

    const auto incoming_summary = incoming->summary();
    const auto incoming_epsilon_min = as_i32(
        incoming_summary.at("epsilon_min"),
        "physical value-hop incoming epsilon minimum");
    const auto incoming_epsilon_max = as_i32(
        incoming_summary.at("epsilon_max"),
        "physical value-hop incoming epsilon maximum");
    const auto incoming_top_valid = parse_validity(
        incoming_summary.at("top_valid"));
    const auto producer_taylor_complete_max = as_u32(
        incoming_summary.at("taylor_complete_max"),
        "physical value-hop producer Taylor maximum");
    if (incoming->source_chart() != producing.handle ||
        incoming->source_operator_identity() != producing.exact_identity)
      throw std::invalid_argument(
          "physical value transport hop incoming local differs from its producing plan owner");
    incoming->require_exact_plan_binding(
        producing.local_geometry, producing.prescriptions,
        "physical value transport hop incoming local");
    const auto incoming_equation_owner = incoming->retained_equation_owner();
    if (!incoming_equation_owner ||
        incoming_equation_owner.get() != producing_equation_owner.get())
      return ineligible(
          "incoming-equation-owner-is-not-the-producing-plan-owner");
    if (incoming_top_valid < incoming_epsilon_max)
      return ineligible(
          "incoming-local-full-epsilon-window-is-not-complete");
    if (incoming_epsilon_max < requested_epsilon.complete_max)
      return ineligible(
          "incoming-local-lacks-complete-epsilon-coverage");
    // EpsilonWindow certifies every row below its stored minimum as
    // structural zero.  Local assembly may therefore trim a requested lower
    // row without losing information.  Preserve any extra stored lower rows,
    // but widen a later stored minimum down to the requested edge with exact
    // zeros before causal physical evolution.
    const auto execution_epsilon_min =
        std::min(incoming_epsilon_min, requested_epsilon.min_power);
    const auto incoming_dimension = as_u32(
        incoming_summary.at("dimension"),
        "physical value-hop incoming dimension");
    (void)physical_ode_detail::checked_physical_evolution_coefficient_count(
        EpsilonWindow{execution_epsilon_min, incoming_epsilon_max},
        incoming_dimension, receiver_taylor_complete_max);

    auto acb_incoming =
        std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(incoming);
    if (!acb_incoming)
      return ineligible(
          "incoming-local-is-not-an-acb-local");
    const bool source_has_framed_tail =
        acb_incoming->tail_model().status == TailMajorantStatus::Certified &&
        acb_incoming->tail_model().model.has_value();
    std::optional<PhysicalRegularTaylorTailModel> physical_source_model;
    std::string physical_source_model_failure;
    const auto& source_equation = acb_incoming->physical_equation();
    if (source_equation &&
        source_equation->payload_identity ==
            producing_equation_owner->physical_payload_identity()) {
      auto prepared_physical =
          prepare_physical_regular_homogeneous_tail_model(
              *source_equation, acb_incoming->solution());
      if (prepared_physical.status == TailMajorantStatus::Certified &&
          prepared_physical.model.has_value())
        physical_source_model = std::move(*prepared_physical.model);
      else
        physical_source_model_failure = std::move(prepared_physical.detail);
    } else {
      physical_source_model_failure =
          "incoming local has no matching retained physical q/C owner";
    }
    if (!physical_source_model.has_value() && !source_has_framed_tail)
      return ineligible(
          "incoming-local-has-no-certified-physical-tail-model",
          physical_source_model_failure);
    const bool use_physical_source_tail =
        physical_source_model.has_value();
    auto equation = std::dynamic_pointer_cast<
        RegularPhysicalEquationOwner<ComplexBall>>(*receiving_owner);
    if (!equation)
      throw std::invalid_argument(
          "physical value transport hop receiver has the wrong Acb equation type");
    std::uint32_t tail_certificate_taylor_complete_max =
        use_physical_source_tail
        ? physical_source_model->taylor_complete_max
        : producer_taylor_complete_max;
    constexpr std::uint32_t kWitnessSearchCap = 16;
    // Reconstructing one physical value is substantially cheaper than
    // solving every basis column.  A coefficient-blind center ratio cannot
    // predict the required Taylor order for a stiff local operator, so retry
    // the exact retained q/C recurrence on a deterministic private-order
    // ladder.  The published local keeps receiver_taylor_complete_max; only
    // this handoff certificate sees the extended prefix.
    constexpr std::uint32_t kPrivateTaylorOrderCap = 800;
    constexpr std::uint32_t kPrivateTaylorOrderMinimumStep = 25;
    std::vector<std::uint32_t> tail_order_attempts{
        tail_certificate_taylor_complete_max};
    std::vector<std::string> tail_order_accuracy_trace;
    std::uint32_t current_accuracy_deficit =
        std::numeric_limits<std::uint32_t>::max();
    std::optional<Rational> current_accuracy_excess;
    std::optional<Rational> witness_radius;
    std::uint32_t witness_dyadic_exponent = 0;
    std::string witness_dyadic_direction;
    auto producing_point = RealEvaluationPoint::rational(
        producing_local.str());
    const auto producing_rim = exact_plan_rim(
        producing.prescriptions, *producing_scale);
    auto expected_rim = producing_point.sign < 0
        ? producing_rim : std::nullopt;
    auto point_modulus = exact_path_detail::abs(producing_local);
    auto witness_gap = producing.geometry.radius - point_modulus;
    std::optional<CertifiedLocalEvaluation> certified;
    EpsilonVector retained_value;
    LocalEvaluation certified_handoff;
    std::string last_accuracy_failure;
    bool certificate_failed_after_disk = false;
    bool used_overlap_recenter = false;
    std::string direct_center_failure;
    std::optional<CertifiedLocalEvaluation> recentered_certified;
    EpsilonVector source_overlap_inflated_value;
    EpsilonVector recentered_retained_value;
    LocalEvaluation recentered_center_handoff;
    std::optional<Rational> recentered_witness_radius;
    std::uint32_t recentered_witness_dyadic_exponent = 0;
    std::string recentered_witness_dyadic_direction;
    std::uint32_t recentered_taylor_complete_max = 0;
    std::vector<std::uint32_t> recentered_order_attempts;
    Rational recentered_chart_radius(0);
    const auto try_witness = [&](const Rational& candidate,
                                 std::string direction,
                                 std::uint32_t exponent) {
      std::optional<CertifiedLocalEvaluation> attempted;
      if (use_physical_source_tail) {
        EvaluationOptions options;
        options.imaginary_sign = producing_rim;
        options.compute_tail_estimate = false;
        attempted =
            evaluate_physical_local_solution_with_certified_tail(
                *physical_source_model, producing_point,
                candidate.str(), options);
      } else {
        attempted =
            incoming->evaluate_retained_point_with_certified_tail(
                producing_point, producing_rim, candidate.str());
      }
      if (!attempted.has_value()) return false;
      if (attempted->tail.status != TailMajorantStatus::Certified) {
        if (attempted->tail.status !=
                TailMajorantStatus::Inconclusive ||
            attempted->tail.disk.status !=
                TailMajorantStatus::Inconclusive)
          certificate_failed_after_disk = true;
        return false;
      }
      if (attempted->evaluation.imaginary_sign != expected_rim)
        throw std::logic_error(
            "center-evaluation rim differs from its producing plan chart");
      if (attempted->evaluation.value.epsilon.min_power <
              incoming_epsilon_min ||
          attempted->evaluation.value.epsilon.complete_max !=
              incoming_epsilon_max ||
          attempted->tail.value.guarantee !=
              ErrorGuarantee::Certified ||
          !tail_majorant_detail::same_epsilon_window(
              EpsilonWindow{incoming_epsilon_min,
                            incoming_epsilon_max},
              attempted->tail.value.frame) ||
          attempted->tail.value.absolute.size() !=
              EpsilonWindow{incoming_epsilon_min,
                            incoming_epsilon_max}.width())
        throw std::logic_error(
            "certified physical value handoff changed its retained epsilon frame");
      auto candidate_retained = attempted->evaluation.value;
      candidate_retained.error = ErrorEnvelope{};
      if (incoming_epsilon_min <
          candidate_retained.epsilon.min_power) {
        auto widened = physical_ode_detail::zero_epsilon_vector(
            EpsilonWindow{incoming_epsilon_min,
                          incoming_epsilon_max},
            candidate_retained.dimension);
        for (std::int64_t raw_power =
                 candidate_retained.epsilon.min_power;
             raw_power <= candidate_retained.epsilon.complete_max;
             ++raw_power) {
          const auto power = static_cast<std::int32_t>(raw_power);
          for (std::uint32_t component = 0;
               component < candidate_retained.dimension; ++component)
            widened.at(power, component) =
                candidate_retained.at(power, component);
        }
        candidate_retained = std::move(widened);
      }
      auto candidate_handoff = attempted->evaluation;
      candidate_handoff.value = candidate_retained;
      for (std::int64_t raw_power =
               candidate_handoff.value.epsilon.min_power;
           raw_power <=
               candidate_handoff.value.epsilon.complete_max;
           ++raw_power) {
        const auto power = static_cast<std::int32_t>(raw_power);
        const auto row = static_cast<std::size_t>(
            raw_power -
            candidate_handoff.value.epsilon.min_power);
        for (std::uint32_t component = 0;
             component < candidate_handoff.value.dimension;
             ++component)
          attempted->tail.value.absolute[row].add_error_to(
              candidate_handoff.value.at(power, component));
      }
      candidate_handoff.value.error = ErrorEnvelope{};
      if (execution_epsilon_min <
          candidate_handoff.value.epsilon.min_power) {
        auto widened = physical_ode_detail::zero_epsilon_vector(
            EpsilonWindow{execution_epsilon_min,
                          incoming_epsilon_max},
            candidate_handoff.value.dimension);
        for (std::int64_t raw_power = incoming_epsilon_min;
             raw_power <= incoming_epsilon_max; ++raw_power) {
          const auto power = static_cast<std::int32_t>(raw_power);
          for (std::uint32_t component = 0;
               component < candidate_handoff.value.dimension;
               ++component)
            widened.at(power, component) =
                candidate_handoff.value.at(power, component);
        }
        candidate_handoff.value = std::move(widened);
      }
      const auto accuracy_deficit =
          value_handoff_accuracy_required_additional_digits(
              candidate_handoff.value, required_complete_max,
              relative_accuracy_max);
      const auto accuracy_excess =
          value_handoff_accuracy_excess_ratio(
              candidate_handoff.value, required_complete_max,
              relative_accuracy_max);
      current_accuracy_deficit =
          std::min(current_accuracy_deficit, accuracy_deficit);
      if (accuracy_excess.has_value() &&
          (!current_accuracy_excess.has_value() ||
           *accuracy_excess < *current_accuracy_excess))
        current_accuracy_excess = *accuracy_excess;
      if (accuracy_deficit != 0) {
        const auto failure = value_handoff_accuracy_failure(
            candidate_handoff.value, required_complete_max,
            relative_accuracy_max);
        if (!failure.has_value())
          throw std::logic_error(
              "value handoff deficit and failure diagnostics disagree");
        last_accuracy_failure =
            *failure + ";required_additional_digits=" +
            std::to_string(accuracy_deficit) +
            ";witness_radius_exact=" +
            candidate.str() + ";witness_direction=" + direction +
            ";witness_dyadic_exponent=" +
            std::to_string(exponent);
        return false;
      }
      witness_radius = candidate;
      witness_dyadic_direction = std::move(direction);
      witness_dyadic_exponent = exponent;
      retained_value = std::move(candidate_retained);
      certified_handoff = std::move(candidate_handoff);
      certified = std::move(attempted);
      return true;
    };
    const auto dyadic_denominator = [](std::uint32_t exponent) {
      Rational denominator(1);
      for (std::uint32_t index = 0; index < exponent; ++index)
        denominator *= Rational(2);
      return denominator;
    };
    const auto search_witness = [&](const std::string& phase) {
      witness_radius.reset();
      witness_dyadic_exponent = 0;
      witness_dyadic_direction.clear();
      certified.reset();
      retained_value = EpsilonVector{};
      certified_handoff = LocalEvaluation{};
      last_accuracy_failure.clear();
      certificate_failed_after_disk = false;
      current_accuracy_deficit =
          std::numeric_limits<std::uint32_t>::max();
      current_accuracy_excess.reset();
      for (std::uint32_t exponent = kWitnessSearchCap;
           exponent >= 1 && !witness_radius.has_value(); --exponent) {
        const auto candidate = producing.geometry.radius -
            witness_gap / dyadic_denominator(exponent);
        (void)try_witness(candidate, "outward", exponent);
      }
      for (std::uint32_t exponent = 2;
           exponent <= kWitnessSearchCap &&
               !witness_radius.has_value(); ++exponent) {
        const auto candidate = point_modulus +
            witness_gap / dyadic_denominator(exponent);
        (void)try_witness(candidate, "inward", exponent);
      }
    };
    // A private Taylor extension is useful only while it tightens the exact
    // worst radius/allowance ratio.  Stop at the first exact-ratio plateau:
    // unlike a rounded decimal-digit count, this comparison cannot hide a
    // sub-bit improvement.  Continuing once the exact ratio is unchanged
    // only rebuilds larger Taylor prefixes against a fixed input/arithmetic
    // radius floor.  Direct-center and overlap evaluation have independent
    // geometry, so each receives its own progress history.
    constexpr std::uint32_t kPrivateTaylorNonprogressCap = 1;
    std::optional<Rational> phase_best_accuracy_excess;
    std::uint32_t phase_nonprogress = 0;
    const auto note_phase_progress = [&] {
      std::string status = "unmeasured";
      if (!current_accuracy_excess.has_value()) return status;
      if (!phase_best_accuracy_excess.has_value() ||
          *current_accuracy_excess < *phase_best_accuracy_excess) {
        status = phase_best_accuracy_excess.has_value()
            ? "improved"
            : "initial";
        phase_best_accuracy_excess = *current_accuracy_excess;
        phase_nonprogress = 0;
      } else {
        status = "stalled";
        ++phase_nonprogress;
      }
      return status;
    };
    const auto record_phase_attempt =
        [&](const std::string& phase, const std::string& progress) {
      tail_order_accuracy_trace.push_back(
          phase + ":" +
          std::to_string(tail_certificate_taylor_complete_max) + ":" +
          (current_accuracy_deficit ==
                  std::numeric_limits<std::uint32_t>::max()
              ? "none"
              : std::to_string(current_accuracy_deficit)) + ":" +
          progress);
    };
    search_witness("direct");
    record_phase_attempt("direct", note_phase_progress());
    while (!witness_radius.has_value() &&
           use_physical_source_tail &&
           !last_accuracy_failure.empty() &&
           phase_nonprogress < kPrivateTaylorNonprogressCap &&
           tail_certificate_taylor_complete_max <
               kPrivateTaylorOrderCap) {
      const auto step = std::max(
          kPrivateTaylorOrderMinimumStep,
          tail_certificate_taylor_complete_max / 2);
      const auto next_order_wide =
          static_cast<std::uint64_t>(
              tail_certificate_taylor_complete_max) + step;
      const auto next_order = static_cast<std::uint32_t>(
          std::min<std::uint64_t>(
              next_order_wide, kPrivateTaylorOrderCap));
      auto extended_physical =
          prepare_physical_regular_homogeneous_tail_model(
              *source_equation, acb_incoming->solution(), next_order);
      if (extended_physical.status != TailMajorantStatus::Certified ||
          !extended_physical.model.has_value())
        return ineligible(
            "physical-tail-private-order-reconstruction-failed",
            "retained_taylor_complete_max=" +
                std::to_string(producer_taylor_complete_max) +
                ";requested_tail_certificate_taylor_complete_max=" +
                std::to_string(next_order) + ";detail=" +
                extended_physical.detail);
      physical_source_model =
          std::move(*extended_physical.model);
      tail_certificate_taylor_complete_max = next_order;
      tail_order_attempts.push_back(next_order);
      search_witness("direct");
      record_phase_attempt("direct", note_phase_progress());
    }
    if (!witness_radius.has_value() &&
        use_physical_source_tail) {
      direct_center_failure =
          last_accuracy_failure.empty()
          ? (certificate_failed_after_disk
                 ? "direct receiver-center certificate failed after disk certification"
                 : "direct receiver-center certificate was inconclusive")
          : last_accuracy_failure;

      // The planner already supplies a common point in the safe overlap of
      // both ordinary charts.  If direct source-to-receiver-center
      // evaluation remains too close to the source theorem's certified q
      // disk even at the private-order cap, certify the source only to that
      // overlap and continue with the receiving physical equation.
      producing_point = RealEvaluationPoint::rational(
          exact_match.producing_local.str());
      expected_rim = producing_point.sign < 0
          ? producing_rim : std::nullopt;
      point_modulus =
          exact_path_detail::abs(exact_match.producing_local);
      witness_gap = producing.geometry.radius - point_modulus;
      phase_best_accuracy_excess.reset();
      phase_nonprogress = 0;
      search_witness("overlap");
      record_phase_attempt("overlap", note_phase_progress());
      while (!witness_radius.has_value() &&
             !last_accuracy_failure.empty() &&
             phase_nonprogress < kPrivateTaylorNonprogressCap &&
             tail_certificate_taylor_complete_max <
                 kPrivateTaylorOrderCap) {
        const auto step = std::max(
            kPrivateTaylorOrderMinimumStep,
            tail_certificate_taylor_complete_max / 2);
        const auto next_order = static_cast<std::uint32_t>(
            std::min<std::uint64_t>(
                static_cast<std::uint64_t>(
                    tail_certificate_taylor_complete_max) + step,
                kPrivateTaylorOrderCap));
        auto extended_physical =
            prepare_physical_regular_homogeneous_tail_model(
                *source_equation, acb_incoming->solution(),
                next_order);
        if (extended_physical.status !=
                TailMajorantStatus::Certified ||
            !extended_physical.model.has_value())
          return ineligible(
              "physical-overlap-source-private-order-reconstruction-failed",
              "requested_tail_certificate_taylor_complete_max=" +
                  std::to_string(next_order) + ";detail=" +
                  extended_physical.detail);
        physical_source_model =
            std::move(*extended_physical.model);
        tail_certificate_taylor_complete_max = next_order;
        tail_order_attempts.push_back(next_order);
        search_witness("overlap");
        record_phase_attempt("overlap", note_phase_progress());
      }

      if (witness_radius.has_value()) {
        source_overlap_inflated_value =
            certified_handoff.value;
        const auto receiving_match_local =
            exact_match.receiving_local;
        const auto recentered =
            recenter_physical_cleared_ode(
                *equation->physical_equation(),
                receiving_match_local);
        if (!recentered.eligible ||
            !recentered.equation.has_value())
          return ineligible(
              "physical-overlap-recenter-ineligible",
              recentered.reason);
        recentered_chart_radius =
            receiving.geometry.radius -
            exact_path_detail::abs(receiving_match_local);
        const auto recentered_point_exact =
            Rational(0) - receiving_match_local;
        const auto recentered_point =
            RealEvaluationPoint::rational(
                recentered_point_exact.str());
        const auto recentered_point_modulus =
            exact_path_detail::abs(recentered_point_exact);
        if (!(recentered_point_modulus <
              recentered_chart_radius))
          return ineligible(
              "physical-overlap-recenter-has-no-certified-return-disk");

        std::optional<PhysicalRegularTaylorTailModel>
            recentered_model;
        std::string recentered_model_failure;
        const auto prepare_recentered_model =
            [&](std::uint32_t order) {
          const auto overlap_evolution =
              evolve_ordinary_center_value<ComplexBall>(
                  *recentered.equation,
                  certified_handoff.value, order);
          if (!overlap_evolution.eligible) {
            recentered_model_failure =
                overlap_evolution.reason;
            return false;
          }
          ChartGeometry chart;
          chart.center_exact = "0";
          chart.scale_exact = "1";
          chart.radius_exact = recentered_chart_radius.str();
          chart.radius = ComplexBall::from_strings(
              recentered_chart_radius.str());
          chart.infinite_radius = false;
          auto overlap_solution =
              ordinary_evolution_local_solution(
                  overlap_evolution, std::move(chart),
                  receiving.prescriptions,
                  incoming->checkpoint_identity() +
                      ":overlap-recenter:" + arm_name + ":" +
                      std::to_string(match_index) + ":" +
                      std::to_string(order));
          auto prepared =
              prepare_physical_regular_homogeneous_tail_model(
                  *recentered.equation, overlap_solution);
          if (prepared.status !=
                  TailMajorantStatus::Certified ||
              !prepared.model.has_value()) {
            recentered_model_failure = prepared.detail;
            return false;
          }
          recentered_model = std::move(*prepared.model);
          recentered_taylor_complete_max = order;
          return true;
        };

        std::string recentered_accuracy_failure;
        bool recentered_certificate_failed_after_disk = false;
        const auto try_recentered_witness =
            [&](const Rational& candidate, std::string direction,
                std::uint32_t exponent) {
          EvaluationOptions options;
          options.compute_tail_estimate = false;
          auto attempted =
              evaluate_physical_local_solution_with_certified_tail(
                  *recentered_model, recentered_point,
                  candidate.str(), options);
          if (attempted.tail.status !=
                  TailMajorantStatus::Certified) {
            if (attempted.tail.status !=
                    TailMajorantStatus::Inconclusive ||
                attempted.tail.disk.status !=
                    TailMajorantStatus::Inconclusive)
              recentered_certificate_failed_after_disk = true;
            return false;
          }
          if (attempted.evaluation.imaginary_sign.has_value() ||
              attempted.evaluation.value.epsilon.min_power <
                  execution_epsilon_min ||
              attempted.evaluation.value.epsilon.complete_max !=
                  incoming_epsilon_max ||
              attempted.tail.value.guarantee !=
                  ErrorGuarantee::Certified ||
              !tail_majorant_detail::same_epsilon_window(
                  attempted.tail.value.frame,
                  EpsilonWindow{execution_epsilon_min,
                                incoming_epsilon_max}))
            throw std::logic_error(
                "recentered physical return changed its epsilon frame or branch");
          auto candidate_retained =
              attempted.evaluation.value;
          candidate_retained.error = ErrorEnvelope{};
          // evaluate_local_solution canonically trims leading epsilon rows
          // that are structurally zero across the reconstructed Taylor
          // tensor.  The transport contract still owns the wider execution
          // frame: restore that exact-zero prefix before attaching the
          // full-frame tail certificate or publishing the receiving value.
          if (execution_epsilon_min <
              candidate_retained.epsilon.min_power) {
            auto widened = physical_ode_detail::zero_epsilon_vector(
                EpsilonWindow{execution_epsilon_min,
                              incoming_epsilon_max},
                candidate_retained.dimension);
            for (std::int64_t raw_power =
                     candidate_retained.epsilon.min_power;
                 raw_power <=
                     candidate_retained.epsilon.complete_max;
                 ++raw_power) {
              const auto power =
                  static_cast<std::int32_t>(raw_power);
              for (std::uint32_t component = 0;
                   component < candidate_retained.dimension;
                   ++component)
                widened.at(power, component) =
                    candidate_retained.at(power, component);
            }
            candidate_retained = std::move(widened);
          }
          auto candidate_handoff = attempted.evaluation;
          candidate_handoff.value = candidate_retained;
          for (std::int64_t raw_power =
                   execution_epsilon_min;
               raw_power <= incoming_epsilon_max;
               ++raw_power) {
            const auto power =
                static_cast<std::int32_t>(raw_power);
            const auto row = static_cast<std::size_t>(
                raw_power - execution_epsilon_min);
            for (std::uint32_t component = 0;
                 component <
                     candidate_handoff.value.dimension;
                 ++component)
              attempted.tail.value.absolute[row].add_error_to(
                  candidate_handoff.value.at(
                      power, component));
          }
          candidate_handoff.value.error = ErrorEnvelope{};
          if (const auto failure =
                  value_handoff_accuracy_failure(
                      candidate_handoff.value,
                      required_complete_max,
                      relative_accuracy_max)) {
            recentered_accuracy_failure =
                *failure + ";witness_radius_exact=" +
                candidate.str() + ";witness_direction=" +
                direction + ";witness_dyadic_exponent=" +
                std::to_string(exponent);
            return false;
          }
          recentered_witness_radius = candidate;
          recentered_witness_dyadic_direction =
              std::move(direction);
          recentered_witness_dyadic_exponent = exponent;
          recentered_retained_value =
              std::move(candidate_retained);
          recentered_center_handoff =
              std::move(candidate_handoff);
          recentered_certified = std::move(attempted);
          return true;
        };
        const auto search_recentered_witness = [&] {
          recentered_witness_radius.reset();
          recentered_witness_dyadic_direction.clear();
          recentered_witness_dyadic_exponent = 0;
          recentered_certified.reset();
          recentered_retained_value = EpsilonVector{};
          recentered_center_handoff = LocalEvaluation{};
          recentered_accuracy_failure.clear();
          recentered_certificate_failed_after_disk = false;
          const auto gap =
              recentered_chart_radius -
              recentered_point_modulus;
          for (std::uint32_t exponent = kWitnessSearchCap;
               exponent >= 1 &&
                   !recentered_witness_radius.has_value();
               --exponent) {
            const auto candidate =
                recentered_chart_radius -
                gap / dyadic_denominator(exponent);
            (void)try_recentered_witness(
                candidate, "outward", exponent);
          }
          for (std::uint32_t exponent = 2;
               exponent <= kWitnessSearchCap &&
                   !recentered_witness_radius.has_value();
               ++exponent) {
            const auto candidate =
                recentered_point_modulus +
                gap / dyadic_denominator(exponent);
            (void)try_recentered_witness(
                candidate, "inward", exponent);
          }
        };

        auto recentered_order =
            receiver_taylor_complete_max;
        recentered_order_attempts.push_back(
            recentered_order);
        if (!prepare_recentered_model(recentered_order))
          return ineligible(
              "physical-overlap-recenter-model-ineligible",
              recentered_model_failure);
        search_recentered_witness();
        while (!recentered_witness_radius.has_value() &&
               !recentered_accuracy_failure.empty() &&
               recentered_order <
                   kPrivateTaylorOrderCap) {
          const auto step = std::max(
              kPrivateTaylorOrderMinimumStep,
              recentered_order / 2);
          recentered_order =
              static_cast<std::uint32_t>(
                  std::min<std::uint64_t>(
                      static_cast<std::uint64_t>(
                          recentered_order) + step,
                      kPrivateTaylorOrderCap));
          recentered_order_attempts.push_back(
              recentered_order);
          if (!prepare_recentered_model(recentered_order))
            return ineligible(
                "physical-overlap-recenter-model-ineligible",
                recentered_model_failure);
          search_recentered_witness();
        }
        if (!recentered_witness_radius.has_value()) {
          std::string detail =
              "direct_center_failure=" +
              direct_center_failure +
              ";recentered_taylor_order_attempts=";
          for (std::size_t index = 0;
               index < recentered_order_attempts.size();
               ++index) {
            if (index != 0) detail += ",";
            detail += std::to_string(
                recentered_order_attempts[index]);
          }
          if (!recentered_accuracy_failure.empty())
            detail += ";recentered_accuracy_failure=" +
                recentered_accuracy_failure;
          else if (recentered_certificate_failed_after_disk)
            detail +=
                ";recentered_certificate_failed_after_disk=true";
          return ineligible(
              "physical-overlap-recenter-return-certificate-failed",
              std::move(detail));
        }
        certified_handoff =
            recentered_center_handoff;
        used_overlap_recenter = true;
      }
    }
    const auto tail_order_attempt_detail = [&] {
      std::string detail =
          "retained_taylor_complete_max=" +
          std::to_string(producer_taylor_complete_max) +
          ";tail_certificate_taylor_order_attempts=";
      for (std::size_t index = 0;
           index < tail_order_attempts.size(); ++index) {
        if (index != 0) detail += ",";
        detail += std::to_string(tail_order_attempts[index]);
      }
      detail += ";tail_certificate_accuracy_deficits=";
      for (std::size_t index = 0;
           index < tail_order_accuracy_trace.size(); ++index) {
        if (index != 0) detail += ",";
        detail += tail_order_accuracy_trace[index];
      }
      return detail;
    };
    if (!witness_radius.has_value()) {
      if (certificate_failed_after_disk)
        return ineligible(
            "receiver-center-tail-certificate-fails-after-disk-certification",
            tail_order_attempt_detail());
      if (!last_accuracy_failure.empty())
        return ineligible(
            "inflated-center-evaluation-fails-relative-accuracy-contract",
            last_accuracy_failure + ";" +
                tail_order_attempt_detail());
      return ineligible(
          "receiver-center-tail-certificate-is-inconclusive",
          tail_order_attempt_detail());
    }

    auto metadata = parse_local_metadata(metadata_prototype);
    const auto same_prescription = [](const Prescription& left,
                                      const Prescription& right) {
      return left.factor_exact == right.factor_exact &&
             left.sign == right.sign &&
             left.multiplicity == right.multiplicity &&
             left.leading_coefficient_sign ==
                 right.leading_coefficient_sign;
    };
    if (!(Rational(metadata.chart.center_exact) ==
              receiving.geometry.center) ||
        !(Rational(metadata.chart.scale_exact) ==
              receiving.geometry.scale) ||
        metadata.chart.infinite_radius ||
        !acb_equal(metadata.chart.radius.raw(),
                   ComplexBall::from_strings(
                       receiving.geometry.radius.str()).raw()) ||
        metadata.prescriptions.size() != receiving.prescriptions.size() ||
        !std::equal(metadata.prescriptions.begin(),
                    metadata.prescriptions.end(),
                    receiving.prescriptions.begin(), same_prescription) ||
        metadata.a.domain != ExactDomain::Rational ||
        metadata.b.domain != ExactDomain::Rational ||
        !(Rational(metadata.a.canonical) == Rational(0)) ||
        !(Rational(metadata.b.canonical) == Rational(0)))
      throw std::invalid_argument(
          "physical value-solver metadata differs from its retained ordinary plan chart");

    const auto kernel_started = std::chrono::steady_clock::now();
    auto evolution = evolve_ordinary_center_value<ComplexBall>(
        *equation->physical_equation(), certified_handoff.value,
        receiver_taylor_complete_max);
    if (!evolution.eligible)
      return ineligible("physical-evolution-ineligible: " +
                        evolution.reason);
    if (evolution.epsilon.min_power != execution_epsilon_min ||
        evolution.epsilon.complete_max != incoming_epsilon_max ||
        evolution.taylor_coefficients.size() !=
            static_cast<std::size_t>(receiver_taylor_complete_max) + 1)
      throw std::logic_error(
          "ordinary physical evolution changed its zero-padded source epsilon window");
    const auto kernel_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - kernel_started).count();

    std::string local_handle;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      if (session->locals.size() + session->pending_local_solves + 1 >
          session->local_capacity)
        throw std::invalid_argument(
            "persistent local capacity is exhausted by physical value transport hop");
      local_handle = "l:" + std::to_string(session->next_local++);
      ++session->pending_local_solves;
    }

    const auto started = std::chrono::steady_clock::now();
    std::shared_ptr<StoredLocalBase> next;
    json::object sealed_lineage;
    try {
      const auto result_checkpoint = arm_checkpoint_identity(
          checkpoint_root, arm_name, "local", match_index + 1);
      metadata.checkpoint_identity = result_checkpoint;
      LocalSolution<ComplexBall> solution;
      solution.chart = std::move(metadata.chart);
      solution.epsilon = evolution.epsilon;
      solution.taylor_complete_max = receiver_taylor_complete_max;
      solution.dimension = evolution.dimension;
      solution.prescriptions = std::move(metadata.prescriptions);
      solution.checkpoint_identity = result_checkpoint;
      LocalSector<ComplexBall> sector;
      sector.a = std::move(metadata.a);
      sector.b = std::move(metadata.b);
      sector.log_power = 0;
      sector.coefficients.reserve(solution.sector_size());
      for (std::int64_t raw_power = solution.epsilon.min_power;
           raw_power <= solution.epsilon.complete_max; ++raw_power) {
        const auto power = static_cast<std::int32_t>(raw_power);
        for (std::uint32_t taylor = 0;
             taylor <= receiver_taylor_complete_max; ++taylor)
          for (std::uint32_t component = 0;
               component < solution.dimension; ++component)
            sector.coefficients.push_back(
                evolution.at(taylor).at(power, component));
      }
      solution.sectors.push_back(std::move(sector));
      validate_local_solution(solution, false);

      const auto owner_reference = [](
          const RetainedPlanChartBinding& binding,
          const PhysicalEquationOwnerBase& owner) {
        return json::object{
            {"kind", owner.equation_owner_kind()},
            {"handle", owner.equation_owner_handle()},
            {"operator_identity", owner.equation_operator_identity()},
            {"plan_exact_identity", binding.exact_identity},
            {"owner_signature_identity", owner.owner_signature_identity()},
            {"physical_payload_identity", owner.physical_payload_identity()}};
      };
      const auto incoming_record = json::object{
          {"local", incoming->handle()},
          {"chart", incoming->source_chart()},
          {"source_operator_identity", incoming->source_operator_identity()},
          {"checkpoint_identity", incoming->checkpoint_identity()},
          {"epsilon", json::object{{"min", incoming_epsilon_min},
                                     {"max", incoming_epsilon_max}}},
          {"taylor_complete_max", producer_taylor_complete_max},
          {"top_valid", incoming_summary.at("top_valid")},
          {"analytic_metadata", incoming->exact_analytic_metadata()}};
      const auto physical_certificate_record =
          [](const CertifiedLocalEvaluation& value) {
        return json::object{
            {"status", tail_majorant_status_name(
                 value.tail.status)},
            {"value", checkpoint_error_envelope_record(
                 value.tail.value)},
            {"theta", checkpoint_error_envelope_record(
                 value.tail.theta)},
            {"disk", json::object{
                 {"witness_radius_exact",
                  value.tail.disk.witness_radius_exact},
                 {"q_lower_exact",
                  value.tail.disk.q_lower.dump_exact()},
                 {"ode_norm_upper_exact",
                  value.tail.disk.ode_norm_upper.dump_exact()},
                 {"cauchy_circle_upper_exact", [&] {
                    json::array values;
                    for (const auto& bound :
                         value.tail.disk.cauchy_circle_upper)
                      values.emplace_back(
                          bound.dump_exact());
                    return values;
                  }()},
                 {"detail", value.tail.disk.detail}}},
            {"detail", value.tail.detail}};
      };
      const auto physical_inflation_record =
          [](const EpsilonVector& retained,
             const EpsilonVector& inflated) {
        return json::object{
            {"gate",
             "each-component-acb-add-error-mag-by-epsilon-row-v1"},
            {"retained_value",
             checkpoint_epsilon_vector_record(retained)},
            {"inflated_value",
             checkpoint_epsilon_vector_record(inflated)}};
      };
      json::object tail_contract;
      if (!use_physical_source_tail) {
        tail_contract = json::object{
            {"mode", "physical-evolution-transient-tail-v2"},
            {"source_model", checkpoint_regular_tail_model_record(
                 *acb_incoming->tail_model().model)},
            {"witness_radius_exact", witness_radius->str()},
            {"witness_search_policy",
             "accuracy-qualified-bidirectional-dyadic-v2"},
            {"witness_dyadic_direction",
             witness_dyadic_direction},
            {"witness_dyadic_exponent", witness_dyadic_exponent},
            {"output_status", "transient-physical-reconstructible"},
            {"output_reason",
             "causal physical evolution retains an exact physical q/C owner "
             "for certified transient tail reconstruction"}};
      } else if (used_overlap_recenter) {
        tail_contract = json::object{
            {"mode",
             "certified-physical-overlap-recenter-point-tail-acb-v1"},
            {"producer_taylor_complete_max",
             tail_certificate_taylor_complete_max},
            {"receiver_taylor_complete_max",
             receiver_taylor_complete_max},
            {"producing_point_exact",
             exact_match.producing_local.str()},
            {"producing_chart_radius_exact",
             producing.geometry.radius.str()},
            {"witness_radius_exact",
             witness_radius->str()},
            {"witness_search_policy",
             "accuracy-qualified-bidirectional-dyadic-v2"},
            {"witness_dyadic_direction",
             witness_dyadic_direction},
            {"witness_dyadic_exponent",
             witness_dyadic_exponent},
            {"source_physical_payload_identity",
             physical_source_model->physical_payload_identity},
            {"source_certificate",
             physical_certificate_record(*certified)},
            {"source_inflation",
             physical_inflation_record(
                 retained_value,
                 source_overlap_inflated_value)},
            {"receiving_match_local_exact",
             exact_match.receiving_local.str()},
            {"recentered_chart_radius_exact",
             recentered_chart_radius.str()},
            {"recentered_taylor_complete_max",
             recentered_taylor_complete_max},
            {"recentered_witness_radius_exact",
             recentered_witness_radius->str()},
            {"recentered_witness_dyadic_direction",
             recentered_witness_dyadic_direction},
            {"recentered_witness_dyadic_exponent",
             recentered_witness_dyadic_exponent},
            {"receiving_physical_payload_identity",
             equation->physical_payload_identity()},
            {"recentered_certificate",
             physical_certificate_record(
                 *recentered_certified)},
            {"recentered_inflation",
             physical_inflation_record(
                 recentered_retained_value,
                 recentered_center_handoff.value)},
            {"direct_center_failure",
             direct_center_failure},
            {"output_status",
             "transient-physical-reconstructible"},
            {"output_reason",
             "causal physical evolution retains an exact physical q/C owner "
             "for certified transient tail reconstruction"}};
      } else {
        tail_contract = json::object{
            {"mode",
             "certified-physical-regular-taylor-point-tail-acb-v1"},
            {"producer_taylor_complete_max",
             tail_certificate_taylor_complete_max},
            {"receiver_taylor_complete_max",
             receiver_taylor_complete_max},
            {"producing_point_exact", producing_local.str()},
            {"producing_chart_radius_exact",
             producing.geometry.radius.str()},
            {"witness_radius_exact", witness_radius->str()},
            {"witness_search_policy",
             "accuracy-qualified-bidirectional-dyadic-v2"},
            {"witness_dyadic_direction",
             witness_dyadic_direction},
            {"witness_dyadic_exponent", witness_dyadic_exponent},
            {"source_physical_payload_identity",
             physical_source_model->physical_payload_identity},
            {"certificate",
             physical_certificate_record(*certified)},
            {"inflation",
             physical_inflation_record(
                 retained_value,
                 certified_handoff.value)},
            {"output_status", "transient-physical-reconstructible"},
            {"output_reason",
             "causal physical evolution retains an exact physical q/C owner "
             "for certified transient tail reconstruction"}};
      }
      json::object derivation{
          {"schema", "diffexp2-retained-plan-value-handoff-v2"},
          {"capability", "retained-native-regular-value-handoff-v2"},
          {"tile_plan", plan->handle()},
          {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
          {"tile_plan_provenance_identity", plan->provenance_identity()},
          {"arm", arm_name}, {"match", match_index},
          {"producing", json::object{
               {"chart", encode_plan_chart(
                    producing, exact_match.producing_chart)},
               {"owner", owner_reference(
                    producing, *producing_equation_owner)}}},
          {"receiving", json::object{
               {"chart", encode_plan_chart(
                    receiving, exact_match.receiving_chart)},
               {"owner", owner_reference(receiving, **receiving_owner)}}},
          {"receiver_center_physical_exact",
           receiving.local_geometry.center_exact},
          {"producing_local_exact", producing_local.str()},
          {"prototype_identity", json::serialize(
               canonical_json_value(value_solver))},
          {"tail_contract", std::move(tail_contract)},
          {"accuracy_contract", json::object{
               {"relative_error_max_exact", relative_accuracy_text},
               {"gate",
                "required-prefix-component-radii-lte-threshold-times-max-one-upper-magnitude-v2"},
               {"acb_preflight_required", true}}},
          {"epsilon", json::object{
               {"min", requested_epsilon.min_power},
               {"max", requested_epsilon.complete_max},
               {"required_complete_max", required_complete_max}}},
          {"incoming", incoming_record},
          {"scope", "single-regular-to-regular-transport-hop"},
          {"coefficient_transport", "native-retained-only"},
          {"whole_arm_complete", false},
          {"evaluated_epsilon", json::object{
               {"min", execution_epsilon_min},
               {"max", incoming_epsilon_max},
               {"required_complete_max", required_complete_max}}},
          {"output", json::object{
               {"checkpoint_identity", solution.checkpoint_identity},
               {"chart", equation->handle()},
               {"source_operator_identity", equation->exact_identity()},
               {"epsilon", json::object{
                    {"min", solution.epsilon.min_power},
                    {"max", solution.epsilon.complete_max}}},
               {"taylor_complete_max", solution.taylor_complete_max},
               {"top_valid", encode_validity(incoming_epsilon_max)},
               {"dimension", solution.dimension}}},
          {"equation_owner_signature_identity",
           equation->owner_signature_identity()},
          {"equation_payload_identity",
           equation->physical_payload_identity()}};
      derivation["provenance_identity"] = json::serialize(
          canonical_json_value(derivation));
      struct PhysicalValueHandoffOwners {
        std::shared_ptr<StoredTilePlan> plan;
        std::shared_ptr<StoredLocalBase> incoming;
        std::shared_ptr<RegularPhysicalEquationOwnerBase> receiving;
      };
      auto derivation_owner = std::make_shared<PhysicalValueHandoffOwners>(
          PhysicalValueHandoffOwners{plan, incoming, *receiving_owner});
      const NativeLocalDiagnostics diagnostics{
          incoming_epsilon_max, 0.0, kernel_ms};
      next = make_retained_typed_shared<ComplexBall,
          StoredLocal<ComplexBall>>(
              local_handle, equation->handle(), equation->exact_identity(),
              std::move(solution), equation->precision_bits(),
              std::vector<PseudoHit<ComplexBall>>{}, diagnostics,
              std::nullopt, std::move(derivation), derivation_owner,
              unavailable_tail_model(
                  "ordinary physical evolution has no framed recurrence tail "
                  "model; its retained physical q/C owner supports certified "
                  "transient reconstruction"),
              std::nullopt, true, true, equation,
              equation->physical_equation());
      sealed_lineage = next->seal_plan_match_lineage();
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "physical value transport hop reservation accounting underflow");
      --session->pending_local_solves;
      throw;
    }

    const auto next_stats = next->stats();
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "physical value transport hop reservation accounting underflow");
      --session->pending_local_solves;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during physical value transport hop");
      const auto found = session->locals.find(incoming_handle);
      if (found == session->locals.end() ||
          found->second.get() != incoming.get())
        throw std::invalid_argument(
            "physical value transport hop incoming owner changed before publication");
      if (!session->locals.emplace(local_handle, next).second)
        throw std::logic_error(
            "physical value transport hop local handle collision");
      ++session->total_local_solves;
      ++session->total_transport_physical_value_hop_successes;
      session->total_local_run_parse_ms += next_stats.create_parse_ms;
      session->total_local_kernel_ms += next_stats.create_kernel_ms;
      plan->note_match_advance(arm_name);
    }
    auto next_summary = compact_transport_local_reference(next);
    next_summary["release_via"] = "local";
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"capability", "ordinary-physical-value-transport-hop-v1"},
        {"execution_mode", "causal-ordinary-physical-evolution"},
        {"native_retained", true}, {"json_coefficients", 0},
        {"used", true},
        {"reason", "eligible-causal-ordinary-physical-evolution"},
        {"arm", arm_name}, {"match", match_index},
        {"retained_taylor_complete_max",
         producer_taylor_complete_max},
        {"tail_certificate_taylor_complete_max",
         tail_certificate_taylor_complete_max},
        {"tail_order_retries",
         tail_order_attempts.size() - 1},
        {"used_overlap_recenter", used_overlap_recenter},
        {"recentered_tail_certificate_taylor_complete_max",
         used_overlap_recenter
             ? json::value(recentered_taylor_complete_max)
             : json::value(nullptr)},
        {"recentered_tail_order_retries",
         used_overlap_recenter
             ? json::value(
                   recentered_order_attempts.size() - 1)
             : json::value(nullptr)},
        {"tail_witness_direction", witness_dyadic_direction},
        {"tail_witness_dyadic_exponent", witness_dyadic_exponent},
        {"value_hops", 1}, {"basis_matches", 0},
        {"output_tail_status", "transient-physical-reconstructible"},
        {"next_hop_policy", "certified-physical-tail-replay"},
        {"next_local", std::move(next_summary)},
        {"sealed_local_lineage", std::move(sealed_lineage)},
        {"elapsed_ms", std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - started).count()}};
  }

  if (operation == "transport.consume_value_hop") {
    require_exact_keys(
        root,
        {"schema", "op", "session", "tile_plan",
         "tile_plan_checkpoint_identity", "arm", "match", "value_solver",
         "incoming", "incoming_checkpoint_identity", "epsilon",
         "checkpoint_policy"},
        "native transport.consume_value_hop request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "transport.consume_value_hop requires Rational or Acb coefficients");
    const auto arm_name = required_string(root, "arm");
    const auto match_index = checkpoint_size_t(
        root.at("match"), "value transport hop match index");
    const auto incoming_handle = required_string(root, "incoming");
    const auto& value_solver = as_object(
        root.at("value_solver"), "regular value-solver prototype");
    require_exact_keys(value_solver,
                       {"schema", "run", "metadata",
                        "tail_proxy_max_exact",
                        "relative_accuracy_max_exact"},
                       "regular value-solver prototype");
    if (required_string(value_solver, "schema") !=
            "diffexp2-native-regular-value-solver-prototype-v1")
      throw std::invalid_argument(
          "unsupported regular value-solver prototype schema");
    const auto& run_prototype = as_object(
        value_solver.at("run"), "regular value-solver run prototype");
    require_exact_keys(
        run_prototype,
        {"nmax", "p", "has_initial", "adaptive_probe", "a_target",
         "b_target", "a_shift_min", "a_shifts", "schedule", "initial",
         "initial_validity", "source", "return_u"},
        "regular value-solver run prototype");
    const auto& metadata_prototype = as_object(
        value_solver.at("metadata"),
        "regular value-solver metadata prototype");
    const Rational tail_proxy_max(
        required_string(value_solver, "tail_proxy_max_exact"));
    const Rational relative_accuracy_max(
        required_string(value_solver, "relative_accuracy_max_exact"));
    if (tail_proxy_max.sign() <= 0 || !(tail_proxy_max < Rational(1)) ||
        relative_accuracy_max.sign() <= 0 ||
        !(relative_accuracy_max < Rational(1)) ||
        required_string(value_solver, "tail_proxy_max_exact") !=
            tail_proxy_max.str() ||
        required_string(value_solver, "relative_accuracy_max_exact") !=
            relative_accuracy_max.str())
      throw std::invalid_argument(
          "regular value-solver accuracy thresholds must be canonical exact rationals strictly between zero and one");

    const auto& epsilon = as_object(
        root.at("epsilon"), "value transport hop epsilon contract");
    require_exact_keys(epsilon, {"min", "max", "required_complete_max"},
                       "value transport hop epsilon contract");
    EpsilonWindow requested_epsilon{
        as_i32(epsilon.at("min"), "value hop epsilon minimum"),
        as_i32(epsilon.at("max"), "value hop epsilon maximum")};
    (void)requested_epsilon.width();
    const auto required_complete_max = as_i32(
        epsilon.at("required_complete_max"),
        "value hop required epsilon maximum");
    if (required_complete_max < requested_epsilon.min_power ||
        required_complete_max > requested_epsilon.complete_max)
      throw std::invalid_argument(
          "value transport hop epsilon contract is inconsistent");
    const auto& checkpoint_policy = as_object(
        root.at("checkpoint_policy"),
        "value transport hop checkpoint policy");
    require_exact_keys(checkpoint_policy, {"schema", "root"},
                       "value transport hop checkpoint policy");
    if (required_string(checkpoint_policy, "schema") !=
            "diffexp2-deterministic-arm-checkpoints-v1")
      throw std::invalid_argument(
          "unsupported value transport hop checkpoint policy");
    const auto checkpoint_root = required_string(checkpoint_policy, "root");

    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> incoming;
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
            "value transport hop tile-plan binding is stale");
      plan = plan_found->second;
      const auto& retained = plan->arm(arm_name);
      if (match_index >= retained.exact.matches.size())
        throw std::invalid_argument(
            "value transport hop match lies outside its exact arm");
      const auto incoming_found = session->locals.find(incoming_handle);
      if (incoming_found == session->locals.end() ||
          required_string(root, "incoming_checkpoint_identity") !=
              incoming_found->second->checkpoint_identity())
        throw std::invalid_argument(
            "value transport hop incoming-local binding is stale");
      incoming = incoming_found->second;
    }

    const auto ineligible = [&](std::string reason,
                                std::string detail = {}) {
      json::object response{
          {"status", "ok"}, {"session", session->handle},
          {"capability", "regular-value-transport-hop-v1"},
          {"used", false}, {"reason", std::move(reason)},
          {"arm", arm_name}, {"match", match_index},
          {"value_hops", 0}, {"basis_matches", 0}};
      if (!detail.empty())
        response["detail"] = std::move(detail);
      return response;
    };
    const auto& retained = plan->arm(arm_name);
    const auto& exact_match = retained.exact.matches[match_index];
    const auto& producing = retained.charts.at(exact_match.producing_chart);
    const auto& receiving = retained.charts.at(exact_match.receiving_chart);
    if (producing.geometry.singular_center ||
        receiving.geometry.singular_center)
      return ineligible("singular-chart-crossing");
    const auto* producing_owner = std::get_if<
        std::shared_ptr<PreparedChartBase>>(&producing.owner);
    const auto* receiving_owner = std::get_if<
        std::shared_ptr<PreparedChartBase>>(&receiving.owner);
    if (!producing_owner || !*producing_owner)
      return ineligible("producing-owner-is-not-a-primitive-regular-chart");
    if (!receiving_owner || !*receiving_owner)
      return ineligible("receiving-owner-is-not-a-primitive-regular-chart");
    const auto& owner_relative_accuracy_max =
        (*receiving_owner)->regular_value_relative_accuracy_max_exact();
    if (!owner_relative_accuracy_max.has_value())
      return ineligible(
          "receiving-owner-has-no-relative-accuracy-contract");
    if (required_string(value_solver, "relative_accuracy_max_exact") !=
        *owner_relative_accuracy_max)
      throw std::invalid_argument(
          "regular value-solver relative-accuracy threshold differs from its prepared owner");
    if (required_string(value_solver, "tail_proxy_max_exact") !=
        (*receiving_owner)->regular_value_tail_proxy_max_exact())
      throw std::invalid_argument(
          "regular value-solver tail threshold differs from its prepared owner");
    // The value shortcut evaluates at the receiving physical centre.  A
    // rational planning surrogate is never authority for that point.  Keep
    // the optimization on exact rational retained maps; algebraic geometry
    // falls back side-effect-free to the full certified basis match.
    std::optional<Rational> producing_center;
    std::optional<Rational> producing_scale;
    std::optional<Rational> receiving_center;
    try {
      producing_center = Rational(
          producing.local_geometry.center_exact);
      producing_scale = Rational(
          producing.local_geometry.scale_exact);
      receiving_center = Rational(
          receiving.local_geometry.center_exact);
    } catch (const std::invalid_argument&) {
      return ineligible(
          "certified-algebraic-chart-requires-basis-match");
    }
    if (!(*producing_center == producing.geometry.center) ||
        !(*producing_scale == producing.geometry.scale) ||
        !(*receiving_center == receiving.geometry.center) ||
        !acb_equal(producing.local_geometry.radius.raw(),
                   ComplexBall::from_strings(
                       producing.geometry.radius.str()).raw()) ||
        !acb_equal(receiving.local_geometry.radius.raw(),
                   ComplexBall::from_strings(
                       receiving.geometry.radius.str()).raw()))
      return ineligible(
          "certified-algebraic-chart-requires-basis-match");
    if (session->domain != "acb")
      return ineligible(
          "rational-value-handoff-has-no-exact-polynomial-tail-zero-certificate");
    if (producing_scale->is_zero())
      throw std::invalid_argument(
          "value transport hop producing chart has a zero exact scale");
    const auto producing_local =
        (*receiving_center - *producing_center) / *producing_scale;
    const auto center_ratio =
        exact_path_detail::abs(producing_local) /
        producing.geometry.radius;
    if (!(center_ratio < Rational(1)))
      return ineligible(
          "receiver-center-lies-outside-producing-certified-disk");
    const auto incoming_summary = incoming->summary();
    const auto incoming_epsilon_min = as_i32(
        incoming_summary.at("epsilon_min"),
        "value-hop incoming epsilon minimum");
    const auto incoming_epsilon_max = as_i32(
        incoming_summary.at("epsilon_max"),
        "value-hop incoming epsilon maximum");
    const auto incoming_top_valid = parse_validity(
        incoming_summary.at("top_valid"));
    const auto incoming_effective_max = std::min(
        incoming_epsilon_max, incoming_top_valid);
    const auto receiver_frame_top = static_cast<std::int64_t>(
        (*receiving_owner)->frame_base()) +
        static_cast<std::int64_t>((*receiving_owner)->frame_width()) - 1;
    const auto producer_taylor_complete_max = as_u32(
        incoming_summary.at("taylor_complete_max"),
        "value-hop producer Taylor maximum");
    const auto receiver_taylor_complete_max = as_u32(
        run_prototype.at("nmax"),
        "regular value-solver receiver Taylor maximum");
    if (incoming->source_chart() != producing.handle ||
        incoming->source_operator_identity() != producing.exact_identity)
      throw std::invalid_argument(
          "value transport hop incoming local differs from its producing plan chart");
    incoming->require_exact_plan_binding(
        producing.local_geometry, producing.prescriptions,
        "value transport hop incoming local");
    const auto incoming_equation_owner = incoming->retained_equation_owner();
    if (!incoming_equation_owner ||
        incoming_equation_owner.get() != producing_owner->get())
      return ineligible(
          "incoming-equation-owner-is-not-the-producing-plan-owner");
    if ((*receiving_owner)->equation_owner_handle() != receiving.handle ||
        (*receiving_owner)->equation_operator_identity() !=
            receiving.exact_identity ||
        (*receiving_owner)->owner_signature_identity().empty() ||
        (*receiving_owner)->physical_payload_identity().empty())
      return ineligible(
          "receiving-owner-has-no-complete-physical-value-contract");

    const auto producing_point = RealEvaluationPoint::rational(
        producing_local.str());
    const auto producing_rim = exact_plan_rim(
        producing.prescriptions, *producing_scale);
    if (incoming_effective_max < requested_epsilon.complete_max)
      return ineligible("incoming-local-lacks-complete-epsilon-coverage");
    if (requested_epsilon.min_power < (*receiving_owner)->frame_base() ||
        incoming_epsilon_min < (*receiving_owner)->frame_base() ||
        static_cast<std::int64_t>(incoming_effective_max) >
            receiver_frame_top)
      return ineligible(
          "incoming-local-does-not-fit-receiver-epsilon-frame");
    auto acb_incoming =
        std::dynamic_pointer_cast<StoredLocal<ComplexBall>>(incoming);
    if (!acb_incoming ||
        acb_incoming->tail_model().status !=
            TailMajorantStatus::Certified ||
        !acb_incoming->tail_model().model.has_value())
      return ineligible(
          "incoming-local-has-no-certified-regular-tail-model");
    constexpr std::uint32_t kWitnessSearchCap = 16;
    std::optional<Rational> witness_radius;
    std::uint32_t witness_dyadic_exponent = 0;
    std::string witness_dyadic_direction;
    const auto point_modulus = exact_path_detail::abs(producing_local);
    const auto witness_gap = producing.geometry.radius - point_modulus;
    const auto expected_rim = producing_point.sign < 0
        ? producing_rim : std::nullopt;
    std::optional<CertifiedLocalEvaluation> certified;
    EpsilonVector retained_value;
    LocalEvaluation certified_handoff;
    std::string last_accuracy_failure;
    bool certificate_failed_after_disk = false;
    const auto try_witness = [&](const Rational& candidate,
                                 std::string direction,
                                 std::uint32_t exponent) {
      auto attempted =
          incoming->evaluate_retained_point_with_certified_tail(
              producing_point, producing_rim, candidate.str());
      if (!attempted.has_value()) return false;
      if (attempted->tail.status != TailMajorantStatus::Certified) {
        if (attempted->tail.status !=
                TailMajorantStatus::Inconclusive ||
            attempted->tail.disk.status !=
                TailMajorantStatus::Inconclusive)
          certificate_failed_after_disk = true;
        return false;
      }
      if (attempted->evaluation.imaginary_sign != expected_rim)
        throw std::logic_error(
            "center-evaluation rim differs from its producing plan chart");
      if (attempted->evaluation.value.epsilon.min_power !=
              incoming_epsilon_min ||
          attempted->evaluation.value.epsilon.complete_max !=
              incoming_epsilon_max ||
          attempted->tail.value.guarantee !=
              ErrorGuarantee::Certified ||
          !tail_majorant_detail::same_epsilon_window(
              attempted->evaluation.value.epsilon,
              attempted->tail.value.frame) ||
          attempted->tail.value.absolute.size() !=
              attempted->evaluation.value.epsilon.width())
        throw std::logic_error(
            "certified value handoff changed its retained epsilon frame");
      auto candidate_retained = attempted->evaluation.value;
      candidate_retained.error = ErrorEnvelope{};
      auto candidate_handoff = attempted->evaluation;
      for (std::int64_t raw_power =
               candidate_handoff.value.epsilon.min_power;
           raw_power <=
               candidate_handoff.value.epsilon.complete_max;
           ++raw_power) {
        const auto power = static_cast<std::int32_t>(raw_power);
        const auto row = static_cast<std::size_t>(
            raw_power -
            candidate_handoff.value.epsilon.min_power);
        for (std::uint32_t component = 0;
             component < candidate_handoff.value.dimension;
             ++component)
          attempted->tail.value.absolute[row].add_error_to(
              candidate_handoff.value.at(power, component));
      }
      candidate_handoff.value.error = ErrorEnvelope{};
      if (const auto failure = value_handoff_accuracy_failure(
              candidate_handoff.value, required_complete_max,
              relative_accuracy_max)) {
        last_accuracy_failure =
            *failure + ";witness_radius_exact=" +
            candidate.str() + ";witness_direction=" + direction +
            ";witness_dyadic_exponent=" +
            std::to_string(exponent);
        return false;
      }
      witness_radius = candidate;
      witness_dyadic_direction = std::move(direction);
      witness_dyadic_exponent = exponent;
      retained_value = std::move(candidate_retained);
      certified_handoff = std::move(candidate_handoff);
      certified = std::move(attempted);
      return true;
    };
    const auto dyadic_denominator = [](std::uint32_t exponent) {
      Rational denominator(1);
      for (std::uint32_t index = 0; index < exponent; ++index)
        denominator *= Rational(2);
      return denominator;
    };
    for (std::uint32_t exponent = kWitnessSearchCap;
         exponent >= 1 && !witness_radius.has_value(); --exponent) {
      const auto candidate = producing.geometry.radius -
          witness_gap / dyadic_denominator(exponent);
      (void)try_witness(candidate, "outward", exponent);
    }
    for (std::uint32_t exponent = 2;
         exponent <= kWitnessSearchCap &&
             !witness_radius.has_value(); ++exponent) {
      const auto candidate = point_modulus +
          witness_gap / dyadic_denominator(exponent);
      (void)try_witness(candidate, "inward", exponent);
    }
    if (!witness_radius.has_value()) {
      if (certificate_failed_after_disk)
        return ineligible(
            "receiver-center-tail-certificate-fails-after-disk-certification");
      if (!last_accuracy_failure.empty())
        return ineligible(
            "inflated-center-evaluation-fails-relative-accuracy-contract",
            last_accuracy_failure);
      return ineligible(
          "receiver-center-tail-certificate-is-inconclusive");
    }

    std::string local_handle;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->closed)
        throw std::invalid_argument("persistent solver session is closed");
      if (session->locals.size() + session->pending_local_solves + 1 >
          session->local_capacity)
        throw std::invalid_argument(
            "persistent local capacity is exhausted by value transport hop");
      local_handle = "l:" + std::to_string(session->next_local++);
      ++session->pending_local_solves;
    }

    const auto started = std::chrono::steady_clock::now();
    std::shared_ptr<StoredLocalBase> next;
    json::object sealed_lineage;
    try {
      const auto owner_reference = [](
          const RetainedPlanChartBinding& binding,
          const PreparedChartBase& owner) {
        return json::object{
            {"kind", owner.equation_owner_kind()},
            {"handle", owner.equation_owner_handle()},
            {"operator_identity", owner.equation_operator_identity()},
            {"plan_exact_identity", binding.exact_identity},
            {"owner_signature_identity", owner.owner_signature_identity()},
            {"physical_payload_identity", owner.physical_payload_identity()}};
      };
      const auto incoming_record = json::object{
          {"local", incoming->handle()},
          {"chart", incoming->source_chart()},
          {"source_operator_identity", incoming->source_operator_identity()},
          {"checkpoint_identity", incoming->checkpoint_identity()},
          {"epsilon", json::object{
               {"min", incoming_summary.at("epsilon_min")},
               {"max", incoming_summary.at("epsilon_max")}}},
          {"taylor_complete_max", producer_taylor_complete_max},
          {"top_valid", incoming_summary.at("top_valid")},
          {"analytic_metadata", incoming->exact_analytic_metadata()}};
      const auto result_checkpoint = arm_checkpoint_identity(
          checkpoint_root, arm_name, "local", match_index + 1);
      json::object derivation{
          {"schema", "diffexp2-retained-plan-value-handoff-v2"},
          {"capability", "retained-native-regular-value-handoff-v2"},
          {"tile_plan", plan->handle()},
          {"tile_plan_checkpoint_identity", plan->checkpoint_identity()},
          {"tile_plan_provenance_identity", plan->provenance_identity()},
          {"arm", arm_name},
          {"match", match_index},
          {"producing", json::object{
               {"chart", encode_plan_chart(
                    producing, exact_match.producing_chart)},
               {"owner", owner_reference(producing, **producing_owner)}}},
          {"receiving", json::object{
               {"chart", encode_plan_chart(
                    receiving, exact_match.receiving_chart)},
               {"owner", owner_reference(receiving, **receiving_owner)}}},
          {"receiver_center_physical_exact",
           receiving.local_geometry.center_exact},
          {"producing_local_exact", producing_local.str()},
          {"prototype_identity", json::serialize(
               canonical_json_value(value_solver))},
          {"tail_contract", json::object{
               {"mode", "certified-regular-taylor-point-tail-acb-v1"},
               {"producer_taylor_complete_max",
                producer_taylor_complete_max},
               {"receiver_taylor_complete_max",
                receiver_taylor_complete_max},
               {"center_ratio_exact", center_ratio.str()},
               {"producing_point_exact", producing_local.str()},
                    {"producing_chart_radius_exact",
                      producing.geometry.radius.str()},
                    {"witness_radius_exact", witness_radius->str()},
                    {"witness_search_policy",
                     "accuracy-qualified-bidirectional-dyadic-v2"},
                    {"witness_dyadic_direction",
                     witness_dyadic_direction},
                    {"witness_dyadic_exponent",
                     witness_dyadic_exponent},
               {"source_model", checkpoint_regular_tail_model_record(
                    *acb_incoming->tail_model().model)},
               {"certificate", json::object{
                    {"status", tail_majorant_status_name(
                         certified->tail.status)},
                    {"value", checkpoint_error_envelope_record(
                         certified->tail.value)},
                    {"theta", checkpoint_error_envelope_record(
                         certified->tail.theta)},
                    {"disk", json::object{
                         {"witness_radius_exact",
                          certified->tail.disk.witness_radius_exact},
                         {"q_lower_exact",
                          certified->tail.disk.q_lower.dump_exact()},
                         {"ode_norm_upper_exact",
                          certified->tail.disk.ode_norm_upper.dump_exact()},
                         {"cauchy_circle_upper_exact", [&] {
                            json::array values;
                            for (const auto& value :
                                 certified->tail.disk.cauchy_circle_upper)
                              values.emplace_back(value.dump_exact());
                            return values;
                          }()},
                         {"detail", certified->tail.disk.detail}}},
                    {"detail", certified->tail.detail}}},
               {"inflation", json::object{
                    {"gate",
                     "each-component-acb-add-error-mag-by-epsilon-row-v1"},
                    {"retained_value",
                     checkpoint_epsilon_vector_record(retained_value)},
                    {"inflated_value",
                     checkpoint_epsilon_vector_record(
                         certified_handoff.value)}}}}},
          {"accuracy_contract", json::object{
               {"relative_error_max_exact", relative_accuracy_max.str()},
               {"gate",
                "required-prefix-component-radii-lte-threshold-times-max-one-upper-magnitude-v2"},
               {"acb_preflight_required", true}}},
          {"epsilon", json::object{
               {"min", requested_epsilon.min_power},
               {"max", requested_epsilon.complete_max},
               {"required_complete_max", required_complete_max}}},
          {"incoming", incoming_record},
          {"scope", "single-regular-to-regular-transport-hop"},
          {"coefficient_transport", "native-retained-only"},
          {"whole_arm_complete", false}};
      struct ValueHandoffOwners {
        std::shared_ptr<StoredTilePlan> plan;
        std::shared_ptr<StoredLocalBase> incoming;
        std::shared_ptr<PreparedChartBase> receiving;
      };
      auto derivation_owner = std::make_shared<ValueHandoffOwners>(
          ValueHandoffOwners{plan, incoming, *receiving_owner});
      auto chart = std::dynamic_pointer_cast<PreparedChart<ComplexBall>>(
          *receiving_owner);
      if (!chart)
        throw std::invalid_argument(
            "value transport hop receiver has the wrong Acb chart type");
      next = chart->solve_regular_value_handoff(
          local_handle, run_prototype, metadata_prototype, incoming,
          certified_handoff.value, receiving.geometry,
          receiving.prescriptions, requested_epsilon,
          required_complete_max, result_checkpoint, std::move(derivation),
          derivation_owner, chart);
      sealed_lineage = next->seal_plan_match_lineage();
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "value transport hop reservation accounting underflow");
      --session->pending_local_solves;
      throw;
    }

    const auto next_stats = next->stats();
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "value transport hop reservation accounting underflow");
      --session->pending_local_solves;
      if (session->closed)
        throw std::invalid_argument(
            "persistent solver session closed during value transport hop");
      const auto found = session->locals.find(incoming_handle);
      if (found == session->locals.end() ||
          found->second.get() != incoming.get())
        throw std::invalid_argument(
            "value transport hop incoming owner changed before publication");
      if (!session->locals.emplace(local_handle, next).second)
        throw std::logic_error(
            "value transport hop local handle collision");
      ++session->total_local_solves;
      session->total_local_run_parse_ms += next_stats.create_parse_ms;
      session->total_local_kernel_ms += next_stats.create_kernel_ms;
      plan->note_match_advance(arm_name);
    }
    auto next_summary = compact_transport_local_reference(next);
    next_summary["release_via"] = "local";
    return json::object{
        {"status", "ok"}, {"session", session->handle},
        {"capability", "regular-value-transport-hop-v1"},
        {"native_retained", true}, {"json_coefficients", 0},
        {"used", true}, {"reason", "eligible-regular-safe-center"},
        {"arm", arm_name}, {"match", match_index},
        {"tail_witness_direction", witness_dyadic_direction},
        {"tail_witness_dyadic_exponent", witness_dyadic_exponent},
        {"value_hops", 1}, {"basis_matches", 0},
        {"next_local", std::move(next_summary)},
        {"sealed_local_lineage", std::move(sealed_lineage)},
        {"elapsed_ms", std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - started).count()}};
  }

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
    bool terminal_basis_match = false;
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
      terminal_basis_match =
          match_index + 1 == retained.exact.matches.size();
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
    auto attempted_epsilon = requested_epsilon;
    try {
      const auto live_epsilon = live_match_epsilon_intersection(
          requested_epsilon, required_complete_max, incoming, basis);
      attempted_epsilon = live_epsilon;
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
            basis_handles, basis, incoming_handle, incoming, true,
            terminal_basis_match);
        if (auto incomplete = match->incomplete_acb_summary();
            incomplete.has_value()) {
          const auto& residual = as_object(
              incomplete->at("residual"),
              "incomplete consuming-hop Acb residual");
          const auto& complete_window = as_object(
              residual.at("complete_window"),
              "incomplete consuming-hop Acb residual window");
          const auto complete_max = as_i32(
              complete_window.at("max"),
              "incomplete consuming-hop Acb residual maximum");
          const auto additional = std::max<std::int32_t>(
              0, required_complete_max - complete_max);
          const auto complete_through_required =
              residual.at("complete_through_required").as_bool();
          const auto& coefficient_verdicts = as_object(
              residual.at("coefficient_verdicts"),
              "incomplete consuming-hop Acb coefficient verdicts");
          const auto inconclusive_coefficients = as_u64(
              coefficient_verdicts.at("inconclusive"),
              "incomplete consuming-hop inconclusive coefficient count");
          const auto& normal_frame_attempt = as_object(
              incomplete->at("normal_frame_attempt"),
              "incomplete consuming-hop normal-frame attempt");
          const bool propagated_enclosure =
              normal_frame_attempt.if_contains(
                  "physical_clearance_source") != nullptr &&
              normal_frame_attempt.at(
                  "physical_clearance_source").is_string() &&
              normal_frame_attempt.at(
                  "physical_clearance_source").as_string() ==
                  "propagated-enclosure";
          const bool retryable_epsilon =
              !complete_through_required && additional > 0;
          const bool retryable_clearance =
              complete_through_required &&
              inconclusive_coefficients > 0 &&
              !propagated_enclosure;
          const auto& lattice = as_object(
              incomplete->at("exact_lattice"),
              "incomplete consuming-hop exact lattice");
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
          {
            std::lock_guard<std::mutex> lock(session->mutex);
            if (session->pending_local_solves == 0)
              throw std::logic_error(
                  "consuming transport hop reservation accounting underflow");
            --session->pending_local_solves;
          }
          return json::object{
              {"status", "error"},
              {"id", "CPP"},
              {"reason", "acb_match_residual_inconclusive"},
              {"retryable_epsilon_reservoir", retryable_epsilon},
              {"retryable_matching_clearance", retryable_clearance},
              {"retryable_propagated_enclosure",
               propagated_enclosure},
              {"required_additional_epsilon_orders", additional},
              {"arm", arm_name},
              {"match", match_index},
              {"geometry", match->handoff().at("geometry")},
              {"residual", residual},
              {"epsilon", incomplete->at("epsilon")},
              {"refinement", incomplete->at("refinement")},
              {"weight_windows", incomplete->at("weight_windows")},
              {"normal_frame_attempt",
               normal_frame_attempt},
              {"exact_lattice", compact_lattice},
              {"detail",
               retryable_epsilon
                   ? "the Acb match needs a wider private epsilon reservoir before materialization"
                   : propagated_enclosure
                   ? "the Acb match reaches the required epsilon order, but uncertainty propagated by its producer encloses the residual tolerance; receiving Taylor-order and precision retries are not applicable"
                   : "the Acb match reaches the required epsilon order but its finite-Taylor overlap is not accurate enough for the residual tolerance"}};
        }
        next = match->materialize(
            local_handle,
            arm_checkpoint_identity(checkpoint_root, arm_name, "local",
                                    match_index + 1),
            session->precision_bits, match, terminal_basis_match);
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
            session->precision_bits, match, terminal_basis_match);
      }
      require_retained_local_complete_max(
          next, required_complete_max, "consuming transport hop");
      sealed_lineage = next->seal_plan_match_lineage();
      if (terminal_basis_match &&
          match->has_terminal_acb_factorization())
        next->attach_terminal_factorized_owner(
            std::static_pointer_cast<void>(match),
            match->checkpoint_identity());
    } catch (const MatchingArithmeticError& error) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "consuming transport hop reservation accounting underflow");
      --session->pending_local_solves;
      if (session->domain == "acb" &&
          error.code ==
              MatchingArithmeticErrorCode::
                  MaterializedContinuityInconclusive) {
        const auto& arm = plan->arm(arm_name);
        const auto failure_power = error.epsilon_power.has_value()
            ? json::value(*error.epsilon_power)
            : json::value(nullptr);
        const auto failure_component = error.row.has_value()
            ? json::value(*error.row)
            : json::value(nullptr);
        const json::object verdicts{
            {"pass", 0}, {"fail", 0}, {"inconclusive", 1}};
        return json::object{
            {"status", "error"},
            {"id", "CPP"},
            {"reason", "acb_match_residual_inconclusive"},
            {"retryable_epsilon_reservoir", false},
            {"retryable_matching_clearance", true},
            {"required_additional_epsilon_orders", 0},
            {"arm", arm_name},
            {"match", match_index},
            {"geometry", encode_plan_match(arm, match_index)},
            {"residual", json::object{
                 {"status", "materialized-continuity-inconclusive"},
                 {"scope", "materialized-continuity-clearance"},
                 {"complete_through_required", true},
                 {"complete_window", json::object{
                      {"min", attempted_epsilon.min_power},
                      {"max", attempted_epsilon.complete_max}}},
                 {"required_complete_max", required_complete_max},
                 {"coefficient_verdicts", verdicts},
                 {"required_coefficient_verdicts", verdicts},
                 {"failure_epsilon", failure_power},
                 {"failure_component", failure_component},
                 {"detail", error.what()}}},
            {"epsilon", json::object{
                 {"min", attempted_epsilon.min_power},
                 {"max", attempted_epsilon.complete_max},
                 {"required_complete_max", required_complete_max}}},
            {"refinement", refinement},
            {"detail",
             "the Acb match residual passed, but the independently materialized receiving local did not retain the required handoff accuracy"}};
      }
      if (session->domain == "acb" &&
          error.code ==
              MatchingArithmeticErrorCode::UnresolvedDeterminantTail) {
        const auto& arm = plan->arm(arm_name);
        return json::object{
            {"status", "error"},
            {"id", "CPP"},
            {"reason", "acb_match_residual_inconclusive"},
            {"retryable_epsilon_reservoir", true},
            {"retryable_matching_clearance", false},
            {"required_additional_epsilon_orders", 1},
            {"arm", arm_name},
            {"match", match_index},
            {"geometry", encode_plan_match(arm, match_index)},
            {"residual", json::object{
                 {"status", "unresolved-determinant-tail"},
                 {"common_complete_max",
                  error.epsilon_power.has_value()
                      ? json::value(*error.epsilon_power)
                      : json::value(nullptr)},
                 {"required_complete_max", required_complete_max},
                 {"detail", error.what()}}},
            {"epsilon", json::object{
                 {"min", attempted_epsilon.min_power},
                 {"max", attempted_epsilon.complete_max},
                 {"required_complete_max", required_complete_max}}},
            {"refinement", refinement},
            {"detail",
             "the exact saturation determinant remains unresolved at the "
             "private epsilon edge; retry with one additional order"}};
      }
      if (session->domain == "acb" &&
          error.code ==
              MatchingArithmeticErrorCode::InsufficientCompleteWindow &&
          error.epsilon_power.has_value()) {
        const auto complete_max = *error.epsilon_power;
        const auto additional = std::max<std::int32_t>(
            0, required_complete_max - complete_max);
        if (additional > 0) {
          const auto& arm = plan->arm(arm_name);
          const bool materialized_source_incomplete =
              std::string(error.what()).find(
                  "materialized physical source") != std::string::npos;
          return json::object{
              {"status", "error"},
              {"id", "CPP"},
              {"reason", "acb_match_residual_inconclusive"},
              {"retryable_epsilon_reservoir", true},
              {"retryable_matching_clearance", false},
              {"required_additional_epsilon_orders", additional},
              {"arm", arm_name},
              {"match", match_index},
              {"geometry", encode_plan_match(arm, match_index)},
              {"residual", json::object{
                   {"status", materialized_source_incomplete
                                  ? "materialized-source-incomplete"
                                  : "no-common-complete-window"},
                   {"common_complete_max", complete_max},
                   {"required_complete_max", required_complete_max},
                   {"detail", error.what()}}},
              {"epsilon", json::object{
                   {"min", requested_epsilon.min_power},
                   {"max", requested_epsilon.complete_max},
                   {"required_complete_max", required_complete_max}}},
              {"refinement", refinement},
              {"detail",
               materialized_source_incomplete
                   ? "the Acb match residual passed, but its Laurent-weighted "
                     "physical source needs a wider private epsilon reservoir"
                   : "the Acb match needs a wider private epsilon reservoir "
                     "before residual certification can begin"}};
        }
      }
      throw;
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "consuming transport hop reservation accounting underflow");
      --session->pending_local_solves;
      throw;
    }

    json::object response;
    try {
      json::object match_reference{
          {"index", match_index},
          {"checkpoint_identity", match->checkpoint_identity()},
          {"provenance_identity", match->provenance_identity()},
          {"planned_hop", match->handoff()},
          {"sealed_local_lineage", sealed_lineage}};
      if (terminal_basis_match && match->has_acb_match())
        match_reference["diagnostic_native_match_summary"] =
            match->compact_terminal_match_diagnostic();
      auto next_summary = compact_transport_local_reference(next);
      next_summary["release_via"] = "local";
      response = json::object{
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
      // Publication consumes the previous basis owners.  Prove that the
      // entire response fits Boost.JSON before making that destructive state
      // change; otherwise a late serialization failure leaves the native arm
      // advanced while the Wolfram caller correctly believes that the hop
      // failed.
      (void)json::serialize(response);
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_local_solves == 0)
        throw std::logic_error(
            "consuming transport hop response reservation accounting underflow");
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
      ++session->total_transport_framed_basis_hops;
      session->total_local_match_ms += match->elapsed_ms();
      plan->note_match_advance(arm_name);
    }
    return response;
  }

  if (operation == "transport.publish_consumed_states" ||
      operation == "transport.publish_consumed_state") {
    const bool single_arm =
        operation == "transport.publish_consumed_state";
    if (single_arm)
      require_exact_keys(
          root,
          {"schema", "op", "session", "tile_plan",
           "tile_plan_checkpoint_identity", "anchor",
           "anchor_checkpoint_identity", "epsilon", "refinement",
           "checkpoint_policy", "arm", "tile_sources"},
          "native transport.publish_consumed_state request");
    else
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
    const auto parse_side = [&](const std::string& arm,
                                const json::value& raw_side) {
      const auto& side = as_object(
          raw_side, "published consumed transport-arm input");
      require_exact_keys(side, {"tile_sources"},
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
    std::vector<ConsumedStateInput> inputs;
    if (single_arm) {
      const auto arm_name = required_string(root, "arm");
      if (arm_name != "lower" && arm_name != "upper")
        throw std::invalid_argument(
            "published consumed transport arm must be lower or upper");
      inputs.push_back(parse_side(
          arm_name, json::object{{"tile_sources",
                                  root.at("tile_sources")}}));
    } else {
      inputs.push_back(parse_side("lower", root.at("lower")));
      inputs.push_back(parse_side("upper", root.at("upper")));
    }
    const auto plan_handle = required_string(root, "tile_plan");
    const auto anchor_handle = required_string(root, "anchor");
    std::shared_ptr<StoredTilePlan> plan;
    std::shared_ptr<StoredLocalBase> anchor;
    std::vector<std::string> state_handles(inputs.size());
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
              session->pending_transport_states + inputs.size() >
          session->transport_state_capacity)
        throw std::invalid_argument(
            "persistent transport-state capacity is exhausted");
      for (auto& handle : state_handles)
        handle = "transport:" +
          std::to_string(session->next_transport_state++);
      session->pending_transport_states += inputs.size();
    }

    const auto started = std::chrono::steady_clock::now();
    std::vector<std::shared_ptr<StoredTransportArmState>> states(
        inputs.size());
    try {
      for (std::size_t index = 0; index < inputs.size(); ++index) {
        auto& input = inputs[index];
        std::shared_ptr<StoredPlannedMatchHop>
            terminal_factorized_match;
        if (const auto erased =
                input.tile_sources.back()
                    ->terminal_factorized_owner();
            erased != nullptr)
          terminal_factorized_match =
              std::static_pointer_cast<StoredPlannedMatchHop>(
                  erased);
        states[index] = std::make_shared<StoredTransportArmState>(
            state_handles[index],
            checkpoint_root + ":" + input.arm + ":state", input.arm,
            plan, anchor, std::move(input.tile_sources), epsilon_contract.work,
            epsilon_contract.public_required_complete_max,
            epsilon_contract.match_required_complete_max, refinement, 0.0,
            std::move(terminal_factorized_match));
      }
    } catch (...) {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_transport_states < inputs.size())
        throw std::logic_error(
            "published consumed-state reservation accounting underflow");
      session->pending_transport_states -= inputs.size();
      throw;
    }

    {
      std::lock_guard<std::mutex> lock(session->mutex);
      if (session->pending_transport_states < inputs.size())
        throw std::logic_error(
            "published consumed-state reservation accounting underflow");
      session->pending_transport_states -= inputs.size();
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
      session->total_transport_arm_marches += inputs.size();
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
        {"capability", single_arm
             ? "published-consumed-transport-arm-state-v1"
             : "published-consumed-transport-arm-states-v1"},
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
            anchor_binding.local_geometry, anchor_binding.prescriptions,
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
            checkpoint_root, true);
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
          anchor_binding.local_geometry, anchor_binding.prescriptions,
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
        checkpoint_root, true);
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

  if (operation == "transport.diagnostic_rebuild_terminal_match") {
    const bool scan_points =
        root.if_contains("receiving_local_points") != nullptr;
    const bool override_precision =
        root.if_contains("precision_bits") != nullptr;
    if (scan_points && override_precision)
      require_exact_keys(
          root, {"schema", "op", "session", "transport_state",
                 "receiving_local_points", "precision_bits"},
          "native terminal-match diagnostic scan request");
    else if (scan_points)
      require_exact_keys(
          root, {"schema", "op", "session", "transport_state",
                 "receiving_local_points"},
          "native terminal-match diagnostic scan request");
    else if (override_precision)
      require_exact_keys(
          root, {"schema", "op", "session", "transport_state",
                 "precision_bits"},
          "native terminal-match diagnostic rebuild request");
    else
      require_exact_keys(
          root, {"schema", "op", "session", "transport_state"},
          "native terminal-match diagnostic rebuild request");
    const auto* diagnostics =
        std::getenv("DE2_DIAGNOSTIC_TERMINAL_STATE");
    if (diagnostics == nullptr || std::string(diagnostics) != "1")
      throw std::invalid_argument(
          "terminal-match diagnostic rebuild requires "
          "DE2_DIAGNOSTIC_TERMINAL_STATE=1");
    const auto diagnostic_precision = override_precision
        ? static_cast<slong>(
              as_u32(root.at("precision_bits"),
                     "terminal diagnostic precision bits"))
        : session->precision_bits;
    if (diagnostic_precision < 64 ||
        diagnostic_precision > 16384)
      throw std::invalid_argument(
          "terminal diagnostic precision must lie in 64..16384 bits");
    std::shared_ptr<StoredTransportArmState> state;
    {
      std::lock_guard<std::mutex> lock(session->mutex);
      const auto found = session->transport_states.find(
          required_string(root, "transport_state"));
      if (found == session->transport_states.end())
        throw std::invalid_argument(
            "unknown native transport-arm state for terminal rematch");
      state = found->second;
    }
    const auto previous = state->terminal_factorized_match();
    if (!previous || !previous->has_acb_match())
      throw std::invalid_argument(
          "terminal rematch requires one retained Acb terminal match");
    const auto plan = state->plan_owner();
    const auto& handoff = previous->handoff();
    const auto arm_name = required_string(handoff, "arm");
    const auto match_index = static_cast<std::size_t>(
        as_u64(handoff.at("match"),
               "terminal diagnostic rematch index"));
    const auto previous_summary =
        previous->compact_terminal_match_diagnostic();
    const auto& epsilon = as_object(
        previous_summary.at("epsilon"),
        "terminal diagnostic rematch epsilon");
    const auto& previous_refinement = as_object(
        previous_summary.at("refinement"),
        "terminal diagnostic rematch refinement");
    const auto checkpoint_identity =
        previous->checkpoint_identity();
    const auto previous_acb =
        std::dynamic_pointer_cast<StoredRefinedAcbMatch>(
            previous->native_match());
    if (!previous_acb)
      throw std::logic_error(
          "terminal diagnostic rematch lost its retained Acb proof");
    const auto retained_singular_request =
        previous_acb->retained_singular_saturation_request();
    if (scan_points) {
      const auto& raw_points = as_array(
          root.at("receiving_local_points"),
          "terminal diagnostic receiving-local scan points");
      if (raw_points.empty() || raw_points.size() > 16)
        throw std::invalid_argument(
            "terminal diagnostic scan requires 1..16 receiving-local points");
      const auto& retained_arm = plan->arm(arm_name);
      if (match_index >= retained_arm.exact.matches.size())
        throw std::logic_error(
            "terminal diagnostic scan match index is outside its plan");
      const auto& exact_match =
          retained_arm.exact.matches[match_index];
      const auto& producing =
          retained_arm.exact.charts.at(
              exact_match.producing_chart);
      const auto& receiving =
          retained_arm.exact.charts.at(
              exact_match.receiving_chart);
      std::vector<std::string> basis_handles;
      json::array basis_checkpoints;
      basis_handles.reserve(previous->basis_owners().size());
      basis_checkpoints.reserve(previous->basis_owners().size());
      for (const auto& owner : previous->basis_owners()) {
        basis_handles.push_back(owner->handle());
        basis_checkpoints.emplace_back(
            owner->checkpoint_identity());
      }
      json::array scans;
      scans.reserve(raw_points.size());
      for (std::size_t index = 0; index < raw_points.size();
           ++index) {
        if (!raw_points[index].is_string())
          throw std::invalid_argument(
              "terminal diagnostic receiving-local points must be exact rational strings");
        const Rational receiving_local(
            std::string(raw_points[index].as_string()));
        const auto physical =
            receiving.center + receiving.scale * receiving_local;
        const auto incoming_local =
            (physical - producing.center) / producing.scale;
        json::object item{
            {"receiving_local", receiving_local.str()},
            {"incoming_local", incoming_local.str()},
            {"physical", physical.str()}};
        try {
          auto singular_request = retained_singular_request;
          const auto scan_checkpoint =
              checkpoint_identity + ":point-scan:" +
              std::to_string(index);
          singular_request["match_checkpoint_identity"] =
              scan_checkpoint;
          singular_request["receiving_basis_point_exact"] =
              receiving_local.str();
          singular_request["receiving_basis_point_sign"] =
              receiving_local.sign();
          singular_request["physical_match_point_exact"] =
              physical.str();
          json::object kernel_request{
              {"basis", [&]() {
                 json::array values;
                 for (const auto& handle : basis_handles)
                   values.emplace_back(handle);
                 return values;
               }()},
              {"incoming", previous->incoming_owner()->handle()},
              {"basis_chart",
               previous->basis_owners().front()->source_chart()},
              {"incoming_chart",
               previous->incoming_owner()->source_chart()},
              {"basis_point",
               json::object{{"exact", receiving_local.str()}}},
              {"incoming_point",
               json::object{{"exact", incoming_local.str()}}},
              {"epsilon", json::object{
                   {"min", epsilon.at("min")},
                   {"max", epsilon.at("max")},
                   {"required_complete_max",
                    epsilon.at("required_complete_max")}}},
              {"basis_checkpoint_identities",
               basis_checkpoints},
              {"incoming_checkpoint_identity",
               previous->incoming_owner()
                   ->checkpoint_identity()},
              {"checkpoint_identity", scan_checkpoint},
              {"refinement", json::object{
                   {"relative_tolerance",
                    previous_refinement.at(
                        "relative_tolerance")},
                   {"max_steps",
                    previous_refinement.at("max_steps")}}},
              {"native_singular_scc_saturation",
               singular_request}};
          if (const auto sign =
                  previous_acb
                      ->effective_basis_imaginary_sign();
              sign.has_value())
            kernel_request["basis_imaginary_sign"] = *sign;
          if (const auto sign =
                  previous_acb
                      ->effective_incoming_imaginary_sign();
              sign.has_value())
            kernel_request["incoming_imaginary_sign"] = *sign;
          const auto rebuilt = build_refined_acb_match(
              previous->native_match()->handle() +
                  ":strict-point-scan:" +
                  std::to_string(index),
              kernel_request, basis_handles,
              previous->basis_owners(),
              previous->incoming_owner()->handle(),
              previous->incoming_owner(), diagnostic_precision,
              checkpoint_configuration_identity(*session),
              singular_request, true);
          item["status"] = "ok";
          item["match"] =
              rebuilt->compact_terminal_diagnostic_summary();
        } catch (const std::exception& error) {
          item["status"] = "error";
          item["detail"] = error.what();
        }
        scans.push_back(std::move(item));
      }
      return json::object{
          {"status", "ok"},
          {"session", session->handle},
          {"scan", std::move(scans)}};
    }
    json::object match_request{
        {"arm", arm_name},
        {"match", match_index},
        {"epsilon", json::object{
             {"min", epsilon.at("min")},
             {"max", epsilon.at("max")},
             {"required_complete_max",
              epsilon.at("required_complete_max")}}},
        {"checkpoint_identity", checkpoint_identity}};
    match_request["native_singular_scc_saturation"] =
        retained_singular_request;
    match_request["refinement"] = json::object{
        {"relative_tolerance",
         previous_refinement.at("relative_tolerance")},
        {"max_steps", previous_refinement.at("max_steps")}};
    std::vector<std::string> basis_handles;
    basis_handles.reserve(previous->basis_owners().size());
    for (const auto& owner : previous->basis_owners())
      basis_handles.push_back(owner->handle());
    const auto rebuilt = build_planned_match_hop(
        previous->native_match()->handle() + ":strict-diagnostic-rematch",
        match_request, session->domain, diagnostic_precision,
        checkpoint_configuration_identity(*session), plan,
        basis_handles, previous->basis_owners(),
        previous->incoming_owner()->handle(),
        previous->incoming_owner(), true, true,
        retained_singular_request);
    return json::object{
        {"status", "ok"},
        {"session", session->handle},
        {"previous", previous_summary},
        {"rebuilt", rebuilt->compact_terminal_match_diagnostic()}};
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
      std::optional<BoundedDivergentCancellation>
          divergent_cancellation;
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
      if (observable.if_contains("divergent_cancellation") != nullptr)
        require_exact_keys(
            observable,
            {"identity", "checkpoint_identity", "integrand_rows",
             "epsilon", "tail_policy", "divergent_cancellation"},
            "native transport observable");
      else
        require_exact_keys(
            observable,
            {"identity", "checkpoint_identity", "integrand_rows",
             "epsilon", "tail_policy"},
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
      if (const auto* policy =
              observable.if_contains("divergent_cancellation")) {
        if (session->domain != "acb")
          throw std::invalid_argument(
              "bounded divergent cancellation is restricted to Acb transport contractions");
        parsed.divergent_cancellation =
            parse_bounded_divergent_cancellation(
                *policy,
                "native transport divergent-cancellation policy");
      }
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
        input.divergent_cancellation =
            observable.divergent_cancellation;
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
        if (observable.divergent_cancellation.has_value())
          input.aggregate_record["divergent_cancellation"] =
              encode_bounded_divergent_cancellation(
                  *observable.divergent_cancellation);
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
    const bool has_divergent_cancellation =
        observable.if_contains("divergent_cancellation") != nullptr;
    const bool has_publication_tolerance =
        observable.if_contains("publication_relative_tolerance") !=
        nullptr;
    if (has_divergent_cancellation && has_publication_tolerance)
      require_exact_keys(
          observable,
          {"identity", "checkpoint_identity", "epsilon", "tail_policy",
           "divergent_cancellation",
           "publication_relative_tolerance"},
          "streamed transport-pair observable");
    else if (has_divergent_cancellation)
      require_exact_keys(
          observable,
          {"identity", "checkpoint_identity", "epsilon", "tail_policy",
           "divergent_cancellation"},
          "streamed transport-pair observable");
    else if (has_publication_tolerance)
      require_exact_keys(
          observable,
          {"identity", "checkpoint_identity", "epsilon", "tail_policy",
           "publication_relative_tolerance"},
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
    std::optional<Magnitude> publication_relative_tolerance;
    if (has_publication_tolerance) {
      const auto text = required_string(
          observable, "publication_relative_tolerance");
      publication_relative_tolerance = Magnitude::decimal(text);
      if (text.empty() ||
          !publication_relative_tolerance->is_finite() ||
          publication_relative_tolerance->is_zero() ||
          Magnitude::one() <= *publication_relative_tolerance)
        throw std::invalid_argument(
            "streamed transport-pair publication tolerance must be one finite decimal strictly between zero and one");
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
          divergent_cancellation, publication_relative_tolerance);
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
    json::object progress;
    try {
      progress = stream->add_tile(side, tile, row);
    } catch (const MatchingArithmeticError& error) {
      if (error.code ==
          MatchingArithmeticErrorCode::TerminalOutputInconclusive)
        return json::object{
            {"status", "error"},
            {"id", "CPP"},
            {"reason", "terminal_output_ball_inconclusive"},
            {"retryable_level_accuracy", true},
            {"request_index", 0},
            {"failure_functional",
             error.row.has_value()
                 ? json::value(*error.row)
                 : json::value(nullptr)},
            {"failure_epsilon",
             error.epsilon_power.has_value()
                 ? json::value(*error.epsilon_power)
                 : json::value(nullptr)},
            {"required_additional_digits",
             error.required_additional_digits.has_value()
                 ? json::value(*error.required_additional_digits)
                 : json::value(nullptr)},
            {"side", side_name},
            {"tile", tile},
            {"detail", error.what()}};
      if (error.code !=
              MatchingArithmeticErrorCode::InsufficientCompleteWindow ||
          !error.row.has_value() || *error.row == 0)
        throw;
      return json::object{
          {"status", "error"},
          {"id", "CPP"},
          {"reason", "acb_match_residual_inconclusive"},
          {"retryable_epsilon_reservoir", true},
          {"retryable_matching_clearance", false},
          {"required_additional_epsilon_orders", *error.row},
          {"side", side_name},
          {"tile", tile},
          {"common_complete_max",
           error.epsilon_power.has_value()
               ? json::value(*error.epsilon_power)
               : json::value(nullptr)},
          {"detail", error.what()}};
    }
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
    TransportPairObservableStream::FinishResult finished;
    try {
      finished = stream->finish();
      if (stream->publication_relative_tolerance().has_value())
        require_transport_line_publication_accuracy(
            finished.line->result(), stream->epsilon(),
            *stream->publication_relative_tolerance(),
            "streamed final transport-pair observable", 0);
    } catch (const MatchingArithmeticError& error) {
      release_reservation();
      if (error.code != MatchingArithmeticErrorCode::
                            TerminalOutputInconclusive)
        throw;
      return json::object{
          {"status", "error"},
          {"id", "CPP"},
          {"reason", "terminal_output_ball_inconclusive"},
          {"retryable_level_accuracy", true},
          {"request_index", 0},
          {"failure_functional",
           error.row.has_value()
               ? json::value(*error.row)
               : json::value(nullptr)},
          {"failure_epsilon",
               error.epsilon_power.has_value()
                   ? json::value(*error.epsilon_power)
                   : json::value(nullptr)},
          {"required_additional_digits",
           error.required_additional_digits.has_value()
               ? json::value(*error.required_additional_digits)
               : json::value(nullptr)},
          {"scope", "final-paired-line"},
          {"conditioning", finished.conditioning},
          {"detail", error.what()}};
    }
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
        {"conditioning", finished.conditioning},
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
    const bool has_threads = root.if_contains("threads") != nullptr;
    if (has_threads)
      require_exact_keys(
          root,
          {"schema", "op", "session", "lower", "upper",
           "checkpoint_policy", "observables", "threads"},
          "native transport.contract_pair request");
    else
      require_exact_keys(
          root,
          {"schema", "op", "session", "lower", "upper",
           "checkpoint_policy", "observables"},
          "native transport.contract_pair request");
    if (session->domain == "symbolic")
      throw std::invalid_argument(
          "native transport-pair contraction requires rational or Acb coefficients");
    const auto requested_threads = has_threads
        ? as_u32(root.at("threads"),
                 "native transport-pair contraction threads")
        : std::uint32_t{1};
    if (requested_threads == 0)
      throw std::invalid_argument(
          "native transport-pair contraction threads must be positive");

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
      std::optional<Magnitude> publication_relative_tolerance;
      std::string publication_relative_tolerance_text;
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
      const bool has_divergent_cancellation =
          observable.if_contains("divergent_cancellation") != nullptr;
      const bool has_publication_tolerance =
          observable.if_contains("publication_relative_tolerance") !=
          nullptr;
      if (has_divergent_cancellation && has_publication_tolerance)
        require_exact_keys(
            observable,
            {"identity", "checkpoint_identity", "lower_integrand_rows",
             "upper_integrand_rows", "epsilon", "tail_policy",
             "divergent_cancellation",
             "publication_relative_tolerance"},
            "native transport-pair observable");
      else if (has_divergent_cancellation)
        require_exact_keys(
            observable,
            {"identity", "checkpoint_identity", "lower_integrand_rows",
             "upper_integrand_rows", "epsilon", "tail_policy",
             "divergent_cancellation"},
            "native transport-pair observable");
      else if (has_publication_tolerance)
        require_exact_keys(
            observable,
            {"identity", "checkpoint_identity", "lower_integrand_rows",
             "upper_integrand_rows", "epsilon", "tail_policy",
             "publication_relative_tolerance"},
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
      if (has_publication_tolerance) {
        parsed.publication_relative_tolerance_text = required_string(
            observable, "publication_relative_tolerance");
        parsed.publication_relative_tolerance = Magnitude::decimal(
            parsed.publication_relative_tolerance_text);
        if (parsed.publication_relative_tolerance_text.empty() ||
            !parsed.publication_relative_tolerance->is_finite() ||
            parsed.publication_relative_tolerance->is_zero() ||
            Magnitude::one() <=
                *parsed.publication_relative_tolerance)
          throw std::invalid_argument(
              "native transport-pair publication tolerance must be one finite decimal strictly between zero and one");
      }
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
    const auto worker_count = std::min<std::size_t>(
        {static_cast<std::size_t>(requested_threads),
         observable_count, kMaxPersistentBatchThreads});
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
          input.publication_relative_tolerance =
              observable.publication_relative_tolerance;
          input.factorize_ordinary_stored_rows =
              !observable.publication_relative_tolerance.has_value();
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
    const auto insufficient_window_response =
        [](const MatchingArithmeticError& error) {
          return json::object{
              {"status", "error"},
              {"id", "CPP"},
              {"reason", "acb_match_residual_inconclusive"},
              {"retryable_epsilon_reservoir", true},
              {"retryable_matching_clearance", false},
              {"required_additional_epsilon_orders",
               error.row.has_value()
                   ? json::value(*error.row)
                   : json::value(nullptr)},
              {"request_index",
               error.column.has_value()
                   ? json::value(*error.column)
                   : json::value(nullptr)},
              {"common_complete_max",
               error.epsilon_power.has_value()
                   ? json::value(*error.epsilon_power)
                   : json::value(nullptr)},
              {"scope", "paired-observable-arm-contraction"},
              {"detail", error.what()}};
        };

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
      // Keep the two arms sequential so their largest projection scratch
      // allocations never overlap.  Within one arm, distinct observables
      // are mathematically and operationally independent; run only the
      // explicitly bounded number in parallel and retain deterministic
      // request-order publication and failure selection.
      try {
        for (std::size_t side = 0; side < 2; ++side) {
          if (session->domain == "acb")
            ComplexBall::set_precision(session->precision_bits);
          const auto started = std::chrono::steady_clock::now();
          if (worker_count == 1) {
            arm_results[side] = contract_transport_arm(
                session->domain, session->precision_bits,
                states[side]->plan_owner(),
                side == 0 ? "lower" : "upper",
                states[side]->tile_sources(), arm_inputs[side],
                states[side]);
          } else {
            arm_results[side].resize(observable_count);
            std::atomic<std::size_t> next_observable{0};
            std::vector<std::exception_ptr> failures(observable_count);
            std::vector<std::jthread> workers;
            workers.reserve(worker_count);
            for (std::size_t worker = 0; worker < worker_count; ++worker)
              workers.emplace_back([&, side] {
                if (session->domain == "acb")
                  ComplexBall::set_precision(session->precision_bits);
                while (true) {
                  const auto index =
                      next_observable.fetch_add(1);
                  if (index >= observable_count) return;
                  try {
                    std::vector<TransportObservableContractionInput>
                        singleton{arm_inputs[side][index]};
                    auto result = contract_transport_arm(
                        session->domain, session->precision_bits,
                        states[side]->plan_owner(),
                        side == 0 ? "lower" : "upper",
                        states[side]->tile_sources(), singleton,
                        states[side]);
                    if (result.size() != 1)
                      throw std::logic_error(
                          "parallel transport-pair observable returned the wrong result count");
                    arm_results[side][index] =
                        std::move(result.front());
                  } catch (const MatchingArithmeticError& error) {
                    if (error.code == MatchingArithmeticErrorCode::
                                          TerminalOutputInconclusive ||
                        error.code == MatchingArithmeticErrorCode::
                                          InsufficientCompleteWindow)
                      failures[index] = std::make_exception_ptr(
                          MatchingArithmeticError(
                              error.code, error.what(), error.row,
                              index, error.epsilon_power,
                              error.required_additional_digits));
                    else
                      failures[index] =
                          std::current_exception();
                  } catch (...) {
                    failures[index] = std::current_exception();
                  }
                }
              });
            for (auto& worker : workers) worker.join();
            const auto failed = std::find_if(
                failures.begin(), failures.end(),
                [](const auto& error) {
                  return error != nullptr;
                });
            if (failed != failures.end())
              std::rethrow_exception(*failed);
          }
          arm_elapsed_ms[side] =
              std::chrono::duration<double, std::milli>(
                  std::chrono::steady_clock::now() - started)
                  .count();
        }
      } catch (const MatchingArithmeticError& error) {
        release_reservation();
        if (error.code == MatchingArithmeticErrorCode::
                              InsufficientCompleteWindow &&
            error.row.has_value() && *error.row > 0)
          return insufficient_window_response(error);
        if (error.code != MatchingArithmeticErrorCode::
                              TerminalOutputInconclusive)
          throw;
        return json::object{
            {"status", "error"},
            {"id", "CPP"},
            {"reason", "terminal_output_ball_inconclusive"},
            {"retryable_level_accuracy", true},
            {"request_index",
             error.column.has_value()
                 ? json::value(*error.column)
                 : json::value(nullptr)},
            {"failure_functional",
             error.row.has_value()
                 ? json::value(*error.row)
                 : json::value(nullptr)},
            {"failure_epsilon",
             error.epsilon_power.has_value()
                 ? json::value(*error.epsilon_power)
                 : json::value(nullptr)},
            {"required_additional_digits",
             error.required_additional_digits.has_value()
                 ? json::value(*error.required_additional_digits)
                 : json::value(nullptr)},
            {"detail", error.what()}};
      } catch (...) {
        release_reservation();
        throw;
      }
    }

    std::vector<std::shared_ptr<StoredLineResult>> combined_lines;
    combined_lines.reserve(observable_count);
    std::vector<json::object> pair_conditioning;
    pair_conditioning.reserve(observable_count);
    std::array<std::uint64_t, 2>
        superseded_tile_integrations{0, 0};
    double superseded_tile_integration_ms = 0.0;
    std::vector<double> direct_probe_elapsed_ms(
        observable_count, 0.0);
    std::size_t ordinary_factorization_retries = 0;
    for (std::size_t index = 0; index < observable_count; ++index) {
      if (arm_results[0].size() != observable_count ||
          arm_results[1].size() != observable_count)
        throw std::logic_error(
            "native transport-pair contraction returned the wrong arm result count");
      const auto build_combined = [&]() {
        auto aggregate_record = pair_aggregate_records[index];
        const auto line_handle = required_string(
            aggregate_record, "line_handle");
        aggregate_record.erase("line_handle");
        return build_compact_transport_pair_observable_line(
            line_handle,
            pending_observables[index].checkpoint_identity,
            std::move(aggregate_record), states[0], states[1],
            arm_results[0][index].aggregate,
            arm_results[1][index].aggregate);
      };
      auto combined = build_combined();
      if (combined->result().value.epsilon.complete_max <
          pending_observables[index].epsilon.required_complete_max)
        throw std::domain_error(
            "native transport-pair aggregate does not cover its required epsilon maximum");
      auto conditioning =
          encode_transport_pair_conditioning_diagnostics(
              arm_results[0][index].aggregate->result(),
              arm_results[1][index].aggregate->result(),
              combined->result(),
              std::move(arm_results[0][index].tile_values),
              std::move(arm_results[1][index].tile_values));
      if (pending_observables[index]
              .publication_relative_tolerance.has_value()) {
        try {
          require_transport_line_publication_accuracy(
              combined->result(), pending_observables[index].epsilon,
              *pending_observables[index]
                   .publication_relative_tolerance,
              "final transport-pair observable", index);
        } catch (const MatchingArithmeticError& error) {
          if (error.code != MatchingArithmeticErrorCode::
                                TerminalOutputInconclusive)
            throw;
          // The direct stored-row integral is rigorous but may lose the
          // correlation of one uncertain center ball reused throughout its
          // Taylor recurrence. Recompute only this failing logical
          // observable with the expensive factorized ordinary-row
          // reassociation. Terminal singular rows retain their established
          // composed/factorized route in both passes.
          direct_probe_elapsed_ms[index] =
              arm_results[0][index].elapsed_ms +
              arm_results[1][index].elapsed_ms;
          for (std::size_t side = 0; side < 2; ++side) {
            superseded_tile_integrations[side] +=
                arm_results[side][index].tile_integrations;
            superseded_tile_integration_ms +=
                arm_results[side][index].tile_integration_ms;
            auto retry_input = arm_inputs[side][index];
            retry_input.factorize_ordinary_stored_rows = true;
            const auto retry_started =
                std::chrono::steady_clock::now();
            try {
              auto retried = contract_transport_arm(
                  session->domain, session->precision_bits,
                  states[side]->plan_owner(),
                  side == 0 ? "lower" : "upper",
                  states[side]->tile_sources(),
                  std::vector<TransportObservableContractionInput>{
                      std::move(retry_input)},
                  states[side]);
              if (retried.size() != 1)
                throw std::logic_error(
                    "factorized transport-pair retry returned the wrong result count");
              arm_results[side][index] =
                  std::move(retried.front());
            } catch (const MatchingArithmeticError& retry_error) {
              if (retry_error.code ==
                      MatchingArithmeticErrorCode::
                          TerminalOutputInconclusive ||
                  retry_error.code ==
                      MatchingArithmeticErrorCode::
                          InsufficientCompleteWindow)
                throw MatchingArithmeticError(
                    retry_error.code, retry_error.what(),
                    retry_error.row, index,
                    retry_error.epsilon_power,
                    retry_error.required_additional_digits);
              throw;
            }
            arm_elapsed_ms[side] +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() -
                    retry_started)
                    .count();
          }
          ++ordinary_factorization_retries;
          combined = build_combined();
          conditioning =
              encode_transport_pair_conditioning_diagnostics(
                  arm_results[0][index].aggregate->result(),
                  arm_results[1][index].aggregate->result(),
                  combined->result(),
                  std::move(
                      arm_results[0][index].tile_values),
                  std::move(
                      arm_results[1][index].tile_values));
          try {
            require_transport_line_publication_accuracy(
                combined->result(),
                pending_observables[index].epsilon,
                *pending_observables[index]
                     .publication_relative_tolerance,
                "final transport-pair observable after ordinary-row "
                "factorization",
                index);
          } catch (const MatchingArithmeticError& retry_error) {
            release_reservation();
            if (retry_error.code !=
                MatchingArithmeticErrorCode::
                    TerminalOutputInconclusive)
              throw;
            return json::object{
                {"status", "error"},
                {"id", "CPP"},
                {"reason",
                 "terminal_output_ball_inconclusive"},
                {"retryable_level_accuracy", true},
                {"request_index", index},
                {"failure_functional",
                 retry_error.row.has_value()
                     ? json::value(*retry_error.row)
                     : json::value(nullptr)},
                {"failure_epsilon",
                 retry_error.epsilon_power.has_value()
                     ? json::value(
                           *retry_error.epsilon_power)
                     : json::value(nullptr)},
                {"required_additional_digits",
                 retry_error.required_additional_digits
                         .has_value()
                     ? json::value(*retry_error
                                        .required_additional_digits)
                     : json::value(nullptr)},
                {"scope", "final-paired-line"},
                {"ordinary_factorization_retry", true},
                {"conditioning", std::move(conditioning)},
                {"detail", retry_error.what()}};
          }
        }
      }
      if (pending_observables[index].tail_policy ==
              TransportTailPolicy::Require &&
          combined->result().scope !=
              LineIntegrationScope::FullLocalWithCertifiedTail)
        throw std::domain_error(
            "native transport-pair contraction requires certified full-local tails on both arms");
      pair_conditioning.push_back(
          std::move(conditioning));
      combined_lines.push_back(std::move(combined));
    }
    const auto pair_wall_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - operation_started).count();

    std::uint64_t tile_integrations =
        superseded_tile_integrations[0] +
        superseded_tile_integrations[1];
    double tile_integration_ms =
        superseded_tile_integration_ms;
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
      for (std::size_t side = 0; side < 2; ++side)
        for (std::uint64_t tile = 0;
             tile < superseded_tile_integrations[side]; ++tile)
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
          {"conditioning", pair_conditioning[index]},
          {"lower_tiles", arm_results[0][index].tile_integrations},
          {"upper_tiles", arm_results[1][index].tile_integrations},
          {"elapsed_ms", arm_results[0][index].elapsed_ms +
                             arm_results[1][index].elapsed_ms +
                             direct_probe_elapsed_ms[index] +
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
        {"requested_observable_threads", requested_threads},
        {"observable_worker_threads", worker_count},
        {"ordinary_factorization_retries",
         ordinary_factorization_retries},
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
      std::string publication_relative_tolerance_text;
      Magnitude publication_relative_tolerance = Magnitude::zero();
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
          {"identity", "checkpoint_identity", "integrand_row", "epsilon",
           "publication_relative_tolerance"},
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
      parsed.publication_relative_tolerance_text = required_string(
          observable, "publication_relative_tolerance");
      if (parsed.publication_relative_tolerance_text.empty())
        throw std::invalid_argument(
            "native transport endpoint publication tolerance cannot be empty");
      parsed.publication_relative_tolerance = Magnitude::decimal(
          parsed.publication_relative_tolerance_text);
      if (!parsed.publication_relative_tolerance.is_finite() ||
          parsed.publication_relative_tolerance.is_zero() ||
          Magnitude::one() <=
              parsed.publication_relative_tolerance)
        throw std::invalid_argument(
            "native transport endpoint publication tolerance must be finite and strictly between zero and one");
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
    for (std::size_t index = 0; index < observable_count; ++index) {
      const auto& observable = pending_observables[index];
      try {
        endpoints.push_back(build_transport_endpoint_row(
            observable.endpoint_handle, observable.checkpoint_identity,
            observable.identity, observable.row, observable.epsilon,
            observable.epsilon_record,
            observable.publication_relative_tolerance,
            observable.publication_relative_tolerance_text,
            observable.projected_handle,
            observable.projected_checkpoint_identity, session->domain,
            session->precision_bits, state, binding));
      } catch (const MatchingArithmeticError& error) {
        if (error.code !=
            MatchingArithmeticErrorCode::
                TerminalOutputInconclusive)
          throw;
        return json::object{
            {"status", "error"},
            {"id", "CPP"},
            {"reason", "terminal_output_ball_inconclusive"},
            {"retryable_level_accuracy", true},
            {"request_index", index},
            {"failure_functional",
             error.row.has_value()
                 ? json::value(*error.row)
                 : json::value(nullptr)},
            {"failure_epsilon",
             error.epsilon_power.has_value()
                 ? json::value(*error.epsilon_power)
                 : json::value(nullptr)},
            {"required_additional_digits",
             error.required_additional_digits.has_value()
                 ? json::value(*error.required_additional_digits)
                 : json::value(nullptr)},
            {"detail", error.what()}};
      }
    }
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
